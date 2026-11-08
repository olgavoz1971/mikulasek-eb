"""Forced second profile component and heavier minima on NSV 807, sector 30, 120 s.

Takes the configuration chosen by ``nsv807_s30_tiered.py`` (read from its
``chosen`` folder) and the same cleaned light curve, and fits four variants:

* second component on the primary only, or on both minima (``force_two_component``
  skips only the rms-gain condition; the other acceptance conditions stay);
* catalogue weights, or weight 3 inside 0.12 cycles of both minima
  (``eclipse_weight``; the points count ninefold).

No condition is loosened. If the package rejects the second component, the figure
shows the single-component fit and the title says "not adopted". Here "adopted"
means that the second component was adopted.

Run ``nsv807_s30_tiered.py`` first, then from ``examples``::

    python nsv807_s30_forced.py [--candidate p1_h23_b0]

Writes ``results/nsv807/s30_forced/``: one figure per variant (``fits``,
``trends``) and ``variants.csv``.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import logging
from pathlib import Path

import mikulasek_eb as eb
from nsv807_s30_tiered import (
    DATA,
    FILES,
    OUT as CHOSEN_RUN,
    figure_title,
    load_sap,
    remove_outliers,
    save_figures,
)

logger = logging.getLogger(__name__)

OUT = CHOSEN_RUN.parent / "s30_forced"
ECLIPSE_WEIGHT = 3.0
ECLIPSE_PHASE_HALF = 0.12


def chosen_config(stem: str, candidate: str | None) -> eb.Config:
    """Read the configuration of the chosen candidate or of a named one.

    Args:
        stem (str): Light curve file stem.
        candidate (str | None): Name of a candidate saved in ``candidates``, for
            example ``p1_h23_b0``. None means the chosen candidate.

    Returns:
        Config: Configuration stored in the one solution file of ``chosen``.

    Raises:
        FileNotFoundError: If ``chosen`` has no solution for this light curve.
    """
    folder = CHOSEN_RUN / ("chosen" if candidate is None else "candidates")
    pattern = f"{stem}__*.json" if candidate is None else f"{stem}__{candidate}.json"
    found = sorted(folder.glob(pattern))
    if len(found) != 1:
        raise FileNotFoundError(f"expected one chosen solution for {stem}: {found}")
    return eb.ModelSolution.load(found[0]).config


def main() -> None:
    """Fit the four variants and write figures and a table."""
    stem = Path(FILES[0]).stem
    payload = remove_outliers(load_sap(DATA / FILES[0]), stem)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--candidate", help="candidate name, for example p1_h23_b0")
    args = parser.parse_args()
    chosen = chosen_config(stem, args.candidate)
    out_dir = OUT if args.candidate is None else OUT / args.candidate
    rows = []
    for eclipses, where in (((1,), "primary"), ((1, 2), "both")):
        for weight, heavy in ((1.0, "w1"), (ECLIPSE_WEIGHT, "w3")):
            policy = dataclasses.replace(
                chosen.policy,
                allow_two_component=True,
                force_two_component=True,
                two_comp_eclipses=eclipses,
                eclipse_weight=weight,
                eclipse_phase_half=ECLIPSE_PHASE_HALF,
            )
            name = f"2c_{where}_{heavy}"
            out = eb.fit(
                t_mjd,
                payload["mag"],
                payload["sig"],
                payload["period"],
                payload["epoch_jd"],
                dataclasses.replace(chosen, policy=policy).validate(),
                data_id=f"{stem}__{name}",
            )
            stats = out.solution.stats
            trial = stats.get("two_component_trial") or {}
            adopted = bool(trial.get("adopted"))
            rms = 1e3 * float(stats["rms_mag"])
            save_figures(
                out_dir,
                f"{stem}__{name}",
                payload,
                out,
                figure_title(stem, name, stats["status"], adopted, rms),
            )
            out.solution.save(out_dir / f"{stem}__{name}.json")
            rows.append(
                {
                    "variant": name,
                    "status": stats["status"],
                    "second_component_adopted": int(adopted),
                    "reason": trial.get("reason", ""),
                    "rms_mmag": f"{rms:.2f}",
                }
            )
            logger.info("%s adopted=%s rms %.2f mmag", name, adopted, rms)
    with open(out_dir / "variants.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    logger.info("wrote %s", OUT)


if __name__ == "__main__":
    main()
