"""Compare the two relative-error rules on the 120 s and 600 s spot curves.

``all`` is the core default: both components must be well determined. ``core`` asks that
only of the narrower component. The width limit still applies to both after adoption.
The default in the package is not changed.

Run from the repository root::

    MPLBACKEND=Agg python experiments/two_component_gate/compare.py
"""

from __future__ import annotations

import csv
import dataclasses
import logging
from pathlib import Path

import mikulasek_eb as eb

logger = logging.getLogger(__name__)
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFEST = GENERATED / "results" / "n01" / "manifest_spotsweep.csv"
OUT = Path(__file__).resolve().parent / "out" / "compare_120_600.csv"
CADENCES = (120, 600)
SCOPES = ("all", "core")


def config_for(payload: dict, t_mjd, scope: str) -> eb.Config:
    """The degraded-spot structure, with one relative-error rule.

    Args:
        payload (dict): Loaded light curve.
        t_mjd (ndarray): Times (MJD).
        scope (str): ``all`` or ``core``.

    Returns:
        Config: Configuration for ``fit``.
    """
    features = eb.extract_features(
        t_mjd, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    data = dataclasses.replace(
        eb.propose_candidates(features)[0].config.data, baseline_order=1
    )
    return eb.Config(
        eb.ModelSpec(profile=1, harmonics=(2,), ocon=True),
        data,
        eb.FitPolicy(allow_two_component=True, two_comp_rel_err_scope=scope),
        eb.TimingRecipe(window=0.12),
    ).validate()


def main() -> None:
    """Fit each curve under both rules and write the comparison table."""
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")
    with open(MANIFEST) as fh:
        rows = [
            r
            for r in csv.DictReader(fh)
            if int(r["cadence_s"]) in CADENCES and int(r["fail"]) == 0
        ]
    table = []
    for row in rows:
        path = GENERATED / row["path"]
        payload = eb.load_tess_votable(path)
        t_mjd = eb.absolute_jd_to_mjd(payload["t"])
        for scope in SCOPES:
            out = eb.fit(
                t_mjd,
                payload["mag"],
                payload["sig"],
                payload["period"],
                payload["epoch_jd"],
                config_for(payload, t_mjd, scope),
                data_id=path.stem,
            )
            stats = out.solution.stats
            trial = stats["two_component_trial"]
            table.append(
                dict(
                    physics=row["physics_id"],
                    cadence_s=row["cadence_s"],
                    scope=scope,
                    adopted=trial["adopted"],
                    status=stats["status"],
                    rms_mmag=round(1e3 * stats["rms_mag"], 2),
                    single_rms_mmag=round(trial["single_rms_mmag"], 2),
                    trial_rms_mmag=(
                        "" if trial["rms_mmag"] is None else round(trial["rms_mmag"], 2)
                    ),
                    failed="; ".join(trial["failed_gates"]),
                    reasons="; ".join(stats["status_reasons"]),
                )
            )
            logger.warning(
                "%s %ss %s: adopted=%s status=%s rms=%.2f",
                row["physics_id"],
                row["cadence_s"],
                scope,
                trial["adopted"],
                stats["status"],
                1e3 * stats["rms_mag"],
            )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    logger.warning("wrote %s", OUT)


if __name__ == "__main__":
    main()
