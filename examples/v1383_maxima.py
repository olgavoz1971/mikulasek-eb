"""Time the maxima of V1383 Tau with the saved Mikulášek template.

The primary maximum is the brightest phase of that sector's template between
phase 0 and phase 0.5, and the secondary maximum is the brightest phase between
phase 0.5 and phase 1. The half-windows are 0.15 cycles for the primary maximum
and 0.14 cycles for the secondary maximum. One shift of the stored template and
a straight-line baseline are fitted, with the same point rule as the minimum
timing. A gap longer than two cadences is flagged and kept. The eclipse
contact-wall check is not applied.

The points are the background-cleaned light curves used for the minimum timing,
and the template is that sector's chosen solution. The SAP errors are the
weights.

Run from ``mikulasek_eb`` with ``MPLBACKEND=Agg``::

    python examples/v1383_maxima.py

Writes ``results/v1383_tau/<sector>_fit_bkg/maxima/``: ``eclipses.csv``, two
``.dat`` files, ``result.json``, a summary figure per kind and one figure per
timed maximum.
"""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path

import numpy as np

import mikulasek_eb as eb
import nsv807_maxima as maxima
import nsv807_s30_tiered as tiered

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "v1383_tau"
CLEANED = RESULTS / "cleaned_bkg"
HALVES = {"primary_maximum": 0.15, "secondary_maximum": 0.14}
RUNS = (
    ("s42_fit_bkg", "TIC388995726_s42_TESS-SPOC_600s.vot"),
    ("s43_fit_bkg", "TIC388995726_s43_TESS-SPOC_600s.vot"),
    ("s44_fit_bkg", "TIC388995726_s44_TESS-SPOC_600s.vot"),
    ("s70_fit_bkg", "TIC388995726_s70_TESS-SPOC_200s.vot"),
    ("s71_fit_bkg", "TIC388995726_s71_TESS-SPOC_200s.vot"),
)


def windows_are_separate(phases: dict) -> None:
    """Require the two maximum windows to leave a gap on both sides.

    Args:
        phases (dict): Template phase of each maximum.

    Raises:
        RuntimeError: If the windows meet or overlap.
    """
    primary = phases["primary_maximum"]
    secondary = phases["secondary_maximum"]
    forward = (secondary - HALVES["secondary_maximum"]) - (
        primary + HALVES["primary_maximum"]
    )
    backward = (primary + 1.0 - HALVES["primary_maximum"]) - (
        secondary + HALVES["secondary_maximum"]
    )
    if forward <= 0.0 or backward <= 0.0:
        raise RuntimeError(
            f"the maximum windows meet: primary {primary:.5f}, "
            f"secondary {secondary:.5f}"
        )


def load_curve(run: Path, name: str):
    """Load one background-cleaned light curve and its chosen solution.

    Args:
        run (Path): Tiered-fit folder.
        name (str): VOTable file name in ``cleaned_bkg``.

    Returns:
        tuple: Payload, times (MJD), solution and cadence (days).

    Raises:
        FileNotFoundError: If the chosen folder does not hold one solution.
        RuntimeError: If the solution did not finish, or the times are not
        strictly increasing.
    """
    tiered.FLUX_ERR = "sap_flux_err"
    payload = tiered.load_sap(CLEANED / name)
    found = sorted((run / "chosen").glob("*.json"))
    if len(found) != 1:
        raise FileNotFoundError(f"expected one solution in {run / 'chosen'}: {found}")
    solution = eb.ModelSolution.load(found[0])
    if solution.stats["status"] != "ok":
        raise RuntimeError(f"{found[0].name}: status {solution.stats['status']}")
    t_mjd = np.asarray(eb.absolute_jd_to_mjd(payload["t"]), float)
    if np.any(np.diff(t_mjd) <= 0):
        raise RuntimeError(f"{name}: times are not strictly increasing")
    return payload, t_mjd, solution, float(np.median(np.diff(t_mjd)))


def run_one(folder: str, name: str) -> None:
    """Time both maxima of one sector and write its figures.

    Args:
        folder (str): Run folder name under ``results/v1383_tau``.
        name (str): VOTable file name.
    """
    run = RESULTS / folder
    payload, t_mjd, solution, cadence = load_curve(run, name)
    ephem = dict(solution.ephemeris)
    model = solution.build_model()
    params = dict(solution.params)
    data = {
        "t_jd": np.asarray(payload["t"], float),
        "t_mjd": t_mjd,
        "mag": np.asarray(payload["mag"], float),
        "sig": np.asarray(payload["sig"], float),
        "th": np.asarray(
            eb.phase(t_mjd, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0)), float
        ),
        "model": model,
        "params": params,
        "ephem": ephem,
        "period": float(ephem["P0"]),
        "cadence": cadence,
    }
    phases = {
        kind: maxima.brightest_phase(model, params, lo, hi)
        for kind, lo, hi in maxima.KINDS
    }
    windows_are_separate(phases)
    out = run / "maxima"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = []
    stem = Path(name).stem
    for kind, _, _ in maxima.KINDS:
        half = HALVES[kind]
        centre = phases[kind]
        logger.info(
            "%s %s: template phase %.5f, half-window %.2f",
            folder,
            kind,
            centre,
            half,
        )
        kind_rows = maxima.time_kind(data, kind, centre, half)
        rows.extend(kind_rows)
        timed = [row for row in kind_rows if row["T_jd"] != ""]
        logger.info(
            "%s %s: %d timed, %d without a time",
            folder,
            kind,
            len(timed),
            len(kind_rows) - len(timed),
        )
        for row in timed:
            maxima.draw_event(out, data, row, stem, half)
        maxima.draw_summary(out, data, kind_rows, stem, half)
    maxima.write_tables(out, stem, rows)
    saved = {
        "file": name,
        "half_cycles": HALVES,
        "template_phase": phases,
        "sigma": "sap_flux_err",
        "n_timed": {
            kind: sum(row["kind"] == kind and row["T_jd"] != "" for row in rows)
            for kind, _, _ in maxima.KINDS
        },
        "n_rejected": {
            kind: sum(row["kind"] == kind and row["T_jd"] == "" for row in rows)
            for kind, _, _ in maxima.KINDS
        },
    }
    (out / "result.json").write_text(json.dumps(saved, indent=2) + "\n")
    logger.info("wrote %s", out)


def main() -> None:
    """Time the maxima of the five V1383 Tau light curves."""
    for folder, name in RUNS:
        run_one(folder, name)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
