"""Record the magnitude of each timed maximum of NSV 807.

The windows, the template and the time fit are those of ``nsv807_maxima.py``.
The half-window is 0.16 cycles for both maxima. Each timed maximum gets the
magnitude of the fitted model at the fitted time: the saved template at its
brightest phase, plus the straight-line baseline of that window. The
uncertainty is the formal error of that magnitude from the same weighted fit.
It is not the brightest observed point.

The files already written by ``nsv807_maxima.py`` are not changed.

Run from ``mikulasek_eb``::

    python examples/nsv807_maxima_mag.py

Writes ``results/nsv807/<run>/maxima_mag/``.
"""

from __future__ import annotations

import csv
import logging
import shutil
from pathlib import Path

import numpy as np

import mikulasek_eb as eb
import nsv807_maxima as maxima
import v1383_maxima_mag as mag

logger = logging.getLogger(__name__)


def run_one(folder: str, name: str, flux_err: str | None) -> None:
    """Record the magnitude of each saved maximum of one sector.

    Args:
        folder (str): Run folder name under ``results/nsv807``.
        name (str): VOTable file name.
        flux_err (str | None): SAP error column, or None for the constant sigma.

    Raises:
        RuntimeError: If this run does not reproduce the saved timing.
    """
    run = maxima.RESULTS / folder
    saved_rows = list(csv.DictReader((run / "maxima" / "eclipses.csv").open()))
    payload, t_mjd, solution, cadence = maxima.load_curve(run, name, flux_err)
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
    rows = []
    for kind, _, _ in maxima.KINDS:
        half = maxima.HALF
        peak = mag.template_peak(model, params, phases[kind])
        kind_rows = maxima.time_kind(data, kind, phases[kind], half)
        previous = [row for row in saved_rows if row["kind"] == kind]
        if len(kind_rows) != len(previous):
            raise RuntimeError(
                f"{folder} {kind}: {len(kind_rows)} events, saved {len(previous)}"
            )
        for row, old in zip(kind_rows, previous):
            if f"{phases[kind]:.6f}" != f"{float(old['template_phase']):.6f}":
                raise RuntimeError(
                    f"{folder} {kind}: the window centre differs from the saved timing"
                )
            if float(old["half_cycles"]) != half:
                raise RuntimeError(
                    f"{folder} {kind} E={row['cycle']}: half-window {old['half_cycles']}"
                )
            if int(row["cycle"]) != int(old["cycle"]) or row["T_jd"] != old["T_jd"]:
                raise RuntimeError(
                    f"{folder} {kind} E={row['cycle']}: time differs from the saved timing"
                )
            if (
                row["reason"] != old["reason"]
                or str(row["gap_in_window"]) != old["gap_in_window"]
            ):
                raise RuntimeError(
                    f"{folder} {kind} E={row['cycle']}: flag differs from the saved timing"
                )
            if row["T_jd"] == "":
                row["mag"] = ""
                row["err_mag"] = ""
                continue
            value, err = mag.maximum_magnitude(data, row, half, peak)
            row["mag"] = f"{value:.6f}"
            row["err_mag"] = f"{err:.6f}"
        rows.extend(kind_rows)
        timed = [row for row in kind_rows if row["T_jd"]]
        logger.info(
            "%s %s: %d times, median magnitude %.4f, median error %.5f mag",
            folder,
            kind,
            len(timed),
            float(np.median([float(row["mag"]) for row in timed])),
            float(np.median([float(row["err_mag"]) for row in timed])),
        )
    out = run / "maxima_mag"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    mag.write_sector(out, Path(name).stem, rows)


def main() -> None:
    """Record the magnitudes of the four NSV 807 light curves."""
    for folder, name, flux_err in maxima.RUNS:
        run_one(folder, name, flux_err)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
