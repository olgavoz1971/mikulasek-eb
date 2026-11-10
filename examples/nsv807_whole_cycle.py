"""Whole-cycle timing of the NSV 807 light curves.

One shift of the whole template per orbital cycle, no eclipse window. The points
of cycle n are those with phase in (n - 0.25, n + 0.75) of the fitted ephemeris,
so the primary at phase n and the secondary at n + 0.5 lie inside the cycle. The
baseline of each cycle is a polynomial (default degree 1). The method is
``time_cycles`` of ``batch_whole_cycle_time.py``. The template and the points are
the ones of the tiered fit (``nsv807_s30_tiered.py``), read through ``prepare`` of
``nsv807_timing.py``.

Point rule, the one of the locked timing protocol applied to a cycle: at least 4
points, with at least 2 on each side of the cycle centre. The primary and the
secondary of one cycle get the same time shift, so their times are not
independent. The flag is ``gap`` when a gap longer than 2 cadences lies inside the
cycle window, and ``ok`` otherwise.

Run from ``examples`` for one light curve::

    python nsv807_whole_cycle.py --run ../results/nsv807/s30_protocol

Writes ``whole_cycle/`` in the run folder: ``eclipses.csv``, the two ``.dat``
files (full JD, error in days, flag, columns separated by two spaces), a summary
figure per kind and ``fits/E_<n>.png`` for every timed cycle.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np

import mikulasek_eb as eb
from batch_whole_cycle_time import time_cycles
from nsv807_timing import KINDS, prepare
from window_scan_mik import DAY_S, window_is_covered, write_csv

logger = logging.getLogger(__name__)

MIN_POINTS = 4
MIN_SIDE = 2


def cycle_arrays(t_mjd, mag, solution, n: int, delta: float, order: int) -> tuple:
    """Points, template and baseline of one timed cycle, for drawing.

    The template is shifted by the fitted ``delta``. The baseline is a polynomial
    found by least squares from the residual of the points to the template.

    Args:
        t_mjd (ndarray): Times, MJD.
        mag (ndarray): Magnitudes.
        solution (ModelSolution): The template.
        n (int): Cycle number.
        delta (float): Fitted time shift of the cycle (days).
        order (int): Baseline degree.

    Returns:
        tuple: Boolean selection, phase of the selected points relative to the
        cycle centre, template at those points, the baseline there (mag), the
        baseline coefficients and the time (MJD) their abscissa is measured from.
    """
    ephem = solution.ephem
    model = solution.build_model()
    phase = eb.phase(t_mjd, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0))
    sel = np.abs(phase - (n + 0.25)) < 0.5
    shifted = phase[sel] - delta / ephem["P0"]
    template = np.asarray(model(shifted, solution.p), float)
    x = t_mjd[sel] - t_mjd[sel].mean()
    design = np.vander(x, order + 1, increasing=True)
    coef = np.linalg.lstsq(design, mag[sel] - template, rcond=None)[0]
    return (
        sel,
        phase[sel] - (n + 0.25),
        template,
        design @ coef,
        coef,
        t_mjd[sel].mean(),
    )


def write_dat(out: Path, stem: str, rows: list[dict], flags: dict) -> None:
    """Write one ``.dat`` file per kind: full JD, error in days, flag.

    Args:
        out (Path): Folder.
        stem (str): Light-curve name.
        rows (list[dict]): Rows of ``eclipses.csv``.
        flags (dict): Flag text by cycle number.
    """
    for kind in KINDS:
        with (out / f"{stem}_{kind}.dat").open("w") as handle:
            for r in rows:
                if r["kind"] == kind and r["T_jd"] != "":
                    n = int(np.floor(float(r["cycle"])))
                    err = float(r["sigma_s"]) / DAY_S
                    handle.write(f"{float(r['T_jd']):.8f}  {err:.8f}  {flags[n]}\n")


def draw_cycles(out, t_mjd, mag, sig, solution, deltas, flags, order, stem) -> None:
    """Draw every timed cycle: points, shifted template with baseline, residual.

    Args:
        out (Path): ``whole_cycle`` folder. Figures go to ``fits/``.
        t_mjd (ndarray): Times, MJD.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        solution (ModelSolution): The template.
        deltas (dict): Fitted shift (days) by cycle number.
        flags (dict): Flag text by cycle number.
        order (int): Baseline degree.
        stem (str): Light-curve name, for the title.
    """
    import matplotlib.pyplot as plt

    folder = out / "fits"
    folder.mkdir(parents=True, exist_ok=True)
    for n, delta in deltas.items():
        sel, phi, template, base, coef, t_ref = cycle_arrays(
            t_mjd, mag, solution, n, delta, order
        )
        fit = template + base
        ephem = solution.ephem
        model = solution.build_model()
        centre_mjd = eb.predicted_time(
            n + 0.25, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0)
        )
        grid = np.linspace(
            centre_mjd - 0.5 * ephem["P0"], centre_mjd + 0.5 * ephem["P0"], 1200
        )
        grid_phase = eb.phase(grid, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0))
        curve = np.asarray(model(grid_phase - delta / ephem["P0"], solution.p), float)
        curve = curve + np.vander(grid - t_ref, order + 1, increasing=True) @ coef
        grid_phi = grid_phase - (n + 0.25)
        fig, axes = plt.subplots(
            2, 1, figsize=(10, 7), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
        )
        axes[0].errorbar(
            phi, mag[sel], yerr=sig[sel], fmt=".", ms=4, color="0.4", label="data"
        )
        axes[0].plot(
            grid_phi, curve, color="red", lw=1.2, label="template and baseline"
        )
        axes[0].invert_yaxis()
        axes[0].set_ylabel("Magnitude")
        axes[0].legend()
        axes[0].set_title(
            f"{stem}\ncycle {n}  shift {delta * DAY_S:+.1f} s  n={int(sel.sum())}  {flags[n]}"
        )
        axes[1].errorbar(
            phi, 1e3 * (mag[sel] - fit), yerr=1e3 * sig[sel], fmt=".", ms=4, color="0.4"
        )
        axes[1].axhline(0.0, color="red", lw=1)
        axes[1].set_xlabel(
            "Phase from the cycle centre (primary at -0.25, secondary at +0.25)"
        )
        axes[1].set_ylabel("Residual (mmag)")
        fig.tight_layout()
        fig.savefig(folder / f"E_{n}.png", dpi=100)
        plt.close(fig)


def draw_summary(out: Path, rows: list[dict], stem: str) -> None:
    """Draw the O-C of the timed cycles about a straight line, one figure per kind.

    Args:
        out (Path): ``whole_cycle`` folder.
        rows (list[dict]): Rows of ``eclipses.csv``.
        stem (str): Light-curve name, for the title.
    """
    import matplotlib.pyplot as plt
    from mikulasek_eb.timing import fit_ephemeris

    for kind in KINDS:
        timed = [r for r in rows if r["kind"] == kind and r["T_jd"] != ""]
        if len(timed) < 3:
            continue
        e = np.array([float(r["cycle"]) for r in timed])
        t = np.array([float(r["T_jd"]) for r in timed])
        s = np.array([float(r["sigma_s"]) for r in timed]) / DAY_S
        _, res, _ = fit_ephemeris(t, s, e, 1)
        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.errorbar(t, DAY_S * res, yerr=DAY_S * s, fmt="o", ms=3, color="tab:blue")
        ax.axhline(0.0, color="grey", lw=0.8)
        ax.set_xlabel("Julian Date")
        ax.set_ylabel("O-C about a straight line (s)")
        ax.set_title(f"{stem}\nwhole-cycle timing, {kind}, {len(timed)} cycles")
        fig.tight_layout()
        fig.savefig(out / f"whole_cycle_{kind}.png", dpi=130)
        plt.close(fig)


def run(run_dir: Path, name: str, flux: str, flux_err, skip: bool, order: int) -> None:
    """Time one light curve cycle by cycle and write tables and figures.

    Args:
        run_dir (Path): Tiered-fit run folder.
        name (str): Light-curve file name.
        flux (str): Flux column name.
        flux_err (str | None): Error column name, None for the constant sigma.
        skip (bool): The tiered fit skipped the outlier step.
        order (int): Baseline degree of each cycle.
    """
    payload, t_mjd, solution, cadence_d = prepare(run_dir, name, flux, flux_err, skip)
    stem = Path(name).stem
    out = run_dir / "whole_cycle"
    out.mkdir(parents=True, exist_ok=True)
    mag, sig = payload["mag"], payload["sig"]
    rows, _, deltas = time_cycles(
        t_mjd,
        mag,
        sig,
        solution,
        float(payload["epoch_jd"]),
        float(payload["period"]),
        order,
        MIN_POINTS,
        MIN_SIDE,
        stem,
    )
    t_jd = np.asarray(eb.mjd_to_jd(t_mjd), float)
    ephem = solution.ephem
    flags = {}
    for n in deltas:
        centre_mjd = eb.predicted_time(
            n + 0.25, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0)
        )
        covered = window_is_covered(
            t_jd, float(eb.mjd_to_jd(centre_mjd)), 0.5 * ephem["P0"], cadence_d
        )
        flags[n] = "ok" if covered else "gap"
    write_csv(out / "eclipses.csv", rows)
    write_dat(out, stem, rows, flags)
    draw_summary(out, rows, stem)
    draw_cycles(out, t_mjd, mag, sig, solution, deltas, flags, order, stem)
    logger.info(
        "%s: %d cycles timed, %d eclipse rows without a time, baseline degree %d",
        stem,
        len(deltas),
        sum(r["T_jd"] == "" for r in rows),
        order,
    )


def main() -> None:
    """Parse the command line and time one light curve."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", required=True, help="tiered-fit run folder")
    parser.add_argument("--file", default="TIC35119266_s30_SPOC_120s.vot")
    parser.add_argument("--flux", default="sap_flux")
    parser.add_argument("--flux-err", default="sap_flux_err")
    parser.add_argument("--no-flux-err", action="store_true")
    parser.add_argument("--skip-outliers", action="store_true")
    parser.add_argument(
        "--order", type=int, default=1, help="baseline degree per cycle"
    )
    args = parser.parse_args()
    flux_err = None if args.no_flux_err else args.flux_err
    run(
        Path(args.run).resolve(),
        args.file,
        args.flux,
        flux_err,
        args.skip_outliers,
        args.order,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
