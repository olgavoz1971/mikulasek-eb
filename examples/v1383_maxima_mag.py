"""Time the maxima of V1383 Tau and record the magnitude of the fitted template.

The windows, the template and the time fit are those of ``v1383_maxima.py``.
Each timed maximum also gets the magnitude of the fitted model at the fitted
time: the saved template at its brightest phase, plus the straight-line
baseline of that window. The uncertainty is the formal error of that
magnitude from the same weighted fit, rescaled as the time error is. It is
not the brightest observed point.

The files already written by ``v1383_maxima.py`` are not changed.

Run from ``mikulasek_eb``::

    python examples/v1383_maxima_mag.py

Writes ``results/v1383_tau/<sector>_fit_bkg/maxima_mag/`` and the collected
files in ``results/v1383_tau/timing_max_mag_dat/``. Each secondary is written
twice there, the second with the times reduced by half the file period. The
magnitudes are the same in both.
"""

from __future__ import annotations

import csv
import logging
import shutil
from decimal import Decimal
from pathlib import Path

import numpy as np

import mikulasek_eb as eb
import nsv807_collect_maxima as files
import nsv807_maxima as maxima
import v1383_collect_dat as minima
import v1383_collect_maxima as collected
import v1383_maxima as base
from nsv807_collect_dat import header
from window_scan_mik import DAY_S

logger = logging.getLogger(__name__)

OUT = base.RESULTS / "timing_max_mag_dat"
SIGMA_AGREE = 1e-3
CHI2_AGREE = 1e-6
PHASE_STEP = 1e-6


def template_peak(model, params: dict, centre: float) -> float:
    """Magnitude of the saved template at one maximum.

    Args:
        model: Template model.
        params (dict): Template parameters.
        centre (float): Template phase of the maximum.

    Returns:
        float: Template magnitude, relative to its out-of-eclipse level.

    Raises:
        RuntimeError: If shifting the phase by whole cycles changes the value.
    """
    peak = float(np.asarray(model(centre, params), float))
    for cycle in (-2, -1, 1, 2):
        shifted = float(np.asarray(model(centre + cycle, params), float))
        if abs(shifted - peak) > 1e-8:
            raise RuntimeError(
                f"the template at phase {centre} changes by a whole cycle"
            )
    return peak


def maximum_magnitude(
    data: dict, row: dict, half: float, peak: float
) -> tuple[float, float]:
    """Magnitude of the fitted template at one timed maximum, and its error.

    Args:
        data (dict): Light curve and template.
        row (dict): Timed event from ``time_kind``.
        half (float): Half-window in cycles.
        peak (float): Template magnitude at the maximum phase.

    Returns:
        tuple[float, float]: Magnitude and its formal 1-sigma error, both in
        magnitudes.

    Raises:
        RuntimeError: If the repeated fit disagrees with the saved time, or if
        the magnitude error is not the error of that same fit.
    """
    event = float(row["cycle"]) + float(row["template_phase"])
    ephem = data["ephem"]
    fit = eb.fit_minimum(
        data["t_mjd"],
        data["mag"],
        data["sig"],
        data["model"],
        data["params"],
        ephem,
        event,
        window=half,
        order=1,
        free=(),
    )
    t_jd = float(eb.mjd_to_jd(fit.T))
    if (
        f"{t_jd:.6f}" != row["T_jd"]
        or f"{float(fit.sigma) * DAY_S:.3f}" != row["sigma_s"]
    ):
        raise RuntimeError(
            f"{row['kind']} E={row['cycle']}: the repeated fit does not match the saved time"
        )
    if int(fit.n) != int(row["n_points"]):
        raise RuntimeError(f"{row['kind']} E={row['cycle']}: point count changed")
    sel = np.abs(data["th"] - event) < half
    if int(sel.sum()) != int(fit.n):
        raise RuntimeError(f"{row['kind']} E={row['cycle']}: the window does not match")
    ts = data["t_mjd"][sel]
    ys = data["mag"][sel]
    ths = data["th"][sel]
    w = 1.0 / data["sig"][sel]
    period = ephem["P0"] * (1.0 + ephem.get("Pdot", 0.0) * event)
    delta = float(fit.delta)
    template = np.asarray(data["model"](ths - delta / period, data["params"]), float)
    spread = float(ts.std())
    if spread <= 0.0:
        raise RuntimeError(
            f"{row['kind']} E={row['cycle']}: the window has no time span"
        )
    scaled = (ts - float(ts.mean())) / spread
    design = np.column_stack([np.ones(len(ts)), scaled])
    coef = np.linalg.lstsq(design * w[:, None], w * (ys - template), rcond=None)[0]
    resid = w * (ys - template - design @ coef)
    dof = int(fit.n) - 3
    rchi2 = float(resid @ resid) / dof
    if abs(rchi2 - float(fit.rchi2)) > CHI2_AGREE * max(float(fit.rchi2), 1e-8):
        raise RuntimeError(
            f"{row['kind']} E={row['cycle']}: chi-square {rchi2} against the fit {fit.rchi2}"
        )
    step = PHASE_STEP * period
    above = np.asarray(
        data["model"](ths - (delta + step) / period, data["params"]), float
    )
    below = np.asarray(
        data["model"](ths - (delta - step) / period, data["params"]), float
    )
    slope = (above - below) / (2.0 * step)
    jac = np.column_stack([-w * slope, -w * design[:, 0], -w * design[:, 1]])
    scale = np.linalg.svd(jac, compute_uv=False)
    if scale[-1] <= 1e-8 * scale[0]:
        raise RuntimeError(
            f"{row['kind']} E={row['cycle']}: the fit for the magnitude is rank deficient"
        )
    cov = float(fit.rchi2) * np.linalg.inv(jac.T @ jac)
    sigma_days = float(np.sqrt(cov[0, 0]))
    if abs(sigma_days - float(fit.sigma)) > SIGMA_AGREE * float(fit.sigma):
        raise RuntimeError(
            f"{row['kind']} E={row['cycle']}: time error {sigma_days} against the fit {fit.sigma}"
        )
    at = (float(fit.T) - float(ts.mean())) / spread
    mag = peak + float(coef[0] + coef[1] * at)
    gradient = np.array([float(coef[1]) / spread, 1.0, at])
    err = float(np.sqrt(gradient @ cov @ gradient))
    if not np.isfinite(mag) or not np.isfinite(err) or err <= 0.0:
        raise RuntimeError(
            f"{row['kind']} E={row['cycle']}: magnitude {mag}, error {err}"
        )
    return mag, err


def write_sector(folder: Path, stem: str, rows: list[dict]) -> None:
    """Write one sector table and its two ``.dat`` files.

    Args:
        folder (Path): ``maxima_mag`` folder.
        stem (str): Light-curve name.
        rows (list[dict]): Every event, with magnitudes on the timed rows.
    """
    fields = [
        "kind",
        "cycle",
        "template_phase",
        "T_jd",
        "sigma_s",
        "mag",
        "err_mag",
        "half_cycles",
        "n_points",
        "dof",
        "gap_in_window",
        "reason",
    ]
    with (folder / "eclipses.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    for kind, _, _ in maxima.KINDS:
        phase = float(
            next(row["template_phase"] for row in rows if row["kind"] == kind)
        )
        half = float(next(row["half_cycles"] for row in rows if row["kind"] == kind))
        with (folder / f"{stem}_{kind}.dat").open("w") as handle:
            handle.write(f"# Maximum: {kind}\n")
            handle.write(
                f"# Window centre: template phase {phase:.6f}, half-window {half:.2f} cycles\n"
            )
            handle.write(
                "# Magnitude: fitted template at the fitted time "
                "(saved template plus the straight-line baseline)\n"
            )
            handle.write("# Columns: jd  err_days  flag  mag  err_mag\n")
            for row in rows:
                if row["kind"] != kind or row["T_jd"] == "":
                    continue
                err = float(row["sigma_s"]) / DAY_S
                handle.write(
                    f"{float(row['T_jd']):.8f}  {err:.8f}  {maxima.flag_text(row)}  "
                    f"{float(row['mag']):.6f}  {float(row['err_mag']):.6f}\n"
                )


def data_line(row: dict, jd: Decimal) -> str:
    """One data line, with the time given separately so a copy can shift it.

    Args:
        row (dict): Timed event.
        jd (Decimal): Time to write, days.

    Returns:
        str: ``jd  err_days  flag  mag  err_mag``.
    """
    err = Decimal(f"{float(row['sigma_s']) / DAY_S:.8f}")
    flag = "gap" if row["gap_in_window"] is True else "ok"
    mag = Decimal(f"{float(row['mag']):.6f}")
    err_mag = Decimal(f"{float(row['err_mag']):.6f}")
    return f"{jd:.8f}  {err:.8f}  {flag}  {mag:.6f}  {err_mag:.6f}"


def write_collected(
    path: Path, head: list[str], rows: list[dict], half: Decimal | None
) -> int:
    """Write one collected file.

    Args:
        path (Path): Destination.
        head (list[str]): Plain header lines.
        rows (list[dict]): Events of one kind.
        half (Decimal or None): Subtract this from each time. None leaves the
            measured times.

    Returns:
        int: Number of rows that hold a time.
    """
    timed = [row for row in rows if row["T_jd"]]
    missed = [row for row in rows if not row["T_jd"]]
    lines = []
    if half is not None:
        lines.append(
            f"# Times in this copy are the source times minus {half} d, half the file period."
        )
        lines.append("# The magnitudes are the magnitudes at the unshifted times.")
    lines += header(head)
    lines.append(
        f"# rows with a time: {len(timed)}, rows without a time: {len(missed)}"
    )
    lines += [f"# no time  E={row['cycle']} {row['reason']}" for row in missed]
    for row in timed:
        jd = Decimal(f"{float(row['T_jd']):.8f}")
        if half is not None:
            jd -= half
        lines.append(data_line(row, jd))
    path.write_text("\n".join(lines) + "\n")
    return len(timed)


def run_one(folder: str, name: str) -> dict:
    """Time one sector and record the magnitude of each maximum.

    Args:
        folder (str): Run folder name under ``results/v1383_tau``.
        name (str): VOTable file name.

    Returns:
        dict: Rows of both kinds and the names used in the collected header.

    Raises:
        RuntimeError: If this run does not reproduce the saved times.
    """
    run = base.RESULTS / folder
    saved_rows = list(csv.DictReader((run / "maxima" / "eclipses.csv").open()))
    payload, t_mjd, solution, cadence = base.load_curve(run, name)
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
    base.windows_are_separate(phases)
    rows = []
    for kind, _, _ in maxima.KINDS:
        half = base.HALVES[kind]
        peak = template_peak(model, params, phases[kind])
        kind_rows = maxima.time_kind(data, kind, phases[kind], half)
        previous = [row for row in saved_rows if row["kind"] == kind]
        if len(kind_rows) != len(previous):
            raise RuntimeError(
                f"{folder} {kind}: {len(kind_rows)} events, saved {len(previous)}"
            )
        for row, old in zip(kind_rows, previous):
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
            mag, err = maximum_magnitude(data, row, half, peak)
            row["mag"] = f"{mag:.6f}"
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
    write_sector(out, Path(name).stem, rows)
    return {
        "rows": rows,
        "solution": next((run / "chosen").glob("*.json")).stem.split("__")[1],
        "stem": Path(name).stem,
    }


def collect(sectors: dict) -> None:
    """Write the shared folder, including the shifted secondary copies.

    Args:
        sectors (dict): Output of ``run_one``, keyed by the sector label.
    """
    half = collected.file_half_period()
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    with open(base.RESULTS / "cleaned_bkg" / "summary.csv", newline="") as handle:
        summaries = {row["file"]: row for row in csv.DictReader(handle)}
    index = []
    for lc, (stem, label) in minima.LIGHT_CURVES.items():
        sector = sectors[lc]
        summary = summaries[f"{stem}.vot"]
        for kind in files.KINDS:
            rows = [row for row in sector["rows"] if row["kind"] == kind]
            phase = float(rows[0]["template_phase"])
            window = float(rows[0]["half_cycles"])
            label_kind = files.LABEL[kind]
            head = [
                "Star: V1383 Tau (TIC 388995726), contact binary with a strong O'Connell effect",
                f"Light curve: {stem}.vot ({label}), SAP flux",
                f"Maximum: {label_kind}",
                "Method: Mikulasek phenomenological template, window timing of the maximum "
                "(shift of the saved template, linear baseline)",
                f"Window: centre at template phase {phase:.6f}, the brightest point of the "
                f"template, half-window {window:.2f} cycles",
                "Cleaning: points with sap_bkg above "
                f"{summary['ratio']} x the median background removed: "
                f"{summary['removed']} of {summary['rows_in']} rows",
                f"Template: tiered fit {sector['solution']}, baseline split at gaps over 0.3 d",
                "Sigma: sap_flux_err of the light curve",
                "Flag: ok, gap (gap longer than 2 cadences in the window). "
                "The eclipse contact-wall flag is not used",
                "Magnitude: the fitted template at the fitted time, saved template plus the "
                "straight-line baseline of the window, in magnitudes. Not the brightest observed point",
                "Time: full JD of the time system of the source file",
                "Columns: jd  err_days  flag  mag  err_mag",
                "err_days is the formal 1 sigma error of the time in days",
                "err_mag is the formal 1 sigma error of the magnitude",
                "Lines starting with # are comments. Rows without a time are comments too.",
            ]
            name = f"V1383Tau_{lc}_mikulasek_maximum_{label_kind}.dat"
            n_times = write_collected(OUT / name, head, rows, None)
            index.append((name, lc, "mikulasek_maximum", label_kind, n_times))
            if kind != "secondary_maximum":
                continue
            shifted = f"V1383Tau_{lc}_mikulasek_maximum_secondary_minus_half_period.dat"
            n_shift = write_collected(OUT / shifted, head, rows, half)
            index.append(
                (
                    shifted,
                    lc,
                    "mikulasek_maximum",
                    "secondary_minus_half_period",
                    n_shift,
                )
            )
    with open(OUT / "index.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "light_curve", "method", "kind", "n_times"])
        writer.writerows(index)
    logger.info("wrote %d files to %s", len(index), OUT)


def main() -> None:
    """Time the five sectors and collect the files."""
    sectors = {}
    for folder, name in base.RUNS:
        lc = folder.split("_")[0]
        sectors[lc] = run_one(folder, name)
    collect(sectors)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
