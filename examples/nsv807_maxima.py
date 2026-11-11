"""Time the maxima of NSV 807 with the saved Mikulášek template.

The primary maximum is the brightest phase of that sector's template between
phase 0 and phase 0.5, and the secondary maximum is the brightest phase between
phase 0.5 and phase 1. The half-window is 0.16 cycles about that phase. One
shift of the stored template and a straight-line baseline are fitted, with the
same point rule as the minimum timing: at least 4 points and at least 2 on each
side of the centre. A gap longer than two cadences is flagged and kept. The
eclipse contact-wall check is not applied, because the window is on the maximum.

The points are the cleaned points of the minimum timing. A file with no SAP
error column keeps the constant sigma implied by the chi-square of its saved
solution.

Run from ``mikulasek_eb`` with ``MPLBACKEND=Agg``::

    python examples/nsv807_maxima.py

Writes ``results/nsv807/<run>/maxima/``: ``eclipses.csv``, two ``.dat`` files,
``result.json``, a summary figure per kind and one figure per timed maximum.
"""

from __future__ import annotations

import csv
import json
import logging
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import mikulasek_eb as eb
import nsv807_s30_tiered as tiered
from window_scan_mik import DAY_S, window_is_covered

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "nsv807"
HALF = 0.16
MIN_POINTS = 4
MIN_EACH_SIDE = 2
MATCH_DAYS = 1e-5
PHASE_GRID = 200001
SECOND_PEAK_CYCLES = 0.05
SECOND_PEAK_MAG = 1e-3
KINDS = (
    ("primary_maximum", 0.0, 0.5),
    ("secondary_maximum", 0.5, 1.0),
)
RUNS = (
    ("s30_protocol", "TIC35119266_s30_SPOC_120s.vot", "sap_flux_err"),
    ("s30_600s_protocol", "TIC35119266_s30_TESS-SPOC_600s.vot", "sap_flux_err"),
    ("s4_protocol", "TIC35119266_s4_TESS-SPOC_1800s.vot", "sap_flux_err"),
    ("s97_protocol", "TIC35119266_s97_QLP_200s.vot", None),
)


def brightest_phase(model, params: dict, lo: float, hi: float) -> float:
    """Phase of the brightest point of the template on an open interval.

    The template is a magnitude, so the brightest point is its lowest value.

    Args:
        model: Template model.
        params (dict): Template parameters.
        lo (float): Left edge, excluded, cycles.
        hi (float): Right edge, excluded, cycles.

    Returns:
        float: Phase in cycles.

    Raises:
        RuntimeError: If a second peak within 1 mmag lies 0.05 cycles away.
    """
    phi = np.linspace(lo, hi, PHASE_GRID)
    phi = phi[(phi > lo) & (phi < hi)]
    curve = np.asarray(model(phi, params), float)
    at = int(np.argmin(curve))
    other = (np.abs(phi - phi[at]) > SECOND_PEAK_CYCLES) & (
        curve <= curve[at] + SECOND_PEAK_MAG
    )
    if np.any(other):
        raise RuntimeError(
            f"two peaks between {lo} and {hi}: {phi[at]:.5f} and "
            f"{phi[other][int(np.argmin(curve[other]))]:.5f}"
        )
    return float(phi[at])


def solution_sigma(solution) -> float:
    """Constant sigma implied by the chi-square and the rms of a solution.

    Args:
        solution (ModelSolution): Saved fit whose weights were one constant.

    Returns:
        float: Sigma in magnitudes.

    Raises:
        RuntimeError: If the saved chi-square cannot give a sigma.
    """
    stats = solution.stats
    chi2 = float(stats["chi2"])
    n_points = float(stats["n_points"])
    if chi2 <= 0.0 or n_points <= 0.0:
        raise RuntimeError("the saved solution has no chi-square for a constant sigma")
    return float(stats["rms_mag"]) * np.sqrt(n_points / chi2)


def drop_listed(payload: dict, removed: Path) -> dict:
    """Drop the points listed by the minimum-timing cleaning.

    Args:
        payload (dict): Light curve from ``load_sap``.
        removed (Path): ``__removed.csv`` of the tiered fit. A missing file means
            that run removed nothing.

    Returns:
        dict: Light curve without the listed times.

    Raises:
        RuntimeError: If a listed time does not match one row.
    """
    if not removed.exists():
        return payload
    with removed.open(newline="") as handle:
        listed = [float(row["t_jd"]) for row in csv.DictReader(handle)]
    if not listed:
        return payload
    drop = []
    for time in listed:
        hit = np.flatnonzero(np.abs(payload["t"] - time) < MATCH_DAYS)
        if hit.size != 1:
            raise RuntimeError(f"{removed.name}: time {time} matches {hit.size} rows")
        drop.append(int(hit[0]))
    if len(set(drop)) != len(drop):
        raise RuntimeError(f"{removed.name}: a time is listed twice")
    keep = np.ones(len(payload["t"]), dtype=bool)
    keep[drop] = False
    out = dict(payload)
    for key in ("t", "mag", "sig"):
        out[key] = payload[key][keep]
    logger.info("%s: dropped %d listed points", removed.name, len(drop))
    return out


def load_curve(run: Path, name: str, flux_err: str | None):
    """Load one cleaned light curve and its chosen solution.

    Args:
        run (Path): Tiered-fit folder.
        name (str): VOTable file name.
        flux_err (str | None): SAP error column. None selects the constant sigma
            of the saved solution.

    Returns:
        tuple: Payload, times (MJD), solution and cadence (days).
    """
    tiered.FLUX_ERR = flux_err
    path = tiered.DATA / name
    payload = drop_listed(tiered.load_sap(path), run / f"{path.stem}__removed.csv")
    found = sorted((run / "chosen").glob("*.json"))
    if len(found) != 1:
        raise FileNotFoundError(f"expected one solution in {run / 'chosen'}: {found}")
    solution = eb.ModelSolution.load(found[0])
    if solution.stats["status"] != "ok":
        raise RuntimeError(f"{found[0].name}: status {solution.stats['status']}")
    if flux_err is None:
        sigma = solution_sigma(solution)
        payload = dict(payload)
        payload["sig"] = np.full_like(payload["mag"], sigma)
        logger.info("%s constant sigma %.2f mmag", name, 1e3 * sigma)
    t_mjd = np.asarray(eb.absolute_jd_to_mjd(payload["t"]), float)
    if np.any(np.diff(t_mjd) <= 0):
        raise RuntimeError(f"{name}: times are not strictly increasing")
    return payload, t_mjd, solution, float(np.median(np.diff(t_mjd)))


def consider(th: np.ndarray, centre: float, half: float = HALF) -> np.ndarray:
    """Integer cycles whose window contains at least one point.

    Args:
        th (ndarray): Phase of the data on the fitted ephemeris.
        centre (float): Template phase of the maximum.
        half (float): Half-window in cycles. The default is the NSV 807 window.

    Returns:
        ndarray: Integers ``n`` for events at phase ``n + centre``.
    """
    counted = np.round(th - centre).astype(int)
    kept = []
    for n in np.unique(counted):
        if np.any(np.abs(th - (int(n) + centre)) < half):
            kept.append(int(n))
    return np.asarray(kept, int)


def time_kind(data: dict, kind: str, centre: float, half: float = HALF) -> list[dict]:
    """Time every maximum of one kind.

    Args:
        data (dict): Arrays and the solution for one light curve.
        kind (str): ``primary_maximum`` or ``secondary_maximum``.
        centre (float): Template phase of this maximum.
        half (float): Half-window in cycles. The default is the NSV 807 window.

    Returns:
        list[dict]: One row per event, timed or rejected.
    """
    th = data["th"]
    ephem = data["ephem"]
    rows = []
    for n in consider(th, centre):
        event = n + centre
        offset = th - event
        inside = np.abs(offset) < half
        n_in = int(inside.sum())
        before = int(np.sum(offset[inside] < 0))
        after = int(np.sum(offset[inside] > 0))
        base = {
            "kind": kind,
            "cycle": n,
            "template_phase": centre,
            "half_cycles": half,
            "T_jd": "",
            "sigma_s": "",
            "n_points": n_in,
            "dof": "",
            "gap_in_window": "",
            "reason": "",
        }
        if n_in < MIN_POINTS:
            base["reason"] = "too_few_points"
            rows.append(base)
            continue
        if before < MIN_EACH_SIDE or after < MIN_EACH_SIDE:
            base["reason"] = "one_side_uncovered"
            rows.append(base)
            continue
        try:
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
        except eb.FitDidNotConverge as ex:
            logger.warning("%s E=%s fit failed: %s", kind, event, ex)
            base["reason"] = "fit_failed"
            rows.append(base)
            continue
        t_jd = float(eb.mjd_to_jd(fit.T))
        base.update(
            {
                "T_jd": f"{t_jd:.6f}",
                "sigma_s": f"{float(fit.sigma) * DAY_S:.3f}",
                "n_points": int(fit.n),
                "dof": int(fit.n) - 3,
                "gap_in_window": not window_is_covered(
                    data["t_jd"], t_jd, half * data["period"], data["cadence"]
                ),
            }
        )
        rows.append(base)
    return rows


def flag_text(row: dict) -> str:
    """Flag of one timed maximum.

    Args:
        row (dict): A row that has a time.

    Returns:
        str: ``ok`` or ``gap``.
    """
    return "gap" if row["gap_in_window"] is True else "ok"


def write_tables(out: Path, stem: str, rows: list[dict]) -> None:
    """Write the table and one ``.dat`` file per kind.

    Args:
        out (Path): Output folder.
        stem (str): Light-curve name.
        rows (list[dict]): Every event.
    """
    fields = [
        "kind",
        "cycle",
        "template_phase",
        "T_jd",
        "sigma_s",
        "half_cycles",
        "n_points",
        "dof",
        "gap_in_window",
        "reason",
    ]
    with (out / "eclipses.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    for kind, _, _ in KINDS:
        phase = next(r["template_phase"] for r in rows if r["kind"] == kind)
        half = float(next(r["half_cycles"] for r in rows if r["kind"] == kind))
        with (out / f"{stem}_{kind}.dat").open("w") as handle:
            handle.write(f"# Maximum: {kind}\n")
            handle.write(
                f"# Window centre: template phase {phase:.6f}, half-window {half:.2f} cycles\n"
            )
            handle.write("# Columns: jd  err_days  flag\n")
            for row in rows:
                if row["kind"] != kind:
                    continue
                if row["T_jd"] == "":
                    handle.write(f"# E={row['cycle']} rejected: {row['reason']}\n")
                    continue
                err = float(row["sigma_s"]) / DAY_S
                handle.write(f"{float(row['T_jd']):.8f}  {err:.8f}  {flag_text(row)}\n")


def draw_event(out: Path, data: dict, row: dict, stem: str, half: float = HALF) -> None:
    """Draw one timed maximum with the shifted template.

    Args:
        out (Path): Output folder.
        data (dict): Light curve and template.
        row (dict): Timed event.
        stem (str): Light-curve name.
        half (float): Half-window in cycles. The default is the NSV 807 window.
    """
    centre = float(row["T_jd"])
    t_mjd = np.asarray(eb.absolute_jd_to_mjd(centre), float)
    event = row["cycle"] + row["template_phase"]
    delta = t_mjd - eb.predicted_time(
        event, data["ephem"]["M0"], data["ephem"]["P0"], data["ephem"].get("Pdot", 0.0)
    )
    period = data["ephem"]["P0"] * (1.0 + data["ephem"].get("Pdot", 0.0) * event)
    inside = np.abs(data["th"] - event) < half
    phase = data["th"][inside] - float(delta) / period
    template = np.asarray(data["model"](phase, data["params"]), float)
    elapsed = data["t_jd"][inside] - centre
    design = np.column_stack([np.ones(int(inside.sum())), elapsed])
    coef = np.linalg.lstsq(design, data["mag"][inside] - template, rcond=None)[0]
    grid = np.linspace(-half * data["period"], half * data["period"], 400)
    grid_phase = row["template_phase"] + grid / period
    curve = np.asarray(data["model"](grid_phase, data["params"]), float)
    folder = out / "fits" / row["kind"]
    folder.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        2, 1, figsize=(8, 7), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
    )
    axes[0].errorbar(
        elapsed,
        data["mag"][inside],
        yerr=data["sig"][inside],
        fmt=".",
        ms=4,
        color="0.4",
        label="data",
    )
    axes[0].plot(
        grid,
        curve + coef[0] + coef[1] * grid,
        color="red",
        lw=1.5,
        label="template and baseline",
    )
    axes[0].axvline(0.0, color="blue", ls=":", lw=1, label="fitted time")
    axes[0].invert_yaxis()
    axes[0].set_ylabel("Magnitude")
    axes[0].legend()
    axes[0].set_title(
        f"{stem}\n{row['kind'].replace('_', ' ')} after cycle {row['cycle']}, "
        f"template phase {row['template_phase']:.5f}\n"
        f"T={centre:.6f}  sigma {row['sigma_s']} s  n={row['n_points']}  {flag_text(row)}"
    )
    fit = template + design @ coef
    axes[1].errorbar(
        elapsed,
        1e3 * (data["mag"][inside] - fit),
        yerr=1e3 * data["sig"][inside],
        fmt=".",
        ms=4,
        color="0.4",
    )
    axes[1].axhline(0.0, color="red", lw=1)
    axes[1].set_xlabel("Time from the fitted maximum (d)")
    axes[1].set_ylabel("Residual (mmag)")
    fig.tight_layout()
    fig.savefig(folder / f"E_{row['cycle']}.png", dpi=100)
    plt.close(fig)


def draw_summary(
    out: Path, data: dict, rows: list[dict], stem: str, half: float = HALF
) -> None:
    """Draw every timed maximum of one kind on the shifted template.

    Args:
        out (Path): Output folder.
        data (dict): Light curve and template.
        rows (list[dict]): Events of one kind.
        stem (str): Light-curve name.
        half (float): Half-window in cycles. The default is the NSV 807 window.
    """
    timed = [r for r in rows if r["T_jd"] != ""]
    if not timed:
        return
    kind = timed[0]["kind"]
    centre_phase = float(timed[0]["template_phase"])
    phases, values = [], []
    for row in timed:
        centre = float(row["T_jd"])
        t_mjd = float(eb.absolute_jd_to_mjd(centre))
        event = row["cycle"] + centre_phase
        delta = t_mjd - eb.predicted_time(
            event,
            data["ephem"]["M0"],
            data["ephem"]["P0"],
            data["ephem"].get("Pdot", 0.0),
        )
        period = data["ephem"]["P0"] * (1.0 + data["ephem"].get("Pdot", 0.0) * event)
        inside = np.abs(data["th"] - event) < half
        phase = data["th"][inside] - event - delta / period
        template = np.asarray(
            data["model"](data["th"][inside] - delta / period, data["params"]), float
        )
        elapsed = data["t_jd"][inside] - centre
        design = np.column_stack([np.ones(int(inside.sum())), elapsed])
        coef = np.linalg.lstsq(design, data["mag"][inside] - template, rcond=None)[0]
        phases.append(phase)
        values.append(data["mag"][inside] - design @ coef)
    phi = np.concatenate(phases)
    mag = np.concatenate(values)
    grid = np.linspace(-half, half, 400)
    line = np.asarray(data["model"](centre_phase + grid, data["params"]), float)
    fig, axis = plt.subplots(figsize=(8, 4.5))
    axis.scatter(phi, mag, s=3, color="0.55", label="data, own baseline removed")
    axis.plot(grid, line, color="red", lw=1.5, label="template")
    axis.invert_yaxis()
    axis.set_xlabel("Phase from the template maximum")
    axis.set_ylabel("Magnitude")
    axis.legend(frameon=False)
    axis.set_title(
        f"{stem}\n{kind.replace('_', ' ')}, phase {centre_phase:.5f}, "
        f"half-window {half:.2f}, {len(timed)} timed"
    )
    fig.tight_layout()
    fig.savefig(out / f"{kind}.png", dpi=120)
    plt.close(fig)


def run_one(folder: str, name: str, flux_err: str | None) -> None:
    """Time both maxima of one sector and write its figures.

    Args:
        folder (str): Run folder name under ``results/nsv807``.
        name (str): VOTable file name.
        flux_err (str | None): SAP error column, or None for the constant sigma.
    """
    run = RESULTS / folder
    payload, t_mjd, solution, cadence = load_curve(run, name, flux_err)
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
    out = run / "maxima"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = []
    phases = {}
    for kind, lo, hi in KINDS:
        centre = brightest_phase(model, params, lo, hi)
        phases[kind] = centre
        logger.info(
            "%s %s: template phase %.5f, window %.5f to %.5f",
            folder,
            kind,
            centre,
            centre - HALF,
            centre + HALF,
        )
        kind_rows = time_kind(data, kind, centre)
        rows.extend(kind_rows)
        timed = [r for r in kind_rows if r["T_jd"] != ""]
        logger.info(
            "%s %s: %d timed, %d without a time",
            folder,
            kind,
            len(timed),
            len(kind_rows) - len(timed),
        )
        for row in timed:
            draw_event(out, data, row, Path(name).stem)
        draw_summary(out, data, kind_rows, Path(name).stem)
    write_tables(out, Path(name).stem, rows)
    saved = {
        "file": name,
        "half_cycles": HALF,
        "template_phase": phases,
        "sigma": (
            "sap_flux_err"
            if flux_err
            else "constant from the saved solution chi-square"
        ),
        "constant_sigma_mmag": None if flux_err else 1e3 * solution_sigma(solution),
        "n_timed": {
            kind: sum(r["kind"] == kind and r["T_jd"] != "" for r in rows)
            for kind, _, _ in KINDS
        },
        "n_rejected": {
            kind: sum(r["kind"] == kind and r["T_jd"] == "" for r in rows)
            for kind, _, _ in KINDS
        },
    }
    (out / "result.json").write_text(json.dumps(saved, indent=2) + "\n")
    logger.info("wrote %s", out)


def main() -> None:
    """Time the maxima of the four NSV 807 light curves."""
    for folder, name, flux_err in RUNS:
        run_one(folder, name, flux_err)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
