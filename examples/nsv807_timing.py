"""Window scan and timing of one real light curve, under the locked timing protocol.

The protocol is the section "Protocol of the timing for a real light curve" of
``docs/talk_outline.md``. The template is the saved ``chosen/`` solution of the
tiered fit (``nsv807_s30_tiered.py``) and the points are the same cleaned points,
made by the same functions.

Two modes, run from ``examples``::

    python nsv807_timing.py scan --run ../results/nsv807/s30_protocol
    python nsv807_timing.py figure --run ../results/nsv807/s30_protocol
    python nsv807_timing.py time --run ../results/nsv807/s30_protocol --half-over-d 4.0 4.0

``scan`` reuses ``scan_kind`` of ``window_scan_mik.py`` with the point rule below,
the eclipse centres taken from the fitted timings (not from the file epoch and
period) and the rms of O-C about a straight line. It writes the scan tables and the
adopted window of each eclipse kind. ``time`` times every eclipse at the window
``half-over-d`` times D (primary, secondary), capped at 0.25 cycles, and keeps flags.
For the 600 s file pass the multiples adopted on the 120 s file.

Point rule of both modes: at least 4 points with at least 2 on each side of the
predicted centre. No other gate: the contact-wall check, the sigma cut and the
clip are off, and the wall check and the gaps become flags.

Writes ``timing_scan/`` or ``timing/`` inside the run folder. ``timing/`` holds
``eclipses.csv``, the two ``.dat`` files (full JD, error in days, flag, columns
separated by two spaces), a summary figure per kind and ``fits/<kind>/E_<n>.png``
for every timed eclipse.
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
from pathlib import Path

import numpy as np

import mikulasek_eb as eb
import nsv807_s30_tiered as tiered
from window_scan_mik import (
    DAY_S,
    MAX_HALF,
    WIDTH_KEY,
    scan_kind,
    window_is_covered,
    draw,
    write_csv,
)

logger = logging.getLogger(__name__)

MIN_POINTS = 4
MIN_EACH_SIDE = 2
KINDS = ("primary", "secondary")
POINT_RULE = {"min_points": MIN_POINTS, "min_each_side": MIN_EACH_SIDE}


def prepare(
    run: Path,
    name: str,
    flux: str,
    flux_err: str | None,
    skip_outliers: bool = False,
    solution_path: Path | None = None,
) -> tuple:
    """Load the cleaned light curve and the chosen solution of a tiered-fit run.

    Args:
        run (Path): Output folder of ``nsv807_s30_tiered.py``.
        name (str): Light-curve file name, or an absolute path.
        flux (str): Flux column name.
        flux_err (str | None): Error column name. None when the file has none: the
            constant sigma is then found by the pilot fit of the tiered fit.
        skip_outliers (bool): True when the tiered fit skipped the outlier step
            (``--skip-outliers``), so that all points are kept here too.
        solution_path (Path | None): Solution to use instead of the one in
            ``chosen``, for example a borrowed solution.

    Returns:
        tuple: ``payload`` (cleaned light curve), ``t_mjd``, the chosen
        ``ModelSolution`` and the median cadence in days.

    Raises:
        FileNotFoundError: If ``chosen`` does not hold exactly one solution.
    """
    tiered.FLUX, tiered.FLUX_ERR = flux, flux_err
    tiered.OUT = run
    path = Path(name) if Path(name).is_absolute() else tiered.DATA / name
    payload = tiered.load_sap(path)
    if not skip_outliers:
        payload = tiered.remove_outliers(payload, path.stem)
    found = (
        [solution_path]
        if solution_path is not None
        else sorted((run / "chosen").glob("*.json"))
    )
    if len(found) != 1:
        raise FileNotFoundError(f"expected one solution in {run / 'chosen'}: {found}")
    solution = eb.ModelSolution.load(found[0])
    if solution.stats["status"] != "ok":
        raise RuntimeError(f"{found[0].name}: status {solution.stats['status']}")
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    if flux_err is None:
        payload = tiered.constant_sigma_from_fit(payload, t_mjd, path.stem)
    return payload, t_mjd, solution, float(np.median(np.diff(t_mjd)))


def adopted_half(rows: list[dict], start_over_d: float = 0.0) -> dict:
    """Choose the window from the scan rows of one eclipse kind.

    Only windows that time the whole fixed set are compared, because an rms taken
    over a few timed eclipses is not comparable with one taken over the whole set.
    The narrowest of them whose rms is within ``rms / sqrt(2 (N - 1))`` of the best
    rms is adopted, N being the number of eclipses of the set.

    Args:
        rows (list[dict]): Scan rows of one kind, in order of window.
        start_over_d (float): Only windows of at least this many D are considered.
            It is the multiple adopted on a light curve of higher cadence, so that a
            scan at lower cadence can only widen the window.

    Returns:
        dict: ``half``, ``half_over_D``, ``rms_oc_s``, ``best_rms_oc_s`` and
        ``n_timed_in_set`` of the adopted row.

    Raises:
        ValueError: If no window at or above ``start_over_d`` times the whole set.
    """
    usable = [
        r
        for r in rows
        if r["rms_oc_s"] != ""
        and int(r["n_timed_in_set"]) == int(r["n_set"])
        and float(r["half_over_D"]) >= start_over_d - 1e-9
    ]
    if not usable:
        raise ValueError("no window of the scan from the start times the whole set")
    best = min(usable, key=lambda r: float(r["rms_oc_s"]))
    n = int(best["n_timed_in_set"])
    tolerance = float(best["rms_oc_s"]) / math.sqrt(2.0 * max(n - 1, 1))
    limit = float(best["rms_oc_s"]) + tolerance
    chosen = next(r for r in usable if float(r["rms_oc_s"]) <= limit)
    return {
        "half": chosen["half"],
        "half_over_D": chosen["half_over_D"],
        "rms_oc_s": chosen["rms_oc_s"],
        "best_rms_oc_s": best["rms_oc_s"],
        "n_timed_in_set": chosen["n_timed_in_set"],
    }


def run_adopt(run: Path, name: str, start_over_d: tuple[float, float]) -> None:
    """Apply the adoption rule again to the saved scan tables of a run.

    The scan itself is not repeated. ``adopted.csv`` and the scan figure are
    rewritten, keeping the D of each kind from the previous ``adopted.csv``.

    Args:
        run (Path): Tiered-fit run folder with ``timing_scan``.
        name (str): Light-curve file name, for the figure title.
        start_over_d (tuple[float, float]): First window considered, in D, for the
            primary and the secondary (the multiples adopted at higher cadence).
    """
    out = run / "timing_scan"
    with (out / "adopted.csv").open(newline="") as handle:
        width = {r["kind"]: r["D_cycles"] for r in csv.DictReader(handle)}
    adopted = []
    for kind in KINDS:
        with (out / f"scan_{kind}.csv").open(newline="") as handle:
            pick = adopted_half(
                list(csv.DictReader(handle)), start_over_d[KINDS.index(kind)]
            )
        pick.update({"kind": kind, "D_cycles": width[kind]})
        adopted.append(pick)
        logger.info(
            "%s adopted half-window %s cycles = %s D (rms %s s, best %s s, N %s)",
            kind,
            pick["half"],
            pick["half_over_D"],
            pick["rms_oc_s"],
            pick["best_rms_oc_s"],
            pick["n_timed_in_set"],
        )
    write_csv(out / "adopted.csv", adopted)
    draw_scan(run, Path(name).stem)


def draw_scan(run: Path, stem: str) -> None:
    """Draw the scan figure from the saved tables of a run.

    Args:
        run (Path): Tiered-fit run folder with ``timing_scan``.
        stem (str): Light-curve name, for the title.
    """
    out = run / "timing_scan"
    summary = []
    for kind in KINDS:
        with (out / f"scan_{kind}.csv").open(newline="") as handle:
            summary.extend(csv.DictReader(handle))
    with (out / "adopted.csv").open(newline="") as handle:
        picked = {r["kind"]: r for r in csv.DictReader(handle)}
    draw(
        summary,
        {k: float(picked[k]["D_cycles"]) for k in KINDS},
        out / "rms_vs_half_width.png",
        f"{stem}\nO-C scan, rms about a straight line, fixed eclipse set",
        adopted={k: float(picked[k]["half"]) for k in KINDS},
    )


def run_scan(
    run: Path, name: str, flux: str, flux_err: str | None, skip_outliers: bool = False
) -> None:
    """Scan the window of both eclipse kinds and write the tables.

    Args:
        run (Path): Tiered-fit run folder.
        name (str): Light-curve file name.
        flux (str): Flux column name.
        flux_err (str | None): Error column name.
    """
    payload, t_mjd, solution, cadence_d = prepare(
        run, name, flux, flux_err, skip_outliers
    )
    t_jd = np.asarray(eb.mjd_to_jd(t_mjd), float)
    out = run / "timing_scan"
    out.mkdir(parents=True, exist_ok=True)
    adopted = []
    for kind in KINDS:
        summary, points, width = scan_kind(
            kind,
            t_mjd,
            payload["mag"],
            payload["sig"],
            solution,
            payload["epoch_jd"],
            payload["period"],
            t_jd,
            cadence_d,
            recipe_kwargs=POINT_RULE,
            reference="linear",
            centres_from_fit=True,
        )
        write_csv(out / f"scan_{kind}.csv", summary)
        write_csv(out / f"oc_points_{kind}.csv", points)
        pick = adopted_half(summary)
        pick.update({"kind": kind, "D_cycles": f"{width:.5f}"})
        adopted.append(pick)
        logger.info(
            "%s adopted half-window %s cycles = %s D (rms %s s, best %s s, N %s)",
            kind,
            pick["half"],
            pick["half_over_D"],
            pick["rms_oc_s"],
            pick["best_rms_oc_s"],
            pick["n_timed_in_set"],
        )
    write_csv(out / "adopted.csv", adopted)
    draw_scan(run, Path(name).stem)
    logger.info("wrote %s", out)


def recipe(half: float, walls: bool) -> eb.TimingRecipe:
    """Timing recipe of the locked protocol.

    Args:
        half (float): Half-window in cycles.
        walls (bool): Whether the contact-wall check is applied. It is used only
            to set the ``walls_ok`` flag.

    Returns:
        eb.TimingRecipe: The point rule, no clip, no sigma cut.
    """
    return eb.TimingRecipe(
        window=half, clip_outliers=False, require_contact_walls=walls, **POINT_RULE
    )


def event_fit(
    t_jd: np.ndarray,
    mag: np.ndarray,
    model,
    params: dict,
    p0: float,
    offset: float,
    centre: float,
    half: float,
) -> tuple:
    """Template and baseline of one timed eclipse, for drawing.

    The baseline is an offset and a slope, found by least squares from the
    residual of the points to the template placed at the fitted time.

    Args:
        t_jd (ndarray): Times of all points, JD.
        mag (ndarray): Magnitudes.
        model: Template model of the solution.
        params (dict): Template parameters.
        p0 (float): Fitted period, days.
        offset (float): 0 for a primary, 0.5 for a secondary.
        centre (float): Fitted time of the eclipse, JD.
        half (float): Half-window, cycles.

    Returns:
        tuple: Boolean selection of the points, their phase from the centre
        (cycles), the template there, and the baseline there (mag).
    """
    sel = np.abs(t_jd - centre) <= half * p0
    phi = (t_jd[sel] - centre) / p0
    template = np.asarray(model(offset + phi, params), float)
    design = np.column_stack([np.ones(int(sel.sum())), t_jd[sel] - centre])
    coef = np.linalg.lstsq(design, mag[sel] - template, rcond=None)[0]
    return sel, phi, template, design @ coef


def flag_text(row: dict) -> str:
    """Flag of one timed eclipse as text.

    Args:
        row (dict): Row of ``eclipses.csv`` with a time.

    Returns:
        str: ``ok``, ``gap``, ``walls`` or ``gap+walls``.
    """
    names = []
    if row["gap_in_window"] is True:
        names.append("gap")
    if row["walls_ok"] is False:
        names.append("walls")
    return "+".join(names) if names else "ok"


def write_dat(out: Path, stem: str, rows: list[dict]) -> None:
    """Write one ``.dat`` file per kind: full JD, error in days, flag.

    The columns are separated by two spaces and there is no header.

    Args:
        out (Path): Folder.
        stem (str): Light-curve name.
        rows (list[dict]): Rows of ``eclipses.csv``.
    """
    for kind in KINDS:
        with (out / f"{stem}_{kind}.dat").open("w") as handle:
            for r in rows:
                if r["kind"] == kind and r["T_jd"] != "":
                    err = float(r["sigma_s"]) / DAY_S
                    handle.write(f"{float(r['T_jd']):.8f}  {err:.8f}  {flag_text(r)}\n")


def draw_each(
    out: Path,
    t_jd: np.ndarray,
    mag: np.ndarray,
    sig: np.ndarray,
    solution,
    rows: list[dict],
    stem: str,
) -> None:
    """Draw every timed eclipse on its own: points, fit and residual.

    Args:
        out (Path): Run ``timing`` folder. Figures go to ``fits/<kind>/``.
        t_jd (ndarray): Times of the points, JD.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        solution (ModelSolution): The template.
        rows (list[dict]): Rows of ``eclipses.csv``.
        stem (str): Light-curve name, for the title.
    """
    import matplotlib.pyplot as plt

    model = solution.build_model()
    params = solution.params
    p0 = float(solution.ephemeris["P0"])
    offset = {"primary": 0.0, "secondary": 0.5}
    for r in rows:
        if r["T_jd"] == "":
            continue
        kind, centre, half = r["kind"], float(r["T_jd"]), float(r["half_cycles"])
        sel, phi, template, base = event_fit(
            t_jd, mag, model, params, p0, offset[kind], centre, half
        )
        folder = out / "fits" / kind
        folder.mkdir(parents=True, exist_ok=True)
        fit = template + base
        fig, axes = plt.subplots(
            2, 1, figsize=(8, 7), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
        )
        axes[0].errorbar(
            t_jd[sel] - centre,
            mag[sel],
            yerr=sig[sel],
            fmt=".",
            ms=4,
            color="0.4",
            label="data",
        )
        grid = np.linspace(-half * p0, half * p0, 400)
        pg = grid / p0
        tmpl = np.asarray(model(offset[kind] + pg, params), float)
        coef = np.linalg.lstsq(
            np.column_stack([np.ones(int(sel.sum())), t_jd[sel] - centre]),
            mag[sel] - template,
            rcond=None,
        )[0]
        axes[0].plot(
            grid,
            tmpl + coef[0] + coef[1] * grid,
            color="red",
            lw=1.5,
            label="template and baseline",
        )
        axes[0].axvline(0.0, color="blue", ls=":", lw=1, label="fitted time")
        axes[0].invert_yaxis()
        axes[0].set_ylabel("Magnitude")
        axes[0].legend()
        axes[0].set_title(
            f"{stem}\n{kind} E={r['cycle']}  T={centre:.6f}  sigma {r['sigma_s']} s  "
            f"n={r['n_points']}  {flag_text(r)}"
        )
        axes[1].errorbar(
            t_jd[sel] - centre,
            1e3 * (mag[sel] - fit),
            yerr=1e3 * sig[sel],
            fmt=".",
            ms=4,
            color="0.4",
        )
        axes[1].axhline(0.0, color="red", lw=1)
        axes[1].set_xlabel("Time from the fitted time (d)")
        axes[1].set_ylabel("Residual (mmag)")
        fig.tight_layout()
        fig.savefig(folder / f"E_{r['cycle']}.png", dpi=100)
        plt.close(fig)


def draw_timing(
    out: Path, t_jd: np.ndarray, mag: np.ndarray, solution, rows: list[dict], stem: str
) -> None:
    """Draw how the template sits in the minima, and the O-C of the timed eclipses.

    One figure per eclipse kind with three panels. Top: the points of every timed
    eclipse, shifted to the fitted time of that eclipse and corrected for its own
    baseline (an offset and a slope, found by least squares from the residual to
    the template), with the template in red. Middle: the binned residual of those
    points, in mmag. Bottom: O-C against time about a straight line fitted to the
    timed eclipses, flagged times in orange.

    Args:
        out (Path): Folder for the PNG files.
        t_jd (ndarray): Times of the points, JD.
        mag (ndarray): Magnitudes.
        solution (ModelSolution): The template.
        rows (list[dict]): Rows of ``eclipses.csv``.
        stem (str): Light-curve name, for the title.
    """
    import matplotlib.pyplot as plt
    from mikulasek_eb.timing import fit_ephemeris

    model = solution.build_model()
    params = solution.params
    p0 = float(solution.ephemeris["P0"])
    offset = {"primary": 0.0, "secondary": 0.5}
    for kind in KINDS:
        timed = [r for r in rows if r["kind"] == kind and r["T_jd"] != ""]
        if not timed:
            continue
        half = float(timed[0]["half_cycles"])
        phis, ys = [], []
        for r in timed:
            centre = float(r["T_jd"])
            sel, phi, _, base = event_fit(
                t_jd, mag, model, params, p0, offset[kind], centre, half
            )
            phis.append(phi)
            ys.append(mag[sel] - base)
        phi_all, y_all = np.concatenate(phis), np.concatenate(ys)
        grid = np.linspace(-half, half, 400)
        line = np.asarray(model(offset[kind] + grid, params), float)
        edges = np.linspace(-half, half, 41)
        index = np.digitize(phi_all, edges) - 1
        resid = y_all - np.asarray(model(offset[kind] + phi_all, params), float)
        centres, medians, scatter = [], [], []
        for k in range(40):
            in_bin = index == k
            if in_bin.sum() >= 3:
                centres.append(0.5 * (edges[k] + edges[k + 1]))
                medians.append(1e3 * float(np.median(resid[in_bin])))
                scatter.append(1e3 * float(np.std(resid[in_bin])))
        fig, axes = plt.subplots(3, 1, figsize=(10, 13))
        axes[0].scatter(
            phi_all, y_all, s=3, color="0.6", label="data, own baseline removed"
        )
        axes[0].plot(
            grid, line, color="red", lw=1.5, label="template at the fitted time"
        )
        axes[0].invert_yaxis()
        axes[0].set_xlabel("Phase from the fitted centre of each eclipse")
        axes[0].set_ylabel("Magnitude (template scale)")
        axes[0].legend()
        axes[0].set_title(
            f"{stem}\n{kind}, {len(timed)} timed eclipses, half-window {half:.4f} cycles"
        )
        axes[1].errorbar(
            centres, medians, yerr=scatter, fmt="o", ms=3, color="tab:blue"
        )
        axes[1].axhline(0.0, color="red", lw=1)
        axes[1].set_xlabel("Phase from the fitted centre of each eclipse")
        axes[1].set_ylabel("Median residual per bin (mmag), bar = scatter")
        e = np.array([float(r["cycle"]) for r in timed])
        t = np.array([float(r["T_jd"]) for r in timed])
        sig = np.array([float(r["sigma_s"]) for r in timed]) / DAY_S
        _, res, _ = fit_ephemeris(t, sig, e, 1)
        flagged = np.array(
            [(r["walls_ok"] is False) or r["gap_in_window"] is True for r in timed]
        )
        axes[2].errorbar(
            t[~flagged],
            DAY_S * res[~flagged],
            yerr=DAY_S * sig[~flagged],
            fmt="o",
            ms=3,
            color="tab:blue",
            label="no flag",
        )
        if flagged.any():
            axes[2].errorbar(
                t[flagged],
                DAY_S * res[flagged],
                yerr=DAY_S * sig[flagged],
                fmt="o",
                ms=3,
                color="orange",
                label="walls or gap flag",
            )
        axes[2].axhline(0.0, color="grey", lw=0.8)
        axes[2].set_xlabel("Julian Date")
        axes[2].set_ylabel("O-C about a straight line (s)")
        axes[2].legend()
        fig.tight_layout()
        fig.savefig(out / f"timing_{kind}.png", dpi=130)
        plt.close(fig)


def run_time(
    run: Path,
    name: str,
    flux: str,
    flux_err: str | None,
    half_over_d: tuple[float, float],
    skip_outliers: bool = False,
    solution_path: Path | None = None,
    out_dir: Path | None = None,
) -> None:
    """Time every eclipse and write the table with its flags.

    Args:
        run (Path): Tiered-fit run folder.
        name (str): Light-curve file name.
        flux (str): Flux column name.
        flux_err (str | None): Error column name.
        half_over_d (tuple[float, float]): Window multiples of D for the primary
            and the secondary.
        skip_outliers (bool): The tiered fit skipped the outlier step.
        solution_path (Path | None): Solution to time with instead of ``chosen``.
        out_dir (Path | None): Output folder instead of ``timing`` in the run folder.
    """
    payload, t_mjd, solution, cadence_d = prepare(
        run, name, flux, flux_err, skip_outliers, solution_path
    )
    t_jd = np.asarray(eb.mjd_to_jd(t_mjd), float)
    period = float(payload["period"])
    out = run / "timing" if out_dir is None else out_dir
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for kind, multiple in zip(KINDS, half_over_d):
        half = min(multiple * float(solution.params[WIDTH_KEY[kind]]), MAX_HALF)
        args = (t_mjd, payload["mag"], payload["sig"], solution)
        timed_all = eb.time(*args, recipe(half, False), which=(kind,))
        walled = eb.time(*args, recipe(half, True), which=(kind,))
        passed = {float(e) for e in walled.timings[kind]["E"]}
        got = timed_all.timings[kind]
        for cycle, t_fit, sigma in zip(got["E"], got["T"], got["sigma"]):
            n_in = int(np.sum(np.abs(t_jd - float(t_fit)) <= half * period))
            rows.append(
                {
                    "kind": kind,
                    "cycle": f"{float(cycle):g}",
                    "T_jd": f"{float(t_fit):.6f}",
                    "sigma_s": f"{float(sigma) * DAY_S:.3f}",
                    "half_cycles": f"{half:.5f}",
                    "n_points": n_in,
                    "dof": n_in - 3,
                    "walls_ok": float(cycle) in passed,
                    "gap_in_window": not window_is_covered(
                        t_jd, float(t_fit), half * period, cadence_d
                    ),
                    "reason": "",
                }
            )
        for item in timed_all.rejected[kind]:
            rows.append(
                {
                    "kind": kind,
                    "cycle": f"{float(item['E']):g}",
                    "T_jd": "",
                    "sigma_s": "",
                    "half_cycles": f"{half:.5f}",
                    "n_points": "",
                    "dof": "",
                    "walls_ok": "",
                    "gap_in_window": "",
                    "reason": item["reason"],
                }
            )
        logger.info(
            "%s half %.5f cycles: %d timed, %d without a time",
            kind,
            half,
            len(got["E"]),
            len(timed_all.rejected[kind]),
        )
    write_csv(out / "eclipses.csv", rows)
    stem = Path(name).stem
    write_dat(out, stem, rows)
    draw_timing(out, t_jd, payload["mag"], solution, rows, stem)
    draw_each(out, t_jd, payload["mag"], payload["sig"], solution, rows, stem)
    logger.info("wrote %s", out)


def main() -> None:
    """Parse the command line and run one mode."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("mode", choices=("scan", "time", "figure", "adopt"))
    parser.add_argument("--run", required=True, help="tiered-fit run folder")
    parser.add_argument("--file", default="TIC35119266_s30_SPOC_120s.vot")
    parser.add_argument("--flux", default="sap_flux")
    parser.add_argument("--flux-err", default="sap_flux_err")
    parser.add_argument(
        "--no-flux-err",
        action="store_true",
        help="the file has no error column: constant sigma from the pilot fit",
    )
    parser.add_argument(
        "--skip-outliers",
        action="store_true",
        help="the tiered fit skipped the outlier step",
    )
    parser.add_argument(
        "--half-over-d",
        type=float,
        nargs=2,
        metavar=("PRIMARY", "SECONDARY"),
        help="window multiples of D for the time mode",
    )
    parser.add_argument(
        "--start-over-d",
        type=float,
        nargs=2,
        default=(0.0, 0.0),
        metavar=("PRIMARY", "SECONDARY"),
        help="adopt mode: first window considered, in D (multiples adopted at higher cadence)",
    )
    parser.add_argument(
        "--out-name",
        help="time mode: write to this folder name in the run folder instead of timing",
    )
    args = parser.parse_args()
    run = Path(args.run).resolve()
    if args.no_flux_err:
        args.flux_err = None
    if args.mode == "figure":
        draw_scan(run, Path(args.file).stem)
    elif args.mode == "adopt":
        run_adopt(run, args.file, tuple(args.start_over_d))
    elif args.mode == "scan":
        run_scan(run, args.file, args.flux, args.flux_err, args.skip_outliers)
    else:
        if args.half_over_d is None:
            parser.error("time mode needs --half-over-d")
        run_time(
            run,
            args.file,
            args.flux,
            args.flux_err,
            tuple(args.half_over_d),
            args.skip_outliers,
            out_dir=None if args.out_name is None else run / args.out_name,
        )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
