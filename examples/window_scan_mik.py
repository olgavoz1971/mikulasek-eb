"""Scan the timing half-window of the Mikulasek template on one light curve.

The solution is the chosen one of ``results/synthetic_scoring/batch_mix_full``.
No fit is run. For each eclipse kind the half-window runs from D/10 in steps of
D/10, where D is the Mikulasek width of that eclipse (``D1`` for the primary,
``D2`` for the secondary, in cycles), up to the largest window ``time`` allows
without overlapping the neighbouring minimum (0.25 cycles).

O is the accepted time and C is the file ephemeris (epoch plus E periods, the
period and epoch stored in the file). No truth column is read. The outlier and
sigma cuts of ``time`` are switched off (``clip_outliers=False``), so that the
only reason for a missing time is the window itself.

The set of eclipses that enters the rms is fixed once, from the data: an eclipse
belongs to it when the widest window holds samples with no gap longer than
``GAP_CADENCES`` cadences and reaches both ends of the window to within the
same tolerance. The same eclipses are used at every width.

Run from the repository root::

    MPLBACKEND=Agg python examples/window_scan_mik.py --physics G1

Writes ``results/synthetic_scoring/window_scan_mik/{physics}_c{cadence}/``:
``scan.csv``, ``oc_points.csv``, ``rms_vs_half_width.png``.
"""

from __future__ import annotations

import argparse
import csv
import logging
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import mikulasek_eb as eb
from mikulasek_eb.timing import fit_ephemeris
from mikulasek_eb.solution import ModelSolution

logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFEST = GENERATED / "results" / "mix" / "manifest.csv"
FULL = ROOT / "results" / "synthetic_scoring" / "batch_mix_full"
OUT = ROOT / "results" / "synthetic_scoring" / "window_scan_mik"
DAY_S = 86400.0
MAX_HALF = 0.25
GAP_CADENCES = 2.0
STEP_FRACTION = 0.1
WIDTH_KEY = {"primary": "D1", "secondary": "D2"}


def find_curve(physics: str, cadence: int) -> dict[str, str]:
    """Manifest row of the light curve to scan (second G1 file excluded).

    Args:
        physics (str): Physics id, for example ``G1``.
        cadence (int): Cadence in seconds.

    Returns:
        dict[str, str]: The manifest row.

    Raises:
        RuntimeError: If there is not exactly one matching written file.
    """
    with MANIFEST.open(newline="") as handle:
        rows = [
            r
            for r in csv.DictReader(handle)
            if r["physics_id"] == physics
            and int(r["cadence_s"]) == cadence
            and int(r["fail"]) == 0
            and r["recipe_id"] != "B22"
        ]
    if len(rows) != 1:
        raise RuntimeError(
            f"expected one {physics} file at {cadence} s, found {len(rows)}"
        )
    return rows[0]


def chosen_solution(curve: Path) -> tuple[ModelSolution, str]:
    """Load the solution chosen by the full comparison for a light curve.

    Args:
        curve (Path): Light-curve file.

    Returns:
        tuple[ModelSolution, str]: The stored solution and its candidate name.

    Raises:
        RuntimeError: If the stored solution is not status ``ok``.
    """
    with (FULL / "chosen.csv").open(newline="") as handle:
        name = next(
            r["chosen"] for r in csv.DictReader(handle) if r["file"] == curve.name
        )
    solution = ModelSolution.load(FULL / "solutions" / f"{curve.stem}__{name}.json")
    if solution.stats["status"] != "ok":
        raise RuntimeError(
            f"{curve.name}: stored solution is {solution.stats['status']}"
        )
    return solution, name


def ephemeris_jd(epoch_jd: float, cycle: float, period: float) -> float:
    """File ephemeris at one eclipse.

    ``cycle`` is the eclipse number from ``time``: an integer for a primary and
    n + 0.5 for a secondary. The centre is the epoch plus that number of periods.
    A secondary is not shifted by a further half period.

    Args:
        epoch_jd (float): File epoch, JD.
        cycle (float): Eclipse number.
        period (float): File period, days.

    Returns:
        float: Predicted time, JD.
    """
    return float(epoch_jd) + float(cycle) * float(period)


def format_cycle(cycle: float) -> str:
    """Write an eclipse number without a negative zero.

    Args:
        cycle (float): Eclipse number.

    Returns:
        str: ``0`` for zero, otherwise the shortest decimal form.
    """
    value = 0.0 if float(cycle) == 0.0 else float(cycle)
    return f"{value:g}"


def integer_cycle(kind: str, cycle: float) -> int:
    """Integer cycle used by the light-curve file.

    A primary keeps its integer. A secondary's integer is the half-integer minus
    one half, so cycle 0.5 is file E 0.

    Args:
        kind (str): ``primary`` or ``secondary``.
        cycle (float): Eclipse number from ``time``.

    Returns:
        int: Integer cycle.

    Raises:
        ValueError: If ``kind`` is unknown.
    """
    if kind == "primary":
        return int(round(float(cycle)))
    if kind == "secondary":
        return int(round(float(cycle) - 0.5))
    raise ValueError(f"Unknown eclipse kind {kind!r}.")


def window_is_covered(
    t_jd: np.ndarray, centre: float, half: float, cadence_d: float
) -> bool:
    """Whether the data cover a window without a long gap.

    Args:
        t_jd (np.ndarray): Sample times, sorted, days.
        centre (float): Predicted time of the eclipse, days.
        half (float): Half-window in days.
        cadence_d (float): Cadence in days.

    Returns:
        bool: True when samples reach both ends and no gap exceeds the tolerance.
    """
    lo, hi = centre - half, centre + half
    inside = t_jd[(t_jd >= lo) & (t_jd <= hi)]
    if len(inside) == 0:
        return False
    tolerance = GAP_CADENCES * cadence_d
    edges = np.concatenate(([lo], inside, [hi]))
    return bool(np.max(np.diff(edges)) <= tolerance)


def scan_kind(
    kind,
    t_mjd,
    mag,
    sig,
    solution,
    epoch_jd,
    period,
    t_jd,
    cadence_d,
    recipe_kwargs=None,
    reference="mean",
    centres_from_fit=False,
):
    """Scan one eclipse kind.

    Args:
        kind (str): ``primary`` or ``secondary``.
        t_mjd (np.ndarray): Times, MJD.
        mag (np.ndarray): Magnitudes.
        sig (np.ndarray): Magnitude errors.
        solution (ModelSolution): Stored template.
        epoch_jd (float): File epoch, JD.
        period (float): File period, days.
        t_jd (np.ndarray): The same times in JD, for the coverage test.
        cadence_d (float): Cadence in days.
        recipe_kwargs (dict | None): Extra ``TimingRecipe`` fields, for example the
            point rule. None keeps the package defaults of the synthetic scans.
        reference (str): ``mean`` takes the rms of O-C about its mean, as in the
            synthetic scans. ``linear`` takes it about a straight line fitted to the
            same times against the eclipse number, so that a stale file period does
            not enter the rms.
        centres_from_fit (bool): False predicts the eclipse centres from the file
            epoch and period, as in the synthetic scans. True takes them from the
            wide-window timing (the fitted time of a timed eclipse, the fitted
            prediction of a rejected one), so that the file ephemeris is not used
            to place an eclipse.

    Returns:
        tuple[list[dict], list[dict], float]: Summary rows, per-eclipse rows and the
        width D in cycles.
    """
    extra = dict(recipe_kwargs or {})
    width = float(solution.params[WIDTH_KEY[kind]])
    step = STEP_FRACTION * width
    halves = np.arange(1, int(MAX_HALF / step) + 1) * step
    halves = halves[halves <= MAX_HALF]
    wide = eb.time(
        t_mjd,
        mag,
        sig,
        solution,
        eb.TimingRecipe(window=float(halves[-1]), clip_outliers=False, **extra),
        which=(kind,),
    )
    predicted = {}
    for item in wide.timings[kind]["E"]:
        predicted[float(item)] = None
    for item in wide.rejected[kind]:
        predicted[float(item["E"])] = float(item["T_pred_jd"])
    cycles = sorted(predicted)
    if centres_from_fit:
        centre = dict(predicted)
        for e, t_fit in zip(wide.timings[kind]["E"], wide.timings[kind]["T"]):
            centre[float(e)] = float(t_fit)
    else:
        centre = {e: ephemeris_jd(epoch_jd, e, period) for e in cycles}
    in_set = {
        e: window_is_covered(t_jd, centre[e], float(halves[-1]) * period, cadence_d)
        for e in cycles
    }
    logger.info(
        "%s: D=%.5f cycles, %d half-windows, %d of %d eclipses in the fixed set",
        kind,
        width,
        len(halves),
        sum(in_set.values()),
        len(cycles),
    )
    summary, points = [], []
    for half in halves:
        out = eb.time(
            t_mjd,
            mag,
            sig,
            solution,
            eb.TimingRecipe(window=float(half), clip_outliers=False, **extra),
            which=(kind,),
        )
        got = out.timings[kind]
        epochs = np.asarray(got["E"], float)
        times = np.asarray(got["T"], float)
        sigma = np.asarray(got["sigma"], float)
        oc = (times - (epoch_jd + epochs * period)) * DAY_S
        in_mask = np.array([in_set.get(float(e), False) for e in epochs], dtype=bool)
        for e, o, s, flag in zip(epochs, oc, sigma, in_mask):
            points.append(
                {
                    "kind": kind,
                    "half": f"{half:.6f}",
                    "cycle": format_cycle(e),
                    "oc_s": f"{o:.3f}",
                    "sigma_s": f"{s * DAY_S:.3f}",
                    "in_set": int(flag),
                }
            )
        n_set = sum(in_set.values())
        used = oc[in_mask]
        used_sigma = sigma[in_mask] * DAY_S
        row = {
            "kind": kind,
            "half": f"{half:.6f}",
            "half_over_D": f"{half / width:.1f}",
            "n_set": n_set,
            "n_timed_in_set": int(in_mask.sum()),
            "n_timed_all": len(epochs),
            "mean_oc_s": "",
            "rms_oc_s": "",
            "rms_over_sigma": "",
            "median_sigma_s": "",
            "n_rejected": str(len(out.rejected[kind])),
            "reasons": "; ".join(
                f"{name} {count}"
                for name, count in Counter(
                    item["reason"] for item in out.rejected[kind]
                ).most_common()
            ),
        }
        if len(used) >= 2:
            mean = float(np.mean(used))
            row["mean_oc_s"] = f"{mean:.3f}"
            if reference == "linear" and len(used) >= 3:
                _, res_d, _ = fit_ephemeris(
                    times[in_mask], sigma[in_mask], epochs[in_mask], 1
                )
                about = res_d * DAY_S
            else:
                about = used - mean
            row["rms_oc_s"] = f"{float(np.sqrt(np.mean(about**2))):.3f}"
            row["rms_over_sigma"] = (
                f"{float(np.sqrt(np.mean((about / used_sigma) ** 2))):.3f}"
            )
            row["median_sigma_s"] = f"{float(np.median(used_sigma)):.3f}"
        summary.append(row)
    return summary, points, width


# Adopted MAVKA half-widths at 120 s, one pair per star, from the batch tables.
MAVKA_HALF_120 = {
    "G1": {"primary": 0.0155, "secondary": 0.0156},
    "G2": {"primary": 0.0187, "secondary": 0.0187},
    "G4": {"primary": 0.0421, "secondary": 0.0454},
    "G5": {"primary": 0.0890, "secondary": 0.0934},
    "G6": {"primary": 0.0790, "secondary": 0.0814},
    "G7": {"primary": 0.0400, "secondary": 0.0400},
    "G8": {"primary": 0.0113, "secondary": 0.0168},
}


def draw(
    summary: list[dict],
    widths: dict[str, float],
    path: Path,
    title: str,
    mavka: dict[str, float] | None = None,
    adopted: dict[str, float] | None = None,
) -> None:
    """Plot rms of O-C and rms of O-C over sigma against half-window over D.

    Args:
        summary (list[dict]): Rows of ``scan.csv``.
        widths (dict[str, float]): D by eclipse kind, cycles.
        path (Path): Output PNG.
        title (str): Figure title.
        mavka (dict[str, float] | None): Adopted MAVKA half-width in cycles, by eclipse
            kind. Drawn as a vertical line at that half-width divided by D.
        adopted (dict[str, float] | None): Adopted half-width in cycles, by eclipse
            kind. Drawn as a dotted vertical line labelled with its multiple of D.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for kind, colour in (("primary", "tab:blue"), ("secondary", "tab:red")):
        rows = [r for r in summary if r["kind"] == kind and r["rms_oc_s"] != ""]
        x = [float(r["half"]) / widths[kind] for r in rows]
        axes[0].plot(
            x, [float(r["rms_oc_s"]) for r in rows], ".-", color=colour, label=kind
        )
        axes[1].plot(
            x,
            [float(r["rms_over_sigma"]) for r in rows],
            ".-",
            color=colour,
            label=kind,
        )
        if mavka is not None and widths[kind] > 0:
            ratio = mavka[kind] / widths[kind]
            for ax in axes:
                ax.axvline(
                    ratio,
                    color=colour,
                    ls="--",
                    lw=1,
                    label=f"MAVKA {kind} {mavka[kind]:.4f} ({ratio:.2f} D)",
                )
    if adopted is not None:
        for kind, colour in (("primary", "tab:blue"), ("secondary", "tab:red")):
            ratio = adopted[kind] / widths[kind]
            for ax in axes:
                ax.axvline(
                    ratio,
                    color=colour,
                    ls=":",
                    lw=1.5,
                    label=f"adopted {kind} {ratio:.1f} D",
                )
    axes[0].set_ylabel("rms of O-C on the fixed set (s)")
    axes[1].set_ylabel("rms of (O-C) / sigma")
    axes[0].set_yscale("log")
    axes[1].axhline(1.0, color="grey", lw=0.8)
    for ax in axes:
        ax.set_xlabel("half-window / D")
        ax.grid(alpha=0.3)
        ax.legend()
    if not any(r["rms_oc_s"] for r in summary):
        axes[0].text(
            0.5, 0.5, "no eclipse timed", transform=axes[0].transAxes, ha="center"
        )
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def write_csv(
    path: Path, rows: list[dict], fieldnames: list[str] | None = None
) -> None:
    """Write a list of dicts as CSV.

    Args:
        path (Path): Output file.
        rows (list[dict]): Rows with identical keys. May be empty.
        fieldnames (list[str] | None): Columns used when ``rows`` is empty.
    """
    fields = fieldnames if fieldnames is not None else list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _fixed_set(
    kind, t_mjd, mag, sig, solution, epoch_jd, period, t_jd, cadence_d, half
):
    """List every eclipse at one half-window and the coverage flag of each.

    Args:
        kind (str): ``primary`` or ``secondary``.
        t_mjd (np.ndarray): Times, MJD.
        mag (np.ndarray): Magnitudes.
        sig (np.ndarray): Magnitude errors.
        solution (ModelSolution): Stored template.
        epoch_jd (float): File epoch, JD.
        period (float): File period, days.
        t_jd (np.ndarray): The same times in JD.
        cadence_d (float): Cadence in days.
        half (float): Half-window in cycles. The coverage test uses this width.

    Returns:
        dict[float, bool]: Coverage flag by eclipse number.
    """
    wide = eb.time(
        t_mjd,
        mag,
        sig,
        solution,
        eb.TimingRecipe(window=float(half), clip_outliers=False),
        which=(kind,),
    )
    cycles = {float(e) for e in wide.timings[kind]["E"]} | {
        float(item["E"]) for item in wide.rejected[kind]
    }
    return {
        cycle: window_is_covered(
            t_jd, ephemeris_jd(epoch_jd, cycle, period), half * period, cadence_d
        )
        for cycle in cycles
    }


def _rows_at(points, kind, half):
    """Return the stored O-C rows of one eclipse kind at one half-window.

    Args:
        points (list[dict]): Rows of ``oc_points.csv``.
        kind (str): ``primary`` or ``secondary``.
        half (str): Half-window text as stored in the table.

    Returns:
        list[dict]: Matching rows.
    """
    return [row for row in points if row["kind"] == kind and row["half"] == half]


def correct_saved_sets() -> None:
    """Recompute the fixed-set flags of every saved scan and redraw the figures.

    The eclipse number from ``time`` is already n + 0.5 for a secondary. An earlier
    version of this script added another half period when it tested coverage, so the
    secondary flags were taken on the following primary. The stored O-C values used
    epoch plus E periods and are left as they are. One ``time`` call per eclipse
    kind recovers the cycle list. The grid is not refitted.
    """
    for folder in sorted(path for path in OUT.iterdir() if path.is_dir()):
        physics, cadence_text = folder.name.split("_c")
        cadence = int(cadence_text)
        with (folder / "scan.csv").open(newline="") as handle:
            summary = list(csv.DictReader(handle))
        with (folder / "oc_points.csv").open(newline="") as handle:
            points = list(csv.DictReader(handle))
        curve = GENERATED / find_curve(physics, cadence)["path"]
        solution, name = chosen_solution(curve)
        payload = eb.load_tess_votable(curve)
        t_jd = np.asarray(payload["t"], float)
        t_mjd = eb.absolute_jd_to_mjd(t_jd)
        period = float(payload["period"])
        epoch_jd = float(payload["epoch_jd"])
        flags = {}
        for kind in ("primary", "secondary"):
            halves = [row["half"] for row in summary if row["kind"] == kind]
            if not halves:
                continue
            flags[kind] = _fixed_set(
                kind,
                t_mjd,
                payload["mag"],
                payload["sig"],
                solution,
                epoch_jd,
                period,
                t_jd,
                cadence / DAY_S,
                float(halves[-1]),
            )
        old_members = {}
        for row in points:
            key = (row["kind"], row["half"])
            old_members.setdefault(key, set())
            if int(row["in_set"]) == 1:
                old_members[key].add(float(row["cycle"]))
        for row in points:
            cycle = float(row["cycle"])
            new_flag = flags[row["kind"]].get(cycle)
            if new_flag is None:
                raise RuntimeError(
                    f"{folder.name} {row['kind']} cycle {row['cycle']} is not in the wide call."
                )
            if row["kind"] == "primary" and int(new_flag) != int(row["in_set"]):
                raise RuntimeError(
                    f"{folder.name} primary cycle {row['cycle']} changed coverage flag."
                )
            row["in_set"] = str(int(new_flag))
            row["cycle"] = format_cycle(cycle)
        for row in summary:
            kind = row["kind"]
            if kind not in flags:
                continue
            n_set = str(sum(flags[kind].values()))
            if kind == "primary" and n_set != row["n_set"]:
                raise RuntimeError(
                    f"{folder.name} primary set size changed from {row['n_set']} to {n_set}."
                )
            used = [
                item
                for item in _rows_at(points, kind, row["half"])
                if int(item["in_set"]) == 1
            ]
            members = {float(item["cycle"]) for item in used}
            if members == old_members.get((kind, row["half"]), set()):
                row["n_set"] = n_set
                continue
            previous = (
                row["n_set"],
                row["n_timed_in_set"],
                row["mean_oc_s"],
                row["rms_oc_s"],
            )
            row["n_set"] = n_set
            row["n_timed_in_set"] = str(len(used))
            row["mean_oc_s"] = ""
            row["rms_oc_s"] = ""
            row["rms_over_sigma"] = ""
            row["median_sigma_s"] = ""
            if len(used) >= 2:
                oc = np.array([float(item["oc_s"]) for item in used])
                sigma = np.array([float(item["sigma_s"]) for item in used])
                mean = float(np.mean(oc))
                row["mean_oc_s"] = f"{mean:.3f}"
                row["rms_oc_s"] = f"{float(np.sqrt(np.mean((oc - mean) ** 2))):.3f}"
                row["rms_over_sigma"] = (
                    f"{float(np.sqrt(np.mean(((oc - mean) / sigma) ** 2))):.3f}"
                )
                row["median_sigma_s"] = f"{float(np.median(sigma)):.3f}"
            logger.info(
                "%s %s half %s: n_set %s -> %s, timed %s -> %s, mean %s -> %s, rms %s -> %s",
                folder.name,
                kind,
                row["half"],
                previous[0],
                row["n_set"],
                previous[1],
                row["n_timed_in_set"],
                previous[2],
                row["mean_oc_s"],
                previous[3],
                row["rms_oc_s"],
            )
        write_csv(folder / "scan.csv", summary)
        write_csv(
            folder / "oc_points.csv",
            points,
            ["kind", "half", "cycle", "oc_s", "sigma_s", "in_set"],
        )
        widths = {
            "primary": float(solution.params["D1"]),
            "secondary": float(solution.params["D2"]),
        }
        draw(
            summary,
            widths,
            folder / "rms_vs_half_width.png",
            f"{curve.stem}, {name}: Mikulasek half-window scan",
        )
        if cadence == 120:
            mark_saved(physics)
        logger.info("corrected %s", folder.name)


def mark_saved(physics: str) -> None:
    """Redraw one 120 s scan from its table, with the MAVKA half-width marked.

    No light curve is timed. D is read from the chosen solution of that file, and the
    MAVKA half-width is the adopted 120 s value, so the mark is ``half / D``.

    Args:
        physics (str): Physics id, one of the first-batch stars.

    Raises:
        KeyError: If the star has no adopted MAVKA half-width.
    """
    folder = OUT / f"{physics}_c120"
    with (folder / "scan.csv").open(newline="") as handle:
        summary = list(csv.DictReader(handle))
    curve = GENERATED / find_curve(physics, 120)["path"]
    solution, name = chosen_solution(curve)
    widths = {
        "primary": float(solution.params["D1"]),
        "secondary": float(solution.params["D2"]),
    }
    mavka = MAVKA_HALF_120[physics]
    for kind in ("primary", "secondary"):
        logger.info(
            "%s %s: MAVKA half %.4f / D %.5f = %.2f D",
            physics,
            kind,
            mavka[kind],
            widths[kind],
            mavka[kind] / widths[kind],
        )
    path = folder / "rms_vs_half_width_mavka.png"
    draw(
        summary,
        widths,
        path,
        f"{curve.stem}, {name}: half-window scan, MAVKA window marked",
        mavka,
    )
    logger.info("wrote %s", path)


def main() -> None:
    """Scan one star at one cadence and write the tables and the figure."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--physics", default="G1")
    parser.add_argument("--cadence", type=int, default=120)
    parser.add_argument(
        "--mark-only",
        action="store_true",
        help="Redraw the 120 s figure from scan.csv, with the MAVKA half-width marked. No timing.",
    )
    parser.add_argument(
        "--correct-sets",
        action="store_true",
        help="Recompute the fixed-set flags from the file ephemeris and redraw. No grid is refitted.",
    )
    args = parser.parse_args()
    if args.correct_sets:
        correct_saved_sets()
        return
    if args.mark_only:
        mark_saved(args.physics)
        return
    manifest_row = find_curve(args.physics, args.cadence)
    curve = GENERATED / manifest_row["path"]
    solution, name = chosen_solution(curve)
    payload = eb.load_tess_votable(curve)
    t_jd = np.asarray(payload["t"], float)
    t_mjd = eb.absolute_jd_to_mjd(t_jd)
    period = float(payload["period"])
    epoch_jd = float(payload["epoch_jd"])
    logger.info(
        "%s, candidate %s, P=%.6f d, epoch JD %.5f", curve.name, name, period, epoch_jd
    )
    summary: list[dict] = []
    points: list[dict] = []
    widths: dict[str, float] = {}
    for kind in ("primary", "secondary"):
        rows, pts, widths[kind] = scan_kind(
            kind,
            t_mjd,
            payload["mag"],
            payload["sig"],
            solution,
            epoch_jd,
            period,
            t_jd,
            args.cadence / DAY_S,
        )
        summary += rows
        points += pts
    folder = OUT / f"{args.physics}_c{args.cadence}"
    folder.mkdir(parents=True, exist_ok=True)
    write_csv(folder / "scan.csv", summary)
    write_csv(
        folder / "oc_points.csv",
        points,
        ["kind", "half", "cycle", "oc_s", "sigma_s", "in_set"],
    )
    draw(
        summary,
        widths,
        folder / "rms_vs_half_width.png",
        f"{curve.stem}, {name}: Mikulasek half-window scan",
    )
    logger.info("wrote %s", folder)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
