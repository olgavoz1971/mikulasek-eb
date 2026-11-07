"""Time the mixed-trend batch with one shift of the whole template per cycle.

The stored solution is not refitted. Profile, depths, widths, Gamma and harmonics
stay as they are. Each orbital cycle gets one time shift, and that shift moves
the primary, the secondary and the harmonics together. There is no eclipse
window. The points of cycle n are those with phase in (n - 0.25, n + 0.75), so
the primary at phase n and the secondary at phase n + 0.5 both lie inside the
cycle, a quarter of a cycle from either edge.

The per-cycle baseline is the timing baseline, an offset and a slope. A cycle
with fewer than 6 points, or fewer than 3 on either side of its centre, is not
timed and keeps its reason. Nothing is filled in.

The second G1 600 s file is timed and left out of the figure, as in the method
comparison. G6 is included.

Run from ``examples``::

    python batch_whole_cycle_time.py
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import mikulasek_eb as eb
from batch_compare_methods import LABELS
from batch_plateau_time import (
    DAY,
    GENERATED,
    ROOT,
    chosen_solution,
    manifest_rows,
    truth_offsets,
    write_csv,
)

logger = logging.getLogger(__name__)
OUT = ROOT / "results" / "synthetic_scoring" / "batch_whole_cycle"
RELAXED = ROOT / "results" / "synthetic_scoring" / "batch_relaxed_mik" / "summary.csv"
MIN_POINTS = 6
MIN_SIDE = 3
STARS = ["G1", "G2", "G4", "G5", "G6", "G7", "G8"]
CADENCES = [120, 600, 1800]


def time_file(row: dict, order: int = 1) -> tuple[list[dict], list[dict]]:
    """Time one file, one shift per cycle.

    Args:
        row (dict): Manifest row.
        order (int): Degree of the straight-line-or-higher baseline fitted to each cycle.

    Returns:
        tuple[list[dict], list[dict]]: One summary row per eclipse kind, and one
        row per eclipse. A cycle that fails keeps two eclipse rows with the reason.
    """
    curve = GENERATED / row["path"]
    solution, name = chosen_solution(curve)
    payload = eb.load_tess_votable(curve)
    t_mjd = eb.absolute_jd_to_mjd(np.asarray(payload["t"], float))
    mag = np.asarray(payload["mag"], float)
    sig = np.asarray(payload["sig"], float)
    epoch, period = float(payload["epoch_jd"]), float(payload["period"])
    ephem = solution.ephem
    model = solution.build_model()
    phase = eb.phase(t_mjd, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0))
    n_min = int(np.floor(float(np.min(phase)) + 0.25))
    n_max = int(np.floor(float(np.max(phase)) + 0.25 - 1e-9))
    measured = {"primary": [], "secondary": []}
    points = []
    for n in range(n_min, n_max + 1):
        centre = n + 0.25
        offset = phase - centre
        inside = np.abs(offset) < 0.5
        n_in = int(inside.sum())
        n_before = int(np.sum(offset[inside] < 0))
        n_after = n_in - n_before
        reason = ""
        if n_in < MIN_POINTS:
            reason = "too_few_points"
        elif n_before < MIN_SIDE or n_after < MIN_SIDE:
            reason = "one_side_uncovered"
        delta = sigma = None
        if not reason:
            try:
                fit = eb.fit_minimum(
                    t_mjd,
                    mag,
                    sig,
                    model,
                    solution.p,
                    ephem,
                    centre,
                    window=0.5,
                    order=order,
                    free=(),
                )
                delta, sigma = fit.delta, fit.sigma
            except Exception as exc:
                reason = f"fit_failed: {type(exc).__name__}"
                logger.warning("%s cycle %s: %s", curve.name, n, reason)
        for kind, eclipse in (("primary", float(n)), ("secondary", n + 0.5)):
            record = {
                "file": curve.name,
                "kind": kind,
                "cycle": f"{eclipse:.1f}",
                "T_jd": "",
                "sigma_s": "",
                "oc_file_s": "",
                "reason": reason,
            }
            if delta is not None:
                t_mjd_min = eb.predicted_time(
                    eclipse, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0)
                )
                t_jd = float(eb.mjd_to_jd(t_mjd_min + delta))
                nearest = round((t_jd - epoch) / period * 2) / 2
                oc_s = (t_jd - (epoch + nearest * period)) * DAY
                record.update(
                    T_jd=f"{t_jd:.6f}",
                    sigma_s=f"{sigma * DAY:.3f}",
                    oc_file_s=f"{oc_s:.3f}",
                )
                measured[kind].append(t_jd)
            points.append(record)
    summary = []
    for kind in ("primary", "secondary"):
        times = measured[kind]
        rejected = [p for p in points if p["kind"] == kind and p["reason"]]
        oc = [float(p["oc_file_s"]) for p in points if p["kind"] == kind and p["T_jd"]]
        sigmas = [
            float(p["sigma_s"]) for p in points if p["kind"] == kind and p["sigma_s"]
        ]
        record = {
            "file": curve.name,
            "physics": row["physics_id"],
            "cadence_s": row["cadence_s"],
            "recipe_id": row["recipe_id"],
            "candidate": name,
            "kind": kind,
            "baseline_order": str(order),
            "n_timed": str(len(times)),
            "n_rejected": str(len(rejected)),
            "mean_oc_s": f"{float(np.mean(oc)):.2f}" if oc else "",
            "rms_oc_s": f"{float(np.std(oc)):.2f}" if len(oc) >= 2 else "",
            "median_sigma_s": f"{float(np.median(sigmas)):.2f}" if sigmas else "",
        }
        record.update(truth_offsets(curve, kind, np.asarray(times, float)))
        summary.append(record)
        logger.info(
            "%s %s n=%d rejected=%d vs conj %s s",
            curve.name,
            kind,
            len(times),
            len(rejected),
            record.get("bias_vs_conj_s", ""),
        )
    return summary, points


def draw(summary: list[dict]) -> None:
    """Scatter about the pure minimum, whole cycle against the relaxed window run.

    Args:
        summary (list[dict]): Rows of this run. The second G1 600 s file is omitted.
    """
    relaxed = list(csv.DictReader(RELAXED.open()))
    fig, axes = plt.subplots(
        2, len(STARS), figsize=(3.0 * len(STARS), 6.2), sharex=True
    )
    for i, kind in enumerate(("primary", "secondary")):
        for j, star in enumerate(STARS):
            ax = axes[i, j]
            for rows, colour, label, ls in (
                (summary, "red", "whole cycle", "-"),
                (relaxed, "tab:purple", "window 0.25, relaxed", "--"),
            ):
                y = []
                for cadence in CADENCES:
                    hit = [
                        r
                        for r in rows
                        if r["physics"] == star
                        and int(r["cadence_s"]) == cadence
                        and r["kind"] == kind
                        and "B22" not in r["file"]
                    ]
                    y.append(
                        float(hit[0]["scatter_vs_pure_s"])
                        if hit and hit[0]["scatter_vs_pure_s"]
                        else np.nan
                    )
                ax.plot(range(3), y, color=colour, ls=ls, marker="o", ms=4, label=label)
            ax.set_xticks(range(3), [str(c) for c in CADENCES])
            ax.set_yscale("log")
            ax.grid(alpha=0.3)
            if i == 0:
                ax.set_title(f"{LABELS[star]} ({star})", fontsize=9)
            if j == 0:
                ax.set_ylabel(f"{kind}: scatter about the pure minimum (s)")
            if i == 1:
                ax.set_xlabel("cadence (s)")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle(
        "One shift of the whole template per cycle, against the relaxed 0.25-cycle windows."
    )
    fig.tight_layout()
    fig.savefig(OUT / "scatter_vs_pure.png", dpi=130)
    plt.close(fig)


def run(order: int, dest: Path, physics: set[str] | None) -> None:
    """Time the batch at one baseline degree and write a new directory.

    Args:
        order (int): Baseline polynomial degree of each cycle.
        dest (Path): Output directory. The straight-line run is not a valid destination.
        physics (set[str] | None): Stars to time, or every star when None.

    Raises:
        RuntimeError: If ``dest`` is the directory of the straight-line run.
    """
    if dest.resolve() == OUT.resolve():
        raise RuntimeError(f"{dest} holds the straight-line run and is not overwritten")
    dest.mkdir(parents=True, exist_ok=True)
    summary, points = [], []
    for row in manifest_rows():
        if physics is not None and row["physics_id"] not in physics:
            continue
        part_summary, part_points = time_file(row, order)
        summary.extend(part_summary)
        points.extend(part_points)
    write_csv(dest / "summary.csv", summary)
    write_csv(dest / "eclipses.csv", points)
    logger.info("wrote %s", dest)


def compare() -> None:
    """Draw the straight-line, quadratic and cubic whole-cycle runs on one figure."""
    folders = {
        "baseline degree 1": OUT / "summary.csv",
        "baseline degree 2": OUT.parent / "batch_whole_cycle_b2" / "summary.csv",
        "baseline degree 3": OUT.parent / "batch_whole_cycle_b3" / "summary.csv",
        "window 0.25, relaxed": RELAXED,
    }
    colours = {
        "baseline degree 1": "red",
        "baseline degree 2": "tab:orange",
        "baseline degree 3": "tab:green",
        "window 0.25, relaxed": "tab:purple",
    }
    styles = ["-", "--", ":", "-."]
    loaded = {name: list(csv.DictReader(path.open())) for name, path in folders.items()}
    fig, axes = plt.subplots(
        2, len(STARS), figsize=(3.0 * len(STARS), 6.4), sharex=True
    )
    for i, kind in enumerate(("primary", "secondary")):
        for j, star in enumerate(STARS):
            ax = axes[i, j]
            for (name, rows), ls in zip(loaded.items(), styles):
                y = []
                for cadence in CADENCES:
                    hit = [
                        r
                        for r in rows
                        if r["physics"] == star
                        and int(r["cadence_s"]) == cadence
                        and r["kind"] == kind
                        and "B22" not in r["file"]
                    ]
                    y.append(
                        float(hit[0]["scatter_vs_pure_s"])
                        if hit and hit[0]["scatter_vs_pure_s"]
                        else np.nan
                    )
                ax.plot(
                    range(3),
                    y,
                    color=colours[name],
                    ls=ls,
                    marker="o",
                    ms=4,
                    label=name,
                )
            ax.set_xticks(range(3), [str(c) for c in CADENCES])
            ax.set_yscale("log")
            ax.grid(alpha=0.3)
            if i == 0:
                ax.set_title(f"{LABELS[star]} ({star})", fontsize=9)
            if j == 0:
                ax.set_ylabel(f"{kind}: scatter about the pure minimum (s)")
            if i == 1:
                ax.set_xlabel("cadence (s)")
    axes[0, 0].legend(fontsize=6)
    fig.suptitle("Whole-cycle shift. The per-cycle baseline degree is the only change.")
    fig.tight_layout()
    path = OUT.parent / "batch_whole_cycle_baseline" / "scatter_vs_pure.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    plt.close(fig)
    logger.info("wrote %s", path)


def main() -> None:
    """Time at a higher baseline degree, or draw the comparison. The degree-1 run is kept."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, choices=(2, 3))
    parser.add_argument("--stars", nargs="*", default=None)
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Draw degree 1, 2 and 3 against the relaxed windows. No timing.",
    )
    args = parser.parse_args()
    if args.compare:
        compare()
        return
    if args.order is None:
        raise SystemExit("Pass --order 2 or 3. The degree-1 run is not repeated.")
    dest = OUT.parent / f"batch_whole_cycle_b{args.order}"
    run(args.order, dest, None if args.stars is None else set(args.stars))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
