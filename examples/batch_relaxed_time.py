"""Time the mixed-trend batch with Mikulasek, relaxed gates, maximal window.

Every eclipse of every file is timed with the stored chosen solution and a
half-window of 0.25 cycles, the widest ``time`` accepts without overlap. Two
gates of the strict run are relaxed: the minimum is 6 points per window (the
floor given to MAVKA, with 3 on each side as in the package) and the
contact-wall coverage check is off.

The wall check is not dropped silently. A second pass with the same window and
6 points but the wall check on says which relaxed times would have failed it.
Their ``walls_ok`` column is False. A time that is still missing keeps a row
with its reason. Nothing is filled in.

The plateau run (``batch_plateau_time.py``) is not changed.

Run from the repository root::

    python examples/batch_relaxed_time.py
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

import mikulasek_eb as eb
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
OUT = ROOT / "results" / "synthetic_scoring" / "batch_relaxed_mik"
HALF = 0.25
MIN_POINTS = 6


def recipe(walls: bool) -> eb.TimingRecipe:
    """Timing recipe of this run.

    Args:
        walls (bool): Whether the contact-wall check is applied.

    Returns:
        eb.TimingRecipe: Window 0.25, 6 points, outlier cut off.
    """
    return eb.TimingRecipe(
        window=HALF,
        clip_outliers=False,
        min_points=MIN_POINTS,
        require_contact_walls=walls,
    )


def time_file(row: dict) -> tuple[list[dict], list[dict]]:
    """Time one file with the relaxed recipe and flag the wall check.

    Args:
        row (dict): Manifest row.

    Returns:
        tuple[list[dict], list[dict]]: Summary rows per kind and per-eclipse rows.
    """
    curve = GENERATED / row["path"]
    solution, name = chosen_solution(curve)
    payload = eb.load_tess_votable(curve)
    t_mjd = eb.absolute_jd_to_mjd(np.asarray(payload["t"], float))
    epoch, period = float(payload["epoch_jd"]), float(payload["period"])
    summary, points = [], []
    for kind in ("primary", "secondary"):
        args = (t_mjd, payload["mag"], payload["sig"], solution)
        relaxed = eb.time(*args, recipe(False), which=(kind,))
        walled = eb.time(*args, recipe(True), which=(kind,))
        passed = {float(e) for e in walled.timings[kind]["E"]}
        timed = relaxed.timings[kind]
        oc, flagged = [], 0
        for cycle, t_jd, sigma in zip(timed["E"], timed["T"], timed["sigma"]):
            ok = float(cycle) in passed
            flagged += not ok
            oc_s = (float(t_jd) - (epoch + float(cycle) * period)) * DAY
            oc.append(oc_s)
            points.append(
                {
                    "file": curve.name,
                    "kind": kind,
                    "cycle": f"{float(cycle):.1f}",
                    "T_jd": f"{float(t_jd):.6f}",
                    "sigma_s": f"{float(sigma) * DAY:.3f}",
                    "oc_file_s": f"{oc_s:.3f}",
                    "walls_ok": str(ok),
                    "reason": "",
                }
            )
        for item in relaxed.rejected[kind]:
            points.append(
                {
                    "file": curve.name,
                    "kind": kind,
                    "cycle": f"{float(item['E']):.1f}",
                    "T_jd": "",
                    "sigma_s": "",
                    "oc_file_s": "",
                    "walls_ok": "",
                    "reason": item["reason"],
                }
            )
        measured = np.asarray(timed["T"], float) if len(oc) else np.array([])
        record = {
            "file": curve.name,
            "physics": row["physics_id"],
            "cadence_s": row["cadence_s"],
            "candidate": name,
            "kind": kind,
            "n_timed": str(len(oc)),
            "n_walls_failed": str(flagged),
            "n_rejected": str(len(relaxed.rejected[kind])),
            "mean_oc_s": f"{float(np.mean(oc)):.2f}" if oc else "",
            "rms_oc_s": f"{float(np.std(oc)):.2f}" if len(oc) >= 2 else "",
            "median_sigma_s": (
                f"{float(np.median(np.asarray(timed['sigma'], float) * DAY)):.2f}"
                if oc
                else ""
            ),
        }
        record.update(truth_offsets(curve, kind, measured))
        summary.append(record)
        logger.info(
            "%s %s n=%d walls_failed=%d rejected=%s vs conj %s s",
            curve.name,
            kind,
            len(oc),
            flagged,
            record["n_rejected"],
            record.get("bias_vs_conj_s", ""),
        )
    return summary, points


def main() -> None:
    """Time every accepted file and write the tables."""
    OUT.mkdir(parents=True, exist_ok=True)
    summary, points = [], []
    for row in manifest_rows():
        part_summary, part_points = time_file(row)
        summary.extend(part_summary)
        points.extend(part_points)
    write_csv(OUT / "summary.csv", summary)
    write_csv(OUT / "eclipses.csv", points)
    logger.info("wrote %s", OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
