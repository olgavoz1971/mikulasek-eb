"""Write the fixed eclipse sets and the half-window grids of the 120 s scan.

The files are the hand-over for a second method that repeats the window scan on
the same eclipses and the same grid. Nothing is fitted for the eclipse sets: the
list of cycles comes from the widest ``time`` call, and the coverage rule is the
one of ``window_scan_mik.py``.

Run from the repository root::

    python examples/window_scan_eclipse_sets.py

Writes ``eclipse_sets_120.csv`` and ``grid_120.csv`` in
``results/synthetic_scoring/window_scan_mik/``.
"""

from __future__ import annotations

import csv
import logging

import numpy as np

import mikulasek_eb as eb
from window_scan_mik import (
    DAY_S,
    GENERATED,
    MAX_HALF,
    OUT,
    STEP_FRACTION,
    WIDTH_KEY,
    chosen_solution,
    ephemeris_jd,
    find_curve,
    format_cycle,
    integer_cycle,
    window_is_covered,
    write_csv,
)

logger = logging.getLogger(__name__)
STARS = ("G1", "G2", "G4", "G5", "G6", "G7", "G8")
CADENCE = 120


def main() -> None:
    """Write the eclipse sets and the half-window grids for every 120 s star."""
    sets: list[dict] = []
    grids: list[dict] = []
    for star in STARS:
        row = find_curve(star, CADENCE)
        curve = GENERATED / row["path"]
        solution, name = chosen_solution(curve)
        payload = eb.load_tess_votable(curve)
        t_jd = np.asarray(payload["t"], float)
        t_mjd = eb.absolute_jd_to_mjd(t_jd)
        period = float(payload["period"])
        epoch_jd = float(payload["epoch_jd"])
        for kind in ("primary", "secondary"):
            width = float(solution.params[WIDTH_KEY[kind]])
            step = STEP_FRACTION * width
            k_max = int(MAX_HALF / step)
            halves = np.arange(1, k_max + 1) * step
            halves = halves[halves <= MAX_HALF]
            wide = eb.time(
                t_mjd,
                payload["mag"],
                payload["sig"],
                solution,
                eb.TimingRecipe(window=float(halves[-1]), clip_outliers=False),
                which=(kind,),
            )
            cycles = sorted(
                {float(e) for e in wide.timings[kind]["E"]}
                | {float(r["E"]) for r in wide.rejected[kind]}
            )
            n_in = 0
            for cycle in cycles:
                centre = ephemeris_jd(epoch_jd, cycle, period)
                flag = window_is_covered(
                    t_jd, centre, float(halves[-1]) * period, CADENCE / DAY_S
                )
                n_in += int(flag)
                sets.append(
                    {
                        "physics_id": star,
                        "kind": kind,
                        "cycle": format_cycle(cycle),
                        "integer_E": integer_cycle(kind, cycle),
                        "predicted_jd": f"{centre:.6f}",
                        "in_set": int(flag),
                    }
                )
            grids.append(
                {
                    "physics_id": star,
                    "file": curve.name,
                    "kind": kind,
                    "D_cycles": f"{width:.6f}",
                    "step_cycles": f"{step:.6f}",
                    "step_s": f"{step * period * DAY_S:.1f}",
                    "step_in_cadences": f"{step * period * DAY_S / CADENCE:.2f}",
                    "n_half_windows": len(halves),
                    "first_half": f"{halves[0]:.6f}",
                    "last_half": f"{halves[-1]:.6f}",
                    "n_cycles": len(cycles),
                    "n_in_set": n_in,
                    "candidate": name,
                }
            )
            logger.info(
                "%s %s: D=%.5f, %d half-windows, %d of %d eclipses in the set",
                star,
                kind,
                width,
                len(halves),
                n_in,
                len(cycles),
            )
    write_csv(OUT / "eclipse_sets_120.csv", sets)
    write_csv(OUT / "grid_120.csv", grids)
    logger.info("wrote %s", OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
