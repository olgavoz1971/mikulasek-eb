"""Magnitude difference of the two maxima of NSV 807 in one cycle.

A cycle is kept only when the magnitude step timed both the primary maximum
and the secondary maximum of that cycle. The windows are the NSV 807 maxima
windows, 0.16 cycles for both. The difference is the fitted-template magnitude
of the primary minus that of the secondary. Its error is the root-sum-square
of the two magnitude errors, because the two maxima are separate fits. The
time is the time of the primary maximum. Cycles with only one maximum are
omitted.

Run from ``mikulasek_eb``::

    python examples/nsv807_max_diff.py

Writes one ``.dat`` file per light curve in ``results/nsv807/timing_max_diff_dat/``.
"""

from __future__ import annotations

import logging
from decimal import Decimal
from pathlib import Path

import nsv807_collect_dat as minima
import nsv807_maxima as maxima
import v1383_max_diff as diff

logger = logging.getLogger(__name__)

OUT = maxima.RESULTS / "timing_max_diff_dat"


def write_one(
    path: Path, label: str, rows: list[tuple[Decimal, Decimal, Decimal]]
) -> None:
    """Write one light-curve file.

    Args:
        path (Path): Destination.
        label (str): Light-curve description.
        rows (list[tuple[Decimal, Decimal, Decimal]]): Primary time, difference
            and its error.
    """
    lines = [
        "Star: NSV 807 (TIC 35119266), W UMa type",
        f"Light curve: {label}",
        "Window: half-window 0.16 cycles for both maxima, centred on the brightest "
        "phase of that sector's template",
        "Difference of the fitted-template magnitudes of the primary maximum and the "
        "secondary maximum of the same cycle. A cycle with only one maximum is omitted.",
        "Columns: jd  mag_prim_minus_mag_sec  err_dmag",
        "jd is the time of the primary maximum",
        "mag_prim_minus_mag_sec is the primary magnitude minus the secondary magnitude",
        "err_dmag is the formal 1 sigma error of that difference, the root-sum-square "
        "of the two magnitude errors",
    ]
    text = [f"# {line}" for line in lines]
    for jd, difference, err in rows:
        text.append(f"{jd:.8f}  {difference:.6f}  {err:.6f}")
    path.write_text("\n".join(text) + "\n")


def main() -> None:
    """Write the difference file of every light curve."""
    OUT.mkdir(parents=True, exist_ok=True)
    for lc, (run, stem, label) in minima.LIGHT_CURVES.items():
        path = maxima.RESULTS / run / "maxima_mag" / "eclipses.csv"
        found = diff.timed(path)
        rows = diff.pairs(found)
        name = f"NSV807_{lc}_mag_prim_minus_mag_sec.dat"
        write_one(OUT / name, f"{stem}.vot ({label})", rows)
        logger.info(
            "%s: %d cycles with both maxima, %d primary and %d secondary timed",
            name,
            len(rows),
            len(found["primary_maximum"]),
            len(found["secondary_maximum"]),
        )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
