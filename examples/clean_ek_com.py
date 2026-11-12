"""Flag folded outliers on EK Com, sector 49, and draw them.

The default of 400 bins raises on this file: phase bin 11 has 4 points,
and a bin needs at least 5 before its scatter can be measured. This example
uses 200 bins and 5 robust sigma, which flags 15 of 2651 points. With 100
bins the same rule flags none. The original file is not modified, and
nothing is written into ``data/``.

Run from the repository root::

    python examples/clean_ek_com.py

Writes ``examples/out/ek_com_clean/folded_outliers.png`` and ``removed.csv``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import mikulasek_eb_tess as tess

from clean_folded import draw_outliers, write_removed

logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "EK_Com_sector_49.vot"
OUT = ROOT / "examples" / "out" / "ek_com_clean"
N_BINS = 200
N_SIGMA = 5.0


def main() -> None:
    """Flag EK Com at 200 bins and save the figure."""
    payload = tess.load_tess_votable(DATA)
    flag, _, phase = tess.flag_folded_outliers(
        payload["t"],
        payload["mag"],
        payload["period"],
        n_bins=N_BINS,
        n_sigma=N_SIGMA,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    n_flag = int(flag.sum())
    title = (
        f"{DATA.name}\n"
        f"{n_flag} of {len(flag)} points beyond {N_SIGMA:.0f} robust sigma "
        f"of their phase bin ({N_BINS} bins).\n"
        "400 bins raises on this file. The original file is not modified."
    )
    figure = OUT / "folded_outliers.png"
    draw_outliers(payload["t"], payload["mag"], phase, flag, figure, title)
    removed = OUT / "removed.csv"
    write_removed(removed, payload["t"], payload["mag"], flag)
    logger.info("flagged %d of %d", n_flag, len(flag))
    logger.info("wrote %s", figure)
    logger.info("wrote %s", removed)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
