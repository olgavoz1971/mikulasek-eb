"""Flag points far from the folded light curve, and draw them.

The file is ``data/v477_Lyr_sector_26.vot``. The rule is the default one:
400 phase bins, 5 robust sigma, at least 5 points in a bin. The original
file is not modified. Run from the repository root::

    python examples/clean_folded.py

Writes ``examples/out/v477_lyr_clean/folded_outliers.png`` and
``removed.csv``.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import mikulasek_eb_tess as tess

logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "v477_Lyr_sector_26.vot"
OUT = ROOT / "examples" / "out" / "v477_lyr_clean"
FLAGGED_COLOUR = "#D55E00"


def draw_outliers(
    t_jd: np.ndarray,
    mag: np.ndarray,
    phase: np.ndarray,
    flag: np.ndarray,
    path: Path,
    title: str,
) -> None:
    """Draw the folded curve and the time series with the flagged points.

    Args:
        t_jd (ndarray): Absolute JD.
        mag (ndarray): Magnitudes.
        phase (ndarray): Phase from the first time (cycles, 0 to 1).
        flag (ndarray): True for a point beyond the robust bin scatter.
        path (Path): PNG file to write.
        title (str): Figure title.
    """
    kept = ~flag
    fig, axes = plt.subplots(2, 1, figsize=(10, 7))
    axes[0].plot(phase[kept], mag[kept], ".", ms=2, color="0.45", label="kept")
    axes[0].plot(
        phase[flag], mag[flag], ".", ms=4, color=FLAGGED_COLOUR, label="flagged"
    )
    axes[0].set_xlabel("Phase from the first time (cycles)")
    axes[0].set_xlim(0.0, 1.0)
    axes[1].plot(t_jd[kept], mag[kept], ".", ms=2, color="0.45", label="kept")
    axes[1].plot(t_jd[flag], mag[flag], ".", ms=4, color=FLAGGED_COLOUR, label="flagged")
    axes[1].set_xlabel("Time (JD)")
    for ax in axes:
        ax.invert_yaxis()
        ax.set_ylabel("Magnitude")
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def write_removed(path: Path, t_jd: np.ndarray, mag: np.ndarray, flag: np.ndarray) -> None:
    """Write the flagged times. The light curve file is not changed.

    Args:
        path (Path): CSV file to write.
        t_jd (ndarray): Absolute JD.
        mag (ndarray): Magnitudes.
        flag (ndarray): True for a flagged point.
    """
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["jd", "mag"])
        for time, value in zip(t_jd[flag], mag[flag]):
            writer.writerow([f"{time:.8f}", f"{value:.6f}"])


def main() -> None:
    """Flag the default folded rule on V477 Lyr and save the figure."""
    payload = tess.load_tess_votable(DATA)
    flag, _, phase = tess.flag_folded_outliers(
        payload["t"], payload["mag"], payload["period"]
    )
    OUT.mkdir(parents=True, exist_ok=True)
    n_flag = int(flag.sum())
    title = (
        f"{DATA.name}\n"
        f"{n_flag} of {len(flag)} points beyond 5 robust sigma "
        "of their phase bin (400 bins).\n"
        "The original file is not modified."
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
