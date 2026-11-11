"""Plot the V1383 Tau magnitude difference of the two maxima.

Every sector is drawn on one figure. The horizontal axis is the MJD of the
primary maximum. The vertical axis is the primary magnitude minus the
secondary magnitude, with its formal error.

Run from ``mikulasek_eb`` with ``MPLBACKEND=Agg``::

    python examples/v1383_max_diff_plot.py

Reads ``results/v1383_tau/timing_max_diff_dat/`` and writes ``dmag_vs_mjd.png``
there.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import mikulasek_eb as eb
import v1383_collect_dat as minima

logger = logging.getLogger(__name__)

FOLDER = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "v1383_tau"
    / "timing_max_diff_dat"
)
STYLE = {
    "s42": ("#0072B2", "Sector 42, 600 s"),
    "s43": ("#E69F00", "Sector 43, 600 s"),
    "s44": ("#009E73", "Sector 44, 600 s"),
    "s70": ("#D55E00", "Sector 70, 200 s"),
    "s71": ("#CC79A7", "Sector 71, 200 s"),
}


def read_sector(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read one difference file.

    Args:
        path (Path): Sector ``.dat`` file.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]: MJD, magnitude difference
        and its error.

    Raises:
        RuntimeError: If a data line does not hold three numbers.
    """
    jd, diff, err = [], [], []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 3:
            raise RuntimeError(f"{path.name}: {line!r}")
        jd.append(float(parts[0]))
        diff.append(float(parts[1]))
        err.append(float(parts[2]))
    if not jd:
        raise RuntimeError(f"{path.name}: no rows")
    mjd = np.asarray(eb.absolute_jd_to_mjd(jd), float)
    return mjd, np.asarray(diff), np.asarray(err)


def main() -> None:
    """Draw every sector on one figure."""
    fig, axis = plt.subplots(figsize=(10, 4.8))
    for lc in minima.LIGHT_CURVES:
        colour, label = STYLE[lc]
        mjd, diff, err = read_sector(
            FOLDER / f"V1383Tau_{lc}_mag_prim_minus_mag_sec.dat"
        )
        axis.errorbar(
            mjd,
            diff,
            yerr=err,
            fmt="o",
            ms=4,
            color=colour,
            ecolor=colour,
            elinewidth=0.6,
            capsize=0,
            linestyle="none",
            label=f"{label}, n={len(mjd)}",
        )
        logger.info(
            "%s: n=%d, median difference %.4f mag, MJD %.1f to %.1f",
            lc,
            len(mjd),
            float(np.median(diff)),
            float(mjd.min()),
            float(mjd.max()),
        )
    axis.axhline(0.0, color="0.5", lw=0.8)
    axis.set_xlabel("MJD of the primary maximum")
    axis.set_ylabel("Primary minus secondary (mag)")
    axis.set_title("V1383 Tau, magnitude difference of the two maxima in one cycle")
    axis.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    dest = FOLDER / "dmag_vs_mjd.png"
    fig.savefig(dest, dpi=140)
    plt.close(fig)
    logger.info("wrote %s", dest)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
