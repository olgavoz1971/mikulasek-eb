"""Diagnostic: scatter of the folded-curve residual against sap_bkg / median(sap_bkg).

For each V1383 Tau file the residual of every point from the median of its phase
bin (the first step of the folded-curve rule, bins 400 for 200 s and 200 for 600 s)
is plotted against the background divided by the median background of that file.
The robust scatter of the residual in log-spaced ratio bins is drawn over the points.
The ratio where that scatter starts to grow is the candidate for a background cut.
Nothing is removed and no rule is changed.

Run from ``examples``::

    python v1383_bkg_diagnostic.py

Writes ``results/v1383_tau/bkg_ratio_diagnostic.png`` and ``bkg_ratio_scatter.csv``.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import nsv807_s30_tiered as tiered
from v1383_outlier_time_plot import background

logger = logging.getLogger(__name__)

DATA = Path("/home/voz/projects/UPJS/lc_shipper/data/V1383_Tau")
RESULTS = Path(__file__).resolve().parents[1] / "results" / "v1383_tau"
EDGES = np.array([0.0, 0.9, 1.0, 1.1, 1.25, 1.5, 2, 3, 5, 10, 30, 100, 1e4])
MIN_PER_BIN = 5


def ratio_and_residual(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Background ratio and folded-curve residual of every point of one file.

    Args:
        path (Path): VOTable.

    Returns:
        tuple[ndarray, ndarray, ndarray]: ``sap_bkg / median``, residual in mmag and
        the flag of the folded-curve rule (True for a flagged point).

    Raises:
        ValueError: If the background times do not match the light curve.
    """
    tiered.OUTLIER_BINS = 200 if "_600s" in path.stem else 400
    payload = tiered.load_sap(path)
    flag, residual, _ = tiered.fold_outliers(payload)
    t_all, bkg_all = background(path)
    order = np.argsort(t_all)
    pos = np.searchsorted(t_all[order], payload["t"])
    if not np.allclose(t_all[order][pos], payload["t"], atol=1e-5):
        raise ValueError(f"{path.stem}: background times do not match")
    bkg = bkg_all[order][pos]
    return bkg / float(np.nanmedian(bkg)), 1e3 * residual, flag


def main() -> None:
    """Draw the diagnostic for every file and write the binned scatter."""
    tiered.FLUX, tiered.FLUX_ERR = "sap_flux", "sap_flux_err"
    paths = sorted(DATA.glob("*.vot"))
    fig, axes = plt.subplots(len(paths), 1, figsize=(12, 3.4 * len(paths)), sharex=True)
    rows = []
    for ax, path in zip(np.atleast_1d(axes), paths):
        ratio, resid, flag = ratio_and_residual(path)
        label = path.stem.split("_")[1]
        good = ratio > 0
        ax.scatter(ratio[good & ~flag], resid[good & ~flag], s=2, color="0.6")
        ax.scatter(
            ratio[good & flag],
            resid[good & flag],
            s=10,
            color="red",
            label="flagged by the folded-curve rule",
        )
        centres, scatter = [], []
        for lo, hi in zip(EDGES[:-1], EDGES[1:]):
            sel = good & (ratio >= lo) & (ratio < hi)
            if sel.sum() < MIN_PER_BIN:
                continue
            r = resid[sel]
            s = tiered.MAD_TO_SIGMA * float(np.median(np.abs(r - np.median(r))))
            centres.append(np.sqrt(max(lo, 0.8) * min(hi, 1e3)))
            scatter.append(s)
            rows.append(
                {
                    "sector": label,
                    "ratio_lo": lo,
                    "ratio_hi": hi,
                    "n": int(sel.sum()),
                    "robust_sigma_mmag": f"{s:.2f}",
                }
            )
        ax.plot(
            centres,
            scatter,
            "o-",
            color="blue",
            label="robust sigma in ratio bin (mmag)",
        )
        ax.plot(centres, [-s for s in scatter], "o-", color="blue")
        ax.axvline(1.0, color="k", lw=0.5)
        ax.set_xscale("log")
        ax.set_yscale("symlog", linthresh=10)
        ax.set_ylabel("residual (mmag)")
        ax.set_title(f"{path.stem}: {len(ratio)} points")
        ax.legend(loc="upper left", fontsize=8)
    np.atleast_1d(axes)[-1].set_xlabel("sap_bkg / median(sap_bkg) of the file")
    fig.tight_layout()
    fig.savefig(RESULTS / "bkg_ratio_diagnostic.png", dpi=110)
    plt.close(fig)
    with open(RESULTS / "bkg_ratio_scatter.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    logger.info("wrote %s", RESULTS / "bkg_ratio_diagnostic.png")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
