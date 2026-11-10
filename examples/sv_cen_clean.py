"""Outlier step of the SV Cen protocol: the folded-curve rule on the QLP SAP flux.

For each QLP sector in ``lc_shipper/data/SV_Cen`` this runs ``remove_outliers`` of
``nsv807_s30_tiered.py`` (local MAD per phase bin, 5 sigma, removed once, file period
only, zero of phase at the first time of the sector). The sap_bkg column is not used:
in these files it scatters about zero (sectors 10 to 64) or sits near 80000
(sectors 99 and 100) and is not a background level. The SPOC file of sector 64 is
ignored. The researcher restored the three flagged points of sector 64 (consecutive,
on the curve): they stay in the cleaned file, ``__removed.csv`` and the figures still
mark them, and ``summary.csv`` gives ``flagged`` and ``removed`` apart. The files have no SAP error column, so the errors are not used here.

The number of phase bins follows the cadence, so that a bin holds about 17 to 40
points: 400 bins at 200 s, 200 bins at 600 s (the choices of V1383 Tau) and 70 bins
at 1800 s. A bin with fewer than 5 points stops the run.

For each sector it writes the removed list, the folded figure with the removed
points, an unfolded figure (SAP magnitude against time, removed points in red) and
the cleaned VOTable (``nsv807_cleaned_files.clean_one``).

Run from ``examples``::

    python sv_cen_clean.py

Writes ``results/sv_cen/s<sector>_protocol/`` and ``results/sv_cen/cleaned/``.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import nsv807_s30_tiered as tiered
from nsv807_cleaned_files import clean_one

logger = logging.getLogger(__name__)

DATA = Path("/home/voz/projects/UPJS/lc_shipper/data/SV_Cen")
RESULTS = Path(__file__).resolve().parents[1] / "results" / "sv_cen"
BINS = {"200s": 400, "600s": 200, "1800s": 70}
RESTORED_SECTORS = {64}


def sector_key(path: Path) -> int:
    """Sector number from a file name such as ``TIC..._s64_QLP_200s``.

    Args:
        path (Path): VOTable.

    Returns:
        int: Sector number.
    """
    return int(path.stem.split("_")[1][1:])


def draw_unfolded(payload: dict, removed_csv: Path, out: Path, title: str) -> None:
    """Draw the SAP magnitude against time with the removed points in red.

    Args:
        payload (dict): Light curve from ``load_sap``.
        removed_csv (Path): Removed-point list of this sector.
        out (Path): Figure file.
        title (str): Figure title.

    Raises:
        ValueError: If a removed time matches no single row.
    """
    with open(removed_csv, newline="") as handle:
        removed = [float(r["t_jd"]) for r in csv.DictReader(handle)]
    flag = np.zeros(len(payload["t"]), bool)
    for t in removed:
        hit = np.flatnonzero(np.abs(payload["t"] - t) < 1e-5)
        if hit.size != 1:
            raise ValueError(f"{title}: removed time {t} matches {hit.size} rows")
        flag[hit[0]] = True
    t0 = float(payload["t"].min())
    x = payload["t"] - t0
    fig, ax = plt.subplots(figsize=(16, 5))
    ax.scatter(x[~flag], payload["mag"][~flag], s=2, color="0.5", label="kept")
    ax.scatter(x[flag], payload["mag"][flag], s=16, color="red", label="removed")
    ax.invert_yaxis()
    ax.set_xlabel(f"JD - {t0:.4f} (days)")
    ax.set_ylabel("SAP magnitude")
    ax.set_title(f"{title}: {int(flag.sum())} of {len(flag)} points removed")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def main() -> None:
    """Flag, draw and write the cleaned files; stop on the first failure."""
    tiered.FLUX, tiered.FLUX_ERR = "sap_flux", None
    out_dir = RESULTS / "cleaned"
    rows = []
    for path in sorted(DATA.glob("*QLP*.vot"), key=sector_key):
        tiered.OUT = RESULTS / f"s{sector_key(path)}_protocol"
        tiered.OUTLIER_BINS = BINS[path.stem.rsplit("_", 1)[1]]
        payload = tiered.load_sap(path)
        tiered.remove_outliers(payload, path.stem)
        removed = tiered.OUT / f"{path.stem}__removed.csv"
        draw_unfolded(
            payload, removed, tiered.OUT / f"{path.stem}__outliers_time.png", path.stem
        )
        applied = removed
        if sector_key(path) in RESTORED_SECTORS:
            applied = tiered.OUT / f"{path.stem}__removed_applied.csv"
            applied.write_text("t_jd\n")
        row = clean_one(path, applied, out_dir)
        row["bins"] = tiered.OUTLIER_BINS
        row["flagged"] = len(removed.read_text().splitlines()) - 1
        rows.append(row)
    with open(out_dir / "summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    logger.info("wrote %s", out_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
