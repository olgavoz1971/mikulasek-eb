"""Outlier step of the locked protocol for the five V1383 Tau light curves.

For each file in ``lc_shipper/data/V1383_Tau`` this runs ``remove_outliers`` of
``nsv807_s30_tiered.py`` (400 phase bins for the 200 s files and 200 for the 600 s files, 5 sigma, local MAD, removed once, file
period only) on the SAP flux, draws the folded light curve with the flagged
points, and writes the original VOTable without those points
(``nsv807_cleaned_files.clean_one``). The fits and the timing use these files.

Run from ``examples``::

    python v1383_clean.py

Writes ``results/v1383_tau/<label>_protocol/`` (``__removed.csv``, ``__outliers.png``)
and ``results/v1383_tau/cleaned/`` (``<stem>.vot``, ``summary.csv``).
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import nsv807_s30_tiered as tiered
from nsv807_cleaned_files import clean_one

logger = logging.getLogger(__name__)

DATA = Path("/home/voz/projects/UPJS/lc_shipper/data/V1383_Tau")
BINS_200S = 400
BINS_600S = 200
RESULTS = Path(__file__).resolve().parents[1] / "results" / "v1383_tau"


def main() -> None:
    """Flag, draw and write the cleaned files; stop on the first failure."""
    tiered.FLUX, tiered.FLUX_ERR = "sap_flux", "sap_flux_err"
    out_dir = RESULTS / "cleaned"
    rows = []
    for path in sorted(DATA.glob("*.vot")):
        label = path.stem.split("_")[1]
        tiered.OUT = RESULTS / f"{label}_protocol"
        tiered.OUTLIER_BINS = BINS_600S if "_600s" in path.stem else BINS_200S
        payload = tiered.load_sap(path)
        tiered.remove_outliers(payload, path.stem)
        row = clean_one(path, tiered.OUT / f"{path.stem}__removed.csv", out_dir)
        rows.append(row)
    with open(out_dir / "summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    logger.info("wrote %s", out_dir)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
