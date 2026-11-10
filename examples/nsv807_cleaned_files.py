"""Write the NSV 807 light curves without the points the tiered-fit protocol removed.

For each light curve the points listed in ``<run>/<stem>__removed.csv`` (the outlier
step of the locked protocol) are dropped from the original VOTable, and the rest is
written unchanged by ``volightcurve.write_lightcurve``. Nothing else is altered.
Sector 4 has no removed list, because the outlier step does not run on it, so it is
not written.

Run from ``examples``::

    python nsv807_cleaned_files.py

Writes ``results/nsv807/cleaned/<stem>.vot`` and ``cleaned/summary.csv``.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import numpy as np
from volightcurve import read_lightcurve, write_lightcurve
from volightcurve.time_reference import (
    resolve_observation_timesys,
    time_offset_to_absolute_jd,
)

logger = logging.getLogger(__name__)

DATA = Path("/home/voz/projects/UPJS/lc_shipper/data")
RESULTS = Path(__file__).resolve().parents[1] / "results" / "nsv807"
RUNS = {
    "TIC35119266_s30_SPOC_120s": "s30_protocol",
    "TIC35119266_s30_TESS-SPOC_600s": "s30_600s_protocol",
    "TIC35119266_s97_QLP_200s": "s97_protocol",
}
MATCH_DAYS = 1e-5


def clean_one(data_file: Path, removed_csv: Path, out_dir: Path) -> dict:
    """Drop the removed points of one light curve and write the VOTable.

    Args:
        data_file (Path): Original VOTable.
        removed_csv (Path): Removed-point list written by the outlier step.
        out_dir (Path): Output folder.

    Returns:
        dict: ``file``, ``rows_in``, ``removed``, ``rows_out``.

    Raises:
        ValueError: If a removed time does not match exactly one row.
    """
    stem = data_file.stem
    volc = read_lightcurve(data_file, filename=data_file.name, strict=True)
    ts = resolve_observation_timesys(volc)
    raw = np.asarray(volc.table[volc.get_time_colnames()[0]], dtype=float)
    t_jd = np.array([time_offset_to_absolute_jd(float(v), ts.timeorigin) for v in raw])
    with open(removed_csv, newline="") as handle:
        removed = [float(r["t_jd"]) for r in csv.DictReader(handle)]
    drop = []
    for t in removed:
        hit = np.flatnonzero(np.abs(t_jd - t) < MATCH_DAYS)
        if hit.size != 1:
            raise ValueError(f"{stem}: removed time {t} matches {hit.size} rows")
        drop.append(int(hit[0]))
    if len(set(drop)) != len(drop):
        raise ValueError(f"{stem}: a row is listed twice")
    rows_in = len(volc.table)
    volc.table.remove_rows(sorted(drop))
    out_dir.mkdir(parents=True, exist_ok=True)
    write_lightcurve(volc, "votable_binary", destination=out_dir / f"{stem}.vot")
    logger.info("%s: %d rows, %d removed", stem, rows_in, len(drop))
    return {
        "file": f"{stem}.vot",
        "rows_in": rows_in,
        "removed": len(drop),
        "rows_out": rows_in - len(drop),
    }


def main() -> None:
    """Write the cleaned files and a summary table."""
    out_dir = RESULTS / "cleaned"
    rows = [
        clean_one(
            DATA / f"{stem}.vot",
            RESULTS / run / f"{stem}__removed.csv",
            out_dir,
        )
        for stem, run in RUNS.items()
    ]
    with open(out_dir / "summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
