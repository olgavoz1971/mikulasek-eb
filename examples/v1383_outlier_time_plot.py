"""Draw the removed outliers of every V1383 Tau file against time, with the background.

Two panels share the time axis. Top: the SAP magnitude, with the points removed by
``v1383_clean.py`` in red. Bottom: the SAP background flux of the same file, with the
same points marked in red, so that a cluster of outliers can be compared with a
rise of the background. Reads the original VOTable and ``__removed.csv``. Nothing is
fitted or removed.

Run from ``examples``::

    python v1383_outlier_time_plot.py

Writes ``results/v1383_tau/<label>_protocol/<stem>__outliers_time.png``.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from volightcurve import read_lightcurve
from volightcurve.time_reference import (
    resolve_observation_timesys,
    time_offset_to_absolute_jd,
)

import nsv807_s30_tiered as tiered

logger = logging.getLogger(__name__)

DATA = Path("/home/voz/projects/UPJS/lc_shipper/data/V1383_Tau")
RESULTS = Path(__file__).resolve().parents[1] / "results" / "v1383_tau"
MATCH_DAYS = 1e-5


def background(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read the time and the SAP background of one file.

    Args:
        path (Path): VOTable.

    Returns:
        tuple[ndarray, ndarray]: Absolute JD and ``sap_bkg``, in file order.

    Raises:
        ValueError: If the file has no ``sap_bkg`` column.
    """
    volc = read_lightcurve(path, filename=path.name, strict=True)
    if "sap_bkg" not in volc.table.colnames:
        raise ValueError(f"{path.name} has no sap_bkg column")
    ts = resolve_observation_timesys(volc)
    raw = np.asarray(volc.table[volc.get_time_colnames()[0]], dtype=float)
    t_jd = np.array([time_offset_to_absolute_jd(float(v), ts.timeorigin) for v in raw])
    return t_jd, np.asarray(volc.table["sap_bkg"], dtype=float)


def draw(path: Path, folder: Path) -> None:
    """Draw the two-panel figure of one file.

    Args:
        path (Path): VOTable.
        folder (Path): Run folder with the ``__removed.csv``.

    Raises:
        ValueError: If a removed time matches no row of the file.
    """
    payload = tiered.load_sap(path)
    t_all, bkg_all = background(path)
    with open(folder / f"{path.stem}__removed.csv", newline="") as handle:
        removed = np.array([float(r["t_jd"]) for r in csv.DictReader(handle)])
    flag = np.zeros(len(payload["t"]), bool)
    for t in removed:
        hit = np.flatnonzero(np.abs(payload["t"] - t) < MATCH_DAYS)
        if hit.size != 1:
            raise ValueError(f"{path.stem}: removed time {t} matches {hit.size} rows")
        flag[hit[0]] = True
    order = np.argsort(t_all)
    pos = np.searchsorted(t_all[order], payload["t"])
    bkg = bkg_all[order][pos]
    if not np.allclose(t_all[order][pos], payload["t"], atol=MATCH_DAYS):
        raise ValueError(f"{path.stem}: background times do not match")
    fig, (top, bottom) = plt.subplots(
        2, 1, figsize=(16, 9), sharex=True, gridspec_kw={"height_ratios": [3, 2]}
    )
    t0 = float(payload["t"].min())
    x = payload["t"] - t0
    top.scatter(x[~flag], payload["mag"][~flag], s=2, color="0.5", label="kept")
    top.scatter(x[flag], payload["mag"][flag], s=16, color="red", label="removed")
    top.invert_yaxis()
    top.set_ylabel("SAP magnitude")
    top.legend(loc="upper right")
    top.set_title(
        f"{path.stem}: {int(flag.sum())} of {len(flag)} points removed, shown in time"
    )
    bottom.scatter(x[~flag], bkg[~flag], s=2, color="0.5")
    bottom.scatter(x[flag], bkg[flag], s=16, color="red")
    bottom.set_ylabel("SAP background flux")
    bottom.set_xlabel(f"JD - {t0:.4f} (days)")
    fig.tight_layout()
    out = folder / f"{path.stem}__outliers_time.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    logger.info("wrote %s", out)


def main() -> None:
    """Draw every file."""
    tiered.FLUX, tiered.FLUX_ERR = "sap_flux", "sap_flux_err"
    for path in sorted(DATA.glob("*.vot")):
        draw(path, RESULTS / f"{path.stem.split('_')[1]}_protocol")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
