"""Optional cleaning of V1383 Tau by the SAP background (not the locked outlier rule).

A point is removed when its ``sap_bkg`` is more than ``RATIO`` times the median
``sap_bkg`` of its own file: ``sap_bkg / median > RATIO``. The quiet level is
measured on the file itself. A sigma of the background was tried and rejected: the
quiet background is so stable that any multiple of its scatter removes harmless
points. The ratio of 10 follows ``v1383_bkg_diagnostic.py``. The rule does not look
at the magnitude, so isolated points on the steep walls of a minimum are not touched. It is an option to compare with the
folded-curve rule of ``v1383_clean.py``, not a replacement for it.

For every file it writes the removed-point list, the folded figure with the removed
points, the time figure with the background below (same time axis), and the cleaned
VOTable. The default ratio is 10; ``--ratio`` changes it.

Run from ``examples``::

    python v1383_clean_background.py [--ratio 10]

Writes ``results/v1383_tau/<label>_bkg/`` and ``results/v1383_tau/cleaned_bkg/``.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

import numpy as np

import nsv807_s30_tiered as tiered
from nsv807_cleaned_files import clean_one
from v1383_outlier_time_plot import background, draw

logger = logging.getLogger(__name__)

DATA = Path("/home/voz/projects/UPJS/lc_shipper/data/V1383_Tau")
RESULTS = Path(__file__).resolve().parents[1] / "results" / "v1383_tau"
RATIO = 10.0


def flag_background(payload: dict, path: Path, ratio: float) -> tuple:
    """Flag the points whose background exceeds ``ratio`` times the file median.

    Args:
        payload (dict): Light curve from ``load_sap``.
        path (Path): VOTable, read again for ``sap_bkg``.
        ratio (float): Threshold as a multiple of the median background.

    Returns:
        tuple: Boolean flag for the points of ``payload``, the background of those
        points and the median background.

    Raises:
        ValueError: If the background times do not match the light curve, or the
            median background is not positive.
    """
    t_all, bkg_all = background(path)
    order = np.argsort(t_all)
    pos = np.searchsorted(t_all[order], payload["t"])
    if not np.allclose(t_all[order][pos], payload["t"], atol=1e-5):
        raise ValueError(f"{path.stem}: background times do not match")
    bkg = bkg_all[order][pos]
    centre = float(np.nanmedian(bkg))
    if not centre > 0:
        raise ValueError(f"{path.stem}: the median background is not positive")
    return bkg > ratio * centre, bkg, centre


def main() -> None:
    """Clean and draw every file by the background rule."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ratio", type=float, default=RATIO)
    ratio = parser.parse_args().ratio
    tiered.FLUX, tiered.FLUX_ERR = "sap_flux", "sap_flux_err"
    out_dir = RESULTS / "cleaned_bkg"
    rows = []
    for path in sorted(DATA.glob("*.vot")):
        folder = RESULTS / f"{path.stem.split('_')[1]}_bkg"
        folder.mkdir(parents=True, exist_ok=True)
        payload = tiered.load_sap(path)
        flag, _, centre = flag_background(payload, path, ratio)
        phase = ((payload["t"] - payload["t"].min()) / payload["period"]) % 1.0
        removed = folder / f"{path.stem}__removed.csv"
        with removed.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["t_jd", "phase", "mag"])
            for i in np.flatnonzero(flag):
                writer.writerow(
                    [
                        f"{payload['t'][i]:.6f}",
                        f"{phase[i]:.5f}",
                        f"{payload['mag'][i]:.5f}",
                    ]
                )
        tiered.outlier_figure(
            payload,
            flag,
            phase,
            folder / f"{path.stem}__outliers.png",
            path.stem,
            rule=f"sap_bkg above {ratio:g} x median ({ratio * centre:.0f} e/s)",
        )
        draw(path, folder)
        row = clean_one(path, removed, out_dir)
        row.update(bkg_median=f"{centre:.1f}", ratio=ratio)
        rows.append(row)
        logger.info(
            "%s: %d of %d removed (median background %.0f e/s, cut %.0f)",
            path.stem,
            int(flag.sum()),
            len(flag),
            centre,
            ratio * centre,
        )
    with open(out_dir / "summary.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
