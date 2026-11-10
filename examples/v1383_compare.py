"""Compare the V1383 Tau timings of all methods from the collected ``.dat`` files.

For each sector, kind and method the table gives the number of times, the median
formal sigma (where finite), the rms of O-C about an unweighted straight line
(``fit_ephemeris``, as in the scan), the same rms after removing times more than
600 s from the line (the number removed is given), and the paired difference to the
Mikulasek first-set times (median and MAD-based scatter, on the eclipses timed by
both). Times are matched to the nearest Mikulasek time of the same kind within
0.1 period. Nothing is fitted to the light curves.

Run from ``examples``::

    python v1383_compare.py

Writes ``results/v1383_tau/method_comparison.csv``.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import numpy as np

from mikulasek_eb.timing import fit_ephemeris

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "v1383_tau"
DAT = RESULTS / "timing_dat"
PERIOD = 0.35064
CLIP_S = 600.0
SECTORS = ("s42", "s43", "s44", "s70", "s71")
KINDS = ("primary", "secondary")
METHODS = (
    "mikulasek_window",
    "mikulasek_window_wide",
    "mavka_AP",
    "mavka_WSAP",
    "mavka_WSAP2",
    "mavka_WSL",
)
MAD_TO_SIGMA = 1.4826


def read_times(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read the data lines of a collected ``.dat`` file.

    Args:
        path (Path): File with ``#`` comment lines and columns ``jd err flag``.

    Returns:
        tuple[ndarray, ndarray]: JD and error in days (NaN where the error is None).
    """
    jd, err = [], []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        jd.append(float(parts[0]))
        err.append(np.nan if parts[1] == "None" else float(parts[1]))
    return np.array(jd), np.array(err)


def line_rms(jd: np.ndarray) -> tuple[float, float, int]:
    """Rms of O-C about an unweighted straight line, and after a clip.

    Args:
        jd (ndarray): Times of minima of one kind (days).

    Returns:
        tuple[float, float, int]: rms in s, rms in s after removing times more than
        ``CLIP_S`` from the line, and the number removed.
    """
    jd = np.sort(jd)
    eclipse = np.round((jd - jd[0]) / PERIOD)
    _, res, _ = fit_ephemeris(jd, np.ones_like(jd), eclipse, 1)
    res_s = res * 86400.0
    rms = float(np.sqrt(np.mean(res_s**2)))
    keep = np.abs(res_s) < CLIP_S
    if keep.sum() < 3:
        return rms, float("nan"), int((~keep).sum())
    _, res2, _ = fit_ephemeris(jd[keep], np.ones(keep.sum()), eclipse[keep], 1)
    return rms, float(np.sqrt(np.mean((res2 * 86400.0) ** 2))), int((~keep).sum())


def paired(jd: np.ndarray, ref: np.ndarray) -> tuple[int, float, float]:
    """Difference of times to the nearest reference time.

    Args:
        jd (ndarray): Times of one method.
        ref (ndarray): Reference times of the same kind.

    Returns:
        tuple[int, float, float]: Number matched within 0.1 period, median
        difference in s and MAD-based scatter in s.
    """
    diffs = []
    for t in jd:
        k = int(np.argmin(np.abs(ref - t)))
        if abs(ref[k] - t) < 0.1 * PERIOD:
            diffs.append((t - ref[k]) * 86400.0)
    if not diffs:
        return 0, float("nan"), float("nan")
    d = np.array(diffs)
    mad = MAD_TO_SIGMA * float(np.median(np.abs(d - np.median(d))))
    return len(d), float(np.median(d)), mad


def main() -> None:
    """Write and log the comparison table."""
    rows = []
    for sector in SECTORS:
        for kind in KINDS:
            ref, _ = read_times(DAT / f"V1383Tau_{sector}_mikulasek_window_{kind}.dat")
            for method in METHODS:
                path = DAT / f"V1383Tau_{sector}_{method}_{kind}.dat"
                if not path.exists():
                    continue
                jd, err = read_times(path)
                rms, rms_clip, n_clip = line_rms(jd)
                n_pair, med, mad = paired(jd, ref)
                finite = err[np.isfinite(err)]
                rows.append(
                    {
                        "sector": sector,
                        "kind": kind,
                        "method": method,
                        "n": len(jd),
                        "median_sigma_s": (
                            f"{86400 * np.median(finite):.1f}"
                            if len(finite)
                            else "None"
                        ),
                        "rms_line_s": f"{rms:.1f}",
                        "rms_line_clipped_s": f"{rms_clip:.1f}",
                        "n_clipped": n_clip,
                        "n_paired": n_pair,
                        "diff_median_s": f"{med:.1f}",
                        "diff_mad_s": f"{mad:.1f}",
                    }
                )
    with open(RESULTS / "method_comparison.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for r in rows:
        logger.info(
            "%s %-9s %-22s n%3d sig %6s rms %6s clipped %6s (%d) diff %6s +/- %5s",
            r["sector"],
            r["kind"],
            r["method"],
            r["n"],
            r["median_sigma_s"],
            r["rms_line_s"],
            r["rms_line_clipped_s"],
            r["n_clipped"],
            r["diff_median_s"],
            r["diff_mad_s"],
        )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    main()
