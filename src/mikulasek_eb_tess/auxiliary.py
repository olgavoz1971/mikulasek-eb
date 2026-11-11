"""Cleaning and the timing-window scan for a TESS light curve.

The folded-outlier rule and the SAP-background rule are the ones used on the
TESS light curves of this project. ``window_scan`` has no TESS number in it, so
its implementation stays in ``mikulasek_eb.tuning``; it is re-exported here so a
TESS workflow can import the cleaning and the scan from one place.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from volightcurve import read_lightcurve
from volightcurve.time_reference import (
    resolve_observation_timesys,
    time_offset_to_absolute_jd,
)

from mikulasek_eb.tuning import wall_band, window_scan

logger = logging.getLogger(__name__)

# Normal-consistent scale of the median absolute deviation, the value used by the
# folded-outlier rule. Kept as written there, so a point on the cut stays on the cut.
MAD_TO_SIGMA = 1.4826

__all__ = [
    "MAD_TO_SIGMA",
    "flag_folded_outliers",
    "flag_sap_background",
    "read_sap_background",
    "wall_band",
    "window_scan",
]


def flag_folded_outliers(
    t,
    mag,
    period_days: float,
    *,
    n_bins: int = 400,
    n_sigma: float = 5.0,
    min_points: int = 5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Flag points far from the median of their own phase bin.

    The light curve is folded on ``period_days``. The catalogue epoch is not
    used: the phase zero is the first time, and any other zero gives the same
    bins. Each point is compared with the median of its bin. The scatter of
    the bin is ``MAD_TO_SIGMA`` times the median absolute deviation of that bin.

    The defaults (400 bins, 5 sigma, at least 5 points) are the values used for
    the TESS light curves in this work. A sparser file needs fewer bins; pass
    ``n_bins``. A bin that is too thin, or that has no scatter, raises.

    Args:
        t (ndarray): Times (days).
        mag (ndarray): Magnitudes.
        period_days (float): Period (days).
        n_bins (int): Number of phase bins.
        n_sigma (float): Flag a point beyond this many robust sigma.
        min_points (int): Minimum number of points in a bin.

    Returns:
        tuple[ndarray, ndarray, ndarray]: Boolean flag (True for an outlier),
        the residual from the bin median (mag) and the phase of every point.

    Raises:
        ValueError: If the inputs are not usable, a bin has fewer than
            ``min_points`` points, or a bin has zero scatter.
    """
    t = np.asarray(t, dtype=float)
    mag = np.asarray(mag, dtype=float)
    if t.ndim != 1 or mag.shape != t.shape:
        raise ValueError("t and mag must be 1-D and the same length.")
    if int(n_bins) != n_bins or n_bins < 1:
        raise ValueError(f"n_bins must be a positive integer, got {n_bins!r}.")
    if not n_sigma > 0:
        raise ValueError(f"n_sigma must be positive, got {n_sigma!r}.")
    if min_points < 2:
        raise ValueError("min_points must be at least 2, so a bin can have a scatter.")
    if not period_days > 0:
        raise ValueError(f"period_days must be positive, got {period_days!r}.")

    n_bins = int(n_bins)
    phase = ((t - np.min(t)) / period_days) % 1.0
    index = np.minimum((phase * n_bins).astype(int), n_bins - 1)
    residual = np.empty_like(mag)
    local = np.empty(n_bins)
    for k in range(n_bins):
        in_bin = index == k
        n = int(in_bin.sum())
        if n < min_points:
            raise ValueError(
                f"phase bin {k} has {n} points, fewer than {min_points}: "
                "its sigma cannot be measured"
            )
        values = mag[in_bin]
        residual[in_bin] = values - np.median(values)
        local[k] = MAD_TO_SIGMA * np.median(
            np.abs(residual[in_bin] - np.median(residual[in_bin]))
        )
    if not np.all(local > 0):
        raise ValueError("a phase bin has zero scatter")
    return np.abs(residual / local[index]) > n_sigma, residual, phase


def read_sap_background(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Read the time and the SAP background of one TESS file.

    Args:
        path (str | Path): VOTable.

    Returns:
        tuple[ndarray, ndarray]: Absolute JD and ``sap_bkg``, in file order.

    Raises:
        ValueError: If the file has no ``sap_bkg`` column, or not exactly one
            time column.
    """
    path = Path(path)
    volc = read_lightcurve(path, filename=path.name, strict=True)
    if "sap_bkg" not in volc.table.colnames:
        raise ValueError(f"{path.name} has no sap_bkg column")
    time_cols = volc.get_time_colnames()
    if len(time_cols) != 1:
        raise ValueError(f"Need exactly one time column, found {time_cols!r}.")
    ts = resolve_observation_timesys(volc)
    raw = np.asarray(volc.table[time_cols[0]], dtype=float)
    t_jd = np.array([time_offset_to_absolute_jd(float(v), ts.timeorigin) for v in raw])
    return t_jd, np.asarray(volc.table["sap_bkg"], dtype=float)


def flag_sap_background(
    t_jd,
    t_background,
    background,
    *,
    ratio: float,
    match_days: float = 1e-5,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Flag points whose SAP background exceeds ``ratio`` times the median.

    The median is taken on the background of the points in ``t_jd``, so the
    quiet level is that of the light curve being cleaned. The comparison does
    not look at the magnitude: a point on the wall of an eclipse is removed
    only when its background is high. A cut in units of the background scatter
    was tried and rejected, because the quiet background is so stable that any
    multiple of its scatter removes ordinary points. There is no default ratio.
    The V1383 Tau light curves were cleaned at 10.

    Args:
        t_jd (ndarray): Times of the light curve (absolute JD).
        t_background (ndarray): Times of the background samples (absolute JD).
        background (ndarray): ``sap_bkg`` at ``t_background``.
        ratio (float): Remove a point when its background is above this many
            times the median.
        match_days (float): A background time matches a light-curve time within
            this tolerance (days).

    Returns:
        tuple[ndarray, ndarray, float]: Boolean flag, the background of the
        matched points, and the median background.

    Raises:
        ValueError: If the times do not match, the median is not positive, or
            ``ratio`` is not positive.
    """
    t_jd = np.asarray(t_jd, dtype=float)
    t_background = np.asarray(t_background, dtype=float)
    background = np.asarray(background, dtype=float)
    if t_background.ndim != 1 or background.shape != t_background.shape:
        raise ValueError("background times and values must be 1-D and the same length.")
    if t_jd.ndim != 1:
        raise ValueError("t_jd must be 1-D.")
    if not ratio > 0:
        raise ValueError(f"ratio must be positive, got {ratio!r}.")
    order = np.argsort(t_background)
    t_sorted = t_background[order]
    pos = np.searchsorted(t_sorted, t_jd)
    if not np.all((pos >= 0) & (pos < len(t_sorted))):
        raise ValueError("background times do not match the light curve")
    if not np.allclose(t_sorted[pos], t_jd, atol=match_days, equal_nan=False):
        raise ValueError("background times do not match the light curve")
    bkg = background[order][pos]
    centre = float(np.nanmedian(bkg))
    if not np.isfinite(centre) or not centre > 0:
        raise ValueError("the median background is not positive")
    logger.info(
        "background cut: median %.3g, ratio %g, %d of %d points flagged",
        centre,
        ratio,
        int(np.sum(bkg > ratio * centre)),
        len(bkg),
    )
    return bkg > ratio * centre, bkg, centre
