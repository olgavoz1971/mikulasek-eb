"""Save the two standard figures of a fitted light curve.

``fit.png`` is the time series, the fold and the two minima.
``trend.png`` is the per-segment baseline next to the template.
"""

from __future__ import annotations

import logging
from pathlib import Path

import mikulasek_eb as eb
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.plotting import plot_global_fit, plot_instrumental_trend

logger = logging.getLogger(__name__)


def save_figures(
    folder: Path,
    title: str,
    *,
    t_mjd,
    mag,
    sig,
    model,
    p_single: dict,
    p_adopted: dict,
    ephem: dict,
    seg,
    order: int,
    timescale: str | None,
) -> None:
    """Write ``fit.png`` and ``trend.png`` for one adopted fit.

    Args:
        folder (Path): Destination. Created if needed.
        title (str): Figure title.
        t_mjd (ndarray): Times the fit used (MJD).
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        model: Fitted ``EBModel``.
        p_single (dict): Single-component parameters. Drawn in grey when they
            differ from ``p_adopted``.
        p_adopted (dict): Adopted template parameters.
        ephem (dict): Fitted ephemeris (M0 in MJD).
        seg (ndarray | None): Baseline segment of each point.
        order (int): Baseline polynomial degree.
        timescale (str | None): TIMESYS timescale, for the axis label.
    """
    folder.mkdir(parents=True, exist_ok=True)
    yhat, template, theta = predicted_magnitude(
        t_mjd, mag, sig, model, p_adopted, ephem, seg, order
    )
    t_jd = eb.mjd_to_jd(t_mjd)
    fit_path = folder / "fit.png"
    trend_path = folder / "trend.png"
    plot_global_fit(
        t_jd,
        mag,
        sig,
        yhat,
        template,
        theta,
        model,
        p_single,
        p_adopted,
        seg,
        title,
        timescale=timescale,
        save_path=str(fit_path),
    )
    plot_instrumental_trend(t_jd, mag, yhat, template, trend_path, title)
    logger.info("wrote %s", fit_path)
    logger.info("wrote %s", trend_path)
