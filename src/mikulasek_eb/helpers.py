"""Time conversion, template completion and fit reporting helpers."""

from __future__ import annotations

import logging

import numpy as np
from astropy.time import Time

from .model import EBModel, phase
from .fitting import FitResult, _Nuisance
from .tuning import binned_residuals

# Phase half-width for out-of-eclipse residual clip (cycles).
OUT_OF_ECLIPSE_PHASE_HALF = 0.12


def absolute_jd_to_mjd(t_jd: np.ndarray | float) -> np.ndarray:
    """Convert absolute Julian Date to Modified Julian Date for LM fitting.

    Args:
        t_jd (ndarray or float): Absolute JD (J0 = 0).

    Returns:
        ndarray: MJD (same shape as input).
    """
    return np.asarray(Time(t_jd, format="jd").mjd, dtype=float)


def complete_template(p: dict, model: EBModel) -> dict:
    """Fill every key the model can use so staged fits never drop a term.

    Args:
        p (dict): Partial template dict.
        model (EBModel): Model that defines harmonics and O'Connell.

    Returns:
        dict: Copy with C, optional D2/G2/C2, H_k and S1 present.
    """
    q = dict(p)
    q.setdefault("C1", 0.0)
    q.setdefault("D2", q["D1"])
    q.setdefault("G2", q["G1"])
    q.setdefault("C2", q.get("C1", 0.0))
    for k in model.harmonics:
        q.setdefault(f"H{k}", 0.0)
    if model.ocon:
        q.setdefault("S1", 0.0)
    return q


def out_of_eclipse_phase_mask(
    theta: np.ndarray, half: float = OUT_OF_ECLIPSE_PHASE_HALF
) -> np.ndarray:
    """Boolean mask of points away from both eclipse centres.

    Args:
        theta (ndarray): Fitted phase (cycles).
        half (float): Half-width around phase 0 and 0.5 (cycles).

    Returns:
        ndarray: True out of eclipse.
    """
    ph = theta - np.round(theta)
    d2 = ph - 0.5
    d2 = d2 - np.round(d2)
    return (np.abs(ph) > half) & (np.abs(d2) > half)


def primary_ripple_mmag(res: np.ndarray, th: np.ndarray) -> float:
    """Peak absolute binned residual around the primary (mmag).

    Args:
        res (ndarray): Residuals (mag).
        th (ndarray): Phase (cycles).

    Returns:
        float: Max |binned mean| in mmag; NaN if empty.
    """
    b = binned_residuals(res, th)
    vals = np.asarray(b["primary"], dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float("nan")
    return float(np.max(np.abs(vals)))


def report_fit(tag: str, r: FitResult, noise: float) -> None:
    """Log one staged-fit attempt: chi-squared, rms vs noise floor, parameters.

    Args:
        tag (str): Attempt label.
        r (FitResult): Global-fit result.
        noise (float): Robust white-noise floor (mag).
    """
    log = logging.getLogger(__name__)
    log.info(
        "[%s] success=%s nfev=%d chi2=%.1f rchi2=%.2f rms=%.5f mag rms/noise=%.2f",
        tag,
        r.success,
        r.nfev,
        r.chi2,
        r.rchi2,
        r.rms,
        r.rms / noise,
    )
    vals = {**r.p, **r.ephem}
    for n in r.names:
        log.info("  %s = %.6f +- %.6f", n, vals[n], r.err[n])


def predicted_magnitude(t, mag, sig, model, p, ephem, seg, order):
    """Return the global-fit prediction, including per-segment polynomial baselines.

    Args:
        t (ndarray): Time in days (MJD for fits in this project).
        mag (ndarray): Observed magnitudes.
        sig (ndarray): Magnitude uncertainties.
        model (EBModel): Phenomenological eclipse model.
        p (dict): Fitted template parameters.
        ephem (dict): Fitted ephemeris (M0, P0, optional Pdot).
        seg (ndarray): Segment indices used in ``fit_global``.
        order (int): Polynomial degree used in ``fit_global``.

    Returns:
        tuple[ndarray, ndarray, ndarray]: Full prediction ``yhat``, template
        ``F`` (no baselines), and fitted phase ``theta``.
    """
    w = 1.0 / sig
    theta = phase(t, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0))
    f_template = model(theta, p)
    nui = _Nuisance(t, w, seg, order)
    yhat = mag - nui.project(w * (mag - f_template)) / w
    return yhat, f_template, theta
