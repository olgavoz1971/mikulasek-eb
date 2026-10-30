"""Exposure time helpers: estimate it from the cadence and choose the sub-sampling.

The exposure integration itself is done by ``EBModel(exposure_days=..., nsub=...,
period_days=...)``: the model is averaged over the exposure, which is what an
integrating detector records (a boxcar in time). The paper does not include it.

TODO (volightcurve): there is no standard way to annotate the exposure time in a
VOTable time series, so ``volightcurve`` has no method that returns it. When the
data carry it, a public method such as ``read_lightcurve(...).get_exposure_time()``
should be added there and used by ``mikulasek_eb.io``. Until then the exposure is a
user parameter, and ``estimate_exposure_days`` gives a cadence-based estimate.
"""

from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger(__name__)

# Fraction of consecutive time steps that must agree with the median cadence.
MIN_REGULAR_FRACTION = 0.5

# Relative tolerance for "agrees with the median cadence".
CADENCE_RELATIVE_TOL = 0.1

# Largest sub-exposure step in units of the eclipse FWHM (midpoint-rule accuracy).
MAX_STEP_OVER_FWHM = 0.25

# Fewest sub-samples. The midpoint rule with n samples captures a fraction 1 - 1/n^2 of the
# leading smearing term (texp^2/24) F'', but with only 5 samples the sharp bottom of a
# profile with G < 1 shows as visible corners in the template (up to 4 mmag on NSV 807 at
# 30 min), so 21 samples are used to keep the template smooth.
MIN_NSUB = 21

# TESS: the onboard cosmic-ray mitigation removes 20 percent of the exposures at every
# cadence except 20 s (480 s of exposure at a 600 s cadence). A mission property, not a
# pipeline one: it holds for SPOC, QLP and full-frame light curves alike.
TESS_DUTY_CYCLE = 0.8
TESS_FAST_CADENCE_S = 20.0
TESS_FAST_CADENCE_TOL = 0.1


def estimate_cadence_days(t: np.ndarray) -> float:
    """Median time step of a light curve.

    Args:
        t (ndarray): Observation times (days).

    Returns:
        float: Cadence (days).

    Raises:
        ValueError: If there are fewer than 3 points, or fewer than half of the steps
            agree with the median (a mixed or irregular cadence: give the exposure
            explicitly instead).
    """
    t = np.sort(np.asarray(t, dtype=float))
    if len(t) < 3:
        raise ValueError("Need at least 3 points to estimate the cadence.")
    steps = np.diff(t)
    cadence = float(np.median(steps))
    if not cadence > 0:
        raise ValueError("Median time step is not positive (duplicate times?).")
    regular = np.mean(np.abs(steps - cadence) < CADENCE_RELATIVE_TOL * cadence)
    if regular < MIN_REGULAR_FRACTION:
        raise ValueError(
            f"Only {regular:.0%} of the time steps agree with the median cadence "
            f"{cadence * 86400:.1f} s: irregular or mixed cadence; set the exposure explicitly."
        )
    return cadence


def estimate_exposure_days(t: np.ndarray, duty_cycle: float = 1.0) -> float:
    """Estimate the exposure time as ``duty_cycle`` times the cadence.

    ``duty_cycle = 1`` is the upper limit (no dead time). Real detectors often integrate
    for less than the cadence (readout, onboard cosmic-ray mitigation; see the mission
    handbook), so pass the mission's value when it is known.

    Args:
        t (ndarray): Observation times (days).
        duty_cycle (float): Exposure divided by cadence, in (0, 1].

    Returns:
        float: Exposure time (days).

    Raises:
        ValueError: If ``duty_cycle`` is outside (0, 1] or the cadence is irregular.
    """
    if not 0.0 < duty_cycle <= 1.0:
        raise ValueError(f"duty_cycle must be in (0, 1], got {duty_cycle!r}.")
    exposure = duty_cycle * estimate_cadence_days(t)
    logger.info(
        "estimated exposure %.1f s (duty cycle %.2f of the cadence)",
        exposure * 86400.0,
        duty_cycle,
    )
    return exposure


def suggest_nsub(exposure_days: float, fwhm_days: float) -> int:
    """Number of sub-exposure samples so the midpoint rule resolves the eclipse.

    The step ``exposure / nsub`` is kept below ``MAX_STEP_OVER_FWHM`` of the FWHM, and
    at least ``MIN_NSUB`` samples are used: even an exposure well below the FWHM biases the
    shape by about (exposure / D)^2 / 24 of the depth, which is not negligible.

    Args:
        exposure_days (float): Exposure time (days).
        fwhm_days (float): Eclipse full width at half depth (days).

    Returns:
        int: Number of sub-samples, at least ``MIN_NSUB``.

    Raises:
        ValueError: If an argument is not positive.
    """
    if not (exposure_days > 0 and fwhm_days > 0):
        raise ValueError("exposure_days and fwhm_days must be positive.")
    return max(MIN_NSUB, int(np.ceil(exposure_days / (MAX_STEP_OVER_FWHM * fwhm_days))))


def tess_exposure_days(t: np.ndarray) -> tuple[float, str]:
    """Exposure time of a TESS light curve from its cadence.

    The effective exposure is 0.8 of the cadence because of the onboard cosmic-ray
    mitigation, except at the 20 s cadence, which has no mitigation (full cadence).

    Args:
        t (ndarray): Observation times (days).

    Returns:
        tuple[float, str]: Exposure (days) and a short text naming the rule used.

    Raises:
        ValueError: If the cadence is irregular (see ``estimate_cadence_days``).
    """
    cadence_s = estimate_cadence_days(t) * 86400.0
    fast = (
        abs(cadence_s - TESS_FAST_CADENCE_S)
        < TESS_FAST_CADENCE_TOL * TESS_FAST_CADENCE_S
    )
    duty = 1.0 if fast else TESS_DUTY_CYCLE
    source = "TESS: full cadence (20 s)" if fast else f"TESS: {duty} x cadence"
    return estimate_exposure_days(t, duty_cycle=duty), source
