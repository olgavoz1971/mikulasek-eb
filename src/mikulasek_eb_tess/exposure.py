"""Exposure time of a TESS light curve.

TODO (volightcurve): there is no standard way to annotate the exposure time in a
VOTable time series, so ``volightcurve`` has no method that returns it. When the
data carry it, a public method such as ``read_lightcurve(...).get_exposure_time()``
should be added there and used by ``mikulasek_eb_tess.io``. Until then the TESS rule
below is used.
"""

from __future__ import annotations

import numpy as np

from mikulasek_eb.exposure import estimate_cadence_days, estimate_exposure_days

# TESS: the onboard cosmic-ray mitigation removes 20 percent of the exposures at every
# cadence except 20 s (480 s of exposure at a 600 s cadence). A mission property, not a
# pipeline one: it holds for SPOC, QLP and full-frame light curves alike.
TESS_DUTY_CYCLE = 0.8
TESS_FAST_CADENCE_S = 20.0
TESS_FAST_CADENCE_TOL = 0.1


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
