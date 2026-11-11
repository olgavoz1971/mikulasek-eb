"""Exposure-time integration of the model and the cadence helpers."""

from __future__ import annotations

import numpy as np
import pytest

from mikulasek_eb import EBModel, phase
from mikulasek_eb.exposure import (
    estimate_cadence_days,
    estimate_exposure_days,
    suggest_nsub,
)
from mikulasek_eb.fitting import fit_global

PERIOD_D = 0.5
M0_D = 0.1
CADENCE_D = 10.0 / 1440.0
TRUE_PARAMS = dict(A1=0.40, A2=0.25, D1=0.012, G1=1.3, C1=0.0)
TRUTH_NSUB = 61
FIT_NSUB = 7
NOISE_MAG = 0.001


def test_model_integrates_a_boxcar_in_time():
    """The exposure model equals a fine numerical average over the exposure window."""
    exposure = CADENCE_D
    model = EBModel(
        profile=1, harmonics=(), exposure_days=exposure, nsub=201, period_days=PERIOD_D
    )
    plain = EBModel(profile=1, harmonics=())
    theta = np.linspace(-0.05, 0.05, 41)
    offsets = np.linspace(-0.5, 0.5, 4001) * exposure / PERIOD_D
    expected = np.array([plain(th + offsets, TRUE_PARAMS).mean() for th in theta])
    assert np.allclose(model(theta, TRUE_PARAMS), expected, atol=2e-5)
    # The eclipse is shallower and the area is conserved (smearing, not shifting).
    assert (
        model(np.array([0.0]), TRUE_PARAMS)[0] < plain(np.array([0.0]), TRUE_PARAMS)[0]
    )


def test_single_sample_and_no_exposure_are_instantaneous():
    """nsub=1 and exposure_days=None both give the instantaneous model."""
    theta = np.linspace(-0.05, 0.05, 11)
    plain = EBModel(profile=1, harmonics=())(theta, TRUE_PARAMS)
    one = EBModel(
        profile=1, harmonics=(), exposure_days=CADENCE_D, nsub=1, period_days=PERIOD_D
    )
    assert np.allclose(one(theta, TRUE_PARAMS), plain)


def test_exposure_settings_are_validated():
    """Inconsistent exposure arguments raise."""
    with pytest.raises(ValueError):
        EBModel(profile=1, exposure_days=CADENCE_D)
    with pytest.raises(ValueError):
        EBModel(profile=1, period_days=PERIOD_D)
    with pytest.raises(ValueError):
        EBModel(profile=1, exposure_days=-1.0, period_days=PERIOD_D)
    with pytest.raises(ValueError):
        EBModel(profile=1, exposure_days=CADENCE_D, nsub=0, period_days=PERIOD_D)


def test_integration_removes_the_width_bias_of_long_cadence():
    """Long-cadence data of a narrow eclipse: only the integrating fit recovers D."""
    rng = np.random.default_rng(3)
    t = np.arange(0.0, 40.0, CADENCE_D)
    truth_model = EBModel(
        profile=1,
        harmonics=(),
        exposure_days=CADENCE_D,
        nsub=TRUTH_NSUB,
        period_days=PERIOD_D,
    )
    mag = truth_model(phase(t, M0_D, PERIOD_D), TRUE_PARAMS) + rng.normal(
        0, NOISE_MAG, len(t)
    )
    sig = np.full(len(t), NOISE_MAG)
    free = ["A1", "A2", "D1", "G1"]
    start = dict(TRUE_PARAMS, A1=0.3, A2=0.2, D1=0.015, G1=1.0)
    ephem = dict(M0=M0_D, P0=PERIOD_D)
    plain = fit_global(
        t, mag, sig, EBModel(profile=1, harmonics=()), start, free, ephem
    )
    integrating = fit_global(
        t,
        mag,
        sig,
        EBModel(
            profile=1,
            harmonics=(),
            exposure_days=CADENCE_D,
            nsub=FIT_NSUB,
            period_days=PERIOD_D,
        ),
        start,
        free,
        ephem,
    )
    assert integrating.rms < plain.rms
    assert abs(integrating.p["D1"] - TRUE_PARAMS["D1"]) < 5 * integrating.err["D1"]
    assert plain.p["D1"] - TRUE_PARAMS["D1"] > 5 * plain.err["D1"]  # biased wide


def test_cadence_and_exposure_estimates():
    """The helpers return the median step and refuse irregular sampling."""
    t = np.arange(0.0, 1.0, 120.0 / 86400.0)
    t = np.delete(t, slice(200, 260))  # a gap does not change the median step
    assert estimate_cadence_days(t) == pytest.approx(120.0 / 86400.0)
    assert estimate_exposure_days(t, duty_cycle=0.8) == pytest.approx(96.0 / 86400.0)
    with pytest.raises(ValueError):
        estimate_exposure_days(t, duty_cycle=1.5)
    rng = np.random.default_rng(1)
    with pytest.raises(ValueError, match="irregular"):
        estimate_cadence_days(np.cumsum(rng.uniform(0.001, 0.05, 500)))
    with pytest.raises(ValueError):
        estimate_cadence_days(np.array([0.0, 1.0]))


def test_suggest_nsub():
    """The sub-step stays below a quarter of the FWHM, with a floor of 21 samples."""
    assert suggest_nsub(10.0 / 1440.0, 0.05) == 21
    assert suggest_nsub(10.0 / 1440.0, 0.0069) == 21
    assert suggest_nsub(10.0 / 1440.0, 0.0010) == 28
    with pytest.raises(ValueError):
        suggest_nsub(0.0, 0.01)
