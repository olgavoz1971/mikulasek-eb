"""The curve drawn for a timing window agrees with the fitted event."""

from __future__ import annotations

import numpy as np
import pytest

from mikulasek_eb import EBModel, phase
from mikulasek_eb.timing import effective_window, event_curve, fit_minimum

PERIOD_D = 0.5
M0_D = 0.123
SHIFT_D = 0.002
NOISE_MAG = 0.002
TRUE_PARAMS = dict(A1=0.40, A2=0.30, D1=0.04, G1=1.5, C1=0.0, H2=0.05)


def test_event_curve_follows_the_fitted_event():
    """Residuals of the drawn curve at the window points are at the noise level."""
    rng = np.random.default_rng(7)
    t = np.arange(0.0, 3.0, 4.0 / 1440.0)
    model = EBModel(profile=1, harmonics=(2,), ocon=False)
    ephem = dict(M0=M0_D, P0=PERIOD_D)
    mag = model(phase(t, M0_D + SHIFT_D, PERIOD_D), TRUE_PARAMS)
    mag = mag + 0.01 * (t - 1.0) + rng.normal(0.0, NOISE_MAG, len(t))
    sig = np.full(len(t), NOISE_MAG)
    res = fit_minimum(t, mag, sig, model, TRUE_PARAMS, ephem, 2.0, window=0.15)
    curve = event_curve(t, mag, sig, model, TRUE_PARAMS, ephem, 2.0, res.delta, 0.15)
    inside = curve["in_window"]
    assert inside.sum() == res.n
    drawn = np.interp(t[inside], curve["t_grid"], curve["y_grid"])
    assert np.std(mag[inside] - drawn) < 2.0 * NOISE_MAG
    assert res.delta == pytest.approx(SHIFT_D, abs=5 * res.sigma)


@pytest.mark.parametrize(
    "requested, overlap, expected",
    [(0.4, False, 0.25), (0.4, True, 0.4), (0.7, True, 0.5), (0.1, False, 0.1)],
)
def test_effective_window(requested, overlap, expected):
    """The window is capped as ``virtual_minima`` caps it."""
    assert effective_window(requested, overlap) == expected
