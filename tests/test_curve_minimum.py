"""Definitions of one minimum read off a curve."""

from __future__ import annotations

import numpy as np
import pytest

from mikulasek_eb import measure_curve_minimum
from mikulasek_eb.errors import InvalidLightCurve


def _asymmetric():
    """A flux dip with a straight left wing twice as wide as the right one.

    Continuum 1, minimum 0 at time 0. The bisector at fraction f is at -f/2, and the
    centroid of the deficit is at -1/3.
    """
    left = np.linspace(-2.0, 0.0, 2001)
    right = np.linspace(0.0, 1.0, 1001)[1:]
    t = np.concatenate([left, right])
    y = np.concatenate([-0.5 * left, right])
    return t, y


def test_symmetric_triangle_puts_every_definition_at_the_centre():
    """A symmetric dip has the same time for every definition."""
    t = np.linspace(-1.0, 1.0, 4001)
    y = np.abs(t)
    m = measure_curve_minimum(t, y, 0.0, 1.0, continuum=1.0)
    assert m.t_lowest == pytest.approx(0.0, abs=1e-6)
    assert m.t_bisector_median == pytest.approx(0.0, abs=1e-4)
    assert all(v == pytest.approx(0.0, abs=1e-4) for v in m.t_bisector.values())
    assert m.t_centroid == pytest.approx(0.0, abs=1e-3)
    assert m.depth == pytest.approx(1.0)


def test_asymmetric_wings_move_the_bisector_and_the_centroid():
    """The bisector at fraction f is -f/2 and the centroid is -1/3."""
    t, y = _asymmetric()
    m = measure_curve_minimum(t, y, 0.0, 2.0, continuum=1.0)
    for frac, time in m.t_bisector.items():
        assert time == pytest.approx(-0.5 * frac, abs=1e-3)
    assert m.t_bisector_median == pytest.approx(-0.25, abs=1e-3)
    assert m.t_centroid == pytest.approx(-1.0 / 3.0, abs=1e-2)
    assert m.t_lowest == pytest.approx(0.0, abs=1e-6)


def test_magnitude_curve_gives_the_same_times_as_its_flux():
    """Turning the curve over and switching the kind does not move the times."""
    t, y = _asymmetric()
    flux = measure_curve_minimum(t, y, 0.0, 2.0, continuum=1.0)
    mag = measure_curve_minimum(t, -y, 0.0, 2.0, kind="mag", continuum=-1.0)
    assert mag.t_bisector_median == pytest.approx(flux.t_bisector_median, abs=1e-6)
    assert mag.t_centroid == pytest.approx(flux.t_centroid, abs=1e-6)
    assert mag.depth == pytest.approx(flux.depth)


def test_a_wing_that_crosses_a_level_twice_is_rejected():
    """A bump on one wing makes that level ambiguous."""
    t = np.linspace(-1.0, 1.0, 401)
    y = np.abs(t)
    y[(t > -0.7) & (t < -0.4)] -= 0.5
    with pytest.raises(InvalidLightCurve, match="more than once"):
        measure_curve_minimum(t, y, 0.0, 1.0, continuum=1.0)


def test_a_level_above_the_wings_is_rejected():
    """A fraction whose level the curve never reaches is an error."""
    t = np.linspace(-1.0, 1.0, 401)
    with pytest.raises(InvalidLightCurve, match="not reached"):
        measure_curve_minimum(t, np.abs(t), 0.0, 1.0, continuum=3.0, fractions=(0.9,))


def test_folded_window_uses_the_period():
    """Points just after time zero belong to the minimum at the end of the series."""
    period = 4.0
    t = np.linspace(0.0, period, 4001, endpoint=False)
    offset = np.mod(t - period + period / 2, period) - period / 2
    y = np.abs(offset)
    m = measure_curve_minimum(t, y, period, 1.0, continuum=1.0, period=period)
    assert m.t_lowest == pytest.approx(period, abs=1e-3)
    assert m.t_bisector_median == pytest.approx(period, abs=1e-3)
