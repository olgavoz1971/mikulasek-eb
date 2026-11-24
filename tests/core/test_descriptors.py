"""Depth and full width at half depth of the fitted eclipses."""

from __future__ import annotations

import numpy as np
import pytest

from mikulasek_eb import EBModel
from mikulasek_eb.descriptors import eclipse_descriptors


def test_gaussian_profile_has_the_analytic_width():
    """Profile 5 is exp(-x^2/2), so FWHM = 2 D sqrt(2 ln 2) and depth = A."""
    model = EBModel(profile=5, harmonics=())
    p = dict(A1=0.4, A2=0.2, D1=0.02)
    out = eclipse_descriptors(model, p)
    expected = 2 * 0.02 * np.sqrt(2 * np.log(2))
    assert out["eclipse_1"]["depth_mag"] == pytest.approx(0.4, abs=1e-6)
    assert out["eclipse_2"]["depth_mag"] == pytest.approx(0.2, abs=1e-6)
    assert out["eclipse_1"]["fwhm_cycles"] == pytest.approx(expected, rel=1e-4)
    assert out["depth_ratio"] == pytest.approx(0.5, abs=1e-6)
    assert out["fwhm_ratio"] == pytest.approx(1.0, rel=1e-4)


def test_triangle_profile_width_is_the_half_width_parameter():
    """Profile 6 is a triangle of half-width D, so its FWHM equals D."""
    out = eclipse_descriptors(EBModel(profile=6, harmonics=()), dict(A1=0.3, A2=0.3, D1=0.03))
    assert out["eclipse_1"]["fwhm_cycles"] == pytest.approx(0.03, rel=1e-4)


def test_harmonics_and_oconnell_do_not_change_the_descriptors():
    """Proximity terms are not part of the eclipse function."""
    base = dict(A1=0.4, A2=0.25, D1=0.012, G1=1.3, C1=0.0)
    plain = eclipse_descriptors(EBModel(profile=1, harmonics=()), base)
    rich = eclipse_descriptors(
        EBModel(profile=1, harmonics=(1, 2, 3), ocon=True),
        dict(base, H1=0.01, H2=0.02, H3=0.005, S1=0.003),
    )
    assert rich == plain


def test_separate_shape_of_the_second_eclipse_is_used():
    """D2 widens only the secondary."""
    model = EBModel(profile=3, harmonics=())
    p = dict(A1=0.4, A2=0.2, D1=0.01, G1=1.0, D2=0.02, G2=1.0)
    out = eclipse_descriptors(model, p)
    assert out["fwhm_ratio"] == pytest.approx(2.0, rel=1e-4)


def test_second_component_enters_the_depth():
    """A flat-bottom pair of copies is measured as one eclipse curve."""
    model = EBModel(profile=3, harmonics=())
    p = dict(A1=0.3, A2=0.1, D1=0.02, G1=1.0)
    single = eclipse_descriptors(model, p)
    p2 = dict(p, A1b=0.1, D1b=0.005, G1b=1.0)
    double = eclipse_descriptors(model, p2)
    assert double["eclipse_1"]["depth_mag"] == pytest.approx(0.4, abs=1e-6)
    assert double["eclipse_2"] == single["eclipse_2"]


def test_exposure_integration_is_not_part_of_the_descriptors():
    """The model keeps its exposure setting, the descriptors ignore it."""
    p = dict(A1=0.4, A2=0.2, D1=0.012, G1=1.3)
    instant = eclipse_descriptors(EBModel(profile=3, harmonics=()), p)
    boxcar = eclipse_descriptors(
        EBModel(profile=3, harmonics=(), exposure_days=0.01, nsub=21, period_days=0.5), p
    )
    assert boxcar == instant


def test_missing_parameters_are_an_error_not_a_default():
    """An EBModel without p raises instead of guessing."""
    with pytest.raises(ValueError, match="p is required"):
        eclipse_descriptors(EBModel(profile=3, harmonics=()))


def test_eclipse_without_depth_is_an_error():
    """A zero-depth eclipse has no width."""
    model = EBModel(profile=3, harmonics=())
    with pytest.raises(ValueError, match="no positive depth"):
        eclipse_descriptors(model, dict(A1=0.0, A2=0.2, D1=0.01, G1=1.0))


def test_curve_that_never_falls_below_half_is_an_error():
    """A width larger than the grid cannot be measured."""
    model = EBModel(profile=5, harmonics=())
    with pytest.raises(ValueError, match="half depth"):
        eclipse_descriptors(model, dict(A1=0.4, A2=0.2, D1=2.0))
