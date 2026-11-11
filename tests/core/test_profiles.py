"""Table 2 profiles (1 to 7) and the Eq. 14 K term."""

from __future__ import annotations

import numpy as np
import pytest

from mikulasek_eb import PROFILES, EBModel, phase
from mikulasek_eb.fitting import fit_global
from mikulasek_eb.model import PROFILE_SHAPE

PHI_GRID = np.linspace(-0.2, 0.2, 801)
WIDTH_D = 0.03
GAMMA = 1.7
CURVATURE_C = -0.2
QUARTIC_K = 0.05
ZERO_TOL = 1e-12


@pytest.mark.parametrize("profile", sorted(PROFILES))
def test_profile_is_one_at_centre_and_symmetric(profile):
    """Every profile is 1 at the centre, symmetric and non-negative."""
    func = PROFILES[profile]
    assert func(np.array([0.0]), WIDTH_D, GAMMA)[0] == pytest.approx(1.0)
    values = func(PHI_GRID, WIDTH_D, GAMMA)
    assert np.allclose(values, values[::-1])
    assert np.all(values >= -ZERO_TOL)


def test_profile_relations_from_table_2():
    """Special cases implied by the formulas of Table 2."""
    x = PHI_GRID
    # Gamma = 1: #2 is exp(1 - cosh x) = #7 and #3 is exp(-x^2/2) = #5.
    assert np.allclose(PROFILES[2](x, WIDTH_D, 1.0), PROFILES[7](x, WIDTH_D))
    assert np.allclose(PROFILES[3](x, WIDTH_D, 1.0), PROFILES[5](x, WIDTH_D))
    # #1 reduces to #2 for C = 0 and K = 0.
    assert np.allclose(PROFILES[1](x, WIDTH_D, GAMMA), PROFILES[2](x, WIDTH_D, GAMMA))


def test_triangle_profile_matches_paper_expression():
    """#6: A(1/2 - |u| + ||u| - 1/2|), u = phi/(2D), is a triangle of half-width D."""
    u = np.abs(PHI_GRID) / (2.0 * WIDTH_D)
    paper = 0.5 - u + np.abs(u - 0.5)
    assert np.allclose(PROFILES[6](PHI_GRID, WIDTH_D), paper)
    assert PROFILES[6](np.array([WIDTH_D]), WIDTH_D)[0] == pytest.approx(0.0)
    assert PROFILES[6](np.array([0.5 * WIDTH_D]), WIDTH_D)[0] == pytest.approx(0.5)


def test_k_term_matches_eq_14_and_is_off_by_default():
    """K enters as (1 + C x^2 + K x^4); K = 0 reproduces Eq. 13 exactly."""
    x = PHI_GRID / WIDTH_D
    base = PROFILES[2](PHI_GRID, WIDTH_D, GAMMA)
    with_k = PROFILES[1](PHI_GRID, WIDTH_D, GAMMA, CURVATURE_C, QUARTIC_K)
    assert np.allclose(with_k, (1 + CURVATURE_C * x**2 + QUARTIC_K * x**4) * base)
    assert np.allclose(
        PROFILES[1](PHI_GRID, WIDTH_D, GAMMA, CURVATURE_C, 0.0),
        PROFILES[1](PHI_GRID, WIDTH_D, GAMMA, CURVATURE_C),
    )


def test_model_uses_only_the_parameters_of_its_profile():
    """Profiles 5 to 7 need no G or C; a missing required key still raises."""
    model = EBModel(profile=5, harmonics=())
    assert model.names() == ["A1", "A2", "D1"]
    flux = model(np.array([0.0, 0.5]), dict(A1=0.3, A2=0.1, D1=WIDTH_D))
    assert flux == pytest.approx([0.3, 0.1])
    with pytest.raises(KeyError):
        EBModel(profile=1, harmonics=())(
            np.array([0.0]), dict(A1=0.3, A2=0.1, D1=WIDTH_D)
        )


def test_k_term_names_and_validation():
    """K is listed only when requested and only for profile 1."""
    assert "K1" in EBModel(profile=1, harmonics=(), kterm=True).names()
    assert "K1" not in EBModel(profile=1, harmonics=()).names()
    with pytest.raises(ValueError):
        EBModel(profile=2, kterm=True)
    with pytest.raises(ValueError):
        EBModel(profile=8)
    assert set(PROFILE_SHAPE) == set(PROFILES)


@pytest.mark.parametrize("profile,truth", [(5, dict(D1=0.03)), (7, dict(D1=0.03))])
def test_one_parameter_profile_fit_recovers_width(profile, truth):
    """A global fit of a Gaussian or exp(1 - cosh) eclipse recovers D within errors."""
    rng = np.random.default_rng(7)
    model = EBModel(profile=profile, harmonics=())
    t = np.arange(0.0, 12.0, 4.0 / 1440.0)
    params = dict(A1=0.35, A2=0.20, **truth)
    mag = model(phase(t, 0.1, 0.9), params) + rng.normal(0.0, 0.002, len(t))
    sig = np.full(len(t), 0.002)
    start = dict(A1=0.3, A2=0.25, D1=0.04)
    fit = fit_global(
        t, mag, sig, model, start, ["A1", "A2", "D1"], dict(M0=0.1, P0=0.9)
    )
    assert not fit.rank_deficient
    assert abs(fit.p["D1"] - truth["D1"]) < 5.0 * fit.err["D1"]


def test_k_term_recovered_on_a_transit():
    """A single transit (A2 = 0) with a quartic bottom: K is recovered within errors."""
    rng = np.random.default_rng(11)
    model = EBModel(profile=1, harmonics=(), kterm=True)
    truth = dict(A1=0.03, A2=0.0, D1=0.012, G1=1.6, C1=-0.15, K1=0.08)
    t = np.arange(0.0, 60.0, 1.0 / 1440.0)
    mag = model(phase(t, 0.2, 3.1), truth) + rng.normal(0.0, 0.0006, len(t))
    sig = np.full(len(t), 0.0006)
    free = ["A1", "D1", "G1", "C1", "K1"]
    start = dict(truth, A1=0.028, D1=0.013, G1=1.4, C1=-0.1, K1=0.0)
    fit = fit_global(
        t, mag, sig, model, start, free, dict(M0=0.2, P0=3.1), free_ephem=("M0",)
    )
    assert fit.success
    assert abs(fit.p["K1"] - truth["K1"]) < 4.0 * fit.err["K1"]
    assert abs(fit.p["D1"] - truth["D1"]) < 4.0 * fit.err["D1"]


def test_pipeline_refuses_profiles_without_gamma():
    """The staged pipeline names the limitation instead of failing deep inside."""
    from mikulasek_eb import FitConfig, run_staged_fit

    t = np.arange(0.0, 3.0, 0.01)
    with pytest.raises(ValueError, match="profiles 1 to 4"):
        run_staged_fit(
            t,
            np.zeros_like(t),
            np.ones_like(t),
            1.0,
            0.0,
            FitConfig(
                profile=5,
                harmonics=(),
                ocon=False,
                baseline_order=0,
                gap_split_days=None,
            ),
        )
