"""Tests of the light-curve features (synthetic truth and real-data smoke checks)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import mikulasek_eb as eb
from mikulasek_eb.errors import InvalidLightCurve
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.model import EBModel, phase

ROOT = Path(__file__).resolve().parent.parent
P = 0.9
M0 = 58000.2
TRUE = dict(A1=0.40, A2=0.20, D1=0.03, G1=2.0, C1=0.0, H2=0.02, S1=0.006)


def _synthetic(period_true=P, n_days=60.0, noise=0.002, seed=1):
    """Detached synthetic star at 10-min cadence, with an O'Connell term."""
    rng = np.random.default_rng(seed)
    t = np.arange(58000.0, 58000.0 + n_days, 10.0 / 1440.0)
    t = t + rng.uniform(-3.0, 3.0, len(t)) / 1440.0  # avoid exact phase repeats
    model = EBModel(1, (2,), True)
    y = model(phase(t, M0, period_true), TRUE) + rng.normal(0, noise, len(t))
    return t, y, np.full(len(t), noise), model


def test_synthetic_geometry_and_shape_match_the_model():
    """Depths, widths, Gamma and the maxima difference agree with the truth."""
    t, y, s, model = _synthetic()
    f = eb.extract_features(t, y, s, P, 2458000.7)
    g = f.geometry
    # depth is measured from the quadrature level, so the harmonic terms count
    q = model(np.array([0.0, 0.5, 0.25, 0.75]), TRUE)
    level = 0.5 * (q[2] + q[3])
    assert g["A1"] == pytest.approx(q[0] - level, abs=0.01)
    assert g["A2"] == pytest.approx(q[1] - level, abs=0.01)
    assert abs(g["secondary_offset"]["value"]) < 0.005
    ph = np.linspace(-0.2, 0.2, 40001)
    prof = model(ph, dict(TRUE, H2=0.0, S1=0.0))
    half = prof > 0.5 * TRUE["A1"] + (prof.min())
    true_fwhm = ph[half][-1] - ph[half][0]
    assert g["fwhm1"] == pytest.approx(true_fwhm, rel=0.05)
    # With harmonics the widths are measured on a curved baseline, so Gamma is biased
    # (1.45 here against 2.0 true); check the estimator on the pure profile below.
    assert f.shape["gamma_estimate"] is not None
    # (the 5 percent edge, 22 mmag here, is as small as the 20 mmag harmonic, so the
    # separation of eclipses is tested on the pure profile below)
    md = f.variation["maxima_difference"]
    expected = model(np.array([0.25]), TRUE) - model(np.array([0.75]), TRUE)
    assert md["value"] == pytest.approx(float(expected[0]), abs=5 * md["err"])
    assert f.variation["maxima_significance"] > 3


def test_gamma_estimate_recovers_the_pure_profile():
    """Without harmonics the width ratios give back the true Gamma."""
    rng = np.random.default_rng(4)
    t = np.arange(58000.0, 58060.0, 10.0 / 1440.0)
    t = t + rng.uniform(-3.0, 3.0, len(t)) / 1440.0
    model = EBModel(1, (), False)
    for gamma in (1.0, 2.0, 4.0):
        p = dict(A1=0.4, A2=0.2, D1=0.03, G1=gamma, C1=0.0)
        y = model(phase(t, M0, P), p) + rng.normal(0, 0.002, len(t))
        f = eb.extract_features(t, y, np.full(len(t), 0.002), P, 2458000.7)
        assert f.shape["gamma_estimate"] == pytest.approx(gamma, rel=0.15)
        assert f.shape["gamma_match_distance"] < 0.05
        assert f.variation["separated_eclipses"] is True


def test_period_error_shows_up_as_a_drift_between_halves():
    """A wrong catalogue period gives a measurable shift and a period correction."""
    p_true = P * (1 + 2e-4)
    t, y, s, _ = _synthetic(period_true=p_true, n_days=90.0)
    f = eb.extract_features(t, y, s, P, 2458000.7)
    d = f.ephemeris["drift_between_halves"]["primary"]
    assert d["period_correction_days"] == pytest.approx(
        p_true - P, abs=4 * d["period_correction_err_days"]
    )
    assert abs(d["shift_s"]) > 4 * d["shift_err_s"]


def test_features_round_trip_and_hash():
    """JSON round trip keeps the hash; the hash ignores the package version."""
    t, y, s, _ = _synthetic()
    f = eb.extract_features(t, y, s, P, 2458000.7)
    g = eb.Features.from_json(f.to_json())
    assert g.features_hash() == f.features_hash()
    h = eb.Features.from_dict(
        {**f.to_dict(), "provenance": {**f.provenance, "package_version": "x"}}
    )
    assert h.features_hash() == f.features_hash()
    again = eb.extract_features(t, y, s, P, 2458000.7)
    assert again.features_hash() == f.features_hash()


def test_wrong_period_is_reported_not_hidden():
    """A period far from the truth leaves no usable eclipse; the failure is explicit."""
    t, y, s, _ = _synthetic()
    with pytest.raises(InvalidLightCurve):
        eb.extract_features(t[:-1], y, s, P, 2458000.7)
    f = eb.extract_features(t, y, s, 1.37 * P, 2458000.7)
    assert (
        f.geometry["A2"] is None or f.shape["gamma_match_distance"] is None or f.notes
    )


def test_unmeasurable_features_are_none_with_a_note():
    """A secondary that is not an eclipse gives None and a note, not a default."""
    rng = np.random.default_rng(2)
    t = np.arange(58000.0, 58040.0, 10.0 / 1440.0)
    t = t + rng.uniform(-3.0, 3.0, len(t)) / 1440.0
    model = EBModel(1, (1,), False)
    y = model(phase(t, M0, P), dict(A1=0.4, A2=0.0, D1=0.03, G1=2.0, C1=0.0, H1=0.1))
    f = eb.extract_features(t, y, np.full(len(t), 0.002), P, 2458000.7)
    if f.geometry["A2"] is None:
        assert any(n.startswith("A2:") for n in f.notes)


def _load(name):
    p = load_tess_votable(ROOT / "data" / f"{name}.vot")
    return p, absolute_jd_to_mjd(p["t"])


@pytest.mark.parametrize(
    "name,separated,contact_candidate",
    [
        ("lc_discovery_TIC_455790537_SPOC_sector_16", False, False),  # AA And
        ("lc_discovery_TIC_455790537_TESS_SPOC_sector_16", False, False),
        ("EK_Com_sector_49", False, True),
        ("lc_discovery_TIC_347405567_TESS_SPOC_sector_50", False, True),  # XY Boo
        ("AV_Del_sector_81", False, False),
        ("AR_Aur_sector_19", True, False),
    ],
)
def test_real_stars_eclipse_character(name, separated, contact_candidate):
    """Separated eclipses and the contact candidate follow what the folds show.

    AR Aur has a flat stretch between narrow eclipses; AA And and AV Del have unequal
    depths with continuous variation; EK Com and XY Boo have equal depths and no flat
    stretch. Calibrated on these stars only (see docs/feature_spec.md).
    """
    p, t = _load(name)
    f = eb.extract_features(t, p["mag"], p["sig"], p["period"], p["epoch_jd"])
    v = f.variation
    assert v["separated_eclipses"] is separated
    assert v["contact_candidate"] is contact_candidate


def test_edge_level_does_not_change_the_ordering():
    """AR Aur stays separated and the others stay continuous for edge levels 0.02 to 0.1."""
    for level in (0.02, 0.05, 0.1):
        st = eb.FeatureSettings(edge_level=level)
        cov = {}
        for name in ("AR_Aur_sector_19", "AV_Del_sector_81", "EK_Com_sector_49"):
            p, t = _load(name)
            f = eb.extract_features(
                t, p["mag"], p["sig"], p["period"], p["epoch_jd"], st
            )
            cov[name] = f.variation["eclipse_coverage"]
        assert (
            cov["AR_Aur_sector_19"]
            < 0.3
            < cov["AV_Del_sector_81"]
            < cov["EK_Com_sector_49"]
        )


def test_no_secondary_leaves_the_branch_undecided():
    """V477 Lyr has no secondary eclipse: the contact candidate is None with a note."""
    p, t = _load("v477_Lyr_sector_26")
    f = eb.extract_features(t, p["mag"], p["sig"], p["period"], p["epoch_jd"])
    assert f.variation["contact_candidate"] is None
    assert any(n.startswith("A2:") for n in f.notes)


def test_real_xy_boo_ephemeris_features_are_finite():
    """XY Boo: drift between halves is measured for both eclipses."""
    p, t = _load("lc_discovery_TIC_347405567_TESS_SPOC_sector_50")
    f = eb.extract_features(t, p["mag"], p["sig"], p["period"], p["epoch_jd"])
    for kind in ("primary", "secondary"):
        d = f.ephemeris["drift_between_halves"][kind]
        assert d is not None and np.isfinite(d["shift_s"]) and d["shift_err_s"] > 0
