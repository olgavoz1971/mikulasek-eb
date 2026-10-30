"""Tests of borrowing a model between light curves (ParameterPlan, fit_borrowed)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import mikulasek_eb as eb
from mikulasek_eb.config import DataSetup, ModelSpec
from mikulasek_eb.errors import ConfigError
from mikulasek_eb.exposure import tess_exposure_days
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.model import EBModel, phase

ROOT = Path(__file__).resolve().parent.parent
TRUE = dict(A1=0.40, A2=0.20, D1=0.03, G1=1.5, C1=0.0, H2=0.02)
EPHEM = dict(M0=58000.1, P0=0.9, Pdot=0.0)


def _reference(params: dict) -> eb.ModelSolution:
    """Hand-built reference solution for a synthetic star."""
    config = eb.Config(
        ModelSpec(profile=1, harmonics=(2,), ocon=False),
        DataSetup(baseline_order=0, gap_split_days=None),
    )
    return eb.ModelSolution(
        config=config,
        params=params,
        errors={},
        ephemeris=EPHEM,
        stats={},
        exposure=None,
        provenance=dict(
            data_id="synthetic",
            epoch_jd=2458000.6,
            period_input_days=0.9,
            package_version="test",
        ),
    )


def test_plan_validation_and_round_trip():
    """Unknown or repeated groups raise ConfigError; plans round-trip."""
    with pytest.raises(ConfigError):
        eb.ParameterPlan(("shapes",)).validate()
    with pytest.raises(ConfigError):
        eb.ParameterPlan(("shape", "shape")).validate()
    plan = eb.ParameterPlan(("amplitudes", "period"))
    assert eb.ParameterPlan.from_dict(plan.to_dict()) == plan


def test_plan_selects_names_by_group():
    """Group names map to the expected parameters, including second components."""
    params = dict(A1=1, A2=1, D1=1, G1=1, C1=0, D2=1, A1b=1, D1b=1, G1b=1, H2=1, S1=1)
    names, eph = eb.ParameterPlan(("shape", "zero_phase")).free_names(params)
    assert sorted(names) == ["C1", "D1", "D1b", "D2", "G1", "G1b"]
    assert eph == ["M0"]
    names, _ = eb.ParameterPlan(("amplitudes", "harmonics")).free_names(params)
    assert sorted(names) == ["A1", "A1b", "A2", "H2", "S1"]


def test_synthetic_borrowing_recovers_amplitudes_and_keeps_fixed_values():
    """Diluted, low-cadence data: free amplitudes are recovered, shape stays fixed."""
    rng = np.random.default_rng(3)
    t = np.arange(58000.0, 58020.0, 30.0 / 1440.0)
    expo = 30.0 / 1440.0
    model = EBModel(1, (2,), False, exposure_days=expo, nsub=9, period_days=EPHEM["P0"])
    dilute = dict(TRUE, A1=0.9 * TRUE["A1"], A2=0.9 * TRUE["A2"])
    sig = np.full(len(t), 0.002)
    y = model(phase(t, EPHEM["M0"] + 0.003, EPHEM["P0"]), dilute)
    y = y + rng.normal(0.0, sig)
    ref = _reference(dict(TRUE))
    out = eb.fit_borrowed(
        t,
        y,
        sig,
        ref,
        DataSetup(baseline_order=0, gap_split_days=None, exposure_days=expo),
        eb.ParameterPlan(("amplitudes", "zero_phase")),
    )
    s = out.solution
    assert s.p["A1"] == pytest.approx(dilute["A1"], abs=5 * s.errors["A1"])
    assert s.p["A2"] == pytest.approx(dilute["A2"], abs=5 * s.errors["A2"])
    assert s.p["D1"] == TRUE["D1"] and s.p["G1"] == TRUE["G1"]
    assert s.ephem["M0"] == pytest.approx(EPHEM["M0"] + 0.003, abs=5 * s.errors["M0"])
    assert s.ephem["P0"] == EPHEM["P0"]
    assert s.provenance["borrowed_from"] == ref.solution_hash()
    assert s.provenance["plan"] == {"free": ["amplitudes", "zero_phase"]}
    assert s.exposure["days"] == expo
    assert out.pipeline is None
    assert eb.ModelSolution.from_json(s.to_json()).solution_hash() == s.solution_hash()


def test_empty_plan_is_an_error():
    """A plan that frees nothing present in the reference is rejected."""
    t = np.arange(58000.0, 58010.0, 0.02)
    ref = _reference(dict(TRUE))
    with pytest.raises(ConfigError):
        eb.fit_borrowed(
            t,
            np.zeros_like(t),
            np.full_like(t, 0.01),
            ref,
            DataSetup(baseline_order=0, gap_split_days=None),
            eb.ParameterPlan(()),
        )


# ---------------------------------------------------------------- real data (AA And)
def _load(name):
    p = load_tess_votable(ROOT / "data" / f"lc_discovery_TIC_455790537_{name}.vot")
    return p, absolute_jd_to_mjd(p["t"])


@pytest.fixture(scope="module")
def aa_and_reference():
    """AA And 2-min sector 16 reference solution (same settings as the per-star script)."""
    import fit_tic455790537_spoc_s16 as star

    p, t = _load("SPOC_sector_16")
    out = eb.fit(
        t,
        p["mag"],
        p["sig"],
        p["period"],
        p["epoch_jd"],
        eb.Config.from_fit_config(star.FIT_CONFIG),
        data_id="s16_2min",
    )
    return out.solution


def _borrow(ref, name, free):
    p, t = _load(name)
    data = DataSetup(
        baseline_order=3,
        gap_split_days=0.5,
        exposure_days=tess_exposure_days(t)[0],
    )
    return (
        t,
        p,
        eb.fit_borrowed(
            t, p["mag"], p["sig"], ref, data, eb.ParameterPlan(free), data_id=name
        ),
    )


def test_aa_and_30min_borrowed_fit_and_timings_agree_with_2min(aa_and_reference):
    """Borrowed 30-min fit fits well, and its timings match the 2-min timings."""
    ref = aa_and_reference
    t, p, out = _borrow(
        ref, "TESS_SPOC_sector_16", ("amplitudes", "harmonics", "zero_phase")
    )
    assert out.solution.stats["rms_mag"] < 0.0020
    recipe = eb.TimingRecipe(require_contact_walls=False, clip_outliers=False)
    slow = eb.time(out.t_mjd, out.mag, out.sig, out.solution, recipe)
    p2, t2 = _load("SPOC_sector_16")
    fast_out = eb.fit_borrowed(
        t2,
        p2["mag"],
        p2["sig"],
        ref,
        DataSetup(baseline_order=3, gap_split_days=0.5),
        eb.ParameterPlan(("amplitudes", "harmonics", "zero_phase")),
    )
    fast = eb.time(
        fast_out.t_mjd, fast_out.mag, fast_out.sig, fast_out.solution, recipe
    )
    for kind in ("primary", "secondary"):
        a, b = fast.timings[kind], slow.timings[kind]
        common = np.intersect1d(a["E"], b["E"])
        assert len(common) >= 20
        ia = [list(a["E"]).index(e) for e in common]
        ib = [list(b["E"]).index(e) for e in common]
        diff = (np.asarray(b["T"])[ib] - np.asarray(a["T"])[ia]) * 86400.0
        assert abs(diff.mean()) < 3.0
        assert diff.std() < 10.0


def test_aa_and_sector_17_borrowed_fit_clips_outliers(aa_and_reference):
    """Sector 17 (scattered-light ramps) fits to a few mmag after clipping."""
    _, _, out = _borrow(
        aa_and_reference, "TESS_SPOC_sector_17", ("amplitudes", "zero_phase")
    )
    assert out.solution.stats["n_clipped"] > 0
    assert out.solution.stats["rms_mag"] < 0.008
