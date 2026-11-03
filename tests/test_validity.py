"""Tests of the fit verdict: a numerical needle must be flagged and must not be timed."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import mikulasek_eb as eb
from mikulasek_eb import validity
from mikulasek_eb.errors import UnreliableSolution
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.model import EBModel, phase

GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
RUNAWAY_1800 = (
    "results/detached/detached_none_nooconnell_SP_L000_SP_L000_1800_c1800.vot",
    "results/detached/detached_none_nooconnell_SP_L180_SP_L180_1800_c1800.vot",
)
SPOT_FREE_120 = (
    "results/detached/detached_none_nooconnell_SP_none_SP_none_120_c120.vot"
)

NORMAL = dict(A1=0.5, A2=0.4, D1=0.02, G1=1.0, C1=0.0, D2=0.02, G2=1.0, C2=0.0)


def test_normal_template_has_no_limit_reason():
    """A typical detached eclipse pair passes the limits."""
    model = EBModel(harmonics=())
    assert validity._limit_reasons(model, NORMAL, 0.5, None) == []


@pytest.mark.parametrize(
    "change",
    [
        dict(G1=2.4e-6, A1=38000.0),
        dict(G1=0.05),
        dict(G1=100.0),
        dict(D1=0.25, G1=2.2),
        dict(A1=5.0),
        dict(A1=-0.2),
    ],
)
def test_runaway_and_unphysical_templates_are_bound_active(change):
    """A needle, a too-wide eclipse and an impossible depth each give a reason."""
    model = EBModel(harmonics=())
    p = {**NORMAL, **change}
    assert validity._limit_reasons(model, p, 0.5, None)


def test_two_component_negative_amplitude_is_allowed():
    """The second component may be negative; only its size is limited."""
    model = EBModel(harmonics=())
    p = {**NORMAL, "A1b": -0.1, "D1b": 0.03, "G1b": 1.5}
    assert validity._limit_reasons(model, p, 0.5, None) == []


def test_fwhm_matches_the_closed_form():
    """The numerical FWHM of profile 1 agrees with 2 D arccosh(1 - ln(1 - 2^(-1/G)))."""
    model = EBModel(harmonics=())
    comp = validity.components(model, NORMAL)[0]
    exact = 2 * NORMAL["D1"] * np.arccosh(1 - np.log(1 - 2 ** (-1 / NORMAL["G1"])))
    assert validity.component_fwhm(model, comp) == pytest.approx(exact, rel=2e-3)


def test_spike_narrower_than_the_sub_step_is_not_converged():
    """A spike narrower than the exposure sub-step gives a comb and is refused."""
    period = 0.87
    model = EBModel(harmonics=(), exposure_days=0.0167, nsub=21, period_days=period)
    spike = {**NORMAL, "A1": 40.0, "D1": 0.055, "G1": 1e-5}
    assert validity._exposure_reason(model, spike, 0.012) is not None
    assert validity._exposure_reason(model, NORMAL, 0.012) is None


def _synthetic(seed=3):
    """Two-eclipse light curve from the model, 20 cycles at 10 minute cadence."""
    rng = np.random.default_rng(seed)
    period = 0.87
    t = 58000.0 + np.arange(0.0, 20 * period, 600.0 / 86400.0)
    model = EBModel(harmonics=())
    mag = 10.0 + model(phase(t, 58000.0, period), NORMAL)
    sig = np.full(len(t), 0.006)
    return t, mag + rng.normal(0.0, 0.006, len(t)), sig, period


def test_clean_synthetic_fit_is_ok_and_can_be_timed():
    """A well-sampled clean curve is ok, and ``time`` runs."""
    t, mag, sig, period = _synthetic()
    cfg = eb.Config(
        eb.ModelSpec(profile=1, harmonics=(), ocon=False),
        eb.DataSetup(baseline_order=1, gap_split_days=None),
        eb.FitPolicy(allow_two_component=False),
        eb.TimingRecipe(window=0.12),
    )
    out = eb.fit(t, mag, sig, period, 58000.0 + 2400000.5, cfg)
    assert out.solution.stats["status"] == "ok"
    result = eb.time(out.t_mjd, out.mag, out.sig, out.solution)
    assert len(result.timings["primary"]["T"]) > 0


def test_time_refuses_a_flagged_solution():
    """``time`` raises for any status but ok, and for a missing status."""
    t, mag, sig, period = _synthetic()
    cfg = eb.Config(
        eb.ModelSpec(profile=1, harmonics=(), ocon=False),
        eb.DataSetup(baseline_order=1, gap_split_days=None),
        eb.FitPolicy(allow_two_component=False),
        eb.TimingRecipe(window=0.12),
    )
    out = eb.fit(t, mag, sig, period, 58000.0 + 2400000.5, cfg)
    solution = out.solution
    for status in ("bound-active", "unresolved"):
        solution.stats["status"] = status
        solution.stats["status_reasons"] = ["test"]
        with pytest.raises(UnreliableSolution):
            eb.time(out.t_mjd, out.mag, out.sig, solution)
    del solution.stats["status"]
    with pytest.raises(UnreliableSolution):
        eb.time(out.t_mjd, out.mag, out.sig, solution)


def _fit_generated(relative: str):
    """Fit one generated curve with the degraded-spot structure (profile 1, H2)."""
    import dataclasses

    from mikulasek_eb.candidates import propose_candidates
    from mikulasek_eb.features import extract_features

    path = GENERATED / relative
    if not path.is_file():
        pytest.skip(f"{path} not available")
    payload = load_tess_votable(path)
    t = absolute_jd_to_mjd(payload["t"])
    features = extract_features(
        t, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    own = propose_candidates(features)[0].config.data
    cfg = eb.Config(
        eb.ModelSpec(profile=1, harmonics=(2,), ocon=True),
        dataclasses.replace(own, baseline_order=1),
        eb.FitPolicy(allow_two_component=True),
        eb.TimingRecipe(window=0.12),
    )
    return eb.fit(
        t, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"], cfg,
        data_id=path.stem,
    )


@pytest.mark.slow
@pytest.mark.parametrize("relative", RUNAWAY_1800)
def test_runaway_1800_curves_are_flagged_and_not_timed(relative):
    """The two 1800 s curves whose primary ran away are not reported as ok."""
    out = _fit_generated(relative)
    assert out.solution.stats["status"] != "ok", out.solution.stats
    with pytest.raises(UnreliableSolution):
        eb.time(out.t_mjd, out.mag, out.sig, out.solution)


@pytest.mark.slow
def test_spot_free_120_stays_ok():
    """The spot-free 120 s curve is ok and is still timed."""
    out = _fit_generated(SPOT_FREE_120)
    assert out.solution.stats["status"] == "ok", out.solution.stats
    result = eb.time(out.t_mjd, out.mag, out.sig, out.solution)
    assert len(result.timings["primary"]["T"]) > 0
