"""Tests of the host-facing API: ``fit``, ``time``, ``ModelSolution`` and typed errors."""

from __future__ import annotations

import json

import numpy as np
import pytest

import mikulasek_eb as eb
from mikulasek_eb_tess.io import load_tess_votable
from mikulasek_eb.errors import (
    ConfigError,
    FitCancelled,
    InvalidLightCurve,
    MikulasekError,
)
from mikulasek_eb.helpers import absolute_jd_to_mjd


@pytest.fixture(scope="module")
def ek_com_api():
    """Run ``eb.fit`` on EK Com with the per-star config; return payload and output."""
    import fit_ek_com_s49 as star

    payload = load_tess_votable(star.VOTABLE_PATH)
    stages = []
    out = eb.fit(
        absolute_jd_to_mjd(payload["t"]),
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        eb.Config.from_fit_config(star.FIT_CONFIG),
        data_id="EK_Com_S49",
        progress=lambda stage, frac: stages.append((stage, frac)),
    )
    return payload, out, stages


def test_exceptions_keep_their_builtin_bases():
    """Typed errors stay catchable as ValueError or RuntimeError."""
    assert issubclass(ConfigError, ValueError)
    assert issubclass(InvalidLightCurve, ValueError)
    assert issubclass(FitCancelled, RuntimeError)
    assert issubclass(FitCancelled, MikulasekError)


@pytest.mark.parametrize(
    "mutate",
    [
        lambda t, m, s: (t[:-1], m, s),
        lambda t, m, s: (t, np.where(np.arange(len(m)) == 3, np.nan, m), s),
        lambda t, m, s: (t, m, np.where(np.arange(len(s)) == 3, -1.0, s)),
    ],
)
def test_validate_lightcurve_rejects_bad_input(mutate):
    """Length mismatch, non-finite values and non-positive errors raise InvalidLightCurve."""
    t = np.linspace(0.0, 5.0, 200)
    m = np.zeros(200)
    s = np.full(200, 0.01)
    t, m, s = mutate(t, m, s)
    with pytest.raises(InvalidLightCurve):
        eb.validate_lightcurve(t, m, s, 1.0, 2459000.0)


def test_solution_json_round_trip_and_hash(ek_com_api):
    """The solution survives JSON, keeps its hash, and rebuilds the same model."""
    _, out, _ = ek_com_api
    sol = out.solution
    again = eb.ModelSolution.from_json(sol.to_json())
    assert again.solution_hash() == sol.solution_hash()
    json.loads(sol.to_json())
    m1 = sol.build_model()
    m2 = again.build_model()
    assert m1.names() == m2.names()
    assert sol.provenance["data_id"] == "EK_Com_S49"
    assert sol.stats["n_points"] > 0


def test_time_reproduces_reference_counts(ek_com_api):
    """``eb.time`` with the solution default recipe returns timings and per-type rejects."""
    payload, out, _ = ek_com_api
    res = eb.time(
        out.t_mjd,
        out.mag,
        out.sig,
        out.solution,
        eb.TimingRecipe(require_contact_walls=False, clip_outliers=False),
    )
    assert set(res.timings) == {"primary", "secondary"}
    assert len(res.timings["primary"]["T"]) == 70
    assert len(res.timings["secondary"]["T"]) == 68
    assert res.solution_hash == out.solution.solution_hash()
    json.dumps(res.to_dict())


def test_time_defaults_to_solution_recipe(ek_com_api):
    """With no recipe, the one stored in the solution config is used."""
    _, out, _ = ek_com_api
    res = eb.time(out.t_mjd, out.mag, out.sig, out.solution, which=("primary",))
    assert res.recipe == out.solution.config.timing
    assert set(res.timings) == {"primary"}


def test_progress_hook_reports_increasing_fractions(ek_com_api):
    """The hook is called with rising fractions that end at 1."""
    _, _, stages = ek_com_api
    fracs = [f for _, f in stages]
    assert fracs == sorted(fracs)
    assert fracs[-1] == 1.0


def test_cancel_from_progress_hook_propagates():
    """An exception raised by the hook (FitCancelled) stops the fit unchanged."""
    import fit_ek_com_s49 as star

    payload = load_tess_votable(star.VOTABLE_PATH)

    def stop(stage, frac):
        raise FitCancelled("user cancelled")

    with pytest.raises(FitCancelled):
        eb.fit(
            absolute_jd_to_mjd(payload["t"]),
            payload["mag"],
            payload["sig"],
            payload["period"],
            payload["epoch_jd"],
            eb.Config.from_fit_config(star.FIT_CONFIG),
            progress=stop,
        )
