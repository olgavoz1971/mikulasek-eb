"""Tests of candidate rules and the evidence stage."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import mikulasek_eb as eb
from mikulasek_eb.config import Config, DataSetup, ModelSpec
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.model import EBModel, phase

ROOT = Path(__file__).resolve().parent.parent


def _load(name):
    p = load_tess_votable(ROOT / "data" / f"{name}.vot")
    return p, absolute_jd_to_mjd(p["t"])


def _features(name):
    p, t = _load(name)
    return p, t, eb.extract_features(t, p["mag"], p["sig"], p["period"], p["epoch_jd"])


def test_contact_stars_get_contact_harmonic_sets_and_detached_get_measured_ones():
    """EK Com: (1,2,3) and (2,3); AR Aur: the significant harmonics only."""
    _, _, f = _features("EK_Com_sector_49")
    sets = {c.config.model.harmonics for c in eb.propose_candidates(f)}
    assert sets == {(1, 2, 3), (2, 3)}
    _, _, f = _features("AR_Aur_sector_19")
    cands = eb.propose_candidates(f)
    assert len({c.config.model.harmonics for c in cands}) == 1
    assert all("separated eclipses" in c.reasons[0] for c in cands)


def test_candidates_are_valid_configs_with_reasons_and_baseline_choices():
    """Each candidate validates, has reasons, and the baseline orders are offered."""
    _, _, f = _features("lc_discovery_TIC_455790537_TESS_SPOC_sector_16")
    cands = eb.propose_candidates(f)
    assert {c.config.data.baseline_order for c in cands} == {0, 1}
    assert all(c.reasons for c in cands)
    assert all(
        c.config.data.exposure_days is not None for c in cands
    )  # 0.16 of the FWHM
    assert len({c.name for c in cands}) == len(cands)


def test_review_flags_odd_even_difference_as_possible_double_period():
    """Two unequal eclipse depths in alternate cycles ask for a human."""
    rng = np.random.default_rng(8)
    P = 0.5
    t = np.arange(58000.0, 58060.0, 10.0 / 1440.0)
    t = t + rng.uniform(-3.0, 3.0, len(t)) / 1440.0
    model = EBModel(1, (2,), False)
    base = dict(A1=0.3, A2=0.28, D1=0.04, G1=1.5, C1=0.0, H2=0.05)
    y = model(phase(t, 58000.1, P), base)
    cycle = np.floor((t - 58000.1) / P + 0.25).astype(int)
    y = y + np.where(cycle % 2 == 1, 0.04, 0.0) * np.exp(
        -0.5 * (((t - 58000.1) / P + 0.5) % 1 - 0.5) ** 2 / 0.03**2
    )
    y = y + rng.normal(0, 0.002, len(t))
    f = eb.extract_features(
        t, y, np.full(len(t), 0.002), P, 2458000.6, eb.FeatureSettings(n_bins=120)
    )
    ask, _ = eb.review_features(f)
    assert any("2P" in a for a in ask)


def test_xy_boo_selection_finds_the_papers_harmonic_set():
    """XY Boo: harmonics (2, 3), as in Eq. 18 of the paper, win on evidence."""
    p, t = _load("lc_discovery_TIC_347405567_TESS_SPOC_sector_50")
    r = eb.select_model(t, p["mag"], p["sig"], p["period"], p["epoch_jd"])
    assert r.best.config.model.harmonics == (2, 3)
    assert r.best.config.model.ocon is True
    assert r.n_eff_measured and r.n_eff < len(t)
    assert next(x for x in r.table if x.name == r.best.name).delta_bic == 0.0
    json.dumps(r.to_dict())
    assert r.output.solution.config == r.best.config


def test_failed_candidates_are_recorded_and_block_a_confident_answer():
    """An unsupported candidate (profile 5) shows in the table and in the reasons."""
    p, t = _load("EK_Com_sector_49")
    f = eb.extract_features(t, p["mag"], p["sig"], p["period"], p["epoch_jd"])
    good = eb.propose_candidates(f)[:1]
    bad = eb.Candidate(
        "profile5",
        Config(ModelSpec(5, (2,), False), DataSetup(0, None)),
        ("test",),
    )
    r = eb.select_model(
        t,
        p["mag"],
        p["sig"],
        p["period"],
        p["epoch_jd"],
        features=f,
        candidates=good + [bad],
    )
    row = next(x for x in r.table if x.name == "profile5")
    assert row.status.startswith("ConfigError")
    assert r.ask_human and any("profile5" in x for x in r.reasons)


def test_v477_lyr_with_a_shallow_secondary_now_fits():
    """The staged pipeline starts from the local-continuum depth, so V477 Lyr is fitted."""
    p, t = _load("v477_Lyr_sector_26")
    g = eb.initial_guess(t, p["mag"], p["period"])
    assert g["A2"] == pytest.approx(0.07, abs=0.03)


def test_two_component_twins_are_offered_only_when_asked():
    """``two_component_candidates`` adds a forced ``_2c`` twin of every candidate."""
    from mikulasek_eb.candidates import RuleSettings

    _, _, f = _features("EK_Com_sector_49")
    plain = eb.propose_candidates(f)
    assert not any(c.name.endswith("_2c") for c in plain)
    both = eb.propose_candidates(f, RuleSettings(two_component_candidates=True))
    assert len(both) == 2 * len(plain)
    twins = [c for c in both if c.name.endswith("_2c")]
    assert all(c.config.policy.force_two_component for c in twins)
    assert not any(c.config.policy.force_two_component for c in plain)
    assert len({c.name for c in both}) == len(both)
