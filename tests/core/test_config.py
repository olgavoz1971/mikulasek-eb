"""Layered configuration: validation, JSON round trip, hash and the flat view."""

from __future__ import annotations

import dataclasses
import json

import pytest

from mikulasek_eb import (
    Config,
    ConfigError,
    DataSetup,
    FitConfig,
    FitPolicy,
    ModelSpec,
    TimingRecipe,
)


def _config(**overrides) -> Config:
    """Small valid configuration for the tests."""
    model = ModelSpec(profile=1, harmonics=(2, 3), ocon=True)
    data = DataSetup(baseline_order=0, gap_split_days=1.0)
    return dataclasses.replace(Config(model, data), **overrides)


def test_json_round_trip_restores_tuples_and_hash():
    """Serialise and parse: identical object, identical hash."""
    cfg = _config()
    again = Config.from_json(cfg.to_json())
    assert again == cfg
    assert isinstance(again.model.harmonics, tuple)
    assert again.config_hash() == cfg.config_hash()
    assert len(cfg.config_hash()) == 12


def test_hash_ignores_notes_but_not_settings():
    """Free-text notes do not change the identity; any setting does."""
    cfg = _config()
    assert dataclasses.replace(cfg, notes="hello").config_hash() == cfg.config_hash()
    changed = dataclasses.replace(cfg, timing=TimingRecipe(window=0.15))
    assert changed.config_hash() != cfg.config_hash()


def test_save_and_load(tmp_path):
    """A file written by ``save`` loads back unchanged."""
    cfg = _config(notes="EK Com sector 49")
    path = tmp_path / "cfg.json"
    cfg.save(path)
    assert Config.load(path) == cfg


def test_policy_and_timing_may_be_omitted_in_json():
    """Defaults fill omitted optional layers; required layers are enforced."""
    doc = dict(
        model=dict(profile=1, harmonics=[2], ocon=False),
        data=dict(baseline_order=1, gap_split_days=None),
    )
    cfg = Config.from_dict(doc)
    assert cfg.policy == FitPolicy() and cfg.timing == TimingRecipe()
    with pytest.raises(ConfigError, match="Missing required layer"):
        Config.from_dict(dict(model=doc["model"]))


@pytest.mark.parametrize(
    "mutate,match",
    [
        (lambda d: d["model"].update(profile=9), "profile"),
        (lambda d: d["model"].update(profile=2, kterm=True), "kterm"),
        (lambda d: d["model"].update(typo=1), "Unknown keys"),
        (lambda d: d.update(extra=1), "Unknown top-level"),
        (lambda d: d["data"].update(exposure_nsub=5), "exposure_days is None"),
        (lambda d: d["data"].update(gap_split_days=-1.0), "gap_split_days"),
        (lambda d: d["policy"].update(gamma_starts=[]), "gamma_starts"),
        (lambda d: d["policy"].update(two_comp_eclipses=[3]), "two_comp_eclipses"),
        (lambda d: d["timing"].update(window=0.4), "window"),
        (lambda d: d["timing"].update(template="reference"), "not implemented"),
        (lambda d: d.update(schema_version=99), "schema_version"),
    ],
)
def test_invalid_configurations_fail_fast(mutate, match):
    """Impossible, unknown or unsupported settings raise ConfigError naming the problem."""
    doc = _config().to_dict()
    mutate(doc)
    with pytest.raises(ConfigError, match=match):
        Config.from_dict(doc)


def test_invalid_json_text():
    """Garbage text raises ConfigError, not a JSON exception."""
    with pytest.raises(ConfigError, match="not valid JSON"):
        Config.from_json("{not json")


def test_flat_view_covers_exactly_the_three_layers():
    """FitConfig fields are the union of model, data and policy fields (no drift)."""
    layers = (ModelSpec, DataSetup, FitPolicy)
    union = [f.name for layer in layers for f in dataclasses.fields(layer)]
    assert len(union) == len(set(union)), "layer fields must have unique names"
    assert {f.name for f in dataclasses.fields(FitConfig)} == set(union)


def test_flat_and_layered_views_round_trip():
    """Config -> FitConfig -> Config keeps the model, data and policy layers."""
    cfg = _config()
    again = Config.from_fit_config(cfg.to_fit_config(), timing=cfg.timing)
    assert again == cfg


def test_defaults_match_the_values_used_before_the_refactor():
    """Guard the numbers that used to be module constants."""
    policy, timing = FitPolicy(), TimingRecipe()
    assert policy.min_points_per_fwhm == 8.0
    assert policy.primary_a2_over_a1_min_ratio == 1.2
    assert policy.two_comp_min_rms_gain_frac == 0.25
    assert policy.two_comp_min_agreeing_starts == 2
    assert policy.two_comp_max_rel_err == 0.5
    assert policy.two_comp_eclipses == (1,)
    assert policy.two_comp_timing_agreement_sigma_max == 1.0
    assert (timing.window, timing.baseline_order, timing.max_sigma_factor) == (
        0.25,
        1,
        5.0,
    )
    assert timing.local_windows == (0.15, 0.10, 0.06)
    json.loads(_config().to_json())  # valid JSON


def test_force_two_component_needs_allow_two_component():
    """Forcing the trial while forbidding it is a configuration error."""
    cfg = _config()
    bad = dataclasses.replace(
        cfg,
        policy=dataclasses.replace(
            cfg.policy, allow_two_component=False, force_two_component=True
        ),
    )
    with pytest.raises(ConfigError):
        bad.validate()


def test_timing_window_above_quarter_needs_overlap_flag():
    """A window up to 0.5 is valid only with ``allow_overlapping_windows``."""
    from mikulasek_eb.config import TimingRecipe

    base = Config(ModelSpec(profile=1, harmonics=(2,), ocon=False), DataSetup(0, None))
    wide = TimingRecipe(window=0.5)
    with pytest.raises(ConfigError):
        dataclasses.replace(base, timing=wide).validate()
    dataclasses.replace(
        base, timing=dataclasses.replace(wide, allow_overlapping_windows=True)
    ).validate()
    with pytest.raises(ConfigError):
        dataclasses.replace(base, timing=TimingRecipe(min_points=4, min_each_side=3)).validate()
