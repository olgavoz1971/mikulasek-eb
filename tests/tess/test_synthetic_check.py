"""Tests of the synthetic-file checker with in-memory files (no disk)."""

from __future__ import annotations

import dataclasses
from pathlib import Path

import numpy as np
import pytest

from mikulasek_eb_tess.synthetic_check import (
    FAIL,
    PASS,
    WARN,
    SyntheticFile,
    check_synthetic,
    load_synthetic,
)

P = 0.34
CADENCE_S = 600.0
EXPOSURE_S = 480.0
GAP = (13.0, 14.0)


def _pure_flux(t):
    """Smooth pure flux with two minima per period (analytic, so tests are exact)."""
    return 1.3 - 0.3 * np.cos(4.0 * np.pi * np.asarray(t) / P)


def _smeared_flux(t, shift=0.0):
    """Average of the pure flux over the exposure window (optionally shifted)."""
    u = shift + (np.arange(201) / 200 - 0.5) * EXPOSURE_S / 86400.0
    return np.mean([_pure_flux(t + x) for x in u], axis=0)


def make(window_shift=0.0):
    """Build a consistent synthetic file and its pure curve in memory.

    Args:
        window_shift (float): Shift of the exposure window used for the truth (days).

    Returns:
        tuple[SyntheticFile, SyntheticFile]: The spoiled file and the pure curve.
    """
    rng = np.random.default_rng(3)
    t = 0.0035 + np.arange(int(27 * 86400 / CADENCE_S)) * CADENCE_S / 86400.0
    t = t[(t < GAP[0]) | (t > GAP[1])]
    tp = _smeared_flux(t, window_shift)
    tr = 0.02 * np.sin(2 * np.pi * t / 13.7)
    sigma = 0.01 * tp
    flux = tp + tr + rng.normal(0.0, sigma)
    E = np.arange(0, 80)
    tc = np.ravel(np.column_stack([E * P, E * P + 0.5 * P]))
    kind = ["primary", "secondary"] * len(E)
    tm = tc + np.where(np.arange(len(tc)) % 2 == 0, 0.05, 0.1) / 86400.0
    depth = float(tp.max() - tp.min())
    meta = dict(
        period=P,
        epoch=0.0,
        product_kind="tess_observe",
        recipe_id="S05",
        physics_id="G5",
        source_file="pure.vot",
        cadence_s=CADENCE_S,
        exposure_s=EXPOSURE_S,
        smear_frac=0.8,
        smear_nsub=201,
        duration_days=27.0,
        gap_start_time=GAP[0],
        gap_end_time=GAP[1],
        noise_type="white",
        rel_sigma=0.01,
        sigma_factor=1.0,
        seed=3,
        P_true_days=P,
        P_meta_days=P,
        period_eps=0.0,
        trend_types="ramps",
        trend_model="additive",
        trend_amp_frac=float(np.ptp(tr) / depth),
        eclipse_depth_mag=float(2.5 * np.log10(tp.max() / tp.min())),
        time_scale="TDB",
        time_refposition="BARYCENTER",
        timeorigin_jd=2458325.0,
        eb_class="contact",
        fillout_factor=0.4,
        has_oconnell="false",
        has_spots="false",
        eccentric="false",
        total_eclipse="false",
        n_eclipses=len(tc),
        eclipse_E=",".join(str(int(e)) for e in np.repeat(E, 2)),
        eclipse_kind=",".join(kind),
        eclipse_T_conj_days=",".join(f"{x:.12g}" for x in tc),
        eclipse_T_min_pure_days=",".join(f"{x:.12g}" for x in tm),
        eclipse_O_C_injected_days=",".join("0" for _ in tc),
        conj_minus_min_pure_s=-0.1,
        eclipse_in_gap="1",
    )
    data = SyntheticFile(
        t,
        flux,
        sigma,
        tp,
        tr,
        meta,
        dict(timescale="TDB", refposition="BARYCENTER", timeorigin=2458325.0),
    )
    grid = np.arange(0.0, P, 40.0 / 86400.0)
    pure = SyntheticFile(grid, _pure_flux(grid), np.array([]), None, None, {}, {})
    return data, pure


def status(findings, name):
    """Status of the finding(s) with this name (the worst one if several)."""
    rank = {PASS: 0, "info": 0, "skip": 0, WARN: 1, FAIL: 2}
    hits = [f for f in findings if f.name == name]
    assert hits, f"no finding named {name!r}: {[f.name for f in findings]}"
    return max((f.status for f in hits), key=rank.get)


def with_meta(data, **changes):
    """Copy of a file with some parameters changed (None removes a key)."""
    meta = dict(data.meta)
    for k, v in changes.items():
        if v is None:
            meta.pop(k, None)
        else:
            meta[k] = v
    return dataclasses.replace(data, meta=meta)


def test_a_consistent_file_has_no_failures_and_no_warnings():
    """The reference in-memory file passes every check."""
    data, pure = make()
    findings = check_synthetic(data, pure)
    bad = [f for f in findings if f.status in (FAIL, WARN)]
    assert not bad, bad
    assert status(findings, "smeared truth") == PASS


@pytest.mark.parametrize(
    "changes,name,expected",
    [
        (dict(cadence_s=300.0), "cadence", FAIL),
        (dict(exposure_s=600.0), "exposure", FAIL),
        (dict(gap_start_time=10.0), "gap", FAIL),
        (dict(seed=None), "parameters", FAIL),
        (dict(eb_class="binary"), "labels", FAIL),
        (dict(fillout_factor=None), "labels", FAIL),
        (dict(has_spots="maybe"), "labels", FAIL),
        (dict(time_scale="UTC"), "time system", FAIL),
        (dict(trend_amp_frac=0.15), "trend amplitude", WARN),
        (dict(sigma_factor=2.0), "noise", WARN),
        (dict(eclipse_depth_mag=0.5), "eclipse depth", WARN),
        (dict(smear_nsub=None), "smear_nsub", WARN),
        (dict(eclipse_in_gap=None), "eclipse flags", WARN),
        (dict(conj_minus_min_pure_s=-0.01), "conj_minus_min_pure_s", WARN),
        (dict(n_eclipses=10), "eclipse truth", FAIL),
        (dict(period=0.35), "period parameter", FAIL),
    ],
)
def test_each_inconsistency_is_reported(changes, name, expected):
    """A single wrong or missing item produces the expected finding."""
    data, pure = make()
    findings = check_synthetic(with_meta(data, **changes), pure)
    assert status(findings, name) == expected


def test_missing_truth_columns_fail():
    """Without the truth columns the check fails and says so."""
    data, pure = make()
    findings = check_synthetic(
        dataclasses.replace(data, truth_flux=None, truth_trend=None), pure
    )
    assert status(findings, "truth columns") == FAIL


@pytest.mark.parametrize(
    "shift", [-0.5 * EXPOSURE_S / 86400.0, 0.5 * EXPOSURE_S / 86400.0]
)
def test_a_wrongly_placed_exposure_window_is_detected(shift):
    """Truth built with a window that starts or ends at the stamp fails."""
    data, pure = make(window_shift=shift)
    assert status(check_synthetic(data, pure), "smeared truth") == FAIL


def test_a_spotted_system_may_have_a_minimum_offset():
    """For an eccentric or spotted system the offset is information, not a warning."""
    data, pure = make()
    tm = np.array([float(x) for x in data.meta["eclipse_T_min_pure_days"].split(",")])
    tm = tm + 5.0 / 86400.0
    spotted = with_meta(
        data,
        has_spots="true",
        eclipse_T_min_pure_days=",".join(f"{x:.12g}" for x in tm),
    )
    assert (
        status(check_synthetic(spotted, pure), "T_min_pure - T_conj (primary)")
        == "info"
    )
    clean = with_meta(
        data, eclipse_T_min_pure_days=spotted.meta["eclipse_T_min_pure_days"]
    )
    assert status(check_synthetic(clean, pure), "T_min_pure - T_conj (primary)") == WARN


def test_without_the_pure_curve_the_smear_check_is_skipped():
    """No pure curve: the check is skipped, not passed."""
    data, _ = make()
    assert status(check_synthetic(data, None), "smeared truth") == "skip"


def test_unreadable_file_is_reported(tmp_path):
    """A file that volightcurve cannot read gives a failing finding, not a crash."""
    bad = tmp_path / "bad.vot"
    bad.write_text("not a votable")
    findings = check_synthetic(load_synthetic(bad))
    assert findings and findings[0].status == FAIL


REAL = Path("/home/voz/projects/UPJS/EB_synthetic_lc/results/n01/contact")
REAL_FILE = REAL / "contact_overcontact_nooconnell_G5_S05_c600.vot"


@pytest.mark.skipif(not REAL_FILE.exists(), reason="generator output not available")
def test_the_real_generator_file_has_no_failures():
    """The first real file passes every hard check (warnings are known and listed)."""
    data = load_synthetic(REAL_FILE)
    pure = load_synthetic(REAL / data.meta["source_file"], with_errors=False)
    findings = check_synthetic(data, pure)
    assert not [f for f in findings if f.status == FAIL], findings


def test_a_wrong_period_in_the_file_does_not_break_the_truth_checks():
    """With period_eps the truth uses P_true while the table carries P_meta."""
    data, pure = make()
    eps = 1e-4
    wrong = with_meta(
        data, period=P * (1 + eps), P_meta_days=P * (1 + eps), period_eps=eps
    )
    findings = check_synthetic(wrong, pure)
    assert not [f for f in findings if f.status in (FAIL, WARN)], findings
    assert status(findings, "conjunction times") == PASS
    assert status(findings, "smeared truth") == PASS
