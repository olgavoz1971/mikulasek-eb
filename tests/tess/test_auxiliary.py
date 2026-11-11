"""Folded-outlier and SAP-background cleaning, and the re-exported window scan."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import mikulasek_eb.tuning as tuning
import mikulasek_eb_tess as tess

ROOT = Path(__file__).resolve().parents[2]


def test_folded_rule_flags_a_planted_outlier_and_refuses_a_thin_bin():
    """One point 0.1 mag off the curve is flagged; a bin with too few points raises."""
    rng = np.random.default_rng(0)
    t = np.linspace(0.0, 10.0, 4000)
    mag = 0.01 * np.sin(2 * np.pi * t) + rng.normal(0.0, 0.001, len(t))
    mag[100] += 0.1
    flag, residual, phase = tess.flag_folded_outliers(
        t, mag, 1.0, n_bins=20, n_sigma=5.0, min_points=5
    )
    assert flag[100]
    assert flag.sum() < 10
    assert residual.shape == t.shape
    assert phase.shape == t.shape
    with pytest.raises(ValueError, match="phase bin"):
        tess.flag_folded_outliers(np.linspace(0, 1, 30), np.zeros(30), 1.0, n_bins=50)


def test_packaged_light_curve_can_be_folded_at_400_bins():
    """V477 Lyr has enough points in every bin for the default rule."""
    payload = tess.load_tess_votable(ROOT / "data" / "v477_Lyr_sector_26.vot")
    flag, _, _ = tess.flag_folded_outliers(
        payload["t"], payload["mag"], payload["period"]
    )
    assert flag.shape == payload["t"].shape
    assert flag.mean() < 0.05


def test_background_cut_is_a_ratio_of_the_median_and_the_times_must_match():
    """A point at 20 times the median is flagged; a missing time raises."""
    background = np.array([1.0, 1.0, 1.0, 1.0, 20.0])
    times = np.arange(5, dtype=float)
    flag, matched, centre = tess.flag_sap_background(
        times, times, background, ratio=10.0
    )
    assert flag.tolist() == [False, False, False, False, True]
    assert matched.tolist() == background.tolist()
    assert centre == 1.0
    with pytest.raises(ValueError, match="not positive"):
        tess.flag_sap_background(
            np.array([0.0]), np.array([0.0]), np.array([0.0]), ratio=10
        )
    with pytest.raises(ValueError, match="do not match"):
        tess.flag_sap_background(np.array([9.0]), times, background, ratio=10)


def test_packaged_files_have_no_sap_background_column():
    """The reduced test files do not carry sap_bkg, and the reader says so."""
    with pytest.raises(ValueError, match="sap_bkg"):
        tess.read_sap_background(ROOT / "data" / "EK_Com_sector_49.vot")


def test_named_flux_column_matches_the_single_column_load():
    """Naming the only flux column gives the same light curve as discovery."""
    path = ROOT / "data" / "EK_Com_sector_49.vot"
    plain = tess.load_tess_votable(path)
    named = tess.load_tess_votable(path, flux="flux", flux_err="flux_err")
    assert np.array_equal(plain["mag"], named["mag"])
    assert np.array_equal(plain["t"], named["t"])
    with pytest.raises(ValueError, match="not a flux column"):
        tess.load_tess_votable(path, flux="sap_flux", flux_err="sap_flux_err")
    with pytest.raises(ValueError, match="together"):
        tess.load_tess_votable(path, flux="flux")


def test_window_scan_is_the_core_function():
    """The adapter re-exports the scan; it does not keep a second copy."""
    assert tess.window_scan is tuning.window_scan
    assert tess.wall_band is tuning.wall_band
