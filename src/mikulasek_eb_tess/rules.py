"""Calibrated settings of the TESS adapter.

These numbers are our own heuristics, set on the test stars with TESS light curves
(cadences of 2, 10 and 30 minutes, TESS noise and systematics). They were the defaults
of ``FeatureSettings`` and ``RuleSettings`` before the survey-neutral layer was
separated. They are not expected to hold for another survey and must not be reused
there without a new calibration.

Change one value with ``dataclasses.replace(FEATURES, gap_days=0.3)``.
"""

from __future__ import annotations

from mikulasek_eb.selection import FeatureSettings, RuleSettings

FEATURES = FeatureSettings(
    n_bins=None,
    min_points_per_bin=3,
    reference_halfwidth=0.05,
    core_halfmax=0.25,
    vertex_core_frac=0.8,
    bottom_core_frac=0.9,
    local_inner=0.03,
    local_outer=0.25,
    secondary_search=0.2,
    min_prominence_sigma=5.0,
    min_secondary_width_bins=3.0,
    edge_level=0.05,
    window_margin=1.0,
    flat_bottom_ratio=0.6,
    separated_coverage=0.3,
    equal_depth_ratio=0.8,
    harmonics_max=4,
    gap_days=0.5,
    step_window_days=1.0,
    min_step_points=10,
    trend_bin_days=1.0,
    outlier_nsig=5.0,
    gamma_match_tol=0.15,
    parity_window_fwhm=0.25,
)

RULES = RuleSettings(
    harmonic_sigma=3.0,
    ocon_sigma=3.0,
    exposure_over_fwhm_min=0.1,
    baseline_orders=(0, 1),
    parity_sigma=3.0,
    offset_sigma=3.0,
    offset_min=0.01,
    shallow_depth_ratio=0.15,
    min_points_per_fwhm=8.0,
    max_outlier_fraction=0.02,
    drift_sigma=3.0,
    two_component_candidates=False,
)
