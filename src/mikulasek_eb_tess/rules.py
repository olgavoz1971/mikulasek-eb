"""Calibrated settings of the TESS adapter.

These numbers are our own heuristics, set on the test stars with TESS light curves
(cadences of 2, 10 and 30 minutes, TESS noise and systematics). They were the defaults
of ``FeatureSettings`` and ``RuleSettings`` before the survey-neutral layer was
separated. They are not expected to hold for another survey and must not be reused
there without a new calibration.

Change one value with ``dataclasses.replace(FEATURES, gap_days=0.3)``.
"""

from __future__ import annotations

import dataclasses

from mikulasek_eb.selection import ContactSettings, FeatureSettings, RuleSettings

CONTACTS = ContactSettings(
    scan_limit_primary=0.25,
    scan_limit_secondary=0.2,
    start_over_fwhm=0.6,
    slope_narrow=0.25,
    slope_wide=0.75,
    min_points_in_window=5,
    min_sigma=10.0,
    max_tail=0.1,
    max_shift=1.2,
    a_k_values=(1.0, 1.5, 2.0),
    a_edge_frac=0.08,
    a_max_out_over_in=0.3,
    a_max_spread=0.10,
    a_min_points=8,
    agree_tol=0.15,
)

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
    contacts=CONTACTS,
)

# The earlier measurement: eclipse windows at the 5 percent edge, no contact search.
FEATURES_EDGE_ONLY = dataclasses.replace(FEATURES, contacts=None)

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
