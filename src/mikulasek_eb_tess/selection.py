"""The selection machinery with the TESS calibration filled in.

Each function has the signature of its survey-neutral counterpart in
``mikulasek_eb.selection``. A missing settings argument is replaced by the calibrated TESS
value, and a missing exposure rule by the TESS rule. The values are in ``rules.py``.
"""

from __future__ import annotations

from mikulasek_eb import selection as neutral
from mikulasek_eb.selection import Features, FeatureSettings, RuleSettings

from .exposure import tess_exposure_days
from .rules import FEATURES, RULES


def extract_features(
    t_mjd,
    mag,
    sig,
    period_days: float,
    epoch_jd: float,
    settings: FeatureSettings | None = None,
    *,
    exposure_days: float | None = None,
) -> Features:
    """Measure the features of one TESS light curve.

    Args:
        t_mjd (ndarray): Times (MJD), strictly increasing.
        mag (ndarray): Magnitudes.
        sig (ndarray): Catalogue magnitude errors.
        period_days (float): Catalogue period (days).
        epoch_jd (float): Catalogue epoch (absolute JD).
        settings (FeatureSettings | None): None uses the calibrated TESS values.
        exposure_days (float | None): Exposure (days). None uses the TESS exposure rule.

    Returns:
        Features: The measurements.

    Raises:
        InvalidLightCurve: If the input is invalid or the fold is too sparse.
    """
    return neutral.extract_features(
        t_mjd,
        mag,
        sig,
        period_days,
        epoch_jd,
        FEATURES if settings is None else settings,
        exposure_days=exposure_days,
        exposure_rule=tess_exposure_days,
    )


def review_features(features: Features, rules: RuleSettings | None = None):
    """Review the features with the calibrated TESS rules.

    Args:
        features (Features): Measured features.
        rules (RuleSettings | None): None uses the calibrated TESS values.

    Returns:
        The result of ``mikulasek_eb.selection.review_features``.
    """
    return neutral.review_features(features, RULES if rules is None else rules)


def propose_candidates(features: Features, rules: RuleSettings | None = None) -> list:
    """Propose candidate configurations with the calibrated TESS rules.

    Args:
        features (Features): Measured features.
        rules (RuleSettings | None): None uses the calibrated TESS values.

    Returns:
        list: Candidates, as ``mikulasek_eb.selection.propose_candidates``.
    """
    return neutral.propose_candidates(features, RULES if rules is None else rules)


def select_model(
    t_mjd,
    mag,
    sig,
    period_days: float,
    epoch_jd: float,
    *,
    features: Features | None = None,
    candidates: list | None = None,
    rules: RuleSettings | None = None,
    feature_settings: FeatureSettings | None = None,
    data_id: str | None = None,
    progress=None,
    retry_not_converged: bool = False,
) -> neutral.SelectionResult:
    """Choose a model for one TESS light curve by fitted evidence.

    Args:
        t_mjd (ndarray): Times (MJD), strictly increasing.
        mag (ndarray): Magnitudes.
        sig (ndarray): Catalogue errors.
        period_days (float): Catalogue period (days).
        epoch_jd (float): Catalogue epoch (absolute JD).
        features (Features | None): Precomputed features; None measures them.
        candidates (list | None): Candidates to fit; None proposes them.
        rules (RuleSettings | None): None uses the calibrated TESS rules.
        feature_settings (FeatureSettings | None): None uses the calibrated TESS values.
        data_id (str | None): Identifier for provenance.
        progress (callable | None): ``progress(stage, fraction)``.
        retry_not_converged (bool): See ``mikulasek_eb.selection.select_model``.

    Returns:
        SelectionResult: Choice, comparison table and flags.

    Raises:
        InvalidLightCurve: If the input is invalid.
    """
    return neutral.select_model(
        t_mjd,
        mag,
        sig,
        period_days,
        epoch_jd,
        features=features,
        candidates=candidates,
        rules=RULES if rules is None else rules,
        feature_settings=FEATURES if feature_settings is None else feature_settings,
        exposure_rule=tess_exposure_days,
        data_id=data_id,
        progress=progress,
        retry_not_converged=retry_not_converged,
    )


def tiered_search(
    t_mjd,
    mag,
    sig,
    period_days: float,
    epoch_jd: float,
    *,
    features: Features,
    profiles,
    start_profile: int,
    rules: RuleSettings | None = None,
    decisive_bic: float = neutral.DELTA_BIC_DECISIVE,
    data_id: str | None = None,
    on_step=None,
) -> neutral.TieredResult:
    """Run the tiered search (baseline, profile, harmonics) with the calibrated TESS rules.

    Args:
        t_mjd (ndarray): Times (MJD), strictly increasing.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        period_days (float): Catalogue period (days).
        epoch_jd (float): Catalogue epoch (absolute JD).
        features (Features): Measured features of this light curve.
        profiles (Sequence[int]): Table 2 profiles that the profile step may try.
        start_profile (int): Profile of the baseline step.
        rules (RuleSettings | None): None uses the calibrated TESS rules.
        decisive_bic (float): Decisive BIC difference of the stop rule.
        data_id (str | None): Identifier for provenance.
        on_step (callable | None): ``on_step(label, result, seconds)`` after each step.

    Returns:
        TieredResult: As ``mikulasek_eb.selection.tiered_search``.

    Raises:
        NoValidSolution: If a step finished without an acceptable candidate.
    """
    return neutral.tiered_search(
        t_mjd,
        mag,
        sig,
        period_days,
        epoch_jd,
        features=features,
        rules=RULES if rules is None else rules,
        profiles=profiles,
        start_profile=start_profile,
        decisive_bic=decisive_bic,
        data_id=data_id,
        on_step=on_step,
    )
