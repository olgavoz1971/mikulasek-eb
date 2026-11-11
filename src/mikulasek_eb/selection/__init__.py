"""Selection machinery: measure features, propose candidates, choose by evidence.

This layer holds no survey numbers. The thresholds of the measurement
(``FeatureSettings``), the thresholds of the rules (``RuleSettings``) and the exposure rule
are arguments that a survey adapter supplies, for example ``mikulasek_eb_tess``. The core
(``mikulasek_eb``) never imports this layer.

``select_model`` fits a given list of candidates. ``tiered_search`` is a named strategy built on it:
baseline, then profile, then harmonics.
"""

from .candidates import Candidate, RuleSettings, propose_candidates, review_features
from .features import (
    FeatureSettings,
    Features,
    extract_features,
    gamma_from_width_ratios,
)
from .search import DELTA_BIC_DECISIVE, CandidateRow, SelectionResult, select_model
from .tiered import TieredResult, TieredStep, make_candidate, tiered_search

__all__ = [
    "Candidate",
    "CandidateRow",
    "DELTA_BIC_DECISIVE",
    "FeatureSettings",
    "Features",
    "RuleSettings",
    "SelectionResult",
    "TieredResult",
    "TieredStep",
    "extract_features",
    "gamma_from_width_ratios",
    "make_candidate",
    "propose_candidates",
    "review_features",
    "select_model",
    "tiered_search",
]
