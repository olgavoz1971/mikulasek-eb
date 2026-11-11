"""TESS adapter of ``mikulasek_eb``.

Everything that is specific to TESS lives here: the VOTable reader (through
``volightcurve``), the TESS exposure rule, the settings calibrated on TESS light curves, the
selection functions that use them, the folded-outlier and SAP-background cleaning, and the
checker of synthetic TESS-like files. The core
``mikulasek_eb`` does not import this package and does not need ``volightcurve``.

Install the extra with ``pip install mikulasek-eb[tess]``.
"""

from mikulasek_eb.selection import (
    DELTA_BIC_DECISIVE,
    Candidate,
    CandidateRow,
    FeatureSettings,
    Features,
    RuleSettings,
    SelectionResult,
    TieredResult,
    TieredStep,
    make_candidate,
)

from .auxiliary import (
    flag_folded_outliers,
    flag_sap_background,
    read_sap_background,
    wall_band,
    window_scan,
)
from .exposure import tess_exposure_days
from .io import load_tess_votable
from .rules import FEATURES, RULES
from .selection import (
    extract_features,
    propose_candidates,
    review_features,
    select_model,
    tiered_search,
)

__all__ = [
    "Candidate",
    "CandidateRow",
    "DELTA_BIC_DECISIVE",
    "FEATURES",
    "FeatureSettings",
    "Features",
    "RULES",
    "RuleSettings",
    "SelectionResult",
    "TieredResult",
    "TieredStep",
    "extract_features",
    "flag_folded_outliers",
    "flag_sap_background",
    "load_tess_votable",
    "make_candidate",
    "propose_candidates",
    "read_sap_background",
    "review_features",
    "select_model",
    "tess_exposure_days",
    "tiered_search",
    "wall_band",
    "window_scan",
]
