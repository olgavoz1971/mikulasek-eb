"""The TESS exposure rule of the adapter."""

from __future__ import annotations

import numpy as np
import pytest

from mikulasek_eb_tess.exposure import tess_exposure_days


def test_tess_exposure_is_80_percent_of_the_cadence_except_at_20_s():
    """480 s at the 600 s cadence; the 20 s cadence has no cosmic-ray mitigation."""
    for cadence, factor in ((600.0, 0.8), (120.0, 0.8), (1800.0, 0.8), (20.0, 1.0)):
        t = np.arange(0.0, 2.0, cadence / 86400.0)
        days, source = tess_exposure_days(t)
        assert days * 86400.0 == pytest.approx(factor * cadence)
        assert source.startswith("TESS")
    with pytest.raises(ValueError, match="irregular"):
        rng = np.random.default_rng(2)
        tess_exposure_days(np.cumsum(rng.uniform(0.001, 0.05, 500)))
