"""Fit AR Aur with the published Table 1 model, as the paper test does.

The light curve is ``data/AR_Aur_sector_19.vot``. The period that is fitted
is the period of the file. The test is ``tests/tess/test_paper_comparison.py``.

Run from the repository root::

    python examples/fit_ar_aur.py

Writes ``examples/out/ar_aur/``.
"""

from __future__ import annotations

import logging

from fit_paper_star import main

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main("AR_Aur")
