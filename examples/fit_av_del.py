"""Fit AV Del with the published Table 1 model, as the paper test does.

The light curve is ``data/AV_Del_sector_81.vot``. The test is
``tests/tess/test_paper_comparison.py``.

Run from the repository root::

    python examples/fit_av_del.py

Writes ``examples/out/av_del/``.
"""

from __future__ import annotations

import logging

from fit_paper_star import main

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main("AV_Del")
