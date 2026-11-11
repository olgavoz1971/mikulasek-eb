"""Fit V477 Lyr with the published Table 1 model, as the paper test does.

The light curve is ``data/v477_Lyr_sector_26.vot``. The test is
``tests/tess/test_paper_comparison.py``.

Run from the repository root::

    python examples/fit_v477_lyr.py

Writes ``examples/out/v477_lyr/``.
"""

from __future__ import annotations

import logging

from fit_paper_star import main

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main("V477_Lyr")
