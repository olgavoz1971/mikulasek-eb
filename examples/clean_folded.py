"""Flag points far from the folded light curve of one packaged TESS file.

The file is ``data/v477_Lyr_sector_26.vot``. The rule is the robust sigma of
each phase bin. Run from the repository root::

    python examples/clean_folded.py
"""

from __future__ import annotations

import logging
from pathlib import Path

import mikulasek_eb_tess as tess

logger = logging.getLogger(__name__)
DATA = Path(__file__).resolve().parents[1] / "data" / "v477_Lyr_sector_26.vot"


def main() -> None:
    """Load the file and report how many points the folded rule flags."""
    payload = tess.load_tess_votable(DATA)
    flag, _, _ = tess.flag_folded_outliers(
        payload["t"], payload["mag"], payload["period"]
    )
    logger.info(
        "%s: %d of %d points beyond 5 robust sigma of their phase bin",
        DATA.name,
        int(flag.sum()),
        len(flag),
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
