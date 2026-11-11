"""Fit one candidate of EK Com and compare timing windows by the O-C scatter.

The file is ``data/EK_Com_sector_49.vot``. The scan is ``window_scan``: several
phase windows, and the O-C scatter of the timings in each. Run from the
repository root::

    python examples/scan_windows.py
"""

from __future__ import annotations

import logging
from pathlib import Path

import mikulasek_eb as eb
import mikulasek_eb_tess as tess

logger = logging.getLogger(__name__)
DATA = Path(__file__).resolve().parents[1] / "data" / "EK_Com_sector_49.vot"


def main() -> None:
    """Fit the first proposed candidate and log the window scan."""
    payload = tess.load_tess_votable(DATA)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    features = tess.extract_features(
        t_mjd, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    candidate = tess.propose_candidates(features)[0]
    out = eb.fit(
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        candidate.config,
        data_id=DATA.stem,
    )
    status = out.solution.stats["status"]
    if status != "ok":
        raise RuntimeError(f"{DATA.name}: fit status {status}")
    tess.window_scan(
        out.t_mjd,
        out.mag,
        out.sig,
        out.pipeline.model,
        out.pipeline.adopted,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
