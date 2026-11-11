"""Fit AA And, SPOC sector 16 at 2 min, the way the regression test does.

The settings are those of ``working_scripts/fit_tic455790537_spoc_s16.py``,
which ``tests/tess/test_regression_real.py`` imports. The light curve is
``data/lc_discovery_TIC_455790537_SPOC_sector_16.vot``.

Run from the repository root::

    python examples/fit_aa_and.py

Writes ``examples/out/aa_and/fit.png`` and ``trend.png``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import mikulasek_eb as eb
import mikulasek_eb_tess as tess
from draw_fit import save_figures

logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "lc_discovery_TIC_455790537_SPOC_sector_16.vot"
OUT = ROOT / "examples" / "out" / "aa_and"

# Total primary, 2 min cadence. Harmonics (1, 2, 3) and the second component
# are allowed; the pipeline keeps the second component only if its checks pass.
FIT_CONFIG = eb.FitConfig(
    profile=1,
    harmonics=(1, 2, 3),
    ocon=True,
    baseline_order=0,
    gap_split_days=None,
    allow_two_component=True,
)


def main() -> None:
    """Fit AA And and write the two figures."""
    payload = tess.load_tess_votable(DATA)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    out = eb.fit(
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        FIT_CONFIG,
        data_id=DATA.stem,
    )
    pipeline = out.pipeline
    adopted = pipeline.adopted
    logger.info(
        "status %s, rms %.2f mmag",
        out.solution.stats["status"],
        1e3 * adopted.rms,
    )
    save_figures(
        OUT,
        f"{payload['title']}  P={payload['period']:.6f} d",
        t_mjd=pipeline.t_mjd,
        mag=pipeline.mag,
        sig=pipeline.sig,
        model=pipeline.model,
        p_single=pipeline.single_component.p,
        p_adopted=adopted.p,
        ephem=adopted.ephem,
        seg=pipeline.seg,
        order=pipeline.order,
        timescale=payload["timescale"],
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
