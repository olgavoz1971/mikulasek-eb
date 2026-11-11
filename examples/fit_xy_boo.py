"""Choose a model for XY Boo the way the selection test does.

``tests/tess/test_selection.py`` fits the candidates and keeps harmonics
(2, 3). The light curve is
``data/lc_discovery_TIC_347405567_TESS_SPOC_sector_50.vot``.

Run from the repository root::

    python examples/fit_xy_boo.py

Writes ``examples/out/xy_boo/fit.png`` and ``trend.png``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import mikulasek_eb as eb
import mikulasek_eb_tess as tess
from draw_fit import save_figures

logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "lc_discovery_TIC_347405567_TESS_SPOC_sector_50.vot"
OUT = ROOT / "examples" / "out" / "xy_boo"


def main() -> None:
    """Select the model for XY Boo and write the two figures of the winner."""
    payload = tess.load_tess_votable(DATA)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    result = tess.select_model(
        t_mjd, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    if result.best is None or result.output is None:
        raise eb.NoValidSolution(f"{DATA.name}: no candidate finished")
    out = result.output
    pipeline = out.pipeline
    adopted = pipeline.adopted
    logger.info(
        "chosen %s, status %s, rms %.2f mmag",
        result.best.name,
        out.solution.stats["status"],
        1e3 * adopted.rms,
    )
    save_figures(
        OUT,
        f"{payload['title']}\n{result.best.name}  {out.solution.stats['status']}",
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
