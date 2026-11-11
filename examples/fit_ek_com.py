"""Fit EK Com, sector 49, the way the regression test does.

The settings are those of ``working_scripts/fit_ek_com_s49.py``, which
``tests/tess/test_regression_real.py`` imports. The light curve is
``data/EK_Com_sector_49.vot``.

Run from the repository root::

    python examples/fit_ek_com.py

Writes ``examples/out/ek_com/fit.png`` and ``trend.png``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import mikulasek_eb as eb
import mikulasek_eb_tess as tess
from draw_fit import save_figures

logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "EK_Com_sector_49.vot"
OUT = ROOT / "examples" / "out" / "ek_com"

# Contact binary at 10 min cadence. The out-of-eclipse scatter is the
# ellipsoidal variation, so the usual rms/noise gate would never try the
# second component. Eclipse points are weighted threefold.
FIT_CONFIG = eb.FitConfig(
    profile=1,
    harmonics=(2, 3),
    ocon=True,
    baseline_order=0,
    gap_split_days=1.0,
    allow_two_component=True,
    rms_noise_two_comp_threshold=0.0,
    eclipse_weight=3.0,
    eclipse_phase_half=0.12,
)


def run() -> tuple[dict, eb.FitOutput]:
    """Load EK Com and fit it.

    Returns:
        tuple[dict, FitOutput]: The VOTable payload and the fit.
    """
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
    return payload, out


def main() -> None:
    """Fit EK Com and write the two figures."""
    payload, out = run()
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
