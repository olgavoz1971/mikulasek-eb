"""Fit one Table 1 star the way ``tests/tess/test_paper_comparison.py`` does.

The published model is laid on the TESS light curve and refitted in steps
(variants A to D). Variant D frees the shape. The comparison figure is the
one the test writes; ``fit.png`` and ``trend.png`` show variant D alone.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import mikulasek_eb_tess as tess
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.tuning import initial_guess

from draw_fit import save_figures

logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]


def main(star_key: str) -> None:
    """Run the published-model fit of one star and write the figures.

    Args:
        star_key (str): Manifest key, one of ``EK_Com``, ``V477_Lyr``,
            ``AV_Del``, ``AR_Aur``.

    Raises:
        KeyError: If the star has no Table 1 template.
    """
    sys.path.insert(0, str(ROOT / "working_scripts"))
    from paper_forward_check import (  # noqa: WPS433 (path is set just above)
        BASELINE_ORDER,
        TEMPLATE_BUILDERS,
        baseline_segments,
        fit_variants,
        load_star,
        plot_forward_check,
    )

    if star_key not in TEMPLATE_BUILDERS:
        raise KeyError(f"No Table 1 template for {star_key}.")
    entry, paper, tess_lambda = load_star(star_key)
    model, params, _perr = TEMPLATE_BUILDERS[star_key](paper, tess_lambda)
    payload = tess.load_tess_votable(ROOT / entry["votable"])
    t_mjd = absolute_jd_to_mjd(payload["t"])
    seg = baseline_segments(t_mjd)
    guess = initial_guess(t_mjd, payload["mag"], payload["period"])
    fits = fit_variants(
        t_mjd,
        payload["mag"],
        payload["sig"],
        model,
        params,
        star_key,
        seg,
        guess,
        payload["period"],
    )
    adopted = fits["D"]
    folder = ROOT / "examples" / "out" / star_key.lower()
    folder.mkdir(parents=True, exist_ok=True)
    comparison = folder / "paper_comparison.png"
    shape = ", ".join(
        f"{k}={params[k]}" for k in ("D1", "G1", "C1", "C2") if k in params
    )
    plot_forward_check(
        t_mjd,
        payload["mag"],
        payload["sig"],
        model,
        fits,
        seg,
        f"{star_key}: Table 1 model on TESS (published shape {shape})",
        comparison,
    )
    logger.info("wrote %s", comparison)
    logger.info("variant D rms %.2f mmag", 1e3 * adopted.rms)
    save_figures(
        folder,
        f"{payload['title']}\nvariant D",
        t_mjd=t_mjd,
        mag=payload["mag"],
        sig=payload["sig"],
        model=model,
        p_single=adopted.p,
        p_adopted=adopted.p,
        ephem=adopted.ephem,
        seg=seg,
        order=BASELINE_ORDER,
        timescale=payload["timescale"],
    )
