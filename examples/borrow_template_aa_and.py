"""Borrow a model between light curves of one star (AA And, TIC 455790537).

Reference: SPOC sector 16, 2-min cadence (full staged fit). The model is then borrowed for
the 30-min cadence light curve of the same sector and for the 30-min light curve of
sector 17: shape fixed, amplitudes, harmonics and the zero phase M0 fitted, own exposure
integration, own baseline. One PNG per light curve goes to ``results/``; the reference
solution and the two borrowed solutions are saved as JSON next to them. The 30-min
timings are compared with the 2-min ones.

Run::

    python examples/borrow_template_aa_and.py
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from astropy.time import Time

import mikulasek_eb as eb
from mikulasek_eb.config import DataSetup, FitPolicy, ModelSpec
from mikulasek_eb.exposure import estimate_exposure_days
from mikulasek_eb.helpers import absolute_jd_to_mjd, predicted_magnitude
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.plotting import plot_global_fit

logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
STEM = "lc_discovery_TIC_455790537_"
REFERENCE = "SPOC_sector_16"  # 2-min cadence
TARGETS = ("TESS_SPOC_sector_16", "TESS_SPOC_sector_17")  # 30-min cadence
FREE_GROUPS = ("amplitudes", "harmonics", "zero_phase")
# Baseline of the borrowed fits: cubic per segment, split at gaps longer than this (days).
BASELINE_ORDER = 3
GAP_SPLIT_DAYS = 0.5
TIMING_RECIPE = eb.TimingRecipe(require_contact_walls=False, clip_outliers=False)


def load(name: str) -> tuple[dict, np.ndarray]:
    """Read a light curve VOTable.

    Args:
        name (str): File name part after ``lc_discovery_TIC_455790537_``.

    Returns:
        tuple: (payload dict from ``load_tess_votable``, times in MJD).
    """
    payload = load_tess_votable(DATA / f"{STEM}{name}.vot")
    return payload, absolute_jd_to_mjd(payload["t"])


def draw(name: str, payload: dict, out: eb.FitOutput, order: int, seg, ref_p: dict):
    """Save the diagnostic figure of one fitted light curve.

    Args:
        name (str): Light-curve name (used in the title and file name).
        payload (dict): VOTable payload (title, time scale, period, epoch).
        out (eb.FitOutput): Fit result with the arrays the fit used.
        order (int): Baseline polynomial degree of that fit.
        seg (ndarray | None): Baseline segments of that fit.
        ref_p (dict): Reference template parameters, drawn as the "earlier template".
    """
    sol = out.solution
    model = sol.build_model()
    yhat, f_template, theta = predicted_magnitude(
        out.t_mjd, out.mag, out.sig, model, sol.p, sol.ephem, seg, order
    )
    t_jd = np.asarray(Time(out.t_mjd, format="mjd").jd, dtype=float)
    path = RESULTS / f"borrow_aa_and_{name}.png"
    plot_global_fit(
        t_jd,
        out.mag,
        out.sig,
        yhat,
        f_template,
        theta,
        model,
        ref_p,
        sol.p,
        seg,
        f"{payload['title']} ({name})  rms {1e3 * sol.stats['rms_mag']:.2f} mmag",
        timescale=payload["timescale"],
        save_path=str(path),
        note="Blue line: reference template (2-min fit). Red: this light curve.",
    )
    logger.info("figure %s", path)


def timing_difference(a: eb.TimingOutput, b: eb.TimingOutput) -> None:
    """Log the differences of the timings of the same events in two timing outputs.

    Args:
        a (eb.TimingOutput): Reference timings (2-min).
        b (eb.TimingOutput): Timings to compare (30-min).
    """
    for kind in ("primary", "secondary"):
        ta, tb = a.timings[kind], b.timings[kind]
        common = np.intersect1d(ta["E"], tb["E"])
        ia = [list(ta["E"]).index(e) for e in common]
        ib = [list(tb["E"]).index(e) for e in common]
        diff = (np.asarray(tb["T"])[ib] - np.asarray(ta["T"])[ia]) * 86400.0
        logger.info(
            "%s: %d common events, 30-min minus 2-min: mean %.1f s, scatter %.1f s",
            kind,
            len(common),
            diff.mean(),
            diff.std(),
        )


def main() -> None:
    """Fit the reference, borrow its model for the other light curves, plot and compare."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    RESULTS.mkdir(exist_ok=True)

    # 1. Reference: full staged fit of the 2-min light curve.
    payload, t = load(REFERENCE)
    config = eb.Config(
        ModelSpec(profile=1, harmonics=(1, 2, 3), ocon=True),
        DataSetup(baseline_order=0, gap_split_days=None),
        FitPolicy(allow_two_component=True),
    )
    ref_out = eb.fit(
        t,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        config,
        data_id=REFERENCE,
    )
    reference = ref_out.solution
    reference.save(RESULTS / f"borrow_aa_and_{REFERENCE}.json")
    draw(
        REFERENCE,
        payload,
        ref_out,
        ref_out.pipeline.order,
        ref_out.pipeline.seg,
        reference.p,
    )
    # Reference timings, with the same cleaned-up baseline setting as the borrowed fits.
    ref_data = DataSetup(baseline_order=BASELINE_ORDER, gap_split_days=GAP_SPLIT_DAYS)
    plan = eb.ParameterPlan(FREE_GROUPS)
    ref_refit = eb.fit_borrowed(
        t, payload["mag"], payload["sig"], reference, ref_data, plan, data_id=REFERENCE
    )
    ref_timing = eb.time(
        ref_refit.t_mjd, ref_refit.mag, ref_refit.sig, ref_refit.solution, TIMING_RECIPE
    )

    # 2. Borrow the model for the 30-min light curves.
    for name in TARGETS:
        payload, t = load(name)
        data = DataSetup(
            baseline_order=BASELINE_ORDER,
            gap_split_days=GAP_SPLIT_DAYS,
            exposure_days=estimate_exposure_days(t),
        )
        out = eb.fit_borrowed(
            t, payload["mag"], payload["sig"], reference, data, plan, data_id=name
        )
        out.solution.save(RESULTS / f"borrow_aa_and_{name}.json")
        seg = np.cumsum(np.r_[0, np.diff(out.t_mjd) > GAP_SPLIT_DAYS])
        draw(name, payload, out, BASELINE_ORDER, seg, reference.p)
        timing = eb.time(out.t_mjd, out.mag, out.sig, out.solution, TIMING_RECIPE)
        logger.info("%s: %s", name, {k: len(v["T"]) for k, v in timing.timings.items()})
        if name == "TESS_SPOC_sector_16":
            timing_difference(ref_timing, timing)


if __name__ == "__main__":
    main()
