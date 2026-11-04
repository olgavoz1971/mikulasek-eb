"""Fit one spotted curve from candidate selection, and save the figures.

The harmonic orders, the O'Connell term and the baseline are chosen on this file.
Nothing is copied from the spot-free curve.

Run from the repository root::

    MPLBACKEND=Agg python examples/fit_from_scratch.py

Figures go to ``examples/out/fit_from_scratch/``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import mikulasek_eb as eb
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.plotting import plot_global_fit

logger = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent
OUT = HERE / "out" / "fit_from_scratch"
CURVE = Path(
    "/home/voz/projects/UPJS/EB_synthetic_lc/results/n01/detached/"
    "detached_none_nooconnell_SP_L315_SP_L315_120_c120.vot"
)


def save_figure(
    out: eb.FitOutput, title: str, path: Path, timescale: str | None
) -> None:
    """Save the four-panel profile figure for one fit.

    Args:
        out (FitOutput): Result of ``fit``.
        title (str): Figure title.
        path (Path): File to write.
        timescale (str | None): Time-scale name from the light curve, for the axis.
    """
    adopted = out.pipeline.adopted
    model = out.solution.build_model()
    yhat, template, theta = predicted_magnitude(
        out.t_mjd,
        out.mag,
        out.sig,
        model,
        adopted.p,
        adopted.ephem,
        out.pipeline.seg,
        out.solution.config.data.baseline_order,
    )
    plot_global_fit(
        eb.mjd_to_jd(out.t_mjd),
        out.mag,
        out.sig,
        yhat,
        template,
        theta,
        model,
        out.pipeline.single_component.p,
        adopted.p,
        out.pipeline.seg,
        title,
        timescale=timescale,
        save_path=str(path),
    )
    logger.info("figure %s", path)


def report(out: eb.FitOutput, name: str) -> None:
    """Log the adopted model, its status, and the times when the status allows them.

    Args:
        out (FitOutput): Result of ``fit``.
        name (str): Name used in the log.
    """
    stats = out.solution.stats
    model = out.solution.config.model
    logger.info(
        "%s harmonics=%s O'Connell=%s status=%s rms=%.2f mmag two-component=%s",
        name,
        model.harmonics,
        model.ocon,
        stats["status"],
        1e3 * stats["rms_mag"],
        out.solution.stats["two_component"],
    )
    for reason in stats["status_reasons"]:
        logger.info("  reason: %s", reason)
    for warning in stats["status_warnings"]:
        logger.info("  warning: %s", warning)
    trial = stats.get("two_component_trial")
    if trial:
        logger.info(
            "  two-component %s (single %.2f mmag, trial %s)",
            trial["reason"],
            trial["single_rms_mmag"],
            "not run" if trial["rms_mmag"] is None else f"{trial['rms_mmag']:.2f} mmag",
        )
        for gate in trial["failed_gates"]:
            logger.info("    failed: %s", gate)
    try:
        timed = eb.time(out.t_mjd, out.mag, out.sig, out.solution)
    except eb.UnreliableSolution as ex:
        logger.info("%s not timed: %s", name, ex)
        return
    logger.info(
        "%s timed: %d primary, %d secondary",
        name,
        len(timed.timings["primary"]["T"]),
        len(timed.timings["secondary"]["T"]),
    )


def main() -> None:
    """Select a model for the 120 s curve at spot longitude 315 degrees."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if not CURVE.is_file():
        raise FileNotFoundError(CURVE)
    OUT.mkdir(parents=True, exist_ok=True)
    payload = eb.load_tess_votable(CURVE)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    result = eb.select_model(
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        data_id=CURVE.stem,
    )
    for cand in result.candidates:
        logger.info("candidate %s: %s", cand.name, "; ".join(cand.reasons))
    for row in result.table:
        fitted = result.fits.get(row.name)
        rms = None if fitted is None else 1e3 * fitted.solution.stats["rms_mag"]
        logger.info(
            "fitted %s: status=%s rms=%s mmag",
            row.name,
            row.status,
            "none" if rms is None else f"{rms:.2f}",
        )
    if result.best is None:
        logger.info("no candidate was chosen")
    else:
        logger.info("chosen %s", result.best.name)
    for name, out in result.fits.items():
        report(out, name)
        save_figure(
            out,
            f"{name}  {out.solution.stats['status']}  "
            f"rms {1e3 * out.solution.stats['rms_mag']:.2f} mmag",
            OUT / f"{name}.png",
            payload.get("timescale"),
        )


if __name__ == "__main__":
    main()
