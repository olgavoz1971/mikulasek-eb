"""Fit a light curve the way a user of the package does, and look at the result.

``fit`` applies the parameter bounds and writes a status on the solution. ``time``
returns eclipse times only when that status is ``ok``. This example uses the same
files as the tests of that behaviour, and saves the profile figure for each one.

Run from the repository root::

    MPLBACKEND=Agg python examples/inspect_fit.py

Figures go to ``examples/out/inspect_fit/``.
"""

from __future__ import annotations

import dataclasses
import logging
from pathlib import Path

import mikulasek_eb as eb
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.plotting import plot_global_fit

logger = logging.getLogger(__name__)
HERE = Path(__file__).resolve().parent
OUT = HERE / "out" / "inspect_fit"
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
SPOT_FILES = (
    GENERATED / "results/detached/detached_none_nooconnell_SP_L000_SP_L000_1800_c1800.vot",
    GENERATED / "results/detached/detached_none_nooconnell_SP_L180_SP_L180_1800_c1800.vot",
)
XY_BOO = HERE.parent / "data" / "lc_discovery_TIC_347405567_TESS_SPOC_sector_50.vot"


def spot_config(payload: dict, t_mjd) -> eb.Config:
    """Profile 1 with H2, as in the test of the 1800 s spot curves.

    The exposure and the gap split are measured on this file. The baseline is a
    straight line and a second component on the primary is allowed.

    Args:
        payload (dict): Light curve as returned by ``load_tess_votable``.
        t_mjd (ndarray): Times in MJD.

    Returns:
        Config: Configuration passed to ``fit``.
    """
    features = eb.extract_features(
        t_mjd, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    data = eb.propose_candidates(features)[0].config.data
    return eb.Config(
        eb.ModelSpec(profile=1, harmonics=(2,), ocon=True),
        dataclasses.replace(data, baseline_order=1),
        # Play, L000 only in spirit: the core default stays 0.5. The collapsed first
        # component there has a relative error of about 2.9, so 3 lets that trial through.
        eb.FitPolicy(allow_two_component=True, two_comp_max_rel_err=3.0),
        eb.TimingRecipe(window=0.12),
    )


def save_figure(out: eb.FitOutput, title: str, path: Path, timescale: str | None) -> None:
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
        out.t_mjd, out.mag, out.sig, model, adopted.p, adopted.ephem,
        out.pipeline.seg, out.solution.config.data.baseline_order,
    )
    plot_global_fit(
        eb.mjd_to_jd(out.t_mjd), out.mag, out.sig, yhat, template, theta, model,
        out.pipeline.single_component.p, adopted.p, out.pipeline.seg, title,
        timescale=timescale, save_path=str(path),
    )
    logger.info("figure %s", path)


def report(out: eb.FitOutput, name: str) -> None:
    """Log the status, and time the curve when the status allows it.

    Args:
        out (FitOutput): Result of ``fit``.
        name (str): Name used in the log.
    """
    stats = out.solution.stats
    logger.info(
        "%s status=%s rms=%.2f mmag", name, stats["status"], 1e3 * stats["rms_mag"]
    )
    for reason in stats["status_reasons"]:
        logger.info("  reason: %s", reason)
    for warning in stats["status_warnings"]:
        logger.info("  warning: %s", warning)
    try:
        timed = eb.time(out.t_mjd, out.mag, out.sig, out.solution)
    except eb.UnreliableSolution as ex:
        logger.info("%s not timed: %s", name, ex)
        return
    logger.info(
        "%s timed: %d primary, %d secondary",
        name, len(timed.timings["primary"]["T"]), len(timed.timings["secondary"]["T"]),
    )


def fit_file(path: Path, config: eb.Config) -> eb.FitOutput:
    """Load one VOTable and fit it.

    Args:
        path (Path): Light curve.
        config (Config): Configuration.

    Returns:
        FitOutput: Result of ``fit``.
    """
    payload = eb.load_tess_votable(path)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    out = eb.fit(
        t_mjd, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"],
        config, data_id=path.stem,
    )
    report(out, path.stem)
    save_figure(
        out,
        f"{path.stem}  {out.solution.stats['status']}  "
        f"rms {1e3 * out.solution.stats['rms_mag']:.2f} mmag",
        OUT / f"{path.stem}.png",
        payload.get("timescale"),
    )
    return out


def inspect_spots() -> None:
    """Fit the two 1800 s curves whose primary used to run away."""
    for path in SPOT_FILES:
        if not path.is_file():
            raise FileNotFoundError(path)
        payload = eb.load_tess_votable(path)
        t_mjd = eb.absolute_jd_to_mjd(payload["t"])
        fit_file(path, spot_config(payload, t_mjd))


def inspect_xy_boo() -> None:
    """Let ``select_model`` choose a model for XY Boo, and draw every candidate.

    A candidate whose status is not ``ok`` stays in the table and is left out of
    the choice.
    """
    if not XY_BOO.is_file():
        raise FileNotFoundError(XY_BOO)
    payload = eb.load_tess_votable(XY_BOO)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    result = eb.select_model(
        t_mjd, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"],
        data_id=XY_BOO.stem,
    )
    for row in result.table:
        logger.info("XY Boo candidate %s: %s", row.name, row.status)
    if result.best is None:
        logger.info("XY Boo: no candidate was chosen")
    else:
        logger.info("XY Boo: chosen %s", result.best.name)
    for name, out in result.fits.items():
        report(out, f"XY Boo {name}")
        save_figure(
            out,
            f"XY Boo {name}  {out.solution.stats['status']}  "
            f"rms {1e3 * out.solution.stats['rms_mag']:.2f} mmag",
            OUT / f"xy_boo_{name}.png",
            payload.get("timescale"),
        )


def main() -> None:
    """Fit the two spot curves and XY Boo, and save a figure for each fit."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    inspect_spots()
    inspect_xy_boo()


if __name__ == "__main__":
    main()
