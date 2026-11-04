"""Fit a light curve the way a user of the package does, and look at the result.

``fit`` uses Levenberg-Marquardt unless ``FitPolicy.bounded_solver`` is set, then
writes a status on the solution. ``time`` returns eclipse times only when that
status is ``ok``. For the two 1800 s curves the example fits twice: first with
the default solver, and, when that fit is refused, again with the bounded solver.
For XY Boo it also forces a second eclipse component on the chosen model, because
the rms cut skips that trial, and refits the eclipse shape with extra weight on
the minima. This example uses the same files as the tests of that behaviour, and
saves a profile figure for each fit.

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
    GENERATED
    / "results/n01/detached/detached_none_nooconnell_SP_L000_SP_L000_1800_c1800.vot",
    GENERATED
    / "results/n01/detached/detached_none_nooconnell_SP_L180_SP_L180_1800_c1800.vot",
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


def save_figure(
    out: eb.FitOutput,
    title: str,
    path: Path,
    timescale: str | None,
    note: str | None = None,
) -> None:
    """Save the four-panel profile figure for one fit.

    Args:
        out (FitOutput): Result of ``fit``.
        title (str): Figure title.
        path (Path): File to write.
        timescale (str | None): Time-scale name from the light curve, for the axis.
        note (str | None): Subtitle. None leaves the plot to choose one.
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
        note=note,
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
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        config,
        data_id=path.stem,
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


def refit_bounded(path: Path, payload: dict, t_mjd, config: eb.Config) -> eb.FitOutput:
    """Refit one refused curve with the eclipse shape kept inside its box.

    The default solver is Levenberg-Marquardt. It can run away (a huge depth and a
    Gamma far below 0.3) and the verdict then refuses the fit, which is the result
    we want to see first. ``bounded_solver`` is a second attempt, not the default:
    the trust-region solver holds the depth, the width and Gamma inside the limits
    while it fits. The verdict is applied again, so a bounded fit can still be
    refused, for example when a component is wider than half a cycle.

    Args:
        path (Path): Light curve, used for the figure name.
        payload (dict): Loaded light curve.
        t_mjd (ndarray): Times in MJD.
        config (Config): The configuration of the refused fit.

    Returns:
        FitOutput: The bounded fit.
    """
    logger.info(
        "%s: default fit refused. Refitting with the bounded solver "
        "(FitPolicy.bounded_solver).",
        path.stem,
    )
    bounded = dataclasses.replace(
        config, policy=dataclasses.replace(config.policy, bounded_solver=True)
    )
    out = eb.fit(
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        bounded,
        data_id=path.stem + "_bounded",
    )
    report(out, path.stem + " bounded")
    save_figure(
        out,
        f"{path.stem} bounded  {out.solution.stats['status']}  "
        f"rms {1e3 * out.solution.stats['rms_mag']:.2f} mmag",
        OUT / f"{path.stem}_bounded.png",
        payload.get("timescale"),
    )
    return out


def inspect_spots() -> None:
    """Fit the two 1800 s curves whose primary used to run away.

    Each curve is fitted with the default solver. A fit that is not ``ok`` is
    fitted again with the bounded solver, and both figures are kept.
    """
    for path in SPOT_FILES:
        if not path.is_file():
            raise FileNotFoundError(path)
        payload = eb.load_tess_votable(path)
        t_mjd = eb.absolute_jd_to_mjd(payload["t"])
        config = spot_config(payload, t_mjd)
        out = fit_file(path, config)
        if out.solution.stats["status"] != "ok":
            refit_bounded(path, payload, t_mjd, config)


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
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
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
    forced_two_component(
        payload, t_mjd, None if result.best is None else result.best.config
    )
    weighted_minima(payload, t_mjd, None if result.best is None else result.best.config)


def weighted_minima(payload: dict, t_mjd, chosen: eb.Config | None) -> None:
    """Refit the eclipse shape with the minima counted more heavily.

    The catalogue-weight fit is done first, so the cosines and the ephemeris
    describe the bright variation. They are then frozen. Points within 0.12
    cycles of each minimum are given weight 3 (they count ninefold) and only
    the eclipse shape is fitted again. The figure shows that shape in red and
    the catalogue-weight shape in grey. The unweighted residuals decide whether
    the bottom got rounder. This is a look, not a change of the package default.

    Args:
        payload (dict): Loaded light curve.
        t_mjd (ndarray): Times in MJD.
        chosen (Config | None): The chosen configuration. None means nothing was chosen.
    """
    if chosen is None:
        logger.info("XY Boo: no model was chosen, so the minima were not reweighted")
        return
    policy = dataclasses.replace(
        chosen.policy, eclipse_weight=3.0, eclipse_phase_half=0.12
    )
    _fit_forced(
        payload,
        t_mjd,
        chosen,
        policy,
        "xy_boo_weighted_minima",
        note="One component. Minima weight 3; cosines and ephemeris frozen.",
    )


def forced_two_component(payload: dict, t_mjd, chosen: eb.Config | None) -> None:
    """Try a second eclipse component on the model ``select_model`` chose.

    The chosen XY Boo fit skips that trial: rms/noise is below 1.5, so the
    pipeline never fits the second copy. ``force_two_component`` fits it anyway.
    It skips the rms-gain gate and no other gate. The eye says both minima are
    rounder than the single component, so the trial is run twice: on the primary
    only, and on both minima. When the remaining gates reject the copy, one more
    fit loosens the relative-error gate so the trial can be drawn. That fit is a
    look, not a change of the package default.

    Args:
        payload (dict): Loaded light curve.
        t_mjd (ndarray): Times in MJD.
        chosen (Config | None): The chosen configuration. None means nothing was chosen.
    """
    if chosen is None:
        logger.info(
            "XY Boo: no model was chosen, so the second component was not forced"
        )
        return
    for eclipses, tag in (((1,), "primary"), ((1, 2), "both")):
        policy = dataclasses.replace(
            chosen.policy,
            allow_two_component=True,
            force_two_component=True,
            two_comp_eclipses=eclipses,
        )
        out = _fit_forced(payload, t_mjd, chosen, policy, f"xy_boo_forced_{tag}")
        adopted = out.solution.stats["two_component_trial"]["adopted"]
        if adopted:
            continue
        # The figure above is the single component. Loosen the determination
        # gates so the trial itself is drawn. This look is not a package default.
        loose = dataclasses.replace(
            policy, two_comp_max_rel_err=50.0, two_comp_min_agreeing_starts=1
        )
        _fit_forced(payload, t_mjd, chosen, loose, f"xy_boo_forced_{tag}_shown")


def _fit_forced(
    payload: dict,
    t_mjd,
    chosen: eb.Config,
    policy: eb.FitPolicy,
    stem: str,
    note: str | None = None,
) -> eb.FitOutput:
    """Fit one forced two-component configuration and save its figure.

    Args:
        payload (dict): Loaded light curve.
        t_mjd (ndarray): Times in MJD.
        chosen (Config): Model and data of the chosen candidate.
        policy (FitPolicy): Policy for this look.
        stem (str): Figure name without the suffix.
        note (str | None): Subtitle of the figure.
    """
    config = dataclasses.replace(chosen, policy=policy)
    out = eb.fit(
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        config,
        data_id=stem,
    )
    report(out, stem)
    save_figure(
        out,
        f"XY Boo {stem}  {out.solution.stats['status']}  "
        f"rms {1e3 * out.solution.stats['rms_mag']:.2f} mmag",
        OUT / f"{stem}.png",
        payload.get("timescale"),
        note=note,
    )
    return out


def main() -> None:
    """Fit the two spot curves and XY Boo, and save a figure for each fit."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    inspect_spots()
    inspect_xy_boo()


if __name__ == "__main__":
    main()
