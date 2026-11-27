"""Host-facing functions: ``fit`` and ``time``.

Both are pure: they read their arguments and return new objects, with no global state, no
file access and no plotting, so a web server, a notebook or an agent can call them. Heavy
arrays stay in ``FitOutput`` and ``TimingOutput``; the light, JSON-serialisable parts are
``ModelSolution`` and ``TimingOutput.to_dict()``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

from .config import (
    MIN_POINTS_EACH_SIDE,
    MIN_POINTS_PER_EVENT,
    Config,
    FitConfig,
    TimingRecipe,
)
from .errors import InvalidLightCurve
from .pipeline import FitPipelineResult, run_staged_fit
from .solution import ModelSolution, package_version, plain
from .helpers import mjd_to_jd
from .model import phase, predicted_time
from .timing import fit_minimum
from .tom import time_with_recipe
from .validity import require_ok

logger = logging.getLogger(__name__)

# Eclipse types timed by default, each with its own output and rejection list.
ECLIPSE_TYPES = ("primary", "secondary")


@dataclass
class FitOutput:
    """Result of ``fit``.

    Attributes:
        solution (ModelSolution): Light, serialisable result.
        t_mjd (ndarray): Times of the points the fit used (after its out-of-eclipse clip).
        mag (ndarray): Magnitudes of those points.
        sig (ndarray): Their errors (after eclipse weighting, if requested).
        pipeline (FitPipelineResult | None): Full pipeline result for plots and diagnostics;
            None for a borrowed fit (``fit_borrowed``), which has no staged search.
    """

    solution: ModelSolution
    t_mjd: np.ndarray
    mag: np.ndarray
    sig: np.ndarray
    pipeline: FitPipelineResult | None


@dataclass
class TimingOutput:
    """Result of ``time``.

    Attributes:
        timings (dict): ``{"primary": {...}, "secondary": {...}}``; each holds arrays
            E (from T0), T (absolute JD), sigma (days), delta, n, rchi2, walls_ok,
            n_before, n_after, balanced.
        rejected (dict): For each eclipse type, the list of events that could not be timed
            (E, T_pred_jd, reason, detail).
        recipe (TimingRecipe): Recipe used.
        solution_hash (str): Hash of the solution whose template was used.
    """

    timings: dict
    rejected: dict
    recipe: TimingRecipe
    solution_hash: str = field(default="")

    def to_dict(self) -> dict:
        """JSON-compatible form (arrays become lists).

        Returns:
            dict: ``timings``, ``rejected``, ``recipe`` and ``solution_hash``.
        """
        import dataclasses

        return plain(
            dict(
                timings={
                    k: {n: np.asarray(v).tolist() for n, v in d.items()}
                    for k, d in self.timings.items()
                },
                rejected=self.rejected,
                recipe=dataclasses.asdict(self.recipe),
                solution_hash=self.solution_hash,
            )
        )


def validate_lightcurve(t_mjd, mag, sig, period_days: float, epoch_jd: float) -> None:
    """Structural checks of the input; no scientific heuristics.

    Args:
        t_mjd (ndarray): Times (MJD).
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        period_days (float): Approximate period (days).
        epoch_jd (float): Folding epoch T0 (absolute JD).

    Raises:
        InvalidLightCurve: If the arrays differ in length or are not 1-D, contain
            non-finite values, errors are not positive, times are not strictly
            increasing, there are fewer than 3 points, or the period or epoch is invalid.
    """
    arrays = [np.asarray(a, dtype=float) for a in (t_mjd, mag, sig)]
    if any(a.ndim != 1 for a in arrays):
        raise InvalidLightCurve("t, mag and sig must be 1-D arrays.")
    if len({a.shape[0] for a in arrays}) != 1:
        raise InvalidLightCurve("t, mag and sig must have the same length.")
    t, m, s = arrays
    if len(t) < 3:
        raise InvalidLightCurve("Need at least 3 points.")
    if not (
        np.all(np.isfinite(t)) and np.all(np.isfinite(m)) and np.all(np.isfinite(s))
    ):
        raise InvalidLightCurve(
            "t, mag and sig must be finite (remove NaN or inf first)."
        )
    if not np.all(s > 0):
        raise InvalidLightCurve("All magnitude errors must be positive.")
    if not np.all(np.diff(t) > 0):
        raise InvalidLightCurve(
            "Times must be strictly increasing (sort and drop duplicates)."
        )
    if not (np.isfinite(period_days) and period_days > 0):
        raise InvalidLightCurve(f"Invalid period {period_days!r}.")
    if not np.isfinite(epoch_jd):
        raise InvalidLightCurve(f"Invalid epoch {epoch_jd!r}.")


def solution_from_pipeline(
    result: FitPipelineResult,
    n_input: int,
    period_days: float,
    epoch_jd: float,
    *,
    data_id: str | None = None,
) -> ModelSolution:
    """Package a staged-pipeline result as a ``ModelSolution``.

    Args:
        result (FitPipelineResult): Result of ``run_staged_fit``.
        n_input (int): Number of input points before clipping.
        period_days (float): Input period (days) given to the pipeline.
        epoch_jd (float): Input epoch T0 (absolute JD) given to the pipeline.
        data_id (str | None): Identifier of the data for provenance.

    Returns:
        ModelSolution: Serialisable solution.
    """
    adopted = result.adopted
    model = result.model
    stats = dict(
        chi2=adopted.chi2,
        dof=adopted.dof,
        rchi2=adopted.rchi2,
        rms_mag=adopted.rms,
        noise_floor_mag=result.noise_floor,
        rms_over_noise=adopted.rms / result.noise_floor,
        n_points=len(result.t_mjd),
        n_input=n_input,
        rank_deficient=adopted.rank_deficient,
        two_component="A1b" in adopted.p,
        two_component_trial=result.two_component_trial,
    )
    if result.verdict is None:
        raise ValueError("The pipeline result has no fit verdict.")
    stats.update(result.verdict.to_stats())
    exposure = (
        None
        if model.exposure_days is None
        else dict(days=model.exposure_days, nsub=model.nsub)
    )
    solution = ModelSolution(
        config=result.config,
        params=plain(adopted.p),
        errors=plain(adopted.err),
        ephemeris=plain(adopted.ephem),
        stats=plain(stats),
        exposure=exposure,
        provenance=dict(
            data_id=data_id,
            epoch_jd=float(epoch_jd),
            period_input_days=float(period_days),
            package_version=package_version(),
        ),
    )
    logger.info(
        "solution %s: rms %.2f mmag, %d of %d points",
        solution.solution_hash(),
        1e3 * adopted.rms,
        len(result.t_mjd),
        n_input,
    )
    return solution


def fit(
    t_mjd,
    mag,
    sig,
    period_days: float,
    epoch_jd: float,
    config: Config | FitConfig,
    *,
    data_id: str | None = None,
    progress=None,
) -> FitOutput:
    """Fit the phenomenological model to one light curve.

    Args:
        t_mjd (ndarray): Times (MJD), strictly increasing.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        period_days (float): Approximate period (days), for example the catalogue value.
        epoch_jd (float): Folding epoch T0 as absolute JD (numbering origin of E).
        config (Config | FitConfig): Configuration.
        data_id (str | None): Identifier of the data (for provenance and cache keys).
        progress (callable | None): ``progress(stage, fraction)``; may raise
            ``FitCancelled`` to stop.

    Returns:
        FitOutput: Solution and the arrays the fit used.

    Raises:
        InvalidLightCurve: If the input is invalid or is not an eclipsing binary.
        ConfigError: If the configuration is invalid or unsupported.
        FitDidNotConverge: If the solver fails.
        NoValidSolution: If no acceptable solution was found.
        FitCancelled: If the progress hook cancelled the fit.
    """
    validate_lightcurve(t_mjd, mag, sig, period_days, epoch_jd)
    t_mjd, mag, sig = (np.asarray(a, dtype=float) for a in (t_mjd, mag, sig))
    result = run_staged_fit(
        t_mjd, mag, sig, period_days, epoch_jd, config, progress=progress
    )
    solution = solution_from_pipeline(
        result, len(t_mjd), period_days, epoch_jd, data_id=data_id
    )
    return FitOutput(solution, result.t_mjd, result.mag, result.sig, result)


def time(
    t_mjd,
    mag,
    sig,
    solution: ModelSolution,
    recipe: TimingRecipe | None = None,
    *,
    which: tuple = ECLIPSE_TYPES,
) -> TimingOutput:
    """Time the eclipses of a light curve with the template of a solution.

    Pass the arrays returned by ``fit`` (``FitOutput.t_mjd`` and so on) to reproduce the
    timings of that fit; raw arrays also contain the points the fit clipped.

    Args:
        t_mjd (ndarray): Times (MJD), strictly increasing.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        solution (ModelSolution): Provides the template and the T0 numbering.
        recipe (TimingRecipe | None): Timing method; None uses the recipe stored in the
            solution's configuration.
        which (tuple): Eclipse types to time, each with its own output and rejections.

    Returns:
        TimingOutput: Timings and rejected events per eclipse type.

    Raises:
        UnreliableSolution: If the solution status is not ``ok`` (``bound-active``,
            ``unresolved`` or missing). Nothing is timed from such a template.
        InvalidLightCurve: If the input is invalid.
        ConfigError: If the recipe is invalid.
    """
    require_ok(solution.stats)
    prov = solution.provenance
    validate_lightcurve(t_mjd, mag, sig, prov["period_input_days"], prov["epoch_jd"])
    recipe = solution.config.timing if recipe is None else recipe
    Config(solution.config.model, solution.config.data, timing=recipe).validate()
    t_mjd, mag, sig = (np.asarray(a, dtype=float) for a in (t_mjd, mag, sig))
    model = solution.build_model()
    timings, rejected = {}, {}
    for kind in which:
        rejected[kind] = []
        timings[kind] = time_with_recipe(
            t_mjd,
            mag,
            sig,
            model,
            solution,
            kind,
            prov["epoch_jd"],
            prov["period_input_days"],
            recipe,
            rejected=rejected[kind],
        )
    return TimingOutput(timings, rejected, recipe, solution.solution_hash())


@dataclass
class CycleTimingOutput:
    """Result of ``time_cycles``: one shift of the whole template per orbital cycle.

    Attributes:
        cycles (dict): Arrays ``E`` (integer cycle), ``delta`` (days), ``T_primary`` and
            ``T_secondary`` (absolute JD), ``sigma`` (days), ``n``, ``rchi2``.
        rejected (list): Cycles that could not be timed (E, reason, detail).
        min_points (int): Point floor used.
        min_each_side (int): Per-side floor used.
        order (int): Per-cycle baseline polynomial degree.
        solution_hash (str): Hash of the solution whose template was used.
    """

    cycles: dict
    rejected: list
    min_points: int
    min_each_side: int
    order: int = 1
    solution_hash: str = field(default="")


def time_cycles(
    t_mjd,
    mag,
    sig,
    solution: ModelSolution,
    min_points: int = MIN_POINTS_PER_EVENT,
    min_each_side: int = MIN_POINTS_EACH_SIDE,
    order: int = 1,
) -> CycleTimingOutput:
    """Time every orbital cycle with one shift of the whole template.

    This is a separate method from ``time``. There is no eclipse window. The points of
    cycle ``n`` are those with phase in ``(n - 0.25, n + 0.75)``, so the primary and the
    secondary of that cycle share one shift. The default floor is 4 points with 2 on
    each side of the cycle centre.

    Args:
        t_mjd (ndarray): Times (MJD), strictly increasing.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        solution (ModelSolution): Provides the template and the T0 numbering.
        min_points (int): Fewest points in a cycle.
        min_each_side (int): Fewest points on each side of the cycle centre.
        order (int): Per-cycle baseline polynomial degree.

    Returns:
        CycleTimingOutput: One timed cycle is one shift, used for both eclipses.

    Raises:
        UnreliableSolution: If the solution status is not ``ok``.
        InvalidLightCurve: If the input is invalid.
    """
    require_ok(solution.stats)
    prov = solution.provenance
    validate_lightcurve(t_mjd, mag, sig, prov["period_input_days"], prov["epoch_jd"])
    t_mjd, mag, sig = (np.asarray(a, dtype=float) for a in (t_mjd, mag, sig))
    ephem = solution.ephem
    model = solution.build_model()
    th = phase(t_mjd, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0))
    n_min = int(np.floor(float(np.min(th)) + 0.25))
    n_max = int(np.floor(float(np.max(th)) + 0.25 - 1e-9))
    rows = []
    rejected = []
    for n in range(n_min, n_max + 1):
        centre = n + 0.25
        offset = th - centre
        inside = np.abs(offset) < 0.5
        n_in = int(inside.sum())
        n_before = int(np.sum(offset[inside] < 0))
        n_after = n_in - n_before
        if n_in < min_points:
            rejected.append(
                dict(E=float(n), reason="too_few_points", detail=f"n={n_in}")
            )
            continue
        if n_before < min_each_side or n_after < min_each_side:
            rejected.append(
                dict(E=float(n), reason="one_side_uncovered", detail=f"n={n_in}")
            )
            continue
        try:
            fit = fit_minimum(
                t_mjd,
                mag,
                sig,
                model,
                solution.p,
                ephem,
                centre,
                window=0.5,
                order=order,
                free=(),
            )
        except Exception as exc:
            logger.warning("cycle %s fit failed: %s: %s", n, type(exc).__name__, exc)
            rejected.append(
                dict(
                    E=float(n),
                    reason="fit_failed",
                    detail=f"{type(exc).__name__}: {exc}",
                )
            )
            continue
        t_pri = predicted_time(
            float(n), ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0)
        )
        t_sec = predicted_time(
            n + 0.5, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0)
        )
        rows.append(
            dict(
                E=float(n),
                delta=float(fit.delta),
                T_primary=float(mjd_to_jd(t_pri + fit.delta)),
                T_secondary=float(mjd_to_jd(t_sec + fit.delta)),
                sigma=float(fit.sigma),
                n=int(fit.n),
                rchi2=float(fit.rchi2),
            )
        )
    if not rows:
        cycles = {
            k: np.array([])
            for k in ("E", "delta", "T_primary", "T_secondary", "sigma", "n", "rchi2")
        }
    else:
        cycles = {k: np.array([row[k] for row in rows]) for k in rows[0]}
    return CycleTimingOutput(
        cycles,
        rejected,
        min_points,
        min_each_side,
        order,
        solution.solution_hash(),
    )
