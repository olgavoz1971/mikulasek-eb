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

from .config import Config, FitConfig, TimingRecipe
from .errors import InvalidLightCurve
from .pipeline import FitPipelineResult, run_staged_fit
from .solution import ModelSolution, package_version, plain
from .tom import time_with_recipe

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
            E (from T0), T (absolute JD), sigma (days), delta, n, rchi2.
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
        n_input=len(t_mjd),
        rank_deficient=adopted.rank_deficient,
        two_component="A1b" in adopted.p,
    )
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
        len(t_mjd),
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
        InvalidLightCurve: If the input is invalid.
        ConfigError: If the recipe is invalid.
    """
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
