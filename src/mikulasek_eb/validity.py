"""Verdict on a finished fit: is the adopted template a believable eclipse model?

A converged least-squares fit is not a good fit. On under-sampled eclipses the solver can
return ``success`` with a template that has no eclipse at the data (a needle of
Gamma -> 0 and A -> infinity), and the rms on the samples stays modest. This module
looks at the finished template and the data and gives one of three statuses:

``ok``
    No check failed.
``bound-active``
    A parameter or a derived width is outside the range of a real eclipse
    (Gamma, FWHM, depth). The numbers describe a numerical artefact, not a star.
``unresolved``
    The checks could not confirm the template: it does not follow the data in an
    eclipse core, the exposure integral has not converged, or there is no data in a core.

Only ``ok`` can be timed (``api.time`` refuses the rest). The limits below are design
choices of this package, not values from Mikulasek (2015); the paper sets no bounds.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

from .errors import UnreliableSolution
from .model import PROFILES, EBModel
from .tuning import initial_guess, noise_floor, out_of_eclipse_mask, residuals

logger = logging.getLogger(__name__)

STATUS_OK = "ok"
STATUS_BOUND_ACTIVE = "bound-active"
STATUS_UNRESOLVED = "unresolved"

# Gamma of a real eclipse. Paper Table 1 has 0.62 to 1.7; Fig. 1 draws 1/2. Total or
# contact eclipses with a flat bottom need a large Gamma, hence the high cap.
GAMMA_MIN = 0.3
GAMMA_MAX = 30.0
# Largest full width at half depth of one component (cycles): a primary eclipse of a
# circular binary cannot be wider than about half an orbit.
FWHM_MAX = 0.5
# Smallest full width at half depth, as a fraction of the exposure in phase.
FWHM_MIN_OVER_EXPOSURE = 0.5
# Depth of one component against the folded depth of the data.
DEPTH_MAX_OVER_FOLDED = 1.5
# Core of an eclipse for the agreement test, as a fraction of the data FWHM (each side).
CORE_HALF_OVER_FWHM = 0.5
# The template may miss the core mean by this fraction of the folded depth.
CORE_MISS_OVER_DEPTH = 0.2
CORE_MISS_OVER_SE = 5.0
CORE_MIN_POINTS = 3
# Curvature range for a warning only (paper Table 1: -0.34 to +0.10).
C_WARN = (-0.6, 0.3)
# rms over the noise floor above which a warning says the residuals are structured.
RMS_NOISE_WARN = 1.5
FINE_GRID_N = 4001
FWHM_GRID_N = 40001
SUBSAMPLE_FACTOR = 4


@dataclass
class FitVerdict:
    """Outcome of ``assess_fit``.

    Attributes:
        status (str): ``ok``, ``bound-active`` or ``unresolved``.
        reasons (list[str]): Why the status is not ``ok``; empty when it is.
        warnings (list[str]): Notes that do not change the status.
    """

    status: str
    reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """True when the status is ``ok``."""
        return self.status == STATUS_OK

    def to_stats(self) -> dict:
        """Entries for ``ModelSolution.stats``.

        Returns:
            dict: ``status``, ``status_reasons``, ``status_warnings``.
        """
        return dict(
            status=self.status,
            status_reasons=list(self.reasons),
            status_warnings=list(self.warnings),
        )


def verdict_from_stats(stats: dict) -> FitVerdict:
    """Rebuild a verdict from ``ModelSolution.stats``.

    Args:
        stats (dict): Statistics of a solution.

    Returns:
        FitVerdict: The stored verdict.

    Raises:
        UnreliableSolution: If the solution carries no status (it was made by an older
            version, or by code that skipped the assessment). Refit it.
    """
    if "status" not in stats:
        raise UnreliableSolution(
            "The solution has no fit status (made before the validity check). Refit it."
        )
    return FitVerdict(
        status=stats["status"],
        reasons=list(stats.get("status_reasons", [])),
        warnings=list(stats.get("status_warnings", [])),
    )


def require_ok(stats: dict) -> None:
    """Refuse a solution whose status is not ``ok``.

    Args:
        stats (dict): Statistics of a solution.

    Raises:
        UnreliableSolution: If the status is missing or is not ``ok``.
    """
    verdict = verdict_from_stats(stats)
    if not verdict.ok:
        raise UnreliableSolution(
            f"Fit status is '{verdict.status}': " + "; ".join(verdict.reasons)
        )


def components(model: EBModel, p: dict) -> list[dict]:
    """List the profile copies of a template, with their own shape parameters.

    Args:
        model (EBModel): Model that defines the profile.
        p (dict): Template parameters.

    Returns:
        list[dict]: One entry per copy: ``name``, ``centre`` (0 or 0.5), ``A``, ``D``,
        ``G``, ``C``.
    """
    out = []
    for k in (1, 2):
        out.append(
            dict(
                name=f"eclipse {k}",
                centre=0.5 * (k - 1),
                A=float(p[f"A{k}"]),
                D=float(p.get(f"D{k}", p["D1"])),
                G=float(p.get(f"G{k}", p["G1"])),
                C=float(p.get(f"C{k}", p.get("C1", 0.0))),
                signed=False,
            )
        )
        if f"A{k}b" in p:
            out.append(
                dict(
                    name=f"eclipse {k} second component",
                    centre=0.5 * (k - 1),
                    A=float(p[f"A{k}b"]),
                    D=float(p[f"D{k}b"]),
                    G=float(p[f"G{k}b"]),
                    C=0.0,
                    signed=True,
                )
            )
    return out


def component_fwhm(model: EBModel, comp: dict) -> float:
    """Full width at half depth of one profile copy, measured on a fine grid.

    Args:
        model (EBModel): Model that defines the profile.
        comp (dict): One entry of ``components``.

    Returns:
        float: FWHM in cycles (grid step 2.5e-5). A needle narrower than the step
        returns the step.
    """
    phi = np.linspace(-0.5, 0.5, FWHM_GRID_N)
    with np.errstate(all="ignore"):
        values = PROFILES[model.profile](phi, D=comp["D"], G=comp["G"], C=comp["C"])
    peak = float(values[len(values) // 2])
    if not np.isfinite(peak) or peak <= 0.0:
        return 0.0
    above = np.isfinite(values) & (values >= 0.5 * peak)
    return max(float(np.count_nonzero(above)), 1.0) * (phi[1] - phi[0])


def _limit_reasons(model: EBModel, p: dict, depth_ref: float, exposure_phase: float | None) -> list[str]:
    """Limits on Gamma, FWHM and depth of each profile copy.

    Args:
        model (EBModel): Model.
        p (dict): Template parameters.
        depth_ref (float): Folded depth of the data (mag).
        exposure_phase (float | None): Exposure in cycles, or None.

    Returns:
        list[str]: One sentence per violated limit.
    """
    reasons = []
    for comp in components(model, p):
        name = comp["name"]
        if not GAMMA_MIN <= comp["G"] <= GAMMA_MAX:
            reasons.append(
                f"{name}: Gamma={comp['G']:.3g} outside [{GAMMA_MIN}, {GAMMA_MAX}]"
            )
        amp = comp["A"]
        limit = DEPTH_MAX_OVER_FOLDED * depth_ref
        if (amp < 0.0 and not comp["signed"]) or abs(amp) > limit:
            reasons.append(
                f"{name}: depth A={amp:.4g} mag outside the range allowed by the "
                f"folded depth {depth_ref:.3f} mag"
            )
        if comp["D"] <= 0.0:
            reasons.append(f"{name}: width D={comp['D']:.3g} is not positive")
            continue
        fwhm = component_fwhm(model, comp)
        if fwhm > FWHM_MAX:
            reasons.append(f"{name}: FWHM={fwhm:.3g} cycles is wider than {FWHM_MAX}")
        if exposure_phase is not None and fwhm < FWHM_MIN_OVER_EXPOSURE * exposure_phase:
            reasons.append(
                f"{name}: FWHM={fwhm:.3g} cycles is narrower than "
                f"{FWHM_MIN_OVER_EXPOSURE} of the exposure ({exposure_phase:.3g} cycles)"
            )
    return reasons


def _core_reasons(res: np.ndarray, th: np.ndarray, g: dict, noise: float) -> list[str]:
    """Check that the template follows the data in each eclipse core.

    A template with the primary removed leaves the whole eclipse depth in the residual of
    the core. This is tested per minimum, so a healthy secondary cannot hide it.

    Args:
        res (ndarray): Residuals after template and baseline (mag).
        th (ndarray): Phase (cycles).
        g (dict): Folded start: ``A1``, ``A2``, ``fwhm``.
        noise (float): Out-of-eclipse noise (mag).

    Returns:
        list[str]: One sentence per minimum that fails.
    """
    reasons = []
    half = CORE_HALF_OVER_FWHM * float(g["fwhm"])
    for k, centre in ((1, 0.0), (2, 0.5)):
        d = th - centre
        d = d - np.round(d)
        core = np.abs(d) < half
        n = int(np.count_nonzero(core))
        if n < CORE_MIN_POINTS:
            reasons.append(
                f"eclipse {k}: only {n} points in the core, the template cannot be checked"
            )
            continue
        mean = float(np.mean(res[core]))
        se = noise / np.sqrt(n)
        depth = max(float(g[f"A{k}"]), 0.0)
        if abs(mean) > max(CORE_MISS_OVER_DEPTH * depth, CORE_MISS_OVER_SE * se):
            reasons.append(
                f"eclipse {k}: the template misses the data in the core by "
                f"{mean:.4g} mag (folded depth {depth:.3f} mag, {n} points)"
            )
    return reasons


def _exposure_reason(model: EBModel, p: dict, noise: float) -> str | None:
    """Check that the exposure integral has converged.

    The template is evaluated with ``nsub`` and with ``4 nsub`` sub-samples on a fine
    phase grid. A spike narrower than the sub-step gives a comb that differs between the
    two.

    Args:
        model (EBModel): Model with exposure integration.
        p (dict): Template parameters.
        noise (float): Out-of-eclipse noise (mag).

    Returns:
        str | None: A sentence when the two differ by more than the noise, else None.
    """
    if model.exposure_days is None or model.nsub == 1:
        return None
    fine = EBModel(
        profile=model.profile,
        harmonics=model.harmonics,
        ocon=model.ocon,
        kterm=model.kterm,
        exposure_days=model.exposure_days,
        nsub=SUBSAMPLE_FACTOR * model.nsub,
        period_days=model.period_days,
    )
    grid = np.linspace(-0.5, 0.5, FINE_GRID_N)
    with np.errstate(all="ignore"):
        diff = float(np.max(np.abs(model(grid, p) - fine(grid, p))))
    if not np.isfinite(diff) or diff > noise:
        return (
            f"the exposure integral has not converged: nsub={model.nsub} and "
            f"{SUBSAMPLE_FACTOR * model.nsub} differ by {diff:.4g} mag (noise {noise:.4g} mag)"
        )
    return None


def assess_fit(
    t_mjd: np.ndarray,
    mag: np.ndarray,
    sig: np.ndarray,
    period_days: float,
    model: EBModel,
    result,
    *,
    seg: np.ndarray | None = None,
    order: int = 1,
    start: dict | None = None,
    noise: float | None = None,
) -> FitVerdict:
    """Decide whether the adopted template is a believable eclipse model.

    Args:
        t_mjd (ndarray): Times used in the fit (MJD).
        mag (ndarray): Magnitudes.
        sig (ndarray): Catalogue errors (mag).
        period_days (float): Period used to fold the data (days).
        model (EBModel): Fitted model.
        result (FitResult): Adopted fit.
        seg (ndarray | None): Baseline segment ids.
        order (int): Baseline polynomial degree.
        start (dict | None): Folded start with ``A1``, ``A2``, ``fwhm``. None measures it
            from the data.
        noise (float | None): Out-of-eclipse noise (mag). None measures it from the data.

    Returns:
        FitVerdict: Status, reasons and warnings.
    """
    g = initial_guess(t_mjd, mag, period_days) if start is None else start
    if noise is None:
        ephem0 = dict(M0=g["M0"], P0=period_days, Pdot=0.0)
        noise = noise_floor(t_mjd, mag, mask=out_of_eclipse_mask(t_mjd, ephem0, g["fwhm"]))
    depth_ref = max(float(g["A1"]), float(g["A2"]))
    exposure_phase = (
        None if model.exposure_days is None else model.exposure_days / model.period_days
    )
    bound_reasons = _limit_reasons(model, result.p, depth_ref, exposure_phase)

    unresolved = []
    with np.errstate(all="ignore"):
        res, th = residuals(t_mjd, mag, sig, model, result, seg=seg, order=order)
    if not np.all(np.isfinite(res)):
        unresolved.append("the residuals are not finite")
    else:
        unresolved += _core_reasons(res, th, g, noise)
    comb = _exposure_reason(model, result.p, noise)
    if comb:
        unresolved.append(comb)

    warnings = []
    for key in ("C1", "C2"):
        if key in result.p and not C_WARN[0] <= result.p[key] <= C_WARN[1]:
            warnings.append(f"{key}={result.p[key]:.3g} outside {C_WARN}")
    if result.rms > RMS_NOISE_WARN * noise:
        warnings.append(
            f"rms/noise={result.rms / noise:.2f} exceeds {RMS_NOISE_WARN}: the template "
            "leaves structure in the residuals (status ok means numerically sound, not a good fit)"
        )
    if result.rank_deficient:
        warnings.append("the Jacobian is rank deficient: parameter errors are unreliable")

    if bound_reasons:
        verdict = FitVerdict(STATUS_BOUND_ACTIVE, bound_reasons + unresolved, warnings)
    elif unresolved:
        verdict = FitVerdict(STATUS_UNRESOLVED, unresolved, warnings)
    else:
        verdict = FitVerdict(STATUS_OK, [], warnings)
    if verdict.ok:
        logger.info("fit status ok%s", f" ({'; '.join(warnings)})" if warnings else "")
    else:
        logger.error("fit status %s: %s", verdict.status, "; ".join(verdict.reasons))
    return verdict
