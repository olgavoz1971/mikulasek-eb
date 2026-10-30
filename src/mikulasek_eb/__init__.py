"""
mikulasek_eb -- phenomenological light-curve model of eclipsing binaries.

Independent re-implementation from the published description (NOT the author's code):
    Mikulasek, Z. 2015, A&A 584, A8  (doi:10.1051/0004-6361/201425244)

Equations used
    Eq. (13)  eclipse profile  A (1 + C phi^2/D^2) {1 - [1 - exp(1 - cosh(phi/D))]^G}
    Eq. (15)  proximity terms  H_k cos(2 pi k theta);  O'Connell term  S1 sin(2 pi theta)
    Eq. (17)  quadratic phase function  theta = v - Pdot v^2/2, v=(t-M0)/P0 ;
              Theta(E) = M0 + P0 E + 0.5 P0 Pdot E^2
    Eq. (9)   per-interval time shifts (O-C) of a fixed template -> "virtual minima"
    Eq. (14)   optional quartic term K x^4 in profile 1 (EBModel(kterm=True))
    Table 2   eclipse profiles 1 to 7 (4 = Andronov-type)

Conventions
    * light curves in magnitudes; eclipse depth A > 0 means fainter
    * E integer = primary minimum, E = n + 0.5 = secondary (circular orbit assumed)
    * use times relative to a reference (t = BJD - T_ref) to avoid float64 loss,
      and put M0 near the middle of the data to decorrelate M0 and P0.
Units (the numerical core uses plain floats on purpose; the boundary adds the units)
    times t, M0, P0 and all timing errors: days (the pipeline passes MJD; the ToM files are
    written in JD, TDB); phase, D, theta: cycles; Gamma, C, Pdot: dimensionless; A, H, S1,
    y, sigma, rms: magnitudes.  astropy.units is applied where the data enter (volightcurve
    reader, time scale conversions), not inside the least-squares residual functions.

Two timing modes
    fit_minimum / virtual_minima : template fixed (optionally scaled in A, D, G),
        free shift delta [d] + per-event polynomial baseline.  window=0.25 -> whole
        phase around the event; small window (e.g. 0.05-0.1) -> local-minimum mode.
    fit_global : all template + ephemeris parameters simultaneously.


This is an independent re-implementation of the published method, not the author's code.
Plotting lives in ``mikulasek_eb.plotting`` (needs matplotlib) and is not imported here.
"""

from __future__ import annotations

from .model import (
    PROFILE_SHAPE,
    PROFILES,
    EBModel,
    phase,
    predicted_time,
    prof1,
    prof2,
    prof3,
    prof4,
    prof5,
    prof6,
    prof7,
)
from .exposure import estimate_cadence_days, estimate_exposure_days, suggest_nsub
from .fitting import (
    FitResult,
    fit_global,
    fit_global_iterative,
    correlation_matrix,
    strong_correlations,
    curve_band,
    ephemeris_at_epoch,
)
from .timing import (
    MinResult,
    ShiftFit,
    fit_minimum,
    virtual_minima,
    fit_ephemeris,
    fit_global_shifts,
)
from .tuning import (
    noise_floor,
    out_of_eclipse_mask,
    D_from_fwhm,
    initial_guess,
    multistart_global,
    residuals,
    binned_residuals,
    clip_timings,
    window_scan,
    wall_band,
    oc_diagnostics,
)
from .extensions import refine_two_component, eclipse_reweighted_sigma
from .api import FitOutput, TimingOutput, fit, time, validate_lightcurve
from .errors import (
    ConfigError,
    FitCancelled,
    FitDidNotConverge,
    InvalidLightCurve,
    MikulasekError,
    NoValidSolution,
)
from .solution import ModelSolution
from .features import Features, FeatureSettings, extract_features
from .borrow import GROUPS, ParameterPlan, fit_borrowed
from .config import (
    Config,
    DataSetup,
    FitConfig,
    FitPolicy,
    ModelSpec,
    TimingRecipe,
)
from .pipeline import FitPipelineResult, run_staged_fit
from .tom import (
    time_both_walls,
    time_with_recipe,
    write_timing_ascii,
    select_eclipses,
    clean_timings,
    choose_local_window,
    mjd_to_jd,
)
from .io import load_tess_votable
from .helpers import absolute_jd_to_mjd

__all__ = [
    "Features",
    "FeatureSettings",
    "extract_features",
    "GROUPS",
    "ParameterPlan",
    "fit_borrowed",
    "fit",
    "time",
    "FitOutput",
    "TimingOutput",
    "validate_lightcurve",
    "ModelSolution",
    "MikulasekError",
    "FitCancelled",
    "InvalidLightCurve",
    "NoValidSolution",
    "Config",
    "ConfigError",
    "DataSetup",
    "FitPolicy",
    "ModelSpec",
    "TimingRecipe",
    "estimate_cadence_days",
    "estimate_exposure_days",
    "suggest_nsub",
    "PROFILES",
    "EBModel",
    "phase",
    "predicted_time",
    "prof1",
    "prof2",
    "prof3",
    "prof4",
    "prof5",
    "prof6",
    "prof7",
    "PROFILE_SHAPE",
    "FitResult",
    "FitDidNotConverge",
    "fit_global",
    "fit_global_iterative",
    "correlation_matrix",
    "strong_correlations",
    "curve_band",
    "ephemeris_at_epoch",
    "MinResult",
    "ShiftFit",
    "fit_minimum",
    "virtual_minima",
    "fit_ephemeris",
    "fit_global_shifts",
    "noise_floor",
    "out_of_eclipse_mask",
    "D_from_fwhm",
    "initial_guess",
    "multistart_global",
    "residuals",
    "binned_residuals",
    "clip_timings",
    "window_scan",
    "wall_band",
    "oc_diagnostics",
    "refine_two_component",
    "eclipse_reweighted_sigma",
    "FitConfig",
    "FitPipelineResult",
    "run_staged_fit",
    "time_both_walls",
    "time_with_recipe",
    "write_timing_ascii",
    "select_eclipses",
    "clean_timings",
    "choose_local_window",
    "mjd_to_jd",
    "load_tess_votable",
    "absolute_jd_to_mjd",
]
