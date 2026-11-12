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
    written in JD). JD and MJD differ by the constant JD_TO_MJD in helpers. phase, D, theta:
    cycles; Gamma, C, Pdot: dimensionless; A, H, S1, y, sigma, rms: magnitudes.

Layers
    mikulasek_eb            the core: model, fits, timing. numpy, scipy, astropy only.
    mikulasek_eb.selection  features, candidates and the choice by evidence, no survey numbers.
    mikulasek_eb.compare    model-free reference measurements of a minimum.
    mikulasek_eb_tess       TESS adapter: reader, exposure rule and calibrated settings.

Two timing modes
    fit_minimum / virtual_minima : template fixed (optionally scaled in A, D, G),
        free shift delta [d] + per-event polynomial baseline.  window=0.25 -> whole
        phase around the event; small window (e.g. 0.05-0.1) -> local-minimum mode.
    fit_global : all template + ephemeris parameters simultaneously.


This is an independent re-implementation of the published method, not the author's code.
The root exports only the everyday surface below. Everything else is public in its own module:
``model`` (profiles), ``fitting``, ``timing``, ``tuning``, ``tom``, ``extensions``, ``exposure``,
``validity``, ``borrow``, ``helpers``.
Plotting lives in ``mikulasek_eb.plotting`` (needs matplotlib) and is not imported here.
"""

from __future__ import annotations

from .api import FitOutput, TimingOutput, fit, time, validate_lightcurve
from .borrow import ParameterPlan, fit_borrowed
from .config import (
    Config,
    DataSetup,
    FitConfig,
    FitPolicy,
    ModelSpec,
    TimingRecipe,
)
from .errors import (
    ConfigError,
    FitCancelled,
    FitDidNotConverge,
    InvalidLightCurve,
    MikulasekError,
    NoValidSolution,
    UnreliableSolution,
)
from .helpers import absolute_jd_to_mjd, mjd_to_jd
from .model import PROFILES, EBModel, phase
from .pipeline import FitPipelineResult, run_staged_fit
from .solution import ModelSolution

__all__ = [
    # run
    "fit",
    "time",
    "validate_lightcurve",
    "fit_borrowed",
    "run_staged_fit",
    # results
    "FitOutput",
    "TimingOutput",
    "ModelSolution",
    "FitPipelineResult",
    "ParameterPlan",
    # configuration
    "Config",
    "ModelSpec",
    "DataSetup",
    "FitConfig",
    "FitPolicy",
    "TimingRecipe",
    # model and time
    "EBModel",
    "PROFILES",
    "phase",
    "mjd_to_jd",
    "absolute_jd_to_mjd",
    # errors
    "MikulasekError",
    "ConfigError",
    "FitCancelled",
    "FitDidNotConverge",
    "InvalidLightCurve",
    "NoValidSolution",
    "UnreliableSolution",
]
