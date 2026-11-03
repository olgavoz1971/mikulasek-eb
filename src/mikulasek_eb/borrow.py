"""Borrowing a model between light curves of the same star.

A reference ``ModelSolution`` (for example from a high-cadence light curve) supplies the
model structure and the starting values. A ``ParameterPlan`` names the parameter groups
that are fitted to the new light curve; all other parameters stay at the reference values.
The new light curve keeps its own exposure, baseline and clipping settings (``DataSetup``).

Groups:
    ``shape``       D, Gamma, C, K of every eclipse and second component.
    ``amplitudes``  A of every eclipse and second component.
    ``harmonics``   H<k> and the O'Connell term S1.
    ``zero_phase``  M0, the epoch of the primary minimum.
    ``period``      P0.

Borrowed values are fixed, not "soft starts": the fit is linear-baseline plus the free
groups only. Parameter errors are reported for the free parameters only.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass

import numpy as np

from .api import FitOutput, validate_lightcurve
from .config import Config, DataSetup
from .errors import ConfigError, InvalidLightCurve
from .exposure import suggest_nsub
from .fitting import fit_global
from .model import PROFILES, EBModel
from .solution import ModelSolution, package_version, plain
from .helpers import out_of_eclipse_phase_mask
from .tuning import residuals
from .validity import assess_fit

logger = logging.getLogger(__name__)

GROUPS = ("shape", "amplitudes", "harmonics", "zero_phase", "period")

_PARAMETER_PATTERNS = {
    "shape": re.compile(r"^[DGCK]\d+b?$"),
    "amplitudes": re.compile(r"^A\d+b?$"),
    "harmonics": re.compile(r"^(H\d+|S1)$"),
}
_EPHEMERIS_GROUPS = {"zero_phase": "M0", "period": "P0"}


@dataclass(frozen=True)
class ParameterPlan:
    """Which parameter groups are fitted; the others keep the reference values.

    Attributes:
        free (tuple of str): Names from ``GROUPS`` that are fitted.
    """

    free: tuple = ("amplitudes", "zero_phase")

    def validate(self) -> ParameterPlan:
        """Check the group names.

        Returns:
            ParameterPlan: ``self``.

        Raises:
            ConfigError: If a group name is unknown or repeated.
        """
        bad = [g for g in self.free if g not in GROUPS]
        if bad:
            raise ConfigError(f"Unknown parameter groups {bad}; use {list(GROUPS)}.")
        if len(set(self.free)) != len(self.free):
            raise ConfigError(f"Repeated parameter groups in {list(self.free)}.")
        return self

    def to_dict(self) -> dict:
        """Plain-data form.

        Returns:
            dict: ``{"free": [...]}``.
        """
        return {"free": list(self.free)}

    @classmethod
    def from_dict(cls, data: dict) -> ParameterPlan:
        """Build a plan from ``to_dict`` output.

        Args:
            data (dict): Dictionary with the key ``free``.

        Returns:
            ParameterPlan: Validated plan.
        """
        return cls(free=tuple(data["free"])).validate()

    def free_names(self, params: dict) -> tuple[list, list]:
        """Split the parameter names into template and ephemeris names to fit.

        Args:
            params (dict): Template parameters of the reference solution.

        Returns:
            tuple: (template names, ephemeris names).
        """
        template = [
            n
            for n in params
            if any(
                g in self.free and rx.match(n) for g, rx in _PARAMETER_PATTERNS.items()
            )
        ]
        ephemeris = [e for g, e in _EPHEMERIS_GROUPS.items() if g in self.free]
        return template, ephemeris


def _narrowest_fwhm_cycles(profile: int, p: dict) -> float:
    """Narrowest full width at half depth of the main eclipse profiles (cycles).

    Args:
        profile (int): Table 2 profile number.
        p (dict): Template parameters (eclipse 2 shape defaults to eclipse 1).

    Returns:
        float: FWHM in cycles.
    """
    fn = PROFILES[profile]
    phi = np.linspace(0.0, 0.25, 20001)
    widths = []
    for k in (1, 2):
        D = p.get(f"D{k}", p["D1"])
        G = p.get(f"G{k}", p.get("G1"))
        C = p.get(f"C{k}", p.get("C1", 0.0))
        K = p.get(f"K{k}", p.get("K1", 0.0))
        y = np.asarray(fn(phi, D, G, C, K), float)
        widths.append(2.0 * phi[np.argmax(y < 0.5 * y[0])])
    return float(min(widths))


def fit_borrowed(
    t_mjd,
    mag,
    sig,
    reference: ModelSolution,
    data: DataSetup,
    plan: ParameterPlan | None = None,
    *,
    data_id: str | None = None,
) -> FitOutput:
    """Fit a light curve with a model borrowed from a reference solution.

    The model structure, ephemeris numbering (T0 and period of the reference input) and
    the policy and timing recipe come from the reference; the exposure, baseline and gap
    settings come from ``data``. Out-of-eclipse outliers are clipped as in the staged
    pipeline (``policy.clip_nsig``), then the fit is repeated once.

    Args:
        t_mjd (ndarray): Times (MJD), strictly increasing.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        reference (ModelSolution): Solution to borrow from.
        data (DataSetup): Properties of this light curve.
        plan (ParameterPlan | None): Groups to fit; default amplitudes and zero phase.
        data_id (str | None): Identifier of the data for provenance.

    Returns:
        FitOutput: New solution (``pipeline`` is None) and the arrays used by the fit.

    Raises:
        InvalidLightCurve: If the input is invalid.
        ConfigError: If the plan or the data setup is invalid.
        FitDidNotConverge: If the solver fails.
    """
    plan = (plan or ParameterPlan()).validate()
    prov = reference.provenance
    validate_lightcurve(t_mjd, mag, sig, prov["period_input_days"], prov["epoch_jd"])
    t, y, s = (np.asarray(a, dtype=float) for a in (t_mjd, mag, sig))
    config = Config(
        reference.config.model,
        data,
        reference.config.policy,
        reference.config.timing,
        notes=reference.config.notes,
    ).validate()
    names, eph_names = plan.free_names(reference.p)
    if not names and not eph_names:
        raise ConfigError("The plan frees no parameter present in the reference.")
    if len(t) < 3 * (len(names) + len(eph_names)) + 10:
        raise InvalidLightCurve("Too few points for the free parameters of the plan.")

    spec = config.model
    period = prov["period_input_days"]
    if data.exposure_days is None:
        model = EBModel(spec.profile, spec.harmonics, spec.ocon, spec.kterm)
        exposure = None
    else:
        fwhm_days = _narrowest_fwhm_cycles(spec.profile, reference.p) * period
        nsub = data.exposure_nsub or suggest_nsub(data.exposure_days, fwhm_days)
        model = EBModel(
            spec.profile,
            spec.harmonics,
            spec.ocon,
            spec.kterm,
            exposure_days=data.exposure_days,
            nsub=nsub,
            period_days=period,
        )
        exposure = dict(days=data.exposure_days, nsub=nsub)
        logger.info(
            "borrowed fit: exposure %.1f s in %d sub-samples",
            86400.0 * data.exposure_days,
            nsub,
        )
    seg = (
        None
        if data.gap_split_days is None
        else np.cumsum(np.r_[0, np.diff(t) > data.gap_split_days])
    )
    order = data.baseline_order
    ephem0 = reference.ephem

    def run(tt, yy, ss, gg):
        """Fit the free groups to the given points."""
        return fit_global(
            tt,
            yy,
            ss,
            model,
            reference.p,
            names,
            ephem0,
            eph_names,
            seg=gg,
            order=order,
        )

    r = run(t, y, s, seg)
    res, th = residuals(t, y, s, model, r, seg=seg, order=order)
    oo = out_of_eclipse_phase_mask(th)
    mad = 1.4826 * np.median(np.abs(res[oo] - np.median(res[oo])))
    keep = np.ones(len(t), dtype=bool)
    keep[oo] = np.abs(res[oo]) < config.policy.clip_nsig * max(mad, 1e-8)
    n_clip = int((~keep).sum())
    logger.info("borrowed fit: clipped %d points (mad %.2f mmag)", n_clip, 1e3 * mad)
    if n_clip:
        t, y, s = t[keep], y[keep], s[keep]
        seg = None if seg is None else seg[keep]
        r = run(t, y, s, seg)

    stats = dict(
        chi2=r.chi2,
        dof=r.dof,
        rchi2=r.rchi2,
        rms_mag=r.rms,
        n_points=len(t),
        n_input=len(t_mjd),
        n_clipped=n_clip,
        rank_deficient=r.rank_deficient,
        two_component="A1b" in r.p,
    )
    stats.update(
        assess_fit(t, y, s, period, model, r, seg=seg, order=order).to_stats()
    )
    solution = ModelSolution(
        config=config,
        params=plain(r.p),
        errors=plain(r.err),
        ephemeris=plain(r.ephem),
        stats=plain(stats),
        exposure=exposure,
        provenance=dict(
            data_id=data_id,
            epoch_jd=prov["epoch_jd"],
            period_input_days=period,
            package_version=package_version(),
            borrowed_from=reference.solution_hash(),
            plan=plan.to_dict(),
        ),
    )
    logger.info(
        "borrowed solution %s: rms %.2f mmag, %d of %d points",
        solution.solution_hash(),
        1e3 * r.rms,
        len(t),
        len(t_mjd),
    )
    return FitOutput(solution, t, y, s, None)
