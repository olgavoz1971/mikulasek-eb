"""Extensions beyond Mikulasek (2015), off by default: second profile component and eclipse weighting, with their acceptance checks."""

from __future__ import annotations

import logging

import numpy as np
from astropy.time import Time

from .model import EBModel, phase
from .tom import HALF_ORBIT_WINDOW, TIMING_BASELINE_ORDER
from .fitting import fit_global
from .timing import virtual_minima
from .tuning import wall_band

logger = logging.getLogger(__name__)


def refine_two_component(
    t,
    mag,
    sig,
    model,
    r,
    seg=None,
    order=1,
    eclipses=(1,),
    frac=0.3,
    widths=(0.5, 1.5),
    gammas=(0.7, 6.0),
    bounds=None,
):
    """EXTENSION (not in the paper): add a second profile component to the chosen eclipses.

    Starts from the single-component fit ``r`` and tries a small grid of starts (width factor
    times Gamma) for the new component. Judge the result by rms/noise, error sizes and
    timing stability, not by chi-square alone.

    Args:
        t (ndarray): Times (days).
        mag (ndarray): Magnitudes (mag).
        sig (ndarray): Magnitude errors (mag).
        model (EBModel): Model.
        r (FitResult): Single-component fit.
        seg (ndarray or None): Baseline segment id per point.
        order (int): Baseline polynomial degree.
        eclipses (tuple of int): Eclipses (1, 2) receiving the second component.
        frac (float): Starting fraction of the amplitude given to the new component.
        widths (tuple of float): Starting width factors for D<k>b relative to D<k>.
        gammas (tuple of float): Starting Gamma values for the new component.
        bounds (dict or None): Bounds passed to ``fit_global``.

    Returns:
        tuple: (best FitResult or None, table of (width factor, Gamma, chi2, rms in mmag, status)).
    """
    free_base = [n for n in r.names if n in r.p]
    free_eph = [n for n in r.names if n not in r.p]
    rows, best = [], None
    for wf in widths:
        for G2 in gammas:
            p = dict(r.p)
            newfree = []
            for k in eclipses:
                Dk = p.get(f"D{k}", p["D1"])
                Gk = p.get(f"G{k}", p["G1"])
                for nm, v in ((f"D{k}", Dk), (f"G{k}", Gk)):
                    p.setdefault(nm, v)  # make shape of this eclipse explicit
                Ak = p[f"A{k}"]
                p[f"A{k}"] = (1 - frac) * Ak
                p[f"A{k}b"] = frac * Ak
                p[f"D{k}b"] = wf * p[f"D{k}"]
                p[f"G{k}b"] = G2
                newfree += [f"D{k}", f"G{k}", f"A{k}b", f"D{k}b", f"G{k}b"]
            free = list(dict.fromkeys(free_base + newfree))
            try:
                rr = fit_global(
                    t, mag, sig, model, p, free, r.ephem, free_eph, seg=seg, order=order,
                    bounds=bounds,
                )
            except Exception as ex:
                logger.warning(
                    "two-component start (width %s, Gamma %s) failed: %s: %s",
                    wf,
                    G2,
                    type(ex).__name__,
                    ex,
                )
                rows.append((wf, G2, np.nan, np.nan, "error"))
                continue
            ok = np.isfinite(rr.chi2) and all(
                rr.p[f"D{k}"] > 0
                and rr.p[f"D{k}b"] > 0
                and rr.p[f"G{k}"] > 0
                and rr.p[f"G{k}b"] > 0
                for k in eclipses
            )
            rows.append((wf, G2, rr.chi2, 1e3 * rr.rms, "ok" if ok else "rejected"))
            if ok and (best is None or rr.chi2 < best.chi2):
                best = rr
    return best, rows


def eclipse_reweighted_sigma(
    t: np.ndarray,
    sig: np.ndarray,
    ephem: dict,
    weight: float,
    half: float,
) -> np.ndarray:
    """Scale catalogue sigma down in both eclipse windows (higher chi-squared weight).

    Args:
        t (ndarray): Times in the optimiser system (MJD in this project).
        sig (ndarray): Catalogue magnitude uncertainties.
        ephem (dict): Ephemeris (M0, P0, optional Pdot).
        weight (float): Factor >= 1. Points in eclipse count as weight squared.
        half (float): Phase half-width around 0 and 0.5 (cycles).

    Returns:
        ndarray: Effective sigma for the shape-only refit.

    Raises:
        ValueError: If ``weight`` is not finite and at least 1.
    """
    if not np.isfinite(weight) or weight < 1.0:
        raise ValueError(f"eclipse_weight must be >= 1, got {weight!r}.")
    th = phase(t, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0))
    ph = th - np.round(th)
    d2 = (ph - 0.5) - np.round(ph - 0.5)
    in_ecl = (np.abs(ph) < half) | (np.abs(d2) < half)
    sig_w = np.asarray(sig, dtype=float).copy()
    sig_w[in_ecl] = sig_w[in_ecl] / weight
    logger.info(
        "eclipse weight=%g on %s/%s points with |phase|<%g of 0 or 0.5",
        weight,
        int(np.count_nonzero(in_ecl)),
        len(t),
        half,
    )
    return sig_w


def two_component_degenerate(p: dict, k: int, depth: float) -> bool:
    """True if A and Ab have opposite sign and each exceeds the eclipse depth.

    Args:
        p (dict): Template including optional ``A{k}b``.
        k (int): Eclipse index (1 or 2).
        depth (float): Single-component depth of that eclipse (mag).

    Returns:
        bool: Degenerate pair.
    """
    a, ab = p[f"A{k}"], p[f"A{k}b"]
    return (a * ab < 0.0) and (abs(a) > abs(depth)) and (abs(ab) > abs(depth))


def _mjd_to_jd(t_mjd: float) -> float:
    """Convert one MJD value to absolute Julian Date.

    Args:
        t_mjd (float): Modified Julian Date.

    Returns:
        float: Julian Date (J0 = 0).
    """
    return float(Time(t_mjd, format="mjd").jd)


def _fold_phase_from_t0(t_jd: float, t0_jd: float, period: float) -> float:
    """Phase of one JD relative to folding epoch T0.

    Args:
        t_jd (float): Absolute Julian Date.
        t0_jd (float): TABLE PARAM epoch (absolute JD).
        period (float): Orbital period (days).

    Returns:
        float: Phase in [0, 1).
    """
    return float(np.mod((t_jd - t0_jd) / period, 1.0))


def _two_component_timing_ok(
    t_mjd: np.ndarray,
    mag: np.ndarray,
    sig: np.ndarray,
    model: EBModel,
    r_single,
    r_two,
    sigma_max: float,
) -> tuple[bool, str]:
    """Compare both-walls virtual minima for single- vs two-component templates.

    Args:
        t_mjd (ndarray): Times in MJD.
        mag (ndarray): Magnitudes.
        sig (ndarray): Uncertainties.
        model (EBModel): Phenomenological model.
        r_single: Single-component ``FitResult``.
        r_two: Two-component ``FitResult``.
        sigma_max (float): Largest accepted median |dT| / combined sigma.

    Returns:
        tuple[bool, str]: ``timing_ok`` and a log line for check 4.
    """
    wb = wall_band(model, r_two.p)
    vm_s = virtual_minima(
        t_mjd,
        mag,
        sig,
        model,
        r_single.p,
        r_single.ephem,
        window=HALF_ORBIT_WINDOW,
        order=TIMING_BASELINE_ORDER,
        wall_band=wb,
    )
    vm_t = virtual_minima(
        t_mjd,
        mag,
        sig,
        model,
        r_two.p,
        r_two.ephem,
        window=HALF_ORBIT_WINDOW,
        order=TIMING_BASELINE_ORDER,
        wall_band=wb,
    )
    common = np.intersect1d(vm_s["E"], vm_t["E"])
    if len(common) == 0:
        return False, "check4 no common both-walls events"
    ds = np.array([vm_s["T"][vm_s["E"] == e][0] for e in common])
    dt = np.array([vm_t["T"][vm_t["E"] == e][0] for e in common])
    comb = np.sqrt(
        np.array(
            [
                vm_s["sigma"][vm_s["E"] == e][0] ** 2
                + vm_t["sigma"][vm_t["E"] == e][0] ** 2
                for e in common
            ]
        )
    )
    d_t = dt - ds
    med_ratio = float(np.median(np.abs(d_t) / comb))
    msg = (
        f"check4 both-walls N={len(common)} "
        f"median dt={86400.0 * np.median(d_t):.1f} s "
        f"rms dt={86400.0 * d_t.std():.1f} s "
        f"median |dt|/comb.sigma={med_ratio:.2f}"
    )
    timing_ok = med_ratio <= sigma_max
    return timing_ok, msg
