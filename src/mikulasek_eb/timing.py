"""Per-eclipse times of minima (virtual minima, Eq. 9), O-C ephemeris fits and the simultaneous template plus shifts solution."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from .fitting import _Nuisance, _check_solution, _cov
from .model import phase, predicted_time

logger = logging.getLogger(__name__)


@dataclass
class MinResult:
    """Outcome of ``fit_minimum``.

    Attributes:
        E (float): Eclipse number.
        T (float): Time of minimum (days).
        sigma (float): Formal 1-sigma error of T (days).
        delta (float): Shift relative to the ephemeris prediction (days).
        n (int): Number of points used.
        rchi2 (float): Reduced chi-square of the event fit.
        scales (dict): Fitted multiplicative shape scales (A, D, G), if requested.
    """

    E: float
    T: float
    sigma: float
    delta: float
    n: int
    rchi2: float
    scales: dict = field(default_factory=dict)


# ------------------------------------------- single event / virtual minima
def fit_minimum(
    t,
    y,
    sigma,
    model,
    p,
    ephem,
    E,
    window=0.1,
    order=1,
    free=(),
    scan=0.1,
    rescale=True,
):
    """Time of one eclipse (number E) from points with |theta - E| < window.

    Model: y = poly_order(t) + F(theta_c(t) - delta/P_loc), Eq. (9). Free: the shift delta
    and optional multiplicative scales of this eclipse's depth 'A', width 'D', kurtosis 'G'
    (use for local-window fits when the shape varies).

    Args:
        t (ndarray): Times (days).
        y (ndarray): Magnitudes (mag).
        sigma (ndarray): Magnitude errors (mag).
        model (EBModel): Template model.
        p (dict): Template parameters held fixed (except the scales).
        ephem (dict): Ephemeris (M0, P0 in days; optional Pdot).
        E (float): Eclipse number.
        window (float): Half-width of the fitted phase window (cycles, at most 0.25).
        order (int): Per-event baseline polynomial degree.
        free (tuple of str): Shape scales to fit, subset of "A", "D", "G".
        scan (float): First scan the shift over +-scan cycles to avoid local minima (0 = off).
        rescale (bool): Rescale the covariance by the reduced chi-square.

    Returns:
        MinResult: Time of minimum with its formal error.

    Raises:
        FitDidNotConverge: If the solver does not converge.
    """
    t, y, sigma = (np.asarray(a, float) for a in (t, y, sigma))
    M0, P0, Pdot = ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0)
    th = phase(t, M0, P0, Pdot)
    sel = np.abs(th - E) < window
    ts, ys, ths = t[sel], y[sel], th[sel]
    w = 1.0 / sigma[sel]
    Pl = P0 * (1.0 + Pdot * E)
    k = 1 if abs(E - round(E)) < 1e-6 else 2
    base = dict(p)
    for kk in (1, 2):
        for nm in "DGCK":
            if f"{nm}1" in p:
                base.setdefault(f"{nm}{kk}", p[f"{nm}1"])
    scales = list(free)
    nui = _Nuisance(ts, w, None, order)

    def build(x):
        """Template parameters with the fitted shape scales applied.

        Args:
            x (ndarray): Solver vector, shift first then the scales.

        Returns:
            dict: Scaled template parameters.
        """
        q = dict(base)
        for nm, sc in zip(scales, x[1:]):
            for key in (f"{nm}{k}", f"{nm}{k}b"):
                if key in q:
                    q[key] = q[key] * sc
        return q

    def resid(x):
        """Weighted, baseline-projected residuals of the event fit.

        Args:
            x (ndarray): Solver vector.

        Returns:
            ndarray: Residuals (dimensionless).
        """
        return nui.project(w * (ys - model(ths - x[0] / Pl, build(x))))

    d0 = 0.0
    if scan > 0:
        grid = np.linspace(-scan, scan, 61) * Pl
        c = [np.sum(resid(np.r_[g, np.ones(len(scales))]) ** 2) for g in grid]
        d0 = grid[int(np.argmin(c))]
    x0 = np.r_[d0, np.ones(len(scales))]
    sol = least_squares(resid, x0, method="lm", x_scale="jac")
    _check_solution(sol, f"fit_minimum E={E}")
    dof = int(sel.sum()) - nui.npar - len(x0)
    cov, chi2, rchi2 = _cov(sol, dof, rescale)
    delta = float(sol.x[0])
    return MinResult(
        E,
        predicted_time(E, M0, P0, Pdot) + delta,
        float(np.sqrt(cov[0, 0])),
        delta,
        int(sel.sum()),
        rchi2,
        dict(zip(scales, sol.x[1:])),
    )


def virtual_minima(
    t,
    y,
    sigma,
    model,
    p,
    ephem,
    window=0.25,
    order=1,
    free=(),
    min_points=10,
    min_each_side=3,
    wall_band=None,
    min_wall=3,
    rejected=None,
    **kw,
):
    """Fit every eclipse covered by the data (window <= 0.25 so events never overlap).

    Events lacking enough points, or points on both sides of the predicted centre, are
    skipped and recorded, never dropped silently. If ``wall_band`` is given an event also
    needs at least ``min_wall`` points on both walls, because partially covered events rest
    on one wall and are biased by any wall-shape mismatch (see ``eb_tuning.wall_band``).
    Call separately for each instrument or filter.

    Args:
        t (ndarray): Times (days).
        y (ndarray): Magnitudes (mag).
        sigma (ndarray): Magnitude errors (mag).
        model (EBModel): Template model.
        p (dict): Template parameters.
        ephem (dict): Ephemeris (M0, P0 in days; optional Pdot).
        window (float): Half-width of the phase window (cycles, capped at 0.25).
        order (int): Per-event baseline polynomial degree.
        free (tuple of str): Shape scales to fit per event.
        min_points (int): Minimum points in the window.
        min_each_side (int): Minimum points on each side of the centre.
        wall_band (tuple or dict or None): (lo, hi) in |phase offset|, or per eclipse
            {1: (lo, hi), 2: (lo, hi)}.
        min_wall (int): Minimum points on each wall when ``wall_band`` is given.
        rejected (list or None): Receives one dict per skipped event:
            ``{"E": float, "reason": "too_few_points" | "one_side_uncovered" |
            "walls_uncovered" | "fit_failed", "detail": str}``.
        **kw: Passed to ``fit_minimum`` (``scan``, ``rescale``).

    Returns:
        dict: Arrays E, T (days), sigma (days), delta (days), n, rchi2.
    """
    window = min(window, 0.25)
    t = np.asarray(t, float)
    th = phase(t, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0))
    En = np.round(2 * th) / 2
    near = np.abs(th - En) < window
    rows = []
    for E in np.unique(En[near]):
        d = th[near & (En == E)] - E
        reason = None
        if len(d) < min_points:
            reason = "too_few_points"
        elif (d < 0).sum() < min_each_side or (d > 0).sum() < min_each_side:
            reason = "one_side_uncovered"
        elif wall_band is not None:
            lo, hi = (
                wall_band[1 if abs(E - round(E)) < 1e-6 else 2]
                if isinstance(wall_band, dict)
                else wall_band
            )
            if ((d > -hi) & (d < -lo)).sum() < min_wall or (
                (d > lo) & (d < hi)
            ).sum() < min_wall:
                reason = "walls_uncovered"
        if reason is not None:
            logger.debug(
                "event E=%s rejected: %s (n=%d, before=%d, after=%d)",
                E,
                reason,
                len(d),
                int((d < 0).sum()),
                int((d > 0).sum()),
            )
            if rejected is not None:
                rejected.append(dict(E=float(E), reason=reason, detail=f"n={len(d)}"))
            continue
        try:
            rows.append(
                fit_minimum(t, y, sigma, model, p, ephem, E, window, order, free, **kw)
            )
        except Exception as ex:
            logger.warning("event E=%s fit failed: %s: %s", E, type(ex).__name__, ex)
            if rejected is not None:
                rejected.append(
                    dict(
                        E=float(E),
                        reason="fit_failed",
                        detail=f"{type(ex).__name__}: {ex}",
                    )
                )
    return {
        k: np.array([getattr(r, k) for r in rows])
        for k in ("E", "T", "sigma", "delta", "n", "rchi2")
    }


# ------------------------------------------------------- O-C helper
def fit_ephemeris(T, sigma, E, deg=2):
    """Weighted polynomial ephemeris T(E).

    Args:
        T (ndarray): Times of minima (days).
        sigma (ndarray): Errors of T (days).
        E (ndarray): Eclipse numbers.
        deg (int): Polynomial degree.

    Returns:
        tuple: (coefficients in increasing order, residuals in days, weighted rms in days).
    """
    T, sigma, E = (np.asarray(a, float) for a in (T, sigma, E))
    V = np.vander(E, deg + 1, increasing=True)
    w = 1 / sigma
    c = np.linalg.lstsq(V * w[:, None], T * w, rcond=None)[0]
    res = T - V @ c
    return c, res, float(np.sqrt(np.sum((res * w) ** 2) / np.sum(w**2)))


@dataclass
class ShiftFit:
    """Outcome of ``fit_global_shifts``.

    Attributes:
        p, p_err (dict): Template parameters and errors.
        E, T, sigma, delta (ndarray): Eclipse numbers, times, joint-fit errors and shifts (days).
        cov_T (ndarray): Full covariance of the timings (days^2).
        sigma_fixed_template (ndarray): Errors of the same events with the template fixed (days).
        names, cov: Parameter names and joint covariance.
        chi2, dof, rchi2, rms: Fit statistics (rms in mag).
        n_unassigned (int): Points outside every accepted event.
    """

    p: dict
    p_err: dict
    E: np.ndarray
    T: np.ndarray
    sigma: np.ndarray
    delta: np.ndarray
    cov_T: np.ndarray  # full covariance matrix of the timings (days^2)
    sigma_fixed_template: (
        np.ndarray
    )  # same events, template held fixed (virtual_minima)
    names: list
    cov: np.ndarray
    chi2: float
    dof: int
    rchi2: float
    rms: float
    n_unassigned: int


def fit_global_shifts(
    t,
    y,
    sigma,
    model,
    p0,
    free,
    ephem,
    seg=None,
    order=1,
    window=0.25,
    min_points=10,
    min_each_side=3,
    wall_band=None,
    min_wall=3,
    rescale=True,
):
    """Simultaneous solution of template parameters and one time shift per eclipse.

    Implements Eq. (9) with Sect. 5.2 principle 1 and Sect. 5.1 error propagation. The
    ephemeris is held fixed (shifts are O-C relative to it; M0, P0 are degenerate with a mean
    or linear shift). The timing errors are the diagonal of the joint covariance, so they
    include the template uncertainty and its correlation with the shifts; ``cov_T`` is the
    full matrix. Points outside accepted events get zero shift. Events are those accepted by
    ``virtual_minima`` with the same criteria.

    Args:
        t (ndarray): Times (days).
        y (ndarray): Magnitudes (mag).
        sigma (ndarray): Magnitude errors (mag).
        model (EBModel): Template model.
        p0 (dict): Starting template parameters.
        free (list of str): Template parameters to fit.
        ephem (dict): Fixed ephemeris (M0, P0 in days; optional Pdot).
        seg (ndarray or None): Baseline segment id per point.
        order (int): Baseline polynomial degree.
        window (float): Half-width of the phase window (cycles, capped at 0.25).
        min_points (int): Minimum points per event.
        min_each_side (int): Minimum points on each side of the centre.
        wall_band (tuple or dict or None): Wall coverage requirement, see ``virtual_minima``.
        min_wall (int): Minimum points per wall.
        rescale (bool): Rescale the covariance by the reduced chi-square.

    Returns:
        ShiftFit: Template, shifts and joint covariance.

    Raises:
        FitDidNotConverge: If the solver does not converge.
    """
    t, y, sigma = (np.asarray(a, float) for a in (t, y, sigma))
    w = 1.0 / sigma
    M0, P0, Pd = ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0)
    th_c = phase(t, M0, P0, Pd)
    En = np.round(2 * th_c) / 2
    window = min(window, 0.25)
    vm = virtual_minima(
        t,
        y,
        sigma,
        model,
        p0,
        ephem,
        window=window,
        order=order,
        min_points=min_points,
        min_each_side=min_each_side,
        wall_band=wall_band,
        min_wall=min_wall,
    )
    ev = vm["E"]
    K = len(ev)
    pos = {float(e): i for i, e in enumerate(ev)}
    idx = np.array(
        [
            pos.get(float(e), -1) if abs(th - e) < window else -1
            for th, e in zip(th_c, En)
        ]
    )
    asg = idx >= 0
    Pl = P0 * (1.0 + Pd * En)
    nui = _Nuisance(t, w, seg, order)
    free = list(free)
    nf = len(free)

    def resid(x):
        """Weighted, baseline-projected residuals with per-event shifts applied.

        Args:
            x (ndarray): Template parameters followed by the shifts (days).

        Returns:
            ndarray: Residuals (dimensionless).
        """
        p = dict(p0)
        for n, v in zip(free, x[:nf]):
            p[n] = v
        th = th_c.copy()
        th[asg] -= x[nf:][idx[asg]] / Pl[asg]
        return nui.project(w * (y - model(th, p)))

    x0 = np.r_[[p0[n] for n in free], vm["delta"]]
    sol = least_squares(resid, x0, method="lm", x_scale="jac")
    _check_solution(sol, "fit_global_shifts")
    dof = len(t) - nui.npar - len(x0)
    cov, chi2, rchi2 = _cov(sol, dof, rescale)
    p = dict(p0)
    for n, v in zip(free, sol.x[:nf]):
        p[n] = v
    d = sol.x[nf:]
    cT = cov[nf:, nf:]
    T = np.array([predicted_time(e, M0, P0, Pd) for e in ev]) + d
    return ShiftFit(
        p,
        dict(zip(free, np.sqrt(np.diag(cov)[:nf]))),
        ev,
        T,
        np.sqrt(np.diag(cT)),
        d,
        cT,
        vm["sigma"],
        list(free) + [f"d[{e}]" for e in ev],
        cov,
        chi2,
        dof,
        rchi2,
        float(np.sqrt(np.mean((sol.fun / w) ** 2))),
        int((~asg).sum()),
    )
