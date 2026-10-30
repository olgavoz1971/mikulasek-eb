"""Global least-squares fits with per-segment baselines, uncertainty propagation and correlations (Sects. 5.1 and 5.2)."""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from .model import phase

logger = logging.getLogger(__name__)


# ------------------------------------------------- linear nuisance (Eq. 11)
class _Nuisance:
    """Per-segment polynomial baselines eliminated by weighted linear least squares.

    This is variable projection: the baselines never enter the nonlinear parameter vector.
    """

    def __init__(self, t, w, seg, order):
        """Build the weighted Vandermonde blocks.

        Args:
            t (ndarray): Times (days).
            w (ndarray): Weights, 1/sigma (1/mag).
            seg (ndarray or None): Segment id per point; None means one segment.
            order (int): Polynomial degree per segment (reduced if a segment is too short).
        """
        seg = np.zeros(len(t), int) if seg is None else np.asarray(seg)
        self.blocks, self.npar = [], 0
        for s in np.unique(seg):
            idx = np.where(seg == s)[0]
            tt = t[idx]
            sd = tt.std()
            x = (tt - tt.mean()) / (sd if sd > 0 else 1.0)
            V = (
                np.vander(x, min(order, len(idx) - 1) + 1, increasing=True)
                * w[idx, None]
            )
            self.blocks.append((idx, V))
            self.npar += V.shape[1]

    def share(self, mask):
        """Number of baseline parameters attributable to a subset of points.

        Args:
            mask (ndarray): Boolean selection of points.

        Returns:
            float: Fractional parameter count.
        """
        return sum(V.shape[1] * mask[idx].mean() for idx, V in self.blocks)

    def project(self, rw):
        """Remove the best-fit baselines from weighted residuals.

        Args:
            rw (ndarray): Weighted residuals (dimensionless).

        Returns:
            ndarray: Weighted residuals orthogonal to the baselines.
        """
        out = rw.copy()
        for idx, V in self.blocks:
            c = np.linalg.lstsq(V, rw[idx], rcond=None)[0]
            out[idx] = rw[idx] - V @ c
        return out


# ------------------------------------------------------------------ results
@dataclass
class FitResult:
    """Outcome of ``fit_global``.

    Attributes:
        p (dict): Fitted template parameters (mag, cycles).
        ephem (dict): Fitted ephemeris (M0, P0 in days, Pdot dimensionless).
        err (dict): 1-sigma errors of every fitted parameter.
        chi2, dof, rchi2 (float, int, float): Chi-square, degrees of freedom, reduced chi-square.
        rms (float): Rms residual (mag).
        names (list): Fitted parameter names in covariance order.
        cov (ndarray): Covariance matrix (rescaled by rchi2 if requested).
        success (bool): Solver success flag.
        resid (ndarray): Residuals (mag) after removing the baselines.
        nuis: Baseline bookkeeping used by ``fit_global_iterative``.
        rank_deficient (bool): True if the Jacobian is rank deficient (degenerate parameters).
    """

    p: dict
    ephem: dict
    err: dict
    chi2: float
    dof: int
    rchi2: float
    rms: float
    names: list
    cov: np.ndarray
    success: bool
    resid: np.ndarray = None  # residuals [mag] after removing the baselines
    nuis: object = None  # baseline bookkeeping (used by fit_global_iterative)
    rank_deficient: bool = (
        False  # True if the Jacobian is rank deficient (degenerate parameters)
    )


class FitDidNotConverge(RuntimeError):
    """Raised when the Levenberg-Marquardt solver does not report convergence."""


def _check_solution(sol, what):
    """Fail fast on a non-converged solver result; report a rank-deficient Jacobian.

    Args:
        sol: ``scipy.optimize.OptimizeResult`` from ``least_squares``.
        what (str): Short description of the fit for the messages.

    Returns:
        bool: True if the Jacobian is rank deficient (some parameters are degenerate, so
        the pseudo-inverse covariance is not a proper covariance).

    Raises:
        FitDidNotConverge: If the solver reports no convergence.
    """
    if not sol.success:
        raise FitDidNotConverge(
            f"{what}: solver did not converge (status={sol.status}: {sol.message})"
        )
    J = sol.jac
    rank = int(np.linalg.matrix_rank(J))
    deficient = rank < J.shape[1]
    if deficient:
        logger.warning(
            "%s: Jacobian rank %d < %d parameters; covariance from pseudo-inverse "
            "(degenerate parameters, errors unreliable)",
            what,
            rank,
            J.shape[1],
        )
    return deficient


def _cov(sol, dof, rescale):
    """Covariance matrix from a least-squares solution.

    Args:
        sol: ``scipy.optimize.OptimizeResult`` with ``jac`` and ``cost``.
        dof (int): Degrees of freedom.
        rescale (bool): Multiply the covariance by the reduced chi-square.

    Returns:
        tuple: (cov, chi2, rchi2).
    """
    J = sol.jac
    cov = np.linalg.pinv(J.T @ J)
    chi2 = 2.0 * sol.cost
    rchi2 = chi2 / max(dof, 1)
    return (cov * rchi2 if rescale else cov), chi2, rchi2


# ------------------------------------------------------------- global fit
def fit_global(
    t,
    y,
    sigma,
    model,
    p0,
    free,
    ephem0,
    free_ephem=("M0", "P0"),
    seg=None,
    order=0,
    rescale=True,
):
    """Simultaneous fit of template parameters and ephemeris to all data.

    Segments ``seg`` (for example night or orbit id) each get a polynomial baseline of
    degree ``order``. With ``rescale=True`` the covariance is multiplied by the reduced
    chi-square (paper, Sect. 5.2: weights from residual scatter).

    Args:
        t (ndarray): Times (days).
        y (ndarray): Magnitudes (mag).
        sigma (ndarray): Magnitude errors (mag).
        model (EBModel): Model to fit.
        p0 (dict): Starting template parameters.
        free (list of str): Names of the template parameters to fit.
        ephem0 (dict): Starting ephemeris (M0, P0 in days; optional Pdot).
        free_ephem (tuple of str): Subset of M0, P0, Pdot to fit.
        seg (ndarray or None): Segment id per point.
        order (int): Baseline polynomial degree per segment.
        rescale (bool): Rescale the covariance by the reduced chi-square.

    Returns:
        FitResult: Fitted parameters, errors and diagnostics.

    Raises:
        FitDidNotConverge: If the solver does not converge.
    """
    t, y, sigma = (np.asarray(a, float) for a in (t, y, sigma))
    w = 1.0 / sigma
    nui = _Nuisance(t, w, seg, order)
    free, fe = list(free), list(free_ephem)
    names = free + fe

    def unpack(x):
        """Split the solver vector into template and ephemeris dicts.

        Args:
            x (ndarray): Parameter vector in the order ``free + free_ephem``.

        Returns:
            tuple: (template dict, ephemeris dict).
        """
        p, e = dict(p0), dict(ephem0)
        for n, v in zip(names, x):
            (p if n in free else e)[n] = v
        return p, e

    def resid(x):
        """Weighted, baseline-projected residuals for the solver.

        Args:
            x (ndarray): Parameter vector.

        Returns:
            ndarray: Residuals (dimensionless).
        """
        p, e = unpack(x)
        F = model(phase(t, e["M0"], e["P0"], e.get("Pdot", 0.0)), p)
        return nui.project(w * (y - F))

    x0 = [p0[n] for n in free] + [ephem0[n] for n in fe]
    sol = least_squares(resid, x0, method="lm", x_scale="jac")
    deficient = _check_solution(sol, "fit_global")
    dof = len(t) - nui.npar - len(x0)
    cov, chi2, rchi2 = _cov(sol, dof, rescale)
    p, e = unpack(sol.x)
    err = dict(zip(names, np.sqrt(np.diag(cov))))
    rms = float(np.sqrt(np.mean((sol.fun / w) ** 2)))
    return FitResult(
        p,
        e,
        err,
        chi2,
        dof,
        rchi2,
        rms,
        names,
        cov,
        sol.success,
        resid=sol.fun / w,
        nuis=nui,
        rank_deficient=deficient,
    )


# =====================================================================================
# Section 5.1 of the paper: uncertainties, correlations, simultaneous solution
# =====================================================================================
def correlation_matrix(r):
    """Correlation matrix of the fitted parameters.

    Args:
        r (FitResult): Fit result.

    Returns:
        tuple: (names, correlation matrix).
    """
    sd = np.sqrt(np.diag(r.cov))
    return list(r.names), r.cov / np.outer(sd, sd)


def strong_correlations(r, thr=0.9):
    """List parameter pairs whose correlation is strong (degeneracy warnings).

    Args:
        r (FitResult): Fit result.
        thr (float): Absolute correlation threshold.

    Returns:
        list: Tuples (name_i, name_j, correlation) with |correlation| >= thr.
    """
    names, c = correlation_matrix(r)
    return [
        (names[i], names[j], float(c[i, j]))
        for i in range(len(names))
        for j in range(i)
        if abs(c[i, j]) >= thr
    ]


def curve_band(model, r, theta, rel_step=1e-5):
    """Template and its 1-sigma band from the full parameter covariance.

    Uses the general law of uncertainty propagation (paper Sect. 5.1).

    Args:
        model (EBModel): Template model.
        r (FitResult): Fit result with covariance.
        theta (ndarray): Phases (cycles).
        rel_step (float): Relative step for the numerical derivatives.

    Returns:
        tuple: (F(theta) in mag, 1-sigma uncertainty in mag).
    """
    theta = np.asarray(theta, float)
    pn = [n for n in r.names if n in r.p]
    ii = [r.names.index(n) for n in pn]
    C = r.cov[np.ix_(ii, ii)]
    J = np.empty((len(theta), len(pn)))
    for k, n in enumerate(pn):
        h = rel_step * max(abs(r.p[n]), 1e-3)
        pp, pm = dict(r.p), dict(r.p)
        pp[n] += h
        pm[n] -= h
        J[:, k] = (model(theta, pp) - model(theta, pm)) / (2 * h)
    var = np.einsum("ij,jk,ik->i", J, C, J)
    return model(theta, r.p), np.sqrt(np.maximum(var, 0.0))


def ephemeris_at_epoch(r, t=None, Ec=None):
    """Re-express the fitted ephemeris at an epoch where M0' and P0' are uncorrelated.

    T(E) = M0' + P0' (E - Ec) + 0.5 P0' Pdot (E - Ec)^2 (paper Sect. 5.1).

    Args:
        r (FitResult): Fit result with covariance.
        t (ndarray or None): Observation times (days) defining the search range for Ec.
        Ec (float or None): Epoch to use; if None, the integer epoch of minimal
            correlation within the data range.

    Returns:
        dict: Ec, M0, M0_err, P0, P0_err (days), Pdot and the correlation before and after.
    """
    eph = r.ephem
    M0, P0, Pd = eph["M0"], eph["P0"], eph.get("Pdot", 0.0)
    fe = [n for n in ("M0", "P0", "Pdot") if n in r.names]
    ii = [r.names.index(n) for n in fe]
    C = r.cov[np.ix_(ii, ii)]

    def stats(E):
        """Covariance and correlation of (M0', P0') at epoch E.

        Args:
            E (float): Epoch.

        Returns:
            tuple: (2x2 covariance, correlation coefficient).
        """
        dM = {"M0": 1.0, "P0": E + 0.5 * Pd * E**2, "Pdot": 0.5 * P0 * E**2}
        dP = {"M0": 0.0, "P0": 1.0 + Pd * E, "Pdot": P0 * E}
        J = np.array([[dM[n] for n in fe], [dP[n] for n in fe]])
        S = J @ C @ J.T
        return S, S[0, 1] / np.sqrt(S[0, 0] * S[1, 1])

    if Ec is None:
        if t is not None:
            E = (np.asarray(t, float) - M0) / P0
            grid = np.arange(np.floor(E.min()), np.ceil(E.max()) + 1)
        else:
            grid = np.arange(-1000.0, 1001.0)
        Ec = float(grid[int(np.argmin([abs(stats(e)[1]) for e in grid]))])
    S, corr = stats(Ec)
    return dict(
        Ec=Ec,
        M0=M0 + P0 * Ec + 0.5 * P0 * Pd * Ec**2,
        M0_err=float(np.sqrt(S[0, 0])),
        P0=P0 * (1 + Pd * Ec),
        P0_err=float(np.sqrt(S[1, 1])),
        Pdot=Pd,
        corr_M0_P0_before=float(stats(0.0)[1]),
        corr_M0_P0_after=float(corr),
    )


def fit_global_iterative(
    t,
    y,
    sigma,
    model,
    p0,
    free,
    ephem0,
    free_ephem=("M0", "P0"),
    seg=None,
    order=0,
    subsets=None,
    n_iter=8,
    tol=0.01,
    rescale=False,
):
    """Iteratively estimate unknown uncertainties from subset residuals (paper Sect. 5.2, principle 2).

    Per subset the sigmas are multiplied by sqrt(chi2_s / dof_s) until all factors are within
    ``tol`` of 1. A systematic model misfit concentrated in part of a subset is absorbed into
    that subset's sigma.

    Args:
        t (ndarray): Times (days).
        y (ndarray): Magnitudes (mag).
        sigma (ndarray): Starting magnitude errors (mag).
        model (EBModel): Model to fit.
        p0 (dict): Starting template parameters.
        free (list of str): Template parameters to fit.
        ephem0 (dict): Starting ephemeris.
        free_ephem (tuple of str): Ephemeris parameters to fit.
        seg (ndarray or None): Baseline segment id per point.
        order (int): Baseline polynomial degree.
        subsets (ndarray or None): Subset id per point; defaults to ``seg``.
        n_iter (int): Maximum iterations.
        tol (float): Convergence tolerance on the scale factors.
        rescale (bool): Rescale the covariance in the fits.

    Returns:
        tuple: (FitResult, rescaled sigma array in mag, history of factors per iteration).
    """
    t, y = np.asarray(t, float), np.asarray(y, float)
    sig = np.asarray(sigma, float).copy()
    sub = (
        np.zeros(len(t), int)
        if (subsets is None and seg is None)
        else np.asarray(seg if subsets is None else subsets)
    )
    p, e, hist = dict(p0), dict(ephem0), []
    for it in range(n_iter):
        r = fit_global(
            t,
            y,
            sig,
            model,
            p,
            free,
            e,
            free_ephem,
            seg=seg,
            order=order,
            rescale=rescale,
        )
        chi = (r.resid / sig) ** 2
        fac = {}
        for u in np.unique(sub):
            m = sub == u
            dof = max(m.sum() - r.nuis.share(m), 1.0)
            fac[u] = float(np.sqrt(chi[m].sum() / dof))
        f_pt = np.array([fac[u] for u in sub])
        hist.append(fac)
        sig = sig * f_pt
        p, e = r.p, r.ephem
        if max(abs(f - 1) for f in fac.values()) < tol:
            break
    r = fit_global(
        t, y, sig, model, p, free, e, free_ephem, seg=seg, order=order, rescale=rescale
    )
    return r, sig, hist
