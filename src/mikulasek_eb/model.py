"""Eclipse profiles, the EBModel light-curve model and the quadratic phase function (Mikulasek 2015, Eqs. 13, 15, 17)."""

from __future__ import annotations


import numpy as np


# ------------------------------------------------------------------ profiles
def _one_minus_pow(u, G):
    """Compute 1 - (1 - u)^G without cancellation for tiny u.

    Args:
        u (ndarray): Values in [0, 1].
        G (float): Exponent (Gamma, dimensionless).

    Returns:
        ndarray: 1 - (1 - u)^G.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return -np.expm1(G * np.log1p(-u))


def prof1(phi, D, G, C=0.0):  # Eq. (13) / Table 2 #1
    """Normalised eclipse profile: Eq. (13) / Table 2 profile 1 with curvature term C.

    Args:
        phi (ndarray): Phase offset from the eclipse centre (cycles).
        D (float): Width parameter (cycles).
        G (float): Kurtosis Gamma (dimensionless).
        C (float): Curvature of the bottom (dimensionless).

    Returns:
        ndarray: Profile, 1 at the centre and 0 far from it (multiply by A in mag).
    """
    x = np.clip(phi / D, -60.0, 60.0)
    return (1.0 + C * x**2) * _one_minus_pow(np.exp(1.0 - np.cosh(x)), G)


def prof2(phi, D, G, C=0.0):  # Table 2 #2 (no C)
    """Normalised eclipse profile: Table 2 profile 2 (the curvature C is ignored).

    Args:
        phi (ndarray): Phase offset from the eclipse centre (cycles).
        D (float): Width parameter (cycles).
        G (float): Kurtosis Gamma (dimensionless).
        C (float): Curvature of the bottom (dimensionless).

    Returns:
        ndarray: Profile, 1 at the centre and 0 far from it (multiply by A in mag).
    """
    return prof1(phi, D, G, 0.0)


def prof3(phi, D, G, C=0.0):  # Table 2 #3
    """Normalised eclipse profile: Table 2 profile 3 (Gaussian-based; the curvature C is ignored).

    Args:
        phi (ndarray): Phase offset from the eclipse centre (cycles).
        D (float): Width parameter (cycles).
        G (float): Kurtosis Gamma (dimensionless).
        C (float): Curvature of the bottom (dimensionless).

    Returns:
        ndarray: Profile, 1 at the centre and 0 far from it (multiply by A in mag).
    """
    x = np.clip(phi / D, -60.0, 60.0)
    return _one_minus_pow(np.exp(-0.5 * x**2), G)


def prof4(phi, D, G, C=0.0):  # Table 2 #4: Re[(1-|x|^G)^(3/2)]
    """Normalised eclipse profile: Table 2 profile 4, Re[(1 - |x|^G)^(3/2)] (Andronov type; C is ignored).

    Args:
        phi (ndarray): Phase offset from the eclipse centre (cycles).
        D (float): Width parameter (cycles).
        G (float): Kurtosis Gamma (dimensionless).
        C (float): Curvature of the bottom (dimensionless).

    Returns:
        ndarray: Profile, 1 at the centre and 0 far from it (multiply by A in mag).
    """
    base = 1.0 - np.abs(phi / D) ** G
    return np.where(base > 0.0, np.abs(base) ** 1.5, 0.0)


PROFILES = {1: prof1, 2: prof2, 3: prof3, 4: prof4}


# ------------------------------------------------------------------ model
class EBModel:
    """Phenomenological eclipsing-binary light curve in magnitudes (Mikulasek 2015).

    F(theta) = sum_k A_k prof(phi_k) + sum H_k cos(2 pi k theta) [+ S1 sin(2 pi theta)].

    Parameter dict keys: A1, A2 (mag), D1, G1, C1 (D2, G2, C2 optional; if absent they
    equal D1, G1, C1), H<k> (mag) for k in ``harmonics``, S1 (mag) if ``ocon``.

    Optional second component per eclipse k (EXTENSION, not in the paper): if the keys
    A<k>b, D<k>b, G<k>b are present, A<k>b * profile(phi; D<k>b, G<k>b, C=0) is added
    to eclipse k. This lets a flat bottom and steep walls coexist; amplitudes may be negative.

    Attributes:
        profile (int): Table 2 profile number (1 to 4).
        harmonics (tuple): Orders k of the cosine terms.
        ocon (bool): Whether the O'Connell sine term S1 is included.
    """

    def __init__(self, profile=1, harmonics=(2,), ocon=False):
        """Store the model configuration.

        Args:
            profile (int): Table 2 profile number (1 to 4).
            harmonics (tuple of int): Orders of the cosine terms H_k.
            ocon (bool): Include the O'Connell sine term S1.
        """
        self.profile, self.harmonics, self.ocon = profile, tuple(harmonics), ocon

    def names(self, shared_shape=True):
        """Names of the template parameters that can be fitted.

        Args:
            shared_shape (bool): If True the two eclipses share D, G, C; otherwise D2, G2, C2
                are listed as well.

        Returns:
            list of str: Parameter names.
        """
        n = ["A1", "A2", "D1", "G1", "C1"]
        if not shared_shape:
            n += ["D2", "G2", "C2"]
        n += [f"H{k}" for k in self.harmonics]
        return n + (["S1"] if self.ocon else [])

    def __call__(self, theta, p):
        """Evaluate the model.

        Args:
            theta (ndarray): Phase (cycles), primary at integer values.
            p (dict): Parameters (A in mag, D in cycles, see class docstring).

        Returns:
            ndarray: Magnitude relative to the out-of-eclipse level.
        """
        theta = np.asarray(theta, float)
        prof = PROFILES[self.profile]
        F = np.zeros_like(theta)
        for k in (1, 2):
            c = 0.5 * (k - 1)
            phi = (theta - c) - np.round(theta - c)
            F += p[f"A{k}"] * prof(
                phi,
                p.get(f"D{k}", p["D1"]),
                p.get(f"G{k}", p["G1"]),
                p.get(f"C{k}", p["C1"]),
            )
            if f"A{k}b" in p:
                F += p[f"A{k}b"] * prof(phi, p[f"D{k}b"], p[f"G{k}b"], 0.0)
        for k in self.harmonics:
            F += p[f"H{k}"] * np.cos(2 * np.pi * k * theta)
        if self.ocon:
            F += p["S1"] * np.sin(2 * np.pi * theta)
        return F


def phase(t, M0, P0, Pdot=0.0):
    """Quadratic phase function theta(t), Eq. (17).

    Args:
        t (ndarray): Times (days; MJD in the pipeline).
        M0 (float): Reference primary minimum (same time system, days).
        P0 (float): Period at M0 (days).
        Pdot (float): Dimensionless period derivative dP/dt.

    Returns:
        ndarray: Phase in cycles.
    """
    v = (np.asarray(t, float) - M0) / P0
    return v - 0.5 * Pdot * v**2


def predicted_time(E, M0, P0, Pdot=0.0):
    """Predicted time of eclipse number E from the quadratic ephemeris.

    Args:
        E (float or ndarray): Eclipse number (integer primary, n + 0.5 secondary).
        M0 (float): Reference primary minimum (days).
        P0 (float): Period at M0 (days).
        Pdot (float): Dimensionless period derivative.

    Returns:
        float or ndarray: Time of minimum (days, same system as M0).
    """
    return M0 + P0 * E + 0.5 * P0 * Pdot * E**2
