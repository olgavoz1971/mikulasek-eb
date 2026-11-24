"""Parametrisation-free descriptors of the fitted eclipses: depth and full width at half depth.

The parameters D, Gamma, C and A of an eclipse profile are correlated and have a different
meaning in each Table 2 profile, so they cannot be compared between fits that use different
profiles, and even within one profile D, Gamma and C trade against each other. The two
numbers returned here are properties of the fitted eclipse curve itself, so two fits can be
compared through them whatever profile, harmonic set or baseline each used.
"""

from __future__ import annotations

import numpy as np

from .model import EBModel
from .solution import ModelSolution

# Grid over one eclipse, in cycles around its centre. The step is 2.5e-5 cycles.
GRID_HALF_WIDTH = 0.25
GRID_N = 20001


def _eclipse_function(model: EBModel, p: dict, k: int, phi: np.ndarray) -> np.ndarray:
    """Evaluate the eclipse function of eclipse ``k`` alone.

    The eclipse function is the sum of the profile copies of that eclipse (the main
    profile and, if present, its second component). Harmonics, the O'Connell term, the
    other eclipse and the exposure integration are all left out.

    Args:
        model (EBModel): Model that defines the profile and the K term.
        p (dict): Template parameters.
        k (int): Eclipse number, 1 (primary, at phase 0) or 2 (secondary, at phase 0.5).
        phi (ndarray): Phase offsets from the eclipse centre (cycles).

    Returns:
        ndarray: Eclipse function (mag, positive means fainter).
    """
    bare = EBModel(model.profile, harmonics=(), ocon=False, kterm=model.kterm)
    other = 3 - k
    q = dict(p)
    q[f"A{other}"] = 0.0
    if f"A{other}b" in q:
        q[f"A{other}b"] = 0.0
    return bare(0.5 * (k - 1) + phi, q)


def _half_depth_width(phi: np.ndarray, values: np.ndarray, k: int) -> float:
    """Full width at half depth of a sampled eclipse function.

    The width is the distance between the two half-depth crossings on either side of the
    deepest point, with linear interpolation between grid points.

    Args:
        phi (ndarray): Phase offsets (cycles), increasing.
        values (ndarray): Eclipse function (mag) on ``phi``.
        k (int): Eclipse number, used only in error messages.

    Returns:
        float: Full width at half depth (cycles).

    Raises:
        ValueError: If the curve does not fall below half depth on both sides of the
            deepest point within the grid, or has no positive depth.
    """
    peak_index = int(np.argmax(values))
    peak = float(values[peak_index])
    if not np.isfinite(peak) or peak <= 0.0:
        raise ValueError(f"Eclipse {k} has no positive depth (peak {peak!r} mag).")
    half = 0.5 * peak
    below = np.flatnonzero(values < half)
    left = below[below < peak_index]
    right = below[below > peak_index]
    if left.size == 0 or right.size == 0:
        raise ValueError(
            f"Eclipse {k} does not fall below half depth within "
            f"+-{GRID_HALF_WIDTH} cycles of its centre; no full width can be measured."
        )
    i0, i1 = int(left[-1]), int(right[0])
    x_left = phi[i0] + (half - values[i0]) * (phi[i0 + 1] - phi[i0]) / (
        values[i0 + 1] - values[i0]
    )
    x_right = phi[i1 - 1] + (half - values[i1 - 1]) * (phi[i1] - phi[i1 - 1]) / (
        values[i1] - values[i1 - 1]
    )
    return float(x_right - x_left)


def eclipse_descriptors(source: EBModel | ModelSolution, p: dict | None = None) -> dict:
    """Depth and full width at half depth of both fitted eclipses.

    The numbers describe the fitted eclipse curve, not its parameters, so they can be
    compared between fits with different profiles, harmonic sets or baselines. Each
    eclipse is taken as the sum of its profile copies (the main profile and, when present,
    the second component of the two-component extension). Proximity terms (harmonics,
    O'Connell), the other eclipse, the baseline and the exposure integration are excluded:
    the descriptors refer to the instantaneous eclipse, measured from its own zero level
    in the far wings. A fit that moves part of the eclipse into low harmonics changes these
    numbers accordingly, and that change is real information about the curve.

    Only the point values are returned. A solution stores parameter errors but not their
    covariance, so no honest uncertainty can be propagated here.

    Args:
        source (EBModel | ModelSolution): A fitted ``ModelSolution`` (then ``p`` is
            ignored and its own parameters are used), or an ``EBModel`` together with
            ``p``.
        p (dict | None): Template parameters (A in mag, D in cycles) when ``source`` is
            an ``EBModel``; required in that case.

    Returns:
        dict: ``{"eclipse_1": {"depth_mag", "fwhm_cycles"}, "eclipse_2": {...},
        "depth_ratio", "fwhm_ratio"}``. Depth is the deepest point of the eclipse
        function (mag, fainter is positive), the width is the full width at half of that
        depth (cycles). The ratios are eclipse 2 divided by eclipse 1.

    Raises:
        ValueError: If ``p`` is missing for an ``EBModel``, or an eclipse has no positive
            depth or does not fall below half depth within +-0.25 cycles.
    """
    if isinstance(source, ModelSolution):
        model, params = source.build_model(), source.params
    else:
        if p is None:
            raise ValueError("p is required when the source is an EBModel.")
        model, params = source, p
    phi = np.linspace(-GRID_HALF_WIDTH, GRID_HALF_WIDTH, GRID_N)
    out = {}
    for k in (1, 2):
        values = _eclipse_function(model, params, k, phi)
        out[f"eclipse_{k}"] = dict(
            depth_mag=float(np.max(values)),
            fwhm_cycles=_half_depth_width(phi, values, k),
        )
    out["depth_ratio"] = out["eclipse_2"]["depth_mag"] / out["eclipse_1"]["depth_mag"]
    out["fwhm_ratio"] = out["eclipse_2"]["fwhm_cycles"] / out["eclipse_1"]["fwhm_cycles"]
    return out
