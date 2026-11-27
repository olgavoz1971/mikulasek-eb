"""Parametrisation-free descriptors of the fitted eclipses: depth and full width at half depth.

The parameters D, Gamma, C and A of an eclipse profile are correlated and have a different
meaning in each Table 2 profile, so they cannot be compared between fits that use different
profiles, and even within one profile D, Gamma and C trade against each other. The two
numbers returned here are properties of the fitted eclipse curve itself, so two fits can be
compared through them whatever profile, harmonic set or baseline each used.
"""

from __future__ import annotations

import numpy as np

from .config import HALF_ORBIT_WINDOW
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


def _level_width(phi: np.ndarray, values: np.ndarray, k: int, level: float) -> float:
    """Full width of a sampled eclipse function at a fraction of its depth.

    The width is the distance between the two crossings of ``level * depth`` on either
    side of the deepest point, with linear interpolation between grid points.

    Args:
        phi (ndarray): Phase offsets (cycles), increasing.
        values (ndarray): Eclipse function (mag) on ``phi``.
        k (int): Eclipse number, used only in error messages.
        level (float): Fraction of the depth, strictly between 0 and 1.

    Returns:
        float: Full width at that level (cycles).

    Raises:
        ValueError: If the curve does not fall below the level on both sides of the
            deepest point within the grid, or has no positive depth.
    """
    peak_index = int(np.argmax(values))
    peak = float(values[peak_index])
    if not np.isfinite(peak) or peak <= 0.0:
        raise ValueError(f"Eclipse {k} has no positive depth (peak {peak!r} mag).")
    target = level * peak
    below = np.flatnonzero(values < target)
    left = below[below < peak_index]
    right = below[below > peak_index]
    if left.size == 0 or right.size == 0:
        label = "half depth" if level == 0.5 else f"{100 * level:g} % of its depth"
        raise ValueError(
            f"Eclipse {k} does not fall below {label} within "
            f"+-{GRID_HALF_WIDTH} cycles of its centre; no full width can be measured."
        )
    i0, i1 = int(left[-1]), int(right[0])
    x_left = phi[i0] + (target - values[i0]) * (phi[i0 + 1] - phi[i0]) / (
        values[i0 + 1] - values[i0]
    )
    x_right = phi[i1 - 1] + (target - values[i1 - 1]) * (phi[i1] - phi[i1 - 1]) / (
        values[i1] - values[i1 - 1]
    )
    return float(x_right - x_left)


def _half_depth_width(phi: np.ndarray, values: np.ndarray, k: int) -> float:
    """Full width at half depth of a sampled eclipse function.

    Args:
        phi (ndarray): Phase offsets (cycles), increasing.
        values (ndarray): Eclipse function (mag) on ``phi``.
        k (int): Eclipse number, used only in error messages.

    Returns:
        float: Full width at half depth (cycles).

    Raises:
        ValueError: As ``_level_width`` with level 0.5.
    """
    return _level_width(phi, values, k, 0.5)


def _unpack_source(source: EBModel | ModelSolution, p: dict | None) -> tuple[EBModel, dict]:
    """Return the model and the parameters of a solution or of a model with ``p``.

    Args:
        source (EBModel | ModelSolution): A fitted solution or a model.
        p (dict | None): Parameters, required when ``source`` is a model.

    Returns:
        tuple[EBModel, dict]: Model and template parameters.

    Raises:
        ValueError: If ``source`` is a model and ``p`` is None.
    """
    if isinstance(source, ModelSolution):
        return source.build_model(), source.params
    if p is None:
        raise ValueError("p is required when the source is an EBModel.")
    return source, p


def eclipse_widths(
    source: EBModel | ModelSolution, p: dict | None = None, level: float = 0.05
) -> dict:
    """Full width of both fitted eclipses at a fraction of their depth.

    Same eclipse function as ``eclipse_descriptors`` (main profile plus second component,
    no harmonics, no O'Connell term, no other eclipse, no exposure integration), so the
    result does not depend on the harmonics fitted to the curve and does not need contacts
    on the light curve. The level 0.5 gives the full width at half depth. A level of 0 is
    not allowed: profiles such as profile 1 only tend to zero.

    Args:
        source (EBModel | ModelSolution): A fitted solution (``p`` is then ignored) or a
            model together with ``p``.
        p (dict | None): Template parameters when ``source`` is a model.
        level (float): Fraction of the depth, strictly between 0 and 1. Default 0.05.

    Returns:
        dict: ``{"primary": width, "secondary": width}`` in cycles.

    Raises:
        ValueError: If ``level`` is not in (0, 1), ``p`` is missing, an eclipse has no
            positive depth, or an eclipse does not fall below the level within
            +-0.25 cycles.
    """
    if not 0.0 < level < 1.0:
        raise ValueError(f"level must be strictly between 0 and 1, got {level!r}.")
    model, params = _unpack_source(source, p)
    phi = np.linspace(-GRID_HALF_WIDTH, GRID_HALF_WIDTH, GRID_N)
    return {
        kind: _level_width(phi, _eclipse_function(model, params, k, phi), k, level)
        for k, kind in ((1, "primary"), (2, "secondary"))
    }


def window_half_widths(
    source: EBModel | ModelSolution,
    p: dict | None = None,
    n_widths: float = 1.5,
    level: float = 0.05,
    max_half: float = HALF_ORBIT_WINDOW,
) -> dict:
    """Half-windows for timing, set from the width of the fitted eclipse.

    The full window of each eclipse is ``n_widths`` times its full width at ``level`` of
    its depth (see ``eclipse_widths``), so the half-window is half of that. The default,
    1.5 widths at 5 percent, is a preliminary choice tested on two synthetic stars, not a
    rule. Nothing is capped or enlarged silently: a window above ``max_half`` is an error.

    Args:
        source (EBModel | ModelSolution): A fitted solution or a model with ``p``.
        p (dict | None): Template parameters when ``source`` is a model.
        n_widths (float): Full window in units of the width at ``level``. Must be
            positive.
        level (float): Fraction of the depth that defines the width, in (0, 1).
        max_half (float): Largest half-window allowed, in cycles.

    Returns:
        dict: ``{"primary": half, "secondary": half}`` in cycles, to be passed as
        ``TimingRecipe(window=half)`` for each kind.

    Raises:
        ValueError: If ``n_widths`` is not positive, a width cannot be measured, or a
            half-window exceeds ``max_half``.
    """
    if not n_widths > 0.0:
        raise ValueError(f"n_widths must be positive, got {n_widths!r}.")
    out = {}
    for kind, width in eclipse_widths(source, p, level).items():
        half = 0.5 * n_widths * width
        if half > max_half:
            raise ValueError(
                f"The {kind} half-window {half:.4f} cycles ({n_widths:g} widths of "
                f"{width:.4f} at {100 * level:g} %) exceeds max_half {max_half:g}. "
                "Choose a smaller n_widths or a larger max_half explicitly."
            )
        out[kind] = half
    return out


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
    model, params = _unpack_source(source, p)
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
