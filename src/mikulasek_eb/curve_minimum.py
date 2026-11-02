"""Times of one minimum read off a light curve, without a fitted model.

Several definitions of "the" minimum, for a pure curve or a smooth estimate of one.
The curve is used as given: nothing is fitted and no point is rejected.

* ``lowest`` is the time of the extreme sample.
* A bisector at a fraction ``f`` (0 at the extreme, 1 at the continuum) is the midpoint of
  the two times where the curve crosses that level, one on each wing.
* The centroid is the first moment of the light missing below the continuum (flux) or of
  the light in excess of it (magnitude).

The same fractions as the generator's ``T_min_pure`` (0.20, 0.35, 0.50, 0.65, 0.80) are the
default, so the median of the bisectors is that definition when the continuum is the 90th
percentile. Each fraction is kept, because their spread measures the asymmetry.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from .errors import InvalidLightCurve

logger = logging.getLogger(__name__)

BISECTOR_FRACTIONS = (0.20, 0.35, 0.50, 0.65, 0.80)


@dataclass(frozen=True)
class CurveMinimum:
    """One minimum under several definitions.

    Attributes:
        t_lowest (float): Time of the extreme sample.
        t_bisector (dict): Fraction from the extreme to the bisector time.
        t_bisector_median (float): Median of the bisector times.
        t_centroid (float): First moment of the eclipse deficit or excess.
        continuum (float): Value used as the top of the depth, in the units of ``y``.
        extreme (float): Value of the extreme sample, in the units of ``y``.
        depth (float): Continuum minus the extreme for flux, the opposite for magnitude.
    """

    t_lowest: float
    t_bisector: dict
    t_bisector_median: float
    t_centroid: float
    continuum: float
    extreme: float
    depth: float


def measure_curve_minimum(
    t,
    y,
    t_guess: float,
    half_width: float,
    *,
    kind: str = "flux",
    fractions: tuple = BISECTOR_FRACTIONS,
    continuum: float | None = None,
    continuum_percentile: float = 90.0,
    period: float | None = None,
) -> CurveMinimum:
    """Measure one minimum inside a window around ``t_guess``.

    Args:
        t (ndarray): Times, one per point. Strictly increasing is not required.
        y (ndarray): Flux or magnitude at those times.
        t_guess (float): Time near the minimum. With ``period`` set, the window is taken
            on the phase folded about this time.
        half_width (float): Half-width of the window (same unit as ``t``).
        kind (str): ``"flux"`` (the minimum is the lowest value) or ``"mag"`` (the highest).
        fractions (tuple of float): Bisector levels, each in (0, 1): 0 at the extreme and
            1 at the continuum.
        continuum (float | None): Top of the depth. None takes the ``continuum_percentile``
            of the values inside the window.
        continuum_percentile (float): Used only when ``continuum`` is None.
        period (float | None): Fold the offsets into one cycle before cutting the window.

    Returns:
        CurveMinimum: The times and the depth they were measured on.

    Raises:
        InvalidLightCurve: If the window, the depth, a wing or a requested level is unusable.
    """
    if kind not in ("flux", "mag"):
        raise InvalidLightCurve("kind must be 'flux' or 'mag'.")
    if half_width <= 0:
        raise InvalidLightCurve("half_width must be positive.")
    if not fractions or any(not 0.0 < f < 1.0 for f in fractions):
        raise InvalidLightCurve("fractions must lie in (0, 1).")
    tt = np.asarray(t, float)
    yy = np.asarray(y, float)
    if tt.shape != yy.shape or tt.ndim != 1 or len(tt) < 5:
        raise InvalidLightCurve("t and y must be 1-d arrays of at least 5 points.")
    if not np.isfinite(tt).all() or not np.isfinite(yy).all():
        raise InvalidLightCurve("t and y must be finite.")

    offset = tt - t_guess
    if period is not None:
        if period <= 0:
            raise InvalidLightCurve("period must be positive.")
        offset = np.mod(offset + 0.5 * period, period) - 0.5 * period
    window = np.abs(offset) <= half_width
    if int(window.sum()) < 5:
        raise InvalidLightCurve("fewer than 5 points inside the window.")
    offset, yy = offset[window], yy[window]
    order = np.argsort(offset)
    offset, yy = offset[order], yy[order]

    dip = yy if kind == "flux" else -yy
    i_ext = int(np.argmin(dip))
    if np.count_nonzero(dip == dip[i_ext]) > 1:
        raise InvalidLightCurve("the extreme value is reached at more than one sample.")
    extreme = float(yy[i_ext])
    if continuum is None:
        # In the dip coordinate the continuum is the bright end, whatever ``kind`` says.
        top_dip = float(np.percentile(dip, continuum_percentile))
    else:
        top_dip = float(continuum) if kind == "flux" else -float(continuum)
    top = top_dip if kind == "flux" else -top_dip
    depth = top_dip - float(dip[i_ext])
    if not np.isfinite(depth) or depth <= 0:
        raise InvalidLightCurve("the continuum is not above the minimum.")

    def wing(side: str) -> tuple[np.ndarray, np.ndarray]:
        """Points of one wing, ordered from the extreme outward.

        Args:
            side (str): ``"left"`` (offsets up to the extreme) or ``"right"``.

        Returns:
            tuple: Offsets and dip values, the extreme first.
        """
        if side == "left":
            sl = slice(i_ext, None, -1)
        else:
            sl = slice(i_ext, None)
        return offset[sl], dip[sl]

    bisectors = {}
    for frac in fractions:
        level = float(dip[i_ext]) + frac * depth
        left = _crossing(*wing("left"), level, frac)
        right = _crossing(*wing("right"), level, frac)
        bisectors[float(frac)] = t_guess + 0.5 * (left + right)
    times = np.array(list(bisectors.values()))

    deficit = np.maximum(top_dip - dip, 0.0)
    span = np.empty_like(offset)
    span[0] = offset[1] - offset[0]
    span[-1] = offset[-1] - offset[-2]
    span[1:-1] = 0.5 * (offset[2:] - offset[:-2])
    if np.any(span <= 0):
        raise InvalidLightCurve("the window times are not strictly increasing.")
    weight = deficit * span
    if weight.sum() <= 0:
        raise InvalidLightCurve("the curve never falls below the continuum.")
    t_centroid = t_guess + float(np.sum(weight * offset) / np.sum(weight))
    logger.debug(
        "curve minimum: depth %.4g, bisector spread %.4g, lowest at offset %.4g",
        depth,
        float(times.max() - times.min()),
        offset[i_ext],
    )
    return CurveMinimum(
        t_lowest=t_guess + float(offset[i_ext]),
        t_bisector=bisectors,
        t_bisector_median=float(np.median(times)),
        t_centroid=t_centroid,
        continuum=top,
        extreme=extreme,
        depth=float(depth if kind == "flux" else depth),
    )


def _crossing(offset: np.ndarray, dip: np.ndarray, level: float, frac: float) -> float:
    """Offset where one wing crosses ``level``, walking outward from the extreme.

    Args:
        offset (ndarray): Offsets of the wing, the extreme first.
        dip (ndarray): Dip values of the wing, the extreme first.
        level (float): Dip value to cross.
        frac (float): Fraction, named in the error.

    Returns:
        float: Offset of the crossing.

    Raises:
        InvalidLightCurve: If the level is not crossed exactly once.
    """
    found = None
    for i in range(len(dip) - 1):
        inner, outer = dip[i], dip[i + 1]
        if outer == level:
            if found is not None:
                raise InvalidLightCurve(
                    f"fraction {frac} is crossed more than once on one wing."
                )
            found = float(offset[i + 1])
            continue
        if inner == outer:
            if inner == level:
                raise InvalidLightCurve(f"fraction {frac} lies on a flat part of a wing.")
            continue
        if (inner - level) * (outer - level) < 0:
            if found is not None:
                raise InvalidLightCurve(
                    f"fraction {frac} is crossed more than once on one wing."
                )
            share = (level - inner) / (outer - inner)
            found = float(offset[i] + share * (offset[i + 1] - offset[i]))
    if found is None:
        raise InvalidLightCurve(f"fraction {frac} is not reached on one wing.")
    return found
