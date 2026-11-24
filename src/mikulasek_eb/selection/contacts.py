"""Outer contacts of the eclipses on the folded light curve, from two finders and a rule.

The outer contact is the corner where an eclipse begins or ends on the median fold. Two
finders look for it on each side of each eclipse, and an agreement rule combines them. The
result replaces the 5 percent edge as the eclipse window when contacts are found, and the
5 percent edge stays, flagged, when they are not (see ``docs/feature_spec.md``).

* **Jump scan.** The slope outside minus the slope inside a candidate point, on two windows
  of ``slope_narrow`` times the full width at half depth (FWHM). The contact is the maximum
  of that jump. It is accepted only if the jump is ``min_sigma`` times its error, falls back
  to (almost) zero beyond the peak, keeps its place when the windows are widened to
  ``slope_wide`` FWHM, and is not at the edge of the scanned range.
* **Method A.** A broken line fitted to the fold in a window that starts at half depth and is
  ``k`` FWHM wide for every ``k`` in ``a_k_values``. The hinge of the line is the contact.
  The broken line is "method A" of Andronov and co-workers as implemented in Mavka:
  Andrych, Andronov and Chinarova (2020, J. Phys. Stud. 24, 1902); code of
  M. Yu. Pyatnytskyy (``lc_approx``, MIT). Here the model is re-implemented and used for a
  different aim, the corner of an eclipse and not its time of minimum, with a window that
  starts at half depth.
* **Agreement rule.** Per contact: both finders accept and differ by at most ``agree_tol`` of
  their mean (``confirmed``, jump-scan position); one finder accepts (``single``, its
  position); both accept and disagree (``disagree``, no contact); none (``none``). An eclipse
  has contacts when both its sides have one. If exactly one eclipse has contacts, the other
  one takes its half-width and is flagged ``borrowed`` (equal duration of both eclipses,
  which assumes a circular orbit). If none has, the flag is ``fallback_5pct``.

Nothing is moved by an empirical factor. The contacts found lie inside the geometric first
contact of a mesh model, by 7 to 19 percent of the half-width on synthetic curves; this is a
documented property of the method (the corner is visible only after the flux has fallen),
and it is not corrected.

All thresholds are in ``ContactSettings`` and have no defaults.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress

logger = logging.getLogger(__name__)

MEDIAN_ERROR_FACTOR = 1.2533  # standard error of the median over that of the mean


@dataclass(frozen=True)
class ContactSettings:
    """Thresholds of the contact finders and of the agreement rule. There are no defaults.

    Attributes:
        scan_limit_primary (float): Largest distance (cycles) from the primary centre used.
        scan_limit_secondary (float): The same for the secondary (the neighbouring eclipse
            is closer, so it is smaller).
        start_over_fwhm (float): The jump scan starts at this many FWHM from the centre.
        slope_narrow (float): Width of each of the two slope windows of the jump, in FWHM,
            for the contact itself.
        slope_wide (float): The same for the test of stability of the peak.
        min_points_in_window (int): Fewest bins in a slope window.
        min_sigma (float): Smallest jump over its error.
        max_tail (float): Largest ratio of the lowest jump beyond the peak to the peak.
        max_shift (float): Largest ratio of the wide-window peak distance to the narrow one.
        a_k_values (tuple[float, ...]): Method A window widths, in FWHM.
        a_edge_frac (float): The hinge must be this fraction of the window away from its ends.
        a_max_out_over_in (float): Largest ratio of the outer to the inner slope.
        a_max_spread (float): Largest range of the hinge over the windows, over its median.
        a_min_points (int): Fewest bins in a method A window.
        agree_tol (float): Largest relative difference (of the mean) for the two finders to agree.
    """

    scan_limit_primary: float
    scan_limit_secondary: float
    start_over_fwhm: float
    slope_narrow: float
    slope_wide: float
    min_points_in_window: int
    min_sigma: float
    max_tail: float
    max_shift: float
    a_k_values: tuple
    a_edge_frac: float
    a_max_out_over_in: float
    a_max_spread: float
    a_min_points: int
    agree_tol: float


def broken_line(t, c1, c2, c3, c4):
    """Broken linear model of method A.

    Args:
        t (ndarray): Distance from the eclipse centre (cycles).
        c1 (float): Magnitude at the hinge.
        c2 (float): Slope inside the hinge (mag per cycle).
        c3 (float): Slope parameter outside the hinge; the outside slope is ``-c3``.
        c4 (float): Hinge (cycles), the outer contact.

    Returns:
        ndarray: Model magnitudes.
    """
    t = np.asarray(t, float)
    return np.where(t < c4, c1 + c2 * (t - c4), c1 - c3 * (t - c4))


def _slope(x, y, s):
    """Ordinary least-squares slope and its error for known per-point errors.

    Args:
        x (ndarray): Positions.
        y (ndarray): Values.
        s (ndarray): One-sigma errors of ``y``.

    Returns:
        tuple[float, float]: Slope and its error.
    """
    fit = linregress(x, y)
    dx = x - x.mean()
    return fit.slope, float(np.sqrt(np.sum((s * dx) ** 2)) / np.sum(dx**2))


def _jump(distance, y, sig, at, width, min_points):
    """Slope outside minus slope inside at one candidate point.

    Args:
        distance (ndarray): Distance of the bins from the eclipse centre (cycles).
        y (ndarray): Bin medians (mag).
        sig (ndarray): Errors of the bin medians (mag).
        at (float): Candidate position (cycles from the centre).
        width (float): Window width on each side (cycles).
        min_points (int): Fewest bins in a window.

    Returns:
        tuple[float, float]: Jump (mag per cycle) and its error; NaN, NaN if a window has too
        few bins.
    """
    out = (distance >= at) & (distance < at + width)
    inn = (distance >= at - width) & (distance < at)
    if out.sum() < min_points or inn.sum() < min_points:
        return np.nan, np.nan
    so, eo = _slope(distance[out], y[out], sig[out])
    si, ei = _slope(distance[inn], y[inn], sig[inn])
    return so - si, float(np.hypot(eo, ei))


def _side_scan(distance, y, sig, fwhm, limit, slope_over_fwhm, st):
    """Jump along one side of an eclipse.

    Args:
        distance (ndarray): Distance of the bins from the eclipse centre (cycles).
        y (ndarray): Bin medians (mag).
        sig (ndarray): Errors of the bin medians (mag).
        fwhm (float): Full width at half depth (cycles).
        limit (float): Largest distance scanned (cycles).
        slope_over_fwhm (float): Half-window of the jump, in FWHM.
        st (ContactSettings): Settings.

    Returns:
        dict: ``d`` (positions), ``half`` (jump) and ``sigma_half`` (jump over its error).
    """
    width = 2.0 * slope_over_fwhm * fwhm
    step = float(np.median(np.diff(np.sort(distance))))
    grid = np.arange(st.start_over_fwhm * fwhm, limit, step)
    half = np.full(len(grid), np.nan)
    sigma = np.full(len(grid), np.nan)
    for i, d in enumerate(grid):
        j, e = _jump(distance, y, sig, d, 0.5 * width, st.min_points_in_window)
        half[i] = j
        sigma[i] = j / e if np.isfinite(e) and e > 0 else np.nan
    return dict(d=grid, half=half, sigma_half=sigma)


def _pick_corner(scan, scan_wide, st):
    """Contact of one side from the jump scan, and the reasons if it is not accepted.

    Args:
        scan (dict): ``_side_scan`` with the narrow window.
        scan_wide (dict): ``_side_scan`` with the wide window.
        st (ContactSettings): Settings.

    Returns:
        tuple[float | None, list[str]]: Accepted position (cycles) or None, and the failed
        tests (empty if accepted; a single reason if the scan has no positive jump).
    """
    half = scan["half"]
    if not np.isfinite(half).any() or np.nanmax(half) <= 0:
        return None, ["no positive jump"]
    i = int(np.nanargmax(half))
    j = int(np.nanargmax(scan_wide["half"]))
    tail = float(np.nanmin(half[i:]) / half[i])
    shift = float(scan_wide["d"][j] / scan["d"][i])
    sigma = float(scan["sigma_half"][i])
    reasons = []
    if sigma < st.min_sigma:
        reasons.append(f"jump only {sigma:.0f} sigma")
    if tail > st.max_tail:
        reasons.append(f"tail {tail:.2f} > {st.max_tail}")
    if shift > st.max_shift:
        reasons.append(f"peak moves outward x{shift:.2f} > {st.max_shift}")
    if i in (0, len(half) - 1):
        reasons.append("peak at the edge of the scan")
    return (None if reasons else float(scan["d"][i])), reasons


def _fit_window(distance, mag, start, stop, st):
    """Method A on one window of one side of an eclipse.

    Args:
        distance (ndarray): Distance of the bins from the eclipse centre (cycles).
        mag (ndarray): Bin medians (mag).
        start (float): Window start (cycles from the centre).
        stop (float): Window end (cycles from the centre).
        st (ContactSettings): Settings.

    Returns:
        dict | None: ``position`` (hinge) and ``reasons`` (failed edge or slope tests); None
        if the window has too few bins or the fit does not converge.
    """
    sel = (distance >= start) & (distance <= stop) & np.isfinite(mag)
    if sel.sum() < st.a_min_points:
        return None
    d, y = distance[sel], mag[sel]
    p0 = (float(np.mean(y)), -1.0, 0.0, 0.5 * (start + stop))
    try:
        p, _ = curve_fit(broken_line, d, y, p0=p0, maxfev=100000)
    except RuntimeError:
        return None
    _, slope_in, c3, hinge = (float(x) for x in p)
    slope_out = -c3
    ratio = abs(slope_out) / abs(slope_in) if slope_in != 0 else np.inf
    span = stop - start
    reasons = []
    if hinge < start + st.a_edge_frac * span:
        reasons.append(f"hinge {hinge:.3f} at the inner end of the window")
    if hinge > stop - st.a_edge_frac * span:
        reasons.append(f"hinge {hinge:.3f} at the outer end of the window")
    if slope_in >= 0:
        reasons.append(f"inside slope {slope_in:.1f} is not a descending wall")
    if ratio > st.a_max_out_over_in:
        reasons.append(f"|out|/|in| {ratio:.2f} > {st.a_max_out_over_in}")
    return dict(position=hinge, reasons=reasons)


def _method_a_side(distance, mag, fwhm, limit, st):
    """Method A on one side for every window width, and the verdict.

    Args:
        distance (ndarray): Distance of the bins from the eclipse centre (cycles).
        mag (ndarray): Bin medians (mag).
        fwhm (float): Full width at half depth (cycles).
        limit (float): Largest distance used (cycles).
        st (ContactSettings): Settings.

    Returns:
        tuple[float | None, list[str]]: Median hinge if all tests pass, else None, and the
        failed tests.
    """
    start = 0.5 * fwhm
    fits = [
        _fit_window(distance, mag, start, min(start + k * fwhm, limit), st)
        for k in st.a_k_values
    ]
    valid = [(k, f) for k, f in zip(st.a_k_values, fits) if f is not None]
    reasons = []
    if len(valid) < len(st.a_k_values):
        reasons.append(f"{len(st.a_k_values) - len(valid)} window(s) could not be fitted")
    for k, f in valid:
        reasons += [f"k={k:g}: {r}" for r in f["reasons"]]
    positions = np.array([f["position"] for _, f in valid])
    position = float(np.median(positions)) if positions.size else np.nan
    spread = (
        float((positions.max() - positions.min()) / position)
        if positions.size and position > 0
        else np.nan
    )
    if not np.isfinite(spread) or spread > st.a_max_spread:
        reasons.append(
            f"hinge moves with the window (spread {spread:.2f} > {st.a_max_spread})"
        )
    return (None if reasons else position), reasons


def _side_arrays(centres, medians, sig, centre, sign):
    """Bins on one side of an eclipse, ordered by distance from the centre.

    Args:
        centres (ndarray): Bin phases (cycles).
        medians (ndarray): Bin medians (mag), NaN for empty bins.
        sig (ndarray): Errors of the bin medians (mag).
        centre (float): Phase of the eclipse centre.
        sign (float): -1 for the left side, +1 for the right side.

    Returns:
        tuple[ndarray, ndarray, ndarray]: Distance, median and error.
    """
    signed = ((centres - centre + 0.5) % 1.0) - 0.5
    keep = np.isfinite(medians) & (sign * signed > 0)
    d = sign * signed[keep]
    order = np.argsort(d)
    return d[order], medians[keep][order], sig[keep][order]


def combine_contact(a_pos, jump_pos, agree_tol: float) -> dict:
    """Combine one contact from the two finders by the agreement rule.

    Args:
        a_pos (float | None): Distance from the eclipse centre found by method A, or None.
        jump_pos (float | None): Distance found by the jump scan, or None.
        agree_tol (float): Largest relative difference (of the mean) for agreement.

    Returns:
        dict: ``status`` (``confirmed``, ``single``, ``disagree`` or ``none``), ``position``
        (cycles or None), ``by`` (finder of a ``single`` contact, else None), ``a`` and ``jump``.
    """
    out = {"a": a_pos, "jump": jump_pos, "position": None, "by": None}
    if a_pos is None and jump_pos is None:
        out["status"] = "none"
    elif a_pos is None or jump_pos is None:
        out["status"] = "single"
        out["by"] = "jump" if a_pos is None else "A"
        out["position"] = jump_pos if a_pos is None else a_pos
    elif abs(a_pos - jump_pos) <= agree_tol * 0.5 * (a_pos + jump_pos):
        out["status"] = "confirmed"
        out["position"] = jump_pos
    else:
        out["status"] = "disagree"
    return out


def find_contacts(
    centres,
    medians,
    counts,
    sigma_point: float,
    centre1: float,
    centre2: float,
    fwhm1: float,
    fwhm2: float,
    settings: ContactSettings,
) -> dict:
    """Outer contacts of both eclipses on the median fold.

    Args:
        centres (ndarray): Phases of the fold bins (cycles).
        medians (ndarray): Bin medians (mag), NaN for empty bins.
        counts (ndarray): Number of points in each bin.
        sigma_point (float): Robust scatter of one point about the fold (mag). The error of a
            bin median is ``1.2533 * sigma_point / sqrt(count)``.
        centre1 (float): Phase of the primary centre.
        centre2 (float): Phase of the secondary centre.
        fwhm1 (float): Full width at half depth of the primary (cycles).
        fwhm2 (float): The same for the secondary (cycles).
        settings (ContactSettings): Thresholds.

    Returns:
        dict: ``primary`` and ``secondary``, each with ``left`` and ``right`` (the output of
        ``combine_contact`` plus ``a_reasons`` and ``jump_reasons``), ``half`` (half-width,
        cycles, or None) and ``flag`` (``measured``, ``borrowed`` or ``fallback_5pct``).
    """
    st = settings
    sig = MEDIAN_ERROR_FACTOR * sigma_point / np.sqrt(np.maximum(counts, 1))
    spec = {
        "primary": (centre1, fwhm1, st.scan_limit_primary),
        "secondary": (centre2, fwhm2, st.scan_limit_secondary),
    }
    out = {}
    for key, (centre, fwhm, limit) in spec.items():
        item = {}
        for side, sign in (("left", -1.0), ("right", 1.0)):
            d, y, e = _side_arrays(centres, medians, sig, centre, sign)
            scan = _side_scan(d, y, e, fwhm, limit, st.slope_narrow, st)
            wide = _side_scan(d, y, e, fwhm, limit, st.slope_wide, st)
            jump_pos, jump_reasons = _pick_corner(scan, wide, st)
            a_pos, a_reasons = _method_a_side(d, y, fwhm, limit, st)
            c = combine_contact(a_pos, jump_pos, st.agree_tol)
            c["a_reasons"], c["jump_reasons"] = a_reasons, jump_reasons
            item[side] = c
        ok = all(item[s]["position"] is not None for s in ("left", "right"))
        item["half"] = (
            0.5 * (item["left"]["position"] + item["right"]["position"]) if ok else None
        )
        item["flag"] = "measured" if ok else None
        out[key] = item
    measured = [k for k in out if out[k]["half"] is not None]
    for key in out:
        if out[key]["flag"] is None:
            if len(measured) == 1:
                out[key]["half"] = out[measured[0]]["half"]
                out[key]["flag"] = "borrowed"
            else:
                out[key]["flag"] = "fallback_5pct"
    logger.debug(
        "contacts: primary %s, secondary %s", out["primary"]["flag"], out["secondary"]["flag"]
    )
    return out
