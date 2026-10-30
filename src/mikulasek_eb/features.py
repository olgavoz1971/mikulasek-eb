"""Light-curve features for automatic model selection (see ``docs/feature_spec.md``).

The features are measured on the folded light curve without fitting a model. Each one is a
number with a unit and, where possible, an uncertainty. A feature that cannot be measured is
``None`` and the reason is added to ``Features.notes``; nothing is replaced by a default.

Conventions: phase in cycles, magnitudes in mag (larger = fainter), depth = how much
fainter than the reference level, widths in cycles, times in days unless the name says
otherwise. Phase 0 is the primary minimum found in the data. The reference level is the
mean of the median fold at the two quadratures (phases 0.25 and 0.75).
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import logging
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .api import validate_lightcurve
from .errors import ConfigError, InvalidLightCurve
from .exposure import estimate_cadence_days, estimate_exposure_days
from .model import prof1
from .solution import package_version, plain
from .tuning import noise_floor

logger = logging.getLogger(__name__)

FEATURES_SCHEMA_VERSION = 1
MAD_TO_SIGMA = 1.4826
MEDIAN_ERROR_FACTOR = 1.2533  # standard error of the median over that of the mean
FINE_GRID_N = 4001


@dataclass(frozen=True)
class FeatureSettings:
    """Measurement choices. All values are our own heuristics, to be tuned on real data.

    Attributes:
        n_bins (int | None): Phase bins of the median fold; None chooses
            ``min(1000, max(100, n_points // 10))`` (about 10 points per bin or more).
        min_points_per_bin (int): Bins with fewer points are empty.
        reference_halfwidth (float): Half-width (cycles) of the quadrature windows that
            define the reference level and the maxima.
        core_halfmax (float): Largest distance (cycles) from the guessed centre of bins used
            for the bottom parabola.
        vertex_core_frac (float): Bins deeper than this fraction of the depth enter the
            parabola that gives the eclipse centre.
        bottom_core_frac (float): Same, for the bottom curvature feature.
        edge_level (float): The eclipse edge is where the eclipse is this fraction of its
            depth (0.05 means the eclipse has fallen to 5 percent of its depth).
        window_margin (float): Eclipse window half-width = margin times half the edge width.
        flat_bottom_ratio (float): ``bottom_ratio`` above this suggests a flat bottom.
        separated_coverage (float): The eclipses are separated (a flat stretch exists between
            them) if the windows cover less than this fraction of the cycle.
        equal_depth_ratio (float): Depths are called equal if A2 / A1 is at least this.
        harmonics_max (int): Highest harmonic of the out-of-eclipse fit.
        gap_days (float): Gaps longer than this split the baseline segments.
        trend_max_order (int): Highest polynomial degree tried per segment.
        trend_window_periods (float): Time-bin width (periods) for the trend amplitude.
        outlier_nsig (float): Threshold of the outlier fraction (robust sigmas).
        gamma_match_tol (float): Largest accepted distance of the Gamma lookup match.
        parity_window_fwhm (float): Half-width (FWHM) of the odd/even depth windows.
    """

    n_bins: int | None = None
    min_points_per_bin: int = 3
    reference_halfwidth: float = 0.05
    core_halfmax: float = 0.25
    vertex_core_frac: float = 0.8
    bottom_core_frac: float = 0.9
    edge_level: float = 0.05
    window_margin: float = 1.0
    flat_bottom_ratio: float = 0.6
    separated_coverage: float = 0.3
    equal_depth_ratio: float = 0.8
    harmonics_max: int = 4
    gap_days: float = 0.5
    trend_max_order: int = 3
    trend_window_periods: float = 2.0
    outlier_nsig: float = 5.0
    gamma_match_tol: float = 0.15
    parity_window_fwhm: float = 0.25


GROUPS = (
    "sampling",
    "geometry",
    "shape",
    "variation",
    "ephemeris",
    "residuals",
)


@dataclass(frozen=True)
class Features:
    """Measured features of one light curve (plain data, JSON round trip).

    Attributes:
        sampling (dict): A group of the spec (counts, cadence, exposure, gaps, noise).
        geometry (dict): B group (depths, centres, widths).
        shape (dict): C group (bottom and wing shape, Gamma estimate).
        variation (dict): D group (out-of-eclipse variation, maxima, harmonics).
        ephemeris (dict): E group (parity test, drift between halves).
        residuals (dict): F group (fold scatter, trend).
        notes (list): ``"<feature>: <reason>"`` for each feature that is None.
        settings (dict): The ``FeatureSettings`` used.
        provenance (dict): period_days, epoch_jd, m0_mjd, package_version.
        schema_version (int): Version of the serialised layout.
    """

    sampling: dict
    geometry: dict
    shape: dict
    variation: dict
    ephemeris: dict
    residuals: dict
    notes: list
    settings: dict
    provenance: dict
    schema_version: int = FEATURES_SCHEMA_VERSION

    def to_dict(self) -> dict:
        """Plain-data form.

        Returns:
            dict: JSON-compatible dictionary.
        """
        return plain(dataclasses.asdict(self))

    @classmethod
    def from_dict(cls, data: dict) -> Features:
        """Rebuild from ``to_dict`` output.

        Args:
            data (dict): Dictionary produced by ``to_dict``.

        Returns:
            Features: The features.

        Raises:
            ConfigError: If keys are missing or unknown, or the schema differs.
        """
        names = {f.name for f in dataclasses.fields(cls)}
        if set(data) != names:
            raise ConfigError(f"Features keys differ: {sorted(set(data) ^ names)}")
        if data["schema_version"] != FEATURES_SCHEMA_VERSION:
            raise ConfigError("Unsupported Features schema version.")
        return cls(**data)

    def to_json(self, indent: int | None = 2) -> str:
        """JSON text of the features.

        Args:
            indent (int | None): JSON indent; None gives one line.

        Returns:
            str: JSON text.
        """
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True)

    @classmethod
    def from_json(cls, text: str) -> Features:
        """Rebuild from ``to_json`` output.

        Args:
            text (str): JSON text.

        Returns:
            Features: The features.
        """
        return cls.from_dict(json.loads(text))

    def features_hash(self) -> str:
        """Deterministic 12-character hash (package version excluded).

        Returns:
            str: Hex digest prefix.
        """
        d = self.to_dict()
        d["provenance"] = {
            k: v for k, v in d["provenance"].items() if k != "package_version"
        }
        text = json.dumps(d, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(text.encode()).hexdigest()[:12]


# ------------------------------------------------------------------ helpers
def _dphi(a, b):
    """Signed phase difference a - b wrapped to [-0.5, 0.5).

    Args:
        a (ndarray or float): Phase(s) in cycles.
        b (ndarray or float): Reference phase(s) in cycles.

    Returns:
        ndarray or float: Wrapped difference.
    """
    return (np.asarray(a) - np.asarray(b) + 0.5) % 1.0 - 0.5


@dataclass
class _Fold:
    """Median fold in phase bins."""

    centres: np.ndarray
    med: np.ndarray
    err: np.ndarray
    n: np.ndarray
    ok: np.ndarray

    def at(self, phi):
        """Periodic linear interpolation of the filled bins.

        Args:
            phi (ndarray or float): Phase(s) in cycles.

        Returns:
            ndarray: Folded magnitude.
        """
        return np.interp(
            np.mod(phi, 1.0), self.centres[self.ok], self.med[self.ok], period=1.0
        )


def _make_fold(ph: np.ndarray, mag: np.ndarray, nb: int, min_n: int) -> _Fold:
    """Median fold of points at phases ``ph`` in [0, 1).

    Args:
        ph (ndarray): Phases in [0, 1).
        mag (ndarray): Magnitudes.
        nb (int): Number of bins.
        min_n (int): Minimum points per bin.

    Returns:
        _Fold: Bin medians with the error of the median (robust scatter).
    """
    idx = np.minimum((ph * nb).astype(int), nb - 1)
    order = np.argsort(idx, kind="stable")
    idx_s, mag_s = idx[order], mag[order]
    edges = np.searchsorted(idx_s, np.arange(nb + 1))
    med = np.full(nb, np.nan)
    err = np.full(nb, np.nan)
    n = np.diff(edges)
    for k in range(nb):
        if n[k] >= min_n:
            v = mag_s[edges[k] : edges[k + 1]]
            m = np.median(v)
            med[k] = m
            mad = MAD_TO_SIGMA * np.median(np.abs(v - m))
            err[k] = MEDIAN_ERROR_FACTOR * mad / np.sqrt(n[k])
    ok = np.isfinite(med)
    if ok.sum() < 0.8 * nb:
        raise InvalidLightCurve(
            f"Only {int(ok.sum())} of {nb} phase bins are filled; too few points "
            f"for the fold (or the period is wrong)."
        )
    fill = np.nanmedian(err)
    err = np.where(np.isfinite(err) & (err > 0), err, fill)
    return _Fold((np.arange(nb) + 0.5) / nb, med, err, n, ok)


def _wls(X: np.ndarray, y: np.ndarray, sig: np.ndarray):
    """Weighted linear least squares by matrix algebra.

    Args:
        X (ndarray): Design matrix.
        y (ndarray): Data.
        sig (ndarray): Data errors.

    Returns:
        tuple: (coefficients, covariance, chi2, dof).
    """
    A = X / sig[:, None]
    b = y / sig
    coef, *_ = np.linalg.lstsq(A, b, rcond=None)
    cov = np.linalg.inv(A.T @ A)
    chi2 = float(np.sum((b - A @ coef) ** 2))
    return coef, cov, chi2, len(y) - X.shape[1]


def _level(fold: _Fold, centre: float, st: FeatureSettings) -> float:
    """Reference level: mean of the median fold at the two quadratures.

    Args:
        fold (_Fold): Median fold.
        centre (float): Phase of the primary minimum.
        st (FeatureSettings): Settings.

    Returns:
        float: Level (mag).
    """
    out = []
    for q in (0.25, 0.75):
        sel = fold.ok & (
            np.abs(_dphi(fold.centres, centre + q)) <= st.reference_halfwidth
        )
        if not sel.any():
            raise InvalidLightCurve("No filled bins near the quadrature phases.")
        out.append(np.median(fold.med[sel]))
    return float(np.mean(out))


def _measure_eclipse(fold: _Fold, nominal: float, level: float, st: FeatureSettings):
    """Depth and bottom parabola of the eclipse expected near ``nominal``.

    Args:
        fold (_Fold): Median fold.
        nominal (float): Expected phase of the minimum.
        level (float): Reference level (mag).
        st (FeatureSettings): Settings.

    Returns:
        dict: depth, peak_phase, and (when a convex-down parabola of at least 5 bins
        could be fitted) vertex, vertex_err, curvature, curvature_err, n_core, else
        those entries are None and ``reason`` explains.
    """
    grid = nominal + np.linspace(-0.15, 0.15, 3001)
    y = fold.at(grid)
    i = int(np.argmax(y))
    peak_phase, depth = float(grid[i]), float(y[i] - level)
    out = dict(
        depth=depth,
        peak_phase=peak_phase,
        vertex=None,
        vertex_err=None,
        curvature=None,
        curvature_err=None,
        n_core=0,
        reason=None,
    )
    if depth <= 0:
        out["reason"] = "no minimum deeper than the reference level"
        return out
    out.update(_core_parabola(fold, peak_phase, level, depth, st.vertex_core_frac, st))
    return out


def _core_parabola(fold, centre, level, depth, frac, st):
    """Fit a parabola to the bins deeper than ``frac`` of the depth around ``centre``.

    Args:
        fold (_Fold): Median fold.
        centre (float): Guess of the minimum phase.
        level (float): Reference level (mag).
        depth (float): Eclipse depth (mag).
        frac (float): Bins fainter than level + frac * depth are used.
        st (FeatureSettings): Settings.

    Returns:
        dict: vertex, vertex_err (phase offset from ``centre`` plus the guess), curvature
        and curvature_err (mag per cycle squared), n_core, reason (None if fitted).
    """
    d = _dphi(fold.centres, centre)
    sel = fold.ok & (np.abs(d) <= st.core_halfmax) & (fold.med >= level + frac * depth)
    res = dict(vertex=None, vertex_err=None, curvature=None, curvature_err=None)
    res.update(n_core=int(sel.sum()), reason=None)
    if sel.sum() < 5:
        res["reason"] = f"only {int(sel.sum())} bins in the bottom core (need 5)"
        return res
    x = d[sel]
    X = np.c_[np.ones_like(x), x, x**2]
    coef, cov, chi2, dof = _wls(X, fold.med[sel], fold.err[sel])
    scale = max(1.0, chi2 / max(dof, 1))
    cov = cov * scale
    c0, b, c = coef
    res["curvature"] = float(c)
    res["curvature_err"] = float(np.sqrt(cov[2, 2]))
    if c >= 0:
        res["reason"] = "bottom parabola is not convex towards fainter magnitudes"
        return res
    v = -b / (2 * c)
    g = np.array([0.0, -1 / (2 * c), b / (2 * c**2)])
    res["vertex"] = float(centre + v)
    res["vertex_err"] = float(np.sqrt(g @ cov @ g))
    return res


def _width(fold: _Fold, centre: float, level_mag: float):
    """Width of the region fainter than ``level_mag`` around ``centre`` (cycles).

    Args:
        fold (_Fold): Median fold.
        centre (float): Phase of the minimum.
        level_mag (float): Magnitude defining the edge.

    Returns:
        float or None: Width, None if the curve does not cross the level within half a cycle.
    """
    x = np.linspace(-0.5, 0.5, FINE_GRID_N)
    y = fold.at(centre + x)
    i0 = FINE_GRID_N // 2
    if y[i0] <= level_mag:
        return None
    edges = []
    for step in (1, -1):
        j = i0
        while 0 < j < FINE_GRID_N - 1 and y[j] > level_mag:
            j += step
        if y[j] > level_mag:
            return None
        k = j - step
        frac = (y[k] - level_mag) / (y[k] - y[j])
        edges.append(x[k] + frac * (x[j] - x[k]))
    return float(edges[0] - edges[1])


@lru_cache(maxsize=1)
def _gamma_table():
    """Width ratios W90/W50 and W10/W50 of profile 1 (C = 0) for a grid of Gamma.

    Returns:
        tuple: (Gamma grid, W90/W50, W10/W50).
    """
    gammas = np.geomspace(0.3, 30.0, 150)
    x = np.linspace(0.0, 60.0, 120001)
    r90, r10 = [], []
    for g in gammas:
        y = np.asarray(prof1(x, 1.0, g), float)
        w = {}
        for q in (0.9, 0.5, 0.1):
            j = int(np.argmax(y < q))
            frac = (y[j - 1] - q) / (y[j - 1] - y[j])
            w[q] = 2.0 * (x[j - 1] + frac * (x[j] - x[j - 1]))
        r90.append(w[0.9] / w[0.5])
        r10.append(w[0.1] / w[0.5])
    return gammas, np.array(r90), np.array(r10)


def gamma_from_width_ratios(bottom_ratio: float, wing_ratio: float):
    """Gamma of profile 1 (C = 0) whose width ratios best match the measured ones.

    Args:
        bottom_ratio (float): W90 / W50.
        wing_ratio (float): W10 / W50.

    Returns:
        tuple: (Gamma, distance of the match in the ratio plane).
    """
    gammas, r90, r10 = _gamma_table()
    dist = np.hypot(r90 - bottom_ratio, r10 - wing_ratio)
    i = int(np.argmin(dist))
    return float(gammas[i]), float(dist[i])


def _median_error(v: np.ndarray) -> float:
    """Standard error of the median from the robust scatter.

    Args:
        v (ndarray): Values.

    Returns:
        float: Error of the median.
    """
    return float(
        MEDIAN_ERROR_FACTOR
        * MAD_TO_SIGMA
        * np.median(np.abs(v - np.median(v)))
        / np.sqrt(len(v))
    )


# ---------------------------------------------------------------- the extractor
def extract_features(
    t_mjd,
    mag,
    sig,
    period_days: float,
    epoch_jd: float,
    settings: FeatureSettings | None = None,
    *,
    exposure_days: float | None = None,
) -> Features:
    """Measure the features of one light curve.

    Args:
        t_mjd (ndarray): Times (MJD), strictly increasing.
        mag (ndarray): Magnitudes.
        sig (ndarray): Catalogue magnitude errors.
        period_days (float): Catalogue period (days).
        epoch_jd (float): Catalogue epoch (absolute JD), recorded in the provenance.
        settings (FeatureSettings | None): Measurement choices; None uses the defaults.
        exposure_days (float | None): Exposure time; None estimates it from the cadence.

    Returns:
        Features: The measurements; unmeasurable ones are None with a note.

    Raises:
        InvalidLightCurve: If the input is invalid or the fold is too sparse.
    """
    validate_lightcurve(t_mjd, mag, sig, period_days, epoch_jd)
    st = settings or FeatureSettings()
    t, y, s = (np.asarray(a, dtype=float) for a in (t_mjd, mag, sig))
    P = float(period_days)
    notes: list[str] = []

    def none(name: str, reason: str):
        """Record that a feature is not measurable and return None."""
        notes.append(f"{name}: {reason}")
        return None

    nb = st.n_bins or int(min(1000, max(100, len(t) // 10)))
    x = (t - t[0]) / P
    # --- fold, centred on the primary minimum
    coarse = _make_fold(x % 1.0, y, nb, st.min_points_per_bin)
    c_abs = float(coarse.centres[np.argmax(np.where(coarse.ok, coarse.med, -np.inf))])
    fold = _make_fold((x - c_abs) % 1.0, y, nb, st.min_points_per_bin)
    level = _level(fold, 0.0, st)
    prim = _measure_eclipse(fold, 0.0, level, st)
    if prim["vertex"] is not None:
        c_abs += prim["vertex"]
        fold = _make_fold((x - c_abs) % 1.0, y, nb, st.min_points_per_bin)
        level = _level(fold, 0.0, st)
        prim = _measure_eclipse(fold, 0.0, level, st)
    else:
        notes.append(f"primary centre refinement: {prim['reason']}")
    A1 = prim["depth"]
    if A1 <= 0:
        raise InvalidLightCurve("The folded curve has no primary eclipse.")
    sec = _measure_eclipse(fold, 0.5, level, st)
    A2 = sec["depth"] if sec["depth"] > 0 else None
    m0_mjd = t[0] + c_abs * P
    m0_mjd += P * round((t.mean() - m0_mjd) / P)

    # --- residuals about the fold, noise
    phi_pts = (x - c_abs) % 1.0
    resid = y - fold.at(phi_pts)
    sigma_noise = noise_floor(t, resid)
    mad_res = MAD_TO_SIGMA * np.median(np.abs(resid - np.median(resid)))

    # --- A: sampling
    try:
        cadence_d = estimate_cadence_days(t)
    except ValueError as ex:
        cadence_d = none("cadence_s", str(ex))
    if exposure_days is not None:
        exposure_d, exposure_src = float(exposure_days), "user"
    else:
        exposure_src = "estimated"
        try:
            exposure_d = estimate_exposure_days(t)
        except ValueError as ex:
            exposure_d = none("exposure_s", str(ex))
    dt = np.diff(t)
    gaps = dt[dt > st.gap_days]
    sampling = dict(
        n_points=len(t),
        span_days=float(t[-1] - t[0]),
        n_cycles=float((t[-1] - t[0]) / P),
        cadence_s=None if cadence_d is None else float(86400 * cadence_d),
        exposure_s=None if exposure_d is None else float(86400 * exposure_d),
        exposure_source=exposure_src,
        points_per_cycle=float(len(t) / ((t[-1] - t[0]) / P)),
        n_gaps=int(len(gaps)),
        gap_days_total=float(gaps.sum()) if len(gaps) else 0.0,
        sigma_noise_mag=float(sigma_noise),
        sigma_catalogue_mag=float(np.median(s)),
        noise_over_catalogue=float(sigma_noise / np.median(s)),
        outlier_fraction=float(
            np.mean(np.abs(resid - np.median(resid)) > st.outlier_nsig * mad_res)
        ),
        n_bins=nb,
        n_empty_bins=int((~fold.ok).sum()),
    )

    # --- B: geometry
    fw1 = _width(fold, 0.0, level + 0.5 * A1)
    fw2 = None
    sec_centre = sec["vertex"] if sec["vertex"] is not None else sec["peak_phase"]
    if A2 is None:
        none("A2", sec["reason"])
    else:
        fw2 = _width(fold, sec_centre, level + 0.5 * A2)
        if fw2 is None:
            none("fwhm2", "secondary does not reach half depth within half a cycle")
    if fw1 is None:
        none("fwhm1", "primary does not reach half depth within half a cycle")
    sec_off = (
        None
        if A2 is None
        else dict(value=float(_dphi(sec_centre, 0.5)), err=sec["vertex_err"])
    )
    if A2 is not None and sec["vertex"] is None:
        none("secondary_offset", f"error not available: {sec['reason']}")
    geometry = dict(
        m0_mjd=float(m0_mjd),
        level_mag=float(level),
        A1=float(A1),
        A2=None if A2 is None else float(A2),
        depth_ratio=None if A2 is None else float(A2 / A1),
        secondary_offset=sec_off,
        fwhm1=fw1,
        fwhm2=fw2,
        fwhm_ratio=None if fw1 is None or fw2 is None else float(fw2 / fw1),
        points_per_fwhm=(
            None if fw1 is None else float(sampling["points_per_cycle"] * fw1)
        ),
        exposure_over_fwhm=(
            None if fw1 is None or exposure_d is None else float(exposure_d / (fw1 * P))
        ),
    )

    # --- C: shape
    w = {q: _width(fold, 0.0, level + q * A1) for q in (0.9, 0.5, 0.1)}
    bottom_ratio = None if w[0.9] is None or w[0.5] is None else float(w[0.9] / w[0.5])
    wing_ratio = None if w[0.1] is None or w[0.5] is None else float(w[0.1] / w[0.5])
    if bottom_ratio is None:
        none("bottom_ratio", "W90 or W50 of the primary not measurable")
    if wing_ratio is None:
        none("wing_ratio", "W10 or W50 of the primary not measurable")
    core = _core_parabola(fold, 0.0, level, A1, st.bottom_core_frac, st)
    if core["curvature"] is None:
        none("bottom_curvature", core["reason"])
    gamma = dist = None
    if bottom_ratio is not None and wing_ratio is not None:
        gamma, dist = gamma_from_width_ratios(bottom_ratio, wing_ratio)
        if dist > st.gamma_match_tol:
            notes.append(
                f"gamma_estimate: poor match to profile 1 (distance {dist:.3f})"
            )
    flat = None
    if bottom_ratio is not None:
        flat = bool(
            bottom_ratio >= st.flat_bottom_ratio
            and (
                core["curvature"] is None
                or abs(core["curvature"]) < 3 * core["curvature_err"]
            )
        )
    shape = dict(
        bottom_ratio=bottom_ratio,
        wing_ratio=wing_ratio,
        bottom_curvature=(
            None
            if core["curvature"] is None
            else dict(
                value=core["curvature"],
                err=core["curvature_err"],
                n_bins=core["n_core"],
            )
        ),
        gamma_estimate=gamma,
        gamma_match_distance=dist,
        flat_bottom=flat,
    )

    # --- D: out-of-eclipse variation
    grid = np.linspace(0.0, 1.0, 2000, endpoint=False)
    coverage = ooe_ratio = harmonics = None
    edge1 = _width(fold, 0.0, level + st.edge_level * A1)
    edge2 = None if A2 is None else _width(fold, sec_centre, level + st.edge_level * A2)
    if A2 is not None:
        # an eclipse whose edge is not reached within half a cycle fills the cycle
        half1 = 0.5 if edge1 is None else st.window_margin * edge1 / 2
        half2 = 0.5 if edge2 is None else st.window_margin * edge2 / 2
        for nm, e in (("edge_width1", edge1), ("edge_width2", edge2)):
            if e is None:
                notes.append(
                    f"{nm}: eclipse does not fall to {st.edge_level:.0%} of its depth "
                    f"within half a cycle; its window is taken as the whole cycle"
                )
        inside = (np.abs(_dphi(grid, 0.0)) < half1) | (
            np.abs(_dphi(grid, sec_centre)) < half2
        )
        coverage = float(inside.mean())
        free = grid[~inside]
        if len(free) >= 0.05 * len(grid):
            ooe_ratio = float(np.ptp(fold.at(free)) / A1)
        else:
            none(
                "ooe_amplitude_over_depth",
                f"eclipse windows cover {coverage:.0%} of the cycle (need at most 95%)",
            )
        d1 = np.abs(_dphi(fold.centres, 0.0))
        d2 = np.abs(_dphi(fold.centres, sec_centre))
        use = fold.ok & (d1 >= half1) & (d2 >= half2)
        K = st.harmonics_max
        if use.sum() >= 2 * (2 * K + 1):
            th = 2 * np.pi * fold.centres[use]
            cols = [np.ones(use.sum())]
            for k in range(1, K + 1):
                cols += [np.cos(k * th), np.sin(k * th)]
            coef, cov, chi2, dof = _wls(np.array(cols).T, fold.med[use], fold.err[use])
            cov = cov * max(1.0, chi2 / max(dof, 1))
            harmonics = {
                str(k): dict(
                    cos=float(coef[2 * k - 1]),
                    cos_err=float(np.sqrt(cov[2 * k - 1, 2 * k - 1])),
                    sin=float(coef[2 * k]),
                    sin_err=float(np.sqrt(cov[2 * k, 2 * k])),
                )
                for k in range(1, K + 1)
            }
        else:
            none("harmonic_amplitudes", "too few bins outside the eclipse windows")
    else:
        none("eclipse_coverage", "no secondary eclipse, windows not defined")
    separated = None if coverage is None else bool(coverage < st.separated_coverage)
    equal_depths = None if A2 is None else bool(A2 / A1 >= st.equal_depth_ratio)
    contact_candidate = (
        None if separated is None else bool(not separated and equal_depths)
    )
    hw = st.reference_halfwidth
    q = []
    for ph0 in (0.25, 0.75):
        sel = fold.ok & (np.abs(_dphi(fold.centres, ph0)) <= hw)
        q.append(
            (np.median(fold.med[sel]), np.sqrt(np.sum(fold.err[sel] ** 2)) / sel.sum())
        )
    mdiff = float(q[0][0] - q[1][0])
    merr = float(np.hypot(q[0][1], q[1][1]) * MEDIAN_ERROR_FACTOR)
    variation = dict(
        eclipse_coverage=coverage,
        ooe_amplitude_over_depth=ooe_ratio,
        edge_width1=edge1,
        edge_width2=edge2,
        separated_eclipses=separated,
        continuous_variation=None if separated is None else not separated,
        equal_depths=equal_depths,
        contact_candidate=contact_candidate,
        maxima_difference=dict(value=mdiff, err=merr),
        maxima_significance=float(abs(mdiff) / merr),
        harmonic_amplitudes=harmonics,
    )

    # --- E: period and ephemeris sanity
    parity = {}
    for name, centre, fw in (("primary", 0.0, fw1), ("secondary", sec_centre, fw2)):
        if fw is None:
            parity[name] = none(f"parity_{name}", "width not measurable")
            continue
        dphi_all = _dphi(phi_pts, centre)
        cyc = np.round((x - c_abs) - centre)
        near = np.abs(dphi_all) < st.parity_window_fwhm * fw
        vals = {}
        for p in (0, 1):
            sel = near & (cyc % 2 == p)
            if sel.sum() < 5:
                vals = None
                break
            vals[p] = (float(np.median(y[sel])), _median_error(y[sel]))
        if vals is None:
            parity[name] = none(
                f"parity_{name}", "fewer than 5 points in a parity class"
            )
            continue
        diff = vals[1][0] - vals[0][0]
        err = float(np.hypot(vals[0][1], vals[1][1]))
        parity[name] = dict(
            even_mag=vals[0][0],
            odd_mag=vals[1][0],
            difference=float(diff),
            significance=float(abs(diff) / err),
        )
    tmed = float(np.median(t))
    halves = {}
    for label, sel in (("first", t <= tmed), ("second", t > tmed)):
        nb_h = max(20, nb // 2)
        try:
            fh = _make_fold(
                ((x[sel] - c_abs) % 1.0), y[sel], nb_h, st.min_points_per_bin
            )
            lev = _level(fh, 0.0, st)
        except InvalidLightCurve as ex:
            halves[label] = none(f"drift_{label}_half", str(ex))
            continue
        res = {}
        for name, nominal in (("primary", 0.0), ("secondary", sec_centre)):
            m = _measure_eclipse(fh, nominal, lev, st)
            res[name] = m if m["vertex"] is not None else None
        halves[label] = dict(res, t_med=float(np.median(t[sel])))
    drift = {}
    for name in ("primary", "secondary"):
        a = halves.get("first")
        b = halves.get("second")
        if not a or not b or a[name] is None or b[name] is None:
            drift[name] = none(f"drift_{name}", "centre not measurable in both halves")
            continue
        shift = (b[name]["vertex"] - a[name]["vertex"]) * P * 86400.0
        err = float(
            np.hypot(a[name]["vertex_err"], b[name]["vertex_err"]) * P * 86400.0
        )
        ncyc = (b["t_med"] - a["t_med"]) / P
        drift[name] = dict(
            shift_s=float(shift),
            shift_err_s=err,
            period_correction_days=float(shift / 86400.0 / ncyc),
            period_correction_err_days=float(err / 86400.0 / ncyc),
        )
    ephemeris = dict(parity=parity, drift_between_halves=drift)

    # --- F: residual structure
    seg = np.cumsum(np.r_[0, dt > st.gap_days])
    width = st.trend_window_periods * P
    amp = []
    for sg in np.unique(seg):
        m = seg == sg
        tb = np.floor((t[m] - t[m][0]) / width).astype(int)
        meds = [
            np.median(resid[m][tb == b]) for b in np.unique(tb) if (tb == b).sum() >= 3
        ]
        if len(meds) >= 3:
            amp.append(float(np.ptp(meds)))
    bic = {}
    for order in range(st.trend_max_order + 1):
        chi2 = 0.0
        for sg in np.unique(seg):
            m = seg == sg
            u = 2 * (t[m] - t[m][0]) / max(np.ptp(t[m]), 1e-9) - 1
            X = np.vander(u, order + 1, increasing=True)
            c, *_ = np.linalg.lstsq(X, resid[m], rcond=None)
            chi2 += float(np.sum(((resid[m] - X @ c) / sigma_noise) ** 2))
        bic[order] = chi2 + (order + 1) * len(np.unique(seg)) * np.log(len(t))
    residuals = dict(
        fold_scatter_ratio=float(mad_res / sigma_noise),
        trend_amplitude_mag=(
            max(amp) if amp else none("trend_amplitude", "segments too short")
        ),
        trend_order_hint=int(min(bic, key=bic.get)),
        trend_bic={str(k): float(v) for k, v in bic.items()},
    )

    prov = dict(
        period_days=P,
        epoch_jd=float(epoch_jd),
        m0_mjd=float(m0_mjd),
        package_version=package_version(),
    )
    feats = Features(
        sampling=sampling,
        geometry=geometry,
        shape=shape,
        variation=variation,
        ephemeris=ephemeris,
        residuals=residuals,
        notes=notes,
        settings=dataclasses.asdict(st),
        provenance=prov,
    )
    return Features.from_dict(feats.to_dict())
