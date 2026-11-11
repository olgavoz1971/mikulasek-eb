"""Checks of the synthetic TESS-like light curves made for calibration.

A synthetic file must carry its spoiling record, the physics labels and the truth as
named parameters (see ``docs/prompt_synthetic_lc_metadata.md``). ``check_synthetic``
verifies that the items are present, that they agree with the data, and that the
truth columns are consistent with the stated noise and trends. It never repairs a
file: every problem is reported as a ``Finding``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

PASS, WARN, FAIL, INFO, SKIP = "pass", "warn", "fail", "info", "skip"

REQUIRED_NUMERIC = (
    "cadence_s",
    "exposure_s",
    "smear_frac",
    "duration_days",
    "gap_start_time",
    "gap_end_time",
    "rel_sigma",
    "sigma_factor",
    "seed",
    "P_true_days",
    "P_meta_days",
    "period_eps",
    "trend_amp_frac",
    "eclipse_depth_mag",
    "n_eclipses",
)
REQUIRED_TEXT = (
    "product_kind",
    "recipe_id",
    "physics_id",
    "source_file",
    "noise_type",
    "trend_types",
    "trend_model",
    "time_scale",
    "time_refposition",
    "eb_class",
    "has_oconnell",
    "has_spots",
    "eccentric",
    "total_eclipse",
    "eclipse_E",
    "eclipse_kind",
    "eclipse_T_conj_days",
    "eclipse_T_min_pure_days",
    "eclipse_O_C_injected_days",
)
EB_CLASSES = ("detached", "semi_detached", "contact")
BOOLEAN_PARAMS = ("has_oconnell", "has_spots", "eccentric", "total_eclipse")

# Tolerances (documented here so that they are not hidden in the checks).
CADENCE_REL_TOL = 1e-6
NOISE_RATIO_TOL = 0.10  # relative difference of the noise rms from its stated value
TREND_AMP_REL_TOL = 0.25
DEPTH_MAG_TOL = 0.005
MINIMUM_OFFSET_WARN_S = 1.0  # |T_conj - T_min_pure| for circular, unspotted systems
SMEAR_TRUTH_REL_TOL = 0.005  # of the pure depth
SMEAR_SUBSAMPLES = 201


@dataclass(frozen=True)
class Finding:
    """One result of a check.

    Attributes:
        name (str): Short name of the check.
        status (str): ``pass``, ``warn``, ``fail``, ``info`` or ``skip``.
        message (str): What was found, with numbers.
    """

    name: str
    status: str
    message: str


@dataclass(frozen=True)
class SyntheticFile:
    """The parts of a synthetic file that the checks need.

    Attributes:
        time (ndarray): Time column (days from the declared origin).
        flux (ndarray): Observed flux.
        flux_err (ndarray): Flux errors.
        truth_flux (ndarray | None): Pure smeared flux at the observed times.
        truth_trend (ndarray | None): Trend alone.
        meta (dict): Table parameters (strings or numbers).
        timesys (dict): ``timescale``, ``refposition``, ``timeorigin`` of the file.
        problems (list): Problems found while reading (volightcurve checks).
        path (Path | None): File the data came from.
    """

    time: np.ndarray
    flux: np.ndarray
    flux_err: np.ndarray
    truth_flux: np.ndarray | None
    truth_trend: np.ndarray | None
    meta: dict
    timesys: dict
    problems: list = field(default_factory=list)
    path: Path | None = None


def load_synthetic(path: str | Path, with_errors: bool = True) -> SyntheticFile:
    """Read a synthetic light curve through ``volightcurve``.

    Args:
        path (str | Path): VOTable written by the generator.
        with_errors (bool): False for the pure curve, which has no error column.

    Returns:
        SyntheticFile: Arrays and parameters; ``problems`` lists what the package
        reported (unreadable file, missing or ambiguous columns, calibration or time
        problems).
    """
    from volightcurve import read_lightcurve

    path = Path(path)
    problems = []
    try:
        volc = read_lightcurve(path, filename=path.name)
    except Exception as ex:  # reported to the user, not hidden
        return SyntheticFile(
            *(np.array([]),) * 3,
            None,
            None,
            {},
            {},
            [f"volightcurve could not read the file: {type(ex).__name__}: {ex}"],
            path,
        )
    tcols = volc.get_time_colnames()
    fcols = volc.get_flux_colnames()
    ecols = volc.get_flux_error_colnames()
    needed_err = 1 if with_errors else 0
    if len(tcols) != 1 or len(fcols) != 1 or len(ecols) != needed_err:
        problems.append(
            f"column discovery: time {tcols}, flux {fcols}, flux error {ecols} "
            f"(exactly one time, one flux and {needed_err} error column needed)"
        )
        return SyntheticFile(
            *(np.array([]),) * 3, None, None, dict(volc.table.meta), {}, problems, path
        )
    problems += [f"conversion: {p}" for p in volc.conversion_problems(fcols[0])]
    problems += [f"time: {p}" for p in volc.time_problems()]
    table = volc.table

    def column(name):
        return np.asarray(table[name], float) if name in table.colnames else None

    ts = volc.timesys
    return SyntheticFile(
        np.asarray(table[tcols[0]], float),
        np.asarray(table[fcols[0]], float),
        np.asarray(table[ecols[0]], float) if ecols else np.array([]),
        column("truth_flux_pure_smeared"),
        column("truth_trend"),
        dict(table.meta),
        dict(
            timescale=ts.timescale, refposition=ts.refposition, timeorigin=ts.timeorigin
        ),
        problems,
        path,
    )


def _num(meta: dict, key: str):
    """Numeric parameter or None.

    Args:
        meta (dict): Parameters.
        key (str): Name.

    Returns:
        float | None: The value, None if absent or not a number.
    """
    try:
        return float(meta[key])
    except (KeyError, TypeError, ValueError):
        return None


def _flag(meta: dict, key: str):
    """Boolean parameter written as text, or None.

    Args:
        meta (dict): Parameters.
        key (str): Name.

    Returns:
        bool | None: True or False for ``true``/``false`` (any case), else None.
    """
    value = str(meta.get(key, "")).strip().lower()
    return {"true": True, "false": False}.get(value)


def _array(meta: dict, key: str, cast=float):
    """Comma-separated parameter as an array, or None.

    Args:
        meta (dict): Parameters.
        key (str): Name.
        cast (callable): Conversion of one item.

    Returns:
        ndarray | None: Values, None if absent or unparsable.
    """
    if key not in meta:
        return None
    try:
        return np.array([cast(x) for x in str(meta[key]).split(",")])
    except ValueError:
        return None


def _pure_curve(pure: SyntheticFile, period: float):
    """Periodic interpolation of a pure one-orbit curve.

    Args:
        pure (SyntheticFile): The high-resolution curve (time, flux).
        period (float): Period (days).

    Returns:
        callable: ``f(t)`` returning the pure flux at times ``t``.
    """
    return lambda t: np.interp(np.mod(t, period), pure.time, pure.flux, period=period)


def _smeared(curve, t: np.ndarray, exposure_days: float, offset_days: float = 0.0):
    """Mean of a pure curve over an exposure window.

    Args:
        curve (callable): Pure flux as a function of time.
        t (ndarray): Time stamps (days).
        exposure_days (float): Window length (days).
        offset_days (float): Shift of the window centre from the stamp (days).

    Returns:
        ndarray: Averaged flux at each stamp.
    """
    u = (
        offset_days
        + (np.arange(SMEAR_SUBSAMPLES) / (SMEAR_SUBSAMPLES - 1) - 0.5) * exposure_days
    )
    return np.mean([curve(t + x) for x in u], axis=0)


def check_synthetic(data: SyntheticFile, pure: SyntheticFile | None = None) -> list:
    """Run all checks on one synthetic file.

    Args:
        data (SyntheticFile): The spoiled light curve and its parameters.
        pure (SyntheticFile | None): The pure high-resolution curve named by
            ``source_file``; enables the independent check of the smeared truth.

    Returns:
        list[Finding]: One or more findings per check, in a fixed order.
    """
    out: list[Finding] = []

    def add(name, status, message):
        out.append(Finding(name, status, message))

    meta = data.meta
    if data.problems:
        for p in data.problems:
            add("volightcurve", FAIL, p)
    else:
        add("volightcurve", PASS, "reads cleanly; one time, flux and error column")
    if len(data.time) == 0:
        return out

    # --- presence of the named parameters
    missing = [k for k in REQUIRED_NUMERIC if _num(meta, k) is None]
    missing += [k for k in REQUIRED_TEXT if str(meta.get(k, "")).strip() == ""]
    if missing:
        add("parameters", FAIL, "missing or unreadable: " + ", ".join(missing))
    else:
        add("parameters", PASS, "all required parameters are present")
    if _num(meta, "smear_nsub") is None:
        add(
            "smear_nsub",
            WARN,
            "smear_nsub is not recorded (how the exposure average was made)",
        )

    t, flux, err = data.time, data.flux, data.flux_err
    dt = np.diff(t)
    cadence = float(np.median(dt)) * 86400.0

    # --- cadence, exposure, duration
    c = _num(meta, "cadence_s")
    if c is not None:
        ok = abs(cadence - c) <= CADENCE_REL_TOL * c
        add(
            "cadence",
            PASS if ok else FAIL,
            f"declared {c:g} s, median time step {cadence:.6f} s",
        )
    e, f = _num(meta, "exposure_s"), _num(meta, "smear_frac")
    if None not in (c, e, f):
        ok = abs(e - f * c) <= CADENCE_REL_TOL * c
        add(
            "exposure",
            PASS if ok else FAIL,
            f"exposure_s {e:g}, smear_frac x cadence_s = {f * c:g}",
        )
    d = _num(meta, "duration_days")
    if d is not None:
        span = float(t[-1] - t[0])
        ok = span <= d and d - span <= 2 * cadence / 86400.0
        add(
            "duration",
            PASS if ok else WARN,
            f"declared {d:g} d, data span {span:.4f} d",
        )

    # --- gap
    gaps = np.where(dt > 5 * cadence / 86400.0)[0]
    g0, g1 = _num(meta, "gap_start_time"), _num(meta, "gap_end_time")
    if g0 is not None and g1 is not None and g1 > g0:
        if len(gaps) != 1:
            add(
                "gap",
                FAIL,
                f"declared one gap [{g0:g}, {g1:g}] d, data have {len(gaps)} gaps",
            )
        else:
            a, b = float(t[gaps[0]]), float(t[gaps[0] + 1])
            tol = 1.5 * cadence / 86400.0
            ok = abs(a - g0) <= tol and abs(b - g1) <= tol
            add(
                "gap",
                PASS if ok else FAIL,
                f"declared [{g0:g}, {g1:g}] d, data [{a:.4f}, {b:.4f}] d",
            )
    elif g0 is not None and g1 is not None:
        add(
            "gap",
            PASS if len(gaps) == 0 else FAIL,
            f"no gap declared, data have {len(gaps)}",
        )

    # --- period and time system
    pt, pm, eps = (
        _num(meta, "P_true_days"),
        _num(meta, "P_meta_days"),
        _num(meta, "period_eps"),
    )
    period = _num(meta, "period")
    # The truth (fluxes, conjunction times) is made with P_true; the table PERIOD is
    # P_meta, which differs by period_eps on purpose.
    truth_period = pt if pt is not None else (period or pm)
    if None not in (pt, pm, eps):
        ok = abs(abs(pm - pt) / pt - abs(eps)) <= 1e-9 + 1e-6 * abs(eps)
        add(
            "period",
            PASS if ok else WARN,
            f"P_true {pt:g}, P_meta {pm:g}, period_eps {eps:g}",
        )
    if period is not None and pm is not None:
        add(
            "period parameter",
            PASS if abs(period - pm) <= 1e-12 else FAIL,
            f"table period {period:g}, P_meta_days {pm:g}",
        )
    ts = data.timesys
    mismatches = []
    if str(meta.get("time_scale", "")) != str(ts.get("timescale")):
        mismatches.append(f"scale {meta.get('time_scale')} vs {ts.get('timescale')}")
    if str(meta.get("time_refposition", "")) != str(ts.get("refposition")):
        mismatches.append(
            f"position {meta.get('time_refposition')} vs {ts.get('refposition')}"
        )
    origin = _num(meta, "timeorigin_jd")
    if (
        origin is not None
        and ts.get("timeorigin") is not None
        and abs(origin - ts["timeorigin"]) > 1e-9
    ):
        mismatches.append(f"origin {origin} vs {ts['timeorigin']}")
    add(
        "time system",
        PASS if not mismatches else FAIL,
        (
            "the parameters repeat the TIMESYS of the file"
            if not mismatches
            else "; ".join(mismatches)
        ),
    )

    # --- labels
    cls = str(meta.get("eb_class", ""))
    if cls not in EB_CLASSES:
        add("labels", FAIL, f"eb_class {cls!r} is not one of {EB_CLASSES}")
    else:
        bad = [k for k in BOOLEAN_PARAMS if _flag(meta, k) is None]
        if bad:
            add("labels", FAIL, "not true/false: " + ", ".join(bad))
        elif cls == "contact" and _num(meta, "fillout_factor") is None:
            add("labels", FAIL, "contact system without fillout_factor")
        else:
            add("labels", PASS, f"class {cls}, flags readable")

    # --- noise, trends and the truth columns
    tp, tr = data.truth_flux, data.truth_trend
    if tp is None or tr is None:
        add(
            "truth columns",
            FAIL,
            "truth_flux_pure_smeared and/or truth_trend are missing",
        )
    else:
        depth_flux = float(tp.max() - tp.min())
        res = flux - tp - tr
        factor = _num(meta, "sigma_factor")
        ratio = float(np.std(res / err))
        if factor:
            expected = 1.0 / factor
            ok = abs(ratio - expected) <= NOISE_RATIO_TOL * expected
            add(
                "noise",
                PASS if ok else WARN,
                f"rms of (flux - truth - trend) / error is {ratio:.3f}, expected {expected:.3f} "
                f"(sigma_factor {factor:g}, additive trend)",
            )
        rel = _num(meta, "rel_sigma")
        if rel is not None:
            r = float(np.std(res / tp))
            add(
                "noise level",
                PASS if abs(r - rel) <= NOISE_RATIO_TOL * rel else WARN,
                f"relative noise {r:.4f}, rel_sigma {rel:g}",
            )
        amp, frac = float(np.ptp(tr)), _num(meta, "trend_amp_frac")
        if frac is not None and depth_flux > 0:
            realised = amp / depth_flux
            ok = abs(realised - frac) <= TREND_AMP_REL_TOL * frac
            add(
                "trend amplitude",
                PASS if ok else WARN,
                f"peak-to-peak trend / flux depth is {realised:.3f}, trend_amp_frac {frac:g}"
                + (
                    ""
                    if ok
                    else " (record the realised value and define the fraction exactly)"
                ),
            )
        dm = _num(meta, "eclipse_depth_mag")
        if dm is not None:
            measured = 2.5 * float(np.log10(tp.max() / tp.min()))
            ok = abs(measured - dm) <= DEPTH_MAG_TOL
            add(
                "eclipse depth",
                PASS if ok else WARN,
                f"eclipse_depth_mag {dm:.4f}, depth of truth_flux_pure_smeared {measured:.4f} mag",
            )

    # --- eclipse truth
    n = _num(meta, "n_eclipses")
    E = _array(meta, "eclipse_E", int)
    kind = (
        None
        if "eclipse_kind" not in meta
        else np.array(str(meta["eclipse_kind"]).split(","))
    )
    tc = _array(meta, "eclipse_T_conj_days")
    tm = _array(meta, "eclipse_T_min_pure_days")
    oc = _array(meta, "eclipse_O_C_injected_days")
    arrays = (E, kind, tc, tm, oc)
    if any(a is None for a in arrays):
        add("eclipse truth", FAIL, "eclipse arrays are missing or unreadable")
    else:
        lengths = {len(a) for a in arrays}
        if lengths != {int(n)} if n is not None else len(lengths) != 1:
            add(
                "eclipse truth",
                FAIL,
                f"array lengths {sorted(lengths)}, n_eclipses {n}",
            )
        elif not set(kind) <= {"primary", "secondary"}:
            add("eclipse truth", FAIL, f"unknown kinds {sorted(set(kind))}")
        else:
            add("eclipse truth", PASS, f"{len(E)} eclipses, arrays consistent")
            _check_eclipses(add, meta, t, cadence, truth_period, E, kind, tc, tm, oc)

    # --- the smeared truth against an independent average of the pure curve
    if pure is None or tp is None or truth_period is None or e is None:
        add("smeared truth", SKIP, "pure curve not available")
    else:
        curve = _pure_curve(pure, truth_period)
        depth = float(tp.max() - tp.min())
        exposure_days = e / 86400.0
        scores = {
            label: float(np.std(tp - _smeared(curve, t, exposure_days, shift)))
            for label, shift in (
                ("centred", 0.0),
                ("ending at the stamp", -0.5 * exposure_days),
                ("starting at the stamp", 0.5 * exposure_days),
            )
        }
        best = min(scores, key=scores.get)
        rel = scores["centred"] / depth
        status = PASS if best == "centred" and rel <= SMEAR_TRUTH_REL_TOL else FAIL
        add(
            "smeared truth",
            status,
            "rms difference from an independent average of the pure curve: "
            + ", ".join(f"{k} {v:.5f}" for k, v in scores.items())
            + f"; centred is {100 * rel:.2f} % of the depth, best window: {best}",
        )
    return out


def _check_eclipses(add, meta, t, cadence_s, period, E, kind, tc, tm, oc):
    """Checks of the per-eclipse truth arrays.

    Args:
        add (callable): Collector ``add(name, status, message)``.
        meta (dict): Parameters.
        t (ndarray): Observation times (days).
        cadence_s (float): Median time step (s).
        period (float | None): Period (days).
        E (ndarray): Cycle numbers.
        kind (ndarray): ``primary`` or ``secondary`` per eclipse.
        tc (ndarray): Geometric conjunction times (days).
        tm (ndarray): Minimum of the pure curve (days).
        oc (ndarray): Injected O-C (days).
    """
    primary = kind == "primary"
    off = (tm - tc) * 86400.0
    circular_clean = (
        _flag(meta, "eccentric") is False and _flag(meta, "has_spots") is False
    )
    for label, sel in (("primary", primary), ("secondary", ~primary)):
        if sel.any():
            add(
                f"T_min_pure - T_conj ({label})",
                (
                    INFO
                    if not circular_clean
                    else (
                        PASS
                        if np.max(np.abs(off[sel])) <= MINIMUM_OFFSET_WARN_S
                        else WARN
                    )
                ),
                f"mean {off[sel].mean():+.3f} s, scatter {off[sel].std():.3f} s, "
                f"largest {np.max(np.abs(off[sel])):.3f} s"
                + (
                    ""
                    if circular_clean
                    else " (eccentric or spotted: an offset is expected)"
                ),
            )
    stated = _num(meta, "conj_minus_min_pure_s")
    if stated is None:
        add("conj_minus_min_pure_s", WARN, "not recorded")
    else:
        worst = max(
            abs(off[primary].mean()) if primary.any() else 0.0,
            abs(off[~primary].mean()) if (~primary).any() else 0.0,
        )
        add(
            "conj_minus_min_pure_s",
            PASS if abs(stated) >= worst - 1e-3 else WARN,
            f"stated {stated:+.3f} s, largest mean offset over both kinds {worst:.3f} s"
            + (
                ""
                if abs(stated) >= worst - 1e-3
                else " (give the secondary or the larger value)"
            ),
        )
    if period is not None and _flag(meta, "eccentric") is False:
        epoch = _num(meta, "epoch") or 0.0
        expected = epoch + E * period + np.where(primary, 0.0, 0.5 * period)
        worst = float(np.max(np.abs(tc - expected))) * 86400.0
        add(
            "conjunction times",
            PASS if worst <= 1e-3 else FAIL,
            f"circular orbit: largest difference from epoch + E P (+ P/2) is {worst:.4f} s",
        )
    if np.any(oc != 0):
        add(
            "injected O-C",
            INFO,
            f"non-zero O-C injected, up to {np.max(np.abs(oc)) * 86400.0:.2f} s",
        )
    else:
        add("injected O-C", INFO, "zero for every eclipse")
    if period is not None:
        expected_n = 2.0 * (t[-1] - t[0]) / period
        add(
            "eclipse count",
            PASS if abs(len(E) - expected_n) <= 3 else WARN,
            f"{len(E)} listed, about {expected_n:.1f} expected from the span",
        )
    flags = ("eclipse_observed", "eclipse_in_gap")
    present = [k for k in flags if k in meta]
    gap_lo, gap_hi = _num(meta, "gap_start_time"), _num(meta, "gap_end_time")
    in_gap = (
        (tc > gap_lo) & (tc < gap_hi)
        if None not in (gap_lo, gap_hi)
        else np.zeros(len(tc), bool)
    )
    before = tc < t[0]
    after = tc > t[-1]
    note = f"{int(in_gap.sum())} in the gap, {int(before.sum() + after.sum())} outside the data range"
    if present:
        add("eclipse flags", PASS, f"{', '.join(present)} present ({note})")
    else:
        add(
            "eclipse flags",
            WARN,
            f"per-eclipse observed or in-gap flag missing ({note})",
        )
