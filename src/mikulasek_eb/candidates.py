"""Candidate configurations from measured features (draft rules of ``docs/feature_spec.md``).

The rules only decide which candidates are worth fitting; the choice between them is made by
fitted evidence (``selection.select_model``). Every threshold is a field of ``RuleSettings``
and is our own heuristic, set on the test stars only.
"""

from __future__ import annotations

import dataclasses
import logging
from dataclasses import dataclass

from .config import Config, DataSetup, FitPolicy, ModelSpec
from .features import Features

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class RuleSettings:
    """Thresholds of the rules.

    Attributes:
        harmonic_sigma (float): A harmonic is kept if its cosine amplitude exceeds this
            many errors (separated eclipses only).
        ocon_sigma (float): O'Connell term (S1) is included above this maxima significance.
        exposure_over_fwhm_min (float): Exposure integration is switched on above this ratio.
        baseline_orders (tuple): Baseline polynomial degrees offered as candidates.
        parity_sigma (float): Odd/even depth difference above this suggests period 2P.
        offset_sigma (float): Secondary offset above this many errors suggests eccentricity.
        offset_min (float): and above this value (cycles).
        shallow_depth_ratio (float): A2 / A1 below this is reported as a shallow secondary.
        min_points_per_fwhm (float): Fewer points per FWHM is reported as unresolved.
        max_outlier_fraction (float): More outliers than this is reported.
        drift_sigma (float): Drift between halves above this many errors is reported.
        two_component_candidates (bool): Also propose a twin of every candidate with the
            second profile component forced (``_2c``), so the evidence decides. Off until
            the synthetic calibration is finished.
    """

    harmonic_sigma: float = 3.0
    ocon_sigma: float = 3.0
    exposure_over_fwhm_min: float = 0.1
    baseline_orders: tuple = (0, 1)
    parity_sigma: float = 3.0
    offset_sigma: float = 3.0
    offset_min: float = 0.01
    shallow_depth_ratio: float = 0.15
    min_points_per_fwhm: float = 8.0
    max_outlier_fraction: float = 0.02
    drift_sigma: float = 3.0
    two_component_candidates: bool = False


@dataclass(frozen=True)
class Candidate:
    """One configuration to fit, with the reasons it was proposed.

    Attributes:
        name (str): Short label, for example ``"h123_b0"``.
        config (Config): The configuration.
        reasons (tuple): Plain-language reasons from the features.
    """

    name: str
    config: Config
    reasons: tuple


def _significant_harmonics(features: Features, rules: RuleSettings):
    """Harmonic orders whose cosine amplitude is significant outside the eclipses.

    Args:
        features (Features): Measured features.
        rules (RuleSettings): Thresholds.

    Returns:
        tuple or None: Sorted orders, None if the amplitudes were not measured.
    """
    amps = features.variation["harmonic_amplitudes"]
    if amps is None:
        return None
    return tuple(
        int(k)
        for k, a in sorted(amps.items(), key=lambda kv: int(kv[0]))
        if abs(a["cos"]) > rules.harmonic_sigma * a["cos_err"]
    )


def review_features(features: Features, rules: RuleSettings | None = None):
    """Reasons why a human should look at this light curve, and informative remarks.

    Args:
        features (Features): Measured features.
        rules (RuleSettings | None): Thresholds; None uses the defaults.

    Returns:
        tuple: (reasons to ask a human, remarks). Both are lists of strings.
    """
    r = rules or RuleSettings()
    ask, info = [], []
    g, e, sm = features.geometry, features.ephemeris, features.sampling
    if g["A2"] is None:
        ask.append("no secondary eclipse was found")
    elif g["depth_ratio"] < r.shallow_depth_ratio:
        info.append(f"shallow secondary (A2/A1 = {g['depth_ratio']:.2f})")
    for kind, p in e["parity"].items():
        if p is not None and p["significance"] > r.parity_sigma:
            ask.append(
                f"odd and even {kind} eclipses differ by {p['significance']:.1f} sigma: "
                f"the true period may be 2P"
            )
    off = g["secondary_offset"]
    if off is not None and off["err"] is not None:
        if (
            abs(off["value"]) > r.offset_min
            and abs(off["value"]) > r.offset_sigma * off["err"]
        ):
            ask.append(
                f"secondary is {off['value']:+.4f} cycles from phase 0.5: eccentric orbit?"
            )
    if (
        g["points_per_fwhm"] is not None
        and g["points_per_fwhm"] < r.min_points_per_fwhm
    ):
        ask.append(
            f"eclipses are poorly resolved ({g['points_per_fwhm']:.1f} points per FWHM)"
        )
    if sm["outlier_fraction"] > r.max_outlier_fraction:
        info.append(f"{100 * sm['outlier_fraction']:.1f}% of the points are outliers")
    shifts = {
        k: d
        for k, d in e["drift_between_halves"].items()
        if d is not None and abs(d["shift_s"]) > r.drift_sigma * d["shift_err_s"]
    }
    if len(shifts) == 2 and (
        shifts["primary"]["shift_s"] * shifts["secondary"]["shift_s"] > 0
    ):
        info.append(
            "both eclipses shift the same way between the halves: period error?"
        )
    elif shifts:
        info.append(
            "eclipse centres shift between the halves ("
            + ", ".join(f"{k} {d['shift_s']:+.0f} s" for k, d in shifts.items())
            + "); opposite signs would point to spots"
        )
    return ask, info


def propose_candidates(features: Features, rules: RuleSettings | None = None) -> list:
    """Candidate configurations for the staged pipeline (profile 1).

    Args:
        features (Features): Measured features.
        rules (RuleSettings | None): Thresholds; None uses the defaults.

    Returns:
        list[Candidate]: Candidates, each with reasons.
    """
    r = rules or RuleSettings()
    v, sm, g = features.variation, features.sampling, features.geometry
    reasons_common = []
    # harmonics
    if v["separated_eclipses"]:
        sig = _significant_harmonics(features, r)
        if sig is None:
            sets = [(2,)]
            reasons_common.append("separated eclipses, harmonics not measured: (2,)")
        else:
            sets = [sig]
            reasons_common.append(f"separated eclipses: significant harmonics {sig}")
    elif v["contact_candidate"]:
        sets = [(1, 2, 3), (2, 3)]
        reasons_common.append("equal depths and no flat stretch: contact-like")
    else:
        sets = [(1, 2, 3), (2,)]
        reasons_common.append("continuous variation, unequal depths")
    # O'Connell
    ocon = v["maxima_significance"] > r.ocon_sigma
    reasons_common.append(
        f"maxima differ by {v['maxima_significance']:.1f} sigma: S1 {'on' if ocon else 'off'}"
    )
    # exposure
    expo_days = None
    if (
        g["exposure_over_fwhm"] is not None
        and g["exposure_over_fwhm"] > r.exposure_over_fwhm_min
    ):
        expo_days = sm["exposure_s"] / 86400.0
        reasons_common.append(
            f"exposure is {g['exposure_over_fwhm']:.2f} of the FWHM: integrate"
        )
    gap = features.settings["gap_days"] if sm["n_gaps"] else None
    out = []
    for harmonics in sets:
        for order in r.baseline_orders:
            name = "h" + "".join(map(str, harmonics)) + f"_b{order}"
            cfg = Config(
                ModelSpec(profile=1, harmonics=harmonics, ocon=ocon),
                DataSetup(
                    baseline_order=order, gap_split_days=gap, exposure_days=expo_days
                ),
                FitPolicy(allow_two_component=True),
            ).validate()
            out.append(Candidate(name, cfg, tuple(reasons_common)))
            if r.two_component_candidates:
                twin = dataclasses.replace(
                    cfg,
                    policy=dataclasses.replace(
                        cfg.policy, allow_two_component=True, force_two_component=True
                    ),
                ).validate()
                out.append(
                    Candidate(
                        name + "_2c",
                        twin,
                        tuple(reasons_common) + ("second profile component forced",),
                    )
                )
    return out
