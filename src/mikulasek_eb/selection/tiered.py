"""Tiered search: baseline, then profile, then harmonics, by BIC.

A named selection strategy. It does not fit all combinations. It fits three small steps and
carries the winner of each step into the next one:

1. Baseline. The richest harmonic set, the proposed baseline orders, the start profile.
2. Profile. The winner again, the other profiles, and (start profile 1 only) profile 1 with
   the curvature ``C`` kept free, for flat-bottomed minima.
3. Harmonics. The winner again, and the smaller harmonic sets proposed from the features.

The search stops before step 3 when step 1 and step 2 each gained less than the decisive BIC
difference. No survey numbers are held here: the rules and the features are arguments.
Every candidate is made with the trust-region retry of the fit policy switched off.
"""

from __future__ import annotations

import dataclasses
import logging
import time
from dataclasses import dataclass, field
from typing import Callable, Sequence

from ..errors import NoValidSolution
from .candidates import Candidate, RuleSettings, propose_candidates
from .features import Features
from .search import DELTA_BIC_DECISIVE, SelectionResult, select_model

logger = logging.getLogger(__name__)

STEP_BASELINE = "1 baseline"
STEP_PROFILE = "2 profile"
STEP_HARMONICS = "3 harmonics"


@dataclass
class TieredStep:
    """One finished step of the tiered search.

    Attributes:
        label (str): ``"1 baseline"``, ``"2 profile"`` or ``"3 harmonics"``.
        result (SelectionResult): The comparison of the candidates of this step.
        seconds (float): Wall time of the step.
        gain (float): BIC of the previous winner (fitted again in this step) minus BIC of
            the winner of this step. For step 1 the reference is the first candidate.
    """

    label: str
    result: SelectionResult
    seconds: float
    gain: float


@dataclass
class TieredResult:
    """Outcome of ``tiered_search``.

    Attributes:
        best (Candidate): The winner of the last step that ran.
        final (SelectionResult): Selection result of the last step. ``final.output`` is the
            fit of the winner.
        steps (list[TieredStep]): The steps that ran, in order.
        stopped_early (str): Why step 3 was skipped, empty when it was not.
    """

    best: Candidate
    final: SelectionResult
    steps: list[TieredStep] = field(default_factory=list)
    stopped_early: str = ""

    @property
    def gains(self) -> list[float]:
        """BIC gain of each step that ran.

        Returns:
            list[float]: One value per step, in order.
        """
        return [step.gain for step in self.steps]

    @property
    def seconds(self) -> float:
        """Total wall time of the steps.

        Returns:
            float: Seconds.
        """
        return float(sum(step.seconds for step in self.steps))


def make_candidate(
    base: Candidate,
    profile: int,
    harmonics: Sequence[int],
    order: int,
    keep_curvature: bool = False,
) -> Candidate:
    """Make one candidate from a proposed one, with the trust-region retry off.

    Args:
        base (Candidate): Proposed candidate, for the data layer and the reasons.
        profile (int): Table 2 profile.
        harmonics (Sequence[int]): Harmonic orders.
        order (int): Baseline polynomial degree.
        keep_curvature (bool): Profile 1 only. Keep C free, so that the bottom of an eclipse
            can be flat.

    Returns:
        Candidate: Named ``p<profile>_h<harmonics>_b<order>``, with ``_C`` added when the
        curvature is kept.

    Raises:
        ConfigError: If the resulting configuration is not valid.
    """
    cfg = dataclasses.replace(
        base.config,
        model=dataclasses.replace(
            base.config.model, profile=profile, harmonics=tuple(harmonics)
        ),
        data=dataclasses.replace(base.config.data, baseline_order=order),
        policy=dataclasses.replace(
            base.config.policy, retry_bound_active=False, keep_curvature=keep_curvature
        ),
    ).validate()
    name = f"p{profile}_h{''.join(map(str, harmonics))}_b{order}"
    if keep_curvature:
        name += "_C"
    return Candidate(name, cfg, base.reasons + (f"Table 2 profile {profile}",))


def _gain(result: SelectionResult, previous_name: str) -> float:
    """BIC of the previous winner minus BIC of the step winner.

    Args:
        result (SelectionResult): Selection result of the step.
        previous_name (str): Previous winner, fitted again in this step.

    Returns:
        float: The gain. Zero when the winner is unchanged.

    Raises:
        NoValidSolution: If the previous winner has no BIC in this step.
    """
    by_name = {row.name: row for row in result.table}
    previous = by_name[previous_name]
    if previous.bic is None:
        raise NoValidSolution(f"{previous_name} did not finish again in this step")
    return float(previous.bic - by_name[result.best.name].bic)


def tiered_search(
    t_mjd,
    mag,
    sig,
    period_days: float,
    epoch_jd: float,
    *,
    features: Features,
    rules: RuleSettings,
    profiles: Sequence[int],
    start_profile: int,
    decisive_bic: float = DELTA_BIC_DECISIVE,
    data_id: str | None = None,
    on_step: Callable[[str, SelectionResult, float], None] | None = None,
) -> TieredResult:
    """Choose a model in three steps: baseline, profile, harmonics.

    Args:
        t_mjd (ndarray): Times (MJD), strictly increasing.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        period_days (float): Catalogue period (days).
        epoch_jd (float): Catalogue epoch (absolute JD).
        features (Features): Measured features of this light curve.
        rules (RuleSettings): Rules that turn the features into harmonic sets, baseline
            orders and a base candidate. Supplied by the survey adapter.
        profiles (Sequence[int]): Table 2 profiles that step 2 may try.
        start_profile (int): Profile of step 1. Profile 1 also adds the curvature trial.
        decisive_bic (float): A step gains when it lowers the BIC by at least this much.
            The search stops before step 3 when steps 1 and 2 both gained less.
        data_id (str | None): Identifier for provenance.
        on_step (callable | None): ``on_step(label, result, seconds)``, called after each
            step in which a candidate finished, for example to draw the candidates.

    Returns:
        TieredResult: The winner, the final comparison, every step and the stop reason.

    Raises:
        NoValidSolution: If no candidate of a step finished with status ``ok``, or the
            previous winner did not finish again.
        ConfigError: If a candidate configuration is not valid.
    """
    proposed = propose_candidates(features, rules)
    sets = sorted(
        {c.config.model.harmonics for c in proposed}, key=lambda h: (-len(h), h)
    )
    orders = sorted({c.config.data.baseline_order for c in proposed})
    base = proposed[0]
    steps: list[TieredStep] = []

    def run(
        label: str, candidates: list[Candidate], previous_name: str
    ) -> SelectionResult:
        """Fit one step, report it, and record its gain.

        Args:
            label (str): Step label.
            candidates (list[Candidate]): Candidates of the step, previous winner first.
            previous_name (str): Name used as the reference of the gain.

        Returns:
            SelectionResult: The comparison of the step.

        Raises:
            NoValidSolution: If no candidate finished with status ``ok``.
        """
        start = time.perf_counter()
        result = select_model(
            t_mjd,
            mag,
            sig,
            period_days,
            epoch_jd,
            features=features,
            candidates=candidates,
            rules=rules,
            data_id=data_id,
            retry_not_converged=False,
        )
        seconds = time.perf_counter() - start
        if result.best is None:
            raise NoValidSolution(
                f"{data_id} step {label}: no candidate finished with status ok"
            )
        if on_step is not None:
            on_step(label, result, seconds)
        steps.append(TieredStep(label, result, seconds, _gain(result, previous_name)))
        logger.info(
            "%s step %s winner %s in %.0f s", data_id, label, result.best.name, seconds
        )
        return result

    step1 = [make_candidate(base, start_profile, sets[0], order) for order in orders]
    result = run(STEP_BASELINE, step1, step1[0].name)
    winner = result.best
    order = winner.config.data.baseline_order

    step2 = (
        [winner]
        + [
            make_candidate(base, profile, sets[0], order)
            for profile in profiles
            if profile != start_profile
        ]
        + (
            [make_candidate(base, 1, sets[0], order, keep_curvature=True)]
            if start_profile == 1
            else []
        )
    )
    result = run(STEP_PROFILE, step2, winner.name)
    winner = result.best

    stop = ""
    if steps[-1].gain < decisive_bic and steps[-2].gain < decisive_bic:
        stop = "two steps gained less than the decisive BIC difference"

    if not stop and len(sets) > 1:
        profile = winner.config.model.profile
        keep = winner.config.policy.keep_curvature
        step3 = [winner] + [
            make_candidate(base, profile, harmonics, order, keep_curvature=keep)
            for harmonics in sets[1:]
        ]
        result = run(STEP_HARMONICS, step3, winner.name)
        winner = result.best

    return TieredResult(best=winner, final=result, steps=steps, stopped_early=stop)
