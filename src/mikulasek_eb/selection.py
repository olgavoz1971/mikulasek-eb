"""Evidence stage: fit the candidates and choose between them (see ``docs/feature_spec.md``).

Each candidate is fitted with the staged pipeline. A candidate whose
Levenberg-Marquardt fit is marked ``bound-active`` is fitted again with the
trust-region solver (``FitPolicy.retry_bound_active``, default true). The
comparison uses a Gaussian likelihood
with a free variance (residual-based sigma, not the catalogue one), evaluated on the points
that every successful candidate kept. Residuals of TESS-like data are strongly correlated,
so the number of points in the likelihood and in the penalty is the effective number of
independent points ``n_eff`` from the features:

    BIC = n_eff * ln(rss / n) + k * ln(n_eff)

with ``k`` the number of fitted parameters including the baselines. If the correlation time
could not be measured, ``n_eff`` is the number of points and the result says so.
"""

from __future__ import annotations

import dataclasses
import logging
from dataclasses import dataclass, field

import numpy as np

from .api import FitOutput, fit, validate_lightcurve
from .candidates import Candidate, RuleSettings, propose_candidates, review_features
from .errors import FitDidNotConverge, MikulasekError
from .features import FeatureSettings, Features, extract_features

logger = logging.getLogger(__name__)

DELTA_BIC_DECISIVE = 10.0


@dataclass
class CandidateRow:
    """One line of the comparison table.

    Attributes:
        name (str): Candidate label.
        config_hash (str): Hash of the candidate configuration.
        status (str): ``"ok"``, the fit status with its reasons when that is not ok
            (such a candidate is not compared), or the exception text.
        rms_mmag (float | None): Rms on the common points, mmag.
        k (int | None): Number of fitted parameters including baselines.
        n_common (int | None): Points used in the comparison.
        bic (float | None): Evidence value.
        delta_bic (float | None): BIC minus the best BIC.
        two_component (bool | None): Whether the fit adopted the second component.
    """

    name: str
    config_hash: str
    status: str
    rms_mmag: float | None = None
    k: int | None = None
    n_common: int | None = None
    bic: float | None = None
    delta_bic: float | None = None
    two_component: bool | None = None


@dataclass
class SelectionResult:
    """Outcome of ``select_model``.

    Attributes:
        best (Candidate | None): Chosen candidate, None if every fit failed.
        output (FitOutput | None): Its fit.
        table (list): ``CandidateRow`` per candidate, in the order tried.
        ask_human (bool): True if the choice is not clear or the features raise a flag.
        reasons (list): Why a human should look (empty if ``ask_human`` is False).
        remarks (list): Informative remarks from the features and the comparison.
        features (Features): The measured features.
        n_eff (float): Effective number of points used in the BIC.
        n_eff_measured (bool): False if ``n_eff`` is just the number of points.
        candidates (list): The candidates that were tried.
        fits (dict): ``FitOutput`` of every candidate that fitted, by name (not serialised).
    """

    best: Candidate | None
    output: FitOutput | None
    table: list
    ask_human: bool
    reasons: list
    remarks: list
    features: Features
    n_eff: float
    n_eff_measured: bool = True
    candidates: list = field(default_factory=list)
    fits: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        """Light, JSON-compatible summary (no arrays, no fit objects).

        Returns:
            dict: Choice, table, flags and the features hash.
        """
        return dict(
            best=None if self.best is None else self.best.name,
            best_config=None if self.best is None else self.best.config.to_dict(),
            table=[row.__dict__ for row in self.table],
            ask_human=self.ask_human,
            reasons=self.reasons,
            remarks=self.remarks,
            n_eff=self.n_eff,
            n_eff_measured=self.n_eff_measured,
            features_hash=self.features.features_hash(),
        )


def _bic(rss: float, n: int, k: int, n_eff: float) -> float:
    """Information criterion with the effective number of independent points.

    Args:
        rss (float): Residual sum of squares on the common points (mag squared).
        n (int): Number of common points.
        k (int): Number of fitted parameters.
        n_eff (float): Effective number of independent points.

    Returns:
        float: ``n_eff ln(rss / n) + k ln(n_eff)``.
    """
    return float(n_eff * np.log(rss / n) + k * np.log(n_eff))


def select_model(
    t_mjd,
    mag,
    sig,
    period_days: float,
    epoch_jd: float,
    *,
    features: Features | None = None,
    candidates: list | None = None,
    rules: RuleSettings | None = None,
    feature_settings: FeatureSettings | None = None,
    data_id: str | None = None,
    progress=None,
    retry_not_converged: bool = False,
) -> SelectionResult:
    """Measure features, fit the proposed candidates and choose by evidence.

    The first fit of each candidate uses the solver in its configuration, which
    is Levenberg-Marquardt unless ``bounded_solver`` is set. When that fit is
    marked ``bound-active`` and ``retry_bound_active`` is true, the same
    candidate is fitted again with the trust-region solver. The verdict of the
    second fit is the one that is compared. A second fit that raises is logged,
    and the bound-active result is kept. ``retry_not_converged`` is off: a
    candidate that raises ``FitDidNotConverge`` is not retried unless the caller
    asks.

    Args:
        t_mjd (ndarray): Times (MJD), strictly increasing.
        mag (ndarray): Magnitudes.
        sig (ndarray): Catalogue errors.
        period_days (float): Catalogue period (days).
        epoch_jd (float): Catalogue epoch (absolute JD).
        features (Features | None): Precomputed features; None measures them.
        candidates (list | None): Candidates to fit; None proposes them from the features.
        rules (RuleSettings | None): Thresholds of the rules.
        feature_settings (FeatureSettings | None): Used if the features are measured here.
        data_id (str | None): Identifier for provenance.
        progress (callable | None): ``progress(stage, fraction)``; may raise ``FitCancelled``.
        retry_not_converged (bool): When a Levenberg-Marquardt fit raises
            ``FitDidNotConverge``, fit that candidate again with the trust-region
            solver. The default is false. This is an experiment switch, not the
            selection default.

    Returns:
        SelectionResult: Choice, comparison table and flags.

    Raises:
        InvalidLightCurve: If the input is invalid.
    """
    validate_lightcurve(t_mjd, mag, sig, period_days, epoch_jd)
    rules = rules or RuleSettings()
    feats = features or extract_features(
        t_mjd, mag, sig, period_days, epoch_jd, feature_settings
    )
    ask, remarks = review_features(feats, rules)
    cands = list(
        candidates if candidates is not None else propose_candidates(feats, rules)
    )
    if progress is not None:
        progress("features", 0.1)
    fits: dict = {}
    table = []
    for i, cand in enumerate(cands):
        try:
            out = fit(
                t_mjd, mag, sig, period_days, epoch_jd, cand.config, data_id=data_id
            )
            verdict = out.solution.stats
            retry_note = ""
            if (
                verdict["status"] == "bound-active"
                and cand.config.policy.retry_bound_active
                and not cand.config.policy.bounded_solver
            ):
                logger.info(
                    "candidate %s is bound-active: refitting with the trust-region solver",
                    cand.name,
                )
                bounded = dataclasses.replace(
                    cand.config,
                    policy=dataclasses.replace(cand.config.policy, bounded_solver=True),
                ).validate()
                try:
                    out = fit(
                        t_mjd,
                        mag,
                        sig,
                        period_days,
                        epoch_jd,
                        bounded,
                        data_id=data_id,
                    )
                    cand = Candidate(
                        cand.name,
                        bounded,
                        cand.reasons
                        + ("trust-region fallback after a bound-active verdict",),
                    )
                    cands[i] = cand
                    verdict = out.solution.stats
                except MikulasekError as ex:
                    retry_note = (
                        f"; trust-region retry failed: {type(ex).__name__}: {ex}"
                    )
                    logger.warning(
                        "candidate %s trust-region retry failed: %s: %s",
                        cand.name,
                        type(ex).__name__,
                        ex,
                    )
            fits[cand.name] = out
            if verdict["status"] == "ok":
                row_status = "ok"
            else:
                row_status = (
                    f"{verdict['status']}: "
                    + "; ".join(verdict["status_reasons"])
                    + retry_note
                )
                logger.warning("candidate %s is %s", cand.name, row_status)
            table.append(CandidateRow(cand.name, cand.config.config_hash(), row_status))
        except MikulasekError as ex:
            if (
                retry_not_converged
                and isinstance(ex, FitDidNotConverge)
                and not cand.config.policy.bounded_solver
            ):
                logger.info(
                    "candidate %s did not converge: refitting with the trust-region solver",
                    cand.name,
                )
                bounded = dataclasses.replace(
                    cand.config,
                    policy=dataclasses.replace(cand.config.policy, bounded_solver=True),
                ).validate()
                try:
                    out = fit(
                        t_mjd,
                        mag,
                        sig,
                        period_days,
                        epoch_jd,
                        bounded,
                        data_id=data_id,
                    )
                except MikulasekError as retry_ex:
                    logger.warning(
                        "candidate %s trust-region retry failed: %s: %s",
                        cand.name,
                        type(retry_ex).__name__,
                        retry_ex,
                    )
                    table.append(
                        CandidateRow(
                            cand.name,
                            cand.config.config_hash(),
                            f"{type(ex).__name__}: {ex}; trust-region retry failed: "
                            f"{type(retry_ex).__name__}: {retry_ex}",
                        )
                    )
                else:
                    cand = Candidate(
                        cand.name,
                        bounded,
                        cand.reasons
                        + ("trust-region fallback after the fit did not converge",),
                    )
                    cands[i] = cand
                    fits[cand.name] = out
                    verdict = out.solution.stats
                    if verdict["status"] == "ok":
                        row_status = "ok"
                    else:
                        row_status = f"{verdict['status']}: " + "; ".join(
                            verdict["status_reasons"]
                        )
                        logger.warning("candidate %s is %s", cand.name, row_status)
                    table.append(
                        CandidateRow(cand.name, cand.config.config_hash(), row_status)
                    )
            else:
                logger.warning(
                    "candidate %s failed: %s: %s", cand.name, type(ex).__name__, ex
                )
                table.append(
                    CandidateRow(
                        cand.name,
                        cand.config.config_hash(),
                        f"{type(ex).__name__}: {ex}",
                    )
                )
        if progress is not None:
            progress(f"candidate {cand.name}", 0.1 + 0.9 * (i + 1) / len(cands))
    ok = [row for row in table if row.status == "ok"]
    n_eff = None
    measured = True
    if not ok:
        reasons = ask + ["every candidate failed to fit"]
        return SelectionResult(
            None,
            None,
            table,
            True,
            reasons,
            remarks,
            feats,
            float("nan"),
            False,
            cands,
            fits,
        )
    # common points
    common = None
    for row in ok:
        t_used = fits[row.name].t_mjd
        common = t_used if common is None else np.intersect1d(common, t_used)
    n = len(common)
    n_eff = feats.residuals["n_effective"]
    if n_eff is None:
        n_eff, measured = float(n), False
        remarks.append(
            "correlation time not measured: BIC uses the number of points as n_eff "
            "(residuals are correlated, so the evidence is overstated)"
        )
        measured_reason = (
            "BIC used n_eff = number of points (correlation time not measured)"
        )
    n_eff = float(min(n_eff, n))
    for row in ok:
        out = fits[row.name]
        r = out.pipeline.adopted
        sel = np.isin(out.t_mjd, common)
        rss = float(np.sum(r.resid[sel] ** 2))
        row.n_common = n
        row.k = int(len(out.t_mjd) - r.dof)
        row.rms_mmag = 1e3 * float(np.sqrt(rss / n))
        row.bic = _bic(rss, n, row.k, n_eff)
        row.two_component = bool(out.solution.stats["two_component"])
    best_bic = min(row.bic for row in ok)
    for row in ok:
        row.delta_bic = row.bic - best_bic
    ordered = sorted(ok, key=lambda row: row.bic)
    best_name = ordered[0].name
    reasons = list(ask)
    if not measured:
        reasons.append(measured_reason)
    if len(ordered) > 1 and ordered[1].delta_bic < DELTA_BIC_DECISIVE:
        reasons.append(
            f"{ordered[1].name} is within {ordered[1].delta_bic:.1f} BIC of {best_name} "
            f"(decisive difference is {DELTA_BIC_DECISIVE:g})"
        )
    if len(ok) < len(table):
        reasons.append(
            "failed candidates: " + ", ".join(r.name for r in table if r.status != "ok")
        )
    best = next(c for c in cands if c.name == best_name)
    return SelectionResult(
        best,
        fits[best_name],
        table,
        bool(reasons),
        reasons,
        remarks,
        feats,
        n_eff,
        measured,
        cands,
        fits,
    )
