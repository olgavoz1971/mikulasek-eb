"""Experiment on the 1800 s spot curves whose primary ran away.

The staged fit is applied to longitude 0 and 180 degrees. The second pass wraps
``fit_global``. A step is rejected when the template range on a fine phase grid
is more than three times the range at the data phases. The parameters of the
previous accepted step are kept. Pictures and a table go to
``experiments/runaway_1800/out/``.

Run from the repository root::

    MPLBACKEND=Agg python experiments/runaway_1800/run.py
"""

from __future__ import annotations

import copy
import csv
import dataclasses
import logging
from pathlib import Path

import numpy as np

import mikulasek_eb.extensions as extensions
import mikulasek_eb.fitting as fitting
import mikulasek_eb.pipeline as pipeline
import mikulasek_eb.tuning as tuning
from mikulasek_eb import fit, mjd_to_jd
from mikulasek_eb.candidates import propose_candidates
from mikulasek_eb.config import Config, FitPolicy, ModelSpec, TimingRecipe
from mikulasek_eb.features import extract_features
from mikulasek_eb.helpers import absolute_jd_to_mjd, predicted_magnitude
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.model import phase
from mikulasek_eb.plotting import plot_global_fit

logger = logging.getLogger(__name__)
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
FILES = (
    "results/detached/detached_none_nooconnell_SP_L000_SP_L000_1800_c1800.vot",
    "results/detached/detached_none_nooconnell_SP_L180_SP_L180_1800_c1800.vot",
)
# A real eclipse has almost the same range on the grid and at the samples.
# The needles on these two curves are hundreds of times larger on the grid.
PEAK_RATIO = 3.0
GRID_N = 8001
_REAL_FIT = fitting.fit_global
_GUARD = False
_STATE = {"last": None, "reverted": 0, "unresolved": False}


def _span(values: np.ndarray) -> float:
    """Difference between the largest and smallest finite value.

    Args:
        values (ndarray): Model values (mag).

    Returns:
        float: Range (mag). Zero when no value is finite.
    """
    finite = np.asarray(values, float)
    finite = finite[np.isfinite(finite)]
    if len(finite) == 0:
        return 0.0
    return float(np.max(finite) - np.min(finite))


def _ranges(model, result, t) -> tuple[float, float]:
    """Template range on a fine phase grid and at the phases of the data.

    Args:
        model (EBModel): Model that was fitted.
        result (FitResult): Parameters of this step.
        t (ndarray): Times of this step (days).

    Returns:
        tuple[float, float]: Grid range and sample range, in magnitudes.
    """
    ephem = result.ephem
    theta = phase(t, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0))
    phi = theta - np.round(theta)
    grid = np.linspace(-0.5, 0.5, GRID_N)
    return _span(model(grid, result.p)), _span(model(phi, result.p))


def _rejected(grid_span: float, sample_span: float) -> bool:
    """Decide whether the curve between the samples outgrows the samples.

    Args:
        grid_span (float): Template range on the fine grid (mag).
        sample_span (float): Template range at the data phases (mag).

    Returns:
        bool: True when the grid range exceeds ``PEAK_RATIO`` times the sample range.
    """
    return grid_span > PEAK_RATIO * max(sample_span, 1e-4)


def _guarded_fit_global(*args, **kwargs):
    """Fit one step, and keep the previous accepted step when the grid test fails.

    Args:
        *args: Positional arguments of ``fit_global``.
        **kwargs: Keyword arguments of ``fit_global``.

    Returns:
        FitResult: This step, or a copy of the last step that passed the grid test.
    """
    result = _REAL_FIT(*args, **kwargs)
    if not _GUARD:
        return result
    t = np.asarray(kwargs.get("t", args[0]), float)
    model = kwargs.get("model", args[3])
    grid_span, sample_span = _ranges(model, result, t)
    rank = bool(result.rank_deficient)
    if not _rejected(grid_span, sample_span):
        if rank:
            logger.warning(
                "rank deficient, grid %.4g mag agrees with samples %.4g mag; step kept",
                grid_span,
                sample_span,
            )
        _STATE["last"] = result
        return result
    logger.warning(
        "step rejected: grid %.4g mag, samples %.4g mag, rank_deficient=%s",
        grid_span,
        sample_span,
        rank,
    )
    if _STATE["last"] is None:
        logger.error("no previous accepted step; this solution stays unresolved")
        _STATE["unresolved"] = True
        return result
    _STATE["reverted"] += 1
    logger.warning("keeping the previous accepted step (%d reverted)", _STATE["reverted"])
    return copy.deepcopy(_STATE["last"])


def _install_guard() -> None:
    """Point every staged-fit import of ``fit_global`` at the wrapper."""
    fitting.fit_global = _guarded_fit_global
    pipeline.fit_global = _guarded_fit_global
    tuning.fit_global = _guarded_fit_global
    extensions.fit_global = _guarded_fit_global


def _reset_state() -> None:
    """Clear the accepted step before a new light curve."""
    _STATE["last"] = None
    _STATE["reverted"] = 0
    _STATE["unresolved"] = False


def _config_for(path: Path) -> Config:
    """Profile 1, H2, O'Connell on, straight baseline, this file's exposure and gaps.

    Args:
        path (Path): Light curve.

    Returns:
        Config: The structure used for the degraded spot run.
    """
    payload = load_tess_votable(path)
    t = absolute_jd_to_mjd(payload["t"])
    features = extract_features(
        t, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    own = propose_candidates(features)[0].config.data
    data = dataclasses.replace(own, baseline_order=1)
    return Config(
        ModelSpec(profile=1, harmonics=(2,), ocon=True),
        data,
        FitPolicy(allow_two_component=True),
        TimingRecipe(window=0.12),
    ).validate()


def _one(path: Path, label: str) -> dict:
    """Fit one curve and save the profile picture.

    Args:
        path (Path): Light curve.
        label (str): ``current`` or ``grid_check``.

    Returns:
        dict: Adopted parameters, ranges, and whether a step was reverted.
    """
    _reset_state()
    payload = load_tess_votable(path)
    t = absolute_jd_to_mjd(payload["t"])
    config = _config_for(path)
    out = fit(
        t, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"], config, data_id=path.stem
    )
    adopted = out.pipeline.adopted
    model = out.solution.build_model()
    yhat, template, theta = predicted_magnitude(
        out.t_mjd, out.mag, out.sig, model, adopted.p, adopted.ephem,
        out.pipeline.seg, config.data.baseline_order,
    )
    grid_span, sample_span = _ranges(model, adopted, out.t_mjd)
    picture = OUT / f"{path.stem}_{label}.png"
    plot_global_fit(
        mjd_to_jd(out.t_mjd), out.mag, out.sig, yhat, template, theta, model,
        out.pipeline.single_component.p, adopted.p, out.pipeline.seg,
        f"{path.stem} {label}  rms {1e3 * float(out.solution.stats['rms_mag']):.2f} mmag",
        timescale=payload.get("timescale"),
        save_path=str(picture),
    )
    row = dict(
        file=path.name,
        mode=label,
        A1=adopted.p.get("A1"),
        G1=adopted.p.get("G1"),
        D1=adopted.p.get("D1"),
        A2=adopted.p.get("A2"),
        G2=adopted.p.get("G2"),
        rms_mmag=1e3 * float(adopted.rms),
        grid_span_mag=grid_span,
        sample_span_mag=sample_span,
        rank_deficient=bool(adopted.rank_deficient),
        n_reverted=_STATE["reverted"],
        unresolved=_STATE["unresolved"],
        flagged=_rejected(grid_span, sample_span),
    )
    logger.info(
        "%s %s A1=%.4g G1=%.4g grid=%.3g sample=%.3g reverted=%d flagged=%s",
        path.stem, label, row["A1"], row["G1"], grid_span, sample_span,
        row["n_reverted"], row["flagged"],
    )
    return row


def main() -> None:
    """Fit both curves with the current solver and with the grid check."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    _install_guard()
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    global _GUARD
    for relative in FILES:
        path = GENERATED / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        _GUARD = False
        rows.append(_one(path, "current"))
        _GUARD = True
        rows.append(_one(path, "grid_check"))
    destination = OUT / "comparison.csv"
    with open(destination, "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    logger.info("wrote %s", destination)


if __name__ == "__main__":
    main()
