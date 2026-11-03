"""Experiment on the 1800 s spot curves whose primary ran away.

The second pass holds every eclipse-shape parameter that leaves a wide physical
interval, and refits the rest. The intervals are stated below. Pictures and a
table go to ``experiments/runaway_1800/out/``.

Run from the repository root::

    MPLBACKEND=Agg python experiments/runaway_1800/run.py
"""

from __future__ import annotations

import csv
import dataclasses
import logging
from pathlib import Path

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
from mikulasek_eb.plotting import plot_global_fit

logger = logging.getLogger(__name__)
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
FILES = (
    "results/detached/detached_none_nooconnell_SP_L000_SP_L000_1800_c1800.vot",
    "results/detached/detached_none_nooconnell_SP_L180_SP_L180_1800_c1800.vot",
)
# Wide intervals for this experiment, in magnitudes and in cycles. Gamma from 0.3
# still includes the paper's cusp (Figure 1 draws 1/2) and excludes the needle.
# The second-component amplitude may be negative. Harmonics are not bounded here.
SHAPE_BOUNDS = {
    "A1": (1e-4, 3.0),
    "A2": (1e-4, 3.0),
    "A1b": (-3.0, 3.0),
    "A2b": (-3.0, 3.0),
    "D1": (0.005, 0.25),
    "D2": (0.005, 0.25),
    "D1b": (0.005, 0.25),
    "D2b": (0.005, 0.25),
    "G1": (0.3, 6.0),
    "G2": (0.3, 6.0),
    "G1b": (0.3, 6.0),
    "G2b": (0.3, 6.0),
    "C1": (-2.0, 2.0),
    "C2": (-2.0, 2.0),
}
_REAL_FIT = fitting.fit_global
_GUARD = False
_STATE = {"held": 0}


def _outside(parameters: dict) -> dict:
    """Parameters that lie outside ``SHAPE_BOUNDS``.

    Args:
        parameters (dict): Template parameters of one step.

    Returns:
        dict: Name to the bound it is held at.
    """
    held = {}
    for name, (low, high) in SHAPE_BOUNDS.items():
        if name not in parameters:
            continue
        value = parameters[name]
        if value < low:
            held[name] = low
        elif value > high:
            held[name] = high
    return held


def _call_fit(args, kwargs, p0, free):
    """Call the real fit with a replaced starting template and free list.

    Args:
        args (tuple): Original positional arguments.
        kwargs (dict): Original keyword arguments.
        p0 (dict): Starting template.
        free (list): Template names still free.

    Returns:
        FitResult: Refit.
    """
    if "p0" in kwargs or "free" in kwargs:
        updated = dict(kwargs)
        updated["p0"] = p0
        updated["free"] = free
        return _REAL_FIT(*args, **updated)
    updated = list(args)
    updated[4] = p0
    updated[5] = free
    return _REAL_FIT(*updated, **kwargs)


def _guarded_fit_global(*args, **kwargs):
    """Fit one step, and repeat it with out-of-range shape parameters held.

    Args:
        *args: Positional arguments of ``fit_global``.
        **kwargs: Keyword arguments of ``fit_global``.

    Returns:
        FitResult: The free fit, or the refit with parameters held on the bound.
    """
    result = _REAL_FIT(*args, **kwargs)
    if not _GUARD:
        return result
    free0 = list(kwargs.get("free", args[5]))
    held_at = {}
    for _ in range(4):
        held = _outside(result.p)
        if not held:
            return result
        held_at.update(held)
        p0 = dict(result.p)
        p0.update(held_at)
        free = [name for name in free0 if name not in held_at]
        logger.warning("holding %s at the bound and refitting", held)
        _STATE["held"] += len(held)
        result = _call_fit(args, kwargs, p0, free)
    still = _outside(result.p)
    if still:
        logger.error("still outside after holding: %s", still)
    return result


def _install_guard() -> None:
    """Point every staged-fit import of ``fit_global`` at the wrapper."""
    fitting.fit_global = _guarded_fit_global
    pipeline.fit_global = _guarded_fit_global
    tuning.fit_global = _guarded_fit_global
    extensions.fit_global = _guarded_fit_global


def _reset_state() -> None:
    """Clear the held-parameter count before a new light curve."""
    _STATE["held"] = 0


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
        label (str): ``current`` or ``bounded``.

    Returns:
        dict: Adopted shape parameters and how many were held at a bound.
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
        C1=adopted.p.get("C1"),
        A2=adopted.p.get("A2"),
        G2=adopted.p.get("G2"),
        rms_mmag=1e3 * float(adopted.rms),
        n_held=_STATE["held"],
        outside=";".join(_outside(adopted.p)),
    )
    logger.info(
        "%s %s A1=%.4g G1=%.4g D1=%.4g rms=%.2f mmag held=%d outside=%s",
        path.stem, label, row["A1"], row["G1"], row["D1"], row["rms_mmag"],
        row["n_held"], row["outside"] or "none",
    )
    return row


def main() -> None:
    """Fit both curves with the current solver and with the shape bounds."""
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
        rows.append(_one(path, "bounded"))
    destination = OUT / "comparison.csv"
    with open(destination, "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    logger.info("wrote %s", destination)


if __name__ == "__main__":
    main()
