"""Mikulasek timings of the degraded spot-sweep curves, with two pictures per file.

One model structure is taken from the spot-free 120 s curve (O'Connell term kept on) and fitted
separately to every requested file, so each file has its own baseline and its own parameters.
Each file gets the instrumental-trend picture and the four-panel profile picture. Run::

    python examples/spot_degraded_timing.py --cadence 120 --window 0.12
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import logging
import sys
from pathlib import Path

import numpy as np

from mikulasek_eb import fit, mjd_to_jd, time as time_eclipses
from mikulasek_eb.config import Config, TimingRecipe
from mikulasek_eb.features import extract_features
from mikulasek_eb.candidates import propose_candidates
from mikulasek_eb.helpers import absolute_jd_to_mjd, predicted_magnitude
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.plotting import plot_global_fit, plot_instrumental_trend
from mikulasek_eb.selection import select_model
from mikulasek_eb.synthetic_check import _array, _num, load_synthetic

logger = logging.getLogger(__name__)
DAY = 86400.0
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFEST = GENERATED / "results" / "manifest_spotsweep.csv"
OUT = Path(__file__).resolve().parents[1] / "results" / "synthetic_scoring" / "spot_degraded"


def injected_trend_mag(pure_flux: np.ndarray, trend_flux: np.ndarray) -> np.ndarray:
    """Additive flux trend as a magnitude change at the out-of-eclipse level.

    Using the flux of each point would put a spike into every eclipse, because the same
    flux trend is a large magnitude where the star is faint. The median of the pure flux
    is the level used here.

    Args:
        pure_flux (ndarray): Noise-free flux without the trend.
        trend_flux (ndarray): Additive trend in flux.

    Returns:
        ndarray: Magnitude change at the median pure flux.
    """
    level = float(np.median(pure_flux))
    return -2.5 * np.log10(level + trend_flux) + 2.5 * np.log10(level)


def rows_for(cadence: int) -> list[dict]:
    """Manifest rows of one cadence.

    Args:
        cadence (int): Cadence in seconds.

    Returns:
        list[dict]: Accepted rows.
    """
    with open(MANIFEST, newline="") as fh:
        return [r for r in csv.DictReader(fh) if int(r["cadence_s"]) == cadence and int(r["fail"]) == 0]


def structure_from_spot_free(window: float) -> Config:
    """Model structure chosen on the spot-free 120 s curve, with the O'Connell term on.

    Args:
        window (float): Timing half-width (cycles).

    Returns:
        Config: Structure. The data layer is replaced per file.
    """
    row = next(r for r in rows_for(120) if r["physics_id"] == "SP_none")
    path = GENERATED / row["path"]
    payload = load_tess_votable(path)
    t = absolute_jd_to_mjd(payload["t"])
    res = select_model(
        t, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"], data_id=path.stem
    )
    if res.best is None:
        raise RuntimeError("no model chosen for the spot-free curve")
    cfg = res.output.solution.config
    logger.info("spot-free selection %s; O'Connell forced on", res.best.name)
    return dataclasses.replace(
        cfg,
        model=dataclasses.replace(cfg.model, ocon=True),
        timing=TimingRecipe(window=window),
    ).validate()


def time_file(path: Path, structure: Config, out_dir: Path) -> list[dict]:
    """Fit one file, save the trend and profile pictures, and time both eclipses.

    Args:
        path (Path): Degraded light curve.
        structure (Config): Model structure to fit.
        out_dir (Path): Directory that receives ``trends/`` and ``fits/``.

    Returns:
        list[dict]: One row per eclipse kind.
    """
    payload = load_tess_votable(path)
    synthetic = load_synthetic(path)
    t = absolute_jd_to_mjd(payload["t"])
    features = extract_features(t, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"])
    own = propose_candidates(features)[0].config.data
    data = dataclasses.replace(own, baseline_order=structure.data.baseline_order)
    config = dataclasses.replace(structure, data=data).validate()
    out = fit(t, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"], config, data_id=path.stem)
    adopted = out.pipeline.adopted
    single = out.pipeline.single_component
    model = out.solution.build_model()
    yhat, template, theta = predicted_magnitude(
        out.t_mjd, out.mag, out.sig, model, adopted.p, adopted.ephem,
        out.pipeline.seg, config.data.baseline_order,
    )
    injected = None
    if synthetic.truth_flux is not None and synthetic.truth_trend is not None:
        injected = injected_trend_mag(synthetic.truth_flux, synthetic.truth_trend)
    trends = out_dir / "trends"
    fits = out_dir / "fits"
    trends.mkdir(parents=True, exist_ok=True)
    fits.mkdir(parents=True, exist_ok=True)
    plot_instrumental_trend(
        out.t_mjd, out.mag, yhat, template, trends / f"{path.stem}.png", path.stem, injected_mag=injected
    )
    rms = 1e3 * float(out.solution.stats["rms_mag"])
    plot_global_fit(
        mjd_to_jd(out.t_mjd), out.mag, out.sig, yhat, template, theta, model,
        single.p, adopted.p, out.pipeline.seg,
        f"{path.stem}  rms {rms:.2f} mmag",
        timescale=payload.get("timescale"),
        save_path=str(fits / f"{path.stem}.png"),
    )
    recipe = dataclasses.replace(config.timing, window=config.timing.window)
    timings = time_eclipses(out.t_mjd, out.mag, out.sig, out.solution, recipe)
    meta = synthetic.meta
    origin = float(synthetic.timesys["timeorigin"])
    kinds = [k.strip() for k in str(meta["eclipse_kind"]).split(",")]
    tmin = _array(meta, "eclipse_T_min_pure_days")
    tconj = _array(meta, "eclipse_T_conj_days")
    period = _num(meta, "P_true_days") or payload["period"]
    rows = []
    for kind in ("primary", "secondary"):
        measured = np.asarray(timings.timings[kind]["T"], float)
        sigma = np.asarray(timings.timings[kind]["sigma"], float)
        mask = np.array([k == kind for k in kinds])
        pure = origin + tmin[mask]
        conj = origin + tconj[mask]
        base = dict(
            file=path.name, physics=meta.get("physics_id"), cadence_s=int(_num(meta, "cadence_s")),
            kind=kind, n=0, n_rejected=len(timings.rejected[kind]),
            window=config.timing.window,
        )
        if len(measured) == 0 or len(pure) == 0:
            rows.append(base)
            continue
        nearest = np.argmin(np.abs(measured[:, None] - pure[None, :]), axis=1)
        keep = np.abs(measured - pure[nearest]) < period / 4
        dp = (measured - pure[nearest])[keep]
        dc = (measured - conj[nearest])[keep]
        sig = sigma[keep]
        base.update(
            n=int(len(dp)),
            bias_vs_pure_s=float(np.mean(dp)) * DAY,
            bias_vs_conj_s=float(np.mean(dc)) * DAY,
            scatter_vs_pure_s=float(np.std(dp)) * DAY,
            median_sigma_s=float(np.median(sig)) * DAY,
            honesty_vs_pure=float(np.sqrt(np.mean(((dp - np.mean(dp)) / sig) ** 2))),
        )
        rows.append(base)
        logger.info(
            "%s %s n=%d bias vs conj %+.1f s, vs pure %+.1f s",
            path.stem, kind, base["n"], base["bias_vs_conj_s"], base["bias_vs_pure_s"],
        )
    return rows


def main() -> None:
    """Time the requested cadences and write the table and both pictures per file."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--cadence", type=int, nargs="+", required=True, help="Cadence in seconds, one or more of 120 600 1800.")
    parser.add_argument("--window", type=float, default=0.12, help="Timing half-width in cycles. Default 0.12.")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "command.txt").write_text(" ".join(sys.argv) + "\n")
    structure = structure_from_spot_free(args.window)
    (OUT / "model.txt").write_text(repr(structure.model) + "\n" + repr(structure.timing) + "\n")
    rows = []
    for cadence in args.cadence:
        for row in rows_for(cadence):
            rows.extend(time_file(GENERATED / row["path"], structure, OUT))
    fields = [
        "file", "physics", "cadence_s", "kind", "n", "n_rejected", "window",
        "bias_vs_conj_s", "bias_vs_pure_s", "scatter_vs_pure_s", "median_sigma_s", "honesty_vs_pure",
    ]
    destination = OUT / ("timing_" + "_".join(str(c) for c in args.cadence) + f"_w{args.window}.csv")
    with open(destination, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    logger.info("wrote %s, %s and %s", destination, OUT / "trends", OUT / "fits")


if __name__ == "__main__":
    main()
