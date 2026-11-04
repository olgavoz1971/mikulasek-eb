"""G4 at 600 s: how the instrumental trend is handled.

The package removes a per-segment polynomial in magnitudes. A trend that is
additive in flux is not a constant magnitude shift: it shifts a deep minimum
by more than a maximum. This experiment tests that on the observed curve only.
No truth column is read, and nothing is timed. The package is not changed.

Runs, all with profile 3, harmonics (1, 2, 3), the gap split and the
trust-region fallback off:

* ``A_mag_b0``, ``A_mag_b1``: the package as it is, trend in magnitudes,
  constant or straight line per segment.
* ``B_flux_o1``, ``B_flux_o2``: back-fitting. The trend is an additive
  polynomial of the given order in flux, per segment, estimated from the
  residual of the current model. It is subtracted from the observed flux, the
  package then fits the result with a constant per segment, and the loop
  repeats until the rms stops changing.

The back-fitting loop is part of this experiment only. The package has no
flux-domain baseline. If it helps, the fix belongs in the package.

Run from the repository root::

    MPLBACKEND=Agg python examples/g4_c600_trend.py

Writes ``results/synthetic_scoring/g4_c600_trend/``.
"""

from __future__ import annotations

import csv
import dataclasses
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u

import mikulasek_eb as eb
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.plotting import plot_global_fit

logger = logging.getLogger(__name__)
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFEST = GENERATED / "results" / "m11" / "manifest.csv"
OUT = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "synthetic_scoring"
    / "g4_c600_trend"
)
MAX_ITERATIONS = 80
RMS_TOLERANCE_MMAG = 0.01
CORE_HALF_WIDTH = 0.03


def g4_600_row() -> dict:
    """The accepted G4 row at 600 s.

    Returns:
        dict: One manifest row.

    Raises:
        RuntimeError: If that row is missing.
    """
    with open(MANIFEST, newline="") as fh:
        rows = [
            r
            for r in csv.DictReader(fh)
            if r["physics_id"] == "G4"
            and int(r["cadence_s"]) == 600
            and int(r["fail"]) == 0
        ]
    if len(rows) != 1:
        raise RuntimeError(f"expected one G4 file at 600 s, found {len(rows)}")
    return rows[0]


def to_flux(mag: np.ndarray) -> np.ndarray:
    """Relative flux of magnitudes.

    Args:
        mag (ndarray): Magnitudes.

    Returns:
        ndarray: Flux, ``10**(-0.4 mag)``, dimensionless.
    """
    return np.asarray(u.Magnitude(mag).physical.to_value(u.dimensionless_unscaled))


def to_mag(flux: np.ndarray) -> np.ndarray:
    """Magnitudes of a relative flux.

    Args:
        flux (ndarray): Positive relative flux.

    Returns:
        ndarray: Magnitudes.
    """
    return np.asarray(u.Magnitude(flux * u.dimensionless_unscaled).value)


def make_config(base, order: int):
    """Profile 3, harmonics (1, 2, 3), the given baseline order, fallback off.

    Args:
        base (Candidate): A proposed candidate for this light curve.
        order (int): Baseline polynomial degree in magnitudes.

    Returns:
        Config: Validated configuration.
    """
    return dataclasses.replace(
        base.config,
        model=dataclasses.replace(base.config.model, profile=3, harmonics=(1, 2, 3)),
        data=dataclasses.replace(base.config.data, baseline_order=order),
        policy=dataclasses.replace(base.config.policy, retry_bound_active=False),
    ).validate()


def model_on_kept(out: eb.FitOutput):
    """Package prediction and fitted phase at the points the fit kept.

    Args:
        out (FitOutput): A finished fit.

    Returns:
        tuple: ``(yhat, template, theta)``.
    """
    adopted = out.pipeline.adopted
    return predicted_magnitude(
        out.t_mjd,
        out.mag,
        out.sig,
        out.solution.build_model(),
        adopted.p,
        adopted.ephem,
        out.pipeline.seg,
        out.solution.config.data.baseline_order,
    )


def kept_index(t_all: np.ndarray, t_kept: np.ndarray) -> np.ndarray:
    """Indices of the kept times in the full time array.

    Args:
        t_all (ndarray): All times, increasing.
        t_kept (ndarray): Times the fit kept.

    Returns:
        ndarray: Index of each kept time.

    Raises:
        RuntimeError: If a kept time is not in the full array.
    """
    index = np.searchsorted(t_all, t_kept)
    if index.max() >= len(t_all) or not np.allclose(t_all[index], t_kept, atol=1e-9):
        raise RuntimeError("kept times do not match the light curve")
    return index


def phase_class(theta: np.ndarray) -> np.ndarray:
    """Label points as primary core, secondary core or out of the cores.

    Args:
        theta (ndarray): Fitted phase in cycles.

    Returns:
        ndarray: Strings ``primary``, ``secondary`` or ``outside``.
    """
    wrapped = np.mod(theta + 0.5, 1.0) - 0.5
    label = np.full(theta.shape, "outside", dtype=object)
    label[np.abs(np.abs(wrapped) - 0.5) < CORE_HALF_WIDTH] = "secondary"
    label[np.abs(wrapped) < CORE_HALF_WIDTH] = "primary"
    return label


def flux_trend(t, seg, resid_flux, weight, order):
    """Per-segment additive polynomial fitted to a flux residual.

    Args:
        t (ndarray): Times at the kept points.
        seg (ndarray): Segment id at the kept points.
        resid_flux (ndarray): Observed minus model flux at the kept points.
        weight (ndarray): Least-squares weights, ``1/sigma``.
        order (int): Polynomial degree.

    Returns:
        dict: Segment id to ``numpy.polynomial.Polynomial``.
    """
    polys = {}
    for s in np.unique(seg):
        sel = seg == s
        polys[int(s)] = np.polynomial.Polynomial.fit(
            t[sel], resid_flux[sel], order, w=weight[sel]
        )
    return polys


def evaluate_trend(polys: dict, t: np.ndarray, seg_all: np.ndarray) -> np.ndarray:
    """Evaluate the per-segment flux trend at every time.

    Args:
        polys (dict): Segment id to polynomial.
        t (ndarray): All times.
        seg_all (ndarray): Segment id of every time.

    Returns:
        ndarray: Trend in relative flux.
    """
    trend = np.zeros_like(t)
    for s, poly in polys.items():
        sel = seg_all == s
        trend[sel] = poly(t[sel])
    return trend


def run_magnitude(label, order, t_mjd, payload, cfg, curve_stem):
    """Package fit with the trend in magnitudes.

    Args:
        label (str): Run label.
        order (int): Baseline degree in magnitudes.
        t_mjd (ndarray): Times (MJD).
        payload (dict): Light curve payload.
        cfg (Config): Configuration.
        curve_stem (str): File stem.

    Returns:
        dict: Fit, kept-point residuals and trend description for plotting.
    """
    out = eb.fit(
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        cfg,
        data_id=f"{curve_stem}_{label}",
    )
    yhat, template, theta = model_on_kept(out)
    resid = out.mag - yhat
    return dict(
        label=label,
        out=out,
        t=out.t_mjd,
        resid=resid,
        theta=theta,
        trend_flux=None,
        iterations=1,
        yhat=yhat,
        template=template,
    )


def run_flux(label, flux_order, t_mjd, payload, cfg, gap, curve_stem):
    """Back-fit an additive flux trend and the package template.

    Args:
        label (str): Run label.
        flux_order (int): Polynomial degree of the additive flux trend.
        t_mjd (ndarray): Times (MJD).
        payload (dict): Light curve payload.
        cfg (Config): Configuration with baseline order 0.
        gap (float | None): Gap that splits the segments, days.
        curve_stem (str): File stem.

    Returns:
        dict: Fit, kept-point residuals, trend and iteration history.
    """
    mag = np.asarray(payload["mag"], float)
    sig = np.asarray(payload["sig"], float)
    flux_obs = to_flux(mag)
    scale = float(np.median(flux_obs))
    flux_rel = flux_obs / scale
    sig_flux = flux_rel * 0.4 * np.log(10.0) * sig
    seg_all = (
        np.zeros(len(t_mjd), int)
        if gap is None
        else np.cumsum(np.r_[0, np.diff(t_mjd) > gap])
    )
    trend = np.zeros_like(flux_rel)
    history = []
    previous = np.inf
    out = None
    for iteration in range(1, MAX_ITERATIONS + 1):
        mag_d = to_mag((flux_rel - trend) * scale)
        out = eb.fit(
            t_mjd,
            mag_d,
            sig,
            payload["period"],
            payload["epoch_jd"],
            cfg,
            data_id=f"{curve_stem}_{label}",
        )
        yhat, template, theta = model_on_kept(out)
        index = kept_index(t_mjd, out.t_mjd)
        model_flux = to_flux(yhat) / scale
        polys = flux_trend(
            out.t_mjd,
            seg_all[index],
            flux_rel[index] - model_flux,
            1.0 / sig_flux[index],
            flux_order,
        )
        trend = evaluate_trend(polys, t_mjd, seg_all)
        full_flux = model_flux + trend[index]
        resid = mag[index] - to_mag(full_flux * scale)
        rms = 1e3 * float(np.sqrt(np.mean(resid**2)))
        history.append((iteration, rms, len(out.t_mjd)))
        logger.info("%s iteration %d rms %.3f mmag", label, iteration, rms)
        if abs(previous - rms) < RMS_TOLERANCE_MMAG:
            break
        previous = rms
    else:
        raise RuntimeError(f"{label}: no convergence in {MAX_ITERATIONS} iterations")
    return dict(
        label=label,
        out=out,
        t=out.t_mjd,
        resid=resid,
        theta=theta,
        trend_flux=trend,
        trend_index=index,
        iterations=iteration,
        yhat=yhat,
        template=template,
        history=history,
        t_all=t_mjd,
        scale=scale,
    )


def plot_run(run: dict, payload: dict, title: str) -> None:
    """Save the profile figure and the residual-versus-time figure of one run.

    Args:
        run (dict): Result of ``run_magnitude`` or ``run_flux``.
        payload (dict): Light curve payload.
        title (str): Figure title.
    """
    out = run["out"]
    adopted = out.pipeline.adopted
    for folder in ("fits", "trends"):
        (OUT / folder).mkdir(parents=True, exist_ok=True)
    plot_global_fit(
        eb.mjd_to_jd(out.t_mjd),
        out.mag,
        out.sig,
        run["yhat"],
        run["template"],
        run["theta"],
        out.solution.build_model(),
        out.pipeline.single_component.p,
        adopted.p,
        out.pipeline.seg,
        title,
        timescale=payload.get("timescale"),
        save_path=str(OUT / "fits" / f"{run['label']}.png"),
    )
    classes = phase_class(run["theta"])
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    colours = {"primary": "tab:red", "secondary": "tab:green", "outside": "0.6"}
    for name, colour in colours.items():
        sel = classes == name
        axes[0].plot(
            run["t"][sel], 1e3 * run["resid"][sel], ".", ms=3, color=colour, label=name
        )
    axes[0].set_ylabel("residual (mmag)")
    axes[0].set_title(title)
    axes[0].legend(loc="upper right")
    if run["trend_flux"] is None:
        axes[1].plot(
            out.t_mjd,
            run["yhat"] - run["template"],
            "b-",
            label="trend in magnitudes",
        )
        axes[1].set_ylabel("trend (mag)")
    else:
        axes[1].plot(
            run["t_all"],
            100 * run["trend_flux"],
            "b-",
            label="additive flux trend",
        )
        axes[1].set_ylabel("trend (% of median flux)")
    axes[1].legend(loc="upper right")
    primary = classes == "primary"
    axes[2].plot(
        run["t"][primary],
        run["template"][primary] + out.mag[primary] - run["yhat"][primary],
        ".",
        ms=3,
        color="tab:red",
        label="primary core, data minus trend",
    )
    axes[2].invert_yaxis()
    axes[2].set_ylabel("magnitude")
    axes[2].set_xlabel("time (MJD)")
    axes[2].legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "trends" / f"{run['label']}.png", dpi=100)
    plt.close(fig)


def summarise(run: dict) -> dict:
    """One table row for a run.

    Args:
        run (dict): Result of ``run_magnitude`` or ``run_flux``.

    Returns:
        dict: Label, rms overall and by phase class, point count, iterations.
    """
    classes = phase_class(run["theta"])
    resid = 1e3 * run["resid"]
    row = dict(
        run=run["label"],
        n_points=len(resid),
        iterations=run["iterations"],
        rms_mmag=f"{np.sqrt(np.mean(resid**2)):.2f}",
    )
    for name in ("primary", "secondary", "outside"):
        sel = classes == name
        row[f"rms_{name}_mmag"] = f"{np.sqrt(np.mean(resid[sel] ** 2)):.2f}"
    return row


def main() -> None:
    """Run the four trend treatments on G4 at 600 s."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "command.txt").write_text(" ".join(sys.argv) + "\n")
    curve = GENERATED / g4_600_row()["path"]
    payload = eb.load_tess_votable(curve)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    features = eb.extract_features(
        t_mjd, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    base = eb.propose_candidates(features)[0]
    gap = base.config.data.gap_split_days
    logger.info("gap split %s days", gap)
    runs = [
        run_magnitude("A_mag_b0", 0, t_mjd, payload, make_config(base, 0), curve.stem),
        run_magnitude("A_mag_b1", 1, t_mjd, payload, make_config(base, 1), curve.stem),
        run_flux("B_flux_o1", 1, t_mjd, payload, make_config(base, 0), gap, curve.stem),
        run_flux("B_flux_o2", 2, t_mjd, payload, make_config(base, 0), gap, curve.stem),
    ]
    rows = []
    for run in runs:
        rows.append(summarise(run))
        plot_run(
            run, payload, f"G4  600 s  {run['label']}  rms {rows[-1]['rms_mmag']} mmag"
        )
    fields = [
        "run",
        "n_points",
        "iterations",
        "rms_mmag",
        "rms_primary_mmag",
        "rms_secondary_mmag",
        "rms_outside_mmag",
    ]
    with open(OUT / "runs.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fields)
        writer.writeheader()
        writer.writerows(rows)
    logger.info("wrote %s", OUT / "runs.csv")


if __name__ == "__main__":
    main()
