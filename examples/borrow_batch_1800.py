"""Borrow the 120 s template for the 1800 s files of the mixed-trend batch.

For each star (G6 is left out) the 1800 s file is described by six solutions:

* ``scratch_old``: the stored from-scratch solution of ``batch_mix_full``, fitted
  without exposure integration;
* ``scratch_exp``: the same configuration refitted from scratch with the exposure
  of the file integrated;
* ``borrow_phase``, ``borrow_amp``, ``borrow_width``, ``borrow_amp_width``: the stored 120 s
  solution of the same star, borrowed with ``fit_borrowed`` and the exposure
  integrated. Profile, Gamma, C, K, harmonics and period are frozen. Free groups:
  zero phase; zero phase and amplitudes; zero phase and widths D; all three.

The baseline and gap settings of each 1800 s file are those of its stored solution.
The exposure is the ``exposure_s`` written in the file (0.8 of the cadence).

Each solution is judged by the residual rms and the residual-against-phase curve,
and by timing the whole file with the relaxed recipe (window 0.25 cycles, 6 points,
wall check off) and scoring the times against the PHOEBE pure minimum and conjunction.
A solution that is not status ``ok`` or a fit that raises keeps a row with its reason.

Writes ``results/synthetic_scoring/borrow_batch_1800/``: ``fits.csv``, ``timing.csv``,
``timing_common.csv``, ``solutions/*.json``, ``residuals_vs_phase.png``, ``timing_scatter.png``, and a profile figure and a
trend figure of every solution in ``fits/`` and ``trends/``.

Run from ``examples``::

    python borrow_batch_1800.py
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import mikulasek_eb as eb
from batch_compare_methods import LABELS, matched, score, truth
from batch_plateau_time import FULL, GENERATED, ROOT, chosen_solution, write_csv
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.plotting import plot_global_fit, plot_instrumental_trend
from mikulasek_eb.solution import ModelSolution
from mikulasek_eb.synthetic_check import _num, load_synthetic

logger = logging.getLogger(__name__)
OUT = ROOT / "results" / "synthetic_scoring" / "borrow_batch_1800"
STARS = ["G1", "G2", "G4", "G5", "G7", "G8"]
CADENCE = 1800
MANIFEST = GENERATED / "results" / "mix" / "manifest.csv"
VARIANTS = {
    "borrow_phase": ("zero_phase",),
    "borrow_amp": ("zero_phase", "amplitudes"),
    "borrow_width": ("zero_phase", "widths"),
    "borrow_amp_width": ("zero_phase", "amplitudes", "widths"),
}
NAMES = ["scratch_old", "scratch_exp", *VARIANTS]
COLOUR = {
    "scratch_old": "black",
    "scratch_exp": "tab:gray",
    "borrow_phase": "tab:blue",
    "borrow_amp": "tab:orange",
    "borrow_width": "tab:green",
    "borrow_amp_width": "tab:red",
}
# What each name fits. The baseline of the 1800 s file is fitted in every one.
FREE_TEXT = {
    "scratch_old": "From scratch on this 1800 s file: profile, shape, depths, widths, harmonics and epoch were all fitted.",
    "scratch_exp": "The same configuration fitted again from scratch. It repeats the stored fit.",
    "borrow_phase": "Borrowed from 120 s. Fitted: the epoch, and the baseline of this file. Frozen: profile, Gamma, C, K, harmonics, period, depths and widths.",
    "borrow_amp": "Borrowed from 120 s. Fitted: the epoch, the eclipse depths, and the baseline. Frozen: profile, Gamma, C, K, harmonics, period and widths.",
    "borrow_width": "Borrowed from 120 s. Fitted: the epoch, the eclipse widths D, and the baseline. Frozen: profile, Gamma, C, K, harmonics, period and depths.",
    "borrow_amp_width": "Borrowed from 120 s. Fitted: the epoch, the depths, the widths D, and the baseline. Frozen: profile, Gamma, C, K, harmonics and period.",
}
RECIPE = eb.TimingRecipe(
    window=0.25, clip_outliers=False, min_points=6, require_contact_walls=False
)
DAY = 86400.0


def files(physics: str) -> dict[int, dict]:
    """Manifest rows of one star at 120 s and 1800 s (second G1 600 s file excluded).

    Args:
        physics (str): Physics id.

    Returns:
        dict[int, dict]: Cadence to manifest row.
    """
    with MANIFEST.open(newline="") as handle:
        rows = [
            r
            for r in csv.DictReader(handle)
            if r["physics_id"] == physics
            and int(r["fail"]) == 0
            and int(r["cadence_s"]) in (120, CADENCE)
        ]
    out = {int(r["cadence_s"]): r for r in rows}
    if set(out) != {120, CADENCE}:
        raise RuntimeError(f"{physics}: need one file at 120 s and at {CADENCE} s")
    return out


def fits_for(row120: dict, row1800: dict, folder: Path) -> dict[str, object]:
    """Fit or borrow every solution of one star at 1800 s.

    Args:
        row120 (dict): Manifest row of the 120 s file (reference).
        row1800 (dict): Manifest row of the 1800 s file.
        folder (Path): Where the solution JSON files are written.

    Returns:
        dict[str, object]: Name to ``FitOutput``, or to the exception raised.
    """
    ref, _ = chosen_solution(GENERATED / row120["path"])
    stored, _ = chosen_solution(GENERATED / row1800["path"])
    curve = GENERATED / row1800["path"]
    payload = eb.load_tess_votable(curve)
    t_mjd = eb.absolute_jd_to_mjd(np.asarray(payload["t"], float))
    args = (t_mjd, payload["mag"], payload["sig"])
    exposure_days = _num(load_synthetic(curve).meta, "exposure_s") / DAY
    data = dataclasses.replace(stored.config.data, exposure_days=exposure_days)
    result: dict[str, object] = {}
    result["scratch_old"] = eb.FitOutput(stored, *args, None)
    try:
        config = eb.Config(
            stored.config.model, data, stored.config.policy, stored.config.timing
        )
        result["scratch_exp"] = eb.fit(
            *args,
            float(payload["period"]),
            float(payload["epoch_jd"]),
            config,
            data_id=curve.stem,
        )
    except Exception as exc:  # recorded, not hidden
        logger.warning("%s scratch_exp: %s: %s", curve.name, type(exc).__name__, exc)
        result["scratch_exp"] = exc
    for name, groups in VARIANTS.items():
        try:
            result[name] = eb.fit_borrowed(
                *args, ref, data, eb.ParameterPlan(groups), data_id=curve.stem
            )
        except Exception as exc:
            logger.warning("%s %s: %s: %s", curve.name, name, type(exc).__name__, exc)
            result[name] = exc
    for name, out in result.items():
        if not isinstance(out, Exception):
            out.solution.save(folder / f"{curve.stem}__{name}.json")
    return result


def fit_row(star: str, name: str, out) -> dict:
    """Summary row of one solution.

    Args:
        star (str): Physics id.
        name (str): Solution name.
        out: ``FitOutput`` or the exception raised.

    Returns:
        dict: rms, status and the main widths and amplitudes.
    """
    row = {"star": star, "solution": name}
    if isinstance(out, Exception):
        row.update(status="raised", reason=f"{type(out).__name__}: {out}")
        return row
    sol = out.solution
    p = sol.params
    row.update(
        status=sol.stats["status"],
        reason="",
        rms_mmag=f"{1e3 * sol.stats['rms_mag']:.3f}",
        n_points=sol.stats["n_points"],
        D1=f"{p['D1']:.6f}",
        D2=f"{p.get('D2', p['D1']):.6f}",
        A1=f"{p['A1']:.4f}",
        A2=f"{p.get('A2', float('nan')):.4f}",
        G1=f"{p['G1']:.3f}",
    )
    return row


def time_all(curve: Path, out) -> dict[str, list[float]]:
    """Time the whole file with a solution and return the times by kind.

    Args:
        curve (Path): 1800 s file.
        out: ``FitOutput``.

    Returns:
        dict[str, list[float]]: Accepted times in absolute JD per kind.
    """
    payload = eb.load_tess_votable(curve)
    t_mjd = eb.absolute_jd_to_mjd(np.asarray(payload["t"], float))
    timed = eb.time(t_mjd, payload["mag"], payload["sig"], out.solution, RECIPE)
    return {
        k: [float(x) for x in timed.timings[k]["T"]] for k in ("primary", "secondary")
    }


def residual_curve(out, bins: int = 100):
    """Binned residual against phase of one solution.

    Args:
        out: ``FitOutput``.
        bins (int): Number of phase bins over one cycle.

    Returns:
        tuple[ndarray, ndarray]: Bin centres (cycles, -0.25 to 0.75) and mean residual (mmag).
    """
    sol = out.solution
    data = sol.config.data
    seg = (
        None
        if data.gap_split_days is None
        else np.cumsum(np.r_[0, np.diff(out.t_mjd) > data.gap_split_days])
    )
    yhat, _, theta = predicted_magnitude(
        out.t_mjd,
        out.mag,
        out.sig,
        sol.build_model(),
        sol.p,
        sol.ephem,
        seg,
        data.baseline_order,
    )
    ph = (theta + 0.25) % 1.0 - 0.25
    edges = np.linspace(-0.25, 0.75, bins + 1)
    idx = np.digitize(ph, edges) - 1
    res = 1e3 * (out.mag - yhat)
    mean = np.array(
        [res[idx == k].mean() if np.any(idx == k) else np.nan for k in range(bins)]
    )
    return 0.5 * (edges[:-1] + edges[1:]), mean


def save_fit_figures(
    star: str, name: str, curve: Path, solution, reference_p: dict
) -> None:
    """Save the profile figure and the trend figure of one stored solution.

    Args:
        star (str): Physics id, for the title.
        name (str): Solution name.
        curve (Path): The 1800 s light curve.
        solution: Stored ``ModelSolution``.
        reference_p (dict): The 120 s template, drawn when it differs from this fit.
    """
    payload = eb.load_tess_votable(curve)
    t_mjd = eb.absolute_jd_to_mjd(np.asarray(payload["t"], float))
    mag = np.asarray(payload["mag"], float)
    sig = np.asarray(payload["sig"], float)
    data = solution.config.data
    seg = (
        None
        if data.gap_split_days is None
        else np.cumsum(np.r_[0, np.diff(t_mjd) > data.gap_split_days])
    )
    model = solution.build_model()
    yhat, template, theta = predicted_magnitude(
        t_mjd, mag, sig, model, solution.p, solution.ephem, seg, data.baseline_order
    )
    title = (
        f"{LABELS[star]} ({star}), {CADENCE} s, {name}, "
        f"rms {1e3 * solution.stats['rms_mag']:.2f} mmag"
    )
    earlier = solution.p if name.startswith("scratch") else reference_p
    (OUT / "fits").mkdir(parents=True, exist_ok=True)
    (OUT / "trends").mkdir(parents=True, exist_ok=True)
    stem = f"{curve.stem}__{name}"
    plot_global_fit(
        eb.mjd_to_jd(t_mjd),
        mag,
        sig,
        yhat,
        template,
        theta,
        model,
        earlier,
        solution.p,
        seg,
        title,
        timescale=payload.get("timescale"),
        save_path=str(OUT / "fits" / f"{stem}.png"),
        note=FREE_TEXT[name],
    )
    plot_instrumental_trend(
        eb.mjd_to_jd(t_mjd),
        mag,
        yhat,
        template,
        OUT / "trends" / f"{stem}.png",
        title,
    )
    logger.info("figures %s", stem)


def save_figures_from_disk() -> None:
    """Draw the profile and trend figures of every solution already saved."""
    for star in STARS:
        rows = files(star)
        curve = GENERATED / rows[CADENCE]["path"]
        reference, _ = chosen_solution(GENERATED / rows[120]["path"])
        for name in NAMES:
            path = OUT / "solutions" / f"{curve.stem}__{name}.json"
            save_fit_figures(star, name, curve, ModelSolution.load(path), reference.p)


def draw_residuals(curves: dict, path: Path) -> None:
    """Residual against phase for every star and solution.

    Args:
        curves (dict): (star, name) to (phase, residual) arrays.
        path (Path): Output PNG.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 7), sharex=True)
    for ax, star in zip(axes.ravel(), STARS):
        for name in NAMES:
            if (star, name) in curves:
                x, y = curves[(star, name)]
                ax.plot(x, y, color=COLOUR[name], lw=1.0, label=name)
        ax.axhline(0, color="0.6", lw=0.6)
        ax.set_title(f"{LABELS[star]} ({star}), {CADENCE} s", fontsize=10)
        ax.grid(alpha=0.3)
    for ax in axes[1]:
        ax.set_xlabel("phase (cycles)")
    for ax in axes[:, 0]:
        ax.set_ylabel("binned residual (mmag)")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle(
        "Residual against phase. Black: stored from-scratch. Gray: from-scratch with exposure. Colours: borrowed from 120 s."
    )
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def draw_timing(rows: list[dict], path: Path) -> None:
    """Scatter of the times about the pure minimum, by star and solution.

    Args:
        rows (list[dict]): Rows of ``timing_common.csv``.
        path (Path): Output PNG.
    """
    shown = [
        ("scratch_exp", "scratch", "tab:purple"),
        ("borrow_phase", "free: epoch", "tab:blue"),
        ("borrow_amp", "free: epoch, depths", "tab:orange"),
        ("borrow_width", "free: epoch, widths", "tab:green"),
        ("borrow_amp_width", "free: epoch, depths, widths", "tab:red"),
    ]
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    width = 0.12
    for ax, kind in zip(axes, ("primary", "secondary")):
        for k, (name, label, colour) in enumerate(shown):
            y = []
            for star in STARS:
                hit = [
                    r
                    for r in rows
                    if (r["star"], r["solution"], r["kind"]) == (star, name, kind)
                ]
                y.append(
                    float(hit[0]["scatter_vs_pure_s"])
                    if hit and hit[0]["scatter_vs_pure_s"] != ""
                    else np.nan
                )
            ax.bar(
                np.arange(len(STARS)) + (k - 2) * width,
                y,
                width,
                color=colour,
                label=label,
            )
        ax.set_yscale("log")
        ax.set_ylabel(f"{kind}: scatter about pure minimum (s)")
        ax.grid(alpha=0.3, axis="y")
    axes[1].set_xticks(range(len(STARS)), [f"{LABELS[s]}\n{s}" for s in STARS])
    axes[0].legend(fontsize=7, ncol=3)
    fig.suptitle(
        f"{CADENCE} s timing, relaxed recipe, eclipses common to the solutions that were timed"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main() -> None:
    """Run every star and write the tables and figures."""
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "solutions").mkdir(exist_ok=True)
    fit_rows, timing_rows, common_rows, curves = [], [], [], {}
    for star in STARS:
        rows = files(star)
        curve = GENERATED / rows[CADENCE]["path"]
        outs = fits_for(rows[120], rows[CADENCE], OUT / "solutions")
        times: dict[str, dict[str, list[float]]] = {}
        for name in NAMES:
            out = outs[name]
            fit_rows.append(fit_row(star, name, out))
            if isinstance(out, Exception):
                continue
            curves[(star, name)] = residual_curve(out)
            try:
                times[name] = time_all(curve, out)
            except Exception as exc:
                logger.warning(
                    "%s %s timing: %s: %s", star, name, type(exc).__name__, exc
                )
                fit_rows[-1]["reason"] = f"timing {type(exc).__name__}: {exc}"
        for kind in ("primary", "secondary"):
            pure, conj, period, _ = truth(curve, kind)
            sets = {n: set(matched(times[n][kind], pure, period)) for n in times}
            common = set.intersection(*sets.values()) if sets else set()
            for name in NAMES:
                for only, bucket in ((None, timing_rows), (common, common_rows)):
                    rec = score(
                        times.get(name, {}).get(kind, []), pure, conj, period, only
                    )
                    rec.update(star=star, solution=name, kind=kind)
                    bucket.append(rec)
    write_csv(OUT / "fits.csv", unify(fit_rows))
    write_csv(OUT / "timing.csv", timing_rows)
    write_csv(OUT / "timing_common.csv", common_rows)
    draw_residuals(curves, OUT / "residuals_vs_phase.png")
    draw_timing(common_rows, OUT / "timing_scatter.png")
    save_figures_from_disk()
    logger.info("wrote %s", OUT)


def unify(rows: list[dict]) -> list[dict]:
    """Give every row the same columns.

    Args:
        rows (list[dict]): Rows with different keys.

    Returns:
        list[dict]: Rows with all keys, empty where missing.
    """
    keys = list(dict.fromkeys(k for r in rows for k in r))
    return [{k: r.get(k, "") for k in keys} for r in rows]


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--figures",
        action="store_true",
        help="Draw the profile and trend figures of the solutions already saved. No fit is run.",
    )
    if parser.parse_args().figures:
        save_figures_from_disk()
    else:
        main()
