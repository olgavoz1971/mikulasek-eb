"""Mikulasek fit and timing of the nine 120 s spot curves, each from its own candidates.

The harmonic orders, the O'Connell term and the baseline are chosen on that file.
The solver, the verdict and the second-component gates are the package defaults.
No extra weight is put on the minima. The timing window is the centre ± 0.12
cycles (full width 0.24), the window of this experiment. Every fitted candidate keeps its profile figure and its
trend figure. Only the chosen candidate is timed, and only when its verdict is ok.

Run from the repository root::

    MPLBACKEND=Agg python examples/spot_120_select.py
    MPLBACKEND=Agg python examples/spot_120_select.py --set m11

``--set n01`` (the default) reads the 1 percent files and writes
``results/synthetic_scoring/spot_120_select/``. ``--set m11`` reads the Tmag 11
files and writes ``results/synthetic_scoring/spot_120_select_m11/``. Each
directory holds ``fits/``, ``trends/`` and ``timing_120_w0.12.csv``.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

import numpy as np
from volightcurve import JD_TO_MJD

import mikulasek_eb as eb
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.plotting import plot_global_fit, plot_instrumental_trend
from mikulasek_eb.synthetic_check import _array, _num, load_synthetic

logger = logging.getLogger(__name__)
DAY = 86400.0
WINDOW = 0.12
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
SCORING = Path(__file__).resolve().parents[1] / "results" / "synthetic_scoring"
SETS = {
    "n01": ("results/n01/manifest_spotsweep.csv", "spot_120_select"),
    "m11": ("results/m11/manifest_spotsweep.csv", "spot_120_select_m11"),
}
MANIFEST = GENERATED / SETS["n01"][0]
OUT = SCORING / SETS["n01"][1]


def rows_120() -> list[dict]:
    """Accepted 120 s rows of the spot-sweep manifest.

    Returns:
        list[dict]: One manifest row per light curve.
    """
    with open(MANIFEST, newline="") as fh:
        return [
            r
            for r in csv.DictReader(fh)
            if int(r["cadence_s"]) == 120 and int(r["fail"]) == 0
        ]


def injected_trend_mag(pure_flux: np.ndarray, trend_flux: np.ndarray) -> np.ndarray:
    """Additive flux trend as a magnitude change at the out-of-eclipse level.

    Args:
        pure_flux (ndarray): Noise-free flux without the trend.
        trend_flux (ndarray): Additive trend in flux.

    Returns:
        ndarray: Magnitude change at the median pure flux.
    """
    level = float(np.median(pure_flux))
    return -2.5 * np.log10(level + trend_flux) + 2.5 * np.log10(level)


def save_figures(
    curve: Path, stem: str, payload: dict, out: eb.FitOutput, title: str
) -> None:
    """Save the profile figure and the trend figure of one candidate.

    Args:
        curve (Path): Light curve, used to read the injected trend.
        stem (str): File name without a suffix.
        payload (dict): Light curve as returned by ``load_tess_votable``.
        out (FitOutput): Fit of this candidate.
        title (str): Figure title.
    """
    adopted = out.pipeline.adopted
    model = out.solution.build_model()
    order = out.solution.config.data.baseline_order
    yhat, template, theta = predicted_magnitude(
        out.t_mjd,
        out.mag,
        out.sig,
        model,
        adopted.p,
        adopted.ephem,
        out.pipeline.seg,
        order,
    )
    fits = OUT / "fits"
    trends = OUT / "trends"
    fits.mkdir(parents=True, exist_ok=True)
    trends.mkdir(parents=True, exist_ok=True)
    plot_global_fit(
        eb.mjd_to_jd(out.t_mjd),
        out.mag,
        out.sig,
        yhat,
        template,
        theta,
        model,
        out.pipeline.single_component.p,
        adopted.p,
        out.pipeline.seg,
        title,
        timescale=payload.get("timescale"),
        save_path=str(fits / f"{stem}.png"),
    )
    synthetic = load_synthetic(curve)
    injected = injected_on_kept_points(synthetic, out.t_mjd)
    plot_instrumental_trend(
        out.t_mjd,
        out.mag,
        yhat,
        template,
        trends / f"{stem}.png",
        title,
        injected_mag=injected,
    )
    logger.info("figures %s", stem)


def injected_on_kept_points(synthetic, t_mjd: np.ndarray) -> np.ndarray | None:
    """Injected trend in magnitudes at the times the fit kept.

    The fit drops out-of-eclipse points beyond 5 sigma, so its time array is shorter
    than the file. The trend is taken from the file rows whose MJD matches.

    Args:
        synthetic: Light curve as returned by ``load_synthetic``.
        t_mjd (ndarray): Times the fit used (MJD).

    Returns:
        ndarray | None: Trend at those times, or None when the file has no truth trend.

    Raises:
        RuntimeError: If a fit time has no matching row in the file.
    """
    if synthetic.truth_flux is None or synthetic.truth_trend is None:
        return None
    injected = injected_trend_mag(synthetic.truth_flux, synthetic.truth_trend)
    if len(injected) == len(t_mjd):
        return injected
    file_mjd = synthetic.time + float(synthetic.timesys["timeorigin"]) - JD_TO_MJD
    if np.any(np.diff(file_mjd) <= 0):
        raise RuntimeError("the file times are not strictly increasing")
    pos = np.searchsorted(file_mjd, t_mjd)
    if np.any(pos >= len(file_mjd)) or not np.allclose(file_mjd[pos], t_mjd, atol=1e-8):
        raise RuntimeError("a fit time has no row in the light curve")
    return injected[pos]


def offsets(curve: Path, out: eb.FitOutput) -> list[dict]:
    """Primary and secondary offsets of one ok solution against the two truths.

    Args:
        curve (Path): Light curve, for its truth metadata.
        out (FitOutput): Chosen fit, already known to be ok.

    Returns:
        list[dict]: One row for the primary and one for the secondary.
    """
    synthetic = load_synthetic(curve)
    meta = synthetic.meta
    recipe = eb.TimingRecipe(window=WINDOW)
    timings = eb.time(out.t_mjd, out.mag, out.sig, out.solution, recipe)
    origin = float(synthetic.timesys["timeorigin"])
    kinds = [k.strip() for k in str(meta["eclipse_kind"]).split(",")]
    tmin = _array(meta, "eclipse_T_min_pure_days")
    tconj = _array(meta, "eclipse_T_conj_days")
    period = _num(meta, "P_true_days")
    rows = []
    for kind in ("primary", "secondary"):
        measured = np.asarray(timings.timings[kind]["T"], float)
        sigma = np.asarray(timings.timings[kind]["sigma"], float)
        mask = np.array([k == kind for k in kinds])
        pure = origin + tmin[mask]
        conj = origin + tconj[mask]
        base = dict(kind=kind, n=0, n_rejected=len(timings.rejected[kind]), reason="")
        if len(measured) == 0 or len(pure) == 0:
            rows.append(base)
            continue
        nearest = np.argmin(np.abs(measured[:, None] - pure[None, :]), axis=1)
        keep = np.abs(measured - pure[nearest]) < period / 4
        dp = (measured - pure[nearest])[keep]
        dc = (measured - conj[nearest])[keep]
        sig = sigma[keep]
        if len(dp) == 0:
            rows.append(base)
            continue
        base.update(
            n=int(len(dp)),
            bias_vs_conj_s=f"{float(np.mean(dc)) * DAY:.2f}",
            bias_vs_pure_s=f"{float(np.mean(dp)) * DAY:.2f}",
            scatter_vs_pure_s=f"{float(np.std(dp)) * DAY:.2f}",
            median_sigma_s=f"{float(np.median(sig)) * DAY:.2f}",
            honesty_vs_pure=f"{float(np.sqrt(np.mean(((dp - np.mean(dp)) / sig) ** 2))):.3f}",
        )
        rows.append(base)
        logger.info(
            "%s %s n=%s vs conj %s s, vs pure %s s",
            curve.stem,
            kind,
            base["n"],
            base["bias_vs_conj_s"],
            base["bias_vs_pure_s"],
        )
    return rows


def blank_offsets(reason: str) -> list[dict]:
    """Two eclipse rows with no time, and the reason in place of a number.

    Args:
        reason (str): Why this curve was not timed.

    Returns:
        list[dict]: Primary and secondary, both empty of offsets.
    """
    return [
        dict(
            kind=kind,
            n=0,
            n_rejected="",
            reason=reason,
            bias_vs_conj_s="",
            bias_vs_pure_s="",
            scatter_vs_pure_s="",
            median_sigma_s="",
            honesty_vs_pure="",
        )
        for kind in ("primary", "secondary")
    ]


def candidate_row(
    curve: Path,
    physics: str,
    name: str,
    chosen: bool,
    out: eb.FitOutput | None,
    status: str,
) -> dict:
    """One line describing a candidate, whether or not it was chosen.

    Args:
        curve (Path): Light curve.
        physics (str): ``physics_id`` from the manifest.
        name (str): Candidate name.
        chosen (bool): True when ``select_model`` chose this candidate.
        out (FitOutput | None): Its fit, or None when the fit raised.
        status (str): Verdict, or the exception text.

    Returns:
        dict: A row of the candidate table.
    """
    row = dict(
        file=curve.name,
        physics=physics,
        candidate=name,
        chosen=int(chosen),
        status=status,
        harmonics="",
        ocon="",
        baseline_order="",
        rms_mmag="",
        two_component="",
        reasons="",
        warnings="",
    )
    if out is None:
        return row
    model = out.solution.config.model
    stats = out.solution.stats
    row.update(
        harmonics=" ".join(str(h) for h in model.harmonics),
        ocon=int(bool(model.ocon)),
        baseline_order=out.solution.config.data.baseline_order,
        rms_mmag=f"{1e3 * float(stats['rms_mag']):.2f}",
        two_component=int(bool(stats["two_component"])),
        reasons="; ".join(stats["status_reasons"]),
        warnings="; ".join(stats["status_warnings"]),
    )
    return row


def identity(
    curve: Path, physics: str, out: eb.FitOutput | None, status: str, reason: str
) -> dict:
    """Fields shared by the two eclipse rows of one curve.

    Args:
        curve (Path): Light curve.
        physics (str): ``physics_id`` from the manifest.
        out (FitOutput | None): Chosen fit, or None.
        status (str): Verdict of the chosen fit, or ``none``.
        reason (str): Why the curve was not timed. Empty when it was.

    Returns:
        dict: Identity columns.
    """
    row = dict(
        file=curve.name,
        physics=physics,
        window=WINDOW,
        status=status,
        reason=reason,
        candidate="",
        harmonics="",
        ocon="",
        baseline_order="",
        rms_mmag="",
        two_component="",
    )
    if out is None:
        return row
    model = out.solution.config.model
    stats = out.solution.stats
    row.update(
        harmonics=" ".join(str(h) for h in model.harmonics),
        ocon=int(bool(model.ocon)),
        baseline_order=out.solution.config.data.baseline_order,
        rms_mmag=f"{1e3 * float(stats['rms_mag']):.2f}",
        two_component=int(bool(stats["two_component"])),
    )
    return row


def one_curve(manifest_row: dict) -> tuple[list[dict], list[dict]]:
    """Select, draw every candidate, and time the chosen fit when it is ok.

    Args:
        manifest_row (dict): One accepted row of the spot-sweep manifest.

    Returns:
        tuple[list[dict], list[dict]]: Candidate rows, then the two timing rows.
    """
    curve = GENERATED / manifest_row["path"]
    physics = manifest_row["physics_id"]
    payload = eb.load_tess_votable(curve)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    result = eb.select_model(
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        data_id=curve.stem,
    )
    chosen_name = None if result.best is None else result.best.name
    logger.info("%s chosen %s", physics, chosen_name if chosen_name else "none")
    cand_rows = []
    for row in result.table:
        out = result.fits.get(row.name)
        cand_rows.append(
            candidate_row(
                curve, physics, row.name, row.name == chosen_name, out, row.status
            )
        )
        if out is None:
            logger.info("%s %s no figure: %s", curve.stem, row.name, row.status)
            continue
        stats = out.solution.stats
        mark = "chosen" if row.name == chosen_name else "not chosen"
        title = (
            f"{physics}  {row.name}  {mark}  {stats['status']}  "
            f"rms {1e3 * float(stats['rms_mag']):.2f} mmag"
        )
        save_figures(curve, f"{curve.stem}__{row.name}", payload, out, title)
    out = result.output
    if out is None:
        reason = "; ".join(result.reasons) if result.reasons else "no candidate chosen"
        timed = blank_offsets(reason)
        shared = identity(curve, physics, None, "none", reason)
    else:
        status = out.solution.stats["status"]
        shared = identity(curve, physics, out, status, "")
        shared["candidate"] = chosen_name
        if status != "ok":
            reason = "; ".join(out.solution.stats["status_reasons"])
            shared["reason"] = reason
            timed = blank_offsets(reason)
            logger.info("%s not timed: %s", physics, reason)
        else:
            timed = offsets(curve, out)
    for row in timed:
        row.update(shared)
    return cand_rows, timed


def write_table(path: Path, rows: list[dict], fields: list[str]) -> None:
    """Write a rectangular table.

    Args:
        path (Path): Destination.
        rows (list[dict]): Rows. Missing keys are written empty.
        fields (list[str]): Column order.
    """
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fields, extrasaction="ignore", restval="")
        writer.writeheader()
        writer.writerows(rows)
    logger.info("wrote %s", path)


def main() -> None:
    """Fit and time the nine 120 s curves, and keep every candidate figure."""
    global MANIFEST, OUT
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--set",
        choices=tuple(SETS),
        default="n01",
        help="n01 is the 1 percent set. m11 is the Tmag 11 set.",
    )
    args = parser.parse_args()
    manifest_rel, out_name = SETS[args.set]
    MANIFEST = GENERATED / manifest_rel
    OUT = SCORING / out_name
    logger.info("noise set %s, manifest %s, output %s", args.set, MANIFEST, OUT)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "command.txt").write_text(" ".join(sys.argv) + "\n")
    manifest = rows_120()
    if len(manifest) != 9:
        raise RuntimeError(f"expected 9 accepted 120 s curves, found {len(manifest)}")
    candidates: list[dict] = []
    timings: list[dict] = []
    for manifest_row in manifest:
        cand_rows, timed = one_curve(manifest_row)
        candidates.extend(cand_rows)
        timings.extend(timed)
    write_table(
        OUT / "candidates_120.csv",
        candidates,
        [
            "file",
            "physics",
            "candidate",
            "chosen",
            "status",
            "harmonics",
            "ocon",
            "baseline_order",
            "rms_mmag",
            "two_component",
            "reasons",
            "warnings",
        ],
    )
    write_table(
        OUT / "timing_120_w0.12.csv",
        timings,
        [
            "file",
            "physics",
            "candidate",
            "harmonics",
            "ocon",
            "baseline_order",
            "rms_mmag",
            "two_component",
            "window",
            "status",
            "reason",
            "kind",
            "n",
            "n_rejected",
            "bias_vs_conj_s",
            "bias_vs_pure_s",
            "scatter_vs_pure_s",
            "median_sigma_s",
            "honesty_vs_pure",
        ],
    )


if __name__ == "__main__":
    main()
