"""Tiered Mikulasek search on NSV 807, sector 30, 120 s, SAP flux only.

This follows the protocol of ``docs/talk_outline.md``, section "Protocol of the
tiered fit for a real light curve":

0. SAP flux only, converted with the photcal of the file. The exposure is 0.8
   of the cadence, estimated from the times, not the ``exptime`` of the file.
1. Outlier control, option A: fold by the file period and epoch, take the
   median of 400 phase bins, flag points more than 5 local robust sigma from
   the median of their bin, and remove them once. The removed points are
   written to ``removed.csv``.
2. Candidates proposed by the package. The exposure integration is whatever the
   package proposes (on only if exposure/FWHM > 0.1). The search starts from the
   richest harmonic set and tries the smaller ones in step 3. One baseline polynomial per orbit, split
   at the gap, with one eclipse template for both orbits.
3. Three steps, baseline, profile, harmonics. BIC, decisive difference 10, stop
   after two steps in a row that each gain less. The trust-region retry is off.
   The profile step also tries profile 1 with the curvature C kept free
   (``keep_curvature``), name suffix ``_C``, for flat-bottomed minima.

Nothing is timed. The second profile component is not a candidate here.

Run from ``examples``::

    python nsv807_s30_tiered.py

Writes ``results/nsv807/s30_protocol/``: ``chosen/`` with the final solution and
its two figures, ``candidates/`` with the figures of every candidate.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import logging
import shutil
import time
from pathlib import Path

import numpy as np
from volightcurve import read_lightcurve
from volightcurve.time_reference import (
    normalise_table_epoch_to_absolute_jd,
    resolve_observation_timesys,
    time_offset_to_absolute_jd,
)

import mikulasek_eb as eb
from mikulasek_eb.candidates import Candidate
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.plotting import plot_global_fit, plot_instrumental_trend
from mikulasek_eb.selection import DELTA_BIC_DECISIVE

logger = logging.getLogger(__name__)
DATA = Path("/home/voz/projects/UPJS/lc_shipper/data")
FILES = ("TIC35119266_s30_SPOC_120s.vot",)
OUT = Path(__file__).resolve().parents[1] / "results" / "nsv807" / "s30_protocol"
PROFILES = (1, 2, 3, 4)
START_PROFILE = 1
SKIP_OUTLIERS = False
GAP_DAYS = None
FLUX = "sap_flux"
FLUX_ERR = "sap_flux_err"
OUTLIER_BINS = 400
SIGMA_FROM_BINS = False
OUTLIER_NSIG = 5.0
MIN_POINTS_PER_BIN = 5
MAD_TO_SIGMA = 1.4826


def load_sap(path: Path) -> dict:
    """Load the SAP flux of one file as magnitudes.

    Args:
        path (Path): VOTable.

    Returns:
        dict: ``t`` absolute JD, ``mag``, ``sig``, ``period``, ``epoch_jd``,
        ``timescale``, ``lightcurve``.

    The error column is ``FLUX_ERR``. When it is None (``--no-flux-err``) the file
    has no error for the SAP flux and ``sig`` is a placeholder of 1 mag, to be
    replaced by ``constant_sigma_from_fit``.

    Raises:
        ValueError: If the SAP columns, the period or the epoch are missing,
            or if no finite row remains.
    """
    volc = read_lightcurve(path, filename=path.name, strict=True)
    if FLUX not in volc.get_flux_colnames() or (
        FLUX_ERR is not None and FLUX_ERR not in volc.get_flux_error_colnames()
    ):
        raise ValueError(f"{path.name} has no {FLUX} and {FLUX_ERR}")
    problems = volc.conversion_problems(FLUX)
    if problems:
        raise ValueError(f"{path.name}: {problems}")
    photcal = volc.photdms[FLUX].photcal
    mag_col = volc.add_mag_column_from_flux(FLUX)
    mag = np.asarray(volc.table[mag_col].to("mag").value, dtype=float)
    if FLUX_ERR is None:
        sig = np.ones_like(mag)
    else:
        sig = np.asarray(
            photcal.flux_err_to_mag_err(volc.table[FLUX], volc.table[FLUX_ERR])
            .to("mag")
            .value,
            dtype=float,
        )
    ts = resolve_observation_timesys(volc)
    raw_t = np.asarray(volc.table[volc.get_time_colnames()[0]], dtype=float)
    t_jd = np.array(
        [time_offset_to_absolute_jd(float(v), ts.timeorigin) for v in raw_t]
    )
    ok = np.isfinite(t_jd) & np.isfinite(mag) & np.isfinite(sig) & (sig > 0)
    if not np.any(ok):
        raise ValueError(f"{path.name}: no finite SAP row")
    period = float(volc.table.meta["period"])
    epoch_jd = float(
        normalise_table_epoch_to_absolute_jd(volc, volc.table.meta["epoch"])
    )
    logger.info(
        "%s SAP %d rows, median %.3f mag, median sigma %.2f mmag",
        path.name,
        int(ok.sum()),
        float(np.median(mag[ok])),
        1e3 * float(np.median(sig[ok])),
    )
    return {
        "t": t_jd[ok],
        "mag": mag[ok],
        "sig": sig[ok],
        "period": period,
        "epoch_jd": epoch_jd,
        "timescale": ts.timescale,
        "lightcurve": path.stem,
    }


def fold_outliers(payload: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Flag points that are far from the folded light curve (option A).

    The light curve is folded with the file period. The epoch of the file is not
    used: the zero of the phase is the first time, and any zero gives the same
    folded shape. Each point is
    compared with the median of its phase bin, and the scatter of that bin is
    its own robust sigma (1.4826 times the median absolute deviation).

    Args:
        payload (dict): Light curve from ``load_sap``.

    Returns:
        tuple[ndarray, ndarray, ndarray]: Boolean flag (True for an outlier),
        the residual from the bin median (mag) and the phase of every point.

    Raises:
        ValueError: If a phase bin has fewer than ``MIN_POINTS_PER_BIN`` points
            or no scatter, so that its sigma cannot be measured.
    """
    t, mag = payload["t"], payload["mag"]
    phase = ((t - t.min()) / payload["period"]) % 1.0
    index = np.minimum((phase * OUTLIER_BINS).astype(int), OUTLIER_BINS - 1)
    residual = np.empty_like(mag)
    local = np.empty(OUTLIER_BINS)
    for k in range(OUTLIER_BINS):
        in_bin = index == k
        if in_bin.sum() < MIN_POINTS_PER_BIN:
            raise ValueError(
                f"phase bin {k} has {int(in_bin.sum())} points, fewer than "
                f"{MIN_POINTS_PER_BIN}: its sigma cannot be measured"
            )
        values = mag[in_bin]
        residual[in_bin] = values - np.median(values)
        local[k] = MAD_TO_SIGMA * np.median(
            np.abs(residual[in_bin] - np.median(residual[in_bin]))
        )
    if not np.all(local > 0):
        raise ValueError("a phase bin has zero scatter")
    z = residual / local[index]
    return np.abs(z) > OUTLIER_NSIG, residual, phase


def outlier_figure(
    payload: dict,
    flag: np.ndarray,
    phase: np.ndarray,
    path: Path,
    stem: str,
    rule: str = "outliers 5 sigma from the folded curve",
) -> None:
    """Draw the folded light curve with the flagged points marked.

    Args:
        payload (dict): Light curve from ``load_sap``.
        flag (ndarray): True for a flagged point.
        phase (ndarray): Phase of every point.
        path (Path): PNG file to write.
        stem (str): File stem, for the title.
        rule (str): Name of the rule that flagged the points, for the title.
    """
    import matplotlib.pyplot as plt

    mag = payload["mag"]
    fig, axes = plt.subplots(2, 1, figsize=(10, 11))
    low, high = float(mag[~flag].min()), float(mag[~flag].max())
    pad = 0.15 * (high - low)
    for ax, zoom in zip(axes, (False, True)):
        ax.scatter(phase[~flag], mag[~flag], s=3, color="0.6", label="kept")
        ax.scatter(phase[flag], mag[flag], s=14, color="red", label="flagged")
        ax.set_ylim(high + pad, low - pad) if zoom else ax.invert_yaxis()
        ax.set_xlabel("Phase")
        ax.set_ylabel("Magnitude")
        ax.legend()
    axes[0].set_title(f"{stem}\n{rule}: {int(flag.sum())} of {len(flag)} flagged")
    axes[1].set_title("Same, zoomed to the range of the kept points")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def remove_outliers(payload: dict, stem: str) -> dict:
    """Remove the flagged points once and write them to ``removed.csv``.

    Args:
        payload (dict): Light curve from ``load_sap``.
        stem (str): File stem, for the table name.

    Returns:
        dict: The light curve without the flagged points.
    """
    flag, residual, phase = fold_outliers(payload)
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"{stem}__removed.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["t_jd", "phase", "mag", "residual_mmag"])
        for i in np.where(flag)[0]:
            writer.writerow(
                [
                    f"{payload['t'][i]:.6f}",
                    f"{phase[i]:.5f}",
                    f"{payload['mag'][i]:.5f}",
                    f"{1e3 * residual[i]:.2f}",
                ]
            )
    outlier_figure(payload, flag, phase, OUT / f"{stem}__outliers.png", stem)
    logger.info(
        "%s outliers: removed %d of %d points (%.2f percent), listed in %s",
        stem,
        int(flag.sum()),
        len(flag),
        100.0 * float(flag.mean()),
        path.name,
    )
    kept = {k: v for k, v in payload.items()}
    for key in ("t", "mag", "sig"):
        kept[key] = payload[key][~flag]
    kept["n_removed"] = int(flag.sum())
    kept["n_input"] = int(len(flag))
    return kept


def constant_sigma_from_fit(payload: dict, t_mjd: np.ndarray, stem: str) -> dict:
    """Set one constant sigma from the residuals of a pilot fit.

    Used only for a file without an error column for the SAP flux. The pilot is
    the first candidate proposed by the package. It is fitted with a constant
    sigma equal to the robust scatter of the points about their phase-bin
    medians. The root mean square residual of that fit is the constant sigma of
    the protocol. Nothing else changes.

    Args:
        payload (dict): Light curve after the outlier step, with placeholder sigma.
        t_mjd (ndarray): Times, MJD.
        stem (str): File stem, for the log.

    Returns:
        dict: Copy of the payload with the constant sigma (mag) in ``sig``.

    Raises:
        RuntimeError: If the pilot fit did not finish with status ``ok``.
    """
    _, residual, _ = fold_outliers(payload)
    provisional = MAD_TO_SIGMA * float(
        np.median(np.abs(residual - np.median(residual)))
    )
    pilot = dict(payload)
    pilot["sig"] = np.full_like(payload["mag"], provisional)
    features = eb.extract_features(
        t_mjd, pilot["mag"], pilot["sig"], pilot["period"], pilot["epoch_jd"]
    )
    config = eb.propose_candidates(features)[0].config
    out = eb.fit(
        t_mjd,
        pilot["mag"],
        pilot["sig"],
        pilot["period"],
        pilot["epoch_jd"],
        config,
        data_id=f"{stem}__pilot",
    )
    if out.solution.stats["status"] != "ok":
        raise RuntimeError(f"{stem}: pilot fit status {out.solution.stats['status']}")
    rms = float(out.solution.stats["rms_mag"])
    logger.info(
        "%s constant sigma: provisional %.2f mmag, pilot fit rms %.2f mmag adopted",
        stem,
        1e3 * provisional,
        1e3 * rms,
    )
    result = dict(payload)
    result["sig"] = np.full_like(payload["mag"], rms)
    return result


def constant_sigma_from_bins(payload: dict, stem: str) -> dict:
    """Set one constant sigma from the scatter about the phase-bin medians.

    Used only for a file without an error column for the SAP flux, when the
    researcher chose this instead of the pilot fit. The sigma is the robust
    scatter (1.4826 times the median absolute deviation) of the residuals of the
    points from the median of their phase bin, with ``OUTLIER_BINS`` bins. No model
    is fitted. The result is smaller than the rms of a model fit by the misfit of
    that model.

    Args:
        payload (dict): Light curve with placeholder sigma.
        stem (str): File stem, for the log.

    Returns:
        dict: Copy of the payload with the constant sigma (mag) in ``sig``.

    Raises:
        ValueError: If a phase bin cannot be measured (see ``fold_outliers``).
    """
    _, residual, _ = fold_outliers(payload)
    sigma = MAD_TO_SIGMA * float(np.median(np.abs(residual - np.median(residual))))
    logger.info(
        "%s constant sigma from %d phase bins: %.2f mmag",
        stem,
        OUTLIER_BINS,
        1e3 * sigma,
    )
    result = dict(payload)
    result["sig"] = np.full_like(payload["mag"], sigma)
    return result


def make_candidate(
    base: Candidate,
    profile: int,
    harmonics: tuple,
    order: int,
    keep_curvature: bool = False,
) -> Candidate:
    """One candidate with the trust-region fallback off.

    Args:
        base (Candidate): Proposed candidate, for the data layer and the reasons.
        profile (int): Table 2 profile.
        harmonics (tuple): Harmonic orders.
        order (int): Baseline polynomial degree.
        keep_curvature (bool): Profile 1 only. Keep C free, so that the bottom
            of an eclipse can be flat.

    Returns:
        Candidate: Named ``p<profile>_h<harmonics>_b<order>``, with ``_C`` added
        when the curvature is kept.
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


def figure_title(
    lightcurve: str, candidate: str, status: str, adopted: bool, rms_mmag: float
) -> str:
    """Title of one fit figure, the same fields on every attempt.

    Args:
        lightcurve (str): File stem.
        candidate (str): Candidate name. It starts with the profile number.
        status (str): Fit status, for example ``ok`` or ``bound-active``.
        adopted (bool): Whether this candidate is the final choice of the protocol.
        rms_mmag (float): Residual root mean square (mmag).

    Returns:
        str: The light curve, then the profile, the status, adopted or not, and the rms.
    """
    kept = "adopted" if adopted else "not adopted"
    return f"{lightcurve}\n{candidate}  {status}  {kept}  rms {rms_mmag:.2f} mmag"


def save_figures(
    folder: Path, stem: str, payload: dict, out: eb.FitOutput, title: str
) -> None:
    """Save the profile figure and the trend figure of one candidate.

    Args:
        folder (Path): Folder that gets ``fits`` and ``trends`` subfolders.
        stem (str): File name without a suffix.
        payload (dict): SAP light curve.
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
    (folder / "fits").mkdir(parents=True, exist_ok=True)
    (folder / "trends").mkdir(parents=True, exist_ok=True)
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
        save_path=str(folder / "fits" / f"{stem}.png"),
        note="",
    )
    plot_instrumental_trend(
        eb.mjd_to_jd(out.t_mjd),
        out.mag,
        yhat,
        template,
        folder / "trends" / f"{stem}.png",
        title,
    )


def title_of(stem: str, name: str, out: eb.FitOutput, adopted: bool) -> str:
    """Title of the figure of one fit.

    Args:
        stem (str): Light curve name.
        name (str): Candidate name.
        out (FitOutput): The fit.
        adopted (bool): True only for the final choice of the protocol.

    Returns:
        str: Title with the fields of ``figure_title``.
    """
    return figure_title(
        stem,
        name,
        out.solution.stats["status"],
        adopted,
        1e3 * float(out.solution.stats["rms_mag"]),
    )


def run_step(step, candidates, t_mjd, payload, features, stem, saved, rows):
    """Fit one step and keep a figure of every finished candidate.

    Args:
        step (str): Step label.
        candidates (list): Candidates, previous winner first when there is one.
        t_mjd (ndarray): Times, MJD.
        payload (dict): SAP light curve after the outlier step.
        features: Measured features.
        stem (str): File stem.
        saved (set): Candidate names already drawn.
        rows (list): Table rows, extended here.

    Returns:
        tuple: Selection result and the step time in seconds.

    Raises:
        RuntimeError: If no candidate finished with status ``ok``.
    """
    start = time.perf_counter()
    result = eb.select_model(
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        features=features,
        candidates=candidates,
        data_id=stem,
        retry_not_converged=False,
    )
    seconds = time.perf_counter() - start
    if result.best is None:
        raise RuntimeError(f"{stem} step {step}: no candidate finished with status ok")
    for row in result.table:
        out = result.fits.get(row.name)
        rows.append(
            {
                "file": stem,
                "step": step,
                "candidate": row.name,
                "status": row.status,
                "rms_mmag": "" if row.rms_mmag is None else f"{row.rms_mmag:.2f}",
                "bic": "" if row.bic is None else f"{row.bic:.1f}",
                "delta_bic": "" if row.delta_bic is None else f"{row.delta_bic:.1f}",
                "step_seconds": f"{seconds:.0f}",
                "step_winner": int(row.name == result.best.name),
            }
        )
        if out is not None and row.name not in saved:
            save_figures(
                OUT / "candidates",
                f"{stem}__{row.name}",
                payload,
                out,
                title_of(stem, row.name, out, False),
            )
            out.solution.save(OUT / "candidates" / f"{stem}__{row.name}.json")
            saved.add(row.name)
    logger.info("%s step %s winner %s in %.0f s", stem, step, result.best.name, seconds)
    return result, seconds


def gain_of(result, previous_name: str) -> float:
    """BIC of the previous winner minus BIC of the step winner.

    Args:
        result: Selection result of the step.
        previous_name (str): Previous winner, fitted again in this step.

    Returns:
        float: The gain. Zero when the winner is unchanged.

    Raises:
        RuntimeError: If the previous winner has no BIC.
    """
    by_name = {row.name: row for row in result.table}
    previous = by_name[previous_name]
    if previous.bic is None:
        raise RuntimeError(f"{previous_name} did not finish again in this step")
    return float(previous.bic - by_name[result.best.name].bic)


def one_file(name: str) -> tuple[list[dict], dict]:
    """Run the protocol on one light curve.

    Args:
        name (str): File name in the shipper data directory.

    Returns:
        tuple[list[dict], dict]: Candidate rows and the summary row.
    """
    stem = Path(name).stem
    payload = load_sap(DATA / name)
    if SKIP_OUTLIERS:
        payload["n_removed"] = 0
        payload["n_input"] = int(len(payload["t"]))
        logger.info("%s outlier step skipped, all points kept", stem)
    else:
        payload = remove_outliers(payload, stem)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    if FLUX_ERR is None and SIGMA_FROM_BINS:
        payload = constant_sigma_from_bins(payload, stem)
    elif FLUX_ERR is None:
        payload = constant_sigma_from_fit(payload, t_mjd, stem)
    settings = None if GAP_DAYS is None else eb.FeatureSettings(gap_days=GAP_DAYS)
    features = eb.extract_features(
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        settings,
    )
    sampling = features.sampling
    logger.info(
        "%s cadence %.1f s, exposure %.1f s (%s)",
        name,
        sampling["cadence_s"],
        sampling["exposure_s"],
        sampling["exposure_source"],
    )
    proposed = eb.propose_candidates(features)
    sets = sorted(
        {c.config.model.harmonics for c in proposed}, key=lambda h: (-len(h), h)
    )
    orders = sorted({c.config.data.baseline_order for c in proposed})
    base = proposed[0]
    rows: list[dict] = []
    saved: set = set()
    total = 0.0
    gains: list[float] = []

    step1 = [make_candidate(base, START_PROFILE, sets[0], order) for order in orders]
    result, sec = run_step(
        "1 baseline", step1, t_mjd, payload, features, stem, saved, rows
    )
    total += sec
    winner = result.best
    order = winner.config.data.baseline_order
    gains.append(gain_of(result, step1[0].name))
    stop = ""

    step2 = (
        [winner]
        + [
            make_candidate(base, profile, sets[0], order)
            for profile in PROFILES
            if profile != START_PROFILE
        ]
        + [
            make_candidate(base, 1, sets[0], order, keep_curvature=True)
            for _ in ([0] if START_PROFILE == 1 else [])
        ]
    )
    result, sec = run_step(
        "2 profile", step2, t_mjd, payload, features, stem, saved, rows
    )
    total += sec
    gains.append(gain_of(result, winner.name))
    winner = result.best
    if gains[-1] < DELTA_BIC_DECISIVE and gains[-2] < DELTA_BIC_DECISIVE:
        stop = "two steps gained less than the decisive BIC difference"

    if not stop and len(sets) > 1:
        profile = winner.config.model.profile
        keep = winner.config.policy.keep_curvature
        step3 = [winner] + [
            make_candidate(base, profile, harmonics, order, keep_curvature=keep)
            for harmonics in sets[1:]
        ]
        result, sec = run_step(
            "3 harmonics", step3, t_mjd, payload, features, stem, saved, rows
        )
        total += sec
        gains.append(gain_of(result, winner.name))
        winner = result.best

    final = result.output
    chosen = OUT / "chosen"
    chosen.mkdir(parents=True, exist_ok=True)
    final.solution.save(chosen / f"{stem}__{winner.name}.json")
    save_figures(chosen, stem, payload, final, title_of(stem, winner.name, final, True))
    solution = final.solution
    summary = {
        "file": name,
        "chosen": winner.name,
        "status": solution.stats["status"],
        "rms_mmag": f"{1e3 * float(solution.stats['rms_mag']):.2f}",
        "n_input": payload["n_input"],
        "n_removed": payload["n_removed"],
        "profile": str(solution.config.model.profile),
        "harmonics": "".join(map(str, solution.config.model.harmonics)),
        "baseline_order": str(solution.config.data.baseline_order),
        "baseline_pieces": str(
            1
            + int(
                np.sum(
                    np.diff(final.t_mjd) > (solution.config.data.gap_split_days or 0.5)
                )
            )
        ),
        "exposure_s": f"{sampling['exposure_s']:.1f}",
        "step_gains_bic": " ".join(f"{g:.1f}" for g in gains),
        "stopped_early": stop or "no",
        "total_seconds": f"{total:.0f}",
    }
    return rows, summary


def write_rows(path: Path, rows: list[dict]) -> None:
    """Write a table.

    Args:
        path (Path): Destination.
        rows (list[dict]): Rows with the same keys.
    """
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    """Run the protocol on the light curves and write the tables.

    Without arguments the NSV 807 file is used. ``--file`` (path), ``--out``
    (folder), ``--flux`` and ``--flux-err`` (column names) and ``--no-flux-err`` (constant sigma from a pilot fit) run the same protocol on another light curve.
    """
    global OUT, FILES, FLUX, FLUX_ERR, SKIP_OUTLIERS, GAP_DAYS
    global SIGMA_FROM_BINS, OUTLIER_BINS, START_PROFILE
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--file", help="VOTable path instead of the NSV 807 file")
    parser.add_argument("--out", help="output folder instead of s30_protocol")
    parser.add_argument("--flux", help="flux column name, default sap_flux")
    parser.add_argument("--flux-err", help="error column name, default sap_flux_err")
    parser.add_argument(
        "--skip-outliers",
        action="store_true",
        help="no outlier step: too few points for the phase bins (for example 1800 s)",
    )
    parser.add_argument(
        "--no-flux-err",
        action="store_true",
        help="the file has no error column for the flux: constant sigma from a pilot fit",
    )
    parser.add_argument(
        "--gap-days",
        type=float,
        help="split the baseline at gaps longer than this (days); default is the package's 0.5",
    )
    parser.add_argument(
        "--sigma-from-bins",
        action="store_true",
        help="with --no-flux-err: constant sigma from the scatter about phase-bin medians, no pilot fit",
    )
    parser.add_argument(
        "--bins", type=int, help="number of phase bins (default %d)" % OUTLIER_BINS
    )
    parser.add_argument(
        "--start-profile",
        type=int,
        default=1,
        help="profile of the baseline step (default 1); the profile step still tries the others",
    )
    args = parser.parse_args()
    GAP_DAYS = args.gap_days
    SIGMA_FROM_BINS = args.sigma_from_bins
    if args.bins:
        OUTLIER_BINS = args.bins
    START_PROFILE = args.start_profile
    if args.flux:
        FLUX = args.flux
    if args.flux_err:
        FLUX_ERR = args.flux_err
    SKIP_OUTLIERS = args.skip_outliers
    if args.no_flux_err:
        FLUX_ERR = None
    if args.file:
        FILES = (args.file,)
    if args.out:
        OUT = Path(args.out).resolve()
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    rows, summaries = [], []
    for name in FILES:
        part_rows, summary = one_file(name)
        rows.extend(part_rows)
        summaries.append(summary)
        logger.info(
            "chosen %s %s rms %s mmag", name, summary["chosen"], summary["rms_mmag"]
        )
    write_rows(OUT / "candidates.csv", rows)
    write_rows(OUT / "summary.csv", summaries)
    logger.info("wrote %s", OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
