"""Tiered search on every 600 s file of the Tmag 11 first batch.

The search is the one tried on G5 (``examples/g5_c600_tiered.py``):

1. Baseline: profile 1 with the smallest proposed harmonic set, every proposed
   baseline order.
2. Profile: that harmonic set and baseline, profiles 1 to 4.
3. Harmonics: the winning profile and baseline, every other proposed harmonic
   set.

A step stops the search for that file when two steps in a row lower BIC by
less than the decisive difference (10). The trust-region fallback is off.
Profiles 5 to 7 are not tried. No truth column is read and nothing is timed.
Every finished candidate keeps its profile figure and its trend figure.
A file whose first step finds no finished candidate is recorded and the batch
continues.

Run from the repository root::

    MPLBACKEND=Agg python examples/batch_m11_tiered_c600.py
    MPLBACKEND=Agg python examples/batch_m11_tiered_c600.py \\
        --manifest /home/voz/projects/UPJS/EB_synthetic_lc/results/mix/manifest.csv \\
        --out results/synthetic_scoring/batch_mix_tiered_c600

The first command writes ``results/synthetic_scoring/batch_m11_tiered_c600/``.
``--manifest`` and ``--out`` repeat the same search on another degradation.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import logging
import shutil
import sys
import time
from pathlib import Path

import mikulasek_eb as eb
from mikulasek_eb.candidates import Candidate
from mikulasek_eb.errors import MikulasekError
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.plotting import plot_global_fit, plot_instrumental_trend
from mikulasek_eb.selection import DELTA_BIC_DECISIVE

logger = logging.getLogger(__name__)
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFEST = GENERATED / "results" / "m11" / "manifest.csv"
OUT = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "synthetic_scoring"
    / "batch_m11_tiered_c600"
)
PROFILES = (1, 2, 3, 4)
N_EXPECTED = 8
RETRY_BOUND_ACTIVE = False
CANDIDATE_FIELDS = [
    "physics_id",
    "recipe_id",
    "file",
    "step",
    "candidate",
    "status",
    "rms_mmag",
    "bic",
    "delta_bic",
    "step_seconds",
    "step_winner",
]
CHOSEN_FIELDS = [
    "physics_id",
    "recipe_id",
    "file",
    "chosen",
    "status",
    "rms_mmag",
    "step_gains_bic",
    "stopped_early",
    "total_seconds",
    "solution",
    "error",
]


def rows_600(physics: str | None = None) -> list[dict]:
    """Accepted 600 s rows of the manifest.

    Args:
        physics (str | None): Keep only this physics id. None keeps every 600 s row.

    Returns:
        list[dict]: One row per file, ordered by physics id then recipe id.

    Raises:
        RuntimeError: If the number of rows is not the expected one.
    """
    with open(MANIFEST, newline="") as fh:
        rows = [
            r
            for r in csv.DictReader(fh)
            if int(r["cadence_s"]) == 600 and int(r["fail"]) == 0
        ]
    rows.sort(key=lambda r: (r["physics_id"], r["recipe_id"]))
    if physics is not None:
        rows = [r for r in rows if r["physics_id"] == physics]
    expected = 1 if physics is not None else N_EXPECTED
    if len(rows) != expected:
        raise RuntimeError(f"expected {expected} files at 600 s, found {len(rows)}")
    return rows


def make_candidate(base: Candidate, profile: int, harmonics: tuple, order: int):
    """Build one candidate from a proposed one.

    The trust-region fallback follows ``RETRY_BOUND_ACTIVE``.

    Args:
        base (Candidate): A proposed candidate, used for data layer and reasons.
        profile (int): Table 2 profile.
        harmonics (tuple): Harmonic orders.
        order (int): Baseline polynomial degree.

    Returns:
        Candidate: The new candidate, named ``p<profile>_h<harmonics>_b<order>``.
    """
    cfg = dataclasses.replace(
        base.config,
        model=dataclasses.replace(
            base.config.model, profile=profile, harmonics=tuple(harmonics)
        ),
        data=dataclasses.replace(base.config.data, baseline_order=order),
        policy=dataclasses.replace(
            base.config.policy, retry_bound_active=RETRY_BOUND_ACTIVE
        ),
    ).validate()
    name = f"p{profile}_h{''.join(map(str, harmonics))}_b{order}"
    return Candidate(name, cfg, base.reasons + (f"Table 2 profile {profile}",))


def save_figures(stem: str, payload: dict, out: eb.FitOutput, title: str) -> None:
    """Save the profile figure and the fitted-trend figure of one candidate.

    Args:
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
    for folder in ("fits", "trends"):
        (OUT / folder).mkdir(parents=True, exist_ok=True)
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
        save_path=str(OUT / "fits" / f"{stem}.png"),
    )
    plot_instrumental_trend(
        eb.mjd_to_jd(out.t_mjd),
        out.mag,
        yhat,
        template,
        OUT / "trends" / f"{stem}.png",
        title,
    )
    logger.info("figures %s", stem)


def run_step(step, candidates, t_mjd, payload, features, label, stem, saved, rows):
    """Fit one step with ``select_model`` and record every candidate.

    Args:
        step (str): Step label for the tables.
        candidates (list): Candidates of this step, previous winner first.
        t_mjd (ndarray): Times (MJD).
        payload (dict): Light curve payload.
        features (Features): Features of this light curve.
        label (str): Physics id, written into the figure title.
        stem (str): File stem of the light curve.
        saved (set): Candidate names that already have pictures; updated here.
        rows (list): Table rows; extended here.

    Returns:
        tuple: ``SelectionResult`` of the step and its wall time in seconds.

    Raises:
        RuntimeError: If no candidate of the step finished with status ok.
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
        for row in result.table:
            rows.append(
                dict(
                    step=step,
                    candidate=row.name,
                    status=row.status,
                    rms_mmag="",
                    bic="",
                    delta_bic="",
                    step_seconds=f"{seconds:.0f}",
                    step_winner=0,
                )
            )
        raise RuntimeError(f"step {step}: no candidate finished with status ok")
    for row in result.table:
        out = result.fits.get(row.name)
        stats = None if out is None else out.solution.stats
        rows.append(
            dict(
                step=step,
                candidate=row.name,
                status=row.status,
                rms_mmag="" if row.rms_mmag is None else f"{row.rms_mmag:.2f}",
                bic="" if row.bic is None else f"{row.bic:.1f}",
                delta_bic="" if row.delta_bic is None else f"{row.delta_bic:.1f}",
                step_seconds=f"{seconds:.0f}",
                step_winner=int(row.name == result.best.name),
            )
        )
        if out is not None and row.name not in saved:
            mark = "ok" if row.status == "ok" else row.status.split(":")[0]
            title = (
                f"{label}  600 s  {row.name}  step {step}  {mark}  "
                f"rms {1e3 * float(stats['rms_mag']):.2f} mmag"
            )
            save_figures(f"{stem}__{row.name}", payload, out, title)
            saved.add(row.name)
    logger.info(
        "%s step %s winner %s in %.0f s", label, step, result.best.name, seconds
    )
    return result, seconds


def gain_of(result, previous_name: str) -> float:
    """BIC gain of the step winner over the previous winner.

    Args:
        result (SelectionResult): Result of the step.
        previous_name (str): Name of the previous winner, fitted in this step.

    Returns:
        float: ``BIC(previous) - BIC(winner)``, zero when the winner is unchanged.

    Raises:
        RuntimeError: If the previous winner has no BIC in this step.
    """
    by_name = {row.name: row for row in result.table}
    previous = by_name[previous_name]
    if previous.bic is None:
        raise RuntimeError(f"{previous_name} did not finish again in this step")
    best = by_name[result.best.name]
    return float(previous.bic - best.bic)


def tiered(manifest_row: dict) -> tuple[list[dict], dict]:
    """Run the three steps with the early stop on one file.

    Args:
        manifest_row (dict): One accepted 600 s row.

    Returns:
        tuple[list[dict], dict]: Candidate rows, then the summary row.
    """
    curve = GENERATED / manifest_row["path"]
    label = manifest_row["physics_id"]
    payload = eb.load_tess_votable(curve)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    features = eb.extract_features(
        t_mjd, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    proposed = eb.propose_candidates(features)
    sets = sorted(
        {c.config.model.harmonics for c in proposed}, key=lambda h: (len(h), h)
    )
    orders = sorted({c.config.data.baseline_order for c in proposed})
    base = proposed[0]
    reference = sets[0]
    rows: list[dict] = []
    saved: set = set()
    total = 0.0
    gains: list[float] = []

    step1 = [make_candidate(base, 1, reference, o) for o in orders]
    result, sec = run_step(
        "1 baseline",
        step1,
        t_mjd,
        payload,
        features,
        label,
        curve.stem,
        saved,
        rows,
    )
    total += sec
    winner = result.best
    order = winner.config.data.baseline_order
    gains.append(gain_of(result, step1[0].name))
    stop = ""

    step2 = [winner] + [
        make_candidate(base, p, reference, order) for p in PROFILES if p != 1
    ]
    result, sec = run_step(
        "2 profile",
        step2,
        t_mjd,
        payload,
        features,
        label,
        curve.stem,
        saved,
        rows,
    )
    total += sec
    gains.append(gain_of(result, winner.name))
    winner = result.best
    if (
        len(gains) >= 2
        and gains[-1] < DELTA_BIC_DECISIVE
        and gains[-2] < DELTA_BIC_DECISIVE
    ):
        stop = "two steps in a row gained less than the decisive BIC difference"

    if not stop and len(sets) > 1:
        profile = winner.config.model.profile
        step3 = [winner] + [make_candidate(base, profile, h, order) for h in sets[1:]]
        result, sec = run_step(
            "3 harmonics",
            step3,
            t_mjd,
            payload,
            features,
            label,
            curve.stem,
            saved,
            rows,
        )
        total += sec
        gains.append(gain_of(result, winner.name))
        winner = result.best
    final = result.output
    solutions = OUT / "solutions"
    solutions.mkdir(parents=True, exist_ok=True)
    destination = solutions / f"{curve.stem}.json"
    final.solution.save(destination)
    summary = dict(
        physics_id=label,
        recipe_id=manifest_row["recipe_id"],
        file=curve.name,
        chosen=winner.name,
        status=final.solution.stats["status"],
        rms_mmag=f"{1e3 * float(final.solution.stats['rms_mag']):.2f}",
        step_gains_bic=" ".join(f"{g:.1f}" for g in gains),
        stopped_early=stop or "no",
        total_seconds=f"{total:.0f}",
        solution=str(destination.relative_to(OUT)),
        error="",
    )
    for row in rows:
        row["physics_id"] = label
        row["recipe_id"] = manifest_row["recipe_id"]
        row["file"] = curve.name
    return rows, summary


def write_table(path: Path, rows: list[dict], fields: list[str]) -> None:
    """Write a rectangular table.

    Args:
        path (Path): Destination.
        rows (list[dict]): Rows.
        fields (list[str]): Column order.
    """
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fields, extrasaction="ignore", restval="")
        writer.writeheader()
        writer.writerows(rows)


def copy_chosen(summary: dict) -> None:
    """Copy the chosen profile figure and trend figure.

    Args:
        summary (dict): Chosen-fit row.

    Raises:
        FileNotFoundError: If a chosen picture is missing.
    """
    if not summary["chosen"]:
        return
    stem = Path(summary["file"]).stem
    picture = f"{stem}__{summary['chosen']}.png"
    for folder in ("fits", "trends"):
        src = OUT / folder / picture
        if not src.is_file():
            raise FileNotFoundError(src)
        dest = OUT / folder / "chosen"
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest / picture)


def main() -> None:
    """Run the tiered search on every 600 s file and save every picture."""
    global MANIFEST, OUT, RETRY_BOUND_ACTIVE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--physics", default=None)
    parser.add_argument("--retry-bound-active", action="store_true")
    args = parser.parse_args()
    MANIFEST = args.manifest
    OUT = args.out if args.out.is_absolute() else Path.cwd() / args.out
    RETRY_BOUND_ACTIVE = args.retry_bound_active
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if not MANIFEST.is_file():
        raise FileNotFoundError(MANIFEST)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "command.txt").write_text(" ".join(sys.argv) + "\n")
    candidates: list[dict] = []
    chosen_rows: list[dict] = []
    for manifest_row in rows_600(args.physics):
        label = f"{manifest_row['physics_id']} {manifest_row['recipe_id']}"
        try:
            rows, summary = tiered(manifest_row)
        except (MikulasekError, RuntimeError) as ex:
            logger.error("%s failed: %s: %s", label, type(ex).__name__, ex)
            summary = dict(
                physics_id=manifest_row["physics_id"],
                recipe_id=manifest_row["recipe_id"],
                file=Path(manifest_row["path"]).name,
                chosen="",
                status="failed",
                rms_mmag="",
                step_gains_bic="",
                stopped_early="",
                total_seconds="",
                solution="",
                error=f"{type(ex).__name__}: {ex}",
            )
            rows = []
        else:
            copy_chosen(summary)
            logger.info(
                "%s chosen %s rms %s mmag in %s s",
                label,
                summary["chosen"],
                summary["rms_mmag"],
                summary["total_seconds"],
            )
        candidates.extend(rows)
        chosen_rows.append(summary)
        write_table(OUT / "candidates.csv", candidates, CANDIDATE_FIELDS)
        write_table(OUT / "chosen.csv", chosen_rows, CHOSEN_FIELDS)
    logger.info("wrote %s", OUT / "chosen.csv")


if __name__ == "__main__":
    main()
