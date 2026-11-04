"""G5 at 600 s: tiered candidate search with an early stop, fallback off.

The exhaustive run (``g5_c600_profiles.py``) fits every profile with every
harmonic set and baseline order. This run fits in steps and keeps the winner
of each step, so far fewer fits are needed:

1. Baseline: profile 1 with the smallest proposed harmonic set, every proposed
   baseline order.
2. Profile: that harmonic set and baseline, profiles 1 to 4.
3. Harmonics: the winning profile and baseline, every other proposed harmonic
   set.

Inside each step ``select_model`` compares the previous winner with the new
candidates, so the BIC values of one step share the same points. A step gains
when the best new candidate lowers BIC by at least the decisive difference of
``select_model`` (10). The search stops after two steps in a row that gain
less. The trust-region fallback is off. A step with no finished candidate stops
the run with an error. No truth column is read and nothing is timed.

Run from the repository root::

    MPLBACKEND=Agg python examples/g5_c600_tiered.py

Writes ``results/synthetic_scoring/g5_c600_tiered/``.
"""

from __future__ import annotations

import csv
import dataclasses
import logging
import shutil
import sys
import time
from pathlib import Path

import mikulasek_eb as eb
from mikulasek_eb.candidates import Candidate
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
    / "g5_c600_tiered"
)
PROFILES = (1, 2, 3, 4)


def g5_600_row() -> dict:
    """The accepted G5 row at 600 s.

    Returns:
        dict: One manifest row.

    Raises:
        RuntimeError: If that row is missing.
    """
    with open(MANIFEST, newline="") as fh:
        rows = [
            r
            for r in csv.DictReader(fh)
            if r["physics_id"] == "G5"
            and int(r["cadence_s"]) == 600
            and int(r["fail"]) == 0
        ]
    if len(rows) != 1:
        raise RuntimeError(f"expected one G5 file at 600 s, found {len(rows)}")
    return rows[0]


def make_candidate(base: Candidate, profile: int, harmonics: tuple, order: int):
    """Build one candidate from a proposed one, fallback off.

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
        policy=dataclasses.replace(base.config.policy, retry_bound_active=False),
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
    plot_instrumental_trend(
        eb.mjd_to_jd(out.t_mjd),
        out.mag,
        yhat,
        template,
        trends / f"{stem}.png",
        title,
    )
    logger.info("figures %s", stem)


def run_step(
    step: str,
    candidates: list,
    t_mjd,
    payload: dict,
    features,
    stem: str,
    saved: set,
    rows: list,
):
    """Fit one step with ``select_model`` and record every candidate.

    Args:
        step (str): Step label for the tables.
        candidates (list): Candidates of this step, previous winner first.
        t_mjd (ndarray): Times (MJD).
        payload (dict): Light curve payload.
        features (Features): Features of this light curve.
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
                f"G5  600 s  {row.name}  step {step}  {mark}  "
                f"rms {1e3 * float(stats['rms_mag']):.2f} mmag"
            )
            save_figures(f"{stem}__{row.name}", payload, out, title)
            saved.add(row.name)
    logger.info("step %s winner %s in %.0f s", step, result.best.name, seconds)
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
    """Run the three steps with the early stop.

    Args:
        manifest_row (dict): The accepted G5 row at 600 s.

    Returns:
        tuple[list[dict], dict]: Candidate rows, then the summary row.
    """
    curve = GENERATED / manifest_row["path"]
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
        "1 baseline", step1, t_mjd, payload, features, curve.stem, saved, rows
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
        "2 profile", step2, t_mjd, payload, features, curve.stem, saved, rows
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
            "3 harmonics", step3, t_mjd, payload, features, curve.stem, saved, rows
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
        file=curve.name,
        chosen=winner.name,
        rms_mmag=f"{1e3 * float(final.solution.stats['rms_mag']):.2f}",
        step_gains_bic=" ".join(f"{g:.1f}" for g in gains),
        stopped_early=stop or "no",
        total_seconds=f"{total:.0f}",
        solution=str(destination.relative_to(OUT)),
    )
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
    logger.info("wrote %s", path)


def main() -> None:
    """Run the tiered search on G5 at 600 s and save every picture."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if not MANIFEST.is_file():
        raise FileNotFoundError(MANIFEST)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "command.txt").write_text(" ".join(sys.argv) + "\n")
    rows, summary = tiered(g5_600_row())
    write_table(
        OUT / "candidates.csv",
        rows,
        [
            "step",
            "candidate",
            "status",
            "rms_mmag",
            "bic",
            "delta_bic",
            "step_seconds",
            "step_winner",
        ],
    )
    write_table(
        OUT / "chosen.csv",
        [summary],
        [
            "file",
            "chosen",
            "rms_mmag",
            "step_gains_bic",
            "stopped_early",
            "total_seconds",
            "solution",
        ],
    )
    stem = Path(summary["file"]).stem
    picture = f"{stem}__{summary['chosen']}.png"
    for folder in ("fits", "trends"):
        dest = OUT / folder / "chosen"
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copy2(OUT / folder / picture, dest / picture)


if __name__ == "__main__":
    main()
