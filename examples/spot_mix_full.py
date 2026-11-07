"""Full comparison on the mixed-trend spot sweep, no timing.

Every accepted file of ``results/mix/manifest_spotsweep.csv`` is fitted from
the degraded light curve only: time, magnitude, magnitude error, period and
epoch. No truth column is read.

The candidates are the sets ``propose_candidates`` would offer, repeated for
Table 2 profiles 1 to 4. Profiles 5 to 7 are not tried. Levenberg-Marquardt
runs first. A candidate marked ``bound-active`` is fitted again with the
trust-region solver, and a candidate that raises ``FitDidNotConverge`` is
fitted again the same way. The reported solution is the best BIC among the
fits whose status is ``ok``.

Run from the repository root::

    MPLBACKEND=Agg python examples/spot_mix_full.py
    MPLBACKEND=Agg python examples/spot_mix_full.py \\
        --manifest /home/voz/projects/UPJS/EB_synthetic_lc/results/mix/manifest.csv \\
        --out results/synthetic_scoring/batch_mix_full --expected 22

The first command writes ``results/synthetic_scoring/spot_mix_full/``.
``--manifest`` repeats the same comparison on another mixed-trend manifest.
The older ``spot_m11_fit/`` and ``batch_m11_fit/`` results are not touched.
Files are fitted from 1800 s down to 120 s, and the tables are rewritten
after each file.

Every candidate that finishes is kept. Its profile figure, its trend figure
and its solution JSON share one stem, ``{light curve}__{candidate}``. The
solution is the fit that entered the comparison, including a successful
trust-region retry. ``selection/{light curve}.json`` is the decision record:
the best name, the comparison table, the reasons, the remarks and ``n_eff``.
``rms_fit_mmag`` is the rms on the points that candidate kept.
``rms_common_mmag`` is the rms on the points common to every candidate with
status ``ok``, which is the rms that entered the BIC. A candidate that raises
has a protocol line and no solution. The light curve is not copied.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import logging
import math
import shutil
import sys
import time
from pathlib import Path

import mikulasek_eb as eb
from mikulasek_eb.candidates import Candidate
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.plotting import plot_global_fit, plot_instrumental_trend

logger = logging.getLogger(__name__)
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFEST = GENERATED / "results" / "mix" / "manifest_spotsweep.csv"
OUT = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "synthetic_scoring"
    / "spot_mix_full"
)
PROFILES = (1, 2, 3, 4)
N_EXPECTED = 27
CANDIDATE_FIELDS = [
    "file",
    "cadence_s",
    "physics",
    "candidate",
    "profile",
    "chosen",
    "status",
    "harmonics",
    "ocon",
    "baseline_order",
    "rms_fit_mmag",
    "rms_common_mmag",
    "bic",
    "delta_bic",
    "k",
    "n_common",
    "two_component",
    "solver",
    "seconds",
    "solution",
]
RETRY_REASONS = (
    "trust-region fallback after a bound-active verdict",
    "trust-region fallback after the fit did not converge",
)
CHOSEN_FIELDS = [
    "file",
    "cadence_s",
    "physics",
    "chosen",
    "status",
    "harmonics",
    "ocon",
    "baseline_order",
    "rms_mmag",
    "solver",
    "seconds",
    "solution",
    "selection_reasons",
]


def manifest_rows() -> list[dict]:
    """Accepted rows of the mixed-trend spot sweep.

    Returns:
        list[dict]: One row per file, 1800 s first, then 600 s, then 120 s.

    Raises:
        RuntimeError: If the number of accepted files is not 27.
    """
    with open(MANIFEST, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if int(r["fail"]) == 0]
    rows.sort(key=lambda r: (-int(r["cadence_s"]), r["physics_id"]))
    if len(rows) != N_EXPECTED:
        raise RuntimeError(f"expected {N_EXPECTED} accepted curves, found {len(rows)}")
    return rows


def full_candidates(features) -> list:
    """Profile 1 to 4 copies of the proposed candidates, fallback on.

    Args:
        features (Features): Features of this light curve.

    Returns:
        list[Candidate]: One candidate per profile, harmonic set and baseline order.
    """
    out = []
    for cand in eb.propose_candidates(features):
        for profile in PROFILES:
            cfg = dataclasses.replace(
                cand.config,
                model=dataclasses.replace(cand.config.model, profile=profile),
                policy=dataclasses.replace(cand.config.policy, retry_bound_active=True),
            ).validate()
            out.append(
                Candidate(
                    f"p{profile}_{cand.name}",
                    cfg,
                    cand.reasons + (f"Table 2 profile {profile}",),
                )
            )
    return out


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
    (OUT / "fits").mkdir(parents=True, exist_ok=True)
    (OUT / "trends").mkdir(parents=True, exist_ok=True)
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


def solver_name(out: eb.FitOutput) -> str:
    """Name the solver that produced a fit.

    Args:
        out (FitOutput): A finished fit.

    Returns:
        str: ``trust-region`` or ``levenberg-marquardt``.
    """
    if out.solution.config.policy.bounded_solver:
        return "trust-region"
    return "levenberg-marquardt"


def _json_number(value):
    """Return a JSON number, or None when the value is missing or not finite.

    Args:
        value: A count, a BIC, or an rms.

    Returns:
        float | int | None: The value, or None.
    """
    if value is None:
        return None
    number = float(value)
    if math.isnan(number) or math.isinf(number):
        return None
    return number


def retry_reason(reasons) -> str | None:
    """Return the trust-region retry note stored on a candidate, if it has one.

    Args:
        reasons (tuple): Proposal reasons, plus a retry note when the stored fit is a retry.

    Returns:
        str | None: The retry note, or None when the stored fit is the first attempt.
    """
    for item in reasons:
        if item in RETRY_REASONS:
            return item
    return None


def save_protocol(stem: str, light_curve: Path, result) -> None:
    """Write the decision record of one light curve.

    Args:
        stem (str): Light-curve stem, shared with the solution files.
        light_curve (Path): Degraded light curve that was fitted. It is recorded, not copied.
        result (SelectionResult): Comparison for this light curve.
    """
    by_name = {cand.name: cand for cand in result.candidates}
    table = []
    for row in result.table:
        cand = by_name[row.name]
        out = result.fits.get(row.name)
        note = retry_reason(cand.reasons)
        fit_rms = None
        solution = None
        solver = None
        if out is not None:
            fit_rms = 1e3 * float(out.solution.stats["rms_mag"])
            solution = f"solutions/{stem}__{row.name}.json"
            solver = solver_name(out)
        table.append(
            dict(
                name=row.name,
                config_hash=row.config_hash,
                status=row.status,
                rms_fit_mmag=_json_number(fit_rms),
                rms_common_mmag=_json_number(row.rms_mmag),
                k=row.k,
                n_common=row.n_common,
                bic=_json_number(row.bic),
                delta_bic=_json_number(row.delta_bic),
                two_component=row.two_component,
                proposal_reasons=list(cand.reasons),
                retried=note is not None,
                retry_reason=note,
                solver=solver,
                solution=solution,
            )
        )
    protocol = dict(
        light_curve=str(light_curve),
        best=None if result.best is None else result.best.name,
        ask_human=result.ask_human,
        reasons=list(result.reasons),
        remarks=list(result.remarks),
        n_eff=_json_number(result.n_eff),
        n_eff_measured=result.n_eff_measured,
        features_hash=result.features.features_hash(),
        notes=dict(
            rms_fit_mmag="rms of this candidate on the points it kept, in mmag",
            rms_common_mmag=(
                "rms on the points common to every candidate with status ok, "
                "in mmag; this is the rms that entered the BIC"
            ),
        ),
        table=table,
    )
    folder = OUT / "selection"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{stem}.json").write_text(
        json.dumps(protocol, indent=2) + "\n", encoding="utf-8"
    )


def one_curve(manifest_row: dict) -> tuple[list[dict], dict]:
    """Fit one spot file with the full comparison.

    Args:
        manifest_row (dict): One accepted row of the spot-sweep manifest.

    Returns:
        tuple[list[dict], dict]: Candidate rows, then the chosen-fit row.
    """
    curve = GENERATED / manifest_row["path"]
    physics = manifest_row["physics_id"]
    cadence_s = manifest_row["cadence_s"]
    payload = eb.load_tess_votable(curve)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    features = eb.extract_features(
        t_mjd, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    started = time.perf_counter()
    result = eb.select_model(
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        features=features,
        candidates=full_candidates(features),
        data_id=curve.stem,
        retry_not_converged=True,
    )
    seconds = f"{time.perf_counter() - started:.0f}"
    chosen_name = None if result.best is None else result.best.name
    logger.info(
        "%s %ss chosen %s in %s s",
        physics,
        cadence_s,
        chosen_name if chosen_name else "none",
        seconds,
    )
    cand_rows = []
    for row in result.table:
        out = result.fits.get(row.name)
        record = dict(
            file=curve.name,
            cadence_s=cadence_s,
            physics=physics,
            candidate=row.name,
            profile=row.name.split("_", 1)[0][1:],
            chosen=int(row.name == chosen_name),
            status=row.status,
            harmonics="",
            ocon="",
            baseline_order="",
            rms_fit_mmag="",
            solver="",
            seconds=seconds,
            solution="",
        )
        if out is not None:
            model = out.solution.config.model
            stats = out.solution.stats
            record.update(
                harmonics=" ".join(str(h) for h in model.harmonics),
                ocon=int(bool(model.ocon)),
                baseline_order=out.solution.config.data.baseline_order,
                rms_fit_mmag=f"{1e3 * float(stats['rms_mag']):.2f}",
                solver=solver_name(out),
            )
            destination = OUT / "solutions" / f"{curve.stem}__{row.name}.json"
            destination.parent.mkdir(parents=True, exist_ok=True)
            out.solution.save(destination)
            record["solution"] = str(destination.relative_to(OUT))
            mark = "chosen" if row.name == chosen_name else "not chosen"
            title = (
                f"{physics}  {cadence_s} s  {row.name}  {mark}  {stats['status']}  "
                f"rms {1e3 * float(stats['rms_mag']):.2f} mmag"
            )
            save_figures(f"{curve.stem}__{row.name}", payload, out, title)
        record["rms_common_mmag"] = (
            "" if row.rms_mmag is None else f"{row.rms_mmag:.2f}"
        )
        record["bic"] = "" if row.bic is None else f"{row.bic:.3f}"
        record["delta_bic"] = "" if row.delta_bic is None else f"{row.delta_bic:.3f}"
        record["k"] = "" if row.k is None else row.k
        record["n_common"] = "" if row.n_common is None else row.n_common
        record["two_component"] = (
            "" if row.two_component is None else int(bool(row.two_component))
        )
        cand_rows.append(record)
    chosen = dict(
        file=curve.name,
        cadence_s=cadence_s,
        physics=physics,
        chosen=chosen_name or "",
        status=(
            "none" if result.output is None else result.output.solution.stats["status"]
        ),
        harmonics="",
        ocon="",
        baseline_order="",
        rms_mmag="",
        solver="",
        seconds=seconds,
        solution="",
        selection_reasons="; ".join(result.reasons),
    )
    save_protocol(curve.stem, curve, result)
    legacy = OUT / "solutions" / f"{curve.stem}.json"
    if legacy.is_file():
        legacy.unlink()
    if result.output is not None:
        model = result.output.solution.config.model
        stats = result.output.solution.stats
        chosen.update(
            harmonics=" ".join(str(h) for h in model.harmonics),
            ocon=int(bool(model.ocon)),
            baseline_order=result.output.solution.config.data.baseline_order,
            rms_mmag=f"{1e3 * float(stats['rms_mag']):.2f}",
            solver=solver_name(result.output),
            solution=f"solutions/{curve.stem}__{chosen_name}.json",
        )
        picture = f"{curve.stem}__{chosen_name}.png"
        for folder in ("fits", "trends"):
            dest = OUT / folder / "chosen"
            dest.mkdir(parents=True, exist_ok=True)
            shutil.copy2(OUT / folder / picture, dest / picture)
    return cand_rows, chosen


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


def main() -> None:
    """Fit a mixed-trend manifest with the full comparison."""
    global MANIFEST, OUT, N_EXPECTED
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--expected", type=int, default=N_EXPECTED)
    args = parser.parse_args()
    MANIFEST = args.manifest
    OUT = args.out if args.out.is_absolute() else Path.cwd() / args.out
    N_EXPECTED = args.expected
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if not MANIFEST.is_file():
        raise FileNotFoundError(MANIFEST)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "command.txt").write_text(" ".join(sys.argv) + "\n")
    candidates: list[dict] = []
    chosen_rows: list[dict] = []
    for manifest_row in manifest_rows():
        cand_rows, chosen = one_curve(manifest_row)
        candidates.extend(cand_rows)
        chosen_rows.append(chosen)
        write_table(OUT / "candidates.csv", candidates, CANDIDATE_FIELDS)
        write_table(OUT / "chosen.csv", chosen_rows, CHOSEN_FIELDS)
    logger.info("wrote %s", OUT / "chosen.csv")


if __name__ == "__main__":
    main()
