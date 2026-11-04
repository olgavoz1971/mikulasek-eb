"""From-scratch Mikulasek fits of the Tmag 11 first batch, no timing.

The manifest is ``results/m11/manifest.csv``. Every accepted file is fitted,
including G6, which carries the spot LC_Generator put on it. The only inputs
are the degraded light curve: time, magnitude, magnitude error, and the
period and epoch stored in that file. No truth column is read, and nothing
is timed.

The chosen solution of each file is written with ``ModelSolution.save``, so
a later timing or borrowing run can load it and does not refit. A file with
no chosen solution saves nothing. The profile figure of each chosen
candidate is copied to ``fits/chosen/``.

Run from the repository root::

    MPLBACKEND=Agg python examples/batch_m11_fit.py
    MPLBACKEND=Agg python examples/batch_m11_fit.py --physics G6

``--physics`` refits only that star and merges it into the existing tables.
Writes ``results/synthetic_scoring/batch_m11_fit/``.
"""

from __future__ import annotations

import argparse
import csv
import logging
import shutil
import sys
from pathlib import Path

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
    / "batch_m11_fit"
)
N_EXPECTED = 22


def manifest_rows(physics: str | None = None) -> list[dict]:
    """Accepted rows of the Tmag 11 first-batch manifest.

    G6 is included. Its spot is the one LC_Generator generated.

    Args:
        physics (str | None): Keep this ``physics_id`` only. None keeps all 22.

    Returns:
        list[dict]: One row per light curve, ordered by physics id then cadence.

    Raises:
        RuntimeError: If a file failed generation, the manifest is not 22
            accepted files, or ``physics`` matches nothing.
    """
    with open(MANIFEST, newline="") as fh:
        rows = list(csv.DictReader(fh))
    failed = [r for r in rows if int(r["fail"]) != 0]
    if failed:
        names = ", ".join(r["path"] for r in failed)
        raise RuntimeError(f"files failed generation: {names}")
    if len(rows) != N_EXPECTED:
        raise RuntimeError(f"expected {N_EXPECTED} accepted files, found {len(rows)}")
    if physics is not None:
        rows = [r for r in rows if r["physics_id"] == physics]
        if not rows:
            raise RuntimeError(f"no accepted files for physics {physics}")
    rows.sort(key=lambda r: (r["physics_id"], int(r["cadence_s"])))
    return rows


def save_figure(stem: str, payload: dict, out: eb.FitOutput, title: str) -> None:
    """Save the profile figure of one candidate.

    The figure is drawn from the fit and from the magnitudes the fit used.

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
    fits.mkdir(parents=True, exist_ok=True)
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
    logger.info("figure %s", stem)


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


def candidate_row(
    curve: Path,
    manifest_row: dict,
    name: str,
    chosen: bool,
    out: eb.FitOutput | None,
    status: str,
) -> dict:
    """One line describing a candidate, whether or not it was chosen.

    Args:
        curve (Path): Light curve.
        manifest_row (dict): Its manifest row, used as labels only.
        name (str): Candidate name.
        chosen (bool): True when ``select_model`` chose this candidate.
        out (FitOutput | None): Its fit, or None when the fit raised.
        status (str): Verdict, or the exception text.

    Returns:
        dict: A row of the candidate table.
    """
    row = dict(
        file=curve.name,
        cadence_s=manifest_row["cadence_s"],
        physics=manifest_row["physics_id"],
        eb_class=manifest_row["eb_class"],
        recipe=manifest_row["recipe_id"],
        candidate=name,
        chosen=int(chosen),
        status=status,
        harmonics="",
        ocon="",
        baseline_order="",
        rms_mmag="",
        two_component="",
        solver="",
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
        solver=solver_name(out),
        reasons="; ".join(stats["status_reasons"]),
        warnings="; ".join(stats["status_warnings"]),
    )
    return row


def one_curve(manifest_row: dict) -> tuple[list[dict], dict]:
    """Select a model for one degraded light curve and save the chosen solution.

    Args:
        manifest_row (dict): One accepted non-spot row of the first-batch manifest.

    Returns:
        tuple[list[dict], dict]: Candidate rows, then the chosen-fit row.
    """
    curve = GENERATED / manifest_row["path"]
    physics = manifest_row["physics_id"]
    cadence_s = manifest_row["cadence_s"]
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
    logger.info(
        "%s %ss chosen %s", physics, cadence_s, chosen_name if chosen_name else "none"
    )
    cand_rows = []
    for row in result.table:
        out = result.fits.get(row.name)
        cand_rows.append(
            candidate_row(
                curve, manifest_row, row.name, row.name == chosen_name, out, row.status
            )
        )
        if out is None:
            logger.info("%s %s no figure: %s", curve.stem, row.name, row.status)
            continue
        stats = out.solution.stats
        mark = "chosen" if row.name == chosen_name else "not chosen"
        title = (
            f"{physics}  {cadence_s} s  {row.name}  {mark}  {stats['status']}  "
            f"rms {1e3 * float(stats['rms_mag']):.2f} mmag"
        )
        save_figure(f"{curve.stem}__{row.name}", payload, out, title)
    solution_path = ""
    if result.output is not None:
        solutions = OUT / "solutions"
        solutions.mkdir(parents=True, exist_ok=True)
        destination = solutions / f"{curve.stem}.json"
        result.output.solution.save(destination)
        solution_path = str(destination.relative_to(OUT))
        logger.info("solution %s", solution_path)
    chosen = dict(
        file=curve.name,
        cadence_s=cadence_s,
        physics=physics,
        eb_class=manifest_row["eb_class"],
        recipe=manifest_row["recipe_id"],
        n_candidates=len(result.table),
        chosen=chosen_name or "",
        status=(
            "none" if result.output is None else result.output.solution.stats["status"]
        ),
        harmonics="",
        ocon="",
        baseline_order="",
        rms_mmag="",
        two_component="",
        solver="",
        solution=solution_path,
        reasons="",
        warnings="",
        selection_reasons="; ".join(result.reasons),
    )
    if result.output is not None:
        model = result.output.solution.config.model
        stats = result.output.solution.stats
        chosen.update(
            harmonics=" ".join(str(h) for h in model.harmonics),
            ocon=int(bool(model.ocon)),
            baseline_order=result.output.solution.config.data.baseline_order,
            rms_mmag=f"{1e3 * float(stats['rms_mag']):.2f}",
            two_component=int(bool(stats["two_component"])),
            solver=solver_name(result.output),
            reasons="; ".join(stats["status_reasons"]),
            warnings="; ".join(stats["status_warnings"]),
        )
    return cand_rows, chosen


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


def copy_chosen_figures(rows: list[dict]) -> None:
    """Copy the profile figure of each chosen candidate into ``fits/chosen``.

    Args:
        rows (list[dict]): Chosen-fit rows. Each must name a chosen candidate.

    Raises:
        RuntimeError: If a row has no chosen candidate.
        FileNotFoundError: If that candidate's figure is missing.
    """
    dest = OUT / "fits" / "chosen"
    dest.mkdir(parents=True, exist_ok=True)
    for row in rows:
        if not row["chosen"]:
            raise RuntimeError(f"{row['file']} has no chosen candidate to copy")
        stem = Path(row["file"]).stem
        src = OUT / "fits" / f"{stem}__{row['chosen']}.png"
        if not src.is_file():
            raise FileNotFoundError(src)
        shutil.copy2(src, dest / src.name)
        logger.info("chosen figure %s", dest / src.name)


def read_table(path: Path) -> list[dict]:
    """Read a CSV table written earlier.

    Args:
        path (Path): Existing table.

    Returns:
        list[dict]: Data rows.
    """
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def main() -> None:
    """Fit the Tmag 11 first batch and save each chosen solution."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--physics",
        help="Refit this physics id only and merge it into the existing tables.",
    )
    args = parser.parse_args()
    if not MANIFEST.is_file():
        raise FileNotFoundError(MANIFEST)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "command.txt").write_text(" ".join(sys.argv) + "\n")
    kept = manifest_rows(args.physics)
    candidates: list[dict] = []
    chosen_rows: list[dict] = []
    for manifest_row in kept:
        cand_rows, chosen = one_curve(manifest_row)
        candidates.extend(cand_rows)
        chosen_rows.append(chosen)
    if args.physics:
        old_candidates = [
            r
            for r in read_table(OUT / "candidates.csv")
            if r["physics"] != args.physics
        ]
        old_chosen = [
            r for r in read_table(OUT / "chosen.csv") if r["physics"] != args.physics
        ]
        candidates = old_candidates + candidates
        chosen_rows = old_chosen + chosen_rows
        candidates.sort(
            key=lambda r: (r["physics"], int(r["cadence_s"]), r["candidate"])
        )
        chosen_rows.sort(key=lambda r: (r["physics"], int(r["cadence_s"])))
    fields_shared = [
        "file",
        "cadence_s",
        "physics",
        "eb_class",
        "recipe",
        "candidate",
        "chosen",
        "status",
        "harmonics",
        "ocon",
        "baseline_order",
        "rms_mmag",
        "two_component",
        "solver",
        "reasons",
        "warnings",
    ]
    write_table(OUT / "candidates.csv", candidates, fields_shared)
    write_table(
        OUT / "chosen.csv",
        chosen_rows,
        [
            "file",
            "cadence_s",
            "physics",
            "eb_class",
            "recipe",
            "n_candidates",
            "chosen",
            "status",
            "harmonics",
            "ocon",
            "baseline_order",
            "rms_mmag",
            "two_component",
            "solver",
            "solution",
            "reasons",
            "warnings",
            "selection_reasons",
        ],
    )
    stale = OUT / "excluded_spotted.csv"
    if stale.is_file():
        stale.unlink()
    copy_chosen_figures(chosen_rows)


if __name__ == "__main__":
    main()
