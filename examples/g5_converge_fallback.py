"""G5 only: trust-region retry of every candidate that did not converge.

The batch selection drops a candidate that raises ``FitDidNotConverge``.
This run asks ``select_model`` to fit that candidate again with the
trust-region solver, then chooses by the same evidence. The three Tmag 11
G5 files are the only input. No truth column is read, and nothing is timed.
Other profiles are not tried.

Run from the repository root::

    MPLBACKEND=Agg python examples/g5_converge_fallback.py

Writes ``results/synthetic_scoring/g5_converge_fallback/``.
"""

from __future__ import annotations

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
    / "g5_converge_fallback"
)
N_EXPECTED = 3


def g5_rows() -> list[dict]:
    """Accepted G5 rows of the Tmag 11 first-batch manifest.

    Returns:
        list[dict]: One row per cadence, ordered by cadence.

    Raises:
        RuntimeError: If the number of accepted G5 files is not 3.
    """
    with open(MANIFEST, newline="") as fh:
        rows = [
            r
            for r in csv.DictReader(fh)
            if r["physics_id"] == "G5" and int(r["fail"]) == 0
        ]
    rows.sort(key=lambda r: int(r["cadence_s"]))
    if len(rows) != N_EXPECTED:
        raise RuntimeError(f"expected {N_EXPECTED} G5 files, found {len(rows)}")
    return rows


def save_figure(stem: str, payload: dict, out: eb.FitOutput, title: str) -> None:
    """Save the profile figure of one candidate.

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


def one_curve(manifest_row: dict) -> tuple[list[dict], dict]:
    """Select a model, retrying any candidate that did not converge.

    Args:
        manifest_row (dict): One accepted G5 row.

    Returns:
        tuple[list[dict], dict]: Candidate rows, then the chosen-fit row.
    """
    curve = GENERATED / manifest_row["path"]
    payload = eb.load_tess_votable(curve)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    result = eb.select_model(
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        data_id=curve.stem,
        retry_not_converged=True,
    )
    chosen_name = None if result.best is None else result.best.name
    logger.info(
        "G5 %ss chosen %s",
        manifest_row["cadence_s"],
        chosen_name if chosen_name else "none",
    )
    cand_rows = []
    for row in result.table:
        out = result.fits.get(row.name)
        record = dict(
            file=curve.name,
            cadence_s=manifest_row["cadence_s"],
            candidate=row.name,
            chosen=int(row.name == chosen_name),
            status=row.status,
            harmonics="",
            baseline_order="",
            rms_mmag="",
            solver="",
        )
        if out is not None:
            model = out.solution.config.model
            stats = out.solution.stats
            record.update(
                harmonics=" ".join(str(h) for h in model.harmonics),
                baseline_order=out.solution.config.data.baseline_order,
                rms_mmag=f"{1e3 * float(stats['rms_mag']):.2f}",
                solver=solver_name(out),
            )
            mark = "chosen" if row.name == chosen_name else "not chosen"
            title = (
                f"G5  {manifest_row['cadence_s']} s  {row.name}  {mark}  "
                f"{stats['status']}  rms {1e3 * float(stats['rms_mag']):.2f} mmag"
            )
            save_figure(f"{curve.stem}__{row.name}", payload, out, title)
        cand_rows.append(record)
    solution_path = ""
    chosen = dict(
        file=curve.name,
        cadence_s=manifest_row["cadence_s"],
        chosen=chosen_name or "",
        status=(
            "none" if result.output is None else result.output.solution.stats["status"]
        ),
        harmonics="",
        baseline_order="",
        rms_mmag="",
        solver="",
        solution="",
        selection_reasons="; ".join(result.reasons),
    )
    if result.output is not None:
        solutions = OUT / "solutions"
        solutions.mkdir(parents=True, exist_ok=True)
        destination = solutions / f"{curve.stem}.json"
        result.output.solution.save(destination)
        solution_path = str(destination.relative_to(OUT))
        model = result.output.solution.config.model
        stats = result.output.solution.stats
        chosen.update(
            harmonics=" ".join(str(h) for h in model.harmonics),
            baseline_order=result.output.solution.config.data.baseline_order,
            rms_mmag=f"{1e3 * float(stats['rms_mag']):.2f}",
            solver=solver_name(result.output),
            solution=solution_path,
        )
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
    logger.info("wrote %s", path)


def main() -> None:
    """Refit G5 with a trust-region retry after every convergence failure."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if not MANIFEST.is_file():
        raise FileNotFoundError(MANIFEST)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "command.txt").write_text(" ".join(sys.argv) + "\n")
    candidates: list[dict] = []
    chosen_rows: list[dict] = []
    for manifest_row in g5_rows():
        cand_rows, chosen = one_curve(manifest_row)
        candidates.extend(cand_rows)
        chosen_rows.append(chosen)
    write_table(
        OUT / "candidates.csv",
        candidates,
        [
            "file",
            "cadence_s",
            "candidate",
            "chosen",
            "status",
            "harmonics",
            "baseline_order",
            "rms_mmag",
            "solver",
        ],
    )
    write_table(
        OUT / "chosen.csv",
        chosen_rows,
        [
            "file",
            "cadence_s",
            "chosen",
            "status",
            "harmonics",
            "baseline_order",
            "rms_mmag",
            "solver",
            "solution",
            "selection_reasons",
        ],
    )
    dest = OUT / "fits" / "chosen"
    dest.mkdir(parents=True, exist_ok=True)
    for row in chosen_rows:
        if not row["chosen"]:
            raise RuntimeError(f"{row['file']} has no chosen candidate")
        stem = Path(row["file"]).stem
        src = OUT / "fits" / f"{stem}__{row['chosen']}.png"
        if not src.is_file():
            raise FileNotFoundError(src)
        shutil.copy2(src, dest / src.name)
        logger.info("chosen figure %s", dest / src.name)


if __name__ == "__main__":
    main()
