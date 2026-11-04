"""From-scratch Mikulasek fits of the Tmag 11 spot curves, at every cadence.

Each file is fitted from its own candidates. The harmonic orders, the
O'Connell term and the baseline are chosen on that file. The only inputs
are the degraded light curve: time, magnitude, magnitude error, and the
period and epoch stored in that file. No truth flux, no conjunction time,
no bisector, and no timing.

The manifest is used only to list the files. Spot longitude is a label of
the file, not an input to the fit.

Run from the repository root::

    MPLBACKEND=Agg python examples/spot_m11_fit.py
    MPLBACKEND=Agg python examples/spot_m11_fit.py --cadence 1800

The full run writes ``results/synthetic_scoring/spot_m11_fit/``.
``--cadence`` writes ``results/synthetic_scoring/spot_m11_fit_c<cadence>/``
and does not overwrite the full run. Selection uses the package default:
Levenberg-Marquardt, then the trust-region solver for a candidate marked
``bound-active``.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

import mikulasek_eb as eb
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.plotting import plot_global_fit

logger = logging.getLogger(__name__)
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFEST = GENERATED / "results" / "m11" / "manifest_spotsweep.csv"
SCORING = Path(__file__).resolve().parents[1] / "results" / "synthetic_scoring"
OUT = SCORING / "spot_m11_fit"
N_EXPECTED = 27


def manifest_rows(cadence_s: int | None) -> list[dict]:
    """Accepted rows of the Tmag 11 spot-sweep manifest.

    Args:
        cadence_s (int | None): Keep this cadence only. None keeps every cadence.

    Returns:
        list[dict]: One row per light curve, ordered by cadence then physics id.

    Raises:
        RuntimeError: If the number of accepted files is not the expected count.
    """
    with open(MANIFEST, newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if int(r["fail"]) == 0]
    if cadence_s is not None:
        rows = [r for r in rows if int(r["cadence_s"]) == cadence_s]
    rows.sort(key=lambda r: (int(r["cadence_s"]), r["physics_id"]))
    expected = 9 if cadence_s is not None else N_EXPECTED
    if len(rows) != expected:
        raise RuntimeError(f"expected {expected} accepted curves, found {len(rows)}")
    return rows


def save_figure(stem: str, payload: dict, out: eb.FitOutput, title: str) -> None:
    """Save the profile figure of one candidate.

    The figure is drawn from the fit and from the magnitudes the fit used.
    Nothing is read from a truth column.

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


def candidate_row(
    curve: Path,
    cadence_s: str,
    physics: str,
    name: str,
    chosen: bool,
    out: eb.FitOutput | None,
    status: str,
) -> dict:
    """One line describing a candidate, whether or not it was chosen.

    Args:
        curve (Path): Light curve.
        cadence_s (str): Cadence in seconds, from the manifest.
        physics (str): ``physics_id`` from the manifest, a file label.
        name (str): Candidate name.
        chosen (bool): True when ``select_model`` chose this candidate.
        out (FitOutput | None): Its fit, or None when the fit raised.
        status (str): Verdict, or the exception text.

    Returns:
        dict: A row of the candidate table.
    """
    row = dict(
        file=curve.name,
        cadence_s=cadence_s,
        physics=physics,
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
        solver=(
            "trust-region"
            if out.solution.config.policy.bounded_solver
            else "levenberg-marquardt"
        ),
        reasons="; ".join(stats["status_reasons"]),
        warnings="; ".join(stats["status_warnings"]),
    )
    return row


def one_curve(manifest_row: dict) -> tuple[list[dict], dict]:
    """Select a model for one degraded light curve and draw every candidate.

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
                curve,
                cadence_s,
                physics,
                row.name,
                row.name == chosen_name,
                out,
                row.status,
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
    chosen = dict(
        file=curve.name,
        cadence_s=cadence_s,
        physics=physics,
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
            solver=(
                "trust-region"
                if result.output.solution.config.policy.bounded_solver
                else "levenberg-marquardt"
            ),
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


def main() -> None:
    """Fit the Tmag 11 spot curves from the degraded light curve alone."""
    global OUT
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cadence",
        type=int,
        choices=(120, 600, 1800),
        help="Fit this cadence only. Writes a separate directory.",
    )
    args = parser.parse_args()
    if args.cadence is not None:
        OUT = SCORING / f"spot_m11_fit_c{args.cadence}"
    if not MANIFEST.is_file():
        raise FileNotFoundError(MANIFEST)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "command.txt").write_text(" ".join(sys.argv) + "\n")
    logger.info("output %s", OUT)
    candidates: list[dict] = []
    chosen_rows: list[dict] = []
    for manifest_row in manifest_rows(args.cadence):
        cand_rows, chosen = one_curve(manifest_row)
        candidates.extend(cand_rows)
        chosen_rows.append(chosen)
    write_table(
        OUT / "candidates.csv",
        candidates,
        [
            "file",
            "cadence_s",
            "physics",
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
        ],
    )
    write_table(
        OUT / "chosen.csv",
        chosen_rows,
        [
            "file",
            "cadence_s",
            "physics",
            "n_candidates",
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
            "selection_reasons",
        ],
    )


if __name__ == "__main__":
    main()
