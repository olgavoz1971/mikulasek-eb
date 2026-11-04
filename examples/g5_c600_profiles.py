"""G5 at 600 s: profiles 1 to 4, trust-region fallback off, every picture kept.

The file is the Tmag 11 G5 light curve at 600 s. Candidates are the usual
profile-1 set, plus the same harmonics and baselines on Table 2 profiles 2, 3
and 4. The staged fit does not run profiles 5 to 7. ``retry_bound_active`` is
false on every candidate, and a fit that does not converge is not retried.
No truth column is read, and nothing is timed. Every finished candidate keeps
its profile figure and its fitted-trend figure.

Run from the repository root::

    MPLBACKEND=Agg python examples/g5_c600_profiles.py

Writes ``results/synthetic_scoring/g5_c600_profiles/``.
"""

from __future__ import annotations

import csv
import dataclasses
import logging
import shutil
import sys
from pathlib import Path

import mikulasek_eb as eb
from mikulasek_eb.candidates import Candidate
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.plotting import plot_global_fit, plot_instrumental_trend

logger = logging.getLogger(__name__)
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFEST = GENERATED / "results" / "m11" / "manifest.csv"
OUT = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "synthetic_scoring"
    / "g5_c600_profiles"
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


def profile_candidates(features) -> list:
    """Profile-1 candidates, repeated for profiles 1 to 4, fallback off.

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
                policy=dataclasses.replace(
                    cand.config.policy, retry_bound_active=False
                ),
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


def one_curve(manifest_row: dict) -> tuple[list[dict], dict]:
    """Fit profiles 1 to 4 on the 600 s G5 file and keep every picture.

    Args:
        manifest_row (dict): The accepted G5 row at 600 s.

    Returns:
        tuple[list[dict], dict]: Candidate rows, then the chosen-fit row.
    """
    curve = GENERATED / manifest_row["path"]
    payload = eb.load_tess_votable(curve)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    features = eb.extract_features(
        t_mjd, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    result = eb.select_model(
        t_mjd,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        features=features,
        candidates=profile_candidates(features),
        data_id=curve.stem,
        retry_not_converged=False,
    )
    chosen_name = None if result.best is None else result.best.name
    logger.info("G5 600 s chosen %s", chosen_name if chosen_name else "none")
    cand_rows = []
    for row in result.table:
        out = result.fits.get(row.name)
        record = dict(
            file=curve.name,
            candidate=row.name,
            profile=row.name.split("_", 1)[0][1:],
            chosen=int(row.name == chosen_name),
            status=row.status,
            harmonics="",
            baseline_order="",
            rms_mmag="",
        )
        if out is not None:
            model = out.solution.config.model
            stats = out.solution.stats
            record.update(
                harmonics=" ".join(str(h) for h in model.harmonics),
                baseline_order=out.solution.config.data.baseline_order,
                rms_mmag=f"{1e3 * float(stats['rms_mag']):.2f}",
            )
            mark = "chosen" if row.name == chosen_name else "not chosen"
            title = (
                f"G5  600 s  {row.name}  {mark}  {stats['status']}  "
                f"rms {1e3 * float(stats['rms_mag']):.2f} mmag"
            )
            save_figures(f"{curve.stem}__{row.name}", payload, out, title)
        cand_rows.append(record)
    chosen = dict(
        file=curve.name,
        chosen=chosen_name or "",
        status=(
            "none" if result.output is None else result.output.solution.stats["status"]
        ),
        harmonics="",
        baseline_order="",
        rms_mmag="",
        solution="",
        selection_reasons="; ".join(result.reasons),
    )
    if result.output is not None:
        solutions = OUT / "solutions"
        solutions.mkdir(parents=True, exist_ok=True)
        destination = solutions / f"{curve.stem}.json"
        result.output.solution.save(destination)
        model = result.output.solution.config.model
        stats = result.output.solution.stats
        chosen.update(
            harmonics=" ".join(str(h) for h in model.harmonics),
            baseline_order=result.output.solution.config.data.baseline_order,
            rms_mmag=f"{1e3 * float(stats['rms_mag']):.2f}",
            solution=str(destination.relative_to(OUT)),
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


def copy_chosen(name: str, folder: str) -> None:
    """Copy one chosen picture into the chosen subfolder.

    Args:
        name (str): Picture file name.
        folder (str): ``fits`` or ``trends``.

    Raises:
        FileNotFoundError: If the source picture is missing.
    """
    src = OUT / folder / name
    if not src.is_file():
        raise FileNotFoundError(src)
    dest = OUT / folder / "chosen"
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest / name)
    logger.info("chosen %s", dest / name)


def main() -> None:
    """Fit G5 at 600 s on profiles 1 to 4 and save every picture."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if not MANIFEST.is_file():
        raise FileNotFoundError(MANIFEST)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "command.txt").write_text(" ".join(sys.argv) + "\n")
    candidates, chosen = one_curve(g5_600_row())
    write_table(
        OUT / "candidates.csv",
        candidates,
        [
            "file",
            "candidate",
            "profile",
            "chosen",
            "status",
            "harmonics",
            "baseline_order",
            "rms_mmag",
        ],
    )
    write_table(
        OUT / "chosen.csv",
        [chosen],
        [
            "file",
            "chosen",
            "status",
            "harmonics",
            "baseline_order",
            "rms_mmag",
            "solution",
            "selection_reasons",
        ],
    )
    if not chosen["chosen"]:
        raise RuntimeError(f"{chosen['file']} has no chosen candidate")
    stem = Path(chosen["file"]).stem
    picture = f"{stem}__{chosen['chosen']}.png"
    copy_chosen(picture, "fits")
    copy_chosen(picture, "trends")


if __name__ == "__main__":
    main()
