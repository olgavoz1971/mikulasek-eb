"""Scoring step 2: fit a library of candidates to synthetic files, score against truth.

For every accepted file of the generator manifest it fits a fixed library of
candidates (harmonics x baseline order x O'Connell on/off; exposure and gap taken
from the rules), lets ``select_model`` choose by evidence, and compares each fit
with the noise-free observation (smeared pure curve plus additive trend, in mag).

Metrics (mmag): ``noise`` = data minus noise-free (the floor), ``rms_truth`` = fit
minus noise-free on the points the fit used, ``regret`` = rms_truth of the chosen
candidate minus the smallest rms_truth in the library.

Run::

    python examples/score_synthetic_fits.py [results_dir] [file_substring ...]
"""

from __future__ import annotations

import csv
import dataclasses
import logging
import sys
import time
from pathlib import Path

import numpy as np
from features_table import ROOT

from mikulasek_eb import select_model
from mikulasek_eb.candidates import Candidate, propose_candidates
from mikulasek_eb.config import Config, ModelSpec
from mikulasek_eb.errors import InvalidLightCurve
from mikulasek_eb.features import extract_features
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.synthetic_check import load_synthetic

logger = logging.getLogger(__name__)
GENERATOR = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
DEFAULT_RESULTS = GENERATOR / "results" / "n01"
OUT = ROOT / "results" / "synthetic_scoring"
HARMONIC_SETS = ((), (2,), (2, 3), (1, 2, 3))
BASELINE_ORDERS = (0, 1)


def build_library(features, rules=None) -> list:
    """Fixed candidate library sharing exposure and gap handling with the rules.

    Args:
        features (Features): Measured features of the light curve.
        rules (RuleSettings | None): Rule thresholds.

    Returns:
        list[Candidate]: Harmonics x baseline order x O'Connell on/off.
    """
    base = propose_candidates(features, rules)[0].config
    out = []
    for harm in HARMONIC_SETS:
        for order in BASELINE_ORDERS:
            for ocon in (False, True):
                name = (
                    "h"
                    + "".join(map(str, harm))
                    + f"_b{order}"
                    + ("_oc" if ocon else "")
                )
                cfg = dataclasses.replace(
                    base,
                    model=ModelSpec(profile=1, harmonics=harm, ocon=ocon),
                    data=dataclasses.replace(base.data, baseline_order=order),
                ).validate()
                out.append(Candidate(name, cfg, ("library",)))
    return out


def noise_free_mag(path: Path, payload: dict):
    """Noise-free observed magnitudes (smeared pure curve plus additive trend).

    The magnitude zero point is taken from the data themselves through
    ``mag = c - 2.5 log10(flux)``; the scatter of ``c`` over the file is checked.

    Args:
        path (Path): Synthetic file.
        payload (dict): ``load_tess_votable`` output of the same file.

    Returns:
        tuple[ndarray, ndarray]: Noise-free magnitudes and the data magnitudes.

    Raises:
        ValueError: If the truth columns are missing, sizes differ or the zero
            point is not constant.
    """
    d = load_synthetic(path)
    if d.truth_flux is None or d.truth_trend is None:
        raise ValueError("truth_flux_pure_smeared or truth_trend column is missing")
    if len(d.flux) != len(payload["mag"]):
        raise ValueError("point counts of the flux and magnitude readers differ")
    c = payload["mag"] + 2.5 * np.log10(d.flux)
    if np.ptp(c) > 1e-6:
        raise ValueError(
            f"mag = c - 2.5 log10(flux) does not hold (ptp c = {np.ptp(c):.2e})"
        )
    return (
        float(np.median(c)) - 2.5 * np.log10(d.truth_flux + d.truth_trend),
        payload["mag"],
    )


def score_file(path: Path) -> list[dict]:
    """Fit the library to one file and score each candidate and the selector.

    Args:
        path (Path): Synthetic file.

    Returns:
        list[dict]: One row per candidate with metrics and selector flags.
    """
    payload = load_tess_votable(path)
    t = absolute_jd_to_mjd(payload["t"])
    truth, mag = noise_free_mag(path, payload)
    feats = extract_features(
        t, mag, payload["sig"], payload["period"], payload["epoch_jd"]
    )
    lib = build_library(feats)
    t0 = time.time()
    res = select_model(
        t,
        mag,
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        features=feats,
        candidates=lib,
        data_id=path.stem,
    )
    seconds = time.time() - t0
    key = lambda c: (
        c.config.model.harmonics,
        c.config.data.baseline_order,
        c.config.model.ocon,
    )  # noqa: E731
    proposed = {key(c) for c in propose_candidates(feats)}
    library_key = {c.name: key(c) for c in lib}
    noise = 1e3 * float(np.std(mag - truth))
    rows = []
    for row in res.table:
        out = res.fits.get(row.name)
        rms_truth = None
        if out is not None:
            idx = np.searchsorted(t, out.t_mjd)
            fit_mag = out.mag - out.pipeline.adopted.resid
            rms_truth = 1e3 * float(np.sqrt(np.mean((fit_mag - truth[idx]) ** 2)))
        rows.append(
            dict(
                file=path.name,
                candidate=row.name,
                status=row.status,
                k=row.k,
                rms_data=row.rms_mmag,
                rms_truth=rms_truth,
                noise=noise,
                bic=row.bic,
                delta_bic=row.delta_bic,
                chosen=bool(res.best is not None and row.name == res.best.name),
                proposed=library_key[row.name] in proposed,
                ask_human=res.ask_human,
                seconds=seconds,
            )
        )
    return rows


def summarise(rows: list[dict]) -> list[str]:
    """Markdown summary: one line per file with chosen, best and regret.

    Args:
        rows (list[dict]): Rows from ``score_file``.

    Returns:
        list[str]: Markdown lines.
    """
    lines = [
        "| File | Noise mmag | Chosen | rms truth chosen | Best in library | rms truth best "
        "| Regret mmag | Best proposed | Proposed regret | ask_human |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for name in dict.fromkeys(r["file"] for r in rows):
        fr = [r for r in rows if r["file"] == name and r["rms_truth"] is not None]
        if not fr:
            lines.append(f"| {name} | - | all failed | | | | | | | |")
            continue
        best = min(fr, key=lambda r: r["rms_truth"])
        chosen = next((r for r in fr if r["chosen"]), None)
        prop = [r for r in fr if r["proposed"]]
        bp = min(prop, key=lambda r: r["rms_truth"]) if prop else None
        if chosen is None:
            lines.append(
                f"| {name} | {fr[0]['noise']:.2f} | none chosen | | {best['candidate']} "
                f"| {best['rms_truth']:.3f} | | | | yes |"
            )
            continue
        proposed_regret = (
            "-" if bp is None else f"{bp['rms_truth'] - best['rms_truth']:.3f}"
        )
        lines.append(
            f"| {name} | {fr[0]['noise']:.2f} | {chosen['candidate']} "
            f"| {chosen['rms_truth']:.3f} | {best['candidate']} | {best['rms_truth']:.3f} "
            f"| {chosen['rms_truth'] - best['rms_truth']:.3f} "
            f"| {bp['candidate'] if bp else 'not in library'} "
            f"| {proposed_regret} "
            f"| {'yes' if fr[0]['ask_human'] else 'no'} |"
        )
    return lines


def resummarise() -> None:
    """Rebuild the markdown summary from the saved CSV without refitting."""
    num = ("k", "rms_data", "rms_truth", "noise", "bic", "delta_bic", "seconds")
    flags = ("chosen", "proposed", "ask_human")
    with open(OUT / "fits_vs_truth.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in num:
            r[k] = float(r[k]) if r[k] != "" else None
        for k in flags:
            r[k] = r[k] == "True"
    (OUT / "fits_vs_truth.md").write_text("\n".join(summarise(rows)) + "\n")


def main() -> None:
    """Score the requested files and write the tables."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = sys.argv[1:]
    if args == ["--resummarise"]:
        resummarise()
        return
    results = Path(args[0]) if args and Path(args[0]).is_dir() else DEFAULT_RESULTS
    filters = [a for a in args if not Path(a).is_dir()]
    with open(results / "manifest.csv", newline="") as fh:
        manifest = [r for r in csv.DictReader(fh) if int(r["fail"]) == 0]
    rows, failures = [], []
    for m in manifest:
        if filters and not any(f in m["path"] for f in filters):
            continue
        path = GENERATOR / m["path"]
        logger.info("Scoring %s", path.name)
        try:
            rows += score_file(path)
        except (InvalidLightCurve, ValueError) as exc:
            logger.error("%s: %s", path.name, exc)
            failures.append((path.name, f"{type(exc).__name__}: {exc}"))
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "fits_vs_truth.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    lines = summarise(rows)
    if failures:
        lines += ["", "Not scored (typed error, not skipped silently):", ""]
        lines += [f"- {n}: {msg}" for n, msg in failures]
    (OUT / "fits_vs_truth.md").write_text("\n".join(lines) + "\n")
    logger.info("Wrote %s", OUT / "fits_vs_truth.md")


if __name__ == "__main__":
    main()
