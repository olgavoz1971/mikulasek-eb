"""Scoring step 2b: shape error alone, with the true trend removed from the data.

The generator observation is ``smeared pure + additive trend + noise`` (flux). Here the
true trend is subtracted from the flux before the magnitude is formed, so the fit sees
only the pure curve plus noise. The library is fitted with baseline order 0 (constant
per segment) and compared with the pure smeared curve in mag. Whatever error remains is
the shape error of the profile and harmonics, not the trend.

Run::

    python examples/score_synthetic_shape.py [file_substring ...]
"""

from __future__ import annotations

import csv
import dataclasses
import logging
import sys
from pathlib import Path

import numpy as np
from score_synthetic_fits import GENERATOR, DEFAULT_RESULTS, OUT, build_library
from features_table import ROOT  # noqa: F401

from mikulasek_eb import select_model
from mikulasek_eb.errors import InvalidLightCurve
from mikulasek_eb.features import extract_features
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.synthetic_check import load_synthetic

logger = logging.getLogger(__name__)


def score_file(path: Path) -> list[dict]:
    """Fit order-0 candidates to a trend-free file and score against the pure curve.

    Args:
        path (Path): Synthetic file.

    Returns:
        list[dict]: One row per candidate.

    Raises:
        ValueError: If truth columns are missing or the magnitude zero point is not constant.
    """
    d = load_synthetic(path)
    payload = load_tess_votable(path)
    if d.truth_flux is None or d.truth_trend is None:
        raise ValueError("truth columns missing")
    c = payload["mag"] + 2.5 * np.log10(d.flux)
    if np.ptp(c) > 1e-6:
        raise ValueError("magnitude zero point is not constant")
    c = float(np.median(c))
    flux_clean = d.flux - d.truth_trend
    mag = c - 2.5 * np.log10(flux_clean)
    truth = c - 2.5 * np.log10(d.truth_flux)
    t = absolute_jd_to_mjd(payload["t"])
    feats = extract_features(
        t, mag, payload["sig"], payload["period"], payload["epoch_jd"]
    )
    lib = [
        dataclasses.replace(cand)
        for cand in build_library(feats)
        if cand.config.data.baseline_order == 0
    ]
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
    noise = 1e3 * float(np.std(mag - truth))
    rows = []
    for row in res.table:
        out = res.fits.get(row.name)
        rms = None
        if out is not None:
            idx = np.searchsorted(t, out.t_mjd)
            fit_mag = out.mag - out.pipeline.adopted.resid
            rms = 1e3 * float(np.sqrt(np.mean((fit_mag - truth[idx]) ** 2)))
        rows.append(
            dict(
                file=path.name,
                candidate=row.name,
                status=row.status,
                k=row.k,
                rms_shape=rms,
                noise=noise,
                bic=row.bic,
                chosen=bool(res.best is not None and row.name == res.best.name),
            )
        )
    return rows


def main() -> None:
    """Score the requested files and write ``shape_vs_truth.md/csv``."""
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")
    filters = sys.argv[1:]
    with open(DEFAULT_RESULTS / "manifest.csv", newline="") as fh:
        manifest = [r for r in csv.DictReader(fh) if int(r["fail"]) == 0]
    rows, failures = [], []
    for m in manifest:
        if filters and not any(f in m["path"] for f in filters):
            continue
        path = GENERATOR / m["path"]
        print_name = path.name
        logger.warning("Scoring %s", print_name)
        try:
            rows += score_file(path)
        except (InvalidLightCurve, ValueError) as exc:
            failures.append((path.name, f"{type(exc).__name__}: {exc}"))
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "shape_vs_truth.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    lines = [
        "| File | Noise mmag | Chosen | rms shape chosen | Best | rms shape best |",
        "|---|---|---|---|---|---|",
    ]
    for name in dict.fromkeys(r["file"] for r in rows):
        fr = [r for r in rows if r["file"] == name and r["rms_shape"] is not None]
        best = min(fr, key=lambda r: r["rms_shape"])
        ch = next((r for r in fr if r["chosen"]), best)
        lines.append(
            f"| {name} | {fr[0]['noise']:.2f} | {ch['candidate']} | {ch['rms_shape']:.3f} "
            f"| {best['candidate']} | {best['rms_shape']:.3f} |"
        )
    if failures:
        lines += ["", "Not scored:", ""] + [f"- {n}: {m}" for n, m in failures]
    (OUT / "shape_vs_truth.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
