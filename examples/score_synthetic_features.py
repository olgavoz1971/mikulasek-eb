"""Scoring step 1: features measured by our code against the generator labels.

Reads every accepted row of the generator ``manifest.csv``, measures the features
of the synthetic light curve (``extract_features``) and compares them with the
labels written in the file: eclipse depth, FWHM of each eclipse, depth ratio.
Writes ``results/synthetic_scoring/features_vs_labels.md`` and ``.csv``.

Run::

    python examples/score_synthetic_features.py
    python examples/score_synthetic_features.py /path/to/EB_synthetic_lc/results
"""

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path

from features_table import ROOT

from mikulasek_eb import extract_features
from mikulasek_eb.errors import InvalidLightCurve
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.synthetic_check import _flag, _num, load_synthetic

logger = logging.getLogger(__name__)
DEFAULT_RESULTS = Path("/home/voz/projects/UPJS/EB_synthetic_lc/results")
OUT = ROOT / "results" / "synthetic_scoring"


def _pct(measured, truth):
    """Relative difference in per cent.

    Args:
        measured (float | None): Measured value.
        truth (float | None): Label value.

    Returns:
        float | None: 100 * (measured - truth) / truth, None if either is missing.
    """
    if measured is None or truth is None or truth == 0:
        return None
    return 100.0 * (measured - truth) / truth


def score_file(path: Path) -> dict:
    """Measure features of one synthetic file and compare with its labels.

    Args:
        path (Path): Synthetic ``*_c<cadence>.vot`` file.

    Returns:
        dict: One row of measured values, labels and differences.
    """
    meta = load_synthetic(path, with_errors=False).meta
    payload = load_tess_votable(path)
    f = extract_features(
        absolute_jd_to_mjd(payload["t"]),
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
    )
    g, s, a = f.geometry, f.shape, f.sampling
    p_true = _num(meta, "P_true_days") or payload["period"]
    fwhm1 = None if g["fwhm1"] is None else g["fwhm1"] * p_true
    fwhm2 = None if g["fwhm2"] is None else g["fwhm2"] * p_true
    lab_f1 = _num(meta, "eclipse_fwhm_primary_days")
    lab_f2 = _num(meta, "eclipse_fwhm_secondary_days")
    lab_depth = _num(meta, "eclipse_depth_mag")
    lab_ratio = (
        None if lab_f1 is None or lab_f2 is None or lab_f1 == 0 else lab_f2 / lab_f1
    )
    return {
        "file": path.name,
        "physics": meta.get("physics_id"),
        "cadence_s": _num(meta, "cadence_s"),
        "A1_meas": g["A1"],
        "A1_label": lab_depth,
        "A1_diff_pct": _pct(g["A1"], lab_depth),
        "fwhm1_meas_d": fwhm1,
        "fwhm1_label_d": lab_f1,
        "fwhm1_diff_pct": _pct(fwhm1, lab_f1),
        "fwhm2_meas_d": fwhm2,
        "fwhm2_label_d": lab_f2,
        "fwhm2_diff_pct": _pct(fwhm2, lab_f2),
        "fwhm_ratio_meas": g["fwhm_ratio"],
        "fwhm_ratio_label": lab_ratio,
        "secondary_offset_phase": (g["secondary_offset"] or {}).get("value"),
        "eccentric_label": _flag(meta, "eccentric"),
        "total_label": _flag(meta, "total_eclipse"),
        "flat_bottom_meas": s["flat_bottom"],
        "gamma_estimate": s["gamma_estimate"],
        "noise_over_catalogue": a["noise_over_catalogue"],
        "points_per_fwhm": g["points_per_fwhm"],
        "n_notes": len(f.notes),
    }


def _fmt(v, digits=3):
    """Format a value for the table.

    Args:
        v: Number, bool, text or None.
        digits (int): Significant digits.

    Returns:
        str: Text.
    """
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.{digits}g}"
    return str(v)


def main() -> None:
    """Score all accepted manifest rows and write the table."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    results = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_RESULTS
    with open(results / "manifest.csv", newline="") as fh:
        manifest = [r for r in csv.DictReader(fh) if int(r["fail"]) == 0]
    rows, failures = [], []
    for m in manifest:
        path = results.parent / m["path"]
        logger.info("Scoring %s", path.name)
        try:
            rows.append(score_file(path))
        except InvalidLightCurve as exc:
            logger.error("%s: features failed: %s", path.name, exc)
            failures.append((path.name, str(exc)))
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "features_vs_labels.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    cols = [
        ("file", "File"),
        ("A1_meas", "A1 meas"),
        ("A1_label", "A1 label"),
        ("A1_diff_pct", "A1 diff %"),
        ("fwhm1_diff_pct", "FWHM1 diff %"),
        ("fwhm2_diff_pct", "FWHM2 diff %"),
        ("fwhm_ratio_meas", "ratio meas"),
        ("fwhm_ratio_label", "ratio label"),
        ("secondary_offset_phase", "sec offset (phase)"),
        ("eccentric_label", "ecc"),
        ("flat_bottom_meas", "flat meas"),
        ("total_label", "total"),
        ("points_per_fwhm", "pts/FWHM"),
    ]
    lines = [
        "| " + " | ".join(c[1] for c in cols) + " |",
        "|" + "---|" * len(cols),
    ]
    for r in rows:
        lines.append("| " + " | ".join(_fmt(r[c[0]]) for c in cols) + " |")
    if failures:
        lines += ["", "Features not measured (typed error, not skipped silently):", ""]
        lines += [f"- {n}: {msg}" for n, msg in failures]
    (OUT / "features_vs_labels.md").write_text("\n".join(lines) + "\n")
    logger.info("Wrote %s", OUT / "features_vs_labels.md")


if __name__ == "__main__":
    main()
