"""Scoring step 3: times of minima against the generator truth.

Per file: select a model from a small library (baseline order 1, harmonics sets,
O'Connell on when the rules switch it on), time every eclipse with the stored recipe
and compare with ``eclipse_T_min_pure_days`` (the independent bisector of the pure
curve). The file has O-C injected as zero, so a constant offset per eclipse kind is
reported separately from the scatter.

Per eclipse kind:
``bias_s``  mean(T_fit - T_min_pure) in seconds;
``scatter_s``  rms about that mean (the method error);
``honesty``  rms of (T_fit - T_min_pure - bias) / sigma (about 1 if errors are honest);
``conj_offset_s``  T_conj - T_min_pure from the file (astrophysical offset).

Run::

    python examples/score_synthetic_timing.py [file_substring ...]
"""

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path

import numpy as np
from features_table import ROOT  # noqa: F401
from score_synthetic_fits import DEFAULT_RESULTS, OUT, build_library

from mikulasek_eb import select_model, time as time_eclipses
from mikulasek_eb.errors import InvalidLightCurve
from mikulasek_eb.features import extract_features
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.synthetic_check import _array, _num, load_synthetic

logger = logging.getLogger(__name__)
DAY = 86400.0


def score_file(path: Path) -> list[dict]:
    """Time one synthetic file and compare with the pure-curve minima.

    Args:
        path (Path): Synthetic file.

    Returns:
        list[dict]: One row per eclipse kind.

    Raises:
        ValueError: If the truth arrays are missing or inconsistent.
    """
    d = load_synthetic(path)
    meta = d.meta
    kinds = [k.strip() for k in str(meta["eclipse_kind"]).split(",")]
    tmin = _array(meta, "eclipse_T_min_pure_days")
    inj = _array(meta, "eclipse_O_C_injected_days")
    if tmin is None or inj is None or len(kinds) != len(tmin):
        raise ValueError("eclipse truth arrays missing or inconsistent")
    if np.any(inj != 0.0):
        raise ValueError("this scorer assumes O-C injected = 0")
    origin = float(d.timesys["timeorigin"])
    payload = load_tess_votable(path)
    t = absolute_jd_to_mjd(payload["t"])
    feats = extract_features(t, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"])
    lib = [
        c
        for c in build_library(feats)
        if c.config.data.baseline_order == 1 and not c.config.model.ocon
    ]
    res = select_model(
        t, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"],
        features=feats, candidates=lib, data_id=path.stem,
    )
    if res.best is None:
        raise ValueError("no candidate fitted")
    out = res.output
    tim = time_eclipses(out.t_mjd, out.mag, out.sig, out.solution)
    period = _num(meta, "P_true_days") or payload["period"]
    rows = []
    for kind in ("primary", "secondary"):
        T = np.asarray(tim.timings[kind]["T"], float)
        sig = np.asarray(tim.timings[kind]["sigma"], float)
        truth = origin + tmin[[k == kind for k in kinds]]
        if len(T) == 0 or len(truth) == 0:
            rows.append(dict(file=path.name, kind=kind, chosen=res.best.name, n=0))
            continue
        j = np.argmin(np.abs(T[:, None] - truth[None, :]), axis=1)
        dt = T - truth[j]
        keep = np.abs(dt) < period / 4
        dt, s = dt[keep], sig[keep]
        bias = float(np.mean(dt))
        conj = _num(meta, f"conj_minus_min_pure_{kind}_s")
        rows.append(
            dict(
                file=path.name,
                kind=kind,
                chosen=res.best.name,
                n=int(len(dt)),
                n_rejected=int(len(tim.rejected[kind])),
                bias_s=bias * DAY,
                scatter_s=float(np.std(dt)) * DAY,
                median_sigma_s=float(np.median(s)) * DAY,
                honesty=float(np.sqrt(np.mean(((dt - bias) / s) ** 2))),
                conj_offset_s=conj,
            )
        )
    return rows


def _fmt(v):
    """Format a table cell.

    Args:
        v: Value.

    Returns:
        str: Text.
    """
    if v is None:
        return "-"
    return f"{v:.3g}" if isinstance(v, float) else str(v)


def main() -> None:
    """Score the requested files and write ``timing_vs_truth.md/csv``."""
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")
    filters = sys.argv[1:]
    with open(DEFAULT_RESULTS / "manifest.csv", newline="") as fh:
        manifest = [r for r in csv.DictReader(fh) if int(r["fail"]) == 0]
    rows, failures = [], []
    for m in manifest:
        if filters and not any(f in m["path"] for f in filters):
            continue
        path = DEFAULT_RESULTS.parent / m["path"]
        logger.warning("Timing %s", path.name)
        try:
            rows += score_file(path)
        except (InvalidLightCurve, ValueError) as exc:
            failures.append((path.name, f"{type(exc).__name__}: {exc}"))
    cols = ["file", "kind", "chosen", "n", "n_rejected", "bias_s", "scatter_s",
            "median_sigma_s", "honesty", "conj_offset_s"]
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "timing_vs_truth.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(_fmt(r.get(c)) for c in cols) + " |" for r in rows]
    if failures:
        lines += ["", "Not scored:", ""] + [f"- {n}: {m}" for n, m in failures]
    (OUT / "timing_vs_truth.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
