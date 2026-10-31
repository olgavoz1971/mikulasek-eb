"""Test: offer the second profile component as a candidate and let BIC decide.

Each harmonic set (baseline order 1, O'Connell on and off) is fitted twice: ``single`` (the
second component never tried) and ``two`` (``FitPolicy.force_two_component``: the trial always runs and
the rms gain criterion is dropped; the other acceptance checks stay). ``select_model`` compares them by BIC.
Reported per file: the BIC winner, whether a second component was really adopted, and the
error against the noise-free observation, so the BIC choice can be judged by the truth.

Run::

    python examples/score_synthetic_twocomp.py [file_substring ...]
"""

from __future__ import annotations

import csv
import dataclasses
import logging
import sys

import numpy as np
from score_synthetic_fits import DEFAULT_RESULTS, OUT, build_library, noise_free_mag

from mikulasek_eb import select_model
from mikulasek_eb.candidates import Candidate
from mikulasek_eb.errors import InvalidLightCurve
from mikulasek_eb.features import extract_features
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable

logger = logging.getLogger(__name__)
HARMONICS = ((2,), (2, 3), (1, 2, 3))


def score_file(path) -> list[dict]:
    """Fit single and forced two-component variants and score them.

    Args:
        path (Path): Synthetic file.

    Returns:
        list[dict]: One row per candidate.
    """
    payload = load_tess_votable(path)
    t = absolute_jd_to_mjd(payload["t"])
    truth, mag = noise_free_mag(path, payload)
    feats = extract_features(t, mag, payload["sig"], payload["period"], payload["epoch_jd"])
    cands = []
    for c in build_library(feats):
        if c.config.data.baseline_order != 1 or c.config.model.harmonics not in HARMONICS:
            continue
        pol = c.config.policy
        single = dataclasses.replace(pol, allow_two_component=False)
        two = dataclasses.replace(pol, allow_two_component=True, force_two_component=True)
        cands.append(Candidate(c.name + "_single", dataclasses.replace(c.config, policy=single).validate(), ()))
        cands.append(Candidate(c.name + "_two", dataclasses.replace(c.config, policy=two).validate(), ()))
    res = select_model(t, mag, payload["sig"], payload["period"], payload["epoch_jd"], features=feats, candidates=cands, data_id=path.stem)
    rows = []
    for row in res.table:
        out = res.fits.get(row.name)
        rms_truth = adopted = None
        if out is not None:
            idx = np.searchsorted(t, out.t_mjd)
            rms_truth = 1e3 * float(np.sqrt(np.mean((out.mag - out.pipeline.adopted.resid - truth[idx]) ** 2)))
            adopted = bool(out.solution.stats["two_component"])
        rows.append(dict(file=path.name, candidate=row.name, status=row.status, k=row.k, bic=row.bic, delta_bic=row.delta_bic, rms_truth=rms_truth, two_adopted=adopted, chosen=bool(res.best and res.best.name == row.name)))
    return rows


def main() -> None:
    """Run the test and write ``twocomp_vs_truth.md/csv``."""
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    filters = sys.argv[1:] or ["_c"]
    with open(DEFAULT_RESULTS / "manifest.csv", newline="") as fh:
        manifest = [r for r in csv.DictReader(fh) if int(r["fail"]) == 0]
    rows, failures = [], []
    for m in manifest:
        if not any(f in m["path"] for f in filters):
            continue
        path = DEFAULT_RESULTS.parent / m["path"]
        logger.warning("Scoring %s", path.name)
        try:
            rows += score_file(path)
        except (InvalidLightCurve, ValueError) as exc:
            failures.append((path.name, f"{type(exc).__name__}: {exc}"))
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "twocomp_vs_truth.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    lines = ["| File | BIC winner | two adopted | rms truth winner | best single (rms) | best two-adopted (rms) | dBIC best two vs best single |", "|---|---|---|---|---|---|---|"]
    for name in dict.fromkeys(r["file"] for r in rows):
        fr = [r for r in rows if r["file"] == name and r["bic"] is not None]
        win = next(r for r in fr if r["chosen"])
        sg = min((r for r in fr if r["candidate"].endswith("_single")), key=lambda r: r["bic"])
        tw = [r for r in fr if r["candidate"].endswith("_two")]
        twb = min(tw, key=lambda r: r["bic"]) if tw else None
        sg_t = min((r["rms_truth"] for r in fr if r["candidate"].endswith("_single")))
        lines.append(
            f"| {name} | {win['candidate']} | {win['two_adopted']} | {win['rms_truth']:.2f} | {sg['candidate']} ({sg_t:.2f}) "
            f"| {twb['candidate'] + ' (' + format(twb['rms_truth'], '.2f') + ')' if twb else '-'}{'' if twb and twb['two_adopted'] else ' not adopted'} "
            f"| {twb['bic'] - sg['bic']:.1f} |"
        )
    if failures:
        lines += ["", "Not scored:", ""] + [f"- {n}: {m}" for n, m in failures]
    (OUT / "twocomp_vs_truth.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
