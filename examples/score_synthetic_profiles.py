"""Offline experiment: Table 2 profiles 1 to 4 on synthetic files, shape error alone.

The true trend is removed (as in ``score_synthetic_shape.py``), baseline order 0, and
each profile is fitted with several harmonic sets (O'Connell on and off). The error is
the rms against the pure smeared curve in mmag. Profiles 5 to 7 are not in the staged
pipeline (they have no Gamma) and are not included here.

Run::

    python examples/score_synthetic_profiles.py [file_substring ...]
"""

from __future__ import annotations

import csv
import dataclasses
import logging
import sys

import numpy as np
from score_synthetic_fits import GENERATOR, DEFAULT_RESULTS, OUT, build_library

from mikulasek_eb import fit
from mikulasek_eb.config import ModelSpec
from mikulasek_eb.errors import MikulasekError
from mikulasek_eb.features import extract_features
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.synthetic_check import load_synthetic

logger = logging.getLogger(__name__)
PROFILES = (1, 2, 3, 4)


def run_file(path) -> list[dict]:
    """Fit profiles 1 to 4 with a few harmonic sets to one trend-free file.

    Args:
        path (Path): Synthetic file.

    Returns:
        list[dict]: One row per fit (profile, harmonics, O'Connell, errors).
    """
    d = load_synthetic(path)
    payload = load_tess_votable(path)
    c = float(np.median(payload["mag"] + 2.5 * np.log10(d.flux)))
    mag = c - 2.5 * np.log10(d.flux - d.truth_trend)
    truth = c - 2.5 * np.log10(d.truth_flux)
    t = absolute_jd_to_mjd(payload["t"])
    feats = extract_features(
        t, mag, payload["sig"], payload["period"], payload["epoch_jd"]
    )
    base = [x for x in build_library(feats) if x.config.data.baseline_order == 0]
    rows = []
    for cand in base:
        m = cand.config.model
        if m.harmonics not in ((2, 3), (1, 2, 3)):
            continue
        for prof in PROFILES:
            cfg = dataclasses.replace(
                cand.config,
                model=ModelSpec(profile=prof, harmonics=m.harmonics, ocon=m.ocon),
            ).validate()
            row = dict(
                file=path.name,
                profile=prof,
                harmonics="".join(map(str, m.harmonics)),
                ocon=m.ocon,
                status="ok",
                rms_truth=None,
                rms_data=None,
                two_component=None,
            )
            try:
                out = fit(
                    t,
                    mag,
                    payload["sig"],
                    payload["period"],
                    payload["epoch_jd"],
                    cfg,
                    data_id=path.stem,
                )
                idx = np.searchsorted(t, out.t_mjd)
                resid = out.pipeline.adopted.resid
                row["rms_truth"] = 1e3 * float(
                    np.sqrt(np.mean((out.mag - resid - truth[idx]) ** 2))
                )
                row["rms_data"] = 1e3 * float(np.sqrt(np.mean(resid**2)))
                row["two_component"] = bool(out.solution.stats["two_component"])
            except MikulasekError as exc:
                row["status"] = f"{type(exc).__name__}"
            rows.append(row)
            logger.warning(
                "%s p%d h%s oc=%s %s %s",
                path.name[-14:],
                prof,
                row["harmonics"],
                m.ocon,
                row["status"],
                row["rms_truth"],
            )
    return rows


def main() -> None:
    """Run the experiment and write ``profiles_vs_truth.md/csv``."""
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    filters = sys.argv[1:]
    with open(DEFAULT_RESULTS / "manifest.csv", newline="") as fh:
        manifest = [r for r in csv.DictReader(fh) if int(r["fail"]) == 0]
    rows = []
    for m in manifest:
        if filters and not any(f in m["path"] for f in filters):
            continue
        rows += run_file(GENERATOR / m["path"])
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "profiles_vs_truth.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    lines = []
    for name in dict.fromkeys(r["file"] for r in rows):
        lines += [
            f"### {name}",
            "",
            "| profile | best rms truth (mmag) | with | two-component | failed fits |",
            "|---|---|---|---|---|",
        ]
        for p in PROFILES:
            fr = [r for r in rows if r["file"] == name and r["profile"] == p]
            ok = [r for r in fr if r["rms_truth"] is not None]
            if not ok:
                lines.append(f"| {p} | - | - | - | {len(fr)} |")
                continue
            b = min(ok, key=lambda r: r["rms_truth"])
            lines.append(
                f"| {p} | {b['rms_truth']:.3f} | h{b['harmonics']}{' oc' if b['ocon'] else ''} | {b['two_component']} | {len(fr) - len(ok)} |"
            )
        lines.append("")
    (OUT / "profiles_vs_truth.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
