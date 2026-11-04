"""Fit every degraded spot curve and report the fit status.

The status (``ok``, ``bound-active``, ``unresolved``) comes from ``mikulasek_eb.validity``.
The sweep is the false-positive check of that rule: curves that fit acceptably at 120 s and
600 s must stay ``ok``, and the 1800 s curves that ran away must not.

Run from the repository root::

    MPLBACKEND=Agg python experiments/runaway_1800/status_sweep.py --cadence 120 600 1800
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import logging
from pathlib import Path

import mikulasek_eb as eb
from mikulasek_eb.candidates import propose_candidates
from mikulasek_eb.features import extract_features
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable

logger = logging.getLogger(__name__)
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFEST = GENERATED / "results" / "n01" / "manifest_spotsweep.csv"
OUT = Path(__file__).resolve().parent / "out" / "status_sweep.csv"


def config_for(payload: dict, t_mjd) -> eb.Config:
    """Profile 1, H2, O'Connell on, straight baseline, this curve's exposure and gaps.

    Args:
        payload (dict): Loaded light curve.
        t_mjd (ndarray): Times (MJD).

    Returns:
        Config: Structure of the degraded spot run.
    """
    features = extract_features(
        t_mjd, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    own = propose_candidates(features)[0].config.data
    return eb.Config(
        eb.ModelSpec(profile=1, harmonics=(2,), ocon=True),
        dataclasses.replace(own, baseline_order=1),
        eb.FitPolicy(allow_two_component=True),
        eb.TimingRecipe(window=0.12),
    ).validate()


def main() -> None:
    """Fit all files of the requested cadences and write the status table."""
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--cadence", type=int, nargs="+", required=True)
    args = parser.parse_args()
    with open(MANIFEST) as fh:
        rows = [
            r
            for r in csv.DictReader(fh)
            if int(r["cadence_s"]) in args.cadence and int(r["fail"]) == 0
        ]
    table = []
    for row in rows:
        path = GENERATED / row["path"]
        payload = load_tess_votable(path)
        t = absolute_jd_to_mjd(payload["t"])
        out = eb.fit(
            t,
            payload["mag"],
            payload["sig"],
            payload["period"],
            payload["epoch_jd"],
            config_for(payload, t),
            data_id=path.stem,
        )
        stats = out.solution.stats
        table.append(
            dict(
                file=path.name,
                cadence_s=row["cadence_s"],
                physics=row["physics_id"],
                status=stats["status"],
                rms_mmag=1e3 * stats["rms_mag"],
                reasons=" | ".join(stats["status_reasons"]),
                warnings=" | ".join(stats["status_warnings"]),
            )
        )
        logger.warning(
            "%s: %s rms=%.2f mmag", path.name, stats["status"], 1e3 * stats["rms_mag"]
        )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    logger.warning("wrote %s", OUT)


if __name__ == "__main__":
    main()
