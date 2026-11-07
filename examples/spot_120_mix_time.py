"""Time the nine mixed-trend 120 s spot curves from the stored solutions.

No fit is run. Each curve uses the solution already chosen for it by the full
comparison in ``results/synthetic_scoring/spot_mix_full/``. The timing window
is the centre ± 0.12 cycles, the window of the spotted experiment. The stored
solution does not keep the fit's clip mask. Four of the nine fits had dropped
1 to 9 points out of 18720, so every finite row of the file is timed.

Run from the repository root::

    MPLBACKEND=Agg python examples/spot_120_mix_time.py

Writes ``results/synthetic_scoring/spot_120_mix_time/timing_120_w0.12.csv``.
The 1 percent timings and the additive Tmag 11 timings are not touched.
"""

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path

import mikulasek_eb as eb
import spot_120_select as spot
from mikulasek_eb.solution import ModelSolution

logger = logging.getLogger(__name__)
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFEST = GENERATED / "results" / "mix" / "manifest_spotsweep.csv"
SOLUTIONS = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "synthetic_scoring"
    / "spot_mix_full"
    / "solutions"
)
CHOSEN = SOLUTIONS.parent / "chosen.csv"
OUT = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "synthetic_scoring"
    / "spot_120_mix_time"
)
FIELDS = [
    "file",
    "physics",
    "candidate",
    "harmonics",
    "ocon",
    "baseline_order",
    "rms_mmag",
    "two_component",
    "window",
    "status",
    "reason",
    "kind",
    "n",
    "n_rejected",
    "bias_vs_conj_s",
    "bias_vs_pure_s",
    "scatter_vs_pure_s",
    "median_sigma_s",
    "honesty_vs_pure",
]


class StoredFit:
    """The arrays and solution ``spot_120_select.offsets`` reads.

    Attributes:
        t_mjd (ndarray): Times (MJD).
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        solution (ModelSolution): Stored template.
    """

    def __init__(self, t_mjd, mag, sig, solution: ModelSolution) -> None:
        self.t_mjd = t_mjd
        self.mag = mag
        self.sig = sig
        self.solution = solution


def rows_120() -> list[dict]:
    """Accepted 120 s rows of the mixed-trend spot sweep.

    Returns:
        list[dict]: One manifest row per light curve, spot-free first.

    Raises:
        RuntimeError: If the number of accepted 120 s curves is not 9.
    """
    with open(MANIFEST, newline="") as fh:
        rows = [
            r
            for r in csv.DictReader(fh)
            if int(r["cadence_s"]) == 120 and int(r["fail"]) == 0
        ]
    rows.sort(key=lambda r: (r["physics_id"] != "SP_none", r["physics_id"]))
    if len(rows) != 9:
        raise RuntimeError(f"expected 9 accepted 120 s curves, found {len(rows)}")
    return rows


def chosen_names() -> dict[str, str]:
    """Candidate name of each stored 120 s solution.

    Returns:
        dict[str, str]: File name to the chosen candidate.
    """
    with open(CHOSEN, newline="") as fh:
        return {
            r["file"]: r["chosen"]
            for r in csv.DictReader(fh)
            if int(r["cadence_s"]) == 120
        }


def one_curve(manifest_row: dict, names: dict[str, str]) -> list[dict]:
    """Time one stored solution.

    Args:
        manifest_row (dict): One accepted 120 s row.
        names (dict[str, str]): Chosen candidate name by file name.

    Returns:
        list[dict]: Primary and secondary timing rows.

    Raises:
        FileNotFoundError: If the stored solution is missing.
        RuntimeError: If the stored solution is not status ``ok``.
    """
    curve = GENERATED / manifest_row["path"]
    solution_path = SOLUTIONS / f"{curve.stem}.json"
    if not solution_path.is_file():
        raise FileNotFoundError(solution_path)
    solution = ModelSolution.load(solution_path)
    status = solution.stats["status"]
    if status != "ok":
        raise RuntimeError(f"{curve.name}: stored solution status is {status}")
    payload = eb.load_tess_votable(curve)
    stored = StoredFit(
        eb.absolute_jd_to_mjd(payload["t"]),
        payload["mag"],
        payload["sig"],
        solution,
    )
    shared = spot.identity(curve, manifest_row["physics_id"], stored, status, "")
    shared["candidate"] = names[curve.name]
    timed = spot.offsets(curve, stored)
    for row in timed:
        row.update(shared)
    return timed


def main() -> None:
    """Time the nine stored solutions and write the offset table."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "command.txt").write_text(" ".join(sys.argv) + "\n")
    names = chosen_names()
    timings: list[dict] = []
    for manifest_row in rows_120():
        timings.extend(one_curve(manifest_row, names))
    spot.write_table(OUT / "timing_120_w0.12.csv", timings, FIELDS)


if __name__ == "__main__":
    main()
