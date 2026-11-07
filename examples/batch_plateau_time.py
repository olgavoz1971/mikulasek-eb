"""Time the mixed-trend batch with Mikulasek at the plateau windows.

The window of each eclipse is the multiple of D adopted in the outline
(section "Timing experiment: how the windows were chosen"). At 120 s, D is
the width of that solution. At 600 s and 1800 s the same multiple is applied
to the D of the solution of that file. A product above 0.25 cycles is cut to
0.25, which is the widest half-window ``time`` accepts without overlapping
the neighbouring minimum, and the cut is recorded.

No fit is repeated. The outlier cut is off, as in the scan that chose the
windows. A rejected eclipse is a row with its reason.

Run from the repository root::

    python examples/batch_plateau_time.py
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import numpy as np

import mikulasek_eb as eb
from mikulasek_eb.solution import ModelSolution
from mikulasek_eb.synthetic_check import _array, _num, load_synthetic

logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]
GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFEST = GENERATED / "results" / "mix" / "manifest.csv"
FULL = ROOT / "results" / "synthetic_scoring" / "batch_mix_full"
OUT = ROOT / "results" / "synthetic_scoring" / "batch_plateau_mik"
DAY = 86400.0
MAX_HALF = 0.25
WIDTH_KEY = {"primary": "D1", "secondary": "D2"}
# Narrowest half-window, in units of D, within one standard error of the best
# rms on the 120 s Mikulasek scan. Steps are tenths of D.
MULTIPLE = {
    "G1": {"primary": 2.8, "secondary": 2.6},
    "G2": {"primary": 2.9, "secondary": 1.5},
    "G4": {"primary": 1.4, "secondary": 1.5},
    "G5": {"primary": 4.7, "secondary": 4.7},
    "G6": {"primary": 4.3, "secondary": 4.7},
    "G7": {"primary": 7.8, "secondary": 6.8},
    "G8": {"primary": 2.6, "secondary": 1.7},
}


def manifest_rows() -> list[dict]:
    """Accepted mixed-trend files, including the second G1 600 s curve.

    Returns:
        list[dict]: Manifest rows with ``fail`` 0 and a physics id G1 to G8.
    """
    with MANIFEST.open(newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle)
            if int(row["fail"]) == 0 and row["physics_id"] in MULTIPLE
        ]
    if len(rows) != 22:
        raise RuntimeError(f"expected 22 files, found {len(rows)}")
    return rows


def chosen_solution(curve: Path) -> tuple[ModelSolution, str]:
    """Load the chosen solution of one file.

    Args:
        curve (Path): Light-curve file.

    Returns:
        tuple[ModelSolution, str]: The solution and its candidate name.

    Raises:
        RuntimeError: If the stored solution is not status ``ok``.
    """
    with (FULL / "chosen.csv").open(newline="") as handle:
        name = next(
            row["chosen"] for row in csv.DictReader(handle) if row["file"] == curve.name
        )
    solution = ModelSolution.load(FULL / "solutions" / f"{curve.stem}__{name}.json")
    if solution.stats["status"] != "ok":
        raise RuntimeError(
            f"{curve.name}: stored solution is {solution.stats['status']}"
        )
    return solution, name


def half_width(physics: str, kind: str, d_cycles: float) -> tuple[float, float, bool]:
    """Half-window for one eclipse, capped at the non-overlap limit.

    Args:
        physics (str): Physics id.
        kind (str): ``primary`` or ``secondary``.
        d_cycles (float): Mikulasek D of this eclipse in this solution, cycles.

    Returns:
        tuple[float, float, bool]: Requested half-width, half-width passed to
        ``time``, and whether the cap was applied.
    """
    requested = MULTIPLE[physics][kind] * d_cycles
    used = min(requested, MAX_HALF)
    return requested, used, used < requested - 1e-12


def truth_offsets(curve: Path, kind: str, measured: np.ndarray) -> dict:
    """Mean offset of the measured times from conjunction and from the pure minimum.

    Args:
        curve (Path): Light curve, for its truth metadata.
        kind (str): ``primary`` or ``secondary``.
        measured (ndarray): Accepted times, absolute JD.

    Returns:
        dict: Bias and scatter in seconds. Empty when nothing matched.

    Raises:
        RuntimeError: If the truth metadata has no eclipse of this kind.
    """
    synthetic = load_synthetic(curve)
    meta = synthetic.meta
    kinds = [item.strip() for item in str(meta["eclipse_kind"]).split(",")]
    mask = np.array([item == kind for item in kinds])
    if not np.any(mask):
        raise RuntimeError(f"{curve.name}: no {kind} in the truth metadata")
    origin = float(synthetic.timesys["timeorigin"])
    pure = origin + _array(meta, "eclipse_T_min_pure_days")[mask]
    conj = origin + _array(meta, "eclipse_T_conj_days")[mask]
    period = _num(meta, "P_true_days")
    if len(measured) == 0:
        return {}
    nearest = np.argmin(np.abs(measured[:, None] - pure[None, :]), axis=1)
    keep = np.abs(measured - pure[nearest]) < period / 4
    if not np.any(keep):
        return {}
    return {
        "bias_vs_conj_s": f"{float(np.mean((measured - conj[nearest])[keep])) * DAY:.2f}",
        "bias_vs_pure_s": f"{float(np.mean((measured - pure[nearest])[keep])) * DAY:.2f}",
        "scatter_vs_pure_s": f"{float(np.std((measured - pure[nearest])[keep])) * DAY:.2f}",
    }


def time_file(row: dict) -> tuple[list[dict], list[dict]]:
    """Time one file at the plateau windows.

    Args:
        row (dict): Manifest row.

    Returns:
        tuple[list[dict], list[dict]]: One summary row per eclipse kind, and one
        row per eclipse.
    """
    curve = GENERATED / row["path"]
    solution, name = chosen_solution(curve)
    payload = eb.load_tess_votable(curve)
    t_mjd = eb.absolute_jd_to_mjd(np.asarray(payload["t"], float))
    epoch = float(payload["epoch_jd"])
    period = float(payload["period"])
    summary, points = [], []
    for kind in ("primary", "secondary"):
        d_cycles = float(solution.params[WIDTH_KEY[kind]])
        requested, used, capped = half_width(row["physics_id"], kind, d_cycles)
        recipe = eb.TimingRecipe(window=used, clip_outliers=False)
        out = eb.time(
            t_mjd, payload["mag"], payload["sig"], solution, recipe, which=(kind,)
        )
        timed = out.timings[kind]
        rejected = out.rejected[kind]
        oc = []
        for cycle, t_jd, sigma in zip(timed["E"], timed["T"], timed["sigma"]):
            oc_s = (float(t_jd) - (epoch + float(cycle) * period)) * DAY
            oc.append(oc_s)
            points.append(
                {
                    "file": curve.name,
                    "kind": kind,
                    "cycle": f"{float(cycle):.1f}",
                    "T_jd": f"{float(t_jd):.6f}",
                    "sigma_s": f"{float(sigma) * DAY:.3f}",
                    "oc_file_s": f"{oc_s:.3f}",
                    "reason": "",
                }
            )
        for item in rejected:
            points.append(
                {
                    "file": curve.name,
                    "kind": kind,
                    "cycle": f"{float(item['E']):.1f}",
                    "T_jd": "",
                    "sigma_s": "",
                    "oc_file_s": "",
                    "reason": item["reason"],
                }
            )
        measured = np.asarray(timed["T"], float) if len(timed["E"]) else np.array([])
        record = {
            "file": curve.name,
            "physics": row["physics_id"],
            "cadence_s": row["cadence_s"],
            "candidate": name,
            "kind": kind,
            "D_cycles": f"{d_cycles:.6f}",
            "multiple": f"{MULTIPLE[row['physics_id']][kind]:.1f}",
            "half_requested": f"{requested:.6f}",
            "half_used": f"{used:.6f}",
            "capped": str(capped),
            "n_timed": str(len(oc)),
            "n_rejected": str(len(rejected)),
            "mean_oc_s": f"{float(np.mean(oc)):.2f}" if oc else "",
            "rms_oc_s": f"{float(np.std(oc)):.2f}" if len(oc) >= 2 else "",
            "median_sigma_s": (
                f"{float(np.median(np.asarray(timed['sigma'], float) * DAY)):.2f}"
                if len(oc)
                else ""
            ),
        }
        record.update(truth_offsets(curve, kind, measured))
        summary.append(record)
        logger.info(
            "%s %s half %.4f%s n=%d rejected=%d vs conj %s s",
            curve.name,
            kind,
            used,
            " capped" if capped else "",
            len(oc),
            len(rejected),
            record.get("bias_vs_conj_s", ""),
        )
    return summary, points


def write_csv(path: Path, rows: list[dict]) -> None:
    """Write rows with the columns of the first row.

    Args:
        path (Path): Destination.
        rows (list[dict]): Table.
    """
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    """Time every accepted mixed-trend file and write the tables."""
    OUT.mkdir(parents=True, exist_ok=True)
    summary, points = [], []
    for row in manifest_rows():
        part_summary, part_points = time_file(row)
        summary.extend(part_summary)
        points.extend(part_points)
    write_csv(OUT / "summary.csv", summary)
    write_csv(OUT / "eclipses.csv", points)
    logger.info("wrote %s", OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
