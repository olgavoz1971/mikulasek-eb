"""Definitions of the minimum on the pure spot-sweep curves.

Reads each one-orbit parent (no noise, no trend) and measures the primary and the secondary
with ``measure_curve_minimum``. The window is 0.12 of the period, as in the generator. The
median of the bisectors is compared with the generator's ``T_conj - T_min_pure`` taken from
the 120 s file of the same system. Writes ``results/synthetic_scoring/spot_pure/``. Run::

    python examples/spot_pure_minima.py
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mikulasek_eb import BISECTOR_FRACTIONS, measure_curve_minimum
from mikulasek_eb.synthetic_check import _num, load_synthetic

logger = logging.getLogger(__name__)
DAY = 86400.0
WINDOW_FRACTION = 0.12
GENERATOR = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
PARENTS = GENERATOR / "results" / "clean" / "detached"
SPOILED = GENERATOR / "results" / "n01" / "detached"
OUT = (
    Path(__file__).resolve().parents[1] / "results" / "synthetic_scoring" / "spot_pure"
)


def parents() -> list[Path]:
    """One-orbit parents of the spot sweep.

    Returns:
        list[Path]: Parent files, the spot-free one first and then by longitude.
    """
    files = [
        p
        for p in PARENTS.glob("*SP_*.vot")
        if "_c120" not in p.name and "_c600" not in p.name and "_c1800" not in p.name
    ]

    def key(path: Path) -> tuple:
        """Sort the spot-free parent before the longitudes."""
        name = path.stem
        if name.endswith("SP_none"):
            return (0, -1.0)
        return (1, float(name.split("SP_L")[-1]))

    return sorted(files, key=key)


def generator_offset(physics: str, kind: str) -> float:
    """``T_conj - T_min_pure`` in seconds from the 120 s file.

    Args:
        physics (str): ``SP_none`` or ``SP_L045`` and so on.
        kind (str): ``primary`` or ``secondary``.

    Returns:
        float: Offset in seconds.
    """
    matches = list(SPOILED.glob(f"*{physics}_{physics}_120_c120.vot"))
    if len(matches) != 1:
        raise FileNotFoundError(f"{physics} 120 s file: found {len(matches)}")
    meta = load_synthetic(matches[0]).meta
    return float(_num(meta, f"conj_minus_min_pure_{kind}_s"))


def main() -> None:
    """Measure every parent and write the table and the figure."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in parents():
        data = load_synthetic(path, with_errors=False)
        if data.problems:
            raise RuntimeError(f"{path.name}: {data.problems}")
        period = float(data.meta["period"])
        physics = (
            "SP_none"
            if path.stem.endswith("SP_none")
            else "SP_L" + path.stem.split("SP_L")[-1]
        )
        longitude = data.meta.get("spot01_long_deg")
        for kind, guess in (("primary", 0.0), ("secondary", 0.5 * period)):
            measured = measure_curve_minimum(
                data.time,
                data.flux,
                guess,
                WINDOW_FRACTION * period,
                kind="flux",
                period=period,
            )
            file_offset = generator_offset(physics, kind)
            ours = (measured.t_bisector_median - guess) * DAY
            row = {
                "physics": physics,
                "longitude_deg": "" if longitude is None else float(longitude),
                "kind": kind,
                "lowest_minus_conj_s": (measured.t_lowest - guess) * DAY,
                "centroid_minus_conj_s": (measured.t_centroid - guess) * DAY,
                "bisector_median_minus_conj_s": ours,
                "bisector_spread_s": (
                    max(measured.t_bisector.values())
                    - min(measured.t_bisector.values())
                )
                * DAY,
                "generator_conj_minus_min_s": file_offset,
                "median_minus_generator_s": ours + file_offset,
            }
            for frac, time in measured.t_bisector.items():
                row[f"bisector_{frac:.2f}_s"] = (time - guess) * DAY
            rows.append(row)
            logger.info(
                "%s %s median %+.1f s, generator %+.1f s, difference %+.2f s",
                physics,
                kind,
                ours,
                -file_offset,
                ours + file_offset,
            )
    fields = list(rows[0])
    with open(OUT / "pure_minima.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fields)
        writer.writeheader()
        writer.writerows(rows)
    _figure(rows)
    logger.info("wrote %s", OUT)


def _figure(rows: list[dict]) -> None:
    """Offset of each definition against the spot longitude.

    Args:
        rows (list[dict]): One row per system and eclipse kind.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, kind in zip(axes, ("primary", "secondary")):
        part = [r for r in rows if r["kind"] == kind]
        x = [(-1 if r["longitude_deg"] == "" else r["longitude_deg"]) for r in part]
        for frac in BISECTOR_FRACTIONS:
            ax.plot(x, [r[f"bisector_{frac:.2f}_s"] for r in part], ":", color="0.6")
        ax.plot(
            x,
            [r["bisector_median_minus_conj_s"] for r in part],
            "o-",
            label="bisector median",
        )
        ax.plot(x, [r["centroid_minus_conj_s"] for r in part], "s-", label="centroid")
        ax.plot(x, [r["lowest_minus_conj_s"] for r in part], "^-", label="lowest point")
        ax.axhline(0, color="k", lw=0.4)
        ax.set_title(kind)
        ax.set_xlabel("spot longitude (deg); -1 is no spot")
        ax.set_xticks([-1, 0, 45, 90, 135, 180, 225, 270, 315])
    axes[0].set_ylabel("photometric minimum minus conjunction (s)")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "pure_minima.png", dpi=100)


if __name__ == "__main__":
    main()
