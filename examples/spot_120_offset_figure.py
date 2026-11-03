"""Offset against spot longitude for the nine 120 s curves.

Two panels, primary and secondary. Each panel shows the bisector median of the
pure curve, the Mikulášek time, and the WSAP time. All three are seconds from
conjunction. Positive means later than conjunction. The centroid is not drawn.
WSL is not read. The numbers come from tables already written. Nothing is refitted.

Run from the repository root::

    MPLBACKEND=Agg python examples/spot_120_offset_figure.py

The figure is ``results/synthetic_scoring/spot_120_select/offset_vs_longitude.png``.
"""

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mikulasek_eb.synthetic_check import _array, _num, load_synthetic

logger = logging.getLogger(__name__)
DAY = 86400.0
ROOT = Path(__file__).resolve().parents[1]
PURE = ROOT / "results" / "synthetic_scoring" / "spot_pure" / "pure_minima.csv"
TIMING = (
    ROOT / "results" / "synthetic_scoring" / "spot_120_select" / "timing_120_w0.12.csv"
)
WSAP = Path(
    "/home/voz/projects/UPJS/Mavka_user_folder/output/detached_120/tables/WSAP.csv"
)
CURVES = Path("/home/voz/projects/UPJS/EB_synthetic_lc/results/detached")
OUT = ROOT / "results" / "synthetic_scoring" / "spot_120_select"
ORDER = (
    "SP_none",
    "SP_L000",
    "SP_L045",
    "SP_L090",
    "SP_L135",
    "SP_L180",
    "SP_L225",
    "SP_L270",
    "SP_L315",
)
LABELS = ("none", "0", "45", "90", "135", "180", "225", "270", "315")


def read_rows(path: Path) -> list[dict]:
    """Read a CSV table.

    Args:
        path (Path): Existing table.

    Returns:
        list[dict]: One dict per row.

    Raises:
        FileNotFoundError: If the table is missing.
    """
    if not path.is_file():
        raise FileNotFoundError(path)
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def one(rows: list[dict], physics: str, kind: str) -> dict:
    """The single row for one system and one eclipse.

    Args:
        rows (list[dict]): Table rows.
        physics (str): ``SP_none`` or ``SP_L045`` and so on.
        kind (str): ``primary`` or ``secondary``.

    Returns:
        dict: That row.

    Raises:
        RuntimeError: If the row is missing or repeated.
    """
    found = [r for r in rows if r["physics"] == physics and r["kind"] == kind]
    if len(found) != 1:
        raise RuntimeError(f"{physics} {kind}: found {len(found)} rows")
    return found[0]


def physics_of(filename: str) -> str:
    """Physics id carried in one light-curve file name.

    Args:
        filename (str): VOTable file name.

    Returns:
        str: One entry of ``ORDER``.

    Raises:
        RuntimeError: If the name matches none or more than one id.
    """
    found = [physics for physics in ORDER if physics in filename]
    if len(found) != 1:
        raise RuntimeError(f"{filename} matches {found}")
    return found[0]


def wsap_offsets() -> dict[tuple[str, str], dict]:
    """Mean WSAP time minus conjunction, per system and eclipse kind.

    ``time_d`` and ``eclipse_T_conj_days`` are both days on the file clock.
    A row is used only when its status is ``accepted``. WSL is not read.

    Returns:
        dict: ``(physics, kind)`` to the mean offset in seconds, the scatter,
        and the number of eclipses.

    Raises:
        RuntimeError: If the table is not WSAP, a time is empty, or a window
        centre is not the conjunction of that eclipse.
    """
    grouped: dict[str, list[dict]] = {}
    for row in read_rows(WSAP):
        if row["method"] != "WSAP":
            raise RuntimeError(f"expected WSAP, found {row['method']}")
        grouped.setdefault(row["file"], []).append(row)
    if set(map(physics_of, grouped)) != set(ORDER):
        raise RuntimeError(f"WSAP files are {sorted(map(physics_of, grouped))}")
    out = {}
    for filename, group in grouped.items():
        physics = physics_of(filename)
        meta = load_synthetic(CURVES / filename).meta
        kinds = [k.strip() for k in str(meta["eclipse_kind"]).split(",")]
        tconj = _array(meta, "eclipse_T_conj_days")
        period = _num(meta, "P_true_days")
        for kind in ("primary", "secondary"):
            accepted = [
                r for r in group if r["kind"] == kind and r["status"] == "accepted"
            ]
            if not accepted or any(r["time_d"] == "" for r in accepted):
                raise RuntimeError(f"{physics} {kind} has no accepted WSAP time")
            conj = tconj[np.array([k == kind for k in kinds])]
            measured = np.array([float(r["time_d"]) for r in accepted])
            predicted = np.array([float(r["t_pred_d"]) for r in accepted])
            nearest = np.argmin(np.abs(measured[:, None] - conj[None, :]), axis=1)
            if np.max(np.abs(predicted - conj[nearest])) > 1e-3:
                raise RuntimeError(
                    f"{physics} {kind}: window centre is not conjunction"
                )
            delta = measured - conj[nearest]
            if np.any(np.abs(delta) > period / 4):
                raise RuntimeError(
                    f"{physics} {kind}: a WSAP time is far from conjunction"
                )
            seconds = delta * DAY
            out[(physics, kind)] = dict(
                wsap_minus_conj_s=float(np.mean(seconds)),
                wsap_scatter_s=float(np.std(seconds)),
                wsap_n=int(len(seconds)),
            )
            logger.info(
                "WSAP %s %s n=%d vs conj %+.2f s",
                physics,
                kind,
                len(seconds),
                float(np.mean(seconds)),
            )
    return out


def series() -> list[dict]:
    """Join the bisector median, the Mikulášek offset and the WSAP offset.

    Returns:
        list[dict]: One row per system and eclipse, with the plotted values.
    """
    pure = read_rows(PURE)
    timing = read_rows(TIMING)
    wsap = wsap_offsets()
    rows = []
    for physics in ORDER:
        for kind in ("primary", "secondary"):
            bisector = one(pure, physics, kind)
            fit = one(timing, physics, kind)
            if fit["status"] != "ok" or fit["bias_vs_conj_s"] == "":
                raise RuntimeError(
                    f"{physics} {kind} has no Mikulášek time: {fit['status']}"
                )
            rows.append(
                dict(
                    physics=physics,
                    kind=kind,
                    bisector_minus_conj_s=float(
                        bisector["bisector_median_minus_conj_s"]
                    ),
                    mikulasek_minus_conj_s=float(fit["bias_vs_conj_s"]),
                    mikulasek_scatter_s=float(fit["scatter_vs_pure_s"]),
                    **wsap[(physics, kind)],
                )
            )
    return rows


def draw(rows: list[dict], path: Path) -> None:
    """Save the two-panel figure.

    Args:
        rows (list[dict]): Plotted values from ``series``.
        path (Path): PNG to write.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    x = list(range(len(ORDER)))
    for ax, kind in zip(axes, ("primary", "secondary")):
        part = [r for r in rows if r["kind"] == kind]
        ax.axhline(0.0, color="0.4", lw=0.6)
        ax.plot(
            x,
            [r["bisector_minus_conj_s"] for r in part],
            "o-",
            color="#1f4e79",
            label="bisector median",
        )
        ax.errorbar(
            x,
            [r["mikulasek_minus_conj_s"] for r in part],
            yerr=[r["mikulasek_scatter_s"] for r in part],
            fmt="s-",
            color="#ff2e2e",
            capsize=3,
            label="Mikulášek",
        )
        ax.errorbar(
            x,
            [r["wsap_minus_conj_s"] for r in part],
            yerr=[r["wsap_scatter_s"] for r in part],
            fmt="^-",
            color="#2e7d32",
            capsize=3,
            label="WSAP",
        )
        ax.set_title(kind.capitalize())
        ax.set_xlabel("Spot longitude (degrees)")
        ax.set_xticks(x, LABELS)
    axes[0].set_ylabel("Offset from conjunction (s)")
    axes[0].legend()
    fig.suptitle(
        "120 s. Positive means later than conjunction. The bar is the eclipse scatter. "
        "WSAP full phase width 0.058."
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    logger.info("figure %s", path)


def main() -> None:
    """Read the two tables and write the longitude figure."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "offset_command.txt").write_text(" ".join(sys.argv) + "\n")
    rows = series()
    fields = list(rows[0])
    table = OUT / "offset_vs_longitude.csv"
    with open(table, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fields)
        writer.writeheader()
        writer.writerows(rows)
    logger.info("wrote %s", table)
    draw(rows, OUT / "offset_vs_longitude.png")


if __name__ == "__main__":
    main()
