"""Offset against spot longitude for the nine 120 s curves.

Two panels, primary and secondary. Each panel shows the bisector median of the
pure curve, the Mikulášek time, the AP time and the WSAP time. All four are
seconds from conjunction. Positive means later than conjunction. The centroid
is not drawn. WSL is not read. The numbers come from tables already written.
Nothing is refitted. AP and WSAP use the centre ± 0.0176 cycles on the primary
and ± 0.0178 cycles on the secondary. Mikulášek uses the centre ± 0.12 cycles.

Run from the repository root::

    MPLBACKEND=Agg python examples/spot_120_offset_figure.py
    MPLBACKEND=Agg python examples/spot_120_offset_figure.py --set m11

``--set n01`` (the default) reads the 1 percent timings and writes
``results/synthetic_scoring/spot_120_select/offset_vs_longitude.png``.
``--set m11`` reads the Tmag 11 timings and writes the same names under
``results/synthetic_scoring/spot_120_select_m11/``.
"""

from __future__ import annotations

import argparse
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
HALF_PHASE = {"primary": 0.0176, "secondary": 0.0178}
GENERATOR = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MAVKA_ROOT = Path("/home/voz/projects/UPJS/Mavka_user_folder/output")
SETS = {
    "n01": (
        ROOT / "results" / "synthetic_scoring" / "spot_120_select",
        MAVKA_ROOT / "detached_120" / "tables",
        GENERATOR / "results" / "n01" / "detached",
        "1 percent of the flux",
    ),
    "m11": (
        ROOT / "results" / "synthetic_scoring" / "spot_120_select_m11",
        MAVKA_ROOT / "detached_120_m11" / "tables",
        GENERATOR / "results" / "m11" / "detached",
        "Tmag 11 (2.6 mmag at 120 s)",
    ),
}
TIMING = SETS["n01"][0] / "timing_120_w0.12.csv"
MAVKA = SETS["n01"][1]
CURVES = SETS["n01"][2]
NOISE = SETS["n01"][3]
OUT = SETS["n01"][0]
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


def method_offsets(method: str) -> dict[tuple[str, str], dict]:
    """Mean time minus conjunction for one MAVKA method.

    ``time_d`` and ``eclipse_T_conj_days`` are both days on the file clock.
    A row is used only when its status is ``accepted``. The half-width on each
    row must be the expert value for that eclipse kind. WSL is not read.

    Args:
        method (str): ``AP`` or ``WSAP``.

    Returns:
        dict: ``(physics, kind)`` to the mean offset in seconds, the scatter,
        and the number of eclipses. Keys are prefixed with the method name
        in lower case.

    Raises:
        RuntimeError: If the table is the wrong method, a window is not the
        expert half-width, a time is empty, or a window centre is not the
        conjunction of that eclipse.
    """
    grouped: dict[str, list[dict]] = {}
    for row in read_rows(MAVKA / f"{method}.csv"):
        if row["method"] != method:
            raise RuntimeError(f"expected {method}, found {row['method']}")
        grouped.setdefault(row["file"], []).append(row)
    if set(map(physics_of, grouped)) != set(ORDER):
        raise RuntimeError(f"{method} files are {sorted(map(physics_of, grouped))}")
    prefix = method.lower()
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
                raise RuntimeError(f"{physics} {kind} has no accepted {method} time")
            windows = {float(r["half_phase_window"]) for r in accepted}
            if windows != {HALF_PHASE[kind]}:
                raise RuntimeError(
                    f"{physics} {kind} {method} half-width is {windows}, "
                    f"expected {HALF_PHASE[kind]}"
                )
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
                    f"{physics} {kind}: a {method} time is far from conjunction"
                )
            seconds = delta * DAY
            out[(physics, kind)] = {
                f"{prefix}_minus_conj_s": float(np.mean(seconds)),
                f"{prefix}_scatter_s": float(np.std(seconds)),
                f"{prefix}_n": int(len(seconds)),
            }
            logger.info(
                "%s %s %s n=%d vs conj %+.2f s",
                method,
                physics,
                kind,
                len(seconds),
                float(np.mean(seconds)),
            )
    return out


def series() -> list[dict]:
    """Join the bisector median, the Mikulášek offset, AP and WSAP.

    Returns:
        list[dict]: One row per system and eclipse, with the plotted values.
    """
    pure = read_rows(PURE)
    timing = read_rows(TIMING)
    ap = method_offsets("AP")
    wsap = method_offsets("WSAP")
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
                    **ap[(physics, kind)],
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
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), sharey=True)
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
            [r["ap_minus_conj_s"] for r in part],
            yerr=[r["ap_scatter_s"] for r in part],
            fmt="D-",
            color="#6b3fa0",
            capsize=3,
            label="AP",
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
        f"120 s, {NOISE}. Positive means later than conjunction. "
        "The bar is the eclipse scatter.\n"
        "Window is the centre ± half-width. AP and WSAP: primary ± 0.0176, "
        "secondary ± 0.0178 cycles. Mikulášek: ± 0.12 cycles (full width 0.24)."
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    logger.info("figure %s", path)


def main() -> None:
    """Read the two tables and write the longitude figure."""
    global TIMING, MAVKA, CURVES, NOISE, OUT
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--set",
        choices=tuple(SETS),
        default="n01",
        help="n01 is the 1 percent set. m11 is the Tmag 11 set.",
    )
    args = parser.parse_args()
    OUT, MAVKA, CURVES, NOISE = SETS[args.set]
    TIMING = OUT / "timing_120_w0.12.csv"
    logger.info("noise set %s, curves %s, mavka %s", args.set, CURVES, MAVKA)
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
