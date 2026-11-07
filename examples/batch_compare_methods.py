"""Compare the timings of the mixed-trend batch: Mikulasek, AP, WSAP and WSL.

Inputs (all already on disk, nothing is timed here):

* Mikulasek at the plateau windows, ``batch_plateau_mik/eclipses.csv``;
* Mikulasek with relaxed gates and 0.25 cycle windows, ``batch_relaxed_mik/eclipses.csv``;
* MAVKA AP, WSAP and WSL at the plateau windows, ``output/batch_mix_plateau/``;
* MAVKA AP, WSAP and WSL at the expert windows, ``output/batch_mix/tables/``.

Each accepted time is matched to the nearest pure-curve minimum of the same kind
(within a quarter period) and compared with it and with the PHOEBE conjunction.
The second G1 600 s file (B22, wrong period) is left out of the figures.
Eclipses not observed, or without a time, are not counted. ``n`` is the number
of times of each row, and the sets of eclipses are not forced to be equal.

Writes ``results/synthetic_scoring/batch_compare/``: ``summary.csv`` and the figures
``scatter_vs_pure.png``, ``bias_vs_conj.png``, ``bias_vs_pure.png``, and the same
scored on the eclipses timed by Mik relaxed, AP, WSAP and WSL together
(``summary_common.csv``, ``*_common.png``).

Run from ``examples``::

    python batch_compare_methods.py
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from batch_plateau_time import GENERATED, ROOT, manifest_rows, write_csv
from mikulasek_eb.synthetic_check import _array, _num, load_synthetic

logger = logging.getLogger(__name__)
DAY = 86400.0
SCORING = ROOT / "results" / "synthetic_scoring"
OUT = SCORING / "batch_compare"
MAVKA = Path("/home/voz/projects/UPJS/Mavka_user_folder/output")
STARS = ["G1", "G2", "G4", "G5", "G6", "G7", "G8"]
CADENCES = [120, 600, 1800]
LABELS = {
    "G1": "Well-detached",
    "G2": "Total eclipse",
    "G4": "Algol",
    "G5": "W UMa",
    "G6": "O'Connell",
    "G7": "Shallow secondary",
    "G8": "Eccentric",
}
# name -> (colour, line style, marker)
COMMON_SET = ("Mik relaxed", "AP plateau", "WSAP plateau", "WSL plateau")
STYLE = {
    "Mik plateau": ("red", "-", "o"),
    "Mik relaxed": ("red", "--", "s"),
    "AP plateau": ("tab:blue", "-", "o"),
    "WSAP plateau": ("tab:orange", "-", "o"),
    "WSL plateau": ("tab:green", "-", "o"),
    "AP expert": ("tab:blue", ":", "x"),
    "WSAP expert": ("tab:orange", ":", "x"),
    "WSL expert": ("tab:green", ":", "x"),
}


def truth(curve: Path, kind: str) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Pure minima and conjunctions of one eclipse kind, absolute JD.

    Args:
        curve (Path): Light-curve file with truth metadata.
        kind (str): ``primary`` or ``secondary``.

    Returns:
        tuple: Pure minima, conjunctions, true period (d) and time origin (JD).

    Raises:
        RuntimeError: If the file has no eclipse of this kind.
    """
    synthetic = load_synthetic(curve)
    meta = synthetic.meta
    kinds = [k.strip() for k in str(meta["eclipse_kind"]).split(",")]
    mask = np.array([k == kind for k in kinds])
    if not mask.any():
        raise RuntimeError(f"{curve.name}: no {kind} in the truth metadata")
    origin = float(synthetic.timesys["timeorigin"])
    return (
        origin + _array(meta, "eclipse_T_min_pure_days")[mask],
        origin + _array(meta, "eclipse_T_conj_days")[mask],
        _num(meta, "P_true_days"),
        origin,
    )


def read_times() -> dict[tuple[str, str, str], list[float]]:
    """All accepted times by (method name, file, kind), absolute JD.

    Returns:
        dict: Times in days of JD, for every method of the comparison.
    """
    origin = {}
    times: dict[tuple[str, str, str], list[float]] = {}
    for row in manifest_rows():
        origin[Path(row["path"]).name] = float(
            load_synthetic(GENERATED / row["path"]).timesys["timeorigin"]
        )
    for name, folder in (
        ("Mik plateau", "batch_plateau_mik"),
        ("Mik relaxed", "batch_relaxed_mik"),
    ):
        with (SCORING / folder / "eclipses.csv").open(newline="") as handle:
            for r in csv.DictReader(handle):
                if r["T_jd"]:
                    times.setdefault((name, r["file"], r["kind"]), []).append(
                        float(r["T_jd"])
                    )
    for tag, folder in (
        ("plateau", MAVKA / "batch_mix_plateau"),
        ("expert", MAVKA / "batch_mix" / "tables"),
    ):
        for method in ("AP", "WSAP", "WSL"):
            with (folder / f"{method}.csv").open(newline="") as handle:
                for r in csv.DictReader(handle):
                    if r["status"] == "accepted" and r["time_d"]:
                        t = origin[r["file"]] + float(r["time_d"])
                        times.setdefault(
                            (f"{method} {tag}", r["file"], r["kind"]), []
                        ).append(t)
    return times


def matched(times: list[float], pure: np.ndarray, period: float) -> dict[int, float]:
    """Map each accepted time to the index of its pure minimum.

    Args:
        times (list[float]): Accepted times, absolute JD.
        pure (ndarray): Pure-curve minima, JD.
        period (float): True period, days.

    Returns:
        dict[int, float]: Pure-minimum index to time, within a quarter period.
    """
    out: dict[int, float] = {}
    for t in times:
        k = int(np.argmin(np.abs(pure - t)))
        if abs(t - pure[k]) < period / 4:
            out[k] = t
    return out


def score(
    times: list[float],
    pure: np.ndarray,
    conj: np.ndarray,
    period: float,
    only: set | None = None,
) -> dict:
    """Offsets of measured times from the truth.

    Args:
        times (list[float]): Accepted times, absolute JD.
        pure (ndarray): Pure-curve minima, JD.
        conj (ndarray): Conjunctions, JD.
        period (float): True period, days.
        only (set | None): Pure-minimum indices to keep, or None for all.

    Returns:
        dict: n, bias and scatter in seconds. Empty values when nothing matched.
    """
    out = {
        "n": 0,
        "bias_vs_pure_s": "",
        "bias_vs_conj_s": "",
        "scatter_vs_pure_s": "",
        "rmse_vs_pure_s": "",
    }
    if not times:
        return out
    pairs = matched(times, pure, period)
    idx = sorted(k for k in pairs if only is None or k in only)
    if not idx:
        return out
    t = np.array([pairs[k] for k in idx])
    dp = (t - pure[idx]) * DAY
    dc = (t - conj[idx]) * DAY
    out.update(
        n=len(idx),
        bias_vs_pure_s=float(np.mean(dp)),
        bias_vs_conj_s=float(np.mean(dc)),
        scatter_vs_pure_s=float(np.std(dp)),
        rmse_vs_pure_s=float(np.sqrt(np.mean(dp**2))),
    )
    return out


def build_table(common: bool = False) -> list[dict]:
    """Score every method on every file and eclipse kind.

    Args:
        common (bool): If True, keep only the eclipses timed by Mik relaxed, AP, WSAP
            and WSL at the plateau windows, so every method scores the same eclipses.

    Returns:
        list[dict]: One row per star, cadence, kind and method.
    """
    times = read_times()
    rows = []
    for m in manifest_rows():
        curve = GENERATED / m["path"]
        if m["recipe_id"] == "B22":
            continue
        for kind in ("primary", "secondary"):
            pure, conj, period, _ = truth(curve, kind)
            only = None
            if common:
                only = set.intersection(
                    *(
                        set(matched(times.get((n, curve.name, kind), []), pure, period))
                        for n in COMMON_SET
                    )
                )
            for name in STYLE:
                rec = score(
                    times.get((name, curve.name, kind), []), pure, conj, period, only
                )
                rec.update(
                    star=m["physics_id"],
                    cadence_s=int(m["cadence_s"]),
                    kind=kind,
                    method=name,
                )
                rows.append(rec)
    return rows


def value(rows: list[dict], star: str, cadence: int, kind: str, method: str, key: str):
    """One table value or NaN.

    Args:
        rows (list[dict]): Table from ``build_table``.
        star (str): Physics id.
        cadence (int): Cadence in seconds.
        kind (str): Eclipse kind.
        method (str): Method name.
        key (str): Column.

    Returns:
        float: The value, NaN where there is no time.
    """
    for r in rows:
        if (r["star"], r["cadence_s"], r["kind"], r["method"]) == (
            star,
            cadence,
            kind,
            method,
        ):
            return float(r[key]) if r[key] != "" else np.nan
    raise RuntimeError(f"no row {star} {cadence} {kind} {method}")


def draw(rows: list[dict], key: str, ylabel: str, path: Path, log: bool) -> None:
    """Draw one metric against cadence for every star and kind.

    Args:
        rows (list[dict]): Table.
        key (str): Metric column.
        ylabel (str): Axis label.
        path (Path): Output PNG.
        log (bool): Log axis for positive metrics, symlog otherwise.
    """
    fig, axes = plt.subplots(
        2, len(STARS), figsize=(3.0 * len(STARS), 6.4), sharex=True
    )
    for i, kind in enumerate(("primary", "secondary")):
        for j, star in enumerate(STARS):
            ax = axes[i, j]
            for name, (colour, ls, marker) in STYLE.items():
                y = [value(rows, star, c, kind, name, key) for c in CADENCES]
                ax.plot(
                    range(3),
                    y,
                    color=colour,
                    ls=ls,
                    marker=marker,
                    lw=1.4 if "expert" not in name else 0.9,
                    ms=4,
                    label=name,
                )
            ax.set_xticks(range(3), [str(c) for c in CADENCES])
            if log:
                ax.set_yscale("log")
            else:
                ax.set_yscale("symlog", linthresh=2)
                ax.axhline(0, color="0.6", lw=0.6)
            ax.grid(alpha=0.3)
            if i == 0:
                ax.set_title(f"{LABELS[star]} ({star})", fontsize=9)
            if j == 0:
                ax.set_ylabel(f"{kind}: {ylabel}")
            if i == 1:
                ax.set_xlabel("cadence (s)")
    axes[0, 0].legend(fontsize=6, loc="best")
    fig.suptitle(
        f"{ylabel}. Solid: plateau windows. Dashed: Mikulasek relaxed. Dotted: MAVKA expert windows."
    )
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main() -> None:
    """Write the comparison table and the figures."""
    OUT.mkdir(parents=True, exist_ok=True)
    rows = build_table()
    write_csv(OUT / "summary.csv", rows)
    draw(
        rows,
        "scatter_vs_pure_s",
        "scatter of time about the pure minimum (s)",
        OUT / "scatter_vs_pure.png",
        True,
    )
    draw(
        rows,
        "bias_vs_conj_s",
        "mean time minus conjunction (s)",
        OUT / "bias_vs_conj.png",
        False,
    )
    draw(
        rows,
        "bias_vs_pure_s",
        "mean time minus pure minimum (s)",
        OUT / "bias_vs_pure.png",
        False,
    )
    paired = build_table(common=True)
    write_csv(OUT / "summary_common.csv", paired)
    draw(
        paired,
        "scatter_vs_pure_s",
        "scatter about the pure minimum, common eclipses (s)",
        OUT / "scatter_vs_pure_common.png",
        True,
    )
    draw(
        paired,
        "bias_vs_conj_s",
        "mean time minus conjunction, common eclipses (s)",
        OUT / "bias_vs_conj_common.png",
        False,
    )
    logger.info("wrote %s", OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
