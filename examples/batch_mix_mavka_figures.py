"""Figures of the first-batch MAVKA timings against the PHOEBE truth.

Offset is the accepted time minus the matched geometric conjunction of the same
kind, in seconds. The pure-curve photometric minimum is drawn as a black dash,
because on the spotted G6 it is not the conjunction. A Mikulasek column can be
added later from the same tables. The second G1 file at 600 s (B22) is left out
of the grid because its catalogue period differs from the true one.

Run from the repository root::

    MPLBACKEND=Agg python examples/batch_mix_mavka_figures.py

Writes ``offset_vs_star.png``, ``scatter_vs_star.png`` and ``figure_data.csv``
in ``results/synthetic_scoring/batch_mix_mavka/``.
"""

from __future__ import annotations

import csv
import logging
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import mad_std

from batch_mix_mavka_offsets import (
    DAY_S,
    GENERATED,
    METHODS,
    OUT,
    PHYSICS,
    TABLES,
    _manifest_rows,
    _match,
    _split,
)
from mikulasek_eb.synthetic_check import _array, _num, load_synthetic

logger = logging.getLogger(__name__)
CADENCES = ("120", "600", "1800")
COLOURS = {"AP": "tab:blue", "WSAP": "tab:orange", "WSL": "tab:green"}
SHIFT = {"AP": -0.22, "WSAP": 0.0, "WSL": 0.22}


def collect() -> list[dict]:
    """Collect median offset, robust scatter and counts per star, cadence, method, kind.

    Returns:
        list[dict]: One record per combination. Offsets are seconds from conjunction.
    """
    rows = _manifest_rows()
    paths = {Path(r["path"]).name: GENERATED / r["path"] for r in rows}
    recipe = {Path(r["path"]).name: r["recipe_id"] for r in rows}
    records: list[dict] = []
    for method in METHODS:
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        with (TABLES / f"{method}.csv").open(newline="") as handle:
            for record in csv.DictReader(handle):
                grouped[record["file"]].append(record)
        for name, group in grouped.items():
            if recipe[name] == "B22":
                continue
            meta = load_synthetic(paths[name]).meta
            kinds = np.array([k == "primary" for k in _split(meta["eclipse_kind"])])
            conj = _array(meta, "eclipse_T_conj_days")
            pure = _array(meta, "eclipse_T_min_pure_days")
            period = _num(meta, "P_true_days")
            for kind in ("primary", "secondary"):
                mask = kinds if kind == "primary" else ~kinds
                accepted = [
                    r for r in group if r["kind"] == kind and r["status"] == "accepted"
                ]
                observed = [
                    r
                    for r in group
                    if r["kind"] == kind and "eclipse_observed" not in r["reason"]
                ]
                measured = np.array([float(r["time_d"]) for r in accepted])
                offset = _match(measured, conj[mask], period)
                records.append(
                    {
                        "physics_id": str(meta["physics_id"]),
                        "cadence_s": group[0]["cadence_s"],
                        "method": method,
                        "kind": kind,
                        "n_observed": len(observed),
                        "n_accepted": len(accepted),
                        "median_s": float(np.median(offset)),
                        "mad_std_s": float(mad_std(offset)),
                        "pure_minus_conj_s": float(np.median(pure[mask] - conj[mask]))
                        * DAY_S,
                    }
                )
    return records


def draw(records: list[dict], field: str, ylabel: str, name: str, signed: bool) -> None:
    """Draw one 2 x 3 panel figure (kind by cadence) of ``field`` against star.

    Args:
        records (list[dict]): Output of ``collect``.
        field (str): Record key to plot, ``median_s`` or ``mad_std_s``.
        ylabel (str): Axis label.
        name (str): File name inside ``OUT``.
        signed (bool): Symmetric log axis with pure-minimum dashes when True,
            ordinary log axis when False.
    """
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, sharey=True)
    for i, kind in enumerate(("primary", "secondary")):
        for j, cadence in enumerate(CADENCES):
            ax = axes[i, j]
            for method in METHODS:
                for k, star in enumerate(PHYSICS):
                    hit = [
                        r
                        for r in records
                        if (r["physics_id"], r["cadence_s"], r["method"], r["kind"])
                        == (star, cadence, method, kind)
                    ]
                    if not hit:
                        continue
                    ax.plot(
                        k + SHIFT[method],
                        hit[0][field],
                        "o",
                        color=COLOURS[method],
                        label=method if (k == 0) else None,
                    )
                    if signed:
                        ax.plot(
                            [k - 0.35, k + 0.35],
                            [hit[0]["pure_minus_conj_s"]] * 2,
                            color="black",
                            lw=1,
                            label=(
                                "pure-curve minimum"
                                if (method == "AP" and k == 0)
                                else None
                            ),
                        )
            ax.set_title(f"{kind}, {cadence} s")
            ax.set_xticks(range(len(PHYSICS)), PHYSICS)
            ax.grid(alpha=0.3)
            if signed:
                ax.set_yscale("symlog", linthresh=10)
                ax.axhline(0, color="grey", lw=0.8)
            else:
                ax.set_yscale("log")
            if j == 0:
                ax.set_ylabel(ylabel)
    axes[0, 0].legend(fontsize=8, loc="best")
    fig.suptitle(
        "MAVKA timings of the mixed-trend first batch, one noise realisation; "
        "1800 s uses the widened windows"
    )
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=140)
    plt.close(fig)
    logger.info("wrote %s", OUT / name)


def main() -> None:
    """Write the data table and both figures."""
    records = collect()
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "figure_data.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    draw(
        records,
        "median_s",
        "median time minus conjunction (s, symlog)",
        "offset_vs_star.png",
        True,
    )
    draw(
        records,
        "mad_std_s",
        "scatter of times, 1.4826 x MAD (s)",
        "scatter_vs_star.png",
        False,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
