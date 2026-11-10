"""Collect every NSV 807 timing ``.dat`` file of every method into one folder.

Each output file starts with ``#`` comment lines (star, light curve, method, window,
cleaning, columns). Every row that holds no time is written as a ``#`` line too, so
all remaining lines are data. Sources:

* Mikulasek window timing: ``results/nsv807/<run>/timing/``.
* Mikulasek whole-cycle timing: ``results/nsv807/<run>/whole_cycle/``.
* Mikulasek borrowed template (sector 4): ``s4_protocol/borrow/<variant>/timing/``.
* MAVKA AP, WSAP, WSL: ``Mavka_user_folder/output/nsv807_cleaned/``.

Run from ``examples``::

    python nsv807_collect_dat.py

Writes ``results/nsv807/timing_dat/``.
"""

from __future__ import annotations

import csv
import logging
import re
from pathlib import Path

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "nsv807"
MAVKA = Path("/home/voz/projects/UPJS/Mavka_user_folder/output/nsv807_cleaned")
OUT = RESULTS / "timing_dat"
LIGHT_CURVES = {
    "s30_120": ("s30_protocol", "TIC35119266_s30_SPOC_120s", "sector 30, SPOC, 120 s"),
    "s30_600": (
        "s30_600s_protocol",
        "TIC35119266_s30_TESS-SPOC_600s",
        "sector 30, TESS-SPOC, 600 s",
    ),
    "s4_1800": (
        "s4_protocol",
        "TIC35119266_s4_TESS-SPOC_1800s",
        "sector 4, TESS-SPOC, 1800 s",
    ),
    "s97_200": ("s97_protocol", "TIC35119266_s97_QLP_200s", "sector 97, QLP, 200 s"),
}
BORROW_VARIANTS = (
    "borrow_phase",
    "borrow_amp",
    "borrow_width",
    "borrow_amp_width",
    "borrow_amp_width_harm",
)
KINDS = ("primary", "secondary")


def cleaning_text(rows: list[dict], lc: str, stem: str) -> str:
    """Describe the outlier cleaning of one light curve.

    Args:
        rows (list[dict]): Rows of ``cleaned/summary.csv``.
        lc (str): Light-curve label.
        stem (str): File stem.

    Returns:
        str: One sentence for the header.
    """
    for row in rows:
        if row["file"] == f"{stem}.vot":
            return (
                f"outliers removed before timing (400 phase bins, 5 sigma, MAD): "
                f"{row['removed']} of {row['rows_in']} rows removed"
            )
    return "no outlier removal (all points kept)"


def header(lines: list[str]) -> list[str]:
    """Prefix every line with ``# ``.

    Args:
        lines (list[str]): Plain header lines.

    Returns:
        list[str]: Comment lines.
    """
    return [f"# {line}".rstrip() for line in lines]


def common(lc: str, kind: str, method: str, window: str, extra: list[str]) -> list[str]:
    """Header lines shared by every file.

    Args:
        lc (str): Light-curve label.
        kind (str): ``primary`` or ``secondary``.
        method (str): Method description.
        window (str): Window description.
        extra (list[str]): Method-specific lines.

    Returns:
        list[str]: Plain header lines.
    """
    _, stem, label = LIGHT_CURVES[lc]
    with open(RESULTS / "cleaned" / "summary.csv", newline="") as handle:
        summary = list(csv.DictReader(handle))
    return [
        "Star: NSV 807 (TIC 35119266), W UMa type",
        f"Light curve: {stem}.vot ({label}), SAP flux",
        f"Minimum: {kind}",
        f"Method: {method}",
        f"Window: {window}",
        f"Cleaning: {cleaning_text(summary, lc, stem)}",
        *extra,
        "Time: full JD of the time system of the source file",
        "Columns: jd  err_days  flag",
        "err_days is the formal 1 sigma error of the time in days",
        "Lines starting with # are comments. Rows without a time are comments too.",
    ]


def read_dat(path: Path) -> list[list[str]]:
    """Read a plain ``.dat`` file.

    Args:
        path (Path): File with columns separated by two spaces.

    Returns:
        list[list[str]]: Split rows.
    """
    with open(path) as handle:
        return [line.split() for line in handle if line.strip()]


def mik_window_text(folder: Path, kind: str) -> str:
    """Describe the half-window of a Mikulasek window timing.

    Args:
        folder (Path): Folder with ``eclipses.csv``.
        kind (str): ``primary`` or ``secondary``.

    Returns:
        str: Window sentence.
    """
    with open(folder / "eclipses.csv", newline="") as handle:
        half = sorted(
            {
                float(r["half_cycles"])
                for r in csv.DictReader(handle)
                if r["kind"] == kind
            }
        )
    return (
        f"centre +/- {half[0]:.5f} cycles"
        if len(half) == 1
        else f"centre +/- {half[0]:.5f} to {half[-1]:.5f} cycles"
    )


def write(
    path: Path, head: list[str], rows: list[list[str]], none_marker: str | None
) -> int:
    """Write one collected file.

    Args:
        path (Path): Destination.
        head (list[str]): Plain header lines.
        rows (list[list[str]]): Source rows, ``jd err flag``.
        none_marker (str | None): First value of a row without a time (``none``) or
            None when the source has no such rows.

    Returns:
        int: Number of data rows written.
    """
    data = [r for r in rows if none_marker is None or r[0] != none_marker]
    skipped = [r for r in rows if none_marker is not None and r[0] == none_marker]
    lines = header(head)
    lines.append(
        f"# rows with a time: {len(data)}, rows without a time: {len(skipped)}"
    )
    lines += [f"# no time  {r[-1]}" for r in skipped]
    lines += ["  ".join(r) for r in data]
    path.write_text("\n".join(lines) + "\n")
    return len(data)


def mavka_halves() -> dict:
    """Read the chosen MAVKA half-widths from ``scan/chosen.txt``.

    Returns:
        dict: ``(method, kind)`` to half-width in cycles and a flag for borrowing.
    """
    text = (MAVKA / "scan" / "chosen.txt").read_text()
    out = {}
    for m in re.finditer(r"^(AP|WSAP) (primary|secondary) ([0-9.]+) rms", text, re.M):
        out[(m[1], m[2])] = (float(m[3]), False)
    for m in re.finditer(
        r"^WSL (primary|secondary) scan , borrowed from WSAP ([0-9.]+)", text, re.M
    ):
        out[("WSL", m[1])] = (float(m[2]), True)
    return out


def main() -> None:
    """Write all collected files and an index."""
    OUT.mkdir(parents=True, exist_ok=True)
    halves = mavka_halves()
    index = []
    for lc, (run, stem, _) in LIGHT_CURVES.items():
        solution = next((RESULTS / run / "chosen").glob("*.json")).stem.split("__")[1]
        for kind in KINDS:
            src = RESULTS / run / "timing" / f"{stem}_{kind}.dat"
            head = common(
                lc,
                kind,
                "Mikulasek phenomenological template, window timing (shift of the saved template, linear baseline)",
                mik_window_text(RESULTS / run / "timing", kind),
                [
                    f"Template: tiered fit {solution}",
                    "Flag: ok, gap (gap in window), walls (contact walls not both covered)",
                ],
            )
            name = f"NSV807_{lc}_mikulasek_window_{kind}.dat"
            n = write(OUT / name, head, read_dat(src), None)
            index.append((name, lc, "mikulasek_window", kind, n))
            src = RESULTS / run / "whole_cycle" / f"{stem}_{kind}.dat"
            head = common(
                lc,
                kind,
                "Mikulasek phenomenological template, whole-cycle timing (one template shift per cycle, linear baseline per cycle)",
                "phase (n - 0.25, n + 0.75) of cycle n; primary and secondary of a cycle share one time",
                [
                    f"Template: tiered fit {solution}",
                    "Flag: ok, gap (gap longer than 2 cadences inside the cycle)",
                ],
            )
            name = f"NSV807_{lc}_mikulasek_whole_cycle_{kind}.dat"
            n = write(OUT / name, head, read_dat(src), None)
            index.append((name, lc, "mikulasek_whole_cycle", kind, n))
            for method in ("AP", "WSAP", "WSL"):
                half, borrowed = halves[(method, kind)]
                head = common(
                    lc,
                    kind,
                    f"MAVKA {method}, window timing of each minimum",
                    f"centre +/- {half:.6f} cycles (chosen by the scan on s30_120"
                    + (
                        "; borrowed from WSAP, WSL was not scanned)"
                        if borrowed
                        else ")"
                    )
                    + "; windows with fewer than 6 points were widened and are flagged extended",
                    [
                        "Flag: ok, or extended (window widened to reach 6 points)",
                        "Time system: TDB, barycentric (as written by MAVKA)",
                    ],
                )
                src = MAVKA / f"{lc}_{method}_{kind}.dat"
                name = f"NSV807_{lc}_mavka_{method}_{kind}.dat"
                n = write(OUT / name, head, read_dat(src), "none")
                index.append((name, lc, f"mavka_{method}", kind, n))
    run = "s4_protocol"
    stem = LIGHT_CURVES["s4_1800"][1]
    for variant in BORROW_VARIANTS:
        folder = RESULTS / run / "borrow" / variant / "timing"
        for kind in KINDS:
            head = common(
                "s4_1800",
                kind,
                f"Mikulasek template borrowed from sector 30 (120 s), variant {variant}, window timing",
                mik_window_text(folder, kind),
                [
                    "Template: p3_h123_b0 of sector 30 (120 s), borrowed; epoch fitted; variant name lists what else was free",
                    "Flag: ok, gap (gap in window), walls (contact walls not both covered)",
                ],
            )
            name = f"NSV807_s4_1800_mikulasek_{variant}_{kind}.dat"
            n = write(OUT / name, head, read_dat(folder / f"{stem}_{kind}.dat"), None)
            index.append((name, "s4_1800", variant, kind, n))
    with open(OUT / "index.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "light_curve", "method", "kind", "n_times"])
        writer.writerows(index)
    logger.info("wrote %d files to %s", len(index), OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
