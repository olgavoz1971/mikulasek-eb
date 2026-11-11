"""Collect the NSV 807 maximum timings into one folder.

The header follows the minimum collection in ``nsv807_collect_dat.py``: star,
light curve, method, window, cleaning, template, flag and columns. A maximum
without a time is a comment, and every other line is ``jd  err_days  flag``.

Run from ``mikulasek_eb``::

    python examples/nsv807_collect_maxima.py

Writes ``results/nsv807/timing_max_dat/`` and ``index.csv``.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path

import nsv807_collect_dat as minima
from window_scan_mik import DAY_S

logger = logging.getLogger(__name__)

RESULTS = minima.RESULTS
OUT = RESULTS / "timing_max_dat"
KINDS = ("primary_maximum", "secondary_maximum")
LABEL = {"primary_maximum": "primary", "secondary_maximum": "secondary"}


def rows_of(folder: Path, kind: str) -> list[dict]:
    """Read one maximum kind from a sector table.

    Args:
        folder (Path): ``maxima`` folder of one sector.
        kind (str): ``primary_maximum`` or ``secondary_maximum``.

    Returns:
        list[dict]: Rows of that kind, in table order.

    Raises:
        RuntimeError: If the template phase or the half-window is not one value.
    """
    with (folder / "eclipses.csv").open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["kind"] == kind]
    if not rows:
        raise RuntimeError(f"{folder} has no {kind}")
    phases = {float(row["template_phase"]) for row in rows}
    halves = {float(row["half_cycles"]) for row in rows}
    if len(phases) != 1 or len(halves) != 1:
        raise RuntimeError(f"{folder} {kind}: phase {phases}, half-window {halves}")
    return rows


def write_one(path: Path, head: list[str], rows: list[dict]) -> int:
    """Write one collected file.

    Args:
        path (Path): Destination.
        head (list[str]): Plain header lines.
        rows (list[dict]): Events of one kind.

    Returns:
        int: Number of rows that hold a time.
    """
    timed = [row for row in rows if row["T_jd"]]
    missed = [row for row in rows if not row["T_jd"]]
    lines = minima.header(head)
    lines.append(
        f"# rows with a time: {len(timed)}, rows without a time: {len(missed)}"
    )
    lines += [f"# no time  E={row['cycle']} {row['reason']}" for row in missed]
    for row in timed:
        err = float(row["sigma_s"]) / DAY_S
        flag = "gap" if row["gap_in_window"] == "True" else "ok"
        lines.append(f"{float(row['T_jd']):.8f}  {err:.8f}  {flag}")
    path.write_text("\n".join(lines) + "\n")
    return len(timed)


def sigma_line(saved: dict) -> str:
    """Describe the magnitude error used for one sector.

    Args:
        saved (dict): ``result.json`` of the maxima run.

    Returns:
        str: One header line.
    """
    if saved["constant_sigma_mmag"] is None:
        return "Sigma: sap_flux_err of the light curve"
    return (
        "Sigma: one constant, "
        f"{saved['constant_sigma_mmag']:.2f} mmag, from the chi-square of the saved solution"
    )


def main() -> None:
    """Write the eight maximum files and the index."""
    OUT.mkdir(parents=True, exist_ok=True)
    index = []
    for lc, (run, stem, _) in minima.LIGHT_CURVES.items():
        folder = RESULTS / run / "maxima"
        saved = json.loads((folder / "result.json").read_text())
        solution = next((RESULTS / run / "chosen").glob("*.json")).stem.split("__")[1]
        with open(RESULTS / "cleaned" / "summary.csv", newline="") as handle:
            summary = list(csv.DictReader(handle))
        for kind in KINDS:
            rows = rows_of(folder, kind)
            phase = float(rows[0]["template_phase"])
            half = float(rows[0]["half_cycles"])
            head = [
                "Star: NSV 807 (TIC 35119266), W UMa type",
                f"Light curve: {stem}.vot ({minima.LIGHT_CURVES[lc][2]}), SAP flux",
                f"Maximum: {LABEL[kind]}",
                "Method: Mikulasek phenomenological template, window timing of the maximum "
                "(shift of the saved template, linear baseline)",
                f"Window: centre at template phase {phase:.6f}, the brightest point of the "
                f"template, half-window {half:.2f} cycles",
                f"Cleaning: {minima.cleaning_text(summary, lc, stem)}",
                f"Template: tiered fit {solution}",
                sigma_line(saved),
                "Flag: ok, gap (gap longer than 2 cadences in the window). "
                "The eclipse contact-wall flag is not used",
                "Time: full JD of the time system of the source file",
                "Columns: jd  err_days  flag",
                "err_days is the formal 1 sigma error of the time in days",
                "Lines starting with # are comments. Rows without a time are comments too.",
            ]
            name = f"NSV807_{lc}_mikulasek_maximum_{LABEL[kind]}.dat"
            n_times = write_one(OUT / name, head, rows)
            index.append((name, lc, "mikulasek_maximum", LABEL[kind], n_times))
            logger.info("%s: %d times", name, n_times)
    with open(OUT / "index.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "light_curve", "method", "kind", "n_times"])
        writer.writerows(index)
    logger.info("wrote %d files to %s", len(index), OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
