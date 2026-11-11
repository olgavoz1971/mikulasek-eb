"""Collect the V1383 Tau maximum timings into one folder.

The header follows the maximum collection of NSV 807. A maximum without a time
is a comment, and every other line is ``jd  err_days  flag``. Each secondary
is written twice: once as measured, and once with every time reduced by half
the file period. The file period is the one stored in the chosen solutions,
and it is the same in every sector.

Run from ``mikulasek_eb``::

    python examples/v1383_collect_maxima.py

Writes ``results/v1383_tau/timing_max_dat/`` and ``index.csv``.
"""

from __future__ import annotations

import csv
import json
import logging
from decimal import Decimal
from pathlib import Path

import nsv807_collect_maxima as collected
import v1383_collect_dat as minima
from nsv807_collect_dat import header
from window_scan_mik import DAY_S

logger = logging.getLogger(__name__)

RESULTS = minima.RESULTS
OUT = RESULTS / "timing_max_dat"
RUNS = tuple(f"{lc}_fit_bkg" for lc in minima.LIGHT_CURVES)


def file_half_period() -> Decimal:
    """Half the file period, required to be the same in every sector.

    Returns:
        Decimal: Half period in days.

    Raises:
        RuntimeError: If the saved file periods are not all equal.
    """
    found = []
    for run in RUNS:
        path = next((RESULTS / run / "chosen").glob("*.json"))
        saved = json.loads(path.read_text())
        found.append(f"{float(saved['provenance']['period_input_days']):.5f}")
    if any(period != found[0] for period in found[1:]):
        raise RuntimeError(f"the file periods differ: {found}")
    half = Decimal(found[0]) / 2
    if half * 2 != Decimal(found[0]):
        raise RuntimeError(f"half of {found[0]} is not exact")
    return half


def write_shifted(path: Path, head: list[str], rows: list[dict], half: Decimal) -> int:
    """Write one secondary file with the times reduced by ``half``.

    Args:
        path (Path): Destination.
        head (list[str]): Plain header lines of the unshifted file.
        rows (list[dict]): Secondary events.
        half (Decimal): Half the file period, days.

    Returns:
        int: Number of rows that hold a time.
    """
    timed = [row for row in rows if row["T_jd"]]
    missed = [row for row in rows if not row["T_jd"]]
    lines = [
        f"# Times in this copy are the source times minus {half} d, half the file period."
    ]
    lines += header(head)
    lines.append(
        f"# rows with a time: {len(timed)}, rows without a time: {len(missed)}"
    )
    lines += [f"# no time  E={row['cycle']} {row['reason']}" for row in missed]
    for row in timed:
        shown = Decimal(f"{float(row['T_jd']):.8f}")
        err = Decimal(f"{float(row['sigma_s']) / DAY_S:.8f}")
        flag = "gap" if row["gap_in_window"] == "True" else "ok"
        lines.append(f"{shown - half:.8f}  {err:.8f}  {flag}")
    path.write_text("\n".join(lines) + "\n")
    return len(timed)


def main() -> None:
    """Write the maximum files, the shifted secondary copies and the index."""
    half = file_half_period()
    OUT.mkdir(parents=True, exist_ok=True)
    with open(RESULTS / "cleaned_bkg" / "summary.csv", newline="") as handle:
        summaries = {row["file"]: row for row in csv.DictReader(handle)}
    index = []
    for lc, (stem, label) in minima.LIGHT_CURVES.items():
        run = RESULTS / f"{lc}_fit_bkg"
        folder = run / "maxima"
        saved = json.loads((folder / "result.json").read_text())
        if saved["sigma"] != "sap_flux_err":
            raise RuntimeError(f"{lc}: sigma {saved['sigma']}")
        solution = next((run / "chosen").glob("*.json")).stem.split("__")[1]
        summary = summaries[f"{stem}.vot"]
        for kind in collected.KINDS:
            rows = collected.rows_of(folder, kind)
            phase = float(rows[0]["template_phase"])
            window = float(rows[0]["half_cycles"])
            head = [
                "Star: V1383 Tau (TIC 388995726), contact binary with a strong O'Connell effect",
                f"Light curve: {stem}.vot ({label}), SAP flux",
                f"Maximum: {collected.LABEL[kind]}",
                "Method: Mikulasek phenomenological template, window timing of the maximum "
                "(shift of the saved template, linear baseline)",
                f"Window: centre at template phase {phase:.6f}, the brightest point of the "
                f"template, half-window {window:.2f} cycles",
                "Cleaning: points with sap_bkg above "
                f"{summary['ratio']} x the median background removed: "
                f"{summary['removed']} of {summary['rows_in']} rows",
                f"Template: tiered fit {solution}, baseline split at gaps over 0.3 d",
                "Sigma: sap_flux_err of the light curve",
                "Flag: ok, gap (gap longer than 2 cadences in the window). "
                "The eclipse contact-wall flag is not used",
                "Time: full JD of the time system of the source file",
                "Columns: jd  err_days  flag",
                "err_days is the formal 1 sigma error of the time in days",
                "Lines starting with # are comments. Rows without a time are comments too.",
            ]
            name = f"V1383Tau_{lc}_mikulasek_maximum_{collected.LABEL[kind]}.dat"
            n_times = collected.write_one(OUT / name, head, rows)
            index.append(
                (name, lc, "mikulasek_maximum", collected.LABEL[kind], n_times)
            )
            logger.info("%s: %d times", name, n_times)
            if kind != "secondary_maximum":
                continue
            shifted = f"V1383Tau_{lc}_mikulasek_maximum_secondary_minus_half_period.dat"
            n_shift = write_shifted(OUT / shifted, head, rows, half)
            if n_shift != n_times:
                raise RuntimeError(f"{shifted}: {n_shift} times, source has {n_times}")
            index.append(
                (
                    shifted,
                    lc,
                    "mikulasek_maximum",
                    "secondary_minus_half_period",
                    n_shift,
                )
            )
            logger.info("%s: %d times", shifted, n_shift)
    with open(OUT / "index.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "light_curve", "method", "kind", "n_times"])
        writer.writerows(index)
    logger.info("wrote %d files to %s", len(index), OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
