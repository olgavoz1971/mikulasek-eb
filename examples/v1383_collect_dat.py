"""Collect the V1383 Tau window-timing ``.dat`` files into one folder with comments.

Each file starts with ``#`` lines (star, light curve, minimum, method, window,
cleaning, template, flags, columns). Rows without a time would be ``#`` lines, so all
other lines are data. Mikulasek window timing, first set (all sectors) and the second set
of the 600 s sectors (``timing_w36_42``). MAVKA files are added when they exist. Reads ``results/v1383_tau/<sector>_fit_bkg/timing/``.

Run from ``examples``::

    python v1383_collect_dat.py

Writes ``results/v1383_tau/timing_dat/`` and ``index.csv``.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

from nsv807_collect_dat import header, mik_window_text, read_dat, write

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "v1383_tau"
OUT = RESULTS / "timing_dat"
LIGHT_CURVES = {
    "s42": ("TIC388995726_s42_TESS-SPOC_600s", "sector 42, TESS-SPOC, 600 s"),
    "s43": ("TIC388995726_s43_TESS-SPOC_600s", "sector 43, TESS-SPOC, 600 s"),
    "s44": ("TIC388995726_s44_TESS-SPOC_600s", "sector 44, TESS-SPOC, 600 s"),
    "s70": ("TIC388995726_s70_TESS-SPOC_200s", "sector 70, TESS-SPOC, 200 s"),
    "s71": ("TIC388995726_s71_TESS-SPOC_200s", "sector 71, TESS-SPOC, 200 s"),
}
KINDS = ("primary", "secondary")


def head_lines(
    lc: str, kind: str, summary: dict, solution: str, window: str, rule: str
) -> list[str]:
    """Header lines of one file.

    Args:
        lc (str): Sector label.
        kind (str): ``primary`` or ``secondary``.
        summary (dict): Row of ``cleaned_bkg/summary.csv`` of this file.
        solution (str): Name of the chosen candidate.
        window (str): Window sentence.
        rule (str): How the multiple of D was chosen.

    Returns:
        list[str]: Plain header lines.
    """
    stem, label = LIGHT_CURVES[lc]
    return [
        "Star: V1383 Tau (TIC 388995726), contact binary with a strong O'Connell effect",
        f"Light curve: {stem}.vot ({label}), SAP flux",
        f"Minimum: {kind}",
        "Method: Mikulasek phenomenological template, window timing (shift of the saved template, linear baseline)",
        f"Window: {window}; {rule}",
        f"Cleaning: points with sap_bkg above {summary['ratio']} x the median background removed: "
        f"{summary['removed']} of {summary['rows_in']} rows",
        f"Template: tiered fit {solution}, baseline split at gaps over 0.3 d",
        "Flag: ok, gap (gap in window), walls (contact walls not both covered)",
        "Time: full JD of the time system of the source file",
        "Columns: jd  err_days  flag",
        "err_days is the formal 1 sigma error of the time in days",
        "Lines starting with # are comments.",
    ]


MAVKA = Path("/home/voz/projects/UPJS/Mavka_user_folder/output/v1383_tau_2")
MAVKA_LABELS = {
    "s42": "s42_600",
    "s43": "s43_600",
    "s44": "s44_600",
    "s70": "s70_200",
    "s71": "s71_200",
}


def collect_mavka() -> list[tuple]:
    """Copy the MAVKA ``.dat`` files, which already carry ``#`` headers.

    AP, WSAP and WSAP2 come from ``v1383_tau_2``. WSL comes from the run that keeps
    times with an infinite sigma (``v1383_tau_2/WSL_sigma_none``); the earlier WSL
    files of ``v1383_tau_2`` hold no time and are not copied. A first line names
    the source file. The rest of each file is copied unchanged.

    Returns:
        list[tuple]: Index rows (file, sector, method, kind, number of times).
    """
    rows = []
    for lc, label in MAVKA_LABELS.items():
        for method, folder in (
            ("AP", MAVKA),
            ("WSAP", MAVKA),
            ("WSAP2", MAVKA),
            ("WSL", MAVKA / "WSL_sigma_none"),
        ):
            for kind in KINDS:
                src = folder / f"{label}_{method}_{kind}.dat"
                if not src.exists():
                    continue
                text = src.read_text().splitlines()
                n = sum(1 for line in text if line.strip() and not line.startswith("#"))
                name = f"V1383Tau_{lc}_mavka_{method}_{kind}.dat"
                (OUT / name).write_text(
                    f"# Collected from {src}\n" + "\n".join(text) + "\n"
                )
                rows.append((name, lc, f"mavka_{method}", kind, n))
    return rows


def main() -> None:
    """Write all collected files and the index."""
    OUT.mkdir(parents=True, exist_ok=True)
    with open(RESULTS / "cleaned_bkg" / "summary.csv", newline="") as handle:
        summaries = {r["file"]: r for r in csv.DictReader(handle)}
    index = []
    for lc, (stem, _) in LIGHT_CURVES.items():
        run = RESULTS / f"{lc}_fit_bkg"
        solution = next((run / "chosen").glob("*.json")).stem.split("__")[1]
        sets = [
            (
                "timing",
                "mikulasek_window",
                "3.3 D (primary) and 3.7 D (secondary) of the D of this sector, from the scan on sector 70",
            )
        ]
        if (run / "timing_w36_42").exists():
            sets.append(
                (
                    "timing_w36_42",
                    "mikulasek_window_wide",
                    "3.6 D (primary) and 4.2 D (secondary) of the D of this sector, from the scan on "
                    "sector 44 started at the windows of sector 70 (second set for the 600 s sectors)",
                )
            )
        for kind in KINDS:
            for sub, method, rule in sets:
                folder = run / sub
                head = head_lines(
                    lc,
                    kind,
                    summaries[f"{stem}.vot"],
                    solution,
                    mik_window_text(folder, kind),
                    rule,
                )
                name = f"V1383Tau_{lc}_{method}_{kind}.dat"
                n = write(
                    OUT / name, head, read_dat(folder / f"{stem}_{kind}.dat"), None
                )
                index.append((name, lc, method, kind, n))
    index.extend(collect_mavka())
    with open(OUT / "index.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "light_curve", "method", "kind", "n_times"])
        writer.writerows(index)
    logger.info("wrote %d files to %s", len(index), OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
