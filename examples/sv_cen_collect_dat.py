"""Collect the SV Cen window timings into one folder with commented headers.

One ``.dat`` file per sector and per minimum. Every line that is not a time starts
with ``#`` (star, light curve, minimum, method, window, the template, the cleaning,
the sigma and the columns), and a rejected eclipse is a ``#`` line with its reason,
so every other line is data: ``jd  err_days  flag``.

The window is 1.6 D for both minima. The scan of sector 64 adopted 1.6 D for the
primary and 0.4 D for the secondary, and the secondary was widened to the primary's
window by decision.

Run from ``examples``::

    python sv_cen_collect_dat.py

Writes ``results/sv_cen/timing_dat/`` and ``index.csv``.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "sv_cen"
OUT = RESULTS / "timing_dat"
RUNS = {
    "10": ("s10_fit", "1800 s", "70"),
    "11": ("s11_fit", "1800 s", "70"),
    "37": ("s37_fit", "600 s", "200"),
    "38": ("s38_fit", "600 s", "200"),
    "64": ("s64_fit", "200 s", "400"),
    "99": ("s99_fit_p4", "600 s", "200"),
    "100": ("s100_fit_p4", "200 s", "400"),
}
KINDS = ("primary", "secondary")
WINDOW = (
    "0.198 cycles for both minima of every sector, the phase width of the sector 64 "
    "secondary at 1.6 D, locked in phase and not as a multiple of D"
)


def main() -> None:
    """Write the collected files and the index."""
    OUT.mkdir(parents=True, exist_ok=True)
    with open(RESULTS / "cleaned" / "summary.csv", newline="") as handle:
        summary = {r["file"]: r for r in csv.DictReader(handle)}
    index = []
    for sector, (run, cadence, bins) in RUNS.items():
        folder = RESULTS / run
        chosen = list((folder / "chosen").glob("*.json"))
        if len(chosen) != 1:
            raise FileNotFoundError(f"{folder}/chosen holds {len(chosen)} solutions")
        solution = chosen[0].stem.split("__")[1]
        rows = list(csv.DictReader(open(folder / "timing" / "eclipses.csv")))
        stem = chosen[0].stem.split("__")[0]
        cleaning = summary[f"{stem}.vot"]
        for kind in KINDS:
            of_kind = [r for r in rows if r["kind"] == kind]
            timed = [r for r in of_kind if r["T_jd"]]
            half = timed[0]["half_cycles"] if timed else "none"
            head = [
                "# Star: SV Cen (TIC 324101484), eclipsing binary with a fast period change",
                f"# Light curve: {stem}.vot (sector {sector}, QLP, {cadence}), SAP flux",
                f"# Minimum: {kind}",
                "# Method: Mikulasek phenomenological template, window timing "
                "(shift of the saved template, linear baseline)",
                f"# Window: {WINDOW}",
                f"# Half-window of this file: {half} cycles",
                f"# Template: {solution}, from the tiered fit of this sector "
                f"(start profile {solution[1]}, Pdot fixed at 0)",
                f"# Cleaning: folded curve, {bins} phase bins, 5 sigma, local MAD, "
                f"removed once; {cleaning['removed']} of {cleaning['rows_in']} rows "
                f"removed, {cleaning['flagged']} flagged",
                "# Sigma: one constant per sector, the robust scatter about the "
                "phase-bin medians; the files have no SAP error column",
                "# Time: TDB, barycentric, JD",
                "# Columns: jd  err_days  flag",
                "# A line without a time is a comment and gives the reason",
            ]
            name = f"SVCen_s{sector}_mikulasek_window_{kind}.dat"
            with (OUT / name).open("w") as handle:
                handle.write("\n".join(head) + "\n")
                for r in of_kind:
                    if r["T_jd"]:
                        flags = []
                        if r["walls_ok"] != "True":
                            flags.append("walls")
                        if r["gap_in_window"] == "True":
                            flags.append("gap")
                        handle.write(
                            f"{float(r['T_jd']):.8f}  "
                            f"{float(r['sigma_s']) / 86400.0:.8f}  "
                            f"{'+'.join(flags) or 'ok'}\n"
                        )
                    else:
                        handle.write(f"# E={r['cycle']} rejected: {r['reason']}\n")
            index.append((name, sector, kind, solution, len(timed)))
            logger.info("%s: %d times", name, len(timed))
    with (OUT / "index.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "sector", "kind", "template", "n_times"])
        writer.writerows(index)
    logger.info("wrote %d files to %s", len(index), OUT)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
