"""Mark the timing windows on the folded light curve of SV Cen sector 64.

The windows are those used for the timing: 1.6 D for both minima, with D1 and D2
from the chosen fit. The scan of this sector adopted 1.6 D for the primary and
0.4 D for the secondary, and the secondary was widened to 1.6 D by decision, so
both are drawn. Phase zero is the fitted M0 of the primary. Nothing is fitted.

Run from ``examples``::

    python sv_cen_window_mark.py

Writes ``results/sv_cen/s64_fit/timing_scan/windows_on_fold.png``.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import nsv807_s30_tiered as tiered

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "sv_cen"
RUNS = {
    "10": "s10_fit",
    "11": "s11_fit",
    "37": "s37_fit",
    "38": "s38_fit",
    "64": "s64_fit",
    "99": "s99_fit_p4",
    "100": "s100_fit_p4",
}


def locked_half() -> float:
    """Phase half-width locked from the sector 64 secondary at 1.6 D.

    Returns:
        float: Half-window in cycles.
    """
    solution = json.loads(
        next((RESULTS / "s64_fit" / "chosen").glob("*.json")).read_text()
    )
    return 1.6 * float(solution["params"]["D2"])


def main() -> None:
    """Draw every sector with the one locked window shaded on both minima."""
    tiered.FLUX_ERR = None
    half = locked_half()
    for sector, run in RUNS.items():
        path = next(
            p
            for p in (RESULTS / "cleaned").glob(f"*_s{sector}_*.vot")
            if "_witherr" not in p.name
        )
        payload = tiered.load_sap(path)
        solution = json.loads(
            next((RESULTS / run / "chosen").glob("*.json")).read_text()
        )
        period = float(solution["ephemeris"]["P0"])
        m0 = float(solution["ephemeris"]["M0"]) + 2400000.5
        phase = ((payload["t"] - m0) / period) % 1.0
        fig, ax = plt.subplots(figsize=(12, 6))
        for shift in (-1, 0, 1):
            ax.scatter(phase + shift, payload["mag"], s=2, color="0.6")
        for centre, colour in ((0.0, "tab:blue"), (0.5, "tab:red")):
            for edge in (-1, 0, 1):
                ax.axvspan(
                    centre + edge - half,
                    centre + edge + half,
                    alpha=0.15,
                    color=colour,
                    label=(
                        f"{'primary' if centre == 0 else 'secondary'} {half:.3f} cycles"
                        if edge == 0
                        else None
                    ),
                )
        ax.set_xlim(-0.5, 1.5)
        ax.invert_yaxis()
        ax.set_xlabel("Phase from the fitted primary minimum")
        ax.set_ylabel("SAP magnitude")
        ax.set_title(f"SV Cen sector {sector}, locked window {half:.3f} cycles")
        ax.legend()
        fig.tight_layout()
        out = RESULTS / run / "windows_on_fold.png"
        fig.savefig(out, dpi=120)
        plt.close(fig)
        logger.info("wrote %s", out)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
