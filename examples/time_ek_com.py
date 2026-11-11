"""Time the eclipses of EK Com, sector 49, the way the regression test does.

The fit is the one in ``examples/fit_ek_com.py``. The timing recipe is the
one in ``tests/tess/test_regression_real.py``: the whole-phase window, both
sides of mid-eclipse required, and the contact-wall cut left off, because a
10 minute cadence does not sample both walls. The recorded counts are 70
primary and 68 secondary timings.

Run from the repository root::

    python examples/time_ek_com.py

Writes ``examples/out/ek_com_timing/``: the two timing files and ``oc.png``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import mikulasek_eb as eb
from mikulasek_eb.timing import fit_ephemeris
from mikulasek_eb.tom import write_timing_ascii

from fit_ek_com import DATA, run

logger = logging.getLogger(__name__)
OUT = Path(__file__).resolve().parents[1] / "examples" / "out" / "ek_com_timing"

# Same choices as the regression test. The default recipe also demands points
# on both contact walls, which this cadence does not provide.
RECIPE = eb.TimingRecipe(require_contact_walls=False, clip_outliers=False)


def draw_oc(timings: dict, path: Path, title: str) -> None:
    """Draw the O-C of the primary and the secondary about a linear ephemeris.

    Args:
        timings (dict): ``primary`` and ``secondary`` blocks from ``time``.
        path (Path): PNG file to write.
        title (str): Figure title.
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    for kind, colour in (("primary", "#0072B2"), ("secondary", "#D55E00")):
        block = timings[kind]
        if len(block["T"]) == 0:
            continue
        _coef, residual, _scale = fit_ephemeris(
            block["T"], block["sigma"], block["E"], 1
        )
        ax.errorbar(
            block["E"],
            residual * 86400.0,
            yerr=np.asarray(block["sigma"]) * 86400.0,
            fmt="o",
            ms=3,
            color=colour,
            label=kind,
        )
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.set_xlabel("Eclipse number E")
    ax.set_ylabel("O-C (seconds)")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def main() -> None:
    """Fit EK Com, time both eclipses, and write the files and the O-C figure."""
    payload, out = run()
    timed = eb.time(out.t_mjd, out.mag, out.sig, out.solution, RECIPE)
    OUT.mkdir(parents=True, exist_ok=True)
    for kind in ("primary", "secondary"):
        block = timed.timings[kind]
        path = OUT / f"{kind}.dat"
        write_timing_ascii(
            path,
            block,
            lightcurve=DATA,
            period=payload["period"],
            t0_jd=payload["epoch_jd"],
            method="mikulasek_eb.time",
            window=RECIPE.window,
            events=kind,
            extra_comments=[
                "require_contact_walls: false",
                "clip_outliers: false",
            ],
            timescale=payload["timescale"],
            rejected=timed.rejected[kind],
            recipe=RECIPE,
        )
        sigma_s = (
            86400.0 * float(np.median(block["sigma"]))
            if len(block["T"])
            else float("nan")
        )
        logger.info(
            "%s: %d timings, %d rejected, median sigma %.1f s",
            kind,
            len(block["T"]),
            len(timed.rejected[kind]),
            sigma_s,
        )
        logger.info("wrote %s", path)
    figure = OUT / "oc.png"
    draw_oc(timed.timings, figure, f"{payload['title']}\nO-C about a linear ephemeris")
    logger.info("wrote %s", figure)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
