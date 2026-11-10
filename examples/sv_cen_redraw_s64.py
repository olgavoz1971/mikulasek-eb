"""Redraw the sector 64 timing figures with the template at the right phase.

``nsv807_timing.draw_each`` evaluates the template at the phase measured from the
fitted time. The eclipse numbers of these sectors are counted from the file epoch,
so a primary carries an odd cycle number and that phase lands on the secondary. The
template is evaluated here at the phase of the fitted ephemeris (``mikulasek_eb.phase``
with the saved M0 and P0), which puts the primary at phase 0. One constant offset,
the mean residual in the window, is added. The template is not changed and nothing
is retimed. The existing figures are removed first.

Run from ``examples``::

    python sv_cen_redraw_s64.py

Rewrites ``results/sv_cen/s64_fit/timing/fits/<kind>/E_<cycle>.png``.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import nsv807_timing as timing

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "sv_cen"
OFFSET = {"primary": 0.0, "secondary": 0.5}
SECTORS = (
    ("10", "s10_fit", 70),
    ("11", "s11_fit", 70),
    ("37", "s37_fit", 200),
    ("38", "s38_fit", 200),
    ("64", "s64_fit", 400),
    ("99", "s99_fit_p4", 200),
    ("100", "s100_fit_p4", 400),
)


def main() -> None:
    """Replace the per-eclipse figures of every sector."""
    import sv_cen_timing as sv

    for sector, run_name, bins in SECTORS:
        run = RESULTS / run_name
        path = next(
            p
            for p in (RESULTS / "cleaned").glob(f"*_s{sector}_*.vot")
            if "_witherr" not in p.name
        )
        one(sv, run, path, bins)


def one(sv, run: Path, path: Path, bins: int) -> None:
    """Redraw the figures of one sector.

    Args:
        sv (module): ``sv_cen_timing``, for the corrected loader.
        run (Path): Fit folder of the sector.
        path (Path): Cleaned VOTable.
        bins (int): Phase bins used for the sigma of this sector.
    """
    payload, t_mjd, solution, _ = sv.prepare(run, str(path), bins)
    t_jd = np.asarray(timing.eb.mjd_to_jd(t_mjd), float)
    rows = list(__import__("csv").DictReader(open(run / "timing" / "eclipses.csv")))
    model, params = solution.build_model(), solution.params
    p0 = float(solution.ephemeris["P0"])
    fits = run / "timing" / "fits"
    if fits.exists():
        shutil.rmtree(fits)
    stem = path.stem
    drawn = 0
    for r in rows:
        if r["T_jd"] == "":
            continue
        kind, centre, half = r["kind"], float(r["T_jd"]), float(r["half_cycles"])
        sel = np.abs(t_jd - centre) <= half * p0
        theta = timing.eb.phase(t_mjd, solution.ephemeris["M0"], p0, 0.0)
        template = np.asarray(model(theta[sel], params), float)
        offset = float(np.mean(payload["mag"][sel] - template))
        at_centre = float(
            timing.eb.phase(
                np.array([centre - 2400000.5]), solution.ephemeris["M0"], p0, 0.0
            )[0]
        )
        grid = np.linspace(-half * p0, half * p0, 400)
        curve = (
            np.asarray(
                model(at_centre + grid / p0, params),
                float,
            )
            + offset
        )
        fig, axes = plt.subplots(
            2, 1, figsize=(8, 7), gridspec_kw={"height_ratios": [3, 1]}, sharex=True
        )
        axes[0].plot(
            t_jd[sel] - centre,
            payload["mag"][sel],
            ".",
            ms=4,
            color="0.4",
            label="data",
        )
        axes[0].plot(
            grid, curve, color="red", lw=1.5, label="template + constant offset"
        )
        axes[0].axvline(0.0, color="blue", ls=":", lw=1, label="fitted time")
        axes[0].invert_yaxis()
        axes[0].set_ylabel("Magnitude")
        axes[0].legend()
        axes[0].set_title(
            f"{stem}\n{kind} E={r['cycle']}  T={centre:.6f}  "
            f"sigma {r['sigma_s']} s  n={int(sel.sum())}  {timing.flag_text(r)}"
        )
        axes[1].plot(
            t_jd[sel] - centre,
            1e3 * (payload["mag"][sel] - template - offset),
            ".",
            ms=4,
            color="0.4",
        )
        axes[1].axhline(0.0, color="red", lw=1)
        axes[1].set_xlabel("Time from the fitted time (d)")
        axes[1].set_ylabel("Residual (mmag)")
        fig.tight_layout()
        folder = fits / kind
        folder.mkdir(parents=True, exist_ok=True)
        fig.savefig(folder / f"E_{r['cycle']}.png", dpi=100)
        plt.close(fig)
        drawn += 1
    logger.info("redrew %d figures in %s", drawn, fits)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
