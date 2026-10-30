"""Supervised check of an automatic selection: numbers and plots.

Runs ``select_model`` on a star, then draws the fold with the fit, residuals
against phase (binned) and against time, and zooms on both eclipses. It prints
the residual rms in eclipse, out of eclipse and after binning in phase.

Run::

    python examples/diagnose_fit.py "AR Aur S19"
    python examples/diagnose_fit.py path/to/lc.vot
    python examples/diagnose_fit.py "NSV 807 S04" --show   # interactive window (zoom, pan)
"""

from __future__ import annotations

import logging
import sys

import matplotlib

if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from select_models import OUT, resolve

from mikulasek_eb import select_model
from mikulasek_eb.helpers import absolute_jd_to_mjd, predicted_magnitude
from mikulasek_eb.io import load_tess_votable
from mikulasek_eb.model import phase
from mikulasek_eb.tuning import residuals

logger = logging.getLogger(__name__)


def binned(
    ph: np.ndarray, y: np.ndarray, n: int = 100
) -> tuple[np.ndarray, np.ndarray]:
    """Median of ``y`` in ``n`` equal phase bins.

    Args:
        ph (ndarray): Phase in [0, 1).
        y (ndarray): Values.
        n (int): Number of bins.

    Returns:
        tuple[ndarray, ndarray]: Bin centres and medians (NaN for empty bins).
    """
    idx = np.minimum((ph * n).astype(int), n - 1)
    med = np.array(
        [np.median(y[idx == k]) if np.any(idx == k) else np.nan for k in range(n)]
    )
    return (np.arange(n) + 0.5) / n, med


def main() -> None:
    """Select a model, print residual statistics and save the diagnostic PNG."""
    logging.basicConfig(level=logging.WARNING)
    args = [a for a in sys.argv[1:] if a != "--show"]
    label, path = resolve(args[0])
    payload = load_tess_votable(path)
    t = absolute_jd_to_mjd(payload["t"])
    r = select_model(
        t,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        data_id=label,
    )
    pl = r.output.pipeline
    ad = pl.adopted
    yhat, tmpl, theta = predicted_magnitude(
        pl.t_mjd, pl.mag, pl.sig, pl.model, ad.p, ad.ephem, pl.seg, pl.order
    )
    res, th = residuals(pl.t_mjd, pl.mag, pl.sig, pl.model, ad, pl.seg, pl.order)
    ph = np.mod(th, 1.0)
    d = np.minimum(ph, 1 - ph)
    d2 = np.abs(ph - 0.5)
    ecl = (d < 0.08) | (d2 < 0.08)
    bx, bres = binned(ph, res)
    sig = pl.sig
    print(f"{label}: best {r.best.name}, n={len(res)}, ask_human={r.ask_human}")
    print(
        f"  rms all {1e3*np.std(res):.2f} mmag, median catalogue sigma {1e3*np.median(sig):.2f}"
    )
    print(
        f"  rms in eclipse {1e3*np.std(res[ecl]):.2f}, out {1e3*np.std(res[~ecl]):.2f} mmag"
    )
    print(
        f"  rms of phase-binned residuals {1e3*np.nanstd(bres):.2f} mmag (structure not noise)"
    )
    print(
        f"  worst bin {1e3*np.nanmax(np.abs(bres)):.2f} mmag at phase {bx[np.nanargmax(np.abs(bres))]:.2f}"
    )
    print(f"  reasons: {r.reasons}")

    fig, ax = plt.subplots(3, 2, figsize=(13, 11))
    o = np.argsort(ph)
    grid = np.linspace(0.0, 1.0, 2001)
    curve = pl.model(grid, ad.p)  # pure template, no baseline
    flat = pl.mag - (yhat - tmpl)  # data with the fitted baseline removed
    ax[0, 0].plot(ph, flat, ".", ms=1.5, color="0.6")
    ax[0, 0].plot(grid, curve, "r-", lw=1)
    ax[0, 0].invert_yaxis()
    ax[0, 0].set_title(f"{label}: {r.best.name} fold")
    ax[0, 1].plot(ph, 1e3 * res, ".", ms=1.5, color="0.6")
    ax[0, 1].plot(bx, 1e3 * bres, "r-o", ms=3)
    ax[0, 1].set_title("residuals (mmag) vs phase, binned in red")
    ax[1, 0].plot(pl.t_mjd, 1e3 * res, ".", ms=1.5)
    ax[1, 0].set_title("residuals (mmag) vs time (MJD)")
    for k, c in enumerate((0.0, 0.5)):
        x = (ph - c + 0.5) % 1 - 0.5
        a = ax[1 + k, 1] if k == 0 else ax[2, 0]
        s = np.abs(x) < 0.15
        a.plot(x[s], flat[s], ".", ms=2, color="0.6")
        g = (grid - c + 0.5) % 1 - 0.5
        gs = np.abs(g) < 0.15
        og = np.argsort(g[gs])
        a.plot(g[gs][og], curve[gs][og], "r-")
        a.invert_yaxis()
        a.set_title("primary zoom" if k == 0 else "secondary zoom")
    ax[2, 1].hist(res / sig, bins=60)
    ax[2, 1].set_title("residual / catalogue sigma")
    fig.tight_layout()
    out = OUT / f"{path.stem}_diagnostic.png"
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=110)
    print(f"  saved {out}")
    if "--show" in sys.argv:
        plt.show()


if __name__ == "__main__":
    main()
