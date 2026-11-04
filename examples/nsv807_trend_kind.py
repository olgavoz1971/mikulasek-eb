"""NSV 807, QLP sector 97: is the instrumental trend additive or multiplicative in flux?

The test uses the eclipses themselves. For every orbital cycle the mean flux is
measured in four windows: the primary minimum, the secondary minimum and the two
quadratures (the maxima). Let ``Fmax`` be the mean of the quadratures and ``Fmin``
the flux of a minimum.

* A trend multiplicative in flux scales the whole cycle: ``Fmin = r * Fmax`` with
  a constant ratio ``r``, so a straight line through the cycles passes through
  the origin (intercept 0).
* A trend additive in flux shifts the whole cycle: ``Fmin = Fmax - A`` with a
  constant amplitude ``A``, so the slope is 1.

The slope and the intercept of ``Fmin`` against ``Fmax`` therefore separate the
two. Real changes of the star (spots, O'Connell effect) change ``r`` or ``A``
as well, and are not separated from the trend by this test. Only the file
is read: time and flux. The flux errors of this file are all NaN, so they are not
used. No fit of the package is involved.

Run from the repository root::

    MPLBACKEND=Agg python examples/nsv807_trend_kind.py
    MPLBACKEND=Agg python examples/nsv807_trend_kind.py <file.vot> <folder>

Writes ``results/nsv807_s97_trend/``.
"""

from __future__ import annotations

import csv
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import fitting, models
from volightcurve import read_lightcurve

logger = logging.getLogger(__name__)
DEFAULT_SOURCE = Path("/home/voz/Downloads/NSV_807_QLP_sector_97.vot")
RESULTS = Path(__file__).resolve().parents[1] / "results"
WINDOW = 0.05
MIN_POINTS = 8
N_BOOTSTRAP = 1000
SEED = 20261104
PHASE_BINS = 40


def load(path: Path):
    """Read time, flux, period and epoch with volightcurve.

    Args:
        path (Path): VOTable path.

    Returns:
        tuple: ``(t, flux, period, epoch)`` as float arrays and floats, in the
        time system of the file.

    Raises:
        FileNotFoundError: If the file is missing.
        ValueError: If the period or epoch are missing or rows are not finite.
    """
    if not path.is_file():
        raise FileNotFoundError(path)
    volc = read_lightcurve(path, filename=path.name, strict=True)
    time_cols = volc.get_time_colnames()
    flux_cols = volc.get_flux_colnames()
    if len(time_cols) != 1 or len(flux_cols) != 1:
        raise ValueError(f"need one time and one flux column: {time_cols} {flux_cols}")
    t = np.asarray(volc.table[time_cols[0]], float)
    flux = np.asarray(volc.table[flux_cols[0]], float)
    ok = np.isfinite(t) & np.isfinite(flux)
    if not np.all(ok):
        raise ValueError(f"{int((~ok).sum())} rows have a non-finite time or flux")
    period = volc.table.meta.get("period")
    epoch = volc.table.meta.get("epoch")
    if period is None or epoch is None:
        raise ValueError("period or epoch missing in the file")
    return t, flux, float(period), float(epoch)


def primary_phase_shift(t, flux, period, epoch) -> float:
    """Phase offset that puts the deeper minimum at phase 0.

    Args:
        t (ndarray): Times.
        flux (ndarray): Flux.
        period (float): Period (d).
        epoch (float): Epoch of the file.

    Returns:
        float: Phase of the deepest bin centre, in cycles.
    """
    phase = np.mod((t - epoch) / period, 1.0)
    edges = np.linspace(0, 1, PHASE_BINS + 1)
    index = np.clip(np.digitize(phase, edges) - 1, 0, PHASE_BINS - 1)
    medians = np.array([np.median(flux[index == i]) for i in range(PHASE_BINS)])
    centres = 0.5 * (edges[:-1] + edges[1:])
    return float(centres[int(np.argmin(medians))])


def per_cycle(t, flux, period, epoch, shift) -> list[dict]:
    """Mean flux of the minima and quadratures in every orbital cycle.

    Args:
        t (ndarray): Times.
        flux (ndarray): Flux.
        period (float): Period (d).
        epoch (float): Epoch of the file.
        shift (float): Phase of the primary minimum.

    Returns:
        list[dict]: One row per cycle with all four windows populated.
    """
    cycle_float = (t - epoch) / period - shift
    centres = {"min1": 0.0, "q1": 0.25, "min2": 0.5, "q2": 0.75}
    windows = {}
    for name, centre in centres.items():
        number = np.floor(cycle_float - centre + 0.5).astype(int)
        offset = cycle_float - centre - number
        inside = np.abs(offset) < WINDOW
        windows[name] = (number, inside)
    rows = []
    for number in np.unique(windows["min1"][0]):
        values, times, ok = {}, {}, True
        for name, (index, inside) in windows.items():
            sel = inside & (index == number)
            if sel.sum() < MIN_POINTS:
                ok = False
                break
            values[name] = float(np.mean(flux[sel]))
            times[name] = float(np.mean(t[sel]))
        if ok:
            rows.append(
                dict(
                    cycle=int(number),
                    t=float(np.mean(list(times.values()))),
                    min1=values["min1"],
                    min2=values["min2"],
                    fmax=0.5 * (values["q1"] + values["q2"]),
                )
            )
    return rows


def line_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Least-squares line ``y = a + b x``.

    Args:
        x (ndarray): Abscissa.
        y (ndarray): Ordinate.

    Returns:
        tuple[float, float]: Intercept and slope.
    """
    fitted = fitting.LinearLSQFitter()(models.Linear1D(1.0, 0.0), x, y)
    return float(fitted.intercept.value), float(fitted.slope.value)


def bootstrap(x: np.ndarray, y: np.ndarray) -> dict:
    """Line fit with bootstrap spread over cycles.

    Args:
        x (ndarray): Abscissa, one value per cycle.
        y (ndarray): Ordinate.

    Returns:
        dict: Intercept, slope and their standard deviations.
    """
    rng = np.random.default_rng(SEED)
    a0, b0 = line_fit(x, y)
    draws = []
    for _ in range(N_BOOTSTRAP):
        pick = rng.integers(0, len(x), len(x))
        draws.append(line_fit(x[pick], y[pick]))
    draws = np.asarray(draws)
    return dict(
        intercept=a0,
        intercept_err=float(draws[:, 0].std()),
        slope=b0,
        slope_err=float(draws[:, 1].std()),
    )


def main() -> None:
    """Measure the cycle fluxes, fit the lines and save the figure and tables.

    Optional arguments: the VOTable path and the output folder name under
    ``results/``. Without them the NSV 807 file is read.
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    SOURCE = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_SOURCE
    OUT = RESULTS / (sys.argv[2] if len(sys.argv) > 2 else "nsv807_s97_trend")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "command.txt").write_text(" ".join(sys.argv) + "\n")
    t, flux, period, epoch = load(SOURCE)
    shift = primary_phase_shift(t, flux, period, epoch)
    rows = per_cycle(t, flux, period, epoch, shift)
    logger.info(
        "%d points, period %.7f d, primary at phase %.3f, %d full cycles",
        len(t),
        period,
        shift,
        len(rows),
    )
    cycles = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    results = []
    for name in ("min1", "min2"):
        fit = bootstrap(cycles["fmax"], cycles[name])
        ratio = float(np.median(cycles[name] / cycles["fmax"]))
        fit.update(
            minimum=name,
            median_ratio=ratio,
            fmax_range_percent=100
            * float(np.ptp(cycles["fmax"]) / np.median(cycles["fmax"])),
            n_cycles=len(rows),
        )
        results.append(fit)
        logger.info(
            "%s: slope %.3f +- %.3f, intercept %.4f +- %.4f, median Fmin/Fmax %.3f, Fmax range %.1f %%",
            name,
            fit["slope"],
            fit["slope_err"],
            fit["intercept"],
            fit["intercept_err"],
            ratio,
            fit["fmax_range_percent"],
        )
    fields = [
        "minimum",
        "n_cycles",
        "slope",
        "slope_err",
        "intercept",
        "intercept_err",
        "median_ratio",
        "fmax_range_percent",
    ]
    with open(OUT / "line_fits.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fields)
        writer.writeheader()
        writer.writerows(results)
    with open(OUT / "cycles.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(3, 1, figsize=(11, 12))
    axes[0].plot(t, flux, ".", ms=1, color="0.7", label="all points")
    axes[0].plot(
        cycles["t"], cycles["fmax"], "o", ms=3, color="tab:green", label="maxima"
    )
    axes[0].plot(
        cycles["t"], cycles["min1"], "o", ms=3, color="tab:red", label="primary"
    )
    axes[0].plot(
        cycles["t"], cycles["min2"], "o", ms=3, color="tab:blue", label="secondary"
    )
    axes[0].set_xlabel("time (file units)")
    axes[0].set_ylabel("flux")
    axes[0].legend(loc="lower right")
    axes[0].set_title("NSV 807, QLP sector 97: flux of the cycles")
    span = np.linspace(cycles["fmax"].min(), cycles["fmax"].max(), 50)
    for axis, name, colour in (
        (axes[1], "min1", "tab:red"),
        (axes[2], "min2", "tab:blue"),
    ):
        y = cycles[name]
        fit = next(r for r in results if r["minimum"] == name)
        axis.plot(cycles["fmax"], y, "o", ms=3, color=colour)
        axis.plot(
            span,
            fit["intercept"] + fit["slope"] * span,
            "-",
            color="k",
            label="fitted line",
        )
        ratio = fit["median_ratio"]
        axis.plot(
            span,
            ratio * span,
            "--",
            color="tab:orange",
            label="multiplicative: through the origin",
        )
        mid = np.median(cycles["fmax"])
        axis.plot(
            span,
            np.median(y) + (span - mid),
            ":",
            color="tab:purple",
            label="additive: slope 1",
        )
        axis.set_xlabel("flux at maxima")
        axis.set_ylabel(
            f"flux at {'primary' if name == 'min1' else 'secondary'} minimum"
        )
        axis.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "cycle_fluxes.png", dpi=100)
    plt.close(fig)
    logger.info("wrote %s", OUT)


if __name__ == "__main__":
    main()
