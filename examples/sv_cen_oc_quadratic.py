"""O-C parabola of SV Cen from the sector-specific window timings.

The observed times are the collected timings in ``results/sv_cen/timing_dat``.
The reference period and the deep-minimum epoch are the joint fit with the
shape frozen. That epoch falls on the secondary of these files, so the primary
epoch used here is half a period earlier. Primary and secondary are fitted
separately:

    O-C = c + b E + a E^2,    Pdot = 2 a / P.

A time is rejected only when its residual from a first parabola, fitted to
every timed eclipse, lies more than five robust standard deviations from the
residual median. Rejected times are marked and left out of the reported fit.
A gap or wall flag is not a rejection.

Run from ``mikulasek_eb`` with ``MPLBACKEND=Agg``::

    python examples/sv_cen_oc_quadratic.py

Writes ``results/sv_cen/oc_quadratic/``: ``oc.png``, ``result.json``, ``oc.csv``.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.modeling import fitting, models
from scipy.stats import median_abs_deviation

import mikulasek_eb as eb

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "sv_cen"
TIMING = RESULTS / "timing_dat"
EPHEMERIS = RESULTS / "pdot_frozen" / "result.json"
OUT = RESULTS / "oc_quadratic"
KINDS = ("primary", "secondary")
SECONDARY_CYCLE = 0.5
OUTLIER_SIGMA = 5.0
PHASE_LIMIT = 0.1


def reference() -> tuple[float, float]:
    """Read the joint-fit period and the primary epoch.

    Returns:
        tuple[float, float]: Primary epoch (MJD) and period (days). The stored
        epoch is a deep minimum, half a cycle from the primary of the timings.
    """
    saved = json.loads(EPHEMERIS.read_text())
    period = float(saved["fitted"]["P0"])
    deep = float(saved["fitted"]["M0"])
    return deep - SECONDARY_CYCLE * period, period


def read_times(kind: str, epoch: float, period: float) -> dict[str, np.ndarray]:
    """Read one minimum and count its cycles on the reference ephemeris.

    Args:
        kind (str): ``primary`` or ``secondary``.
        epoch (float): Primary epoch, MJD.
        period (float): Reference period, days.

    Returns:
        dict: ``sector``, ``cycle``, ``oc`` (days), ``sigma`` (days), ``flag``.

    Raises:
        RuntimeError: A time is more than ``PHASE_LIMIT`` cycles from its grid
        point, so the cycle count is not safe.
    """
    half = 0.0 if kind == "primary" else SECONDARY_CYCLE
    sector, cycle, oc, sigma, flag = [], [], [], [], []
    for path in sorted(TIMING.glob(f"SVCen_*_{kind}.dat")):
        number = int(path.name.split("_")[1][1:])
        for line in path.read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            jd_text, err_text, mark = line.split()
            mjd = float(eb.absolute_jd_to_mjd(float(jd_text)))
            raw = (mjd - epoch) / period - half
            count = round(raw) + half
            if abs(raw - round(raw)) > PHASE_LIMIT:
                raise RuntimeError(
                    f"{kind} sector {number} is {raw - round(raw):.3f} cycles "
                    "from the nearest grid point"
                )
            sector.append(number)
            cycle.append(count)
            oc.append(mjd - (epoch + period * count))
            sigma.append(float(err_text))
            flag.append(mark)
    if not cycle:
        raise RuntimeError(f"no {kind} times in {TIMING}")
    return {
        "sector": np.asarray(sector, int),
        "cycle": np.asarray(cycle, float),
        "oc": np.asarray(oc, float),
        "sigma": np.asarray(sigma, float),
        "flag": np.asarray(flag),
    }


def fit_parabola(cycle: np.ndarray, oc: np.ndarray, sigma: np.ndarray) -> dict:
    """Fit ``c + b u + a u^2`` with weights ``1/sigma``.

    Args:
        cycle (ndarray): Cycle number.
        oc (ndarray): O-C in days.
        sigma (ndarray): Timing error in days.

    Returns:
        dict: Quadratic coefficient ``a`` (days per cycle squared), its formal
        and residual-rescaled errors, the reduced chi-square, and a callable
        of the parabola on the original cycle number.

    Raises:
        RuntimeError: The weighted normal matrix cannot be solved, or the
        astropy solution and the normal equation disagree.
    """
    origin = float(np.mean(cycle))
    u = cycle - origin
    design = np.column_stack((np.ones_like(u), u, u**2))
    weight = 1.0 / sigma**2
    normal = design.T @ (design * weight[:, None])
    if np.linalg.matrix_rank(normal) < 3:
        raise RuntimeError("the parabola normal matrix is rank deficient")
    coeff = np.linalg.solve(normal, design.T @ (weight * oc))
    fitted = fitting.LinearLSQFitter()(
        models.Polynomial1D(degree=2), u, oc, weights=1.0 / sigma
    )
    astropy_coeff = np.array([fitted.c0.value, fitted.c1.value, fitted.c2.value])
    if not np.allclose(coeff, astropy_coeff, rtol=1e-8, atol=1e-12):
        raise RuntimeError("the parabola solution and the astropy fit disagree")
    residual = oc - design @ coeff
    chi2 = float(np.sum(weight * residual**2))
    dof = len(cycle) - 3
    if dof < 1:
        raise RuntimeError("fewer than four times remain for the parabola")
    scale = chi2 / dof
    covariance = np.linalg.inv(normal)
    return {
        "a": float(coeff[2]),
        "error_formal": float(np.sqrt(covariance[2, 2])),
        "error_rescaled": float(np.sqrt(covariance[2, 2] * scale)),
        "chi2_reduced": float(scale),
        "origin": origin,
        "coeff": coeff,
        "residual": residual,
    }


def parabola(cycle: np.ndarray, solution: dict) -> np.ndarray:
    """Evaluate a fitted parabola at cycle numbers.

    Args:
        cycle (ndarray): Cycle number.
        solution (dict): Result of ``fit_parabola``.

    Returns:
        ndarray: O-C in days.
    """
    u = cycle - solution["origin"]
    c0, c1, c2 = solution["coeff"]
    return c0 + c1 * u + c2 * u**2


def per_year(rate: float) -> float:
    """Convert a dimensionless ``dP/dt`` to seconds per year.

    Args:
        rate (float): ``dP/dt`` with both ``P`` and ``t`` in days.

    Returns:
        float: The same rate in seconds per year.
    """
    return float((rate * (u.day / u.day)).to(u.s / u.yr).value)


def to_seconds(days: np.ndarray | float) -> np.ndarray | float:
    """Convert days to seconds.

    Args:
        days (ndarray or float): Interval in days.

    Returns:
        ndarray or float: The same interval in seconds.
    """
    return (np.asarray(days, float) * u.day).to(u.s).value


def outliers(residual: np.ndarray) -> np.ndarray:
    """Mark residuals that stand away from the first parabola.

    Args:
        residual (ndarray): Residual of the first parabola, days.

    Returns:
        ndarray: Boolean mask, true where a time is rejected.
    """
    seconds = to_seconds(residual)
    centre = float(np.median(seconds))
    scatter = float(median_abs_deviation(seconds, scale="normal"))
    if scatter == 0.0:
        raise RuntimeError("the parabola residuals have no scatter")
    return np.abs(seconds - centre) > OUTLIER_SIGMA * scatter


def summarise(
    kind: str, data: dict, kept: np.ndarray, solution: dict, period: float
) -> dict:
    """Collect the reported period change for one minimum.

    Args:
        kind (str): ``primary`` or ``secondary``.
        data (dict): Times from ``read_times``.
        kept (ndarray): Boolean mask of the times used in the reported fit.
        solution (dict): Parabola fitted to the kept times.
        period (float): Reference period, days.

    Returns:
        dict: Counts, residual size and the period change.
    """
    residual_s = to_seconds(solution["residual"])
    factor = 2.0 / period
    pdot = factor * solution["a"]
    return {
        "kind": kind,
        "n": int(kept.sum()),
        "n_rejected": int((~kept).sum()),
        "rejected_cycles": [float(c) for c in data["cycle"][~kept]],
        "n_flagged_kept": int(np.sum((data["flag"] != "ok") & kept)),
        "rms_s": float(np.sqrt(np.mean(residual_s**2))),
        "robust_sigma_s": float(median_abs_deviation(residual_s, scale="normal")),
        "largest_residual_s": float(np.max(np.abs(residual_s))),
        "chi2_reduced": solution["chi2_reduced"],
        "quadratic_days": solution["a"],
        "pdot_s_per_year": per_year(pdot),
        "pdot_error_formal_s_per_year": per_year(factor * solution["error_formal"]),
        "pdot_error_rescaled_s_per_year": per_year(factor * solution["error_rescaled"]),
    }


def draw(panels: list[dict], path: Path) -> None:
    """Draw the two O-C diagrams with the reported parabola.

    Args:
        panels (list[dict]): One entry per minimum, with the times, the kept
        mask and the reported solution.
        path (Path): Figure path.
    """
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True, sharey=True)
    for axis, panel in zip(axes, panels):
        data, kept, solution = panel["data"], panel["kept"], panel["solution"]
        cycle = data["cycle"]
        oc_s = to_seconds(data["oc"])
        err_s = to_seconds(data["sigma"])
        flagged = kept & (data["flag"] != "ok")
        plain = kept & ~flagged
        axis.errorbar(
            cycle[plain],
            oc_s[plain],
            err_s[plain],
            fmt="o",
            ms=3.5,
            color="0.2",
            ecolor="0.6",
            elinewidth=0.6,
            label="timed",
            zorder=2,
        )
        if np.any(flagged):
            axis.errorbar(
                cycle[flagged],
                oc_s[flagged],
                err_s[flagged],
                fmt="s",
                ms=4,
                color="0.2",
                markerfacecolor="none",
                ecolor="0.6",
                elinewidth=0.6,
                label="gap or wall, kept",
                zorder=3,
            )
        if np.any(~kept):
            axis.scatter(
                cycle[~kept],
                oc_s[~kept],
                marker="x",
                s=36,
                color="C3",
                label="rejected",
                zorder=4,
            )
        grid = np.linspace(cycle.min(), cycle.max(), 400)
        axis.plot(
            grid,
            to_seconds(parabola(grid, solution)),
            color="C3",
            lw=1.4,
            label="parabola",
            zorder=1,
        )
        axis.axhline(0.0, color="0.7", lw=0.6)
        rate = panel["summary"]["pdot_s_per_year"]
        error = panel["summary"]["pdot_error_rescaled_s_per_year"]
        axis.set_ylabel("O-C (s)")
        axis.set_title(
            f"{panel['kind'].capitalize()}, period change {rate:.2f} ± {error:.2f} s/yr"
        )
        axis.legend(frameon=False, fontsize=8)
    axes[1].set_xlabel("Cycle number")
    fig.suptitle("SV Cen, sector timings against the joint-fit period")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def write_table(panels: list[dict], path: Path) -> None:
    """Write one row per timed eclipse.

    Args:
        panels (list[dict]): Fitted minima.
        path (Path): CSV path.
    """
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["kind", "sector", "cycle", "oc_s", "sigma_s", "flag", "used", "residual_s"]
        )
        for panel in panels:
            data, kept = panel["data"], panel["kept"]
            residual = np.full(len(kept), np.nan)
            residual[kept] = to_seconds(panel["solution"]["residual"])
            for i in range(len(kept)):
                writer.writerow(
                    [
                        panel["kind"],
                        int(data["sector"][i]),
                        data["cycle"][i],
                        f"{to_seconds(data['oc'][i]):.3f}",
                        f"{to_seconds(data['sigma'][i]):.3f}",
                        data["flag"][i],
                        int(kept[i]),
                        "" if not kept[i] else f"{residual[i]:.3f}",
                    ]
                )


def main() -> None:
    """Fit both minima and write the figure and the numbers."""
    epoch, period = reference()
    panels = []
    for kind in KINDS:
        data = read_times(kind, epoch, period)
        first = fit_parabola(data["cycle"], data["oc"], data["sigma"])
        kept = ~outliers(first["residual"])
        solution = fit_parabola(
            data["cycle"][kept], data["oc"][kept], data["sigma"][kept]
        )
        summary = summarise(kind, data, kept, solution, period)
        panels.append(
            {
                "kind": kind,
                "data": data,
                "kept": kept,
                "solution": solution,
                "summary": summary,
            }
        )
        logger.info(
            "%s: Pdot %.3f (%.3f rescaled, %.4f formal) s/yr, n %d, rejected %d, "
            "residual rms %.1f s, largest %.1f s",
            kind,
            summary["pdot_s_per_year"],
            summary["pdot_error_rescaled_s_per_year"],
            summary["pdot_error_formal_s_per_year"],
            summary["n"],
            summary["n_rejected"],
            summary["rms_s"],
            summary["largest_residual_s"],
        )
    OUT.mkdir(parents=True, exist_ok=True)
    draw(panels, OUT / "oc.png")
    write_table(panels, OUT / "oc.csv")
    saved = {
        "primary_epoch_mjd": epoch,
        "period_days": period,
        "ephemeris": str(EPHEMERIS),
        "pdot": "2 a / P from O-C = c + b E + a E^2, a in days",
        "error_quoted": "rescaled by the reduced chi-square of the timing errors",
        "outlier_rule": f"{OUTLIER_SIGMA:g} robust standard deviations of the first parabola",
        "primary": panels[0]["summary"],
        "secondary": panels[1]["summary"],
    }
    (OUT / "result.json").write_text(json.dumps(saved, indent=2) + "\n")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
