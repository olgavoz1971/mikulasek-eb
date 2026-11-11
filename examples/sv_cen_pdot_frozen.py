"""One joint fit of SV Cen with the shape frozen and ``Pdot`` free.

The template is the sector 64 profile 3 solution, held fixed, and the exposure
integration is off, because the sectors have different cadences. The free parameters
are the epoch, the period and ``Pdot``, plus one constant baseline per stretch of
data (a new stretch starts at a gap longer than 0.5 d). The start comes from the
linear O-C elements, and ``Pdot`` starts at 0.

Run from ``examples``::

    python sv_cen_pdot_frozen.py

Writes ``results/sv_cen/pdot_frozen/``: ``result.json``, ``fit.png``.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import mikulasek_eb as eb
import nsv807_s30_tiered as tiered

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "sv_cen"
DATA = RESULTS / "cleaned"
SOLUTION = RESULTS / "s64_fit" / "chosen"
OUT = RESULTS / "pdot_frozen"
BINS = {"1800s": 70, "600s": 200, "200s": 400}
GAP_DAYS = 0.5
P0 = 1.6573037839324687
# The linear O-C epoch falls on the shallow minimum. The template is deeper at
# phase 0, so the start is moved by half a period to put that side on the deep
# minimum. The period and the zero of Pdot are the linear O-C values.
M0 = 51870.08157841349 + 0.5 * P0
PDOT = 0.0
SECONDS_PER_YEAR = 86400.0 * 365.25


def sectors() -> list[Path]:
    """Cleaned QLP files, one per sector, in sector order.

    Returns:
        list[Path]: VOTables, without the stand-in error copies.
    """
    files = [p for p in DATA.glob("*QLP*.vot") if "_witherr" not in p.name]
    return sorted(files, key=lambda p: int(p.stem.split("_")[1][1:]))


def load() -> dict:
    """Load every sector with its own constant sigma.

    Returns:
        dict: ``t`` (MJD), ``mag``, ``sig`` and ``seg`` (one id per stretch).
    """
    t_all, mag_all, sig_all = [], [], []
    for path in sectors():
        tiered.FLUX_ERR = None
        tiered.OUTLIER_BINS = BINS[path.stem.rsplit("_", 1)[1]]
        payload = tiered.constant_sigma_from_bins(tiered.load_sap(path), path.stem)
        t_all.append(eb.absolute_jd_to_mjd(payload["t"]))
        mag_all.append(payload["mag"])
        sig_all.append(payload["sig"])
    t = np.concatenate(t_all)
    mag = np.concatenate(mag_all)
    sig = np.concatenate(sig_all)
    order = np.argsort(t)
    t, mag, sig = t[order], mag[order], sig[order]
    seg = np.cumsum(np.r_[0, np.diff(t) > GAP_DAYS])
    for s in np.unique(seg):
        sel = seg == s
        logger.info(
            "stretch %d: %d points, MJD %.2f to %.2f",
            int(s),
            int(sel.sum()),
            float(t[sel][0]),
            float(t[sel][-1]),
        )
    return {"t": t, "mag": mag, "sig": sig, "seg": seg}


def frozen_model():
    """Sector 64 template with the exposure integration turned off.

    Returns:
        tuple: ``EBModel`` and its parameter dict. The exposure is omitted because
        one exposure cannot describe every cadence.
    """
    solution = eb.ModelSolution.load(next(SOLUTION.glob("*.json")))
    spec = solution.config.model
    model = eb.EBModel(
        profile=spec.profile,
        harmonics=spec.harmonics,
        ocon=spec.ocon,
        kterm=spec.kterm,
    )
    return model, dict(solution.params)


def per_year(pdot: float) -> float:
    """Convert a dimensionless ``Pdot`` to seconds per year.

    Args:
        pdot (float): ``dP/dt`` with ``P`` and ``t`` in days.

    Returns:
        float: The same rate in seconds per year.
    """
    return pdot * SECONDS_PER_YEAR


def draw(data: dict, model, params: dict, ephem: dict, path: Path) -> None:
    """Draw the folded fit and the residual against time.

    Args:
        data (dict): Arrays from ``load``.
        model: Frozen template.
        params (dict): Frozen template parameters.
        ephem (dict): Fitted ``M0``, ``P0`` and ``Pdot``.
        path (Path): Figure file.
    """
    theta = eb.phase(data["t"], ephem["M0"], ephem["P0"], ephem["Pdot"])
    template = np.asarray(model(theta, params), float)
    baseline = data["mag"] - template - data["resid"]
    phase = theta % 1.0
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    axes[0].scatter(phase, data["mag"] - baseline, s=1, color="0.6", label="data")
    grid = np.linspace(0, 1, 1000)
    axes[0].plot(
        grid,
        np.asarray(model(grid, params), float),
        color="red",
        lw=1,
        label="template",
    )
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Phase")
    axes[0].set_ylabel("Magnitude, baseline removed")
    axes[0].legend()
    axes[0].set_title(f"SV Cen, shape frozen, Pdot {per_year(ephem['Pdot']):.2f} s/yr")
    axes[1].scatter(data["t"], 1e3 * data["resid"], s=1, color="0.5")
    axes[1].axhline(0.0, color="red", lw=1)
    axes[1].set_xlabel("MJD")
    axes[1].set_ylabel("Residual (mmag)")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def main() -> None:
    """Fit and write the result."""
    data = load()
    model, params = frozen_model()
    ephem0 = {"M0": M0, "P0": P0, "Pdot": PDOT}
    result = eb.fit_global(
        data["t"],
        data["mag"],
        data["sig"],
        model,
        params,
        [],
        ephem0,
        free_ephem=("M0", "P0", "Pdot"),
        seg=data["seg"],
        order=0,
    )
    data["resid"] = result.resid
    OUT.mkdir(parents=True, exist_ok=True)
    draw(data, model, params, result.ephem, OUT / "fit.png")
    saved = {
        "start": ephem0,
        "fitted": result.ephem,
        "error": {name: result.err[name] for name in result.ephem},
        "pdot_s_per_year": per_year(result.ephem["Pdot"]),
        "pdot_error_s_per_year": per_year(result.err["Pdot"]),
        "rms_mmag": 1e3 * result.rms,
        "n_points": int(len(data["t"])),
        "n_stretches": int(len(np.unique(data["seg"]))),
        "nfev": result.nfev,
        "rank_deficient": bool(result.rank_deficient),
        "template": "sector 64 p3_h123_b1, shape frozen, no exposure integration",
    }
    (OUT / "result.json").write_text(json.dumps(saved, indent=2) + "\n")
    logger.info(
        "M0 %.5f (%.5f) d, P0 %.8f (%.8f) d, Pdot %.3f (%.3f) s/yr, rms %.2f mmag",
        result.ephem["M0"],
        result.err["M0"],
        result.ephem["P0"],
        result.err["P0"],
        saved["pdot_s_per_year"],
        saved["pdot_error_s_per_year"],
        saved["rms_mmag"],
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
