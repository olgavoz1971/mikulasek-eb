"""Measure the white noise of the real TESS light curves in ``data/`` against cadence and magnitude.

For each file the script reports the cadence, the magnitude of the star from the file's own
photometric calibration, the white-noise sigma per point, that sigma scaled to one hour, and the
sigma of one hour expected from the published TESS noise model at that magnitude
(sigma(1 hr) = a + b 10^(0.2 (T - 10)) + c 10^(0.4 (T - 10)) ppm, median star, a = 50.2,
b = 97.4, c = 92.9, arXiv:2202.03656, as quoted in ``docs/synthetic_noise.md``).

The noise is ``sampling["sigma_noise_mag"]`` of ``extract_features``: the robust sigma of
consecutive differences of the residuals about the folded curve, so the eclipses are removed.
Magnitudes come from ``load_tess_votable`` (volightcurve photometric calibration). The magnitude of
a star is taken as the 10th percentile of its magnitudes (the out-of-eclipse level). This is the
magnitude of the flux column in the file, not the TIC magnitude.

Run::

    cd /home/voz/projects/UPJS/mikulasek_eb
    MPLBACKEND=Agg python examples/real_noise_vs_cadence.py
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mikulasek_eb import extract_features
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable

logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "real_noise"

STARS = {
    "AA And S16 2-min": "lc_discovery_TIC_455790537_SPOC_sector_16",
    "AA And S16 30-min": "lc_discovery_TIC_455790537_TESS_SPOC_sector_16",
    "AA And S17 30-min": "lc_discovery_TIC_455790537_TESS_SPOC_sector_17",
    "XY Boo S50 10-min": "lc_discovery_TIC_347405567_TESS_SPOC_sector_50",
    "EK Com S49": "EK_Com_sector_49",
    "V477 Lyr S26": "v477_Lyr_sector_26",
    "AV Del S81": "AV_Del_sector_81",
    "AR Aur S19": "AR_Aur_sector_19",
    "NSV 807 S04": "NSV_807_sector_04",
}

MAG_TO_FLUX = 1.0 / 1.0857  # sigma_mag * this = relative flux sigma
MODEL_A, MODEL_B, MODEL_C = 50.2, 97.4, 92.9  # ppm, median star, arXiv:2202.03656


def model_sigma_1hr_ppm(tmag):
    """Published sigma of one hour for a 2-min TESS target of the median star.

    Args:
        tmag (float or ndarray): TESS magnitude.

    Returns:
        float or ndarray: Noise of one hour (ppm).
    """
    x = np.asarray(tmag, float) - 10.0
    return MODEL_A + MODEL_B * 10 ** (0.2 * x) + MODEL_C * 10 ** (0.4 * x)


def measure(label: str, stem: str) -> dict:
    """Measure cadence, magnitude and white noise of one light curve.

    Args:
        label (str): Short name of the star and sector.
        stem (str): File stem in ``data/``.

    Returns:
        dict: One row of the table.
    """
    payload = load_tess_votable(ROOT / "data" / f"{stem}.vot")
    t = absolute_jd_to_mjd(payload["t"])
    f = extract_features(
        t, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    s = f.sampling
    cadence_s = float(s["cadence_s"])
    tmag = float(np.percentile(payload["mag"], 10))
    sigma_mag = float(s["sigma_noise_mag"])
    sigma_ppm = 1e6 * sigma_mag * MAG_TO_FLUX
    sigma_1hr = sigma_ppm * np.sqrt(cadence_s / 3600.0)
    expected = float(model_sigma_1hr_ppm(tmag))
    expected_pt = expected * np.sqrt(3600.0 / cadence_s)
    return {
        "star": label,
        "cadence_s": round(cadence_s, 1),
        "mag_file": round(tmag, 3),
        "sigma_mmag": round(1e3 * sigma_mag, 3),
        "sigma_ppm": round(sigma_ppm, 0),
        "sigma_1hr_ppm": round(sigma_1hr, 0),
        "model_1hr_ppm": round(expected, 0),
        "model_per_point_ppm": round(expected_pt, 0),
        "measured_over_model": round(sigma_1hr / expected, 2),
        "catalogue_sigma_mmag": round(1e3 * float(s["sigma_catalogue_mag"]), 3),
        "noise_over_catalogue": round(float(s["noise_over_catalogue"]), 2),
    }


def draw(rows: list[dict], path: Path) -> None:
    """Plot the measured sigma of one hour against the magnitude and the published model.

    Args:
        rows (list[dict]): Rows from ``measure``.
        path (Path): Output PNG.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    grid = np.linspace(5, 14, 200)
    ax.plot(grid, model_sigma_1hr_ppm(grid), color="#444444", lw=2, label="TESS model, median star")
    ax.fill_between(
        grid,
        model_sigma_1hr_ppm(grid) / 2,
        model_sigma_1hr_ppm(grid) * 2,
        color="#999999",
        alpha=0.25,
        label="factor 2 either side",
    )
    colours = {120: "#d62728", 200: "#9467bd", 600: "#2ca02c", 1800: "#1f77b4"}
    for r in rows:
        c = colours.get(int(round(r["cadence_s"])), "#000000")
        ax.scatter(r["mag_file"], r["sigma_1hr_ppm"], s=70, color=c, zorder=3)
        ax.annotate(r["star"], (r["mag_file"], r["sigma_1hr_ppm"]), fontsize=8,
                    xytext=(5, 4), textcoords="offset points")
    for cad, c in colours.items():
        ax.scatter([], [], color=c, label=f"{cad} s")
    ax.axhline(438, color="#d62728", ls=":", lw=1)
    ax.text(5.1, 460, "Tmag 11 in docs/synthetic_noise.md: 438 ppm", fontsize=8, color="#d62728")
    ax.set_yscale("log")
    ax.set_xlabel("Magnitude of the flux column in the file (10th percentile)")
    ax.set_ylabel("White noise scaled to one hour (ppm)")
    ax.set_title("Real TESS light curves: measured noise against the published model")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    """Measure every file, write the table and the figure."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [measure(label, stem) for label, stem in STARS.items()]
    with (OUT / "noise_vs_cadence.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    draw(rows, OUT / "noise_vs_cadence.png")
    for r in rows:
        logger.info("%s", r)
    logger.info("wrote %s", OUT)


if __name__ == "__main__":
    main()
