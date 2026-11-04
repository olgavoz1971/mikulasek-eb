"""Check the Tmag 11 synthetic light curves of LC_Generator against what was requested.

Reads every file of ``results/m11/manifest.csv`` and ``results/m11/manifest_spotsweep.csv`` of
the generator repository and compares it with the request (the cadence table, the Tmag 11 flux
scale, white flux-relative noise, the error column equal to the noise). For each file:

- the calibrated magnitude of the bright state (90th percentile of the noise-free smeared curve)
  through the file's own PhotCal (``volightcurve``), and the 10th percentile of the observed
  magnitudes, the estimate used for the real stars in ``real_noise_vs_cadence.py``;
- the noise actually present: ``flux - truth_flux_pure_smeared - truth_trend`` divided by the
  trended truth, its standard deviation against ``rel_sigma``, the median error column divided
  by the trended truth, and the lag-1 autocorrelation (white noise: about 0);
- the same white-noise estimate that was applied to the real files (``extract_features``,
  ``sigma_noise_mag``), scaled to one hour, against the published model at magnitude 11.

The spot sweep also gets one check that the noise realisation is shared across longitudes. A file
that cannot be processed is written to the table with the error text; nothing is skipped.

Run::

    cd /home/voz/projects/UPJS/mikulasek_eb
    MPLBACKEND=Agg python examples/check_m11_files.py
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from volightcurve import read_lightcurve

from mikulasek_eb import extract_features
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable

from real_noise_vs_cadence import MAG_TO_FLUX, model_sigma_1hr_ppm

logger = logging.getLogger(__name__)

GEN = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFESTS = ("results/m11/manifest.csv", "results/m11/manifest_spotsweep.csv")
OUT = Path(__file__).resolve().parent.parent / "results" / "real_noise"
REQUESTED = {120: 0.00240, 600: 0.00140, 1800: 0.00125}
REQUESTED_MMAG = {120: 2.6, 600: 1.5, 1800: 1.4}


def lag1(x: np.ndarray) -> float:
    """Lag-1 autocorrelation of a series.

    Args:
        x (ndarray): Values in time order.

    Returns:
        float: Correlation of consecutive points.
    """
    return float(np.corrcoef(x[:-1], x[1:])[0, 1])


def check(path: Path) -> tuple[dict, np.ndarray | None]:
    """Compare one file with the request.

    Args:
        path (Path): A VOTable in ``results/m11``.

    Returns:
        tuple: The table row (dict) and the noise in units of its own sigma (ndarray), or None.
    """
    volc = read_lightcurve(path, filename=path.name, strict=True)
    tab = volc.table
    meta = tab.meta
    photcal = volc.photdms["flux"].photcal
    flux, err = tab["flux"], tab["flux_err"]
    truth = np.asarray(tab["truth_flux_pure_smeared"], float)
    trend = np.asarray(tab["truth_trend"], float)
    base = truth + trend
    resid = np.asarray(flux, float) - base
    rel = resid / base
    cadence = float(meta["cadence_s"])
    key = int(round(cadence))
    bright = photcal.flux_to_mag(np.percentile(truth, 90) * flux.unit)
    mag_obs = np.percentile(photcal.flux_to_mag(flux).value, 10)
    err_rel = float(np.median(np.asarray(err, float) / base))
    sig_mmag = float(
        np.median(photcal.flux_err_to_mag_err(flux, err).to("mmag").value)
    )
    dt = np.diff(np.asarray(tab["time"], float))
    row = {
        "file": path.name,
        "cadence_s": cadence,
        "n": len(tab),
        "median_dt_s": round(float(np.median(dt)) * 86400, 1),
        "span_d": round(float(tab["time"][-1] - tab["time"][0]), 2),
        "bright_mag_truth": round(float(bright.value), 4),
        "mag_10pct_observed": round(float(mag_obs), 3),
        "rel_sigma_meta": meta["rel_sigma"],
        "rel_sigma_requested": REQUESTED.get(key),
        "rel_noise_measured": round(float(np.std(rel)), 5),
        "err_over_trended_truth": round(err_rel, 5),
        "lag1": round(lag1(resid), 3),
        "err_mmag": round(sig_mmag, 3),
        "requested_mmag": REQUESTED_MMAG.get(key),
        "noise_set": meta.get("noise_set"),
        "seed": meta.get("seed"),
    }
    try:
        payload = load_tess_votable(path)
        t = absolute_jd_to_mjd(payload["t"])
        f = extract_features(
            t, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
        )
        s = f.sampling
        sigma_ppm = 1e6 * float(s["sigma_noise_mag"]) * MAG_TO_FLUX
        sig_1hr = sigma_ppm * np.sqrt(cadence / 3600.0)
        row["estimator_mmag"] = round(1e3 * float(s["sigma_noise_mag"]), 3)
        row["estimator_1hr_ppm"] = round(sig_1hr, 0)
        row["model_1hr_ppm_tmag11"] = round(float(model_sigma_1hr_ppm(11.0)), 0)
        row["estimator_over_model"] = round(sig_1hr / float(model_sigma_1hr_ppm(11.0)), 2)
        row["error"] = ""
    except Exception as exc:  # reported in the table, not hidden
        row["error"] = f"{type(exc).__name__}: {exc}"
    return row, resid / (err_rel * base)


def main() -> None:
    """Check all files, write the table, the same-seed test and the figure."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    rows, noise = [], {}
    for man in MANIFESTS:
        with (GEN / man).open() as fh:
            for rec in csv.DictReader(fh):
                p = GEN / rec["path"]
                row, z = check(p)
                row["batch"] = "spotsweep" if "spotsweep" in man else "first"
                rows.append(row)
                noise[p.name] = z
    cols = list(dict.fromkeys(k for r in rows for k in r))
    with (OUT / "m11_check.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    # shared noise realisation across the spot sweep, 120 s
    names = sorted(n for n in noise if "_SP_" in n and n.endswith("_120_c120.vot"))
    ref = noise[names[0]]
    shared = [(n, float(np.corrcoef(ref, noise[n])[0, 1])) for n in names[1:]]
    with (OUT / "m11_same_noise.txt").open("w") as fh:
        fh.write(f"reference {names[0]}\n")
        for n, c in shared:
            fh.write(f"{c:.4f}  {n}\n")

    # figure: noise against the real files, and one light curve against its truth
    fig, (a, b) = plt.subplots(1, 2, figsize=(14, 5.5))
    grid = np.linspace(9, 13, 100)
    a.plot(grid, model_sigma_1hr_ppm(grid), color="#444", label="TESS model, median star")
    a.fill_between(grid, model_sigma_1hr_ppm(grid) / 2, model_sigma_1hr_ppm(grid) * 2,
                   color="#999", alpha=0.25, label="factor 2 either side")
    col = {120: "#d62728", 600: "#2ca02c", 1800: "#1f77b4"}
    for r in rows:
        if r["error"]:
            continue
        a.scatter(r["mag_10pct_observed"], r["estimator_1hr_ppm"], facecolors="none",
                  edgecolors=col[int(round(r["cadence_s"]))], s=40)
    real = list(csv.DictReader((OUT / "noise_vs_cadence.csv").open()))
    for r in real:
        a.scatter(float(r["mag_file"]), float(r["sigma_1hr_ppm"]), marker="s", color="k", s=30)
    for k, c in col.items():
        a.scatter([], [], facecolors="none", edgecolors=c, label=f"synthetic {k} s")
    a.scatter([], [], marker="s", color="k", label="real files")
    a.set_xlim(9, 13)
    a.set_yscale("log")
    a.set_xlabel("Magnitude (10th percentile of the observed magnitudes)")
    a.set_ylabel("White noise scaled to one hour (ppm)")
    a.set_title("Tmag 11 synthetic files against the real files")
    a.legend(fontsize=8)
    a.grid(alpha=0.3, which="both")

    volc = read_lightcurve(
        GEN / "results/m11/detached/detached_none_nooconnell_SP_L045_SP_L045_120_c120.vot",
        filename="x.vot", strict=True)
    tb = volc.table
    t = np.asarray(tb["time"], float)
    b.plot(t, np.asarray(tb["flux"], float), ".", ms=1.5, color="#888", label="observed flux")
    b.plot(t, np.asarray(tb["truth_flux_pure_smeared"], float) + np.asarray(tb["truth_trend"], float),
           "-", lw=0.7, color="#d62728", label="truth + trend")
    b.set_xlim(2, 6)
    b.set_xlabel("Time (d, file clock)")
    b.set_ylabel("Flux (ct/s)")
    b.set_title("SP_L045, 120 s, m11")
    b.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "m11_check.png", dpi=150)
    plt.close(fig)
    logger.info("wrote %s", OUT)


if __name__ == "__main__":
    main()
