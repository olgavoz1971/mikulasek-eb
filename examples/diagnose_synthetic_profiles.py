"""Where does the shape error sit? Pictures of profile fits against the pure truth.

Fits one synthetic file (true trend removed, baseline order 0) with the chosen profiles,
draws for each: the phase-binned data, the fit and the pure smeared truth, the error
(fit minus truth) over the whole phase and zoomed on both eclipses, and prints the rms of
the error in phase zones (primary core, primary wings, secondary core, secondary wings,
out of eclipse). Pictures go to ``results/synthetic_scoring/profile_pictures/``.

Run::

    python examples/diagnose_synthetic_profiles.py G4_B08 [profile:harmonics:ocon ...]
    # optional 4th field ':2' forces the two-component trial (default variants: 1:123:0 3:123:0 3:123:1)
"""

from __future__ import annotations

import csv
import dataclasses
import logging
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from score_synthetic_fits import DEFAULT_RESULTS, OUT, build_library  # noqa: E402

from mikulasek_eb import fit  # noqa: E402
from mikulasek_eb.config import ModelSpec  # noqa: E402
from mikulasek_eb.features import extract_features  # noqa: E402
from mikulasek_eb.helpers import absolute_jd_to_mjd  # noqa: E402
from mikulasek_eb.io import load_tess_votable  # noqa: E402
from mikulasek_eb.synthetic_check import _num, load_synthetic  # noqa: E402

logger = logging.getLogger(__name__)
PICTURES = OUT / "profile_pictures"


def binned(ph, y, nb=200):
    """Mean of ``y`` in phase bins.

    Args:
        ph (ndarray): Phase in [0, 1).
        y (ndarray): Values.
        nb (int): Number of bins.

    Returns:
        tuple[ndarray, ndarray]: Bin centres and means (NaN for empty bins).
    """
    idx = np.minimum((ph * nb).astype(int), nb - 1)
    n = np.bincount(idx, minlength=nb)
    s = np.bincount(idx, weights=y, minlength=nb)
    return (np.arange(nb) + 0.5) / nb, np.where(n > 0, s / np.maximum(n, 1), np.nan)


def zone_rms(ph_c, err, half):
    """Rms of the error in phase zones around the eclipses.

    Args:
        ph_c (ndarray): Phase centred on the primary, in [-0.5, 0.5).
        err (ndarray): Fit minus truth (mag).
        half (float): Half-width of the eclipse core (phase).

    Returns:
        dict: Rms in mmag per zone.
    """
    d1 = np.abs(ph_c)
    d2 = np.abs(np.abs(ph_c) - 0.5)
    zones = {
        "primary core": d1 < half,
        "primary wings": (d1 >= half) & (d1 < 2.5 * half),
        "secondary core": d2 < half,
        "secondary wings": (d2 >= half) & (d2 < 2.5 * half),
    }
    zones["out of eclipse"] = ~(zones["primary core"] | zones["primary wings"] | zones["secondary core"] | zones["secondary wings"])
    return {k: 1e3 * float(np.sqrt(np.mean(err[m] ** 2))) for k, m in zones.items() if m.any()}


def main() -> None:
    """Fit the variants and write pictures and a zone table."""
    logging.basicConfig(level=logging.INFO if "-v" in sys.argv else logging.WARNING, format="%(message)s")
    key = sys.argv[1]
    sys.argv = [a for a in sys.argv if a != "-v"] if False else sys.argv
    specs = [a for a in sys.argv[2:] if a != "-v"] or ["1:123:0", "3:123:0", "3:123:1"]
    with open(DEFAULT_RESULTS / "manifest.csv", newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if key in r["path"]]
    path = DEFAULT_RESULTS.parent / rows[0]["path"]
    d = load_synthetic(path)
    payload = load_tess_votable(path)
    c = float(np.median(payload["mag"] + 2.5 * np.log10(d.flux)))
    mag = c - 2.5 * np.log10(d.flux - d.truth_trend)
    truth = c - 2.5 * np.log10(d.truth_flux)
    t = absolute_jd_to_mjd(payload["t"])
    feats = extract_features(t, mag, payload["sig"], payload["period"], payload["epoch_jd"])
    base = next(x for x in build_library(feats) if x.config.data.baseline_order == 0)
    PICTURES.mkdir(parents=True, exist_ok=True)
    half = 0.5 * max(feats.geometry["fwhm1"] or 0.03, feats.geometry["fwhm2"] or 0.03)
    results = []
    fig, axes = plt.subplots(3, len(specs), figsize=(6 * len(specs), 11), squeeze=False)
    for j, spec in enumerate(specs):
        prof, harm, oc, *force = spec.split(":")
        model = ModelSpec(profile=int(prof), harmonics=tuple(int(h) for h in harm), ocon=oc == "1")
        policy = base.config.policy
        if force and force[0] == "2":  # force the two-component trial
            policy = dataclasses.replace(policy, rms_noise_two_comp_threshold=0.0, two_comp_min_rms_gain_frac=0.01)
        cfg = dataclasses.replace(base.config, model=model, policy=policy).validate()
        out = fit(t, mag, payload["sig"], payload["period"], payload["epoch_jd"], cfg, data_id=path.stem)
        idx = np.searchsorted(t, out.t_mjd)
        fit_mag = out.mag - out.pipeline.adopted.resid
        err = fit_mag - truth[idx]
        e = out.solution.ephemeris
        ph = ((out.t_mjd - e["M0"]) / e["P0"]) % 1.0
        ph_c = (ph + 0.5) % 1.0 - 0.5
        order = np.argsort(ph_c)
        z = zone_rms(ph_c, err, half)
        results.append((spec, z, bool(out.solution.stats["two_component"])))
        bx, by = binned(ph_c + 0.5, out.mag)
        ax = axes[0, j]
        ax.plot(bx - 0.5, by, ".", ms=3, color="0.6", label="data (binned)")
        ax.plot(ph_c[order], truth[idx][order], "k-", lw=1, label="pure truth")
        ax.plot(ph_c[order], fit_mag[order], "r-", lw=1, label="fit")
        ax.invert_yaxis()
        ax.set_title(f"profile {prof}, harmonics {harm}, ocon {oc}")
        ax.set_xlabel("phase")
        ax.set_ylabel("mag")
        ax.legend(fontsize=8)
        ex, ey = binned(ph_c + 0.5, err)
        ax = axes[1, j]
        ax.plot(ex - 0.5, 1e3 * ey, "r.-", lw=0.8, ms=3)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_title("fit minus truth, binned (mmag)")
        ax.set_xlabel("phase")
        ax = axes[2, j]
        for centre, col in ((0.0, "C0"), (0.5, "C1")):
            m = np.abs(((ph_c - centre + 0.5) % 1.0) - 0.5) < 3 * half
            x = ((ph_c[m] - centre + 0.5) % 1.0) - 0.5
            o = np.argsort(x)
            ax.plot(x[o], 1e3 * err[m][o], ".", ms=2, color=col, label="primary" if centre == 0 else "secondary")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_title("error zoomed on the eclipses (mmag)")
        ax.set_xlabel("phase from eclipse centre")
        ax.legend(fontsize=8)
    fig.suptitle(path.name)
    fig.tight_layout()
    png = PICTURES / f"{path.stem}_profiles{'_' + '_'.join(s.replace(':', '') for s in specs) if sys.argv[2:] else ''}.png"
    fig.savefig(png, dpi=110)
    lines = [f"File {path.name}, eclipse core half-width {half:.4f} phase", "",
             "| variant | two-comp | " + " | ".join(results[0][1]) + " |", "|---|---|" + "---|" * len(results[0][1])]
    for spec, z, two in results:
        lines.append(f"| {spec} | {two} | " + " | ".join(f"{v:.2f}" for v in z.values()) + " |")
    (PICTURES / f"{path.stem}_zones.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(png)


if __name__ == "__main__":
    main()
