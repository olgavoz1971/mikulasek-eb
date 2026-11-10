"""Look at what the package fits to SV Cen sector 64, before choosing a model.

Diagnostic only, no result. It fits the first candidate of the search (profile 1,
harmonics 1, 2, 3, baseline order 0) to the cleaned sector 64 file, with the constant
sigma from the phase-bin scatter, twice: (a) with the package defaults (the run that
ends with status ``bound-active``), and (b) with ``separate_shape_chi2_frac`` set to
1.0 so that the shared eclipse shape of the earlier stage is kept. Both fits are
drawn with the figure function of the tiered script.

Run from ``examples``::

    python sv_cen_model_look.py

Option ``--max-harmonic N`` uses harmonics 1 to N and ``--tag`` renames the folder.

Writes ``results/sv_cen/s64_look<tag>/{a_defaults,b_shared_shape}/fits/*.png``.
"""

from __future__ import annotations

import argparse
import dataclasses
import logging
from pathlib import Path

import mikulasek_eb as eb

import nsv807_s30_tiered as tiered

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "sv_cen"
FILE = RESULTS / "cleaned" / "TIC324101484_s64_QLP_200s.vot"
BINS = 400
HARMONICS = (1, 2, 3)


def main() -> None:
    """Fit the two variants and draw them; the figures are written even when the
    status is not ``ok``."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--max-harmonic", type=int, default=max(HARMONICS))
    parser.add_argument("--profiles", default="1,2,3,4")
    parser.add_argument("--tag", default="")
    args = parser.parse_args()
    harmonics = tuple(range(1, args.max_harmonic + 1))
    profiles = tuple(int(v) for v in args.profiles.split(","))
    unsupported = [p for p in profiles if p > 4]
    if unsupported:
        raise SystemExit(
            f"run_staged_fit supports profiles 1 to 4 (they have Gamma); got {unsupported}"
        )
    tiered.FLUX, tiered.FLUX_ERR, tiered.OUTLIER_BINS = "sap_flux", None, BINS
    payload = tiered.constant_sigma_from_bins(tiered.load_sap(FILE), FILE.stem)
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    features = eb.extract_features(
        t_mjd, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
    )
    base = eb.propose_candidates(features)[0]
    variants = {}
    for profile in profiles:
        cand = tiered.make_candidate(base, profile, harmonics, 0)
        variants.update(
            {
                "a_defaults": cand.config,
                "b_shared_shape": dataclasses.replace(
                    cand.config,
                    policy=dataclasses.replace(
                        cand.config.policy, separate_shape_chi2_frac=1.0
                    ),
                ).validate(),
            }
        )
    for name, config in variants.items():
        out = eb.fit(
            t_mjd,
            payload["mag"],
            payload["sig"],
            payload["period"],
            payload["epoch_jd"],
            config,
            data_id=f"{FILE.stem}__p{profile}_{name}",
        )
        status = out.solution.stats["status"]
        rms = 1e3 * float(out.solution.stats["rms_mag"])
        logger.info("%s: status %s rms %.2f mmag", name, status, rms)
        tiered.save_figures(
            RESULTS / f"s64_look{args.tag}" / f"p{profile}_{name}",
            name,
            payload,
            out,
            f"SV Cen s64 p{profile} {name}: status {status}, rms {rms:.1f} mmag",
        )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
