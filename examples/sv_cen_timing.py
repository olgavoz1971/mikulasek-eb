"""Window timing of SV Cen, the timing of ``nsv807_timing.py`` with one change.

The SV Cen QLP files have no SAP error column, and the pilot fit that
``nsv807_timing.prepare`` uses to get a constant sigma ends ``bound-active`` on them.
The tiered fits of these sectors (``nsv807_s30_tiered.py --sigma-from-bins``) took the
constant sigma from the scatter about the phase-bin medians instead. This script does
the same for the timing, and nothing else: every mode, the window rule and the point
rule are those of ``nsv807_timing.py``, which is not modified.

``--bins`` must be the number used for that sector's fit (400 at 200 s, 200 at 600 s,
70 at 1800 s), so that the timing sigma equals the fit sigma.

Run from ``examples``::

    python sv_cen_timing.py scan --run <fit folder> --file <cleaned.vot> --bins 400
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np

import nsv807_s30_tiered as tiered
import nsv807_timing as timing

logger = logging.getLogger(__name__)

# Kept before ``main`` replaces ``timing.prepare``, so the replacement can still
# call the original.
PREPARE = timing.prepare


def prepare(run: Path, name: str, bins: int, skip_outliers: bool = True):
    """Load one cleaned sector with the constant sigma from the phase-bin scatter.

    The cleaned files have no SAP error column. A copy with a stand-in error of one
    percent of the flux lets the loader run. The stand-in is thrown away: the sigma
    below replaces it, and no fit uses it. The pilot fit of ``nsv807_timing.prepare``
    is not used: it ends ``bound-active`` on these files.

    Args:
        run (Path): Tiered-fit folder, with the chosen solution in ``chosen``.
        name (str): Cleaned VOTable, absolute path or file name.
        bins (int): Phase bins, the same number as in the tiered fit of this sector.
        skip_outliers (bool): True when the tiered fit skipped the outlier step.

    Returns:
        tuple: ``payload``, ``t_mjd``, the chosen ``ModelSolution`` and the median
        cadence in days, as from ``nsv807_timing.prepare``.
    """
    src = Path(name)
    if src.stem.endswith("_witherr"):
        src = src.with_name(src.stem[: -len("_witherr")] + src.suffix)
    stand_in = src.with_name(src.stem + "_witherr.vot")
    if not stand_in.exists():
        from astropy.table import Column
        from volightcurve import read_lightcurve, write_lightcurve

        volc = read_lightcurve(src, filename=src.name, strict=True)
        flux = np.asarray(volc.table["sap_flux"], float)
        volc.table.add_column(
            Column(
                np.full_like(flux, np.nanmedian(flux) * 0.01),
                name="sap_flux_err",
                unit=volc.table["sap_flux"].unit,
            )
        )
        volc.table["sap_flux_err"].meta["ucd"] = "stat.error;phot.flux"
        write_lightcurve(volc, "votable_binary", destination=stand_in)
        logger.info("wrote %s with a stand-in sap_flux_err", stand_in.name)
    tiered.FLUX_ERR = "sap_flux_err"
    payload, t_mjd, solution, cadence = PREPARE(
        run, str(stand_in), "sap_flux", "sap_flux_err", skip_outliers
    )
    tiered.OUTLIER_BINS = bins
    tiered.FLUX_ERR = None
    measured = tiered.constant_sigma_from_bins(
        tiered.load_sap(Path(name)), Path(name).stem
    )
    payload["sig"] = measured["sig"]
    return payload, t_mjd, solution, cadence


def main() -> None:
    """Run one mode of ``nsv807_timing`` on an SV Cen sector."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("mode", choices=("scan", "time", "figure", "adopt"))
    parser.add_argument("--run", required=True, help="tiered-fit folder of the sector")
    parser.add_argument("--file", required=True, help="cleaned VOTable")
    parser.add_argument("--bins", type=int, required=True, help="phase bins of the fit")
    parser.add_argument(
        "--half-over-d",
        type=float,
        nargs=2,
        metavar=("PRIMARY", "SECONDARY"),
        help="window multiples of D for the time mode",
    )
    parser.add_argument(
        "--start-over-d",
        type=float,
        nargs=2,
        default=(0.0, 0.0),
        metavar=("PRIMARY", "SECONDARY"),
        help="adopt mode: first window considered, in D",
    )
    parser.add_argument("--out-name", help="time mode: folder name inside the run")
    args = parser.parse_args()
    run = Path(args.run).resolve()

    def prepared(run, name, flux, flux_err, skip_outliers=True, solution_path=None):
        return prepare(run, name, args.bins, skip_outliers)

    timing.prepare = prepared
    if args.mode == "figure":
        timing.draw_scan(run, Path(args.file).stem)
    elif args.mode == "adopt":
        timing.run_adopt(run, args.file, tuple(args.start_over_d))
    elif args.mode == "scan":
        timing.run_scan(run, args.file, "sap_flux", "sap_flux_err", True)
    else:
        if args.half_over_d is None:
            parser.error("time mode needs --half-over-d")
        timing.run_time(
            run,
            args.file,
            "sap_flux",
            "sap_flux_err",
            tuple(args.half_over_d),
            True,
            out_dir=None if args.out_name is None else run / args.out_name,
        )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
