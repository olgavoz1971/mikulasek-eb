"""Borrow the sector 30 (120 s) template for sector 4 (1800 s) and time sector 4.

The reference is the chosen solution of the tiered fit of sector 30 at 120 s
(``results/nsv807/s30_protocol/chosen``). Sector 4 keeps its own data layer: the
baseline order and gap split of its own tiered fit (``results/nsv807/s4_protocol``)
and the exposure of 1440 s. Profile, Gamma, C, K and the period are frozen in every
variant, and the epoch is fitted. The variants differ by what else is fitted:

* ``borrow_phase``: the epoch only.
* ``borrow_amp``: epoch and depths.
* ``borrow_width``: epoch and widths D.
* ``borrow_amp_width``: epoch, depths and widths.
* ``borrow_amp_width_harm``: as above, with the harmonics and the O'Connell term
  free (the 2018 maxima differ from 2020, an earlier trial needed them free).

Every variant is timed with the window timing of ``nsv807_timing.py`` at 3.0 D and
3.2 D. Nothing is selected here: all variants are reported.

Run from ``examples``::

    python nsv807_borrow.py

Writes ``results/nsv807/s4_protocol/borrow/<variant>/``: the solution, ``fits/`` and
``trends/`` figures, and ``timing/`` (``eclipses.csv``, ``.dat`` files, figures),
and ``borrow/fits.csv``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

import mikulasek_eb as eb
import nsv807_s30_tiered as tiered
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.plotting import plot_global_fit, plot_instrumental_trend
from nsv807_timing import run_time
from window_scan_mik import write_csv

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "nsv807"
REFERENCE_RUN = RESULTS / "s30_protocol"
TARGET_RUN = RESULTS / "s4_protocol"
TARGET_FILE = "TIC35119266_s4_TESS-SPOC_1800s.vot"
VARIANTS = {
    "borrow_phase": ("zero_phase",),
    "borrow_amp": ("zero_phase", "amplitudes"),
    "borrow_width": ("zero_phase", "widths"),
    "borrow_amp_width": ("zero_phase", "amplitudes", "widths"),
    "borrow_amp_width_harm": ("zero_phase", "amplitudes", "widths", "harmonics"),
}
HALF_OVER_D = (3.0, 3.2)


def single_solution(folder: Path) -> eb.ModelSolution:
    """Load the one solution of a ``chosen`` folder.

    Args:
        folder (Path): ``chosen`` folder of a tiered-fit run.

    Returns:
        ModelSolution: The stored solution.

    Raises:
        FileNotFoundError: If the folder does not hold exactly one solution.
    """
    found = sorted(folder.glob("*.json"))
    if len(found) != 1:
        raise FileNotFoundError(f"expected one solution in {folder}: {found}")
    return eb.ModelSolution.load(found[0])


def save_figures(folder: Path, stem: str, out, reference_p: dict, title: str) -> None:
    """Draw the profile figure and the trend figure of one borrowed fit.

    Args:
        folder (Path): Variant folder.
        stem (str): File stem.
        out (FitOutput): The borrowed fit.
        reference_p (dict): The sector 30 template, drawn in grey.
        title (str): Figure title.
    """
    solution = out.solution
    data = solution.config.data
    seg = (
        None
        if data.gap_split_days is None
        else np.cumsum(np.r_[0, np.diff(out.t_mjd) > data.gap_split_days])
    )
    model = solution.build_model()
    yhat, template, theta = predicted_magnitude(
        out.t_mjd,
        out.mag,
        out.sig,
        model,
        solution.p,
        solution.ephem,
        seg,
        data.baseline_order,
    )
    (folder / "fits").mkdir(parents=True, exist_ok=True)
    (folder / "trends").mkdir(parents=True, exist_ok=True)
    plot_global_fit(
        eb.mjd_to_jd(out.t_mjd),
        out.mag,
        out.sig,
        yhat,
        template,
        theta,
        model,
        reference_p,
        solution.p,
        seg,
        title,
        save_path=str(folder / "fits" / f"{stem}.png"),
        note="",
    )
    plot_instrumental_trend(
        eb.mjd_to_jd(out.t_mjd),
        out.mag,
        yhat,
        template,
        folder / "trends" / f"{stem}.png",
        title,
    )


def main() -> None:
    """Borrow, fit and time sector 4 under every variant."""
    reference = single_solution(REFERENCE_RUN / "chosen")
    stored = single_solution(TARGET_RUN / "chosen")
    tiered.FLUX, tiered.FLUX_ERR = "sap_flux", "sap_flux_err"
    tiered.OUT = TARGET_RUN
    payload = tiered.load_sap(tiered.DATA / TARGET_FILE)
    stem = Path(TARGET_FILE).stem
    t_mjd = eb.absolute_jd_to_mjd(payload["t"])
    rows = []
    for name, groups in VARIANTS.items():
        folder = TARGET_RUN / "borrow" / name
        folder.mkdir(parents=True, exist_ok=True)
        out = eb.fit_borrowed(
            t_mjd,
            payload["mag"],
            payload["sig"],
            reference,
            stored.config.data,
            eb.ParameterPlan(groups),
            data_id=stem,
        )
        stats = out.solution.stats
        rms = 1e3 * float(stats["rms_mag"])
        path = folder / f"{stem}__{name}.json"
        out.solution.save(path)
        title = tiered.figure_title(stem, name, stats["status"], False, rms)
        save_figures(folder, stem, out, reference.p, title)
        rows.append(
            {
                "variant": name,
                "free": " ".join(groups),
                "status": stats["status"],
                "rms_mmag": f"{rms:.2f}",
                "n_points": stats["n_points"],
            }
        )
        logger.info("%s status %s rms %.2f mmag", name, stats["status"], rms)
        if stats["status"] != "ok":
            logger.warning("%s is not status ok, so it is not timed", name)
            continue
        run_time(
            TARGET_RUN,
            TARGET_FILE,
            "sap_flux",
            "sap_flux_err",
            HALF_OVER_D,
            skip_outliers=True,
            solution_path=path,
            out_dir=folder / "timing",
        )
    write_csv(TARGET_RUN / "borrow" / "fits.csv", rows)
    logger.info("wrote %s", TARGET_RUN / "borrow")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
