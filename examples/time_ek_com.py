"""Time the eclipses of EK Com, sector 49, the way the regression test does.

The fit is the one in ``examples/fit_ek_com.py``. The timing recipe is the
one in ``tests/tess/test_regression_real.py``: half-width 0.25 cycles (phase
from -0.25 to +0.25, quadrature to quadrature), both sides of the predicted
centre required, and the contact-wall cut left off, because a 10 minute
cadence does not sample both walls. The recorded counts are 70 primary and
68 secondary timings.

``window`` on the recipe is that half-width. The points fitted for eclipse
number E are those with ``|phase - E| < window``.

Run from the repository root::

    python examples/time_ek_com.py

Writes ``examples/out/ek_com_timing/``: ``primary.dat``, ``secondary.dat``,
``fit.png``, ``trend.png``, ``oc.png``, ``windows_on_fold.png``, and one
figure per timed eclipse in ``eclipses/``.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import mikulasek_eb as eb
from mikulasek_eb.helpers import predicted_magnitude
from mikulasek_eb.model import phase
from mikulasek_eb.plotting import LINE_ADOPTED_COLOUR
from mikulasek_eb.timing import fit_ephemeris
from mikulasek_eb.tom import write_timing_ascii

from draw_fit import save_figures
from fit_ek_com import DATA, run

logger = logging.getLogger(__name__)
OUT = Path(__file__).resolve().parents[1] / "examples" / "out" / "ek_com_timing"

# Same choices as the regression test. The default recipe also demands points
# on both contact walls, which this cadence does not provide.
RECIPE = eb.TimingRecipe(require_contact_walls=False, clip_outliers=False)
PRIMARY_COLOUR = "#0072B2"
SECONDARY_COLOUR = "#D55E00"


def half_width_sentence(half: float) -> str:
    """State the fitted phase interval without calling it a full width.

    Args:
        half (float): ``TimingRecipe.window``, a half-width in cycles.

    Returns:
        str: One sentence naming the half-width and the inequality.
    """
    return (
        f"Half-width {half:.2f} cycles.\n"
        f"Points with |phase - E| < {half:.2f} are the ones fitted."
    )


def fitted_event(
    t_mjd: np.ndarray,
    mag: np.ndarray,
    sig: np.ndarray,
    model,
    p: dict,
    ephem: dict,
    t_jd: float,
    delta_days: float,
    half: float,
    order: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Rebuild one eclipse the way ``fit_minimum`` fitted it.

    Args:
        t_mjd (ndarray): Times used in the fit (MJD).
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        model: Template model.
        p (dict): Template parameters.
        ephem (dict): Fitted ephemeris, M0 in MJD.
        t_jd (float): Fitted time of minimum (absolute JD).
        delta_days (float): Shift from the predicted centre (days).
        half (float): Half-width in cycles.
        order (int): Per-event baseline degree.

    Returns:
        tuple: Phase offset from the predicted centre, magnitudes, the
        template plus the per-event baseline, and the shift in cycles.

    Raises:
        RuntimeError: If the half-width contains no points.
    """
    period = float(ephem["P0"])
    pdot = float(ephem.get("Pdot", 0.0))
    theta = phase(t_mjd, ephem["M0"], period, pdot)
    predicted_mjd = float(eb.absolute_jd_to_mjd(t_jd - delta_days))
    centre = float(phase(np.array([predicted_mjd]), ephem["M0"], period, pdot)[0])
    offset = theta - centre
    selected = np.abs(offset) < half
    if not np.any(selected):
        raise RuntimeError(
            f"no points with |phase - E| < {half} around MJD {predicted_mjd}"
        )
    times = t_mjd[selected]
    values = mag[selected]
    weights = 1.0 / sig[selected]
    local_period = period * (1.0 + pdot * centre)
    shift_cycles = delta_days / local_period
    template = np.asarray(model(theta[selected] - shift_cycles, p), float)
    scale = float(times.std())
    x = (times - times.mean()) / (scale if scale > 0.0 else 1.0)
    degree = min(order, len(times) - 1)
    basis = np.vander(x, degree + 1, increasing=True)
    coef = np.linalg.lstsq(
        basis * weights[:, None], weights * (values - template), rcond=None
    )[0]
    return offset[selected], values, template + basis @ coef, shift_cycles


def draw_windows_on_fold(
    theta: np.ndarray,
    detrended: np.ndarray,
    template: np.ndarray,
    half: float,
    path: Path,
    title: str,
) -> None:
    """Shade the fitted half-width of each eclipse on the folded curve.

    Args:
        theta (ndarray): Fitted phase (cycles).
        detrended (ndarray): Magnitudes with the global baseline removed.
        template (ndarray): Template at ``theta``, without the baseline.
        half (float): Half-width in cycles.
        path (Path): PNG file to write.
        title (str): Figure title.
    """
    phi = theta - np.round(theta)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.axvspan(-half, half, color=PRIMARY_COLOUR, alpha=0.15, label="primary half-width")
    ax.axvspan(-0.5, -0.5 + half, color=SECONDARY_COLOUR, alpha=0.15, label="secondary half-width")
    ax.axvspan(0.5 - half, 0.5, color=SECONDARY_COLOUR, alpha=0.15)
    ax.plot(phi, detrended, ".", ms=2, color="0.45")
    order = np.argsort(phi)
    ax.plot(phi[order], template[order], color=LINE_ADOPTED_COLOUR, lw=1.2, label="template")
    ax.axvline(0.0, color=PRIMARY_COLOUR, ls="--", lw=0.8)
    ax.axvline(0.5, color=SECONDARY_COLOUR, ls="--", lw=0.8)
    ax.axvline(-0.5, color=SECONDARY_COLOUR, ls="--", lw=0.8)
    ax.set_xlim(-0.5, 0.5)
    ax.invert_yaxis()
    ax.set_xlabel("Phase from the primary (cycles)")
    ax.set_ylabel("Magnitude, baseline removed")
    ax.set_title(title)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)


def draw_eclipse(
    offset: np.ndarray,
    mag: np.ndarray,
    curve: np.ndarray,
    shift_cycles: float,
    half: float,
    path: Path,
    title: str,
) -> None:
    """Draw one timed eclipse and the shift of its template.

    Args:
        offset (ndarray): Phase minus the predicted centre (cycles).
        mag (ndarray): Magnitudes in the half-width.
        curve (ndarray): Template plus the per-event baseline.
        shift_cycles (float): Fitted shift (cycles).
        half (float): Half-width in cycles.
        path (Path): PNG file to write.
        title (str): Figure title.
    """
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(offset, mag, ".", ms=4, color="0.3", label="data")
    order = np.argsort(offset)
    ax.plot(offset[order], curve[order], color=LINE_ADOPTED_COLOUR, lw=1.4, label="template + baseline")
    ax.axvline(0.0, color="0.2", ls="--", lw=0.8, label="predicted centre")
    ax.axvline(shift_cycles, color=LINE_ADOPTED_COLOUR, ls=":", lw=1.0, label="fitted shift")
    ax.set_xlim(-half, half)
    ax.invert_yaxis()
    ax.set_xlabel("Phase from the predicted centre (cycles)")
    ax.set_ylabel("Magnitude")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def eclipse_path(folder: Path, kind: str, eclipse: float) -> Path:
    """Name one eclipse figure from its kind and eclipse number.

    Args:
        folder (Path): Directory for the eclipse figures.
        kind (str): ``primary`` or ``secondary``.
        eclipse (float): Eclipse number. Secondaries are half-integers.

    Returns:
        Path: ``primary_E12.png`` or ``secondary_E12.5.png``.
    """
    if abs(eclipse - round(eclipse)) < 1e-6:
        name = f"{kind}_E{int(round(eclipse))}.png"
    else:
        name = f"{kind}_E{eclipse:.1f}.png"
    return folder / name


def draw_oc(timings: dict, path: Path, title: str) -> None:
    """Draw the O-C of the primary and the secondary about a linear ephemeris.

    Args:
        timings (dict): ``primary`` and ``secondary`` blocks from ``time``.
        path (Path): PNG file to write.
        title (str): Figure title.
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    for kind, colour in (("primary", "#0072B2"), ("secondary", "#D55E00")):
        block = timings[kind]
        if len(block["T"]) == 0:
            continue
        _coef, residual, _scale = fit_ephemeris(
            block["T"], block["sigma"], block["E"], 1
        )
        ax.errorbar(
            block["E"],
            residual * 86400.0,
            yerr=np.asarray(block["sigma"]) * 86400.0,
            fmt="o",
            ms=3,
            color=colour,
            label=kind,
        )
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.set_xlabel("Eclipse number E")
    ax.set_ylabel("O-C (seconds)")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def main() -> None:
    """Fit EK Com, time both eclipses, and write the files and the figures."""
    payload, out = run()
    pipeline = out.pipeline
    adopted = pipeline.adopted
    timed = eb.time(out.t_mjd, out.mag, out.sig, out.solution, RECIPE)
    OUT.mkdir(parents=True, exist_ok=True)
    title = f"{payload['title']}  P={payload['period']:.6f} d"
    save_figures(
        OUT,
        title,
        t_mjd=pipeline.t_mjd,
        mag=pipeline.mag,
        sig=pipeline.sig,
        model=pipeline.model,
        p_single=pipeline.single_component.p,
        p_adopted=adopted.p,
        ephem=adopted.ephem,
        seg=pipeline.seg,
        order=pipeline.order,
        timescale=payload["timescale"],
    )
    yhat, template, theta = predicted_magnitude(
        pipeline.t_mjd,
        pipeline.mag,
        pipeline.sig,
        pipeline.model,
        adopted.p,
        adopted.ephem,
        pipeline.seg,
        pipeline.order,
    )
    fold = OUT / "windows_on_fold.png"
    sentence = half_width_sentence(RECIPE.window)
    draw_windows_on_fold(
        theta,
        pipeline.mag - (yhat - template),
        template,
        RECIPE.window,
        fold,
        f"{payload['title']}\n{sentence}",
    )
    logger.info("wrote %s", fold)
    eclipses = OUT / "eclipses"
    eclipses.mkdir(parents=True, exist_ok=True)
    model = out.solution.build_model()
    ephem = out.solution.ephem
    n_figures = 0
    for kind in ("primary", "secondary"):
        block = timed.timings[kind]
        for index, eclipse in enumerate(block["E"]):
            offset, values, curve, shift = fitted_event(
                out.t_mjd,
                out.mag,
                out.sig,
                model,
                out.solution.p,
                ephem,
                float(block["T"][index]),
                float(block["delta"][index]),
                RECIPE.window,
                RECIPE.baseline_order,
            )
            if len(offset) != int(block["n"][index]):
                raise RuntimeError(
                    f"{kind} E={eclipse}: the figure has {len(offset)} points, "
                    f"the timing used {int(block['n'][index])}"
                )
            path = eclipse_path(eclipses, kind, float(eclipse))
            draw_eclipse(
                offset,
                values,
                curve,
                shift,
                RECIPE.window,
                path,
                f"{kind} E={float(eclipse):g}\n{sentence}",
            )
            n_figures += 1
    logger.info("wrote %d eclipse figures in %s", n_figures, eclipses)
    for kind in ("primary", "secondary"):
        block = timed.timings[kind]
        path = OUT / f"{kind}.dat"
        write_timing_ascii(
            path,
            block,
            lightcurve=DATA,
            period=payload["period"],
            t0_jd=payload["epoch_jd"],
            method="mikulasek_eb.time",
            window=RECIPE.window,
            events=kind,
            extra_comments=[
                "require_contact_walls: false",
                "clip_outliers: false",
            ],
            timescale=payload["timescale"],
            rejected=timed.rejected[kind],
            recipe=RECIPE,
        )
        sigma_s = (
            86400.0 * float(np.median(block["sigma"]))
            if len(block["T"])
            else float("nan")
        )
        logger.info(
            "%s: %d timings, %d rejected, median sigma %.1f s",
            kind,
            len(block["T"]),
            len(timed.rejected[kind]),
            sigma_s,
        )
        logger.info("wrote %s", path)
    figure = OUT / "oc.png"
    draw_oc(timed.timings, figure, f"{payload['title']}\nO-C about a linear ephemeris")
    logger.info("wrote %s", figure)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
