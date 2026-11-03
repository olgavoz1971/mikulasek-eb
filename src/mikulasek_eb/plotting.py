"""Diagnostic figures of the global fit (requires matplotlib)."""

from __future__ import annotations


import numpy as np
import matplotlib.pyplot as plt

from .model import EBModel

# Primary zoom panel: |phase| < this (cycles).
PRIMARY_PHASE_ZOOM_HALF = 0.12


# Plot styling: single-component vs adopted template.
LINE_SINGLE_COLOUR = "#5b6b7a"


LINE_SINGLE_WIDTH = 1.5


LINE_ADOPTED_COLOUR = "#ff2e2e"


# Fold/zoom overlay can stay a bit thicker than the time-series line.
LINE_ADOPTED_WIDTH = 2.0


# Time-series fit: thin enough not to hide 10-min / 2-min points.
LINE_TIME_SERIES_WIDTH = 0.8


# Zoomed minima: primary points stay blue, secondary points are orange.
POINT_PRIMARY_COLOUR = "C0"
POINT_SECONDARY_COLOUR = "C1"


# Full folded panel (not the primary zoom): two cycles so Min II at ±0.5 is interior.
FOLD_PHASE_MIN = -1.0


FOLD_PHASE_MAX = 1.0


FOLD_PHASE_GRID_N = 4000


def replicate_folded_points(
    phi: np.ndarray,
    *cols: np.ndarray,
    ph_min: float = FOLD_PHASE_MIN,
    ph_max: float = FOLD_PHASE_MAX,
) -> tuple[np.ndarray, ...]:
    """Repeat one orbital cycle of folded points until they fill ``[ph_min, ph_max]``.

    ``phi`` is expected in ``[-0.5, 0.5]``. Integer cycle shifts are appended
    so a secondary minimum that sits on the edge of one cycle appears in the
    interior of a two-cycle window.

    Args:
        phi (ndarray): Folded phase in one cycle.
        *cols (ndarray): Parallel arrays (magnitude, sigma, residuals, ...).
        ph_min (float): Left edge of the display window (cycles).
        ph_max (float): Right edge of the display window (cycles).

    Returns:
        tuple[ndarray, ...]: ``(phi_out, *cols_out)`` restricted to the window.
    """
    k_lo = int(np.floor(ph_min - 0.5))
    k_hi = int(np.ceil(ph_max + 0.5))
    phi_parts = []
    col_parts = [[] for _ in cols]
    for k in range(k_lo, k_hi + 1):
        shifted = phi + k
        keep = (shifted >= ph_min) & (shifted <= ph_max)
        if not np.any(keep):
            continue
        phi_parts.append(shifted[keep])
        for i, col in enumerate(cols):
            col_parts[i].append(col[keep])
    phi_out = np.concatenate(phi_parts)
    extras = tuple(np.concatenate(p) for p in col_parts)
    return (phi_out, *extras)


def _secondary_offset(phi: np.ndarray) -> np.ndarray:
    """Phase distance from the secondary minimum.

    ``phi`` is folded into ``[-0.5, 0.5]``. The secondary sits at both edges.

    Args:
        phi (ndarray): Folded phase (cycles).

    Returns:
        ndarray: Offset from phase 0.5, in ``[-0.5, 0.5]``.
    """
    return np.where(phi >= 0.0, phi - 0.5, phi + 0.5)


def _draw_template(ax, x, y, *, colour, width, label):
    """Draw one template curve.

    Args:
        ax: Matplotlib axes.
        x (ndarray): Horizontal coordinate.
        y (ndarray): Template magnitude.
        colour (str): Line colour. The same function keeps this colour in every panel.
        width (float): Line width.
        label (str | None): Legend entry. None draws the line without a new name.
    """
    ax.plot(x, y, color=colour, lw=width, zorder=5 if colour == LINE_ADOPTED_COLOUR else 4, label=label)


def plot_global_fit(
    t_jd: np.ndarray,
    mag: np.ndarray,
    sig: np.ndarray,
    yhat: np.ndarray,
    f_template: np.ndarray,
    theta: np.ndarray,
    model: EBModel,
    p_single: dict,
    p_adopted: dict,
    seg: np.ndarray | None,
    title: str,
    timescale: str | None = None,
    save_path: str | None = None,
    note: str | None = None,
) -> None:
    """Plot time series, fold, both minima, and residuals (matplotlib).

    The adopted template is the red line in every panel. The single-component
    template is the grey line, and only when its parameters differ. One
    component therefore has no grey line. The zoom panel shows the primary and
    the secondary against phase measured from each minimum's own centre, with
    those same two lines.

    Args:
        t_jd (ndarray): Absolute Julian Date for the time axis.
        mag (ndarray): Observed magnitudes.
        sig (ndarray): Magnitude uncertainties.
        yhat (ndarray): Full prediction from the adopted fit.
        f_template (ndarray): Adopted template without baselines.
        theta (ndarray): Fitted phase (cycles).
        model (EBModel): Fitted model.
        p_single (dict): Single-component template parameters.
        p_adopted (dict): Adopted template parameters.
        seg (ndarray | None): Baseline segment indices.
        title (str): Figure title.
        timescale (str, optional): TIMESYS timescale for the axis label.
        save_path (str, optional): If set, also save the figure to this path.
        note (str, optional): Figure subtitle; the default describes single- or
            two-component templates.
    """
    phi = theta - np.round(theta)
    phi_grid = np.linspace(FOLD_PHASE_MIN, FOLD_PHASE_MAX, FOLD_PHASE_GRID_N)
    zoom = np.linspace(-PRIMARY_PHASE_ZOOM_HALF, PRIMARY_PHASE_ZOOM_HALF, 2000)
    mag_corr = mag - (yhat - f_template)
    same = p_single is p_adopted or p_single == p_adopted
    # Default matplotlib size is 10 pt. Twice that, for this figure only.
    style = {
        "font.size": 20,
        "axes.labelsize": 20,
        "axes.titlesize": 22,
        "legend.fontsize": 18,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
    }

    with plt.rc_context(style):
        _draw_global_fit(
            t_jd, mag, sig, yhat, phi, phi_grid, zoom, mag_corr,
            same, model, p_single, p_adopted, seg, title, timescale, save_path, note,
        )


def _draw_global_fit(
    t_jd, mag, sig, yhat, phi, phi_grid, zoom, mag_corr,
    same, model, p_single, p_adopted, seg, title, timescale, save_path, note,
) -> None:
    """Draw the four panels. Called with the enlarged font already set.

    Args:
        t_jd (ndarray): Absolute Julian Date.
        mag (ndarray): Observed magnitudes.
        sig (ndarray): Magnitude uncertainties.
        yhat (ndarray): Full adopted prediction.
        phi (ndarray): Phase folded into one cycle.
        phi_grid (ndarray): Fine phase grid for the folded panel.
        zoom (ndarray): Phase offset grid for a minimum, centred on zero.
        mag_corr (ndarray): Magnitudes with the baseline removed.
        same (bool): True when the single-component and adopted parameters match.
        model (EBModel): Fitted model.
        p_single (dict): Single-component parameters.
        p_adopted (dict): Adopted parameters.
        seg (ndarray | None): Baseline segment indices.
        title (str): Figure title.
        timescale (str | None): Time-scale name for the axis.
        save_path (str | None): File to write.
        note (str | None): Figure subtitle.
    """
    fig, axes = plt.subplots(4, 1, figsize=(16, 24), sharex=False)
    ax = axes[0]
    ax.errorbar(
        t_jd,
        mag,
        yerr=sig,
        fmt="o",
        ms=2,
        alpha=0.35,
        elinewidth=0.4,
        capsize=0,
        color="C0",
        zorder=1,
        label="data",
    )
    if seg is None:
        ax.plot(
            t_jd,
            yhat,
            color=LINE_ADOPTED_COLOUR,
            lw=LINE_TIME_SERIES_WIDTH,
            zorder=5,
            label="fit",
        )
    else:
        for i, s in enumerate(np.unique(seg)):
            m = seg == s
            ax.plot(
                t_jd[m],
                yhat[m],
                color=LINE_ADOPTED_COLOUR,
                lw=LINE_TIME_SERIES_WIDTH,
                zorder=5,
                label="fit" if i == 0 else None,
            )
    axis = "Julian Date (d)" if timescale is None else f"Julian Date ({timescale}, d)"
    ax.set_xlabel(axis)
    ax.set_ylabel("Magnitude")
    ax.set_title(title)
    ax.invert_yaxis()
    ax.legend()

    ax = axes[1]
    phi_fold, mag_fold, sig_fold = replicate_folded_points(phi, mag_corr, sig)
    ord_ = np.argsort(phi_fold)
    ax.errorbar(
        phi_fold[ord_],
        mag_fold[ord_],
        yerr=sig_fold[ord_],
        fmt="o",
        ms=2,
        alpha=0.35,
        elinewidth=0.4,
        capsize=0,
        color="C0",
        zorder=1,
        label="baseline-corrected",
    )
    if not same:
        ax.plot(
            phi_grid,
            model(phi_grid, p_single),
            color=LINE_SINGLE_COLOUR,
            lw=LINE_SINGLE_WIDTH,
            zorder=4,
            label="single-component",
        )
    ax.plot(
        phi_grid,
        model(phi_grid, p_adopted),
        color=LINE_ADOPTED_COLOUR,
        lw=LINE_ADOPTED_WIDTH,
        zorder=5,
        label="adopted fit",
    )
    ax.set_xlabel("Phase")
    ax.set_xlim(FOLD_PHASE_MIN, FOLD_PHASE_MAX)
    ax.set_ylabel("Magnitude (template scale)")
    ax.invert_yaxis()
    ax.legend()

    ax = axes[2]
    in_p = np.abs(phi) < PRIMARY_PHASE_ZOOM_HALF
    sec_off = _secondary_offset(phi)
    in_s = np.abs(sec_off) < PRIMARY_PHASE_ZOOM_HALF
    ax.errorbar(
        phi[in_p],
        mag_corr[in_p],
        yerr=sig[in_p],
        fmt="o",
        ms=4,
        alpha=0.5,
        elinewidth=0.4,
        capsize=0,
        color=POINT_PRIMARY_COLOUR,
        zorder=1,
        label="primary",
    )
    ax.errorbar(
        sec_off[in_s],
        mag_corr[in_s],
        yerr=sig[in_s],
        fmt="s",
        ms=4,
        alpha=0.5,
        elinewidth=0.4,
        capsize=0,
        color=POINT_SECONDARY_COLOUR,
        zorder=1,
        label="secondary",
    )
    if not same:
        _draw_template(
            ax, zoom, model(zoom, p_single),
            colour=LINE_SINGLE_COLOUR, width=LINE_SINGLE_WIDTH, label="single-component",
        )
        _draw_template(
            ax, zoom, model(0.5 + zoom, p_single),
            colour=LINE_SINGLE_COLOUR, width=LINE_SINGLE_WIDTH, label=None,
        )
    _draw_template(
        ax, zoom, model(zoom, p_adopted),
        colour=LINE_ADOPTED_COLOUR, width=LINE_ADOPTED_WIDTH, label="adopted fit",
    )
    _draw_template(
        ax, zoom, model(0.5 + zoom, p_adopted),
        colour=LINE_ADOPTED_COLOUR, width=LINE_ADOPTED_WIDTH, label=None,
    )
    ax.set_xlabel("Phase from minimum centre")
    ax.set_ylabel("Magnitude (template scale)")
    ax.invert_yaxis()
    ax.legend()

    ax = axes[3]
    res = mag - yhat
    phi_res, res_fold = replicate_folded_points(phi, res)
    ax.scatter(phi_res, 1e3 * res_fold, s=4, alpha=0.35, zorder=1, label="residual")
    ax.axhline(0.0, color=LINE_ADOPTED_COLOUR, lw=1.2, zorder=5)
    ax.set_xlabel("Phase")
    ax.set_xlim(FOLD_PHASE_MIN, FOLD_PHASE_MAX)
    ax.set_ylabel("Residual (mmag)")
    ax.legend()
    if note is None:
        note = (
            "Two-component primary: extension, not Mikulasek (2015)."
            if not same
            else "Single-component template."
        )
    fig.suptitle(note, y=0.995)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=120, bbox_inches="tight")
    import matplotlib

    if matplotlib.get_backend().lower() != "agg":
        plt.show()
    else:
        plt.close(fig)


def plot_instrumental_trend(
    t,
    mag,
    yhat,
    template,
    save_path,
    title: str,
    *,
    injected_mag=None,
) -> None:
    """Save the fitted instrumental trend next to the light curve and the template.

    The trend is the per-segment baseline of the fit: the full prediction minus the
    eclipse template. When ``injected_mag`` is given it is drawn on the same panel,
    shifted to the median of the fitted trend, so the shapes can be compared; a
    constant between them is the magnitude zero point.

    Args:
        t (ndarray): Times (days).
        mag (ndarray): Observed magnitudes.
        yhat (ndarray): Full prediction, template plus baseline.
        template (ndarray): Eclipse template at the same times, without the baseline.
        save_path (str | Path): File to write.
        title (str): Figure title.
        injected_mag (ndarray | None): Instrumental trend known independently, in magnitudes.
    """
    t = np.asarray(t, float)
    mag = np.asarray(mag, float)
    trend = np.asarray(yhat, float) - np.asarray(template, float)
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    axes[0].plot(t, mag, ".", ms=2, color="0.6", label="data")
    axes[0].plot(t, yhat, "r-", lw=0.8, label="template + trend")
    axes[0].invert_yaxis()
    axes[0].set_ylabel("magnitude")
    axes[0].legend(loc="upper right")
    axes[0].set_title(title)
    axes[1].plot(t, trend, "b-", lw=0.8, label="fitted trend")
    if injected_mag is not None:
        injected = np.asarray(injected_mag, float)
        injected = injected - np.median(injected) + np.median(trend)
        axes[1].plot(t, injected, "k-", lw=0.8, alpha=0.7, label="injected trend, median aligned")
    axes[1].set_ylabel("trend (mag)")
    axes[1].legend(loc="upper right")
    axes[2].plot(t, mag - trend, ".", ms=2, color="0.6", label="data minus trend")
    axes[2].plot(t, template, "r-", lw=0.8, label="template")
    axes[2].invert_yaxis()
    axes[2].set_ylabel("magnitude")
    axes[2].set_xlabel("time (days)")
    axes[2].legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(save_path, dpi=100)
    plt.close(fig)
