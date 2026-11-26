"""Diagnostic figures of the global fit (requires matplotlib)."""

from __future__ import annotations


from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from .helpers import absolute_jd_to_mjd
from .model import EBModel, phase, predicted_time
from .timing import effective_window, event_curve

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
    ax.plot(
        x,
        y,
        color=colour,
        lw=width,
        zorder=5 if colour == LINE_ADOPTED_COLOUR else 4,
        label=label,
    )


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
            t_jd,
            mag,
            sig,
            yhat,
            phi,
            phi_grid,
            zoom,
            mag_corr,
            same,
            model,
            p_single,
            p_adopted,
            seg,
            title,
            timescale,
            save_path,
            note,
        )


def _draw_global_fit(
    t_jd,
    mag,
    sig,
    yhat,
    phi,
    phi_grid,
    zoom,
    mag_corr,
    same,
    model,
    p_single,
    p_adopted,
    seg,
    title,
    timescale,
    save_path,
    note,
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
            ax,
            zoom,
            model(zoom, p_single),
            colour=LINE_SINGLE_COLOUR,
            width=LINE_SINGLE_WIDTH,
            label="single-component",
        )
        _draw_template(
            ax,
            zoom,
            model(0.5 + zoom, p_single),
            colour=LINE_SINGLE_COLOUR,
            width=LINE_SINGLE_WIDTH,
            label=None,
        )
    _draw_template(
        ax,
        zoom,
        model(zoom, p_adopted),
        colour=LINE_ADOPTED_COLOUR,
        width=LINE_ADOPTED_WIDTH,
        label="adopted fit",
    )
    _draw_template(
        ax,
        zoom,
        model(0.5 + zoom, p_adopted),
        colour=LINE_ADOPTED_COLOUR,
        width=LINE_ADOPTED_WIDTH,
        label=None,
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
        axes[1].plot(
            t, injected, "k-", lw=0.8, alpha=0.7, label="injected trend, median aligned"
        )
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


# ------------------------------------------------------------ timing windows
DAY_S = 86400.0
WINDOW_REJECTED_COLOUR = "0.55"


def eclipse_picture_stem(kind: str, E: float) -> str:
    """File stem of one eclipse picture, in the MAVKA layout ``E00_secondary``.

    Eclipse number ``E`` is counted from T0: an integer is a primary, ``n + 0.5``
    is the secondary of the same cycle ``n``. The stem uses the cycle number, so
    ``E=1`` is ``E01_primary`` and ``E=0.5`` is ``E00_secondary``.

    Args:
        kind (str): ``primary`` or ``secondary``.
        E (float): Eclipse number from T0.

    Returns:
        str: Stem without the ``.png`` suffix.
    """
    cycle = int(np.floor(float(E) + 1e-9))
    return f"E{cycle:02d}_{kind}"


def _plot_shifted_template(ax, t_rel, y_grid, in_grid) -> None:
    """Draw the template in red inside the window and in grey outside it.

    Args:
        ax (matplotlib.axes.Axes): Target axes.
        t_rel (ndarray): Times from the predicted minimum (days).
        y_grid (ndarray): Template plus baseline.
        in_grid (ndarray): True on the fitted window.
    """
    inside = np.where(in_grid, y_grid, np.nan)
    outside = np.where(~in_grid, y_grid, np.nan)
    ax.plot(
        t_rel,
        outside,
        color=LINE_SINGLE_COLOUR,
        lw=LINE_SINGLE_WIDTH,
        zorder=4,
        label="template outside the window",
    )
    ax.plot(
        t_rel,
        inside,
        color=LINE_ADOPTED_COLOUR,
        lw=LINE_ADOPTED_WIDTH,
        zorder=5,
        label="template in the window",
    )


def _timing_events(t_mjd, solution, timing, kind: str) -> list[dict]:
    """List the timed and the rejected eclipses of one kind, numbered from M0 and from T0.

    Args:
        t_mjd (ndarray): Times (MJD) that were timed.
        solution (ModelSolution): Solution whose template was used.
        timing (TimingOutput): Result of ``mikulasek_eb.time``.
        kind (str): ``primary`` or ``secondary``.

    Returns:
        list[dict]: One dict per eclipse, sorted by ``E``: ``E`` (from T0), ``E_m0``
        (engine number), ``timed`` (bool), ``T_mjd``, ``sigma`` and ``delta`` (days),
        ``n``, ``reason``, ``detail``.

    Raises:
        ValueError: If ``timing`` was made with another solution or has no ``kind``.
    """
    if timing.solution_hash != solution.solution_hash():
        raise ValueError("timing was made with a different solution")
    if kind not in timing.timings:
        raise ValueError(f"kind {kind!r} was not timed; timed: {list(timing.timings)}")
    ephem = solution.ephem
    events = []
    timed = timing.timings[kind]
    for i, e in enumerate(np.asarray(timed.get("E", []), float)):
        t_mjd_found = float(absolute_jd_to_mjd(float(timed["T"][i])))
        th = float(phase(t_mjd_found, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0)))
        events.append(
            dict(
                E=float(e),
                E_m0=float(np.round(2 * th) / 2),
                timed=True,
                T_mjd=t_mjd_found,
                sigma=float(timed["sigma"][i]),
                delta=float(timed["delta"][i]),
                n=int(timed["n"][i]),
                reason="",
                detail="",
            )
        )
    for item in timing.rejected[kind]:
        t_pred = float(absolute_jd_to_mjd(float(item["T_pred_jd"])))
        th = float(phase(t_pred, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0)))
        events.append(
            dict(
                E=float(item["E"]),
                E_m0=float(np.round(2 * th) / 2),
                timed=False,
                reason=item["reason"],
                detail=item["detail"],
            )
        )
    return sorted(events, key=lambda r: r["E"])


def _eclipse_geometry(t, y, sig, model, solution, timing, ev, kind) -> dict:
    """Arrays and labels needed to draw one eclipse, in MJD of the observations.

    Args:
        t (ndarray): Times (MJD).
        y (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        model (EBModel): Template model.
        solution (ModelSolution): Solution whose template was used.
        timing (TimingOutput): Result of ``mikulasek_eb.time``.
        ev (dict): One element of ``_timing_events``.
        kind (str): ``primary`` or ``secondary``.

    Returns:
        dict: Times in MJD, masks, the optional template curve, window edges and title.
    """
    ephem = solution.ephem
    window = effective_window(
        timing.recipe.window, timing.recipe.allow_overlapping_windows
    )
    colour = POINT_PRIMARY_COLOUR if kind == "primary" else POINT_SECONDARY_COLOUR
    E = ev["E_m0"]
    th = phase(t, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0))
    t_pred = float(predicted_time(E, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0)))
    cycle_label = int(np.floor(ev["E"] + 1e-9))
    delta = float(ev["delta"]) if ev["timed"] else 0.0
    try:
        curve = event_curve(
            t,
            y,
            sig,
            model,
            solution.p,
            ephem,
            E,
            delta,
            window,
            timing.recipe.baseline_order,
        )
        edge = curve["half_days"]
        t_grid = curve["t_grid"]
        y_grid = curve["y_grid"]
        in_grid = curve["in_grid"]
    except ValueError:
        t_grid = y_grid = in_grid = None
        Pl = ephem["P0"] * (1.0 + ephem.get("Pdot", 0.0) * E)
        edge = window * Pl
    if ev["timed"]:
        t_found = float(ev["T_mjd"])
        title = (
            f"{kind}, cycle {cycle_label:d}, half-window {window:.3f}, n={ev['n']}\n"
            f"T = {t_found:.6f} MJD, sigma {ev['sigma'] * DAY_S:.1f} s"
        )
    else:
        t_found = None
        title = (
            f"{kind}, cycle {cycle_label:d}, half-window {window:.3f}\n"
            f"NOT TIMED: {ev['reason']} ({ev['detail']})"
        )
    return dict(
        t=t,
        y=y,
        cycle=np.abs(th - E) < 0.5,
        inside=np.abs(th - E) < window,
        t_grid=t_grid,
        y_grid=y_grid,
        in_grid=in_grid,
        t_lo=t_pred - edge,
        t_hi=t_pred + edge,
        colour=colour,
        timed=ev["timed"],
        t_found=t_found,
        sigma=ev.get("sigma"),
        title=title,
    )


def _mark_found_time(ax, geo) -> None:
    """Draw a vertical line at the found time of minimum, in MJD.

    Args:
        ax (matplotlib.axes.Axes): Target axes.
        geo (dict): Result of ``_eclipse_geometry``.
    """
    if geo["t_found"] is None:
        return
    ax.axvline(geo["t_found"], color=LINE_ADOPTED_COLOUR, lw=1.2, zorder=4)
    ax.axvspan(
        geo["t_found"] - geo["sigma"],
        geo["t_found"] + geo["sigma"],
        color=LINE_ADOPTED_COLOUR,
        alpha=0.25,
        zorder=1,
    )


def _scatter_window_points(ax, geo, size: float) -> None:
    """Draw the observed points that fall in the timing window.

    Args:
        ax (matplotlib.axes.Axes): Target axes.
        geo (dict): Result of ``_eclipse_geometry``.
        size (float): Marker size.
    """
    if geo["timed"]:
        ax.scatter(
            geo["t"][geo["inside"]],
            geo["y"][geo["inside"]],
            s=size,
            color=geo["colour"],
            zorder=3,
        )
    else:
        ax.scatter(
            geo["t"][geo["inside"]],
            geo["y"][geo["inside"]],
            s=size,
            facecolor="none",
            edgecolor=WINDOW_REJECTED_COLOUR,
            zorder=3,
        )


def _draw_cycle(ax, geo) -> None:
    """Draw the whole orbital cycle around one eclipse, in MJD.

    Args:
        ax (matplotlib.axes.Axes): Target axes.
        geo (dict): Result of ``_eclipse_geometry``.
    """
    t, y, inside, cycle = geo["t"], geo["y"], geo["inside"], geo["cycle"]
    span = geo["colour"] if geo["timed"] else WINDOW_REJECTED_COLOUR
    ax.axvspan(geo["t_lo"], geo["t_hi"], color=span, alpha=0.12, zorder=0)
    ax.scatter(t[cycle & ~inside], y[cycle & ~inside], s=5, color="0.65", zorder=2)
    if geo["t_grid"] is not None:
        _plot_shifted_template(ax, geo["t_grid"], geo["y_grid"], geo["in_grid"])
    _scatter_window_points(ax, geo, 12)
    _mark_found_time(ax, geo)
    ax.set_title(geo["title"], fontsize=9)
    ax.set_xlabel("Time (MJD)")
    ax.set_ylabel("Magnitude")
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.invert_yaxis()


def _draw_zoom(ax, geo) -> None:
    """Draw only the timing window of one eclipse, in MJD.

    Args:
        ax (matplotlib.axes.Axes): Target axes.
        geo (dict): Result of ``_eclipse_geometry``.
    """
    if geo["t_grid"] is not None:
        ax.plot(
            geo["t_grid"][geo["in_grid"]],
            geo["y_grid"][geo["in_grid"]],
            color=LINE_ADOPTED_COLOUR,
            lw=LINE_ADOPTED_WIDTH,
            zorder=5,
        )
    _scatter_window_points(ax, geo, 18)
    _mark_found_time(ax, geo)
    ax.set_xlim(geo["t_lo"], geo["t_hi"])
    ax.set_title(geo["title"], fontsize=9)
    ax.set_xlabel("Time (MJD)")
    ax.set_ylabel("Magnitude")
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=4, prune=None))
    ax.invert_yaxis()


def _timing_context(t_mjd, mag, sig, solution, timing, kind):
    """Common inputs of the timing pictures, validated.

    Args:
        t_mjd (ndarray): Times (MJD) that were timed.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        solution (ModelSolution): Solution whose template was used.
        timing (TimingOutput): Result of ``mikulasek_eb.time``.
        kind (str): ``primary`` or ``secondary``.

    Returns:
        tuple: ``(t, y, s, model, events)``.
    """
    t, y, s = (np.asarray(a, float) for a in (t_mjd, mag, sig))
    return t, y, s, solution.build_model(), _timing_events(t, solution, timing, kind)


def _one_eclipse_geo(t_mjd, mag, sig, solution, timing, kind, E):
    """Geometry of one named eclipse.

    Args:
        t_mjd (ndarray): Times (MJD).
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        solution (ModelSolution): Solution whose template was used.
        timing (TimingOutput): Result of ``mikulasek_eb.time``.
        kind (str): ``primary`` or ``secondary``.
        E (float): Eclipse number from T0.

    Returns:
        dict: Result of ``_eclipse_geometry``.

    Raises:
        ValueError: If the eclipse is not in the result.
    """
    t, y, s, model, events = _timing_context(t_mjd, mag, sig, solution, timing, kind)
    found = [ev for ev in events if ev["E"] == float(E)]
    if not found:
        raise ValueError(f"eclipse {E} is not in the {kind} timing result")
    return _eclipse_geometry(t, y, s, model, solution, timing, found[0], kind)


def plot_timing_window(
    t_mjd, mag, sig, solution, timing, kind: str, E: float, save_path=None
):
    """Draw the whole orbital cycle around one timed eclipse (requires matplotlib).

    The x-axis is MJD of the observations. Grey points are the rest of the cycle, coloured
    points are the fitted window, the template is at the fitted position (red in the
    window, grey outside), and a red line marks the found time. An eclipse that could
    not be timed is drawn in grey with the reason. Pass the same arrays and the same
    solution that ``mikulasek_eb.time`` received.

    Args:
        t_mjd (ndarray): Times (MJD) that were timed.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        solution (ModelSolution): Solution whose template was used.
        timing (TimingOutput): Result of ``mikulasek_eb.time``.
        kind (str): ``primary`` or ``secondary``.
        E (float): Eclipse number counted from T0, as in ``timing``.
        save_path (str | None): If set, also save the figure here.

    Returns:
        matplotlib.figure.Figure: The figure.

    Raises:
        ValueError: If ``timing`` does not belong to ``solution``, ``kind`` was not
            timed or ``E`` is not in the result.
    """
    geo = _one_eclipse_geo(t_mjd, mag, sig, solution, timing, kind, E)
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    _draw_cycle(ax, geo)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=130)
    return fig


def plot_timing_window_zoom(
    t_mjd, mag, sig, solution, timing, kind: str, E: float, save_path=None
):
    """Draw only the fitted window of one eclipse, in MJD (requires matplotlib).

    This is a separate picture, not a second panel of ``plot_timing_window``.

    Args:
        t_mjd (ndarray): Times (MJD) that were timed.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        solution (ModelSolution): Solution whose template was used.
        timing (TimingOutput): Result of ``mikulasek_eb.time``.
        kind (str): ``primary`` or ``secondary``.
        E (float): Eclipse number counted from T0.
        save_path (str | None): If set, also save the figure here.

    Returns:
        matplotlib.figure.Figure: The figure.

    Raises:
        ValueError: If the eclipse is not in the result.
    """
    geo = _one_eclipse_geo(t_mjd, mag, sig, solution, timing, kind, E)
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    _draw_zoom(ax, geo)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=130)
    return fig


def save_timing_windows(
    t_mjd, mag, sig, solution, timing, kind: str, out_dir, eclipses=None
) -> list:
    """Write one whole-cycle picture per eclipse, and a zoom in ``out_dir/zoom/``.

    Names follow MAVKA: ``E01_primary.png``. The zoom uses the same name in the
    subfolder ``zoom``. Each file is closed after saving.

    Args:
        t_mjd (ndarray): Times (MJD) that were timed.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        solution (ModelSolution): Solution whose template was used.
        timing (TimingOutput): Result of ``mikulasek_eb.time``.
        kind (str): ``primary`` or ``secondary``.
        out_dir (str | Path): Folder for the whole-cycle pictures; created if missing.
        eclipses (iterable of float | None): Eclipse numbers E (from T0) to draw;
            None draws all, timed and rejected.

    Returns:
        list[Path]: The whole-cycle files written.

    Raises:
        ValueError: If ``timing`` does not belong to ``solution``, ``kind`` was not
            timed, an eclipse in ``eclipses`` is not in the result or there is
            nothing to draw.
    """
    t, y, s, model, events = _timing_context(t_mjd, mag, sig, solution, timing, kind)
    if eclipses is not None:
        wanted = [float(e) for e in eclipses]
        missing = [e for e in wanted if e not in {ev["E"] for ev in events}]
        if missing:
            raise ValueError(f"eclipses not in the timing result: {missing}")
        events = [ev for ev in events if ev["E"] in wanted]
    if not events:
        raise ValueError(f"no {kind} eclipses to draw")
    out_dir = Path(out_dir)
    zoom_dir = out_dir / "zoom"
    out_dir.mkdir(parents=True, exist_ok=True)
    zoom_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ev in events:
        geo = _eclipse_geometry(t, y, s, model, solution, timing, ev, kind)
        stem = eclipse_picture_stem(kind, ev["E"])
        fig, ax = plt.subplots(figsize=(8.5, 4.6))
        _draw_cycle(ax, geo)
        fig.tight_layout()
        path = out_dir / f"{stem}.png"
        fig.savefig(path, dpi=130)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(5.5, 4.2))
        _draw_zoom(ax, geo)
        fig.tight_layout()
        fig.savefig(zoom_dir / f"{stem}.png", dpi=130)
        plt.close(fig)
        paths.append(path)
    return paths


def plot_timing_overview(t_mjd, mag, sig, solution, timings, save_path=None):
    """Draw the folded light curve with the timing windows of all eclipses (matplotlib).

    The points are folded on the template ephemeris over two cycles. The windows of the
    timed eclipses of every kind present in ``timings`` are shaded in one picture
    (primary and secondary in the colours of the package) and their points are coloured.
    The windows of rejected eclipses are not drawn on the fold, because the same phase
    holds both timed and rejected eclipses.

    Args:
        t_mjd (ndarray): Times (MJD) that were timed.
        mag (ndarray): Magnitudes.
        sig (ndarray): Magnitude errors.
        solution (ModelSolution): Solution whose template was used.
        timings (TimingOutput | sequence of TimingOutput): Results of ``time``. Pass
            one per recipe when the kinds were timed with different windows.
        save_path (str | None): If set, also save the figure here.

    Returns:
        matplotlib.figure.Figure: The figure.
    """
    if hasattr(timings, "timings"):
        timings = [timings]
    t, y = np.asarray(t_mjd, float), np.asarray(mag, float)
    ephem = solution.ephem
    th = phase(t, ephem["M0"], ephem["P0"], ephem.get("Pdot", 0.0))
    fold = th - np.floor(th + 0.5)
    fig, ax = plt.subplots(figsize=(11, 4.8))
    ax.scatter(fold, y, s=3, color="0.6", zorder=2)
    ax.scatter(fold + 1.0, y, s=3, color="0.6", zorder=2)
    labels = []
    for timing in timings:
        window = effective_window(
            timing.recipe.window, timing.recipe.allow_overlapping_windows
        )
        for kind in timing.timings:
            colour = (
                POINT_PRIMARY_COLOUR if kind == "primary" else POINT_SECONDARY_COLOUR
            )
            events = [
                ev for ev in _timing_events(t, solution, timing, kind) if ev["timed"]
            ]
            used = np.zeros(len(t), bool)
            for ev in events:
                used |= np.abs(th - ev["E_m0"]) < window
            centres = (0.0, 1.0) if kind == "primary" else (-0.5, 0.5, 1.5)
            for c in centres:
                ax.axvspan(c - window, c + window, color=colour, alpha=0.15, zorder=0)
            for shift in (0.0, 1.0):
                ax.scatter(fold[used] + shift, y[used], s=4, color=colour, zorder=3)
            ax.plot(
                [],
                [],
                color=colour,
                lw=6,
                alpha=0.4,
                label=f"{kind}, half-window {window:.3f}",
            )
            labels.append(kind)
    ax.set_xlim(-0.5, 1.5)
    ax.set_xlabel("Phase of the template ephemeris")
    ax.set_ylabel("Magnitude")
    ax.invert_yaxis()
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    ax.set_title("Folded light curve and timing windows")
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=130)
    return fig
