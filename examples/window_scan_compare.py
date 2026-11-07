"""Plot the MAVKA 120 s window scan against the Mikulasek scan.

Reads the tables Mavka_user has already written. Does not fit anything. MAVKA
curves are drawn only up to 3 D, the half-window at which the template of these
stars is down to about five percent of its central depth. The Mikulasek curve
is the full scan. One figure per star, and one overview with the pure
(noise-free) light curve.

Run from the repository root::

    MPLBACKEND=Agg python examples/window_scan_compare.py
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from volightcurve import read_lightcurve

from window_scan_mik import GENERATED, MAVKA_HALF_120, MAX_HALF, OUT

logger = logging.getLogger(__name__)
MAVKA_ROOT = Path(
    "/home/voz/projects/UPJS/Mavka_user_folder/output/window_scan_mavka_120"
)
COMPARE = OUT.parent / "window_scan_compare"
CLEAN = GENERATED / "results" / "clean"
MAVKA_LIMIT_D = 3.0
# G7 primary: the pure-curve shoulder, where the steep wall has given way to
# the slow outer rise. The two profile-1 components are down to about one
# percent of their combined depth here.
G7_PRIMARY_DRAW_D = 4.2
METHODS = ("AP", "WSAP", "WSL")
COLOUR = {
    "Mikulasek": "red",
    "AP": "tab:blue",
    "WSAP": "tab:orange",
    "WSL": "tab:green",
}
KIND_TITLE = {"primary": "Primary", "secondary": "Secondary"}
STAR_LABEL = {
    "G1": "Well-detached",
    "G2": "Total eclipse",
    "G4": "Algol",
    "G5": "W UMa",
    "G6": "O'Connell",
    "G7": "Shallow secondary",
    "G8": "Eccentric",
}
PURE_CURVE = {
    "G1": CLEAN / "detached" / "detached_none_nooconnell_G1.vot",
    "G2": CLEAN / "detached" / "detached_none_nooconnell_G2.vot",
    "G4": CLEAN / "semidetached" / "semidetached_none_nooconnell_G4.vot",
    "G5": CLEAN / "contact" / "contact_overcontact_nooconnell_G5.vot",
    "G6": CLEAN / "contact" / "contact_overcontact_oconnell_G6.vot",
    "G7": CLEAN / "semidetached" / "semidetached_none_nooconnell_G7.vot",
    "G8": CLEAN / "detached" / "detached_none_nooconnell_G8.vot",
}


def star_title(star: str) -> str:
    """Short description of a physics slot, then its code.

    Args:
        star (str): Physics id.

    Returns:
        str: For example ``Eccentric (G8)``.
    """
    return f"{STAR_LABEL[star]} ({star})"


def read_csv(path: Path) -> list[dict]:
    """Read a scan table.

    Args:
        path (Path): CSV file.

    Returns:
        list[dict]: Rows. Empty if the file is missing.
    """
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def mavka_rows(star: str) -> list[dict]:
    """MAVKA scan rows, including the later WSAP and WSL file.

    Args:
        star (str): Physics id.

    Returns:
        list[dict]: Rows from ``scan.csv`` and, when present,
        ``wsap_wsl_scan.csv``.
    """
    folder = MAVKA_ROOT / star
    return read_csv(folder / "scan.csv") + read_csv(folder / "wsap_wsl_scan.csv")


def mavka_limit(star: str, kind: str) -> float | None:
    """Largest half-window, in units of D, drawn for one MAVKA curve.

    G5 and G6 are contact stars, so there is no continuum to stop at. Their
    saved scans already reach 0.25 cycles and are drawn in full. The G7
    primary is drawn to 4.2 D1, the shoulder past the wall. Every other
    eclipse stays at 3 D.

    Args:
        star (str): Physics id.
        kind (str): ``primary`` or ``secondary``.

    Returns:
        float | None: Limit in units of D. ``None`` draws the whole saved scan.
    """
    if star in ("G5", "G6"):
        return None
    if star == "G7" and kind == "primary":
        return G7_PRIMARY_DRAW_D
    return MAVKA_LIMIT_D


def series(
    rows: list[dict],
    kind: str,
    method: str | None,
    column: str,
    limit_d: float | None = None,
) -> tuple[list[float], list[float]]:
    """Finite (half / D, value) pairs for one curve.

    Args:
        rows (list[dict]): Scan rows.
        kind (str): ``primary`` or ``secondary``.
        method (str | None): Method name. ``None`` selects the Mikulasek table,
            which has no method column.
        column (str): ``rms_oc_s`` or ``mean_oc_s``.
        limit_d (float | None): Largest half-window, in units of D, to draw.

    Returns:
        tuple[list[float], list[float]]: Abscissa in units of D, and the column.
    """
    xs, ys = [], []
    for row in rows:
        if row["kind"] != kind:
            continue
        if method is not None and row.get("method") != method:
            continue
        text = row.get(column, "")
        if text == "":
            continue
        x = float(row["half_over_D"])
        if limit_d is not None and x > limit_d + 1e-6:
            continue
        xs.append(x)
        ys.append(float(text))
    return xs, ys


def widths(star: str) -> dict[str, float]:
    """Mikulasek D of each eclipse at 120 s.

    Args:
        star (str): Physics id.

    Returns:
        dict[str, float]: D in cycles, by eclipse kind.

    Raises:
        KeyError: If the grid has no row for that star.
    """
    found = {}
    with (OUT / "grid_120.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["physics_id"] == star:
                found[row["kind"]] = float(row["D_cycles"])
    if set(found) != {"primary", "secondary"}:
        raise KeyError(f"No 120 s D for {star}.")
    return found


def pure_curve(star: str) -> tuple[list[np.ndarray], list[np.ndarray], float]:
    """Noise-free one-orbit curve, repeated over phase -1.5 to 1.5.

    The clean parent has no error column, so the flux and the time are read
    from the VOTable directly. The time column is days from the file epoch.

    Args:
        star (str): Physics id.

    Returns:
        tuple[list[ndarray], list[ndarray], float]: One phase array and one
        flux array per cycle copy, and the period in days.
    """
    path = PURE_CURVE[star]
    curve = read_lightcurve(path, filename=path.name)
    time = np.asarray(curve.table["time"], dtype=float)
    flux = np.asarray(curve.table["flux"], dtype=float)
    period = float(curve.table.meta["period"])
    epoch = float(curve.table.meta.get("epoch", 0.0))
    phase = (time - epoch) / period
    phases, fluxes = [], []
    for cycle in (-2, -1, 0, 1):
        shifted = phase + cycle
        keep = (shifted >= -1.5) & (shifted <= 1.5)
        order = np.argsort(shifted[keep])
        phases.append(shifted[keep][order])
        fluxes.append(flux[keep][order])
    return phases, fluxes, period


def draw_scans(ax, mik_rows, mavka, kind, column, d_of, adopted, star: str) -> None:
    """Draw one rms or mean panel.

    Args:
        ax: Matplotlib axes.
        mik_rows (list[dict]): Mikulasek scan rows.
        mavka (list[dict]): MAVKA scan rows.
        kind (str): ``primary`` or ``secondary``.
        column (str): ``rms_oc_s`` or ``mean_oc_s``.
        d_of (dict[str, float]): D by eclipse kind.
        adopted (dict[str, float]): Adopted MAVKA half-width in cycles.
        star (str): Physics id. Chooses how far the MAVKA curves are drawn.
    """
    x, y = series(mik_rows, kind, None, column)
    if x:
        ax.plot(x, y, "-", color=COLOUR["Mikulasek"], label="Mikulasek", lw=1.4)
    limit = mavka_limit(star, kind)
    for name in METHODS:
        x, y = series(mavka, kind, name, column, limit)
        if x:
            ax.plot(x, y, "-", color=COLOUR[name], label=name, lw=1.2)
    mark = adopted[kind] / d_of[kind]
    ax.axvline(mark, color="0.35", ls="--", lw=1, label=f"adopted {mark:.2f} D")
    ax.axvline(MAVKA_LIMIT_D, color="0.55", ls=":", lw=1, label="3 D")
    ax.axvline(MAX_HALF / d_of[kind], color="0.15", ls="-.", lw=1, label="0.25 cycles")
    ax.grid(alpha=0.3)


def draw_mean_twin(ax, mik_rows, mavka, kind, star: str, show_label: bool):
    """Draw the mean O-C of every method on a right-hand axis.

    The lines use the colours of the rms curves, thinner and fainter, and
    the same drawing limits. The right axis is symmetric logarithmic, linear
    inside +-10 s, so that a drift of hundreds of seconds and a few seconds
    can be read on one axis.

    Args:
        ax: Left axes, which holds the rms.
        mik_rows (list[dict]): Mikulasek scan rows.
        mavka (list[dict]): MAVKA scan rows.
        kind (str): ``primary`` or ``secondary``.
        star (str): Physics id. Chooses how far the MAVKA curves are drawn.
        show_label (bool): Write the axis label. Done on the secondary panel.

    Returns:
        The twin axes.
    """
    twin = ax.twinx()
    x, y = series(mik_rows, kind, None, "mean_oc_s")
    if x:
        twin.plot(x, y, "-", color=COLOUR["Mikulasek"], lw=0.8, alpha=0.5)
    limit = mavka_limit(star, kind)
    for name in METHODS:
        x, y = series(mavka, kind, name, "mean_oc_s", limit)
        if x:
            twin.plot(x, y, "-", color=COLOUR[name], lw=0.8, alpha=0.5)
    twin.axhline(0.0, color="0.6", lw=0.5)
    twin.set_yscale("symlog", linthresh=10.0)
    if show_label:
        twin.set_ylabel("mean O-C (s), thin lines")
    return twin


def draw_star(star: str, mavka: list[dict], mik_rows: list[dict], path: Path) -> None:
    """Draw rms and mean O-C against half-window over D for one star.

    Args:
        star (str): Physics id.
        mavka (list[dict]): MAVKA scan rows.
        mik_rows (list[dict]): Mikulasek scan rows.
        path (Path): Output PNG.
    """
    d_of = widths(star)
    adopted = MAVKA_HALF_120[star]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex="row")
    for i, kind in enumerate(("primary", "secondary")):
        draw_scans(axes[i, 0], mik_rows, mavka, kind, "rms_oc_s", d_of, adopted, star)
        draw_scans(axes[i, 1], mik_rows, mavka, kind, "mean_oc_s", d_of, adopted, star)
        for ax in axes[i]:
            ax.set_xlabel("half-window / D")
            ax.legend(fontsize=8)
        axes[i, 0].set_yscale("log")
        axes[i, 0].set_ylabel("rms of O-C (s)")
        axes[i, 1].set_ylabel("mean O-C (s)")
        axes[i, 0].set_title(KIND_TITLE[kind])
        if not any(
            line.get_label() in METHODS or line.get_label() == "Mikulasek"
            for line in axes[i, 1].lines
        ):
            axes[i, 1].text(
                0.5, 0.5, "mean not saved", transform=axes[i, 1].transAxes, ha="center"
            )
    fig.suptitle(f"{star_title(star)}, 120 s")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
    logger.info("wrote %s", path)


def _flux_at(
    phases: list[np.ndarray], fluxes: list[np.ndarray], target: float
) -> float:
    """Flux of the pure curve at one phase.

    Args:
        phases (list[ndarray]): Phase of each plotted copy, each increasing.
        fluxes (list[ndarray]): Flux of each plotted copy.
        target (float): Phase. A phase outside the plotted range is folded
            by whole cycles, because the curve repeats.

    Returns:
        float: Interpolated flux.

    Raises:
        RuntimeError: If the plotted curve covers less than one cycle.
    """
    phase = np.concatenate(phases)
    flux = np.concatenate(fluxes)
    order = np.argsort(phase)
    phase, flux = phase[order], flux[order]
    if phase[-1] - phase[0] < 1.0:
        raise RuntimeError("Pure curve covers less than one cycle.")
    folded = target
    while folded < phase[0]:
        folded += 1.0
    while folded > phase[-1]:
        folded -= 1.0
    return float(np.interp(folded, phase, flux))


def mark_chord(
    ax,
    phases: list[np.ndarray],
    fluxes: list[np.ndarray],
    centre: float,
    half: float,
    color: str,
    label: str | None,
) -> None:
    """Draw one horizontal chord of length ``2 * half`` on the pure curve.

    The chord is centred on ``centre``. Its flux is the mean of the pure curve
    at the two ends. ``half`` is a half-width, so a Mikulasek D is drawn with
    length 2 D, and the MAVKA half-window of 3 D1 is drawn with length 6 D1.

    Args:
        ax: Matplotlib axes.
        phases (list[ndarray]): Phase of each plotted copy.
        fluxes (list[ndarray]): Flux of each plotted copy.
        centre (float): Phase of the chord centre.
        half (float): Half-length of the chord, in cycles.
        color (str): Matplotlib colour.
        label (str | None): Legend text. ``None`` draws the chord with no entry.
    """
    left, right = centre - half, centre + half
    level = 0.5 * (_flux_at(phases, fluxes, left) + _flux_at(phases, fluxes, right))
    ax.plot(
        [left, right],
        [level, level],
        color=color,
        lw=1.8,
        solid_capstyle="butt",
        zorder=4,
        label=label,
    )


def draw_pure(ax, star: str, d_of: dict[str, float], show_legend: bool) -> None:
    """Draw the noise-free light curve over phase -1 to 1, with three chords.

    The red chord is 2 D1 long and is centred on phase 0. The green chord is
    2 D2 long and is centred on phase 0.5. The blue chord is 2 x 3 D1 long,
    the full width of the MAVKA half-window of 3 primary D, centred on phase 0.
    D1 and D2 are the Mikulasek widths of this star.

    Args:
        ax: Matplotlib axes.
        star (str): Physics id.
        d_of (dict[str, float]): Mikulasek D in cycles, by eclipse kind.
        show_legend (bool): Draw the one legend shared by every light curve.
    """
    phases, fluxes, period = pure_curve(star)
    for phase, flux in zip(phases, fluxes):
        ax.plot(phase, flux, color="0.15", lw=0.7, label="_nolegend_")
    d1 = d_of["primary"]
    d2 = d_of["secondary"]
    mark_chord(
        ax, phases, fluxes, 0.0, 3.0 * d1, "blue", "2 x 3 D1" if show_legend else None
    )
    mark_chord(ax, phases, fluxes, 0.0, d1, "red", "2 D1" if show_legend else None)
    mark_chord(ax, phases, fluxes, 0.5, d2, "green", "2 D2" if show_legend else None)
    ax.set_xlim(-1.0, 1.0)
    ax.set_xlabel("phase")
    ax.set_ylabel("flux")
    ax.set_title(star_title(star))
    ax.grid(alpha=0.3)
    ax.text(
        0.98,
        0.04,
        f"P = {period:g} d",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
    )
    if show_legend:
        ax.legend(fontsize=7, loc="upper right")


def report(star: str, rows: list[dict], method: str | None, label: str) -> None:
    """Log the rms at 1 D, 2 D, 3 D and at the last drawn window.

    Args:
        star (str): Physics id.
        rows (list[dict]): Scan rows.
        method (str | None): Method, or ``None`` for Mikulasek.
        label (str): Name used in the log line.
    """
    for kind in ("primary", "secondary"):
        limit = None if method is None else mavka_limit(star, kind)
        x, y = series(rows, kind, method, "rms_oc_s", limit)
        if not x:
            logger.info("%s %s %s: no rms", star, label, kind)
            continue
        pairs = list(zip(x, y))
        bits = []
        for target in (1.0, 2.0, 3.0):
            near = min(pairs, key=lambda item: abs(item[0] - target))
            if abs(near[0] - target) <= 0.25:
                bits.append(f"{near[0]:.1f} D {near[1]:.2f} s")
        last = pairs[-1]
        bits.append(f"last {last[0]:.1f} D {last[1]:.2f} s")
        logger.info("%s %s %s: %s", star, label, kind, "; ".join(bits))


def main() -> None:
    """Write one comparison figure per star and one overview of the rms."""
    COMPARE.mkdir(parents=True, exist_ok=True)
    stars = [
        path.name
        for path in sorted(MAVKA_ROOT.iterdir())
        if (path / "scan.csv").is_file()
    ]
    overview, overview_axes = plt.subplots(
        len(stars), 3, figsize=(15, 2.35 * len(stars)), sharex=False
    )
    for row, star in enumerate(stars):
        mavka = mavka_rows(star)
        mik_rows = read_csv(OUT / f"{star}_c120" / "scan.csv")
        draw_star(star, mavka, mik_rows, COMPARE / f"{star}_rms_and_mean.png")
        report(star, mik_rows, None, "Mikulasek")
        for method in METHODS:
            report(star, mavka, method, method)
        d_of = widths(star)
        adopted = MAVKA_HALF_120[star]
        for col, kind in enumerate(("primary", "secondary")):
            ax = overview_axes[row][col]
            draw_scans(ax, mik_rows, mavka, kind, "rms_oc_s", d_of, adopted, star)
            ax.set_yscale("log")
            draw_mean_twin(ax, mik_rows, mavka, kind, star, show_label=(col == 1))
            ax.set_title(f"{star_title(star)}, {kind}")
            if row == len(stars) - 1:
                ax.set_xlabel("half-window / D")
            if col == 0:
                ax.set_ylabel("rms (s)")
            if row == 0 and col == 0:
                ax.legend(fontsize=7)
        draw_pure(overview_axes[row][2], star, d_of, show_legend=(row == 0))
    overview.suptitle(
        "120 s. Red: Mikulasek, full scan. "
        "Dashed: adopted window. Dotted: 3 D. Dash-dot: 0.25 cycles. "
        "Thick: rms (left axis). Thin: mean O-C (right axis, symlog). "
        "Right: noise-free curve."
    )
    overview.tight_layout()
    overview_path = COMPARE / "overview_rms.png"
    overview.savefig(overview_path, dpi=140)
    plt.close(overview)
    logger.info("wrote %s", overview_path)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
