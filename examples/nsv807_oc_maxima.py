"""Common O-C of the NSV 807 minima and maxima.

The calculated times use the file ephemeris, which is the same in every light
curve. Minima are measured from the primary minimum and maxima from the
quadrature, phase 0.25. A secondary time is moved back by half a period before
the O-C is formed, so primary and secondary share a diagram. Sector 30 at
600 s is left out. Minima and maxima have different colours, and each sector
has its own symbol.

Run from ``mikulasek_eb`` with ``MPLBACKEND=Agg``::

    python examples/nsv807_oc_maxima.py

Writes ``results/nsv807/oc_maxima/oc.png`` and ``oc.csv``.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "nsv807"
MINIMA = RESULTS / "timing_dat"
MAXIMA = RESULTS / "timing_max_dat"
OUT = RESULTS / "oc_maxima"
PHASE_LIMIT = 0.2
# light curve, feature, kind, file, reference phase of that feature
FILES = (
    (
        "s30_120",
        "minimum",
        "primary",
        MINIMA / "NSV807_s30_120_mikulasek_window_primary.dat",
        0.0,
    ),
    (
        "s30_120",
        "minimum",
        "secondary",
        MINIMA / "NSV807_s30_120_mikulasek_window_secondary.dat",
        0.0,
    ),
    (
        "s4_1800",
        "minimum",
        "primary",
        MINIMA / "NSV807_s4_1800_mikulasek_window_primary.dat",
        0.0,
    ),
    (
        "s4_1800",
        "minimum",
        "secondary",
        MINIMA / "NSV807_s4_1800_mikulasek_window_secondary.dat",
        0.0,
    ),
    (
        "s97_200",
        "minimum",
        "primary",
        MINIMA / "NSV807_s97_200_mikulasek_window_primary.dat",
        0.0,
    ),
    (
        "s97_200",
        "minimum",
        "secondary",
        MINIMA / "NSV807_s97_200_mikulasek_window_secondary.dat",
        0.0,
    ),
    (
        "s30_120",
        "maximum",
        "primary",
        MAXIMA / "NSV807_s30_120_mikulasek_maximum_primary.dat",
        0.25,
    ),
    (
        "s30_120",
        "maximum",
        "secondary",
        MAXIMA / "NSV807_s30_120_mikulasek_maximum_secondary.dat",
        0.25,
    ),
    (
        "s4_1800",
        "maximum",
        "primary",
        MAXIMA / "NSV807_s4_1800_mikulasek_maximum_primary.dat",
        0.25,
    ),
    (
        "s4_1800",
        "maximum",
        "secondary",
        MAXIMA / "NSV807_s4_1800_mikulasek_maximum_secondary.dat",
        0.25,
    ),
    (
        "s97_200",
        "maximum",
        "primary",
        MAXIMA / "NSV807_s97_200_mikulasek_maximum_primary.dat",
        0.25,
    ),
    (
        "s97_200",
        "maximum",
        "secondary",
        MAXIMA / "NSV807_s97_200_mikulasek_maximum_secondary.dat",
        0.25,
    ),
)
COLOUR = {"minimum": "#0072B2", "maximum": "#D55E00"}
MARKER = {"s30_120": "o", "s4_1800": "s", "s97_200": "^"}
LABEL = {
    "s30_120": "Sector 30, 120 s",
    "s4_1800": "Sector 4, 1800 s",
    "s97_200": "Sector 97, 200 s",
}


def file_ephemeris() -> tuple[float, float]:
    """Read the file epoch and period, and require every sector to share them.

    Returns:
        tuple[float, float]: Epoch of primary minimum (JD) and period (days).

    Raises:
        RuntimeError: If a sector carries a different file ephemeris.
    """
    found = []
    for run in ("s30_protocol", "s30_600s_protocol", "s4_protocol", "s97_protocol"):
        path = next((RESULTS / run / "chosen").glob("*.json"))
        saved = json.loads(path.read_text())
        found.append(
            (
                float(saved["provenance"]["epoch_jd"]),
                float(saved["provenance"]["period_input_days"]),
            )
        )
    if any(item != found[0] for item in found[1:]):
        raise RuntimeError(f"the file ephemerides differ: {found}")
    return found[0]


def read_times(path: Path) -> list[tuple[float, float, str]]:
    """Read one collected maximum file.

    Args:
        path (Path): ``.dat`` file.

    Returns:
        list[tuple[float, float, str]]: JD, error in days, flag.

    Raises:
        RuntimeError: If a data line does not have three columns.
    """
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 3:
            raise RuntimeError(f"{path.name}: {line!r}")
        rows.append((float(parts[0]), float(parts[1]), parts[2]))
    if not rows:
        raise RuntimeError(f"{path.name} has no times")
    return rows


def to_seconds(days: np.ndarray) -> np.ndarray:
    """Convert days to seconds.

    Args:
        days (ndarray): Interval in days.

    Returns:
        ndarray: The same interval in seconds.
    """
    return (np.asarray(days, float) * u.day).to(u.s).value


def assemble(epoch: float, period: float) -> list[dict]:
    """Build the O-C of every kept minimum and maximum.

    Args:
        epoch (float): File epoch of primary minimum, JD.
        period (float): File period, days.

    Returns:
        list[dict]: One entry per time, with the cycle and the O-C in seconds.

    Raises:
        RuntimeError: If a time is far from the phase of its feature after the
        secondary shift.
    """
    points = []
    for light, feature, kind, path, reference in FILES:
        shift = 0.0 if kind == "primary" else -0.5 * period
        for jd, err, flag in read_times(path):
            used = jd + shift
            cycle = round((used - epoch) / period - reference)
            phase = (used - epoch) / period - cycle
            if abs(phase - reference) > PHASE_LIMIT:
                raise RuntimeError(f"{path.name}: phase {phase:.4f} is not a {feature}")
            oc_days = used - (epoch + (cycle + reference) * period)
            points.append(
                {
                    "light": light,
                    "feature": feature,
                    "kind": kind,
                    "cycle": cycle,
                    "oc_s": float(to_seconds(oc_days)),
                    "sigma_s": float(to_seconds(err)),
                    "flag": flag,
                }
            )
    return points


def draw(points: list[dict], epoch: float, period: float, path: Path) -> None:
    """Draw minima and maxima on one O-C diagram.

    Args:
        points (list[dict]): Rows from ``assemble``.
        epoch (float): File epoch, JD, for the title.
        period (float): File period, days, for the title.
        path (Path): Figure path.
    """
    fig, axis = plt.subplots(figsize=(9, 5.2))
    for feature in ("minimum", "maximum"):
        for light in ("s30_120", "s4_1800", "s97_200"):
            chosen = [
                p for p in points if p["feature"] == feature and p["light"] == light
            ]
            cycle = np.array([p["cycle"] for p in chosen])
            oc = np.array([p["oc_s"] for p in chosen])
            sigma = np.array([p["sigma_s"] for p in chosen])
            axis.errorbar(
                cycle,
                oc,
                sigma,
                fmt=MARKER[light],
                ms=4.5,
                color=COLOUR[feature],
                markerfacecolor=COLOUR[feature],
                ecolor=COLOUR[feature],
                elinewidth=0.5,
                label=f"{feature.capitalize()}, {LABEL[light]}",
            )
    axis.axhline(0.0, color="0.6", lw=0.7)
    axis.set_xlabel("Cycle number")
    axis.set_ylabel("O-C (s)")
    axis.set_title(
        "NSV 807\n"
        f"File epoch JD {epoch:.2f}, P = {period:.7f} d. "
        "Secondary times shifted back by half a cycle."
    )
    axis.legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def write_table(points: list[dict], path: Path) -> None:
    """Write the plotted points.

    Args:
        points (list[dict]): Rows from ``assemble``.
        path (Path): CSV path.
    """
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["light", "feature", "kind", "cycle", "oc_s", "sigma_s", "flag"],
        )
        writer.writeheader()
        for point in points:
            writer.writerow(
                {
                    "light": point["light"],
                    "feature": point["feature"],
                    "kind": point["kind"],
                    "cycle": point["cycle"],
                    "oc_s": f"{point['oc_s']:.3f}",
                    "sigma_s": f"{point['sigma_s']:.3f}",
                    "flag": point["flag"],
                }
            )


def main() -> None:
    """Plot the common O-C and write the points."""
    epoch, period = file_ephemeris()
    points = assemble(epoch, period)
    OUT.mkdir(parents=True, exist_ok=True)
    draw(points, epoch, period, OUT / "oc.png")
    write_table(points, OUT / "oc.csv")
    logger.info(
        "epoch JD %.5f, P %.7f d, %d points, wrote %s",
        epoch,
        period,
        len(points),
        OUT / "oc.png",
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
