"""Magnitude difference of the two maxima of V1383 Tau in one cycle.

A cycle is kept only when the previous step timed both the primary maximum
and the secondary maximum of that cycle. The difference is the fitted-template
magnitude of the primary minus that of the secondary. Its error is the
root-sum-square of the two magnitude errors, because the two maxima are
separate fits. The time is the time of the primary maximum. Cycles with only
one maximum are omitted.

Run from ``mikulasek_eb``::

    python examples/v1383_max_diff.py

Writes one ``.dat`` file per sector in ``results/v1383_tau/timing_max_diff_dat/``.
"""

from __future__ import annotations

import csv
import logging
from decimal import Decimal
from pathlib import Path

import v1383_collect_dat as minima
import v1383_maxima as base

logger = logging.getLogger(__name__)

OUT = base.RESULTS / "timing_max_diff_dat"
KINDS = ("primary_maximum", "secondary_maximum")


def timed(path: Path) -> dict[str, dict[int, dict]]:
    """Timed maxima of one sector, keyed by kind and cycle.

    Args:
        path (Path): ``eclipses.csv`` from the magnitude step.

    Returns:
        dict[str, dict[int, dict]]: Primary and secondary rows that have a time.

    Raises:
        RuntimeError: If one cycle is timed twice for the same maximum.
    """
    found = {kind: {} for kind in KINDS}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["kind"] not in KINDS or not row["T_jd"]:
                continue
            if not row["mag"] or not row["err_mag"]:
                raise RuntimeError(
                    f"{path}: cycle {row['cycle']} has a time and no magnitude"
                )
            cycle = int(row["cycle"])
            if cycle in found[row["kind"]]:
                raise RuntimeError(
                    f"{path}: cycle {cycle} is repeated for {row['kind']}"
                )
            found[row["kind"]][cycle] = row
    return found


def pairs(found: dict[str, dict[int, dict]]) -> list[tuple[Decimal, Decimal, Decimal]]:
    """Cycles that hold both maxima, in order of the primary time.

    Args:
        found (dict[str, dict[int, dict]]): Timed rows from ``timed``.

    Returns:
        list[tuple[Decimal, Decimal, Decimal]]: Primary time, magnitude
        difference and its formal error.

    Raises:
        RuntimeError: If the secondary maximum of a cycle is not the later one,
        or a magnitude error is not positive.
    """
    both = sorted(set(found["primary_maximum"]) & set(found["secondary_maximum"]))
    rows = []
    for cycle in both:
        primary = found["primary_maximum"][cycle]
        secondary = found["secondary_maximum"][cycle]
        jd_primary = Decimal(f"{float(primary['T_jd']):.8f}")
        jd_secondary = Decimal(f"{float(secondary['T_jd']):.8f}")
        if jd_secondary <= jd_primary:
            raise RuntimeError(
                f"cycle {cycle}: the secondary maximum is not after the primary"
            )
        difference = Decimal(primary["mag"]) - Decimal(secondary["mag"])
        err_primary = Decimal(primary["err_mag"])
        err_secondary = Decimal(secondary["err_mag"])
        if err_primary <= 0 or err_secondary <= 0:
            raise RuntimeError(f"cycle {cycle}: a magnitude error is not positive")
        err = (err_primary**2 + err_secondary**2).sqrt()
        rows.append((jd_primary, difference, err))
    rows.sort()
    return rows


def write_one(
    path: Path, label: str, rows: list[tuple[Decimal, Decimal, Decimal]]
) -> None:
    """Write one sector file.

    Args:
        path (Path): Destination.
        label (str): Sector description.
        rows (list[tuple[Decimal, Decimal, Decimal]]): Primary time, difference
            and its error.
    """
    lines = [
        "Star: V1383 Tau (TIC 388995726), contact binary with a strong O'Connell effect",
        f"Light curve: {label}",
        "Difference of the fitted-template magnitudes of the primary maximum and the "
        "secondary maximum of the same cycle. A cycle with only one maximum is omitted.",
        "Columns: jd  mag_prim_minus_mag_sec  err_dmag",
        "jd is the time of the primary maximum",
        "mag_prim_minus_mag_sec is the primary magnitude minus the secondary magnitude",
        "err_dmag is the formal 1 sigma error of that difference, the root-sum-square "
        "of the two magnitude errors",
    ]
    text = [f"# {line}" for line in lines]
    for jd, difference, err in rows:
        text.append(f"{jd:.8f}  {difference:.6f}  {err:.6f}")
    path.write_text("\n".join(text) + "\n")


def main() -> None:
    """Write the difference file of every sector."""
    OUT.mkdir(parents=True, exist_ok=True)
    for lc, (stem, label) in minima.LIGHT_CURVES.items():
        path = base.RESULTS / f"{lc}_fit_bkg" / "maxima_mag" / "eclipses.csv"
        found = timed(path)
        rows = pairs(found)
        primary = len(found["primary_maximum"])
        secondary = len(found["secondary_maximum"])
        name = f"V1383Tau_{lc}_mag_prim_minus_mag_sec.dat"
        write_one(OUT / name, f"{stem}.vot ({label})", rows)
        logger.info(
            "%s: %d cycles with both maxima, %d primary and %d secondary timed",
            name,
            len(rows),
            primary,
            secondary,
        )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
