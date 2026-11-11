"""Magnitude difference of the AP maxima of V1383 Tau in one cycle.

The times and magnitudes are the AP results. A cycle is kept only when both
the primary maximum and the secondary maximum of that cycle have a time. The
cycle is the one of the supplied window centre. The difference is the primary
magnitude minus the secondary magnitude. Its error is the root-sum-square of
the two magnitude errors. The time is the AP time of the primary maximum.

Run from ``mikulasek_eb``::

    python examples/v1383_mavka_max_diff.py

Writes one ``.dat`` file per sector in ``results/v1383_tau/mavka_ap_max_diff/``.
"""

from __future__ import annotations

import csv
import logging
from decimal import Decimal
from pathlib import Path

import v1383_maxima as base

logger = logging.getLogger(__name__)

MAVKA = Path("/home/voz/projects/UPJS/Mavka_user_folder/output/v1383_tau_maxima")
CENTRES = Path("/home/voz/projects/UPJS/agent_mailroom/v1383_tau_maxima_centres")
OUT = base.RESULTS / "mavka_ap_max_diff"
PERIOD = Decimal("0.35064")
HALF = {"primary": Decimal("0.145"), "secondary": Decimal("0.128")}
RUNS = (
    ("s42", "s42_600", "s42_fit_bkg"),
    ("s43", "s43_600", "s43_fit_bkg"),
    ("s44", "s44_600", "s44_fit_bkg"),
    ("s70", "s70_200", "s70_fit_bkg"),
    ("s71", "s71_200", "s71_fit_bkg"),
)


def centres_of(path: Path) -> list[str]:
    """Window centres, in file order.

    Args:
        path (Path): Centre list.

    Returns:
        list[str]: Full JD strings.

    Raises:
        RuntimeError: If the file holds no centre.
    """
    found = [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if not found:
        raise RuntimeError(f"{path.name}: no centres")
    return found


def read_mavka(path: Path) -> list[tuple]:
    """Accepted rows and rejected centres, in file order.

    Args:
        path (Path): AP ``.dat`` file.

    Returns:
        list[tuple]: ``("ok", jd, err, mag, err_mag)`` or ``("reject", jd)``.

    Raises:
        RuntimeError: If a data line does not hold four columns.
    """
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            parts = line[1:].split()
            if len(parts) >= 2 and parts[1] == "rejected:":
                rows.append(("reject", parts[0]))
            continue
        parts = line.split()
        if len(parts) != 4:
            raise RuntimeError(f"{path.name}: {line!r}")
        rows.append(("ok", *parts))
    return rows


def assign(
    centres: list[str], rows: list[tuple], half: Decimal, label: str
) -> dict[str, tuple]:
    """Give each centre its AP row.

    Args:
        centres (list[str]): Window centres.
        rows (list[tuple]): AP rows in the same order.
        half (Decimal): Half-window in cycles.
        label (str): Sector and maximum, for the error text.

    Returns:
        dict[str, tuple]: Accepted row by centre. A rejected centre is absent.

    Raises:
        RuntimeError: If the file does not follow the centre list, or a fitted
        time lies outside the window of its centre.
    """
    if len(rows) != len(centres):
        raise RuntimeError(f"{label}: {len(rows)} rows for {len(centres)} centres")
    limit = half * PERIOD
    accepted = {}
    for centre, row in zip(centres, rows):
        if row[0] == "reject":
            if row[1] != centre:
                raise RuntimeError(f"{label}: rejected {row[1]}, centre {centre}")
            continue
        shift = abs(Decimal(row[1]) - Decimal(centre))
        if shift > limit:
            raise RuntimeError(
                f"{label}: fitted time {row[1]} is outside the window of {centre}"
            )
        accepted[centre] = row
    return accepted


def cycles_of(folder: Path) -> dict[tuple[str, str], int]:
    """Cycle number of each supplied centre.

    Args:
        folder (Path): Sector ``maxima`` folder.

    Returns:
        dict[tuple[str, str], int]: Keyed by maximum and centre JD.

    Raises:
        RuntimeError: If one centre belongs to two cycles.
    """
    found = {}
    with (folder / "eclipses.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            if not row["T_jd"]:
                continue
            kind = "primary" if row["kind"] == "primary_maximum" else "secondary"
            if row["kind"] not in ("primary_maximum", "secondary_maximum"):
                continue
            centre = f"{float(row['T_jd']):.8f}"
            key = (kind, centre)
            cycle = int(row["cycle"])
            if key in found and found[key] != cycle:
                raise RuntimeError(f"{folder.name}: centre {centre} has two cycles")
            found[key] = cycle
    return found


def write_one(path: Path, label: str, rows: list[tuple[str, Decimal, Decimal]]) -> None:
    """Write one sector file.

    Args:
        path (Path): Destination.
        label (str): Light-curve description.
        rows (list[tuple[str, Decimal, Decimal]]): Primary time, difference
            and its error.
    """
    lines = [
        "Star: V1383 Tau (TIC 388995726), contact binary with a strong O'Connell effect",
        f"Light curve: {label}",
        "Method: MAVKA AP",
        "Difference of the AP magnitudes of the primary maximum and the secondary "
        "maximum of the same cycle. A cycle with only one maximum is omitted.",
        "Columns: jd  mag_prim_minus_mag_sec  err_dmag",
        "jd is the AP time of the primary maximum",
        "mag_prim_minus_mag_sec is the primary magnitude minus the secondary magnitude",
        "err_dmag is the formal 1 sigma error of that difference, the root-sum-square "
        "of the two magnitude errors",
    ]
    text = [f"# {line}" for line in lines]
    for jd, difference, err in rows:
        text.append(f"{jd}  {difference:.8f}  {err:.8f}")
    path.write_text("\n".join(text) + "\n")


def one_sector(lc: str, label: str, folder: str) -> None:
    """Pair one sector and write its difference file.

    Args:
        lc (str): Sector label used in the file name.
        label (str): MAVKA file label.
        folder (str): Fit folder under ``results/v1383_tau``.
    """
    cycle = cycles_of(base.RESULTS / folder / "maxima")
    accepted = {}
    for kind in ("primary", "secondary"):
        accepted[kind] = assign(
            centres_of(CENTRES / f"{label}_{kind}.dat"),
            read_mavka(MAVKA / f"{label}_AP_{kind}.dat"),
            HALF[kind],
            f"{label} {kind}",
        )
    by_cycle = {"primary": {}, "secondary": {}}
    for kind in ("primary", "secondary"):
        for centre, row in accepted[kind].items():
            if (kind, centre) not in cycle:
                raise RuntimeError(
                    f"{label}: centre {centre} is not in the saved timing"
                )
            by_cycle[kind][cycle[(kind, centre)]] = row
    both = sorted(set(by_cycle["primary"]) & set(by_cycle["secondary"]))
    rows = []
    for n in both:
        primary = by_cycle["primary"][n]
        secondary = by_cycle["secondary"][n]
        if Decimal(secondary[1]) <= Decimal(primary[1]):
            raise RuntimeError(
                f"{label} cycle {n}: the secondary is not after the primary"
            )
        difference = Decimal(primary[3]) - Decimal(secondary[3])
        err_primary = Decimal(primary[4])
        err_secondary = Decimal(secondary[4])
        if err_primary <= 0 or err_secondary <= 0:
            raise RuntimeError(f"{label} cycle {n}: a magnitude error is not positive")
        rows.append(
            (primary[1], difference, (err_primary**2 + err_secondary**2).sqrt())
        )
    rows.sort()
    name = f"V1383Tau_{lc}_AP_mag_prim_minus_mag_sec.dat"
    write_one(OUT / name, label, rows)
    logger.info("%s: %d cycles with both maxima", name, len(rows))


def main() -> None:
    """Write the AP difference file of every sector."""
    OUT.mkdir(parents=True, exist_ok=True)
    for lc, label, folder in RUNS:
        one_sector(lc, label, folder)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
