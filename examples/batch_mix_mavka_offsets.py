"""Collate the first-batch MAVKA times with the PHOEBE parameters on the files.

No fit is run and no light curve is timed. Each accepted time in the MAVKA
tables is matched to the nearest geometric conjunction of the same kind, and
to the pure-curve minimum stored beside it. Both clocks are the file clock.
Rejected rows are counted and are not given a time.

Run from the repository root::

    python examples/batch_mix_mavka_offsets.py

Writes ``results/synthetic_scoring/batch_mix_mavka/physics.csv`` and
``offsets.csv``. The MAVKA tables are not modified.
"""

from __future__ import annotations

import csv
import logging
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from mikulasek_eb.synthetic_check import _array, _num, load_synthetic

logger = logging.getLogger(__name__)

GENERATED = Path("/home/voz/projects/UPJS/EB_synthetic_lc")
MANIFEST = GENERATED / "results" / "mix" / "manifest.csv"
TABLES = Path("/home/voz/projects/UPJS/Mavka_user_folder/output/batch_mix/tables")
OUT = (
    Path(__file__).resolve().parents[1]
    / "results"
    / "synthetic_scoring"
    / "batch_mix_mavka"
)
DAY_S = 86400.0
METHODS = ("AP", "WSAP", "WSL")
PHYSICS = ("G1", "G2", "G4", "G5", "G6", "G7", "G8")


def _manifest_rows() -> list[dict[str, str]]:
    """Return the manifest rows that the generator marked as written.

    Returns:
        list[dict[str, str]]: One dict per light curve, with ``fail`` equal to 0.
    """
    with MANIFEST.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if int(row["fail"]) == 0]
    if len(rows) != 22:
        raise RuntimeError(f"expected 22 written light curves, found {len(rows)}")
    return rows


def _split(value: object) -> list[str]:
    """Split a comma-separated metadata string.

    Args:
        value (object): The stored parameter.

    Returns:
        list[str]: The stripped items.
    """
    return [part.strip() for part in str(value).split(",")]


def _physics_row(meta: dict) -> dict[str, str]:
    """Collect the PHOEBE and eclipse parameters used in the collation.

    Args:
        meta (dict): Metadata from ``load_synthetic``.

    Returns:
        dict[str, str]: One row for ``physics.csv``.
    """
    period = _num(meta, "P_true_days")
    fwhm_p = _num(meta, "eclipse_fwhm_primary_days")
    fwhm_s = _num(meta, "eclipse_fwhm_secondary_days")
    requiv_p = _num(meta, "requiv.primary.star.component")
    requiv_s = _num(meta, "requiv.secondary.star.component")
    lobe_p = _num(meta, "requiv_max.primary.star.component")
    lobe_s = _num(meta, "requiv_max.secondary.star.component")
    return {
        "physics_id": str(meta["physics_id"]),
        "eb_class": str(meta["eb_class"]),
        "P_true_days": f"{period:.6g}",
        "q": f"{_num(meta, 'q.binary.orbit.component'):.3f}",
        "incl_deg": f"{_num(meta, 'incl.binary.orbit.component'):.1f}",
        "ecc": f"{_num(meta, 'ecc.binary.orbit.component'):.3f}",
        "per0_deg": f"{_num(meta, 'per0.binary.orbit.component'):.1f}",
        "teff_primary_K": f"{_num(meta, 'teff.primary.star.component'):.0f}",
        "teff_secondary_K": f"{_num(meta, 'teff.secondary.star.component'):.0f}",
        "fill_primary": f"{requiv_p / lobe_p:.3f}",
        "fill_secondary": f"{requiv_s / lobe_s:.3f}",
        "fillout_factor": str(meta.get("fillout_factor", "")),
        "total_eclipse": str(meta["total_eclipse"]),
        "has_spots": str(meta["has_spots"]),
        "has_oconnell": str(meta["has_oconnell"]),
        "eclipse_depth_mag": f"{_num(meta, 'eclipse_depth_mag'):.3f}",
        "fwhm_primary_days": f"{fwhm_p:.6g}",
        "fwhm_secondary_days": f"{fwhm_s:.6g}",
        "half_fwhm_primary_phase": f"{0.5 * fwhm_p / period:.4f}",
        "half_fwhm_secondary_phase": f"{0.5 * fwhm_s / period:.4f}",
        "conj_minus_min_pure_primary_s": f"{_num(meta, 'conj_minus_min_pure_primary_s'):.2f}",
        "conj_minus_min_pure_secondary_s": f"{_num(meta, 'conj_minus_min_pure_secondary_s'):.2f}",
    }


def _match(measured: np.ndarray, truth: np.ndarray, period_days: float) -> np.ndarray:
    """Offset of each measured time from the nearest truth time, in seconds.

    Args:
        measured (np.ndarray): Accepted times, days on the file clock.
        truth (np.ndarray): Conjunction or pure-minimum times of one kind.
        period_days (float): True period. A match farther than half a period is a
            pairing error, not a large timing residual.

    Returns:
        np.ndarray: ``measured`` minus the nearest truth time, in seconds.

    Raises:
        RuntimeError: If the nearest truth time is more than half a period away.
    """
    nearest = np.argmin(np.abs(measured[:, None] - truth[None, :]), axis=1)
    delta_days = measured - truth[nearest]
    if np.any(np.abs(delta_days) > 0.5 * period_days):
        raise RuntimeError(
            "a measured time is more than half a period from every truth time"
        )
    return delta_days * DAY_S


def main() -> None:
    """Write the physics table and the per-file offset summary."""
    rows = _manifest_rows()
    by_name = {Path(row["path"]).name: GENERATED / row["path"] for row in rows}
    loaded: dict[str, object] = {}

    def synthetic(name: str):
        if name not in loaded:
            loaded[name] = load_synthetic(by_name[name])
        return loaded[name]

    OUT.mkdir(parents=True, exist_ok=True)
    physics_path = OUT / "physics.csv"
    seen: set[str] = set()
    physics_rows: list[dict[str, str]] = []
    for row in rows:
        if row["physics_id"] in seen or row["recipe_id"] == "B22":
            continue
        seen.add(row["physics_id"])
        physics_rows.append(_physics_row(synthetic(Path(row["path"]).name).meta))
    physics_rows.sort(key=lambda item: PHYSICS.index(item["physics_id"]))
    with physics_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(physics_rows[0]))
        writer.writeheader()
        writer.writerows(physics_rows)
    logger.info("wrote %s (%d stars)", physics_path, len(physics_rows))

    offset_rows: list[dict[str, str]] = []
    for method in METHODS:
        grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
        reasons: Counter[str] = Counter()
        with (TABLES / f"{method}.csv").open(newline="") as handle:
            for record in csv.DictReader(handle):
                grouped[record["file"]].append(record)
                if record["status"] != "accepted":
                    reasons[record["reason"]] += 1
        logger.info("%s rejected reasons: %s", method, dict(reasons))
        for name, group in grouped.items():
            meta = synthetic(name).meta
            kinds = _split(meta["eclipse_kind"])
            conj_all = _array(meta, "eclipse_T_conj_days")
            pure_all = _array(meta, "eclipse_T_min_pure_days")
            period = _num(meta, "P_true_days")
            for kind in ("primary", "secondary"):
                accepted = [
                    record
                    for record in group
                    if record["kind"] == kind and record["status"] == "accepted"
                ]
                fit_fail = [
                    record
                    for record in group
                    if record["kind"] == kind
                    and record["status"] != "accepted"
                    and "eclipse_observed" not in record["reason"]
                ]
                unobserved = [
                    record
                    for record in group
                    if record["kind"] == kind and "eclipse_observed" in record["reason"]
                ]
                widened = [
                    record
                    for record in accepted
                    if abs(
                        float(record["half_phase_window"])
                        - float(record["adopted_half_phase_window"])
                    )
                    > 1e-12
                ]
                if accepted:
                    measured = np.array(
                        [float(record["time_d"]) for record in accepted]
                    )
                    mask = np.array([item == kind for item in kinds])
                    versus_conj = _match(measured, conj_all[mask], period)
                    versus_pure = _match(measured, pure_all[mask], period)
                    median_conj = f"{float(np.median(versus_conj)):.2f}"
                    std_conj = f"{float(np.std(versus_conj)):.2f}"
                    median_pure = f"{float(np.median(versus_pure)):.2f}"
                    std_pure = f"{float(np.std(versus_pure)):.2f}"
                else:
                    median_conj = std_conj = median_pure = std_pure = ""
                offset_rows.append(
                    {
                        "physics_id": str(meta["physics_id"]),
                        "recipe_id": next(
                            row["recipe_id"]
                            for row in rows
                            if Path(row["path"]).name == name
                        ),
                        "cadence_s": group[0]["cadence_s"],
                        "method": method,
                        "kind": kind,
                        "n_accepted": str(len(accepted)),
                        "n_fit_fail": str(len(fit_fail)),
                        "n_unobserved": str(len(unobserved)),
                        "n_widened": str(len(widened)),
                        "median_minus_conj_s": median_conj,
                        "std_minus_conj_s": std_conj,
                        "median_minus_pure_s": median_pure,
                        "std_minus_pure_s": std_pure,
                    }
                )
    offset_path = OUT / "offsets.csv"
    with offset_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(offset_rows[0]))
        writer.writeheader()
        writer.writerows(offset_rows)
    logger.info("wrote %s (%d rows)", offset_path, len(offset_rows))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
