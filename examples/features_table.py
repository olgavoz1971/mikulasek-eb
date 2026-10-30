"""Measure the features of every light curve in ``data/`` and write ``results/features_summary.md``.

Run::

    python examples/features_table.py
"""

from __future__ import annotations

import logging
from pathlib import Path

from mikulasek_eb import extract_features
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable

logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
STARS = {
    "AA And S16 2-min": "lc_discovery_TIC_455790537_SPOC_sector_16",
    "AA And S16 30-min": "lc_discovery_TIC_455790537_TESS_SPOC_sector_16",
    "XY Boo S50 10-min": "lc_discovery_TIC_347405567_TESS_SPOC_sector_50",
    "EK Com S49": "EK_Com_sector_49",
    "V477 Lyr S26": "v477_Lyr_sector_26",
    "AV Del S81": "AV_Del_sector_81",
    "AR Aur S19": "AR_Aur_sector_19",
}


def fmt(value, digits: int = 3) -> str:
    """Format a number, None or bool for the table.

    Args:
        value: Value to show.
        digits (int): Significant digits for floats.

    Returns:
        str: Text.
    """
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def main() -> None:
    """Write the summary table and the notes of each star."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    rows, notes = [], []
    for label, stem in STARS.items():
        payload = load_tess_votable(ROOT / "data" / f"{stem}.vot")
        t = absolute_jd_to_mjd(payload["t"])
        f = extract_features(
            t, payload["mag"], payload["sig"], payload["period"], payload["epoch_jd"]
        )
        g, sh, v, e, r = f.geometry, f.shape, f.variation, f.ephemeris, f.residuals
        d = e["drift_between_halves"]
        shift = lambda k: (
            "-"
            if d[k] is None
            else f"{d[k]['shift_s']:.0f} +- {d[k]['shift_err_s']:.0f}"
        )
        md = v["maxima_difference"]
        rows.append(
            [
                label,
                fmt(g["A1"]),
                fmt(g["A2"]),
                fmt(g["fwhm1"]),
                fmt(g["fwhm2"]),
                fmt(g["points_per_fwhm"]),
                fmt(g["exposure_over_fwhm"]),
                fmt(sh["bottom_ratio"]),
                fmt(sh["wing_ratio"]),
                fmt(sh["gamma_estimate"]),
                fmt(v["eclipse_coverage"]),
                fmt(v["continuous_variation"]),
                f"{1e3 * md['value']:.1f} ({fmt(v['maxima_significance'], 2)} sigma)",
                shift("primary"),
                shift("secondary"),
                fmt(r["fold_scatter_ratio"]),
                fmt(r["trend_order_hint"]),
            ]
        )
        notes.append(
            f"- **{label}** (hash {f.features_hash()}): "
            + ("; ".join(f.notes) or "none")
        )
    head = [
        "star",
        "A1",
        "A2",
        "fwhm1",
        "fwhm2",
        "pts/fwhm",
        "expo/fwhm",
        "W90/W50",
        "W10/W50",
        "Gamma est.",
        "coverage",
        "contact",
        "maxima diff (mmag)",
        "drift prim. (s)",
        "drift sec. (s)",
        "fold scatter/noise",
        "trend order",
    ]
    lines = ["# Features of the light curves in data/", ""]
    lines.append("| " + " | ".join(head) + " |")
    lines.append("|" + "---|" * len(head))
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    lines += ["", "Notes (features that could not be measured):", ""] + notes
    out = ROOT / "results" / "features_summary.md"
    out.write_text("\n".join(lines) + "\n")
    logger.info("wrote %s", out)


if __name__ == "__main__":
    main()
