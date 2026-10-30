"""Rebuild a QLP VOTable whose flux zero point is declared dimensionless.

The QLP export declares ``zeroPointFlux = 1`` with the unit ``---`` while the
flux column is in ``s-1``. The researcher confirmed that the zero point is
1 s-1, so the light curve is re-assembled through the public volightcurve API
with that explicit calibration. Everything else is copied from the file.

Run::

    python working_scripts/fix_qlp_zero_point_unit.py data/NSV_807_sector_04.vot
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from volightcurve import assemble_volightcurve, read_lightcurve, write_lightcurve

logger = logging.getLogger(__name__)


def fix(path: Path) -> Path:
    """Re-assemble a light curve with the flux zero point declared in s-1.

    Args:
        path (Path): VOTable to fix; it is overwritten.

    Returns:
        Path: The written file.

    Raises:
        ValueError: If the declared calibration is not the expected QLP one.
    """
    volc = read_lightcurve(path, filename=path.name)
    hub = next(iter(volc.photdms.values()))
    pc = hub.photcal
    if float(pc.zp_flux.value) != 1.0 or not pc.zp_flux.unit.is_equivalent(""):
        raise ValueError(f"Unexpected declared flux zero point {pc.zp_flux!r}")
    logger.info("Declared: zp_flux=%r, zp_mag=%r", pc.zp_flux, pc.zp_mag)
    meta = volc.table.meta
    fixed = assemble_volightcurve(
        volc.table,
        timeorigin=volc.timesys.timeorigin,
        period=meta.get("period"),
        epoch=meta.get("epoch"),
        filter_id=hub.filter.filter_id,
        filter_name=hub.filter.name,
        zp_flux=float(pc.zp_flux.value),
        zp_flux_unit="1 / s",
        zp_mag=float(pc.zp_mag.value),
        zp_mag_unit="mag",
        mag_sys=pc.mag_sys,
        effective_wavelength=float(hub.filter.spectral_location.value),
        effective_wavelength_unit="m",
        free_comments=[
            "Flux zero point unit set to 1/s on the researcher's confirmation "
            "(the QLP export declared it dimensionless)."
        ],
    )
    # assemble_volightcurve has no time scale argument; carry over what the file declared.
    fixed.timesys.timescale = volc.timesys.timescale
    fixed.timesys.refposition = volc.timesys.refposition
    write_lightcurve(fixed, "votable_binary", destination=path)
    return path


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    print(fix(Path(sys.argv[1])))
