"""Copy every NSV 807 secondary timing file with the times moved back.

Each copy sits next to its source. The file name gains ``_minus_half_period``
before ``.dat``. The first column of every data line is the original JD minus
half the file period. The file period is the one stored in the chosen
solutions, and it is the same in every sector. Comment lines are kept, and one
comment is added to say what was done. The source files are not changed.

Run from ``mikulasek_eb``::

    python examples/nsv807_secondary_half_period.py
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "nsv807"
RUNS = ("s30_protocol", "s30_600s_protocol", "s4_protocol", "s97_protocol")
SUFFIX = "_minus_half_period"


def file_period() -> float:
    """Read the file period and require every sector to carry the same value.

    Returns:
        float: Period in days.

    Raises:
        RuntimeError: If the saved file periods are not all equal.
    """
    found = []
    for run in RUNS:
        path = next((RESULTS / run / "chosen").glob("*.json"))
        saved = json.loads(path.read_text())
        found.append(float(saved["provenance"]["period_input_days"]))
    if any(period != found[0] for period in found[1:]):
        raise RuntimeError(f"the file periods differ: {found}")
    return found[0]


def sources() -> list[Path]:
    """Secondary timing files that are not already shifted copies.

    Returns:
        list[Path]: Paths, in sorted order.

    Raises:
        RuntimeError: If no secondary file is found.
    """
    found = [
        path
        for path in sorted(RESULTS.rglob("*secondary*.dat"))
        if SUFFIX not in path.stem
    ]
    if not found:
        raise RuntimeError(f"no secondary files in {RESULTS}")
    return found


def shift_file(path: Path, half: float) -> Path:
    """Write one copy with the times reduced by ``half``.

    Args:
        path (Path): Source ``.dat`` file.
        half (float): Half the file period, days.

    Returns:
        Path: The copy.

    Raises:
        RuntimeError: If a data line does not start with a JD and two further
        columns.
    """
    lines = [
        f"# Times in this copy are the source times minus {half:.7f} d, half the file period."
    ]
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            lines.append(line)
            continue
        parts = line.split()
        if len(parts) != 3:
            raise RuntimeError(f"{path}: {line!r}")
        try:
            jd = float(parts[0])
        except ValueError as ex:
            raise RuntimeError(f"{path}: {line!r}") from ex
        lines.append(f"{jd - half:.8f}  {parts[1]}  {parts[2]}")
    dest = path.with_name(f"{path.stem}{SUFFIX}.dat")
    dest.write_text("\n".join(lines) + "\n")
    return dest


def main() -> None:
    """Write a shifted copy of every secondary file."""
    half = file_period() / 2.0
    written = [shift_file(path, half) for path in sources()]
    logger.info("half period %.7f d, wrote %d copies", half, len(written))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    main()
