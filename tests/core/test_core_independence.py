"""The core and the selection layer stay free of survey packages.

``mikulasek_eb`` needs numpy, scipy and astropy. ``mikulasek_eb.selection`` and
``mikulasek_eb.compare`` need the same. Nothing there may import ``volightcurve``,
``lightkurve`` or the TESS adapter. The adapter imports the core, never the reverse.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

CORE_SOURCE = Path(__file__).resolve().parents[2] / "src" / "mikulasek_eb"
BLOCKED = ("volightcurve", "lightkurve", "lc_discovery", "matplotlib")
SURVEY_IMPORT = re.compile(
    r"^\s*(?:import|from)\s+(volightcurve|lightkurve|lc_discovery|mikulasek_eb_tess)\b",
    re.M,
)


def test_core_imports_with_survey_packages_blocked():
    """A fresh interpreter imports the three survey-neutral packages without them."""
    blocking = "".join(f"sys.modules[{name!r}] = None\n" for name in BLOCKED)
    code = (
        "import sys\n"
        + blocking
        + "import mikulasek_eb\n"
        + "import mikulasek_eb.selection\n"
        + "import mikulasek_eb.compare\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr


def test_no_core_source_file_imports_a_survey_package():
    """No file under ``src/mikulasek_eb`` imports a reader or the adapter."""
    offenders = []
    for path in sorted(CORE_SOURCE.rglob("*.py")):
        for match in SURVEY_IMPORT.finditer(path.read_text()):
            offenders.append(
                f"{path.relative_to(CORE_SOURCE)}: {match.group(0).strip()}"
            )
    assert not offenders, offenders
