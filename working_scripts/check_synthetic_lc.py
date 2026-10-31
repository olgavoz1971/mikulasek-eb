"""Check synthetic light curves against the metadata requirements.

Run::

    python working_scripts/check_synthetic_lc.py file1.vot [file2.vot ...]

The pure high-resolution curve named by ``source_file`` is read from the same
folder when it exists, which enables the independent check of the smeared truth.
The exit code is 1 if any check failed.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from mikulasek_eb.synthetic_check import (
    FAIL,
    check_synthetic,
    load_synthetic,
)

logger = logging.getLogger(__name__)

MARK = {"pass": "PASS", "warn": "WARN", "fail": "FAIL", "info": "info", "skip": "skip"}


def main(paths: list[str]) -> int:
    """Check each file and write a report to standard output.

    Args:
        paths (list[str]): Synthetic VOTables.

    Returns:
        int: 1 if any check failed, otherwise 0.
    """
    failed = False
    for name in paths:
        path = Path(name)
        data = load_synthetic(path)
        pure = None
        source = str(data.meta.get("source_file", "")).strip()
        if source and (path.parent / source).is_file():
            pure = load_synthetic(path.parent / source, with_errors=False)
            if pure.problems:
                logger.warning("pure curve %s: %s", source, pure.problems)
                pure = None
        findings = check_synthetic(data, pure)
        sys.stdout.write(f"\n{path.name}\n")
        for f in findings:
            sys.stdout.write(f"  [{MARK[f.status]}] {f.name}: {f.message}\n")
        counts = {s: sum(f.status == s for f in findings) for s in MARK}
        sys.stdout.write(
            f"  summary: {counts['pass']} pass, {counts['warn']} warn, "
            f"{counts['fail']} fail, {counts['skip']} skipped\n"
        )
        failed = failed or any(f.status == FAIL for f in findings)
    return 1 if failed else 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    raise SystemExit(main(sys.argv[1:]))
