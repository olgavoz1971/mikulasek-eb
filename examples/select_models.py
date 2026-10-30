"""Run automatic model selection on light curves and save the evidence.

For each star it writes, into ``results/selection/``:
``<stem>_selection.json`` (candidates, BIC table, reasons, features hash),
``<stem>_solution.json`` (the chosen fit) and one row of ``summary.md``.

Run::

    python examples/select_models.py                 # all stars in features_table.STARS
    python examples/select_models.py "EK Com S49"    # by label
    python examples/select_models.py path/to/lc.vot  # any VOTable (new stars)
"""

from __future__ import annotations

import json
import logging
import sys
import time
from pathlib import Path

from features_table import ROOT, STARS

from mikulasek_eb import select_model
from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb.io import load_tess_votable

logger = logging.getLogger(__name__)
OUT = ROOT / "results" / "selection"


def resolve(arg: str) -> tuple[str, Path]:
    """Map a command-line argument to a label and a VOTable path.

    Args:
        arg (str): Star label from ``STARS`` or a path to a ``.vot`` file.

    Returns:
        tuple[str, Path]: Label and file path.

    Raises:
        SystemExit: If the argument is neither a known label nor an existing file.
    """
    if arg in STARS:
        return arg, ROOT / "data" / f"{STARS[arg]}.vot"
    path = Path(arg)
    if path.is_file():
        return path.stem, path
    raise SystemExit(f"Unknown star label or missing file: {arg}")


def run_one(label: str, path: Path) -> str:
    """Select a model for one star and save the results.

    Args:
        label (str): Name used in the summary.
        path (Path): VOTable with the light curve.

    Returns:
        str: One markdown table row.
    """
    payload = load_tess_votable(path)
    t = absolute_jd_to_mjd(payload["t"])
    t0 = time.time()
    result = select_model(
        t,
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        data_id=label,
    )
    OUT.mkdir(parents=True, exist_ok=True)
    stem = path.stem
    (OUT / f"{stem}_selection.json").write_text(json.dumps(result.to_dict(), indent=2))
    best = "-"
    if result.best is not None:
        result.output.solution.save(OUT / f"{stem}_solution.json")
        best = result.best.name
    rows = [r for r in result.table if r.delta_bic is not None]
    second = sorted(r.delta_bic for r in rows)[1] if len(rows) > 1 else None
    logger.info("%s: best %s, ask_human %s", label, best, result.ask_human)
    return (
        f"| {label} | {best} | {'-' if second is None else f'{second:.1f}'} | "
        f"{'yes' if result.ask_human else 'no'} | {result.n_eff:.0f}"
        f"{'' if result.n_eff_measured else ' (not measured)'} | "
        f"{time.time() - t0:.0f} s | {'; '.join(result.reasons) or '-'} |"
    )


def main() -> None:
    """Process the requested stars and write ``summary.md``."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    args = sys.argv[1:] or list(STARS)
    lines = [
        "| Star | Best | dBIC to next | Ask human | n_eff | Time | Reasons |",
        "|---|---|---|---|---|---|---|",
    ]
    for arg in args:
        lines.append(run_one(*resolve(arg)))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")
    logger.info("Wrote %s", OUT / "summary.md")


if __name__ == "__main__":
    main()
