"""The plotting function saves a file and closes the figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np

from mikulasek_eb.plotting import plot_instrumental_trend


def test_plot_instrumental_trend_writes_the_file(tmp_path: Path):
    """A known baseline is accepted and the picture file is created."""
    t = np.linspace(0.0, 2.0, 50)
    template = 0.2 * np.abs((t % 1.0) - 0.5)
    trend = 0.01 * t
    path = tmp_path / "trend.png"
    plot_instrumental_trend(t, template + trend, template + trend, template, path, "test")
    assert path.is_file()
    assert path.stat().st_size > 0
