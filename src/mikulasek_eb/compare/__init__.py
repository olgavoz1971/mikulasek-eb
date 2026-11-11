"""Model-free reference measurements, for comparison with the template timing.

Nothing here uses the Mikulasek model. These tools measure a minimum straight from the
points, so a template time can be checked against them.
"""

from .curve_minimum import BISECTOR_FRACTIONS, CurveMinimum, measure_curve_minimum

__all__ = ["BISECTOR_FRACTIONS", "CurveMinimum", "measure_curve_minimum"]
