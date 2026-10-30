"""Typed exceptions, so a host application can tell a user error from a numerical failure.

All derive from ``MikulasekError``. Those that replace earlier built-in exceptions
also keep that base (``ValueError`` or ``RuntimeError``), so existing handlers still work.
"""

from __future__ import annotations


class MikulasekError(Exception):
    """Base class of every error raised on purpose by this package."""


class ConfigError(MikulasekError, ValueError):
    """The configuration is incomplete, inconsistent, unsupported or has unknown keys."""


class InvalidLightCurve(MikulasekError, ValueError):
    """The input light curve is structurally invalid or is not an eclipsing binary."""


class FitDidNotConverge(MikulasekError, RuntimeError):
    """The Levenberg-Marquardt solver did not report convergence."""


class NoValidSolution(MikulasekError, RuntimeError):
    """The search finished without any acceptable (finite, positive-width) solution."""


class FitCancelled(MikulasekError, RuntimeError):
    """A progress hook asked to stop the fit (raise this from the hook to cancel)."""
