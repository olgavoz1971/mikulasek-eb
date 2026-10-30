"""``ModelSolution``: the light, JSON-serialisable result of a fit.

It holds numbers only (parameters, errors, ephemeris, statistics, provenance), never
arrays, so it can be stored in the browser, a database or a cache, sent over REST, saved
next to results, and later used to time eclipses or to seed a fit of another light curve.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path

import numpy as np

from .config import Config, ModelSpec
from .errors import ConfigError
from .model import EBModel

SOLUTION_SCHEMA_VERSION = 1

# Keys of the plain-data form, in order.
_KEYS = (
    "schema_version",
    "config",
    "params",
    "errors",
    "ephemeris",
    "stats",
    "exposure",
    "provenance",
)


def plain(obj):
    """Convert numpy scalars and containers to JSON-compatible Python objects.

    Args:
        obj: A number, string, bool, None, numpy scalar, or a dict, list or tuple of those.

    Returns:
        The same structure with Python ``float``, ``int``, ``bool`` and ``list``.
    """
    if isinstance(obj, dict):
        return {str(k): plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [plain(v) for v in obj]
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    return obj


@dataclass(frozen=True)
class ModelSolution:
    """Result of one fit, as plain data.

    Attributes:
        config (Config): Full configuration that produced it.
        params (dict): Template parameters (A in mag, D in cycles, H and S1 in mag, ...).
        errors (dict): 1-sigma errors of the fitted template and ephemeris parameters.
        ephemeris (dict): M0 (MJD), P0 (days), Pdot (dimensionless).
        stats (dict): chi2, dof, rchi2, rms_mag, noise_floor_mag, rms_over_noise, n_points,
            n_input, rank_deficient, two_component.
        exposure (dict | None): ``{"days", "nsub"}`` actually used, or None.
        provenance (dict): data_id, epoch_jd, period_input_days, package_version.
        schema_version (int): Version of the serialised layout.
    """

    config: Config
    params: dict
    errors: dict
    ephemeris: dict
    stats: dict
    exposure: dict | None
    provenance: dict
    schema_version: int = SOLUTION_SCHEMA_VERSION

    @property
    def spec(self) -> ModelSpec:
        """Model structure of this solution."""
        return self.config.model

    @property
    def p(self) -> dict:
        """Template parameters (read access compatible with ``FitResult.p``)."""
        return dict(self.params)

    @property
    def ephem(self) -> dict:
        """Ephemeris (read access compatible with ``FitResult.ephem``)."""
        return dict(self.ephemeris)

    def build_model(self) -> EBModel:
        """Rebuild the ``EBModel`` that was fitted, including exposure integration.

        Returns:
            EBModel: Model with the stored profile, harmonics, O'Connell, K and exposure.
        """
        spec = self.config.model
        extra = {}
        if self.exposure is not None:
            extra = dict(
                exposure_days=self.exposure["days"],
                nsub=self.exposure["nsub"],
                period_days=self.provenance["period_input_days"],
            )
        return EBModel(
            profile=spec.profile,
            harmonics=spec.harmonics,
            ocon=spec.ocon,
            kterm=spec.kterm,
            **extra,
        )

    def to_dict(self) -> dict:
        """Plain-data form.

        Returns:
            dict: JSON-compatible dict with the keys of ``_KEYS``.
        """
        return plain(
            dict(
                schema_version=self.schema_version,
                config=self.config.to_dict(),
                params=self.params,
                errors=self.errors,
                ephemeris=self.ephemeris,
                stats=self.stats,
                exposure=self.exposure,
                provenance=self.provenance,
            )
        )

    @classmethod
    def from_dict(cls, data: dict) -> ModelSolution:
        """Rebuild a solution from its plain-data form.

        Args:
            data (dict): Output of ``to_dict`` or parsed JSON.

        Returns:
            ModelSolution: The solution.

        Raises:
            ConfigError: On unknown or missing keys or a wrong schema version.
        """
        if not isinstance(data, dict):
            raise ConfigError("A solution must be a JSON object.")
        if set(data) != set(_KEYS):
            raise ConfigError(
                f"Solution keys must be {sorted(_KEYS)}, got {sorted(data)}."
            )
        if data["schema_version"] != SOLUTION_SCHEMA_VERSION:
            raise ConfigError(
                f"Unsupported solution schema_version {data['schema_version']}."
            )
        return cls(
            config=Config.from_dict(data["config"]),
            params=dict(data["params"]),
            errors=dict(data["errors"]),
            ephemeris=dict(data["ephemeris"]),
            stats=dict(data["stats"]),
            exposure=None if data["exposure"] is None else dict(data["exposure"]),
            provenance=dict(data["provenance"]),
        )

    def to_json(self, indent: int | None = 2) -> str:
        """Serialise to JSON text.

        Args:
            indent (int | None): JSON indent; None gives a compact string.

        Returns:
            str: JSON document.
        """
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_json(cls, text: str) -> ModelSolution:
        """Parse a JSON document.

        Args:
            text (str): JSON text.

        Returns:
            ModelSolution: The solution.

        Raises:
            ConfigError: If the text is not valid JSON or the content is invalid.
        """
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ConfigError(f"Solution is not valid JSON: {exc}") from exc
        return cls.from_dict(parsed)

    def save(self, path: str | Path) -> None:
        """Write the solution as JSON.

        Args:
            path (str | Path): Output file.
        """
        Path(path).write_text(self.to_json() + "\n", encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> ModelSolution:
        """Read a solution file.

        Args:
            path (str | Path): JSON file.

        Returns:
            ModelSolution: The solution.
        """
        return cls.from_json(Path(path).read_text(encoding="utf-8"))

    def solution_hash(self) -> str:
        """Deterministic identity of the solution (config, numbers and provenance).

        Returns:
            str: First 12 hex digits of the SHA-256 of the canonical JSON.
        """
        canon = json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canon.encode("utf-8")).hexdigest()[:12]


def package_version() -> str:
    """Installed version of this package, for provenance.

    Returns:
        str: Version string.

    Raises:
        importlib.metadata.PackageNotFoundError: If the package is not installed.
    """
    return metadata.version("mikulasek-eb")
