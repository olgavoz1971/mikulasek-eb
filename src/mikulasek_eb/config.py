"""Layered, serialisable configuration: model spec, data setup, fit policy and timing recipe.

The four layers have different lifetimes (see ``docs/mikulasek_implementation_notes.md``):

* ``ModelSpec``: what the star needs (profile, harmonics, ...). Transferable between
  light curves of the same star.
* ``DataSetup``: properties of one light curve (baseline, exposure). Not transferable.
* ``FitPolicy``: how the search runs and what is accepted. Global defaults.
* ``TimingRecipe``: how times of minima are measured. Kept separately so that the same
  recipe can be applied to several light curves for homogeneity.

``Config`` bundles them. It is plain data: it round-trips through JSON (the canonical
form for a host application, a notebook or an agent) and has a deterministic hash for
cache keys and provenance. ``FitConfig`` is the flat, resolved view that the staged
pipeline reads.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
from dataclasses import dataclass, field, fields
from pathlib import Path

from .errors import ConfigError

SCHEMA_VERSION = 1

# Fraction-of-phase limit for any timing window: whole phase around one event.
MAX_WINDOW = 0.25

# Defaults of the timing recipe (also exported by ``mikulasek_eb.tom`` under these names).
HALF_ORBIT_WINDOW = 0.25
LOCAL_WINDOWS = (0.15, 0.10, 0.06)
TIMING_BASELINE_ORDER = 1
MAX_SIGMA_FACTOR = 5.0

# Fields that are tuples in Python and lists in JSON.
_TUPLE_FIELDS = {"harmonics", "gamma_starts", "two_comp_eclipses", "local_windows"}


@dataclass(frozen=True)
class ModelSpec:
    """What the star needs: the structure of the phenomenological model.

    Attributes:
        profile (int): Mikulasek Table 2 profile number (1 to 7).
        harmonics (tuple[int, ...]): Orders k of the cosine terms H_k.
        ocon (bool): Include the O'Connell sine term S1.
        kterm (bool): Include the Eq. (14) quartic term K (profile 1 only).
    """

    profile: int
    harmonics: tuple[int, ...]
    ocon: bool
    kterm: bool = False


@dataclass(frozen=True)
class DataSetup:
    """Properties of one light curve.

    Attributes:
        baseline_order (int): Per-segment polynomial degree of the instrumental term.
        gap_split_days (float | None): Split baselines at gaps longer than this
            (days); None means one baseline for the whole series.
        exposure_days (float | None): Exposure time (days) for exposure integration of
            the model (extension, not in the paper); None means instantaneous samples.
        exposure_nsub (int | None): Sub-exposure samples; None chooses from the exposure
            and the folded FWHM (``exposure.suggest_nsub``).
    """

    baseline_order: int
    gap_split_days: float | None
    exposure_days: float | None = None
    exposure_nsub: int | None = None


@dataclass(frozen=True)
class FitPolicy:
    """How the staged search runs and what it accepts.

    The two-component and eclipse-weight settings are extensions, not in the paper.
    Values are our heuristics unless stated otherwise.

    Attributes:
        gamma_starts (tuple): Multi-start Gamma values (includes Gamma < 1, as in Table 1).
        separate_shape_chi2_frac (float): Unlock D2/G2 if chi2 drops by more than this fraction.
        clip_nsig (float): Out-of-eclipse MAD clip threshold (sigma).
        allow_two_component (bool): Try the second profile component after a poor fit.
        rms_noise_two_comp_threshold (float): Try the two-component fit above this rms/noise.
        eclipse_weight (float): Extra least-squares weight inside both eclipse windows
            (``sigma_eff = sigma / weight``); 1.0 means catalogue sigma only.
        eclipse_phase_half (float): Phase half-width around 0 and 0.5 for that weight (cycles).
        min_points_per_fwhm (float): Warn below this many samples per eclipse FWHM.
        primary_a2_over_a1_min_ratio (float): Warn that P may be P/2 if A1 < ratio * A2.
        two_comp_min_rms_gain_frac (float): Acceptance: fractional rms improvement.
        two_comp_min_agreeing_starts (int): Acceptance: starts agreeing with the best chi2.
        two_comp_max_rel_err (float): Acceptance: largest relative parameter error.
        two_comp_eclipses (tuple[int, ...]): Eclipses (1, 2) that get the second component.
        two_comp_timing_agreement_sigma_max (float): Acceptance: median |dT|/combined sigma.
    """

    gamma_starts: tuple = (0.5, 1, 2, 4, 8)
    separate_shape_chi2_frac: float = 0.01
    clip_nsig: float = 5.0
    allow_two_component: bool = True
    rms_noise_two_comp_threshold: float = 1.5
    eclipse_weight: float = 1.0
    eclipse_phase_half: float = 0.12
    min_points_per_fwhm: float = 8.0
    primary_a2_over_a1_min_ratio: float = 1.2
    two_comp_min_rms_gain_frac: float = 0.25
    two_comp_min_agreeing_starts: int = 2
    two_comp_max_rel_err: float = 0.5
    two_comp_eclipses: tuple = (1,)
    two_comp_timing_agreement_sigma_max: float = 1.0


@dataclass(frozen=True)
class TimingRecipe:
    """How times of minima are measured (kept separately for homogeneity).

    Attributes:
        window (float): Phase half-width of the fit around each eclipse (cycles, at most 0.25).
        baseline_order (int): Per-event polynomial degree.
        max_sigma_factor (float): Drop events whose sigma exceeds this times the median.
        require_contact_walls (bool): Require points on both walls (90 to 10 percent band).
        clip_outliers (bool): Apply the sigma cut and the 4-sigma O-C clip.
        local_windows (tuple[float, ...]): Windows tried for local primary/secondary-only files.
        template (str): Which template the timings use. Only ``"own"`` (the solution of the
            same light curve) is implemented.
    """

    window: float = HALF_ORBIT_WINDOW
    baseline_order: int = TIMING_BASELINE_ORDER
    max_sigma_factor: float = MAX_SIGMA_FACTOR
    require_contact_walls: bool = True
    clip_outliers: bool = True
    local_windows: tuple = LOCAL_WINDOWS
    template: str = "own"


@dataclass(frozen=True)
class FitConfig:
    """Flat, resolved view of the model, data and policy layers read by the pipeline.

    Use ``Config`` as the public, layered form; ``Config.to_fit_config()`` builds this.
    Field names are the union of ``ModelSpec``, ``DataSetup`` and ``FitPolicy``.
    """

    profile: int
    harmonics: tuple[int, ...]
    ocon: bool
    baseline_order: int
    gap_split_days: float | None
    kterm: bool = False
    exposure_days: float | None = None
    exposure_nsub: int | None = None
    gamma_starts: tuple = (0.5, 1, 2, 4, 8)
    separate_shape_chi2_frac: float = 0.01
    clip_nsig: float = 5.0
    allow_two_component: bool = True
    rms_noise_two_comp_threshold: float = 1.5
    eclipse_weight: float = 1.0
    eclipse_phase_half: float = 0.12
    min_points_per_fwhm: float = 8.0
    primary_a2_over_a1_min_ratio: float = 1.2
    two_comp_min_rms_gain_frac: float = 0.25
    two_comp_min_agreeing_starts: int = 2
    two_comp_max_rel_err: float = 0.5
    two_comp_eclipses: tuple = (1,)
    two_comp_timing_agreement_sigma_max: float = 1.0


def _layer_from_dict(cls, data: dict, layer: str):
    """Build one layer from a dict, rejecting unknown keys and restoring tuples.

    Args:
        cls: Layer dataclass.
        data (dict): Parsed JSON object.
        layer (str): Layer name for the messages.

    Returns:
        An instance of ``cls``.

    Raises:
        ConfigError: On unknown or missing keys.
    """
    if not isinstance(data, dict):
        raise ConfigError(
            f"Layer {layer!r} must be an object, got {type(data).__name__}."
        )
    known = {f.name for f in fields(cls)}
    unknown = set(data) - known
    if unknown:
        raise ConfigError(f"Unknown keys in {layer!r}: {sorted(unknown)}.")
    kwargs = {
        k: (tuple(v) if k in _TUPLE_FIELDS and v is not None else v)
        for k, v in data.items()
    }
    try:
        return cls(**kwargs)
    except TypeError as exc:
        raise ConfigError(f"Invalid {layer!r} layer: {exc}") from exc


@dataclass(frozen=True)
class Config:
    """Complete configuration of one analysis: four layers plus free-text notes.

    Attributes:
        model (ModelSpec): Model structure.
        data (DataSetup): Light-curve properties.
        policy (FitPolicy): Search and acceptance rules.
        timing (TimingRecipe): Timing method.
        notes (str): Free text (not part of the hash).
        schema_version (int): Version of the serialised layout.
    """

    model: ModelSpec
    data: DataSetup
    policy: FitPolicy = field(default_factory=FitPolicy)
    timing: TimingRecipe = field(default_factory=TimingRecipe)
    notes: str = ""
    schema_version: int = SCHEMA_VERSION

    def validate(self) -> Config:
        """Check the configuration for impossible or unsupported combinations.

        Returns:
            Config: ``self``, so calls can be chained.

        Raises:
            ConfigError: With a message naming the first problem found.
        """
        m, d, p, t = self.model, self.data, self.policy, self.timing
        if self.schema_version != SCHEMA_VERSION:
            raise ConfigError(
                f"Unsupported schema_version {self.schema_version} (this code reads {SCHEMA_VERSION})."
            )
        if m.profile not in range(1, 8):
            raise ConfigError(f"profile must be 1 to 7, got {m.profile!r}.")
        if m.kterm and m.profile != 1:
            raise ConfigError("kterm (Eq. 14) is defined for profile 1 only.")
        if any((not isinstance(k, int)) or k < 1 for k in m.harmonics):
            raise ConfigError(
                f"harmonics must be positive integers, got {m.harmonics!r}."
            )
        if d.baseline_order < 0:
            raise ConfigError("baseline_order must be >= 0.")
        if d.gap_split_days is not None and not d.gap_split_days > 0:
            raise ConfigError("gap_split_days must be positive or None.")
        if d.exposure_days is not None and not d.exposure_days > 0:
            raise ConfigError("exposure_days must be positive or None.")
        if d.exposure_nsub is not None:
            if d.exposure_days is None:
                raise ConfigError("exposure_nsub is set but exposure_days is None.")
            if int(d.exposure_nsub) != d.exposure_nsub or d.exposure_nsub < 1:
                raise ConfigError("exposure_nsub must be a positive integer or None.")
        if not p.gamma_starts or any(not g > 0 for g in p.gamma_starts):
            raise ConfigError(
                "gamma_starts must be a non-empty list of positive numbers."
            )
        for name in (
            "clip_nsig",
            "eclipse_weight",
            "min_points_per_fwhm",
            "primary_a2_over_a1_min_ratio",
            "two_comp_max_rel_err",
            "two_comp_timing_agreement_sigma_max",
        ):
            if not getattr(p, name) > 0:
                raise ConfigError(f"{name} must be positive.")
        if not 0 < p.eclipse_phase_half <= MAX_WINDOW:
            raise ConfigError(f"eclipse_phase_half must be in (0, {MAX_WINDOW}].")
        if not 0 < p.two_comp_min_rms_gain_frac < 1:
            raise ConfigError("two_comp_min_rms_gain_frac must be in (0, 1).")
        if p.two_comp_min_agreeing_starts < 1:
            raise ConfigError("two_comp_min_agreeing_starts must be >= 1.")
        if not p.two_comp_eclipses or not set(p.two_comp_eclipses) <= {1, 2}:
            raise ConfigError("two_comp_eclipses must be a non-empty subset of {1, 2}.")
        if not 0 < t.window <= MAX_WINDOW:
            raise ConfigError(f"timing window must be in (0, {MAX_WINDOW}].")
        if any(not 0 < w <= MAX_WINDOW for w in t.local_windows):
            raise ConfigError(f"timing local_windows must be in (0, {MAX_WINDOW}].")
        if t.baseline_order < 0 or not t.max_sigma_factor > 0:
            raise ConfigError(
                "timing baseline_order must be >= 0 and max_sigma_factor > 0."
            )
        if t.template != "own":
            raise ConfigError(
                f"timing template {t.template!r} is not implemented; only 'own' is."
            )
        return self

    def to_dict(self) -> dict:
        """Plain-data form (lists instead of tuples), ready for JSON.

        Returns:
            dict: ``{"schema_version", "model", "data", "policy", "timing", "notes"}``.
        """

        def plain(obj):
            out = dataclasses.asdict(obj)
            return {k: (list(v) if isinstance(v, tuple) else v) for k, v in out.items()}

        return dict(
            schema_version=self.schema_version,
            model=plain(self.model),
            data=plain(self.data),
            policy=plain(self.policy),
            timing=plain(self.timing),
            notes=self.notes,
        )

    @classmethod
    def from_dict(cls, data: dict) -> Config:
        """Build and validate a configuration from its plain-data form.

        Args:
            data (dict): Output of ``to_dict`` or an equivalent parsed JSON object.
                ``policy`` and ``timing`` may be omitted (defaults); ``model`` and
                ``data`` are required.

        Returns:
            Config: Validated configuration.

        Raises:
            ConfigError: On unknown or missing keys, or invalid values.
        """
        if not isinstance(data, dict):
            raise ConfigError("Configuration must be a JSON object.")
        unknown = set(data) - {
            "schema_version",
            "model",
            "data",
            "policy",
            "timing",
            "notes",
        }
        if unknown:
            raise ConfigError(f"Unknown top-level keys: {sorted(unknown)}.")
        for required in ("model", "data"):
            if required not in data:
                raise ConfigError(f"Missing required layer {required!r}.")
        cfg = cls(
            model=_layer_from_dict(ModelSpec, data["model"], "model"),
            data=_layer_from_dict(DataSetup, data["data"], "data"),
            policy=_layer_from_dict(FitPolicy, data.get("policy", {}), "policy"),
            timing=_layer_from_dict(TimingRecipe, data.get("timing", {}), "timing"),
            notes=data.get("notes", ""),
            schema_version=data.get("schema_version", SCHEMA_VERSION),
        )
        return cfg.validate()

    def to_json(self, indent: int | None = 2) -> str:
        """Serialise to JSON text.

        Args:
            indent (int | None): JSON indent; None gives a compact string.

        Returns:
            str: JSON document.
        """
        return json.dumps(self.to_dict(), indent=indent)

    @classmethod
    def from_json(cls, text: str) -> Config:
        """Parse and validate a JSON document.

        Args:
            text (str): JSON text.

        Returns:
            Config: Validated configuration.

        Raises:
            ConfigError: If the text is not valid JSON or the content is invalid.
        """
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ConfigError(f"Configuration is not valid JSON: {exc}") from exc
        return cls.from_dict(parsed)

    def save(self, path: str | Path) -> None:
        """Write the configuration as JSON.

        Args:
            path (str | Path): Output file.
        """
        Path(path).write_text(self.to_json() + "\n", encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> Config:
        """Read and validate a configuration file.

        Args:
            path (str | Path): JSON file.

        Returns:
            Config: Validated configuration.
        """
        return cls.from_json(Path(path).read_text(encoding="utf-8"))

    def config_hash(self) -> str:
        """Deterministic identity of the configuration (notes excluded).

        Returns:
            str: First 12 hex digits of the SHA-256 of the canonical JSON.
        """
        body = self.to_dict()
        body.pop("notes")
        canon = json.dumps(body, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canon.encode("utf-8")).hexdigest()[:12]

    def to_fit_config(self) -> FitConfig:
        """Flat view of the model, data and policy layers for the pipeline.

        Returns:
            FitConfig: Resolved flat configuration.
        """
        flat = {}
        for layer in (self.model, self.data, self.policy):
            flat.update({f.name: getattr(layer, f.name) for f in fields(layer)})
        return FitConfig(**flat)

    @classmethod
    def from_fit_config(
        cls, flat: FitConfig, timing: TimingRecipe | None = None, notes: str = ""
    ) -> Config:
        """Layered view of a flat configuration.

        Args:
            flat (FitConfig): Flat configuration.
            timing (TimingRecipe | None): Timing recipe; None means the default recipe.
            notes (str): Free text.

        Returns:
            Config: Validated layered configuration.
        """

        def pick(layer):
            return layer(**{f.name: getattr(flat, f.name) for f in fields(layer)})

        return cls(
            model=pick(ModelSpec),
            data=pick(DataSetup),
            policy=pick(FitPolicy),
            timing=timing if timing is not None else TimingRecipe(),
            notes=notes,
        ).validate()
