"""Validated parameters and paper presets for the compact resonator.

All geometry dimensions are expressed in micrometres.  The values of ``g_c``,
``g_l``, ``g_r``, ``w`` and ``Z0`` in the presets are transcribed from
Geerlings et al., Applied Physics Letters 100, 192601 (2012).  The paper does
not report capacitor-finger or inductor-turn counts; the preset counts are
therefore explicitly documented, deterministic reconstruction assumptions.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
from math import isfinite
from os import PathLike
from pathlib import Path
from types import MappingProxyType
from typing import Any, ClassVar, Mapping
import tomllib


COUNT_INFERENCE_NOTE = (
    "The paper text does not specify capacitor-finger or inductor-turn counts. "
    "Design A uses the counts visible in the Fig. 1 schematic (14 IDC fingers "
    "and 10 half-width spans per inductor side); B/C counts are reconstruction "
    "assumptions constrained by their reported envelopes. Change the counts "
    "when mask data or a target L/C is available."
)


def _positive_number(name: str, value: object) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a real number, got {value!r}")
    if not isfinite(float(value)) or float(value) <= 0.0:
        raise ValueError(f"{name} must be finite and greater than zero")


def _layer_number(name: str, value: object) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer GDS layer number")
    if not 0 <= value <= 65535:
        raise ValueError(f"{name} must be between 0 and 65535")


@dataclass(frozen=True, slots=True)
class ResonatorParameters:
    """Parameters controlling the reconstructed Fig. 1 layout.

    ``w_um`` is the common trace width used by designs A and B.  Design C uses
    the paper's split widths through ``capacitor_w_um=40`` and
    ``inductor_w_um=10``.  A value of ``None`` for either override means
    ``w_um`` is used.

    The default constructor is the design-A reconstruction requested by the
    paper: ``gC=10 um``, ``gL=20 um``, ``gR=10 um``, ``w=5 um``, ``Z0=100
    ohm`` and 300 um nominal resonator width.
    """

    COUNT_ASSUMPTION: ClassVar[str] = COUNT_INFERENCE_NOTE

    design_name: str = "A"
    g_c_um: float = 10.0
    g_l_um: float = 20.0
    g_r_um: float = 10.0
    w_um: float = 5.0
    z0_ohm: float = 100.0
    nominal_width_um: float = 300.0

    # Counts are reconstruction assumptions, not values reported in the paper.
    capacitor_finger_count: int = 14
    # Per side: the generated path has 2*inductor_turns half-width spans plus
    # one full-width bridge, i.e. 21 horizontal spans for the Fig. 1 default.
    inductor_turns: int = 10
    capacitor_w_um: float | None = None
    inductor_w_um: float | None = None

    # CPW and finite layout window.  A None length is resolved from the
    # resonator width and the ground margin by the geometry generator.
    include_feedline: bool = False
    feedline_width_um: float = 10.0
    feedline_gap_um: float = 6.0
    feedline_length_um: float | None = None
    feedline_end_clearance_um: float = 10.0
    ground_bridge_width_um: float = 20.0
    coupling_ground_width_um: float = 20.0
    ground_margin_um: float = 100.0

    # setupEM-compatible driven-verification marker geometry.
    include_ports: bool = False
    port_depth_um: float = 2.0
    metal_layer: int = 1
    metal_datatype: int = 0
    left_port_layer: int = 201
    right_port_layer: int = 202
    port_datatype: int = 0

    # Physical thicknesses used only by the optional build123d STEP exporter.
    metal_thickness_um: float = 0.2
    substrate_thickness_um: float = 430.0

    def __post_init__(self) -> None:
        if not isinstance(self.design_name, str) or not self.design_name.strip():
            raise ValueError("design_name must be a non-empty string")

        for name in (
            "g_c_um",
            "g_l_um",
            "g_r_um",
            "w_um",
            "z0_ohm",
            "nominal_width_um",
            "feedline_width_um",
            "feedline_gap_um",
            "feedline_end_clearance_um",
            "ground_bridge_width_um",
            "coupling_ground_width_um",
            "ground_margin_um",
            "port_depth_um",
            "metal_thickness_um",
            "substrate_thickness_um",
        ):
            _positive_number(name, getattr(self, name))

        for name in ("capacitor_w_um", "inductor_w_um", "feedline_length_um"):
            value = getattr(self, name)
            if value is not None:
                _positive_number(name, value)

        for name in ("capacitor_finger_count", "inductor_turns"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"{name} must be an integer")
            if value < 2:
                raise ValueError(f"{name} must be at least 2")

        for name in ("include_feedline", "include_ports"):
            if not isinstance(getattr(self, name), bool):
                raise TypeError(f"{name} must be true or false")
        if self.include_ports and not self.include_feedline:
            raise ValueError("include_ports=true requires include_feedline=true")

        for name in (
            "metal_layer",
            "metal_datatype",
            "left_port_layer",
            "right_port_layer",
            "port_datatype",
        ):
            _layer_number(name, getattr(self, name))
        if self.include_ports:
            if self.left_port_layer == self.right_port_layer:
                raise ValueError("left_port_layer and right_port_layer must differ")
            if self.metal_layer in {self.left_port_layer, self.right_port_layer}:
                raise ValueError("port marker layers must differ from metal_layer")

        minimum_width = (
            2.0 * self.g_r_um
            + 2.0 * self.capacitor_trace_width_um
            + self.g_c_um
            + 2.0 * self.inductor_trace_width_um
        )
        if self.nominal_width_um <= minimum_width:
            raise ValueError(
                "nominal_width_um is too small for g_r, g_c and the trace "
                f"widths; it must exceed {minimum_width:g} um"
            )

        if self.feedline_length_um is not None:
            required = self.nominal_width_um + 2.0 * self.ground_margin_um
            if self.feedline_length_um < required:
                raise ValueError(
                    "feedline_length_um must be at least nominal_width_um + "
                    f"2*ground_margin_um ({required:g} um)"
                )
        if 2.0 * self.port_depth_um >= self.resolved_feedline_length_um:
            raise ValueError("port_depth_um is too large for the feedline")

    @property
    def capacitor_trace_width_um(self) -> float:
        """Effective capacitor trace width in micrometres."""

        return self.w_um if self.capacitor_w_um is None else self.capacitor_w_um

    @property
    def inductor_trace_width_um(self) -> float:
        """Effective inductor trace width in micrometres."""

        return self.w_um if self.inductor_w_um is None else self.inductor_w_um

    @property
    def resolved_feedline_length_um(self) -> float:
        """Finite CPW/layout-window length used by geometry generation."""

        if self.feedline_length_um is not None:
            return self.feedline_length_um
        return self.nominal_width_um + 2.0 * self.ground_margin_um

    @property
    def inferred_counts(self) -> Mapping[str, int]:
        """Return the two reconstruction counts, explicitly labelled inferred."""

        return MappingProxyType(
            {
                "capacitor_finger_count": self.capacitor_finger_count,
                "inductor_turns": self.inductor_turns,
            }
        )

    def with_updates(self, **changes: object) -> "ResonatorParameters":
        """Return a validated copy with selected values replaced."""

        return replace(self, **changes)


# Paper values.  Counts and widths used only to reconstruct an actual mask are
# noted explicitly above.  The design-B/C nominal widths follow the paper's
# stated ~700 x 500 um and <1000 x 1000 um envelope constraints.
DESIGN_A = ResonatorParameters()
DESIGN_B = ResonatorParameters(
    design_name="B",
    g_c_um=20.0,
    g_l_um=5.0,
    g_r_um=10.0,
    w_um=10.0,
    z0_ohm=200.0,
    nominal_width_um=700.0,
    capacitor_finger_count=10,
    inductor_turns=10,
)
DESIGN_C = ResonatorParameters(
    design_name="C",
    g_c_um=80.0,
    g_l_um=10.0,
    g_r_um=10.0,
    w_um=10.0,
    z0_ohm=300.0,
    nominal_width_um=1000.0,
    capacitor_finger_count=6,
    inductor_turns=10,
    capacitor_w_um=40.0,
    inductor_w_um=10.0,
)

PRESETS: Mapping[str, ResonatorParameters] = MappingProxyType(
    {"A": DESIGN_A, "B": DESIGN_B, "C": DESIGN_C}
)


def get_preset(name: str = "A") -> ResonatorParameters:
    """Return immutable paper preset A, B or C."""

    key = name.strip().upper()
    try:
        return PRESETS[key]
    except KeyError as exc:
        choices = ", ".join(PRESETS)
        raise ValueError(f"unknown design preset {name!r}; choose {choices}") from exc


_ALIASES: Mapping[str, str] = MappingProxyType(
    {
        "name": "design_name",
        "design": "design_name",
        "gC": "g_c_um",
        "gc": "g_c_um",
        "g_c": "g_c_um",
        "gL": "g_l_um",
        "gl": "g_l_um",
        "g_l": "g_l_um",
        "gR": "g_r_um",
        "gr": "g_r_um",
        "g_r": "g_r_um",
        "w": "w_um",
        "trace_width_um": "w_um",
        "Z0": "z0_ohm",
        "z0": "z0_ohm",
        "nominal_width": "nominal_width_um",
        "finger_count": "capacitor_finger_count",
        "turn_count": "inductor_turns",
        "cap_w_um": "capacitor_w_um",
        "ind_w_um": "inductor_w_um",
        "left_layer": "left_port_layer",
        "right_layer": "right_port_layer",
    }
)

_FIELD_NAMES = frozenset(field.name for field in fields(ResonatorParameters))
_TABLE_NAMES = ("resonator", "geometry", "cpw", "ports", "layers", "stack")


def _flatten_config(data: Mapping[str, Any]) -> dict[str, Any]:
    flattened: dict[str, Any] = {}
    for key, value in data.items():
        if key in _TABLE_NAMES:
            if not isinstance(value, Mapping):
                raise TypeError(f"TOML table [{key}] must contain key/value pairs")
            for nested_key, nested_value in value.items():
                if nested_key in flattened:
                    raise ValueError(f"duplicate configuration key {nested_key!r}")
                flattened[nested_key] = nested_value
        else:
            if key in flattened:
                raise ValueError(f"duplicate configuration key {key!r}")
            flattened[key] = value
    return flattened


def parameters_from_mapping(
    values: Mapping[str, Any], *, preset: str | None = None
) -> ResonatorParameters:
    """Build parameters from a mapping, optionally updating a paper preset.

    Both Python field names and paper-style aliases such as ``gC``, ``gL``,
    ``gR``, ``w`` and ``Z0`` are accepted.  Unknown keys are rejected so a
    misspelling cannot silently change the simulated structure.
    """

    if not isinstance(values, Mapping):
        raise TypeError("parameter values must be a mapping")
    flat = _flatten_config(values)
    inline_preset = flat.pop("preset", None)
    if inline_preset is not None and not isinstance(inline_preset, str):
        raise TypeError("preset must be a string")
    if preset is not None and inline_preset is not None:
        if preset.strip().upper() != inline_preset.strip().upper():
            raise ValueError("preset argument conflicts with the TOML preset value")
    selected = inline_preset or preset or "A"
    base = get_preset(selected)

    updates: dict[str, Any] = {}
    for raw_name, value in flat.items():
        name = _ALIASES.get(raw_name, raw_name)
        if name not in _FIELD_NAMES:
            allowed = ", ".join(sorted(_FIELD_NAMES))
            raise ValueError(
                f"unknown resonator parameter {raw_name!r}; valid fields: {allowed}"
            )
        if name in updates:
            raise ValueError(f"parameter {name!r} was specified more than once")
        updates[name] = value
    return replace(base, **updates)


def load_parameters(
    path: str | PathLike[str], *, preset: str | None = None
) -> ResonatorParameters:
    """Load a validated configuration from TOML using the standard library."""

    config_path = Path(path)
    try:
        with config_path.open("rb") as stream:
            values = tomllib.load(stream)
    except tomllib.TOMLDecodeError as exc:
        raise ValueError(f"invalid TOML in {config_path}: {exc}") from exc
    return parameters_from_mapping(values, preset=preset)


__all__ = [
    "COUNT_INFERENCE_NOTE",
    "DESIGN_A",
    "DESIGN_B",
    "DESIGN_C",
    "PRESETS",
    "ResonatorParameters",
    "get_preset",
    "load_parameters",
    "parameters_from_mapping",
]
