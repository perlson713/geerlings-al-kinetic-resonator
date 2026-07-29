"""Unified TOML configuration for geometry, mesh, and Palace stages."""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Iterable, Mapping
import tomllib

from .em import (
    EigenmodeParameters,
    MeshParameters,
    StackupParameters,
    dataclass_overrides,
)
from .errors import ConfigurationError
from .parameters import ResonatorParameters, parameters_from_mapping


@dataclass(frozen=True, slots=True)
class ProjectParameters:
    name: str = "geerlings_design_a"
    output_directory: str = "build/design_a"

    def validate(self) -> None:
        if not self.name.strip():
            raise ConfigurationError("project.name cannot be empty")
        if not self.output_directory.strip():
            raise ConfigurationError("project.output_directory cannot be empty")


@dataclass(frozen=True, slots=True)
class ProjectConfig:
    project: ProjectParameters
    resonator: ResonatorParameters
    stackup: StackupParameters
    mesh: MeshParameters
    eigenmode: EigenmodeParameters
    source: Path | None = None

    def manifest(self) -> dict[str, Any]:
        return {
            "project": asdict(self.project),
            "resonator": asdict(self.resonator),
            "stackup": asdict(self.stackup),
            "mesh": asdict(self.mesh),
            "eigenmode": asdict(self.eigenmode),
            "source": str(self.source) if self.source else None,
        }


_GEOMETRY_TABLES = {"resonator", "geometry", "cpw", "ports", "layers", "stack"}
_TOP_LEVEL = {"preset", "project", "stackup", "mesh", "eigenmode", *_GEOMETRY_TABLES}


def _parse_override_value(raw: str) -> Any:
    try:
        return tomllib.loads(f"value = {raw}")["value"]
    except tomllib.TOMLDecodeError:
        # Unquoted CLI strings such as mesh.boundary=PEC are ergonomic.
        return raw


def apply_overrides(
    data: Mapping[str, Any], overrides: Iterable[str] = ()
) -> dict[str, Any]:
    """Apply dotted ``section.key=value`` overrides to a copied mapping."""

    copied: dict[str, Any] = {
        key: dict(value) if isinstance(value, Mapping) else value
        for key, value in data.items()
    }
    for expression in overrides:
        if "=" not in expression:
            raise ConfigurationError(
                f"Override {expression!r} must use section.key=value"
            )
        dotted, raw_value = expression.split("=", 1)
        path = [part.strip() for part in dotted.split(".") if part.strip()]
        if len(path) != 2:
            raise ConfigurationError(
                f"Override {expression!r} must have exactly section.key"
            )
        section, key = path
        table = copied.setdefault(section, {})
        if not isinstance(table, dict):
            raise ConfigurationError(f"{section!r} is not a TOML table")
        table[key] = _parse_override_value(raw_value.strip())
    return copied


def project_from_mapping(
    values: Mapping[str, Any], *, source: Path | None = None
) -> ProjectConfig:
    unknown = set(values) - _TOP_LEVEL
    if unknown:
        raise ConfigurationError(
            "Unknown top-level TOML key(s): " + ", ".join(sorted(unknown))
        )
    geometry_values: dict[str, Any] = {}
    if "preset" in values:
        geometry_values["preset"] = values["preset"]
    for table_name in _GEOMETRY_TABLES:
        if table_name in values:
            geometry_values[table_name] = values[table_name]
    try:
        resonator = parameters_from_mapping(geometry_values)
        project = dataclass_overrides(ProjectParameters, values.get("project"))
        stackup = dataclass_overrides(StackupParameters, values.get("stackup"))
        mesh = dataclass_overrides(MeshParameters, values.get("mesh"))
        eigenmode = dataclass_overrides(EigenmodeParameters, values.get("eigenmode"))
    except (TypeError, ValueError) as exc:
        if isinstance(exc, ConfigurationError):
            raise
        raise ConfigurationError(str(exc)) from exc
    if resonator.metal_layer != stackup.metal_layer:
        raise ConfigurationError(
            "resonator metal_layer and stackup.metal_layer must match "
            f"({resonator.metal_layer} != {stackup.metal_layer})"
        )
    # STEP preview and simulation stackup should share physical thicknesses.
    geometry_stack = values.get("stack", {})
    if not isinstance(geometry_stack, Mapping):
        raise ConfigurationError("[stack] must be a TOML table")
    changes: dict[str, float] = {}
    if "metal_thickness_um" not in geometry_stack:
        changes["metal_thickness_um"] = stackup.nb_thickness_um
    if "substrate_thickness_um" not in geometry_stack:
        changes["substrate_thickness_um"] = stackup.substrate_thickness_um
    if changes:
        resonator = replace(resonator, **changes)
    return ProjectConfig(project, resonator, stackup, mesh, eigenmode, source)


def load_project(
    path: str | Path, *, overrides: Iterable[str] = ()
) -> ProjectConfig:
    source = Path(path).resolve()
    try:
        with source.open("rb") as handle:
            values = tomllib.load(handle)
    except OSError as exc:
        raise ConfigurationError(f"Unable to read configuration {source}: {exc}") from exc
    except tomllib.TOMLDecodeError as exc:
        raise ConfigurationError(f"Invalid TOML in {source}: {exc}") from exc
    return project_from_mapping(apply_overrides(values, overrides), source=source)
