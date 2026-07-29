"""Electromagnetic stackup, meshing, and solver parameter models."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from pathlib import Path
from typing import Any
import xml.etree.ElementTree as ET

from .errors import ConfigurationError


@dataclass(frozen=True, slots=True)
class StackupParameters:
    """Fabrication assumptions not specified completely in the paper.

    Coordinates are micrometres.  The scalar sapphire permittivity is a
    practical gds2palace input; a tensor can be substituted in Palace later.
    """

    substrate_thickness_um: float = 430.0
    air_thickness_um: float = 300.0
    sapphire_permittivity: float = 9.4
    sapphire_loss_tangent: float = 1.0e-6
    nb_thickness_um: float = 0.2
    nb_conductivity_s_per_m: float = 1.0e10
    metal_layer: int = 1

    def validate(self) -> None:
        positive = {
            "substrate_thickness_um": self.substrate_thickness_um,
            "air_thickness_um": self.air_thickness_um,
            "sapphire_permittivity": self.sapphire_permittivity,
            "nb_thickness_um": self.nb_thickness_um,
            "nb_conductivity_s_per_m": self.nb_conductivity_s_per_m,
        }
        for name, value in positive.items():
            if value <= 0:
                raise ConfigurationError(f"{name} must be positive, got {value!r}")
        if self.sapphire_loss_tangent < 0:
            raise ConfigurationError("sapphire_loss_tangent cannot be negative")
        if self.metal_layer <= 0:
            raise ConfigurationError("metal_layer must be a positive GDS layer")


@dataclass(frozen=True, slots=True)
class MeshParameters:
    dielectric_margin_um: float = 150.0
    air_margin_xy_um: float = 250.0
    air_below_um: float = 100.0
    air_above_um: float = 300.0
    refined_cell_size_um: float = 4.0
    max_cell_size_um: float = 80.0
    cells_per_wavelength: int = 12
    adaptive_iterations: int = 0
    boundary: str = "ABC"
    preprocess_gds: bool = True

    def validate(self) -> None:
        for name in (
            "dielectric_margin_um",
            "air_margin_xy_um",
            "air_below_um",
            "air_above_um",
            "refined_cell_size_um",
            "max_cell_size_um",
        ):
            if getattr(self, name) <= 0:
                raise ConfigurationError(f"{name} must be positive")
        if self.cells_per_wavelength < 10:
            raise ConfigurationError("cells_per_wavelength must be at least 10")
        if self.adaptive_iterations < 0:
            raise ConfigurationError("adaptive_iterations cannot be negative")
        if self.boundary.upper() not in {"ABC", "PEC"}:
            raise ConfigurationError("boundary must be ABC or PEC")


@dataclass(frozen=True, slots=True)
class EigenmodeParameters:
    target_frequency_ghz: float = 3.0
    search_max_frequency_ghz: float = 12.0
    modes: int = 8
    save_modes: int = 8
    order: int = 2
    tolerance: float = 1.0e-8
    linear_tolerance: float = 1.0e-8
    linear_max_iterations: int = 500
    mpi_processes: int = 4
    metal_model: str = "pec"
    keep_passive_ports: bool = False

    def validate(self) -> None:
        if self.target_frequency_ghz <= 0:
            raise ConfigurationError("target_frequency_ghz must be nonzero and positive")
        if self.search_max_frequency_ghz <= self.target_frequency_ghz:
            raise ConfigurationError(
                "search_max_frequency_ghz must exceed target_frequency_ghz"
            )
        if self.modes <= 0 or self.save_modes < 0:
            raise ConfigurationError("modes must be positive and save_modes nonnegative")
        if self.save_modes > self.modes:
            raise ConfigurationError("save_modes cannot exceed modes")
        if self.order not in {1, 2, 3}:
            raise ConfigurationError("Palace basis order must be 1, 2, or 3")
        if self.tolerance <= 0 or self.linear_tolerance <= 0:
            raise ConfigurationError("solver tolerances must be positive")
        if self.linear_max_iterations <= 0 or self.mpi_processes <= 0:
            raise ConfigurationError("iteration and MPI counts must be positive")
        if self.metal_model.lower() not in {"pec", "conductivity"}:
            raise ConfigurationError("metal_model must be 'pec' or 'conductivity'")


@dataclass(frozen=True, slots=True)
class LumpedTargets:
    """Ideal LC targets associated with a modal impedance and frequency."""

    frequency_ghz: float
    impedance_ohm: float
    inductance_nh: float
    capacitance_ff: float


def ideal_lumped_targets(frequency_ghz: float, impedance_ohm: float) -> LumpedTargets:
    """Compute L=Z/omega and C=1/(omega*Z).

    These equations are reported in the author's thesis and are useful for
    synthesis, but the paper notes that a full resonator simulation is needed
    because isolated ideal L/C estimates overpredict the measured frequency.
    """

    if frequency_ghz <= 0 or impedance_ohm <= 0:
        raise ConfigurationError("frequency and impedance must be positive")
    omega = 2.0 * math.pi * frequency_ghz * 1.0e9
    inductance_h = impedance_ohm / omega
    capacitance_f = 1.0 / (omega * impedance_ohm)
    return LumpedTargets(
        frequency_ghz=frequency_ghz,
        impedance_ohm=impedance_ohm,
        inductance_nh=inductance_h * 1.0e9,
        capacitance_ff=capacitance_f * 1.0e15,
    )


def dataclass_overrides(cls: type[Any], values: dict[str, Any] | None = None) -> Any:
    """Construct a parameter dataclass while rejecting misspelled keys."""

    values = values or {}
    allowed = set(cls.__dataclass_fields__)
    unknown = set(values) - allowed
    if unknown:
        names = ", ".join(sorted(unknown))
        raise ConfigurationError(f"Unknown {cls.__name__} setting(s): {names}")
    result = cls(**values)
    result.validate()
    return result


def stackup_xml(parameters: StackupParameters) -> str:
    """Return a gds2palace-compatible sapphire/Nb stackup XML document."""

    parameters.validate()
    root = ET.Element("Stackup", {"schemaVersion": "2.0"})
    materials = ET.SubElement(root, "Materials")
    ET.SubElement(
        materials,
        "Material",
        {
            "Name": "Nb",
            "Type": "Conductor",
            "Permittivity": "1",
            "DielectricLossTangent": "0",
            "Conductivity": f"{parameters.nb_conductivity_s_per_m:.12g}",
            "Color": "d9e7ff",
        },
    )
    ET.SubElement(
        materials,
        "Material",
        {
            "Name": "Sapphire",
            "Type": "Dielectric",
            "Permittivity": f"{parameters.sapphire_permittivity:.12g}",
            "DielectricLossTangent": f"{parameters.sapphire_loss_tangent:.12g}",
            "Conductivity": "0",
            "Color": "6baed6",
        },
    )
    ET.SubElement(
        materials,
        "Material",
        {
            "Name": "AIR",
            "Type": "Dielectric",
            "Permittivity": "1",
            "DielectricLossTangent": "0",
            "Conductivity": "0",
            "Color": "d0d0d0",
        },
    )
    elayers = ET.SubElement(root, "ELayers", {"LengthUnit": "um"})
    dielectrics = ET.SubElement(elayers, "Dielectrics")
    # gds2palace expects top-to-bottom order and computes z in reverse.
    ET.SubElement(
        dielectrics,
        "Dielectric",
        {
            "Name": "AIR",
            "Material": "AIR",
            "Thickness": f"{parameters.air_thickness_um:.12g}",
        },
    )
    ET.SubElement(
        dielectrics,
        "Dielectric",
        {
            "Name": "Sapphire",
            "Material": "Sapphire",
            "Thickness": f"{parameters.substrate_thickness_um:.12g}",
        },
    )
    layers = ET.SubElement(elayers, "Layers")
    ET.SubElement(
        layers,
        "Substrate",
        {"Offset": f"{parameters.substrate_thickness_um:.12g}"},
    )
    ET.SubElement(
        layers,
        "Layer",
        {
            "Name": "Nb",
            "Type": "conductor",
            "Zmin": "0",
            "Zmax": f"{parameters.nb_thickness_um:.12g}",
            "Material": "Nb",
            "Layer": str(parameters.metal_layer),
        },
    )
    ET.indent(root, space="  ")
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + ET.tostring(
        root, encoding="unicode"
    ) + "\n"


def write_stackup_xml(path: str | Path, parameters: StackupParameters) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(stackup_xml(parameters), encoding="utf-8")
    return target


def settings_manifest(
    stackup: StackupParameters,
    mesh: MeshParameters,
    eigenmode: EigenmodeParameters,
) -> dict[str, Any]:
    return {
        "stackup": asdict(stackup),
        "mesh": asdict(mesh),
        "eigenmode": asdict(eigenmode),
    }
