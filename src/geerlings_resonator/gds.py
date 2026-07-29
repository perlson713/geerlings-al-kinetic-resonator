"""Optional gdsfactory export for the dependency-free planar geometry."""

from __future__ import annotations

from dataclasses import dataclass
from os import PathLike
from pathlib import Path
import re
from typing import Any

from .exceptions import MissingOptionalDependencyError
from .geometry import LayoutGeometry, generate_layout
from .parameters import ResonatorParameters


LayerSpec = tuple[int, int]


@dataclass(frozen=True, slots=True)
class GDSLayerMap:
    """GDSII layer/datatype pairs used by the exporter."""

    metal: LayerSpec = (1, 0)
    left_port: LayerSpec = (201, 0)
    right_port: LayerSpec = (202, 0)

    @classmethod
    def from_parameters(cls, parameters: ResonatorParameters) -> "GDSLayerMap":
        return cls(
            metal=(parameters.metal_layer, parameters.metal_datatype),
            left_port=(parameters.left_port_layer, parameters.port_datatype),
            right_port=(parameters.right_port_layer, parameters.port_datatype),
        )

    def __post_init__(self) -> None:
        for name in ("metal", "left_port", "right_port"):
            spec = getattr(self, name)
            if (
                not isinstance(spec, tuple)
                or len(spec) != 2
                or any(isinstance(value, bool) or not isinstance(value, int) for value in spec)
                or any(not 0 <= value <= 65535 for value in spec)
            ):
                raise ValueError(
                    f"{name} must be a (layer, datatype) pair in the range 0..65535"
                )


def _coerce_layout(
    value: LayoutGeometry | ResonatorParameters | None,
) -> LayoutGeometry:
    if value is None:
        return generate_layout()
    if isinstance(value, LayoutGeometry):
        return value
    if isinstance(value, ResonatorParameters):
        return generate_layout(value)
    raise TypeError("expected LayoutGeometry, ResonatorParameters, or None")


def _gdsfactory() -> Any:
    try:
        import gdsfactory as gf
    except ImportError as exc:
        raise MissingOptionalDependencyError(
            "GDS export requires the optional 'gdsfactory' package. Install "
            "the project with `python -m pip install -e .` or run "
            "`python -m pip install gdsfactory`."
        ) from exc
    return gf


def _cell_name(layout: LayoutGeometry, requested: str | None) -> str:
    raw = requested or f"geerlings_design_{layout.parameters.design_name}"
    cleaned = re.sub(r"[^A-Za-z0-9_$?]", "_", raw).strip("_")
    if not cleaned:
        cleaned = "geerlings_resonator"
    # GDS cell names are commonly limited to 32 characters by downstream CAD.
    return cleaned[:32]


def build_gds_component(
    value: LayoutGeometry | ResonatorParameters | None = None,
    *,
    layers: GDSLayerMap | None = None,
    cell_name: str | None = None,
    include_ports: bool = True,
) -> Any:
    """Build and return a gdsfactory Component without writing it.

    Positive niobium metal, including resonator, CPW signal and surrounding
    ground, is placed on layer 1/datatype 0 by default.  Left/right setupEM
    driven-port markers are placed on layers 201/202.
    """

    layout = _coerce_layout(value)
    gf = _gdsfactory()
    layer_map = layers or GDSLayerMap.from_parameters(layout.parameters)
    base_name = _cell_name(layout, cell_name)
    # gdsfactory's default KCLayout is process-global and rejects a duplicate
    # cell name before write_gds() can deduplicate it.  This matters for the
    # multi-point sweep command, which builds several Design-A cells in one
    # Python process.
    try:
        unique_name = gf.kcl.layout.unique_cell_name(base_name)
    except AttributeError:
        # Compatibility fallback for older gdsfactory 9 releases.
        unique_name = base_name
    try:
        component = gf.Component(name=unique_name)
    except ValueError as first_error:
        for suffix in range(1, 10_000):
            retry_name = f"{base_name[:27]}_{suffix:04d}"
            try:
                component = gf.Component(name=retry_name)
                break
            except ValueError:
                continue
        else:
            raise first_error
    for polygon in layout.metal_polygons:
        component.add_polygon(polygon.tuples(), layer=layer_map.metal)
    if include_ports:
        for port in layout.ports:
            layer = layer_map.left_port if port.name == "left" else layer_map.right_port
            component.add_polygon(port.polygon.tuples(), layer=layer)

    # Component.info values must remain JSON serialisable across gdsfactory
    # releases.  They make paper assumptions discoverable in a saved project.
    try:
        component.info.update(
            {
                "paper_design": layout.parameters.design_name,
                "z0_ohm": layout.parameters.z0_ohm,
                "capacitor_finger_count": layout.parameters.capacitor_finger_count,
                "inductor_half_spans_per_side": layout.parameters.inductor_turns,
                "counts_source": "Fig.1/inferred; see Python parameters",
            }
        )
    except (AttributeError, TypeError, ValueError):
        # Metadata APIs have changed between gdsfactory releases; geometry is
        # the required result and is deliberately version-agnostic.
        pass
    return component


def write_gds(
    value: LayoutGeometry | ResonatorParameters | None,
    path: str | PathLike[str],
    *,
    layers: GDSLayerMap | None = None,
    cell_name: str | None = None,
    include_ports: bool = True,
) -> Path:
    """Write a gdsfactory GDS file and return its path."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    component = build_gds_component(
        value,
        layers=layers,
        cell_name=cell_name,
        include_ports=include_ports,
    )
    component.write_gds(str(output))
    return output


export_gds = write_gds


__all__ = [
    "GDSLayerMap",
    "LayerSpec",
    "build_gds_component",
    "export_gds",
    "write_gds",
]
