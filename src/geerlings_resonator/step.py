"""Optional build123d STEP export of substrate and niobium geometry."""

from __future__ import annotations

from os import PathLike
from pathlib import Path
from typing import Any

from .exceptions import MissingOptionalDependencyError
from .geometry import LayoutGeometry, Point2D, Polygon2D, generate_layout
from .parameters import ResonatorParameters


_UM_TO_MM = 1.0e-3


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


def _build123d() -> tuple[Any, Any, Any, Any, Any, Any]:
    try:
        from build123d import Compound, Face, Solid, Vector, Wire, export_step
    except ImportError as exc:
        raise MissingOptionalDependencyError(
            "STEP export requires the optional 'build123d' package. Install "
            "the project with `python -m pip install -e '.[cad]'` or run "
            "`python -m pip install build123d`."
        ) from exc
    return Compound, Face, Solid, Vector, Wire, export_step


def _extrude_polygon(
    polygon: Polygon2D,
    *,
    z_base_mm: float,
    height_mm: float,
    face_type: Any,
    solid_type: Any,
    vector_type: Any,
    wire_type: Any,
) -> Any:
    vertices = [
        vector_type(point.x * _UM_TO_MM, point.y * _UM_TO_MM, z_base_mm)
        for point in polygon.points
    ]
    wire = wire_type.make_polygon(vertices, close=True)
    face = face_type(wire)
    return solid_type.extrude(face, vector_type(0.0, 0.0, height_mm))


def _substrate_polygon(layout: LayoutGeometry) -> Polygon2D:
    bounds = layout.domain_bounds
    return Polygon2D(
        name="sapphire_substrate",
        role="substrate",
        points=(
            Point2D(bounds.x_min, bounds.y_min),
            Point2D(bounds.x_max, bounds.y_min),
            Point2D(bounds.x_max, bounds.y_max),
            Point2D(bounds.x_min, bounds.y_max),
        ),
    )


def build_step_model(
    value: LayoutGeometry | ResonatorParameters | None = None,
    *,
    include_substrate: bool = True,
    include_ground: bool = True,
    include_feedline: bool = True,
) -> Any:
    """Build a build123d Compound in millimetre STEP units.

    The substrate occupies negative z and its top surface is z=0.  The 200 nm
    default niobium layer occupies positive z, matching the paper.  Port marker
    polygons are simulation annotations and are intentionally not extruded.
    """

    layout = _coerce_layout(value)
    Compound, Face, Solid, Vector, Wire, _ = _build123d()
    metal_height = layout.parameters.metal_thickness_um * _UM_TO_MM
    substrate_height = layout.parameters.substrate_thickness_um * _UM_TO_MM

    assembly_children: list[Any] = []
    if include_substrate:
        substrate = _extrude_polygon(
            _substrate_polygon(layout),
            z_base_mm=-substrate_height,
            height_mm=substrate_height,
            face_type=Face,
            solid_type=Solid,
            vector_type=Vector,
            wire_type=Wire,
        )
        substrate.label = "sapphire"
        assembly_children.append(substrate)

    metal_solids: list[Any] = []
    for polygon in layout.metal_polygons:
        if polygon.role == "ground" and not include_ground:
            continue
        if polygon.role == "feedline" and not include_feedline:
            continue
        metal_solids.append(
            _extrude_polygon(
                polygon,
                z_base_mm=0.0,
                height_mm=metal_height,
                face_type=Face,
                solid_type=Solid,
                vector_type=Vector,
                wire_type=Wire,
            )
        )
    if metal_solids:
        # The dependency-free layout intentionally represents rounded traces
        # as overlapping rectangles/disks.  Fuse those overlaps so STEP
        # contains a valid Nb body (or a compound of disconnected Nb bodies),
        # rather than an assembly with coincident/interpenetrating solids.
        niobium = (
            metal_solids[0]
            if len(metal_solids) == 1
            else metal_solids[0].fuse(*metal_solids[1:])
        )
        niobium.label = "niobium"
        assembly_children.append(niobium)
    if not assembly_children:
        raise ValueError("STEP model would contain no solids")
    return Compound(children=assembly_children)


def write_step(
    value: LayoutGeometry | ResonatorParameters | None,
    path: str | PathLike[str],
    *,
    include_substrate: bool = True,
    include_ground: bool = True,
    include_feedline: bool = True,
) -> Path:
    """Write a build123d STEP assembly and return its path."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    model = build_step_model(
        value,
        include_substrate=include_substrate,
        include_ground=include_ground,
        include_feedline=include_feedline,
    )
    *_, backend_export_step = _build123d()
    backend_export_step(model, str(output))
    return output


export_step = write_step


__all__ = ["build_step_model", "export_step", "write_step"]
