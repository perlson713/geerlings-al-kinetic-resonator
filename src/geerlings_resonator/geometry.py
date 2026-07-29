"""Dependency-free topology and planar geometry for the Fig. 1 resonator.

The module deliberately uses small immutable Python objects rather than a CAD
kernel.  Geometry can therefore be inspected and unit-tested without
gdsfactory, build123d, shapely, or a mesher.  Coordinates and widths are in
micrometres.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, hypot, isfinite, pi, sin
from types import MappingProxyType
from typing import Iterable, Iterator, Mapping, Sequence

from .parameters import ResonatorParameters


_EPSILON = 1.0e-12


def _finite_coordinate(name: str, value: object) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a real number")
    if not isfinite(float(value)):
        raise ValueError(f"{name} must be finite")


@dataclass(frozen=True, slots=True, order=True)
class Point2D:
    """A point in the layout plane, in micrometres."""

    x: float
    y: float

    def __post_init__(self) -> None:
        _finite_coordinate("x", self.x)
        _finite_coordinate("y", self.y)

    def as_tuple(self) -> tuple[float, float]:
        return (float(self.x), float(self.y))


@dataclass(frozen=True, slots=True)
class Bounds2D:
    """An axis-aligned bounding box in micrometres."""

    x_min: float
    y_min: float
    x_max: float
    y_max: float

    def __post_init__(self) -> None:
        for name in ("x_min", "y_min", "x_max", "y_max"):
            _finite_coordinate(name, getattr(self, name))
        if self.x_max <= self.x_min or self.y_max <= self.y_min:
            raise ValueError("bounds must have positive width and height")

    @property
    def width(self) -> float:
        return self.x_max - self.x_min

    @property
    def height(self) -> float:
        return self.y_max - self.y_min

    @property
    def center(self) -> Point2D:
        return Point2D(
            0.5 * (self.x_min + self.x_max),
            0.5 * (self.y_min + self.y_max),
        )

    def expanded(self, margin_um: float) -> "Bounds2D":
        _finite_coordinate("margin_um", margin_um)
        if margin_um < 0.0:
            raise ValueError("margin_um cannot be negative")
        return Bounds2D(
            self.x_min - margin_um,
            self.y_min - margin_um,
            self.x_max + margin_um,
            self.y_max + margin_um,
        )


def _signed_area(points: Sequence[Point2D]) -> float:
    return 0.5 * sum(
        first.x * second.y - second.x * first.y
        for first, second in zip(points, (*points[1:], points[0]), strict=True)
    )


@dataclass(frozen=True, slots=True)
class Polygon2D:
    """A simple, counter-clockwise polygon with a semantic role."""

    name: str
    points: tuple[Point2D, ...]
    role: str

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("polygon name cannot be empty")
        if not self.role:
            raise ValueError("polygon role cannot be empty")
        points = tuple(self.points)
        if len(points) >= 2 and points[0] == points[-1]:
            points = points[:-1]
        if len(points) < 3:
            raise ValueError("a polygon requires at least three distinct points")
        if len(set(points)) < 3:
            raise ValueError("a polygon requires at least three distinct points")
        area = _signed_area(points)
        if abs(area) <= _EPSILON:
            raise ValueError("polygon area must be non-zero")
        if area < 0.0:
            points = tuple(reversed(points))
        object.__setattr__(self, "points", points)

    @property
    def signed_area(self) -> float:
        return _signed_area(self.points)

    @property
    def area(self) -> float:
        return abs(self.signed_area)

    @property
    def bounds(self) -> Bounds2D:
        xs = tuple(point.x for point in self.points)
        ys = tuple(point.y for point in self.points)
        return Bounds2D(min(xs), min(ys), max(xs), max(ys))

    def tuples(self) -> tuple[tuple[float, float], ...]:
        return tuple(point.as_tuple() for point in self.points)


@dataclass(frozen=True, slots=True)
class Polyline2D:
    """An ordered centerline with a physical width and semantic role."""

    name: str
    points: tuple[Point2D, ...]
    width_um: float
    role: str

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("polyline name cannot be empty")
        if not self.role:
            raise ValueError("polyline role cannot be empty")
        points = tuple(self.points)
        if len(points) < 2:
            raise ValueError("a polyline requires at least two points")
        if any(first == second for first, second in zip(points, points[1:])):
            raise ValueError("a polyline cannot contain zero-length segments")
        _finite_coordinate("width_um", self.width_um)
        if self.width_um <= 0.0:
            raise ValueError("width_um must be greater than zero")
        object.__setattr__(self, "points", points)

    @property
    def length(self) -> float:
        return sum(
            hypot(second.x - first.x, second.y - first.y)
            for first, second in zip(self.points, self.points[1:])
        )

    @property
    def bounds(self) -> Bounds2D:
        half = 0.5 * self.width_um
        xs = tuple(point.x for point in self.points)
        ys = tuple(point.y for point in self.points)
        return Bounds2D(
            min(xs) - half,
            min(ys) - half,
            max(xs) + half,
            max(ys) + half,
        )


@dataclass(frozen=True, slots=True)
class Port2D:
    """A finite port marker plus its transverse reference line."""

    name: str
    polygon: Polygon2D
    line: Polyline2D
    layer: int
    datatype: int
    propagation_direction: tuple[float, float]

    def __post_init__(self) -> None:
        if self.name not in {"left", "right"}:
            raise ValueError("port name must be 'left' or 'right'")
        for field_name in ("layer", "datatype"):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"port {field_name} must be an integer")
            if not 0 <= value <= 65535:
                raise ValueError(f"port {field_name} must be between 0 and 65535")
        dx, dy = self.propagation_direction
        if abs(hypot(dx, dy) - 1.0) > 1.0e-9:
            raise ValueError("port propagation_direction must be a unit vector")


@dataclass(frozen=True, slots=True)
class CircuitElement:
    """A simple netlist-level element used to expose the LC topology."""

    name: str
    kind: str
    node_a: str
    node_b: str


@dataclass(frozen=True, slots=True)
class ResonatorTopology:
    """Topology of a meander inductor in parallel with an IDC capacitor."""

    nodes: tuple[str, str]
    elements: tuple[CircuitElement, CircuitElement]

    @property
    def is_parallel_lc(self) -> bool:
        endpoints = {(element.node_a, element.node_b) for element in self.elements}
        kinds = {element.kind for element in self.elements}
        return len(endpoints) == 1 and kinds == {"inductor", "capacitor"}


def compact_resonator_topology() -> ResonatorTopology:
    """Return the dependency-free two-node parallel-LC topology."""

    nodes = ("left_node", "right_node")
    return ResonatorTopology(
        nodes=nodes,
        elements=(
            CircuitElement("L_meander", "inductor", *nodes),
            CircuitElement("C_interdigitated", "capacitor", *nodes),
        ),
    )


@dataclass(frozen=True, slots=True)
class CenterlineGeometry:
    """Named centerlines used to build the conductor polygons."""

    inductor: Polyline2D
    feedline: Polyline2D | None
    capacitor_fingers: tuple[Polyline2D, ...]
    capacitor_buses: tuple[Polyline2D, Polyline2D]

    def __iter__(self) -> Iterator[Polyline2D]:
        yield self.inductor
        if self.feedline is not None:
            yield self.feedline
        yield from self.capacitor_buses
        yield from self.capacitor_fingers


@dataclass(frozen=True, slots=True)
class LayoutGeometry:
    """Complete positive-metal, etch-guide, centerline, and port geometry."""

    parameters: ResonatorParameters
    topology: ResonatorTopology
    centerlines: CenterlineGeometry
    metal_polygons: tuple[Polygon2D, ...]
    etch_polygons: tuple[Polygon2D, ...]
    ports: tuple[Port2D, ...]
    domain_bounds: Bounds2D

    @property
    def bounds(self) -> Bounds2D:
        return self.domain_bounds

    @property
    def polygons(self) -> tuple[Polygon2D, ...]:
        """Alias for all positive metal polygons."""

        return self.metal_polygons

    def metal_by_role(self) -> Mapping[str, tuple[Polygon2D, ...]]:
        grouped: dict[str, list[Polygon2D]] = {}
        for polygon in self.metal_polygons:
            grouped.setdefault(polygon.role, []).append(polygon)
        return MappingProxyType(
            {role: tuple(polygons) for role, polygons in grouped.items()}
        )

    def get_port(self, name: str) -> Port2D:
        for port in self.ports:
            if port.name == name:
                return port
        raise KeyError(f"layout does not contain a {name!r} port")


@dataclass(frozen=True, slots=True)
class _Dimensions:
    cap_w: float
    ind_w: float
    cutout_x_min: float
    cutout_x_max: float
    bus_left_x: float
    bus_right_x: float
    cap_bottom: float
    cap_top: float
    inductor_y0: float
    bridge_y: float
    cutout_y_min: float
    cutout_y_max: float
    feedline_y: float
    lower_slot_y_min: float
    lower_slot_y_max: float
    upper_slot_y_min: float
    upper_slot_y_max: float
    signal_x_min: float
    signal_x_max: float
    slot_x_min: float
    slot_x_max: float
    domain: Bounds2D


def _derive(parameters: ResonatorParameters) -> _Dimensions:
    cap_w = parameters.capacitor_trace_width_um
    ind_w = parameters.inductor_trace_width_um
    half_width = 0.5 * parameters.nominal_width_um
    cutout_x_min = -half_width
    cutout_x_max = half_width
    bus_left_x = cutout_x_min + parameters.g_r_um + 0.5 * cap_w
    bus_right_x = cutout_x_max - parameters.g_r_um - 0.5 * cap_w

    cap_bottom = 0.0
    cap_top = (
        parameters.capacitor_finger_count * cap_w
        + (parameters.capacitor_finger_count - 1) * parameters.g_c_um
    )
    inductor_y0 = cap_top + parameters.g_l_um + 0.5 * ind_w
    bridge_y = inductor_y0 + parameters.inductor_turns * (
        ind_w + parameters.g_l_um
    )
    cutout_y_min = cap_bottom - parameters.g_r_um
    cutout_y_max = bridge_y + 0.5 * ind_w + parameters.g_r_um

    lower_slot_y_min = cutout_y_max + parameters.coupling_ground_width_um
    lower_slot_y_max = lower_slot_y_min + parameters.feedline_gap_um
    feedline_y = lower_slot_y_max + 0.5 * parameters.feedline_width_um
    upper_slot_y_min = feedline_y + 0.5 * parameters.feedline_width_um
    upper_slot_y_max = upper_slot_y_min + parameters.feedline_gap_um

    signal_half_length = 0.5 * parameters.resolved_feedline_length_um
    signal_x_min = -signal_half_length
    signal_x_max = signal_half_length
    slot_x_min = signal_x_min - parameters.feedline_end_clearance_um
    slot_x_max = signal_x_max + parameters.feedline_end_clearance_um
    domain_half_length = (
        signal_half_length
        + parameters.feedline_end_clearance_um
        + parameters.ground_bridge_width_um
        if parameters.include_feedline
        else half_width + parameters.ground_margin_um
    )
    domain_y_max = (
        upper_slot_y_max + parameters.ground_margin_um
        if parameters.include_feedline
        else cutout_y_max + parameters.ground_margin_um
    )
    domain = Bounds2D(
        -domain_half_length,
        cutout_y_min - parameters.ground_margin_um,
        domain_half_length,
        domain_y_max,
    )
    return _Dimensions(
        cap_w=cap_w,
        ind_w=ind_w,
        cutout_x_min=cutout_x_min,
        cutout_x_max=cutout_x_max,
        bus_left_x=bus_left_x,
        bus_right_x=bus_right_x,
        cap_bottom=cap_bottom,
        cap_top=cap_top,
        inductor_y0=inductor_y0,
        bridge_y=bridge_y,
        cutout_y_min=cutout_y_min,
        cutout_y_max=cutout_y_max,
        feedline_y=feedline_y,
        lower_slot_y_min=lower_slot_y_min,
        lower_slot_y_max=lower_slot_y_max,
        upper_slot_y_min=upper_slot_y_min,
        upper_slot_y_max=upper_slot_y_max,
        signal_x_min=signal_x_min,
        signal_x_max=signal_x_max,
        slot_x_min=slot_x_min,
        slot_x_max=slot_x_max,
        domain=domain,
    )


def _rectangle(
    name: str,
    x_min: float,
    y_min: float,
    x_max: float,
    y_max: float,
    role: str,
) -> Polygon2D:
    if x_max <= x_min or y_max <= y_min:
        raise ValueError(f"rectangle {name!r} must have positive width and height")
    return Polygon2D(
        name=name,
        points=(
            Point2D(x_min, y_min),
            Point2D(x_max, y_min),
            Point2D(x_max, y_max),
            Point2D(x_min, y_max),
        ),
        role=role,
    )


def _disk(
    name: str,
    center: Point2D,
    radius: float,
    role: str,
    segments: int,
) -> Polygon2D:
    if segments < 8:
        raise ValueError("round-stroke approximation requires at least 8 segments")
    return Polygon2D(
        name=name,
        points=tuple(
            Point2D(
                center.x + radius * cos(2.0 * pi * index / segments),
                center.y + radius * sin(2.0 * pi * index / segments),
            )
            for index in range(segments)
        ),
        role=role,
    )


def stroke_polyline(
    centerline: Polyline2D, *, round_segments: int = 16
) -> tuple[Polygon2D, ...]:
    """Convert a centerline to overlapping polygons with round joins/caps.

    The operation is local and dependency-free.  Polygon overlap is
    intentional; GDS and STEP exporters preserve the same connected conductor.
    """

    half = 0.5 * centerline.width_um
    polygons: list[Polygon2D] = []
    for index, (first, second) in enumerate(
        zip(centerline.points, centerline.points[1:])
    ):
        dx = second.x - first.x
        dy = second.y - first.y
        length = hypot(dx, dy)
        nx = -dy * half / length
        ny = dx * half / length
        polygons.append(
            Polygon2D(
                name=f"{centerline.name}_segment_{index:03d}",
                points=(
                    Point2D(first.x - nx, first.y - ny),
                    Point2D(second.x - nx, second.y - ny),
                    Point2D(second.x + nx, second.y + ny),
                    Point2D(first.x + nx, first.y + ny),
                ),
                role=centerline.role,
            )
        )
    polygons.extend(
        _disk(
            f"{centerline.name}_join_{index:03d}",
            point,
            half,
            centerline.role,
            round_segments,
        )
        for index, point in enumerate(centerline.points)
    )
    return tuple(polygons)


def _inductor_centerline(
    parameters: ResonatorParameters, dimensions: _Dimensions
) -> Polyline2D:
    width = dimensions.ind_w
    pitch = width + parameters.g_l_um
    bend_radius = 0.5 * pitch
    start_left = dimensions.bus_left_x
    # Turn centerlines are inset by the semicircle radius so the outside metal
    # edges retain gR at the cutout and gL across the center opening.
    outer_turn_left = (
        dimensions.cutout_x_min
        + parameters.g_r_um
        + 0.5 * width
        + bend_radius
    )
    inner_turn_left = -(0.5 * parameters.g_l_um + 0.5 * width + bend_radius)
    if inner_turn_left <= outer_turn_left:
        raise ValueError(
            "nominal_width_um leaves no horizontal run for the inductor; "
            "increase it or reduce g_r/trace widths"
        )

    def append_turn(
        points: list[Point2D], x: float, y_bottom: float, *, bulge_right: bool
    ) -> None:
        segments = 12
        center_y = y_bottom + bend_radius
        direction = 1.0 if bulge_right else -1.0
        for segment in range(1, segments + 1):
            if segment == segments:
                points.append(Point2D(x, y_bottom + pitch))
                continue
            angle = -0.5 * pi + direction * pi * segment / segments
            points.append(
                Point2D(
                    x + bend_radius * cos(angle),
                    center_y + bend_radius * sin(angle),
                )
            )

    # Build the left ascending serpentine with the rounded hairpin turns shown
    # in Fig. 1. Mirroring it in reverse produces the right descending half;
    # the top bridge makes one continuous conductor between the IDC buses.
    left: list[Point2D] = [Point2D(start_left, dimensions.inductor_y0)]
    current_x = start_left
    for row in range(parameters.inductor_turns):
        target_x = inner_turn_left if row % 2 == 0 else outer_turn_left
        if abs(current_x - target_x) > _EPSILON:
            left.append(Point2D(target_x, dimensions.inductor_y0 + row * pitch))
            current_x = target_x
        append_turn(
            left,
            current_x,
            dimensions.inductor_y0 + row * pitch,
            bulge_right=(current_x == inner_turn_left),
        )
    assert abs(left[-1].y - dimensions.bridge_y) <= _EPSILON

    right_descending = tuple(Point2D(-point.x, point.y) for point in reversed(left))
    points = tuple(left) + right_descending
    return Polyline2D(
        name="inductor",
        points=points,
        width_um=width,
        role="resonator_inductor",
    )


def _capacitor_centerlines(
    parameters: ResonatorParameters, dimensions: _Dimensions
) -> tuple[tuple[Polyline2D, ...], tuple[Polyline2D, Polyline2D]]:
    width = dimensions.cap_w
    left_outer = dimensions.bus_left_x - 0.5 * width
    left_inner = dimensions.bus_left_x + 0.5 * width
    right_inner = dimensions.bus_right_x - 0.5 * width
    right_outer = dimensions.bus_right_x + 0.5 * width
    left_tip_edge = right_inner - parameters.g_c_um
    right_tip_edge = left_inner + parameters.g_c_um

    fingers: list[Polyline2D] = []
    pitch = width + parameters.g_c_um
    for index in range(parameters.capacitor_finger_count):
        y = dimensions.cap_bottom + 0.5 * width + index * pitch
        if index % 2 == 0:
            points = (
                Point2D(dimensions.bus_left_x, y),
                Point2D(left_tip_edge - 0.5 * width, y),
            )
            side = "left"
        else:
            points = (
                Point2D(dimensions.bus_right_x, y),
                Point2D(right_tip_edge + 0.5 * width, y),
            )
            side = "right"
        fingers.append(
            Polyline2D(
                name=f"capacitor_finger_{index:02d}_{side}",
                points=points,
                width_um=width,
                role="resonator_capacitor",
            )
        )

    buses = (
        Polyline2D(
            name="capacitor_bus_left",
            points=(
                Point2D(dimensions.bus_left_x, dimensions.cap_bottom),
                Point2D(dimensions.bus_left_x, dimensions.cap_top),
            ),
            width_um=width,
            role="resonator_capacitor",
        ),
        Polyline2D(
            name="capacitor_bus_right",
            points=(
                Point2D(dimensions.bus_right_x, dimensions.cap_bottom),
                Point2D(dimensions.bus_right_x, dimensions.cap_top),
            ),
            width_um=width,
            role="resonator_capacitor",
        ),
    )
    # The variables document the exact metal edges used for gC; geometry is
    # built as rectangles below, while these lines are exposed for inspection.
    assert left_outer < left_inner < right_inner < right_outer
    return tuple(fingers), buses


def generate_centerlines(
    parameters: ResonatorParameters | None = None,
) -> CenterlineGeometry:
    """Generate inductor, IDC, and optional CPW centerlines."""

    parameters = ResonatorParameters() if parameters is None else parameters
    dimensions = _derive(parameters)
    fingers, buses = _capacitor_centerlines(parameters, dimensions)
    feedline = None
    if parameters.include_feedline:
        feedline = Polyline2D(
            name="cpw_feedline",
            points=(
                Point2D(dimensions.signal_x_min, dimensions.feedline_y),
                Point2D(dimensions.signal_x_max, dimensions.feedline_y),
            ),
            width_um=parameters.feedline_width_um,
            role="feedline",
        )
    return CenterlineGeometry(
        inductor=_inductor_centerline(parameters, dimensions),
        feedline=feedline,
        capacitor_fingers=fingers,
        capacitor_buses=buses,
    )


def _capacitor_polygons(
    parameters: ResonatorParameters, dimensions: _Dimensions
) -> tuple[Polygon2D, ...]:
    width = dimensions.cap_w
    half = 0.5 * width
    left_outer = dimensions.bus_left_x - half
    left_inner = dimensions.bus_left_x + half
    right_inner = dimensions.bus_right_x - half
    right_outer = dimensions.bus_right_x + half
    left_tip = right_inner - parameters.g_c_um
    right_tip = left_inner + parameters.g_c_um

    polygons = [
        _rectangle(
            "capacitor_bus_left",
            left_outer,
            dimensions.cap_bottom,
            left_inner,
            dimensions.cap_top,
            "resonator_capacitor",
        ),
        _rectangle(
            "capacitor_bus_right",
            right_inner,
            dimensions.cap_bottom,
            right_outer,
            dimensions.cap_top,
            "resonator_capacitor",
        ),
    ]
    pitch = width + parameters.g_c_um
    for index in range(parameters.capacitor_finger_count):
        y_min = dimensions.cap_bottom + index * pitch
        y_max = y_min + width
        if index % 2 == 0:
            x_min, x_max, side = left_outer, left_tip, "left"
        else:
            x_min, x_max, side = right_tip, right_outer, "right"
        polygons.append(
            _rectangle(
                f"capacitor_finger_{index:02d}_{side}",
                x_min,
                y_min,
                x_max,
                y_max,
                "resonator_capacitor",
            )
        )
    return tuple(polygons)


def _connector_polygons(
    dimensions: _Dimensions,
) -> tuple[Polygon2D, Polygon2D]:
    half = 0.5 * dimensions.ind_w
    y_min = dimensions.cap_top - min(half, 0.5 * dimensions.cap_w)
    y_max = dimensions.inductor_y0
    return (
        _rectangle(
            "node_connector_left",
            dimensions.bus_left_x - half,
            y_min,
            dimensions.bus_left_x + half,
            y_max,
            "resonator_connector",
        ),
        _rectangle(
            "node_connector_right",
            dimensions.bus_right_x - half,
            y_min,
            dimensions.bus_right_x + half,
            y_max,
            "resonator_connector",
        ),
    )


def _etch_polygons(
    parameters: ResonatorParameters, dimensions: _Dimensions
) -> tuple[Polygon2D, ...]:
    polygons = [
        _rectangle(
            "resonator_ground_cutout",
            dimensions.cutout_x_min,
            dimensions.cutout_y_min,
            dimensions.cutout_x_max,
            dimensions.cutout_y_max,
            "ground_cutout",
        )
    ]
    if not parameters.include_feedline:
        return tuple(polygons)
    polygons.extend(
        (
            _rectangle(
                "cpw_gap_lower",
                dimensions.slot_x_min,
                dimensions.lower_slot_y_min,
                dimensions.slot_x_max,
                dimensions.lower_slot_y_max,
                "cpw_gap",
            ),
            _rectangle(
                "cpw_gap_upper",
                dimensions.slot_x_min,
                dimensions.upper_slot_y_min,
                dimensions.slot_x_max,
                dimensions.upper_slot_y_max,
                "cpw_gap",
            ),
        )
    )
    return tuple(polygons)


def _ground_polygons(
    parameters: ResonatorParameters, dimensions: _Dimensions
) -> tuple[Polygon2D, ...]:
    domain = dimensions.domain
    polygons = [
        _rectangle(
            "ground_bottom",
            domain.x_min,
            domain.y_min,
            domain.x_max,
            dimensions.cutout_y_min,
            "ground",
        ),
        _rectangle(
            "ground_cutout_left",
            domain.x_min,
            dimensions.cutout_y_min,
            dimensions.cutout_x_min,
            dimensions.cutout_y_max,
            "ground",
        ),
        _rectangle(
            "ground_cutout_right",
            dimensions.cutout_x_max,
            dimensions.cutout_y_min,
            domain.x_max,
            dimensions.cutout_y_max,
            "ground",
        ),
    ]
    if not parameters.include_feedline:
        polygons.append(
            _rectangle(
                "ground_top",
                domain.x_min,
                dimensions.cutout_y_max,
                domain.x_max,
                domain.y_max,
                "ground",
            )
        )
        return tuple(polygons)
    polygons.extend(
        (
            _rectangle(
                "ground_coupling_strip",
                domain.x_min,
                dimensions.cutout_y_max,
                domain.x_max,
                dimensions.lower_slot_y_min,
                "ground",
            ),
            _rectangle(
                "ground_top",
                domain.x_min,
                dimensions.upper_slot_y_max,
                domain.x_max,
                domain.y_max,
                "ground",
            ),
            _rectangle(
                "ground_cpw_bridge_left",
                domain.x_min,
                dimensions.lower_slot_y_min,
                dimensions.slot_x_min,
                dimensions.upper_slot_y_max,
                "ground",
            ),
            _rectangle(
                "ground_cpw_bridge_right",
                dimensions.slot_x_max,
                dimensions.lower_slot_y_min,
                domain.x_max,
                dimensions.upper_slot_y_max,
                "ground",
            ),
        )
    )
    return tuple(polygons)


def _ports(
    parameters: ResonatorParameters, dimensions: _Dimensions
) -> tuple[Port2D, ...]:
    if not parameters.include_feedline or not parameters.include_ports:
        return ()
    # A gds2palace in-plane lumped-port sheet spans dielectric between the
    # CPW signal and one ground terminal.  It must not cross the signal metal
    # and continue through the opposite slot, which would overlap a conductor.
    y_min = dimensions.lower_slot_y_min
    y_max = dimensions.lower_slot_y_max
    center_y = 0.5 * (y_min + y_max)
    line_width = min(parameters.port_depth_um, parameters.feedline_width_um) / 4.0
    left_line_x = dimensions.signal_x_min + 0.5 * parameters.port_depth_um
    right_line_x = dimensions.signal_x_max - 0.5 * parameters.port_depth_um
    left_polygon = _rectangle(
        "port_left_rectangle",
        dimensions.signal_x_min,
        y_min,
        dimensions.signal_x_min + parameters.port_depth_um,
        y_max,
        "port",
    )
    right_polygon = _rectangle(
        "port_right_rectangle",
        dimensions.signal_x_max - parameters.port_depth_um,
        y_min,
        dimensions.signal_x_max,
        y_max,
        "port",
    )
    # center_y is retained to make the transverse construction explicit.
    assert y_min < center_y < y_max
    return (
        Port2D(
            name="left",
            polygon=left_polygon,
            line=Polyline2D(
                "port_left_line",
                (Point2D(left_line_x, y_min), Point2D(left_line_x, y_max)),
                line_width,
                "port",
            ),
            layer=parameters.left_port_layer,
            datatype=parameters.port_datatype,
            propagation_direction=(1.0, 0.0),
        ),
        Port2D(
            name="right",
            polygon=right_polygon,
            line=Polyline2D(
                "port_right_line",
                (Point2D(right_line_x, y_min), Point2D(right_line_x, y_max)),
                line_width,
                "port",
            ),
            layer=parameters.right_port_layer,
            datatype=parameters.port_datatype,
            propagation_direction=(-1.0, 0.0),
        ),
    )


def generate_layout(
    parameters: ResonatorParameters | None = None,
) -> LayoutGeometry:
    """Generate the finite layout, with optional CPW feedline and ports."""

    parameters = ResonatorParameters() if parameters is None else parameters
    if not isinstance(parameters, ResonatorParameters):
        raise TypeError("parameters must be a ResonatorParameters instance")
    dimensions = _derive(parameters)
    centerlines = generate_centerlines(parameters)
    metal: list[Polygon2D] = list(_ground_polygons(parameters, dimensions))
    if parameters.include_feedline:
        metal.append(
            _rectangle(
                "cpw_feedline",
                dimensions.signal_x_min,
                dimensions.feedline_y - 0.5 * parameters.feedline_width_um,
                dimensions.signal_x_max,
                dimensions.feedline_y + 0.5 * parameters.feedline_width_um,
                "feedline",
            )
        )
    metal.extend(_capacitor_polygons(parameters, dimensions))
    metal.extend(_connector_polygons(dimensions))
    metal.extend(stroke_polyline(centerlines.inductor))
    return LayoutGeometry(
        parameters=parameters,
        topology=compact_resonator_topology(),
        centerlines=centerlines,
        metal_polygons=tuple(metal),
        etch_polygons=_etch_polygons(parameters, dimensions),
        ports=_ports(parameters, dimensions),
        domain_bounds=dimensions.domain,
    )


def generate_polygons(
    parameters: ResonatorParameters | None = None,
) -> tuple[Polygon2D, ...]:
    """Return all positive metal polygons for a parameter set."""

    return generate_layout(parameters).metal_polygons


def generate_etch_polygons(
    parameters: ResonatorParameters | None = None,
) -> tuple[Polygon2D, ...]:
    """Return non-metal guide polygons for the ground cutout and CPW slots."""

    return generate_layout(parameters).etch_polygons


def generate_ports(
    parameters: ResonatorParameters | None = None,
) -> tuple[Port2D, ...]:
    """Return optional left/right CPW port rectangles and reference lines."""

    return generate_layout(parameters).ports


# Concise aliases useful in notebooks and parametric sweeps.
Point = Point2D
Polygon = Polygon2D
Polyline = Polyline2D
Layout = LayoutGeometry
make_layout = generate_layout


__all__ = [
    "Bounds2D",
    "CenterlineGeometry",
    "CircuitElement",
    "Layout",
    "LayoutGeometry",
    "Point",
    "Point2D",
    "Polygon",
    "Polygon2D",
    "Polyline",
    "Polyline2D",
    "Port2D",
    "ResonatorTopology",
    "compact_resonator_topology",
    "generate_centerlines",
    "generate_etch_polygons",
    "generate_layout",
    "generate_polygons",
    "generate_ports",
    "make_layout",
    "stroke_polyline",
]
