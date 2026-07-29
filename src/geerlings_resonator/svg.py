"""Deterministic, dependency-free SVG previews of generated layouts."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from os import PathLike
from pathlib import Path

from .geometry import LayoutGeometry, Point2D, generate_layout
from .parameters import COUNT_INFERENCE_NOTE, ResonatorParameters


@dataclass(frozen=True, slots=True)
class SVGStyle:
    """Colours used by :func:`render_svg`."""

    substrate: str = "#111820"
    ground: str = "#d8dde3"
    feedline: str = "#f8fafc"
    resonator: str = "#f6c453"
    port_left: str = "#24a8ff"
    port_right: str = "#f97391"
    etch_guide: str = "#73808c"


def _number(value: float) -> str:
    if abs(value) < 5.0e-13:
        value = 0.0
    return f"{value:.6f}".rstrip("0").rstrip(".") or "0"


def _layout(value: LayoutGeometry | ResonatorParameters | None) -> LayoutGeometry:
    if value is None:
        return generate_layout()
    if isinstance(value, LayoutGeometry):
        return value
    if isinstance(value, ResonatorParameters):
        return generate_layout(value)
    raise TypeError("expected LayoutGeometry, ResonatorParameters, or None")


def _svg_point(point: Point2D, layout: LayoutGeometry) -> tuple[float, float]:
    bounds = layout.domain_bounds
    return (point.x - bounds.x_min, bounds.y_max - point.y)


def _points_attribute(points: tuple[Point2D, ...], layout: LayoutGeometry) -> str:
    return " ".join(
        f"{_number(x)},{_number(y)}" for x, y in (_svg_point(p, layout) for p in points)
    )


def render_svg(
    value: LayoutGeometry | ResonatorParameters | None = None,
    *,
    style: SVGStyle = SVGStyle(),
    show_ports: bool = True,
    show_etch_guides: bool = False,
) -> str:
    """Return a stable SVG string for a layout or parameter set.

    No timestamp, random identifier, or CAD-library metadata is emitted, so
    identical parameters produce byte-identical output.
    """

    layout = _layout(value)
    bounds = layout.domain_bounds
    width = _number(bounds.width)
    height = _number(bounds.height)
    design = escape(layout.parameters.design_name, quote=True)
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} '
            f'{height}" width="{width}" height="{height}" role="img" '
            f'aria-label="Geerlings compact resonator design {design}">'
        ),
        f"  <title>Geerlings compact resonator - design {design}</title>",
        f"  <desc>{escape(COUNT_INFERENCE_NOTE)}</desc>",
        (
            f'  <rect id="substrate" x="0" y="0" width="{width}" '
            f'height="{height}" fill="{escape(style.substrate, quote=True)}"/>'
        ),
        '  <g id="metal">',
    ]

    for polygon in layout.metal_polygons:
        if polygon.role == "ground":
            colour = style.ground
        elif polygon.role == "feedline":
            colour = style.feedline
        else:
            colour = style.resonator
        lines.append(
            "    "
            f'<polygon data-name="{escape(polygon.name, quote=True)}" '
            f'data-role="{escape(polygon.role, quote=True)}" '
            f'points="{_points_attribute(polygon.points, layout)}" '
            f'fill="{escape(colour, quote=True)}"/>'
        )
    lines.append("  </g>")

    if show_etch_guides:
        lines.append(
            '  <g id="etch-guides" fill="none" stroke-width="0.75" '
            f'stroke="{escape(style.etch_guide, quote=True)}" '
            'stroke-dasharray="3 3">'
        )
        for polygon in layout.etch_polygons:
            lines.append(
                "    "
                f'<polygon data-name="{escape(polygon.name, quote=True)}" '
                f'points="{_points_attribute(polygon.points, layout)}"/>'
            )
        lines.append("  </g>")

    if show_ports and layout.ports:
        lines.append('  <g id="ports" fill="none" stroke-width="1">')
        for port in layout.ports:
            colour = style.port_left if port.name == "left" else style.port_right
            lines.append(
                "    "
                f'<polygon data-name="port-{port.name}" data-layer="{port.layer}" '
                f'points="{_points_attribute(port.polygon.points, layout)}" '
                f'stroke="{escape(colour, quote=True)}" stroke-dasharray="2 2"/>'
            )
            first, second = port.line.points
            x1, y1 = _svg_point(first, layout)
            x2, y2 = _svg_point(second, layout)
            lines.append(
                "    "
                f'<line x1="{_number(x1)}" y1="{_number(y1)}" '
                f'x2="{_number(x2)}" y2="{_number(y2)}" '
                f'stroke="{escape(colour, quote=True)}"/>'
            )
        lines.append("  </g>")

    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def write_svg(
    value: LayoutGeometry | ResonatorParameters | None,
    path: str | PathLike[str],
    *,
    style: SVGStyle = SVGStyle(),
    show_ports: bool = True,
    show_etch_guides: bool = False,
) -> Path:
    """Write a deterministic SVG preview and return its resolved input path."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        render_svg(
            value,
            style=style,
            show_ports=show_ports,
            show_etch_guides=show_etch_guides,
        ),
        encoding="utf-8",
        newline="\n",
    )
    return output


export_svg = write_svg


__all__ = ["SVGStyle", "export_svg", "render_svg", "write_svg"]
