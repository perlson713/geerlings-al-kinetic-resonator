#!/usr/bin/env python3
"""Render a clean top-view PNG and deterministic SVG for a project config."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as PlotPolygon, Rectangle

from geerlings_resonator.config import load_project
from geerlings_resonator.geometry import generate_layout
from geerlings_resonator.svg import write_svg


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    project = load_project(args.config)
    layout = generate_layout(project.resonator)
    cutout = next(poly for poly in layout.etch_polygons if poly.role == "ground_cutout")
    bounds = cutout.bounds
    args.output.mkdir(parents=True, exist_ok=True)

    write_svg(layout, args.output / "layout.svg", show_ports=False, show_etch_guides=True)
    figure, axis = plt.subplots(figsize=(6.4, 6.4), constrained_layout=True)
    axis.set_facecolor("#dbe2ea")
    axis.add_patch(
        Rectangle(
            (bounds.x_min, bounds.y_min),
            bounds.width,
            bounds.height,
            facecolor="#101820",
            edgecolor="#5f7180",
            linewidth=1.2,
        )
    )
    for polygon in layout.metal_polygons:
        if polygon.role == "ground":
            continue
        axis.add_patch(
            PlotPolygon(
                polygon.tuples(),
                closed=True,
                facecolor="#f6c453",
                edgecolor="#f6c453",
                linewidth=0.25,
            )
        )
    margin = 0.035 * max(bounds.width, bounds.height)
    axis.set_xlim(bounds.x_min - margin, bounds.x_max + margin)
    axis.set_ylim(bounds.y_min - margin, bounds.y_max + margin)
    axis.set_aspect("equal")
    axis.set_axis_off()
    axis.set_title(
        f"Near-square 9 GHz resonator  |  {bounds.width:.2f} x {bounds.height:.2f} um",
        color="#193549",
        fontsize=12,
        pad=10,
    )
    figure.savefig(args.output / "layout_preview.png", dpi=220, facecolor="#dbe2ea")
    plt.close(figure)
    print(
        f"rendered {bounds.width:.6f} x {bounds.height:.6f} um, "
        f"aspect={bounds.width / bounds.height:.9f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
