#!/usr/bin/env python3
"""Generate the cavity-C one-cell, two-layer four-resonator masks."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as PlotPolygon, Rectangle

from geerlings_resonator.config import load_project
from geerlings_resonator.geometry import generate_layout


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "four_resonator_chip.toml"
DEFAULT_PLACEMENTS = (
    PROJECT_ROOT
    / "results"
    / "four_resonator_chip"
    / "centered_square_placements.csv"
)


def _load_toml(path: Path) -> dict:
    with path.resolve().open("rb") as handle:
        return tomllib.load(handle)


def _load_json(path: Path) -> dict:
    return json.loads(path.resolve().read_text(encoding="utf-8"))


def _project_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def _load_cavity_c_placements(path: Path) -> list[dict[str, str]]:
    with path.resolve().open(newline="", encoding="utf-8") as handle:
        rows = [
            row
            for row in csv.DictReader(handle)
            if row["cavity"] == "cavity_C"
        ]
    if len(rows) != 8:
        raise RuntimeError("expected eight cavity-C placement rows")
    for chip in ("left", "right"):
        if sum(row["chip"] == chip for row in rows) != 4:
            raise RuntimeError(f"expected four cavity-C/{chip} placement rows")
    return rows


def _build_pattern_layouts(config: dict, pattern_ids: set[str]) -> dict[str, object]:
    bank_path = _project_path(config["frequency_bank_results"])
    bank_config_path = _project_path(config["frequency_bank_config"])
    bank = _load_json(bank_path)
    bank_config = _load_toml(bank_config_path)
    base_config = _project_path(bank_config["base_design_config"])
    rows = {
        str(row["pattern_id"]): row
        for row in bank["patterns"]
        if str(row["pattern_id"]) in pattern_ids
    }
    if set(rows) != pattern_ids:
        raise RuntimeError("frequency-bank patterns do not match cavity-C placements")

    layouts = {}
    for pattern_id, row in rows.items():
        project = load_project(
            base_config,
            overrides=[
                f"resonator.nominal_width_um={float(row['nominal_width_um']):.12f}"
            ],
        )
        layout = generate_layout(project.resonator)
        if layout.centerlines.feedline is not None or layout.ports:
            raise RuntimeError(f"{pattern_id}: expected an isolated resonator")
        layouts[pattern_id] = layout
    return layouts


def _translated_points(points, center_x, center_y, target_x, target_y):
    return [
        (target_x + point.x - center_x, target_y + point.y - center_y)
        for point in points
    ]


def _ipolygon(kdb, points, database_unit):
    return kdb.DPolygon(
        [kdb.DPoint(x, y) for x, y in points]
    ).to_itype(database_unit)


def _corner_l_points(
    corner_x_um: float,
    corner_y_um: float,
    *,
    outer_um: float,
    cutout_um: float,
) -> list[tuple[float, float]]:
    """Return outer-square minus inward 75-square as one six-vertex L."""

    leg_um = outer_um - cutout_um
    inward_x = -math.copysign(1.0, corner_x_um)
    inward_y = -math.copysign(1.0, corner_y_um)
    local_points = (
        (0.0, 0.0),
        (outer_um, 0.0),
        (outer_um, leg_um),
        (leg_um, leg_um),
        (leg_um, outer_um),
        (0.0, outer_um),
    )
    return [
        (
            corner_x_um + inward_x * local_x,
            corner_y_um + inward_y * local_y,
        )
        for local_x, local_y in local_points
    ]


def _band_filename(chip: str, selected: list[dict[str, str]]) -> str:
    below = sorted(
        row["pattern_id"]
        for row in selected
        if row["frequency_band"] == "below_9ghz"
    )
    above = sorted(
        row["pattern_id"]
        for row in selected
        if row["frequency_band"] == "above_9ghz"
    )
    if len(below) != 2 or len(above) != 2:
        raise RuntimeError(f"cavity-C/{chip}: expected two patterns in each band")
    return (
        f"cavity_C_{chip}_centered_square_"
        f"below9_{'-'.join(below)}_above9_{'-'.join(above)}_4res.gds"
    )


def _verify_gds(
    path: Path,
    *,
    chip_size_um: float,
    marker_outer_um: float,
    marker_cutout_um: float,
) -> dict:
    from klayout import db as kdb

    readback = kdb.Layout()
    readback.read(str(path))
    top_cells = list(readback.top_cells())
    if readback.cells() != 1 or len(top_cells) != 1:
        raise RuntimeError(f"{path.name}: expected exactly one cell")
    top = top_cells[0]
    if list(top.each_inst()):
        raise RuntimeError(f"{path.name}: hierarchy is not flat")
    populated_layers = {
        (info.layer, info.datatype)
        for info in readback.layer_infos()
        if not top.shapes(readback.find_layer(info.layer, info.datatype)).is_empty()
    }
    if populated_layers != {(1, 0), (2, 0)}:
        raise RuntimeError(f"{path.name}: expected only layers 1/0 and 2/0")

    resonator_layer = readback.find_layer(1, 0)
    marker_layer = readback.find_layer(2, 0)
    resonator_shapes = list(top.each_shape(resonator_layer))
    marker_shapes = list(top.each_shape(marker_layer))
    if not resonator_shapes or any(shape.is_text() for shape in resonator_shapes):
        raise RuntimeError(f"{path.name}: invalid resonator geometry")
    if len(marker_shapes) != 4 or any(not shape.is_polygon() for shape in marker_shapes):
        raise RuntimeError(f"{path.name}: expected four polygonal L marks")

    expected_marker_area = marker_outer_um**2 - marker_cutout_um**2
    marker_areas = []
    marker_bounding_boxes = []
    for shape in marker_shapes:
        polygon = shape.polygon
        bbox = polygon.bbox()
        width_um = bbox.width() * readback.dbu
        height_um = bbox.height() * readback.dbu
        area_um2 = polygon.area() * readback.dbu**2
        if not math.isclose(width_um, marker_outer_um, abs_tol=readback.dbu):
            raise RuntimeError(f"{path.name}: L-mark width mismatch")
        if not math.isclose(height_um, marker_outer_um, abs_tol=readback.dbu):
            raise RuntimeError(f"{path.name}: L-mark height mismatch")
        if not math.isclose(area_um2, expected_marker_area, abs_tol=1.0e-6):
            raise RuntimeError(f"{path.name}: L-mark area mismatch")
        marker_areas.append(area_um2)
        marker_bounding_boxes.append([width_um, height_um])

    bbox = top.dbbox()
    if not math.isclose(bbox.width(), chip_size_um, abs_tol=readback.dbu):
        raise RuntimeError(f"{path.name}: chip-width extent mismatch")
    if not math.isclose(bbox.height(), chip_size_um, abs_tol=readback.dbu):
        raise RuntimeError(f"{path.name}: chip-height extent mismatch")
    return {
        "cell_count": readback.cells(),
        "top_cell_count": len(top_cells),
        "instance_count": 0,
        "populated_layers": ["1/0", "2/0"],
        "shape_counts": {
            "1/0_resonator": len(resonator_shapes),
            "2/0_corner_L": len(marker_shapes),
        },
        "bounding_box_um": [bbox.width(), bbox.height()],
        "corner_L_bounding_boxes_um": marker_bounding_boxes,
        "corner_L_areas_um2": marker_areas,
    }


def _coupling_screen(spacing_mm: float) -> dict | None:
    path = (
        PROJECT_ROOT
        / "results"
        / "four_resonator_chip"
        / "centered_square_coupling_sweep.csv"
    )
    if not path.exists():
        return None
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected = min(
        rows,
        key=lambda row: abs(float(row["square_side_mm"]) - spacing_mm),
    )
    if not math.isclose(
        float(selected["square_side_mm"]), spacing_mm, abs_tol=1.0e-9
    ):
        return None
    return {
        "method": "existing reduced-order electric-plus-magnetic dipole sweep",
        "maximum_pair_coupling_khz": float(
            selected["maximum_pair_coupling_khz"]
        ),
        "maximum_mixing_amplitude_percent": float(
            selected["maximum_mixing_amplitude_percent"]
        ),
        "maximum_dispersive_shift_khz": float(
            selected["maximum_dispersive_shift_khz"]
        ),
        "worst_pair": selected["worst_pair"],
        "passes_previous_negligible_coupling_thresholds": False,
    }


def _shared_readback_entry(details: dict, payload: dict) -> dict:
    placements = details["placements"]
    return {
        "cell_count": details["cell_count"],
        "top_cells": details["top_cell_count"],
        "instance_count": details["instance_count"],
        "bounding_box_um": details["bounding_box_um"],
        "populated_layers": details["populated_layers"],
        "shape_counts": {
            "1/0": details["shape_counts"]["1/0_resonator"],
            "2/0": details["shape_counts"]["2/0_corner_L"],
        },
        "below_9ghz_patterns": sorted(
            row["pattern_id"]
            for row in placements
            if row["frequency_band"] == "below_9ghz"
        ),
        "above_9ghz_patterns": sorted(
            row["pattern_id"]
            for row in placements
            if row["frequency_band"] == "above_9ghz"
        ),
        "resonator_center_spacing_mm": payload["resonator_center_spacing_mm"],
        "corner_L_outer_square_um": payload["corner_L"]["outer_square_um"],
        "corner_L_cutout_square_um": payload["corner_L"]["cutout_square_um"],
    }


def _update_existing_metadata(output: Path, payload: dict) -> None:
    readback_path = output / "centered_square_gds_readback_verification.json"
    shared_readback = None
    if readback_path.exists():
        shared_readback = _load_json(readback_path)
        for filename, details in payload["verification"].items():
            shared_readback[filename] = _shared_readback_entry(details, payload)
        readback_path.write_text(
            json.dumps(shared_readback, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    artifact_names = [payload["preview"], payload["verification_artifact"]]
    override = {
        "cavity": "cavity_C",
        "cavity_height_mm": 15.2,
        "gds_files": sorted(payload["verification"]),
        "cell_count_per_gds": 1,
        "instance_count_per_gds": 0,
        "populated_layers": payload["layers"],
        "square_side_and_nearest_center_spacing_mm": payload[
            "resonator_center_spacing_mm"
        ],
        "resonator_centers_mm": payload["resonator_centers_mm"],
        "corner_L": payload["corner_L"],
        "preview": payload["preview"],
        "verification": payload["verification_artifact"],
        "analysis_scope": payload["analysis_scope"],
        "prior_reduced_order_coupling_screen": payload[
            "prior_reduced_order_coupling_screen"
        ],
    }

    manifest_path = output / "centered_square_layout_manifest.json"
    if manifest_path.exists():
        manifest = _load_json(manifest_path)
        manifest["cavity_C_two_layer_gds_override"] = override
        if shared_readback is not None:
            manifest["gds_readback"] = shared_readback
        manifest["artifacts"] = list(
            dict.fromkeys([*manifest.get("artifacts", []), *artifact_names])
        )
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    summary_path = output / "centered_square_optimization_summary.json"
    if summary_path.exists():
        summary = _load_json(summary_path)
        summary["cavity_C_two_layer_gds_override"] = override
        summary["artifacts"] = list(
            dict.fromkeys([*summary.get("artifacts", []), *artifact_names])
        )
        summary_path.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )


def _draw_preview(
    path: Path,
    masks: list[dict],
    *,
    chip_size_um: float,
    spacing_mm: float,
    marker_outer_um: float,
    marker_cutout_um: float,
) -> None:
    figure = plt.figure(figsize=(13.0, 5.2), dpi=210, layout="constrained")
    grid = figure.add_gridspec(1, 3, width_ratios=(1.0, 1.0, 0.68))
    layer_colors = {1: "#145a7a", 2: "#ef8a24"}
    half_mm = 0.5 * chip_size_um / 1000.0

    for index, mask in enumerate(masks):
        axis = figure.add_subplot(grid[0, index])
        axis.add_patch(
            Rectangle(
                (-half_mm, -half_mm),
                2.0 * half_mm,
                2.0 * half_mm,
                facecolor="none",
                edgecolor="#777777",
                linestyle="--",
                linewidth=0.9,
            )
        )
        for polygon in mask["resonator_polygons_um"]:
            axis.add_patch(
                PlotPolygon(
                    [(x / 1000.0, y / 1000.0) for x, y in polygon],
                    closed=True,
                    facecolor=layer_colors[1],
                    edgecolor="none",
                )
            )
        for polygon in mask["marker_polygons_um"]:
            axis.add_patch(
                PlotPolygon(
                    [(x / 1000.0, y / 1000.0) for x, y in polygon],
                    closed=True,
                    facecolor=layer_colors[2],
                    edgecolor="none",
                )
            )
        half_spacing = 0.5 * spacing_mm
        dimension_y = half_spacing + 0.55
        axis.annotate(
            "",
            xy=(-half_spacing, dimension_y),
            xytext=(half_spacing, dimension_y),
            arrowprops={"arrowstyle": "<->", "color": "#333333", "lw": 1.0},
        )
        axis.text(
            0.0,
            dimension_y + 0.08,
            f"{spacing_mm:.2f} mm center-to-center",
            ha="center",
            va="bottom",
            fontsize=8,
        )
        axis.set_title(f"cavity C — {mask['chip']} chip")
        axis.set_xlabel("x [mm]")
        if index == 0:
            axis.set_ylabel("y [mm]")
        axis.set_xlim(-half_mm - 0.12, half_mm + 0.12)
        axis.set_ylim(-half_mm - 0.12, half_mm + 0.12)
        axis.set_aspect("equal")
        axis.grid(alpha=0.12)

    detail = figure.add_subplot(grid[0, 2])
    leg_um = marker_outer_um - marker_cutout_um
    canonical = (
        (0.0, 0.0),
        (marker_outer_um, 0.0),
        (marker_outer_um, leg_um),
        (leg_um, leg_um),
        (leg_um, marker_outer_um),
        (0.0, marker_outer_um),
    )
    detail.add_patch(
        PlotPolygon(
            canonical,
            closed=True,
            facecolor=layer_colors[2],
            edgecolor="#8a4b10",
            linewidth=0.9,
        )
    )
    detail.add_patch(
        Rectangle(
            (leg_um, leg_um),
            marker_cutout_um,
            marker_cutout_um,
            facecolor="none",
            edgecolor="#555555",
            linestyle="--",
            linewidth=0.9,
        )
    )
    detail.text(
        leg_um + 0.5 * marker_cutout_um,
        leg_um + 0.5 * marker_cutout_um,
        f"{marker_cutout_um:.0f} × {marker_cutout_um:.0f}\ncut-out",
        ha="center",
        va="center",
        fontsize=8,
    )
    detail.set_title("Layer 2 corner-L detail")
    detail.set_xlabel("µm")
    detail.set_ylabel("µm")
    detail.set_xlim(-8.0, marker_outer_um + 8.0)
    detail.set_ylim(-8.0, marker_outer_um + 8.0)
    detail.set_aspect("equal")
    detail.grid(alpha=0.15)
    figure.suptitle(
        "Cavity-C preview: one flat cell, layer 1 resonators, layer 2 corner L marks",
        fontsize=12,
    )
    figure.savefig(path)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--placements", type=Path, default=DEFAULT_PLACEMENTS)
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=PROJECT_ROOT / "build" / "cavity_c_two_layer_preview",
    )
    parser.add_argument("--spacing-mm", type=float)
    parser.add_argument("--marker-outer-um", type=float)
    parser.add_argument("--marker-cutout-um", type=float)
    args = parser.parse_args()

    config = _load_toml(args.config)
    mask_config = config["cavity_c_two_layer_gds"]
    spacing_mm = float(
        mask_config["spacing_mm"] if args.spacing_mm is None else args.spacing_mm
    )
    marker_outer_um = float(
        mask_config["marker_outer_um"]
        if args.marker_outer_um is None
        else args.marker_outer_um
    )
    marker_cutout_um = float(
        mask_config["marker_cutout_um"]
        if args.marker_cutout_um is None
        else args.marker_cutout_um
    )
    if int(mask_config["resonator_layer"]) != 1 or int(
        mask_config["corner_l_layer"]
    ) != 2:
        raise ValueError("cavity-C mask layers must remain 1 and 2")
    if spacing_mm <= 0.0:
        raise ValueError("spacing must be positive")
    if not 0.0 < marker_cutout_um < marker_outer_um:
        raise ValueError("marker cutout must be positive and smaller than its outer square")

    placements = _load_cavity_c_placements(args.placements)
    pattern_ids = {row["pattern_id"] for row in placements}
    layouts = _build_pattern_layouts(config, pattern_ids)
    output = args.output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)
    chip_size_um = float(config["chip"]["size_mm"]) * 1000.0
    half_chip_um = 0.5 * chip_size_um
    half_spacing_um = 0.5 * spacing_mm * 1000.0

    from klayout import db as kdb

    masks = []
    verification = {}
    for chip in ("left", "right"):
        selected = [row for row in placements if row["chip"] == chip]
        filename = _band_filename(chip, selected)
        library = kdb.Layout()
        library.dbu = 0.001
        cell = library.create_cell(f"CAV_C_{'L' if chip == 'left' else 'R'}_4RES_2L")
        resonator_layer = library.layer(1, 0)
        marker_layer = library.layer(2, 0)
        resonator_polygons_um = []
        resolved_placements = []

        for row in selected:
            source_x = float(row["chip_local_x_mm"])
            source_y = float(row["y_mm"])
            target_x = math.copysign(half_spacing_um, source_x)
            target_y = math.copysign(half_spacing_um, source_y)
            layout = layouts[row["pattern_id"]]
            cutout = next(
                polygon
                for polygon in layout.etch_polygons
                if polygon.role == "ground_cutout"
            )
            center_x = 0.5 * (cutout.bounds.x_min + cutout.bounds.x_max)
            center_y = 0.5 * (cutout.bounds.y_min + cutout.bounds.y_max)
            for polygon in layout.metal_polygons:
                if polygon.role == "ground":
                    continue
                points = _translated_points(
                    polygon.points,
                    center_x,
                    center_y,
                    target_x,
                    target_y,
                )
                cell.shapes(resonator_layer).insert(
                    _ipolygon(kdb, points, library.dbu)
                )
                resonator_polygons_um.append(points)
            resolved_placements.append(
                {
                    "pattern_id": row["pattern_id"],
                    "frequency_band": row["frequency_band"],
                    "center_x_mm": target_x / 1000.0,
                    "center_y_mm": target_y / 1000.0,
                }
            )

        marker_polygons_um = []
        for sign_x, sign_y in ((-1, -1), (-1, 1), (1, -1), (1, 1)):
            points = _corner_l_points(
                sign_x * half_chip_um,
                sign_y * half_chip_um,
                outer_um=marker_outer_um,
                cutout_um=marker_cutout_um,
            )
            cell.shapes(marker_layer).insert(_ipolygon(kdb, points, library.dbu))
            marker_polygons_um.append(points)

        path = output / filename
        save_options = kdb.SaveLayoutOptions()
        save_options.gds2_write_timestamps = False
        library.write(str(path), save_options)
        verification[filename] = {
            **_verify_gds(
                path,
                chip_size_um=chip_size_um,
                marker_outer_um=marker_outer_um,
                marker_cutout_um=marker_cutout_um,
            ),
            "placements": resolved_placements,
        }
        masks.append(
            {
                "chip": chip,
                "filename": filename,
                "resonator_polygons_um": resonator_polygons_um,
                "marker_polygons_um": marker_polygons_um,
            }
        )

    preview_name = "cavity_C_two_layer_1p5mm_preview.png"
    _draw_preview(
        output / preview_name,
        masks,
        chip_size_um=chip_size_um,
        spacing_mm=spacing_mm,
        marker_outer_um=marker_outer_um,
        marker_cutout_um=marker_cutout_um,
    )
    payload = {
        "schema_id": "geerlings.cavity_C.two_layer_mask.v1",
        "cavity": "cavity_C",
        "cell_structure": "one flat cell per GDS; no instances",
        "layers": {
            "1/0": "four isolated resonator metal patterns",
            "2/0": "four corner L marks",
        },
        "resonator_center_spacing_mm": spacing_mm,
        "resonator_centers_mm": [
            [-0.5 * spacing_mm, -0.5 * spacing_mm],
            [-0.5 * spacing_mm, 0.5 * spacing_mm],
            [0.5 * spacing_mm, -0.5 * spacing_mm],
            [0.5 * spacing_mm, 0.5 * spacing_mm],
        ],
        "corner_L": {
            "outer_square_um": [marker_outer_um, marker_outer_um],
            "cutout_square_um": [marker_cutout_um, marker_cutout_um],
            "leg_width_um": marker_outer_um - marker_cutout_um,
            "count": 4,
        },
        "preview": preview_name,
        "verification_artifact": "cavity_C_two_layer_1p5mm_verification.json",
        "analysis_scope": {
            "kind": "mask-only geometry override",
            "qc_recomputed": False,
            "direct_coupling_recomputed_fullwave": False,
            "warning": (
                "the 4.00 mm centered-square Qc and negligible-coupling results "
                "do not describe these 1.50 mm cavity-C GDS files"
            ),
        },
        "prior_reduced_order_coupling_screen": _coupling_screen(spacing_mm),
        "verification": verification,
    }
    verification_name = payload["verification_artifact"]
    (output / verification_name).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _update_existing_metadata(output, payload)
    print(
        json.dumps(
            {
                "output": str(output),
                "gds": [mask["filename"] for mask in masks],
                "preview": preview_name,
                "verification": verification_name,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
