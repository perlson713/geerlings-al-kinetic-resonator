#!/usr/bin/env python3
"""Optimize and build four-resonator chips for the two-cavity assembly."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon as PlotPolygon, Rectangle

from geerlings_resonator.chip_layout import (
    frequency_map,
    load_json,
    merge_current_cavity_field_calibration,
    optimize_rectangle_qc_positions,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _load_toml(path: Path) -> dict:
    with path.resolve().open("rb") as handle:
        return tomllib.load(handle)


def _project_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _load_bank_generator():
    script = PROJECT_ROOT / "scripts" / "generate_frequency_bank.py"
    spec = importlib.util.spec_from_file_location("generate_frequency_bank", script)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {script}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _build_pattern_layouts(
    bank_summary: dict,
    bank_config_path: Path,
) -> tuple[dict[str, object], dict[str, dict]]:
    generator = _load_bank_generator()
    bank_config = generator._load(bank_config_path)
    base_config = PROJECT_ROOT / bank_config["base_design_config"]
    layouts: dict[str, object] = {}
    resolved: dict[str, dict] = {}
    for row in bank_summary["patterns"]:
        pattern = str(row["pattern_id"])
        width_um = float(row["nominal_width_um"])
        evaluated = generator._evaluate(width_um, bank_config, base_config)
        layouts[pattern] = evaluated["layout"]
        resolved[pattern] = {
            "frequency_ghz": float(row["predicted_al200_ghz"]),
            "nominal_width_um": width_um,
            "cutout_width_um": float(evaluated["topology"]["cutout_width_um"]),
            "cutout_height_um": float(evaluated["topology"]["cutout_height_um"]),
            "fig1_idc_bottom": bool(evaluated["topology"]["fig1_idc_bottom"]),
            "full_width_top_return": bool(
                evaluated["topology"]["full_width_top_return"]
            ),
            "constant_radius_folds": bool(
                evaluated["topology"]["constant_radius_folds"]
            ),
        }
    return layouts, resolved


def _translated_points(points, center_x, center_y, target_x, target_y):
    return [
        (target_x + point.x - center_x, target_y + point.y - center_y)
        for point in points
    ]


def _ipolygon(kdb, points, database_unit):
    return kdb.DPolygon(
        [kdb.DPoint(x, y) for x, y in points]
    ).to_itype(database_unit)


def _write_chip_gds(
    output: Path,
    layouts: dict[str, object],
    placements: list[dict],
    *,
    chip_size_mm: float,
    chip_edge_clearance_mm: float,
    case_codes: dict[str, str],
    variant_tag: str = "",
) -> tuple[list[str], dict[str, dict]]:
    from klayout import db as kdb

    chip_size_um = chip_size_mm * 1000.0
    edge_um = chip_edge_clearance_mm * 1000.0
    artifacts: list[str] = []
    verification: dict[str, dict] = {}
    for cavity, cavity_code in case_codes.items():
        for chip in ("left", "right"):
            selected = [
                row
                for row in placements
                if row["cavity"] == cavity and row["chip"] == chip
            ]
            if len(selected) != 4:
                raise RuntimeError(f"expected four placements for {cavity}/{chip}")
            below = sorted(
                str(row["pattern_id"])
                for row in selected
                if row["frequency_band"] == "below_9ghz"
            )
            above = sorted(
                str(row["pattern_id"])
                for row in selected
                if row["frequency_band"] == "above_9ghz"
            )
            if len(below) + len(above) != 4:
                raise RuntimeError(f"invalid 9 GHz band labels for {cavity}/{chip}")
            band_name = "_".join(
                part
                for part in (
                    f"below9_{'-'.join(below)}" if below else "",
                    f"above9_{'-'.join(above)}" if above else "",
                )
                if part
            )
            variant_part = f"{variant_tag}_" if variant_tag else ""
            stem = f"cavity_{cavity_code}_{chip}_{variant_part}{band_name}_4res"
            library = kdb.Layout()
            library.dbu = 0.001
            cell = library.create_cell(
                f"CAV_{cavity_code}_{'L' if chip == 'left' else 'R'}_4RES"
            )
            metal_layer = library.layer(1, 0)
            outline_layer = library.layer(100, 0)
            label_layer = library.layer(101, 0)
            scale = 1.0 / library.dbu
            half = 0.5 * chip_size_um
            ground_half = half - edge_um
            ground = kdb.Region(
                kdb.Box(
                    round(-ground_half * scale),
                    round(-ground_half * scale),
                    round(ground_half * scale),
                    round(ground_half * scale),
                )
            )
            resonator_polygons = []
            for row in selected:
                target_x = 1000.0 * float(row["chip_local_x_mm"])
                target_y = 1000.0 * float(row["y_mm"])
                layout = layouts[row["pattern_id"]]
                cutout = next(
                    polygon
                    for polygon in layout.etch_polygons
                    if polygon.role == "ground_cutout"
                )
                center_x = 0.5 * (cutout.bounds.x_min + cutout.bounds.x_max)
                center_y = 0.5 * (cutout.bounds.y_min + cutout.bounds.y_max)
                ground -= kdb.Region(
                    _ipolygon(
                        kdb,
                        _translated_points(
                            cutout.points,
                            center_x,
                            center_y,
                            target_x,
                            target_y,
                        ),
                        library.dbu,
                    )
                )
                for polygon in layout.metal_polygons:
                    if polygon.role == "ground":
                        continue
                    resonator_polygons.append(
                        _ipolygon(
                            kdb,
                            _translated_points(
                                polygon.points,
                                center_x,
                                center_y,
                                target_x,
                                target_y,
                            ),
                            library.dbu,
                        )
                    )
                cell.shapes(label_layer).insert(
                    kdb.Text(
                        f"{row['pattern_id']}_{'BELOW9' if row['frequency_band'] == 'below_9ghz' else 'ABOVE9'}_ROT0",
                        round(target_x * scale),
                        round(target_y * scale),
                    )
                )
            cell.shapes(metal_layer).insert(ground)
            for polygon in resonator_polygons:
                cell.shapes(metal_layer).insert(polygon)
            cell.shapes(outline_layer).insert(
                kdb.Box(
                    round(-half * scale),
                    round(-half * scale),
                    round(half * scale),
                    round(half * scale),
                )
            )
            path = output / f"{stem}.gds"
            library.write(str(path))
            artifacts.append(path.name)

            readback = kdb.Layout()
            readback.read(str(path))
            top_cells = list(readback.top_cells())
            if len(top_cells) != 1:
                raise RuntimeError(f"{path.name}: expected one top cell")
            top = top_cells[0]
            bbox = top.dbbox()
            if abs(bbox.width() - chip_size_um) > readback.dbu:
                raise RuntimeError(f"{path.name}: unexpected chip width")
            if abs(bbox.height() - chip_size_um) > readback.dbu:
                raise RuntimeError(f"{path.name}: unexpected chip height")
            layer_counts = {}
            for layer_number, datatype in ((1, 0), (100, 0), (101, 0)):
                layer_index = readback.find_layer(layer_number, datatype)
                count = (
                    sum(1 for _ in top.each_shape(layer_index))
                    if layer_index is not None
                    else 0
                )
                layer_counts[f"{layer_number}/{datatype}"] = count
            if layer_counts["1/0"] < 1000:
                raise RuntimeError(f"{path.name}: incomplete metal geometry")
            if layer_counts["100/0"] != 1 or layer_counts["101/0"] != 4:
                raise RuntimeError(f"{path.name}: outline or label count mismatch")
            verification[path.name] = {
                "top_cells": 1,
                "bounding_box_um": [bbox.width(), bbox.height()],
                "shape_counts": layer_counts,
                "rotation_labels_all_zero": True,
                "below_9ghz_patterns": below,
                "above_9ghz_patterns": above,
            }
    return artifacts, verification


def _add_resonator_to_axis(
    axis,
    layout,
    row: dict,
    *,
    global_coordinates: bool,
) -> None:
    cutout = next(
        polygon for polygon in layout.etch_polygons if polygon.role == "ground_cutout"
    )
    center_x = 0.5 * (cutout.bounds.x_min + cutout.bounds.x_max)
    center_y = 0.5 * (cutout.bounds.y_min + cutout.bounds.y_max)
    x_mm = float(row["x_mm"] if global_coordinates else row["chip_local_x_mm"])
    y_mm = float(row["y_mm"])
    target_x = x_mm * 1000.0
    target_y = y_mm * 1000.0
    cutout_points = _translated_points(
        cutout.points, center_x, center_y, target_x, target_y
    )
    axis.add_patch(
        PlotPolygon(
            [(x / 1000.0, y / 1000.0) for x, y in cutout_points],
            closed=True,
            facecolor="#152530",
            edgecolor="#4ba3d2",
            linewidth=0.5,
        )
    )
    for polygon in layout.metal_polygons:
        if polygon.role == "ground":
            continue
        points = _translated_points(
            polygon.points, center_x, center_y, target_x, target_y
        )
        axis.add_patch(
            PlotPolygon(
                [(x / 1000.0, y / 1000.0) for x, y in points],
                closed=True,
                facecolor="#78c8f3",
                edgecolor="none",
            )
        )


def _render_layout_preview(
    path: Path,
    layouts: dict[str, object],
    placements: list[dict],
    config: dict,
    case_codes: dict[str, str],
    layout_description: str = "Qc-optimized rectangle corners",
) -> None:
    chip = config["chip"]
    cavity = config["cavity"]
    chip_size = float(chip["size_mm"])
    chip_center = float(chip["center_abs_x_mm"])
    straight = float(cavity["straight_length_mm"])
    width = float(cavity["width_mm"])
    radius = float(cavity["end_radius_mm"])
    total_length = straight + 2.0 * radius
    figure, axes = plt.subplots(2, 1, figsize=(14.2, 8.1), dpi=190)
    figure.subplots_adjust(
        left=0.075, right=0.985, bottom=0.08, top=0.95, hspace=0.28
    )
    for axis, (case, code) in zip(axes, case_codes.items()):
        axis.add_patch(
            FancyBboxPatch(
                (-0.5 * total_length, -0.5 * width),
                total_length,
                width,
                boxstyle=f"round,pad=0,rounding_size={radius}",
                fill=False,
                edgecolor="#315f7d",
                linewidth=1.5,
            )
        )
        for center in (-chip_center, chip_center):
            axis.add_patch(
                Rectangle(
                    (center - 0.5 * chip_size, -0.5 * chip_size),
                    chip_size,
                    chip_size,
                    facecolor="#8f72b5",
                    edgecolor="#78599b",
                    alpha=0.12,
                    linewidth=1.0,
                )
            )
        selected = [row for row in placements if row["cavity"] == case]
        for row in selected:
            _add_resonator_to_axis(
                axis, layouts[row["pattern_id"]], row, global_coordinates=True
            )
            is_upper = row["row"] == "upper"
            frequency_mark = "<9" if row["frequency_band"] == "below_9ghz" else ">9"
            axis.text(
                float(row["x_mm"]),
                float(row["y_mm"]) + (0.25 if is_upper else -0.25),
                f"{row['pattern_id']}\n{frequency_mark}",
                fontsize=6.2,
                ha="center",
                va="bottom" if is_upper else "top",
                linespacing=0.85,
            )
        axis.set_xlim(-0.5 * total_length - 0.5, 0.5 * total_length + 0.5)
        axis.set_ylim(-3.0, 3.0)
        axis.set_aspect("equal")
        axis.set_ylabel(f"Cavity {code}: y [mm]")
        axis.grid(alpha=0.12)
        axis.set_title(f"Cavity {code}: {layout_description}, below/above 9 GHz")
    axes[-1].set_xlabel("cavity long axis x [mm]")
    figure.savefig(path)
    plt.close(figure)


def _render_qc_preview(
    path: Path,
    summary: dict,
    nominal_target: float,
    case_codes: dict[str, str],
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(11.8, 4.6), dpi=190, layout="constrained")
    colors = {key: color for key, color in zip(case_codes, ("#2878a8", "#b25b24"))}
    for case, code in case_codes.items():
        rows = [row for row in summary["pin_calibration"] if row["cavity"] == case]
        axes[0].loglog(
            [row["pin_length_mm"] for row in rows],
            [row["target_resonator_external_qc"] for row in rows],
            "o-",
            color=colors[case],
            label=f"Cavity {code}",
        )
    axes[0].set_xlabel("SMA pin length [mm]")
    axes[0].set_ylabel("resonator-to-port external Qc")
    axes[0].grid(alpha=0.2, which="both")
    axes[0].legend(frameon=False)

    for case, code in case_codes.items():
        rows = [
            row
            for row in summary["qc_predictions"]
            if row["cavity"] == case
            and float(row["target_resonator_external_qc"]) == nominal_target
        ]
        rows.sort(key=lambda row: row["frequency_ghz"])
        axes[1].plot(
            [row["frequency_ghz"] for row in rows],
            [100.0 * (row["predicted_external_qc"] / nominal_target - 1.0) for row in rows],
            "o",
            color=colors[case],
            label=f"Cavity {code}",
        )
    axes[1].axhline(0.0, color="#555555", linestyle="--", linewidth=1.0)
    axes[1].set_xlabel("resonator frequency [GHz]")
    axes[1].set_ylabel(f"Qc error from {nominal_target:.0e} [%]")
    axes[1].grid(alpha=0.2)
    axes[1].legend(frameon=False)
    figure.savefig(path, bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=PROJECT_ROOT / "configs" / "four_resonator_chip.toml",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    config = _load_toml(args.config)
    bank_path = _project_path(config["frequency_bank_results"])
    bank_config_path = _project_path(config["frequency_bank_config"])
    calibration_path = _project_path(config["fem_calibration"])
    field_fem_path = _project_path(config["current_field_fem_results"])
    bank = load_json(bank_path)
    calibration = merge_current_cavity_field_calibration(
        load_json(calibration_path), load_json(field_fem_path)
    )
    chip = config["chip"]
    rectangle = config["rectangle"]
    centered_square = config["centered_square"]
    cavity = config["cavity"]
    optimization = config["optimization"]
    summary = optimize_rectangle_qc_positions(
        frequency_map(bank),
        calibration,
        config["sites"],
        target_qcs=optimization["target_resonator_external_qcs"],
        mean_coupling_mhz=float(optimization["mean_resonator_cavity_g_mhz"]),
        chip_center_abs_x_mm=float(chip["center_abs_x_mm"]),
        chip_size_mm=float(chip["size_mm"]),
        chip_edge_clearance_mm=float(chip["edge_clearance_mm"]),
        frequency_split_ghz=float(chip["frequency_split_ghz"]),
        coordinate_limit_mm=float(rectangle["coordinate_limit_mm"]),
        minimum_center_separation_mm=float(
            rectangle["minimum_center_separation_mm"]
        ),
        maximum_ground_cutout_width_mm=float(
            chip["maximum_ground_cutout_width_mm"]
        ),
        maximum_ground_cutout_height_mm=float(
            chip["maximum_ground_cutout_height_mm"]
        ),
        pin_wall_thickness_mm=float(cavity["pin_wall_thickness_mm"]),
        fixed_point_iterations=int(optimization["fixed_point_iterations"]),
        optimizer_seed=int(optimization["optimizer_seed"]),
        optimizer_maxiter=int(optimization["optimizer_maxiter"]),
        optimizer_popsize=int(optimization["optimizer_popsize"]),
        position_tolerances_mm=[
            float(value) / 1000.0
            for value in optimization["position_tolerances_um"]
        ],
    )
    square_summary = optimize_rectangle_qc_positions(
        frequency_map(bank),
        calibration,
        config["sites"],
        target_qcs=optimization["target_resonator_external_qcs"],
        mean_coupling_mhz=float(optimization["mean_resonator_cavity_g_mhz"]),
        chip_center_abs_x_mm=float(chip["center_abs_x_mm"]),
        chip_size_mm=float(chip["size_mm"]),
        chip_edge_clearance_mm=float(chip["edge_clearance_mm"]),
        frequency_split_ghz=float(chip["frequency_split_ghz"]),
        coordinate_limit_mm=float(rectangle["coordinate_limit_mm"]),
        minimum_center_separation_mm=float(
            rectangle["minimum_center_separation_mm"]
        ),
        maximum_ground_cutout_width_mm=float(
            chip["maximum_ground_cutout_width_mm"]
        ),
        maximum_ground_cutout_height_mm=float(
            chip["maximum_ground_cutout_height_mm"]
        ),
        pin_wall_thickness_mm=float(cavity["pin_wall_thickness_mm"]),
        fixed_point_iterations=int(optimization["fixed_point_iterations"]),
        optimizer_seed=int(optimization["optimizer_seed"]),
        optimizer_maxiter=int(optimization["optimizer_maxiter"]),
        optimizer_popsize=int(optimization["optimizer_popsize"]),
        centered_square_side_mm=float(centered_square["side_mm"]),
        position_tolerances_mm=[
            float(value) / 1000.0
            for value in optimization["position_tolerances_um"]
        ],
    )
    output = args.output or _project_path(config["output_directory"])
    if not output.is_absolute():
        output = PROJECT_ROOT / output
    output.mkdir(parents=True, exist_ok=True)

    case_codes = {
        str(cavity["case_a"]): "A",
        str(cavity["case_b"]): "B",
    }
    if set(case_codes) != {
        row["cavity"] for row in summary["placements"]
    }:
        raise ValueError("configured cavity cases do not match FEM calibration")
    if set(case_codes) != {
        row["cavity"] for row in square_summary["placements"]
    }:
        raise ValueError("centered-square cavity cases do not match FEM calibration")

    layouts, patterns = _build_pattern_layouts(bank, bank_config_path)
    gds_artifacts, verification = _write_chip_gds(
        output,
        layouts,
        summary["placements"],
        chip_size_mm=float(chip["size_mm"]),
        chip_edge_clearance_mm=float(chip["edge_clearance_mm"]),
        case_codes=case_codes,
    )
    square_gds_artifacts, square_verification = _write_chip_gds(
        output,
        layouts,
        square_summary["placements"],
        chip_size_mm=float(chip["size_mm"]),
        chip_edge_clearance_mm=float(chip["edge_clearance_mm"]),
        case_codes=case_codes,
        variant_tag="centered_square",
    )
    layout_preview = "four_per_chip_layout_preview.png"
    qc_preview = "qc_uniformity.png"
    square_layout_preview = "centered_square_layout_preview.png"
    square_qc_preview = "centered_square_qc_uniformity.png"
    _render_layout_preview(
        output / layout_preview,
        layouts,
        summary["placements"],
        config,
        case_codes,
    )
    _render_layout_preview(
        output / square_layout_preview,
        layouts,
        square_summary["placements"],
        config,
        case_codes,
        layout_description=(
            f"centered {float(centered_square['side_mm']):.2f} mm square"
        ),
    )
    _render_qc_preview(
        output / qc_preview,
        summary,
        float(optimization["nominal_plot_target_qc"]),
        case_codes,
    )
    _render_qc_preview(
        output / square_qc_preview,
        square_summary,
        float(optimization["nominal_plot_target_qc"]),
        case_codes,
    )

    _write_csv(output / "placements.csv", summary["placements"])
    _write_csv(output / "pin_calibration.csv", summary["pin_calibration"])
    _write_csv(output / "qc_predictions.csv", summary["qc_predictions"])
    _write_csv(
        output / "centered_square_placements.csv", square_summary["placements"]
    )
    _write_csv(
        output / "centered_square_pin_calibration.csv",
        square_summary["pin_calibration"],
    )
    _write_csv(
        output / "centered_square_qc_predictions.csv",
        square_summary["qc_predictions"],
    )
    summary["inputs"] = {
        "config": args.config.resolve().relative_to(PROJECT_ROOT).as_posix(),
        "frequency_bank_results": bank_path.relative_to(PROJECT_ROOT).as_posix(),
        "fem_calibration": calibration_path.relative_to(PROJECT_ROOT).as_posix(),
        "current_field_fem_results": field_fem_path.relative_to(PROJECT_ROOT).as_posix(),
    }
    summary["artifacts"] = gds_artifacts + [
        layout_preview,
        qc_preview,
        "placements.csv",
        "pin_calibration.csv",
        "qc_predictions.csv",
        "current_cavity_field_fem.json",
        "optimization_summary.json",
        "layout_manifest.json",
        "gds_readback_verification.json",
    ]
    (output / "optimization_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "gds_readback_verification.json").write_text(
        json.dumps(verification, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest = {
        "schema_id": "geerlings.four_resonator_chip.gds.v1",
        "layout_variant": "qc_optimized_rectangle",
        "chip_size_mm": [float(chip["size_mm"]), float(chip["size_mm"])],
        "resonators_per_chip": 4,
        "chips_per_cavity": 2,
        "total_resonators": 16,
        "rotation_deg": 0.0,
        "centers_form_axis_aligned_rectangle": True,
        "frequency_partition": summary["constraints"]["frequency_partition"],
        "patterns": patterns,
        "rectangle_geometry": summary["rectangle_geometry"],
        "placements": summary["placements"],
        "pin_calibration": summary["pin_calibration"],
        "gds_layers": {
            "positive_aluminium_metal": [1, 0],
            "chip_outline": [100, 0],
            "labels": [101, 0],
        },
        "gds_readback": verification,
        "artifacts": summary["artifacts"],
    }
    (output / "layout_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    square_summary["inputs"] = dict(summary["inputs"])
    square_summary["artifacts"] = square_gds_artifacts + [
        square_layout_preview,
        square_qc_preview,
        "centered_square_placements.csv",
        "centered_square_pin_calibration.csv",
        "centered_square_qc_predictions.csv",
        "current_cavity_field_fem.json",
        "centered_square_optimization_summary.json",
        "centered_square_layout_manifest.json",
        "centered_square_gds_readback_verification.json",
    ]
    (output / "centered_square_optimization_summary.json").write_text(
        json.dumps(square_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (output / "centered_square_gds_readback_verification.json").write_text(
        json.dumps(square_verification, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    square_manifest = {
        "schema_id": "geerlings.four_resonator_chip.gds.v1",
        "layout_variant": "centered_square",
        "chip_size_mm": [float(chip["size_mm"]), float(chip["size_mm"])],
        "resonators_per_chip": 4,
        "chips_per_cavity": 2,
        "total_resonators": 16,
        "rotation_deg": 0.0,
        "centers_form_axis_aligned_rectangle": True,
        "centered_on_chip": True,
        "square_side_mm": float(centered_square["side_mm"]),
        "frequency_partition": square_summary["constraints"]["frequency_partition"],
        "patterns": patterns,
        "rectangle_geometry": square_summary["rectangle_geometry"],
        "placements": square_summary["placements"],
        "pin_calibration": square_summary["pin_calibration"],
        "gds_layers": {
            "positive_aluminium_metal": [1, 0],
            "chip_outline": [100, 0],
            "labels": [101, 0],
        },
        "gds_readback": square_verification,
        "artifacts": square_summary["artifacts"],
    }
    (output / "centered_square_layout_manifest.json").write_text(
        json.dumps(square_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(output),
                "gds": gds_artifacts,
                "centered_square_gds": square_gds_artifacts,
                "global_qc_max_to_min": summary[
                    "global_16_resonator_qc_max_to_min"
                ],
                "minimum_gap_mm": summary[
                    "minimum_nominal_footprint_edge_gap_mm"
                ],
                "centered_square_global_qc_max_to_min": square_summary[
                    "global_16_resonator_qc_max_to_min"
                ],
                "centered_square_minimum_gap_mm": square_summary[
                    "minimum_nominal_footprint_edge_gap_mm"
                ],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
