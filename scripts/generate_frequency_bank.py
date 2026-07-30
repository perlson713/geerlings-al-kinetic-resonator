#!/usr/bin/env python3
"""Generate eight Fig. 1-topology layouts on a nominal 10 MHz frequency grid."""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import json
import math
from pathlib import Path
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as PlotPolygon, Rectangle

from geerlings_resonator.config import load_project
from geerlings_resonator.em import write_stackup_xml
from geerlings_resonator.gds import write_gds
from geerlings_resonator.geometry import generate_layout
from geerlings_resonator.kinetic import (
    AluminiumLondonModel,
    ReducedOrderResonator,
    evaluate_point,
)
from geerlings_resonator.svg import write_svg


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def symmetric_targets(center_ghz: float, spacing_mhz: float, count: int) -> list[float]:
    """Return an exactly centered frequency grid, including for an even count."""

    if count < 1:
        raise ValueError("count must be positive")
    if spacing_mhz <= 0:
        raise ValueError("spacing_mhz must be positive")
    step_ghz = spacing_mhz * 1.0e-3
    midpoint = 0.5 * (count - 1)
    return [center_ghz + (index - midpoint) * step_ghz for index in range(count)]


def _load(path: Path) -> dict:
    with path.resolve().open("rb") as handle:
        return tomllib.load(handle)


def _pec_frequency(width_um: float, calibration: dict) -> float:
    width_delta = float(calibration["second_width_um"]) - float(
        calibration["anchor_width_um"]
    )
    if width_delta == 0:
        raise ValueError("PEC calibration widths must differ")
    slope = (
        float(calibration["second_pec_frequency_ghz"])
        - float(calibration["anchor_pec_frequency_ghz"])
    ) / width_delta
    return float(calibration["anchor_pec_frequency_ghz"]) + slope * (
        width_um - float(calibration["anchor_width_um"])
    )


def _topology_guard(layout) -> dict:
    parameters = layout.parameters
    if parameters.inductor_turns != 6 or parameters.inductor_turns % 2:
        raise RuntimeError("frequency bank requires the approved even six-turn meander")
    if parameters.capacitor_finger_count != 14:
        raise RuntimeError("frequency bank requires the approved 14-finger IDC")
    if layout.centerlines.feedline is not None or layout.ports:
        raise RuntimeError("frequency bank requires the isolated no-port layout")

    fingers = layout.centerlines.capacitor_fingers
    left_bus, right_bus = layout.centerlines.capacitor_buses
    if fingers[0].name != "capacitor_finger_00_left":
        raise RuntimeError("lowest IDC finger must remain connected to the left bus")
    if fingers[1].name != "capacitor_finger_01_right":
        raise RuntimeError("second IDC finger must remain connected to the right bus")
    if not math.isclose(fingers[0].points[0].y, 2.5, abs_tol=1.0e-12):
        raise RuntimeError("lowest IDC centerline must remain at y=2.5 um")
    if not math.isclose(fingers[1].points[0].y, 17.5, abs_tol=1.0e-12):
        raise RuntimeError("first right-bus finger must remain one pitch higher")
    if not math.isclose(fingers[0].points[0].x, left_bus.points[0].x, abs_tol=1.0e-12):
        raise RuntimeError("lowest IDC finger is detached from the left bus")
    if not math.isclose(fingers[1].points[0].x, right_bus.points[0].x, abs_tol=1.0e-12):
        raise RuntimeError("second IDC finger is detached from the right bus")
    if left_bus.points[0].y != 0.0 or right_bus.points[0].y != 0.0:
        raise RuntimeError("IDC buses must retain their common y=0 lower edge")

    inductor = layout.centerlines.inductor
    if len(inductor.points) != 158:
        raise RuntimeError("approved six-turn inductor must have 158 sampled points")
    left_points = inductor.points[:79]
    right_points = inductor.points[79:]
    for left, mirrored_right in zip(left_points, reversed(right_points)):
        if not math.isclose(left.x, -mirrored_right.x, abs_tol=1.0e-9):
            raise RuntimeError("meander lost left/right mirror symmetry")
        if not math.isclose(left.y, mirrored_right.y, abs_tol=1.0e-9):
            raise RuntimeError("meander lost left/right mirror symmetry")
    for row in range(6):
        radius_um = 12.5
        center_x = (
            -25.0
            if row % 2 == 0
            else -(parameters.nominal_width_um / 2.0 - 25.0)
        )
        center_y = 235.0 + row * 25.0
        turn_start = 1 + row * 13
        for point in left_points[turn_start : turn_start + 13]:
            radius = math.hypot(point.x - center_x, point.y - center_y)
            if not math.isclose(radius, radius_um, abs_tol=1.0e-9):
                raise RuntimeError("meander fold radius changed")
    left_top, right_top = inductor.points[78], inductor.points[79]
    if not math.isclose(left_top.y, right_top.y, abs_tol=1.0e-12):
        raise RuntimeError("top return is not horizontal")
    if not math.isclose(left_top.x, -right_top.x, abs_tol=1.0e-12):
        raise RuntimeError("top return is not full-width and symmetric")

    cutout = next(poly for poly in layout.etch_polygons if poly.role == "ground_cutout")
    return {
        "fig1_idc_bottom": True,
        "full_width_top_return": True,
        "constant_radius_folds": True,
        "left_right_mirror_symmetry": True,
        "cutout_width_um": cutout.bounds.width,
        "cutout_height_um": cutout.bounds.height,
    }


def _evaluate(width_um: float, data: dict, base_config: Path) -> dict:
    project = load_project(
        base_config,
        overrides=[f"resonator.nominal_width_um={width_um:.12f}"],
    )
    layout = generate_layout(project.resonator)
    topology = _topology_guard(layout)
    calibration = data["calibration"]
    trace_width_um = layout.centerlines.inductor.width_um
    centerline_length_um = layout.centerlines.inductor.length
    branch_squares = (
        centerline_length_um
        / trace_width_um
        * float(calibration["branch_square_scale"])
    )
    pec_frequency_ghz = _pec_frequency(width_um, calibration)
    resonator = ReducedOrderResonator(
        pec_frequency_ghz=pec_frequency_ghz,
        modal_impedance_ohm=float(calibration["modal_impedance_ohm"]),
        branch_squares=branch_squares,
        meander_energy_fraction=float(calibration["meander_energy_fraction"]),
    )
    aluminium = data["aluminium"]
    material = AluminiumLondonModel(
        london_penetration_nm=float(aluminium["london_penetration_nm"]),
        clean_coherence_length_nm=float(aluminium["clean_coherence_length_nm"]),
        scattering_factor=float(aluminium["scattering_factor"]),
    )
    points = {
        str(int(float(thickness))): evaluate_point(
            float(thickness), material=material, resonator=resonator
        ).as_dict()
        for thickness in data["grid"]["reported_thicknesses_nm"]
    }
    return {
        "project": project,
        "layout": layout,
        "topology": topology,
        "width_um": width_um,
        "centerline_length_um": centerline_length_um,
        "branch_squares": branch_squares,
        "pec_frequency_ghz": pec_frequency_ghz,
        "thickness_points": points,
    }


def _solve_width(target_ghz: float, data: dict, base_config: Path) -> float:
    geometry = data["geometry"]
    reference_key = str(int(float(data["grid"]["reference_thickness_nm"])))
    low = float(geometry["solve_min_width_um"])
    high = float(geometry["solve_max_width_um"])

    def frequency(width: float) -> float:
        return float(_evaluate(width, data, base_config)["thickness_points"][reference_key]["frequency_ghz"])

    if not frequency(low) > target_ghz > frequency(high):
        raise RuntimeError(f"target {target_ghz} GHz is outside the width bracket")
    for _ in range(60):
        middle = 0.5 * (low + high)
        if frequency(middle) > target_ghz:
            low = middle
        else:
            high = middle
    raw_width = 0.5 * (low + high)
    quantum = float(geometry["width_quantization_um"])
    decimal_places = max(0, -Decimal(str(quantum)).as_tuple().exponent)
    return round(round(raw_width / quantum) * quantum, decimal_places)


def _write_pattern(output: Path, pattern: dict, data: dict) -> None:
    output.mkdir(parents=True, exist_ok=True)
    layout = pattern.pop("layout")
    project = pattern.pop("project")
    write_gds(
        layout,
        output / "layout.gds",
        cell_name=f"fig1_{pattern['pattern_id']}",
        include_ports=False,
    )
    write_svg(layout, output / "layout.svg", show_ports=False, show_etch_guides=True)
    write_stackup_xml(output / "stackup.xml", project.stackup)
    manifest = {
        "schema_id": "geerlings.frequency_bank.pattern.v1",
        **pattern,
        "material_reference": {
            "reference_thickness_nm": data["grid"]["reference_thickness_nm"],
            "reported_thicknesses_nm": data["grid"]["reported_thicknesses_nm"],
            "temperature_mk": data["aluminium"]["temperature_mk"],
        },
        "artifacts": ["layout.gds", "layout.svg", "stackup.xml", "manifest.json"],
        "frequency_semantics": "calibrated_PEC_surrogate_plus_reduced_order_kinetic_inductance",
        "not_independent_full_wave_solve": True,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _write_preview(path: Path, patterns: list[dict]) -> None:
    figure, axes = plt.subplots(2, 4, figsize=(13.2, 7.3), constrained_layout=True)
    figure.patch.set_facecolor("#dbe2ea")
    for axis, pattern in zip(axes.flat, patterns):
        layout = pattern["layout"]
        bounds = pattern["topology"]
        x_min = -0.5 * bounds["cutout_width_um"]
        axis.set_facecolor("#dbe2ea")
        axis.add_patch(
            Rectangle(
                (x_min, 0.0),
                bounds["cutout_width_um"],
                bounds["cutout_height_um"],
                facecolor="#101820",
                edgecolor="#5f7180",
                linewidth=0.8,
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
                    linewidth=0.1,
                )
            )
        margin = 12.0
        axis.set_xlim(x_min - margin, -x_min + margin)
        axis.set_ylim(-margin, bounds["cutout_height_um"] + margin)
        axis.set_aspect("equal")
        axis.set_axis_off()
        axis.set_title(
            f"{pattern['pattern_id']}  {pattern['target_frequency_ghz']:.3f} GHz\n"
            f"width {pattern['width_um']:.3f} um",
            fontsize=10,
            color="#193549",
        )
    figure.suptitle(
        "Eight Fig. 1-topology patterns | Al 200 nm nominal | 10 MHz spacing",
        fontsize=15,
        color="#193549",
    )
    figure.savefig(path, dpi=180, facecolor="#dbe2ea")
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=PROJECT_ROOT / "configs" / "frequency_bank_8x.toml",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = _load(args.config)
    base_config = PROJECT_ROOT / data["base_design_config"]
    output = args.output or PROJECT_ROOT / data["output_directory"]
    if not output.is_absolute():
        output = PROJECT_ROOT / output
    output.mkdir(parents=True, exist_ok=True)

    grid = data["grid"]
    targets = symmetric_targets(
        float(grid["center_frequency_ghz"]),
        float(grid["spacing_mhz"]),
        int(grid["count"]),
    )
    reference_key = str(int(float(grid["reference_thickness_nm"])))
    patterns: list[dict] = []
    for index, target in enumerate(targets, start=1):
        width_um = _solve_width(target, data, base_config)
        evaluated = _evaluate(width_um, data, base_config)
        reference_frequency = float(
            evaluated["thickness_points"][reference_key]["frequency_ghz"]
        )
        evaluated.update(
            pattern_id=f"P{index:02d}",
            target_frequency_ghz=target,
            reference_frequency_ghz=reference_frequency,
            target_error_mhz=(reference_frequency - target) * 1.0e3,
        )
        patterns.append(evaluated)

    for pattern in patterns:
        serializable = {
            key: value
            for key, value in pattern.items()
            if key not in {"project", "layout"}
        }
        _write_pattern(
            output / "patterns" / pattern["pattern_id"],
            {**serializable, "project": pattern["project"], "layout": pattern["layout"]},
            data,
        )

    rows = []
    for pattern in patterns:
        points = pattern["thickness_points"]
        rows.append(
            {
                "pattern_id": pattern["pattern_id"],
                "target_al200_ghz": pattern["target_frequency_ghz"],
                "nominal_width_um": pattern["width_um"],
                "predicted_pec_ghz": pattern["pec_frequency_ghz"],
                "predicted_al150_ghz": points["150"]["frequency_ghz"],
                "predicted_al200_ghz": points["200"]["frequency_ghz"],
                "target_error_al200_mhz": pattern["target_error_mhz"],
                "inductor_centerline_um": pattern["centerline_length_um"],
                "branch_squares": pattern["branch_squares"],
                "cutout_width_um": pattern["topology"]["cutout_width_um"],
                "cutout_height_um": pattern["topology"]["cutout_height_um"],
            }
        )
    with (output / "frequency_table.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    slope_mhz_per_um = (
        float(data["calibration"]["second_pec_frequency_ghz"])
        - float(data["calibration"]["anchor_pec_frequency_ghz"])
    ) / (
        float(data["calibration"]["second_width_um"])
        - float(data["calibration"]["anchor_width_um"])
    ) * 1.0e3
    summary = {
        "schema_id": "geerlings.frequency_bank.v1",
        "count": len(patterns),
        "center_frequency_ghz": grid["center_frequency_ghz"],
        "nominal_spacing_mhz": grid["spacing_mhz"],
        "reference_thickness_nm": grid["reference_thickness_nm"],
        "calibration": {
            **data["calibration"],
            "pec_slope_mhz_per_um": slope_mhz_per_um,
            "method": "linear interpolation of two 5 um local-metal PEC solves",
        },
        "patterns": rows,
        "uncertainty": {
            "absolute_pec_mesh_envelope_ghz": data["calibration"]["pec_mesh_uncertainty_ghz"],
            "spacing_is_nominal_not_resolved_by_absolute_mesh_envelope": True,
        },
    }
    (output / "results.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _write_preview(output / "frequency_bank_preview.png", patterns)
    print(json.dumps({"output": str(output), "patterns": rows}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
