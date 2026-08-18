"""Position-only external-Qc equalization for four-resonator chips.

The cavity inputs are compact samples from direct 3-D FEM calculations.  This
module interpolates those samples and applies a coupled-mode conversion for the
micron-scale resonators; it does not relabel that conversion as a conformal
full-wave solve.
"""

from __future__ import annotations

from copy import deepcopy
from itertools import permutations
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq, differential_evolution, least_squares


def load_json(path: Path) -> dict[str, Any]:
    """Load a UTF-8 JSON object."""

    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return payload


def resonator_external_qc(
    resonator_frequency_ghz: float,
    cavity_frequency_ghz: float,
    cavity_port_qc: float,
    coupling_ghz: float,
) -> float:
    """Convert cavity-port coupling to resonator-to-port external Qc."""

    if min(resonator_frequency_ghz, cavity_frequency_ghz, cavity_port_qc, coupling_ghz) <= 0:
        raise ValueError("frequencies, Qc, and coupling must be positive")
    detuning_ghz = cavity_frequency_ghz - resonator_frequency_ghz
    return (
        resonator_frequency_ghz
        * cavity_port_qc
        * detuning_ghz**2
        / (cavity_frequency_ghz * coupling_ghz**2)
    )


def frequency_map(bank: Mapping[str, Any]) -> dict[str, float]:
    """Extract the Al-200-nm frequency assigned to each bank pattern."""

    patterns = bank.get("patterns")
    if not isinstance(patterns, list) or not patterns:
        raise ValueError("frequency bank has no patterns")
    result = {
        str(row["pattern_id"]): float(row["predicted_al200_ghz"])
        for row in patterns
    }
    if len(result) != len(patterns):
        raise ValueError("frequency bank contains duplicate pattern IDs")
    return result


def merge_current_cavity_field_calibration(
    calibration: Mapping[str, Any],
    field_result: Mapping[str, Any],
) -> dict[str, Any]:
    """Combine current-depth FEM fields with historical direct pin sweeps.

    The current closed-cavity frequency replaces the weak-pin anchor frequency;
    the historical pin-induced frequency shifts and cavity-port Qc values are
    retained.  This makes the hybrid nature of the absolute-Qc estimate explicit
    while using the current cavity-height field shapes for position balancing.
    """

    merged = deepcopy(dict(calibration))
    raw_pin_anchors = merged.get("pin_anchors")
    profiles = field_result.get("profiles")
    cases = field_result.get("cases")
    if not isinstance(raw_pin_anchors, Mapping):
        raise ValueError("calibration has no pin_anchors object")
    pin_anchors = deepcopy(dict(raw_pin_anchors))
    raw_aliases = merged.get("pin_anchor_aliases", {})
    if not isinstance(raw_aliases, Mapping):
        raise ValueError("pin_anchor_aliases must be an object")
    pin_anchor_aliases = {
        str(cavity): str(source) for cavity, source in raw_aliases.items()
    }
    for cavity, source in pin_anchor_aliases.items():
        if cavity in pin_anchors:
            raise ValueError(f"pin anchor alias target {cavity} already exists")
        if source not in pin_anchors:
            raise ValueError(f"pin anchor alias source {source} does not exist")
        pin_anchors[cavity] = deepcopy(pin_anchors[source])
    merged["pin_anchors"] = pin_anchors
    if not isinstance(profiles, Mapping) or not isinstance(cases, Mapping):
        raise ValueError("current cavity FEM result has no profiles/cases objects")
    if set(pin_anchors) != set(profiles) or set(pin_anchors) != set(cases):
        raise ValueError("pin and current-field cavity cases do not match")

    profiles_by_y: dict[str, dict[str, list[dict[str, float]]]] = {}
    field_grids: dict[str, dict[str, dict[str, Any]]] = {}
    frequency_alignment: dict[str, dict[str, float]] = {}
    for cavity in pin_anchors:
        samples = cases[cavity].get("samples")
        if not isinstance(samples, list) or not samples:
            raise ValueError(f"current FEM case {cavity} has no field samples")
        by_y: dict[float, list[Mapping[str, Any]]] = {}
        for row in samples:
            by_y.setdefault(round(float(row["y_mm"]), 6), []).append(row)
        profiles_by_y[cavity] = {}
        for y_mm, y_rows in sorted(by_y.items()):
            by_x: dict[float, list[float]] = {}
            for row in y_rows:
                x_mm = round(abs(float(row["x_mm"])), 6)
                by_x.setdefault(x_mm, []).append(
                    float(row["e_in_plane_rms_normalized"])
                )
            profiles_by_y[cavity][f"{y_mm:.6g}"] = [
                {
                    "x_mm": x_mm,
                    "e_in_plane_rms_normalized_mean": float(np.mean(values)),
                    "contributing_samples": len(values),
                }
                for x_mm, values in sorted(by_x.items())
            ]

        field_grids[cavity] = {}
        for chip, sign in (("left", -1.0), ("right", 1.0)):
            chip_rows = [
                row for row in samples if math.copysign(1.0, float(row["x_mm"])) == sign
            ]
            local_x_values = sorted(
                {
                    round(
                        float(row["x_mm"])
                        - (-3.525 if chip == "left" else 3.525),
                        6,
                    )
                    for row in chip_rows
                }
            )
            y_values = sorted({round(float(row["y_mm"]), 6) for row in chip_rows})
            values = []
            for local_x in local_x_values:
                global_x = local_x + (-3.525 if chip == "left" else 3.525)
                value_row = []
                for y_mm in y_values:
                    matches = [
                        row
                        for row in chip_rows
                        if math.isclose(float(row["x_mm"]), global_x, abs_tol=1.0e-9)
                        and math.isclose(float(row["y_mm"]), y_mm, abs_tol=1.0e-9)
                    ]
                    if len(matches) != 1:
                        raise ValueError(
                            f"current FEM field grid is incomplete at {cavity}/{chip} "
                            f"x={local_x}, y={y_mm}"
                        )
                    value_row.append(float(matches[0]["e_in_plane_rms_normalized"]))
                values.append(value_row)
            field_grids[cavity][chip] = {
                "local_x_mm": local_x_values,
                "y_mm": y_values,
                "e_in_plane_rms_normalized": values,
            }

        rows = sorted(
            pin_anchors[cavity], key=lambda row: float(row["pin_length_mm"])
        )
        weak_anchor_frequency = float(rows[0]["frequency_ghz"])
        current_frequency = float(cases[cavity]["fundamental"]["frequency_ghz"])
        offset = current_frequency - weak_anchor_frequency
        for row in rows:
            row["historical_frequency_ghz"] = float(row["frequency_ghz"])
            row["frequency_ghz"] = float(row["frequency_ghz"]) + offset
        pin_anchors[cavity] = rows
        frequency_alignment[cavity] = {
            "reference_pin_length_mm": float(rows[0]["pin_length_mm"]),
            "historical_reference_frequency_ghz": weak_anchor_frequency,
            "current_closed_cavity_frequency_ghz": current_frequency,
            "applied_frequency_offset_ghz": offset,
        }

    merged["profiles"] = deepcopy(dict(profiles))
    merged["profiles_by_y"] = profiles_by_y
    merged["field_grids"] = field_grids
    provenance = dict(merged.get("provenance", {}))
    provenance["current_field_fem"] = {
        "schema_id": field_result.get("schema_id"),
        "solver": field_result.get("solver"),
        "validity": field_result.get("validity"),
        "frequency_alignment": frequency_alignment,
    }
    provenance["absolute_qc_semantics"] = (
        "historical direct pin-Qc anchors with current closed-cavity frequency "
        "alignment; not a current-geometry port full-wave sweep"
    )
    provenance["pin_anchor_aliases_applied"] = pin_anchor_aliases
    merged["provenance"] = provenance
    return merged


def _profile(rows: Sequence[Mapping[str, Any]]) -> tuple[PchipInterpolator, float, float]:
    xs = np.asarray([float(row["x_mm"]) for row in rows], dtype=float)
    fields = np.asarray(
        [float(row["e_in_plane_rms_normalized_mean"]) for row in rows],
        dtype=float,
    )
    if len(xs) < 3 or np.any(np.diff(xs) <= 0):
        raise ValueError("field profile x coordinates must be strictly increasing")
    if np.any(fields <= 0) or np.any(np.diff(fields) >= 0):
        raise ValueError("field profile must be positive and strictly decreasing")
    return PchipInterpolator(xs, fields), float(xs[0]), float(xs[-1])


def _pin_curves(
    rows: Sequence[Mapping[str, Any]],
) -> tuple[np.ndarray, PchipInterpolator, PchipInterpolator]:
    ordered = sorted(rows, key=lambda row: float(row["pin_length_mm"]))
    lengths = np.asarray([float(row["pin_length_mm"]) for row in ordered])
    frequencies = np.asarray([float(row["frequency_ghz"]) for row in ordered])
    qcs = np.asarray([float(row["cavity_port_qc"]) for row in ordered])
    if len(lengths) < 3 or np.any(np.diff(lengths) <= 0):
        raise ValueError("pin lengths must be strictly increasing")
    if np.any(frequencies <= 0) or np.any(qcs <= 0):
        raise ValueError("pin FEM anchors must be positive")
    return (
        lengths,
        PchipInterpolator(lengths, frequencies),
        PchipInterpolator(lengths, np.log(qcs)),
    )


def _inverse_field(
    profile: PchipInterpolator,
    target: float,
    x_min: float,
    x_max: float,
) -> float:
    low_field = float(profile(x_max))
    high_field = float(profile(x_min))
    if not low_field <= target <= high_field:
        raise ValueError(
            f"required field {target:.6g} is outside calibrated range "
            f"[{low_field:.6g}, {high_field:.6g}]"
        )
    return float(brentq(lambda x: float(profile(x)) - target, x_min, x_max))


def optimize_qc_positions(
    frequencies_ghz: Mapping[str, float],
    calibration: Mapping[str, Any],
    sites: Sequence[Mapping[str, Any]],
    *,
    target_qcs: Sequence[float] = (1.0e5, 1.0e6, 1.0e7),
    mean_coupling_mhz: float = 25.0,
    chip_center_abs_x_mm: float = 3.525,
    chip_size_mm: float = 5.05,
    footprint_size_mm: float = 0.626,
    chip_edge_clearance_mm: float = 0.1,
    pin_wall_thickness_mm: float = 7.5,
    fixed_point_iterations: int = 5,
    position_tolerances_mm: Sequence[float] = (0.002, 0.005, 0.010),
) -> dict[str, Any]:
    """Balance Qc using x position only while keeping every rotation at zero."""

    if mean_coupling_mhz <= 0:
        raise ValueError("mean coupling must be positive")
    if fixed_point_iterations < 1:
        raise ValueError("fixed_point_iterations must be positive")
    targets = sorted(float(value) for value in target_qcs)
    if not targets or any(value <= 0 for value in targets):
        raise ValueError("target Qc values must be positive")

    frequencies = {str(key): float(value) for key, value in frequencies_ghz.items()}
    site_patterns = [str(site["pattern_id"]) for site in sites]
    if sorted(site_patterns) != sorted(frequencies):
        raise ValueError("site assignment must contain every bank pattern exactly once")
    for site in sites:
        if site["chip"] not in {"left", "right"}:
            raise ValueError("site chip must be 'left' or 'right'")

    profile_payload = calibration.get("profiles")
    profiles_by_y_payload = calibration.get("profiles_by_y")
    pin_payload = calibration.get("pin_anchors")
    if not isinstance(profile_payload, Mapping) or not isinstance(pin_payload, Mapping):
        raise ValueError("calibration requires profiles and pin_anchors objects")
    cavity_keys = list(profile_payload)
    if len(cavity_keys) != 2 or set(cavity_keys) != set(pin_payload):
        raise ValueError("calibration must contain the same two cavity cases")

    site_y = {str(site["pattern_id"]): float(site["y_mm"]) for site in sites}
    profiles: dict[
        str, dict[str, tuple[PchipInterpolator, float, float]]
    ] = {}
    curves: dict[str, dict[str, Any]] = {}
    for cavity in cavity_keys:
        common_profile = _profile(profile_payload[cavity])
        profiles[cavity] = {}
        for pattern in frequencies:
            selected_profile = common_profile
            if isinstance(profiles_by_y_payload, Mapping):
                by_y = profiles_by_y_payload.get(cavity)
                if not isinstance(by_y, Mapping):
                    raise ValueError(f"profiles_by_y has no case {cavity}")
                matches = [
                    rows
                    for key, rows in by_y.items()
                    if math.isclose(float(key), site_y[pattern], abs_tol=1.0e-9)
                ]
                if len(matches) != 1:
                    raise ValueError(
                        f"profiles_by_y has no unique y={site_y[pattern]} for {cavity}"
                    )
                selected_profile = _profile(matches[0])
            profiles[cavity][pattern] = selected_profile
        lengths, fc_curve, log_qc_curve = _pin_curves(pin_payload[cavity])
        curves[cavity] = {
            "lengths": lengths,
            "frequency": fc_curve,
            "log_qc": log_qc_curve,
        }

    positions = {
        cavity: {pattern: chip_center_abs_x_mm for pattern in frequencies}
        for cavity in cavity_keys
    }

    def q_values(cavity: str, length: float) -> dict[str, float]:
        curve = curves[cavity]
        cavity_frequency = float(curve["frequency"](length))
        cavity_qc = math.exp(float(curve["log_qc"](length)))
        sampled = {
            pattern: float(
                profiles[cavity][pattern][0](positions[cavity][pattern])
            )
            for pattern in frequencies
        }
        field_reference = float(np.mean(list(sampled.values())))
        return {
            pattern: resonator_external_qc(
                frequency,
                cavity_frequency,
                cavity_qc,
                mean_coupling_mhz * sampled[pattern] / field_reference / 1000.0,
            )
            for pattern, frequency in frequencies.items()
        }

    def solve_pin_length(cavity: str, target_qc: float) -> float:
        lengths = curves[cavity]["lengths"]

        def residual(length: float) -> float:
            values = list(q_values(cavity, length).values())
            return float(np.mean(np.log(values)) - math.log(target_qc))

        low = float(lengths.min())
        high = float(lengths.max())
        if residual(low) * residual(high) > 0:
            raise ValueError(
                f"target Qc {target_qc:g} is outside pin calibration for {cavity}"
            )
        return float(brentq(residual, low, high))

    settings: dict[str, dict[float, float]] = {}
    for _ in range(fixed_point_iterations):
        settings = {
            cavity: {target: solve_pin_length(cavity, target) for target in targets}
            for cavity in cavity_keys
        }
        for cavity in cavity_keys:
            relative_logs = {pattern: [] for pattern in frequencies}
            for target in targets:
                cavity_frequency = float(
                    curves[cavity]["frequency"](settings[cavity][target])
                )
                raw = {
                    pattern: math.sqrt(frequency / cavity_frequency)
                    * abs(cavity_frequency - frequency)
                    for pattern, frequency in frequencies.items()
                }
                geometric_mean = math.exp(float(np.mean(np.log(list(raw.values())))))
                for pattern in frequencies:
                    relative_logs[pattern].append(
                        math.log(raw[pattern] / geometric_mean)
                    )
            center_field = float(
                np.mean(
                    [
                        float(profiles[cavity][pattern][0](chip_center_abs_x_mm))
                        for pattern in frequencies
                    ]
                )
            )
            for pattern in frequencies:
                relative = math.exp(float(np.mean(relative_logs[pattern])))
                profile, x_min, x_max = profiles[cavity][pattern]
                positions[cavity][pattern] = _inverse_field(
                    profile,
                    center_field * relative,
                    x_min,
                    x_max,
                )

    settings = {
        cavity: {target: solve_pin_length(cavity, target) for target in targets}
        for cavity in cavity_keys
    }

    placements: list[dict[str, Any]] = []
    for cavity in cavity_keys:
        for site in sites:
            pattern = str(site["pattern_id"])
            chip = str(site["chip"])
            x_abs = positions[cavity][pattern]
            x_global = -x_abs if chip == "left" else x_abs
            chip_center = -chip_center_abs_x_mm if chip == "left" else chip_center_abs_x_mm
            placements.append(
                {
                    "cavity": cavity,
                    "pattern_id": pattern,
                    "chip": chip,
                    "column": str(site.get("column", "")),
                    "x_mm": x_global,
                    "chip_local_x_mm": x_global - chip_center,
                    "y_mm": float(site["y_mm"]),
                    "rotation_deg": 0.0,
                    "frequency_ghz": frequencies[pattern],
                    "sampled_field": float(profiles[cavity][pattern][0](x_abs)),
                }
            )

    usable_half = 0.5 * chip_size_mm - chip_edge_clearance_mm
    footprint_half = 0.5 * footprint_size_mm
    minimum_gap = math.inf
    for cavity in cavity_keys:
        for chip in ("left", "right"):
            selected = [
                row
                for row in placements
                if row["cavity"] == cavity and row["chip"] == chip
            ]
            for row in selected:
                if abs(row["chip_local_x_mm"]) + footprint_half > usable_half:
                    raise ValueError(f"{row['pattern_id']} exceeds chip x clearance")
                if abs(row["y_mm"]) + footprint_half > usable_half:
                    raise ValueError(f"{row['pattern_id']} exceeds chip y clearance")
            for index, first in enumerate(selected):
                for second in selected[index + 1 :]:
                    dx = abs(first["chip_local_x_mm"] - second["chip_local_x_mm"])
                    dy = abs(first["y_mm"] - second["y_mm"])
                    if dx < footprint_size_mm and dy < footprint_size_mm:
                        raise ValueError(
                            f"nominal footprint overlap: {first['pattern_id']} / "
                            f"{second['pattern_id']}"
                        )
                    minimum_gap = min(
                        minimum_gap,
                        max(dx - footprint_size_mm, dy - footprint_size_mm),
                    )

    position_tolerance_ratios: dict[str, float] = {}
    for tolerance_mm in sorted(float(value) for value in position_tolerances_mm):
        if tolerance_mm <= 0:
            raise ValueError("position tolerances must be positive")
        multipliers = []
        for row in placements:
            cavity = str(row["cavity"])
            pattern = str(row["pattern_id"])
            x_nominal = abs(float(row["x_mm"]))
            profile, x_min, x_max = profiles[cavity][pattern]
            nominal_field = float(profile(x_nominal))
            for shifted_x in (x_nominal - tolerance_mm, x_nominal + tolerance_mm):
                if not x_min <= shifted_x <= x_max:
                    raise ValueError("position tolerance leaves calibrated field range")
                shifted_field = float(profile(shifted_x))
                multipliers.append((nominal_field / shifted_field) ** 2)
        position_tolerance_ratios[f"{1000.0 * tolerance_mm:g}"] = (
            max(multipliers) / min(multipliers)
        )

    pin_calibration: list[dict[str, Any]] = []
    qc_predictions: list[dict[str, Any]] = []
    global_ratios: dict[str, float] = {}
    for target in targets:
        all_values: list[float] = []
        target_rows: list[dict[str, Any]] = []
        for cavity in cavity_keys:
            length = settings[cavity][target]
            curve = curves[cavity]
            cavity_frequency = float(curve["frequency"](length))
            cavity_qc = math.exp(float(curve["log_qc"](length)))
            values = q_values(cavity, length)
            all_values.extend(values.values())
            target_rows.append(
                {
                    "target_resonator_external_qc": target,
                    "cavity": cavity,
                    "pin_length_mm": length,
                    "pin_protrusion_into_cavity_mm": length - pin_wall_thickness_mm,
                    "interpolated_cavity_frequency_ghz": cavity_frequency,
                    "interpolated_cavity_port_qc": cavity_qc,
                    "resonator_qc_min": min(values.values()),
                    "resonator_qc_max": max(values.values()),
                    "resonator_qc_max_to_min": max(values.values()) / min(values.values()),
                }
            )
            for pattern, value in values.items():
                qc_predictions.append(
                    {
                        "target_resonator_external_qc": target,
                        "cavity": cavity,
                        "pattern_id": pattern,
                        "frequency_ghz": frequencies[pattern],
                        "predicted_external_qc": value,
                    }
                )
        global_ratio = max(all_values) / min(all_values)
        global_ratios[f"{target:.0e}"] = global_ratio
        for row in target_rows:
            row["global_resonator_qc_max_to_min"] = global_ratio
            pin_calibration.append(row)

    return {
        "schema_id": "geerlings.four_resonator_chip.position_only_qc.v1",
        "qc_definition": (
            "external coupling Q from each chip resonator to its SMA port "
            "through the detuned cavity mode"
        ),
        "constraints": {
            "resonator_rotation_deg": 0.0,
            "resonators_per_chip": 4,
            "chips_per_cavity": 2,
            "chip_size_mm": [chip_size_mm, chip_size_mm],
            "position_may_change": True,
        },
        "coupled_mode_model": "Qc_res = fr*Qc_cavity*(fc-fr)^2/(fc*g^2)",
        "mean_resonator_cavity_g_mhz_assumption": mean_coupling_mhz,
        "frequency_bank_ghz": frequencies,
        "target_resonator_external_qcs": targets,
        "pin_calibration": pin_calibration,
        "global_resonator_qc_max_to_min": global_ratios,
        "minimum_nominal_footprint_edge_gap_mm": minimum_gap,
        "position_tolerance_worst_case_qc_max_to_min": {
            "tolerance_unit": "um",
            "independent_x_error_envelope": position_tolerance_ratios,
            "semantics": (
                "deterministic worst case from local FEM field gradients; excludes "
                "linewidth, film, cavity-machining, and model-form uncertainty"
            ),
        },
        "placements": placements,
        "qc_predictions": qc_predictions,
        "fem_calibration_provenance": calibration.get("provenance", {}),
        "validity": {
            "cavity_pin_response": "direct 3-D FEM anchors with PCHIP interpolation",
            "field_position_response": (
                "direct 3-D FEM field samples reduced to an x profile with PCHIP interpolation"
            ),
            "resonator_external_qc": (
                "coupled-mode conversion, not a conformal micron-to-centimeter full-wave solve; "
                "absolute values scale as 1/g^2"
            ),
            "relative_position_equalization": (
                "more robust than the absolute Qc calibration but still inherits cavity-model "
                "and fabrication uncertainty"
            ),
        },
    }


def optimize_rectangle_qc_positions(
    frequencies_ghz: Mapping[str, float],
    calibration: Mapping[str, Any],
    sites: Sequence[Mapping[str, Any]],
    *,
    target_qcs: Sequence[float] = (1.0e6,),
    mean_coupling_mhz: float = 25.0,
    chip_center_abs_x_mm: float = 3.525,
    chip_size_mm: float = 5.05,
    chip_edge_clearance_mm: float = 0.10,
    frequency_split_ghz: float = 9.0,
    coordinate_limit_mm: float = 1.18,
    minimum_center_separation_mm: float = 0.60,
    maximum_ground_cutout_width_mm: float = 0.425968,
    maximum_ground_cutout_height_mm: float = 0.395,
    pin_wall_thickness_mm: float = 7.5,
    fixed_point_iterations: int = 2,
    optimizer_seed: int = 17,
    optimizer_maxiter: int = 60,
    optimizer_popsize: int = 8,
    centered_square_side_mm: float | None = None,
    position_tolerances_mm: Sequence[float] = (0.002, 0.005, 0.010, 0.020),
) -> dict[str, Any]:
    """Equalize Qc while constraining each chip to four rectangle corners.

    The rectangle center, width, and height are continuous variables.  For each
    trial geometry, all 4! assignments of the four chip frequencies to its four
    corners are evaluated.  The two chips are optimized together for each
    cavity case, so the reported spread covers all eight resonators in that
    cavity.
    """

    frequencies = {str(key): float(value) for key, value in frequencies_ghz.items()}
    targets = sorted({float(value) for value in target_qcs})
    if len(frequencies) != 8 or any(value <= 0 for value in frequencies.values()):
        raise ValueError("rectangle optimization requires eight positive frequencies")
    if not targets or any(value <= 0 for value in targets):
        raise ValueError("target Qc values must be positive")
    if mean_coupling_mhz <= 0:
        raise ValueError("mean coupling must be positive")
    if frequency_split_ghz <= 0:
        raise ValueError("frequency split must be positive")
    if fixed_point_iterations < 1:
        raise ValueError("fixed_point_iterations must be at least one")
    half_minimum = 0.5 * minimum_center_separation_mm
    if not 0 < half_minimum <= coordinate_limit_mm:
        raise ValueError("invalid rectangle separation or coordinate limit")
    if centered_square_side_mm is not None:
        centered_square_side_mm = float(centered_square_side_mm)
        if not minimum_center_separation_mm <= centered_square_side_mm <= 2.0 * coordinate_limit_mm:
            raise ValueError(
                "centered square side must satisfy separation and coordinate limits"
            )

    patterns_by_chip: dict[str, list[str]] = {"left": [], "right": []}
    site_by_pattern: dict[str, Mapping[str, Any]] = {}
    for site in sites:
        pattern = str(site["pattern_id"])
        chip = str(site["chip"])
        if chip not in patterns_by_chip:
            raise ValueError(f"unknown chip {chip!r}")
        if pattern in site_by_pattern:
            raise ValueError(f"duplicate site for {pattern}")
        site_by_pattern[pattern] = site
        patterns_by_chip[chip].append(pattern)
    if set(site_by_pattern) != set(frequencies):
        raise ValueError("site patterns and frequency-bank patterns do not match")
    for chip, patterns in patterns_by_chip.items():
        if len(patterns) != 4:
            raise ValueError(f"{chip} chip must contain four patterns")
        patterns.sort(key=lambda pattern: (frequencies[pattern], pattern))
    frequency_partition = {
        chip: {
            "below_9ghz": [
                pattern
                for pattern in patterns_by_chip[chip]
                if frequencies[pattern] < frequency_split_ghz
            ],
            "above_9ghz": [
                pattern
                for pattern in patterns_by_chip[chip]
                if frequencies[pattern] > frequency_split_ghz
            ],
        }
        for chip in ("left", "right")
    }
    if any(
        len(frequency_partition[chip]["below_9ghz"])
        + len(frequency_partition[chip]["above_9ghz"])
        != 4
        for chip in ("left", "right")
    ):
        raise ValueError("a resonator lies exactly on the frequency split")

    pin_anchors = calibration.get("pin_anchors")
    field_grids = calibration.get("field_grids")
    if not isinstance(pin_anchors, Mapping) or not isinstance(field_grids, Mapping):
        raise ValueError("calibration must contain pin_anchors and field_grids")
    if set(pin_anchors) != set(field_grids):
        raise ValueError("pin and field cavity cases do not match")
    cavity_keys = list(pin_anchors)

    curves: dict[str, dict[str, Any]] = {}
    interpolators: dict[str, dict[str, RegularGridInterpolator]] = {}
    grid_bounds: dict[str, dict[str, tuple[float, float, float, float]]] = {}
    maximum_tolerance = max((float(value) for value in position_tolerances_mm), default=0.0)
    for cavity in cavity_keys:
        lengths, frequency_curve, log_qc_curve = _pin_curves(pin_anchors[cavity])
        curves[cavity] = {
            "lengths": lengths,
            "frequency": frequency_curve,
            "log_qc": log_qc_curve,
        }
        interpolators[cavity] = {}
        grid_bounds[cavity] = {}
        for chip in ("left", "right"):
            grid = field_grids[cavity].get(chip)
            if not isinstance(grid, Mapping):
                raise ValueError(f"missing field grid for {cavity}/{chip}")
            xs = np.asarray(grid["local_x_mm"], dtype=float)
            ys = np.asarray(grid["y_mm"], dtype=float)
            values = np.asarray(grid["e_in_plane_rms_normalized"], dtype=float)
            if (
                len(xs) < 2
                or len(ys) < 2
                or np.any(np.diff(xs) <= 0)
                or np.any(np.diff(ys) <= 0)
                or values.shape != (len(xs), len(ys))
                or np.any(values <= 0)
            ):
                raise ValueError(f"invalid field grid for {cavity}/{chip}")
            required = coordinate_limit_mm + maximum_tolerance
            if xs[0] > -required or xs[-1] < required or ys[0] > -required or ys[-1] < required:
                raise ValueError(
                    f"field grid for {cavity}/{chip} does not cover rectangle plus tolerance"
                )
            interpolators[cavity][chip] = RegularGridInterpolator(
                (xs, ys), values, method="linear", bounds_error=True
            )
            grid_bounds[cavity][chip] = (
                float(xs[0]),
                float(xs[-1]),
                float(ys[0]),
                float(ys[-1]),
            )

    slot_names = {
        "left": (
            ("outer", "lower"),
            ("outer", "upper"),
            ("inner", "lower"),
            ("inner", "upper"),
        ),
        "right": (
            ("inner", "lower"),
            ("inner", "upper"),
            ("outer", "lower"),
            ("outer", "upper"),
        ),
    }
    all_permutations = np.asarray(list(permutations(range(4))), dtype=int)
    permutation_indices_by_chip: dict[str, np.ndarray] = {}
    for chip in ("left", "right"):
        configured_columns = [
            str(site_by_pattern[pattern].get("column", ""))
            for pattern in patterns_by_chip[chip]
        ]
        if all(column in {"inner", "outer"} for column in configured_columns):
            slot_columns = [column for column, _ in slot_names[chip]]
            allowed = [
                permutation
                for permutation in all_permutations
                if all(
                    configured_columns[int(pattern_index)] == slot_columns[slot_index]
                    for slot_index, pattern_index in enumerate(permutation)
                )
            ]
            if not allowed:
                raise ValueError(f"column constraints have no assignment for {chip}")
            permutation_indices_by_chip[chip] = np.asarray(allowed, dtype=int)
        elif any(configured_columns):
            raise ValueError(f"column constraints are incomplete for {chip}")
        else:
            permutation_indices_by_chip[chip] = all_permutations

    def decode_rectangle(vector: Sequence[float]) -> dict[str, dict[str, float]]:
        rectangles: dict[str, dict[str, float]] = {}
        for chip_index, chip in enumerate(("left", "right")):
            half_x, fraction_x, half_y, fraction_y = (
                float(value) for value in vector[4 * chip_index : 4 * chip_index + 4]
            )
            center_x = fraction_x * (coordinate_limit_mm - half_x)
            center_y = fraction_y * (coordinate_limit_mm - half_y)
            rectangles[chip] = {
                "center_x_mm": center_x,
                "center_y_mm": center_y,
                "half_width_mm": half_x,
                "half_height_mm": half_y,
                "x_low_mm": center_x - half_x,
                "x_high_mm": center_x + half_x,
                "y_low_mm": center_y - half_y,
                "y_high_mm": center_y + half_y,
                "width_mm": 2.0 * half_x,
                "height_mm": 2.0 * half_y,
            }
        return rectangles

    def rectangle_slots(
        rectangles: Mapping[str, Mapping[str, float]],
    ) -> dict[str, list[tuple[float, float]]]:
        result: dict[str, list[tuple[float, float]]] = {}
        for chip in ("left", "right"):
            rectangle = rectangles[chip]
            x_by_column = (
                {
                    "outer": float(rectangle["x_low_mm"]),
                    "inner": float(rectangle["x_high_mm"]),
                }
                if chip == "left"
                else {
                    "inner": float(rectangle["x_low_mm"]),
                    "outer": float(rectangle["x_high_mm"]),
                }
            )
            y_by_row = {
                "lower": float(rectangle["y_low_mm"]),
                "upper": float(rectangle["y_high_mm"]),
            }
            result[chip] = [
                (x_by_column[column], y_by_row[row])
                for column, row in slot_names[chip]
            ]
        return result

    def sample_fields(
        cavity: str, slots: Mapping[str, Sequence[tuple[float, float]]]
    ) -> dict[str, np.ndarray]:
        sampled: dict[str, np.ndarray] = {}
        for chip in ("left", "right"):
            points = np.asarray(slots[chip], dtype=float)
            values = np.asarray(interpolators[cavity][chip](points), dtype=float)
            if np.any(values <= 0):
                raise ValueError(f"non-positive interpolated field in {cavity}/{chip}")
            sampled[chip] = values
        return sampled

    def best_assignment(
        cavity: str,
        pin_length_mm: float,
        vector: Sequence[float],
    ) -> tuple[float, dict[str, int], dict[str, Any]]:
        rectangles = decode_rectangle(vector)
        slots = rectangle_slots(rectangles)
        fields = sample_fields(cavity, slots)
        cavity_frequency = float(curves[cavity]["frequency"](pin_length_mm))
        log_q_by_chip: dict[str, np.ndarray] = {}
        for chip in ("left", "right"):
            permutation_indices = permutation_indices_by_chip[chip]
            pattern_base = np.asarray(
                [
                    math.log(frequencies[pattern])
                    + 2.0 * math.log(abs(cavity_frequency - frequencies[pattern]))
                    for pattern in patterns_by_chip[chip]
                ],
                dtype=float,
            )
            log_q_by_chip[chip] = (
                pattern_base[permutation_indices]
                - 2.0 * np.log(fields[chip])[None, :]
            )
        left = log_q_by_chip["left"]
        right = log_q_by_chip["right"]
        minima = np.minimum(left.min(axis=1)[:, None], right.min(axis=1)[None, :])
        maxima = np.maximum(left.max(axis=1)[:, None], right.max(axis=1)[None, :])
        sums = left.sum(axis=1)[:, None] + right.sum(axis=1)[None, :]
        sums_of_squares = (
            np.square(left).sum(axis=1)[:, None]
            + np.square(right).sum(axis=1)[None, :]
        )
        variance = np.maximum(sums_of_squares / 8.0 - np.square(sums / 8.0), 0.0)
        scores = maxima - minima + 0.05 * np.sqrt(variance)
        left_index, right_index = np.unravel_index(int(np.argmin(scores)), scores.shape)
        return (
            float(scores[left_index, right_index]),
            {"left": int(left_index), "right": int(right_index)},
            {"rectangles": rectangles, "slots": slots, "fields": fields},
        )

    def assignment_map(
        assignment_indices: Mapping[str, int],
    ) -> dict[str, tuple[str, int]]:
        result: dict[str, tuple[str, int]] = {}
        for chip in ("left", "right"):
            permutation = permutation_indices_by_chip[chip][assignment_indices[chip]]
            for slot_index, pattern_index in enumerate(permutation):
                pattern = patterns_by_chip[chip][int(pattern_index)]
                result[pattern] = (chip, slot_index)
        return result

    def q_values(
        cavity: str,
        pin_length_mm: float,
        vector: Sequence[float],
        assignment_indices: Mapping[str, int],
    ) -> tuple[dict[str, float], dict[str, Any]]:
        rectangles = decode_rectangle(vector)
        slots = rectangle_slots(rectangles)
        fields = sample_fields(cavity, slots)
        assignment = assignment_map(assignment_indices)
        field_reference = float(np.mean(np.concatenate(list(fields.values()))))
        cavity_frequency = float(curves[cavity]["frequency"](pin_length_mm))
        cavity_qc = math.exp(float(curves[cavity]["log_qc"](pin_length_mm)))
        values: dict[str, float] = {}
        for pattern, (chip, slot_index) in assignment.items():
            coupling_ghz = (
                mean_coupling_mhz
                * float(fields[chip][slot_index])
                / field_reference
                / 1000.0
            )
            values[pattern] = resonator_external_qc(
                frequencies[pattern], cavity_frequency, cavity_qc, coupling_ghz
            )
        return values, {
            "rectangles": rectangles,
            "slots": slots,
            "fields": fields,
            "field_reference": field_reference,
            "cavity_frequency_ghz": cavity_frequency,
            "cavity_port_qc": cavity_qc,
        }

    def solve_pin_length(
        cavity: str,
        target_qc: float,
        vector: Sequence[float],
        assignment_indices: Mapping[str, int],
    ) -> float:
        lengths = curves[cavity]["lengths"]

        def residual(length: float) -> float:
            values, _ = q_values(cavity, length, vector, assignment_indices)
            return float(np.mean(np.log(list(values.values()))) - math.log(target_qc))

        low = float(lengths.min())
        high = float(lengths.max())
        if residual(low) * residual(high) > 0:
            raise ValueError(f"target Qc {target_qc:g} is outside pin calibration for {cavity}")
        return float(brentq(residual, low, high))

    bounds = [
        (half_minimum, coordinate_limit_mm),
        (-1.0, 1.0),
        (half_minimum, coordinate_limit_mm),
        (-1.0, 1.0),
    ] * 2
    nominal_target = targets[0]
    optimized: dict[str, dict[str, Any]] = {}
    for cavity in cavity_keys:
        pin_length = float(np.median(curves[cavity]["lengths"]))
        previous_vector: np.ndarray | None = None
        assignment_indices: dict[str, int] = {}
        details: dict[str, Any] = {}
        score = math.inf
        if centered_square_side_mm is not None:
            half_side = 0.5 * centered_square_side_mm
            previous_vector = np.asarray(
                [half_side, 0.0, half_side, 0.0] * 2, dtype=float
            )
            for _ in range(fixed_point_iterations):
                score, assignment_indices, details = best_assignment(
                    cavity, pin_length, previous_vector
                )
                pin_length = solve_pin_length(
                    cavity, nominal_target, previous_vector, assignment_indices
                )
            optimizer_metadata = {
                "mode": "fixed_centered_square_with_corner_assignment_enumeration",
                "objective": score,
                "success": True,
                "message": "fixed geometry evaluated",
                "iterations": fixed_point_iterations,
                "function_evaluations": fixed_point_iterations,
                "seed": None,
                "local_refinement_success": None,
                "local_refinement_function_evaluations": 0,
            }
        else:
            result = None
            refinement = None
            for iteration in range(fixed_point_iterations):
                kwargs: dict[str, Any] = {}
                if previous_vector is not None:
                    kwargs["x0"] = previous_vector
                result = differential_evolution(
                    lambda vector: best_assignment(cavity, pin_length, vector)[0],
                    bounds,
                    seed=optimizer_seed,
                    maxiter=optimizer_maxiter,
                    popsize=optimizer_popsize,
                    polish=True,
                    tol=1.0e-9,
                    atol=1.0e-10,
                    workers=1,
                    updating="immediate",
                    **kwargs,
                )
                previous_vector = np.asarray(result.x, dtype=float)
                score, assignment_indices, details = best_assignment(
                    cavity, pin_length, previous_vector
                )
                for _ in range(3):
                    fixed_assignment = dict(assignment_indices)

                    def equal_q_residual(vector: Sequence[float]) -> np.ndarray:
                        values, _ = q_values(
                            cavity, pin_length, vector, fixed_assignment
                        )
                        logs = np.log(
                            [values[pattern] for pattern in sorted(values)]
                        )
                        return logs - float(np.mean(logs))

                    refinement = least_squares(
                        equal_q_residual,
                        previous_vector,
                        bounds=(
                            np.asarray([low for low, _ in bounds], dtype=float),
                            np.asarray([high for _, high in bounds], dtype=float),
                        ),
                        max_nfev=2000,
                        ftol=1.0e-13,
                        xtol=1.0e-13,
                        gtol=1.0e-13,
                    )
                    refined_vector = np.asarray(refinement.x, dtype=float)
                    refined_score, refined_assignment, refined_details = best_assignment(
                        cavity, pin_length, refined_vector
                    )
                    if refined_score > score + 1.0e-12:
                        break
                    previous_vector = refined_vector
                    score = refined_score
                    details = refined_details
                    if refined_assignment == assignment_indices:
                        assignment_indices = refined_assignment
                        break
                    assignment_indices = refined_assignment
                pin_length = solve_pin_length(
                    cavity, nominal_target, previous_vector, assignment_indices
                )
            if result is None or previous_vector is None:
                raise RuntimeError("rectangle optimizer did not run")
            optimizer_metadata = {
                "mode": "differential_evolution_plus_local_least_squares",
                "objective": score,
                "success": bool(result.success),
                "message": str(result.message),
                "iterations": int(result.nit),
                "function_evaluations": int(result.nfev),
                "seed": optimizer_seed,
                "local_refinement_success": bool(
                    refinement is not None and refinement.success
                ),
                "local_refinement_function_evaluations": (
                    int(refinement.nfev) if refinement is not None else 0
                ),
            }
        if previous_vector is None:
            raise RuntimeError("placement evaluation did not run")
        nominal_values, final_details = q_values(
            cavity, pin_length, previous_vector, assignment_indices
        )
        optimized[cavity] = {
            "vector": previous_vector,
            "assignment_indices": assignment_indices,
            "details": final_details,
            "nominal_pin_length_mm": pin_length,
            "nominal_q_values": nominal_values,
            "optimizer": optimizer_metadata,
        }

    placements: list[dict[str, Any]] = []
    rectangle_geometry: dict[str, dict[str, dict[str, float]]] = {}
    for cavity in cavity_keys:
        solution = optimized[cavity]
        assignment = assignment_map(solution["assignment_indices"])
        details = solution["details"]
        rectangle_geometry[cavity] = deepcopy(details["rectangles"])
        for pattern in sorted(frequencies, key=lambda value: (frequencies[value], value)):
            chip, slot_index = assignment[pattern]
            local_x, y_mm = details["slots"][chip][slot_index]
            column, row_name = slot_names[chip][slot_index]
            chip_center = -chip_center_abs_x_mm if chip == "left" else chip_center_abs_x_mm
            placements.append(
                {
                    "cavity": cavity,
                    "pattern_id": pattern,
                    "frequency_band": (
                        "below_9ghz"
                        if frequencies[pattern] < frequency_split_ghz
                        else "above_9ghz"
                    ),
                    "chip": chip,
                    "column": column,
                    "row": row_name,
                    "x_mm": chip_center + local_x,
                    "chip_local_x_mm": local_x,
                    "y_mm": y_mm,
                    "rotation_deg": 0.0,
                    "frequency_ghz": frequencies[pattern],
                    "sampled_field": float(details["fields"][chip][slot_index]),
                    "predicted_external_qc": float(solution["nominal_q_values"][pattern]),
                }
            )

    minimum_center_separation = math.inf
    minimum_cutout_gap = math.inf
    usable_half = 0.5 * chip_size_mm - chip_edge_clearance_mm
    for cavity in cavity_keys:
        for chip in ("left", "right"):
            selected = [
                row for row in placements if row["cavity"] == cavity and row["chip"] == chip
            ]
            for row in selected:
                if abs(float(row["chip_local_x_mm"])) + 0.5 * maximum_ground_cutout_width_mm > usable_half:
                    raise ValueError(f"{row['pattern_id']} exceeds chip x clearance")
                if abs(float(row["y_mm"])) + 0.5 * maximum_ground_cutout_height_mm > usable_half:
                    raise ValueError(f"{row['pattern_id']} exceeds chip y clearance")
            for index, first in enumerate(selected):
                for second in selected[index + 1 :]:
                    dx = abs(float(first["chip_local_x_mm"]) - float(second["chip_local_x_mm"]))
                    dy = abs(float(first["y_mm"]) - float(second["y_mm"]))
                    minimum_center_separation = min(
                        minimum_center_separation, math.hypot(dx, dy)
                    )
                    gap_x = dx - maximum_ground_cutout_width_mm
                    gap_y = dy - maximum_ground_cutout_height_mm
                    if gap_x < 0 and gap_y < 0:
                        raise ValueError(
                            f"ground cutout overlap: {first['pattern_id']} / {second['pattern_id']}"
                        )
                    minimum_cutout_gap = min(minimum_cutout_gap, max(gap_x, gap_y))

    pin_calibration: list[dict[str, Any]] = []
    qc_predictions: list[dict[str, Any]] = []
    global_ratios: dict[str, float] = {}
    target_q_by_cavity: dict[float, dict[str, dict[str, float]]] = {}
    for target in targets:
        target_q_by_cavity[target] = {}
        all_values: list[float] = []
        target_rows: list[dict[str, Any]] = []
        for cavity in cavity_keys:
            solution = optimized[cavity]
            length = solve_pin_length(
                cavity,
                target,
                solution["vector"],
                solution["assignment_indices"],
            )
            values, details = q_values(
                cavity, length, solution["vector"], solution["assignment_indices"]
            )
            target_q_by_cavity[target][cavity] = values
            all_values.extend(values.values())
            target_rows.append(
                {
                    "target_resonator_external_qc": target,
                    "cavity": cavity,
                    "pin_length_mm": length,
                    "pin_protrusion_into_cavity_mm": length - pin_wall_thickness_mm,
                    "interpolated_cavity_frequency_ghz": details["cavity_frequency_ghz"],
                    "interpolated_cavity_port_qc": details["cavity_port_qc"],
                    "resonator_qc_min": min(values.values()),
                    "resonator_qc_max": max(values.values()),
                    "resonator_qc_max_to_min": max(values.values()) / min(values.values()),
                }
            )
            placement_by_pattern = {
                row["pattern_id"]: row
                for row in placements
                if row["cavity"] == cavity
            }
            for pattern in sorted(values, key=lambda value: (frequencies[value], value)):
                placement = placement_by_pattern[pattern]
                qc_predictions.append(
                    {
                        "target_resonator_external_qc": target,
                        "cavity": cavity,
                        "chip": placement["chip"],
                        "frequency_band": placement["frequency_band"],
                        "pattern_id": pattern,
                        "frequency_ghz": frequencies[pattern],
                        "predicted_external_qc": values[pattern],
                    }
                )
        global_ratio = max(all_values) / min(all_values)
        global_ratios[f"{target:.0e}"] = global_ratio
        for row in target_rows:
            row["global_resonator_qc_max_to_min"] = global_ratio
            pin_calibration.append(row)

    tolerance_ratios: dict[str, float] = {}
    nominal_q = target_q_by_cavity[nominal_target]
    for tolerance_mm in sorted(float(value) for value in position_tolerances_mm):
        if tolerance_mm <= 0:
            raise ValueError("position tolerances must be positive")
        candidates: list[float] = []
        for row in placements:
            cavity = str(row["cavity"])
            chip = str(row["chip"])
            pattern = str(row["pattern_id"])
            nominal_field = float(row["sampled_field"])
            x_nominal = float(row["chip_local_x_mm"])
            y_nominal = float(row["y_mm"])
            x_min, x_max, y_min, y_max = grid_bounds[cavity][chip]
            for dx in (-tolerance_mm, 0.0, tolerance_mm):
                for dy in (-tolerance_mm, 0.0, tolerance_mm):
                    shifted_x = x_nominal + dx
                    shifted_y = y_nominal + dy
                    if not (x_min <= shifted_x <= x_max and y_min <= shifted_y <= y_max):
                        raise ValueError("position tolerance leaves calibrated field grid")
                    shifted_field = float(
                        np.asarray(
                            interpolators[cavity][chip]((shifted_x, shifted_y))
                        ).item()
                    )
                    candidates.append(
                        nominal_q[cavity][pattern]
                        * (nominal_field / shifted_field) ** 2
                    )
        tolerance_ratios[f"{1000.0 * tolerance_mm:g}"] = max(candidates) / min(candidates)

    return {
        "schema_id": (
            "geerlings.four_resonator_chip.centered_square_qc.v1"
            if centered_square_side_mm is not None
            else "geerlings.four_resonator_chip.rectangle_qc.v1"
        ),
        "qc_definition": (
            "external coupling Q from each chip resonator to its SMA port "
            "through the detuned cavity mode"
        ),
        "constraints": {
            "centers_form_axis_aligned_rectangle": True,
            "resonator_rotation_deg": 0.0,
            "resonators_per_chip": 4,
            "chips_per_cavity": 2,
            "chip_size_mm": [chip_size_mm, chip_size_mm],
            "rectangle_coordinate_limit_mm": coordinate_limit_mm,
            "minimum_center_separation_mm": minimum_center_separation_mm,
            "centered_on_chip": centered_square_side_mm is not None,
            "square_side_mm": centered_square_side_mm,
            "frequency_split_ghz": frequency_split_ghz,
            "frequency_partition": frequency_partition,
        },
        "coupled_mode_model": "Qc_res = fr*Qc_cavity*(fc-fr)^2/(fc*g^2)",
        "mean_resonator_cavity_g_mhz_assumption": mean_coupling_mhz,
        "frequency_bank_ghz": frequencies,
        "target_resonator_external_qcs": targets,
        "rectangle_geometry": rectangle_geometry,
        "optimizer": {cavity: optimized[cavity]["optimizer"] for cavity in cavity_keys},
        "pin_calibration": pin_calibration,
        "global_resonator_qc_max_to_min": global_ratios,
        "minimum_resonator_center_separation_mm": minimum_center_separation,
        "minimum_ground_cutout_edge_gap_mm": minimum_cutout_gap,
        "minimum_nominal_footprint_edge_gap_mm": minimum_cutout_gap,
        "position_tolerance_worst_case_qc_max_to_min": {
            "tolerance_unit": "um",
            "independent_xy_error_envelope": tolerance_ratios,
            "semantics": (
                "deterministic corner envelope from the local 2-D FEM field interpolation; "
                "excludes linewidth, film, cavity-machining, and model-form uncertainty"
            ),
        },
        "placements": placements,
        "qc_predictions": qc_predictions,
        "fem_calibration_provenance": calibration.get("provenance", {}),
        "validity": {
            "cavity_pin_response": (
                "historical direct 3-D FEM pin-Qc anchors with current-depth frequency alignment"
            ),
            "field_position_response": (
                "bilinear interpolation of direct current-depth 3-D FEM field samples"
            ),
            "resonator_external_qc": (
                "coupled-mode conversion, not a conformal micron-to-centimeter full-wave solve; "
                "absolute values scale as 1/g^2"
            ),
            "relative_position_equalization": (
                "more robust than the absolute Qc calibration but still inherits cavity-model "
                "and fabrication uncertainty"
            ),
        },
    }
