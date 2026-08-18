#!/usr/bin/env python3
"""Evaluate direct-coupling suppression for the centered-square chip layout."""

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

from geerlings_resonator.config import load_project
from geerlings_resonator.coupling import pair_coupling_estimate, resonator_dipole


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _load_toml(path: Path) -> dict:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=PROJECT_ROOT / "configs" / "four_resonator_chip.toml",
    )
    parser.add_argument("--output-directory", type=Path)
    args = parser.parse_args()

    config = _load_toml(args.config)
    square = config["centered_square"]
    coupling = square["coupling"]
    chip = config["chip"]
    bank_path = PROJECT_ROOT / config["frequency_bank_results"]
    bank_config_path = PROJECT_ROOT / config["frequency_bank_config"]
    bank = json.loads(bank_path.read_text(encoding="utf-8"))
    bank_config = _load_toml(bank_config_path)
    design = load_project(PROJECT_ROOT / bank_config["base_design_config"])

    modal_impedance = float(coupling["modal_impedance_ohm"])
    descriptors = []
    for row in bank["patterns"]:
        parameters = design.resonator.with_updates(
            nominal_width_um=float(row["nominal_width_um"])
        )
        descriptors.append(
            resonator_dipole(
                str(row["pattern_id"]),
                parameters,
                frequency_ghz=float(row["predicted_al200_ghz"]),
                impedance_ohm=modal_impedance,
            )
        )

    selected_side_mm = float(square["side_mm"])
    effective_epsr = float(coupling["effective_relative_permittivity"])
    maximum_coupling_hz = float(coupling["maximum_pair_coupling_khz"]) * 1.0e3
    minimum_detuning_hz = min(
        abs(first.frequency_hz - second.frequency_hz)
        for index, first in enumerate(descriptors)
        for second in descriptors[index + 1 :]
    )

    maximum_width_mm = max(row["cutout_width_um"] for row in bank["patterns"]) / 1000.0
    maximum_height_mm = max(row["cutout_height_um"] for row in bank["patterns"]) / 1000.0
    usable_ground_half_mm = 0.5 * float(chip["size_mm"]) - float(
        chip["edge_clearance_mm"]
    )
    ground_strip_mm = float(coupling["minimum_outer_ground_strip_mm"])
    maximum_side_mm = 2.0 * min(
        usable_ground_half_mm - 0.5 * maximum_width_mm - ground_strip_mm,
        usable_ground_half_mm - 0.5 * maximum_height_mm - ground_strip_mm,
    )
    if selected_side_mm > maximum_side_mm + 1.0e-12:
        raise ValueError("selected centered-square side exceeds the ground-clearance limit")

    sweep_start = float(coupling["sweep_start_side_mm"])
    sweep_step = float(coupling["sweep_step_mm"])
    side_values = []
    value = sweep_start
    while value < selected_side_mm - 1.0e-12:
        side_values.append(value)
        value += sweep_step
    side_values.append(selected_side_mm)

    sweep_rows = []
    selected_worst = None
    for side_mm in side_values:
        estimates = [
            pair_coupling_estimate(
                first,
                second,
                separation_m=side_mm * 1.0e-3,
                effective_relative_permittivity=effective_epsr,
                electric_orientation_factor=float(
                    coupling["electric_orientation_factor"]
                ),
                magnetic_orientation_factor=float(
                    coupling["magnetic_orientation_factor"]
                ),
            )
            for index, first in enumerate(descriptors)
            for second in descriptors[index + 1 :]
        ]
        worst = max(estimates, key=lambda estimate: estimate.total_coupling_hz)
        mixing = worst.total_coupling_hz / minimum_detuning_hz
        dispersive_shift_hz = worst.total_coupling_hz**2 / minimum_detuning_hz
        row = {
            "square_side_mm": side_mm,
            "maximum_pair_coupling_khz": worst.total_coupling_hz / 1.0e3,
            "electric_component_khz": worst.electric_coupling_hz / 1.0e3,
            "magnetic_component_khz": worst.magnetic_coupling_hz / 1.0e3,
            "maximum_mixing_amplitude_percent": 100.0 * mixing,
            "maximum_dispersive_shift_khz": dispersive_shift_hz / 1.0e3,
            "worst_pair": f"{worst.first_pattern}-{worst.second_pattern}",
        }
        sweep_rows.append(row)
        if math.isclose(side_mm, selected_side_mm, abs_tol=1.0e-12):
            selected_worst = (worst, mixing, dispersive_shift_hz)
    if selected_worst is None:
        raise RuntimeError("selected side was not evaluated")
    worst, mixing, dispersive_shift_hz = selected_worst

    if worst.total_coupling_hz > maximum_coupling_hz:
        raise RuntimeError("selected square side does not meet the coupling threshold")
    if mixing > float(coupling["maximum_mixing_amplitude_percent"]) / 100.0:
        raise RuntimeError("selected square side does not meet the mixing threshold")
    if dispersive_shift_hz > float(coupling["maximum_dispersive_shift_khz"]) * 1.0e3:
        raise RuntimeError("selected square side does not meet the dispersive-shift threshold")

    output = args.output_directory or PROJECT_ROOT / config["output_directory"]
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "centered_square_coupling_sweep.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sweep_rows[0]))
        writer.writeheader()
        writer.writerows(sweep_rows)

    figure, axis = plt.subplots(figsize=(7.2, 4.5), dpi=190, layout="constrained")
    axis.semilogy(
        [row["square_side_mm"] for row in sweep_rows],
        [row["maximum_pair_coupling_khz"] for row in sweep_rows],
        color="#2878a8",
        linewidth=1.8,
    )
    axis.axhline(maximum_coupling_hz / 1.0e3, color="#b25b24", linestyle="--")
    axis.scatter(
        [selected_side_mm],
        [worst.total_coupling_hz / 1.0e3],
        color="#238b45",
        zorder=3,
        label=f"selected {selected_side_mm:.2f} mm",
    )
    axis.set_xlabel("centered-square side / nearest center spacing [mm]")
    axis.set_ylabel("conservative pair coupling estimate |J|/2pi [kHz]")
    axis.grid(alpha=0.2, which="both")
    axis.legend(frameon=False)
    figure.savefig(output / "centered_square_coupling_sweep.png")
    plt.close(figure)

    payload = {
        "schema_id": "geerlings.four_resonator_chip.centered_square_coupling.v1",
        "layout_variant": "centered_square",
        "selected_square_side_mm": selected_side_mm,
        "chip_local_centers_mm": [
            [-0.5 * selected_side_mm, -0.5 * selected_side_mm],
            [-0.5 * selected_side_mm, 0.5 * selected_side_mm],
            [0.5 * selected_side_mm, -0.5 * selected_side_mm],
            [0.5 * selected_side_mm, 0.5 * selected_side_mm],
        ],
        "geometry_clearance": {
            "maximum_allowed_side_mm": maximum_side_mm,
            "minimum_outer_ground_strip_mm": usable_ground_half_mm
            - 0.5 * selected_side_mm
            - 0.5 * max(maximum_width_mm, maximum_height_mm),
            "minimum_ground_cutout_edge_gap_mm": selected_side_mm
            - max(maximum_width_mm, maximum_height_mm),
        },
        "frequency_bank_minimum_detuning_mhz": minimum_detuning_hz / 1.0e6,
        "thresholds": {
            "maximum_pair_coupling_khz": maximum_coupling_hz / 1.0e3,
            "maximum_mixing_amplitude_percent": float(
                coupling["maximum_mixing_amplitude_percent"]
            ),
            "maximum_dispersive_shift_khz": float(
                coupling["maximum_dispersive_shift_khz"]
            ),
        },
        "selected_result": {
            **worst.to_dict(),
            "maximum_pair_coupling_khz": worst.total_coupling_hz / 1.0e3,
            "maximum_mixing_amplitude_percent": 100.0 * mixing,
            "maximum_dispersive_shift_khz": dispersive_shift_hz / 1.0e3,
        },
        "resonator_dipoles": [descriptor.to_dict() for descriptor in descriptors],
        "model": {
            "method": "far-field electric plus magnetic dipoles; magnitudes added",
            "effective_relative_permittivity": effective_epsr,
            "electric_orientation_factor": float(
                coupling["electric_orientation_factor"]
            ),
            "magnetic_orientation_factor": float(
                coupling["magnetic_orientation_factor"]
            ),
            "cancellation_credit": False,
            "exact_zero_claimed": False,
            "semantics": (
                "spacing-screening estimate for negligible hybridization; "
                "a simultaneous multi-resonator full-wave solve or measurement "
                "is required for sign-off"
            ),
        },
        "artifacts": [
            "centered_square_coupling_analysis.json",
            "centered_square_coupling_sweep.csv",
            "centered_square_coupling_sweep.png",
        ],
    }
    (output / "centered_square_coupling_analysis.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload["selected_result"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
