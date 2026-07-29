#!/usr/bin/env python3
"""Apply a low-temperature Al kinetic-inductance correction to a PEC mode."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import json
from pathlib import Path
import platform
import sys
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geerlings_resonator.config import load_project
from geerlings_resonator.geometry import generate_layout
from geerlings_resonator.kinetic import (
    AluminiumLondonModel,
    ReducedOrderResonator,
    evaluate_point,
    inclusive_thicknesses,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=PROJECT_ROOT / "configs" / "al_kinetic_sweep.toml",
        help="kinetic-inductance sweep TOML",
    )
    parser.add_argument(
        "--design-config",
        type=Path,
        default=PROJECT_ROOT / "configs" / "design_a.toml",
        help="layout TOML used to validate the frozen geometry",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="override the output directory from the sweep TOML",
    )
    return parser.parse_args()


def _load_toml(path: Path) -> dict:
    with path.resolve().open("rb") as handle:
        return tomllib.load(handle)


def _resonator(baseline: dict, reduction: dict, *, frequency_ghz: float, scale: float):
    return ReducedOrderResonator(
        pec_frequency_ghz=frequency_ghz,
        modal_impedance_ohm=float(baseline["modal_impedance_ohm"]),
        branch_squares=float(reduction["branch_squares"]) * scale,
        meander_energy_fraction=float(reduction["meander_energy_fraction"]),
    )


def _material(material: dict, *, scattering_factor: float):
    return AluminiumLondonModel(
        london_penetration_nm=float(material["london_penetration_nm"]),
        clean_coherence_length_nm=float(material["clean_coherence_length_nm"]),
        scattering_factor=scattering_factor,
    )


def _evaluate_grid(data: dict) -> tuple[list[dict], dict]:
    baseline = data["baseline"]
    reduction = data["geometry_reduction"]
    aluminium = data["aluminium_london"]
    sweep = data["sweep"]

    pec_ghz = float(baseline["pec_frequency_ghz"])
    pec_uncertainty_ghz = float(baseline["pec_mesh_uncertainty_ghz"])
    square_uncertainty = float(reduction["squares_relative_uncertainty"])
    factor = float(aluminium["scattering_factor"])
    factor_sigma = float(aluminium["scattering_factor_sigma"])
    thicknesses = inclusive_thicknesses(
        sweep["start_nm"], sweep["stop_nm"], sweep["step_nm"]
    )

    nominal_material = _material(aluminium, scattering_factor=factor)
    nominal_resonator = _resonator(
        baseline, reduction, frequency_ghz=pec_ghz, scale=1.0
    )
    model_low_material = _material(
        aluminium, scattering_factor=factor + factor_sigma
    )
    model_low_resonator = _resonator(
        baseline,
        reduction,
        frequency_ghz=pec_ghz,
        scale=1.0 + square_uncertainty,
    )
    model_high_material = _material(
        aluminium, scattering_factor=factor - factor_sigma
    )
    model_high_resonator = _resonator(
        baseline,
        reduction,
        frequency_ghz=pec_ghz,
        scale=1.0 - square_uncertainty,
    )
    total_low_resonator = _resonator(
        baseline,
        reduction,
        frequency_ghz=pec_ghz - pec_uncertainty_ghz,
        scale=1.0 + square_uncertainty,
    )
    total_high_resonator = _resonator(
        baseline,
        reduction,
        frequency_ghz=pec_ghz + pec_uncertainty_ghz,
        scale=1.0 - square_uncertainty,
    )

    rows: list[dict] = []
    for thickness in thicknesses:
        point = evaluate_point(
            thickness, material=nominal_material, resonator=nominal_resonator
        )
        model_low = evaluate_point(
            thickness,
            material=model_low_material,
            resonator=model_low_resonator,
        )
        model_high = evaluate_point(
            thickness,
            material=model_high_material,
            resonator=model_high_resonator,
        )
        total_low = evaluate_point(
            thickness,
            material=model_low_material,
            resonator=total_low_resonator,
        )
        total_high = evaluate_point(
            thickness,
            material=model_high_material,
            resonator=total_high_resonator,
        )
        row = point.as_dict()
        row.update(
            model_frequency_low_ghz=model_low.frequency_ghz,
            model_frequency_high_ghz=model_high.frequency_ghz,
            total_frequency_low_ghz=total_low.frequency_ghz,
            total_frequency_high_ghz=total_high.frequency_ghz,
        )
        rows.append(row)

    metadata = {
        "analysis_type": "reduced_order_kinetic_inductance_correction",
        "not_full_wave_kinetic_inductance_eigenanalysis": True,
        "pec_current_distribution_frozen": True,
        "finite_geometric_metal_thickness_included": False,
        "temperature_mk": float(aluminium["temperature_mk"]),
        "runtime": {
            "python": platform.python_version(),
            "matplotlib": matplotlib.__version__,
            "platform": platform.platform(),
        },
        "equations": {
            "penetration_depth": "lambda(d)=a*lambda_L*sqrt(xi_0/d)",
            "sheet_inductance": "Lk_square=mu0*lambda*coth(d/lambda)",
            "frequency_correction": "f=f_PEC/sqrt(1+Lk/Lg)",
            "geometric_inductance": "Lg=Z_mode/(2*pi*f_PEC)",
        },
        "baseline": baseline,
        "geometry_reduction": reduction,
        "aluminium_london": aluminium,
        "derived": {
            "effective_squares": nominal_resonator.effective_squares,
            "geometric_inductance_ph": nominal_resonator.geometric_inductance_h
            * 1.0e12,
            "capacitance_ff": nominal_resonator.capacitance_f * 1.0e15,
        },
        "uncertainty_semantics": {
            "model_band": "corner envelope of a +/- sigma and N_eff +/- relative uncertainty",
            "total_band": "model corners plus converged PEC frequency +/- mesh uncertainty",
            "not_statistical_confidence_interval": True,
        },
        "sources": [
            {
                "id": "Lopez-Nunez-2025",
                "url": "https://doi.org/10.1088/1361-6668/adf360",
                "use": "Al thickness-dependent penetration depth and sheet inductance",
            },
            {
                "id": "Geerlings-2012",
                "url": "https://doi.org/10.1063/1.4710520",
                "use": "resonator geometry and modal-impedance target",
            },
        ],
    }
    return rows, metadata


def _validate_geometry(design_config: Path, reduction: dict) -> dict:
    resolved_config = design_config.resolve()
    project = load_project(resolved_config)
    layout = generate_layout(project.resonator)
    centerline = layout.centerlines.inductor
    expected_length = float(reduction["expected_centerline_length_um"])
    expected_width = float(reduction["expected_trace_width_um"])
    length_error = centerline.length - expected_length
    width_error = centerline.width_um - expected_width
    if abs(length_error) > 1.0e-9 or abs(width_error) > 1.0e-12:
        raise RuntimeError(
            "design geometry differs from kinetic-reduction reference: "
            f"length error={length_error:.12g} um, width error={width_error:.12g} um"
        )
    if project.resonator.include_feedline or project.resonator.include_ports:
        raise RuntimeError("kinetic sweep requires the isolated no-feedline layout")
    try:
        config_label = resolved_config.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        config_label = resolved_config.name
    return {
        "design_config": config_label,
        "inductor_centerline_length_um": centerline.length,
        "inductor_trace_width_um": centerline.width_um,
        "feedline_present": layout.centerlines.feedline is not None,
        "port_count": len(layout.ports),
    }


def _write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_plot(path: Path, rows: list[dict], pec_frequency_ghz: float) -> None:
    thickness = [row["thickness_nm"] for row in rows]
    nominal = [row["frequency_ghz"] for row in rows]
    model_low = [row["model_frequency_low_ghz"] for row in rows]
    model_high = [row["model_frequency_high_ghz"] for row in rows]
    total_low = [row["total_frequency_low_ghz"] for row in rows]
    total_high = [row["total_frequency_high_ghz"] for row in rows]

    figure, axis = plt.subplots(figsize=(7.2, 4.5), constrained_layout=True)
    axis.fill_between(
        thickness,
        total_low,
        total_high,
        color="#b9d5ef",
        alpha=0.42,
        linewidth=0,
        label="model + PEC-mesh envelope",
    )
    axis.fill_between(
        thickness,
        model_low,
        model_high,
        color="#4c92ce",
        alpha=0.35,
        linewidth=0,
        label="kinetic-model envelope",
    )
    axis.plot(thickness, nominal, "o-", color="#075985", lw=2.0, label="nominal")
    axis.axhline(
        pec_frequency_ghz,
        color="#9a3412",
        ls="--",
        lw=1.4,
        label="zero-thickness PEC baseline",
    )
    axis.set_xlabel("Al thickness (nm)")
    axis.set_ylabel("Corrected resonance frequency (GHz)")
    axis.set_title("Kinetic-inductance thickness sweep")
    axis.grid(True, alpha=0.22)
    axis.legend(frameon=False, fontsize=8, loc="lower right")
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> int:
    args = _arguments()
    data = _load_toml(args.config)
    output = args.output or PROJECT_ROOT / data["sweep"]["output_directory"]
    if not output.is_absolute():
        output = PROJECT_ROOT / output
    output.mkdir(parents=True, exist_ok=True)

    geometry_validation = _validate_geometry(
        args.design_config, data["geometry_reduction"]
    )
    rows, metadata = _evaluate_grid(data)
    metadata["geometry_validation"] = geometry_validation
    metadata["sweep"] = {
        "count": len(rows),
        "start_nm": rows[0]["thickness_nm"],
        "stop_nm": rows[-1]["thickness_nm"],
        "frequency_change_100_to_200_mhz": (
            rows[-1]["frequency_ghz"] - rows[0]["frequency_ghz"]
        )
        * 1.0e3,
    }

    _write_csv(output / "sweep.csv", rows)
    with (output / "results.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {"metadata": metadata, "points": rows},
            handle,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")
    _write_plot(
        output / "frequency_vs_thickness.png",
        rows,
        float(data["baseline"]["pec_frequency_ghz"]),
    )
    print(json.dumps({"output": str(output), "points": len(rows)}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
