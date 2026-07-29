"""Collect the Design-A NGSolve mesh-convergence and boundary checks."""

from __future__ import annotations

import csv
import json
from pathlib import Path


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    root = Path(__file__).resolve().parents[1]
    output = root / "build" / "design_a" / "ngsolve"
    cases = [
        ("coarse", 15.0),
        ("medium", 7.5),
        ("fine", 5.0),
        ("finer", 4.0),
    ]
    rows: list[dict[str, object]] = []
    previous: float | None = None
    for name, edge_um in cases:
        result = json.loads((output / name / "results.json").read_text(encoding="utf-8"))
        physical = next(mode for mode in result["modes"] if mode["classification"] == "physical")
        frequency = float(physical["frequency_ghz"])
        mesh = result["mesh"]
        rows.append(
            {
                "case": name,
                "metal_edge_um": edge_um,
                "nodes": mesh["nodes"],
                "tetrahedra": mesh["volume_elements"],
                "frequency_ghz": frequency,
                "change_from_previous_percent": (
                    "" if previous is None else 100.0 * (frequency - previous) / previous
                ),
            }
        )
        previous = frequency

    pmc = json.loads((output / "coarse_pmc" / "results.json").read_text(encoding="utf-8"))
    pmc_mode = next(mode for mode in pmc["modes"] if mode["classification"] == "physical")
    pec_coarse = float(rows[0]["frequency_ghz"])
    pmc_frequency = float(pmc_mode["frequency_ghz"])
    summary = {
        "recommended_baseline_frequency_ghz": rows[-1]["frequency_ghz"],
        "conservative_mesh_uncertainty_ghz": 0.1,
        "model": "Al 175 nm metadata; zero-thickness PEC sheet on eps_r=9.4 sapphire",
        "aluminium_range_nm": [150, 200],
        "thickness_dependence_in_pec_model": "none",
        "outer_boundary": "PEC baseline",
        "coarse_pec_frequency_ghz": pec_coarse,
        "coarse_pmc_frequency_ghz": pmc_frequency,
        "pec_pmc_difference_percent": 100.0 * (pmc_frequency - pec_coarse) / pec_coarse,
        "convergence": rows,
    }
    with (output / "convergence.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    figure, axis = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    axis.plot(
        [float(row["metal_edge_um"]) for row in rows],
        [float(row["frequency_ghz"]) for row in rows],
        "o-",
        color="#2369bd",
        linewidth=1.8,
    )
    axis.invert_xaxis()
    axis.grid(True, alpha=0.3)
    axis.set(
        xlabel="Local metal mesh size (um; finer to the right)",
        ylabel="Resonance frequency (GHz)",
        title="Design A eigenfrequency mesh convergence",
    )
    for row in rows:
        axis.annotate(
            f"{float(row['frequency_ghz']):.4f}",
            (float(row["metal_edge_um"]), float(row["frequency_ghz"])),
            xytext=(4, 7),
            textcoords="offset points",
            fontsize=8,
        )
    figure.savefig(output / "convergence.png", dpi=180)
    plt.close(figure)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
