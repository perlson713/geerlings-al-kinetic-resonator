"""Command-line entry point for layout through Palace eigenmode results."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import re
import sys
from typing import Sequence

from .config import ProjectConfig, load_project
from .em import ideal_lumped_targets, write_stackup_xml
from .errors import ConfigurationError, GeerlingsError
from .exceptions import MissingOptionalDependencyError
from .gds import write_gds
from .geometry import LayoutGeometry, generate_layout
from .mesh import build_mesh_and_eigenmode_config
from .palace import run_palace
from .results import (
    filter_frequency_window,
    format_modes,
    read_eigenmodes,
    write_results_json,
)
from .setupem import launch_setupem, write_setupem_model
from .step import write_step
from .svg import write_svg


DEFAULT_CONFIG = Path("configs/design_a.toml")


def _safe_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_.")
    return label or "value"


def _project_from_args(args: argparse.Namespace) -> ProjectConfig:
    return load_project(args.config, overrides=getattr(args, "set_values", ()))


def _output_path(config: ProjectConfig, requested: str | None) -> Path:
    return Path(requested or config.project.output_directory).resolve()


def build_artifacts(
    config: ProjectConfig,
    output_dir: str | Path,
    *,
    include_step: bool = False,
) -> tuple[LayoutGeometry, dict[str, Path]]:
    """Create deterministic layout artifacts shared by every workflow stage."""

    target = Path(output_dir).resolve()
    target.mkdir(parents=True, exist_ok=True)
    layout = generate_layout(config.resonator)
    artifacts = {
        "gds": write_gds(layout, target / "layout.gds"),
        "svg": write_svg(layout, target / "layout.svg", show_etch_guides=True),
        "stackup": write_stackup_xml(target / "stackup.xml", config.stackup),
    }
    if config.resonator.include_ports:
        artifacts["setupem_model"] = write_setupem_model(
            target / "setupem_driven_model.py",
            gds_filename=artifacts["gds"].name,
            stackup_filename=artifacts["stackup"].name,
            mesh=config.mesh,
            eigenmode=config.eigenmode,
            model_name=f"{config.project.name}_driven",
            left_port_layer=config.resonator.left_port_layer,
            right_port_layer=config.resonator.right_port_layer,
            metal_datatype=config.resonator.metal_datatype,
            port_datatype=config.resonator.port_datatype,
        )
    if include_step:
        artifacts["step"] = write_step(layout, target / "model.step")

    lumped = ideal_lumped_targets(
        config.eigenmode.target_frequency_ghz, config.resonator.z0_ohm
    )
    manifest = config.manifest()
    manifest.update(
        {
            "paper": {
                "citation": "Geerlings et al., Appl. Phys. Lett. 100, 192601 (2012)",
                "doi": "10.1063/1.4710520",
                "figure": 1,
                "reported_values": ["gC", "gL", "gR", "w", "Z0", "Nb thickness"],
                "inferred_values": [
                    "IDC finger count",
                    "inductor span count",
                    "finger/span lengths",
                    "CPW and coupling dimensions",
                    "substrate and enclosure properties",
                ],
            },
            "topology": "meander inductor in parallel with a two-terminal IDC",
            "ideal_lumped_target": asdict(lumped),
            "layout_bounds_um": asdict(layout.domain_bounds),
            "artifacts": {name: path.name for name, path in artifacts.items()},
        }
    )
    manifest_path = target / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    artifacts["manifest"] = manifest_path
    return layout, artifacts


def mesh_artifacts(
    config: ProjectConfig,
    output_dir: str | Path,
    *,
    include_step: bool = False,
) -> dict[str, Path]:
    _, artifacts = build_artifacts(config, output_dir, include_step=include_step)
    mesh_path, eigen_path = build_mesh_and_eigenmode_config(
        artifacts["gds"],
        artifacts["stackup"],
        Path(output_dir) / "palace",
        config.mesh,
        config.eigenmode,
        model_name=config.project.name,
        left_port_layer=(
            config.resonator.left_port_layer
            if config.resonator.include_ports
            else None
        ),
        right_port_layer=(
            config.resonator.right_port_layer
            if config.resonator.include_ports
            else None
        ),
        metal_datatype=config.resonator.metal_datatype,
        port_datatype=config.resonator.port_datatype,
    )
    artifacts["mesh"] = mesh_path
    artifacts["eigenmode_config"] = eigen_path
    return artifacts


def simulate_artifacts(
    config: ProjectConfig,
    output_dir: str | Path,
    *,
    include_step: bool = False,
    palace_executable: str = "palace",
    dry_run_only: bool = False,
) -> dict[str, Path]:
    artifacts = mesh_artifacts(config, output_dir, include_step=include_step)
    eigen_config = artifacts["eigenmode_config"]
    run_palace(
        eigen_config.name,
        executable=palace_executable,
        dry_run=True,
        cwd=eigen_config.parent,
    )
    if dry_run_only:
        return artifacts
    run_palace(
        eigen_config.name,
        executable=palace_executable,
        mpi_processes=config.eigenmode.mpi_processes,
        cwd=eigen_config.parent,
    )
    result_path = (
        eigen_config.parent
        / "output"
        / f"{config.project.name}_eigenmode"
        / "eig.csv"
    )
    modes = read_eigenmodes(result_path)
    artifacts["eigenmodes"] = result_path
    artifacts["results_json"] = write_results_json(
        Path(output_dir) / "eigenmodes.json", modes
    )
    print(format_modes(modes))
    return artifacts


def _print_artifacts(artifacts: dict[str, Path]) -> None:
    width = max(len(name) for name in artifacts)
    for name, path in artifacts.items():
        print(f"{name.ljust(width)}  {path}")


def _cmd_build(args: argparse.Namespace) -> int:
    config = _project_from_args(args)
    output = _output_path(config, args.output)
    _, artifacts = build_artifacts(config, output, include_step=args.step)
    _print_artifacts(artifacts)
    return 0


def _cmd_mesh(args: argparse.Namespace) -> int:
    config = _project_from_args(args)
    output = _output_path(config, args.output)
    artifacts = mesh_artifacts(config, output, include_step=args.step)
    _print_artifacts(artifacts)
    return 0


def _cmd_simulate(args: argparse.Namespace) -> int:
    config = _project_from_args(args)
    output = _output_path(config, args.output)
    artifacts = simulate_artifacts(
        config,
        output,
        include_step=args.step,
        palace_executable=args.palace_executable,
        dry_run_only=args.dry_run_only,
    )
    _print_artifacts(artifacts)
    return 0


def _cmd_results(args: argparse.Namespace) -> int:
    modes = read_eigenmodes(args.eig_csv)
    modes = filter_frequency_window(modes, args.minimum_ghz, args.maximum_ghz)
    if not modes:
        print("No modes are inside the requested frequency window.")
        return 1
    print(format_modes(modes))
    if args.json:
        write_results_json(args.json, modes)
    return 0


def _cmd_lc(args: argparse.Namespace) -> int:
    config = _project_from_args(args)
    frequency = args.frequency_ghz or config.eigenmode.target_frequency_ghz
    target = ideal_lumped_targets(frequency, config.resonator.z0_ohm)
    print(f"f0 = {target.frequency_ghz:g} GHz")
    print(f"Z0 = {target.impedance_ohm:g} ohm")
    print(f"L  = {target.inductance_nh:.6g} nH")
    print(f"C  = {target.capacitance_ff:.6g} fF")
    return 0


def _cmd_setupem(args: argparse.Namespace) -> int:
    config = _project_from_args(args)
    output = _output_path(config, args.output)
    _, artifacts = build_artifacts(config, output, include_step=False)
    if "setupem_model" not in artifacts:
        raise ConfigurationError(
            "setupEM verification requires ports.include_ports=true"
        )
    launch_setupem(
        artifacts["gds"], artifacts["stackup"], executable=args.setupem_executable
    )
    return 0


def _cmd_sweep(args: argparse.Namespace) -> int:
    root_config = _project_from_args(args)
    root_output = _output_path(root_config, args.output)
    summaries: list[tuple[str, Path]] = []
    for raw_value in args.values:
        overrides = [*args.set_values, f"{args.parameter}={raw_value}"]
        config = load_project(args.config, overrides=overrides)
        label = f"{args.parameter.replace('.', '-')}-{_safe_label(raw_value)}"
        output = root_output / label
        if args.simulate:
            simulate_artifacts(
                config,
                output,
                include_step=False,
                palace_executable=args.palace_executable,
            )
        elif args.mesh:
            mesh_artifacts(config, output)
        else:
            build_artifacts(config, output)
        summaries.append((raw_value, output))
    for value, output in summaries:
        print(f"{args.parameter}={value}  {output}")
    return 0


def _add_config_options(parser: argparse.ArgumentParser, *, output: bool = True) -> None:
    parser.add_argument(
        "-c", "--config", type=Path, default=DEFAULT_CONFIG, help="Project TOML"
    )
    parser.add_argument(
        "--set",
        dest="set_values",
        action="append",
        default=[],
        metavar="SECTION.KEY=VALUE",
        help="Override a TOML value; repeatable",
    )
    if output:
        parser.add_argument("-o", "--output", help="Artifact directory override")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="geerlings-em",
        description="Parametric Geerlings Fig. 1 layout and Palace eigenmode workflow",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build", help="Write GDS, SVG, stackup and setupEM model")
    _add_config_options(build)
    build.add_argument("--step", action="store_true", help="Also export build123d STEP")
    build.set_defaults(func=_cmd_build)

    mesh = subparsers.add_parser("mesh", help="Build and create gds2palace Gmsh mesh")
    _add_config_options(mesh)
    mesh.add_argument("--step", action="store_true")
    mesh.set_defaults(func=_cmd_mesh)

    simulate = subparsers.add_parser(
        "simulate", help="Build, mesh, Palace dry-run, solve, and parse eig.csv"
    )
    _add_config_options(simulate)
    simulate.add_argument("--step", action="store_true")
    simulate.add_argument("--palace-executable", default="palace")
    simulate.add_argument(
        "--dry-run-only", action="store_true", help="Validate JSON but do not solve"
    )
    simulate.set_defaults(func=_cmd_simulate)

    results = subparsers.add_parser("results", help="Read a Palace eig.csv")
    results.add_argument("eig_csv", type=Path)
    results.add_argument("--minimum-ghz", type=float)
    results.add_argument("--maximum-ghz", type=float)
    results.add_argument("--json", type=Path)
    results.set_defaults(func=_cmd_results)

    lc = subparsers.add_parser("lc", help="Show ideal L/C synthesis targets")
    _add_config_options(lc, output=False)
    lc.add_argument("--frequency-ghz", type=float)
    lc.set_defaults(func=_cmd_lc)

    setupem = subparsers.add_parser("setupem", help="Build and open artifacts in setupEM")
    _add_config_options(setupem)
    setupem.add_argument("--setupem-executable", default="setupEM")
    setupem.set_defaults(func=_cmd_setupem)

    sweep = subparsers.add_parser("sweep", help="Generate or simulate a parameter sweep")
    _add_config_options(sweep)
    sweep.add_argument("parameter", help="Dotted parameter, e.g. resonator.gC")
    sweep.add_argument("values", nargs="+", help="TOML values")
    mode = sweep.add_mutually_exclusive_group()
    mode.add_argument("--mesh", action="store_true")
    mode.add_argument("--simulate", action="store_true")
    sweep.add_argument("--palace-executable", default="palace")
    sweep.set_defaults(func=_cmd_sweep)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (GeerlingsError, MissingOptionalDependencyError, OSError, ValueError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
