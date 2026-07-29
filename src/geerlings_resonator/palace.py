"""Create, validate, and execute Palace v0.17 eigenmode configurations."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import shutil
import subprocess
from typing import Any, Iterable

from .em import EigenmodeParameters
from .errors import ConfigurationError, ExternalToolError


JsonObject = dict[str, Any]


def _attribute_values(items: Iterable[dict[str, Any]]) -> list[int]:
    values: list[int] = []
    for item in items:
        for value in item.get("Attributes", []):
            ivalue = int(value)
            if ivalue not in values:
                values.append(ivalue)
    return values


def driven_to_eigenmode(
    driven: JsonObject,
    parameters: EigenmodeParameters,
    *,
    output: str = "output/eigenmode",
) -> JsonObject:
    """Convert a gds2palace Driven config without altering mesh attributes.

    gds2palace is deliberately used only for geometry/mesh generation.  This
    function replaces its driven solver block with the Palace v0.17 Eigenmode
    block and optionally promotes all conductor surfaces to ideal PEC.
    """

    parameters.validate()
    config = deepcopy(driven)
    required = {"Problem", "Model", "Domains", "Boundaries", "Solver"}
    missing = required - set(config)
    if missing:
        raise ConfigurationError(
            "Generated Palace config is missing: " + ", ".join(sorted(missing))
        )

    config["Problem"] = {
        "Type": "Eigenmode",
        "Verbose": int(config.get("Problem", {}).get("Verbose", 2)),
        "Output": output,
    }

    solver = config["Solver"]
    solver.pop("Driven", None)
    solver["Order"] = parameters.order
    solver.setdefault("Device", "CPU")
    solver["Eigenmode"] = {
        "N": parameters.modes,
        "Tol": parameters.tolerance,
        "Target": parameters.target_frequency_ghz,
        "TargetUpper": parameters.search_max_frequency_ghz,
        "Save": parameters.save_modes,
    }
    linear = solver.setdefault("Linear", {})
    linear.setdefault("Type", "Default")
    linear.setdefault("KSPType", "GMRES")
    linear["Tol"] = parameters.linear_tolerance
    linear["MaxIts"] = parameters.linear_max_iterations

    boundaries = config["Boundaries"]
    if not parameters.keep_passive_ports:
        boundaries.pop("LumpedPort", None)
        boundaries.pop("WavePort", None)

    if parameters.metal_model.lower() == "pec":
        conductor_items = boundaries.pop("Conductivity", [])
        conductor_attributes = _attribute_values(conductor_items)
        existing_pec = boundaries.get("PEC", {}).get("Attributes", [])
        pec_attributes: list[int] = []
        for value in [*existing_pec, *conductor_attributes]:
            ivalue = int(value)
            if ivalue not in pec_attributes:
                pec_attributes.append(ivalue)
        if not pec_attributes:
            raise ConfigurationError(
                "No metal surface attributes were found in Boundaries.Conductivity"
            )
        boundaries["PEC"] = {"Attributes": pec_attributes}
    elif not boundaries.get("Conductivity"):
        raise ConfigurationError(
            "metal_model='conductivity' requested, but no conductivity boundary exists"
        )

    validate_eigenmode_config(config)
    return config


def validate_eigenmode_config(config: JsonObject) -> None:
    """Perform fast local invariants before Palace's authoritative dry-run."""

    for section in ("Problem", "Model", "Domains", "Boundaries", "Solver"):
        if not isinstance(config.get(section), dict):
            raise ConfigurationError(f"Palace section {section!r} is required")
    if config["Problem"].get("Type") != "Eigenmode":
        raise ConfigurationError("Problem.Type must be Eigenmode")
    if not config["Model"].get("Mesh"):
        raise ConfigurationError("Model.Mesh is required")
    materials = config["Domains"].get("Materials", [])
    if not materials:
        raise ConfigurationError("At least one domain material is required")
    for material in materials:
        if not material.get("Attributes"):
            raise ConfigurationError("Every domain material needs Attributes")
    eigen = config["Solver"].get("Eigenmode")
    if not isinstance(eigen, dict):
        raise ConfigurationError("Solver.Eigenmode is required")
    if float(eigen.get("Target", 0)) <= 0:
        raise ConfigurationError("Solver.Eigenmode.Target must be positive")
    if int(eigen.get("N", 0)) <= 0:
        raise ConfigurationError("Solver.Eigenmode.N must be positive")
    if not any(
        key in config["Boundaries"]
        for key in ("PEC", "Conductivity", "Impedance")
    ):
        raise ConfigurationError("No conductor boundary condition is configured")


def load_json(path: str | Path) -> JsonObject:
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ConfigurationError(f"Unable to read Palace JSON {path}: {exc}") from exc


def write_json(path: str | Path, config: JsonObject) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return target


def palace_argv(
    config_path: str | Path,
    *,
    executable: str = "palace",
    mpi_processes: int = 1,
    dry_run: bool = False,
) -> list[str]:
    if mpi_processes <= 0:
        raise ConfigurationError("mpi_processes must be positive")
    args = [executable]
    if dry_run:
        args.append("--dry-run")
    else:
        args.extend(["-np", str(mpi_processes)])
    args.append(str(Path(config_path)))
    return args


def run_palace(
    config_path: str | Path,
    *,
    executable: str = "palace",
    mpi_processes: int = 1,
    dry_run: bool = False,
    cwd: str | Path | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run Palace without a shell and surface a concise actionable failure."""

    resolved = shutil.which(executable)
    if resolved is None and not Path(executable).is_file():
        raise ExternalToolError(
            f"Palace executable {executable!r} was not found. Install Palace v0.17 "
            "inside Linux/WSL or pass --palace-executable."
        )
    args = palace_argv(
        config_path,
        executable=resolved or executable,
        mpi_processes=mpi_processes,
        dry_run=dry_run,
    )
    completed = subprocess.run(
        args,
        cwd=cwd,
        check=False,
        text=True,
        stdout=None,
        stderr=None,
    )
    if completed.returncode:
        mode = "configuration dry-run" if dry_run else "eigenmode solve"
        raise ExternalToolError(
            f"Palace {mode} failed with exit code {completed.returncode}"
        )
    return completed
