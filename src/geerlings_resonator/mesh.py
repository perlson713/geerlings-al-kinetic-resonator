"""gds2palace bridge: GDS + XML stackup -> Gmsh mesh + Palace JSON."""

from __future__ import annotations

from contextlib import contextmanager
import os
from pathlib import Path
from typing import Iterator

from .em import EigenmodeParameters, MeshParameters
from .errors import ConfigurationError, ExternalToolError, OptionalDependencyError
from .palace import driven_to_eigenmode, load_json, write_json


@contextmanager
def _working_directory(path: Path) -> Iterator[None]:
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def build_mesh_and_eigenmode_config(
    gds_path: str | Path,
    stackup_path: str | Path,
    output_dir: str | Path,
    mesh: MeshParameters,
    eigenmode: EigenmodeParameters,
    *,
    model_name: str = "geerlings_fig1",
    left_port_layer: int | None = 201,
    right_port_layer: int | None = 202,
    metal_datatype: int = 0,
    port_datatype: int = 0,
) -> tuple[Path, Path]:
    """Build a tetrahedral mesh with gds2palace, then emit Eigenmode JSON.

    The generated Driven JSON remains next to the eigenmode JSON for setupEM
    inspection and for optional S-parameter validation.
    """

    mesh.validate()
    eigenmode.validate()
    try:
        from gds2palace import (  # type: ignore[import-not-found]
            gds_reader,
            simulation_setup,
            stackup_reader,
        )
    except ImportError as exc:
        raise OptionalDependencyError(
            "Meshing requires gds2palace and gmsh. Install with "
            "`pip install -e '.[mesh]'`."
        ) from exc

    gds = Path(gds_path).resolve()
    stackup = Path(stackup_path).resolve()
    target = Path(output_dir).resolve()
    if not gds.is_file():
        raise ExternalToolError(f"GDS file does not exist: {gds}")
    if not stackup.is_file():
        raise ExternalToolError(f"Stackup XML does not exist: {stackup}")
    target.mkdir(parents=True, exist_ok=True)

    materials, dielectrics, metals = stackup_reader.read_substrate(str(stackup))
    ports = simulation_setup.all_simulation_ports()
    if eigenmode.keep_passive_ports:
        if left_port_layer is None or right_port_layer is None:
            raise ConfigurationError(
                "keep_passive_ports=true requires GDS port marker layers"
            )
        # Resistive, unexcited feedline terminations contribute external Q in
        # an Eigenmode solve.  When disabled, port layers are not meshed and
        # the finite feedline is left open for a bare-frequency comparison.
        ports.add_port(
            simulation_setup.simulation_port(
                portnumber=1,
                voltage=0,
                port_Z0=50,
                source_layernum=left_port_layer,
                target_layername="Nb",
                direction="+y",
            )
        )
        ports.add_port(
            simulation_setup.simulation_port(
                portnumber=2,
                voltage=0,
                port_Z0=50,
                source_layernum=right_port_layer,
                target_layername="Nb",
                direction="+y",
            )
        )
    layer_numbers = metals.getlayernumbers()
    layer_numbers.extend(ports.portlayers)
    polygons = gds_reader.read_gds(
        str(gds),
        layer_numbers,
        purposelist=list(dict.fromkeys((metal_datatype, port_datatype))),
        metals_list=metals,
        preprocess=mesh.preprocess_gds,
        merge_polygon_size=0,
    )
    boundary = mesh.boundary.upper()
    settings = {
        "unit": 1.0e-6,
        "margin": mesh.dielectric_margin_um,
        "air_around": [
            mesh.air_margin_xy_um,
            mesh.air_margin_xy_um,
            mesh.air_margin_xy_um,
            mesh.air_margin_xy_um,
            mesh.air_below_um,
            mesh.air_above_um,
        ],
        # gds2palace uses fmax only to derive a wavelength mesh constraint.
        "fstart": eigenmode.target_frequency_ghz * 1.0e9,
        "fstop": eigenmode.search_max_frequency_ghz * 1.0e9,
        "fstep": 0.1e9,
        "refined_cellsize": mesh.refined_cell_size_um,
        "cells_per_wavelength": mesh.cells_per_wavelength,
        "meshsize_max": mesh.max_cell_size_um,
        "adaptive_mesh_iterations": mesh.adaptive_iterations,
        "order": eigenmode.order,
        "boundary": [boundary] * 6,
        "nogui": True,
        "simulation_ports": ports,
        "materials_list": materials,
        "dielectrics_list": dielectrics,
        "metals_list": metals,
        "layernumbers": layer_numbers,
        "allpolygons": polygons,
        "sim_path": str(target),
        "model_basename": model_name,
    }

    try:
        # Some gds2palace versions resolve auxiliary files against cwd.
        with _working_directory(target):
            driven_path_raw, data_dir_raw = simulation_setup.create_palace([], settings)
    except SystemExit as exc:
        raise ExternalToolError(
            f"gds2palace stopped while building the model (code {exc.code})"
        ) from exc
    except Exception as exc:
        raise ExternalToolError(f"gds2palace meshing failed: {exc}") from exc

    raw_config = Path(driven_path_raw)
    raw_data_dir = Path(data_dir_raw)
    candidates = [raw_config]
    if not raw_config.is_absolute():
        candidates.extend((target / raw_config, target / raw_data_dir / raw_config))
    driven_path = next((path for path in candidates if path.is_file()), candidates[0])
    if not driven_path.is_file():
        searched = ", ".join(str(path) for path in candidates)
        raise ExternalToolError(
            "gds2palace returned without a readable Driven config; checked "
            + searched
        )
    generated = load_json(driven_path)
    eigen_config = driven_to_eigenmode(
        generated, eigenmode, output=f"output/{model_name}_eigenmode"
    )
    eigen_path = write_json(target / "config_eigenmode.json", eigen_config)
    mesh_reference = Path(str(generated.get("Model", {}).get("Mesh", "")))
    mesh_path = (
        mesh_reference
        if mesh_reference.is_absolute()
        else driven_path.parent / mesh_reference
    )
    if not mesh_path.is_file():
        raise ExternalToolError(
            "gds2palace returned successfully but its Model.Mesh does not exist: "
            f"{mesh_path}"
        )
    return mesh_path, eigen_path
