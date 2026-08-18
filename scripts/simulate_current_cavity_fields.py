#!/usr/bin/env python3
"""Solve the configured silicon-loaded cavity field profiles.

This is a closed-PEC 3-D NGSolve eigenmode calculation used for relative
resonator placement.  SMA bores and pins are intentionally excluded here;
their external-Q response is calibrated separately.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import time

import numpy as np


C0_MM_PER_NS = 299.792458


@dataclass(frozen=True)
class Geometry:
    straight_length_mm: float = 15.0
    cavity_width_mm: float = 4.8
    end_radius_mm: float = 2.4
    cavity_height_mm: float = 18.0
    split_from_floor_mm: float = 6.7
    chip_size_mm: float = 5.05
    chip_thickness_mm: float = 0.50
    chip_center_abs_x_mm: float = 3.525
    silicon_epsr: float = 11.68

    @property
    def chip_z_min_mm(self) -> float:
        return self.split_from_floor_mm - self.chip_thickness_mm


@dataclass(frozen=True)
class Settings:
    maxh_mm: float = 0.9
    chip_maxh_mm: float = 0.20
    order: int = 1
    curve_order: int = 3
    iterations: int = 20
    sample_z_mm: float = 6.71
    sample_patch_offset_mm: float = 0.12
    sample_extent_mm: float = 1.20
    sample_step_mm: float = 0.40


def build_mesh(geometry: Geometry, settings: Settings, mesh_path: Path):
    from netgen.occ import Box, Cylinder, Glue, OCCGeometry, Z
    from ngsolve import Mesh

    half_straight = 0.5 * geometry.straight_length_mm
    half_width = 0.5 * geometry.cavity_width_mm
    central = Box(
        (-half_straight, -half_width, 0.0),
        (half_straight, half_width, geometry.cavity_height_mm),
    )
    left = Cylinder(
        (-half_straight, 0.0, 0.0),
        Z,
        r=geometry.end_radius_mm,
        h=geometry.cavity_height_mm,
    )
    right = Cylinder(
        (half_straight, 0.0, 0.0),
        Z,
        r=geometry.end_radius_mm,
        h=geometry.cavity_height_mm,
    )
    cavity = central + left + right
    cavity.faces.name = "pec"

    chips = []
    for center_x in (-geometry.chip_center_abs_x_mm, geometry.chip_center_abs_x_mm):
        chip = Box(
            (
                center_x - 0.5 * geometry.chip_size_mm,
                -half_width,
                geometry.chip_z_min_mm,
            ),
            (
                center_x + 0.5 * geometry.chip_size_mm,
                half_width,
                geometry.split_from_floor_mm,
            ),
        )
        chip.mat("silicon")
        chip.faces.name = "silicon_interface"
        chip.faces.maxh = settings.chip_maxh_mm
        chips.append(chip)
    air = cavity - chips[0] - chips[1]
    air.mat("air")
    shape = Glue([air, *chips])
    ngmesh = OCCGeometry(shape).GenerateMesh(
        maxh=settings.maxh_mm,
        curvaturesafety=2.0,
        grading=0.25,
    )
    ngmesh.Curve(settings.curve_order)
    mesh_path.parent.mkdir(parents=True, exist_ok=True)
    ngmesh.Save(str(mesh_path))
    return Mesh(ngmesh)


def _pinvit(stiffness, mass, preconditioner, *, iterations: int):
    from ngsolve import InnerProduct, Matrix, MultiVector, Vector

    work = stiffness.CreateRowVector()
    eigenvectors = MultiVector(work, 1)
    trial = MultiVector(work, 2)
    trial[0].SetRandom()
    eigenvectors[:] = preconditioner * trial[:1]
    eigenvalues = Vector([1.0])
    history = []
    for iteration in range(iterations):
        trial[:1] = stiffness * eigenvectors - (mass * eigenvectors).Scale(eigenvalues)
        trial[1:2] = preconditioner * trial[:1]
        trial[:1] = eigenvectors
        trial.Orthogonalize(mass)
        small_a = np.asarray(InnerProduct(trial, stiffness * trial).NumPy()).copy()
        small_m = np.asarray(InnerProduct(trial, mass * trial).NumPy()).copy()
        small_a = 0.5 * (small_a + small_a.T)
        small_m = 0.5 * (small_m + small_m.T)
        try:
            chol = np.linalg.cholesky(small_m)
        except np.linalg.LinAlgError:
            small_m += np.eye(2) * 1.0e-12
            chol = np.linalg.cholesky(small_m)
        reduced = np.linalg.solve(chol, small_a)
        reduced = np.linalg.solve(chol, reduced.T).T
        reduced = 0.5 * (reduced + reduced.T)
        values, basis = np.linalg.eigh(reduced)
        generalized = np.linalg.solve(chol.T, basis)
        eigenvalues = Vector([float(values[0])])
        eigenvectors[:] = trial * Matrix(generalized[:, :1])
        history.append(float(values[0]))
        print(f"PINVIT {iteration + 1:02d}: {values[0]:.10g}", flush=True)
    return float(eigenvalues[0]), eigenvectors[0], history


def solve_fundamental(mesh, geometry: Geometry, settings: Settings):
    from ngsolve import (
        BilinearForm,
        GridFunction,
        HCurl,
        IdentityMatrix,
        InnerProduct,
        Integrate,
        Preconditioner,
        TaskManager,
        curl,
        dx,
    )

    epsr = mesh.MaterialCF({"air": 1.0, "silicon": geometry.silicon_epsr})
    space = HCurl(mesh, order=settings.order, dirichlet="pec")
    u, v = space.TnT()
    stiffness = BilinearForm(InnerProduct(curl(u), curl(v)) * dx)
    mass = BilinearForm(epsr * InnerProduct(u, v) * dx)
    augmented = BilinearForm(
        (InnerProduct(curl(u), curl(v)) + epsr * InnerProduct(u, v)) * dx
    )
    preconditioner = Preconditioner(
        augmented, "direct", inverse="sparsecholesky"
    )
    started = time.perf_counter()
    with TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        augmented.Assemble()
        gradient, h1_space = space.CreateGradient()
        gradient_transpose = gradient.CreateTranspose()
        h1_mass = gradient_transpose @ mass.mat @ gradient
        h1_inverse = h1_mass.Inverse(
            freedofs=h1_space.FreeDofs(), inverse="sparsecholesky"
        )
        projector = IdentityMatrix() - gradient @ h1_inverse @ gradient_transpose @ mass.mat
        projected_preconditioner = projector @ preconditioner.mat
        eigenvalue, vector, history = _pinvit(
            stiffness.mat,
            mass.mat,
            projected_preconditioner,
            iterations=settings.iterations,
        )
    solve_seconds = time.perf_counter() - started

    field = GridFunction(space)
    field.vec.data = vector
    energy = float(np.real(Integrate(epsr * InnerProduct(field, field), mesh)))
    silicon_energy = float(
        np.real(
            Integrate(
                epsr * InnerProduct(field, field),
                mesh,
                definedon=mesh.Materials("silicon"),
            )
        )
    )
    field.vec.data *= 1.0 / math.sqrt(max(energy, 1.0e-30))
    frequency_ghz = C0_MM_PER_NS * math.sqrt(max(eigenvalue, 0.0)) / (2.0 * math.pi)
    return field, {
        "eigenvalue_per_mm2": eigenvalue,
        "frequency_ghz": frequency_ghz,
        "electric_energy_silicon_fraction": silicon_energy / energy,
        "ndof": space.ndof,
        "solve_seconds": solve_seconds,
        "eigenvalue_last_relative_change": abs(history[-1] - history[-2])
        / max(abs(history[-1]), 1.0e-30),
        "pinvit_history": history,
    }


def sample_profile(mesh, field, geometry: Geometry, settings: Settings):
    count = round(2.0 * settings.sample_extent_mm / settings.sample_step_mm)
    local_values = [
        -settings.sample_extent_mm + index * settings.sample_step_mm
        for index in range(count + 1)
    ]
    if not any(abs(value) < 1.0e-12 for value in local_values):
        local_values.append(0.0)
    local_values = sorted(round(value, 12) for value in local_values)
    x_values = tuple(
        geometry.chip_center_abs_x_mm + value for value in local_values
    )
    y_values = tuple(local_values)
    offsets = (
        -settings.sample_patch_offset_mm,
        0.0,
        settings.sample_patch_offset_mm,
    )
    samples = []
    profile = []
    for x_abs in x_values:
        aggregate = []
        for x_sign in (-1.0, 1.0):
            for y_mm in y_values:
                ex2 = ey2 = ez2 = 0.0
                count = 0
                for dx_mm in offsets:
                    for dy_mm in offsets:
                        vector = np.asarray(
                            field(
                                mesh(
                                    x_sign * (x_abs + dx_mm),
                                    y_mm + dy_mm,
                                    settings.sample_z_mm,
                                )
                            ),
                            dtype=complex,
                        )
                        ex2 += float(abs(vector[0]) ** 2)
                        ey2 += float(abs(vector[1]) ** 2)
                        ez2 += float(abs(vector[2]) ** 2)
                        count += 1
                ex = math.sqrt(ex2 / count)
                ey = math.sqrt(ey2 / count)
                ez = math.sqrt(ez2 / count)
                in_plane = math.hypot(ex, ey)
                aggregate.append(in_plane)
                samples.append(
                    {
                        "x_mm": x_sign * x_abs,
                        "y_mm": y_mm,
                        "z_mm": settings.sample_z_mm,
                        "ex_rms_normalized": ex,
                        "ey_rms_normalized": ey,
                        "ez_rms_normalized": ez,
                        "e_in_plane_rms_normalized": in_plane,
                    }
                )
        profile.append(
            {
                "x_mm": x_abs,
                "e_in_plane_rms_normalized_mean": float(np.mean(aggregate)),
                "e_in_plane_rms_normalized_std": float(np.std(aggregate)),
                "contributing_samples": len(aggregate),
            }
        )
    return profile, samples


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--work-directory",
        type=Path,
        default=Path("build/current_cavity_field_fem"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/four_resonator_chip/current_cavity_field_fem.json"),
    )
    parser.add_argument("--maxh", type=float, default=0.9)
    parser.add_argument("--chip-maxh", type=float, default=0.20)
    parser.add_argument("--iterations", type=int, default=20)
    parser.add_argument(
        "--sample-extent",
        type=float,
        default=1.20,
        help="symmetric chip-local x/y field-grid extent in mm",
    )
    parser.add_argument(
        "--case",
        action="append",
        metavar="NAME=HEIGHT_MM",
        help=(
            "solve only the named case at the supplied cavity height; repeat for "
            "multiple cases (default: cavity_A, cavity_B, and cavity_C)"
        ),
    )
    parser.add_argument(
        "--update-existing",
        action="store_true",
        help="merge selected cases into an existing output JSON",
    )
    args = parser.parse_args()
    settings = Settings(
        maxh_mm=args.maxh,
        chip_maxh_mm=args.chip_maxh,
        iterations=args.iterations,
        sample_extent_mm=args.sample_extent,
    )
    default_cases = {
        "cavity_A": Geometry(cavity_height_mm=18.00),
        "cavity_B": Geometry(cavity_height_mm=17.77),
        "cavity_C": Geometry(cavity_height_mm=15.20),
    }
    if args.case:
        cases = {}
        for specification in args.case:
            name, separator, raw_height = specification.partition("=")
            name = name.strip()
            if not separator or not name:
                parser.error(f"invalid --case {specification!r}; expected NAME=HEIGHT_MM")
            try:
                height_mm = float(raw_height)
            except ValueError:
                parser.error(f"invalid cavity height in --case {specification!r}")
            if height_mm <= Geometry().split_from_floor_mm:
                parser.error("cavity height must exceed the 6.7 mm split plane")
            if name in cases:
                parser.error(f"duplicate --case name {name!r}")
            cases[name] = Geometry(cavity_height_mm=height_mm)
    else:
        cases = default_cases

    if args.update_existing:
        if not args.output.is_file():
            parser.error("--update-existing requires an existing --output JSON")
        payload = json.loads(args.output.read_text(encoding="utf-8"))
        if payload.get("schema_id") != "geerlings.current_cavity_field_fem.v1":
            parser.error("existing output has an incompatible schema")
        if not isinstance(payload.get("profiles"), dict) or not isinstance(
            payload.get("cases"), dict
        ):
            parser.error("existing output has no profiles/cases objects")
    else:
        payload = {
            "schema_id": "geerlings.current_cavity_field_fem.v1",
            "solver": "NGSolve H(curl) closed-PEC eigenproblem",
            "purpose": "relative in-plane field profile for position-only Qc equalization",
            "profiles": {},
            "cases": {},
            "validity": {
                "current_cavity_depths_included": True,
                "two_silicon_chips_included": True,
                "sma_bore_and_pin_included": False,
                "absolute_external_qc_computed": False,
            },
        }
    for name, geometry in cases.items():
        print(f"=== {name}: height {geometry.cavity_height_mm:.2f} mm ===", flush=True)
        mesh_path = args.work_directory / name / "mesh.vol"
        mesh = build_mesh(geometry, settings, mesh_path)
        field, result = solve_fundamental(mesh, geometry, settings)
        profile, samples = sample_profile(mesh, field, geometry, settings)
        payload["profiles"][name] = profile
        payload["cases"][name] = {
            "geometry": asdict(geometry),
            "settings": asdict(settings),
            "mesh": {
                "elements_3d": mesh.ne,
                "materials": list(mesh.GetMaterials()),
                "boundaries": sorted(set(mesh.GetBoundaries())),
                "raw_mesh": mesh_path.as_posix(),
            },
            "fundamental": result,
            "samples": samples,
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                name: {
                    "frequency_ghz": row["fundamental"]["frequency_ghz"],
                    "elements_3d": row["mesh"]["elements_3d"],
                    "ndof": row["fundamental"]["ndof"],
                    "last_relative_change": row["fundamental"][
                        "eigenvalue_last_relative_change"
                    ],
                }
                for name, row in payload["cases"].items()
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
