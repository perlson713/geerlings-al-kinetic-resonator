"""Windows-native Gmsh/NGSolve eigenmode baseline for the Fig. 1 resonator.

The generated mesh uses millimetres, a zero-thickness PEC aluminium sheet at
z=0, scalar sapphire permittivity, and a finite PEC outer box.  Aluminium film
thickness is recorded as metadata because it does not enter a PEC-sheet solve.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


C0_M_PER_S = 299_792_458.0
UM_TO_MM = 1.0e-3


@dataclass(frozen=True)
class RunSettings:
    config: Path
    output: Path
    mesh_path: Path
    metal_edge_um: float = 12.0
    max_cell_um: float = 120.0
    transition_um: float = 160.0
    simplify_um: float = 0.75
    modes: int = 6
    iterations: int = 24
    order: int = 1
    preconditioner: str = "direct"
    outer_boundary: str = "pec"
    aluminium_nominal_nm: float = 175.0
    aluminium_min_nm: float = 150.0
    aluminium_max_nm: float = 200.0
    overrides: tuple[str, ...] = ()
    gradient_inverse: str = "sparsecholesky"


def _metal_components(layout, simplify_um: float):
    from shapely.geometry import Polygon
    from shapely.ops import unary_union

    merged = unary_union([Polygon(poly.tuples()) for poly in layout.metal_polygons])
    if not merged.is_valid:
        merged = merged.buffer(0)
    components = list(merged.geoms) if merged.geom_type == "MultiPolygon" else [merged]
    components = [
        item.simplify(simplify_um, preserve_topology=True) for item in components
    ]
    components.sort(key=lambda item: item.area, reverse=True)
    if len(components) != 2:
        raise RuntimeError(
            "Expected exactly two metal components (ground and resonator), "
            f"found {len(components)}"
        )
    return components


def _add_occ_polygon(gmsh, polygon) -> int:
    loops: list[int] = []
    rings = [polygon.exterior, *polygon.interiors]
    for ring in rings:
        coords = list(ring.coords)[:-1]
        point_tags = [
            gmsh.model.occ.addPoint(float(x) * UM_TO_MM, float(y) * UM_TO_MM, 0.0)
            for x, y in coords
        ]
        line_tags = [
            gmsh.model.occ.addLine(point_tags[index], point_tags[(index + 1) % len(point_tags)])
            for index in range(len(point_tags))
        ]
        loops.append(gmsh.model.occ.addCurveLoop(line_tags))
    return gmsh.model.occ.addPlaneSurface(loops)


def generate_mesh(config, settings: RunSettings) -> dict[str, object]:
    import gmsh
    from geerlings_resonator.geometry import generate_layout

    layout = generate_layout(config.resonator)
    components = _metal_components(layout, settings.simplify_um)
    bounds = layout.domain_bounds.expanded(config.mesh.air_margin_xy_um)
    xmin = bounds.x_min * UM_TO_MM
    ymin = bounds.y_min * UM_TO_MM
    lx = bounds.width * UM_TO_MM
    ly = bounds.height * UM_TO_MM
    sub_t = config.stackup.substrate_thickness_um * UM_TO_MM
    air_t = config.stackup.air_thickness_um * UM_TO_MM

    settings.output.mkdir(parents=True, exist_ok=True)
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("General.Verbosity", 2)
        gmsh.model.add("geerlings_design_a")
        sapphire = gmsh.model.occ.addBox(xmin, ymin, -sub_t, lx, ly, sub_t)
        air = gmsh.model.occ.addBox(xmin, ymin, 0.0, lx, ly, air_t)
        metal_tools = [(2, _add_occ_polygon(gmsh, polygon)) for polygon in components]
        gmsh.model.occ.fragment(
            [(3, sapphire), (3, air)],
            metal_tools,
            removeObject=True,
            # Keep the sheet tools so Gmsh exports conformal internal PEC faces.
            removeTool=False,
        )
        gmsh.model.occ.synchronize()

        volumes = [tag for _, tag in gmsh.model.getEntities(3)]
        sapphire_tags: list[int] = []
        air_tags: list[int] = []
        for tag in volumes:
            z = gmsh.model.occ.getCenterOfMass(3, tag)[2]
            (sapphire_tags if z < -1.0e-9 else air_tags).append(tag)
        if not sapphire_tags or not air_tags:
            raise RuntimeError(
                f"Unable to classify volumes: sapphire={sapphire_tags}, air={air_tags}"
            )

        surface_info: list[tuple[int, float, tuple[int, ...], bool]] = []
        for _, tag in gmsh.model.getEntities(2):
            bbox = gmsh.model.getBoundingBox(2, tag)
            # OCC expands planar bounding boxes by roughly 1e-7 model units.
            at_z0 = abs(bbox[2]) < 1.0e-6 and abs(bbox[5]) < 1.0e-6
            upward, _ = gmsh.model.getAdjacencies(2, tag)
            area = gmsh.model.occ.getMass(2, tag)
            surface_info.append((tag, area, tuple(int(x) for x in upward), at_z0))

        candidates = {
            tag: area
            for tag, area, upward, at_z0 in surface_info
            if at_z0 and len(upward) >= 2
        }
        metal_tags: list[int] = []
        for component in components:
            expected = component.area * UM_TO_MM * UM_TO_MM
            if not candidates:
                raise RuntimeError("No conformal z=0 surface remains for a metal component")
            tag = min(candidates, key=lambda item: abs(candidates[item] - expected))
            error = abs(candidates[tag] - expected)
            if error > max(2.0e-8, 2.0e-4 * expected):
                raise RuntimeError(
                    f"Metal surface area mismatch: expected {expected:g} mm^2, "
                    f"closest {candidates[tag]:g} mm^2"
                )
            metal_tags.append(tag)
            del candidates[tag]

        outer_tags = [
            tag for tag, _, upward, _ in surface_info if len(upward) == 1
        ]
        interface_tags = [
            tag
            for tag, _, upward, at_z0 in surface_info
            if at_z0 and len(upward) >= 2 and tag not in metal_tags
        ]
        if not outer_tags or not interface_tags:
            raise RuntimeError(
                f"Surface classification failed: outer={outer_tags}, interface={interface_tags}"
            )

        groups = [
            (3, sapphire_tags, "sapphire"),
            (3, air_tags, "air"),
            (2, metal_tags, "metal_pec"),
            (2, interface_tags, "interface"),
            (2, outer_tags, "outer_pec"),
        ]
        for dimension, tags, name in groups:
            group = gmsh.model.addPhysicalGroup(dimension, tags)
            gmsh.model.setPhysicalName(dimension, group, name)

        distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(distance, "FacesList", metal_tags)
        gmsh.model.mesh.field.setNumber(distance, "Sampling", 200)
        threshold = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(threshold, "InField", distance)
        gmsh.model.mesh.field.setNumber(
            threshold, "SizeMin", settings.metal_edge_um * UM_TO_MM
        )
        gmsh.model.mesh.field.setNumber(
            threshold, "SizeMax", settings.max_cell_um * UM_TO_MM
        )
        gmsh.model.mesh.field.setNumber(threshold, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(
            threshold, "DistMax", settings.transition_um * UM_TO_MM
        )
        gmsh.model.mesh.field.setAsBackgroundMesh(threshold)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.5 * settings.metal_edge_um * UM_TO_MM)
        gmsh.option.setNumber("Mesh.MeshSizeMax", settings.max_cell_um * UM_TO_MM)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 0)
        start = time.perf_counter()
        gmsh.model.mesh.generate(3)
        mesh_seconds = time.perf_counter() - start
        node_count = len(gmsh.model.mesh.getNodes()[0])
        element_types, element_tags, _ = gmsh.model.mesh.getElements(3)
        tetrahedra = sum(len(tags) for tags in element_tags)
        gmsh.write(str(settings.mesh_path))
    finally:
        gmsh.finalize()

    summary = {
        "mesh_path": str(settings.mesh_path),
        "units": "mm",
        "nodes": node_count,
        "volume_elements": tetrahedra,
        "mesh_seconds": mesh_seconds,
        "materials": {"sapphire": sapphire_tags, "air": air_tags},
        "boundaries": {
            "metal_pec": metal_tags,
            "interface": interface_tags,
            "outer_pec": outer_tags,
        },
        "metal_components": [
            {
                "role": role,
                "area_um2": component.area,
                "vertices": len(component.exterior.coords) - 1,
            }
            for role, component in zip(("ground", "resonator"), components, strict=True)
        ],
        "domain_mm": {
            "xmin": xmin,
            "xmax": xmin + lx,
            "ymin": ymin,
            "ymax": ymin + ly,
            "zmin": -sub_t,
            "zmax": air_t,
        },
    }
    (settings.output / "mesh_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary


def pinvit_numpy(mata, matm, pre, *, num: int, maxit: int):
    """NGSolve PINVIT with NumPy for the tiny projected Ritz problems."""
    from ngsolve import InnerProduct, Matrix, MultiVector, Vector

    work = mata.CreateRowVector()
    uvecs = MultiVector(work, num)
    vecs = MultiVector(work, 2 * num)
    for vector in vecs[:num]:
        vector.SetRandom()
    uvecs[:] = pre * vecs[:num]
    eigenvalues = Vector(num * [1.0])
    history: list[list[float]] = []
    for iteration in range(maxit):
        vecs[:num] = mata * uvecs - (matm * uvecs).Scale(eigenvalues)
        vecs[num : 2 * num] = pre * vecs[:num]
        vecs[:num] = uvecs
        vecs.Orthogonalize(matm)
        small_a = np.asarray(InnerProduct(vecs, mata * vecs).NumPy()).copy()
        small_m = np.asarray(InnerProduct(vecs, matm * vecs).NumPy()).copy()
        small_a = 0.5 * (small_a + small_a.T)
        small_m = 0.5 * (small_m + small_m.T)
        try:
            chol = np.linalg.cholesky(small_m)
        except np.linalg.LinAlgError:
            small_m += np.eye(small_m.shape[0]) * 1.0e-12
            chol = np.linalg.cholesky(small_m)
        reduced = np.linalg.solve(chol, small_a)
        reduced = np.linalg.solve(chol, reduced.T).T
        reduced = 0.5 * (reduced + reduced.T)
        values, basis = np.linalg.eigh(reduced)
        generalized = np.linalg.solve(chol.T, basis)
        eigenvalues = Vector(values[:num])
        uvecs[:] = vecs * Matrix(generalized[:, :num])
        rates = [float(value) for value in values[:num]]
        history.append(rates)
        print(f"PINVIT {iteration + 1:02d}: " + ", ".join(f"{x:.7g}" for x in rates))
    return eigenvalues, uvecs, history


def save_xy_field_preview(mesh, field, layout, path: Path, title: str) -> None:
    """Sample |E| just above the chip and save a directly viewable PNG."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    bounds = layout.domain_bounds
    xs_mm = np.linspace(bounds.x_min * UM_TO_MM, bounds.x_max * UM_TO_MM, 181)
    ys_mm = np.linspace(bounds.y_min * UM_TO_MM, bounds.y_max * UM_TO_MM, 241)
    values = np.full((len(ys_mm), len(xs_mm)), np.nan)
    for iy, y_mm in enumerate(ys_mm):
        for ix, x_mm in enumerate(xs_mm):
            try:
                vector = np.asarray(field(mesh(float(x_mm), float(y_mm), 0.002)))
                values[iy, ix] = float(np.linalg.norm(vector))
            except Exception:
                continue
    vmax = float(np.nanmax(values))
    if not math.isfinite(vmax) or vmax <= 0.0:
        return
    vmin = max(vmax * 1.0e-4, float(np.nanmin(values[values > 0.0])))
    figure, axis = plt.subplots(figsize=(7.2, 8.0), constrained_layout=True)
    image = axis.imshow(
        values,
        origin="lower",
        extent=(bounds.x_min, bounds.x_max, bounds.y_min, bounds.y_max),
        cmap="inferno",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        interpolation="bilinear",
        aspect="equal",
    )
    for polygon in layout.metal_polygons:
        points = [*polygon.tuples(), polygon.tuples()[0]]
        axis.plot(
            [point[0] for point in points],
            [point[1] for point in points],
            color="cyan",
            linewidth=0.25,
            alpha=0.65,
        )
    axis.set(title=title, xlabel="x (um)", ylabel="y (um)")
    figure.colorbar(image, ax=axis, label="|E| (normalized eigenvector)")
    figure.savefig(path, dpi=180)
    plt.close(figure)


def solve_mesh(config, settings: RunSettings, mesh_summary: dict[str, object]):
    from netgen.read_gmsh import ReadGmsh
    from ngsolve import (
        BilinearForm,
        GridFunction,
        HCurl,
        IdentityMatrix,
        InnerProduct,
        Integrate,
        Mesh,
        Norm,
        Preconditioner,
        TaskManager,
        VTKOutput,
        curl,
        dx,
    )
    from geerlings_resonator.geometry import generate_layout

    mesh = Mesh(ReadGmsh(str(settings.mesh_path)))
    materials = set(mesh.GetMaterials())
    boundaries = set(mesh.GetBoundaries())
    if materials != {"air", "sapphire"}:
        raise RuntimeError(f"Unexpected materials from MSH2: {sorted(materials)}")
    required_boundaries = {"metal_pec", "interface", "outer_pec"}
    if not required_boundaries.issubset(boundaries):
        raise RuntimeError(
            f"Missing boundary names: {sorted(required_boundaries - boundaries)}; "
            f"found {sorted(boundaries)}"
        )
    epsr = mesh.MaterialCF(
        {"air": 1.0, "sapphire": config.stackup.sapphire_permittivity}
    )
    fes = HCurl(
        mesh,
        order=settings.order,
        dirichlet=(
            "metal_pec|outer_pec"
            if settings.outer_boundary == "pec"
            else "metal_pec"
        ),
    )
    u, v = fes.TnT()
    stiffness = BilinearForm(InnerProduct(curl(u), curl(v)) * dx)
    mass = BilinearForm(epsr * InnerProduct(u, v) * dx)
    augmented = BilinearForm(
        (InnerProduct(curl(u), curl(v)) + epsr * InnerProduct(u, v)) * dx
    )
    if settings.preconditioner == "direct":
        preconditioner = Preconditioner(
            augmented, "direct", inverse="sparsecholesky"
        )
    else:
        preconditioner = Preconditioner(augmented, settings.preconditioner)
    start = time.perf_counter()
    with TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        augmented.Assemble()
        gradient, h1_space = fes.CreateGradient()
        gradient_t = gradient.CreateTranspose()
        h1_mass = gradient_t @ mass.mat @ gradient
        h1_inverse = h1_mass.Inverse(
            freedofs=h1_space.FreeDofs(), inverse=settings.gradient_inverse
        )
        projector = (
            IdentityMatrix()
            - gradient @ h1_inverse @ gradient_t @ mass.mat
        )
        projected_preconditioner = projector @ preconditioner.mat
        eigenvalues, eigenvectors, history = pinvit_numpy(
            stiffness.mat,
            mass.mat,
            projected_preconditioner,
            num=settings.modes,
            maxit=settings.iterations,
        )
    solve_seconds = time.perf_counter() - start

    rows: list[dict[str, object]] = []
    layout = generate_layout(config.resonator)
    preview_written = False
    for index, (eigenvalue, eigenvector) in enumerate(
        zip(eigenvalues, eigenvectors, strict=True), start=1
    ):
        lam = float(eigenvalue)
        field = GridFunction(fes)
        field.vec.data = eigenvector
        frequency_ghz = (
            C0_M_PER_S * 1.0e3 * math.sqrt(max(lam, 0.0)) / (2.0 * math.pi) / 1.0e9
        )
        raw_residual = stiffness.mat * eigenvector - lam * (mass.mat * eigenvector)
        correction = projected_preconditioner * raw_residual
        correction_rel = float(correction.Norm()) / max(
            float(eigenvector.Norm()), 1.0e-30
        )
        eigenvalue_change = abs(history[-1][index - 1] - history[-2][index - 1]) / max(
            abs(history[-1][index - 1]), 1.0e-30
        )
        total_energy = float(
            Integrate(epsr * InnerProduct(field, field), mesh)
        )
        sapphire_energy = float(
            Integrate(
                epsr
                * InnerProduct(field, field),
                mesh,
                definedon=mesh.Materials("sapphire"),
            )
        )
        curl_energy = float(Integrate(InnerProduct(curl(field), curl(field)), mesh))
        rows.append(
            {
                "mode": index,
                "eigenvalue_per_mm2": lam,
                "frequency_ghz": frequency_ghz,
                "pinvit_correction_relative": correction_rel,
                "eigenvalue_last_relative_change": eigenvalue_change,
                "rayleigh_quotient_per_mm2": curl_energy / total_energy,
                "electric_energy_sapphire_fraction": sapphire_energy / total_energy,
                "classification": "physical" if lam > 1.0e-6 else "harmonic_zero",
            }
        )
        VTKOutput(
            ma=mesh,
            coefs=[Norm(field)],
            names=["E_norm"],
            filename=str(settings.output / f"mode_{index:02d}"),
            subdivision=0,
        ).Do()
        if not preview_written and 1.0 < frequency_ghz < 30.0:
            save_xy_field_preview(
                mesh,
                field,
                layout,
                settings.output / f"mode_{index:02d}_E_xy.png",
                f"Mode {index}: |E| at z=+2 um, {frequency_ghz:.4f} GHz",
            )
            preview_written = True

    with (settings.output / "modes.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "solver": "NGSolve 6.2.2606 HCurl / NumPy PINVIT",
        "formulation": "lossless Maxwell eigenproblem, zero-thickness PEC sheet",
        "outer_boundary": settings.outer_boundary.upper(),
        "coordinate_unit": "mm",
        "frequency_conversion": "f = c0*1000*sqrt(lambda)/(2*pi)",
        "aluminium": {
            "nominal_thickness_nm": settings.aluminium_nominal_nm,
            "provided_range_nm": [
                settings.aluminium_min_nm,
                settings.aluminium_max_nm,
            ],
            "model": "PEC sheet",
            "thickness_effect_in_this_model": "none",
        },
        "sapphire_relative_permittivity": config.stackup.sapphire_permittivity,
        "finite_box_warning": "PEC box; not an open/radiating boundary model",
        "configuration": {
            "source": settings.config.name,
            "overrides": list(settings.overrides),
            "resonator": asdict(config.resonator),
        },
        "basis_order": settings.order,
        "preconditioner": settings.preconditioner,
        "gradient_inverse": settings.gradient_inverse,
        "hcurl_dofs": fes.ndof,
        "solve_seconds": solve_seconds,
        "mesh": mesh_summary,
        "modes": rows,
        "pinvit_history": history,
    }
    (settings.output / "results.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("configs/design_a.toml"))
    parser.add_argument("--output", type=Path, default=Path("build/design_a/ngsolve"))
    parser.add_argument("--mesh", type=Path)
    parser.add_argument("--metal-edge-um", type=float, default=12.0)
    parser.add_argument("--max-cell-um", type=float, default=120.0)
    parser.add_argument("--transition-um", type=float, default=160.0)
    parser.add_argument("--simplify-um", type=float, default=0.75)
    parser.add_argument("--modes", type=int, default=6)
    parser.add_argument("--iterations", type=int, default=24)
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument(
        "--preconditioner", choices=("direct", "bddc", "local"), default="direct"
    )
    parser.add_argument(
        "--gradient-inverse",
        choices=("sparsecholesky", "pardisospd", "pardiso", "umfpack", "superlu"),
        default="sparsecholesky",
        help="direct inverse used by the H1 gradient-nullspace projector",
    )
    parser.add_argument("--outer-boundary", choices=("pec", "pmc"), default="pec")
    parser.add_argument("--reuse-mesh", action="store_true")
    parser.add_argument("--mesh-only", action="store_true")
    parser.add_argument(
        "--set",
        dest="overrides",
        action="append",
        default=[],
        metavar="SECTION.KEY=VALUE",
        help="repeatable typed TOML override, for example resonator.finger_count=10",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    from geerlings_resonator.config import load_project

    config_path = args.config.resolve()
    output = args.output.resolve()
    mesh_path = (args.mesh or output / "model.msh").resolve()
    settings = RunSettings(
        config=config_path,
        output=output,
        mesh_path=mesh_path,
        metal_edge_um=args.metal_edge_um,
        max_cell_um=args.max_cell_um,
        transition_um=args.transition_um,
        simplify_um=args.simplify_um,
        modes=args.modes,
        iterations=args.iterations,
        order=args.order,
        preconditioner=args.preconditioner,
        outer_boundary=args.outer_boundary,
        overrides=tuple(args.overrides),
        gradient_inverse=args.gradient_inverse,
    )
    config = load_project(config_path, overrides=args.overrides)
    if config.resonator.include_feedline:
        raise RuntimeError("This baseline expects include_feedline=false")
    if args.reuse_mesh:
        mesh_summary = json.loads(
            (mesh_path.parent / "mesh_summary.json").read_text(encoding="utf-8")
        )
    else:
        mesh_summary = generate_mesh(config, settings)
    print(json.dumps(mesh_summary, indent=2))
    if not args.mesh_only:
        result = solve_mesh(config, settings, mesh_summary)
        print(json.dumps(result["modes"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
