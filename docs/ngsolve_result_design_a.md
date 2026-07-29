# Design A NGSolve eigenmode result (Al variant)

Result date: 2026-07-23

This is the Windows-native electromagnetic baseline for the isolated Fig. 1
reconstruction in `configs/design_a.toml`.

## Model

- Feedline and ports: omitted
- Film: aluminium, nominal 175 nm; zero-thickness PEC electromagnetic model
- Substrate: 430 um sapphire, scalar relative permittivity 9.4
- Air above substrate: 300 um
- Outer box: 1.000 x 1.195 x 0.730 mm
- Discretization: first-order H(curl) tetrahedral finite elements
- Mesher/solver: Gmsh 4.15.2 and NGSolve 6.2.2606
- Eigensolver: gradient-nullspace Poisson projection plus PINVIT

The two near-zero eigenpairs are electrostatic harmonic fields caused by the
disconnected PEC components. They are not resonances. The first positive
eigenpair is the resonator mode; its electric field is localized in the IDC
gaps and around the meander.

## Mesh convergence

| Case | Local metal size | Nodes | Tetrahedra | Resonance | Change |
|---|---:|---:|---:|---:|---:|
| coarse | 15.0 um | 9,371 | 38,792 | 10.749728 GHz | -- |
| medium | 7.5 um | 45,523 | 212,684 | 10.563558 GHz | -1.732% |
| fine | 5.0 um | 68,934 | 288,488 | 10.474939 GHz | -0.839% |
| finer | 4.0 um | 105,368 | 452,993 | **10.416985 GHz** | -0.553% |

Recommended baseline: **10.42 +/- 0.10 GHz numerical mesh envelope**. This
does not include reconstruction or material-model uncertainty.

Changing the coarse outer wall from PEC to PMC shifted the selected mode from
10.749728 to 10.747728 GHz (-0.0186%), indicating low wall sensitivity for
the selected box.

## Thickness interpretation

The PEC-sheet Maxwell operator has no film-thickness parameter. Thickness
dependence is introduced separately through the low-temperature sheet kinetic
inductance in `scripts/kinetic_inductance_sweep.py`. The PEC solve cannot
predict conductor-loss Q.

## Reproduce

From the repository root:

```powershell
$ProjectRoot=(Resolve-Path .).Path
$env:PYTHONPATH=(Join-Path $ProjectRoot 'src')
& (Join-Path $ProjectRoot '.venv-ngsolve\Scripts\python.exe') -B `
  (Join-Path $ProjectRoot 'scripts\ngsolve_quick.py') `
  --config (Join-Path $ProjectRoot 'configs\design_a.toml') `
  --output (Join-Path $ProjectRoot 'build\design_a\ngsolve\finer') `
  --metal-edge-um 4 --max-cell-um 100 --transition-um 60 `
  --simplify-um 0.05 --modes 3 --iterations 25 --preconditioner bddc
```

Compact outputs are committed under `results/baseline_pec/`. Raw meshes and
VTU fields are attached to GitHub Release `em-results-v1`.
