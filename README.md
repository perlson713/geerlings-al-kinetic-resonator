# Geerlings compact resonator: parametric layout and eigenfrequency

This repository reconstructs the compact resonator in Fig. 1 of Geerlings
*et al.*, *Applied Physics Letters* **100**, 192601 (2012), and provides a
reproducible path from parameterized layout to electromagnetic eigenfrequency.
The analyzed variant is an isolated aluminium resonator; the upper feedline and
all ports are omitted.

```text
TOML geometry -> GDS/SVG/STEP -> Gmsh tetrahedra -> NGSolve H(curl) PEC mode
                                                     |
                                                     v
                         low-temperature Al sheet kinetic inductance
                                                     |
                                                     v
                                  thickness-corrected resonance sweep
```

![Approved Fig. 1-topology 9 GHz layout](results/fig1_square_9ghz/layout_preview.png)

## Approved Fig. 1-topology near-square 9 GHz variant

The accepted isolated layout is configured in
[`configs/design_fig1_square_9ghz_al200.toml`](configs/design_fig1_square_9ghz_al200.toml),
with a 150 nm alternative in the adjacent config. It has 14 IDC fingers and
six meander turns. The even turn count preserves the paper's full-width top
return. The IDC bottom also follows Fig. 1 exactly: the lowest finger is tied
to the left bus, the first right-bus finger is one 15 um pitch above it, and
both buses extend to the common lower edge. The feedline and ports are absent.

The ground-cutout rectangle is **424 x 395 um** (width/height 1.0734; 29 um
side difference). The 5 um H(curl) PEC solve gives **9.111777769 GHz** on
270,428 tetrahedra and 632,082 H(curl) degrees of freedom. Applying the
documented low-temperature Al kinetic-inductance correction gives
**8.988895988 GHz at 150 nm** and **9.006908779 GHz at 200 nm**. These are
-11.104 MHz (-0.123%) and +6.909 MHz (+0.077%) from 9 GHz, respectively.

The preview, field plot, compact solver output, and thickness results are under
[`results/fig1_square_9ghz`](results/fig1_square_9ghz); the numerical details
are in [`docs/fig1_square_9ghz_result.md`](docs/fig1_square_9ghz_result.md).
The previous seven-turn search result remains under `results/square_9ghz` as
superseded history. Numerical mesh uncertainty remains larger than the
residual target errors, so the displayed digits are not a fabrication-accuracy
claim. The complete mesh and VTU fields are attached to GitHub Release
[`fig1-square-9ghz-v2`](https://github.com/perlson713/geerlings-al-kinetic-resonator/releases/tag/fig1-square-9ghz-v2).

## Main result

The converged zero-thickness PEC baseline is **10.416985 GHz**. Applying the
documented reduced-order aluminium kinetic-inductance model gives:

| Al thickness | Sheet inductance | Kinetic-corrected frequency |
|---:|---:|---:|
| 100 nm | 116.694 fH/sq | 10.174467 GHz |
| 150 nm | 82.767 fH/sq | 10.243224 GHz |
| 200 nm | 70.422 fH/sq | 10.268590 GHz |

The nominal change from 100 to 200 nm is **+94.123 MHz**. The complete 10 nm
grid, uncertainty envelopes, plot, and machine-readable metadata are under
[`results/kinetic_inductance_sweep_100_200nm`](results/kinetic_inductance_sweep_100_200nm).

This is a frozen-current reduced-order correction of a full-wave PEC mode, not
a surface-impedance full-wave eigensolve. The distinction, equations, and
uncertainty semantics are recorded in
[`docs/kinetic_inductance_sweep.md`](docs/kinetic_inductance_sweep.md).

## Reproduce

Python 3.12--3.14 is supported. Install the layout and analysis dependencies:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[cad,mesh,analysis,dev]"
```

Generate the layout:

```powershell
.\.venv\Scripts\geerlings-em.exe build -c configs\design_a.toml --step
```

Run the kinetic-inductance sweep after a PEC baseline has been selected:

```powershell
$env:PYTHONPATH="$PWD\src"
.\.venv\Scripts\python.exe -B scripts\kinetic_inductance_sweep.py
```

The Windows-native Gmsh/NGSolve baseline command and convergence study are in
[`docs/ngsolve_result_design_a.md`](docs/ngsolve_result_design_a.md). Palace
input generation remains available, but Palace itself is supported on Linux;
the results committed here use NGSolve on Windows.

## Parameterization

Edit [`configs/design_a.toml`](configs/design_a.toml), or apply dotted command
line overrides, for example:

```powershell
.\.venv\Scripts\geerlings-em.exe build -c configs\design_a.toml `
  --set resonator.gC=20 --set resonator.gL=5 --set resonator.w=10
```

The paper specifies Design-A values `gC=10 um`, `gL=20 um`, `gR=10 um`,
`w=5 um`, and modal target `Z0=100 ohm`. Finger length, meander span length,
substrate dimensions, dielectric constant, and enclosure dimensions were not
reported and are explicit reconstruction assumptions. The visible defaults
use 14 IDC fingers and 21 horizontal inductor spans.

The original paper used 200 nm niobium. Aluminium, the 100--200 nm sweep, and
the material model are choices for this follow-up analysis.

## Artifacts

- `results/baseline_pec/`: compact PEC convergence data and selected field plot.
- `results/kinetic_inductance_sweep_100_200nm/`: CSV, JSON, and sweep plot.
- `build/`: regenerable raw solver output; excluded from Git history.
- GitHub Release `em-results-v1`: archived raw NGSolve meshes and VTU fields.
- `provenance/`: LLM-oriented context, event log, and SHA-256 manifest.

Run the tests with:

```powershell
$env:PYTHONPATH="$PWD\src"
.\.venv\Scripts\python.exe -B -m unittest discover -s tests -v
```

## Interpretation limits

- The geometry is reconstructed from a schematic, not from the original mask.
- The PEC frequency has a conservative `+/-0.10 GHz` mesh envelope.
- The modal impedance is the configured 100 ohm target, not a field-extracted
  impedance; its uncertainty is represented by an effective-squares envelope.
- The kinetic model assumes a 20 mK low-temperature aluminium film and the
  published thickness fit. Process-specific `Tc`, sheet resistance, or a
  measured penetration depth should replace it for fabrication prediction.
- Finite geometric metal thickness and conductor-loss Q are not solved here.

## References

- [Geerlings et al. (2012)](https://doi.org/10.1063/1.4710520)
- [López-Núñez et al. aluminium thin-film characterization (2025)](https://doi.org/10.1088/1361-6668/adf360)
- [setupEM](https://github.com/VolkerMuehlhaus/setupEM)
- [Palace](https://awslabs.github.io/palace/stable/)
- [gdsfactory](https://gdsfactory.github.io/gdsfactory/)
- [build123d](https://github.com/gumyr/build123d)
- [NGSolve Maxwell eigenproblem tutorial](https://docu.ngsolve.org/latest/i-tutorials/unit-2.4-Maxwell/Maxwellevp.html)
