# Near-square 9 GHz aluminium resonator

Result date: 2026-07-29

## Selected geometry

```yaml
feedline: false
ports: false
trace_width_um: 5.0
gC_um: 10.0
gL_um: 20.0
gR_um: 10.0
IDC_fingers: 14
meander_turns: 7
cutout_width_um: 413.0
cutout_height_um: 420.0
aspect_width_over_height: 0.9833333333
side_mismatch_percent: 1.6666667
```

The preceding Design-A reconstruction was 300 x 495 um. The selected layout
widens the horizontal span and reduces the meander from 10 to 7 turns, while
preserving the 14-finger IDC and all paper-specified gaps and trace width.

## Electromagnetic result

The selected first positive resonator mode was solved with a 5 um local-metal
Gmsh mesh and first-order NGSolve H(curl) elements.

| Quantity | Value |
|---|---:|
| PEC eigenfrequency | 9.098137315 GHz |
| tetrahedra | 316,634 |
| H(curl) DOFs | 739,652 |
| sapphire electric-energy fraction | 0.902952 |

The two preceding near-zero modes are disconnected-conductor harmonic fields.
The selected mode is localized in the IDC gaps and meander, as shown in
`results/square_9ghz/pec_fine/mode_03_E_xy.png`.

## Aluminium kinetic-inductance correction

The effective meander reduction uses 569.588 branch squares, a 0.98 meander
energy fraction, and the same 20 mK literature Al model as the preceding sweep.

| Al thickness | sheet Lk | corrected frequency | error from 9 GHz |
|---:|---:|---:|---:|
| 150 nm | 82.7665 fH/sq | 8.975563322 GHz | -24.437 MHz (-0.272%) |
| 200 nm | 70.4219 fH/sq | 8.993531312 GHz | -6.469 MHz (-0.072%) |

The 200 nm film is the primary design point. Both allowed thicknesses are
within 0.3% of 9 GHz. A conservative +/-0.10 GHz mesh envelope dominates the
remaining target error.

## Numerical scope

The full-wave calculation is a zero-thickness PEC eigenmode followed by a
frozen-current sheet-kinetic-inductance correction. It is not a coupled
surface-impedance eigensolve. A 4 um mesh containing 502,971 tetrahedra was
generated, but the available Windows memory could not factor its gradient
projector; therefore the validated 5 um result and conservative mesh envelope
are retained.

Reproduction inputs:

- `configs/design_square_9ghz_al150.toml`
- `configs/design_square_9ghz_al200.toml`
- `configs/square_9ghz_kinetic.toml`

The complete selected solver output (mesh plus three VTU mode fields) is
published in [GitHub Release `square-9ghz-v1`](https://github.com/perlson713/geerlings-al-kinetic-resonator/releases/tag/square-9ghz-v1)
as asset `raw-square-9ghz-ngsolve.tar.gz` with
SHA-256 `83F91FF8FA03D05B36EA7EFDDBC3ABD30D3951EB4ECD466999C8747D1FB9FF90`.
