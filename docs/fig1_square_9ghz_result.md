# Approved Fig. 1-topology near-square 9 GHz aluminium resonator

Result date: 2026-07-30

## Frozen topology

```yaml
feedline: false
ports: false
trace_width_um: 5.0
gC_um: 10.0
gL_um: 20.0
gR_um: 10.0
IDC_fingers: 14
meander_turns: 6
cutout_width_um: 424.0
cutout_height_um: 395.0
aspect_width_over_height: 1.0734177215
inductor_centerline_length_um: 2812.8942919921856
```

The six-turn meander is the accepted replacement for the earlier seven-turn
search result. An even turn count preserves the full-width top return visible
in Geerlings et al. Fig. 1. Horizontal span length remains parametric; the
constant-radius semicircular folds, pitch, end transitions, and left/right
mirror symmetry are frozen.

The IDC lower edge was checked directly against Fig. 1 and is regression
locked. In centerline coordinates, the lowest finger is
`capacitor_finger_00_left` at `y=2.5 um`; it starts at the left bus. The first
right-bus finger is `capacitor_finger_01_right` at `y=17.5 um`, exactly one
15 um pitch higher. Both vertical buses extend to `y=0`. In conductor-outline
coordinates the two finger rectangles span `y=[0,5] um` and `y=[15,20] um`.

## Electromagnetic result

The selected first positive resonator mode was solved using a 5 um local-metal
Gmsh mesh and first-order NGSolve H(curl) elements. The final reproducibility
run used the 200 nm configuration as metadata; the electromagnetic conductor
is a zero-thickness PEC sheet, so 150 versus 200 nm does not change this PEC
eigenvalue.

| Quantity | Value |
|---|---:|
| PEC eigenfrequency | 9.111777769 GHz |
| tetrahedra | 270,428 |
| H(curl) DOFs | 632,082 |
| sapphire electric-energy fraction | 0.903081 |
| final PINVIT physical-mode relative change | 1.85e-12 |

The two preceding near-zero solutions are disconnected-conductor harmonic
fields. The selected physical mode is shown in
`results/fig1_square_9ghz/pec_fine/mode_03_E_xy.png`.

## Aluminium kinetic-inductance correction

The reduced-order correction uses 569.321903 branch squares, a 0.98 meander
energy fraction, a 100 ohm modal-impedance target, and the documented 20 mK Al
film model.

| Al thickness | sheet Lk | kinetic L | corrected frequency | error from 9 GHz |
|---:|---:|---:|---:|---:|
| 150 nm | 82.7665 fH/sq | 48.0824 pH | 8.988895988 GHz | -11.104 MHz (-0.123%) |
| 200 nm | 70.4219 fH/sq | 40.9109 pH | 9.006908779 GHz | +6.909 MHz (+0.077%) |

Both allowed thicknesses are within 0.13% of 9 GHz and straddle the target.
The conservative +/-0.10 GHz PEC mesh envelope dominates the residual target
errors.

## Numerical scope and reproduction

This is a full-wave zero-thickness PEC eigensolve followed by a frozen-current
sheet-kinetic-inductance correction. It is not a coupled surface-impedance
eigensolve, does not compute conductor-loss Q, and does not claim the displayed
digits as fabrication accuracy.

Inputs and compact artifacts:

- `configs/design_fig1_square_9ghz_al150.toml`
- `configs/design_fig1_square_9ghz_al200.toml`
- `configs/fig1_square_9ghz_kinetic.toml`
- `results/fig1_square_9ghz/tuning_candidates.csv`
- `results/fig1_square_9ghz/pec_fine/`
- `results/fig1_square_9ghz/kinetic/`

The complete selected solver output is published in GitHub Release
[`fig1-square-9ghz-v2`](https://github.com/perlson713/geerlings-al-kinetic-resonator/releases/tag/fig1-square-9ghz-v2)
as `raw-fig1-square-9ghz-ngsolve.tar.gz` (95,613,636 bytes), SHA-256
`F9FBD7BFF1CA26491ECEB6E17E841B0E96C7D09E0B7AAA365D29E60BFF25F177`.
