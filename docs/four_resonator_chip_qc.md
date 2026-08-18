# Four-resonator rectangle chips with Qc equalization

## Result

Each 5.05 mm chip carries four unrotated Geerlings resonators whose centers are
the four corners of an axis-aligned rectangle.  The layout keeps two modes
below 9 GHz and two above 9 GHz on each chip.  This mixed assignment is
intentional: putting all four sub-9-GHz modes on one chip and all four
over-9-GHz modes on the other left at least about 3.3% Qc spread under the same
rectangle and non-overlap constraints.

At the nominal resonator external-Qc target of `1.0e6`, the optimized
reduced-order result is:

| Cavity | SMA pin length (mm) | Pin protrusion (mm) | Resonator Qc min | Resonator Qc max | max/min |
|---|---:|---:|---:|---:|---:|
| A | 8.463572 | 0.963572 | 999,494.8 | 1,000,397.9 | 1.000904 |
| B | 8.930918 | 1.430918 | 998,436.4 | 1,001,049.5 | 1.002617 |
| C (15.20 mm) | 9.563236 | 2.063236 | 987,839.2 | 1,011,936.0 | 1.024393 |

Cases A and B retain the **0.262%** global deterministic model spread.  The
added 15.20 mm case is **2.439%** by itself and sets the max/min ratio over all
24 generated placements to **1.024393**.  Its solution reaches the 0.600 mm
minimum center separation; removing the pattern-column restriction produced
the same optimum.  The conservative minimum edge gap between ground cutouts
remains 0.174032 mm.

![Rectangle layout preview](../results/four_resonator_chip/four_per_chip_layout_preview.png)

## Centered-square companion version

A second modeled set fixes every chip to a square centered at the chip origin.
Its reference side is 4.00 mm, so the four chip-local centers are exactly
`(-2.00, -2.00)`, `(-2.00, +2.00)`, `(+2.00, -2.00)`, and
`(+2.00, +2.00)` mm.  Only the pattern-to-corner assignment is selected from
the extended 3-D field data.  The current cavity-A/B GDS files retain this
geometry; cavity C has the mask-only override documented below.

| Cavity | SMA pin length (mm) | Resonator Qc min | Resonator Qc max | max/min |
|---|---:|---:|---:|---:|
| A | 8.4844 | 669,372.9 | 1,604,804.7 | 2.3975 |
| B | 8.9573 | 616,509.2 | 1,602,432.1 | 2.5992 |
| C (15.20 mm, 4.00 mm reference) | 9.5932 | 663,438.4 | 1,525,164.7 | 2.2989 |

The centered-square global max/min ratio is **2.603051**.  This spread is the
tradeoff for moving the resonators far apart to suppress direct coupling while
keeping the square centered.  Its minimum ground-cutout gap is 3.574032 mm,
and the minimum cutout-to-outer-ground boundary strip is 0.212016 mm.

![Centered-square layout preview](../results/four_resonator_chip/centered_square_layout_preview.png)

The six companion GDS files contain `centered_square` in their names, for
example
`cavity_A_left_centered_square_below9_P01-P02_above9_P05-P06_4res.gds`.
The associated CSV, JSON, preview, manifest, and readback files use the same
`centered_square_` prefix in `results/four_resonator_chip/`.

### Cavity-C one-cell/two-layer fabrication override

The two current cavity-C `centered_square` GDS files replace only the modeled
4.00 mm cavity-C masks.  They preserve the same below/above-9-GHz pattern
partition, but use the following fabrication geometry:

- four centers at `(-0.75, -0.75)`, `(-0.75, +0.75)`, `(+0.75, -0.75)`, and
  `(+0.75, +0.75)` mm: 1.50 mm nearest center-to-center spacing
- exactly one flat cell per GDS, with no instances
- layer 1/0: the four isolated resonator metal patterns only
- layer 2/0: four corner L marks; each is a 100 x 100 um square minus an
  inward 75 x 75 um square, giving a 25 um leg width
- no chip-outline layer and no text-label layer

![Cavity-C two-layer preview](../results/four_resonator_chip/cavity_C_two_layer_1p5mm_preview.png)

`cavity_C_two_layer_1p5mm_verification.json` records direct KLayout readback of
the cell count, populated layers, shape counts, marker bounding boxes, areas,
and resolved resonator centers.  This is a mask-only override: the 4.00 mm Qc
table above and the 4.00 mm coupling result below do not describe these two
cavity-C GDS files.  For orientation, the existing reduced-order sweep at
1.50 mm estimates 1.646 MHz worst-pair coupling and 16.52% mixing; it therefore
does not meet the previous negligible-coupling thresholds.

### Direct resonator-resonator coupling criterion

The spacing selection uses a reduced-order far-field model.  Each reconstructed
resonator is reduced to its IDC electric dipole and closed meander-loop magnetic
dipole.  The worst electric orientation factor is used, electric and magnetic
coupling magnitudes are added, and no cancellation or ground-screening benefit
is credited.

At 4.00 mm nearest-neighbor spacing for the cavity-A/B masks and the reference
model, the worst pair is P01--P02:

- electric contribution: 82.77 kHz
- magnetic contribution: 4.03 kHz
- conservative total `|J|/2pi`: **86.79 kHz**
- mixing amplitude at the 9.96 MHz minimum detuning: **0.871%**
- dispersive frequency shift estimate: **0.756 kHz**

This meets the configured 100 kHz, 1%, and 1 kHz thresholds, respectively.
It means negligible hybridization under this model, not mathematically zero
coupling.  A simultaneous four-resonator full-wave solve or measurement is
still required for fabrication sign-off.  The complete sweep is in
`centered_square_coupling_analysis.json`,
`centered_square_coupling_sweep.csv`, and
`centered_square_coupling_sweep.png`.

## Optimized rectangle coordinates

Coordinates are relative to each chip center.  Every chip uses exactly the
Cartesian product `{x_low, x_high} x {y_low, y_high}`.

| Cavity/chip | x low (mm) | x high (mm) | y low (mm) | y high (mm) |
|---|---:|---:|---:|---:|
| A / left | 0.462089 | 1.170965 | -1.144980 | 0.190405 |
| A / right | -0.818393 | -0.213233 | -0.001215 | 1.015194 |
| B / left | 0.503932 | 1.180000 | -0.003953 | 0.993915 |
| B / right | -0.819876 | -0.219876 | -1.176560 | 0.502291 |
| C / left | 0.579176 | 1.179219 | -0.777434 | 0.015693 |
| C / right | -1.111399 | -0.511191 | -1.179140 | 0.866611 |

The low-frequency P01--P04 patterns occupy the stronger inner columns and the
high-frequency P05--P08 patterns occupy the weaker outer columns.  The
optimizer also chooses the upper/lower assignment in each column.

## GDS outputs and naming

The optimized-rectangle filenames and GDS text labels explicitly distinguish
modes below and above 9 GHz:

- `cavity_A_left_below9_P01-P02_above9_P05-P06_4res.gds`
- `cavity_A_right_below9_P03-P04_above9_P07-P08_4res.gds`
- `cavity_B_left_below9_P01-P02_above9_P05-P06_4res.gds`
- `cavity_B_right_below9_P03-P04_above9_P07-P08_4res.gds`
- `cavity_C_left_below9_P01-P02_above9_P05-P06_4res.gds`
- `cavity_C_right_below9_P03-P04_above9_P07-P08_4res.gds`

Layer 1/0 is positive aluminium metal, 100/0 is the 5.05 mm chip outline, and
101/0 contains labels such as `P01_BELOW9_ROT0` and `P05_ABOVE9_ROT0`.  All
six GDS files were read back with KLayout; each has one top cell, one outline,
four labels, and more than 1,000 metal shapes.

## Electromagnetic calculation and Qc semantics

The position response comes from a current-geometry, silicon-loaded, closed-PEC
3-D NGSolve H(curl) eigenmode calculation.  The sampled in-plane electric
field is interpolated in both chip-local x and y.

| Case | Cavity height (mm) | Closed-cavity mode (GHz) | Tetrahedra | H(curl) DOF | last eigenvalue change |
|---|---:|---:|---:|---:|---:|
| A | 18.00 | 9.759925 | 68,506 | 166,858 | 1.81e-4 |
| B | 17.77 | 9.797195 | 65,689 | 160,022 | 3.34e-4 |
| C | 15.20 | 10.328000 | 49,822 | 122,416 | 5.04e-6 |

The absolute resonator Qc is a hybrid calibration, not a conformal
micron-to-centimeter full-wave solution.  Historical direct 3-D FEM pin/Qc
anchors are frequency-aligned to the current closed-cavity solves, then
converted with

`Qc_res = fr * Qc_cavity * (fc - fr)^2 / (fc * g^2)`.

The nominal mean coupling is `g/2pi = 25 MHz`; absolute Qc scales as `1/g^2`.
The robust output is the relative position equalization.  The current field
solve does not include the SMA bore or pin, and the current
18.00/17.77/15.20 mm geometries have not received a new direct port-inclusive
pin sweep.  Case C explicitly uses the historical case-B pin-Qc curve as a
proxy, frequency-aligned to the current 15.20 mm closed-cavity solve.
Measurement or a conformal port-inclusive model is required before fabrication
sign-off.

## Reproduce

Generate the current-depth field calibration (requires NGSolve):

```powershell
$env:PYTHONPATH="$PWD\src"
python -B scripts\simulate_current_cavity_fields.py --maxh 0.9 --chip-maxh 0.20 --iterations 40
```

The coupling-suppressed square uses the separately committed extended field
grid `centered_square_cavity_field_fem.json`, which covers chip-local
`x,y = +/-2.2 mm`.  Reproduce it without replacing the rectangle calibration:

```powershell
$env:PYTHONPATH="$PWD\src"
python -B scripts\simulate_current_cavity_fields.py `
  --output results\four_resonator_chip\centered_square_cavity_field_fem.json `
  --case cavity_A=18.0 --case cavity_B=17.77 `
  --maxh 0.9 --chip-maxh 0.20 --iterations 40 --sample-extent 2.2
python -B scripts\simulate_current_cavity_fields.py `
  --output results\four_resonator_chip\centered_square_cavity_field_fem.json `
  --case cavity_C=15.2 --update-existing `
  --maxh 0.9 --chip-maxh 0.20 --iterations 80 --sample-extent 2.2
```

Optimize the rectangles and regenerate GDS/CSV/JSON/previews:

```powershell
$env:PYTHONPATH="$PWD\src"
python -B scripts\generate_four_resonator_chip.py
```

Then apply the approved cavity-C one-cell/two-layer override.  This command
overwrites only the two cavity-C `centered_square` GDS files and updates their
preview/readback metadata:

```powershell
$env:PYTHONPATH="$PWD\src"
python -B scripts\generate_cavity_c_two_layer_gds.py `
  --output-directory results\four_resonator_chip
```

Reproduce the direct-coupling spacing check:

```powershell
$env:PYTHONPATH="$PWD\src"
python -B scripts\analyze_centered_square_coupling.py
```

Machine-readable coordinates, Qc values, provenance, sensitivity envelopes,
and GDS readback are under `results/four_resonator_chip/`.
