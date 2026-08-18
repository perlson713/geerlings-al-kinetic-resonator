# Eight Fig. 1-topology patterns at nominal 15 MHz spacing

Result date: 2026-08-18

## Frequency grid and geometry

This bank is the 15 MHz counterpart of the 10 MHz bank in
[`frequency_bank_8x.md`](frequency_bank_8x.md). Eight layouts are placed on a
uniform 15 MHz grid that is exactly centered on 9 GHz at Al 200 nm: four below
and four above the center. Because the count is even, 9.000 GHz itself is not
occupied and lies midway between P04 and P05, so the innermost pair sits at
9 GHz -/+ 7.5 MHz and the outermost pair at 9 GHz -/+ 52.5 MHz.

| ID | offset from 9 GHz | target at Al 200 nm | cutout width | PEC surrogate | predicted Al 150 nm | predicted Al 200 nm |
|---|---:|---:|---:|---:|---:|---:|
| P01 | -52.5 MHz | 8.9475 GHz | 426.790 um | 9.051725775 GHz | 8.929617174 GHz | 8.947516533 GHz |
| P02 | -37.5 MHz | 8.9625 GHz | 426.086 um | 9.066878680 GHz | 8.944573902 GHz | 8.962502039 GHz |
| P03 | -22.5 MHz | 8.9775 GHz | 425.382 um | 9.082031584 GHz | 8.959531357 GHz | 8.977488167 GHz |
| P04 | -7.5 MHz | 8.9925 GHz | 424.676 um | 9.097227537 GHz | 8.974532038 GHz | 8.992517497 GHz |
| P05 | +7.5 MHz | 9.0075 GHz | 423.972 um | 9.112380442 GHz | 8.989490959 GHz | 9.007504880 GHz |
| P06 | +22.5 MHz | 9.0225 GHz | 423.268 um | 9.127533346 GHz | 9.004450618 GHz | 9.022492894 GHz |
| P07 | +37.5 MHz | 9.0375 GHz | 422.564 um | 9.142686251 GHz | 9.019411017 GHz | 9.037481542 GHz |
| P08 | +52.5 MHz | 9.0525 GHz | 421.858 um | 9.157882204 GHz | 9.034414665 GHz | 9.052513412 GHz |

The width is quantized to 2 nm so that both half-width coordinates remain on
the 1 nm GDS grid while preserving exact left/right symmetry. Quantization
leaves every Al 200 nm nominal frequency within 0.019 MHz of its target. The
total span 8.9475 to 9.0525 GHz stays inside the 420 to 430 um width bracket
of the two-point PEC calibration, so no extrapolation is introduced relative
to the 10 MHz bank.

## Frozen topology

The generator enforces the same accepted constraints as the 10 MHz bank, and
the run aborts if any of them changes:

- six constant-radius semicircular meander folds;
- full-width horizontal top return;
- unchanged 25 um fold pitch and 12.5 um fold radius;
- 14 IDC fingers with 5 um trace and 10 um capacitor gap;
- lowest IDC finger connected to the left bus at centerline `y=2.5 um`;
- first right-bus finger one 15 um pitch above at `y=17.5 um`;
- both IDC buses extended to the common `y=0` lower edge;
- no feedline and no ports.

Only the horizontal straight-span length and therefore the cutout width vary.
The cutout height remains 395 um.

## Frequency model

The frequency model is unchanged: the linear two-point PEC calibration

```text
(width, PEC frequency) = (424 um, 9.111777769 GHz)
(width, PEC frequency) = (431 um, 8.961109683 GHz)
slope = -21.524012329 MHz/um
```

followed by the documented 20 mK Al kinetic-inductance correction evaluated
independently at 150 and 200 nm, with the target widths solved against the
200 nm prediction.

This process is a calibrated PEC surrogate plus reduced-order kinetic
correction, not eight independent full-wave surface-impedance eigensolves.
The nominal 15 MHz relative spacing is smaller than the conservative
+/-0.10 GHz absolute PEC mesh envelope. Fabrication targeting at this spacing
requires a process-specific frequency calibration or measured test structures.

## Artifacts and reproduction

Each `results/frequency_bank_8x_15mhz/patterns/Pxx/` directory contains:

- `layout.gds`: positive metal on layer 1/datatype 0;
- `layout.svg`: deterministic visual/vector layout;
- `stackup.xml`: setupEM/gds2palace stackup;
- `manifest.json`: geometry, topology, and frequency metadata.

The aggregate outputs are `frequency_table.csv`, `results.json`, and
`frequency_bank_preview.png`. The grid center, spacing, and count are read
from the config, so the same script produces both banks:

```bash
PYTHONPATH=src python -B scripts/generate_frequency_bank.py \
  --config configs/frequency_bank_8x_15mhz.toml
```

Re-running `configs/frequency_bank_8x.toml` with this code reproduces the
committed 10 MHz widths and frequencies; only inductor-length summation noise
below 1e-9 um differs.

The independent GDS read-back check used KLayout 0.30.10. All eight files have
one top cell, 337 metal polygons, no layers 201/202, and a bounding width equal
to the requested design width plus the 200 um ground margin. The GDS writer was
gdsfactory 9.45.0.

The chip-level placement of this bank, four equally pitched resonators on each
of two 5.05 mm square chips, is generated in the companion repository
`aluminum-cavity-chip-layouts` under `chip_layouts/equal_pitch_15mhz`.
