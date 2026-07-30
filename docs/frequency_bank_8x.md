# Eight Fig. 1-topology patterns at nominal 10 MHz spacing

Result date: 2026-07-30

## Frequency grid and geometry

The bank contains eight layouts whose Al 200 nm nominal frequencies are
exactly centered on 9 GHz with 10 MHz target spacing. Because the count is
even, 9.000 GHz lies midway between P04 and P05.

| ID | target at Al 200 nm | cutout width | PEC surrogate | predicted Al 150 nm | predicted Al 200 nm |
|---|---:|---:|---:|---:|---:|
| P01 | 8.965 GHz | 425.968 um | 9.069418513 GHz | 8.947080925 GHz | 8.965013875 GHz |
| P02 | 8.975 GHz | 425.498 um | 9.079534799 GHz | 8.957066726 GHz | 8.975018819 GHz |
| P03 | 8.985 GHz | 425.030 um | 9.089608037 GHz | 8.967010358 GHz | 8.984981465 GHz |
| P04 | 8.995 GHz | 424.560 um | 9.099724322 GHz | 8.976996810 GHz | 8.994986966 GHz |
| P05 | 9.005 GHz | 424.090 um | 9.109840608 GHz | 8.986983589 GHz | 9.004992746 GHz |
| P06 | 9.015 GHz | 423.620 um | 9.119956894 GHz | 8.996970696 GHz | 9.014998808 GHz |
| P07 | 9.025 GHz | 423.150 um | 9.130073180 GHz | 9.006958133 GHz | 9.025005152 GHz |
| P08 | 9.035 GHz | 422.680 um | 9.140189466 GHz | 9.016945900 GHz | 9.035011778 GHz |

The width is quantized to 2 nm so that both half-width coordinates remain on
the 1 nm GDS grid while preserving exact left/right symmetry. Quantization
leaves every Al 200 nm nominal frequency within 0.019 MHz of its target.

## Frozen topology

Every pattern retains the accepted constraints:

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

The PEC frequency-versus-width rule is a linear local calibration using two
5 um local-metal NGSolve runs:

```text
(width, PEC frequency) = (424 um, 9.111777769 GHz)
(width, PEC frequency) = (431 um, 8.961109683 GHz)
slope = -21.524012329 MHz/um
```

For each generated geometry, the exact inductor centerline length determines
the effective branch-square count. The documented 20 mK Al kinetic-inductance
model is then applied independently at 150 and 200 nm. The eight target widths
are solved against the 200 nm prediction and quantized to the GDS grid.

This process is a calibrated PEC surrogate plus reduced-order kinetic
correction, not eight independent full-wave surface-impedance eigensolves.
The nominal 10 MHz relative spacing is smaller than the conservative
+/-0.10 GHz absolute PEC mesh envelope. Fabrication targeting at this spacing
requires a process-specific frequency calibration or measured test structures.

## Artifacts and reproduction

Each `results/frequency_bank_8x/patterns/Pxx/` directory contains:

- `layout.gds`: positive metal on layer 1/datatype 0;
- `layout.svg`: deterministic visual/vector layout;
- `stackup.xml`: setupEM/gds2palace stackup;
- `manifest.json`: geometry, topology, and frequency metadata.

The aggregate outputs are `frequency_table.csv`, `results.json`, and
`frequency_bank_preview.png`. Regenerate them with:

```powershell
$env:PYTHONPATH="$PWD\src"
.\.venv-ngsolve\Scripts\python.exe -B scripts\generate_frequency_bank.py `
  --config configs\frequency_bank_8x.toml
```

The independent GDS read-back check used KLayout 0.30.10. All eight files have
one top cell, 337 metal polygons, no layers 201/202, and a bounding width equal
to the requested design width plus the 200 um ground margin.

The complete CAD bundle is published in GitHub Release
[`frequency-bank-8x-v1`](https://github.com/perlson713/geerlings-al-kinetic-resonator/releases/tag/frequency-bank-8x-v1)
as `frequency-bank-8x-cad.zip` (424,603 bytes), SHA-256
`E6937EBBED0519AFB0376F4583BB086A9B651FB59B90E3A8276E7E8D17FBDC1F`.
