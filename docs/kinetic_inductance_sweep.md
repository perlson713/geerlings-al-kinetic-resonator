# Aluminium kinetic-inductance thickness sweep

Result date: 2026-07-29

## Analysis contract

`analysis_type = reduced_order_kinetic_inductance_correction`

The converged NGSolve PEC eigenfrequency is corrected using a low-temperature
aluminium sheet inductance while freezing the PEC current distribution. This
is not a full-wave surface-impedance eigenanalysis.

The material equations are

```text
lambda(d) = a * lambda_L * sqrt(xi_0 / d)
Lk_square(d) = mu0 * lambda(d) * coth(d / lambda(d))
f(d) = f_PEC / sqrt(1 + Lk(d) / Lg)
```

with `lambda_L=15.7 nm`, `xi_0=1600 nm`, `a=1.26`, and `T=20 mK`, following
the thin-film aluminium characterization of
[López-Núñez et al.](https://doi.org/10.1088/1361-6668/adf360). The fit is
used over 100--200 nm, inside its reported thin-film range.

The layout reduction uses 618.962 branch squares, including terminal and bend
corrections. Dividing by the 0.98 meander-energy fraction gives 631.594
effective squares. The configured 100 ohm modal target and the PEC frequency
give `Lg=1527.841 pH` and `C=152.784 fF`.

## Nominal results

| Thickness (nm) | lambda (nm) | Lk sheet (fH/sq) | Lk mode (pH) | kinetic fraction | f (GHz) |
|---:|---:|---:|---:|---:|---:|
| 100 | 79.128 | 116.694 | 73.703 | 0.046020 | 10.174467 |
| 110 | 75.446 | 105.663 | 66.736 | 0.041852 | 10.196670 |
| 120 | 72.234 | 97.563 | 61.620 | 0.038768 | 10.213067 |
| 130 | 69.400 | 91.427 | 57.745 | 0.036419 | 10.225541 |
| 140 | 66.875 | 86.631 | 54.716 | 0.034574 | 10.235322 |
| 150 | 64.608 | 82.767 | 52.275 | 0.033083 | 10.243224 |
| 160 | 62.556 | 79.560 | 50.250 | 0.031842 | 10.249795 |
| 170 | 60.688 | 76.828 | 48.524 | 0.030782 | 10.255403 |
| 180 | 58.979 | 74.447 | 47.020 | 0.029857 | 10.260299 |
| 190 | 57.406 | 72.331 | 45.684 | 0.029033 | 10.264655 |
| 200 | 55.952 | 70.422 | 44.478 | 0.028288 | 10.268590 |

The 100--200 nm change is **+94.123 MHz**. Every corrected value remains below
the 10.416985 GHz PEC baseline because kinetic inductance adds positive modal
inductance.

## Envelope and limitations

The kinetic-model envelope combines `a=1.26 +/- 0.03` and effective squares
`+/-5%` at their corner values. The total envelope additionally combines the
PEC convergence estimate `+/-0.10 GHz`. These are deterministic sensitivity
envelopes, not statistical confidence intervals.

At 100 nm the kinetic-only envelope is 10.1545--10.1936 GHz and the total
envelope is 10.0594--10.2894 GHz. At 200 nm they are respectively
10.2577--10.2792 GHz and 10.1607--10.3766 GHz.

Finite geometric metal thickness is not included in the full-wave mesh. The
100 ohm impedance is a paper/configuration target rather than an impedance
extracted from the field solution. Fabrication prediction should replace the
generic Al fit with measured film `Tc`, normal sheet resistance, or penetration
depth and should rerun a coupled surface-impedance eigenanalysis.

Machine-readable results and exact assumptions are in
`results/kinetic_inductance_sweep_100_200nm/results.json`.
