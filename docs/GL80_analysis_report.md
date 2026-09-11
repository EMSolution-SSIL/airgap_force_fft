# GL80 Air-gap Force FFT Analysis Report

## 1. Purpose

This report summarizes the validation of the air-gap force FFT postprocessor
using the GL80 sample data.

The analysis starts from air-gap flux density:

- `Br(theta, t)`: radial air-gap flux density
- `Btheta(theta, t)`: circumferential air-gap flux density

and evaluates:

- Maxwell radial stress `Fr = (Br^2 - Btheta^2) / (2 mu0)`
- Maxwell tangential stress `Ftheta = Br * Btheta / mu0`
- torque from the circumferential stress
- time-order and space-order components by 2D FFT
- torque ripple DFT-bin orders
- dominant `Br/Btheta` mode pairs contributing to `Fr/Ftheta`

In this implementation, `Fr` and `Ftheta` are stress-like quantities in Pa, not
integrated force in N.

## 2. Run Configuration

The run conditions are recorded in [config.yaml](../config.yaml). The same merged
configuration is saved in [analysis_summary.json](../data/GL80/fft_output/analysis_summary.json).

Key conditions:

| Item | Value |
|---|---:|
| Input file | `GL80/motorGapB.json` |
| Torque reference | `GL80/transient_results.json` |
| Original angular span | 60 deg |
| Symmetricity | 1, anti-periodic |
| Sector repetition | 6 |
| Reconstructed theta points | 1512 |
| Time samples for torque comparison | 49 |
| Time samples for FFT | 48 |
| Air-gap radius | 0.04 m |
| Stack length | 0.01 m |

The JSON field `postData.gapB.symmetricity` controls sector reconstruction:

- `0`: periodic symmetry
- `1`: anti-periodic symmetry

For this GL80 case, the 60-degree model is anti-periodic. The reconstructed
360-degree flux density is therefore:

```text
sector 0:  +Br, +Btheta
sector 1:  -Br, -Btheta
sector 2:  +Br, +Btheta
...
```

This is important for flux-density spectra. A simple periodic copy would shift
the dominant flux density harmonics and obscure the correct relationship between
`Br/Btheta` and `Fr/Ftheta`.

The stack length is inferred as `2 * mean(position.z)`. For this sample,
`position.z = 0.005 m`, so the stack length is `0.01 m`.

The time series contains both the start point and the periodic end point. Torque
comparison uses all 49 samples. FFT-based outputs exclude the final time sample
to avoid counting the periodic endpoint twice. This is controlled by:

```yaml
analysis:
  exclude_time_endpoint_for_fft: true
```

## 3. Torque Validation

The integrated `Ftheta` torque is opposite in sign to the usual positive rotor
torque direction for this data. Therefore, the comparison plot uses:

- `-Ftheta integrated`
- `-Torque_Stator`
- `Torque_Rotor`

![Torque comparison](../data/GL80/fft_output/torque_comparison.png)

The agreement is excellent after applying the correct anti-periodic sector
reconstruction.

| Reference | Bias [Nm] | MAE [Nm] | Max abs. error [Nm] | Correlation |
|---|---:|---:|---:|---:|
| `-Torque_Stator` | -1.52e-6 | 1.52e-6 | 1.57e-6 | 1.0000 |
| `Torque_Rotor` | 2.93e-6 | 2.93e-6 | 3.62e-6 | 1.0000 |

Rotor-positive torque summary from the FFT evaluation samples:

| Quantity | Value |
|---|---:|
| Mean torque | 1.008339 Nm |
| Min torque | 0.993625 Nm |
| Max torque | 1.022243 Nm |
| Peak-to-peak ripple | 0.028618 Nm |
| RMS ripple | 0.009243 Nm |
| Peak-to-peak ripple / mean | 2.84 % |
| RMS ripple / mean | 0.92 % |

Detailed values are written to:

- [integrated_metrics.csv](../data/GL80/fft_output/integrated_metrics.csv)
- [torque_reference_comparison.csv](../data/GL80/fft_output/torque_reference_comparison.csv)
- [torque_ripple_summary.json](../data/GL80/fft_output/torque_ripple_summary.json)
- [torque_ripple_orders.csv](../data/GL80/fft_output/torque_ripple_orders.csv)

## 4. Harmonic Mode Maps

The following 4-panel plot compares the dominant 2D FFT components of the input
flux density and the resulting Maxwell stress.

The horizontal axis is DFT-bin time order `n`. The vertical axis is space
harmonic order `m`. Positive and negative space orders are folded into the
positive side for readability.

![4-panel harmonic map](../data/GL80/fft_output/airgap_2d_fft_4panel.png)

Single-field figures are also available:

- [Br_2d_fft.png](../data/GL80/fft_output/Br_2d_fft.png)
- [Btheta_2d_fft.png](../data/GL80/fft_output/Btheta_2d_fft.png)
- [Fr_2d_fft.png](../data/GL80/fft_output/Fr_2d_fft.png)
- [Ftheta_2d_fft.png](../data/GL80/fft_output/Ftheta_2d_fft.png)

## 5. Dominant Modes

### Br

| Frequency [Hz] | Space order `m` | Coefficient [T] |
|---:|---:|---:|
| -241.490 | 21 | 0.317623 |
| 241.490 | -21 | 0.317623 |
| 0 | -21 | 0.158845 |
| 0 | 21 | 0.158845 |
| -482.980 | 21 | 0.079503 |
| 482.980 | -21 | 0.079503 |
| 241.490 | 21 | 0.063641 |
| -241.490 | -21 | 0.063641 |

### Btheta

| Frequency [Hz] | Space order `m` | Coefficient [T] |
|---:|---:|---:|
| -241.490 | 21 | 0.085319 |
| 241.490 | -21 | 0.085319 |
| 0 | 21 | 0.042673 |
| 0 | -21 | 0.042673 |
| -241.490 | 57 | 0.037451 |
| 241.490 | -57 | 0.037451 |
| 482.980 | -63 | 0.029268 |
| -482.980 | 63 | 0.029268 |

### Fr

| Frequency [Hz] | Space order `m` | Coefficient [Pa] |
|---:|---:|---:|
| -241.490 | 42 | 38615.0 |
| 241.490 | -42 | 38615.0 |
| 0 | -36 | 25914.0 |
| 0 | 36 | 25914.0 |
| -482.980 | 42 | 19315.1 |
| 482.980 | -42 | 19315.1 |
| 724.470 | -84 | 16491.1 |
| -724.470 | 84 | 16491.1 |

### Ftheta

| Frequency [Hz] | Space order `m` | Coefficient [Pa] |
|---:|---:|---:|
| 0 | 36 | 22820.5 |
| 0 | -36 | 22820.5 |
| 241.490 | -42 | 22628.0 |
| -241.490 | 42 | 22628.0 |
| -482.980 | 42 | 11313.7 |
| 482.980 | -42 | 11313.7 |
| 0 | 72 | 11045.8 |
| 0 | -72 | 11045.8 |

The complete tables are written to:

- [Br_dominant_modes.csv](../data/GL80/fft_output/Br_dominant_modes.csv)
- [Btheta_dominant_modes.csv](../data/GL80/fft_output/Btheta_dominant_modes.csv)
- [Fr_dominant_modes.csv](../data/GL80/fft_output/Fr_dominant_modes.csv)
- [Ftheta_dominant_modes.csv](../data/GL80/fft_output/Ftheta_dominant_modes.csv)

## 6. Torque Ripple DFT-bin Orders

The table below is not a mechanical or electrical order-tracking result. It is a
DFT-bin order table for the sampled analysis window. Because the FFT uses 48
samples after excluding the periodic endpoint, the base frequency is:

```text
f_bin = 1 / (48 * dt) = 241.489896 Hz
```

The dominant positive-frequency torque ripple components are:

| DFT-bin order `n` | Frequency [Hz] | Single-sided amplitude [Nm] |
|---:|---:|---:|
| 8 | 1931.919 | 0.012563 |
| 4 | 965.960 | 0.003412 |
| 16 | 3863.838 | 0.001114 |
| 12 | 2897.879 | 0.000309 |
| 20 | 4829.798 | 0.000217 |
| 1 | 241.490 | 5.09e-7 |
| 2 | 482.980 | 2.46e-7 |
| 5 | 1207.449 | 2.04e-7 |
| 6 | 1448.939 | 1.64e-7 |
| 7 | 1690.429 | 1.61e-7 |

The largest torque ripple component is DFT-bin order `n=8`. The full table is
available in [torque_ripple_orders.csv](../data/GL80/fft_output/torque_ripple_orders.csv).

If a true mechanical or electrical order is needed, the next implementation step
should add rotor-angle-based order tracking or use a known rotational speed to
convert frequency into mechanical/electrical order.

## 7. Mode Pair Contribution Examples

Maxwell stress is generated by products of flux density components. Therefore,
large `Fr/Ftheta` modes can often be traced back to pairs of `Br/Btheta` modes.

For example, the dominant `Fr(n=1, m=-42)` mode is mainly produced by `Br * Br`
pairs:

| Rank | Product | Pair `(n, m)` | Contribution [Pa] | Fraction of target |
|---:|---|---|---:|---:|
| 1 | `Br*Br/(2mu0)` | `(0, -21) + (1, -21)` | 20074.5 | 52.0 % |
| 2 | `Br*Br/(2mu0)` | `(1, -21) + (0, -21)` | 20074.5 | 52.0 % |
| 3 | `Br*Br/(2mu0)` | `(2, -63) + (-1, 21)` | 7662.9 | 19.8 % |
| 4 | `Br*Br/(2mu0)` | `(-1, 21) + (2, -63)` | 7662.9 | 19.8 % |

The dominant `Ftheta(n=0, m=36)` mode is mainly produced by `Br * Btheta`
pairs:

| Rank | Product | Pair `(n, m)` | Contribution [Pa] | Fraction of target |
|---:|---|---|---:|---:|
| 1 | `Br*Btheta/mu0` | `(1, -21) + (-1, 57)` | 9465.9 | 41.5 % |
| 2 | `Br*Btheta/mu0` | `(-1, 21) + (1, 15)` | 5132.8 | 22.5 % |
| 3 | `Br*Btheta/mu0` | `(-1, 57) + (1, -21)` | 2669.2 | 11.7 % |
| 4 | `Br*Btheta/mu0` | `(0, -21) + (0, 57)` | 2367.0 | 10.4 % |

The full contribution table is written to
[mode_pair_contributions.csv](../data/GL80/fft_output/mode_pair_contributions.csv).

## 8. Interpretation

The result is consistent with the discussion memo, with one important
refinement: because the GL80 60-degree model is anti-periodic, the flux density
must be sign-reversed in adjacent sectors before the full 360-degree FFT is
evaluated.

After anti-periodic reconstruction, the input flux density has strong `m=21`
components. This is physically natural for a 42-pole machine because the pole
pair count is 21. Since Maxwell stress depends on quadratic products of flux
density, these input modes generate sum and difference components. This explains
why the stress spectra contain strong `m=42`, `m=36`, `m=72`, and `m=84`
components.

The torque validation supports the sign and scaling convention of the
implementation. The rotor-positive torque from `-Ftheta` agrees with the FEM
rotor torque to within a few micro-newton-meters for this sample.

The present GL80 sample is still reconstructed from a 60-degree sector. This is
suitable for checking periodic and anti-periodic harmonic content, torque, and
major force modes. It is not sufficient for evaluating true non-periodic
asymmetry such as eccentricity-induced `m=1` UMP. For that, the electromagnetic
model itself must include the asymmetry.

## 9. Recommended Next Steps

1. Add rotor-angle-based order tracking so torque ripple can be reported as
   mechanical or electrical order, not only DFT-bin order.
2. Apply the same workflow to another operating point, such as no-load, rated
   load, and high-current operation.
3. Compare how `Br/Btheta` dominant modes shift with current condition, and how
   that changes `Fr/Ftheta`.
4. Add a multi-case comparison report that overlays torque ripple and dominant
   force modes across cases.
5. For eccentricity or manufacturing variation studies, use a model that
   explicitly contains the non-periodic asymmetry rather than only sector
   reconstruction.
6. Before OSS publication, add a small synthetic sample dataset and a tutorial
   config so new users can reproduce the expected outputs without proprietary
   simulation data.
