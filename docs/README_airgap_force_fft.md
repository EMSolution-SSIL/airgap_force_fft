# Air-gap force 2D FFT

`src/airgap_force_fft.py` reads air-gap flux density and evaluates Maxwell stress,
integrated torque/UMP, and 2D FFT modes.

## Inputs

Supported formats:

- EMSolution `gapB` JSON, such as `data/data/GL80/motorGapB.json`
- long-form CSV with default columns:

```csv
time,theta_deg,Br,Btheta
0.0,0.0,0.7,0.02
0.0,1.0,0.71,0.01
```

For EMSolution JSON, the script reads:

- `timeStep.time`
- `postData.gapB.symmetricity`
- `postData.gapB.position.theta`
- `postData.gapB.position.r`
- `postData.gapB.magneticDensity[*].br`
- `postData.gapB.magneticDensity[*].btheta`

If the spatial model is a periodic sector, the script can reconstruct the full
360 degrees. The default `--sector-periods auto` infers this from the angular
span. The GL80 sample is inferred as a 60-degree sector repeated 6 times.

`symmetricity` controls how the sector is reconstructed:

- `0`: periodic copy
- `1`: anti-periodic copy, where adjacent sectors use `-Br` and `-Btheta`

## CSV Metadata

CSV input contains only sampled values, so information that EMSolution JSON can
carry as metadata must be supplied in `config.yaml`.

Required CSV columns by default:

| Column | Meaning |
|---|---|
| `time` | time [s] or another uniformly sampled step axis |
| `theta_deg` | spatial angle [deg] |
| `Br` | radial flux density [T] |
| `Btheta` | circumferential flux density [T] |

For CSV input, choose one of the following workflows:

1. **Full 360-degree CSV**
   - Prepare the CSV as a full 360-degree dataset.
   - Set `input.sector.periods: 1`.
   - Set `input.sector.symmetricity: 0`.
   - Specify `geometry.radius_m` and `geometry.stack_length_m`.

2. **Sector CSV with metadata in config**
   - Prepare the CSV as one periodic or anti-periodic sector.
   - Set `input.sector.periods`, for example `6` for a 60-degree model.
   - Set `input.sector.symmetricity`:
     - `0` for periodic symmetry
     - `1` for anti-periodic symmetry
   - Specify `geometry.radius_m` and `geometry.stack_length_m`.

If the symmetry is uncertain, prefer creating a full 360-degree CSV outside this
tool and using `input.sector.periods: 1`. This avoids accidental use of periodic
copying for anti-periodic data.

A sample CSV is available at `../data/GL80/motorGapB_sample.csv`. It is the original
60-degree GL80 sector converted from JSON, so use `../config_csv.yaml` or set
`input.sector.periods: 6`, `input.sector.symmetricity: 1`, `geometry.radius_m`,
and `geometry.stack_length_m` explicitly.

## Example

```powershell
python .\src\airgap_force_fft.py --config .\config.yaml
```

The run conditions are recorded in `config.yaml`, and the merged effective
configuration is also written to `analysis_summary.json`.

For EMSolution `gapB` JSON, the script uses the mean evaluation radius from
`position.r`. It also uses `2 * mean(position.z)` as the stack length when
`geometry.stack_length_m` is omitted. For the GL80 sample, `z = 0.005 m`, so the
inferred stack length is `0.01 m`.

To create a template config:

```powershell
python .\src\airgap_force_fft.py --write-default-config --config .\config.yaml
```

## Outputs

- `integrated_metrics.csv`: `Fx`, `Fy`, `UMP`, and integrated torque
- `Torque_temporal_spectrum.csv`: 1D FFT of integrated torque
- `torque_ripple_summary.json`: mean torque, peak-to-peak ripple, and RMS ripple
  using the configured positive torque sign
- `torque_ripple_orders.csv`: dominant positive-frequency torque ripple orders
- `Br_spectrum.npz`, `Btheta_spectrum.npz`, `Fr_spectrum.npz`,
  `Ftheta_spectrum.npz`: complex 2D FFT coefficients
- `Br_dominant_modes.csv`, `Btheta_dominant_modes.csv`,
  `Fr_dominant_modes.csv`, `Ftheta_dominant_modes.csv`: largest spectral bins
- `Br_2d_fft.png`, `Btheta_2d_fft.png`, `Fr_2d_fft.png`,
  `Ftheta_2d_fft.png`: sparse mode maps with time harmonic order on x, space
  harmonic order on y, and coefficient magnitude as marker size/color
- `airgap_2d_fft_4panel.png`: 2x2 sparse mode map for `Br`, `Btheta`, `Fr`,
  and `Ftheta` with the same displayed harmonic range
- `torque_comparison.png`: calculated torque and FEM reference torque, plus
  error traces
- `Br_m0_temporal_spectrum.csv`, `Btheta_m0_temporal_spectrum.csv`,
  `Fr_m0_temporal_spectrum.csv`, `Ftheta_m0_temporal_spectrum.csv`: temporal
  FFT of the spatial 0th component
- `torque_reference_comparison.csv`: comparison against transient torque JSON
- `mode_pair_contributions.csv`: leading `Br/Btheta` Fourier coefficient pairs
  contributing to dominant `Fr/Ftheta` modes through Maxwell stress products
- `analysis_summary.json`: run metadata

Torque ripple and mode-pair analysis are controlled by:

```yaml
analysis:
  exclude_time_endpoint_for_fft: true
  torque_ripple:
    dominant_count: 12
    positive_sign: -1.0
  mode_pair_contributions:
    enabled: true
    target_count: 8
    pairs_per_target: 12
```

`torque_ripple.positive_sign: -1.0` reports torque in the usual rotor-positive
direction for the GL80 sample, because the integrated `Ftheta` torque is opposite
in sign. `mode_pair_contributions.csv` is intended for tracing statements such
as `Br(m=18)` and `Btheta(m=24)` combining into `Ftheta(m=42)`.

`exclude_time_endpoint_for_fft: true` is useful when the time series contains
both the start point and the periodic endpoint. Torque comparison still uses all
samples, but FFT outputs exclude the final endpoint to avoid double-counting it.
The current torque ripple order table reports DFT-bin order over the sampled
window, not rotor mechanical order or electrical order.

The sparse mode maps are controlled by:

```yaml
plotting:
  sparse_modes:
    enabled: true
    top_n: 80
    min_relative: 0.02
    max_time_order:
    time_order_min:
    time_order_max:
    space_order_min: 0
    space_order_max:
    positive_space_only: true
    label_top_n: 0
  heatmap:
    enabled: false
  mode_grid:
    enabled: true
    filename: airgap_2d_fft_4panel.png
  torque_comparison:
    enabled: true
    filename: torque_comparison.png
```

Set `heatmap.enabled: true` to also write `Fr_2d_fft_heatmap.png` and
`Ftheta_2d_fft_heatmap.png`.

`positive_space_only: true` folds `+m` and `-m` into the positive space-order
side for readability. Use `time_order_min`, `time_order_max`,
`space_order_min`, and `space_order_max` to limit the displayed harmonic range.

## Notes

If `torque_reference_comparison.csv` differs strongly from the FEM torque,
check the coordinate conversion into `Br` and `Btheta`, the stack length, the
evaluation radius, and the sign convention of `Btheta`.
