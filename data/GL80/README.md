# GL80 Sample Data

This directory contains the healthy-position GL80 sample input and generated
analysis outputs under `healthy_position`.

The motor model is based on the IEEJ public-use GL80 outer-rotor PMSM:
<http://www2.iee.or.jp/~drm/PublicUseMotors/actualMotors/PMSM/GL80/detail.html>.
The JSON and derived CSV/figure files are EMSolution analysis and
postprocessing results prepared for reproducible examples in this repository.

Run from the repository root with `config.yaml`. The sample contains the
full 360-degree air-gap field, its pyemsol torque reference, and selected FFT
CSV/figure outputs. Large reproducible spectrum arrays (`*.npz`) are omitted.
