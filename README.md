# Air-Gap Force FFT

Python postprocessor for motor air-gap flux-density results. It reads `Br` and `Btheta`, calculates Maxwell stress, evaluates integrated force and torque, and extracts time-space harmonic modes with 2-D FFT.

The included sample is the healthy-position result for the GL80 PMSM.

The project is open source under the [MIT License](LICENSE).

## Repository Layout

```text
.
|-- src/                         # Processing scripts and reusable modules
|   |-- airgap_force_fft.py       # CLI entry point
|   |-- airgap_io.py              # JSON/CSV readers and sector expansion
|   |-- airgap_fft.py             # Maxwell stress, FFT, torque comparison
|   `-- airgap_plotting.py        # Plotting utilities
|-- tests/                       # Unit tests
|-- data/
|   `-- GL80/healthy_position/   # Public sample input and generated outputs
|-- docs/                        # Discussion notes and analysis reports
|-- config.yaml                  # Default GL80 example configuration
|-- pyproject.toml               # Packaging metadata
`-- requirements.txt             # Runtime dependencies
```

## Quick Start

Install the published package from PyPI:

```powershell
python -m pip install airgap-force-fft
```

For development, create an environment and install the repository dependencies:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

Run the default GL80 example:

```powershell
python .\src\airgap_force_fft.py --config .\config.yaml
```

The main outputs are written under `data/GL80/healthy_position/fft_output` by
default.

## Input Formats

Supported air-gap flux-density inputs:

- EMSolution/pyemsol `gapB` JSON
- Long-form CSV with columns such as `time`, `theta_deg`, `Br`, and `Btheta`

For EMSolution JSON, geometry and symmetry metadata can be read from the file when available. For CSV input, set the missing metadata explicitly in `config.yaml`, especially:

- `input.sector.periods`
- `input.sector.symmetricity`
- `geometry.radius_m`
- `geometry.stack_length_m`

For full 360-degree CSV data, use `input.sector.periods: 1` and `input.sector.symmetricity: 0`.

## Torque Reference Files

Torque comparison supports two reference formats:

- The earlier `transient_results.json` style with `Time.data` and `Torque.*.data`
- pyemsol `output_transient.json`, where torque is read from `postData.forceNodal.forceNodalData` entries whose `propertyNum` is `stator` or `rotor`, using `forceMZ`

The comparison plot uses rotor-positive convention: `-Ftheta integrated`, `-stator`, and `rotor`.

## Reports

Public documentation is in `docs/`:

- [Air-gap force 2D FFT guide](docs/README_airgap_force_fft.md)
- [Two-dimensional FFT analysis method](docs/airgap_fft_2d_method.md)
- [Reference prototype](docs/airgap_force_fft_prototype.py)

## Tests

Run the unit tests with:

```powershell
python -m unittest discover -s tests -v
```

## Packaging

For editable development installation:

```powershell
pip install -e .
```

After installation, the CLI can also be invoked as:

```powershell
airgap-force-fft --config config.yaml
```

