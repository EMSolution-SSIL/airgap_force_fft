from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MU0 = 4.0e-7 * np.pi


@dataclass
class AirgapData:
    step: np.ndarray
    theta: np.ndarray
    br: np.ndarray
    bt: np.ndarray
    step_name: str
    step_unit: str


def _check_uniform(x: np.ndarray, name: str, rtol: float = 1e-6) -> float:
    dx = np.diff(x)
    if len(dx) == 0:
        raise ValueError(f"{name}: at least two samples are required.")
    atol = rtol * max(1.0, abs(float(dx[0])))
    if not np.allclose(dx, dx[0], rtol=rtol, atol=atol):
        raise ValueError(f"{name}: samples must be uniformly spaced for FFT.")
    return float(dx[0])


def load_airgap_csv(
    path: Path,
    step_col: str = "time",
    theta_col: str = "theta_deg",
    br_col: str = "Br",
    bt_col: str = "Btheta",
) -> AirgapData:
    """
    Expected long-form CSV columns:
        step_col, theta_col, br_col, bt_col

    Example:
        time,theta_deg,Br,Btheta
        0.0000,0.0,...
        0.0000,1.0,...
        ...
        0.0001,0.0,...
    """
    df = pd.read_csv(path)
    required = [step_col, theta_col, br_col, bt_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    # 0 deg and 360 deg are the same spatial point.
    # Grouping makes the prototype tolerant of a duplicated endpoint.
    theta_mod = np.mod(df[theta_col].to_numpy(float), 360.0)
    df = df.assign(_theta_deg=theta_mod)
    df = (
        df.groupby([step_col, "_theta_deg"], as_index=False)[[br_col, bt_col]]
        .mean()
        .sort_values([step_col, "_theta_deg"])
    )

    step = np.sort(df[step_col].unique().astype(float))
    theta_deg = np.sort(df["_theta_deg"].unique().astype(float))

    if len(df) != len(step) * len(theta_deg):
        raise ValueError(
            "Input is not a complete rectangular grid of step x theta samples."
        )

    br = (
        df.pivot(index=step_col, columns="_theta_deg", values=br_col)
        .reindex(index=step, columns=theta_deg)
        .to_numpy()
    )
    bt = (
        df.pivot(index=step_col, columns="_theta_deg", values=bt_col)
        .reindex(index=step, columns=theta_deg)
        .to_numpy()
    )

    if np.isnan(br).any() or np.isnan(bt).any():
        raise ValueError("NaN detected after pivoting.")

    theta = np.deg2rad(theta_deg)

    _check_uniform(step, step_col)
    dtheta = _check_uniform(theta, theta_col)

    # Initial implementation requirement:
    # full 360-degree air-gap sampling, endpoint excluded.
    span = dtheta * len(theta)
    if not np.isclose(span, 2.0 * np.pi, rtol=1e-5, atol=1e-8):
        raise ValueError(
            "Prototype requires one full 360-degree spatial period "
            "(endpoint 360 deg excluded)."
        )

    if step_col == "time":
        step_unit = "Hz"
    elif "angle" in step_col.lower():
        step_unit = "order"
    else:
        step_unit = "cycles_per_step_unit"

    return AirgapData(
        step=step,
        theta=theta,
        br=br,
        bt=bt,
        step_name=step_col,
        step_unit=step_unit,
    )


def maxwell_stress(
    br: np.ndarray, bt: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Air-gap Maxwell stress [Pa].

    sigma_r:
        radial normal stress
        = (Br^2 - Btheta^2) / (2 mu0)

    sigma_t:
        tangential shear stress
        = Br * Btheta / mu0

    Sign conventions depend on the definition of Br, Btheta and surface normal.
    """
    sigma_r = (br**2 - bt**2) / (2.0 * MU0)
    sigma_t = br * bt / MU0
    return sigma_r, sigma_t


def integrated_quantities(
    theta: np.ndarray,
    sigma_r: np.ndarray,
    sigma_t: np.ndarray,
    radius: float,
    stack_length: float,
) -> pd.DataFrame:
    """
    Integrate stress over the cylindrical air-gap surface.

    dA = radius * stack_length * dtheta

    Fx, Fy:
        resultant radial electromagnetic force [N]

    UMP:
        sqrt(Fx^2 + Fy^2) [N]

    Torque:
        radius * integral(sigma_t dA) [N m]
    """
    dtheta = 2.0 * np.pi / len(theta)

    fx = radius * stack_length * dtheta * np.sum(
        sigma_r * np.cos(theta)[None, :], axis=1
    )
    fy = radius * stack_length * dtheta * np.sum(
        sigma_r * np.sin(theta)[None, :], axis=1
    )
    ump = np.hypot(fx, fy)

    torque = (
        radius**2
        * stack_length
        * dtheta
        * np.sum(sigma_t, axis=1)
    )

    return pd.DataFrame(
        {
            "Fx_N": fx,
            "Fy_N": fy,
            "UMP_N": ump,
            "Torque_Nm": torque,
        }
    )


def fft2_spectrum(
    x: np.ndarray,
    step: np.ndarray,
    theta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    2-D DFT of x(step, theta).

    Array shape:
        x.shape == (Nstep, Ntheta)

    Returns:
        step_frequency:
            signed temporal frequency when step is time [Hz]
            or signed order when step is rotor mechanical angle [rev]

        spatial_order:
            signed circumferential spatial order m

        coeff:
            shifted complex Fourier coefficient,
            normalized by Nstep * Ntheta

    NumPy DFT convention:
        cos(m*theta - 2*pi*f*t)
    produces conjugate peaks around
        (m, -f) and (-m, +f).

    Keep signed axes in the first implementation so that traveling-wave
    direction information is not lost.
    """
    dstep = _check_uniform(step, "step")
    dtheta = _check_uniform(theta, "theta")

    coeff = np.fft.fft2(x) / x.size

    step_frequency = np.fft.fftfreq(len(step), d=dstep)

    # theta/(2*pi) is mechanical revolutions.
    # For full 360-deg data, this returns integer spatial order m.
    dtheta_rev = dtheta / (2.0 * np.pi)
    spatial_order = np.fft.fftfreq(len(theta), d=dtheta_rev)

    return (
        np.fft.fftshift(step_frequency),
        np.fft.fftshift(spatial_order),
        np.fft.fftshift(coeff, axes=(0, 1)),
    )


def plot_spectrum(
    step_frequency: np.ndarray,
    spatial_order: np.ndarray,
    coeff: np.ndarray,
    title: str,
    output: Path,
    max_order: int | None = None,
    max_step_frequency: float | None = None,
) -> None:
    mag = np.abs(coeff)

    row_mask = np.ones_like(step_frequency, dtype=bool)
    col_mask = np.ones_like(spatial_order, dtype=bool)

    if max_step_frequency is not None:
        row_mask &= np.abs(step_frequency) <= max_step_frequency
    if max_order is not None:
        col_mask &= np.abs(spatial_order) <= max_order

    z = mag[np.ix_(row_mask, col_mask)]
    xf = step_frequency[row_mask]
    ym = spatial_order[col_mask]

    # Display in dB. Saved spectrum stays as linear complex coefficients.
    z_db = 20.0 * np.log10(np.maximum(z, np.finfo(float).tiny))

    fig, ax = plt.subplots(figsize=(10, 6))
    mesh = ax.pcolormesh(xf, ym, z_db.T, shading="auto")
    ax.set_xlabel("Temporal frequency / rotor order")
    ax.set_ylabel("Spatial order m")
    ax.set_title(title)
    fig.colorbar(mesh, ax=ax, label="|complex coefficient| [dB re 1 Pa]")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def dominant_modes(
    step_frequency: np.ndarray,
    spatial_order: np.ndarray,
    coeff: np.ndarray,
    n: int = 30,
) -> pd.DataFrame:
    """Return the largest complex spectral bins for quick inspection."""
    mag = np.abs(coeff)
    flat = np.argsort(mag.ravel())[::-1][:n]
    rows, cols = np.unravel_index(flat, mag.shape)

    return pd.DataFrame(
        {
            "step_frequency_or_order": step_frequency[rows],
            "spatial_order_m": spatial_order[cols],
            "complex_real": coeff[rows, cols].real,
            "complex_imag": coeff[rows, cols].imag,
            "coefficient_abs": mag[rows, cols],
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="2-D FFT of air-gap Maxwell radial/tangential stress."
    )
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("--step-col", default="time")
    parser.add_argument("--theta-col", default="theta_deg")
    parser.add_argument("--br-col", default="Br")
    parser.add_argument("--bt-col", default="Btheta")
    parser.add_argument(
        "--radius", type=float, required=True, help="Air-gap evaluation radius [m]"
    )
    parser.add_argument(
        "--stack-length", type=float, required=True, help="Stack length [m]"
    )
    parser.add_argument("--output-dir", type=Path, default=Path("fft_output"))
    parser.add_argument("--max-order", type=int, default=30)
    parser.add_argument("--max-step-frequency", type=float, default=None)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    data = load_airgap_csv(
        args.input_csv,
        step_col=args.step_col,
        theta_col=args.theta_col,
        br_col=args.br_col,
        bt_col=args.bt_col,
    )

    sigma_r, sigma_t = maxwell_stress(data.br, data.bt)

    metrics = integrated_quantities(
        data.theta,
        sigma_r,
        sigma_t,
        args.radius,
        args.stack_length,
    )
    metrics.insert(0, data.step_name, data.step)
    metrics.to_csv(args.output_dir / "integrated_metrics.csv", index=False)

    for name, field in [("Fr", sigma_r), ("Ftheta", sigma_t)]:
        sf, sm, coeff = fft2_spectrum(field, data.step, data.theta)

        np.savez_compressed(
            args.output_dir / f"{name}_spectrum.npz",
            step_frequency_or_order=sf,
            spatial_order=sm,
            coefficient=coeff,
        )

        dominant_modes(sf, sm, coeff).to_csv(
            args.output_dir / f"{name}_dominant_modes.csv",
            index=False,
        )

        plot_spectrum(
            sf,
            sm,
            coeff,
            title=f"{name} 2-D FFT",
            output=args.output_dir / f"{name}_2d_fft.png",
            max_order=args.max_order,
            max_step_frequency=args.max_step_frequency,
        )


if __name__ == "__main__":
    main()
