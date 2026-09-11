from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from airgap_io import AirgapData

MU0 = 4.0e-7 * np.pi


def check_uniform(x: np.ndarray, name: str, rtol: float = 1e-6) -> float:
    dx = np.diff(x)
    if len(dx) == 0:
        raise ValueError(f"{name}: at least two samples are required.")
    atol = rtol * max(1.0, abs(float(dx[0])))
    if not np.allclose(dx, dx[0], rtol=rtol, atol=atol):
        raise ValueError(f"{name}: samples must be uniformly spaced for FFT.")
    return float(dx[0])


def maxwell_stress(br: np.ndarray, bt: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
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
    dtheta = check_uniform(theta, "theta")
    fx = radius * stack_length * dtheta * np.sum(
        sigma_r * np.cos(theta)[None, :], axis=1
    )
    fy = radius * stack_length * dtheta * np.sum(
        sigma_r * np.sin(theta)[None, :], axis=1
    )
    ump = np.hypot(fx, fy)
    torque = radius**2 * stack_length * dtheta * np.sum(sigma_t, axis=1)
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
    dstep = check_uniform(step, "step")
    dtheta = check_uniform(theta, "theta")
    coeff = np.fft.fft2(x) / x.size
    step_frequency = np.fft.fftfreq(len(step), d=dstep)
    spatial_order = np.fft.fftfreq(len(theta), d=dtheta / (2.0 * np.pi))
    return (
        np.fft.fftshift(step_frequency),
        np.fft.fftshift(spatial_order),
        np.fft.fftshift(coeff, axes=(0, 1)),
    )


def fft1_spectrum(
    x: np.ndarray,
    step: np.ndarray,
    rotor_speed_rpm: float | None = None,
) -> pd.DataFrame:
    dstep = check_uniform(step, "step")
    freq = np.fft.fftshift(np.fft.fftfreq(len(step), d=dstep))
    coeff = np.fft.fftshift(np.fft.fft(x) / len(x))
    data: dict[str, Any] = {
        "frequency_Hz": freq,
        "complex_real": coeff.real,
        "complex_imag": coeff.imag,
        "coefficient_abs": np.abs(coeff),
    }
    if rotor_speed_rpm:
        data["rotor_order"] = freq / (rotor_speed_rpm / 60.0)
    return pd.DataFrame(data)


def time_harmonic_order_from_frequency(frequency: np.ndarray) -> np.ndarray:
    positive = np.sort(np.abs(frequency[np.abs(frequency) > 0.0]))
    if len(positive) == 0:
        return frequency.copy()
    fundamental = float(positive[0])
    order = frequency / fundamental
    rounded = np.round(order)
    return np.where(np.isclose(order, rounded, rtol=1e-6, atol=1e-6), rounded, order)


def dominant_modes(
    step_frequency: np.ndarray,
    spatial_order: np.ndarray,
    coeff: np.ndarray,
    n: int = 30,
    include_dc: bool = True,
) -> pd.DataFrame:
    mag = np.abs(coeff).copy()
    if not include_dc:
        zero_row = int(np.argmin(np.abs(step_frequency)))
        zero_col = int(np.argmin(np.abs(spatial_order)))
        mag[zero_row, zero_col] = -np.inf
    flat = np.argsort(mag.ravel())[::-1][:n]
    rows, cols = np.unravel_index(flat, mag.shape)
    return pd.DataFrame(
        {
            "step_frequency_Hz_or_order": step_frequency[rows],
            "spatial_order_m": spatial_order[cols],
            "complex_real": coeff[rows, cols].real,
            "complex_imag": coeff[rows, cols].imag,
            "coefficient_abs": np.abs(coeff[rows, cols]),
        }
    )


def _extract_pyemsol_force_nodal_torque(obj: dict[str, Any], path: Path) -> pd.DataFrame | None:
    force_nodal = obj.get("postData", {}).get("forceNodal", {})
    force_data = force_nodal.get("forceNodalData")
    time = obj.get("timeStep", {}).get("time")
    if not isinstance(force_data, list) or time is None:
        return None

    cols: dict[str, Any] = {"time": np.asarray(time, dtype=float)}
    names = {"stator": "reference_Torque_Stator_Nm", "rotor": "reference_Torque_Rotor_Nm"}
    for item in force_data:
        if not isinstance(item, dict):
            continue
        property_num = str(item.get("propertyNum", "")).lower()
        if property_num not in names or "forceMZ" not in item:
            continue
        torque = np.asarray(item["forceMZ"], dtype=float)
        if len(torque) != len(cols["time"]):
            raise ValueError(
                f"{path}: forceMZ length for {property_num!r} does not match timeStep.time."
            )
        cols[names[property_num]] = torque

    if len(cols) == 1:
        return None
    return pd.DataFrame(cols)


def extract_torque_references(path: Path) -> pd.DataFrame:
    with path.open("r", encoding="utf-8") as f:
        obj = json.load(f)

    pyemsol = _extract_pyemsol_force_nodal_torque(obj, path)
    if pyemsol is not None:
        return pyemsol

    time = np.asarray(obj["Time"]["data"], dtype=float)
    cols: dict[str, Any] = {"time": time}

    def walk(node: Any, parts: list[str]) -> None:
        if isinstance(node, dict) and "data" in node:
            name = "_".join(p for p in parts if p) or "Torque"
            cols[f"reference_{name}_Nm"] = np.asarray(node["data"], dtype=float)
            return
        if isinstance(node, dict):
            for key, value in node.items():
                walk(value, [*parts, str(key)])

    walk(obj.get("Torque", {}), ["Torque"])
    return pd.DataFrame(cols)


def compare_torque(
    step: np.ndarray,
    calculated_torque: np.ndarray,
    reference: pd.DataFrame,
) -> pd.DataFrame:
    out = pd.DataFrame({"time": step, "calculated_Torque_Nm": calculated_torque})
    ref_time = reference["time"].to_numpy(float)
    for col in reference.columns:
        if col == "time":
            continue
        ref = reference[col].to_numpy(float)
        if len(ref) == len(step) and np.allclose(ref_time, step, rtol=1e-8, atol=1e-12):
            aligned = ref
        else:
            aligned = np.interp(step, ref_time, ref)
        err = calculated_torque - aligned
        denom = np.maximum(np.abs(aligned), np.finfo(float).eps)
        out[col] = aligned
        out[f"error_vs_{col}_Nm"] = err
        out[f"relative_error_vs_{col}"] = err / denom
    return out


def calculate_airgap_fields(data: AirgapData) -> dict[str, np.ndarray]:
    sigma_r, sigma_t = maxwell_stress(data.br, data.bt)
    return {"Fr": sigma_r, "Ftheta": sigma_t}


def torque_ripple_analysis(
    step: np.ndarray,
    torque: np.ndarray,
    dominant_count: int = 12,
    positive_sign: float = -1.0,
) -> tuple[dict[str, float], pd.DataFrame]:
    torque_positive = positive_sign * torque
    mean = float(np.mean(torque_positive))
    ripple = torque_positive - mean
    peak_to_peak = float(np.max(torque_positive) - np.min(torque_positive))
    rms_ripple = float(np.sqrt(np.mean(ripple**2)))
    denom = max(abs(mean), np.finfo(float).eps)
    summary = {
        "mean_torque_Nm": mean,
        "min_torque_Nm": float(np.min(torque_positive)),
        "max_torque_Nm": float(np.max(torque_positive)),
        "peak_to_peak_ripple_Nm": peak_to_peak,
        "rms_ripple_Nm": rms_ripple,
        "peak_to_peak_ripple_percent_of_mean": 100.0 * peak_to_peak / denom,
        "rms_ripple_percent_of_mean": 100.0 * rms_ripple / denom,
    }

    spectrum = fft1_spectrum(torque_positive, step)
    spectrum["time_harmonic_order_n"] = time_harmonic_order_from_frequency(
        spectrum["frequency_Hz"].to_numpy(float)
    )
    spectrum["single_sided_amplitude_Nm"] = np.where(
        np.isclose(spectrum["frequency_Hz"], 0.0),
        spectrum["coefficient_abs"],
        2.0 * spectrum["coefficient_abs"],
    )
    dominant = spectrum.loc[~np.isclose(spectrum["frequency_Hz"], 0.0)].copy()
    dominant = dominant.loc[dominant["frequency_Hz"] > 0.0]
    dominant = dominant.sort_values("coefficient_abs", ascending=False).head(dominant_count)
    return summary, dominant.reset_index(drop=True)


def mode_pair_contributions(
    br: np.ndarray,
    bt: np.ndarray,
    step: np.ndarray,
    theta: np.ndarray,
    target_count: int = 8,
    pairs_per_target: int = 12,
    include_dc: bool = False,
) -> pd.DataFrame:
    dstep = check_uniform(step, "step")
    dtheta = check_uniform(theta, "theta")
    nt, nth = br.shape
    size = br.size
    time_frequency = np.fft.fftfreq(nt, d=dstep)
    time_order = time_harmonic_order_from_frequency(time_frequency)
    spatial_order = np.fft.fftfreq(nth, d=dtheta / (2.0 * np.pi))

    br_c = np.fft.fft2(br) / size
    bt_c = np.fft.fft2(bt) / size
    sigma_r, sigma_t = maxwell_stress(br, bt)
    stress_coeffs = {
        "Fr": np.fft.fft2(sigma_r) / size,
        "Ftheta": np.fft.fft2(sigma_t) / size,
    }
    products = {
        "Fr": [
            ("Br*Br/(2mu0)", "Br", br_c, "Br", br_c, 1.0 / (2.0 * MU0)),
            ("-Btheta*Btheta/(2mu0)", "Btheta", bt_c, "Btheta", bt_c, -1.0 / (2.0 * MU0)),
        ],
        "Ftheta": [
            ("Br*Btheta/mu0", "Br", br_c, "Btheta", bt_c, 1.0 / MU0),
        ],
    }

    rows: list[dict[str, Any]] = []
    for stress_name, stress_c in stress_coeffs.items():
        mag = np.abs(stress_c).copy()
        if not include_dc:
            mag[0, 0] = -np.inf
        target_flat = np.argsort(mag.ravel())[::-1]
        target_flat = target_flat[np.isfinite(mag.ravel()[target_flat])][:target_count]
        target_indices = np.unravel_index(target_flat, mag.shape)

        for target_t, target_m in zip(*target_indices):
            target_coeff = stress_c[target_t, target_m]
            target_abs = abs(target_coeff)
            contributions: list[dict[str, Any]] = []
            for product_name, a_name, a_coeff, b_name, b_coeff, scale in products[stress_name]:
                for ia in range(nt):
                    ib = (target_t - ia) % nt
                    for ja in range(nth):
                        jb = (target_m - ja) % nth
                        contribution = scale * a_coeff[ia, ja] * b_coeff[ib, jb]
                        contributions.append(
                            {
                                "product": product_name,
                                "source_a": a_name,
                                "time_frequency_a_Hz": time_frequency[ia],
                                "time_order_a_n": time_order[ia],
                                "spatial_order_a_m": spatial_order[ja],
                                "source_b": b_name,
                                "time_frequency_b_Hz": time_frequency[ib],
                                "time_order_b_n": time_order[ib],
                                "spatial_order_b_m": spatial_order[jb],
                                "contribution_real": contribution.real,
                                "contribution_imag": contribution.imag,
                                "contribution_abs": abs(contribution),
                            }
                        )

            contributions.sort(key=lambda row: row["contribution_abs"], reverse=True)
            for rank, row in enumerate(contributions[:pairs_per_target], start=1):
                row.update(
                    {
                        "stress": stress_name,
                        "target_time_frequency_Hz": time_frequency[target_t],
                        "target_time_order_n": time_order[target_t],
                        "target_spatial_order_m": spatial_order[target_m],
                        "target_coeff_real": target_coeff.real,
                        "target_coeff_imag": target_coeff.imag,
                        "target_coeff_abs": target_abs,
                        "pair_rank": rank,
                        "abs_fraction_of_target": row["contribution_abs"]
                        / max(target_abs, np.finfo(float).eps),
                    }
                )
                rows.append(row)

    return pd.DataFrame(rows)
