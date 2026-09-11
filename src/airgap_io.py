from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


@dataclass
class AirgapData:
    step: np.ndarray
    theta: np.ndarray
    br: np.ndarray
    bt: np.ndarray
    step_name: str
    step_unit: str
    radius: float | None = None
    stack_length: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


def _deduplicate_endpoint(theta_deg: np.ndarray, *arrays: np.ndarray) -> tuple[np.ndarray, ...]:
    theta_mod = np.mod(theta_deg, 360.0)
    order = np.argsort(theta_mod)
    theta_sorted = theta_mod[order]
    arrays_sorted = [a[..., order] for a in arrays]

    if len(theta_sorted) > 1 and np.isclose(theta_sorted[-1], 360.0):
        theta_sorted = theta_sorted[:-1]
        arrays_sorted = [a[..., :-1] for a in arrays_sorted]

    return (theta_sorted, *arrays_sorted)


def _infer_periodic_span(theta_deg: np.ndarray) -> tuple[float, int]:
    theta_sorted = np.sort(theta_deg)
    dtheta = np.diff(theta_sorted)
    if len(dtheta) == 0:
        raise ValueError("theta: at least two samples are required.")
    nominal_step = float(np.median(dtheta))
    span_deg = float(theta_sorted[-1] - theta_sorted[0] + nominal_step)
    periods = max(1, int(round(360.0 / span_deg)))
    if not np.isclose(span_deg * periods, 360.0, rtol=2e-3, atol=0.05):
        periods = 1
    return 360.0 / periods, periods


def _periodic_interp_rows(
    theta_deg: np.ndarray,
    values: np.ndarray,
    target_theta_deg: np.ndarray,
    period_deg: float,
) -> np.ndarray:
    order = np.argsort(theta_deg)
    x = theta_deg[order]
    y = values[:, order]
    x_aug = np.concatenate([x - period_deg, x, x + period_deg])
    out = np.empty((values.shape[0], len(target_theta_deg)), dtype=float)
    for i in range(values.shape[0]):
        y_aug = np.concatenate([y[i], y[i], y[i]])
        out[i] = np.interp(target_theta_deg, x_aug, y_aug)
    return out


def _resample_to_uniform_sector(
    theta_deg: np.ndarray,
    br: np.ndarray,
    bt: np.ndarray,
    sector_angle_deg: float,
    spatial_samples: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ntheta = int(spatial_samples or len(theta_deg))
    target = (np.arange(ntheta, dtype=float) + 0.5) * sector_angle_deg / ntheta
    br_u = _periodic_interp_rows(theta_deg, br, target, sector_angle_deg)
    bt_u = _periodic_interp_rows(theta_deg, bt, target, sector_angle_deg)
    return target, br_u, bt_u


def expand_sector(
    theta_deg: np.ndarray,
    br: np.ndarray,
    bt: np.ndarray,
    sector_angle_deg: float,
    sector_periods: int,
    symmetricity: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if sector_periods <= 1:
        return theta_deg, br, bt
    if not np.isclose(sector_angle_deg * sector_periods, 360.0, rtol=1e-8, atol=1e-8):
        raise ValueError("sector_angle_deg * sector_periods must be 360 deg.")
    theta_parts = [theta_deg + k * sector_angle_deg for k in range(sector_periods)]
    theta_full = np.concatenate(theta_parts)
    if symmetricity == 0:
        signs = np.ones(sector_periods)
    elif symmetricity == 1:
        if sector_periods % 2 != 0:
            raise ValueError("Anti-periodic symmetricity requires an even number of sectors.")
        signs = np.array([1.0 if k % 2 == 0 else -1.0 for k in range(sector_periods)])
    else:
        raise ValueError(f"Unsupported symmetricity: {symmetricity}")
    br_full = np.concatenate([sign * br for sign in signs], axis=1)
    bt_full = np.concatenate([sign * bt for sign in signs], axis=1)
    return theta_full, br_full, bt_full


def _resolve_sector(
    theta_deg: np.ndarray,
    sector_angle_deg: float | None,
    sector_periods: int | str,
) -> tuple[float, int]:
    inferred_angle, inferred_periods = _infer_periodic_span(theta_deg)
    if sector_periods == "auto":
        periods = inferred_periods
    else:
        periods = int(sector_periods)
    if sector_angle_deg is None:
        sector_angle = 360.0 / periods if periods > 1 else inferred_angle
    else:
        sector_angle = float(sector_angle_deg)
    return sector_angle, periods


def load_airgap_csv(
    path: Path,
    step_col: str = "time",
    theta_col: str = "theta_deg",
    br_col: str = "Br",
    bt_col: str = "Btheta",
    sector_periods: int | str = "auto",
    sector_angle_deg: float | None = None,
    symmetricity: int = 0,
    resample_spatial: bool = False,
    spatial_samples: int | None = None,
) -> AirgapData:
    df = pd.read_csv(path)
    required = [step_col, theta_col, br_col, bt_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

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
        raise ValueError("Input is not a complete rectangular grid of step x theta samples.")

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

    sector_angle, periods = _resolve_sector(theta_deg, sector_angle_deg, sector_periods)
    if resample_spatial:
        theta_deg, br, bt = _resample_to_uniform_sector(
            theta_deg, br, bt, sector_angle, spatial_samples
        )
    theta_deg, br, bt = expand_sector(theta_deg, br, bt, sector_angle, periods, symmetricity)

    step_unit = "Hz" if step_col == "time" else "cycles_per_step_unit"
    return AirgapData(
        step=step,
        theta=np.deg2rad(theta_deg),
        br=br,
        bt=bt,
        step_name=step_col,
        step_unit=step_unit,
        metadata={
            "input_format": "csv",
            "sector_angle_deg": sector_angle,
            "sector_periods": periods,
            "symmetricity": symmetricity,
        },
    )


def load_airgap_json(
    path: Path,
    sector_periods: int | str = "auto",
    sector_angle_deg: float | None = None,
    resample_spatial: bool = True,
    spatial_samples: int | None = None,
) -> AirgapData:
    with path.open("r", encoding="utf-8") as f:
        obj = json.load(f)

    try:
        gap = obj["postData"]["gapB"]
        position = gap["position"]
        magnetic_density = gap["magneticDensity"]
    except KeyError as exc:
        raise ValueError("JSON does not look like EMSolution gapB post data.") from exc

    step = np.asarray(obj["timeStep"]["time"], dtype=float)
    theta_deg = np.asarray(position["theta"], dtype=float)
    radius = float(np.mean(np.asarray(position.get("r", []), dtype=float)))
    z = np.asarray(position.get("z", []), dtype=float)
    stack_length = float(2.0 * np.mean(z)) if len(z) else None

    br = np.asarray([row["br"] for row in magnetic_density], dtype=float)
    bt = np.asarray([row["btheta"] for row in magnetic_density], dtype=float)
    if br.shape != bt.shape:
        raise ValueError("br and btheta arrays have different shapes.")
    if br.shape != (len(step), len(theta_deg)):
        raise ValueError("magneticDensity shape does not match time/theta axes.")

    theta_deg, br, bt = _deduplicate_endpoint(theta_deg, br, bt)
    sector_angle, periods = _resolve_sector(theta_deg, sector_angle_deg, sector_periods)
    if resample_spatial:
        theta_deg, br, bt = _resample_to_uniform_sector(
            theta_deg, br, bt, sector_angle, spatial_samples
        )
    symmetricity = int(gap.get("symmetricity", 0))
    theta_deg, br, bt = expand_sector(
        theta_deg,
        br,
        bt,
        sector_angle,
        periods,
        symmetricity=symmetricity,
    )

    metadata = {
        "input_format": "emsolution_gapb_json",
        "time_unit": obj.get("timeStep", {}).get("timeUnit"),
        "magnetic_density_unit": gap.get("magneticDensityUnit"),
        "sector_angle_deg": sector_angle,
        "sector_periods": periods,
        "symmetricity": symmetricity,
        "original_theta_count": int(len(position["theta"])),
        "theta_count": int(len(theta_deg)),
        "radius_from_json_m": radius,
    }
    if len(z):
        metadata["z_min_m"] = float(np.min(z))
        metadata["z_max_m"] = float(np.max(z))
        metadata["stack_length_from_2z_m"] = stack_length

    return AirgapData(
        step=step,
        theta=np.deg2rad(theta_deg),
        br=br,
        bt=bt,
        step_name="time",
        step_unit="Hz",
        radius=radius,
        stack_length=stack_length,
        metadata=metadata,
    )


def load_airgap(path: Path, input_format: str = "auto", **kwargs: Any) -> AirgapData:
    fmt = input_format.lower()
    if fmt == "auto":
        fmt = "json" if path.suffix.lower() == ".json" else "csv"
    if fmt == "csv":
        return load_airgap_csv(path, **kwargs)
    if fmt in {"json", "emsolution-json", "emsolution_gapb_json"}:
        json_kwargs = {
            k: v
            for k, v in kwargs.items()
            if k in {"sector_periods", "sector_angle_deg", "resample_spatial", "spatial_samples"}
        }
        return load_airgap_json(path, **json_kwargs)
    raise ValueError(f"Unsupported input format: {input_format}")
