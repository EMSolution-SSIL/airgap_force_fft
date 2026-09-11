from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from airgap_fft import (
    calculate_airgap_fields,
    compare_torque,
    dominant_modes,
    extract_torque_references,
    fft1_spectrum,
    fft2_spectrum,
    integrated_quantities,
    mode_pair_contributions,
    torque_ripple_analysis,
)
from airgap_io import AirgapData, load_airgap
from airgap_plotting import plot_sparse_mode_grid, plot_sparse_modes, plot_spectrum, plot_torque_comparison


DEFAULT_CONFIG: dict[str, Any] = {
    "input": {
        "path": "data/GL80/motorGapB.json",
        "format": "auto",
        "csv": {
            "step_col": "time",
            "theta_col": "theta_deg",
            "br_col": "Br",
            "bt_col": "Btheta",
        },
        "sector": {
            "periods": "auto",
            "angle_deg": None,
            "symmetricity": None,
            "resample_spatial": True,
            "spatial_samples": None,
        },
    },
    "geometry": {
        "radius_m": None,
        "stack_length_m": None,
    },
    "analysis": {
        "max_order": 120,
        "max_step_frequency": None,
        "dominant_count": 40,
        "rotor_speed_rpm": None,
        "include_dc": False,
        "exclude_time_endpoint_for_fft": True,
        "torque_reference_json": "data/GL80/transient_results.json",
        "torque_ripple": {
            "dominant_count": 12,
            "positive_sign": -1.0,
        },
        "mode_pair_contributions": {
            "enabled": True,
            "target_count": 8,
            "pairs_per_target": 12,
        },
    },
    "plotting": {
        "sparse_modes": {
            "enabled": True,
            "top_n": 80,
            "min_relative": 0.02,
            "max_time_order": None,
            "time_order_min": None,
            "time_order_max": None,
            "space_order_min": 0,
            "space_order_max": None,
            "positive_space_only": True,
            "label_top_n": 0,
        },
        "heatmap": {
            "enabled": False,
        },
        "mode_grid": {
            "enabled": True,
            "filename": "airgap_2d_fft_4panel.png",
        },
        "torque_comparison": {
            "enabled": True,
            "filename": "torque_comparison.png",
        },
    },
    "output": {
        "dir": "data/GL80/fft_output",
    },
}


def deep_update(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy(base)
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_update(out[key], value)
        else:
            out[key] = value
    return out


def load_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        user_config = yaml.safe_load(f) or {}
    if not isinstance(user_config, dict):
        raise ValueError(f"{path} must contain a YAML mapping.")
    return deep_update(DEFAULT_CONFIG, user_config)


def resolve_path(value: str | Path | None, base_dir: Path) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value)
    return path if path.is_absolute() else base_dir / path


def resolve_geometry(data: AirgapData, config: dict[str, Any]) -> tuple[float, float]:
    geometry = config["geometry"]
    radius = geometry.get("radius_m") if geometry.get("radius_m") is not None else data.radius
    if radius is None:
        raise ValueError("geometry.radius_m is required when the input file has no radius data.")

    stack_length = (
        geometry.get("stack_length_m")
        if geometry.get("stack_length_m") is not None
        else data.stack_length
    )
    if stack_length is None:
        raise ValueError(
            "geometry.stack_length_m is required when the input file has no position.z data."
        )
    return float(radius), float(stack_length)


def write_analysis(
    data: AirgapData,
    config: dict[str, Any],
    config_path: Path,
) -> None:
    base_dir = config_path.parent
    output_dir = resolve_path(config["output"]["dir"], base_dir)
    if output_dir is None:
        raise ValueError("output.dir is required.")
    output_dir.mkdir(parents=True, exist_ok=True)

    radius, stack_length = resolve_geometry(data, config)
    analysis = config["analysis"]
    plotting = config["plotting"]
    sparse_cfg = plotting["sparse_modes"]
    heatmap_cfg = plotting["heatmap"]

    fields = calculate_airgap_fields(data)
    sigma_r = fields["Fr"]
    sigma_t = fields["Ftheta"]
    metrics = integrated_quantities(data.theta, sigma_r, sigma_t, radius, stack_length)
    metrics.insert(0, data.step_name, data.step)
    metrics.to_csv(output_dir / "integrated_metrics.csv", index=False)

    fft_slice = slice(None, -1) if analysis.get("exclude_time_endpoint_for_fft", False) else slice(None)
    fft_step = data.step[fft_slice]
    fft_br = data.br[fft_slice, :]
    fft_bt = data.bt[fft_slice, :]
    fft_fields = {name: field[fft_slice, :] for name, field in fields.items()}
    fft_metrics = metrics.iloc[fft_slice].reset_index(drop=True)

    fft1_spectrum(
        fft_metrics["Torque_Nm"].to_numpy(float),
        fft_step,
        analysis.get("rotor_speed_rpm"),
    ).to_csv(output_dir / "Torque_temporal_spectrum.csv", index=False)

    ripple_cfg = analysis.get("torque_ripple", {})
    ripple_summary, ripple_orders = torque_ripple_analysis(
        fft_step,
        fft_metrics["Torque_Nm"].to_numpy(float),
        dominant_count=int(ripple_cfg.get("dominant_count", 12)),
        positive_sign=float(ripple_cfg.get("positive_sign", -1.0)),
    )
    with (output_dir / "torque_ripple_summary.json").open("w", encoding="utf-8") as f:
        json.dump(ripple_summary, f, indent=2, ensure_ascii=False)
    ripple_orders.to_csv(output_dir / "torque_ripple_orders.csv", index=False)

    spectral_fields = {
        "Br": {"values": fft_br, "unit": "T"},
        "Btheta": {"values": fft_bt, "unit": "T"},
        "Fr": {"values": fft_fields["Fr"], "unit": "Pa"},
        "Ftheta": {"values": fft_fields["Ftheta"], "unit": "Pa"},
    }
    spectra: dict[str, dict[str, object]] = {}

    for name, item in spectral_fields.items():
        field = item["values"]
        unit = item["unit"]
        sf, sm, coeff = fft2_spectrum(field, fft_step, data.theta)
        spectra[name] = {
            "step_frequency": sf,
            "spatial_order": sm,
            "coefficient": coeff,
            "unit": unit,
        }
        np.savez_compressed(
            output_dir / f"{name}_spectrum.npz",
            step_frequency_Hz_or_order=sf,
            spatial_order=sm,
            coefficient=coeff,
        )
        dominant_modes(
            sf,
            sm,
            coeff,
            n=int(analysis["dominant_count"]),
            include_dc=bool(analysis["include_dc"]),
        ).to_csv(output_dir / f"{name}_dominant_modes.csv", index=False)
        if sparse_cfg.get("enabled", True):
            plot_sparse_modes(
                sf,
                sm,
                coeff,
                title=f"{name} 2-D FFT sparse modes",
                output=output_dir / f"{name}_2d_fft.png",
                colorbar_label=f"|complex coefficient| [{unit}]",
                max_order=analysis.get("max_order"),
                max_step_frequency=analysis.get("max_step_frequency"),
                max_time_order=sparse_cfg.get("max_time_order"),
                time_order_min=sparse_cfg.get("time_order_min"),
                time_order_max=sparse_cfg.get("time_order_max"),
                space_order_min=sparse_cfg.get("space_order_min"),
                space_order_max=sparse_cfg.get("space_order_max"),
                positive_space_only=bool(sparse_cfg.get("positive_space_only", True)),
                top_n=int(sparse_cfg.get("top_n", 80)),
                min_relative=float(sparse_cfg.get("min_relative", 0.02)),
                include_dc=bool(analysis["include_dc"]),
                label_top_n=int(sparse_cfg.get("label_top_n", 0)),
            )
        if heatmap_cfg.get("enabled", False):
            plot_spectrum(
                sf,
                sm,
                coeff,
                title=f"{name} 2-D FFT heatmap",
                output=output_dir / f"{name}_2d_fft_heatmap.png",
                max_order=analysis.get("max_order"),
                max_step_frequency=analysis.get("max_step_frequency"),
            )

        spatial_m0 = np.mean(field, axis=1)
        fft1_spectrum(spatial_m0, fft_step, analysis.get("rotor_speed_rpm")).to_csv(
            output_dir / f"{name}_m0_temporal_spectrum.csv",
            index=False,
        )

    pair_cfg = analysis.get("mode_pair_contributions", {})
    if pair_cfg.get("enabled", True):
        mode_pair_contributions(
            fft_br,
            fft_bt,
            fft_step,
            data.theta,
            target_count=int(pair_cfg.get("target_count", 8)),
            pairs_per_target=int(pair_cfg.get("pairs_per_target", 12)),
            include_dc=bool(analysis["include_dc"]),
        ).to_csv(output_dir / "mode_pair_contributions.csv", index=False)

    mode_grid_cfg = plotting.get("mode_grid", {})
    if mode_grid_cfg.get("enabled", True):
        plot_sparse_mode_grid(
            spectra,
            output=output_dir / mode_grid_cfg.get("filename", "airgap_2d_fft_4panel.png"),
            max_order=analysis.get("max_order"),
            max_step_frequency=analysis.get("max_step_frequency"),
            max_time_order=sparse_cfg.get("max_time_order"),
            time_order_min=sparse_cfg.get("time_order_min"),
            time_order_max=sparse_cfg.get("time_order_max"),
            space_order_min=sparse_cfg.get("space_order_min"),
            space_order_max=sparse_cfg.get("space_order_max"),
            positive_space_only=bool(sparse_cfg.get("positive_space_only", True)),
            top_n=int(sparse_cfg.get("top_n", 80)),
            min_relative=float(sparse_cfg.get("min_relative", 0.02)),
            include_dc=bool(analysis["include_dc"]),
            label_top_n=int(sparse_cfg.get("label_top_n", 0)),
        )

    torque_reference = resolve_path(analysis.get("torque_reference_json"), base_dir)
    if torque_reference:
        reference = extract_torque_references(torque_reference)
        comparison = compare_torque(data.step, metrics["Torque_Nm"].to_numpy(float), reference)
        comparison.to_csv(output_dir / "torque_reference_comparison.csv", index=False)
        torque_plot_cfg = plotting.get("torque_comparison", {})
        if torque_plot_cfg.get("enabled", True):
            plot_torque_comparison(
                comparison,
                output=output_dir / torque_plot_cfg.get("filename", "torque_comparison.png"),
            )

    summary = {
        "config_path": str(config_path),
        "effective_config": config,
        "step_count": int(len(data.step)),
        "fft_step_count": int(len(fft_step)),
        "exclude_time_endpoint_for_fft": bool(analysis.get("exclude_time_endpoint_for_fft", False)),
        "theta_count": int(len(data.theta)),
        "radius_m": radius,
        "stack_length_m": stack_length,
        "step_name": data.step_name,
        "step_unit": data.step_unit,
        "metadata": data.metadata,
        "outputs": sorted(p.name for p in output_dir.iterdir()),
    }
    with (output_dir / "analysis_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)


def run_from_config(config_path: Path) -> None:
    config = load_config(config_path)
    base_dir = config_path.parent
    input_cfg = config["input"]
    csv_cfg = input_cfg.get("csv", {})
    sector_cfg = input_cfg.get("sector", {})
    input_path = resolve_path(input_cfg["path"], base_dir)
    if input_path is None:
        raise ValueError("input.path is required.")

    data = load_airgap(
        input_path,
        input_format=input_cfg.get("format", "auto"),
        step_col=csv_cfg.get("step_col", "time"),
        theta_col=csv_cfg.get("theta_col", "theta_deg"),
        br_col=csv_cfg.get("br_col", "Br"),
        bt_col=csv_cfg.get("bt_col", "Btheta"),
        sector_periods=sector_cfg.get("periods", "auto"),
        sector_angle_deg=sector_cfg.get("angle_deg"),
        symmetricity=sector_cfg.get("symmetricity") or 0,
        resample_spatial=bool(sector_cfg.get("resample_spatial", True)),
        spatial_samples=sector_cfg.get("spatial_samples"),
    )
    write_analysis(data, config, config_path)


def write_default_config(path: Path) -> None:
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(DEFAULT_CONFIG, f, sort_keys=False, allow_unicode=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="2-D FFT postprocessor for air-gap Maxwell force from Br/Btheta."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config.yaml"),
        help="YAML file that records input, geometry, analysis, and output settings.",
    )
    parser.add_argument(
        "--write-default-config",
        action="store_true",
        help="Write a template config YAML and exit.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.write_default_config:
        write_default_config(args.config)
        return
    run_from_config(args.config)


if __name__ == "__main__":
    main()
