from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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
    z_db = 20.0 * np.log10(np.maximum(z, np.finfo(float).tiny))

    fig, ax = plt.subplots(figsize=(10, 6))
    mesh = ax.pcolormesh(xf, ym, z_db.T, shading="auto")
    ax.set_xlabel("Temporal frequency [Hz] or rotor order")
    ax.set_ylabel("Spatial order m")
    ax.set_title(title)
    fig.colorbar(mesh, ax=ax, label="|complex coefficient| [dB re 1 Pa]")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def _time_harmonic_order(step_frequency: np.ndarray) -> np.ndarray:
    positive = np.sort(np.abs(step_frequency[np.abs(step_frequency) > 0.0]))
    if len(positive) == 0:
        return step_frequency.copy()
    fundamental = float(positive[0])
    order = step_frequency / fundamental
    rounded = np.round(order)
    return np.where(np.isclose(order, rounded, rtol=1e-6, atol=1e-6), rounded, order)


def _sparse_mode_points(
    step_frequency: np.ndarray,
    spatial_order: np.ndarray,
    coeff: np.ndarray,
    max_order: int | None = None,
    max_step_frequency: float | None = None,
    max_time_order: int | None = None,
    time_order_min: float | None = None,
    time_order_max: float | None = None,
    space_order_min: float | None = None,
    space_order_max: float | None = None,
    positive_space_only: bool = True,
    top_n: int = 80,
    min_relative: float = 0.02,
    include_dc: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mag = np.abs(coeff)
    time_order = _time_harmonic_order(step_frequency)
    if positive_space_only:
        plot_spatial_order = np.unique(np.abs(spatial_order))
        folded_mag = np.zeros((mag.shape[0], len(plot_spatial_order)))
        for i, order in enumerate(plot_spatial_order):
            folded_mag[:, i] = np.max(mag[:, np.isclose(np.abs(spatial_order), order)], axis=1)
        mag = folded_mag
    else:
        plot_spatial_order = spatial_order

    mask = np.ones_like(mag, dtype=bool)
    if max_step_frequency is not None:
        mask &= np.abs(step_frequency)[:, None] <= max_step_frequency
    if max_time_order is not None:
        mask &= np.abs(time_order)[:, None] <= max_time_order
    if time_order_min is not None:
        mask &= time_order[:, None] >= time_order_min
    if time_order_max is not None:
        mask &= time_order[:, None] <= time_order_max
    if max_order is not None:
        mask &= plot_spatial_order[None, :] <= max_order
    if space_order_min is not None:
        mask &= plot_spatial_order[None, :] >= space_order_min
    if space_order_max is not None:
        mask &= plot_spatial_order[None, :] <= space_order_max
    if not include_dc:
        mask &= ~(
            np.isclose(time_order, 0.0)[:, None]
            & np.isclose(plot_spatial_order, 0.0)[None, :]
        )

    candidate = np.where(mask, mag, -np.inf)
    finite = candidate[np.isfinite(candidate)]
    finite = finite[finite > 0.0]
    if len(finite) == 0:
        return np.array([]), np.array([]), np.array([])

    threshold = float(np.max(finite) * min_relative)
    candidate = np.where(candidate >= threshold, candidate, -np.inf)
    flat = np.argsort(candidate.ravel())[::-1]
    flat = flat[np.isfinite(candidate.ravel()[flat])][:top_n]
    rows, cols = np.unravel_index(flat, mag.shape)

    x = time_order[rows]
    y = plot_spatial_order[cols]
    z = mag[rows, cols]
    return x, y, z


def _apply_sparse_axes_limits(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    max_order: int | None = None,
    max_time_order: int | None = None,
    time_order_min: float | None = None,
    time_order_max: float | None = None,
    space_order_min: float | None = None,
    space_order_max: float | None = None,
    positive_space_only: bool = True,
) -> None:
    if len(x):
        x_pad = max(1.0, 0.08 * (float(np.max(x)) - float(np.min(x)) + 1.0))
        y_pad = max(1.0, 0.08 * (float(np.max(y)) - float(np.min(y)) + 1.0))
        ax.set_xlim(float(np.min(x)) - x_pad, float(np.max(x)) + x_pad)
        ax.set_ylim(float(np.min(y)) - y_pad, float(np.max(y)) + y_pad)

    if time_order_min is not None or time_order_max is not None:
        xmin, xmax = ax.get_xlim()
        ax.set_xlim(
            time_order_min if time_order_min is not None else xmin,
            time_order_max if time_order_max is not None else xmax,
        )
    elif max_time_order is not None:
        ax.set_xlim(-max_time_order, max_time_order)

    if space_order_min is not None or space_order_max is not None:
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(
            space_order_min if space_order_min is not None else ymin,
            space_order_max if space_order_max is not None else ymax,
        )
    elif max_order is not None:
        if positive_space_only:
            ax.set_ylim(0.0, max_order)
        else:
            ax.set_ylim(-max_order, max_order)
    elif positive_space_only:
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(max(0.0, ymin), ymax)

    if len(x) and np.allclose(x, np.round(x), rtol=1e-6, atol=1e-6):
        xmin, xmax = ax.get_xlim()
        xticks = np.arange(np.ceil(xmin), np.floor(xmax) + 1.0)
        if len(xticks) <= 45:
            ax.set_xticks(xticks)
    if len(y) and np.allclose(y, np.round(y), rtol=1e-6, atol=1e-6):
        ymin, ymax = ax.get_ylim()
        yticks = np.arange(np.ceil(ymin), np.floor(ymax) + 1.0)
        if len(yticks) <= 55:
            ax.set_yticks(yticks)


def _draw_sparse_points(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    title: str,
    label_top_n: int = 0,
    size_scale: float = 260.0,
) -> plt.Collection | None:
    ax.axhline(0.0, color="0.35", linewidth=0.8)
    ax.axvline(0.0, color="0.35", linewidth=0.8)
    ax.set_xlabel("Time harmonic order n")
    ax.set_ylabel("Space harmonic order m")
    ax.set_title(title)
    ax.grid(True, which="both", color="0.55", linewidth=0.7, alpha=0.75)

    if len(z) == 0:
        ax.text(0.5, 0.5, "No modes above threshold", ha="center", va="center")
        return None

    z_norm = z / np.max(z)
    size = 35.0 + size_scale * np.sqrt(z_norm)
    sc = ax.scatter(
        x,
        y,
        c=z,
        s=size,
        cmap="viridis",
        edgecolors="black",
        linewidths=0.45,
    )

    if label_top_n:
        for idx in range(min(label_top_n, len(x))):
            ax.annotate(
                f"{z[idx]:.3g}",
                (x[idx], y[idx]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=8,
            )
    return sc


def plot_sparse_modes(
    step_frequency: np.ndarray,
    spatial_order: np.ndarray,
    coeff: np.ndarray,
    title: str,
    output: Path,
    colorbar_label: str = "|complex coefficient|",
    max_order: int | None = None,
    max_step_frequency: float | None = None,
    max_time_order: int | None = None,
    time_order_min: float | None = None,
    time_order_max: float | None = None,
    space_order_min: float | None = None,
    space_order_max: float | None = None,
    positive_space_only: bool = True,
    top_n: int = 80,
    min_relative: float = 0.02,
    include_dc: bool = False,
    label_top_n: int = 0,
) -> None:
    x, y, z = _sparse_mode_points(
        step_frequency,
        spatial_order,
        coeff,
        max_order=max_order,
        max_step_frequency=max_step_frequency,
        max_time_order=max_time_order,
        time_order_min=time_order_min,
        time_order_max=time_order_max,
        space_order_min=space_order_min,
        space_order_max=space_order_max,
        positive_space_only=positive_space_only,
        top_n=top_n,
        min_relative=min_relative,
        include_dc=include_dc,
    )

    fig, ax = plt.subplots(figsize=(10, 6))
    sc = _draw_sparse_points(ax, x, y, z, title=title, label_top_n=label_top_n)
    _apply_sparse_axes_limits(
        ax,
        x,
        y,
        max_order=max_order,
        max_time_order=max_time_order,
        time_order_min=time_order_min,
        time_order_max=time_order_max,
        space_order_min=space_order_min,
        space_order_max=space_order_max,
        positive_space_only=positive_space_only,
    )
    if sc is not None:
        fig.colorbar(sc, ax=ax, label=colorbar_label)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_sparse_mode_grid(
    spectra: dict[str, dict[str, object]],
    output: Path,
    max_order: int | None = None,
    max_step_frequency: float | None = None,
    max_time_order: int | None = None,
    time_order_min: float | None = None,
    time_order_max: float | None = None,
    space_order_min: float | None = None,
    space_order_max: float | None = None,
    positive_space_only: bool = True,
    top_n: int = 80,
    min_relative: float = 0.02,
    include_dc: bool = False,
    label_top_n: int = 0,
) -> None:
    names = ["Br", "Btheta", "Fr", "Ftheta"]
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharex=True, sharey=True)

    for ax, name in zip(axes.ravel(), names):
        item = spectra[name]
        x, y, z = _sparse_mode_points(
            item["step_frequency"],
            item["spatial_order"],
            item["coefficient"],
            max_order=max_order,
            max_step_frequency=max_step_frequency,
            max_time_order=max_time_order,
            time_order_min=time_order_min,
            time_order_max=time_order_max,
            space_order_min=space_order_min,
            space_order_max=space_order_max,
            positive_space_only=positive_space_only,
            top_n=top_n,
            min_relative=min_relative,
            include_dc=include_dc,
        )
        sc = _draw_sparse_points(
            ax,
            x,
            y,
            z,
            title=f"{name} 2-D FFT",
            label_top_n=label_top_n,
            size_scale=170.0,
        )
        _apply_sparse_axes_limits(
            ax,
            x,
            y,
            max_order=max_order,
            max_time_order=max_time_order,
            time_order_min=time_order_min,
            time_order_max=time_order_max,
            space_order_min=space_order_min,
            space_order_max=space_order_max,
            positive_space_only=positive_space_only,
        )
        if sc is not None:
            fig.colorbar(sc, ax=ax, label=f"|complex coefficient| [{item['unit']}]", shrink=0.82)

    fig.suptitle("Air-gap flux density and Maxwell stress harmonic modes", fontsize=15)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def plot_torque_comparison(
    comparison: pd.DataFrame,
    output: Path,
    title: str = "Torque comparison",
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, height_ratios=[2.0, 1.0])
    ax_wave, ax_err = axes
    time = comparison["time"].to_numpy(float)
    calculated_raw = comparison["calculated_Torque_Nm"].to_numpy(float)
    calculated = -calculated_raw

    ax_wave.plot(time, calculated, color="black", linewidth=2.0, label="-Ftheta integrated")
    reference_cols = [c for c in comparison.columns if c.startswith("reference_")]
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
    for i, col in enumerate(reference_cols):
        ref = comparison[col].to_numpy(float)
        label = col.replace("reference_", "").replace("_Nm", "")
        if "Stator" in label:
            sign = -1.0
        elif "Rotor" in label:
            sign = 1.0
        else:
            corr = np.corrcoef(calculated, ref)[0, 1] if len(ref) > 1 else 1.0
            sign = -1.0 if np.isfinite(corr) and corr < 0.0 else 1.0
        ref_plot = sign * ref
        plot_label = f"-{label}" if sign < 0.0 else label
        ax_wave.plot(
            time,
            ref_plot,
            linewidth=1.4,
            linestyle="--",
            color=colors[i % len(colors)],
            label=plot_label,
        )
        ax_err.plot(
            time,
            calculated - ref_plot,
            linewidth=1.2,
            color=colors[i % len(colors)],
            label=f"calc - {plot_label}",
        )

    ax_wave.set_ylabel("Torque [Nm]")
    ax_wave.set_title(title)
    ax_wave.grid(True, color="0.7", alpha=0.8)
    ax_wave.legend(loc="best")

    ax_err.axhline(0.0, color="0.35", linewidth=0.8)
    ax_err.set_xlabel("Time [s]")
    ax_err.set_ylabel("Error [Nm]")
    ax_err.grid(True, color="0.7", alpha=0.8)
    ax_err.legend(loc="best")

    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)
