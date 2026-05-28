#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
from pathlib import Path

import numpy as np

_MPL_CACHE = Path(tempfile.gettempdir()) / "simple-matplotlib"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))
os.environ.setdefault("XDG_CACHE_HOME", str(_MPL_CACHE))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot SIMPLE fsc_area_score CSV outputs in a cryoSPARC-like summary style."
    )
    parser.add_argument(
        "fbody",
        nargs="?",
        default="fsc_area_score",
        help="Output file body used by fsc_area_score, default: fsc_area_score",
    )
    parser.add_argument("--curves", type=Path, help="Path to *_curves.csv")
    parser.add_argument("--directions", type=Path, help="Path to *_directions.csv")
    parser.add_argument("--summary", type=Path, help="Path to *_summary.txt")
    parser.add_argument("--output", type=Path, help="Output image path, default: <fbody>_plot.png")
    parser.add_argument("--threshold", type=float, help="FSC crossing threshold; defaults to summary value or 0.143")
    parser.add_argument("--title", default="Directional FSC area score", help="Plot title")
    parser.add_argument("--dpi", type=int, default=180, help="Raster output DPI")
    return parser.parse_args()


def paths_from_args(args: argparse.Namespace) -> tuple[Path, Path, Path, Path]:
    fbody = Path(args.fbody)
    curves = args.curves or fbody.with_name(fbody.name + "_curves.csv")
    directions = args.directions or fbody.with_name(fbody.name + "_directions.csv")
    summary = args.summary or fbody.with_name(fbody.name + "_summary.txt")
    output = args.output or fbody.with_name(fbody.name + "_plot.png")
    return curves, directions, summary, output


def read_threshold(summary_path: Path, override: float | None) -> float:
    if override is not None:
        return override
    if summary_path.exists():
        with summary_path.open() as handle:
            for line in handle:
                fields = line.split()
                if len(fields) >= 2 and fields[0] == "threshold":
                    return float(fields[1])
    return 0.143


def read_curves(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if len(header) < 3:
            raise ValueError(f"{path} does not contain direction columns")

        wave_numbers: list[float] = []
        resolutions: list[float] = []
        rows: list[list[float]] = []
        for row in reader:
            if not row:
                continue
            wave_numbers.append(float(row[0]))
            resolutions.append(float(row[1]))
            vals = []
            for field in row[2:]:
                field = field.strip()
                vals.append(float(field) if field and field.lower() != "nan" else math.nan)
            rows.append(vals)

    curves = np.asarray(rows, dtype=float)
    return np.asarray(wave_numbers, dtype=float), np.asarray(resolutions, dtype=float), curves


def read_directions(path: Path) -> tuple[np.ndarray, np.ndarray]:
    crossing_shells: list[float] = []
    crossing_resolutions: list[float] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            shell = float(row["crossing_wave_number"])
            resolution = float(row["crossing_resolution_A"])
            if shell > 0.0 and resolution > 0.0 and math.isfinite(resolution):
                crossing_shells.append(shell)
                crossing_resolutions.append(resolution)
    return np.asarray(crossing_shells, dtype=float), np.asarray(crossing_resolutions, dtype=float)


def curve_stats(curves: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    valid = np.isfinite(curves)
    counts = valid.sum(axis=1)
    values = np.where(valid, curves, 0.0)

    mean = np.full(curves.shape[0], np.nan)
    std = np.full(curves.shape[0], np.nan)
    minv = np.full(curves.shape[0], np.nan)
    maxv = np.full(curves.shape[0], np.nan)

    usable = counts > 0
    mean[usable] = values[usable].sum(axis=1) / counts[usable]
    centered = np.where(valid, curves - mean[:, None], 0.0)
    std[usable] = np.sqrt((centered[usable] ** 2).sum(axis=1) / counts[usable])

    for idx in np.flatnonzero(usable):
        row = curves[idx, valid[idx]]
        minv[idx] = row.min()
        maxv[idx] = row.max()

    return mean, std, minv, maxv


def frequency_axis(resolutions: np.ndarray) -> np.ndarray:
    freq = np.zeros_like(resolutions, dtype=float)
    valid = np.isfinite(resolutions) & (resolutions > 0.0)
    freq[valid] = 1.0 / resolutions[valid]
    return freq


def shell_widths(freq: np.ndarray) -> np.ndarray:
    if len(freq) == 1:
        return np.ones(1)
    edges = np.empty(len(freq) + 1, dtype=float)
    edges[1:-1] = 0.5 * (freq[:-1] + freq[1:])
    edges[0] = max(0.0, freq[0] - (edges[1] - freq[0]))
    edges[-1] = freq[-1] + (freq[-1] - edges[-2])
    return np.maximum(np.diff(edges), 1.0e-6)


def resolution_label(freq: float) -> str:
    if freq <= 0.0:
        return "DC"
    resolution = 1.0 / freq
    if resolution >= 10.0:
        return f"{resolution:.0f}A"
    return f"{resolution:.1f}A"


def add_resolution_ticks(ax: plt.Axes, freq: np.ndarray) -> None:
    max_freq = float(np.nanmax(freq))
    ticks = np.linspace(0.0, max_freq, 7)
    ax.set_xticks(ticks)
    ax.set_xticklabels([resolution_label(tick) for tick in ticks])


def plot(
    wave_numbers: np.ndarray,
    resolutions: np.ndarray,
    curves: np.ndarray,
    crossing_shells: np.ndarray,
    crossing_resolutions: np.ndarray,
    threshold: float,
    output: Path,
    title: str,
    dpi: int,
) -> None:
    freq = frequency_axis(resolutions)
    mean, std, minv, maxv = curve_stats(curves)

    fig, ax = plt.subplots(figsize=(7.0, 3.4))
    ax2 = ax.twinx()

    lower_sigma = mean - std
    upper_sigma = mean + std
    minmax = ax.fill_between(freq, minv, maxv, color="#9ecae1", alpha=0.45, linewidth=0, label="[min, max]")
    sigma = ax.fill_between(freq, lower_sigma, upper_sigma, color="#2b8cbe", alpha=0.55, linewidth=0, label="[+sigma, -sigma]")
    (mean_line,) = ax.plot(freq, mean, color="#0072bc", linewidth=1.9, label="mean")
    threshold_label = f"{threshold:.3g} threshold" if not crossing_shells.size else "_nolegend_"
    threshold_line = ax.axhline(threshold, color="#777777", linewidth=1.0, label=threshold_label)

    bar_handle = None
    if crossing_shells.size:
        shell_indices = np.asarray(
            [int(np.argmin(np.abs(wave_numbers - crossing))) for crossing in crossing_shells],
            dtype=int,
        )
        counts = np.bincount(shell_indices, minlength=len(freq)).astype(float)
        relative = counts / max(1.0, counts.sum())
        widths = shell_widths(freq) * 0.86
        bar_label = (
            f"{threshold:.3g} crossings\n"
            f"worst: {np.max(crossing_resolutions):.2f} A\n"
            f"best: {np.min(crossing_resolutions):.2f} A"
        )
        bar_handle = ax2.bar(
            freq,
            relative,
            width=widths,
            align="center",
            color="#41ab5d",
            alpha=0.55,
            edgecolor="#41ab5d",
            linewidth=0.4,
            label=bar_label,
        )

    ax.set_title(title)
    ax.set_xlabel("Resolution")
    ax.set_ylabel("FSC")
    ax.set_ylim(0.0, 1.0)
    ax.set_xlim(0.0, float(np.nanmax(freq)) * 1.02)
    ax.grid(True, color="#dddddd", linewidth=0.6, alpha=0.7)
    add_resolution_ticks(ax, freq)

    ax2.set_ylabel("Relative occurrence", color="#008a22")
    ax2.tick_params(axis="y", colors="#008a22")
    ax2.spines["right"].set_color("#008a22")
    ax2.set_ylim(0.0, 1.0)

    handles: list[object] = [mean_line, sigma, minmax]
    if threshold_line.get_label() != "_nolegend_":
        handles.append(threshold_line)
    labels = [handle.get_label() for handle in handles]
    if bar_handle is not None:
        handles.append(bar_handle)
        labels.append(bar_handle.get_label())
    ax.legend(handles, labels, loc="upper right", fontsize=8, frameon=True, framealpha=0.9)

    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    curves_path, directions_path, summary_path, output_path = paths_from_args(args)
    threshold = read_threshold(summary_path, args.threshold)
    wave_numbers, resolutions, curves = read_curves(curves_path)
    crossing_shells, crossing_resolutions = read_directions(directions_path)
    plot(
        wave_numbers,
        resolutions,
        curves,
        crossing_shells,
        crossing_resolutions,
        threshold,
        output_path,
        args.title,
        args.dpi,
    )
    print(output_path)


if __name__ == "__main__":
    main()
