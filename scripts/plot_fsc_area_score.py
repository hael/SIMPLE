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
 
    # ── Palette ───────────────────────────────────────────────────────────────
    C_BLUE_DARK   = "#185FA5"   # mean line
    C_BLUE_MID    = "#378ADD"   # ±σ fill
    C_BLUE_LIGHT  = "#B5D4F4"   # min/max fill
    C_GREEN_BAR   = "#3B6D11"   # crossing bars
    C_GREEN_AXIS  = "#3B6D11"   # right-axis labels/spine
    C_THRESHOLD   = "#888780"   # dashed threshold line
    C_GRID        = "#E8E7E2"   # subtle grid
 
    # ── Figure & axes ─────────────────────────────────────────────────────────
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Inter", "Helvetica Neue", "Arial", "DejaVu Sans", "sans-serif"],
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
 
    fig, ax = plt.subplots(figsize=(8.0, 3.8))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax2 = ax.twinx()
    ax2.set_facecolor("white")
    ax2.set_prop_cycle(color=[C_GREEN_BAR])
 
    # ── Bands & mean ──────────────────────────────────────────────────────────
    lower_sigma = mean - std
    upper_sigma = mean + std
 
    minmax = ax.fill_between(
        freq, minv, maxv,
        color=C_BLUE_LIGHT, alpha=0.50, linewidth=0,
        label="[min, max]",
    )
    sigma = ax.fill_between(
        freq, lower_sigma, upper_sigma,
        color=C_BLUE_MID, alpha=0.40, linewidth=0,
        label="[+sigma, -sigma]",
    )
    (mean_line,) = ax.plot(
        freq, mean,
        color=C_BLUE_DARK, linewidth=2.2,
        solid_capstyle="round", label="mean",
    )
    threshold_label = f"{threshold:.3g} threshold" if not crossing_shells.size else "_nolegend_"
    threshold_line = ax.axhline(
        threshold,
        color=C_THRESHOLD, linewidth=0.9,
        linestyle=(0, (5, 4)),
        label=threshold_label,
    )
 
    # ── Crossing bars ─────────────────────────────────────────────────────────
    bar_handle = None
    if crossing_shells.size:
        shell_indices = np.asarray(
            [int(np.argmin(np.abs(wave_numbers - c))) for c in crossing_shells],
            dtype=int,
        )
        counts = np.bincount(shell_indices, minlength=len(freq)).astype(float)
        relative = counts / max(1.0, counts.sum())
        widths = shell_widths(freq) * 0.86
        bar_label = (
            f"{threshold:.3g} crossings\n"
            f"worst: {np.max(crossing_resolutions):.2f} Å\n"
            f"best:  {np.min(crossing_resolutions):.2f} Å"
        )
        bar_handle = ax2.bar(
            freq, relative,
            width=widths, align="center",
            color=C_GREEN_BAR, alpha=0.50,
            edgecolor="none",
            label=bar_label,
        )
 
    # ── Axes cosmetics ────────────────────────────────────────────────────────
    max_freq = float(np.nanmax(freq))
 
    ax.set_xlim(0.0, max_freq * 1.02)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("Resolution", fontsize=11, labelpad=6, color="#444441")
    ax.set_ylabel("FSC", fontsize=11, labelpad=6, color="#444441")
    ax.tick_params(axis="both", labelsize=9.5, colors="#5F5E5A", length=3, width=0.6)
    add_resolution_ticks(ax, freq)
 
    ax.grid(True, axis="y", color=C_GRID, linewidth=0.7, zorder=0)
    ax.grid(False, axis="x")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#D3D1C7")
    ax.spines["bottom"].set_color("#D3D1C7")
    ax.spines["left"].set_linewidth(0.7)
    ax.spines["bottom"].set_linewidth(0.7)
 
    ax2.set_ylim(0.0, 1.0)
    ax2.set_ylabel("Relative occurrence", fontsize=11, labelpad=6, color=C_GREEN_AXIS)
    ax2.tick_params(axis="y", labelsize=9.5, colors=C_GREEN_AXIS, length=3, width=0.6)
    ax2.spines["right"].set_color(C_GREEN_AXIS)
    ax2.spines["right"].set_linewidth(0.7)
    ax2.spines["top"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.spines["bottom"].set_visible(False)
    ax2.grid(False)
    ax2.set_zorder(ax.get_zorder() - 1)
 
    # ── Layout ────────────────────────────────────────────────────────────────
    fig.subplots_adjust(top=0.68, bottom=0.13, left=0.09, right=0.91)
    ax_pos = ax.get_position()
    x0 = ax_pos.x0
    fig_w_in, fig_h_in = fig.get_size_inches()
 
    # ── Title (two-line header matching widget design) ────────────────────────
    parts = title.split(" ", 1)
    label_text = parts[0].upper()
    subtitle_text = parts[1] if len(parts) > 1 else ""
 
    fig.text(
        x0, ax_pos.y1 + 0.28,
        label_text,
        fontsize=8.5, fontweight="normal",
        color="#888780",
        transform=fig.transFigure,
        va="bottom",
    )
    fig.text(
        x0, ax_pos.y1 + 0.18,
        subtitle_text,
        fontsize=14, fontweight="medium",
        color="#2C2C2A",
        transform=fig.transFigure,
        va="bottom",
    )
 
    # ── Legend (horizontal row in header, swatches via patches) ──────────────
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D
 
    # Absolute sizes (inches) converted to figure fraction
    SW_IN  = 0.20   # swatch width
    SH_IN  = 0.11   # swatch height
    GAP_IN = 0.06   # swatch → label gap
    SEP_IN = 0.18   # label → next swatch gap
    sw = SW_IN / fig_w_in
    sh = SH_IN / fig_h_in
    gap = GAP_IN / fig_w_in
    sep = SEP_IN / fig_w_in
 
    # Legend row sits just above the axes top edge
    row_mid  = ax_pos.y1 + 0.055
    row_bot  = row_mid - sh / 2.0
 
    items = [
        ("line", C_BLUE_DARK,   1.00, "mean"),
        ("rect", C_BLUE_MID,    0.45, "±\u03c3"),
        ("rect", C_BLUE_LIGHT,  0.55, "min / max"),
    ]
    if crossing_shells.size:
        items.append(("rect", C_GREEN_BAR, 0.55,
                      f"{threshold:.3g} crossings"))
 
    # Pass 1 — draw placeholder labels at x=0 to measure their widths
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
 
    def label_width_frac(txt: str) -> float:
        t = fig.text(0, 0, txt, fontsize=8.5, va="center", ha="left",
                     transform=fig.transFigure)
        bb = t.get_window_extent(renderer=renderer)
        t.remove()
        return bb.width / (fig_w_in * fig.dpi)
 
    # Pre-compute widths so we can place accurately
    widths = [label_width_frac(label) for _, _, _, label in items]
 
    # Pass 2 — draw swatches + labels at correct positions
    x = ax_pos.x0
    for (kind, color, alpha, label), lw in zip(items, widths):
        if kind == "line":
            fig.add_artist(Line2D(
                [x, x + sw], [row_mid, row_mid],
                transform=fig.transFigure,
                color=color, linewidth=2.2,
                solid_capstyle="round", clip_on=False,
            ))
        else:
            fig.add_artist(mpatches.FancyBboxPatch(
                (x, row_bot), sw, sh,
                boxstyle="square,pad=0",
                transform=fig.transFigure,
                facecolor=color, alpha=alpha,
                edgecolor="none", clip_on=False,
            ))
        fig.text(
            x + sw + gap, row_mid, label,
            fontsize=8.5, color="#5F5E5A",
            va="center", ha="left",
            transform=fig.transFigure,
        )
        x += sw + gap + lw + sep
 
    # Worst / best stats — right-aligned
    if crossing_shells.size:
        stats = (
            f"worst  {np.max(crossing_resolutions):.2f} Å"
            f"   ·   best  {np.min(crossing_resolutions):.2f} Å"
        )
        fig.text(
            ax_pos.x1, row_mid, stats,
            fontsize=8.5, color="#5F5E5A",
            va="center", ha="right",
            transform=fig.transFigure,
        )
 
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, facecolor="white")
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