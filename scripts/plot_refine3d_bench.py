#!/usr/bin/env python3
"""Create a publication-ready summary from refine3D benchmark reports.

Usage:
    python3 scripts/plot_refine3d_bench.py BENCHMARK_DIRECTORY [OUTPUT_DIRECTORY]

The input directory is expected to contain the text reports written with
L_BENCH_GLOB: REFINE3D_STRATEGY_BENCH_ITER*.txt,
REFINE3D_STAGE_BENCH_ITER*.txt, REFINE3D_BENCH_ITER*.txt, and
VOLASSEMBLE_BENCH_ITER*.txt.  The script uses only the Python standard
library and writes:

    refine3d_benchmark_summary.csv  one row per iteration
    refine3d_benchmark.svg          vector figure suitable for publication
    refine3d_stage_entry.svg        stage-entry figure, when stage reports exist

The SVG intentionally distinguishes the strategy-level wall-clock phases from
worker and volassemble totals.  A strategy ``setup/init`` sample is not a
sigma-estimation measurement: it combines reprojection-model preparation with
the group-sigma consolidation performed by the strategy.
"""

from __future__ import annotations

import csv
import math
import re
import statistics
import sys
from collections import defaultdict
from html import escape
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


FAMILIES = {
    "REFINE3D_STRATEGY_BENCH_ITER": "strategy",
    "REFINE3D_STAGE_BENCH_ITER": "stage",
    "REFINE3D_BENCH_ITER": "worker",
    "VOLASSEMBLE_BENCH_ITER": "assembly",
}

STRATEGY_COMPONENTS = (
    ("refine3D strategy setup/init", "setup / init (inclusive)", "#6f7d8c"),
    ("refine3D probabilistic pre-step", "probability-table preparation", "#2b6cb0"),
    ("refine3D matcher/scheduler", "matching / scheduling", "#c05621"),
    ("refine3D assembly/postprocess", "assembly / postprocessing", "#2f855a"),
)

NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
ITER_RE = re.compile(r"ITER(\d+)", re.IGNORECASE)
HEADER_RE = re.compile(r"^\s*\*{3}\s*(.*?)\s*\*{3}\s*$")
VALUE_RE = re.compile(rf"^\s*(.+?)\s*:\s*({NUMBER})\s*$")


def parse_file(path: Path) -> Tuple[Dict[str, str], Dict[str, float]]:
    """Return the context and seconds blocks without treating context as time."""
    context: Dict[str, str] = {}
    seconds: Dict[str, float] = {}
    section = ""
    for raw in path.read_text().splitlines():
        header = HEADER_RE.match(raw)
        if header:
            title = " ".join(header.group(1).upper().split())
            section = "context" if title == "BENCHMARK CONTEXT" else "seconds" if title.endswith("TIMINGS (S)") else ""
            continue
        match = VALUE_RE.match(raw)
        if not match:
            continue
        key, value = " ".join(match.group(1).split()), match.group(2)
        if section == "context":
            context[key] = value
        elif section == "seconds" and "%" not in key:
            seconds[key] = float(value)
    return context, seconds


def load_reports(indir: Path) -> Dict[int, Dict[str, object]]:
    rows: Dict[int, Dict[str, object]] = defaultdict(dict)
    for path in sorted(indir.glob("*_BENCH_ITER*.txt")):
        family = next((name for prefix, name in FAMILIES.items() if path.name.startswith(prefix)), None)
        if family is None:
            continue
        match = ITER_RE.search(path.name)
        if match is None:
            continue
        iteration = int(match.group(1))
        context, seconds = parse_file(path)
        row = rows[iteration]
        row[f"{family}_file"] = path.name
        for key, value in context.items():
            row[f"{family}:{key}"] = value
        for key, value in seconds.items():
            row[key] = value
    if not rows:
        raise ValueError(f"No recognized refine3D benchmark reports in {indir}")
    return dict(rows)


def as_float(row: Dict[str, object], key: str) -> float:
    value = row.get(key, 0.0)
    return float(value) if value not in (None, "") else 0.0


def context_value(row: Dict[str, object], family: str, key: str) -> str:
    return str(row.get(f"{family}:{key}", ""))


def write_summary(rows: Dict[int, Dict[str, object]], outpath: Path) -> None:
    columns = [
        "iteration", "nspace", "kto", "strategy_total_s", "strategy_setup_init_s",
        "reprojection_model_s", "group_sigma_consolidation_s", "probabilistic_prestep_s",
        "matcher_scheduler_s", "assembly_postprocess_s",
        "worker_total_s", "volassemble_total_s", "volassemble_nonuniform_filter_s",
        "stage_initialization_total_s", "stage_entry_calc_pspec_s",
    ]
    with outpath.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for iteration in sorted(rows):
            row = rows[iteration]
            writer.writerow({
                "iteration": iteration,
                "nspace": context_value(row, "strategy", "refine3D nspace") or
                          context_value(row, "stage", "refine3D nspace"),
                "kto": context_value(row, "strategy", "refine3D kto") or
                       context_value(row, "stage", "refine3D kto"),
                "strategy_total_s": as_float(row, "refine3D total time"),
                "strategy_setup_init_s": as_float(row, "refine3D strategy setup/init"),
                "reprojection_model_s": as_float(row, "refine3D reprojection model"),
                "group_sigma_consolidation_s": as_float(row, "refine3D group-sigma consolidation"),
                "probabilistic_prestep_s": as_float(row, "refine3D probabilistic pre-step"),
                "matcher_scheduler_s": as_float(row, "refine3D matcher/scheduler"),
                "assembly_postprocess_s": as_float(row, "refine3D assembly/postprocess"),
                "worker_total_s": as_float(row, "match3D total time"),
                "volassemble_total_s": as_float(row, "volassemble total time"),
                "volassemble_nonuniform_filter_s": as_float(row, "volassemble nonuniform_filter"),
                "stage_initialization_total_s": as_float(row, "refine3D stage initialization total"),
                "stage_entry_calc_pspec_s": as_float(row, "refine3D stage-entry calc_pspec"),
            })


def linear(value: float, low: float, high: float, pix_low: float, pix_high: float) -> float:
    if high <= low:
        return (pix_low + pix_high) / 2.0
    return pix_low + (value - low) * (pix_high - pix_low) / (high - low)


def nice_limit(values: Iterable[float]) -> float:
    maximum = max(values, default=1.0)
    if maximum <= 0:
        return 1.0
    step = 10 ** math.floor(math.log10(maximum))
    return math.ceil(maximum / step * 5.0) / 5.0 * step


def svg_text(x: float, y: float, text: str, size: int = 16, anchor: str = "start", weight: int = 400) -> str:
    return (f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" font-family="Arial,Helvetica,sans-serif" '
            f'fill="#17202a" font-size="{size}" font-weight="{weight}">{escape(text)}</text>')


def svg_vertical_text(x: float, y: float, text: str, size: int = 16) -> str:
    return (f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="middle" font-size="{size}" '
            f'font-family="Arial,Helvetica,sans-serif" fill="#17202a" font-weight="400" '
            f'transform="rotate(-90 {x:.1f} {y:.1f})">{escape(text)}</text>')


def draw_legend(items: Iterable[Tuple[str, str]], x: float, y: float) -> List[str]:
    out: List[str] = []
    xpos = x
    for label, color in items:
        width = 28 + len(label) * 12.2
        out.append(f'<rect x="{xpos:.1f}" y="{y - 17:.1f}" width="18" height="18" fill="{color}"/>')
        out.append(svg_text(xpos + 26, y, label, 20))
        xpos += width
    return out


def draw_line_legend(items: Iterable[Tuple[str, str]], x: float, y: float) -> List[str]:
    """Draw a legend whose samples match the line traces used in a panel."""
    out: List[str] = []
    xpos = x
    for label, color in items:
        width = 38 + len(label) * 12.2
        out.append(f'<line x1="{xpos:.1f}" y1="{y - 7:.1f}" x2="{xpos + 26:.1f}" y2="{y - 7:.1f}" '
                   f'stroke="{color}" stroke-width="2.667"/>')
        out.append(svg_text(xpos + 34, y, label, 20))
        xpos += width
    return out


def stage_boundaries(rows: Dict[int, Dict[str, object]]) -> List[Tuple[int, str]]:
    prior = None
    result = []
    for iteration in sorted(rows):
        row = rows[iteration]
        nspace = context_value(row, "strategy", "refine3D nspace")
        if nspace and nspace != prior:
            result.append((iteration, nspace))
            prior = nspace
    return result


def write_figure(rows: Dict[int, Dict[str, object]], outpath: Path) -> None:
    its = sorted(rows)
    first, last = its[0], its[-1]
    width, height = 1800, 1400
    margin_l, margin_r = 120, 70
    top_a, bottom_a = 175, 670
    top_b, bottom_b = 930, 1255
    plot_l, plot_r = margin_l, width - margin_r
    stack_max = nice_limit([sum(as_float(rows[i], metric) for metric, _, _ in STRATEGY_COMPONENTS) for i in its])
    total_max = nice_limit([max(as_float(rows[i], "match3D total time"), as_float(rows[i], "volassemble total time")) for i in its])

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<title>refine3D benchmark profile</title>',
        '<desc>Stacked per-iteration strategy wall times and independently timed worker and volassemble wall times.</desc>',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#17202a}.axis{stroke:#4a5568;stroke-width:2.667}.grid{stroke:#cbd5e0;stroke-width:1}.muted{fill:#4a5568}</style>',
        svg_text(plot_l, 58, "refine3D benchmark profile", 38, weight=500),
        svg_text(plot_l, 96, "Wall time by refinement iteration. Strategy components are stacked; worker and volassemble totals are shown separately.", 22),
    ]

    # Top panel: stacked strategy components.
    out += [svg_text(plot_l, top_a - 25, "A  Strategy-level wall time", 26, weight=500)]
    for tick in range(0, int(stack_max) + 1, max(1, int(math.ceil(stack_max / 5)))):
        y = linear(tick, 0, stack_max, bottom_a, top_a)
        out += [f'<line class="grid" x1="{plot_l}" y1="{y:.1f}" x2="{plot_r}" y2="{y:.1f}"/>', svg_text(plot_l - 16, y + 7, f"{tick:g}", 20, "end")]
    out += [
        f'<line class="axis" x1="{plot_l}" y1="{top_a}" x2="{plot_l}" y2="{bottom_a}"/>',
        f'<line class="axis" x1="{plot_l}" y1="{bottom_a}" x2="{plot_r}" y2="{bottom_a}"/>',
        svg_vertical_text(42, (top_a + bottom_a) / 2, "wall time (s)", 24),
    ]
    bar_w = max(2.0, (plot_r - plot_l) / max(1, last - first + 1) * 0.74)
    for iteration in its:
        x = linear(iteration, first, last, plot_l, plot_r)
        y = bottom_a
        for metric, _, color in STRATEGY_COMPONENTS:
            value = as_float(rows[iteration], metric)
            h = bottom_a - linear(value, 0, stack_max, bottom_a, top_a)
            y -= h
            out.append(f'<rect x="{x - bar_w/2:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{color}"/>')
    for iteration, nspace in stage_boundaries(rows):
        if iteration == first:
            continue
        x = linear(iteration, first, last, plot_l, plot_r)
        out += [f'<line x1="{x:.1f}" y1="{top_a}" x2="{x:.1f}" y2="{bottom_a}" stroke="#1a202c" stroke-width="2.667" stroke-dasharray="4 4"/>',
                svg_text(x + 9, top_a + 23, f"nspace {nspace}", 18)]
    for iteration in range(first, last + 1):
        if iteration not in rows or (iteration - first) % 10 != 0:
            continue
        x = linear(iteration, first, last, plot_l, plot_r)
        out += [f'<line class="axis" x1="{x:.1f}" y1="{bottom_a}" x2="{x:.1f}" y2="{bottom_a + 6}"/>',
                svg_text(x, bottom_a + 31, str(iteration), 20, "middle")]
    out += [
        svg_text((plot_l + plot_r) / 2, bottom_a + 62, "refinement iteration", 24, "middle"),
        svg_text(plot_l, bottom_a + 94, "Colour key — stacked contributions to strategy-level wall time:", 20, weight=500),
    ]
    out += draw_legend([(label, color) for _, label, color in STRATEGY_COMPONENTS], plot_l, bottom_a + 124)

    # Bottom panel: totals from independently timed executables.
    out += [
        svg_text(plot_l, top_b - 70, "B  Independently timed executables", 26, weight=500),
        svg_text(plot_l, top_b - 43,
                 "These lines are separate executable wall-clock totals; they are not added to the stacked bars above.", 19),
    ]
    for tick in range(0, int(total_max) + 1, max(1, int(math.ceil(total_max / 5)))):
        y = linear(tick, 0, total_max, bottom_b, top_b)
        out += [f'<line class="grid" x1="{plot_l}" y1="{y:.1f}" x2="{plot_r}" y2="{y:.1f}"/>', svg_text(plot_l - 16, y + 7, f"{tick:g}", 20, "end")]
    out += [
        f'<line class="axis" x1="{plot_l}" y1="{top_b}" x2="{plot_l}" y2="{bottom_b}"/>',
        f'<line class="axis" x1="{plot_l}" y1="{bottom_b}" x2="{plot_r}" y2="{bottom_b}"/>',
        svg_vertical_text(42, (top_b + bottom_b) / 2, "wall time per process (s)", 24),
    ]
    lines = (("match3D total time", "match3D worker executable", "#805ad5"),
             ("volassemble total time", "volume-assembly executable", "#d69e2e"))
    for metric, _, color in lines:
        points = []
        for iteration in its:
            value = as_float(rows[iteration], metric)
            x = linear(iteration, first, last, plot_l, plot_r)
            y = linear(value, 0, total_max, bottom_b, top_b)
            points.append(f"{x:.1f},{y:.1f}")
        out.append(f'<polyline fill="none" stroke="{color}" stroke-width="2.667" points="{" ".join(points)}"/>')
    for iteration in range(first, last + 1):
        if iteration not in rows or (iteration - first) % 10 != 0:
            continue
        x = linear(iteration, first, last, plot_l, plot_r)
        out += [f'<line class="axis" x1="{x:.1f}" y1="{bottom_b}" x2="{x:.1f}" y2="{bottom_b + 8}"/>', svg_text(x, bottom_b + 34, str(iteration), 20, "middle")]
    out.append(svg_text((plot_l + plot_r) / 2, height - 35, "refinement iteration", 24, "middle"))
    out += draw_line_legend([(label, color) for _, label, color in lines], plot_l, top_b - 18)
    n_strategy = len([i for i in its if "refine3D total time" in rows[i]])
    n_worker = len([i for i in its if "match3D total time" in rows[i]])
    n_assembly = len([i for i in its if "volassemble total time" in rows[i]])
    out += [
        svg_text(plot_r, height - 82, f"strategy n={n_strategy}; worker n={n_worker}; volassemble n={n_assembly}", 18, "end"),
        '</svg>',
    ]
    outpath.write_text("\n".join(out) + "\n")


def write_stage_figure(rows: Dict[int, Dict[str, object]], outpath: Path) -> None:
    """Draw stage initialization separately so calc_pspec remains legible."""
    samples = [(iteration, row) for iteration, row in sorted(rows.items())
               if "refine3D stage initialization total" in row]
    if not samples:
        return
    width, height = 1200, 680
    left, right, top, bottom = 120, 70, 115, 550
    values = [as_float(row, "refine3D stage initialization total") for _, row in samples]
    values += [as_float(row, "refine3D stage-entry calc_pspec") for _, row in samples]
    ymax = nice_limit(values)
    plot_w = width - left - right
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<title>refine3D stage-entry benchmark</title>',
        '<desc>Stage initialization wall time and the independently measured calc_pspec wall time.</desc>',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#17202a}.axis{stroke:#4a5568;stroke-width:2.667}.grid{stroke:#cbd5e0;stroke-width:1}</style>',
        svg_text(left, 45, "refine3D stage-entry benchmark", 28, weight=500),
        svg_text(left, 72, "Power-spectrum estimation is measured independently from the total strategy initialization path.", 16),
    ]
    for tick in range(0, int(ymax) + 1, max(1, int(math.ceil(ymax / 5)))):
        y = linear(tick, 0, ymax, bottom, top)
        out += [f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{width-right}" y2="{y:.1f}"/>',
                svg_text(left - 13, y + 5, f"{tick:g}", 13, "end")]
    out += [
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{bottom}"/>',
        f'<line class="axis" x1="{left}" y1="{bottom}" x2="{width-right}" y2="{bottom}"/>',
        svg_vertical_text(38, (top + bottom) / 2, "wall time (s)", 15),
    ]
    step = plot_w / max(1, len(samples))
    bar_w = min(70.0, step * 0.32)
    for index, (iteration, row) in enumerate(samples):
        center = left + step * (index + 0.5)
        total = as_float(row, "refine3D stage initialization total")
        pspec = as_float(row, "refine3D stage-entry calc_pspec")
        for offset, value, color in ((-bar_w * 0.6, total, "#6f7d8c"), (bar_w * 0.6, pspec, "#2b6cb0")):
            y = linear(value, 0, ymax, bottom, top)
            out.append(f'<rect x="{center+offset-bar_w/2:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{bottom-y:.1f}" fill="{color}"/>')
        nspace = context_value(row, "stage", "refine3D nspace")
        out.append(svg_text(center, bottom + 25, f"iter {iteration}", 13, "middle"))
        if nspace:
            out.append(svg_text(center, bottom + 44, f"nspace {nspace}", 12, "middle"))
    out += draw_legend((("stage initialization total", "#6f7d8c"), ("stage-entry calc_pspec", "#2b6cb0")), left, bottom + 83)
    out += [svg_text((left + width - right) / 2, height - 25, "stage first iteration", 15, "middle"), '</svg>']
    outpath.write_text("\n".join(out) + "\n")


def main(argv: List[str]) -> int:
    if not 2 <= len(argv) <= 3:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    indir = Path(argv[1]).expanduser().resolve()
    outdir = Path(argv[2]).expanduser().resolve() if len(argv) == 3 else indir
    if not indir.is_dir():
        print(f"Benchmark directory does not exist: {indir}", file=sys.stderr)
        return 2
    outdir.mkdir(parents=True, exist_ok=True)
    rows = load_reports(indir)
    summary = outdir / "refine3d_benchmark_summary.csv"
    figure = outdir / "refine3d_benchmark.svg"
    stage_figure = outdir / "refine3d_stage_entry.svg"
    write_summary(rows, summary)
    write_figure(rows, figure)
    write_stage_figure(rows, stage_figure)
    print(f"Parsed {len(rows)} benchmark iterations from {indir}")
    print(f"Wrote {summary}")
    print(f"Wrote {figure}")
    if stage_figure.exists():
        print(f"Wrote {stage_figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
