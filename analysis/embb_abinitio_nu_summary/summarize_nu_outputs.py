#!/usr/bin/env python3
"""Parse abinitio3D nonuniform-filter logs and render validation diagrams.

The script intentionally uses only the Python standard library so the figures
can be regenerated in a lean SIMPLE development environment.
"""

from __future__ import annotations

import csv
import html
import math
import os
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path


INPUT_DIR = Path("/Users/elmlundho/abinitio_outputs/embb_ab3D_nu_outputs")
OUT_DIR = Path(__file__).resolve().parent


def parse_float(text: str) -> float:
    return float(text.replace("%", ""))


def parse_int(text: str) -> int:
    return int(text.replace(",", ""))


def run_sort_key(path: Path) -> tuple[int, str]:
    match = re.search(r"RESTART(\d+)", path.name)
    return (int(match.group(1)) if match else 10**9, path.name)


def run_label(path: Path) -> str:
    match = re.search(r"RESTART(\d+)", path.name)
    return f"R{int(match.group(1)):02d}" if match else path.name


def base_record(path: Path, stage: int | None, stage_lp: float | None, iteration: int | None) -> dict:
    return {
        "run": run_label(path),
        "file": path.name,
        "stage": stage,
        "stage_lp": stage_lp,
        "iteration": iteration,
        "nu_iter": None,
        "mask_voxels": None,
        "before_base_voxels": None,
        "before_base_pct": None,
        "before_aux_voxels": None,
        "before_aux_pct": None,
        "after_base_voxels": None,
        "after_base_pct": None,
        "after_aux_voxels": None,
        "after_aux_pct": None,
        "aux_resolution": None,
        "aux_nearest_lp": None,
        "aux_margin_mean": None,
        "aux_margin_win": None,
        "aux_margin_wins": None,
        "smoothing_beta": None,
        "smoothing_maxits": None,
        "smoothing_candidates": None,
        "smoothing_auxiliary": None,
        "smoothing_step_tolerance": None,
        "smoothing_neighborhood": None,
        "smoothing_color_passes": None,
        "smoothing_quad_frac": None,
        "smoothing_initial_energy": None,
        "smoothing_changed_total": None,
        "smoothing_final_energy": None,
        "candidate_coords": "",
        "voxels_analyzed": None,
        "local_mean_A": None,
        "local_median_A": None,
        "local_sigma_A": None,
        "local_min_A": None,
        "local_max_A": None,
        "base_bank_voxels": None,
        "disc_voxels": None,
        "disc_voxels_pct": None,
        "neighbor_pairs": None,
        "identical_pairs": None,
        "identical_pairs_pct": None,
        "tolerated_pairs": None,
        "tolerated_pairs_pct": None,
        "discontinuous_pairs": None,
        "discontinuous_pairs_pct": None,
        "max_step_diff": None,
        "matching_lp_A": None,
        "fsc143_A": None,
        "filter_mode": None,
        "orientation_overlap": None,
    }


def parse_logs(input_dir: Path) -> tuple[list[dict], list[dict], list[dict]]:
    records_by_key: dict[tuple[str, int | None, int | None], dict] = {}
    bank_rows: list[dict] = []
    step_rows: list[dict] = []

    for path in sorted(input_dir.glob("*"), key=run_sort_key):
        if not path.is_file():
            continue
        stage: int | None = None
        stage_lp: float | None = None
        iteration: int | None = None
        current: dict | None = None
        section: str | None = None
        bank_idx = 0
        changed_total = 0

        def rec() -> dict:
            nonlocal current
            key = (path.name, stage, iteration)
            if key not in records_by_key:
                records_by_key[key] = base_record(path, stage, stage_lp, iteration)
            current = records_by_key[key]
            return current

        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                stripped = line.strip()

                m = re.match(r">>> STAGE\s+(\d+)\s+WITH LP\s+([0-9.]+)", stripped)
                if m:
                    stage = int(m.group(1))
                    stage_lp = float(m.group(2))
                    section = None
                    current = None
                    continue

                m = re.match(r">>> ITERATION\s+(\d+)", stripped)
                if m:
                    iteration = int(m.group(1))
                    section = None
                    current = None
                    changed_total = 0
                    continue

                if "NU auxiliary unary margins" in stripped:
                    rec()
                    section = "aux_margins"
                    continue
                if "NU candidate source assignments before" in stripped:
                    rec()
                    section = "candidate_before"
                    continue
                if "NU candidate source assignments after" in stripped:
                    rec()
                    section = "candidate_after"
                    continue
                if "NU FILTER LOCAL RESOLUTION SUMMARY" in stripped:
                    rec()
                    section = "resolution_summary"
                    continue
                if "NU LOW-PASS ASSIGNMENTS" in stripped:
                    rec()
                    section = "bank"
                    bank_idx = 0
                    continue
                if "NU AUXILIARY SOURCE ASSIGNMENTS" in stripped:
                    rec()
                    section = "aux_summary"
                    continue
                if "NU NEIGHBOR CONTINUITY" in stripped:
                    rec()
                    section = "continuity"
                    continue
                if "LP-step difference distribution" in stripped:
                    rec()
                    section = "step_distribution"
                    continue
                if stripped.startswith("**** SIMPLE_VOLASSEMBLE") or stripped.startswith("**** SIMPLE_POSTPROCESS"):
                    section = None

                if "NU ordered-label smoothing:" in stripped:
                    r = rec()
                    m = re.search(
                        r"beta=\s*([0-9.Ee+-]+), max iterations=(\d+), candidates=(\d+), auxiliary=(\d+), step tolerance=(\d+)",
                        stripped,
                    )
                    if m:
                        r["smoothing_beta"] = float(m.group(1))
                        r["smoothing_maxits"] = int(m.group(2))
                        r["smoothing_candidates"] = int(m.group(3))
                        r["smoothing_auxiliary"] = int(m.group(4))
                        r["smoothing_step_tolerance"] = int(m.group(5))
                    section = "smoothing"
                    continue

                if current is not None:
                    if "NU ordered-label smoothing neighborhood" in stripped:
                        m = re.search(r"neighborhood:\s*([^,]+), color passes=(\d+)", stripped)
                        if m:
                            current["smoothing_neighborhood"] = m.group(1).strip()
                            current["smoothing_color_passes"] = int(m.group(2))
                        continue
                    if "NU ordered-label smoothing quadratic jump fraction" in stripped:
                        current["smoothing_quad_frac"] = float(stripped.rsplit(":", 1)[1])
                        continue
                    if "NU ordered-label smoothing candidate coordinates" in stripped:
                        current["candidate_coords"] = " ".join(stripped.rsplit(":", 1)[1].split())
                        continue
                    if "NU ordered-label smoothing initial mean site energy" in stripped:
                        current["smoothing_initial_energy"] = float(stripped.rsplit(":", 1)[1])
                        continue
                    m = re.search(r"NU ordered-label smoothing iteration\s+(\d+) changed voxels:\s+(\d+), mean site energy:\s*([0-9.Ee+-]+)", stripped)
                    if m:
                        changed = int(m.group(2))
                        changed_total += changed
                        current[f"smoothing_changed_{m.group(1)}"] = changed
                        current[f"smoothing_energy_{m.group(1)}"] = float(m.group(3))
                        current["smoothing_changed_total"] = changed_total
                        current["smoothing_final_energy"] = float(m.group(3))
                        continue

                    m = re.search(r"MATCHING\s+LOW-PASS LIMIT AVG/SDEV/MIN/MAX:\s*([0-9.]+)", stripped)
                    if m:
                        current["matching_lp_A"] = float(m.group(1))
                        continue
                    m = re.search(r"RESOLUTION @ FSC=0\.143\s+AVG/SDEV/MIN/MAX:\s*([0-9.]+)", stripped)
                    if m:
                        current["fsc143_A"] = float(m.group(1))
                        continue
                    m = re.search(r"ORIENTATION OVERLAP:\s*([0-9.]+)", stripped)
                    if m:
                        current["orientation_overlap"] = float(m.group(1))
                        continue
                    m = re.search(r"\|\s*FILTER MODE\s*\|\s*([^|]+)", stripped)
                    if m:
                        current["filter_mode"] = m.group(1).strip()
                        continue

                if current is None:
                    continue

                if section in ("candidate_before", "candidate_after"):
                    prefix = "before" if section == "candidate_before" else "after"
                    m = re.match(r"Mask voxels:\s+(\d+)", stripped)
                    if m:
                        current["mask_voxels"] = int(m.group(1))
                        continue
                    m = re.match(r"Base-bank voxels:\s+(\d+)\s+\(\s*([0-9.]+)%\)", stripped)
                    if m:
                        current[f"{prefix}_base_voxels"] = int(m.group(1))
                        current[f"{prefix}_base_pct"] = float(m.group(2))
                        continue
                    m = re.match(r"Auxiliary voxels:\s+(\d+)\s+\(\s*([0-9.]+)%\)", stripped)
                    if m:
                        current[f"{prefix}_aux_voxels"] = int(m.group(1))
                        current[f"{prefix}_aux_pct"] = float(m.group(2))
                        continue
                    m = re.match(r"Aux\S*\s+([0-9.]+)\s+([0-9.]+)\s+(\d+)\s+([0-9.]+)%", stripped)
                    if m:
                        current["aux_resolution"] = float(m.group(1))
                        current["aux_nearest_lp"] = float(m.group(2))
                        current[f"{prefix}_aux_voxels"] = int(m.group(3))
                        current[f"{prefix}_aux_pct"] = float(m.group(4))
                        continue

                elif section == "aux_margins":
                    m = re.match(r"Aux\S*\s+([0-9.]+)\s+([0-9.]+)\s+(\d+)\s+([0-9.]+)%\s+([-0-9.Ee+]+)\s+([-0-9.Ee+]+)", stripped)
                    if m:
                        current["aux_resolution"] = float(m.group(1))
                        current["aux_nearest_lp"] = float(m.group(2))
                        current["aux_margin_wins"] = int(m.group(3))
                        current["aux_margin_mean"] = float(m.group(5))
                        current["aux_margin_win"] = float(m.group(6))
                        continue

                elif section == "resolution_summary":
                    m = re.match(r"Voxels analyzed:\s+(\d+)", stripped)
                    if m:
                        current["voxels_analyzed"] = int(m.group(1))
                        continue
                    m = re.match(r"Angstrom\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)", stripped)
                    if m:
                        current["local_mean_A"] = float(m.group(1))
                        current["local_median_A"] = float(m.group(2))
                        current["local_sigma_A"] = float(m.group(3))
                        current["local_min_A"] = float(m.group(4))
                        current["local_max_A"] = float(m.group(5))
                        continue

                elif section == "bank":
                    m = re.match(r"Base-bank voxels:\s+(\d+)", stripped)
                    if m:
                        current["base_bank_voxels"] = int(m.group(1))
                        continue
                    if not stripped or stripped.startswith("Bank") or stripped.startswith("LP limit"):
                        continue
                    parts = stripped.replace("%", "").split()
                    try:
                        if len(parts) >= 5 and re.match(r"^\d+$", parts[0]):
                            row_bank = int(parts[0])
                            fourier_k = int(parts[1])
                            lp = float(parts[2])
                            voxels = int(parts[3])
                            pct = float(parts[4])
                        elif len(parts) >= 3:
                            bank_idx += 1
                            row_bank = bank_idx
                            fourier_k = None
                            lp = float(parts[0])
                            voxels = int(parts[1])
                            pct = float(parts[2])
                        else:
                            continue
                    except ValueError:
                        continue
                    bank_rows.append(
                        {
                            "run": current["run"],
                            "file": current["file"],
                            "stage": current["stage"],
                            "stage_lp": current["stage_lp"],
                            "iteration": current["iteration"],
                            "bank": row_bank,
                            "fourier_k": fourier_k,
                            "lp_A": lp,
                            "voxels": voxels,
                            "pct_base": pct,
                        }
                    )
                    continue

                elif section == "aux_summary":
                    m = re.match(r"Mask voxels:\s+(\d+)", stripped)
                    if m:
                        current["mask_voxels"] = int(m.group(1))
                        continue
                    m = re.match(r"Aux\S*@?\s+([0-9.]+)\s+(\d+)\s+([0-9.]+)%", stripped)
                    if m:
                        current["aux_resolution"] = float(m.group(1))
                        current["after_aux_voxels"] = int(m.group(2))
                        current["after_aux_pct"] = float(m.group(3))
                        continue

                elif section == "continuity":
                    m = re.match(r"Voxels analyzed:\s+(\d+)", stripped)
                    if m:
                        current["voxels_analyzed"] = int(m.group(1))
                        continue
                    m = re.match(r"Voxels with discontinuous neighbors:\s+(\d+)\s+\(\s*([0-9.]+)%\)", stripped)
                    if m:
                        current["disc_voxels"] = int(m.group(1))
                        current["disc_voxels_pct"] = float(m.group(2))
                        continue
                    m = re.match(r"Neighbor pairs examined:\s+(\d+)", stripped)
                    if m:
                        current["neighbor_pairs"] = int(m.group(1))
                        continue
                    m = re.match(r"identical:\s+(\d+)\s+\(\s*([0-9.]+)%\)", stripped)
                    if m:
                        current["identical_pairs"] = int(m.group(1))
                        current["identical_pairs_pct"] = float(m.group(2))
                        continue
                    m = re.match(r"tolerated:\s+(\d+)\s+\(\s*([0-9.]+)%\)", stripped)
                    if m:
                        current["tolerated_pairs"] = int(m.group(1))
                        current["tolerated_pairs_pct"] = float(m.group(2))
                        continue
                    m = re.match(r"discontinuous:\s+(\d+)\s+\(\s*([0-9.]+)%\)", stripped)
                    if m:
                        current["discontinuous_pairs"] = int(m.group(1))
                        current["discontinuous_pairs_pct"] = float(m.group(2))
                        continue

                elif section == "step_distribution":
                    parts = stripped.split()
                    if len(parts) >= 4 and parts[0].isdigit():
                        step = int(parts[0])
                        step_rows.append(
                            {
                                "run": current["run"],
                                "file": current["file"],
                                "stage": current["stage"],
                                "stage_lp": current["stage_lp"],
                                "iteration": current["iteration"],
                                "step": step,
                                "pairs": int(parts[1]),
                                "pct": float(parts[2]),
                                "class": parts[3],
                            }
                        )
                        current["max_step_diff"] = max(current["max_step_diff"] or 0, step)
                        continue

    records = [r for r in records_by_key.values() if r.get("voxels_analyzed") is not None]
    records.sort(key=lambda r: (r["run"], r["iteration"] or -1))

    by_run: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        by_run[r["run"]].append(r)
    for rows in by_run.values():
        rows.sort(key=lambda r: r["iteration"] or -1)
        for i, r in enumerate(rows, start=1):
            r["nu_iter"] = i

    return records, bank_rows, step_rows


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        keys: list[str] = []
        seen = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    keys.append(key)
                    seen.add(key)
        fieldnames = keys
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: object, ndigits: int = 2) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float):
        return f"{value:.{ndigits}f}"
    return str(value)


def svg_text(x: float, y: float, text: object, size: int = 13, fill: str = "#1f2933", anchor: str = "start", weight: str = "400", rotate: float | None = None) -> str:
    transform = f' transform="rotate({rotate:.1f} {x:.1f} {y:.1f})"' if rotate is not None else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{fill}" '
        f'text-anchor="{anchor}" font-weight="{weight}"{transform}>{html.escape(str(text))}</text>'
    )


def svg_rect(x: float, y: float, w: float, h: float, fill: str, stroke: str = "none", rx: float = 0, opacity: float = 1.0) -> str:
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx:.1f}" fill="{fill}" stroke="{stroke}" opacity="{opacity:.3f}"/>'


def color_lerp(c1: str, c2: str, frac: float) -> str:
    frac = max(0.0, min(1.0, frac))
    a = tuple(int(c1[i : i + 2], 16) for i in (1, 3, 5))
    b = tuple(int(c2[i : i + 2], 16) for i in (1, 3, 5))
    c = tuple(round(a[i] + (b[i] - a[i]) * frac) for i in range(3))
    return f"#{c[0]:02x}{c[1]:02x}{c[2]:02x}"


def pct_color(value: float | None, vmax: float, base: str = "#f3f7f4", high: str = "#116149") -> str:
    if value is None:
        return "#eef2f7"
    return color_lerp(base, high, math.sqrt(max(0.0, min(1.0, value / vmax))))


def path_from_points(points: list[tuple[float, float]]) -> str:
    if not points:
        return ""
    return "M " + " L ".join(f"{x:.1f},{y:.1f}" for x, y in points)


def grouped(records: list[dict]) -> dict[str, list[dict]]:
    by_run: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        by_run[r["run"]].append(r)
    for rows in by_run.values():
        rows.sort(key=lambda r: r["nu_iter"] or 0)
    return dict(sorted(by_run.items()))


def final_records(records: list[dict]) -> list[dict]:
    return [rows[-1] for rows in grouped(records).values() if rows]


def bank_for_record(bank_rows: list[dict], rec: dict) -> list[dict]:
    return [
        b
        for b in bank_rows
        if b["run"] == rec["run"] and b["iteration"] == rec["iteration"] and b["stage"] == rec["stage"]
    ]


def selected_finest_lp(bank_rows: list[dict], rec: dict) -> float | None:
    rows = [b for b in bank_for_record(bank_rows, rec) if b["voxels"] > 0]
    return min((b["lp_A"] for b in rows), default=None)


def dominant_lp(bank_rows: list[dict], rec: dict) -> tuple[float | None, float | None]:
    rows = bank_for_record(bank_rows, rec)
    if not rows:
        return None, None
    best = max(rows, key=lambda b: b["pct_base"])
    return best["lp_A"], best["pct_base"]


def make_dashboard(records: list[dict], bank_rows: list[dict], path: Path) -> None:
    finals = final_records(records)
    runs = [r["run"] for r in finals]

    cols = sorted({round(b["lp_A"], 1) for rec in finals for b in bank_for_record(bank_rows, rec)}, reverse=True)
    cols_with_aux: list[str] = [f"{c:.1f}A" for c in cols] + ["Aux"]
    values: dict[tuple[str, str], float] = {}
    for rec in finals:
        base_to_pct: dict[float, float] = defaultdict(float)
        for b in bank_for_record(bank_rows, rec):
            base_to_pct[round(b["lp_A"], 1)] += b["pct_base"] * ((rec.get("after_base_pct") or 100.0) / 100.0)
        for c in cols:
            values[(rec["run"], f"{c:.1f}A")] = base_to_pct.get(c, 0.0)
        values[(rec["run"], "Aux")] = rec.get("after_aux_pct") or 0.0

    w, h = 1600, 1120
    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
        '<style>text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif}.small{font-size:12px}.label{font-size:13px}.axis{stroke:#9aa5b1;stroke-width:1}.grid{stroke:#d9e2ec;stroke-width:1}.panel{fill:#ffffff;stroke:#d9e2ec}.muted{fill:#627d98}</style>',
        svg_rect(0, 0, w, h, "#f7fafc"),
        svg_text(40, 52, "abinitio3D NU filtering validation across repeated runs", 28, "#102a43", weight="700"),
        svg_text(40, 80, f"{len(finals)} repeated runs; {len(records)} NU-filtered iterations parsed from SIMPLE abinitio3D logs", 15, "#52606d"),
    ]

    # Process strip.
    parts.append(svg_rect(40, 110, 1520, 130, "#ffffff", "#d9e2ec", 8))
    steps = [
        ("Even/Odd maps", "base Fourier volumes"),
        ("Discrete LP bank", "coarse-to-fine candidates"),
        ("Aux ML map", "competes in every voxel"),
        ("26-neighbor Potts", "step<=1 tolerated"),
        ("NU assignment", "local filter + continuity audit"),
    ]
    x0, y0, boxw, boxh, gap = 70, 138, 245, 62, 38
    for i, (title, subtitle) in enumerate(steps):
        x = x0 + i * (boxw + gap)
        fill = ["#e0f2f1", "#d8ecff", "#fcebc8", "#e9ddff", "#d7f5df"][i]
        parts.append(svg_rect(x, y0, boxw, boxh, fill, "#bcccdc", 8))
        parts.append(svg_text(x + boxw / 2, y0 + 27, title, 17, "#102a43", "middle", "700"))
        parts.append(svg_text(x + boxw / 2, y0 + 49, subtitle, 12, "#52606d", "middle"))
        if i < len(steps) - 1:
            ax = x + boxw + 8
            ay = y0 + boxh / 2
            parts.append(f'<path d="M {ax:.1f},{ay:.1f} L {ax+22:.1f},{ay:.1f}" stroke="#486581" stroke-width="2"/>')
            parts.append(f'<path d="M {ax+22:.1f},{ay:.1f} L {ax+15:.1f},{ay-6:.1f} M {ax+22:.1f},{ay:.1f} L {ax+15:.1f},{ay+6:.1f}" stroke="#486581" stroke-width="2" fill="none"/>')

    # Assignment heatmap.
    heat_x, heat_y = 90, 305
    cell_w, cell_h = 92, 32
    parts.append(svg_text(40, 275, "Final voxel-source assignment per run", 20, "#102a43", weight="700"))
    parts.append(svg_text(40, 296, "Base-bank values are normalized to mask voxels; Aux is the ML-regularized auxiliary assignment after Potts smoothing.", 13, "#52606d"))
    for j, col in enumerate(cols_with_aux):
        parts.append(svg_text(heat_x + j * cell_w + cell_w / 2, heat_y - 12, col, 12, "#52606d", "middle", rotate=-35))
    for i, run in enumerate(runs):
        y = heat_y + i * cell_h
        parts.append(svg_text(heat_x - 14, y + 21, run, 13, "#102a43", "end", "700"))
        for j, col in enumerate(cols_with_aux):
            value = values.get((run, col), 0.0)
            x = heat_x + j * cell_w
            color = pct_color(value, 45.0, high="#006d77" if col != "Aux" else "#9a3412")
            parts.append(svg_rect(x, y, cell_w - 2, cell_h - 2, color, "#ffffff", 2))
            txt_color = "#ffffff" if value > 22 else "#102a43"
            parts.append(svg_text(x + cell_w / 2, y + 20, f"{value:.1f}", 11, txt_color, "middle", "700" if value > 0.05 else "400"))
    parts.append(svg_text(heat_x + len(cols_with_aux) * cell_w + 12, heat_y + 18, "% of mask", 12, "#52606d"))

    # Final continuity and local resolution summary.
    panel_x, panel_y = 40, 690
    parts.append(svg_text(panel_x, panel_y - 22, "Continuity and local-resolution end state", 20, "#102a43", weight="700"))
    chart_x, chart_y, chart_w, chart_h = panel_x + 50, panel_y + 35, 650, 250
    max_disc = max((r.get("discontinuous_pairs_pct") or 0.0 for r in finals), default=1.0)
    max_vox = max((r.get("disc_voxels_pct") or 0.0 for r in finals), default=1.0)
    max_y = max(2.0, math.ceil(max(max_disc, max_vox) * 1.2))
    parts.append(svg_rect(panel_x, panel_y, 735, 325, "#ffffff", "#d9e2ec", 8))
    parts.append(svg_text(chart_x, chart_y - 12, "Discontinuity after Potts smoothing", 14, "#102a43", weight="700"))
    for k in range(0, int(max_y) + 1, max(1, int(max_y / 4) or 1)):
        yy = chart_y + chart_h - (k / max_y) * chart_h
        parts.append(f'<line x1="{chart_x:.1f}" y1="{yy:.1f}" x2="{chart_x+chart_w:.1f}" y2="{yy:.1f}" class="grid"/>')
        parts.append(svg_text(chart_x - 8, yy + 4, f"{k}%", 11, "#627d98", "end"))
    bw = chart_w / max(1, len(finals)) * 0.34
    for i, r in enumerate(finals):
        cx = chart_x + (i + 0.5) * chart_w / len(finals)
        pair = r.get("discontinuous_pairs_pct") or 0.0
        vox = r.get("disc_voxels_pct") or 0.0
        ph = (pair / max_y) * chart_h
        vh = (vox / max_y) * chart_h
        parts.append(svg_rect(cx - bw - 1, chart_y + chart_h - ph, bw, ph, "#006d77"))
        parts.append(svg_rect(cx + 1, chart_y + chart_h - vh, bw, vh, "#d97706"))
        parts.append(svg_text(cx, chart_y + chart_h + 18, r["run"], 11, "#52606d", "middle"))
    parts.append(svg_rect(chart_x + chart_w - 160, chart_y - 28, 12, 12, "#006d77"))
    parts.append(svg_text(chart_x + chart_w - 142, chart_y - 18, "pairs >1 LP step", 11, "#52606d"))
    parts.append(svg_rect(chart_x + chart_w - 160, chart_y - 10, 12, 12, "#d97706"))
    parts.append(svg_text(chart_x + chart_w - 142, chart_y, "voxels with discontinuity", 11, "#52606d"))

    table_x, table_y = 825, 690
    parts.append(svg_rect(table_x, table_y, 735, 325, "#ffffff", "#d9e2ec", 8))
    headers = ["Run", "Mean A", "Median A", "Finest base A", "Aux %", "Disc pairs %"]
    xs = [table_x + 28, table_x + 120, table_x + 225, table_x + 345, table_x + 505, table_x + 610]
    for x, head in zip(xs, headers):
        parts.append(svg_text(x, table_y + 34, head, 12, "#52606d", weight="700"))
    parts.append(f'<line x1="{table_x+20}" y1="{table_y+45}" x2="{table_x+715}" y2="{table_y+45}" stroke="#d9e2ec"/>')
    for i, r in enumerate(finals[:10]):
        y = table_y + 68 + i * 24
        finest = selected_finest_lp(bank_rows, r)
        vals = [
            r["run"],
            fmt(r.get("local_mean_A"), 2),
            fmt(r.get("local_median_A"), 2),
            fmt(finest, 1),
            fmt(r.get("after_aux_pct") or 0.0, 2),
            fmt(r.get("discontinuous_pairs_pct"), 2),
        ]
        for x, val in zip(xs, vals):
            parts.append(svg_text(x, y, val, 12, "#102a43", weight="700" if x == xs[0] else "400"))

    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def make_trajectories(records: list[dict], path: Path) -> None:
    by_run = grouped(records)
    w, h = 1500, 950
    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
        '<style>text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif}.grid{stroke:#d9e2ec;stroke-width:1}.axis{stroke:#9aa5b1;stroke-width:1}</style>',
        svg_rect(0, 0, w, h, "#f7fafc"),
        svg_text(40, 52, "NU filtering trajectories within abinitio3D", 28, "#102a43", weight="700"),
        svg_text(40, 78, "Each line is one repeated run; x-axis is the NU-filtered iteration count within the run.", 14, "#52606d"),
    ]
    palette = ["#006d77", "#d97706", "#3b82f6", "#7c3aed", "#047857", "#be123c", "#0f766e", "#9333ea", "#ca8a04", "#334155"]

    def line_panel(x: float, y: float, width: float, height: float, title: str, field: str, y_label: str, invert: bool = False) -> None:
        vals = [r.get(field) for rows in by_run.values() for r in rows if r.get(field) is not None]
        if not vals:
            return
        ymin, ymax = min(vals), max(vals)
        pad = (ymax - ymin) * 0.08 if ymax > ymin else 1.0
        ymin -= pad
        ymax += pad
        if field.endswith("_pct"):
            ymin = min(0.0, ymin)
        xmax = max((len(rows) for rows in by_run.values()), default=1)
        parts.append(svg_rect(x - 10, y - 38, width + 40, height + 76, "#ffffff", "#d9e2ec", 8))
        parts.append(svg_text(x, y - 14, title, 16, "#102a43", weight="700"))
        for i in range(5):
            frac = i / 4
            yy = y + height - frac * height
            val = ymin + frac * (ymax - ymin)
            if invert:
                val = ymax - frac * (ymax - ymin)
            parts.append(f'<line x1="{x:.1f}" y1="{yy:.1f}" x2="{x+width:.1f}" y2="{yy:.1f}" class="grid"/>')
            parts.append(svg_text(x - 8, yy + 4, f"{val:.2f}", 11, "#627d98", "end"))
        parts.append(f'<line x1="{x:.1f}" y1="{y+height:.1f}" x2="{x+width:.1f}" y2="{y+height:.1f}" class="axis"/>')
        parts.append(f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x:.1f}" y2="{y+height:.1f}" class="axis"/>')
        parts.append(svg_text(x + width / 2, y + height + 38, "NU iteration", 12, "#52606d", "middle"))
        parts.append(svg_text(x - 52, y + height / 2, y_label, 12, "#52606d", "middle", rotate=-90))
        for idx, (run, rows) in enumerate(by_run.items()):
            pts: list[tuple[float, float]] = []
            for r in rows:
                val = r.get(field)
                if val is None:
                    continue
                xi = x + ((r["nu_iter"] or 1) - 1) / max(1, xmax - 1) * width
                f = (val - ymin) / (ymax - ymin) if ymax > ymin else 0.5
                if invert:
                    f = 1.0 - f
                yi = y + height - f * height
                pts.append((xi, yi))
            color = palette[idx % len(palette)]
            parts.append(f'<path d="{path_from_points(pts)}" fill="none" stroke="{color}" stroke-width="2.0" opacity="0.82"/>')
            if pts:
                parts.append(f'<circle cx="{pts[-1][0]:.1f}" cy="{pts[-1][1]:.1f}" r="3.2" fill="{color}"/>')
        parts.append(svg_text(x + width - 5, y + 16, f"n={len(by_run)}", 11, "#627d98", "end"))

    line_panel(95, 150, 590, 260, "Local median resolution", "local_median_A", "A", invert=True)
    line_panel(845, 150, 590, 260, "Auxiliary ML-volume assignment", "after_aux_pct", "% mask")
    line_panel(95, 570, 590, 260, "Discontinuous neighbor pairs", "discontinuous_pairs_pct", "% pairs")
    line_panel(845, 570, 590, 260, "Potts smoothing changed voxels", "smoothing_changed_total", "voxels")

    # Legend.
    lx, ly = 40, 895
    for idx, run in enumerate(by_run):
        x = lx + idx * 140
        color = palette[idx % len(palette)]
        parts.append(f'<line x1="{x:.1f}" y1="{ly:.1f}" x2="{x+24:.1f}" y2="{ly:.1f}" stroke="{color}" stroke-width="3"/>')
        parts.append(svg_text(x + 32, ly + 4, run, 12, "#52606d"))

    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def make_report(records: list[dict], bank_rows: list[dict], path: Path) -> None:
    finals = final_records(records)
    mean_disc_pair = statistics.mean([r["discontinuous_pairs_pct"] for r in finals if r.get("discontinuous_pairs_pct") is not None])
    mean_disc_vox = statistics.mean([r["disc_voxels_pct"] for r in finals if r.get("disc_voxels_pct") is not None])
    mean_aux_before = statistics.mean([(r.get("before_aux_pct") or 0.0) for r in finals])
    mean_aux = statistics.mean([(r.get("after_aux_pct") or 0.0) for r in finals])
    mean_local = statistics.mean([r["local_mean_A"] for r in finals if r.get("local_mean_A") is not None])
    median_local = statistics.mean([r["local_median_A"] for r in finals if r.get("local_median_A") is not None])
    mean_changed = statistics.mean([r["smoothing_changed_total"] for r in finals if r.get("smoothing_changed_total") is not None])

    lines = [
        "# EMBB abinitio3D NU Filtering Summary",
        "",
        f"Parsed `{len(records)}` NU-filtered iterations from `{len(finals)}` repeated abinitio3D runs.",
        "",
        "## Final-Iteration Aggregate",
        "",
        f"- Mean base-bank local resolution across final maps: `{mean_local:.2f} A`",
        f"- Mean of final base-bank local-resolution medians: `{median_local:.2f} A`",
        f"- Mean auxiliary assignment before/after Potts smoothing: `{mean_aux_before:.2f}%` -> `{mean_aux:.2f}%` of mask voxels",
        f"- Mean Potts smoothing changed voxels in final iterations: `{mean_changed:.0f}`",
        f"- Mean discontinuous neighbor-pair rate: `{mean_disc_pair:.2f}%`",
        f"- Mean voxels with any discontinuous neighbor: `{mean_disc_vox:.2f}%`",
        "",
        "Note: SIMPLE reports the NU local-resolution summary over base-bank voxels. Auxiliary ML-regularized assignments are reported separately as a percentage of the support mask.",
        "",
        "## Final-Iteration Runs",
        "",
        "| Run | Stage | Iter | Mean A | Median A | Finest selected base A | Dominant base LP | Aux % | Disc pair % | Disc voxel % | FSC143 A | Matching LP A |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in finals:
        finest = selected_finest_lp(bank_rows, r)
        dominant, dominant_pct = dominant_lp(bank_rows, r)
        dominant_txt = "n/a" if dominant is None else f"{dominant:.1f} A ({dominant_pct:.1f}%)"
        lines.append(
            "| "
            + " | ".join(
                [
                    r["run"],
                    fmt(r.get("stage"), 0),
                    fmt(r.get("iteration"), 0),
                    fmt(r.get("local_mean_A"), 2),
                    fmt(r.get("local_median_A"), 2),
                    fmt(finest, 1),
                    dominant_txt,
                    fmt(r.get("after_aux_pct") or 0.0, 2),
                    fmt(r.get("discontinuous_pairs_pct"), 2),
                    fmt(r.get("disc_voxels_pct"), 2),
                    fmt(r.get("fsc143_A"), 3),
                    fmt(r.get("matching_lp_A"), 3),
                ]
            )
            + " |"
        )
    lines += [
        "",
        "## Generated Figures",
        "",
        "- `nu_validation_dashboard.svg`: final assignment heatmap plus final continuity/local-resolution summary.",
        "- `nu_trajectories.svg`: per-run NU-stage trajectories for median resolution, auxiliary assignment, discontinuity, and Potts changed voxels.",
        "",
        "## Parsed Data",
        "",
        "- `nu_iteration_metrics.csv`: one row per NU-filtered iteration.",
        "- `nu_bank_assignments.csv`: one row per low-pass bank member per NU-filtered iteration.",
        "- `nu_step_distribution.csv`: one row per LP-step-difference class per NU-filtered iteration.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    input_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else INPUT_DIR
    if not input_dir.exists():
        raise SystemExit(f"Input directory does not exist: {input_dir}")

    records, bank_rows, step_rows = parse_logs(input_dir)
    if not records:
        raise SystemExit(f"No NU filter records parsed from: {input_dir}")

    write_csv(OUT_DIR / "nu_iteration_metrics.csv", records)
    write_csv(OUT_DIR / "nu_bank_assignments.csv", bank_rows)
    write_csv(OUT_DIR / "nu_step_distribution.csv", step_rows)
    make_dashboard(records, bank_rows, OUT_DIR / "nu_validation_dashboard.svg")
    make_trajectories(records, OUT_DIR / "nu_trajectories.svg")
    make_report(records, bank_rows, OUT_DIR / "nu_validation_summary.md")

    finals = final_records(records)
    print(f"Parsed {len(records)} NU-filtered iterations from {len(finals)} repeated runs.")
    print(f"Wrote outputs to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
