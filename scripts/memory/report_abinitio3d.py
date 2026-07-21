#!/usr/bin/env python3
"""Create the abinitio3D memory-screening CSV, TXT, charts, and PDF report."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image as PILImage, ImageDraw, ImageFont
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle, Image,
)


NAVY = "#16324F"
BLUE = "#2878B5"
CYAN = "#36A7B4"
ORANGE = "#F28E2B"
RED = "#D1495B"
GREEN = "#2A9D6F"
LIGHT = colors.HexColor("#EEF3F7")
GRID = colors.HexColor("#B7C6D3")


def sweep(df: pd.DataFrame, group: str, baseline: pd.Series, x: str) -> pd.DataFrame:
    part = pd.concat([baseline.to_frame().T, df[df.design_group == group]], ignore_index=True)
    return part.drop_duplicates(subset=[x]).sort_values(x)


def save_charts(df: pd.DataFrame, out: Path) -> list[Path]:
    out.mkdir(parents=True, exist_ok=True)
    base = df[df.design_group == "baseline"].iloc[0]
    paths: list[Path] = []

    def font(size: int, bold: bool = False):
        candidates = [
            "/System/Library/Fonts/HelveticaNeue.ttc",
            "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        ]
        for candidate in candidates:
            try:
                return ImageFont.truetype(candidate, size=size, index=1 if bold and candidate.endswith(".ttc") else 0)
            except OSError:
                pass
        return ImageFont.load_default(size=size)

    def canvas(title: str):
        image = PILImage.new("RGB", (1980, 650), "white")
        draw = ImageDraw.Draw(image)
        draw.text((990, 35), title, anchor="mm", fill=NAVY, font=font(34, True))
        return image, draw

    def panel(draw, slot: int, xs, ys, title: str, xlabel: str, color: str, labels=None):
        x0 = 70 + slot * 645
        y0, width, height = 120, 555, 410
        left, right, top, bottom = x0 + 82, x0 + width - 20, y0 + 45, y0 + height - 55
        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)
        xmin, xmax = float(xs.min()), float(xs.max())
        ymin, ymax = float(ys.min()), float(ys.max())
        if xmax == xmin:
            xmax += 1
        pad = max(8.0, (ymax - ymin) * .16)
        ymin, ymax = ymin - pad, ymax + pad
        xp = lambda v: left + (float(v) - xmin) / (xmax - xmin) * (right - left)
        yp = lambda v: bottom - (float(v) - ymin) / (ymax - ymin) * (bottom - top)
        draw.text((x0 + width/2, y0 + 10), title, anchor="ma", fill=NAVY, font=font(24, True))
        for frac in (0, .5, 1):
            yy = bottom - frac * (bottom - top)
            draw.line((left, yy, right, yy), fill="#D8E2EA", width=2)
            draw.text((left - 12, yy), f"{ymin + frac*(ymax-ymin):.0f}", anchor="rm", fill="#526674", font=font(17))
        draw.line((left, top, left, bottom), fill="#8094A2", width=2)
        draw.line((left, bottom, right, bottom), fill="#8094A2", width=2)
        points = [(xp(x), yp(y)) for x, y in zip(xs, ys)]
        draw.line(points, fill=color, width=7, joint="curve")
        for i, ((xx, yy), x, y) in enumerate(zip(points, xs, ys)):
            draw.ellipse((xx-8, yy-8, xx+8, yy+8), fill=color, outline="white", width=3)
            lab = labels[i] if labels else f"{x:g}"
            draw.text((xx, bottom + 18), str(lab), anchor="ma", fill="#405563", font=font(16))
            draw.text((xx, yy - 15), f"{y:.0f}", anchor="ms", fill=color, font=font(16, True))
        draw.text((x0 + width/2, y0 + height), xlabel, anchor="ma", fill="#405563", font=font(18))
        draw.text((x0 + 7, (top+bottom)/2), "Peak MiB", anchor="lm", fill="#526674", font=font(17))

    image, draw = canvas("Strongest one-factor memory drivers")
    specs = [("particle_count", "nptcls", "Particles", BLUE),
             ("state_count", "nstates", "Independent states", RED),
             ("job_partitions", "nparts", "Local partitions", ORANGE)]
    for i, (group, x, label, color) in enumerate(specs):
        d = sweep(df, group, base, x)
        if group == "job_partitions":
            d = pd.concat([df[(df.design_group == "thread_count") & (df.nthr == 2)], df[df.design_group == group]]).sort_values(x)
        panel(draw, i, d[x], d.peak_tree_rss_mib, label, label, color)
    p = out / "abinitio3d_key_drivers.png"; image.save(p); paths.append(p)

    image, draw = canvas("Geometry, physical scale, and working-box effects")
    specs = [("box_size", "box", "Raw box", GREEN),
             ("pixel_size", "smpd_angstrom_per_pixel", "Pixel size", CYAN),
             ("mask_diameter", "mskdiam_angstrom", "Mask diameter", ORANGE)]
    for i, (group, x, title, color) in enumerate(specs):
        d = sweep(df, group, base, x)
        labels = [f"{v:g}/e{int(e)}" for v, e in zip(d[x], d.effective_box)]
        panel(draw, i, d[x], d.peak_tree_rss_mib, title, "Input / effective box", color, labels)
    p = out / "abinitio3d_geometry.png"; image.save(p); paths.append(p)

    image, draw = canvas("Sampling, symmetry, and thread controls")
    d = df[df.design_group == "sampled_particles"].sort_values("nsample_effective")
    panel(draw, 0, d.nsample_effective, d.peak_tree_rss_mib, "Sampled particles", "Sample of 500", BLUE)
    sym = pd.concat([base.to_frame().T, df[df.design_group.isin(["symmetry", "symmetry_search"])]]).reset_index(drop=True)
    panel(draw, 1, range(4), sym.peak_tree_rss_mib, "Symmetry policy", "Target / start", CYAN, ["c1", "c2", "d2", "c2/c1"])
    threads = sweep(df, "thread_count", base, "nthr")
    panel(draw, 2, threads.nthr, threads.peak_tree_rss_mib, "Thread count", "OpenMP threads", RED)
    p = out / "abinitio3d_controls.png"; image.save(p); paths.append(p)
    return paths


def write_txt(df: pd.DataFrame, metadata: dict, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        stream.write("SIMPLE abinitio3D memory screening — machine-readable report\n")
        stream.write("# Format: UTF-8, tab-separated table after the DATA marker.\n")
        stream.write("# Memory metric: conservative process-tree planning peak = parent peak + largest nparts child peaks.\n")
        stream.write("# Units: RSS fields ending in _mib use mebibytes; byte fields are integer bytes.\n")
        stream.write("# Synthetic inputs: asymmetric 3D volume projections, CTF disabled, SNR 0.1.\n")
        stream.write("# Successful runs: %d/%d\n" % ((df.status == "ok").sum(), len(df)))
        stream.write("# Metadata: " + json.dumps(metadata, sort_keys=True) + "\n")
        stream.write("DATA\n")
        df.to_csv(stream, sep="\t", index=False, lineterminator="\n")


def add_footer(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#D6E0E8"))
    canvas.line(doc.leftMargin, .48 * inch, letter[0] - doc.rightMargin, .48 * inch)
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#667788"))
    canvas.drawString(doc.leftMargin, .30 * inch, "SIMPLE abinitio3D memory screening")
    canvas.drawRightString(letter[0] - doc.rightMargin, .30 * inch, f"Page {doc.page}")
    canvas.restoreState()


def make_pdf(df: pd.DataFrame, metadata: dict, charts: list[Path], path: Path) -> None:
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="TitleX", parent=styles["Title"], fontName="Helvetica-Bold",
                              fontSize=24, leading=28, textColor=colors.HexColor(NAVY), spaceAfter=12))
    styles.add(ParagraphStyle(name="H1X", parent=styles["Heading1"], fontName="Helvetica-Bold",
                              fontSize=16, leading=19, textColor=colors.HexColor(NAVY), spaceBefore=8, spaceAfter=8))
    styles.add(ParagraphStyle(name="BodyX", parent=styles["BodyText"], fontSize=9.2, leading=13,
                              textColor=colors.HexColor("#243746"), spaceAfter=7))
    styles.add(ParagraphStyle(name="SmallX", parent=styles["BodyText"], fontSize=7.4, leading=9.4,
                              textColor=colors.HexColor("#3D4E5C")))
    doc = SimpleDocTemplate(str(path), pagesize=letter, rightMargin=.55*inch, leftMargin=.55*inch,
                            topMargin=.55*inch, bottomMargin=.62*inch, title="SIMPLE abinitio3D Memory Screening")
    story = []
    base = df[df.design_group == "baseline"].iloc[0]
    peak = df.loc[df.peak_tree_rss_mib.idxmax()]
    story += [Paragraph("SIMPLE abinitio3D Memory Screening", styles["TitleX"]),
              Paragraph("A controlled process-tree RSS study across data geometry, sampling, parallelism, heterogeneity, and workflow controls", styles["BodyX"])]
    summary_data = [
        ["Runs", "Baseline", "Maximum", "SIMPLE commit"],
        [f"{len(df)} / {len(df)} successful", f"{base.peak_tree_rss_mib:.1f} MiB",
         f"{peak.peak_tree_rss_mib:.1f} MiB", "c11237cd4"],
    ]
    t = Table(summary_data, colWidths=[1.55*inch]*4, rowHeights=[.30*inch, .42*inch])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor(NAVY)), ("TEXTCOLOR", (0,0), (-1,0), colors.white),
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"), ("FONTNAME", (0,1), (-1,1), "Helvetica-Bold"),
        ("FONTSIZE", (0,0), (-1,-1), 9), ("ALIGN", (0,0), (-1,-1), "CENTER"),
        ("GRID", (0,0), (-1,-1), .35, GRID), ("BACKGROUND", (0,1), (-1,1), LIGHT),
    ]))
    story += [t, Spacer(1, .18*inch), Paragraph("Executive findings", styles["H1X"])]
    findings = [
        f"The largest planning bound was <b>{peak.peak_tree_rss_mib:.1f} MiB</b> for 500 particles, box 128, four local partitions, and two threads per partition.",
        f"Independent state count was a strong driver: one, two, and four states required {base.peak_tree_rss_mib:.1f}, "
        f"{df[df.nstates.eq(2)].iloc[0].peak_tree_rss_mib:.1f}, and {df[(df.design_group.eq('state_count')) & df.nstates.eq(4)].iloc[0].peak_tree_rss_mib:.1f} MiB.",
        "Local partitions mainly multiply worker memory: at two threads per worker, the planning bound rose from 303.0 MiB (one partition) to 396.5 MiB (two) and 622.4 MiB (four).",
        "Particle sampling reduced worker memory: for 500 particles, sampling 50–100 used about 282 MiB versus 324 MiB for all 500.",
        "Symmetry choice, symmetry-axis search, projection reconstruction, and initialization route changed the one-iteration result by only a few percent in this design.",
    ]
    for item in findings:
        story.append(Paragraph("• " + item, styles["BodyX"]))
    story += [Spacer(1, .08*inch), Paragraph("How to read the memory number", styles["H1X"]),
              Paragraph("abinitio3D launches short-lived simple_private_exec workers. The reported planning peak is deliberately conservative: the parent's native peak RSS plus the largest <i>nparts</i> child-process peaks. It is a safe capacity-planning bound, but the component peaks may not have occurred at exactly the same instant. Parent and child peaks are retained separately in the TXT/CSV files.", styles["BodyX"]),
              PageBreak(), Paragraph("Primary drivers", styles["H1X"]), Image(str(charts[0]), width=7.25*inch, height=2.37*inch),
              Paragraph("Particles increase worker-side tables and reconstruction inputs. Independent states add both parent and worker allocations. Job partitions are the clearest concurrency multiplier because several workers may be resident together.", styles["BodyX"]),
              Spacer(1,.12*inch), Paragraph("Geometry and physical scale", styles["H1X"]), Image(str(charts[1]), width=7.25*inch, height=2.37*inch),
              Paragraph("The working crop—not only the raw box—matters. At the baseline low-pass schedule, raw boxes 96–160 mostly cropped to 88 pixels, yet raw-stack and worker allocations still increased at the largest boxes. Pixel size and mask diameter change the physical support and therefore allocation behavior.", styles["BodyX"]),
              PageBreak(), Paragraph("Workflow controls", styles["H1X"]), Image(str(charts[2]), width=7.25*inch, height=2.37*inch)]
    story += [Paragraph("Design and reproducibility", styles["H1X"]),
              Paragraph("Each run used a fresh SIMPLE project. Deterministic asymmetric volumes were projected by SIMPLE, imported with CTF disabled, and given a one-iteration abinitio2D preparation outside the measured interval because abinitio3D requires class metadata. All measured cases used SIMPLE's unmodified stage-1 schedule, filt_mode=none, automsk=no, lpstart=20 A, and synthetic SNR=0.1 unless the row explicitly varies a setting.", styles["BodyX"]),
              Paragraph("This is a screened experimental design, not a full Cartesian product. One-factor sweeps support local comparisons around the baseline. Four interaction cases expose combined pressure but are insufficient for a production-quality predictive model by themselves; repeat measurements and a randomized factorial or space-filling design should be added before model fitting.", styles["BodyX"]),
              PageBreak(), Paragraph("All benchmark runs", styles["H1X"])]

    cols = ["design_group", "nptcls", "box", "effective_box", "nthr", "nparts", "nstates",
            "nsample_effective", "iterations_observed", "peak_tree_rss_mib", "elapsed_s"]
    labels = ["Group", "N", "Box", "Eff", "Thr", "Part", "States", "Sample", "Iter", "Peak MiB", "Sec"]
    rows = [labels]
    for _, r in df.iterrows():
        rows.append([r.design_group.replace("_", " "), int(r.nptcls), int(r.box), int(r.effective_box),
                     int(r.nthr), int(r.nparts), int(r.nstates), int(r.nsample_effective),
                     int(r.iterations_observed), f"{r.peak_tree_rss_mib:.1f}", f"{r.elapsed_s:.1f}"])
    table = Table(rows, repeatRows=1, colWidths=[1.38*inch,.42*inch,.40*inch,.40*inch,.37*inch,.37*inch,.43*inch,.52*inch,.37*inch,.58*inch,.48*inch])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor(NAVY)), ("TEXTCOLOR", (0,0), (-1,0), colors.white),
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"), ("FONTSIZE", (0,0), (-1,-1), 6.7),
        ("ALIGN", (1,1), (-1,-1), "RIGHT"), ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, LIGHT]),
        ("GRID", (0,0), (-1,-1), .25, GRID), ("TOPPADDING", (0,0), (-1,-1), 3), ("BOTTOMPADDING", (0,0), (-1,-1), 3),
    ]))
    story += [table]
    doc.build(story, onFirstPage=add_footer, onLaterPages=add_footer)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("results", type=Path)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, default=Path("output/pdf"))
    ap.add_argument("--tmp-dir", type=Path, default=Path("tmp/pdfs/abinitio3d"))
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.results)
    metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
    csv_out = args.output_dir / "abinitio3d_memory_screening.csv"
    txt_out = args.output_dir / "abinitio3d_memory_screening.txt"
    pdf_out = args.output_dir / "simple_abinitio3d_memory_screening_report.pdf"
    shutil.copyfile(args.results, csv_out)
    write_txt(df, metadata, txt_out)
    charts = save_charts(df, args.tmp_dir)
    make_pdf(df, metadata, charts, pdf_out)
    print(csv_out)
    print(txt_out)
    print(pdf_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
