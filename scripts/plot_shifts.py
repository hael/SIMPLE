#!/usr/bin/env python3
"""PDF report of the particle origin shifts in a STAR file or a SIMPLE project.

Page 1: the summary statistics, a histogram of the x shifts, a histogram of the y shifts and a
2D histogram of x against y. The bins are 1 A wide (--binwidth), centred on whole Angstroms,
and the plots are symmetric: +/- the largest |shift| (or --limit), rounded out to a bin edge.
The x and y histograms share that range and one count axis, and the 2D histogram is square,
so the two directions compare at a glance.

For a SIMPLE project, the following pages repeat the three plots per group, the most
populated first (--max-groups), each row led by a thumbnail of the group:
  - --2d: per 2D class, with its class average (the project's cavg stack);
  - --3d: per projection direction (ptcl3D proj, per state), with a projection of the state's
    volume (the project's registered volumes, or --vols) along the direction of the group's
    particle closest to the group's mean viewing direction, made by SIMPLE's reproject.
The shift axes are the page-1 axes on every row; each row's counts are scaled to its group.

Inputs:
  - a STAR file: the first loop block (or --block) holding rlnOriginXAngst/rlnOriginYAngst
    (Angstroms; RELION 3.1+ and SIMPLE's export), else rlnOriginX/rlnOriginY (pixels;
    converted to Angstroms with --angpix, otherwise plotted in pixels with 1-pixel bins);
  - a SIMPLE project (.simple): the x/y shifts of ptcl2D (--2d) or ptcl3D (--3d), stored in
    pixels of each particle's stack and converted with that stack's smpd, read through
    SIMPLE's print_project_vals. Only selected particles (state > 0) unless
    --include-deselected. SIMPLE's executables are found on PATH, or in --simple-bin.

Needs numpy and matplotlib for the plots; the reading, the statistics and the thumbnails'
images are plain Python, and the statistics are printed before the plotting starts.

Examples:
  python scripts/plot_shifts.py particles.star
  python scripts/plot_shifts.py project.simple --2d
  python scripts/plot_shifts.py project.simple --3d --max-groups 100 --log2d
  python scripts/plot_shifts.py project.simple --3d --vols sharp_state01.mrc sharp_state02.mrc
"""

from __future__ import annotations

import argparse
import array
import math
import os
import shutil
import statistics
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

ANGST = "Å"
# the STAR shift columns, the preferred pair first, with their units
SHIFT_COLUMNS = (
    ("rlnOriginXAngst", "rlnOriginYAngst", ANGST),
    ("rlnOriginX", "rlnOriginY", "px"),
)
ROWS_PER_PAGE = 5

# the reference palette's light-mode roles: one blue for a single series, its ramp for
# magnitude (light = near zero), text and hairlines in neutral inks
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
HAIRLINE = "#e4e3df"
SERIES = "#2a78d6"
SEQUENTIAL = ("#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
              "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="PDF report of the origin shifts in a STAR file or a SIMPLE project.")
    parser.add_argument("input", type=Path, help="STAR file, or SIMPLE project file (.simple)")
    field = parser.add_mutually_exclusive_group()
    field.add_argument("--2d", dest="field", action="store_const", const="ptcl2D",
                       help="project: the shifts of ptcl2D, and per 2D class")
    field.add_argument("--3d", dest="field", action="store_const", const="ptcl3D",
                       help="project: the shifts of ptcl3D, and per projection direction")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="output PDF (default: <input stem>[_ptcl2D|_ptcl3D]_shifts.pdf beside the input)")
    parser.add_argument("--binwidth", type=float, default=1.0,
                        help="bin width, in the shifts' unit (default: 1, i.e. 1 A)")
    parser.add_argument("--limit", type=float, default=None,
                        help="plot shifts within +/- this value (default: the largest |shift|)")
    parser.add_argument("--log2d", action="store_true",
                        help="logarithmic colour scale for the 2D histograms (shows sparse tails)")
    parser.add_argument("--block", default=None,
                        help="STAR: data block to read, without 'data_' (default: the first loop block with shifts)")
    parser.add_argument("--angpix", type=float, default=None,
                        help="STAR: pixel size (A) to convert pixel shifts (rlnOriginX/Y) to Angstroms")
    parser.add_argument("--include-deselected", action="store_true",
                        help="project: include particles with state <= 0 (default: excluded)")
    parser.add_argument("--max-groups", type=int, default=50,
                        help="project: plot this many classes/directions, the most populated first "
                             "(default: 50; 0: none, only page 1)")
    parser.add_argument("--all-groups", action="store_true",
                        help="project: plot every class/direction")
    parser.add_argument("--vols", nargs="+", type=Path, default=None,
                        help="project, --3d: the volumes for the thumbnails, one per state in order "
                             "(default: the project's registered volumes)")
    parser.add_argument("--simple-bin", type=Path, default=None,
                        help="project: the directory with simple_exec and simple_private_exec (default: PATH)")
    return parser.parse_args()


# ---------------- STAR files ----------------

def read_star_shifts(path: Path, block: str | None) -> tuple[str, list[float], list[float], str]:
    """The block name, the x and y shifts and their unit, read in one pass from the requested
    loop block, or from the first loop block that holds a shift column pair (SHIFT_COLUMNS).

    Only the two shift columns of that block are kept, so a large particle file costs two
    floats per particle. Labels lose their leading underscore and any trailing '#N'; fields
    are split on whitespace; key-value pairs outside loops are skipped.
    """
    name = None          # the current data block
    labels: list[str] = []
    in_loop = False
    in_header = False
    target = None        # (ix, iy, unit) while reading the chosen loop
    xs: list[float] = []
    ys: list[float] = []
    blocks: list[str] = []

    def pick() -> tuple[int, int, str] | None:
        if block is not None and name != block:
            return None
        for xlab, ylab, unit in SHIFT_COLUMNS:
            if xlab in labels and ylab in labels:
                return labels.index(xlab), labels.index(ylab), unit
        return None

    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for nline, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("data_") or line.startswith("loop_") or (in_loop and not in_header and line.startswith("_")):
                if target is not None:
                    break  # the chosen loop has ended
                if line.startswith("data_"):
                    name = line[5:]
                    blocks.append(name)
                    in_loop = in_header = False
                else:
                    in_loop = in_header = True
                labels = [] if line.startswith(("data_", "loop_")) else [line[1:].split()[0]]
                continue
            if in_loop and in_header and line.startswith("_"):
                labels.append(line[1:].split()[0])
                continue
            if in_loop:
                if in_header:
                    in_header = False
                    target = pick()
                if target is None:
                    continue
                ix, iy, _ = target
                fields = line.split()
                if len(fields) <= max(ix, iy):
                    raise SystemExit(f"line {nline}: {len(fields)} fields, the shifts are in columns {ix + 1} and {iy + 1}")
                xs.append(float(fields[ix]))
                ys.append(float(fields[iy]))
    if target is not None:
        return name, xs, ys, target[2]
    if block is not None and block not in blocks:
        raise SystemExit(f"no data block 'data_{block}' in the STAR file; blocks: {', '.join(blocks) or 'none'}")
    wanted = " or ".join(f"{x}/{y}" for x, y, _ in SHIFT_COLUMNS)
    where = f"data_{block}" if block is not None else "any loop block"
    raise SystemExit(f"no shift columns ({wanted}) in {where}")


def load_star(args: argparse.Namespace) -> dict:
    project_only = [("--2d/--3d", args.field is not None), ("--vols", args.vols is not None),
                    ("--simple-bin", args.simple_bin is not None), ("--all-groups", args.all_groups),
                    ("--include-deselected", args.include_deselected)]
    given = [name for name, on in project_only if on]
    if given:
        raise SystemExit(f"{', '.join(given)} apply to a SIMPLE project, not a STAR file")
    block, xs, ys, unit = read_star_shifts(args.input, args.block)
    if args.angpix is not None:
        if unit == "px":
            xs = [v * args.angpix for v in xs]
            ys = [v * args.angpix for v in ys]
            unit = ANGST
        else:
            print(f"note: the shifts are already in {ANGST}; --angpix is not used", file=sys.stderr)
    elif unit == "px":
        print("note: pixel shifts (rlnOriginX/Y): plotted in pixels; give --angpix for Angstroms", file=sys.stderr)
    return {"source": f"{args.input}   data_{block}", "xs": xs, "ys": ys, "unit": unit, "groups": []}


# ---------------- SIMPLE projects ----------------

def simple_tool(args: argparse.Namespace, name: str) -> str:
    """The path of SIMPLE's executable @name, from --simple-bin or PATH."""
    if args.simple_bin is not None:
        path = (args.simple_bin / name).resolve()  # absolute: reproject runs in a temporary folder
        if not path.is_file():
            raise SystemExit(f"no {name} in {args.simple_bin}")
        return str(path)
    found = shutil.which(name)
    if found is None:
        raise SystemExit(f"{name} is not on PATH; give SIMPLE's bin directory with --simple-bin")
    return found


def run(cmd: list[str], cwd: str | None = None) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(cmd, text=True, capture_output=True, check=False, cwd=cwd)
    except FileNotFoundError:
        raise SystemExit(f"cannot run '{cmd[0]}'")


def print_project_vals(exe: str, projfile: Path, oritype: str, keys: list[str],
                       required: bool = True) -> list[tuple[int, int, list[float]]] | None:
    """The rows print_project_vals writes for @oritype: (index, state, [key values]), 'eulnorm'
    giving three. None, when not @required, if the segment lacks a key."""
    cmd = [exe, "prg=print_project_vals", f"projfile={projfile}", f"oritype={oritype}", f"keys={','.join(keys)}"]
    out = run(cmd)
    nvals = sum(3 if k == "eulnorm" else 1 for k in keys)
    rows = []
    for line in out.stdout.splitlines():
        fields = line.split()
        if len(fields) != 2 + nvals:
            continue
        try:
            rows.append((int(fields[0]), int(fields[1]), [float(v) for v in fields[2:]]))
        except ValueError:
            continue
    if not rows:
        said = [l.strip() for l in (out.stdout + out.stderr).splitlines() if "missing" in l or "ERROR" in l]
        if not required:
            return None
        detail = "; ".join(said) if said else f"exit status {out.returncode}"
        raise SystemExit(f"no {oritype} rows with {','.join(keys)} in {projfile} ({detail})")
    return rows


def project_out_entries(exe: str, projfile: Path) -> list[dict[str, str]]:
    """The entries of the project's out segment (print_project_field), as key -> value."""
    out = run([exe, "prg=print_project_field", f"projfile={projfile}", "oritype=out"])
    entries = []
    for line in out.stdout.splitlines():
        pairs = [tok.split("=", 1) for tok in line.split() if "=" in tok]
        if pairs:
            entries.append({k: v for k, v in pairs})
    return entries


def load_project(args: argparse.Namespace) -> dict:
    if args.field is None:
        raise SystemExit("a SIMPLE project needs --2d (ptcl2D shifts) or --3d (ptcl3D shifts)")
    if args.block is not None or args.angpix is not None:
        raise SystemExit("--block and --angpix apply to STAR files; a project's pixel sizes are its stacks'")
    if args.vols is not None and args.field != "ptcl3D":
        raise SystemExit("--vols are the --3d thumbnails' volumes")
    private = simple_tool(args, "simple_private_exec")
    smpd = {idx: vals[0] for idx, _, vals in print_project_vals(private, args.input, "stk", ["smpd"])}
    plot_groups = args.all_groups or args.max_groups > 0
    keys = ["x", "y", "stkind"]
    if plot_groups:
        keys += ["class"] if args.field == "ptcl2D" else ["proj", "e1", "e2", "eulnorm"]
    rows = print_project_vals(private, args.input, args.field, keys, required=not plot_groups)
    if rows is None:
        what = "class" if args.field == "ptcl2D" else "proj (a 3D refinement's directions)"
        print(f"note: {args.field} has no {what}: only page 1 is plotted", file=sys.stderr)
        plot_groups = False
        keys = ["x", "y", "stkind"]
        rows = print_project_vals(private, args.input, args.field, keys)
    xs, ys, gkeys, eulers, normals = [], [], [], [], []
    for idx, state, vals in rows:
        if state <= 0 and not args.include_deselected:
            continue
        x, y, stkind = vals[:3]
        istk = int(round(stkind))
        if istk not in smpd:
            raise SystemExit(f"{args.field} particle {idx}: stack {istk} is not in the project's stk segment")
        xs.append(x * smpd[istk])
        ys.append(y * smpd[istk])
        if plot_groups:
            if args.field == "ptcl2D":
                gkeys.append(int(round(vals[3])))
            else:
                gkeys.append((state, int(round(vals[3]))))
                eulers.append((vals[4], vals[5]))
                normals.append(tuple(vals[6:9]))
    which = "all particles" if args.include_deselected else "state > 0"
    data = {"source": f"{args.input}   {args.field} ({which})", "xs": xs, "ys": ys, "unit": ANGST, "groups": []}
    if plot_groups and xs:
        data["groups"] = select_groups(args, gkeys)
        if data["groups"]:
            out = project_out_entries(simple_tool(args, "simple_exec"), args.input)
            if args.field == "ptcl2D":
                thumbs_2d(data["groups"], out)
            else:
                thumbs_3d(args, data["groups"], eulers, normals, out)
    return data


def select_groups(args: argparse.Namespace, gkeys: list) -> list[dict]:
    """The groups to plot: classes (2D) or (state, direction)s (3D) with a particle, the most
    populated first, unassigned ones (class or direction 0) left out."""
    members: dict = {}
    for i, key in enumerate(gkeys):
        if (key if isinstance(key, int) else key[1]) <= 0:
            continue
        members.setdefault(key, []).append(i)
    order = sorted(members, key=lambda k: (-len(members[k]), k))
    nstates = len({k[0] for k in members}) if order and not isinstance(order[0], int) else 1
    if not args.all_groups:
        order = order[:args.max_groups]
    groups = []
    for key in order:
        if isinstance(key, int):
            label = f"class {key}"
        elif nstates > 1:
            label = f"state {key[0]} · direction {key[1]}"
        else:
            label = f"direction {key[1]}"
        groups.append({"key": key, "label": label, "members": members[key], "thumb": None, "ntotal": len(members)})
    return groups


# ---------------- thumbnails ----------------

def read_mrc_image(path: str, index: int) -> tuple[int, int, list[float]] | None:
    """Image @index (0-based) of an MRC stack: (nx, ny, values row by row); None when the file
    is missing, the index outside the stack, or the mode not one of 0, 1, 2, 6, 12."""
    try:
        with open(path, "rb") as f:
            head = f.read(1024)
            if len(head) < 1024:
                return None
            little = head[212] != 0x11
            end = "<" if little else ">"
            nx, ny, nz, mode = struct.unpack(end + "4i", head[:16])
            nsymbt = struct.unpack(end + "i", head[92:96])[0]
            types = {0: ("b", 1), 1: ("h", 2), 2: ("f", 4), 6: ("H", 2), 12: ("e", 2)}
            if mode not in types or not 0 <= index < nz:
                return None
            code, size = types[mode]
            f.seek(1024 + nsymbt + index * nx * ny * size)
            raw = f.read(nx * ny * size)
            if len(raw) < nx * ny * size:
                return None
    except OSError:
        return None
    vals = list(struct.unpack(f"{end}{nx * ny}{code}", raw))
    return nx, ny, vals


def mrc_pixel_size(path: str) -> float | None:
    """The pixel size (A) of an MRC file's header: cella x / mx."""
    try:
        with open(path, "rb") as f:
            head = f.read(1024)
    except OSError:
        return None
    end = ">" if len(head) > 212 and head[212] == 0x11 else "<"
    mx = struct.unpack(end + "i", head[28:32])[0]
    cella = struct.unpack(end + "f", head[40:44])[0]
    return cella / mx if mx > 0 and cella > 0 else None


def thumbs_2d(groups: list[dict], out: list[dict[str, str]]) -> None:
    """Each class's average from the project's cavg stack (class n: image n)."""
    stk = next((e["stk"] for e in out if e.get("imgkind") == "cavg" and "stk" in e), None)
    if stk is None:
        print("note: the project registers no class averages (cavg): no thumbnails", file=sys.stderr)
        return
    for g in groups:
        g["thumb"] = read_mrc_image(stk, g["key"] - 1)
    if any(g["thumb"] is None for g in groups):
        print(f"note: some classes have no image in {stk}", file=sys.stderr)


def thumbs_3d(args: argparse.Namespace, groups: list[dict], eulers: list, normals: list,
              out: list[dict[str, str]]) -> None:
    """Each direction's projection of its state's volume (SIMPLE's reproject), along the direction
    of the group's particle closest to the group's mean viewing direction, in-plane angle 0."""
    if args.vols is not None:
        vols = {s: str(p.resolve()) for s, p in enumerate(args.vols, start=1)}
        smpds = {s: mrc_pixel_size(p) for s, p in vols.items()}
    else:
        vols = {int(float(e["state"])): e["vol"] for e in out if e.get("imgkind") == "vol" and "vol" in e and "state" in e}
        smpds = {int(float(e["state"])): float(e["smpd"]) for e in out
                 if e.get("imgkind") == "vol" and "smpd" in e and "state" in e}
    vols = {s: p for s, p in vols.items() if Path(p).is_file()}
    wanted = [g for g in groups if g["key"][0] in vols]
    if len(wanted) < len(groups):
        missing = sorted({g["key"][0] for g in groups} - set(vols))
        print(f"note: no volume for state(s) {missing}: those directions have no thumbnail "
              "(give the volumes with --vols)", file=sys.stderr)
    if not wanted:
        return
    smpd = next((v for v in smpds.values() if v), None)
    if smpd is None:
        print("note: the volumes' pixel size is unknown: no thumbnails", file=sys.stderr)
        return
    lines = []
    for g in wanted:
        nmean = [sum(normals[i][k] for i in g["members"]) for k in range(3)]
        rep = max(g["members"], key=lambda i: sum(normals[i][k] * nmean[k] for k in range(3)))
        e1, e2 = eulers[rep]
        g["angles"] = (e1, e2)
        lines.append(f"e1={e1:.4f} e2={e2:.4f} e3=0.0 state={g['key'][0]}")
    nstates = max(vols)
    with tempfile.TemporaryDirectory(prefix="simple_shift_thumbs_") as tmp:
        Path(tmp, "oris.txt").write_text("\n".join(lines) + "\n")
        cmd = [simple_tool(args, "simple_exec"), "prg=reproject", "oritab=oris.txt", "outstk=thumbs.mrc",
               f"smpd={smpd}", f"nstates={nstates}", "nspace=2", "pgrp=c1", "nthr=4", "mkdir=no"]
        cmd += [f"vol{s}={p}" for s, p in sorted(vols.items())]
        res = run(cmd, cwd=tmp)
        stk = str(Path(tmp, "thumbs.mrc"))
        if res.returncode != 0 or not Path(stk).is_file():
            print("note: reproject failed: no thumbnails", file=sys.stderr)
            return
        for i, g in enumerate(wanted):
            g["thumb"] = read_mrc_image(stk, i)


# ---------------- statistics ----------------

def percentile(sorted_vals: list[float], q: float) -> float:
    """Linear-interpolated percentile (0-100) of an ascending list."""
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = (len(sorted_vals) - 1) * q / 100.0
    lo = math.floor(pos)
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (pos - lo)


def summarise(vals: list[float]) -> dict[str, float]:
    s = sorted(vals)
    return {
        "mean": statistics.fmean(vals),
        "sd": statistics.pstdev(vals) if len(vals) > 1 else 0.0,
        "median": percentile(s, 50.0),
        "p95": percentile(s, 95.0),
        "min": s[0],
        "max": s[-1],
    }


def stats_lines(xs: list[float], ys: list[float]) -> list[tuple[str, dict[str, float]]]:
    radial = [math.hypot(x, y) for x, y in zip(xs, ys)]
    return [("x", summarise(xs)), ("y", summarise(ys)), ("|shift|", summarise(radial))]


def print_summary(data: dict, stats) -> None:
    print(f"{data['source']}   {len(data['xs'])} particles   shifts in {data['unit']}")
    print(f"{'':>8} {'mean':>9} {'sd':>9} {'median':>9} {'95%':>9} {'min':>9} {'max':>9}")
    for name, st in stats:
        print(f"{name:>8} {st['mean']:9.3f} {st['sd']:9.3f} {st['median']:9.3f} {st['p95']:9.3f} "
              f"{st['min']:9.3f} {st['max']:9.3f}")
    if data["groups"]:
        g0 = data["groups"]
        nthumbs = sum(g["thumb"] is not None for g in g0)
        print(f"groups: {len(g0)} of {g0[0]['ntotal']} plotted, the most populated first "
              f"({len(g0[0]['members'])} to {len(g0[-1]['members'])} particles); {nthumbs} with a thumbnail")


def bin_edges(halfrange: float, width: float) -> list[float]:
    """Bin edges of @width centred on whole multiples of @width, symmetric about 0 and covering
    +/- @halfrange: the outermost bins hold the largest |shift|."""
    k = max(0, math.ceil(halfrange / width - 0.5))
    return [(i - 0.5) * width for i in range(-k, k + 2)]


# ---------------- the report ----------------

def plot_report(out: Path, data: dict, stats, edges: list[float], log2d: bool, field: str | None) -> None:
    cache = Path(tempfile.gettempdir()) / "simple-matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.colors import ListedColormap, LogNorm, Normalize
    except ImportError as err:
        raise SystemExit(f"the PDF needs numpy and matplotlib ({err}); the statistics are above")

    unit = data["unit"]
    x = np.asarray(data["xs"])
    y = np.asarray(data["ys"])
    e = np.asarray(edges)
    lim = float(e[-1])
    width = float(e[1] - e[0])
    nout = int(np.count_nonzero((np.abs(x) > lim) | (np.abs(y) > lim)))
    cmap = ListedColormap(SEQUENTIAL)

    plt.rcParams.update({
        "font.size": 9,
        "text.color": TEXT_PRIMARY,
        "axes.labelcolor": TEXT_SECONDARY,
        "axes.edgecolor": HAIRLINE,
        "axes.linewidth": 0.6,
        "axes.facecolor": SURFACE,
        "figure.facecolor": SURFACE,
        "xtick.color": TEXT_SECONDARY,
        "ytick.color": TEXT_SECONDARY,
        "xtick.major.size": 0,
        "ytick.major.size": 0,
        "grid.color": HAIRLINE,
        "grid.linewidth": 0.5,
        "grid.linestyle": "-",
    })

    def histogram(ax, vals, label, ymax, ylabel=True, small=False):
        ax.hist(vals, bins=e, color=SERIES, rwidth=0.82, zorder=2)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(0, ymax if ymax > 0 else 1)
        ax.set_xlabel(f"{label} shift ({unit})", fontsize=7 if small else None)
        if ylabel:
            ax.set_ylabel("particles", fontsize=7 if small else None)
        ax.grid(axis="y", zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        if small:
            ax.tick_params(labelsize=6)

    def histogram2d(ax, xv, yv, small=False):
        counts, xe, ye = np.histogram2d(xv, yv, bins=[e, e])
        vmax = max(counts.max(), 1)
        norm = LogNorm(vmin=1, vmax=vmax) if log2d else Normalize(vmin=0, vmax=vmax)
        mesh = ax.pcolormesh(xe, ye, np.ma.masked_equal(counts.T, 0), cmap=cmap, norm=norm, shading="flat")
        ax.set_aspect("equal")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.axhline(0, color=HAIRLINE, linewidth=0.5, zorder=0)
        ax.axvline(0, color=HAIRLINE, linewidth=0.5, zorder=0)
        ax.set_xlabel(f"x shift ({unit})", fontsize=7 if small else None)
        ax.set_ylabel(f"y shift ({unit})", fontsize=7 if small else None)
        ax.spines[["top", "right"]].set_visible(False)
        if small:
            ax.tick_params(labelsize=6)
        return mesh

    pages = []
    # ---- page 1: the source, the statistics and the three plots of every particle
    fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait
    grid = fig.add_gridspec(3, 2, height_ratios=[0.55, 1.0, 1.55], hspace=0.42, wspace=0.28,
                            left=0.09, right=0.93, top=0.95, bottom=0.05)
    head = fig.add_subplot(grid[0, :])
    head.axis("off")
    head.text(0.0, 1.0, "Particle origin shifts", fontsize=14, weight="bold", va="top")
    detail = f"{len(x):,} particles   shifts in {unit}   bins {width:g} {unit}"
    if nout > 0:
        detail += f"   {nout:,} outside ±{lim:g} {unit} not plotted"
    head.text(0.0, 0.84, data["source"], color=TEXT_SECONDARY, va="top")
    head.text(0.0, 0.72, detail, color=TEXT_SECONDARY, va="top")
    cols = ["", "mean", "sd", "median", "95%", "min", "max"]
    colx = [0.0, 0.16, 0.29, 0.42, 0.55, 0.68, 0.81]
    rowy = 0.50
    for cx, c in zip(colx, cols):
        head.text(cx, rowy, c, color=TEXT_SECONDARY, va="top", ha="left", family="monospace")
    for name, st in stats:
        rowy -= 0.15
        vals = [name] + [f"{st[k]:8.3f}" for k in ("mean", "sd", "median", "p95", "min", "max")]
        for cx, v in zip(colx, vals):
            head.text(cx, rowy, v, va="top", ha="left", family="monospace")
    counts_x, _ = np.histogram(x, bins=e)
    counts_y, _ = np.histogram(y, bins=e)
    ymax = max(counts_x.max(), counts_y.max()) * 1.05
    for col, (vals, label) in enumerate(((x, "x"), (y, "y"))):
        ax = fig.add_subplot(grid[1, col])
        histogram(ax, vals, label, ymax)
        ax.set_title(f"{label} shift", loc="left", fontsize=10, color=TEXT_PRIMARY)
    ax = fig.add_subplot(grid[2, :])
    mesh = histogram2d(ax, x, y)
    ax.set_title("x against y shift", loc="left", fontsize=10, color=TEXT_PRIMARY)
    bar = fig.colorbar(mesh, ax=ax, fraction=0.04, pad=0.03)
    bar.set_label("particles per bin" + (" (log scale)" if log2d else ""), color=TEXT_SECONDARY)
    bar.outline.set_visible(False)
    bar.ax.tick_params(colors=TEXT_SECONDARY, size=0)
    pages.append(fig)

    # ---- the group pages: one row per class/direction, led by its thumbnail
    groups = data["groups"]
    if groups:
        kind = "2D class" if field == "ptcl2D" else "projection direction"
        thumb = "class average" if field == "ptcl2D" else "projection of the state's volume"
        npages = math.ceil(len(groups) / ROWS_PER_PAGE)
        for p in range(npages):
            rows = groups[p * ROWS_PER_PAGE:(p + 1) * ROWS_PER_PAGE]
            fig = plt.figure(figsize=(8.27, 11.69))
            grid = fig.add_gridspec(ROWS_PER_PAGE, 4, width_ratios=[1.0, 1.25, 1.25, 1.0], hspace=0.62,
                                    wspace=0.42, left=0.06, right=0.97, top=0.91, bottom=0.05)
            fig.text(0.06, 0.975, f"Shifts per {kind}", fontsize=13, weight="bold", va="top")
            fig.text(0.06, 0.952, f"the {len(groups)} most populated of {groups[0]['ntotal']}, page {p + 1} of {npages}; "
                     f"each row led by its {thumb}", color=TEXT_SECONDARY, va="top")
            fig.text(0.06, 0.935, f"shift axes as on page 1 (±{lim:g} {unit}, {width:g} {unit} bins); "
                     "each row's counts and colours scaled to its own group", color=TEXT_SECONDARY, va="top")
            for r, g in enumerate(rows):
                idx = np.asarray(g["members"])
                gx, gy = x[idx], y[idx]
                ax = fig.add_subplot(grid[r, 0])
                ax.set_xticks([])
                ax.set_yticks([])
                for side in ax.spines.values():
                    side.set_visible(False)
                if g["thumb"] is not None:
                    nx, ny, vals = g["thumb"]
                    img = np.asarray(vals, dtype=float).reshape(ny, nx)
                    lo, hi = np.percentile(img, [1, 99])
                    ax.imshow(img, cmap="gray", origin="lower", vmin=lo, vmax=hi if hi > lo else lo + 1)
                else:
                    ax.set_facecolor(HAIRLINE)
                    ax.text(0.5, 0.5, "no image", ha="center", va="center", color=TEXT_SECONDARY, fontsize=7,
                            transform=ax.transAxes)
                ax.set_title(g["label"], loc="left", fontsize=9, color=TEXT_PRIMARY)
                sub = f"{len(idx):,} particles"
                if "angles" in g:
                    sub += f"\ne1 {g['angles'][0]:.0f}°, e2 {g['angles'][1]:.0f}°"
                ax.set_xlabel(sub, fontsize=7, color=TEXT_SECONDARY)
                cx, _ = np.histogram(gx, bins=e)
                cy, _ = np.histogram(gy, bins=e)
                gmax = max(cx.max(), cy.max()) * 1.05
                for col, (vals, label) in enumerate(((gx, "x"), (gy, "y")), start=1):
                    ax = fig.add_subplot(grid[r, col])
                    histogram(ax, vals, label, gmax, ylabel=(col == 1), small=True)
                    ax.set_title(f"mean {vals.mean():.2f}, sd {vals.std():.2f}", loc="left", fontsize=7,
                                 color=TEXT_SECONDARY)
                ax = fig.add_subplot(grid[r, 3])
                histogram2d(ax, gx, gy, small=True)
            pages.append(fig)

    with PdfPages(out) as pdf:
        for fig in pages:
            pdf.savefig(fig)
            plt.close(fig)


def main() -> int:
    args = parse_args()
    if not args.input.is_file():
        raise SystemExit(f"no such file: {args.input}")
    if args.binwidth <= 0.0:
        raise SystemExit("--binwidth must be positive")
    if args.limit is not None and args.limit <= 0.0:
        raise SystemExit("--limit must be positive")
    if args.max_groups < 0:
        raise SystemExit("--max-groups must be 0 or more")
    is_project = args.input.suffix == ".simple"
    data = load_project(args) if is_project else load_star(args)
    if not data["xs"]:
        raise SystemExit("no particles to plot")
    stats = stats_lines(data["xs"], data["ys"])
    print_summary(data, stats)
    xs, ys = data["xs"], data["ys"]
    halfrange = args.limit if args.limit is not None else max(max(abs(v) for v in xs), max(abs(v) for v in ys))
    edges = bin_edges(halfrange, args.binwidth)
    tag = f"_{args.field}" if is_project else ""
    out = args.output if args.output is not None else args.input.with_name(f"{args.input.stem}{tag}_shifts.pdf")
    plot_report(out, data, stats, edges, args.log2d, args.field)
    print(f"written: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
