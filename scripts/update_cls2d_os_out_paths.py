#!/usr/bin/env python3
"""Update os_out cls2D/cavg stack paths for SIMPLE project files.

For each discovered *.simple project, this script points the os_out cavg
(entry identified by imgkind, default: "cavg") to a sibling .mrc file that
shares the same stem as the project file.

Implementation strategy:
1) Query os_out via:
    simple_private_exec prg=print_project_vals oritype=out keys=imgkind,smpd
2) Extract smpd for the target imgkind.
3) Apply update via:
     simple_exec prg=import_cavgs projfile=<proj> stk=<proj_stem>.mrc smpd=<smpd>

Notes:
- import_cavgs updates (or creates) the os_out entry for the imgkind.
- The target .mrc file must exist (simple_exec validates dimensions).
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable, Optional


def find_simple_files(root: Path) -> Iterable[Path]:
    for dirpath, _, filenames in os.walk(root):
        for name in filenames:
            if name.endswith(".simple"):
                yield Path(dirpath) / name


def run_cmd(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, text=True, capture_output=True, check=False)


def parse_smpd(print_vals_stdout: str, imgkind: str) -> Optional[float]:
    """Parse smpd from print_project_vals output.

    Expected data rows from SIMPLE look like:
      <index> <state> <imgkind> <smpd>
    We ignore non-data lines and logs.
    """
    print(print_vals_stdout)
    for raw in print_vals_stdout.splitlines():
        line = raw.strip()
        if not line:
            continue
        cols = line.split()
        if len(cols) < 4:
            continue
        # Data rows begin with index and state integers.
        if not (cols[0].isdigit() and cols[1].lstrip("+-").isdigit()):
            continue
        if cols[2] != imgkind:
            continue
        try:
            return float(cols[3])
        except ValueError:
            continue
    return None


def parse_ncls(print_vals_stdout: str) -> int:
    """Count cls2D rows from print_project_vals output."""
    ncls = 0
    for raw in print_vals_stdout.splitlines():
        line = raw.strip()
        if not line:
            continue
        cols = line.split()
        if len(cols) < 3:
            continue
        if not (cols[0].isdigit() and cols[1].lstrip("+-").isdigit()):
            continue
        ncls += 1
    return ncls


def parse_selection_file(selection_file: Path, ncls: int) -> tuple[Optional[list[int]], str]:
    """Parse sibling .txt into a cls2D state vector of length ncls.

    Supported formats:
    1) 0/1 states with exactly ncls entries.
    2) 1-based selected class indices (any length).
    """
    tokens: list[str] = []
    for raw in selection_file.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        tokens.extend(line.replace(",", " ").split())

    if not tokens:
        return None, f"skip selection (empty selection file): {selection_file}"

    values: list[int] = []
    for tok in tokens:
        try:
            values.append(int(tok))
        except ValueError:
            return None, f"skip selection (non-integer token '{tok}'): {selection_file}"

    if len(values) == ncls and all(v in (0, 1) for v in values):
        return values, ""

    states = [0] * ncls
    for idx in values:
        if idx < 1 or idx > ncls:
            return None, f"skip selection (index out of range 1..{ncls}: {idx}): {selection_file}"
        states[idx - 1] = 1
    return states, ""


def apply_cls2d_selection(
    projfile: Path,
    simple_exec: str,
    states: list[int],
    dry_run: bool,
) -> tuple[bool, str]:
    if dry_run:
        return True, f"dry-run: apply cls2D selection to {projfile} ({sum(states)} selected)"

    tmp_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False, encoding="utf-8") as tf:
            tmp_path = Path(tf.name)
            for s in states:
                tf.write(f"{s}\n")

        cmd = [
            simple_exec,
            "prg=selection",
            f"projfile={projfile}",
            "oritype=cls2D",
            f"infile={tmp_path}",
            "prune=no",
            "mkdir=no",
        ]
        out = run_cmd(cmd)
        if out.returncode != 0:
            return False, f"selection failed: {projfile}\n{out.stderr.strip()}"
        return True, f"selection applied: {projfile} ({sum(states)} selected)"
    finally:
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink()


def get_ncls(projfile: Path, simple_exec: str) -> tuple[int, str]:
    cmd = [
        simple_exec,
        "prg=print_project_vals",
        f"projfile={projfile}",
        "oritype=cls2D",
        "keys=state",
    ]
    out = run_cmd(cmd)
    if out.returncode != 0:
        return 0, f"skip selection (print_project_vals cls2D failed): {projfile}\n{out.stderr.strip()}"
    ncls = parse_ncls(out.stdout)
    if ncls <= 0:
        return 0, f"skip selection (no cls2D rows found): {projfile}"
    return ncls, ""


def update_one_project(
    projfile: Path,
    simple_exec: str,
    simple_private_exec: str,
    imgkind: str,
    dry_run: bool,
) -> tuple[bool, str]:
    messages: list[str] = []
    target_stk = projfile.with_suffix(".mrc")
    if not target_stk.exists():
        return False, f"skip (missing target mrc): {target_stk}"

    q_cmd = [
        simple_private_exec,
        "prg=print_project_vals",
        f"projfile={projfile}",
        "oritype=out",
        "keys=imgkind,smpd",
    ]
    q = run_cmd(q_cmd)
    if q.returncode != 0:
        return False, f"skip (print_project_vals failed): {projfile}\n{q.stderr.strip()}"

    smpd = parse_smpd(q.stdout, imgkind)
    if smpd is None:
        return False, f"skip (imgkind={imgkind} with smpd not found): {projfile}"

    u_cmd = [
        simple_exec,
        "prg=import_cavgs",
        f"projfile={projfile}",
        f"stk={target_stk}",
        f"smpd={smpd}",
        "mkdir=no",
    ]

    if dry_run:
        messages.append(f"dry-run: {' '.join(u_cmd)}")
    else:
        u = run_cmd(u_cmd)
        if u.returncode != 0:
            return False, f"failed: {projfile}\n{u.stderr.strip()}"
        messages.append(f"updated: {projfile} -> {target_stk}")

    selection_file = projfile.with_suffix(".txt")
    if selection_file.exists():
        ncls, err = get_ncls(projfile, simple_private_exec)
        if ncls <= 0:
            messages.append(err)
            return False, "\n".join(messages)

        states, perr = parse_selection_file(selection_file, ncls)
        if states is None:
            messages.append(perr)
            return False, "\n".join(messages)

        sel_ok, sel_msg = apply_cls2d_selection(projfile, simple_exec, states, dry_run)
        messages.append(sel_msg)
        if not sel_ok:
            return False, "\n".join(messages)
    else:
        messages.append(f"selection file not found (skip): {selection_file}")

    return True, "\n".join(messages)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Find *.simple files and update os_out cls2D/cavg stack path "
            "to sibling .mrc (same stem), and apply cls2D selection from "
            "sibling .txt when present."
        )
    )
    parser.add_argument(
        "root",
        nargs="?",
        default=".",
        help="Root directory to scan recursively (default: current directory).",
    )
    parser.add_argument(
        "--simple-exec",
        default="simple_exec",
        help="Path to simple_exec binary (default: simple_exec).",
    )
    parser.add_argument(
        "--simple-private-exec",
        default="simple_private_exec",
        help="Path to simple_private_exec binary for print_project_vals (default: simple_private_exec).",
    )
    parser.add_argument(
        "--imgkind",
        default="cavg",
        help="os_out imgkind to update (default: cavg).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print actions without modifying projects.",
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    if not root.exists():
        print(f"error: root does not exist: {root}", file=sys.stderr)
        return 2

    files = list(find_simple_files(root))
    if not files:
        print(f"no .simple files found under {root}")
        return 0

    ok = 0
    fail = 0
    for proj in files:
        success, msg = update_one_project(
            projfile=proj,
            simple_exec=args.simple_exec,
            simple_private_exec=args.simple_private_exec,
            imgkind=args.imgkind,
            dry_run=args.dry_run,
        )
        print(msg)
        if success:
            ok += 1
        else:
            fail += 1

    print(f"summary: updated={ok} skipped_or_failed={fail} total={len(files)}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
