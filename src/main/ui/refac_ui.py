#!/usr/bin/env python3
"""
refactor_ui_api_modules.py

Refactors SIMPLE Fortran UI API modules of the form:

module simple_ui_api_<cat>
use simple_ui_program, only: ui_program
implicit none
type(ui_program), target :: foo
type(ui_program), target :: bar_
end module simple_ui_api_<cat>

into the concrete "register" pattern:

module simple_ui_api_<cat>
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none
public :: register_ui_<cat>

type(ui_program), target :: foo
type(ui_program), target :: bar_

contains
  subroutine register_ui_<cat>(prgtab)
    class(ui_hash), intent(inout) :: prgtab
    call add_ui_program('foo', foo, prgtab)
    call add_ui_program('bar', bar_, prgtab)
  end subroutine
end module

Notes:
- By default, key names are variable names with a trailing "_" stripped (mkdir_ -> "mkdir").
- Modules that already contain "contains" are skipped unless --force is passed.
- Creates .bak backups before writing.

Usage:
  python3 refactor_ui_api_modules.py /path/to/simple_ui_api_*.f90
  python3 refactor_ui_api_modules.py -r /path/to/src --glob "simple_ui_api_*.f90"
  python3 refactor_ui_api_modules.py -r /path/to/src --glob "simple_ui_api_*.f90" --force
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Tuple


RE_MOD = re.compile(r'^\s*module\s+([a-zA-Z_]\w*)\s*$', re.IGNORECASE | re.MULTILINE)
RE_ENDMOD = re.compile(r'^\s*end\s+module\s+([a-zA-Z_]\w*)\s*$', re.IGNORECASE | re.MULTILINE)
RE_USE = re.compile(r'^\s*use\s+([a-zA-Z_]\w*)\b(.*)$', re.IGNORECASE | re.MULTILINE)
RE_IMPLICIT = re.compile(r'^\s*implicit\s+none\b', re.IGNORECASE | re.MULTILINE)
RE_CONTAINS = re.compile(r'^\s*contains\b', re.IGNORECASE | re.MULTILINE)

# Matches lines like: type(ui_program), target :: a, b, c_
RE_UI_DECL = re.compile(
    r'^\s*type\s*\(\s*ui_program\s*\)\s*,\s*target\s*::\s*(.+?)\s*$',
    re.IGNORECASE | re.MULTILINE
)


def strip_inline_comment(line: str) -> str:
    # Fortran inline comments start with !
    # Keep everything before !
    return line.split('!')[0].rstrip()


def split_decl_list(rhs: str) -> List[str]:
    # rhs may contain comma separated vars; remove comments already
    parts = [p.strip() for p in rhs.split(',')]
    return [p for p in parts if p]


def derive_suffix(modname: str) -> str:
    # expected: simple_ui_api_<suffix>
    lower = modname.lower()
    prefix = "simple_ui_api_"
    if lower.startswith(prefix):
        return modname[len(prefix):]  # preserve original case
    # fallback
    return modname


def key_from_var(var: str) -> str:
    # strip trailing underscore ONLY if it's a final underscore
    return var[:-1] if var.endswith('_') else var


def normalize_use_lines(existing: str) -> Tuple[List[str], bool]:
    """
    Return a minimal list of use lines:
      - keep 'use simple_ui_program, only: ui_program' if present; otherwise add it.
      - ensure 'use simple_ui_hash, only: ui_hash'
      - ensure 'use simple_ui_utils, only: add_ui_program'
    Also return whether we found any use lines at all (for formatting).
    """
    uses = {}
    found_any = False

    for m in RE_USE.finditer(existing):
        found_any = True
        mod = m.group(1)
        rest = m.group(2).rstrip()
        # store last occurrence
        uses[mod.lower()] = f"use {mod}{rest}"

    def set_use(mod: str, line: str):
        uses.setdefault(mod.lower(), line)

    # Ensure required use lines exist; prefer canonical formatting
    set_use("simple_ui_program", "use simple_ui_program, only: ui_program")
    set_use("simple_ui_hash",    "use simple_ui_hash,    only: ui_hash")
    set_use("simple_ui_utils",   "use simple_ui_utils,   only: add_ui_program")

    # Output in a stable order: program, hash, utils, then any other existing uses
    out = []
    for m in ["simple_ui_program", "simple_ui_hash", "simple_ui_utils"]:
        out.append(uses[m])

    # Append any other existing uses (preserve original spelling as stored)
    for k, v in uses.items():
        if k in ("simple_ui_program", "simple_ui_hash", "simple_ui_utils"):
            continue
        out.append(v)

    return out, found_any


def extract_ui_program_vars(text: str) -> List[str]:
    vars_: List[str] = []
    for m in RE_UI_DECL.finditer(text):
        rhs = strip_inline_comment(m.group(1))
        vars_.extend(split_decl_list(rhs))
    # de-dup while preserving order
    seen = set()
    out = []
    for v in vars_:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def build_refactored_module(modname: str, vars_: List[str]) -> str:
    suffix = derive_suffix(modname)
    reg = f"register_ui_{suffix}"

    use_lines, _ = normalize_use_lines("")  # canonical only
    lines: List[str] = []
    lines.append(f"!@descr: \"{suffix}\" UI api (concrete implementation)")
    lines.append(f"module {modname}")
    lines.extend(use_lines)
    lines.append("implicit none")
    lines.append(f"public :: {reg}")
    lines.append("")
    for v in vars_:
        lines.append(f"type(ui_program), target :: {v}")
    lines.append("")
    lines.append("contains")
    lines.append("")
    lines.append(f"    subroutine {reg}(prgtab)")
    lines.append("        class(ui_hash), intent(inout) :: prgtab")
    for v in vars_:
        key = key_from_var(v)
        lines.append(f"        call add_ui_program('{key}', {v}, prgtab)")
    lines.append(f"    end subroutine {reg}")
    lines.append("")
    lines.append(f"end module {modname}")
    lines.append("")
    return "\n".join(lines)


def can_safely_rewrite(text: str, force: bool) -> Tuple[bool, str]:
    if not RE_MOD.search(text):
        return False, "no MODULE statement found"
    if not RE_ENDMOD.search(text):
        return False, "no END MODULE statement found"
    if RE_CONTAINS.search(text) and not force:
        return False, "already has CONTAINS (use --force to overwrite)"
    if not RE_IMPLICIT.search(text):
        # not a dealbreaker, but indicates non-standard structure
        if not force:
            return False, "no IMPLICIT NONE found (use --force to overwrite)"
    if not extract_ui_program_vars(text):
        return False, "no type(ui_program), target :: ... declarations found"
    return True, ""


def process_file(path: Path, force: bool, dry_run: bool) -> Tuple[bool, str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    ok, reason = can_safely_rewrite(text, force=force)
    if not ok:
        return False, reason

    m = RE_MOD.search(text)
    assert m
    modname = m.group(1)

    vars_ = extract_ui_program_vars(text)
    new_text = build_refactored_module(modname, vars_)

    if dry_run:
        return True, f"DRY RUN: would rewrite ({len(vars_)} programs)"
    # backup then write
    bak = path.with_suffix(path.suffix + ".bak")
    bak.write_text(text, encoding="utf-8")
    path.write_text(new_text, encoding="utf-8")
    return True, f"rewrote ({len(vars_)} programs), backup -> {bak.name}"


def discover_files(root: Path, glob_pat: str) -> List[Path]:
    return sorted(root.rglob(glob_pat))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="*", help="Fortran files to rewrite (e.g. simple_ui_api_*.f90)")
    ap.add_argument("-r", "--root", type=str, default=None, help="Root directory to search recursively")
    ap.add_argument("--glob", type=str, default="simple_ui_api_*.f90", help="Glob pattern for recursive search")
    ap.add_argument("--force", action="store_true", help="Overwrite modules that already have CONTAINS / nonstandard layout")
    ap.add_argument("--dry-run", action="store_true", help="Do not write files; just report what would change")
    args = ap.parse_args()

    paths: List[Path] = []
    if args.root:
        paths.extend(discover_files(Path(args.root), args.glob))
    paths.extend(Path(f) for f in args.files)

    # de-dup
    uniq = []
    seen = set()
    for p in paths:
        p = p.resolve()
        if p not in seen and p.exists():
            seen.add(p)
            uniq.append(p)

    if not uniq:
        print("No files found.")
        return

    changed = 0
    skipped = 0
    for p in uniq:
        ok, msg = process_file(p, force=args.force, dry_run=args.dry_run)
        if ok:
            changed += 1
            print(f"[OK]   {p}: {msg}")
        else:
            skipped += 1
            print(f"[SKIP] {p}: {msg}")

    print(f"\nSummary: ok={changed} skipped={skipped} total={changed+skipped}")


if __name__ == "__main__":
    main()
