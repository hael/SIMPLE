#!/usr/bin/env python3
"""Fix mechanical gfortran warnings from a build log.

Adapted from SIMPLE's scripts/delete_unused_variables.pl. It handles unused
locals, genuinely unused private module variables, and stale `use ..., only:`
names. It runs dry by default and never compiles anything itself -- the owner
rebuilds to verify, which is also how the result is checked.

Files containing conditional-preprocessor branches are deliberately report-only:
a symbol unused in the current build can be required by another configuration.

Usage:
    ./compile_debug.sh 2>&1 | tee build.log
    python3 scripts/fix_warnings.py build.log            # report what would change
    python3 scripts/fix_warnings.py build.log --apply    # edit the files
    python3 scripts/fix_warnings.py build.log --show-manual
    (rebuild; the log should come back clean)

Handled classes:
  - "Unused variable 'x'"                        -> removed from its declaration
  - "Unused PRIVATE module variable 'x'"         -> removed unless submodules use the parent
  - "Unused parameter 'x' ... explicitly imported" -> removed from the use-only list
Everything else is summarized by warning flag; use --show-manual for details.
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

LOC_RE          = re.compile(r"^(.+\.[Ff](?:90|08|95|03)):(\d+):(\d+):")
UNUSED_VAR_RE   = re.compile(r"Warning: Unused variable '(\w+)'")
UNUSED_MOD_RE   = re.compile(r"Warning: Unused PRIVATE module variable '(\w+)' declared")
UNUSED_IMP_RE   = re.compile(r"Warning: Unused parameter '(\w+)' which has been explicitly imported")
OTHER_WARN_RE   = re.compile(r"Warning: (.+)$")
WARN_FLAG_RE    = re.compile(r"\[(-W[^]]+)\]$")
CPP_COND_RE     = re.compile(r"^\s*#\s*(?:if|ifdef|ifndef|elif|else|endif)\b", re.IGNORECASE)
MODULE_RE       = re.compile(
    r"^\s*module\s+(?!procedure\b|subroutine\b|function\b)(\w+)",
    re.IGNORECASE | re.MULTILINE,
)
SUBMODULE_RE    = re.compile(r"^\s*submodule\s*\(\s*([^)]+)\)", re.IGNORECASE | re.MULTILINE)


def discover_submodule_ancestors(repo_root: Path) -> set[str]:
    """Return module names whose state can be consumed by a submodule."""
    ancestors: set[str] = set()
    for path in (repo_root / "src").rglob("*"):
        if not path.is_file() or path.suffix.lower() not in {".f90", ".f08", ".f95", ".f03"}:
            continue
        if "build" in path.relative_to(repo_root).parts:
            continue
        text = path.read_text(errors="replace")
        for match in SUBMODULE_RE.finditer(text):
            ancestors.update(part.strip().lower() for part in match.group(1).split(":") if part.strip())
    return ancestors


def declares_submodule_parent(lines: list[str], submodule_ancestors: set[str]) -> bool:
    """Whether this source defines a module named as a submodule ancestor."""
    match = MODULE_RE.search("".join(lines))
    return bool(match and match.group(1).lower() in submodule_ancestors)


def split_decl_items(text: str) -> list[str]:
    """Split a declaration's entity list on top-level commas (paren-aware)."""
    items, cur, depth = [], "", 0
    for ch in text:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth = max(0, depth - 1)
        if ch == "," and depth == 0:
            items.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        items.append(cur.strip())
    return items


def entity_name(item: str) -> str:
    """The bare name of one declared entity: strip initializer and dimensions."""
    name = re.split(r"=", item, maxsplit=1)[0]
    name = re.split(r"\(", name, maxsplit=1)[0]
    return name.strip().lower()


def remove_from_declaration(line: str, var: str) -> str | None:
    """Remove `var` from a declaration line; None means delete the line."""
    if "::" not in line:
        return line
    lhs, rhs = line.split("::", 1)
    comment = ""
    m = re.search(r"(\s*!.*)$", rhs)
    if m:
        comment = m.group(1)
        rhs = rhs[: m.start()]
    kept = [it for it in split_decl_items(rhs) if entity_name(it) != var.lower()]
    if not kept:
        return None
    return lhs + ":: " + ", ".join(kept) + comment.rstrip() + "\n"


def statement_extent(lines: list[str], start: int) -> int:
    """Last line index (inclusive) of the free-form statement starting at start."""
    i = start
    while True:
        code = re.sub(r"!.*$", "", lines[i]).rstrip()
        if code.endswith("&") and i + 1 < len(lines):
            i += 1
        else:
            return i


def remove_from_use_only(lines: list[str], start: int, name: str) -> list[str] | None:
    """Rewrite the use-statement spanning lines[start..] without `name`.

    Returns the replacement lines, or None if the name was not found."""
    end = statement_extent(lines, start)
    joined = ""
    for i in range(start, end + 1):
        code = re.sub(r"!.*$", "", lines[i]).strip()
        code = code.rstrip("&").strip()
        code = code.lstrip("&").strip()
        joined += (" " if joined else "") + code
    m = re.match(r"(use\s+\w+\s*,\s*only\s*:)(.*)$", joined, re.IGNORECASE)
    if not m:
        return None
    head, taillist = m.group(1), m.group(2)
    names = [n.strip() for n in taillist.split(",") if n.strip()]
    kept = [n for n in names if n.lower() != name.lower()]
    if len(kept) == len(names):
        return None
    indent = re.match(r"^(\s*)", lines[start]).group(1)
    if not kept:
        return []                                   # the whole import is stale
    return [f"{indent}{head.strip().replace('only :', 'only:')} " + ", ".join(kept) + "\n"]


def main() -> int:
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    log_path = Path(sys.argv[1])
    apply = "--apply" in sys.argv[2:]
    show_manual = "--show-manual" in sys.argv[2:]
    repo_root = Path(__file__).resolve().parent.parent
    submodule_ancestors = discover_submodule_ancestors(repo_root)

    # collect (file, line, kind, name) from the log
    edits: list[tuple[str, int, str, str]] = []
    manual: list[tuple[str, str]] = []
    loc: tuple[str, int] | None = None
    for raw in log_path.read_text(errors="replace").splitlines():
        m = LOC_RE.match(raw)
        if m:
            loc = (m.group(1), int(m.group(2)))
            continue
        if loc is None:
            continue
        if m := UNUSED_VAR_RE.search(raw):
            edits.append((loc[0], loc[1], "var", m.group(1)))
        elif m := UNUSED_MOD_RE.search(raw):
            edits.append((loc[0], loc[1], "module_var", m.group(1)))
        elif m := UNUSED_IMP_RE.search(raw):
            edits.append((loc[0], loc[1], "import", m.group(1)))
        elif m := OTHER_WARN_RE.search(raw):
            flag_match = WARN_FLAG_RE.search(raw)
            flag = flag_match.group(1) if flag_match else "unclassified"
            manual.append((flag, f"{loc[0]}:{loc[1]}: {m.group(1)}"))

    # apply per file, bottom-up so line numbers stay valid
    by_file: dict[str, list[tuple[int, str, str]]] = {}
    for f, ln, kind, name in dict.fromkeys(edits):
        by_file.setdefault(f, []).append((ln, kind, name))
    nchanges = 0
    nskipped = 0
    for f, items in sorted(by_file.items()):
        path = Path(f)
        if not path.is_file():
            print(f"SKIP (missing): {f}")
            nskipped += len(items)
            continue
        lines = path.read_text().splitlines(keepends=True)
        if any(CPP_COND_RE.match(line) for line in lines):
            print(f"SKIP (conditional preprocessing; manual review required): {f}")
            nskipped += len(items)
            continue
        is_submodule_parent = declares_submodule_parent(lines, submodule_ancestors)
        for ln, kind, name in sorted(items, reverse=True):
            if kind == "module_var" and is_submodule_parent:
                print(f"SKIP (parent-module state; may be used by a submodule): {f}:{ln} ({name})")
                nskipped += 1
                continue
            idx = ln - 1
            if idx >= len(lines):
                print(f"SKIP (line out of range): {f}:{ln}")
                nskipped += 1
                continue
            if kind in {"var", "module_var"}:
                new = remove_from_declaration(lines[idx], name)
                old = lines[idx].rstrip()
                if new is None:
                    print(f"{f}:{ln}: delete declaration: {old.strip()}")
                    del lines[idx]
                elif new != lines[idx]:
                    print(f"{f}:{ln}: {old.strip()}  ->  {new.strip()}")
                    lines[idx] = new
                else:
                    print(f"SKIP (no match for {name!r}): {f}:{ln}")
                    nskipped += 1
                    continue
            else:
                repl = remove_from_use_only(lines, idx, name)
                if repl is None:
                    print(f"SKIP (import {name!r} not found): {f}:{ln}")
                    nskipped += 1
                    continue
                end = statement_extent(lines, idx)
                print(f"{f}:{ln}: drop import {name!r}" + ("" if repl else " (whole use statement)"))
                lines[idx : end + 1] = repl
            nchanges += 1
        if apply:
            path.write_text("".join(lines))

    print(f"\n{nchanges} mechanical fix(es) " + ("applied" if apply else "previewed (rerun with --apply)"))
    if nskipped:
        print(f"{nskipped} recognized warning(s) left unchanged by safety/staleness checks")
    if manual:
        counts = Counter(flag for flag, _ in manual)
        print(f"{len(manual)} warning(s) need human judgment:")
        for flag, count in counts.most_common():
            print(f"  {count:5d}  {flag}")
        if show_manual:
            print("Manual warning details:")
            for _, message in manual:
                print("  " + message)
    print("Rebuild to verify; the log should come back clean.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
