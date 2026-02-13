#!/usr/bin/env python3
"""
Format-preserving refactor for legacy ui_program constructors.

Transforms:

1) call prg%new( name, short, long, exec, nimg,nparm,nalt,nsrch,nfilt,nmask,ncomp, sp_required, [kw...] )
   -> removes the 7 count integers only, leaving existing line breaks intact:
      call prg%new( name, short, long, exec, sp_required, [kw...] )

   NOTE: We DO NOT reflow lines. We surgically delete the "counts" segment including its commas.

2) call prg%set_input('parm_ios', i, ...)
   -> call prg%add_input(UI_PARM, ...)
   (same for other sections; drops the explicit index argument)
   NOTE: Preserves original multiline formatting and continuations.

3) Post-hoc list field edits (which no longer compile after switching to linked_list):
     prg%filt_ctrls(6)%descr_long = '...'
   are rewritten to:

     block
         type(ui_param) :: ptmp
         class(*), allocatable :: tmp
         call prg%filt_ctrls%at(6, tmp)
         select type(t => tmp)
         type is (ui_param)
             ptmp = t
         class default
             THROW_HARD('... not ui_param')
         end select
         if (allocated(tmp)) deallocate(tmp)
         ptmp%descr_long = '...'
         call prg%filt_ctrls%replace_at(6, ptmp)
     end block

   NOTE: The RHS is preserved verbatim (including quotes, etc.). Inline trailing comments are preserved.

Usage:
  python refactor_ui_programs.py path/to/src --write
  python refactor_ui_programs.py file.f90

Default is dry-run with unified diff output.

Assumptions (as you confirmed):
- Post-hoc edits are of the form: <obj>%<section>(<int>)%<field> = <rhs>
- <section> is one of img_ios/parm_ios/alt_ios/srch_ctrls/filt_ctrls/mask_ctrls/comp_ctrls
"""

from __future__ import annotations
import argparse
import difflib
import glob
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Optional


SECTION_UI = {
    "img_ios": "UI_IMG",
    "parm_ios": "UI_PARM",
    "alt_ios": "UI_ALT",
    "srch_ctrls": "UI_SRCH",
    "filt_ctrls": "UI_FILT",
    "mask_ctrls": "UI_MASK",
    "comp_ctrls": "UI_COMP",
}

SECTION_NAMES = "|".join(map(re.escape, SECTION_UI.keys()))

CALL_START_RE = re.compile(r"^(\s*)call\b", re.IGNORECASE)

# Match a set_input(...) where first arg is a quoted section name and second arg is an integer index
# We will ONLY replace the "set_input(" token and the leading "'section', idx," segment.
SET_INPUT_HEAD_RE = re.compile(
    rf"""(?is)
    (?P<prefix>\b%)(?P<meth>set_input)\s*\(\s*
    (?P<q>['"])(?P<section>{SECTION_NAMES})(?P=q)\s*,\s*
    (?P<idx>[+-]?\d+)\s*,\s*
    """,
    re.VERBOSE,
)

# Used to locate %new( ... ) within a call statement text
NEW_CALL_RE = re.compile(r"(?is)\b(?P<obj>[A-Za-z_]\w*(?:%\w+)*)%new\s*\(")

# Post-hoc edit line (single line assignment)
POSTHOC_RE = re.compile(
    rf"""^
    (?P<indent>\s*)
    (?P<obj>[A-Za-z_]\w*(?:%\w+)*)
    %(?P<section>{SECTION_NAMES})
    \(\s*(?P<idx>\d+)\s*\)
    %(?P<field>[A-Za-z_]\w*)
    \s*=\s*
    (?P<rhs>.*?)
    (?P<trail>\s*(?:!.*)?)$
    """,
    re.VERBOSE,
)

@dataclass
class CallBlock:
    start: int
    end: int
    indent: str
    text: str  # including embedded newlines, WITHOUT trailing newline guarantee


def _strip_comment_for_scan(line: str) -> str:
    """Remove Fortran trailing comment (!...) unless inside string. Used for scanning only."""
    out = []
    in_str = False
    quote = ""
    i = 0
    while i < len(line):
        ch = line[i]
        if in_str:
            out.append(ch)
            if ch == quote:
                # doubled quote inside string
                if i + 1 < len(line) and line[i + 1] == quote:
                    out.append(line[i + 1])
                    i += 2
                    continue
                in_str = False
            i += 1
            continue
        else:
            if ch in ("'", '"'):
                in_str = True
                quote = ch
                out.append(ch)
                i += 1
                continue
            if ch == "!":
                break
            out.append(ch)
            i += 1
    return "".join(out)


def iter_call_blocks(lines: List[str]) -> List[CallBlock]:
    """
    Collect multi-line 'call ...' blocks using typical Fortran continuation rules:
    - A line ending with '&' continues
    - Or next line beginning with '&' continues
    We keep original newlines exactly.
    """
    blocks: List[CallBlock] = []
    i = 0
    n = len(lines)

    def ends_with_amp(s: str) -> bool:
        s2 = _strip_comment_for_scan(s).rstrip()
        return s2.endswith("&")

    def starts_with_amp(s: str) -> bool:
        return s.lstrip().startswith("&")

    while i < n:
        m = CALL_START_RE.match(lines[i])
        if not m:
            i += 1
            continue

        indent = m.group(1)
        start = i
        buf = [lines[i]]
        i += 1
        while i < n:
            if ends_with_amp(buf[-1]) or starts_with_amp(lines[i]):
                buf.append(lines[i])
                i += 1
                continue
            break
        blocks.append(CallBlock(start=start, end=i, indent=indent, text="".join(buf)))
    return blocks


def find_matching_paren(s: str, open_pos: int) -> Optional[int]:
    """Find the matching ')' for '(' at open_pos in s, ignoring strings and comments."""
    depth = 0
    in_str = False
    quote = ""
    i = open_pos
    while i < len(s):
        ch = s[i]
        if in_str:
            if ch == quote:
                # doubled quote
                if i + 1 < len(s) and s[i + 1] == quote:
                    i += 2
                    continue
                in_str = False
            i += 1
            continue

        # comments: if we hit '!' outside string, skip to end of line
        if ch == "!":
            nl = s.find("\n", i)
            if nl == -1:
                return None
            i = nl + 1
            continue

        if ch in ("'", '"'):
            in_str = True
            quote = ch
            i += 1
            continue

        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return None


def split_args_with_spans(argstr: str) -> List[Tuple[str, int, int]]:
    """
    Split argument list on commas at depth 0, ignoring strings/comments.
    Returns list of (token_text, start_idx, end_idx) spans in *argstr*.
    Spans include the raw token region (no surrounding commas).
    """
    items: List[Tuple[str, int, int]] = []
    depth = 0
    in_str = False
    quote = ""
    i = 0
    tok_start = 0

    def emit(end_i: int):
        nonlocal tok_start
        tok = argstr[tok_start:end_i]
        items.append((tok, tok_start, end_i))
        tok_start = end_i + 1  # skip comma

    while i < len(argstr):
        ch = argstr[i]

        if in_str:
            if ch == quote:
                if i + 1 < len(argstr) and argstr[i + 1] == quote:
                    i += 2
                    continue
                in_str = False
            i += 1
            continue

        if ch == "!":
            # skip to end of line (comment)
            nl = argstr.find("\n", i)
            if nl == -1:
                # comment to end
                i = len(argstr)
                continue
            i = nl + 1
            continue

        if ch in ("'", '"'):
            in_str = True
            quote = ch
            i += 1
            continue

        if ch == "(":
            depth += 1
            i += 1
            continue
        if ch == ")":
            depth = max(0, depth - 1)
            i += 1
            continue

        if ch == "," and depth == 0:
            emit(i)
            i += 1
            continue

        i += 1

    # last token
    if tok_start <= len(argstr):
        items.append((argstr[tok_start:], tok_start, len(argstr)))

    # Trim pure-empty trailing tokens (rare)
    out = []
    for tok, a, b in items:
        if tok.strip() == "":
            continue
        out.append((tok, a, b))
    return out


def _is_int_literal(tok: str) -> bool:
    t = tok.replace("&", " ").strip()
    return re.fullmatch(r"[+-]?\d+", t) is not None


def refactor_new_call_in_block(block_text: str) -> Tuple[str, bool]:
    """
    Remove legacy count args from prg%new(...) calls inside a call block.
    Edits are done by deleting a contiguous slice of the argstring.
    """
    changed = False
    s = block_text
    pos = 0

    while True:
        m = NEW_CALL_RE.search(s, pos)
        if not m:
            break

        open_paren = s.find("(", m.end() - 1)
        if open_paren == -1:
            pos = m.end()
            continue

        close_paren = find_matching_paren(s, open_paren)
        if close_paren is None:
            pos = m.end()
            continue

        argstr = s[open_paren + 1:close_paren]
        args = split_args_with_spans(argstr)

        # legacy signature needs at least 12 args; args[4..10] are 7 integers
        if len(args) >= 12 and all(_is_int_literal(args[i][0]) for i in range(4, 11)):
            # spans in argstr coordinates
            _, start4, _ = args[4]
            _, _, end10 = args[10]

            # Expand removal start leftwards to include comma before arg #5 if present
            start_rm = start4
            j = start_rm - 1
            # include whitespace and '&' immediately before
            while j >= 0 and argstr[j] in " \t\r\n&":
                j -= 1
            if j >= 0 and argstr[j] == ",":
                start_rm = j  # include the comma separating arg4 and arg5

            # Expand removal end rightwards to include comma after arg #11 (the 7th count)
            end_rm = end10
            j = end_rm
            while j < len(argstr) and argstr[j] in " \t\r\n&":
                j += 1
            if j < len(argstr) and argstr[j] == ",":
                # include the comma after the last count, and any trailing whitespace/& after it
                j2 = j + 1
                while j2 < len(argstr) and argstr[j2] in " \t\r\n&":
                    j2 += 1
                end_rm = j2
            else:
                end_rm = end10

            new_argstr = argstr[:start_rm] + argstr[end_rm:]
            s = s[:open_paren + 1] + new_argstr + s[close_paren:]
            changed = True
            # continue scanning after this close paren
            pos = open_paren + 1
        else:
            pos = close_paren + 1

    return s, changed


def refactor_set_input_in_block(block_text: str) -> Tuple[str, bool]:
    """
    Convert %set_input('section', idx, ...) -> %add_input(UI_*, ...)
    preserving line breaks exactly.
    """
    def _sub(m: re.Match) -> str:
        section = m.group("section")
        ui = SECTION_UI[section]
        return f"%add_input({ui}, "

    new_s, n = SET_INPUT_HEAD_RE.subn(_sub, block_text)
    return new_s, (n > 0)


def refactor_posthoc_edits(lines: List[str]) -> Tuple[List[str], bool]:
    """
    Rewrite single-line post-hoc edits:
      obj%section(i)%field = rhs
    into a block using at()/replace_at().
    """
    changed = False
    out: List[str] = []

    for line in lines:
        m = POSTHOC_RE.match(line.rstrip("\n"))
        if not m:
            out.append(line)
            continue

        indent = m.group("indent")
        obj = m.group("obj")
        section = m.group("section")
        idx = m.group("idx")
        field = m.group("field")
        rhs = m.group("rhs")
        trail = m.group("trail") or ""

        # Avoid double-refactoring: if line already contains replace_at/at blocks, skip
        if "%replace_at" in line or "%at(" in line:
            out.append(line)
            continue

        # Build block; keep RHS verbatim (including any inline comment captured in trail)
        # Put trailing comment on the assignment line only.
        block_lines = [
            f"{indent}block\n",
            f"{indent}    type(ui_param) :: ptmp\n",
            f"{indent}    class(*), allocatable :: tmp\n",
            f"{indent}    call {obj}%{section}%at({idx}, tmp)\n",
            f"{indent}    select type(t => tmp)\n",
            f"{indent}    type is (ui_param)\n",
            f"{indent}        ptmp = t\n",
            f"{indent}    class default\n",
            f"{indent}        THROW_HARD('{obj}%{section}({idx}) not ui_param')\n",
            f"{indent}    end select\n",
            f"{indent}    if (allocated(tmp)) deallocate(tmp)\n",
            f"{indent}    ptmp%{field} = {rhs}{trail}\n",
            f"{indent}    call {obj}%{section}%replace_at({idx}, ptmp)\n",
            f"{indent}end block\n",
        ]
        out.extend(block_lines)
        changed = True

    return out, changed


def refactor_text(text: str) -> Tuple[str, bool]:
    lines = text.splitlines(keepends=True)

    # Pass 1: post-hoc edits (line-based)
    lines2, ch_post = refactor_posthoc_edits(lines)

    # Pass 2: call blocks (multiline-preserving)
    blocks = iter_call_blocks(lines2)
    if not blocks:
        return "".join(lines2), ch_post

    changed = ch_post
    new_lines = lines2[:]

    patches: List[Tuple[int, int, str]] = []
    for b in blocks:
        s = b.text

        s1, ch1 = refactor_new_call_in_block(s)
        s2, ch2 = refactor_set_input_in_block(s1)

        if ch1 or ch2:
            patches.append((b.start, b.end, s2))
            changed = True

    for start, end, repl in reversed(patches):
        new_lines[start:end] = [repl]

    return "".join(new_lines), changed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="Root directory or a single Fortran file")
    ap.add_argument("--write", action="store_true", help="Write changes back to files")
    ap.add_argument("--glob", default="*.f90,*.F90,*.f,*.for,*.F,*.fpp",
                    help="Comma-separated globs when path is a directory")
    ap.add_argument("--backup-ext", default=".bak", help="Backup extension when writing")
    args = ap.parse_args()

    root = Path(args.path)
    globs = [g.strip() for g in args.glob.split(",") if g.strip()]

    files: List[Path] = []
    if root.is_file():
        files = [root]
    else:
        for g in globs:
            files.extend(Path(p) for p in glob.glob(str(root / "**" / g), recursive=True))
        files = sorted(set(files))

    if not files:
        raise SystemExit("No matching files found.")

    any_changes = False
    for fp in files:
        try:
            orig = fp.read_text(encoding="utf-8", errors="replace")
        except Exception as e:
            print(f"[skip] {fp}: {e}")
            continue

        new_text, changed = refactor_text(orig)
        if not changed:
            continue

        any_changes = True
        if args.write:
            backup = fp.with_suffix(fp.suffix + args.backup_ext)
            backup.write_text(orig, encoding="utf-8")
            fp.write_text(new_text, encoding="utf-8")
            print(f"[write] {fp} (backup: {backup.name})")
        else:
            print(f"\n=== {fp} ===")
            diff = difflib.unified_diff(
                orig.splitlines(),
                new_text.splitlines(),
                fromfile=str(fp) + " (old)",
                tofile=str(fp) + " (new)",
                lineterm=""
            )
            for line in diff:
                print(line)

    if not any_changes:
        print("No changes detected.")


if __name__ == "__main__":
    main()
