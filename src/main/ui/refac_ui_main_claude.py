#!/usr/bin/env python3
"""
Format-preserving refactor for legacy ui_program constructors.

Transforms:
1) call prg%new(name, short, long, exec, n1,n2,n3,n4,n5,n6,n7, sp_required, [kw...])
   -> call prg%new(name, short, long, exec, sp_required, [kw...])
   Removes the 7 count integers, preserving line breaks.

2) call prg%set_input('parm_ios', idx, param, ...)
   -> call prg%add_input_param(UI_PARM, param, ...)
   Preserves multiline formatting.

3) Post-hoc edits like prg%filt_ctrls(6)%descr_long = '...'
   -> Commented out with TODO for manual handling

Usage:
  python refactor.py input.f90 --write
  python refactor.py src/ --write  # process directory
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

# Match set_input with quoted section name and integer index
SET_INPUT_RE = re.compile(
    rf"""(?ix)
    (?P<prefix>%)set_input\s*\(\s*
    (?P<q>['"])(?P<section>{SECTION_NAMES})(?P=q)\s*,\s*
    (?P<idx>\d+)\s*,\s*
    """,
    re.VERBOSE,
)

# Match %new( ... ) 
NEW_CALL_RE = re.compile(r"(?i)%new\s*\(")

# Post-hoc edit: obj%section(idx)%field = rhs
POSTHOC_RE = re.compile(
    rf"""^
    (?P<indent>\s*)
    (?P<obj>\w+)
    %(?P<section>{SECTION_NAMES})
    \(\s*(?P<idx>\d+)\s*\)
    %(?P<field>\w+)
    \s*=\s*
    (?P<rhs>.+?)
    (?P<comment>\s*!.*)?$
    """,
    re.VERBOSE,
)


@dataclass
class CallBlock:
    start: int
    end: int  
    text: str


def strip_fortran_comment(line: str) -> str:
    """Remove trailing ! comment, respecting strings."""
    in_str = False
    quote = None
    for i, ch in enumerate(line):
        if in_str:
            if ch == quote:
                # Check for doubled quote
                if i + 1 < len(line) and line[i + 1] == quote:
                    continue
                in_str = False
        elif ch in ('"', "'"):
            in_str = True
            quote = ch
        elif ch == '!':
            return line[:i]
    return line


def iter_call_blocks(lines: List[str]) -> List[CallBlock]:
    """Collect multi-line call statements using & continuation."""
    blocks = []
    i = 0
    
    def has_continuation(s: str) -> bool:
        clean = strip_fortran_comment(s).rstrip()
        return clean.endswith('&')
    
    def is_continuation(s: str) -> bool:
        return s.lstrip().startswith('&')
    
    while i < len(lines):
        if not CALL_START_RE.match(lines[i]):
            i += 1
            continue
            
        start = i
        buf = [lines[i]]
        i += 1
        
        while i < len(lines):
            if has_continuation(buf[-1]) or is_continuation(lines[i]):
                buf.append(lines[i])
                i += 1
            else:
                break
                
        blocks.append(CallBlock(start=start, end=i, text=''.join(buf)))
    
    return blocks


def find_matching_paren(text: str, start: int) -> Optional[int]:
    """Find matching ) for ( at start position, ignoring strings/comments."""
    depth = 0
    in_str = False
    quote = None
    
    i = start
    while i < len(text):
        ch = text[i]
        
        if in_str:
            if ch == quote:
                if i + 1 < len(text) and text[i + 1] == quote:
                    i += 2
                    continue
                in_str = False
            i += 1
            continue
        
        # Handle comments
        if ch == '!':
            nl = text.find('\n', i)
            if nl == -1:
                return None
            i = nl + 1
            continue
            
        if ch in ('"', "'"):
            in_str = True
            quote = ch
            i += 1
            continue
            
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
            if depth == 0:
                return i
        i += 1
    
    return None


def tokenize_args(argstr: str) -> List[Tuple[str, int, int]]:
    """
    Split arguments on commas at depth 0, ignoring strings/comments.
    Returns [(token, start_pos, end_pos), ...]
    """
    tokens = []
    depth = 0
    in_str = False
    quote = None
    tok_start = 0
    i = 0
    
    def emit_token(end: int):
        nonlocal tok_start
        if tok_start < end:
            tokens.append((argstr[tok_start:end], tok_start, end))
        tok_start = end + 1
    
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
            i += 1
            continue
        
        if ch == '!':
            nl = argstr.find('\n', i)
            if nl == -1:
                i = len(argstr)
            else:
                i = nl + 1
            continue
            
        if ch in ('"', "'"):
            in_str = True
            quote = ch
            i += 1
            continue
            
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth = max(0, depth - 1)
        elif ch == ',' and depth == 0:
            emit_token(i)
            i += 1
            continue
            
        i += 1
    
    emit_token(len(argstr))
    
    # Filter empty tokens
    result = []
    for tok, start, end in tokens:
        clean = tok.replace('&', '').strip()
        if clean:
            result.append((tok, start, end))
    
    return result


def is_integer_literal(token: str) -> bool:
    """Check if token is an integer literal."""
    cleaned = token.replace('&', '').replace('\n', '').replace('\r', '').strip()
    int_pattern = r'[+-]?\d+'
    return bool(re.fullmatch(int_pattern, cleaned))


def refactor_new_call(text: str, debug: bool = False) -> Tuple[str, bool]:
    """Remove 7 count arguments from %new() calls."""
    changed = False
    pos = 0
    
    while True:
        m = NEW_CALL_RE.search(text, pos)
        if not m:
            break
            
        open_pos = m.end() - 1
        close_pos = find_matching_paren(text, open_pos)
        if close_pos is None:
            pos = m.end()
            continue
            
        argstr = text[open_pos + 1:close_pos]
        args = tokenize_args(argstr)
        
        if debug:
            print(f"\nFound %new() with {len(args)} args")
            for i, (tok, _, _) in enumerate(args[:15]):
                raw = tok.replace('\n', '\\n')[:60]
                clean = tok.replace('&', '').replace('\n', '').replace('\r', '').strip()[:60]
                print(f"  arg[{i}]: '{clean}' (raw: {raw})")
        
        # Need at least 12 args
        if len(args) >= 12:
            counts = [args[i] for i in range(4, 11)]
            is_all_int = all(is_integer_literal(tok) for tok, _, _ in counts)
            
            if debug:
                print(f"  Checking args[4:11] for integers: {is_all_int}")
                for i, (tok, _, _) in enumerate(counts, 4):
                    print(f"    arg[{i}]: '{tok.strip()}' -> {is_integer_literal(tok)}")
            
            if is_all_int:
                _, start_rm, _ = args[4]
                _, _, end_rm = args[10]
                
                if debug:
                    print(f"  Initial removal span: [{start_rm}:{end_rm}]")
                    print(f"  Text to remove initially: '{argstr[start_rm:end_rm]}'")
                
                # Extend left to include comma before arg[4]
                temp = start_rm - 1
                while temp >= 0 and argstr[temp] in ' \t\n\r&':
                    temp -= 1
                if temp >= 0 and argstr[temp] == ',':
                    start_rm = temp
                    if debug:
                        print(f"  Extended left to position {start_rm} (included comma)")
                
                # Extend right to include comma after arg[10]
                temp = end_rm
                while temp < len(argstr) and argstr[temp] in ' \t\n\r&':
                    temp += 1
                if temp < len(argstr) and argstr[temp] == ',':
                    temp += 1
                    while temp < len(argstr) and argstr[temp] in ' \t\n\r&':
                        temp += 1
                    end_rm = temp
                    if debug:
                        print(f"  Extended right to position {end_rm} (included comma+whitespace)")
                
                # Check for comment about "entries in each group"
                line_start = argstr.rfind('\n', 0, start_rm) + 1
                line_end = argstr.find('\n', end_rm)
                if line_end == -1:
                    line_end = len(argstr)
                
                line_with_counts = argstr[line_start:line_end]
                comment_match = re.search(r'!\s*#\s*entries\s+in\s+each\s+group[^!\n]*', line_with_counts)
                if comment_match:
                    comment_start = line_start + comment_match.start()
                    comment_end = line_start + comment_match.end()
                    if comment_start >= end_rm:
                        end_rm = comment_end
                        if debug:
                            print(f"  Extended to remove comment at [{comment_start}:{comment_end}]")
                
                if debug:
                    print(f"  Final removal span: [{start_rm}:{end_rm}]")
                    print(f"  Will remove: '{argstr[start_rm:end_rm]}'")
                    print(f"  Result will be: '{argstr[:start_rm]}' + '{argstr[end_rm:]}'")
                
                new_argstr = argstr[:start_rm] + argstr[end_rm:]
                new_full_text = text[:open_pos + 1] + new_argstr + text[close_pos:]
                
                if debug:
                    print(f"  Original argstr length: {len(argstr)}")
                    print(f"  New argstr length: {len(new_argstr)}")
                    print(f"  ✓ Removed 7 integer args")
                    # Show a snippet of the new result
                    snippet_start = max(0, open_pos - 50)
                    snippet_end = min(len(new_full_text), close_pos + 50)
                    print(f"  New text snippet: {new_full_text[snippet_start:snippet_end]}")
                
                text = new_full_text
                changed = True
                pos = open_pos + 1
            else:
                pos = close_pos + 1
        else:
            if debug:
                print(f"  Not enough args ({len(args)} < 12), skipping")
            pos = close_pos + 1
    
    return text, changed


def refactor_set_input(text: str) -> Tuple[str, bool]:
    """Convert set_input('section', idx, ...) to add_input_param(UI_*, ...)."""
    def replace(m: re.Match) -> str:
        section = m.group('section')
        ui_const = SECTION_UI[section]
        return f"%add_input_param({ui_const}, "
    
    new_text, count = SET_INPUT_RE.subn(replace, text)
    return new_text, count > 0


def refactor_posthoc_edits(lines: List[str]) -> Tuple[List[str], bool]:
    """Comment out post-hoc edits with TODO markers."""
    changed = False
    out = []
    
    for line in lines:
        m = POSTHOC_RE.match(line.rstrip('\n'))
        if not m:
            out.append(line)
            continue
        
        if line.strip().startswith('!'):
            out.append(line)
            continue
            
        indent = m.group('indent')
        obj = m.group('obj')
        section = m.group('section')
        idx = m.group('idx')
        field = m.group('field')
        rhs = m.group('rhs')
        
        out.append(f"{indent}! TODO: Manual fix needed - was: {obj}%{section}({idx})%{field} = {rhs}\n")
        out.append(f"{indent}! The new API doesn't support array indexing. Add {field}={rhs} to add_input_param() call.\n")
        out.append(f"{indent}! {line.rstrip()}\n")
        changed = True
    
    return out, changed


def refactor_text(text: str, debug: bool = False) -> Tuple[str, bool]:
    """Main refactoring function."""
    lines = text.splitlines(keepends=True)
    
    # Pass 1: Handle post-hoc edits
    lines, changed_post = refactor_posthoc_edits(lines)
    
    # Pass 2: Handle call blocks
    blocks = iter_call_blocks(lines)
    changed_calls = False
    
    if blocks:
        patches = []
        for block in blocks:
            old_text = block.text
            new_text = block.text
            new_text, ch1 = refactor_new_call(new_text, debug=debug)
            new_text, ch2 = refactor_set_input(new_text)
            
            if ch1 or ch2:
                if debug:
                    print(f"\nBlock changed at lines {block.start}-{block.end}")
                    if ch1:
                        print(f"  - Refactored %new() call")
                    if ch2:
                        print(f"  - Refactored set_input() calls")
                    if old_text != new_text:
                        print(f"  - Text actually changed: {len(old_text)} -> {len(new_text)} chars")
                    else:
                        print(f"  - WARNING: Flags say changed but text is identical!")
                patches.append((block.start, block.end, new_text))
                changed_calls = True
        
        if debug:
            print(f"\nApplying {len(patches)} patches to lines array")
        
        for start, end, repl in reversed(patches):
            if debug:
                print(f"  Replacing lines[{start}:{end}] with new content")
            lines[start:end] = [repl]
    
    return ''.join(lines), changed_post or changed_calls


def main():
    ap = argparse.ArgumentParser(description="Refactor Fortran ui_program constructors")
    ap.add_argument('path', help='File or directory to process')
    ap.add_argument('--write', action='store_true', help='Write changes to files')
    ap.add_argument('--backup', default='.bak', help='Backup extension')
    ap.add_argument('--glob', default='*.f90,*.F90', help='File patterns for directories')
    ap.add_argument('--debug', action='store_true', help='Print debug info')
    args = ap.parse_args()
    
    root = Path(args.path)
    patterns = [p.strip() for p in args.glob.split(',')]
    
    files = []
    if root.is_file():
        files = [root]
    else:
        for pattern in patterns:
            files.extend(root.rglob(pattern))
        files = sorted(set(files))
    
    if not files:
        print("No files found")
        return
    
    any_changed = False
    for fp in files:
        try:
            orig = fp.read_text(encoding='utf-8')
        except Exception as e:
            print(f"[skip] {fp}: {e}")
            continue
        
        new_text, changed = refactor_text(orig, debug=args.debug)
        if not changed:
            continue
        
        any_changed = True
        
        if args.write:
            backup = fp.with_suffix(fp.suffix + args.backup)
            backup.write_text(orig, encoding='utf-8')
            fp.write_text(new_text, encoding='utf-8')
            print(f"[wrote] {fp} (backup: {backup.name})")
        else:
            print(f"\n{'='*60}\n{fp}\n{'='*60}")
            diff = difflib.unified_diff(
                orig.splitlines(),
                new_text.splitlines(),
                fromfile=f'{fp} (old)',
                tofile=f'{fp} (new)',
                lineterm=''
            )
            for line in diff:
                print(line)
    
    if not any_changed:
        print("No changes needed")


if __name__ == '__main__':
    main()