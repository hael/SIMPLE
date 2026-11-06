#!/usr/bin/env python3
"""check_tester_linkage.py -- a submodule tester calls no private plain procedure of its module.

gfortran before 16 gives a private module procedure internal linkage, so a submodule tester (another
object file) that calls one by name compiles and then fails to link: on the Linux builds, not on a
macOS build with gfortran 16. Private type-bound procedures (also called by name: they are bound),
private constants and variables, and separate module procedures link. A procedure a submodule tester
calls is public, bound to a type, or a separate module procedure with its body in a production submodule
(doc/policies/test_environment_policy.md, section 4.3).
usage: scripts/check_tester_linkage.py [ROOT]
exit status: 0 none found, 1 a tester calls a private plain procedure
"""

import os
import re
import subprocess
import sys

KEYWORDS = {"write", "read", "open", "close", "print", "inquire", "flush", "rewind", "backspace", "format"}
# prefixes in any order: pure, elemental, recursive, impure, a type (real(dp), type(x)), module
PROC_RE = re.compile(r"\s*(?:(?:pure|elemental|recursive|impure|non_recursive|integer|real|logical|complex|character"
                     r"|double\s+precision|type\s*\(\w+\)|class\s*\(\w+\))(?:\s*\([^)]*\))?\s+)*"
                     r"(module\s+)?(?:(?:pure|elemental|recursive|impure)\s+)*(subroutine|function)\s+(\w+)")


def code_lines(path):
    """lower-case source lines without strings and comments"""
    out = []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = re.sub(r"'[^']*'|\"[^\"]*\"", "''", line)
            out.append(re.sub(r"!.*", "", line).rstrip().lower())
    return out


def module_contains(lines):
    """index of the module-level contains (not a type's), or None"""
    depth = 0
    for i, l in enumerate(lines):
        if re.match(r"\s*type\s*(,[^:]*)?::|\s*type\s+\w+\s*$", l):
            depth += 1
        elif re.match(r"\s*end\s*type\b", l):
            depth -= 1
        elif depth == 0 and re.match(r"\s*contains\s*$", l):
            return i
    return None


def private_plain_procedures(lines):
    ic = module_contains(lines)
    if ic is None:
        return set()
    spec = lines[:ic]
    if not any(l.strip() == "private" for l in spec):
        return set()
    public, bound = set(), set()
    i = 0
    while i < ic:
        l = spec[i]
        m = re.match(r"\s*public\s*::(.*)", l)
        if m:
            text = m.group(1)
            while text.rstrip().endswith("&") and i + 1 < ic:
                i += 1
                text = text.rstrip()[:-1] + " " + spec[i]
            public |= {x.strip().strip("&").strip() for x in text.split(",") if x.strip().strip("&").strip()}
        m = re.match(r"\s*procedure\b[^:]*::(.*)", l)
        if m:
            bound |= {p.split("=>")[-1].strip() for p in m.group(1).split(",")}
        i += 1
    plain, depth = set(), 0
    for l in lines[ic + 1:]:
        if re.match(r"\s*end\s+(subroutine|function)\b", l):
            depth -= 1
            continue
        m = PROC_RE.match(l)
        if m and not re.match(r"\s*module\s+procedure\b", l):
            if depth == 0 and not m.group(1):
                plain.add(m.group(3))
            depth += 1
    return {n for n in plain if n not in public and n not in bound and n not in KEYWORDS}


def main():
    root = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), ".."))
    files = subprocess.run(["rg", "--files", "-g", "*.f90", "src"], cwd=root, capture_output=True,
                           text=True).stdout.split()
    modfile = {}
    for f in files:
        for l in code_lines(os.path.join(root, f)):
            m = re.match(r"\s*module\s+(\w+)\s*$", l)
            if m and m.group(1) not in ("procedure", "subroutine", "function"):
                modfile[m.group(1)] = f
    bad = []
    for f in files:
        if not f.endswith("_tester.f90"):
            continue
        lines = code_lines(os.path.join(root, f))
        head = next((l for l in lines if re.match(r"\s*submodule\s*\(", l)), None)
        if head is None:
            continue
        host = re.match(r"\s*submodule\s*\(\s*(\w+)", head).group(1)
        if host not in modfile:
            continue
        for name in sorted(private_plain_procedures(code_lines(os.path.join(root, modfile[host])))):
            call = re.compile(r"(?<![%\w])" + re.escape(name) + r"\s*\(|\bcall\s+" + re.escape(name) + r"\b")
            for n, l in enumerate(lines, 1):
                if call.search(l) and not re.search(r"\b(subroutine|function)\s+" + re.escape(name) + r"\b", l):
                    bad.append("%s:%d calls %s, a private module procedure of %s" % (f, n, name, host))
                    break
    for b in bad:
        print("check_tester_linkage:", b)
    if bad:
        print("check_tester_linkage: make each public, bind it to a type, or give it a separate module procedure "
              "body in a production submodule (gfortran before 16 does not link it)")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
