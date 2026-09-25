#!/usr/bin/env python3
"""check_test_registry.py -- the CTest registrations, the test UI and the test routers agree.

Registry consistency (doc/refactoring_notes/completed/uniform_test_environment_refactoring.md,
section 7, item 6). A static check of the source tree, run by scripts/run_fast_gate.sh
before the gate, so a build with tests fails on a mismatch:

  * every selector production/CMakeLists.txt registers (test=<id>) is a program of the
    test UI (src/main/ui/simple_test) with exactly one router case (src/main/exec/
    simple_test_exec_*.f90);
  * every program of the test UI has exactly one router case, and every router case is
    a program of the test UI;
  * every area suite (unit_<area>) and library suite (lib_<area>) of the test UI is
    registered; `units`, the umbrella of the unit suites, is not, by design;
  * every program that runs sub-suites (the tables of simple_commanders_test_class) lists
    them, in table order, in the help of its suite= input in the test UI, spelled as the
    runner matches them (suite_id: lower case, blanks and hyphens as underscores, commas
    and slashes dropped); `units`, which runs every table, has no list.

Programs that are neither suites nor registered (manual tools, the focused routes of
lib_single) are allowed. Nothing is compiled or run.

usage: scripts/check_test_registry.py [ROOT] [--verbose]
exit status: 0 consistent, 1 inconsistent, 2 a source file is missing
"""
import glob
import os
import re
import sys

SUITE_RE = re.compile(r'^(unit|lib)_\w+$')
UMBRELLA = {'units'}


def read(path):
    with open(path, encoding='utf-8', errors='replace') as f:
        return f.read()


def strip_fortran_comments(text):
    return '\n'.join(line.split('!', 1)[0] for line in text.split('\n'))


def registered_selectors(cmake_text):
    """test=<id> selectors of the CTest registrations, with foreach(<var> IN ITEMS ...) expanded"""
    lines = [l.split('#', 1)[0] for l in cmake_text.split('\n')]
    text = '\n'.join(lines)
    loops = {m.group(1): m.group(2).split()
             for m in re.finditer(r'foreach\(\s*(\w+)\s+IN\s+ITEMS\s+([^)]*)\)', text)}
    selectors = set()
    for m in re.finditer(r'test=([A-Za-z0-9_]*(?:\$\{\w+\}[A-Za-z0-9_]*)*)', text):
        sel = m.group(1)
        var = re.search(r'\$\{(\w+)\}', sel)
        if var:
            if var.group(1) not in loops:
                raise SystemExit('check_test_registry: cannot expand %s in production/CMakeLists.txt' % sel)
            for item in loops[var.group(1)]:
                selectors.add(sel.replace('${%s}' % var.group(1), item))
        elif sel:
            selectors.add(sel)
    return selectors


def ui_programs(ui_dir):
    progs = {}
    for f in sorted(glob.glob(os.path.join(ui_dir, '*.f90'))):
        for m in re.finditer(r"call\s+add_ui_program\(\s*'([A-Za-z0-9_]+)'", strip_fortran_comments(read(f))):
            progs.setdefault(m.group(1), []).append(os.path.basename(f))
    return progs


def router_cases(exec_dir):
    cases = {}
    for f in sorted(glob.glob(os.path.join(exec_dir, 'simple_test_exec_*.f90'))):
        for m in re.finditer(r"\bcase\s*\(([^)]*)\)", strip_fortran_comments(read(f)), re.I):
            for name in re.findall(r"'([A-Za-z0-9_]+)'", m.group(1)):
                cases.setdefault(name, []).append(os.path.basename(f))
    return cases


def suite_id(name):
    """the selector of a sub-suite, as suite_id in simple_commanders_test_class"""
    out = []
    for ch in name.lower().rstrip():
        if ch in ',/':
            continue
        out.append('_' if ch in ' -' else ch)
    return ''.join(out)


def suite_tables(class_text):
    """sub-suite selectors per test program: the suites_<x> tables (unit_<area> runs
    suites_<area>, lib_<area> runs suites_lib_<area>) and the add_suite calls made directly
    in an exec_test_<program> routine"""
    text = strip_fortran_comments(class_text)
    add = r"add_suite\(\s*s\s*,\s*n\s*,\s*'([^']+)'"
    programs = {}
    for m in re.finditer(r'subroutine\s+suites_(\w+)\s*\((.*?)end\s+subroutine\s+suites_\1\b', text, re.S | re.I):
        key = m.group(1)
        prog = key if key.startswith('lib_') else 'unit_' + key
        programs[prog] = [suite_id(x) for x in re.findall(add, m.group(2))]
    for m in re.finditer(r'subroutine\s+exec_test_(\w+)\s*\((.*?)end\s+subroutine\s+exec_test_\1\b', text, re.S | re.I):
        names = re.findall(add, m.group(2))
        if names:
            programs[m.group(1)] = [suite_id(x) for x in names]
    return programs


def ui_suite_lists(ui_dir):
    """the sub-suites each test UI program lists in the help of its suite= input"""
    lists = {}
    for f in sorted(glob.glob(os.path.join(ui_dir, '*.f90'))):
        text = strip_fortran_comments(read(f)).replace('&\n', '').replace('&', '')
        for m in re.finditer(r"(\w+)%add_input\(\s*UI_PARM\s*,\s*'suite'.*?to run alone \(([^)]*)\)", text, re.S):
            lists[m.group(1)] = [x.strip() for x in m.group(2).split(',') if x.strip()]
    return lists


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    verbose = '--verbose' in sys.argv[1:]
    root = os.path.abspath(args[0] if args else os.path.join(os.path.dirname(__file__), '..'))
    cmake = os.path.join(root, 'production', 'CMakeLists.txt')
    test_class = os.path.join(root, 'src', 'main', 'commanders', 'test', 'simple_commanders_test_class.f90')
    ui_dir = os.path.join(root, 'src', 'main', 'ui', 'simple_test')
    exec_dir = os.path.join(root, 'src', 'main', 'exec')
    for p in (cmake, ui_dir, exec_dir, test_class):
        if not os.path.exists(p):
            print('check_test_registry: missing %s' % p)
            return 2
    registered = registered_selectors(read(cmake))
    progs = ui_programs(ui_dir)
    cases = router_cases(exec_dir)
    problems = []
    for name, files in sorted(progs.items()):
        if len(files) > 1:
            problems.append("test program '%s' is defined twice in the test UI (%s)" % (name, ', '.join(files)))
    for sel in sorted(registered):
        if sel not in progs:
            problems.append("CTest registers test=%s, which is not a program of the test UI" % sel)
    for name in sorted(progs):
        n = len(cases.get(name, []))
        if n == 0:
            problems.append("test program '%s' has no router case in src/main/exec/simple_test_exec_*.f90" % name)
        elif n > 1:
            problems.append("test program '%s' has %d router cases (%s)" % (name, n, ', '.join(cases[name])))
    for name in sorted(cases):
        if name not in progs:
            problems.append("router case '%s' (%s) is not a program of the test UI" % (name, ', '.join(cases[name])))
    for name in sorted(progs):
        if SUITE_RE.match(name) and name not in UMBRELLA and name not in registered:
            problems.append("suite '%s' is in the test UI but not registered with CTest" % name)
    tables = suite_tables(read(test_class))
    ui_lists = ui_suite_lists(ui_dir)
    for name in sorted(set(tables) | set(ui_lists)):
        if name in UMBRELLA or name not in progs:
            continue
        if name not in ui_lists:
            problems.append("test program '%s' runs sub-suites but its test UI has no suite= list" % name)
        elif name not in tables:
            problems.append("test program '%s' lists sub-suites in its test UI but runs no suite table" % name)
        elif ui_lists[name] != tables[name]:
            problems.append("test program '%s': the suite= list of the test UI (%s) is not the suite table (%s)"
                            % (name, ', '.join(ui_lists[name]), ', '.join(tables[name])))
    if problems:
        print('check_test_registry: %d problem(s):' % len(problems))
        for p in problems:
            print('   ' + p)
        return 1
    if verbose:
        print('check_test_registry: %d registered selectors, %d test programs, %d router cases, '
              '%d suite tables: consistent' % (len(registered), len(progs), len(cases), len(tables)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
