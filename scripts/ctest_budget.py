#!/usr/bin/env python3
"""ctest_budget.py -- enforce the fast-gate time budget.

Reads the stdout of a `ctest -L fast` run (tee it to a file; see
scripts/run_fast_gate.sh) and checks the one rule of section 5.1 of
doc/refactoring_notes/uniform_test_environment_refactoring.md: the whole run
finishes within --budget seconds of real time ("Total Test time (real)" as
ctest reports it). Any failed entry also fails the check.

Writes every entry sorted by time next to the log as <log>.timing.txt, so the
numbers are kept from build to build and a suite that grows is visible, and
prints the same table unless --quiet is given. With --quiet (the compile
scripts) a passing run within budget prints nothing, so what the developer
sees is ctest's own report, as in X; the table and the problems are printed
only when there is something to fix. Exits 1 on a failed entry or a run over
budget, 2 if the log holds no ctest results.

usage: scripts/ctest_budget.py LOG [--budget 30] [--no-budget] [--quiet]
"""
import argparse
import re
import sys

RE_TEST = re.compile(r'^\s*\d+/\d+\s+Test\s+#\d+:\s+(\S+)\s+\.+\s*(?:\*+)?(Passed|Failed|Timeout|Not Run|Exception|Subprocess aborted|\*\*\*[^\d]*)\s+([\d.]+)\s+sec', re.M)
RE_TOTAL = re.compile(r'^Total Test time \(real\)\s*=\s*([\d.]+)\s+sec', re.M)
RE_LABEL = re.compile(r'^(\w+)\s*=\s*([\d.]+)\s+sec\*proc', re.M)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('log')
    ap.add_argument('--budget', type=float, default=30.0, help='seconds of real time for the whole run (default 30)')
    ap.add_argument('--no-budget', action='store_true', help='report times but do not fail on the budget (only failed entries count)')
    ap.add_argument('--quiet', action='store_true', help='print nothing when the run passes within budget (the table still goes to <log>.timing.txt)')
    a = ap.parse_args()
    text = open(a.log, errors='replace').read()
    tests = [(name, status.strip(), float(secs)) for name, status, secs in RE_TEST.findall(text)]
    if not tests:
        print('ctest_budget: no ctest results in %s' % a.log)
        return 2
    total = RE_TOTAL.search(text)
    total = float(total.group(1)) if total else sum(t[2] for t in tests)
    labels = dict((k, float(v)) for k, v in RE_LABEL.findall(text))
    tests.sort(key=lambda t: -t[2])
    lines = ['%8s  %-8s  %s' % ('sec', 'status', 'test')]
    problems = []
    for name, status, secs in tests:
        flag = ''
        if not status.startswith('Passed'):
            flag = '  <-- ' + status
            problems.append('%s: %s' % (name, status))
        lines.append('%8.2f  %-8s  %s%s' % (secs, status.split()[0], name, flag))
    lines.append('')
    over = total > a.budget and not a.no_budget
    lines.append('%d entries; total real time %.1f s (budget %.0f s)%s' % (
        len(tests), total, a.budget, '  <-- OVER BUDGET' if over else ''))
    if labels:
        lines.append('label proc-seconds: ' + ', '.join('%s %.1f' % kv for kv in sorted(labels.items())))
    if over:
        problems.append('total real time %.1f s exceeds the %.0f s budget' % (total, a.budget))
    out = '\n'.join(lines)
    with open(a.log + '.timing.txt', 'w') as fh:
        fh.write(out + '\n')
    if problems:
        print(out)
        print('\nctest_budget: %d problem(s):' % len(problems))
        for p in problems:
            print('   ' + p)
        return 1
    if not a.quiet:
        print(out)
        print('\nctest_budget: within budget')
    return 0


if __name__ == '__main__':
    sys.exit(main())
