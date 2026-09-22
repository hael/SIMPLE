#!/bin/bash
# scripts/run_fast_gate.sh — run the build-time test gate on an installed
# --compile-tests build and check it against the budget.
#
#   scripts/run_fast_gate.sh [BUILD_DIR]        (default: build)
#
# Runs `ctest -L "fast|provisional"` in BUILD_DIR with half the cores, tees
# the output to BUILD_DIR/test_runs/ctest_fast.log and hands it to
# scripts/ctest_budget.py. Exit status is the budget checker's: non-zero on
# a failed entry or, the gate being declared (Phase 2 of
# doc/refactoring_notes/uniform_test_environment_refactoring.md, 2026-09-22),
# on a run over 30 s of real time. The checker keeps the per-entry table
# beside the log, so a suite that grows is visible from build to build.
#
# Called by every compile_*.sh after `make install` when --compile-tests is
# given, and by `make check`.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="${1:-$ROOT/build}"
case "$BUILD" in /*) ;; *) BUILD="$ROOT/$BUILD" ;; esac
[ -f "$BUILD/CTestTestfile.cmake" ] || { echo "run_fast_gate: no CTest configuration in $BUILD (build with --compile-tests)" >&2; exit 2; }
if command -v nproc >/dev/null 2>&1; then ncpu=$(nproc); else ncpu=$(sysctl -n hw.ncpu); fi
jobs=$(( ncpu / 2 )); [ "$jobs" -lt 1 ] && jobs=1
mkdir -p "$BUILD/test_runs"
LOG="$BUILD/test_runs/ctest_fast.log"
# GATE_DECLARED went to yes in Phase 2 (2026-09-22): the fast area suites
# exist and SIMPLE_CTEST_BUDGET is armed in production/CMakeLists.txt.
GATE_DECLARED=yes
echo "== fast gate: ctest -L 'fast|provisional' -j$jobs in $BUILD"
( cd "$BUILD" && ctest -L "fast|provisional" --output-on-failure --parallel "$jobs" --timeout 600 2>&1 ) | tee "$LOG"
if [ "$GATE_DECLARED" = yes ]; then
    python3 "$ROOT/scripts/ctest_budget.py" "$LOG" --budget 30
else
    python3 "$ROOT/scripts/ctest_budget.py" "$LOG" --no-budget
fi
