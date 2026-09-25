#!/bin/bash
# scripts/run_fast_gate.sh — run the build-time test gate on a build with tests
# (BUILD_TESTS=ON) and check it against the budget.
#
#   scripts/run_fast_gate.sh [BUILD_DIR]        (default: build)
#
# Runs `ctest -L "fast|provisional"` in BUILD_DIR with half the cores, tees
# the output to BUILD_DIR/test_runs/ctest_fast.log and hands it to
# scripts/ctest_budget.py. What is printed is ctest's own report, as in X;
# the budget checker is silent when the gate passes within budget and speaks
# only when an entry failed or the run went over 30 s (Phase 2 of
# doc/refactoring_notes/completed/uniform_test_environment_refactoring.md, 2026-09-22).
# The per-entry timing table is always written beside the log
# (ctest_fast.log.timing.txt), so a suite that grows is visible from build to
# build. Exit status is the budget checker's.
#
# Before ctest, scripts/check_test_registry.py checks that the CTest
# registrations, the test UI and the test routers agree (plan, section 7);
# a mismatch fails the gate with status 1 before any test runs.
#
# Called by every compile_*.sh between `make` and `make install` unless
# --exclude-tests is given (a failed gate installs nothing), and by
# `make check`. The fast suites run in-process from the build tree and need
# nothing from the install tree.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="${1:-$ROOT/build}"
case "$BUILD" in /*) ;; *) BUILD="$ROOT/$BUILD" ;; esac
[ -f "$BUILD/CTestTestfile.cmake" ] || { echo "run_fast_gate: no CTest configuration in $BUILD (built with --exclude-tests? rebuild without it)" >&2; exit 2; }
if command -v nproc >/dev/null 2>&1; then ncpu=$(nproc); else ncpu=$(sysctl -n hw.ncpu); fi
jobs=$(( ncpu / 2 )); [ "$jobs" -lt 1 ] && jobs=1
python3 "$ROOT/scripts/check_test_registry.py" "$ROOT" || exit 1
mkdir -p "$BUILD/test_runs"
LOG="$BUILD/test_runs/ctest_fast.log"
# GATE_DECLARED went to yes in Phase 2 (2026-09-22): the fast area suites
# exist and SIMPLE_CTEST_BUDGET is armed in production/CMakeLists.txt.
GATE_DECLARED=yes
( cd "$BUILD" && ctest -L "fast|provisional" --output-on-failure --parallel "$jobs" --timeout 600 2>&1 ) | tee "$LOG"
if [ "$GATE_DECLARED" = yes ]; then
    python3 "$ROOT/scripts/ctest_budget.py" "$LOG" --budget 30 --quiet
else
    python3 "$ROOT/scripts/ctest_budget.py" "$LOG" --no-budget --quiet
fi
