#!/bin/bash
# Tests are built by default (simple_test_exec and the *_tester modules) and the fast
# gate runs before installation; --exclude-tests builds the library and executables only.
BUILD_TESTS=ON
ROOT="$(cd "$(dirname "$0")" && pwd)"
source "$ROOT/scripts/simple_build_jobs.sh" || exit $?
for arg in "$@"; do
    case "$arg" in
        --exclude-tests) BUILD_TESTS=OFF ;;
        -h|--help)
            echo "usage: $(basename "$0") [--exclude-tests]"
            echo "Default jobs: available CPUs, or CMAKE_BUILD_PARALLEL_LEVEL when set."
            exit 0
            ;;
        *) echo "$(basename "$0"): unknown option: $arg (see --help)" >&2; exit 1 ;;
    esac
done
JOBS=$(simple_build_jobs) || exit $?

# NICE pins NumPy 2.2.6, whose binary wheels support CPython 3.10-3.13.
# Prefer 3.12 for a stable GUI environment instead of letting CMake select a
# newer, incompatible python3 from PATH (for example Miniforge Python 3.14).
if [ -n "${NICE_PYTHON3:-}" ]; then
    NICE_PYTHON_CANDIDATES=("$NICE_PYTHON3")
else
    NICE_PYTHON_CANDIDATES=(python3.12 python3.13 python3.11 python3.10 python3)
fi
NICE_PYTHON=""
for candidate in "${NICE_PYTHON_CANDIDATES[@]}"; do
    candidate_path=$(command -v "$candidate" 2>/dev/null) || continue
    if "$candidate_path" -c \
        'import sqlite3, sys, venv; raise SystemExit(0 if (3, 10) <= sys.version_info[:2] < (3, 14) else 1)' \
        >/dev/null 2>&1; then
        NICE_PYTHON="$candidate_path"
        break
    fi
done
if [ -z "$NICE_PYTHON" ]; then
    echo "$(basename "$0"): NICE requires Python 3.10-3.13 with sqlite3 and venv." >&2
    echo "Install Python 3.12 or set NICE_PYTHON3 to a compatible interpreter." >&2
    exit 1
fi
NICE_PYTHON_VERSION=$("$NICE_PYTHON" -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
echo "NICE Python: $NICE_PYTHON ($NICE_PYTHON_VERSION)"

export TMPDIR="$HOME/tmp"
mkdir -p "$TMPDIR"
chmod -R u+rwx "$TMPDIR"
rm -rf build
mkdir build
cd build
cmake -DBUILD_TESTS=${BUILD_TESTS} -DNICE=ON \
    -DPython3_EXECUTABLE:FILEPATH="$NICE_PYTHON" ..
make -j"$JOBS" || exit $?
# Unless --exclude-tests is given, the build-time test gate (scripts/run_fast_gate.sh)
# runs between build and install: a failed gate is a failed build and
# nothing is installed; its status is the script's status.
if [ "$BUILD_TESTS" = ON ]; then "$ROOT/scripts/run_fast_gate.sh" "$PWD" || GATE_RC=$?; fi
[ "${GATE_RC:-0}" = 0 ] && { make -j"$JOBS" install || exit $?; }
exit ${GATE_RC:-0}
