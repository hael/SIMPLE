#!/bin/bash
# Clean release build and install into build/.
#
#   ./compile_clean.sh                  library, executables and tests; the fast
#                                       test gate runs before installation (default)
#   ./compile_clean.sh --exclude-tests  library and executables only, no gate
#   ./compile_clean.sh --jobs 4         limit concurrent compiler processes
# High-level tests are registered but never run here; invoke them explicitly
# after building with: (cd build && ctest -L highlevel --output-on-failure)
#
# The test code is about a fifth of the summed compile time; --exclude-tests skips
# it when only the executables are needed (doc/policies/compile_time_policy.md).
BUILD_TESTS=ON
JOBS=""
ROOT="$(cd "$(dirname "$0")" && pwd)"
source "$ROOT/scripts/simple_build_jobs.sh" || exit $?
while [[ $# -gt 0 ]]; do
    case "$1" in
        --exclude-tests) BUILD_TESTS=OFF ;;
        --jobs|-j)
            shift
            if [[ $# -eq 0 || -z "$1" ]]; then echo "compile_clean.sh: --jobs requires a positive integer" >&2; exit 1; fi
            JOBS="$1"
            ;;
        --jobs=) echo "compile_clean.sh: --jobs requires a positive integer" >&2; exit 1 ;;
        --jobs=*) JOBS="${1#*=}" ;;
        -h|--help)
            echo "usage: $(basename "$0") [--exclude-tests] [--jobs N]"
            echo "Default jobs: available CPUs, or CMAKE_BUILD_PARALLEL_LEVEL when set."
            exit 0
            ;;
        *) echo "compile_clean.sh: unknown option: $1 (see --help)" >&2; exit 1 ;;
    esac
    shift
done
JOBS=$(simple_build_jobs "$JOBS") || exit $?
rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=${BUILD_TESTS} || exit $?
make -j"$JOBS" || exit $?
# Unless --exclude-tests is given, the build-time test gate (scripts/run_fast_gate.sh)
# runs between build and install: a failed gate is a failed build and
# nothing is installed; its status is the script's status.
if [ "$BUILD_TESTS" = ON ]; then "$ROOT/scripts/run_fast_gate.sh" "$PWD" || GATE_RC=$?; fi
[ "${GATE_RC:-0}" = 0 ] && { make -j"$JOBS" install || exit $?; }
exit ${GATE_RC:-0}
