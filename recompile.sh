#!/bin/bash
# Incremental build of an existing build/ (made by a compile_*.sh script) and install.
#
#   ./recompile.sh                 rebuild what changed; when the build has tests, the
#                                  fast test gate runs before installation (default)
#   ./recompile.sh --no-gate       rebuild and install without the gate
#   ./recompile.sh --jobs 4        limit concurrent compiler processes
#   ./recompile.sh --reconfigure   rerun cmake first (it reruns by itself when a source
#                                  file under src/ was added or removed)
#   ./recompile.sh --build-dir D   another configured build tree (or SIMPLE_BUILD_DIR)
# The build type and the test setting stay as the compile script configured them.
GATE=yes
RECONF=no
JOBS=""
ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD="${SIMPLE_BUILD_DIR:-$ROOT/build}"
source "$ROOT/scripts/simple_build_jobs.sh" || exit $?
while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-gate) GATE=no ;;
        --reconfigure) RECONF=yes ;;
        --jobs|-j)
            shift
            if [[ $# -eq 0 || -z "$1" ]]; then echo "recompile.sh: --jobs requires a positive integer" >&2; exit 1; fi
            JOBS="$1"
            ;;
        --jobs=) echo "recompile.sh: --jobs requires a positive integer" >&2; exit 1 ;;
        --jobs=*) JOBS="${1#*=}" ;;
        --build-dir)
            shift
            if [[ $# -eq 0 || -z "$1" ]]; then echo "recompile.sh: --build-dir requires a directory" >&2; exit 1; fi
            BUILD="$1"
            ;;
        -h|--help)
            echo "usage: $(basename "$0") [--no-gate] [--jobs N] [--reconfigure] [--build-dir DIR]"
            echo "Default jobs: available CPUs, or CMAKE_BUILD_PARALLEL_LEVEL when set."
            exit 0
            ;;
        *) echo "recompile.sh: unknown option: $1 (see --help)" >&2; exit 1 ;;
    esac
    shift
done
JOBS=$(simple_build_jobs "$JOBS") || exit $?
case "$BUILD" in /*) ;; *) BUILD="$PWD/$BUILD" ;; esac
if [[ ! -f "$BUILD/CMakeCache.txt" ]]; then
    echo "recompile.sh: no configured build in $BUILD; run ./compile_clean.sh or ./compile_debug.sh first" >&2
    exit 1
fi
# The library globs its sources when cmake runs, so a file added or removed under src/ is
# invisible to make: compare the source list with the one of the last configure.
SOURCES=$(cd "$ROOT/src" && find . -type f \( -name '*.[Ff][90][08]' -o -name '*.c' -o -name '*.cpp' \) | LC_ALL=C sort)
STAMP="$BUILD/recompile_sources.txt"
if [[ $RECONF = yes || ! -f "$STAMP" || "$SOURCES" != "$(cat "$STAMP")" ]]; then
    cmake -S "$ROOT" -B "$BUILD" || exit $?
    printf '%s\n' "$SOURCES" > "$STAMP"
fi
make -C "$BUILD" -j"$JOBS" || exit $?
# As in the compile scripts: a failed gate installs nothing, and its status is the script's.
if [[ $GATE = yes ]] && grep -Eq '^BUILD_TESTS:[A-Z]+=(ON|1|YES|TRUE)$' "$BUILD/CMakeCache.txt"; then
    "$ROOT/scripts/run_fast_gate.sh" "$BUILD" || GATE_RC=$?
fi
[ "${GATE_RC:-0}" = 0 ] && { make -C "$BUILD" -j"$JOBS" install || exit $?; }
exit ${GATE_RC:-0}
