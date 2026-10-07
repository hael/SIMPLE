#!/bin/bash
# Rerun tests from an existing SIMPLE build without configuring or compiling.
#
# The default is the same registry-checked, budgeted fast gate used by the
# compile scripts. Set SIMPLE_BUILD_DIR or pass --build-dir to use another
# configured build tree.
set -eu

ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD="${SIMPLE_BUILD_DIR:-$ROOT/build}"

usage() {
    cat <<EOF
usage: $(basename "$0") [--build-dir DIR] [MODE|TEST] [CTEST_OPTIONS...]

Rerun tests from an existing build; this script never configures or compiles.

  (no mode), fast  run the registry-checked, 30-second fast gate
  library          run tests labelled library
  highlevel        run tests labelled highlevel
  platform         run tests labelled platform
  all              run every registered CTest entry
  list             list registered CTest entries without running them
  ctest ARGS...    pass ARGS directly to CTest
  TEST             run one exact CTest entry, for example unit_heterogeneity

The build directory defaults to ./build. SIMPLE_BUILD_DIR provides the same
override as --build-dir.
EOF
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --build-dir)
            [ "$#" -ge 2 ] || { echo "run_tests: --build-dir requires a directory" >&2; exit 2; }
            BUILD="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            break
            ;;
    esac
done

case "$BUILD" in
    /*) ;;
    *) BUILD="$ROOT/$BUILD" ;;
esac

[ -f "$BUILD/CTestTestfile.cmake" ] || {
    echo "run_tests: no CTest configuration in $BUILD" >&2
    echo "run_tests: compile once with tests enabled, then rerun this script" >&2
    exit 2
}
command -v ctest >/dev/null 2>&1 || {
    echo "run_tests: ctest is not available on PATH" >&2
    exit 2
}

mode="${1:-fast}"
if [ "$#" -gt 0 ]; then shift; fi

case "$mode" in
    fast)
        [ "$#" -eq 0 ] || { echo "run_tests: fast mode does not accept CTest options; use 'ctest' mode" >&2; exit 2; }
        exec "$ROOT/scripts/run_fast_gate.sh" "$BUILD"
        ;;
    library|highlevel|platform)
        cd "$BUILD"
        exec ctest -L "$mode" --no-tests=error --output-on-failure "$@"
        ;;
    all)
        cd "$BUILD"
        exec ctest --no-tests=error --output-on-failure "$@"
        ;;
    list)
        cd "$BUILD"
        exec ctest -N "$@"
        ;;
    ctest)
        cd "$BUILD"
        exec ctest "$@"
        ;;
    *)
        cd "$BUILD"
        exec ctest -R "^${mode}$" --no-tests=error --output-on-failure "$@"
        ;;
esac
