#!/usr/bin/env bash
set -euo pipefail

# ============ How to use the script: ============
# Expect MSYS UCRT64 environment inside Windows
# Run this script from an MSYS2 UCRT64 shell.  
#
# export PATH=/ucrt64/bin:/usr/bin:$PATH
# cd ~/SIMPLE
#
# Keep build products outside the source checkout so generated modules and CMake metadata cannot pollute Git.
# Must Keep build products isolated from the source files.
src_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${SIMPLE_BUILD_DIR:-${src_dir}/../SIMPLE_build_debug}"
build_type="${SIMPLE_BUILD_TYPE:-Debug}"
requested_jobs="${SIMPLE_BUILD_JOBS:-12}"
clean_build="${SIMPLE_CLEAN_BUILD:-no}"

# ============ For a Debug build: ==============
# SIMPLE_BUILD_DIR=~/SIMPLE_build_debug SIMPLE_BUILD_TYPE=Debug SIMPLE_BUILD_JOBS=12 ./compile_windows.sh
#
# ============ For a Release build: ============
# SIMPLE_BUILD_DIR=~/SIMPLE_build_release SIMPLE_BUILD_TYPE=Release SIMPLE_BUILD_JOBS=12 ./compile_windows.sh
#
# After completion, executables should be under: ============
# export SIMPLE_EMAIL="my.name@uni.edu"
# export SIMPLE_QSYS="local"
# export SIMPLE_PATH=~/SIMPLE_build_debug
# export PATH=${SIMPLE_PATH}/scripts:${SIMPLE_PATH}/bin:${PATH}
#
# Test a local project: ============
# mkdir -p "$HOME/betagal"
# simple_exec prg=new_project projname=betagal dir="$HOME/betagal"

if [[ "$build_dir" != /* ]]; then
    build_dir="${src_dir}/${build_dir}"
fi
src_dir="$(realpath -m "$src_dir")"
build_dir="$(realpath -m "$build_dir")"

case "$build_dir" in
    "$src_dir"|"/")
        printf 'Refusing unsafe build directory: %s\n' "$build_dir" >&2
        exit 2
        ;;
    "$src_dir"/build*)
        ;;
    "$src_dir"/*)
        printf 'Build directory must be outside the source tree: %s\n' "$build_dir" >&2
        exit 2
        ;;
esac

if ! [[ "$requested_jobs" =~ ^[1-9][0-9]*$ ]]; then
    printf 'SIMPLE_BUILD_JOBS must be a positive integer: %s\n' "$requested_jobs" >&2
    exit 2
fi

case "$clean_build" in
    yes|no) ;;
    *)
        printf 'SIMPLE_CLEAN_BUILD must be yes or no: %s\n' "$clean_build" >&2
        exit 2
        ;;
esac

available_jobs="$(nproc)"
build_jobs="$requested_jobs"
if (( build_jobs > available_jobs )); then
    build_jobs="$available_jobs"
fi

printf 'Source:      %s\n' "$src_dir"
printf 'Build:       %s\n' "$build_dir"
printf 'Build type:  %s\n' "$build_type"
printf 'Build jobs:  %s\n' "$build_jobs"
printf 'Clean build: %s\n' "$clean_build"

if [[ "$clean_build" == yes ]]; then
    rm --recursive --force -- "$build_dir"
fi
mkdir -p "$build_dir"

cmake -S "$src_dir" -B "$build_dir" \
    -G "MSYS Makefiles" \
    -DCMAKE_BUILD_TYPE="$build_type" \
    -DCMAKE_INSTALL_PREFIX="$build_dir"

cmake --build "$build_dir" --parallel "$build_jobs"
cmake --install "$build_dir"
