#!/usr/bin/env bash

simple_build_jobs() {
    local jobs="${1:-${CMAKE_BUILD_PARALLEL_LEVEL:-}}"
    if [[ -z "$jobs" ]]; then
        jobs=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc 2>/dev/null || \
            getconf _NPROCESSORS_ONLN 2>/dev/null || \
            sysctl -n hw.ncpu 2>/dev/null || printf '%s\n' "${NUMBER_OF_PROCESSORS:-1}")
    fi
    if [[ ! "$jobs" =~ ^[1-9][0-9]*$ ]]; then
        printf '%s: jobs must be a positive integer (got %s)\n' "${0##*/}" "$jobs" >&2
        return 1
    fi
    printf '%s\n' "$jobs"
}