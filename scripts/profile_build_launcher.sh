#!/bin/bash
# Compiler launcher used by scripts/profile_build.sh (CMAKE_<LANG>_COMPILER_LAUNCHER).
# Invoked as:  profile_build_launcher.sh <logfile> <compiler> <args...>
# Appends one tab-separated line per compile:  start  end  seconds  source
# and forwards the compiler's exit status untouched.
log="$1"; shift
src=""
for a in "$@"; do
    case "$a" in
        *.f90|*.F90|*.f08|*.F08|*.f|*.F|*.c|*.cpp|*.cc|*.cu) src="$a" ;;
    esac
done
t0=$(python3 -c 'import time; print("%.3f" % time.monotonic())') || exit $?
"$@"
rc=$?
t1=$(python3 -c 'import time; print("%.3f" % time.monotonic())') || exit $?
dt=$(awk -v start="$t0" -v end="$t1" 'BEGIN {printf "%.3f", end - start}')
printf '%s\t%s\t%s\t%s\n' "$t0" "$t1" "$dt" "$src" >> "$log"
exit $rc
