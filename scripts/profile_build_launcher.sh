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
t0=$(perl -MTime::HiRes=time -e 'printf "%.3f", time')
"$@"
rc=$?
t1=$(perl -MTime::HiRes=time -e 'printf "%.3f", time')
dt=$(perl -e "printf '%.3f', $t1 - $t0")
printf '%s\t%s\t%s\t%s\n' "$t0" "$t1" "$dt" "$src" >> "$log"
exit $rc
