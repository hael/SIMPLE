#!/bin/bash
# scripts/test_timing_run.sh — Phase 0 of doc/refactoring_notes/
# uniform_test_environment_refactoring.md: run every Fortran test on both
# routes, each in its own directory under a timeout, and record wall time,
# exit status and the tail of its output.
#
# Needs an installed --compile-tests build (build/bin holds simple_test_exec
# and the standalone simple_test_* binaries). Nothing is compiled here.
#
# usage: scripts/test_timing_run.sh [--label NAME] [--timeout SECS] [--omp N]
#                                   [--route standalone|exec|both] [--only a,b,c]
#                                   [--args FILE]
#
#   --label    run name; results go to build_test_runs/<label>/ at the repo root
#              (default: timing). Outside build/, so compile_*.sh's rm -rf build
#              does not erase them; covered by .gitignore's build* pattern.
#   --timeout  per-test wall-clock limit in seconds (default: 120; a fast-tier
#              candidate that needs more is not fast, and the extensive tier is
#              timed separately with a longer limit)
#   --omp      OMP_NUM_THREADS for every test (default: 1, the fast-tier setting)
#   --route    which route(s) to run (default: both)
#   --only     comma-separated test names to run (default: all)
#   --args     per-test argument table, "name<TAB>args" per line
#              (default: scripts/test_args.tsv if present)
#
# Output: build_test_runs/<label>/results.tsv with the columns
#   route  name  status  wall_s  exit  args
# where status is pass | fail | timeout | crash (exit >= 128) | noargs
# (the test refused for lack of required arguments), plus one log per test
# in build_test_runs/<label>/<route>_<name>/run.log. Never edits the tree.

set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/build"
BIN="$BUILD/bin"
label="timing"; tmo=120; omp=1; route="both"; only=""; argsfile="$ROOT/scripts/test_args.tsv"
while [ $# -gt 0 ]; do
    case "$1" in
        --label)   label="$2"; shift 2 ;;
        --timeout) tmo="$2"; shift 2 ;;
        --omp)     omp="$2"; shift 2 ;;
        --route)   route="$2"; shift 2 ;;
        --only)    only="$2"; shift 2 ;;
        --args)    argsfile="$2"; shift 2 ;;
        -h|--help) awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; exit 0 ;;
        *) echo "test_timing_run: unknown option $1" >&2; exit 1 ;;
    esac
done
[ -x "$BIN/simple_test_exec" ] || { echo "test_timing_run: $BIN/simple_test_exec not found; build with --compile-tests first" >&2; exit 1; }

OUT="$ROOT/build_test_runs/$label"
mkdir -p "$OUT"
RES="$OUT/results.tsv"
printf 'route\tname\tstatus\twall_s\texit\targs\n' > "$RES"

# environment the tests expect (as in CI)
export SIMPLE_PATH="$BUILD"
export SIMPLE_QSYS="${SIMPLE_QSYS:-local}"
export SIMPLE_EMAIL="${SIMPLE_EMAIL:-test@localhost}"
export PATH="$BIN:$BUILD/scripts:$PATH"
export OMP_NUM_THREADS="$omp"
export GFORTRAN_UNBUFFERED_ALL=y   # so run.log shows progress while a test runs

now() { perl -MTime::HiRes=time -e 'printf "%.3f", time'; }
wanted() {   # name -> 0 if selected
    [ -z "$only" ] && return 0
    case ",$only," in *",$1,"*) return 0 ;; esac
    return 1
}
args_for() { # name -> args from the table, if any
    [ -f "$argsfile" ] || return 0
    awk -F'\t' -v n="$1" '$1 == n { print $2; exit }' "$argsfile"
}

# run one command in its own directory under a timeout; sets STATUS EXIT WALL
run_one() {
    local dir="$1"; shift
    rm -rf "$dir"; mkdir -p "$dir"
    local t0 t1
    t0=$(now)
    ( cd "$dir" && perl -e 'my $t = shift; alarm $t; exec @ARGV' "$tmo" "$@" > run.log 2>&1 ) 2>/dev/null
    EXIT=$?
    t1=$(now)
    WALL=$(perl -e "printf '%.2f', $t1 - $t0")
    if   [ "$EXIT" -eq 142 ]; then STATUS=timeout          # SIGALRM
    elif [ "$EXIT" -eq 0 ];   then STATUS=pass
    elif [ "$EXIT" -ge 128 ]; then STATUS=crash
    elif grep -q -i -E "not defined|required|missing|usage:|ERROR! .*(defined|given|specified)" "$dir/run.log"; then STATUS=noargs
    else STATUS=fail
    fi
}

ntot=0
# ---- route 1: standalone binaries
if [ "$route" = "both" ] || [ "$route" = "standalone" ]; then
    for exe in "$BIN"/simple_test_*; do
        [ -x "$exe" ] || continue
        name=$(basename "$exe"); name=${name#simple_test_}
        [ "$name" = "exec" ] && continue
        wanted "$name" || continue
        a=$(args_for "$name")
        printf '%-12s %-45s ' standalone "$name"
        run_one "$OUT/standalone_$name" "$exe" $a  # tail -f the run.log to watch a slow one
        printf '%-8s %8s s\n' "$STATUS" "$WALL"
        printf 'standalone\t%s\t%s\t%s\t%s\t%s\n' "$name" "$STATUS" "$WALL" "$EXIT" "$a" >> "$RES"
        ntot=$((ntot + 1))
    done
fi
# ---- route 2: simple_test_exec cases, names taken from the test UI sources
if [ "$route" = "both" ] || [ "$route" = "exec" ]; then
    cases=$(perl -ne 'print "$1\n" if /call\s+add_ui_program\(\s*'"'"'([A-Za-z0-9_]+)'"'"'/' "$ROOT"/src/main/ui/simple_test/*.f90 | sort -u)
    for name in $cases; do
        wanted "$name" || continue
        a=$(args_for "$name")
        printf '%-12s %-45s ' exec "$name"
        run_one "$OUT/exec_$name" "$BIN/simple_test_exec" "test=$name" $a
        printf '%-8s %8s s\n' "$STATUS" "$WALL"
        printf 'exec\t%s\t%s\t%s\t%s\t%s\n' "$name" "$STATUS" "$WALL" "$EXIT" "$a" >> "$RES"
        ntot=$((ntot + 1))
    done
fi
echo
echo "$ntot tests run; results in $RES"
awk -F'\t' 'NR > 1 { c[$3]++ } END { for (k in c) printf "  %-8s %d\n", k, c[k] }' "$RES"
echo "slowest:"
tail -n +2 "$RES" | sort -t"$(printf '\t')" -k4,4gr | head -15 | awk -F'\t' '{ printf "  %8s s  %-8s %-11s %s\n", $4, $3, $1, $2 }'
