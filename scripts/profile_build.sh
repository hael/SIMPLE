#!/bin/bash
# scripts/profile_build.sh — compile-time profiling harness for SIMPLE (Make generator).
#
# Records per-file compile times through a compiler launcher, plus the
# whole-build numbers listed in doc/refactoring_notes/planned/
# contract_submodule_architecture.md section 8. Works on macOS (bash 3.2) and
# Linux. Never runs `make install` (two worktrees would clobber one prefix).
#
# Usage (run from anywhere; the tree is located from the script's own path):
#
#   scripts/profile_build.sh clean  [--label NAME] [--debug] [--no-tests] [-j N]
#       rm -rf build; cmake; make -jN. Full clean-build profile.
#
#   scripts/profile_build.sh touch  <file> [--label NAME] [-j N]
#       touch <file> in an existing profiled build tree, make -jN, and report
#       exactly which sources recompiled and how long it took.
#
#   scripts/profile_build.sh probe  <file> [--label NAME] [-j N]
#       Interface change: insert a public named constant into the module that
#       <file> declares (after its first `implicit none`), rebuild, and report
#       exactly which sources recompiled; then restore <file> and rebuild
#       quietly so the tree is back in its original state. A plain `touch`
#       never changes a .mod and so never cascades; this does.
#
#   scripts/profile_build.sh report <profile_dir>
#       Recompute summary.txt from a profile directory's compile.tsv.
#
#   scripts/profile_build.sh compare <profile_dir_A> <profile_dir_B>
#       Per-file and total deltas between two profiles (typically master vs branch).
#
# Every run writes build_profile/<label>_<timestamp>/ in the tree root with:
#   meta.txt            commit, branch, host, compiler, -j, wall time
#   make.log            full make output
#   compile.tsv         start  end  seconds  source   (one line per compile job)
#   compile_sorted.tsv  same, slowest first, source paths relative to the tree
#   summary.txt         the headline numbers
# build_profile/ is covered by the `build*` pattern in .gitignore.

set -u

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LAUNCHER="$ROOT/scripts/profile_build_launcher.sh"
BUILD="$ROOT/build"
LIVELOG="$BUILD/build_profile_compile.tsv"   # fixed path baked into the build tree at configure
PROFDIR_BASE="$ROOT/build_profile"

now() { perl -MTime::HiRes=time -e 'printf "%.3f", time'; }
ncpu() {
    if command -v nproc >/dev/null 2>&1; then nproc
    else sysctl -n hw.ncpu
    fi
}
die() { echo "build_profile: $*" >&2; exit 1; }

# ------------------------------------------------------------------------------
# summary: compile.tsv + build dir -> summary.txt
# ------------------------------------------------------------------------------
write_summary() {
    local prof="$1"
    local tsv="$prof/compile.tsv"
    [ -f "$tsv" ] || die "no compile.tsv in $prof"
    awk -F'\t' -v root="$ROOT/" 'BEGIN{OFS="\t"} { s=$4; if (index(s,root)==1) s=substr(s,length(root)+1); print $1,$2,$3,s }' "$tsv" \
        | sort -t"$(printf '\t')" -k3,3gr > "$prof/compile_sorted.tsv"
    local njobs cpu
    njobs=$(wc -l < "$tsv" | tr -d ' ')
    cpu=$(awk -F'\t' '{s+=$3} END{printf "%.1f", s}' "$tsv")
    {
        echo "== build_profile summary =="
        [ -f "$prof/meta.txt" ] && cat "$prof/meta.txt"
        echo
        echo "compile jobs:            $njobs"
        echo "compile CPU (sum, s):    $cpu"
        if [ -d "$BUILD/modules" ]; then
            echo "modules dir size (KB):   $(du -sk "$BUILD/modules" | awk '{print $1}')"
            echo "modules dir files:       $(find "$BUILD/modules" -type f | wc -l | tr -d ' ')"
            echo ".mod files:              $(find "$BUILD/modules" -name '*.mod' | wc -l | tr -d ' ')"
            echo ".smod files:             $(find "$BUILD/modules" -name '*.smod' | wc -l | tr -d ' ')"
            echo
            echo "largest .mod files (bytes):"
            find "$BUILD/modules" -name '*.mod' -exec wc -c {} + | grep -v ' total$' | sort -k1,1nr | head -10 \
                | awk -v root="$BUILD/" '{ p=$2; if (index(p,root)==1) p=substr(p,length(root)+1); printf "  %10d  %s\n", $1, p }'
        fi
        echo
        echo "slowest compiles (s):"
        head -25 "$prof/compile_sorted.tsv" | awk -F'\t' '{ printf "  %8.2f  %s\n", $3, $4 }'
        echo
        echo "all recompiled sources are listed in compile_sorted.tsv"
    } > "$prof/summary.txt"
    cat "$prof/summary.txt"
}

write_meta() {
    local prof="$1" mode="$2" wall="$3" jobs="$4" extra="$5"
    local fc fcver=""
    fc=$(grep '^CMAKE_Fortran_COMPILER:' "$BUILD/CMakeCache.txt" 2>/dev/null | cut -d= -f2)
    [ -n "$fc" ] && fcver=$("$fc" --version 2>/dev/null | head -1)
    {
        echo "mode:        $mode"
        echo "tree:        $ROOT"
        echo "commit:      $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null) ($(git -C "$ROOT" rev-parse --abbrev-ref HEAD 2>/dev/null))"
        echo "dirty files: $(git -C "$ROOT" status --porcelain 2>/dev/null | grep -v '^??' | wc -l | tr -d ' ')"
        echo "host:        $(hostname)  ($(uname -sm))"
        echo "date:        $(date '+%Y-%m-%d %H:%M:%S')"
        echo "compiler:    $fcver"
        echo "build type:  $(grep '^CMAKE_BUILD_TYPE:' "$BUILD/CMakeCache.txt" 2>/dev/null | cut -d= -f2)"
        echo "build tests: $(grep '^BUILD_TESTS:' "$BUILD/CMakeCache.txt" 2>/dev/null | cut -d= -f2)"
        echo "make -j:     $jobs"
        [ -n "$extra" ] && echo "$extra"
        echo "wall time (s): $wall"
    } > "$prof/meta.txt"
}

new_profdir() {
    local d="$PROFDIR_BASE/${1}_$(date '+%Y%m%d_%H%M%S')"
    mkdir -p "$d" || die "cannot create $d"
    echo "$d"
}

# run make, tee the full log, show progress lines only, return make's status
timed_make() {
    local jobs="$1" prof="$2" filter="$3"
    local t0 t1
    t0=$(now)
    ( cd "$BUILD" && { make -j"$jobs" 2>&1; echo $? > "$prof/make.rc"; } | tee "$prof/make.log" | grep -E --line-buffered "$filter" )
    t1=$(now)
    WALL=$(perl -e "printf '%.1f', $t1 - $t0")
    RC=$(cat "$prof/make.rc" 2>/dev/null || echo 1)
}

# ------------------------------------------------------------------------------
# clean
# ------------------------------------------------------------------------------
cmd_clean() {
    local label="clean" jobs="" debug=0 notests=0
    while [ $# -gt 0 ]; do
        case "$1" in
            --label)    label="$2"; shift 2 ;;
            --debug)    debug=1; shift ;;
            --no-tests) notests=1; shift ;;
            -j)         jobs="$2"; shift 2 ;;
            -j*)        jobs="${1#-j}"; shift ;;
            *) die "unknown option for clean: $1" ;;
        esac
    done
    [ -n "$jobs" ] || jobs=$(ncpu)
    local prof; prof=$(new_profdir "$label")

    echo "build_profile: clean build of $ROOT -> $prof (make -j$jobs)"
    rm -rf "$BUILD"
    mkdir -p "$BUILD"
    : > "$LIVELOG"
    ( cd "$BUILD" && cmake .. \
        "-DCMAKE_Fortran_COMPILER_LAUNCHER=$LAUNCHER;$LIVELOG" \
        "-DCMAKE_C_COMPILER_LAUNCHER=$LAUNCHER;$LIVELOG" \
        "-DCMAKE_CXX_COMPILER_LAUNCHER=$LAUNCHER;$LIVELOG" \
        $( [ $debug -eq 1 ] && echo -DCMAKE_BUILD_TYPE=debug ) \
        $( [ $notests -eq 1 ] && echo -DBUILD_TESTS=OFF ) \
        > "$prof/cmake.log" 2>&1 ) || { cat "$prof/cmake.log"; die "cmake failed"; }
    : > "$LIVELOG"    # drop anything the configure probes logged

    timed_make "$jobs" "$prof" '^\[|[Ee]rror'
    cp "$LIVELOG" "$prof/compile.tsv"
    write_meta "$prof" "clean" "$WALL" "$jobs" ""
    write_summary "$prof"
    echo
    echo "profile written to: $prof"
    [ "$RC" = "0" ] || echo "build_profile: WARNING make exited with status $RC (see make.log)"
    return "$RC"
}

# ------------------------------------------------------------------------------
# touch
# ------------------------------------------------------------------------------
cmd_touch() {
    [ $# -ge 1 ] || die "touch needs a file"
    local file="$1"; shift
    local label="" jobs=""
    while [ $# -gt 0 ]; do
        case "$1" in
            --label) label="$2"; shift 2 ;;
            -j)      jobs="$2"; shift 2 ;;
            -j*)     jobs="${1#-j}"; shift ;;
            *) die "unknown option for touch: $1" ;;
        esac
    done
    [ -n "$jobs" ] || jobs=$(ncpu)
    [ -f "$file" ] || file="$ROOT/$file"
    [ -f "$file" ] || die "no such file: $1"
    [ -f "$BUILD/CMakeCache.txt" ] || die "no build tree; run 'clean' first"
    grep -q "CMAKE_Fortran_COMPILER_LAUNCHER.*profile_build_launcher" "$BUILD/CMakeCache.txt" \
        || die "build tree was not configured by this script; run 'clean' first"
    [ -n "$label" ] || label="touch_$(basename "$file" | sed 's/\.[fF]90$//')"
    local prof; prof=$(new_profdir "$label")

    # bring the tree up to date first so only the touch is measured
    ( cd "$BUILD" && make -j"$jobs" > "$prof/make_pre.log" 2>&1 ) || die "pre-build failed (see $prof/make_pre.log)"
    : > "$LIVELOG"

    echo "build_profile: touch ${file#$ROOT/} -> $prof (make -j$jobs)"
    touch "$file"
    timed_make "$jobs" "$prof" 'Building|Linking|[Ee]rror'
    cp "$LIVELOG" "$prof/compile.tsv"
    write_meta "$prof" "touch" "$WALL" "$jobs" "touched:     ${file#$ROOT/}"
    write_summary "$prof"
    echo
    echo "profile written to: $prof"
    return "$RC"
}

# ------------------------------------------------------------------------------
# probe (interface change)
# ------------------------------------------------------------------------------
cmd_probe() {
    [ $# -ge 1 ] || die "probe needs a file"
    local file="$1"; shift
    local label="" jobs=""
    while [ $# -gt 0 ]; do
        case "$1" in
            --label) label="$2"; shift 2 ;;
            -j)      jobs="$2"; shift 2 ;;
            -j*)     jobs="${1#-j}"; shift ;;
            *) die "unknown option for probe: $1" ;;
        esac
    done
    [ -n "$jobs" ] || jobs=$(ncpu)
    [ -f "$file" ] || file="$ROOT/$file"
    [ -f "$file" ] || die "no such file: $1"
    [ -f "$BUILD/CMakeCache.txt" ] || die "no build tree; run 'clean' first"
    grep -q "CMAKE_Fortran_COMPILER_LAUNCHER.*profile_build_launcher" "$BUILD/CMakeCache.txt" \
        || die "build tree was not configured by this script; run 'clean' first"
    grep -qi '^[[:space:]]*implicit[[:space:]]*none' "$file" || die "no implicit none in $file"
    [ -n "$label" ] || label="probe_$(basename "$file" | sed 's/\.[fF]90$//')"
    local prof; prof=$(new_profdir "$label")

    ( cd "$BUILD" && make -j"$jobs" > "$prof/make_pre.log" 2>&1 ) || die "pre-build failed (see $prof/make_pre.log)"
    cp -p "$file" "$prof/probe_original.f90"
    perl -i -pe 'if (!$done && /^\s*implicit\s+none/i) { $_ .= "integer, parameter, public :: PROFILE_BUILD_PROBE = 1\n"; $done = 1 }' "$file"
    : > "$LIVELOG"

    echo "build_profile: probe ${file#$ROOT/} -> $prof (make -j$jobs)"
    timed_make "$jobs" "$prof" 'Building|Linking|[Ee]rror'
    cp "$LIVELOG" "$prof/compile.tsv"
    cp -p "$prof/probe_original.f90" "$file"
    touch "$file"
    echo "build_profile: restoring ${file#$ROOT/} and rebuilding"
    ( cd "$BUILD" && make -j"$jobs" > "$prof/make_restore.log" 2>&1 ) || echo "build_profile: WARNING restore build failed (see $prof/make_restore.log)"
    write_meta "$prof" "probe" "$WALL" "$jobs" "probed:      ${file#$ROOT/}"
    write_summary "$prof"
    echo
    echo "profile written to: $prof"
    return "$RC"
}

# ------------------------------------------------------------------------------
# compare
# ------------------------------------------------------------------------------
cmd_compare() {
    [ $# -eq 2 ] || die "compare needs two profile directories"
    local a="$1" b="$2"
    [ -f "$a/compile_sorted.tsv" ] || die "$a has no compile_sorted.tsv (run report on it)"
    [ -f "$b/compile_sorted.tsv" ] || die "$b has no compile_sorted.tsv (run report on it)"
    local out="$b/compare_vs_$(basename "$a").tsv"
    echo "A: $a"; grep -E '^(tree|commit|wall time|make -j|compile)' "$a/meta.txt" | sed 's/^/   /'
    echo "B: $b"; grep -E '^(tree|commit|wall time|make -j|compile)' "$b/meta.txt" | sed 's/^/   /'
    echo
    awk -F'\t' -v A="$a/compile_sorted.tsv" -v B="$b/compile_sorted.tsv" -v out="$out" '
        BEGIN{OFS="\t"}
        FILENAME==A { ta[$4]+=$3; sa+=$3; na++; next }
        FILENAME==B { tb[$4]+=$3; sb+=$3; nb++; next }
        END {
            printf "compile jobs:   A %d   B %d\n", na, nb
            printf "compile CPU:    A %.1f s   B %.1f s   (B-A %+.1f s, %+.1f%%)\n", sa, sb, sb-sa, (sa>0 ? 100*(sb-sa)/sa : 0)
            for (f in ta) seen[f]=1
            for (f in tb) seen[f]=1
            for (f in seen) printf "%.2f\t%.2f\t%.2f\t%s\n", tb[f]-ta[f], ta[f], tb[f], f > out
        }' "$a/compile_sorted.tsv" "$b/compile_sorted.tsv"
    local tmp="$out.sorted"
    sort -t"$(printf '\t')" -k1,1g "$out" > "$tmp" && mv "$tmp" "$out"
    echo
    echo "per-file delta (B-A s, A s, B s, file) — largest speed-ups:"
    head -20 "$out" | awk -F'\t' '{ printf "  %+8.2f  %8.2f  %8.2f  %s\n", $1, $2, $3, $4 }'
    echo "  ..."
    echo "largest slow-downs / new files:"
    tail -20 "$out" | awk -F'\t' '{ printf "  %+8.2f  %8.2f  %8.2f  %s\n", $1, $2, $3, $4 }'
    echo
    echo "full per-file table: $out"
}

# ------------------------------------------------------------------------------
WALL=""; RC=1
case "${1:-}" in
    clean)   shift; cmd_clean "$@" ;;
    touch)   shift; cmd_touch "$@" ;;
    probe)   shift; cmd_probe "$@" ;;
    report)  shift; [ $# -eq 1 ] || die "report needs a profile directory"; write_summary "$1" ;;
    compare) shift; cmd_compare "$@" ;;
    *)  awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; exit 1 ;;
esac
