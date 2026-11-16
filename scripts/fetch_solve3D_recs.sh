#!/usr/bin/env bash
# fetch_solve3D_rec.sh -- fetch rec_final_state01_pproc.mrc from every
# N_solve3D directory of a remote project, and nothing else.
#
# Usage:
#   fetch_solve3D_rec.sh [-f] [-n] [-m FILENAME] HOST:/remote/project/dir [LOCAL_DEST]
#
#   HOST:/remote/project/dir  directory containing 1_solve3D ... N_solve3D
#   LOCAL_DEST                local target (default: current directory)
#
# Options:
#   -f           flatten: save as LOCAL_DEST/N_solve3D_rec_final_state01_pproc.mrc
#                instead of keeping LOCAL_DEST/N_solve3D/rec_final_state01_pproc.mrc
#   -n           dry run (list what would be transferred)
#   -m FILENAME  fetch a different file name (default rec_final_state01_pproc.mrc)
#
# Examples:
#   fetch_solve3D_rec.sh dell:/data/proj/ab_initio ./recs
#   fetch_solve3D_rec.sh -f dell:/data/proj/ab_initio ./recs
#
# Re-running only transfers new/changed files (rsync).

set -euo pipefail

FILE="rec_final_state01_pproc.mrc"
FLAT=0
DRY=()

usage() { sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-1}"; }

while getopts "fnm:h" opt; do
    case "$opt" in
        f) FLAT=1 ;;
        n) DRY=(--dry-run) ;;
        m) FILE="$OPTARG" ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done
shift $((OPTIND - 1))

[ $# -ge 1 ] && [ $# -le 2 ] || usage 1
SRC="${1%/}"
DEST="${2:-.}"
case "$SRC" in
    *:*) ;;
    *) echo "error: source must be HOST:/path" >&2; exit 1 ;;
esac

mkdir -p "$DEST"

# Only [0-9]*_solve3D directories, only the one file, prune everything else.
# -m drops directories that end up empty (e.g. a run that has no such file).
rsync -avhm --progress ${DRY[@]+"${DRY[@]}"} \
    --include='[0-9]*_solve3D/' \
    --include="$FILE" \
    --exclude='*' \
    "$SRC/" "$DEST/"

if [ "$FLAT" -eq 1 ] && [ "${#DRY[@]}" -eq 0 ]; then
    for d in "$DEST"/[0-9]*_solve3D; do
        [ -f "$d/$FILE" ] || continue
        mv -f "$d/$FILE" "$DEST/$(basename "$d")_$FILE"
        rmdir "$d" 2>/dev/null || true
    done
fi