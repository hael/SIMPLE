#!/usr/bin/env bash
# Clean all *_solve3D directories in the current working directory: everything is removed
# except the patterns in KEEP below and the project files (*.simple).
# Dry run unless -x is given.   usage: clean_solve3D_dirs.sh [-x]
set -euo pipefail
shopt -s nullglob

KEEP=(
  '*.simple'
  'fsc_state*_iter*.pdf'
  'RESOLUTION_STATE*_ITER*'
  'RESOLUTION_FINAL_STATE*'
  'orientations_distribution_state*.jpg'
  'orthogonal_reprojs_state*.jpg'
  'prob_posterior_iter*_part*.txt'
  'recvol_state??.mrc'
  'recvol_state??_even.mrc'
  'recvol_state??_odd.mrc'
  'recvol_state??_even_unfil.mrc'
  'recvol_state??_odd_unfil.mrc'
  'recvol_state??_even_nu_filt.mrc'
  'recvol_state??_odd_nu_filt.mrc'
  'recvol_state??_nu_filt.mrc'
  'recvol_state??_nu_locres.mrc'
  'recvol_state??_lp.mrc'
  'recvol_state??_pproc.mrc'
  'recvol_state??_pproc_mirr.mrc'
  'recvol_state??_pcg_support.txt'
  'recvol_state??_iter*_lp.mrc'
  'recvol_state??_iter*_pproc.mrc'
  'recvol_state??_stage*_lp.mrc'
  'rec_final_state??.mrc'
  'rec_final_state??_lp.mrc'
  'rec_final_state??_pproc.mrc'
  'rec_final_state??_pproc_mirr.mrc'
  'rec_final_state??_pcg_support.txt'
)

execute=0
while getopts "x" o; do
  case $o in x) execute=1;; *) echo "usage: $0 [-x]" >&2; exit 1;; esac
done

dirs=(*_solve3D/)
[ ${#dirs[@]} -gt 0 ] || { echo "no *_solve3D directories in $PWD" >&2; exit 1; }
for d in "${dirs[@]%/}"; do
  n=0; kb=0
  while IFS= read -r -d '' f; do
    name=${f##*/}
    keep=0
    for g in "${KEEP[@]}"; do [[ $name == $g ]] && { keep=1; break; }; done
    [ $keep = 1 ] && continue
    kb=$((kb + $(du -sk -- "$f" | cut -f1))); n=$((n+1))
    if [ $execute = 1 ]; then rm -rf -- "$f"; else echo "would remove: $f"; fi
  done < <(find "$d" -mindepth 1 -maxdepth 1 -print0)
  echo "$d: $n entries, $((kb/1024)) MB $([ $execute = 1 ] && echo removed || echo 'to remove (dry run, use -x)')"
done
