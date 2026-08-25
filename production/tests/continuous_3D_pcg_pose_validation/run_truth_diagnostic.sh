#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_truth_diagnostic.sh \
    --truth-volume /absolute/path/to/truth.mrc \
    --truth-oris /absolute/path/to/truth_oris.txt \
    --noisy-stack /absolute/path/to/noisy.mrcs \
    --box 144 --smpd 1.3 --mskdiam 120 --pgrp c1 \
    --kv 300 --cs 2.7 --fraca 0.1 \
    --lp-fsc05 4.255 --lp-fsc0143 3.671 \
    [--nthr 64] [--output-parent "$PWD"]

The runner executes all 16 matched arms after preflight. It collects clean and
noisy, exact and perturbed, and full/FSC-limited evidence in one timestamped
directory. An arm failure is recorded without stopping later independent arms.
EOF
}

TRUTH_VOLUME=""
TRUTH_ORIS=""
NOISY_STACK=""
BOX=""
SMPD=""
MSKDIAM=""
PGRP=""
KV=""
CS=""
FRACA=""
LP_FSC05=""
LP_FSC0143=""
NTHR=64
OUTPUT_PARENT="$PWD"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --truth-volume) TRUTH_VOLUME=$2; shift 2 ;;
        --truth-oris) TRUTH_ORIS=$2; shift 2 ;;
        --noisy-stack) NOISY_STACK=$2; shift 2 ;;
        --box) BOX=$2; shift 2 ;;
        --smpd) SMPD=$2; shift 2 ;;
        --mskdiam) MSKDIAM=$2; shift 2 ;;
        --pgrp) PGRP=$2; shift 2 ;;
        --kv) KV=$2; shift 2 ;;
        --cs) CS=$2; shift 2 ;;
        --fraca) FRACA=$2; shift 2 ;;
        --lp-fsc05) LP_FSC05=$2; shift 2 ;;
        --lp-fsc0143) LP_FSC0143=$2; shift 2 ;;
        --nthr) NTHR=$2; shift 2 ;;
        --output-parent) OUTPUT_PARENT=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

for value_name in TRUTH_VOLUME TRUTH_ORIS NOISY_STACK BOX SMPD MSKDIAM PGRP KV CS FRACA LP_FSC05 LP_FSC0143; do
    if [[ -z ${!value_name} ]]; then
        echo "ERROR: required input is missing: $value_name" >&2
        usage >&2
        exit 2
    fi
done
for input_file in "$TRUTH_VOLUME" "$TRUTH_ORIS" "$NOISY_STACK"; do
    if [[ ! -f $input_file ]]; then
        echo "ERROR: input file does not exist: $input_file" >&2
        exit 2
    fi
done

# Require the Python syntax level used by the metadata and result analyzers.
python_supported() {
    command -v python3 >/dev/null 2>&1 && \
        python3 -c 'import sys; raise SystemExit(sys.version_info < (3, 7))'
}
if ! python_supported && type module >/dev/null 2>&1; then
    echo "Loading miniconda/25.9.1 for the diagnostic analyzer"
    module load miniconda/25.9.1
fi
for executable in simple_exec simple_private_exec python3; do
    if ! command -v "$executable" >/dev/null 2>&1; then
        echo "ERROR: required executable is not on PATH: $executable" >&2
        exit 2
    fi
done
if ! python_supported; then
    echo "ERROR: Python 3.7 or newer is required" >&2
    exit 2
fi

TRUTH_VOLUME=$(readlink -f "$TRUTH_VOLUME")
TRUTH_ORIS=$(readlink -f "$TRUTH_ORIS")
NOISY_STACK=$(readlink -f "$NOISY_STACK")
mkdir -p "$OUTPUT_PARENT"
OUTPUT_PARENT=$(readlink -f "$OUTPUT_PARENT")
RESULT_ROOT="$OUTPUT_PARENT/continuous_3D_pose_truth_diagnostic_$(date +%Y%m%d_%H%M%S)"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE_ROOT=$(readlink -f "$SCRIPT_DIR/../../..")
mkdir -p "$RESULT_ROOT"/{input,generated,bases,arms,analysis}
mkdir -p "$RESULT_ROOT/input/validation_source"
cp "$SCRIPT_DIR/run_truth_diagnostic.sh" "$RESULT_ROOT/input/validation_source/"
cp "$SCRIPT_DIR/analyze_truth_matrix.py" "$RESULT_ROOT/input/validation_source/"
cp "$SCRIPT_DIR/analyze_pose_ab.py" "$RESULT_ROOT/input/validation_source/"
cp "$SCRIPT_DIR/prepare_truth_oritab.py" "$RESULT_ROOT/input/validation_source/"
cp "$SCRIPT_DIR/README.md" "$RESULT_ROOT/input/validation_source/"
# Capture the current draft contract and consolidated historical handoff with the evidence.
cp "$SOURCE_ROOT/doc/implementation_notes/continuous_3D_pose_end_polishing_spec.md" \
    "$RESULT_ROOT/input/"
cp "$SOURCE_ROOT/doc/implementation_notes/continuous_3D_pose_end_polishing_plan.md" \
    "$RESULT_ROOT/input/"
cp "$SOURCE_ROOT/doc/implementation_notes/completed/continuous_3D_pose_end_polishing_history_and_handoff.md" \
    "$RESULT_ROOT/input/"
cp "$SOURCE_ROOT/doc/implementation_notes/continuous_3D_refinement_on_pcg_operator.md" \
    "$RESULT_ROOT/input/"
cp "$SOURCE_ROOT/doc/policies/reconstruct3D_pcg_policy.md" "$RESULT_ROOT/input/"
exec > >(tee -a "$RESULT_ROOT/diagnostic.log") 2>&1
echo "Diagnostic directory: $RESULT_ROOT"

FAILURE_COUNT=0
FAILURE_LOG="$RESULT_ROOT/failures.txt"
set +e

record_failure() {
    FAILURE_COUNT=$((FAILURE_COUNT+1))
    printf 'FAIL [%02d] %s\n' "$FAILURE_COUNT" "$*" | tee -a "$FAILURE_LOG" >&2
}

# Package status, checksums, and the complete log after all scheduled arms finish.
finalize() {
    local status=$?
    set +e
    if [[ ! -f $RESULT_ROOT/STATUS.txt ]]; then
        echo FAIL >"$RESULT_ROOT/STATUS.txt"
    fi
    (cd "$RESULT_ROOT" && find . -type f ! -name MANIFEST.sha256 ! -name diagnostic.log -print0 | \
        sort -z | xargs -0 sha256sum) >"$RESULT_ROOT/MANIFEST.sha256"
    if (( status != 0 )); then
        echo "Diagnostic stopped with status $status. Download: $RESULT_ROOT" >&2
    fi
}
trap finalize EXIT

# Run one required command and record a failure without hiding later independent arms.
run_success() {
    local workdir=$1
    local logfile=$2
    shift 2
    printf '%q ' "$@" >"${logfile%.log}.command.txt"
    printf '\n' >>"${logfile%.log}.command.txt"
    (cd "$workdir" && "$@") 2>&1 | tee "$logfile"
    local status=${PIPESTATUS[0]}
    if (( status != 0 )); then
        record_failure "command returned status $status; see $logfile"
        return 1
    fi
    if grep -Eqi 'Missing key on command line|not enough input variables defined' "$logfile"; then
        record_failure "command stopped during UI input validation; see $logfile"
        return 1
    fi
    return 0
}

# Export the final particle poses used by the truth-error analyzer.
export_poses() {
    local project=$1
    local output=$2
    simple_private_exec prg=print_project_vals projfile="$project" oritype=ptcl3D \
        keys=e1,e2,e3,x,y,eo >"$output" 2>&1
}

# Record exact source hashes so the diagnostic is tied to one implementation.
write_source_manifest() {
    local list="$RESULT_ROOT/input/source_files.txt"
    local manifest="$RESULT_ROOT/input/source_manifest.sha256"
    cat >"$list" <<EOF
src/main/commanders/simple/simple_commanders_rec.f90
src/main/strategies/search/simple_pcg_pose_polisher.f90
src/main/volume/simple_reconstructor_pcg.f90
production/tests/continuous_3D_pcg_pose_validation/run_truth_diagnostic.sh
production/tests/continuous_3D_pcg_pose_validation/analyze_truth_matrix.py
production/tests/continuous_3D_pcg_pose_validation/analyze_pose_ab.py
production/tests/continuous_3D_pcg_pose_validation/prepare_truth_oritab.py
doc/implementation_notes/continuous_3D_pose_end_polishing_spec.md
doc/implementation_notes/continuous_3D_pose_end_polishing_plan.md
doc/implementation_notes/completed/continuous_3D_pose_end_polishing_history_and_handoff.md
doc/implementation_notes/continuous_3D_refinement_on_pcg_operator.md
doc/policies/reconstruct3D_pcg_policy.md
EOF
    : >"$manifest"
    while IFS= read -r relative; do
        if [[ -f $SOURCE_ROOT/$relative ]]; then
            sha256sum "$SOURCE_ROOT/$relative" >>"$manifest"
        else
            printf 'MISSING  %s\n' "$SOURCE_ROOT/$relative" >>"$manifest"
        fi
    done <"$list"
}

NPTCLS=$(grep -cve '^[[:space:]]*$' "$TRUTH_ORIS")
if (( NPTCLS < 1 )); then
    echo "ERROR: truth orientation table has no rows" >&2
    exit 2
fi
cat >"$RESULT_ROOT/run_config.env" <<EOF
TRUTH_VOLUME=$TRUTH_VOLUME
TRUTH_ORIS=$TRUTH_ORIS
NOISY_STACK=$NOISY_STACK
NPTCLS=$NPTCLS
BOX=$BOX
SMPD=$SMPD
MSKDIAM=$MSKDIAM
PGRP=$PGRP
KV=$KV
CS=$CS
FRACA=$FRACA
LP_FSC05=$LP_FSC05
LP_FSC0143=$LP_FSC0143
NTHR=$NTHR
EOF
{
    date --iso-8601=seconds
    uname -a
    python3 --version
    command -v simple_exec
    command -v simple_private_exec
} >"$RESULT_ROOT/environment.txt"
git -C "$SOURCE_ROOT" rev-parse HEAD >"$RESULT_ROOT/input/source_commit.txt" 2>/dev/null || true
write_source_manifest
sha256sum "$TRUTH_VOLUME" "$TRUTH_ORIS" "$NOISY_STACK" >"$RESULT_ROOT/input/fixture_manifest.sha256"

PERTURBED_ORIS="$RESULT_ROOT/generated/perturbed_oris.txt"
EXACT_ORIS="$RESULT_ROOT/generated/exact_oris.txt"
if ! run_success "$RESULT_ROOT/generated" "$RESULT_ROOT/generated/prepare_truth_oritab.log" \
    python3 "$SCRIPT_DIR/prepare_truth_oritab.py" --input "$TRUTH_ORIS" \
    --exact-output "$EXACT_ORIS" --output "$PERTURBED_ORIS" \
    --kv "$KV" --cs "$CS" --fraca "$FRACA"; then
    echo "ERROR: cannot create the perturbed orientation table" >&2
    exit 2
fi

CLEAN_STACK="$RESULT_ROOT/generated/clean_particles.mrcs"
CLEAN_ORIS="$RESULT_ROOT/generated/clean_truth_oris.txt"
run_success "$RESULT_ROOT/generated" "$RESULT_ROOT/generated/simulate_clean.log" simple_exec \
    prg=simulate_particles vol1="$TRUTH_VOLUME" oritab="$EXACT_ORIS" \
    outstk="$CLEAN_STACK" outfile="$CLEAN_ORIS" smpd="$SMPD" mskdiam="$MSKDIAM" \
    nptcls="$NPTCLS" snr=10 ctf=yes pgrp="$PGRP" nthr="$NTHR" \
    bfac=0 sherr=0 dferr=0 astigerr=0 bfacerr=0 mkdir=no || true
if [[ ! -f $CLEAN_STACK ]]; then
    echo "ERROR: clean independent particle generation failed" >&2
    exit 2
fi

# Import one observation stack and create matched exact and perturbed base projects.
create_base() {
    local observation=$1
    local start=$2
    local base="$RESULT_ROOT/bases/${observation}_${start}"
    local stack=$NOISY_STACK
    local oris=$EXACT_ORIS
    [[ $observation == clean ]] && stack=$CLEAN_STACK
    [[ $start == perturbed ]] && oris=$PERTURBED_ORIS
    mkdir -p "$base"
    if ! run_success "$base" "$base/new_project.log" simple_exec \
        prg=new_project projname=diagnostic dir=./ qsys_name=local; then
        return
    fi
    run_success "$base" "$base/import_particles.log" simple_exec \
        prg=import_particles projfile="$base/diagnostic.simple" stk="$stack" \
        oritab="$oris" smpd="$SMPD" ctf=yes kv="$KV" cs="$CS" fraca="$FRACA" mkdir=no || return
    if ! export_poses "$base/diagnostic.simple" "$base/imported_poses.txt"; then
        record_failure "base ${observation}_${start} could not export imported particles"
        return
    fi
    local imported_count even_count odd_count
    imported_count=$(awk 'NF >= 8 && $1 ~ /^[0-9]+$/ {n++} END {print n+0}' "$base/imported_poses.txt")
    even_count=$(awk 'NF >= 8 && $1 ~ /^[0-9]+$/ && $8 == 0 {n++} END {print n+0}' "$base/imported_poses.txt")
    odd_count=$(awk 'NF >= 8 && $1 ~ /^[0-9]+$/ && $8 == 1 {n++} END {print n+0}' "$base/imported_poses.txt")
    if (( imported_count != NPTCLS || even_count == 0 || odd_count == 0 )); then
        record_failure "base ${observation}_${start} has invalid particle or half-set counts: total=$imported_count even=$even_count odd=$odd_count"
        return
    fi
    printf 'BASE_IMPORT_PASS total=%d even=%d odd=%d\n' \
        "$imported_count" "$even_count" "$odd_count" | tee "$base/import_validation.txt"
}

for observation in clean noisy; do
    for start in exact perturbed; do
        create_base "$observation" "$start"
    done
done

# Calculate half-map and truth-map FSC curves for one completed arm.
calculate_fsc() {
    local name=$1
    local vol1=$2
    local vol2=$3
    local output=$4
    local workdir="$5"
    mkdir -p "$workdir"
    if ! run_success "$workdir" "$workdir/fsc.log" simple_exec prg=fsc \
        vol1="$vol1" vol2="$vol2" smpd="$SMPD" mskdiam="$MSKDIAM" \
        automsk=no nthr="$NTHR" mkdir=no; then
        record_failure "$name FSC calculation failed"
        return 1
    fi
    if [[ ! -f $workdir/fsc_state01.bin ]]; then
        record_failure "$name FSC did not produce fsc_state01.bin"
        return 1
    fi
    run_success "$workdir" "$output" simple_exec prg=print_fsc \
        fsc="$workdir/fsc_state01.bin" box="$BOX" smpd="$SMPD" mkdir=no
}

# Run one off/on reconstruction arm and collect pose, map, and FSC evidence.
run_arm() {
    local case_name=$1
    local observation=$2
    local start=$3
    local band=$4
    local enabled=$5
    local base="$RESULT_ROOT/bases/${observation}_${start}"
    local arm="$RESULT_ROOT/arms/${case_name}_${enabled}"
    if [[ ! -f $base/diagnostic.simple ]]; then
        record_failure "arm ${case_name}_${enabled} has no imported base project"
        return
    fi
    mkdir -p "$arm"
    cp -a "$base/." "$arm/" || {
        record_failure "arm ${case_name}_${enabled} could not copy its base project"
        return
    }
    local project="$arm/diagnostic.simple"
    local polish_args=(pcg_pose_polish=no)
    [[ $enabled == on ]] && polish_args=(pcg_pose_polish=yes)
    local lp_args=()
    [[ $band == fsc05 ]] && lp_args=(lp="$LP_FSC05")
    [[ $band == fsc0143 ]] && lp_args=(lp="$LP_FSC0143")
    mkdir -p "$arm/evidence"
    export_poses "$project" "$arm/evidence/poses_before.txt" || \
        record_failure "arm ${case_name}_${enabled} could not export initial poses"
    if ! run_success "$arm" "$arm/reconstruct3D.log" simple_exec \
        prg=reconstruct3D projfile="$project" rec_backend=pcg "${polish_args[@]}" \
        "${lp_args[@]}" combine_eo=no objfun=cc ml_reg=no envfsc=no projrec=no \
        pgrp="$PGRP" mskdiam="$MSKDIAM" nthr="$NTHR" nparts=1 \
        maxits=2 rtol=0 postprocess=no mkdir=no; then
        return
    fi
    export_poses "$project" "$arm/evidence/poses_after.txt" || \
        record_failure "arm ${case_name}_${enabled} could not export final poses"
    if [[ ! -f $arm/fsc_state01.bin || ! -f $arm/recvol_state01.mrc || \
          ! -f $arm/recvol_state01_even.mrc || ! -f $arm/recvol_state01_odd.mrc ]]; then
        record_failure "arm ${case_name}_${enabled} is missing reconstruction evidence"
        return
    fi
    run_success "$arm" "$arm/evidence/half_fsc.txt" simple_exec prg=print_fsc \
        fsc="$arm/fsc_state01.bin" box="$BOX" smpd="$SMPD" mkdir=no || true
    calculate_fsc "${case_name}_${enabled} truth average" "$TRUTH_VOLUME" \
        "$arm/recvol_state01.mrc" "$arm/evidence/truth_fsc_avg.txt" \
        "$arm/evidence/truth_avg_calc" || true
    calculate_fsc "${case_name}_${enabled} truth even" "$TRUTH_VOLUME" \
        "$arm/recvol_state01_even.mrc" "$arm/evidence/truth_fsc_even.txt" \
        "$arm/evidence/truth_even_calc" || true
    calculate_fsc "${case_name}_${enabled} truth odd" "$TRUTH_VOLUME" \
        "$arm/recvol_state01_odd.mrc" "$arm/evidence/truth_fsc_odd.txt" \
        "$arm/evidence/truth_odd_calc" || true
}

cat >"$RESULT_ROOT/case_matrix.tsv" <<EOF
clean_exact_full clean exact full
clean_perturbed_full clean perturbed full
noisy_exact_full noisy exact full
noisy_perturbed_full noisy perturbed full
noisy_exact_fsc05 noisy exact fsc05
noisy_perturbed_fsc05 noisy perturbed fsc05
noisy_exact_fsc0143 noisy exact fsc0143
noisy_perturbed_fsc0143 noisy perturbed fsc0143
EOF

while read -r case_name observation start band; do
    run_arm "$case_name" "$observation" "$start" "$band" off
    run_arm "$case_name" "$observation" "$start" "$band" on
done <"$RESULT_ROOT/case_matrix.tsv"

python3 "$SCRIPT_DIR/analyze_truth_matrix.py" --root "$RESULT_ROOT" \
    --truth-oris "$TRUTH_ORIS" --pgrp "$PGRP" 2>&1 | tee "$RESULT_ROOT/analysis/analyzer.log"
ANALYSIS_STATUS=${PIPESTATUS[0]}
if (( ANALYSIS_STATUS != 0 )); then
    record_failure "truth-matrix analyzer returned status $ANALYSIS_STATUS"
fi
if (( FAILURE_COUNT == 0 )); then
    echo PASS >"$RESULT_ROOT/STATUS.txt"
else
    echo FAIL >"$RESULT_ROOT/STATUS.txt"
fi
echo "Diagnostic directory: $RESULT_ROOT"
echo "Download this one directory. Review analysis/truth_matrix.md first."
if (( FAILURE_COUNT == 0 )); then
    exit 0
fi
echo "Completed all scheduled arms with $FAILURE_COUNT recorded failure(s)." >&2
exit 1
