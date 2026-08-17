#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_oracle_validation.sh \
    --projfile /absolute/path/to/checkpoint/project.simple \
    --box 144 --smpd 1.3 --mskdiam 120 --pgrp c1 \
    [--objfun cc] [--nthr 8] [--distributed-parts 2] \
    [--output-parent "$PWD"]

The source checkpoint is read-only. Five complete arm directories and all
logs, metadata exports, FSC tables, analysis files, and checksums are written
under one timestamped output directory. After preflight, a failed test or arm
is recorded and the remaining independent work continues.
EOF
}

PROJFILE=""
BOX=""
SMPD=""
MSKDIAM=""
PGRP=""
OBJFUN=cc
NTHR=8
DISTRIBUTED_PARTS=2
OUTPUT_PARENT="$PWD"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --projfile) PROJFILE=$2; shift 2 ;;
        --box) BOX=$2; shift 2 ;;
        --smpd) SMPD=$2; shift 2 ;;
        --mskdiam) MSKDIAM=$2; shift 2 ;;
        --pgrp) PGRP=$2; shift 2 ;;
        --objfun) OBJFUN=$2; shift 2 ;;
        --nthr) NTHR=$2; shift 2 ;;
        --distributed-parts) DISTRIBUTED_PARTS=$2; shift 2 ;;
        --output-parent) OUTPUT_PARENT=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

for value_name in PROJFILE BOX SMPD MSKDIAM PGRP; do
    if [[ -z ${!value_name} ]]; then
        echo "ERROR: --${value_name,,} is required" >&2
        usage >&2
        exit 2
    fi
done
if [[ ! -f $PROJFILE ]]; then
    echo "ERROR: project file does not exist: $PROJFILE" >&2
    exit 2
fi
if (( DISTRIBUTED_PARTS < 2 )); then
    echo "ERROR: --distributed-parts must be at least 2" >&2
    exit 2
fi
if [[ $OBJFUN != cc && $OBJFUN != euclid ]]; then
    echo "ERROR: --objfun must be cc or euclid" >&2
    exit 2
fi
# Require the Python syntax level used by the evidence analyzers.
python_supported() {
    command -v python3 >/dev/null 2>&1 && \
        python3 -c 'import sys; raise SystemExit(sys.version_info < (3, 7))'
}

if ! python_supported && type module >/dev/null 2>&1; then
    echo "Loading miniconda/25.9.1 for the validation analyzer"
    module load miniconda/25.9.1
fi
for executable in simple_exec simple_private_exec simple_test_continuous_3D_pcg_refinement python3; do
    if ! command -v "$executable" >/dev/null 2>&1; then
        echo "ERROR: required executable is not on PATH: $executable" >&2
        exit 2
    fi
done
if ! python_supported; then
    echo "ERROR: Python 3.7 or newer is required; load miniconda/25.9.1 on Oracle Linux" >&2
    exit 2
fi

PROJFILE=$(readlink -f "$PROJFILE")
CHECKPOINT_DIR=$(dirname "$PROJFILE")
PROJECT_NAME=$(basename "$PROJFILE")
mkdir -p "$OUTPUT_PARENT"
OUTPUT_PARENT=$(readlink -f "$OUTPUT_PARENT")
RESULT_ROOT="$OUTPUT_PARENT/continuous_3D_pose_validation_$(date +%Y%m%d_%H%M%S)"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
case "$RESULT_ROOT/" in
    "$CHECKPOINT_DIR/"*)
        echo "ERROR: output directory cannot be inside the source checkpoint" >&2
        exit 2
        ;;
esac
mkdir -p "$RESULT_ROOT"/{unit,policy,input,arms,analysis}
mkdir -p "$RESULT_ROOT/input/validation_source"
cp "$SCRIPT_DIR/run_oracle_validation.sh" "$RESULT_ROOT/input/validation_source/"
cp "$SCRIPT_DIR/analyze_pose_ab.py" "$RESULT_ROOT/input/validation_source/"
cp "$SCRIPT_DIR/README.md" "$RESULT_ROOT/input/validation_source/"
exec > >(tee -a "$RESULT_ROOT/validation.log") 2>&1
echo "Validation directory: $RESULT_ROOT"

record_failed_run() {
    local status=$?
    if (( status == 0 )); then
        return
    fi
    if [[ -f $RESULT_ROOT/STATUS.txt ]]; then
        echo "Validation completed with failures. Download: $RESULT_ROOT" >&2
        return
    fi
    set +e
    echo FAIL >"$RESULT_ROOT/STATUS.txt"
    (cd "$RESULT_ROOT" && find . -type f ! -name MANIFEST.sha256 ! -name validation.log -print0 | \
        sort -z | xargs -0 sha256sum) >"$RESULT_ROOT/MANIFEST.sha256"
    echo "Validation stopped unexpectedly with status $status. Download: $RESULT_ROOT" >&2
}
trap record_failed_run EXIT
FAILURE_COUNT=0
FAILURE_LOG="$RESULT_ROOT/failures.txt"
set +e

record_failure() {
    FAILURE_COUNT=$((FAILURE_COUNT+1))
    printf 'FAIL [%02d] %s\n' "$FAILURE_COUNT" "$*" | tee -a "$FAILURE_LOG" >&2
}

MASK_RADIUS_PIXELS=$(python3 -c 'import sys; print(float(sys.argv[1])/(2.0*float(sys.argv[2])))' "$MSKDIAM" "$SMPD")
ANGULAR_BOUND=$(python3 -c 'import sys; print(8.0/max(float(sys.argv[1]),1.0))' "$MASK_RADIUS_PIXELS")
cat >"$RESULT_ROOT/run_config.env" <<EOF
PROJFILE=$PROJFILE
BOX=$BOX
SMPD=$SMPD
MSKDIAM=$MSKDIAM
PGRP=$PGRP
OBJFUN=$OBJFUN
NTHR=$NTHR
DISTRIBUTED_PARTS=$DISTRIBUTED_PARTS
MASK_RADIUS_PIXELS=$MASK_RADIUS_PIXELS
CUMULATIVE_ANGULAR_BOUND_RAD=$ANGULAR_BOUND
EOF
{
    date --iso-8601=seconds
    uname -a
    command -v simple_exec
    command -v simple_private_exec
    command -v simple_test_continuous_3D_pcg_refinement
    python3 --version
} >"$RESULT_ROOT/environment.txt"
git -C "$SCRIPT_DIR/../../.." rev-parse HEAD >"$RESULT_ROOT/input/source_commit.txt" 2>/dev/null || true
(cd "$CHECKPOINT_DIR" && find . -type f -print0 | sort -z | xargs -0 sha256sum) \
    >"$RESULT_ROOT/input/checkpoint_manifest.sha256"

# Run one required command and record a failure without stopping independent arms.
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
    return 0
}

# Require a typed policy rejection and reject an earlier UI-input failure.
expect_failure() {
    local name=$1
    local expected=$2
    shift 2
    local workdir="$RESULT_ROOT/policy/$name"
    local logfile="$workdir/test.log"
    local command=("$@" "projfile=$workdir/$PROJECT_NAME" "pgrp=$PGRP" \
        "mskdiam=$MSKDIAM" "nthr=$NTHR")
    local failed=0
    mkdir -p "$workdir"
    if ! cp "$PROJFILE" "$workdir/$PROJECT_NAME"; then
        record_failure "policy $name could not copy its project"
        return
    fi
    sha256sum "$workdir/$PROJECT_NAME" >"$workdir/project_before.sha256"
    printf '%q ' "${command[@]}" >"$workdir/command.txt"
    printf '\n' >>"$workdir/command.txt"
    (cd "$workdir" && "${command[@]}") >"$logfile" 2>&1
    local status=$?
    if grep -Eqi 'Missing key on command line|not enough input variables defined' "$logfile"; then
        record_failure "policy $name did not reach parameter validation"
        failed=1
    fi
    if (( status == 0 )); then
        record_failure "policy $name succeeded but failure was required"
        failed=1
    fi
    if ! grep -Eqi "$expected" "$logfile"; then
        record_failure "policy $name did not report the expected reason"
        failed=1
    fi
    sha256sum "$workdir/$PROJECT_NAME" >"$workdir/project_after.sha256"
    if [[ $(cut -d' ' -f1 "$workdir/project_before.sha256") != \
          $(cut -d' ' -f1 "$workdir/project_after.sha256") ]]; then
        record_failure "policy $name changed its project before failing"
        failed=1
    fi
    if (( failed == 0 )); then
        echo "POLICY [$name]: PASS" | tee "$workdir/result.txt"
    else
        echo "POLICY [$name]: FAIL" | tee "$workdir/result.txt"
    fi
}

# Export pose and immutable metadata for matched A/B comparison.
export_project_evidence() {
    local project=$1
    local evidence_dir=$2
    local suffix=$3
    local status=0
    simple_private_exec prg=print_project_vals projfile="$project" oritype=ptcl3D \
        keys=e1,e2,e3,x,y,eo >"$evidence_dir/poses_${suffix}.txt" 2>&1 || status=1
    simple_private_exec prg=print_project_vals projfile="$project" oritype=ptcl3D \
        keys=dfx,dfy,angast,eo >"$evidence_dir/invariants_${suffix}.txt" 2>&1 || status=1
    return "$status"
}

# Run one reconstruction arm from the frozen checkpoint and collect final evidence.
run_arm() {
    local route=$1
    local enabled=$2
    local nparts=$3
    local arm_name="${route}_${enabled}"
    local arm_dir="$RESULT_ROOT/arms/$arm_name"
    mkdir -p "$arm_dir"
    if ! cp -a "$CHECKPOINT_DIR/." "$arm_dir/"; then
        record_failure "arm $arm_name could not copy the frozen checkpoint"
        return
    fi
    local arm_project="$arm_dir/$PROJECT_NAME"
    mkdir -p "$arm_dir/evidence"
    sha256sum "$arm_project" >"$arm_dir/evidence/project_before.sha256"
    if ! export_project_evidence "$arm_project" "$arm_dir/evidence" before; then
        record_failure "arm $arm_name could not export its initial project evidence"
    fi
    local polish_args=()
    if [[ $enabled == on ]]; then
        polish_args=(pcg_pose_polish=yes)
    elif [[ $enabled == off ]]; then
        polish_args=(pcg_pose_polish=no)
    elif [[ $enabled != default ]]; then
        record_failure "arm $arm_name has an unknown mode"
        return
    fi
    if ! run_success "$arm_dir" "$arm_dir/reconstruct3D.log" simple_exec \
        prg=reconstruct3D projfile="$arm_project" rec_backend=pcg "${polish_args[@]}" \
        combine_eo=no objfun="$OBJFUN" ml_reg=no envfsc=no projrec=no \
        pgrp="$PGRP" mskdiam="$MSKDIAM" nthr="$NTHR" nparts="$nparts" \
        maxits=2 rtol=0 mkdir=no; then
        return
    fi
    sha256sum "$arm_project" >"$arm_dir/evidence/project_after.sha256"
    if ! export_project_evidence "$arm_project" "$arm_dir/evidence" after; then
        record_failure "arm $arm_name could not export its final project evidence"
    fi
    if [[ ! -f $arm_dir/fsc_state01.bin ]]; then
        record_failure "arm $arm_name did not produce fsc_state01.bin"
        return
    fi
    run_success "$arm_dir" "$arm_dir/evidence/fsc.txt" simple_exec prg=print_fsc \
        fsc="$arm_dir/fsc_state01.bin" box="$BOX" smpd="$SMPD" mkdir=no || true
}

run_success "$RESULT_ROOT/unit" "$RESULT_ROOT/unit/pose_contract.log" \
    simple_test_continuous_3D_pcg_refinement case=pose_contract || true
run_success "$RESULT_ROOT/unit" "$RESULT_ROOT/unit/mother.log" \
    simple_test_continuous_3D_pcg_refinement || true

expect_failure invalid_backend 'requires rec_backend=pcg' simple_exec prg=reconstruct3D \
    rec_backend=gridding pcg_pose_polish=yes combine_eo=no mkdir=no
expect_failure invalid_combined_halves 'requires independent even and odd' simple_exec prg=reconstruct3D \
    rec_backend=pcg pcg_pose_polish=yes combine_eo=yes mkdir=no
expect_failure invalid_value 'must be yes or no' simple_exec prg=reconstruct3D \
    rec_backend=pcg pcg_pose_polish=maybe combine_eo=no mkdir=no
expect_failure invalid_program 'not allowed|accepted only by reconstruct3D' simple_exec prg=refine3D \
    rec_backend=pcg pcg_pose_polish=yes combine_eo=no mkdir=no

run_arm shared default 1
run_arm shared off 1
run_arm shared on 1
run_arm distributed off "$DISTRIBUTED_PARTS"
run_arm distributed on "$DISTRIBUTED_PARTS"

python3 "$SCRIPT_DIR/analyze_pose_ab.py" --root "$RESULT_ROOT" --pgrp "$PGRP" \
    --angular-bound "$ANGULAR_BOUND" 2>&1 | tee "$RESULT_ROOT/analysis/analyzer.log"
ANALYSIS_STATUS=${PIPESTATUS[0]}
if (( ANALYSIS_STATUS != 0 )); then
    record_failure "aggregate A/B analyzer returned status $ANALYSIS_STATUS"
fi
if (( FAILURE_COUNT == 0 )); then
    echo PASS >"$RESULT_ROOT/STATUS.txt"
else
    echo FAIL >"$RESULT_ROOT/STATUS.txt"
fi
(cd "$RESULT_ROOT" && find . -type f ! -name MANIFEST.sha256 ! -name validation.log -print0 | \
    sort -z | xargs -0 sha256sum) \
    >"$RESULT_ROOT/MANIFEST.sha256"
echo "Validation directory: $RESULT_ROOT"
echo "Download this one directory with rsync. Review analysis/summary.md first."
if (( FAILURE_COUNT == 0 )); then
    exit 0
fi
echo "Completed every scheduled test and arm with $FAILURE_COUNT recorded failure(s)." >&2
exit 1
