#!/usr/bin/env bash
set -uo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_oracle_validation.sh [--output-parent "$HOME/Projects"]

Run the complete continuous 3-D pose-refinement validation matrix after the
authorized Oracle build. Every test runs independently from the output parent.
The runner retains all statuses, source snapshots, generated comparisons, and
analysis in one timestamped package.
EOF
}

OUTPUT_PARENT="${HOME}/Projects"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-parent) OUTPUT_PARENT=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE_ROOT=$(readlink -f "${SCRIPT_DIR}/../../..")
mkdir -p "$OUTPUT_PARENT"
OUTPUT_PARENT=$(readlink -f "$OUTPUT_PARENT")
RESULT_ROOT="${OUTPUT_PARENT}/continuous_3D_pose_validation_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULT_ROOT"/{analysis,capture/artifacts,capture/logs,generated_matrices,input,logs,source,tests}
exec > >(tee -a "$RESULT_ROOT/run.log") 2>&1

echo "Continuous 3-D pose validation directory: $RESULT_ROOT"
echo "Runtime working directory: $OUTPUT_PARENT"

FAILURE_COUNT=0
for required_tool in python3 sha256sum find sort awk sed grep xargs comm cp readlink; do
    if ! command -v "$required_tool" >/dev/null 2>&1; then
        echo "ERROR: required tool is not available: $required_tool" >&2
        FAILURE_COUNT=$((FAILURE_COUNT+1))
    fi
done

resolve_executable() {
    local requested=$1
    local resolved

    resolved=$(command -v "$requested" 2>/dev/null || true)
    if [[ -z $resolved || ! -x $resolved ]]; then
        echo ""
        return 1
    fi
    readlink -f "$resolved"
}

POSE_EXEC=$(resolve_executable simple_test_pose_cont_refinement || true)
CART_EXEC=$(resolve_executable simple_test_cartesian_fourier || true)
PCG_EXEC=$(resolve_executable simple_test_continuous_3D_pcg_reconstruction || true)
SIMPLE_TEST_EXEC=$(resolve_executable simple_test_exec || true)
for executable_name in POSE_EXEC CART_EXEC PCG_EXEC SIMPLE_TEST_EXEC; do
    executable_value=${!executable_name}
    if [[ -z $executable_value ]]; then
        echo "ERROR: required executable is unavailable: $executable_name" >&2
        FAILURE_COUNT=$((FAILURE_COUNT+1))
    fi
done

{
    echo "date=$(date --iso-8601=seconds)"
    echo "host=$(hostname)"
    echo "runtime_working_directory=$OUTPUT_PARENT"
    echo "source_root=$SOURCE_ROOT"
    echo "simple_path=${SIMPLE_PATH:-unset}"
    echo "path=$PATH"
    echo "omp_num_threads=${OMP_NUM_THREADS:-unset}"
    echo "gcc_version=$(gcc --version 2>/dev/null | sed -n '1p')"
    echo "gfortran_version=$(gfortran --version 2>/dev/null | sed -n '1p')"
    echo "python_version=$(python3 --version 2>&1 || true)"
    echo "pose_executable=$POSE_EXEC"
    echo "cartesian_executable=$CART_EXEC"
    echo "pcg_executable=$PCG_EXEC"
    echo "simple_test_exec=$SIMPLE_TEST_EXEC"
    echo "module_list_begin"
    module list 2>&1 || true
    echo "module_list_end"
    uname -a
} >"$RESULT_ROOT/input/environment.txt"

snapshot_sources() {
    local list_path="$RESULT_ROOT/input/source_files.txt"
    local relative_path

    {
        echo "doc/implementation_notes/continuous_3D_pose_end_polishing.md"
        echo "src/main/volume/simple_cartesian_pose_refiner.f90"
        echo "src/main/interp/simple_cartesian_fourier.f90"
        echo "src/main/interp/simple_gridding.f90"
        echo "src/main/volume/simple_reconstructor_pcg.f90"
        echo "production/simple_test_exec.f90"
        find "$SOURCE_ROOT/production/tests" -maxdepth 1 -type f \
            \( -name 'simple_pose_cont_refinement*.f90' \
            -o -name 'simple_test_pose_cont_refinement.f90' \
            -o -name 'simple_cartesian_fourier*.f90' \
            -o -name 'simple_test_cartesian_fourier.f90' \
            -o -name 'simple_continuous_3D_pcg_reconstruction*.f90' \
            -o -name 'simple_test_continuous_3D_pcg_reconstruction.f90' \) \
            -printf '%P\n' | sed 's|^|production/tests/|'
        find "$SCRIPT_DIR" -maxdepth 1 -type f -printf '%P\n' \
            | sed 's|^|production/tests/pose_cont_validation/|'
        find "$SOURCE_ROOT/production/tests/pose_cont_capture" -maxdepth 1 -type f \
            -printf '%P\n' | sed 's|^|production/tests/pose_cont_capture/|'
    } | sort -u >"$list_path"

    while IFS= read -r relative_path; do
        if [[ ! -f "$SOURCE_ROOT/$relative_path" ]]; then
            echo "ERROR: required source snapshot is missing: $relative_path" >&2
            FAILURE_COUNT=$((FAILURE_COUNT+1))
            continue
        fi
        mkdir -p "$RESULT_ROOT/source/$(dirname "$relative_path")"
        cp "$SOURCE_ROOT/$relative_path" "$RESULT_ROOT/source/$relative_path"
    done <"$list_path"

    (cd "$RESULT_ROOT/source" && find . -type f -print0 | sort -z \
        | xargs -0 sha256sum) >"$RESULT_ROOT/input/source.sha256"
}

snapshot_sources
if [[ -f "$SOURCE_ROOT/build/build_debug.log" ]]; then
    cp "$SOURCE_ROOT/build/build_debug.log" "$RESULT_ROOT/logs/build_debug.log"
fi
if [[ -f "$SOURCE_ROOT/build_debug_incremental_failed.log" ]]; then
    cp "$SOURCE_ROOT/build_debug_incremental_failed.log" \
        "$RESULT_ROOT/logs/build_debug_incremental_failed.log"
fi

printf 'test_id\texecutable\targuments\tworking_directory\texit_status\n' \
    >"$RESULT_ROOT/statuses.tsv"
printf 'test_id\texecutable\tsha256\n' >"$RESULT_ROOT/executables.tsv"

record_executable() {
    local test_id=$1
    local executable=$2
    local digest
    digest=$(sha256sum "$executable" | awk '{print $1}')
    printf '%s\t%s\t%s\n' "$test_id" "$executable" "$digest" \
        >>"$RESULT_ROOT/executables.tsv"
}

run_test() {
    local test_id=$1
    local executable=$2
    shift 2
    local test_root="$RESULT_ROOT/tests/$test_id"
    local status
    local arg

    mkdir -p "$test_root/evidence"
    record_executable "$test_id" "$executable"
    {
        printf 'executable=%q\n' "$executable"
        printf 'working_directory=%q\n' "$OUTPUT_PARENT"
        printf 'arguments='
        for arg in "$@"; do printf '%q ' "$arg"; done
        printf '\n'
    } >"$test_root/command.txt"

    echo ">>> TEST [$test_id] START"
    (
        cd "$OUTPUT_PARENT" || exit 125
        "$executable" "$@"
    ) 2>&1 | tee "$test_root/run.log"
    status=${PIPESTATUS[0]}
    printf '%s\n' "$status" >"$test_root/exit_status.txt"
    printf '%s\t%s\t' "$test_id" "$executable" >>"$RESULT_ROOT/statuses.tsv"
    printf '%q ' "$@" >>"$RESULT_ROOT/statuses.tsv"
    printf '\t%s\t%s\n' "$OUTPUT_PARENT" "$status" >>"$RESULT_ROOT/statuses.tsv"
    echo ">>> TEST [$test_id] STATUS $status"
}

run_capture_test() {
    local test_id=$1
    local case_name=$2
    local log_name=$3
    local status

    record_executable "$test_id" "$POSE_EXEC"
    {
        printf 'executable=%q\n' "$POSE_EXEC"
        printf 'working_directory=%q\n' "$OUTPUT_PARENT"
        printf 'environment=CONTINUOUS_3D_POSE_CAPTURE_DIR=%q\n' "$RESULT_ROOT/capture/artifacts"
        printf 'arguments=%q\n' "case=$case_name"
    } >"$RESULT_ROOT/capture/${test_id}_command.txt"
    echo ">>> TEST [$test_id] START"
    (
        cd "$OUTPUT_PARENT" || exit 125
        CONTINUOUS_3D_POSE_CAPTURE_DIR="$RESULT_ROOT/capture/artifacts" \
            "$POSE_EXEC" "case=$case_name"
    ) 2>&1 | tee "$RESULT_ROOT/capture/logs/$log_name"
    status=${PIPESTATUS[0]}
    printf '%s\n' "$status" >"$RESULT_ROOT/capture/${test_id}_exit_status.txt"
    printf '%s\t%s\t%s\t%s\t%s\n' "$test_id" "$POSE_EXEC" "case=$case_name" \
        "$OUTPUT_PARENT" "$status" >>"$RESULT_ROOT/statuses.tsv"
    echo ">>> TEST [$test_id] STATUS $status"
}

find "$OUTPUT_PARENT" -maxdepth 1 -type d -name 'continuous_3D_matrix_volumes_*' \
    -printf '%f\n' | sort >"$RESULT_ROOT/input/matrix_roots_before.txt"

if (( FAILURE_COUNT == 0 )); then
    run_test configuration_matrix "$POSE_EXEC" case=configuration_matrix \
        "evidence_dir=$RESULT_ROOT/tests/configuration_matrix/evidence"
    run_test pose_mother "$POSE_EXEC"
    run_test fixed_reference "$POSE_EXEC" case=fixed_reference
    run_test forward_path "$POSE_EXEC" case=forward_path
    run_test matched_window "$POSE_EXEC" case=matched_window
    run_test reference_bias "$POSE_EXEC" case=reference_bias
    run_test operator_contract "$POSE_EXEC" case=operator_contract
    run_capture_test pose_capture_range pose_capture_range pose_capture.log
    run_capture_test pose_capture_mechanism pose_capture_mechanism pose_mechanism.log
    run_test tolerance_calibration "$POSE_EXEC" case=tolerance_calibration \
        "evidence_dir=$RESULT_ROOT/tests/tolerance_calibration/evidence"
    run_test objective_normals "$POSE_EXEC" case=objective_normals \
        "evidence_dir=$RESULT_ROOT/tests/objective_normals/evidence"
    run_test lm_transactions "$POSE_EXEC" case=lm_transactions \
        "evidence_dir=$RESULT_ROOT/tests/lm_transactions/evidence"
    run_test ctf_sigma "$POSE_EXEC" case=ctf_sigma \
        "evidence_dir=$RESULT_ROOT/tests/ctf_sigma/evidence"
    run_test forward_hierarchy "$POSE_EXEC" case=forward_hierarchy \
        "evidence_dir=$RESULT_ROOT/tests/forward_hierarchy/evidence"
    run_test cartesian_mother "$CART_EXEC"
    run_test pcg_mother "$PCG_EXEC"
    run_test pcg_recon "$SIMPLE_TEST_EXEC" test=pcg_recon
    run_test pcg_priors "$SIMPLE_TEST_EXEC" test=pcg_priors
fi

find "$OUTPUT_PARENT" -maxdepth 1 -type d -name 'continuous_3D_matrix_volumes_*' \
    -printf '%f\n' | sort >"$RESULT_ROOT/input/matrix_roots_after.txt"
comm -13 "$RESULT_ROOT/input/matrix_roots_before.txt" \
    "$RESULT_ROOT/input/matrix_roots_after.txt" >"$RESULT_ROOT/input/generated_matrix_roots.txt"
while IFS= read -r matrix_root; do
    [[ -z $matrix_root ]] && continue
    cp -a "$OUTPUT_PARENT/$matrix_root" "$RESULT_ROOT/generated_matrices/"
done <"$RESULT_ROOT/input/generated_matrix_roots.txt"

ANALYSIS_STATUS=2
if (( FAILURE_COUNT == 0 )); then
    python3 "$SCRIPT_DIR/analyze_oracle_validation.py" --result-root "$RESULT_ROOT" \
        2>&1 | tee "$RESULT_ROOT/analysis/analysis.log"
    ANALYSIS_STATUS=${PIPESTATUS[0]}
fi
printf '%s\n' "$ANALYSIS_STATUS" >"$RESULT_ROOT/analysis/exit_status.txt"

if grep -R -Fq 'pcg_pose_polish' "$RESULT_ROOT"/tests/*/command.txt \
    "$RESULT_ROOT"/capture/*_command.txt 2>/dev/null; then
    echo "ERROR: command manifest invokes removed pcg_pose_polish option" >&2
    FAILURE_COUNT=$((FAILURE_COUNT+1))
fi

if (( FAILURE_COUNT == 0 && ANALYSIS_STATUS == 0 )); then
    echo PASS >"$RESULT_ROOT/STATUS.txt"
else
    echo FAIL >"$RESULT_ROOT/STATUS.txt"
fi

(cd "$RESULT_ROOT" && find . -type f ! -name MANIFEST.sha256 ! -name run.log -print0 | sort -z \
    | xargs -0 sha256sum) >"$RESULT_ROOT/MANIFEST.sha256"

echo "Package status: $(cat "$RESULT_ROOT/STATUS.txt")"
echo "Package path: $RESULT_ROOT"
exit $(( FAILURE_COUNT == 0 && ANALYSIS_STATUS == 0 ? 0 : 1 ))
