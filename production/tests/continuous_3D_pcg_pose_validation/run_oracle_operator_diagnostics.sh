#!/usr/bin/env bash
set -uo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_oracle_operator_diagnostics.sh [--output-parent "$HOME/Projects"] \
    [--executable simple_test_continuous_3D_pcg_refinement]

The runner collects the five operator/reference diagnostics and the unchanged
mother suite in one timestamped directory. A failed case is recorded without
hiding later independent evidence.
EOF
}

OUTPUT_PARENT="$PWD"
EXECUTABLE="simple_test_continuous_3D_pcg_refinement"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-parent) OUTPUT_PARENT=$2; shift 2 ;;
        --executable) EXECUTABLE=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if ! command -v "$EXECUTABLE" >/dev/null 2>&1 && [[ ! -x $EXECUTABLE ]]; then
    echo "ERROR: test executable is not available: $EXECUTABLE" >&2
    exit 2
fi
mkdir -p "$OUTPUT_PARENT"
OUTPUT_PARENT=$(readlink -f "$OUTPUT_PARENT")
RESULT_ROOT="$OUTPUT_PARENT/continuous_3D_operator_diagnostics_$(date +%Y%m%d_%H%M%S)"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE_ROOT=$(readlink -f "$SCRIPT_DIR/../../..")
mkdir -p "$RESULT_ROOT"/{logs,input}
cp "$SCRIPT_DIR/run_oracle_operator_diagnostics.sh" "$RESULT_ROOT/input/"
cp "$SCRIPT_DIR/README.md" "$RESULT_ROOT/input/"
cp "$SOURCE_ROOT/doc/implementation_notes/continuous_3D_pose_end_polishing.md" \
    "$RESULT_ROOT/input/"
cp "$SOURCE_ROOT/doc/implementation_notes/continuous_3D_refinement_on_pcg_operator.md" \
    "$RESULT_ROOT/input/"
cp "$SOURCE_ROOT/doc/policies/reconstruct3D_pcg_policy.md" "$RESULT_ROOT/input/"
exec > >(tee -a "$RESULT_ROOT/diagnostics.log") 2>&1
echo "Diagnostic directory: $RESULT_ROOT"

FAILURE_COUNT=0
run_case() {
    local label=$1
    local marker=$2
    local logfile="$RESULT_ROOT/logs/${label}.log"
    echo ">>> DIAGNOSTIC [$label] START"
    if [[ $label == mother_suite ]]; then
        "$EXECUTABLE" 2>&1 | tee "$logfile"
    else
        "$EXECUTABLE" "case=$label" 2>&1 | tee "$logfile"
    fi
    local status=${PIPESTATUS[0]}
    if (( status != 0 )); then
        FAILURE_COUNT=$((FAILURE_COUNT+1))
        echo "FAIL: $label returned status $status" >&2
    elif ! grep -Fq "$marker" "$logfile"; then
        FAILURE_COUNT=$((FAILURE_COUNT+1))
        echo "FAIL: $label did not emit its completion marker" >&2
    else
        echo ">>> DIAGNOSTIC [$label] COMPLETE"
    fi
}

run_case fixed_reference 'CONTINUOUS_3D_FIXED_REFERENCE: EVIDENCE COMPLETE'
run_case forward_path 'CONTINUOUS_3D_FORWARD_PATH: EVIDENCE COMPLETE'
run_case matched_window 'CONTINUOUS_3D_MATCHED_WINDOW: EVIDENCE COMPLETE'
run_case reference_bias 'CONTINUOUS_3D_REFERENCE_BIAS: EVIDENCE COMPLETE'
run_case operator_contract 'CONTINUOUS_3D_OPERATOR_CONTRACT: EVIDENCE COMPLETE'
run_case mother_suite 'Continuous 3D PCG refinement test suite: PASS'

{
    echo '# Continuous 3D operator diagnostics'
    echo
    echo "Recorded failures: $FAILURE_COUNT"
    echo
    echo '## Diagnosis markers'
    grep -hE 'DIAGNOSIS:|components (gauge/leave-in|clip/taper)|reconstructed-reference deapod|test suite: (PASS|FAIL)' \
        "$RESULT_ROOT"/logs/*.log || true
    echo
    echo 'Review the full reference-bias arm metrics in `logs/reference_bias.log`.'
    echo 'Review the four-model and box-margin matrix in `logs/operator_contract.log`.'
} >"$RESULT_ROOT/summary.md"

cat >"$RESULT_ROOT/input/source_files.txt" <<'EOF'
production/tests/simple_continuous_3D_pcg_refinement_fixed_reference_test.f90
production/tests/simple_continuous_3D_pcg_refinement_forward_path_test.f90
production/tests/simple_continuous_3D_pcg_refinement_matched_window_test.f90
production/tests/simple_continuous_3D_pcg_refinement_reference_bias_test.f90
production/tests/simple_continuous_3D_pcg_refinement_operator_contract_support.f90
production/tests/simple_continuous_3D_pcg_refinement_operator_contract_test.f90
production/tests/simple_continuous_3D_pcg_refinement_halfset_support.f90
production/tests/simple_test_continuous_3D_pcg_refinement.f90
src/main/strategies/search/simple_pcg_pose_polisher.f90
src/main/volume/simple_reconstructor_pcg.f90
src/main/image/simple_projector.f90
src/main/commanders/test/simple_commanders_test_highlevel.f90
doc/implementation_notes/continuous_3D_pose_end_polishing.md
doc/implementation_notes/continuous_3D_refinement_on_pcg_operator.md
doc/policies/reconstruct3D_pcg_policy.md
EOF
while IFS= read -r relative; do
    if [[ -f $SOURCE_ROOT/$relative ]]; then
        sha256sum "$SOURCE_ROOT/$relative"
    else
        printf 'MISSING  %s\n' "$SOURCE_ROOT/$relative"
    fi
done <"$RESULT_ROOT/input/source_files.txt" >"$RESULT_ROOT/input/source_manifest.sha256"

if (( FAILURE_COUNT == 0 )); then
    echo PASS >"$RESULT_ROOT/STATUS.txt"
    echo "CONTINUOUS_3D_OPERATOR_DIAGNOSTICS: COMPLETE"
else
    echo FAIL >"$RESULT_ROOT/STATUS.txt"
    echo "CONTINUOUS_3D_OPERATOR_DIAGNOSTICS: FAIL ($FAILURE_COUNT recorded failures)" >&2
fi
(cd "$RESULT_ROOT" && find . -type f ! -name MANIFEST.sha256 ! -name diagnostics.log -print0 | \
    sort -z | xargs -0 sha256sum) >"$RESULT_ROOT/MANIFEST.sha256"
echo "Download this directory: $RESULT_ROOT"
(( FAILURE_COUNT == 0 ))
