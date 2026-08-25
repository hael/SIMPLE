#!/usr/bin/env bash
set -uo pipefail

usage() {
    cat <<'EOF'
Usage:
  run_oracle_pose_capture.sh [--output-parent "$HOME/Projects"] \
    [--executable simple_test_continuous_3D_pcg_refinement]

Run the isolated PCG Cartesian pose-capture sweep. The script stores the raw
log, CSV metrics, MRC images, analysis, source provenance, and checksums in one
timestamped directory. It does not run refine3D or the full mother suite.
EOF
}

OUTPUT_PARENT="${HOME}/Projects"
EXECUTABLE="simple_test_continuous_3D_pcg_refinement"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output-parent) OUTPUT_PARENT=$2; shift 2 ;;
        --executable) EXECUTABLE=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

mkdir -p "$OUTPUT_PARENT"
OUTPUT_PARENT=$(readlink -f "$OUTPUT_PARENT")
RESULT_ROOT="$OUTPUT_PARENT/continuous_3D_pose_capture_$(date +%Y%m%d_%H%M%S)"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE_ROOT=$(readlink -f "$SCRIPT_DIR/../../..")
mkdir -p "$RESULT_ROOT"/{analysis,artifacts,input,logs}
exec > >(tee -a "$RESULT_ROOT/run.log") 2>&1
echo "Pose-capture directory: $RESULT_ROOT"

FAILURE_COUNT=0
if ! command -v "$EXECUTABLE" >/dev/null 2>&1 && [[ ! -x $EXECUTABLE ]]; then
    echo "ERROR: test executable is not available: $EXECUTABLE" >&2
    FAILURE_COUNT=$((FAILURE_COUNT+1))
fi
if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: python3 is required for the result analysis" >&2
    FAILURE_COUNT=$((FAILURE_COUNT+1))
fi
if ! command -v sha256sum >/dev/null 2>&1; then
    echo "ERROR: sha256sum is required for the evidence manifest" >&2
    FAILURE_COUNT=$((FAILURE_COUNT+1))
fi

cp "$SCRIPT_DIR/run_oracle_pose_capture.sh" "$RESULT_ROOT/input/"
cp "$SCRIPT_DIR/analyze_pose_capture.py" "$RESULT_ROOT/input/"
cp "$SCRIPT_DIR/reanalyze_oracle_pose_capture.sh" "$RESULT_ROOT/input/"
cp "$SCRIPT_DIR/README.md" "$RESULT_ROOT/input/"
cp "$SOURCE_ROOT/production/tests/simple_continuous_3D_pcg_refinement_pose_capture_test.f90" \
    "$RESULT_ROOT/input/"
cp "$SOURCE_ROOT/production/tests/simple_continuous_3D_pcg_refinement_pose_mechanism_test.f90" \
    "$RESULT_ROOT/input/"
cp "$SOURCE_ROOT/production/tests/simple_test_continuous_3D_pcg_refinement.f90" \
    "$RESULT_ROOT/input/"
cp "$SOURCE_ROOT/doc/implementation_notes/completed/continuous_3D_pose_end_polishing_validation_evidence.md" \
    "$RESULT_ROOT/input/"
cp "$SOURCE_ROOT/doc/implementation_notes/continuous_3D_refinement_on_pcg_operator.md" \
    "$RESULT_ROOT/input/"

{
    echo "date=$(date --iso-8601=seconds)"
    echo "host=$(hostname)"
    echo "source_root=$SOURCE_ROOT"
    echo "executable=$EXECUTABLE"
    echo "resolved_executable=$(command -v "$EXECUTABLE" 2>/dev/null || true)"
    echo "command=$EXECUTABLE case=pose_capture_range"
    echo "gcc_version=$(gcc --version 2>/dev/null | sed -n '1p')"
    echo "python_version=$(python3 --version 2>&1 || true)"
    echo "git_commit=$(git -C "$SOURCE_ROOT" rev-parse HEAD 2>/dev/null || echo unavailable)"
    echo "git_status_begin"
    git -C "$SOURCE_ROOT" status --short --branch 2>/dev/null || true
    echo "git_status_end"
    uname -a
} >"$RESULT_ROOT/input/environment.txt"

TEST_STATUS=2
MECHANISM_STATUS=2
REGRESSION_STATUS=2
if (( FAILURE_COUNT == 0 )); then
    echo ">>> POSE CAPTURE TEST START"
    export CONTINUOUS_3D_POSE_CAPTURE_DIR="$RESULT_ROOT/artifacts"
    "$EXECUTABLE" case=pose_capture_range 2>&1 | tee "$RESULT_ROOT/logs/pose_capture.log"
    TEST_STATUS=${PIPESTATUS[0]}
    if (( TEST_STATUS != 0 )); then
        echo "FAIL: pose-capture test returned status $TEST_STATUS" >&2
        FAILURE_COUNT=$((FAILURE_COUNT+1))
    elif ! grep -Fq 'CONTINUOUS_3D_POSE_CAPTURE: EVIDENCE COMPLETE' \
        "$RESULT_ROOT/logs/pose_capture.log"; then
        echo "FAIL: pose-capture completion marker is absent" >&2
        FAILURE_COUNT=$((FAILURE_COUNT+1))
    else
        echo ">>> POSE CAPTURE TEST COMPLETE"
    fi
fi

if (( FAILURE_COUNT == 0 )); then
    echo ">>> DEFAULT LM REGRESSION START"
    "$EXECUTABLE" case=pose_recovery 2>&1 | tee "$RESULT_ROOT/logs/pose_recovery.log"
    REGRESSION_STATUS=${PIPESTATUS[0]}
    if (( REGRESSION_STATUS != 0 )); then
        echo "FAIL: default LM regression returned status $REGRESSION_STATUS" >&2
        FAILURE_COUNT=$((FAILURE_COUNT+1))
    else
        echo ">>> DEFAULT LM REGRESSION COMPLETE"
    fi
fi

if (( FAILURE_COUNT == 0 )); then
    echo ">>> POSE MECHANISM TEST START"
    "$EXECUTABLE" case=pose_capture_mechanism 2>&1 | tee "$RESULT_ROOT/logs/pose_mechanism.log"
    MECHANISM_STATUS=${PIPESTATUS[0]}
    if (( MECHANISM_STATUS != 0 )); then
        echo "FAIL: pose-mechanism test returned status $MECHANISM_STATUS" >&2
        FAILURE_COUNT=$((FAILURE_COUNT+1))
    elif ! grep -Fq 'CONTINUOUS_3D_POSE_MECHANISM: EVIDENCE COMPLETE' \
        "$RESULT_ROOT/logs/pose_mechanism.log"; then
        echo "FAIL: pose-mechanism completion marker is absent" >&2
        FAILURE_COUNT=$((FAILURE_COUNT+1))
    else
        echo ">>> POSE MECHANISM TEST COMPLETE"
    fi
fi

ANALYSIS_STATUS=2
if command -v python3 >/dev/null 2>&1; then
    python3 "$SCRIPT_DIR/analyze_pose_capture.py" --result-root "$RESULT_ROOT" \
        2>&1 | tee "$RESULT_ROOT/logs/analysis.log"
    ANALYSIS_STATUS=${PIPESTATUS[0]}
    if (( ANALYSIS_STATUS != 0 )); then
        echo "FAIL: pose-capture analysis returned status $ANALYSIS_STATUS" >&2
        FAILURE_COUNT=$((FAILURE_COUNT+1))
    fi
fi

{
    echo "test_status=$TEST_STATUS"
    echo "mechanism_status=$MECHANISM_STATUS"
    echo "regression_status=$REGRESSION_STATUS"
    echo "analysis_status=$ANALYSIS_STATUS"
    echo "recorded_failures=$FAILURE_COUNT"
} >"$RESULT_ROOT/status_details.txt"

cat >"$RESULT_ROOT/input/source_files.txt" <<'EOF'
production/tests/simple_continuous_3D_pcg_refinement_pose_capture_test.f90
production/tests/simple_continuous_3D_pcg_refinement_pose_mechanism_test.f90
production/tests/simple_test_continuous_3D_pcg_refinement.f90
src/main/volume/simple_reconstructor_pcg.f90
doc/implementation_notes/completed/continuous_3D_pose_end_polishing_validation_evidence.md
doc/implementation_notes/continuous_3D_refinement_on_pcg_operator.md
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
    echo "CONTINUOUS_3D_POSE_CAPTURE_PACKAGE: PASS"
else
    echo FAIL >"$RESULT_ROOT/STATUS.txt"
    echo "CONTINUOUS_3D_POSE_CAPTURE_PACKAGE: FAIL ($FAILURE_COUNT recorded failures)" >&2
fi
(cd "$RESULT_ROOT" && find . -type f ! -name MANIFEST.sha256 ! -name run.log -print0 | \
    sort -z | xargs -0 sha256sum) >"$RESULT_ROOT/MANIFEST.sha256"
echo "Download this one directory: $RESULT_ROOT"
echo "Review analysis/summary.md first."
(( FAILURE_COUNT == 0 ))
