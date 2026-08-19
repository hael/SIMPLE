#!/usr/bin/env bash
set -uo pipefail

usage() {
    cat <<'EOF'
Usage:
  reanalyze_oracle_pose_capture.sh /absolute/path/to/continuous_3D_pose_capture_TIMESTAMP

Reanalyze a completed pose-capture directory without rerunning the Fortran
test. The script requires the original numerical completion marker. It updates
the analysis, package status, copied analyzer, and checksums in that directory.
EOF
}

if [[ $# -ne 1 ]]; then
    usage >&2
    exit 2
fi
RESULT_ROOT=$(readlink -f "$1")
if [[ ! -d $RESULT_ROOT ]]; then
    echo "ERROR: result directory does not exist: $RESULT_ROOT" >&2
    exit 2
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$RESULT_ROOT"/{analysis,input,logs}
exec > >(tee -a "$RESULT_ROOT/reanalysis.log") 2>&1
echo "Reanalyzing pose-capture directory: $RESULT_ROOT"

FAILURE_COUNT=0
TEST_STATUS=1
if [[ -f $RESULT_ROOT/logs/pose_capture.log ]] && \
    grep -Fq 'CONTINUOUS_3D_POSE_CAPTURE: EVIDENCE COMPLETE' \
        "$RESULT_ROOT/logs/pose_capture.log"; then
    TEST_STATUS=0
    echo "PASS: original numerical completion marker is present"
else
    echo "FAIL: original numerical completion marker is absent" >&2
    FAILURE_COUNT=$((FAILURE_COUNT+1))
fi

ANALYSIS_STATUS=2
if ! command -v python3 >/dev/null 2>&1; then
    echo "FAIL: python3 is not available" >&2
    FAILURE_COUNT=$((FAILURE_COUNT+1))
else
    python3 "$SCRIPT_DIR/analyze_pose_capture.py" --result-root "$RESULT_ROOT" \
        2>&1 | tee "$RESULT_ROOT/logs/analysis_recheck.log"
    ANALYSIS_STATUS=${PIPESTATUS[0]}
    if (( ANALYSIS_STATUS != 0 )); then
        echo "FAIL: pose-capture analysis returned status $ANALYSIS_STATUS" >&2
        FAILURE_COUNT=$((FAILURE_COUNT+1))
    fi
fi

HAS_SHA256=yes
if ! command -v sha256sum >/dev/null 2>&1; then
    echo "FAIL: sha256sum is unavailable; MANIFEST.sha256 cannot be refreshed" >&2
    FAILURE_COUNT=$((FAILURE_COUNT+1))
    HAS_SHA256=no
fi

cp "$SCRIPT_DIR/analyze_pose_capture.py" "$RESULT_ROOT/input/"
cp "$SCRIPT_DIR/reanalyze_oracle_pose_capture.sh" "$RESULT_ROOT/input/"
{
    echo "reanalysis=yes"
    echo "test_status=$TEST_STATUS"
    echo "analysis_status=$ANALYSIS_STATUS"
    echo "recorded_failures=$FAILURE_COUNT"
    echo "python_version=$(python3 --version 2>&1 || true)"
} >"$RESULT_ROOT/status_details.txt"

if (( FAILURE_COUNT == 0 )); then
    echo PASS >"$RESULT_ROOT/STATUS.txt"
    echo "CONTINUOUS_3D_POSE_CAPTURE_PACKAGE: PASS"
else
    echo FAIL >"$RESULT_ROOT/STATUS.txt"
    echo "CONTINUOUS_3D_POSE_CAPTURE_PACKAGE: FAIL ($FAILURE_COUNT recorded failures)" >&2
fi

if [[ $HAS_SHA256 == yes ]]; then
    (cd "$RESULT_ROOT" && find . -type f ! -name MANIFEST.sha256 \
        ! -name run.log ! -name reanalysis.log -print0 | \
        sort -z | xargs -0 sha256sum) >"$RESULT_ROOT/MANIFEST.sha256"
fi

echo "Review analysis/summary.md first."
(( FAILURE_COUNT == 0 ))
