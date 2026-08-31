#!/usr/bin/env python3
"""Analyze the complete continuous 3-D pose-refinement Oracle package."""

import argparse
import csv
import json
import os
import re
import subprocess
import sys


REQUIRED_TESTS = [
    "configuration_matrix", "pose_mother", "fixed_reference", "forward_path",
    "matched_window", "reference_bias", "operator_contract",
    "pose_capture_range", "pose_capture_mechanism", "tolerance_calibration",
    "objective_normals", "lm_transactions", "ctf_sigma", "forward_hierarchy",
    "cartesian_mother", "pcg_mother", "pcg_recon", "pcg_priors",
]

EXPECTED_MARKERS = {
    "configuration_matrix": "POSE_CONT_REFINEMENT_CONFIGURATION_MATRIX: PASS",
    "pose_mother": "Continuous 3D pose refinement test suite: PASS",
    "fixed_reference": "CONTINUOUS_3D_FIXED_REFERENCE: EVIDENCE COMPLETE",
    "forward_path": "CONTINUOUS_3D_FORWARD_PATH: EVIDENCE COMPLETE",
    "matched_window": "CONTINUOUS_3D_MATCHED_WINDOW: EVIDENCE COMPLETE",
    "reference_bias": "CONTINUOUS_3D_REFERENCE_BIAS: EVIDENCE COMPLETE",
    "operator_contract": "CONTINUOUS_3D_OPERATOR_CONTRACT: EVIDENCE COMPLETE",
    "tolerance_calibration": "POSE_CONT_REFINEMENT_FORWARD_AMENDMENT: FROZEN",
    "objective_normals": "POSE_CONT_REFINEMENT_OBJECTIVE_NORMALS: PASS",
    "lm_transactions": "POSE_CONT_REFINEMENT_LM_TRANSACTIONS: PASS",
    "ctf_sigma": "POSE_CONT_REFINEMENT_CTF_SIGMA: PASS",
    "forward_hierarchy": "POSE_CONT_REFINEMENT_FORWARD_HIERARCHY: PASS",
    "cartesian_mother": "Neutral Cartesian Fourier test suite: PASS",
}

KNOWN_PCG = {
    "lambda": [10.0, 10.0],
    "pcg_raw_l2": [0.6612755, 0.6608928],
    "grid_raw_l2": [0.6518738, 0.6388823],
}


def read_text(path):
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def close_pair(actual, expected, tolerance=5.0e-7):
    return len(actual) == 2 and all(
        abs(got - want) <= tolerance for got, want in zip(actual, expected)
    )


def extract_pair(text, label):
    pattern = re.compile(re.escape(label) + r"\s*([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)")
    matches = pattern.findall(text)
    if not matches:
        return []
    return [float(value) for value in matches[-1]]


def count_tsv_rows(path):
    with open(path, "r", encoding="utf-8", newline="") as handle:
        return max(sum(1 for _ in csv.reader(handle, delimiter="\t")) - 1, 0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-root", required=True)
    args = parser.parse_args()
    root = os.path.realpath(args.result_root)
    checks = []

    def add(name, passed, detail):
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    statuses_path = os.path.join(root, "statuses.tsv")
    statuses = {}
    if os.path.isfile(statuses_path):
        with open(statuses_path, "r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                statuses[row["test_id"]] = int(row["exit_status"])
    add("required_status_rows", sorted(statuses) == sorted(REQUIRED_TESTS),
        "observed={}".format(",".join(sorted(statuses))))

    for test_id in REQUIRED_TESTS:
        if test_id == "pcg_mother":
            continue
        expected_status = 0
        actual_status = statuses.get(test_id)
        add("{}_exit".format(test_id), actual_status == expected_status,
            "expected={} actual={}".format(expected_status, actual_status))

    for test_id, marker in EXPECTED_MARKERS.items():
        path = os.path.join(root, "tests", test_id, "run.log")
        text = read_text(path) if os.path.isfile(path) else ""
        add("{}_marker".format(test_id), marker in text, marker)

    capture_log = os.path.join(root, "capture", "logs", "pose_capture.log")
    mechanism_log = os.path.join(root, "capture", "logs", "pose_mechanism.log")
    add("pose_capture_marker", os.path.isfile(capture_log) and
        "CONTINUOUS_3D_POSE_CAPTURE: EVIDENCE COMPLETE" in read_text(capture_log),
        "pose capture evidence marker")
    add("pose_mechanism_marker", os.path.isfile(mechanism_log) and
        "CONTINUOUS_3D_POSE_MECHANISM: EVIDENCE COMPLETE" in read_text(mechanism_log),
        "pose mechanism evidence marker")

    capture_analyzer = os.path.join(root, "source", "production", "tests",
                                    "pose_cont_capture", "analyze_pose_capture.py")
    capture_status = 2
    if os.path.isfile(capture_analyzer):
        completed = subprocess.run(
            [sys.executable, capture_analyzer, "--result-root", os.path.join(root, "capture")],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        capture_status = completed.returncode
        with open(os.path.join(root, "analysis", "pose_capture_analysis.log"), "w",
                  encoding="utf-8") as handle:
            handle.write(completed.stdout)
    add("pose_capture_analysis", capture_status == 0,
        "analyzer_exit={}".format(capture_status))

    config_path = os.path.join(root, "tests", "configuration_matrix", "evidence",
                               "configuration_matrix.tsv")
    config_rows = count_tsv_rows(config_path) if os.path.isfile(config_path) else -1
    add("configuration_matrix_rows", config_rows == 1296,
        "expected=1296 actual={}".format(config_rows))

    frozen_path = os.path.join(root, "tests", "tolerance_calibration", "evidence",
                               "FORWARD_AMENDMENT.FROZEN")
    add("forward_amendment_frozen", os.path.isfile(frozen_path), frozen_path)

    required_artifacts = {
        "configuration_matrix": ["configuration_matrix.tsv", "configuration_matrix_manifest.tsv"],
        "tolerance_calibration": ["calibration_raw_errors.tsv", "frozen_tolerances.tsv",
                                  "fixture_manifest.tsv", "FORWARD_AMENDMENT.FROZEN"],
        "objective_normals": ["residual_jacobian_components.tsv",
                              "normal_equation_components.tsv",
                              "objective_normal_summary.tsv", "blas_normal_equations.tsv"],
        "lm_transactions": ["lm_systems.tsv", "pose_transactions.tsv"],
        "ctf_sigma": ["ctf_sigma_cases.tsv", "variance_outcomes.tsv", "pcg_fallback.tsv"],
        "forward_hierarchy": ["forward_hierarchy.tsv", "forward_hierarchy_components.tsv",
                              "forward_hierarchy_manifest.tsv", "forward_holdout_fixtures.tsv"],
    }
    for test_id, filenames in required_artifacts.items():
        missing = []
        for filename in filenames:
            path = os.path.join(root, "tests", test_id, "evidence", filename)
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                missing.append(filename)
        add("{}_artifacts".format(test_id), not missing,
            "missing={}".format(",".join(missing) if missing else "none"))

    pcg_path = os.path.join(root, "tests", "pcg_mother", "run.log")
    pcg_text = read_text(pcg_path) if os.path.isfile(pcg_path) else ""
    pcg_status = statuses.get("pcg_mother")
    pcg_lambda = extract_pair(pcg_text, "CONTINUOUS_3D_MATRIX best noisy lambda even/odd:")
    pcg_raw = extract_pair(pcg_text, "CONTINUOUS_3D_MATRIX best noisy raw L2 even/odd:")
    grid_rows = re.findall(
        r"CONTINUOUS_3D_MATRIX conventional raw L2 errors:\s*"
        r"([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)",
        pcg_text)
    grid_raw = [float(grid_rows[-1][2]), float(grid_rows[-1][3])] if grid_rows else []
    known_failure = (
        pcg_status == 1
        and close_pair(pcg_lambda, KNOWN_PCG["lambda"], 1.0e-8)
        and close_pair(pcg_raw, KNOWN_PCG["pcg_raw_l2"])
        and close_pair(grid_raw, KNOWN_PCG["grid_raw_l2"])
        and "the best scale-sensitive PCG reconstruction does not beat conventional gridding" in pcg_text
    )
    add("pcg_mother_phase2_disposition", known_failure,
        "status={} lambda={} pcg_raw={} grid_raw={}".format(
            pcg_status, pcg_lambda, pcg_raw, grid_raw))

    source_hash_path = os.path.join(root, "input", "source.sha256")
    source_root = os.path.join(root, "source")
    source_hash_ok = os.path.isfile(source_hash_path)
    if source_hash_ok:
        completed = subprocess.run(["sha256sum", "-c", source_hash_path], cwd=source_root,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   universal_newlines=True)
        source_hash_ok = completed.returncode == 0
        with open(os.path.join(root, "analysis", "source_hash_check.log"), "w",
                  encoding="utf-8") as handle:
            handle.write(completed.stdout)
    add("source_snapshot_integrity", source_hash_ok, source_hash_path)

    command_text = ""
    for directory, _, files in os.walk(root):
        for filename in files:
            if filename.endswith("command.txt"):
                command_text += read_text(os.path.join(directory, filename))
    add("removed_route_absent", "pcg_pose_polish" not in command_text and
        "simple_pcg_pose_polisher" not in command_text,
        "command manifests contain no deleted route")

    environment = read_text(os.path.join(root, "input", "environment.txt"))
    working_match = re.search(r"^runtime_working_directory=(.+)$", environment, re.MULTILINE)
    working_directory = working_match.group(1) if working_match else ""
    add("runtime_parent_is_projects", os.path.basename(working_directory) == "Projects",
        working_directory)

    passed = all(check["passed"] for check in checks)
    summary = {
        "schema_version": 1,
        "status": "PASS" if passed else "FAIL",
        "known_pcg_reconstruction_failure_retained": known_failure,
        "checks": checks,
    }
    with open(os.path.join(root, "analysis", "summary.json"), "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    with open(os.path.join(root, "analysis", "summary.tsv"), "w", encoding="utf-8") as handle:
        handle.write("check\tstatus\tdetail\n")
        for check in checks:
            detail = check["detail"].replace("\t", " ").replace("\n", " ")
            handle.write("{}\t{}\t{}\n".format(
                check["name"], "PASS" if check["passed"] else "FAIL", detail))

    print("Continuous 3-D pose package analysis: {}".format(summary["status"]))
    for check in checks:
        print("{}: {} ({})".format(
            check["name"], "PASS" if check["passed"] else "FAIL", check["detail"]))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
