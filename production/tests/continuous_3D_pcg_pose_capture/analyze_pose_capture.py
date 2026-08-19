#!/usr/bin/env python3
"""Validate and summarize one isolated PCG Cartesian pose-capture sweep."""

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Tuple


LEGACY_SHIFT_MAGNITUDES = (0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0)
V2_SHIFT_MAGNITUDES = (0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 12.5, 15.0)
V3_SHIFT_MAGNITUDES = (
    0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 10.0, 12.5, 15.0
)
ROTATION_MAGNITUDES = (0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 12.5, 15.0)
V2_JOINT_CASES = ((2.0, 0.5), (5.0, 1.0), (10.0, 3.0), (15.0, 5.0))
V3_SMALL_JOINT_CASES = ((2.0, 0.5), (5.0, 1.0))
V3_JOINT_ROTATIONS = (10.0, 12.5, 15.0)
V3_JOINT_SHIFTS = (3.0, 4.0, 5.0)
V3_MIXED_CASES = (
    (5.0, 0.0, 0.0, -1.0, 0.0),
    (-10.0, 0.0, 0.0, 3.0, 0.0),
    (-15.0, 0.0, 0.0, 5.0, 0.0),
)
V3_MULTI_AXIS_CASES = (
    (2.0, -3.0, 4.0, 0.5, -1.0),
    (5.0, -7.5, 10.0, 2.0, -3.0),
    (-5.0, 7.5, -10.0, -2.0, 3.0),
)
MECHANISM_VOLUMES = ("gaussian_blobs", "shifted_mixture", "permuted_texture")
MECHANISM_CASES = {
    "joint_p10_p5": (10.0, 5.0),
    "joint_p12p5_p5": (12.5, 5.0),
    "joint_p15_p5": (15.0, 5.0),
    "joint_p15_p4": (15.0, 4.0),
    "joint_m15_p5": (15.0, 5.0),
}
MECHANISM_ROUTES = (
    "joint",
    "guarded_joint",
    "shift_then_joint",
    "rotation_then_joint",
)
TRUTH_PATH_ROTATION_TOL_DEG = 1.0e-6
ENDPOINT_OBJECTIVE_ABS_TOL = 1.0e-11
ENDPOINT_OBJECTIVE_REL_TOL = 1.0e-6
MECHANISM_FLOAT_COLUMNS = {
    "final_rotation_x_deg",
    "final_rotation_y_deg",
    "final_rotation_z_deg",
    "final_shift_x_px",
    "final_shift_y_px",
    "final_rotation_norm_deg",
    "final_shift_norm_px",
    "objective_before",
    "objective_after",
}
MECHANISM_INT_COLUMNS = {
    "accepted_steps",
    "attempted_steps",
    "within_15deg_5px_from_seed",
}
VECTOR_COLUMNS = {
    "injected_rotation_x_deg",
    "injected_rotation_y_deg",
    "injected_rotation_z_deg",
    "injected_shift_x_px",
    "injected_shift_y_px",
}
REQUIRED_COLUMNS = {
    "trial_id",
    "family",
    "axis",
    "sign",
    "injected_magnitude",
    "initial_rotation_deg",
    "final_rotation_deg",
    "initial_shift_px",
    "final_shift_px",
    "objective_before",
    "objective_after",
    "accepted_steps",
    "attempted_steps",
    "status_name",
    "max_rotation_step_deg",
    "max_shift_step_px",
    "stencil_switches",
    "monotone",
}
FLOAT_COLUMNS = {
    "injected_magnitude",
    "injected_rotation_deg",
    "injected_shift_px",
    "initial_rotation_deg",
    "final_rotation_deg",
    "initial_shift_px",
    "final_shift_px",
    "objective_before",
    "objective_after",
    "max_rotation_step_deg",
    "max_shift_step_px",
}
INT_COLUMNS = {
    "axis",
    "sign",
    "accepted_steps",
    "attempted_steps",
    "status",
    "stencil_switches",
    "monotone",
    "rotation_improved",
    "shift_improved",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-root", type=Path, required=True)
    return parser.parse_args()


def read_configuration(path: Path) -> Dict[str, str]:
    values = {}  # type: Dict[str, str]
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def read_rows(
    path: Path, errors: List[str], scenario_version: str
) -> List[Dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = set(REQUIRED_COLUMNS)
        if scenario_version in ("3", "4"):
            required.update(VECTOR_COLUMNS)
        missing = required.difference(reader.fieldnames or [])
        if missing:
            errors.append("CSV is missing columns: " + ", ".join(sorted(missing)))
        rows = []  # type: List[Dict[str, Any]]
        for line_number, raw in enumerate(reader, start=2):
            row = dict(raw)  # type: Dict[str, Any]
            try:
                float_columns = set(FLOAT_COLUMNS)
                if scenario_version in ("3", "4"):
                    float_columns.update(VECTOR_COLUMNS)
                for name in float_columns:
                    row[name] = float(raw[name])
                for name in INT_COLUMNS:
                    row[name] = int(raw[name])
            except (KeyError, TypeError, ValueError) as exc:
                errors.append(f"CSV line {line_number} is not parseable: {exc}")
                continue
            if scenario_version in ("1", "2"):
                add_legacy_vectors(row)
            rows.append(row)
    return rows


def read_typed_csv(
    path: Path, float_columns: set, int_columns: set, errors: List[str], label: str
) -> List[Dict[str, Any]]:
    rows = []  # type: List[Dict[str, Any]]
    if not path.is_file():
        errors.append(f"missing required artifact: {path.name}")
        return rows
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = set(float_columns).union(int_columns)
        missing = required.difference(reader.fieldnames or [])
        if missing:
            errors.append(f"{label} CSV is missing columns: " + ", ".join(sorted(missing)))
        for line_number, raw in enumerate(reader, start=2):
            row = dict(raw)  # type: Dict[str, Any]
            try:
                for name in float_columns:
                    row[name] = float(raw[name])
                for name in int_columns:
                    row[name] = int(raw[name])
            except (KeyError, TypeError, ValueError) as exc:
                errors.append(f"{label} line {line_number} is not parseable: {exc}")
                continue
            rows.append(row)
    return rows


def add_legacy_vectors(row: Dict[str, Any]) -> None:
    rotation = [0.0, 0.0, 0.0]
    shift = [0.0, 0.0]
    axis = int(row["axis"])
    sign = int(row["sign"])
    if row["family"] == "rotation" and 1 <= axis <= 3:
        rotation[axis - 1] = sign * float(row["injected_rotation_deg"])
    elif row["family"] == "shift" and 1 <= axis <= 2:
        shift[axis - 1] = sign * float(row["injected_shift_px"])
    elif row["family"] == "joint":
        rotation[0] = float(row["injected_rotation_deg"])
        shift[0] = float(row["injected_shift_px"])
    row["injected_rotation_x_deg"] = rotation[0]
    row["injected_rotation_y_deg"] = rotation[1]
    row["injected_rotation_z_deg"] = rotation[2]
    row["injected_shift_x_px"] = shift[0]
    row["injected_shift_y_px"] = shift[1]


def scenario(
    family: str, axis: int, sign: int, rotation: Tuple[float, float, float],
    shift: Tuple[float, float]
) -> Tuple[Any, ...]:
    return (
        family,
        axis,
        sign,
        round(rotation[0], 6),
        round(rotation[1], 6),
        round(rotation[2], 6),
        round(shift[0], 6),
        round(shift[1], 6),
    )


def expected_scenario_matrix(configuration: Dict[str, str]) -> set:
    version = configuration.get("scenario_version", "1")
    if version in ("3", "4"):
        shift_magnitudes = V3_SHIFT_MAGNITUDES
    elif version == "2":
        shift_magnitudes = V2_SHIFT_MAGNITUDES
    else:
        shift_magnitudes = LEGACY_SHIFT_MAGNITUDES
    scenarios = {scenario("exact", 0, 0, (0.0, 0.0, 0.0), (0.0, 0.0))}
    for axis in (1, 2):
        for sign in (-1, 1):
            for magnitude in shift_magnitudes:
                vector = [0.0, 0.0]
                vector[axis - 1] = sign * magnitude
                scenarios.add(scenario("shift", axis, sign, (0.0, 0.0, 0.0), tuple(vector)))
    for axis in (1, 2, 3):
        for sign in (-1, 1):
            for magnitude in ROTATION_MAGNITUDES:
                vector = [0.0, 0.0, 0.0]
                vector[axis - 1] = sign * magnitude
                scenarios.add(scenario("rotation", axis, sign, tuple(vector), (0.0, 0.0)))
    if version == "2":
        for rotation_deg, shift_px in V2_JOINT_CASES:
            scenarios.add(scenario("joint", 1, 1, (rotation_deg, 0.0, 0.0), (shift_px, 0.0)))
    elif version in ("3", "4"):
        for rotation_deg, shift_px in V3_SMALL_JOINT_CASES:
            scenarios.add(scenario("joint", 1, 1, (rotation_deg, 0.0, 0.0), (shift_px, 0.0)))
        for rotation_deg in V3_JOINT_ROTATIONS:
            for shift_px in V3_JOINT_SHIFTS:
                scenarios.add(scenario("joint", 1, 1, (rotation_deg, 0.0, 0.0), (shift_px, 0.0)))
        for values in V3_MIXED_CASES:
            scenarios.add(scenario("joint", 1, 0, values[:3], values[3:]))
        for values in V3_MULTI_AXIS_CASES:
            scenarios.add(scenario("joint", 0, 0, values[:3], values[3:]))
    return scenarios


def validate_rows(
    rows: List[Dict[str, Any]], configuration: Dict[str, str], artifacts: Path
) -> List[str]:
    errors = []  # type: List[str]
    if configuration.get("scenario_version", "1") not in ("1", "2", "3", "4"):
        errors.append("configuration has an unsupported scenario version")
    expected_scenarios = expected_scenario_matrix(configuration)
    expected_trials = len(expected_scenarios)
    if len(rows) != expected_trials:
        errors.append(f"expected {expected_trials} trials, found {len(rows)}")
    trial_ids = [str(row["trial_id"]) for row in rows]
    if len(set(trial_ids)) != len(trial_ids):
        errors.append("trial identifiers are not unique")
    actual_scenarios = {
        scenario(
            str(row["family"]),
            int(row["axis"]),
            int(row["sign"]),
            (
                float(row["injected_rotation_x_deg"]),
                float(row["injected_rotation_y_deg"]),
                float(row["injected_rotation_z_deg"]),
            ),
            (float(row["injected_shift_x_px"]), float(row["injected_shift_y_px"])),
        )
        for row in rows
    }
    if actual_scenarios != expected_scenarios:
        missing = expected_scenarios.difference(actual_scenarios)
        extra = actual_scenarios.difference(expected_scenarios)
        if missing:
            errors.append(f"scenario matrix is missing {len(missing)} trial(s)")
        if extra:
            errors.append(f"scenario matrix has {len(extra)} unexpected trial(s)")

    rotation_bound = float(configuration.get("rotation_step_bound_deg", "nan"))
    shift_bound = float(configuration.get("shift_step_bound_px", "nan"))
    if not math.isfinite(rotation_bound) or not math.isfinite(shift_bound):
        errors.append("configuration does not contain finite step bounds")

    exact_rows = [row for row in rows if row["family"] == "exact"]
    if len(exact_rows) != 1:
        errors.append(f"expected one exact-pose control, found {len(exact_rows)}")
    elif (
        exact_rows[0]["status_name"] != "finite_no_improvement"
        or int(exact_rows[0]["accepted_steps"]) != 0
        or float(exact_rows[0]["final_rotation_deg"]) > 1.0e-8
        or float(exact_rows[0]["final_shift_px"]) > 1.0e-10
    ):
        errors.append("exact-pose control was not stationary")

    for row in rows:
        trial_id = str(row["trial_id"])
        numeric = [float(row[name]) for name in FLOAT_COLUMNS]
        if not all(math.isfinite(value) for value in numeric):
            errors.append(f"{trial_id}: non-finite metric")
        if int(row["monotone"]) != 1:
            errors.append(f"{trial_id}: accepted objective trace is not monotone")
        if int(row["attempted_steps"]) < int(row["accepted_steps"]):
            errors.append(f"{trial_id}: attempted steps are less than accepted steps")
        if str(row["status_name"]) == "invalid_numerics":
            errors.append(f"{trial_id}: optimizer returned invalid_numerics")
        if float(row["max_rotation_step_deg"]) > rotation_bound + 1.0e-6:
            errors.append(f"{trial_id}: rotation step exceeded its bound")
        if float(row["max_shift_step_px"]) > shift_bound + 1.0e-8:
            errors.append(f"{trial_id}: shift step exceeded its bound")
        before = float(row["objective_before"])
        after = float(row["objective_after"])
        if after > before + 1.0e-10 * max(1.0, abs(before)):
            errors.append(f"{trial_id}: final objective exceeds initial objective")
        for suffix in (
            "initial_prediction.mrc",
            "initial_residual.mrc",
            "final_prediction.mrc",
            "final_residual.mrc",
        ):
            artifact = artifacts / f"{trial_id}_{suffix}"
            if not artifact.is_file() or artifact.stat().st_size <= 1024:
                errors.append(f"{trial_id}: missing or empty {suffix}")
    return errors


def read_mechanism_evidence(
    artifacts: Path, errors: List[str]
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    summary_rows = read_typed_csv(
        artifacts / "pose_mechanism_summary.csv",
        MECHANISM_FLOAT_COLUMNS,
        MECHANISM_INT_COLUMNS,
        errors,
        "mechanism summary",
    )
    trajectory_float = {
        "rotation_x_deg", "rotation_y_deg", "rotation_z_deg", "shift_x_px",
        "shift_y_px", "rotation_norm_deg", "shift_norm_px", "objective",
    }
    trajectory_rows = read_typed_csv(
        artifacts / "pose_mechanism_trajectory.csv",
        trajectory_float,
        {"accepted_step"},
        errors,
        "mechanism trajectory",
    )
    path_float = {
        "alpha", "rotation_x_deg", "rotation_y_deg", "rotation_z_deg",
        "shift_x_px", "shift_y_px", "rotation_norm_deg", "shift_norm_px", "objective",
    }
    path_rows = read_typed_csv(
        artifacts / "pose_objective_paths.csv",
        path_float,
        set(),
        errors,
        "objective path",
    )
    return summary_rows, trajectory_rows, path_rows


def validate_mechanism_evidence(
    summary_rows: List[Dict[str, Any]], trajectory_rows: List[Dict[str, Any]],
    path_rows: List[Dict[str, Any]], artifacts: Path
) -> List[str]:
    errors = []  # type: List[str]
    expected_summary = {
        (volume, case, route)
        for volume in MECHANISM_VOLUMES
        for case in MECHANISM_CASES
        for route in MECHANISM_ROUTES
    }
    actual_summary = {
        (str(row.get("volume")), str(row.get("case")), str(row.get("route")))
        for row in summary_rows
    }
    if actual_summary != expected_summary or len(summary_rows) != len(expected_summary):
        errors.append(
            f"mechanism summary expected {len(expected_summary)} unique rows, "
            f"found {len(summary_rows)}"
        )

    for row in summary_rows:
        trial = "{}/{}/{}".format(row.get("volume"), row.get("case"), row.get("route"))
        numeric = [float(row[name]) for name in MECHANISM_FLOAT_COLUMNS]
        if not all(math.isfinite(value) for value in numeric):
            errors.append(f"{trial}: mechanism summary contains a non-finite value")
        if int(row["attempted_steps"]) < int(row["accepted_steps"]):
            errors.append(f"{trial}: attempted steps are less than accepted steps")
        if float(row["objective_after"]) > float(row["objective_before"]) + 1.0e-12:
            errors.append(f"{trial}: staged route increased the objective")
        if str(row.get("status_name")) in ("invalid_numerics", "unknown"):
            errors.append(f"{trial}: invalid terminal status")
        if row.get("route") == "guarded_joint" and int(row["within_15deg_5px_from_seed"]) != 1:
            errors.append(f"{trial}: guarded route escaped its cumulative bound")

    trajectory_groups = defaultdict(list)
    for row in trajectory_rows:
        key = (str(row.get("volume")), str(row.get("case")), str(row.get("route")))
        trajectory_groups[key].append(row)
        values = [
            float(row[name])
            for name in (
                "rotation_x_deg", "rotation_y_deg", "rotation_z_deg", "shift_x_px",
                "shift_y_px", "rotation_norm_deg", "shift_norm_px", "objective",
            )
        ]
        if not all(math.isfinite(value) for value in values):
            errors.append("{}/{}/{}: non-finite trajectory value".format(*key))
    for key in expected_summary:
        if key not in trajectory_groups:
            errors.append("{}/{}/{}: missing accepted-pose trajectory".format(*key))
            continue
        expected_stages = {
            "joint": {"joint"},
            "guarded_joint": {"joint"},
            "shift_then_joint": {"shift_only", "joint"},
            "rotation_then_joint": {"rotation_only", "joint"},
        }[key[2]]
        group = trajectory_groups[key]
        if {str(row.get("stage")) for row in group} != expected_stages:
            errors.append("{}/{}/{}: staged trajectory is incomplete".format(*key))
        for stage in expected_stages:
            stage_rows = [row for row in group if row.get("stage") == stage]
            steps = [int(row["accepted_step"]) for row in stage_rows]
            if steps != list(range(len(stage_rows))):
                errors.append("{}/{}/{}: {} step indices are not contiguous".format(*key, stage))
            objectives = [float(row["objective"]) for row in stage_rows]
            if any(right > left + 1.0e-12 for left, right in zip(objectives, objectives[1:])):
                errors.append("{}/{}/{}: {} objective trace is not monotone".format(*key, stage))

    summary_by_key = {
        (str(row["volume"]), str(row["case"]), str(row["route"])): row
        for row in summary_rows
    }
    for key, group in trajectory_groups.items():
        if key not in summary_by_key or not group:
            continue
        final = group[-1]
        summary = summary_by_key[key]
        if abs(float(final["rotation_norm_deg"]) - float(summary["final_rotation_norm_deg"])) > 1.0e-8:
            errors.append("{}/{}/{}: final rotation does not match trajectory".format(*key))
        if abs(float(final["shift_norm_px"]) - float(summary["final_shift_norm_px"])) > 1.0e-8:
            errors.append("{}/{}/{}: final shift does not match trajectory".format(*key))
        if abs(float(final["objective"]) - float(summary["objective_after"])) > 1.0e-12:
            errors.append("{}/{}/{}: final objective does not match trajectory".format(*key))

    expected_path_groups = {
        (volume, case, target)
        for volume in MECHANISM_VOLUMES
        for case in MECHANISM_CASES
        for target in ("truth", "joint_endpoint")
    }
    path_groups = defaultdict(list)
    for row in path_rows:
        key = (str(row.get("volume")), str(row.get("case")), str(row.get("target")))
        path_groups[key].append(row)
        if not all(
            math.isfinite(float(row[name]))
            for name in (
                "alpha", "rotation_x_deg", "rotation_y_deg", "rotation_z_deg",
                "shift_x_px", "shift_y_px", "rotation_norm_deg", "shift_norm_px", "objective",
            )
        ):
            errors.append("{}/{}/{}: non-finite objective-path value".format(*key))
    if set(path_groups) != expected_path_groups:
        errors.append("objective-path scenario matrix is incomplete")
    for key, group in path_groups.items():
        ordered = sorted(group, key=lambda row: float(row["alpha"]))
        if len(ordered) != 41 or abs(float(ordered[0]["alpha"])) > 1.0e-12 or abs(
            float(ordered[-1]["alpha"]) - 1.0
        ) > 1.0e-12:
            errors.append("{}/{}/{}: objective path does not contain 41 endpoints".format(*key))
        if key[2] == "truth" and ordered and float(ordered[-1]["objective"]) > 1.0e-12:
            errors.append("{}/{}/{}: truth endpoint objective is not zero".format(*key))
        if key[2] == "truth" and ordered and (
            float(ordered[-1]["rotation_norm_deg"]) > TRUTH_PATH_ROTATION_TOL_DEG
            or float(ordered[-1]["shift_norm_px"]) > 1.0e-10
        ):
            errors.append("{}/{}/{}: truth path does not end at the truth pose".format(*key))
        if key[2] == "joint_endpoint" and ordered:
            direct = summary_by_key.get((key[0], key[1], "joint"))
            if direct is not None:
                direct_objective = float(direct["objective_after"])
                objective_tolerance = max(
                    ENDPOINT_OBJECTIVE_ABS_TOL,
                    ENDPOINT_OBJECTIVE_REL_TOL * abs(direct_objective),
                )
                if abs(float(ordered[-1]["objective"]) - direct_objective) > objective_tolerance:
                    errors.append("{}/{}/{}: endpoint path does not match direct joint LM".format(*key))

    for volume in MECHANISM_VOLUMES:
        for suffix in ("volume.mrc", "observation.mrc"):
            path = artifacts / f"mechanism_{volume}_{suffix}"
            if not path.is_file() or path.stat().st_size <= 1024:
                errors.append(f"missing or empty mechanism artifact: {path.name}")
    return errors


def target_values(row: Dict[str, Any]) -> Tuple[float, float]:
    if row["family"] == "rotation":
        return float(row["initial_rotation_deg"]), float(row["final_rotation_deg"])
    return float(row["initial_shift_px"]), float(row["final_shift_px"])


def scientific_recovery_warnings(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    warnings = []  # type: List[Dict[str, Any]]
    for row in rows:
        if row["family"] == "exact":
            continue
        reasons = []  # type: List[str]
        initial_rotation = float(row["initial_rotation_deg"])
        final_rotation = float(row["final_rotation_deg"])
        initial_shift = float(row["initial_shift_px"])
        final_shift = float(row["final_shift_px"])
        if initial_rotation > 1.0e-12:
            if final_rotation / initial_rotation > 0.1:
                reasons.append("rotation retained more than 10% of its initial error")
        elif final_rotation > 0.1:
            reasons.append("rotation coupling exceeded 0.1 degree")
        if initial_shift > 1.0e-12:
            if final_shift / initial_shift > 0.1:
                reasons.append("shift retained more than 10% of its initial error")
        elif final_shift > 0.1:
            reasons.append("shift coupling exceeded 0.1 pixel")
        if reasons:
            warnings.append(
                {
                    "trial_id": str(row["trial_id"]),
                    "initial_rotation_deg": initial_rotation,
                    "final_rotation_deg": final_rotation,
                    "initial_shift_px": initial_shift,
                    "final_shift_px": final_shift,
                    "reasons": reasons,
                }
            )
    return warnings


def vector_text(row: Dict[str, Any]) -> Tuple[str, str]:
    rotation = (
        float(row["injected_rotation_x_deg"]),
        float(row["injected_rotation_y_deg"]),
        float(row["injected_rotation_z_deg"]),
    )
    shift = (
        float(row["injected_shift_x_px"]),
        float(row["injected_shift_y_px"]),
    )
    return (
        "({:.4g}, {:.4g}, {:.4g})".format(*rotation),
        "({:.4g}, {:.4g})".format(*shift),
    )


def markdown_summary(
    rows: List[Dict[str, Any]], configuration: Dict[str, str], errors: List[str],
    recovery_warnings: List[Dict[str, Any]], mechanism_rows: List[Dict[str, Any]],
    path_rows: List[Dict[str, Any]]
) -> str:
    expected_trials = len(expected_scenario_matrix(configuration))
    version = configuration.get("scenario_version", "1")
    if version in ("3", "4"):
        scenario_description = (
            "- Separate sweeps bracket each coordinate. Joint trials include a 3 by 3 boundary grid, "
            "mixed signs, and simultaneous multi-axis rotations."
        )
    elif version == "2":
        scenario_description = (
            "- Separate sweeps change one coordinate. Four joint trials change one x rotation and "
            "x shift together."
        )
    else:
        scenario_description = (
            "- Each trial changes one coordinate while the LM solve can update all five coordinates."
        )
    lines = [
        "# PCG Cartesian pose-capture diagnostic",
        "",
        f"**Package integrity:** {'FAIL' if errors else 'PASS'}",
        f"**Scientific recovery warnings:** {len(recovery_warnings)}",
        "",
        "This result describes one clean, matched-operator numerical experiment. "
        "It does not validate `refine3D` integration or real-data improvement.",
        "",
        "## Fixed conditions",
        "",
        f"- Volume box: {configuration.get('box', 'unknown')}",
        f"- Sampling distance: {configuration.get('sampling_distance_A', 'unknown')} A/pixel",
        f"- Fourier shells: {configuration.get('shell_min', 'unknown')} to "
        f"{configuration.get('shell_max', 'unknown')}",
        f"- Nominal low-pass: {configuration.get('nominal_lowpass_A', 'unknown')} A",
        f"- LM rotation-step bound: {configuration.get('rotation_step_bound_deg', 'unknown')} degrees",
        f"- LM shift-step bound: {configuration.get('shift_step_bound_px', 'unknown')} pixels",
        f"- Maximum LM iterations: {configuration.get('max_lm_iterations', 'unknown')}",
        scenario_description,
        "",
        "## Integrity checks",
        "",
        f"- Trial rows: {len(rows)} of {expected_trials} expected",
        "- Exact-pose stationarity, finite metrics, monotone accepted objectives, and step bounds "
        "are package gates.",
    ]
    if errors:
        lines.extend(["", "### Integrity errors", ""])
        lines.extend(f"- {error}" for error in errors)

    lines.extend(
        [
            "",
            "## Scientific recovery warnings",
            "",
            "These warnings are not package-integrity failures. They identify gross capture loss. "
            "A warning occurs when a nonzero truth error retains more than 10 percent of its "
            "initial value, or when an initially zero coupled error exceeds 0.1 degree or 0.1 pixel.",
            "The limits are diagnostic thresholds, not production acceptance criteria.",
        ]
    )
    if recovery_warnings:
        lines.extend(
            [
                "",
                "| Trial | Initial rotation | Final rotation | Initial shift | Final shift | Reason |",
                "|:---|---:|---:|---:|---:|:---|",
            ]
        )
        for warning in recovery_warnings:
            lines.append(
                "| `{}` | {:.6g} | {:.6g} | {:.6g} | {:.6g} | {} |".format(
                    warning["trial_id"],
                    warning["initial_rotation_deg"],
                    warning["final_rotation_deg"],
                    warning["initial_shift_px"],
                    warning["final_shift_px"],
                    "; ".join(warning["reasons"]),
                )
            )
    else:
        lines.extend(["", "- None."])

    status_counts = Counter(str(row["status_name"]) for row in rows)
    lines.extend(["", "## Terminal outcomes", ""])
    for name, count in sorted(status_counts.items()):
        lines.append(f"- `{name}`: {count}")

    exact = next((row for row in rows if row["family"] == "exact"), None)
    if exact is not None:
        lines.extend(
            [
                "",
                "## Exact-pose control",
                "",
                f"- Status: `{exact['status_name']}`",
                f"- Accepted/attempted steps: {exact['accepted_steps']}/{exact['attempted_steps']}",
                f"- Final rotation error: {float(exact['final_rotation_deg']):.6g} degrees",
                f"- Final shift error: {float(exact['final_shift_px']):.6g} pixels",
            ]
        )

    groups = defaultdict(list)  # type: Dict[Tuple[str, int, int], List[Dict[str, Any]]]
    for row in rows:
        if row["family"] in ("rotation", "shift"):
            groups[(str(row["family"]), int(row["axis"]), int(row["sign"]))].append(row)

    lines.extend(
        [
            "",
            "## One-coordinate sweeps",
            "",
            "`Final/initial` is the error ratio for the injected parameter family. "
            "A value below 1 means that the final pose moved closer to truth. "
            "No automatic capture threshold is assigned.",
        ]
    )
    axis_names = {1: "x", 2: "y", 3: "z"}
    for key in sorted(groups):
        family, axis, sign = key
        group = sorted(groups[key], key=lambda row: float(row["injected_magnitude"]))
        unit = "degrees" if family == "rotation" else "pixels"
        lines.extend(
            [
                "",
                f"### {family} {axis_names[axis]}, {'positive' if sign > 0 else 'negative'}",
                "",
                f"| Injected ({unit}) | Initial error | Final error | Final/initial | "
                "Coupled final error | Objective reduction | Accepted/attempted | Switches | Status |",
                "|---:|---:|---:|---:|---:|---:|---:|---:|:---|",
            ]
        )
        for row in group:
            initial, final = target_values(row)
            ratio = final / initial if initial > 0.0 else math.nan
            before = float(row["objective_before"])
            after = float(row["objective_after"])
            reduction = (before - after) / before if before > 0.0 else 0.0
            coupled = (
                float(row["final_shift_px"])
                if family == "rotation"
                else float(row["final_rotation_deg"])
            )
            lines.append(
                f"| {float(row['injected_magnitude']):.4g} | {initial:.6g} | {final:.6g} | "
                f"{ratio:.6g} | {coupled:.6g} | {100.0 * reduction:.3f}% | "
                f"{row['accepted_steps']}/{row['attempted_steps']} | "
                f"{row['stencil_switches']} | `{row['status_name']}` |"
            )

    joint_rows = [row for row in rows if row["family"] == "joint"]
    if joint_rows:
        lines.extend(
            [
                "",
                "## Joint trials",
                "",
                "| Trial | Rotation vector (degrees) | Shift vector (pixels) | "
                "Final rotation error | Final shift error | Objective reduction | "
                "Accepted/attempted | Switches | Status |",
                "|:---|:---|:---|---:|---:|---:|---:|---:|:---|",
            ]
        )
        for row in joint_rows:
            before = float(row["objective_before"])
            after = float(row["objective_after"])
            reduction = (before - after) / before if before > 0.0 else 0.0
            rotation_vector, shift_vector = vector_text(row)
            lines.append(
                f"| `{row['trial_id']}` | {rotation_vector} | {shift_vector} | "
                f"{float(row['final_rotation_deg']):.6g} | "
                f"{float(row['final_shift_px']):.6g} | {100.0 * reduction:.3f}% | "
                f"{row['accepted_steps']}/{row['attempted_steps']} | "
                f"{row['stencil_switches']} | `{row['status_name']}` |"
            )

    if mechanism_rows:
        lines.extend(
            [
                "",
                "## Multi-volume mechanism diagnostic",
                "",
                "The table compares direct joint LM, the optional cumulative 15-degree/5-pixel "
                "guard, shift-first staging, and rotation-first staging. Rotation and shift "
                "errors are measured against truth. The guard is centered on the input seed.",
                "",
                "| Volume | Seed case | Route | Final rotation | Final shift | Objective reduction | "
                "Accepted/attempted | Within guard | Status |",
                "|:---|:---|:---|---:|---:|---:|---:|:---:|:---|",
            ]
        )
        for row in sorted(
            mechanism_rows,
            key=lambda item: (str(item["volume"]), str(item["case"]), str(item["route"])),
        ):
            before = float(row["objective_before"])
            after = float(row["objective_after"])
            reduction = (before - after) / before if before > 0.0 else 0.0
            lines.append(
                f"| {row['volume']} | {row['case']} | `{row['route']}` | "
                f"{float(row['final_rotation_norm_deg']):.6g} | "
                f"{float(row['final_shift_norm_px']):.6g} | {100.0 * reduction:.3f}% | "
                f"{row['accepted_steps']}/{row['attempted_steps']} | "
                f"{'yes' if int(row['within_15deg_5px_from_seed']) else 'no'} | "
                f"`{row['status_name']}` |"
            )

    if path_rows:
        path_groups = defaultdict(list)
        for row in path_rows:
            path_groups[(str(row["volume"]), str(row["case"]), str(row["target"]))].append(row)
        lines.extend(
            [
                "",
                "## Objective-path shape",
                "",
                "Each path uses a straight right-tangent rotation interpolation and a linear "
                "shift interpolation. `Uphill segments` counts sampled moves that increase the "
                "objective. A value above zero identifies a barrier on that straight path.",
                "",
                "| Volume | Seed case | Target | Peak/start | End/start | Uphill segments |",
                "|:---|:---|:---|---:|---:|---:|",
            ]
        )
        for key in sorted(path_groups):
            ordered = sorted(path_groups[key], key=lambda item: float(item["alpha"]))
            objectives = [float(item["objective"]) for item in ordered]
            start = objectives[0]
            scale = start if start > 0.0 else 1.0
            uphill = sum(
                right > left + 1.0e-12 * max(1.0, abs(left))
                for left, right in zip(objectives, objectives[1:])
            )
            lines.append(
                f"| {key[0]} | {key[1]} | `{key[2]}` | "
                f"{max(objectives) / scale:.6g} | {objectives[-1] / scale:.6g} | {uphill} |"
            )

    lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "- Review axes and signs separately. Do not hide a failed trial in an average.",
            "- Objective reduction alone is not pose recovery. Compare the final truth errors.",
            "- Stencil-switch counts identify trials that cross local interpolation cells.",
            "- The MRC files permit visual comparison of the observation, predictions, and residuals.",
            "- Discuss any asymmetric or non-monotonic capture behavior with Hans before integration work.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    root = args.result_root.resolve()
    artifacts = root / "artifacts"
    analysis = root / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    errors = []  # type: List[str]

    csv_path = artifacts / "pose_capture.csv"
    config_path = artifacts / "configuration.txt"
    for required in (csv_path, config_path, artifacts / "truth_volume.mrc", artifacts / "truth_observation.mrc"):
        if not required.is_file() or required.stat().st_size == 0:
            errors.append(f"missing required artifact: {required.relative_to(root)}")
    configuration = read_configuration(config_path) if config_path.is_file() else {}
    scenario_version = configuration.get("scenario_version", "1")
    rows = read_rows(csv_path, errors, scenario_version) if csv_path.is_file() else []
    if rows and configuration:
        errors.extend(validate_rows(rows, configuration, artifacts))
    recovery_warnings = scientific_recovery_warnings(rows)
    mechanism_rows = []  # type: List[Dict[str, Any]]
    trajectory_rows = []  # type: List[Dict[str, Any]]
    path_rows = []  # type: List[Dict[str, Any]]
    if scenario_version == "4":
        mechanism_rows, trajectory_rows, path_rows = read_mechanism_evidence(artifacts, errors)
        errors.extend(
            validate_mechanism_evidence(mechanism_rows, trajectory_rows, path_rows, artifacts)
        )

    (analysis / "summary.md").write_text(
        markdown_summary(
            rows, configuration, errors, recovery_warnings, mechanism_rows, path_rows
        ),
        encoding="utf-8",
    )
    (analysis / "summary.json").write_text(
        json.dumps(
            {
                "integrity": "FAIL" if errors else "PASS",
                "errors": errors,
                "scientific_recovery_warning_count": len(recovery_warnings),
                "scientific_recovery_warnings": recovery_warnings,
                "mechanism_summary_rows": mechanism_rows,
                "mechanism_trajectory_row_count": len(trajectory_rows),
                "objective_path_row_count": len(path_rows),
                "configuration": configuration,
                "rows": rows,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"POSE_CAPTURE_ANALYSIS: {'FAIL' if errors else 'PASS'}")
    print(f"Review: {analysis / 'summary.md'}")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
