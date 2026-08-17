#!/usr/bin/env python3
"""Analyze the matched continuous 3-D PCG pose-polishing validation arms."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path


POSE_KEYS = 6
INVARIANT_KEYS = 4
FSC_AREA_TOLERANCE = 0.01
ROUTE_ROTATION_RMS_TOLERANCE = 0.002
ROUTE_SHIFT_RMS_TOLERANCE = 0.02


def parse_rows(path: Path, payload_count: int) -> dict[int, tuple[int, tuple[float, ...]]]:
    rows: dict[int, tuple[int, tuple[float, ...]]] = {}
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = raw.split()
        if len(fields) < 2 + payload_count:
            continue
        try:
            index = int(fields[0])
            state = int(fields[1])
            payload = tuple(float(value) for value in fields[2 : 2 + payload_count])
        except ValueError:
            continue
        rows[index] = (state, payload)
    if not rows:
        raise ValueError(f"no project rows were parsed from {path}")
    return rows


def matmul(left: list[list[float]], right: list[list[float]]) -> list[list[float]]:
    return [
        [sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3)]
        for i in range(3)
    ]


def transpose(matrix: list[list[float]]) -> list[list[float]]:
    return [[matrix[j][i] for j in range(3)] for i in range(3)]


def simple_rotmat(angle_degrees: float, axis: int) -> list[list[float]]:
    angle = math.radians(angle_degrees)
    cosine = math.cos(angle)
    sine = math.sin(angle)
    if axis == 2:
        return [[cosine, 0.0, -sine], [0.0, 1.0, 0.0], [sine, 0.0, cosine]]
    if axis == 3:
        return [[cosine, sine, 0.0], [-sine, cosine, 0.0], [0.0, 0.0, 1.0]]
    raise ValueError(f"unsupported SIMPLE rotation axis: {axis}")


def euler2m(eulers: tuple[float, float, float]) -> list[list[float]]:
    first = simple_rotmat(eulers[0], 3)
    tilt = simple_rotmat(eulers[1], 2)
    third = simple_rotmat(eulers[2], 3)
    return matmul(matmul(third, tilt), first)


def symmetry_matrices(point_group: str) -> list[list[list[float]]]:
    identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    normalized = point_group.lower()
    if normalized == "c1":
        return [identity]
    if normalized == "d2":
        return [
            identity,
            [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
            [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
            [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
        ]
    raise ValueError("the validation checker currently supports pgrp=c1 or pgrp=d2")


def rotation_distance(left: list[list[float]], right: list[list[float]]) -> float:
    relative = matmul(transpose(left), right)
    cosine = 0.5 * (relative[0][0] + relative[1][1] + relative[2][2] - 1.0)
    return math.acos(max(-1.0, min(1.0, cosine)))


def point_group_distance(
    left: list[list[float]], right: list[list[float]], symmetries: list[list[list[float]]]
) -> float:
    return min(rotation_distance(left, matmul(right, symmetry)) for symmetry in symmetries)


def parse_polish_summary(log_path: Path) -> dict[str, int]:
    text = log_path.read_text(encoding="utf-8", errors="replace")
    patterns = {
        "primary": re.compile(
            r"PCG POSE POLISH: PARTICLES\s+(\d+)\s+IMPROVED\s+(\d+)\s+"
            r"UNCHANGED\s+(\d+)\s+UNRELIABLE\s+(\d+)"
        ),
        "terminal": re.compile(
            r"PCG POSE POLISH: STEP-BOUND\s+(\d+)\s+INVALID\s+(\d+)\s+"
            r"ITERATION-LIMIT\s+(\d+)"
        ),
        "steps": re.compile(
            r"PCG POSE POLISH: ATTEMPTED STEPS\s+(\d+)\s+ACCEPTED STEPS\s+(\d+)\s+"
            r"STENCIL SWITCHES\s+(\d+)"
        ),
    }
    matches = {name: pattern.findall(text) for name, pattern in patterns.items()}
    if any(len(items) != 1 for items in matches.values()):
        raise ValueError(f"{log_path} does not contain exactly one complete polish summary")
    primary = tuple(int(value) for value in matches["primary"][0])
    terminal = tuple(int(value) for value in matches["terminal"][0])
    steps = tuple(int(value) for value in matches["steps"][0])
    return {
        "particles": primary[0],
        "improved": primary[1],
        "unchanged": primary[2],
        "unreliable": primary[3],
        "step_bound": terminal[0],
        "invalid": terminal[1],
        "iteration_limit": terminal[2],
        "attempted_steps": steps[0],
        "accepted_steps": steps[1],
        "stencil_switches": steps[2],
        "final_reconstruction_markers": text.count(
            "PCG POSE POLISH: RUNNING FINAL PCG RECONSTRUCTION"
        ),
    }


def parse_fsc(path: Path) -> list[float]:
    pattern = re.compile(r">>> FSC:\s+([-+0-9.Ee]+)")
    values = [float(match.group(1)) for match in pattern.finditer(path.read_text(errors="replace"))]
    if len(values) < 3:
        raise ValueError(f"too few FSC shells were parsed from {path}")
    return values


def fsc_area(values: list[float]) -> float:
    shells = values[1:]
    return sum(max(-1.0, min(1.0, value)) for value in shells) / len(shells)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_arm(root: Path, name: str) -> dict[str, object]:
    evidence = root / "arms" / name / "evidence"
    arm_dir = root / "arms" / name
    volumes = [arm_dir / "recvol_state01_even.mrc", arm_dir / "recvol_state01_odd.mrc"]
    if not all(path.is_file() for path in volumes):
        raise ValueError(f"{name} is missing a final even or odd volume")
    return {
        "pose": parse_rows(evidence / "poses_after.txt", POSE_KEYS),
        "invariant": parse_rows(evidence / "invariants_after.txt", INVARIANT_KEYS),
        "fsc": parse_fsc(evidence / "fsc.txt"),
        "log": arm_dir / "reconstruct3D.log",
        "volume_hashes": [sha256_file(path) for path in volumes],
    }


def compare_metadata(
    left: dict[int, tuple[int, tuple[float, ...]]],
    right: dict[int, tuple[int, tuple[float, ...]]],
    label: str,
    failures: list[str],
) -> None:
    if left.keys() != right.keys():
        failures.append(f"{label}: particle index sets differ")
        return
    for index in left:
        if left[index] != right[index]:
            failures.append(f"{label}: immutable metadata differs at particle {index}")
            return


def pose_deltas(
    left: dict[int, tuple[int, tuple[float, ...]]],
    right: dict[int, tuple[int, tuple[float, ...]]],
    symmetries: list[list[list[float]]],
) -> tuple[list[float], list[float]]:
    if left.keys() != right.keys():
        raise ValueError("pose exports contain different particle index sets")
    rotations: list[float] = []
    shifts: list[float] = []
    for index in left:
        left_state, left_values = left[index]
        right_state, right_values = right[index]
        if left_state != right_state:
            raise ValueError(f"particle {index} changed state")
        left_matrix = euler2m(left_values[0:3])
        right_matrix = euler2m(right_values[0:3])
        rotations.append(point_group_distance(left_matrix, right_matrix, symmetries))
        shifts.append(math.hypot(right_values[3] - left_values[3], right_values[4] - left_values[4]))
    return rotations, shifts


def rms(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values) / len(values))


def analyze_pair(
    root: Path,
    route: str,
    arms: dict[str, dict[str, object]],
    symmetries: list[list[list[float]]],
    angular_bound: float,
    failures: list[str],
) -> dict[str, object]:
    off = arms[f"{route}_off"]
    on = arms[f"{route}_on"]
    compare_metadata(off["invariant"], on["invariant"], f"{route} A/B", failures)
    rotations, shifts = pose_deltas(off["pose"], on["pose"], symmetries)
    summary = parse_polish_summary(on["log"])
    off_text = Path(off["log"]).read_text(encoding="utf-8", errors="replace")
    accounted = sum(
        summary[key]
        for key in ("improved", "unchanged", "unreliable", "step_bound", "invalid", "iteration_limit")
    )
    changed = sum(
        rotation > 2.0e-4 or shift > 2.0e-4 for rotation, shift in zip(rotations, shifts)
    )
    active_particles = sum(state > 0 for state, _ in off["pose"].values())
    off_markers = off_text.count("PCG POSE POLISH:")
    off_area = fsc_area(off["fsc"])
    on_area = fsc_area(on["fsc"])
    if off_markers != 0:
        failures.append(f"{route}: disabled arm contains pose-polish markers")
    if summary["final_reconstruction_markers"] != 1:
        failures.append(f"{route}: enabled arm did not run exactly one final PCG reconstruction")
    if accounted != summary["particles"]:
        failures.append(f"{route}: terminal outcome counts do not balance")
    if summary["particles"] != active_particles:
        failures.append(f"{route}: polished-particle count does not match active project particles")
    if summary["improved"] < 1 or changed < 1:
        failures.append(f"{route}: no measurable accepted pose change occurred")
    if changed > summary["improved"]:
        failures.append(f"{route}: more poses changed than were reported improved")
    if summary["invalid"] != 0:
        failures.append(f"{route}: invalid numerical outcomes occurred")
    if summary["accepted_steps"] < summary["improved"]:
        failures.append(f"{route}: accepted-step count is smaller than improved-particle count")
    if max(rotations) > angular_bound + 2.0e-3:
        failures.append(f"{route}: a cumulative angular change exceeded the declared LM bound")
    if max(shifts) > 8.001:
        failures.append(f"{route}: a cumulative shift change exceeded eight one-pixel LM steps")
    if on_area < off_area - FSC_AREA_TOLERANCE:
        failures.append(f"{route}: enabled FSC area declined by more than {FSC_AREA_TOLERANCE}")
    changed_volumes = sum(
        left != right for left, right in zip(off["volume_hashes"], on["volume_hashes"])
    )
    if summary["improved"] > 0 and changed_volumes < 1:
        failures.append(f"{route}: accepted poses did not change either final half-map artifact")
    return {
        "summary": summary,
        "changed_particles_at_export_precision": changed,
        "rotation_max_rad": max(rotations),
        "rotation_mean_rad": sum(rotations) / len(rotations),
        "shift_max_px": max(shifts),
        "shift_mean_px": sum(shifts) / len(shifts),
        "fsc_area_off": off_area,
        "fsc_area_on": on_area,
        "fsc_area_change": on_area - off_area,
        "changed_final_halfmaps": changed_volumes,
    }


def analyze_default_behavior(
    arms: dict[str, dict[str, object]], failures: list[str]
) -> dict[str, object]:
    default = arms["shared_default"]
    explicit_off = arms["shared_off"]
    default_log = Path(default["log"]).read_text(encoding="utf-8", errors="replace")
    if default_log.count("PCG POSE POLISH:") != 0:
        failures.append("default arm contains pose-polish markers")
    if default["pose"] != explicit_off["pose"]:
        failures.append("omitted pcg_pose_polish does not match explicit no poses")
    if default["invariant"] != explicit_off["invariant"]:
        failures.append("omitted pcg_pose_polish does not match explicit no metadata")
    if default["fsc"] != explicit_off["fsc"]:
        failures.append("omitted pcg_pose_polish does not match explicit no FSC")
    return {
        "poses_identical": default["pose"] == explicit_off["pose"],
        "metadata_identical": default["invariant"] == explicit_off["invariant"],
        "fsc_identical": default["fsc"] == explicit_off["fsc"],
    }


def analyze_route_equivalence(
    arms: dict[str, dict[str, object]],
    symmetries: list[list[list[float]]],
    failures: list[str],
) -> dict[str, object]:
    result: dict[str, object] = {}
    for mode in ("off", "on"):
        shared = arms[f"shared_{mode}"]
        distributed = arms[f"distributed_{mode}"]
        compare_metadata(
            shared["invariant"], distributed["invariant"], f"{mode} shared/distributed", failures
        )
        rotations, shifts = pose_deltas(shared["pose"], distributed["pose"], symmetries)
        shared_area = fsc_area(shared["fsc"])
        distributed_area = fsc_area(distributed["fsc"])
        rotation_rms = rms(rotations)
        shift_rms = rms(shifts)
        if rotation_rms > ROUTE_ROTATION_RMS_TOLERANCE:
            failures.append(f"{mode}: shared/distributed rotation RMS exceeds tolerance")
        if shift_rms > ROUTE_SHIFT_RMS_TOLERANCE:
            failures.append(f"{mode}: shared/distributed shift RMS exceeds tolerance")
        if abs(shared_area - distributed_area) > FSC_AREA_TOLERANCE:
            failures.append(f"{mode}: shared/distributed FSC areas differ beyond tolerance")
        result[mode] = {
            "rotation_rms_rad": rotation_rms,
            "shift_rms_px": shift_rms,
            "shared_fsc_area": shared_area,
            "distributed_fsc_area": distributed_area,
        }
    return result


def write_markdown(path: Path, report: dict[str, object]) -> None:
    lines = [
        "# Continuous 3-D PCG pose validation result",
        "",
        f"**Conclusion:** {report['conclusion']}",
        "",
        "## A/B results",
        "",
    ]
    for route in ("shared", "distributed"):
        pair = report["pairs"].get(route)
        if pair is None:
            continue
        lines.extend(
            [
                f"### {route}",
                "",
                f"- Improved particles: {pair['summary']['improved']} / {pair['summary']['particles']}",
                f"- Maximum angular change: {pair['rotation_max_rad']:.6g} rad",
                f"- Maximum shift change: {pair['shift_max_px']:.6g} pixels",
                f"- FSC area change: {pair['fsc_area_change']:.6g}",
                "",
            ]
        )
    lines.extend(["## Failures", ""])
    if report["failures"]:
        lines.extend(f"- {failure}" for failure in report["failures"])
    else:
        lines.append("- None")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--pgrp", required=True)
    parser.add_argument("--angular-bound", type=float, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    failures: list[str] = []
    try:
        symmetries = symmetry_matrices(args.pgrp)
        arms = {
            name: load_arm(root, name)
            for name in (
                "shared_default",
                "shared_off",
                "shared_on",
                "distributed_off",
                "distributed_on",
            )
        }
        default_behavior = analyze_default_behavior(arms, failures)
        pairs = {
            route: analyze_pair(root, route, arms, symmetries, args.angular_bound, failures)
            for route in ("shared", "distributed")
        }
        route_equivalence = analyze_route_equivalence(arms, symmetries, failures)
    except (OSError, ValueError, KeyError) as error:
        failures.append(str(error))
        pairs = {}
        default_behavior = {}
        route_equivalence = {}
    report = {
        "conclusion": "PASS" if not failures else "FAIL",
        "fsc_area_tolerance": FSC_AREA_TOLERANCE,
        "route_rotation_rms_tolerance_rad": ROUTE_ROTATION_RMS_TOLERANCE,
        "route_shift_rms_tolerance_px": ROUTE_SHIFT_RMS_TOLERANCE,
        "cumulative_angular_bound_rad": args.angular_bound,
        "pairs": pairs,
        "default_behavior": default_behavior,
        "route_equivalence": route_equivalence,
        "failures": failures,
    }
    analysis_dir = root / "analysis"
    analysis_dir.mkdir(exist_ok=True)
    (analysis_dir / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    write_markdown(analysis_dir / "summary.md", report)
    print(f"CONTINUOUS_3D_POSE_VALIDATION: {report['conclusion']}")
    for failure in failures:
        print(f"FAIL: {failure}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
