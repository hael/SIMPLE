#!/usr/bin/env python3
"""Analyze clean/noisy exact/perturbed PCG pose-polishing diagnostics."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

from analyze_pose_ab import (
    euler2m,
    fsc_area,
    parse_fsc,
    parse_polish_summary,
    parse_rows,
    point_group_distance,
    rms,
    symmetry_matrices,
)


POSE_KEYS = 6
EXACT_ROTATION_TOLERANCE_RAD = 0.001745
EXACT_SHIFT_TOLERANCE_PX = 0.05
FSC_AREA_TOLERANCE = 0.005
CLEAN_RECOVERY_FRACTION = 0.50
NOISY_RECOVERY_FRACTION = 0.10
KEY_PATTERN = re.compile(r"(?<![A-Za-z0-9_])([A-Za-z0-9_]+)=([-+0-9.Ee]+)")


def parse_truth(path: Path) -> dict[int, tuple[int, tuple[float, ...]]]:
    """Read the simulator truth poses in the common project-row representation."""
    rows: dict[int, tuple[int, tuple[float, ...]]] = {}
    for index, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        if not raw.strip():
            continue
        values = {key: float(value) for key, value in KEY_PATTERN.findall(raw)}
        missing = [key for key in ("e1", "e2", "e3", "x", "y") if key not in values]
        if missing:
            raise ValueError(f"truth row {index} is missing {','.join(missing)}")
        rows[len(rows) + 1] = (
            1,
            (values["e1"], values["e2"], values["e3"], values["x"], values["y"], 0.0),
        )
    if not rows:
        raise ValueError(f"no truth rows were parsed from {path}")
    return rows


def pose_errors(
    poses: dict[int, tuple[int, tuple[float, ...]]],
    truth: dict[int, tuple[int, tuple[float, ...]]],
    symmetries: list[list[list[float]]],
) -> tuple[list[float], list[float]]:
    """Calculate per-particle symmetry-aware rotation and image-shift errors."""
    if poses.keys() != truth.keys():
        raise ValueError("project and truth orientation tables have different particle index sets")
    rotations: list[float] = []
    shifts: list[float] = []
    for index in poses:
        state, values = poses[index]
        _, expected = truth[index]
        if state <= 0:
            continue
        rotations.append(
            point_group_distance(euler2m(values[0:3]), euler2m(expected[0:3]), symmetries)
        )
        shifts.append(math.hypot(values[3] - expected[3], values[4] - expected[4]))
    if not rotations:
        raise ValueError("no active project poses were available for truth comparison")
    return rotations, shifts


def load_arm(root: Path, name: str) -> dict[str, object]:
    """Load poses, FSC curves, and optional polish telemetry for one arm."""
    arm = root / "arms" / name
    evidence = arm / "evidence"
    return {
        "poses": parse_rows(evidence / "poses_after.txt", POSE_KEYS),
        "half_fsc": parse_fsc(evidence / "half_fsc.txt"),
        "truth_avg_fsc": parse_fsc(evidence / "truth_fsc_avg.txt"),
        "truth_even_fsc": parse_fsc(evidence / "truth_fsc_even.txt"),
        "truth_odd_fsc": parse_fsc(evidence / "truth_fsc_odd.txt"),
        "polish": parse_polish_summary(arm / "reconstruct3D.log") if name.endswith("_on") else None,
    }


def area_changes(off: dict[str, object], on: dict[str, object]) -> dict[str, float]:
    """Calculate enabled-minus-disabled FSC-area changes for all map pairs."""
    changes: dict[str, float] = {}
    for key in ("half_fsc", "truth_avg_fsc", "truth_even_fsc", "truth_odd_fsc"):
        changes[key] = fsc_area(on[key]) - fsc_area(off[key])
    return changes


def analyze_case(
    root: Path,
    case: str,
    start: str,
    truth: dict[int, tuple[int, tuple[float, ...]]],
    symmetries: list[list[list[float]]],
) -> dict[str, object]:
    """Apply the predeclared pose and FSC gates to one matched off/on case."""
    off = load_arm(root, f"{case}_off")
    on = load_arm(root, f"{case}_on")
    off_rotation, off_shift = pose_errors(off["poses"], truth, symmetries)
    on_rotation, on_shift = pose_errors(on["poses"], truth, symmetries)
    rotation_before = rms(off_rotation)
    rotation_after = rms(on_rotation)
    shift_before = rms(off_shift)
    shift_after = rms(on_shift)
    changes = area_changes(off, on)
    fsc_stable = min(changes.values()) >= -FSC_AREA_TOLERANCE
    reasons: list[str] = []
    if start == "exact":
        pose_pass = (
            rotation_after <= rotation_before + EXACT_ROTATION_TOLERANCE_RAD
            and shift_after <= shift_before + EXACT_SHIFT_TOLERANCE_PX
        )
        if not pose_pass:
            reasons.append("exact truth poses moved beyond the stability tolerance")
    else:
        required = CLEAN_RECOVERY_FRACTION if case.startswith("clean_") else NOISY_RECOVERY_FRACTION
        pose_pass = (
            rotation_after <= (1.0 - required) * rotation_before
            and shift_after <= (1.0 - required) * shift_before
        )
        if not pose_pass:
            reasons.append(f"both pose RMS errors did not decrease by at least {100.0 * required:.0f}%")
    if not fsc_stable:
        reasons.append("an FSC area declined beyond the predeclared tolerance")
    polish = on["polish"]
    if polish["invalid"] != 0:
        reasons.append("the enabled arm reported invalid numerical outcomes")
    return {
        "conclusion": "PASS" if pose_pass and fsc_stable and not reasons else "FAIL",
        "rotation_rms_before_rad": rotation_before,
        "rotation_rms_after_rad": rotation_after,
        "shift_rms_before_px": shift_before,
        "shift_rms_after_px": shift_after,
        "fsc_area_changes": changes,
        "polish_summary": polish,
        "reasons": reasons,
    }


def write_markdown(path: Path, report: dict[str, object]) -> None:
    """Write a concise human-readable report from the complete result object."""
    lines = [
        "# Continuous 3-D pose truth-diagnostic result",
        "",
        f"**Feature gate:** {report['conclusion']}",
        "",
        "This matrix separates clean operator consistency, clean local recovery,",
        "noisy exact-pose stability, and noisy recovery. A lower internal objective",
        "does not override a truth-pose or FSC failure.",
        "",
        "## Case results",
        "",
    ]
    for name, result in report["cases"].items():
        lines.extend(
            [
                f"### {name}: {result['conclusion']}",
                "",
                f"- Rotation RMS: {result['rotation_rms_before_rad']:.6g} -> "
                f"{result['rotation_rms_after_rad']:.6g} rad",
                f"- Shift RMS: {result['shift_rms_before_px']:.6g} -> "
                f"{result['shift_rms_after_px']:.6g} pixels",
                f"- Half-map FSC-area change: {result['fsc_area_changes']['half_fsc']:.6g}",
                f"- Truth-average FSC-area change: "
                f"{result['fsc_area_changes']['truth_avg_fsc']:.6g}",
            ]
        )
        if result["reasons"]:
            lines.extend(f"- Failure: {reason}" for reason in result["reasons"])
        lines.append("")
    lines.extend(["## Interpretation", ""])
    lines.extend(f"- {item}" for item in report["interpretation"])
    lines.extend(["", "## Feature-gate rule", ""])
    lines.append(
        "Clean exact and clean perturbed must pass. At least one noisy frequency "
        "band must pass both its exact and perturbed cases."
    )
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    """Analyze all eight cases and apply the aggregate scientific feature gate."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--truth-oris", type=Path, required=True)
    parser.add_argument("--pgrp", required=True)
    args = parser.parse_args()
    cases = {
        "clean_exact_full": "exact",
        "clean_perturbed_full": "perturbed",
        "noisy_exact_full": "exact",
        "noisy_perturbed_full": "perturbed",
        "noisy_exact_fsc05": "exact",
        "noisy_perturbed_fsc05": "perturbed",
        "noisy_exact_fsc0143": "exact",
        "noisy_perturbed_fsc0143": "perturbed",
    }
    errors: list[str] = []
    results: dict[str, object] = {}
    try:
        truth = parse_truth(args.truth_oris)
        symmetries = symmetry_matrices(args.pgrp)
        for name, start in cases.items():
            try:
                results[name] = analyze_case(args.root, name, start, truth, symmetries)
            except (OSError, ValueError, KeyError) as error:
                errors.append(f"{name}: {error}")
    except (OSError, ValueError) as error:
        errors.append(str(error))
    interpretation: list[str] = []
    if results:
        if results.get("clean_exact_full", {}).get("conclusion") != "PASS":
            interpretation.append(
                "Clean exact failure supports a mismatch between the truth generator and the "
                "fixed-map objective, including operator inconsistency or reconstruction-map bias."
            )
        elif results.get("clean_perturbed_full", {}).get("conclusion") != "PASS":
            interpretation.append("Clean exact is stable, but clean recovery fails; inspect pose numerics and solver policy.")
        else:
            interpretation.append("Clean exact stability and clean local recovery pass.")
        if results.get("noisy_exact_full", {}).get("conclusion") != "PASS":
            interpretation.append("Full-band noisy exact failure supports noise overfitting.")
        if any(
            results.get(f"noisy_exact_{band}", {}).get("conclusion") == "PASS"
            and results.get(f"noisy_perturbed_{band}", {}).get("conclusion") == "PASS"
            for band in ("fsc05", "fsc0143")
        ):
            interpretation.append("A data-supported frequency cutoff provides stable and useful noisy recovery.")
    clean_gate = all(results.get(name, {}).get("conclusion") == "PASS" for name in (
        "clean_exact_full", "clean_perturbed_full"
    ))
    noisy_gate = any(
        results.get(f"noisy_exact_{band}", {}).get("conclusion") == "PASS"
        and results.get(f"noisy_perturbed_{band}", {}).get("conclusion") == "PASS"
        for band in ("full", "fsc05", "fsc0143")
    )
    conclusion = "PASS" if clean_gate and noisy_gate and not errors else "FAIL"
    report = {
        "conclusion": conclusion,
        "tolerances": {
            "exact_rotation_rad": EXACT_ROTATION_TOLERANCE_RAD,
            "exact_shift_px": EXACT_SHIFT_TOLERANCE_PX,
            "fsc_area_decline": FSC_AREA_TOLERANCE,
            "clean_recovery_fraction": CLEAN_RECOVERY_FRACTION,
            "noisy_recovery_fraction": NOISY_RECOVERY_FRACTION,
        },
        "cases": results,
        "interpretation": interpretation,
        "errors": errors,
    }
    analysis = args.root / "analysis"
    analysis.mkdir(exist_ok=True)
    (analysis / "truth_matrix.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    write_markdown(analysis / "truth_matrix.md", report)
    print(f"CONTINUOUS_3D_POSE_TRUTH_DIAGNOSTIC: {conclusion}")
    for error in errors:
        print(f"ERROR: {error}")
    return 0 if conclusion == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
