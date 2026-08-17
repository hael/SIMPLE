#!/usr/bin/env python3
"""Create a deterministic joint-pose perturbation of a SIMPLE orientation table."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


ROTATION_E1_DEGREES = 0.8
ROTATION_E3_DEGREES = 0.6
SHIFT_X_PIXELS = 0.25
SHIFT_Y_PIXELS = 0.20


def replace_value(line: str, key: str, value: float) -> str:
    """Replace exactly one numeric key in a SIMPLE orientation-table row."""
    pattern = re.compile(rf"(?<![A-Za-z0-9_]){re.escape(key)}=([-+0-9.Ee]+)")
    updated, count = pattern.subn(f"{key}={value:.9g}", line, count=1)
    if count != 1:
        raise ValueError(f"orientation row does not contain exactly one {key} value")
    return updated


def read_value(line: str, key: str) -> float:
    """Read one required numeric key from a SIMPLE orientation-table row."""
    match = re.search(rf"(?<![A-Za-z0-9_]){re.escape(key)}=([-+0-9.Ee]+)", line)
    if match is None:
        raise ValueError(f"orientation row is missing {key}")
    return float(match.group(1))


def perturb_line(line: str, index: int) -> str:
    """Apply the deterministic alternating joint-pose perturbation to one row."""
    sign = 1.0 if index % 2 else -1.0
    updates = {
        "e1": (read_value(line, "e1") + sign * ROTATION_E1_DEGREES) % 360.0,
        "e3": (read_value(line, "e3") - sign * ROTATION_E3_DEGREES) % 360.0,
        "x": read_value(line, "x") + sign * SHIFT_X_PIXELS,
        "y": read_value(line, "y") - sign * SHIFT_Y_PIXELS,
    }
    updated = line
    for key, value in updates.items():
        updated = replace_value(updated, key, value)
    return updated


def normalize_even_odd(lines: list[str]) -> tuple[list[str], bool]:
    """Keep a valid split or create deterministic alternating half ownership."""
    values = {int(round(read_value(line, "eo"))) for line in lines}
    if values == {0, 1}:
        return lines, False
    normalized = [
        replace_value(line, "eo", 0.0 if index % 2 else 1.0)
        for index, line in enumerate(lines, start=1)
    ]
    return normalized, True


def validate_constant(lines: list[str], key: str, expected: float) -> None:
    """Require one microscope value to be uniform and equal to the run input."""
    values = [read_value(line, key) for line in lines]
    if any(abs(value - expected) > 1.0e-5 * max(abs(expected), 1.0) for value in values):
        raise ValueError(f"orientation table {key} values do not match {expected}")


def main() -> int:
    """Validate truth metadata and write matched exact and perturbed tables."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--exact-output", type=Path)
    parser.add_argument("--kv", type=float)
    parser.add_argument("--cs", type=float)
    parser.add_argument("--fraca", type=float)
    args = parser.parse_args()
    lines = [line for line in args.input.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise ValueError("truth orientation table is empty")
    for key in ("kv", "cs", "fraca"):
        expected = getattr(args, key)
        if expected is not None:
            validate_constant(lines, key, expected)
    exact, assigned_even_odd = normalize_even_odd(lines)
    perturbed = [perturb_line(line, index) for index, line in enumerate(exact, start=1)]
    if args.exact_output is not None:
        args.exact_output.write_text("\n".join(exact) + "\n", encoding="utf-8")
    args.output.write_text("\n".join(perturbed) + "\n", encoding="utf-8")
    print(f"WROTE_PERTURBED_ORITAB rows={len(lines)} output={args.output}")
    if args.exact_output is not None:
        print(f"WROTE_EXACT_ORITAB rows={len(lines)} output={args.exact_output}")
    print(f"EVEN_ODD_ASSIGNED value={'yes' if assigned_even_odd else 'no'}")
    print(
        "PERTURBATION "
        f"e1=+/-{ROTATION_E1_DEGREES}deg "
        f"e3=-/+{ROTATION_E3_DEGREES}deg "
        f"x=+/-{SHIFT_X_PIXELS}px y=-/+{SHIFT_Y_PIXELS}px"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
