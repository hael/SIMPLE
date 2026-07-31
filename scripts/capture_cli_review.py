#!/usr/bin/env python3
"""Write a fast, reviewable snapshot of the generated SIMPLE UI registry.

The registry is emitted once by ``simple_private_exec prg=print_ui_json``.
This script saves that exact JSON stdout and also writes a readable text
inventory of the same program, input, visibility, default, option, activation,
group, and requirement metadata. With ``--capture-cli``, it also writes one
plain-text CLI usage review for every registered program. It deliberately does
not run CTest or workflow programs.

Example:

  python3 scripts/capture_cli_review.py --build-dir build
  python3 scripts/capture_cli_review.py --build-dir build \\
      --capture-cli
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


INPUT_SECTIONS = (
    ("image input/output", "IMAGE INPUT/OUTPUT"),
    ("file input/output", "FILE INPUT/OUTPUT"),
    ("parameter input/output", "PARAMETER INPUT/OUTPUT"),
    ("search controls", "SEARCH CONTROLS"),
    ("filter controls", "FILTER CONTROLS"),
    ("mask controls", "MASK CONTROLS"),
    ("computer controls", "COMPUTER CONTROLS"),
)
ANSI_ESCAPE = re.compile(rb"\x1b\[[0-?]*[ -/]*[@-~]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build-dir",
        type=Path,
        required=True,
        help="CMake build directory containing production/simple_private_exec.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Readable review text file (default: <build-dir>/simple_ui_review.txt).",
    )
    parser.add_argument(
        "--registry-output",
        type=Path,
        help="Exact JSON stdout file (default: alongside --output with a .json suffix).",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Seconds allowed for the single registry query (default: 30).",
    )
    parser.add_argument(
        "--capture-cli",
        action="store_true",
        help=(
            "Capture exact usage output from every registered program. This invokes the "
            "parser-level usage=yes mode, which stops before workflow execution."
        ),
    )
    parser.add_argument(
        "--cli-output",
        type=Path,
        help="Combined plain-text --capture-cli review file (default: <build-dir>/simple_cli_instructions.txt).",
    )
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be greater than 0")
    args.build_dir = args.build_dir.resolve()
    if args.output is None:
        args.output = args.build_dir / "simple_ui_review.txt"
    else:
        args.output = args.output.resolve()
    if args.registry_output is None:
        args.registry_output = args.output.with_suffix(".json")
    else:
        args.registry_output = args.registry_output.resolve()
    if args.cli_output is None:
        args.cli_output = args.build_dir / "simple_cli_instructions.txt"
    else:
        args.cli_output = args.cli_output.resolve()
    return args


def private_executable(build_dir: Path) -> Path:
    suffix = ".exe" if sys.platform.startswith("win") else ""
    return build_dir / "production" / f"simple_private_exec{suffix}"


def production_executable(build_dir: Path, executable_name: str) -> Path:
    suffix = ".exe" if sys.platform.startswith("win") else ""
    if executable_name == "all":
        executable_name = "simple_exec"
    return build_dir / "production" / f"{executable_name}{suffix}"


def load_registry(executable: Path, timeout: float) -> tuple[str, dict[str, Any]]:
    command = (str(executable), "prg=print_ui_json")
    try:
        completed = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(f"registry query timed out after {timeout:g} seconds") from error
    if completed.returncode != 0:
        detail = (completed.stdout + completed.stderr).strip()
        raise RuntimeError(f"registry query failed with exit code {completed.returncode}: {detail}")
    try:
        registry = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise RuntimeError(f"registry query emitted invalid JSON: {error}") from error
    if not isinstance(registry, dict):
        raise RuntimeError("registry query emitted a non-object JSON value")
    return completed.stdout, registry


def capture_cli_instructions(
    build_dir: Path, name: str, program: dict[str, Any], timeout: float
) -> bytes:
    executable_name = program["program"]["executable"]
    executable = production_executable(build_dir, executable_name)
    if not executable.is_file():
        raise RuntimeError(f"CLI executable not found for {name}: {executable}")
    command = (str(executable), f"prg={name}", "usage=yes")
    try:
        completed = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(f"CLI capture for {name} timed out after {timeout:g} seconds") from error
    if completed.returncode != 0:
        raise RuntimeError(f"CLI capture for {name} exited with {completed.returncode}")
    return completed.stdout


def sorted_programs(registry: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    return sorted(
        registry.items(),
        key=lambda item: (
            item[1]["program"]["executable"],
            item[1]["program"]["category_order"],
            item[1]["program"]["category"],
            item[0],
        ),
    )


def render_cli_instructions(captures: list[tuple[str, dict[str, Any], bytes]]) -> str:
    lines = [
        "SIMPLE CLI instructions",
        f"Generated: {dt.datetime.now().astimezone().isoformat(timespec='seconds')}",
        "Each section is the program's usage output with terminal ANSI styling removed.",
        "",
    ]
    for name, program, output in captures:
        executable_name = program["program"]["executable"]
        lines.extend(("=" * 96, f"{executable_name}:{name}", "-" * 96))
        text = ANSI_ESCAPE.sub(b"", output).decode("utf-8", errors="replace")
        lines.append(text.rstrip("\n"))
        lines.append("")
    return "\n".join(lines)


def value_text(value: Any) -> str:
    if isinstance(value, str):
        return value
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def input_by_key(program: dict[str, Any]) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for section, _ in INPUT_SECTIONS:
        for input_value in program.get(section, []):
            result[input_value["key"]] = input_value
    return result


def write_input(lines: list[str], input_value: dict[str, Any]) -> None:
    qualifiers = [input_value["keytype"], "required" if input_value["required"] else "optional"]
    qualifiers.append(f"visibility={input_value['visibility']}")
    if input_value.get("has_default"):
        qualifiers.append(f"default={value_text(input_value.get('default', ''))}")
    lines.append(f"  {input_value['key']} [{'; '.join(qualifiers)}]")
    lines.append(f"    Label: {input_value['label']}")
    lines.append(f"    Help: {input_value['help']}")
    if input_value.get("placeholder"):
        lines.append(f"    Placeholder: {input_value['placeholder']}")
    if input_value.get("units"):
        lines.append(f"    Units: {input_value['units']}")
    if input_value.get("options"):
        lines.append("    Options: " + ", ".join(input_value["options"]))
    if group := input_value.get("group"):
        lines.append(f"    Group: {group['label']} ({group['id']}, order {group['order']})")
    if activation := input_value.get("activation"):
        lines.append(
            "    Active when: " + activation["key"] + " in {" + ", ".join(activation["equals_any"]) + "}"
        )


def cardinality_text(requirement: dict[str, Any]) -> str:
    minimum = requirement["min_selected"]
    maximum = requirement["max_selected"]
    keys = requirement["keys"]
    if minimum == maximum:
        return f"exactly {minimum}"
    if maximum == len(keys):
        return f"at least {minimum}"
    return f"{minimum} to {maximum}"


def render_program(lines: list[str], name: str, program: dict[str, Any]) -> None:
    metadata = program["program"]
    lines.extend(("=" * 96, f"{metadata['executable']}:{name}"))
    lines.append(f"Title: {metadata['display_name']}")
    lines.append(f"Category: {metadata['category_display_name']} ({metadata['category']}, order {metadata['category_order']})")
    lines.append(f"Visibility: {metadata['visibility']}")
    lines.append(f"Summary: {metadata['summary']}")
    lines.append(f"Help: {metadata['help']}")
    groups = metadata.get("groups", [])
    if groups:
        lines.append("Declared groups: " + ", ".join(f"{group['label']} ({group['id']}, order {group['order']})" for group in groups))
    inputs = input_by_key(program)
    for section, section_title in INPUT_SECTIONS:
        entries = program.get(section, [])
        if not entries:
            continue
        lines.extend(("", section_title))
        for input_value in entries:
            write_input(lines, input_value)
    requirements = metadata.get("requirements", [])
    if requirements:
        lines.extend(("", "INPUT REQUIREMENTS"))
        for requirement in requirements:
            lines.append(f"{requirement['label']}: {requirement['help']}")
            lines.append(f"  Required: {cardinality_text(requirement)}")
            lines.append("  Accepted inputs:")
            for key in requirement["keys"]:
                input_value = inputs[key]
                lines.append(f"    {key}: {input_value['label']}")
    lines.append("")


def render_review(registry: dict[str, Any], executable: Path, registry_output: Path) -> str:
    programs = sorted_programs(registry)
    input_count = sum(
        len(program.get(section, []))
        for program in registry.values()
        for section, _ in INPUT_SECTIONS
    )
    lines = [
        "SIMPLE UI registry review",
        f"Generated: {dt.datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"Registry command: {executable} prg=print_ui_json",
        f"Exact registry JSON: {registry_output}",
        f"Programs: {len(programs)}",
        f"Program-input instances: {input_count}",
        "",
        "This is a read-only UI metadata report. It does not execute workflows or CTest cases.",
        "",
    ]
    for name, program in programs:
        render_program(lines, name, program)
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    executable = private_executable(args.build_dir)
    if not args.build_dir.is_dir():
        print(f"error: build directory does not exist: {args.build_dir}", file=sys.stderr)
        return 2
    if not executable.is_file():
        print(f"error: registry executable not found: {executable}", file=sys.stderr)
        return 2
    try:
        raw_json, registry = load_registry(executable, args.timeout)
    except RuntimeError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    try:
        captures = (
            [
                (name, program, capture_cli_instructions(args.build_dir, name, program, args.timeout))
                for name, program in sorted_programs(registry)
            ]
            if args.capture_cli
            else []
        )
    except RuntimeError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.registry_output.parent.mkdir(parents=True, exist_ok=True)
    if captures:
        args.cli_output.parent.mkdir(parents=True, exist_ok=True)
        args.cli_output.write_text(render_cli_instructions(captures), encoding="utf-8")
    args.registry_output.write_text(raw_json, encoding="utf-8")
    args.output.write_text(render_review(registry, executable, args.registry_output), encoding="utf-8")
    print(f"wrote {args.output} and {args.registry_output} ({len(registry)} programs)")
    if captures:
        print(f"wrote ANSI-free CLI output for {len(captures)} programs to {args.cli_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
