#!/usr/bin/env python3
"""Split the legacy UI description catalog into suite and group TOML files.

This is a one-time migration tool. It consumes the current monolithic catalog
and the executable's ``prg=list``/``test=list`` output, preserving every
existing description record while making program-group membership explicit.
It must not be used to refresh a catalog after cutover.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tomllib
from collections import defaultdict
from pathlib import Path
from typing import Any


FENCE = re.compile(r"```toml\s*\n(.*?)```", re.DOTALL)
ANSI_ESCAPE = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")
SUITES = {
    "simple_exec": {"display_name": "SIMPLE", "layout": "grouped", "order": 10, "argument": "prg=list"},
    "single_exec": {"display_name": "SINGLE", "layout": "grouped", "order": 20, "argument": "prg=list"},
    "simple_test_exec": {"display_name": "SIMPLE tests", "layout": "grouped", "order": 30, "argument": "test=list"},
    "simple_stream": {"display_name": "SIMPLE stream", "layout": "flat", "order": 40, "argument": "prg=list"},
}


def quote(value: str) -> str:
    """Return a TOML basic string using JSON's compatible escaping."""
    return json.dumps(value, ensure_ascii=False)


def toml_scalar(value: Any) -> str:
    if isinstance(value, str):
        return quote(value)
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return repr(value)
    raise ValueError(f"unsupported legacy catalog value: {value!r}")


def slug(value: str) -> str:
    result = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    if not result:
        raise ValueError(f"cannot derive a group ID from {value!r}")
    return result


def parse_legacy_catalog(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for match in FENCE.finditer(path.read_text(encoding="utf-8")):
        record = tomllib.loads(match.group(1))
        if "name" not in record or "executable" not in record:
            continue
        if not isinstance(record["name"], str) or not isinstance(record["executable"], str):
            raise ValueError("legacy program name and executable must be strings")
        records.append(record)
    return records


def parse_program_list(output: str) -> list[tuple[str, list[str]]]:
    groups: list[tuple[str, list[str]]] = []
    title: str | None = None
    programs: list[str] = []
    for raw_line in output.splitlines():
        line = ANSI_ESCAPE.sub("", raw_line).strip()
        if not line:
            continue
        if line.endswith(":"):
            if title is not None:
                groups.append((title, programs))
            title = line[:-1]
            programs = []
            continue
        if title is None:
            raise ValueError(f"program list contains an entry before its first heading: {line!r}")
        programs.append(line)
    if title is not None:
        groups.append((title, programs))
    if not groups:
        raise ValueError("program list did not contain any groups")
    return groups


def read_program_groups(executable: Path, argument: str) -> list[tuple[str, list[str]]]:
    completed = subprocess.run(
        [str(executable), argument],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if completed.returncode != 0:
        raise ValueError(f"{executable} {argument} failed:\n{completed.stdout}{completed.stderr}")
    return parse_program_list(completed.stdout)


def write_suite(path: Path, suite_id: str) -> None:
    suite = SUITES[suite_id]
    lines = [
        f'id = {quote(suite_id)}',
        f'executable = {quote(suite_id)}',
        f'display_name = {quote(suite["display_name"])}',
        f'layout = {quote(suite["layout"])}',
        f'order = {suite["order"]}',
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def write_field(lines: list[str], record: dict[str, Any], key: str, default: Any = None) -> None:
    if key in record and record[key] != default:
        lines.append(f"{key} = {toml_scalar(record[key])}")


def write_choices(lines: list[str], input_record: dict[str, Any]) -> None:
    choices = input_record.get("choice", [])
    if not choices:
        return
    if not all(isinstance(choice, dict) for choice in choices):
        raise ValueError("legacy input choice is not a table")
    if all(choice.get("label") == choice.get("value") and choice.get("help", "") == "" for choice in choices):
        lines.append("choices = [" + ", ".join(quote(choice["value"]) for choice in choices) + "]")
        return
    lines.append("choices = [")
    for choice in choices:
        if not isinstance(choice.get("value"), str):
            raise ValueError("legacy choice value is not a string")
        fields = [f"value = {quote(choice['value'])}"]
        if choice.get("label", choice["value"]) != choice["value"]:
            fields.append(f"label = {quote(choice['label'])}")
        if choice.get("help", ""):
            fields.append(f"help = {quote(choice['help'])}")
        lines.append("    { " + ", ".join(fields) + " },")
    lines.append("]")


def write_program(lines: list[str], program: dict[str, Any]) -> None:
    lines.extend(("", "[[program]]"))
    for key in ("name", "executable", "display_name", "summary"):
        write_field(lines, program, key)
    write_field(lines, program, "help", "")
    write_field(lines, program, "visibility", "standard")
    for input_record in program.get("input", []):
        if not isinstance(input_record, dict):
            raise ValueError(f"{program['executable']}/{program['name']}: input is not a table")
        lines.extend(("", "[[program.input]]"))
        for key in ("key", "label"):
            write_field(lines, input_record, key)
        for key in ("help", "placeholder", "units"):
            write_field(lines, input_record, key, "")
        write_field(lines, input_record, "visibility", "standard")
        write_choices(lines, input_record)


def write_group(path: Path, suite_id: str, title: str, order: int, programs: list[dict[str, Any]]) -> None:
    group_id = slug(title)
    lines = [
        f"# {title}",
        f"suite_id = {quote(suite_id)}",
        f"group_id = {quote(group_id)}",
        f"group_title = {quote(title)}",
        f"group_order = {order}",
    ]
    for program in programs:
        write_program(lines, program)
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def build_membership(program_groups: dict[str, list[tuple[str, list[str]]]]) -> dict[tuple[str, str], tuple[str, str, int, int]]:
    membership: dict[tuple[str, str], tuple[str, str, int, int]] = {}
    for suite_id, groups in program_groups.items():
        for group_index, (title, names) in enumerate(groups, start=1):
            for program_index, name in enumerate(names, start=1):
                key = (suite_id, name)
                if key in membership:
                    raise ValueError(f"{suite_id}: duplicate program in list output: {name}")
                membership[key] = (suite_id, title, group_index, program_index)
    return membership


def suite_for_legacy_record(record: dict[str, Any], membership: dict[tuple[str, str], tuple[str, str, int, int]]) -> str:
    executable = record["executable"]
    name = record["name"]
    if executable == "all":
        candidate = "simple_exec"
    else:
        candidate = executable
    if (candidate, name) not in membership:
        raise ValueError(f"unresolved legacy catalog program: {executable}/{name}")
    return candidate


def import_catalog(legacy_path: Path, output_root: Path, executables: dict[str, Path]) -> tuple[int, int]:
    if output_root.exists() and any(output_root.iterdir()):
        raise ValueError(f"output directory must be empty: {output_root}")
    records = parse_legacy_catalog(legacy_path)
    program_groups = {
        suite_id: read_program_groups(executables[suite_id], suite["argument"])
        for suite_id, suite in SUITES.items()
    }
    membership = build_membership(program_groups)
    grouped_records: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        suite_id = suite_for_legacy_record(record, membership)
        _, title, _, _ = membership[(suite_id, record["name"])]
        grouped_records[(suite_id, title)].append(record)

    output_root.mkdir(parents=True, exist_ok=True)
    for suite_id, groups in program_groups.items():
        suite_records = [record for (suite, _), values in grouped_records.items() if suite == suite_id for record in values]
        if not suite_records:
            continue
        suite_dir = output_root / suite_id
        suite_dir.mkdir()
        write_suite(suite_dir / "suite.toml", suite_id)
        for group_index, (title, names) in enumerate(groups, start=1):
            records_in_group = grouped_records.get((suite_id, title), [])
            if not records_in_group:
                continue
            record_by_name = {record["name"]: record for record in records_in_group}
            ordered_records = [record_by_name[name] for name in names if name in record_by_name]
            write_group(suite_dir / f"{slug(title)}.toml", suite_id, title, group_index * 10, ordered_records)
    return len(records), len(grouped_records)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-catalog", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--simple-exec", type=Path, required=True)
    parser.add_argument("--single-exec", type=Path, required=True)
    parser.add_argument("--simple-test-exec", type=Path, required=True)
    parser.add_argument("--simple-stream", type=Path, required=True)
    args = parser.parse_args()
    executables = {
        "simple_exec": args.simple_exec,
        "single_exec": args.single_exec,
        "simple_test_exec": args.simple_test_exec,
        "simple_stream": args.simple_stream,
    }
    try:
        records, groups = import_catalog(args.legacy_catalog, args.output_root, executables)
    except (OSError, ValueError, subprocess.SubprocessError, tomllib.TOMLDecodeError) as error:
        print(f"legacy UI catalog import failed: {error}", file=sys.stderr)
        return 1
    print(f"imported {records} programs into {groups} group files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
