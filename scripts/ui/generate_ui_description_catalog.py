#!/usr/bin/env python3
"""Validate the UI description catalog and generate its Fortran lookup module."""

from __future__ import annotations

import argparse
import json
import sys
import tomllib
from pathlib import Path
from typing import Any


VISIBILITIES = {"standard": "UI_VIS_STANDARD", "advanced": "UI_VIS_ADVANCED", "developer": "UI_VIS_DEVELOPER"}
LAYOUTS = {"grouped", "flat"}
TEXT_LIMITS = {"display_name": 80, "summary": 120, "label": 80, "help": 500, "placeholder": 48, "units": 24}
SUITE_FIELDS = {"id", "executable", "display_name", "layout", "order"}
GROUP_FIELDS = {"suite_id", "group_id", "group_title", "group_order", "review_status", "program"}
PROGRAM_FIELDS = {"name", "executable", "display_name", "summary", "help", "visibility", "input"}
INPUT_FIELDS = {"key", "label", "help", "placeholder", "units", "visibility", "choices"}


def reject_unknown_fields(record: dict[str, Any], allowed: set[str], context: str) -> None:
    unknown = sorted(set(record) - allowed)
    if unknown:
        raise ValueError(f"{context}: unrecognized fields: {', '.join(unknown)}")


def normalize_choices(input_record: dict[str, Any], context: str) -> list[dict[str, str]]:
    choices = input_record.pop("choices", [])
    if not isinstance(choices, list):
        raise ValueError(f"{context}: choices must be an array")
    normalized: list[dict[str, str]] = []
    for choice in choices:
        if isinstance(choice, str):
            normalized.append({"value": choice, "label": choice, "help": ""})
        elif isinstance(choice, dict):
            unknown = sorted(set(choice) - {"value", "label", "help"})
            if unknown:
                raise ValueError(f"{context}: choice has unrecognized fields: {', '.join(unknown)}")
            value = choice.get("value")
            if not isinstance(value, str):
                raise ValueError(f"{context}: expanded choice value must be a string")
            label = choice.get("label", value)
            help_text = choice.get("help", "")
            if not isinstance(label, str) or not isinstance(help_text, str):
                raise ValueError(f"{context}: expanded choice label and help must be strings")
            normalized.append({"value": value, "label": label, "help": help_text})
        else:
            raise ValueError(f"{context}: each choice must be a string or inline table")
    return normalized


def normalize_program(program: Any, context: str) -> dict[str, Any]:
    if not isinstance(program, dict):
        raise ValueError(f"{context}: program is not a TOML table")
    normalized = dict(program)
    reject_unknown_fields(normalized, PROGRAM_FIELDS, context)
    normalized.setdefault("help", "")
    normalized.setdefault("visibility", "standard")
    inputs = normalized.get("input", [])
    if not isinstance(inputs, list):
        raise ValueError(f"{context}: input must be an array of tables")
    normalized_inputs: list[dict[str, Any]] = []
    for input_index, input_record in enumerate(inputs, start=1):
        input_context = f"{context} input {input_index}"
        if not isinstance(input_record, dict):
            raise ValueError(f"{input_context}: input is not a TOML table")
        normalized_input = dict(input_record)
        reject_unknown_fields(normalized_input, INPUT_FIELDS, input_context)
        normalized_input.setdefault("help", "")
        normalized_input.setdefault("placeholder", "")
        normalized_input.setdefault("units", "")
        normalized_input.setdefault("visibility", "standard")
        normalized_input["choice"] = normalize_choices(normalized_input, input_context)
        normalized_inputs.append(normalized_input)
    normalized["input"] = normalized_inputs
    return normalized


def parse_catalog(path: Path) -> list[dict[str, Any]]:
    if not path.is_dir():
        raise ValueError(f"catalog path must be a directory: {path}")

    records: list[dict[str, Any]] = []
    suite_files = sorted(path.glob("*/suite.toml"))
    if not suite_files:
        raise ValueError(f"catalog directory has no suite.toml files: {path}")
    suite_ids: set[str] = set()
    suite_orders: set[int] = set()
    for suite_file in suite_files:
        suite = tomllib.loads(suite_file.read_text(encoding="utf-8"))
        reject_unknown_fields(suite, SUITE_FIELDS, str(suite_file))
        suite_id = suite.get("id")
        if not isinstance(suite_id, str) or not suite_id:
            raise ValueError(f"{suite_file}: suite id must be a non-empty string")
        if suite.get("executable") != suite_id:
            raise ValueError(f"{suite_file}: executable must match suite id")
        if not isinstance(suite.get("display_name"), str) or not suite["display_name"]:
            raise ValueError(f"{suite_file}: display_name must be a non-empty string")
        if suite.get("layout") not in LAYOUTS:
            raise ValueError(f"{suite_file}: layout must be grouped or flat")
        if not isinstance(suite.get("order"), int) or suite["order"] < 0:
            raise ValueError(f"{suite_file}: order must be a non-negative integer")
        if suite_id in suite_ids:
            raise ValueError(f"{suite_file}: duplicate suite id {suite_id}")
        if suite["order"] in suite_orders:
            raise ValueError(f"{suite_file}: duplicate suite order {suite['order']}")
        suite_ids.add(suite_id)
        suite_orders.add(suite["order"])

        group_files = sorted(group_file for group_file in suite_file.parent.glob("*.toml") if group_file.name != "suite.toml")
        group_ids: set[str] = set()
        group_orders: set[int] = set()
        for group_file in group_files:
            group = tomllib.loads(group_file.read_text(encoding="utf-8"))
            reject_unknown_fields(group, GROUP_FIELDS, str(group_file))
            if group.get("suite_id") != suite_id:
                raise ValueError(f"{group_file}: suite_id does not match {suite_file}")
            if not isinstance(group.get("group_id"), str) or not group["group_id"]:
                raise ValueError(f"{group_file}: group_id must be a non-empty string")
            if not isinstance(group.get("group_title"), str) or not group["group_title"]:
                raise ValueError(f"{group_file}: group_title must be a non-empty string")
            if not isinstance(group.get("group_order"), int) or group["group_order"] < 0:
                raise ValueError(f"{group_file}: group_order must be a non-negative integer")
            if group["group_id"] in group_ids:
                raise ValueError(f"{group_file}: duplicate group_id {group['group_id']}")
            if group["group_order"] in group_orders:
                raise ValueError(f"{group_file}: duplicate group_order {group['group_order']}")
            group_ids.add(group["group_id"])
            group_orders.add(group["group_order"])
            group.setdefault("review_status", "legacy")
            if group.get("review_status") not in {"legacy", "reviewed"}:
                raise ValueError(f"{group_file}: review_status must be legacy or reviewed")
            programs = group.get("program")
            if not isinstance(programs, list) or not programs:
                raise ValueError(f"{group_file}: group must contain at least one program")
            for program_index, program in enumerate(programs, start=1):
                program = normalize_program(program, f"{group_file} program {program_index}")
                program["_review_status"] = group["review_status"]
                records.append(program)
    return records


def strings(record: dict[str, Any], fields: tuple[str, ...], context: str, errors: list[str], strict_text: bool,
            allow_empty: set[str] | None = None) -> None:
    allow_empty = allow_empty or set()
    for field in fields:
        value = record.get(field)
        if not isinstance(value, str):
            errors.append(f"{context}: {field} must be a string")
            continue
        if field not in allow_empty and not value:
            errors.append(f"{context}: {field} must not be empty")
        if "\n" in value:
            errors.append(f"{context}: {field} must be one line")
        if strict_text and field in TEXT_LIMITS and len(value) > TEXT_LIMITS[field]:
            errors.append(f"{context}: {field} exceeds {TEXT_LIMITS[field]} characters")


def validate(records: list[dict[str, Any]], strict_text: bool) -> list[str]:
    errors: list[str] = []
    identities: set[tuple[str, str]] = set()
    for program in records:
        context = f"program {program.get('executable', '?')}/{program.get('name', '?')}"
        reviewed = strict_text or program.get("_review_status") == "reviewed"
        program_allow_empty = set() if reviewed else {"help"}
        strings(program, ("name", "executable", "display_name", "summary", "help", "visibility"), context, errors, reviewed,
                program_allow_empty)
        identity = (program.get("executable", ""), program.get("name", ""))
        if identity in identities:
            errors.append(f"{context}: duplicate program identity")
        identities.add(identity)
        if program.get("visibility") not in VISIBILITIES:
            errors.append(f"{context}: visibility must be standard, advanced, or developer")
        input_keys: set[str] = set()
        for input_record in program.get("input", []):
            input_context = f"{context} input {input_record.get('key', '?')}"
            input_allow_empty = {"placeholder", "units"}
            if not reviewed:
                input_allow_empty.add("help")
            strings(input_record, ("key", "label", "help", "placeholder", "units", "visibility"), input_context, errors,
                    reviewed, input_allow_empty)
            key = input_record.get("key", "")
            if key in input_keys:
                errors.append(f"{input_context}: duplicate input key")
            input_keys.add(key)
            if input_record.get("visibility") not in VISIBILITIES:
                errors.append(f"{input_context}: visibility must be standard, advanced, or developer")
            if reviewed and ("{default}" in input_record.get("placeholder", "") or "(" in input_record.get("placeholder", "")):
                errors.append(f"{input_context}: placeholder contains legacy default or choice encoding")
            choice_values: set[str] = set()
            for choice in input_record.get("choice", []):
                choice_context = f"{input_context} choice {choice.get('value', '?')}"
                strings(choice, ("value", "label", "help"), choice_context, errors, reviewed, {"help"})
                value = choice.get("value", "")
                if value in choice_values:
                    errors.append(f"{choice_context}: duplicate value")
                choice_values.add(value)
    return errors


def fortran_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def generate(records: list[dict[str, Any]]) -> str:
    lines = [
        "! This file is generated by scripts/ui/generate_ui_description_catalog.py.",
        "! Do not edit manually; edit src/main/ui/catalog/ and regenerate.",
        "module simple_ui_description_catalog",
        "use simple_ui_descriptor_types, only: ui_choice",
        "use simple_ui_visibility, only: UI_VIS_STANDARD, UI_VIS_ADVANCED, UI_VIS_DEVELOPER",
        "implicit none",
        "private",
        "public :: get_ui_program_presentation, get_ui_input_presentation, apply_ui_input_choice_presentation",
        "contains",
        "",
        "subroutine get_ui_program_presentation(executable, name, display_name, summary, help, visibility, found)",
        "    character(len=*),              intent(in)  :: executable, name",
        "    character(len=:), allocatable, intent(out) :: display_name, summary, help",
        "    integer,                       intent(out) :: visibility",
        "    logical,                       intent(out) :: found",
        "    display_name = ''",
        "    summary = ''",
        "    help = ''",
        "    visibility = UI_VIS_DEVELOPER",
        "    found = .false.",
        "    select case(trim(executable)//'/'//trim(name))",
    ]
    for program in records:
        lines.extend((f"    case({fortran_string(program['executable'] + '/' + program['name'])})", "        found = .true.",
                      f"        display_name = {fortran_string(program['display_name'])}",
                      f"        summary = {fortran_string(program['summary'])}",
                      f"        help = {fortran_string(program['help'])}",
                      f"        visibility = {VISIBILITIES[program['visibility']]}"))
    lines.extend(("    end select", "end subroutine get_ui_program_presentation", "",
                  "subroutine get_ui_input_presentation(executable, program_name, input_key, label, help, placeholder, units, visibility, found)",
                  "    character(len=*),              intent(in)  :: executable, program_name, input_key",
                  "    character(len=:), allocatable, intent(out) :: label, help, placeholder, units",
                  "    integer,                       intent(out) :: visibility",
                  "    logical,                       intent(out) :: found",
                  "    label = ''", "    help = ''", "    placeholder = ''", "    units = ''",
                  "    visibility = UI_VIS_DEVELOPER", "    found = .false.", "    select case(trim(executable)//'/'//trim(program_name))"))
    for program in records:
        lines.append(f"    case({fortran_string(program['executable'] + '/' + program['name'])})")
        lines.append("        select case(trim(input_key))")
        for input_record in program.get("input", []):
            lines.extend((f"        case({fortran_string(input_record['key'])})", "            found = .true.",
                          f"            label = {fortran_string(input_record['label'])}",
                          f"            help = {fortran_string(input_record['help'])}",
                          f"            placeholder = {fortran_string(input_record['placeholder'])}",
                          f"            units = {fortran_string(input_record['units'])}",
                          f"            visibility = {VISIBILITIES[input_record['visibility']]}"))
        lines.extend(("        end select",))
    lines.extend(("    end select", "end subroutine get_ui_input_presentation", "",
                  "subroutine apply_ui_input_choice_presentation(executable, program_name, input_key, choices, found)",
                  "    character(len=*), intent(in) :: executable, program_name, input_key",
                  "    type(ui_choice),  intent(inout) :: choices(:)",
                  "    logical,          intent(out) :: found",
                  "    integer :: i", "    found = .false.", "    select case(trim(executable)//'/'//trim(program_name))"))
    for program in records:
        choice_inputs = [entry for entry in program.get("input", []) if entry.get("choice")]
        if not choice_inputs:
            continue
        lines.append(f"    case({fortran_string(program['executable'] + '/' + program['name'])})")
        lines.append("        select case(trim(input_key))")
        for input_record in choice_inputs:
            lines.extend((f"        case({fortran_string(input_record['key'])})", "            found = .true.",
                          "            do i = 1, size(choices)", "                select case(choices(i)%value%to_char())"))
            for choice in input_record["choice"]:
                lines.extend((f"                case({fortran_string(choice['value'])})",
                              f"                    choices(i)%label = {fortran_string(choice['label'])}",
                              f"                    choices(i)%help = {fortran_string(choice['help'])}"))
            lines.extend(("                end select", "            enddo"))
        lines.append("        end select")
    lines.extend(("    end select", "end subroutine apply_ui_input_choice_presentation", "", "end module simple_ui_description_catalog", ""))
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true", help="fail when generated output is missing or stale")
    parser.add_argument("--strict-text", action="store_true", help="enforce user-facing text length and placeholder rules")
    args = parser.parse_args()
    try:
        records = parse_catalog(args.catalog)
    except (OSError, ValueError, tomllib.TOMLDecodeError) as error:
        print(f"UI description catalog validation failed: {error}", file=sys.stderr)
        return 1
    errors = validate(records, args.strict_text)
    if errors:
        print("UI description catalog validation failed:", file=sys.stderr)
        print("\n".join(f"- {error}" for error in errors), file=sys.stderr)
        return 1
    output = generate(records)
    if args.check:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != output:
            print(f"generated catalog is stale: {args.output}", file=sys.stderr)
            return 1
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
