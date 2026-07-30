#!/usr/bin/env python3
"""Validate the full UI JSON registry emitted by ``simple_private_exec``.

This is intentionally standard-library-only.  CMake's ``string(JSON)`` reparses
the complete registry for every field, which makes a complete contract check
unacceptably slow in the normal build path.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


SECTIONS = (
    "image input/output",
    "file input/output",
    "parameter input/output",
    "search controls",
    "filter controls",
    "mask controls",
    "computer controls",
)
VISIBILITIES = {"standard", "advanced", "developer"}


def fail(context: str, message: str) -> None:
    raise ValueError(f"{context}: {message}")


def no_duplicate_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            fail("JSON", f"duplicate object member '{key}'")
        result[key] = value
    return result


def object_value(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        fail(context, "must be an object")
    return value


def array_value(value: Any, context: str) -> list[Any]:
    if not isinstance(value, list):
        fail(context, "must be an array")
    return value


def string_value(value: Any, context: str, *, nonempty: bool = True) -> str:
    if not isinstance(value, str):
        fail(context, "must be a string")
    if nonempty and not value:
        fail(context, "must not be empty")
    return value


def bool_value(value: Any, context: str) -> bool:
    if not isinstance(value, bool):
        fail(context, "must be a boolean")
    return value


def number_value(value: Any, context: str) -> int | float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        fail(context, "must be a number")
    return value


def required_string(entry: dict[str, Any], key: str, context: str, *, nonempty: bool = True) -> str:
    if key not in entry:
        fail(context, f"missing '{key}'")
    return string_value(entry[key], f"{context}.{key}", nonempty=nonempty)


def validate_visibility(entry: dict[str, Any], context: str) -> str:
    visibility = required_string(entry, "visibility", context)
    if visibility not in VISIBILITIES:
        fail(context, f"invalid visibility '{visibility}'")
    return visibility


def validate_group(group: Any, context: str) -> str:
    group_object = object_value(group, context)
    group_id = required_string(group_object, "id", context)
    required_string(group_object, "label", context)
    if "order" not in group_object:
        fail(context, "missing 'order'")
    number_value(group_object["order"], f"{context}.order")
    return group_id


def validate_program(program_name: str, entry: Any) -> None:
    context = f"program '{program_name}'"
    program = object_value(entry, context)
    descriptor = object_value(program.get("program"), f"{context}.program")
    for key in ("name", "category", "category_display_name", "display_name", "summary", "help", "executable"):
        required_string(descriptor, key, f"{context}.program")
    if required_string(descriptor, "name", f"{context}.program") != program_name:
        fail(context, "descriptor name does not match its registry key")
    validate_visibility(descriptor, f"{context}.program")
    if "category_order" not in descriptor:
        fail(f"{context}.program", "missing 'category_order'")
    number_value(descriptor["category_order"], f"{context}.program.category_order")

    group_ids: set[str] = set()
    if "groups" in descriptor:
        for index, group in enumerate(array_value(descriptor["groups"], f"{context}.program.groups")):
            group_id = validate_group(group, f"{context}.program.groups[{index}]")
            if group_id in group_ids:
                fail(context, f"duplicate group id '{group_id}'")
            group_ids.add(group_id)

    input_keys: set[str] = set()
    activation_keys: list[str] = []
    for section in SECTIONS:
        for index, input_entry in enumerate(array_value(program.get(section), f"{context}.{section}")):
            input_context = f"{context}.{section}[{index}]"
            input_object = object_value(input_entry, input_context)
            for key in ("key", "keytype", "label", "help", "visibility"):
                required_string(input_object, key, input_context)
            required_string(input_object, "placeholder", input_context, nonempty=False)
            required = bool_value(input_object.get("required"), f"{input_context}.required")
            bool_value(input_object.get("has_default"), f"{input_context}.has_default")
            array_value(input_object.get("options"), f"{input_context}.options")
            visibility = validate_visibility(input_object, input_context)
            input_key = required_string(input_object, "key", input_context)
            if required and visibility != "standard":
                fail(input_context, "required input must have standard visibility")
            if input_key in input_keys:
                fail(context, f"duplicate input key '{input_key}'")
            input_keys.add(input_key)

            if "group" in input_object:
                group_id = validate_group(input_object["group"], f"{input_context}.group")
                if group_id not in group_ids:
                    fail(input_context, f"references unknown group '{group_id}'")
            if "activation" in input_object:
                activation = object_value(input_object["activation"], f"{input_context}.activation")
                activation_keys.append(required_string(activation, "key", f"{input_context}.activation"))
                array_value(activation.get("equals_any"), f"{input_context}.activation.equals_any")

    for activation_key in activation_keys:
        if activation_key not in input_keys:
            fail(context, f"activation controller '{activation_key}' is not an input")

    requirements = array_value(descriptor.get("requirements"), f"{context}.program.requirements")
    requirement_ids: set[str] = set()
    for index, requirement_entry in enumerate(requirements):
        requirement_context = f"{context}.program.requirements[{index}]"
        requirement = object_value(requirement_entry, requirement_context)
        requirement_id = required_string(requirement, "id", requirement_context)
        if requirement_id in requirement_ids:
            fail(context, f"duplicate requirement id '{requirement_id}'")
        requirement_ids.add(requirement_id)
        required_string(requirement, "label", requirement_context)
        required_string(requirement, "help", requirement_context)
        keys = array_value(requirement.get("keys"), f"{requirement_context}.keys")
        if not keys:
            fail(requirement_context, "key list must not be empty")
        for key in keys:
            key = string_value(key, f"{requirement_context}.keys item")
            if key not in input_keys:
                fail(requirement_context, f"references unknown input '{key}'")
        number_value(requirement.get("min_selected"), f"{requirement_context}.min_selected")
        number_value(requirement.get("max_selected"), f"{requirement_context}.max_selected")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--stamp", required=True, type=Path)
    args = parser.parse_args()

    completed = subprocess.run(
        [str(args.executable), "prg=print_ui_json"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode:
        sys.stderr.buffer.write(completed.stderr)
        raise RuntimeError(f"print_ui_json failed with exit code {completed.returncode}")

    try:
        registry = json.loads(completed.stdout, object_pairs_hook=no_duplicate_object)
        registry = object_value(registry, "complete UI registry")
        if not registry:
            fail("complete UI registry", "contains no programs")
        for program_name, program in registry.items():
            string_value(program_name, "complete UI registry key")
            validate_program(program_name, program)
    except (json.JSONDecodeError, ValueError) as error:
        raise SystemExit(f"UI JSON validation failed: {error}") from error

    args.output.write_bytes(completed.stdout)
    args.stamp.touch()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
