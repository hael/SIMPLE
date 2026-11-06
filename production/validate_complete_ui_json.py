#!/usr/bin/env python3
"""Validate the full UI JSON registry emitted by a SIMPLE executable.

The stream launcher's JSON (``prg=print_ui_stream``, written by hand in
``print_stream_ui_json``) is validated against the registry too: its processes are
simple_stream programs, the user inputs a process takes are defined, are inputs of its
program and carry that input's default. With the stream contract it is also compared
with its reviewed reference, ``stream_launcher_contract.json``.

This is intentionally standard-library-only.  CMake's ``string(JSON)`` reparses
the complete registry for every field, which makes a complete contract check
unacceptably slow in the normal build path.
"""

from __future__ import annotations

import argparse
import difflib
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
UI_DEFAULT_MAX_SIGNIFICANT_DIGITS = 6
PATH_PLACEHOLDER_PREFIX = "e.g. "
FORBIDDEN_PATH_PLACEHOLDERS = {"e.g. input.file", "e.g. input.mrc"}
PATH_PLACEHOLDER_SUFFIXES = {
    "blocktree": (".bin",),
    "boxfile": (".box",),
    "boxtab": (".txt",),
    "ciffile": (".cif",),
    "deftab": (".txt", ".simple"),
    "deselfile": (".txt",),
    "filetab": (".txt",),
    "fsc": (".bin",),
    "gainref": (".mrc", ".mrcs"),
    "infile": (".txt", ".bin", ".star"),
    "oritab": (".txt", ".simple"),
    "oritab2": (".txt", ".simple"),
    "outstk": (".mrc", ".mrcs", ".spi"),
    "outvol": (".mrc", ".spi"),
    "pdbfile": (".pdb",),
    "pdbfile2": (".pdb",),
    "pdbfiles": (".txt",),
    "pdbout": (".pdb",),
    "pickrefs": (".mrc", ".mrcs", ".spi"),
    "plaintexttab": (".txt",),
    "projtab": (".txt",),
    "refs": (".mrc", ".mrcs", ".spi"),
    "star_mic": (".star",),
    "star_model": (".star",),
    "star_ptcl": (".star",),
    "starfile": (".star",),
    "stk": (".mrc", ".mrcs", ".spi"),
    "stk2": (".mrc", ".mrcs", ".spi"),
    "stk_backgr": (".mrc", ".mrcs"),
    "stktab": (".txt",),
    "vol1": (".mrc", ".spi"),
    "vol2": (".mrc", ".spi"),
    "vol3": (".mrc", ".spi"),
    "vol_even": (".mrc", ".spi"),
    "vol_odd": (".mrc", ".spi"),
}
STREAM_UI_CONTRACT_REFERENCE = Path(__file__).resolve().with_name("stream_ui_contract.json")
STREAM_LAUNCHER_CONTRACT_REFERENCE = STREAM_UI_CONTRACT_REFERENCE.with_name("stream_launcher_contract.json")
STREAM_LAUNCHER_CONTRACT_DISPLAY_NAME = str(
    STREAM_LAUNCHER_CONTRACT_REFERENCE.relative_to(STREAM_LAUNCHER_CONTRACT_REFERENCE.parent.parent)
)
STREAM_UI_CONTRACT_DISPLAY_NAME = str(
    STREAM_UI_CONTRACT_REFERENCE.relative_to(STREAM_UI_CONTRACT_REFERENCE.parent.parent)
)


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


def validate_placeholder(input_entry: dict[str, Any], context: str) -> None:
    key = required_string(input_entry, "key", context)
    keytype = required_string(input_entry, "keytype", context)
    placeholder = required_string(input_entry, "placeholder", context, nonempty=False)
    if keytype not in {"file", "dir"}:
        return
    if not placeholder.startswith(PATH_PLACEHOLDER_PREFIX):
        fail(context, f"path placeholder must start with '{PATH_PLACEHOLDER_PREFIX}'")
    if placeholder in FORBIDDEN_PATH_PLACEHOLDERS:
        fail(context, f"path placeholder '{placeholder}' is generic and does not identify the accepted input")
    if any(character in placeholder for character in "{}|"):
        fail(context, "path placeholder must be a single example without defaults or alternatives")

    expected_suffixes = PATH_PLACEHOLDER_SUFFIXES.get(key)
    if key.startswith("projfile"):
        expected_suffixes = (".simple",)
    if expected_suffixes is None:
        return
    if not any(suffix in placeholder.lower() for suffix in expected_suffixes):
        expected = " or ".join(expected_suffixes)
        fail(context, f"placeholder for '{key}' must show {expected}")


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
            validate_placeholder(input_object, input_context)
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


def canonical_stream_contract(registry: dict[str, Any]) -> str:
    """Return the reviewable, deterministic stream-GUI descriptor JSON."""
    stream_registry = {
        name: entry
        for name, entry in registry.items()
        if entry["program"]["executable"] == "simple_stream"
    }
    if not stream_registry:
        fail("stream UI registry", "contains no programs")
    return json.dumps(stream_registry, indent=2, sort_keys=True) + "\n"


def program_inputs(entry: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """The inputs of one registry program, by key."""
    return {
        input_entry["key"]: input_entry
        for section in SECTIONS
        for input_entry in entry.get(section, [])
    }


def same_default(launcher_default: Any, ui_default: Any) -> bool:
    """A launcher default (often a string) against the descriptor's (often a number)."""
    try:
        a, b = float(launcher_default), float(ui_default)
    except (TypeError, ValueError):
        return str(launcher_default) == str(ui_default)
    return abs(a - b) <= 1e-6 * max(1.0, abs(a), abs(b))


def validate_stream_launcher(launcher: Any, registry: dict[str, Any]) -> None:
    """The hand-written stream launcher (print_ui_stream) against the program registry."""
    context = "stream launcher (print_ui_stream)"
    launcher = object_value(launcher, context)
    if set(launcher) != {"user_inputs", "processes"}:
        fail(context, f"has keys {sorted(launcher)}, expected user_inputs and processes")
    user_inputs: dict[str, dict[str, Any]] = {}
    for index, item in enumerate(array_value(launcher["user_inputs"], f"{context}.user_inputs")):
        item_context = f"{context}.user_inputs[{index}]"
        item = object_value(item, item_context)
        key = required_string(item, "key", item_context)
        for field in ("keytype", "label", "help"):
            required_string(item, field, f"{item_context}.{key}")
        bool_value(item.get("required"), f"{item_context}.{key}.required")
        if key in user_inputs:
            fail(item_context, f"user input {key!r} is defined twice")
        user_inputs[key] = item
    used: set[str] = set()
    names: set[str] = set()
    for index, process in enumerate(array_value(launcher["processes"], f"{context}.processes")):
        process_context = f"{context}.processes[{index}]"
        process = object_value(process, process_context)
        name = required_string(process, "name", process_context)
        prg = required_string(process, "prg", f"{process_context}.{name}")
        number_value(process.get("nthr_master"), f"{process_context}.{name}.nthr_master")
        if name in names:
            fail(process_context, f"process name {name!r} is used twice")
        names.add(name)
        entry = registry.get(prg)
        if entry is None or entry["program"]["executable"] != "simple_stream":
            fail(f"{process_context}.{name}", f"prg={prg} is no simple_stream program of the registry")
        inputs = program_inputs(entry)
        for key in array_value(process.get("user_inputs", []), f"{process_context}.{name}.user_inputs"):
            key = string_value(key, f"{process_context}.{name}.user_inputs")
            if key not in user_inputs:
                fail(f"{process_context}.{name}", f"user input {key!r} is not defined in user_inputs")
            if key not in inputs:
                fail(f"{process_context}.{name}", f"user input {key!r} is no input of program {prg}")
            used.add(key)
            launcher_default = user_inputs[key].get("default")
            ui_default = inputs[key].get("default")
            if launcher_default is not None and inputs[key].get("has_default") and \
                    not same_default(launcher_default, ui_default):
                fail(f"{process_context}.{name}",
                     f"user input {key!r} defaults to {launcher_default!r}, program {prg} to {ui_default!r}")
        for item in array_value(process.get("static_inputs", []), f"{process_context}.{name}.static_inputs"):
            item = string_value(item, f"{process_context}.{name}.static_inputs")
            key, sep, _ = item.partition("=")
            if not sep or not key:
                fail(f"{process_context}.{name}", f"static input {item!r} is not key=value")
    unused = sorted(set(user_inputs) - used)
    if unused:
        fail(context, f"user inputs no process takes: {unused}")


def canonical_stream_launcher(launcher: Any) -> str:
    """Return the reviewable, deterministic stream launcher JSON."""
    return json.dumps(launcher, indent=2, sort_keys=True) + "\n"


def validate_stream_contract(
    actual: str,
    reference: Path = STREAM_UI_CONTRACT_REFERENCE,
    display_name: str = STREAM_UI_CONTRACT_DISPLAY_NAME,
) -> None:
    """Reject an unreviewed stream-GUI change and show its JSON diff."""
    try:
        expected = reference.read_text(encoding="utf-8")
    except FileNotFoundError as error:
        raise ValueError(
            "stream UI contract: missing review reference "
            f"'{display_name}'"
        ) from error

    if expected == actual:
        return

    diff = "".join(
        difflib.unified_diff(
            expected.splitlines(keepends=True),
            actual.splitlines(keepends=True),
            fromfile=display_name,
            tofile="generated stream UI JSON",
        )
    )
    if not diff:
        diff = "(The reference and generated JSON differ only in line endings.)\n"
    fail(
        "stream UI contract",
        "The generated stream GUI differs from its reviewed JSON reference.\n"
        "The ordinary build has already confirmed that the complete UI JSON is valid.\n"
        "Review the readable diff below. If the GUI change is intended, run:\n"
        "  cmake --build <build-dir> --target update_stream_ui_contract\n"
        "Then commit the descriptor and updated reference file together.\n"
        "Do not update the reference for an unintended change.\n\n"
        f"{diff}",
    )


def numeric_significant_digits(value: str) -> int:
    """Return the decimal significant digits in an already-valid JSON number."""
    mantissa = value.lstrip("+-").lower().split("e", maxsplit=1)[0]
    digits = mantissa.replace(".", "").lstrip("0")
    return len(digits) if digits else 1


def validate_numeric_default_precision(registry: dict[str, Any]) -> None:
    """Keep JSON default values concise instead of exposing binary residue."""
    for program_name, program in registry.items():
        for section in SECTIONS:
            for index, input_entry in enumerate(program[section]):
                if not input_entry.get("has_default"):
                    continue
                if input_entry["keytype"] not in {"num", "int", "float"}:
                    continue
                default = input_entry.get("default")
                if not isinstance(default, str):
                    fail(f"program '{program_name}'.{section}[{index}].default", "must be numeric text")
                digits = numeric_significant_digits(default)
                if digits > UI_DEFAULT_MAX_SIGNIFICANT_DIGITS:
                    fail(
                        f"program '{program_name}'.{section}[{index}].default",
                        f"has {digits} significant digits; UI defaults allow at most "
                        f"{UI_DEFAULT_MAX_SIGNIFICANT_DIGITS}",
                    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument(
        "--launcher-executable",
        type=Path,
        help="the executable that prints the stream launcher (prg=print_ui_stream); default --executable",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--stamp", required=True, type=Path)
    stream_contract_action = parser.add_mutually_exclusive_group()
    stream_contract_action.add_argument(
        "--check-stream-contract",
        action="store_true",
        help="also require the reviewed stream-GUI snapshot to match",
    )
    stream_contract_action.add_argument(
        "--update-stream-contract",
        action="store_true",
        help="replace the reviewed stream-GUI reference after intentional review",
    )
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
    launcher_run = subprocess.run(
        [str(args.launcher_executable or args.executable), "prg=print_ui_stream"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if launcher_run.returncode:
        sys.stderr.buffer.write(launcher_run.stderr)
        raise RuntimeError(f"print_ui_stream failed with exit code {launcher_run.returncode}")

    try:
        registry = json.loads(completed.stdout, object_pairs_hook=no_duplicate_object)
        registry = object_value(registry, "complete UI registry")
        if not registry:
            fail("complete UI registry", "contains no programs")
        for program_name, program in registry.items():
            string_value(program_name, "complete UI registry key")
            validate_program(program_name, program)
        launcher = json.loads(launcher_run.stdout, object_pairs_hook=no_duplicate_object)
        validate_stream_launcher(launcher, registry)
        stream_contract = None
        if args.check_stream_contract or args.update_stream_contract:
            stream_contract = canonical_stream_contract(registry)
            launcher_contract = canonical_stream_launcher(launcher)
            if args.update_stream_contract:
                STREAM_UI_CONTRACT_REFERENCE.write_text(stream_contract, encoding="utf-8")
                STREAM_LAUNCHER_CONTRACT_REFERENCE.write_text(launcher_contract, encoding="utf-8")
            else:
                # both references are compared, so a stale one does not hide the other's diff
                failures = []
                for actual, reference, display_name in (
                    (stream_contract, STREAM_UI_CONTRACT_REFERENCE, STREAM_UI_CONTRACT_DISPLAY_NAME),
                    (launcher_contract, STREAM_LAUNCHER_CONTRACT_REFERENCE, STREAM_LAUNCHER_CONTRACT_DISPLAY_NAME),
                ):
                    try:
                        validate_stream_contract(actual, reference, display_name)
                    except ValueError as error:
                        failures.append(str(error))
                if failures:
                    raise ValueError("\n\n".join(failures))
        textual_registry = json.loads(
            completed.stdout,
            object_pairs_hook=no_duplicate_object,
            parse_float=lambda value: value,
            parse_int=lambda value: value,
        )
        validate_numeric_default_precision(object_value(textual_registry, "complete UI registry"))
    except (json.JSONDecodeError, ValueError) as error:
        raise SystemExit(f"UI JSON validation failed: {error}") from error

    if stream_contract is None:
        args.output.write_bytes(completed.stdout)
    else:
        args.output.write_text(stream_contract, encoding="utf-8")
    args.stamp.touch()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
