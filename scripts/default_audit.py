#!/usr/bin/env python3
"""Generate a reviewable audit of SIMPLE command-line defaults.

This tool is intentionally conservative.  It reads Fortran source and the
Fortran-emitted UI registry, but is never read by a SIMPLE executable at
runtime.  Unsupported control flow is reported as unresolved rather than
being interpreted as a default.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SECTION_NAMES = (
    "image input/output",
    "file input/output",
    "parameter input/output",
    "search controls",
    "filter controls",
    "mask controls",
    "computer controls",
)

MISSING_KEY = re.compile(
    r"^\s*if\s*\(\s*\.not\.\s*cline%defined\s*\(\s*'([^']+)'\s*\)\s*\)\s*"
    r"call\s+cline%set\s*\(\s*'([^']+)'\s*,\s*(.+?)\s*\)\s*$",
    re.IGNORECASE,
)
NUMERIC_LITERAL = re.compile(
    r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[ed][+-]?\d+)?(?:_[a-z0-9_]+)?$", re.IGNORECASE
)
SET_CALL = re.compile(
    r"^\s*call\s+cline%set\s*\(\s*'([^']+)'\s*,\s*(.+?)\s*\)\s*$", re.IGNORECASE
)
NEW_CALL = re.compile(r"\b[a-z_][a-z0-9_]*%new\s*\(\s*cline\s*\)", re.IGNORECASE)
PARAMS_INITIALIZER_CALL = re.compile(
    r"\b[a-z_][a-z0-9_]*%init_params_and_build_[a-z0-9_]*\s*\(\s*cline\b", re.IGNORECASE
)
ROUTER_CALL = re.compile(r"\bcall\s+(exec_[a-z0-9_]+_commander)\s*\(", re.IGNORECASE)
EXECUTE_CALL = re.compile(r"\bcall\s+([a-z_][a-z0-9_]*)%execute\s*\(\s*cline\s*\)", re.IGNORECASE)
CASE_STMT = re.compile(r"^\s*case\s*\((.+)\)\s*$", re.IGNORECASE)
SUBROUTINE_START = re.compile(r"^\s*(?:module\s+)?subroutine\s+([a-z_][a-z0-9_]*)\b", re.IGNORECASE)
SUBROUTINE_END = re.compile(r"^\s*end\s*(?:module\s+)?subroutine\s*([a-z_][a-z0-9_]*)?", re.IGNORECASE)
CLINE_HELPER_CALL = re.compile(r"^\s*call\s+([a-z_][a-z0-9_]*)\s*\(.*\bcline\b", re.IGNORECASE)

# These helpers are deliberately small command-default bundles.  Their direct
# missing-key assignments are part of the caller's display-default evidence.
TRACED_DEFAULT_HELPERS = {"set_automask2d_defaults"}

# These helpers inspect or report the command line but do not assign defaults.
READ_ONLY_CLINE_HELPERS = {"warn_for_forced_bootstrap_overrides"}

# Stream project creation only supplies internal project identity and scheduler
# metadata when a caller did not provide a project file.
INTERNAL_CLINE_HELPERS = {
    "create_stream_project": "stream-generated project metadata is not a GUI default",
    "run_make_cavgs_workflow": "runtime-selected strategy derives defaults from project and execution mode",
}

# These defaults are deliberate execution contracts, not GUI choices.  Keep
# them in the audit as explicit informational evidence rather than creating
# descriptors that could override project- or workflow-derived state.
INTENTIONAL_NON_UI_DEFAULTS = {
    ("simple_exec", "abinitio3D_cavgs", "imgkind"):
        "class-average routing is internal to abinitio3D_cavgs",
    ("simple_exec", "abinitio3D_cavgs", "noise_norm"):
        "noise normalization is exposed only by the normalize workflow",
    ("simple_exec", "ctfops", "box"):
        "image geometry is inferred from supplied data; the no-stack fallback is internal setup",
    ("simple_exec", "reconstruct3D_pcg", "nspace"):
        "the PCG operator ignores the setup grid; this is an internal even placeholder",
    ("simple_exec", "volcluster", "box"):
        "box is derived from the first input volume",
    ("simple_stream", "abinitio2D_stream", "projfile_optics"):
        "the optics project filename is a stream-generated artifact",
}


@dataclass(frozen=True)
class Statement:
    text: str
    line: int


@dataclass(frozen=True)
class Procedure:
    name: str
    path: Path
    start_line: int
    statements: tuple[Statement, ...]


def strip_comment(line: str) -> str:
    """Remove a Fortran comment without treating an exclamation in a string as one."""
    quote = ""
    index = 0
    while index < len(line):
        char = line[index]
        if quote:
            if char == quote:
                if index + 1 < len(line) and line[index + 1] == quote:
                    index += 2
                    continue
                quote = ""
        elif char in "'\"":
            quote = char
        elif char == "!":
            return line[:index]
        index += 1
    return line


def logical_lines(path: Path) -> list[Statement]:
    """Return continued, comment-free Fortran statements with their first line."""
    result: list[Statement] = []
    pending = ""
    pending_line = 0
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        text = strip_comment(raw).rstrip()
        if not text.strip():
            continue
        if not pending:
            pending_line = number
        text = text.lstrip()
        if pending and text.startswith("&"):
            text = text[1:].lstrip()
        continued = text.endswith("&")
        if continued:
            text = text[:-1].rstrip()
        pending = f"{pending} {text}".strip()
        if not continued:
            result.append(Statement(pending, pending_line))
            pending = ""
    if pending:
        result.append(Statement(pending, pending_line))
    return result


def source_location(root: Path, path: Path, line: int) -> str:
    return f"{path.relative_to(root)}:{line}"


def is_literal(value: str, constants: dict[str, str]) -> tuple[str | None, str]:
    """Return canonical CLI text and a conservative source classification."""
    value = value.strip()
    value = re.sub(r"\s*!.*$", "", value).strip()
    if re.fullmatch(r"'([^']|'')*'", value):
        return value[1:-1].replace("''", "'"), "literal"
    if re.fullmatch(r'"([^"]|"")*"', value):
        return value[1:-1].replace('""', '"'), "literal"
    if re.fullmatch(r"\.(?:true|false)\.", value, re.IGNORECASE):
        return value.lower(), "literal"
    if NUMERIC_LITERAL.fullmatch(value):
        return value, "literal"
    constant = constants.get(value.lower())
    if constant is not None:
        return constant, "named_constant"
    return None, "expression"


def declaration_binding_kind(declaration: str) -> str:
    """Map a parameter-component declaration to its command-line value kind."""
    declaration = declaration.strip().lower()
    if declaration.startswith("integer"):
        return "int"
    if declaration.startswith("real"):
        return "real"
    if declaration.startswith("logical"):
        return "logical"
    return "char"


def extract_named_constants(paths: Iterable[Path]) -> dict[str, str]:
    constants: dict[str, str] = {}
    pattern = re.compile(r"\bparameter\b.*::\s*([a-z_][a-z0-9_]*)\s*=\s*(.+)$", re.IGNORECASE)
    for path in paths:
        for statement in logical_lines(path):
            match = pattern.search(statement.text)
            if not match:
                continue
            value, kind = is_literal(match.group(2), constants)
            if kind == "literal" and value is not None:
                constants[match.group(1).lower()] = value
    return constants


def extract_parameter_defaults(root: Path, constants: dict[str, str]) -> dict[str, dict]:
    path = root / "src/main/params/simple_parameters.f90"
    defaults: dict[str, dict] = {}
    in_parameters = False
    declaration = re.compile(
        r"^\s*(character(?:\s*\([^)]*\))?|integer(?:\s*\([^)]*\))?|real(?:\s*\([^)]*\))?|"
        r"logical(?:\s*\([^)]*\))?|type\s*\(\s*string\s*\))[^:]*::\s*"
        r"([a-z_][a-z0-9_]*)(?:\s*\([^)]*\))?\s*=\s*(.+)$",
        re.IGNORECASE,
    )
    for statement in logical_lines(path):
        text = statement.text
        if re.match(r"^\s*type\s*::\s*parameters\b", text, re.IGNORECASE):
            in_parameters = True
            continue
        if in_parameters and re.match(r"^\s*contains\b", text, re.IGNORECASE):
            break
        if not in_parameters:
            continue
        match = declaration.match(text)
        if not match:
            continue
        value, value_kind = is_literal(match.group(3), constants)
        if value is None:
            continue
        defaults[match.group(2).lower()] = {
            "value": value,
            "value_kind": value_kind,
            "binding_kind": declaration_binding_kind(match.group(1)),
            "source": source_location(root, path, statement.line),
            "source_kind": "declaration",
        }
    core = root / "src/main/params/simple_parameters_core.f90"
    assignment = re.compile(r"^\s*self%([a-z_][a-z0-9_]*)(?:\s*\([^)]*\))?\s*=\s*(.+)$", re.IGNORECASE)
    in_defaults = False
    for statement in logical_lines(core):
        text = statement.text
        if re.match(r"^\s*(?:module\s+)?subroutine\s+init_dynamic_defaults\b", text, re.IGNORECASE):
            in_defaults = True
            continue
        if in_defaults and re.match(r"^\s*end\s*(?:module\s+)?subroutine\s+init_dynamic_defaults\b", text, re.IGNORECASE):
            break
        if not in_defaults:
            continue
        match = assignment.match(text)
        if not match:
            continue
        value, value_kind = is_literal(match.group(2), constants)
        if value is None:
            continue
        component = match.group(1).lower()
        previous = defaults.get(component, {})
        defaults[component] = {
            "value": value,
            "value_kind": value_kind,
            "binding_kind": previous.get("binding_kind", "char"),
            "source": source_location(root, core, statement.line),
            "source_kind": "dynamic_initialization",
        }
    return defaults


def extract_parameter_bindings(root: Path) -> dict[str, dict]:
    path = root / "src/main/params/simple_parameters_parse.f90"
    bindings: dict[str, dict] = {}
    pattern = re.compile(
        r"\bcall\s+reg%add_(char|int|real|file|dir)\s*\(\s*'([^']+)'\s*,\s*self%"
        r"([a-z_][a-z0-9_]*)(?:\s*\([^)]*\))?",
        re.IGNORECASE,
    )
    for statement in logical_lines(path):
        match = pattern.search(statement.text)
        if not match:
            continue
        key = match.group(2).lower()
        component = match.group(3).lower()
        bindings[key] = {
            "component": component,
            "binding_kind": match.group(1).lower(),
            "source": source_location(root, path, statement.line),
        }
    return bindings


def extract_procedures(paths: Iterable[Path]) -> dict[str, list[Procedure]]:
    """Return every procedure, retaining same-named procedures in distinct modules."""
    procedures: dict[str, list[Procedure]] = {}
    for path in paths:
        statements = logical_lines(path)
        start_index = 0
        while start_index < len(statements):
            start = SUBROUTINE_START.match(statements[start_index].text)
            if not start:
                start_index += 1
                continue
            name = start.group(1).lower()
            end_index = start_index + 1
            while end_index < len(statements):
                end = SUBROUTINE_END.match(statements[end_index].text)
                if end and (end.group(1) is None or end.group(1).lower() == name):
                    break
                end_index += 1
            if end_index == len(statements):
                start_index += 1
                continue
            procedures.setdefault(name, []).append(
                Procedure(name, path, statements[start_index].line, tuple(statements[start_index + 1:end_index]))
            )
            start_index = end_index + 1
    return procedures


def extract_commander_methods(commander_paths: Iterable[Path]) -> dict[str, str]:
    methods: dict[str, str] = {}
    type_start = re.compile(r"^\s*type(?:\s*,[^:]*)?\s*::\s*((?:commander|stream)_[a-z0-9_]+)\b", re.IGNORECASE)
    execute = re.compile(r"\bprocedure\b.*::\s*execute\s*=>\s*([a-z_][a-z0-9_]*)", re.IGNORECASE)
    for path in commander_paths:
        statements = logical_lines(path)
        current_type: str | None = None
        for statement in statements:
            start = type_start.match(statement.text)
            if start:
                current_type = start.group(1).lower()
                continue
            if current_type and re.match(r"^\s*end\s*type\b", statement.text, re.IGNORECASE):
                current_type = None
                continue
            if current_type:
                binding = execute.search(statement.text)
                if binding:
                    methods[current_type] = binding.group(1).lower()
    return methods


def parse_ui_payload(raw: str) -> dict:
    first = raw.find("{")
    last = raw.rfind("}")
    if first < 0 or last < first:
        raise ValueError("simple_private_exec did not print JSON")
    return json.loads(raw[first:last + 1])


def load_ui_registry(executable: Path) -> dict[tuple[str, str], set[str]]:
    completed = subprocess.run(
        [str(executable), "prg=print_ui_json"], text=True, capture_output=True, check=False
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip() or completed.stdout.strip() or "print_ui_json failed")
    payload = parse_ui_payload(completed.stdout)
    registry: dict[tuple[str, str], set[str]] = {}
    for program_name, record in payload.items():
        program = record.get("program", {})
        executable_name = program.get("executable")
        if not executable_name:
            continue
        keys: set[str] = set()
        for section in SECTION_NAMES:
            for entry in record.get(section, []):
                key = entry.get("key")
                if isinstance(key, str):
                    keys.add(key.lower())
        registry[(executable_name, program_name)] = keys
    return registry


def extract_entry_rules(root: Path, executable: str, constants: dict[str, str]) -> tuple[list[str], list[dict]]:
    path = root / "production" / f"{executable}.f90"
    statements = logical_lines(path)
    routers: list[str] = []
    rules: list[dict] = []
    first_router = len(statements)
    for index, statement in enumerate(statements):
        if ROUTER_CALL.search(statement.text):
            first_router = index
            break
    for statement in statements[:first_router]:
        rule = parse_set_rule(root, path, statement, "executable", constants)
        if rule:
            rules.append(rule)
    for statement in statements[first_router:]:
        match = ROUTER_CALL.search(statement.text)
        if match:
            routers.append(match.group(1).lower())
    return routers, rules


def parse_set_rule(root: Path, path: Path, statement: Statement, phase: str, constants: dict[str, str]) -> dict | None:
    missing = MISSING_KEY.match(statement.text)
    if missing:
        if missing.group(1).lower() != missing.group(2).lower():
            return {
                "key": missing.group(2).lower(), "value": None, "value_kind": "expression", "rule": "unresolved",
                "phase": phase, "source": source_location(root, path, statement.line),
                "note": "missing-key condition and assigned key differ",
            }
        return make_rule(root, path, statement.line, missing.group(1), missing.group(3), "missing_key", phase, constants)
    forced = SET_CALL.match(statement.text)
    if forced:
        return make_rule(root, path, statement.line, forced.group(1), forced.group(2), "forced", phase, constants)
    return None


def make_rule(root: Path, path: Path, line: int, key: str, raw_value: str, rule: str, phase: str, constants: dict[str, str]) -> dict:
    value, value_kind = is_literal(raw_value, constants)
    if value is None and rule == "forced":
        # An unconditional computed assignment is command orchestration, not a
        # default presented to a user.  Keep computed missing-key assignments
        # unresolved: those can still require a context-dependent UI default.
        rule = "derived"
    return {
        "key": key.lower(), "value": value, "value_kind": value_kind,
        "rule": rule if value is not None or rule == "derived" else "unresolved", "phase": phase,
        "source": source_location(root, path, line),
        "note": "" if value is not None else f"unsupported assigned value: {raw_value.strip()}",
    }


def extract_router_routes(root: Path, procedures: dict[str, list[Procedure]], executable: str,
                          routers: Iterable[str], constants: dict[str, str]) -> dict[str, list[dict]]:
    routes: dict[str, list[dict]] = {}
    for router_name in routers:
        for router in procedures.get(router_name, []):
            if not router.path.stem.startswith(f"{executable}_"):
                continue
            variable_types: dict[str, str] = {}
            for statement in logical_lines(router.path):
                match = re.match(r"^\s*type\s*\(\s*((?:commander|stream)_[a-z0-9_]+)\s*\)\s*::\s*([a-z_][a-z0-9_]*)", statement.text, re.IGNORECASE)
                if match:
                    variable_types[match.group(2).lower()] = match.group(1).lower()
            current_cases: list[str] = []
            current_statements: list[Statement] = []
            for statement in router.statements:
                case = CASE_STMT.match(statement.text)
                if case:
                    if current_cases:
                        add_router_case(root, router, current_cases, current_statements, variable_types, routes, constants)
                    current_cases = re.findall(r"'([^']+)'", case.group(1))
                    current_statements = []
                    continue
                if current_cases:
                    current_statements.append(statement)
            if current_cases:
                add_router_case(root, router, current_cases, current_statements, variable_types, routes, constants)
    return routes


def extract_stream_routes(root: Path) -> dict[str, list[dict]]:
    """Extract direct simple_stream program-to-stage execute routes."""
    path = root / "production" / "simple_stream.f90"
    statements = logical_lines(path)
    variable_types: dict[str, str] = {}
    for statement in statements:
        match = re.match(r"^\s*type\s*\(\s*(stream_[a-z0-9_]+)\s*\)\s*::\s*([a-z_][a-z0-9_]*)", statement.text, re.IGNORECASE)
        if match:
            variable_types[match.group(2).lower()] = match.group(1).lower()
    routes: dict[str, list[dict]] = {}
    current_cases: list[str] = []
    current_statements: list[Statement] = []

    def add_case() -> None:
        execute = next((EXECUTE_CALL.search(statement.text) for statement in current_statements
                        if EXECUTE_CALL.search(statement.text)), None)
        if execute is None:
            return
        commander_type = variable_types.get(execute.group(1).lower())
        for name in current_cases:
            routes.setdefault(name, []).append({
                "router": "simple_stream",
                "commander_type": commander_type,
                "source": source_location(root, path, current_statements[0].line),
                "rules": [],
            })

    for statement in statements:
        case = CASE_STMT.match(statement.text)
        if case:
            if current_cases:
                add_case()
            current_cases = re.findall(r"'([^']+)'", case.group(1))
            current_statements = []
        elif current_cases:
            current_statements.append(statement)
    if current_cases:
        add_case()
    return routes


def add_router_case(root: Path, router: Procedure, names: list[str], statements: list[Statement], variable_types: dict[str, str], routes: dict[str, list[dict]], constants: dict[str, str]) -> None:
    execute = next((EXECUTE_CALL.search(statement.text) for statement in statements if EXECUTE_CALL.search(statement.text)), None)
    forced = [rule for statement in statements if (rule := parse_set_rule(root, router.path, statement, "router", constants))]
    target = None
    if execute:
        target = variable_types.get(execute.group(1).lower())
    for name in names:
        routes.setdefault(name, []).append({
            "router": router.name,
            "commander_type": target,
            "source": source_location(root, router.path, router.start_line),
            "rules": forced,
        })


def commander_rules(root: Path, procedure: Procedure, constants: dict[str, str],
                    helper_procedures: dict[str, list[Procedure]]) -> list[dict]:
    rules: list[dict] = []
    before_constructor = True
    forced_keys: set[str] = set()
    for statement in procedure.statements:
        if NEW_CALL.search(statement.text) or PARAMS_INITIALIZER_CALL.search(statement.text):
            before_constructor = False
            continue
        direct = parse_set_rule(root, procedure.path, statement, "commander" if before_constructor else "runtime", constants)
        if direct:
            if not before_constructor:
                direct["rule"] = "runtime"
            elif direct["value_kind"] == "expression" and direct["rule"] != "derived":
                direct["rule"] = "unresolved"
            elif direct["rule"] == "forced":
                forced_keys.add(direct["key"])
            elif direct["rule"] == "missing_key" and direct["key"] in forced_keys:
                # A direct forced assignment above this fallback makes the
                # missing-key branch unreachable on this command path.
                direct["rule"] = "shadowed"
            rules.append(direct)
            continue
        helper_call = CLINE_HELPER_CALL.match(statement.text)
        if before_constructor and helper_call:
            helper_name = helper_call.group(1).lower()
            if helper_name in READ_ONLY_CLINE_HELPERS:
                continue
            if helper_name in INTERNAL_CLINE_HELPERS:
                rules.append({
                    "key": None, "value": None, "value_kind": "internal", "rule": "internal_helper",
                    "phase": "commander", "source": source_location(root, procedure.path, statement.line),
                    "note": INTERNAL_CLINE_HELPERS[helper_name],
                })
                continue
            helper_candidates = helper_procedures.get(helper_name, [])
            if helper_name in TRACED_DEFAULT_HELPERS and len(helper_candidates) == 1:
                for helper_statement in helper_candidates[0].statements:
                    helper_rule = parse_set_rule(root, helper_candidates[0].path, helper_statement,
                                                 "helper", constants)
                    if helper_rule:
                        rules.append(helper_rule)
                continue
        if before_constructor and re.search(r"\bcall\s+[a-z_][a-z0-9_]*\s*\(.*\bcline\b", statement.text, re.IGNORECASE):
            rules.append({
                "key": None, "value": None, "value_kind": "expression", "rule": "unresolved",
                "phase": "commander", "source": source_location(root, procedure.path, statement.line),
                "note": "pre-construction helper call requires manual audit",
            })
    return rules


def resolve_execute_method(procedures: dict[str, list[Procedure]], method_name: str | None) -> Procedure | None:
    """Return a uniquely identifiable commander execute method."""
    candidates = procedures.get(method_name or "", [])
    return candidates[0] if len(candidates) == 1 else None


def make_records(root: Path, registry: dict[tuple[str, str], set[str]], bindings: dict[str, dict], baselines: dict[str, dict],
                 procedures: dict[str, list[Procedure]], helper_procedures: dict[str, list[Procedure]],
                 commander_methods: dict[str, str], constants: dict[str, str]) -> tuple[list[dict], list[dict]]:
    records: list[dict] = []
    diagnostics: list[dict] = []
    ui_keys_anywhere = set().union(*registry.values()) if registry else set()
    entry_cache: dict[str, tuple[list[str], list[dict]]] = {}
    router_cache: dict[str, dict[str, list[dict]]] = {}
    for (executable, program), ui_keys in sorted(registry.items()):
        if executable not in {"simple_exec", "single_exec", "simple_stream"}:
            diagnostics.append({"severity": "info", "program": program, "message": f"{executable} is outside the command-default audit scope"})
            continue
        if executable == "simple_stream":
            if executable not in router_cache:
                router_cache[executable] = extract_stream_routes(root)
            candidates = router_cache[executable].get(program, [])
            entry_rules = []
        else:
            if executable not in entry_cache:
                entry_cache[executable] = extract_entry_rules(root, executable, constants)
            routers, entry_rules = entry_cache[executable]
            if executable not in router_cache:
                router_cache[executable] = extract_router_routes(root, procedures, executable, routers, constants)
            candidates = router_cache[executable].get(program, [])
        if len(candidates) != 1:
            diagnostics.append({
                "severity": "error", "program": f"{executable}:{program}",
                "message": "no unique executable-router-command route" if not candidates else "multiple executable-router-command routes",
            })
            continue
        route = candidates[0]
        method_name = commander_methods.get(route["commander_type"] or "")
        execute_method = resolve_execute_method(procedures, method_name)
        if execute_method is None:
            diagnostics.append({
                "severity": "error", "program": f"{executable}:{program}",
                "message": f"could not locate execute method for {route['commander_type'] or 'router case'}",
            })
            continue
        evidence = list(entry_rules) + list(route["rules"]) + commander_rules(
            root, execute_method, constants, helper_procedures
        )
        evidence_by_key: dict[str, list[dict]] = {}
        for item in evidence:
            if item["key"] is None:
                severity = "info" if item["rule"] == "internal_helper" else "warning"
                diagnostics.append({"severity": severity, "program": f"{executable}:{program}", "source": item["source"], "message": item["note"]})
                continue
            evidence_by_key.setdefault(item["key"], []).append(item)
        for key in sorted(ui_keys):
            binding = bindings.get(key)
            baseline = baselines.get(binding["component"]) if binding else None
            key_evidence = evidence_by_key.pop(key, [])
            verified = [item for item in key_evidence if item["rule"] == "missing_key" and item["value"] is not None]
            conflicts = {item["value"] for item in verified}
            if len(conflicts) > 1:
                diagnostics.append({"severity": "error", "program": f"{executable}:{program}", "key": key,
                                    "message": "conflicting statically-known missing-key defaults"})
            override = verified[-1] if verified and len(conflicts) <= 1 else None
            records.append({
                "executable": executable,
                "program": program,
                "key": key,
                "parameter_component": binding["component"] if binding else None,
                "baseline": baseline,
                "commander_override": override,
                "effective_display_default": override["value"] if override else (baseline or {}).get("value"),
                "confidence": "verified" if override or baseline else "unresolved",
                "route": {"router": route["router"], "commander_type": route["commander_type"],
                          "execute": method_name, "source": route["source"]},
                "evidence": key_evidence,
            })
        for key, items in sorted(evidence_by_key.items()):
            rules = {item["rule"] for item in items}
            intentional_reason = INTENTIONAL_NON_UI_DEFAULTS.get((executable, program, key))
            if intentional_reason:
                severity = "info"
                message = f"intentional non-UI default: {intentional_reason}"
            elif "missing_key" in rules:
                if key not in bindings:
                    severity = "error"
                    message = "default rule has no parameter binding"
                elif key in ui_keys_anywhere:
                    severity = "warning"
                    message = "default rule is not exposed by this UI program"
                else:
                    severity = "info"
                    message = "default key is not represented in the UI registry"
            elif "unresolved" in rules:
                severity = "warning"
                message = "unresolved command-line assignment requires manual audit"
            elif rules == {"derived"}:
                severity = "info"
                message = "derived command setting is not a UI display default"
            elif rules == {"shadowed"}:
                severity = "info"
                message = "missing-key fallback is shadowed by a forced command setting"
            elif rules == {"runtime"}:
                severity = "info"
                message = "runtime command-line output is not a UI input"
            else:
                severity = "info"
                message = "forced command setting is not a UI display default"
            diagnostics.append({
                "severity": severity, "program": f"{executable}:{program}", "key": key,
                "source": items[0]["source"],
                "message": message,
            })
    return records, diagnostics


def write_markdown(path: Path, records: list[dict], diagnostics: list[dict]) -> None:
    lines = ["# SIMPLE command-default audit", "", "Generated review evidence; this file is not runtime metadata.", "",
             "## Effective display defaults", "", "| Executable | Program | Key | Baseline | Override | Display default | Confidence |", "| --- | --- | --- | --- | --- | --- | --- |"]
    for record in records:
        baseline = (record["baseline"] or {}).get("value", "")
        override = (record["commander_override"] or {}).get("value", "")
        lines.append("| {executable} | {program} | {key} | {baseline} | {override} | {effective_display_default} | {confidence} |".format(
            executable=record["executable"], program=record["program"], key=record["key"], baseline=baseline,
            override=override, effective_display_default=record["effective_display_default"] or "", confidence=record["confidence"]))
    lines += ["", "## Diagnostics", ""]
    if diagnostics:
        lines += ["| Severity | Program | Key | Source | Message |", "| --- | --- | --- | --- | --- |"]
        for item in diagnostics:
            lines.append("| {severity} | {program} | {key} | {source} | {message} |".format(
                severity=item.get("severity", ""), program=item.get("program", ""), key=item.get("key", ""),
                source=item.get("source", ""), message=item.get("message", "")))
    else:
        lines.append("No diagnostics.")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def fortran_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def key_binding_kind(key: str, bindings: dict[str, dict], baselines: dict[str, dict]) -> str:
    """Return the scalar parser kind for a CLI key when source establishes one."""
    binding = bindings.get(key)
    if binding is not None:
        return binding["binding_kind"]
    baseline = baselines.get(key)
    if baseline is not None:
        return baseline.get("binding_kind", "char")
    return "char"


def as_ui_real(value: str, binding_kind: str) -> str | None:
    """Represent a known numeric CLI value in the UI's real-valued storage."""
    text = value.strip()
    if text.lower() == "no":
        # Legacy numeric switches sometimes spell zero as ``no``.  cmdline
        # stores that token as a character argument, so get_rarg/get_iarg
        # expose their initialized numeric value, zero.
        return "0.0"
    if not NUMERIC_LITERAL.fullmatch(text):
        return None
    text = re.sub(r"_[a-z0-9_]+$", "", text, flags=re.IGNORECASE)
    text = re.sub(r"[dD]", "e", text)
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text


def normalize_ui_default(key: str, value: str, bindings: dict[str, dict], baselines: dict[str, dict], source: str) -> str:
    """Return a type-valid display value or fail before emitting invalid Fortran."""
    binding_kind = key_binding_kind(key, bindings, baselines)
    if binding_kind not in ("int", "real"):
        return value
    normalized = as_ui_real(value, binding_kind)
    if normalized is None:
        raise ValueError(
            f"numeric UI default for {key} is not a valid {binding_kind} value: {value} ({source})"
        )
    return normalized


def collect_program_overrides(root: Path, procedures: dict[str, list[Procedure]], commander_methods: dict[str, str],
                              constants: dict[str, str], bindings: dict[str, dict],
                              baselines: dict[str, dict], helper_procedures: dict[str, list[Procedure]]) -> dict[tuple[str, str], dict[str, str]]:
    """Collect only exact missing-key rules; unsupported routes stay absent."""
    result: dict[tuple[str, str], dict[str, str]] = {}
    for executable in ("simple_exec", "single_exec", "simple_stream"):
        if executable == "simple_stream":
            entry_rules = []
            routes = extract_stream_routes(root)
        else:
            routers, entry_rules = extract_entry_rules(root, executable, constants)
            routes = extract_router_routes(root, procedures, executable, routers, constants)
        for program, candidates in routes.items():
            if len(candidates) != 1:
                continue
            route = candidates[0]
            method_name = commander_methods.get(route["commander_type"] or "")
            execute_method = resolve_execute_method(procedures, method_name)
            if execute_method is None:
                continue
            evidence = list(entry_rules) + list(route["rules"]) + commander_rules(
                root, execute_method, constants, helper_procedures
            )
            overrides: dict[str, str] = {}
            for rule in evidence:
                if rule["rule"] == "missing_key" and rule["key"] and rule["value"] is not None:
                    overrides[rule["key"]] = normalize_ui_default(
                        rule["key"], rule["value"], bindings, baselines, rule["source"]
                    )
            if overrides:
                result[(executable, program)] = overrides
    return result


def write_fortran_defaults(path: Path, baselines: dict[str, dict], bindings: dict[str, dict],
                           overrides: dict[tuple[str, str], dict[str, str]]) -> None:
    """Write a read-only lookup module for UI construction only."""
    baseline_by_key: dict[str, str] = {}
    for component, baseline in baselines.items():
        if baseline["value"] != "":
            baseline_by_key[component] = normalize_ui_default(
                component, baseline["value"], bindings, baselines, baseline["source"]
            )
    for key, binding in bindings.items():
        baseline = baselines.get(binding["component"])
        if baseline is not None and baseline["value"] != "":
            baseline_by_key[key] = normalize_ui_default(
                key, baseline["value"], bindings, baselines, baseline["source"]
            )
    lines = [
        "! Generated by scripts/default_audit.py. Do not edit.",
        "! This module is a UI display-default lookup; commanders do not use it.",
        "module simple_ui_default_values",
        "implicit none",
        "private",
        "public :: get_ui_default",
        "contains",
        "",
        "subroutine get_ui_default(executable, program, key, found, value)",
        "    character(len=*), intent(in)  :: executable, program, key",
        "    logical,          intent(out) :: found",
        "    character(len=*), intent(out) :: value",
        "    found = .false.",
        "    value = ''",
        "    select case(trim(key))",
    ]
    for key, value in sorted(baseline_by_key.items()):
        lines += [f"        case({fortran_string(key)})", "            found = .true.", f"            value = {fortran_string(value)}"]
    lines += ["    end select", "    select case(trim(executable))"]
    for executable in ("simple_exec", "single_exec", "simple_stream"):
        programs = [(program, values) for (exe, program), values in overrides.items() if exe == executable]
        if not programs:
            continue
        lines += [f"        case({fortran_string(executable)})", "            select case(trim(program))"]
        for program, values in sorted(programs):
            lines += [f"                case({fortran_string(program)})", "                    select case(trim(key))"]
            for key, value in sorted(values.items()):
                lines += [f"                        case({fortran_string(key)})", "                            found = .true.", f"                            value = {fortran_string(value)}"]
            lines += ["                    end select"]
        lines += ["            end select"]
    lines += ["    end select", "end subroutine get_ui_default", "", "end module simple_ui_default_values", ""]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--ui-executable", type=Path)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--output-markdown", type=Path)
    parser.add_argument("--output-fortran", type=Path)
    parser.add_argument("--strict", action="store_true", help="fail when the audit finds an error diagnostic")
    args = parser.parse_args()
    root = args.source_root.resolve()
    commander_paths = sorted((root / "src/main/commanders").rglob("*.f90"))
    exec_paths = sorted((root / "src/main/exec").rglob("*.f90"))
    stream_paths = sorted((root / "src/main/stream").rglob("*.f90"))
    helper_paths = [root / "src/defs/simple_default_clines.f90"]
    parameter_paths = [root / "src/main/params/simple_parameters.f90", root / "src/main/params/simple_parameters_core.f90"]
    constants = extract_named_constants(parameter_paths + stream_paths)
    baselines = extract_parameter_defaults(root, constants)
    bindings = extract_parameter_bindings(root)
    procedures = extract_procedures(exec_paths + commander_paths + stream_paths)
    helper_procedures = extract_procedures(helper_paths)
    commander_methods = extract_commander_methods(commander_paths + stream_paths)
    if args.output_fortran:
        overrides = collect_program_overrides(root, procedures, commander_methods, constants, bindings, baselines, helper_procedures)
        write_fortran_defaults(args.output_fortran, baselines, bindings, overrides)
        print(f"generated UI defaults: {len(baselines)} baselines, {len(overrides)} program override sets")
    audit_args = (args.ui_executable, args.output_json, args.output_markdown)
    if any(audit_args) and not all(audit_args):
        parser.error("--ui-executable, --output-json, and --output-markdown must be supplied together")
    if all(audit_args):
        registry = load_ui_registry(args.ui_executable.resolve())
        records, diagnostics = make_records(
            root, registry, bindings, baselines, procedures, helper_procedures, commander_methods, constants
        )
        output = {
            "schema_version": 1,
            "scope": "simple_exec, single_exec, and simple_stream programs emitted by the Fortran UI registry",
            "records": records,
            "diagnostics": diagnostics,
        }
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_markdown.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        write_markdown(args.output_markdown, records, diagnostics)
        errors = sum(item.get("severity") == "error" for item in diagnostics)
        print(f"default audit: {len(records)} UI inputs, {len(diagnostics)} diagnostics, {errors} errors")
        if args.strict and errors:
            return 1
    elif args.strict:
        parser.error("--strict requires the audit JSON and Markdown outputs")
    elif not args.output_fortran:
        parser.error("supply --output-fortran or the three audit output arguments")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
