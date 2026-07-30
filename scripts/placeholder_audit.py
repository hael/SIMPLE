#!/usr/bin/env python3
"""Generate an editable placeholder review inventory from a SIMPLE UI registry.

The output is review evidence only.  Approved placeholder wording must be
applied back to the Fortran UI descriptors; this file never feeds the build.

Examples:

  python3 scripts/placeholder_audit.py \
    --ui-executable build/production/simple_private_exec \
    --output-markdown doc/policies/ui_placeholder_review.md

  python3 scripts/placeholder_audit.py \
    --registry-json build/simple_ui_review.json \
    --output-markdown doc/policies/ui_placeholder_review.md
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any


INPUT_SECTIONS = (
    "image input/output",
    "file input/output",
    "parameter input/output",
    "search controls",
    "filter controls",
    "mask controls",
    "computer controls",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--ui-executable", type=Path, help="Built simple_private_exec executable.")
    source.add_argument("--registry-json", type=Path, help="Existing JSON emitted by prg=print_ui_json.")
    parser.add_argument("--output-markdown", type=Path, required=True, help="Review Markdown to write.")
    parser.add_argument("--timeout", type=float, default=30.0, help="Registry query timeout in seconds (default: 30).")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be greater than 0")
    args.output_markdown = args.output_markdown.resolve()
    if args.ui_executable is not None:
        args.ui_executable = args.ui_executable.resolve()
    if args.registry_json is not None:
        args.registry_json = args.registry_json.resolve()
    return args


def load_registry(args: argparse.Namespace) -> tuple[dict[str, Any], str]:
    if args.registry_json is not None:
        try:
            raw = args.registry_json.read_text(encoding="utf-8")
        except OSError as error:
            raise RuntimeError(f"cannot read registry JSON: {error}") from error
        source = str(args.registry_json)
    else:
        if not args.ui_executable.is_file():
            raise RuntimeError(f"UI executable not found: {args.ui_executable}")
        try:
            completed = subprocess.run(
                (str(args.ui_executable), "prg=print_ui_json"),
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=args.timeout,
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            raise RuntimeError(f"print_ui_json timed out after {args.timeout:g} seconds") from error
        if completed.returncode != 0:
            raise RuntimeError(f"print_ui_json failed: {(completed.stdout + completed.stderr).strip()}")
        raw = completed.stdout
        source = f"{args.ui_executable} prg=print_ui_json"
    try:
        registry = json.loads(raw)
    except json.JSONDecodeError as error:
        raise RuntimeError(f"registry is not valid JSON: {error}") from error
    if not isinstance(registry, dict):
        raise RuntimeError("registry root must be a JSON object")
    return registry, source


def markdown_cell(value: Any) -> str:
    text = str(value).replace("|", r"\|").replace("\n", " ").strip()
    return text if text else "*(empty)*"


def build_rows(registry: dict[str, Any]) -> dict[str, list[dict[str, str]]]:
    rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for program_name, program in registry.items():
        for section in INPUT_SECTIONS:
            for input_value in program.get(section, []):
                rows[input_value["key"]].append(
                    {
                        "program": program_name,
                        "section": section,
                        "type": input_value["keytype"],
                        "label": input_value["label"],
                        "placeholder": input_value["placeholder"],
                    }
                )
    return rows


def render_markdown(rows: dict[str, list[dict[str, str]]], source: str) -> str:
    instance_count = sum(len(contexts) for contexts in rows.values())
    lines = [
        "# UI Placeholder Review",
        "",
        "This document inventories placeholder values by actual parameter key and program context.",
        "It is generated review evidence from the registered Fortran UI, not a source of runtime metadata.",
        "After review, apply approved wording to the owning Fortran descriptor and regenerate this inventory.",
        "",
        f"Registry source: `{source}`.",
        "",
        "Each heading is one parameter key. Each row is a separate `ui_program_input` instance, so a key may have context-dependent wording.",
        "Edit **Placeholder** directly. Leave choice, binary, and hidden inputs empty unless the Fortran descriptor rule changes.",
        "",
        f"Inventory: {len(rows)} unique parameter keys and {instance_count} parameter instances.",
    ]
    for key in sorted(rows):
        lines.extend(
            (
                "",
                f"## Parameter `{key}`",
                "",
                "| Used by program | Section | Type | Parameter label | Placeholder |",
                "| --- | --- | --- | --- | --- |",
            )
        )
        for row in sorted(rows[key], key=lambda item: (item["program"], item["section"], item["label"])):
            lines.append(
                "| "
                + " | ".join(
                    markdown_cell(row[column]) for column in ("program", "section", "type", "label", "placeholder")
                )
                + " |"
            )
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    try:
        registry, source = load_registry(args)
    except RuntimeError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    rows = build_rows(registry)
    args.output_markdown.parent.mkdir(parents=True, exist_ok=True)
    args.output_markdown.write_text(render_markdown(rows, source), encoding="utf-8")
    print(f"wrote {args.output_markdown} ({len(rows)} parameter keys)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
