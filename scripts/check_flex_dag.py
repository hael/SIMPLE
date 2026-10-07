#!/usr/bin/env python3
"""Check FLEX source-unit boundaries, dependency layers and cycles."""

from pathlib import Path
import re
import sys


LAYERS = {
    # L0: contracts and shared definitions.
    "simple_defs_flex": 0,
    "simple_flex_pca_records": 0,
    "simple_flex_pca_stages": 0,
    "simple_flex_pca_artifacts": 0,
    "simple_flex_pca_rounds": 0,
    "simple_flex_pca_run_types": 0,
    # L1: leaves, codecs and project I/O.
    "simple_flex_pca_plane_cache": 1,
    "simple_flex_pca_planes": 1,
    "simple_flex_reconstructor_latent_ops": 1,
    "simple_flex_pca_pcg": 1,
    "simple_flex_pca_polar": 1,
    "simple_flex_pca_crossfsc": 1,
    "simple_flex_pca_deconv": 1,
    "simple_flex_pca_gmm": 1,
    "simple_flex_pca_targets": 1,
    "simple_flex_pca_util": 1,
    "simple_flex_pca_plot": 1,
    "simple_umap": 1,
    "simple_flex_pca_project_gateway": 1,
    "simple_flex_pca_embedding_io": 1,
    "simple_flex_pca_state_parts": 1,
    "simple_flex_weights_state": 1,
    "simple_flex_weights_file": 1,
    # L2: probe-fit state and services.
    "simple_flex_pca_fit_types": 2,
    "simple_flex_pca_mstep": 2,
    "simple_flex_pca_posterior": 2,
    "simple_flex_probe_fit": 2,
    "simple_flex_pca_basis": 2,
    "simple_flex_pca_embed": 2,
    "simple_flex_pca_pairmerge": 2,
    "simple_flex_pca_fit_driver": 2,
    # L3: state inference and reconstruction.
    "simple_flex_pca_states_backend": 3,
    "simple_flex_pca_states_gridding": 3,
    "simple_flex_pca_states_pcg": 3,
    "simple_flex_pca_state_delivery": 3,
    "simple_flex_pca_rec3d": 3,
    "simple_flex_pca_weights": 3,
    "simple_flex_pca_merge": 3,
    # L4-L6: application services, application and strategy.
    "simple_flex_pca_state_service": 4,
    "simple_flex_pca_delivery_3d": 4,
    "simple_flex_pca_application": 5,
    "simple_flex_pca_strategy": 6,
}

IGNORED = {
    "simple_flex_cls_expansion",  # A different program; outside the FLEX-PCA refactor.
}

MODULE_RE = re.compile(r"^\s*module\s+(?!procedure\b|subroutine\b|function\b)(\w+)", re.I)
SUBMODULE_RE = re.compile(r"^\s*submodule\s*\(\s*(\w+)\s*\)", re.I)
USE_RE = re.compile(
    r"^\s*use\b\s*(?:,\s*(?:intrinsic|non_intrinsic)\s*)?(?:::)?\s*([a-z]\w*)\b(.*)$",
    re.I,
)
ONLY_RE = re.compile(r"(?:^|,)\s*only\s*:", re.I)
END_MODULE_RE = re.compile(r"^\s*end\s+(?:sub)?module\b", re.I)


def strip_fortran_comment(line):
    """Remove a free-form Fortran comment while preserving exclamation marks in strings."""
    quote = None
    i = 0
    while i < len(line):
        char = line[i]
        if quote is None:
            if char in ("'", '"'):
                quote = char
            elif char == "!":
                return line[:i]
        elif char == quote:
            if i + 1 < len(line) and line[i + 1] == quote:
                i += 1
            else:
                quote = None
        i += 1
    return line


def fortran_statements(lines):
    """Yield comment-free logical statements with free-form continuations joined."""
    pending = ""
    for line in lines:
        code = strip_fortran_comment(line).strip()
        if not code:
            continue
        if pending and code.startswith("&"):
            code = code[1:].lstrip()
        pending = f"{pending} {code}".strip()
        if pending.endswith("&"):
            pending = pending[:-1].rstrip()
            continue
        yield pending
        pending = ""
    if pending:
        yield pending


def source_units(lines):
    """Return the module or submodule unit declared by one source file."""
    units = []
    owner = None
    body = []
    for line in lines:
        code = line.split("!", 1)[0]
        if owner is None:
            match = SUBMODULE_RE.match(code) or MODULE_RE.match(code)
            if not match:
                continue
            owner = match.group(1).lower()
            body = [line]
        else:
            body.append(line)
        if owner is not None and END_MODULE_RE.match(code):
            units.append((owner, body))
            owner = None
            body = []
    if owner is not None:
        units.append((owner, body))
    return units


def module_imports(lines):
    imports = set()
    for statement in fortran_statements(lines):
        match = USE_RE.match(statement)
        if match:
            imports.add(match.group(1).lower())
    return imports


def broad_core_imports(lines):
    """Return every core-API USE that lacks an explicit ONLY list."""
    broad = []
    for statement in fortran_statements(lines):
        match = USE_RE.match(statement)
        if not match or match.group(1).lower() != "simple_core_module_api":
            continue
        if not ONLY_RE.search(match.group(2)):
            broad.append(statement)
    return broad


def find_cycle(edges):
    state = {}
    stack = []

    def visit(node):
        state[node] = 1
        stack.append(node)
        for dep in sorted(edges.get(node, ())):
            if state.get(dep, 0) == 0:
                cycle = visit(dep)
                if cycle:
                    return cycle
            elif state.get(dep) == 1:
                start = stack.index(dep)
                return stack[start:] + [dep]
        stack.pop()
        state[node] = 2
        return None

    for node in sorted(edges):
        if state.get(node, 0) == 0:
            cycle = visit(node)
            if cycle:
                return cycle
    return None


def main():
    args = [arg for arg in sys.argv[1:] if not arg.startswith("--")]
    verbose = "--verbose" in sys.argv[1:]
    root = Path(args[0] if args else Path(__file__).resolve().parents[1]).resolve()
    paths = sorted((root / "src/main/flex").rglob("*.f90"))
    paths += [
        root / "src/defs/simple_defs_flex.f90",
        root / "src/fileio/simple_flex_weights_file.f90",
        root / "src/main/strategies/parallelization/simple_flex_pca_strategy.f90",
    ]
    problems = []
    edges = {name: set() for name in LAYERS}
    seen = set()
    for path in paths:
        if not path.is_file():
            problems.append(f"missing source: {path.relative_to(root)}")
            continue
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        units = source_units(lines)
        if not units:
            problems.append(f"cannot identify module owner: {path.relative_to(root)}")
            continue
        if len(units) != 1:
            problems.append(
                f"multiple module units in one source: {path.relative_to(root)} ({len(units)})"
            )
        for owner, unit_lines in units:
            if owner in IGNORED or owner.endswith("_tester"):
                continue
            if owner not in LAYERS:
                problems.append(f"unmapped FLEX module {owner}: {path.relative_to(root)}")
                continue
            seen.add(owner)
            if broad_core_imports(unit_lines):
                problems.append(
                    f"broad simple_core_module_api import in {owner}: {path.relative_to(root)}"
                )
            for dep in module_imports(unit_lines):
                if dep in IGNORED:
                    continue
                if dep.startswith("simple_flex") and dep not in LAYERS:
                    problems.append(f"unmapped FLEX dependency {owner} -> {dep}")
                    continue
                if dep not in LAYERS:
                    continue
                edges[owner].add(dep)
                if LAYERS[dep] > LAYERS[owner]:
                    problems.append(
                        f"upward edge L{LAYERS[owner]} -> L{LAYERS[dep]}: {owner} -> {dep}"
                    )
    missing = sorted(set(LAYERS) - seen)
    if missing:
        problems.append("mapped modules without a source: " + ", ".join(missing))
    cycle = find_cycle(edges)
    if cycle:
        problems.append("dependency cycle: " + " -> ".join(cycle))
    if problems:
        print(f"check_flex_dag: {len(problems)} problem(s):")
        for problem in problems:
            print("   " + problem)
        return 1
    if verbose:
        print(
            f"check_flex_dag: {len(seen)} production modules, "
            f"{sum(map(len, edges.values()))} edges: one unit per source, layered and acyclic"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
