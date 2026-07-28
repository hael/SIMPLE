# Codex Repository References

- `doc/code_overview/code_base_map.md` — generated source tree and module overview.
- `origin-push-workflow.md` — canonical-repository commit and push workflow.
- `windows-msys-install-workflow.md` — Windows/MSYS build and install workflow.
- `windows-conda-run-workflow.md` — Windows/CONDA build and install workflow.
- `Branch-Context-Notes.md` — local branch overview; verify claims against source.
- `branch-context/` — local branch notes when present; verify claims against source.

Regenerate navigation artifacts after relevant source changes:

```text
perl scripts/generate_codeoverview.pl --root . --out doc/code_overview/code_base_map.md
perl scripts/gen_fortran_indexes.pl --root src --out doc/code_overview/fortran-indexes
```

The Fortran index generator writes:

- `doc/code_overview/fortran-indexes/modules.md` — compact module summary.
- `doc/code_overview/fortran-indexes/api_index.md` — procedure/API index.
- `doc/code_overview/fortran-indexes/module_index.md` — human-readable modules, dependencies, and symbols.
- `doc/code_overview/fortran-indexes/symbol_index.csv` — machine-readable symbol index.
- `doc/code_overview/fortran-indexes/module_graph.dot` — module dependency graph for Graphviz-compatible tools.

The generators produce navigational artifacts; source code remains authoritative.
