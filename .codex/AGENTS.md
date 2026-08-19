# SIMPLE Codex agent Instructions

## Feature safety rule

Every new feature must be protected by an explicit opt-in key, with the default
set to disabled. Existing behavior must remain unchanged unless the user
enables the feature deliberately. Add validation for the key and test both the
default-off and opt-in paths before release.

## Codex Repository References

- `local-history.md` - mandatory local-only history workflow for every `.codex`
  edit. Agents must verify and snapshot before and after changes, must use the
  guarded helper, and must never reinitialize or rewrite the history database.

- `SPEC-First Workflow Policy for Codex.md` - mandatory SPEC-before-PLAN gate
  for nontrivial features, behavior changes, significant fixes, integrations,
  workflow changes, and multi-component work. Resolve blocking requirement
  questions and finalize acceptance criteria before detailed planning or
  implementation; trivial fully specified work may use a lightweight SPEC.
  Store SPEC and PLAN files in `doc/implementation_notes/` by default and link
  them to each other; use OneDrive or another external vault only when the user
  explicitly requests it.

- `ASD-STE100 Technical Communication and Test Evidence Policy.md` - mandatory
  communication policy for implementations, test plans, result reviews, and
  handoffs. Explain the algorithm and purpose before commands, cite important
  changes with current file and line numbers, define test hypotheses and gates,
  and explain why supplied evidence passes, fails, or needs more testing.

- `msys-ucrt64-rsync-oracle-workflow.md` - controlled Windows-to-Oracle rsync,
  remote line-ending normalization, and debug-compilation procedure. These are
  high-risk operations and require the user's explicit approval for each dry
  run, real transfer, remote rewrite, and compilation action.

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
