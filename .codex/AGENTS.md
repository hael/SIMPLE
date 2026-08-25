# SIMPLE Codex agent instructions

## Scope

This file contains portable repository guidance that is useful in every SIMPLE
checkout. Machine-specific paths, private infrastructure, user communication
preferences, and local operating procedures do not belong here.

When `.codex/AGENTS-local.md` exists, read it after this file. It is ignored by
the main repository and can contain checkout-specific guidance. Do not assume
that it exists in another clone.

## Feature safety rule

Every new feature must be protected by an explicit opt-in key, with the default
set to disabled. Existing behavior must remain unchanged unless the user
enables the feature deliberately. Add validation for the key and test both the
default-off and opt-in paths before release.

## Generated repository navigation

- `doc/code_overview/code_base_map.md` - generated source tree and module
  overview.
- `doc/code_overview/fortran-indexes/modules.md` - compact module summary.
- `doc/code_overview/fortran-indexes/api_index.md` - procedure and API index.
- `doc/code_overview/fortran-indexes/module_index.md` - human-readable modules,
  dependencies, and symbols.
- `doc/code_overview/fortran-indexes/symbol_index.csv` - machine-readable
  symbol index.
- `doc/code_overview/fortran-indexes/module_graph.dot` - module dependency
  graph for Graphviz-compatible tools.

Regenerate navigation artifacts after relevant source changes:

```text
perl scripts/generate_codeoverview.pl --root . --out doc/code_overview/code_base_map.md
perl scripts/gen_fortran_indexes.pl --root src --out doc/code_overview/fortran-indexes
```

Generated navigation is not source authority. Verify important claims against
the current source, production caller, and Git history.
