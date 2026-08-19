# Project data relocation implementation plan

This plan implements the acceptance contract in
[project_data_relocation_spec.md](project_data_relocation_spec.md).

1. Register the relocation arguments in the typed `parameters` lifecycle and
   expose them in the `update_project` UI metadata.
2. Add a project-domain helper that maps supported segment fields, validates
   all proposed targets, and applies changes only after validation succeeds.
3. Add an opt-in relocation branch to the `update_project` commander while
   retaining the existing metadata-only branch unchanged.
4. Add focused tests for global and independently scoped mappings, plus user
   documentation for the CLI contract and safety behavior.
5. Run non-compiling validation: inspect the diff, scan Fortran macro and
   conflict-marker hazards, run formatting checks, and record unrun build and
   runtime tests for handoff.
