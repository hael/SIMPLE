# Project data relocation specification

## Purpose

Extend `update_project` so a user can create a new SIMPLE project whose stored
dataset paths point from an old storage root to an existing new storage root.
The ordinary metadata-only `update_project` behavior must remain unchanged
unless relocation arguments are supplied.

## Inputs

- `old_root` and `new_root` define a global mapping.
- `mic_old_root` and `mic_new_root` override the mapping for micrograph paths.
- `ptcl_old_root` and `ptcl_new_root` override the mapping for particle-stack paths.
- `cavg_old_root` and `cavg_new_root` override the mapping for class-average paths.
- `vol_old_root` and `vol_new_root` override the mapping for volume paths.
- `projfile_out` optionally names the new project. The default is
  `<input>_remapped.simple` beside the input project.

Each old/new pair is atomic: specifying only one member is an error. A scoped
pair overrides the global pair for that scope.

## Path ownership

- Micrographs: `mic.movie`, `mic.intg`, and `mic.boxfile`.
- Particles: `stk.stk`, `stk.stk_den`, and `stk.boxfile`.
- Class averages: `out.stk`, `out.stkpath`, `out.frcs` for `frc2D`, and
  `out.sigma2`.
- Volumes: `out.vol`, `out.fsc`, and `out.frcs` for `frc3D`.

## Safety requirements

- Match the old root only on a complete path-component boundary.
- Refuse an empty old or new root, a filesystem root as the old root, or equal
  old and new roots.
- Validate that the new root and every proposed target exist before writing.
- Require an explicitly supplied scoped mapping to match at least one path.
- Never modify the input project or overwrite an existing output project.
- Do not write any output project when validation fails.

## Acceptance criteria

1. With no relocation pair, `update_project` retains its existing metadata-only
   read/update/write path.
2. A global mapping applies to all four scopes that contain matching paths.
3. Scoped mappings can be used independently and override a global mapping.
4. Supported file and directory fields are rewritten and persisted in a new
   complete project file.
5. Prefix lookalikes such as `/old/data_backup` remain unchanged when the old
   root is `/old/data`.
6. Focused project tests cover global and per-scope remapping.

Implementation steps are recorded in
[project_data_relocation_plan.md](project_data_relocation_plan.md).
