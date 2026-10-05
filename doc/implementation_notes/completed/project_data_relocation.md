# Project data relocation

**Status:** completed. `update_project` exposes typed global and per-scope root
mappings, validates the complete relocation before writing, preserves the
ordinary metadata-only route when no mapping is supplied, and writes a new
project rather than modifying the input. The user-facing procedure is in
[`doc/relocating_project_data.md`](../../relocating_project_data.md).

## Contract

Create a new SIMPLE project whose stored dataset paths point from an old
storage root to an existing new storage root.

- `old_root` and `new_root` define a global mapping.
- `mic_old_root` and `mic_new_root` override it for micrograph paths.
- `ptcl_old_root` and `ptcl_new_root` override it for particle-stack paths.
- `cavg_old_root` and `cavg_new_root` override it for class-average paths.
- `vol_old_root` and `vol_new_root` override it for volume paths.
- `projfile_out` optionally names the new project; by default the command writes
  `<input>_remapped.simple` beside the input project.

Each old/new pair is atomic. A scoped pair overrides the global pair for its
scope.

## Path ownership

- Micrographs: `mic.movie`, `mic.intg`, and `mic.boxfile`.
- Particles: `stk.stk` and `stk.boxfile`.
- Class averages: `out.stk`, `out.stkpath`, `out.frcs` for `frc2D`, and
  `out.sigma2`.
- Volumes: `out.vol`, `out.fsc`, and `out.frcs` for `frc3D`.

## Safety requirements

- Match the old root only on a complete path-component boundary.
- Refuse an empty old or new root, a filesystem root as the old root, or equal
  old and new roots.
- Validate the new root and every proposed target before writing.
- Require an explicitly supplied scoped mapping to match at least one path.
- Never modify the input project or overwrite an existing output project.
- Do not write any output project when validation fails.

## Landed implementation

1. The arguments are registered in the typed `parameters` lifecycle and the
   `update_project` UI metadata.
2. The project-domain helper maps supported segment fields and validates all
   targets before the commander publishes a result.
3. The opt-in relocation branch leaves the existing metadata-only branch
   unchanged.
4. Focused project tests cover global and independently scoped mappings, and
   the user guide documents the CLI contract and safety behavior.
