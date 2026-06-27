---
name: simple-refine3d
description: Use when working on SIMPLE's refine3D subsystem, including CLI/UI entrypoints, commanders, probabilistic pre-alignment, search and matcher strategies, volume assembly, automasking, nonuniform filtering, distributed execution, and related policy docs.
---

# SIMPLE Refine3D

This is a cross-folder workflow. Read code and policy together.

## Start Here

1. `doc/policies/refine3D_policy.md`
2. `src/main/ui/simple/simple_ui_refine3D.f90`
3. `src/main/exec/simple_exec_refine3D.f90`
4. `src/main/commanders/simple/simple_commanders_refine3D.f90`
5. `src/main/strategies/parallelization/simple_refine3D_strategy.f90`
6. `src/main/strategies/search/simple_strategy3D_matcher.f90`
7. `src/main/commanders/simple/simple_commanders_rec_distr.f90`
8. `src/main/volume/simple_volume_postprocess_policy.f90`

## Ownership Map

- UI/CLI contract: `ui`, `exec`, `simple_cmdline`, `simple_parameters`
- Top-level command ownership: `simple_commanders_refine3D.f90`
- Iteration/execution mode policy: `simple_refine3D_strategy.f90`
- Particle-domain update/search: `simple_strategy3D_matcher.f90` and neighboring `search/` modules
- Probabilistic pre-alignment: `simple_commanders_prob.f90`, `simple_eul_prob_tab*.f90`
- Volume-domain assembly/postprocess: `commander_volassemble` and `volume_postprocess_policy`

## Conceptual Boundary

- Particle-domain work: sampling, candidate generation, matching, pose/state updates, partial reconstructions
- Volume-domain work: assembly, even/odd handling, gridding correction, automask, nonuniform filtering, merged outputs

Preserve that split when refactoring.

## Adjacent Files Worth Reading

- `doc/policies/automasking_policy.md`
- `doc/policies/nonuniform_filtering_policy.md`
- `src/main/volume/simple_reconstructor_eo.f90`
- `src/main/pftc/*`
- `src/main/ori/*`
- `src/main/project/*`

## Common Traps

- `refine3D` file references in docs may use logical names; verify real paths under `strategies/parallelization`, `strategies/search`, and `volume/`.
- Shared-memory and distributed routes should differ in launch/orchestration, not scientific ownership boundaries.
- `volassemble` is intentionally the execution site for expensive volume-domain work.
