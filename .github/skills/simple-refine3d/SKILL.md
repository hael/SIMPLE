---
name: simple-refine3d
description: Use when working on SIMPLE's refine3D workflow, including CLI/UI entrypoints, commanders, probabilistic pre-alignment, 3D search and matcher strategies, online partial reconstruction, volassemble, automasking, nonuniform filtering, multi-state behavior, distributed execution, and related abinitio3D refinement code paths.
---

# SIMPLE Refine3D

This is the single 3D refinement skill. Treat `refine3D` as a layered
fixed-point workflow: commanders choose the route, strategies control iteration
and execution mode, matcher/search modules update particle poses and write
partial reconstructions, and explicit assembly commanders own assembled-reference
work.

## Start Here

1. `doc/policies/refine3D_policy.md`
2. `src/main/ui/simple/simple_ui_refine3D.f90`
3. `src/main/exec/simple_exec_refine3D.f90`
4. `src/main/commanders/simple/simple_commanders_refine3D.f90`
5. `src/main/strategies/parallelization/simple_refine3D_strategy.f90`
6. `src/main/strategies/search/simple_strategy3D_matcher.f90`
7. `src/main/commanders/simple/simple_commanders_rec_distr.f90`
8. `src/main/volume/simple_vol_pproc_policy.f90`

For deeper navigation, read:

- `references/refine3d-code-map.md`
- `references/bayesian-model.md`
- `references/volume-postprocessing.md`

## Ownership Map

- UI/CLI contract: `ui`, `exec`, `simple_cmdline`, `src/main/params/simple_parameters.f90`
- Top-level command ownership: `simple_commanders_refine3D.f90`
- Iteration/execution mode policy: `simple_refine3D_strategy.f90`
- Particle-domain update/search: `simple_strategy3D_matcher.f90` and neighboring `search/` modules
- Probabilistic pre-alignment: `simple_commanders_prob.f90`, `simple_eul_prob_tab*.f90`
- Volume-domain assembly/postprocess: assembly commanders and `simple_vol_pproc_policy`

## Conceptual Boundary

- Particle-domain work: sampling, candidate generation, matching, pose/state updates, partial reconstructions
- Volume-domain work: assembly, even/odd handling, gridding correction, automask, nonuniform filtering, merged outputs

Preserve that split when refactoring.

## Working Rules

- Keep commanders orchestration-heavy, strategies iteration-focused,
  matcher/search modules search-focused, and numerical modules
  implementation-focused.
- Keep shared-memory and distributed refine3D behavior aligned.
- Use `parameters` as parsed workflow state. After `params%new(cline)` or a
  builder initializer, read parsed values from `params`; use `cmdline` mainly
  for pre-parse defaults and child command construction.
- Let project metadata own native sampling. Child command lines should not set
  native `box`/`smpd` or redundant `smpd_crop`; pass `box_crop` only for a
  deliberate crop.
- Keep command-line contracts registered in `src/main/params` before passing new
  keys on SIMPLE command lines.
- Use "probabilistic" precisely. SIMPLE uses probabilistic candidate generation
  and assignment with hard orientation/state updates in key paths; do not call a
  path full soft-assignment EM unless the code actually does that.
- Preserve matcher single-read particle I/O when reconstruction is active.
  Matching/reconstruction should read each particle batch once and update
  partial reconstructions from those images after assignment.
- Preserve artifact handoffs: assignment maps, partition-local partial
  reconstructions, `POLAR_REFS*`, state volumes, even/odd volumes, FSC outputs,
  and state automasks are workflow contracts.
- For `abinitio3D` stage outputs, distinguish planned stage LP from saved
  low-pass diagnostic volumes. `_stageNN_lp.mrc` snapshots use measured state
  FSC resolution when available, with planned stage LP only as fallback.
- For multi-state `abinitio3D_cavgs` with `pgrp_start != pgrp`, keep the
  symmetry-axis search state-local.

## Adjacent Files Worth Reading

- `doc/policies/automasking_policy.md`
- `doc/policies/nonuniform_filtering_policy.md`
- `src/main/volume/simple_reconstructor_eo.f90`
- `src/main/volume/simple_vol_pproc_policy.f90`
- `src/main/nu_filt/*`
- `src/main/pftc/*`
- `src/main/ori/*`
- `src/main/project/*`

## Focused Companion Skills

- For abinitio3D `update_frac`, `nsample*`, `sampled`/`updatecnt`,
  `prob_align`/`prob_tab` sampling reuse, or trailing-reconstruction weighting,
  read `.github/skills/simple-abinitio3d-importance-sampling/SKILL.md`.
- For fractional-update contracts, current partial reconstruction handoffs,
  previous even/odd/rho compatibility, or obsfield mirrors, read
  `.github/skills/simple-frac-update-trailing/SKILL.md`.
- For `filt_mode=nonuniform`, `_nu_filt` reference products, automask/mask
  precedence, or `simple_nu_filter` changes, read
  `.github/skills/simple-nonuniform-regularization/SKILL.md`.

## Common Traps

- `refine3D` file references in docs may use logical names; verify real paths under `strategies/parallelization`, `strategies/search`, and `volume/`.
- Shared-memory and distributed routes should differ in launch/orchestration, not scientific ownership boundaries.
- `volassemble` is intentionally the execution site for expensive volume-domain work.
- Avoid moving assembled-reference postprocessing back into strategy or matcher layers.
- Avoid trading the matcher single-read particle I/O contract for lower peak memory without an explicit policy change.
