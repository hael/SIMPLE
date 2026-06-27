---
name: simple-bayesian-3d-refinement
description: Use when working on SIMPLE's Bayesian-style and probabilistic 3D refinement stack, including refine3D, probabilistic alignment, 3D search strategies, matcher particle I/O, online partial reconstruction, volassemble, automasking, nonuniform filtering, multi-state 3D workflows, or related ab initio and refinement code paths.
---

# SIMPLE Bayesian-Style 3D Refinement

Use this skill as the default mental model for SIMPLE's probabilistic 3D
refinement code. Treat the implementation as a layered fixed-point workflow:
commanders choose the route, strategies control iteration and execution mode,
matcher/search modules update particle poses and write partial Cartesian
reconstructions or polar partial sums, and explicit assembly commanders own
assembled-reference work.

## Quick Start

1. Identify which layer the task belongs to.
   Use [references/refine3d-code-map.md](./references/refine3d-code-map.md) for code navigation.
2. Reconstruct the scientific intent before changing code.
   Use [references/bayesian-model.md](./references/bayesian-model.md) for the probabilistic model and terminology.
3. Treat volume-domain policy as a separate concern from particle-search policy.
   Use [references/volume-postprocessing.md](./references/volume-postprocessing.md) for `volassemble`, automasking, FSC-mask consumption, and nonuniform filtering.
4. For current public workflow policy, read `doc/policies/refine3D_policy.md`.

## Working Rules

- Distinguish clearly between particle-domain work and volume-domain work.
  Particle-domain work lives in probabilistic alignment, search strategies,
  matcher preparation, and partial reconstruction. Volume-domain work lives in
  `volassemble`, reconstruction postprocessing, and FSC-mask consumption.
- Preserve the current SIMPLE architecture. Keep commanders orchestration-heavy,
  strategies iteration-focused, matcher/search modules search-focused, and
  numerical modules implementation-focused.
- Keep shared-memory and distributed refine3D behavior aligned. When changing
  workflow logic, inspect both execution modes and avoid introducing drift
  between them.
- Preserve explicit assembly ownership. `simple_refine3D_strategy.f90`
  dispatches assembly, but Cartesian and polar assembly commanders execute
  assembled-reference work.
- Use `parameters` as the parsed workflow state. After `params%new(cline)` or a
  builder initializer, read parsed values from `params`; use `cmdline` mainly
  for pre-parse defaults and constructing child command lines.
- Let project metadata own native sampling. For SIMPLE commands driven by
  `projfile`, child command lines should not set native `box`/`smpd` or
  redundant `smpd_crop`. Pass `box_crop` only when a deliberate crop is
  requested, and let the parameter parser derive the effective cropped sampling.
- Keep child command lines sparse. Do not restate parameters that already have
  the defaults needed by the child commander; set only values that change
  behavior or preserve the parent execution mode.
- Keep command-line contracts registered. Any new key passed on a SIMPLE command
  line must first be represented in `src/main/params/simple_parameters.f90` and
  registered through the `src/main/params` parsing/checking path.
- Use "Bayesian-style" or "probabilistic" precisely. SIMPLE uses probabilistic
  candidate generation and assignment with hard orientation/state updates in key
  paths. Do not describe the implementation as full soft-assignment EM unless
  the code path actually does that.
- Treat handoff artifacts as part of the workflow contract. Assignment maps,
  partition-local partial reconstructions or polar sums, `POLAR_REFS*`, state
  volumes, even/odd volumes, FSC outputs, and `automask3D_stateNN.mrc` are
  conventions that downstream code depends on.
- Preserve matcher particle I/O performance. Online refine3D
  matching/reconstruction should read each particle batch once, retain the raw
  batch images needed for Cartesian reconstruction, and update partial
  reconstructions from those images after assignment.
- For `abinitio3D` stage outputs, distinguish the planned stage LP from the
  saved low-pass diagnostic volume. `lpinfo(istage)%lp` controls the staged
  search/reference schedule, while `_stageNN_lp.mrc` snapshots are filtered to
  the measured state FSC resolution when available, with the planned stage LP
  only as fallback.
- For multi-state `abinitio3D_cavgs` with `pgrp_start != pgrp`, keep the
  symmetry-axis search state-local. Each state must search its own current map,
  apply the resulting axis only to orientations assigned to that state, and only
  then continue under the target point-group symmetry.

## Decision Workflow

### Explaining or writing policy

Read:

- [references/bayesian-model.md](./references/bayesian-model.md)
- [references/volume-postprocessing.md](./references/volume-postprocessing.md)

Then:

- explain the scientific objective first
- separate public policy from implementation details
- call out the current owner of each responsibility
- name the exact files that implement the policy
- when the question concerns the current top-level workflow, treat
  `doc/policies/refine3D_policy.md` as authoritative

### Modifying probabilistic alignment or search behavior

Read:

- [references/bayesian-model.md](./references/bayesian-model.md)
- [references/refine3d-code-map.md](./references/refine3d-code-map.md)

Then inspect the nearest files among:

- `src/main/commanders/simple/simple_commanders_prob.f90`
- `src/main/simple_eul_prob_tab.f90`
- `src/main/simple_eul_prob_tab_neigh.f90`
- `src/main/strategies/search/simple_strategy3D_matcher.f90`
- `src/main/strategies/search/simple_strategy3D_prob.f90`
- `src/main/strategies/search/simple_strategy3D_*.f90`

Preserve sampled-particle bookkeeping, separation between table generation and
assignment, single-read matcher batch preparation when reconstruction is active,
state/projection/in-plane/shift ownership, and parity between `prob`,
`prob_state`, and `prob_neigh` where intended.

### Modifying reconstruction postprocessing

Read:

- [references/volume-postprocessing.md](./references/volume-postprocessing.md)
- [references/refine3d-code-map.md](./references/refine3d-code-map.md)

Inspect first:

- `src/main/commanders/simple/simple_commanders_rec_distr.f90`
- `src/main/volume/simple_vol_pproc_policy.f90`
- `src/main/volume/simple_reconstructor_eo.f90`
- `doc/policies/refine3D_policy.md`
- `doc/policies/automasking_policy.md`
- `doc/policies/nonuniform_filtering_policy.md`

Preserve Cartesian and polar assembly commanders as the execution point for
assembled-reference work, state-specific mask compatibility rules,
nonuniform-mask precedence, and shared-memory/distributed workflow parity.

## Non-Obvious Project Facts

- `commander_refine3D` is only the entrypoint/default layer. The real execution
  path is strategy-driven through `create_refine3D_strategy`.
- `simple_refine3D_strategy` owns iteration state, shared-memory versus
  distributed orchestration, probabilistic pre-step dispatch, bootstrap
  polar-reference projection, matcher calls, and assembly-command dispatch.
- `prob_align` is a pre-step that samples particles, generates partition-local
  probability tables, aggregates them, and writes a single assignment map
  consumed by `strategy3D_prob`.
- `refine3D_exec` is the core particle-update loop. It prepares reference
  sections, dispatches search strategies, updates orientations/states/shifts,
  and writes Cartesian partial reconstructions or polar partial sums when
  instructed by the strategy.
- `commander_cartesian_volassemble` assembles Cartesian state volumes, refreshes
  `POLAR_REFS_even.bin` / `POLAR_REFS_odd.bin` for subsequent 3D matching, and
  runs volume postprocessing.
- `commander_polar_volassemble` assembles state-major polar references for
  `polar=yes` and `polar=obsfield`, writing the `POLAR_REFS.bin`,
  `POLAR_REFS_even.bin`, and `POLAR_REFS_odd.bin` handoff set.
- `POLAR_REFS*` are file handoffs for probability-table consumers and for polar
  matcher paths without a fresh complete volume source. Non-probabilistic
  current-volume matching derives references directly from `vol1..volN`; when a
  file handoff is needed, `ensure_polar_refs_on_disk` or
  `materialize_polar_refs_from_volume_source` materialize it after callers
  settle the k-range with `set_bp_range3D`.
- In `abinitio3D_cavgs`, multi-state point-group activation cannot use a single
  state as the symmetry-axis source for all states. The symmetry-search stage
  must determine/apply axes independently per state before the next refinement
  stage enforces `pgrp`.
- The newer volume-free ab initio work is related but distinct. It reuses the
  fixed-point/probabilistic philosophy while replacing explicit intermediate
  volumes with polar-Fourier reprojection estimation for ab initio orientation
  recovery.

## Language To Prefer

- "probabilistic pre-alignment followed by hard assignment"
- "particle-domain search versus volume-domain postprocessing"
- "explicit Cartesian/polar assembly pathway"
- "state-major polar reference blocks"
- "state-specific automask artifact"
- "shared-memory and distributed parity"
- "volume-free ab initio path" for the polar common-line manuscript

Avoid:

- conflating `abinitio3D` and `refine3D`
- describing assembly commanders as mere I/O combiners or distributed-only helpers
- moving assembled-reference postprocessing back into strategy or matcher layers
- trading the matcher single-read particle I/O contract for lower peak memory without an explicit policy change
- describing the current refine3D implementation as purely Bayesian without acknowledging the hard-assignment workflow
