# Refine3D Code Map

## Scope

Use this note to find the right SIMPLE module before making changes.

## Top-Level Flow

For normal 3D refinement, the high-level call chain is:

1. `src/main/commanders/simple/simple_commanders_refine3D.f90`
2. `src/main/strategies/parallelization/simple_refine3D_strategy.f90`
3. `src/main/strategies/search/simple_strategy3D_matcher.f90`
4. `src/main/commanders/simple/simple_commanders_rec_distr.f90`

Read them in that order unless the task is tightly scoped.

For public workflow policy, read `doc/policies/refine3D_policy.md` first.

## Ownership By Layer

### Commanders

- `simple_commanders_refine3D.f90`
  Entry point, defaults, automatic workflow, and execution-strategy selection
  through `create_refine3D_strategy`.
- `simple_commanders_prob.f90`
  Probabilistic pre-alignment and distributed table orchestration.
- `simple_commanders_rec_distr.f90`
  Explicit Cartesian and polar assembly pathways.
  `commander_cartesian_volassemble` reduces Cartesian partial reconstructions,
  writes state/even/odd volumes, updates FSC-derived resolution metadata, runs
  automask/nonuniform postprocessing, and refreshes `POLAR_REFS_even.bin` /
  `POLAR_REFS_odd.bin`.
  `commander_polar_volassemble` reduces polar partial sums for
  `polar=yes|obsfield` and writes the `POLAR_REFS.bin` /
  `POLAR_REFS_even.bin` / `POLAR_REFS_odd.bin` handoff set.

### Strategy Layer

- `simple_refine3D_strategy.f90`
  Shared-memory versus distributed master execution, iteration loop, scheduler
  interaction, probabilistic pre-step orchestration, bootstrap polar-reference
  projection, matcher invocation, and assembly-command dispatch.

### Matcher/Search Layer

- `simple_strategy3D_matcher.f90`
  Core per-iteration particle-domain execution, reference preparation from a
  current volume source or `POLAR_REFS*` handoff, batch preparation, strategy
  dispatch, orientation/state/shift updates, and Cartesian partial
  reconstruction or polar partial-sum writing. In online Cartesian
  reconstruction, batch preparation must preserve the already-read raw particle
  images for reconstruction so the matcher does not re-read stacks in a second pass.
- `simple_strategy3D_prob.f90`
  Consumer of the precomputed assignment map.
- `simple_strategy3D_alloc.f90`
  Search-space allocation and probabilistic thresholds.
- `simple_strategy3D_utils.f90`
  Shared 3D search helpers used by concrete strategy implementations.
- `simple_strategy3D_*.f90`
  Concrete search strategies such as greedy, SHC, SNHC, ptree, and neighborhood variants.

### Probability Tables

- `simple_eul_prob_tab.f90`
  Core probabilistic 3D table generation and assignment.
- `simple_eul_prob_tab_neigh.f90`
  Neighborhood-aware extension.

### Reconstruction/Postprocessing

- `simple_matcher_3Drec.f90`
  Partial reconstruction bookkeeping written by the matcher.
- `simple_matcher_pftc_prep.f90`
  `POLAR_REFS*` availability checks and polar central-section loading for
  matcher/probability-table consumers.
- `simple_vol_pproc_policy.f90`
  Plan object for automask and nonuniform-filter mask selection.
- `simple_reconstructor_eo.f90`
  FSC consumer of state masks.

## Current Artifact Contracts

- `ASSIGNMENT.dat` is produced by probabilistic pre-alignment and consumed by `strategy3D_prob`.
- Cartesian matcher work writes partition-local reconstruction updates;
  Cartesian assembly reduces them into state/even/odd volumes and projects
  updated polar central sections for the next iteration.
- Online Cartesian matcher work preserves one particle-stack read per active
  batch. When partial reconstruction is enabled, reconstruction preparation
  consumes the raw images retained during matcher batch construction after the
  particle assignment has been made.
- Polar matcher work writes partition-local polar partial sums; polar assembly
  reduces them into state-major polar references.
- Multi-state `abinitio3D_cavgs` point-group activation is state-local at the
  symmetry-search stage: each current state map supplies its own axis search,
  and the resulting orientation transform applies only to that state.
- `POLAR_REFS.bin`, `POLAR_REFS_even.bin`, and `POLAR_REFS_odd.bin` are file
  handoffs for `prob_tab`, `prob_tab_neigh`, and polar matcher preparation when
  no fresh complete volume source is active. Non-probabilistic current-volume
  matching derives references directly and may materialize these files as cache
  or handoff artifacts.
- Multi-state polar references are state-major: state `s` occupies
  `(s-1)*nspace+1 : s*nspace`.
- If `POLAR_REFS*` are missing or stale but all `vol1..volN` inputs exist,
  `ensure_polar_refs_on_disk` materializes refs through
  `materialize_polar_refs_from_volume_source`; callers must settle the k-range
  with `set_bp_range3D` first.
- If `FRCS_FILE` is absent in the bootstrap pass, matcher preparation creates
  neutral in-memory FRCs rather than using silent zeros.

## Questions To Ask Before Editing

- Is this change about candidate generation, assignment, or final orientation stamping?
- Is this change iteration control or search logic?
- Does this affect both shared-memory and distributed refine3D?
- Does this change touch particle-domain work or assembled-volume work?
- Does a policy doc in `doc/policies/` need to change with the code?
- Does the change preserve the source-vs-handoff `POLAR_REFS*` contract between assembly, matcher, and probability-table code?
- Does online Cartesian reconstruction reuse matcher batch images instead of adding a second image-stack read?
- If multi-state `abinitio3D_cavgs` activates a point group, does the symmetry-axis search remain state-local?
