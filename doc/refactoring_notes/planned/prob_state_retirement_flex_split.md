# Retire `prob_state`: flex PCA initializes every particle multi-state split

Date: 2026-10-06. Status: planned.

## Why

`refine=prob_state` keeps each particle's previous projection direction and
in-plane angle, scores that pose against the reference of every state (with
shift refinement), and draws the state label by a balanced assignment over
states. It does not work (maintainer, 2026-10-06). It is today the only engine
of four places that turn a single-state refinement into a multi-state one:

1. the split stage of docked `solve3D` on particles (stage 6 by default);
2. the split stage of docked `solve3D_cavgs`;
3. the initialization phase of `refine3D_states` for state-0/1 input with
   `flex=no`;
4. `refine3D_states pose_policy=fixed` (`multivol_mode=input_oris_fixed`), whose
   whole refinement and missing-update pass are `prob_state`.

Projection-aware flex PCA (`flex_pca`, the current default initializer of
`refine3D_states`) becomes the only way to split particles into states.
Class averages carry much more signal per image and keep a random balanced
partition.

## Rulings (2026-10-06)

1. Particle workflows initialize states with `flex_pca`: the docked `solve3D`
   split and state-0/1 input to `refine3D_states`.
2. `refine3D_states pose_policy=fixed` is retired, together with
   `multivol_mode=input_oris_fixed`. Fixed-pose classification is `flex_pca`
   itself.
3. `refine3D_states flex=no` (stochastic initializer plus the `prob_state`
   initialization phase) is retired; the `flex` key goes with it (release 4
   keeps no compatibility keys).
4. Docked `solve3D_cavgs` splits the class averages into `nstates` random
   balanced partitions (equal populations), reconstructs the state volumes, and
   runs its split stage with the ordinary post-split search mode. The
   maintainer's note was cut off after "There is much more signal so n";
   read as "no flex or `prob_state` is needed for class averages".
5. Flex-initialized runs require `nstates >= 3` (docked `solve3D` on particles,
   state-0/1 input to `refine3D_states`); `nstates` is a ceiling, because flex
   merges indistinct states.
6. Executed as an arun run on the Dell; the working tree is committed before
   launch.
7. Flex inside `solve3D` runs at the crop and low-pass of the pre-split
   stage (stage 5 by default) and on the cohort only.
8. When flex returns fewer states than requested, the run continues with the
   flex count, records it in the manifest, and stops with an error below 2.
9. The class-average split gives every state partition the same number of
   class averages (within one).
10. No large-dataset validation in this run; the maintainer tests on big
    datasets afterwards.

## What changes

### `refine3D_states`

- State-0/1 input always runs `flex_pca` (the present default path). Inputs
  with project multi-state labels and docked checkpoints are unchanged.
- `pose_policy` becomes `local|global` (default `global`); `fixed` is rejected
  as an unknown value. `input_oris_fixed` disappears from the code.
- The `prob_state` initialization phase (`l_init_state_assignment`,
  `init_sweep_iters`, the init frequency march and stage) is removed: every
  run starts at the `prob_neigh` frequency march.
- The stochastic state-volume startup for state-0/1 input in
  `initialize_state_volumes` is removed; the branch that uses or reconstructs
  project state volumes for labelled input stays.
- The final missing-update pass uses `refine=greedy` only (the present
  non-fixed path).

### Docked `solve3D` on particles

Today, at the split (`build_solve3D_split_checkpoint`): a cohort-forming
`refine=prob` pass, then random balanced labels (`randomize_states`), then a
post-split-sized subset restricted to the cohort reconstructs the state
volumes, and the split stage runs `prob_state`.

New: the random labels are replaced by `flex_pca` on the cohort with the
pre-split consensus map, at the crop and low-pass of the pre-split stage
(ruling 7). Flex sets the state labels of the cohort and seeds the state
volumes. If flex returns fewer states than requested, the run continues with
the flex count, writes it to the manifest and to the stage command lines, and
stops with an error below 2 (ruling 8). The split stage runs the ordinary
post-split mode (`prob_neigh`, `prob_neigh_mode=geom`). The cohort, the
sticky class sampling after the split, trailing from `TRAILREC_STAGE_SINGLE`,
and the terminal missing-update pass (`refine=greedy`) stay as they are.
`nstates < 3` is rejected for docked particle runs before stage 1.

### Docked `solve3D_cavgs`

`randomize_states` gives every state the same number of class averages
(within one) by a random balanced partition (ruling 9); the split stage emits
the post-split mode instead of `prob_state`. `nstates >= 2` stays allowed.

### `prob_state` removal

- `solve/simple_solve3D_controller.f90`: the `docked_split_stage` branch that
  emits `prob_state`.
- `strategies/parallelization/simple_refine3D_strategy.f90`:
  `l_prob_state_mode` and its branches.
- `strategies/search/simple_strategy3D_matcher.f90`,
  `simple_strategy3D_alloc.f90`: the `prob_state` cases.
- `strategies/search/probabilistic/simple_strategy3D_prob.f90`: the fixed
  projection (`l_fixed_projection`).
- `commanders/simple/simple_commanders_prob.f90`: the `l_state_only` branches.
- `strategies/search/probabilistic/simple_eul_prob_tab.f90`:
  `fill_tab_state_only`, `fill_tab_state_only_range`, `state_assign`, and the
  state-only table storage if nothing else uses it;
  `simple_eul_prob_tab_neigh.f90`: the `prob_state` reference.
- `params/`: `refine` and `pose_policy` choices; `ui/simple/simple_ui_refine3D.f90`,
  `ui/single/single_ui_nano3D.f90`, `ui/simple/simple_ui_heterogeneity.f90`:
  `prob_state`, `fixed`, and `flex` removed from the choices and help.

## Phases

The run uses an out-of-tree Release build for workflow runs and the existing
test executables for unit tests. Data: the beta-galactosidase project after
2D clustering, `/data2/bgal_addon/control/bgal_addon.simple` (5,513
particles, `mskdiam=180 pgrp=d2`), a fresh copy per run. Every `solve3D` run
uses `nsample=2500` and starts from one common `solve3D_cavgs` result with
`cavg_ini_ext=yes`. Beta-galactosidase is homogeneous, so these runs test that
the workflows complete and stay sane, not that they find real states.

### Phase 0: baseline (no source change)

- S: `solve3D_cavgs` once; keep its project under `scratch/keep/` as the
  common start.
- A: `solve3D multivol_mode=docked nstates=3 nsample=2500 cavg_ini_ext=yes`
  from S, twice.
- B: `refine3D_states nstates=3 nsample=1500 lpstart=12 lpstop=6 pgrp=d2
  mskdiam=180` (flex default, `pose_policy=global`) on the final project of a
  single-state `solve3D` from S, once.
- C: `solve3D_cavgs multivol_mode=docked nstates=3` once.

Record command lines, per-state populations, per-state FSC 0.143
resolutions, and wall times in Progress.

### Phase 1: implementation, tests, documents

Implement "What changes". Tests: the balanced partition gives every state the
same number of class averages within one; a docked particle run with
`nstates=2` stops with the documented message; existing unit suites of the
touched areas pass. `rg` gates: no `prob_state`,
`input_oris_fixed`, `fill_tab_state_only`, or `state_assign` in `src/`; no
`flex=no` path; `pose_policy` choices are `local|global`.

Update `doc/policies/heterogeneity/refine3D_states_policy.md`,
`doc/policies/3D/solve3D_policy.md`, `doc/policies/3D/solve3D_cavgs_policy.md`,
`doc/policies/importance_sampling_fractional_update_policy.md`, the
algorithms chapters `refine3d.md`, `solve3d.md`,
`sampling_and_fractional_updates.md`, `heterogeneity_analysis/README.md`,
`heterogeneity_analysis/refine3d_states.md`, the UI help, and the skills
`simple-refine3d` and `simple-solve3d-importance-sampling` with their
references.

Exit: Release and Debug builds without new warnings; the fast test gate and
the unit suites of the touched areas pass; `scripts/check_test_registry.py`
passes; the `rg` gates hold.

### Phase 2: comparison runs and report

Repeat A, B, and C on the new build (A twice). Acceptance: every run
completes; no state falls below the split population floor; in A, the flex
split reports its state count and the run continues with it; the best-state
FSC 0.143 resolution of A and B is no worse than the Phase 0 value by more
than the larger of the Phase 0 spread and 3 %; C completes with every state
populated. Then move this plan to `doc/refactoring_notes/completed/` and
write a self-contained report beside it.

### File table

| File | Change | Phases |
|---|---|---|
| `doc/refactoring_notes/planned/prob_state_retirement_flex_split.md` | plan and progress | each |
| `doc/refactoring_notes/completed/prob_state_retirement_flex_split.md` | plan moved, report beside it | 2 |
| `src/main/commanders/simple/simple_commanders_refine3D.f90` | `refine3D_states`: flex only, no fixed policy, no init phase | 1 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90`, `src/main/solve/simple_solve3D_split_checkpoint.f90`, `src/main/solve/simple_solve3D_utils.f90`, `src/main/solve/simple_solve3D_controller.f90`, `src/main/solve/simple_solve3D_manifest.f90` | flex split, cavgs balanced split, `nstates` floor | 1 |
| `src/main/commanders/simple/simple_commanders_flex_pca.f90`, `src/main/flex/run/*.f90` | only if running flex at the pre-split crop on the cohort needs an entry point or option | 1 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | `l_state_only` removed | 1 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | `l_prob_state_mode` removed | 1 |
| `src/main/strategies/search/simple_strategy3D_matcher.f90`, `simple_strategy3D_alloc.f90`, `simple_strategy3D_srch.f90` | `prob_state` cases removed | 1 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90`, `simple_strategy3D_prob.f90` | state-only table and assignment removed | 1 |
| `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`, `simple_parameters_phases.f90` | choices, `flex` key | 1 |
| `src/main/ui/simple/simple_ui_refine3D.f90`, `simple_ui_heterogeneity.f90`, `src/main/ui/single/single_ui_nano3D.f90` | UI metadata | 1 |
| unit testers of the touched modules, `src/main/ui/simple_test/*.f90` | tests and registration | 1 |
| `doc/policies/heterogeneity/refine3D_states_policy.md`, `doc/policies/3D/solve3D_policy.md`, `doc/policies/3D/solve3D_cavgs_policy.md`, `doc/policies/importance_sampling_fractional_update_policy.md` | policies | 1 2 |
| `doc/algorithms/refine3d.md`, `solve3d.md`, `sampling_and_fractional_updates.md`, `doc/algorithms/heterogeneity_analysis/README.md`, `refine3d_states.md` | algorithms docs | 1 |
| `.github/skills/simple-refine3d/`, `.github/skills/simple-solve3d-importance-sampling/` | skills and references | 1 |

## Progress

Not started.
