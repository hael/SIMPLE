# 3-D fractional-update balance refactoring

Date: 2026-10-04. Status: planned, not started.

## What changes

The sampling interface becomes `balance=none|class|cavg` and the sampling
`partition` variable disappears. `cavg` replaces `balance=yes` as the default
everywhere it is set today (`solve3D`, the split checkpoint, `refine3D_states`,
`classify3D_refs`); `refine3D`, `refine3D_auto` and external-reference pose
initialization use `none`. `refine3D_states` and `classify3D_refs` fall back
to `none` when the project has no selected class averages.

- `none`: no groups; global lowest-`updatecnt` tiers (current `sample4update_cnt`).
- `class`: one sampling unit per selected 2-D class, equal quota per unit
  (current `balance=yes` without partition).
- `cavg`: the same units, each assigned to one of `nclust` groups from the
  existing CC average-linkage clustering of the aligned class averages. The
  quota is nested: equal over groups, then equal over the units inside a group,
  then lowest `updatecnt` first within a unit. This applies to single- and
  multi-state runs alike. `clust_crit` is not a parameter of this path.

The nested quota exists because several conformational states hide inside one
view and 2-D classification separates them into classes that average linkage
merges first; the group level keeps the view balance, the unit level keeps the
hidden classes on an equal footing.

Nothing in the library may use projection-direction assignments, or anything
else derived from the 3-D maps or poses, to decide which particles are
sampled. The maps carry bias. `make_projdir_class_samples` and
`prepare_refine3D_states_class_sampling` (`refine3D_states`), the unused
`oris%get_proj_sample_stats`, the `cls2D.cluster` branch of
`classify3D_refs`'s sampling setup, and the sampling use of `partition` and
`clust_crit` are deleted. The `proj` field, `set_projs` and `proj2class` stay
for search, convergence and 2-D averaging.

The `selection` program also reads `balance` (balanced particle selection
across classes). Under the enum, its balanced mode is `balance=class`, off is
`none`, and `cavg` is rejected there. Its own `partition` option (states as
partitions of classes, written to `cls2D.cluster`) is not sampling and stays.

Only active particles (state > 0) are sampled; they always belong to selected
classes, so nothing else needs handling.

## On disk

`clssmp.bin` is the single sampling file for every mode that needs one.
`class_sample` gains an integer `group` (0 = the unit is its own group),
serialized as a fourth header value before `pinds`/`ccs`. The producer of
`cavg` (`simple_view_partition_sampling`) writes one record per selected
class with its group index instead of interleaving the particles of a group
into one record. `class` writes the same records with `group = 0`. Stage
command lines carry `balance` and, for `cavg`, `nclust`; the stages only read
the file. The manifest records `balance`, `nclust` and `nsample`;
continuation, add-on and split handoff rebuild the file from the current
`ptcl2D.class` and never persist particle indices.

## Coverage planning and report

With equal quotas, the number of iterations that visits every particle once is
`sweep = max over units of ceil(pop/quota)`, not `ceil(1/update_frac)`. The
`refine3D_states` `prob_state` init phase (which must label every active
particle before `prob_neigh`), the automatic `maxits` of `refine3D_states` and
`refine3D_auto` (`target updates × sweep`), and any `solve3D` stage rule that
assumes one sweep per `ceil(1/update_frac)` iterations use `sweep` computed
from the unit table before the first stage.

Under the cohort schedule a new set is drawn once per frequency block, so a
full visit takes `sweep` blocks. Nothing is adjusted automatically: the report
states the blocks needed against the blocks planned, and the existing terminal
missing-update pass labels any particle the march did not reach.

The sampler prints the unit table once before the first stage: per group and
unit, population, quota and expected visits per particle over the planned
iterations (blocks under cohorts), with the minimum and maximum over units and
a warning when the minimum is below one or the maximum exceeds ten times the
target.

## Cohort schedule (`refine3D_states` only)

State labelling and pose refinement converge as fixed-point iterations; in
`refine3D_states` they run on one particle set for a whole frequency block.
`solve3D` and `sticky_class_sampling` are not changed by this refactoring;
`sticky_class_sampling` is removed from the `refine3D_states` UI, where
cohorts replace it.

A new internal parameter `cohort_sampling=yes|no` (registered in
`simple_parameters`, not in the UI) is set by `refine3D_states` on the command
lines of its `prob_neigh` frequency blocks. No new project field is needed:
the cohort is the latest sampling round. With `cohort_sampling=yes`,
`sample_ptcls4update3D` draws with the grouped allocator at the first
iteration of a stage (`which_iter == startit`), which stamps a new round
marker on the cohort (lowest `updatecnt` first deprioritizes the previous
cohort). At the other iterations a new `oris%sample4update_rescore` returns
exactly the rows whose `sampled` equals the latest round marker — what
`sample4update_reprod` returns — and then advances `sampled` and `updatecnt`
on those rows once, so the next iteration finds the same set under the new
marker. Earlier cohorts carry older markers and cannot be confused with the
current one; `fromto` and `allow_empty` behave as in `sample4update_reprod`.
`prob_tab` and the matcher keep reproducing the round with
`sample4update_reprod`, which stays marker-neutral. The first-iteration
clean-up of `sampled`/`updatecnt` in `prob_align` (`startit == 1`) is
unchanged. The `prob_state` init phase and the terminal missing-update pass do
not use cohorts.

Trailing reconstruction is unchanged: `current *= u/f`, `previous *=
(1−u)·N/M`, current-map coefficient `u` (`u = f = n/N` per state unless
`ufrac_trec`), chain seeded at full mass by `1/f` and reversed for that
iteration's restoration. Holding one cohort for `k` iterations therefore
gives the cohort a cumulative coefficient `1 − (1−u)^k` with earlier cohorts
fading geometrically. The equal quota biases the composition of the selected
partial reconstruction (about `min(pop, quota)` particles per unit per
iteration); the state-level blend computes no per-unit weight. The policy
states this and that the composition moves with `nsample`, `nclust` and the
class selection.

## Convergence

The 3-D convergence statistics are computed over the current round as today;
under the cohort schedule that is the rescored cohort. The multi-state block
additionally prints, per state, the share of active particles with
`updatecnt == 1`; this diagnostic is not part of the convergence mask.

## Phases

The run uses an out-of-tree Release build for workflow runs and the existing
test executables for unit tests. The beta-galactosidase project after 2-D
clustering is `/data2/bgal_addon/control/bgal_addon.simple` (5,513 particles,
`mskdiam=180 pgrp=d2`); every run starts from a fresh copy. `nsample=1000` is
required in `solve3D`, otherwise sampling is full and balancing never engages.

### Phase 0: baseline (no source change)

On a build of the base state:

- A: `solve3D nsample=1000` with the current default (`balance=yes
  partition=no`), twice.
- B: the same with `partition=yes`, twice. Record the number of selected
  classes; if it is not above 20, use `nclust=8` in B and in every later
  `cavg` run.
- C: `refine3D_states` on the final project of one A run, `nstates=3
  nsample=1500 lpstart=12 lpstop=6 pgrp=d2 mskdiam=180`, once. Record wall
  time and, from the final project, the per-state share of `updatecnt == 1`.

Record exact command lines, final FSC 0.143 resolutions and wall times in
Progress. Keep the final project of one A run under `scratch/keep/`.

Exit: all five runs complete and are recorded.

### Phase 1: implementation, tests, documents

Implement the sections above. Tests:

- `simple_oris_tester`: `class` gives every unit the same count capped at its
  population, lowest `updatecnt` first; `cavg` gives every group the same
  count capped at its population and, inside a group, every unit the same
  count capped at its population; a fixed synthetic skewed population reaches
  every particle in `sweep` iterations; `sample4update_rescore` returns the
  identical subset on repeated calls, advances `sampled` and `updatecnt`
  exactly once per call, returns empty under `allow_empty` on a range without
  cohort rows, and a following draw returns no row of the previous cohort
  while lower-`updatecnt` rows remain.
- `simple_class_sample_io_tester`: `group` round-trips.
- `simple_accum_blend_tester` (`trailing_reconstruction_blend` sub-suite): a
  fixed synthetic cohort reconstructed for `k` iterations matches the closed
  form of the chain (even, odd, rho).
- `rg` gates: no sampling path reads `partition` or `clust_crit`;
  `get_proj_sample_stats`, `make_projdir_class_samples` and
  `prepare_refine3D_states_class_sampling` no longer exist.

Update `doc/policies/importance_sampling_fractional_update_policy.md`,
`doc/policies/3D/solve3D_policy.md`,
`doc/policies/heterogeneity/refine3D_states_policy.md`,
`doc/policies/3D/classify3D_refs_policy.md`, the UI help, and the skills
`simple-solve3d-importance-sampling` and `simple-frac-update-trailing` with
their references.

Exit: Release and Debug builds without new warnings; the fast test gate and
the unit suites of the touched areas pass; `scripts/check_test_registry.py`
passes; the `rg` gates hold.

### Phase 2: comparison runs and report

- A′: `solve3D nsample=1000 balance=class`, twice.
- B′: `solve3D nsample=1000 balance=cavg` (with the Phase 0 `nclust`), twice.
- C′: the Phase 0 C command line on a copy of the same starting project.

Acceptance: A′ and B′ final FSC 0.143 resolutions are no worse than the mean
of A and B respectively by more than the larger of the Phase 0 spread and 3 %
of that mean. C′ completes; its log shows the coverage report; within every
frequency block `% PARTICLES UPDATED SO FAR` stays constant after the block's
first iteration and grows at the next block's first iteration; the per-state
`updatecnt == 1` share is printed.

Then move this plan to `doc/refactoring_notes/completed/` and write the
report beside it.

### File table

| File | Change | Phases |
|---|---|---|
| `doc/refactoring_notes/planned/view_balanced_fractional_sampling_refactoring.md` | plan and progress | each |
| `doc/refactoring_notes/completed/view_balanced_fractional_sampling_refactoring.md` | plan moved, report beside it `view_balanced_fractional_sampling_refactoring_report.md` | 2 |
| `src/defs/simple_type_defs.f90` | `class_sample%group` | 1 |
| `src/fileio/simple_class_sample_io.f90` | serialize `group` | 1 |
| `ori/simple_oris.f90`, `simple_oris_sampling.f90`, `simple_oris_getters.f90`, `simple_oris_reshape.f90` | nested allocator, rescore, deletions | 1 |
| `strategies/search/simple_matcher_smpl_and_lplims.f90`, `simple_view_partition_sampling.f90`, `simple_strategy3D_matcher.f90` | sampling dispatch, `cavg` producer | 1 |
| `commanders/simple/simple_commanders_prob.f90`, `simple_commanders_refine3D.f90`, `simple_commanders_solve3D.f90`, `simple_commanders_project_cls.f90`, `simple_commanders_project_core.f90`, `simple_commanders_mkcavgs.f90` | workflow defaults, deletions, `selection` | 1 |
| `solve/simple_solve3D_controller.f90`, `simple_solve3D_utils.f90`, `simple_solve3D_manifest.f90`, `simple_solve3D_split_checkpoint.f90` | enum, manifest, handoff | 1 |
| `simple_external_reference_pose_initialization.f90`, `simple_convergence.f90`, `sigma2/simple_sigma2_bootstrap.f90`, `apis/simple_core_api.f90` | `none`, diagnostic, key stripping, exports | 1 |
| `params/simple_parameters.f90`, `simple_parameters_parse.f90`, `simple_parameters_phases.f90` | enum, `cohort_sampling` | 1 |
| `ui/simple/simple_ui_solve3D.f90`, `simple_ui_heterogeneity.f90`, `simple_ui_project.f90`, `simple_ui_refine3D.f90`, `ui/simple_ui_params_common.f90` | UI metadata | 1 |
| `src/main/ui/simple_test/*.f90` | test registration | 1 |
| `image/simple_accum_blend_tester.f90`, `src/fileio/simple_class_sample_io_tester.f90` | tests | 1 |
| `doc/policies/importance_sampling_fractional_update_policy.md`, `doc/policies/3D/solve3D_policy.md`, `doc/policies/heterogeneity/refine3D_states_policy.md`, `doc/policies/3D/classify3D_refs_policy.md` | policies | 1 2 |
| `.github/skills/simple-solve3d-importance-sampling/`, `.github/skills/simple-frac-update-trailing/` | skills and references | 1 2 |

## Progress

(none yet)
