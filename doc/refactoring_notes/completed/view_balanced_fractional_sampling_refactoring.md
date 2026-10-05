# 3-D fractional-update balance refactoring

Date: 2026-10-04. Status: completed 2026-10-05 (Phases 0–2; Phase 2 finished under the maintainer's ruling of 2026-10-05; Phase 3, a same-seed reproducibility control added by the maintainer, done; see Progress). Report: `view_balanced_fractional_sampling_refactoring_report.md` beside this plan.

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
under the cohort schedule that is the rescored cohort. By the maintainer's
ruling of 2026-10-05 there is no per-state diagnostic of the share of particles
with `updatecnt == 1`: under cohorts a rescore advances `updatecnt` every
iteration, so the share always reads 0. The coverage report and
`% PARTICLES UPDATED SO FAR` show whether every particle is reached.

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
  time.

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
first iteration and grows at the next block's first iteration.

By the maintainer's ruling of 2026-10-05, the equivalence of `balance=class`
with the old `balance=yes` is settled by a fixed-seed paired test rather than by
more unseeded repeats: case A once on the base build (`balance=yes`) and once on
the new build (`balance=class`), in sequence, with the same `SIMPLE_SEED`
inherited by the distributed workers. The sampled particle set is compared at
every iteration of every stage. Identical throughout: `class` is accepted as
equivalent, both final FSC 0.143 resolutions are recorded, and the de novo
divergence rate of the unseeded runs is reported as an open observation. Not
identical: the phase stops at the first stage and iteration that differ.

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
| `src/main/simple_convergence.f90` | remove the per-state `updatecnt == 1` diagnostic (maintainer's ruling of 2026-10-05) | 2 |
| `ori/simple_oris_tester.f90` | tests of the nested quota, the sweep and cohort rescoring (named in Phase 1's tests, missing from this table) | 1 |
| `doc/policies/importance_sampling_fractional_update_policy.md`, `doc/policies/3D/solve3D_policy.md`, `doc/policies/heterogeneity/refine3D_states_policy.md`, `doc/policies/3D/classify3D_refs_policy.md` | policies | 1 2 |
| `.github/skills/simple-solve3d-importance-sampling/`, `.github/skills/simple-frac-update-trailing/` | skills and references | 1 2 |

## Progress

### Phase 0: baseline, 2026-10-04, done

Files changed: this plan only (this Progress entry). No source file changed.

Build. Every Phase 0 number comes from an out-of-tree Release build of the
base state (HEAD 9ac9517d0, clean working tree, no launch differences) in
`/home/elmlundho/agent_runs/balance_sampling/scratch/build_release`, made by
`scratch/build_release.sh`: CMake 3.26.5, GNU Fortran 15.2.1 from
gcc-toolset-15, `CMAKE_BUILD_TYPE=Release`, `BUILD_TESTS=ON`, OpenMP on,
installed into the build tree. The base build prints five linker warnings
("requires executable stack", from `json_value_module` and
`simple_commanders_cavgs`); they are the warning baseline for Phase 1. Every
run sets `SIMPLE_PATH` and `PATH` to this build, because the machine's login
shell points them at another checkout.

Data and settings. Each run starts from a fresh copy of
`/data2/bgal_addon/control/bgal_addon.simple`: beta-galactosidase, 5,513
particles, box 256 at 1.275 Å per pixel, after 2-D classification into 50
classes. All 50 class averages are selected (state 1), so the number of
selected classes is above 20 and the view-partition runs keep the default
`nclust=20`; Phase 2's `cavg` runs use `nclust=20` too (the default, not set
on the command line). Execution uses `nparts=4 nthr=6` (four distributed parts
with six threads each) as in earlier runs on this machine. The cases ran one
after another on an otherwise idle machine (load 9 to 11, which is these runs'
own load), so wall times are comparable. Runs are started by
`scratch/run_case.sh NAME KIND [extra]`, which writes
`scratch/logs/case_NAME.log` with every line time-stamped and the
`/usr/bin/time -v` summary at the end. Per-iteration tables are in
`scratch/logs/table_NAME.md` (made by `scratch/extract_iters.py`). Paths below
are relative to `/home/elmlundho/agent_runs/balance_sampling`.

"Final FSC 0.143" below is the resolution where the Fourier shell correlation
(FSC) between the two independently reconstructed half maps drops to 0.143.
It is taken from the last `POSTPROCESS: FSC ... 0.5/0.143` line of the log,
which belongs to the final reconstruction from all particles at the original
pixel size. The last refinement iteration of every `solve3D` run reports
4.03 Å, the limit of its downscaled box, so it cannot separate the runs.

Exact command lines, each run in its own fresh directory under
`scratch/runs/NAME` holding the copy `bgal.simple`:

- A (current default, class-balanced sampling, `balance=yes partition=no`):
  `simple_exec prg=solve3D projfile=bgal.simple mskdiam=180 pgrp=d2 nsample=1000 nparts=4 nthr=6`
- B (view-partition sampling over 20 groups of similar class averages):
  `simple_exec prg=solve3D projfile=bgal.simple mskdiam=180 pgrp=d2 nsample=1000 nparts=4 nthr=6 partition=yes`
- C (`refine3D_states` on the final project of run A1):
  `simple_exec prg=refine3D_states projfile=bgal.simple nstates=3 nsample=1500 lpstart=12 lpstop=6 pgrp=d2 mskdiam=180 nparts=4 nthr=6`

C's input project. `solve3D` writes its final project and maps in the
execution directory `scratch/runs/A1/1_solve3D`. The top-level `bgal.simple`
is left unchanged. The final project is kept with its maps in
`scratch/keep/caseA_final` (`bgal.simple`, `recvol_state01.mrc` with its even
and odd halves, `rec_final_state01.mrc`, `fsc_state01.bin`,
`sigma2_state.bin`, `RESOLUTION_FINAL_STATE01`, `solve3D_manifest.txt`). The
project's output segment points at the consensus map and FSC file by absolute
path inside `scratch/runs/A1/1_solve3D`. A first C attempt with that directory
removed stopped at once, because flex principal component analysis (flex PCA,
the default state initialization of `refine3D_states`) found no consensus map
(`scratch/logs/case_C_attempt0_hidden_A1.log`). That is a fault in my test
setup, not in SIMPLE. Before C, `run_case.sh` therefore recreates
`scratch/runs/A1/1_solve3D` from `scratch/keep/caseA_final` and copies
`bgal.simple` from there, so C depends only on `scratch/keep/`. Phase 2's C′
must use the same `run_case.sh C refine3D_states` step.

Results:

| Run | Final FSC 0.5 / 0.143 (Å) | Wall time | Peak memory | Log |
|---|---|---|---|---|
| A1 | 4.41 / 3.84 | 28:17 | 3.2 GB | `scratch/logs/case_A1.log` |
| A2 | 4.35 / 3.93 | 28:12 | 3.2 GB | `scratch/logs/case_A2.log` |
| B1 | 4.47 / 3.98 | 28:27 | 3.2 GB | `scratch/logs/case_B1.log` |
| B2 | 4.41 / 3.84 | 28:28 | 3.2 GB | `scratch/logs/case_B2.log` |
| C | state 1: 4.95 / 4.03; state 2: 13.06 / 10.20 | 10:59 | 4.9 GB | `scratch/logs/case_C.log` |

Acceptance limits for Phase 2, from the rule in Phase 2 (no worse than the
mean by more than the larger of the Phase 0 spread and 3 % of the mean):

- A: mean 3.885 Å, spread 0.09 Å, 3 % of the mean 0.117 Å, so A′ ≤ 4.00 Å.
- B: mean 3.91 Å, spread 0.14 Å, 3 % of the mean 0.117 Å, so B′ ≤ 4.05 Å.

All four `solve3D` runs ran eight stages (low-pass limits 23.3, 21.8, 20.4 and
10.2 Å, then about 9.3, 8.4, 7.8 and 4.5 Å, with stages 5 to 7 sometimes moved
finer early when FSC 0.5 allowed it), sampled about 1,000 particles (18 %) per
iteration, and reached 100 % of particles updated in stage 8. B's log prints
the group table (`scratch/logs/B1_groups.txt`): 50 classes in 20 groups of 1 to
4 classes, holding 57 to 701 particles. The largest group holds 26 times as
many particles as the smallest, but every group gets the same 52 particles
per iteration (1.9 times between largest and smallest sampled share after
capping at the population).

C in detail (`scratch/logs/table_C.md`, `scratch/logs/updatecnt_C.md`):

- Flex PCA asked for three states but its population floor removed one. The
  third state's 1,309 particles had fewer than its 2,000 effective particles
  and were set to state 0 (excluded). The run continued with two states:
  `NSTATES FROM PROJECT: 2`, 4,204 active particles. The run did not fail, so
  this command line is the C that Phase 2 repeats. A supplementary run with
  `flex=no` is reported below.
- The `prob_state` initialization phase did not run, because flex PCA had
  already labelled every particle. Five frequency blocks of three iterations
  each (`prob_neigh`, matching low-pass 12.1, 9.1, 7.3, 6.0 and 6.0 Å), 15
  iterations, state-overlap targets reached in blocks 4 and 5 only.
- Sampling today follows projection directions
  (`PROJDIR-BALANCED SAMPLING BINS 2339`) and drew 55.6 % (2,339 of 4,204)
  every iteration, although `nsample=1500`. The update fraction 1500/5513 =
  0.272 was computed before flex PCA dropped the third state.
  `% PARTICLES UPDATED SO FAR` rose every iteration (55.6, 78.3, 88.6, ...)
  and reached 100 % at iteration 11. This is the behavior the cohort schedule
  replaces.
- Final project (`scratch/runs/C/1_refine3D_states/bgal.simple`): share of
  active particles with update count 1 (`updatecnt == 1`) is 0.016 in
  state 1 (3,377 particles, mean update count 7.85) and 0.006 in state 2
  (827 particles, mean 10.39). Update counts carry over from `solve3D`.

Supplementary C with `flex=no` (C_flexno, the C command line plus `flex=no`,
`scratch/logs/case_C_flexno.log`, `table_C_flexno.md`,
`updatecnt_C_flexno.md`). Not an exit item. It shows the three-state path,
including the `prob_state` initialization phase, which the cohort schedule
and the coverage plan also change. It completed in 15:36 (peak memory
5.4 GB) with three states. The initialization phase ran its cap of 10
iterations (state overlap 0.885 against the 0.95 target), followed by five
frequency blocks of three `prob_neigh` iterations, 25 iterations in all.
Projection-direction sampling used 2,771 bins and drew 50.3 % of particles
every iteration. 100 % had been updated by iteration 12. Final FSC 0.5 / 0.143:
state 1 11.66 / 9.89 Å (1,134 particles), states 2 and 3 6.04 / 4.08 Å (2,166
and 2,213 particles). No particle had update count 1 in the final project
(share 0.000 in all three states; mean update counts 14.9, 12.4 and 11.6).
Phase 2 should repeat C as recorded above and may repeat C_flexno for the
three-state comparison.

Exit: all five runs (A1, A2, B1, B2, C) completed with exit status 0 and are
recorded above with command lines, final FSC 0.143 resolutions and wall times.
The final project of A1 is kept in `scratch/keep/caseA_final`. The bulky run
outputs were deleted after their numbers were recorded; the logs remain.

### Phase 1: implementation, tests, documents, 2026-10-04, done

Files changed: `src/defs/simple_type_defs.f90`,
`src/fileio/simple_class_sample_io.f90` and its tester,
`src/main/ori/simple_oris.f90`, `simple_oris_sampling.f90`,
`simple_oris_getters.f90`, `simple_oris_tester.f90`,
`src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`,
`simple_view_partition_sampling.f90`,
`src/main/commanders/simple/simple_commanders_solve3D.f90`,
`simple_commanders_refine3D.f90`, `simple_commanders_project_core.f90`,
`src/main/solve/simple_solve3D_controller.f90`, `simple_solve3D_utils.f90`,
`simple_solve3D_manifest.f90`, `simple_solve3D_split_checkpoint.f90`,
`src/main/simple_external_reference_pose_initialization.f90`,
`src/main/simple_convergence.f90`, `src/main/sigma2/simple_sigma2_bootstrap.f90`,
`src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`,
`simple_parameters_phases.f90`, `src/main/ui/simple/simple_ui_solve3D.f90`,
`simple_ui_heterogeneity.f90`, `simple_ui_project.f90`,
`src/main/image/simple_accum_blend_tester.f90`, the four policies and the two
skills with their references named in this phase, and this plan (one file-table
row added for `simple_oris_tester.f90`). Files in the table that needed no
change: `simple_strategy3D_matcher.f90`, `simple_commanders_prob.f90`,
`simple_commanders_mkcavgs.f90`, `simple_oris_reshape.f90`,
`simple_core_api.f90`, `simple_ui_refine3D.f90`, `simple_ui_params_common.f90`
and the test registration under `src/main/ui/simple_test/` (the new tests live
in existing sub-suites).

Inventory of the base code, as the plan asked (callers and reads before the
change):

- `sample4update_class` (the per-class quota allocator in `simple_oris`):
  the 3D dispatch `sample_ptcls4update3D` (used by `prob_align`,
  `prob_align_neigh` and the non-probabilistic path of the 3D matcher), the
  initial greedy sample of `solve3D` when it starts from class averages, and the
  reconstruction sample of the docked split checkpoint (restricted to the
  sticky cohort).
- `read_class_samples`: the 3D dispatch and the split checkpoint.
  `write_class_samples`: `solve3D`, `refine3D_states` (projection-direction
  bins) and `classify3D_refs`.
- `get_class_sample_stats` (one record per class, particles ranked by 2D
  score): `solve3D`, `classify3D_refs` (also over `cls2D` `cluster` labels),
  and three non-sampling users that stay unchanged: `sample_classes`
  (balanced particle selection, also over `cluster`), `bootstrap_cavgs` and
  `reseed_classes`. `get_proj_sample_stats` had no caller.
- `params%balance` was read by the 3D dispatch, by the sampling setup of
  `refine3D_states` and `classify3D_refs`, and by the `solve3D` stage
  controller. It was set to `yes` by `solve3D` (default), the split checkpoint,
  `refine3D_states` and `classify3D_refs`; to `no` by `refine3D_auto`, the three
  terminal missing-update passes and external-reference pose initialization;
  `solve2D` deletes it for class averaging. The `selection` program lists
  `balance=yes|no` in its UI, but its commander never reads it.
- `params%partition` was read by `solve3D` (view groups; rejection with input
  volumes), by `classify3D_refs` (`cls2D` `cluster` groups), by the parameter
  check that forbade `partition=yes` with `balance=no`, and by the
  `selection` and `sample_classes` programs. The last two are not sampling and
  keep it.
- `params%clust_crit` was read for sampling only by the view-group producer.
  The class-average clustering programs keep their own use of it.
- `sticky_class_sampling` is set internally by the `solve3D` docked split
  handoff to `refine3D_states`. That use stays, with its meaning unchanged.
- No `solve3D` stage rule assumed one sweep per `ceil(1/update_frac)`
  iterations. `refine3D_states` (initialization sweep and `maxits`),
  `refine3D_auto` and `classify3D_refs` (`maxits`) used `ceiling(4N/n)`.

What was implemented (the design sections above, with these decisions where the
plan left a choice open):

- `balance=none|class|cavg` replaces `balance=yes|no` and the sampling use of
  `partition`. It is held in a five-character `parameters` field, defaults to
  `none`, and is validated as an enum. `cohort_sampling=yes|no` is an internal
  parameter, registered in `parameters` but not in the user interface.
- The `class_sample` type gains `group`, written as the fourth header value of
  every record of `clssmp.bin`. Files from before this change are not read.
- The nested quota of `sample4update_class` works at three levels:
  - groups: equal shares, increasing in equal steps as before, capped at each
    group's population; under `class` every unit is its own group, so `class`
    reproduces the old allocation exactly;
  - units inside a group: equal shares capped at population. The remainder of
    the equal split goes to the units whose particles have the lowest mean
    update count. This makes it the same in every distributed partition and
    rotates it between units over iterations;
  - inside a unit: lowest update count first.

  `class_sample_quotas` and `class_sample_sweep` give the expected quotas and
  the draws one full sweep takes. `sample4update_rescore` serves the cohort
  schedule. `get_proj_sample_stats`, `make_projdir_class_samples`, the old
  `prepare_refine3D_states_class_sampling`, the `cls2D` `cluster` branch of
  `classify3D_refs` and every sampling use of `partition` and `clust_crit`
  are deleted.
- `simple_view_partition_sampling` now produces the units for both modes
  (`make_class_samples`). Under `cavg` it clusters the selected class averages
  by average linkage on their aligned correlation into `nclust` groups. It also
  prints the coverage report before the first stage
  (`report_class_sample_coverage`). `refine3D_states` and `classify3D_refs`
  leave rows that are inactive in `ptcl3D` out of their units; `solve3D` does
  not, because it resets its `ptcl3D` states from the 2D selection later.
- `solve3D` defaults to `cavg`. With input volumes it defaults to `class`
  (random classes, no averages to group) and rejects `cavg`, as it rejected
  `partition=yes` before. It writes `clssmp.bin` for every fractional run,
  with class units when `balance=none`, because the initial greedy sample and
  the split checkpoint always drew from class units. The split checkpoint's
  assignment pass uses the run's `balance`, or `class` when the run uses
  `none` (before, it forced `yes`). Child command lines no longer inherit
  `balance` and `nclust` (`strip_sampling_keys`, which replaces
  `strip_view_partition_keys`). The stage controller sets `balance` on every
  stage line, and `nclust` under `cavg`. The run manifest records `balance`
  and `nclust` when given on the command line (replacing `partition`, `nclust`
  and `clust_crit`), plus the effective `nsample` as before. A manifest
  written before this change that recorded `partition` no longer reads, as
  release 4 allows.
- `refine3D_states` and `classify3D_refs` use `cavg`, and fall back to `none`
  without selected class averages. `refine3D_states` computes the sweep from
  the unit table before the first stage. The `prob_state` initialization phase
  runs at least that many iterations, and the automatic `maxits` is 4 × sweep
  (between 10 and 50). The prob_neigh frequency blocks carry
  `cohort_sampling=yes`; the flag is removed after the march, on the
  missing-update pass and in the sigma2 bootstrap. `refine3D_auto` uses `none`,
  and its `maxits` is 4 × `ceil(active/nsample)`. The `classify3D_refs`
  `maxits` rule is not among the rules the plan names and stays as it was.
- The coverage report is printed once per run, before the first stage, by
  `solve3D` (over the planned stage iterations), `refine3D_states` (over the
  planned frequency blocks) and `classify3D_refs` (over the planned
  iterations). It warns when the least-visited unit gets fewer than one visit,
  or the most-visited more than ten times the target.
- The multi-state convergence block prints, per state, the share of active
  particles with update count 1.
- `selection` lists `balance=none|class` and rejects `cavg` on its command
  line. `sticky_class_sampling` is removed from the `refine3D_states` UI.

Tests and checks (logs under `/home/elmlundho/agent_runs/balance_sampling/scratch/logs`):

- New tests:
  - `simple_oris_tester`: `class` gives every unit the same count capped at its
    population, lowest update count first. `cavg` gives every group, then every
    unit of a group, the same count capped at its population, and the group
    remainder goes to the least-updated unit. A fixed skewed population (units
    of 2, 6, 6, 3 and 20 particles in three groups, 19 drawn) has a sweep of 3,
    and three draws visit every particle. `sample4update_rescore` returns the
    same three rows on repeated calls and advances `sampled` and `updatecnt`
    once per call. It returns an empty set under `allow_empty` on a range
    without cohort rows. The next global draw, and the next unit draw, take no
    row of the previous cohort.
  - `simple_class_sample_io_tester`: `group` round-trips, including on an
    empty class.
  - `simple_accum_blend_tester`, `trailing_reconstruction_blend` sub-suite: a
    fixed cohort held for four iterations follows the closed form
    `1 − (1−u)^k` on the even and odd chains, with the density at full mass.
- Release build (`scratch/build_release`; logs `build_release_p1.log`,
  `build_release/make_p1.out`): the same five linker "requires executable
  stack" warnings as the base, no compiler warnings.
- Debug build by `compile_debug.sh` (`build_debug_p1.log`): six linker
  warnings, the same six as a Debug build of the base source
  (`scratch/build_debug_base`). The fast test gate passed 14 of 14.
- Unit areas of the touched code, all passing: `unit_core`, `unit_ori`,
  `unit_image`, `unit_project`, `unit_ui`, `unit_heterogeneity`,
  `unit_reconstruction` and `unit_parallel` (logs `p1_unit_*.log`, and
  `p1b_unit_*.log` after the review fix).
- `scripts/check_test_registry.py` passes.
- The searches the plan requires hold. No sampling path reads `partition` or
  `clust_crit`: the remaining reads are `selection`, `sample_classes`,
  `cluster_cavgs` and stack matching. `get_proj_sample_stats`,
  `make_projdir_class_samples` and `prepare_refine3D_states_class_sampling` no
  longer exist.
- An independent review of the diff found one defect, now fixed: the
  `balance` field was first declared with four characters, which would have
  cut `class` to `clas` and rejected it. The unit tests call the samplers
  directly and could not see it; the `balance=class` smoke run below covers it.
- Smoke runs on beta-galactosidase (not acceptance runs; Phase 2 does those):
  - `solve3D nsample=1000 nstages=1` with the default `cavg` (log
    `case_smoke_s1.log`): 50 class averages in 20 groups (1.7 s), the unit
    table with a 16-iteration sweep, 1,015 particles per iteration, normal
    stop in 1:46.
  - The Phase 0 C command line (log `case_smoke_C.log`, table
    `table_smoke_C.md`): sweep 9 blocks against 12 planned, `maxits` 36.
    Within every frequency block `% PARTICLES UPDATED SO FAR` stays constant
    after the block's first iteration and grows at the next block's first
    iteration. The per-state `updatecnt == 1` share is printed, every active
    particle was updated by the march (no missing-update pass), normal stop
    in 21:40.
  - `solve3D balance=class nstages=1` after the fix (log
    `case_smoke_class.log`): 50 units in 50 groups, an 11-iteration sweep,
    normal stop in 1:43.

Exit: the Release and Debug builds show no new warnings; the fast gate and the
unit suites of the touched areas pass; `check_test_registry.py` passes; the
required searches hold.

### Phase 2: comparison runs, 2026-10-04/05, done after the maintainer's ruling

Files changed: this plan (status line, the "Convergence" section, the Phase 0
C item, the Phase 2 acceptance, one file-table row and this entry), moved to
`doc/refactoring_notes/completed/`; `src/main/simple_convergence.f90` (the
per-state `updatecnt == 1` diagnostic removed again, so the file is back to its
base content); `doc/policies/heterogeneity/refine3D_states_policy.md` (the
diagnostic's sentence replaced); the report
`view_balanced_fractional_sampling_refactoring_report.md` beside the plan.

The phase first stopped on the A′ acceptance (first pass below). The
maintainer's ruling of 2026-10-05 then removed the diagnostic and settled the
equivalence of `balance=class` by a fixed-seed paired test (second pass,
further below).

#### First pass (stopped)

All runs used the Release build of the Phase 1 code
(`scratch/build_release`), `nparts=4 nthr=6`, one after another on an
otherwise idle machine. Paths are under
`/home/elmlundho/agent_runs/balance_sampling/scratch`; the per-run logs are
`logs/case_NAME.log` and the per-iteration tables `logs/table_NAME.md`. The
queue is `queue_phase2.sh` (log `logs/queue_phase2.log`). The command lines
are the Phase 0 ones, with `balance` added as the plan prescribes:

- A′: `simple_exec prg=solve3D projfile=bgal.simple mskdiam=180 pgrp=d2 nsample=1000 nparts=4 nthr=6 balance=class`
- B′: the same with `balance=cavg` (`nclust` at its default of 20)
- C′: `simple_exec prg=refine3D_states projfile=bgal.simple nstates=3 nsample=1500 lpstart=12 lpstop=6 pgrp=d2 mskdiam=180 nparts=4 nthr=6`
  (`refine3D_states` sets `balance=cavg` itself), on a copy of the kept A1
  final project prepared by `run_case.sh` as in Phase 0

Results. "Final FSC 0.143" is the resolution where the two half-map
reconstructions of the final all-particle reconstruction stop agreeing, as in
Phase 0.

| Run | Final FSC 0.5 / 0.143 (Å) | Wall time | Limit (Phase 0) | Verdict |
|---|---|---|---|---|
| A′1 | 4.35 / 3.89 | 27:56 | ≤ 4.00 | pass |
| A′2 | 6.40 / 4.29 | 28:14 | ≤ 4.00 | fail |
| B′1 | 4.41 / 3.93 | 28:48 | ≤ 4.05 | pass |
| B′2 | 4.41 / 3.89 | 28:27 | ≤ 4.05 | pass |
| C′ | state 1: 5.10 / 4.08; state 2: 14.19 / 12.09 | 21:50 | none | completes |

The A′ acceptance fails: A′2 is 0.29 Å above the limit, and the A′ mean
(4.09 Å) is too. B′ passes. C′ meets every item of its acceptance:

- It completed.
- Its log shows the coverage report: 20 groups of 50 classes; one sweep needs
  9 frequency blocks and 12 are planned; visits per particle range from 1.41 to
  12.0 against a target of 4.27, so no warning was raised.
- `% PARTICLES UPDATED SO FAR`, read from the refine3D_states log at the three
  iterations of each block (`logs/case_Cp.log`), is constant within every
  block and grows at the next block's first iteration: 27.374, 48.243,
  62.369, 74.240, 83.998, 90.646, 94.444, 97.270, 98.789, 99.359, 99.739,
  100.000. Every active particle (4,212) was updated by the march; no
  missing-update pass ran.
- The per-state share of particles with update count 1 was printed at every
  iteration, for example `>>> SHARE UPDATED ONCE (UPDATECNT=1), STATE  1:
  0.245` and `STATE  2: 0.328` at the first iteration. The rescoring advances
  the update count once per iteration, as this plan specifies, so the share
  read 0 after a block's first iteration. In the final project it is 0.000 in
  both states (Phase 0: 0.016 and 0.006). The ruling below removed this
  diagnostic and this acceptance item.

C′ ran 36 iterations (4 × a sweep of 9), against 15 in Phase 0's C, because
the automatic iteration cap now comes from the sweep. Flex initialization
again kept 2 of the 3 requested states.

Why A′2 failed. `balance=class` gives every class its own group, so it
allocates exactly as `balance=yes` did. Three checks show this:

- The code paths are the same (Phase 1 unit test: `class` gives every class
  the same count capped at its population, which is the old rule).
- The early `solve3D` stages draw greedily, the same top-scoring particles of
  every class each iteration. That sample was identical, with no differing
  particle, in base run A4 and in `class` run A′4 (`cmp/A4/sampled_rows.txt`
  and `cmp/Ap4/sampled_rows.txt`, 1,000 particles each).
- The child command lines of those two runs differ only by `balance=yes`
  against `balance=class` and the dropped `partition=no`
  (`cmp/*/child_cline.txt`).

The failed run lags from the ab initio stages on. Its FSC 0.143 at the end of
stage 4 is 8.8 Å, against 5.5–6.0 Å in the good runs, before the stages in
which the sampling is stochastic. To measure the spread, two more baseline A
runs on a Release build of the base source (`scratch/build_base`) and two more
A′ runs were made, alternating (`queue_phase2_spread.sh`). They are
supplementary evidence, not acceptance runs:

| Code | balance | Final FSC 0.143 (Å), all runs |
|---|---|---|
| base | `yes` (A) | 3.84, 3.93, 3.89, 3.98 |
| new | `class` (A′) | 3.89, 4.29, 5.10, 3.89 |
| base | `yes partition=yes` (B) | 3.98, 3.84 |
| new | `cavg` (B′) | 3.93, 3.89 |

A′3 (5.10 Å) also diverged in the ab initio stages (stage 4: 9.1 Å). Counting
both modes, the new code gave 2 poor de novo runs out of 6 and the base code 0
out of 6. With identical sampling in the stages where the runs diverge, that
difference is not explained by this change; six runs per side cannot tell
chance from a real effect. Phase 0's two-run spread (0.09 Å), and the 3 %
floor, do not cover the de novo failure rate that `solve3D` shows on this data
set.

What the maintainer needs to decide: whether the A′ acceptance stands as
written (then this phase fails, and the de novo divergence needs its own
investigation), or is restated, for example on the median of more runs, or by
counting only runs whose ab initio stages converged. After a decision, Phase 2
needs only the plan move and the report; the C′ and B′ evidence is complete.

#### Second pass, under the maintainer's ruling of 2026-10-05

The ruling (recorded in `spec/rulings.md` of the run) has two parts.

First, the per-state `updatecnt == 1` diagnostic is removed: under cohorts it
always reads 0, and the coverage report and `% PARTICLES UPDATED SO FAR`
already show whether every particle is reached. It is gone from
`simple_convergence.f90`, from this plan's "Convergence" section, the Phase 0
C item and the Phase 2 acceptance, and from the `refine3D_states` policy. Checks
after the change:

- The Release build was rebuilt.
- The fast test gate on the Release tree passed 14 of 14
  (`logs/p2_fast_gate.log`).
- The unit areas `unit_ori`, `unit_heterogeneity`, `unit_parallel`,
  `unit_reconstruction` and `unit_image` passed (`logs/p2_unit_*.log`;
  driver `p2_ruling1.sh`).

Second, the equivalence of `balance=class` with the old `balance=yes` is
settled by a fixed-seed paired test. Case A was run once on the base build
(`scratch/build_base`, `balance=yes`) and once on the Phase 1 build
(`balance=class`), in sequence, both with `SIMPLE_SEED=20261005` exported to the
`solve3D` process. The distributed workers inherit it: their process
environments carry the same value. Driver `seed_pair.sh`, logs
`logs/case_Aseed_base.log` and `logs/case_Aseed_class.log`.

The comparison runs through `snap_watch.sh`, `seed_compare.sh` and
`seed_rounds.py`, with results in `cmp/seed/`:

- `snap_watch.sh` copied the project after every iteration's assignment,
  when the log prints `% PARTICLES SAMPLED THIS ITERATION`.
- Every particle's round marker (`sampled`) and update count (`updatecnt`)
  were compared across the two runs.
- A few copies landed one round late, or caught a file being written. Each
  round's particle set was therefore taken either from a snapshot holding that
  round's marker, or from the snapshots of the two neighbouring rounds. The
  update-count rise over two rounds, minus membership in the later round,
  leaves exactly the middle round's set.

Results:

- **Rounds:** all 133 rounds of all eight stages are known in both runs, and
  every round's particle set is identical in the two runs.
- **Cross-check:** for every interval between consecutive snapshots of one
  run, each particle's update-count rise equals the number of the other run's
  rounds that contain it (116 and 117 intervals, no mismatch).
- **Final projects:** they agree on `sampled` and `updatecnt` for all 5,513
  particles.
- **Stages:** both runs entered the same stages at the same iterations (1, 21,
  41, 58, 73, 85, 97, 109).

So `balance=class` reproduces the old `balance=yes` sampling particle for
particle, and is accepted as equivalent.

The two seeded runs nevertheless ended at different resolutions: final FSC
0.5 / 0.143 of 4.35 / 3.89 Å (base, 28:16) and 7.10 / 4.41 Å (`class`, 28:20).
Their alignments differ from the first iteration on (mean orientation change
45.86° against 44.96°), with identical particles. A fixed seed does not make
the multithreaded alignment reproducible, and that is where the outcomes part.

Open observation for the maintainer, not settled by this run: unseeded de
novo divergence was seen in 2 of 4 `class` runs and 0 of 4 base runs here, and
in 0 of 6 identical `solve3D` runs in the earlier trailing_halfmap run. The
seeded pair shows identical sampling ending in different resolutions, so the
spread arises in alignment, outside sampling.

Exit:
- **A′:** accepted under the ruling, by identical sampling (unseeded results
  3.89 and 4.29 Å; seeded pair 3.89 and 4.41 Å).
- **B′:** passes, at 3.93 and 3.89 Å against a limit of 4.05 Å.
- **C′:** completes; its log shows the coverage report; `% PARTICLES UPDATED SO
  FAR` stays constant within every frequency block and grows at the next
  block's first iteration.
- **Plan and report:** the plan is moved to `completed/` and the report is
  written beside it.

### Phase 3: same-seed reproducibility control, 2026-10-05, done

Added by the maintainer after Phase 2. It is not among the plan's phases: the
run's phase list defines it, and the ruling of 2026-10-05 records it. No source
change. Files changed: this plan (status line and this entry) and the report.

Question. In Phase 2 the seeded pair (base build with `balance=yes`, new build
with `balance=class`, both with `SIMPLE_SEED=20261005`) sampled identical
particles in all 133 rounds, yet ended at 3.89 Å and 4.41 Å, and their
alignments parted at iteration 1. Two explanations were possible: `solve3D` is
not reproducible under a fixed seed even on one binary, or the new binary
changes the alignment numerics.

Runs. Case A was run twice more with `SIMPLE_SEED=20261005`, one after the
other, exactly as the Phase 2 pair: the same scripts (`seed_pair2.sh`, a copy of
`seed_pair.sh` with new run names, plus `run_case.sh` and `snap_watch.sh`), the
same execution settings (`nparts=4 nthr=6`), the same command line, and an
idle machine.

- `Aseed_base2`: base build (`scratch/build_base`), `balance=yes`. Final FSC
  0.5 / 0.143 of 6.66 / 4.47 Å, wall time 28:07.
- `Aseed_class2`: Phase 1 build (`scratch/build_release`, the same binary as
  `Aseed_class`), `balance=class`. Final FSC 0.5 / 0.143 of 4.41 / 3.93 Å,
  wall time 28:23.

Logs: `/home/elmlundho/agent_runs/balance_sampling/scratch/logs/case_Aseed_base2.log`
and `case_Aseed_class2.log`.

Comparison. `seed_logcmp.py` compares the quantities each iteration prints:
the stage, the conical FSC area ratio (cFAR, a map-quality figure), the
orientation overlap, the orientation change, the shift increment, the FSC
0.143 resolution and the score, each with its full
average/deviation/minimum/maximum line. Output:
`scratch/logs/p3_seed_logcmp.txt`. "Iteration 1" below means the first
iteration of stage 1.

| Pair | First difference | Final FSC 0.143 |
|---|---|---|
| `Aseed_base` / `Aseed_base2` (same binary) | iteration 1: cFAR 0.4410 / 0.4314, mean orientation change 45.864° / 44.518°, score maximum 0.425 / 0.418 | 3.89 / 4.47 Å |
| `Aseed_class` / `Aseed_class2` (same binary) | iteration 1: cFAR 0.2579 / 0.0046, mean orientation change 44.964° / 45.629°, score minimum 0.355 / 0.358 | 4.41 / 3.93 Å |
| `Aseed_base` / `Aseed_class` (Phase 2, two binaries) | iteration 1: cFAR 0.4410 / 0.2579, mean orientation change 45.864° / 44.964° | 3.89 / 4.41 Å |

All three pairs differ from iteration 1 in the same quantities, and in every
one of the 134 iterations. `solve3D` is therefore not reproducible under a
fixed `SIMPLE_SEED`, even on one binary: the seed fixes the random draws,
including the particle sampling, but not the multithreaded alignment and
reconstruction numerics. The base/`class` outcome difference in Phase 2 is
not evidence against the change.

All case-A runs of this refactoring, counting as poor a final FSC 0.143 above
the Phase 0 limit of 4.00 Å:

| Code | Runs | Final FSC 0.143 (Å) | Poor |
|---|---|---|---|
| base, `balance=yes` | Phase 0 A1 and A2, Phase 2 A3 and A4, `Aseed_base`, `Aseed_base2` | 3.84, 3.93, 3.89, 3.98, 3.89, 4.47 | 1 of 6 |
| new, `balance=class` | A′1 to A′4, `Aseed_class`, `Aseed_class2` | 3.89, 4.29, 5.10, 3.89, 4.41, 3.93 | 3 of 6 |

`balance=class` samples exactly as `balance=yes`, so both rows run the same
sampling. The base code also produces poor de novo runs (`Aseed_base2`). The
earlier trailing_halfmap run had 0 poor runs in 6 identical `solve3D` runs, and
the view-partition runs here (B and B′) had 0 of 4. With these numbers, the
difference between 1 of 6 and 3 of 6 does not separate from chance.

Bulky outputs were deleted; the logs, the scripts and the text comparisons
remain. No job is left running.

Exit:
- **Runs:** both seeded control runs completed.
- **Comparisons:** both same-binary pairs were compared from iteration 1, and
  the first difference is recorded above.
- **Outcome:** both same-binary pairs differ from iteration 1, so the
  prescribed outcome applies. `solve3D` is not reproducible under a fixed seed,
  and the report's open item on de novo divergence is updated.
