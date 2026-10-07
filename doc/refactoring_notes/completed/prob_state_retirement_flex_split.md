# Retire `prob_state`: flex PCA initializes every particle multi-state split

Date: 2026-10-06. Status: completed 2026-10-07 (Phases 0, 1 and 2 done; report beside this plan).

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

### Ruling of 2026-10-06 during the run: `solve3D` keeps only the independent multi-state mode

Taken after the Phase 0 baseline showed that docked particle `solve3D` never
runs `prob_state`: at its split stage it gives a cohort of particles random
balanced state labels, reconstructs the state volumes, and hands the rest of
the run to `refine3D_states` (see Progress, Phase 0). This ruling replaces
rulings 1, 4, 7 and 8 above and the sections "Docked `solve3D` on particles"
and "Docked `solve3D_cavgs`" below; where they differ, it wins.

1. `solve3D` and `solve3D_cavgs` lose the `multivol_mode` key. `nstates=1` is
   the single-state run; `nstates>1` is today's independent mode (every state
   refined from its own start), with its current defaults. Docked mode (one
   consensus model up to a split stage, then states) is removed with all of
   its code: the split checkpoint
   (`src/main/solve/simple_solve3D_split_checkpoint.f90`), the cohort pass,
   the handoff to `refine3D_states`, the docked stage rules in the solve3D
   controller (split stage, widened pre-split sample,
   `PROB_NEIGH_MODE_DOCKED`), the `split_stage` parameter, the docked fields
   of the solve3D manifest, and the docked-checkpoint input of
   `refine3D_states`. `sticky_class_sampling` and the `sampled_only` path of
   `sample4update_class` lose their only producer and go too. Only docked
   multi-state code goes: "dock" and "docked" elsewhere (volume docking,
   `solve3D_addon` docking, `dock_vols`) are unrelated and stay.
2. Flex PCA is therefore not used inside `solve3D`. It stays the initializer
   of `refine3D_states` for input without state labels (all particles in state
   0 or 1). `prob_state`, `pose_policy=fixed` and `flex=no` are still retired
   as planned; `nstates>=3` applies to flex-initialized `refine3D_states`. No
   balanced class-average partition is needed.
3. `multivol_mode` stays as an internal parameter only where
   `refine3D_states` and `classify3D_refs` set it. Once `input_oris_fixed` is
   gone, Progress records whether the parameter has become removable; it is
   not removed in this run.
4. Validation is reduced to what the change can break; no new baselines are
   run. The Phase 0 docked cases (A1, A2, C) stay in Progress as context only.
   Phase 2 runs, on the new build with the Phase 0 execution settings and
   inputs: the single-state `solve3D nsample=2500 cavg_ini_ext=yes` from S
   once (it exercises the `prob`, `prob_neigh` and `balance=cavg` paths the
   removals touch; final FSC 0.143 within 3 % of the Phase 0 value, 3.84 Å);
   B (`refine3D_states`) as in Phase 0 once (best-state FSC 0.143 within 3 %
   of its Phase 0 value); `solve3D nsample=2500 cavg_ini_ext=yes nstates=3`
   from S and `solve3D_cavgs nstates=3`, once each, which must complete with
   every state populated (smoke tests of the now implicit independent mode,
   no resolution criterion). Plus the fast test gate, the unit suites of the
   touched areas and `scripts/check_test_registry.py`, as in Phase 1's exit.
5. File-table rows are added before editing files the removal touches beyond
   the table (for example `src/main/ui/simple/simple_ui_solve3D.f90`,
   `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`,
   `src/main/ori/simple_oris_sampling.f90`,
   `src/main/sigma2/simple_sigma2_bootstrap.f90`,
   `doc/policies/3D/solve3D_addon_policy.md`). Documents: the "Multi-state
   modes" section of `doc/algorithms/solve3d.md` keeps single and independent
   only; the docked sections of `doc/policies/3D/solve3D_policy.md`,
   `doc/policies/3D/solve3D_cavgs_policy.md` and
   `doc/policies/importance_sampling_fractional_update_policy.md` go.

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

Replaced by the ruling of 2026-10-06 during the run: docked mode is removed,
so this section no longer applies. Its description of today's code is also
incomplete: after the steps below, the split stage does not run `prob_state`
but hands off to `refine3D_states` (`flex=no`, `pose_policy=local`), which
skips its `prob_state` initialization phase because the particles already
carry state labels.

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

Replaced by the ruling of 2026-10-06 during the run: docked mode is removed
and `nstates>1` runs the independent mode, so no balanced partition is needed.

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

Replaced by item 4 of the ruling of 2026-10-06 during the run. Phase 2 runs,
with the command lines recorded in Progress (Phase 0) and the same
`scratch/keep/` inputs: the single-state `solve3D` from S (final FSC 0.143 at
most 3.96 Å, i.e. within 3 % of 3.84 Å); B (best-state FSC 0.143 at most
4.15 Å, within 3 % of 4.03 Å); `solve3D_cavgs nstates=3` and
`solve3D nstates=3` (started from the `solve3D_cavgs nstates=3` result, see
the open point in Progress, Phase 0) as smoke tests (complete, every state populated);
the fast test gate, the unit suites of the touched areas and
`scripts/check_test_registry.py`. Then the plan moves to
`doc/refactoring_notes/completed/` with a self-contained report beside it.
The original text follows for reference.

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
| `src/main/ui/simple/simple_ui_solve3D.f90` | `solve3D` and `solve3D_cavgs` lose the `multivol_mode` and `split_stage` inputs; help text says "nstates>1" instead of "independent" (ruling of 2026-10-06) | 1 |
| `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` | the matcher stops passing the docked sticky-cohort flag to the sampler | 1 |
| `src/main/ori/simple_oris_sampling.f90`, `src/main/ori/simple_oris.f90`, `src/main/ori/simple_oris_tester.f90` | `sample4update_class` loses its `sampled_only` (sticky cohort) argument, its only producer being docked mode; its unit test goes with it | 1 |
| `src/main/sigma2/simple_sigma2_bootstrap.f90` | stops deleting the removed `sticky_class_sampling` key from its child command line | 1 |
| `src/main/solve/simple_solve3D_manifest_tester.f90` | manifest unit test follows the removed docked fields | 1 |
| `src/main/flex/run/simple_flex_pca_state_service.f90`, `src/main/flex/run/simple_flex_pca_application.f90` | comments that name the retired `flex=yes` key | 1 |
| `scripts/memory/benchmark_solve3d.py` | the solve3D memory benchmark stops passing `multivol_mode` on the `solve3D` command line | 1 |
| `doc/policies/3D/solve3D_addon_policy.md`, `doc/how2s/how2_process_heterogeneous_datasets.md`, `doc/policies/NU/nonuniform_filtering_policy.md` | living documents that mention docked mode, `multivol_mode` on `solve3D`, `prob_state`, `pose_policy=fixed` or `flex=no` | 1 |

## Progress

### Phase 0: baseline (2026-10-06, done)

No source file changed; this plan is the only file edited (the ruling of
2026-10-06 recorded above, the replaced sections marked, this entry). All
numbers come from an out-of-tree Release build with `BUILD_TESTS=ON` of the
unchanged base (commit `6914834d3`, clean tree), built with GCC 15.2
(gcc-toolset-15), CMake 3.26.5 and `USE_ARCHOPT=ON`, the same configuration
as the maintainer's `~/src/SIMPLE/build`. The build has no compiler
warnings and five linker warnings "requires executable stack"
(`json_value_module.f90.o` three times, `simple_commanders_cavgs.f90.o`,
`simple_stream_stage_pool2D.f90.o`); that is the warning baseline for
Phase 1.

Evidence is under `/home/elmlundho/agent_runs/prob_state_flex/scratch`:
`build_release.sh`, `run_case.sh` (the exact command lines), `queue_phase0.sh`
(the order), and in `logs/` per run `case_<name>.log` (full output, each line
prefixed by its time in seconds, `/usr/bin/time -v` at the end),
`summary_<name>.txt` (populations, final resolutions, wall time, memory) and
`table_<name>.md` (per-iteration table). The runs ran one at a time on the
otherwise idle 24-core machine.

**Inputs.** Each run starts from a fresh copy of
`/data2/bgal_addon/control/bgal_addon.simple`: beta-galactosidase after 2D
classification, 5,513 particles, box 256 at 1.275 Å per pixel, 50 selected
class averages. Run S (`solve3D_cavgs`, single state) is the common start of
every `solve3D` run; its run directory is kept as `scratch/keep/S/1_solve3D_cavgs`.
The single-state `solve3D` from S is kept as `scratch/keep/single/1_solve3D`;
it is the input of case B. A project's output segment refers to its maps by
absolute path inside the directory that wrote them, so `run_case.sh` recreates
`scratch/runs/S/1_solve3D_cavgs` (and, for B, `scratch/runs/single/1_solve3D`)
from `scratch/keep/` before each run that reads a kept project; a run depends
only on `scratch/keep/`.

**Command lines.** Every run sets `SIMPLE_PATH` to the Release build and puts
its `bin` first on `PATH`; the project file is the copy named `bgal.simple`
in the run directory `scratch/runs/<name>`. Phase 2 repeats these exactly.

| Run | Starts from | Command line |
|---|---|---|
| S | fresh copy of the input project | `simple_exec prg=solve3D_cavgs projfile=bgal.simple mskdiam=180 pgrp=d2 nparts=4 nthr=6` |
| single | `keep/S/1_solve3D_cavgs/bgal.simple` | `simple_exec prg=solve3D projfile=bgal.simple mskdiam=180 pgrp=d2 nsample=2500 cavg_ini_ext=yes nparts=4 nthr=6` |
| A1, A2 | `keep/S/1_solve3D_cavgs/bgal.simple` | `simple_exec prg=solve3D projfile=bgal.simple mskdiam=180 pgrp=d2 nsample=2500 cavg_ini_ext=yes nparts=4 nthr=6 multivol_mode=docked nstates=3` |
| B | `keep/single/1_solve3D/bgal.simple` | `simple_exec prg=refine3D_states projfile=bgal.simple nstates=3 nsample=1500 lpstart=12 lpstop=6 pgrp=d2 mskdiam=180 nparts=4 nthr=6` |
| C | fresh copy of the input project | `simple_exec prg=solve3D_cavgs projfile=bgal.simple mskdiam=180 pgrp=d2 nparts=4 nthr=6 multivol_mode=docked nstates=3` |

**Results.** Population is the number of particles with that state label in
the final project (state 0 means excluded). Resolutions are the Fourier shell
correlation (FSC) thresholds 0.5 and 0.143 of the final reconstruction of
each state, from its `RESOLUTION_FINAL_STATE<nn>` file; for S and C the half
sets are the even and odd class averages, not particles. Every run exited
with status 0.

| Run | Wall time | Peak memory | State: population, FSC 0.5 / 0.143 (Å) |
|---|---|---|---|
| S | 10:56 | 3.0 GB | 1: 5,513 particles (50 class averages), 7.77 / 6.80 |
| single | 14:38 | 3.2 GB | 1: 5,513, 4.35 / 3.84 |
| A1 | 55:56 | 5.4 GB | 1: 1,330, 7.25 / 4.35; 2: 1,146, 10.20 / 7.77; 3: 3,037, 4.47 / 3.80 |
| A2 | 56:01 | 5.4 GB | 1: 868, 9.07 / 6.40; 2: 3,777, 4.41 / 3.80; 3: 868, 10.20 / 7.59 |
| B | 18:55 | 4.8 GB | 1: 3,460, 5.02 / 4.03; 2: 726, 14.19 / 12.09; 0: 1,327 |
| C | 11:18 | 5.3 GB | 1: 1,316 (9 class averages), 8.82 / 5.26; 2: 2,688 (19), 8.82 / 6.94; 3: 1,509 (22), 20.40 / 9.07 |

What the runs did, read from their logs:

- single: entered at stage 4 (the stage after the symmetry search, because
  `cavg_ini_ext=yes` takes the orientations of S), stages 4 to 8 at low-pass
  limits 9.1, 7.4, 7.3, 5.5 and 4.5 Å.
- A1 and A2 (docked `solve3D`, context only under the ruling): stage 5, the
  pre-split stage, ran `refine=prob` at a cropped box of 130 pixels and
  7.6 / 7.4 Å. At stage 6 the split checkpoint formed a cohort with one
  `refine=prob` pass (4,962 of 5,513 particles), gave it random balanced
  labels, reconstructed three state volumes from the cohort, and handed the
  remaining 49 iterations to `refine3D_states` (`pose_policy=local`,
  `flex=no`), which skipped its `prob_state` initialization phase because the
  particles already had labels and ran 17 frequency blocks of
  `refine=prob_neigh` from 8.4 to 4.5 Å, then the final missing-update pass
  on the 549 particles never sampled. `prob_state` never ran. This is the
  observation that led to the ruling of 2026-10-06 during the run.
- B (`refine3D_states` on the single-state project): flex PCA initialized the
  states; it delivered 3, then its occupancy floor (at least 2,000 effective
  particles per state) dropped one state of 1,327 particles to state 0, so
  the run continued with 2 states (`REFINE3D_STATES FLEX_PCA INITIALIZED
  NSTATES: 2`). No initialization phase; 11 frequency blocks from 12 to 6 Å,
  32 iterations at an update fraction of 0.27 (1,500 of 5,513 particles per
  iteration).
- C (docked `solve3D_cavgs`, context only): stages 1 to 5 single-state, then at
  stage 6 random labels and three iterations of `refine=prob_state`, then 12
  iterations of `refine=prob_neigh`. It is the only Phase 0 run that ran
  `prob_state`.

**Phase 2 limits** (item 4 of the ruling): the single-state `solve3D` final
FSC 0.143 at most 3.96 Å (3.84 Å plus 3 %); B's best-state FSC 0.143 at most
4.15 Å (4.03 Å plus 3 %). For the record, the A pair's best states both
reached 3.80 Å.

**Open point for Phase 2.** The ruling's smoke test `solve3D nsample=2500
cavg_ini_ext=yes nstates=3` from S cannot run as written: S carries one
state, and `solve3D` with `cavg_ini_ext=yes` and `nstates>1` requires the
project's particles to carry exactly `nstates` populated states. On the base
build it stops at once with "cavg_ini_ext=yes with nstates>1 requires
matching existing ptcl3D state assignments" (`logs/case_probe_indep.log`).
Phase 1 does not change that check (it is not docked code), so Phase 2 starts
this smoke test from the final project of the `solve3D_cavgs nstates=3` smoke
test, whose particles carry three populated states, with the rest of the
command line unchanged. That keeps the test on the independent mode with
class-average initialization, as the ruling intends.

### Phase 1: implementation (2026-10-06)

#### Inventory of the base code, taken before editing

The base code was read in full for every name the run retires. The
design that applies is the ruling of 2026-10-06 during the run (docked mode
removed), so flex PCA inside `solve3D` and the balanced class-average
partition are not built; the questions the plan asked about them (can flex
run on a cohort at a crop, does `gen_labelling(..., 'uniform')` balance) no
longer arise.

- `refine=prob_state` (keep the projection direction, score the pose against
  every state, draw the state). Emitted by three producers: the docked split
  stage of the solve3D stage controller (`set_refine3D_mode_policy` in
  `simple_solve3D_controller.f90`, through `docked_split_stage`), the
  initialization stage of `refine3D_states` for input without state labels
  when `flex=no`, and every stage plus the final missing-update pass of
  `refine3D_states pose_policy=fixed`. Consumed by: the parallel refine3D
  strategy (`l_prob_state_mode`, which only ever guards the `prob_neigh`
  branches, so removing it changes no other mode), the matcher and strategy
  allocator (`prob_state` cases), the probabilistic search
  (`l_fixed_projection` in `simple_strategy3D_prob.f90`, restoring the
  projection direction), the probability-table commanders
  (`l_state_only` in `simple_commanders_prob.f90`), and the state-only table
  in `simple_eul_prob_tab.f90` (`state_tab` storage, `new_state`,
  `fill_tab_state_only`, `fill_tab_state_only_range`, `read_state_tab`,
  `write_state_tab`, already without callers, and `state_assign`); nothing
  else uses that storage. `simple_eul_prob_tab_neigh.f90` mentions
  `prob_state` in a comment only. The parameter phases accept
  `refine=prob_state` as a probabilistic-alignment mode; the UI of
  `refine3D` and of SINGLE's `nano3D` list it as a choice.
- `input_oris_fixed` (the internal `multivol_mode` of `pose_policy=fixed`):
  set and checked only in `refine3D_states`; read in the shared
  `prep4prob` of `simple_strategy3D_srch.f90` (which skips recording the
  projection index) and in `simple_strategy3D_prob.f90`. Retiring it with
  `prob_state` is safe because it only ever ran with `prob_state`.
- `pose_policy`: parameter, parser, validation (`fixed|local|global`), the
  `refine3D_states` UI, and the `refine3D_states` commander, which maps it to
  the internal `multivol_mode` and `prob_neigh_mode`. Also set to `local` by
  the docked handoff in `simple_commanders_solve3D.f90`.
- The `flex` key: parameter, parser, yes/no validation, the
  `refine3D_states` UI, the `refine3D_states` commander (explicit
  `flex=yes|no`, otherwise derived from whether the project carries state
  labels and written back to the command line), the docked handoff
  (`flex=no`), and comments in two flex PCA files. The rule "flex
  initialization needs `nstates>=3`" already exists in the commander and in
  `run_flex_pca`; only its messages, which suggest `flex=no`, change.
- `l_init_state_assignment` (the `prob_state` initialization phase): only in
  `refine3D_states`, with `init_stage_cap`, `init_stage_minits`, the init
  stage and its frequency march, the stochastic startup branch of
  `initialize_state_volumes`, and the state-0/1 branch of
  `set_refine3D_states_nstates`. The variable `init_sweep_iters` also sets
  the automatic `maxits` (four updates per particle times one sweep), so it
  stays under a neutral name.
- `fill_tab_state_only`, `state_assign`, `l_state_only`,
  `l_prob_state_mode`: as listed under `prob_state`; no other callers.
- `docked_split_stage`, `randomize_states`,
  `build_solve3D_split_checkpoint` and the rest of docked mode: the docked
  branches of `exec_solve3D_cavgs` and `exec_solve3D`, the handoff to
  `refine3D_states`, the split checkpoint module, the controller's docked
  rules (`HET_DOCKED_STAGE`, `NSAMPLE_HET_SPLIT_CAP`,
  `PROB_NEIGH_MODE_DOCKED`, `solve3D_het_docked_stage`,
  `solve3D_docked_cohort_active`, `calc_docked_multistate_max_sampling`, the
  docked case of the mode, trailing and stage-control policies, the 0.98
  `frac_best` of multi-state stages 7 and 8), the `current_sample_only`
  option of `calc_rec` (used only by the checkpoint), the `split_stage`
  parameter and UI entries, and the manifest fields `base_multivol_mode` and
  `split_stage` (plus the input keys `multivol_mode` and `split_stage`).
  `randomize_states` serves only docked mode; `gen_labelling(...,
  'uniform')` stays, it labels independent multi-state starts.
- `sticky_class_sampling`: produced only by docked mode (the controller and
  the handoff); consumed by the matcher's sampling call
  (`sampled_only` of `sample4update_class`) and deleted from child command
  lines in `calc_rec`, `calc_frozen_rec` and the sigma bootstrap.
- `multivol_mode` outside `solve3D`: `refine3D_states` sets
  `input_oris_refine` (after this phase always), `classify3D_refs` sets
  `independent`, and `simple_convergence.f90` tests the `input_oris` prefix
  to pick the state-overlap convergence rule. Nothing else reads its value.
  The solve3D memory benchmark script passes `multivol_mode` on the
  `solve3D` command line.

Decisions taken from the inventory:

- With docked mode gone, `solve3D` and `solve3D_cavgs` decide by `nstates`
  alone: one state is the single-state run, more is today's independent mode
  with its defaults (five stages and a 6 Å final low-pass for particles
  unless given, `refine=shc` in stages 1 and 2 for particles, stochastic
  sampling from stage 4, trailing reconstruction at the last stage only, the
  greedy missing-update pass and the final reconstruction always). Both
  commanders reject `multivol_mode` on their command line with a message
  that says so, because the parser accepts every known key for every program
  and a silently ignored mode would mislead.
- The manifest reader keeps rejecting fields it does not know, as its unit
  test requires for removed features under release 4; a manifest written by
  an earlier docked or independent run therefore no longer reads, which
  affects only `solve3D_addon` replays of runs made before this change.
  Project files (`.simple`) are untouched and read as before.
- The generated Fortran indexes under `doc/code_overview/fortran-indexes/`
  are not regenerated in this run; they are produced by
  `scripts/gen_fortran_indexes.pl` and carry no hand-written content.

#### What changed

- `refine=prob_state` is gone from every layer: the solve3D stage controller
  no longer emits it, the parallel refine3D strategy, matcher and strategy
  allocator have no case for it, the probabilistic search no longer restores
  a fixed projection direction, the probability-table commanders always fill
  and read the full table, and the state-only table of
  `simple_eul_prob_tab.f90` (storage, constructor, fill, read, write and the
  balanced state draw) is deleted. The parameter phases no longer count it
  as a probabilistic-alignment mode; the `refine` choices of `refine3D` and
  SINGLE's `nano3D` no longer list it. The multi-state minimum-population
  constant of the strategy, which applies to every probabilistic multi-state
  run, is renamed from `PROB_STATE_MIN_POP` to `PROB_MULTISTATE_MIN_POP`.
- `refine3D_states`: input without state labels (state 0 or 1) is always
  initialized by flex PCA and requires `nstates >= 3` ("state=0/1 input
  requires nstates >= 3 for flex_pca initialization"); labelled input refines
  its states and uses or reconstructs its state maps. The `flex` key, the
  `prob_state` initialization phase (its stage cap, minimum iterations, init
  frequency march and stochastic startup volumes) and `pose_policy=fixed`
  with `multivol_mode=input_oris_fixed` are removed; `pose_policy` is
  `local|global` in the parameter validation, the commander and the UI. Every
  run starts at the `prob_neigh` frequency march, trailing reconstruction is
  always requested, and the final missing-update pass is `refine=greedy`.
  The `startit` offset that only the docked handoff used is dropped. The
  variable that counts the draws of one sweep keeps setting the automatic
  `maxits` under the name `sweep_iters`.
- Docked mode is removed from `solve3D` and `solve3D_cavgs` (ruling of
  2026-10-06): the split checkpoint module is deleted, together with the
  cohort pass, the handoff to `refine3D_states`, the controller's docked stage
  rules and constants, `randomize_states`, the current-sample option of
  `calc_rec`, the `split_stage` parameter and UI entries, the manifest fields
  `base_multivol_mode` and `split_stage` and the manifest input keys
  `multivol_mode` and `split_stage`, `sticky_class_sampling` (parameter,
  derived flag, validation, controller emission and the deletes in
  `calc_rec`, `calc_frozen_rec` and the sigma bootstrap) and the
  `sampled_only` path of `sample4update_class`. The controller's mode,
  trailing and stage-control policies and the commanders now decide by
  `params%nstates` (one state: single-state rules; more: the former
  independent rules, unchanged); both commanders refuse `multivol_mode` on
  their command line ("solve3D takes no multivol_mode: nstates>1 refines
  independent states"); `solve3D` state continuation refuses `nstates>1`, as
  it did through the forced single mode before. The manifest replay of
  `solve3D_addon` no longer sets a mode. The UI help of `nstages` and
  `lpstop` says "nstates>1" instead of "multivol_mode=independent". The
  solve3D memory benchmark script no longer passes `multivol_mode`.
- Unit tests: the manifest test now checks that a replay sets no
  `multivol_mode` and that a manifest carrying the removed `split_stage` field
  is refused; the `sampled_only` case of the sampler test is removed with the
  argument. The plan's two planned tests no longer apply under the ruling
  (there is no balanced class-average partition and no docked particle run);
  their place is taken by the command-line checks below.
- `multivol_mode` after the change: `refine3D_states` always sets
  `input_oris_refine` and refuses a user value, `classify3D_refs` sets
  `independent`, and only `simple_convergence.f90` reads the value (the
  `input_oris` prefix selects the state-overlap convergence rule of
  `refine3D_states`). The parameter has therefore become removable, but only
  together with a replacement for that convergence discriminator; it is kept
  in this run, as the ruling asks.
- Documents: the `refine3D_states`, `solve3D`, `solve3D_cavgs`,
  importance-sampling and fractional-update, NU-filtering and `solve3D_addon`
  policies, the algorithm chapters `refine3d.md`, `solve3d.md`,
  `sampling_and_fractional_updates.md`, the heterogeneity-analysis overview
  and `refine3d_states.md`, the heterogeneous-datasets how-to, and the skills
  `simple-refine3d` and `simple-solve3d-importance-sampling` with their
  references describe the retired names as gone, or no longer mention them.
  The generated Fortran indexes are not regenerated.

#### Evidence

All paths are under `/home/elmlundho/agent_runs/prob_state_flex/scratch`.

- Builds. Base Debug build and fast gate (`logs/debug_base.log`): 14 of 14
  fast suites pass; no compiler warnings; six linker warnings "requires
  executable stack" (`json_value_module.f90.o` three times,
  `simple_stream_stage_pool2D.f90.o`, `simple_commanders_cavgs.f90.o`,
  `simple_persistent_worker.f90.o`). Phase 1 Debug build with the repository's
  `./compile_debug.sh` (`logs/debug_p1_a.log`, and the final tree in
  `logs/debug_p1_final.log`): compiles, no compiler warnings with `-Wunused`
  on, the same six linker warnings, fast gate 14 of 14. Phase 1 Release build
  (`build_release.sh`, `logs/build_release_p1.log`): the same five linker
  warnings as the Phase 0 base build, no compiler warnings.
- Tests (`logs/tests_p1.log`, Release tree): fast gate 14 of 14, including
  the ori suite (sampler) and the project suite (solve3D manifest, 68 checks);
  library test `lib_heterogeneity` passes. The high-level test
  `solve3D_addon` fails at its first base `solve3D` with "balance=cavg needs
  selected class averages": its simulated particles carry no class averages.
  The same failure occurs on a Release build of the unchanged base
  (`logs/build_base_addon.log`), so it predates this run and is not caused
  by it; it is left for the maintainer.
- `scripts/check_test_registry.py` passes (test registrations unchanged).
- Text checks: no `prob_state`, `input_oris_fixed`, `fill_tab_state_only` or
  `state_assign` in `src/` (an unrelated helper whose name contained
  `state_assign`, `read_multistate_assignment_coverage`, is renamed
  `read_multistate_update_coverage`); no `flex` key or `flex=no` path; the
  `pose_policy` choices are `local|global`; no `THROW_HARD` continued with
  `//&` among the changed lines; no file mode changes.
- Command-line checks on the Release build with the kept single-state project
  (`logs/checks_p1/`): `refine3D_states nstates=2` stops with "state=0/1
  input requires nstates >= 3 for flex_pca initialization";
  `pose_policy=fixed` stops with "supports pose_policy=local|global";
  `flex=no` and `split_stage=6` are unknown arguments; `solve3D
  multivol_mode=docked` stops with "solve3D takes no multivol_mode".
  `refine3D refine=prob_state` is refused by the workers ("refinement mode:
  prob_state unsupported") but the distributed master keeps waiting; any
  unknown `refine` value behaves the same way in the base code, and refusing
  it in the master would need a list of refine modes kept beside the matcher,
  a new contract this run does not add. Recorded for the maintainer.

- Independent review of the diff (a separate reviewer reading the code
  against this plan): no regression found. Replacing the mode tests by
  `nstates` tests is exact, because `nstates` is set only from the command
  line and the old code refused every mode that disagreed with it. One
  behaviour changes for the better: `solve3D cavg_ini=yes nstates>1` used to
  fail, because the nested `solve3D_cavgs` inherited
  `multivol_mode=independent` without `nstates`; it now runs a single-state
  class-average initialization and then labels the states, as intended.
- Final trees: Debug (`logs/debug_p1_final.log`) and Release
  (`logs/build_release_p1_final.log`, fast gate in
  `logs/fast_gate_release_final.log`) built from the final working tree; the
  Release tree is the one Phase 2 uses.

#### Exit items

- Release and Debug builds without new warnings: met (see Builds).
- Fast test gate and the unit suites of the touched areas pass: met (Debug
  and Release, 14 of 14; the touched areas are covered by the ori, project,
  heterogeneity, parallel and alignment suites).
- `scripts/check_test_registry.py` passes: met.
- Text checks hold: met.
- Open for Phase 2 (from Phase 0): the smoke test `solve3D nstates=3
  cavg_ini_ext=yes` needs a start whose particles already carry three
  populated states; it starts from the `solve3D_cavgs nstates=3` result.

#### Files changed in Phase 1

- `.github/skills/simple-refine3d/SKILL.md`
- `.github/skills/simple-refine3d/references/bayesian-model.md`
- `.github/skills/simple-solve3d-importance-sampling/SKILL.md`
- `.github/skills/simple-solve3d-importance-sampling/references/coupling-map.md`
- `doc/algorithms/heterogeneity_analysis/README.md`
- `doc/algorithms/heterogeneity_analysis/refine3d_states.md`
- `doc/algorithms/refine3d.md`
- `doc/algorithms/sampling_and_fractional_updates.md`
- `doc/algorithms/solve3d.md`
- `doc/how2s/how2_process_heterogeneous_datasets.md`
- `doc/policies/3D/solve3D_addon_policy.md`
- `doc/policies/3D/solve3D_cavgs_policy.md`
- `doc/policies/3D/solve3D_policy.md`
- `doc/policies/NU/nonuniform_filtering_policy.md`
- `doc/policies/heterogeneity/refine3D_states_policy.md`
- `doc/policies/importance_sampling_fractional_update_policy.md`
- `doc/refactoring_notes/planned/prob_state_retirement_flex_split.md`
- `scripts/memory/benchmark_solve3d.py`
- `src/main/commanders/simple/simple_commanders_prob.f90`
- `src/main/commanders/simple/simple_commanders_refine3D.f90`
- `src/main/commanders/simple/simple_commanders_solve3D.f90`
- `src/main/flex/run/simple_flex_pca_application.f90`
- `src/main/flex/run/simple_flex_pca_state_service.f90`
- `src/main/ori/simple_oris.f90`
- `src/main/ori/simple_oris_sampling.f90`
- `src/main/ori/simple_oris_tester.f90`
- `src/main/params/simple_parameters.f90`
- `src/main/params/simple_parameters_parse.f90`
- `src/main/params/simple_parameters_phases.f90`
- `src/main/sigma2/simple_sigma2_bootstrap.f90`
- `src/main/solve/simple_solve3D_controller.f90`
- `src/main/solve/simple_solve3D_manifest.f90`
- `src/main/solve/simple_solve3D_manifest_tester.f90`
- `src/main/solve/simple_solve3D_split_checkpoint.f90` (deleted)
- `src/main/solve/simple_solve3D_utils.f90`
- `src/main/strategies/parallelization/simple_refine3D_strategy.f90`
- `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`
- `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90`
- `src/main/strategies/search/probabilistic/simple_strategy3D_prob.f90`
- `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`
- `src/main/strategies/search/simple_strategy3D_alloc.f90`
- `src/main/strategies/search/simple_strategy3D_matcher.f90`
- `src/main/strategies/search/simple_strategy3D_srch.f90`
- `src/main/ui/simple/simple_ui_heterogeneity.f90`
- `src/main/ui/simple/simple_ui_refine3D.f90`
- `src/main/ui/simple/simple_ui_solve3D.f90`
- `src/main/ui/single/single_ui_nano3D.f90`

### Phase 2: comparison runs (2026-10-07, done)

The runs follow item 4 of the ruling of 2026-10-06 during the run: no new
baselines, the single-state `solve3D` and case B compared with Phase 0, and
two smoke tests of the now implicit independent multi-state mode. They ran one
at a time on the otherwise idle machine (load 0.05 at the start), on the
Release build of the final Phase 1 tree, with the Phase 0 command lines and
the `scratch/keep/` inputs unchanged. Evidence is under
`/home/elmlundho/agent_runs/prob_state_flex/scratch`: `queue_phase2.sh` (the
order), `run_case.sh` (the command lines; it gained the kind `solve3D_C3` for
the second smoke test), and in `logs/` the files `case_<name>.log`,
`summary_<name>.txt` and `table_<name>.md` for `single_p2`, `B_p2`, `C3` and
`I3`, plus `queue_phase2.log`. Every run exited with status 0.

| Run | Command line (Phase 0 table, or as noted) | Wall time | State: population, FSC 0.5 / 0.143 (Å) | Acceptance |
|---|---|---|---|---|
| single (single-state `solve3D` from S) | as in Phase 0 | 14:31 | 1: 5,513, 4.35 / 3.80 | at most 3.96 Å: met (Phase 0: 3.84 Å) |
| B (`refine3D_states` on the kept single-state project) | as in Phase 0 | 18:58 | 1: 3,490, 5.02 / 4.03; 2: 723, 15.54 / 12.55; 0: 1,300 | best state at most 4.15 Å: met (Phase 0: 4.03 Å) |
| C3 (`solve3D_cavgs nstates=3`) | the S line plus `nstates=3` | 14:49 | 1: 4,191 (30 class averages), 7.77 / 6.80; 2: 272 (5), 36.27 / 23.31; 3: 1,050 (15), 25.11 / 19.20 | completes with every state populated: met |
| I3 (`solve3D nstates=3` from the C3 result) | the single line plus `nstates=3`, project from C3 | 37:19 | 1: 4,189, 4.66 / 3.93; 2: 471, 21.76 / 13.60; 3: 853, 14.84 / 10.53 | completes with every state populated: met |

What the runs did, read from their logs:

- single ran the same stages and low-pass limits as in Phase 0 (stages 4 to
  8 at 9.1, 7.4, 7.3, 5.5 and 4.5 Å).
- B: "REFINE3D_STATES STATE INITIALIZATION BY FLEX PCA"; flex delivered three
  states, its occupancy floor (at least 2,000 effective particles per state)
  dropped one of 1,300 particles to state 0, and the run continued with
  "REFINE3D_STATES FLEX_PCA INITIALIZED NSTATES: 2", then 32 `prob_neigh`
  iterations in 11 frequency blocks from 12 to 6 Å ("STAGE ITERATIONS
  PROB_NEIGH/TOTAL: 32/32"), as in Phase 0 but without the empty
  initialization stage.
- C3 ran stages 1 to 7 with `prob_neigh`, `prob` and `prob_neigh` and no
  `prob_state` stage (Phase 0's docked C had three `prob_state` iterations at
  its split).
- I3 was started from the C3 result rather than from S, as decided in Phase 0,
  because `cavg_ini_ext=yes` with three states needs particles that already
  carry three populated states. It reported "SOLVE3D INDEPENDENT MULTI-STATE
  DEFAULT NSTAGES: 5" and "DEFAULT LPSTOP: 6.0 A", ran stages 4 and 5 with
  `refine=prob`, and its coverage check found every active particle updated
  ("MULTISTATE ASSIGNMENT COVERAGE NSTATES=3 UPDATED/ACTIVE/MISSING:
  5513/5513/0"), so no missing-update pass was needed.
- Tests on the same tree: the fast gate passes 14 of 14
  (`logs/fast_gate_p2.log`) and `scripts/check_test_registry.py` passes; no
  file mode changed.
- Observation, not caused by this run: in B the sampling line "ACTIVE
  PARTICLES/SAMPLE TARGET/UPDATE_FRAC: 5513/1500/0.2721" counts the 1,300
  particles flex dropped to state 0, so the realized update fraction over the
  4,190 remaining particles is about 0.36, not 0.27; Phase 0 shows the same.

Bulky outputs were deleted after the numbers were recorded (`runs/`, and the
C3 result kept only for I3); `scratch/keep/S` and `scratch/keep/single` stay as
the inputs that reproduce these runs.

Files changed in Phase 2: this plan (status line and this entry), moved from
`doc/refactoring_notes/planned/` to `doc/refactoring_notes/completed/`, and
the new report `doc/refactoring_notes/completed/prob_state_retirement_flex_split_report.md`.
