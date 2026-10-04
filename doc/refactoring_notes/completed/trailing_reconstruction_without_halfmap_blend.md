# Trailing reconstruction without the finished-halfmap blend

Date: 2026-10-03 (revised the same day).

Status: completed 2026-10-04. Implemented in Phase 1 (2026-10-03) and
validated in Phase 2 against the Phase 0 baseline on the old code; see
section 12 and the review report beside this file,
`trailing_reconstruction_without_halfmap_blend_report.md`. It was the last
open design item of the release 4 legacy cleanup.

Validation level: static source inspection only. No source code was changed,
compiled or executed while preparing this note. Line references are to the
working tree of 2026-10-03 and will drift.

This is the single living design record for this refactor. Update it as the
implementation and validation land rather than creating companion plans.

## 1. Background

**Fractional update.** With `update_frac` below 1, each 3D refinement
iteration aligns and reconstructs only a sample of the particles: a fraction
`f` of the pool that has been updated so far.

**Trailing reconstruction** (`trail_rec=yes`) blends that sample with what came
before, so the map does not rest on the sample alone.

**The accumulator chain.** Since the recent rework, the "before" is kept as an
accumulator chain: per state and per even/odd half, the unregularized Fourier
sums and sampling densities of the particles, stored on disk between
iterations at the mass of the full updated pool. Each iteration blends the
current sample's sums with the chain's sums and then restores (density
correction, regularization, deapodization) once. Blending sums before
restoration is what makes the result a single consistent estimator.

**Finished half maps.** "Finished" half maps are maps after that restoration.
Blending them is a different and inconsistent estimator, because restoration
is not linear.

**Gridding and PCG.** The two 3D reconstruction backends. Gridding is the
default; PCG solves the reconstruction iteratively by the preconditioned
conjugate-gradient method. Each keeps its own chain.

## 2. Ruling

The maintainer's ruling of 2026-10-03:

> Trailing needs to be done consistently across using the newly implemented
> approach. No trailing should ever happen on finished halfmaps. If the
> solution to correct weighting and consistent application is one
> reconstruction, I am ok with it.

A first version of this note satisfied the ruling with one extra full
reconstruction of every particle whenever a blend was due and no chain existed.
For a data set of 2.5 million particles, that pass costs about as much as several
ordinary iterations. The revised design below needs no extra reconstruction.
It removes the finished-map blend by starting the chain from the current
sample.

## 3. Current behaviour

### 3.1 Normal iterations (kept)

Notation:

- `N`: the number of particles updated so far in a state;
- `n`: the number in the current sample, so `f = n/N`;
- `u`: the applied update weight (`ufrac_trec` for a single state, otherwise
  `f`);
- `M`: the mass stored in the chain.

The population rule (`population_blend_weights`) scales the current partial
sums by `u/f` and the chain by `(1-u)N/M`. One restoration then yields halves in
which the current sample has weight exactly `u`, and the FSC is estimated from
those halves. When `u >= 0.99`, nothing is blended: the chain is rewritten from
the current sums scaled to full mass.

Code:

- gridding: `simple_commanders_rec_distr.f90`, `restore_state_from_parts` and
  `blend_trailing_accumulators`;
- PCG: `simple_rec3D_pcg_strategy.f90`, `set_chain_blend_weights` and the
  chain branch of the half reduction (around lines 1125-1150).

### 3.2 The start without a chain (to change)

When a blend is due (`u < 0.99`) but no valid chain exists, both backends
already do the consistent thing with the accumulators:

- they write the current sample's sums, scaled by `1/f`, as the new chain at
  full mass;
- they scale the current sums back by `f`, so that this iteration is
  reconstructed from the current sample alone.

They then add a step on finished maps:

- **Gridding.**
  - `restore_eos_and_write_fsc` reads the previous half maps from `vol<state>`
    on the command line (`read_previous_halfmaps`), and uses their FSC to
    regularize the new halves.
  - `trail_restored_halves_if_needed` then writes
    `u * new + (1-u) * previous` for the restored halves, for the inputs of
    the nonuniform filter, and, when `lp` is set, for the merged volume.
- **PCG** (`execute_rec3D_pcg_distributed_master`, with `l_bootstrap` set).
  - The FSC is computed from the previous shipped half maps
    (`load_previous_state_halves`).
  - `blend_bootstrap_half` blends the solved halves, the regularized pair and
    the solvent-weighted pair with the previous pair.
  - The support provenance of the blend is combined from both contributions.
  - `trail_bootstrap_states` reports this per state to
    `filter_pcg_nonuniform_maps`, which only checks the array's size.

So the chain is already built consistently. Only the map of this one
iteration is a blend of finished maps, and that is the part the ruling
removes.

### 3.3 When the start without a chain happens

A blend is due with no chain in these cases:

- **solve3D, once per run**, when trailing switches on. Single-state runs
  switch it on at stage 5 (`TRAILREC_STAGE_SINGLE`), multi-state runs at
  `TRAILREC_STAGE_MULTI`. solve3D's full stage reconstruction (`calc_rec`)
  runs only at the start stage, at state splits and at add-on boundaries. It
  seeds a chain only if the stage it feeds trails, which the start stage does
  not. Once a chain exists, it survives later stage boundaries: a chain from a
  smaller grid is zero-padded to the current one.
- **refine3D with `trail_rec=yes`**, started on a project with update history
  but without a chain in its directory (for example probabilistic modes,
  which keep the update counts, or `startit > 1`).
- **refine3D_auto and other commanders** that switch trailing on mid-run.
- **A chain discarded on validation** (particle population, state layout, a
  larger grid, physical extent, or corrupt or mixed-generation files).

A state that received no sample at all (`f` about 0) also takes this path
today. Its output is then the previous maps, because `u` is about 0.

In a fresh run the first trailing iteration has `f = 1` (the updated pool is
the sample itself), so no blend is due and the chain is simply seeded.

## 4. Proposed design

### 4.1 Principle

Every blend happens in the accumulator domain. When a blend is due but the
chain is missing, the iteration does not blend: it seeds the chain from the
current sample at full mass and ships the map of the current sample alone.
From the next iteration on, every blend uses the chain.

### 4.2 Why this is enough

- **No extra I/O or reconstruction.** The current sample has already been
  reconstructed by the matcher, and the chain seed is written today.
- **No regression in the transition.** The first trailing iteration's map
  rests on the current sample only. That is exactly what every iteration did
  before trailing switched on: with `update_frac` below 1 and trailing off,
  each map is reconstructed from the current sample alone.
- **The chain fills over about `1/f` iterations** (10 at `f = 0.1`), the same
  way it does in a run that starts trailing from a fresh project.
- **The FSC always describes the data that was restored**, and nothing needs
  the previous volumes on the command line.

### 4.3 The trade-off

Compared with today's blend, the first trailing iteration's map is noisier
for one iteration, because the previous maps no longer contribute. Compared
with a full-reconstruction seed (section 5.1), the chain reaches full-dataset
statistics only after about `1/f` iterations instead of at once. The full
reconstruction is still used where one runs anyway: solve3D's add-on
boundaries and state splits already seed the chain when the next stage
trails.

## 5. Alternatives considered

### 5.1 Seed the chain by one full reconstruction (the first version)

When a blend is due and no chain exists, run reconstruct3D with
`trail_seed=yes` over every updated particle instead of assembling the
partials. That gives a full-dataset chain at once, but it reads and grids
every particle again, and it discards the matcher's partials of that
iteration. At 2.5 million particles this costs about as much as several
ordinary iterations, for a gain of one or a few iterations of smoothing.
Rejected on cost.

### 5.2 Seed before matching

Run the full reconstruction at the start of the iteration instead, so that
the first trailing iteration can blend with weight `u` against it. This has
the cost of 5.1. Worse, reconstruct3D writes the same volume, half-map and FSC
files that the iteration is about to use as matching references, so it would
silently replace them. Rejected.

## 6. Changes

### 6.1 Gridding (`simple_commanders_rec_distr.f90`)

**Remove:**

- `read_previous_halfmaps`;
- the start-without-chain branch of `restore_eos_and_write_fsc`; the FSC is
  then always estimated from the restored halves;
- `trail_restored_halves_if_needed` and its call;
- the `vol_prev_even`, `vol_prev_odd` and `vol_merged` arguments of
  `restore_state_from_parts`, and their declarations in `exec_volassemble`;
- the requirement that `vol<state>` be on the command line under
  `trail_rec=yes`.

**Keep:**

- in `blend_trailing_accumulators`, the branch that seeds or rewrites the chain
  at full mass and scales the current sums back by `f`; it now serves both the
  start without a chain and `u >= 0.99`;
- the log line, reworded: "SEEDED FULL-MASS TRAILING CHAIN ...; THIS
  ITERATION USES THE CURRENT SAMPLE ONLY".

**A state with no sample** (`f` below 0.001):

- with a valid chain, carry the chain unchanged and restore this iteration
  from it, as the solve3D_addon cohort path already does;
- without a chain, there is nothing to restore. Carry the previous volume
  forward and skip the state, as volassemble already does for a state without
  partial reconstructions (`determine_dropped_states`).
- *Decision taken in Phase 1:* carrying forward needs the previous state
  volume, its half maps and its FSC file in the run directory under their
  standard names, because the refinement strategy reads them after assembly.
  That is not always so: at the first iteration of a refine3D run the start
  volumes carry other names, or live in another directory when `vol1` is
  given. When they are missing, the run stops with an error that says what is
  missing and to raise `update_frac`. Seeding from the sample instead was
  tried. It ships a map of almost no particles: with one sampled particle, one
  half set is empty, so the FSC is zero and postprocess crashed
  (`scratch/logs/nosample_gridding_seedfallback.log`); PCG stopped with "requires both
  halfsets". Reading the command-line `vol<state>` instead would bring back
  the dependence this plan removes. The case needs a state whose sample is
  below 0.1 % of its updated pool at the start of a run, so in practice
  `update_frac` well below 0.001. The PCG master applies the same rule.

### 6.2 PCG (`simple_rec3D_pcg_strategy.f90`)

**Remove:**

- the `l_bootstrap` branches: the previous FSC pair
  (`load_previous_state_halves`), `blend_bootstrap_half`, and the combination
  of support provenance for blended pairs. The FSC pair is always the current
  base pair.
- the `trail_bootstrap_states` output of
  `execute_rec3D_pcg_distributed_master`;
- the `l_trail_bootstrap` argument of `filter_pcg_nonuniform_maps`, which only
  checks its size;
- the matching arguments in the callers, `simple_refine3D_strategy.f90`
  (`assemble_refine3D_pcg`) and `simple_rec3D_strategy.f90`.

**Keep** the existing chain seeding at full mass and the scale-back by `f` for
the start without a chain. A state with no sample follows the gridding rules.

### 6.3 What does not change

- The population rule, the `u/f` scaling, the chain manifest, the integrity
  checks and the carry-over of chains on `continue=yes`.
- solve3D's full reconstructions and their chain seeding (`trail_seed`), and
  the solve3D_addon paths, which already require a seeded chain.
- The matcher, the strategies and the single-read particle I/O.

## 7. Documents and skills to update

- `.github/skills/simple-frac-update-trailing/SKILL.md` and its
  `references/frac-update-contract.md`: replace the bootstrap paragraph, which
  describes the finished-map blend and the mandatory previous half maps, with
  the start-without-chain rule of section 4.1.
- `doc/policies/importance_sampling_fractional_update_policy.md` and
  `doc/policies/3D/reconstruct3D_pcg_policy.md`: remove the bootstrap blend and
  the lag-one FSC pair.
- `doc/refactoring_notes/planned/release4_legacy_cleanup_inventory.md`: remove
  the item when this lands.

## 8. Tests

- `src/main/image/simple_accum_blend_tester.f90`. Keep
  `run_bootstrap_then_update` ("post-bootstrap effective update weight equals
  realized fraction f"): it checks the full-mass seed, which this design keeps.
  Rename it to describe a chain start rather than a bootstrap. Add one case:
  the seeding iteration's own restored map equals the current sample's map.
- A gridding check that a trailing iteration without a chain neither reads
  nor needs `vol<state>` and ships the current-sample map.
- Run `scripts/check_test_registry.py` if test names change.

## 9. Validation plan (maintainer runs)

On the Oracle Linux test machine, with ground truth: orientation errors
against truth for small noisy sets, FSC for larger ones.

1. **solve3D end to end, single state.** At the stage where trailing switches
   on, expect one "seeded ... current sample only" log line per state and no
   blend log line; accumulator blends follow. Compare orientation error and
   FSC per iteration with the same run before the change. A small dip in that
   one iteration is expected; no lasting difference should be.
2. **The same run with `rec_backend=pcg`.**
3. **refine3D, `trail_rec=yes`, `update_frac=0.2`**, started from a project
   with update history and no chain.
4. **A multi-state run in which one state receives no sample** in some
   iteration, with and without an existing chain.

## 10. Open questions

- **Marking the seeding iteration** in the convergence output
  (`TRAIL_REC_UPDATE_FRACTION`), so the one-iteration change of behaviour is
  visible in run reports.
- **Seeding solve3D's chain earlier.** If the first-iteration dip turns out to
  matter, solve3D could seed the chain at the full reconstruction it already
  runs at its start stage, when a later stage will trail. That costs nothing
  extra, but the chain would then carry early, low-resolution alignments that
  decay only slowly (by `1-u` per iteration). The current code deliberately
  avoids this ("seeding earlier would park stale full-weight alignments in the
  chain"). Decide on evidence from run 1.

## 11. Execution

This plan is carried out as an unattended run on the Oracle Linux test machine
(arun run `trailing_halfmap`), in three phases. Nothing is committed; the
maintainer reviews the complete diff when all phases are done. Every phase
records itself in section 12.

### 11.1 Phase 0: baseline on the current code

No source file changes in this phase.

1. Build the tree as it is (`./compile_debug.sh`, which also runs the fast
   gate) and record how SIMPLE is built on this machine.
2. Find the beta-galactosidase single-state set (5,513 particles) on this
   machine and copy or link what is needed into the run's scratch directory.
   Never write to the original.
3. Run the cases below. For each, record per iteration from the last
   non-trailing stage to the end:
   - FSC 0.5 and 0.143 resolution;
   - mean orientation and shift change;
   - the trailing log lines of volassemble, or of the PCG master for case B;
   - the wall time.

   In every case, identify the iteration where trailing starts without a
   chain: today it logs "USING LEGACY PREVIOUS-HALFMAP BLEND". Record which
   states take that path and how often.

   The cases:
   - **A.** `solve3D` single state on beta-gal, gridding backend, default
     stage plan (trailing switches on at stage 5). Run it twice, so that the
     run-to-run spread is known; a third time if the machine allows.
   - **B.** The same with `rec_backend=pcg`. Twice if time allows, otherwise
     once.
   - **C.** `refine3D` with `trail_rec=yes` and `update_frac=0.2`, five
     iterations, starting from the final project of a case A run. Use a
     setting that keeps the update history, so that the first iteration
     starts trailing without a chain; a probabilistic refine mode keeps it.
     Confirm from the log that the start-without-chain path is taken.
   - **D, optional.** A multi-state case in which one state receives no sample
     in some iteration, if one can be set up at reasonable cost (for example
     with small `update_frac` on a two-state split). If not, record why.
4. If the simulated high-level workflows (`simulated_workflow_1jxy`,
   `simulated_workflow_6vxx`) reach a trailing stage, record their
   ground-truth metrics (orientation and shift errors against the simulated
   truth). If they do not, record that.

Exit criteria:

- the build and the fast gate pass;
- cases A to C are recorded, with the start-without-chain iteration
  identified in each;
- the spread of case A is known from at least two runs;
- all numbers are in the "Phase 0 baseline" subsection of section 12.

### 11.2 Phase 1: implementation

Make the changes of section 6, the tests of section 8 and the document and
skill updates of section 7. Do not change the matcher, the strategies'
iteration logic or the population rule.

Exit criteria:

- the build and the fast gate pass, including the `unit_image` trailing-blend
  tests;
- `python3 scripts/check_test_registry.py .` passes;
- the source contains none of `read_previous_halfmaps`,
  `trail_restored_halves_if_needed`, `blend_bootstrap_half`,
  `load_previous_state_halves`, `trail_bootstrap_states` and
  `l_trail_bootstrap`;
- volassemble no longer reads `vol<state>` from its command line under
  `trail_rec=yes`.

### 11.3 Phase 2: validation and close

1. Repeat cases A to C (and D if it was run) on the Phase 1 build, with the
   same settings and data.
   - At the iteration where trailing starts without a chain, expect the new
     log line ("seeded full-mass trailing chain ... current sample only") and
     no blend.
   - From the next iteration on, expect accumulator blends.
   - For case D, expect the unsampled state's chain to be carried, or its
     volume carried forward when it has no chain.
2. Acceptance: the final FSC 0.143 resolution of each case lies within the
   larger of the Phase 0 spread and 3 % of its Phase 0 value. The resolution
   at the start-without-chain iteration itself may dip; report it, but it is
   not a criterion. Where the simulated workflows reach a trailing stage,
   their ground-truth errors must lie within their Phase 0 values plus the
   same margin.
3. Remove this item from
   `doc/refactoring_notes/planned/release4_legacy_cleanup_inventory.md`, and
   update the trailing-reconstruction paragraph under "What remains" in
   `doc/refactoring_notes/completed/release4_legacy_cleanup_report.md` to say
   that the change is done.
4. Do not regenerate `doc/code_overview/fortran-indexes`; the maintainer does
   that.
5. Move this plan to `doc/refactoring_notes/completed/` with `mv`, and write
   the review report beside it. The report must be readable on its own:
   - say what changed and why, in plain words;
   - give the evidence against the Phase 0 baseline in tables;
   - list any decision taken during the run and where it is recorded;
   - define every SIMPLE term and spell out every acronym on first use;
   - use no references to run-internal labels.

Exit criteria:

- the acceptance of step 2 holds;
- steps 3 to 5 are done.

### File table

| File | Change | Phases |
| --- | --- | --- |
| `src/main/commanders/simple/simple_commanders_rec_distr.f90` | Section 6.1. | 1 |
| `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90` | Section 6.2. | 1 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | Section 6.2: caller arguments only, plus the comment of `carry_over_trail_rec_chains` that named the previous-halfmap bootstrap (Phase 1). | 1 |
| `src/main/strategies/parallelization/simple_rec3D_strategy.f90` | Section 6.2: caller arguments only. | 1 |
| `src/main/volume/simple_halfmap_diagnostics.f90` | Section 6.2: the support-provenance kind `mixed` was written only by the removed PCG bootstrap blend; drop it from the accepted kinds and fix the comment that names the bootstrap as a reader (added in Phase 1). | 1 |
| `src/main/image/simple_accum_blend_tester.f90` | Section 8. | 1 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | Only if a test name changes (section 8). | 1 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | Only if a sub-suite name changes (section 8). | 1 |
| `.github/skills/simple-frac-update-trailing/SKILL.md` | Section 7. | 1 |
| `.github/skills/simple-frac-update-trailing/references/frac-update-contract.md` | Section 7. | 1 |
| `doc/policies/importance_sampling_fractional_update_policy.md` | Section 7. | 1 |
| `doc/policies/3D/reconstruct3D_pcg_policy.md` | Section 7. | 1 |
| `doc/refactoring_notes/planned/release4_legacy_cleanup_inventory.md` | Section 11.3, step 3. | 2 |
| `doc/refactoring_notes/completed/release4_legacy_cleanup_report.md` | Section 11.3, step 3. | 2 |

## 12. Progress

(One row per phase, newest last, written by the run. Phase 0 adds a "Phase 0
baseline" subsection with its tables.)

| Phase | Date | Files changed | Evidence | Exit criteria |
| --- | --- | --- | --- | --- |
| 0: baseline on the current code | 2026-10-03 | This plan only (section 12, and the status line at the top). No source file. | Under `/home/elmlundho/agent_runs/trailing_halfmap/scratch/`: `logs/build_debug_p0.log` (Debug build and fast gate), `logs/build_release_p0.log` (Release build), `logs/case_<run>.log` (every run, each line prefixed by its Unix time), `logs/table_<run>.md` (per-iteration tables), `logs/compact_{A,B,C}.md`, `logs/ctest_sw1jxy_LastTest.log`, `logs/ctest_sw6vxx_r1_LastTest.log`, `logs/ctest_sw6vxx_r2_LastTest.log`; scripts `run_case.sh`, `run_D.sh`, `queue_AB.sh`, `queue_B3.sh`, `extract_iters.py`, `compact_tables.py`. Kept for Phase 2: `scratch/keep/`. | Build and fast gate pass (13 of 13 entries, 3.25 s against a 30 s budget). Cases A, B and C recorded below, the start-without-chain iteration identified in each. Spread of case A known from three runs, of case B from three. Optional case D tried three ways, not obtained (reason below). All numbers are in "Phase 0 baseline". |
| 1: implementation | 2026-10-03 | `simple_commanders_rec_distr.f90`, `simple_rec3D_pcg_strategy.f90`, `simple_refine3D_strategy.f90`, `simple_rec3D_strategy.f90`, `simple_halfmap_diagnostics.f90` (added to the file table in this phase), `simple_accum_blend_tester.f90`, the skill `simple-frac-update-trailing` and its `references/frac-update-contract.md`, `importance_sampling_fractional_update_policy.md`, `reconstruct3D_pcg_policy.md`, and this plan (file table, section 6.1, section 12). | Under `/home/elmlundho/agent_runs/trailing_halfmap/scratch/logs/`: `build_debug_p1e.log` (final Debug build and fast gate), `build_release_p1e.log` (Release build used by Phase 2), `gridding_check_p1e.out` with `gcheck_G0.log`, `gcheck_Ga.log`, `gcheck_Gb.log` (script `scratch/gridding_check.sh`, comparison `scratch/mrc_compare.py`), `nosample_check_p1e.out` with `ns_<backend>_<test>.log` (script `scratch/nosample_check.sh`); the rejected seed fallback for the no-sample corner: `nosample_check_seedfallback.out`, `nosample_gridding_seedfallback.log`, `nosample_pcg_seedfallback.log`. | (1) Build and fast gate pass: 13 of 13 entries, no new compiler warning (11, as on the base); the `unit_image` trailing-blend sub-suite passes 60 of 60 checks (40 on the base; the new chain-start case adds 20). (2) `python3 scripts/check_test_registry.py .` passes (no test or sub-suite name changed). (3) The source contains none of `read_previous_halfmaps`, `trail_restored_halves_if_needed`, `blend_bootstrap_half`, `load_previous_state_halves`, `trail_bootstrap_states`, `l_trail_bootstrap` (`rg` over `src` and `production` finds nothing). (4) volassemble no longer reads `vol<state>` from its command line: no `vol<state>` access is left in the commander, and the gridding check ran it without `vol1`. |
| 2: validation and close | 2026-10-04 | This plan (section 12, status line; moved to `doc/refactoring_notes/completed/`), `doc/refactoring_notes/planned/release4_legacy_cleanup_inventory.md` (item C6 and its ruling removed), `doc/refactoring_notes/completed/release4_legacy_cleanup_report.md` ("What remains" paragraph), and the new report `doc/refactoring_notes/completed/trailing_reconstruction_without_halfmap_blend_report.md`. No source file. | Under `/home/elmlundho/agent_runs/trailing_halfmap/scratch/logs/`: `case_<run>p2.log` and `table_<run>p2.md` for every Phase 2 run, `queue_p2.log`, `queue_c.log`, `case_C2base.log`, `case_C3base.log`, `build_base_p2.log`, `compact_A_p0p2.md`, `compact_B_p0p2.md`, `compact_C_all.md`; scripts `scratch/queue_p2.sh`, `scratch/queue_c.sh`, `scratch/build_base.sh`, and the Phase 0 runners `run_case.sh` and `run_D.sh`, unchanged except that `run_case.sh` can select another build. | (1) Acceptance holds for every case, on the reading recorded below: case A (gridding), all three runs at most 3.98 Å against a bound of 4.02 Å; case B (PCG), at most 3.93 Å against 5.72 Å; case C (refine3D), at most 4.30 Å against 4.33 Å, with the Phase 0 spread of case C completed by two base-code runs. Every start-without-chain iteration logs the new "seeded ... current sample only" line and no blend; from the next iteration on, every iteration blends from the chain. No legacy line appears anywhere. Case D behaves as in Phase 0 (no unsampled state can be produced) with the new line in place of the legacy one. The simulated workflows do not trail, so they set no criterion. (2) Step 3: the inventory item and the report paragraph are updated. (3) Step 4: `doc/code_overview/fortran-indexes` is not regenerated. (4) Step 5: this plan is moved with `mv` and the report is written beside it. |

### Phase 0 baseline

Measured on 2026-10-03 on the Oracle Linux test machine (Intel Core Ultra 9
285, 24 cores, 62 GB memory), on the unchanged base `d32aa9d16` (no extra
launch differences).

**Terms used below.** *FSC* is the Fourier shell correlation between the two
independently refined half maps ("even" and "odd" halves of the particles);
the resolution where it falls to 0.5 and to 0.143 is the usual measure of map
quality, in ångström (Å, smaller is better). The *orientation change* and
*shift change* are the mean change of each sampled particle's projection
direction (degrees) and in-plane shift (pixels) from the previous iteration,
as refine3D reports them; they measure how far the alignment still moves.
The *start without a chain* is the iteration where trailing is due (`u` below
0.99) but no accumulator chain exists. Today volassemble logs it as "SEEDED
FULL-MASS TRAILING CHAIN ...; USING LEGACY PREVIOUS-HALFMAP BLEND THIS
ITERATION". The PCG master has no such line: there the iteration shows "PCG
TRAIL" lines without a "PCG TRAILING BLEND" line. *Gridding* and *PCG* are
the two reconstruction backends (section 1).

**Build.** gfortran, gcc and g++ 15.2.1 from gcc-toolset-15, CMake 3.26.5,
system FFTW 3, OpenMP. `./compile_debug.sh` (Debug, tests on) built the tree
and passed the fast gate. The long runs use a separate Release build of the
same tree (`cmake -DBUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release
-DCMAKE_INSTALL_PREFIX=<build dir>`, `make -j24`, `make install`; script
`scratch/build_release.sh`), as the earlier runs on this machine did, because
Debug is far too slow for solve3D. `~/.bashrc` points `SIMPLE_PATH` at another
checkout, so every run sets `SIMPLE_PATH` and `PATH` to the run's Release build
(see `run_case.sh`).

**Data.** The beta-galactosidase set: 5,513 particles, box 256 pixels,
1.275 Å per pixel, point group D2, mask diameter 180 Å, stack
`/data/bgal/DAT/sumstack.mrc`. Every solve3D run starts from a fresh copy of
the project after 2D clustering, `/data2/bgal_addon/control/bgal_addon.simple`
(read only).

**Decisions taken in this phase.**

- *`nsample=1000` in cases A and B.* With the default `nsample` (10,000),
  solve3D samples more particles than the 5,513 active ones and forces full
  sampling ("FORCING FULL ACTIVE SAMPLING (NO FRACTIONAL OR TRAILING
  UPDATE)"), so trailing never switches on. `nsample=1000` gives a sampled
  fraction of 1000/5513 = 0.181 per iteration and keeps the default stage plan
  otherwise. Trailing then switches on at stage 5, as section 3.3 describes.
- *`startit=2` in case C.* refine3D with `refine=prob` clears the update
  counts at the first iteration of a fresh run (`startit=1`). Then the updated
  pool equals the sample (`f = 1`), no blend is due, and the start without a
  chain is not exercised in a meaningful way (run C0 below). With `startit=2`,
  the update history of the solve3D project is kept (every particle has been
  updated, `N = 5513`), so the first iteration has `f = 0.2` and takes the
  start-without-chain path. `maxits` counts iterations from `startit`, so the
  run still has five iterations (2 to 6).
- *Three runs of case B.* The two first PCG runs ended far apart, so a third
  run was added to show the spread.
- *Nothing else heavy ran during the timed runs.* Cases A, B and C ran one at
  a time; the load average was about 8 to 13 during runs, all of it from the
  run itself (4 parts of 6 threads each).

**Exact command lines** (run from an empty directory under `scratch/runs/`;
`run_case.sh` and `run_D.sh` do exactly this):

- Case A (gridding), runs A1 to A3: copy
  `/data2/bgal_addon/control/bgal_addon.simple` to `bgal.simple`, then
  `simple_exec prg=solve3D projfile=bgal.simple mskdiam=180 pgrp=d2 nsample=1000 nparts=4 nthr=6`
- Case B (PCG), runs B1 to B3: the same with `rec_backend=pcg` appended.
- Case C, run C1: copy `scratch/keep/caseA_final/bgal.simple` (the final
  project of run A1) to `bgal.simple`, then
  `simple_exec prg=refine3D projfile=bgal.simple vol1=/home/elmlundho/agent_runs/trailing_halfmap/scratch/keep/caseA_final/recvol_state01.mrc mskdiam=180 pgrp=d2 refine=prob trail_rec=yes update_frac=0.2 startit=2 maxits=5 nparts=4 nthr=6`
- Case D attempts, runs D1 to D3: copy the same kept project, then
  `simple_exec prg=selection projfile=bgal.simple oritype=ptcl3D infile=/home/elmlundho/agent_runs/trailing_halfmap/scratch/keep/caseD_states.txt mkdir=no`
  (20 particles, drawn with seed 20261003, go to state 2), then
  `simple_exec prg=refine3D projfile=bgal.simple nstates=2 vol1=<V> vol2=<V2> mskdiam=180 pgrp=d2 refine=shc trail_rec=yes update_frac=<UF> startit=2 maxits=5 nparts=4 nthr=6`
  with `V` = `scratch/keep/caseA_final/recvol_state01.mrc`; D1: `V2 = V`,
  `UF = 0.05`; D2: `V2 = V`, `UF = 0.01`; D3: `UF = 0.01`, `V2` =
  `scratch/keep/caseD_vol2_lp40.mrc` (V low-pass filtered to 40 Å by
  `simple_exec prg=filter vol1=V lp=40 smpd=1.275 nthr=8 outvol=...`).

**Summary.** "Final map" is the solve3D final reconstruction of all particles
at full box (for case C, the last iteration). "Last iteration" is the last
refinement iteration. Wall time and peak resident memory are from
`/usr/bin/time -v` around the whole program.

| Run | Wall time | Peak memory | Final map FSC 0.5 / 0.143 (Å) | Last iteration: FSC 0.5 / 0.143 (Å) | Start without a chain |
| --- | --- | --- | --- | --- | --- |
| A1 | 28:16 | 3.1 GiB | 4.35 / 3.93 | 133: 4.53 / 4.03 | iteration 73 (first of stage 5), state 1, `N` 1832, `n` 1000 |
| A2 | 27:51 | 3.1 GiB | 4.41 / 3.89 | 131: 4.60 / 4.03 | iteration 71 (first of stage 5), same counts |
| A3 | 28:13 | 3.1 GiB | 4.41 / 3.89 | 133: 4.53 / 4.03 | iteration 73 (first of stage 5), same counts |
| B1 | 39:43 | 6.3 GiB | 4.41 / 3.93 | 132: 4.47 / 4.03 | iteration 72 (first of stage 5), `F = U = 0.5459` |
| B2 | 39:41 | 6.6 GiB | 7.59 / 5.18 | 133: 7.77 / 5.27 | iteration 73 (first of stage 5), same |
| B3 | 40:03 | 6.0 GiB | 6.53 / 4.29 | 133: 6.66 / 4.41 | iteration 73 (first of stage 5), same |
| C1 | 4:35 | 2.9 GiB | (refine3D) 7.77 / 4.13 | 6: 7.77 / 4.13 | iteration 2 (the first), state 1, `N` 5513, `n` 1103 |

Spread of the final map at FSC 0.143: case A 3.89 to 3.93 Å; case B 3.93 to
5.18 Å.

**What the runs show.**

- In every solve3D run the start without a chain happens exactly once, at the
  first iteration of stage 5, for the one state. The updated pool then holds
  1832 particles and the sample 1000, so `f = u = 0.546`. Every later
  iteration, including all of stages 6 to 8, blends in the accumulator domain
  ("TRAILING ACCUMULATOR BLEND" or "PCG TRAILING BLEND"); no chain is
  discarded at the stage boundaries.
- The FSC printed at the start-without-chain iteration is the FSC of the
  previous iteration's half maps, not of the map shipped (the lag-one pair of
  section 3.2). This is clearest in case C: iteration 2 reports 4.35 / 3.89 Å,
  exactly the solve3D map it started from, and the first accumulator blend
  (iteration 3) reports 7.96 / 4.53 Å, after which the FSC 0.143 recovers to
  4.13 Å by iteration 6 as the chain fills (`f = 0.2`, so about five
  iterations). Phase 2 should compare case C over all five iterations, not at
  iteration 2.
- The log line "USING LEGACY PREVIOUS-HALFMAP BLEND" is printed whenever no
  chain exists, also when `u` is about 1 and no blend is applied (run C0:
  `refine=prob`, `startit=1`, `N = n = 1103`).
- The wide PCG spread arises before trailing switches on. At the end of
  stage 4 (no trailing yet), B2 and B3 were at 8.59 and 7.77 Å (FSC 0.143)
  against 5.53 to 6.16 Å for B1 and all gridding runs, and the gap persisted to
  the end. It is run-to-run variation of de novo solving with the PCG
  backend, not an effect of trailing.
- At the switch-on iteration (stage 5.1) the gridding runs show 8.82 to
  9.07 / 5.53 Å, the iteration after it 9.07 to 9.33 / 6.40 to 6.53 Å. These
  are the values Phase 2 compares the "dip" against.

**Case D (optional): not obtained.** A state that receives no sample cannot
be produced at reasonable cost with refine3D. Sampling is drawn over all
states together, and the state search moves sampled particles between
states. Runs D1 to D3 start with 20 particles in state 2, but every
iteration some sampled particles moved into state 2 (it grew to about 320,
128 and 53 particles), so state 2 always had a sample. This held even with a
40 Å low-pass reference for state 2. The runs do show the multi-state start
without a chain: at the first iteration, both states take the legacy path
(state 2 with `N` = 77, 35 and 24); from the next iteration on, both blend
in the accumulator domain. The code shows what an unsampled state does today.
- Gridding: a state with no sampled particle writes no partial
  reconstruction. volassemble drops it ("HAS NO RECONSTRUCTED PARTICLES;
  DROPPING AND CARRYING PREVIOUS VOLUME FORWARD") before any trailing code
  runs.
- PCG: the state is skipped ("HAS NO SELECTED PARTICLES; SKIPPING") and its
  chain is left untouched.
- The "realized update fraction ~0" warning in `blend_trailing_accumulators`
  needs at least one sampled particle and `n/N` below 0.001.

Logs: `logs/case_D1.log` to `logs/case_D3.log`; the first, misconfigured
attempt is `logs/case_D0_badsetup.log` (the selection wrote into its own
directory, and refine3D relabelled the states uniformly).

**Simulated high-level workflows.** Neither reaches a trailing stage. Their
solve3D runs on about 70 particles with the default `nsample`, so solve3D
forces full sampling. On this base, `simulated_workflow_1jxy` (Release)
passed in 228 s: "SOLVE3D NSAMPLE/ACTIVE FRACTION 142.8571 > 0.9000 ->
FORCING FULL ACTIVE SAMPLING", every iteration sampled 70 of 70, and the
masked FSC against the truth reached 0.143 at 4.46 Å.
`simulated_workflow_6vxx` was run twice. Each run simulates its own
particles: 76 in the first run and 85 in the second, both with full sampling
forced. The first run failed its final-volume correlation check: the masked
FSC against the truth reached 0.143 only at 27.7 Å, a failed de novo solution.
The second run passed: correlation 0.9834, masked FSC 0.143 at 4.71 Å. This
workflow never trails, so the failure has nothing to do with this plan; it
is recorded as run-to-run variation of the unchanged base. It passed three
times on this machine on 2026-10-01, at an older base.
There are no ground-truth trailing metrics to record, so Phase 2 has no
simulated-workflow criterion.

**Kept for Phase 2** in `scratch/keep/`:
- `caseA_final/`: the final project and maps of run A1, the start of case C;
- `caseD_states.txt`;
- `caseD_vol2_lp40.mrc`.

**Per-iteration tables.** Each cell holds, for one run: the iteration
number; a code (`-` no trailing, `L` start without a chain, `B` accumulator
blend); FSC 0.5 / 0.143 (Å); orientation change (degrees) / shift change
(pixels); and the iteration's wall time (seconds). Rows are aligned by stage
and by the position of the iteration within the stage (`5.1` is the first
iteration of stage 5), because the stages do not have the same number of
iterations in every run. The tables start at stage 4, the last stage without
trailing. The last refinement iteration's wall time includes the final
reconstruction.

Case A (gridding), runs A1 to A3:

| stage.k | A1 | A2 | A3 |
|---|---|---|---|
| 4.1 | 58 - 9.60/6.40 22.7/1.24 10 | 58 - 9.89/7.77 20.1/1.69 10 | 58 - 9.60/6.66 22.5/1.48 10 |
| 4.2 | 59 - 9.60/6.04 14.1/0.62 9 | 59 - 9.89/7.59 13.6/0.76 10 | 59 - 9.33/6.16 14.1/0.67 9 |
| 4.3 | 60 - 9.33/6.53 9.6/0.50 10 | 60 - 9.60/7.25 9.4/0.55 9 | 60 - 9.33/6.53 9.6/0.50 9 |
| 4.4 | 61 - 9.33/5.73 7.0/0.42 9 | 61 - 9.33/7.10 7.1/0.46 9 | 61 - 9.33/6.16 7.0/0.46 10 |
| 4.5 | 62 - 9.33/5.73 5.4/0.37 10 | 62 - 9.33/6.16 5.5/0.39 10 | 62 - 9.33/6.28 5.6/0.41 9 |
| 4.6 | 63 - 9.33/5.53 4.5/0.36 9 | 63 - 9.07/6.04 4.5/0.34 9 | 63 - 9.33/6.04 4.4/0.36 10 |
| 4.7 | 64 - 9.33/5.63 3.8/0.32 10 | 64 - 9.33/5.63 3.7/0.29 10 | 64 - 9.33/5.53 3.9/0.33 9 |
| 4.8 | 65 - 9.33/5.63 3.5/0.29 9 | 65 - 9.07/6.16 3.1/0.26 8 | 65 - 9.33/6.04 3.4/0.28 9 |
| 4.9 | 66 - 9.33/5.93 3.1/0.26 9 | 66 - 9.33/5.53 2.6/0.20 8 | 66 - 9.33/5.53 3.1/0.28 10 |
| 4.10 | 67 - 9.33/5.93 2.9/0.25 9 | 67 - 9.07/5.53 2.0/0.16 9 | 67 - 9.07/6.04 2.7/0.23 8 |
| 4.11 | 68 - 9.33/5.63 2.8/0.24 8 | 68 - 9.07/5.53 1.4/0.13 8 | 68 - 9.33/5.53 2.3/0.20 9 |
| 4.12 | 69 - 9.07/5.53 2.5/0.23 9 | 69 - 9.07/5.53 0.9/0.09 9 | 69 - 9.07/6.04 1.7/0.14 8 |
| 4.13 | 70 - 9.07/5.53 2.2/0.20 8 | 70 - 9.33/5.53 0.6/0.06 8 | 70 - 9.07/5.53 1.2/0.10 9 |
| 4.14 | 71 - 9.07/5.53 1.9/0.18 9 |  | 71 - 9.07/6.04 0.8/0.07 8 |
| 4.15 | 72 - 9.07/6.16 1.5/0.13 8 |  | 72 - 9.33/6.04 0.5/0.04 8 |
| 5.1 | 73 L 8.82/5.53 18.1/2.06 13 | 71 L 9.07/5.53 18.1/2.14 13 | 73 L 8.82/5.53 17.9/2.46 13 |
| 5.2 | 74 B 9.07/6.40 17.5/1.47 14 | 72 B 9.33/6.53 17.8/1.45 13 | 74 B 9.07/6.53 18.0/1.78 14 |
| 5.3 | 75 B 8.59/5.53 15.3/1.13 16 | 73 B 9.07/5.63 14.9/1.15 15 | 75 B 8.82/6.40 15.5/1.40 15 |
| 5.4 | 76 B 8.16/5.10 12.0/0.78 15 | 74 B 8.37/5.10 12.3/0.83 16 | 76 B 8.37/5.10 12.9/0.91 14 |
| 5.5 | 77 B 7.96/5.02 10.0/0.61 15 | 75 B 8.16/5.02 10.1/0.63 15 | 77 B 7.96/5.10 10.3/0.64 15 |
| 5.6 | 78 B 7.96/5.02 9.0/0.53 15 | 76 B 7.77/5.02 8.7/0.50 15 | 78 B 7.96/5.10 9.6/0.57 15 |
| 5.7 | 79 B 7.77/5.02 7.7/0.45 15 | 77 B 7.77/5.02 7.5/0.47 15 | 79 B 7.96/5.02 7.8/0.49 15 |
| 5.8 | 80 B 7.77/5.02 7.1/0.44 15 | 78 B 7.77/5.02 7.0/0.42 15 | 80 B 7.77/5.02 7.1/0.42 16 |
| 5.9 | 81 B 7.77/5.02 6.0/0.36 15 | 79 B 7.77/5.02 6.1/0.41 16 | 81 B 7.77/5.10 6.4/0.41 15 |
| 5.10 | 82 B 7.77/5.02 5.5/0.34 15 | 80 B 7.77/5.02 5.7/0.38 14 | 82 B 7.77/5.02 5.8/0.35 15 |
| 5.11 | 83 B 7.77/5.02 4.9/0.31 15 | 81 B 7.77/5.02 5.0/0.34 14 | 83 B 7.59/5.10 5.1/0.34 14 |
| 5.12 | 84 B 7.77/5.02 4.7/0.31 14 | 82 B 7.77/5.02 4.6/0.31 14 | 84 B 7.77/5.10 4.5/0.32 13 |
| 6.1 | 85 B 7.77/5.02 4.0/0.24 14 | 83 B 7.59/4.95 4.0/0.27 14 | 85 B 7.77/5.02 4.0/0.25 15 |
| 6.2 | 86 B 7.59/4.66 3.7/0.40 14 | 84 B 7.59/4.73 3.4/0.37 14 | 86 B 7.59/4.66 3.5/0.37 14 |
| 6.3 | 87 B 7.10/4.66 3.2/0.31 14 | 85 B 7.10/4.66 3.2/0.31 14 | 87 B 7.42/4.66 3.1/0.30 13 |
| 6.4 | 88 B 7.10/4.66 3.0/0.26 13 | 86 B 7.10/4.66 2.8/0.25 13 | 88 B 7.10/4.66 2.8/0.26 14 |
| 6.5 | 89 B 6.80/4.66 2.8/0.23 14 | 87 B 6.95/4.66 2.4/0.21 14 | 89 B 6.80/4.66 2.6/0.21 14 |
| 6.6 | 90 B 6.80/4.66 2.5/0.20 13 | 88 B 6.80/4.66 2.0/0.15 13 | 90 B 6.80/4.66 2.3/0.19 13 |
| 6.7 | 91 B 6.80/4.66 2.2/0.17 14 | 89 B 6.95/4.66 1.7/0.13 14 | 91 B 6.80/4.66 1.9/0.15 14 |
| 6.8 | 92 B 6.66/4.66 2.0/0.15 13 | 90 B 6.80/4.66 1.6/0.13 14 | 92 B 6.80/4.66 1.8/0.15 13 |
| 6.9 | 93 B 6.53/4.66 1.8/0.15 14 | 91 B 6.80/4.66 1.4/0.10 13 | 93 B 6.80/4.66 1.6/0.12 14 |
| 6.10 | 94 B 6.80/4.66 1.5/0.11 13 | 92 B 6.80/4.66 1.2/0.10 14 | 94 B 6.80/4.66 1.4/0.11 13 |
| 6.11 | 95 B 6.80/4.66 1.3/0.10 14 | 93 B 6.80/4.66 1.0/0.10 13 | 95 B 6.80/4.66 1.2/0.10 14 |
| 6.12 | 96 B 6.80/4.66 1.2/0.09 13 | 94 B 6.80/4.66 0.9/0.08 14 | 96 B 6.80/4.66 1.0/0.08 13 |
| 7.1 | 97 B 6.80/4.66 19.7/3.03 21 | 95 B 6.80/4.66 19.1/3.00 21 | 97 B 6.80/4.66 19.6/3.20 21 |
| 7.2 | 98 B 6.66/4.66 18.0/2.18 24 | 96 B 6.66/4.66 18.4/2.21 23 | 98 B 6.80/4.66 18.2/2.42 23 |
| 7.3 | 99 B 6.66/4.66 15.6/1.44 23 | 97 B 6.53/4.66 15.7/1.40 23 | 99 B 6.53/4.66 15.8/1.58 23 |
| 7.4 | 100 B 6.40/4.60 12.7/0.97 23 | 98 B 6.53/4.60 13.0/0.96 23 | 100 B 6.53/4.60 13.1/1.09 23 |
| 7.5 | 101 B 6.40/4.47 11.1/0.81 22 | 99 B 6.53/4.53 11.0/0.80 23 | 101 B 6.53/4.47 10.9/0.79 22 |
| 7.6 | 102 B 6.40/4.41 10.0/0.70 23 | 100 B 6.53/4.53 9.5/0.64 23 | 102 B 6.53/4.41 9.9/0.69 23 |
| 7.7 | 103 B 6.40/4.41 8.3/0.55 23 | 101 B 6.66/4.41 8.2/0.52 22 | 103 B 6.53/4.41 8.2/0.53 23 |
| 7.8 | 104 B 6.53/4.35 7.3/0.49 23 | 102 B 6.53/4.41 7.2/0.48 23 | 104 B 6.53/4.41 7.4/0.51 23 |
| 7.9 | 105 B 6.53/4.35 6.5/0.44 23 | 103 B 6.53/4.41 6.2/0.44 23 | 105 B 6.66/4.41 6.5/0.45 23 |
| 7.10 | 106 B 6.53/4.35 5.7/0.35 23 | 104 B 6.53/4.35 5.5/0.35 22 | 106 B 6.53/4.41 5.6/0.39 23 |
| 7.11 | 107 B 6.53/4.35 5.0/0.34 23 | 105 B 6.53/4.35 4.9/0.31 23 | 107 B 6.53/4.35 4.8/0.29 23 |
| 7.12 | 108 B 6.53/4.35 4.3/0.28 22 | 106 B 6.66/4.35 4.1/0.26 23 | 108 B 6.66/4.35 4.4/0.28 23 |
| 8.1 | 109 B 6.40/4.35 3.8/0.27 22 | 107 B 6.53/4.35 3.8/0.26 22 | 109 B 6.53/4.35 3.8/0.24 23 |
| 8.2 | 110 B 5.26/4.35 16.8/3.14 23 | 108 B 6.53/4.35 3.2/0.23 22 | 110 B 6.28/4.35 18.4/3.26 22 |
| 8.3 | 111 B 5.53/4.24 3.4/0.24 21 | 109 B 6.53/4.24 2.9/0.18 22 | 111 B 6.40/4.24 3.0/0.21 22 |
| 8.4 | 112 B 6.28/4.24 2.9/0.20 22 | 110 B 6.40/4.18 17.6/3.22 22 | 112 B 6.40/4.13 2.8/0.20 22 |
| 8.5 | 113 B 6.40/4.18 2.6/0.18 22 | 111 B 6.40/4.08 2.5/0.17 21 | 113 B 6.40/4.08 2.5/0.17 23 |
| 8.6 | 114 B 6.40/4.18 2.2/0.15 22 | 112 B 6.40/4.08 2.2/0.16 22 | 114 B 6.40/4.08 2.2/0.15 22 |
| 8.7 | 115 B 5.18/4.08 13.6/0.92 22 | 113 B 6.40/4.08 2.0/0.15 22 | 115 B 5.26/4.08 14.4/1.01 23 |
| 8.8 | 116 B 5.18/4.08 2.0/0.17 22 | 114 B 6.53/4.08 1.9/0.14 22 | 116 B 6.28/4.08 1.9/0.15 21 |
| 8.9 | 117 B 6.04/4.08 2.0/0.16 22 | 115 B 5.10/4.08 13.8/0.94 21 | 117 B 6.28/4.08 1.8/0.15 22 |
| 8.10 | 118 B 6.28/4.08 1.8/0.14 22 | 116 B 6.40/4.08 1.6/0.12 22 | 118 B 6.40/4.08 1.7/0.14 22 |
| 8.11 | 119 B 6.40/4.08 1.7/0.14 22 | 117 B 6.40/4.08 1.4/0.12 22 | 119 B 6.40/4.08 1.5/0.13 22 |
| 8.12 | 120 B 5.10/4.08 10.0/0.68 21 | 118 B 6.53/4.08 1.3/0.11 22 | 120 B 4.95/4.08 10.5/0.70 22 |
| 8.13 | 121 B 5.18/4.08 1.6/0.14 23 | 119 B 6.53/4.08 1.2/0.11 21 | 121 B 5.18/4.08 1.3/0.12 22 |
| 8.14 | 122 B 5.18/4.08 1.5/0.13 22 | 120 B 5.02/4.03 10.5/0.68 22 | 122 B 6.28/4.08 1.3/0.12 22 |
| 8.15 | 123 B 6.04/4.08 1.5/0.14 22 | 121 B 6.40/4.08 1.2/0.12 23 | 123 B 6.28/4.08 1.2/0.10 21 |
| 8.16 | 124 B 6.40/4.08 1.4/0.12 22 | 122 B 6.40/4.08 1.1/0.11 22 | 124 B 6.28/4.08 1.2/0.10 22 |
| 8.17 | 125 B 4.80/4.08 8.0/0.59 21 | 123 B 6.40/4.08 1.1/0.12 22 | 125 B 4.66/4.03 8.3/0.60 23 |
| 8.18 | 126 B 4.87/4.03 1.5/0.15 23 | 124 B 6.40/4.08 1.1/0.12 22 | 126 B 4.73/4.03 1.2/0.14 22 |
| 8.19 | 127 B 5.10/4.03 1.3/0.13 22 | 125 B 4.73/4.03 8.3/0.65 22 | 127 B 4.87/4.03 1.1/0.12 22 |
| 8.20 | 128 B 4.87/4.03 1.2/0.12 22 | 126 B 4.73/4.03 1.0/0.12 23 | 128 B 4.66/4.08 1.2/0.14 22 |
| 8.21 | 129 B 5.18/4.03 1.3/0.13 21 | 127 B 4.66/4.03 1.0/0.12 22 | 129 B 4.66/4.03 1.1/0.11 22 |
| 8.22 | 130 B 4.60/4.03 6.8/0.55 22 | 128 B 4.66/4.03 1.1/0.12 22 | 130 B 4.60/4.03 6.5/0.56 23 |
| 8.23 | 131 B 4.60/4.03 1.2/0.13 22 | 129 B 4.60/4.03 1.0/0.12 22 | 131 B 4.60/4.03 1.2/0.14 22 |
| 8.24 | 132 B 4.60/4.03 1.4/0.14 23 | 130 B 4.60/4.03 7.0/0.62 22 | 132 B 4.60/4.03 1.2/0.14 22 |
| 8.25 | 133 B 4.53/4.03 1.4/0.15 56 | 131 B 4.60/4.03 1.0/0.12 54 | 133 B 4.53/4.03 1.0/0.13 55 |

Case B (PCG), runs B1 to B3:

| stage.k | B1 | B2 | B3 |
|---|---|---|---|
| 4.1 | 58 - 9.89/7.42 21.1/1.66 13 | 58 - 15.54/9.89 23.0/2.29 13 | 58 - 12.09/8.37 19.6/1.72 13 |
| 4.2 | 59 - 9.33/6.66 13.9/0.80 13 | 59 - 14.19/9.07 16.2/1.17 12 | 59 - 10.20/7.96 13.5/0.88 13 |
| 4.3 | 60 - 9.33/6.66 9.6/0.54 12 | 60 - 13.06/9.60 11.4/0.86 13 | 60 - 10.20/8.82 9.6/0.65 12 |
| 4.4 | 61 - 9.07/6.04 6.9/0.45 13 | 61 - 10.88/9.33 8.8/0.69 12 | 61 - 10.20/8.82 7.5/0.57 13 |
| 4.5 | 62 - 9.07/5.53 5.4/0.37 12 | 62 - 11.26/9.07 7.2/0.63 13 | 62 - 10.20/7.96 6.1/0.53 13 |
| 4.6 | 63 - 9.07/6.04 4.3/0.31 13 | 63 - 10.53/9.07 6.1/0.56 12 | 63 - 10.20/8.82 5.1/0.45 12 |
| 4.7 | 64 - 9.07/5.63 3.6/0.29 12 | 64 - 10.88/9.07 5.4/0.49 13 | 64 - 10.20/7.77 4.2/0.38 13 |
| 4.8 | 65 - 9.07/6.04 3.2/0.26 12 | 65 - 10.20/9.33 4.7/0.45 13 | 65 - 10.20/8.82 3.6/0.35 12 |
| 4.9 | 66 - 9.07/5.44 2.7/0.22 12 | 66 - 10.20/9.33 4.2/0.42 12 | 66 - 10.20/7.77 2.9/0.28 13 |
| 4.10 | 67 - 9.07/5.63 2.0/0.16 11 | 67 - 10.20/9.07 3.9/0.41 12 | 67 - 10.20/7.96 2.4/0.25 12 |
| 4.11 | 68 - 9.07/5.53 1.4/0.12 12 | 68 - 10.20/9.07 3.6/0.39 12 | 68 - 9.89/7.96 2.1/0.24 11 |
| 4.12 | 69 - 9.07/5.44 1.0/0.09 11 | 69 - 10.20/9.07 3.5/0.36 13 | 69 - 9.89/7.77 1.6/0.20 12 |
| 4.13 | 70 - 9.07/5.53 0.7/0.06 12 | 70 - 10.20/8.59 3.3/0.34 12 | 70 - 9.89/7.77 1.2/0.14 11 |
| 4.14 | 71 - 9.07/5.53 0.4/0.05 11 | 71 - 10.20/8.82 3.3/0.34 11 | 71 - 9.89/7.77 0.9/0.12 12 |
| 4.15 |  | 72 - 10.20/8.59 3.2/0.35 13 | 72 - 9.89/7.77 0.6/0.09 12 |
| 5.1 | 72 L 8.82/5.53 18.0/2.30 24 | 73 L 9.89/7.25 17.6/3.02 23 | 73 L 9.89/7.10 17.3/2.62 23 |
| 5.2 | 73 B 9.07/6.53 17.8/1.58 25 | 74 B 10.20/9.33 17.1/2.43 25 | 74 B 9.89/7.96 18.4/1.89 25 |
| 5.3 | 74 B 8.59/5.18 15.0/1.17 26 | 75 B 9.89/9.07 16.2/1.80 25 | 75 B 9.60/7.77 16.9/1.48 26 |
| 5.4 | 75 B 7.96/5.18 12.5/0.83 28 | 76 B 9.89/8.59 14.7/1.18 26 | 76 B 9.33/7.59 14.6/1.09 26 |
| 5.5 | 76 B 7.77/5.10 10.1/0.62 26 | 77 B 9.89/7.77 12.8/0.93 27 | 77 B 9.07/7.42 12.1/0.82 27 |
| 5.6 | 77 B 7.77/5.02 9.1/0.51 25 | 78 B 9.60/7.77 11.8/0.82 27 | 78 B 9.07/6.80 10.5/0.73 26 |
| 5.7 | 78 B 7.77/5.02 7.8/0.47 27 | 79 B 9.60/7.59 10.2/0.70 26 | 79 B 9.07/6.66 9.2/0.64 26 |
| 5.8 | 79 B 7.77/5.02 6.9/0.41 26 | 80 B 9.60/7.59 8.9/0.60 26 | 80 B 9.07/6.53 8.1/0.55 26 |
| 5.9 | 80 B 7.77/5.02 6.2/0.36 25 | 81 B 9.33/7.59 8.2/0.59 26 | 81 B 9.07/6.53 7.1/0.52 26 |
| 5.10 | 81 B 7.77/5.02 5.5/0.33 25 | 82 B 9.33/8.16 7.8/0.57 27 | 82 B 9.07/6.53 6.4/0.49 26 |
| 5.11 | 82 B 7.59/5.02 4.9/0.31 24 | 83 B 9.60/7.59 7.0/0.49 26 | 83 B 9.07/6.53 5.6/0.46 25 |
| 5.12 | 83 B 7.59/5.02 4.5/0.28 24 | 84 B 9.60/8.16 6.5/0.49 25 | 84 B 9.07/6.28 5.4/0.45 24 |
| 6.1 | 84 B 7.59/5.10 3.9/0.27 19 | 85 B 9.33/7.59 5.6/0.41 18 | 85 B 9.07/6.28 4.5/0.37 18 |
| 6.2 | 85 B 7.59/4.66 3.5/0.36 17 | 86 B 9.33/7.59 5.3/0.46 17 | 86 B 8.82/6.28 3.9/0.37 19 |
| 6.3 | 86 B 7.10/4.66 3.1/0.30 17 | 87 B 9.33/7.59 4.6/0.38 17 | 87 B 8.82/5.93 3.6/0.34 18 |
| 6.4 | 87 B 7.10/4.66 2.7/0.24 18 | 88 B 9.33/7.10 4.8/0.37 16 | 88 B 8.82/5.93 3.2/0.29 18 |
| 6.5 | 88 B 7.10/4.66 2.3/0.20 18 | 89 B 9.33/6.95 4.7/0.36 18 | 89 B 8.59/5.93 2.8/0.26 18 |
| 6.6 | 89 B 6.53/4.66 2.0/0.15 17 | 90 B 9.07/6.95 4.0/0.33 18 | 90 B 8.37/5.53 2.5/0.23 18 |
| 6.7 | 90 B 6.53/4.66 1.7/0.13 17 | 91 B 9.07/6.95 3.9/0.32 17 | 91 B 8.37/5.53 2.3/0.21 18 |
| 6.8 | 91 B 6.40/4.66 1.6/0.13 18 | 92 B 9.07/6.95 3.8/0.29 17 | 92 B 8.37/5.53 2.1/0.21 17 |
| 6.9 | 92 B 6.40/4.66 1.3/0.10 18 | 93 B 9.07/6.95 3.7/0.29 18 | 93 B 8.37/5.53 1.8/0.17 19 |
| 6.10 | 93 B 6.53/4.66 1.2/0.10 18 | 94 B 9.07/6.95 3.4/0.28 17 | 94 B 8.37/5.10 1.7/0.17 18 |
| 6.11 | 94 B 6.40/4.66 1.1/0.09 18 | 95 B 8.82/6.80 3.2/0.27 18 | 95 B 8.37/5.10 1.6/0.16 18 |
| 6.12 | 95 B 6.53/4.66 1.0/0.09 18 | 96 B 8.82/6.66 3.2/0.28 17 | 96 B 8.37/4.95 1.4/0.16 17 |
| 7.1 | 96 B 6.53/4.73 19.6/2.96 30 | 97 B 8.59/6.66 19.0/3.76 30 | 97 B 7.77/5.93 19.1/3.46 30 |
| 7.2 | 97 B 6.40/4.66 18.2/2.37 32 | 98 B 8.16/6.80 17.9/2.80 31 | 98 B 7.77/5.93 17.9/2.47 31 |
| 7.3 | 98 B 6.40/4.66 15.6/1.46 32 | 99 B 8.16/6.66 15.5/1.70 31 | 99 B 7.77/5.53 15.9/1.64 32 |
| 7.4 | 99 B 6.40/4.60 12.9/0.97 32 | 100 B 8.16/6.66 14.5/0.97 31 | 100 B 7.77/5.53 14.0/0.97 31 |
| 7.5 | 100 B 6.40/4.35 11.4/0.86 32 | 101 B 8.37/6.66 12.5/0.83 32 | 101 B 7.77/5.53 12.5/0.88 32 |
| 7.6 | 101 B 6.40/4.35 9.9/0.71 31 | 102 B 8.37/6.66 11.3/0.72 31 | 102 B 7.77/5.53 10.4/0.78 32 |
| 7.7 | 102 B 6.40/4.41 8.3/0.55 32 | 103 B 8.37/6.66 10.2/0.64 30 | 103 B 7.77/5.18 9.0/0.61 31 |
| 7.8 | 103 B 6.40/4.35 7.4/0.46 32 | 104 B 8.59/6.66 9.4/0.61 30 | 104 B 7.77/5.18 8.5/0.60 32 |
| 7.9 | 104 B 6.40/4.35 6.5/0.40 32 | 105 B 8.59/6.66 8.3/0.57 31 | 105 B 7.77/5.18 7.8/0.55 30 |
| 7.10 | 105 B 6.40/4.35 5.8/0.41 32 | 106 B 8.37/6.66 7.8/0.56 31 | 106 B 7.77/5.18 6.8/0.49 32 |
| 7.11 | 106 B 6.40/4.35 4.9/0.30 31 | 107 B 8.37/6.66 6.8/0.48 30 | 107 B 7.77/5.18 6.5/0.49 32 |
| 7.12 | 107 B 6.40/4.35 4.3/0.27 32 | 108 B 8.37/6.66 6.4/0.49 30 | 108 B 7.77/4.95 5.7/0.47 31 |
| 8.1 | 108 B 6.40/4.35 3.8/0.24 31 | 109 B 8.37/6.40 6.0/0.43 31 | 109 B 7.77/4.73 5.3/0.40 31 |
| 8.2 | 109 B 6.40/4.35 3.1/0.20 31 | 110 B 7.77/6.28 17.4/3.70 31 | 110 B 7.59/4.73 17.4/3.28 31 |
| 8.3 | 110 B 5.10/4.29 17.7/3.26 30 | 111 B 7.77/6.28 5.7/0.47 30 | 111 B 7.59/4.73 4.7/0.37 31 |
| 8.4 | 111 B 6.28/4.24 2.7/0.18 31 | 112 B 7.96/6.28 5.1/0.38 30 | 112 B 7.59/4.73 4.4/0.33 32 |
| 8.5 | 112 B 6.28/4.13 2.4/0.15 30 | 113 B 8.16/6.04 5.0/0.37 30 | 113 B 7.59/4.73 4.0/0.31 31 |
| 8.6 | 113 B 6.40/4.13 2.0/0.13 32 | 114 B 8.16/6.28 4.3/0.34 31 | 114 B 7.77/4.73 3.4/0.26 31 |
| 8.7 | 114 B 6.40/4.13 1.8/0.13 30 | 115 B 7.77/5.93 15.6/0.95 31 | 115 B 7.59/4.66 13.6/0.92 30 |
| 8.8 | 115 B 5.10/4.08 13.5/1.01 31 | 116 B 7.77/5.93 4.5/0.33 31 | 116 B 7.59/4.66 3.3/0.28 31 |
| 8.9 | 116 B 6.28/4.08 1.6/0.12 31 | 117 B 7.96/5.93 4.1/0.32 31 | 117 B 7.59/4.66 3.1/0.25 31 |
| 8.10 | 117 B 6.28/4.08 1.5/0.11 30 | 118 B 7.96/5.93 3.9/0.30 31 | 118 B 7.59/4.66 2.7/0.23 31 |
| 8.11 | 118 B 6.40/4.08 1.4/0.12 31 | 119 B 7.96/5.93 3.6/0.30 31 | 119 B 7.59/4.66 2.5/0.21 31 |
| 8.12 | 119 B 6.40/4.08 1.4/0.12 31 | 120 B 7.77/5.53 12.5/0.75 31 | 120 B 7.59/4.60 10.4/0.75 30 |
| 8.13 | 120 B 5.10/4.03 10.0/0.71 30 | 121 B 7.77/5.53 3.3/0.29 31 | 121 B 7.59/4.53 2.3/0.22 31 |
| 8.14 | 121 B 5.10/4.03 1.2/0.12 32 | 122 B 7.77/5.93 3.5/0.28 31 | 122 B 7.59/4.66 2.1/0.21 31 |
| 8.15 | 122 B 5.10/4.03 1.2/0.11 31 | 123 B 7.77/5.93 3.3/0.29 30 | 123 B 7.59/4.66 2.2/0.21 31 |
| 8.16 | 123 B 5.10/4.03 1.2/0.11 31 | 124 B 7.96/5.93 3.2/0.28 31 | 124 B 7.59/4.66 2.0/0.21 31 |
| 8.17 | 124 B 5.10/4.03 1.1/0.10 31 | 125 B 7.77/5.44 10.0/0.61 32 | 125 B 7.25/4.60 8.1/0.61 31 |
| 8.18 | 125 B 4.66/4.08 8.1/0.57 31 | 126 B 7.77/5.35 3.3/0.30 31 | 126 B 7.42/4.47 1.9/0.21 31 |
| 8.19 | 126 B 4.66/4.03 1.1/0.13 31 | 127 B 7.77/5.44 3.2/0.31 30 | 127 B 7.25/4.47 1.8/0.21 31 |
| 8.20 | 127 B 4.66/4.03 1.0/0.12 32 | 128 B 7.77/5.35 3.2/0.30 31 | 128 B 7.10/4.47 1.8/0.21 31 |
| 8.21 | 128 B 4.66/4.03 1.0/0.11 31 | 129 B 7.77/5.44 3.0/0.30 30 | 129 B 7.59/4.47 1.8/0.22 32 |
| 8.22 | 129 B 4.66/4.03 1.1/0.12 31 | 130 B 7.77/5.35 8.3/0.57 31 | 130 B 6.80/4.41 6.6/0.59 31 |
| 8.23 | 130 B 4.53/4.03 6.9/0.58 32 | 131 B 7.77/5.26 3.2/0.32 31 | 131 B 6.80/4.41 1.9/0.24 31 |
| 8.24 | 131 B 4.53/4.03 1.1/0.14 31 | 132 B 7.77/5.26 2.8/0.29 31 | 132 B 6.66/4.41 1.9/0.24 31 |
| 8.25 | 132 B 4.47/4.03 0.9/0.10 64 | 133 B 7.77/5.26 3.1/0.29 62 | 133 B 6.66/4.41 1.9/0.22 61 |

Case C (refine3D, `trail_rec=yes`, `update_frac=0.2`, `startit=2`), run C1; the row label is the iteration's position in the run:

| stage.k | C1 |
|---|---|
| 1 | 2 L 4.35/3.89 6.2/0.69 58 |
| 2 | 3 B 7.96/4.53 5.8/0.56 53 |
| 3 | 4 B 7.77/4.29 5.6/0.47 55 |
| 4 | 5 B 7.77/4.24 5.1/0.42 55 |
| 5 | 6 B 7.77/4.13 4.9/0.42 53 |

### Phase 1 implementation notes

What changed, in plain words:

- **Gridding (volassemble).** When trailing is due and no accumulator chain
  exists, volassemble still seeds the chain from the current sample at full
  mass. It no longer reads the previous half maps, no longer takes the FSC
  (the resolution measure) from them, and no longer blends the restored half
  maps with them. The iteration restores and ships the current sample's map,
  and its FSC describes that map. The log line now reads "SEEDED FULL-MASS
  TRAILING CHAIN, STATE s, REPRESENTED POPULATION N; THIS ITERATION USES THE
  CURRENT SAMPLE ONLY". The optional input FSC of `restore_gridding_pair`
  existed only for the removed path and went with it, as did a benchmark
  timer.
- **PCG master.** The same, for the PCG backend:
  - the previous-pair FSC, the bootstrap blend of the solved, regularized
    and solvent pairs, and the "mixed" support provenance are gone;
  - the `trail_bootstrap_states` output is gone, and with it the matching
    argument of `filter_pcg_nonuniform_maps` and of its two callers;
  - the internal flag is now called `l_sample_seed`;
  - the master logs the same "SEEDED ... CURRENT SAMPLE ONLY" line as
    gridding.
- **A state with no sample** (realized fraction `f` below 0.001, with at least
  one sampled particle). Both backends apply the same rules:
  - with a valid chain, the chain carries unchanged and the iteration is
    restored from it;
  - without a chain, the previous volume is carried forward when the
    directory holds it;
  - otherwise the run stops with an error (decision in section 6.1).

  A state with no sampled particle at all is dropped or skipped before any
  trailing code runs, as before.
- **Tests.** The trailing-blend tester's bootstrap case is renamed to a chain
  start. A new case checks that the chain-start iteration writes a
  full-mass chain holding the current map, and restores with the current
  sample's own sums and density, so the previous map has weight 0.

Checks beyond the fast gate, all on beta-galactosidase with the Phase 1
Release build:

- **Gridding check (section 8).** One refine3D iteration from the kept case A
  project (`refine=prob trail_rec=yes update_frac=0.2 startit=2 maxits=1`)
  starts without a chain. On copies of its directory, volassemble was then run
  by hand twice on the same partial reconstructions, without `vol1` and with
  the chain files removed: once with `trail_rec=yes` and once with
  `trail_rec=no`, which reconstructs the current sample alone. The results:
  - the `trail_rec=yes` run seeded the chain and logged the new line;
  - its even, odd and merged maps equal those of the `trail_rec=no` run to a
    relative 4e-7, and its FSC to 2e-10 (the `1/f` then `f` rounding);
  - they are bitwise identical to what the refine3D iteration itself shipped.
- **No-sample check.** A single state with `update_frac=0.0001` samples 1
  particle of 5513, so `f` is about 2e-4. The same three tests ran on both
  backends, with the same outcome:
  - in a fresh directory (no chain, no previous volume) the run stops with
    the new error;
  - after one ordinary iteration has written a chain, the next iteration
    logs "HAS NO SAMPLE; TRAILING CHAIN CARRIED UNCHANGED" and completes;
  - with the chain removed, it logs "CARRYING PREVIOUS VOLUME FORWARD",
    completes, and leaves the volume, half map and FSC files byte-identical.
- **Review.** A review subagent checked the diff against this plan. It found
  two faults in the first version of the no-sample rule, both fixed:
  - with a single state, the PCG master would stop after carrying the state
    forward;
  - the previous volume may be missing under its standard name.

  It found nothing outside the plan. The matcher, the population rule, the
  solve3D seeding and the strategies' iteration logic are unchanged.

Left as they were (outside this plan, for the maintainer):

- **Resolution fields of a carried-forward state.** In a multi-state gridding
  run such a state gets resolution 0 in the project, exactly like an existing
  dropped state.
- **solve3D_addon with PCG and a near-empty cohort sample.** It still blends
  by the population rule, while gridding carries the cohort chain. This is
  pre-existing.
- **An unused NU evidence source.** The NU filter still accepts the evidence
  source `previous_shipped` (`NU_EVIDENCE_SOURCE_PREV`), which nothing has
  produced since 2026-09-06.

### Phase 2 validation

All runs used the Phase 1 Release build and the Phase 0 command lines unchanged (`scratch/run_case.sh`, `scratch/run_D.sh`), one at a time, on the same machine. The load average was 7 to 15, all of it from the run itself.

**How acceptance was read.** The plan asks that each case's final FSC 0.143 resolution lies within the larger of the Phase 0 spread and 3 % of its Phase 0 value. Here, "Phase 0 value" is the mean over the Phase 0 runs and the spread is their range (maximum minus minimum). A Phase 2 run fails only if it is worse (a larger number in ångström) than that mean by more than the margin; a better resolution never fails. Each run is judged on its own, and so is the Phase 2 mean.

**Case C needed a spread.** Phase 0 ran case C once, so its spread was unknown. The first Phase 2 run ended at 4.30 Å against 4.13 Å, just outside 3 % (4.25 Å). FSC values move in steps of one Fourier shell, here about 1.3 to 2.5 %, and `refine=prob` draws a random sample every iteration. So case C was run twice more on the base code, from an out-of-tree Release build of the unchanged base commit (`git archive`, read-only; `scratch/build_base.sh`, deleted afterwards), alternating with two more runs on the new code. The base-code runs are part of the Phase 0 spread; no tolerance was changed.

| Case | Phase 0 runs: final FSC 0.143 (Å) | Mean / spread / margin | Worse-side bound (Å) | Phase 2 runs: final FSC 0.143 (Å) | Result |
| --- | --- | --- | --- | --- | --- |
| A, solve3D, gridding | 3.93, 3.89, 3.89 | 3.90 / 0.04 / 0.12 (3 %) | 4.02 | 3.80, 3.93, 3.98 (mean 3.90) | pass |
| B, solve3D, PCG | 3.93, 5.18, 4.29 | 4.47 / 1.25 / 1.25 (spread) | 5.72 | 3.93, 3.89, 3.93 (mean 3.92) | pass |
| C, refine3D | 4.13 (Phase 0), 4.24, 4.24 (base code, run in Phase 2) | 4.20 / 0.11 / 0.13 (3 %) | 4.33 | 4.30, 4.24, 4.24 (mean 4.26) | pass |

The FSC 0.5 of the final maps is 4.41 Å in all six solve3D runs of Phase 2 (Phase 0: 4.35 to 4.41 Å for case A, and 4.41, 7.59 and 6.53 Å for case B). In this Phase 2 sample, none of the PCG runs fell into the poor de novo solutions that two of the Phase 0 PCG runs reached before trailing switched on. That is chance, not an effect of this change: the divergence starts in stage 4, before trailing switches on.

**Wall time and memory** (whole program, `/usr/bin/time -v`):

| Run | Phase 0 | Phase 2 |
| --- | --- | --- |
| A1 / A2 / A3 | 28:16 / 27:51 / 28:13, 3.1 GiB | 28:09 / 28:14 / 28:14, 3.0 to 3.1 GiB |
| B1 / B2 / B3 | 39:43 / 39:41 / 40:03, 6.0 to 6.6 GiB | 39:56 / 40:19 / 40:15, 6.0 to 6.5 GiB |
| C1 | 4:35, 2.9 GiB | 4:12, 2.9 GiB (C2, C3: 4:30, 4:38) |

No reconstruction was added and no file of previous maps is read, so time and memory are unchanged within run-to-run variation.

**The start-without-chain iteration.** It is the first iteration of stage 5 in every solve3D run (iteration 72 or 73) and the first iteration (2) of case C. It logs "SEEDED FULL-MASS TRAILING CHAIN, STATE 1, REPRESENTED POPULATION N; THIS ITERATION USES THE CURRENT SAMPLE ONLY", with `N` = 1832 in solve3D and 5513 in case C. The PCG master now logs the same line. No blend line appears in that iteration, and accumulator blends follow in every later iteration.

**The dip at the switch-on iteration** (open question of section 10, recorded, not decided). Phase 0 printed, for that iteration, the FSC of the previous iteration's half maps (the lag-one pair). Phase 2 prints the FSC of the map it ships, the current sample of 1000 particles. The printed FSC 0.143 therefore drops at that iteration:
- solve3D: from 5.53 Å (both backends, every run on the good path) to 7.59 to 8.16 Å;
- case C: from 3.89 Å to 4.53 to 4.60 Å.

The alignment does not suffer. Orientation and shift changes in the next iteration match Phase 0 (about 18° and 1.5 to 1.8 pixels in solve3D). The printed FSC is back at Phase 0 values after one or two iterations (rows 5.2 and 5.3 below), and the final maps are the same. In case C the iterations after the switch-on lie in the base-code range (iteration 3: 4.35 to 4.53 Å on both codes; iteration 4: 4.29 Å on all six runs), and the finals overlap. The dip is cosmetic in the run reports, so marking that iteration in the convergence output (section 10) would help readers, but nothing calls for seeding the chain earlier.

Switch-on rows, case A (columns: Phase 0 runs A1 to A3, then Phase 2 runs A1p2 to A3p2; cell format as in the Phase 0 tables, with `S` for the new start without a chain):

| stage.k | A1 | A2 | A3 | A1p2 | A2p2 | A3p2 |
| 4.1 | 58 - 9.60/6.40 22.7/1.24 10 | 58 - 9.89/7.77 20.1/1.69 10 | 58 - 9.60/6.66 22.5/1.48 10 | 58 - 9.60/7.10 21.0/1.43 9 | 58 - 9.60/6.80 18.8/1.41 10 | 58 - 9.89/6.80 20.8/1.27 10 |
| 4.13 | 70 - 9.07/5.53 2.2/0.20 8 | 70 - 9.33/5.53 0.6/0.06 8 | 70 - 9.07/5.53 1.2/0.10 9 | 70 - 9.07/5.53 0.9/0.09 8 | 70 - 9.07/5.53 1.8/0.17 8 | 70 - 9.07/6.40 2.6/0.25 9 |
| 4.14 | 71 - 9.07/5.53 1.9/0.18 9 |  | 71 - 9.07/6.04 0.8/0.07 8 | 71 - 9.07/6.04 0.6/0.07 9 | 71 - 9.07/5.53 1.4/0.13 9 | 71 - 9.33/6.40 2.5/0.23 8 |
| 4.15 | 72 - 9.07/6.16 1.5/0.13 8 |  | 72 - 9.33/6.04 0.5/0.04 8 |  | 72 - 9.07/5.53 0.9/0.07 8 | 72 - 9.07/5.53 2.3/0.22 9 |
| 5.1 | 73 L 8.82/5.53 18.1/2.06 13 | 71 L 9.07/5.53 18.1/2.14 13 | 73 L 8.82/5.53 17.9/2.46 13 | 72 S 9.89/8.16 18.4/2.24 12 | 73 S 9.89/7.96 18.1/2.17 13 | 73 S 9.89/7.77 18.7/2.13 12 |
| 5.2 | 74 B 9.07/6.40 17.5/1.47 14 | 72 B 9.33/6.53 17.8/1.45 13 | 74 B 9.07/6.53 18.0/1.78 14 | 73 B 9.07/6.80 18.0/1.69 14 | 74 B 9.07/6.53 17.9/1.57 14 | 74 B 9.07/6.95 18.0/1.54 14 |
| 5.3 | 75 B 8.59/5.53 15.3/1.13 16 | 73 B 9.07/5.63 14.9/1.15 15 | 75 B 8.82/6.40 15.5/1.40 15 | 74 B 8.82/5.18 15.7/1.30 16 | 75 B 8.82/5.18 15.1/1.18 15 | 75 B 8.82/6.40 15.1/1.20 15 |
| 5.4 | 76 B 8.16/5.10 12.0/0.78 15 | 74 B 8.37/5.10 12.3/0.83 16 | 76 B 8.37/5.10 12.9/0.91 14 | 75 B 8.37/5.10 12.7/0.84 16 | 76 B 8.37/5.10 12.1/0.85 15 | 76 B 7.96/5.10 12.4/0.83 17 |
| 5.5 | 77 B 7.96/5.02 10.0/0.61 15 | 75 B 8.16/5.02 10.1/0.63 15 | 77 B 7.96/5.10 10.3/0.64 15 | 76 B 8.16/5.02 10.3/0.67 15 | 77 B 7.96/5.10 10.1/0.61 15 | 77 B 7.77/5.10 9.8/0.54 15 |
| 6.1 | 85 B 7.77/5.02 4.0/0.24 14 | 83 B 7.59/4.95 4.0/0.27 14 | 85 B 7.77/5.02 4.0/0.25 15 | 84 B 7.77/5.02 4.0/0.27 14 | 85 B 7.77/4.95 3.9/0.25 14 | 85 B 7.77/5.02 3.8/0.27 14 |

Switch-on rows, case B (PCG), in the same layout:

| stage.k | B1 | B2 | B3 | B1p2 | B2p2 | B3p2 |
| 4.1 | 58 - 9.89/7.42 21.1/1.66 13 | 58 - 15.54/9.89 23.0/2.29 13 | 58 - 12.09/8.37 19.6/1.72 13 | 58 - 10.20/7.77 20.8/1.47 13 | 58 - 12.55/8.82 24.2/1.77 13 | 58 - 9.33/6.80 18.7/1.44 13 |
| 4.13 | 70 - 9.07/5.53 0.7/0.06 12 | 70 - 10.20/8.59 3.3/0.34 12 | 70 - 9.89/7.77 1.2/0.14 11 | 70 - 9.07/6.04 2.0/0.18 12 | 70 - 9.07/5.53 1.3/0.11 12 | 70 - 9.07/5.73 2.7/0.25 12 |
| 4.14 | 71 - 9.07/5.53 0.4/0.05 11 | 71 - 10.20/8.82 3.3/0.34 11 | 71 - 9.89/7.77 0.9/0.12 12 | 71 - 9.07/5.53 1.5/0.13 11 | 71 - 9.07/6.04 0.9/0.09 11 | 71 - 9.33/5.53 2.6/0.24 11 |
| 4.15 |  | 72 - 10.20/8.59 3.2/0.35 13 | 72 - 9.89/7.77 0.6/0.09 12 | 72 - 9.07/5.44 1.0/0.09 12 | 72 - 9.07/6.04 0.7/0.07 12 | 72 - 9.07/5.63 2.3/0.21 12 |
| 5.1 | 72 L 8.82/5.53 18.0/2.30 24 | 73 L 9.89/7.25 17.6/3.02 23 | 73 L 9.89/7.10 17.3/2.62 23 | 73 S 9.89/8.16 18.0/2.41 23 | 73 S 9.89/7.96 17.3/2.24 23 | 73 S 9.60/7.59 18.4/2.17 23 |
| 5.2 | 73 B 9.07/6.53 17.8/1.58 25 | 74 B 10.20/9.33 17.1/2.43 25 | 74 B 9.89/7.96 18.4/1.89 25 | 74 B 9.07/6.53 18.3/1.82 25 | 74 B 9.07/6.40 17.9/1.63 25 | 74 B 9.07/6.40 18.2/1.64 25 |
| 5.3 | 74 B 8.59/5.18 15.0/1.17 26 | 75 B 9.89/9.07 16.2/1.80 25 | 75 B 9.60/7.77 16.9/1.48 26 | 75 B 8.82/5.26 15.3/1.31 26 | 75 B 8.82/5.18 15.2/1.25 27 | 75 B 8.82/5.18 15.4/1.19 26 |
| 5.4 | 75 B 7.96/5.18 12.5/0.83 28 | 76 B 9.89/8.59 14.7/1.18 26 | 76 B 9.33/7.59 14.6/1.09 26 | 76 B 7.96/5.18 12.9/0.88 27 | 76 B 7.96/5.10 12.5/0.89 26 | 76 B 7.77/5.10 12.2/0.82 27 |
| 5.5 | 76 B 7.77/5.10 10.1/0.62 26 | 77 B 9.89/7.77 12.8/0.93 27 | 77 B 9.07/7.42 12.1/0.82 27 | 77 B 7.77/5.10 10.3/0.63 26 | 77 B 7.77/5.02 10.2/0.61 27 | 77 B 7.77/5.10 10.2/0.56 26 |
| 6.1 | 84 B 7.59/5.10 3.9/0.27 19 | 85 B 9.33/7.59 5.6/0.41 18 | 85 B 9.07/6.28 4.5/0.37 18 | 85 B 7.59/5.02 3.7/0.24 18 | 85 B 7.59/4.95 3.7/0.21 19 | 85 B 7.59/4.95 3.9/0.24 18 |

Case C, all six runs (Phase 0 run C1, base-code runs C2base and C3base, Phase 2 runs C1p2 to C3p2; the row label is the iteration's position in the run):

| stage.k | C1 | C2base | C3base | C1p2 | C2p2 | C3p2 |
|---|---|---|---|---|---|---|
| 0.1 | 2 L 4.35/3.89 6.2/0.69 58 | 2 L 4.35/3.89 6.3/0.71 59 | 2 L 4.35/3.89 6.3/0.70 60 | 2 S 8.82/4.60 6.5/0.72 55 | 2 S 8.82/4.53 6.5/0.71 60 | 2 S 8.82/4.60 6.2/0.74 62 |
| 0.2 | 3 B 7.96/4.53 5.8/0.56 53 | 3 B 8.16/4.35 6.0/0.57 52 | 3 B 8.59/4.35 6.0/0.55 53 | 3 B 8.16/4.35 6.1/0.63 52 | 3 B 8.82/4.41 6.1/0.58 54 | 3 B 8.59/4.53 5.9/0.56 54 |
| 0.3 | 4 B 7.77/4.29 5.6/0.47 55 | 4 B 7.77/4.29 5.6/0.49 50 | 4 B 7.77/4.29 5.6/0.49 54 | 4 B 7.77/4.29 5.6/0.49 52 | 4 B 8.16/4.29 5.8/0.48 54 | 4 B 7.96/4.29 5.5/0.48 54 |
| 0.4 | 5 B 7.77/4.24 5.1/0.42 55 | 5 B 7.77/4.24 5.1/0.43 53 | 5 B 7.77/4.29 5.3/0.43 53 | 5 B 7.77/4.29 5.4/0.41 49 | 5 B 7.77/4.24 5.1/0.40 56 | 5 B 7.77/4.24 5.1/0.46 53 |
| 0.5 | 6 B 7.77/4.13 4.9/0.42 53 | 6 B 7.77/4.24 4.8/0.41 51 | 6 B 7.59/4.24 4.9/0.46 47 | 6 B 7.77/4.29 4.9/0.40 44 | 6 B 7.77/4.24 4.8/0.41 52 | 6 B 7.77/4.24 5.0/0.42 54 |

**Case D.** Repeated with the Phase 0 command lines (D1 to D3). Just as in Phase 0, refine3D moves sampled particles into the small second state every iteration, so no state is ever left without a sample. At the first iteration both states log the new seed line, including a state of 28 to 89 particles, and from then on both blend from their chains. The no-sample paths themselves were exercised in Phase 1 by a dedicated check (`scratch/nosample_check.sh`, see the Phase 1 notes).

**Simulated workflows.** They run on 70 to 85 particles and never trail (Phase 0), so they set no Phase 2 criterion and were not rerun.

Per-iteration tables of the Phase 2 runs, from stage 4 on. Case A (gridding):

| stage.k | A1p2 | A2p2 | A3p2 |
|---|---|---|---|
| 4.1 | 58 - 9.60/7.10 21.0/1.43 9 | 58 - 9.60/6.80 18.8/1.41 10 | 58 - 9.89/6.80 20.8/1.27 10 |
| 4.2 | 59 - 9.33/6.40 13.4/0.70 10 | 59 - 9.60/6.66 12.3/0.76 10 | 59 - 9.60/6.40 13.2/0.69 9 |
| 4.3 | 60 - 9.60/6.40 9.2/0.50 9 | 60 - 9.33/6.66 8.9/0.53 9 | 60 - 9.33/6.66 9.4/0.54 10 |
| 4.4 | 61 - 9.33/6.40 6.7/0.42 9 | 61 - 9.60/6.16 6.9/0.46 10 | 61 - 9.07/6.53 6.9/0.47 9 |
| 4.5 | 62 - 9.33/6.53 5.3/0.40 10 | 62 - 9.60/6.66 5.6/0.41 9 | 62 - 9.33/6.66 5.5/0.41 10 |
| 4.6 | 63 - 9.07/6.16 4.6/0.36 9 | 63 - 9.33/5.63 4.6/0.37 9 | 63 - 9.07/6.40 4.6/0.35 9 |
| 4.7 | 64 - 9.33/6.53 3.9/0.31 10 | 64 - 9.33/5.53 3.9/0.32 10 | 64 - 9.07/6.16 4.0/0.33 9 |
| 4.8 | 65 - 9.33/6.66 3.3/0.30 9 | 65 - 9.33/5.53 3.3/0.28 9 | 65 - 9.07/5.63 3.5/0.31 10 |
| 4.9 | 66 - 9.33/6.66 2.9/0.26 8 | 66 - 9.33/5.53 2.9/0.26 9 | 66 - 9.33/5.63 3.3/0.31 9 |
| 4.10 | 67 - 9.07/5.53 2.5/0.21 9 | 67 - 9.33/5.53 2.6/0.24 8 | 67 - 9.07/5.63 3.1/0.30 9 |
| 4.11 | 68 - 9.33/6.04 1.8/0.16 8 | 68 - 9.07/5.53 2.1/0.18 8 | 68 - 9.07/5.53 3.0/0.28 8 |
| 4.12 | 69 - 9.07/6.04 1.3/0.12 9 | 69 - 9.07/6.28 2.0/0.18 9 | 69 - 9.07/5.63 2.8/0.26 8 |
| 4.13 | 70 - 9.07/5.53 0.9/0.09 8 | 70 - 9.07/5.53 1.8/0.17 8 | 70 - 9.07/6.40 2.6/0.25 9 |
| 4.14 | 71 - 9.07/6.04 0.6/0.07 9 | 71 - 9.07/5.53 1.4/0.13 9 | 71 - 9.33/6.40 2.5/0.23 8 |
| 4.15 |  | 72 - 9.07/5.53 0.9/0.07 8 | 72 - 9.07/5.53 2.3/0.22 9 |
| 5.1 | 72 S 9.89/8.16 18.4/2.24 12 | 73 S 9.89/7.96 18.1/2.17 13 | 73 S 9.89/7.77 18.7/2.13 12 |
| 5.2 | 73 B 9.07/6.80 18.0/1.69 14 | 74 B 9.07/6.53 17.9/1.57 14 | 74 B 9.07/6.95 18.0/1.54 14 |
| 5.3 | 74 B 8.82/5.18 15.7/1.30 16 | 75 B 8.82/5.18 15.1/1.18 15 | 75 B 8.82/6.40 15.1/1.20 15 |
| 5.4 | 75 B 8.37/5.10 12.7/0.84 16 | 76 B 8.37/5.10 12.1/0.85 15 | 76 B 7.96/5.10 12.4/0.83 17 |
| 5.5 | 76 B 8.16/5.02 10.3/0.67 15 | 77 B 7.96/5.10 10.1/0.61 15 | 77 B 7.77/5.10 9.8/0.54 15 |
| 5.6 | 77 B 7.96/5.02 9.4/0.56 15 | 78 B 7.77/5.02 8.8/0.55 15 | 78 B 7.77/5.10 9.0/0.52 15 |
| 5.7 | 78 B 7.96/5.02 7.9/0.48 15 | 79 B 7.77/5.02 7.6/0.45 15 | 79 B 7.77/5.02 7.7/0.47 16 |
| 5.8 | 79 B 7.77/5.02 7.1/0.46 16 | 80 B 7.77/5.02 6.7/0.41 15 | 80 B 7.77/5.10 6.7/0.43 15 |
| 5.9 | 80 B 7.77/5.02 6.5/0.39 15 | 81 B 7.77/5.02 6.1/0.39 15 | 81 B 7.77/5.10 5.9/0.37 16 |
| 5.10 | 81 B 7.77/5.02 5.9/0.35 14 | 82 B 7.77/5.02 5.3/0.36 14 | 82 B 7.77/5.02 5.4/0.37 14 |
| 5.11 | 82 B 7.77/5.02 5.1/0.33 14 | 83 B 7.77/5.02 4.9/0.35 15 | 83 B 7.77/5.02 5.0/0.32 13 |
| 5.12 | 83 B 7.77/5.02 4.6/0.30 14 | 84 B 7.77/5.02 4.5/0.30 14 | 84 B 7.77/5.02 4.2/0.30 14 |
| 6.1 | 84 B 7.77/5.02 4.0/0.27 14 | 85 B 7.77/4.95 3.9/0.25 14 | 85 B 7.77/5.02 3.8/0.27 14 |
| 6.2 | 85 B 7.59/4.66 3.7/0.37 14 | 86 B 7.59/4.66 3.4/0.36 14 | 86 B 7.59/4.66 3.5/0.39 14 |
| 6.3 | 86 B 7.59/4.66 3.0/0.28 14 | 87 B 7.59/4.66 3.1/0.30 13 | 87 B 7.42/4.66 3.0/0.30 14 |
| 6.4 | 87 B 7.59/4.66 2.7/0.23 13 | 88 B 7.10/4.66 3.0/0.26 14 | 88 B 7.25/4.66 2.8/0.26 13 |
| 6.5 | 88 B 7.10/4.66 2.5/0.20 14 | 89 B 7.10/4.66 2.7/0.23 13 | 89 B 7.10/4.66 2.5/0.23 13 |
| 6.6 | 89 B 7.10/4.66 2.2/0.18 13 | 90 B 7.10/4.66 2.5/0.21 14 | 90 B 7.10/4.66 2.4/0.20 14 |
| 6.7 | 90 B 7.10/4.66 2.1/0.18 14 | 91 B 6.80/4.66 2.3/0.19 14 | 91 B 7.10/4.66 2.1/0.17 13 |
| 6.8 | 91 B 7.10/4.66 1.8/0.15 14 | 92 B 6.80/4.66 2.2/0.19 13 | 92 B 6.95/4.66 2.0/0.15 14 |
| 6.9 | 92 B 7.10/4.66 1.5/0.12 13 | 93 B 6.80/4.66 2.0/0.16 14 | 93 B 7.10/4.66 1.6/0.12 13 |
| 6.10 | 93 B 7.10/4.66 1.4/0.11 14 | 94 B 6.80/4.66 1.9/0.14 13 | 94 B 7.10/4.66 1.7/0.14 14 |
| 6.11 | 94 B 6.80/4.66 1.2/0.10 13 | 95 B 6.80/4.66 1.7/0.14 13 | 95 B 7.10/4.66 1.4/0.12 14 |
| 6.12 | 95 B 7.59/4.66 1.0/0.09 14 | 96 B 6.80/4.66 1.7/0.14 14 | 96 B 7.10/4.66 1.3/0.11 13 |
| 7.1 | 96 B 7.59/4.66 19.5/3.27 21 | 97 B 6.80/4.66 19.8/3.07 21 | 97 B 6.95/4.66 19.3/3.15 21 |
| 7.2 | 97 B 6.66/4.66 18.5/2.36 23 | 98 B 6.80/4.66 18.0/2.29 24 | 98 B 6.53/4.66 18.4/2.28 23 |
| 7.3 | 98 B 6.53/4.66 15.6/1.56 23 | 99 B 6.66/4.66 15.4/1.41 23 | 99 B 6.53/4.60 15.3/1.50 23 |
| 7.4 | 99 B 6.40/4.60 13.1/0.97 23 | 100 B 6.53/4.60 12.7/1.00 22 | 100 B 6.40/4.47 13.2/1.02 23 |
| 7.5 | 100 B 6.40/4.53 11.4/0.79 23 | 101 B 6.40/4.47 11.4/0.85 23 | 101 B 6.40/4.47 10.8/0.80 23 |
| 7.6 | 101 B 6.53/4.53 9.8/0.69 23 | 102 B 6.40/4.41 10.0/0.69 23 | 102 B 6.40/4.41 9.4/0.69 23 |
| 7.7 | 102 B 6.66/4.41 8.2/0.58 23 | 103 B 6.40/4.41 8.2/0.56 23 | 103 B 6.40/4.41 8.1/0.58 23 |
| 7.8 | 103 B 6.66/4.41 7.9/0.52 24 | 104 B 6.40/4.41 7.4/0.49 23 | 104 B 6.80/4.41 7.1/0.51 23 |
| 7.9 | 104 B 6.66/4.41 6.6/0.45 22 | 105 B 6.53/4.41 6.6/0.44 23 | 105 B 6.80/4.41 6.6/0.44 23 |
| 7.10 | 105 B 6.66/4.35 5.5/0.39 24 | 106 B 6.53/4.35 5.7/0.36 23 | 106 B 6.80/4.35 5.5/0.37 23 |
| 7.11 | 106 B 6.66/4.35 5.1/0.34 22 | 107 B 6.53/4.35 5.1/0.34 24 | 107 B 6.80/4.35 5.0/0.34 22 |
| 7.12 | 107 B 6.80/4.35 4.3/0.29 23 | 108 B 6.53/4.35 4.5/0.30 23 | 108 B 6.80/4.35 4.5/0.29 23 |
| 8.1 | 108 B 6.53/4.35 3.9/0.25 22 | 109 B 6.53/4.35 3.9/0.26 22 | 109 B 6.40/4.35 3.8/0.26 22 |
| 8.2 | 109 B 6.53/4.35 3.5/0.25 22 | 110 B 5.18/4.29 17.9/3.15 22 | 110 B 6.40/4.35 17.7/3.12 22 |
| 8.3 | 110 B 6.40/4.29 17.5/3.20 22 | 111 B 5.53/4.29 3.5/0.25 22 | 111 B 6.40/4.24 3.3/0.23 22 |
| 8.4 | 111 B 6.40/4.29 3.1/0.20 23 | 112 B 6.40/4.24 3.1/0.24 22 | 112 B 6.40/4.13 3.0/0.21 22 |
| 8.5 | 112 B 6.40/4.24 2.7/0.20 21 | 113 B 6.40/4.13 2.8/0.20 22 | 113 B 6.40/4.08 2.7/0.19 22 |
| 8.6 | 113 B 6.53/4.08 2.4/0.17 22 | 114 B 6.40/4.08 2.6/0.20 21 | 114 B 6.40/4.08 2.4/0.16 22 |
| 8.7 | 114 B 6.53/4.08 2.2/0.15 22 | 115 B 5.53/4.08 13.7/1.00 22 | 115 B 5.18/4.08 13.5/0.91 22 |
| 8.8 | 115 B 5.10/4.03 13.6/0.98 22 | 116 B 5.53/4.08 2.4/0.18 23 | 116 B 5.18/4.08 2.4/0.17 22 |
| 8.9 | 116 B 5.53/4.03 1.9/0.14 22 | 117 B 5.53/4.08 2.4/0.19 22 | 117 B 5.18/4.08 2.1/0.15 22 |
| 8.10 | 117 B 6.40/4.03 1.7/0.13 22 | 118 B 5.63/4.08 2.1/0.16 22 | 118 B 6.40/4.08 2.0/0.16 22 |
| 8.11 | 118 B 6.40/4.03 1.5/0.12 22 | 119 B 6.40/4.08 2.0/0.15 22 | 119 B 6.40/4.08 1.8/0.13 22 |
| 8.12 | 119 B 6.40/4.03 1.4/0.12 23 | 120 B 4.95/4.08 10.1/0.70 22 | 120 B 5.10/4.08 10.0/0.69 23 |
| 8.13 | 120 B 5.10/4.03 10.2/0.69 22 | 121 B 5.18/4.08 1.9/0.16 22 | 121 B 5.18/4.08 1.7/0.14 22 |
| 8.14 | 121 B 5.10/4.03 1.3/0.12 21 | 122 B 5.18/4.08 1.9/0.16 21 | 122 B 5.18/4.08 1.7/0.14 21 |
| 8.15 | 122 B 5.10/4.03 1.2/0.12 22 | 123 B 5.18/4.08 1.8/0.15 22 | 123 B 6.40/4.08 1.6/0.14 22 |
| 8.16 | 123 B 6.40/4.03 1.2/0.11 23 | 124 B 5.18/4.08 1.7/0.14 22 | 124 B 6.40/4.08 1.6/0.14 22 |
| 8.17 | 124 B 6.40/4.03 1.1/0.11 22 | 125 B 4.73/4.03 7.6/0.57 22 | 125 B 5.10/4.08 8.0/0.56 22 |
| 8.18 | 125 B 4.66/4.03 7.9/0.60 22 | 126 B 4.73/4.03 1.6/0.15 22 | 126 B 4.95/4.03 1.4/0.14 22 |
| 8.19 | 126 B 4.66/4.03 1.1/0.13 22 | 127 B 4.73/4.03 1.6/0.17 22 | 127 B 5.10/4.03 1.3/0.13 22 |
| 8.20 | 127 B 4.66/4.03 1.1/0.12 22 | 128 B 4.73/4.03 1.7/0.16 22 | 128 B 4.95/4.03 1.5/0.15 22 |
| 8.21 | 128 B 4.66/4.03 1.1/0.12 22 | 129 B 4.66/4.03 1.6/0.16 22 | 129 B 4.95/4.03 1.4/0.14 22 |
| 8.22 | 129 B 4.66/4.03 1.1/0.14 22 | 130 B 4.60/4.03 6.4/0.57 22 | 130 B 4.60/4.03 6.7/0.56 22 |
| 8.23 | 130 B 4.60/4.03 6.3/0.55 22 | 131 B 4.53/4.03 1.6/0.15 22 | 131 B 4.60/4.03 1.3/0.15 22 |
| 8.24 | 131 B 4.53/4.03 1.2/0.13 22 | 132 B 4.53/4.03 1.7/0.15 21 | 132 B 4.60/4.03 1.3/0.15 22 |
| 8.25 | 132 B 4.53/4.03 1.1/0.11 55 | 133 B 4.53/4.03 1.6/0.15 56 | 133 B 4.60/4.03 1.4/0.15 57 |

Case B (PCG):

| stage.k | B1p2 | B2p2 | B3p2 |
|---|---|---|---|
| 4.1 | 58 - 10.20/7.77 20.8/1.47 13 | 58 - 12.55/8.82 24.2/1.77 13 | 58 - 9.33/6.80 18.7/1.44 13 |
| 4.2 | 59 - 9.89/7.59 13.5/0.74 12 | 59 - 10.20/8.82 17.1/1.04 12 | 59 - 9.07/6.28 12.3/0.70 12 |
| 4.3 | 60 - 9.89/7.25 9.6/0.59 13 | 60 - 9.89/7.59 13.1/0.87 13 | 60 - 9.07/6.53 8.7/0.50 13 |
| 4.4 | 61 - 9.60/7.10 7.3/0.52 12 | 61 - 9.60/6.80 9.9/0.79 12 | 61 - 9.07/6.53 6.6/0.45 12 |
| 4.5 | 62 - 9.33/6.95 6.0/0.47 13 | 62 - 9.07/6.40 7.7/0.57 13 | 62 - 9.33/6.53 5.3/0.40 13 |
| 4.6 | 63 - 9.07/6.53 5.1/0.44 12 | 63 - 9.07/6.04 6.0/0.45 12 | 63 - 9.33/6.66 4.3/0.33 12 |
| 4.7 | 64 - 9.07/6.16 4.4/0.37 12 | 64 - 9.33/6.04 5.0/0.36 13 | 64 - 9.33/5.73 3.9/0.31 13 |
| 4.8 | 65 - 9.07/6.16 4.0/0.33 13 | 65 - 9.33/6.04 4.1/0.32 13 | 65 - 9.33/6.53 3.4/0.30 13 |
| 4.9 | 66 - 9.07/5.63 3.5/0.31 12 | 66 - 9.07/6.04 3.6/0.30 12 | 66 - 9.33/6.53 3.2/0.28 12 |
| 4.10 | 67 - 9.07/6.04 2.9/0.24 12 | 67 - 9.07/5.44 3.1/0.26 12 | 67 - 9.07/5.63 3.1/0.28 12 |
| 4.11 | 68 - 9.07/6.53 2.5/0.21 12 | 68 - 9.07/5.53 2.5/0.21 11 | 68 - 9.33/5.73 3.0/0.27 12 |
| 4.12 | 69 - 9.07/6.04 2.2/0.18 11 | 69 - 9.07/5.53 2.0/0.17 12 | 69 - 9.07/6.04 2.8/0.24 11 |
| 4.13 | 70 - 9.07/6.04 2.0/0.18 12 | 70 - 9.07/5.53 1.3/0.11 12 | 70 - 9.07/5.73 2.7/0.25 12 |
| 4.14 | 71 - 9.07/5.53 1.5/0.13 11 | 71 - 9.07/6.04 0.9/0.09 11 | 71 - 9.33/5.53 2.6/0.24 11 |
| 4.15 | 72 - 9.07/5.44 1.0/0.09 12 | 72 - 9.07/6.04 0.7/0.07 12 | 72 - 9.07/5.63 2.3/0.21 12 |
| 5.1 | 73 S 9.89/8.16 18.0/2.41 23 | 73 S 9.89/7.96 17.3/2.24 23 | 73 S 9.60/7.59 18.4/2.17 23 |
| 5.2 | 74 B 9.07/6.53 18.3/1.82 25 | 74 B 9.07/6.40 17.9/1.63 25 | 74 B 9.07/6.40 18.2/1.64 25 |
| 5.3 | 75 B 8.82/5.26 15.3/1.31 26 | 75 B 8.82/5.18 15.2/1.25 27 | 75 B 8.82/5.18 15.4/1.19 26 |
| 5.4 | 76 B 7.96/5.18 12.9/0.88 27 | 76 B 7.96/5.10 12.5/0.89 26 | 76 B 7.77/5.10 12.2/0.82 27 |
| 5.5 | 77 B 7.77/5.10 10.3/0.63 26 | 77 B 7.77/5.02 10.2/0.61 27 | 77 B 7.77/5.10 10.2/0.56 26 |
| 5.6 | 78 B 7.77/5.10 9.3/0.56 27 | 78 B 7.59/5.02 9.1/0.53 27 | 78 B 7.59/5.10 9.2/0.49 26 |
| 5.7 | 79 B 7.77/5.10 8.0/0.48 27 | 79 B 7.59/5.02 7.6/0.43 26 | 79 B 7.59/5.02 7.8/0.46 27 |
| 5.8 | 80 B 7.59/5.10 7.1/0.43 25 | 80 B 7.59/5.02 6.8/0.41 26 | 80 B 7.59/5.02 7.2/0.43 26 |
| 5.9 | 81 B 7.59/5.02 6.3/0.38 26 | 81 B 7.77/5.02 6.1/0.39 25 | 81 B 7.59/5.02 6.1/0.33 26 |
| 5.10 | 82 B 7.77/5.02 5.5/0.33 25 | 82 B 7.77/5.10 5.5/0.32 25 | 82 B 7.59/5.02 5.3/0.33 25 |
| 5.11 | 83 B 7.59/5.02 4.9/0.31 24 | 83 B 7.59/5.10 4.8/0.29 25 | 83 B 7.77/5.10 5.1/0.33 25 |
| 5.12 | 84 B 7.59/5.02 4.5/0.30 25 | 84 B 7.59/5.10 4.5/0.28 24 | 84 B 7.59/5.02 4.5/0.29 25 |
| 6.1 | 85 B 7.59/5.02 3.7/0.24 18 | 85 B 7.59/4.95 3.7/0.21 19 | 85 B 7.59/4.95 3.9/0.24 18 |
| 6.2 | 86 B 7.25/4.66 3.3/0.35 19 | 86 B 7.42/4.66 3.4/0.35 18 | 86 B 7.59/4.66 3.4/0.36 19 |
| 6.3 | 87 B 7.10/4.66 2.9/0.30 18 | 87 B 7.10/4.66 2.9/0.28 18 | 87 B 7.10/4.66 3.2/0.33 18 |
| 6.4 | 88 B 7.10/4.66 2.5/0.22 17 | 88 B 7.10/4.66 2.6/0.23 18 | 88 B 6.80/4.66 2.9/0.26 17 |
| 6.5 | 89 B 6.80/4.66 2.2/0.20 17 | 89 B 6.80/4.66 2.4/0.20 18 | 89 B 6.80/4.66 2.6/0.23 18 |
| 6.6 | 90 B 6.95/4.66 1.9/0.16 17 | 90 B 6.80/4.66 2.0/0.15 18 | 90 B 6.53/4.66 2.5/0.18 18 |
| 6.7 | 91 B 6.40/4.66 1.6/0.14 17 | 91 B 6.80/4.66 1.6/0.13 18 | 91 B 6.53/4.66 2.3/0.19 18 |
| 6.8 | 92 B 6.40/4.66 1.5/0.13 18 | 92 B 6.80/4.66 1.5/0.13 18 | 92 B 6.40/4.66 2.2/0.17 18 |
| 6.9 | 93 B 6.40/4.66 1.3/0.10 17 | 93 B 6.80/4.66 1.3/0.10 18 | 93 B 6.40/4.66 1.9/0.15 17 |
| 6.10 | 94 B 6.40/4.66 1.2/0.10 17 | 94 B 6.80/4.66 1.1/0.09 18 | 94 B 6.40/4.66 1.8/0.14 17 |
| 6.11 | 95 B 6.80/4.66 1.0/0.09 18 | 95 B 6.80/4.66 1.0/0.10 18 | 95 B 6.40/4.66 1.5/0.12 19 |
| 6.12 | 96 B 6.80/4.66 0.9/0.08 17 | 96 B 6.80/4.66 0.9/0.08 18 | 96 B 6.53/4.66 1.3/0.10 17 |
| 7.1 | 97 B 6.80/4.73 19.2/3.29 31 | 97 B 6.66/4.73 19.0/3.11 30 | 97 B 6.53/4.66 18.9/3.04 30 |
| 7.2 | 98 B 6.40/4.66 17.4/2.40 32 | 98 B 6.40/4.66 18.2/2.27 32 | 98 B 6.40/4.66 17.5/2.32 31 |
| 7.3 | 99 B 6.40/4.60 14.9/1.53 31 | 99 B 6.40/4.66 15.4/1.53 32 | 99 B 6.40/4.60 14.7/1.49 32 |
| 7.4 | 100 B 6.28/4.60 13.0/1.03 32 | 100 B 6.40/4.60 13.2/0.99 32 | 100 B 6.28/4.60 12.6/1.03 32 |
| 7.5 | 101 B 6.28/4.41 11.2/0.85 31 | 101 B 6.40/4.60 11.0/0.79 32 | 101 B 6.40/4.53 11.1/0.81 32 |
| 7.6 | 102 B 6.40/4.60 9.7/0.64 32 | 102 B 6.40/4.47 9.6/0.64 32 | 102 B 6.40/4.53 9.4/0.66 31 |
| 7.7 | 103 B 6.40/4.41 8.3/0.58 32 | 103 B 6.40/4.41 8.0/0.57 31 | 103 B 6.40/4.35 7.9/0.55 32 |
| 7.8 | 104 B 6.40/4.41 7.5/0.49 31 | 104 B 6.40/4.41 7.3/0.51 32 | 104 B 6.40/4.47 7.4/0.51 31 |
| 7.9 | 105 B 6.40/4.41 6.4/0.41 32 | 105 B 6.40/4.35 6.4/0.44 32 | 105 B 6.40/4.47 6.4/0.45 32 |
| 7.10 | 106 B 6.40/4.41 5.8/0.38 32 | 106 B 6.40/4.35 5.7/0.37 31 | 106 B 6.40/4.47 5.8/0.39 32 |
| 7.11 | 107 B 6.40/4.41 4.9/0.30 31 | 107 B 6.40/4.35 5.0/0.31 32 | 107 B 6.53/4.35 5.1/0.35 32 |
| 7.12 | 108 B 6.40/4.41 4.3/0.28 32 | 108 B 6.53/4.35 4.2/0.26 31 | 108 B 6.53/4.35 4.4/0.32 31 |
| 8.1 | 109 B 6.28/4.41 3.7/0.21 31 | 109 B 6.40/4.35 3.9/0.23 31 | 109 B 6.40/4.41 4.2/0.30 31 |
| 8.2 | 110 B 5.18/4.35 16.7/3.29 32 | 110 B 5.18/4.35 16.9/3.15 31 | 110 B 5.53/4.35 17.2/3.10 32 |
| 8.3 | 111 B 5.18/4.18 3.2/0.21 31 | 111 B 5.53/4.29 3.0/0.21 32 | 111 B 5.53/4.24 3.7/0.26 31 |
| 8.4 | 112 B 5.18/4.13 2.9/0.18 31 | 112 B 6.28/4.24 2.8/0.17 31 | 112 B 5.53/4.24 3.4/0.24 31 |
| 8.5 | 113 B 5.18/4.13 2.3/0.15 31 | 113 B 6.28/4.08 2.5/0.17 30 | 113 B 5.53/4.08 3.2/0.22 31 |
| 8.6 | 114 B 6.28/4.13 2.1/0.15 31 | 114 B 6.40/4.08 2.2/0.14 31 | 114 B 5.53/4.03 2.9/0.20 30 |
| 8.7 | 115 B 5.10/4.08 13.3/0.95 30 | 115 B 5.18/4.08 13.5/1.02 31 | 115 B 5.02/4.03 13.0/0.99 32 |
| 8.8 | 116 B 5.18/4.08 1.9/0.14 32 | 116 B 5.53/4.08 1.9/0.14 30 | 116 B 5.02/4.03 2.6/0.19 31 |
| 8.9 | 117 B 5.18/4.08 1.7/0.13 31 | 117 B 5.53/4.08 1.8/0.13 31 | 117 B 5.18/4.03 2.4/0.19 31 |
| 8.10 | 118 B 6.28/4.08 1.6/0.12 31 | 118 B 5.63/4.08 1.5/0.12 31 | 118 B 5.53/4.03 2.2/0.17 32 |
| 8.11 | 119 B 6.28/4.08 1.4/0.11 31 | 119 B 6.40/4.08 1.4/0.11 32 | 119 B 5.53/4.03 2.1/0.15 31 |
| 8.12 | 120 B 4.87/4.08 10.0/0.68 31 | 120 B 5.10/4.08 9.7/0.67 31 | 120 B 4.95/4.03 9.5/0.66 30 |
| 8.13 | 121 B 5.10/4.08 1.3/0.11 30 | 121 B 5.18/4.03 1.4/0.12 31 | 121 B 5.02/4.03 2.0/0.16 32 |
| 8.14 | 122 B 5.18/4.08 1.3/0.11 31 | 122 B 5.53/4.08 1.2/0.11 31 | 122 B 5.02/4.08 1.9/0.15 31 |
| 8.15 | 123 B 5.10/4.08 1.1/0.10 31 | 123 B 5.53/4.08 1.3/0.11 31 | 123 B 5.02/4.08 1.8/0.15 31 |
| 8.16 | 124 B 6.28/4.08 1.2/0.12 30 | 124 B 5.53/4.03 1.2/0.11 31 | 124 B 5.02/4.03 1.9/0.16 30 |
| 8.17 | 125 B 4.66/4.08 7.9/0.57 32 | 125 B 4.73/4.08 7.6/0.57 31 | 125 B 4.95/4.08 7.7/0.56 32 |
| 8.18 | 126 B 4.66/4.08 1.1/0.11 31 | 126 B 4.73/4.08 1.1/0.11 31 | 126 B 4.95/4.08 1.7/0.15 31 |
| 8.19 | 127 B 4.66/4.08 1.1/0.10 31 | 127 B 4.73/4.08 1.1/0.11 31 | 127 B 4.95/4.08 1.7/0.14 30 |
| 8.20 | 128 B 4.66/4.08 1.1/0.11 30 | 128 B 4.66/4.08 1.0/0.11 31 | 128 B 4.95/4.03 1.7/0.17 31 |
| 8.21 | 129 B 4.66/4.08 1.0/0.12 31 | 129 B 4.66/4.08 1.0/0.11 31 | 129 B 4.73/4.03 1.7/0.16 32 |
| 8.22 | 130 B 4.53/4.08 6.5/0.54 31 | 130 B 4.53/4.08 6.0/0.48 31 | 130 B 4.53/4.03 6.3/0.50 31 |
| 8.23 | 131 B 4.53/4.08 1.0/0.11 31 | 131 B 4.53/4.03 1.1/0.13 31 | 131 B 4.53/4.08 1.6/0.14 31 |
| 8.24 | 132 B 4.53/4.08 1.0/0.11 31 | 132 B 4.53/4.03 1.0/0.12 32 | 132 B 4.53/4.08 1.6/0.15 31 |
| 8.25 | 133 B 4.53/4.08 1.0/0.12 63 | 133 B 4.53/4.08 1.0/0.11 64 | 133 B 4.53/4.03 1.7/0.16 65 |
