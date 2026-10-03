# Trailing reconstruction without the finished-halfmap blend

Date: 2026-10-03 (revised the same day).

Status: proposed; implementation not started. It is the last open
design item of the release 4 legacy cleanup.

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
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | Section 6.2: caller arguments only. | 1 |
| `src/main/strategies/parallelization/simple_rec3D_strategy.f90` | Section 6.2: caller arguments only. | 1 |
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
