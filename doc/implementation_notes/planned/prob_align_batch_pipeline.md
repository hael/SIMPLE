# Probabilistic alignment as a batch pipeline: per-particle draws, keyed random streams, no global table

Planned 2026-11-04, on master `a723f4919`. **Nothing here is implemented.** Ruled
by Hans the same day: this refactoring is done ("we have reached the point of no
return; this refactoring is obvious"). It is the design and the plan, with the
gates, for the architecture review of 2026-11-04 in
`speed_architecture_proposals.md` (its "Architecture review" section is the
frame; this note is the execution). It supersedes the Part 2 rows of that note on
the worker fold, the streaming assignment, the read-ahead thread and the 2D table
pre-step, and it closes the worker fold that was implemented and dropped on
2026-11-04 (that note's status of 2026-11-03).

## 1. What is proposed, in one paragraph

Since 2026-10-24 every probabilistic assignment is one draw from one particle's
own row; no population balance, no other particle's row, nothing enters a draw
that is not fixed before the iteration's particle loop. The global table, the
part files, the master's assignment step and the second worker pass exist only
to serialise those draws through one random stream. Give every draw its own
stream, keyed on the iteration, the particle and the draw, and the table becomes
unnecessary: one pipeline prepares a batch of particles, scores its rows, draws
each particle's pose and extras, polishes the drawn pose, computes the noise
residual and commits the decision, then releases the batch. The same pipeline
runs in shared memory on the caller's `(params, build)` and in every distributed
part; the master draws the sample, waits, merges the decisions and assembles, as
it does today. The stochastic trajectory changes when the draws move to keyed
streams. The draws then become independent of thread and partition placement.
Complete-run byte identity additionally needs a numerical execution/reduction
contract; it does not follow from keyed RNG alone (review, section 11).

## 2. Where the code stands

### 2.1 The two-phase structure

A probabilistic iteration runs each particle through two programs or phases:

| phase | 3D dense (`refine=prob`) | 3D neighbourhood (`prob_neigh`) | 2D dense (`prob`) | 2D sparse (`prob_neigh`) |
| --- | --- | --- | --- | --- |
| table worker | `exec_prob_tab` → `eul_prob_tab%fill_tab_range` (`process_particle_with_map` or `refine_best_refs`) | `exec_prob_tab_neigh` → `eul_prob_tab_neigh%fill_tab_range` | `exec_prob_tab2D` → `eul_prob_tab2D%fill_tab_range` | `exec_prob_tab2D` → `fill_tab_prob_neigh_range` |
| part file | dense rows, format -4 | neighbourhood candidates | dense rows, format -7 | candidate stream, format -6 |
| master | `prob_align_exec`: `read_tabs_to_glob`, `ref_assign`, `write_assignment` (`write_insert_table` under `topk > 1`) | `prob_align_neigh_exec`, `ref_assign_neigh` | `prob_align2D_exec`, `ref_assign` | `ref_assign_sparse` |
| search worker | `refine3D_exec` → `strategy3D_prob%srch_prob`: `read_assignment`, extras, continuous polish, commit; sigma2 residual; partial reconstruction | the same | `refine2D_exec` → `strategy2D_prob%srch_prob`: polish, shift state, commit; class sums | the same |

The commanders are in `src/main/commanders/simple/simple_commanders_prob.f90`
(550 lines), the masters in
`src/main/strategies/parallelization/simple_prob_align_strategy.f90` (522), the
tables in `src/main/strategies/search/probabilistic/` (`simple_eul_prob_tab.f90`
1,733, `_2D` 1,901, `_neigh` 1,743, `_utils` 1,592), the matchers in
`src/main/strategies/search/simple_strategy3D_matcher.f90` (641) and
`simple_strategy2D_matcher.f90` (516).

Shared memory runs the same two phases through the file system: the strategy
writes the particle field, runs the `prob_align` commander (which builds a second
`parameters` and builder from the file, reseeds the generator, fills the table in
process, writes the part file, reads it back, assigns and writes the assignment
file), then rereads the field (`simple_refine3D_strategy.f90` lines 1267-1276;
the 2D strategy's `run_prob_align2D`). The matcher then reads the assignment
file. The layering policy lists this as the one accepted exception.

### 2.2 What couples particles, and what does not

Fixed before the particle loop, the same for every particle (the iteration's
immutable model):

- the sample (`sample4update_cnt`, `_rnd`, `_fillin`, drawn by the master or
  the shared-memory strategy, written to the project, reproduced by
  `sample4update_reprod`);
- the references (the reprojection model, or the class averages), the halves,
  the band, the matching low-pass limit, the shift window and map geometry;
- the committed sigma2 spectra;
- the shift prior's width (`shift_prior_sigma2`, from the previous iteration's
  sampled `shincarg`), the particle temperatures' formula;
- the per-state FSC cones of the weighted insertion (`insert_cone_from_fsc`);
- the 2D class populations of the previous assignment (`prep_strategy2D_glob`);
- `prob_beta`, the extremal fraction, the refined-reference budget
  (`calc_num2sample`).

Row-local, reading only the particle's own row and the model above:

- the assignment draw: `draw_posterior_row` over the row, in
  `simple_eul_prob_tab.f90` line 1115 (dense 3D: the closest reference's state,
  then the direction within its block), `simple_eul_prob_tab2D.f90` lines 1083
  and 1223, `simple_eul_prob_tab_neigh.f90` line 828;
- the extras (`select_extras`): the row and the state's cone;
- the fallbacks: a row without a valid distance, a particle without a
  temperature;
- the polish (`refine_assignment_continuously`,
  `refine_selected_continuously`): deterministic L-BFGS from the drawn pose;
- the sigma2 residual of the committed pose; the 2D shift state
  (`record_shift_state`).

Nothing else. The multi-state proposal (`multistate_posterior_state_draw.md`,
the joint draw of state and direction over the whole row) is row-local as well.
The `n_eff` and row-mode readouts (`prob_neff_readout`, `POSTERIOR DRAW N_EFF`,
`% at the row mode`) and the disk-set line (`print_disk_set`) are sums over
particles, not inputs to any draw.

### 2.3 The random draws

Every draw on the probabilistic path, with where it runs today:

| draw | where | stream today |
| --- | --- | --- |
| the sample | master or shared-memory strategy, before the loop | the process's intrinsic stream |
| the in-plane draw per scored reference | the fill (`draw_posterior_crvec`, `simple_polarft_corr.f90` lines 672-740) | the fill thread's stream (gfortran keeps one per thread); each fill calls `seed_rnd` at its start (`simple_eul_prob_tab.f90` 302, `_2D` 282 and 528, `_neigh` 197 and 451) |
| the disk draw of the polar shift map | the fill (`shiftmap_draw`, `simple_polarft_shiftmap.f90` line 607) | the fill thread's |
| the candidate picks of the sparse fills | the fill (`ran_tabu` per thread, `_2D` 539, `_neigh` 205) | the fill thread's |
| the assignment draw | the master, in particle order | the master's |

The fills run under a dynamic schedule, so which draws a particle meets depends
on the thread it lands on: byte identity holds at one thread only (the policies
say so). The classic strategies' samplers (`sample_bounded_dists`,
`sample_power_dists`, `shcloc`, the 2D search order) are not on this path.

### 2.4 Measured on the current code (dell, 2026-11-03/04)

trpv1, 35,645 particles, c2, four parts of six threads, `refine=prob`:

| | 3D, crop 160, 2,500 directions | 2D, crop 120, 100 classes | 2D, box 256, 100 classes |
| --- | --- | --- | --- |
| table job per part | 264-333 s (fill 264-274 s, preparation 3.0 s, read 0.5 s, references 0.4 s) | 24 s (fill 20 s, preparation 2.6 s, read 0.4 s) | 62 s (fill 56 s, preparation 5.1 s) |
| master assignment | 0.05 s (reads 0.02 s) | 0.05 s | 0.05 s |
| search job per part | 21-23 s: preparation 3.0, references 0.4, read 0.55, polish 8-9.4, partial reconstruction 9.3 | 4.2-4.7 s: preparation 1.7-2.0, read 0.5, polish 0.4, class sums 1.3 | 14.8 s: preparation 5.4, read 0.55, class sums 7.5 |
| iteration | 309-362 s | 30 s | 78 s |
| duplicated per part (preparation, read, references, launch) | about 4 s, 1.3 % | about 2.4 s, 8 % | about 6.2 s, 8 % |
| master table | 713 MB (`dists` + `inpls`, 8 bytes per entry); part files 199 MB each | | |

The worker fold of 2026-11-03 (two programs in one job, a marker-file protocol,
the table's buffer handed to the search) recovered the read and the full-box
preparation in cropped 2D (search phase 5.1 to 4.1 s) and nothing in 3D, lost
it in the table phase's noise, and needed its own scheduling path. Dropped.

At 750,000 particles and 2,500 directions the dense master table is 15 GB per
state. That, not the 1-8 % of duplicated preparation, is what the two-phase
structure costs.

## 3. The contract

The execution contract of the architecture review holds; what follows makes it
concrete for this code.

**Lifetimes.**

| lifetime | holds |
| --- | --- |
| the run (`refine3D`, `refine2D` command; a distributed part's program) | partition identity, the project owner's project, control metadata; nothing of the alignment |
| the iteration's model, immutable from before the particle loop to after it | references (polar model or class averages, both halves), band and grid, shift-map geometry, committed sigma2 spectra, sample, incoming orientations/shifts/widths, temperature formula, `prob_beta`, shift-prior width, refined-reference budgets, insertion cones and 2D class populations; particle temperatures are computed from the prepared particles |
| the alignment phase's mutable workspace | compatible FFT plans and per-thread optimizer/map/numerical scratch; reused between batches and released before 3D reconstruction |
| one batch (capacity bounded by the memory budget, initially capped by `nthr x BATCHTHRSZ`) | the canonical observations, prepared polar particles and CTF terms, complete candidate rows, particle temperatures, drawn poses/extras, polished poses and residuals; released when its decisions are committed |
| across the 3D alignment/reconstruction barrier | the final poses, shift states, the frozen insertion poses and weights, the sigma2 contributions, the record buffer of cropped observations (as now); no alignment calculator or reference memo |

**The batch path.**

```text
batch of sampled particles (the part's share, in order)
  -> read, or take from the record buffer      (ptcl_read_session, as now)
  -> prepare the alignment batch               (build_batch_particles3D / 2D, as now)
  -> fill the complete rows                    (the route's fill on the batch, as now)
  -> draw each particle's pose and extras      (the row selector; the route's rule)
  -> polish the drawn pose                     (the strategies' continuous refinement, as now)
  -> sigma2 residual of the committed pose     (as now)
  -> commit: orientation, shift state, extras  (as now)
  -> 2D: the batch's class sums; 3D: the record for the reconstruction (as now)
  -> release the rows
```

The row never leaves the process, so the sparse routes keep their candidate
lists in memory as the row (an adapter over the candidate stream's records) and
the dense routes keep the batch's `batchsz x nrefs` entries (at 2,500
references and 300 particles, 6 MB for distances and in-plane indices alone).
References include all active states; refined entries, diagnostics, scratch and
observations are additional. The memory contract and budget are in section 11.

**One selector.** `draw_posterior_row` stays the primitive. The route adapters
give it the row (the dense block of the drawn state, the neighbourhood's pooled
candidates, the 2D class row or candidate list) and take back the draw, as
`ref_assign`, `ref_assign_neigh`, `ref_assign_sparse` and the 2D `ref_assign`
do today; those four become the adapters' selection step, called per particle
inside the batch instead of over the global table. The multi-state joint draw
lands in the same place.

**Consumers decide the preparation.** A batch prepares what its consumers need:
scoring and polish need the polar particles; the sigma2 residual needs the
committed pose against the reference; the 3D reconstruction needs the cropped
record; the 2D class sums need the cropped (or raw) particle. Nothing is prepared
twice, and nothing is prepared for a consumer that is off (`volrec=no`,
`objfun=cc`).

**What stays as it is.** The master's sample draw and project ownership
(`layering_and_project_ownership_policy.md`); `merge_algndocs`; volume assembly
and the sigma2 commit; the 3D partial reconstruction after the loop from the
record buffer (`separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`
stays a hard boundary); the 2D class-sum order (batch, class, half); the
`pose_cont` Cartesian polish pass, which stays a pass of its own; the classic
modes (`shc`, `snhc_smpl`, ...), untouched.

## 4. The random-stream contract

**The generator.** A keyed, stateless draw: `rnd_keyed(key) -> real in [0,1)`,
with the key the tuple (`SIMPLE_SEED` or the run's seed, iteration, particle
index, phase, index), hashed by a 64-bit mixer (SplitMix64 or equivalent; the
mixer is fixed by the contract, since changing it changes results) and mapped
to the unit interval. One hash per draw; the kernels dwarf it. Thread-safe with
no state, so a particle draws the same numbers on any thread, in any batch, in
any part.

**The sites and their keys.** This table identifies draw families, not a complete
key encoding. Section 11 requires the command/stage/pass namespace, canonical
particle/reference identity, repeated-evaluation identity, generator version
and restart contract to be fixed before phase 1.

| draw | phase | index |
| --- | --- | --- |
| the sample (`sample4update_cnt`, `_rnd`, `_fillin`) | `SAMPLE` | the draw's ordinal in the sampler; particle index 0 |
| the in-plane draw per scored reference (`draw_posterior_crvec`) | `INPL` | the reference |
| the disk draw (`shiftmap_draw`) | `DISK` | the screen's ordinal for the particle |
| the sparse fills' candidate picks (`ran_tabu`) | `CAND` | the pick's ordinal for the particle |
| the assignment draw (`draw_posterior_row` in the selector) | `ASSIGN` | 0 (the joint multi-state draw is one draw) |

The `seed_rnd` calls at the fills' starts and the `ran_tabu` per thread go. A
run without `SIMPLE_SEED` takes a seed drawn once by the project owner and
written to the job description, so the parts share it; the reseed counter of
`simple_rnd.f90` (`nseed_rnd_calls`) stays for the classic path and the rest of
the code.

**What changes, once.** The draws' values change (a different generator, keyed),
so the probabilistic trajectory changes in phase 1. The same logical draw must
then return identical bits on any thread or part. Particle decisions must be
identical for identical prepared inputs and arithmetic. Reconstruction and
class-sum grouping, FFT plans and other workflow randomness remain separate
sources of differences in a complete run; section 11 defines their acceptance.
Keyed RNG alone does not make the existing identity suites pass across thread
or part counts.

**Validation of the change** (the truth suites, not identity): the simulated
workflows (`simulated_workflow` 1jxy and 6vxx, shared memory and distributed),
the multi-state truth cases (`multistate_workflows`, once its `labels` suite
exists), `mini_stream`, and the quality readouts on real data of the 2026-11-03
status (bgal solve2D over four seeds; trpv1 refine3D's FSC curves over three
iterations) against the current code. Accept when the distributions match
within the seeds' spread, as the master-drawn sample was accepted.

## 5. The pipeline, per route

**3D dense.** The batch loop of `refine3D_exec` takes over the fill: for each
batch, `fill_tab_range` on the batch's rows (the map's screen and polish per
particle, as now), then per particle the dense selection (the state by the
closest reference, the direction by the draw in its block; or the joint draw of
the multi-state proposal), the extras under `topk > 1` from the row and the
cone (`select_extras`), then `strategy3D_prob%srch_prob` as it is (polish,
commit, extras resolved), the residual, the record. The `eul_prob_tab` object
shrinks to a per-batch row store plus the assignment map of the part.

**3D neighbourhood.** The same with `eul_prob_tab_neigh%fill_tab_range` and the
pooled-candidate selection of `ref_assign_neigh`.

**2D dense and sparse.** The batch loop of `refine2D_exec` takes over the fill
(`fill_tab_range` or `fill_tab_prob_neigh_range` on the batch), then the class
draw (`ref_assign` or `ref_assign_sparse` per particle), then
`strategy2D_prob%srch_prob` as it is, then the batch's class sums. The particle
buffer's hand-over between the table and the search becomes moot: the batch's
prepared particles serve both. The stream's pool runs the same `refine2D`.

**Empty parts, fractional updates, `assign_all`, `update_missing`.** A part
with no sampled particles runs no batches and writes its empty documents, as
now. Sampling is unchanged (the master's draw, keyed). `update_missing` keeps
its refusal under the probabilistic modes.

**The readouts.** Global diagnostic row identities are chosen before
partitioning. Each part publishes sufficient statistics and counts; the master
merges them using the section 11 contract rather than summing per-part means.
The disk-set line moves to the master's sample draw, where its inputs are.

## 6. Shared memory and distributed, one pipeline

- `refine3D_exec` and `refine2D_exec` run the pipeline on the caller's
  `(params, build)`. The shared-memory strategies stop calling the `prob_align`
  commanders: no second toolbox, no reseed, no particle-field write and reread
  around the step, no table or assignment file. The layering policy's accepted
  exception ends.
- A distributed part runs the same `refine3D_exec` / `refine2D_exec` as the
  worker program it is today. The master draws the sample, writes the field
  (as `master_draws_sample3D` and `master_draws_sample2D` do), schedules the
  parts with the existing queue (one phase per iteration instead of two; no
  change to the scheduler), merges the documents, assembles, commits sigma2.
- The job description loses nothing and gains the run's seed.

**Retired** (after the gates of section 8): the private programs `prob_tab`,
`prob_tab_neigh`, `prob_tab2D`, `prob_align`, `prob_align_neigh`,
`prob_align2D` (`simple_private_exec_driver.f90` lines 322-333, the API and UI
entries) and their commanders; `prob_align_exec`, `prob_align_neigh_exec`,
`prob_align2D_exec` and `run_prob_align*`; the part-file writers and readers
(`begin_write`, `write_tab`, `read_tabs_to_glob`, the dense, neighbourhood and
candidate-stream formats), the assignment and insert-table files
(`write_assignment`, `read_assignment`, `write_insert_table`,
`read_insert_table`), the global table constructors (`new`, `new_neigh_global`)
and the master-side `ref_assign*` loops; `DIST_FBODY`, `ASSIGNMENT_FBODY`,
`INSERT_FBODY`; the `PROB_TAB*` and `PROB_ALIGN*` bench families (their
counters move into the matchers' reports). Roughly a third of the 9,000 lines
of section 2.1, and the layering policy's "worker frames in commanders"
deviation with them.

**Kept**: the fills, the selection rules, `draw_posterior_row`, the extras and
cones, the temperatures, the shift prior, the bench profile of the fill
(`prob_fill_profile`), the tests of the selection rules (section 8).

## 7. What it buys

- No master table and no part files: 713 MB and 4 x 199 MB per iteration on
  trpv1; 15 GB per state at 750,000 particles. A worker's alignment memory is
  the references, the thread scratch and one batch.
- One compatible preparation per particle per probabilistic alignment pass;
  one process start per part per iteration. Separate Cartesian polish and
  reconstruction still perform their required preparation. The historical
  duplicated-work estimate is 1.3 % of a 3D iteration and about 8 % of a 2D iteration at trpv1's
  size, more on runs with small shares per part (PfCRT: 1,700 particles per
  part).
- One code path for shared memory and distributed alignment; identical logical
  random draws across thread and part counts, with the numerical acceptance
  levels of section 11.
- The batch is the unit the batched and GPU kernels of the next release need
  (`speed_architecture_proposals.md`, Part 2) and the unit the coarray SPMD
  proposal (`coarray_spmd_execution.md`) runs on an image: the pipeline is their
  prerequisite, and neither is part of this work.

What it does not buy: the fill (264-333 s of a 3D iteration on trpv1) is
untouched. Wall time comes later, from the kernels, inside this batch.

## 8. Order of work and gates

Each phase ends with its gate on the Dell; nothing moves to the next phase
without it. The fast gate (`scripts/run_fast_gate.sh`) and `unit_pftc_align2D3D`
run at every phase.

Before phase 1, freeze the section 11 RNG schema/test vectors, immutable-input
and decision interfaces, restart metadata and numerical acceptance levels.
This is contract hardening of the target pipeline, not an intermediate worker
architecture. Keep evidence in this note. The proposed `labels` truth route
does not exist at this review's HEAD: implement its asserted gate or explicitly
use the existing asserted multi-state routes until it lands; never count a
missing route as a pass. Check declared artifacts before using a high-level
entry as a byte-identity target (`solve3D_addon` currently declares none).

**Phase 1: keyed draws.** `rnd_keyed` in `simple_rnd.f90` with its unit test
(determinism, independence of order and thread, the unit-interval mapping); the
five sites of section 4 keyed; the fills' `seed_rnd` and `ran_tabu` removed;
the seed, generator version and logical scope in the job/checkpoint metadata.
Gate: fixed-vector RNG tests, sparse permutation assertions, and identical
per-particle decisions on frozen prepared inputs at one and eight threads and
different part assignments; complete-run comparisons at different layouts use
the numerical/truth gates of section 11. Run the truth validation of section 4
and obtain Hans's acceptance of the trajectory change. Scope and line count
must be reassessed after the complete draw-site audit.

**Phase 2: the pipeline in shared memory.** The batch selection inside
`refine3D_exec` and `refine2D_exec` (section 5), the route adapters over the
existing fills and selection rules, the strategies calling the pipeline instead
of the commanders. The distributed path still runs the old two phases (the table
programs stay for it) and is checked against the shared-memory result. Gate:
shared memory byte-identical to phase 1 at one and eight threads
(`simulated_workflow` both systems, `solve3D_addon`, the 2D cases); the
distributed path unchanged. About 800 lines, mostly moved.

**Phase 3: the distributed parts on the pipeline.** The master's iteration
drops the probabilistic step and the second scheduling phase; the parts run the
pipeline; the readouts merged with counts. Gate: frozen-input alignment
decisions identical across shared/distributed layouts (1, 2, 3, 4 parts),
complete-run numerical/truth comparisons under section 11, `mini_stream` (the
stream's pool), `solve3D_addon` distributed, the nightly entries. About 300 lines, mostly
deleted.

**Phase 4: retirement and documents.** The programs, commanders, files and
objects of section 6; the policies and skills of section 10; the layering
policy's section 8; the bench families; the generated indexes. Gate: the full
nightly green on dell2; `scripts/check_test_registry.py`. About 3,000 lines
deleted.

**Phase 5 (its own proposal and ruling): batched kernels** inside the batch.
Results change (FFT plans, precision); not part of this work.

The measurements of Part 3 of `speed_architecture_proposals.md` (PfCRT, msp1,
the master's peak memory) are taken on phase 3's code and recorded in this
note.

## 9. Risks and open points

- **The trajectory change is isolated in phase 1**, before the structural
  move. Later phases preserve the phase-1 science and use the separate
  frozen-input identity and complete-run numerical/truth contracts of section 11.
- **The 2D sparse route's candidate stream** was a file format with a coverage
  check (each particle in exactly one part). In memory, the batch's candidate
  lists are the row; the coverage check becomes a check on the part's share.
- **`topk > 1` extras** need the row and the cone at draw time; both are in the
  batch. The extras' weights stay the table's; the polish does not reweight
  them (the review's rule).
- **The pose-posterior diagnostic** (`prob_post_diag`, the own- and other-half
  columns) is row-local and keeps its columns in the part's documents.
- **Multi-state.** The joint state draw of `multistate_posterior_state_draw.md`
  is one more selection rule in the same place; the refined-reference cap
  across states is a fill budget and orthogonal. The two developments should
  land in the order the maintainer sets; neither blocks the other.
- **The `solve2D` chunks and the stream's pool** run `refine2D` through its
  commander and must stay unchanged in behaviour; `mini_stream` and the chunk
  tests are in every gate from phase 2.
- **Memory per worker** must be measured. Existing workers already batch dense
  rows; eliminating the master's global table does not by itself reduce worker
  peak RSS. Fused preparation/polish workspaces must fit the section 11 budget.
- **The coarray backend** keeps working through phase 3 unchanged (it
  dispatches the same programs); the SPMD proposal builds on the pipeline
  afterwards.

## 10. Documents and skills to update (phase 4)

- `doc/policies/3D/refine3D_policy.md` section 8 (the pre-step, the part
  files, the assignment), `doc/policies/2D/solve2D_policy.md` (the table
  contract, the part files, the project-segment table),
  `doc/policies/importance_sampling_fractional_update_policy.md`,
  `doc/policies/layering_and_project_ownership_policy.md` (sections 2, 3 and
  8), `doc/policies/particle_read_path_policy.md` (the buffer's hand-over),
  `doc/policies/test_environment_policy.md` if the identity entries gain thread
  counts.
- `doc/algorithms/sampling_and_fractional_updates.md`,
  `doc/algorithms/refine2d_class_averaging.md`, `doc/algorithms/solve2d.md`.
- Skills: `simple-refine3d` (and its code map), `simple-solve2d`,
  `simple-solve3d-importance-sampling` (and its coupling map),
  `simple-frac-update-trailing`, `simple-main-commanders`,
  `simple-main-strategies`, `simple-main-apis`.
- `speed_architecture_proposals.md`: the status of this work; Part 2 reduced to
  the batched kernels, the persistent workers and the SPMD route.

## 11. Source review, 2026-11-04: contracts to settle before implementation

Reviewed against `6b924ccf5`, including the current fill, matcher, sampling and
reduction code and the earlier architecture findings. The batch pipeline is
the right final direction. No intermediate marker-file worker fold or retained
global assignment service is needed. The corrections below are design review
requirements, not evidence that the new pipeline has been implemented or tested.

### 11.1 P1: make the iteration snapshot an enforced input boundary

The current `fill_tab_range` is not merely a numerical row evaluator. At each
call, dense 3D derives `projs_ns` through `calc_num2sample` on the live project
field, and the shifted fills derive the prior through `shift_prior_sigma2`.
The former reads sampled distances/state labels; the latter reads sampled
`shincarg` (`simple_eul_prob_tab_utils.f90`, `calc_athres` and
`shift_prior_sigma2`). Finalization updates those fields. Calling the old fill
unchanged after committing an earlier batch would therefore change later
batches' budgets and priors. It would also repeat full-field reductions per
batch. Changing batch size could then change the scientific calculation.

Compute these values once in `begin_alignment`, from the complete incoming
field, before any decisions are published. Pass them with the immutable
reference/noise/class-population model into every fill and polish consumer.
Read incoming per-particle metadata from the input view, and write decisions
to a separate output view. Batch routines must not derive iteration statistics
from the evolving output project. Do not copy an entire `parameters` or builder
to achieve this; retain only required scalars/arrays and explicitly owned or
borrowed domain resources.

Add a regression that changes earlier batches' output distances, states and
shift increments while leaving the iteration input fixed: later candidate rows
must remain identical. Repeat with different batch boundaries and part layouts.
Test row selectors with hand-constructed rows and explicit variates against an
independent cumulative-posterior oracle, including ties, invalid/empty rows,
state blocks, sparse eligibility and frozen extras; full workflow truth tests
alone do not establish this contract.

### 11.2 P1: finish the RNG namespace, permutation and restart contract

The five-family table in section 4 leaves several identities unspecified.
Iteration numbers and reference indices alone do not distinguish repeated
refine commands, solve stages, assignment-only/final passes or diagnostic
evaluations of the same reference. A screen ordinal also changes meaning when
the candidate ordering changes. Define a logical key schema before coding:

- Seed and algorithm/schema version; deterministic command/stage/pass identity;
  iteration identity within that scope; canonical global particle identity;
  draw family; canonical state/reference identity where applicable; and an
  explicit subdraw/evaluation ordinal. Specify intended reuse of a variate
  versus an independent repeated evaluation. Diagnostics get their own family.
- Never key on a thread, part, image, batch slot, arrival order, absolute path
  or process identifier. Global particle IDs come from the authoritative input
  manifest and must not be replaced by partition-local indices.
- Pin the exact generator, key packing, integer arithmetic and float mapping,
  including endpoint behavior, with fixed expected bit vectors in supported
  Debug/Release and Mac/Linux builds. Do not rely on unspecified signed integer
  overflow. A generic "64-bit mixer or equivalent" is not an executable contract.
- Sparse candidate generation must remain sampling without replacement.
  Removing `ran_tabu` needs a defined keyed permutation/partial permutation,
  including bounded-integer mapping and rejection subdraws; independent class
  index hashes do not preserve that distribution. Preserve the previous-class
  ordering rule and existing candidate-count/stopping rules.
- Persist the seed, RNG version, logical command/stage/pass and next iteration,
  plus input/sample/model identity, in restart metadata as well as the job
  description. A retry of a logical iteration reuses its keys; a distinct pass
  gets a distinct scope. Reject inconsistent restart metadata. Transport the
  seed losslessly as an integer/encoded bit string, never through a real value.

Audit stochastic inputs outside the five named sites too, including reference
initialization and classic stages following probabilistic stages. Removing
command/toolbox initialization can change the intrinsic RNG's reseed counter
even when the new probabilistic draws are keyed. Phase 2's complete-workflow
identity needs an explicit caller-RNG contract; it is not automatic.

### 11.3 P1: distinguish draw identity from complete-run numerical identity

The present 3D master adds partial reconstructions in part order
(`simple_volassemble_strategy`, `reduce_partials`); 2D adds part contributions
(`cavg_sums%accumulate`). Changing `nparts` changes the grouping of floating
point additions. Thread count also changes batch sizes, and FFT/preparation
paths can change rounding. Different reconstructed references feed subsequent
iterations, where an arbitrarily small score difference can change a draw.
These dependencies remain even with perfect keyed RNG.

The acceptance levels are:

1. Identical logical keys produce identical random bits regardless of execution
   order, threads, batch size, partitioning and restart.
2. Identical frozen prepared inputs and numerical operations produce identical
   candidate rows, decisions, fallbacks, shift metadata and frozen extras under
   different schedules. Compare by canonical particle/reference IDs.
3. Structural phases preserve outputs against phase 1 at each fixed execution
   configuration. Compare shared and distributed numerical handoffs explicitly.
4. Complete runs at different numerical layouts have measured numerical
   tolerances on fixed-input accumulators/maps and unchanged multi-seed truth
   gates; stochastic trajectories need not remain byte-identical. Freeze the
   tolerances and truth floors before evaluating the change.

This review recommends these levels for this refactor. If byte identity of
complete runs across all part counts is itself required, add a separately
costed reproducible-reduction/FFT/preparation and initialization contract before
phase 1. It cannot be promised while leaving those operations unchanged.

### 11.4 P1: freeze scientific interfaces independently of execution adapters

Use one strategy-owned alignment context and one bounded batch workspace, with
explicit operations corresponding to:

```text
begin_alignment(input model, incoming metadata, sampling/RNG scope, consumers)
prepare_batch(global particle IDs, canonical observations, workspace)
evaluate_rows(model, incoming metadata, RNG scope, workspace)
select_rows(model, RNG scope, workspace)
finalize_batch(model, workspace, decisions, noise, canonical records)
publish_results(decisions, noise, diagnostics)
end_alignment() -> release alignment resources before 3D reconstruction
```

These are ownership/API boundaries, not a request for seven new modules or a
second implementation of the existing algorithms. The outer entry still runs
on `(params, build)`; inner operations receive the specific inputs they need.
Geometry-compatible plans, map setup, optimizers and scratch belong to the
alignment context and are reused where safe, rather than rebuilt on every
batch call. Mutable numerical scratch remains per thread. Model generation,
geometry/band, noise, masks, halves and preprocessing invalidate caches
explicitly. Incompatible input/model generations are rejected at the handoff.

Define the decision schema: canonical particle/state/reference IDs, pose and
score, absolute and incremental shifts with units/frame, shift width/state,
fallback status, frozen insertion poses/weights, and validity flags for optional
consumers. Define input/output ownership and exact-once coverage. Only the
publication layer writes changed project rows or partition artifacts. Domain
evaluation and selection know neither commanders, file names, queues nor images.
GPU batching can replace evaluation through this boundary; SPMD can replace
publication and iteration transport without changing the scientific pipeline.

Keep canonical observations distinct from normalized/masked alignment PFTs.
At the 3D barrier no decision or diagnostic may keep a pointer to a released
batch/PFTC/reference allocation. Keep the current one-state/one-half
reconstruction and budgeted record-buffer fallback. The proposed context owns
alignment resources only during alignment, not through reconstruction.

### 11.5 P2: budget the full working set and preserve diagnostic meaning

Existing worker tables already bound their dense rows by batch size, so the
guaranteed memory saving is removal of the master's global table and its
transport buffers. Fusing contexts may increase worker peak memory if old and
new polish workspaces coexist. Budget all-state references, per-thread map and
optimizer scratch, canonical buffers, prepared batches, complete/refined rows,
optional diagnostics and result metadata using checked 64-bit byte arithmetic.
Choose batch capacity within the remaining budget, not solely `nthr * 50`.
Measure worker/root and aggregate-node peaks, including the 3D barrier.

The current `n_eff` diagnostic selects a bounded subsample by global row
position. Per-part subsampling followed by summing means is not equivalent and
changes with `nparts`. Choose diagnostic rows from the global sampled manifest
before partitioning and publish sufficient sums/counts (and defined merge order
where exact diagnostic reproducibility is required). Integer row-mode/coverage
counts can be summed directly. Preserve the future multi-state readout's
population, mass and normalization meanings through the same result boundary.

Before retiring the old path, exercise all four routes, topk, shift prior and
map/window/disk cases, polish/noise on and off, fractional sampling, empty parts,
raw/cropped and over-budget observations, gridding/PCG, restart/retry, chunk and
stream callers. Against the keyed phase-1 baseline, measure unchanged
candidate/map work, one compatible preparation
per particle in the fused alignment phase, no shared-memory table transport,
bounded memory and actual end-to-end wall time. The reported 1-8 percent is a
case-specific historical estimate, not a new speed claim.

The durable decision is the pipeline with immutable inputs, explicit randomness,
bounded workspaces and typed publication. Numerical kernels, selection policies
and transports can evolve behind those boundaries without reopening ownership
or reinstating a global probabilistic pre-step. No architecture can guarantee
that numerical methods or hardware will never need further changes.
