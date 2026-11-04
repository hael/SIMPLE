# Probabilistic alignment as a batch pipeline: per-particle draws, keyed random streams, no global table

Planned 2026-11-04 on master `a723f4919`, ruled done by Hans the same day ("we
have reached the point of no return; this refactoring is obvious"), reviewed by
Hans against `6b924ccf5` the same day (the review's requirements are folded into
sections 3, 4, 5 and 8), and planned as the autonomous run `prob_pipeline` on the
Dell (sections 10-14). **Nothing in sections 1-9 is implemented.** The frame is
the architecture review of 2026-11-04 in `speed_architecture_proposals.md`; this
note is its execution. It supersedes that note's Part 2 rows on the worker fold,
the streaming assignment, the read-ahead thread and the 2D table pre-step, and it
closes the worker fold implemented and dropped on 2026-11-04.

## 1. What is proposed, in one paragraph

Since 2026-10-24 every probabilistic assignment is one draw from one particle's
own row; no population balance, no other particle's row, nothing enters a draw
that is not fixed before the iteration's particle loop. The global table, the
part files, the master's assignment step and the second worker pass exist only
to serialise those draws through one random stream. Give every draw its own
stream, keyed on the run, the iteration, the particle and the draw, and the
table becomes unnecessary: one pipeline prepares a batch of particles, scores
its rows, draws each particle's pose and extras, polishes the drawn pose,
computes the noise residual and commits the decision, then releases the batch.
The same pipeline runs in shared memory on the caller's `(params, build)` and
in every distributed part; the master draws the sample, waits, merges the
decisions and assembles, as it does today. The stochastic trajectory changes
once, when the draws move to keyed streams; the draws are then independent of
thread and partition placement. Byte identity of a complete run across thread
or part counts is not promised: the reductions, the FFT plans and the batch
sizes keep their rounding (section 8 states what is required instead).

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
model):

- the sample (`sample4update_cnt`, `_rnd`, `_fillin`, drawn by the master or
  the shared-memory strategy, written to the project, reproduced by
  `sample4update_reprod`);
- the references (the reprojection model, or the class averages), the halves,
  the band, the matching low-pass limit, the shift window and map geometry;
- the committed sigma2 spectra;
- the shift prior's width (`shift_prior_sigma2`, from the previous iteration's
  sampled `shincarg`); the particle temperatures' formula;
- the per-state FSC cones of the weighted insertion (`insert_cone_from_fsc`);
- the 2D class populations of the previous assignment (`prep_strategy2D_glob`);
- `prob_beta`, the extremal fraction, the refined-reference budget `projs_ns`
  (`calc_num2sample`).

Two of these are computed *inside* the fills today, from the live particle
field: dense 3D derives `projs_ns` through `calc_num2sample` at every
`fill_tab_range` call (`simple_eul_prob_tab.f90` line 306, reading the sampled
particles' distances and labels through `calc_athres`), and every shifted fill
derives the prior's width through `shift_prior_sigma2` (line 340; `_2D` 324 and
554; `_neigh` 218 and 470, reading the sampled `shincarg`). The finalisation
writes those same fields. A fill called after an earlier batch's decisions were
committed would therefore see a different budget and prior, and the result would
depend on the batch size. Section 3 moves both to the iteration's model.

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
`% at the row mode`) and the disk-set line (`print_disk_set`) are diagnostics
over particles, not inputs to any draw; `n_eff` reads a bounded subsample chosen
by global row position, which section 5 keeps meaningful across parts.

### 2.3 The random draws

Every draw on the probabilistic path, with where it runs today:

| draw | where | stream today |
| --- | --- | --- |
| the sample | master or shared-memory strategy, before the loop | the process's intrinsic stream |
| the in-plane draw per scored reference | the fill (`draw_posterior_crvec`, `simple_polarft_corr.f90` lines 672-740) | the fill thread's stream (gfortran keeps one per thread); each fill calls `seed_rnd` at its start (`simple_eul_prob_tab.f90` 302, `_2D` 282 and 528, `_neigh` 197 and 451) |
| the disk draw of the polar shift map | the fill (`shiftmap_draw`, `simple_polarft_shiftmap.f90` line 607) | the fill thread's |
| the candidate picks of the sparse fills | the fill (`ran_tabu` per thread, `_2D` 539, `_neigh` 205): sampling without replacement | the fill thread's |
| the assignment draw | the master, in particle order | the master's |

The fills run under a dynamic schedule, so which draws a particle meets depends
on the thread it lands on: byte identity holds at one thread only (the policies
say so). The classic strategies' samplers (`sample_bounded_dists`,
`sample_power_dists`, `shcloc`, the 2D search order) are not on this path, and
neither is the reference initialisation of a workflow's first stage; the run
audits every other `ran3`, `irnd_uni`, `ran_tabu` and `seed_rnd` on the
probabilistic path before phase 1 (section 13).

### 2.4 The reductions

The 3D master adds the parts' partial reconstructions in part order
(`simple_volassemble_strategy.f90`, `reduce_partials`, `do part = 1,
params%nparts`); 2D adds the parts' class-sum contributions
(`cavg_sums%accumulate`). The grouping of those floating-point sums changes
with `nparts`; the batch size (`nthr x BATCHTHRSZ`) changes with the thread
count; the FFT plans and the preparation paths carry their own rounding. A map
that differs in the last bit feeds the next iteration, where an arbitrarily
small score difference can move a draw. None of that changes with keyed
streams, which is why section 8 separates draw identity from complete-run
identity.

### 2.5 Measured on the current code (Dell, 2026-11-03/04)

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
structure costs. The percentages are this case's; they are not a speed claim
for the pipeline.

## 3. The contract

### 3.1 Lifetimes

| lifetime | holds |
| --- | --- |
| the run (`refine3D`, `refine2D` command; a distributed part's program) | partition identity, the project owner's project, control metadata; nothing of the alignment |
| the iteration's model, immutable from before the particle loop to after it | references (polar model or class averages, both halves), band and grid, shift-map geometry, committed sigma2 spectra, the sample, the incoming orientations, shifts and widths (the input view), the temperature formula, `prob_beta`, the shift prior's width, the refined-reference budget, the insertion cones, the 2D class populations |
| the alignment phase's workspace | compatible FFT plans, the map's setup, the optimisers, per-thread numerical scratch; reused between batches, released before the 3D reconstruction |
| one batch (capacity from the memory budget of section 3.5, at most `nthr x BATCHTHRSZ` to begin with) | the canonical observations, the prepared polar particles and CTF terms, the particle temperatures (computed from the prepared particles), the complete candidate rows, the drawn poses and extras, the polished poses, the residuals; released when its decisions are committed |
| across the 3D alignment/reconstruction barrier | the final poses, shift states, the frozen insertion poses and weights, the sigma2 contributions, the record buffer of cropped observations (as now); no alignment calculator, reference memo or row |

### 3.2 The iteration's model is an enforced input boundary

Everything in the model row above is computed once, in `begin_alignment`, from
the complete incoming particle field before any decision is published, and
passed to every fill and polish consumer. In particular `projs_ns`
(`calc_num2sample`) and the shift prior's width (`shift_prior_sigma2`) leave
the fills (section 2.2). The fills and the selector read per-particle metadata
from the input view and write decisions to an output view; no batch routine
derives an iteration statistic from the evolving output. No whole `parameters`
or builder is copied for this: the model holds the scalars and arrays it needs
and borrows the domain objects (the calculator, the references) explicitly.

The regression that guards it (section 12): with the iteration's input fixed,
change earlier batches' output distances, states and shift increments; later
batches' rows must be identical, at several batch boundaries and part layouts.

### 3.3 The operations

One strategy-owned alignment context and one bounded batch workspace, with
these operations (ownership boundaries, not seven modules; the existing fills,
selection rules, polish and residual routines are reused behind them, not
rewritten):

```text
begin_alignment(input model, incoming metadata, sampling/RNG scope, consumers)
prepare_batch(global particle IDs, canonical observations, workspace)
evaluate_rows(model, incoming metadata, RNG scope, workspace)
select_rows(model, RNG scope, workspace)
finalize_batch(model, workspace, decisions, noise, canonical records)
publish_results(decisions, noise, diagnostics)
end_alignment() -> release alignment resources before 3D reconstruction
```

The outer entry runs on `(params, build)`; the inner operations receive the
specific inputs they need. Plans, map setup, optimisers and scratch belong to
the context and are reused where safe; mutable numerical scratch stays per
thread. A new model generation, band, noise, mask, half or preprocessing
invalidates the context's caches explicitly, and an input of another generation
is rejected at the hand-off. Domain evaluation and selection know neither
commanders, file names, queues nor images; only `publish_results` writes
changed project rows or partition artefacts. GPU batching later replaces
`evaluate_rows` through this boundary; SPMD replaces publication and the
iteration's transport without touching the science.

The batch path, with the existing routine behind each step:

```text
batch of sampled particles (the part's share, in order)
  -> read, or take from the record buffer      (ptcl_read_session)
  -> prepare the alignment batch               (build_batch_particles3D / 2D)
  -> fill the complete rows                    (the route's fill on the batch)
  -> draw each particle's pose and extras      (the route's selection rule over draw_posterior_row)
  -> polish the drawn pose                     (refine_assignment_continuously / refine_selected_continuously)
  -> sigma2 residual of the committed pose     (as now)
  -> commit: orientation, shift state, extras  (as now)
  -> 2D: the batch's class sums; 3D: the record for the reconstruction
  -> release the rows
```

The row never leaves the process: the sparse routes keep their candidate lists
as the row (an adapter over the candidate stream's records), the dense routes
the batch's `batchsz x nrefs` entries (at 2,500 references over all active
states and 300 particles, 6 MB for distances and in-plane indices; refined
entries, diagnostics and scratch are additional).

**One selector.** `draw_posterior_row` stays the primitive. The route adapters
give it the row (the dense block of the drawn state, the neighbourhood's pooled
candidates, the 2D class row or candidate list) and take back the draw, as
`ref_assign`, `ref_assign_neigh`, `ref_assign_sparse` and the 2D `ref_assign`
do today; those four become the adapters' selection step, called per particle
inside the batch instead of over the global table. The multi-state joint draw
lands in the same place.

**The decision schema**: canonical particle, state and reference identities;
pose and score; absolute and incremental shifts with their unit and frame;
shift width and state; fallback status; the frozen insertion poses and weights;
validity flags for the optional consumers. Exact-once coverage of the part's
share is checked at publication.

**Consumers decide the preparation.** Scoring and polish need the polar
particles; the sigma2 residual needs the committed pose against the reference;
the 3D reconstruction needs the cropped record; the 2D class sums need the
cropped (or raw) particle. Nothing is prepared twice for the alignment, and
nothing is prepared for a consumer that is off (`volrec=no`, `objfun=cc`). The
canonical observations stay distinct from the normalised, masked alignment
transforms: reusing one as the other needs a proven preparation contract (CTF,
band, half, noise, centring, frame), and none is assumed here.

### 3.4 What stays as it is

The master's sample draw and project ownership
(`layering_and_project_ownership_policy.md`); `merge_algndocs`; volume assembly
and the sigma2 commit; the 3D partial reconstruction after the loop from the
record buffer, one state and one half at a time
(`separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`
stays a hard boundary; at the barrier no decision or diagnostic keeps a pointer
into a released batch, calculator or reference allocation); the 2D class-sum
order (batch, class, half); the `pose_cont` Cartesian polish pass, which stays
a pass of its own and keeps its preparation; the classic modes (`shc`,
`snhc_smpl`, ...), untouched; the stream's pool and the `solve2D` chunks, which
run `refine2D` through its commander and are not edited.

### 3.5 The memory budget

The guaranteed saving is the master's global table and its transport. A worker
already bounds its dense rows by the batch; fusing the table's and the matcher's
workspaces can raise its peak if both coexist. The context budgets, in checked
64-bit byte arithmetic, the references of all active states, the per-thread
map and optimiser scratch, the canonical buffers, the prepared batch, the
complete and refined rows, the optional diagnostics and the result metadata,
and chooses the batch capacity inside the remaining budget
(`ptcl_buffer_budget` is the model: half the task memory the queue grants).
Worker, master and whole-node peaks are measured, including at the 3D barrier
(section 12).

## 4. The random-stream contract

### 4.1 The generator

A keyed, stateless draw `rnd_keyed(key) -> real(dp) in [0,1)`: the key's words
are mixed by a multiply-xor-shift hash on 32-bit lanes with 64-bit products
(a 32-bit lane times a 32-bit constant fits a 64-bit integer; no step relies on
integer overflow, which Fortran leaves undefined), finalised to 64 bits, and the
top 53 bits divided by 2^53 (so 0 is possible and 1 is not). The exact
constants, the lane order, the finaliser and the float mapping are pinned by
fixed test vectors in `simple_rnd_tester.f90`, checked in the Debug and Release
builds of this run and recorded in this note for the Mac build. The generator
carries a version in its key; changing any of this changes results and bumps
the version. One hash per draw; the kernels dwarf it.

### 4.2 The key schema

Every draw's key is the tuple

| field | what |
| --- | --- |
| seed | the run's 64-bit seed: `SIMPLE_SEED` when set, else drawn once by the project owner at the run's start and exported to the parts as `SIMPLE_SEED` through the environment the queue scripts inherit (`put_simple_seed`, as today), never through a real value or a new key |
| version | the generator's schema version |
| scope | the logical pass: the program (`refine3D`, `refine2D`), the pass kind (discrete alignment; the `pose_cont` polish; the final full-population assignment; a diagnostic evaluation), which together with the iteration distinguish repeated commands and solve stages (a solve stage's `refine3D` continues `which_iter`, so stage and iteration are one number) |
| iteration | `which_iter` within the scope |
| particle | the particle's row in the run's particle field (global, the same on the master and every part, as `fromp`..`top` are); never a batch slot, a part-local index, a thread or a process |
| family | `SAMPLE`, `INPL`, `DISK`, `CAND`, `ASSIGN`, `DIAG` |
| reference | the canonical reference (state, direction or class) where the family has one |
| ordinal | the sub-draw's ordinal within the family for the particle (the screen's ordinal, the pick's ordinal, the sampler's draw count) |

Reuse versus independence is explicit: an evaluation that is repeated on the
same reference in the same pass reuses the variate (same key); a distinct pass
or a diagnostic gets a distinct scope or family. The sample's draws belong to
the project owner (particle 0, the sampler's ordinal). Nothing is keyed on a
thread, a part, an image, an arrival order, a path or a process.

### 4.3 The sites

| draw | family | keyed on |
| --- | --- | --- |
| the sample (`sample4update_cnt`, `_rnd`, `_fillin`) | `SAMPLE` | particle 0, the sampler's draw ordinal |
| the in-plane draw per scored reference (`draw_posterior_crvec`) | `INPL` | the reference |
| the disk draw (`shiftmap_draw`) | `DISK` | the screen's ordinal for the particle |
| the sparse fills' candidate picks | `CAND` | a keyed partial Fisher-Yates permutation of the candidate set (sampling without replacement, as `ran_tabu` gives today), one keyed uniform per swap (the swap's ordinal), the previous class or reference first as today, the candidate count and stopping rules unchanged |
| the assignment draw (`draw_posterior_row` in the selector) | `ASSIGN` | ordinal 0 (the joint multi-state draw is one draw) |
| the pose-posterior diagnostic and any other diagnostic draw | `DIAG` | its own ordinals |

The `seed_rnd` calls at the fills' starts and the per-thread `ran_tabu` go. The
reseed counter of `simple_rnd.f90` (`nseed_rnd_calls`) stays for the classic
path and the rest of the code, and the shared-memory strategies' removal of the
commander detour (phase 2) changes that counter's count: the classic stages that
follow a probabilistic stage in a workflow therefore keep their own
reproducibility only through the existing `SIMPLE_SEED` mechanism, which the
run checks (section 12, level 3) rather than assumes.

### 4.4 Restart

The seed, the generator version, the scope and the next iteration are written
with the run's restart metadata (the job description and the project's
`projinfo`, where `refine3D continue=yes` and the stream's restart read their
state) as well as the environment. A retry of a logical iteration reuses its
keys; a distinct pass gets a distinct scope; inconsistent restart metadata is
rejected.

### 4.5 What changes, once

The draws' values change (a different generator, keyed), so every probabilistic
trajectory changes once, in phase 1, before any structural move. The same
logical draw then returns the same bits on any thread, in any batch, in any
part; identical prepared inputs and arithmetic give identical rows and
decisions wherever they run. A complete run still differs across layouts by the
rounding of section 2.4; section 8 says what is accepted.

## 5. The pipeline, per route

**3D dense.** The batch loop of `refine3D_exec` takes over the fill: for each
batch, `fill_tab_range` on the batch's rows (the map's screen and polish per
particle, as now, with `projs_ns` and the prior's width from the model), then
per particle the dense selection (the state by the closest reference, the
direction by the draw in its block; or the joint draw of the multi-state
proposal), the extras under `topk > 1` from the row and the cone
(`select_extras`), then `strategy3D_prob%srch_prob` as it is (polish, commit,
extras resolved), the residual, the record. The `eul_prob_tab` object shrinks
to a per-batch row store plus the assignment map of the part.

**3D neighbourhood.** The same with `eul_prob_tab_neigh%fill_tab_range` and the
pooled-candidate selection of `ref_assign_neigh`.

**2D dense and sparse.** The batch loop of `refine2D_exec` takes over the fill
(`fill_tab_range` or `fill_tab_prob_neigh_range` on the batch), then the class
draw (`ref_assign` or `ref_assign_sparse` per particle), then
`strategy2D_prob%srch_prob` as it is, then the batch's class sums in the order
of today. The particle buffer's hand-over between the table and the search
becomes moot: the batch's prepared particles serve both. The stream's pool runs
the same `refine2D`.

**Empty parts, fractional updates, `assign_all`, `update_missing`.** A part
with no sampled particles runs no batches and writes its empty documents, as
now. Sampling is unchanged (the master's draw, keyed). `update_missing` keeps
its refusal under the probabilistic modes.

**The readouts.** The rows of the `n_eff` readout are chosen from the global
sampled manifest before partitioning (the same rows whatever `nparts`); each
part publishes sums and counts over its share of them, and the master merges
them in part order (so the readout is exactly reproducible), never means of
means. Integer counts (row mode, coverage, disk set) are summed. The disk-set
line moves to the master's sample draw, where its inputs are. The multi-state
readouts (population, mass, normalisation) go through the same result boundary.

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
- The job description is unchanged; the seed travels in the environment.

**Retired** (phase 4): the private programs `prob_tab`, `prob_tab_neigh`,
`prob_tab2D`, `prob_align`, `prob_align_neigh`, `prob_align2D`
(`simple_private_exec_driver.f90` lines 322-333, the API
module) and
their commanders; `prob_align_exec`, `prob_align_neigh_exec`,
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
(`prob_fill_profile`), the tests of the selection rules.

## 7. What it buys

- No master table and no part files: 713 MB and 4 x 199 MB per iteration on
  trpv1; 15 GB per state at 750,000 particles. A worker's alignment memory is
  the references, the thread scratch and one batch, within the budget of
  section 3.5 (measured, not assumed).
- One compatible preparation per particle per probabilistic alignment pass and
  one process start per part per iteration. The Cartesian polish pass and the
  reconstruction keep theirs. On trpv1 that is 1.3 % of a 3D iteration and
  about 8 % of a 2D one; more where shares per part are small (PfCRT: 1,700
  particles per part). Not a speed claim: the run measures it (section 12).
- One code path for shared memory and distributed alignment; identical logical
  draws across thread and part counts, with the acceptance levels of section 8.
- The batch is the unit the batched and GPU kernels of the next release need
  (`speed_architecture_proposals.md`, Part 2) and the unit the coarray SPMD
  proposal (`coarray_spmd_execution.md`) runs on an image: the pipeline is their
  prerequisite, and neither is part of this work.

What it does not buy: the fill (264-333 s of a 3D iteration on trpv1) is
untouched. Wall time comes later, from the kernels, inside this batch.

## 8. Acceptance levels

The run holds every phase to these; none is relaxed.

1. **Draw identity.** Identical logical keys give identical random bits
   whatever the execution order, thread count, batch size, partitioning or
   restart (fixed test vectors; the permutation test of the sparse picks).
2. **Frozen-input identity.** Identical prepared inputs and arithmetic give
   identical candidate rows, decisions, fallbacks, shift metadata and frozen
   extras under different schedules (one and eight threads, two batch sizes,
   one and three parts), compared by canonical particle and reference
   identities, on the four routes.
3. **Structural identity at a fixed configuration.** Each structural phase
   (2, 3, 4) reproduces the previous phase's outputs byte for byte at each
   fixed execution configuration (the declared outputs of the registry's
   entries, `test=identity`, at one thread and at eight; shared memory and
   distributed each against itself), including the classic stages that follow
   a probabilistic stage in a workflow.
4. **Complete-run numerical tolerance across layouts.** Runs of one build at
   different layouts (shared memory; 2, 3 and 4 parts; one and eight threads)
   agree on fixed-input accumulators (a reconstruction and a class-average set
   from the same committed poses) within a tolerance frozen in phase 0 from the
   *base code's* own layout-to-layout differences (twice the base's largest
   observed difference, recorded in Progress), and the multi-seed truth gates
   (the suites' floors) stay unchanged and met. Stochastic trajectories need
   not match.

Byte identity of complete runs across part counts is not required by this
note. If the maintainer requires it, a reproducible-reduction, FFT and
preparation contract is costed separately before phase 1; it cannot be had
while those operations are left as they are.

## 9. Risks and open points

- The trajectory change is isolated in phase 1, before the structural move.
- The 2D sparse route's candidate stream was a file format with a coverage
  check (each particle in exactly one part); in memory the batch's candidate
  lists are the row, and the coverage check becomes exact-once coverage of the
  part's share at publication.
- `topk > 1` extras need the row and the cone at draw time; both are in the
  batch. The extras' weights stay the table's; the polish does not reweight
  them.
- The pose-posterior diagnostic (`prob_post_diag`, the own- and other-half
  columns) is row-local, draws under the `DIAG` family and keeps its columns
  in the part's documents.
- Multi-state: the joint state draw of `multistate_posterior_state_draw.md` is
  one more selection rule in the same place; the refined-reference cap across
  states is a fill budget and orthogonal. The two developments touch the same
  files (ruling R1).
- The coarray backend keeps working through phase 3 unchanged (it dispatches
  the same programs); the SPMD proposal builds on the pipeline afterwards.
- `solve3D_addon` declares no outputs in the registry, so it is a pass/fail
  gate, not an identity target; the identity targets are the entries that
  declare outputs (`simulated_workflow_*`, `mini_stream_*`, `stream_*`).
- The run can test the Linux gfortran builds only; the generator's Mac vectors
  are recorded for Hans to check on the Mac before phase 1's diff is accepted.

## 10. Rulings fixed before the launch

Fixed by Hans on 2026-11-04 unless marked **(needed)**, which the launch waits
on:

1. The refactoring is done (the pipeline of sections 3-6); no marker-file
   worker fold, no retained global assignment service.
2. No parameter is added, internal or public; the seed travels as
   `SIMPLE_SEED` in the environment. If the work seems to need a key, the run
   stops and writes the question up.
3. The acceptance levels of section 8; complete-run byte identity across part
   counts is not required.
4. No floor, tolerance or timeout of an existing test is relaxed; the layout
   tolerance of level 4 is frozen in phase 0 from the base code's own
   differences and never moved.
5. The stream's code (`src/main/stream`, `src/main/commanders/stream`) and the
   chunks' code are not edited; their behaviour is gated (`mini_stream`,
   `stream_*`, the chunk tests).
6. The classic modes, the `pose_cont` pass, the reconstruction after the loop
   and the 2D class-sum order are unchanged.
7. **(needed) R1, the order against `ms_posterior2`.** Both developments edit
   `simple_eul_prob_tab.f90`, `_utils`, `_neigh`, `simple_prob_align_strategy.f90`,
   `simple_refine3D_strategy.f90` and `simple_commanders_prob.f90`. Proposed:
   this run goes first; `ms_posterior2`'s Phase 0 (the `labels` route, P2S:
   test and data files only) may run on dell2 meanwhile, and its Phase 1
   onward starts on the tree this run leaves, with the joint draw as a
   selection rule of section 3.3.
8. **(needed) R2, the generator and key schema** of section 4 as written
   (lanes, families, the restart fields), with the exact constants and test
   vectors produced in phase 0 and reviewed in the diff.
9. **(needed) R3, the `labels` route**: it does not exist at this HEAD; the
   truth gates of this run are the existing asserted routes of
   `multistate_workflows` (`states`, `cavgs`, `shc`, `consensus`) and the
   simulated workflows. A missing route is never counted as a pass.

## 11. The machine and the data

**The Dell** (Oracle Linux 8.10, user elmlundho): Core Ultra 9 285, 24 cores
(8 performance cores, CPUs 0-7, and 16 efficiency cores; no hyperthreading),
62 GB, `/data2` 1.2 TB free, `/home` 16 GB free. The run root is
`/data2/agent_runs/prob_pipeline`, reached as `~/agent_runs/prob_pipeline`
through a symbolic link; nothing of the run goes under `/home`. GCC 15.2 from
gcc-toolset-15, CMake 3.26, `python3.12` for CMake (`python3` is 3.6, no
numpy); `OPENBLAS_CORETYPE=Haswell` is exported; `TMPDIR=$HOME/tmp`. The
maintainer uses the machine interactively during the run: timed measurements
record the load, and a measurement under load is repeated.

**Builds on the Dell**: `./compile_clean.sh` (Release, tests on, the fast gate,
install into `build/`) in about 4 minutes; increments by `make -C build -j20`,
`scripts/run_fast_gate.sh build` (15 suites, 1-2 minutes), `make -C build
install`. The base HEAD is exported with `git archive` into
`scratch/keep/src_base` and built in `scratch/keep/build_base` for the
identity and tolerance baselines; HEAD is never moved.

**Timings on the Dell** (2026-11-03, for planning): `simulated_workflow_1jxy`
about 5 minutes at one thread to its solve3D stage, `simulated_workflow_6vxx`
18 minutes at one thread, `mini_stream_1jxy` 1-2 minutes at 8 threads,
`multistate_workflows` about 55 minutes at 8 threads (one route at a time),
distributed `solve2D` of the 1jxy project 30 s, distributed `solve3D` of the
6vxx project 7-10 minutes at one thread per part, `refine3D` of trpv1 three
iterations at 4 x 6 threads 17 minutes, `solve2D` of the bgal selection 35 s,
`test=nightly` about 2 hours.

**Data** (read-only; work on copies under `scratch`):

- the simulated fixtures the tests build themselves (1jxy, 6vxx, the
  multi-state beta-galactosidase fixture, the stream's simulated movies);
- `/data2/trpv1_c2`: 35,645 particles, c2, `2_solve2D/trpv1_c2.simple` with
  its class averages, `5_solve3D/trpv1_c2.simple` with its map
  (`recvol_state01.mrc`); the stack `/data/trpv1/flipnew.mrc` (9.3 GB, box
  256 at 1.2156 A); `mskdiam=175`;
- `/data2/bgal_addon/3_selection/bgal_addon.simple`: 5,513 particles, 2,756
  active, `mskdiam=180`;
- the scaled references `/data2/agent_runs/speed_items_2026-11-03/runs/nofold/refs120/cavgs120.mrc`
  (trpv1's class averages at box 120) and the run scripts of 2026-11-03 under
  `/data2/agent_runs/speed_items_2026-11-03/tools` (adapt copies).

PfCRT and msp1 are not on this machine; Part 3's measurements on them stay with
the maintainer.

## 12. How the work is tested and judged

**The ladder.** A change goes to a slower rung only when it passes the faster
ones: the unit suites of the touched areas and the fast gate (minutes); the
frozen-input identity harness and the identity entries at one thread (half an
hour); the identity entries at eight threads and the distributed cases (an
hour); the truth suites and the real-data readouts (hours).

**Level 1, draw identity**: `unit_numerics` (`random draws`): the generator's
fixed vectors (keys and expected bits, in both builds), the unit-interval
mapping's endpoints, independence of order (a permuted key sequence gives the
same values), and the keyed partial permutation's distribution (every subset
of size k equally likely over many keys, and the previous-class-first rule).

**Level 2, frozen-input identity**: a harness on the synthetic references and
particles of `polar_matching_bench` (which already compares the production
kernels at one thread against eight): the fill, the selection and the extras of
each route on the same prepared batch at one and eight threads, two batch
sizes and one and three parts, compared by particle and reference identity;
the selectors on hand-built rows with explicit variates against an independent
cumulative-posterior oracle (ties, invalid and empty rows, state blocks, sparse
eligibility, frozen extras); the input-boundary regression of section 3.2.
These are unit tests in the tables' and the bench's testers.

**Level 3, structural identity**: `test=identity` between a run directory of
the previous phase's build and one of this phase's, same entry, same layout,
same seed, at one thread and at eight, for `simulated_workflow_1jxy`,
`simulated_workflow_6vxx` and `mini_stream_1jxy` (shared memory), and for the
distributed `solve2D` and `solve3D` cases of 2026-11-03 (two parts, one thread
per part; compared with `compare_dirs.py`, paths normalised). Phase 1's build
is the reference for every later phase.

**Level 4, layouts and truth**: on one build, the fixed-input accumulators
(`reconstruct3D` from one committed `ptcl3D` field at `nparts` 1, 2 and 4; the
class averages from one `ptcl2D` field at `nparts` 1, 2 and 4) within the
frozen tolerance; the truth suites unchanged and passing: the two simulated
workflows (shared memory and the 1jxy distributed step), `solve3D_addon`,
`cont_refine3D_1jxy`, `multistate_workflows` (the existing routes), the stream
entries, and `test=nightly` at the end of phases 3 and 4; the real-data
readouts of the 2026-11-03 status against the base build: bgal `solve2D` over
four seeds (median class resolution, best ten, count under 10 A) and trpv1
`refine3D` three iterations (FSC=0.143 and 0.5 shells per iteration), judged
within the seeds' spread as the master-drawn sample was.

**Measured, not gated** (phase 4, on the final build, machine quiet): wall time
per iteration and per phase of the 2D and 3D trpv1 cases of section 2.5 against
the base build at the same layout, three runs each; preparations per particle
(the bench counters: reference loads, canonical and polar preparations, table
bytes written, map and optimiser work); peak RSS of a worker, of the master
and of the node, including at the 3D barrier; the fill's own cost unchanged.

**Base failures.** At this HEAD `simulated_workflow_1jxy` may stop at solve3D
stage 3 (`new_shiftmap`, the map's last shell outside the band, seen on
2026-11-03) and `simulated_workflow_6vxx` fails its final-volume correlation
(0.54 against 0.80). Phase 0 records what the base does. A base failure is not
this run's to fix unless its cause is in this run's files; the identity
comparisons use whatever outputs both runs produce, and a gate that the base
already fails is reported as such, not counted against the phase.

## 13. Phases

### Phase 0: contracts, baselines, the harness

1. Build the tree (`./compile_clean.sh`) and the base export
   (`scratch/keep/src_base`, `scratch/keep/build_base`).
2. The draw-site audit: every `ran3`, `irnd_uni`, `ran_tabu`, `seed_rnd` and
   `random_number` reached from the probabilistic routes (the fills, the
   selectors, the matchers' prob strategies, the masters, the sample draw, the
   diagnostics), listed in Progress with its family of section 4.3 or the
   reason it is outside the contract.
3. The generator: `rnd_keyed`, the key type and packing, the keyed partial
   permutation, in `simple_rnd.f90`; the level-1 tests; the test vectors
   recorded in Progress (Debug and Release).
4. The baselines on the base build: the level-3 run directories (one and eight
   threads; shared memory and the distributed cases); the level-4 layout
   differences of the base code (the fixed-input reconstruction and class
   averages at `nparts` 1, 2, 4) and the tolerance frozen from them; the
   real-data readouts; the bench counters and timings of section 12.
5. The frozen-input harness of level 2 on the current fills and selectors (it
   must show the thread dependence of today's draws, and the identity of the
   deterministic parts).

Exit: both builds pass the fast gate; the audit, the vectors, the tolerance and
the baselines are in Progress; the harness runs on the four routes.

### Phase 1: keyed draws

1. The six site families of section 4.3 keyed; the fills' `seed_rnd` and
   `ran_tabu` removed; the seed drawn once by the project owner when
   `SIMPLE_SEED` is unset and exported; the restart fields of section 4.4.
2. `projs_ns` and the shift prior's width computed once per iteration by the
   caller and passed in (section 3.2), with the input-boundary regression.
3. Level 1 and level 2 pass (the harness now identical across threads, batch
   sizes and parts); the fast gate and `unit_pftc_align2D3D`.
4. The level-4 truth suites and the real-data readouts against the base
   (the trajectory change's validation); the level-3 run directories of this
   build made (one and eight threads; they are the reference for phases 2-4).

Exit: levels 1 and 2 hold; the truth suites pass with their floors unchanged;
the readouts are within the seeds' spread; the reference run directories are
kept under `scratch/keep/identity/phase1/`. Progress records the audit's
closure, the vectors and the readouts.

### Phase 2: the pipeline in shared memory

1. The alignment context and the batch operations of section 3.3 inside
   `refine3D_exec` and `refine2D_exec`; the route adapters over the existing
   fills and selection rules; the decision schema and the output view; the
   budgeted batch capacity.
2. The shared-memory strategies call the pipeline on their own `(params,
   build)`; the commander detour, the second toolbox, the segment write and
   reread and the table and assignment files go from the shared-memory path.
   The distributed path still runs the old two phases on the old objects.
3. Level 3 against phase 1 in shared memory at one and eight threads; level 2;
   the fast gate; `mini_stream_1jxy`; the chunk tests; the distributed cases
   unchanged against phase 1.

Exit: shared memory byte-identical to phase 1 at both thread counts on the
declared outputs; the distributed path unchanged; memory per process measured
and within the budget.

### Phase 3: the distributed parts on the pipeline

1. The masters' iterations drop the probabilistic step and the second
   scheduling phase; the parts run the pipeline; the readouts of section 5
   merged in part order; the sample draw keyed on the master.
2. Level 3 against phase 1 for the distributed cases at one thread per part;
   level 4 across layouts (shared memory; 2, 3, 4 parts); `mini_stream_1jxy`,
   `solve3D_addon`, `cont_refine3D_1jxy`, the stream entries,
   `multistate_workflows`, `test=nightly`.

Exit: the distributed cases byte-identical to phase 1 at their layout; level 4
within the frozen tolerance with the truth gates met; the nightly green.

### Phase 4: retirement, measurements, documents

1. The programs, commanders, files, objects and bench families of section 6
   retired; the generated indexes; `scripts/check_test_registry.py`.
2. The measurements of section 12 on the final build against the base build.
3. The documents of section 14; the layering policy's section 8; the status
   section of `speed_architecture_proposals.md`; skill text proposed in the
   report, not edited.
4. `test=nightly` on the final tree.

Exit: the nightly green; the measurements and the Mac vectors to check are in
Progress; the documents updated.

### Phase 5: the report

The report beside the completed note
(`doc/implementation_notes/completed/prob_align_batch_pipeline_report.md`),
leading with what was learned, then the evidence, then the measurements, then
the skill text and the Mac check for the maintainer, ending with a commit
message; the note moved to `completed/` with `mv`.

## 14. Common rules

- Nothing is committed; git is read-only; the base is a `git archive` export.
- Fix as you go: a defect found is fixed when found and recorded with its
  evidence; a base failure outside this run's files is recorded, not fixed.
- Thresholds are frozen before the runs they judge, never loosened; no run is
  repeated until it passes.
- New or changed OpenMP code gets a tester at 16 threads, repeated, compared
  bit for bit with one thread; no string temporaries inside parallel regions;
  work arrays sized by the team in use.
- No `goto`; THROW macros never continued with `//&`; no `.and.` short circuit
  on allocatable arrays; `use ..., only:`; no large type as a plain component
  (`doc/policies/compile_time_policy.md`); no new parameter (ruling 2); no
  generator or hash that relies on integer overflow (section 4.1).
- Timed work runs with the load recorded; a measurement under load is repeated.
- Bulky outputs are deleted once their numbers are in Progress;
  `scratch/keep/` holds the base export and build, the identity run
  directories of phase 1 and the measurements' logs.
- Write prose for people: Progress entries and the report read without other
  documents, with SIMPLE terms explained and acronyms spelled out on first use;
  name things by what they do.

Documents updated in phase 4: `doc/policies/3D/refine3D_policy.md` (section 8:
the pre-step, the part files, the assignment), `doc/policies/2D/solve2D_policy.md`
(the table contract, the part files, the project-segment table),
`doc/policies/importance_sampling_fractional_update_policy.md`,
`doc/policies/layering_and_project_ownership_policy.md` (sections 2, 3, 8),
`doc/policies/particle_read_path_policy.md` (the buffer's hand-over),
`doc/policies/test_environment_policy.md` (the identity entries at two thread
counts), `doc/algorithms/sampling_and_fractional_updates.md`,
`doc/algorithms/refine2d_class_averaging.md`, `doc/algorithms/solve2d.md`,
`speed_architecture_proposals.md`. Skill text (`simple-refine3d` and its code
map, `simple-solve2d`, `simple-solve3d-importance-sampling` and its coupling
map, `simple-frac-update-trailing`, `simple-main-commanders`,
`simple-main-strategies`, `simple-main-apis`) is proposed in the report.

### File table

| file | why | phase |
| --- | --- | --- |
| `doc/implementation_notes/planned/prob_align_batch_pipeline.md` | this plan and its progress | each |
| `src/utils/math/simple_rnd.f90` | the keyed generator, the key type, the keyed partial permutation | 0 1 |
| `src/utils/math/simple_rnd_tester.f90` | the level-1 tests and vectors | 0 1 |
| `src/main/commanders/test/simple_polar_matching_bench.f90` | the frozen-input harness (fill, selection, extras per route; threads, batch sizes, parts) | 0 1 2 |
| `src/main/commanders/test/simple_commanders_test_highlevel_matchbench.f90` | the harness's suites and reports | 0 1 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | the selector's keyed draw; the readout's global rows; `calc_num2sample` and `shift_prior_sigma2` as model inputs | 1 2 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | the keyed fill; the per-batch row store and selection; the retired global table | 1 2 3 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | the same for the neighbourhood route | 1 2 3 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90` | the same for the 2D routes | 1 2 3 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90` | the selectors against the oracle; the input-boundary regression | 0 1 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D_tester.f90` | the same, 2D | 0 1 2 |
| `src/main/pftc/simple_polarft_corr.f90` | the in-plane draw keyed | 1 |
| `src/main/pftc/simple_polarft_shiftmap.f90` | the disk draw keyed | 1 |
| `src/main/pftc/simple_polarft_calc.f90` | the key scope carried to the draw sites | 1 |
| `src/main/ori/simple_oris_sampling.f90` | the sample draw keyed | 1 |
| `src/main/ori/simple_oris.f90` | its interface | 1 |
| `src/main/ori/simple_oris_tester.f90` | its tests | 1 |
| `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` | the sample draw's scope; the master-drawn sample | 1 3 |
| `src/main/strategies/search/simple_strategy3D_matcher.f90` | the alignment context and batch loop (3D) | 2 3 4 |
| `src/main/strategies/search/simple_strategy2D_matcher.f90` | the alignment context and batch loop (2D) | 2 3 4 |
| `src/main/strategies/search/probabilistic/simple_strategy3D_prob.f90` | the decision schema's consumer (3D) | 2 |
| `src/main/strategies/search/probabilistic/simple_strategy2D_prob.f90` | the decision schema's consumer (2D) | 2 |
| `src/main/strategies/search/simple_strategy3D_srch.f90` | `prep4prob` from the model, not the field | 2 |
| `src/main/strategies/search/simple_strategy2D_srch.f90` | the same, 2D | 2 |
| `src/main/strategies/search/simple_matcher_ptcl_batch.f90` | the batch preparation by consumer | 2 |
| `src/main/strategies/search/simple_matcher_ptcl_io.f90` | the buffer's hand-over retired | 2 4 |
| `src/main/strategies/parallelization/simple_prob_align_strategy.f90` | the masters' step retired from the shared-memory path (2) and the distributed path (3); the module retired (4) | 2 3 4 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | the shared-memory strategy on the pipeline (2); the distributed master's one phase (3) | 1 2 3 |
| `src/main/strategies/parallelization/simple_refine2D_strategy.f90` | the same, 2D | 1 2 3 |
| `src/main/strategies/parallelization/simple_strategy_sigma2.f90` | the residual's consumer flag | 2 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | the thin wrappers (2, 3); retired (4) | 2 3 4 |
| `src/main/commanders/simple/simple_commanders_refine3D.f90` | the worker program on the pipeline; the seed | 1 3 |
| `src/main/commanders/simple/simple_commanders_refine2D.f90` | the same, 2D | 1 3 |
| `src/main/apis/simple_private_exec_api.f90` | the retired programs | 4 |
| `production/simple_private_exec_driver.f90` | the retired programs | 4 |
| `src/defs/simple_defs_fname.f90` | the retired file names | 4 |
| `src/main/params/simple_parameters_phases.f90` | the restart fields read with `continue`; no new key | 1 |
| `src/main/project/simple_sp_project.f90` | the restart fields in `projinfo` | 1 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | suites registered | 0 1 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | their help | 0 1 |
| `src/main/commanders/test/simple_test_registry.f90` | an entry's row if its outputs change | 3 4 |
| `doc/policies/3D/refine3D_policy.md` | documents | 4 |
| `doc/policies/2D/solve2D_policy.md` | documents | 4 |
| `doc/policies/importance_sampling_fractional_update_policy.md` | documents | 4 |
| `doc/policies/layering_and_project_ownership_policy.md` | documents | 4 |
| `doc/policies/particle_read_path_policy.md` | documents | 4 |
| `doc/policies/test_environment_policy.md` | documents | 4 |
| `doc/algorithms/sampling_and_fractional_updates.md` | documents | 4 |
| `doc/algorithms/refine2d_class_averaging.md` | documents | 4 |
| `doc/algorithms/solve2d.md` | documents | 4 |
| `doc/implementation_notes/planned/speed_architecture_proposals.md` | the status of this work | 4 |
| `doc/code_overview/fortran-indexes/` | regenerated | 4 |
| `doc/implementation_notes/completed/prob_align_batch_pipeline_report.md` | the report | 5 |
| `doc/implementation_notes/completed/prob_align_batch_pipeline.md` | the plan, moved when the run completes | 5 |

Other files a phase needs are added to this table, with the reason, by the
phase that needs them, before they are edited.

## Progress

(none yet)
