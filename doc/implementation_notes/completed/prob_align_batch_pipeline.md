# Probabilistic alignment as a batch pipeline: per-particle draws, keyed random streams, no global table

Planned 2026-11-04 on master `a723f4919`, ruled done by Hans the same day ("we
have reached the point of no return; this refactoring is obvious"), reviewed by
Hans against `6b924ccf5` the same day (the review's requirements are folded into
sections 3, 4, 5 and 8), and planned as the autonomous run `prob_pipeline` on the
Dell (sections 10-14). **Implemented** by the run `prob_pipeline` (2026-11-04/05; the Progress entries below and `prob_align_batch_pipeline_report.md` beside this note). The frame is
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
| the sample (`sample4update_cnt`, `_rnd`, `_class`, `_fillin`) | `SAMPLE` | particle 0; the sampler's draw ordinal (its shuffle count) in the reference lane, the swap's ordinal in the ordinal lane (phase 1: the keyed partial permutation keys each swap on the ordinal lane, so a sampler with several shuffles counts them in the reference lane) |
| the in-plane draw per scored reference (`draw_posterior_crvec`) | `INPL` | the reference |
| the disk draw (`shiftmap_draw`) | `DISK` | the screen's ordinal for the particle |
| the sparse fills' candidate picks | `CAND` | a keyed partial Fisher-Yates permutation of the candidate set (sampling without replacement, as `ran_tabu` gives today), one keyed uniform per swap (the swap's ordinal), the previous class or reference last as today (`put_last`; the 3D stochastic fill then forces it in as the incumbent), the candidate count and stopping rules unchanged |
| the assignment draw (`draw_posterior_row` in the selector) | `ASSIGN` | ordinal 0 (the joint multi-state draw is one draw); the 2D search's fallback class pick (a populated class when the previous one is empty or the run starts fresh) draws under this family with the pick's ordinal from 1 |
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

Fixed by Hans on 2026-11-04 (the three the launch waited on, R1-R3, were
answered the same day; the run's `spec/rulings.md` holds the rulings' text):

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
7. **R1, the order against `ms_posterior2`** (fixed by Hans on 2026-11-04,
   "I'm good with all this", as proposed). Both developments edit
   `simple_eul_prob_tab.f90`, `_utils`, `_neigh`, `simple_prob_align_strategy.f90`,
   `simple_refine3D_strategy.f90` and `simple_commanders_prob.f90`. This run
   goes first; `ms_posterior2`'s Phase 0 (the `labels` route, P2S: test and
   data files only) may run on dell2 meanwhile, and its Phase 1 onward starts
   on the tree this run leaves, with the joint draw as a selection rule of
   section 3.3. Both machines are available to the maintainer's runs; this
   run uses the Dell only.
8. **R2, the generator and key schema** of section 4 as written (lanes,
   families, the restart fields), fixed by Hans on 2026-11-04; the exact
   constants and test vectors were produced in phase 0 (its Progress entry
   records them), are pinned by `simple_rnd_tester` and are never changed
   afterwards in this run.
9. **R3, the `labels` route**, fixed by Hans on 2026-11-04: it does not exist
   at this HEAD; the truth gates of this run are the existing asserted routes
   of `multistate_workflows` (`states`, `cavgs`, `shc`, `consensus`) and the
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
| `doc/implementation_notes/completed/prob_align_batch_pipeline.md` | this plan and its progress | each |
| `src/utils/math/simple_rnd.f90` | the keyed generator, the key type, the keyed partial permutation | 0 1 |
| `src/utils/math/simple_rnd_tester.f90` | the level-1 tests and vectors | 0 1 |
| `src/main/commanders/test/simple_polar_matching_bench.f90` | the frozen-input harness (fill, selection, extras per route; threads, batch sizes, parts); its reference route (worker tables, part files, the master's merge and assignment) retired with that code, the pipeline route taking over the part and perturbation dimensions (added to phase 4) | 0 1 2 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_matchbench.f90` | the harness's suites and reports; the configurations of the pipeline route alone (added to phase 4) | 0 1 2 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | the selector's keyed draw; the readout's global rows; `calc_num2sample` and `shift_prior_sigma2` as model inputs | 1 2 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | the keyed fill; the per-batch row store and selection; the retired global table | 1 2 3 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | the same for the neighbourhood route | 1 2 3 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90` | the same for the 2D routes | 1 2 3 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90` | the selectors against the oracle; the input-boundary regression; the part-file and global-table tests retired with that code, the selection tests on the pipeline's table (added to phase 4) | 0 1 2 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D_tester.f90` | the same, 2D | 0 1 2 4 |
| `src/main/pftc/simple_polarft_corr.f90` | the in-plane draw keyed; a comment naming the retired table worker (4) | 1 4 |
| `src/main/pftc/simple_polarft_shiftmap.f90` | the disk draw keyed | 1 |
| `src/main/pftc/simple_polarft_calc.f90` | the key scope carried to the draw sites | 1 |
| `src/main/pftc/simple_polarft_core.f90` | the scope's setter and its reset in `kill` (the calculator's bodies live in its submodules; added by phase 1) | 1 |
| `src/main/pftc/simple_polarft_corr_tester.f90` | its posterior-draw tests set the calculator's scope and vary it per draw (added by phase 1) | 1 |
| `src/main/pftc/simple_polar_shiftmap_tester.f90` | its disk-draw test sets the scope and gives each draw its ordinal (added by phase 1) | 1 |
| `src/main/pftc/simple_pftc_inplane_tester.f90` | its fixture sets the scope (added by phase 1) | 1 |
| `src/main/ori/simple_oris_sampling.f90` | the sample draw keyed | 1 |
| `src/main/ori/simple_oris.f90` | its interface | 1 |
| `src/main/ori/simple_oris_tester.f90` | its tests | 1 |
| `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` | the sample draw's scope; the master-drawn sample | 1 3 |
| `src/main/strategies/search/simple_strategy3D_matcher.f90` | the alignment context and batch loop (3D) | 2 3 4 |
| `src/main/strategies/search/simple_strategy2D_matcher.f90` | the alignment context and batch loop (2D) | 2 3 4 |
| `src/main/strategies/search/probabilistic/simple_strategy3D_prob.f90` | the decision schema's consumer (3D) | 2 |
| `src/main/strategies/search/probabilistic/simple_strategy2D_prob.f90` | the decision schema's consumer (2D); its header's mention of the retired table (4) | 2 4 |
| `src/main/strategies/search/simple_strategy3D_srch.f90` | `prep4prob` from the model, not the field | 2 |
| `src/main/strategies/search/simple_strategy2D_srch.f90` | the fallback class pick of the 2D search keyed (1, found by the audit); `prep4prob` from the model, not the field (2) | 1 2 |
| `src/main/strategies/search/simple_matcher_ptcl_batch.f90` | the batch preparation by consumer | 2 |
| `src/main/strategies/search/simple_strategy2D_alloc.f90` | the classic search order of the 2D batch preparation not drawn under the probabilistic routes, which never read it (added by phase 2: the audit's matcher item) | 2 |
| `src/main/strategies/search/simple_matcher_ptcl_io.f90` | the buffer's hand-over retired | 2 4 |
| `src/main/strategies/parallelization/simple_prob_align_strategy.f90` | the masters' step retired from the shared-memory path (2) and the distributed path (3); the module retired (4) | 2 3 4 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | the shared-memory strategy on the pipeline (2); the distributed master's one phase (3) | 1 2 3 |
| `src/main/strategies/parallelization/simple_refine2D_strategy.f90` | the same, 2D | 1 2 3 |
| `src/main/strategies/parallelization/simple_strategy_sigma2.f90` | the residual's consumer flag | 2 |
| `src/main/commanders/simple/simple_commanders_prob.f90` | the table workers compute the iteration's model once and pass it to the fills (1, the fills' caller today); the thin wrappers (2, 3); retired (4) | 1 2 3 4 |
| `src/main/commanders/simple/simple_commanders_refine3D.f90` | the worker program on the pipeline; the seed; the call to the retired part-file cleanup (4) | 1 3 4 |
| `src/main/commanders/simple/simple_commanders_refine2D.f90` | the same, 2D | 1 3 |
| `src/main/apis/simple_private_exec_api.f90` | the retired programs | 4 |
| `production/simple_private_exec_driver.f90` | the retired programs | 4 |
| `src/defs/simple_defs_fname.f90` | the retired file names | 4 |
| `src/main/params/simple_parameters_phases.f90` | the restart fields read with `continue`; no new key | 1 |
| `src/main/project/simple_sp_project.f90` | the restart fields in `projinfo` | 1 |
| `src/main/project/simple_sp_project_core.f90` | their bodies (the `projinfo` accessors live in this submodule; added by phase 1) | 1 |
| `src/main/project/simple_sp_project_tester.f90` | the restart fields' write-and-read round trip (added by phase 1) | 1 |
| `src/utils/qsys/simple_qsys_ctrl.f90` | the queue scripts export `SIMPLE_SEED` when it is set, so a part run by the stream's persistent worker, whose environment is not the master's, reads the run's seed (added by phase 1) | 1 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | suites registered; the 2D table suite's name without the retired part-file I/O (phase 4) | 0 1 4 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | their help; the 2D table sub-suite's name without the retired part-file I/O (4) | 0 1 4 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | the help of `polar_matching_bench`'s `suite=` names the fill harness (added by phase 0: the bench's help lives here, not in the class file); the harness's description without the retired route (phase 4) | 0 4 |
| `src/main/commanders/test/simple_test_registry.f90` | an entry's row if its outputs change | 3 4 |
| `src/main/simple_refine3D_states_controller.f90` | the call to the retired part-file cleanup `del_prob_align_files` (added by phase 4) | 4 |
| `src/main/solve/simple_solve3D_utils.f90` | the same call after each refine3D stage (added by phase 4) | 4 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | the same call in the multi-state missing update (added by phase 4) | 4 |
| `src/main/simple_external_reference_pose_initialization.f90` | its deleter of the retired part and assignment files (added by phase 4) | 4 |
| `src/main/commanders/simple/simple_commanders_solve2D.f90` | its final cleanup of the retired part and assignment files (added by phase 4) | 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_cont.f90` | its check that a Cartesian pass leaves no table: the retired file names replaced by the readout records (added by phase 4) | 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_workflow.f90` | a comment naming the retired master step (added by phase 4) | 4 |
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
| `doc/policies/2D/particle_cache_policy.md` | its mentions of the retired 2D table worker and step (added by phase 4) | 4 |
| `doc/policies/3D/solve3D_addon_policy.md` | its mentions of the retired table workers and step (added by phase 4) | 4 |
| `doc/algorithms/refine3d.md` | its candidate-preparation step, now inside the search (added by phase 4) | 4 |
| `doc/implementation_notes/README.md` | the pipeline row's status (added by phase 4) | 4 |
| `doc/code_overview/fortran-indexes/` | regenerated | 4 |
| `doc/implementation_notes/completed/prob_align_batch_pipeline_report.md` | the report | 5 |
| `doc/implementation_notes/completed/prob_align_batch_pipeline.md` | the plan, moved when the run completes | 5 |

Other files a phase needs are added to this table, with the reason, by the
phase that needs them, before they are edited.

## Progress

### Phase 0: contracts, baselines, the harness (2026-11-04, the Dell)

Done on 2026-11-04 by the autonomous run `prob_pipeline` on the Dell (Oracle Linux 8.10, GCC
15.2.1, CMake 3.26.5; run root `/data2/agent_runs/prob_pipeline`, linked as
`~/agent_runs/prob_pipeline`; every path below that is not absolute is relative to it). Base HEAD
`22a643642` (clean), exported with `git archive` to `scratch/keep/src_base` and built in
`scratch/keep/build_base`; the tree built in `repo/build` (Release) and `scratch/keep/build_debug`
(Debug, `-O0 -fcheck=all`), all three with tests on, `-DPython3_EXECUTABLE=/usr/bin/python3.12`,
NICE off and the architecture optimisation on, as the maintainer's checkout is configured. The
journal of the phase is `state/journal_phase_0.md`.

**Files changed** (the file table's rows of phase 0, plus one row added by it):
`src/utils/math/simple_rnd.f90` (the keyed generator), `src/utils/math/simple_rnd_tester.f90`
(its level-1 tests and vectors), `src/main/commanders/test/simple_polar_matching_bench.f90` (the
frozen-input harness), `src/main/commanders/test/simple_commanders_test_highlevel_matchbench.f90`
(its suite `fills`), `src/main/ui/simple_test/simple_test_ui_highlevel.f90` (the bench's help
names the suite; the row was added to the table because the bench's help lives there, not in the
class file), and this plan (section 10's rulings recorded as fixed, the file-table row, this
entry). The table testers `simple_eul_prob_tab_tester.f90` and `_2D_tester.f90` and the class
registration files were not changed: the selectors' oracle tests need the explicit variate that
phase 1 gives `draw_posterior_row`, and the harness is a suite of an existing test program, so
nothing new is registered.

**The draw-site audit** (item 2). Every random draw reachable from the probabilistic routes
(`refine=prob` and `refine=prob_neigh` in `refine3D` and `refine2D`), traced from the primitive
to the route through its callers. "Route" names which of the four routes reach the site; a
site reached by none of them is outside the contract and says why. The primitives: `ran3` is
the intrinsic `random_number` on the calling thread's stream; `irnd_uni` draws one double on
that stream; `ran_tabu%ne_ran_iarr` draws a full permutation by rejection (a varying number of
draws); `seed_rnd` reseeds the calling thread's stream (`SIMPLE_SEED` plus 7919 times the
process's call count when the seed is set).

| file | line | routine | primitive | route(s) | family of section 4.3, or outside |
| --- | --- | --- | --- | --- | --- |
| `src/main/ori/simple_oris_sampling.f90` | 227 | `sample4update_cnt` | `partial_shuffle` at the cutoff tier when the fill-in is a subset | all four (the default fractional-update sampler; the master's draw or the shared-memory strategy's) | `SAMPLE` |
| same | 177 | `sample4update_rnd` | `partial_shuffle` | 2D dense and 2D sparse (the first sample of a 2D run with `update_frac`; also classic 2D) | `SAMPLE` |
| same | 287 | `sample4update_class` | `shuffle` per update-count tier and class unit | 3D dense and 3D neighbourhood (`balance=class` or `cavg`) | `SAMPLE`; it also shuffles tiers taken whole, draws that change nothing |
| same | 578 | `sample4update_fillin` | `shuffle` per tier | 3D dense and 3D neighbourhood (every fifth iteration under the fill-in) | `SAMPLE`; the same wasted draws |
| same | 21 | `select_particles_set` | `partial_shuffle` | no caller in `src/` | outside: dead code, to delete in phase 1 (its file) |
| same | 657, 687, 719 | `sample_balanced_1/_2`, `sample_balanced_inv` | `shuffle` | the `sample_classes` program only | outside: not on a refinement path |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | 302 | `fill_tab_range` | `seed_rnd` at every call, that is every batch | 3D dense | goes in phase 1 (section 4.3) |
| `.../simple_eul_prob_tab_neigh.f90` | 197, 451 | `fill_tab_stoch_range`, `fill_tab_subspace_range` | `seed_rnd` per call | 3D neighbourhood | goes in phase 1 |
| `.../simple_eul_prob_tab2D.f90` | 282, 528 | `fill_tab_range`, `fill_tab_prob_neigh_range` | `seed_rnd` per call | 2D dense, 2D sparse | goes in phase 1 |
| `.../simple_eul_prob_tab_neigh.f90` | 205, 319 | `fill_tab_stoch_range` (`direct_rts` per thread; `ne_ran_iarr` per particle inside the dynamic-schedule region) | `ran_tabu` | 3D neighbourhood, modes `shc` and `snhc` | `CAND`: the keyed partial Fisher-Yates permutation of the full reference set; today the previous reference is put last of the order (`put_last`, line 324) and forced in afterwards as the incumbent (365-376); the count bound and the stopping rules stay |
| `.../simple_eul_prob_tab2D.f90` | 539, 609 | `fill_tab_prob_neigh_range` (`direct_rts`; `ne_ran_iarr` per particle) | `ran_tabu` | 2D sparse | `CAND`: the previous class put last (`put_last`, 613), the first `nrefs_bound` of the order scored, no incumbent |
| `.../simple_eul_prob_tab_utils.f90` | 987 | `draw_posterior_row`, through `draw_posterior_crvec` (`simple_polarft_corr.f90` 735) from `gen_posterior_val` (662) and `gen_posterior_val_bank` (710) | `ran3` | all four fills: `simple_eul_prob_tab.f90` 863 and 871, `_neigh` 1463 and 1466, `_2D` 471, 863 and 866, every scored reference or class of every particle, inside the dynamic-schedule regions | `INPL`, keyed on the reference |
| `.../simple_eul_prob_tab_neigh.f90` | 1066 | `seed_fallback_if_empty` inside `ref_assign_neigh`, through `gen_posterior_val` | `ran3` | 3D neighbourhood master, a particle without a valid candidate | `INPL` on the fallback reference |
| `src/main/pftc/simple_polarft_shiftmap.f90` | 607 | `shiftmap_draw` with `l_greedy=.false.` | `ran3` | none: production calls `shiftmap_mode` (652-663, greedy) from `simple_eul_prob_tab.f90` 557, `_neigh` 1541 and 1580, `_2D` 416, 726 and 763 and `simple_shiftmap_srch.f90` 152; the draw is reached by `simple_polar_shiftmap_tester.f90` only | `DISK`: no live production site at this HEAD; phase 1 keys the site so that the family is ready and records that it is tester-only |
| `.../simple_eul_prob_tab_utils.f90` | 987 | `draw_posterior_row` from `ref_assign` (`simple_eul_prob_tab.f90` 1115) | `ran3` | 3D dense master, particle order | `ASSIGN` |
| same | 987 | from `ref_assign_neigh` (`_neigh` 828) | `ran3` | 3D neighbourhood master | `ASSIGN` |
| same | 987 | from the 2D `ref_assign` (`_2D` 1223) | `ran3` | 2D dense master | `ASSIGN` |
| same | 987 | from `ref_assign_sparse` (`_2D` 1083) | `ran3` | 2D sparse master | `ASSIGN` |
| `src/main/strategies/search/simple_strategy2D_srch.f90` | 188-197 | `prep4srch`: a class drawn until a populated one when the previous class is empty (population below 2) or the run starts fresh | `irnd_uni` | 2D dense and 2D sparse (`strategy2D_prob` 43; the drawn class sets `prev_corr` and is the fallback of `assign_ori` when the assigned class is empty), also every classic 2D strategy | inside the contract: a per-particle pick on the probabilistic path, `CAND` with the pick's ordinal; its file joins the table in phase 1 |
| `src/main/strategies/search/simple_strategy2D_alloc.f90` | 144-160 | `prep_strategy2D_batch`: a search order per particle, every batch | `ran_tabu` | the 2D matcher's every batch (171), so the 2D prob routes draw it, but only the `snhc*` strategies read it | outside: a classic-mode order drawn and unread on the probabilistic routes (it only advances the worker's intrinsic stream, which the keyed draws no longer use); skipping it under `prob` is a matcher item of phase 2 |
| `src/main/strategies/search/simple_strategy3D_matcher.f90` | 519 | `choose_and_run_strategy` | `ran3 < GREEDY_FREQ` | `refine=shc` only | outside: classic |
| `src/main/strategies/search/simple_strategy3D_alloc.f90`, `simple_strategy3D_srch.f90` | 88-90; 221-226 | `prep_strategy3D` (the per-thread tabu samplers), `prep4srch` | `ran_tabu` | the classic 3D strategies (`prep_strategy3D` skips them for `prob` and `prob_neigh`, 78-79; `strategy3D_prob` calls `prep4prob`, which draws nothing) | outside: classic |
| `.../simple_eul_prob_tab_utils.f90` | 1353, 1386, 1433 | `angle_sampling_2`, `sample_bounded_dists`, `sample_power_dists` | `ran3` | `refine=neigh`, `greedy_sub`, the `snhc_smpl*` and `inpl_smpl` strategies through `gen_prob_euclid_val` and `gen_prob_power_euclid_val` | outside: classic |
| `src/main/strategies/parallelization/simple_refine2D_strategy.f90` | 855-856 | `init_standard_refs`, `cls_init=randcls` | `irnd_uni`, `ran3` | a fresh `refine2D` without references (`init_refine2D_refs` 171 and 326) | outside: a workflow's first-stage initialisation |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | 908-931 | `seed_multistate_startup_labels` through `gen_labelling` (`simple_cluster_seed.f90` 26) | reseed and uniform labels | a multi-state start (`nstates > 1`, `startit = 1`, one state bin), reached on the probabilistic routes (1159-1161, 1598) | outside: first-stage initialisation |
| `src/main/ori/simple_oris_setters.f90` | 690 | `rnd_cls` | `irnd_uni` | `solve3D` with a supplied volume and a class balance (`simple_commanders_solve3D.f90` 1145) | outside: first-stage initialisation |
| `src/main/ori/simple_ori.f90`, `simple_ori_utils.f90` | 660-662; 203-205 | `rnd_ori`, `rnd_romat` through `oris%rnd_oris` | `ran3` | a virgin orientation field at a `refine3D` start (`simple_refine3D_strategy.f90` 1168, 1615), `solve3D` starts (484, 1275, 1721), simulation, `make_oris` | outside: first-stage initialisation and simulation |
| `src/main/ori/simple_sym.f90` | 812 | `symrandomize` | `irnd_uni` | `solve3D` between its stages (`symmetrize`), `oriops` | outside: a `solve3D` step, not the refinement |
| `src/main/strategies/parallelization/simple_cls_expansion_strategy.f90` | 1033, 1049 | `split_class_flex` | `seed_rnd_fixed` | the `cls_expansion` program | outside |
| `src/main/ori/simple_ori.f90`, `simple_oris_setters.f90`, `simple_oris_transform.f90`, `simple_oris_reshape.f90`, `simple_oris_getters.f90` | 595-618, 679-705; 544-636; 24-60; 49, 76; 519 | `rnd_euler_3/4`, `rnd_inpl`, `rnd_shift`, `gau_rnd_shift`, `rnd_inpls`, `rnd_ctf`, `rnd_states`, `rnd_lps`, `rnd_corrs`, `rnd_proj_space`, `introd_alig_err`, `introd_ctf_err`, `split_state`, `split_class`, `get_pinds` with a shuffle | `ran3`, `gasdev`, `irnd_uni`, `ran_tabu` | testers, simulation, `make_oris`, `oriops`, `selection`, `make_cavgs` with `remap_cls`; `rnd_inpl`, `rnd_shift`, `rnd_proj_space` and `split_state` have no caller | outside; the four without a caller are dead code (their files are not in this run's table; recorded for the maintainer) |
| `src/main/strategies/search/probabilistic/simple_strategy3D_prob.f90`, `simple_strategy2D_prob.f90`; `src/main/strategies/parallelization/simple_prob_align_strategy.f90`; `src/main/commanders/simple/simple_commanders_prob.f90` | | the search strategies, the masters' steps, the table commanders | none | | no draw of their own: `strategy3D_prob` runs `prep4prob`, the deterministic polish and the commit; `strategy2D_prob` draws only through `prep4srch` and the matcher's batch order above; the masters draw only through the sample and `ref_assign`; the commanders reseed only through `parameters%new` (`simple_parameters_phases.f90` 28), which every commander start does |
| the diagnostics `prob_neff_readout` (`simple_eul_prob_tab_utils.f90` 560-680), `prob_post_diag` (`simple_eul_prob_tab.f90` 620-771, 930-980), `print_disk_set` (`simple_prob_align_strategy.f90` 489-520) | | | none | all four | no draw: the readout's rows are chosen by a fixed stride, the posterior diagnostic evaluates the deterministic kernels, the disk-set line counts; `DIAG` is reserved with no live site |

Two observations the audit adds to section 2.3. The reseed of `parameters%new` runs at every
commander start, so the shared-memory detour of today reseeds the process between the sample's
draws and the assignment's; in a distributed run the master's stream carries from the one to the
other. And the previous reference or class enters the sparse fills' order last, not first: the
plan's `CAND` row ("the previous class or reference first as today") is to be read as the code's
rule, which phase 1 keeps.

**The generator** (item 3; ruling R2 applied here, section 10). `rnd_keyed` in
`src/utils/math/simple_rnd.f90` beside `seed_rnd`: a draw is a pure function of its key, the type
`rnd_key` with the fields of section 4.2 (`seed` as a 64-bit integer; `version`, `scope`,
`iteration`, `particle`, `family`, `reference`, `ordinal` as default integers). The key's nine
32-bit lanes, in this pinned order: the seed's low word, the seed's high word, the version, the
scope, the iteration, the particle, the family, the reference, the ordinal (a negative field
contributes its two's-complement low word). The mixer is lowbias32, two multiply-xor-shift rounds
on a 32-bit lane held in a signed 64-bit integer: `x ^= x >> 16; x = (x * 0x21f0aaad) mod 2^32;
x ^= x >> 15; x = (x * 0x735a2d97) mod 2^32; x ^= x >> 15` (the constants 569420461 and
1935289751, both below 2^31, so a lane times a constant stays below 2^63: no step relies on
integer overflow, which the Fortran standard leaves undefined). Two lanes a and b start at
0x6a09e667 and 0x510e527f (the fractional parts of the square roots of 2 and 11, 1779033703 and
1359893119); lane i of the key is absorbed as `a = mix32(a xor w_i)`, `b = mix32(b xor a xor i)`;
the finaliser is `a = mix32(a xor b)`, `b = mix32(b xor a)`; the 64-bit word is a in the high and b
in the low 32 bits; the unit variate is the word's top 53 bits over 2^53 (0 is possible, 1 is not).
`RND_KEYED_VERSION = 1`; the families `RND_FAM_SAMPLE` to `RND_FAM_DIAG` are 1 to 6 in the order
of section 4.3; a scope is `rnd_scope(program, pass) = 16 * program + pass` with the programs
`refine3D = 1`, `refine2D = 2` and the pass kinds alignment 1, `pose_cont` polish 2, final
full-population assignment 3, diagnostic 4. Beside the variate: `rnd_keyed_bits` (the word),
`rnd_keyed_int(key, n)` (the uniform integer in 1..n by the variate's equal-width bin, as
`irnd_uni`), and `rnd_keyed_partial_shuffle(key, arr, nshuffle)`, the keyed partial Fisher-Yates
permutation whose swap i draws its uniform under the key with `ordinal = i` (so the permutation
depends on the key alone, never on the order of other draws).

The ten pinned vectors, produced by an independent Python reference (unbounded integers, every step
exact; `scratch/keep/rnd/rnd_keyed_ref.py`) and reproduced bit for bit by the Release
(`-O3 -march=native`) and the Debug (`-O0 -fcheck=all`) builds of GCC 15.2.1 on the Dell
(`scratch/keep/rnd/vectors_{python,release,debug}.txt`); the Mac check runs `unit_numerics`:

| seed | version | scope | iteration | particle | family | reference | ordinal | 64-bit word (signed) | unit variate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | -4131892048949355934 | 7.76009683202676270E-01 |
| 20261103 | 1 | 17 | 1 | 1 | 2 | 1 | 0 | -3402735607585271042 | 8.15537333093113270E-01 |
| 20261103 | 1 | 17 | 1 | 2 | 2 | 1 | 0 | 4682640351440764835 | 2.53846442100018055E-01 |
| 20261103 | 1 | 17 | 1 | 1 | 2 | 2 | 0 | -2904353384885242822 | 8.42554687522089507E-01 |
| 20261103 | 1 | 17 | 2 | 1 | 2 | 1 | 0 | 6790283504847523914 | 3.68102006387408531E-01 |
| 20261103 | 1 | 17 | 1 | 750000 | 5 | 0 | 0 | -4401584271307322743 | 7.61389638533528701E-01 |
| -1 | 1 | 33 | 3 | 12345 | 3 | 0 | 7 | -2048942294402738840 | 8.88926615655555752E-01 |
| 4611686018427400249 | 1 | 18 | 1 | 1 | 4 | 20000 | 2147483647 | 3588558501286480495 | 1.94536146159306367E-01 |
| 123456789012345 | 1 | 17 | 25 | 35645 | 1 | 0 | 35645 | -3713091948123483349 | 7.98712882160304227E-01 |
| 20261103 | 1 | 17 | 1 | 1 | 6 | 0 | 1 | -5700306898616924037 | 6.90985743834271138E-01 |

The level-1 tests in `simple_rnd_tester.f90` (`unit_numerics`, suite "random draws", 123
assertions, both builds): the vectors; the mapping's endpoints (the zero word gives 0, the all-one
word 1 - 2^-53, the low eleven bits are dropped, the lanes follow the fields); independence of
order (a key sequence drawn in reverse, or interleaved with intrinsic draws and reseeds, gives the
same variates, and the keyed draws leave the intrinsic stream untouched); sensitivity (changing any
field of the key, either seed word included, changes the word; 100,000 consecutive particle keys
give distinct words, a uniform distribution over 64 bins and no lag-one serial correlation); the
integer draw's bins around 2^24 and at the largest integer; and the keyed partial permutation (a
duplicate-free prefix of a permutation of the input, the same under the same key whatever was drawn
before, every subset of size three of six equally likely over 20,000 keys and every order of a
permutation of four over 24,000, the key's own ordinal overridden by the swap index). The first
transcription of the two multipliers to decimal was wrong and the vectors caught it: the Python
reference used the hexadecimal constants, the Fortran decimals 569682605 and 1935093143 were
corrected to 569420461 and 1935289751, and every other keyed test had passed on the wrong
constants, which is what the pinned vectors are for.

**The level-3 baselines on the base build** (item 4, first part). Run directories under
`scratch/base_runs/` of the run root `/data2/agent_runs/prob_pipeline` (the symbolic link
`~/agent_runs/prob_pipeline`), every run with `SIMPLE_SEED=20261103` in the environment, the
timeline in `scratch/base_runs/timeline.txt`; the entries ran as `simple_test_exec` runs with
`nthr` on the command line and `OMP_NUM_THREADS` set to it, the distributed cases as `simple_exec`
runs on copies of the 2026-11-03 base run's projects (`run_base_level3.sh`).

| run | layout | result | time |
| --- | --- | --- | --- |
| `simulated_workflow_1jxy` | one thread; "eight threads" | both stop at solve3D stage 3 with `the map's last shell must lie in the band; new_shiftmap` (the base failure of section 12), exit 1, no `FAIL:` line before it | 242 s; 257 s |
| `simulated_workflow_6vxx` | one thread; "eight threads" | both complete with the one known `FAIL: final-volume correlation is below the required minimum` (the base failure of section 12), exit 0 | 787 s; 840 s |
| `mini_stream_1jxy` | eight threads | passes (`PASS: mini_stream 1jxy validated 6 micrographs, 21 particles, and 2 classes`) | 56 s |
| `solve2D` of the 1jxy workflow's extracted project, `nparts=2 nthr=1 ncls=4` | two parts | normal stop | 12 s |
| the same with `nsample=30` (the fractional variant) | two parts | normal stop | 9 s |
| `solve3D` of the 6vxx workflow's solve2D project, `pgrp=c3 nparts=2 nthr=1` | two parts | normal stop | 350 s |

The "eight-thread" runs of the simulated workflows exercise one thread: the test fixes its
thread count (`NTHR = 1`, `simple_commanders_test_highlevel_workflow.f90` line 54) and calls
`omp_set_num_threads(1)` before its first commander (line 152), and its registry row says it does
not scale with `nthr`. `test=identity` between the one- and the "eight-thread" directories
confirms it: `simulated_workflow_6vxx` 184 declared files the same, none differing;
`simulated_workflow_1jxy` 143 the same, and the one item reported as differing or missing is the
`solve3D_manifest.txt` pattern, which neither run wrote because the base failure stops solve3D
before the manifest (`scratch/base_runs/identity_*_t1_vs_t8/identity_report.txt`). The level-3
rung at eight threads therefore checks the stream entry only, where `mini_stream` passes
`params%nthr` to its simulation and reprojection commanders and runs its stream commanders at one
thread. A multi-threaded probabilistic alignment is covered at level 3 by nothing with declared
outputs at this HEAD; phase 1 keeps, as its reference directories, the shared-memory bgal
`solve2D` and a shared-memory trpv1 `refine3D` at eight threads beside the entries, compared with
`compare_dirs.py` (section 12's tooling for the cases that are not entries).

**Base failures** (section 12): the two the plan names, and no other. `simulated_workflow_1jxy`
stops at solve3D stage 3 in `new_shiftmap` ("the map's last shell must lie in the band") at both
thread counts; `simulated_workflow_6vxx` fails its final-volume correlation gate at both. The
identity comparisons of the later phases use the outputs both runs produce before the stop, and
the two gates are reported as base failures, not counted against a phase.

**The level-4 layout differences of the base code and the frozen tolerance** (item 4, second
part). The fixed-input accumulators of trpv1 (35,645 particles, c2, box 256 at 1.2156 A, mask
diameter 175 A) on the base build with `SIMPLE_SEED=20261103`: `reconstruct3D` from the committed
`ptcl3D` field of the solve3D project (`/data2/trpv1_c2/5_solve3D`), and `make_cavgs` (100
classes) from the `ptcl2D` field of the solve2D project (`/data2/trpv1_c2/2_solve2D`), each in
shared memory at six threads and at one thread and at `nparts` 1, 2 and 4 with six threads per
part (`make_cavgs` requires `nparts`; its one-part run executes in shared memory and needs
`refs=`, so its shared-memory layouts are the one-part runs at six and at one thread). Run
directories `scratch/base_runs/level4/`, the measurement in `measure.txt` there (the scratch
tool `scratch/keep/tools/mrcdiff.py` reads the maps for the largest absolute voxel or pixel
difference inside the spherical mask and the correlation of the masked maps; the FSC shells are
`simple_exec prg=fsc` on each pair of maps).

| pair of layouts | reconstruction: largest voxel difference in the mask (the maps' largest value is 10.10) | map correlation | FSC=0.5 and 0.143 shells between the pair | class averages: largest pixel difference over the 100 classes (relative to the class's largest value) | lowest class correlation |
| --- | --- | --- | --- | --- | --- |
| shared memory, six threads against one thread | 0 | 1 | 2.43 A (Nyquist), 2.43 A | 0 (the one-part runs at six and one thread) | 1 |
| shared memory against one part, six threads | 0 | 1 | 2.43, 2.43 | (the same run) | |
| one part against two parts | 6.25e-5 | 1.00000000 | 2.43, 2.43 | 7.79e-7 (1.17e-6) | 1.00000000 |
| one part against four parts | 8.58e-5 | 1.00000000 | 2.43, 2.43 | 6.33e-7 (9.02e-7) | 1.00000000 |
| two parts against four parts | 2.57e-5 | 1.00000000 | 2.43, 2.43 | 5.96e-7 (7.34e-7) | 1.00000000 |

The thread count changes nothing in either accumulator (the reductions are ordered within a
process); the part count changes the order of the partial sums and the maps by less than one part
in a hundred thousand of their largest value. The FSC between any two layouts' maps stays at the
Nyquist shell at both thresholds, so the shell is a pass/fail criterion rather than a measured
spread. **The tolerance of level 4, frozen here and never moved** (ruling 4): a reconstruction
from the same committed poses at two layouts of one build must agree to a largest absolute voxel
difference inside the mask of at most **1.72e-4** (twice 8.58e-5; maps on the scale of this
measurement, whose largest value is about 10), a correlation of the masked maps of at least
**0.99999998** (the observed correlations round to 1 at eight decimals; twice the printing
resolution), and FSC=0.5 and 0.143 shells between the two maps at the Nyquist shell; a class-average
set from the same `ptcl2D` field must agree to a largest absolute pixel difference of at most
**1.56e-6** per class (twice 7.79e-7) and a class correlation of at least 0.99999998. Phase 3
measures its layouts (shared memory; 2, 3 and 4 parts; one and eight threads) on the same inputs
with the same tool and holds them to these numbers.

**The frozen-input harness** (item 5; level 2 of section 8). `polar_matching_bench` gets the
suite `fills` (`simple_test_exec test=polar_matching_bench suite=fills nthr=8`, logs and the table
`polar_fill_identity.tsv` in the working directory). The bench already builds synthetic references
(directions times states, both halves) and CTF-weighted particles at production polar sizes in its
own calculator; with a route given, its constructor now builds a second copy of the same data in a
builder's calculator (the tables reach the calculator through `build%pftc`), a particle field with
the CTF records, the state (round robin, so every state is populated above the tables' threshold
of five), a previous pose of random direction, a previous class, the sampling bookkeeping
`calc_num2sample` reads (`sampled`, `dist`) and the centring status of the polar shift map (odd
particles not yet centred take the disk, even ones the window), the point group `c1` and its
direction space; the parameters come from the bench's own command line with the route's keys
(`refine=prob` or `prob_neigh`, `oritype`, `prob_neigh_mode=snhc` for the 3D neighbourhood,
`ncls` equal to the reference count for 2D, `trs=3`; the polar shift map on, as in production).
One run of the harness (`run_fill(nthr, batchsz, nparts, result)`) is the production two-phase
path in process: one worker table per part (the part's particles in order, as the distributed
workers hold them), each filled in batches of `batchsz` on a team of `nthr` threads and written to
its part file by the production writers; the parts merged into the global table by the production
reader; the insertion cone set (3D, 300 A, above the mask diameter, so that every refined entry
is an eligible extra on random rows); the assignment drawn by the production `ref_assign` of the
route. The seed is reset before the fills (so that a one-thread run reproduces under today's
`seed_rnd`) and before the assignment. The result is gathered by particle and full reference
identity: the rows' distances and in-plane indices (the dense routes' rows; the sparse routes'
candidate lists placed at their references), the refined entries' shifts, the particle
temperatures, the seeded shifts and their flags, the assignment (reference, in-plane index,
distance, shift) and the extras of the weighted insertion (reference, in-plane index, shift,
weight). `compare_fill_results` counts the entries differing between two results in those seven
categories (a distance differing in value or in presence counts once).

The suite runs each route at 200 directions, one and two states and 64 particles (unless
`nspace`, `nstates`, `nptcls` are given): the reference at one thread, one batch and one part,
then five configurations against it: the same again (the repeat), the full team (`nthr`), batches
of a third at one thread, three parts at one thread, and all three together. The run fails when
the repeat differs anywhere or the temperatures, which no draw enters, differ in any
configuration; every other count is reported. On the code of this phase the table below is what
the harness must show (section 13, item 5): the drawn entries differ wherever the thread count,
the batch size or the part count changes (the fills reseed at every call and draw on the calling
thread's stream, so a particle's draws depend on its batch, its part and its thread), the
deterministic parts agree everywhere, and the one-thread repeat is identical. Phase 1's keyed
draws turn every row of the table to zeros, and the suite's checks are strengthened to assert it.

The final run (`scratch/keep/phase0/harness/phase0_release_final3.tsv`, the tree build, Release,
`nthr=8`; two earlier runs with narrower cones, 8 and 60 A, in the same directory): the entries
differing from the one-thread, one-batch, one-part reference, by category. "Entries" are the
row entries compared (references times particles for the dense routes; the candidates for the
sparse ones, where a candidate present in one result and absent in the other counts in every
category it carries).

| route | states | threads | batch | parts | entries | distances | in-plane | refined shifts | temperatures | seeds | assignments (of 64) | extras |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3D_dense | 1 | 1 | 64 | 1 | 12800 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 3D_dense | 1 | 8 | 64 | 1 | 12800 | 7820 | 7818 | 256 | 0 | 0 | 20 | 12 |
| 3D_dense | 1 | 1 | 22 | 1 | 12800 | 5126 | 5124 | 173 | 0 | 0 | 16 | 11 |
| 3D_dense | 1 | 1 | 64 | 3 | 12800 | 5157 | 5156 | 180 | 0 | 0 | 15 | 9 |
| 3D_dense | 1 | 8 | 22 | 3 | 12800 | 7722 | 7722 | 265 | 0 | 0 | 23 | 12 |
| 3D_neigh | 1 | 1 | 64 | 1 | 704 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 3D_neigh | 1 | 8 | 64 | 1 | 704 | 1222 | 1222 | 1214 | 0 | 0 | 61 | 0 |
| 3D_neigh | 1 | 1 | 22 | 1 | 704 | 802 | 802 | 800 | 0 | 0 | 41 | 0 |
| 3D_neigh | 1 | 1 | 64 | 3 | 704 | 811 | 811 | 807 | 0 | 0 | 40 | 0 |
| 3D_neigh | 1 | 8 | 22 | 3 | 704 | 1215 | 1215 | 1206 | 0 | 0 | 60 | 0 |
| 2D_dense | 1 | 1 | 64 | 1 | 12800 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 2D_dense | 1 | 8 | 64 | 1 | 12800 | 7733 | 7731 | 268 | 0 | 0 | 22 | 0 |
| 2D_dense | 1 | 1 | 22 | 1 | 12800 | 5126 | 5124 | 173 | 0 | 0 | 16 | 0 |
| 2D_dense | 1 | 1 | 64 | 3 | 12800 | 5157 | 5156 | 180 | 0 | 0 | 15 | 0 |
| 2D_dense | 1 | 8 | 22 | 3 | 12800 | 7750 | 7749 | 255 | 0 | 0 | 24 | 0 |
| 2D_sparse | 1 | 1 | 64 | 1 | 9600 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 2D_sparse | 1 | 8 | 64 | 1 | 9600 | 7901 | 7901 | 3764 | 0 | 0 | 52 | 0 |
| 2D_sparse | 1 | 1 | 22 | 1 | 9600 | 5300 | 5300 | 2495 | 0 | 0 | 31 | 0 |
| 2D_sparse | 1 | 1 | 64 | 3 | 9600 | 5264 | 5264 | 2441 | 0 | 0 | 31 | 0 |
| 2D_sparse | 1 | 8 | 22 | 3 | 9600 | 7856 | 7855 | 3719 | 0 | 0 | 51 | 0 |
| 3D_dense | 2 | 1 | 64 | 1 | 25600 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 3D_dense | 2 | 8 | 64 | 1 | 25600 | 15597 | 15597 | 265 | 0 | 0 | 23 | 11 |
| 3D_dense | 2 | 1 | 22 | 1 | 25600 | 10363 | 10361 | 169 | 0 | 0 | 11 | 7 |
| 3D_dense | 2 | 1 | 64 | 3 | 25600 | 10387 | 10385 | 175 | 0 | 0 | 13 | 5 |
| 3D_dense | 2 | 8 | 22 | 3 | 25600 | 15487 | 15486 | 248 | 0 | 0 | 23 | 8 |
| 3D_neigh | 2 | 1 | 64 | 1 | 1344 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 3D_neigh | 2 | 8 | 64 | 1 | 1344 | 2455 | 2455 | 2446 | 0 | 0 | 62 | 0 |
| 3D_neigh | 2 | 1 | 22 | 1 | 1344 | 1597 | 1597 | 1595 | 0 | 0 | 42 | 0 |
| 3D_neigh | 2 | 1 | 64 | 3 | 1344 | 1609 | 1609 | 1604 | 0 | 0 | 42 | 0 |
| 3D_neigh | 2 | 8 | 22 | 3 | 1344 | 2420 | 2420 | 2414 | 0 | 0 | 61 | 0 |
| 2D_dense | 2 | 1 | 64 | 1 | 25600 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 2D_dense | 2 | 8 | 64 | 1 | 25600 | 15603 | 15601 | 236 | 0 | 0 | 30 | 0 |
| 2D_dense | 2 | 1 | 22 | 1 | 25600 | 10362 | 10361 | 148 | 0 | 0 | 17 | 0 |
| 2D_dense | 2 | 1 | 64 | 3 | 25600 | 10394 | 10392 | 157 | 0 | 0 | 18 | 0 |
| 2D_dense | 2 | 8 | 22 | 3 | 25600 | 15662 | 15662 | 230 | 0 | 0 | 23 | 0 |
| 2D_sparse | 2 | 1 | 64 | 1 | 19200 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 2D_sparse | 2 | 8 | 64 | 1 | 19200 | 15805 | 15801 | 7533 | 0 | 0 | 55 | 0 |
| 2D_sparse | 2 | 1 | 22 | 1 | 19200 | 10496 | 10496 | 4963 | 0 | 0 | 35 | 0 |
| 2D_sparse | 2 | 1 | 64 | 3 | 19200 | 10505 | 10503 | 4963 | 0 | 0 | 37 | 0 |
| 2D_sparse | 2 | 8 | 22 | 3 | 19200 | 15727 | 15727 | 7478 | 0 | 0 | 54 | 0 |

Read: the one-thread repeat of every route and state count is identical in every category (the
policies' promise that a one-thread run reproduces holds for the fills and the assignment);
every change of thread count, batch size or part count changes most drawn entries (the in-plane
draw of every scored reference, hence the distances, the refined shifts that start from the drawn
rotation, the sparse routes' candidate sets, and so the assignments; the 3D neighbourhood's
candidate lists differ almost entirely, as their picks come from the thread's tabu sampler); the
temperatures and the seeded shifts never differ (the temperature is the particle's own noise
scale, the seed under the polar shift map is zero with its flag set). The extras of the weighted
insertion are exercised on the 3D dense route (about 1.2 poses per particle at the 300 A cone)
and differ with the rows they are drawn from; the 3D neighbourhood's `shc`/`snhc` modes insert
the committed pose alone by design (`ref_assign_neigh`, line 788), and 2D has no extras, so those
columns are empty. Between two processes (the runs with the 8 A and the 60 A cone, same seed),
every one-thread row of the table is identical and every eight-thread row differs, the thread
dependence at the process level. Phase 1's exit is this table at zero in every column.

**Exit criteria of phase 0.** Both builds pass the fast gate: the base build 15 of 15 suites
(`scratch/keep/phase0/logs/build_base_gate.log`), the tree with the phase's changes 15 of 15
(`gate_tree_phase0.log`, repeated as `scratch/logs/gate_tree_phase0_final.log` after the last
edit; the Debug tree passes it too, `gate_debug_phase0_final.log`). The audit (above), the
vectors (above, Release and Debug), the tolerance (above) and the baselines (above and below) are
in this entry. The harness runs on the four routes (above).

**The real-data readouts and the bench counters of section 2.5 on this machine** (item 4, third
part; base build, `SIMPLE_SEED` set, the machine otherwise quiet with the load recorded in
`scratch/base_runs/readouts/timeline.txt`; the extraction in `extract.txt` there, by
`scratch/extract_readouts.sh`). These are the numbers phase 1 judges its trajectory change against
(section 12, level 4: within the seeds' spread) and phase 4 measures the final build against.

*bgal `solve2D`* (`/data2/bgal_addon/3_selection`, 2,756 active of 5,513 particles, 40 classes,
`nsample=1000`, four parts of six threads; 26-27 s per run): the class resolutions of the final
class set over the populated classes.

| seed | populated classes | median resolution (A) | mean of the ten best (A), and the ten | classes under 10 A |
| --- | --- | --- | --- | --- |
| 11 | 37 | 9.33 | 6.94 (4.73 6.40 6.40 6.66 7.25 7.25 7.59 7.59 7.77 7.77) | 23 |
| 22 | 37 | 9.07 | 6.92 (5.26 6.40 6.66 6.80 6.80 7.10 7.25 7.42 7.77 7.77) | 21 |
| 33 | 37 | 9.60 | 6.53 (4.87 5.18 5.53 6.40 6.80 6.80 7.25 7.42 7.42 7.59) | 21 |
| 20261103 | 39 | 11.66 | 6.87 (4.24 6.40 6.40 6.40 6.94 7.42 7.59 7.77 7.77 7.77) | 19 |

The spread over the four seeds: median 9.07 to 11.66 A, the ten best 6.53 to 6.94 A, under 10 A
19 to 23 classes. A run at six threads is a realisation, not a reproduction: the same seed on
2026-11-03 gave seed 11 a median of 9.60 A, ten best 6.54 A and 21 classes under 10 A, which is
the thread dependence of today's draws (the harness above) seen on real data, and the reason the
readouts are judged within the seeds' spread.

*trpv1 `refine3D`* (`/data2/trpv1_c2/5_solve3D`, 35,645 particles, c2, from its map; `refine=prob
objfun=euclid box_crop=160 shift_sigma=0 maxits=3`, four parts of six threads, seed 20261103;
1,000 s in all): the FSC shells of the assembled map after each iteration.

| iteration | FSC=0.143 (A) | FSC=0.5 (A) |
| --- | --- | --- |
| 1 | 4.041 | 6.102 |
| 2 | 4.095 | 6.621 |
| 3 | 4.095 | 6.765 |

They are the shells the same command gave on 2026-11-03 to the last digit printed (its run
`speed_items_2026-11-03/runs/item2/trpv1_r3d`): the resolution readout of this case is robust to
the draws' thread dependence, so phase 1 holds its first iteration to 4.04 and 6.10 A and its
third to about 4.1 and 6.8 A, within the spread the seeds give (to be measured on the new code
with the four seeds if the first comparison is not within a few hundredths of an angstrom).

*The counters and timings of section 2.5 on the Dell* (the same runs; the 2D cases at 100
classes, two iterations, four parts of six threads: `trpv1_r2d_box256` from the solve2D class
averages at box 256, `trpv1_r2d_crop120` at `box_crop=120` from the scaled class averages of
2026-11-03):

| | 3D, crop 160, 2,500 directions | 2D, crop 120, 100 classes | 2D, box 256, 100 classes |
| --- | --- | --- | --- |
| table job per part | 240-306 s (fill 236-302 s, preparation 2.3-2.9 s, read 0.4-0.7 s, references 0.5-0.6 s) | 10-12 s (fill 7.6-9.7 s, preparation 1.9-2.2 s, read 0.4-0.5 s, references 0.1 s) | 24-28 s (fill 18-21 s, preparation 5.4-5.6 s, read 0.4-0.6 s, references 0.3 s) |
| master assignment | 1.2 s (reads 0.5 s, `ref_assign` 0.75 s with the extras of `topk=3`) | 0.05 s (reads 0.02 s, `ref_assign` 0.03 s) | 0.05 s (reads 0.02 s, `ref_assign` 0.03 s) |
| search job per part | 20.5-22.8 s: preparation 2.7-3.2, references 0.4, read 0.4-0.6, polish 7.8-10.5, partial reconstruction 8.3-9.1 | 4.7-5.2 s: preparation 1.7-1.9, read 0.4-0.5, polish 0.6-0.7, class sums 1.6-1.9 | 14.9-16.0 s: preparation 5.0, read 0.6, polish 1.1-1.3, class sums 7.7-8.6 |
| iteration | 319-353 s (the probabilistic pre-step 292-326 s, the search phase 22-23 s, assembly 2.1 s) | 18 s (the pre-step 12.5 s, the search phase 5.3 s, assembly 0.15 s) | 45 s (the pre-step 27.8 s, the search phase 16.2 s, assembly 0.5 s) |
| master table and part files | 35,504 sampled particles x 2,500 references x 8 bytes = 710 MB; part files 198 MB each | | |

The Dell's numbers sit where section 2.5's do (its 3D iteration 309-362 s, its 2D ones 30 and
78 s; the 2D runs here are the first iteration of two): the fill is 80-85 % of a 3D iteration;
the preparation, read and reference load duplicated by the two phases are about 1 % of a 3D
iteration, about 12 % of the cropped 2D one (2.2 s of 18) and about 13 % of the 2D one at box
256 (6 s of 45).


### Phase 1: keyed draws (2026-11-04, the Dell)

Done on 2026-11-04 by the autonomous run `prob_pipeline` on the Dell, on the tree phase 0 left (base HEAD
`22a643642`; the run root `/data2/agent_runs/prob_pipeline`, linked as `~/agent_runs/prob_pipeline`, and paths
below that are not absolute are relative to it). The tree built in `repo/build` (Release) and
`scratch/keep/build_debug` (Debug); the base build of phase 0 (`scratch/keep/build_base`) is the comparison
point of the readouts. The journal is `state/journal_phase_1.md`.

**What changed, in one paragraph.** Every random draw of the probabilistic routes now comes from the keyed
generator of phase 0: the sample the project owner draws, the in-plane rotation drawn for every scored
reference, the sparse fills' candidate order, the assignment draw of every particle, the 2D search's fallback
class pick and the shift map's disk draw (which production does not reach). A draw's key is the run's seed
(`SIMPLE_SEED`, drawn and exported by the project owner when unset, recorded in the project with the
generator's version), the program and pass, the iteration, the particle's row in the field, the family, the
reference and the draw's ordinal: the same logical draw returns the same bits on any thread, in any batch, in
any part, in any process. The fills no longer reseed at every batch and no longer read the field for the
iteration's statistics: the refined-reference budget and the shift prior's width are computed once per
iteration by the fill's caller and passed in. The harness of phase 0 now gives identical rows, selections and
extras in every configuration, including one that perturbs the earlier parts' committed outputs before the
later parts' fills. The trajectory of every probabilistic run changed once, here, as section 4.5 says; the
classic modes, the Cartesian polish pass, the reconstruction and the 2D class-sum order are untouched.

**Files changed** (the file table's rows of phase 1, plus the rows it added): `src/utils/math/simple_rnd.f90` (the run seed, the scope and draw keys, the cached seed), `simple_rnd_tester.f90`; `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` (the selector takes its variate; the scope; the keyed candidate order), `simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90`, `simple_eul_prob_tab2D.f90` (the keyed fills and assignments, the iteration model, the reseeds and tabu samplers removed), `simple_eul_prob_tab_tester.f90`, `simple_eul_prob_tab2D_tester.f90`; `src/main/pftc/simple_polarft_calc.f90`, `simple_polarft_core.f90`, `simple_polarft_corr.f90`, `simple_polarft_shiftmap.f90` (the calculator's scope and the keyed in-plane and disk draws), `simple_polarft_corr_tester.f90`, `simple_pftc_inplane_tester.f90`, `simple_polar_shiftmap_tester.f90`; `src/main/ori/simple_oris.f90`, `simple_oris_sampling.f90`, `simple_oris_tester.f90` (the keyed samplers); `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` (the sample's key, the project owner's seed), `simple_strategy2D_srch.f90` (the fallback pick); `src/main/strategies/parallelization/simple_refine3D_strategy.f90`, `simple_refine2D_strategy.f90` (the seed adopted at initialisation); `src/main/project/simple_sp_project.f90`, `simple_sp_project_core.f90`, `simple_sp_project_tester.f90` (the restart fields); `src/utils/qsys/simple_qsys_ctrl.f90` (the seed exported into the queue scripts); `src/main/commanders/simple/simple_commanders_prob.f90` (the workers compute the iteration model); `src/main/commanders/test/simple_polar_matching_bench.f90`, `simple_commanders_test_highlevel_matchbench.f90`, `src/main/ui/simple_test/simple_test_ui_highlevel.f90` (the harness); and this plan (section 4.3's rows for the sample's lanes, the previous-last rule and the fallback pick; the file-table rows this phase added; this entry). Thirty-one files, about 2,300 lines added and 300 removed over phases 0 and 1. The rows of `simple_parameters_phases.f90`, `simple_commanders_refine3D.f90` and `simple_commanders_refine2D.f90` were not needed (the project owner is the strategy, where the project is read; nothing new is parsed).

**The audit's closure** (phase 0's table, site by site). Every site of the six families now draws on the keyed
generator through one scope per iteration, `rnd_scope_key(seed, program, pass, iteration)`: the seed is the
run's `SIMPLE_SEED`, the program `refine3D` or `refine2D`, the pass the discrete alignment, the iteration
`which_iter`. The tables take the scope at every fill and assignment (`prob_rnd_scope` in the utilities module,
from `params`) and set it on the calculator they fill through (`set_rnd_scope`), so no caller outside the table
modules threads a key through; the sampling helpers of the matcher build the sample's key from `params` under
the probabilistic modes.

| family | site | keyed on | what went |
| --- | --- | --- | --- |
| `SAMPLE` | `sample4update_cnt`, `_rnd`, `_class`, `_fillin` (`simple_oris_sampling.f90`), called from `sample_ptcls4update3D`, `sample_ptcls4fillin` and `sample_ptcls4update2D` (`simple_matcher_smpl_and_lplims.f90`) | an optional key: particle 0, family `SAMPLE`, the sampler's shuffle count in the reference lane, the swap in the ordinal lane (`rnd_keyed_partial_shuffle`); the matcher's helpers pass it under `refine=prob` and `prob_neigh` (`sample_rnd_key`), the classic modes pass none and draw on the intrinsic stream as before | nothing: the classic samplers are untouched |
| `INPL` | `gen_posterior_val` and `gen_posterior_val_bank` (`simple_polarft_corr.f90`), through `draw_posterior_crvec` and `draw_posterior_row` | the calculator's scope, the particle (the project index the fills pass), the reference (the calculator's index); the same pair in the same pass draws the same variate on any thread, in any batch or part | `draw_posterior_row` takes the variate as an argument instead of drawing it; `gen_posterior_evid` (no caller) deleted; the calculator stops when a draw is asked without a scope |
| `DISK` | `shiftmap_draw` with `l_greedy=.false.` (`simple_polarft_shiftmap.f90`) | the scope, the pair the map holds (`wk%iptcl`, `wk%iref`), an optional draw ordinal (the screen's; 0 by default) | the intrinsic draw; still reached only by the shift-map tester at this HEAD, which passes the draw's ordinal |
| `CAND` | the 3D neighbourhood fill's candidate order (`evaluate_direct_stochastic_refs`) and the 2D sparse fill's (`process_particle` of `fill_tab_prob_neigh_range`) | `keyed_candidate_order(scope, particle, n, previous, order)`: the previous reference or class last (as `put_last` arranged it), the others a keyed partial Fisher-Yates permutation of the particle (family `CAND`, one keyed uniform per swap); the candidate counts, the stopping rules and the forced incumbent of the 3D stochastic mode unchanged | the per-thread tabu samplers (`direct_rts`, `ran_tabu`), `ne_ran_iarr` and `put_last` from the fills |
| `ASSIGN` | the 2D search's fallback class pick (`prep4srch`, `simple_strategy2D_srch.f90`), a populated class drawn until one is found when the previous class is empty or the run starts fresh | the scope, the particle, family `ASSIGN`, the pick's ordinal from 1 (the assignment draw is ordinal 0; a first draft keyed it under `CAND`, where it collided with the sparse fill's candidate permutation of the same particle, which the review caught), under the probabilistic modes; the intrinsic stream under the classic ones | nothing else |
| `ASSIGN` | `ref_assign` (3D dense), `ref_assign_neigh`, the 2D `ref_assign` and `ref_assign_sparse` | the scope, the particle (`pinds(i)`), family `ASSIGN`, ordinal 0; the neighbourhood assignment's fallback in-plane draw (`seed_fallback_if_empty`) goes through the calculator's `INPL` key | the master's stream |
| `DIAG` | none at this HEAD (the readouts and the posterior diagnostic draw nothing) | reserved | |

The fills' `seed_rnd` calls (`simple_eul_prob_tab.f90` 302, `_neigh` 197 and 451, `_2D` 282 and 528) are gone; the
reseed counter of `simple_rnd.f90` stays for `parameters%new` and the classic path. The three classic samplers
of the utilities module (`angle_sampling_2`, `sample_bounded_dists`, `sample_power_dists`) and the 2D matcher's
per-batch search order (`prep_strategy2D_batch`), unread on the probabilistic routes, keep the intrinsic
stream. Of the dead code the audit found, `select_particles_set` (`simple_oris.f90`, a file of this phase) is
deleted; `rnd_inpl`, `rnd_shift`, `rnd_proj_space` and `split_state` are left for the maintainer: their files
are not in this run's table and the rule applies to what the run touches.

**The run's seed and the restart metadata** (section 4.2 and 4.4). The project owner is the refine strategy
on its own parameters and builder, where the project is read: `adopt_run_seed(params, build)`
(`simple_matcher_smpl_and_lplims.f90`) runs in the shared-memory and the distributed initialisation of
`refine3D` and `refine2D` under the probabilistic modes (`refine=prob`, `prob_neigh`; the classic modes are
left as they are, ruling 6), before any part is scheduled and before any draw of the run. A run keeps a set
`SIMPLE_SEED`; a run without one draws a 64-bit seed from the operating system's generator (the clock and the
process id when that fails), non-negative so that it reads back as written, and exports it with
`put_simple_seed` (`ensure_run_seed`, `simple_rnd.f90`), so every part, every nested commander and every keyed
draw of the run reads the same seed; the queue scripts of every backend now carry `export SIMPLE_SEED=<seed>`
whenever the seed is set (`simple_qsys_ctrl.f90`), because a part run by the stream's persistent worker inherits
that worker's environment, not the master's; the line stands after the scheduler's directive block, since SLURM,
PBS and LSF stop reading their directives at the first non-comment line (the review's second round caught it
placed before them, where a cluster run would have lost its partition, time and memory directives; the local
backend writes none, so no test showed it). The seed is cached by the serial code that sets or reads it
(`seed_rnd` at every commander's start, `put_simple_seed`, `ensure_run_seed`, `set_run_seed`), so a keyed draw
inside a parallel region reads no environment and takes no lock; an uncached read inside a region is a
programming error that stops the run. The
seed and the generator's schema version are recorded in the project's info segment (`set_rnd_seed`,
`get_rnd_seed` of `sp_project`: the 64-bit seed as four 16-bit words `rnd_seed_w1`..`w4`, numbers the project's
text records keep exactly, and `rnd_version`; a text field came back as a number from the record parser, which
the review caught and the project tester's round-trip test now guards with a 63-bit seed with all four words set
and a negative seed) and the segment written once at
initialisation; a `solve3D` or `solve2D` stage finds the seed set by its first probabilistic stage and records
the same values again. A run that continues (`continue=yes`, the distributed `refine3D`) reuses the recorded
seed (`set_run_seed`): a `SIMPLE_SEED` set to another value, or a recorded generator version other than this
build's, is inconsistent restart metadata and stops the run with a message naming both values. The scope (the
program and the discrete alignment pass) and the next iteration are the `prg` and `startit` of the job
description the run already carries, so nothing else is written; the generator's version lives in every key.
No parameter was added (ruling 2): the planned `rnd_seed` key of
`doc/refactoring_notes/planned/rnd_seed_parameter.md` stays a plan for the maintainer, and `parameters` reads
nothing new (the table row for `simple_parameters_phases.f90` was not needed: the project owner validates the
restart metadata where it reads the project).

One side effect to know: `has_global_seed()`, true whenever `SIMPLE_SEED` is set, makes the polar memo, the
shift map, the PCG reconstructor and the wave projector plan their FFTs with `FFTW_ESTIMATE` instead of
`FFTW_MEASURE`, so that a seeded run replays exactly. Since the project owner now exports a seed into the
environment of every probabilistic run, those runs and their parts plan that way, including runs the user did
not seed. The plans are then reproducible; whether the estimated plans cost time is for phase 4's wall-time
measurements against the base build, which runs unseeded production the other way.

**Level 1 (draw identity)** holds unchanged: the generator's constants and vectors are phase 0's (nothing in
`rnd_keyed` moved; the vectors of phase 0's table are reproduced by the Release and the Debug build of this
phase, `unit_numerics`, suite "random draws", 130 checks with the scope tests added: the run seed read
back, the scope key's fields, a draw key completing it, the family and the program separating one particle's
draws). The selector `draw_posterior_row` takes its unit variate and is tested against a cumulative-posterior
oracle of its own arithmetic (single-precision weights summed in double): a variate below an entry's cumulative
bound draws that entry, at the bound the entry itself, 0 draws the first valid entry, the largest variate the
last, a huge temperature the argmin at every variate and either of two equal minima by the variate's half,
temperature 0 the valid entries uniformly over 20,000 keyed variates, the row [0, 1, 2] at temperature 1 with
masses 0.665, 0.245, 0.090 within 0.02, the same variate the same draw, and a row without a valid entry 0
("3D probability table", `unit_pftc_align2D3D`).

**Level 2 (frozen-input identity)** holds: the harness's first run on this build
(`scratch/keep/phase1/harness/phase1_run1.tsv`, before the suite's checks were strengthened) gave every one of
its 40 rows zero in every category (distances, in-plane indices, refined shifts, temperatures, seeds,
assignments, extras) at eight threads, in batches of a third, in three parts and all three together, on the
four routes at one and two states: the entries phase 0's table showed differing by the thousands. The suite
now fails on any differing entry in any configuration and runs a sixth configuration, three parts at one
thread with the earlier parts' committed outputs perturbed (their distances, states, classes and shift
increments changed as a commit would) before the later parts' fills, the input-boundary regression of section
3.2: the second run (`phase1_run2.tsv`), 48 rows: every count zero in every configuration, the perturbed ones included. The table, with the six configurations of each route and state count collapsed (every entry of every row is 0):

| route | states | entries per configuration (dense: references x particles; sparse: the candidates) | configurations | largest count in any category |
| --- | --- | --- | --- | --- |
| 3D_dense | 1 | 12,800 | 6 (1 and 8 threads; batches of 64 and 22; 1 and 3 parts; 3 parts with the earlier parts' outputs perturbed) | 0 |
| 3D_neigh | 1 | 704 | 6 (1 and 8 threads; batches of 64 and 22; 1 and 3 parts; 3 parts with the earlier parts' outputs perturbed) | 0 |
| 2D_dense | 1 | 12,800 | 6 (1 and 8 threads; batches of 64 and 22; 1 and 3 parts; 3 parts with the earlier parts' outputs perturbed) | 0 |
| 2D_sparse | 1 | 9,600 | 6 (1 and 8 threads; batches of 64 and 22; 1 and 3 parts; 3 parts with the earlier parts' outputs perturbed) | 0 |
| 3D_dense | 2 | 25,600 | 6 (1 and 8 threads; batches of 64 and 22; 1 and 3 parts; 3 parts with the earlier parts' outputs perturbed) | 0 |
| 3D_neigh | 2 | 1,344 | 6 (1 and 8 threads; batches of 64 and 22; 1 and 3 parts; 3 parts with the earlier parts' outputs perturbed) | 0 |
| 2D_dense | 2 | 25,600 | 6 (1 and 8 threads; batches of 64 and 22; 1 and 3 parts; 3 parts with the earlier parts' outputs perturbed) | 0 |
| 2D_sparse | 2 | 19,200 | 6 (1 and 8 threads; batches of 64 and 22; 1 and 3 parts; 3 parts with the earlier parts' outputs perturbed) | 0 |

The second run's table printed the entries column wrongly (a formatting defect of the suite, fixed; the difference counts are what the check reads); the third run (`phase1_run3.tsv`, on the tree after the review fixes, eight threads, 384 s beside the level-3 runs) printed it right and gave the same 48 rows with every count zero, the entries per configuration as tabled. The fills' times in that run (one thread, batches of 64, one part): 3D dense 6.0 s at one state and 10.4 s at two, 3D neigh 1.9 and 4.0 s, 2D dense 6.2 and 11.9 s, 2D sparse 9.3 and 21.2 s; at eight threads 0.77 and 1.55 s, 0.23 and 0.54 s, 0.78 and 1.66 s, 1.13 and 3.65 s, the same in the three runs to a few percent and inside the spread of phase 0's three runs (3D dense at one state 4.8 to 6.1 s, 2D sparse 7.0 to 9.3 s, the eight-thread rows 0.82 and 1.24 s): the keyed draws cost the fills nothing the harness can see (the fill's cost is measured against the base in phase 4).

The OpenMP rule (new or changed parallel code tested at sixteen threads, repeated, against one thread): the
suite run twice at sixteen threads on the final tree (`phase1_nthr16_rep1.tsv`, `rep2.tsv`; 335 s each), 48 rows
each (32 at one thread, 16 at sixteen), every difference count zero in both, and the two tables agree in every
count. The fills at sixteen threads take 0.46 s (3D dense, one state, batches of 64), 0.13 s (3D neigh), 0.45 s
(2D dense) and 0.68 s (2D sparse), the same to a few percent in the repetition.


**The level-3 reference directories of this build** (item 4; the reference of phases 2-4). Under
`scratch/keep/identity/phase1/` of the run root, made by `run_phase1_level3.sh` as phase 0's baselines were
(`SIMPLE_SEED=20261103`, the entries as `simple_test_exec` runs with `nthr` on the command line and
`OMP_NUM_THREADS` set to it, the distributed cases as `simple_exec` runs on copies of the 2026-11-03 base run's
projects), with two shared-memory probabilistic cases at eight threads that the entries do not cover (phase 0's
finding that the simulated workflows pin one thread): the bgal `solve2D` of the readouts in shared memory and a
trpv1 `refine3D` of two iterations (`refine=prob objfun=euclid box_crop=160 shift_sigma=0`). The timeline is in
the directory; the machine carried the truth suites and the readouts at the same time, so the times are not
measurements.

| run | layout | result | time |
| --- | --- | --- | --- |
| `simulated_workflow_1jxy` | one thread; "eight threads" | both stop at solve3D stage 3, iteration 41, with `the map's last shell must lie in the band; new_shiftmap` (the base failure), exit 1, no `FAIL:` line before it | 479 s; 498 s |
| `simulated_workflow_6vxx` | one thread; "eight threads" | both complete with the one known `FAIL: final-volume correlation is below the required minimum` (the base failure), exit 0 | 1,689 s; 1,774 s |
| `mini_stream_1jxy` | eight threads | passes (`PASS: mini_stream 1jxy validated 6 micrographs, 21 particles, and 2 classes`) | 68 s |
| `solve2D` of the 1jxy workflow's extracted project, `nparts=2 nthr=1 ncls=4` | two parts | normal stop | 19 s |
| the same with `nsample=30` (the fractional variant) | two parts | normal stop | 18 s |
| `solve3D` of the 6vxx workflow's solve2D project, `pgrp=c3 nparts=2 nthr=1` | two parts | normal stop | 1,295 s |
| bgal `solve2D` in shared memory, `ncls=40 nsample=1000 nthr=8` | eight threads | normal stop | 33 s |
| trpv1 `refine3D` in shared memory, two iterations, `nthr=8` | eight threads | normal stop; FSC=0.143 at 4.041 then 4.095 A, FSC=0.5 at 6.224 then 6.621 A, the readout's values | 1,331 s |

The base failures reproduce at the same points as on the base build: `simulated_workflow_1jxy` stops at solve3D
stage 3, iteration 41, with `the map's last shell must lie in the band; new_shiftmap` at both thread counts (its
log the same length as the base run's, 2,578 lines); `simulated_workflow_6vxx` completes with its one
`FAIL: final-volume correlation is below the required minimum` (whole-volume correlation 0.5415 against the base
run's 0.5492, floor 0.80; its docking correlation passes at 0.870, the base's 0.874). `mini_stream_1jxy` passes
as before. The trajectories differ from the base build's, as section 4.5 says they must: compared with
`compare_dirs.py`, the base `solve2D_distr2` and this build's agree in 27 files and differ in 117 (every class
average, the project, the statistics).

The reference directories were made twice: once by the build of 11:16 (the installed binaries, used by the
truth suites and the readouts) and again by the final tree (a staged install) after the review's second round,
whose fixes change no output; the first set was kept beside the second under a name of the same length
(`ph1old`, so that binary offsets of embedded paths agree) and the two compared with `compare_dirs.py`
(`scratch/keep/tools/`; byte identity per file, logs and markers skipped; results in
`scratch/phase1_runs/level3_compare/`). Every map, stack, star and binary output agrees byte for byte:
`simulated_workflow_1jxy` 130 mrc, 21 mrcs, 37 star and all 8 projects at both thread counts;
`simulated_workflow_6vxx` 175 mrc, 30 mrcs, 43 star and 7 of 8 projects (the solve3D stage's project carries
the run id); `mini_stream_1jxy` 107 mrc, 7 mrcs and 27 star (compared below its pid-named test directory);
`solve2D_distr2` 82 mrc and 27 star, its fractional variant 85 and 28; `solve3D_distr2` 41 mrc, 4 bin and 2
dat (its 104 differing files are the FSC plots' PDFs, which carry their creation time); the shared-memory trpv1
`refine3D` at eight threads 9 mrc. What differs is what carries a time, a run id or a process id: the job
records' `time` field in a project, `metadata.json`'s heartbeat, the PDFs, the pid-named sample images. The
registry's `test=identity` on the entries agrees for every declared map and stack but reports the projects as
`DIFF`: it strips each run directory's own prefix from the records, and the first set's records name the
directory it was moved out of; `cmp` shows those project files byte-identical.

The one case that does not agree is the shared-memory bgal `solve2D` at eight threads, whose runs diverge from
the second class averages on. It is not the build: two runs of the final tree diverge the same way (iteration
1's mean in-plane distance 88.9 against 91.1 deg; `scratch/phase1_runs/bgal_shmem_repro/`), and the log shows
why: `solve2D`'s first stage runs in refinement mode `snhc_smpl`, a classic mode whose draws stay on the
intrinsic per-thread stream (unchanged by design; the probabilistic stages `prob` and `prob_neigh` follow it),
so at eight threads the stage realises a different trajectory each run, as the base did. At one thread the same command reproduces itself (two runs: all 94 class-average stacks and every star file
byte-identical, only the project's time field and the heartbeat differ), and differs from the eight-thread run
the same way the two eight-thread runs differ from each other. The
distributed `solve2D` cases reproduce because each part runs one thread. The shared-memory 2D reference at
eight threads is therefore a purely probabilistic case instead: the trpv1 `refine2D` of the readouts (`refine=prob objfun=euclid`, crop 120, 100 classes, two iterations) in
shared memory at eight threads, run twice (`scratch/phase1_runs/r2d_prob_shmem/`, 195 and 282 s under load):
all six class-average stacks, both `clusters2D` star files, the three binary states and the assignment file
byte-identical; the first run is kept as `trpv1_r2d_prob_shmem_nthr8` in the reference set (its line in the
timeline). The one file of that pair that differs is the table's part file, in its tail only: the refined
entries (particle, class, shift, distance, deviation; 354,488 of them, 24 bytes each) are appended from the
threads' candidate buffers in thread order, so two runs hold the same entries (the same multiset, checked) in a
different order, while the dense rows before them agree byte for byte. The readers key those entries by
particle and class, which is why nothing downstream differs; phase 2, which rewrites the part files (section
3.3), should write that section in particle order so the file itself is reproducible.

**The truth suites** (level 4 of section 8, ruling R3; `scratch/phase1_runs/truth/`, `run_phase1_truth.sh`,
one at a time at eight threads with `SIMPLE_SEED=20261103`, beside the level-3 runs). `multistate_workflows`
passes with its floors unchanged: exit 0 in 2,345 s, no `FAIL:` line, 126 `PASS:` lines, every metric of its
four asserted routes at or inside its floor (`states` 32 of 32, `cavgs` 32 of 32, `shc` 32 of 32, `consensus`
33 of 33; the `labels` route does not exist at this HEAD and is not counted). The gated values on the keyed
draws, with the floors:

| route | metric | value | floor |
| --- | --- | --- | --- |
| `states` | adjusted Rand index; minimum state recall; accuracy | 1.000; 1.000; 1.000 | at least 0.99 each |
| `states` | largest population error; discrimination margin | 0.000; 0.450 | at most 0.01; at least 0.40 |
| `states` | median rotation error (deg); median shift error (px) | 1.24; 0.047 | at most 2.0; 0.10 |
| `cavgs` | median rotation error (deg); median shift error (px) | 2.39; 0.237 | at most 4.5; 0.40 |
| `shc` | median rotation error (deg); median shift error (px) | 1.63; 0.059 | at most 2.5; 0.15 |
| `consensus` | adjusted Rand index; accuracy; largest population error | 0.666; 0.879; 0.121 | at least 0.40; at least 0.70; at most 0.30 |
| `consensus` | median rotation error (deg); median shift error (px) | 1.66; 0.063 | at most 2.5; 0.15 |

The suite's `success` markers are reported, not gated: `cavgs` resolved every true state (marker 1.0, its
shuffled control 0.0), `shc` did not (0.0, its control 0.0), as the suite prints them (`SUCCEEDED`, `MISSED`).

`solve3D_addon` passes with its floors unchanged: exit 0 in 178 s, no `FAIL:` line, 47 `PASS:` lines, 46 of
46 metrics inside their floors; the gated ones: the union map's correlation with the truth 0.963 (at least
0.90), its FSC=0.143 shell 4.56 A (at most 5.65 A, the base run's shell), the union minus the base correlation
+0.003 (at least -0.03), the cohort's coverage 1.00 (at least 0.90), the cohort's pose-error excess -0.10 deg
(at most 10 deg), the add-on report's union-base correlation 0.9998 (at least 0.90).

`cont_refine3D_1jxy` passes with its floors unchanged: exit 0 in 396 s, no `FAIL:` line, 138 `PASS:` lines,
135 of 135 metrics inside their floors (its gates: median rotation errors 0.23 to 0.92 deg against floors of
1.5 deg, median shift errors 0.05 to 0.07 px against 0.2 px, fractions in the 8 A basin 0.997 to 1.000 against
0.95, the in-place continuation identical to its fresh counterpart to 0.01 deg, the symmetric pose errors of
the polish, polar and continued routes 0.20 to 0.63 deg inside their floors). The three suites ran on the build
of 11:16 (the installed binaries); the fixes after it (the queue scripts' export line moved below the directives,
the cached seed read, a tester, a comment) change no output, which the level-3 comparison below confirms.

**The real-data readouts against the base** (item 4; `scratch/phase1_runs/readouts/`, `run_phase1_readouts.sh`,
the same commands, layouts and seeds as phase 0's baselines, four parts of six threads; the extraction in
`extract.txt` there by `scratch/extract_readouts.sh`; the machine carried the level-3 re-run at the same time,
so the times are not measurements and section 2.5's counters are left to phase 4).

*bgal `solve2D`* (2,756 active of 5,513 particles, 40 classes, `nsample=1000`): the class resolutions of the
final class set over the populated classes, the base build's values of phase 0 beside them.

| seed | populated classes | median resolution (A) | mean of the ten best (A), and the ten | classes under 10 A |
| --- | --- | --- | --- | --- |
| 11 | 37 (base 37) | 9.89 (base 9.33) | 6.91 (4.41 5.93 6.80 6.80 7.42 7.42 7.59 7.59 7.59 7.59) (base 6.94) | 20 (base 23) |
| 22 | 39 (base 37) | 9.60 (base 9.07) | 6.40 (4.18 4.29 6.28 6.40 6.53 6.80 7.10 7.42 7.42 7.59) (base 6.92) | 22 (base 21) |
| 33 | 37 (base 37) | 9.89 (base 9.60) | 6.72 (4.41 6.40 6.53 6.66 6.66 6.80 6.80 7.59 7.59 7.77) (base 6.53) | 19 (base 21) |
| 20261103 | 40 (base 39) | 13.70 (base 11.66) | 6.56 (4.80 5.10 6.16 6.28 6.28 6.80 7.25 7.59 7.59 7.77) (base 6.87) | 18 (base 19) |

Over the four seeds the ten best classes are as good or better (6.40 to 6.91 A against the base's 6.53 to
6.94 A), the counts under 10 A overlap (18 to 22 against 19 to 23), and the medians reach past the base's (9.60
to 13.70 A against 9.07 to 11.66 A): the median counts every populated class, and the run with the highest
median populates all 40 classes (the base's highest, 11.66 A, populates 39), so the extra weakest classes move
it. Because the medians alone did not settle it, the readout was taken over four more seeds (44, 55, 66, 77) on
both builds (`scratch/phase1_runs/bgal_extra/`, the same command, the base build from `scratch/keep/build_base`;
`readouts_bgal_extra_extract.txt` under `scratch/keep/phase1/`). Over eight seeds per build:

| statistic over eight seeds | base build | this build |
| --- | --- | --- |
| median class resolution (A): range, mean | 9.07 to 11.66, 9.75 | 8.82 to 13.70, 9.99 |
| mean of the ten best (A): range, mean | 6.53 to 6.97, 6.78 | 6.35 to 6.91, 6.61 |
| classes under 10 A: range, mean | 18 to 23, 20.4 | 18 to 24, 20.9 |
| populated classes: range | 37 to 39 | 36 to 40 |

The two distributions overlap in every statistic; this build's ten best are better by 0.17 A on average and its
count under 10 A higher by half a class, its medians wider on both sides (8.82 A at the low end, the one
13.70 A at the high end with 40 populated classes). The bgal readout is within the seeds' spread.

*trpv1 `refine3D`* (35,645 particles, c2, from its map; `refine=prob objfun=euclid box_crop=160 shift_sigma=0
maxits=3`, seed 20261103): the FSC shells of the assembled map after each iteration, the base's beside them.

| iteration | FSC=0.143 (A) | FSC=0.5 (A) |
| --- | --- | --- |
| 1 | 4.041 (base 4.041) | 6.224 (base 6.102) |
| 2 | 4.095 (base 4.095) | 6.621 (base 6.621) |
| 3 | 4.095 (base 4.095) | 6.765 (base 6.765) |

The second and third iterations' shells are the base's to the last digit printed; the first iteration's
FSC=0.5 shell is 0.12 A behind the base's, more than the few hundredths phase 0 allowed before asking for the
seeds' spread. The first iteration was therefore run over three more seeds on both builds
(`maxits=1`, the same command otherwise; the base build from `scratch/keep/build_base`;
`scratch/phase1_runs/trpv1_extra/`, the shells in `readouts_trpv1_extra.txt` under `scratch/keep/phase1/`):

| seed | base: FSC=0.143, FSC=0.5 (A) | this build: FSC=0.143, FSC=0.5 (A) |
| --- | --- | --- |
| 20261103 | 4.041, 6.102 | 4.041, 6.224 |
| 11 | 4.041, 6.102 | 4.041, 6.102 |
| 22 | 4.041, 6.102 | 4.041, 6.102 |
| 33 | 3.990, 5.872 | 3.990, 6.483 |

At FSC=0.143, the resolution readout, this build gives the base's shell at every seed, seed for seed (the two
builds move together to 3.990 A at seed 33). At FSC=0.5 the shells are quantised by the Fourier shell spacing
(about 0.12 A at 6 A): this build gives the base's shell at two seeds and lies one and three shells behind it at
the other two, while the base's own seeds move it two shells the other way (5.872 A); the two builds' spreads
(5.87 to 6.10 A and 6.10 to 6.48 A) meet at 6.10 A, and by the second and third iterations the shells are the
base's exactly. The trpv1 readout is within the seeds' spread at FSC=0.143 and at the later iterations; the
first iteration's FSC=0.5 shell is the one measure that sits at the spread's edge rather than inside it, on four
seeds per build, and phase 4's measurements on the final build (three runs each) are where it is settled. The shared-memory trpv1 run of level 3 at eight threads (two iterations)
gives the same shells as this distributed readout (4.041 and 6.224 A, then 4.095 and 6.621 A): the keyed draws
make the readout the same across layouts, which the base's thread-dependent draws did not promise.

**Exit criteria of phase 1.** Levels 1 and 2 hold: the generator's vectors are phase 0's in both builds, and
the harness gives identical rows, selections and extras across threads (one, eight and sixteen), batch sizes,
parts and perturbed earlier parts, in five runs. The fast gate passes on both trees (15 of 15 suites each, after
every review round) and `unit_pftc_align2D3D`, `unit_numerics`, `unit_ori` and `unit_project` pass on both. The
truth suites pass with their floors unchanged (`multistate_workflows` on its four existing routes,
`solve3D_addon`, `cont_refine3D_1jxy`); the simulated workflows reproduce their two base failures at the same
points and nothing else. The readouts are within the seeds' spread: the bgal class resolutions over eight seeds per build, the trpv1
FSC=0.143 shells seed for seed and the second and third iterations' shells to the last digit; the one measure at
the spread's edge is the first trpv1 iteration's FSC=0.5 shell (the base's at two of four seeds, one and three
shells behind at the other two), recorded for phase 4's measurements. The reference run directories of this build are kept under
`scratch/keep/identity/phase1/`. Nothing is committed; the working tree holds the phase as one diff. The side
effect to carry forward: every probabilistic run now plans its FFTs with `FFTW_ESTIMATE` (the seeded-run rule of
`has_global_seed()`), to be weighed in phase 4's wall-time measurements. Phase 2 (the pipeline in shared
memory) can start on this tree.

### Phase 2: the pipeline in shared memory (2026-11-04, the Dell)

Done on 2026-11-04 by the autonomous run `prob_pipeline` on the Dell, on the tree phase 1 left (base HEAD
`22a643642`; the run root `/data2/agent_runs/prob_pipeline`, linked as `~/agent_runs/prob_pipeline`, and paths
below that are not absolute are relative to it). The tree built in `repo/build` (Release) and
`scratch/keep/build_debug` (Debug). Phase 1's binaries are kept in `scratch/keep/bin_phase1` and its run
directories under `scratch/keep/identity/phase1/`; both are this phase's comparison points. The journal is
`state/journal_phase_2.md`.

**What changed, in one paragraph.** In shared memory a probabilistic iteration of `refine3D` or `refine2D`
(`refine=prob`, `prob_neigh`) is one batch loop inside the matcher, `refine3D_exec` or `refine2D_exec`, on the
caller's parameters and builder. The matcher draws the iteration's sample and builds the part's table in its
pipeline form, which holds the assignment map and, batch by batch, the rows the fill keeps in memory. It computes
the iteration's model once from the complete incoming field. For each batch it prepares the particles as before and
fills the batch's rows with the existing fill. It draws each particle's pose or class from its own row with the
existing selection rule and lays out the weighted insertion's extras. Then it runs the existing polish, residual
and commit, and in 2D the batch's class sums. The strategies no longer write the particle field, run the
`prob_align` commander and reread the field. That commander built a second parameters object and builder from the
file, ran the table commander as a third, wrote a part file, read it back, assigned and wrote an assignment file.
The layering policy's accepted exception ends on the shared-memory path, and no table or assignment file is
written there. The selection rules are the old `ref_assign` bodies restructured into an opening, a per-row core
over a range of rows and a closing, which the distributed masters run over their whole table exactly as before.
The distributed path is untouched: the table workers write their part files, the masters merge and assign, the
workers read the assignment.

**Files of this phase.** The matchers `src/main/strategies/search/simple_strategy3D_matcher.f90` and
`simple_strategy2D_matcher.f90` hold the alignment context: the sample draw, the references with the diagnostic
halves and their memo, `begin_pipeline`, `pipeline_batch`, `close_assignment`, the context's memory estimate and
the bench lines. The 2D matcher also lost its `release_pass_buffer` call, as nothing fills a 2D pass buffer now.
The tables `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90`
and `simple_eul_prob_tab2D.f90` gained the pipeline constructors, `set_model_from_field`, the fills keeping the
batch's rows, the per-row cores and `open_assignment`, `assign_batch`, `close_assignment`. Their utilities
`simple_eul_prob_tab_utils.f90` took `print_disk_set` and the cone maker `insert_cone_from_fsc` from the strategy
module, and `prob_insert_poses%setup_range` for the batch-by-batch extras layout. The strategies
`src/main/strategies/parallelization/simple_refine3D_strategy.f90` and `simple_refine2D_strategy.f90` lost the
shared-memory iteration's pre-step detour (the field write, the commander, the segment reread, the
`run_prob_align2D` routine), and `simple_prob_align_strategy.f90` lost its two moved routines. The workers'
commander `src/main/commanders/simple/simple_commanders_prob.f90` calls the tables' `set_model_from_field`, and
lost the private model helpers and the 2D worker's shared-memory buffer hand-over. The batch module
`src/main/strategies/search/simple_matcher_ptcl_batch.f90` lost the 2D record paths, as no 2D pass opens the
record buffer now. The 2D batch preparation `simple_strategy2D_alloc.f90` (a row added to the file table) no
longer draws the classic search order under the probabilistic routes. The harness
`src/main/commanders/test/simple_polar_matching_bench.f90` and `simple_commanders_test_highlevel_matchbench.f90`
gained the pipeline route, its three configurations and the assignment's width in the comparison. The decision
schema's consumers (`simple_strategy3D_prob.f90`, `simple_strategy2D_prob.f90`), the sigma2 strategy and the
particle I/O module listed for this phase were not touched: the consumers read the same assignment map and extras
as before. The searches' `prep4prob` item of the file table needed no change either, as the decisions below say.
The working tree's diff over phases 0 to 2 touches 36 files, about 3,990 lines added and 900 removed.

**The pipeline's shape** (section 3.3's operations, as placed). `refine3D_exec` and `refine2D_exec` own the
alignment context under a probabilistic mode in the process that owns the project (`do_pipeline`, `l_pipeline`:
not a distributed worker). The matcher draws the iteration's sample where the classic modes draw theirs: the
sampling history is cleared at the run's first iteration, the fill-in sample is drawn every fifth iteration under
`fillin`, the keyed `sample_ptcls4update3D` or `sample_ptcls4update2D` otherwise, and the disk-set readout
`print_disk_set` follows the draw. It reads the references as the table worker did, with the diagnostic halves
under `prob_post_diag`, and memoizes them: the fill scores with them, and no polar kernel branches on their
presence, so the polish is unchanged. Then `begin_pipeline` opens the context. The part's table takes its
pipeline form (`new_pipeline`, `new_neigh_pipeline`): the assignment map and the thread buffers or the row store,
no part file. The iteration's model comes from the complete incoming field (`set_model_from_field`): the
refined-reference budget over the active states and the shift prior's width, the arithmetic the table workers'
helpers had. The weighted insertion's cone (`topk > 1`) comes from the state FSC files, bounded by the matching
low-pass limit when the command line sets one and by nothing otherwise. The assignment opens
(`open_assignment`): the insertion table, the cone's geometry, the n_eff readout, the draws' scope, the counts.
Each batch is then prepared once (`build_batch_particles3D`, `build_batch_particles2D`) for the fill, the polish,
the residual and the 2D sums. `fill_tab_range` runs on the batch, and in pipeline mode the fill keeps the batch's
rows: the dense 3D fill moves its row and refined-entry arrays into the table, the neighbourhood and the 2D fills
leave their candidates in the thread buffers. `assign_batch` follows. The 3D dense core runs over the kept rows. The
neighbourhood and 2D sparse cores run over a candidate store built from the thread buffers in thread order, the
order the part file gave the master. The 2D dense core runs over rows packed from the buffers as the flush packed
them. The
extras are laid out for the batch (`prob_insert_poses%setup_range`), with the slots `setup` would give from the
complete table, as the insertion table is filled in particle order. The existing search loop runs
(`strategy3D_prob%srch_prob`, `strategy2D_prob%srch_prob` read the same assignment map and extras), then the
sigma2 residual and, in 2D, the batch's class sums in the order of before. After the loop `close_assignment`
prints what the master printed: the n_eff readout, the extras' report or the "no FSC yet" line, the row mode's
share and the particles without a temperature. In the pipeline it also reports the rows without a valid or with a
flat distance, the counts and stops of `warn_invalid_and_flat` taken per row as the batches pass. The strategies' iterations lost
the detour and nothing else. The matching low-pass reset to the command line's value after the pre-step stays. The
canonical sigma2 preparation stays before the matcher: it validates the store and clears stale transaction files,
and the fill and the residual read the same committed spectra. The distributed masters still run
`prob_align_exec`, `prob_align_neigh_exec` and `prob_align2D_exec` on their own projects.

**One selector, two callers.** The master's whole-table `ref_assign`, `ref_assign_neigh` and the 2D `ref_assign`
(the distributed path) are now `open_assignment`, the per-row core over rows 1..n and `close_assignment`. The
pipeline's `assign_batch` is the same core over the batch's rows. The dense 3D core reads a particle's row
through one accessor (`get_row`: the kept batch rows, or the master's `dists` and `inpls` columns and refined
list). It draws the state as the row's closest reference and the direction within its block, commits the drawn
entry with its refined shift or the seed shift, and selects the extras over the same rows: the parallel
eligibility pass per particle, then the insertion table filled in particle order. The neighbourhood core keeps its
fallback, a particle without a valid candidate re-seeded at its previous orientation and scored on the
calculator, which in the pipeline holds the particle. It keeps its state-row gathering and its extras. The 2D
cores keep their class order, their first-minimum tie rule and their peak counts. The n_eff readout offers the
rows it wants in particle order, during the per-row pass instead of before the draws, which selects the same rows.

**The batch capacity** (section 3.5) stays `nthr x BATCHTHRSZ`. The 2D class sums are accumulated per batch in
batch order, which section 3.4 fixes, and the level-3 identity depends on the same batches. The context's memory
is estimated in checked 64-bit arithmetic from the calculator's dimensions: the references of both halves and
their memoized transforms, the batch's polar particles with their CTF terms and memos, the 2D batch's images for
the sums, the batch's rows and refined entries, the assignment map. The matchers' bench files report it against
the particle buffer's budget, which is half the task memory a queue grants, else the node's memory over twice its
processes. A context over budget warns. The measured peaks are below.

**Decisions and findings of the phase.**

*The reseed counter* (section 4.3). The pipeline does not replay the two `seed_rnd` calls the retired commanders
made per iteration. The counter changes as section 4.3 says it would. No probabilistic search draws from the
intrinsic stream: the polish's shift search is L-BFGS-B from a coarse grid. A replay was built and tested while
the second-iteration deviation below was being isolated, changed nothing, and was removed. Whether a classic
stage after a probabilistic one moves with the count could not be tested by the entries. No classic stage follows a
probabilistic one in these single-state runs: stage 1 of `solve2D` precedes them, and the single-state `solve3D`
schedule is probabilistic throughout. A multi-state `solve3D` runs its classic stages before its probabilistic ones
as well. A workflow with a classic stage after a probabilistic one draws that stage from a
differently advanced stream than phase 1, as the plan accepts. This is a side effect to carry forward, beside
phase 1's `FFTW_ESTIMATE` planning rule, which stands and still waits for phase 4's wall-time measurements.

*The classic order draw.* The audit's matcher item of this phase (section 9) is done: the 2D batch preparation
drew a classic search order per particle on every batch, read by the `snhc` strategies only. Under the
probabilistic routes the preparation now leaves the order in class order and draws nothing
(`prep_strategy2D_batch`'s `l_order`, false from the 2D matcher under `prob`). Nothing probabilistic read it, so
no output moves, and the intrinsic stream no longer advances per batch there. Checked on the install with the order
draw skipped, before the cone's bound, against the phase-1 binaries: the band experiment with `lp` and the
two-iteration continuation with `trs=5` are byte-identical in every class average and class document
(`scratch/keep/phase2/final_d/`).

*`prep4prob` from the model.* The file table's item needed no change. The polish's prior width,
`s3D%shift_prior_sigma2` and `s2D%shift_prior_sigma2`, is computed once from the incoming field before the batch
loop (`prep_strategy3D`, `prep_strategy2D_glob`, after the sample draw and the band's adoption). It is the same
`shift_prior_sigma2(field, lp/smpd)` the table's `set_model_from_field` computes from the same field. Nothing of
the polish is derived from the evolving output, which is section 2.2's rule.

*The insertion cone* (`topk > 1`). The cone's rim limit is bounded by the matching low-pass limit when the
command line sets one, which the stage entries do, and by nothing otherwise. That reproduces the old pre-step,
whose fresh parameters held the command line's `lp`, else the parameters' default of 20 A, a bound no map of the
level-3 runs reaches at FSC=0.5. A first version read the strategy's `lp` after the model's band was adopted, and
the review round moved the capture to the matcher's entry. That was still wrong at a later iteration, where the
strategy's `lp` without a command-line value is the previous iteration's adopted band, and the 3D case of the
third deviation found it. The cone maker moved to the tables' utilities with scalar inputs
(`insert_cone_from_fsc(nstates, box, smpd, mskdiam, lp)`), so the utilities take no parameters dependency.

*An empty sample.* A 2D partition without a sampled particle runs no batch and opens no pipeline. The 2D matcher's
context estimate read the map's first element for the entry's size, out of bounds in that case; both matchers now
take `storage_size` of a local entry.

*The harness's blind spot.* The fills harness compared the assignment's reference, in-plane index, distance and
shift but not its posterior shift width (`sdev`), which the polish records into the field's `shsig` for the next
iteration's disk set. It compares the width as well since 15:44 (`final_sdev_nthr8.tsv` and every later run).

*One accessor kept.* The dense 3D table's `dense_candidate` stays as the 3D tester's accessor over the per-row
core (`get_row` and `entry_candidate`). It goes with the global table in phase 4.

*The review round.* An agent read the diff and its findings were applied: the cone's low-pass capture, the empty
2D sample, the cone maker's home, a guard in the extras' layout (`setup_range` before the first batch throws),
stale comments and a shadowed local in the neighbourhood selector. One observation stands: the stream's
`stream2d` override of the record path is unreachable in shared memory now, as the record path is gone.

*The OpenMP rule.* New or changed parallel code was tested at sixteen threads, repeated, against one thread,
through the harness below (the pipeline's per-thread buffers and stores, the per-row cores' parallel passes).

*The old path's three toolbox defects* (below). Each is a value the pre-step's second toolbox did not share with
the strategy. The band without an explicit `lp` is latent from iteration 5 on and surfaces at the first iteration
without a limit. The fresh start's inherited shifts are zeroed in the matcher only. The dynamic shift window is
widened in the strategy's parameters only. The pipeline's own cone bound, found through the 3D case, was the
pipeline's defect and is fixed. The pipeline shares everything by construction. The distributed parts
of phase 3 will have to take the band, the reset field and the window from the master's step rather than from
their command line.

**Level 1 (draw identity)** holds unchanged: nothing in the generator moved, and `unit_numerics` passes in Release
and Debug (1,480 assertions, its "random draws" suite 130 checks).

**Level 2 (frozen-input identity)** holds for the pipeline route. The harness (`polar_matching_bench
suite=fills`) runs each route's fill, merge and whole-table selection as phase 1 did, and now the shared-memory
pipeline route as well (`run_fill(..., l_pipeline=.true.)`). That route builds one table over the particles, the
model from the field as the matchers compute it and the cone, then runs batch by batch the fill with the rows
kept in memory and the batch's selection. The pipeline runs at one thread with one batch, at the run's thread count with three batches and
at one thread with three batches. Every row is compared with the reference route's at one thread, one batch, one
part, by particle and reference identity: the rows' distances and in-plane indices, the refined shifts, the
temperatures, the seeds, the assignment and the extras. The assignment's comparison covers its reference, in-plane
index, distance and shift, and since 15:44 its posterior shift width as well. The table below is the final tree's
eight-thread run (`scratch/keep/phase2/harness/final_e_nthr8.tsv`). The nine configurations are phase 1's six (one
and eight threads, batches of 64 and 22, one and three parts, the earlier parts perturbed) and the pipeline's three.

| route | states | entries per configuration | configurations | largest count in any category |
| --- | --- | --- | --- | --- |
| 3D_dense | 1 | 12,800 | 9 | 0 |
| 3D_dense | 2 | 25,600 | 9 | 0 |
| 3D_neigh | 1 | 704 | 9 | 0 |
| 3D_neigh | 2 | 1,344 | 9 | 0 |
| 2D_dense | 1 | 12,800 | 9 | 0 |
| 2D_dense | 2 | 25,600 | 9 | 0 |
| 2D_sparse | 1 | 9,600 | 9 | 0 |
| 2D_sparse | 2 | 19,200 | 9 | 0 |

Every run of the phase gives 72 rows with every count zero. The runs, by tree:

| tree | threads | table |
| --- | --- | --- |
| pre-review sources | 8 | `phase2_run2.tsv` |
| pre-review sources | 16, twice | `phase2_nthr16_rep1.tsv`, `phase2_nthr16_rep2.tsv` |
| pre-review sources, Debug | 4 | `phase2_debug_prereview.tsv` |
| review-fixed sources | 8 | `final_nthr8.tsv` |
| review-fixed sources | 16, twice | `final_nthr16_rep1.tsv`, `final_nthr16_rep2.tsv`, identical in every column but the wall time |
| review-fixed sources, the width compared | 8 | `final_sdev_nthr8.tsv` |
| review-fixed sources, Debug | 4 | `final_debug_nthr4.tsv` (before the width), `final_c_debug_nthr4.tsv` (with it) |
| the order draw skipped | 8 | `final_d_nthr8.tsv` |
| final sources (the cone's bound as well) | 8 | `final_e_nthr8.tsv` |
| final sources, Debug | 4 | `final_d_debug_nthr4.tsv` |

The first run on the pre-review tree (`phase2_run1.tsv`) had shown 177 differing entries in the 3D neighbourhood
route and 8,272 in the 2D sparse route, in the refined-shift column alone. The assignment, rows, temperatures,
seeds and extras of those configurations were identical. The harness's new gather took a candidate's shift only
when it was refined, while the reference gather takes every candidate's shift from the master's store, zero for
the unrefined. The gather now mirrors each route (`gather_buffers(l_all_shifts)`). The pipeline's fills take what
the reference fills take: 3D dense with one state 5.1 s at one thread and 0.92 s at eight, 2D sparse with two
states 13.7 and 2.5 s (`phase2_run2.tsv`).

**Level 3 (structural identity with phase 1 at a fixed configuration).** `scratch/run_phase2_level3.sh` ran the
same commands, layouts and seed as phase 1's reference directories. It ran twice: on the staged install of the
review-fixed sources (`stage_p2`, sources of 14:57, the table below) and, after the two later source changes (the
order draw skipped, the cone's bound), on the final install (`stage_p2c`, `run_phase2_level3_final.sh`).
On the final install every comparison gives the same counts as in the table, with one change: the trpv1
`refine3D` case is byte-identical in every map and binary on disk (9 mrc, 4 bin), only its project's path
differing, so with the cone's bound corrected that case matches phase 1 at both iterations
(`scratch/keep/phase2/level3_compare_final/`). Each run was compared with its reference directory under `scratch/keep/identity/phase1/`. The
registry's `test=identity` checks the entries' declared outputs, and `compare_dirs.py` every file of every case.
The run roots are normalised and the job records' time lines dropped, because the two sets of directories sit
under roots of different length and the binary project files embed their paths. The projects of the distributed and stream cases were compared segment by segment through
`compare_projects.sh`: every segment the same, except the `solve3D` case's `projinfo`, which carries the run's
time-stamped id (`level3_compare_final/segments.txt`). Both tools are under
`scratch/keep/tools/`. The comparisons and the timelines are kept under `scratch/keep/phase2/level3_compare/` (the
review-fixed install) and `level3_compare_final/` (the final install), with the registry's reports; the run
directories themselves were not kept (14 GB each). The simulated-workflow entries pin one thread, as phase 0 found,
so their "eight threads" layout is nominal.

| case (this build against phase 1's reference) | layout, time, exit | result |
| --- | --- | --- |
| `simulated_workflow` 1jxy, 1 thread | 430 s, exit 1 | identical through iteration 22 of `6_solve2D` (stage 1 classic, stages 2 to 5 probabilistic): 98 mrc, 21 mrcs and 24 star the same, with 10 more star the same after path normalisation; differing from `cavgs_iter023` (stage 6) on, including the whole `7_solve3D`; the registry's identity test: 102 files the same, 42 differing or missing |
| `simulated_workflow` 1jxy, 8 threads | 465 s, exit 1 | the same counts as at one thread |
| `simulated_workflow` 6vxx, 1 thread | 1,461 s, exit 0, 1 FAIL line | identical through iteration 22: 116 mrc, 30 mrcs and 27 star the same, with 16 more star the same after path normalisation; differing from stage 6 on; identity test 116 the same, 70 differing or missing |
| `simulated_workflow` 6vxx, 8 threads | 1,429 s, exit 0, 1 FAIL line | the same counts as at one thread |
| `mini_stream_1jxy`, 8 threads | 62 s, exit 0 | identical below its pid-named directory: 107 mrc, 7 mrcs, 27 star; the project the same in every segment; the sigma2 state's 8 bytes of embedded path and the metadata's times differ; the registry's identity test matches none of its 13 patterns, because the directory is pid-named |
| `solve2D` 2 parts; `solve2D nsample=30` 2 parts | 15 s and 13 s, exit 0 | identical: 82 and 85 mrc, 27 and 28 star, 3 bin each; the projects the same in every segment |
| `solve3D` 2 parts (6vxx, c3) | 1,085 s, exit 0 | identical: 41 mrc, 4 bin, 2 dat; the project the same in every segment; the 104 FSC plots and the manifest differ by their creation times, as in phase 1 |
| trpv1 `refine3D refine=prob`, shared memory, 8 threads, 2 iterations | 1,913 s, exit 0 | on the review-fixed install: iteration 1's resolution file identical (4.04 and 6.22 A, printed to two decimals) and iteration 2 differing, for the cone's bound; on the final install (1,374 s): byte-identical in every map and binary |
| trpv1 `refine2D refine=prob`, shared memory, 8 threads, 2 iterations | 74 s, exit 0 | differs from iteration 1 |

The end-of-stage binaries of the workflows (`cavg_state.bin`, `frcs.bin`, `sigma2_state.bin` of `6_solve2D`, the
four binaries of `7_solve3D`) carry the stages' final state and differ, as stage 6 differs. So the identity holds
on the distributed and stream cases and on iterations 1 to 22 of both workflows at both thread counts, and it does
not hold from stage 6 of `6_solve2D` on, including the `7_solve3D` that follows, nor on the two trpv1 cases. The
simulated workflows end as the base does. 1jxy stops at iteration 41, solve3D stage 3, with `the map's last shell
must lie in the band; new_shiftmap`, exit 1, after 2,449 log lines. 6vxx completes with the one `FAIL: final-volume
correlation is below the required minimum`, exit 0, its whole-volume correlation 0.5611 and its docking correlation
0.8764 (phase 1: 0.5415 and 0.870). These are the entries rerun at eight threads on the review-fixed install
(`scratch/keep/phase2/evidence/`). The three places where identity does not hold are explained below, each
isolated to its mechanism. The band and the window were pinned and the runs made byte-identical; the shifts were
shown by the toolbox prints and by a continuation's identical first iteration.

**The first deviation, and its cause.** The shared-memory workflows are byte-identical to phase 1 through the
first five stages of `6_solve2D`. The first differing file of every one of the four runs is `cavgs_iter023.mrc`.
That is the first iteration of stage 6, the first probabilistic iteration of the workflow without an explicit
`lp`. There the matching band comes from the class FRCs (`frcs.bin`), 7.564 A in the phase-1 and the phase-2
logs. The cause is a defect of the old two-phase shared-memory 2D path, not of the pipeline. The
table worker `exec_prob_tab2D` built its own builder from the project and never read the class FRCs, so its
`clsfrcs` object stayed unconstructed. From iteration 5 on, `set_bp_range2D` without an explicit limit takes
`estimate_lp_for_align` of that object, a huge value, which the band logic turns into the `lpstart` fallback. The
old table was scored at the lpstart band while the matcher of the same iteration polished at the FRC band. The
pipeline scores and polishes in one process on one band, the band the matcher always used. The mechanism was
isolated on the bgal solve2D project after its stage 5 (`scratch/run_phase2_band_experiment.sh`: `refine2D
refine=prob objfun=euclid startit=6 maxits=1 nthr=1` with the stage's `frcs.bin` and class averages, the phase-1
binaries against the pre-review install; the comparisons and logs are kept under
`scratch/keep/phase2/band_experiment/`). The case with `lp` was repeated on the install with the order draw skipped
(`scratch/keep/phase2/final_d/compare_band.txt`), identical again; the case without `lp` was not repeated.
Without `lp` the two differ in every output: the class averages and their halves, the star file, the class state,
the sigma2 state, the stats. With `lp=` set to the matcher's band (4.730 A, printed by all four runs) every output
is byte-identical and the project differs only in its absolute run path. From stage 6 on the shared-memory 2D
trajectory is therefore a new one, and so is the `solve3D` that follows it in the workflows. The distributed 2D
worker has the same defect, as it also builds its own builder without the FRCs. Its fix changes the distributed
path, which this phase leaves untouched, and belongs to phase 3.

**The second deviation: a fresh 2D start.** The purely probabilistic shared-memory `refine2D` case (trpv1, 100
classes, `refine=prob objfun=euclid box_crop=120`, two iterations, from a project whose particles carry the shifts
of an earlier run) differs from phase 1 from its first iteration. The distributed one-part run reproduces the old
shared-memory numbers to the digit on either build (the posterior's n_eff 2.79, 66.1 per cent at the row mode,
`scratch/keep/phase2/r2d_diag/diag4_lines.txt`), so the old shared-memory result is the distributed table worker's.
The two toolboxes were instrumented with temporary prints, applied and reverted (`scratch/diag_toolbox_patch.py`;
the prints are `scratch/keep/phase2/r2d_diag/diag5_distr.txt` and `diag5_shmem.txt`). The old table worker and the
pipeline's matcher agree in every parameter of the fill: box and crop, sampling, the band 14.819 A, the window,
the masks, centring, the prior, the Fourier range, lpstart and the 2D limits. They agree in the calculator's
dimensions, in the first reference's polar transform and in the sigma2 spectra. They differ in one input only, the
particle field's shifts: particle 1 has x -2.710 and y -1.991 px in the worker and 0 and 0 in the matcher. A fresh
2D start zeroes the inherited shifts in memory before the search (`prep_strategy2D_glob`: "zeroed prior to and
consistently with search"), inside the matcher. The old shared-memory path ran the table worker on the project
file as written before that reset. Its table was scored on particles still carrying the earlier run's shifts,
while the polish of the same iteration started from zero: the fill and the polish disagreed on the particles'
frame at every fresh start from a project with shifts. The distributed path has the same defect, as its search
workers zero in their own memory and its table workers never do. The pipeline fills and polishes on the same reset
field, and its first iteration is the better one: the class resolution averages 10.8 A against 20.0 A, and the
shift increments average 6.65 px against 1.84 px, as the particles are actually centred. Neither the crop nor the
centring is involved: the uncropped case and the case without centring differ the same way
(`scratch/keep/phase2/r2d_diag/timeline.txt`). A continuation from an existing iteration (`startit=2 maxits=1`)
is byte-identical in its class averages and class documents, with the particle field identical in every key for
all 35,645 particles.

**The third deviation: the dynamic shift window.** The same continuation run for two iterations (`startit=2
maxits=2`) differs at its second iteration. Every file committed after the identical first iteration is byte-identical between the two paths: the sigma2
state, the FRCs, the class state, the class averages, the class documents. The particle fields at the matcher's
entry are identical at both iterations in the dumped keys, x, y, shsig, shincarg, updatecnt and sampled (a
temporary dump, `scratch/diag_shsig_patch.py`).
The assignment maps of the two paths, dumped after the batch loop (`scratch/diag_assgn_patch.py`), are identical
at the first iteration and differ at the second for 11,760 of 35,645 particles: the distances at the seventh
digit, the refined shifts at the fifth, the posterior widths more
(`scratch/keep/phase2/r2d_diag/DIAG_ASSGN_3_distr.txt` and `_shmem.txt`). With the two toolboxes' parameters
printed per iteration (`diag11_iter3_diff.txt`) they differ in one parameter at the second iteration: the shift
window `trs`, 5.000 in the old table worker and 5.040 in the matcher. The convergence check widens the window
after an iteration whose search fraction reaches its limit, when `trs` is not on the command line and not above
the minimum (`simple_convergence.f90`: `MSK_FRAC` times the mask radius, between `MINSHIFT` and `MAXSHIFT`; 0.07
times 72 px in this 2D case). It does so in the strategy's parameters. The old shared-memory table worker was
rebuilt from the command line every iteration and kept the initial window, while the polish of the same process
used the widened one. From the second iteration on, the old path's tables were scored with a narrower shift
window than its polish. The distributed parts keep the initial window too, as the master's job description sets
`trs` once, so the distributed path is self-consistent but narrower than the shared-memory polish, a layout
difference of the base that phase 3 will meet. The pipeline fills and polishes with the strategy's window. The
solve2D and solve3D stage controllers put `trs` on their command lines, so the workflows' probabilistic stages
never widen it, which is why the entries are identical across every iteration of stages 2 to 5. After the
polish of the second iteration 10,402 particles' shifts differ, the median difference 4.8e-4 px, a third below
1e-4 px, 211 above half a pixel (a different basin after a last-bit difference). The control pins the window:
with `trs=5` on the command line the same two-iteration continuation is byte-identical between the phase-1
binaries and this build in every class average and class document
(`scratch/keep/phase2/r2d_diag/compare_cont2_trs5.txt`). The trpv1
`refine3D` case, which differed at its second iteration on the review-fixed install, had another cause. Its
convergence check widens the window the same way, and the direct `refine3D` call has no `trs`, but with `trs=5`
pinned the case still differed at its second iteration (`scratch/keep/phase2/evidence/compare_r3d_trs5.txt`),
with phase 1's shells to the last digit printed (4.095 A at FSC=0.143 and 6.621 A at FSC=0.5,
`scratch/keep/phase2/memory/r3d_phase2.readout.txt`) and the volumes differing at the bit level. The cause was the
insertion cone's bound. The matcher bounded the cone by the strategy's `lp` at its entry. At the second iteration
without a command-line `lp` that is the previous iteration's adopted band, 4.041 A, against the old pre-step's
default of 20 A, so the pipeline's cone was the narrower one. The cone is now bounded by the command line's limit
only. With that correction and `trs=5` pinned the case is byte-identical to the phase-1 binaries in every map and
binary on disk (9 mrc, 4 bin) and in both iterations' resolution files
(`scratch/keep/phase2/evidence/compare_r3d_trs5_final.txt`), and with nothing pinned, on the final install's
level-3 run, it is byte-identical as well. The 3D case's second-iteration difference was therefore the cone's bound
alone. Why the widened window leaves this case's searches unchanged, when it changes the 2D case's, was not
pursued; the case runs without centring (`shift_sigma=0`).

**Validation of the new trajectory.** The truth suites of ruling R3 on the final install: `multistate_workflows` on its four existing routes passes 129 of 129 metrics (126 PASS lines, no FAIL,
2,099 s), `solve3D_addon` 46 of 46 (47 PASS, 150 s) and `cont_refine3D_1jxy` 135 of 135 (138 PASS, 358 s), the same
counts as phase 1 (`scratch/keep/phase2/truth_final/`). The same three suites had passed with the same counts on the
review-fixed install (`scratch/keep/phase2/truth/`). Every fixed floor is phase 1's. Two run-derived comparison
targets of `cont_refine3D_1jxy`, a control's own symmetric pose error that bounds its continued run, moved with
the trajectory (0.7343 degrees in phase 1, 0.7349 on the review-fixed install, 0.7383 on the final one), and the
multistate file lists its report rows in another order.

*bgal `solve2D` in shared memory* (`mskdiam=180 ncls=40 nsample=1000 nthr=8`, eight seeds, the phase-1 binaries
against the review-fixed install; `scratch/run_phase2_bgal_shmem_seeds.sh`, the readout in
`scratch/keep/phase2/readouts_bgal_shmem_seeds.txt`): the class resolutions of the final class set, with phase
1's extraction. The review-fixed install stands for the final tree here, as 2D has no cone and never read the
order. The layout differs from phase 1's own bgal readout, which ran four parts of six threads, and the
eight-thread shared-memory run is not reproducible run to run (phase 1's record: stage 1's classic draws), so the
pairs compare the two builds' seed spreads rather than one trajectory.

| seed | phase 1: classes, median, best ten, under 10 A | this build: classes, median, best ten, under 10 A |
| --- | --- | --- |
| 11 | 37, 9.33 A, 6.76 A, 20 | 37, 9.89 A, 6.89 A, 22 |
| 22 | 38, 9.33 A, 6.61 A, 22 | 38, 9.20 A, 6.82 A, 21 |
| 33 | 38, 9.33 A, 6.69 A, 20 | 39, 9.89 A, 6.68 A, 21 |
| 20261103 | 39, 9.89 A, 6.41 A, 20 | 39, 9.07 A, 7.03 A, 21 |
| 44 | 37, 9.07 A, 6.94 A, 22 | 36, 9.61 A, 6.72 A, 19 |
| 55 | 40, 11.14 A, 6.72 A, 18 | 39, 9.89 A, 6.67 A, 21 |
| 66 | 38, 9.90 A, 6.89 A, 19 | 38, 9.60 A, 6.94 A, 20 |
| 77 | 39, 9.33 A, 6.78 A, 21 | 38, 9.33 A, 6.68 A, 21 |

Over the eight seeds the medians average 9.67 A for phase 1 (9.07 to 11.14) and 9.56 A for this build (9.07 to
9.89). The ten best average 6.725 A (6.41 to 6.94) and 6.80 A (6.67 to 7.03). The counts under 10 A average 20.25
(18 to 22) and 20.75 (19 to 22). The two distributions overlap in every statistic. This build's medians are
tighter at the high end, its ten best 0.08 A behind on average with a narrower spread, its count under 10 A
higher by half a class: within the seeds' spread.

**Memory per process** (`/usr/bin/time -v`, maximum resident set, the shared-memory trpv1 cases at eight threads,
`scratch/run_phase2_memory.sh`). The base is phase 0's export and phase 1 its kept binaries. This phase is the
review-fixed install, whose matchers and tables are the final tree's but for the order draw and the cone's bound,
neither of which allocates. The budget is the particle buffer's: half the task memory a queue grants, else the node's
memory over twice its processes, which on this node is 62 GB over 2 for one process, about 31 GB. The matchers'
bench lines carry the context's estimate against it.

| case | base | phase 1 | this phase |
| --- | --- | --- | --- |
| `refine3D refine=prob` (35,645 particles, crop 160, 2 iterations) | 5.52 GB | 5.54 GB | 5.47 GB |
| `refine2D refine=prob` (100 classes, crop 120, 2 iterations) | 3.48 GB | 3.46 GB | 1.55 GB |

The 3D peak of this phase is 1 per cent below phase 1's and the base's: the old shared-memory path held the second
toolbox's field, references and table rows beside the strategy's, and the pipeline holds one batch's rows. The 2D peak
is 45 per cent of phase 1's and the base's: the old 2D path kept every sampled particle's record for the matcher in
the pass buffer (35,645 cropped particles of 120 by 120 pixels, about 2 GB), and the pipeline reads each batch once.
Both cases are far inside the budget of about 31 GB per process (`scratch/keep/phase2/memory/`).

**Exit criteria of phase 2.** The plan's first exit point, shared memory byte-identical to phase 1 at both
thread counts on the declared outputs, holds where the old path ran on one set of parameters and fails where it
did not, and the boundary is exact. In shared memory the identity holds, at one and eight threads, on iterations
1 to 22 of both workflows: every class average, class-average stack and star file of stage 1 and of the four
probabilistic stages 2 to 5. It does not hold from stage 6 of `6_solve2D` on, including the `7_solve3D` that
follows, where the registry's identity test fails with 42 and 70 differing or missing files. Nor does it hold on
the two trpv1 shared-memory cases beyond the 3D one's first iteration. On the final install the 3D trpv1 case is byte-identical at both iterations as well, so in
shared memory only the 2D trpv1 case and the workflows from stage 6 on part from phase 1. Each place where it parts is a value
the old pre-step's second toolbox did not share with its strategy: the band without an explicit limit, the fresh
start's shifts, the dynamic shift window. The band and the window were pinned and the runs made byte-identical.
The shifts were shown by the toolbox prints and by a continuation's identical first iteration. The 3D trpv1 case's
second-iteration difference on the review-fixed install was the pipeline's own cone bound, corrected, and the case
is byte-identical on the final install. The run's own ruling was that the stop rule for a structural phase
does not apply, because the cause is a defect of the replaced path, removed by construction within this phase's
files, and not a defect of the pipeline. The distributed path still carries all three defects. The pipeline's
trajectory from stage 6 on is validated as phase 1's was: the truth suites pass with phase 1's counts on both installs (129, 46 and 135 metrics), the simulated
workflows' two base failures reproduced at the same points, the bgal readout within the seeds' spread. The second
exit point, the distributed path unchanged, holds. Its three cases are byte-identical to phase 1 in every map,
class average, star file and binary, and in the 2 dat of the `solve3D` case, whose FSC plots and manifest differ
by their creation times as in phase 1. The stream's mini pipeline is identical in every map, class average and
star file, with 8 bytes of embedded path in its sigma2 state. The projects of the distributed and stream cases are the same in every segment, the
`solve3D` case's time-stamped run id apart. The third exit point, memory per
process measured and within the budget, holds: 5.47 GB for the 3D case and 1.55 GB for the 2D case, against
5.54 and 3.46 GB in phase 1 and a budget of about 31 GB. The shared-memory strategies call the pipeline on their
own parameters and builder. There is no commander detour, no second toolbox, no segment write and reread, and no
table or assignment file on the shared-memory path. The three files of a run directory, `dist_part*.dat`,
`assignment_part.dat` and `insert_poses.dat`, are the distributed path's alone now. Level 2 holds in Release at
one, eight and sixteen threads and in Debug at four, with the assignment's width compared in the runs since
15:44. The fast gate passes on both trees after every round, 15 of 15. The four unit suites pass on both trees:
`unit_pftc_align2D3D` (1,516 assertions), `unit_numerics` (1,480), `unit_ori` (1,047) and `unit_project` (2,979),
which holds the `solve3D` chunks plan suite; the stream's chunk suites are in the gate's `unit_stream`. Two side
effects carry forward to phase 4 beside phase 1's `FFTW_ESTIMATE` planning rule: the reseed counter's change,
which moves a classic stage that follows a probabilistic one, and the classic order no longer drawn under
`prob`. Phase 1's note that the part files' seed-shift section should be written in particle order moves to
phase 3 with the part files, which this phase did not touch. Nothing is committed; the working tree holds phases
0 to 2 as one diff. Phase 3 (the distributed parts on the pipeline) can start on this tree, with the three toolbox
defects to remove from the distributed path as its first items.

### Phase 3: the distributed parts on the pipeline (2026-11-04, the Dell)

Done on 2026-11-04 by the autonomous run `prob_pipeline` on the Dell, on the tree phase 2 left (base HEAD
`22a643642`; the run root `/data2/agent_runs/prob_pipeline`, linked as `~/agent_runs/prob_pipeline`, and paths
below that are not absolute are relative to it). The tree built in `repo/build` (Release) and
`scratch/keep/build_debug` (Debug). Phase 1's run directories under `scratch/keep/identity/phase1/` are the level-3
references, the phase-1 and phase-2 binaries (`scratch/keep/bin_phase1`, `bin_phase2`) the controls' comparison
points. The journal is `state/journal_phase_3.md`.

**What changed, in one paragraph.** A distributed probabilistic iteration of `refine3D` or `refine2D` schedules
its parts once. The master draws the iteration's sample itself, keyed, in its particle field, prints the disk set
and writes the field, as the classic modes' master does. Every part runs the same batch pipeline the shared-memory
matcher runs: it reproduces its share of the master's sample, reads the references with their halves, fills its
rows batch by batch, draws each particle's pose or class from its own row, polishes, computes the residual, commits
and writes its documents, and publishes a small record of the readouts an assignment prints. The master reads the
parts' records in part order after the scheduling phase, sums them and prints the readouts the owner prints in
shared memory. The table workers, their part files, the master's merge and assignment step, the assignment file and
the insertion-pose file are no longer run on any path; their routines and programs stay in the tree until phase 4
retires them. The n_eff readout chooses its rows by rank in the global sample, the same rows whatever the
partitioning. The 2D master's job description now carries the iteration's tempering factor and, when the command
line sets no shift window, the master's current window, which the parts fill and polish with.

**Files of this phase.** The tables' utilities
`src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` hold the readout record (`prob_readout`, its
file name, writer, reader, the part-order sum, the report and the master's merge), the n_eff readout's rows by global
rank and its sums, the insertion report from counts, and the cone's report switch. The three tables
`simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90` and `simple_eul_prob_tab2D.f90` open their assignment with
the global sample's size and the part's rank offset and close it into the record, which a part writes and the owner
reports. The sampling module `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` makes the masters draw
under the probabilistic modes and holds the two draw helpers the owner's matcher and the master share. The matchers
`simple_strategy3D_matcher.f90` and `simple_strategy2D_matcher.f90` run the pipeline on a part too, reproduce the
master's sample there, open the cone without reporting it, publish the record (a zero one for an empty part) and lost
their assignment-file branches. The strategies `src/main/strategies/parallelization/simple_refine3D_strategy.f90` and
`simple_refine2D_strategy.f90` lost the pre-step and its second scheduling phase, draw the sample, report the cone
(3D), set the 2D job description's keys and merge the records. The commanders of the table workers and the masters'
step modules (`simple_commanders_prob.f90`, `simple_prob_align_strategy.f90`, the `refine3D`/`refine2D` commanders,
the registry) needed no change: the programs are unreachable from production now and go in phase 4. The phase's diff
over the nine files is about 870 lines.

**The masters' one phase.** The distributed `refine3D` master's iteration ran, under `refine=prob` or `prob_neigh`,
the `prob_align` commander on its own parameters and builder: that drew the sample, wrote the field, scheduled the
table workers, merged their part files into a global table, drew every particle's assignment, wrote the assignment
file and the insertion poses, then the master scheduled the search workers, which read those files. The 2D master
did the same through `prob_align2D`. Now the master's iteration calls the draw helper the shared-memory matcher
calls (`sample_ptcls4prob3D`: the sampling history cleared at the run's first iteration, the fill-in sample every
fifth iteration under `fillin`, else the fractional-update draw, then the disk-set readout; `sample_ptcls4prob2D`:
the same with the 2D draw, every particle on a full update), writes the field and schedules its parts once with the
existing queue, job description, scripts and markers. The `master_draws_sample3D` and `master_draws_sample2D`
predicates return true under the probabilistic modes whatever the update fraction, because the parts reproduce the
master's sample from the field's stamps, which must then be the master's. The 3D master already set the iteration's
`prob_beta` in the job description; the 2D master sets it now, and, under the probabilistic modes when `trs` is not
on the command line, its current shift window every iteration, so the parts fill and polish with the window the
shared-memory matcher has (the classic modes keep their one-time rule). The 2D full-update override of the pre-step
is applied to the draw (`l_updatefrac` false on a full update; the parts' `update_frac=1.0` key stays).

**The parts on the pipeline.** One switch kept a worker off the pipeline (`do_pipeline`, `l_pipeline`: the
probabilistic mode and not a distributed worker); it is the mode alone now. A part reproduces its share of the
master's sample over its range (`sample4update_reprod`, an empty share allowed), prints neither the disk set nor the
insertion cone's rim limits (the master prints both after its draw: the disk set from the sample, the cone from the
same FSC files the parts open it from, bounded by the command line's limit as the matcher bounds it), reads the model
with the diagnostic halves and memoizes the references as the owner does, and runs `begin_pipeline`, the batch loop
and `close_assignment` unchanged. The worker's old branch, the assignment map read from the file and the insertion
table read and laid out from the master's file, is deleted; the readers stay in the tables until phase 4 with the
writers and the programs.

**The readouts.** An assignment printed the posterior's n_eff readout, the insertion report (or the "no FSC yet"
line), the row-mode share, the warnings on missing temperatures, on invalid and on flat rows (and stopped when every
row was invalid or flat), and in the neighbourhood route the particles kept at their orientation. A part cannot print
them (its output is wiped at the next schedule) and per-part means would be means of means. `close_assignment` now
fills a record of counts and sums: the particles and references, the row-mode count, the missing temperatures, the
invalid and flat rows, the kept particles, `topk` and the extras flag, the insertion table's counts (particles, extra
poses, the committed weights' sum, the eligible sum, the capped count) and the n_eff readout's kept rows with its
sums of exp(entropy) at the draw's factor and at factor 1. A distributed part writes the record as
`prob_readout_part<NN>.txt` (an empty part a zero record); the owner reports its own record. After the scheduling
phase the master reads the records in part order, sums them (`add_prob_readout`), prints the same lines from the
sums (`report_prob_readout`), stops on the global all-invalid or all-flat counts of the dense 3D route (a part only
counts), and deletes the records. The n_eff readout's rows are chosen by rank in the global sample (plan section 5):
with the stride the global sample's size over 1,000, every stride-th rank among the first 1,000 strides is wanted,
valid or not, and a wanted row is kept when it has enough valid distances and a positive temperature. A part computes
its rank offset from the field (the master's sample reproduced over the whole project, the indices below its range
counted). In shared memory the rank is the local index, and the rows are the old rule's whenever every wanted row is
valid; the old rule skipped on to later rows after an invalid one, so the readout's row set changes only there. The
insertion report's committed-weight sum stays the single-precision sum the report always took, so the owner's line is
unchanged.

**Decisions of the phase.** *The readout's rows.* The n_eff readout chooses its rows by rank in the global sample,
valid or not; the old rule skipped on to later rows after an invalid one, so the readout's row set changes only where
a wanted row is invalid, and the printed values are the owner's in every level-3 case (the readout lines of the
control are identical). *The 2D job description.* `prob_beta` is set per iteration as the 3D master always did; the
parts draw now, with the iteration's factor. Under the probabilistic modes the master's current shift window is set
per iteration when the command line sets none, so a part fills and polishes with the window the shared-memory matcher
has; the classic modes keep their one-time rule. *The all-invalid and all-flat stops.* A part counts its invalid and
flat rows; the master stops on the global counts of the dense 3D route, as its whole-table assignment did; a part
alone no longer stops. *An empty part.* It publishes a zero record with the type's defaults, and the merge takes the
reference count and `topk` from the parts that assigned (the review round found the first version carrying the part's
own defaults into the merge's maximum). *The stream's pool.* It runs the classic 2D routes only; a probabilistic part
of it would find no master sample, and the 2D matcher now stops on that combination. *The cone's report.* The level-4
readout comparison found the insertion cone's rim-limit line gone from the distributed logs (the parts printed it
where they open the cone, into logs the next schedule wipes); the master now opens the same cone from the same FSC
files after its draw and reports it, and a part opens it silently. The line is a report, not a number: the third build
differs from the second by it alone. *The retired programs.* The table workers, the masters' step and their files are
unreachable from production; their code, the readers and writers and the private programs go in phase 4, with the
`DIST_FBODY`, `ASSIGNMENT_FBODY` and `INSERT_FBODY` names. *The review round* (an agent reading the diff): the empty
record's defaults, the record reader's key check (a stale or foreign record is refused), the stream guard, an unused
import, stale comments and the bench labels of the retired pre-step; it verified the rank arithmetic, the record's
fields, the merge, the draws' equivalence with the retired pre-step, the empty paths, the job description's keys, the
run's Fortran rules and that no classic path changed.

**Level 1 (draw identity)** holds unchanged: nothing in the generator moved, and `unit_numerics` passes in Release
and Debug (1,480 assertions) on every build of the phase.

**Level 2 (frozen-input identity)** holds unchanged for the four routes: the harness's 72 rows (the reference and the
pipeline routes at one and eight threads, one and three batches, one and three parts) give every difference count zero
on the three builds of the phase (the first, the one with the review fixes, the final one with the cone report), in
Release at eight threads and in Debug at four (`scratch/keep/phase3/harness/`). The phase changes no OpenMP directive
and no code that runs inside a parallel region (the n_eff readout's offer is sequential), so the 16-thread rule for
new or changed OpenMP code has nothing to apply to; phase 2's 16-thread runs stand. The fast gate passes 15 of 15 on
both trees after every round and the four unit suites pass on both after every round: `unit_pftc_align2D3D` (1,516
assertions), `unit_numerics` (1,480), `unit_ori` (1,047) and `unit_project` (2,979). The readout's unit tests in
`unit_pftc_align2D3D`'s table suite cover the rows by global rank (every row wanted at stride 1, the stride's rows,
the rank offset's effect) as before.

**Level 3 (structural identity).** `scratch/run_phase3_level3.sh` ran the same commands, layouts and seed as the
phase-1 reference directories on the staged install of this phase's sources (`stage_p3`; the review fixes that
followed touch no production output, as the controls below show), compared with `compare_dirs.py`, the registry's
`test=identity` on the entries and `compare_projects.sh` on the projects (the run roots normalised, the time lines
dropped). The comparisons are kept under `scratch/keep/phase3/level3_compare/`, the run directories were not.

| case | layout, time | against phase 1 |
| --- | --- | --- |
| `solve3D` (6vxx, c3), 2 parts, 1 thread per part | 539 s | byte-identical: 41 mrc, 4 bin, 2 dat the same; the project the same in every segment but the run id; the two neighbourhood part files absent |
| `solve2D` (1jxy, 4 classes), 2 parts, 1 thread per part | 8 s | identical through iteration 22 (66 mrc, the 27 star of stages 1 to 5); parting at `cavgs_iter023`, stage 6 |
| `solve2D nsample=30`, 2 parts, 1 thread per part | 7 s | the same: identical through iteration 22, parting at stage 6 |
| `simulated_workflow` 1jxy, 1 and 8 threads | 229 s, 244 s | the counts of phase 2: 98 mrc, 21 mrcs, 24 star the same, parting at stage 6; identity test 102 the same, 42 differing or missing; the base failure at iteration 41 |
| `simulated_workflow` 6vxx, 1 and 8 threads | 857 s, 902 s | the counts of phase 2: 116 mrc, 30 mrcs, 27 star the same, parting at stage 6; identity test 116 the same, 70 differing or missing; the base's final-volume FAIL |
| `mini_stream_1jxy`, 8 threads | 46 s | identical below its pid-named directory: 107 mrc, 7 mrcs, 27 star |
| trpv1 `refine3D refine=prob`, shared memory, 8 threads, 2 iterations | 1,292 s | byte-identical: 9 mrc, 4 bin |
| trpv1 `refine2D refine=prob`, shared memory, 8 threads, 2 iterations | 118 s | differs from phase 1 as phase 2 did (the fresh start's shifts); byte-identical to the phase-2 binaries run beside it: 6 mrc, 2 star, 3 bin |

The exit asks the distributed cases to be byte-identical to phase 1 at their layout. The `solve3D` case is. The two
`solve2D` cases are through every stage whose command line sets the low-pass limit (stages 1 to 5, iterations 1
to 22) and part at the first iteration of stage 6, which sets none: there the retired table worker scored its table
at the `lpstart` fallback, never having read the class FRCs, while the search worker polished at the FRC band
(phase 2's first deviation). The parts on the pipeline score and polish at one band, so the distributed trajectory
from stage 6 on is the shared-memory one of phase 2, and the deviation is at the same point, for the same cause,
in both layouts. The control that pins the mechanism in the distributed layout: a trpv1 `refine2D` continuation at
two parts with `lp=14.819 trs=5` (`startit=2 maxits=1`, one phase-1 iteration as the start), where none of the
three toolbox defects of phase 2 applies, is byte-identical between the phase-1 binaries and this build in every
map, binary and document (`scratch/keep/phase3/control2D/`), and the readout lines the two paths print (the n_eff
readout, the row-mode share, the disk set) are identical to the character. The old two-phase distributed path and
the parts on the pipeline agree bit for bit wherever the old path agreed with itself.

*The final build.* The three distributed cases were rerun on the final build's staged install
(`scratch/run_phase3_level3_final.sh`, kept under `scratch/keep/phase3/level3_final/` with the masters' logs) after
the review fixes and the cone report: `solve3D` byte-identical to phase 1 (41 mrc, 4 bin; the project the same in
every segment but `projinfo`'s run id), the two `solve2D` cases identical through iteration 22 and parting at stage 6
as above. The masters' readout lines (the n_eff readout, the insertion reports and the cone's rim limits, the
assignment share, the disk set, the warnings, the resolution lines) were compared with phase 1's logs to the
character: every value is the same, including the 103 cone lines of the `solve3D` case; what differs is the 104
schedule lines of the retired table workers (`PROB_TAB: 2 parts completed`, `PROB_TAB_NEIGH: ...`), gone, and the
place of the warnings on invalid and flat rows inside the iteration's block, printed with the merged report now rather
than during the master's assignment.

**The merged readouts and the iteration's time.** A distributed probabilistic iteration prints what the shared-memory
one prints, from the parts' records summed in part order: the n_eff readout at the draw's and at the physical
tempering (`POSTERIOR DRAW N_EFF: 49.54 at beta 0.30, 44.37 at beta 1.00 (103 rows)` in the `solve3D` case's first
iteration, the owner's line to the character), the insertion report (`PROB_TAB_INSERT (topk 3)`: the pose counts, the
capped and the eligible, or `no FSC yet, committed pose only`), the assignment report (`PROB_TAB_ASSIGN: 103 particles
x 500 references, 9.7% at the row mode`; in 2D `PROB_TAB2D_ASSIGN`), the warnings on invalid, missing-temperature and
flat rows, and the dense 3D route's stop on a table with no valid or no peaked row, now on the global counts. The
master prints the disk set (`% PARTICLES SAMPLED THIS ITERATION`) from its own draw and, under `topk > 1`, the
insertion cone's rim limits (`PROB_TAB_INSERT CONE (state 1): rim limit 124.80 A, i.e. 79.45 deg of rotation or 96.00
px of shift`) from the FSC files of the previous iteration, which it wrote; a part opens the same cone without
reporting it (the layouts' readout comparison of level 4 found the line missing from the first two builds' distributed
logs, where the parts had printed it into logs the next schedule wipes). Nothing is left in the run directory: no
`dist_part*.dat`, no `assignment_part.dat`, no `insert_poses.dat`, and the records are deleted once merged. The
retired second scheduling phase cost the distributed cases of level 3 about half their time: `solve2D` 8 s against 15
s in phase 1 at the same layout, `solve3D` 539 s against 1,085 s (the plan's section 2.5 measured the trpv1 `refine3D`
pre-step at 310 to 354 s of each 344 to 389 s iteration at four parts of six threads; phase 4 measures the final build
with the machine quiet). The bench of a distributed master now records the draw (`refine3D sample draw`, `refine2D
sample draw`) where it recorded the pre-step, and the `PROB_ALIGN*` and `PROB_TAB*` bench files end on the distributed
path.

**Level 4 (layouts).** *The fixed-input accumulators within the frozen tolerance.* `scratch/run_phase3_level4.sh`
(part a) ran phase 0's accumulator runs on this build's staged install with `SIMPLE_SEED=20261103`: `reconstruct3D`
from the committed `ptcl3D` field of the trpv1 solve3D project and `make_cavgs` (100 classes) from the `ptcl2D`
field of its solve2D project, in shared memory at six and at one thread and at `nparts` 1, 2 and 4 with six threads
per part (`make_cavgs` in shared memory needs `refs=`, as phase 0 found: its shared-memory layouts are the one-part
runs at six and at one thread; the bisect of this phase showed the refusal on the base, phase-1, phase-2 and phase-3
builds alike). `scratch/measure_level4.sh` compared the pairs with phase 0's tool (`scratch/keep/phase3/level4/`).

| pair of layouts | reconstruction: largest voxel difference in the mask (the maps' largest value 10.10) | correlation | FSC=0.5 and 0.143 shells | class averages: largest pixel difference over the 100 classes | lowest class correlation |
| --- | --- | --- | --- | --- | --- |
| shared memory, six threads against one thread | 0 | 1 | 2.43 A, 2.43 A | 0 (the one-part runs at six and one thread) | 1 |
| shared memory against one part, six threads | 0 | 1 | 2.43, 2.43 | (the same run) | |
| one part against two parts | 6.25e-5 | 1.00000000 | 2.43, 2.43 | 7.79e-7 | 1.00000000 |
| one part against four parts | 8.58e-5 | 1.00000000 | 2.43, 2.43 | 6.33e-7 | 1.00000000 |
| two parts against four parts | 2.57e-5 | 1.00000000 | 2.43, 2.43 | 5.96e-7 | 1.00000000 |

Every number is phase 0's base measurement to the digit (`reconstruct3D` and `make_cavgs` are untouched by the run,
so the base's layout-to-layout differences recur exactly), inside the frozen tolerance of 1.72e-4 and 1.56e-6 with
the correlations at 1 and the FSC shells at Nyquist.

*The probabilistic runs across layouts* (part b of the script, `scratch/phase3_runs/level4`, measured by
`scratch/measure_phase3_layouts.sh`, kept under
`scratch/keep/phase3/level4/layouts_*.txt`): the three distributed
cases of level 3 (`solve2D`, `solve2D nsample=30` on the 1jxy project; `solve3D` on the 6vxx one) and the trpv1
`refine3D refine=prob` of level 3 (two iterations, crop 160, `shift_sigma=0`), each in shared memory and at 2, 3 and
4 parts, at one and at eight threads per process, on the staged install of the second build (the final build's
change is a report line). Section 8.4 asks these runs for the accumulators' tolerance above and the truth gates, not
for matching trajectories.

*The 2D cases* (`solve2D` on 82 particles into 4 classes, 27 or 28 iterations over 6 stages, stage 1 classic). Every
layout and thread count is its own trajectory from iteration 1: the classic stage's draws are not keyed and fall in a
different order per part and per thread (the base's behaviour; phase 1 recorded it for the bgal shared-memory run).
The final class sets of the eight runs per case correlate at 0.93 to 0.97 class for class across the part counts
at either thread count, the readout lines of the probabilistic stages differ in 2 to 8 of 22 between a layout and
shared memory, and every layout ends stage 6 at the same band (124.8 A, the tiny case's last band). One thread
against eight at the same layout differs from the first class average on.

*`solve3D`* (103 particles, c3, 4 stages, every stage probabilistic from a random start). Across the part counts and
across the thread counts the search is the same through stage 2: the readout lines (n_eff, the insertion report,
the row-mode share) are identical over the first 42 iterations (126 of 376 lines between shared memory and two
parts; 86 of 273 for two parts at one and eight threads), and the stage-1 and stage-2 maps agree to the accumulators'
order-of-summation differences (stage 1: 1.8e-6 shared memory against two parts, 8.9e-7 two against three parts,
4.2e-7 one against eight threads; stage 2: 2.0e-5, 1.1e-5, 1.1e-5; the maps' largest value 0.77 and 14.3). The
ab initio trajectory then amplifies them: the stage-3 maps correlate at 0.9999990 to 0.993, the stage-4 maps at 0.80
to 0.99, and the final maps of the eight runs at 0.72 to 0.97 with the FSC between two layouts' maps anywhere from
4.1 A to the lowest shell (one eight-thread three-part run found a different solution). Section 8.4 asks no more
of these runs. On the base build the same `solve3D` in shared memory at one and at eight threads
(`scratch/run_phase3_base_threads.sh`, `scratch/keep/phase3/level4/base_thread_control_timeline.txt`) parts at its
first readout line and its stage-1 maps correlate at 0.996 (the base's thread-dependent draws); this build's search is
thread- and part-independent, and what remains is the accumulators' floating-point order, which the plan leaves as
it is (section 8, the reproducible-reduction contract not required).

*trpv1 `refine3D`* (35,645 particles, c2, from its map, two iterations). This is the case with no classic stage and no
random start, and the layouts agree as the keyed draws and the master's sample promise. At eight threads the four
layouts print the same FSC shells after each iteration (4.041 and 6.224 A, then 4.095 and 6.621 A: phase 1's
distributed readout and the shared-memory shells of level 3), the same n_eff readout (3.27 at the draw's tempering and
1.37 at factor 1 over 1,000 rows, then 2.62 and 1.28), the same row-mode share (67.4 then 74.9 per cent) and the same
insertion report (1.17 or 1.18 poses per particle, committed weight 0.92, 28.7 per cent capped). Their final maps
agree to a largest voxel difference of 2.2e-2 to 2.7e-2 inside the mask (the maps' largest value 8.28), a correlation
of 0.99998 and an FSC between any two layouts' maps at the Nyquist shell (3.89 A): the first iteration's search is the
same at every layout, and the second's differs where the assembled map and the FSC it reads differ by the
accumulators' order. At one thread per process the picture across the part counts is the same: the shells phase 1's at
every layout, the final maps of the 2, 3 and 4-part runs agreeing to 2.1e-2 to 2.5e-2 (correlation 0.99998 to 0.99999)
and the shared-memory run's to 2.5e-2 to 2.9e-2 (0.99998), the FSC between any two at Nyquist. Between one thread and
eight at the same layout the maps differ more, by 0.13 of 8.3 at every layout (correlation 0.9996), and the difference
begins in the first iteration: its row-mode share is 67.3 per cent at one thread against 67.4 at eight, the n_eff
readout and the shells the same. A probe pins where (`scratch/run_phase3_thread_probe.sh`: the two-part layout for one
iteration at one and at eight threads per part, `scratch/keep/phase3/level4/thread_probe_timeline.txt` and the two
fields beside it): after one iteration the discrete decisions are all but identical (the projection direction differs
in 25 of 35,645 particles, the in-plane index in 27), the continuous values differ in their last digits for most
particles (the shifts in about 18,000 rows, the in-plane angle in 6,409) and by more than 1e-3 degrees or pixels for
about 1,200 particles (3 per cent), the parts' partial reconstructions differ, and the assembled map by 2.6e-2
(correlation 0.99998), the assembly's FSC-area diagnostic with it. Level 2 holds the fill, the selection and the
extras identical across threads on frozen prepared input, and the sample and the n_eff rows are the master's; what the
thread count moves is the prepared input itself, through the preparations' threaded arithmetic, and from it the
polish's optimiser. That is the reproducible-reduction, FFT and preparation contract section 8 leaves as it is, which
the base's thread-dependent draws kept out of sight; this phase records it and changes nothing there.

**The truth suites and the gate entries** (ruling R3; `scratch/chain_phase3_truth.sh` on the second build's staged
install, one entry at a time at its registry thread count, under the layouts' load; the logs and metric files under
`scratch/keep/phase3/truth/`): `multistate_workflows` on its four existing routes passes 129 of 129 metrics (126 PASS
lines, no FAIL; 5,384 s), the floors the same as phase 2's file metric for metric, 36 of the 129 values moved with the
distributed trajectory; `solve3D_addon` 46 of 46 (47 PASS lines; 2,453 s), its floors and its 46 values phase 2's to
the digit (the addon runs in shared memory, whose path this phase does not change); `cont_refine3D_1jxy` 135 of 135
(138 PASS lines; 410 s), every fixed floor phase 2's, the two run-derived targets (a control's own symmetric pose
error bounding its continued run) at 0.7370 degrees against 0.7383 in phase 2, nine values moved with the distributed
trajectory; the distributed gate entry `flex_pca_blobs` (two parts) 94 of 94 metrics (99 PASS lines; 511 s); the
stream entries `stream_preproc` (95 PASS lines, 5 movies validated against the simulation truth; 82 s) and
`stream_chains` (`PASS: sieve to 3D (11 checks)`, `PASS: movies to 3D (16 checks)`; 790 s), whose chains run the
pool's classic 2D routes and the initial `solve2D` at one part through the stream commanders. Every fixed floor is
phase 1's. The two simulated workflows reproduced their base failures at the same points in level 3 (solve3D stage 3,
iteration 41 of 1jxy; the final-volume FAIL of 6vxx), and `mini_stream_1jxy` passed there.

**The real-data readouts at the distributed layout** (section 12's level 4; `scratch/run_phase3_readouts.sh` on the
final build's staged install, the same commands, layout and seeds as phases 0 and 1, four parts of six threads, with
the one-thread trpv1 layouts running beside; the extraction in `scratch/keep/phase3/readouts/extract.txt`).

*bgal `solve2D`* (2,756 active of 5,513 particles, 40 classes, `nsample=1000`), the class resolutions of the final
class set over eight seeds, beside the base's values of phase 0 and phase 1's at the same layout (its four seeds and
its four extra ones; `scratch/keep/phase1/`):

| seed | base: classes, median, best ten, under 10 A | phase 1 | this build |
| --- | --- | --- | --- |
| 11 | 37, 9.33 A, 6.94 A, 23 | 37, 9.89 A, 6.91 A, 20 | 38, 9.76 A, 6.69 A, 19 |
| 22 | 37, 9.07 A, 6.92 A, 21 | 39, 9.60 A, 6.40 A, 22 | 39, 10.88 A, 6.72 A, 19 |
| 33 | 37, 9.60 A, 6.53 A, 21 | 37, 9.89 A, 6.72 A, 19 | 37, 9.60 A, 6.61 A, 22 |
| 20261103 | 39, 11.66 A, 6.87 A, 19 | 40, 13.70 A, 6.56 A, 18 | 39, 11.66 A, 6.59 A, 19 |
| 44 | 37, 9.33 A, 6.78 A, 20 | 36, 9.20 A, 6.59 A, 23 | 36, 9.46 A, 7.07 A, 19 |
| 55 | 39, 9.60 A, 6.62 A, 20 | 38, 9.46 A, 6.61 A, 24 | 39, 9.60 A, 7.07 A, 20 |
| 66 | 38, 10.36 A, 6.97 A, 18 | 36, 8.82 A, 6.35 A, 20 | 39, 9.60 A, 6.50 A, 22 |
| 77 | 37, 9.07 A, 6.57 A, 21 | 39, 9.33 A, 6.75 A, 21 | 37, 9.60 A, 6.69 A, 22 |

| statistic over eight seeds | base build | phase 1 | this build |
| --- | --- | --- | --- |
| median class resolution (A): range, mean | 9.07 to 11.66, 9.75 | 8.82 to 13.70, 9.99 | 9.46 to 11.66, 10.02 |
| mean of the ten best (A): range, mean | 6.53 to 6.97, 6.78 | 6.35 to 6.91, 6.61 | 6.50 to 7.07, 6.74 |
| classes under 10 A: range, mean | 18 to 23, 20.4 | 18 to 24, 20.9 | 19 to 22, 20.3 |
| populated classes: range | 37 to 39 | 36 to 40 | 36 to 39 |

The three distributions overlap in every statistic: this build's medians sit inside the base's range with a mean
0.27 A past it, its ten best are 0.03 A behind the base's mean and 0.13 A behind phase 1's, its counts under 10 A
are the base's to within a tenth of a class, and its populated classes are phase 1's. The bgal readout is within
the seeds' spread.

*trpv1 `refine3D`* (35,645 particles, c2, from its map; `refine=prob objfun=euclid box_crop=160 shift_sigma=0
maxits=3`, seed 20261103, four parts of six threads; 1,439 s beside the one-thread layouts): the FSC shells of the
assembled map after each iteration, phase 1's and the base's beside them.

| iteration | FSC=0.143 (A) | FSC=0.5 (A) |
| --- | --- | --- |
| 1 | 4.041 (phase 1 4.041, base 4.041) | 6.224 (phase 1 6.224, base 6.102) |
| 2 | 4.095 (phase 1 4.095, base 4.095) | 6.621 (phase 1 6.621, base 6.621) |
| 3 | 4.095 (phase 1 4.095, base 4.095) | 6.765 (phase 1 6.765, base 6.765) |

Every shell is phase 1's to the last digit printed, at every iteration: the parts on the pipeline give the readout
the retired two-phase path gave at this layout, and the shared-memory run of level 3 and the eight-thread layouts
of level 4 give the same first two iterations. Against the base the record is phase 1's: the resolution readout
(FSC=0.143) the base's at every iteration, the first iteration's FSC=0.5 shell one Fourier shell behind it
(6.224 against 6.102 A, the measure phase 1 found at the seeds' spread's edge over four seeds per build and left
to phase 4's measurements), the later iterations' shells the base's exactly. The master's merged readouts over the
three iterations: n_eff 3.27, 2.62 and 2.51 at the draw's tempering (1.37, 1.28, 1.24 at factor 1) over 1,000
rows; 67.4, 74.9 and 76.8 per cent at the row mode; the cone's rim 6.22 then 6.62 A; 1.18 then 1.17 poses per
particle, 28.7 then 28.2 per cent capped.

**`test=nightly`** on the final build's staged install (`scratch/phase3_runs/nightly`, the report kept under
`scratch/keep/phase3/nightly/`): 19 executions in 3,397 s (started while the one-thread trpv1 layouts and the
thread probe ran), 18 passed with every metric at or past its floor (1,067 metric rows), one failed:
`simulated_workflow_1jxy` stops at solve3D stage 3, iteration 41, in `new_shiftmap` ("the map's last shell must lie
in the band"), the base failure section 12 records and every phase's level 3 reproduced at the same point.
`simulated_workflow_6vxx` passes its gate with exit 0 and prints the recorded `FAIL: final-volume correlation is below
the required minimum` of the base. The entries that pass include the distributed ones (`flex_pca_blobs`,
`cont_refine3D_1jxy`, the simulated workflow's distributed step inside the failed 1jxy run up to its stop), the stream
entries (`stream_preproc`, `stream_chains`, `mini_stream_1jxy`, `mini_stream_6vxx`), the truth suites
(`multistate_workflows`, `solve3D_addon`) and the single-particle, species and atoms entries. The nightly is green
but for the base failure the plan names, which is reported as such and not counted against the phase (section 12).

**Exit criteria of phase 3.** The first exit point, the distributed cases byte-identical to phase 1 at their layout,
holds for the `solve3D` case in every map, binary and segment on the first build and on the final one, and for the two
`solve2D` cases on every stage whose command line sets the band (iterations 1 to 22); from stage 6 on they follow the
shared-memory trajectory phase 2 validated, parting at the same iteration for the same cause (the retired table
worker's band without an explicit limit), and the control with the band and the window pinned is byte-identical. The
run's ruling is phase 2's: the stop rule for a structural phase does not apply, because the cause is a defect of the
replaced path, removed by construction within this phase's files; the three toolbox defects phase 2 left on the
distributed path are gone with the path. The second exit point, level 4 within the frozen tolerance with the truth
gates met, holds: the fixed-input accumulators give phase 0's numbers to the digit (at most 8.58e-5 and 7.79e-7
against 1.72e-4 and 1.56e-6), the probabilistic runs across the layouts are the same search wherever the base's
accumulators and classic draws let them be and stochastic trajectories where they do not, and the truth suites pass
with their floors unchanged (`multistate_workflows` 129 of 129, `solve3D_addon` 46 of 46, `cont_refine3D_1jxy` 135 of
135, `flex_pca_blobs` 94 of 94, the two stream entries, the simulated workflows' two base failures at their points).
The third exit point, the nightly green, holds (18 of 19 executions, the nineteenth the base failure of
`simulated_workflow_1jxy` at its recorded point). Level 2 holds on the three builds in Release at eight threads and in
Debug at four; the fast gate 15 of 15 and the four unit suites on both trees after every round. The real-data readouts
at the distributed layout are within the seeds' spread (the bgal class resolutions over eight seeds overlapping the
base's and phase 1's in every statistic; the trpv1 shells phase 1's to the last digit at every iteration). One
scheduling phase per iteration on both masters; the readouts merged in part order; no table, assignment or
insertion-pose file on any path. The table workers, the masters' step, their readers, writers, file names and private
programs are unreachable from production and go in phase 4 with the measurements. Phase 1's note on the part files'
seed-shift section is closed with the part files. Phase 2's two side effects (the reseed counter, the classic order
not drawn under `prob`) and phase 1's `FFTW_ESTIMATE` rule carry forward to phase 4's measurements. Nothing is
committed; the working tree holds phases 0 to 3 as one diff. Phase 4 (retirement, measurements, documents) can start
on this tree.

### Phase 4: retirement, measurements, documents (2026-11-05, the Dell)

Done on 2026-11-05 by the autonomous run `prob_pipeline` on the Dell, on the tree phase 3 left (base HEAD
`22a643642`; the run root `/data2/agent_runs/prob_pipeline`, linked as `~/agent_runs/prob_pipeline`, and paths
below that are not absolute are relative to it). The tree built in `repo/build` (Release) and
`scratch/keep/build_debug` (Debug); the phase-3 sources rebuilt in `scratch/build_p3` for the identity comparison
(the driver had deleted the staged phase-3 binaries at phase 3's end, as its cleanup rule says); the base build
`scratch/keep/build_base` for the measurements. The journal is `state/journal_phase_4.md`.

**What changed, in one paragraph.** The old two-phase path is gone. The probability tables keep the fills, the
selection rules, the extras and cones, the temperatures, the n_eff readout and the fill's cost record, and lost
their global constructors, their part-file writers and readers, their whole-table assignment loops and the
storage of the master's table; the `prob_align` step module and the table workers' commander module are deleted
with the private programs `prob_tab`, `prob_tab_neigh`, `prob_tab2D`, `prob_align`, `prob_align_neigh` and
`prob_align2D`, their entries in the private driver and its API module, the file-name constants of the part,
assignment and insertion-pose files, the deleters of those files at every site that called them, and the
`PROB_TAB*` and `PROB_ALIGN*` bench families, whose counters now print in the matchers' own bench records. The
fill harness runs the pipeline in parts and the table testers test the selection rules on the pipeline's batch.
The net change is 3,313 source lines fewer (the plan's estimate: a third of the 9,000 lines of section 2.1). The
documents of section 14 describe the pipeline; the generated indexes are regenerated. The measurements of
section 12 and the nightly are below.

**The retirement, in the order of the notes** (one build after the tables, the modules, the entries, the
constants and their users together, because each step's only callers were the next step's: the master-side
readers' caller was the `prob_align` step module, the writers' the table workers' module, the programs' the
driver and the API module, so the tree compiles only with all of them gone; a second build after the bench
counters moved; the fast gate 15 of 15 and the four unit suites after each). *The master-side table objects and
readers.* `simple_eul_prob_tab.f90` lost `new` (the global constructor), `new_worker`, `new_compact_global`,
`new_assignment`, `read_tabs_to_glob` and `read_dense_part`, `ref_assign` (the whole-table loop) and
`warn_invalid_and_flat`, `read_assignment`, `read_insert_table`, `dense_candidate` (the tester's view of a row),
the master's storage (`dists`, `inpls`, `refined`, `refined_offsets`), the part-file state, the format and status
constants, and the master branches of `get_row`, `max_refined`, `assign_rows` and `close_assignment`; the
pipeline flag went with them, there being one constructor left (`new_pipeline`). `simple_eul_prob_tab_neigh.f90`
lost `new_neigh`, `new_neigh_worker`, `new_neigh_global`, `ref_assign_neigh`, `read_tabs_to_glob` and
`read_sparse_tab_to_glob`, the candidate counts and chunk counters of the stream. `simple_eul_prob_tab2D.f90`
lost `new`, `new_worker`, `new_assignment`, `ref_assign`, `read_tabs_to_glob`, `read_part_head`,
`read_part_cands`, `part_read_error`, `read_assignment`, the master's rows and refined lists, the part type and
the formats, the assignment profile of the retired `PROB_ALIGN2D` record. *The part-file writers.* The three
tables lost `begin_write`, `write_tab`, `write_rows`, `write_dense_head`, `write_stream_head`,
`flush_candidate_buffers`, `flush_dense_rows`, `write_assignment` and `write_insert_table`; the utilities lost the
seed-shift table's writer and reader, the particle-index lookup of the readers, the insertion table's file
writer, reader and subset and the whole-table layout of the insertion poses (the batch layout, `setup_range`,
stays). *The programs and their commanders.* `simple_commanders_prob.f90` (the six commanders, the table workers'
frames, `write_prob_tab_bench`, `write_empty_prob_tab`) and `simple_prob_align_strategy.f90` (`prob_align_exec`,
`prob_align_neigh_exec`, `prob_align2D_exec`, their `run_*` entries, `write_prob_align_bench`) are deleted; the 2D
matcher's `set_b_p_ptrs2D`, which only the 2D step called, with them. *The API and driver entries.*
`simple_private_exec_api.f90` no longer imports the commanders and `simple_private_exec_driver.f90` no longer
dispatches the six programs; no exec router, UI table or generated list named them. *The file-name constants.*
`DIST_FBODY`, `ASSIGNMENT_FBODY` and `INSERT_FBODY` are gone from `simple_defs_fname.f90`, and with them the
deleter `del_prob_align_files` and its calls after each refine3D stage (`simple_solve3D_utils.f90`), in the
missing-update passes of the multi-state controller and of `solve3D`, in frequency marching
(`simple_commanders_refine3D.f90`), the deleter of the external-reference pose initialisation, the 2D master's
iteration cleanup (which still deletes the readout records) and `solve2D`'s final cleanup; the `cont_refine3D_1jxy`
test's check that a Cartesian pass leaves no table now looks for a part's readout record. *The bench families.*
The fill's cost record (`prob_fill_profile`: the particles filled, the references scored and refined, the shift
searches, the coarse maps, the optimiser's evaluations, the shifted particle memos; the thread-seconds of the seed
estimation, the score loop and the refinement, the fill region's wall time and the load balance) is printed by
`write_fill_profile_bench` (the utilities) in each part's `REFINE3D_BENCH_PART` and `REFINE2D_BENCH_PART` record
as the `match3D fill ...` and `match2D fill ...` lines, where the `PROB_TAB*_BENCH_PART` records carried it; the
strategies' records already timed the master's draw where they timed the pre-step. The generated Fortran indexes
(`perl scripts/gen_fortran_indexes.pl --root src --out doc/code_overview/fortran-indexes`) no longer list the two
modules; `scripts/check_test_registry.py .` and `python3 scripts/check_descr.py` pass.

**The harness and the testers.** The fill harness's reference route was the retired code: per-part worker tables
written to part files, merged into the master's table and assigned at once. Its pipeline route takes over the
part and perturbation dimensions (D1 of the journal): `run_fill` now runs one pipeline table per part over the
part's index range, as a distributed part does, with the iteration's model computed once from the unperturbed
field and set on every part (as each part computes it from the field the master wrote), gathers the rows,
temperatures, seeds, assignment and extras by particle identity, and perturbs the earlier parts' committed rows
before the later parts' fills as before; the configurations are the six of the reference route (one and eight
threads, one batch and a third, one and three parts, the earlier parts perturbed), the reference one thread,
one batch, one part. The harness's `polar_fill_identity.tsv` loses its `pipeline` column. The 3D table tester
lost the dense part files' round trip and the parallel part read (50 assertions), and its selection tests run on
the pipeline's batch (the rows set on `batch_rows`, the candidates appended to the thread buffers, then
`open_assignment`, `assign_batch`, `close_assignment`); the insertion test keeps the in-memory table and drops
the file round trip, the subset and the whole-table layout; the multi-state extras test builds its rows in memory.
The 2D tester's stream tests became batch tests through the thread buffers, its parallel-read test went, and
its sub-suite is `2D probability table` (the `I/O` of its name retired with the part files). The harness's help
and the sub-suite list of `unit_pftc_align2D3D` name the pipeline.

**Decisions of the phase.** *The order of retirement and the builds* (D3). The notes retire in six steps with a
build after each; each step's only callers were the next step's (the readers' the `prob_align` module, the
writers' the table workers' module, the programs' the driver and the API module), so the tree compiles only with
the first four steps gone together, and the constants' users and the bench counters are two more edits that
compile alone. Two builds were made: the retirement with the constants, and the counters moved into the matchers'
records; each passed the fast gate and the four unit suites. *The harness* (D1) and *the testers* (D2) as above:
the reference route retired with the code it ran, the pipeline route in parts holding level 2's schedules; the
selection tests ported to the pipeline's batch, the part-file tests retired. *The phase-3 binaries* (D4). The
staged phase-3 install under `scratch/` lost its binaries to the driver's end-of-phase cleanup (files over 1 MB
outside build trees and `scratch/keep/`); the phase-3 sources were rebuilt from `HEAD` and the driver's
`review/full_after_3.diff` in `scratch/build_p3`, and level 3 compares the two builds' runs of the same cases
(a lesson for the run's records: binaries a later phase needs go under `scratch/keep/` or stay in a build tree).
*The scripts outside the file table.* `scripts/plot_refine3d_bench.py` reads the retired bench families and the
strategies' old pre-step label; it was already out of step with the labels phase 3 changed, is not in the file
table, and is left for the maintainer (the report's follow-ups). `doc/code_overview/code_base_map.md`, generated
by another script, still names the two modules; the plan regenerates the Fortran indexes only.

**Level 1 (draw identity)** holds unchanged: nothing in the generator moved, and `unit_numerics` passes in Release
and Debug (1,480 assertions) on both builds of the phase.

**Level 2 (frozen-input identity)** holds on the pipeline route's six configurations, the only route the harness runs
now: 48 rows (four routes, two state counts, the configurations at one and eight threads, one batch and a third, one
and three parts, the earlier parts perturbed), every difference count zero and the run's own verdict (`the fill
harness gives identical rows, selections and extras in every configuration`) in Release at eight threads on the final
build (`scratch/keep/phase4/harness/a_nthr8.tsv`) and in Debug at four threads
(`scratch/keep/phase4/harness/a_debug_nthr4.tsv`, 48 rows, every count zero). The phase changes no OpenMP directive
and no code inside a parallel region (deletions, a bench writer and the harness's driving loop), so the 16-thread rule
has nothing to apply to. The fast gate passes 15 of 15 on both trees after both builds and the four unit suites pass
on both: `unit_pftc_align2D3D` (1,466 assertions; 1,516 before, the part-file tests' 50 retired with the code),
`unit_numerics` (1,480), `unit_ori` (1,047) and `unit_project` (2,979).

**Level 3 (structural identity with phase 3 at a fixed configuration).** The staged phase-3 install had lost its
binaries to the driver's end-of-phase cleanup, so the phase-3 sources were rebuilt (`scratch/build_p3`, from `HEAD`
and `review/full_after_3.diff`) and `scratch/run_phase4_level3.sh` ran the same cases from both builds with the same
seed, the phase-3 set after the final one, and compared the pairs with `compare_dirs.py`, the registry's
`test=identity` on the entries and `compare_projects.sh` on the projects (`scratch/keep/phase4/level3_compare/`).
Every declared output is byte-identical. The three distributed cases at two parts and one thread per part: the
`solve3D` case 41 mrc and 4 bin the same and the project the same in every segment but `projinfo`'s run id; the two
`solve2D` cases 82 and 85 mrc, 27 and 28 star and 3 bin the same, every segment the same. The shared-memory entries:
`simulated_workflow_1jxy` at one and at eight threads 126 mrc, 21 mrcs and the eight projects the same (the identity
test: 139 files the same, 1 differing, the time-stamped sigma2 provenance; the base failure at solve3D stage 3,
iteration 41, on both builds), `simulated_workflow_6vxx` at both thread counts 174 mrc and 30 mrcs the same (184 the
same, 0 differing; the base's final-volume FAIL on both), `mini_stream_1jxy` 107 mrc, 7 mrcs and 27 star the same; the
trpv1 `refine3D` case at eight threads 9 mrc and 4 bin the same, the trpv1 `refine2D` case 6 mrc, 3 bin and 2 star the
same. What differs carries a time or a process id: `metadata.json`, the shift and FSC pdf plots, the manifest, 8 bytes
of the sigma2 state's header, the provenance files, the pid-named jpgs. The masters' readout lines are identical to
the character in every case (none of 110, 113, 828, 15 and 10 lines differ), the cone lines included. The retirement
changed no number.

**The measurements of section 12** (`scratch/run_phase4_measure.sh`; the machine otherwise quiet, the load before
each run in the timelines; the base build `scratch/keep/build_base` and the final install, three runs each of the
trpv1 cases of section 2.5 at four parts of six threads: `refine3D refine=prob objfun=euclid box_crop=160
shift_sigma=0` for three iterations, `refine2D refine=prob objfun=euclid` with 100 classes for two iterations at
crop 120 and at box 256; the master under `/usr/bin/time -v`, the node's memory in use and every worker's resident
set sampled every second from `/proc` by `scratch/sample_mem.sh`; the records, timelines and the summary by
`scratch/summarise_phase4_measure.py` under `scratch/keep/phase4/measure/`). Medians over the three runs; the
parts' numbers are means over the four parts.

| trpv1 `refine3D`, 35,645 particles, 2,500 directions, crop 160 | base | final |
| --- | --- | --- |
| wall time of the run, three iterations (s) | 1,015 (990, 1,015, 1,026) | 1,001 (989, 1,001, 1,012) |
| the iteration, 1 / 2 / 3 (s) | 342 / 353 / 312 | 340 / 349 / 309 |
| of which the probabilistic pre-step (base) or the parts' pipeline (final), 1 / 2 / 3 (s) | 315 / 326 / 286 | 335 / 345 / 304 |
| the matcher's scheduling phase (base), 1 / 2 / 3 (s) | 23 / 23 / 22 | (one phase) |
| the master's sample draw (s) | inside the pre-step | 0.02 |
| the part's work per iteration, 1 / 2 / 3 (s): table job + search job (base), the pipeline job (final) | 297 + 22 / 304 + 22 / 269 + 22 | 318 / 324 / 288 |
| the fill, 1 / 2 / 3 (s) | 293 / 300 / 265 | 295 / 299 / 263 |
| particle preparation + reading + references, per job (s) | 3.0 + 0.6 + 0.4 (table) and 3.0 + 0.5 + 0.3 (search) | 3.3 + 0.6 + 0.3 |
| the search's polish, the partial reconstruction (s) | 9.1, 9.0 | 10, 8 |
| the assignment (s) | 0.05 on the master | 0.3 per part |
| references refined per particle, 1 / 2 / 3 | 112 / 120 / 99 | 112 / 119 / 100 |
| shift searches, objective and gradient evaluations, memos (iteration 1) | 1.00e6, 1.29e9, 1.34e7, 8,770 | 1.00e6, 1.29e9, 1.34e7, 8,770 |
| the fill's thread-seconds: seed, score loop, refinement (iteration 1) | 31, 570, 1,130 | 31, 567, 1,140 |
| the fill region's load balance | 0.984 | 0.984 |
| master maximum resident set (GB) | 1.26 (1.22 to 1.31) | 1.10 |
| worker peak resident set (GB) | 2.43 the search worker, 1.56 the table worker | 2.42 (the alignment context 1.38 GB within its 6.7 GB budget) |
| node peak memory in use (GB), over a 6.2 GB baseline | 15.8 to 16.1 | 15.9 to 16.0 |

| trpv1 `refine2D`, 100 classes, two iterations | base, crop 120 | final, crop 120 | base, box 256 | final, box 256 |
| --- | --- | --- | --- | --- |
| wall time of the run (s) | 41.6 (41.4 to 42.0) | 34.6 (34.3 to 34.8) | 96.0 (96.0 to 97.3) | 81.5 (81.3 to 82.0) |
| the iteration, 1 / 2 (s) | 17.8 / 20.7 | 15.6 / 16.1 | 44.3 / 49.1 | 39.6 / 39.0 |
| the pre-step + the matcher phase (base), the parts' pipeline (final), 1 / 2 (s) | 12.3 + 5.3 / 15.3 + 5.2 | 15.3 / 15.8 | 27.5 + 16.2 / 32.6 + 15.9 | 39.0 / 38.5 |
| the part's work per iteration, 1 / 2 (s): table job + search job, or the pipeline job | 11.1 + 5.0 / 14.1 + 4.9 | 14.0 / 14.4 | 25.9 + 15.5 / 31.1 + 15.3 | 37.1 / 36.4 |
| the fill, 1 / 2 (s) | 8.5 / 11.3 | 9.1 / 9.1 | 19.5 / 24.8 | 21.0 / 20.2 |
| particle preparation per job (s) | 2.0 (table) and 1.85 (search) | 2.0 | 5.6 (table) and 5.1 (search) | 5.5 |
| the search, the class sums (s) | 0.6, 1.7 | 0.6, 1.8 | 1.2, 8.2 | 1.1, 8.7 |
| shift searches, objective evaluations, coarse maps (iteration 1) | 8.91e4, 1.10e8, 1.8e5 | 8.91e4, 1.11e8, 1.9e5 | 8.91e4, 2.39e8, 1.8e5 | 8.91e4, 2.40e8, 1.9e5 |
| master, worker, node peak (GB) | 0.32, 0.48 (0.35 the table worker), 8.2 | 0.32, 0.49, 8.3 | 0.46, 0.83 (0.62), 9.7 | 0.44, 0.84, 9.7 |

What the numbers say. *Wall time.* In 3D the fill is 93 per cent of the part's job on both builds and is
unchanged (295 against 293 s at the first iteration, the same shift searches, evaluations and memos, the same
thread-seconds within 1 per cent); the pipeline saves the second scheduling phase and the second preparation, read
and reference load of each part (about 4 s a part, 1.3 per cent, as section 2.5 estimated), and the run is 1.4 per
cent faster (1,001 against 1,015 s). In 2D the duplicated preparation was 8 per cent (crop 120) and the full-box
preparation 13 per cent (box 256) of the part's job, and the runs are 17 and 15 per cent faster (34.6 against 41.6 s,
81.5 against 96.0 s); the 2D fill's own wall time is 7 to 8 per cent higher in the pipeline at the first iteration
(9.1 against 8.5 s, 21.0 against 19.5 s, with the same shift searches and 1 to 9 per cent more objective evaluations
and refinement thread-seconds), which the batch structure explains best (the pipeline fills the matcher's batches,
smaller than the table worker's, and a batch builds its search objects), and which the preparation saved outweighs
by three to four times. The second iterations are not comparable across builds: the keyed draws give another
trajectory, and the base's second iteration refined more (its coarse maps 5.3e5 against 2.4e5). *Memory.* The
master's peak falls by 13 per cent in 3D (1.26 to 1.10 GB) with the global table, and stays in 2D (0.32, 0.44 GB);
the worker's peak is the search worker's (2.42 GB in 3D: the particle buffer and the reconstruction, the alignment
context at 1.38 GB within its budget; 0.49 and 0.84 GB in 2D), where the table worker's 1.56 GB (0.35, 0.62 GB in
2D) is gone from the node; the node's peak is unchanged (16 GB in 3D, 8.3 and 9.7 GB in 2D, over a 6.2 GB baseline),
being the four search workers' peaks beside the master's; at the 3D barrier the master assembles while the parts
have exited. At 750,000 particles the master's table would have been 15 GB per state (section 2.5): that, not the
wall time, is what the pipeline removes. *Preparations per particle.* One compatible preparation per particle per
iteration on the final build (3.3 s a part in 3D, 2.0 and 5.5 s in 2D) against two on the base; the reference loads
once a part (0.3 s) against twice; no table bytes written (199 MB a part per iteration before); the map and
optimiser work of the fill unchanged.

**The documents of section 14.** `doc/policies/3D/refine3D_policy.md`: the workflow's steps (the search with the
pipeline inside), the ownership (the matcher owns the pipeline, the tables the fills and selection rules, no
`prob_align` module), the master's project (its sample draw instead of the step, the parts' readout records),
section 8 (the pipeline in the matcher, the keyed draws against the thread-count paragraph, the dense and
neighbourhood assignments per batch, the n_eff rows by rank in the global sample merged in part order, the cone
opened by each part and reported by the master, the extras kept in the part's insertion table, the part files'
subsection replaced by the rows' lifetime in the batch, the fill's counters in the matcher's record), the bench
file table without the retired families, the refactor rules with the readout records in the place of the
assignment maps. `doc/policies/2D/solve2D_policy.md`: the steps, the ownership (the master's draw and one
scheduling phase, the matcher's pipeline), the table contract (batch-bounded, no partition file, no global object),
the part files' retirement, the assignment per batch, the empty part, the one preparation per batch, the
project-segment table (no write-and-reread in shared memory; the worker's readout record), the artifacts, the
checklist and the rules. `doc/policies/importance_sampling_fractional_update_policy.md`: the master's keyed draw
in the place of `prob_align_exec`, the parts' reproduction, section 6 rewritten as the pipeline's sampling
contract (the four steps and the keyed inner draws), the invariants and the checklist.
`doc/policies/layering_and_project_ownership_policy.md`: the three master steps and the exception that ended for
the probabilistic step (section 2), the step table and the worker frame held once in the matchers (section 3),
the readout record among the part files (5), the keyed draws under the reseed rule (6), the red flag without its
retired example (7), section 8 without the three entries the pipeline closed and with the write counts recounted,
the audit without the worker frame (9). `doc/policies/particle_read_path_policy.md`: the hand-over retired, the
budget shared with the 3D pipeline's context. `doc/policies/test_environment_policy.md`: identity decided at a
fixed thread count on both builds (one, and eight where a run checks it), the keyed draws against the
thread-reproducibility statements, the sub-suite's name. `doc/algorithms/sampling_and_fractional_updates.md`,
`refine2d_class_averaging.md` and `solve2d.md`: the rows filled per batch and drawn from, the keyed sample, the
implementation pointers (and `solve2d.md`'s stale untempered-draw and balanced-assignment paragraphs brought to
the tempered draw of 2026-10-23/24). `speed_architecture_proposals.md`: the header and the review marked
implemented, the RNG section resolved, a status entry of 2026-11-05. Four other live documents that named the
retired path were added to the file table and corrected: `doc/policies/2D/particle_cache_policy.md` (one read per
iteration), `doc/policies/3D/solve3D_addon_policy.md` (the partitioning and the empty part),
`doc/algorithms/refine3d.md` (the candidate preparation inside the search) and the implementation notes' README
(the pipeline row's status). The skill text is proposed in the report, not edited (section 14).

**`test=nightly`** on the final install (`scratch/phase4_runs/nightly`, the report kept under
`scratch/keep/phase4/nightly/`): 19 executions in 2,221 s on the quiet machine, 18 passed with every metric at or past
its floor (1,078 metric rows, none below), one failed: `simulated_workflow_1jxy` stops at solve3D stage 3, iteration
41, in `new_shiftmap` ("the map's last shell must lie in the band"), the base failure section 12 records and every
phase's level 3 reproduced at the same point on both builds. `simulated_workflow_6vxx` passes its gate with exit 0 and
prints the recorded `FAIL: final-volume correlation is below the required minimum` of the base. The passing entries
include the distributed ones (`flex_pca_blobs`, `cont_refine3D_1jxy` with its no-table check on the readout records),
the stream entries (`stream_preproc`, `stream_chains`, `mini_stream_1jxy`, `mini_stream_6vxx`), the truth suites
(`multistate_workflows`, `solve3D_addon`) and the single-particle, species and atoms entries. The nightly is green but
for the base failure the plan names, reported as such and not counted against the phase (section 12).

**The Mac check for the maintainer.** The run tested the Linux gfortran builds only (GCC 15.2.1, Release
`-O3 -march=native` and Debug `-O0 -fcheck=all`); the generator's vectors of phase 0 are what the Mac build has to
reproduce before the diff is accepted (section 9). The check is `simple_test_exec test=unit_numerics` on the Mac
build, whose keyed-generator suite pins the ten vectors below (the same table as phase 0's entry; the Python
reference `scratch/keep/rnd/rnd_keyed_ref.py` and the three vector files `vectors_{python,release,debug}.txt`
are kept beside it). A vector that differs on the Mac is a difference in the Mac compiler's integer or
floating-point arithmetic of the hash or of the 53-bit mapping, not a difference in SIMPLE's draws, and would mean
that the keyed draws of a Mac run are not those of a Linux run of the same seed.

| seed | version | scope | iteration | particle | family | reference | ordinal | 64-bit word (signed) | unit variate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | -4131892048949355934 | 7.76009683202676270E-01 |
| 20261103 | 1 | 17 | 1 | 1 | 2 | 1 | 0 | -3402735607585271042 | 8.15537333093113270E-01 |
| 20261103 | 1 | 17 | 1 | 2 | 2 | 1 | 0 | 4682640351440764835 | 2.53846442100018055E-01 |
| 20261103 | 1 | 17 | 1 | 1 | 2 | 2 | 0 | -2904353384885242822 | 8.42554687522089507E-01 |
| 20261103 | 1 | 17 | 2 | 1 | 2 | 1 | 0 | 6790283504847523914 | 3.68102006387408531E-01 |
| 20261103 | 1 | 17 | 1 | 750000 | 5 | 0 | 0 | -4401584271307322743 | 7.61389638533528701E-01 |
| -1 | 1 | 33 | 3 | 12345 | 3 | 0 | 7 | -2048942294402738840 | 8.88926615655555752E-01 |
| 4611686018427400249 | 1 | 18 | 1 | 1 | 4 | 20000 | 2147483647 | 3588558501286480495 | 1.94536146159306367E-01 |
| 123456789012345 | 1 | 17 | 25 | 35645 | 1 | 0 | 35645 | -3713091948123483349 | 7.98712882160304227E-01 |
| 20261103 | 1 | 17 | 1 | 1 | 6 | 0 | 1 | -5700306898616924037 | 6.90985743834271138E-01 |

**Exit criteria of phase 4.** The first exit point, the nightly green, holds (18 of 19 executions, the nineteenth the
base failure of `simulated_workflow_1jxy` at its recorded point). The second, the measurements and the Mac vectors to
check in Progress, holds: the wall times per iteration and per phase of the 2D and 3D trpv1 cases against the base
build, the preparations per particle from the bench counters, and the peak resident sets of a worker, of the master
and of the node are above, and the Mac check with phase 0's ten vectors is the section before this one. The third, the
documents updated, holds: the ten documents of section 14, the layering policy's section 8, the status of
`speed_architecture_proposals.md`, four other live documents and the regenerated Fortran indexes. The programs,
commanders, files, objects and bench families of section 6 are retired; the fast gate, the registry check and the
description check pass on the final tree; level 2 holds on the pipeline route's configurations in Release and Debug;
level 3 against phase 3 holds at every fixed configuration (every declared output of the three distributed cases, the
two simulated workflows at one and eight threads, the mini stream and the two trpv1 shared-memory cases byte-identical
between the phase-3 and the final build, the masters' readout lines identical to the character). The skill text is
proposed in the report. Nothing is committed; the working tree holds phases 0 to 4 as one diff. Phase 5 (the report)
can start on this tree.
