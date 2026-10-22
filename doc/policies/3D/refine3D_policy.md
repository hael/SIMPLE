# Refine3D Policy

This document records the current policy for the base `refine3D` command. It
does not describe the staged `solve3D` controller, the class-average
`solve3D_cavgs` initializer, or the automated wrapper `refine3D_auto`. The one
exception is the per-state M-estimation of `refine3D_auto state=X` (section 8),
because its weighted reconstruction runs in the base matcher and assembly.

Related workflow policies:

- [solve3D_policy.md](solve3D_policy.md)
- [solve3D_cavgs_policy.md](solve3D_cavgs_policy.md)
- [class-average rejection policy](../sieving_and_rejection/model_cavgs_rejection.md)
- [refine3D_auto_policy.md](refine3D_auto_policy.md)
- [refine3D states policy](../heterogeneity/refine3D_states_policy.md)
- [classify3D_refs_policy.md](classify3D_refs_policy.md)
- [automasking_policy.md](automasking_policy.md)
- [nonuniform filtering policy](../NU/nonuniform_filtering_policy.md)
- [reconstruct3D_pcg_policy.md](reconstruct3D_pcg_policy.md)
- [particle cache policy](../2D/particle_cache_policy.md)
- [separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md](separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md)

## 1. Scope

`refine3D` is an iterative projection-matching workflow over an existing
project and one or more starting 3D references. Its durable contract is:

1. initialize project state and execution mode
2. materialize the iteration reprojection model from current Cartesian volumes
3. optionally run probabilistic pre-alignment
4. run particle-domain search and hard assignment
5. write partition-local Cartesian reconstruction inputs
6. run volume assembly when volume reconstruction is enabled
7. persist project state for the next iteration or caller

The supported reconstruction path is Cartesian. Cartesian matching still uses
polar Fourier central sections internally, but those sections are generated
from the current Cartesian volumes. `refine3D` does not own a separate
non-Cartesian reconstruction branch.

Particle-domain work stays in the probabilistic pre-step, matcher preparation,
search strategies, pose/state/shift updates, sigma updates during Euclidean
matching, and partial reconstruction writing. Volume-domain work stays in
`volassemble` and the volume postprocessing helpers it calls.

Shared-memory and distributed `refine3D` must preserve the same scientific
workflow and artifact contracts. Only toolbox lifetime, process launch, and
scheduling differ. Three rules belong to one strategy: only the distributed
strategy accepts convergence from the third iteration of the call on
(`which_iter >= startit + 2`) and schedules the final `combine_eo` iteration
(section 10); only the shared-memory run finalizer removes the `os_out`
entries of empty states (section 11.2).

## 2. Entry Points and Execution Mode

The public command is `refine3D`, registered in `simple_ui_refine3D.f90` and
routed by `simple_exec_refine3D.f90`.

`simple_commanders_refine3D.f90` performs the thin command-level setup:

- validate multi-volume input against `nstates`
- set local defaults such as `mkdir`, `cenlp`, `oritype`, and `prg`
- select a strategy through `create_refine3D_strategy`
- run the iteration loop
- delegate finalization and cleanup to the selected strategy

Strategy selection is command-line shaped:

- `nparts` without `part` selects the distributed master strategy
- otherwise the shared-memory strategy is used
- distributed workers run the partition command path with `part` and `outfile`
- the distributed master splits the rows into contiguous partitions that
  balance the particles with state > 0 (`qsys_env%new(..., l_active)`,
  `split_nobjs_active`); with every row active this is the even split. The
  master merges the partition documents by their own ranges (`merge_algndocs`),
  so no consumer assumes a particular split
- a worker whose partition samples nothing (fewer active particles than
  partitions, or none of a class-balanced or probabilistic sample drawn over
  the whole project) completes an empty transaction: its samplers accept an
  empty sample only when `part` is set (`allow_empty`), and it writes its
  unchanged sigma2 slice, its range's orientations, zero PCG accumulators and
  `JOB_FINISHED`

`maxits` is the number of iterations to run in the current invocation.
`which_iter` starts at `startit`, and `extr_iter` follows the same per-call
counter unless supplied by a caller.

Shared-memory `refine3D` currently rejects `continue=yes`. Distributed
`continue=yes` resumes from project-carried output volumes, FSC files, the
trailing chains and the sigma mode. It does not depend on the previous
partitioning: the partial reconstructions of the previous run are transient
(removed before every round) and are never read on continuation, so no
partial-file count is checked; a state without particles in a part writes no
partial, which is a zero contribution.

## 3. Ownership

`simple_commanders_refine3D.f90` owns the base command entry point and the
shared iteration loop. It should stay thin. The module also hosts the commanders
of `refine3D_auto`, `refine3D_states`, `refine3D_states_cavgs`,
`classify3D_refs` and `bootstrap_rec3D`.

`simple_refine3D_strategy.f90` owns:

- shared-memory versus distributed orchestration
- scheduler interaction
- iteration counters and convergence checks
- probabilistic pre-step dispatch
- matcher dispatch
- volume assembly dispatch
- per-iteration reprojection-model materialization
- run finalization

It must not absorb numerical volume postprocessing or assembled-reference
construction.

`simple_strategy3D_matcher.f90` owns `refine3D_exec`: consuming the
driver-generated reprojection model, selecting and running the search strategy,
updating orientations/states/shifts, updating Euclidean sigma estimates during
search, and writing partition-local Cartesian partial reconstructions.

`simple_prob_align_strategy.f90` (`prob_align_exec`, `prob_align_neigh_exec`)
plus `simple_eul_prob_tab*.f90` own probabilistic pre-alignment and the
assignment artifacts consumed by the matcher.

`simple_volassemble_strategy.f90` owns volume assembly (`volassemble_exec`):
reducing partials, restoring even/odd volumes, calculating FSCs, running
volume postprocessing, writing state volumes, and updating the particles'
resolution fields. It returns the names of the volumes each state ends with
and leaves `params%vols` alone.

`simple_strategy_sigma2.f90` owns the canonical sigma2 commit
(`commit_sigma2_exec`).

These are strategy routines on a caller's `(params, build)` with explicit
arguments for what differs from the caller's own settings (the thread count,
`which_iter`, a write flag). The commanders `prob_align`, `prob_align_neigh`,
`volassemble` and `commit_sigma2` are thin wrappers: their `execute(cline)`
builds parameters and a builder from the command line and calls the routine.
Shared-memory `refine3D` calls the `prob_align`, `prob_align_neigh` and
`commit_sigma2` commanders; it and `reconstruct3D` call `run_volassemble`, the
routine's own entry that builds from the command line. The other command-line
callers use the commanders. A routine that runs on the caller's parameters does
not run `params%new`: it sets the thread count it needs and the caller restores
its own afterwards, and it draws from the random generator without a reseed
(ruling of 2026-10-20).

`simple_vol_pproc_policy.f90` owns one check: whether a state's density-mask
file exists and matches the box and sampling of the volumes it would mask
(`state_mask_is_compatible`, shared by volume assembly, postprocessing and the
solve3D final reconstruction). The NU setup API owns its spherical `mskdiam`
support invariant.

## 4. Project and Builder State

The project is durable workflow state. The `builder` owns derived execution
state for one parsed command line and must not be treated as durable policy
state.

Shared-memory `refine3D` rebuilds its strategy toolbox each iteration after
the iteration command line has settled. Persist any initial project mutations,
such as random orientations, sampling-counter cleanup, or state initialization,
before the next rebuild reads the project.

Distributed `refine3D` keeps prototype command lines for the workers, the start
reconstruction, postprocessing and the thread count of the sigma commit; the
probability and assembly command lines are made each iteration from the
master's own. These command lines are execution templates, not independent
workflow policies.

The distributed master is the one owner of its project. It materialises the
reprojection model on its own builder
(`materialize_reprojection_model_master`): the general toolbox once, at the
first iteration, which also reads the segments the model uses (the
non-data segments, `out`, `optics` and the particle field, then
`update_projinfo`, so that writes go to the run's own project and not to a
project the run was copied from), because the start reconstruction may
have written them; on later iterations the builder's volumes again
(`new_vols`), the search toolbox, the band and the model. `prob_align` (or
`prob_align_neigh`), volume assembly (gridding) and the sigma commit run on
the master's `(params, build)`; the PCG master already did. The particle
field therefore holds this iteration's resolution fields (`cfar`,
`lp_est`) in memory, and the matching low-pass refresh reads the master's
own project.

The master writes the particle field once before the workers start (with
the alignment's sample, or as the model left it), once more after the
merge only when a polish pass follows (its workers read the merged field),
and once at the end of the iteration; `merge_algndocs(..., write=.false.)`
merges the workers' documents in memory, and gridding assembly does not write
the project on this path. The PCG master (`rec_backend=pcg`) writes the
particle field once, at the end of its call. The `out` segment is written once,
before postprocessing, which reads it from the file. The final full write
is unchanged. Workers read only the file; the master reads only what the
workers write.

Fresh-start checks that decide whether to consume project-carried matching
metadata must use the stage interval: `which_iter <= startit` with
`continue != yes` is fresh, even when `startit` is greater than one. A fresh
start still consumes the project's matching low-pass when the `_nu_filt`
references it was derived from are on disk (the `refine3D_auto` startup).

Before a matcher iteration writes partition-local Cartesian partials, stale
partial reconstruction files for the active states and partitions must be
removed.

## 5. Startup and State Initialization

If multiple starting volumes are supplied, `nstates` must be defined and
greater than one. Multi-state startup requires either all `vol1..volN` inputs
or none.

When orientations are absent, `refine3D` randomizes 3D orientations. When
projection indices are absent, distributed initialization sets them from the
Euler sampling.

For multi-state startup (the distributed master always, shared memory when the
field holds at most one state bin):

- existing state assignments must match the requested `nstates`
- fresh state labels are randomized uniformly, whether complete starting
  volumes are supplied or the project references are used
- probabilistic multi-state startup requires each state to exceed the minimum
  population guard

Even/odd partitioning is required for consistent half-map reconstruction. The
partition belongs to the project, not to a field
(`ensure_field_halves`, `src/main/project/simple_eo_partition.f90`). The
master makes it once, before any worker runs, and workers stop rather than make
one. The same labels go to `ptcl2D` and `ptcl3D`, and existing labels never
change. A stack row that maps to a micrograph is a group: all its particles
share one half, so overlapping neighbouring picks never sit in opposite halves
and a reconstruction opens a stack once per state rather than once per half.
Groups go, largest first, to the half with fewer active particles. Their order
comes from the partition's own random stream, which leaves the particle sample
unchanged. A field with stack rows but no micrograph rows takes the per-index
route (`partition_eo`), as do class-average fields. The stack rows record which
route labelled them, so a missing, partial or inconsistent partition is told
apart, and an inconsistent one stops the run.

The offline `select_particles eo_rebalance=yes` export is an
exception to label preservation and micrograph grouping. It pairs selected
particles with adjacent scores within each symmetry-aware viewing bin into
opposite halves, balances odd leftovers, updates both particle fields and
replaces affected `eo_group` markers with `eo_byindex`. This is the selector's
default; `eo_rebalance=no` preserves the original labels.
The exported project no longer preserves the original gold-standard
independence: it mixes previously refined halves and may split correlated
neighbouring picks. Original half-maps/FSC must not validate the new split;
subsequent reconstruction uses the newly recorded labels without reassigning them.

### 5.1 Sigma2 bootstrap contract

Work that has alignments, or nothing at all, but no noise-power estimate is
one problem: a fresh solve3D start, the ini3D routes, external starting
volumes, the refine3D_auto startup, a refine3D started at a later iteration in
an empty directory, and every final reconstruction at a new sampling. Since
2026-09-06 all of them follow one rule, owned by `simple_sigma2_bootstrap`:

- Seed canonical state from the particle image power spectra (`calc_pspec`),
  and let the first euclid pass replace the seed with residual sigmas. Image
  power needs neither alignment nor volume, and it is the basis every
  refinement already starts from. Half-map power (the former
  `bootstrap_rec3D` estimator) sat on a different basis than the residual
  sigmas a refinement then computes and conditioned the euclid system
  markedly worse; it is retired.
- `ensure_sigma2_for_iteration` is a no-op when the project already owns a
  compatible committed state and seeds otherwise. The consuming refine3D
  reads that state directly; the external-reference pose initialization uses
  the same transaction path after its CC residual pass.
- Refine3D is self-healing: a missing, invalid, wrong-grid, wrong-layout, or
  wrong-grouping state is rebuilt from particle power before workers launch.
- Where no refinement iteration follows (final reconstructions at original
  sampling in solve3D and refine3D_auto), the seed is upgraded by one
  residual pass: `refine=sigma` against the seeded map at the final sampling
  (no search, no volume assembly, no orientation output, alignment docs are
  not merged, no state labels seeded or written, so a multi-state project
  collapsed into one state keeps it (`labels_may_change`, as for
  `refine=eval`; until 13 October 2026 the distributed start relabelled such a
  project at random and wrote it), and the particle field left exactly as it was: shared-memory
  initialisation does not reset `updatecnt`/`sampled` for it, the run
  finalizers do not write the field, and the distributed master refreshes it
  from disk before its whole-project write), committed as the next canonical
  generation, then the shipped euclid
  ML reconstruction runs on the residual sigmas. Since 2026-09-07
  `bootstrap_rec3D` (module `simple_commanders_refine3D`) owns this whole
  sequence: seed, bootstrap map, residual pass, commit, final map; the shared
  ending `calc_final_rec` (module `simple_final_rec`, 2026-09-12), which
  closes solve3D, refine3D_auto, refine3D_states, refine3D_states_cavgs and
  classify3D_refs,
  calls it and carries no copy of the sequence. It runs standalone on any project with 3D orientations
  and is the test entry point for the final-reconstruction stage. The
  bootstrap map only serves as the residual reference, so it is always a
  gridding assembly with ML regularization (one particle pass, no
  postprocessing) that keeps the workflow's `filt_mode` and
  `automsk`: the residual sigmas depend on the regularization of the
  reference they are scored against, so that reference is regularized
  exactly as the refinement's matching references were; the shipped map keeps the
  caller's backend and `automsk` (so on PCG it is estimated on the same
  density-envelope support as the refinement, 2026-09-09; `filt_mode=none`
  keeps it classical) and, on PCG, starts from nothing at the native box and
  therefore gets the cold-solve budget of at least `FINAL_PCG_MAXITS_FLOOR`
  (5) iterations whoever the caller is (2026-09-07). The base plus ML solve
  pair of that final PCG reconstruction is inherent to ML regularization:
  the prior is built from the base pair's independent-half FSC, which is
  also the reported FSC and the unfiltered pair postprocessing uses.
- The residual pass is not an update (2026-09-27). Before this, the
  shared-memory run finalizer wrote the particle field with the pass's
  sampling bookkeeping: every active particle came back with `updatecnt` + 1,
  a new `sampled` round and cleared search statistics (poses untouched). So
  particles never updated by the refinement entered the shipped map, unlike
  the bootstrap map (`sample4rec` takes `updatecnt > 0` rows), and counted as
  updated for `update_missing` and for the frozen membership of
  `solve3D_addon`; a standalone `bootstrap_rec3D` at `which_iter=1`
  (`startit=1`) even reset `updatecnt`/`sampled` for every particle first.
  Nothing reads the bump: the pass only closes final reconstructions, and a
  later refinement starts its own sampling round.
- The unfiltered pair belongs to the assembly that wrote the map it
  accompanies (2026-09-12): every gridding `volassemble` writes the `_unfil`
  halves, the unregularized pair under ML regularization and copies of the
  halves otherwise, so `postprocess_nu` always finds a current pair; and
  postprocessing ignores a pair whose box differs from the map. A final
  reconstruction at native sampling therefore never reads the cropped pair
  left behind by the last refinement block.
- The registration-box rule is the same for both stores (2026-09-07): a
  final reconstruction whose registration (crop) box differs from the native
  box refreshes the sigmas at native sampling (image-power seed, bootstrap
  map, residual pass) before the shipped map. Before this the canonical store
  reused its committed state whenever it was structurally consumable, which
  is always true at the native box, and so skipped the refresh legacy
  performs; the two stores could not produce the same final map.

## 6. Reference Preparation

The strategy materializes the iteration reprojection model before
probabilistic pre-alignment or matcher work starts. Reference preparation:

1. reads the current Cartesian reference volumes
2. applies reference masking/filtering/centering policy
3. determines the active high-pass/low-pass shell range
4. writes even/odd PFTC reference sections for the matcher

With `automsk=yes` or `nu`, the matching references carry the selected envelope
after filtering (`automasking_policy.md`, "Reconstruction and reference use").
NU assembly multiplies the `_nu_filt` references by it: the density envelope
under `yes`; under `nu` the NU-evidence envelope when valid, otherwise the
density envelope. Outside it the NU filter field takes the coarsest bank
candidate. The matcher multiplies references that are not NU products by the
density envelope under `yes`, and under `nu` while no `_nu_filt` references
exist by the lagged NU envelope with density fallback. `automsk=yes` also
implies `envfsc=yes` (policy 2026-09-09), so FSC evaluation follows the density
envelope (post hoc with solvent correction on gridding; inside the estimator on
PCG, reported as `>>> FSC MODE`). Particles and matching-bandwidth selection
are unchanged.

`automsk=nu` requires `filt_mode=nonuniform|nonuniform_lpset`; `automsk=yes` is
accepted with any `filt_mode`. `automsk=tight` is rejected together with NU
filtering because the NU-evidence envelope has no Otsu tightness mode.

The lifecycle is lag-by-one: volume assembly from iteration N writes the NU
envelope used when iteration N+1 materializes its reprojection model. A missing
mask at the initial bootstrap is therefore expected.

`volassemble` must not generate or promote matcher reprojections. The next
iteration's reprojection model is generated by the strategy from the current
Cartesian outputs.

The model file header is the authority for the matching shell range consumed
by matcher workers and probabilistic table workers.

Reference preparation must project only the active matching shell range. It
must not silently project to the crop interpolation limit unless matching
policy explicitly asks for that range.

The 3D low-pass range comes from one of these sources:

- explicit LP-set `lp`
- project-carried NU-selected handoff in a continuing nonuniform refinement
- current FSC files
- an existing project `lp` field

`lpstop` caps the selected matching bandwidth. The crop Nyquist limit caps the
final Fourier index after the matching policy has selected a candidate range.

## 7. Matching Topology

Gold-standard matching keeps even and odd references independent.

If low-resolution even/odd docking is needed for registration, it is applied
only after reference read/mask/filter and immediately before PFTC generation.
Those blended references are not written as half-map outputs and do not feed
FSC, automasking, NU filtering, or ordinary half-map handoff.

LP-set matching uses a merged registration reference. State count alone must
not force merged-reference matching; topology is controlled by LP-set policy.
`filt_mode=uniform` is disabled for multi-state search until it has a
state-specific low-pass contract.

In nonuniform mode, reference loading follows the NU policy:

- plain `nonuniform` prefers independent `_nu_filt` even/odd references and
  falls back to regular even/odd references before using a merged map
- `nonuniform_lpset` with active LP-set matching uses the merged reference and
  prefers the merged `_nu_filt` product when present. In both modes the
  matching band is the finest member of the NU bank, the finest rung of the
  ladder cut at `fsc/1.5` or the regularized pair once it is at or beyond
  the ladder's finest rung (`nonuniform_filtering_policy.md` sections 8 and
  12, 2026-09-19)

The ordinary reference low-pass filter is not applied on top of a NU reference
path. NU and ML-regularized references are treated as already filtered during
volume assembly.

## 8. Probabilistic and Matcher Work

The current workflow is probabilistic pre-alignment followed by hard particle
assignment. It is not a monolithic soft-assignment volume-integrated EM
implementation. Every particle gets one pose and one state label, and
reconstruction follows the labels. The one exception is the per-state
M-estimation below: there reconstruction weighs each particle by a frozen
weight that `flex_pca` computed before the run. The refinement never computes
these weights from its own scores, and alignment never reads them.

Use the terms "probabilistic pre-alignment" and "hard-assignment particle
update" unless the ownership, artifacts, and update model actually change.

`refine=prob*` modes run the probability pre-step before the main matcher.
`prob_neigh` dispatches to the neighborhood probability command. The matcher
then consumes the generated assignment artifact.

The matcher's search loop (3D and 2D) runs under a dynamic OpenMP schedule: a
particle is searched whole by one thread, so at more than one thread only
which random draws a particle meets changes; one thread gives the same
results as a static schedule.

`prob_neigh_mode` controls how `prob_neigh` chooses candidate neighborhoods
before evaluating candidates:

- `state`: score coarse subspace representatives independently per state, pool
  the selected neighborhoods across states, and evaluate the same pooled
  projection search space for every active state.
- `geom`: select the subspace containing the particle's current projection for
  every active state, with no coarse scoring or pooled peaks.
- `shc`: use the direct stochastic candidate path. When shifts are enabled, it
  first estimates a shift seed before candidate scoring.
- `snhc`: use the same direct stochastic candidate path without the initial
  shift seed. It is consequently also the zero-shift default for this mode.

Probabilistic-table assignments use the calibrated likelihood path. The stored
distances are noise-normalized negative log-likelihoods for the Euclidean
objective, and evaluated candidates are sampled with weights proportional to
`exp(-dist)` over an explicit top-K support. The implementation uses a
per-particle minimum shift before exponentiation for numerical stability; this
does not change normalized weights.

The top-K truncation is deliberate. It defines the local discrete support
actually evaluated by the pre-alignment step; it is not meant to represent a
full posterior over all SO(3) grid points. If the CC objective is enabled, its
existing likelihood-compatible transformation remains the source of the
probability weights, but should not be described as a calibrated Gaussian
likelihood.

Likelihood-weighted probability-table modes may still profile or MAP-refine
shifts, and sometimes in-plane rotation, after stochastic candidate selection.
The stored assignment distance is then the refined/profiled objective value.
This is intended current behavior, not a full soft-assignment EM update.

### Dense-table assignment (`ref_assign`)

The `prob` driver assigns the particles from the dense table
(`simple_eul_prob_tab.f90`, `ref_assign`). Every reference sorts the
particles by distance, and its head is its closest particle not yet assigned.

- **One state.** Each step draws a reference from the heads' distances with
  the likelihood sampler above and assigns that reference's head particle.
  Only the heads that pointed at that particle move on.
- **Several states.** A particle's state is the state of its closest
  reference over all states. An exact tie between references goes to the
  lowest reference index, which is the label the earlier global greedy
  frontier gave. Then the same draw runs within each state, over that state's
  references and particles.

The routine keeps a count of the unassigned particles and, for every
particle, the list of references whose head it is. After an assignment it
advances only that particle's list, instead of testing every particle and
every reference again. The driver holds the table by reference
(`dists(nptcls,nrefs)`, `inpls(nptcls,nrefs)`), so each reference's
distances are contiguous for its sort. The states' labels and the check for
invalid or flat particles run with the references outside, over blocks of
256 particles in parallel, each particle keeping its best so far (strict `<`,
references in ascending order). The draws and their order are those of the earlier implementation,
so the assignment is identical (since 2026-10-16). The unit suite of the
table keeps the earlier implementation as its reference and compares the
two.

### Probability-table part files

Every table part writes one file that the driver reads after the parts have
finished.

- **The dense table (`prob`).** The header holds the number of references
  and particles, a format word (`PROB_TAB_DENSE_FORMAT`, -3, private to the
  dense table) and the number of refined entries; then the number and list
  of the active states, the particle indices, the seed table, one row per
  particle of (distance, in-plane index) for every reference
  (`prob_dense_entry`, 8 bytes; the reference and particle follow from the
  position), and a list of the refined entries (particle, reference,
  shift; `prob_refined_entry`). An unrefined entry has no shift of its own.
  The worker writes a batch's rows at the particles' positions and checks at
  the end that every row was written. The driver reads the parts in
  parallel: consecutive rows of a part in one statement (up to 8 MiB),
  scattered to the driver's columns, and each particle's refined entries
  counted during the read; the refined lists are then laid out by a prefix
  sum and every part places its entries in parallel (a particle's entries
  come from its one part, in file order). Every particle
  of the driver must be in exactly one part, and the parts must agree on the
  rotation grid and on the active states (a row's reference index counts the
  active references). Before reading past the header the file size must be
  the size the header implies, and every read checks its status. Inside the
  parallel read nothing stops: each part returns a status, and the first
  failing part stops the run after the loop, with its file name and the
  reason. `dense_candidate` rebuilds the entry `assign_candidate` reads. The
  driver's table is not initialised before the read; its refined-entry
  offsets are default integers behind a size guard.
- **The neighbourhood table (`prob_neigh`)** keeps its candidate stream
  (particle index and full candidate per entry). The neighbourhood part
  writes the number of candidates of each particle after the particle
  indices; the driver sizes its store from these counts and then places the
  candidates in file order (the order decides ties). The 2D table's parts
  are described in `solve2D_policy.md`.
- An empty part writes the header of its own table type. The neighbourhood
  table refuses `write_rows`, the dense table's row writer.

The shift search of the tables keeps a double-precision copy of the
reference for one minimisation (`load_ref_aln_8`, which records the
reference index and half per thread; the `_aln` evaluators take the index and
stop on a mismatch), and reuses the
all-rotation objective it computed at a shift when the same shift (bit for
bit) comes back, as L-BFGS-B's last accepted point does. The likelihood
scores, and the power and bounded samplers of the 2D and greedy searches,
work in per-thread buffers over arrays of distances (no internal procedures,
so no heap trampolines).

### Continuous in-plane policy

`inpl_cont` is internal since release 4 (2026-10-03), not a command-line key.
It is `yes` everywhere except under `refine=cont`, which sets it to `no` (the
Cartesian pass solves the in-plane angle itself); `no` is the historical
alternating search: continuous shift optimization with the discrete in-plane
angle callback. The convergence output keeps reporting the continuous in-plane
runs. `yes` follows the **polish-only
principle**: continuous refinement changes pose *precision*, never search
*behavior*. Selection sees exactly what the legacy route sees -- candidate
scoring, probability-table profiling, shift-seed estimation, and multi-peak
shift refinement all use the legacy discrete machinery unchanged -- and the
joint raw-Euclidean `(sx,sy,rotind_frac)` optimizer runs exactly once per
particle per iteration, after selection, to polish the committed assignment.
This makes the search trajectory identical to the validated legacy trajectory
by construction; the continuous route is a strict refinement of it and cannot
alter exploration dynamics on any sample. Objective type is a capability
check and must never activate continuous-angle behavior by itself.

Raw `objfun=euclid` and `objfun=cc` provide the analytic joint
`(sx,sy,rotind_frac)` gradient. The cc route minimizes `-cc` with a
quotient-rule angular derivative and maps scores as the clamped correlation
rather than `exp(-loss)`. Parameter validation does not couple `inpl_cont` to
`objfun` or a specific program. Objective capability
is owned by the PFTC/search implementation. The opt-in is not restricted to `refine=shc`:
deterministic, neighborhood, evaluation, and probabilistic matcher routes use
the same policy wherever they commit a pose. A mode with no pose search,
such as sigma-only setup, naturally invokes neither optimizer. Unsupported
capability combinations fail validation rather than silently reverting to the
callback. The joint cc route is validated on nanoparticle data as well, so
the nano 2D/3D workflows follow the `inpl_cont=yes` default.

The `continuous in-plane` sub-suite of the fast gate
(`simple_test_exec test=unit_pftc_align2D3D`, `simple_pftc_inplane_tester`; it
replaced `simple_test_continuous_inplane_hybrid_grad` in September 2026) guards
the hybrid route with an integer-grid score identity, finite-difference checks
of the derivative components, and construction of the production joint
optimizer.

Probability tables are pure legacy under both `inpl_cont` values: candidate
scoring, shift-seed estimation, and per-candidate shift refinement use the
legacy discrete optimizer, and the table schema remains discrete (rounded
`inpl`, score, shift in the rounded-index frame). The joint optimizer never
touches table construction; polishing candidates would sharpen the score
contrast the stochastic sampler sees and collapse exploration on difficult
samples (statistics parity below addresses the same failure through the
convergence side).

After hard state/projection/in-plane assignment, the matcher runs the joint
three-parameter optimizer locally for the selected pose -- the single joint
solve of the iteration. The assignment index is authoritative: it converts
the stored shift back to the native particle frame and seeds a solve bounded
to plus or minus two in-plane cells, without another global all-angle
selection. An accepted result persists one coupled pose: fractional Euler
`e3`, nearest integer `inpl`, shift, and score. A finite non-improving
evaluation retains the incoming discrete pose with its re-scored objective
value. A numerically invalid evaluation leaves the incoming assignment
untouched. No `inpl_cont=yes` path falls back to the legacy callback. State
and projection identity remain fixed during the final solve. A persisted
fractional `e3` is rounded to the canonical integer cell before the next
search.

**Statistics parity.** Every statistic that steers search control --
`dist_inpl`, angular distances, class/projection overlap, and anything the
convergence checker reads -- is computed from the DISCRETE cells on both
sides of the comparison, never from fractional poses. Sub-grid pose precision
would otherwise read as premature convergence and truncate the annealing
schedule that difficult samples require. The fractional `e3` and polished
shift are carried for reconstruction and persistence only.

Joint acceptance is guarded twice, in the shared optimizer. An improvement
must be material -- it must beat the seed cost by a relative tolerance far
above solver roundoff -- otherwise the solve counts as non-improving and the
seed pose stands; this keeps the improved fraction an honest diagnostic and
prevents noise-level "improvements" from committing fractional poses. A
solution pinned to any joint search bound (the shift box or the plus-or-minus
two-cell angular window) is likewise demoted to non-improving: an
angular-bound-pinned solution contradicts the exhaustive seed selection, and a
shift-corner solution is the signature of descent into truncated-series
artifacts rather than signal.

Convergence `dist_inpl` reporting is symmetry-aware. When multiple symmetry
operations have the same Euler distance within numerical tolerance, including
projection-preserving operations in dihedral groups, the reported in-plane
distance uses the equivalent operation with the smallest in-plane change. A
symmetry-equivalent branch switch must not appear as a spurious 180-degree
in-plane move.

Reconstruction, sigma, assembly, and postprocessing children must not receive
matcher-only search policy (`pose_cont`); `inpl_cont` never appears on a
command line.

For multi-state alignment (`nstates > 1`), shift-first candidate scoring is
disabled. The matcher and probability-table paths may still refine shifts after
candidate/state selection, but they must not use a shift seed estimated from a
particle's previous state to score candidates in other states.

A matcher pass reads each selected particle from the stacks once
(`doc/policies/particle_read_path_policy.md`, "The particle buffer"). The
reconstruction after the search-phase teardown takes its images from the
pass's particle buffer: the raw images, or, when it reconstructs at a cropped
box from a polar search, the cropped observations the search kept as records.
Only when the buffer does not fit its memory budget does the reconstruction
read the stacks again, once per particle. Do not add further reads unless the
performance contract is explicitly changed.

The downscaled particle cache is a 2D-only feature: `refine3D` rejects
`cache=yes` and always reads the original full-size stacks (see
`doc/policies/2D/particle_cache_policy.md`).

### Continuous Cartesian pose policy

`refine=cont` is a refinement mode of `refine3D` (contract of the pose_cont
refactoring, `doc/refactoring_notes/completed/pose_cont_refactoring.md`). It
refines the five pose parameters of every sampled particle (three rotation
components, two shifts) continuously against Cartesian central sections of
the reference, from the poses stored in `ptcl3D`. It runs none of the polar
machinery: no projection grid, no polar reprojection model, no `prob_align`,
no in-plane table. It is a continuation mode: every active particle must hold
a pose from a discrete workflow (positive state, half-set, projection
direction, finite angles and shifts); `check_cont_seeds` stops the run at
entry otherwise, before the random-orientation fallback of an unsearched
field. It searches no state: each pose is refined against its own state's
reference, and a multi-state start seeds and writes no labels for it
(`labels_may_change`, as for `refine=sigma` and `eval`), so a state left
empty, or a project in one state, is kept as it is (since 14 October 2026;
stream p07's local refinement relies on it).

Ownership:

- `src/main/cftc`: `cartft_calc` (reference slots per state and half,
  particle slots of the batch, objective, gradient, normal terms, residual)
  and `cartft_pose_opt` (one bounded Levenberg-Marquardt transaction per
  particle: the joint five-parameter stage; `cont_route=shift_then_joint` is an
  internal test seam, C16).
- `simple_strategy3D_cont`: `strategy3D_cont`, the matcher's only
  continuous-pose search object; it extends `strategy3D` beside the polar
  parent `strategy3D_pftc`. The matcher also sums `cartft_pose_work`, the work
  record of the solves, per thread for its bench report (section 11.1).
- The iteration strategy materializes the prepared real-space reference
  volumes once per iteration (`cart_refvols_{even,odd}.bin`, section 6 of this
  policy applies unchanged up to the polar projection). Every matcher reads
  them, pads them by 2 and keeps the packed padded transform per state and
  half, scaled by `padf**3` and weighted per axis by the KB stencil's response
  over the B-spline stencil's (`ref_weight`).
- Central sections are gathered on the fly with the separable quadratic
  B-spline on three taps per axis, continuous in value and derivative. On grid
  it reproduces the polar branch's KB model (KB interpolation policy,
  section 6a).
- Particle preparation follows `prepimg4align` up to the soft mask (centred
  on the stored shift, phase-flipped for `CTFFLAG_YES`) and multiplies by the
  KB stencil's taper (`prepimg4align_cart`), so the two representations share
  one sigma2. `build_batch_particles3D` fills the Cartesian slots of the batch
  from the same read.

Evaluation: one pass (`pose_normal_terms`) returns the objective, the
five-vector gradient and the 5x5 Gauss-Newton matrix, and under `euclid` the
whitened `|r|^2` per shell. It runs over the Friedel half plane of the
particle's band: for a real image every summand is Hermitian between `(h,k)`
and `(-h,-k)`, so the half-plane sums are doubled and DC is added once. The
shift phase comes from two 1-D tables per pass. The full-disk kernels are
kept (`*_full_disk`) as the tests' reference.

Objectives follow `objfun`: `cc` is `1 - cc` with uniform per-pixel weight and
no sigma2 (a Cartesian `cc` pass reads, requires and writes none); `euclid` is
the polar loss normalized by the particle's whitened power, stored as
`exp(-L)`. Its sigma2 contribution is the per-shell residual of the evaluation
that chose the committed pose (the seed, or the last accepted step), which the
transaction keeps (`get_shell_residual`) and the sigma owner unwhitens
(`calc_sigma2`); it costs no pass of its own.
The bounds are `trs` (total shift from the stored shift) and `athres_cont`
(default 10 degrees: the halfwidth of the total rotation from the stored pose,
the capture range for wrongly assigned orientations); `athres` and
`prob_athres`, which size the polar searches, do not apply. The solve is the
joint five-parameter Levenberg-Marquardt stage alone, about a third faster than
a shift-only stage before it, which the internal `cont_route=shift_then_joint`
restores for the single-pass capture tests (E25, N20 a-c); repeated passes
supply what one joint pass misses. `inpl_cont` is overridden to `no` (the
in-plane angle is part of the solve).

Steps and stopping (`cartft_pose_opt`):

- **Curvature scale.** Both objectives' Gauss-Newton matrices are multiplied
  by the regression slope of the data on the model, weighted by `|f|^2` as the
  matrices are and clamped to [0.05, 20] (`curvature_scale`). For `cc` the
  slope stands in for the correlation, by which the plain `J^T J` of a
  normalized objective overstates the curvature. For `euclid` it corrects the
  over-step of a model brighter than the data.
- **Stopping.** A stage stops:
  - when an accepted step reduces the objective by less than `1e-6` of it;
  - when a proposal predicts less than that (it is not evaluated);
  - after 8 consecutive rejections;
  - at 40 iterations.
- **Passes per solve.** The joint stage reuses the seed's evaluation and
  returns the objective at its endpoint, so a solve with N trial evaluations
  makes N+1 passes.

The rationale and the measured gains are in
`doc/implementation_notes/completed/cont_lm_acceleration.md`. The work
counters are the `CONTINUOUS POSE WORK` section of the bench report (section
11.1); `test=cont_refine3D_bench` times the solves on simulated 1JYX particles
(test environment policy).

Scores: a Cartesian pass writes `corr_cart` (slot 51 of the particle record),
the improved flag (slot 52) and the transaction status (`pose_cont_status`,
slot 42, a `CARTFT_*` code) for every particle it processes and never writes
`corr`, which keeps the value of the last polar pass. The two are not
comparable (C9). The convergence report of a Cartesian pass shows the mean
`corr_cart`, the attempts (the sampled particles) and the improved fraction.

Convergence under `refine=cont` (`check_cont_conv` in `simple_convergence`,
constants in `simple_defs_conv`): a sampled particle is stable when its
transaction ended accepted or with a finite no-improvement and it moved by at
most 0.5 degrees (`dist + dist_inpl`, an upper bound on the rotation angle)
and 1 A (`shincarg * smpd`). An iteration passes when 90% of the sampled
particles of every populated state are stable. The run converges when
consecutive passing iterations have together sampled 90% of the active
particles (`sampled` at or above the generation the streak began with); a
failing iteration resets the streak, and `minits` and `maxits` apply as for
every mode. A rejected transaction keeps the seed, so zero motion alone never
counts as convergence. The three thresholds are parameters with these values as defaults:
`cont_conv_rot` (degrees), `cont_conv_shift` (A) and `cont_conv_frac` (the
stable fraction and the streak coverage).

Polish: `pose_cont=yes` in a discrete mode makes the iteration strategy follow
every discrete pass with a `refine=cont` pass over the same particle sample
(`sample4update_reprod`), in shared memory and distributed. The discrete pass
then writes no partial reconstructions; the polish writes them, so the
iteration reconstructs once, from the polished poses. The polish writes pose,
`corr_cart`, the improved flag and, under `euclid`, sigma2 into the sigma2
transaction the discrete pass opened: it starts from the discrete pass's
residual rows and replaces those of the particles it evaluates, so a particle
whose Cartesian slot is invalid keeps the residual of its discrete pose. It
leaves `corr` and the convergence fields of the discrete search (`dist`,
`dist_inpl`, `shincarg`, `mi_proj`, `mi_state`, `frac`, `frac_greedy`), so
convergence keeps measuring the discrete search. `refine=cont` with `pose_cont=yes` is the
internal name of the polish pass and is refused on the command line.

Tests: `unit_cart_align3D` (`cart calculator`, `pose optimizer`,
`pose strategy`, `pose statistics`), `lib_cart_align3D` (1 000 simulated 1JYX
particles through `strategy3D_cont`), and the high-level gate
`cont_refine3D_1jxy` (shared memory, distributed, `refine3D_auto`, the polish),
all gated on the simulation's ground-truth orientations.

### Per-state M-estimation

`refine3D_auto state=X` refines state X of a multi-state project alone, as a
single-state refinement. `m_estimator=no|flex {no}` selects how particles enter
its reconstructions. With `m_estimator=flex` the run is a per-state
M-estimation: state X's map is estimated from every particle with a positive
weight for X in the project's state weight set, and each particle contributes
with that weight. The state weight set holds one weight per particle and state;
`flex_pca` publishes it (refine3D states policy). The weights are frozen:
nothing in the run recomputes them. With `m_estimator=no` the run refines the
particles labelled X from their hard labels, the hard comparison.
`m_estimator=flex` without `state=` is refused, and so is `m_estimator=flex`
on a multi-state `refine3D` or on `refine3D_states`: a search over states moves
the labels the frozen weights were computed for.

**Work project.** `prepare_state_work_project`
(`simple_commanders_refine3D.f90`) builds it right after the parameters are
built. `state=` requires `mkdir=yes`, so the work project and its weight set
live in the run's own directory.

- The project is copied to `refine3D_auto_stateXX.simple`, and the selection
  program runs on the copy with pruning. The selected rows are those with a
  positive weight for X (`m_estimator=flex`) or those labelled X
  (`m_estimator=no`); they reach the selection program as its per-row selection
  file. The copy holds them alone, as state 1. Each work row is checked against
  its parent row by stack and index in the stack.
- State X's map and FSC, read from the parent, become state 1's in the copy,
  where `refine3D_auto` finds them as its project volume. Every other state's
  map and FSC entries are removed, as are the parent's weight-set entry and its
  canonical sigma2 registration. The copy's sigmas are seeded in this run, as
  for any project.
- With `m_estimator=flex` the copy gets its own single-state weight set
  (`publish_work_state` in `simple_state_weight_set`): X's column at the work
  rows, the parent's hard-label flags for X, and the parent's producer and kind.
  Its manifest records the parent set's generation, layout digest and state as
  provenance. Publishing it over the parent's manifest is refused.
- The run continues on the copy as an ordinary single-state `refine3D_auto`.
  `state=` leaves the command line; `m_estimator` stays. The run directory's
  copy of the parent is deleted, so the work project is the directory's only
  project. The project given on the command line is never written.

**Reconstruction.** With `m_estimator=flex`, every reconstruction of the run
reads the work project's weight set, from the startup reconstruction to the
final native one. In each iteration the matcher's `write_partial_recs`
(`simple_strategy3D_matcher`) passes the set to `calc_3Drec` on gridding, or to
the PCG (preconditioned conjugate gradients) worker. A particle enters with its
weight when the weight is above zero on gridding, and above 0.01
(`PCG_WEIGHT_THRESHOLD`) on PCG. Volume assembly, update fractions and trailing
reconstruction count applied mass (the sum of the weights as the backend
applies them) instead of rows; see the
[importance-sampling and fractional-update policy](../importance_sampling_fractional_update_policy.md).

**Alignment and sampling are unchanged.** Every selected particle is aligned
against state 1's reference and gets one pose, as in a hard run. The sampler
never reads the weights; it chooses among the work project's rows by its usual
rule.

**Population gates** use the effective sample size (ESS: the squared sum of a
state's weights divided by the sum of their squares), not row counts. The
probability table treats a state as present when the ESS of its weights exceeds
the table's minimum population (`simple_eul_prob_tab`). The per-iteration
registration and carry-forward decision (`get_state_pops` in
`simple_refine3D_strategy`) treats a weighted state as populated when its
weights carry applied mass, which is exactly when it has partial
reconstructions.

**Final native reconstruction.** `calc_final_rec` (`simple_final_rec`) puts
`m_estimator`, and `state` when the refinement command line has one, on the
child command line it builds. That line drives both of its routes, the direct
`reconstruct3D` and `bootstrap_rec3D`. It records the weight set's identity
(generation and layout digest) before reconstructing, and stops if the identity
has changed when it registers. The `vol` entry of each weighted map carries
three separate numbers: the hard population `pop` (the parent's label count of
X, for reporting), the applied mass `mass` and the ESS `ess`. Mass and ESS are
computed over the reconstructed rows, from the weights as the backend applied
them.

The refined map of X, its FSC and the refined poses live in the work project.

### Probabilistic pre-step benchmark files and particle-frame shifts

The pre-step has benchmark reports of its own, so its time is split into named
parts instead of being one number in the strategy report (the files are in
`BENCH/`, section 11.1):

- `PROB_TAB_BENCH_PARTppp.txt`, one record per table part and iteration:
  wall-clock entries measured outside the parallel regions (setup, reference
  preparation, particle reading and preparation, table fill, table write), and
  a separate thread-seconds block (seed estimation, score loop, refinement,
  their sum, the fill region's wall time and the load balance). Thread-seconds
  are never part of `% accounted for`. The context lists the references, the
  references refined per particle, the shift searches, their alternations, the
  coarse-grid evaluations and the objective and gradient evaluations.
- `PROB_ALIGN_BENCH.txt` (`prob`) and `PROB_ALIGN_NEIGH_BENCH.txt`
  (`prob_neigh`), one record per iteration, written by the driver, all wall
  time: setup, sampling and project write, scheduling (launch to last part
  finished), table reads, `ref_assign` and the assignment write.

Every table fill with a shift search (the single-state and multi-state `prob`
tables and all `prob_neigh` fills) runs its particle loop with dynamic OpenMP
scheduling: the particles' costs differ (the shift searches, the stopping rule
of `prob_neigh`), and static blocks left threads idle at the tail. A particle
is processed whole by one thread with that thread's search object, so the
order only changes which thread's random draws a particle meets; a one-thread
run is unchanged.

Every probabilistic table holds its shifts in the particle frame (the only
path since 2026-10-15; a developer key used to switch it on). The kernels shift the reference in its own polar frame and correlate
over every in-plane rotation, so a reference-frame shift s is the
particle-frame shift `s * R(angle of the rotation)` at that rotation (what
`sh_rot` already converts). A particle-frame shift is the same at every
rotation, so the particle can be shifted once, transformed over rotation like
`ft_ptcl_ctf` (`memoize_shifted_ptcl`), and scored against the unshifted
memoised references. The kernels that score such an entry
(`gen_objfun_vals_bank`, `gen_likelihood_val_bank`) sum the angular Nyquist
bin too, so a score is exact against the true shifted loss; the
shifted-reference kernels drop that bin and differ from it slightly.

- Particle-frame seed, single-state `prob` table and every `prob_neigh` table
  that seeds (`shc` mode, and the `geom` and `state` modes in single-state
  runs) with a shift search: the seed of `estimate_shift_seed` is converted
  once at the rotation where its search ended and held as one particle-frame
  shift. Every scored reference, in the coarse and in the fine pass, is
  scored from one shifted particle (a seed below `SHERRSQ` uses the plain
  memo), a refined reference starts from the seed converted at its own
  rotation and clamped into the box, and an assigned candidate that was not
  refined records the seed as it is (`record_seed_shift`). The seed table,
  the table file and `prob_align` keep their layout.
- Particle-frame grid, multi-state `prob` table: the 5 x 5 lattice of the
  coarse shift grid is a lattice of particle-frame shifts, one entry per
  lattice point built once per particle and shared by its refined references.
  The grid minimum starts the shift search at the reference-frame shift it
  corresponds to at its rotation, clamped into the box.
- Quadratic step, multi-state `prob` table: the L-BFGS-B alternations of a
  refined reference are replaced by a quadratic fit to the lattice, where the
  grid minimum has a full 3 x 3 stencil, the fitted curvature is positive
  definite and well conditioned, the step stays inside the centre cell and the
  search box, and one exact evaluation at the fitted shift (one more lattice
  entry) beats the grid minimum. Otherwise the alternations run from the grid
  minimum.

All three are approximations of shifting each reference by the seed: away
from the seed rotation, and on the lattice, the scores, starts and recorded
shifts differ from it. The table parts' bench files record the shifted
particle memos built and (multi-state) the fits kept and the fallbacks per
reason.

## 9. Volume Assembly

When `volrec=yes`, matcher workers write partition-local Cartesian partials and
the strategy dispatches `volassemble`.

`volassemble`:

- reduces partition-local partial reconstructions
- restores dense even/odd half-volumes, deapodizes them and applies the soft
  spherical support at `msk_crop` (identical to the PCG solve support) before
  writing; the merged volume gets the same support; the support is recorded
  in the `<vol>_pcg_support.txt` sidecar so downstream consumers never mask
  again (2026-09-09)
- calculates FSC curves and state resolutions on the halves as shipped
- calculates conical FSC and cFAR from copies of the shipped halves; with
  `envfsc=yes` the copies are additionally masked by the on-the-fly density
  envelope low-pass filtered at `envmsklp`, and the phase-randomized radial FSC
  uses the shipped half-volumes to determine its 0.8 randomization onset. No
  second spherical mask is applied anywhere in the FSC evaluation. `envmsklp` defaults to `ENVMSKLP_DEFAULT` (20 A) and is separate from
  the `amsklp` NU-evidence smoothing scale
- writes `automask3D_stateNN.mrc` when `envfsc=yes`; the same density envelope
  is available to compatible non-PCG final postprocessing
- restores merged state volumes
- derives the NU-evidence envelope from the completed accepted NU bank before
  unary storage is released
- writes derived NU reference products (`_nu_filt`, `_nu_locres`) on both
  backends through the shared `simple_nu_state_filter` competition
- records diagnostics and NU matching metadata in the project: per-particle
  `res` (state FSC=0.143 resolution), `res05` (state FSC=0.5 resolution), and
  `cfar` (state conical FSC area ratio),
  the raw NU matching handoff in `lp` (clipped against `lpstop` by the
  matcher before use) and, when NU filtering is active, the same NU-estimated
  limit in `lp_est`. The convergence readout reports `RESOLUTION @ FSC=0.143`,
  `RESOLUTION @ FSC=0.5`, the band this iteration actually matched at
  (`MATCHING  LOW-PASS LIMIT (THIS ITERATION)`, from `params%lp`/`kfromto(2)`
  as fixed by `set_bp_range3D` before the search -- an explicit `lp` or the
  previous handoff) and the handoff for the next iteration (`NU HANDOFF LP
  (NEXT)`, the project `lp` field assembly writes after the search;
  2026-09-14, after an explicit `lp=3.6` run reported the handoff as the
  matching band), per state when `nstates > 1`, and
  persists them as `RESOLUTION`, `RESOLUTION_FSC05`, `LP_MATCHING` (the
  matched band), `LP_NU_HANDOFF`, `LP_ESTIMATED` (plus `_STATEnn` variants)
  in the iteration stats. Both
  backends write the same fields. `cfar` occupies particle-record slot 53;
  older, narrower binary records zero-fill appended diagnostic slots when read

Volume assembly does not refresh matcher PFTC references. It only produces
Cartesian volumes and metadata for the next iteration.

The shared helper `restore_state_from_parts` owns the restoration sequence.
Changes to ML ordering, sampling-density correction, FSC handling, trailing
reconstruction, low-resolution even/odd insertion, or LP-set behavior must go
through that shared helper rather than duplicated assembly branches.

FSCs and FSC-derived diagnostics are calculated from dense Cartesian
half-volumes, not from sparse, intermediate, or low-resolution-blended
registration-reference representations.

The nonuniform filter's source halves are the restored halves kept in memory
as they are written (the box only, as a read of the file gives them). The
trailing counts use the particle field in memory. Assembly back-fills
missing projection labels before it; the labels do not enter the counts and
are kept in the field it writes (since 2026-10-18; before, the counts reread
the field and so dropped the labels on trailing runs).

### 9.1 PCG reconstruction backend

`rec_backend=gridding|pcg` selects the reconstruction backend; `gridding` is
the unchanged default and protects all existing workflows. The PCG solver
contract, output convention, ML two-map behavior, warm starts, and
diagnostics are owned by
[reconstruct3D_pcg_policy.md](reconstruct3D_pcg_policy.md); this section
records only the refine3D-side integration contract:

- **Solver controls are distinct from refinement controls.** `maxits` counts
  refinement iterations; the PCG solve budget is `maxits_pcg` (default 2;
  more than 8 warns, more than 100 is refused) and `rtol` is PCG-specific
  (`rtol <= 0` demands exactly `maxits_pcg` iterations).
- **The backend seam is dense even/odd half maps.** `volassemble` dispatches
  to a backend-specific state restorer: gridding reduces `(cmat,rho)` and
  sampling-density-corrects; PCG reduces raw `(B,D)`, finalizes the kernel
  and solves. FSC/cFAR diagnostics, filenames, and project updates share
  policy helpers, and both backends synthesize `_nu_filt` references through
  the same assembly-owned NU competition (base pair = the unregularized
  `_unfil` halves, auxiliary member = the `P_tau`-regularized halves). PCG
  maps never receive gridding correction, a second sampling-density
  correction, or a post-hoc mask.
- **Raw statistics boundary.** Workers accumulate and atomically publish raw,
  unregularized `(B,D)` per `(state,half,part)`; only the master folds,
  finalizes, regularizes, and solves, reducing parts in ascending order.
- **Sigma ordering.** For `objfun=euclid`, iteration `n` scores and reconstructs
  with the sigma model committed before the iteration began. It writes new
  residual estimates into a pending transaction, completes assembly, and only
  then publishes them for iteration `n+1`. A stage-owned reconstruction such
  as solve3D symmetry consumes the same lagged model before the final
  transaction is published. `objfun=cc` is unweighted.
- **Weights versus priors.** Particle/data weights (including `1/sigma2`)
  multiply both `B` and `D`. The zero-mean ML prior adds precision to the
  normal operator and preconditioner only; it never weights `B`, and nothing
  prior-related is persisted in raw statistics.
- **FSC ownership.** The FSC comes from the unregularized `_unfil` base pair
  and remains the resolution authority. The regularized pair is the
  closed-form `P_tau` optimum of the base pair on the replayed operator
  (voxelwise, `shrink_by_ml_prior`, 2026-09-14; the replay solve is retired);
  NU filtering is applied afterwards by the shared competition, never inside
  the solve (policy 2026-09-06).
- **Solve support.** `automsk=yes` uses the conservative density envelope;
  `automsk=nu` prefers the lagged NU-evidence envelope and falls back to
  density; an explicit `pcg_mskfile` remains an override. The selected mask
  constrains both base and replay and `envfsc=yes` is implied. PCG reports FSC
  on the constrained pair without post-hoc masking or phase randomization,
  logged as `>>> FSC MODE`. Before either mask source exists, the base
  bootstraps on the sphere and its current pair supplies replay density.
- **Current exclusions** (hard-errored, not approximated): matrix-free
  workflow execution. Fractional/trailing
  reconstruction is implemented in the distributed master path.
- New regularization is research, tracked in
  `doc/implementation_notes/completed/pcg_priors_history.md`; it cannot be used to close
  integration gates.

## 10. Trailing and Combined Even/Odd

Trailing reconstruction blends the even accumulators (sums and sampling
density) with the previous even chain and the odd with the previous odd chain,
in the accumulator domain (importance-sampling policy, section 7); finished
half maps are never blended. It runs before restoration, NU-evidence envelope
generation and derived reference filtering.

In LP-set mode, if trailing is active, assembly blends the even and odd
accumulators first and restores the merged volume from the blended sums.
Assembly must not use low-resolution even/odd docking insertion for LP-set
outputs.

When `combine_eo=yes`, distributed `refine3D` schedules one additional final
iteration after convergence or run-length termination. That final iteration:

- disables fractional update
- forces full update
- sets `combine_eo=yes`
- tightens the low-pass criterion to at most 0.143

The combined even/odd iteration is part of base `refine3D`, not a terminal
`refine3D_auto` or `solve3D` reconstruction step.

## 11. Finalization and Artifacts

### 11.1 Benchmark reports

Benchmarking (`L_BENCH_GLOB`, on in every build) writes its reports into one
subdirectory of the run directory, `BENCH/`, never into the run directory
itself. Each benchmark family has one file there, and a family that every
partition reports has one file per partition. A report is appended to its
file as one record, so a run of any length writes a fixed number of files: for
a distributed solve3D run with 10 partitions and both probability-table modes,
about 45 instead of one file per iteration and partition (3,104 for a solve3D run of 94 refine3D iterations).

| file | written by | records |
| --- | --- | --- |
| `REFINE3D_STRATEGY_BENCH.txt` | the refine3D strategy | one per iteration |
| `REFINE3D_BENCH_PARTppp.txt` | each 3D search partition | one per iteration |
| `REC3D_BENCH_PARTppp.txt` | each partition's reconstruction (gridding or PCG) | one per pass |
| `PROB_TAB_BENCH_PARTppp.txt`, `PROB_TAB_NEIGH_BENCH_PARTppp.txt`, `PROB_TAB2D_BENCH_PARTppp.txt` | each probability-table part | one per iteration |
| `PROB_ALIGN_BENCH.txt`, `PROB_ALIGN_NEIGH_BENCH.txt`, `PROB_ALIGN2D_BENCH.txt` | the probabilistic pre-step driver of `prob`, of `prob_neigh`, of the 2D tables | one per iteration |
| `VOLASSEMBLE_BENCH.txt` | volume assembly | one per iteration |
| `PCG_ASSEMBLE_BENCH.txt` | the distributed PCG master (`rec_backend=pcg`) | one per iteration |
| `SIGMA2_COMMIT_BENCH.txt` | the sigma2 commit | one per commit |
| `SOLVE3D_BENCH.txt` | `solve3D`, `solve3D_addon` | one per run |
| `REFINE2D_STRATEGY_BENCH.txt`, `REFINE2D_BENCH_PARTppp.txt` | the refine2D strategy; each 2D search partition | one per iteration |
| `SOLVE2D_BENCH.txt` | `solve2D` | one per report |
| `CAVGASSEMBLE_BENCH.txt`, `CAVGASSEMBLE_SUMS_BENCH.txt` | class-average assembly; its sum of the parts' contributions | one per assembly |
| `CAVGSUMS_BENCH_PARTppp.txt` | each partition's offline class sums | one per call |

A record starts with a `*** BENCHMARK RECORD ***` header that names the family
and, where they apply, the iteration (`which_iter`) and the partition. The
writer's own sections follow, as before: `BENCHMARK CONTEXT` and `TIMINGS (s)`,
plus `THREAD-SECONDS` for the table parts, the 2D and 3D searches and the class sums
(busy time summed over the team, the region's wall time, the team size and
the load balance), and `RELATIVE TIMINGS (%)` where a
writer reports them. A Cartesian pass (`refine=cont`, the polish) adds
`CONTINUOUS POSE WORK` to its `REFINE3D` record: the solves by `CARTFT_*`
status, LM iterations, trial evaluations, accepted steps, bound hits,
calculator passes (the sigma2 residual is the committed evaluation's, with no
pass of its own), pixel-evaluations (passes
times the pixels of the particle's shells in the full disk), the solves that
reached the iteration cap, a histogram of actual over predicted reduction at
accepted steps, and the search
thread-seconds per pixel-evaluation. Its context gives the bytes of the
Cartesian references, and their read and transform count as reference
preparation. Each partition appends only to its own file, so
concurrent workers never write to the same file. The shared helper is
`open_bench_report` in `src/fileio/simple_bench_io.f90`; a new benchmark
writer uses it rather than opening a file of its own. SIMPLE never truncates
these files, and a record does not identify its run: clear `BENCH/` before a
second run in the same directory.

`scripts/plot_refine3d_bench.py` reads a `BENCH/` directory, or a directory of
the per-iteration files of earlier builds, and writes the per-iteration
summary table, a long table of every record, a text report and the figure;
for refine2D iterations also their critical path (`bench_iterations2D.csv`,
`bench_profile2D.svg`).

`VOLASSEMBLE_BENCH.txt` has, after its `% accounted for` line, a section
`SPLIT TIMINGS (s)`: the nonuniform filter split into candidate filtering,
objective, label optimisation and composition, and per state the restore and
the filter. These are parts of the timings above and are not added to the
accounted total. The summary script adds the filter split as columns of
`bench_iterations.csv`.

`PCG_ASSEMBLE_BENCH.txt` is the PCG backend's counterpart, written by
`assemble_refine3D_pcg` around `execute_rec3D_pcg_distributed_master` and the
NU low-pass handoff. Its `TIMINGS (s)` are disjoint wall-clock phases summed
over states and halves: init, solve support, operator setup, raw read (I/O,
headers and validation of the partitions' accumulators, or of a carried
chain), raw summation, trailing chain (scaling, reading and blending the
previous chain, writing the new one and its manifest), operator finalize
(for the ml pair only `P_tau` and the preconditioner, since it reuses the base
pair's operator and reads nothing), the base, prior and ml solves (wall time
of the concurrent pair), solve finish (validation and the diagnostics files),
FSC and cFAR, map output, the nonuniform filter, project write, cleanup, the
handoff, total and `% accounted for`. Its context gives the raw files and
bytes read and the solve iterations per half; its `SPLIT TIMINGS (s)` give
the operator finalize per solve kind, each half's solve time and the NU
filter split. The summary script
splits the strategy's assembly into volassemble, four PCG master groups (raw
reduction; operator and solves; FSC, maps and NU filter; other) and the rest,
and reports the PCG phases over the run and by stage.

On each iteration, the strategy report should stay simple: context plus one
`TIMINGS (s)` section. Labels should include the coarse operation buckets
setup, probabilistic pre-step, matcher/scheduler, assembly/postprocess, and
total time, together with separately timed reprojection-model materialization
and the sigma2 transaction (the candidate preparation and the commit). These
measurements keep a per-iteration sigma update distinguishable from reference
preparation.

### 11.2 Run artifacts

Base `objfun=cc` refinement still does not bootstrap or consume sigma2. The
external-reference wrappers are the explicit exception: their shared service
runs `calc_pspec` before the fixed-reference CC pass, preserves that native-grid
partition state outside the initialized cohort, and replaces active cohort
shells with residual sigma2 after committed assignments.

Distributed matching writes partition alignment documents and merges them into
the project after worker completion. Shared-memory matching writes the project
directly after each iteration.

After the current iteration's committed orientations are available, both
strategies attempt to write `orientations_distribution_stateNN.jpg` for each
nonempty state. This projection-direction distribution is a secondary,
latest-view visualization, generated and overwritten at every iteration. It is not
authoritative scientific state, is not registered in `os_out`, and is not a
restart or workflow-handoff artifact. An empty state may produce no
new image. Failure to write this JPEG warns but does not invalidate the
refinement iteration.

On finalization:

- `endit` is written to the command line
- `startit` is removed from the command line
- active state volumes and FSCs are registered in `os_out`
- empty states are removed from `os_out` by the shared-memory finalizer; the
  distributed run leaves the entries of such a state as its iterations
  registered them. Under per-state M-estimation a state is empty when its
  weights carry no applied mass
- `cls3D` distributed runs map class-orientation output back to particles
- `JOB_FINISHED` is touched by the shared-memory path

The original-sampling final reconstruction is distinct from an ordinary
refinement iteration. On the PCG backend, every `calc_final_rec` caller
applies the shared minimum five-iteration budget to this cold solve. An explicit
positive residual tolerance may still stop convergence earlier. Final
automatic sharpening follows the isotropic postprocess protocol
(2026-09-21, the postprocess_nu v2 recipe with one cutoff): the cutoff is
the FSC=0.143 of the reconstruction's FSC file (the base pair's curve,
envfsc-corrected where that applies), or of the `_even_unfil/_odd_unfil`
pair beside the map computed by postprocess when no file is given, one
Guinier B-factor of the density-windowed unfiltered pair average
(windowing restored 2026-10-16; see the PCG policy's final-sharpening
section for why a support-wide fit under-estimates |B|) between `HPLIM_GUINIER`
(10 A since 2026-09-26, RELION's `--autob_lowres` default and Rosenthal &
Henderson's lower limit; 20 A before, which put the steep envelope/micelle
region into the fit) and that cutoff (never of the shipped regularized
map: its prior's amplitude suppression steepens the slope, -150 on
streptavidin against -77 to -83 from unregularized maps), sharpen, the
FSC weighting exactly once, then the Butterworth low-pass at the cutoff.
Exactly once (2026-09-26): an ML-regularized map (support-provenance
`solve_kind=regularized`, or `gridding_regularized` from a gridding
reconstruction with `ml_reg=yes`) already carries it -- its prior shrank
every Fourier component by `rho/(rho + <rho>/(tau SSNR_half))`, about FSC
per shell -- and is only B-sharpened and Butterworth-filtered; every
other map (the `imgkind=unfil|solvent` pair averages, base or mixed PCG
solves, unregularized gridding maps, foreign maps) gets RELION's
weighting `sqrt(2FSC/(1+FSC))` (Rosenthal & Henderson's C_ref,
`postprocessing.cpp` `applyFscWeighting`, zero from the first shell with
FSC < 1e-4). The 2026-09-22 recipe applied the Wiener `2FSC/(1+FSC)` on
top of the regularized map's own shrinkage, leaving 0.036 of the
amplitude at FSC=0.143 (0.33 at 0.5): every map over-smoothed. Its
trigger -- exp_gate a cloud of structured noise under the Butterworth
alone at an automatic B of -108 while -50 looked right -- was a defect in
the closing filter (2026-09-27): the Butterworth array was sized to the
box, so `apply_filter` reached the corners of the Fourier cube (out to
sqrt(3) x Nyquist), where `apply_bfac`'s exp(-B s^2/4) is largest, with
only the Butterworth's k^-8 tail against it. At exp_gate's 0.822 A/pixel
the corner gain is ~1e8 at B -108, ~1e5 at -80, ~4e2 at -60 (msp1 at
1.073 A: ~4e4 at -115); the noise sits outside the molecule because
there is no density there to hide it. Any filter sized to the FSC zeroes
the corners, which is why the FSC weighting and the NU sharpening never
showed it. The Butterworth now covers the Fourier shells up to Nyquist
only, and the sharpening is capped at the cutoff shell (2026-09-28): every
shell beyond the cutoff keeps the cutoff's gain exp(-B s_c^2/4) and the
Butterworth alone rolls it off (`bfac_cap_filter`), since there is no
signal beyond the cutoff to restore and at fine pixels exp(-B s^2/4)
otherwise outgrows the order-8 Butterworth from ~1.5 x the cutoff on
(exp_gate's Nyquist shell: x45 at B -108 uncapped). RELION differs in two
places kept on purpose: it fits B on the FSC-weighted map out to where
the FSC reaches zero (steeper) and applies no low-pass by default; SIMPLE
fits the unweighted unfiltered pair up to the FSC=0.143 cutoff (a
conservative B) and closes with the Butterworth there. No
density-windowed pair estimate. With
`pcg_solvent=yes` the `_unfil` pair is the prior-free base pair and the
sharpened map is the shipped (prior'd, replayed) map; `imgkind=unfil` and
`imgkind=solvent` postprocess the respective pair averages with the same
cutoff.

The transaction files of the sigma2 store (the candidate and the workers' range
files) are run-local. They are written and consumed inside a running
refinement and are not registered as `os_out` handoff artifacts; only the
committed store is registered, in `projinfo`.

## 12. Refactor Rules

- Preserve the particle-domain versus volume-domain boundary.
- Keep shared-memory and distributed workflows behaviorally aligned.
- Keep volume assembly as the owner of assembled-reference work.
- Do not move assembled-volume postprocessing into `refine3D_strategy`.
- Do not merge probabilistic particle-update logic with volume postprocessing.
- Do not introduce a second ambiguous source of matching references.
- Preserve the single image-stack read per matcher batch.
- Treat assignment maps, partition-local partials, state volumes, even/odd
  volumes, FSC files, automasks, and NU outputs as explicit workflow
  contracts.
- Preserve the two-mode in-plane boundary: only `new_alternating` may attach the
  callback, and no `inpl_cont=yes` path may invoke it.
- Initialize candidate-profiling joint solves with one callback-equivalent
  all-angle selection at the supplied shift; do not run the 5-by-5 coarse shift
  initializer.
- Treat the final 3D hard assignment as authoritative: refine locally around
  its rounded in-plane cell and never perform a second global angle selection.
- Keep probabilistic in-plane artifacts rounded and run durable fractional
  refinement only after final hard assignment.
- Keep `inpl_cont` internal: no UI input, no command-line key, `no` only under
  `refine=cont`.
- Keep state weights out of alignment and sampling: they enter only
  reconstruction, its mass bookkeeping and its population gates, and no
  refinement step recomputes or changes them.
