# Solve2D Policy

## 1. Purpose and Scope

This document defines the current architectural policy for `solve2D` and the `refine2D` workflow it drives.

`solve2D` is de novo 2D class-average determination, the 2D equivalent of
`solve3D` (see `doc/policies/3D/solve3D_policy.md`), and mirrors its
methodology: one staged schedule couples ab initio 2D classification from
random references with initial 2D refinement of the resulting classes. It was
called `abinitio2D` before 2026-10-03; the old name still runs. Its stages are
`refine2D` runs (named `cluster2D` before 2026-10-03), the 2D counterpart of
`refine3D`, as the stages of `solve3D` are `refine3D` runs.

It mirrors the recent `refine3D` cleanup where the same design pressure exists:

- keep stage policy separate from execution mechanics
- keep particle/class assignment work separate from class-average assembly
- preserve shared-memory and distributed parity
- make sampled-update and probabilistic handoffs explicit
- treat class-average files, assignment files, FRCs, and partial sums as workflow contracts rather than incidental scratch files

The 2D workflow is intentionally Cartesian. The old `polar` command-line branch
selector has been removed for `solve2D` and `refine2D`.

## 2. Architectural Policy

`solve2D` is a staged 2D classification workflow:

1. set run defaults and read project state
2. determine stage geometry, low-pass limits, sampling policy, and refine2D command lines
3. initialize references when needed
4. run staged `refine2D` iterations
5. optionally run probabilistic pre-alignment in later stages
6. update particle class, in-plane, shift, sampled, and update-count state
7. restore class averages through shared-memory or distributed class-average pathways
8. run a final fill-in assignment pass for active particles that were never updated
9. when sampled updates were active, run a terminal dense greedy all-particle pass
10. generate final class averages, FRC metadata, and ranked outputs

The main policy boundary is:

- particle-domain work owns particle sampling, probabilistic assignment tables, search, class assignment, shift/in-plane updates, and partition-local outputs
- class-average assembly/restoration owns class-average sums, even/odd outputs, merged class averages, FRC/class documents, and project output metadata

Command-line `lp` is a fixed low-pass override for the ML-regularized
`solve2D` stages: the pre-ML Gaussian-reference stage keeps the automatic
starting low-pass so Gaussian regularization still acts, and the fixed `lp`
takes effect once `ml_reg` is active. `solve2D_chunks` must preserve that
behavior when constructing child `solve2D` command lines, applying only the
chunk-local Nyquist floor.

The staged `solve2D` controller uses sampled SNHC (`refine=snhc_smpl`) for
stages 1-2. From stage 3 onward, `refine=prob` requests dense probabilistic
assignment at every stage. For `refine=prob_snhc`, intermediate stages use
sparse probabilistic SNHC and the final staged invocation uses dense
`refine=prob` so the previous class remains a valid assignment candidate and
class-overlap convergence reporting can recover. The separate terminal
all-particle coverage pass after sampled staged updates also uses dense
`refine=greedy`. `solve2D_chunks` must preserve this policy when constructing
child `solve2D` command lines.

### Seeded restart (`cls_init=prev`)

`cls_init=prev` re-enters the workflow from a previous 2D clustering held in
the project instead of a random start. It is a `solve2D`-only mode
(`refine2D` rejects it; `solve2D_chunks` and the stream keep
`cls_init=rand`). Design record:
`doc/implementation_notes/completed/abinitio2D_seeded_restart.md`.

- The seed partition is built from metadata only (`ptcl2D` `class`, `state`,
  `corr`; `cls2D` `state`): no image is read, registered or split before
  the search. Accepted parents are classes with `cls2D%state > 0` (every
  labelled class when no selection state exists) and at least
  `MINCLSPOPLIM` active particles. Seed classes are allocated to parents by
  largest remainder proportional to population (at least one per parent
  when `ncls >= nparents`), so every seed class holds about
  `nptcls/ncls` particles and the seed set represents the previous view
  distribution; a parent with several seed classes is split by rank
  interleaving on `corr`. When `ncls < nparents` the least populous parents
  are dropped. Particles of dropped, rejected, under-populated or
  unlabelled classes get `class=0` and are assigned in the seed pass
  (`oris%reseed_classes`).
- The only hard error is the absence of a previous clustering (virgin
  `ptcl2D`, or no active particle with a class label). Missing `cls2D`
  state or labels beyond `cls2D` are repaired with a counted warning; a
  project with no even/odd partition gets one (`partition_eo`, before the
  sigma2 state is built). Existing `eo` values are never touched, and no
  per-particle `eo` check is made: `isthere('eo')` is false for `eo=0`.
  `delete_2Dclustering` is never called.
- Seed references are `make_cavgs` from the seed labels at the working
  `box_crop` (`start2Drefs*`), made after the canonical sigma2 state has
  been validated or rebuilt (`ensure_resume_sigma_state`).
- The run is entered at `PROBREFINE_STAGE` through the same path as a
  stream checkpoint resume, preceded by one seed pass: a single dense
  `refine=prob` `refine2D` iteration of every active particle
  (`update_frac`, `nsample` and `fillin` deleted, `extr_iter=extr_lim+1`)
  at the low-pass limit of `PROBREFINE_STAGE`. The pass runs as the
  iteration at which stage `PROBREFINE_STAGE-1` would have ended
  (`solve2D_seed_pass_iter`), so the stages from `PROBREFINE_STAGE` to
  the terminal greedy pass run exactly as on an unseeded run: same limits,
  refine policy, iteration counts, sampling and fractional restore. The
  pass must not be a fresh start (`startit > 1`): a fresh start zeroes the
  shifts after `prob_tab2D` has built its table against them.
- Diagnostics: `>>> SOLVE2D SEED` (parents, seed classes, dropped and
  unassigned counts, seed populations), `>>> SOLVE2D SEED PASS
  REASSIGNED` (% class changes, % shifts moved > 1 px, % seed classes
  retaining at least half their members) and `seed_lineage.txt` (seed
  class, parent, seed population, final population).

## 3. Ownership Policy

`simple_commanders_solve2D.f90` owns:

- the `solve2D` entry point
- top-level defaults
- run orchestration across stages
- initial reference handling, including the `cls_init=prev` seed
  (validation/repair of the previous clustering, seed references, seed
  pass and its diagnostics); the seed partition itself is an `oris`
  operation (`reseed_classes`, `simple_oris_reshape.f90`)
- final fill-in dispatch
- terminal dense greedy all-particle dispatch after sampled staged updates
- final class-average generation/ranking

This layer should stay thin enough that stage rules are readable elsewhere.

`simple_solve2D_controller.f90` owns:

- stage counts and stage constants
- low-pass limit helpers
- stage-local `refine2D` command construction
- search-mode policy by stage
- sampled-update policy, including `NSAMPLE_DEFAULT_2D` and `nsample` override handling
- the rule that stage 1 may sample particles but does not fractionally restore previous class averages
- the seed-pass iteration index and command line of a `cls_init=prev` run
  (`solve2D_seed_pass_iter`, `set_cline_refine2D_seed_pass`)

`simple_refine2D_strategy.f90` owns:

- shared-memory versus distributed execution selection
- iteration control inside one `refine2D` invocation
- scheduler interaction
- probabilistic pre-alignment dispatch
- distributed worker scheduling
- distributed class-average assembly dispatch
- convergence and run-finalization bookkeeping

`simple_strategy2D_matcher.f90` owns:

- particle-domain alignment/search
- reproduction of the probabilistic sampled subset when `prob_align2D` is active
- strategy-object selection
- sigma updates during Euclidean search
- writing orientation updates
- writing distributed partial class-average sums when running as a worker

The 2D matcher must preserve a single particle-stack read per batch in the
online alignment/restoration path. Batch construction should keep the
already-read raw particle images for class-average restoration, and restoration
should consume those in-memory images after assignment in the same batch.

Do not split online class-average restoration into a second full particle pass
that re-reads image stacks to lower peak memory. Offline or terminal
class-average assembly commands may have their own explicit reads, but that is
separate from the matcher worker's online single-read contract.

Probabilistic table construction has a separate bounded-memory contract:

- `prob_tab2D` workers retain only thread-local compact candidates for the
  current particle batch and stream them to one partition file;
- `prob_align2D` constructs the global object only after workers complete;
- dense `refine=prob` uses a compact rectangular candidate table, while sparse
  `refine=prob_snhc` uses a particle-oriented ragged candidate store;
- compact storage is materialized as the established `ptcl_ref` assignment only
  after global selection, so assignment-file semantics do not change.

`simple_commanders_mkcavgs.f90` and the classaverager modules own:

- explicit class-average assembly from partial sums
- merged/even/odd class-average output
- class-document generation
- class FRC output and project output metadata

## 4. Sampling and Fractional-Update Policy

`solve2D` uses a fixed run-local target sample size:

- default: `NSAMPLE_DEFAULT_2D = 200000`
- override: `nsample=<integer>`

The effective update fraction is:

```text
update_frac_2D = min(1.0, real(min(nptcls_eff, nsample_target_2D)) / real(nptcls_eff))
```

where `nptcls_eff` is the number of active particles with `state > 0`.

Stage policy:

- stage 1 uses a random sampled subset but disables fractional carry-over of previous class-average sums
- stages 2 and later use sampled update with fractional class-average restoration when the sample is smaller than the active set
- a `cls_init=prev` run has no sticky stage: the seed pass updates every
  active particle (and sets `updatecnt=1` on all of them), then the stages
  from `PROBREFINE_STAGE` on sample and restore as they do on any run
- probabilistic stages preserve sample-once-and-reuse: `prob_align2D` chooses the subset, and `prob_tab2D`/`refine2D_exec` reproduce that subset rather than resampling
- solve2D uses likelihood-weighted probabilistic assignment: raw objective
  distances are kept and evaluated class/in-plane candidates are sampled with
  weights proportional to `exp(-dist)` over the explicit top-K support selected
  by the current probabilistic mode
- variance-normalized Euclidean distances are
  likelihood-like negative log weights; `objfun=cc` is supported as a monotone
  pseudo-likelihood with `dist = 1 - clamp(cc, 0, 1)` before applying
  `exp(-dist)`
- likelihood-weighted 2D modes may still profile/MAP-refine shifts, and
  sometimes in-plane rotation, after stochastic candidate selection; the
  assignment table then stores the refined/profiled distance, which is intended
  current behavior rather than a full soft-assignment EM update
- staged `fillin=yes` currently acts as a full-assignment coverage guard: it
  keeps iterating until active particles have assignments, while particle
  selection still follows the normal sampled-update path
- if any staged update used `update_frac`, `solve2D` runs a terminal
  dense `refine=greedy` all-particle `refine2D` pass with `update_frac` and
  `fillin` disabled, refreshing class, in-plane, and shift parameters before
  final class-average generation

The restoration model is class-local: each class carries forward its previous sums under the population rule of the [class-average state note](../../refactoring_notes/completed/class_average_and_reconstruct3d_partials_refactoring.md) (Section 4.1). The carried set records `M(c)`, the population its sums represent; the assembly owner counts `N(c)` (active, updated rows of the class) and `n(c)` (those sampled this round) on the merged project and blends `new = current + w * shift(previous)` with `w = (N - n) / M`, recording `M <- N`. The carried mass therefore equals the represented population whatever joined or left the class; without a population change `w = 1 - f`. The rule keeps the mass right, not the membership: the carried sums are an approximation (a stochastic recurrence), not an exact sum over the current class members, and old contributions are removed in proportion, not particle by particle. Shared memory is the reference behaviour; distributed runs give the same result, independent of `nparts`.

### Continuous in-plane policy

`inpl_cont=no|yes` is propagated by the `solve2D` controller to every
`refine2D` child, including probabilistic staged calls, the final staged
invocation, and the terminal dense all-particle refresh. `no` preserves the
historical alternating shift/discrete-angle callback route. `yes` is the
default and replaces every callback-based angle/shift optimization with the
joint raw-Euclidean `(sx,sy,rotind_frac)` optimizer.

During candidate profiling, each joint invocation keeps the selected class
fixed but discards every incoming in-plane index or fractional coordinate. At
the caller-supplied native `(x,y)` shift seed, it performs exactly one
all-angle discrete evaluation and selects the best grid index, matching one
invocation of the legacy callback. The joint solve then starts from that index
with a local plus-or-minus-two-cell angular bound. This initialization never
scans alternative shifts: the legacy 5-by-5 shift/all-angle coarse initializer
is not part of `inpl_cont=yes`.

The active joint route is supported by non-streaming, non-time-series search
under the raw Euclidean and cc objectives (the cc route minimizes `-cc` with
a quotient-rule angular derivative and reports the clamped correlation as its
score). The hybrid/denoised blend is not a continuous-angle capability and
fails rather than silently selecting the legacy callback. Time-series
shift-only search uses its fixed-angle optimizer and invokes neither angle
route.

Probabilistic particle and class/reference sampling remain discrete. During
candidate profiling, the joint optimizer may evaluate a fractional angle, but
the probability artifact carries the rounded `inpl`, score, and shift only.
The shift is stored in the rounded-index frame. Once the final class/in-plane
assignment has been selected, the matcher reruns the joint optimizer LOCALLY
for the selected pose: the assignment is authoritative, so the solve seeds at
the table `inpl` and its recovered native-frame shift, bounded to plus or
minus two in-plane cells, with no further global all-angle reselection --
under near-degenerate in-plane branches (rotationally self-similar class
averages) a second global selection at a slightly different shift hops
branches on floating-point noise. An accepted result persists the fractional
`e3`, nearest integer `inpl`, shift, and score as one pose. A valid
non-improving run retains the incoming pose with its re-scored objective
value; a numerically invalid run retains the incoming assignment untouched.
Neither outcome falls back to the callback. Joint acceptance guards
(material-improvement tolerance, bound-pinning demotion) are shared with
refine3D; see the continuous in-plane section of
[refine3D_policy.md](refine3D_policy.md).

## 5. Iteration Semantics

For `refine2D`:

- `startit` is the stage/invocation start
- `which_iter` is the current iteration
- `extr_iter` tracks the 2D extrapolation/search schedule
- `endit` is written after an invocation finishes and is consumed by the next stage setup

Do not collapse these counters into one another. Child command lines, including probabilistic pre-alignment and fill-in, must preserve the distinction between stage start and current iteration.

## 6. Artifact and Handoff Policy

Stable 2D workflow artifacts include:

- `assignment_part*.dat` and `assignment.dat`
- `dist_part*.dat` and `dist.dat`
- `cavg_contrib_part<N>.bin`: one worker's current-iteration class sums (four
  unregularized accumulators, centering offsets, populations), accumulated from zero
- `cavg_state.bin`: the carried class sums of the run, written only by the
  assembly owner, with `M(c)` per class; no part number
- `cavgs_iterNNN.mrc`, `cavgs_iterNNN_even.mrc`, `cavgs_iterNNN_odd.mrc`
- `FRCS_FILE`
- `sigma2` iteration files
- `ptcl2D`, `cls2D`, `cls3D`, and `out` project segments

Partition-local probabilistic assignment/dist files and the class-sum contributions are per-iteration artifacts: the master removes them before the next distributed iteration writes new ones, and the assembly owner deletes the contributions after the blend. Workers never read carried state. The carried set `cavg_state.bin` is the only durable source of class-average carry-over: the owner (the `cavgassemble` step in distributed runs, the matcher's restoration finalisation in shared memory) reads it, applies the class-centering shifts the workers report once, blends, and publishes the new set through a temporary name and a rename. When an iteration would blend but the carried set is missing, unreadable, or disagrees with the run on class count, `box_crop` or `smpd_crop`, the master runs that iteration as a full update (every particle, no carry-over). Legacy `cavgs_*_partN.mrc` / `ctfsqsums_*_partN.mrc` files are never read; a run continued from such a directory starts with one full update. Changing `nparts` or the execution mode across `continue=yes` needs nothing else.

## 7. Review Checklist

For any `solve2D` or `refine2D` change, check:

- Does the command layer remain mostly orchestration?
- Is stage policy in the controller rather than scattered through matcher or strategy code?
- Does probabilistic 2D sample once and then reproduce the same subset?
- Do shared-memory and distributed paths use the same scientific workflow?
- Are class-average assembly/restoration responsibilities explicit?
- Does online class-average restoration reuse the matcher batch images instead
  of introducing a second particle-stack read?
- Does probabilistic table construction stay batch-bounded in workers, avoid
  worker/global overlap, and keep sparse global storage proportional to the
  evaluated candidate count?
- Are stale distributed handoffs (assignment, dist and class-sum contribution files) removed, while the carried set `cavg_state.bin` is left to its owner?
- Are `startit`, `which_iter`, `extr_iter`, and `endit` semantics preserved?
- Is `fillin=yes` treated as a full-assignment coverage guard unless the
  implementation is deliberately changed to missing-only assignment?
- When staged updates are sampled, does terminal dense greedy assignment refresh all active
  particles before final class-average generation?
- Does the change preserve Cartesian-only `solve2D`?
- Does `inpl_cont=no` retain the callback route and `inpl_cont=yes` avoid it
  throughout deterministic and probabilistic search?
- Does candidate profiling reselect its discrete seed at the supplied shift
  (no previous `inpl`, fractional restart coordinate, or 5-by-5 shift scan),
  while the post-assignment durable pass refines locally around the
  authoritative table `inpl` without another global reselection?
- Do probability artifacts remain rounded while final assignment alone owns
  durable fractional `e3`?
- Does `cls_init=prev` keep existing `eo` values, never call
  `delete_2Dclustering`, build its seed from metadata only, and leave the
  stages from `PROBREFINE_STAGE` on identical to an unseeded run?

## 8. Rules to Preserve During Refactors

- Do not reintroduce a `polar` branch selector into `solve2D`.
- Do not bury stage-policy tables in the matcher.
- Do not let probabilistic pre-alignment and matcher update sample different particle subsets.
- Do not make distributed-only class-average assembly semantics diverge from shared-memory scientific behavior.
- Do not describe `fillin=yes` as missing-only assignment while it still uses
  the normal sampled-update path.
- Do not use final fill-in as a substitute for the terminal dense all-particle
  refresh when sampled solve2D updates were active.
- Do not reuse stale assignment files as valid current-iteration inputs.
- Do not re-read particle stacks in the online matcher/restoration path when the
  raw batch images are already available.
- Do not let workers read or write the carried class sums, and do not reintroduce partition-shaped carry-over: workers accumulate the current sample from zero; only the assembly owner blends and publishes `cavg_state.bin`.
- Do not add fractional in-plane coordinates to probabilistic assignment
  artifacts; rerun the joint optimizer after final assignment instead.
- Do not invoke the legacy angle callback from any `inpl_cont=yes` failure or
  no-improvement path.
