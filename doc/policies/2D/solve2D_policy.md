# Solve2D Policy

## 1. Purpose and Scope

This document defines the current architectural policy for `solve2D` and the `refine2D` workflow it drives.

`solve2D` is de novo 2D class-average determination, the 2D equivalent of
`solve3D` (see `doc/policies/3D/solve3D_policy.md`), and mirrors its
methodology: one staged schedule couples ab initio 2D classification from
random references with initial 2D refinement of the resulting classes. It was
called `abinitio2D` before 2026-10-03 (no alias is kept). Its stages are
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
8. when updates are sampled, keep the last stages iterating (`fillin=yes`) until every active particle has a class
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

### Defaults

`solve2D` sets `oritype=ptcl2D`. When unset, it supplies:

- `autoscale=yes`, `mkdir=yes`
- `center=yes`, `center_type=seg`
- `cls_init=rand`
- `gauref=yes`, `ml_reg=yes`
- `refine=prob_snhc`
- `extr_lim=20`, `nits_per_stage=5`, `eo_stage=yes`
- `rank_cavgs=yes`, `stats=no`

The controller derives the rest:

- the crop: `box_crop` from a target pixel size of 2.67 A (never below a box
  of 88) when `autoscale=yes`; the class averages of the stages live at it
- `lpstart` = mask diameter / 12, kept within 15-20 A; `lpstop` = mask
  diameter / 22, kept within 6-8 A; `cenlp` = mask diameter / 6, kept within
  20-30 A; each at least twice the cropped pixel size. A command-line `lp`
  fixes every stage after the first at `lp`
- `nsample=200000` (`NSAMPLE_DEFAULT_2D`): the sample of each iteration, every
  active particle when there are fewer
- `nstages=6` (5 with `eo_stage=no`); a smaller `nstages` re-spaces the
  low-pass ladder from `lpstart` to `lpstop` over the stages before the last,
  which stays the gold-standard stage

The interface shows these values (`scripts/default_audit.py` follows
`exec_solve2D` into `exec_solve2D_workflow`). The shown numeric defaults of
`lpstart`, `lpstop` and `cenlp` are the parameter baselines, since the
interface has no way to show a value derived from the mask diameter; their
help texts give the rule.

### Stage controller

`simple_solve2D_controller.f90` emits one `refine2D` command line per stage.
Six stages, in two phases (stages 1-4 and 5-6):

| stage | search (`refine=prob_snhc`) | low-pass | references |
| --- | --- | --- | --- |
| 1 | `snhc_smpl` | `lpstart` | random, Gaussian-filtered (`gauref`), no centering, no shift search |
| 2 | `snhc_smpl` | ladder | the previous iteration's class averages, `ml_reg`, centering and shift search |
| 3-4 | sparse `prob_snhc` table | ladder | as stage 2 |
| 5 | sparse `prob_snhc` table | `lpstop` | as stage 2 |
| 6 | dense `prob` table | gold standard (even/odd) | as stage 2 |

- The ladder halves the distance to `lpstop` at each stage.
- `refine=prob` uses the dense table from stage 3 on; `refine=snhc_smpl` keeps
  the sparse table in the last stage too. Stages 1 and 2 are always
  `snhc_smpl`.
- Phase 1 stages end at fixed iterations (`minits = maxits`), spread over
  `extr_lim`; phase 2 stages run up to `nits_per_stage` iterations and end on
  convergence from the second on; a sample that needs more than
  `nits_per_stage` iterations to cover every particle fixes the count
  (`minits = maxits`).
- When updates are sampled, stages 5 and 6 (and the last stage) take
  `fillin=yes`: they keep iterating until every active particle has a class.
  Stage 1 samples but does not carry previous class sums forward.
- After the stages, a sampled run ends with one dense `refine=greedy` pass
  over every particle (the terminal pass) before the final class averages.

The log marks each stage (`>>> STAGE n WITH LP x A`, or `WITH GOLD STANDARD
E/O`) and its call (`>>> SOLVE2D ENTERING REFINE2D STAGE n`, `>>> SOLVE2D
RETURNED FROM REFINE2D STAGE n`). Each iteration prints one `>>> CONVERGED:
yes|no`, the decision after `minits`, `maxits` and the full-assignment
coverage, as 3D does.

### Seeded restart (`cls_init=prev`)

`cls_init=prev` re-enters the workflow from a previous 2D clustering held in
the project instead of a random start. It is a `solve2D`-only mode
(`refine2D` rejects it when it makes the starting references itself;
`solve2D_chunks` and the stream pass `cls_init` on and default it to `rand`).
Design record:
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
  project with no even/odd partition gets the project's partition
  (`ensure_field_halves`, before the sigma2 state is built; per micrograph
  group where the stack rows map to micrographs, by particle index otherwise;
  [refine3D_policy.md](../3D/refine3D_policy.md), section 5). Existing `eo`
  values are never touched. The partition's provenance lives on the stack rows,
  since `isthere('eo')` is false for `eo=0`.
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

`refine2D` is a public `simple_exec` program, as `refine3D` is: it runs in
shared memory, or as the distributed master when `nparts` is given (the
private `simple_private_exec prg=refine2D_distr` is the same commander). The
private `simple_private_exec prg=refine2D` is the worker of one partition,
which the master schedules; the memory estimator keys on that name. 2D has
the 3D structure: the matcher searches and accumulates, the strategy restores
the class averages and decides convergence.

`simple_commanders_solve2D.f90` owns:

- the `solve2D` entry point
- top-level defaults
- run orchestration across stages
- initial reference handling, including the `cls_init=prev` seed
  (validation/repair of the previous clustering, seed references, seed
  pass and its diagnostics); the seed partition itself is an `oris`
  operation (`reseed_classes`, `simple_oris_reshape.f90`)
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
- the full-assignment coverage request (`fillin=yes`) on the staged calls
- the seed-pass iteration index and command line of a `cls_init=prev` run
  (`solve2D_seed_pass_iter`, `set_cline_refine2D_seed_pass`)

`simple_refine2D_strategy.f90` owns:

- shared-memory versus distributed execution selection
- iteration control inside one `refine2D` invocation
- scheduler interaction
- probabilistic pre-alignment dispatch
- distributed worker scheduling
- the distributed master's own project: it runs `prob_align2D`
  (`prob_align2D_exec`), class-average assembly (`cavgassemble_exec`) and the
  sigma2 commit (`commit_sigma2_exec`) on its own `(params, build)`, each at
  the thread count it needs, the master's own restored after; a full-update
  iteration sets the master's `update_frac` to 1 for the alignment and restores
  it. Assembly updates the master's `cls2D`, `cls3D`, `out` and the particles'
  `res` in memory; the master writes those segments and the particle segment
  (Project segments, below) and does not reread the particle segment
- shared-memory restoration of the class averages from the sums the matcher
  accumulated (`restore_cavgs_inmem`: the carried sums committed, the averages
  restored and written, the class document, the `cls2D` and `cls3D` segments)
- the convergence decision of both modes (`check_refine2D_convergence`:
  `check_conv2D`, `minits`, `maxits`, the full-assignment coverage, the
  verdict and the progress file); it runs after the restoration, whose class
  document writes the particles' `res` that `check_conv2D` averages, and,
  in shared memory, before a full-update iteration restores `update_frac`
- run-finalization bookkeeping

`simple_strategy2D_matcher.f90` owns:

- particle-domain alignment/search
- reproduction of the probabilistic sampled subset when `prob_align2D` is active
- strategy-object selection
- sigma updates during Euclidean search
- accumulating the class sums from the batch images
- writing orientation updates
- writing distributed partial class-average sums when running as a worker;
  in shared memory it returns with the sums in the class averager and tells
  the strategy whether the iteration blends carried sums

The 2D matcher must preserve a single particle-stack read per batch in the
online alignment/restoration path. Batch construction should keep the
already-read particle images (or the Fourier-cropped images made from them,
below) for class-average restoration, and restoration should consume those
in-memory images after assignment in the same batch.

Do not split online class-average restoration into a second full particle pass
that re-reads image stacks to lower peak memory. Offline or terminal
class-average assembly commands may have their own explicit reads, but that is
separate from the matcher worker's online single-read contract.

Probabilistic table construction has a separate bounded-memory contract:

- `prob_tab2D` workers retain only thread-local compact candidates for the
  current particle batch and stream them to one partition file;
- `prob_align2D` constructs the global object only after workers complete;
- dense `refine=prob` uses a compact rectangular table, while sparse
  `refine=prob_snhc` uses a particle-oriented ragged candidate store;
- compact storage is materialized as the established `ptcl_ref` assignment only
  after global selection, so assignment-file semantics do not change.

Both table fills (`fill_tab_range` for `refine=prob`: the last stage and the
seed pass; `fill_tab_prob_snhc_range` for the `prob_snhc` stages) hold the
shift in the particle frame, as the 3D tables do (`refine3D_policy.md`,
probabilistic pre-step): the shift estimated at the particle's previous class
is converted once at the rotation its search ended at, the particle is
shifted once by it (`memoize_shifted_ptcl`, unless the shift is below
`SHERRSQ`), every class is scored against that memo with the bank kernel
(`gen_likelihood_val_bank`, which includes the angular Nyquist bin), a
refined class starts its shift search from the seed converted at its own
rotation and clamped into the search box, and an assigned, unrefined
candidate records the seed as it is (`record_seed_shift`). The particle loops
run with dynamic OpenMP scheduling (a particle is processed whole by one
thread, so only its random draws change). Each table part's bench record
(`BENCH/PROB_TAB2D_BENCH_PARTppp.txt`) carries the 3D tables' fill profile: the
thread-seconds of the seed, the score loop and the refinement, the fill
regions' wall time, the load balance and the shifted particle memos built.

The 2D part files (`dist_partN.dat`; one iteration's artifacts) have two
layouts, told apart by a format word in the header:

- the candidate stream (`refine=prob_snhc`, format -4): header (classes,
  particles, format word, candidates, chunks), particle indices, the
  candidate count of every particle, the seed table, then the chunks the
  worker threads flushed (particle indices and 24-byte candidates);
- the dense rows (`refine=prob`, format -5): header (classes, particles,
  format word, refined entries), particle indices, the seed table, one row
  of (distance, in-plane index) per particle and class (8 bytes per entry),
  then the refined entries (particle, class, shift) of the candidates whose
  shift was refined. A worker turns each thread's buffered candidates into
  rows when it flushes them: a particle is processed whole by one thread,
  so its candidates are consecutive in one buffer.

`prob_align2D` reads the parts in parallel (one task per part file, file
names built before the loop, a status per part raised after it): first every
part's head, which sizes the sparse store from the counts, then every part's
candidates or rows, once. Each part's size must be the size its header
implies, and every particle of the master must be in exactly one part
(counted with an atomic capture). The dense master keeps `rows(nclasses,
nptcls)` of 8 bytes and the refined entries in compressed-row form; the
assignment rebuilds a full candidate from them (`dense_candidate`). The
`PROB_ALIGN2D` report splits the reads (heads, candidates) and the
assignment (class lists: counts, fill, sort; the frontier loop; the leftover
particles).

Dense assignment (`ref_assign_likelihood`): every active class sorts the
particles by distance; each step draws among the classes' closest available
particles (the frontier). The frontier is kept, not rescanned: after an
assignment only the classes whose head was the assigned particle move on.

Sparse assignment (`ref_assign_sparse_likelihood`): each class's candidates
are gathered with their particles and distances during the fill pass and
sorted in parallel over the classes (the same unstable `hpsort` on the same
input, so ties fall as before); the head walk reads those arrays in order,
and each particle's peak count is computed once. The frontier loop keeps its
order, its one random draw per step and the swap-remove that decides ties
within a particle.

A distributed part with none of the sampled particles is valid: its table
part writes an empty table (header, no candidates; the readers skip its
rotation grid), and its matcher pass runs without batches, writing its
orientations, its sigma2 range and an empty class-sum contribution carrying
the references' centering offsets.

The search loop of the matcher runs under a dynamic OpenMP schedule (a
particle is searched whole by one thread, so only its random draws change).
In shared memory, without the particle cache, the probability table fills the
particle buffer (`doc/policies/particle_read_path_policy.md`) and the matcher
of the same iteration copies its batches from it, so a table iteration reads
each particle once. When the particles are cropped (below), the buffer holds
their records, the cropped particles, so the matcher neither reads nor
prepares them at the full box; otherwise it holds raw images and the
preparation runs twice.

`simple_cavgassemble_strategy.f90` (`cavgassemble_exec`), the classaverager
modules and the `make_cavgs`/`cavgassemble` commanders of
`simple_commanders_mkcavgs.f90` (command-line wrappers; what each writes is in
Project segments, below) own:

- explicit class-average assembly from partial sums (restoration reads the
  class FRCs for the even/odd low-resolution averaging before storing the new
  ones: in `refine2D` the previous iteration's, from the file, as the
  shared-memory matcher reads them, so the two modes restore alike;
  `make_cavgs` starts from none in both modes)
- merged/even/odd class-average output
- class-document generation
- class FRC output and project output metadata

The class sums of a batch (`cavger_update_sums`, the matcher's online sums and
the offline `make_cavgs` sums) are bucketed by class and half; the buckets run
under a dynamic schedule, the largest first, and each adds its particles in
batch order. The sums are therefore the same at any thread count. The work
images are sized by the thread team.

The matcher's class averages live at `box_crop`, so whenever `box_crop < box`
(and outside the stream's 2D pool; `crop_particles2D` decides it for the
matcher and the table alike) its online sums take the Fourier-cropped
particle the particle cache would serve (`particle_cache_policy.md`, section
3), made in memory from the full-size read of the batch: noise-normalized at
the full box, cropped, then tapered, padded and transformed at `box_crop`.
That replaces a transform of the padded full box of which only the
`box_crop` region was used. The alignment then prepares the particle exactly
as with `cache=yes`, and so does the probability table (`prob_tab2D`), so an
uncached run is byte-identical to a cached one. Both differ from the old
uncached full-box preparation by rounding: the cropped particle goes back to
real space and is transformed again, which is not bitwise exact. The class
sums change as the cached restoration does (statistical parity, not bit
equality). The final class averages of `solve2D` are restored at the full
box by `make_cavgs`.

Benchmark reports (`refine3D_policy.md`, section 11.1, lists the files; all in
`BENCH/`): `REFINE2D_BENCH_PARTppp.txt` for every partition, with the
thread-seconds and load balance of the search loop and of the class sums;
`REFINE2D_STRATEGY_BENCH.txt` for the whole iteration, setup to the
class-average JPEG (pre-step, matcher or scheduler, merge, class-average
assembly, the master's project writes and rereads, sigma2, convergence, star
file); `PROB_ALIGN2D_BENCH.txt` for the pre-step driver;
`CAVGSUMS_BENCH_PARTppp.txt` with the class sums' load balance.
`scripts/plot_refine3d_bench.py` turns them into the 2D critical path per
iteration.

### Project segments

Every process builds its own project from the file. `solve2D` runs `refine2D`
and `make_cavgs` as whole programs: it writes the segments they read before the
call and rereads the segments they wrote after it. What writes what, as built:

| Process | Writes |
| --- | --- |
| `solve2D` | `ptcl2D` at the start (fresh start: after the clustering is deleted; `cls_init=prev`: after the seed partition); `ptcl3D` at the end unless the project holds a valid 3D alignment; `out` after the final class averages. It rereads `cls2D` and `out` after every stage, and `ptcl2D` and `ptcl3D` before the end. |
| `refine2D`, shared memory | `cls2D` and `cls3D` every iteration (`restore_cavgs_inmem`); `ptcl2D` and `out` at the end of the invocation. With `prob_align2D`, `ptcl2D` is written before the table commander and again once the sample is drawn, and reread after. |
| `refine2D`, distributed master | `ptcl2D` once the sample is drawn (probabilistic modes), by the merge of the workers' documents, after class-average assembly (with `cls2D`, `cls3D` and `out`), and in the iteration's finalization; `out` at the end of the run. |
| `refine2D` worker | no project segment: its alignment document, sigma2 part, `cavg_contrib_part<N>.bin`, `dist_part*` and `assignment` files. |
| `make_cavgs`, shared memory | `ptcl2D` (class remap, halves), then `cls2D` and `out`. |
| `make_cavgs`, distributed | the master, before the workers start: `ptcl2D` (class remap, expansion, time-series split, halves); the workers: no project segment; the assembly (`run_cavgassemble`): the whole project. |
| `cavgassemble` | the whole project. |
| the sigma2 rebuild (`ensure_sigma2_for_iteration`) | `projinfo` (the registered sigma2 state); the caller rereads `projinfo`. |

`ensure_field_halves` runs on the master only (a worker checks that the halves
exist); when it labels the halves it writes `stk`, `ptcl2D` and `ptcl3D`.

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
  selection still follows the normal sampled-update path; on a fresh start
  (`startit=1`) it selects the update-count sampler instead of the sticky
  subset of `solve2D`'s first stage, which could never cover every particle
- if any staged update used `update_frac`, `solve2D` runs a terminal
  dense `refine=greedy` all-particle `refine2D` pass with `update_frac` and
  `fillin` disabled, refreshing class, in-plane, and shift parameters before
  final class-average generation

The restoration model is class-local: each class carries forward its previous sums under the population rule of the [class-average state note](../../refactoring_notes/completed/class_average_and_reconstruct3d_partials_refactoring.md) (Section 4.1). The carried set records `M(c)`, the population its sums represent; the assembly owner counts `N(c)` (active, updated rows of the class) and `n(c)` (those sampled this round) on the merged project and blends `new = current + w * shift(previous)` with `w = (N - n) / M`, recording `M <- N`. The carried mass therefore equals the represented population whatever joined or left the class; without a population change `w = 1 - f`. The rule keeps the mass right, not the membership: the carried sums are an approximation (a stochastic recurrence), not an exact sum over the current class members, and old contributions are removed in proportion, not particle by particle. Shared memory is the reference behaviour; distributed runs give the same result, independent of `nparts`.

### Continuous in-plane policy

`inpl_cont` is internal since release 4 (2026-10-03): it is not a
command-line key, and every `refine2D` child (probabilistic staged calls, the
final staged invocation, and the terminal dense all-particle refresh) runs
with the default `yes`. The historical alternating shift/discrete-angle
callback route (`no`) is no longer selectable in 2D. `yes` replaces every callback-based angle/shift optimization with the
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
[refine3D policy](../3D/refine3D_policy.md).

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
- `cavg_state.bin`: the carried class sums of the run, with `M(c)` per class; no
  part number. Written by an assembly owner (class-average assembly on the
  distributed master, the strategy's shared-memory restoration); the
  distributed `make_cavgs` assembly also publishes its sums there, at its own
  box, and the stream pool pads the set to a larger box. Nothing removes it:
  the queue cleanup (`qsys_cleanup`) deletes the 2D class sums only when asked
  to, and no caller asks
- `cavgs_iterNNN.mrc`, `cavgs_iterNNN_even.mrc`, `cavgs_iterNNN_odd.mrc`
- `FRCS_FILE`
- `sigma2_state.bin`: the project's one canonical sigma2 store. Before the
  parts of an iteration write, the master validates the committed store and
  clears what another transaction of the next generation left
  (`prepare_canonical_sigma_update`); the parts write their ranges as a
  candidate, and `commit_sigma2` publishes it as the next generation
- `ptcl2D`, `cls2D`, `cls3D`, and `out` project segments

Partition-local probabilistic assignment/dist files and the class-sum contributions are per-iteration artifacts: the master removes them before the next distributed iteration writes new ones, and the assembly owner deletes the contributions after the blend. Workers never read carried state. The carried set `cavg_state.bin` is the only durable source of class-average carry-over: the owner (class-average assembly on the distributed master, the strategy's shared-memory restoration) reads it, applies the class-centering shifts the workers report once, blends, and publishes the new set through a temporary name and a rename. When an iteration would blend but the carried set is missing, unreadable, or disagrees with the run on class count, `box_crop` or `smpd_crop`, the master runs that iteration as a full update (every particle, no carry-over). Legacy `cavgs_*_partN.mrc` / `ctfsqsums_*_partN.mrc` files are never read; a run continued from such a directory starts with one full update. Changing `nparts` or the execution mode across `continue=yes` needs nothing else.

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
- Does `inpl_cont=yes`, the value every 2D run uses, avoid the callback route
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
- Do not use the staged full-assignment coverage (`fillin=yes`) as a substitute
  for the terminal dense all-particle refresh when sampled solve2D updates
  were active.
- Do not reuse stale assignment files as valid current-iteration inputs.
- Do not re-read particle stacks in the online matcher/restoration path when the
  raw batch images are already available.
- Do not let workers read or write the carried class sums, and do not reintroduce partition-shaped carry-over: workers accumulate the current sample from zero; only the assembly owner blends and publishes `cavg_state.bin`.
- Do not add fractional in-plane coordinates to probabilistic assignment
  artifacts; rerun the joint optimizer after final assignment instead.
- Do not invoke the legacy angle callback from any `inpl_cont=yes` failure or
  no-improvement path.
