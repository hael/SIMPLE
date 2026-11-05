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
4. run staged `refine2D` iterations (the later stages probabilistic: each
   particle's class drawn from its own row inside the matcher's batch pipeline)
5. update particle class, in-plane, shift, sampled, and update-count state
6. restore class averages through shared-memory or distributed class-average pathways
8. when updates are sampled, keep the last stages iterating (`assign_all=yes`) until every active particle has a class
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
assignment at every stage. For `refine=prob_neigh` (2D has only
`prob_neigh_mode=snhc`; the mode was `refine=prob_snhc` before 2026-11-03,
when 2D took the 3D names), intermediate stages use sparse probabilistic SNHC and the final staged invocation uses dense
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
- `refine=prob_neigh` (`prob_neigh_mode=snhc`)
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

| stage | search (`refine=prob_neigh`) | low-pass | references |
| --- | --- | --- | --- |
| 1 | `snhc_smpl` (dense `prob` on the evidence over the disk under the polar shift map) | `lpstart` | random, Gaussian-filtered (`gauref`), no centering, no shift search |
| 2 | `snhc_smpl` | ladder | the previous iteration's class averages, `ml_reg`, centering and shift search |
| 3-4 | sparse `prob_neigh` table | ladder | as stage 2 |
| 5 | sparse `prob_neigh` table | `lpstop` | as stage 2 |
| 6 | dense `prob` table | gold standard (even/odd) | as stage 2 |

- The ladder halves the distance to `lpstop` at each stage.
- `refine=prob` uses the dense table from stage 3 on; `refine=snhc_smpl` keeps
  the sparse table in the last stage too. Stages 1 and 2 are always
  `snhc_smpl`, except stage 2 under the polar shift map (below).
- Phase 1 stages end at fixed iterations (`minits = maxits`), spread over
  `extr_lim`; phase 2 stages run up to `nits_per_stage` iterations and end on
  convergence from the second on; a sample that needs more than
  `nits_per_stage` iterations to cover every particle fixes the count
  (`minits = maxits`). The controller plans the stages in absolute
  iterations and emits each as `refine2D` run lengths from the stage's
  `startit`: `refine2D` counts `minits` and `maxits` from `startit`, as
  `refine3D` does (2026-11-03; before, the 2D keys were absolute iterations).
- When updates are sampled, stages 5 and 6 (and the last stage) take
  `assign_all=yes` (`fillin=yes` before 2026-11-03; `fillin` is 3D's
  periodic sampling of the never-updated particles): they keep iterating until
  every active particle has a class.
  Stage 1 samples but does not carry previous class sums forward.
- After the stages, a sampled run ends with one dense `refine=greedy` pass
  over every particle (the terminal pass) before the final class averages.

**Stage omission (2026-11-03).** A stage whose planned limit the class
averages already resolve past is omitted, as `solve3D` omits a stage behind
its FSC=0.5 resolution (`doc/policies/3D/solve3D_policy.md`). The evidence is
the class FRCs (`frcs.bin`, the out segment's `frc2D`) at the end of every
executed stage: per selected populated class the FRC=0.5 crossing
(`class_frcs%estimate_res`), and from them the resolution the classes of
`FRC05_SKIP_PTCL_FRAC` (0.5) of the particles reach, the classes best first
and their populations cumulated (`frc05_ptcl_quantile`,
`solve2D_class_frc05_quantile`): the median particle's class resolution,
which a few overfit or empty classes cannot move, where the 3D rule has one
map and one FSC. The log prints it after every stage (`>>> STAGE n CLASS
FRC=0.5 RESOLUTION REACHED BY 50% OF THE PARTICLES`) with the fraction of the
particles in classes past the next planned limit. Stage k is omitted
(`solve2D_stage_behind_frc05`, a `cycle` in the marching loop, `>>> SOLVE2D
OMITTING STAGE`) when that resolution at the ends of the last two executed
stages is finer than `stage_parms(k)%lp`: the estimate must have held across
a stage, the stability guard the many per-class estimates need. Never
omitted: the entry stage, the first `SKIP_MIN_STAGE` (2) stages (the second
is the centring stage, whose class averages estimate an automatic run's
mask), the last stage of the invocation (the gold-standard stage has no
planned limit; a checkpoint segment's last stage too) and a stage without a
planned limit; the first two stages after a checkpoint resume or a
`cls_init=prev` seed pass only record the evidence. Phase 1 ends at fixed
iterations because the extremal schedule runs over `extr_lim`: an omitted
phase-1 stage before the last advances the ladder, the phase-1 stage that
follows taking its iterations at its own limit, while an omitted last phase-1
stage (`PHASES(1)`) drops them and ends the extremal schedule there, at the
small neighbourhood fraction its decay has reached (logged either way). A
phase-2 stage's iterations are dropped. On the beta-galactosidase fixture
(4463 particles, 20 classes, run 1428 of 2026-11-03) the median particle's
class resolution tracked the ladder within one rung and no stage was omitted;
the rule is meant for the easy data whose classes outrun the march.

The log marks each stage (`>>> STAGE n WITH LP x A`, or `WITH GOLD STANDARD
E/O`) and its call (`>>> SOLVE2D ENTERING REFINE2D STAGE n`, `>>> SOLVE2D
RETURNED FROM REFINE2D STAGE n`). Each iteration prints one `>>> CONVERGED:
yes|no`, the decision after `minits`, `maxits` and the full-assignment
coverage, as 3D does.

#### The polar shift map's dynamic centring (`shift_sigma` > 0)

`shift_sigma` (the width in A of the Gaussian prior on a particle's displacement
from its pick, `doc/implementation_notes/completed/polar_shift_search.md` and
`polar_shift_search_dynamic_centring.md`) is the parameter's default, 10 A
(Hans, 2026-11-02, across the board), so every workflow has the key where a mode
can carry it: a polar mode with a shift search under objfun=euclid, a
probability table (`prob`, `prob_neigh`) or, since ruling 20
(2026-11-03), a classic mode (`snhc`, `snhc_smpl`, `snhc_smpl_many`, `greedy`;
3D `shc`, `neigh`, `greedy`, `greedy_inpl`); elsewhere (`tseries`, the `inpl`
and `greedy_smpl` modes, a stage without shifts) the default is silently off,
and under objfun=cc an explicit key is refused.
`shift_sigma=0` is the window search in every stage, the control. Two flags
carry this: `l_shift_map`, the key's presence (a workflow's master passes it to
its stages), and `l_centring`, a run that centres (the key in such a mode with
a shift search), which the searches, the preparation and the consumers read.

The classic modes (ruling 20) take the table's structure, not a seed: under the
key no shift search at the previous class (`inpl_srch_first`) in any state count,
the candidates scored at the stored shift as the mode has it, and the screen and
polish of the tables (`simple_shiftmap_srch`: a map per candidate of the pool,
the disk or the window, the best cells polished with `minimize_in_cell`) over
every scored class and the previous one, which competes at its stored-shift
score: one polished for the committed pose, the peaks for `snhc_smpl`. Under the
key `snhc` keeps every class it scores. (A first build mapped only the class SHC
committed to for a particle in its window: a wrong class's map moved the particle
to that class's optimum, and the 2D fixture fell well below the window search.)
A particle that needs the disk is searched greedily (every class scored), as the
sparse table scores every class for one; the written class carries its own
mapped shift, and `assign_ori` records the shift state. The matcher builds the
maps once per run (`setup_shiftmap_srch`). The stream's pool (`snhc_smpl`, its
appended particles greedy) is centred this way.
The continuous in-plane route's joint domain contains the particle's box and the
incoming shift (the refine3D policy), so a disk shift is polished, not clipped.

The centring is a per-particle state, not a stage (ruling 19, 2026-11-02; the
sampling strategies inject particles no stage would otherwise centre). Each
particle carries two numbers the consumer writes after its polish
(`record_shift_state`): `shincr`, the move of the iteration, and `shsig`, its
prior width for its next search, in pixels. In every table fill from stage 2 on
(stage 1 is the SNHC window stage, which has no table):

- a sampled particle takes the disk when `shsig` is the sentinel (below zero:
  not yet centred in this workflow) or when three of its prior widths exceed
  the stage's window `trs` (`shift_needs_disk`), and the window otherwise;
- the disk: the coarse map over the disk of 3 sigma (`disk_radius_px`, bounded
  by the box; `SHIFTMAP_SCREEN`) under the prior at the particle's first
  centring, and without it after (a re-centring is the data's alone: under the
  draw's tempering, `prob_beta` 0.3, the untempered prior weighs more than three
  times its physical weight, and when it sat at the pick it pulled correctly
  centred particles back towards their picks at every re-centring, which cost
  the solve2D fixture 56 points). The prior moves with the particle (ruling 22,
  2026-11-03): it is on the cell's shift, the increment from the stored shift
  the particle is prepared at, which the class centring moves too, at the classes
  holding 99 % of the particle's posterior at the
  origin (`screen_pool_size`, at least twice the refined count; in the sparse
  `prob_neigh` fill a disk particle scores every class first), each keeping its
  coarse optimum, and the best classes polished at full resolution in a box of
  one map step about that optimum in the disk's own search object
  (`pftc_shsrch_grad%minimize_in_cell`: the L-BFGS from the cell; the better of
  the polish and the cell at its best rotation is the pose if it beats the
  origin's row, else that row stands). The rows are polished losses;
- the window: the same screen and polish over the coarse window map
  (`SHIFTMAP_WINDOW`, radius `trs`, no prior), the stage's search object;
- a disk particle is prepared under the box's mask (`ptcl_on_disk`,
  `prepimg4align`: the box's noise mask and soft mask radius, `lmsk_box`,
  `msk_box_crop`), since the molecular mask at the pick would cut it; the
  others under the stage's mask. The preparation records its choice per particle
  (`builder%note_obs_mask`), and every later normalisation of the particle in the
  same process reads it (`obs_on_box`: the class sums from full-size images, the
  3D rereads `prep_rec_observation_ptcl`), so a class sum or a reconstruction is
  the same whether the particle's record was kept or it was read again (the
  review's M2, 2026-11-03);
- the width: every coarse map's draw measures the standard deviation of the
  shift under its posterior at the physical temperature, about the mean, in the
  particle's axes (`shiftmap_draw`); it travels with the candidate
  (`prob_candidate%sdev`, `prob_refined_entry%sdev`, `ptcl_ref%sdev`) to the
  consumer, which stores as the next prior width the width of the assigned
  pose (the previous width when the pose carries none), at least the radius
  searched when the move reached that radius's bound, within a tenth of it or
  half a map step where the map is finer (`next_shift_width`: a posterior
  truncated at the bound says nothing of the mass beyond). The width chooses
  the next box; it never enters the likelihood, whose only prior is the pick's;
- the sentinel is set by the refine2D masters at a run's first iteration
  (every workflow's first stage starts there) and by the stream's pool for the
  particles it appends (`append_sets`); the master prints the disk set it is
  about to search by cause (not yet centred, wide) and the widths' median and
  90th percentile over the centred particles after its sample draw
  (`print_disk_set`; a distributed part prints none).

No stage ends on the centring: the stages exit by their own rules. An
automatic run (`mskdiam=0`) runs stage 2 under the box's mask so that its class
averages can estimate the mask diameter at its exit
(`solve2D_mskdiam_estimate`, adopted at the exits of stage 1 and stage 2 with
the low-pass ladder re-derived, `rederive_lplims`); a user's `mskdiam` is the
stages' mask from stage 2 on, the disk particles excepted.

Stage 1 is the plain SNHC window stage under the key too. Two dense-table
variants of stage 1 over the disk of 3 sigma (the class from the evidence, the
drawn cell assigned or not) were measured and dropped 2026-11-01: both made
the centring stage collapse the populated classes at its entry (90 to 50 with
no shift assigned, 38 at the end; the SNHC stage 1 gives 79 at the entry and
57 at the end, the control 65), the per-iteration counts from the runs'
`clusters2D_iterNNN.star`.
Under the key every table scores the classes at the stored shift, without the
previous class's seed and its shifted-particle bank, and polishes from a coarse
map's optimum (Phase 4 of the note), so the 5 x 5 lattice at every rotation
and the seed retire wherever the key runs; `shift_sigma=0` keeps the seed and
the lattice. The terminal greedy pass runs the window search (`shift_sigma=0` on
its command line; since ruling 20 the greedy mode could carry the key, which this
pass was validated without). `cls_init=prev` enters at stage 3 and runs no
centring stage; the stream's checkpointing is refused with the key. With the
key off every command line is as before.

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
  (`update_frac`, `nsample`, `fillin` and `assign_all` deleted,
  `extr_iter=extr_lim+1`)
  at the low-pass limit of `PROBREFINE_STAGE`. It uses that stage's ordinary
  shift window: `shift_sigma` and `trs_window` are deleted even when the
  parent solve2D uses the default shift map, since seeded restarts skip the
  centring stage. The pass runs as the
  iteration at which stage `PROBREFINE_STAGE-1` would have ended
  (`solve2D_seed_pass_iter`), so the stages from `PROBREFINE_STAGE` to
  the terminal greedy pass run exactly as on an unseeded run: same limits,
  refine policy, iteration counts, sampling and fractional restore. The
  pass must not be a fresh start (`startit > 1`): a fresh start zeroes the
  shifts the matcher's table fill would seed from.
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
which the master schedules; the memory estimator keys on that name. The
partitions balance the active particles (`state > 0`), as in 3D (2026-11-03;
before, equal particle ranges). 2D has the
3D structure: the matcher searches and accumulates, the strategy restores
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
- the full-assignment coverage request (`assign_all=yes`) on the staged calls
- the seed-pass iteration index and command line of a `cls_init=prev` run
  (`solve2D_seed_pass_iter`, `set_cline_refine2D_seed_pass`)

`simple_refine2D_strategy.f90` owns:

- shared-memory versus distributed execution selection
- iteration control inside one `refine2D` invocation
- scheduler interaction
- the iteration's sample draw under the probabilistic and the sampled modes
  (`sample_ptcls4prob2D`, `sample_ptcls4iter2D`, keyed on the run's seed, the
  pass, the iteration and the particle; written to the project before the
  parts, which reproduce it)
- distributed worker scheduling: one scheduling phase per iteration; the
  probabilistic parts publish readout records that the master merges in part
  order and prints
- the distributed master's own project: it draws the sample, runs
  class-average assembly (`cavgassemble_exec`) and the sigma2 commit
  (`commit_sigma2_exec`) on its own `(params, build)`, each at the thread count
  it needs, the master's own restored after; a full-update iteration draws
  every particle. Assembly updates the master's `cls2D`, `cls3D`, `out` and the
  particles' `res` in memory; the master writes those segments and the
  particle segment (Project segments, below) and does not reread the particle
  segment
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

- particle-domain alignment/search; under `refine=prob` and `prob_neigh` the
  batch pipeline (`prob_align_batch_pipeline.md`, 2026-11-04): the batch's
  rows or candidates filled, each particle's class drawn from its own row,
  then the polish, residual and class sums, the same body in shared memory and
  in a distributed part
- reproduction of the master's sampled subset
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

The probability table has a bounded-memory contract: the fill keeps
thread-local compact candidates for the current batch alone, the assignment
turns them into the batch's dense rows (`refine=prob`, one entry per class) or
ragged candidate store (`refine=prob_neigh`) and releases them after the
batch's draws; the assignment map (`ptcl_ref`) is the part's only state across
batches. No partition file and no global table exist (retired 2026-11-05).

Both table fills (`fill_tab_range` for `refine=prob`: the last stage and the
seed pass; `fill_tab_prob_neigh_range` for the `prob_neigh` stages) hold the
shift in the particle frame, as the 3D tables do (`refine3D_policy.md`,
section 8): the shift estimated at the particle's previous class
is converted once at the rotation its search ended at, the particle is
shifted once by it (`memoize_shifted_ptcl`, unless the shift is below
`SHERRSQ`), every class is scored against that memo with the bank kernel
(`gen_posterior_val_bank`, which includes the angular Nyquist bin), a
refined class starts its shift search from the seed converted at its own
rotation and clamped into the search box, and an assigned, unrefined
candidate records the seed as it is (`record_seed_shift`). The particle loops
run with dynamic OpenMP scheduling (a particle is processed whole by one
thread; its draws are keyed on the particle, so the schedule changes none).
Each part's bench record (`BENCH/REFINE2D_BENCH_PARTppp.txt`, the `match2D fill
...` lines) carries the fill profile: the thread-seconds of the seed, the score
loop and the refinement, the fill regions' wall time, the load balance and the
shifted particle memos built, beside the wall time of the table fill and the
assignment.

The 2D part files (`dist_partN.dat`, the candidate stream of `prob_neigh` and
the dense rows of `prob`), their parallel read by `prob_align2D` and the
`PROB_ALIGN2D` report were retired on 2026-11-05 with the `prob_tab2D` worker:
a particle's candidates stay in the thread buffer that filled them until the
batch's assignment.

The assignment is one draw per particle from its tempered row posterior
(2026-10-23; the tables' only assignment since 2026-10-24, when the
class-balanced frontier was deleted; the constructor requires
`objfun=euclid`, the loss the posterior is defined by). The rule of the 3D
tables (`refine3D_policy.md`,
`doc/implementation_notes/completed/pose_marginalization.md`) on the 2D tables,
for `refine=prob` and `refine=prob_neigh`: no class is forced to be
populated. At the fill, every class's in-plane rotation is one draw from the
posterior over all grid rotations of the class against the seeded particle at
the particle's temperature T_i = pi wsqsum_i / (2 pftsz) (`gen_posterior_val`,
its bank form; `score_seeded_class`, the shift-free direct call), no angular
threshold; the fill keeps T_i with the row. At the batch's assignment
(`assign_batch`), each particle's class is one draw from the posterior over its
row, p proportional to exp(-T_i (L - L_min)): over the active classes of its
dense row (`assign_rows_dense`, the drawn class's `row_candidate`, `frac` 100),
or over its valid evaluated candidates (`assign_rows_sparse`, an active class
with an in-plane index; the drawn candidate with its shift, the peak count
the candidates with an in-plane index and the fraction their share of the
classes). A particle without a valid distance takes the first active class
(dense) or its best candidate (sparse); one without a temperature its argmin,
with a warning. The draw is keyed on the particle, so shared-memory and
distributed runs draw alike; the master reports the fraction of particles
drawn at their row mode from the parts' readout records. The tempering factor
`PROB_BETA` (0.3, `simple_defs`; for small membrane proteins; below one
hotter) multiplies the temperature of every draw, the fills' in-plane draws
and the class draws as of the 3D direction draws, so below one they draw
from the tempered posterior p^beta (`refine3D_policy.md`). The workflows
carry the iteration's value as `prob_beta` and the tail moves it; since
2026-10-29 neither `prob_beta` nor `prob_beta_end` is a UI input of
`refine2D`, `solve2D` or `solve2D_chunks`. The 2D
tables log the n_eff readout of the 3D ones (`POSTERIOR DRAW N_EFF`, the mean
effective number of classes a particle draws among at the draw's factor and
at beta 1; `refine3D_policy.md`); the n_eff control that briefly replaced
`prob_beta` (2026-10-27 to 2026-10-29) is gone from both routes, and the
weighted insertion (`topk`) is 3D only.

**The annealing tail (2026-10-27, as the 3D one).** A `refine=prob|prob_neigh`
run with `prob_beta_end` searches at `prob_beta` until its class overlap
converges past `minits` (`check_conv2D`'s own rule), then `exec_refine2D`
ramps the draws geometrically to `prob_beta_end` over the next
`PROB_BETA_TAIL` (5) iterations and the run ends (`>>> CLASS OVERLAP
CONVERGED: annealing ...`, `>>> ANNEALING TAIL COMPLETE`); a budget too short
for that starts the tail so that it completes at the run's last iteration,
`startit + maxits - 1`, as in 3D (`>>> BUDGET END`). The tail's
state lives in the strategy's convergence object (the shared-memory strategy
keeps one across iterations since this step, as the 3D one). `solve2D` sets
`prob_beta_end=1` on its last stage unless the user set it, so the last
probabilistic stage ends with draws from the physical posterior before the
terminal greedy pass and the final class averages; the earlier stages run
flat, as before. Until this step the 2D runs never annealed: a run at 0.3
wrote its last class averages from draws at 0.3, in expectation a soft
membership among the best few classes.

The posterior sampling is the tables' (`refine=prob|prob_neigh`, the modes
of the later stages); the SHC strategies of the early stages (`snhc_smpl`,
`snhc_smpl_many`, `inpl_smpl`, `greedy_smpl`) keep their own samplers. A
day's migration of those onto the posterior primitive (2026-10-23) was
taken back on 2026-10-24.

A distributed part with none of the sampled particles is valid: its matcher
pass runs without batches, writing its orientations, its sigma2 range, an
empty readout record and an empty class-sum contribution carrying the
references' centering offsets.

The search loop of the matcher runs under a dynamic OpenMP schedule (a
particle is searched whole by one thread; the draws are keyed on the particle,
so the schedule changes none). The pipeline reads and prepares each batch once
for the fill, the draw and the search; the particle buffer's hand-over from a
table pass to the search pass went with the table pass
(`doc/policies/particle_read_path_policy.md`).

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
particle, made in memory from the full-size read of the batch: noise-normalized
at the full box, cropped, then tapered, padded and transformed at `box_crop`.
That replaces a transform of the padded full box of which only the
`box_crop` region was used. The alignment then prepares the cropped particle
(`prepimg4align_cropped`), the table fill's and the search's one preparation.
It differs from the old uncached full-box preparation by rounding: the cropped
particle goes back to real space and is transformed again, which is not
bitwise exact (statistical parity, not bit equality). The final class averages
of `solve2D` are restored at the full box by `make_cavgs`. The downscaled
particle disk cache was retired 2026-11-02 (`particle_cache_policy.md`).

Benchmark reports (`refine3D_policy.md`, section 11.1, lists the files; all in
`BENCH/`): `REFINE2D_BENCH_PARTppp.txt` for every partition, with the
thread-seconds and load balance of the search loop and of the class sums;
`REFINE2D_STRATEGY_BENCH.txt` for the whole iteration, setup to the
class-average JPEG (the sample draw, matcher or scheduler, merge, class-average
assembly, the master's project writes and rereads, sigma2, convergence, star
file); `CAVGSUMS_BENCH_PARTppp.txt` with the class sums' load balance.
`scripts/plot_refine3d_bench.py` turns them into the 2D critical path per
iteration.

### Project segments

Every process builds its own project from the file. `solve2D` runs `refine2D`
and `make_cavgs` as whole programs: it writes the segments they read before the
call and rereads the segments they wrote after it. What writes what, as built:

| Process | Writes |
| --- | --- |
| `solve2D` | `ptcl2D` at the start (fresh start: after the clustering is deleted; `cls_init=prev`: after the seed partition); `ptcl3D` at the end unless the project holds a valid 3D alignment; `out` after the final class averages. It rereads `cls2D` and `out` after every stage, and `ptcl2D` and `ptcl3D` before the end. |
| `refine2D`, shared memory | `cls2D` and `cls3D` every iteration (`restore_cavgs_inmem`); `ptcl2D` and `out` at the end of the invocation. |
| `refine2D`, distributed master | `ptcl2D` once the sample is drawn (the probabilistic modes and the sampled modes under fractional updates), by the merge of the workers' documents, after class-average assembly (with `cls2D`, `cls3D` and `out`), and in the iteration's finalization; `out` at the end of the run. |
| `refine2D` worker | no project segment: its alignment document, sigma2 part, `cavg_contrib_part<N>.bin` and, under the probabilistic modes, its readout record `prob_readout_part<N>.txt`. |
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
- the stages preserve sample-once-and-reuse: the master (or the shared-memory
  matcher) chooses the subset with a keyed draw, and `refine2D_exec`
  reproduces it rather than resampling
- every sampled stage, probabilistic or not, draws on the distributed master
  (`master_draws_sample2D`, 2026-11-03 for the classic sampled stages, 2026-11-04
  for the probabilistic ones, as 3D's `master_draws_sample3D`): over every
  particle, written to the project before the workers start, which reproduce
  it (`sample4update_reprod`), so the sample does not depend on `nparts`; the
  stream's pool (`stream2d=yes`) runs the classic routes and keeps its own draw
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
- staged `assign_all=yes` (`fillin=yes` before 2026-11-03) acts as a full-assignment coverage guard: it
  keeps iterating until active particles have assignments, while particle
  selection still follows the normal sampled-update path; on a fresh start
  (`startit=1`) it selects the update-count sampler instead of the sticky
  subset of `solve2D`'s first stage, which could never cover every particle
- if any staged update used `update_frac`, `solve2D` runs a terminal
  dense `refine=greedy` all-particle `refine2D` pass with `update_frac`,
  `fillin` and `assign_all` disabled, refreshing class, in-plane, and shift parameters before
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

- `startit` is the stage/invocation start; `minits` and `maxits` count from it
- `which_iter` is the current iteration
- `extr_iter` tracks the 2D extrapolation/search schedule
- `endit` is written after an invocation finishes and is consumed by the next stage setup

Do not collapse these counters into one another. Child command lines, including fill-in, must preserve the distinction between stage start and current iteration.

## 6. Artifact and Handoff Policy

Stable 2D workflow artifacts include:

- `prob_readout_part<N>.txt`: one probabilistic part's readout record (counts
  and sums), merged in part order and deleted by the master after the schedule
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

The parts' readout records and the class-sum contributions are per-iteration artifacts: the master removes them before the next distributed iteration writes new ones, and the assembly owner deletes the contributions after the blend. Workers never read carried state. The carried set `cavg_state.bin` is the only durable source of class-average carry-over: the owner (class-average assembly on the distributed master, the strategy's shared-memory restoration) reads it, applies the class-centering shifts the workers report once, blends, and publishes the new set through a temporary name and a rename. When an iteration would blend but the carried set is missing, unreadable, or disagrees with the run on class count, `box_crop` or `smpd_crop`, the master runs that iteration as a full update (every particle, no carry-over). Legacy `cavgs_*_partN.mrc` / `ctfsqsums_*_partN.mrc` files are never read; a run continued from such a directory starts with one full update. Changing `nparts` or the execution mode across `continue=yes` needs nothing else.

## 7. Review Checklist

For any `solve2D` or `refine2D` change, check:

- Does the command layer remain mostly orchestration?
- Is stage policy in the controller rather than scattered through matcher or strategy code?
- Does probabilistic 2D sample once and then reproduce the same subset?
- Do shared-memory and distributed paths use the same scientific workflow?
- Are class-average assembly/restoration responsibilities explicit?
- Does online class-average restoration reuse the matcher batch images instead
  of introducing a second particle-stack read?
- Does the probability table stay batch-bounded in the matcher (the batch's
  rows or candidates released after its draws, no table across batches, no
  file)?
- Are stale distributed handoffs (the readout records and the class-sum contribution files) removed, while the carried set `cavg_state.bin` is left to its owner?
- Are `startit`, `which_iter`, `extr_iter`, and `endit` semantics preserved?
- Is `assign_all=yes` treated as a full-assignment coverage guard unless the
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
- Do not let the master's sample draw and the parts' matcher update use different particle subsets: the parts reproduce, never redraw.
- Do not make distributed-only class-average assembly semantics diverge from shared-memory scientific behavior.
- Do not describe `assign_all=yes` as missing-only assignment while it still uses
  the normal sampled-update path.
- Do not use the staged full-assignment coverage (`assign_all=yes`) as a substitute
  for the terminal dense all-particle refresh when sampled solve2D updates
  were active.
- Do not re-read particle stacks in the online matcher/restoration path when the
  raw batch images are already available.
- Do not let workers read or write the carried class sums, and do not reintroduce partition-shaped carry-over: workers accumulate the current sample from zero; only the assembly owner blends and publishes `cavg_state.bin`.
- Do not add fractional in-plane coordinates to the probability tables' rows;
  rerun the joint optimizer after final assignment instead.
- Do not invoke the legacy angle callback from any `inpl_cont=yes` failure or
  no-improvement path.
