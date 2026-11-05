# Importance Sampling and Fractional Update Policy

This document records durable workflow contracts for sampled particle updates,
probabilistic candidate sampling, fractional class-average restoration, and
trailing reconstruction in `solve2D`, `refine2D`, `solve3D`, and
`refine3D`. It is policy, not a line-by-line implementation map.

## 1. Core Model

SIMPLE has two sampling layers that must remain separate:

1. outer fractional-update sampling chooses which particles participate in the
   current iteration
2. inner importance sampling chooses which reference, orientation, or in-plane
   candidates are explored for those participating particles

The outer subset is recorded in the project through `sampled` and `updatecnt`.
Downstream restoration and reconstruction consume that recorded state. They must
not infer participation from the nominal command-line `update_frac` alone.

Probabilistic alignment is a sample-once-and-reuse path: the master (or the
shared-memory matcher) chooses the outer subset with a keyed draw and writes
it, and every part's matcher reproduces it once for its batch pipeline (the
table fill, the draw and the hard particle update in one pass).

## 2. Ownership

`simple_commanders_solve2D.f90` owns `solve2D` orchestration: defaults,
stage execution and final class-average generation; fill-in is requested on
the staged `refine2D` calls (`fillin=yes`, set by the controller).

`simple_solve2D_controller.f90` owns the 2D stage policy: `NSAMPLE_DEFAULT_2D`,
`nsample` override handling, stage-local `update_frac`, search-mode transitions,
and the rule that stage 1 may sample particles without fractionally restoring
previous class averages.

`simple_commanders_solve3D.f90` and `simple_solve3D_controller.f90` own 3D
stage scheduling: dynamic `update_frac`, `fillin`, `frac_best`, `balance`,
`trail_rec`, and transitions between early `prob_neigh` modes, `prob`, and
late `prob_neigh`.

`simple_matcher_smpl_and_lplims.f90` owns the shared outer subset-selection
helpers for 2D and 3D. This is where full update, random sampling,
update-count-biased sampling, nested-quota sampling over class units
(`balance=class|cavg`), cohort rescoring, fill-in sampling, and subset
reproduction are dispatched.

`simple_view_partition_sampling.f90` owns the sampling units: one unit per
selected 2D class, the class-average groups of `balance=cavg`, and the
coverage report printed before the first stage.

`simple_oris_sampling.f90` and `simple_oris_getters.f90` own the bookkeeping:
`sampled`, `updatecnt`, the nested quota (`sample4update_class`,
`class_sample_quotas`, `class_sample_sweep`), exact subset reproduction and
cohort rescoring, global realized update fraction, and class-local realized
update fractions.

The outer subset of a probabilistic iteration is drawn once, by the master on
its own project (`sample_ptcls4prob3D`, `sample_ptcls4prob2D` in
`simple_matcher_smpl_and_lplims.f90`; the shared-memory matcher calls the same
helpers), keyed on the run's seed, the pass, the iteration and the particle,
and written before the parts start. The `prob_align` step that sampled,
scheduled table workers, aggregated their outputs and wrote an assignment
artifact was retired on 2026-11-05.

`simple_eul_prob_tab*.f90` owns inner candidate importance sampling. These
modules may sample references, orientations, neighbors, or in-plane candidates
inside the active particle subset, but they must not choose a new particle
subset.

`simple_strategy2D_matcher.f90` and `simple_strategy3D_matcher.f90` own
particle-domain search on the active subset: the probabilistic routes' batch
pipeline (the rows filled, each particle's pose or class drawn from its own
row), the pose or class updates, sigma updates during search, and writing
partition-local reconstruction or class-average inputs.

The classaverager modules own the 2D class-average sums and restoration;
the 2D strategy restores in shared memory and `cavgassemble_exec`
(`simple_cavgassemble_strategy.f90`) assembles distributed parts.
`volassemble_exec` (`simple_volassemble_strategy.f90`) owns 3D volume
assembly and trailing reconstruction; `commander_volassemble` is its
command-line wrapper. These layers consume sampled-update state; they do not
own particle selection.

## 3. Bookkeeping Contracts

`sampled` marks the current sampling round. All particles with the latest
`sampled` value belong to the current active subset.

`updatecnt` tracks cumulative update history. Count-biased and fill-in paths use
it to prefer under-updated or never-updated active particles.

`sample4update_reprod` is the only correct way to reuse a previously selected
probabilistic subset. A part's matcher must not silently resample when the
master has sampled the subset.

`get_state_update_fracs` returns state-local realized update fractions for 3D
trailing from the current `sampled` round, state labels, and particles with
`updatecnt > 0`. `get_update_frac` remains available for callers that need the
legacy global summary.

`get_group_update_counts` returns the counts behind the realized fraction
`f = n/N` of each group, shared by 2D (label `class`) and 3D (label `state`):
`N(g)` active rows of the group with `updatecnt > 0`, and `n(g)` those carrying
the latest `sampled` marker. `get_state_update_fracs` is `n/N` per state. The
assembly owners (2D class-average owner, `volassemble`, the distributed PCG
(preconditioned conjugate gradients) master) compute them from the merged
project; workers never do.

A 3D reconstruction with `m_estimator=flex` reads the project's state weight
set: a weight for every particle and state, published by `flex_pca`. A particle
can then contribute to several states, so the per-state quantities are applied
masses, not row counts. `get_update_rows` returns the rows behind
`get_group_update_counts`: the active updated rows, those of the current sample,
and the rows `sample4rec` reconstructs. `N(s)` and `n(s)` are the sums of the
state's weights over the first two. Only weights above the backend's membership
threshold count: zero for gridding, `PCG_WEIGHT_THRESHOLD`
(`simple_reconstructor_pcg`, 0.01) for PCG. The realized fractions `f = n/N`,
the population a seed represents, the counts the PCG master checks its parts
against and records for each chain half, and the mass a trailing chain
represents are all such masses. Under hard labels every mass is the former
integer count.

`sample4rec` decides "nothing updated yet" over the whole project: when any
active row has `updatecnt > 0`, every range reconstructs only its updated rows.

The nominal `update_frac` is a target used by sampling. The realized fraction in
`simple_oris` is the downstream restoration and trailing contract.

### In-plane representation during probabilistic search

`inpl_cont` does not change either sampling layer or the assignment map's
schema. With `inpl_cont=yes`, callback-style local angle/shift profiling is
replaced by joint `(sx,sy,rotind_frac)` optimization, but inner importance
sampling still carries a canonical rounded in-plane index. Its stored shift is
expressed in that rounded-index frame, and no fractional angle is persisted in
the probability table.

When a sampled class or state/projection candidate is chosen for joint
profiling, its sampled `inpl` is not the continuous seed. The selected
class/state/projection remains fixed, the shift is converted to the native
particle frame, and one all-angle discrete evaluation at that exact shift
selects the profiling seed. The continuous route does not search a 5-by-5 grid
of alternative shifts.

The hard-assignment matcher owns the durable continuous result. In 3D, the
rounded probability-table assignment is authoritative: the matcher recovers
its native shift and reruns the joint optimizer locally within plus or minus two
cells of that in-plane index, without another global all-angle selection. It
persists fractional `e3`, integer `inpl`, shift, and score for the same final
pose. Valid non-improving work retains the incoming discrete pose with a
consistent re-scored objective; invalid work leaves the assignment untouched.
The 2D durable path retains its global all-angle seed selection. Neither path
may enter the legacy callback route. These policies do not alter `sampled`,
`updatecnt`, top-K support, assignment probabilities, or fractional-update
weighting.

## 4. Solve2D and Refine2D

`solve2D` uses a fixed run-local target sample size:

- default: `NSAMPLE_DEFAULT_2D = 200000`
- override: `nsample=<integer>`

The stage controller converts that target into:

```text
update_frac_2D = min(1.0, real(min(nptcls_eff, nsample_target_2D)) / real(nptcls_eff))
```

where `nptcls_eff` is the number of active particles with `state > 0`. If the
target covers almost all active particles, the stage command omits
`update_frac` and naturally becomes a full update.

Current stage policy:

- stage 1 uses the sampled-update machinery when needed, but fractional
  class-average carry-over is disabled
- while `startit == 1` and without `assign_all`, `sample_ptcls4update2D` keeps the
  initial subset sticky by reproducing it after the first random draw; a fresh
  start with `assign_all=yes` samples by update count
- later non-probabilistic iterations use `sample4update_cnt`, which is
  stochastic but biased toward particles with lower `updatecnt`
- the distributed master draws the sample of every sampled iteration,
  probabilistic or not (`master_draws_sample2D`, as 3D's
  `master_draws_sample3D`; not in the stream's pool), keyed, and the workers
  reproduce it (`refine2D_exec`), so the sample does not depend on `nparts`
- staged `assign_all=yes` (`fillin=yes` before 2026-11-03; `fillin` is 3D's
  periodic sampling of the never-updated particles) acts as a full-assignment coverage guard. It
  requires active particles to have assignments before convergence, while
  particle selection still follows the normal sampled-update path
- staged `solve2D` refinement uses sampled SNHC (`refine=snhc_smpl`) for
  stages 1-2. From stage 3 onward, `refine=prob` uses dense probabilistic
  assignment; `refine=prob_neigh` (`prob_neigh_mode=snhc`) uses sparse probabilistic SNHC until the
  final staged invocation, which uses dense `refine=prob`
- when staged updates were sampled, `solve2D` then runs a separate terminal
  dense greedy all-particle pass with `update_frac`, `fillin` and
  `assign_all` disabled, refreshing class, in-plane, and shift parameters before final
  class-average generation

Fractional 2D restoration is class-local and happens once, at the assembly
owner. Workers and the shared-memory matcher accumulate the current sample from
zero; distributed workers write `cavg_contrib_part<N>.bin` with their sums and
the class-centering offsets they applied. The owner sums the contributions in
ascending part order, reads the one partless carried set `cavg_state.bin`,
shifts it once by the centering offsets, and blends each class with the
population rule (Section 7), independently for the even/odd numerators and
CTF-squared sums. This is the 2D analogue of respecting independently updated
objects in 3D.

Distributed cleanup removes the parts' readout records and the class-sum
contribution files before each iteration; the owner deletes the contributions
after the blend. The carried set is never partition-shaped and is owned by the assembly
step only. When an iteration would blend but the carried set is missing,
unreadable or disagrees with the run (class count, `box_crop`, `smpd_crop`),
the master runs that iteration as a full update.

## 5. Solve3D and Refine3D

The 3D controller derives the solve3D outer update policy from `nsample`.
The resulting update fraction is capped by `UPDATE_FRAC_MAX`.

Current high-level solve3D stage policy (single state; the schedule is in
`3D/solve3D_policy.md`):

- stages 1 and 2 use `prob_neigh` with `prob_neigh_mode=shc`
- stages 3-5 use `prob`
- final neighborhood stages (6-8) use `prob_neigh`
- stage 1 uses `nspace=500`; stages 2-4 use `nspace=1000`
- every stage gets its low-pass and crop information independently from the
  normal schedule
- final active stages may switch to `fillin`, except where the multi-state
  policy disables it

For multi-state `solve3D` (`nstates > 1`, which refines independent states;
there is no docked mode), the default policy is an inspection-first run:
`nstages=5` and `lpstop=6.0 A` unless the user overrides them. This stops
after the `prob` phase (stages 3-5) and before `prob_neigh` (stage 6),
multi-state trailing reconstruction and staged automasking; its last stage
is NU-filtered (since 2026-10-18). The workflow still runs the final reconstruction step
at the configured last stage so it writes inspectable final state volumes. To
increase the chance that all active particles receive assignments before that
exit, the multi-state run starts stochastic balanced sampling at stage 4: the child
`refine3D` stages use `greedy_sampling=no` with `frac_best=1.0` from stage 4
onward. This keeps class-balanced quotas but draws from the whole class, not a
top-ranked fraction. The outer particle target remains the fixed
`nsample`-derived update fraction at every stage.

`sample_ptcls4update3D` applies the normal 3D subset policy, chosen by
`balance=none|class|cavg`:

- if fractional update is off, select all active particles
- `none`: no sampling units; update-count-biased sampling over the whole range
  (`sample4update_cnt`, lowest `updatecnt` tiers first)
- `class` and `cavg`: the nested square-root quota of `sample4update_class` over the
  sampling units of the class sampling file `clssmp.bin` (below)
- with `cohort_sampling=yes` (an internal flag `refine3D_states` sets on its
  `prob_neigh` frequency blocks) the stage draws by the rule above at its first
  iteration (`which_iter == startit`) and calls `sample4update_rescore` at the
  others (see "Cohort schedule")

`class` is the default wherever balanced sampling is the workflow default
(`solve3D`, `refine3D_states`, `classify3D_refs`);
`refine3D`, `refine3D_auto` and external-reference pose initialization use
`none`. `refine3D_states` and `classify3D_refs` fall back to `none` when the
project carries no selected class averages. `solve3D` with input volumes
(`vol1`) samples over the random classes it assigns and rejects `cavg`, since
those classes have no averages.

### Sampling units and the nested quota

A sampling unit is one selected 2D class (`cls2D` state > 0): its active
particles ordered by their 2D score (`corr`, best first), so the greedy and
`frac_best` selections take the best of every class. Each unit carries an
integer `group`:

- `class`: `group = 0`, every unit is its own group
- `cavg`: the selected class averages are aligned pairwise under their own
  parameters (`objfun=cc`, no CTF, `hp=30`, `lp=6`, `trs=10`) and their
  in-plane, shift and mirror invariant correlation is clustered into `nclust`
  groups (default 20) by average linkage; with no more selected classes than
  `nclust` every class is its own group. Average linkage merges the most
  similar class averages first, so a tight preferred view becomes one group
  however many classes 2D classification split it into. Cut at a fixed count
  it also leaves the outliers, odd or junk classes far from everything, as
  groups of their own, so a group holding fewer than 2.5% of the particles in
  the units (`VIEW_PARTITION_MIN_FRAC`) keeps merging into its nearest group by
  the same linkage (`hac` with weights and a minimum weight), and fewer than
  `nclust` groups may result; with no more classes than `nclust` the classes
  below that size are merged the same way. The 30 A high-pass keeps the
  overall outline out of the comparison: for a small protein in a micelle or
  nanodisc the outline dominates the low frequencies and grouped side views
  against top views rather than by the protein. `cls_expansion` groups its
  parent classes with the same alignment and linkage but without the minimum
  group size. The groups are
  written as `view_partitionNN_cavgs` stacks with the class-to-group table
  `view_partition.txt` for inspection; the project's `cluster` labels are not
  used or changed. On local execution the clustering in `solve3D` takes the
  idle workers' cores (`nparts*nthr`, capped at the cores the process owns)

The quota of an iteration, `nint(update_frac * active)`, is nested:

1. over groups in proportion to the square root of their populations, each
   capped at its population, the excess shared again over the others;
2. over the units of a group the same way (`class_sample_quotas` gives the
   real-valued quotas); a unit draws the floor of its quota, and the
   particles left to reach the target go one each to the units with a
   fractional quota whose particles have the lowest mean `updatecnt` (ties:
   larger population, then lower class index), which the project state fixes,
   so the master's draw and every partition's reproduction agree on the
   quotas and the remainder rotates over iterations;
3. inside a unit, lowest `updatecnt` first, drawn uniformly within a tier.

Under `class` the first level is the whole rule. The square root sits between
proportional sampling (`none`, which leaves the view skew of the data) and
equal shares, which correct it fully but hand a small class, often an odd or
junk one, many times its share: a unit's visits per particle relative to
another's are the inverse square root of their population ratio, so a class
100 times smaller is visited 10 times as often, not 100 times. The nested level
exists because several conformational states can hide inside one view: 2D
classification separates them into classes that average linkage merges into
one group first. The group level keeps the view balance, the unit level keeps
the hidden classes on an equal footing. The quota applies to single- and
multi-state runs alike.

Nothing in the sampler uses 3D maps, poses, projection directions or any other
quantity derived from them to decide which particles are sampled: the maps
carry bias. The `proj` field and `proj2class` serve search, convergence and 2D
averaging only.

`clssmp.bin` is the single sampling file. Every record is one unit: class
index, population, quota, `group`, then the particle indices and their scores.
The producer (`make_class_samples` in `simple_view_partition_sampling`) runs
once before the first stage; the stages only read the file. Stage command
lines carry `balance` and, for `cavg`, `nclust`. `solve3D` writes the file for
every fractional run: the run's units under `class` or `cavg`, class units
under `none`, because its initial greedy sample draws from class units.
Continuation and add-on rebuild the file from
the current `ptcl2D` class labels; no particle indices are persisted in the run
manifest, which records `balance` and `nclust` as given on the entry command
line and the effective `nsample`. `refine3D_states` and `classify3D_refs` leave
particles that are inactive in `ptcl3D` out of their units.

### Coverage planning and report

With unequal quotas the iterations that visit every particle once are
`sweep = max over units of ceil(pop / quota)`, not `ceil(1 / update_frac)`:
a unit larger than its share takes longer. `class_sample_sweep` computes it
from the unit table with the expected (real-valued) quotas of
`class_sample_quotas`. The automatic `maxits` of `refine3D_states` is the target
updates per particle times `sweep`, and that of `refine3D_auto` (balance
`none`) the target times `ceil(active / nsample)`. No `solve3D` stage rule
assumes a sweep length.

Before the first stage the producer's caller prints the unit table
(`>>> FRACTIONAL-UPDATE SAMPLING UNITS`): per group and unit the population,
the quota and the expected visits per particle over the planned draws
(iterations; frequency blocks under the cohort schedule), the draws one sweep
needs, the minimum and maximum visits over units, the view imbalance the
balancing corrected (the ratio of the two, against the visits an unbalanced
draw would give every particle) and the number of units smaller than their
quota, which are drawn whole every draw. A large ratio is a property of the
data, not a fault. The one warning is coverage short: the minimum below one,
some particles not reached in the planned draws. Nothing is adjusted
automatically: neither `nsample` nor the frequency march changes, and the
terminal missing-update pass labels any particle the march did not reach.

### Cohort schedule (`refine3D_states`)

State labelling and pose refinement converge as fixed-point iterations, so in
`refine3D_states` each `prob_neigh` frequency block refines one particle set.
`refine3D_states` sets `cohort_sampling=yes` on the command lines of those
blocks. `sample_ptcls4update3D` then draws with the grouped allocator at the
block's first iteration, which stamps a new round marker on the cohort; lowest
`updatecnt` first deprioritizes the previous cohort. At the other iterations
`sample4update_rescore` returns exactly the rows whose `sampled` equals the
latest marker, which is what `sample4update_reprod` returns, and advances
`sampled` and `updatecnt` on those rows once, so the next iteration finds the
same set under the new marker. Earlier cohorts carry older markers; `fromto`
and `allow_empty` behave as in `sample4update_reprod`. The parts' matchers
keep reproducing the round with `sample4update_reprod`. The first-iteration
clean-up of `sampled` and `updatecnt` in the master's draw (`startit == 1`,
`sample_ptcls4prob3D`) is unchanged. The terminal missing-update pass does not use
cohorts. A full visit takes `sweep`
blocks; the report states the blocks needed against those planned.
`solve3D` is not part of the cohort schedule.

Trailing reconstruction is unchanged by the cohort: the current partials are
scaled by `u/f` and the chain by `(1−u)·N/M`, with current-map coefficient
`u` (`u = f = n/N` per state unless `ufrac_trec`). Holding one cohort for `k`
iterations therefore gives it a cumulative coefficient `1 − (1−u)^k` while
earlier cohorts fade geometrically (covered by the cohort case of the
`trailing_reconstruction_blend` sub-suite). The quota biases the
composition of the selected partial reconstruction towards about
`min(pop, quota)` particles per unit per iteration; the state-level blend
computes no per-unit weight, so the composition moves with `nsample`, `nclust`
and the class selection.

`sample_ptcls4fillin` is a separate late-stage coverage policy. Its purpose is
to update particles with insufficient history, not to preserve the normal
balanced or count-biased exploration distribution.

The 3D matcher writes partial reconstructions from the active subset. Volume
assembly then restores volumes, calculates FSCs, postprocesses references, and
applies trailing reconstruction when requested. Trailing uses an explicit
`ufrac_trec` only when the parsed `params%l_ufrac_trec_defined` flag is true
and the run is single-state; otherwise it consumes realized per-state fractions
from `get_state_update_fracs`, or under a state weight set the mass fractions of
Section 3. The numeric `params%ufrac_trec` field has a
default value and must not be interpreted as an active override by itself.
Multi-state convergence reporting records those effective state-local fractions
as `TRAIL_REC_UPDATE_FRAC_STATE01`, `TRAIL_REC_UPDATE_FRAC_STATE02`, and so on.

## 6. Probabilistic Alignment

Probabilistic alignment is not a second outer sampler.

The workflow is (`doc/implementation_notes/completed/prob_align_batch_pipeline.md`,
2026-11-04):

1. the master (or the shared-memory matcher) chooses the outer subset through
   the normal 2D or 3D sampling helper, with a draw keyed on the run's seed,
   the pass, the iteration and the particle
2. the master writes the sampled project state and schedules its parts once
3. each part's matcher reproduces its share of the subset
   (`sample4update_reprod`), and for each batch fills the probability-table
   rows, draws each particle's pose or class from its own row (keyed on the
   particle), polishes, computes the residual and commits: the hard particle
   update in the same pass
4. a part publishes a readout record of its counts and sums; the master merges
   the records in part order and prints the readouts

`simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90`, and
`simple_eul_prob_tab2D.f90` perform candidate-level importance sampling inside
that subset. They may use score-derived candidate distributions,
`angle_sampling`, `greedy_sampling`, or neighborhood sampling, but the selected
particle set is already fixed before they run. Every draw of theirs is keyed
(the in-plane draw per scored reference, the shift map's disk draw, the sparse
fills' candidate order, the assignment draw), so it depends on the seed, the
pass, the iteration, the particle and the draw, never on a thread, a batch or
a part.

## 7. Restoration and Assembly

Both the 2D and the 3D blends follow the population rule below. Each stored set records `M(g)`, the
population its sums represent. With `N(g)` and `n(g)` from the merged project, `f = n/N`, and `u = f`
or an applied override:

- current contribution scaled by `s = u / f` (1 by default)
- previous contribution scaled by `w = (1 - u) * N / M` (`(N - n) / M` by default)
- the new set records `M <- s*n + w*M`, which is `N` whenever `M > 0`

The carried mass after the blend is therefore the represented population `N`
whatever joined or left the group: first-time particles add mass without
displacing old mass, deactivated rows take their share of old mass with them,
and returning rows get theirs back (`w > 1` is allowed). With no population
change `w = 1 - u`, the former recurrence. The rule keeps the mass right, not
the membership: the aggregate is the state of a stochastic recurrence, not an
exact sum over the current class or state members. Classes or states with no
active updated particles carry nothing; full sampled participation replaces the
previous sums.

`population_blend_weights` (`simple_oris`) computes these weights. It is
generic over integer counts and real masses; integer counts go through the mass
version with identical arithmetic. A state weight set (Section 3) therefore
changes the inputs `N`, `n` and `M`, not the rule.

2D class-average restoration applies the rule per class at the assembly owner
(Section 4); the restored averages use the represented even/odd populations.

3D volume assembly performs trailing in the accumulator domain, mirroring the
2D scheme. The persistent per-state chain (`trailrec_stateNN_{even,odd}` plus
rho files and a `trailrec_stateNN.txt` manifest) holds blended, unregularized
e/o Fourier sums and sampling densities at the mass of the population it
represents, `M(s)`, which the manifest records. Two fractions govern the blend:

- `f` — the realized state-local fraction that produced the current partials
  (`get_state_update_fracs`, or the mass fraction under a state weight set);
  always computed
- `u` — the applied map-update weight; equals `f` unless a single-state
  `ufrac_trec` override is provided

The recurrence keeps the chain at the mass of the represented population `N`
and makes `u` the restored current-map coefficient, preserving the historical
`ufrac_trec` meaning:

- current contribution: partial sums and rho scaled by `s = u / f`
  (mass `u * N`)
- previous chain contribution: sums and rho scaled by `w = (1 - u) * N / M`
  (mass `(1 - u) * N`); with an unchanged population `w = 1 - u`
- the chain records `M <- N`
- a single sampling-density correction after the blend restores the trailed
  halves, so each Fourier component is weighted by its accumulated sampling
  density; the FSC is estimated post-blend and describes the on-disk artifact

The chain is written before restoration (regularization mutates rho in place).
Every blend happens in the accumulator domain; restored (finished) half maps
are never blended. When a blend is due but the chain does not exist yet, the
iteration does not blend: volassemble seeds the chain with the current partials
scaled by `1/f`, so the stored chain carries full-dataset mass and the next
iteration's effective update weight is the requested fraction (an unnormalized
fractional seed would make a 10 percent request act like a ~53 percent update),
and restores and ships the current sample alone (partials scaled back by `f`);
the FSC describes those halves. The previous volumes are neither read nor
needed and no extra reconstruction runs; the chain fills over about `1/f`
iterations. A state whose realized fraction is below 0.001 carries a valid
chain unchanged and is restored from it, or, without a chain, keeps its
previous volume like a dropped state (when the directory holds that volume, its
half maps and FSC; with neither a chain nor such a volume the run stops with an
error, because a map of the sample alone would be degenerate).
Stage-boundary full reconstructions
seed the chain at full-dataset weight through the internal `trail_seed`
handshake, but only when the consuming stage actually trails. Every seed records
the population it represents: `N` for the chain-start seed (`1/f` scaling), and the
rows `sample4rec` reconstructs for a stage-boundary seed (their applied mass
under a state weight set). The distributed PCG chain applies the same weights.
Its two halves share a manifest of the same kind, `pcg_trail_stateNN.txt`, and
the represented population is read from it. The raw header of each half also
carries its contributor count, applied mass and weight-set identity (chain
identity `pcgtrail-v3`).

The four accumulator files plus manifest form one artifact set. The manifest is
deleted before and rewritten after the data files with per-component byte
sizes, generation counter, provenance (box, sampling, particle population,
state layout), a format version, the represented population `M(s)`, the
contributor count and the state weight set identity (its generation and layout
digest; zero under hard labels), so interrupted writes never validate. The PCG
manifest follows the same order: it is deleted before the first half is
written, written after both, and checked (identity, sizes, state layout) before
use. A chain whose manifest has an older format version is discarded and
re-seeded. A chain built under another weight-set identity is re-seeded, never
blended; this holds for both backends. Readers accept a chain
only when the manifest parses, provenance matches the current project, every
component size matches, and the grid is not larger than the current one with
the same physical extent; smaller grids are zero-padded on read (downsampling
ramp). Any validation failure discards the complete set and re-seeds.
Cross-directory continuation carries chains over as complete sets only,
manifest last.

Neither class-average restoration nor volume assembly should make new particle
sampling decisions. If a restoration or assembly change requires a different
subset policy, that policy belongs in the commander/controller/sampling-helper
layer and must be reflected in `sampled` and `updatecnt`.

Online matcher restoration/reconstruction paths must read active particle
images from disk once per batch and reuse those batch images for both matching
and restoration/reconstruction. Do not introduce a trailing full reconstruction
or class-average restoration pass that re-reads image stacks as a memory
optimization unless the single-read performance contract is explicitly changed.
Probabilistic table-generation programs and explicit offline assembly commands
are separate workflow stages and may perform their own reads.

## 8. Invariants

- Outer particle sampling happens on the master before the parts fill their rows.
- Every part's matcher reproduces the master's subset.
- Candidate importance sampling never changes the particle subset.
- `sampled` remains the current-round marker.
- `updatecnt` remains cumulative update history.
- Downstream restoration uses realized update state, not only nominal
  `update_frac`.
- Stage 1 of `solve2D` may be sampled but must not fractionally carry over
  previous class-average sums.
- Independent multi-state `solve3D` defaults to a five-stage,
  `lpstop=6.0 A` inspection run, starts stochastic balanced sampling at stage
  4, and still writes final reconstruction outputs.
- 2D fractional class-average restoration remains class-local.
- Staged `solve2D` `fillin=yes` remains a full-assignment coverage guard
  unless the implementation is deliberately changed to missing-only assignment.
- Sampled `solve2D` runs a terminal dense greedy all-particle refresh before
  final class-average generation.
- `volassemble` and the classaverager remain consumers of sampled-update state,
  not producers of particle-selection policy.
- Before an assembly that rereads the project file (gridding `volassemble`),
  the iteration strategy writes the particle field the matcher left: pose,
  state, `sampled`, `updatecnt` and the partial reconstructions describe one
  iteration (`publish_before_file_assembly` in `simple_refine3D_strategy`).
- Online matcher restoration/reconstruction reuses the particle images already
  read for the current batch.

## 9. Review Checklist

For sampling, probabilistic alignment, class-average restoration, or volume
assembly changes, check:

- Does the outer subset get selected exactly once, by the master, for a
  probabilistic iteration?
- Do the parts' matchers reuse the recorded subset through
  `sample4update_reprod`?
- Is candidate-level importance sampling kept separate from particle-level
  subset selection?
- With `inpl_cont=yes`, does candidate profiling retain only rounded in-plane
  metadata and defer durable fractional `e3` to the final hard assignment?
- Does candidate profiling perform one all-angle selection at the supplied
  shift without a coarse shift scan?
- Does final 3D refinement retain the authoritative rounded assignment and
  avoid a second global angle selection?
- Can any joint no-improvement or invalid-result path accidentally enter the
  legacy callback?
- Are `sampled` and `updatecnt` updated consistently before downstream
  restoration or trailing consumes them?
- Does 2D restoration use class-local realized fractions?
- Does 3D trailing consume the realized or explicit trailing fraction?
- Does the online matcher path preserve one image-stack read per particle batch?
- Are shared-memory and distributed paths preserving the same scientific
  workflow and artifact contracts?
