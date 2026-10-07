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

Probabilistic pre-alignment is a sample-once-and-reuse path: the pre-alignment
commander chooses the outer subset, probability-table workers reuse it, and the
matcher reuses it again for the hard particle update.

## 2. Ownership

`simple_commanders_solve2D.f90` owns `solve2D` orchestration: defaults,
stage execution, final fill-in, and final class-average generation.

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

`simple_commanders_prob.f90` owns probabilistic pre-alignment orchestration:
sampling the outer subset once, writing it to the project, running table
generation, aggregating probability-table outputs, and writing the assignment
artifact.

`simple_eul_prob_tab*.f90` owns inner candidate importance sampling. These
modules may sample references, orientations, neighbors, or in-plane candidates
inside the active particle subset, but they must not choose a new particle
subset.

`simple_strategy2D_matcher.f90` and `simple_strategy3D_matcher.f90` own
particle-domain search on the active subset, assignment consumption, pose or
class updates, sigma updates during search, and writing partition-local
reconstruction or class-average inputs.

The classaverager modules own 2D class-average restoration and assembly.
`commander_volassemble` owns 3D volume assembly and trailing reconstruction.
These layers consume sampled-update state; they do not own particle selection.

## 3. Bookkeeping Contracts

`sampled` marks the current sampling round. All particles with the latest
`sampled` value belong to the current active subset.

`updatecnt` tracks cumulative update history. Count-biased and fill-in paths use
it to prefer under-updated or never-updated active particles.

`sample4update_reprod` is the only correct way to reuse a previously selected
probabilistic subset. A probability-table worker or downstream matcher must not
silently resample when a probabilistic pre-step has already sampled the subset.

`get_state_update_fracs` returns state-local realized update fractions for 3D
trailing from the current `sampled` round, state labels, and particles with
`updatecnt > 0`. `get_update_frac` remains available for callers that need the
legacy global summary.

`get_group_update_counts` returns the counts behind the realized fraction
`f = n/N` of each group, shared by 2D (label `class`) and 3D (label `state`):
`N(g)` active rows of the group with `updatecnt > 0`, and `n(g)` those carrying
the latest `sampled` marker. `get_state_update_fracs` is `n/N` per state. The
assembly owners (2D class-average owner, `volassemble`, the distributed PCG
master) compute them from the merged project; workers never do.

`sample4rec` decides "nothing updated yet" over the whole project: when any
active row has `updatecnt > 0`, every range reconstructs only its updated rows.

The nominal `update_frac` is a target used by sampling. The realized fraction in
`simple_oris` is the downstream restoration and trailing contract.

### In-plane representation during probabilistic search

`inpl_cont` does not change either sampling layer or the assignment-file
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
- while `startit == 1`, `sample_ptcls4update2D` keeps the initial subset sticky
  by reproducing it after the first random draw
- later non-probabilistic iterations use `sample4update_cnt`, which is
  stochastic but biased toward particles with lower `updatecnt`
- probabilistic stages use `prob_align2D` to sample once, then `prob_tab2D` and
  `refine2D_exec` reproduce the same subset
- staged `fillin=yes` currently acts as a full-assignment coverage guard. It
  requires active particles to have assignments before convergence, while
  particle selection still follows the normal sampled-update path
- staged `solve2D` refinement uses sampled SNHC (`refine=snhc_smpl`) for
  stages 1-2. From stage 3 onward, `refine=prob` uses dense probabilistic
  assignment; `refine=prob_snhc` uses sparse probabilistic SNHC until the
  final staged invocation, which uses dense `refine=prob`
- when staged updates were sampled, `solve2D` then runs a separate terminal
  dense greedy all-particle pass with `update_frac` and `fillin`
  disabled, refreshing class, in-plane, and shift parameters before final
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

Distributed cleanup removes assignment, distance and class-sum contribution
files before each iteration; the owner deletes the contributions after the
blend. The carried set is never partition-shaped and is owned by the assembly
step only. When an iteration would blend but the carried set is missing,
unreadable or disagrees with the run (class count, `box_crop`, `smpd_crop`),
the master runs that iteration as a full update.

## 5. Solve3D and Refine3D

The 3D controller derives the solve3D outer update policy from `nsample`.
The resulting update fraction is capped by `UPDATE_FRAC_MAX`.

Current high-level solve3D stage policy:

- stages 1 and 2 use `prob_neigh` with `prob_neigh_mode=shc`
- stages 3-5 use `prob`
- final neighborhood stages use `prob_neigh`
- stage 1 uses `nspace=500`; stages 2-4 use `nspace=1000`
- every stage gets its low-pass and crop information independently from the
  normal schedule
- final active stages may switch to `fillin`, except where the multi-state
  policy disables it

For multi-state `solve3D` (`nstates > 1`, which refines independent states;
there is no docked mode), the default policy is an inspection-first run:
`nstages=5` and `lpstop=6.0 A` unless the user overrides them. This stops
after the `prob` phase and before
`prob_neigh`, staged NU filtering, multi-state trailing reconstruction,
and staged automasking. The workflow still runs the final reconstruction step
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
- `class` and `cavg`: the nested equal quota of `sample4update_class` over the
  sampling units of the class sampling file `clssmp.bin` (below)
- with `cohort_sampling=yes` (an internal flag `refine3D_states` sets on its
  `prob_neigh` frequency blocks) the stage draws by the rule above at its first
  iteration (`which_iter == startit`) and calls `sample4update_rescore` at the
  others (see "Cohort schedule")

`cavg` is the default wherever balanced sampling is the workflow default
(`solve3D`, `refine3D_states`, `classify3D_refs`);
`refine3D`, `refine3D_auto` and external-reference pose initialization use
`none`. `refine3D_states` and `classify3D_refs` fall back to `none` when the
project carries no selected class averages. `solve3D` with input volumes
(`vol1`) defaults to `class` over the random classes it assigns and rejects
`cavg`, since those classes have no averages.

### Sampling units and the nested quota

A sampling unit is one selected 2D class (`cls2D` state > 0): its active
particles ordered by their 2D score (`corr`, best first), so the greedy and
`frac_best` selections take the best of every class. Each unit carries an
integer `group`:

- `class`: `group = 0`, every unit is its own group
- `cavg`: the selected class averages are aligned pairwise under their own
  parameters (`objfun=cc`, no CTF, `lp=6`, `trs=10`, as `cluster_cavgs`) and
  their in-plane, shift and mirror invariant correlation is clustered into
  `nclust` groups (default 20) by average linkage; with no more selected
  classes than `nclust` every class is its own group. Average linkage merges
  the most similar class averages first, so a tight preferred view becomes one
  group however many classes 2D classification split it into. The groups are
  written as `view_partitionNN_cavgs` stacks with the class-to-group table
  `view_partition.txt` for inspection; the project's `cluster` labels are not
  used or changed. On local execution the clustering in `solve3D` takes the
  idle workers' cores (`nparts*nthr`, capped at the cores the process owns)

The quota of an iteration, `nint(update_frac * active)`, is nested:

1. equal over groups: every group gets one more particle per round until the
   target is reached, capped at its population (the total may exceed the
   target by fewer particles than there are groups);
2. equal over the units of a group, capped at their populations; the remainder
   of the equal split goes to the units whose particles have the lowest mean
   `updatecnt` (ties: larger population, then lower class index), which the
   project state fixes, so every partition computes the same quotas and the
   remainder rotates over iterations;
3. inside a unit, lowest `updatecnt` first, drawn uniformly within a tier.

Under `class` the first level is the whole rule, as before. The nested level
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

With equal quotas the iterations that visit every particle once are
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
and `allow_empty` behave as in `sample4update_reprod`. `prob_tab` and the
matcher keep reproducing the round with `sample4update_reprod`. The
first-iteration clean-up of `sampled` and `updatecnt` in `prob_align`
(`startit == 1`) is unchanged. The terminal missing-update pass does not use
cohorts. A full visit takes `sweep`
blocks; the report states the blocks needed against those planned.
`solve3D` is not part of the cohort schedule.

Trailing reconstruction is unchanged by the cohort: the current partials are
scaled by `u/f` and the chain by `(1−u)·N/M`, with current-map coefficient
`u` (`u = f = n/N` per state unless `ufrac_trec`). Holding one cohort for `k`
iterations therefore gives it a cumulative coefficient `1 − (1−u)^k` while
earlier cohorts fade geometrically (covered by the cohort case of the
`trailing_reconstruction_blend` sub-suite). The equal quota biases the
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
from `get_state_update_fracs`. The numeric `params%ufrac_trec` field has a
default value and must not be interpreted as an active override by itself.
Multi-state convergence reporting records those effective state-local fractions
as `TRAIL_REC_UPDATE_FRAC_STATE01`, `TRAIL_REC_UPDATE_FRAC_STATE02`, and so on.

## 6. Probabilistic Pre-Alignment

Probabilistic pre-alignment is not a second outer sampler.

The workflow is:

1. choose the outer subset through the normal 2D or 3D sampling helper
2. write the sampled project state
3. run probability-table generation only for that subset
4. aggregate table outputs into one assignment artifact
5. reproduce the same subset in the matcher
6. perform the hard particle update

`simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90`, and
`simple_eul_prob_tab2D.f90` perform candidate-level importance sampling inside
that subset. They may use score-derived candidate distributions,
`angle_sampling`, `greedy_sampling`, or neighborhood sampling, but the selected
particle set is already fixed before they run.

## 7. Restoration and Assembly

Both the 2D and the 3D blends follow the population rule (class-average state
note, Section 4.1). Each stored set records `M(g)`, the population its sums
represent. With `N(g)` and `n(g)` from the merged project, `f = n/N`, and `u = f`
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

2D class-average restoration applies the rule per class at the assembly owner
(Section 4); the restored averages use the represented even/odd populations.

3D volume assembly performs trailing in the accumulator domain, mirroring the
2D scheme. The persistent per-state chain (`trailrec_stateNN_{even,odd}` plus
rho files and a `trailrec_stateNN.txt` manifest) holds blended, unregularized
e/o Fourier sums and sampling densities at the mass of the population it
represents, `M(s)`, which the manifest records. Two fractions govern the blend:

- `f` — the realized state-local fraction that produced the current partials
  (`get_state_update_fracs`); always computed
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
rows `sample4rec` reconstructs for a stage-boundary seed. The distributed PCG
chain applies the same weights; its represented population is the particle
count in the raw header of each half (chain identity `pcgtrail-v3`).

The four accumulator files plus manifest form one artifact set. The manifest is
deleted before and rewritten after the data files with per-component byte
sizes, generation counter, provenance (box, sampling, particle population,
state layout), a format version and the represented population `M(s)`, so
interrupted writes never validate. A chain written by an older build, whose
manifest has no represented population, is discarded and re-seeded. Readers accept a chain
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

- Outer particle sampling happens before probabilistic table generation.
- Probabilistic table workers and downstream matchers reproduce the same subset.
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
- Online matcher restoration/reconstruction reuses the particle images already
  read for the current batch.

## 9. Review Checklist

For sampling, probabilistic alignment, class-average restoration, or volume
assembly changes, check:

- Does the outer subset get selected exactly once for a probabilistic
  pre-alignment iteration?
- Do table workers and matchers reuse the recorded subset through
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
