# Solve3D Policy

This document records the current policy for `solve3D`, the staged
particle-based de novo map determination workflow: it couples ab initio 3D
reconstruction from random orientations with initial 3D refinement of the
resulting map in one stage schedule. `solve2D` is its 2D equivalent. The
programs were called `abinitio3D` and `abinitio2D` before 2026-10-03 (release
4 keeps no aliases for the old names). The base `refine3D` contracts are in
[refine3D_policy.md](refine3D_policy.md); this document describes how
`solve3D` configures and chains those stages.

Multi-state matcher and reconstruction lifetimes follow
[separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md](separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md).

## 1. Scope

`solve3D` determines 3D maps de novo from particles by preparing starting
orientations/states, marching through staged `refine3D` runs (global-search
ab initio reconstruction stages, then local-search initial refinement stages),
optionally performing symmetry-axis search, and reconstructing final
original-sampling maps.

It owns stage scheduling. It does not own a separate particle matcher or volume
assembly implementation.

## 2. Defaults

`solve3D` sets:

- `objfun=euclid`
- `bfac=0`

When unset, it supplies:

- `mkdir=yes`
- `overlap=0.95`
- `prob_athres=10`
- `center=no`
- `cenlp` from the solve3D controller default
- `oritype=ptcl3D`
- `pgrp=c1`
- `pgrp_start=c1`
- `filt_mode=nonuniform`
- `automsk=no`
- `gauref=yes`
- `balance=cavg` (`class` with input volumes, which bring random classes and
  no class averages, `cavg` being rejected there; `none` for a project
  without a 2D solution)

For a multi-state run (`nstates > 1`), it also supplies conservative
inspection defaults when the user has not overridden them:

- `nstages=5`
- `lpstop=6.0 A`

The public `filt_mode` values are `none`, `nonuniform`, and
`nonuniform_lpset`. Automatic low-pass modes `uniform` and `fsc` are rejected
for `solve3D`.

## 3. Stage Controller

The stage controller in `simple_solve3D_controller.f90` emits a concrete
`refine3D` command line for each stage. The full particle workflow has eight
stages. Independent multi-state startup defaults to five stages, so it stops
after its `prob` stages 3-5, before the first `prob_neigh` stage (6); stage 5 is
also its first NU-filtered stage.

Stage policy includes:

- stage-specific `nspace` and `maxits`
- cropped box and sampling from the low-pass plan
- stage 1 starts with `nspace=500`; stages 2 through 4 use `nspace=1000`
- every stage gets its low-pass and crop information independently from the
  normal FRC/input schedule
- staged search by state count:
  - single-state particle runs use `prob_neigh` with `prob_neigh_mode=shc`
    in stages 1-2, middle `prob`, and late `prob_neigh`
  - multi-state particle runs use direct `shc` in stages 1-2, `prob` in
    stages 3-5, and `prob_neigh` with `prob_neigh_mode=state` from stage 6
    (`PROB_NEIGH_REFINE_STAGE`, the single-state switch), so the default
    five-stage multi-state run uses no `prob_neigh`; the class-average route
    switches at the same stage (section 6)
- the multi-state direction cap: with `nstates > 1` no stage exceeds 2,500
  directions, so stages 7 and 8 (5,000 in single state) run at 2,500
- the multi-state stage-2 exit: stage 2 of a multi-state particle run (`shc`)
  also ends once the search scans more than 95% of the search space on
  average, whatever the state overlap (section 6)
- `nspace_sub` for `prob_neigh`
- staged point-group policy between `pgrp_start` and `pgrp`
- staged translation limits
- staged ML regularization
- staged fractional update with a fixed `nsample` target while `nsample/active_particles <= 0.9`
- fractional-update selection over the sampling units of `balance`: `cavg` by
  default (groups of similar class averages, then their classes), `class` (one
  unit per selected class) or `none` (global lowest `updatecnt` tiers)
- stochastic sampling start by state count
- early Gaussian reference filtering
- optional trailing reconstruction by stage and state count
- staged NU filtering from `NU_FILTER_STAGE`
- staged automasking only from `AUTOMSK_STAGE`
- an explicitly requested PCG solvent prior from `PCG_SOLVENT_START_STAGE`
  (stage 8, the final stage), so the NU label field the prior'd pair receives
  has settled over stages 5-7 without the prior

The downscaled particle cache is a 2D-only feature: `solve3D` rejects
`cache=yes`, and each stage uses its own crop from the low-pass/downscaling
ladder.

The effective crop and its physically equivalent pixel size are shared by
starting-volume generation, matcher reconstruction, stage-boundary
reconstruction, symmetry-map commands, FSC diagnostics, low-pass snapshots, and project volume
registration. This preserves `box * smpd == box_crop * smpd_crop` across every
particle-to-volume handoff.

### Full-Sampling Switch

`solve3D` now applies a global sampling override when:

`nsample / active_particles > 0.9`

In that regime, the staged controller forces full active-particle updates for
each iteration by suppressing fractional controls (`update_frac`, `nsample`,
`fillin`) in emitted child `refine3D` commands. This also disables trailing
reconstruction for staged `solve3D` commands, and startup class-biased
sampling setup is bypassed in favor of all-active sampling.

### Sampling units

Every stage runs with the run's `balance` (default `cavg`). Below the
full-sampling switch, `solve3D` writes the class sampling file `clssmp.bin`
once at startup: one unit per selected 2D class, and under `cavg` each unit's
group among `nclust` (default 20) groups of similar class averages, formed by
average linkage on their aligned correlation, so a preferred view spread over
many classes no longer dominates the sample. Every fractional update shares
its target equally over the groups, then equally over the classes of a group,
each capped at its population, then lowest `updatecnt` first inside a class.
Under `class` every class is its own group; under `none` the stages sample the
global lowest `updatecnt` tiers, and the file holds class units only for the
initial greedy sample. The groups are written
as `view_partitionNN_cavgs` stacks with the table `view_partition.txt`. The
units are formed once, before the first stage, and serve the whole run: the
stages only read the file, and every stage command line carries `balance` and,
for `cavg`, `nclust`. Before the first stage `solve3D` prints the unit table
with the expected visits per particle over the planned stage iterations, the
view imbalance corrected, and a warning on short coverage only; nothing is
adjusted automatically. The
sampler never uses 3D maps, poses or projection directions. The details are in
`doc/policies/importance_sampling_fractional_update_policy.md`.

The emitted child command line owns `startit` and `which_iter` for the current
stage. `refine3D` then treats `maxits` as the run length for that stage.

### Without a 2D solution

`solve3D` and `solve3D_addon` run on a project without a usable 2D solution
(`ptcl2D` never searched, or no classes in `cls2D`, or no class FRCs
registered, as in a project built from particles alone with per-particle
labels) under `balance=none` only; `class` and `cavg` are refused with a message naming `balance=none`, as
are `cavg_ini` and `cavg_ini_ext`. `balance` then defaults to `none`. The
selection is unchanged: `ptcl3D%state` is reset from the `ptcl2D` state flags
(set at import or by a selection). The 2D shifts are transferred when
`ptcl2D` was searched and otherwise the shifts `ptcl3D` holds stay. No class sampling file is written; the initial
greedy draw and every stage sample the global lowest `updatecnt` tiers
(`sample4update_cnt`), as every `balance=none` stage does. The ladder cannot
be planned from class FRCs: `lpstart`/`lpstop` are taken from the command line
and otherwise `lpstart` from `mskdiam` (`mskdiam2lplimits`) and `lpstop` is the
ladder bound `LPSTOP_BOUNDS(1)`, 4.5 A, as on the input-volume
route; an add-on takes its ladder from the frozen run's manifest regardless.
Input volumes on such a project no longer need the random classes the
`class` default fabricates.

### Stage-start sigma bootstrap (ini3D routes)

The `cavg_ini` and `cavg_ini_ext` routes enter the particle stages at a stage
whose starting reconstruction is ML-regularized (`objfun=euclid`) before any
refine3D iteration has estimated particle sigmas.
`calc_rec` applies the single sigma2 bootstrap rule (`simple_sigma2_bootstrap`,
`doc/policies/3D/refine3D_policy.md` section 5): for a euclid starting
reconstruction it calls `ensure_sigma2_for_iteration`, which is a no-op when
the project owns a compatible committed state and otherwise seeds canonical
particle spectra from image power with `calc_pspec`. The starting reconstruction
then runs as planned, and the stage's first euclid iteration replaces the seed
with residual spectra in the next canonical transaction.
When the final reconstruction at original sampling must rebuild its sigmas
(section 9) it is one program call, `bootstrap_rec3D`, which owns the complete
sequence: the same image-power seed, a euclid ML bootstrap map on it, one
residual sigma pass (`refine=sigma`) against that map, canonical commit of the
residual generation, and the shipped euclid ML reconstruction on it. Because the
program takes any project with 3D orientations, it is also the standalone
test for this stage (`simple_exec prg=bootstrap_rec3D projfile=... pgrp=...
mskdiam=... nparts=... nthr=... rec_backend=...`), so failures in the final
reconstruction can be reproduced in minutes rather than after a full run.
Its bootstrap map is a gridding assembly that carries the last stage's
`filt_mode` and `automsk` (the residual sigmas depend on the
regularization of the reference they are scored against, so that reference is
regularized like the last stage's matching references); the shipped map is
classical (`filt_mode=none`, except that `automsk=nu` keeps the caller's NU
`filt_mode`) and runs on the workflow's backend with the PCG cold-solve budget
applied inside `bootstrap_rec3D` (2026-09-07).
Whether the final reconstruction refreshes its sigmas at native sampling is
decided by the registration-box rule: a registration box different from the
native box forces refresh (2026-09-07). Solve3D drops any state
registration inherited with its input project before the first stage, so
every run seeds and owns its sigmas in its own directory.

## 4. Low-Pass and Cropping

`lpinfo(istage)%lp` controls staged search/reference scheduling. Stage limits
are derived from class FRCs by default, with `lpstart`/`lpstop` overrides and a
`force_lp_range=yes` path that uses the requested range directly.

Where the ladder is planned by `lpstages` or `lpstages_fast`, stage 1 of
`solve3D` never runs with a low-pass limit finer than 20 A: both apply this
floor before generating the remaining low-pass and crop ladder, preserving
gradual frequency marching. `force_lp_range=yes` and the routes planned by
`lpstages_setlims` (input volumes, no 2D solution) take `lpstart` as given.

**Aggressive frequency marching (`freq_march_aggressive=yes`, 2026-10-21;
experimental, default `no`).** For a small membrane target the class FRCs
measure micelle self-consistency at the low resolutions where the default
ladder parks stages 1-3 (PfCRT: 20.0, 19.4, 18.3 A), and stage 3's band does
not discriminate orientations: on the `latest6` reference its overlap peaks
at the fourth iteration and declines for the remaining thirteen (the
iteration-budget note below). The signal that fixes a small particle's
orientation sits in the 5-10 A band, which the default schedule withholds
until stage 4. With the key on, the particle route plans its ladder from the
class FRCs as usual and then `lpstages_aggressive` overwrites the rungs from
`AGGRESSIVE_MARCH_STAGE` (= `PROB_REFINE_STAGE`, 3) on: linear in spatial
frequency from `LPSTART_AGGRESSIVE` (12 A, no coarser than the kept stage-2
rung) to the standard final limit `lpfinal`, so every march stage admits the
same number of new Fourier shells and the steps are smallest in Angstrom at
the fine end (12.0, 10.0, 8.6, 7.5, 6.7, 6.0 A for a 6 A final limit); each
march stage is cropped on its own at the particle-route target (smpd lp/3, at
least 2 A), never below the preceding stage's crop. The finer bands need
finer angular sampling, so `nspace` follows `NSPACE_AGGRESSIVE` = [500, 1000,
2500, 2500, 2500, 2500, 5000, 5000] (a spiral of n directions resolves about
223/sqrt(n) degrees; a band lp on a particle of diameter D needs about lp/D
radians, 5.7 degrees at 12 A on 120 A). Everything else stands: `MAXITS`, the
stage modes, the NU start stage (5; NU before that reinforces a bad
alignment, see section 8), the non-NU matching ceiling at the rung, the
FSC=0.5 promotion, stage omission and mid-stage exit, which now act on the
closer rungs.

*Why the march starts at stage 3 and not at stage 1.* The first version
(PfCRT `latest10`, 2026-10-21) marched from stage 1: 12.0, 10.5, 9.3 A with
`nspace` 1000/2500. The search from random orientations never consolidated:
stage 1's overlap plateaued at 0.22 (0.34 at 19.4 A in `latest6` and
`latest9`), stage 2's at 0.17 (0.23), the FSC=0.5 crossing sat at the band
(12.9 A) in stage 1 and then *coarsened* to 13.5-15.5 A in stage 2 (noise
fitted inside a band the reference has no content in), the cFAR of 0.6-0.7
was that of isotropic noise, and stage 3 consolidated to 0.32 instead of
0.50: garbage after two stages where the default gives a smooth, correct
volume. The reference built from random poses carries nothing at 12 A, so the
stochastic hill climb has no landscape to climb; the low-resolution start is
what makes the global search converge. In the default the large
consolidation (overlap 0.23 to 0.50) happens in stage 3 at nearly stage 2's
band, so it is the switch to the full `prob` search with `ml_reg`, not the
band, that does it; the march therefore begins with that switch. 12 A at
stage 3 is also about res05/1.5 after stage 2 on PfCRT (res05 18.3 A), the
bounded lead the NU stages are known to climb under.

The key applies to the particle route with a 2D solution only: the
input-volume and no-2D routes and the add-on have explicit limits and refuse
it, the class-average route keeps its own ladder and refuses it, and
`force_lp_range=yes` is refused beside it. Test: the ten-restart PfCRT
regression against the default schedule, scored on stage 3's overlap
trajectory (does it still churn?), the stage-5 FSC, the iterations to the
stage-6 Nyquist, overlap through stage 6 and above all the gold-standard
bootstrap of the final poses, plus the same restarts on a target that fails
under the default schedule.

The single-state external-volume route (`vol1`, `nstates=1`) starts at
`VOL_SUPPLIED_STAGE`, the stage after the no-axis-search block
(`NOSYMSRCH_STAGE + 1`, stage 4 in the standard schedule). That first
executable stage is anchored at exactly `lpstart`; the unused preceding stages
carry the same low-pass and crop metadata but are omitted from schedule
diagnostics. The remaining stages march linearly from `lpstart` to `lpstop`.
Those delayed rungs are the planned limits used by FSC=0.5 promotion, stage
shift and mid-stage exit; the entry stage remains ineligible for shift or exit.
Multi-state external-volume schedules and other callers of
`lpstages_setlims` retain the ordinary full-ladder interpolation.

The controller passes an `lpstop` ceiling to each staged `refine3D` child
alongside the effective planned matching limit. In the non-NU stages the
ceiling is the stage's matching limit, so matching never exceeds the
printed stage limit.

Stage-boundary FSC=0.5 advance (2026-09-06; one ladder rung since
2026-10-18, particle route only): past stage 2, when the project's FSC=0.5
resolution of the best resolved populated state (the per-particle `res05`
field written by the reconstruction) is finer than the planned
`lpinfo(istage)%lp`, the stage runs at the *next* rung
`lpinfo(istage+1)%lp` instead, bounded by the ladder cap
(`LPSTOP_BOUNDS(1)`, 4.5 A, or the coarser explicit `lpstop`; the
class-average route keeps its final limit). `res05` gates the advance; the
ladder sets its size.

**Stage omission (2026-10-18).** Where the stage can be omitted, it is omitted
rather than run at an advanced limit. When `res05` is finer than
`lpinfo(istage)%lp` the stage's own crop cannot represent the resolution the
halves already agree on -- `set_bp_range3D` clips the band to the crop
Nyquist, so the stage would print an advanced limit while matching at the
crop limit and spend its whole iteration budget there. `solve3D_stage_behind_fsc05`
gates a `cycle` in the frequency-marching loop, logging `>>> SOLVE3D OMITTING
STAGE`. Excluded, because they carry one-time work: the first stage
(`start_stage`, the add-on registration pass), the last stage (automasking,
envelope FSC, trailing reconstruction) and, **only when an axis is actually
searched**, `SYMSRCH_STAGE`. `NOSYMSRCH_STAGE` can be omitted: with
`pgrp == pgrp_start` the `symmetrize` call is a complete no-op, since
`l_symran` is set only when `pgrp_start` has more symmetry than `pgrp` and
therefore implies `l_srch4symaxis`. Omitted
stages leave `emitted_lp` at its `-1` sentinel, which the stage report and
the add-on manifest already treat as "not run". The one-rung advance above
therefore applies to the stages that cannot be omitted: the entry stage when it
lies past stage 2 and, when an axis is searched, `SYMSRCH_STAGE`; the last stage
has no next rung and is never advanced. Single-state
runs only since 2026-10-21: a multi-state run omits no stage (section 6).

**Stage iteration budget (`MAXITS`): measured slack, not yet reallocated.**
The budget stands at `[20,20,17,15,12,12,12,25]`. A reallocation of stage 3's
budget to stage 6 was made on 2026-10-18 and reverted the same day, unrun,
because moving `NU_FILTER_STAGE` to 5 removes the premise it rested on: stage
6 was starved as the *first* NU stage, and it no longer is. Revisit with
per-stage trajectories measured after the NU-at-5 run. The evidence stands on
its own and is recorded here.

From the ten-restart PfCRT reference `latest6`, which exhausts the budget of
every stage 1-7 in every restart, so no stage before the last converges and
all are cut off rather than finished:

- Stage 3's orientation overlap peaks at iteration 4 (0.50-0.52) and then
  declines monotonically for the remaining 13 (to 0.43-0.45) with resolution
  flat at 11.089 A and cFAR trendless. Churn, not refinement: the band is
  too coarse to discriminate, so assignments wander.
- Stage 6, where NU filtering starts, carries the run from 9.70 to 4.50 A in
  12 iterations and was still gaining resolution at its final iteration in 4
  of the 10 restarts (gains of 0.32-0.58 A at cutoff).

A stage-3-to-6 move would have been roughly time-neutral: -9 iterations at
9.01 s and +8 at 13.26 s, measured per-iteration partition times. Note that
flat FSC=0.143 alone does not mean
a stage is wasteful, because the resolution estimate is quantised to Fourier
shells: stages 4, 5 and the 19-iteration tail of stage 8 all hold resolution
flat while cFAR climbs (0.48-0.57, 0.56-0.68 and 0.59-0.68), so their budgets
stand. `MAXITS_CAVGS_EARLY` still overrides stages 1-3 on the class-average
route, so the stage-3 change applies to the particle route only.

**Mid-stage resolution exit (2026-10-18).** The stage-shift test is also
evaluated during the stage, as a convergence criterion: the internal `lpexit`
key carries `lpinfo(istage+1)%lp`, the next rung's planned limit, and
`check_conv3D` ends the stage once the best populated state's `res05`
resolves past it **and has stopped gaining**: `res05` must be no better than
the previous iteration's, read back from `STATS_FILE` before the iteration
overwrites it (the best of its `RESOLUTION_FSC05_STATEnn`, the same quantity),
and the comparison is within the stage, so the exit never fires at a stage's
first iteration nor against stats whose `ITERATION` precedes the stage.
Logged as `>>> STAGE RESOLUTION EXIT`. Until 2026-10-19 the previous value
was never read back, so the guard could not fire: PfCRT stage 5 ran all 12
`MAXITS` iterations (about 108 s each) where the `latest9` restarts had left
after 3-8.

The gain guard was added after PfCRT `latest9` (9/10 converged, best
3.981 A, run at `8c8d8b42`): the threshold alone ended the NU stage mid-climb
in **all ten** restarts, with last-iteration gains of 0.2-0.8 A, and where it
cut predicted the outcome. The six restarts that left stage 5 at 7.573 A
finished at 4.03-4.14, the two that left at 7.762 A at 3.981 and 4.032, and
the two that left coarsest at 8.171 A gave the two worst maps, 4.140 and the
run's only failure at 5.353. A stage still gaining resolution is productive
whatever the next rung's limit says. Emitted for the stages past stage 2 other
than the entry stage, the last stage and, when an axis is searched,
`SYMSRCH_STAGE`, and it only ever adds convergence, so
`minits` and the Cartesian rule still hold. Single-state runs only since
2026-10-21: the controller emits no `lpexit` for a multi-state run (section 6).

A first version tested against `params%lp`, the band in use, and was
reverted the same day: PfCRT `latest8` converged in 2 of 6 restarts. **Do not
test the band in force.** A band-limited map's FSC=0.5 crossing sits at its
band, because the map resolves to what it was matched at, so the test fires
at the first opportunity in nearly every stage -- the run logged exits at
12.937 A against a band of 12.938 A, and 11.089 against 11.089. Stages 3 and
4 got 3 of their 17 and 15 iterations in every restart, and the final
FSC=0.143 ranged over 4.03-7.96 A against a 4.05 A ten-restart baseline.

The next rung's limit carries information because it has not been applied
yet. On PfCRT's own ladder (20.0, 19.4, 18.3, 12.9, 11.1, 9.7, 9.1, 6.0 A)
against its measured `res05` trajectory (14.8, 12.9, 12.4, 11.5, 9.7 A at
stages 3-7) the test declines at every stage, which is correct: that map does
not run ahead of its schedule. It is expected to fire on datasets whose maps
do, which is where the saving is.

Before 2026-10-18 the band was set to `res05` itself,
which made it a function of a single crossing: on PfCRT, whose FSC runs
nearly parallel to 0.5 over 7.0-4.7 A, a 0.03 FSC difference moves that
crossing 1.5 A, so the mid-stage band varied over 4.6-7.1 A between
restarts for reasons unrelated to the map. One rung bounds what a noisy
crossing can cost. The advance does not accumulate: each stage starts from
its own planned rung, so the lead is at most one rung at any boundary. The
resolution fields are cleared when the particles are reset from `ptcl2D`, so
only resolutions measured in the run itself advance the band. The advanced
value
is the printed stage limit and, in non-NU stages, the `lpstop` ceiling. The
decision is taken once per stage boundary and never per iteration:
`solve3D` runs without gold-standard halves, so an FSC crossing is
trustworthy only where it lies beyond the band that produced the alignments.
The crossing at the end of the previous stage lies beyond that stage's band
and is clean; a per-iteration rule would ratchet on noise fitted inside the
newly opened band. An explicit command-line `lp` (with `ml_reg=yes`) disables
the advance for that stage. Multi-state follows the standard single-band
rule: the best resolved populated state sets the band for all states.
(Streptavidin log set 2026-09-06: the plan sat at 8.6/7.6 A in stages 4/5
while the half maps agreed to 4.3 A at FSC=0.5.)

In the NU stages there is no ceiling (July 2026 policy, restored
2026-09-08): matching runs at the finest member of the NU bank, the finest
rung of the static ladder cut at `fsc/1.5`, or the ML-regularized pair once
its FSC=0.143 is at or beyond the ladder's 4 A rung (2026-09-19, the same
rule as every other workflow; `doc/policies/NU/nonuniform_filtering_policy.md`
sections 8 and 12). The band therefore leads the FSC by 1.25-1.5x, which is
what carries the climb: with the band at the previous FSC=0.143 crossing the
alignment was starved of the shells where the merged reference still holds
signal and PfCRT converged in 0/4 and 1/9 restarts (2026-09-16/18), with a
1.5x lead in 10/10 (2026-09-16), and a ceiling only pins the map. Two
ceilings were tried and
retired: the class-FRC final limit `lpfinal` (6.0 A on PfCRT, whose 2D
classes stop at 6 A while the 3D map reaches 4 A) pinned the NU stages at
5.97 A on 2026-09-07; the ladder's hard bound of 4.5 A pinned the 2026-09-08
run at exactly 4.50 A for 30 iterations while the handoff asked for 4.14 A
(stage 7) and 3.98 A (stage 8), the values the July runs matched at on their
way to 4.1-4.3 A with side chains. Only an explicit command-line `lpstop`
remains a ceiling in NU stages.

Early stopping applies in every stage (the NU-stage `minits=maxits` of
2026-09-07 was retired on 2026-09-08). The early stops that motivated it
happened under a matching-band ceiling: pinned at 5.97 A the sampler
converged on a coarse solution within a few iterations while the FSC could
still improve. Without a ceiling the overlap tracks the map: the 2026-09-08
PfCRT run sat at 0.11-0.25 through stages 6-7 and passed 0.9 in stage 8 only
after the FSC had been flat at the band for ten iterations, where the forced
remainder of the budget changed nothing; on streptavidin the forced budget
cost 700 s of a 2000 s run.
An explicitly supplied, coarser command-line `lpstop` is folded into the
ladder and is also retained as an independent ceiling when the staged child
command is rebuilt; the effective ceiling is the coarser of the two limits.
The workflow logs the acknowledged command-line ceiling before entering the
stage loop. (Record 2026-09-06: capping the NU stages at the per-stage value
stalled every NU stage on streptavidin and msp1 on the pcg path; see
`doc/implementation_notes/completed/pcg_priors_history.md`, dev item 2.)

Saved `_stageNN_lp.mrc` diagnostic volumes are filtered to the current state
FSC resolution when an FSC exists. The planned stage LP is only a fallback.

## 5. Initialization Modes

`solve3D` supports these model-start routes:

- random starting volumes
- user-supplied input volumes
- initialization from `solve3D_cavgs` inside the workflow through
  `cavg_ini=yes`
- externally supplied class-average initialization through `cavg_ini_ext=yes`

Volume input is allowed for single- and multi-state runs (one volume per
state). It cannot be combined with class-average initialization or
partitioned startup. User-supplied input volumes are assumed to be aligned to
the target symmetry axis, so `pgrp_start` is set to `pgrp` and the particle
workflow does not run symmetry-axis search on them.

Input volumes are not assumed to share the particles' alignment provenance.
Before the Euclidean stage ladder, this route calls the shared fixed-reference
CC pose-initialization service. It first runs the normal native-grid
particle-image sigma2 bootstrap, then performs one pass over at most 100,000
active particles at a common 15 A limit. In selected particle records, the
pass replaces active matching shells with reference-conditioned residual sigma2,
consolidates those residuals for the first Euclidean stage, and reconstructs a
data-derived checkpoint. Particles outside the capped cohort retain their
image-bootstrap partition records until later refinement updates them. The
external maps remain fixed during the CC pass and are not blended into the
checkpoint.

Normal particle-based starts treat `solve3D` as the producer of new
`ptcl3D` orientation and multi-state information. The workflow resets `ptcl3D`
sampling, deletes previous 3D alignment while preserving shifts, transfers 2D
shifts from `ptcl2D`, and initializes `ptcl3D%state` only from the 2D
selection state: selected particles become state 1 and unselected particles
become state 0. Fresh multi-state runs then randomize active particles into
the requested 3D states with balanced uniform labels.

Class-average initialization and external class-average initialization both
skip the random-volume start. With `cavg_ini=yes`, the nested
`solve3D_cavgs` run owns any `pgrp_start` to `pgrp` symmetry-axis search.
When control returns to the particle workflow, `pgrp_start` is set to `pgrp` so
the axis search is not repeated. The nested run has as many stages as the
symmetry-search stage (`solve3D_nstages_ini3D`, four), so it ends with the axis
search on the class averages whenever `pgrp_start != pgrp` (section 7). The
particles then enter at that last class-average stage, stage 4, which runs
its full budget (`overlap` 0.99, `fracsrch` 99) on the particles, whatever the
point group (2026-10-17; they entered at stage 3 before). `cavg_ini_ext=yes` is the explicit
exception to the fresh-start rule: it requires prior `ptcl3D` alignment,
preserves the external orientation/state information needed by that route,
assumes the input orientations are already symmetrized, and starts at stage 4,
after the stages a run without an axis search keeps in `pgrp_start`. A prior
`solve3D_cavgs` run with `pgrp_start != pgrp` and fewer than 4 stages, which
would stop before its axis search, is refused by `solve3D_cavgs`. If `nstates > 1`, every requested prior `ptcl3D` state must exist and be
populated.

## 6. Multi-State Policy

`nstates` alone selects the mode; `solve3D` has no `multivol_mode` input and
refuses it on the command line. `nstates=1` is the single-state run.
`nstates > 1` refines independent states from the start: every state has its
own starting model and its own refinement. There is no docked mode (one
consensus model up to a split stage, then states); it was removed on
2026-10-06. To split a finished single-state solution into states, run
`refine3D_states` on its project, where flex PCA initializes the states.

A multi-state run preserves the staged point-group policy between
`pgrp_start` and `pgrp`. When `pgrp_start != pgrp`, the symmetry-search stage
searches the symmetry axis independently for each state, matching the
state-wise behavior used by direct `solve3D_cavgs` runs. This applies to
fresh particle starts only; `cavg_ini=yes`, `cavg_ini_ext=yes`, and
user-supplied input volumes are already in the target symmetry frame before
the parent particle workflow resumes. The multi-state run is intended for
severe heterogeneity where early inspection is more valuable than committing
to a longer refinement immediately. Its particle stages 1 and 2 use direct
`refine=shc` rather than the probabilistic-neighborhood startup of the
single-state run; stages 3 to 5 use `refine=prob`, and `refine=prob_neigh`
with `prob_neigh_mode=state` starts at stage 6, the stage single-state runs
switch at (`PROB_NEIGH_REFINE_STAGE`), so the default five-stage run never
reaches it. From 2026-10-16 to 2026-10-21 stage 5 ran `prob_neigh`, because
the dense `prob` table over 2,500 directions of every state was the most
expensive stage of a multi-state run. That was a mistake (ruling of
2026-10-21): stage 5 is where a multi-state run makes its final labels, and
there the neighbourhood search lost states. On a simulated two-state set the
builds with `prob_neigh` at stage 5 drained a state during it or split the
states too slowly, while the build of 2026-10-15 with `prob` at stage 5
completed the separation there (accuracy 0.86 to 0.97). The class-average
route (`solve3D_cavgs` with `nstates > 1`) shares this mode policy and
switches at the same stage; its early stages keep their own class-average
search schedule. The cap of 2,500 directions (`NSPACE_MULTI_MAX`) also limits every
multi-state stage to 2,500 directions: stages 7 and 8, which run 5,000
directions in single state, run 2,500 when a user extends a multi-state run
to them.

Stage 2 of a multi-state particle run ends by the multi-state convergence
rule (every state's orientation overlap above 0.96 and the mean percentage of
the search space scanned above `fracsrch`) or, since the rulings of
2026-10-15, once `shc` scans more than 95% of the search space on average,
whatever the state overlap. The measure is the mean over the sampled
particles of the per-particle `frac`, which the log prints as "% SEARCH SPACE
SCANNED". The controller sets `fracsrch=95` together with the internal key
`scan_exit=yes` on stage 2's `refine3D` command line (an internal parameter
with no interface entry); the convergence check then ends the stage on the
scan alone. Without the key, stage 1, every later stage and every other user
of the multi-state rule keep the rule as it was. An `shc` search scans more
of the references the longer the states settle, so a stage that kept running
until the overlap converged spent most of its time scanning nearly the whole
space. A `solve3D_addon` pass keeps its own early stopping and does not set
the key.

No stage of a multi-state run is skipped or cut short by a resolution rule
(ruling of 2026-10-21): the stage omission and the mid-stage resolution exit of
section 4 apply to single-state runs only (`solve3D_stage_skipping` in the stage
controller), so every multi-state stage runs until its own convergence rule or
its `MAXITS`. Both rules read the best state's FSC=0.5 resolution. At the start
of a multi-state run every state holds a mixture of nearly the same structure,
so the best state resolves early, and the rules removed exactly the `prob`
stages in which the states split. On a simulated two-state set (30,000
particles, two RNA polymerase conformations) the build of 2026-10-20 omitted
stage 3, ended stage 4 after 3 of its 15 iterations and lost the second state
in stage 5 (accuracy 0.57), while the build of 2026-10-15, which ran stages 3
and 4 in full, separated the states during stage 4 (accuracy 0.97). The
stage-boundary FSC=0.5 advance of a stage's limit still applies.

Unless the user overrides them, the commander sets `nstages=5` and
`lpstop=6.0 A`. Stage 5 is NU-filtered (since 2026-10-18, section 8) but runs
no multi-state trailing reconstruction or staged automasking. After
stage 5, the workflow still runs the final original-sampling reconstruction
so the run produces inspectable `rec_final_stateNN` volumes. To improve
particle coverage before that early exit, the multi-state run starts
stochastic balanced sampling at stage 4: the child `refine3D` stages switch
to `greedy_sampling=no` with `frac_best=1.0` from stage 4 onward. This
samples each unit's quota from the full class rather than from a top-ranked
fraction of that class. The outer particle target remains the fixed
`nsample`-derived update fraction while `nsample/active_particles <= 0.9`;
above that threshold the workflow runs full active-particle updates each
stage.

Prior-orientation multi-state refinement belongs in `refine3D_states`, not in
particle-based solve3D startup.

## 7. Symmetry

`pgrp_start` and `pgrp` must be compatible symmetry groups. If the workflow
raises symmetry, the start group must be a subgroup of the target group. If it
lowers symmetry, the target must be a subgroup of the start group and symmetry
randomization may be applied.

The axis is searched after stage 4 (`SYMSRCH_STAGE`; after stage 3 until
2026-10-17, when the stage-3 map was too poor and the search failed too often).
The stages through it reconstruct in `pgrp_start` and run their full iteration
budget (`overlap` 0.99, `fracsrch` 99); the symmetry-search stage defers its
canonical sigma2 commit to after the symmetric reconstruction. A run that
searches no axis (`pgrp_start = pgrp`: every c1 run and the input-volume,
`cavg_ini`, `cavg_ini_ext` and add-on routes) ends that block at stage 3
(`NOSYMSRCH_STAGE`) as before, so its stage 4 early-stops (`overlap` 0.9,
`fracsrch` 90) and its entry stage is unchanged, except that the particle
stages of a `cavg_ini` run enter at stage 4 and run it to its full budget
(section 5). `solve3D_symsrch_stage`
returns the block's last stage for either case; the final-reconstruction rule
(section 9) and the `solve3D_chunks` refusal follow it.

At the symmetry-search stage, symmetry-axis search is state-local for
multi-state runs. Each active state determines its own axis from its current
map and applies that transform only to orientations assigned to that state.

The axis search scores the Fourier shells between its high-pass and its
low-pass. The low-pass is the planned limit of the symmetry-search stage, at
least 12 Å; the high-pass is the run's `hp` (by default half the box), but
never finer than that low-pass, so the band always holds at least the shell
of the low-pass. A stage limit as coarse as half the box (planned from the
class Fourier ring correlations of very few particles) once left the band empty, and the search
stopped on undefined scores; the polar volume correlator now refuses a band
without a shell above the origin (2026-10-15).

After symmetry handling, the selected maps are injected back into the staged
`refine3D` command line and reference-section files are invalidated.

## 8. Filtering and Automasking

Staged `solve3D` uses the one generated-ladder nonuniform competition
when `filt_mode` is NU-enabled (2026-09-16; the same competition as
`refine3D_auto`, the `nu_refine` shell walk is retired). NU always computes its objective over the full spherical
`mskdiam` support. With the default `automsk=no`, no envelope constrains the
static local-resolution field. `automsk=yes` fixes the filter field outside the
density envelope and masks `_nu_filt` references with it. `automsk=nu` uses a
valid current NU-evidence envelope for those roles and density as fallback;
early FSC/PCG consumers use the previous iteration's evidence artifact.

`solve3D` has no gold-standard stage; gold-standard refinement belongs to
`refine3D_auto`. `envfsc` defaults to `no`. `automsk=yes` implies
`envfsc=yes` (policy 2026-09-09), so the two engage together at
`AUTOMSK_STAGE`/`ENVFSC_STAGE`; the controller forces `envfsc=no` before
`ENVFSC_STAGE`, forces it off for the cavgs route, and forwards the requested
or implied value at and after that boundary. With `envfsc=yes`, volume
assembly generates a density envelope from the current merged half maps for
phase-randomized FSC correction and cFAR on gridding (PCG applies no phase
randomization), using `envmsklp` as its smoothing low-pass. With an active
`automsk` the PCG backend uses the same envelope as the support of both solves;
with `automsk=no` it stays on the `mskdiam` sphere whatever `envfsc` says. The
controller forwards `envmsklp` to staged and final reconstruction commands; its default is `ENVMSKLP_DEFAULT` (20 A). The
density-mask dilation is no longer injected as a layer count: every program
applies the shared physical minimum `ENVMSKWIDTH_A_MIN` (7.5 A) at the
sampling the envelope is built at, so the envelope keeps the same physical
width across the autoscaled stages. The controller keeps a scheduled
`lp` on the refine3D command line. From `NU_FILTER_STAGE`, staged
`nonuniform` is promoted to `nonuniform_lpset`, so the NU frontier can feed an
explicit merged-reference LP-set matching run, bounded at the high-resolution
end by the current `lpstages` limit.

Phase randomization corrects mask-induced correlation; it does not make the
non-independent solve3D half maps a gold-standard validation pair. FSC-derived
resolution metadata in this workflow must be interpreted accordingly.

Automasking is opt-in at the public interface and defaults to `no`. Even when
enabled, staged NU-evidence envelope generation starts only once both
`AUTOMSK_STAGE` and the NU filtering stage are active.
Selecting `automsk=nu` therefore requires an NU `filt_mode`; `automsk=tight`
is rejected. Once staged automasking is active, `yes` uses density and `nu`
uses current evidence for assembled references plus the lagged evidence
artifact for early consumers, with density fallback. There is no separate
reference-mask control, so early stages remain spherical.

NU filtering starts at `NU_FILTER_STAGE` = 5 (6 until 2026-10-18),
overridable from the command line by `nu_stage` for testing (0 = the built-in
value; 5 is the earliest accepted). NU filtering needs an alignment good
enough for its local-resolution estimate: below that the estimate is
meaningless, the filtered reference reinforces the bad alignment rather than
correcting it, and the run goes into a hole. `nu_stage` 3 and 4 were tried on
PfCRT and streptavidin and did not work. Stage 5 is one
stage earlier than the `prob` -> `prob_neigh` switch at `PROB_NEIGH_REFINE_STAGE`.
The two were the same stage before, so a reference becoming NU-filtered and a
global search becoming a neighbourhood search arrived together and could not
be attributed separately. Stage 5 also had capacity to spare on the PfCRT
reference: resolution flat at 9.703 A from its third of twelve iterations,
with only cFAR still moving (0.56-0.68). `l_state_continue_mode` already
enabled NU from stage 5 through `TRAILREC_STAGE_SINGLE`, so this is an
exercised path rather than new ground.

The default multi-state (`nstates > 1`) stage limit is stage 5
(`NSTAGES_INDEPENDENT`), so **that final stage is now an NU stage**; before
2026-10-18 default multi-state ran no NU-filtered stage at all. Applied to
multi-state deliberately, for the quality the filtered reference is expected
to give. The envfsc and automasking boundaries are unaffected: both are
`NSTAGES` (8) and remain outside the default multi-state schedule, as does
`TRAILREC_STAGE_MULTI`. Users who override `nstages` past stage 5 re-enter the
rest of the staged NU policy described here.

Detailed NU behavior belongs to
[nonuniform filtering policy](../NU/nonuniform_filtering_policy.md); detailed
automasking behavior belongs to [automasking_policy.md](automasking_policy.md).

## 9. Final Reconstruction

`solve3D` runs a fresh original-sampling reconstruction from selected
particles for full schedules, for every multi-state schedule, and for a
single-state schedule stopped explicitly (`nstages`) at the symmetry-search
stage or later (stage 3; stage 4 when the run searches an axis, section 7); the
run manifest follows it, so such a run can be the frozen input of
`solve3D_addon`, which then runs that stage. Single-state schedules stopped
before it skip this final all-particle reconstruction.
Before a multi-state final reconstruction, a greedy missing-update pass
(`refine=greedy`) assigns every active particle that no stage updated.

The final reconstruction inherits only the scientific reconstruction policy it
needs. It preserves the parent `envfsc` request so the original-sampling half
maps use the same radial FSC/cFAR masking policy, but it does not inherit staged
search or reference-filter controls such as `refine`, `lp` or `gauref`. It
inherits `automsk` from the last stage, so the shipped map is estimated on the
same density or NU support as the refinement iterations. There is no backend
exception: the final reconstruction runs with `filt_mode=none` on both
backends, except that `automsk=nu` keeps the caller's NU `filt_mode`.

If the final stage used `objfun=euclid` and `ml_reg=yes`, final reconstruction
uses the committed sigma2 store of the workflow's particles: one noise spectrum
per even/odd half. A missing or incompatible store is rebuilt through the one
rebuild entry before the regularized map is produced.

On the PCG backend, a final ML-regularized stage uses the ordinary `P_tau`
replay in `bootstrap_rec3D`; the `Q_NU` prior, its calibration pass and its
controllers were removed on 2026-09-06.
When `pcg_solvent=yes`, the staged controller withholds the soft solvent prior
from every stage before stage 8 and enables it only in stage 8. A workflow
that stops before stage 8 does not use the prior.

The final reconstruction does not apply fractional-update sampling or trailing
average blending. Final-map postprocessing is classical, even when staged
refinement used NU-filtered references.

## 10. Outputs

Stage snapshots are written by `simple_solve3D_utils.f90` with `_stageNN`
suffixes and companion `_lp` diagnostics.

Final outputs use the `rec_final_stateNN` naming convention and include raw and
low-pass diagnostic volumes. Final low-pass diagnostic maps use the state FSC
resolution when available, otherwise the supplied final fallback LP.

Every completed run also writes `solve3D_manifest.txt` last, atomically and
never fatally, and registers it in `projinfo` by bare name with its run
identifier (`simple_solve3D_manifest`): the project identity, the particle
source, the solution, the ladder as planned and emitted, the entry inputs,
content digests of the state maps, halves and FSCs, and a cheap committed-header
and global-section identity for sigma2. Sigma particle spectra are not scanned
to write or validate the manifest. A failed manifest write warns and leaves the
run as it is.

## 11. solve3D_addon

`solve3D_addon projfile=<current> projfile_frozen=<solve3D or solve3D_addon run project>`
extends a completed solution with the particles the current project adds (the
cohort); the frozen particles contribute their accumulators unsearched. Its
policy is [solve3D_addon_policy.md](solve3D_addon_policy.md); the design
history is in
`doc/implementation_notes/completed/abinitio3D_addon_mode_proposal.md`. In
short:

- Only a `solve3D` or `solve3D_addon` output whose manifest validates
  against its project is accepted, so add-ons chain; the solution's settings
  are replayed from the manifest and refused on the command line, which
  carries only the 11 add-on inputs.
- Both projects must share one particle index space: row `i` names the same
  image in both wherever both hold a row `i`, and they may differ in size. Rows
  the current project appends past the frozen project's last row (a stream's
  later particle sets) join the cohort; appended rows from a stack the frozen
  project holds are refused. Every frozen particle must lie within the current
  project's rows and be active there. The cohort needs at least 5
  particles per inherited state (a warning below 5 % of the frozen population).
- The execution-environment keys pass through, the stream's persistent-worker
  keys (`worker_server`, `worker_priority`) included.
- The run searches the cohort in the last stage of the base run's planned
  ladder (read from its manifest), after one registration pass, with
  `center=no` and `overlap` 0.95 by default (it applies when that stage is
  stage 3); the stage limits follow the rule above (FSC=0.5 promotion from the
  union FSC, the NU handoff in the NU stages), and the stage logs them next to
  the base run's.
- The run ends with a validation report against the base solution
  (`solve3D_addon_report.txt` in the run directory, and the log): per state
  the FSC verdict (improved or regressed beyond one FSC=0.143 shell, else
  unchanged), the union map's correlation with the base map up to the base
  resolution (docked below 0.9), the cohort-only map against the base map
  with `addon_diag=yes`, and both runs' stage limits. A regression is warned
  about; the result is published all the same.
- The final reconstruction is incremental (`calc_final_rec_addon`): it reads
  the cohort's particles and the carried frozen sums, never the frozen
  particles, seeds no sigma from image power, commits the cohort's new sigmas
  in a residual pass and ships the union map with the frozen rows restored
  (solve3D_addon_policy.md section 8).
- `mkdir=yes` by default; `mkdir=no` is for NICE. When the run
  has completed, the finished project (every particle posed, the union's
  sigma2 state) replaces the current project file; the frozen project is never
  written. The add-on's own manifest is eligible: the output is the frozen
  input of the next add-on.
