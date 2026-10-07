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
- `sigma_est=global`
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
  no class averages; `cavg` is rejected there)

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
stages. Independent multi-state startup defaults to five stages so it stops
before the `prob_neigh` and NU-filtering stages.

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
    stages 3-5, and `prob_neigh` in later user-enabled stages
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
  has settled over stages 6-7 without the prior

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
The final reconstruction at original sampling is one program call,
`bootstrap_rec3D`, which owns the complete sequence: the same image-power
seed, a euclid ML bootstrap map on it, one residual sigma pass
(`refine=sigma`) against that map, canonical commit of the residual generation,
and the shipped euclid ML reconstruction on it. Because the
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

Stage 1 of `solve3D` never runs with a low-pass limit finer than 20 A.
`lpstages` and `lpstages_fast` apply this floor before generating the remaining
low-pass and crop ladder, preserving gradual frequency marching. Explicit
external-volume schedules use `lpstages_setlims` and are unchanged.

The controller passes an `lpstop` ceiling to each staged `refine3D` child
alongside the effective planned matching limit. In the non-NU stages the
ceiling is the stage's matching limit, so matching never exceeds the
printed stage limit.

Stage-boundary FSC=0.5 promotion (2026-09-06, particle route only): past
stage 2, the planned `lpinfo(istage)%lp` is replaced by the project's FSC=0.5
resolution of the best resolved populated state (the per-particle `res05`
field written by the reconstruction) when that is finer, bounded by the
ladder cap (`LPSTOP_BOUNDS(1)`, 4.5 A, or the coarser explicit `lpstop`; the
class-average route keeps its final limit). The resolution fields are
cleared when the particles are reset from `ptcl2D`, so only resolutions
measured in the run itself promote. The promoted value
is the printed stage limit and, in non-NU stages, the `lpstop` ceiling. The
decision is taken once per stage boundary and never per iteration:
`solve3D` runs without gold-standard halves, so an FSC crossing is
trustworthy only where it lies beyond the band that produced the alignments.
The crossing at the end of the previous stage lies beyond that stage's band
and is clean; a per-iteration rule would ratchet on noise fitted inside the
newly opened band. An explicit command-line `lp` (with `ml_reg=yes`) disables
the promotion for that stage. Multi-state follows the standard single-band
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
the axis search is not repeated. `cavg_ini_ext=yes` is the explicit exception to
the fresh-start rule: it requires prior `ptcl3D` alignment, preserves the
external orientation/state information needed by that route, assumes the input
orientations are already symmetrized, and starts after the symmetry-search
stage. If `nstates > 1`, every requested prior `ptcl3D` state must exist and be
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
single-state run; stages 3-5 use `refine=prob`. Direct `solve3D_cavgs`
initialization retains its own class-average search schedule. Unless the
user overrides them, the commander sets `nstages=5` and `lpstop=6.0 A`.
Stage 5 remains in `prob`; it does not enter `prob_neigh`, static NU
filtering, multi-state trailing reconstruction, or staged automasking. After
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

At the symmetry-search stage, symmetry-axis search is state-local for
multi-state runs. Each active state determines its own axis from its current
map and applies that transform only to orientations assigned to that state.

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
phase-randomized FSC correction and cFAR on gridding, and the PCG backend uses
the same envelope as the support of both solves, using `envmsklp` as its
smoothing low-pass. The controller forwards `envmsklp` to staged and final
reconstruction commands; its default is `ENVMSKLP_DEFAULT` (20 A). The
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

The default multi-state (`nstates > 1`) stage limit stops at stage 5, before
this NU-filtering policy is activated. Users who override `nstages` past that
point re-enter the staged NU policy described here.

Detailed NU behavior belongs to
[nonuniform filtering policy](../NU/nonuniform_filtering_policy.md); detailed
automasking behavior belongs to [automasking_policy.md](automasking_policy.md).

## 9. Final Reconstruction

`solve3D` runs a fresh original-sampling reconstruction from selected
particles for full schedules and for every multi-state schedule. Other
explicit early-stop schedules skip this final all-particle reconstruction.
Before a multi-state final reconstruction, a greedy missing-update pass
(`refine=greedy`) assigns every active particle that no stage updated.

The final reconstruction inherits only the scientific reconstruction policy it
needs. It preserves the parent `envfsc` request so the original-sampling half
maps use the same radial FSC/cFAR masking policy, but it does not inherit staged
search, matching-reference automasking, or reference-filter controls such as
`refine`, `lp`, `automsk`, or `gauref`. There is no backend exception: the
final reconstruction runs with `filt_mode=none` on both backends.

If the final stage used `objfun=euclid` and `ml_reg=yes`, final reconstruction
uses compatible grouped sigma estimates when they are local to the workflow.
If needed, it bootstraps sigmas locally before producing the regularized map.

On the PCG backend, a final ML-regularized stage uses the ordinary `P_tau`
replay in `bootstrap_rec3D`; the `Q_NU` prior, its calibration pass and its
controllers were removed on 2026-09-06.
When `pcg_solvent=yes`, the staged controller withholds the soft solvent prior
from stages 3-7 and enables it only in stage 8. A workflow that stops before
stage 8 does not use the prior.

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
source, the solution, the ladder as planned and emitted, the entry inputs and
the digests of the state maps, halves, FSCs and the committed sigma2 state. A
failed manifest write warns and leaves the run as it is.

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
- The run enters at stage 3 of the base run's planned ladder with
  `center=no` and `overlap` 0.95 by default; the stage limits follow the rule
  above (FSC=0.5 promotion from the union FSC, the NU handoff in the NU
  stages), and each stage boundary logs them next to the base run's.
- The run ends with a validation report against the base solution
  (`solve3D_addon_report.txt` in the run directory, and the log): per state
  the FSC verdict (improved or regressed beyond one FSC=0.143 shell, else
  unchanged), the union map's correlation with the base map up to the base
  resolution (docked below 0.9), the cohort-only map against the base map
  with `addon_diag=yes`, and both runs' stage limits. A regression is warned
  about; the result is published all the same.
- The frozen rows are restored before the final reconstruction and the
  cohort-only sigma2 registration is dropped, so the final reconstruction
  bootstraps the union's sigma2 state over every particle exactly as
  `bootstrap_rec3D` does, and ships the union map on it.
- `mkdir=yes` by default; `mkdir=no` is for NICE and the stream. When the run
  has completed, the finished project (every particle posed, the union's
  sigma2 state) replaces the current project file; the frozen project is never
  written. The add-on's own manifest is eligible: the output is the frozen
  input of the next add-on.
