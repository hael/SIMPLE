# Abinitio3D Policy

This document records the current policy for `abinitio3D`, the staged
particle-based ab initio 3D workflow. The base `refine3D` contracts are in
[refine3D_policy.md](refine3D_policy.md); this document describes how
`abinitio3D` configures and chains those stages.

Multi-state matcher and reconstruction lifetimes follow
[separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md](separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md).

## 1. Scope

`abinitio3D` builds initial 3D models from particles by preparing starting
orientations/states, marching through staged `refine3D` runs, optionally
performing symmetry-axis search, and reconstructing final original-sampling
maps.

It owns stage scheduling. It does not own a separate particle matcher or volume
assembly implementation.

### Projection-direction reconstruction

`projrec=yes` enables an experimental compact reconstruction path for the
particle refinement stages. It first assembles raw Fourier numerator and
CTF-squared sums for each discrete projection direction, state, and even/odd
half using the same native-grid 2D Kaiser-Bessel interpolation machinery as
class averaging. Those un-restored sums are inserted directly into the 3D
partial reconstructions with the 3D Kaiser-Bessel kernel. They are never
CTF-density corrected and never transformed through real space between the 2D
and 3D assembly steps. The default is `projrec=no`.

## 2. Defaults

`abinitio3D` sets:

- `objfun=euclid`
- `sigma_est=global`
- `bfac=0`
- `nu_refine=no`

When unset, it supplies:

- `mkdir=yes`
- `overlap=0.95`
- `prob_athres=10`
- `center=no`
- `cenlp` from the ab initio controller default
- `oritype=ptcl3D`
- `pgrp=c1`
- `pgrp_start=c1`
- `filt_mode=nonuniform`
- `automsk=no`
- `gauref=yes`

For `multivol_mode=independent`, it also supplies conservative inspection
defaults when the user has not overridden them:

- `nstages=5`
- `lpstop=6.0 A`

The public `filt_mode` values are `none`, `nonuniform`, and
`nonuniform_lpset`. Automatic low-pass modes `uniform` and `fsc` are rejected
for `abinitio3D`.

## 3. Stage Controller

The stage controller in `simple_abinitio_controller.f90` emits a concrete
`refine3D` command line for each stage. The full particle workflow has eight
stages. Independent multi-state startup defaults to five stages so it stops
before the `prob_neigh` and NU-filtering stages.

Stage policy includes:

- stage-specific `nspace` and `maxits`
- cropped box and sampling from the low-pass plan
- stage 1 starts with `nspace=500`; stages 2 through 4 use `nspace=1000`
- every stage gets its low-pass and crop information independently from the
  normal FRC/input schedule
- mode-specific staged search:
  - `single` and `docked` particle runs use `prob_neigh` with
    `prob_neigh_mode=shc` in stages 1-2, middle `prob`, and late `prob_neigh`
  - `independent` particle runs use direct `shc` in stages 1-2, `prob` in
    stages 3-5, and `prob_neigh` in later user-enabled stages
- `nspace_sub` for `prob_neigh`
- staged point-group policy between `pgrp_start` and `pgrp`
- staged translation limits
- staged ML regularization
- conical FSC regularization by default only while ML regularization is active
- staged fractional update with a fixed `nsample` target while `nsample/active_particles <= 0.9`
- mode-specific stochastic sampling start
- early Gaussian reference filtering
- optional trailing reconstruction by stage and multivol mode
- staged NU filtering from `NU_FILTER_STAGE`
- staged automasking only from `AUTOMSK_STAGE`

When the downscaled particle cache is requested (`cache=yes`), the stage
controller uses the final active entry in the low-pass/downscaling ladder as
`box_crop` for every emitted `refine3D` command. The stable crop keeps the cache
key unchanged, so all eligible stages reuse one cache instead of rebuilding it
at stage boundaries. Stage low-pass limits remain stage-specific. With
`cache=no`, each stage retains its own cheaper crop from the ladder.

The effective crop and its physically equivalent pixel size are shared by
starting-volume generation, matcher reconstruction, split-stage reconstruction,
symmetry-map commands, FSC diagnostics, low-pass snapshots, and project volume
registration. This preserves `box * smpd == box_crop * smpd_crop` across every
particle-to-volume handoff.

The controller still re-stamps the user's cache request onto every child, so a
stage-local uniform fallback does not permanently disable later stages. A stage
that cannot use the cache releases it; a later eligible stage may rebuild it at
the same final crop. If the final crop is the native box, the existing
`box_crop >= box` eligibility rule disables the cache. See
[particle_cache_policy.md](particle_cache_policy.md) for the full cache
contract.

### Full-Sampling Switch

`abinitio3D` now applies a global sampling override when:

`nsample / active_particles > 0.9`

In that regime, the staged controller forces full active-particle updates for
each iteration by suppressing fractional controls (`update_frac`, `nsample`,
`fillin`) in emitted child `refine3D` commands. This also disables trailing
reconstruction for staged `abinitio3D` commands, and startup class-biased
sampling setup is bypassed in favor of all-active sampling.

The emitted child command line owns `startit` and `which_iter` for the current
stage. `refine3D` then treats `maxits` as the run length for that stage.

## 4. Low-Pass and Cropping

`lpinfo(istage)%lp` controls staged search/reference scheduling. Stage limits
are derived from class FRCs by default, with `lpstart`/`lpstop` overrides and a
`force_lp_range=yes` path that uses the requested range directly.

Stage 1 of `abinitio3D` never runs with a low-pass limit finer than 20 A.
`lpstages` and `lpstages_fast` apply this floor before generating the remaining
low-pass and crop ladder, preserving gradual frequency marching. Explicit
external-volume schedules use `lpstages_setlims` and are unchanged.

Saved `_stageNN_lp.mrc` diagnostic volumes are filtered to the current state
FSC resolution when an FSC exists. The planned stage LP is only a fallback.

## 5. Initialization Modes

`abinitio3D` supports these model-start routes:

- random starting volumes
- user-supplied input volumes
- initialization from `abinitio3D_cavgs` inside the workflow through
  `cavg_ini=yes`
- externally supplied class-average initialization through `cavg_ini_ext=yes`

Volume input is allowed for `single`, `independent`, and `docked`
multi-volume modes. It cannot be combined with class-average initialization or
partitioned startup. User-supplied input volumes are assumed to be aligned to
the target symmetry axis, so `pgrp_start` is set to `pgrp` and the particle
workflow does not run symmetry-axis search on them.

Normal particle-based starts treat `abinitio3D` as the producer of new
`ptcl3D` orientation and multi-state information. The workflow resets `ptcl3D`
sampling, deletes previous 3D alignment while preserving shifts, transfers 2D
shifts from `ptcl2D`, and initializes `ptcl3D%state` only from the 2D
selection state: selected particles become state 1 and unselected particles
become state 0. Fresh `independent` runs then randomize active particles into
the requested 3D states with balanced uniform labels.

Class-average initialization and external class-average initialization both
skip the random-volume start. With `cavg_ini=yes`, the nested
`abinitio3D_cavgs` run owns any `pgrp_start` to `pgrp` symmetry-axis search.
When control returns to the particle workflow, `pgrp_start` is set to `pgrp` so
the axis search is not repeated. `cavg_ini_ext=yes` is the explicit exception to
the fresh-start rule: it requires prior `ptcl3D` alignment, preserves the
external orientation/state information needed by that route, assumes the input
orientations are already symmetrized, and starts after the symmetry-search
stage. If `nstates > 1`, every requested prior `ptcl3D` state must exist and be
populated.

## 6. Multi-Volume Policy

Supported `multivol_mode` values are:

- `single`
- `independent`
- `docked`

`single` requires `nstates=1`. `independent` and `docked` require more than
one state.
When the user gives `nstates > 1` and no `multivol_mode`, the commander
defaults to `independent`.

In `independent` mode, the workflow preserves the staged point-group policy
between `pgrp_start` and `pgrp`. When `pgrp_start != pgrp`, the symmetry-search
stage searches the symmetry axis independently for each state, matching the
state-wise behavior used by direct `abinitio3D_cavgs` runs. This applies to
fresh particle starts only; `cavg_ini=yes`, `cavg_ini_ext=yes`, and user-supplied
input volumes are already in the target symmetry frame before the parent
particle workflow resumes. This mode is intended for severe heterogeneity where
early inspection is more valuable than committing to a longer refinement
immediately. Its particle stages 1 and 2 use direct `refine=shc` rather than
the probabilistic-neighborhood startup used by `single` and `docked`; stages
3-5 use `refine=prob`. Direct `abinitio3D_cavgs` initialization retains its
own class-average search schedule. Unless the user overrides them, the
commander sets `nstages=5` and `lpstop=6.0 A`. Stage 5 remains in `prob`; it
does not enter `prob_neigh`, static NU filtering,
independent-mode trailing reconstruction, or staged automasking. After stage 5,
the workflow still runs the final original-sampling reconstruction so the run
produces inspectable `rec_final_stateNN` volumes. To improve particle coverage
before that early exit, independent mode starts stochastic balanced sampling at
stage 4: the child `refine3D` stages switch to `greedy_sampling=no` with
`frac_best=1.0` from stage 4 onward. This samples each class-balanced quota from
the full class rather than from a top-ranked fraction of that class. The outer
particle target remains the fixed `nsample`-derived update fraction while
`nsample/active_particles <= 0.9`; above that threshold the workflow runs full
active-particle updates each stage.

In `docked` mode, the controller starts as one state, runs stages 1-5 as a
single-state ab initio model, then expands to the requested number of states at
the docked split stage. The default split stage is 6, meaning the split occurs
after stage 5. Docked schedules must reach the configured split stage; an
ordinary `nstages` early stop before `split_stage` is rejected rather than
silently producing a single-state result.

The docked split starts a new multi-state update epoch:

- ordinary pre-split stages use the one-state target
  `min(UPDATE_FRAC_MAX, nsample / active_particles)`
- the stage immediately before the split increases the one-state target to
  `min(UPDATE_FRAC_MAX, nstates * nsample / active_particles)` to broaden pose
  coverage before relabeling
- immediately before the pre-split cohort pass, clear `ptcl3D%sampled` and
  `ptcl3D%updatecnt`
- run one class-balanced `refine=prob` assignment pass with `frac_best=1.0`,
  `fillin=no`, `trail_rec=no`, `volrec=no`, and
  `sticky_class_sampling=no`; this pass
  defines a persistent cohort through `sampled > 0`
- in the fractional regime, use a nominal cohort target of
  `round(2.5 * nstates * nsample)`, subject to
  `NSAMPLE_HET_SPLIT_CAP = 100000`, the active population, and the 90-percent
  `UPDATE_FRAC_MAX` limit; the cohort target is never allowed below the
  effective post-split update target
- in the global full-sampling regime, make the auxiliary pass target all active
  particles
- restore `nstates` to the requested value and recompute the fixed post-split
  target `min(UPDATE_FRAC_MAX, nstates * nsample / active_particles)`
- randomize active particles into balanced uniform state labels without
  clearing the cohort metadata
- require each randomized split state to exceed the probabilistic-table minimum
  population threshold
- select one post-split-sized, class-balanced subset from the persistent cohort
- reconstruct state-specific split volumes and halfmaps from exactly that latest
  sampled subset without trailing volume averaging

The cohort pass prepares one-state pose coverage and establishes the only
ordinary post-split sampling pool. Because its reset happens before the pass,
particles outside the cohort retain `sampled == 0` and `updatecnt == 0`.
After that pass succeeds and the requested state count is restored,
`set_cline_refine3D` obtains the cohort state from
`abinitio_docked_cohort_active` and emits the explicit internal child flag
`sticky_class_sampling=yes`. The flag applies only to the class-balanced
`sample4update_class` path, where it maps to `sampled_only`; particles outside
the cohort are ineligible while cohort members are rotated by increasing
`updatecnt`. It does not alter unbalanced or full particle sampling. The matcher
does not infer stickiness from `nstates` or `multivol_mode`. The flag remains
off in the full-sampling regime.
`sampled == max(sampled)` continues to identify the exact current update, while
`sampled > 0` identifies the persistent cohort.

The split stage is a fractional stabilization stage. It uses
`refine=prob_state`, retains the fixed post-split `update_frac` and `nsample`,
and keeps `fillin=no`. From `TRAILREC_STAGE_SINGLE` onward, trailing
reconstruction is enabled in this fractional regime and uses realized per-state
fractions. The default split stage 6 is after that threshold and therefore
trails against the state-specific split halfmaps reconstructed immediately
after relabeling from the first post-split-sized subset; it must not reuse or
blend the pre-split one-state halfmaps. With an uncapped 2.5-times cohort, the
initial realized fraction is approximately `1 / 2.5`, with state-local and
class-quota rounding. An earlier custom split stage does not enable trailing
before `TRAILREC_STAGE_SINGLE`.

The remaining post-split docked neighborhood stages use `refine=prob_neigh`
with `prob_neigh_mode=geom`. This mode selects the geometric neighborhood
containing each particle's previous best projection and evaluates that same
neighborhood for every state. These stages keep the fixed post-split fractional
target and trailing reconstruction.

When `nsample/active_particles > 0.9`, the global full-sampling switch remains
authoritative throughout docked refinement: emitted stage commands omit
`update_frac`, `nsample`, and `fillin`, and trailing is disabled, including at
the split stage.

`input_oris_start` and `input_oris_fixed` are no longer supported by
`abinitio3D`. Prior-orientation multi-state refinement belongs in the explicit
multi-state refinement workflows, not in particle-based ab initio startup.

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

Staged `abinitio3D` uses static discrete-bank nonuniform filtering when
`filt_mode` is NU-enabled. It always emits `nu_refine=no`; high-resolution NU
shell extension is reserved for `refine3D_auto` and explicit base
`refine3D` use.

Because `abinitio3D` currently keeps gold-standard refinement disabled,
`GOLD_STD_STAGE` is off and `envfsc` defaults to `no`. The advanced `envfsc`
control is nevertheless independent of `automsk` and NU filtering. The
controller forces `envfsc=no` before `ENVFSC_STAGE`, forces it off for the cavgs
route, and forwards the requested value at and after that boundary. With
`envfsc=yes`, volume assembly generates a density envelope from the current
merged half maps for phase-randomized FSC correction and cFAR, using
`envmsklp` as its smoothing low-pass. The controller forwards `envmsklp` to
staged and final reconstruction commands; its default is
`ENVMSKLP_DEFAULT` (20 A). The density-mask growth width defaults to
`ENVMSKWIDTH_DEFAULT` (7 voxels) through `binwidth` and is likewise preserved
in staged and final reconstruction commands. The controller keeps a scheduled
`lp` on the refine3D command line. From `NU_FILTER_STAGE`, staged
`nonuniform` is promoted to `nonuniform_lpset`, so the NU frontier can feed an
explicit merged-reference LP-set matching run.

Phase randomization corrects mask-induced correlation; it does not make the
non-independent ab initio half maps a gold-standard validation pair. FSC-derived
resolution metadata in this workflow must be interpreted accordingly.

Automasking is opt-in at the public interface and defaults to `no`. Even when
enabled, staged NU-evidence envelope generation starts only once both
`AUTOMSK_STAGE` and the NU filtering stage are active.
Selecting `automsk=yes` therefore requires an NU `filt_mode`; `automsk=tight`
and non-NU filtering with automasking are rejected rather than producing a
refinement configuration that cannot generate the requested NU envelope.
Once staged automasking is active, matching references use the lagged state NU
envelope; there is no separate reference-mask control (the former `envref`
parameter has been removed), so early stages without automasking simply remain
spherical.

The default `multivol_mode=independent` stage limit stops at stage 5, before
this NU-filtering policy is activated. Users who override `nstages` past that
point re-enter the staged NU policy described here.

Detailed NU behavior belongs to
[nonuniform_filtering_policy.md](nonuniform_filtering_policy.md); detailed
automasking behavior belongs to [automasking_policy.md](automasking_policy.md).

## 9. Final Reconstruction

After the staged `refine3D` loop, `abinitio3D` runs a fresh
original-sampling reconstruction from selected particles for full schedules
and for `multivol_mode=independent` schedules. Other explicit early-stop
schedules skip this final all-particle reconstruction.

For `docked` mode, the workflow first verifies post-split coverage: every
active particle must have `updatecnt > 0` in the multi-state epoch before final
reconstruction is allowed. This prevents final maps from being produced from
state labels that were never refreshed after the split.

The final reconstruction inherits only the scientific reconstruction policy it
needs. It preserves the parent `envfsc` request so the original-sampling half
maps use the same radial FSC/cFAR masking policy, but it does not inherit staged
search, matching-reference automasking, or NU controls such as `refine`, `lp`,
`automsk`, `gauref`, or NU `filt_mode`.

If the final stage used `objfun=euclid` and `ml_reg=yes`, final reconstruction
uses compatible grouped sigma estimates when they are local to the workflow.
If needed, it bootstraps sigmas locally before producing the regularized map.
For the final ML-regularized stage, final reconstruction preserves the
`conical_fsc` policy selected by the parent workflow.

The final reconstruction does not apply fractional-update sampling or trailing
average blending. Final-map postprocessing is classical, even when staged
refinement used NU-filtered references.

## 10. Outputs

Stage snapshots are written by `simple_abinitio_utils.f90` with `_stageNN`
suffixes and companion `_lp` diagnostics.

Final outputs use the `rec_final_stateNN` naming convention and include raw and
low-pass diagnostic volumes. Final low-pass diagnostic maps use the state FSC
resolution when available, otherwise the supplied final fallback LP.
