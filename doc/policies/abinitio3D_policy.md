# Abinitio3D Policy

This document records the current policy for `abinitio3D`, the staged
particle-based ab initio 3D workflow. The base `refine3D` contracts are in
[refine3D_policy.md](refine3D_policy.md); this document describes how
`abinitio3D` configures and chains those stages.

## 1. Scope

`abinitio3D` builds initial 3D models from particles by preparing starting
orientations/states, marching through staged `refine3D` runs, optionally
performing symmetry-axis search, and reconstructing final original-sampling
maps.

It owns stage scheduling. It does not own a separate particle matcher or volume
assembly implementation.

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
- staged search mode: early `shc_smpl`, middle `prob`, late `prob_neigh`
- `nspace_sub` for `prob_neigh`
- staged point-group policy between `pgrp_start` and `pgrp`
- staged translation limits
- staged ML regularization
- conical FSC regularization by default only while ML regularization is active
- staged fractional update and optional `nsample` ramp
- mode-specific stochastic sampling start
- early Gaussian reference filtering
- optional trailing reconstruction by stage and multivol mode
- staged NU filtering from `NU_FILTER_STAGE`
- staged automasking only from `AUTOMSK_STAGE`

The emitted child command line owns `startit` and `which_iter` for the current
stage. `refine3D` then treats `maxits` as the run length for that stage.

## 4. Low-Pass and Cropping

`lpinfo(istage)%lp` controls staged search/reference scheduling. Stage limits
are derived from class FRCs by default, with `lpstart`/`lpstop` overrides and a
`force_lp_range=yes` path that uses the requested range directly.

For external starting volumes, the low-pass plan is derived from the input
volume dimensions and mask diameter.

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
the requested 3D states.

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
immediately. Unless the user overrides them, the commander sets `nstages=5` and
`lpstop=6.0 A`. Stage 5 is still in the `prob` phase; it does not enter
`prob_neigh`, static NU filtering,
independent-mode trailing reconstruction, or staged automasking. After stage 5,
the workflow still runs the final original-sampling reconstruction so the run
produces inspectable `rec_final_stateNN` volumes. To improve particle coverage
before that early exit, independent mode starts stochastic balanced sampling at
stage 4:
`nsample_start` reaches the full `nsample` target by stage 4, and the child
`refine3D` stages switch to `greedy_sampling=no` with `frac_best=1.0` from
stage 4 onward. This samples each class-balanced quota from the full class
rather than from a top-ranked fraction of that class.

In `docked` mode, the controller starts as one state, runs stages 1-5 as a
single-state ab initio model, then expands to the requested number of states at
the docked split stage. The default split stage is 6, meaning the split occurs
after stage 5.

The docked split starts a new multi-state update epoch:

- restore `nstates` to the requested value
- force the split-stage particle update target to `UPDATE_FRAC_MAX`
- clear `ptcl3D%sampled` and `ptcl3D%updatecnt`
- randomize active particles into the requested state labels
- reconstruct split state volumes from the randomized labels without trailing
  volume averaging

The first post-split stage is a stabilization stage. It uses `refine=shc_smpl`
with fractional particle updates still active, but with fractional volume
averaging disabled. The remaining post-split docked stages use
`refine=prob_neigh` and restore trailing reconstruction within the new
multi-state update epoch.

For docked mode, trailing reconstruction is allowed before the split and after
the split stage. The split stage itself must not blend current state volumes
with previous mixed single-state or pre-split volumes.

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
`GOLD_STD_STAGE` is off, `envfsc=no`, and the controller keeps a scheduled
`lp` on the refine3D command line. From `NU_FILTER_STAGE`, staged
`nonuniform` is promoted to `nonuniform_lpset`, so the NU frontier can feed an
explicit merged-reference LP-set matching run.

Automasking is opt-in at the public interface and defaults to `no`. Even when
enabled, staged automasking starts only from `AUTOMSK_STAGE`.

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
needs from the final stage. It must not inherit staged search or mask-generation
controls such as `refine`, `lp`, `automsk`, `envfsc`, `gauref`, or NU
`filt_mode`.

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
