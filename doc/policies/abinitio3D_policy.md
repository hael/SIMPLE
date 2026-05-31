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

The public `filt_mode` values are `none`, `nonuniform`, and
`nonuniform_lpset`. Automatic low-pass modes `uniform` and `fsc` are rejected
for `abinitio3D`.

## 3. Stage Controller

The stage controller in `simple_abinitio_controller.f90` emits a concrete
`refine3D` command line for each stage. The current particle workflow has eight
stages.

Stage policy includes:

- stage-specific `nspace` and `maxits`
- cropped box and sampling from the low-pass plan
- staged search mode: early `shc_smpl`, middle `prob`, late `prob_neigh`
- `nspace_sub` for `prob_neigh`
- staged point-group policy between `pgrp_start` and `pgrp`
- staged translation limits
- staged ML regularization
- staged fractional update and optional `nsample` ramp
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
- prior-orientation modes through `input_oris_start` and `input_oris_fixed`

Volume input is allowed for `single`, `independent`, and `docked`
multi-volume modes. It cannot be combined with class-average initialization or
partitioned startup.

Class-average initialization and external class-average initialization both
skip the random-volume start, but only the external route assumes the input
orientations are already symmetrized and starts after the symmetry-search
stage.

## 6. Multi-Volume Policy

Supported `multivol_mode` values are:

- `single`
- `independent`
- `docked`
- `input_oris_start`
- `input_oris_fixed`

`single` requires `nstates=1`. The other modes require more than one state.
When the user gives `nstates > 1` and no `multivol_mode`, the commander
defaults to `independent`.

In `independent` mode, the workflow uses the target point group from the start
and skips the staged symmetry-axis search.

In `docked` mode, the controller starts as one state and expands to the
requested number of states at the docked split stage.

Input-orientation modes require a prior `ptcl3D` alignment. They seed or
preserve states according to the selected mode before entering the staged
refinement schedule.

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

Detailed NU behavior belongs to
[nonuniform_filtering_policy.md](nonuniform_filtering_policy.md); detailed
automasking behavior belongs to [automasking_policy.md](automasking_policy.md).

## 9. Final Reconstruction

After the staged `refine3D` loop, `abinitio3D` runs a fresh
original-sampling reconstruction from selected particles.

The final reconstruction inherits only the scientific reconstruction policy it
needs from the final stage. It must not inherit staged search or mask-generation
controls such as `refine`, `lp`, `automsk`, `envfsc`, `gauref`, or NU
`filt_mode`.

If the final stage used `objfun=euclid` and `ml_reg=yes`, final reconstruction
uses compatible grouped sigma estimates when they are local to the workflow.
If needed, it bootstraps sigmas locally before producing the regularized map.

The final reconstruction does not apply fractional-update sampling or trailing
average blending. Final-map postprocessing is classical, even when staged
refinement used NU-filtered references.

## 10. Outputs

Stage snapshots are written by `simple_abinitio_utils.f90` with `_stageNN`
suffixes and companion `_lp` diagnostics.

Final outputs use the `rec_final_stateNN` naming convention and include raw and
low-pass diagnostic volumes. Final low-pass diagnostic maps use the state FSC
resolution when available, otherwise the supplied final fallback LP.
