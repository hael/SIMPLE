# Refine3D Multi Policy

This document records the current policy for `refine3D_multi`. It describes
the automated multi-state wrapper around base `refine3D`; the low-level
projection-matching, probabilistic pre-step, partial reconstruction, volume
assembly, automasking, and nonuniform-filtering contracts remain in the
related policy documents.

Related workflow policies:

- [refine3D_policy.md](refine3D_policy.md)
- [refine3D_auto_policy.md](refine3D_auto_policy.md)
- [automasking_policy.md](automasking_policy.md)
- [nonuniform_filtering_policy.md](nonuniform_filtering_policy.md)
- [sigma_calculation_policy.md](sigma_calculation_policy.md)

Primary implementation files:

- `src/main/ui/simple/simple_ui_refine3D.f90`
- `src/main/exec/simple_exec_refine3D.f90`
- `src/main/commanders/simple/simple_commanders_refine3D.f90`
- `src/main/strategies/parallelization/simple_refine3D_strategy.f90`
- `src/main/strategies/search/simple_strategy3D_matcher.f90`
- `src/main/commanders/simple/simple_commanders_rec_distr.f90`

## 1. Scope

`refine3D_multi` is an automated multi-state refinement workflow for projects
with existing input orientations. It chooses multi-state defaults, establishes
the number of states by inquiring either the command line or the project state
labels, prepares or validates state volumes, maps `multivol_mode` onto a staged
base-`refine3D` workflow, and then runs a final all-particle `reconstruct3D`
pass at native project sampling.

It is not a separate matcher implementation. Once the wrapper has prepared the
command line for a stage, all particle-domain and volume-domain iteration work
is delegated to `commander_refine3D`.

Its durable contract is:

1. establish `nstates` from project labels or command-line input
2. apply multi-state defaults
3. derive sampling, update fraction, and per-stage iteration caps
4. derive default autoscaling and translation limits, unless `autoscale=no`
5. validate or initialize per-state starting volumes
6. run the `multivol_mode` stage plan
7. verify that every active particle has been updated at least once
8. reconstruct and postprocess final maps from all particles at native sampling
9. write final reconstruction products

`refine3D_multi` should be reviewed as a staged policy wrapper. Changes to
the base matcher, probability tables, volume assembly, automasking, or NU
filtering should be checked against their own policy documents as well.

The old `multivol_assign` public application has been removed from the
library. Multi-volume assignment and refinement policy lives in
`refine3D_multi` through `multivol_mode`.

## 2. Entry Points and Ownership

The public command is `refine3D_multi`, registered in
`simple_ui_refine3D.f90` and routed by `simple_exec_refine3D.f90`.

`simple_exec_refine3D.f90` owns command dispatch:

- ordinary runs call `commander_refine3D_multi%execute`
- runs with `nrestarts` are delegated to `restarted_exec`

`simple_commanders_refine3D.f90` owns the wrapper policy:

- multi-state default selection
- project-label versus command-line `nstates` policy
- stage planning and stage-level `refine` mode selection
- final `reconstruct3D` command-line preparation
- final output copying through `write_final_rec_outputs`

Base `refine3D` owns each iteration after the wrapper calls it. The wrapper
must not duplicate or bypass:

- reprojection-model materialization
- probabilistic alignment execution
- hard-assignment matcher execution
- partial reconstruction writing
- `volassemble`
- FSC/resolution registration
- per-state automask and NU filtering decisions

## 3. Public Command and Defaults

The UI describes `refine3D_multi` as an automated multi-state 3D refinement
workflow that requires an `sp_project`. It exposes search, filter, mask, and
compute controls:

- search: `maxits`, `nstates`, `pgrp`, `autoscale`, `continue`
- mode: `multivol_mode`
- filter: `filt_mode`, `lpstop`, `ml_reg`
- mask: `mskdiam`, `automsk`
- compute: `nparts`, `nthr`

The commander sets these hard workflow defaults:

- `balance=yes`
- `greedy_sampling=no`
- `frac_best=1.0`
- `trail_rec=yes`, except `input_oris_fixed` sets `trail_rec=no`
- `objfun=euclid`
- `envfsc=no`
- `lplim_crit=0.5`
- `incrreslim=no`
- `nu_refine=no`
- `combine_eo=no`

The commander also sets these defaults only when the user did not provide a
value:

- `mkdir=yes`
- `center=no`
- `sigma_est=global`
- `prob_inpl=yes`
- `prob_neigh_mode=sum`
- `nsample=min(100000, 10000 * nstates)`
- `autoscale=yes`
- `multivol_mode=input_oris_start`
- `ml_reg=yes`
- `filt_mode=nonuniform_lpset`
- `lpstop=6.0`
- `automsk=no`
- `overlap=0.99` for the `prob_neigh` phase
- `keepvol=no`

The public `filt_mode` values are `fsc`, `nonuniform_lpset`, and `none`;
`nonuniform_lpset` is the default. Plain `nonuniform` is not a
`refine3D_multi` mode because the workflow is LP-set driven, and `uniform` is
disabled for multi-state search.

`nu_refine` is not a public `refine3D_multi` option. The UI does not expose it,
and the commander forces `nu_refine=no` even if a command line attempts to pass
another value.

`combine_eo` is not a public `refine3D_multi` option. The UI does not expose it,
and the commander forces `combine_eo=no`; `combine_eo=yes` is rejected because
combined even/odd terminal alignment is gold-standard-style base `refine3D`
policy, not multi-state LP-set refinement policy.

The UI defaults should stay aligned with those commander defaults. If a UI
default changes without a commander default change, batch and GUI behavior can
diverge.

## 4. Mode and State Model

`refine3D_multi` is always multi-state. The effective `nstates` is established
before the wrapper parses the full parameter object, and `multivol_mode`
decides how project state labels are interpreted.

Supported `multivol_mode` values are:

- `input_oris_start`
- `input_oris_refine`
- `input_oris_fixed`

`input_oris_start` is the default. It is the full workflow:

- state 0/1 input with command-line `nstates`: run `prob_state`, then
  `shc_smpl`, then `prob_neigh`
- existing project state labels greater than one with no command-line
  `nstates`: reuse those labels and skip `prob_state`
- existing project state labels greater than one with command-line `nstates`:
  fail, because command-line `nstates` means the caller expects a state 0/1
  project and a fresh state split

`input_oris_refine` is the refinement-from-input-orientations route:

- state 0/1 input requires command-line `nstates`
- state 0/1 input runs `prob_state`, then goes directly to `prob_neigh`
- existing project state labels greater than one are reused as the starting
  point for `prob_neigh`, with no `prob_state` initialization
- command-line `nstates`, when supplied for an already multi-state project,
  must match the existing project state count

`input_oris_fixed` is the fixed-orientation state-separation route:

- input is expected to contain only `state=0` and `state=1` particles
- command-line `nstates` is required and must be greater than one
- existing project state labels greater than one are rejected; the command
  fails before refinement because `input_oris_fixed` is only valid for state
  0/1 input
- the only refinement stage is `prob_state`
- the assignment path updates state/correlation/accounting fields but does not
  update Euler angles, projection direction, in-plane angle, or shifts

For all modes, every accepted populated project state from `1..nstates` must
have at least one active particle. `state=0` remains inactive. Active-particle
counts use `state > 0` when labels exist. In state 0/1 initialization mode, the
active count honors existing selection state when possible and otherwise falls
back to all `ptcl3D` rows.

## 5. Starting State Volumes

State volume policy is all-or-none.

If `vol1..volN` are all supplied:

- every referenced file must exist
- every volume must match the native project particle box
- every volume must match the native project particle sampling
- those volumes are used as state starting references

If some, but not all, `vol1..volN` keys are supplied, the wrapper fails. Partial
state-volume input is ambiguous and must not be accepted.

If no complete command-line volume set is supplied, compatible project output
volumes are preferred:

- `os_out` must contain a `vol` entry for every state
- every referenced file must exist
- every output volume must match the native project box
- every output volume must have positive sampling matching the native project
  sampling
- compatible project volumes are injected back as `vol1..volN`

If no compatible volume source is available and the project already has
multi-state labels, the wrapper runs a startup `reconstruct3D` pass with
`postprocess=no`, `nu_refine=no`, `objfun=cc`, and no crop overrides. The
resulting `vol_stateNN.mrc` products become the starting `volN` references.

Startup volume selection does not silently re-filter user or project volumes.
Compatible `os_out` volumes are accepted as-is and then become the references
for the first staged base-`refine3D` call. If startup volumes are generated
from pre-existing project labels, the startup reconstruction is native,
unpostprocessed, and non-NU. The first real refinement stage is where
`volassemble` may apply the normal multi-state reference policy:
`ml_reg=yes`, `filt_mode=nonuniform_lpset`, `nu_refine=no`, and any selected
LP-set handoff owned by volume assembly.

When no compatible volume source is available for state 0/1 initialization and
distributed execution is requested with `nparts`, the wrapper leaves
`vol1..volN` absent and delegates startup reference preparation to base
`refine3D`. The first `prob_state` stage seeds fresh state labels in the base
strategy and reconstructs state-specific `startvol_stateNN.mrc` references
from those random state groups before particle-domain probability-table work
begins. This path requires previous objective-function values when no complete
state-volume source is supplied, following the base `refine3D` startup guard.

Shared-memory base `refine3D` still requires explicit starting volumes. Until
that path gains the same internal startup reconstruction parity, state 0/1
initialization without `vol1..volN` must run distributed by setting `nparts`.

## 6. Sampling and Update Fraction

The default sample target is `10000` active particles per state, capped at
`100000` active particles per iteration:

```text
nsample = min(100000, 10000 * nstates)
```

`nsample` must be at least one.

When active particles do not exceed the sample target:

- `update_frac` is removed from the command line
- effective update fraction is `1.0`
- fractional-update mode is disabled
- trailing reconstruction is disabled in the effective parameter state

When active particles exceed the sample target:

- `update_frac = nsample / active_particles`
- fractional-update mode is enabled only if the fraction is at most `0.99`
- trailing reconstruction follows the `trail_rec` setting only when
  fractional-update mode is enabled
- near-full fractions above `0.99` are treated as full updates

`input_oris_fixed` always runs full updates. The wrapper removes
`update_frac`, disables fractional-update mode, and disables trailing
reconstruction for this mode even when the active particle count exceeds the
ordinary sample target.

Fractional multi-state sampling is class-balanced by default:

- `balance=yes`
- `greedy_sampling=no`
- `frac_best=1.0`

This uses the prior 2D class analysis only to define class-biased sampling
quotas. It does not sample greedily with respect to 2D objective-function
values, and `frac_best=1.0` leaves all particles in each selected class eligible
for balanced coverage sampling. Within each class quota, the sampler prefers
the lowest `updatecnt` tier first and randomizes only among particles tied at
the same update count. When fractional balanced updates are active, the wrapper
writes the same class-sampling sidecar consumed by base `refine3D`. If the
project lacks selected `cls2D` entries, the run fails early with an explicit
class-balanced sampling error.

The wrapper does not set `ufrac_trec`; trailing reconstruction consumes the
realized per-state update fractions recorded in `sampled` and `updatecnt`.

Automatic stage budget planning targets roughly four updates per active
particle:

```text
ceil(4.0 * active_particles / particles_per_iteration)
```

The automatically planned stage cap is clamped to:

- at least `10` iterations
- at least the stage-2 minimum of `5` iterations
- at most `50` iterations

If the user supplies `maxits`, that value is interpreted as a stage cap for
`refine3D_multi`, not as the total wrapper iteration count.
`input_oris_fixed` accepts `maxits >= 1` because it has no `prob_neigh` stage.
The other modes require `maxits >= 5` because `prob_neigh` has a minimum of
five iterations. In `input_oris_start` and `input_oris_refine`, the
state-initialization `prob_state` phase remains capped at five iterations.

After all stages finish, the wrapper reads `ptcl3D` and verifies update
coverage. The run is not allowed to continue to final reconstruction unless
every active particle has `updatecnt > 0`, meaning every active particle has
seen the multi-state model and had its state assignment and, for non-fixed
modes, orientation-accounting path updated at least once.

## 7. Autoscaling and Matching Bandwidth

`autoscale=yes` is the default. It uses a two-stage low-pass planning helper
with:

- low-pass start: `10.0 A`
- low-pass stop: `lpstop`, default `6.0 A`

The second planned stage controls the active translation limit and any crop:

- `trs` is set from the stage-2 planning result
- if autoscaling is selected, `box_crop` and `smpd_crop` are set from the
  stage-2 crop
- if autoscaling is not selected, crop overrides are removed

The wrapper leaves final reconstruction uncropped by deleting `box_crop` and
`smpd_crop` before the final all-particle `reconstruct3D` pass, and it writes
final output snapshots with native project box and sampling.

The default `filt_mode=nonuniform_lpset` means volume assembly may run ordinary
NU filtering and may promote the selected NU bandwidth into LP-set matching.
`nu_refine` is always `no`, so high-resolution NU shell expansion is not part
of multi-state refinement.

## 8. Stage Plan

After startup, the wrapper changes `prg` from `refine3D_multi` to `refine3D`
and invokes the base commander once per enabled stage. Each call uses base
`refine3D`'s normal staged iteration policy: `maxits` is the run length for
that stage, `minits` is the stage minimum, and the base strategy reports the
actual final iteration through `endit`.

The possible stages are:

1. state initialization
   - `refine=prob_state`
   - `nspace=2500`
   - minimum iterations: `1`
   - maximum iterations: `5` for initialization, or the stage cap when
     `input_oris_fixed` makes this the only stage
2. stochastic orientation refinement
   - `refine=shc_smpl`
   - `nspace=2500`
   - minimum iterations: `1`
   - maximum iterations: the planned or user-supplied stage cap
3. probabilistic neighborhood refinement
   - `refine=prob_neigh`
   - `nspace=5000`
   - `nspace_sub=500`
   - `prob_neigh_mode=sum`
   - minimum iterations: `5`
   - maximum iterations: the planned or user-supplied stage cap

The `multivol_mode` stage map is:

- `input_oris_start` with state 0/1 input: `prob_state`, `shc_smpl`,
  `prob_neigh`
- `input_oris_start` with existing multi-state labels: `shc_smpl`,
  `prob_neigh`
- `input_oris_refine` with state 0/1 input: `prob_state`, `prob_neigh`
- `input_oris_refine` with existing multi-state labels: `prob_neigh`
- `input_oris_fixed`: `prob_state`

For every stage-level base-`refine3D` call, the wrapper sets:

- `maxits` to the stage cap
- `minits` to the stage minimum
- `startit` to the next global wrapper iteration number
- `which_iter` to the same stage-start global iteration number
- `extr_iter` to the same stage-start global iteration number
- `refine` to the stage mode
- `nspace` to the stage projection count
- `nspace_sub=500` for `prob_neigh`, and no `nspace_sub` for other stages
- `prob_neigh_mode=sum` by default, unless explicitly overridden
- `overlap` to the stage state-overlap target

The wrapper deletes stale `endit` before each staged call and reads the new
`endit` afterward to keep global wrapper iteration numbers monotonic across the
whole run. Iteration numbering does not reset per stage.

`maxits_glob` is set before the stages so annealing and regularization logic
can see the intended global wrapper horizon. It is the sum of the caps for the
stages that will actually run.

## 9. Stage Stop Criterion

For `input_oris_*` multi-state modes, the base convergence signal is
state-label overlap (`mi_state`). Other multi-state base `refine3D` workflows
keep their existing joint-orientation-overlap convergence.

Each stage passes an overlap target to base `refine3D`. After the stage call
returns, the wrapper also reads the project `ptcl3D` segment and reports the
final mean `mi_state` over active particles:

- if `mi_state` is absent, overlap is `0`
- if `sampled` exists, only the latest sampled subset is counted
- otherwise all active state-labelled particles are counted
- particles with state less than or equal to zero are excluded

The stage targets are:

- `prob_state`: `0.95`
- `shc_smpl`: `0.95`
- `prob_neigh`: `overlap`, default `0.99`

A stage may stop only after its minimum iteration count has been reached and
the measured state overlap is greater than or equal to the stage target.

If the stage reaches its cap before the target, the wrapper logs that the cap
was reached and continues to the next stage.

State overlap measures state-label stability. It is not a substitute for FSC,
map quality, pose overlap, or resolution convergence, and final workflow exit
still requires the active-particle update coverage guard described above.

## 10. Delegation to Base Refine3D

Every stage is a normal base `refine3D` invocation after the wrapper has set
the per-stage command-line values. The base command then selects shared-memory
or distributed execution:

- `nparts` without `part` selects the distributed master strategy
- otherwise the shared-memory strategy is used
- distributed workers run through the private worker route with `part` and
  `outfile`

Base `refine3D` owns:

- reference section materialization from current state volumes
- probability-table generation for `prob_state` and `prob_neigh`
- particle sampling for each update
- hard-assignment search and state updates
- `mi_state`, `mi_proj`, `corr`, shift, and orientation metadata updates
- Euclidean sigma updates when `objfun=euclid`
- partition-local partial reconstruction output
- `volassemble`
- per-state FSC and resolution metadata
- state volume registration in `os_out`

All multi-state base calls run with shift-first candidate scoring disabled.
Shifts may still be refined after the candidate/state choice, but a shift seed
estimated from a particle's previous state must not be used to rank another
state.

The wrapper must preserve the base `refine3D` particle-domain versus
volume-domain boundary. It should not call matcher internals or volume
postprocessing helpers directly.

## 11. Filtering, Masking, and Volume Assembly

The default filtering policy is optimized for multi-state refinement
references:

- `filt_mode=nonuniform_lpset`
- `nu_refine=no`, forced and not exposed as a `refine3D_multi` option
- `ml_reg=yes`
- `automsk=no`
- `envfsc=no`
- `lplim_crit=0.5`

`volassemble` remains the execution site for volume-domain work. In a default
run it may restore per-state half maps, calculate FSCs, run ML regularization,
run ordinary NU filtering, write NU-derived reference products, and hand an
LP-set matching bandwidth back through the project. It must not perform
`nu_refine` shell expansion for `refine3D_multi`.

`automsk=no` means the wrapper does not request automatic state-mask
generation by default. If the user enables `automsk=yes` or `automsk=tight`,
state-specific mask production and compatibility rules follow
[automasking_policy.md](automasking_policy.md).

When NU filtering is active, support-mask selection and matching-bandwidth
handoff follow [nonuniform_filtering_policy.md](nonuniform_filtering_policy.md).
`refine3D_multi` must not add a second source of NU handoff state.

The workflow is LP-set multi-state refinement, not gold-standard-style terminal
alignment. `combine_eo` is therefore not exposed and `combine_eo=yes` is
rejected before the wrapper calls base `refine3D`.

## 12. Final Reconstruction

After all staged refinement iterations complete, the wrapper runs
`reconstruct3D` from all particle images. This pass is always at native project
box and sampling, mirroring the final native reconstruction policy used by
`refine3D_auto`.

The final reconstruction command line:

- sets `prg=reconstruct3D`
- sets `outfile=RESOLUTION_FINAL.txt`
- sets `postprocess=yes`
- deletes `trail_rec`
- deletes `refine`
- deletes `sigma_est`
- deletes `update_frac`
- deletes `ufrac_trec`
- deletes the internal stage handoff key `endit`
- deletes `box_crop`
- deletes `smpd_crop`
- sets `objfun=cc`
- sets `nu_refine=no`
- sets `filt_mode=none` when the staged refinement used NU filtering
- sets final-output bookkeeping to the native project `box` and `smpd`

This final pass is an ordinary global final reconstruction/postprocess route.
NU filtering is a refinement-reference feature here, not a separate final-map
postprocessing path.

Before this final pass, the wrapper verifies that every active particle has
`updatecnt > 0`. After final reconstruction, `write_final_rec_outputs` copies
each populated state's `vol_stateNN.mrc` to the final reconstruction naming
convention, writes a low-pass snapshot using `lpstop` at native sampling, and
copies available postprocessed and mirrored products.

## 13. Continue, Restarts, and Execution Mode

The UI exposes `continue`, and the router supports `nrestarts` through the
same restarted-execution mechanism used by `refine3D` and `refine3D_auto`.

The wrapper has no independent resume implementation. It passes the staged
command line into base `refine3D`, so resume behavior is constrained by the
selected base strategy:

- shared-memory base `refine3D` rejects `continue=yes`
- distributed base `refine3D` owns continue checks for previous volumes, FSCs,
  partition counts, and sigma files

A reviewer should treat `continue=yes` as a base-`refine3D` contract that the
wrapper must not weaken. Wrapper startup volume selection should not hide
missing or incompatible resume artifacts that the base distributed strategy is
supposed to validate.

## 14. Artifacts

The staged refinement uses the ordinary `refine3D` artifact names:

- state volumes: `vol_stateNN.mrc`
- half maps: `vol_stateNN_even.mrc`, `vol_stateNN_odd.mrc`
- unfiltered half maps: `vol_stateNN_even_unfil.mrc`,
  `vol_stateNN_odd_unfil.mrc`
- FSCs: `fsc_stateNN.bin`
- per-iteration FSC plots: `fsc_stateNN_iterIII.*`
- partial reconstructions: `vol_stateNN_partPP_even.mrc`,
  `vol_stateNN_partPP_odd.mrc`
- partial density weights: `rho_vol_stateNN_partPP_even.mrc`,
  `rho_vol_stateNN_partPP_odd.mrc`
- reprojection models: `reprojection_model_even.bin`,
  `reprojection_model_odd.bin`
- matcher benchmarks: `REFINE3D_BENCH_ITERIII.txt`
- strategy benchmarks: `REFINE3D_STRATEGY_BENCH_ITERIII.txt`
- assembly benchmarks: `VOLASSEMBLE_BENCH_ITERIII.txt`

The final reconstruction writes ordinary reconstruct3D products and
`RESOLUTION_FINAL.txt`, then `write_final_rec_outputs` writes final state-map
copies and low-pass snapshots for populated states.

The wrapper's `cleanup_init_vols` releases local string objects only. It does
not delete user-supplied or project-derived starting volumes.

## 15. Future Implementation Plan

`input_oris_fixed` currently performs full reconstruction updates and never
uses fractional updates. A future optimization should make fixed-orientation
state separation faster without changing the scientific contract:

- keep unnormalized per-state Fourier volumes resident in memory
- cache each particle representative in the same Fourier/gridding convention
  used for reconstruction
- when a particle changes state, subtract its cached representative from the
  old state accumulator and add it to the new state accumulator
- update density weights and even/odd accumulators consistently with the
  existing reconstruction path
- preserve the final native all-particle `reconstruct3D` pass as the
  authoritative final map generation step

This optimization should be implemented behind the same `input_oris_fixed`
workflow contract: state assignment may change, orientations must remain
fixed, and the run must still verify that every active particle was updated.

## 16. Reviewer Checklist

When reviewing `refine3D_multi`, check these policy points first:

- UI defaults and commander defaults still match.
- `multivol_mode` accepts only `input_oris_start`, `input_oris_refine`, and
  `input_oris_fixed`.
- `nstates` is established before full parameter parsing.
- `input_oris_start` with command-line `nstates` rejects existing project
  multi-state labels.
- `input_oris_start` without command-line `nstates` reuses existing
  multi-state labels.
- `input_oris_refine` skips `prob_state` when project multi-state labels
  already exist.
- `input_oris_fixed` rejects existing project multi-state labels.
- state 0/1 initialization modes require command-line `nstates > 1`.
- state 0/1 initialization modes use complete `vol1..volN` inputs, compatible
  project state volumes, or distributed base-`refine3D` startup reconstruction
  from freshly seeded state labels.
- volume-free state 0/1 initialization requires distributed execution
  (`nparts`) because shared-memory base `refine3D` still requires starting
  volume input.
- fixed-orientation assignment updates state/correlation/accounting fields
  without changing Euler angles, projection direction, in-plane angle, or
  shifts.
- partial `vol1..volN` input is rejected.
- input and project state volumes are checked against native box and sampling.
- default `nsample` is `10000 * nstates` capped at `100000`.
- public `filt_mode` values are `fsc`, `nonuniform_lpset`, and `none`, with
  `nonuniform_lpset` as the default.
- fractional updates use `balance=yes`, `greedy_sampling=no`, and
  `frac_best=1.0`.
- class-balanced fractional updates prefer the lowest `updatecnt` tier within
  each class quota.
- `input_oris_fixed` disables fractional updates and trailing reconstruction.
- `nsample`, `update_frac`, and trailing reconstruction behave consistently
  for full-update versus fractional-update runs.
- user `maxits` is treated as a stage cap, not a total wrapper cap.
- `prob_state` initialization remains capped at five iterations except in
  `input_oris_fixed`, where it is the only stage and uses the stage cap.
- `startit`, `which_iter`, and `extr_iter` remain monotonic global wrapper
  counters.
- stale `endit` is deleted before each stage-level base call, and the returned
  `endit` is used to advance the wrapper counter.
- `prob_neigh` sets `nspace_sub=500` and uses `prob_neigh_mode=sum`
  unless the caller overrides it.
- stage stop logic uses `mi_state` over the latest sampled subset when
  available.
- `prob_state` and `shc_smpl` use state-overlap target `0.95`.
- `prob_neigh` uses state-overlap target `0.99` by default.
- stage stop logic does not claim FSC or pose convergence.
- every active particle must have `updatecnt > 0` before final reconstruction.
- final reconstruction deletes crop and fractional-update controls.
- final reconstruction is at native project box and sampling.
- final reconstruction disables NU filtering and `nu_refine`.
- `nu_refine` remains forced to `no` and absent from the `refine3D_multi` UI.
- `combine_eo` remains absent from the `refine3D_multi` UI and
  `combine_eo=yes` is rejected.
- `lplim_crit=0.5` remains the multi-state default.
- wrapper code does not bypass base `refine3D` strategy ownership.
- NU handoff remains project/volume-assembly owned.
- per-state masks remain `volassemble` artifacts, not wrapper artifacts.
- `continue=yes` remains governed by base strategy support.

## 17. Refactor Rules

- Keep `refine3D_multi` as a wrapper around base `refine3D`.
- Do not move matcher behavior into `exec_refine3D_multi`.
- Do not move volume assembly or NU filtering into `exec_refine3D_multi`.
- Preserve all-or-none state-volume input.
- Preserve mode-specific project-label authority:
  `input_oris_start` reuses labels only when `nstates` is omitted,
  `input_oris_refine` reuses labels as the neighborhood-refinement start, and
  `input_oris_fixed` rejects labels greater than one.
- Preserve explicit state 0/1 initialization mode for projects without
  multi-state labels or for `input_oris_start` calls that supply `nstates`.
- Keep stage counters global and monotonic.
- Keep `maxits` semantics documented if changing them; total-wrapper and
  per-stage budgets are easy to confuse.
- Keep final reconstruction independent of fractional-update and crop state.
- Keep `nu_refine` out of `refine3D_multi`; ordinary NU filtering is only
  available through `filt_mode=nonuniform_lpset`, and NU shell expansion is not
  part of multi-state analysis.
- Preserve the active-particle update coverage guard.
- Treat state volumes, half maps, FSCs, NU products, automasks, partial
  reconstructions, and final maps as explicit workflow contracts.
