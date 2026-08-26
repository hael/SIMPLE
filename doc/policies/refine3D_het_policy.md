# Refine3D Het Policy

> **Development status:** `refine3D_het` is in development. Its public
> interface, stage plan, defaults, artifacts, and scientific policy are
> provisional and may change before the workflow is declared production-ready.

This document records the current implementation policy and the intended
boundaries for `refine3D_het`. It describes the frequency-marched heterogeneous
wrapper around base `refine3D`; low-level projection matching, probabilistic
search, partial reconstruction, volume assembly, automasking, and nonuniform
filtering remain governed by the related policy documents.

Related workflow policies:

- [refine3D_policy.md](refine3D_policy.md)
- [refine3D_multi_policy.md](refine3D_multi_policy.md)
- [refine3D_auto_policy.md](refine3D_auto_policy.md)
- [automasking_policy.md](automasking_policy.md)
- [nonuniform_filtering_policy.md](nonuniform_filtering_policy.md)
- [sigma_calculation_policy.md](sigma_calculation_policy.md)
- [particle_cache_policy.md](particle_cache_policy.md)

Primary implementation files:

- `src/main/ui/simple/simple_ui_refine3D.f90`
- `src/main/exec/simple_exec_refine3D.f90`
- `src/main/commanders/simple/simple_commanders_refine3D.f90`
- `src/main/params/simple_parameters_phases.f90`
- `src/main/strategies/parallelization/simple_refine3D_strategy.f90`
- `src/main/strategies/search/simple_strategy3D_matcher.f90`
- `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`
- `src/main/commanders/simple/simple_commanders_rec_distr.f90`

## 1. Status and Scope

`refine3D_het` is an experimental heterogeneous 3D refinement workflow. It
starts from a requested number of particle states, prepares a fractional or
full-update sampling policy, establishes state reference volumes, advances the
matching resolution in short fixed-length stages, ensures that every active
particle has been updated, and performs a final all-particle reconstruction at
native project sampling.

Unlike `refine3D_multi`, the current `refine3D_het` implementation does not have
public `multivol_mode`, `flex`, or fixed-orientation branches. It is a separate
wrapper with its own UI program and commander. It may reuse base `refine3D`
services and shared parameter definitions, but it must not derive its public UI
by copying the `refine3D_multi` UI object at runtime.

The provisional workflow contract is:

1. require `nstates >= 2`
2. reconcile the requested count with project state labels
3. initialize missing multi-state labels and virgin orientations
4. derive the per-iteration sample target and update fraction
5. prepare balanced class or cluster sampling when fractional updates use it
6. derive translation and crop settings when autoscaling is enabled
7. validate, recover, or reconstruct a complete set of state references
8. run frequency-marched base-`refine3D` stages
9. verify that every active particle has `updatecnt > 0`
10. assign any still-missing particles without reconstructing intermediate maps
11. reconstruct and postprocess final maps from all particles at native sampling
12. write final reconstruction products

`refine3D_het` is not a matcher, reconstructor, volume assembler, or
nonuniform-filter implementation. It is an orchestration layer. Once it calls
base `refine3D` or `reconstruct3D`, the called workflow owns the corresponding
particle-domain or volume-domain operation.

## 2. Entry Points and Ownership

The public command is `prg=refine3D_het`.

`simple_ui_refine3D.f90` owns the independent `refine3D_het` UI definition and
registers it under the display name `Heterogeneous 3D Refinement`. The UI record
requires an `sp_project` and is registered with `simple_exec`.

`simple_exec_refine3D.f90` owns dispatch:

- ordinary calls execute `commander_refine3D_het`
- calls with `nrestarts` use `restarted_exec` with the `refine3D_het` program
  name

`simple_commanders_refine3D.f90` owns wrapper policy:

- hard and overridable workflow defaults
- state-label validation and initialization
- sampling and automatic or explicit iteration-budget selection
- class/cluster sampling-sidecar preparation
- autoscaling and frequency-marching plans
- starting-volume selection
- initial mapping to supplied references when trailing reconstruction requires
  it
- update-coverage enforcement
- final `reconstruct3D` command preparation
- final output copying through `write_final_rec_outputs`

Base `refine3D` owns each refinement call after the wrapper has prepared its
command line. The wrapper must not duplicate or bypass:

- reprojection-model materialization
- probabilistic table generation
- hard-assignment search
- pose, state, correlation, and update-count writes
- partition-local reconstruction output
- trailing-reconstruction bookkeeping
- `volassemble`
- FSC and resolution registration
- per-state reference filtering and masking

`reconstruct3D` owns startup and final reconstruction passes. The wrapper owns
only their command-line policy and the choice of when they run.

## 3. Public Command and Defaults

The current UI exposes these options directly:

- image input: `vol1..volN`
- parameter and I/O: `cache`, `cache_dir`
- search: `maxits`, `nstates`, `nsample`, `prob_neigh_mode`, `inpl_cont`,
  `pgrp`, `ptcl_src`, `center`, `autoscale`
- filter: `filt_mode`, `lpstart`, `lpstop`, `ml_reg`, `envfsc`, `envmsklp`
- mask: `mskdiam`, `automsk`, `nu_msk_sig`
- compute: `nparts`, `nthr`

The commander forces these hard workflow values even if a caller supplies a
different value:

- `balance=yes` initially
- `greedy_sampling=no`
- `frac_best=1.0`
- `trail_rec=yes` initially
- `objfun=euclid`
- `lplim_crit=0.5`
- `incrreslim=no`
- `nu_refine=no`
- `combine_eo=no`
- `multivol_mode=independent`

The effective `balance` and trailing-reconstruction state can later be disabled
when the wrapper selects full-particle updates or cannot derive class sampling.

The commander supplies these values only when the user did not define them:

- `filt_mode=nonuniform_lpset`
- `envfsc=no`
- `mkdir=yes`
- `center=no`
- `sigma_est=global`
- `prob_inpl=yes`
- `refine=prob_neigh`
- `prob_neigh_mode=state`
- `autoscale=yes`
- `ml_reg=yes`
- `lpstart=10.0`
- `lpstop=6.0`
- `automsk=no`
- `nsample=10000` before state multiplication
- `keepvol=no`
- `nspace=5000`
- `nspace_sub=500`

The independent UI definition must remain aligned with the commander. In
particular, `lpstart` and `lpstop` are both public frequency-marching controls,
and the public filtering choices are currently only `nonuniform_lpset` and
`none`.

After the first `parameters%new` call has established the execution context and
private project, the wrapper sets `mkdir=no` for its child commands. Child
refinement and reconstruction calls must remain in that execution context.

## 4. Project and State Model

The parsed command line must provide an effective `nstates` of at least two.
The current wrapper does not infer a missing requested state count from a
multi-state project before parameter parsing.

Active-particle counts use `ptcl3D` rows with `state > 0`. A project with no
active particles is rejected.

The current state-label branches are:

- if the project reports one positive state label, the wrapper generates a
  uniform random labeling over `1..nstates`
- if the project reports exactly `nstates` positive labels, those labels are
  reused
- any other project label count is rejected as inconsistent with the command
  line

If the `ptcl3D` orientation field is virgin, orientations are randomized before
refinement. This applies both when the one-state input is expanded and when a
project already contains the requested number of state labels.

After initialization, every state in `1..nstates` must contain at least one
active particle. Empty states cause a hard failure and their populations are
reported before exit. `state=0` is not part of the active population.

The random state initialization is a starting condition, not evidence that the
states represent distinct structures. Scientific interpretation begins only
after state assignments and maps have been refined and independently assessed.

## 5. Particle Sampling and Update Fraction

When `nsample` is absent, the commander starts from `10000` particles per state.
The current implementation then computes:

```text
nsample_total = min(active_particles,
                    min(100000, nsample_input * nstates))
```

Consequently, an explicitly supplied positive `nsample` is also currently
treated as a per-state value before the `100000` and active-particle caps are
applied. The UI describes it as particles sampled per iteration per state.

The resulting total target must be at least one.

When the active-particle count does not exceed the target:

- the effective update fraction is `1.0`
- `update_frac` is removed from the parent command line
- fractional-update mode is disabled
- trailing reconstruction is disabled in the effective parameter state

When active particles exceed the target:

```text
update_frac = nsample_total / active_particles
```

- fractional updates are enabled only when the ratio is at most `0.99`
- trailing reconstruction follows `trail_rec`, which defaults to `yes`
- ratios greater than `0.99` are treated as full updates

The wrapper does not set `ufrac_trec`. Realized sampled/update-count metadata
and the base trailing-reconstruction contract remain authoritative.

Before staged refinement, previously sampled `ptcl3D` metadata is cleared by
removing `sampled` and `updatecnt` when the project reports prior sampling. The
new run then measures its own update coverage.

## 6. Class- and Cluster-Balanced Sampling

Fractional updates start with `balance=yes`, `greedy_sampling=no`, and
`frac_best=1.0`. The wrapper prepares the class-sampling sidecar consumed by the
base matcher when usable `ptcl2D` metadata exists.

For ordinary, non-partitioned input:

- selected 2D class indices are obtained from the project
- particle membership statistics are collected for those classes
- the resulting class samples are written to `CLASS_SAMPLING_FILE`

For `partition=yes` input:

- `cls2D` must contain `cluster` metadata
- clusters belonging to inactive class rows are excluded
- only populated cluster indices are retained
- particle membership statistics are collected using the `cluster` label
- the resulting samples are written to the same sidecar

If the `ptcl2D` field is virgin, the wrapper disables balancing. If the run is a
full update, it also disables balancing because no fractional subset has to be
distributed. Base `refine3D` then follows its unbalanced/full-update path.

The sidecar format and the matcher behavior that consumes it are shared
infrastructure. `refine3D_het` owns only the choice of class or cluster bins and
the timing of sidecar generation.

## 8. Autoscaling and Translation Limits

The master low-pass/crop record starts at native project box and sampling with:

- translation limit `min(8.0, max(2.0, AHELIX_WIDTH / smpd))`
- autoscaling disabled in the record
- native `box_crop`
- native `smpd_crop`

With `autoscale=no`, the wrapper removes `box_crop` and `smpd_crop` from the
child command line and retains native sampling.

With `autoscale=yes`, the wrapper calls the shared low-pass planning helper for
one stage using the parsed `lpstart` and `lpstop`. It then:

- sets `trs` from the returned translation limit
- sets `box_crop` and `smpd_crop` when the helper selects a crop
- removes crop overrides when the helper retains native sampling

All frequency-marched refinement blocks use the same master translation and
crop plan. The frequency limit changes per block; the crop does not.

The final reconstruction always removes crop overrides so its authoritative
maps are regenerated at native project sampling.

The same parsed `lpstart` therefore controls autoscale planning, the explicit
input-volume mapping pass, and the frequency-marching start.

## 9. Starting State Volumes

State-volume input is all-or-none.

If every `vol1..volN` is supplied:

- every file must exist
- every volume must be cubic with the native project particle box
- every volume must have sampling matching the native project particles within
  the implementation tolerance
- those volumes become the initial state references

If some, but not all, state-volume keys are supplied, the wrapper fails. Partial
multi-state input is ambiguous and must not be accepted.

If state labels were expanded from a single state, or an existing multi-state
project had virgin orientations, explicit input volumes are required. The
wrapper fails when no complete command-line volume set is present; it does not
use project volumes or reconstruct starting volumes in these branches.

Otherwise, if no command-line volumes are supplied, the wrapper first checks
project output volumes:

- `os_out` must contain a `vol` entry for every requested state
- every referenced file must exist
- every output volume must have the native project box
- sampling must be positive and match the native project sampling

If the complete project set is compatible, its paths are injected into
`vol1..volN` for the staged base calls.

If no compatible command-line or project set is available, the wrapper runs a
startup `reconstruct3D` pass. That pass uses:

- `objfun=cc`
- `postprocess=no`
- `nu_refine=no`
- `nsample` equal to the derived total sample target
- the requested `nstates`
- no `trail_rec`, `refine`, sigma-estimation, fractional-update, or stale
  `endit` controls

The resulting `vol_stateNN.mrc` files become `vol1..volN`. The startup
reconstruction inherits the current crop command line when autoscaling selected
a crop; unlike the final reconstruction, it does not explicitly delete
`box_crop` and `smpd_crop`.

The wrapper does not silently accept incompatible project volumes. It falls
back to reconstruction instead. Incompatible explicit input volumes are a hard
error because the user explicitly selected them.

## 10. Initial Mapping to Supplied Volumes

When a complete explicit input-volume set is supplied and the effective run
uses trailing reconstruction, the wrapper performs one initial greedy mapping
pass before the frequency march. This establishes assignments and reconstructed
state references compatible with the fractional/trailing update path.

The initial mapping pass:

- calls base `refine3D`
- samples at most `100000` active particles
- sets `update_frac` to that sample count divided by the active count
- sets `balance=no`
- sets `frac_best=1.0`
- sets `fillin=no`
- sets `trail_rec=yes`
- sets `volrec=yes`
- sets `filt_mode=none`
- runs exactly one iteration at iteration number `1`
- uses `refine=greedy` and `greedy_sampling=yes`
- uses `nspace=2500`
- uses `lp=lpstart`
- removes `partition` and stale `endit`

The mapping command currently inherits `automsk` from the parent. Therefore an
explicit `automsk=yes` is incompatible with its forced `filt_mode=none` and is
rejected by base refine3D parameter validation. The bootstrap command must also
disable `automsk`, or retain a compatible nonuniform filtering mode, before the
public HET automasking option works with this path.

After this pass, the ordinary `vol_stateNN.mrc` products replace the supplied
paths as the references for staged refinement.

If the effective run is a full update, trailing reconstruction is disabled and
the supplied volumes are used directly; this extra mapping pass does not run.

The initial pass is preparation and uses iteration number one. It also sets the
wrapper's global iteration counter to one, so the frequency march begins at
iteration two. When the initial mapping pass does not run, the global counter
remains zero and frequency marching begins at iteration one.

## 11. Search Modes

The default search mode is `refine=prob_neigh`. The commander currently accepts
these base search modes from a direct command line:

- `prob_neigh`
- `prob`
- `shc`

Other values are rejected by the HET wrapper even if base `refine3D` supports
them in another workflow.

For `prob_neigh`, supported neighborhood modes are:

- `state`, the HET default
- `geom`

`prob_neigh_mode` is validated and forwarded only for `prob_neigh`. The
`state` mode pools neighborhood policy by state; the underlying table and
candidate behavior remain owned by the base refinement strategy.

The wrapper defaults to `nspace=5000` and `nspace_sub=500`. It also defaults
`prob_inpl=yes`; the public UI independently exposes continuous in-plane
refinement through `inpl_cont=yes`.

Search-mode changes must preserve the base division between probabilistic
candidate preparation and hard assignment. The HET wrapper must not implement
candidate scoring directly.

## 12. Frequency-Marching Stage Plan

Staged refinement is divided into blocks of three iterations:

```text
STEP_IT = 3
nstages = ceiling(maxits / 3)
```

The final block contains the remainder when `maxits` is not divisible by three.
For example, with no initial mapping pass, `maxits=10` produces planned ranges
`1..3`, `4..6`, `7..9`, and `10`. With an initial mapping at iteration one,
the same frequency-march budget produces `2..4`, `5..7`, `8..10`, and `11`.
These ranges assume no block converges early.

The native Fourier-index limits are derived as:

```text
find_start = max(5, fourier_index(lpstart))
find_stop  = min(box / 2 - 2, fourier_index(lpstop))
```

For more than two stages, the index increment is divided by `nstages - 2`.
This plans the first stage at the starting limit and gives the final two stages
the stop limit after low-pass clamping. With two stages, the plan is start then
stop. With one stage, it remains at the start.

Each block receives:

- `maxits` equal to the number of iterations in that block
- `minits=1`
- `startit`, `which_iter`, and `extr_iter` equal to the wrapper's last completed
  global iteration plus one
- `trs` from the master autoscaling plan
- `lp` from the block's frequency plan, never numerically below `lpstop`
- the parent `refine`, `nspace`, `nspace_sub`, sampling, filter, and mask policy

Each block is an ordinary base-`refine3D` execution. The base commander writes
`endit`; the HET wrapper stores it as the global iteration counter, derives the
next block's start from it, and deletes it after each call. With `minits=1`, a
block may finish before its planned cap when base `refine3D` reports
convergence; `maxits` is therefore a planning horizon rather than a guarantee
that every planned iteration executes. HET forces
`multivol_mode=independent` and does not expose or supply an `overlap` control.

After all blocks, the wrapper deletes ordinary distance and assignment
temporary files for both partitioned and non-partitioned naming forms.

## 13. Delegation to Base Refine3D

For each frequency block, the wrapper changes `prg` to `refine3D` and invokes
`commander_refine3D`. Base strategy selection then follows the ordinary
refinement policy:

- `nparts` without `part` selects the distributed master
- otherwise the shared-memory strategy is selected
- distributed workers use the private worker route with their partition
  arguments

Base `refine3D` owns:

- reference section generation from current state volumes
- probabilistic pre-alignment where the selected mode requires it
- per-iteration particle sampling
- hard-assignment search and pose/state writes
- Euclidean sigma handling
- `sampled` and `updatecnt` bookkeeping
- partial even/odd reconstruction output
- trailing-update reconstruction handoff
- volume assembly and state-volume publication
- FSC and resolution metadata

The wrapper's forced `objfun=euclid`, `sigma_est=global`, `frac_best=1.0`, and
`lplim_crit=0.5` become inputs to those inherited contracts. Any change in
their low-level scientific meaning belongs in the base policies first.

## 14. Trailing Reconstruction and Update Coverage

Trailing reconstruction is enabled only for a genuinely fractional run and
defaults to `yes` there. It remains part of the base `refine3D` partial
reconstruction and `volassemble` contract. The HET wrapper must not synthesize
or merge partial Fourier volumes itself.

After frequency marching, the wrapper reads the persisted `ptcl3D` segment and
requires:

```text
state > 0  implies  updatecnt > 0
```

If `updatecnt` is absent, the workflow fails. If some active particles have not
been updated, the wrapper performs a dedicated missing-update assignment pass
and checks the project again.

The missing-update pass:

- calls base `refine3D`
- sets `balance=no`
- sets `frac_best=1.0`
- sets `fillin=no`
- sets `update_frac=1.0`
- sets `trail_rec=no`
- sets `volrec=no`
- runs one greedy iteration
- sets `greedy_sampling=yes`
- sets `update_missing=yes`
- removes `partition` and stale `endit`

This pass updates particle assignments and accounting only. It must not replace
the last staged volumes because volume reconstruction is explicitly disabled.

Final reconstruction is prohibited until a second coverage read confirms that
no active particle is missing an update.

## 15. Filtering, Masking, and Volume Assembly

The default staged-reference policy is:

- `filt_mode=nonuniform_lpset`
- `nu_refine=no`, forced and not exposed
- `ml_reg=yes`
- `automsk=no`
- `envfsc=no`
- `lplim_crit=0.5`

The other public HET filtering mode is `none`. `fsc`, `uniform`, and plain
`nonuniform` are not currently accepted by the wrapper.

`volassemble` remains the owner of volume-domain processing. Depending on the
selected mode and inherited base policy, it may restore half maps, calculate
FSCs, run ML regularization, run ordinary NU filtering, publish reference
products, and return an LP-set matching bandwidth through the project.

`nu_refine` shell expansion is not part of `refine3D_het`; the commander forces
it off. `combine_eo` is also forced off because this workflow is a multi-state
LP-set refinement rather than a combined even/odd terminal-alignment mode.

If `automsk=yes`, parameter validation requires a nonuniform filtering mode.
With the current HET choices, that means `filt_mode=nonuniform_lpset`.
State-specific NU envelope production and matching-reference flattening follow
[automasking_policy.md](automasking_policy.md). The wrapper must not introduce a
second mask artifact or mask-selection path.

With `envfsc=yes`, envelope-corrected FSC behavior follows the base volume
assembly policy. `envmsklp` controls that envelope scale and does not directly
replace the frequency-marching `lp` value.

## 16. Final Reconstruction and Outputs

After update coverage is complete, the wrapper runs `reconstruct3D` using all
particles and the final assignments. The command line:

- sets `prg=reconstruct3D`
- sets `outfile=RESOLUTION_FINAL.txt`
- sets `postprocess=yes`
- removes `trail_rec`
- removes `refine`
- removes objective-function denominator controls
- removes `sigma_est`
- removes `update_frac` and `ufrac_trec`
- removes stale `endit`
- removes `box_crop` and `smpd_crop`
- sets `objfun=cc`
- sets `nu_refine=no`

When staged refinement used nonuniform filtering, the final pass additionally
sets `filt_mode=none` and `automsk=no`. Ordinary NU filtering is therefore a
reference-generation feature during refinement, not a second final-map
postprocessing pass.

Deleting the crop overrides makes this final reconstruction authoritative at
native project box and sampling. After reconstruction, the wrapper reads the
project output segment and calls `write_final_rec_outputs` with `lpstop` to
publish the final per-state products and low-pass snapshots.

Final-map scientific interpretation must use these all-particle products, not
an arbitrary fractional or trailing intermediate map.

## 17. Restarts, Execution Mode, and Artifacts

The execution router supports `nrestarts` through the common restarted-command
mechanism. The HET UI does not currently expose a separate `continue` control,
and the wrapper has no independent resume implementation.

Every child base-`refine3D` call still inherits the selected shared-memory or
distributed strategy. Resume and worker validation must not be weakened by the
wrapper; any direct `continue` behavior is limited to what the selected base
strategy supports.

The staged workflow uses ordinary base-refinement artifacts, including:

- state volumes: `vol_stateNN.mrc`
- half maps: `vol_stateNN_even.mrc`, `vol_stateNN_odd.mrc`
- unfiltered half maps: `vol_stateNN_even_unfil.mrc`,
  `vol_stateNN_odd_unfil.mrc`
- FSC records: `fsc_stateNN.bin`
- per-iteration FSC plots and benchmark files
- partition-local even/odd partial reconstructions
- partition-local density weights
- reprojection models
- the class/cluster sampling sidecar for balanced fractional updates

The final reconstruction writes ordinary `reconstruct3D` products and
`RESOLUTION_FINAL.txt`; `write_final_rec_outputs` then writes the established
final state-map copies and low-pass snapshots.

The wrapper deletes distance and assignment temporary files after frequency
marching. It must not delete user-supplied reference volumes or compatible
project volumes selected at startup.
