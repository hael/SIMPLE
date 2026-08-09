# Refine3D Policy

This document records the current policy for the base `refine3D` command. It
does not describe the staged `abinitio3D` controller, the class-average
`abinitio3D_cavgs` initializer, or the automated wrapper `refine3D_auto`.

Related workflow policies:

- [abinitio3D_policy.md](abinitio3D_policy.md)
- [abinitio3D_cavgs_policy.md](abinitio3D_cavgs_policy.md)
- [abinitio3D_cavgs_reject_policy.md](abinitio3D_cavgs_reject_policy.md)
- [refine3D_auto_policy.md](refine3D_auto_policy.md)
- [refine3D_multi_policy.md](refine3D_multi_policy.md)
- [automasking_policy.md](automasking_policy.md)
- [nonuniform_filtering_policy.md](nonuniform_filtering_policy.md)
- [sigma_calculation_policy.md](sigma_calculation_policy.md)
- [particle_cache_policy.md](particle_cache_policy.md)
- [separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md](separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md)

## 1. Scope

`refine3D` is an iterative projection-matching workflow over an existing
project and one or more starting 3D references. Its durable contract is:

1. initialize project state and execution mode
2. materialize the iteration reprojection model from current Cartesian volumes
3. optionally run probabilistic pre-alignment
4. run particle-domain search and hard assignment
5. write partition-local Cartesian reconstruction inputs
6. run volume assembly when volume reconstruction is enabled
7. persist project state for the next iteration or caller

The supported reconstruction path is Cartesian. Cartesian matching still uses
polar Fourier central sections internally, but those sections are generated
from the current Cartesian volumes. `refine3D` does not own a separate
non-Cartesian reconstruction branch.

Particle-domain work stays in the probabilistic pre-step, matcher preparation,
search strategies, pose/state/shift updates, sigma updates during Euclidean
matching, and partial reconstruction writing. Volume-domain work stays in
`volassemble` and the volume postprocessing helpers it calls.

Shared-memory and distributed `refine3D` must preserve the same scientific
workflow and artifact contracts. Only toolbox lifetime, process launch, and
scheduling differ.

## 2. Entry Points and Execution Mode

The public command is `refine3D`, registered in `simple_ui_refine3D.f90` and
routed by `simple_exec_refine3D.f90`.

`simple_commanders_refine3D.f90` performs the thin command-level setup:

- validate multi-volume input against `nstates`
- set local defaults such as `mkdir`, `cenlp`, `oritype`, and `prg`
- select a strategy through `create_refine3D_strategy`
- run the iteration loop
- delegate finalization and cleanup to the selected strategy

Strategy selection is command-line shaped:

- `nparts` without `part` selects the distributed master strategy
- otherwise the shared-memory strategy is used
- distributed workers run the partition command path with `part` and `outfile`

`maxits` is the number of iterations to run in the current invocation.
`which_iter` starts at `startit`, and `extr_iter` follows the same per-call
counter unless supplied by a caller.

Shared-memory `refine3D` currently rejects `continue=yes`. Distributed
`continue=yes` resumes from project-carried output volumes, FSC files, and
run-local artifacts that match the requested partitioning and sigma mode.

## 3. Ownership

`simple_commanders_refine3D.f90` owns the base command entry point and the
shared iteration loop. It should stay thin.

`simple_refine3D_strategy.f90` owns:

- shared-memory versus distributed orchestration
- scheduler interaction
- iteration counters and convergence checks
- probabilistic pre-step dispatch
- matcher dispatch
- volume assembly dispatch
- per-iteration reprojection-model materialization
- run finalization

It must not absorb numerical volume postprocessing or assembled-reference
construction.

`simple_strategy3D_matcher.f90` owns `refine3D_exec`: consuming the
driver-generated reprojection model, selecting and running the search strategy,
updating orientations/states/shifts, updating Euclidean sigma estimates during
search, and writing partition-local Cartesian partial reconstructions.

`simple_commanders_prob.f90` plus `simple_eul_prob_tab*.f90` own
probabilistic pre-alignment and assignment artifacts consumed by the matcher.

`simple_commanders_rec_distr.f90` owns `commander_volassemble`: reducing
partials, restoring even/odd volumes, calculating FSCs, running volume
postprocessing, writing state volumes, and updating resolution metadata.

`simple_vol_pproc_policy.f90` owns NU-envelope cadence, artifact compatibility,
and legacy-density-mask compatibility. The NU setup API owns its spherical
`mskdiam` support invariant.

## 4. Project and Builder State

The project is durable workflow state. The `builder` owns derived execution
state for one parsed command line and must not be treated as durable policy
state.

Shared-memory `refine3D` rebuilds its strategy toolbox each iteration after
the iteration command line has settled. Persist any initial project mutations,
such as random orientations, sampling-counter cleanup, or state initialization,
before the next rebuild reads the project.

Distributed `refine3D` keeps prototype command lines for worker, probability,
sigma, assembly, and postprocess steps. These command lines are execution
templates, not independent workflow policies.

Fresh-start checks that decide whether to consume project-carried matching
metadata must use the stage interval: `which_iter <= startit` with
`continue != yes` is fresh, even when `startit` is greater than one.

Before a matcher iteration writes partition-local Cartesian partials, stale
partial reconstruction files for the active states and partitions must be
removed.

## 5. Startup and State Initialization

If multiple starting volumes are supplied, `nstates` must be defined and
greater than one. Multi-state startup requires either all `vol1..volN` inputs
or none.

When orientations are absent, `refine3D` randomizes 3D orientations. When
projection indices are absent, distributed initialization sets them from the
Euler sampling.

For multi-state distributed startup:

- existing state assignments must match the requested `nstates`
- fresh state labels are randomized uniformly, whether complete starting
  volumes are supplied or the project references are used
- probabilistic multi-state startup requires each state to exceed the minimum
  population guard

Even/odd partitioning is required for consistent half-map reconstruction. If it
is absent, initialization partitions the active orientation segment and writes
that state back to the project.

## 6. Reference Preparation

The strategy materializes the iteration reprojection model before
probabilistic pre-alignment or matcher work starts. Reference preparation:

1. reads the current Cartesian reference volumes
2. applies reference masking/filtering/centering policy
3. determines the active high-pass/low-pass shell range
4. writes even/odd PFTC reference sections for the matcher

When `automsk != no`, reference preparation replaces the ordinary spherical
reference mask with the current compatible `nu_envmask3D_stateNN.mrc`. The
selected regular or NU-derived reference is background-zeroed using the mask
transition and multiplied by the mask before projection. A compatible legacy
density mask is the transition fallback, followed by the sphere. This does not
modify particles, reconstruction volumes on disk, FSC curves, NU filtering, or
matching-bandwidth selection, and there is no separate `envref` control.

This path is valid only with `filt_mode=nonuniform|nonuniform_lpset`.
`automsk=yes|tight` with `none`, `uniform`, or `fsc` is rejected during
parameter validation.

The lifecycle is lag-by-one: volume assembly from iteration N writes the NU
envelope used when iteration N+1 materializes its reprojection model. A missing
mask at the initial bootstrap is therefore expected.

`volassemble` must not generate or promote matcher reprojections. The next
iteration's reprojection model is generated by the strategy from the current
Cartesian outputs.

The model file header is the authority for the matching shell range consumed
by matcher workers and probabilistic table workers.

Reference preparation must project only the active matching shell range. It
must not silently project to the crop interpolation limit unless matching
policy explicitly asks for that range.

The 3D low-pass range comes from one of these sources:

- explicit LP-set `lp`
- project-carried NU-selected handoff in a continuing nonuniform refinement
- current FSC files
- an existing project `lp` field

`lpstop` caps the selected matching bandwidth. The crop Nyquist limit caps the
final Fourier index after the matching policy has selected a candidate range.

## 7. Matching Topology

Gold-standard matching keeps even and odd references independent.

If low-resolution even/odd docking is needed for registration, it is applied
only after reference read/mask/filter and immediately before PFTC generation.
Those blended references are not written as half-map outputs and do not feed
FSC, automasking, NU filtering, or ordinary half-map handoff.

LP-set matching uses a merged registration reference. State count alone must
not force merged-reference matching; topology is controlled by LP-set policy.
`filt_mode=uniform` is disabled for multi-state search until it has a
state-specific low-pass contract.

In nonuniform mode, reference loading follows the NU policy:

- plain `nonuniform` prefers independent `_nu_filt` even/odd references and
  falls back to regular even/odd references before using a merged map
- `nonuniform_lpset` with active LP-set matching uses the merged reference and
  prefers the merged `_nu_filt` product when present

The ordinary reference low-pass filter is not applied on top of a NU reference
path. NU and ML-regularized references are treated as already filtered during
volume assembly.

## 8. Probabilistic and Matcher Work

The current workflow is probabilistic pre-alignment followed by hard particle
assignment. It is not a monolithic soft-assignment volume-integrated EM
implementation.

Use the terms "probabilistic pre-alignment" and "hard-assignment particle
update" unless the ownership, artifacts, and update model actually change.

`refine=prob*` modes run the probability pre-step before the main matcher.
`prob_neigh` dispatches to the neighborhood probability command. The matcher
then consumes the generated assignment artifact.

`prob_neigh_mode` controls how `prob_neigh` chooses candidate neighborhoods
before evaluating candidates:

- `state`: score coarse subspace representatives independently per state, pool
  the selected neighborhoods across states, and evaluate the same pooled
  projection search space for every active state.
- `geom`: select the subspace containing the particle's current projection for
  every active state, with no coarse scoring or pooled peaks.
- `shc`: use the direct stochastic candidate path. When shifts are enabled, it
  first estimates a shift seed before candidate scoring.
- `snhc`: use the same direct stochastic candidate path without the initial
  shift seed. It is consequently also the zero-shift default for this mode.

Probabilistic-table assignments use the calibrated likelihood path. The stored
distances are noise-normalized negative log-likelihoods for the Euclidean
objective, and evaluated candidates are sampled with weights proportional to
`exp(-dist)` over an explicit top-K support. The implementation uses a
per-particle minimum shift before exponentiation for numerical stability; this
does not change normalized weights.

The top-K truncation is deliberate. It defines the local discrete support
actually evaluated by the pre-alignment step; it is not meant to represent a
full posterior over all SO(3) grid points. If the CC objective is enabled, its
existing likelihood-compatible transformation remains the source of the
probability weights, but should not be described as a calibrated Gaussian
likelihood.

Likelihood-weighted probability-table modes may still profile or MAP-refine
shifts, and sometimes in-plane rotation, after stochastic candidate selection.
The stored assignment distance is then the refined/profiled objective value.
This is intended current behavior, not a full soft-assignment EM update.

### Continuous in-plane policy

`inpl_cont` has exactly two values. `no` is the default and preserves the
historical alternating search: continuous shift optimization with the
discrete in-plane angle callback. `yes` replaces every use of that callback
with the joint raw-Euclidean `(sx,sy,rotind_frac)` optimizer. The numerical
route is selected before shift-search objects are constructed; objective type
is a capability check and must never activate continuous-angle behavior by
itself.

The refine3D opt-in requires `objfun=euclid`, `objfun_den=no`,
`ptcl_src=raw`, and `projrec=no`. It is not restricted to `refine=shc`:
deterministic, neighborhood, evaluation, and probabilistic matcher routes use
the same policy wherever they perform pose search. A mode with no pose search,
such as sigma-only setup, naturally invokes neither optimizer. Unsupported
capability combinations fail validation rather than silently reverting to the
callback.

Joint candidate profiling fixes the already selected state and projection. It
evaluates every discrete in-plane cell once at the caller-supplied native
`(x,y)` shift and chooses the best index, exactly matching one legacy callback
selection. The continuous three-parameter solve starts from that index and is
bounded to plus or minus two cells around it. The `inpl_cont=yes` initializer
never searches other shifts; in particular, it does not use the legacy 5-by-5
shift/all-angle coarse initializer.

For probabilistic refinement, the probability-table schema remains discrete.
Joint optimization may evaluate a fractional in-plane coordinate while
profiling a candidate, but the assignment artifact carries only the canonical
rounded `inpl`, its score, and its shift. That shift is expressed in the
rounded-index frame so the matcher can recover the exact native shift seed
without storing a fractional coordinate in the table.

After hard state/projection/in-plane assignment, the matcher reruns the joint
three-parameter optimizer locally for the selected pose. The rounded assignment
index is authoritative: it converts the stored shift back to the native
particle frame and seeds a solve bounded to plus or minus two in-plane cells,
without another global all-angle selection. An accepted result persists one
coupled pose: fractional Euler `e3`, nearest integer `inpl`, shift, and score. A
finite non-improving evaluation retains the incoming discrete pose with its
re-scored objective value. A numerically invalid evaluation leaves the incoming
assignment untouched. No `inpl_cont=yes` path falls back to the legacy callback.
State and projection identity remain fixed during the final solve. A persisted
fractional `e3` is rounded to the canonical integer cell before the next search.

Joint acceptance is guarded twice, in the shared optimizer. An improvement
must be material -- it must beat the seed cost by a relative tolerance far
above solver roundoff -- otherwise the solve counts as non-improving and the
seed pose stands; this keeps the improved fraction an honest diagnostic and
prevents noise-level "improvements" from committing fractional poses. A
solution pinned to any joint search bound (the shift box or the plus-or-minus
two-cell angular window) is likewise demoted to non-improving: an
angular-bound-pinned solution contradicts the exhaustive seed selection, and a
shift-corner solution is the signature of descent into truncated-series
artifacts rather than signal.

Convergence `dist_inpl` reporting is symmetry-aware. When multiple symmetry
operations have the same Euler distance within numerical tolerance, including
projection-preserving operations in dihedral groups, the reported in-plane
distance uses the equivalent operation with the smallest in-plane change. A
symmetry-equivalent branch switch must not appear as a spurious 180-degree
in-plane move.

Shared-memory and distributed probabilistic child commands must retain
`inpl_cont`; reconstruction, sigma, assembly, and postprocessing children must
not receive matcher-only search policy.

For multi-state alignment (`nstates > 1`), shift-first candidate scoring is
disabled. The matcher and probability-table paths may still refine shifts after
candidate/state selection, but they must not use a shift seed estimated from a
particle's previous state to score candidates in other states.

The matcher must preserve a single particle-stack read per batch *within each
phase*: batch construction reads each particle once for search, and the
reconstruction phase reads each selected particle once after the search-phase
teardown (the PFTC memory phase boundary makes retaining the raw images across
phases prohibitive). Do not add further reads to either phase unless the
performance contract is explicitly changed.

When the downscaled particle cache is active (`cache=yes`), both phases read
the cache instead of the original full-size stacks; the single-read-per-phase
rule then applies to cache records. See
`doc/policies/particle_cache_policy.md` for the cache contract, including why
cached reconstruction is a deliberate, uniform-across-ranks numerics change.

## 9. Volume Assembly

When `volrec=yes`, matcher workers write partition-local Cartesian partials and
the strategy dispatches `volassemble`.

`volassemble`:

- reduces partition-local partial reconstructions
- restores dense even/odd half-volumes
- calculates FSC curves and state resolutions
- applies conical FSC curves for directional ML regularization when
  `ml_reg=yes` and `conical_fsc=yes`; this is opt-in
- restores merged state volumes
- derives the NU-evidence envelope from the completed accepted NU bank before
  unary storage is released
- writes derived reference products
- records resolution and NU matching metadata in the project

Volume assembly does not refresh matcher PFTC references. It only produces
Cartesian volumes and metadata for the next iteration.

The shared helper `restore_state_from_parts` owns the restoration sequence.
Changes to ML ordering, sampling-density correction, FSC handling, trailing
reconstruction, low-resolution even/odd insertion, or LP-set behavior must go
through that shared helper rather than duplicated assembly branches.

FSCs and FSC-derived diagnostics are calculated from dense Cartesian
half-volumes, not from sparse, intermediate, or low-resolution-blended
registration-reference representations.

## 10. Trailing and Combined Even/Odd

Trailing reconstruction blends even with previous even and odd with previous
odd. It runs before NU-evidence envelope generation and derived reference
filtering.

In LP-set mode, if trailing is active, assembly trails the even and odd
half-volumes first and then merges the trailed halves. Assembly must not use
low-resolution even/odd docking insertion for LP-set outputs.

When `combine_eo=yes`, distributed `refine3D` schedules one additional final
iteration after convergence or run-length termination. That final iteration:

- disables fractional update
- forces full update
- sets `combine_eo=yes`
- tightens the low-pass criterion to at most 0.143

The combined even/odd iteration is part of base `refine3D`, not a terminal
`refine3D_auto` or ab initio reconstruction step.

## 11. Finalization and Artifacts

On each iteration, strategy benchmark files should stay simple: context plus
one `TIMINGS (s)` section. Labels should include the coarse operation buckets
setup, probabilistic pre-step, matcher/scheduler, assembly/postprocess, and
total time, together with separately timed reprojection-model materialization
and group-sigma consolidation. These two measurements keep a per-iteration
sigma update distinguishable from reference preparation.

Each refine3D stage also writes one stage-entry benchmark at its first
iteration. It records the total stage initialization wall time and the nested
`calc_pspec` wall time. The latter is the per-particle sigma estimation that
may be reused across compatible stage changes; it must not be inferred from a
per-iteration setup bucket.

Distributed matching writes partition alignment documents and merges them into
the project after worker completion. Shared-memory matching writes the project
directly after each iteration.

On finalization:

- `endit` is written to the command line
- `startit` is removed from the command line
- active state volumes and FSCs are registered in `os_out`
- empty states are removed from `os_out`
- `cls3D` distributed runs map class-orientation output back to particles
- `JOB_FINISHED` is touched by the shared-memory path

Grouped sigma files are run-local noise-model state. They may be written and
consumed inside a running refinement, but they are not registered as `os_out`
handoff artifacts for a later `refine3D` execution.

## 12. Refactor Rules

- Preserve the particle-domain versus volume-domain boundary.
- Keep shared-memory and distributed workflows behaviorally aligned.
- Keep volume assembly as the owner of assembled-reference work.
- Do not move assembled-volume postprocessing into `refine3D_strategy`.
- Do not merge probabilistic particle-update logic with volume postprocessing.
- Do not introduce a second ambiguous source of matching references.
- Preserve the single image-stack read per matcher batch.
- Treat assignment maps, partition-local partials, state volumes, even/odd
  volumes, FSC files, automasks, and NU outputs as explicit workflow
  contracts.
- Preserve the two-mode in-plane boundary: only `new_legacy` may attach the
  callback, and no `inpl_cont=yes` path may invoke it.
- Initialize candidate-profiling joint solves with one callback-equivalent
  all-angle selection at the supplied shift; do not run the 5-by-5 coarse shift
  initializer.
- Treat the final 3D hard assignment as authoritative: refine locally around
  its rounded in-plane cell and never perform a second global angle selection.
- Keep probabilistic in-plane artifacts rounded and run durable fractional
  refinement only after final hard assignment.
- Keep `inpl_cont` on shared-memory and distributed matcher/probability child
  commands while stripping it from non-matcher children.
