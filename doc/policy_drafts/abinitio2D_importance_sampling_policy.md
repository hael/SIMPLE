# Proposed Abinitio2D Importance-Sampling and Fractional-Update Scheme

## Status

This document is a proposal, not a description of current behavior.

It defines a 2D scheme that mirrors the current 3D importance-sampling and fractional-update policy while respecting the main differences of `abinitio2D`:

- there is no reliable prior class information at stage 1
- 2D classification typically tolerates a larger fixed particle sample per iteration
- class-average restoration should use class-local fractional-update weights rather than one global carry-over weight

In this document, "stage" always means an `abinitio2D` workflow stage as controlled by `simple_abinitio2D_controller.f90`, not a single iteration and not a two-step mini-workflow.

## Design goals

The proposed 2D design should preserve the same high-level invariants as 3D:

- outer sampling chooses which particles participate in the current iteration
- probabilistic alignment reuses that exact sampled subset rather than resampling
- the restoration step consumes the realized update fraction rather than a nominal target
- stage control decides when the workflow changes search behavior without hardwiring the sampling scheme to any one refine mode

At the same time, it should differ from 3D where the problem structure requires it:

- stage 1 must start from random particle sampling because no class statistics are available yet
- the sampled subset size should be controlled by a fixed large target count instead of class-balance-driven sampling
- the default target count should still be overridable with `nsample` for testing
- fractional-update weighting should be computed per class, because classes in 2D play the same role that states play in the multi-state 3D workflow

## Proposed constant

Introduce:

```fortran
integer, parameter :: NPTCLS2SAMPLE_2D = 200000
```

This constant is the target outer sample size for `abinitio2D`.

It should be treated as the default target, not as a hardwired mandatory value.

## Proposed command-line override

For testing, benchmarking, and controlled comparisons, the workflow should allow:

```text
nsample=<integer>
```

to override the default `NPTCLS2SAMPLE_2D` target for the whole run.

The policy should be:

- if `nsample` is not provided, use `NPTCLS2SAMPLE_2D`
- if `nsample` is provided, use `nsample` as the run-local sample-count target
- the rest of the sampled-update and class-local restoration scheme remains unchanged

So the run-local target becomes:

```text
nsample_target_2D = nsample if provided, else NPTCLS2SAMPLE_2D
```

The effective global update fraction for an iteration becomes:

```text
update_frac_2D = min(1.0, real(nsample_target_2D) / real(nptcls_eff))
```

where `nptcls_eff` is the number of active particles with `state > 0`.

This mirrors the 3D policy in spirit:

- use a fixed target sample count
- convert it into a realized global update fraction at runtime
- keep the actual subset encoded through `sampled` and `updatecnt`

## Proposed public policy

The proposed `abinitio2D` workflow is:

1. choose a random outer particle subset of size `nsample_target_2D` or all active particles if fewer exist
2. record that subset in `ptcl2D` through `sampled` and `updatecnt`
3. apply a stage-controlled update policy across the full `abinitio2D` stage schedule, just as `abinitio3D` does across its stages
4. in stage 1, use the sampled subset to generate the first class assignments and class averages, but do not fractionally mix previous class-average sums
5. in every later stage, continue sampling a random subset of size `nsample_target_2D` each iteration
6. if probabilistic alignment is enabled in the later probabilistic stages, sample once in `prob_align2D` and reuse that subset in `prob_tab2D` and `cluster2D_exec`
7. restore class averages with class-local carry-over weights derived from the realized sampled fraction in each class

This keeps the 2D workflow conceptually aligned with 3D while avoiding fake early class-balance assumptions.

## Refine-mode independence

The sampling and fractional-update scheme should be orthogonal to the chosen 2D refine mode.

That means:

- outer sampled-update policy should not require `snhc_smpl`
- outer sampled-update policy should not require `prob`
- class-local fractional restoration should not require any one search strategy
- the same sampled-update contract should work with whatever 2D refine mode the stage controller selects

The only refine-mode-specific exception is probabilistic pre-alignment:

- when `l_prob_align_mode` is off, the current iteration samples directly in the matcher path
- when `l_prob_align_mode` is on, the current iteration must sample once in `prob_align2D` and reproduce that subset later

So the correct abstraction is:

- sampled outer update is a workflow policy
- refine mode is a search policy
- probabilistic sample-once-and-reuse is a special execution rule only for the probabilistic pre-alignment path

## Why stage 1 must differ from 3D

In 3D ab initio, the workflow can lean on prior `ptcl2D` clustering information and class-sampling statistics.

In 2D ab initio, stage 1 has no trustworthy class statistics yet. Therefore:

- class-balanced outer sampling is not available at stage 1
- class-dependent `frac_best` logic is not meaningful at stage 1
- the correct stage-1 policy is random outer subset selection with no fractional carry-over of previous class averages

This should be treated as a deliberate design choice, not as a missing optimization.

## Proposed stage policy

The intended 2D mirror of 3D is not a special "stage 1 versus stage 2" toy schedule. It is a full stage-controlled policy over the complete `abinitio2D` workflow.

Using the current 2D controller constants, the natural stage groups are:

- stage 1
- stages 2 to `PROBREFINE_STAGE - 1`
- stages `PROBREFINE_STAGE` to `NSTAGES_CLS`

with the EO terminal stage continuing the policy of its owning stage group.

These stage groups are sampling-policy groups, not hard locks on refine mode.

### Stage 1: random discovery stage

Policy:

- select a random subset of size `nsample_target_2D`
- set `sampled` and `updatecnt` for that subset
- run normal stage-1 2D search and clustering on only that subset
- do not fractionally update or trail the previous class averages
- write the resulting class averages as the new stage output with full weight

Rationale:

- there is no prior class information worth preserving
- this stage is analogous to initial 3D seeding, except the 2D version must discover class structure from scratch

Implementation consequence:

- outer sampling is active
- class-average carry-over is disabled
- `cavger_init_online` must behave as if `do_frac_update = .false.` for stage 1, even though only a subset of particles is processed

### Stages 2 to `PROBREFINE_STAGE - 1`: sampled update before the probabilistic stage group

Policy:

- keep using a random subset of size `nsample_target_2D`
- continue updating `sampled` and `updatecnt`
- enable fractional class-average carry-over
- compute restoration weights independently for each class
- allow any non-probabilistic refine mode selected by the stage controller

Rationale:

- these stages mirror the early sampled-update regime of 3D
- but the outer sampling remains random rather than class-balanced because the main goal is broad, fast 2D coverage
- the sampled-update policy should not force a particular 2D search engine

Implementation consequence:

- unlike the current code, the 2D stage controller should not simply delete `update_frac`
- instead it should set a stage-local `update_frac_2D` for every stage after stage 1

### Stages `PROBREFINE_STAGE` to `NSTAGES_CLS`: sampled update with optional probabilistic pre-alignment

Policy:

- keep the same outer sample size target
- if `l_prob_align_mode` is enabled, `prob_align2D` chooses the random outer subset once for the iteration
- if `l_prob_align_mode` is enabled, `prob_tab2D` reproduces that exact subset through `sample4update_reprod`
- if `l_prob_align_mode` is enabled, `cluster2D_exec` also reproduces that exact subset
- if `l_prob_align_mode` is disabled, the matcher path samples directly but still follows the same outer sampled-update policy
- class-average carry-over remains class-local

Rationale:

- this preserves the direct 2D analogue of the 3D sample-once-and-reuse policy when probabilistic alignment is active
- but it does not tie the whole workflow to `prob` as the only legal refine mode

### Final EO stage

Recommended initial policy:

- keep the same sampled outer-update contract
- continue class-local fractional carry-over
- let EO resolution estimation operate on the currently restored class averages

Optional later extension:

- add a final fill-in or full-update stage if benchmarking shows a consistent benefit

That extension should be treated as a separate policy change, not as part of the initial mirroring effort.

## Proposed stage-control table

This is the clearest direct mirror of the 3D workflow:

| `abinitio2D` stage group | Search-mode policy | Outer sampling policy | Fractional restoration policy |
| --- | --- | --- | --- |
| stage 1 | whichever stage-1 refine mode is selected | random subset of size `nsample_target_2D` | off |
| stages 2 to `PROBREFINE_STAGE - 1` | whichever non-probabilistic refine mode is selected | random subset of size `nsample_target_2D` | on, class-local |
| stages `PROBREFINE_STAGE` to `NSTAGES_CLS` | whichever refine mode is selected; if probabilistic pre-alignment is enabled then sample-once-and-reuse applies | random subset of size `nsample_target_2D` | on, class-local |
| EO terminal stage if enabled | current EO mode | same as owning stage group | on, class-local |

So the mirror should be read as:

- stage 1 is special only because restoration is non-fractional
- every later stage participates in the sampled-update policy
- probabilistic stages preserve sample-once-and-reuse semantics exactly as in 3D when that path is enabled
- refine mode remains a separate choice from the sampled-update scheme

## Outer subset-selection policy

### Random sampling helper

The current 2D code already has a random outer-sampling helper:

- `sample_ptcls4update2D`
- `sample4update_rnd`

That is the correct basis for the proposal.

The policy change is not to invent a new type of 2D sampler, but to make sampled update a first-class stage-controlled behavior in `abinitio2D`.

### Proposed effective update fraction

For every iteration in stages 2 and later:

```text
nsample_2D_eff   = min(nsample_target_2D, nptcls_eff)
update_frac_2D   = real(nsample_2D_eff) / real(nptcls_eff)
```

If `nptcls_eff <= nsample_target_2D`, the iteration naturally degenerates to a full update.

This makes the proposal robust for both small and large datasets without adding extra control branches.

For stage 1:

- the same random sampled subset is used
- but restoration still behaves as non-fractional because there is no meaningful previous class-average state to trail

## Proposed class-local fractional-update policy

### Problem with current behavior

Current 2D class restoration applies a single global carry-over weight:

```text
previous_weight = 1 - update_frac
```

through `apply_weights2cavgs(1.0 - p_ptr%update_frac)`.

That is too coarse for the proposed sampled workflow because different classes may receive very different numbers of sampled particles in a given iteration.

### Proposed class-local weight

For each class `k`, define:

```text
n_active(k)   = number of active particles currently assigned to class k
n_sampled(k)  = number of particles from the latest sampled round currently assigned to class k
rho(k)        = n_sampled(k) / n_active(k)
```

with the conventions:

- if `n_active(k) == 0`, keep the class empty
- if `n_active(k) > 0` and `n_sampled(k) == 0`, preserve the previous class-average sums for that class
- if `rho(k) >= 1`, fully replace the previous class-average sums for that class

Then restore class `k` using:

```text
previous_class_sum_weight(k) = 1 - rho(k)
current_class_sum_weight(k)  = rho(k)
```

This is the direct 2D analogue of the multi-state 3D idea:

- do not treat all classes as though they received the same fractional update
- let each class carry forward exactly the fraction implied by its realized sampled participation

### Why this mirrors 3D correctly

In multi-state 3D, the important policy idea is not literally "one global `update_frac` for everyone".

The important idea is:

- outer selection defines the updated subset
- downstream restoration respects how much of each independently restored object was updated

For 2D, the independently restored objects are class averages rather than volumes or states. Therefore class-local `rho(k)` is the right mirrored quantity.

## Required bookkeeping extensions

To support class-local restoration cleanly, the 2D workflow should expose a helper that computes class-local realized update fractions from `ptcl2D`.

Suggested helper semantics:

```text
get_class_update_fracs(self, ncls, rho)
```

where `rho(k)` is computed from:

- particles with active `state > 0`
- current class assignment `class == k`
- latest sampled round in `sampled`

This should live alongside the existing `simple_oris` sampling and getter logic rather than inside the classaverager.

That keeps the bookkeeping contract in one place.

## Ownership proposal

### `simple_commanders_abinitio2D`

Should own:

- `NPTCLS2SAMPLE_2D` policy
- optional `nsample` override for testing
- stage-specific activation of sampled update
- the rule that stage 1 is random sampled but non-fractional in restoration

### `simple_abinitio2D_controller`

Should own:

- conversion from `NPTCLS2SAMPLE_2D` to stage-local `update_frac`
- application of the optional `nsample` override before stage-local `update_frac` is computed
- stage-specific toggling of `update_frac`
- stage-local selection of search behavior without hardwiring the sampled-update scheme to a single refine mode

Current gap:

- `set_cline_cluster2D_stage` always deletes `update_frac`

Proposed change:

- thread stage-specific `update_frac` into the `cluster2D` command line for every stage after stage 1
- keep stage 1 restoration non-fractional even though stage-1 sampling is already active
- preserve refine-mode freedom and only activate sample-once-and-reuse when the chosen stage configuration actually enables probabilistic pre-alignment

### `simple_commanders_prob`

Should continue to own:

- sample-once behavior in `prob_align2D`
- writing the sampled subset before probability-table generation

This part of the current architecture is already close to the desired design.

### `simple_strategy2D_matcher`

Should continue to own:

- reuse of `sample4update_reprod` when `l_prob_align_mode` is enabled
- restore-time use of partial sums

Needed change:

- distinguish "sampled outer subset" from "fractionally mix previous class averages"
- stage 1 needs the first behavior without the second

### `simple_classaverager_restore`

Should own:

- execution of class-local carry-over weighting

Needed change:

- replace the current scalar `apply_weights2cavgs(w)` semantics with a class-local weighting path
- preserve the current scalar path as a fallback for non-sampled or legacy workflows if desired

## Minimal implementation sketch

### 1. Add a 2D sampling constant and runtime override

In the 2D control layer, add:

```fortran
integer, parameter :: NPTCLS2SAMPLE_2D = 200000
```

Then define the run-local target sample count as:

```text
nsample_target_2D = nsample if provided, else NPTCLS2SAMPLE_2D
```

### 2. Compute stage-local 2D update fraction

In `exec_abinitio2D` or `simple_abinitio2D_controller`, compute:

```fortran
update_frac_2D = min(1.0, real(nsample_target_2D) / real(nptcls_eff))
```

and thread it into the stage command line for every stage after stage 1.

### 3. Keep stage 1 sampled but non-fractional in restoration

Stage 1 should:

- sample randomly
- update `sampled` and `updatecnt`
- run clustering on the sampled subset
- disable carry-over weighting of previous class-average sums

### 4. Add a class-local update-fraction getter

In `simple_oris`, add a helper that returns `rho(k)` for all classes based on the latest sampled round.

### 5. Add class-local carry-over in class restoration

In `simple_classaverager_restore`, replace the scalar carry-over multiplier with a per-class multiplier:

```text
prev_class_weight(k) = 1 - rho(k)
```

applied independently to:

- even class sums
- odd class sums
- corresponding `ctfsq` accumulators

### 6. Keep probabilistic sample-once-and-reuse behavior

No policy change is needed here beyond ensuring the staged sampled-update logic is active when `prob_align2D` runs.

### 7. Keep the scheme mode-agnostic

Do not encode the scheme as:

- "sampled update means `snhc_smpl`"
- "fractional restoration means `prob`"

Instead encode it as:

- sampled update is enabled by stage policy
- class-local restoration is enabled by stage policy
- refine mode is chosen independently
- probabilistic sample-once-and-reuse is activated only when the selected mode uses `prob_align2D`

## Invariants to preserve

- The outer sampled subset must be chosen before `prob_tab2D` work begins.
- `prob_align2D`, `prob_tab2D`, and `cluster2D_exec` must operate on the same current sampled round.
- Stage 1 must not pretend to have reliable class-balance information.
- Fractional carry-over of class averages must be class-local, not only global.
- `sampled` must remain the current-round marker.
- `updatecnt` must remain the cumulative update history.
- The sampled-update scheme must remain usable with different 2D refine modes.

## Recommended first implementation cut

The first implementation should do only the following:

1. add `NPTCLS2SAMPLE_2D = 200000`
2. allow `nsample` to override the default target for testing
3. activate random outer sampled update for all `abinitio2D` stages
4. keep stage 1 restoration non-fractional
5. activate class-local fractional restoration in every later stage
6. compute class-local carry-over weights from the realized sampled fraction in each class
7. preserve sample-once-and-reuse behavior in the probabilistic stages
8. keep refine-mode selection orthogonal to the sampled-update scheme
9. leave fill-in and class-biased outer sampling for later

That gives a clean 2D mirror of the 3D design without overfitting early to heuristics that require mature class statistics.

## Summary

The proposed mirrored 2D policy is:

- fixed random outer sample size of 200,000 particles
- `nsample` command-line override for testing and benchmarking
- stage-controlled sampled update across the full `abinitio2D` workflow
- stage 1 sampled but with no fractional carry-over of class averages
- every later stage sampled with random outer selection and class-local fractional carry-over
- probabilistic stages reuse the exact sampled subset rather than resampling when that path is enabled
- the scheme does not lock `abinitio2D` into any one refine mode
- class averages use per-class realized update fractions in the same spirit that multi-state 3D restoration respects independently updated objects

This is the simplest 2D design that mirrors the 3D architecture without importing assumptions that only become valid once stable class structure already exists.
