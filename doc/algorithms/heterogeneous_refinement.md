# Multi-State and Heterogeneous 3D Refinement

## Problem

Estimate several 3D states and particle-state assignments while preserving
pose coverage, preventing large views/classes from monopolizing fractional
updates, and ensuring that final scientific maps use every active particle.

SIMPLE exposes two production wrappers over base `refine3D`; neither implements
a separate matcher or reconstructor.

## Prior-orientation multi-state refinement

`refine3D_multi` starts from existing 3D orientations. It obtains `nstates`
from compatible project labels or the command line, validates an all-or-none
set of state volumes, and supports:

- `input_oris_refine`: refine state and pose; state-0/1 input first runs
  `prob_state`, then `prob_neigh`, while an existing multi-state project starts
  directly in `prob_neigh`;
- `input_oris_fixed`: hold Euler angles, projection direction, in-plane angle,
  and shifts fixed while `prob_state` updates state/correlation/accounting.

The default outer target is

```text
N_sample = min(100000, 10000*nstates).
```

Fractional subsets are balanced across current projection directions and favor
the lowest `updatecnt` tier within each direction. The default neighborhood is
geometric around each particle's prior projection; pooled per-state
neighborhoods are optional.

## Competitive heterogeneous refinement

`refine3D_het` can begin without trusted state labels. One-state input is
uniformly randomized across the requested states; an already labeled project
must contain exactly that state count. All states must remain populated.

The public sample target is interpreted per state, then capped at 100000 and
the active population. Fractional updates are class- or cluster-balanced when
usable `ptcl2D` metadata exists. The default `prob_neigh_mode=state` scores
coarse subspaces per state, pools selected neighborhoods, and evaluates the
same pooled support for each state so state comparison uses a common candidate
geometry.

Frequency marching is divided into blocks of at most three iterations. The
first block uses `lpstart`, later blocks advance in Fourier index toward
`lpstop`, and the final two planned blocks share the stop limit. Crop and
translation policy stay fixed across the blocks; only the active frequency
limit changes.

## State volumes and initialization

Explicit `vol1..volN` input is all-or-none and must match native project box and
sampling. Otherwise compatible project volumes are used or, where the wrapper
allows it, a startup reconstruction builds state references. Fractional HET
runs with explicit volumes perform one coarse greedy mapping pass to establish
assignments and trailing-compatible reconstructed state references.

`refine3D_multi flex=yes` is an additional state-0/1 initialization route. It
runs projection-aware `flex_pca`, adopts its hard labels and state volumes, and
then enters the ordinary multi-state stage plan.

## Inherited iteration

Every stage/block is a normal base `refine3D` call:

```text
state/pose candidate table -> hard particle update -> partial halfmaps
-> volume assembly -> filtered next references.
```

The wrappers own `nstates`, starting references, frequency schedules, outer
sampling, balancing, and stage termination. Base refinement owns candidate
scores, state/pose writes, sigma updates, partial reconstruction, FSC, trailing,
NU filtering, and masks.

## Filtering and final maps

Both wrappers default to ML-regularized `nonuniform_lpset` references with
static NU filtering (`nu_refine=no`). Automasking is optional and, when enabled,
uses state-specific NU-evidence envelopes through the base lag-by-one contract.
Combined even/odd terminal alignment is disabled.

After staged work, a missing-update pass assigns any active particle with
`updatecnt=0` without replacing the last staged volumes. Only after coverage is
complete does a fresh all-particle `reconstruct3D` run at native sampling.
Those all-particle per-state products—not an arbitrary fractional/trailing
intermediate—are the authoritative maps for interpretation.

## Implementation

- Wrappers: `src/main/commanders/simple/simple_commanders_refine3D.f90`.
- Base estimator: `src/main/strategies/parallelization/simple_refine3D_strategy.f90`
  and `src/main/strategies/search/simple_strategy3D_matcher.f90`.
- Policies: `doc/policies/refine3D_multi_policy.md` and
  `doc/policies/refine3D_het_policy.md`.

