# Iterative 3D Refinement

## Problem

Given particle images and one or more starting volumes, update particle poses
and state assignments and reconstruct the next references while preserving
even/odd half-map topology, objective/noise semantics, and a clear
particle-domain/volume-domain boundary.

## One iteration

Base `refine3D` performs:

1. materialize the reprojection model from current Cartesian volumes;
2. optionally run probabilistic pre-alignment;
3. search and commit one hard pose/state assignment per active particle;
4. update Euclidean noise estimates when required;
5. write partition-local Cartesian reconstruction statistics;
6. reduce and restore even/odd volumes in `volassemble`;
7. calculate FSC/cFAR and resolution metadata;
8. create filtered/masked references and project state for the next iteration.

Shared-memory and distributed implementations differ in launch and reduction
lifetime, not in this scientific sequence.

## Reference preparation

The strategy reads current Cartesian state volumes, applies the selected
reference mask/filter/centering policy, chooses the active Fourier shell range,
and writes polar Fourier central sections. Gold-standard topology keeps even
and odd references separate. Merged-reference matching is activated by LP-set
policy, not merely by the number of states.

With NU automasking, iteration `n` uses the compatible state envelope generated
by volume assembly at iteration `n-1`. The selected reference is
background-zeroed and multiplied by the mask immediately before reprojection;
the on-disk reconstruction and FSC inputs are unchanged.

## Candidate preparation and hard update

`refine=prob*` modes generate a bounded probability table before the matcher.
Neighborhood policy is one of:

- `state`: score coarse neighborhoods per state, pool them, then evaluate the
  same pooled support for every active state;
- `geom`: use the neighborhood containing the previous projection for each
  state;
- `shc`: stochastic direct candidates with a shift seed;
- `snhc`: stochastic direct candidates without that seed.

Euclidean table costs are noise-normalized negative log weights and are sampled
proportional to `exp(-dist)` over explicit top-K support. The matcher consumes
the selected candidate and commits one state, projection, in-plane cell, and
shift. Multi-state search does not use a shift seed from the old state to rank a
different state.

## Continuous in-plane precision

`inpl_cont=yes` is polish-only. Candidate scoring, probability-table
construction, and selection remain discrete. After the hard assignment, one
bounded joint `(sx,sy,theta)` solve may refine the committed pose while holding
state and projection fixed. Material improvement and non-bound-pinning guards
decide whether the fractional angle is persisted. Convergence statistics remain
discrete so polishing cannot change the exploration schedule.

## Partial reconstruction

Selected particles are reread once in the reconstruction phase and contribute
to per-state/per-half Cartesian statistics under their committed poses, CTFs,
and data weights. The search phase and reconstruction phase are separate memory
lifetime regions. A downscaled cache, when active, is used consistently by both.

## Volume assembly

`volassemble` reduces partials, restores dense halves through the selected
backend, calculates FSCs from dense Cartesian halfmaps, restores merged maps,
and runs ML regularization, nonuniform filtering, masking, and reference-product
generation in the prescribed order. It does not build the next iteration's
PFTC model; the strategy does that from the published Cartesian products.

The backend seam is the dense even/odd half pair. Downstream FSC, merged maps,
NU filtering, masks, filenames, and project updates are common to gridding and
PCG reconstruction.

## Completion

Convergence and run-length policy are owned by the caller/strategy. A
`combine_eo=yes` base run schedules one final full-update iteration with a
tighter matching limit; multi-state wrappers intentionally disable that mode.
Finalization registers active state volumes and FSCs and persists project
orientation/state metadata.

## Implementation

- Commander: `src/main/commanders/simple/simple_commanders_refine3D.f90`.
- Iteration strategy: `src/main/strategies/parallelization/simple_refine3D_strategy.f90`.
- Matcher: `src/main/strategies/search/simple_strategy3D_matcher.f90`.
- Volume assembly: `src/main/commanders/simple/simple_commanders_rec_distr.f90`.
- Policy: `doc/policies/refine3D_policy.md`.

