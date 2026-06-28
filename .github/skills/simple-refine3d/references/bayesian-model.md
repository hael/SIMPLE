# Probabilistic 3D Refinement Model

## Scope

Use this note when the task is about the scientific or algorithmic meaning of
SIMPLE's probabilistic 3D refinement implementation.

## Core Model

SIMPLE's 3D refinement is a fixed-point workflow that alternates between:

- updating particle orientations, in-plane angles, shifts, and sometimes states
- writing particle-derived Cartesian partial reconstructions or polar partial sums from those assignments
- assembling and postprocessing Cartesian state volumes or state-major polar references through explicit assembly pathways

The implementation uses probability-like candidate generation, thresholded
sampling, regularization, even/odd cross-validation, and likelihood-related
objectives. In the main `refine3D` path, however, the update that finally lands
in the orientation table is a hard assignment, not a full soft responsibility
vector.

## Probabilistic Alignment Split

The probabilistic path is split into two phases:

1. `prob_align` or `prob_align_neigh`
   Sample particles for update, distribute probability-table work, aggregate all partitions, and write `ASSIGNMENT.dat`.
2. `strategy3D_prob`
   Consume the precomputed assignment map and stamp the chosen `(state, projection, in-plane angle, shift)` back into the orientation model.

This split matters. If a task changes how candidates are scored, sampled,
normalized, or selected, inspect both the table-building phase and the consumer
phase.

The probabilistic phase is particle-domain work. It may depend on prepared
`POLAR_REFS*` central-section files, but it should not take over
assembled-volume or polar-reference postprocessing.

## Probability Tables

`simple_eul_prob_tab.f90` is the core 3D probability-table implementation.

Its job is to build candidate tables over:

- state
- projection direction
- in-plane rotation
- optional shift

Important details:

- `fill_tab` evaluates per-reference candidates and can refine shifts with BFGS
  after identifying a probabilistic neighborhood.
- `fill_tab_state_only` handles `prob_state` mode, where the table only resolves state.
- `ref_assign` and `state_assign` convert aggregated tables into a single assignment map.
- `prob_athres` governs thresholded probabilistic sampling rather than exact posterior integration.

## Search Modes

The main refine modes relevant here are:

- `prob`
- `prob_state`
- `prob_neigh`
- greedy and SHC variants with probabilistic in-plane sampling
- tree-based variants that use probabilistic descent without the same exhaustive table shape

Do not assume all "probabilistic" modes share the same control flow. Some go
through `prob_align`; some perform probabilistic decisions inside search strategies.

## Regularization and Statistical Controls

Common controls that matter to the model:

- `objfun=euclid`
  Noise-normalized Euclidean objective with sigma estimation.
- `sigma_est`
  Controls sigma calculation policy.
- `ml_reg`
  Enables ML-style regularization in reconstruction-related phases.
- `lam_anneal`
  Applies lambda annealing for connectivity regularization.
- `prob_inpl`
  Enables probabilistic in-plane sampling in relevant search modes.
- even/odd partitioning
  Supports FSC- and SSNR-style estimates and guards against overfitting.

## Polar Reference Model

Non-probabilistic current-volume 3D matching can derive central sections
directly from the supplied `vol1..volN` source. `prob_tab` / `prob_tab_neigh`
and polar matcher preparation without a fresh complete volume source consume
`POLAR_REFS*` handoff files. Cartesian assembly refreshes
`POLAR_REFS_even.bin` / `POLAR_REFS_odd.bin` after volume assembly; polar
assembly writes the merged/even/odd `POLAR_REFS*` triplet after reducing polar
partial sums.

For multi-state runs, polar references are state-major. State 1 occupies
projections `1:nspace`, state 2 occupies `nspace+1:2*nspace`, and so on.
Common-line and obsfield restoration are intra-state operations; cross-state
common lines are not physical.

There is one bootstrap/materialization exception: when `POLAR_REFS*` are
missing or stale but all state volumes are present, `ensure_polar_refs_on_disk`
materializes them from volumes before probabilistic alignment or matching.
Callers must settle the PFTC k-range with `set_bp_range3D` before invoking the
materialization path.

## Relation To The Papers

The recent SIMPLE manuscripts describe two closely related but distinct ideas:

- the earlier probabilistic projection-matching framework for ab initio/refinement-style orientation recovery
- the newer volume-free polar-Fourier formulation, which removes intermediate volumes from the reprojection-estimation step in ab initio 3D reconstruction

For code work, treat the newer manuscript as conceptually adjacent rather than
identical to `refine3D`. It shares the fixed-point spirit and probabilistic
search philosophy, but its reprojection model is different.
