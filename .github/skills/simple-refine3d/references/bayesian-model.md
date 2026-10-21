# Probabilistic 3D Refinement Model

## Scope

Use this note when the task is about the scientific or algorithmic meaning of
SIMPLE's probabilistic 3D refinement implementation.

## Core Model

SIMPLE's 3D refinement is a fixed-point workflow that alternates between:

- updating particle orientations, in-plane angles, shifts, and sometimes states
- writing particle-derived Cartesian partial reconstructions from those assignments
- assembling and postprocessing Cartesian state volumes in volume assembly
- projecting the next iteration's reprojection model from those volumes

The implementation uses probability-like candidate generation, thresholded
sampling, regularization, even/odd cross-validation, and likelihood-related
objectives. In the main `refine3D` path, however, the update that finally lands
in the orientation table is a hard assignment, not a full soft responsibility
vector. The one exception is reconstruction under `m_estimator=flex`, which
weighs particles by frozen `flex_pca` state weights that alignment never reads.

## Probabilistic Alignment Split

The probabilistic path is split into two phases:

1. `prob_align_exec` or `prob_align_neigh_exec` (`simple_prob_align_strategy.f90`)
   Sample particles for update once, run the table workers (`prob_tab`,
   `prob_tab_neigh`) on that sample, read all partitions, assign, and write
   `assignment_part.dat`.
2. `strategy3D_prob`
   Consume the precomputed assignment map and stamp the chosen `(state, projection, in-plane angle, shift)` back into the orientation model.

This split matters. If a task changes how candidates are scored, sampled,
normalized, or selected, inspect both the table-building phase and the consumer
phase.

The probabilistic phase is particle-domain work. It reads the strategy's
reprojection model and must not take over assembled-volume postprocessing.

## Probability Tables

`src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` is the
core 3D probability-table implementation.

Its job is to build candidate tables over:

- state
- projection direction
- in-plane rotation
- optional shift

Important details:

- `fill_tab_range` evaluates per-reference candidates for a batch of the
  sample and can refine shifts (L-BFGS-B, or a quadratic fit to the coarse
  lattice in multi-state tables) after identifying a probabilistic
  neighborhood. Shifts are held in the particle frame.
- `ref_assign` converts the aggregated table into a single assignment map; the
  state label is the deterministic argmin over states, and only the projection
  within the state is drawn (likelihood weights `exp(-dist)` over a top-K
  support, `sample_likelihood_dist`).
- `prob_athres` bounds the angular range of the thresholded candidate draws
  (`angle_sampling`, `sample_bounded_dists`); it is not exact posterior
  integration.

## Search Modes

The `refine` modes of the 3D matcher:

- `prob` (dense table) and `prob_neigh` (neighbourhood table;
  `prob_neigh_mode=state|geom|shc|snhc`), both consumed by `strategy3D_prob`
- `shc` (greedy for unsearched particles and at `GREEDY_FREQ`), `greedy`,
  `greedy_inpl`, `neigh` (`strategy3D_greedy_sub`, which can sample in-plane
  angles under `prob_inpl`)
- `eval` and `sigma` (no pose search)
- `cont` (Cartesian continuous pose refinement)

Do not assume all "probabilistic" modes share the same control flow. Some go
through `prob_align`; some perform probabilistic decisions inside search strategies.

## Regularization and Statistical Controls

Common controls that matter to the model:

- `objfun=euclid`
  Noise-normalized Euclidean objective with a canonical sigma2 store; each
  iteration scores with the committed sigmas and commits new residual sigmas
  after assembly (`commit_sigma2_exec`).
- `ml_reg`
  Enables ML-style regularization in reconstruction-related phases.
- `lam_anneal`
  Applies lambda annealing for connectivity regularization.
- `prob_inpl`
  Enables probabilistic in-plane sampling in relevant search modes.
- even/odd partitioning
  Supports FSC- and SSNR-style estimates and guards against overfitting.

## Reprojection Model

The iteration strategy materializes one reprojection model per iteration from
the current Cartesian volumes (`materialize_reprojection_model_from_volumes`,
`simple_matcher_refvol_utils.f90`): reference read, masking/filtering/centering
policy, band selection (`set_bp_range3D`), and even/odd polar sections written
to `reprojection_model_{even,odd}.bin`. The matcher and the table workers read
it (`read_reprojection_model`) and adopt its band from the header. Volume
assembly produces Cartesian volumes only; it does not generate or refresh
matcher references.

For multi-state runs, references are state-major. State 1 occupies
projections `1:nspace`, state 2 occupies `nspace+1:2*nspace`, and so on.

## Relation To The Papers

The recent SIMPLE manuscripts describe two closely related but distinct ideas:

- the earlier probabilistic projection-matching framework for ab initio/refinement-style orientation recovery
- the newer volume-free polar-Fourier formulation, which removes intermediate volumes from the reprojection-estimation step in ab initio 3D reconstruction

For code work, treat the newer manuscript as conceptually adjacent rather than
identical to `refine3D`. It shares the fixed-point spirit and probabilistic
search philosophy, but its reprojection model is different: `refine3D` always
projects its references from Cartesian volumes.
