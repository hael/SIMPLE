# Abinitio3D Importance Sampling and Fractional Updates

## Scope

This note explains how SIMPLE couples:

- outer particle-subset selection driven by `update_frac` or `nsample*`
- inner importance-style sampling used by probabilistic orientation assignment
- downstream reconstruction behavior that depends on the realized subset

Keep these layers separate. The code is easier to understand when split into
policy, bookkeeping, assignment, and reconstruction.

## Core Contract

### 1. Ab initio owns the stage policy

`src/main/commanders/simple/simple_commanders_abinitio.f90` computes the base
`update_frac` before the stage loop.

Important details:

- `nsample` is converted to `update_frac` by scaling with `nstates`.
- explicit `update_frac` is accepted directly.
- `nsample_start` and optional `nsample_stop` trigger dynamic scheduling through
  `calc_update_frac_dyn`.
- otherwise the default min and max sample targets are used.
- the result is capped by `UPDATE_FRAC_MAX = 0.9`, so ab initio keeps
  fractional update enabled.

That same command writes class-sampling statistics to `CLASS_SAMPLING_FILE`,
which later balanced sampling reuses.

### 2. The stage controller changes search mode and update policy together

`src/main/simple_abinitio_controller.f90` emits the per-stage `refine3D`
command line.

Important stage transitions:

- stages 1 to 4 use `refine='shc_smpl'`
- stages 5 to 6 use `refine='prob'`
- stages 7 to 8 use `refine='prob_neigh'`
- stage 8 can turn on `fillin='yes'` for single-state ab initio

Important coupling facts:

- dynamic `update_frac` is recomputed per stage when the run uses dynamic sampling
- `frac_best` becomes less greedy once probabilistic refinement starts
- `trail_rec` is introduced only in later stages

This is the main place where the code couples search strategy, sampling size,
and reconstruction behavior.

## Bookkeeping Layer: `simple_oris`

### 3. `sampled` and `updatecnt` are the handshake

`src/main/ori/simple_oris_sampling.f90` defines the subset-management routines.

Key routines:

- `sample4update_class`: balanced class-aware sampling using `CLASS_SAMPLING_FILE`
- `sample4update_cnt`: prefers particles with lower `updatecnt`
- `sample4update_fillin`: fills in particles with low update history late in the run
- `sample4update_reprod`: reproduces the most recent sampled subset exactly
- `sample4update_updated`: retrieves all particles updated in a previous step

Interpretation:

- `sampled` identifies the current sampling round
- `updatecnt` counts how often a particle has been selected for update
- `incr_sampled_updatecnt` advances both the round marker and per-particle
  history when requested

### 4. Realized update fraction is recovered from bookkeeping

`src/main/ori/simple_oris_getters.f90` implements `get_update_frac`.

It does not return the command-line target directly. Instead it computes:

- numerator: active particles in the latest `sampled` round
- denominator: active particles with `updatecnt > 0`

So `get_update_frac` is the realized fraction of the currently updated subset
within the active updated pool.

## Update Selection During Refinement

### 5. Non-probabilistic refine samples directly

`src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` provides
`sample_ptcls4update3D`.

Behavior:

- if `l_update_frac` is false, sample all active particles
- if `balance=yes`, use `sample4update_class`
- otherwise use `sample4update_cnt`
- if late-stage fill-in is enabled, `sample_ptcls4fillin` redirects to
  `sample4update_fillin`

### 6. Probabilistic refine samples once, then reproduces

`src/main/commanders/simple/simple_commanders_prob.f90` is the key coupling point.

In `exec_prob_align` and `exec_prob_align_neigh`:

- first iteration can clean `updatecnt` and `sampled`
- the active subset is chosen by `sample_ptcls4update3D` or `sample_ptcls4fillin`
- the chosen subset is written back to the project file
- partition jobs then compute probability tables only for that subset

In `exec_prob_tab` and `exec_prob_tab_neigh`:

- no new sampling happens
- `sample4update_reprod` reproduces the subset created by `prob_align`

In `src/main/strategies/search/simple_strategy3D_matcher.f90`:

- when `l_prob_align_mode` is true, `refine3D_exec` also uses `sample4update_reprod`
- this ensures search and assignment operate on the exact same particles
- when current partial reconstruction is active, the matcher should reuse the
  particle images read during batch construction instead of re-reading stacks in
  a later reconstruction pass

This is the most important contract in the whole workflow.

## Inner Importance Sampling

### 7. Probabilistic tables do their own candidate sampling

`src/main/simple_eul_prob_tab.f90` and
`src/main/simple_eul_prob_tab_neigh.f90` implement the inner
importance-style sampling logic.

Important points:

- `angle_sampling` wraps `greedy_sampling`
- distances or correlations are transformed into a positive proxy
- the code samples likely in-plane or reference candidates instead of
  exhaustively committing to a single best candidate immediately
- `prob_neigh` tightens this around neighborhood-style candidate sets

This sampling is not the same as the outer `update_frac` selection:

- outer selection chooses which particles participate this iteration
- inner sampling chooses which candidate orientations or references are explored
  for each participating particle

## Reconstruction Consequences

### 8. `volassemble` consumes the realized fraction

`src/main/commanders/simple/simple_commanders_rec_distr.f90` does not choose particles.

Its coupling to sampling appears when `trail_rec=yes`:

- if `ufrac_trec` is provided, that value is used
- otherwise it reads `ptcl3D` and calls `get_update_frac`
- previous even and odd volumes are weighted by `1 - update_frac_trail_rec`
- current even and odd volumes are weighted by `update_frac_trail_rec`

So the realized particle subset affects how strongly the new reconstruction
replaces the previous one.

Online current reconstruction is still produced by the matcher path. Preserving
one particle-stack read per active batch is part of the performance contract:
batch preparation keeps the raw images needed for reconstruction, assignment is
made, and reconstruction preparation consumes those same images. Do not move
this into a second image-reading pass unless the repository policy changes.

## Practical Reading Order

For explanations:

1. `simple_commanders_abinitio.f90`
2. `simple_abinitio_controller.f90`
3. `simple_matcher_smpl_and_lplims.f90`
4. `simple_oris_sampling.f90`
5. `simple_commanders_prob.f90`
6. `simple_strategy3D_matcher.f90`
7. `simple_eul_prob_tab*.f90`
8. `simple_commanders_rec_distr.f90`

For modifications:

1. change policy in ab initio or stage-controller code
2. confirm `simple_oris` bookkeeping still matches the intended subset semantics
3. confirm probabilistic paths reproduce rather than resample
4. confirm trailing reconstruction still uses the intended realized fraction
5. confirm online reconstruction/restoration still reuses matcher batch images
