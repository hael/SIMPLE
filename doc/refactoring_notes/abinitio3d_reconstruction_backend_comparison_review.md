# Abinitio3D Reconstruction Backend Comparison Review

Date: 2026-09-09

Status: review complete; two comparison blockers remain before gridding and
PCG quality or cost results should be treated as decisive.

Purpose: single living record of the review findings, required refactoring,
comparison contract, and validation gates for `rec_backend=gridding` versus
`rec_backend=pcg` in `abinitio3D`.

## 1. Executive Verdict

The recent changes materially improve comparability:

- both backends evaluate the reconstructed even/odd half maps through the
  shared half-map diagnostics path;
- gridding now evaluates the deapodized, shipped half maps rather than a
  different intermediate representation;
- the gridding and PCG ML-prior formulas use the same FSC-to-SSNR convention
  and compatible scaling;
- both master reconstruction paths use the same master thread-budget policy;
- the default spherical-support and no-automask paths are closer in output and
  postprocessing behavior.

The comparison is not yet controlled, however. Two issues can change the
scientific result independently of the reconstruction backend:

1. cropped particle observations are prepared differently for gridding and
   PCG;
2. support-provenance sidecars do not follow reconstructed volumes through
   every rename, copy, publication, and restart transition.

These are blockers for a direct backend comparison. Cost instrumentation and
the mathematical meaning of the support also need clarification, but those
can be addressed without changing the core reconstruction algorithms.

## 2. What Is Already Aligned

The following conditions are sufficiently shared for the proposed comparison:

- stages 1 and 2 use gridding for both runs; PCG begins at stage 3;
- particle half assignment, orientation metadata, CTF metadata, sampling and
  fractional-update bookkeeping use the same workflow state;
- both paths consume the same sigma model selected by the workflow;
- radial FSC and cFAR are calculated through `evaluate_halfmap_pair`;
- gridding halves are deapodized before they are shipped and evaluated;
- the ML regularization expressions implement the same shell-average
  `rho / (tau * SSNR)` model, including the same FSC clamp and high-pass limit;
- the master reconstruction thread budget is selected through
  `rec3D_master_nthr` for both backends.

This is a sound base. The remaining differences must either be removed or be
explicitly classified as part of the backend being measured.

## 3. Blocking Findings

### 3.1 P1: cropped particle preparation is not backend-neutral

The gridding matcher and the PCG reconstruction strategy do not currently
construct the same cropped observation when `box_crop < box`.

The gridding path normalizes at the native box, Fourier-crops the particle,
returns to real space, applies the edge taper at the cropped box, and then
pads/transforms the observation for reconstruction. The shared and distributed
PCG paths normalize and taper at the native box, transform, and select the
native Fourier plane needed by the cropped solve.

Cropping and tapering do not commute. The two backends therefore receive
different pixel values before CTF, sigma, interpolation, or solver behavior
can be compared. This matters particularly in staged `abinitio3D`, where
`box_crop` changes across stages.

Required correction:

- introduce or reuse one backend-neutral particle-preparation operation that
  produces the reconstruction-ready cropped observation;
- call it from the gridding matcher and from both shared and worker PCG
  accumulation;
- preserve the current matcher single-read contract: the comparison fix must
  not introduce a second particle-stack read;
- add a diagnostic that compares prepared observations before backend-specific
  weighting or accumulation.

Until this is fixed, an observed FSC, success-rate, or map-quality difference
cannot be attributed solely to gridding versus PCG.

### 3.2 P1: support provenance is not managed as part of the volume artifact

The support sidecar is now shared by both backends, but the lifecycle is still
volume-file-centric. Several `abinitio3D` stage and final-output transitions
rename or copy the MRC file without applying the equivalent operation to its
support-provenance file.

Consequences include:

- a valid PCG half map can lose its provenance at a PCG-to-PCG stage boundary,
  forcing an unintended cold base solve;
- an old destination sidecar can survive beside a newly copied volume and
  describe the wrong artifact;
- a final copied gridding map can lack provenance, allowing later standalone
  postprocessing to apply the spherical support a second time;
- writing the sidecar before the MRC is successfully published can leave an
  orphan sidecar after interruption.

With the current two-iteration PCG configuration, an unintended cold start is
large enough to affect both quality and runtime.

Required correction:

- treat the volume and support sidecar as one reconstruction artifact;
- provide narrow copy, rename, remove, and publish operations for that
  artifact, rather than duplicating filename logic in controllers;
- use those operations for random-start setup, stage snapshots, final copies,
  restart handling, and cleanup;
- publish the MRC first and its sidecar second;
- remove a stale destination sidecar when the source artifact has none;
- retain the current sidecar filename for compatibility.

## 4. Additional Review Findings

### 4.1 P2: current timing answers wall time, not full compute cost

The new reconstruction-master timer usefully isolates master assembly. The
existing iteration benchmark also reports total refinement, matcher/scheduler,
and assembly/postprocessing wall time. Together these are adequate for a
wall-clock comparison on the same hardware.

They do not yet provide a defensible compute-cost or memory comparison:

- worker reconstruction timing is emitted only by partition 1;
- load imbalance across partitions is invisible;
- process peak RSS is lifetime-wide and is sampled for only that worker;
- master peak RSS is not recorded;
- elapsed seconds are not converted to thread-seconds or core-hours.

Required reporting for a cost claim:

- collision-free per-part worker timing, or equivalent scheduler accounting;
- maximum worker wall time and the worker-time distribution;
- sum of `worker seconds * worker threads`;
- `master seconds * master threads`;
- master and maximum-worker peak RSS;
- stage, iteration, backend, `box`, `box_crop`, `nparts`, worker/master thread
  counts, PCG `maxits` and `rtol`, and hardware identity.

Wall time, compute cost, and memory must be reported as separate quantities.

### 4.2 P2: the two support formulations are not proven mathematically equal

The gridding result is explicitly multiplied by the cosine-edged spherical
support after restoration. PCG solves a support-constrained system of the form
`(P H P)u = P b` and returns `x = P u`. Its warm-start conversion divides by
`P` above a safety threshold and sets lower-support voxels to zero, preventing
the soft edge from being multiplied repeatedly across iterations.

These constructions share a zero-support boundary, but they do not by
themselves prove identical amplitudes through the cosine edge. Where `P` is
nonzero, the solved variable can compensate for its magnitude. It is therefore
too strong to claim that the backends have an identical soft-support operator.

The current no-remasking warm-start direction is appropriate. The remaining
choice is policy:

- accept support behavior as part of each backend's reconstruction bundle; or
- use a hard solve-domain support and apply one common soft output/evaluation
  window after both reconstructions.

Whichever policy is selected must be documented and checked with a radial
edge-profile regression and a repeated-warm-start no-compounding regression.

### 4.3 P2: the experiment must lock all non-backend choices

The default sphere-only path is the cleanest first comparison. Automasking,
explicit PCG masks, and conical FSC can select backend-specific support or
regularization behavior and should not be mixed into the baseline experiment.

## 5. Required Baseline Comparison Profile

For the first quality and cost study, use:

```text
automsk=no
envfsc=no
conical_fsc=no
pcg_mskfile=<unset>
```

Lock the following between paired runs:

- starting project and initial maps;
- random seed;
- even/odd assignment;
- sampled particle cohort and fractional-update state;
- sigma state;
- stage schedule and all resolution limits;
- `nparts`, worker threads, master threads, queue, host type, and hardware;
- all reconstruction-independent parameters.

Fix and record PCG `maxits`, `rtol`, and operator mode. Run each backend in an
independent clone of the starting project; do not run them sequentially in one
mutable project directory.

## 6. Implementation Plan

1. **Unify observation preparation.** Add the backend-neutral cropped-particle
   preparation boundary and route gridding plus shared/distributed PCG through
   it.
2. **Make provenance transactional.** Add artifact-level lifecycle helpers and
   update all stage, final-output, restart, and cleanup transitions.
3. **Define support semantics.** Record whether support is an intentional
   backend difference or a common post-solve operation, then add the relevant
   regressions.
4. **Complete measurement.** Add per-part or scheduler-derived cost data,
   master/worker memory data, and a machine-readable comparison manifest.
5. **Run a reconstruction-only gate.** With fixed orientations and sigma,
   compare one `reconstruct3D` execution per backend before the full
   `abinitio3D` study.
6. **Run paired abinitio3D replicas.** Use identical starting states and seeds;
   evaluate quality, success rate, wall time, compute cost, and memory.

## 7. Validation and Acceptance Gates

Before interpreting backend results, verify:

1. **Prepared-observation parity:** for `box_crop < box`, compare the complex
   Fourier observations before backend-specific weighting or accumulation.
2. **Raw-input parity:** with fixed particles and orientations, compare CTF,
   sigma, shell weights, symmetry replication, and the observation cohort.
3. **Shared/distributed parity:** for each backend, compare `nparts=1` with a
   fixed multi-part execution using the same input and reduction order where
   practical.
4. **Support behavior:** compare radial RMS through the taper and verify that
   repeated warm starts do not progressively shrink the edge.
5. **Artifact lifecycle:** exercise stage rename, snapshot copy, final copy,
   restart, missing-sidecar, and stale-destination-sidecar cases.
6. **Reproducible cost:** confirm that reported wall time and thread-seconds can
   be recomputed from the emitted records.

The comparison is ready when:

- prepared input parity passes at an agreed numerical tolerance;
- no valid reconstruction loses provenance and no stale provenance survives;
- the configuration differs only by `rec_backend` and explicitly recorded PCG
  solver parameters;
- wall time, compute cost, and peak memory are independently reproducible;
- quality conclusions are based on paired replicas rather than a single run.

For the production quality study, use at least 10 paired seeds; 20 is preferred
if the success-rate difference is small. Record final FSC 0.5 and 0.143
crossings, cFAR, map correlation or ground-truth metrics when available,
success/failure classification, per-stage timings, core-hours, and peak RSS.

## 8. Recommended Decision

Proceed with the comparison refactoring now, but do not start the definitive
quality/cost campaign until Sections 3.1 and 3.2 pass their validation gates.
The first implementation should be deliberately narrow: shared observation
preparation, transactional provenance handling, and complete measurement.
Changes to the numerical gridding or PCG algorithms are not required to make
the initial experiment meaningful.
