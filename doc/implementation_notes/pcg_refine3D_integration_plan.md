# PCG integration plan for `refine3D`

## Status and objective

The kernel PCG backend is implemented for shared and distributed
`reconstruct3D`. It is not yet a `refine3D` backend: parameter validation rejects
`rec_backend=pcg` outside `reconstruct3D`, the matcher writes gridding partials,
and `volassemble` only knows the gridding `(cmat,rho)` contract.

The target is an opt-in `rec_backend=pcg` path for `refine3D` and then
`refine3D_auto`, without changing the default gridding path. Production PCG
continues to require `pcgop=kernel`, fixed poses during reconstruction, a valid
`mskdiam`, and one to eight solver iterations (default two). Matrix-free remains
a small-fixture oracle and is never a workflow fallback.

The practical target is `refine3D_auto`. Consequently, cropped reconstruction,
the particle cache, fractional updates, trailing reconstruction, distribution,
and ML regularization are requirements rather than optional follow-up features.

| Phase | Deliverable |
| --- | --- |
| P0 | define and validate base-lambda scaling against effective data mass |
| 0 | repair the standalone gridding Euclidean/ML baseline and expose `box_crop` |
| 1 | prove cropped PCG scale equivalence in standalone `reconstruct3D` |
| 2 | generate and validate PCG base plus ML map pairs in standalone `reconstruct3D` |
| 3 | close distributed crop/ML manifest and restart gates; extract reusable services |
| 4 | implement raw `(B,D)` trailing chains and restart equivalence |
| 5 | add backend-aware `volassemble` restoration |
| 6 | wire and validate base `refine3D` |
| 7 | enable and benchmark the normal `refine3D_auto` lifecycle |

## Fixed design decisions

### Raw statistics and worker boundary

Each worker accumulates and atomically publishes raw, unregularized `(B,D)` for
each `(iteration,state,half,part)`:

```text
B = weighted adjoint-data accumulator
D = weighted Gram/sampling-density accumulator
```

Workers never fold, finalize a kernel, apply a prior, or solve. The master
validates the complete artifact set, reduces parts in ascending part order, and
is the only process that folds, finalizes, regularizes, and solves.

### `volassemble` boundary

`volassemble` remains the volume-domain workflow owner. It dispatches to a
backend-specific state restorer:

- gridding reduces `(cmat,rho)` and performs sampling-density correction;
- PCG reduces/blends `(B,D)`, finalizes the kernel and solves.

The backend seam is the dense even/odd halfmaps, not an attempt to disguise
`(B,D)` as gridding partials. After that seam, the existing FSC/cFAR, merged-map,
automask, nonuniform filtering, filenames, project updates, and next-iteration
handoffs remain common. PCG maps must not receive gridding correction or a
second sampling-density correction.

### Trailing domain

Trailing state is raw, unregularized PCG statistics. For realized sampled
fraction `f` and applied update weight `u`:

```text
B_chain_new = (u/f) B_current + (1-u) B_chain_previous
D_chain_new = (u/f) D_current + (1-u) D_chain_previous
```

The chain is written before kernel finalization or regularization. It is kept at
full-dataset mass. Scalar lambda and all FSC-derived prior terms are applied
only after the blend and are never persisted in the chain.

At a full-reconstruction stage boundary, PCG seeds the consuming trailing chain
at full mass. If no valid chain exists, the first fractional iteration mirrors
the existing gridding bootstrap: use the previous halfmaps for that iteration's
legacy volume-domain handoff and seed the new PCG chain with
`(B_current/f,D_current/f)` for subsequent iterations.

### Statistical weights versus priors

Particle/data weights, including `1/sigma2`, multiply both data statistics:

```text
B = sum_i w_i A_i^H y_i
D represents sum_i w_i A_i^H A_i
```

The zero-mean ML prior is different. It does not weight `B`; it adds prior
precision to the normal operator:

```text
(H_data + P_tau + lambda I) x = B
```

Only a nonzero prior mean would add a prior term to the right-hand side. For
parity with the current gridding path, `P_tau` is derived from the post-trailing
unregularized halfmap FSC and the data-density scale, then included in both the
PCG operator and its preconditioner. Persisted `D` remains data-only.

### Prior sequencing and the P0 lambda contract

The small base Tikhonov term and the FSC/SSNR ML prior have different status:

- `lambda I` is part of every PCG solve, including the `_unfil` base pair. Its
  scale must therefore be invariant to how the same data are partitioned,
  reduced, fractionally sampled, or blended into a trailing chain.
- `P_tau` is required gridding-parity behavior. It is derived only after the
  unregularized halfmaps and is applied to the second, ML-regularized solve.
- solvent-flatness and smoothness penalties are new research priors. They are
  excluded until hard-state `refine3D` parity is established.

P0 defines a deterministic effective-data scale `s_data(D)` from the
master-reduced or trailing-blended, data-only `D`:

```text
lambda_eff = lambda_rel * s_data(D_chain)
s_data(c D) = c s_data(D), for c > 0
s_data(D1 + D2) = s_data(D1) + s_data(D2)
```

The scale functional must use a declared, crop-aware support/band of `D` and
must reflect the statistical weights and symmetry expansion already present in
the data operator. It cannot depend on the number of files, workers, or parts.
Linearity makes partition reduction and trailing blends compose exactly, while
homogeneity makes duplication, uniform weight scaling, and `u/f` mass changes
scale lambda with the data term. Equivalently, the
solve may divide `B` and `H_data` by `s_data` and use `lambda_rel I`; the two
forms must be numerically equivalent. Workers continue to publish raw `(B,D)`
only. Lambda is derived by the master after reduction or trailing blend and is
never accumulated into `D` or stored in the chain.

P0 is closed only when duplication/data-weight rescaling, shared versus fixed-
order distributed reduction, repartitioning, `u/f` trailing updates, and crop
versus full-box matched-band tests preserve the same relative regularization
and reconstructed map. The current fixed lambda is retained only as a baseline
for choosing the dimensionless `lambda_rel`; it is not the production scaling
contract.

Soft state responsibilities do not block this work. Hard states close the base
parity gates first. When soft states are added during refinement semantics,
their responsibilities scale both `B` and `D`; `s_data(D)` then follows without
a separate state-weight path.

### ML/NU two-map contract

PCG preserves the existing two-pair `volassemble` contract. Both pairs come
from the same reduced and trailed raw `(B,D)`; they differ only by the ML prior:

| Product | PCG solve | Consumer |
| --- | --- | --- |
| `*_even_unfil.mrc`, `*_odd_unfil.mrc` | `H_data + lambda I`, with `P_tau=0` | FSC/cFAR and the base NU low-pass bank |
| standard even/odd halfmaps | `H_data + P_tau + lambda I` | normal state-map handoff and, when `nu_refine=no`, the ML auxiliary competitor |

Here `_unfil` means no FSC/SSNR ML prior. It still includes the declared small
PCG base lambda and the common real-space support constraint; those are part of
the PCG base estimator and are identical in the base and auxiliary solves.

Assembly for one state is therefore:

1. reduce and trailing-blend raw even/odd `(B,D)`;
2. persist the validated raw statistics until both output pairs are complete;
3. finalize and solve each half with `P_tau=0`, one half at a time, and write
   the `_unfil` pair;
4. calculate FSC/cFAR from that pair and retain this FSC for resolution/project
   metadata;
5. derive separate even/odd isotropic `P_tau` arrays from the shared FSC and
   each half's data-density scale, matching the current `add_invtausq2rho`
   policy and its low-resolution exclusion;
6. reload and deterministically re-finalize the same raw statistics for each
   half, add `P_tau` as a Fourier-diagonal positive precision term to the PCG
   normal operator and corresponding preconditioner approximation, and leave
   `B`, raw `D`, and the trailing chain unchanged;
7. warm-start from the corresponding `_unfil` map and solve each half again,
   writing the standard ML-regularized even/odd maps;
8. average the standard pair for the ordinary merged state map and pass both
   pairs to the unchanged NU orchestration.

This is one particle accumulation and two solve phases, not two particle
reconstructions. Particle I/O and raw reduction are performed once. Production
processes one half at a time to preserve the memory budget, so it cannot retain
both finalized half operators while waiting for the FSC. The correctness path
therefore replays deterministic kernel finalization from persisted raw `(B,D)`
for the ML solves. Retaining both operators is not an acceptable default; a
smaller finalized-operator cache may be considered only after measurement. Both
solves initially receive the same `pcg_maxits` bound; the warm-started ML solve
may be shortened only after a measured equivalence gate.

When `filt_mode=nonuniform*`, the NU bank is constructed from the `_unfil`
pair. With `ml_reg=yes` and `nu_refine=no`, the standard ML pair is supplied as
the auxiliary replacement candidate for the finest static bank member, exactly
as in the gridding path. With `nu_refine=yes`, the ML pair is still generated
as the standard reconstruction output but is not entered as an auxiliary
competitor; NU high-resolution shell extension owns that experiment. With
`ml_reg=no`, there is one effective pair, no auxiliary competitor, and the
standard maps equal the base PCG maps.

PCG outputs are already corrected according to the PCG operator. Neither the
base nor ML pair is multiplied by the gridding correction image or subjected
to gridding sampling-density correction before NU filtering.

### Sigma ordering

For `objfun=euclid`, iteration `n` has the explicit order:

1. search and commit iteration-`n` poses/states/shifts;
2. update and persist iteration-`n` per-particle sigma estimates;
3. accumulate iteration-`n` `(B,D)` using those in-memory sigma estimates;
4. reduce, trail, finalize, and solve.

The matcher path must not reload a previous sigma generation through the
standalone `reconstruct3D` loader. Artifact provenance records the sigma
generation/iteration. `objfun=cc` remains unweighted.

### Solver controls

`refine3D` already uses `maxits` for the number of refinement iterations. PCG
therefore receives distinct typed controls, `pcg_maxits` and `pcg_rtol`, with
defaults 2 and 0 respectively. The existing `reconstruct3D maxits/rtol` CLI can
remain as a compatibility interface.

## Implementation phases

Each phase has a test gate. A later phase does not begin until the preceding
gate passes.

### P0 prerequisite — data-mass-scaled base lambda

Implement the lambda contract above before crop/ML distribution is extended or
PCG is wired into `refine3D`. Freeze the exact `s_data(D)` support/band and
normalization in an operator-level test; derive `lambda_eff` only on the master
after fixed-order reduction or trailing blend.

Gate:

- one dataset reconstructed shared, with multiple partition counts, and after
  an equivalent rescaling of all data statistics has the same solution;
- fractional/trailing representations of the same full-mass statistics give
  the same `lambda_eff` and map;
- crop/full matched-band comparisons do not change the relative lambda merely
  because lattice size changed;
- diagnostics report `s_data`, `lambda_rel`, and the applied `lambda_eff` for
  every state and half.

### Phase 0 — repair the standalone gridding reference path

This establishes a fair Euclidean/ML and crop baseline before changing PCG.
The current gridding strategy attempts to read a partition sigma file directly,
whereas PCG uses the shared prior-run discovery/group loader. Gridding also
couples particle sigma weighting to `l_ml_reg`; the data objective must instead
be controlled by `objfun` and the tau prior independently by `ml_reg`.

Work:

- expose `box_crop` on the `reconstruct3D` CLI and derive `smpd_crop` from
  project geometry; do not expose an independent native-sampling override;
- use `load_sigma2_groups` for shared and distributed gridding whenever
  `objfun=euclid`, with the same missing-sigma hard error as PCG;
- remove the direct `SIGMA2_FBODY<part>` loading branches;
- apply `1/sigma2` particle weights whenever `objfun=euclid`, including
  `ml_reg=no`; let `ml_reg` control only FSC/SSNR prior addition during
  assembly;
- retain `objfun=cc` as the unweighted path;
- validate the existing gridding crop machinery and correct scale/CTF/shift
  handling if the acceptance test exposes a defect.

Gate:

- existing full-box CC results remain unchanged;
- `objfun=euclid,ml_reg=no` consumes prior sigmas and produces weighted base
  maps without an ML prior;
- `objfun=euclid,ml_reg=yes` produces `_unfil` and standard ML-regularized
  halfmaps without requiring `objfun=cc`;
- shared and distributed gridding agree;
- cropped and full gridding maps agree in their matched Fourier band and scale.

### Phase 1 — standalone PCG crop equivalence

Cropping is tested in `reconstruct3D` before any cache or refinement coupling.

Work:

- allow PCG to run at `box_crop,smpd_crop` while project metadata remains the
  authority for native `box,smpd`;
- keep standalone particle reads at native size, but limit the PCG plane,
  sigma, operator, mask, and output geometry to the cropped representation;
- record both native and cropped geometry in raw-artifact provenance;
- remove the crop hard error only after the acceptance gate passes.

Gate:

- cropped raw `B` and `D` agree with the overlapping band of a full-box build;
- cropped kernel agrees with matrix-free on a deterministic fixture;
- full-box and cropped PCG maps agree after matching resolution and scale;
- CC and Euclidean weights, shifts, astigmatic CTF, phase shift, C1, and C2/D2
  are covered;
- PCG and repaired gridding crop tests use the same fixture and comparison
  convention.

### Phase 2 — standalone PCG ML two-map path

This implements and validates ML regularization without `refine3D` or NU
integration.

Work:

- add the isotropic FSC/SSNR `P_tau` builder and positive Fourier-diagonal prior
  term to `reconstructor_pcg` and its preconditioner approximation;
- generate the `_unfil` pair from base PCG solves and calculate FSC/cFAR from
  that pair;
- deterministically replay each half from persisted raw `(B,D)`, add `P_tau`,
  warm-start from `_unfil`, and write the standard ML-regularized pair;
- keep the raw artifacts until both pairs and the merged map are published;
- make `ml_reg=no` execute only the base solve;
- retain a hard error for conical ML regularization.

Gate:

- `ml_reg=no` reproduces the current PCG result;
- the prior is symmetric positive semidefinite and passes matrix-free/kernel
  tests;
- small direct-system tests verify `(H_data+P_tau+lambda I)x=B`;
- shared `reconstruct3D` writes correctly named `_unfil`, regularized halfmaps,
  raw-FSC metadata, and the regularized merged map;
- the PCG and gridding two-map pairs are compared on the same Euclidean sigma
  input at full and cropped boxes.

### Phase 3 — close the distributed standalone boundary

The additive fixed-order reduction is already implemented. This phase extends
it through crop and ML operation, closes restart safety, and extracts the
services later used by refinement.

Work:

- add a manifest for the complete expected `(state,half,part)` artifact set;
- bind artifacts to run/iteration, native/cropped geometry, objective, symmetry,
  source, particle count, and sigma provenance;
- validate the complete set before reducing or solving anything;
- define restart behavior for valid, missing, stale, truncated, and mixed-
  generation artifacts;
- retain ascending-part association and master-only finalization;
- retain reduced raw statistics through both ML solve phases;
- extract reusable raw accumulation, reduction, base solve, prior, and ML solve
  services from standalone orchestration.

Gate:

- shared, `nparts=2`, and `nparts=3` produce equivalent raw statistics, base
  maps, FSC, and ML maps at full and cropped boxes;
- interrupted/restarted execution is equivalent to uninterrupted execution;
- corruption and mixed generations fail before a halfmap is published;
- disk, memory, replay-finalization cost, and cleanup are measured at a
  production-sized box.

### Phase 4 — PCG trailing-chain algebra and restart

Work:

- add raw-accumulator scale, add, import, export, and geometry-conversion APIs
  to `reconstructor_pcg` without exposing its internal arrays to workflow code;
- add one manifest-validated chain per state containing even and odd `(B,D)`;
- write the manifest last and invalidate it first on replacement;
- implement `u/f` plus `1-u` blending and full-mass seeding;
- carry complete chain sets between continued refinement directories;
- support a smaller previous grid by centered zero-padding when physical extent
  matches; discard and reseed on larger previous grids or extent mismatch;
- keep lambda and all prior precision outside the stored chain;
- optionally warm-start from the previous dense solution, but never use that
  solution as authoritative state.

Gate:

- deterministic algebra tests cover `f=1`, `u=f`, explicit `u`, repeated
  rounds, and convergence to the expected weighted history;
- artifact round trips and restart reproduce uninterrupted results;
- stage-boundary seeding and bootstrap normalization preserve full mass;
- box growth reproduces the overlapping-band result;
- shared and partitioned association orders reproduce the same chain.

### Phase 5 — backend-aware `volassemble`

Work:

- add a PCG state-restoration branch in `commander_volassemble` using the
  already validated standalone crop/ML services;
- detect populated/dropped states from PCG raw artifacts rather than gridding
  partial filenames;
- reduce and trail raw statistics, produce both map pairs, and retain raw FSC
  as the resolution authority;
- supply `_unfil` as the NU base and standard ML halfmaps as the auxiliary
  competitor when `ml_reg=yes,nu_refine=no`;
- average the standard halfmaps for the merged map; do not run a full-data
  solve;
- join existing automask, nonuniform filtering, project, and next-iteration
  handoffs without gridding correction on PCG maps.

Gate:

- standalone and `volassemble` invocation of the same PCG services produce
  identical map pairs and FSC;
- NU receives the correct base and auxiliary pairs in every `ml_reg`/
  `nu_refine` combination;
- standard filenames, merged maps, project fields, and postprocessing products
  follow the existing lifecycle contract.

### Phase 6 — wire base `refine3D`

Work:

- register `rec_backend`, `pcg_maxits`, and `pcg_rtol` for `refine3D`;
- allow `rec_backend=pcg` in typed validation with a precise supported-mode
  matrix;
- dispatch the post-search reconstruction boundary to the reusable PCG raw
  producer after poses and sigma are final;
- add a cache-aware PCG reader matching the existing gridding contract: no
  second normalization, taper/FFT at `box_crop`, cropped shifts/CTF sampling,
  and sigma shells capped at crop Nyquist;
- preserve one particle read in the reconstruction phase and the same producer
  for shared/distributed workers;
- update partial cleanup, dropped states, continue/carry-over, benchmarking,
  and failure cleanup for PCG artifacts;
- keep gridding defaults unchanged.

Gate:

- at least two iterations complete in shared and distributed modes;
- cached/uncached cropped accumulation agrees for CC and Euclidean weights;
- current-iteration sigma perturbations affect current-iteration `(B,D)` and
  never lag by one iteration;
- multi-state hard assignments, empty-state carry-forward, standard project
  outputs, and shared/distributed parity pass.

### Phase 7 — enable `refine3D_auto`

Work:

- expose and propagate backend and PCG solver controls through
  `refine3D_auto`;
- run the normal autoscale, cache, fractional-update, trailing, ML, automask,
  and nonuniform-filter lifecycle with PCG;
- retain hard errors for `projrec=yes`, conical regularization, matrix-free
  workflow execution, and unsupported high-order symmetry performance cases;
- add phase/resource diagnostics suitable for multi-iteration runs.

Gate:

- an autoscaled run with `update_frac<1`, `trail_rec=yes`, `ml_reg=yes`, cache
  enabled, and `nparts>1` completes and restarts cleanly;
- half independence, FSC history, resolution metadata, automasks, NU products,
  and final reconstruction are valid across iterations;
- matched PCG/gridding experiments compare FSC, cFAR, scale, anisotropy,
  runtime, peak RSS, and disk use;
- macOS/Linux and debug/release coverage passes before production promotion.

## Primary code ownership

| Area | Owning code |
| --- | --- |
| typed controls and validation | `src/main/params/simple_parameters*.f90` |
| `refine3D`/`refine3D_auto` UI | `src/main/ui/simple/simple_ui_refine3D.f90` |
| post-search backend dispatch | `src/main/strategies/search/simple_strategy3D_matcher.f90` |
| cache-aware partial production | `src/main/strategies/search/simple_matcher_3Drec.f90` plus reusable PCG accumulation services |
| raw PCG algebra and solver | `src/main/volume/simple_reconstructor_pcg.f90` |
| raw filenames and manifests | `src/defs/simple_refine3D_fnames.f90` |
| shared/distributed iteration policy | `src/main/strategies/parallelization/simple_refine3D_strategy.f90` |
| backend restoration and common handoff | `src/main/commanders/simple/simple_commanders_rec_distr.f90` |
| operator and workflow gates | `src/main/commanders/test/simple_commanders_test_highlevel.f90` and end-to-end fixtures |

## Initial support and deferred work

The first production-usable PCG refinement path supports hard state/pose
assignments, raw or validated cached particles, isotropic ML regularization,
multiple states, C1/C2/D2, cropped reconstruction, fractional updates, trailing,
and shared/distributed execution.

Deferred behind that gate:

- conical/anisotropic FSC regularization;
- `projrec=yes` compression;
- soft state responsibilities (Phase 4 refinement semantics, after hard-state
  parity; responsibilities must weight both `B` and `D`);
- solvent-flatness and smoothness quadratic priors;
- symmetry projection/permutation optimization pending distributed profiling;
- high-order symmetry promotion beyond measured performance limits;
- nonlinear priors, non-negativity, TV, or wavelet penalties;
- matrix-free workflow fallback;
- GPU/offload implementation.

## Approval checkpoint

Coding should begin only after approval of these five contracts:

1. `objfun=euclid` controls sigma-weighted data accumulation in both backends;
   `ml_reg` independently controls addition of the FSC/SSNR prior.
2. Crop and ML behavior are first implemented and gated in standalone
   `reconstruct3D`, using repaired gridding as the reference workflow.
3. `volassemble` orchestrates backend restoration and common postprocessing;
   the PCG-specific seam is dense halfmaps.
4. PCG trailing state is manifest-validated raw `(B,D)` at full mass, with all
   priors applied only after blending.
5. Isotropic ML regularization is a separate positive precision operator,
   derived from unregularized post-trailing FSC, not a modification of stored
   `B` or `D`. `pcg_maxits`/`pcg_rtol` remain distinct from refinement-loop
   `maxits`.
6. The only prior-related prerequisite before integration is P0 base-lambda
   scaling. `P_tau` is required gridding parity; solvent-flatness and smoothness
   priors are post-parity research and cannot be used to close integration
   gates.
