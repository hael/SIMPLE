# `flex_analysis` Policy: Sparse Diffusion-Manifold 3D Pre-Images

This policy records the implemented scientific and software contract for
`flex_analysis`. Forward-facing developments are separated explicitly from
the production behavior.

## 1. Production model

`flex_analysis` is a diffusion-map workflow. It does not fit PCA, probabilistic
PCA, alternating E/M updates, a Gaussian latent prior, an external sigma array,
or a model noise variance.

The implemented input assumptions are:

- one supplied and fixed mean volume (`vol1`);
- one active `ptcl3D` state;
- fixed 3D orientations, in-plane rotations, shifts, and CTF parameters;
- an explicit projection-grid size (`nspace`);
- `ctf=yes` or already phase-flipped (`ctf=flip`) particle images;
- one low-pass limit `lp` used only for graph-feature construction.

Before feature preparation, the implementation requires at least three
selected active particles, matching ptcl2D/ptcl3D field sizes, valid and unique
selected stack slots, and a uniform supported CTF mode. Mixed
`ctf=yes`/`ctf=flip` selected particles are rejected. `nspace` must be supplied
on the command line: the maximum populated `ptcl3D%proj` value cannot recover
the grid size because the final projection directions may be empty.

The standalone workflow returns a continuous low-dimensional coordinate system
and representative nonlinear 3D pre-image state volumes. The default number
of state volumes is `npreimages=8`; callers may override it explicitly. This is
separate from `nstates`, which remains one because the input is a single-state
`ptcl3D` project. Callers that only consume the embedding skip pre-image
selection and reconstruction.

## 2. Canonical registered frame

Each selected active particle is geometrically registered once for feature
preparation by the inverse of its stored in-plane shift and in-plane angle.
Feature preparation includes the standard particle normalization performed by
`transform_ptcls`. For the current diagnostic workflow, the residual 3D fit
then reads `flex_registered_particles.simple` directly through the
projected-model batch path. No registered-to-native project is created and no
orientation, particle-index, or image mapping is performed. The implementation extends and
uses `transform_ptcls` with an explicit
particle-index batch, following the matcher batch policy
`min(nptcls,nthr*BATCHTHRSZ)` and the established project-aware particle
reader. The explicit `nspace` value defines the builder grid and the stored
project `proj` field is authoritative. Every selected `proj` must lie in
`1:nspace`; an out-of-range value is a hard error and is never remapped or
clamped. After registration the particle is directly comparable with the
corresponding discrete mean-volume reprojection; no rotation or shift search
is performed during graph construction.

The registered particles are phase-flipped when the input project is
`ctf=yes`. Already phase-flipped particles are not flipped again. The durable
outputs are:

- `flex_registered_particles.mrcs` in shared-memory mode, or
  `flex_registered_particles_partNNN.mrcs` stacks in distributed mode;
- `flex_registered_particles.simple`;
- `flex_registered_particle_map.txt`.

Each stack is written incrementally with one process-level `stack_io` open and
close; the implementation does not open an MRC stack per particle. In a
distributed run, the master generates the discrete mean-volume projections
once with `simple_simple_volinterp::reproject`, writes one temporary projection
stack, and feature workers read it with buffered `stack_io`; workers do not
repeat the reprojection calculation.

The registered stack or part stacks are compact and collectively contain only
the selected particles in feature-row order. The copied project maps selected
`ptcl3D` rows to compact stack indices, marks unselected `ptcl3D` rows inactive,
records every stack as `ctf=flip`, and stores zero in-plane angle and shift in
the registered `ptcl3D` field. The `ptcl2D` field remains untouched. Defocus,
phase shift, and other CTF
metadata remain available. For
astigmatic CTFs, `angast` is rotated into the canonical registered frame by
adding the removed particle `e3` modulo 180 degrees; leaving the raw
astigmatism axis unchanged would be incorrect. An existing `frc3D` record is
preserved when available.

The diagnostic pre-image reconstruction uses the registered project unchanged.
Native-image reconstruction is deliberately disabled until the pre-image path
has passed its reconstruction identity checks.

## 3. Graph feature and sparse candidate policy

For particle `i`, the graph feature is the masked, band-limited registered
residual

```text
r_i = register(phaseflip(y_i)) - |CTF_i^c| P^c_proj(i)(mean)
```

Here the phase-flip operation is omitted for input already marked `ctf=flip`,
and superscript `c` denotes the canonical registered frame. The resulting
residual is serialized in Cartesian image space. When active, the configured
mean-volume mask is applied before the discrete projections are generated.
Each residual is then restricted by the active `lp` and 2D soft mask before
serialization. The comparison is therefore CTF-consistent in the phase-flipped
observation model. This graph-feature limit is not applied to the generative
volumes.

Particles are binned by the validated discrete projection index from the
builder's symmetry-reduced `eulspace`. Projection directions are taken from
that same `eulspace`. For every occupied source direction, projection bins are
traversed in dot-product angular order. Particles belonging to a source bin are
OpenMP-parallel owner rows: each thread updates only its particle's top-k list.
At most `nang_nbrs` candidate particles are compared and only the `k_nn`
closest registered-residual features are retained. The default public
parameters are:

```text
k_nn      = 10
nang_nbrs = 100
```

`nang_nbrs` is a candidate cap, not the number of retained graph edges. It is
never represented as a dense `N x nang_nbrs` table. The retained directed
neighbors are symmetrized into the existing CSR `diffmap_graph` and normalized
with the existing diffusion-map policy.

The default candidate pool is ten times the default retained-neighbor count.
For a roughly one-particle-per-direction data set with 5000 sampled directions,
this is about 2% of the direction space rather than 20%. Because the cap counts
particles rather than direction bins, callers with unusually high or low bin
occupancy may override it explicitly.

Angular proximity gates which comparisons are meaningful; registered image
distance decides which gated particles are neighbors.

### 3.1 Graph kernel bandwidth and density normalization

The Gaussian graph kernel uses a single fixed convention,
`w_ij = exp(-d_ij^2 / eps)`, where `eps` is the squared-distance bandwidth. The
kernel, both bandwidth estimators, and the Ferguson log-bandwidth scan all use
this same convention, so the estimated optimum is directly the kernel
denominator. The kernel and normalization live in the shared engine
`src/main/pca/simple_diff_map_graphs.f90` and are used by every diffusion-map
program (`ppca_denoise`, `cls_split`, `denoise_project`, `flex_analysis`).

`bandwidth_mode` selects the estimator for `eps`:

- `median` — `eps` is the median of the k-th-nearest-neighbor squared
  distances (`median_positive(kth_d2)`), with a mean fallback on underflow.
- `ferguson` — the Coifman–Lafon / Ferguson tanh selector: a log-bandwidth
  scan over `[d2med·10^-3, d2med·10^3]` computes `log sum_ij exp(-d_ij^2/eps)`
  at each trial `eps`, a `d + c·tanh(a x + b)` model is fitted, and the
  inflection `logeps* = -b/a` (clamped to the scanned range) gives
  `eps = bandwidth_tune · exp(logeps*)`. A non-finite, non-positive, or
  off-scale result (outside `[d2med/100, 100·d2med]`) is rejected and the
  estimator falls back to `median`, so a degenerate fit degrades gracefully
  instead of collapsing the embedding.

`flex_analysis` defaults to `bandwidth_mode=ferguson`; the other diffusion-map
programs default to `median`. `bandwidth_tune` is a linear multiplier of the
Ferguson optimum (default `1.0`; larger broadens, smaller sharpens) and has no
effect under `median`.

Density normalization is controlled by `dm_alpha` (Coifman–Lafon α), applied in
`normalize_diffmap_graph` before the symmetric degree normalization:

- `dm_alpha = 0.0` (default) — plain graph Laplacian; sampling density retained.
- `dm_alpha = 0.5` — Fokker–Planck diffusion.
- `dm_alpha = 1.0` — Laplace–Beltrami; sampling density divided out, leaving
  intrinsic manifold geometry.

For `dm_alpha > 0` the weights are rescaled
`w_ij <- w_ij / (d_i^alpha d_j^alpha)` on the raw degree and then symmetrically
renormalized; `dm_alpha = 0` reproduces the original degree normalization
exactly. Because cryo-EM orientation coverage is non-uniform (preferred
orientations, angular gaps), `dm_alpha = 1` is recommended when intrinsic
geometry is desired; the default remains `0.0` for backward compatibility and
must be enabled explicitly.

Example invocations:

```bash
# Ferguson data-adaptive bandwidth (flex default), explicit tune
simple_exec prg=flex_analysis ... bandwidth_mode=ferguson bandwidth_tune=1
# Laplace–Beltrami density normalization
simple_exec prg=flex_analysis ... dm_alpha=1
```

## 4. Embedding and feature selection

The normalized sparse graph is embedded by the existing sparse diffusion-map
eigensolver. The trivial stationary vector is omitted. Nontrivial eigenpairs
are returned in descending eigenvalue order, and the coordinate convention of
`embed_graph` is retained.

With the default `icm=yes`, ICM selects a prefix of the spectrum with a minimum
rank of one. With `icm=no`, rank selection is bypassed and every nontrivial
eigenpair returned by the scan is retained for target selection, residual-basis
fitting, and representative synthesis. `neigs` must be positive and has no
fixed implementation maximum. The eigensolver limits the requested scan only
to the nontrivial modes available for the selected-particle graph (at most
`N-2` for the current embedding). Thus “all” means all eigenpairs returned
under the configured `neigs` scan limit, not necessarily every possible graph
eigenpair. Zero retained modes is invalid.

The coordinate and spectrum outputs are:

- `flex_diffmap_coordinates.txt`;
- `flex_diffmap_spectrum.txt`;
- `flex_diffmap_graph.txt`.

## 5. Manifold targets and kernel-weighted 3D pre-images

The normalized coordinates written for users and embedding consumers are not
used to place reconstruction targets. `embed_graph` also returns the raw
diffusion coordinates `lambda_q psi_q`; these coordinates define the manifold
metric used for state-descriptor spacing and kernel weighting.

Representative targets are selected from actual particles with the same
`cluster_dmat(...,'kmed',...)` path used by `cls_split`. The distance matrix is
formed from the retained raw diffusion coordinates. The requested count is
controlled by `npreimages` and defaults to eight. It is clamped only when more
representatives than particles are requested.

The medoids are **manifold descriptors only**. They are not interpreted as
single-particle states and are not used as hard reconstruction exemplars.

For each descriptor state `s`, the workflow builds a soft particle-weight
vector with a diffusion-kernel model centered at the medoid coordinate:

```text
w_is ∝ exp(-||x_i - x_medoid(s)||^2 / eps_s) / q_i
```

where `x_i` are raw diffusion coordinates, `eps_s` is a local bandwidth, and
`q_i` is a kernel-density proxy. Weights are normalized per state so
`sum_i w_is = 1`.

To avoid sparse-state collapse, the implementation applies adaptive global
bandwidth inflation: per-state local bandwidth seeds are estimated first, then
a shared scale factor is increased until the minimum state effective sample
size reaches a floor (or an iteration cap is reached). This preserves relative
state locality while preventing extremely peaky weights that produce unstable
"noisy ball" reconstructions.

The 2D-to-3D translation then performs direct weighted reconstruction for each
state (independent state volumes, no residual-basis synthesis):

```text
V_s = argmin_V  sum_i w_is ||A_i V - y_i||^2
```

Operationally, each particle contributes to every state with the same scalar
weight in both the Fourier numerator and CTF/sampling density denominator.
This keeps the inversion mathematically consistent with weighted least-squares
reconstruction while preserving SIMPLE's established reconstruction operators.

This policy is aligned with state-of-the-art diffusion-map pre-image practice:
manifold points are sampled by diverse descriptors (medoids), but state volume
estimation uses soft neighborhoods on the manifold rather than hard exemplar
coefficients.

The weighted reconstruction reads `flex_registered_particles.simple` through the
projected-model matcher batch reader. Its selected `ptcl3D` orientations have
zero in-plane angle and shift, matching the registered images used to train the
diffusion model. Kernel weights remain caller-owned tables and are not written
to project metadata. The original input project is never edited.

There is no alternating latent-coordinate refinement, trajectory-volume fan-out,
or sigma-derived latent scaling. Diffusion coordinates remain fixed after the
graph eigensolve; pre-image states are obtained by one weighted reconstruction
pass from those coordinates.

### Reconstruction identity diagnostic

Before native-image mapping is reconsidered, the registered-project path has a
supported regression diagnostic:

```text
simple_test_exec test=flex_preimage_identity \
  projfile=<flex-run>/flex_registered_particles.simple \
  vol1=<fixed-mean-volume> nspace=<projection-grid-size> nthr=<threads>
```

The diagnostic disables ML regularization, output low-pass filtering, and
output masking. First it prepares every selected registered particle both with
the standard `prep_imgs4rec` path and with the projected-model observation /
forward-transfer path. It requires `observation * transfer` to reproduce the
standard CTF-weighted numerator and requires matching CTF-squared denominators.
It records the result in
`flex_preimage_plane_preparation_metrics_test.txt`. Only after that contract
passes does it run a one-component constant-latent test (`z_i=1`, target `1`):
an ordinary reconstruction and `mean + fitted residual` must agree. It writes
the reference, flex result, difference map, Fourier-shell correlation table,
and scalar metrics as `flex_*_preimage_*_test` outputs. A failure is therefore
localized to plane preparation or to the residual density correction / synthesis;
it is not interpreted as a native-to-registered mapping failure.

The separate basis-conditioning diagnostic is:

```text
simple_test_exec test=flex_preimage_basis_ab \
  projfile=<input-project> vol1=<fixed-mean-volume> \
  nspace=<projection-grid-size> neigs=<rank> nthr=<threads>
```

It runs the graph and registered-frame setup once, then reconstructs identical
medoid targets twice: arm A uses the raw graph eigenfunctions; arm B applies a
full-rank, non-centering whitening to the training eigenfunctions and the
inverse transform to the Nyström targets. The no-centering rule matters because
there is no intercept residual mode. The test records transform/Gram/target
invariants and each raw-versus-canonical state correlation in
`flex_preimage_basis_ab_metrics.txt`, with both volume sets retained for visual
inspection. Because the two arms describe the same states, a material map
difference is a numerical/model-covariance defect, not evidence of a different
representative selection.

## 6. Shared- and distributed-memory execution

The executable uses the same strategy lifecycle as `cls_split` and
`denoise_project`: an abstract `flex_analysis_strategy` exposes `initialize`,
`execute`, `finalize_run`, and `cleanup`; a factory selects its `shmem`,
`worker`, or `master` extension from `part` and `nparts`. The public
`flex_analysis` is the single workflow program name for all three strategy
modes, but distributed workers do use SIMPLE's normal private-executable path.
The user-facing invocation is dispatched by `simple_exec`; `qsys_env` launches
`simple_private_exec prg=flex_analysis ...` for each part. Both executable
dispatch tables invoke `commander_flex_analysis`, and the strategy factory
selects master, worker, or shared-memory execution from `nparts` and `part`.
Thus flex has a private-dispatch registration, but no second worker-specific
program name or reconstruction-only commander such as the removed
`flex_analysis_reconstruct` entrypoint. This is the same single-program-name
factory pattern currently used by `cls_split` and `denoise_project`.

The shared strategy performs registered feature preparation, sparse graph
construction, embedding, ICM selection, representative selection, one modal
projected-residual basis fit, and Nyström state synthesis in the main process
using the established OpenMP kernels. The nano3D latent
caller uses a separate in-process embedding API and deliberately skips
reconstruction; this is not a second executable path.

The distributed path has the same scientific result and uses three worker
barriers behind the same `flex_analysis` program:

1. Feature workers register assigned particles and write registered-particle
   and residual-feature part stacks. Each output stack has one worker-level
   `stack_io` open/close lifetime and is written in assignment order.
2. Graph workers read the residual part stacks through buffered `stack_io`,
   calculate only their assigned source rows of the angularly gated kNN table,
   and write those sparse neighbor rows for master assembly.
3. Residual-model workers receive registered project row indices together with
   their spectral vectors in the standard text assignment files, read registered
   particle images, and accumulate their assigned modal residual numerators and
   one self-density volume per mode with canonical registered-frame
   orientations.

The master constructs the normal `qsys_env` job description and writes one
text particle-assignment file per worker. Feature and graph assignments are
written with `arr2txtfile`; stage-3 assignments additionally carry the spectral
vector needed by that worker. After feature
workers finish, it follows the distributed `denoise_project` convention by
building `flex_registered_particles.simple` over the registered part stacks;
it does not concatenate or rewrite them. After graph workers finish, the
master validates complete, nonoverlapping row coverage, assembles and
normalizes the sparse graph, and owns the global sparse eigensolve, ICM
selection, k-medoid selection, and Nyström target evaluation. It launches
residual-model workers on `flex_registered_particles.simple` with standard text
assignments containing each registered project row index and its spectral
vector. Particle auxiliary metadata is not used as a transport because it is
not persistent project state. There is no custom flex worker-state file.

The workflow-specific implementation is collected in `src/main/flex`; its
`README.md` maps the stages and their dependencies. The thin commander, UI and
executable dispatch registrations remain in their normal architectural layers.
The generic diffusion-graph engines remain in `src/main/pca`, and the generic
reconstructor remains in `src/main/volume`, because those components are shared
with established SIMPLE workflows.

The deliberate scheduling unit is a contiguous range of selected-particle
rows, not class. Feature preparation, graph source-row evaluation, and modal
residual accumulation are particle-separable. The fitted residual basis
and all representative volumes are global outputs; there are no independent
per-class outputs to schedule. The effective worker count is capped by the
selected particle count, and the `qsys_env` ranges are materialized exactly in
the assignment files.

Each worker processes its assigned particles in the established projected
latent-model matcher batches and writes one
`write_mstep_stats_part_file` result. The master calls
`update_basis_from_mstep_stats_part_files`, which validates and sums the modal
right-hand sides plus the diagonal self-density for every mode, then applies
the normal `reconstructor%sampl_dens_correct` path independently to every
residual mode.
The master alone adds the fixed mean and synthesizes the
requested Nyström representatives. Successful reduction removes the temporary
statistics files; standard strategy cleanup removes assignment files and
invokes queue-system cleanup.

Shared and worker M-step batches use the same optimized execution policy. Mean
reprojections are particle-parallel with one projection workspace per OpenMP
thread, following the established projected-model E-step pattern. Mean-only and
fused mean-plus-basis projection use the same normalized three-point Cartesian
Kaiser-Bessel gather and are regression-tested against
`reconstructor%project_fplane`; the corresponding gather/splat weights are also
subject to an adjointness test. Modal 3D insertion keeps one OpenMP team alive
for the complete particle batch and uses the reconstructor's scalarized
three-point Kaiser-Bessel gridding pattern.
The expanded reconstruction lattice retains SIMPLE's nonredundant `h >= 0`
Friedel half. A gather whose rotated first coordinate is negative is evaluated
at the reflected coordinate and complex-conjugated; negative-`h` cells are
interpolation halo only and never receive independent projection samples.
OpenMP insertion colours source lines by at least
`ceil(sqrt(3)*interpolation_width)`, so concurrently updated rotated 3D
interpolation windows cannot overlap.
For every particle, the residual numerator is added once per retained mode
with its `z` coefficient, while that mode's density receives `z^2` times the
CTF-squared contribution. This diagonal regression is required to preserve
the amplitude of L2-normalized diffusion eigenfunctions; a shared density
would suppress every synthesized residual by approximately the particle count.
Shared and distributed paths use the same batched Cartesian Kaiser-Bessel
insertion and the same per-mode density correction.

The former `flex_analysis_mstep` and `flex_analysis_estep` programs and their
iterative PCA state are not part of this policy.

## 7. Output and logging policy

Standard output is stage-oriented and concise. It reports:

- particle count, `k_nn`, `nang_nbrs`, `lp`, and maximum rank;
- registered-feature size and preparation time;
- candidate-count range/mean, graph edge count, and graph time;
- retained rank, whether ICM was enabled, and its convergence summary when used;
- representative count, medoid particle, and hard population;
- cumulative read/preparation, mean-projection, and modal-accumulation times
  after every M-step batch;
- modal projected-residual density correction, finalization, and distributed reduction
  progress when state reconstruction is enabled;
- an explicit reconstruction-skipped message for embedding-only callers;
- registered stack/project names.

When reconstruction is enabled, the requested `outvol` is pre-image state 1.
Further states use the same stem with three-digit suffixes. The default is
`flex_state_001.mrc`. `flex_diffmap_preimages.txt` records state number, medoid
row and particle, and hard population.
Embedding-only execution writes no state volumes.

The nano3D latent trajectory integration is an embedding-only caller. It uses
the rank already selected by `flex_analysis`; it does not run a second ICM rank
selection and does not reconstruct otherwise-unused flex state volumes before
constructing trajectory chunks.

## 8. Implemented versus forward-facing work

| Capability | Status |
| --- | --- |
| Registered Cartesian particle stack and project | Implemented |
| Discrete mean reprojection by `proj` | Implemented |
| Phase-flipped, CTF-matched registered residual | Implemented |
| Angular candidate cap and sparse residual kNN graph | Implemented |
| Ferguson/median kernel bandwidth (`bandwidth_mode`, `bandwidth_tune`) | Implemented |
| Coifman–Lafon α density normalization (`dm_alpha`) | Implemented |
| Sparse diffusion embedding | Implemented |
| ICM prefix rank selection, rank at least one | Implemented |
| `icm=no` retention of every eigenpair returned by the `neigs` scan | Implemented |
| Existing `cluster_dmat` k-medoids in raw diffusion distance | Implemented |
| Denoise-style graph eigenfunction/Nyström coefficients | Implemented |
| Independent modal projection-aware residual 3D basis fit | Implemented |
| Fixed-mean Nyström representative synthesis | Implemented |
| Generative volumes without FSC or low-pass filtering | Implemented |
| Shared-memory execution | Implemented |
| Single shmem/worker/master strategy factory and executable | Implemented |
| qsys assignment files plus project-backed spectral coordinates | Implemented |
| Distributed registered-feature preparation with part stacks | Implemented |
| Distributed angular-gated graph-row calculation and sparse assembly | Implemented |
| Projected-latent distributed sufficient-statistics reduction | Implemented |
| Nano3D embedding-only integration without duplicate ICM/reconstruction | Implemented |
| Iterative PCA/PPCA and sigma estimation | Removed from flex workflow |
| Even/odd maps and FSC per pre-image state | Forward-facing development |
| Posterior covariance/model-variance map | Forward-facing development |
| Project-FSC low-pass filtering of pre-image states | Implemented for weighted-state reconstruction |

A spatial variance map derived from posterior covariance or model variance
would also be valuable, but it requires a clearly specified likelihood and
calibrated uncertainty model. It must not be inferred from diffusion
eigenvalues alone.

Project FSCs are already two-fold cross-validated. The weighted-state
reconstruction path now applies FSC-based low-pass filtering in Fourier space
before writing each state volume, using the original input project's
state-1 FSC (from its `out` segment) for all states. Filtering uses the
existing `fsc2optlp_sub(..., merged=.false.)`
mapping rather than a hard cutoff.

If no compatible FSC metadata/file is available, or if FSC shell count is
incompatible with the reconstruction Nyquist, low-pass filtering is skipped with
an explicit log message and unfiltered volumes are written.
