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
`transform_ptcls`. The residual 3D fit separately reads native images through
the projected-model batch path; it does not interpolate them into the
registered 2D frame. The registered project is a row-preserving copy of the
native `ptcl3D` field with its selected in-plane angle and shift set to zero.
Before the native fit, each canonical registered view is composed only with
the corresponding native `ptcl3D` in-plane transform. `ptcl2D` is not read,
modified, or used for this mapping. The implementation hard-fails if the two
projects do not have the same `ptcl3D` row count, if a selected row is inactive,
if the registered row is not canonical, or if their viewing directions
disagree. The implementation extends and
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
- `flex_native_model.simple`, which retains native stack/CTF metadata and
  carries the native-frame `ptcl3D` geometry required for reconstruction;
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

Parameters obtained in the registered frame are mapped back to native images by
restoring only the corresponding source `ptcl3D` in-plane transform.

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

## 5. Manifold targets and residual Nyström 3D pre-images

The normalized coordinates written for users and embedding consumers are not
used to place reconstruction targets. `embed_graph` also returns the raw
diffusion coordinates `lambda_q psi_q`; this preserves diffusion-eigenvalue
scaling when distances are evaluated.

Representative targets are selected from actual particles with the same
`cluster_dmat(...,'kmed',...)` path used by `cls_split`. The distance matrix is
formed from the retained raw diffusion coordinates. The requested count is
controlled by `npreimages` and defaults to eight. It is clamped only when more
representatives than particles are requested; the workflow does not use
silhouette selection or claim that the requested representatives are physical
discrete conformational classes.

The eigensolve also retains the nontrivial graph eigenfunctions `psi_q`. The
training design for the 3D residual model is `z_iq = psi_q(i)`, matching the
residual-mode coefficients used by the `denoise_project` Nyström pre-image.
For every graph node, the corresponding evaluation coefficient is calculated
with the same normalized graph operator:

```text
phi_q(i) = sum_j Wnorm(i,j) psi_q(j) / lambda_q
```

The representative target coefficient is `phi_q` at its medoid node. There is
no Gaussian distance kernel, fitted bandwidth, soft cluster-membership table,
or independent reconstruction of common signal for each representative.

The 2D-to-3D translation uses `simple_flex_projected_latent_model`. It projects the
fixed supplied mean, subtracts that prediction from every prepared particle
Fourier plane, and accumulates the coupled residual normal equations. The
right-hand side is weighted by `z_iq`; the density contains every cross-term
`z_iq*z_ir`. The established coupled block solver therefore accounts for
correlation between retained graph eigenfunctions instead of reconstructing
each component independently. The fitted residual basis volumes retain the
existing spatial-mask policy. As in refine3D restoration, each solved basis is
transformed to real space, normalized for the reconstruction-box scaling,
multiplied once by the inverse 3D Kaiser-Bessel instrument function, and
soft-masked once before it is returned to Fourier space for projection and
state synthesis. Flex explicitly disables the projected model's `lp`-controlled
Fourier-plane truncation and final `bp` operation. The public `lp` parameter
affects graph features only. No FSC-derived weighting, cutoff, or filtering is
applied to the residual basis or synthesized states.

Pre-image state `s` is synthesized as

```text
V_s = V_mean + sum_q phi_q(s) B_q
```

where the original supplied mean is fixed and `B_q` are the fitted residual
3D basis volumes. The basis volumes are an internal representation; the
representative state volumes, not linear one-mode trajectories, are the
diagnostic outputs. Setting every target coefficient to zero reproduces the
supplied mean by construction.

The residual basis fit reads `flex_native_model.simple` through the established
projected-model matcher batch reader. That project is a copy of the native
project, so its images, stack addresses, CTF mode, defocus parameters, and
astigmatism axes remain native. Its selected `ptcl3D` orientations and shifts
are formed by restoring the source `ptcl3D` in-plane transform onto the
canonical registered view. Spectral coordinates remain in the caller-owned
row table (or the distributed worker assignment), rather than project
metadata. This avoids reconstructing interpolated `transform_ptcls` output
while still using the geometry in the frame in which the model was trained. Feeding the
registered stack directly to the kernel would instead normalize and
interpolate those images a second time. The original input project is never
edited.

There is no alternating latent-coordinate refinement, PCA canonicalization,
trajectory-volume fan-out, duplicate difference-map output, or sigma-derived
scaling. Diffusion eigenfunctions remain fixed after the graph eigensolve; the
projected residual basis is fitted once.

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
construction, embedding, ICM selection, representative selection, one coupled
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
3. The master maps registered geometry onto `flex_native_model.simple`.
   Residual-model workers receive native project row indices together with their
   spectral vectors in the standard text assignment files, read native particle
   images, and accumulate their assigned coupled projected-latent M-step
   statistics with the composed native-frame orientations.

The master constructs the normal `qsys_env` job description and writes one
text particle-assignment file per worker. Feature and graph assignments are
written with `arr2txtfile`; stage-3 assignments additionally carry the spectral
vector needed by that worker. After feature
workers finish, it follows the distributed `denoise_project` convention by
building `flex_registered_particles.simple` over the registered part stacks;
it does not concatenate or rewrite them. After graph workers finish, the
master validates complete, nonoverlapping row coverage, assembles and
normalizes the sparse graph, and owns the global sparse eigensolve, ICM
selection, k-medoid selection, and Nyström target evaluation. It maps
registered geometry by `ptcl3D` row into `flex_native_model.simple` and launches
residual-model workers with standard text assignments containing each native
project row index and its spectral vector. Particle auxiliary metadata is not
used as a transport because it is not persistent project state. There is no
custom flex worker-state file.

The workflow-specific implementation is collected in `src/main/flex`; its
`README.md` maps the stages and their dependencies. The thin commander, UI and
executable dispatch registrations remain in their normal architectural layers.
The generic diffusion-graph engines remain in `src/main/pca`, and the generic
reconstructor remains in `src/main/volume`, because those components are shared
with established SIMPLE workflows.

The deliberate scheduling unit is a contiguous range of selected-particle
rows, not class. Feature preparation, graph source-row evaluation, and coupled
normal-equation accumulation are particle-separable. The fitted residual basis
and all representative volumes are global outputs; there are no independent
per-class outputs to schedule. The effective worker count is capped by the
selected particle count, and the `qsys_env` ranges are materialized exactly in
the assignment files.

Each worker processes its assigned particles in the established projected
latent-model matcher batches and writes one
`write_mstep_stats_part_file` result. The master calls
`update_basis_from_mstep_stats_part_files`, which validates and sums those
coupled sufficient statistics, solves the global basis system, and finalizes
the residual basis. The master alone adds the fixed mean and synthesizes the
requested Nyström representatives. Successful reduction removes the temporary
statistics files; standard strategy cleanup removes assignment files and
invokes queue-system cleanup.

Shared and worker M-step batches use the same optimized execution policy. Mean
reprojections are particle-parallel with one projection workspace per OpenMP
thread, following the established projected-model E-step pattern. Mean-only and
fused mean-plus-basis projection use the same normalized three-point Cartesian
Kaiser-Bessel gather and are regression-tested against
`reconstructor%project_fplane`; the corresponding gather/splat weights are also
subject to an adjointness test. Coupled 3D insertion keeps one OpenMP team alive
for the complete particle batch and uses the reconstructor's scalarized
three-point Kaiser-Bessel gridding pattern.
The expanded reconstruction lattice retains SIMPLE's nonredundant `h >= 0`
Friedel half. A gather whose rotated first coordinate is negative is evaluated
at the reflected coordinate and complex-conjugated; negative-`h` cells are
interpolation halo only and never receive independent projection samples.
OpenMP insertion colours source lines by at least
`ceil(sqrt(3)*interpolation_width)`, so concurrently updated rotated 3D
interpolation windows cannot overlap.
Each particle's symmetric `z*z^T` coefficients are packed once in contiguous
pair order; the full packed cross-density vector is updated without dropping
off-diagonal terms. The voxel solve first rejects empty expanded cells from
their diagonal density and right-hand side before materializing the full local
normal matrix. These are execution optimizations only: shared and distributed
paths retain the same complete coupled normal equations.

The former `flex_analysis_mstep` and `flex_analysis_estep` programs and their
iterative PCA state are not part of this policy.

## 7. Output and logging policy

Standard output is stage-oriented and concise. It reports:

- particle count, `k_nn`, `nang_nbrs`, `lp`, and maximum rank;
- registered-feature size and preparation time;
- candidate-count range/mean, graph edge count, and graph time;
- retained rank, whether ICM was enabled, and its convergence summary when used;
- representative count, medoid particle, and hard population;
- cumulative read/preparation, mean-projection, and coupled-accumulation times
  after every M-step batch;
- coupled projected-residual solve, finalization, and distributed reduction
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
| Sparse diffusion embedding | Implemented |
| ICM prefix rank selection, rank at least one | Implemented |
| `icm=no` retention of every eigenpair returned by the `neigs` scan | Implemented |
| Existing `cluster_dmat` k-medoids in raw diffusion distance | Implemented |
| Denoise-style graph eigenfunction/Nyström coefficients | Implemented |
| Coupled projection-aware residual 3D basis fit | Implemented |
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
| Project-FSC shrinkage of pre-image states | Forward-facing development |

A spatial variance map derived from posterior covariance or model variance
would also be valuable, but it requires a clearly specified likelihood and
calibrated uncertainty model. It must not be inferred from diffusion
eigenvalues alone.

Project FSCs are already two-fold cross-validated. Once the pre-image model and
representative generation are otherwise settled, a future shrinkage method can
use those independently validated FSC estimates directly—preferably
state-specific FSCs when available—to shrink Fourier coefficients continuously.
This is distinct from applying a hard or soft FSC cutoff as an output filter.
Until that work is implemented and validated, the generative volumes remain
free of FSC-derived weighting, cutoff, and filtering.
