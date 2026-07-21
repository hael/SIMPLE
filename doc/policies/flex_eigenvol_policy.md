# `flex_eigenvol` Policy: Sparse Diffusion-Map 3D Variability

This policy records the implemented scientific and software contract for
`flex_eigenvol`. Forward-facing developments are separated explicitly from
the production behavior.

## 1. Production model

`flex_eigenvol` is a diffusion-map workflow. It does not fit PCA, probabilistic
PCA, alternating E/M updates, a Gaussian latent prior, an external sigma array,
or a model noise variance.

The implemented input assumptions are:

- one supplied and fixed mean volume (`vol1`);
- one active `ptcl3D` state;
- fixed 3D orientations, in-plane rotations, shifts, and CTF parameters;
- an explicit projection-grid size (`nspace`);
- `ctf=yes` or already phase-flipped (`ctf=flip`) particle images;
- one low-pass limit `lp` used for both graph features and 3D mode
  reconstruction.

Before feature preparation, the implementation requires at least three
selected active particles, matching ptcl2D/ptcl3D field sizes, valid and unique
selected stack slots, and a uniform supported CTF mode. Mixed
`ctf=yes`/`ctf=flip` selected particles are rejected. `nspace` must be supplied
on the command line: the maximum populated `ptcl3D%proj` value cannot recover
the grid size because the final projection directions may be empty.

The standalone workflow returns a continuous low-dimensional coordinate system
and up to 20 corresponding 3D mode volumes. Callers that only consume the
embedding may explicitly disable mode reconstruction. A mode volume is a
linear reconstruction of one selected diffusion coordinate; it is not a PCA
eigenvector and its diffusion eigenvalue is not a particle-image variance.

## 2. Canonical registered frame

Each selected active particle is geometrically registered once by the inverse
of its stored in-plane shift and in-plane angle. Feature preparation includes
the standard particle normalization performed by `transform_ptcls`.
Reconstruction later applies the normal reconstruction batch preprocessing to
the persisted registered image, but it does not repeat the geometric
registration. The implementation extends and uses `transform_ptcls` with an
explicit particle-index batch, following the matcher batch policy
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
the selected particles in feature-row order. The copied project maps those
particles to compact stack indices, marks unselected rows inactive, validates
the resulting ptcl2D and ptcl3D mappings, records every stack as `ctf=flip`,
and stores zero in-plane angle and shift. Defocus, phase shift, and other CTF
metadata remain available. For
astigmatic CTFs, `angast` is rotated into the canonical registered frame by
adding the removed particle `e3` modulo 180 degrees; leaving the raw
astigmatism axis unchanged would be incorrect. An existing `frc3D` record is
preserved when available.

Parameters obtained later in the registered frame must be mapped back to raw
particles using the established `denoise_project`/`map_params_from_den`
composition convention.

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
observation model and uses the same low-pass setting as reconstruction.

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

ICM selects a prefix of the spectrum with a minimum rank of one. The requested
scan limit is clamped to `min(max(neigs,1),20)` and is further limited by the
number of nontrivial graph modes available. Zero selected modes is invalid. No
scientific assumption is made about a preferred ordering beyond the descending
order required by the current prefix-based ICM selector.

The coordinate and spectrum outputs are:

- `flex_diffmap_coordinates.txt`;
- `flex_diffmap_spectrum.txt`;
- `flex_diffmap_graph.txt`.

## 5. Hermitian 3D reconstruction

The Hermitian Fourier reconstruction formulation is required and must be
retained. After subtracting the fixed mean projection, diffusion mode `q`
accumulates particle `i` with

```text
numerator scale = z_iq
density scale   = z_iq**2
```

through `insert_plane_oversamp_multi_scaled`. The numerator is complex and
Hermitian; the density is real and uses the particle CTF-squared sampling
density. Modes are independent, so reconstruction cost and memory are linear
in the selected rank rather than quadratic in the number of modes.

For phase-flipped input, `gen_fplane4rec` applies `abs(CTF)` to the numerator
and `CTF**2` to the density. The resulting maps are CTF-density-corrected
Hermitian least-squares reconstructions. They must not be described as
unit-transfer backprojections or as statistically Wiener-restored maps.

Reconstruction initializes a new builder from
`flex_registered_particles.simple` and reads the preserved registered stack
through the established matcher batch reader. It does not reread or
re-register the raw particle stack. The canonical project has zero in-plane
geometry and canonical astigmatism metadata, while `ctf=flip` causes
`gen_fplane4rec` to use the correct amplitude and density model.

Reconstructors are processed in ordered blocks of at most four modes. The
implementation measures and logs the expanded bytes per mode and the estimated
mode-plus-mean allocation. A smaller final block is used when fewer than four
modes remain.

The configured mask and `lp` are applied during finalization. There is no
post-reconstruction Gram rotation, PCA orthogonalization, trajectory-volume
fan-out, duplicate difference-map output, or sigma-derived scaling.

## 6. Shared- and distributed-memory execution

The executable uses the same strategy lifecycle as `cls_split` and
`denoise_project`: an abstract `flex_eigenvol_strategy` exposes `initialize`,
`execute`, `finalize_run`, and `cleanup`; a factory selects its `shmem`,
`worker`, or `master` extension from `part` and `nparts`. The public
`flex_eigenvol` is the single workflow program name for all three strategy
modes, but distributed workers do use SIMPLE's normal private-executable path.
The user-facing invocation is dispatched by `simple_exec`; `qsys_env` launches
`simple_private_exec prg=flex_eigenvol ...` for each part. Both executable
dispatch tables invoke `commander_flex_eigenvol`, and the strategy factory
selects master, worker, or shared-memory execution from `nparts` and `part`.
Thus flex has a private-dispatch registration, but no second worker-specific
program name or reconstruction-only commander such as the removed
`flex_eigenvol_reconstruct` entrypoint. This is the same single-program-name
factory pattern currently used by `cls_split` and `denoise_project`.

The shared strategy performs registered feature preparation, sparse graph
construction, embedding, ICM selection, and blocked Hermitian reconstruction
in the main process using the established OpenMP kernels. The nano3D latent
caller uses a separate in-process embedding API and deliberately skips
reconstruction; this is not a second executable path.

The distributed path has the same scientific result and uses three worker
barriers behind the same `flex_eigenvol` program:

1. Feature workers register assigned particles and write registered-particle
   and residual-feature part stacks. Each output stack has one worker-level
   `stack_io` open/close lifetime and is written in assignment order.
2. Graph workers read the residual part stacks through buffered `stack_io`,
   calculate only their assigned source rows of the angularly gated kNN table,
   and write those sparse neighbor rows for master assembly.
3. Reconstruction workers reread the registered project, obtain their selected
   coordinates from `ptcl3D` metadata, and accumulate their mode statistics.

The master constructs the normal `qsys_env` job description and writes one
text particle-assignment file per worker with `arr2txtfile`. After feature
workers finish, it follows the distributed `denoise_project` convention by
building `flex_registered_particles.simple` over the registered part stacks;
it does not concatenate or rewrite them. After graph workers finish, the
master validates complete, nonoverlapping row coverage, assembles and
normalizes the sparse graph, and owns only the global sparse eigensolve and ICM
selection. It writes the resulting coordinates into the registered project's
`ptcl3D` metadata before launching reconstruction. There is no custom binary
worker state file.

The deliberate scheduling unit is a contiguous range of selected-particle
rows, not class. Feature preparation and graph source-row evaluation are
particle-separable. Reconstruction produces global per-mode volumes, so every
worker contributes additive Fourier numerator and density statistics to every
mode; there are no independent per-class outputs to schedule. The effective
worker count is capped by the selected particle count, and the `qsys_env`
ranges are materialized exactly in the assignment files.

Each worker processes its assigned particles in the established matcher
batches and modes in blocks of at most four. For every mode and part it calls
`reconstructor%compress_exp`, writes the Fourier numerator with the inherited
image writer, and writes density with `reconstructor%write_rho`. The master
reads those volume/rho pairs into regular reconstructors and combines them
with `reconstructor%sum_reduce`, one mode block at a time. The distributed
handoff and reduction do not serialize, seek through, or manually add
`cmat_exp`/`rho_exp`. Successful
reduction removes the temporary volume/rho pairs and the standard strategy
cleanup removes assignment files and invokes queue-system cleanup.

The former `flex_eigenvol_mstep` and `flex_eigenvol_estep` programs and their
iterative PCA state are not part of this policy.

## 7. Output and logging policy

Standard output is stage-oriented and concise. It reports:

- particle count, `k_nn`, `nang_nbrs`, `lp`, and maximum rank;
- registered-feature size and preparation time;
- candidate-count range/mean, graph edge count, and graph time;
- selected ICM rank and convergence summary;
- reconstruction block size/allocation and progress when mode reconstruction
  is enabled;
- an explicit reconstruction-skipped message for embedding-only callers;
- registered stack/project names.

When reconstruction is enabled, the requested `outvol` is mode 1. Further
selected modes use the same stem with three-digit suffixes. Only selected mode
volumes are written. Embedding-only execution writes no mode volumes.

The nano3D latent trajectory integration is an embedding-only caller. It uses
the rank already selected by `flex_eigenvol`; it does not run a second ICM rank
selection and does not reconstruct otherwise-unused flex mode volumes before
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
| Hermitian CTF-density 3D mode reconstruction | Implemented |
| Shared-memory execution | Implemented |
| Single shmem/worker/master strategy factory and executable | Implemented |
| qsys assignment files plus project-backed worker coordinates | Implemented |
| Distributed registered-feature preparation with part stacks | Implemented |
| Distributed angular-gated graph-row calculation and sparse assembly | Implemented |
| Reconstructor-owned distributed volume/rho I/O and `sum_reduce` | Implemented |
| Bounded four-mode reconstruction and blocked reduction | Implemented |
| Nano3D embedding-only integration without duplicate ICM/reconstruction | Implemented |
| Iterative PCA/PPCA and sigma estimation | Removed from flex workflow |
| Kernel- or neighborhood-weighted state reconstructions along the diffusion manifold | Forward-facing development |
| Posterior covariance/model-variance map | Forward-facing development |
| Project-FSC shrinkage of mode volumes | Forward-facing development |

Kernel- or neighborhood-weighted state reconstructions along the diffusion
manifold are a valuable later development but are outside this first
replacement.

A spatial variance map derived from posterior covariance or model variance
would also be valuable, but it requires a clearly specified likelihood and
calibrated uncertainty model. It must not be inferred from diffusion
eigenvalues alone.

Project FSCs are already two-fold cross-validated. A future shrinkage method
should use those independently validated FSC estimates directly—preferably
mode-specific FSCs when available—to shrink Fourier coefficients continuously.
This is distinct from applying a hard or soft FSC cutoff as an output filter.
