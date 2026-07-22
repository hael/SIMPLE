# `flex_analysis`: Cartesian diffusion-map replacement

**Status:** implemented initial replacement; validation and performance benchmarking in progress  
**Date:** July 21, 2026  
**Scope:** replace the iterative PCA machinery in the shared- and
distributed-memory `flex_analysis` workflows.  
**Policy:** `doc/policies/flex_analysis_policy.md` describes the implemented
diffusion-map workflow and separates current behavior from later developments.

The initial implementation keeps registered preparation, sparse graph
construction, embedding, and ICM on the master. Its distributed path divides
Hermitian reconstruction statistics by particle interval and reduces them on
the master. Distributed feature/edge preparation remains the next performance
extension, not part of the initial cutover.

## 1. Executive decision

Replace projected EM-PCA with a sparse diffusion-map workflow built from
particles that have first been placed in the in-plane frame supplied by the
fixed 3D alignment.

This refactor must assemble existing SIMPLE components. It must not introduce a
parallel particle-I/O, reprojection, image-registration, batching, or worker
framework. The mature `abinitio2D` and `refine3D` workflows provide the
execution template; the existing Cartesian `projector`, `reproject`,
`transform_ptcls`, matcher particle-I/O, and reconstruction machinery provide
the numerical operations.

PFTC is not part of the scientific representation proposed here. The graph is
defined on Cartesian registered images. Polar code may be consulted for
optimization patterns such as memoization, preallocated thread workspaces, and
batched ownership, but `flex_analysis` must not convert the problem to polar
coordinates.

The essential observation is that SIMPLE's search space is discrete and almost
all expensive geometry can be materialized once:

- roughly 5,000 zero-shift, zero-in-plane-angle reprojections of the supplied
  mean volume, one for each sampled projection direction;
- one registered, band-limited particle record per active particle, roughly
  100,000 records for a representative large data set;
- one symmetry-aware angular-neighborhood table over the approximately 5,000
  projection directions;
- one sparse particle graph with approximately 10 retained neighbors per
  particle.

After those products exist, graph construction, sparse eigendecomposition,
ICM rank selection, and mode reconstruction require no alternating fit and no
repeated projection of evolving eigenvolumes.

The proposed production pipeline is:

```text
fixed mean volume + fixed ptcl3D poses
                    |
                    v
      discrete Cartesian mean-reprojection stack
                    |
raw particles --> one-time in-plane registration --> registered particle stack/features
                    |                                  |
                    +------ CTF-matched residuals <----+
                                      |
                         angular candidate gating
                                      |
                      registered-image k-nearest graph
                                      |
                       sparse diffusion eigensolver
                                      |
                            ICM rank selection
                                      |
             signed Hermitian residual backprojection, O(Q)
                                      |
                   at most 20 reconstructed eigenvolumes
```

This is a deliberate algorithm replacement, not another modification of the
probabilistic or iterative PCA model.

## 2. Scientific contract

### 2.1 Fixed information

The workflow takes the following as fixed:

- the supplied mean volume;
- the active `ptcl3D` particle set;
- each particle's 3D orientation, in-plane angle, and 2D shift;
- the SIMPLE projection-direction grid and point-group symmetry;
- particle CTF metadata;
- the configured mask and low-pass support.

`flex_analysis` does not refine poses, shifts, CTF parameters, or the mean map.

### 2.2 Canonical particle frame

For particle `i`, let its stored pose be decomposed into:

- a projection direction `p_i`;
- an in-plane rotation `e3_i`;
- a 2D shift `s_i`.

The canonical image is obtained by applying `-s_i` and `-e3_i`. The result has
zero shift and zero in-plane rotation and is associated with the discrete
projection direction `p_i`.

This is already the behavior of the `ptcl3D` branch of `transform_ptcls` in
`simple_classaverager_restore.f90`: it groups by `proj`, applies the negative
stored shift, and rotates by the negative stored `e3`. The new workflow must
reuse that convention exactly; it must not create a second interpretation of
the Euler angles or shift sign.

This creates the central implementation invariant:

```text
registered_particle(i) is in the same 2D frame as mean_reprojs(proj_i)
```

Once `registered_particle(i)` has been produced, its geometric preparation is
complete. Comparing it with the model requires only the integer lookup

```text
mean_reprojs(proj_i)
```

There is no shift application, in-plane rotation, Euler-matrix construction,
or volume projection in the comparison loop. Those operations have already
been completed once during particle registration or once during construction
of the approximately 5,000-entry reprojection stack.

### 2.3 Why angularly distant particles are excluded

The graph is not an all-particle similarity graph. Orientation supplies a
physical gate. Only particles whose projection directions lie within an
allowed angular neighborhood are candidates for an image comparison.

The neighborhood may be limited by both:

1. a resolution/diameter-derived angular cutoff; and
2. a computational cap on the number of candidate particles.

A useful starting heuristic is

```text
angular_cutoff_radians ~= target_resolution / particle_diameter
```

The exact constant should be validated rather than embedded as an unexplained
formula. The implementation must report the resulting angular cutoff, the
median candidate count, and the number of particles that reached the cap.

The first implementation uses the following computational cap:

```text
nang_nbrs = 100
```

Here `nang_nbrs` means the maximum number of candidate **particles** considered
for one particle after expanding through nearby projection-direction bins. It
does not mean 1,000 retained graph edges. The final graph default remains:

```text
k_nn = 10
```

For the first release, `nang_nbrs=100` is a hard per-particle candidate cap
applied after the resolution-derived angular cutoff. It is also a normal public
UI parameter so larger or smaller data sets can change the computational bound
without changing code. The implementation must report how often the cap is
reached; experience can later determine whether a different adaptive policy is
worth introducing.

Thus orientation decides which comparisons are meaningful; registered
residual-image distance decides which approximately 10 edges survive.

### 2.4 Residual feature used for graph distances

Let `R_p` be the cached reprojection of the supplied mean map at discrete
direction `p`, with zero in-plane angle and zero shift. Let `T_i` be the
canonical in-plane transform for particle `i`, let `C_i^c` be its real signed
CTF expressed in that canonical frame, and let `S_i^c = sign(C_i^c)`.

The graph feature should represent a CTF-consistent, registered residual:

```text
x_i = graph_preprocess(S_i^c T_i y_i)
    - graph_preprocess(abs(C_i^c) R_p_i)
```

For an input stack already marked `ctf=flip`, the first term is already
phase-flipped and `S_i^c` must not be applied again. For astigmatic CTFs, the
astigmatism axis must rotate with the particle image; merely changing `e3`
while leaving the CTF frame unchanged is incorrect.

Any CTF preparation in this expression belongs to the one-time particle
feature-construction pass. It does not alter the geometric invariant above and
must not introduce a rotation or shift into pairwise graph comparisons.

The correctness-first baseline should mirror the successful
`denoise_project` convention:

- normalize the particle noise using the existing image preprocessing;
- phase-flip consistently;
- apply the established soft mask;
- restrict both particle and mean projection to the same low-pass support;
- compare serialized registered residuals with squared Euclidean distance.

When phase-flipped particles are used, the corresponding mean reprojection
must receive the particle-specific `abs(CTF)` amplitude before subtraction.
The model reprojection itself is cached once; applying the particle transfer
is an inexpensive per-particle multiplication, not a new reprojection.

This phase-flipped representation defines the graph feature and the persistent
registered stack, but it does not imply unit-transfer 3D reconstruction. The
registered project must retain the particle CTF parameters and mark the stack
as already phase-flipped (`ctf=flip`). The reconstruction path must therefore
use `abs(CTF)` as the effective forward transfer and `CTF^2` as its density. It
must never relabel the registered stack as `ctf=no`, apply the signed CTF a
second time, or backproject the phase-flipped residual with unit density.

No external sigma array, fitted noise variance, posterior covariance, or
probabilistic latent prior is part of this definition. The diffusion-kernel
bandwidth is a graph scale and must not be described as a model noise sigma.

## 3. Existing SIMPLE machinery to reuse

### 3.1 Cartesian mean-reprojection model

Materialize the supplied mean map once on the exact discrete orientation grid
used by the particle `proj` labels:

```text
mean_reprojs(p), p = 1 ... nspace
```

Use the established Cartesian path:

- `simple_projector` owns padded Fourier-volume expansion and Cartesian central
  section extraction through `fproject`/`fproject_serial`;
- `simple_simple_volinterp::reproject` already creates a complete stack of
  reprojections for an `oris` object with thread-local padded planes;
- `commander_reproject` demonstrates project/eulspace construction, masking,
  stack writing, and lifecycle cleanup.

The flex workflow should call or narrowly extend those routines. It should not
implement another central-section interpolator or introduce a flex-specific
reprojection file format.

Required invariants are:

- `nspace` and projection ordering match the project `proj` labels exactly;
- the point-group symmetry and reference spiral are the same as those used to
  assign `proj`;
- every reprojection has `e3=0` and zero shift;
- box, sampling, mask, and low-pass support match the registered particles;
- the reprojection stack is generated once from the supplied `vol1`.

For approximately 5,000 projections it is acceptable to retain the Cartesian
reprojection stack in memory at the analysis box, or write/read it using the
normal image-stack machinery if the memory budget requires that. This is an
ordinary SIMPLE image stack, not a new cache abstraction.

### 3.2 Existing particle batching and I/O

Follow the matcher batch policy already used by both `cluster2D` and
`refine3D`:

```fortran
batchsz_max = min(nptcls, params%nthr * BATCHTHRSZ)
nbatches    = ceiling(real(nptcls) / real(batchsz_max))
batches     = split_nobjs_even(nptcls, nbatches)
```

Reuse:

- `prepimgbatch`/`killimgbatch` for builder-owned batch images;
- `discrete_read_imgbatch` and `discrete_read_imgbatch_source` for project-aware
  batched reads, sorted stack runs, bounded file handles, and alternate particle
  sources;
- `alloc_ptcl_imgs` and the matcher pattern of one workspace per thread rather
  than one workspace per particle;
- `build%imgbatch` as the owner of the current raw batch.

The important policy inherited from abinitio2D/refine3D is a single-read batch:
read a particle once, perform all graph preparation that needs the raw image,
and release the batch. Do not add a new whole-particle reader inside the flex
strategy.

### 3.3 Existing Cartesian in-plane registration

Use the `ptcl3D` transform convention already implemented by
`transform_ptcls`:

- obtain the stored `ptcl3D` shift and `e3`;
- apply the negative shift;
- rotate by negative `e3` with the established Cartesian Kaiser-Bessel
  interpolation;
- apply the same normalization, phase-flip, masking, and grid-correction
  conventions as the validated `denoise_project` path.

The required refactor is small: split the current routine into a group-selection
wrapper and an index-driven batch kernel. Existing callers keep the group
wrapper; `flex_analysis` supplies the current batch's particle indices. The
mathematics and interpolation code remain single-source.

The registered Cartesian particle images are a durable workflow product, not
disposable graph scratch. Write them to ordinary project-owned particle stacks
as soon as each batch has been transformed. There is no requirement to keep
100,000 `image` objects resident simultaneously. Reconstruction and later
developments must consume the preserved registered stacks through the same
`discrete_read_imgbatch` pattern; they must not reread and re-register the raw
particles.

The registered project representation must follow the established
`denoise_project`/`map_params_from_den` contract:

- the registered particle image has zero stored in-plane angle and zero stored
  shift;
- the registered stack is recorded as `ctf=flip` when phase flipping was
  applied, while the per-particle defocus, astigmatism, phase-shift, and other
  CTF metadata remain available;
- its projection direction (`proj`, equivalently the fixed out-of-plane part
  of the `ptcl3D` orientation) remains attached to the particle;
- project particle index and source-particle identity are preserved through a
  strict registered-to-raw mapping;
- the original raw project remains the authority for the transform that was
  removed during registration;
- parameters learned later in the registered frame are mapped back by
  composing them with the original raw-particle 2D transform, as
  `map_params_from_den` does with `compose3d2d`, rather than by copying the
  registered-frame angle and shift directly onto the raw particle;
- if a later workflow writes transformed particle images back in the raw
  frame, it applies the inverse of the canonical registration transform.

This preservation and inverse-mapping contract is required even though the
first diffusion-map replacement reconstructs only 3D modes. It allows later
state assignment, local reconstruction, or refinement developments to build on
the registered representation without repeating interpolation or losing the
route back to the raw images.

### 3.4 Reuse first, optimize second

The initial implementation should consist almost entirely of existing calls:

1. build the discrete eulspace through the builder/symmetry machinery;
2. call the Cartesian reprojection machinery once for `vol1`;
3. split active particle indices with the normal matcher batch policy;
4. read each batch through the existing discrete stack I/O;
5. call the index-driven form of `transform_ptcls`;
6. subtract the matching Cartesian mean reprojection and create graph features;
7. write/release the batch.

After step 6, graph construction consumes only:

```text
particle feature
projection-direction index (`proj`)
precomputed angular-neighborhood table indexed by `proj`
```

It must not call `transform_ptcls`, `reproject`, `fproject`, `shift2Dserial`, or
`rotmat2d`. This is the hot-path boundary that makes the comparison stage fast.

Only profiling may justify a new optimized kernel. If needed, borrow the mature
polar implementation patterns—precomputed address maps, reusable thread
workspaces, fused normalization/FFT/shift operations, and explicit memoization—
while retaining Cartesian images and the exact Cartesian metric. An optimization
must be numerically compared against `transform_ptcls`; it must not change the
scientific representation to PFTs.

### 3.5 Working-set estimates

For single-precision real square images, the uncropped storage estimates are:

| analysis box | 100,000 particles | 5,000 reprojections |
| ---: | ---: | ---: |
| 64 | 1.64 GB | 81.9 MB |
| 96 | 3.69 GB | 184 MB |
| 128 | 6.55 GB | 328 MB |
| 300 | 36.0 GB | 1.80 GB |

These figures explain why particle images are processed in the existing bounded
batches and written incrementally to ordinary persistent registered-particle
stacks rather than retained as 100,000 resident `image` objects.
The approximately 5,000 model reprojections are the smaller working set and can
usually stay resident at the analysis box.

Every run must log:

- bytes for the reprojection stack;
- bytes for resident graph features;
- batch size and bytes for registered particle images;
- bytes for graph CSR arrays and diffusion vectors;
- bytes written/read for the registered-particle stack, if present.

## 4. Sparse angular candidate index

### 4.1 Direction-level table

Compute the symmetry-aware nearest projection directions once for the
approximately 5,000-element eulspace. Use the existing point-group-aware
`nearest_proj_neighbors` behavior rather than raw Euler-coordinate distance.

Store the result as a small direction-level adjacency table ordered by angular
distance. This table is independent of particles and can be reused for every
particle assigned to the same `proj` bin.

### 4.2 Particle bins

Build a counting-sort/prefix-sum index:

```text
proj_counts(p)
proj_rowptr(p:p+1)
proj_particles(:)
```

This makes all particles for one projection direction contiguous without
changing their global project indices.

For each owner particle, expand through its ordered direction neighbors until
one of the following occurs:

- the resolution-derived angular cutoff is exceeded;
- `nang_nbrs` candidate particles have been collected;
- all permitted neighboring bins have been exhausted.

The particle itself is excluded. Candidate ordering and cap behavior must be
deterministic. If one projection bin alone contains more than the cap, use a
documented deterministic subsampling policy rather than dependence on worker
or thread order.

### 4.3 Do not materialize `N * nang_nbrs`

For `N=100,000` and `nang_nbrs=1,000`, a dense candidate matrix contains 100
million integers before any distances are stored. That is avoidable.

Graph construction should stream direction-bin pairs and maintain only the
best `k_nn` distances for each particle. Candidate bin pairs should be visited
once where possible, updating both endpoint neighbor lists. This also permits
blocked or BLAS-backed distance evaluation for populated bin pairs.

The production graph builder must therefore be gated from the start. It must
not call the current all-pairs `find_euclidean_neighbors`, which evaluates
every particle pair before retaining a sparse graph.

## 5. Registered-residual k-nearest graph

Add a narrow graph API, conceptually:

```fortran
subroutine build_gated_euclidean_knn_graph(features, proj_index, &
                                           proj_neighbors, nang_nbrs, &
                                           k_nn, graph, diagnostics)
```

The exact signature should use derived types where that clarifies ownership,
but the graph routine must receive all dependencies explicitly.

For each particle:

1. obtain angularly permitted candidate bins;
2. evaluate squared Euclidean distances between registered residual features;
3. retain the nearest `k_nn` candidates;
4. union-symmetrize directed neighbor choices;
5. convert once to `diffmap_graph` CSR;
6. normalize using the existing diffusion-map convention.

The initial kernel weighting should preserve the currently validated scalar
graph convention:

```text
w_ij = exp(-d_ij^2 / epsilon)
```

with `epsilon` derived from retained-neighbor distances. A self-tuning local
kernel can be considered later, but should not be mixed into the first
scientific replacement.

Required graph diagnostics are:

- minimum, median, 95th percentile, and maximum candidate counts;
- minimum, median, and maximum retained degree before and after
  symmetrization;
- total `nnz`;
- number of isolated vertices;
- number and sizes of connected components;
- distance and weight quantiles;
- fraction of edges crossing discrete projection bins;
- angular-distance quantiles of retained edges.

An isolated vertex is a hard construction error unless the workflow explicitly
removes that particle and records why. A graph split into major disconnected
components must be reported before eigendecomposition; it must not silently
produce component-indicator eigenvectors interpreted as conformational modes.

For `N=100,000` and `k_nn=10`, the symmetrized graph should be on the order of
one to two million directed CSR entries, not `N^2`. Integer columns plus weight
and normalized-weight arrays are only on the order of a few tens of megabytes.

## 6. Diffusion embedding and ICM selection

Reuse `diffmap_graph`, `embed_graph`, `sparse_eigh`, and the existing ICM
spectral-rank selection. Do not create a second diffusion eigensolver inside
the flex strategy.

The eigensolver must:

- consume the sparse normalized graph;
- exclude the trivial stationary eigenvector;
- scan enough of the spectrum for ICM to distinguish retained modes from the
  tail;
- reconstruct no more than 20 modes;
- present nontrivial modes to the existing ICM selector in descending
  diffusion-eigenvalue order;
- impose a deterministic sign based on the largest-magnitude coordinate, with
  project particle index as the tie breaker.

The descending order is an internal precondition of the existing ICM
implementation: it uses adjacent eigengaps and returns a retained prefix. It is
not a claim that mode 1 has the greatest PCA variance or is biologically more
important than mode 2. Outside satisfying ICM and providing stable file
indices, the workflow assigns no scientific significance to mode order.

The public `neigs` parameter should retain a simple meaning:

```text
neigs = maximum number of reconstructed diffusion modes, clamped to 1..20
```

The internal spectral scan may exceed `neigs` so that ICM sees a noise tail.
The scan limit should be a named constant or a separately reviewed typed
parameter, not another interpretation of `neigs`.

For retained mode `q`, standardize the coordinate over active particles:

```text
z_iq = (phi_iq - mean(phi_q)) / sdev(phi_q)
```

The diffusion eigenvalue remains spectral evidence used by ICM; it is not an
explained-variance fraction and must not be labeled as one in the output table.

ICM is the denoising/model-selection step. The flex call must set a minimum
retained rank of one: returning zero modes is not a meaningful public
`flex_analysis` result. A homogeneous control may produce a weak diagnostic
mode and a spectrum without a convincing separation, but the command still
emits at least one mode. No posterior variance or PPCA noise term is required.

## 7. Hermitian 3D mode reconstruction

Removing PCA does not remove the Hermitian Fourier formulation.

For canonical particle `i`, let:

- `r_i^pf` be its phase-flipped canonical Fourier residual after subtracting
  the cached mean reprojection multiplied by `abs(C_i^c)`;
- `C_i^c` be its real signed CTF in the canonical frame;
- `H_i^pf = sign(C_i^c) C_i^c = abs(C_i^c)` be the effective transfer of the
  phase-flipped observation;
- `P_p_i` insert/backproject the zero-in-plane-angle plane for projection bin
  `p_i`;
- `z_iq` be standardized diffusion coordinate `q`.

Reconstruct mode `q` from:

```text
numerator_q = sum_i z_iq * P_p_i^H * conjugate(H_i^pf) * r_i^pf

density_q   = sum_i z_iq^2 * P_p_i^H * |H_i^pf|^2 * P_p_i

B_q         = restore(numerator_q, density_q)
```

For SIMPLE's real signed CTF, this is algebraically identical to constructing
the normal equations from the unflipped residual
`r_i = y_i - C_i^c R_p_i`, because
`r_i^pf = sign(C_i^c) r_i`:

```text
conjugate(H_i^pf) * r_i^pf
    = abs(C_i^c) * sign(C_i^c) * r_i
    = C_i^c * r_i

|H_i^pf|^2 = |C_i^c|^2
```

Phase flipping therefore does not prevent CTF-amplitude restoration. The
adjoint still uses the conjugate effective transfer and the density still uses
its squared magnitude. This is the Hermitian normal-equation convention and
must be tested with an adjoint check.

Use the established `gen_fplane4rec` handling rather than reimplementing this
logic. Its `CTFFLAG_FLIP` branch already multiplies the phase-flipped Fourier
signal by `abs(CTF)` for the numerator and stores `CTF^2` for the denominator.
This also prevents a duplicate raw registered stack: the preserved
phase-flipped registered images are sufficient for both graph construction and
proper CTF-density reconstruction when their project metadata are correct.

The resulting mode maps are CTF-density-corrected Hermitian least-squares
reconstructions. They should not be called strictly Wiener-restored maps unless
a defined signal/noise or FSC reliability term is also applied. The first
replacement deliberately has no external sigma model; project-FSC-driven,
already cross-validated shrinkage remains the forward-facing regularization
development described in Section 17.

There are only `Q` numerators and `Q` densities. There are no `Q(Q+1)/2`
cross-density volumes and no alternating E/M loop. ICM runs before this phase,
so only retained modes are reconstructed.

Use the existing reconstructor, gridding, masking, and low-pass restoration
machinery. Do not reproduce interpolation or Fourier indexing inside the
diffusion module.

### 7.1 Coordinate and map semantics

Diffusion eigenvectors are orthogonal over graph vertices, but their separately
reconstructed 3D maps need not be exactly orthogonal under a 3D Fourier metric.

The first implementation does not reorder the reconstructed maps and does not
apply a post-reconstruction Gram rotation or orthogonalization. Each output map
remains paired one-to-one with the standardized diffusion coordinate from
which it was reconstructed. This preserves the learned continuous state
coordinates as the low-dimensional particle features rather than replacing
them with a second, Fourier-map-derived coordinate system.

Map norms and the small Hermitian Fourier Gram matrix may be reported as
diagnostics. Non-orthogonality of the restored 3D maps is not itself a failure.
A nonfinite or numerically empty reconstruction is a hard diagnostic problem
and must not be hidden by rotating it into other modes.

### 7.2 Linear interpretation

Because each retained coordinate is standardized, a user-requested linear
diagnostic trajectory has a direct scale:

```text
mean + a * B_q,  a = -3, -2, -1, 0, 1, 2, 3
```

Such maps are linear visualizations of a nonlinear embedding direction. They
are not written by default and must not be described as a complete
parameterization of the manifold. Kernel- or neighborhood-weighted state
reconstructions along the diffusion manifold are a valuable later development
but are outside this first replacement.

## 8. Shared-memory execution

The shared-memory workflow should have explicit, noniterative phases:

1. validate inputs and select active particles;
2. construct the exact discrete eulspace and symmetry-aware direction table;
3. generate the Cartesian mean-reprojection stack once;
4. use the matcher batch policy and existing discrete I/O to register each
   particle once;
5. form residual graph features;
6. build the gated sparse kNN graph;
7. run the sparse diffusion eigensolver;
8. select the retained rank with ICM;
9. reconstruct retained 3D modes in memory-budgeted mode blocks;
10. write the paired maps, coordinates, project outputs, and diagnostics.

The particle loop should mirror `refine3D_exec`/`cluster2D_exec`: derive
`batchsz_max` from `nthr * BATCHTHRSZ`, allocate builder/thread workspaces once,
iterate `split_nobjs_even` batches, and clean the batch toolbox once at the end.
Particle preprocessing should be parallel over the current batch. Graph distance work
should be parallel over owner bins or rows with thread-local top-`k_nn` state.
No two threads should update one mutable top-k list without explicit ownership.

Mode reconstruction can process all retained modes in one particle pass if
memory permits. Otherwise use a small mode block selected from an explicit
memory budget. The run must log the block size and expected reconstructor
bytes. It must never silently allocate one expanded reconstruction object per
mode per thread.

## 9. Distributed-memory execution

The distributed path must implement the same graph and reconstruction model as
the shared path. Distribution changes ownership, not scientific semantics.

### 9.1 Stage A: registered-particle/feature parts

Use the normal SIMPLE distributed particle partition and `qsys_env` scheduling
pattern. Each worker owns its `fromp:top` range, applies the same local matcher
batch loop used in shared memory, and writes ordinary part-local products:

- registered Cartesian particle stack part, retained for reconstruction and
  later registered-frame workflows;
- a compact graph-feature part;
- particle-index and projection-bin mapping using the established part-file
  header/version pattern.

Workers use `discrete_read_imgbatch`, `build%imgbatch`, thread-local transform
workspaces, and standard stack I/O. They do not own a new cache service. The
master generates the common Cartesian mean-reprojection stack once; workers
read it as a normal reference stack or receive the same in-process product when
running locally.

### 9.2 Stage B: graph edge parts

Assign each graph row to exactly one worker. That worker reads only feature
parts needed by the angular bins for its rows, computes its nearest neighbors,
and writes a versioned edge part keyed by global graph row.

The master must validate:

- every expected row appears exactly once;
- no row is duplicated or omitted;
- all particle indices and projection bins match the feature-part manifest;
- every column index is in range;
- edge distances are finite;
- part headers are mutually compatible.

The master then symmetrizes, normalizes, and constructs the global CSR graph.
The sparse eigensolve remains a shared-memory operation on the assembled graph;
for approximately 100,000 vertices and roughly one to two million entries this
is the simpler and expected design.

### 9.3 Stage C: reconstruction parts

After ICM rank selection, distribute particle subsets again. Each worker
accumulates Hermitian numerator and density parts for a bounded block of
retained modes and writes versioned reconstruction statistics. The master
reduces and restores one mode block at a time.

This stage replaces both current `flex_analysis_mstep` and
`flex_analysis_estep` workers. Proposed private stage names should describe the
new responsibilities, for example:

```text
flex_analysis_prepare
flex_analysis_graph
flex_analysis_reconstruct
```

The final names should follow existing private-program conventions. The old
M-step/E-step programs, private-program registrations, APIs, state files, and
temporary basis-volume protocol must be removed after the new distributed path
passes equivalence tests. Do not leave dead PCA workers as an undocumented
fallback.

### 9.4 Artifact lifecycle

Any new compact feature/edge/statistics part format requires:

- magic and format version;
- stage and part number;
- particle count and expected global row interval;
- box, sampling, Fourier limits, `nspace`, symmetry, and mode count;
- source project, registered stack, and mean-reprojection identity;
- write-to-temporary followed by atomic rename.

Successful stages remove obsolete temporary artifacts. Failed stages preserve
them and print their paths for diagnosis.
The registered particle stacks, their project records, and the
registered-to-raw map are successful-run products and are never classified as
temporary stage artifacts.

## 10. Proposed module boundaries and reuse map

### `simple_classaverager` / `simple_classaverager_restore`

Extend the registration utility without changing existing callers:

- retain the current group-based `transform_ptcls` wrapper;
- add an arbitrary-index registration entry point;
- make both delegate to one implementation of shift, `e3`, interpolation, CTF,
  and gridding conventions.

This layer owns particle-image registration, not diffusion policy.

### `simple_projector` and `simple_simple_volinterp`

Continue to own Cartesian volume projection. Flex may add a narrow batched or
preallocated entry point only if the existing `reproject` allocation pattern is
measured as limiting. Central-section interpolation remains here and is not
copied into a flex module.

### `simple_matcher_ptcl_io` and `simple_matcher_ptcl_batch`

Continue to own project-aware particle reads, builder batch images, alternate
particle sources, batch allocation, and thread-local image workspaces. Flex
should reuse these public helpers or extract a genuinely common Cartesian batch
helper. It must not create its own stack-reader implementation.

### Small new `simple_flex_diffmap_features` under `src/main/pca`

Owns only flex-specific operations that do not already exist:

- mapping registered particles to their cached Cartesian mean reprojection;
- forming the graph residual feature;
- serializing/deserializing compact feature parts;
- feature diagnostics.

It does not own particle I/O, image registration, projection, distributed
scheduling, or reconstruction.

### `simple_diff_map_graphs`

Owns:

- gated Euclidean kNN construction;
- top-k retention;
- CSR symmetrization and normalization;
- graph structural diagnostics.

Keep the existing all-pairs graph builder for small class workflows. Do not
silently change `cls_split` or `denoise_project` behavior as part of this
refactor.

### `simple_diffusion_maps` and `simple_diff_map_denoise`

Continue to own sparse embedding and ICM spectral selection. Add only narrow
reusable interfaces required to avoid recomputing the same eigensystem twice.
In particular, the flex path should obtain coordinates, eigenvalues, and the
selected rank from one sparse solve rather than calling a rank helper that
discards coordinates and then solving again.

### New `simple_flex_diffmap_rec3D` under `src/main/pca`

Owns:

- diffusion-coordinate standardization;
- Hermitian signed residual accumulation;
- mode-block sufficient-statistics I/O;
- final map norms and Hermitian Gram diagnostics;
- deterministic mode signs and output spectral metadata.

It uses reconstruction domain objects but does not own workflow scheduling.
The `rec3D` suffix is intentional: this module converts selected diffusion
coordinates into Hermitian 3D mode-reconstruction statistics and restored
volumes. It does not own the generic diffusion embedding or imply that the
graph machinery itself is specific to 3D data.

### `simple_flex_analysis_strategy`

Remains the orchestration boundary for the public command. It sequences
Cartesian reprojection, matcher-style batch registration, graph construction,
embedding, ICM, reconstruction, and output stages and selects shared or
distributed execution. It must not contain low-level rotation, graph, or
gridding kernels.

The public routine should be renamed from `run_flex_analysis_linear` to a
model-neutral name such as `run_flex_analysis`, with all in-process callers
updated.

### Result and trajectory consumers

The implementation uses the model-neutral `flex_embedding_result`, containing:

- project particle indices;
- retained standardized coordinates;
- diffusion eigenvalues;
- a method identifier naming diffusion maps.

`simple_trajectory_chunker` and the `trajectory_make_projavgs` caller consume
this model-neutral result. A diffusion eigenvalue is not labeled as explained
variance.

### Commanders and private dispatch

`simple_commanders_flex_analysis` remains orchestration-only. Replace M/E worker
commanders with prepare/graph/reconstruction worker commanders, then update:

- `simple_private_exec_api.f90`;
- `production/simple_private_exec_driver.f90`;
- `src/main/exec/simple_exec_denoise.f90`;
- `src/main/ui/simple/simple_ui_denoise.f90`;
- private-program tests and generated dispatch surfaces as required.

## 11. Parameter policy

Proposed public parameters are:

| parameter | proposed meaning | default |
| --- | --- | ---: |
| `neigs` | maximum diffusion modes requested from the graph eigensolve | 20, with no workflow-specific hard maximum |
| `k_nn` | retained particle neighbors in sparse graph | 10 |
| `nang_nbrs` | maximum angularly gated candidate particles per particle | 100 |
| `lp` | graph-feature low-pass limit; it does not filter generated volumes | existing flex default |
| `mskdiam` | particle/model support used during feature construction | existing project convention |
| `nparts`, `nthr` | execution controls only | existing conventions |

`nang_nbrs` must be added to `simple_parameters.f90`, initialized through the
normal parameter defaults, parsed through the typed parameter path, and
registered as an input of `new_flex_analysis` in
`src/main/ui/simple/simple_ui_denoise.f90`. Its UI default is 100. It
must not be read ad hoc from `cmdline` after parameter initialization.

The first implementation uses the same `lp` for graph features and Hermitian
mode reconstruction. It does not introduce a second graph-specific low-pass
parameter.

The first implementation should avoid exposing analysis box, feature dimension,
kernel epsilon, eigensolver tolerance, or mode block size as casual user knobs.
Derive them from physical support and memory budgets, log the derivation, and
promote them to typed parameters only if experiments demonstrate a real need.

`maxits` has no meaning in the diffusion workflow and should be removed from
the public UI when the replacement becomes production. It must not remain as a
no-op compatibility parameter.

## 12. Output contract

No backward compatibility with the PCA output set is required. The replacement
should produce a small, explicit output set derived from one `outvol` stem:

- `<prefix>_mode001.mrc` through `<prefix>_modeQQQ.mrc`, exactly one final 3D
  map for each ICM-retained diffusion coordinate;
- `<prefix>_coordinates.txt`, one row per active particle;
- `<prefix>_spectrum.txt`, one row per scanned nontrivial diffusion mode;
- `<prefix>_graph.txt`, one compact graph/connectivity summary;
- the normal execution-directory SIMPLE project, updated to reference the
  preserved registered particle stack or distributed stack parts and the
  registered-to-raw particle mapping.

The supplied mean volume is not copied. The command does not write automatic
`mean + a * mode` trajectory volumes, duplicate difference volumes, a second
alias for mode 1, PCA residual volumes, or temporary reprojection volumes as
successful-run products. A user can form a linear diagnostic trajectory from
the mean and a mode when needed without multiplying the default volume output
count by seven.

The spectral table contains:

```text
MODE  DIFFUSION_EIGENVALUE  ICM_RETAINED  COORDINATE_SD  MAP_FOURIER_NORM
```

Do not report diffusion eigenvalues as PCA variance fractions or cumulative
explained variance.

The coordinate table should include at least:

```text
PROJECT_PARTICLE_INDEX  PROJECTION_BIN  Z1 ... ZQ
```

Optional diagnostics may include graph degree, distance to the `k_nn`th
neighbor, and component identifier. These are preferable to carrying obsolete
PCA residual-energy fields whose semantics no longer apply.

### 12.1 Standard output

Only the master process writes normal progress to STDOUT. Use one short startup
summary, one completion line per major phase, one retained-mode table, and one
final output summary. The intended shape is:

```text
FLEX_ANALYSIS DIFFUSION MAP
particles=... projections=... lp=... k_nn=10 nang_nbrs=100 max_modes=20
[1/6] reprojections and particle registration ... done (... s)
[2/6] residual features                    ... done (... s)
[3/6] sparse graph                         ... done (... s)
      candidates=... nnz=... degree[min/median/max]=... components=...
[4/6] diffusion embedding and ICM          ... done (... s)
      scanned=... retained=...
[5/6] Hermitian 3D reconstruction          ... done (... s)
      input=phase-flipped transfer=abs(CTF) density=CTF^2
      mode  eigenvalue  map
         1  ...         <prefix>_mode001.mrc
[6/6] outputs                              ... done (... s)
      project=... coordinates=... spectrum=... graph=...
```

Do not print per-particle, per-candidate, or per-edge messages. Distributed
workers write at most one concise part-completion summary during normal
operation; detailed worker timings belong in benchmark/diagnostic files.
Warnings must identify an actionable condition such as a capped candidate
list, isolated vertex, disconnected graph, failed eigensolve, or invalid mode
reconstruction. Internal temporary filenames are printed only on failure when
they are intentionally preserved for diagnosis.

## 13. Work explicitly removed

Once the replacement is validated, remove:

- iterative latent initialization and orthonormalization;
- alternating M- and E-step loops;
- coupled eigenvolume M-step statistics;
- `Q(Q+1)/2` cross-density arrays and files;
- temporary basis maps passed from M-step to E-step workers;
- E-step per-particle Hermitian solves;
- PCA canonicalization based on the observation metric and latent covariance;
- `flex_analysis_mstep` and `flex_analysis_estep` private programs;
- PCA-specific tests that have no remaining consumer.

`simple_flex_projected_latent_model.f90` is currently consumed only by
`flex_analysis` and its dedicated test. It can be removed from production once
the new workflow and tests replace those references. General Hermitian
reconstruction helpers that are useful to the new path should be moved to the
nearest reconstruction or flex-diffusion module rather than retaining the
entire obsolete model.

## 14. Refactoring sequence

### Phase 0: freeze evidence and benchmark

- preserve the failing current output statistics and representative volume
  comparisons;
- record wall time, peak RSS, worker count, box, sampling, and particle count;
- define one small synthetic data set and one representative PfCRT-scale run;
- do not optimize before phase timings exist.

### Phase 1: registration and reprojection reference path

- split `transform_ptcls` into its existing group wrapper and an index-driven
  Cartesian batch kernel;
- materialize the exact discrete Cartesian mean-reprojection stack through
  `reproject`/`projector`;
- drive particle reads with `prepimgbatch`, `discrete_read_imgbatch`,
  `build%imgbatch`, and the matcher batch sizing policy;
- produce CTF-consistent registered residuals in real space;
- write persistent registered particle stacks and the registered-to-raw project
  mapping while each batch is resident;
- compare them against direct on-demand transforms and projections;
- add unit tests for shift sign, `e3` sign, symmetry, and astigmatic CTF frame.

### Phase 2: gated shared-memory graph

- build direction bins and the symmetry-aware direction-neighbor table;
- implement streamed gated Euclidean top-k selection;
- construct and validate CSR without an all-pairs particle loop;
- compare against brute force on small data using the same angular gate;
- add graph diagnostics and connectivity checks.

### Phase 3: diffusion embedding and ICM

- run one sparse eigensolve;
- sort nontrivial eigenpairs into the descending spectral order required by the
  existing ICM prefix selector;
- select at least one retained coordinate with ICM;
- standardize and sign-fix coordinates;
- write graph, spectrum, and coordinate diagnostics before reconstructing any
  volume;
- validate repeatability across thread counts.

### Phase 4: Hermitian mode reconstruction

- implement one-mode signed numerator/density accumulation first;
- add the adjoint and synthetic one-mode tests;
- extend to at most 20 modes with memory-budgeted blocks;
- retain the one-to-one pairing between each diffusion coordinate and its 3D
  map without a post-reconstruction Gram rotation;
- validate map scale, signs, Gram diagnostics, and shared-memory repeatability.

### Phase 5: measured Cartesian optimization

- profile reprojection, particle registration, stack I/O, and feature
  serialization independently;
- first reuse existing fused image operations and memoized address maps;
- if a new kernel is still required, borrow workspace/memoization templates
  from mature polar code without introducing a PFTC dependency;
- compare registered images, residuals, distances, edges, and final maps against
  the phase-1 Cartesian reference before enabling the optimized path.

### Phase 6: distributed execution

- use normal `qsys_env` particle partitions and matcher-style local batch loops;
- write ordinary registered-stack/feature parts with strict versioned metadata;
- add owner-row graph edge workers and strict reduction checks;
- add mode-block reconstruction workers and reduction;
- compare graph edges, diffusion subspaces, coordinates, and maps against the
  shared path;
- add cleanup and failed-artifact preservation policy.

### Phase 7: cutover and cleanup

- switch the public `flex_analysis` command to diffusion maps;
- remove the public `maxits` input;
- remove PCA M/E workers, formats, APIs, private-program entries, and tests;
- update in-process trajectory consumers to the model-neutral result;
- rewrite `flex_analysis_policy.md` to describe the implemented diffusion
  contract and clearly retain forward-facing variance-map/FSC developments;
- regenerate any affected documentation or dispatch products.

The implementation branch may temporarily contain both paths for comparison,
but production should not expose a permanent `pca_mode` switch for this
command. The purpose of the refactor is to replace the failing machinery, not
maintain two scientific workflows indefinitely.

## 15. Validation matrix

### 15.1 Registration and reprojection tests

- zero shift and zero `e3` are identity operations;
- known shifts and in-plane rotations are removed with the existing sign
  convention;
- arbitrary-index and group-based `transform_ptcls` agree;
- precomputed and on-demand model reprojections agree;
- raw `ctf=yes` and phase-flipped `ctf=flip` preparation of the same synthetic
  observations produce equivalent Hermitian numerators, `CTF^2` densities, and
  restored maps within tolerance;
- the registered project retains particle CTF metadata and advertises
  phase-flipped stacks as `ctf=flip`, never `ctf=no`;
- C1 and nontrivial point-group direction lookup agree with existing symmetry
  helpers;
- astigmatic CTF orientation is correct after canonical rotation;
- any optimized Cartesian transform agrees with `transform_ptcls` within the
  defined tolerance;
- part headers reject altered model, box, sampling, low-pass, symmetry,
  eulspace, or particle assignment.

### 15.2 Graph tests

- gated graph equals a brute-force gated reference for a small set;
- final neighbor lists contain no self edges, invalid indices, duplicates, or
  nonfinite distances;
- symmetrization is deterministic;
- `k_nn` and `nang_nbrs` limits are respected;
- results are stable across thread counts and distributed row partitions;
- graph connectivity and component diagnostics are correct on constructed
  examples.

### 15.3 Diffusion tests

- known clustered and continuous synthetic manifolds recover expected
  nontrivial coordinates up to sign or degenerate-subspace rotation;
- the stationary eigenvector is never emitted as a flex coordinate;
- eigenpairs are supplied to ICM in descending diffusion-eigenvalue order;
- ICM retains the same rank when the same graph is supplied through shared and
  distributed assembly;
- at least one coordinate is retained;
- coordinate signs and file indices are deterministic outside degenerate
  eigenspaces, without assigning biological importance to their order.

### 15.4 Reconstruction tests

- projection/backprojection adjoint test with complex Fourier data;
- homogeneous particles produce at least one diagnostic mode but no falsely
  well-separated diffusion spectrum or strong structural interpretation;
- a simulated one-coordinate deformation recovers the mode and coordinates up
  to the documented normalization/sign convention;
- two independent simulated modes recover the correct subspace;
- nonzero shifts, nonzero `e3`, varied defocus, astigmatism, and symmetry are
  included;
- each reconstructed map remains paired with the diffusion coordinate used in
  its Hermitian accumulation;
- shared and distributed numerator/density parts reduce to the same maps.

### 15.5 Workflow tests

- `nparts=1` and `nparts>1` produce equivalent graph spectra and mode
  subspaces;
- raw and supported alternate particle sources follow their established I/O
  conventions;
- the registered particle stacks and registered-to-raw map persist after a
  successful run;
- registered-frame test parameters composed back through the
  `map_params_from_den` convention recover the corresponding raw-frame
  transforms;
- trajectory chunking consumes diffusion coordinates without PCA variance
  assumptions;
- output maps are visibly distinct when the data contain heterogeneity;
- repeated runs with the same inputs are deterministic within documented
  eigenspace degeneracy.

## 16. Performance measurements and acceptance criteria

Measure every phase separately:

```text
mean reprojections
particle read/normalize/register/write
residual feature construction
angular candidate traversal
distance/top-k selection
CSR assembly/normalization
sparse eigensolve
ICM selection
Hermitian reconstruction per retained mode
output
```

Record wall time, CPU time, peak RSS, bytes read/written, candidate comparisons,
and graph `nnz` for shared and distributed runs.

The replacement is accepted only when:

1. every particle is read and registered once during particle preparation;
2. every discrete Cartesian mean reprojection is generated once;
3. comparison/graph construction performs no particle rotation, particle
   shift, or model reprojection and uses `proj` only for reference lookup and
   angular candidate gating;
4. graph construction evaluates only angularly permitted candidate pairs and
   contains no `N^2` particle loop;
5. the final graph has approximately `k_nn` local support and no unexplained
   isolated particles;
6. only the ICM-retained modes, at most 20, are reconstructed;
7. reconstruction storage is `O(Q)` and contains no coupled `Q^2` density
   volumes;
8. no alternating PCA iterations or per-particle latent solves remain;
9. the Hermitian adjoint and shared/distributed equivalence tests pass;
10. synthetic heterogeneous data produce distinct, recoverable modes; a
    homogeneous control still emits at least one diagnostic mode but does not
    show a misleadingly separated diffusion spectrum;
11. the representative PfCRT-scale run is substantially faster and uses less
    peak memory than the current iterative PCA workflow;
12. output semantics and policy documentation no longer describe diffusion
    eigenvalues as PCA explained variance.

No fixed speedup factor should be promised before the phase-0 benchmark is
captured. The structural expectation is nevertheless strong: the expensive
registration and reprojection work becomes one-time, graph storage is sparse,
the eigensolver operates on CSR, and 3D reconstruction is performed once for
only the selected modes.

## 17. Forward-facing items, not part of this refactor

The following remain worthwhile but must not delay or silently enter the first
replacement:

- posterior covariance or probabilistic uncertainty estimates;
- a spatial model-variance or predictive-variance map;
- project-FSC-based, already two-fold cross-validated shrinkage/filtering of
  reconstructed modes;
- local/self-tuning diffusion kernels;
- kernel-weighted nonlinear state maps along the diffusion manifold;
- joint pose or shift refinement;
- multi-state mean models;
- accelerator-specific graph or reconstruction kernels.

The existing policy's forward-facing variance-map note should be retained when
the policy is rewritten. A deterministic structural variance map derived from
retained modes must be distinguished from posterior uncertainty.

## 18. Resolved implementation decisions

1. `nang_nbrs` is a public typed UI parameter with an initial default of 100.
   For the first implementation it is a hard per-particle candidate cap applied
   after the resolution-derived angular cutoff. Cap-hit statistics are
   reported so this policy can be revisited with evidence.
2. The initial graph metric is the phase-flipped, CTF-matched registered
   residual described in Section 2.4.
3. Graph feature construction and 3D mode reconstruction use the same `lp`.
   There is no second graph-specific low-pass parameter.
4. Registered Cartesian particle stacks are preserved as project products and
   reused directly. The project retains a strict registered-to-raw mapping and
   later registered-frame parameters are composed back to the raw frame using
   the established `denoise_project`/`map_params_from_den` convention.
5. The first release does not perform ordered map orthogonalization or a
   compensating coordinate rotation. It preserves the diffusion coordinates as
   the continuous low-dimensional state features and reconstructs one paired 3D
   map per retained coordinate. Descending eigenvalue order is maintained only
   because the existing ICM selector assumes an ordered spectrum and returns a
   retained prefix.
6. Zero retained modes is not a valid public result. Flex invokes ICM with a
   minimum rank of one and always emits at least one diagnostic mode.
7. The sparse eigensolve is centralized on the master after distributed graph
   assembly for the initial approximately 100,000-particle target. Workers
   distribute particle preparation, graph-row comparison, and reconstruction
   statistics; the master owns the assembled sparse CSR graph and its
   eigensolve. This avoids adding a distributed eigensolver while the expected
   one-to-two-million-entry graph fits comfortably in shared memory. Revisit
   only if measured graph size or eigensolve memory demonstrates a need.
8. No PCA-output filename or column compatibility is required. Section 12 is
   the authoritative minimal output and STDOUT contract; obsolete trajectory,
   difference, residual, and explained-variance products are removed rather
   than renamed.

The review questions are therefore resolved. Implementation starts with the
phase-0 benchmark and validation data set so the replacement has a measured
correctness and performance baseline.
