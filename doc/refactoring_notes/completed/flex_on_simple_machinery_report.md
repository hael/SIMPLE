# FLEX on SIMPLE machinery: report

Date: 2026-10-08. Written for Afan, who developed FLEX, and for anyone who works on
`flex_pca` next. It describes what changed between master `92d684c6f` and the working tree of
2026-10-08, why, where every piece of FLEX now lives, what SIMPLE gained, what was measured and
what is still to be tested. It reads on its own; the plan it implements is
`doc/refactoring_notes/completed/flex_on_simple_machinery.md` (the plan), for anyone who wants the
day-by-day record.

## 1. Summary

FLEX now runs on SIMPLE's machinery. It no longer reconstructs maps, stores weights, writes its
own part files for the state stage, publishes a special volume kind or carries clustering code
of its own. What stays in `src/main/flex` is the method: the probe fit with its E-step and
coupled M-step, the posterior, the latent deconvolution's noise model, the state placement
decisions, the merge gates, UMAP and the readouts.

SIMPLE gained four capabilities on the way:

1. **A state weight set.** Per-particle, per-state weights in [0, 1], published as one file per
   state plus a manifest, registered in the project and validated against it.
2. **One reconstruction service.** `reconstruct3D` and FLEX reconstruct state maps through the
   same service, and so do the refinement's startup and final reconstructions, which run
   `reconstruct3D`. The refinement's per-iteration partial reconstructions call the gridding and
   PCG reconstruction routines directly, with the same weight handling.
3. **Fractional reconstruction.** Both backends (gridding and PCG), in shared memory and
   distributed, with trailing reconstruction and the final native reconstruction, reconstruct a
   state from weighted particles.
4. **Per-state M-estimation in the refinement.** `refine3D_auto state=X m_estimator=flex` refines
   one FLEX state with every particle weighted by its frozen FLEX weight for that state.

FLEX's non-test source went from 21,051 to 16,532 lines (44 to 37 files). Its six state
reconstruction modules and its weight store with the store's file codec are gone; SIMPLE does
their work. Five general clustering modules (each with a tester), a shared PCG lattice layer, a
sample-set gather and a multi-volume insertion entered the SIMPLE library.

Status. Everything is written and was read by three independent reviews. On the maintainer's
first builds of 2026-10-08 the compiler found three name clashes, the fast test gate failed one
check, and the first `flex_pca_blobs` run stopped on a per-iteration buffer that was allocated
twice; all are fixed. The library and high-level test tiers have not yet been run to the end on
the final code; section 11 lists what to test first.

## 2. Terms

- **Project file** (`*.simple`). SIMPLE's binary project: segments of records (micrographs,
  stacks, particles in 2D and 3D, class averages, and `out`). A particle record has fixed named
  slots, so per-particle vectors cannot live in it.
- **`out` segment.** The project's list of outputs: maps (`vol`), Fourier shell correlation
  files (`fsc`), sigma2 files and other artifacts, each an entry with an `imgkind`. Programs find
  their inputs here.
- **Particle layout digest** (`ptcl_layout_digest`). A 64-bit digest of the project's name and
  particle layout (stacks and particle order). A side file that stores one value per particle
  records it, so a file built for another project or layout is refused. Because the name is part of
  it, a renamed copy of a project cannot open its parent's side files.
- **State weight set.** Weights `w_is` in [0, 1] for particle `i` and state `s`, one file per state,
  with a manifest. Kind **partition**: every weighted row sums to one (mixture
  responsibilities). Kind **kernel**: rows need not sum to one (kernel weights on a path or an
  axis).
- **Applied mass, effective sample size, hard population.** Three different populations of a
  weighted state. Applied mass is the sum of the weights the backend actually uses; it conserves
  accumulators. The effective sample size (ESS) is `(sum w)^2 / sum w^2`; it says how many
  particles' worth of information a state carries and gates support quality. The hard population
  is the number of particles whose largest weight is the state; it is for reporting.
- **Symmetric positive-definite (SPD)** matrices and the **median absolute deviation (MAD)** appear
  in FLEX's linear algebra and robust statistics.
- **Gridding and PCG.** SIMPLE's two reconstruction backends. Gridding inserts every particle's
  Fourier plane into a 3D lattice with a Kaiser-Bessel (KB) window and divides by the sampling
  density. PCG (preconditioned conjugate gradients) solves the least-squares reconstruction
  problem iteratively.
- **Trailing reconstruction** (trailing chain). Under fractional updates only a sample of the
  particles is aligned per iteration. The trailing chain keeps the accumulated reconstruction of
  earlier iterations and blends it with the current one by the population rule: the weights that
  keep the blended map representing the state's current population.
- **M-estimator.** An estimator that weights each observation. Here `m_estimator=flex` means each
  particle enters a state's reconstruction with its FLEX weight for that state.
- **Work project.** A copy of a project holding one state as state 1, refined on its own, so the
  original project is never written.
- **Sigma2.** The per-particle noise power spectra SIMPLE's likelihood objective uses.
- **Plane cache.** FLEX's disk cache of prepared particle Fourier planes at the covariance box. Its
  contract record (the header that says what the cache holds) is written last; it is built under a
  name that carries the run token, a string unique to the run.
- **GMM AUTO.** FLEX's automatic placement: over-fit a tied-covariance Gaussian mixture, join
  components whose pairwise density is unimodal into macro-clusters, give every macro-cluster at
  least one state and share the rest of the states by mass.
- **Hybrid radius.** The shell below which FLEX's E-step uses exact Cartesian statistics; above it,
  the polar rings of a shared direction bank.
- **Recoded label accuracy, truth-matched map.** On a phantom with known conformers, a delivered
  state is matched to the truth conformer its map correlates with best; the accuracy is the fraction
  of particles whose delivered state is their truth's state.
- **Embedding artifact.** The file `flex_pca_embedding.bin`: the latent coordinates of every
  particle with their precisions, plus the calibrated noise scale and the deconvolved
  coordinates. A states-only rerun (`infile=`) starts from it.

## 3. Why

On master before this work, FLEX lived next to SIMPLE rather than on it. It had its own
gridding and PCG state backends, its own part files and delivery for state maps, a weight store
and file codec of 1,265 lines, copied from the sigma2 store, that no other code read, its own `vol_flex` volume
kind with special cases in project I/O and in the `refine3D_states` handoff, a nonuniform filter
behind a key no caller set, a plane cache and an embedding cache that a crash could leave half
written, and its own copies of a Gaussian mixture, k-means (twice), a diffusion k-center, extreme
deconvolution, complete linkage with a union-find, and three copies each of farthest-point
seeding and log-sum-exp.

Two scientific needs pointed the same way. FLEX's weights are fractional, and SIMPLE's
reconstruction was hard-state only. And the weights are most useful in the refinement, as a
different M-estimator, which only works if SIMPLE's own reconstruction reads them.

The rulings that shaped the work: FLEX reuses SIMPLE's machinery wherever it feasibly can; the
weights live in one file per state outside the particle records; FLEX computes the weights and
the refinement uses them frozen, never recomputing them; the refinement is a single-state
refinement of one state on a work project; and every clustering algorithm FLEX uses becomes a
general library tool.

## 4. What FLEX does now, end to end

The `flex_pca` application runs the same five steps as before. What changed is who does the
work inside them.

1. **Prepare.** The project gateway validates the input and the canonical sigma2 state. The
   shared-memory and distributed-master strategies prepare sigma2 before they generate the
   workers' job description, so the global-sigma fallback for an empty half reaches the workers.
   With `cache=yes` the plane cache is built in the particle cache's directory, under a name with
   the run token, within a quarter of the free disk space, with its contract record written last
   and an atomic rename into place.
2. **Embedding.** The paired probe fits, the merge of the two half fits and the embedding run as
   before. The mean and basis projections go through `simple_reconstructor` (a sample-set gather
   that reads any number of volumes with one KB window per sample). The M-step insertion is
   `simple_reconstructor`'s `insert_planes_multi`. The coupled PCG operator shares its lattice
   geometry and support with the production PCG solver. Particle preparation is SIMPLE's
   `prep_imgs4rec`. Once the master has the embedding, the resident planes are freed and the plane
   cache is deleted: nothing after the embedding reads planes.
3. **Infer states.** The latent deconvolution calibrates the noise from the even and odd half
   solutions and fits the mixture prior by extreme deconvolution (`simple_xd_gmm`), choosing the
   number of components by held-out likelihood (this latent mixture prior is not the probe fit's
   mixture prior, which is started by `mcfa_init`). Placement uses the diffusion embedding with
   farthest-point k-center (`simple_kcenter`), k-means (`simple_kmeans`) or the tied-covariance
   mixture (`simple_gmm`), as selected. States below the occupancy floor are pruned and the weights
   decide the labels (section 6.5). The bandwidth cross-validation scores trial maps from the
   reconstruction service.
4. **Reconstruct states.** With `preimage_auto=yes` the service reconstructs raw trial half maps,
   and the two-gate merge joins states whose maps agree within their own reproducibility, by
   complete linkage (`simple_hac`), and folds view-clustered states into their most similar map.
5. **Publish.** FLEX writes its tables, publishes the weights as the project's state weight set,
   registers the embedding artifact in the `out` segment (imgkind `flex_embedding`) and writes each
   particle's largest-weight state into the project (state 0 when it has no weight). With
   `rec_states=yes`, the default, it then runs `reconstruct3D m_estimator=flex` from the published
   set at the native box, and each state's map and FSC are registered as ordinary `vol` and `fsc`
   entries.

The distribution model is unchanged: one strategy per mode, a master that plans the rounds and
workers that run one part per stage. Workers now run under their own worker command of the
private executable, and the public `flex_pca` command refuses `part=`.

## 5. Where everything went

| What FLEX had | Where it is now | Notes |
|---|---|---|
| State backends (`simple_flex_pca_states_backend`, `_states_gridding`, `_states_pcg`), the state driver (`simple_flex_pca_rec3D`), state parts and delivery (`_state_parts`, `_state_delivery`) | The reconstruction service, `src/main/strategies/parallelization/simple_rec3D_service.f90`, and `reconstruct3D m_estimator=flex` | Deleted. Trial maps: the service, in the FLEX process, with the weights as a caller-owned table, at `box_crop`, raw. Delivered maps: `reconstruct3D`, at the native box. |
| Round-weights broadcast and the worker state stage | Gone | Workers no longer take part in the state stage. |
| Weight store `simple_flex_weights_state` and its codec `simple_flex_weights_file` (`flex_weights_state_NNN.bin`) | `src/main/project/simple_state_weight_set.f90` and `src/fileio/simple_state_weights_file.f90` | Deleted. One general weight set for any producer (section 6.1). |
| `vol_flex` volume kind, `flex_pca_state_NNN.mrc` names, the delivery Butterworth filter | Ordinary `vol` and `fsc` entries, `recvol_stateNN.mrc` | A `vol` lookup still falls back to `vol_flex` in projects of earlier releases, read-only. |
| `box_rec`, `smpd_rec`, `flex_rec_box`, `flex_rec_smpd` | The service's box rules | Trial maps at `box_crop`, delivered maps at the native box. |
| Consensus nonuniform filter (`apply_consensus_nu_filter`, `nufilt`) | Removed | FLEX does no nonuniform filtering. |
| Multi-target gridding insert and the coupled insertion | `insert_planes_multi` in `simple_reconstructor` | K target volumes, density shared, diagonal or as packed pairs. |
| The mean and basis projection, the banded mean projection and FLEX's own window gathers (latent operations, polar bank, hybrid E-step statistics) | `exp_samples` and `project_fplane` in `simple_reconstructor` | One KB window per sample, read from every volume. SIMPLE's own unused `project_polar`, `interp_cmat_exp`, `interp_rho_exp` and `insert_plane_oversamp_opt` were deleted with them. |
| The body of `prep_imgs4projected_model` and its masking option | `prep_imgs4rec` and `gen_rec_plane` in `simple_matcher_3Drec` | `prep_imgs4projected_model` remains as a thin FLEX adapter. Optional band, observation-model and whitening arguments. The unused masking path is deleted. |
| PCG operator geometry, wrap table, deposition envelope, support window | `pcg_lattice`, `src/main/volume/simple_pcg_lattice.f90` | `reconstructor_pcg` and `flex_pcg_t` both extend it. The coupled kernels and the per-voxel preconditioner stay in FLEX. |
| The PCG operator's self-test inside the production module | Submodule `src/main/flex/fit/simple_flex_pca_pcg_tester.f90` | Built only with tests. |
| Gaussian mixture in `simple_flex_pca_gmm` | `simple_gmm` | FLEX keeps `gmm_state_weights` (the standardised frame, floors and readouts) and the GMM AUTO placement. |
| Extreme deconvolution in `simple_flex_pca_deconv` | `simple_xd_gmm` | FLEX keeps the noise calibration and the projection and noise matrices. |
| k-means (two copies) | `simple_kmeans` | `kmeans_latent_targets` is a thin adapter over it; the probe fit's mixture-prior start (`mcfa_init`) calls `simple_kmeans` directly. |
| Diffusion k-center | `simple_kcenter` | FLEX keeps the diffusion embedding in front of it. |
| Complete linkage and union-find in the merge | `simple_hac` (formerly `simple_avg_linkage`) | FLEX keeps both gates. |
| Farthest-point seeding, equal-mass quantile start, log-sum-exp | `simple_clustering_utils`, `simple_stat` | One copy each. |
| SPD solves, inverse and log-determinant | `cholesky`, `chol_forward`, `chol_backward`, `spd_inverse`, `spd_logdet` in `simple_linalg` | FLEX keeps its rescale-and-ridge retry around them. |
| Median, MAD, k-th smallest, effective sample size, correlation | `median`, `mad`, `selec` (double precision), `kish_ess`, `pearsn_serial_8` in `simple_stat` and `simple_srch_sort_loc` | |
| Power iteration for the equal-occupancy axis | `jacobi` with a pinned sign | |
| Mean-scale broadcast | `arr2file` | |
| Cross-fit FSC record fields nobody read, the series mirror file | Removed | The record keeps what the ridge uses (version 2). |
| Embedding cache with the deconvolution appended in place, `flex_pca_noise_scale.txt` | One embedding artifact, written under a temporary name and atomically renamed | Registered in the `out` segment. |
| Sigma2 preflight and worker-key injection in the commander | `ensure_canonical_sigma_state` in the project gateway, a worker command of the private executable | |
| `simple_finch` | Deleted | It had no caller. |

## 6. What SIMPLE gained

### 6.1 The state weight set

`state_weight_set` (`src/main/project/simple_state_weight_set.f90`) is the only way weights reach
SIMPLE code.

- **Files.** One binary file per state, `state_weights_gGGGGGG_sNNN.bin` (six-digit generation,
  three-digit state): a magic, a header (version, particle count, state count, state, generation,
  layout digest, kind, size), the state's weight for every project row, and a hard-label flag per
  row. One text manifest, `state_weights.txt`, records the generation, kind, producer, state count,
  particle count and layout digest, and per state the file name, size, checksum, the sum of its raw
  weights, its ESS and its hard population. Consumers compute the applied mass after their own
  threshold. A work project's set also records its parent set's generation, digest and
  state.
- **Publication.** `publish` writes the generation's files, then the manifest under a temporary
  name with a flush and an atomic rename, then the project's `out` entry (imgkind
  `state_weights`), then deletes the superseded generation. A crash at any point leaves a
  complete set selected; a mixed set cannot be read.
- **Validation.** `new` validates the whole set once: manifest, file sizes and checksums, file
  headers, the layout digest and particle count against the current project, values in [0, 1],
  labels only on positive weights, the populations, and row sums for a partition set.
- **Use.** One state column is held at a time. Consumers read raw weights, applied weights (after
  their own threshold), the hard labels, the selection, and applied mass and ESS over a row subset.
  The set's identity is its generation plus its layout digest.
- **Derived and withdrawn sets.** `publish_work_state` writes one state of a parent set for a work
  project. `discard_state_weight_set` withdraws a set: the `out` entry first, then the files.

### 6.2 The reconstruction service

`simple_rec3D_service.f90` takes a request: the particle rows, the weight source (hard labels, the
project's weight set, or a caller-owned table), the backend, the output kind (raw half maps, or
delivered maps with assembly, FSC and post-processing), an optional file-name prefix for raw
output, an explicit publication policy and the dispatch (in-process, or the queue rounds). The
`reconstruct3D` strategies now only parse, partition, prepare sigma2, sample and send the request.
With hard labels the outputs of all four paths (gridding and PCG, each in shared memory and with
two parts) were byte-identical to the references before the extraction.

### 6.3 Fractional reconstruction

`reconstruct3D m_estimator=flex` reconstructs every state of the project's weight set;
`state=X` restricts it to one state. Without `m_estimator=flex` a project that carries a weight
set reconstructs exactly as before.

- **Gridding.** `insert_plane_oversamp` takes an optional weight that scales the data term and
  the density term alike. A particle enters a state when its weight is above zero, and is read
  once per state it weighs into. Because fractional weights make the sampling density sparse and
  irregular, the assembly applies a shell-relative density floor before the division, but only
  when some weight of the state lies strictly between 0 and 1. Hard and 0/1 runs are untouched.
- **PCG.** A particle enters a state when its weight is above 0.01 (`PCG_WEIGHT_THRESHOLD`). Its
  noise spectrum is divided by its weight, which weights both the right-hand side and the
  density. The raw accumulator files (version 2) carry the applied mass and the weight-set
  identity, and the master refuses a part built under another set.
- **Trailing and update fractions.** Under a weight set the update counts, the realized
  fractions, the seed population and a chain's represented population are applied masses. The
  population rule (`population_blend_weights`) is generic over counts and masses with identical
  arithmetic for counts. A chain built under one weight-set identity and read under another is
  re-seeded, never blended. The PCG chain now has the same manifest contract as the gridding chain.
- **Population gates** that skip empty states use the ESS.
- **0/1 weights equal to the hard labels reproduce the hard maps and half maps byte for byte**,
  for both backends.

### 6.4 Per-state M-estimation in the refinement

`refine3D_auto state=X m_estimator=flex`:

1. Copies the project into the run's own directory as `refine3D_auto_stateXX.simple` and selects
   into it every particle with a positive weight for state X. State X's map and FSC become state
   1's. The copy gets its own single-state weight set, derived from the parent's with the parent's
   identity recorded. The parent project is never written. `state=` needs `mkdir=yes`.
2. Runs the normal single-state `refine3D_auto` on the copy. Alignment and sampling are unchanged
   and never read the weights: every selected particle is aligned against the state's reference.
   Every reconstruction (the startup reconstruction, the iterations' partial reconstructions in
   shared memory and distributed, trailing, and the final native reconstruction) is weighted.
3. Registers the final map with its hard population, applied mass and ESS, and checks that the
   weight set did not change during the run.

`refine3D_auto state=X` without `m_estimator=flex` is the hard comparison on the particles
labelled X. `m_estimator=flex` without `state=` is refused, and so is `m_estimator=flex` on
`refine3D_states` or on a multi-state `refine3D`: a search over states moves the labels the
frozen weights were computed for. `refine3D_states` turns the FLEX solution into a hard
assignment for its initialization and withdraws the weight set after reading the labels and
maps.

### 6.5 Pruning

When FLEX's occupancy floor drops a state, the weights decide: a row with weight in a kept state
is labelled by its largest kept weight, partition rows are renormalized over the kept states,
and a row without kept weight is unassigned. Before, the particles of a dropped state got label 0
but kept their weights in the surviving columns.

### 6.6 Library modules

All in `src/utils/clustering` unless named otherwise. Each follows the library's standard: one
module per algorithm with a type of the same name, `new` with plain arrays and the algorithm's
parameters, one clustering call, getters for labels and representatives, `kill`, deterministic
seeding and ties, labels by decreasing population, and a tester. Extreme deconvolution takes its
data per call instead of in `new`, because its projection and noise matrices are d x d per
particle.

- `simple_hac`: agglomerative clustering of a distance matrix, average or complete linkage,
  stopped at a cluster count or a distance threshold, with a mask of entries that stay singletons
  and the merge history.
- `simple_kmeans`: seeded at the point nearest the mean, then farthest-point; Lloyd iterations to
  convergence; an empty cluster is reseeded at the worst-fitted point; optional per-dimension
  weights.
- `simple_kcenter`: greedy farthest-point k-center from the point farthest from the mean.
- `simple_gmm`: a tied-covariance mixture by expectation maximization, with a mixing-proportion
  floor, a responsibility floor, the respawn of redundant components, the Bayesian and integrated
  completed likelihoods (BIC, ICL) and the pairwise Mahalanobis separation.
- `simple_xd_gmm`: extreme deconvolution with per-point projections and noise covariances, the
  posterior per point, and `xd_select_k` for the number of components by held-out likelihood on two
  strided halves.
- In `src/main/volume`: `pcg_lattice`, `exp_samples`, `insert_planes_multi`.
- Numerics: the Cholesky family in `simple_linalg`; `logsumexp`, `kish_ess`, double-precision
  `median` and `mad` in `simple_stat`; double-precision `selec` in `simple_srch_sort_loc`.

## 7. What FLEX keeps

The method is FLEX's and stays in `src/main/flex`:

- the probe fit (`simple_flex_probe_fit` and its four submodules: E-step, update, engine,
  cross-fit), the posterior, the coupled M-step and its PCG operator's coupled kernels and per-voxel
  preconditioner, the polar E-step bank, the cross-fit FSC ridge, the paired merge of the half fits
  and the embedding;
- the latent deconvolution's noise calibration and its projection and noise matrices;
- the state placement decisions (diffusion embedding, paths, equal occupancy, kernel bandwidths,
  the population floor, GMM AUTO's over-fit, unimodality merge and seat apportioning), the pruning,
  the bandwidth cross-validation and the merge gates;
- UMAP, the readouts and figures;
- the probe and embedding part files of the distributed fits, which carry sufficient statistics no
  SIMPLE file holds.

The folders are `run/` (the application, its contracts and the services it composes), `fit/`
(the probe fit and its kernels), `states/` (`weights`, `gmm`, `targets`, `deconv`, `merge`) and the
folder root (helpers, figures with UMAP, the class-expansion program and the testers).
`scripts/check_flex_dag.py` layers the modules themselves, not the folders: contracts first, then
leaves and codecs, then the probe-fit services, then state inference, then the application
services, the application and the strategy. It runs before the fast gate.
`src/main/flex/README.md` describes every module.

## 8. What a user sees differently

- **Keys removed:** `box_rec` (the service's box rules replace it; its derived `smpd_rec` went
  with it) and `nufilt`.
- **Keys added:** `m_estimator=no|flex` and `state=` on `reconstruct3D` and `refine3D_auto`.
- **Outputs.** State maps are `recvol_stateNN.mrc` with their FSC files, registered as `vol` and
  `fsc`; there is no `vol_flex` and no `flex_pca_state_NNN.mrc`. The weights are
  `state_weights_g*_s*.bin` with `state_weights.txt`. The embedding artifact is registered in the
  `out` segment. `flex_pca_noise_scale.txt` and the cross-fit series file are gone.
- **The plane cache** no longer outlives the run.
- **Logs.** The mixture and extreme-deconvolution libraries write a few lines under their own
  prefixes (`>>> GMM components ... respawning`, `>>> XD K=...`); FLEX's own k-means line is
  unchanged.
- **refine3D_states** refuses `m_estimator=flex` and leaves no weight set in its output project.

## 9. Behaviour changes to expect

Most of the reuse is exact. These parts are not, by design:

- **State numbering.** Clusters and mixture components are now labelled by decreasing
  population (extreme-deconvolution components by mixing proportion). FLEX's states can come out
  in a different order than before. The class expansion program also uses the k-means adapter,
  so the order of its two children can change.
- **Merge ties.** The merge's ties follow the library's order. This differs only on exact ties.
- **Probe-fit mixture prior start.** `mcfa_init` now seeds k-means at the point nearest the mean and runs to
  convergence, instead of starting from the point farthest from the origin and running twelve
  iterations; an empty cluster is reseeded instead of collapsing to zero.
- **Rounding-level changes.** The coupled insertion computes sample locations directly instead of
  incrementally; the mean projection uses the separable sample-set gather; the Cholesky routines
  subtract term by term; log-sum-exp normalizes by `exp(x - lse)`; the split-half reliability is
  summed in double precision but returned in single precision by `pearsn_serial_8`.
- **Inference maps.** The bandwidth cross-validation and the merge now read raw half maps from the
  service instead of maps that went through FLEX's delivery filter. This changed inference when
  it was introduced and was measured against truth (section 10).
- **Cross-fit log.** The reported resolution is the last shell above the threshold instead of the
  first below it.

## 10. What was measured

The foundations, the fractional reconstruction, FLEX on the service and the per-state
refinement were implemented and measured by an autonomous run on the group's Linux test
workstation on 2026-10-07, starting from master `bd3997758`. The reuse work (the sample-set gather
and multi-volume insertion, the PCG lattice layer, the clustering modules, the numerical helpers,
the plane cache and the embedding artifact) and the documents were written by hand afterwards
and have not been measured yet.

**Hard reconstruction unchanged.** After the service extraction, all four `reconstruct3D` paths
(gridding and PCG, shared memory and two parts) produced byte-identical maps to the references
recorded before any change.

**Fractional phantom** (2,000 particles, box 64, two conformers A and B differing by one lobe
moved 8 pixels; state 1 weighted 0.75 on A and 0.30 on B with a small modulation, state 2 the
rest). Correlation of each weighted map with its expected mixture, against the correlation of the
hard map with its truth:

| Backend | State | Mass from A / B | Weighted map vs mixture | Hard map vs truth |
|---|---|---|---:|---:|
| gridding | 1 | 749.9 / 299.8 | 0.9856 | 0.9857 |
| gridding | 2 | 250.1 / 700.2 | 0.9853 | 0.9852 |
| PCG | 1 | 749.9 / 299.8 | 0.9981 | 0.9977 |
| PCG | 2 | 250.1 / 700.2 | 0.9981 | 0.9978 |

0/1 weights equal to the labels reproduced the hard maps byte for byte, and both trailing chains
represented the state's applied mass (expected 1049.693; gridding and PCG 1049.694).

**FLEX on the service** (application phantom, four arms: gridding or PCG for the basis and for
the state maps). Every arm delivers 3 states with recoded label accuracy 0.9995, as before. The
truth-matched maps correlate with their truth at least as well as before (before in brackets):

| Arm | Matched A | Matched B |
|---|---:|---:|
| grid_basis_grid | 0.9857 (0.9851) | 0.9852 (0.9846) |
| pcg_basis_grid | 0.9857 (0.9851) | 0.9852 (0.9846) |
| grid_basis_pcg | 0.9985 (0.9976) | 0.9987 (0.9979) |
| pcg_basis_pcg | 0.9985 (0.9976) | 0.9987 (0.9979) |

The workflow gate (`flex_pca_blobs`) gave accuracy 0.99950 in shared memory and 0.99750
distributed, as at the start; the two modes' labels agree on 99.80 % of the particles. (On master
the distributed mode failed because the phantom project lacked the computing-environment segment;
the start values came from a build with that test-only fix, which is part of this work.) A fifth
arm with the bandwidth cross-validation and a covariance box of 48 of 64 reached accuracy 0.8655
on the Release build and 0.8775 on Debug; a Debug run without the cross-validation gave 0.9165, so
most of the drop comes from the smaller box. It has no earlier baseline.

**Per-state refinement phantom** (inside `flex_pca_blobs`): each truth-matched FLEX state refined
by `refine3D_auto state=X`, once with `m_estimator=flex` and once on the hard labels, five
iterations at most. Correlation of the refined map with its truth:

| State (truth) | Work project with `m_estimator=flex` | Weighted | Hard |
|---|---|---:|---:|
| 3 (A) | 1,332 particles, applied mass 1016.7, ESS 1050.5, hard population 996 | 0.98979 | 0.98723 |
| 1 (B) | 1,332 particles, applied mass 1029.9, ESS 1071.0, hard population 1001 | 0.99082 | 0.98622 |

The weighted refinement includes particles of the other conformer with small weights and
correlates slightly better with the truth. The parent project was byte-identical after all four
refinements.

**Real data** (beta-galactosidase): the weighted and the hard refinement of the most populated
state reached 3.98 Å and 3.89 Å (FSC 0.143). FLEX published 0/1 weights on that data set, so this
checks the weighted machinery, not FLEX's weights. Whether FLEX's weights help on real data is to
be tested once the implementation is validated on simulated data.

## 11. What to test first

In this order; each step guards the next. Commands run from the build directory with tests
(`BUILD_TESTS=ON`); `simple_test_exec` runs one entry, `ctest -R <name>` runs it as registered.

1. **The fast gate**, on a Debug build (bounds checking): `scripts/run_fast_gate.sh build` from
   the repository root, or `ctest -L fast` in the build directory. It covers the clustering testers (`unit_numerics`), the state weight set
   with its withdrawal and the segment clearing (`unit_project`), the weighted insertion
   (`unit_reconstruction`), and FLEX's testers with the PCG operator self-test
   (`unit_heterogeneity`).
2. **The hard paths byte for byte.** The reuse work touched two production paths that must not
   change: particle preparation for every gridding reconstruction (`prep_imgs4rec`) and the PCG
   solver's geometry (`pcg_lattice`).
   - `simple_test_exec test=pcg_recon`.
   - The four reference reconstructions on the test workstation, from a fresh copy of the
     beta-galactosidase result project (`keep/B/1_refine3D_states/bgal.simple` under
     `~/agent_runs/flex_on_simple/scratch`), compared with the kept references in `keep/refs/`:
     `simple_exec prg=reconstruct3D projfile=bgal.simple mskdiam=180 pgrp=d2 nstates=2` with
     `rec_backend=gridding nthr=24`, `rec_backend=gridding nparts=2 nthr=12`,
     `rec_backend=pcg nthr=24` and `rec_backend=pcg nparts=2 nthr=12`. Every `.mrc` file should be
     byte-identical (24, 32, 16 and 16 files).
3. **The reconstruction library tier:** `simple_test_exec test=lib_reconstruction`, in particular
   the `fractional_reconstruction` sub-suite (the table in section 10, and 0/1 weights
   byte-identical to hard).
4. **The FLEX application phantom:** `simple_test_exec test=lib_heterogeneity`. Every arm against
   the numbers above: accuracy at least the earlier value minus 0.02 and at least 0.95, matched
   maps within 0.02 of the earlier correlation, each matched map preferring its own truth. State
   numbers may permute (section 9).
5. **The FLEX workflow gate:** `simple_test_exec test=flex_pca_blobs nthr=4` (as CTest runs it; the
   per-state refinement at its end uses twice the threads). The distributed mode exercises the
   worker command, the plane cache adopted by the workers and deleted by the master, and the
   embedding artifact.
6. **The refinement entries** that passed before, for regressions from the reconstruction
   plumbing: `simple_test_exec test=simulated_workflow suite=1jxy`, the same with `suite=6vxx`,
   `simple_test_exec test=solve3D_addon nthr=8` and `simple_test_exec test=cont_refine3D_1jxy
   nthr=8`.

## 12. Open points

- **Polar bank and symmetry.** The E-step's polar bank snaps each particle to the bank direction
  nearest its raw plane normal without mapping it into the asymmetric unit first. Particles that
  SIMPLE refined already lie in the bank's coverage (on a d2 project all 5,513 particles snapped
  within 3.05°), but imported or symmetry-expanded orientations would not: over random d2
  orientations the mean snap error was 31.7° without the mapping and 1.25° with it. Only the polar
  rings above the hybrid radius use the snapped direction. The fix (apply the symmetry mapping
  first) is still to be made.
- **UMAP default.** UMAP was ruled off by default; the code still defaults to `umap=yes`.
- **The registered embedding entry** has no reader yet; a states-only rerun still takes
  `infile=`.
- **Rows with weights between 0 and 0.01** are aligned in a PCG per-state refinement but not
  reconstructed.
- **Library dependencies.** `simple_kmeans`, `simple_kcenter` and `simple_xd_gmm` import
  `simple_clustering_utils` for the seeding, which brings the distance-matrix algorithms into their
  compile chain.
- **Emptying a project segment.** Writing a segment inside a project file does nothing for an
  empty table, because an unread segment looks the same in memory. Emptying a segment on disk is
  now an explicit call, `sp_project%clear_segment_inside`; the weight-set withdrawal uses it. Other
  code that removes the last entry of a segment must call it too; existing callers were not
  audited.

## 13. Working on FLEX from here

- In FLEX, reconstruct through the service, read weights through `state_weight_set`, and publish
  maps as `vol` and `fsc`. Do not add a state backend, a weight file or a volume kind.
- A clustering or numerical routine that is not specific to FLEX goes into the library with a
  tester, not into `src/main/flex`.
- Keep the module layering that `scripts/check_flex_dag.py` checks. One module or submodule per
  file, `use ..., only:` imports, no run state at module scope, and
  self-tests in `*_tester.f90` files.
- Modern Fortran only: no `goto`, no numbered labels for control flow. Fortran names are not case
  sensitive, so a local `r` clashes with an argument `R`, and a function `mean_scale_fname` with a
  constant `MEAN_SCALE_FNAME`.
- Fortran does not short-circuit `.and.`: never test an element of an array that may be
  unallocated on the same line as the flag that guards it.
- Measure changes against truth on the phantoms (`lib_heterogeneity`, `flex_pca_blobs`), not by
  FSC alone, and record the numbers.

## 14. How the work was done

The plan was written and ruled by the maintainer on 2026-10-07. An autonomous run carried out
the foundations, the fractional reconstruction, FLEX on the service and the per-state refinement
on the test workstation the same day, one step at a time with builds, the test tiers and the
measurements above; the maintainer stopped it after the per-state refinement. Its changes were
applied to master `92d684c6f`. The reuse of SIMPLE machinery, the clustering modules and the
documents were then written by hand on 2026-10-07 and 08 as one reviewable change per item, read
by three independent reviews (which found two build errors and one failing test, fixed in their
items), and followed by the maintainer's rulings of 2026-10-08: the PCG membership threshold of
0.01, the plane cache deleted after the embedding, no weighted multi-state refinement, and no
weights kept by `refine3D_states`.
