# FLEX on SIMPLE machinery: fractional state weights and reuse

Date: 2026-10-07. Revision 4: the maintainer's review of revision 1 added the foundations of
section 4 and the contracts of sections 5 and 6; revision 4 makes the refinement a per-state
M-estimation in `refine3D_auto` with frozen FLEX weights (rulings 7 to 9) and adds the clustering
extraction (section 7.1). Every decision of section 3 is ruled. Status: planned. Owner: the
maintainer.

This note is the maintainer's target. Implementing agents record progress in section 12, and add
a row to the file table of section 11 (with a one-line reason) before editing a file it does not
list; nothing else. The goal, rulings, contracts, tolerances and exits change only by a maintainer
ruling. The plan is carried out by an autonomous run, phases 0 to 6.

## 1. Goal

1. **Fractional state weights become a SIMPLE capability.** FLEX computes, per particle and state, a
   weight `w_is` in [0,1] and publishes them as a transactional set of per-state files registered in
   the project. The standard reconstruction (gridding and PCG, shared memory and distributed,
   trailing, final native reconstruction) reads them.
2. **FLEX weights drive the refinement as an M-estimator.** `refine3D_auto state=X m_estimator=flex`
   refines state `X` as a single-state refinement in which every particle contributes with its frozen
   FLEX weight for `X`. The weights are never recomputed during the refinement.
3. **One reconstruction service.** `reconstruct3D`, the refinement and FLEX reconstruct state maps
   through the same service, extracted from today's `reconstruct3D` strategy.
4. **FLEX reuses SIMPLE machinery wherever it feasibly can.** What stays in `src/main/flex` is the
   method itself: the probe fit and its E-step, the coupled M-step, the latent deconvolution's noise
   model, state placement decisions and the merge gates. FLEX stops reconstructing maps, storing
   weights, writing part formats, publishing special volume kinds and carrying clustering code of its
   own: its clustering is extracted into `src/utils/clustering` as general tools.

Not a goal: new scientific features, and recomputing state weights from the refinement's
likelihoods.

## 2. Where we are

On master: the `prob_state` retirement (`fbffebf73`) and the first FLEX refactoring pass
(`878adf2d4`).

**History.** Soft multi-state reconstruction existed in 2020-2021: each particle carried several
weighted orientations with their own states (`insert_planes_2(state=)`, `grid_ptcl_2`,
`reconstruct3D rec_soft=yes`). The search never produced mixed-state weights, because
`states_reweight` kept only the state with the largest summed weight. The code went in
`e8779a629`, `3245e505d`, `4887ae235` (2021-07) and `7da020aac` (2022-02). A per-particle insertion
weight (`pwght`, keys `ptclw`, `cavgw`) multiplied both the backprojection and the density until
`83571b305` (2026-05-24); its project slot is now `I_RETIRED_W`. Soft state weights were planned
for the PCG operator ("weight both B and D", `doc/implementation_notes/completed/pcg_priors_history.md`)
but never built.

**Standard reconstruction** is hard-state only. `insert_plane_oversamp` takes no weight
(`simple_reconstructor.f90`); `update_state_half_rec` stops on a mixed-state batch
(`simple_matcher_3Drec.f90`); the PCG worker selects by hard state, and its raw accumulators carry
an integer particle count that the PCG trailing chain uses as the represented population
(`write_raw_accum` in `simple_reconstructor_pcg.f90`, the chain weights in
`simple_rec3D_pcg_strategy.f90`). The gridding trailing chain has a manifest with a real
represented population (`simple_trail_chain_manifest.f90`). The `reconstruct3D` strategy owns
particle sampling, sigma2 setup, assembly, canonical file names, project registration and
post-processing in one lifecycle (`simple_rec3D_strategy.f90`); registration uses `mkdir=yes` as
its publication switch, which a nested internal run cannot use. `calc_final_rec` builds its child
command line from scratch and derives `pop` from hard state counts (`simple_final_rec.f90`).

**Single-state refinement of one state** exists today only as `solve3D`'s state continuation
(`prepare_state_continue_project` in `simple_commanders_solve3D.f90`): it copies the project,
selects the rows labelled with the state into the copy and refines that copy as a single-state
project. `refine3D_auto` has no `state` input.

**FLEX** reconstructs its weighted state maps itself: a multi-target gridding insert used only by
its gridding backend, a PCG backend that re-reads particles with `sig2/w` above a 1e-3 floor, its
own part files, a shell density floor and its own delivery filter. Its bandwidth cross-validation
runs two reconstructions per bin (even, then odd) and scores maps that went through that delivery.
It publishes `vol_flex`, special-cased in project I/O and in a parallel `refine3D_states` handoff.
It can apply one nonuniform filter from the consensus half maps to every state
(`apply_consensus_nu_filter`, behind `nufilt=yes`); the key defaults to `no`, is not in the UI and
no caller sets it. It stores the weights in `flex_weights_state_NNN.bin`, a 1,265-line copy of the
sigma2 state store; no production code reads them. When its occupancy floor drops a state, the
particles labelled with it get label 0 but keep their weights in the surviving columns
(`prune_underpopulated_states`). Its plane cache is written in place, and its embedding cache
appends the deconvolved block in place, so a crash can leave a torn resume artifact. It carries its
own clustering: a tied-covariance Gaussian mixture with a hierarchical choice of the number of
components (`simple_flex_pca_gmm.f90`), k-means and a diffusion k-center (`simple_flex_pca_targets.f90`,
with a second k-means copy in the mixture prior of `simple_flex_pca_posterior.f90`), extreme
deconvolution (`simple_flex_pca_deconv.f90`), complete-linkage agglomeration with its own union-find
(`simple_flex_pca_merge.f90`), and three copies each of farthest-point seeding and log-sum-exp.

**The project file cannot hold per-particle vectors.** Particle records are fixed 64-slot real
records; keys outside the named slots are dropped on write (`simple_binoris_tester.f90` asserts
it). Per-particle vectors live in side files registered in the `out` segment and tied to the
particle layout by a digest; that digest is `sigma2_state_project_layout_digest`, already used by
twelve files outside the sigma2 store.

**The library's clustering standard** is set by `simple_avg_linkage` (hierarchical agglomerative
clustering, average linkage), `simple_kmedoids` and `simple_aff_prop` in `src/utils/clustering`: one
module per algorithm, `simple_<algorithm>` holding the type `<algorithm>`, with `new` taking plain
arrays and the algorithm's parameters, one clustering call returning labels and representatives,
`kill`, deterministic tie-breaking, labels ordered by decreasing population, a tester per module, and
`simple_clustering_utils` dispatching the distance-matrix algorithms (`cluster_dmat`).

## 3. Rulings

Rulings (2026-10-07):

1. FLEX reuses existing SIMPLE machinery wherever it feasibly can.
2. Fractional state weights are needed, for a different M-estimator and for using FLEX state
   weights in the refinement.
3. Weights are written as one weight distribution per state, outside the particle records.
4. The foundations of section 4 (particle-layout identity, the transactional `state_weight_set`,
   the reconstruction service) come before any weighted reconstruction.
5. `vol_flex` is retired: FLEX publishes ordinary `vol` and FSC entries. Old projects keep a
   read-only `vol_flex` fallback.
6. Final native reconstruction carries the weights.
7. FLEX computes the state weights; the refinement uses them frozen. Recomputing them during the
   refinement is not part of this plan.
8. The refinement is `refine3D_auto state=X m_estimator=flex`: a single-state refinement of state
   `X` on a work project, as `solve3D`'s state continuation does (section 6).
9. Whatever clustering FLEX uses is extracted, generalized and made to follow the library's
   clustering standard (section 7.1).

Decisions, all ruled 2026-10-07:

| # | Decision | Ruling |
|---|---|---|
| D1 | Poses under fractional weights | Within one run every state shares the particle's pose. Each per-state refinement (ruling 8) refines its own work project, so its poses belong to that run; the parent project's poses are untouched. |
| D2 | Where the weights come from during refinement | From FLEX, frozen (ruling 7). |
| D3 | Weight kinds | PARTITION (rows sum to 1, mixture responsibilities) and KERNEL (FLEX's equal-occupancy and single-axis bumps). A per-state M-estimation uses one column, so `refine3D_auto` and `reconstruct3D` accept either kind; the kind is recorded with every result. |
| D4 | Gridding density floor for fractional weights | `floor_rho_shellwise` whenever the weights are fractional; never on hard runs. PCG keeps its shell-floored preconditioner. |
| D5 | Activation | `m_estimator=no|flex {no}` with `state=X` on `refine3D_auto` and `reconstruct3D`. With `no`, a project that carries a weight set reconstructs and refines exactly as today; `state=X` alone selects the rows labelled `X` (the hard comparison). `refine3D_states` gets no weighted mode. |
| D6 | FLEX probe and embedding part codecs | Keep them: they carry FLEX sufficient statistics that no SIMPLE file holds, and SIMPLE has no generic summed-array part file. |
| D7 | FLEX plane cache | Keep its payload (prepped Fourier planes, persistent across runs for `infile=` resumes). Adopt the standard cache machinery of `simple_ptcl_cache`: cache directory, run token, space budget and key-last publication. |
| D8 | UMAP and its figures | Keep in FLEX (no SIMPLE equivalent). Default `umap=no`. |
| D9 | FLEX state-map box | The service's box rules: `box_crop` for trial maps, native for delivered maps, as `reconstruct3D` does. `flex_rec_box`, `flex_rec_smpd` and the temporary `smpd_crop` override go. |
| D10 | Rows whose argmax state is pruned | The weight set's selection is authoritative: a row is active when its applied weights sum above zero; its label is the argmax over the kept states; PARTITION rows are renormalized over the kept states; a row with no kept weight is deselected. `state > 0` then agrees with the set. |
| D11 | Maps FLEX infers from (bandwidth CV, merge) | Raw half maps from the service, no delivery filter, FLEX's own masks for its statistics. This changes inference relative to today (CV now scores delivered maps); Phase 3 measures it against truth. |
| D12 | Nonuniform filtering of FLEX state maps | None. FLEX does no nonuniform filtering: the opt-in consensus filter, the `nufilt` parameter and its parser entry go. |

Ruling of 2026-10-07 during the run, the base of the run: the run starts from master at `bd3997758`
(the commit after `878adf2d4` that carries this plan and one stream change), not from `878adf2d4`.
The base is the commit the run's working copy started from; every Phase 0 number is measured on a
build of it.

Rulings of 2026-10-07 after Phase 4. The autonomous run stopped after Phase 4 by the maintainer's
decision; its diff was applied to master at `92d684c6f`, and Phases 5 and 6 are done by hand.
- *Cross-mode tolerance (Phase 2).* Accepted as judged: the tolerance of 8.2 applies to the
  reconstruction's half maps before the ML prior. The prior is computed from each run's own FSC and
  amplifies summation-order differences, with hard labels as well as with weights.
- *Long tests.* A test of the per-state refinement's length (463 s) is never part of the build; it
  belongs to the nightly run. `flex_pca_blobs` is a `highlevel` entry: the build runs only the `fast`
  tier, and the nightly test workflow runs `highlevel`. The per-state refinement stays inside
  `flex_pca_blobs`; the CTest budget is raised if a later test needs its own entry.
- *Real data.* Whether FLEX's weights help on real data is tested after this implementation is
  finished and validated on simulated data. Until then the real-data exit of Phase 4, which ran with
  the 0/1 weights FLEX publishes on beta-galactosidase, counts as a check of the weighted machinery,
  not of FLEX's weights.

Rulings of 2026-10-07 on the Phase 5 interfaces:
- *Dead code.* Code found dead is deleted: the multi-target gridding insert, `project_polar`,
  `insert_plane_oversamp_opt`, the image-masking option and whatever else turns up.
- *The embedding artifact* is registered in the project's `out` segment through the `oris`
  getters and setters of `os_out`; the container is the interface consumers read.
- *Gather and insertion.* One sample-set gather in `simple_reconstructor` (window geometry once per
  point, read from any volume on the lattice); `project_fplane` and FLEX's gathers use it. FLEX's
  coupled insertion moves to `simple_reconstructor` with shared, diagonal and pair densities;
  `insert_plane_oversamp` stays separate.
- *PCG lattice layer.* Only the code identical today is shared (lattice geometry, support window,
  deposition envelope, kernel-finalisation core). The FLEX PCG self-test moves into a submodule
  file and the test environment policy is updated.
- *Clustering.* The five modules of 7.1 with tied-covariance mixtures only; the number of
  components is chosen by held-out likelihood in the extreme deconvolution; FLEX's over-fit and
  merge stays FLEX. The mixture-prior k-means of the probe fit switches to `simple_kmeans`
  (seeding at the mean, convergence, empty-cluster recovery).
- *Numerical helpers.* A double-precision SPD Cholesky family in `simple_linalg` (FLEX keeps its
  rescale-and-ridge retry); double-precision `median` and `mad` in `simple_stat`; `jacobi` with a
  pinned sign for the equal-occupancy power iteration; the cross-fit FSC artifact loses its
  write-only fields; the basis Wiener filter keeps its 0.999 cap; the mean-scale broadcast uses
  `arr2file`.

Rulings of 2026-10-08 after Phase 5:
- *PCG membership threshold.* `PCG_WEIGHT_THRESHOLD` is 0.01 (the code's value), not the 1e-3 of
  section 5 and the Phase 2 record.
- *Plane cache lifetime.* The run that built the plane cache deletes it once the master has its
  embedding.
- *No weighted multi-state refinement.* `refine3D_states` has no weighted mode (revision 4); the code
  now refuses `m_estimator=flex` there and on a multi-state `refine3D`.
- *No persisted weights in `refine3D_states`.* It turns the flex solution into a hard assignment for
  its initialization, so it withdraws the weight set after reading the labels and state maps.

## 4. Foundations

### 4.1 Particle-layout identity

`sigma2_state_project_layout_digest` becomes `ptcl_layout_digest` in a small module of its own
(`src/main/project/simple_ptcl_layout.f90`). The sigma2 store, the weight set, the trailing
manifests and the solve3D manifest call it. Same digest, same values: a rename and a move.

### 4.2 `state_weight_set`

`src/main/project/simple_state_weight_set.f90` holds the type `state_weight_set`;
`src/fileio/simple_state_weights_file.f90` holds the file codec. They replace
`src/main/flex/states/simple_flex_weights_state.f90` and `src/fileio/simple_flex_weights_file.f90`.

**Storage.** One file per state, `state_weights_gGGGGG_sNNN.bin` (generation `G`, state `N`), and
one manifest `state_weights.txt` published last. The manifest records the generation, the kind,
the producer, the state count, the particle count, the layout digest, each file's byte size and
checksum, and per state the applied mass, the effective sample size and the hard population. The
project's `out` segment points at the manifest (`imgkind='state_weights'`). Publication writes the
generation-stamped files, then the manifest by temporary name and rename, then the project
pointer, then deletes older generations. A crash between steps leaves the previous generation
intact and selected; a mixed delivery cannot be read. This is the pattern of
`simple_trail_chain_manifest.f90`, generalized.

**Object.**
- `new` validates the complete set once (manifest, sizes, checksums, digest against the current
  project) and fails loudly.
- One state column is loaded at a time, on demand; memory does not grow with `nstates`.
- It exposes the selected rows, the raw and the applied weights (after the consumer's threshold),
  and per state the applied mass, the effective sample size and the hard population, kept apart.
- It exposes the kind and the set identity (generation plus digest).
- It provides the reducers the reconstruction needs: weighted population and update counts per
  state for a given row subset.
- It writes a single-state set for a work project: one column, rows remapped to the work project's
  layout, the parent set's identity recorded as provenance (section 6).
- `oris` never opens weight files; consumers hold a `state_weight_set`.

**The three populations are not interchangeable.** Applied mass conserves accumulators (trailing
population rule, update fractions). Effective sample size gates support quality (population gates,
FLEX's occupancy floor). Hard population is label reporting. Each consumer names which one it uses.

`ptcl3D` `state` stays the selection flag and the maximum-weight label, consistent with the set
by D10. Old projects with `flex_weights` entries still load; the entries are ignored.

### 4.3 The reconstruction service

A lower-level request and service extracted from the `reconstruct3D` strategy
(`src/main/strategies/parallelization/simple_rec3D_service.f90`, name to settle at extraction).
The request carries:

- explicit particle indices;
- a weight source: hard labels, a registered `state_weight_set`, or a transient table owned by the
  caller (FLEX's trial weights);
- the backend (gridding or PCG);
- the output kind: raw half maps, or delivered maps (assembly, FSC, post-processing);
- a caller-owned file-name prefix;
- an explicit publication policy (register in the project or not), replacing the `mkdir=yes`
  proxy;
- the dispatch: in-process, or the standard queue rounds.

`reconstruct3D`'s strategy becomes sampling, sigma2 setup and canonical names around the service;
with hard labels its outputs are bitwise unchanged. FLEX calls the service for its trial and final
maps. Its bandwidth CV gets the even and odd maps of a bin from one call instead of two
reconstructions.

## 5. Fractional reconstruction

Through the service, for both backends, shared memory and distributed.

**Gridding.**
- An optional weight on `insert_plane_oversamp` scales both the data term and the density term,
  as before `83571b305`.
- Each `(state, half)` membership list holds the rows with applied weight above the gridding
  threshold, which is zero: every positive weight is inserted, as FLEX's gridding does today
  (FLEX's mixture model already floors responsibilities at 1e-3 and renormalizes).
- One half-map reconstructor at a time, as now: memory independent of `nstates`. A particle is
  read once per state it weighs into.
- Partial names and payloads unchanged; the master's summation unchanged.

**PCG.**
- Membership above the PCG threshold, 1e-3 (the `sig2/w` conditioning floor FLEX uses); the
  particle's noise spectrum divided by `w_is` before `prep_particles` weights both B and D.
- The raw accumulator header carries the real applied mass beside the integer contributor count,
  and the weight-set identity in its provenance. The master's summation is unchanged.

**Thresholds are per backend; every mass is computed from the weights as applied, after the
threshold.**

**Trailing reconstruction and update fractions.** The represented population of a chain becomes
the applied mass: the gridding manifest already holds a real population; the PCG chain moves to the
same manifest contract (real mass, contributor count, weight-set identity) instead of reading
integer header counts; the frozen and add-on accumulators record the same provenance.
`population_blend_weights`, the state update counts and fractions, the PCG master's representation
counts and the seed population count applied mass. The weights are frozen for a run, so a chain
built under one weight-set identity and read under another (FLEX re-run in between) is re-seeded,
never blended. The sampler never reads weights (sampling is driven by class averages only).

**Population gates** that skip empty states use the effective sample size.

**`reconstruct3D state=X m_estimator=flex`** reconstructs state `X` from every row with a positive
applied weight for `X`; without `state`, every state of the set.

**Hard runs are bitwise unchanged** with `m_estimator=no`.

## 6. Refinement: per-state M-estimation in `refine3D_auto`

`refine3D_auto state=X m_estimator=flex`:

1. **Inputs.** The project carries a published weight set (FLEX's) with state `X`. The starting map
   is state `X`'s `vol` entry, which FLEX publishes (ruling 5), unless `vol1` is given.
2. **Work project.** As `solve3D`'s state continuation does: copy the project; select the rows with
   a positive applied weight for `X` (with `m_estimator=no`, the rows labelled `X`) into the copy as
   a single-state project; write the copy's single-state weight set (column `X`, rows remapped, the
   parent set's identity as provenance). The copy is refined; the parent project is not written.
3. **Refinement.** The normal single-state `refine3D_auto` on the copy. Alignment is unweighted:
   every selected particle is aligned against state `X`'s reference. Every reconstruction (the
   iterations, trailing, the final native reconstruction) goes through the service with the
   column-`X` weights. Sampling chooses among the selected rows and never reads the weights; update
   fractions and trailing use applied mass; population gates use effective sample size. Sigma2 is
   bootstrapped for the copy as `refine3D_auto` does for any project.
4. **Results.** The refined map of `X`, its FSC and the refined poses live in the run's work
   project. The registered volume keeps the hard population (`pop`) for reporting and adds the
   applied mass and effective sample size as separate real fields.

**Final native reconstruction.** `calc_final_rec` propagates `m_estimator`, `state` and the
weight-set identity through both its direct and its bootstrap routes.

## 7. What FLEX reuses

| FLEX piece (lines) | SIMPLE target | Lines away | Risk | Phase |
|---|---|---:|---|---|
| State backends, state parts, rec3D driver, delivery, round-weights broadcast (~1,350) | The reconstruction service (4.3) with weights (5); delivery through the standard assembly, FSC, filter and mask | 1,050-1,350 | Delivered maps and inference change (D11): the old filter is a squared 4th-order Butterworth, not the "8th-order" its comment says; the PCG ridge becomes the production one | 3 |
| `vol_flex` publication and the refinement handoff | `vol` and FSC entries; normal discovery; read-only fallback | ~60 | Old projects must load | 3 |
| Weights store (1,265) | `state_weight_set` (4.2) | ~700 net | Low | 1 |
| Opt-in consensus nonuniform filter (`apply_consensus_nu_filter`, `nufilt`) | Removed (D12); FLEX does no nonuniform filtering | ~70 | None: off by default and set by no caller | 1 |
| Clustering: Gaussian mixture, extreme deconvolution, k-means (two copies), diffusion k-center, complete linkage and union-find, farthest-point seeding and log-sum-exp (three copies each) | General modules in `src/utils/clustering` (7.1) | 900-1,300 from FLEX | State order and labels may change (seeding, tie-breaking): judged against truth | 5 |
| Multi-target and coupled insertion, three "one window, many volumes" gathers (~700 in latent ops, polar, E-step) | One batch insertion into K accumulators with per-target weights and a density mode (own, shared, diagonal, packed pairs), and one multi-volume point gather, in `simple_reconstructor.f90`; `project_fplane` becomes the K=1 case; the unused `project_polar` goes | 550-700 (adds ~350 in the reconstructor) | refine3D hot path: keep `insert_plane_oversamp` separate until measured. The coupled insert moves to direct sample locations: rounding-level change | 5 |
| `prep_imgs4projected_model` (~75) | `prep_imgs4rec` with optional band, transfer-plane and observation-model arguments | 70-90 | Optional arguments only on the refine3D path | 5 |
| Coupled PCG operator: geometry, support window, kernel and right-hand-side accumulation, scatter, fold, kernel finalisation (~800 duplicated) | A shared lattice layer composed by `reconstructor_pcg` and the coupled operator; the K×K pair kernels and per-voxel K×K preconditioner stay FLEX | 350-400 each side | Share only code that is identical today | 5 |
| Embedded PCG self-test (~600) | `simple_flex_pcg_tester.f90` | 0 (moved) | None | 5 |
| Sigma2 preflight in the commander (~80) | `canonical_sigma2_consumable` and `simple_sigma2_bootstrap` (keep per-stack grouping), before `gen_job_descr` | ~80 | The global-sigma fallback for an empty half moves with it | 5 |
| Commander and strategy (worker strategy, worker-key injection, thread boost) | A separate worker commander routed by the private executable, `set_master_num_threads`; the rounds callback stays for FLEX's nested EM loops | 110-130 | Thread sizing differs; FLEX stops overwriting `params%nthr` | 5 |
| Plane cache housekeeping | `simple_ptcl_cache`'s cache directory, run token, space budget and key-last publication (D7) | ~50 | None numerical | 5 |
| Embedding cache (184) | One atomic artifact (raw embedding, deconvolution, labels, noise scale) registered in the `out` segment and path-remapped with the project | ~30 | Low; fixes torn resumes | 5 |
| Mean-scale broadcast | `arr2file` | ~45 | Low | 5 |
| Robust median/MAD, power iteration, the SPD Cholesky family | `simple_stat`, `jacobi_dp`, `simple_linalg` | ~150 | Pin eigenvector signs | 5 |
| Cross-fit FSC helpers | `get_find_at_crit`, `fsc2optlp_sub`, `get_resolution_at_fsc`; drop write-only artifact fields | ~100 | The reported band moves (first failing versus last passing shell) | 5 |

### 7.1 Clustering

Every clustering algorithm FLEX uses moves to `src/utils/clustering` as a general tool that follows
the library's clustering standard (section 2): module `simple_<algorithm>` holding type
`<algorithm>`; `new` takes plain arrays (feature matrices, per-point covariances or distance
matrices) and the algorithm's parameters, never a FLEX type, `parameters` or `builder`; one
clustering call returns labels and the representatives (centroids, medoids, component means and
covariances, responsibilities); getters for what callers read; `kill`; seeding and tie-breaking
deterministic for given inputs (a seed is an argument where randomness is needed); labels ordered
by decreasing population, ties by the smallest member index; one tester per module with an
independent oracle.

| Module | Generalizes | Replaces in FLEX |
|---|---|---|
| `simple_hac` (renamed from `simple_avg_linkage`, whose callers move) | linkage average or complete; stop at a cluster count or at a distance threshold; union-find inside | the complete-linkage loop and union-find of the two-gate merge, which keeps its gates and passes `1 - R` as the distance |
| `simple_kmeans` | k-means on feature vectors with optional per-dimension weights, deterministic farthest-point seeding, empty-cluster recovery | `kmeans_latent_targets` and the k-means initialisation of the mixture prior in `simple_flex_pca_posterior.f90` |
| `simple_gmm` | Gaussian mixture by expectation maximization, tied or full covariance, responsibility floor and renormalization, log-sum-exp, the number of components chosen by held-out log-likelihood or by hierarchical splitting | the mixture fits of `simple_flex_pca_gmm.f90`; FLEX keeps its placement decisions (ceiling, merge, how responsibilities become state weights) |
| `simple_xd_gmm` | extreme deconvolution: a Gaussian mixture fitted to points with known per-point covariances, held-out choice of the number of components | the fit, ladder and posterior of `simple_flex_pca_deconv.f90`; FLEX keeps the noise calibration and the projection |
| `simple_kcenter` | farthest-point k-center on feature vectors | the clustering half of `diffusion_kcenter_targets`; the diffusion coordinates come from `simple_diffusion_maps` and `simple_diff_map_graphs` (extension: per-node bandwidth, raw degree) |

Shared numerics they need (log-sum-exp, the SPD Cholesky family, farthest-point seeding where two
modules use it) go to `simple_stat`, `simple_linalg` or `simple_clustering_utils`, once.
`simple_clustering_utils` dispatches `simple_hac` beside k-medoids and affinity propagation in
`cluster_dmat`. Stays in FLEX: path and reliability targets, the mixture prior's coupling to the probe
fit (it calls `simple_gmm` where its update is the same), and every decision of how cluster output
becomes state weights.

Two findings to check in Phase 0: the polar bank snaps particles to bank directions without
symmetry copies, so a non-c1 particle outside the asymmetric unit may snap to a distant direction;
the resident plane store allocates planes on the padded lattice but fills every second sample, so
it holds about four times the memory it uses.

Stays in FLEX, genuinely new: the probe fit and its E-step, the posterior, the coupled solve, the
latent deconvolution's noise model, state placement decisions and the two-gate merge's gates, the
polar bank's variable-length rings and their quadrature contract with the exact low-k part, UMAP,
and the cross-fit component matching. Kept by decision: the probe and embedding codecs (D6), the
plane-cache payload (D7).

## 8. Phases

Each phase is one or a few commits by the maintainer and lands only with its exit met, the
checks of 8.1 included. The measures and tolerances named in the exits are defined in 8.2.

| # | Phase | Scope | Exit |
|---|---|---|---|
| 0 | Preconditions and baselines | `bd3997758` (the base of the run, by the base ruling of section 3; the plan first named `878adf2d4`) builds, and its fast gate, library tier and the high-level entries of 8.1 run; which entries pass is recorded. Record the references of 8.2; the FLEX application phantom (`lib_heterogeneity`) and the workflow gate (`flex_pca_blobs`): per arm the recoded label accuracy, the state count and each truth-matched state map's correlation with its truth; `refine3D_states` on the same data set as the `prob_state` retirement run (best state 4.03 Å there). Check the two findings of section 7 and record what is true; fix nothing. | References and observations in section 12; base build and tier results recorded |
| 1 | Foundations | 4.1, 4.2 (FLEX writes the set; nothing reads it yet), 4.3 with hard labels only (`reconstruct3D` migrated onto the service), D10 in FLEX's pruning, D12 (remove FLEX's nonuniform filter) | `reconstruct3D` same-path equal (8.2) to its Phase 0 references on all four paths; the weight-set tests of section 9 green; FLEX application phantom and workflow gate within the FLEX margins of 8.2 |
| 2 | Fractional reconstruction | Section 5 through the service; `m_estimator` and `state` on `reconstruct3D` with UI and parser | `m_estimator=no` same-path equal to the references; `m_estimator=flex` with 0/1 weights equal to the labels same-path equal on all four paths; fractional shared memory versus two parts within the cross-mode tolerance; the linearity, trailing and fractional-phantom tests of section 9 green |
| 3 | FLEX on the service | Delete FLEX's state backends, parts, delivery and round-weights broadcast; trial and final maps through the service (D9, D11); retire `vol_flex` | Application phantom and workflow gate within the FLEX margins of 8.2 |
| 4 | Per-state M-estimation in `refine3D_auto` | Section 6: `state` and `m_estimator` on `refine3D_auto` with UI and parser, the work project, weighted reconstruction throughout, the final native reconstruction | On the Phase 0 data set: FLEX publishes the weights from the kept single-state project; for the state with the most applied mass, `refine3D_auto state=X m_estimator=flex` reaches an FSC 0.143 resolution at most 1.03 times that of `refine3D_auto state=X` on the same state's hard labels; the per-state refinement phantom of section 9 passes; the parent project is byte-identical after both runs; the final native reconstruction logs and registers the weighted map with applied mass and effective sample size |
| 5 | FLEX reuse | The section 7 rows marked 5, one small diff each, in the order: sigma2 and commander shape; plane-cache housekeeping and the atomic embedding artifact; `prep_imgs4rec`; insertion and gather; PCG lattice layer; the clustering modules of 7.1 (`simple_hac`, `simple_kmeans`, `simple_kcenter`, `simple_gmm`, `simple_xd_gmm`, in that order); statistics, linear algebra and FSC helpers | After each diff: its suites green, the module testers of 7.1 green, and the FLEX A/B of 8.2 against the previous diff |
| 6 | Documents | `doc/algorithms/reconstruction.md` (weighted variant; the density floor no longer FLEX-only), `doc/algorithms/heterogeneity_analysis/flex_pca.md`, `src/main/flex/README.md`, the policies below, the skills | Text matches code; the stale-reference scan of 8.1 is empty |

### 8.1 Checks at the end of every phase

- Debug and Release builds without new compiler warnings; the fast gate green.
- The library tier, and those of the high-level entries `pcg_recon`, `simulate_particles`,
  `simulated_workflow_1jxy`, `simulated_workflow_6vxx`, `solve3D_addon`, `cont_refine3D_1jxy` and
  `flex_pca_blobs` that passed at Phase 0, still pass (the stream and single-particle entries are
  outside this refactor). An entry that failed at Phase 0 is reported, not fixed, unless the phase's
  scope covers it.
- `scripts/check_test_registry.py`, `scripts/check_descr.py` and `scripts/check_flex_dag.py` pass
  (the checker's module table follows modules a phase adds or deletes); no file-mode changes.
- No `THROW_HARD`/`THROW_WARN` argument continued with `//&`; no `flag .and. array(i)` on one
  line where the array may be unallocated.
- Phase 6 only: no living document, skill or README names a deleted file, module, routine or key.

### 8.2 References, measures and tolerances

- **References (Phase 0).** Hard `reconstruct3D` on four paths (gridding and PCG, each in shared
  memory and with two parts), each run twice, from the three-state beta-galactosidase project that
  Phase 0's `refine3D_states` run produces (hard labels, state 0 present).
- **Same-path equal.** Bit for bit, on a path whose two Phase 0 runs were bitwise identical. On a
  path whose two runs differed, the relative L2 difference of every restored half map is at most
  ten times the larger Phase 0 run-to-run difference; Phase 0 records which paths are
  deterministic and their run-to-run differences.
- **Cross-mode tolerance (fractional).** Shared memory against two parts: relative L2 difference of
  each state's even and odd half maps at most 1e-5, and their correlation at least 0.99999.
- **FLEX margins.** Per arm of the application phantom and per mode of the workflow gate: the
  recoded label accuracy at least the Phase 0 value minus 0.02 and at least 0.95; the delivered
  state count between 2 and the provision ceiling; each truth-matched state map's correlation with
  its truth (same frame, low-pass 8 Å) at least the Phase 0 value minus 0.02; each matched map
  prefers its own truth; in the workflow gate, shared and distributed recoded labels agree on at
  least 98 % of the particles.
- **FLEX A/B (Phase 5).** On the application phantom, before and after the diff. For a diff that
  does not touch clustering: identical recoded hard labels, and the latent embedding equal bit for
  bit where the section 7 row says bitwise, otherwise within a maximum relative difference of 1e-4
  after matching component order and sign. For a clustering diff, whose seeding and tie-breaking
  may renumber or move states: the FLEX margins instead, with the embedding equal as above.

Policies to read before Phase 1 and amend in Phase 6:
`doc/policies/3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`
(invariant 3, "every valid selected particle belongs to exactly one (state, half) group";
invariant 4 stands), `doc/policies/3D/refine3D_policy.md` (hard assignment only; the per-state
M-estimation of section 6), `doc/policies/importance_sampling_fractional_update_policy.md` (row
counts per state), `doc/policies/heterogeneity/refine3D_states_policy.md`,
`doc/policies/3D/reconstruct3D_pcg_policy.md`, `doc/policies/2D/particle_cache_policy.md` (shared
cache machinery).

## 9. Test net

Independent oracles, in the fast or library tier unless stated:

- **Same path, 0/1 weights equal hard labels.** Gridding and PCG, shared memory and two parts,
  each same-path equal (8.2) to its own hard reference.
- **Cross-mode, fractional.** Shared memory against two parts within the cross-mode tolerance of
  8.2 (the reduction trees differ, so bitwise is not expected).
- **Linearity.** Inserting a particle with weight `w` equals inserting it with weight 1 and scaling
  the accumulators by `w`, to round-off.
- **Fractional phantom** (library tier): particles simulated from two conformers with known
  poses, given fractional weights that need not match their conformer. The expected map of state
  `s` is the mixture of the truth maps with the fractions of state `s`'s applied mass that come
  from each conformer. Each weighted state map's correlation with its expected mixture (same frame,
  low-pass 8 Å) is at least the correlation of the hard reconstruction of the same particles with
  its truth, minus 0.01.
- **Per-state refinement phantom** (high-level tier, Phase 4): FLEX on the two-conformer phantom
  publishes the weights; for each state, `refine3D_auto state=X m_estimator=flex` and
  `refine3D_auto state=X` (hard labels) refine their work projects; each weighted map correlates with
  its truth at least as well as the hard one, minus 0.01; the parent project is unchanged.
- **Weight set.** Round trip exact; digest, count, size or checksum mismatch refused; a crash
  simulated between file and manifest publication leaves the previous generation selected; an old
  project with `flex_weights` entries loads; pruning leaves `state > 0` and the set's selection in
  agreement; a work project's single-state set reads back the parent's column for the selected rows.
- **Trailing.** A chain built under one weight-set identity is re-seeded, never blended, under
  another; PCG and gridding chains report the same applied mass for the same input.
- **Clustering modules** (fast tier, one tester each): well-separated synthetic clusters recovered
  exactly, with labels ordered by population; a mixture with known parameters recovered within its
  sampling error; extreme deconvolution recovers the noise-free mixture from points with known
  covariances; average and complete linkage against hand-computed merges; determinism (two calls,
  same labels).
- **Regression:** the Phase 0 `refine3D_states` run and the FLEX application phantom, against the
  recorded observations with stated margins.

## 10. Risks

- **Cost.** Gridding re-reads a particle once per state it weighs into: 1-2× for sparse
  responsibilities, `nstates`× for flat ones. PCG re-reads per `(state, half)`, as FLEX does today.
  A per-state refinement refines the particles of one state; refining every state costs one
  `refine3D_auto` run per state.
- **Inference moves in Phase 3 and in the clustering diffs of Phase 5.** Raw trial maps (D11) and
  new seeding or tie-breaking change bandwidths, merges, labels and possibly the state count. Those
  exits are against truth.
- **refine3D hot path.** Phases 1, 2, 4 and 5 touch `insert_plane_oversamp`, `calc_3Drec`,
  `prep_imgs4rec`, the `reconstruct3D` strategy and `refine3D_auto`. Every change is optional or
  gated, and the same-path oracle runs after each.
- **Scope creep.** Recomputed weights, per-state poses inside one run and a generic part-file
  facility are out of this plan.
- **Plan drift.** Earlier rounds rewrote their own plan to match what was delivered. Here the plan
  is fixed: a phase that cannot meet its exit stops and reports.

Out of scope, noted: `maybe_postprocess_reconstruct3D` copies a whole `parameters` per state
(`params_pp = params` in `simple_rec3D_strategy.f90`), which `doc/policies/compile_time_policy.md`
forbids; the service extraction should remove it.

## 11. File table

One row per file or folder; a folder (trailing `/`) covers every file under it. Unit testers
(`*_tester.f90`) of a listed file are covered by its row. Before editing a file no row lists, add
a row with the phase and a one-line reason.

| File | Change | Phases |
|---|---|---|
| `doc/refactoring_notes/planned/flex_on_simple_machinery.md` | this plan: progress, added rows | each |
| `doc/refactoring_notes/completed/flex_on_simple_machinery.md` | the plan, moved when done, with its report beside it | 6 |
| `doc/code_overview/` | indexes and maps the build regenerates | each |
| `production/CMakeLists.txt` | CTest registration of new tests | 1 2 3 4 5 |
| `src/main/commanders/test/`, `src/main/ui/simple_test/`, `src/main/exec/simple_test_exec_highlevel.f90`, `src/main/apis/simple_test_exec_api.f90` | test commanders, gates, registration and test UI | 1 2 3 4 5 |
| `scripts/check_flex_dag.py` | module table follows modules added or deleted | 1 3 5 |
| `src/main/project/simple_ptcl_layout.f90` | new: the particle-layout digest (4.1) | 1 |
| `src/main/sigma2/simple_sigma2_state.f90`, `src/fileio/simple_sigma2_files.f90`, `src/fileio/simple_projfile_utils.f90`, `src/main/commanders/simple/simple_commanders_euclid.f90`, `src/main/commanders/simple/simple_commanders_euclid_distr.f90`, `src/main/commanders/simple/simple_commanders_solve2D.f90`, `src/main/solve/simple_solve3D_manifest.f90`, `src/main/solve/simple_solve3D_utils.f90` | callers of the renamed digest | 1 |
| `src/main/project/simple_state_weight_set.f90` | new: `state_weight_set` (4.2) | 1 2 4 |
| `src/fileio/simple_state_weights_file.f90` | new: the weight-set file codec (4.2) | 1 4 |
| `src/fileio/simple_flex_weights_file.f90`, `src/main/flex/states/simple_flex_weights_state.f90` | deleted, replaced by the two above | 1 |
| `src/main/project/simple_sp_project_out.f90`, `src/main/project/simple_sp_project.f90` | `state_weights` entry; ignored `flex_weights` entries; `vol_flex` read-only fallback; mass and effective size fields | 1 3 4 |
| `src/main/strategies/parallelization/simple_rec3D_service.f90` | new: the reconstruction service (4.3) | 1 2 3 4 |
| `src/main/strategies/parallelization/simple_rec3D_strategy.f90`, `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`, `src/main/commanders/simple/simple_commanders_rec_distr.f90`, `src/main/commanders/simple/simple_commanders_rec.f90` | `reconstruct3D` onto the service; weighted PCG workers and master; trailing counts; `state` and `m_estimator` | 1 2 |
| `src/main/ui/simple/simple_ui_reconstruct3D.f90` | `state` and `m_estimator` inputs | 2 |
| `src/main/ui/simple/simple_ui_refine3D.f90` | `refine3D_auto`: `state` and `m_estimator` inputs | 4 |
| `src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_parse.f90` | `nufilt` removed (D12); `m_estimator` added; `state` for `refine3D_auto` | 1 2 4 |
| `src/main/params/simple_parameters_phases.f90` | derived flags of the new keys: `m_estimator` validated, whether `state` was given | 2 4 |
| `src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_parse.f90`, `src/main/params/simple_parameters_phases.f90`, `src/main/ui/simple/simple_ui_heterogeneity.f90` | FLEX's `box_rec` and `smpd_rec` removed: its maps follow the service's box rules (D9) | 3 |
| `src/main/strategies/search/simple_matcher_3Drec.f90`, `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`, `src/main/commanders/simple/simple_commanders_rec_distr.f90` | the caller-owned weight table of the service (FLEX's trial maps): membership, weights and the density floor | 3 |
| `src/main/nu_filt/simple_nu_filter.f90` | the comment that names FLEX's `nufilt` report | 1 |
| `src/main/volume/simple_reconstructor.f90` | weighted insertion (5); batch insertion and multi-volume gather (section 7) | 2 5 |
| `src/main/volume/simple_reconstructor_pcg.f90` | real mass and weight-set identity in raw accumulators; shared lattice layer | 2 5 |
| `src/main/volume/simple_trail_chain_manifest.f90`, `src/main/volume/simple_frozen_accum.f90` | the manifest contract for both backends; provenance in frozen accumulators | 2 |
| `src/main/strategies/search/simple_matcher_3Drec.f90` | weighted membership lists; `prep_imgs4rec` options | 2 5 |
| `src/main/ori/simple_oris_sampling.f90`, `src/main/ori/simple_oris_getters.f90` | population rule and update counts on applied mass | 2 |
| `src/main/ori/simple_oris.f90` | interfaces of the population rule on real masses and of the row lists the weighted counts reduce over | 2 |
| `src/defs/simple_refine3D_fnames.f90` | name of the PCG trailing chain's manifest | 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`, `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | population gates on effective sample size; weight set threaded through assembly | 2 4 |
| `src/main/strategies/search/simple_strategy3D_matcher.f90` | reconstruction through the service with the weight set | 2 4 |
| `src/main/commanders/simple/simple_commanders_refine3D.f90` | `vol_flex` handoff replaced by normal discovery; `refine3D_auto state=X` work project and `m_estimator` | 3 4 |
| `src/main/commanders/simple/simple_commanders_project_core.f90` | selection of a work project's rows by weight (beside selection by state) | 4 |
| `src/main/simple_final_rec.f90` | weights through the direct and bootstrap routes; mass and effective size registered | 4 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | the state continuation drops `state=` once its work project holds that state alone (since Phase 2, `state=` restricts reconstruction) | 4 |
| `src/main/flex/` | FLEX (D10, D12, Phase 3 deletions, Phase 5 reuse) | 1 3 5 |
| `src/main/commanders/simple/simple_commanders_flex_pca.f90`, `src/main/strategies/parallelization/simple_flex_pca_strategy.f90`, `production/simple_private_exec_driver.f90` | FLEX writes the weight set; worker commander and standard shape | 1 3 5 |
| `src/utils/clustering/` | the clustering modules of 7.1; `simple_avg_linkage` renamed to `simple_hac` | 5 |
| `src/main/apis/simple_private_exec_api.f90` | the private executable's import of FLEX's worker commander | 5 |
| `src/main/strategies/search/simple_ptcl_cache.f90`, `src/main/sigma2/simple_sigma2_bootstrap.f90`, `src/main/pca/simple_diff_map_graphs.f90`, `src/main/pca/simple_diffusion_maps.f90`, `src/utils/math/` | shared machinery FLEX reuses (section 7) | 5 |
| `doc/algorithms/reconstruction.md`, `doc/algorithms/heterogeneity_analysis/`, `doc/policies/3D/`, `doc/policies/heterogeneity/`, `doc/policies/importance_sampling_fractional_update_policy.md`, `doc/policies/2D/particle_cache_policy.md`, `doc/how2s/how2_process_heterogeneous_datasets.md` | documents and policies | 6 |
| `src/main/flex/README.md` | the FLEX layout after the phases | 6 |
| `.github/skills/` | skills that describe the changed contracts | 6 |

## 12. Progress

### Phase 0: preconditions and baselines (2026-10-07, done)

No source file changed. The only edited file is this plan: the base ruling recorded in section 3
and in the Phase 0 row of section 8, and this entry. Paths below are under
`/home/elmlundho/agent_runs/flex_on_simple/scratch` unless absolute.

**Base and builds.** The base is commit `bd3997758` with a clean tree (the base ruling of section 3;
the copy's launch status lists no local differences). The compiler is GNU Fortran 15.2.1
(gcc-toolset-15), with CMake 3.26.5 and Python 3.12 for the build scripts. The configuration matches
the maintainer's `~/src/SIMPLE/build`: OpenMP on, `USE_ARCHOPT=ON`, no MPI. The login shell points
`SIMPLE_PATH` and `PATH` at that other checkout, so every script here strips those entries and sets
both to the build it uses.
- Release: out of tree in `build_release` (`build_release.sh`: `CMAKE_BUILD_TYPE=Release`,
  `BUILD_TESTS=ON`, `USE_ARCHOPT=ON`, installed into the build tree). Warning baseline: no compiler
  warnings; five linker warnings "requires executable stack" (`json_value_module.f90.o` three times,
  `simple_commanders_cavgs.f90.o`, `simple_stream_stage_pool2D.f90.o`).
- Debug: `./compile_debug.sh` in the repository (`build_debug.sh`; output `logs/debug_build.out`).
  Warning baseline: no compiler warnings; six linker warnings, the five above plus
  `simple_persistent_worker.f90.o`. The fast gate passed 15 of 15 in 6.6 s and the build installed.

**Cases.** These ran one at a time, in the background, on the otherwise idle 24-core machine. The
load average was 6.8 at the start of S, from the tail of the builds; after that the load came from
the cases themselves. During S, the two single-threaded check programs of the findings below ran
for a few seconds; nothing else ran. The execution settings
were `nparts=4 nthr=6`, as in the earlier runs on this machine. Every case starts from
`/data2/bgal_addon/control/bgal_addon.simple`: beta-galactosidase after 2D classification, 5,513
particles, box 256 at 1.275 Å per pixel, 50 selected class averages.
- Scripts: `run_case.sh` (the exact command lines), `queue_phase0.sh` (the order),
  `case_summary.sh`.
- Logs per case, in `logs/`: `case_<name>.log` (each line prefixed by its time in seconds,
  `/usr/bin/time -v` at the end), `summary_<name>.txt` and `table_<name>.md`.
- A project's output segment refers to its maps by absolute path inside the run directory that
  wrote them. Before a case that reads a kept project, `run_case.sh` therefore recreates
  `runs/S/1_solve3D_cavgs` (and, for B, `runs/single/1_solve3D`) from `keep/`, so a case depends only
  on `keep/`.
- Kept for later phases: `keep/S/1_solve3D_cavgs`, `keep/single/1_solve3D` (Phase 4 starts from
  it), and `keep/B` (the whole run directory of B). B's result project is
  `keep/B/1_refine3D_states/bgal.simple`, because `refine3D_states` writes into `1_refine3D_states/`;
  `keep/B/bgal.simple` is only the input copy.

| Case | Starts from | Command line (in `runs/<case>`, project copy `bgal.simple`) |
|---|---|---|
| S | fresh copy of the input project | `simple_exec prg=solve3D_cavgs projfile=bgal.simple mskdiam=180 pgrp=d2 nparts=4 nthr=6` |
| single | `keep/S/1_solve3D_cavgs/bgal.simple` | `simple_exec prg=solve3D projfile=bgal.simple mskdiam=180 pgrp=d2 nsample=2500 cavg_ini_ext=yes nparts=4 nthr=6` |
| B | `keep/single/1_solve3D/bgal.simple` | `simple_exec prg=refine3D_states projfile=bgal.simple nstates=3 nsample=1500 lpstart=12 lpstop=6 pgrp=d2 mskdiam=180 nparts=4 nthr=6` |

Population is the number of particles carrying that state label in the final project; state 0 means
excluded. Resolutions are the Fourier shell correlation (FSC) at the 0.5 and 0.143 thresholds of
each state's final reconstruction. For S the two half sets are the even and odd class averages.
Every case exited with status 0.

| Case | Wall time | Peak memory | State: population, FSC 0.5 / 0.143 (Å) | `prob_state` retirement run, Phase 0 |
|---|---|---|---|---|
| S | 10:58 | 3.1 GB | 1: 5,513 (50 class averages), 7.77 / 6.80 | 10:56, 7.77 / 6.80 |
| single | 15:44 | 3.2 GB | 1: 5,513, 4.41 / 3.89 | 14:38, 4.35 / 3.84 |
| B | 19:27 | 4.8 GB | 1: 3,052, 5.44 / 4.03; 2: 738, 14.84 / 12.09; 0: 1,723 | 18:55; 1: 3,460, 5.02 / 4.03; 2: 726, 14.19 / 12.09; 0: 1,327 |

What B did, read from its log:
- FLEX (`flex_pca`) initialized the states. Its population floor delivered 3 states, then its
  occupancy floor, which requires at least 2,000 effective particles per state, dropped one of them.
  The 1,723 particles of that state (31.25 %) went to state 0, and the run continued with
  "REFINE3D_STATES FLEX_PCA INITIALIZED NSTATES: 2".
- It ran 32 `prob_neigh` iterations, with sampling "ACTIVE PARTICLES/SAMPLE TARGET/UPDATE_FRAC:
  3790/1500/0.3958". On this base, unlike the earlier run, the active count excludes state 0.
- So the reference project has two populated states plus state 0, not three: hard labels, with
  state 0 present, as 8.2 requires. The best state reaches 4.03 Å, the value the `prob_state`
  retirement run recorded.

**References (8.2).** These are hard `reconstruct3D` runs on B's result project, from fresh copies,
with the gridding backend and with the preconditioned conjugate gradients (PCG) backend.
- Scripts: `run_ref.sh` (one run: it recreates `runs/S/1_solve3D_cavgs`, `runs/single/1_solve3D`
  and `runs/B` from `keep/`, then copies `keep/B/1_refine3D_states/bgal.simple` to
  `runs/ref_<path>_<rep>/bgal.simple`), `queue_refs.sh`, and `cmp_maps.py` (the comparison).
- Logs: `logs/ref_<path>_<rep>.log`, `logs/queue_refs.log`, `logs/cmp_ref_<path>.txt`.
- The first run of each path is kept under `keep/refs/<path>` (run directory with
  `1_reconstruct3D/`).
- The `reconstruct3D` user interface lists `rec_backend`, `nparts`, `nthr`, `mskdiam`, `pgrp` and
  others, but not `nstates`. Without `nstates=2` the program reconstructs state 1 only: it skips the
  738 state-2 particles and prints "your input doc has multiple states but NSTATES is not given"
  (first attempt kept as `logs/ref_grid_shm_1_no_nstates.log`). The parser accepts `nstates`, so
  the reference lines carry `nstates=2`.

Later phases repeat these lines unchanged, in `runs/ref_<path>_<rep>`:

| Path | Command line |
|---|---|
| gridding, shared memory | `simple_exec prg=reconstruct3D projfile=bgal.simple mskdiam=180 pgrp=d2 nstates=2 rec_backend=gridding nthr=24` |
| gridding, two parts | `simple_exec prg=reconstruct3D projfile=bgal.simple mskdiam=180 pgrp=d2 nstates=2 rec_backend=gridding nparts=2 nthr=12` |
| PCG, shared memory | `simple_exec prg=reconstruct3D projfile=bgal.simple mskdiam=180 pgrp=d2 nstates=2 rec_backend=pcg nthr=24` |
| PCG, two parts | `simple_exec prg=reconstruct3D projfile=bgal.simple mskdiam=180 pgrp=d2 nstates=2 rec_backend=pcg nparts=2 nthr=12` |

Results:
- Each path ran twice, and all four paths are deterministic. Every `.mrc` file the two runs of a
  path wrote is byte-identical: the merged, even and odd maps of both states, their unfiltered,
  low-pass and post-processed variants, and for gridding the partial and density files (24 files
  for gridding in shared memory, 32 for gridding with two parts, 16 for each PCG path). Same-path
  equality in later phases is therefore bit for bit on all four paths.
- Wall times: gridding 14 to 16 s, PCG 2:20 to 2:29.
- Resolutions on every path: state 1 FSC 0.143 4.03 Å, state 2 12.09 Å. At FSC 0.5, state 1 is
  5.44 Å with gridding and 5.10 Å with PCG; state 2 is 14.84 Å on all paths.
- For context (no tolerance attached), shared memory against two parts: the relative L2 difference
  of the restored maps is at most 7.3e-6 for gridding and 7.8e-7 for PCG
  (`logs/cmp_ref_shm_vs_np2.txt`).
- Two mistakes of the script were caught and redone. The first queue started on B's input copy:
  it was stopped after 11 s and its output deleted. The first comparison script failed on the
  complex partial files after the queue had already deleted the second runs: the second runs of the
  two gridding paths were repeated (`redo_grid_rep2.sh`) and compared against the kept first runs.

**FLEX observations.**
- Application phantom (`lib_heterogeneity`, Release build): FLEX on a simulated two-conformer data
  set of 2,000 particles (box 64; conformers A and B differ by one lobe moved 8 pixels), with the
  requested state ceiling `npreimages=3`. There are four arms, the basis backend times the
  state-map backend, gridding or preconditioned conjugate gradients (PCG). The table gives, per
  arm, the delivered states, the recoded label accuracy (fraction of particles whose state maps to
  their true conformer), `pc_corr` (the correlation of FLEX's leading eigenvolume with the true
  difference map A minus B) and the smallest correlation of a truth-matched state map with
  its truth. From `logs/tiers/LastTest_lib_heterogeneity.log`, after "FLEX_PCA TWO-STATE PHANTOM
  PHASE-0 OBSERVATIONS":

  | Arm | States | Label accuracy | pc_corr | Min matched corr | Seconds |
  |---|---:|---:|---:|---:|---:|
  | grid_basis_grid | 3 | 0.9995 | 0.9605 | 0.9846 | 1.9 |
  | pcg_basis_grid | 3 | 0.9995 | 0.9647 | 0.9845 | 5.1 |
  | grid_basis_pcg | 3 | 0.9995 | 0.9605 | 0.9976 | 3.3 |
  | pcg_basis_pcg | 3 | 0.9995 | 0.9647 | 0.9976 | 6.5 |

  Per-state truth counts (particles of conformer A and B in each state) and correlations of each
  state map with the two truth maps:

  | Arm | State: nA, nB, corr A, corr B | Matched |
  |---|---|---|
  | grid_basis_grid | 1: 1, 1000, 0.9310, 0.9846; 2: 3, 0, 0.9681, 0.9716; 3: 996, 0, 0.9851, 0.9354 | A = state 3, B = state 1 |
  | pcg_basis_grid | 1: 1, 1000, 0.9311, 0.9846; 2: 3, 0, 0.9677, 0.9720; 3: 996, 0, 0.9851, 0.9353 | A = 3, B = 1 |
  | grid_basis_pcg | 1: 1, 1000, 0.9469, 0.9979; 2: 3, 0, 0.9825, 0.9846; 3: 996, 0, 0.9976, 0.9463 | A = 3, B = 1 |
  | pcg_basis_pcg | 1: 1, 1000, 0.9469, 0.9979; 2: 3, 0, 0.9817, 0.9848; 3: 996, 0, 0.9976, 0.9461 | A = 3, B = 1 |

  Each matched map prefers its own truth. The unmatched state 2 holds only three particles, all of
  conformer A.
- Workflow gate (`flex_pca_blobs`, Release build): the same phantom, run through the public
  `flex_pca` program with PCG for both the basis and the state maps, once in shared memory and once
  distributed over two worker processes (`nparts=2`).
  - **It fails on the base.** The shared-memory mode passes every check. The distributed mode stops
    in the master before any worker starts, with "oris object does not exist; get_ori"
    (`logs/tiers/flex_pca_blobs.out`; Debug backtrace in `logs/flex_pca_blobs_debug.log`). The
    chain is `master_plan_partitions` in `simple_flex_pca_strategy.f90`, then `qsys_env%new`
    (`simple_qsys_env.f90`), which reads the project's `compenv` segment (the queue settings, which
    `new_project` writes). The project of the gate, written by `build_flex_pca_phantom_project` in
    `simple_flex_pca_application_tester.f90`, has no such segment.
  - Nothing on this path changed between `878adf2d4` and the base, so the failure is deterministic
    on the base. `doc/refactoring_notes/completed/flex_refactoring_plan_2026_10_06.md`
    records the gate as green in the maintainer's own runs (its phase ledger, phases 0 and 1).
  - Under 8.1 an entry that fails at Phase 0 is reported, not fixed; nothing in the repository was
    changed.
  - To give the later FLEX margins a distributed baseline anyway, the gate was also run once on an
    evidence build outside the repository: the base from `git archive bd3997758`, plus a three-line
    change that makes the phantom fixture call `project%update_compenv` before writing its project
    (`logs/evidence_patch.diff`). This build is not part of the run's diff and was deleted after
    the run. The gate passes there in 22.8 s (`logs/tiers/LastTest_flex_pca_blobs_evidence.log`).
    Its `metrics.tsv` rows, as the gate logs them:

  | Row | Shared | Distributed |
  |---|---:|---:|
  | truth_label_accuracy (floor 0.6) | 0.99950 | 0.99750 |
  | leading_eigenvolume_truth_correlation | 0.96468 | 0.92030 |
  | deconvolution_noise_scale | 0.75683 | 0.77554 |
  | delivered states (ceiling 3) | 3 | 3 |
  | State: nA, nB, corr A, corr B | 1: 1, 1000, 0.9469, 0.9979; 2: 3, 0, 0.9817, 0.9848; 3: 996, 0, 0.9976, 0.9461 | 1: 0, 983, 0.9467, 0.9979; 2: 5, 17, 0.9876, 0.9791; 3: 995, 0, 0.9975, 0.9480 |

  Between the two modes: published_state_count_difference 0, shared_distributed_truth_label_agreement
  0.99800 (the 98 % margin of 8.2 is met on the evidence build), shared_distributed_latent_correlation
  0.56132, and mode-A and mode-B map correlations 0.99995 and 0.99998. The gate's other checks pass
  in both modes: the weight set validates, both truth modes are retained, labels agree with the
  weights, the eigenvolume and embedding are published, and each matched map prefers its own truth.
  The shared-memory rows of the failing base run are identical to the shared rows above.
  - Consequence for later phases: on the base itself the distributed mode of the workflow gate has
    no Phase 0 value, so the exits of Phases 1 and 3 can hold the gate's distributed mode only
    against the evidence-build numbers above, after the fixture writes a `compenv` segment. The
    fixture lives in `src/main/flex/`, which the file table assigns to Phases 1, 3 and 5. The
    phase that first needs the distributed gate must decide whether to make that fixture change and
    record it.

**Tiers (Release build, one entry at a time; `queue_tiers.sh`, `logs/queue_tiers.log`, outputs and
CTest logs in `logs/tiers/`).**
The commands: `scripts/run_fast_gate.sh build_release` for the fast gate, then
`ctest -R '^<entry>$'` per entry, each with `SIMPLE_PATH` set to the Release build.

| Entry | Label | Result | Wall time |
|---|---|---|---:|
| fast gate (15 suites) | fast | pass, 15 of 15 (Debug build: 15 of 15 too) | 3 s |
| `lib_reconstruction` | library | pass | 8 s |
| `lib_cart_align3D` | library | pass | 25 s |
| `lib_heterogeneity` | library | pass | 23 s |
| `lib_single` | library | pass | 14 s |
| `lib_stream` | library | pass | 617 s |
| `pcg_recon` | highlevel | pass | 1 s |
| `simulate_particles` | highlevel | pass | 8 s |
| `simulated_workflow_1jxy` | highlevel | pass | 256 s |
| `simulated_workflow_6vxx` | highlevel | failed, then passed on one rerun (flaky) | 690 s; rerun 786 s |
| `solve3D_addon` | highlevel | pass | 853 s |
| `cont_refine3D_1jxy` | highlevel | pass | 545 s |
| `flex_pca_blobs` | highlevel | **fails** (distributed mode; see the FLEX observations) | 9 s |

- `simulated_workflow_6vxx`: the first run failed its final-map validation. Its de novo solution
  correlated 0.5420 with the truth (minimum 0.80), with the truth FSC 0.143 at 49.92 Å.
- The rerun passed, with correlation 0.8821 and truth FSC 0.143 at 5.43 Å
  (`logs/tiers/LastTest_simulated_workflow_6vxx_rerun.log`). The rerun overlapped the evidence build
  of the next paragraph, so its wall time is not comparable.
- `doc/refactoring_notes/completed/trailing_reconstruction_without_halfmap_blend_report.md` records
  the same flakiness on its unchanged base. Later phases should treat one failure of this entry
  that passes on a rerun as this known flakiness, not as a regression.
- `scripts/check_test_registry.py`, `scripts/check_descr.py` and `scripts/check_flex_dag.py` pass, and
  no file mode changed.

**The two findings of section 7.** Checked with two small programs in `findings/`, linked against
the Release library (`build_check.sh`). Nothing was fixed.
1. *Polar-bank snap: true, but inert for particles SIMPLE refined.*
   - How it works: the E-step builds the bank of shared directions with
     `pgrpsyms%build_refspiral`, which covers the point group's asymmetric unit with its mirror
     (`simple_flex_probe_fit_estep.f90`, `cov_polar_enabled`). It then assigns every particle the
     bank direction whose plane normal has the largest dot product with the particle's raw normal
     (`polar_assign_directions` in `simple_flex_pca_polar.f90`). No symmetry operator and no
     `rot_to_asym` is applied first.
   - `check_polar_snap.f90` (output `check_polar_snap.out`), point group d2 with 1,000 directions:
     the orientation (30°, 40°, 0°) and its three symmetry mates snap with errors of 1.35°,
     19.63°, 59.67° and 50.14°, and all four snap with 1.35° after `rot_to_asym`.
   - Over 20,000 random orientations the raw snap error has mean 31.7°, median 27.3° and maximum
     93.5°, with 71 % above 5°. After `rot_to_asym` the mean is 1.25° and the maximum 3.5°. The c1
     control has a maximum of 5.1°.
   - On the real d2 project `keep/single` (`check_polar_snap_project.out`), all 5,513 particles snap
     within 3.05° (mean 1.20°), identical with and without `rot_to_asym`. Particles that SIMPLE
     refined already lie in the bank's coverage. Imported or symmetry-expanded orientations would
     not.
   - Only the polar rings, the high-frequency part above the hybrid radius, use the snapped
     direction. The mean projection and the exact low-frequency part use the particle's true pose.
2. *Resident plane store: true, about four times.*
   - How it works: the store (`plane_store_store` in `simple_flex_pca_planes.f90`) keeps whole
     `fplane_type` planes from `gen_fplane4rec` (`simple_image_ctf.f90`). That routine allocates
     h in [-B, B] and k in [-B, 0] for B = `box_crop` (the padded box is 2B), then fills only every
     second sample in both directions (stride `OSMPL_PAD_FAC` = 2). With the observation model, a
     lattice point costs 20 bytes (complex data, real squared transfer, complex transfer).
   - `check_plane_bytes.f90` (output `check_plane_bytes.out`) measured, for `box_crop` 256: 131,841
     entries, 2,636,820 bytes per plane, of which 33,153 entries are on the written lattice (ratio
     3.98) and 25,973 are nonzero inside Nyquist (ratio 5.08). For `box_crop` 128: 663,060 bytes,
     ratio 3.95.
   - The store is enabled only with `cache=yes`. The over-allocation itself is in
     `gen_fplane4rec`, which every caller shares, not in FLEX.

**Exit items.**
- The base builds, Debug and Release, with the warning baselines above. Its fast gate, library tier
  and the high-level entries of 8.1 ran, and which pass is recorded above: everything passes except
  `flex_pca_blobs`, and `simulated_workflow_6vxx` is flaky.
- References of 8.2: recorded, with the command lines; all four paths are deterministic.
- FLEX application phantom and workflow gate: the phantom is recorded per arm. The workflow gate
  fails on the base in its distributed mode (reported under 8.1). Its per-mode values are recorded
  from an evidence build outside the repository that differs only in the test fixture.
- `refine3D_states` on the `prob_state` retirement data set: recorded (best state 4.03 Å, as
  there).
- The two findings: checked and recorded; nothing fixed.

Bulky outputs were deleted after their numbers were recorded: the recreated run directories,
the evidence build and the test outputs in the build tree. `scratch/keep/` holds S, single, B and
the references.

Files changed in Phase 0: `doc/refactoring_notes/planned/flex_on_simple_machinery.md` only.

### Phase 1: foundations (2026-10-07, done)

Paths below are under `/home/elmlundho/agent_runs/flex_on_simple/scratch` unless absolute. The
builds are those of Phase 0: Release out of tree in `build_release`, and Debug with
`./compile_debug.sh` in `repo/build`. `inc_debug.sh` is an incremental Debug build followed by
the fast gate.

**4.1, the particle-layout identity.** The digest moved from the sigma2 store to a module of its
own, `src/main/project/simple_ptcl_layout.f90`. Its generic name `ptcl_layout_digest` has two
specifics:
- from a project and its particle field: the former `sigma2_state_project_layout_digest`;
- from plain arrays: the former `sigma2_state_layout_digest`, moved with it so the new module
  does not depend on the sigma2 store.

The bodies are unchanged, so the digest takes the same values. The callers switched: the sigma2
files and store tester, project-file utilities, the Euclidean-objective commanders, `solve2D`, the
`solve3D` manifest and utilities, the `reconstruct3D` strategy, the FLEX commander and the FLEX
weight store (deleted in 4.2).

One correction to the plan's text: the gridding trailing-chain manifest
(`simple_trail_chain_manifest.f90`) holds no layout digest, so it is not a caller.

After this step the Debug build was clean and the fast gate passed 15 of 15
(`logs/inc_debug_p1_41.log`).

**4.2, the state weight set.**
- *Codec:* `src/fileio/simple_state_weights_file.f90`.
  - One binary file per state, `state_weights_gGGGGGG_sNNN.bin` (six-digit generation, three-digit state). It holds a 16-byte magic and
    eight 64-bit header words (version, particle count, state count, state index, generation,
    layout digest, kind, byte size), then the state's weight (32-bit real) for every project row,
    then a hard-label flag (8-bit integer) per row.
  - One text manifest, `state_weights.txt`, written under a temporary name, flushed and renamed
    atomically. It records the format and version, generation, kind (`partition` or `kernel`),
    producer, state count, particle count and layout digest. Per state it records the file name
    (relative to the manifest), byte size, whole-file checksum (the 64-bit FNV-1a, Fowler–Noll–Vo,
    digest the sigma2 store uses), applied mass, effective sample size and hard population.
- *Object:* `state_weight_set` in `src/main/project/simple_state_weight_set.f90`.
  - `new` opens the set the project registers and validates all of it once: the manifest, each
    file's existence, size and checksum, and each file header against the manifest. It also checks
    the layout digest and the particle count against the current project, values in [0,1], at most
    one label per row, labels only on positive weights, the populations recomputed from the files,
    and for PARTITION sets the rows summing to one. Without a status argument it stops the program
    on any failure.
  - One state column is held at a time.
  - Getters: kind, producer, generation and layout digest (the set identity), and per state the
    mass, effective sample size and hard population. Readers: raw weights, applied weights (values
    at or below a consumer's threshold set to zero), hard labels, the selection (rows whose applied
    weights sum above zero), and applied mass and effective sample size over a row subset.
  - `publish` writes a new generation in the working directory, in the order of section 4.2: state
    files, then the manifest, then the project pointer (the `out` segment entry
    `imgkind=state_weights`, written to the project file), then deletion of the superseded
    generation's files in the same directory. The kind is inferred: PARTITION when every weighted
    row sums to one within 1e-3, otherwise KERNEL. The producer must be one word, because the
    manifest stores it as a single token.
  - The update counts of section 4.2 and the single-state set for a work project are left to
    Phases 2 and 4, where they are used.
- *Project:* `add_state_weights2os_out` and `get_state_weights` replace the per-state `flex_weights`
  helpers. `out` entries of earlier releases with `imgkind=flex_weights` still load and are ignored.
- *FLEX:* `write_state_weight_set` in the project gateway replaces `write_flex_weights_store`. The
  producer is `flex_pca`, or `flex_pca_merged` after the automatic merge. The FLEX-only per-state
  latent targets and bandwidths are no longer stored with the weights: nothing read them.
  `simple_flex_weights_state.f90` and `simple_flex_weights_file.f90` are deleted, and
  `scripts/check_flex_dag.py` no longer lists them. The FLEX testers read the set through a helper
  of the application tester (`read_state_weight_table`).
- *Tests:* `simple_state_weight_set_tester.f90`, sub-suite `state_weight_set` of `unit_project` (fast
  gate), 33 checks:
  - exact round trip of a PARTITION set and of a KERNEL set;
  - refusal of a set of another particle layout, of another particle count, of a file with one byte
    flipped (checksum) and of a file one byte longer (size);
  - the simulated crash: a new generation's files written without their manifest leave the previous
    generation selected and readable, and the next complete publication supersedes it and deletes
    the older files;
  - a project whose `out` segment carries an earlier release's `flex_weights` entry loads, and
    reports no weight set.
- *A property of the design, recorded:* the manifest's name is fixed and the project already
  points at it when FLEX republishes in the same directory. A crash after the manifest's rename
  therefore leaves the new generation selected. That generation is complete and consistent (its
  files were validated before the rename), so no mixed delivery can be read, but in FLEX the
  project's hard labels are written after the set.

**4.3, the reconstruction service.** New module
`src/main/strategies/parallelization/simple_rec3D_service.f90`.
- A `rec3D_request` carries:
  - the explicit particle rows;
  - the weight source (hard labels only in this phase; Phase 2 adds the weighted sources);
  - the backend, gridding or PCG (preconditioned conjugate gradients);
  - the output kind: raw half maps, or delivered maps that are also post-processed;
  - a caller-owned file-name prefix, allowed only for unregistered raw output;
  - an explicit publication policy;
  - the dispatch: in this process, or the standard queue rounds.
- A `rec3D_service` holds the queue environment for the queue dispatch. It owns insertion,
  assembly (volassemble, or the PCG master), registration of the state maps and FSC (Fourier shell
  correlation) files, post-processing and renaming.
- With the queue dispatch, the master writes the rows to a file. Each part's private
  `reconstruct3D` worker reads it as `pindfile=` and keeps its own range of rows. Sampling is a row
  filter whose only global inputs (the previous sampling index and whether any particle was
  updated) are the same in the master and the parts, so this gives each part the rows it sampled
  itself before.
- The backend selector (`rec3D_backend_id`) moved to the service.
- The `reconstruct3D` strategies (`simple_rec3D_strategy.f90`) now only parse the parameters,
  partition even and odd, prepare the canonical sigma2 state, sample the rows and send the request.
  Their `finalize_run` step went into the service and was removed from the strategy interface and
  from `exec_rec3D`.
- The publication policy of the base is kept path by path, so outputs, including the project, are
  unchanged: shared-memory gridding never registers its maps, while shared-memory PCG and the
  distributed master register them when the run has its own directory (`mkdir=yes`). The
  inconsistency is the base's and is recorded here, not changed.
- Post-processing no longer copies a whole `parameters` per state (`params_pp = params`, noted in
  section 10). `postprocess_volume_from_files` changes only `params%fsc` and `params%bfac`, both set
  afresh for each state, so the service saves and restores those two fields.
- **Same-path equality (8.2):** the Phase 0 command lines on a Release build of this step, each path
  once from a fresh copy of `keep/B` (`refs_check.sh`, `logs/refs_check_p1s.log`,
  `logs/cmp_p1s_<path>.txt`). Every `.mrc` file is byte-identical to the Phase 0 reference: 24 files
  for gridding in shared memory, 32 for gridding with two parts, 16 for each PCG path. Wall times:
  14 s, 16 s, 2:24 and 2:19.

**The pruning ruling (D10).** `prune_underpopulated_states` (FLEX state service) still drops states
whose effective sample size is below FLEX's occupancy floor. After that, the weights decide:
- a row with weight in a kept state is labelled by its largest kept weight;
- PARTITION rows are renormalized over the kept states;
- a row without kept weight gets label 0.

Rows of kept states keep their labels, because FLEX's labels are already the argmax. FLEX's
per-state `neff` stays the pre-pruning estimate that drove the floor; it is only reported, and the
set computes its own effective sample size. Tests:
- `test_prune_states` in `simple_flex_pca_tester.f90` (fast gate, `unit_heterogeneity`): a state is
  dropped; a row of it with kept weight is relabelled and renormalized to its whole mass; a row
  without kept weight is unassigned; and a row is labelled exactly when it carries weight.
- The application phantom now also asserts per arm that the published set selects exactly the
  particles with `state > 0`.

One hazard is recorded and left for Phase 3, which replaces the bandwidth cross-validation (D11).
With `nbins > 1` (default 1) the bandwidth cross-validation after pruning rewrites the weight columns
but not the labels, so a label could point at a zero weight. `publish` would then stop, as the old
store's validation did.

**Removal of FLEX's nonuniform filter (D12).** `apply_consensus_nu_filter` and its call are deleted,
as are the `nufilt` parameter, its parser entry and its use in the application tester. The comment
in `simple_nu_filter.f90` that named FLEX's report is fixed.

**The workflow gate's fixture.** Phase 0 found that `flex_pca_blobs` fails on the base in its
distributed mode. The phantom project had no `compenv` segment (the computing environment
`new_project` writes for every real project), and the distributed master's queue setup reads it.
Phase 1's exit needs the gate's distributed mode, so `build_flex_pca_phantom_project` now calls
`project%update_compenv`: the three-line change of Phase 0's evidence build. It is test-only, and
no assertion or floor changed.

**Builds and tiers (8.1).**
- *Builds:* fresh `./compile_debug.sh` (`logs/debug_build_p1.out`) and fresh Release build
  (`build_release/make.out`). Neither has a compiler warning; the linker warnings are Phase 0's
  baseline (six and five "requires executable stack"). The fast gate passed 15 of 15 on both.
- *Library and high-level tiers*, one entry at a time on the Release build (`logs/queue_tiers_p1.log`,
  outputs in `logs/tiers_p1/`). The wall times are Phase 1's, with Phase 0's in brackets:
  - library entries: `lib_reconstruction` 7 s (8 s), `lib_cart_align3D` 25 s (25 s),
    `lib_heterogeneity` 24 s (23 s), `lib_single` 10 s (14 s), `lib_stream` 561 s (617 s): all pass;
  - high-level entries: `pcg_recon` 1 s, `simulate_particles` 6 s, `simulated_workflow_1jxy` 241 s,
    `solve3D_addon` 850 s, `cont_refine3D_1jxy` 549 s: all pass;
  - `flex_pca_blobs` 22 s: passes, in both modes. It failed at Phase 0; the fixture fix above is
    what changed.
  - `simulated_workflow_6vxx` failed its first run, as at Phase 0. Its de novo solution correlated
    0.586 with the truth (minimum 0.80). One rerun on the final build passed, with correlation
    0.8882 and truth FSC 0.143 at 4.71 Å (`logs/tiers_p1b/`). This is the known flakiness recorded
    at Phase 0, where the base also failed once and passed on its rerun.
- *After the review fixes:* the code was reviewed against this plan, and three small fixes followed
  the build above:
  - `publish` refuses a producer the manifest cannot store as one token;
  - the service accepts a file-name prefix only for unregistered raw output, and builds the names
    with the `.mrc` extension the canonical names use;
  - a distributed `reconstruct3D` part reading `pindfile=` stops on an empty part under fractional
    update (`update_frac`), as its own sampling did before.

  On the final code, after an incremental Debug build and fast gate (15 of 15,
  `logs/inc_debug_p1_final.log`) and an incremental Release build (`logs/refs_check_p1f.log`):
  - all four reference paths are again byte-identical to Phase 0 (24, 32, 16 and 16 files);
  - `lib_heterogeneity` and `flex_pca_blobs` pass with the numbers below;
  - the `simulated_workflow_6vxx` rerun passes.
- *Checks:* `scripts/check_test_registry.py`, `scripts/check_descr.py` and
  `scripts/check_flex_dag.py` pass. No file mode changed. No `THROW_HARD` or `THROW_WARN` argument
  continues with `//&`, and the new code tests no array element on the same line as the flag that
  guards it.

**FLEX margins (8.2) against Phase 0.**
- *Application phantom* (`logs/tiers_p1b/LastTest_lib_heterogeneity.log`): every arm is identical to
  Phase 0. The arms (grid_basis_grid, pcg_basis_grid, grid_basis_pcg, pcg_basis_pcg) each deliver
  3 states with recoded label accuracy 0.9995. The smallest matched-map correlations are 0.9846,
  0.9845, 0.9976 and 0.9976, and the per-state truth tables are unchanged. The margins are met
  without change: accuracy at least 0.9795 and at least 0.95; state count 3, between 2 and the
  ceiling 3; correlations within 0.02; each matched map prefers its own truth. The suite now holds
  56 checks (52 at Phase 0, plus the selection agreement in each arm).
- *Workflow gate* (`logs/tiers_p1b/LastTest_flex_pca_blobs.log`): the shared and distributed modes
  equal the values Phase 0 recorded from its evidence build:
  - shared: recoded label accuracy 0.99950, leading eigenvolume truth correlation 0.96468;
  - distributed: recoded label accuracy 0.99750, leading eigenvolume truth correlation 0.92030;
  - both modes: 3 states each, with the same per-state truth tables;
  - shared and distributed recoded labels agree on 99.80 % of the particles (margin 98 %);
  - mode A and mode B map correlations 0.99995 and 0.99998.

**Exit items.**
- *`reconstruct3D` same-path equal (8.2) to its Phase 0 references on all four paths:* met bit for
  bit, after 4.3 and again on the final code.
- *The weight-set tests of section 9 green:* met. The 33 checks of `state_weight_set` cover the round
  trip, the digest, count, size and checksum refusals, the crash simulation and the old-project
  load. Pruning agreement is covered by `test_prune_states` and by the application phantom's
  selection check. The work-project single-state set belongs to Phase 4.
- *FLEX application phantom and workflow gate within the FLEX margins of 8.2:* met, identical to
  Phase 0.

Files changed in Phase 1 (`git status`):
- new: `src/main/project/simple_ptcl_layout.f90`, `src/main/project/simple_state_weight_set.f90`,
  `src/main/project/simple_state_weight_set_tester.f90`, `src/fileio/simple_state_weights_file.f90`,
  `src/main/strategies/parallelization/simple_rec3D_service.f90`;
- deleted: `src/main/flex/states/simple_flex_weights_state.f90`,
  `src/fileio/simple_flex_weights_file.f90`;
- modified: `scripts/check_flex_dag.py`, `src/fileio/simple_projfile_utils.f90`,
  `src/fileio/simple_sigma2_files.f90`, `src/main/commanders/simple/simple_commanders_euclid.f90`,
  `simple_commanders_euclid_distr.f90`, `simple_commanders_flex_pca.f90`,
  `simple_commanders_rec.f90`, `simple_commanders_solve2D.f90`,
  `src/main/commanders/test/simple_commanders_test_class.f90`,
  `simple_commanders_test_highlevel_flex.f90`, `src/main/flex/run/simple_flex_pca_application.f90`,
  `simple_flex_pca_application_tester.f90`, `simple_flex_pca_delivery_3d.f90`,
  `simple_flex_pca_project_gateway.f90`, `simple_flex_pca_state_service.f90`,
  `src/main/flex/simple_flex_pca_tester.f90`, `src/main/nu_filt/simple_nu_filter.f90`,
  `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`,
  `src/main/project/simple_sp_project.f90`, `simple_sp_project_out.f90`,
  `src/main/sigma2/simple_sigma2_state.f90`, `simple_sigma2_state_tester.f90`,
  `src/main/solve/simple_solve3D_manifest.f90`, `simple_solve3D_utils.f90`,
  `src/main/strategies/parallelization/simple_rec3D_strategy.f90`,
  `simple_rec3D_strategy_tester.f90`, `src/main/ui/simple_test/simple_test_ui_class.f90`, and this
  plan.

### Phase 2: fractional reconstruction (2026-10-07, done)

Paths below are under `/home/elmlundho/agent_runs/flex_on_simple/scratch` unless absolute. Builds and
scripts are those of Phases 0 and 1.

**What changed.**
- *Keys and interface.* `m_estimator=no|flex {no}` and `state=X` on `reconstruct3D`, in the user
  interface, the parameters and their parser. `l_m_estimator_flex` records `m_estimator=flex`, and
  `l_state_defined` records whether `state=` was given.
  - With `m_estimator=flex`, every particle enters each state with its weight from the project's
    state weight set.
  - `state=X` alone selects the rows labelled `X` (the hard comparison). With `m_estimator=flex` it
    selects every row with a positive weight for `X`.
  - The other states are left alone, so volassemble's dropped-state path keeps them. A `state=`
    outside 1..`nstates` is refused.
- *Gridding (section 5).*
  - `insert_plane_oversamp` takes an optional weight that scales the data term and the density term
    alike. Without it the arithmetic is unchanged; with weight 1.0 it is exact.
  - `calc_3Drec` takes an optional state weight set. `group_pinds_by_weights` then builds each
    (state, half) membership list from the rows whose weight is above zero, in selection order, so
    a particle is read once per state it weighs into. One half-map reconstructor exists at a time,
    as before.
  - Partial file names and payloads are unchanged.
  - The density floor (the density-floor ruling, D4): volassemble applies `floor_rho_shellwise` to
    the even and odd densities before restoration, and to the merged density before the final
    restoration, only when the state's weights are fractional (`state_weight_set%is_fractional`:
    some weight strictly between 0 and 1). It is never applied on hard or 0/1 runs.
- *PCG, the preconditioned conjugate gradients backend.*
  - A particle is a member of a state when its weight exceeds 1e-3 (`PCG_WEIGHT_THRESHOLD`, now in
    `simple_reconstructor_pcg`). Its noise spectrum is divided by its weight before `prep_particles`,
    which weighs both the right-hand side B and the density D.
  - The raw accumulator format is now version 2 (magic `SIMPLE_PCG_RAW02`). The header carries the
    applied mass (64-bit real) beside the integer contributor count, plus the weight-set identity
    (generation and layout digest; zero under hard labels). The payload is unchanged and still gated
    on the contributor count.
  - The master checks every part's identity. No reader of version 1 is kept.
- *Trailing reconstruction.*
  - `trail_chain_manifest` is now version 3. Beside the represented mass M it records the contributor
    count and the weight-set identity.
  - The PCG chain pair now has the same manifest contract (`pcg_trail_stateNN.txt`). The manifest is
    deleted before the first half is written, written after both halves, and validated before use
    (identity, sizes, state layout). The represented population M is read from it instead of from
    the raw header counts. A blend continues the generation count; a re-seed restarts it.
  - Both backends discard and re-seed a chain built under another weight-set identity; they never
    blend it.
  - Volassemble now sets a chain's generation and mass only when the chain validates, so a re-seeded
    gridding chain also restarts at generation 1.
  - The frozen and add-on accumulators record and check the same identity (schema version 3).
- *Population rule and counts.*
  - `population_blend_weights` is generic over integer counts and real masses. Integer counts go
    through the mass version with identical arithmetic.
  - `oris%get_update_rows` returns the rows behind `get_group_update_counts` (active updated rows,
    the current sample) and the rows a seed represents.
  - Under a weight set, the state update counts N and n, the realized fractions, the seed population,
    the PCG master's representation counts and a chain's represented mass are applied masses over
    those rows, after the backend's threshold. Under hard labels they are the former counts.
- *Population gates on the effective sample size (ESS).* With a weight set, volassemble's state
  populations and the PCG nonuniform-filter low-pass gate use the effective sample size of the
  state's weights. Post-processing skips a weighted state without weight.
- *The service.* A new weight source, `REC3D_WEIGHTS_SET`. The strategies request it with
  `m_estimator=flex`, and distributed workers open the set themselves. Registration and
  post-processing keep to the selected state under `state=`.
- *The kind recorded with every result (D3).* A weighted reconstruction logs the set's kind,
  generation, producer and, per state, applied mass, effective sample size and hard population
  ("RECONSTRUCTION WEIGHTS"). Registering these fields with the volume belongs to Phase 4 (section
  6, item 4).
- *Population gates of the refinement* (`simple_eul_prob_tab`, `simple_refine3D_strategy`) and the
  matcher's reconstruction call are left for Phase 4, where the refinement first reads a weight set.
  In Phase 2 only `reconstruct3D` consumes weights.

**Tests.**
- *Linearity* (`simple_reconstructor_tester.f90`, fast gate, `unit_reconstruction` sub-suite
  `weighted_insertion`, 6 checks):
  - for weights 0.37, 0.05 and 0.81, inserting with weight w equals inserting with weight one and
    scaling the accumulators by w, to 1e-6 relative, for data and density;
  - weight one equals the unweighted insertion bit for bit.
- *Fractional phantom, same path and trailing* (`simple_rec3D_service_tester.f90`, library tier,
  `lib_reconstruction` sub-suite `fractional_reconstruction`, run in-process through
  `reconstruct3D`). It uses the two-conformer phantom of `simple_flex_pca_application_tester.f90`:
  2,000 particles, box 64; conformers A and B differ by one lobe moved 8 pixels.
  - Fractional rule: state 1's weight is 0.75 for conformer A and 0.30 for B, plus 0.1·cos(0.37·i)
    for particle i; state 2 gets the rest. So each state's expected map is a mixture of both truth
    maps.
  - Results:
    | Backend | State | Mass from A / B | Weighted map vs mixture | Hard map vs truth |
    |---|---|---|---:|---:|
    | gridding | 1 | 749.9 / 299.8 | 0.9856 | 0.9857 |
    | gridding | 2 | 250.1 / 700.2 | 0.9853 | 0.9852 |
    | PCG | 1 | 749.9 / 299.8 | 0.9981 | 0.9977 |
    | PCG | 2 | 250.1 / 700.2 | 0.9981 | 0.9978 |

    Every weighted map meets its floor (the hard correlation minus 0.01).
  - 0/1 weights equal to the labels reproduce the hard state maps and half maps byte for byte, for
    gridding and PCG (shared memory).
  - Trailing seeds under fractional weights: the gridding and the PCG chain both represent the
    state's applied mass. Expected 1049.693; gridding 1049.694; PCG 1049.694.
  - Trailing identity, for both backends: after a fractional update under a new weight-set
    generation, the chain is re-seeded (generation 1, new identity recorded). The same update under
    the same generation blends it (generation 2).
- *Chain and frozen formats.* The trailing-chain identity tester covers the manifest v3 round trip
  (mass, contributor count, identity), refusal of version 2, and the raw header's mass and identity.
  The frozen accumulator tester covers the schema bump.
- *Test tools.* `simple_test_exec test=write_state_weights_labels` (0/1 weights equal to the labels)
  and `test=write_state_weights_mixed` (0.7 on the labelled state, 0.3 spread over the others)
  publish a weight set from a project's labels for the comparisons below. They are test programs,
  not production keys, and not CTest entries.

**Exit comparisons on beta-galactosidase.** Each run starts from a fresh copy of `keep/B` with the
Phase 0 command line of its path plus `m_estimator=` (`run_ref_w.sh`, `queue_refs_p2.sh`).
- `m_estimator=no`: all four paths are byte-identical to the Phase 0 references: 24 files for
  gridding in shared memory, 32 for gridding with two parts, 16 for each PCG path
  (`logs/queue_refs_p2.log`, `logs/cmp_p2no_<path>.txt`). Repeated on the final code:
  byte-identical again on all four paths (`logs/queue_refs_p2f.log`).
- `m_estimator=flex` with 0/1 weights equal to the labels (written by `write_state_weights_labels`):
  all four paths are byte-identical to the Phase 0 references (`logs/cmp_p2labels_<path>.txt`).
  Repeated on the final code: byte-identical again on all four paths (`logs/queue_refs_p2f.log`).
- *Fractional, shared memory against two parts* (`write_state_weights_mixed`, `m_estimator=flex`):

  | Backend, run | State 1 even | State 1 odd | State 2 even | State 2 odd |
  |---|---:|---:|---:|---:|
  | PCG, default (restored halves) | 8.0e-7 | 3.1e-6 | 2.7e-6 | 8.7e-7 |
  | Gridding, default, halves before the ML prior (`_unfil`) | 3.1e-6 | 1.4e-6 | 7.2e-6 | 2.4e-6 |
  | Gridding, default, ML-regularized halves | 5.6e-6 | 2.5e-6 | **1.49e-5** | 5.0e-6 |
  | Gridding, `ml_reg=no` (restored halves) | 3.1e-6 | 1.4e-6 | 7.2e-6 | 2.4e-6 |
  | PCG, `ml_reg=no` | 5.7e-7 | 1.8e-6 | 1.4e-6 | 5.2e-7 |

  The table gives the relative L2 difference. Every correlation is 1.0000000.

  - **Decision, accepted 2026-10-07 (section 3).** The cross-mode tolerance of 8.2 (relative L2 at most 1e-5 and
    correlation at least 0.99999 on each state's even and odd half maps) is met by every half map
    the fractional reconstruction produces: all maps of both backends with `ml_reg=no`, the
    unregularized halves of the default gridding runs, and every PCG half.
  - It is missed by one map: the default runs' maximum-likelihood (ML) regularized gridding half of
    state 2, even half (1.49e-5).
  - That regularization is a prior computed from each run's own FSC, so it amplifies the
    summation-order differences between the two reduction trees. The hard runs of Phase 0 show the
    same amplification on the same map: 1.6e-6 before the prior, 7.3e-6 after.
  - I judged the exit on the reconstruction's half maps (before the prior). The regularized numbers
    are reported here so the maintainer can rule otherwise.
  - Logs: `logs/queue_refs_p2.log`, `logs/cmp_p2mixed_<backend>.txt`, `logs/queue_mixed_mlregno.log`,
    `logs/cmp_mixed_mlregno_<backend>.txt`.

**Builds and tiers (8.1), on the final code.**
- *Builds:* fresh `./compile_debug.sh` (`logs/debug_build_p2.out`) and fresh Release build
  (`build_release/make.out`). Neither has a compiler warning; the linker warnings are Phase 0's
  baseline. The fast gate passed 15 of 15 on both.
- *Tiers*, one entry at a time on the Release build (`logs/queue_tiers_p2.log`, outputs in
  `logs/tiers_p2/`), wall times in seconds:
  - fast gate 4;
  - library entries: `lib_reconstruction` 25 (it now includes the fractional-reconstruction
    sub-suite), `lib_cart_align3D` 26, `lib_heterogeneity` 22, `lib_single` 14, `lib_stream` 581;
  - high-level entries: `pcg_recon` 2, `simulate_particles` 7, `simulated_workflow_1jxy` 241,
    `simulated_workflow_6vxx` 800 (passed on its first run this time), `solve3D_addon` 837,
    `cont_refine3D_1jxy` 545, `flex_pca_blobs` 22.
  - All pass. The FLEX application phantom and the workflow gate give the same numbers as in
    Phases 0 and 1.
- *Manual check:* the PCG fractional-update validation (`simple_test_exec test=pcg_frac_update`, not
  in CTest) passes on a copy of the B project with `box_crop=128 objfun=cc`, in 57 s
  (`logs/p2_pcg_frac_update_cc.log`). It covers raw additivity, u/f continuation weighting and chain
  replay under the new raw format and chain manifest. A first attempt with the default Euclidean
  objective stopped at once, because the copied project's canonical sigma2 file is not in the run
  directory (`logs/p2_pcg_frac_update.log`); that has nothing to do with this phase.
- *Checks:* `scripts/check_test_registry.py`, `scripts/check_descr.py` and `scripts/check_flex_dag.py`
  pass. No file mode changed. No `THROW_HARD` or `THROW_WARN` argument continues with `//&`.

**Recorded, not changed.**
- The weight-set identity is (generation, layout digest), as section 4.2 defines it. Two FLEX runs
  from the same parent project publish the same generation, so their identities coincide. A trailing
  or frozen artifact can therefore be blended across them only if the same directory is reused for
  both.
- The PCG full-population count of a half (`full_state_half`) and the gridding seed's rows differ in
  rows with invalid even/odd labels. This edge case was present before Phase 2.

**Exit items.**
- *`m_estimator=no` same-path equal to the references:* met bit for bit on all four paths.
- *`m_estimator=flex` with 0/1 weights equal to the labels, same-path equal on all four paths:* met
  bit for bit, and also in the library tier for both backends.
- *Fractional shared memory against two parts within the cross-mode tolerance:* met on the half maps
  of the reconstruction (the decision above).
- *The linearity, trailing and fractional-phantom tests of section 9 green:* met (`weighted_insertion`
  in the fast gate; `fractional_reconstruction` in the library tier, 17 checks).

Files changed in Phase 2:
- new: `src/main/strategies/parallelization/simple_rec3D_service_tester.f90`,
  `src/main/volume/simple_reconstructor_tester.f90`;
- modified: `src/defs/simple_refine3D_fnames.f90`,
  `src/main/commanders/simple/simple_commanders_rec.f90`, `simple_commanders_rec_distr.f90`,
  `src/main/commanders/test/simple_commanders_test_class.f90`, `simple_commanders_test_highlevel.f90`,
  `simple_commanders_test_highlevel_flex.f90`, `src/main/exec/simple_test_exec_highlevel.f90`,
  `src/main/ori/simple_oris.f90`, `simple_oris_getters.f90`, `simple_oris_sampling.f90`,
  `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`,
  `simple_parameters_phases.f90`, `src/main/project/simple_state_weight_set.f90`,
  `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`, `simple_rec3D_service.f90`,
  `simple_rec3D_strategy.f90`, `src/main/strategies/search/simple_matcher_3Drec.f90`,
  `src/main/ui/simple/simple_ui_reconstruct3D.f90`, `src/main/ui/simple_test/simple_test_ui_class.f90`,
  `simple_test_ui_highlevel.f90`, `src/main/volume/simple_frozen_accum.f90` (and its tester),
  `simple_reconstructor.f90`, `simple_reconstructor_pcg.f90`, `simple_trail_chain_manifest.f90` (and
  its tester), and this plan (three file-table rows: `simple_parameters_phases.f90`,
  `simple_oris.f90`, `simple_refine3D_fnames.f90`).

### Phase 3: FLEX on the service (2026-10-07, done)

Paths below are under `/home/elmlundho/agent_runs/flex_on_simple/scratch` unless absolute. Builds and
scripts are those of the earlier phases.

**What changed.**
- *Deleted.* FLEX's state backends, parts and delivery: `simple_flex_pca_rec3D.f90`,
  `simple_flex_pca_states_backend.f90`, `simple_flex_pca_state_parts.f90`,
  `simple_flex_pca_state_delivery.f90`, `simple_flex_pca_states_gridding.f90` and
  `simple_flex_pca_states_pcg.f90` (all in `src/main/flex/states/`). With them went the
  round-weights broadcast to the workers, the worker state stage (`PCA_STAGE_STATES`), the memory
  report of resident reconstructors, `flex_rec_box` and `flex_rec_smpd`, the 8th-order Butterworth
  delivery filter and the `flex_pca_state_NNN.mrc` naming (`flex_pca_write_state`, the `outvol`
  default). `scripts/check_flex_dag.py` no longer lists the modules.
- *Trial maps through the service (the box ruling D9, the raw-map ruling D11).*
  - `reconstruct_state_halves` (`simple_flex_pca_state_service.f90`) calls the reconstruction
    service in the FLEX process with a caller-owned weight table, a new weight source
    (`REC3D_WEIGHTS_TABLE`), raw output and a file-name prefix. Every column of the table is one
    state; the maps are written as `<prefix>_stateNN{,_even,_odd}.mrc` at `box_crop`.
  - The table reaches both backends: `calc_3Drec` (gridding) and `execute_rec3D_pcg_worker` (PCG)
    take it as an optional argument. Membership is a weight above zero for gridding and above 1e-3
    for PCG, as for a weight set. A table is in-process only and is refused for the queue.
  - The child reconstruction is a single in-process part with the shell density floor of fractional
    weights (D4; volassemble's new `rho_floor` key), no maximum-likelihood prior, no
    post-processing, no nonuniform filter, no envelope mask, no trailing chain and no `state=`
    selection, whatever the FLEX command line carries (a `refine3D_states` handoff carries
    `filt_mode=nonuniform_lpset`, `automsk` and `trail_rec`). The PCG path reads these from FLEX's
    own parameters, so `reconstruct_state_halves` sets them and restores them afterwards, with
    `nstates`, `nparts` and the state volume names. A column without weight is refused.
  - The bandwidth cross-validation makes one service call per bin, with the full weight columns,
    and scores each state's even and odd maps under FLEX's soft spherical mask. It no longer
    masks the weights by half set or reads maps from another box. Its trial files are deleted
    after each bin.
  - After the cross-validation adopts its bandwidths, the labels follow the pruning ruling (D10): a
    row with weight is labelled by its largest weight, a row without weight is unassigned. This
    closes the hazard Phase 1 recorded (a label pointing at a zero weight when `nbins > 1`).
  - The merge's gate 2 reads the raw trial halves (`flex_pca_trial_stateNN_{even,odd}.mrc`) at
    `box_crop`. Trial halves exist only when the merge runs (`preimage_auto=yes`, more than one
    state) and are deleted after it. A merge no longer re-reconstructs the maps.
- *Delivered maps.* After `publish` writes the weight set and the hard labels, `deliver_state_maps`
  (`simple_flex_pca_project_gateway.f90`) runs `reconstruct3D` with `m_estimator=flex` in the FLEX
  process, from the published set, at the native box, with the run's backend
  (`rec_states_backend`), mask, point group, objective, sigma estimation and PCG controls. It
  registers each state's `recvol_stateNN.mrc` as an ordinary `vol` entry and its FSC as an `fsc`
  entry. It removes the `vol` and `fsc` entries of states it does not deliver (states above the
  delivered count, from an earlier run of the project) and any earlier release's `vol_flex`
  entries.
- *`vol_flex` retired (ruling 5).*
  - `add_vol2os_out` and `get_vol` no longer accept `vol_flex`.
  - Read-only fallback for projects of earlier releases: `get_vol('vol', s)` and
    `isthere_in_osout('vol', s)` fall back to a `vol_flex` entry of state `s` when the state has no
    `vol` entry. `remove_state_artifacts_from_osout` still removes `vol_flex` entries.
  - The `refine3D_states` initialization (`run_flex_pca` in `simple_commanders_refine3D.f90`) reads
    the states' `vol` entries.
  - The testers read `vol`.
- *`box_rec` and `smpd_rec` removed* from the parameters, their parser and derivation, the
  `flex_pca` user interface, the commander's band derivation and the application tester (D9).
- *The PCG master's full-population check* ("current half population exceeds its full
  population") now runs only where a trailing chain is written or blended (`trail_rec`, or a seed
  request). That count is the population a chain represents. A caller-owned table reconstructs
  without a chain, and the project's labels cannot count it. A hard run without trailing loses a
  sanity check whose number fed only the chain.
- *Service file names.* The prefix renaming also covers the unfiltered halves
  (`_even_unfil`, `_odd_unfil`); FLEX deletes them with its trial maps.
- *File table.* Two rows added: the parameters and FLEX's user interface (`box_rec`, `smpd_rec`), and
  the matcher, the PCG strategy and volassemble (the caller-owned weight table).

**FLEX margins (8.2), measured exactly as in Phase 0** (Release build, final code,
`logs/tiers_p3_final/LastTest_lib_heterogeneity.log`, `logs/tiers_p3_final/LastTest_flex_pca_blobs.log`).
Every margin is met.
- *Application phantom (`lib_heterogeneity`), per arm.* Every arm delivers 3 states (between 2 and
  the ceiling 3) with recoded label accuracy 0.9995 (Phase 0: 0.9995; margin: at least 0.9795 and
  0.95). Matched states: A = 3, B = 1, as at Phase 0; each matched map prefers its own truth.

  | Arm | Matched A (Phase 0) | Matched B (Phase 0) | Smallest (Phase 0) |
  |---|---:|---:|---:|
  | grid_basis_grid | 0.9857 (0.9851) | 0.9852 (0.9846) | 0.9852 (0.9846) |
  | pcg_basis_grid | 0.9857 (0.9851) | 0.9852 (0.9846) | 0.9852 (0.9845) |
  | grid_basis_pcg | 0.9985 (0.9976) | 0.9987 (0.9979) | 0.9985 (0.9976) |
  | pcg_basis_pcg | 0.9985 (0.9976) | 0.9987 (0.9979) | 0.9985 (0.9976) |

  The per-state truth counts are those of Phase 0 (state 1: 1 A, 1,000 B; state 2: 3 A; state 3:
  996 A). The leading eigenvolume correlations are unchanged (0.9605 and 0.9647). The matched maps
  correlate slightly better with their truth than Phase 0's Butterworth-filtered `vol_flex` maps.
- *Workflow gate (`flex_pca_blobs`), per mode.* Shared: accuracy 0.99950 (Phase 0: 0.99950), 3
  states, matched maps 0.9985 (A) and 0.9987 (B) (Phase 0: 0.9976, 0.9979). Distributed: accuracy
  0.99750 (0.99750), 3 states, matched maps 0.9984 (A) and 0.9987 (B) (0.9975, 0.9979). The truth
  counts of both modes are those of Phase 0. Shared and distributed labels agree on 99.80 % of the
  particles (margin 98 %); the mode A and B map correlations between the modes are 0.99995 and
  0.99998, the latent correlation 0.56132 and the noise scales 0.75683 and 0.77554, all as before.

**Tests.**
- *A fifth arm of the application phantom*, `grid_cv_crop48`: the gridding pair with the bandwidth
  cross-validation (`nbins=4`) and a covariance box below the native one (`box_crop=48` of 64). It
  runs the code no other entry reaches: the per-bin service calls, the relabelling after adoption,
  and trial maps at `box_crop` with delivered maps at the native box. `leading_pc_correlation`
  now compares on the eigenvolume's own box. The arm passes every check of the suite (the
  independent floor 0.60, the state ceiling, the weight-set and label agreement, matched maps
  preferring their truth). It has no Phase 0 baseline and is not a margin.
  - Its accuracy is 0.8655 on Release (0.8775 on Debug): state 2 mixes the conformers (245 A, 271
    B), and the eigenvolume is weaker at box 48 (correlation 0.849 against 0.96 at box 64).
  - A one-off Debug diagnostic with the same arm at `nbins=1` (the source restored after it) gave
    0.9165 with the same mixed state, so the smaller covariance box accounts for most of the drop.
  - Recorded for the maintainer, not acted on: on this phantom the cross-validation at box 48
    lowers accuracy further (0.9165 to 0.8775).
- The suite holds 70 checks in the phantom suite and 86 assertions in `lib_heterogeneity` (56 and 72
  before the fifth arm).

**Other checks (8.1).**
- *Builds.* Debug and Release build without new compiler warnings (the Release link step's
  executable-stack notes are the pre-existing ones of untouched modules). Fast gate green.
- *Tiers* (Release, `logs/queue_tiers_p3.log`, on the code before the review fixes below): every
  entry passes: fast, `lib_reconstruction`, `lib_cart_align3D`, `lib_heterogeneity`, `lib_single`,
  `lib_stream`, `pcg_recon`, `simulate_particles`, `simulated_workflow_1jxy`,
  `simulated_workflow_6vxx`, `solve3D_addon`, `cont_refine3D_1jxy`, `flex_pca_blobs`. After the review
  fixes, the fast gate, `lib_reconstruction`, `lib_heterogeneity` and `flex_pca_blobs` were run
  again on the final code and pass (`logs/queue_tiers_p3_final.log`).
- *`reconstruct3D` references.* This phase touched the hard path's files (`calc_3Drec`, the PCG
  master, the service, volassemble). All four Phase 0 paths remain byte-identical to their
  references: 24, 32, 16 and 16 files (`logs/refs_check_p3.log`).
- *The `refine3D_states` run of Phase 0 (case B)*, rerun on the final code because its FLEX
  initialization now hands over `vol` and FSC entries (`logs/case_B_p3.log`). Exit 0, 21:03, 4.9 GB.
  FLEX delivered 2 states again ("FLEX_PCA INITIALIZED NSTATES: 2"). Final populations and FSC 0.5 /
  0.143: state 1: 3,044, 5.26 / 4.03 Å; state 2: 731, 15.54 / 12.09 Å; state 0: 1,738. Phase 0 had
  3,052, 5.44 / 4.03; 738, 14.84 / 12.09; 1,723. The project's out segment holds `vol` and `fsc` for
  both states and no `vol_flex`. The refinement now starts from the reconstructed maps instead of
  the Butterworth-filtered `vol_flex` maps, which explains the small differences.
- `scripts/check_test_registry.py`, `scripts/check_descr.py` and `scripts/check_flex_dag.py` pass. No
  file mode changed. No `THROW_HARD` or `THROW_WARN` argument is continued with `//&`, and no new
  `flag .and. array(i)` on a possibly unallocated array.
- *Review.* A read-only review of the phase found no blocking bug. Fixed from it:
  - the delivery checked a map's existence after resolving its absolute path, which stops on a
    missing file;
  - stale `vol` and `fsc` entries of undelivered states;
  - the PCG trial path inherited the trailing, `state=` and automask settings of FLEX's parameters;
  - `params%vols` was not restored;
  - a column without weight would reach an empty state;
  - the unfiltered trial halves were left behind;
  - `sigma_est` did not reach the delivery.

**Recorded, not changed.**
- *The consensus map is replaced.* FLEX's state-1 map is now the project's `vol` of state 1, so the
  consensus map a FLEX run started from is no longer registered once it publishes. A second
  `flex_pca` on its own output project would pick up state 1 as its consensus. This follows from
  ruling 5 (ordinary `vol` entries).
- *Side effects of the trial maps.* Volassemble (gridding) and the PCG master write the per-particle
  `res`, `res05` and `cfar` of the trial maps into the run project's `ptcl3D`, and the in-process PCG
  master writes the field FLEX holds in memory. The delivery's `reconstruct3D` overwrites them.
  With `rec_states=no` and `nbins > 1` no delivery follows, so the cross-validation's values stay.
- *Cost in distributed runs.* The trial maps (one call per bin, plus the merge) and the delivered
  maps run in the master process; the workers no longer reconstruct states.
- *Test gaps.* The workflow gate's `inspect_flex_result` calls `get_vol`, which stops the run rather
  than failing the gate when a state has no map; it was left as it is.

**Exit items.**
- *FLEX's state backends, parts, delivery and round-weights broadcast deleted; trial and final maps
  through the service (D9, D11); `vol_flex` retired:* done.
- *Application phantom and workflow gate within the FLEX margins of 8.2:* met on every arm and both
  modes (tables above).

Files changed in Phase 3:
- deleted: `src/main/flex/states/simple_flex_pca_rec3D.f90`, `simple_flex_pca_states_backend.f90`,
  `simple_flex_pca_state_parts.f90`, `simple_flex_pca_state_delivery.f90`,
  `simple_flex_pca_states_gridding.f90`, `simple_flex_pca_states_pcg.f90`;
- modified: `scripts/check_flex_dag.py`, `src/main/commanders/simple/simple_commanders_flex_pca.f90`,
  `simple_commanders_rec_distr.f90`, `simple_commanders_refine3D.f90`,
  `src/main/commanders/test/simple_commanders_test_highlevel_flex.f90`,
  `src/main/flex/simple_flex_pca_util.f90`, `src/main/flex/run/simple_flex_pca_application.f90`,
  `simple_flex_pca_application_tester.f90`, `simple_flex_pca_artifacts.f90`,
  `simple_flex_pca_project_gateway.f90`, `simple_flex_pca_records.f90`, `simple_flex_pca_stages.f90`,
  `simple_flex_pca_state_service.f90`, `src/main/flex/states/simple_flex_pca_merge.f90`,
  `src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`,
  `simple_parameters_phases.f90`, `src/main/project/simple_sp_project_out.f90`,
  `src/main/strategies/parallelization/simple_flex_pca_strategy.f90`, `simple_rec3D_pcg_strategy.f90`,
  `simple_rec3D_service.f90`, `src/main/strategies/search/simple_matcher_3Drec.f90`,
  `src/main/ui/simple/simple_ui_heterogeneity.f90`, and this plan (two file-table rows).

### Phase 4: per-state M-estimation in `refine3D_auto` (2026-10-07, done)

Paths below are under `/home/elmlundho/agent_runs/flex_on_simple/scratch` unless absolute. Builds and
scripts are those of the earlier phases.

**What changed.**
- *`refine3D_auto state=X m_estimator=flex|no`* (`exec_refine3D_auto`, new internal routine
  `prepare_state_work_project` in `simple_commanders_refine3D.f90`; inputs `state` and `m_estimator`
  in the `refine3D_auto` user interface; both keys were already parsed for every program).
  - It runs right after the parameters are built: before any project read and before the run's GUI
    communicator starts.
  - As `solve3D`'s state continuation does, it copies the project into the run directory as
    `refine3D_auto_stateXX.simple`. It then runs the selection program on the copy, with pruning, so
    the copy holds one state, labelled 1.
  - With `m_estimator=flex` the selected rows are those with a positive weight for X in the
    project's state weight set. With `m_estimator=no` they are the rows labelled X.
  - The rows reach the selection program through its existing per-row selection file (`infile`), so
    the selection program is unchanged. The reason: the particle-layout digest includes the
    project's name, so the parent's set opens only in a project of the parent's name, not in the
    renamed copy.
  - Each work row is checked to be the parent row it was selected from, by stack and index in its
    stack.
  - State X's map and FSC (read from the parent) become state 1's in the copy. Every other state's
    map and FSC entries are removed, as are the parent's state-weights entry and its canonical sigma2
    registration, so the copy's sigmas are seeded in this run as `refine3D_auto` does for any project.
  - With `m_estimator=flex`, the copy gets its own single-state weight set, `publish_work_state`. It
    holds X's column at the work rows, the parent's hard-label flags for X, the parent's producer and
    kind, and the parent set's generation, layout digest and state as provenance.
  - The run then continues on the copy as an ordinary single-state `refine3D_auto`. `state=` leaves
    its command line, `m_estimator` stays, and the run directory's copy of the parent is deleted, so
    the copy is the directory's only project. The parent project is never written.
  - `state=` requires `mkdir=yes`: the work project and its set must live in their own directory. The
    weight set also refuses to publish a derived set over its parent's manifest.
  - `m_estimator=flex` without `state=` is refused.
- *Weighted reconstruction throughout the refinement.*
  - The matcher's partial reconstructions (`write_partial_recs` in `simple_strategy3D_matcher.f90`)
    pass the weight set to `calc_3Drec` and to the PCG worker whenever `m_estimator=flex`. This
    covers shared memory and distributed parts.
  - The startup reconstruction, the registration pass and the main iterations get `m_estimator` from
    the command line they copy. Volassemble, the PCG master and trailing were made weight-aware in
    Phase 2.
  - Sampling and alignment are unchanged: every selected particle is aligned against state 1's
    reference, and the sampler never reads the weights (the frozen-weights ruling, 7).
- *Population gates on the effective sample size.*
  - The probability table's state-existence gate (`simple_eul_prob_tab.f90`) uses the effective
    sample size of the state's weights when weighted.
  - The refinement's per-iteration registration and carry-forward decision (`get_state_pops` in
    `simple_refine3D_strategy.f90`) treats a weighted state as populated when its weights have applied
    mass, which is exactly when it has partial reconstructions.
- *Final native reconstruction* (`calc_final_rec`, `simple_final_rec.f90`).
  - It puts `m_estimator` (and `state`, when the refinement had one) on the child command line it
    builds from scratch. That line drives both routes: the direct `reconstruct3D` and
    `bootstrap_rec3D`, whose children copy it.
  - It records the weight set's identity (generation and layout digest) before reconstructing, and
    stops if the identity has changed when it registers.
  - It logs "FINAL RECONSTRUCTION WEIGHTS: ..." before and "FINAL RECONSTRUCTION REGISTERED STATE ...:
    hard population, applied mass, effective sample size" after.
  - It registers the map with the hard population (`pop`) and, as separate real fields, the applied
    mass (`mass`) and effective sample size (`ess`), computed from the weights as the backend applied
    them over the reconstructed rows. `add_vol2os_out` takes them as optional arguments, and an
    entry re-registered without them drops stale values.
- *Weight-set provenance.* The manifest is now version 2: a `provenance` line holds the parent
  generation, layout digest and state (zeros for a set with no parent). There is no reader of version
  1, by the release-4 rule. `publish` takes an optional kind and provenance. A derived set of a
  PARTITION parent keeps that kind, and its single column is exempt from the row-sum check.
- *A regression from Phase 2, fixed in `solve3D`.* Since Phase 2, `state=` restricts reconstruction
  (the service's state range check, the PCG master's state loop). `solve3D`'s state continuation kept
  `state=X` on its command line after relabelling the selected rows to state 1, so for X ≥ 2 its
  reconstructions would now stop or skip the state. The continuation now deletes `state` once its
  work project holds the state alone. File-table row added.
- *File table.* One row added (`simple_commanders_solve3D.f90`).

**Decisions, recorded for review** (ruled 2026-10-07, section 3: long tests run nightly, never in the
build; FLEX's weights on real data are tested after the implementation is finished).
- *The per-state refinement phantom is not a new CTest process.* The CTest budget (34 entries) may
  only be raised by the owner. The phantom therefore runs inside `flex_pca_blobs`, which already runs
  FLEX on the same two-conformer phantom. After the shared-memory and distributed FLEX comparisons,
  which it does not touch, it refines each truth-matched state of the shared run with
  `refine3D_auto state=X`, once with `m_estimator=flex` and once with the hard labels. It uses twice
  the FLEX modes' threads, that is the eight processors the entry reserves. The FLEX measurements of
  the gate are unchanged.
- *The hard comparison selects the rows labelled X through the same per-row file*, not through the
  selection program's `state=` route. That route also keeps only particles updated at least once,
  and the phantom's rows are never updated. Both modes therefore build their work projects in the
  same way.
- *The data-set refinements ran with two parts of twelve threads*, not case B's four parts of six.
  With four parts, `refine3D_auto`'s probability-table workers (`prob_tab_neigh`) took about 19 GB
  each on this data set and exhausted the machine's 62 GB. Both modes ran with identical settings.

**Tests.**
- *Derived weight set* (`test_work_state` in `simple_state_weight_set_tester.f90`, fast gate, sub-suite
  `state_weight_set` of `unit_project`, now 43 checks). One state of a PARTITION parent is published
  for a work project of six of its rows, in its own directory. The test checks:
  - the derived set validates, holds one state and keeps the PARTITION kind;
  - it records the parent's generation, layout digest and state;
  - its weights are the parent's column at the mapped rows, bit for bit;
  - its hard labels and population are the parent's for that state;
  - the parent's set still validates at its own generation.
- *Per-state refinement phantom* (plan section 9; in `flex_pca_blobs`, as decided above;
  `logs/tiers_p4/LastTest_flex_pca_blobs.log`, Release). FLEX's shared-memory run on the
  two-conformer phantom delivers 3 states. The truth-matched ones are state 3 (conformer A) and state
  1 (conformer B). Each is refined with `refine3D_auto state=X`, five iterations at most, from the
  same parent project. Correlation of the refined map with its truth (same frame, low-pass 8 Å):

  | State (truth) | Work project, `m_estimator=flex` | Weighted map | Hard-label map | Floor (hard − 0.01) |
  |---|---|---:|---:|---:|
  | 3 (A) | 1,332 particles, applied mass 1016.727, effective sample size 1050.455, hard population 996 | 0.98979 | 0.98723 | 0.97723 |
  | 1 (B) | 1,332 particles, applied mass 1029.948, effective sample size 1071.020, hard population 1001 | 0.99082 | 0.98622 | 0.97622 |

  - The weighted refinement includes particles of the other conformer with small weights. Its maps
    correlate slightly better with their truth than the hard ones.
  - Each weighted work project's single-state set validates and records its parent state.
  - Each final map is registered with applied mass and effective sample size, for example "FINAL
    RECONSTRUCTION REGISTERED STATE 1: hard population 996, applied mass 1016.727, effective sample
    size 1050.455".
  - The parent project is byte-identical after the four refinements.
  - The Debug build gives the same picture (`logs/p4dbg3_blobs.log`: 0.98980 against 0.98722, and
    0.99080 against 0.98627).
  - The FLEX measurements of the gate are those of Phase 3: shared and distributed accuracy 0.99950
    and 0.99750, agreement 99.80 %, mode map correlations 0.99995 and 0.99998, and the same per-state
    truth tables.

**Data-set exit** (`p4_dataset.sh`, `p4_refine.sh`, Release build; `logs/p4_dataset.log`,
`logs/p4_flex.log`, `logs/p4_refine_ds.log`, `logs/p4_ds_flex.log`, `logs/p4_ds_hard.log`).
- *FLEX.* The Phase 0 data set is the kept single-state beta-galactosidase project,
  `keep/single/1_solve3D`. On a fresh copy, `flex_pca` ran with the settings `refine3D_states` gives
  it when it initializes from FLEX:
  - `npreimages=3` and `rec_backend=gridding`, from `run_flex_pca`;
  - `refine3D_states`' own defaults that reach the in-process `flex_pca`: `min_state_frac=0.1`,
    `objfun=euclid`, `sigma_est=global`, `ml_reg=yes`, `filt_mode=nonuniform_lpset`, `automsk=no`,
    `envfsc=no`, `lpstart=12`, `lpstop=6`, `nsample=1500`;
  - case B's `pgrp=d2 mskdiam=180 nparts=4 nthr=6`.

  It ran in 5:38 and published a 2-state PARTITION weight set: state 1 with applied mass 2023 and
  state 2 with 1748 (state 0: 1,742 particles).
  - **The weights are 0/1 on this data set.** The population floor (`min_state_frac`) places exactly
    the requested number of hard-labelled states, so applied mass, effective sample size and hard
    population coincide.
  - A one-off run without the population floor also published effectively hard weights (applied
    mass equal to hard population, 4705 and 405; `logs/p4_flexk.log`). It was stopped and not
    pursued.
  - On this data set the comparison therefore tests the weighted machinery with 0/1 weights.
    Fractional weights are tested by the phantom above.
- *Refinements.* State 1 has the most applied mass. `refine3D_auto state=1 m_estimator=flex` and
  `refine3D_auto state=1` ran one after the other, each from a fresh copy of FLEX's project, with
  `pgrp=d2 mskdiam=180 nparts=2 nthr=12` (the memory decision above). Machine load: one other light
  test in parallel during the first.

  | Run | Wall time | FSC 0.5 | FSC 0.143 |
  |---|---:|---:|---:|
  | `m_estimator=flex` | 13:24 | 4.53 Å | 3.98 Å |
  | hard labels | 13:47 | 4.53 Å | 3.89 Å |

  - The weighted resolution is 1.023 times the hard one (bound 1.03).
  - The inputs are identical here (0/1 weights), so the 0.09 Å difference is the refinement's own
    run-to-run variation.
- *Parent unchanged.* The md5 of FLEX's project (`e98e42a91f00e7b3b6732a0a55af44e2`) and of each
  fresh copy is the same before and after both runs.
- *The final native reconstruction*, quoted from `logs/p4_ds_flex.log`:
  - "REFINE3D_AUTO STATE 1 WEIGHT SET: particles 2023, applied mass 2023.000, effective sample size
    2023.000, hard population 2023";
  - "FINAL RECONSTRUCTION WEIGHTS: state weight set generation 1, layout digest
    -1191944246994513173, producer flex_pca, parent state 1";
  - "RECONSTRUCTION WEIGHTS: state weight set of kind partition, generation 1, producer flex_pca";
  - "STATE 1 APPLIED MASS 2023.000 EFFECTIVE SAMPLE SIZE 2023.000 HARD POPULATION 2023";
  - "FINAL RECONSTRUCTION REGISTERED STATE 1: hard population 2023, applied mass 2023.000, effective
    sample size 2023.000".
- *Route coverage.* Both data-set runs and the phantom took the final reconstruction's direct route
  ("reusing committed canonical sigmas"; box 256 is not downscaled). The bootstrap route receives
  `m_estimator` through the same child command line, and its children copy it. That was checked in
  review, not at run time.

**Other checks (8.1).**
- *Builds.* Debug and Release build without compiler warnings. Fast gate green.
- *Tiers* (Release, `logs/queue_tiers_p4.log`): every entry passes: fast, `lib_reconstruction`,
  `lib_cart_align3D`, `lib_heterogeneity`, `lib_single`, `lib_stream`, `pcg_recon`,
  `simulate_particles`, `simulated_workflow_1jxy`, `simulated_workflow_6vxx`, `solve3D_addon`,
  `cont_refine3D_1jxy`, `flex_pca_blobs` (now 463 s, with the per-state refinement).
- *`reconstruct3D` references.* All four Phase 0 paths are byte-identical to their references: 24, 32,
  16 and 16 files (`logs/refs_check_p4.log`).
- `scripts/check_test_registry.py`, `scripts/check_descr.py` and `scripts/check_flex_dag.py` pass. No
  file mode changed. No `THROW_HARD` or `THROW_WARN` argument is continued with `//&`.
- *Review.* A read-only review found no blocking bug. Fixed from it:
  - `state=` now requires `mkdir=yes`; with `mkdir=no` a second state's run in the same directory
    would have superseded the first state's set;
  - state X's map and FSC are read from the parent, because the selection drops the entries of states
    none of the selected rows is labelled with;
  - a missing map file falls back to reconstruction;
  - the run directory's copy of the parent is removed;
  - uninitialised identity variables in `calc_final_rec`;
  - a non-short-circuit index in the gate.

**Recorded, not changed.**
- *GUI communicator in nested commanders.* The GUI communicator's mutexes are process-wide, and a
  nested commander that runs its own communicator (the selection program) destroys them when it ends.
  `refine3D_auto` therefore builds its work project before its communicator starts. `solve3D`'s
  state continuation runs the selection program after its communicator has started; this predates
  the run, and its later metadata calls may fail the same way. Not changed: outside this phase's
  scope beyond the `state=` fix.
- *Other state-associated entries in the work project.* The work project keeps the project-wide
  `vol_msk` entry, and earlier releases' `flex_weights` entries, which are ignored.
- *Manifest version 1 is no longer read* (release 4 rule). The weight sets published by Phases 1 to 3
  in scratch runs are not readable by the new code; none is kept.
- *Within the work project*, ptcl3D `state` is 1 for every selected row, while the work set's
  hard-label flags mark only the rows labelled X in the parent. The work set's hard population is
  therefore the parent's label count of X, as section 6 asks (`pop` for reporting).
- *The data-set comparison uses two parts of twelve threads.* Four parts overran memory in
  `prob_tab_neigh` (about 19 GB per worker). This is a property of `refine3D_auto` on this data set,
  not of this phase's code; it was not investigated further.

**Exit items.**
- *FLEX publishes the weights from the kept single-state project:* done (2 states).
- *For the state with the most applied mass, `refine3D_auto state=X m_estimator=flex` reaches an FSC
  0.143 resolution at most 1.03 times that of `refine3D_auto state=X` on the hard labels:* met, 3.98 Å
  against 3.89 Å (1.023), with the 0/1 weights FLEX publishes on this data set.
- *The per-state refinement phantom of section 9 passes:* met, in `flex_pca_blobs` (Release and
  Debug).
- *The parent project is byte-identical after both runs:* met (md5 before and after; checksum check in
  the phantom).
- *The final native reconstruction logs and registers the weighted map with applied mass and
  effective sample size:* met (log lines above; `mass` and `ess` on the `vol` entry, checked by the
  phantom).

Files changed in Phase 4:
- modified: `production/CMakeLists.txt` (comment of `flex_pca_blobs`),
  `src/fileio/simple_state_weights_file.f90`,
  `src/main/commanders/simple/simple_commanders_refine3D.f90`, `simple_commanders_solve3D.f90`,
  `src/main/commanders/test/simple_commanders_test_highlevel_flex.f90`,
  `src/main/project/simple_sp_project.f90`, `simple_sp_project_out.f90`,
  `simple_state_weight_set.f90` (and its tester), `src/main/simple_final_rec.f90`,
  `src/main/strategies/parallelization/simple_refine3D_strategy.f90`,
  `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`,
  `src/main/strategies/search/simple_strategy3D_matcher.f90`,
  `src/main/ui/simple/simple_ui_refine3D.f90`, `src/main/ui/simple_test/simple_test_ui_highlevel.f90`
  (description of `flex_pca_blobs`), and this plan (one file-table row).

### Phase 5: FLEX reuse (2026-10-07 and 08, written; builds and tiers pending)

The run was stopped after Phase 4 (section 3, rulings of 2026-10-07). Phases 5 and 6 were written by
hand on the Mac, on top of the Phase 0–4 diff applied to master `92d684c6f`, one diff per item for the
maintainer to build and bisect. Nothing was compiled or run by the agent: the builds, the fast gate,
the tiers and the FLEX A/B of 8.2 are the maintainer's.

**What changed, per item.**
- *Sigma2 preflight and commander shape.* `flex_pca` workers run under a worker commander of the
  private executable (`commander_flex_pca_worker`); the strategy factory refuses `part=` on the public
  command. The shared-memory and master strategies prepare the canonical sigma2 state before
  `gen_job_descr` through the project gateway's `ensure_canonical_sigma_state`, which calls
  `ensure_sigma2_for_iteration` with the run's `sigma_est` (so the global-sigma fallback for an empty
  half reaches the workers). The master reports its budget with `set_master_num_threads`; FLEX no
  longer boosts or overwrites `params%nthr`, and `set_worker_key` is gone.
- *Plane cache and embedding artifact.* `simple_ptcl_cache` exports `ptcl_cache_dir`,
  `ptcl_cache_run_token` and `ptcl_cache_space_ok`; the plane cache uses them (cache directory, a
  quarter of the free space, a build name with the run token, the contract record written last,
  then an atomic rename; stale caches deleted; version 3, the project path hashed with FNV-1a). The
  embedding artifact is one file written under a temporary name, flushed and atomically renamed, with
  a trailing block (magic `SIMPLFXD`) holding the noise scale and, after the deconvolution, the
  deconvolved coordinates, precisions and labels (format version 5). A states-only resume reads the
  noise scale from the artifact. The artifact is registered in the out segment of the run's project
  (imgkind `flex_embedding`) through the `oris` setters.
- *`prep_imgs4rec`.* FLEX's particle preparation is `prep_imgs4rec` (in `simple_matcher_3Drec`) with
  optional band, observation-model and whitening arguments, plus the public `gen_rec_plane` for the
  plane-cache path. The unused masking path of the FLEX preparation was deleted with its constants.
- *Insertion and gather.* `simple_reconstructor` gains `exp_samples` (the KB window corners,
  normalized separable weights and Friedel flags of a sample set, with `new`, `new_plane`, `gather`
  and `put_plane`; buffers grow only, so a set kept per thread allocates nothing per particle) and
  `insert_planes_multi` (K target volumes, a shared, diagonal or packed-pair density, direct sample
  locations). `project_fplane` is the K = 1 case. Deleted: `project_polar`, `interp_cmat_exp`,
  `interp_rho_exp`, `insert_plane_oversamp_opt` and FLEX's banded projection. The FLEX mean and basis
  projection, the polar bank's volume reads, the hybrid exact statistics and the M-step insertion use
  them; the E-step and the mean-scale estimate keep one sample set per thread.
- *PCG lattice layer.* `src/main/volume/simple_pcg_lattice.f90` holds what both PCG operators had
  identical: padded-lattice geometry, KB window, colouring stride, wrap table, deposition envelope and
  the solve support (`set_window_sphere`, `set_window_volume`, formerly `set_mask` and
  `set_mask_volume` in `reconstructor_pcg`). `reconstructor_pcg` and `flex_pcg_t` extend it (a
  deviation from the composed component the plan named: extension shares the code without
  forwarding routines, and the type is small). FLEX's padded deposition array and its experiment
  toggle `PCG_HARD_SOLVE_SUPPORT` (always on) are gone. FLEX's white-box self-test moved into the
  submodule `src/main/flex/fit/simple_flex_pca_pcg_tester.f90` (built only with tests); the test
  environment policy, section 4.3, now says so.
- *Clustering modules* (`src/utils/clustering`, section 7.1). `simple_hac` replaces
  `simple_avg_linkage` (average or complete linkage, stop at a count or a distance threshold, a mask
  of entries that stay singletons, the merge history; distances in double precision); FLEX's merge
  uses it on one minus the map ratio, and its union-find is gone. `simple_kmeans` (seeded at the
  point nearest the mean, then farthest-point; Lloyd to convergence; empty clusters reseeded at the
  worst-fitted point), `simple_kcenter` (farthest-point from the point farthest from the mean),
  `simple_gmm` (tied covariance only, by the ruling; floors, respawn, BIC, ICL, pairwise separation)
  and `simple_xd_gmm` (extreme deconvolution with the data passed per call rather than copied,
  because R and N are d x d per particle; `xd_select_k` is the held-out ladder). `logsumexp` is in
  `simple_stat`; `farthest_point_seeds` and `equal_mass_quantile_start` in `simple_clustering_utils`,
  whose unused `labels2smat`, `aggregate`, `silhouette_score`, `DBIndex` and `DunnIndex` were
  deleted. FLEX keeps thin adapters (`kmeans_latent_targets`, the diffusion embedding in front of
  k-center, `gmm_state_weights`, the deconvolution's calibration and projection). The mixture prior's
  initialisation (`mcfa_init`) uses `simple_kmeans`. Testers: `hierarchical_clustering`, `k_means`,
  `k_center`, `gaussian_mixture` and `extreme_deconvolution` in `unit_numerics`, plus a log-sum-exp
  check in `statistics`.
- *Statistics, linear algebra and FSC helpers.* A double-precision Cholesky family in `simple_linalg`
  (`cholesky`, `chol_forward`, `chol_backward`, `spd_inverse`, `spd_logdet`, with a tester); generic
  double-precision `median`, `mad` and `selec`; `kish_ess` in `simple_stat`; `jacobi` with a pinned
  sign for the equal-occupancy axis; the cross-fit FSC artifact (version 2) without its write-only
  fields, and the resolution reported by `get_resolution_at_fsc`; the mean-scale broadcast through
  `arr2file`. `simple_finch` was unused and is deleted.

**Expected behaviour changes, for the FLEX A/B of 8.2.**
- Rounding level: the coupled insertion samples at direct locations (as the plan's row says), the
  mean projection uses the separable gather, the Cholesky routines subtract term by term,
  `logsumexp` normalizes by `exp(x - lse)`, and the split-half reliability is computed by
  `pearsn_serial_8` (single-precision result).
- State numbering: clusters and mixture components are labelled by decreasing population (mixing
  proportion for extreme deconvolution), so states can permute. The merge's ties follow the library's
  order. `mcfa_init` now seeds at the point nearest the mean and runs to convergence instead of
  twelve iterations. The class expansion program also calls `kmeans_latent_targets`, so the order of
  its two children can change.
- The cross-fit log reports the last passing shell instead of the first failing one.

**Review.** Three read-only reviews of the item diffs found two build breakers and one failing test,
fixed in their items: the pair merge still passed the removed count argument of
`crossfsc_harvest_h` (helpers item); the targets module dropped the `matinv` import it still needs
(k-means item); the log-sum-exp shift test was tighter than double precision allows (mixture item).
Also fixed in their items: the mixture tester starts on the x axis (a start offset in y fell into
the y-split optimum in about one draw in ten of a reimplementation), the dead quantile start of
GMM AUTO's multi-seat macro-clusters, the unused `spd_solve` and the unused distance output of the
seeding, the aliasing of `native_deposition_envelope` (now a function), a strided window argument in
the insertion's inner loop and the per-particle allocations of the mean projection. A separate
review-fix diff on top of everything deletes `image%norm_noise_mask_pad_fft` (its only caller was
the deleted FLEX masking path), publishes the embedding artifact and the plane cache with
`simple_atomic_replace` (`simple_rename` deleted the destination first), and corrects two comments.

**Checks.** `scripts/check_test_registry.py`, `scripts/check_descr.py` and
`scripts/check_flex_dag.py` pass. No file mode changed. No `THROW_HARD`/`THROW_WARN` argument
continues with `//&`, and no new condition tests an array element beside the flag that guards it.

**Rulings after the review (2026-10-08), implemented.** `PCG_WEIGHT_THRESHOLD` stays 0.01; on PCG,
rows with a weight between 0 and 0.01 are aligned by `refine3D_auto state=X m_estimator=flex` but not
reconstructed. The plane cache is deleted by the run that built it once the master has its embedding
(`flex_plane_store%release`, which also frees the resident planes). `refine3D_states` refuses
`m_estimator=flex`, `refine3D` refuses it with more than one state, and `refine3D_states` withdraws
the weight set its `flex_pca` initialization published (`discard_state_weight_set` in
`simple_state_weight_set`: the `out` entry first, then the files).

**First builds and runs by the maintainer (2026-10-08).** The compiler found three names that
differ only in case (an argument `R` and a local `r` in the extreme-deconvolution tester; the
function `mean_scale_fname` and the constant `MEAN_SCALE_FNAME`), fixed in their items. The fast
gate failed one check of the weight-set withdrawal: writing a segment inside a project file does
nothing for an empty table, so a withdrawal that emptied the `out` segment left the entry on disk.
Emptying a segment is now an explicit call (`binoris%empty_segment_inside`,
`sp_project%clear_segment_inside`, with a test in the `binoris` suite). The first `flex_pca_blobs`
run stopped at the second probe iteration: the per-thread sample sets of the E-step were allocated
every iteration but not released with the other per-iteration buffers; fixed in the insertion and
gather item. The remaining `goto` statements of the touched code (FLEX's merge and embedding, the
copied quickselect, one in the STAR export) were replaced by structured control flow.
The report for the FLEX developer is `flex_on_simple_machinery_report.md`, beside this plan.

**Open, for the maintainer.**
- The registered embedding entry has no reader yet; the resume still reads `infile=`.
- `simple_kmeans`, `simple_kcenter` and `simple_xd_gmm` import `simple_clustering_utils` for the
  seeding, which pulls the distance-matrix algorithms into their compile chain.

### Phase 6: documents (2026-10-08, written)

- Amended to match the code: `doc/algorithms/reconstruction.md` (weighted state reconstruction; the
  density floor applies to any fractional weight set, not only FLEX's),
  `doc/algorithms/heterogeneity_analysis/flex_pca.md` (state maps through the service, the merge
  gates as the code runs them, the clustering modules), `src/main/flex/README.md`,
  `doc/policies/3D/reconstruct3D_pcg_policy.md`,
  `doc/policies/3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`
  (invariant 3 under a weight set), `doc/policies/3D/refine3D_policy.md` (per-state M-estimation),
  `doc/policies/heterogeneity/refine3D_states_policy.md` (the handoff),
  `doc/policies/importance_sampling_fractional_update_policy.md` (masses instead of counts),
  `doc/policies/2D/particle_cache_policy.md` (the shared cache machinery),
  `doc/how2s/how2_process_heterogeneous_datasets.md` and `doc/policies/test_environment_policy.md`;
  stale code comments in `simple_rec3D_pcg_strategy.f90`, `simple_nu_filter.f90` and
  `simple_commanders_flex_pca.f90`.
- Stale-reference scan (8.1): the procedures, modules, types, parameters and keys that exist at
  `HEAD` and no longer exist anywhere in `src/` or `production/` were searched for in `doc/`
  (excluding completed, rejected and history notes and the generated indexes), `.github/skills`,
  READMEs and scripts. No living text names one. `vol_flex` remains only as the documented
  read-only fallback for projects of earlier releases. The code map is regenerated by the build;
  the tracked Fortran indexes under `doc/code_overview/fortran-indexes/` still list deleted modules
  until `scripts/gen_fortran_indexes.pl` is run.
- Skills: not edited. `simple-main-nu-filt` still lists `flex_pca` as a caller of the nonuniform
  filter, which it no longer is; a proposed update is with the maintainer.
- This plan moves to `completed/` with its report once the maintainer's builds, tiers and A/B pass.
