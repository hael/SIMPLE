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
| 0 | Preconditions and baselines | `878adf2d4` (the base of the run) builds, and its fast gate, library tier and the high-level entries of 8.1 run; which entries pass is recorded. Record the references of 8.2; the FLEX application phantom (`lib_heterogeneity`) and the workflow gate (`flex_pca_blobs`): per arm the recoded label accuracy, the state count and each truth-matched state map's correlation with its truth; `refine3D_states` on the same data set as the `prob_state` retirement run (best state 4.03 Å there). Check the two findings of section 7 and record what is true; fix nothing. | References and observations in section 12; base build and tier results recorded |
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
| `src/main/nu_filt/simple_nu_filter.f90` | the comment that names FLEX's `nufilt` report | 1 |
| `src/main/volume/simple_reconstructor.f90` | weighted insertion (5); batch insertion and multi-volume gather (section 7) | 2 5 |
| `src/main/volume/simple_reconstructor_pcg.f90` | real mass and weight-set identity in raw accumulators; shared lattice layer | 2 5 |
| `src/main/volume/simple_trail_chain_manifest.f90`, `src/main/volume/simple_frozen_accum.f90` | the manifest contract for both backends; provenance in frozen accumulators | 2 |
| `src/main/strategies/search/simple_matcher_3Drec.f90` | weighted membership lists; `prep_imgs4rec` options | 2 5 |
| `src/main/ori/simple_oris_sampling.f90`, `src/main/ori/simple_oris_getters.f90` | population rule and update counts on applied mass | 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90`, `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | population gates on effective sample size; weight set threaded through assembly | 2 4 |
| `src/main/strategies/search/simple_strategy3D_matcher.f90` | reconstruction through the service with the weight set | 2 4 |
| `src/main/commanders/simple/simple_commanders_refine3D.f90` | `vol_flex` handoff replaced by normal discovery; `refine3D_auto state=X` work project and `m_estimator` | 3 4 |
| `src/main/commanders/simple/simple_commanders_project_core.f90` | selection of a work project's rows by weight (beside selection by state) | 4 |
| `src/main/simple_final_rec.f90` | weights through the direct and bootstrap routes; mass and effective size registered | 4 |
| `src/main/flex/` | FLEX (D10, D12, Phase 3 deletions, Phase 5 reuse) | 1 3 5 |
| `src/main/commanders/simple/simple_commanders_flex_pca.f90`, `src/main/strategies/parallelization/simple_flex_pca_strategy.f90`, `production/simple_private_exec_driver.f90` | FLEX writes the weight set; worker commander and standard shape | 1 3 5 |
| `src/utils/clustering/` | the clustering modules of 7.1; `simple_avg_linkage` renamed to `simple_hac` | 5 |
| `src/main/strategies/search/simple_ptcl_cache.f90`, `src/main/sigma2/simple_sigma2_bootstrap.f90`, `src/main/pca/simple_diff_map_graphs.f90`, `src/main/pca/simple_diffusion_maps.f90`, `src/utils/math/` | shared machinery FLEX reuses (section 7) | 5 |
| `doc/algorithms/reconstruction.md`, `doc/algorithms/heterogeneity_analysis/`, `doc/policies/3D/`, `doc/policies/heterogeneity/`, `doc/policies/importance_sampling_fractional_update_policy.md`, `doc/policies/2D/particle_cache_policy.md`, `doc/how2s/how2_process_heterogeneous_datasets.md` | documents and policies | 6 |
| `src/main/flex/README.md` | the FLEX layout after the phases | 6 |
| `.github/skills/` | skills that describe the changed contracts | 6 |

## 12. Progress

(empty)
