# FLEX on SIMPLE machinery: fractional state weights and reuse

Date: 2026-10-07. Status: planned. Owner: the maintainer.

This note is the maintainer's target. Implementing agents record progress in section 10 and
nowhere else; the goal, rulings, contracts and exits change only by a maintainer ruling.

## 1. Goal

1. **Fractional state weights become a SIMPLE capability.** A particle may contribute to every
   state with a weight `w_is` in [0,1]. The weights live in one file per state, registered in
   the project. The standard reconstruction (gridding and PCG, shared memory and distributed,
   trailing) and `refine3D_states` read them. This is what a different M-estimator needs,
   starting with FLEX state weights driving the refinement.
2. **FLEX reuses SIMPLE machinery wherever it feasibly can.** What stays in `src/main/flex` is the
   method itself: the probe fit and its E-step, the coupled M-step, the latent deconvolution,
   state placement and the merge. FLEX stops reconstructing maps, storing weights, writing part
   formats and running clustering code of its own where SIMPLE has an equivalent.

Not a goal: new scientific features. Recomputing responsibilities during refinement is
designed here (section 4.4) but lands only by a separate ruling.

## 2. Where we are

**History.** Soft multi-state reconstruction existed in 2020-2021: each particle carried several
weighted orientations with their own states (`insert_planes_2(state=)`, `grid_ptcl_2`,
`reconstruct3D rec_soft=yes`). The search never produced mixed-state weights, because
`states_reweight` kept only the state with the largest summed weight. The code went in
`e8779a629`, `3245e505d`, `4887ae235` (2021-07, "removal of probabilistic orientation assignment
for improved performance and cleaner code") and `7da020aac` (2022-02). A per-particle insertion
weight (`pwght`, keys `ptclw`, `cavgw`) multiplied both the backprojection and the density until
`83571b305` (2026-05-24); its project slot is now `I_RETIRED_W`. Soft state weights were planned
for the PCG operator ("weight both B and D", `doc/implementation_notes/completed/pcg_priors_history.md`)
but never built.

**Standard reconstruction today** is hard-state only. `insert_plane_oversamp` takes no weight
(`src/main/volume/simple_reconstructor.f90`). `calc_3Drec` groups particles by `(state, half)` and
`update_state_half_rec` stops on a mixed-state batch (`src/main/strategies/search/simple_matcher_3Drec.f90`).
The PCG worker selects by hard state (`simple_rec3D_pcg_strategy.f90`), but `prep_particles`
already takes a per-particle noise spectrum (`simple_reconstructor_pcg.f90`), so dividing it by
`w_is` weights both the data and the density. Partial accumulators are linear, so weighted parts
sum exactly on the master.

**FLEX today** reconstructs its weighted state maps itself: a multi-target gridding insert
(`insert_planes_oversamp_multi_scaled_batch` in `simple_flex_reconstructor_latent_ops.f90`, used
only by its gridding backend), a PCG backend that re-reads particles with `sig2/w` above a 1e-3
floor, its own part files, a shell density floor, and its own delivery filter. It stores the
weights in `flex_weights_state_NNN.bin`, a 1,265-line clone of the sigma2 state store with
candidate/commit transactions that only a single master ever writes. `refine3D_states` reads
FLEX's hard labels and `vol_flex` maps; no production code reads the weights.

**The project file cannot hold per-particle vectors.** Particle records are fixed 64-slot real
records; keys outside the named slots are dropped on write (`simple_binoris_tester.f90` asserts
it). The precedent for per-particle vectors is a side file registered in the `out` segment and
tied to the particle layout by a digest (the sigma2 state, the FLEX weights).

## 3. Rulings and open decisions

Rulings (2026-10-07):

1. FLEX reuses existing SIMPLE machinery wherever it feasibly can.
2. Fractional state weights are needed, for a different M-estimator and for using FLEX state
   weights in the refinement.
3. Weights are written as one weight distribution per state, outside the particle records.

Open decisions, with proposals:

| # | Decision | Proposal |
|---|---|---|
| D1 | Pose model under fractional weights | One shared pose per particle for every state, the `refine3D_states` scaffold. Per-state poses would need a per-particle, per-state pose store; not now. |
| D2 | Where the weights come from during refinement | Phase 4 holds them fixed (the FLEX weights are the M-step weights). Recomputing them each iteration (section 4.4) is Phase 5, by separate ruling. |
| D3 | Weight kinds | Keep the two kinds the store has: PARTITION (rows sum to 1, responsibilities) and KERNEL (FLEX's equal-occupancy and single-axis bumps). Refinement accepts PARTITION only; reconstruction accepts both. |
| D4 | Density floor for fractional weights | Apply `floor_rho_shellwise` in the gridding path whenever weights are fractional; never on hard runs, which stay bitwise. PCG keeps its own shell-floored preconditioner. |
| D5 | Activation | A new parameter `state_weights=yes|no {no}` on `reconstruct3D` and `refine3D_states`. With `no`, a project that carries weight files reconstructs exactly as today. |
| D6 | FLEX probe and embedding part codecs | Keep them: they carry FLEX sufficient statistics that no SIMPLE file holds, and SIMPLE has no generic summed-array part file. No new generic facility in this round. |
| D7 | Plane cache versus `simple_ptcl_cache` | Keep FLEX's plane cache. The 2D cache would move the taper from the full box to `box_crop` (a numerical change), needs a change to `doc/policies/2D/particle_cache_policy.md`, and is deleted on exit, which breaks `infile=` resumes. Revisit with measurements. |
| D8 | UMAP and its figures | Keep in FLEX (no SIMPLE equivalent, no other user). Default `umap=no`; it is a figure, not a product. |
| D9 | FLEX state-map box | Reconstruct state maps at the box the standard path uses for the run (`box_crop` for trial maps, native for delivered maps, as rec3D does), dropping `flex_rec_box`/`flex_rec_smpd` and the temporary `smpd_crop` override. |

## 4. Fractional state weights in SIMPLE

### 4.1 The weights file

`src/fileio/simple_state_weights_file.f90` and `src/main/project/simple_state_weights.f90`
replace `src/fileio/simple_flex_weights_file.f90` and `src/main/flex/states/simple_flex_weights_state.f90`.

- One file per state, `state_weights_NNN.bin`, registered in the project's `out` segment with
  `imgkind='state_weights'` and the state index (today `add_flex_weights2os_out` with
  `flex_weights`).
- Header: magic, version, kind (PARTITION or KERNEL), producer (FLEX or EXTERNAL), particle
  count, the `ptcl3D` layout digest the sigma2 state already uses, and the per-state scalars
  (weight mass, effective sample size, hard population).
- Payload: one real per particle, zero outside the active selection.
- Written whole by one process to a temporary name and renamed. No candidates, generations or
  commits: there is never more than one writer.
- Readers check magic, version, count and digest against the current project and fail loudly.
- `state` in `ptcl3D` stays the selection flag and the maximum-weight label.
- Old projects with `flex_weights` entries must still load (`.simple` files are the one
  backwards-compatibility exception); the entries are ignored, their files not read.

### 4.2 Reconstruction

**Gridding.**
- Restore an optional weight on `insert_plane_oversamp` that scales both the data term and the
  density term, as before `83571b305`.
- `calc_3Drec` builds each `(state, half)` membership list from `w_is > W_FLOOR` and inserts with
  `w_is`. One half-map reconstructor at a time, as now: memory stays independent of `nstates`.
  A particle is read once per state it weighs into.
- Workers write the same partial names and payloads; the master sums them unchanged.

**PCG.**
- The worker selects `(state, half)` members by `w_is > W_FLOOR` and divides the particle's noise
  spectrum by `w_is` before `prep_particles`, which weights both B and D. This is FLEX's existing
  method, moved into `simple_rec3D_pcg_strategy.f90`.
- Master reduction and raw-accumulator files are unchanged. The ridge stays the production one.

**Trailing reconstruction and update fractions.** Every count that today counts rows of a state
counts summed weights instead: `population_blend_weights` and the trail counts in
`simple_commanders_rec_distr.f90`, the state update counts and fractions in `simple_oris_getters.f90`,
the PCG master's representation counts, the seed population, and the add-on mode's frozen
accumulators. The sampler is untouched: it never reads state weights (sampling is driven by
class averages only).

**Population gates** that skip empty states compare the summed weight, not the hard count.

**Hard runs are bitwise unchanged.** With `state_weights=no`, no code path changes. With
`state_weights=yes` and 0/1 weights equal to the hard labels, results equal the hard run bit for
bit; this is the main oracle (section 7).

`W_FLOOR` is 1e-3, FLEX's current value.

### 4.3 Refinement with fixed weights (Phase 4)

`refine3D_states state_weights=yes` aligns as today, with the hard label selecting the
particle's reference, and reconstructs every state with the fixed weights. The weights come from
FLEX (PARTITION kind). FLEX already initializes `refine3D_states` when the project carries only
state labels 0 and 1; it now also writes the weights the refinement uses.

### 4.4 Recomputed responsibilities (Phase 5, separate ruling)

Under the Euclidean objective the probability tables hold noise-normalized negative
log-likelihoods for every active state (`loc_tab(nspace*nstates, nptcls)` in
`simple_eul_prob_tab.f90`, merged by `read_tab_to_glob`). With `d_s` the best distance of state
`s` over the evaluated candidates, the responsibilities are `r_is ∝ π_s exp(-(d_s - d_min))`. The
master writes them to the weights files after the merge and the next reconstruction uses them.
The `prob_state` retirement removes the fixed-pose state table; this phase uses the pose-search
table instead and must confirm every active state is evaluated in the chosen `prob_neigh_mode`.

## 5. What FLEX reuses

| FLEX piece (lines) | SIMPLE target | Lines away | Risk | Phase |
|---|---|---:|---|---|
| State backends, state parts, rec3D driver, delivery (~1,250) | Standard weighted reconstruction (4.2); delivery through `evaluate_halfmap_pair`, `simple_butterworth`, `fsc2optlp_sub`, `mask3D_soft` | 1,000-1,250 | Delivered maps change: FLEX's "8th-order" filter is really a squared 4th-order one, and the PCG ridge becomes the production one | 3 |
| Weights store (1,265) | `simple_state_weights` (4.1) | ~900 | Low | 1 |
| Multi-target and coupled insertion, three "one window, many volumes" gathers (~700 in latent ops, polar, E-step) | One batch insertion into K accumulators with per-target weights and a density mode (own, shared, diagonal, packed pairs), and one multi-volume point gather, in `simple_reconstructor.f90`; `project_fplane` becomes the K=1 case; the unused `project_polar` goes | 550-700 (adds ~350 in the reconstructor) | refine3D hot path: keep `insert_plane_oversamp` separate until measured. The coupled insert moves from incremental to direct sample locations: rounding-level change | 6 |
| `prep_imgs4projected_model` (~75) | `prep_imgs4rec` with optional band, transfer-plane and observation-model arguments | 70-90 | Optional arguments only on the refine3D path | 6 |
| Coupled PCG operator: geometry, support window, kernel and right-hand-side accumulation, scatter, fold, kernel finalisation (~800 duplicated) | A shared lattice layer composed by `reconstructor_pcg` and the coupled operator; the K×K pair kernels and the per-voxel K×K preconditioner stay FLEX | 350-400 each side | Bitwise changes if loop order moves; share only code that is identical today | 6 |
| Embedded PCG self-test (~600) | `simple_flex_pcg_tester.f90` | 0 (moved) | None | 6 |
| Sigma2 preflight in the commander (~80) | `canonical_sigma2_consumable` and `simple_sigma2_bootstrap` (small extension: keep per-stack grouping), run before `gen_job_descr` as reconstruct3D does | ~80 | The global-sigma fallback for an empty half must move with it | 6 |
| Commander and strategy (worker strategy, worker-key injection, thread boost) | The standard shape: a separate worker commander routed by the private executable, `set_master_num_threads`. The rounds callback stays, because FLEX has nested EM loops | 110-130 | Thread sizing differs; FLEX must stop overwriting `params%nthr` | 6 |
| Embedding cache (184) | Registered in the `out` segment, written atomically, path remapped with the project | ~30 | Low (gains resume correctness) | 6 |
| Mean-scale and round-weights broadcasts | `arr2file`; the round weights vanish with the state backends | ~100 | Low | 3, 6 |
| Clustering and statistics: robust median/MAD, power iteration, complete linkage, diffusion k-center, k-means, GMM core, the SPD Cholesky family, log-sum-exp, farthest-point seeding, union-find | `simple_stat`, `jacobi_dp`, `simple_avg_linkage` (extension), `simple_diff_map_graphs` (extension); the generic ones move to `src/utils` as shared tools | 300-500 | Tie-breaking and eigenvector sign can reorder states: pin the sign | 6 |
| Cross-fit FSC helpers | `get_find_at_crit`, `fsc2optlp_sub`, `get_resolution_at_fsc`; drop the write-only artifact fields | ~100 | The reported band moves (first failing versus last passing shell) | 6 |

Stays in FLEX, genuinely new: the probe fit and its E-step, the posterior, the coupled solve,
extreme deconvolution, state placement and the two-gate merge, the polar bank's variable-length
rings and their quadrature contract with the exact low-k part, UMAP, and the cross-fit
component matching. Kept by decision: the probe and embedding codecs (D6), the plane cache (D7).

Two findings to check in Phase 0, not refactor items: the polar bank snaps particles to bank
directions without symmetry copies, so a non-c1 particle outside the asymmetric unit may snap
to a distant direction; and the resident plane store allocates planes on the padded lattice but
fills every second sample, so it holds about four times the memory it uses.

## 6. Phases

Each phase is one or a few commits by the maintainer and lands only with its exit met.

| # | Phase | Scope | Exit |
|---|---|---|---|
| 0 | Preconditions and baselines | The current FLEX working tree compiles and its tiers are green (including the open review fixes: the commander's missing `simple_abspath`, the repeated envelope builds, the state-stage cache adoption). The `prob_state` retirement diff (Dell run `prob_state_flex`, finished 2026-10-07) is merged. Record: hard `reconstruct3D` outputs (gridding and PCG, shared memory and two parts) as bitwise references; the FLEX application phantom and workflow-gate state maps against truth; `refine3D_states` on the `prob_state_flex` Phase 0 data set (best state 4.03 Å). Check the two findings above. | References and observations in section 10 |
| 1 | Weights file in SIMPLE | 4.1. FLEX writes the new file; nothing reads it yet | Round trip exact; a digest mismatch, wrong count or truncated file is refused; an old project with `flex_weights` entries loads |
| 2 | Weighted reconstruction | 4.2, gridding and PCG, shared memory and distributed, trailing counts, population gates, `state_weights` parameter with UI and parser | `state_weights=no` bitwise equal to the Phase 0 references; `state_weights=yes` with 0/1 weights equal to the labels bitwise equal too; two parts equal shared memory bitwise; a fractional two-state phantom recovers both truth maps |
| 3 | FLEX on the standard reconstruction | Delete FLEX's state backends, state parts, delivery and round-weights broadcast; its trial maps (bandwidth selection, merge) call the standard weighted reconstruction; D9 | Application phantom and workflow gate: labels and state maps against truth no worse than the Phase 0 observations by the margins recorded there |
| 4 | Refinement with fixed weights | 4.3 | On the Phase 0 data set the best state is within 3 % of the hard run; on a phantom with known fractional membership each state map correlates with its truth no worse than the hard run |
| 5 | Recomputed responsibilities | 4.4, by separate ruling | Defined in the ruling |
| 6 | FLEX reuse items | The section 5 rows marked 6, one small diff each, in the order: sigma2 and commander shape; `prep_imgs4rec`; insertion and gather; PCG lattice layer; clustering, statistics and FSC helpers; embedding cache | Each diff: its suites green and FLEX latents A/B against the previous step (bitwise where the row says so, otherwise within the PCG solver tolerance) |
| 7 | Documents | `doc/algorithms/reconstruction.md` (weighted variant; the density floor is no longer FLEX-only), `doc/algorithms/heterogeneity_analysis/flex_pca.md`, `src/main/flex/README.md`, the policies below, the skills | Text matches code |

Policies to amend in Phase 7 (and to read before Phase 2):
`doc/policies/3D/separate_alignment_and_reconstruction_for_multistate_peak_mem_reduction.md`
(invariant 3, "every valid selected particle belongs to exactly one (state, half) group";
invariant 4 stands), `doc/policies/3D/refine3D_policy.md` (hard assignment only),
`doc/policies/importance_sampling_fractional_update_policy.md` (row counts per state),
`doc/policies/heterogeneity/refine3D_states_policy.md`, `doc/policies/3D/reconstruct3D_pcg_policy.md`.

## 7. Test net

Independent oracles, in the fast or library tier unless stated:

- **0/1 weights equal hard labels.** For gridding and PCG, shared memory and distributed, a weight
  file of 0/1 values equal to the labels reproduces the hard reconstruction bit for bit.
- **Linearity.** Inserting a particle with weight `w` equals inserting it with weight 1 and
  scaling the accumulators by `w`, to round-off.
- **Distributed equals shared memory** for fractional weights, bit for bit (fixed part order).
- **Fractional phantom** (library tier): two conformers, each particle a known mixture; the
  weighted reconstruction recovers both truth maps; the hard reconstruction from the argmax
  labels is the comparison.
- **Weights file:** round trip, digest and count mismatch refused, old-project load.
- **Regression:** the Phase 0 `refine3D_states` run and the FLEX application phantom, against
  the recorded observations with stated margins.

## 8. Risks

- **Cost.** Gridding re-reads a particle once per state it weighs into: 1-2× for sparse
  responsibilities, `nstates`× for flat ones. PCG re-reads per `(state, half)`, as FLEX does today.
- **Numbers.** FLEX state maps change in Phase 3 (filter, ridge, box). Phase 3's exit is against
  truth, not against the old maps.
- **refine3D hot path.** Phases 2 and 6 touch `insert_plane_oversamp`, `calc_3Drec` and
  `prep_imgs4rec`. Every change is optional-argument or gated, and the hard-run bitwise oracle
  runs after each.
- **Scope creep.** Recomputed responsibilities, per-state poses and a generic part-file facility
  are out of this round unless ruled in.
- **Plan drift.** Earlier rounds rewrote their own plan to match what was delivered. Here the
  plan is fixed: a phase that cannot meet its exit stops and reports.

Out of scope, noted: `maybe_postprocess_reconstruct3D` copies a whole `parameters` per state
(`params_pp = params` in `simple_rec3D_strategy.f90`), which `doc/policies/compile_time_policy.md`
forbids.

## 9. File table

| Phase | Files |
|---|---|
| 1 | `src/fileio/simple_state_weights_file.f90` (new, replaces `simple_flex_weights_file.f90`), `src/main/project/simple_state_weights.f90` (new, replaces `src/main/flex/states/simple_flex_weights_state.f90`), `src/main/project/simple_sp_project_out.f90`, `src/main/project/simple_sp_project.f90`, the FLEX writer `src/main/flex/run/simple_flex_pca_project_gateway.f90`, testers of the store |
| 2 | `src/main/volume/simple_reconstructor.f90`, `src/main/strategies/search/simple_matcher_3Drec.f90`, `src/main/strategies/parallelization/simple_rec3D_strategy.f90`, `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`, `src/main/commanders/simple/simple_commanders_rec_distr.f90`, `src/main/ori/simple_oris_sampling.f90`, `src/main/ori/simple_oris_getters.f90`, the frozen accumulator module, `src/main/params/simple_parameters.f90` and its parser, `src/main/ui/simple/` for `reconstruct3D` and `refine3D_states`, new tests |
| 3 | `src/main/flex/states/` (state backends, gridding, PCG, delivery, rec3D, state parts deleted), `src/main/flex/run/simple_flex_pca_application.f90`, `src/main/flex/run/simple_flex_pca_state_service.f90`, `src/main/flex/states/simple_flex_pca_merge.f90`, `src/main/flex/fit/simple_flex_reconstructor_latent_ops.f90`, FLEX testers |
| 4 | `src/main/commanders/simple/simple_commanders_refine3D.f90`, `src/main/strategies/search/simple_strategy3D_matcher.f90`, the `refine3D_states` UI, tests |
| 6 | `src/main/volume/simple_reconstructor.f90`, `src/main/volume/simple_reconstructor_pcg.f90`, `src/main/strategies/search/simple_matcher_3Drec.f90`, `src/main/commanders/simple/simple_commanders_flex_pca.f90`, `src/main/strategies/parallelization/simple_flex_pca_strategy.f90`, `production/simple_private_exec_driver.f90`, `src/main/sigma2/simple_sigma2_bootstrap.f90`, `src/utils/clustering/`, `src/utils/math/`, `src/main/pca/simple_diff_map_graphs.f90`, `src/main/flex/**` |
| 7 | the documents and policies named in section 6, `.github/skills/` |

## 10. Progress

(empty)
