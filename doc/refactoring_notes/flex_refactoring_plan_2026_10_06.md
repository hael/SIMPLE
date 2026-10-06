# FLEX refactoring plan: reduction, reuse, a numerical test net, and `inpl_refine`

Date: 2026-10-06. Revision 3, after the review of the same day (section 12 lists what it
corrected). Status: planned, nothing implemented.

Baseline: `master` `b818cd7cc`. Scope: the 45 Fortran files of `src/main/flex` that belong
to `flex_pca` (21,422 lines), its commander (`simple_commanders_flex_pca.f90`, 344) and
strategy (`simple_flex_pca_strategy.f90`, 405). `simple_flex_cls_expansion.f90` (1,109
lines) lives in the same folder but is a different program: out of scope except for the
style pass. Every count below is over the scoped set.

This note is the single living design record for this round: contract, measurements,
contracts, test net, phases, ledger. Update it as phases land; do not create companion
documents. It supersedes the architectural plan of 2026-09-18
(`completed/flex_pca_architecture_audit_and_refactoring_plan_2026_09_18.md`), whose
byte-identity contract is dropped here.

Validation level: static source inspection (`rg`/`awk`/python over the tree, the module
`use` graph extracted from the sources). No compilation or run. Timing statements are
structural until phase 0 measures them with the instrumentation of section 4.5.

## 1. Why another round

The September round reorganised ownership under a byte-identical contract: 36 → 46 files,
22,678 → 22,531 lines (both counts including the class-expansion file; −0.6%). It named the
boundaries and by design deleted, merged and replaced nothing. Its fixture harness (phase
0a) was never built, so no test today guards the E-step, the M-step, the probe parts or the
delivered maps: 349 tester lines cover 21,073 scoped production lines (volume: 1,438 over
8,245).

Measured against the rest of the library (`.github/skills/simple-modern-fortran`,
`doc/policies/compile_time_policy.md`):

| Measure | flex (scoped) | volume | strategies/parallelization |
|---|---:|---:|---:|
| module-level `parameter` lines | 119 | 77 | -- |
| `public ::` entries | ~250 | 62 | 77 |
| procedures / type-bound | 510 / 22% | 293 / 56% | 406 / 59% |
| files importing the umbrella `simple_core_module_api` | 41 of 45 | -- | -- |
| comment blocks of 6+ lines | 48 | 19 | -- |
| lines over 132 characters | 120 | 20 | -- |

### 1.1 Findings the phases act on

- **The E-step that runs today stays exactly as it is.** It is a hybrid: for every
  particle the Fourier shells `0..rhyb`, with `rhyb = nint(0.72·band)` derived at bank
  build (`fit_polar_bank_build`), are accumulated exactly on the Cartesian lattice
  (`polar_hybrid_exact_accum`: one KB stencil per lattice position shared across the mean
  and the `ncomp` basis volumes, no allocation, double accumulation), and the rings
  `rhyb+1..band` come from the polar shared-direction bank (`fit_estep_former_polar`).
  `l_pol_hyb` is `.false.` at stage begin and `.true.` after bank build in every run;
  `l_pol_grid` and `l_pol_bank_it` are lazy-initialisation flags. None is a selector; none
  is touched.
- **Dead code**, each item verified at zero call sites or at a pinned selector:
  `fit_estep_former_cart` (selected only by `l_pol_es=.false.`, which nothing sets); the
  retired split overrides `rhyb_req` (pinned 0) and `l_rhyb_off` (pinned false), so `rhyb`
  is always the derived value; the angular oversampling argument `osamp_pol` (pinned 1);
  the ECM contrast alternation (`n_probe_cm=0`, `nml_plain=0`, `l_probe_mls=.false.`);
  `flex_weights_commit`, `_load_all`, `_merge_local_ranges`, `_prepare_update`,
  `_range_path`; the unfused `polar_sample_particle`, `polar_sample_at_pose`,
  `polar_dir_neighbours`, `polar_apply_shift`; `invert_lower`, `punit`;
  `flex_pcg_t%get_cdim`, `%get_pijk`, `%list_signature`. Removing them changes no number.
- **Configuration outside `parameters`**: `flex_run_settings` reads 20 `SIMPLE_COV_*`
  environment keys. Every behaviour reachable only through a non-default key value is
  unsupported after this round.
- **Two live codecs for one payload**: `PROBE_PART_VERSION=12` (single fit) and
  `PROBE_PART_VERSION5=10` (paired) both serialise `flex_probe_part`; stage ids 4–7 keep
  1–3 unused; seven magic/version pairs.
- **Duplicated constants**: `MERGED_PC_FBODY`, `MERGED_META`, `MERGED_EIG_FNAME`,
  `PAIRED_MANIFEST` in both `pairmerge` and `probe_fit_update`; `FLEX_PCG_LAMBDA_REL` =
  `FLEX_PCG_LAMBDA_REL_DEFAULT`; 0.143 twice; three `reconstructor_pcg` constants copied.
  Six prefix families; 41 identifiers carry the covariance-era `cov`.
- **Ownership faults the use graph shows** (section 5.1 draws the graph): `weights`, a
  domain module, calls `reconstruct_flex_weighted_states`, mutates `params%outvol`, reads
  and deletes map files (`cv_select_bandwidths`); `fit_types`, the typed payload, imports
  `crossfsc`, `polar`, `mstep`, `run_types`; `latent_ops` imports `plane_cache` while
  `planes` imports `latent_ops`; `state_delivery` imports `project_gateway`.
- **Records and singletons instead of objects**: `planes` (10 module variables),
  `plane_cache` (8), `pcg` (a `save`d envelope image and flags), `artifacts`
  (`part_dir_prefix`); free-procedure modules over records (`crossfsc`, `weights_state`,
  `deconv`, `umap`). Four `_t` suffixes. `simple_flex_pca_pcg` takes `parameters` in six
  signatures.
- **Plane-cache lifecycle**: the header carries project rows, geometry, selection count
  and hash (`header_key`), but a worker adopts a cache on `cache=yes` and geometry alone
  (`plane_cache_in_use`: "does not need the selection"); `plane_cache_close` releases the
  unit and buffer but leaves `l_probed`, `l_avail`, `fname_glob` set, so a second session
  in the same process inherits a stale probe.
- **Parallel implementations of library machinery**: `simple_flex_pca_pcg` (2,427 lines)
  mirrors `reconstructor_pcg` with its own 164-line `cg_core`, outcome type and stopping
  constants; seven private SPD/Cholesky solvers; two particle-plane caches; k-means and
  diffusion k-centre beside `utils/clustering`.
- **Instrumentation is mixed**: the per-pass `read/prep/project+solve/insert` line is
  wall time (`tic`/`toc`), `sec_proj_thr`/`sec_gram_thr` are per-thread sums. No number
  from today's log can be compared across the two.
- **Too many files** (Afan): 45 scoped files in three folders, median 340 lines, ten under
  160. Section 5.3 sets a cohesion-driven target.
- **Style**: seven extraction scars (bodies at their old indentation, four in
  `probe_fit_estep`), seven file headers over the five-line limit, 17 comments citing
  measurements, 11 "used to" narrations.

### 1.2 Code inquiry: the M-step insertion and the E-step (Afan's questions)

**Insertion.** `insert_planes_oversamp_coupled_batch_scaled` (`latent_ops`, 198 lines)
and `flex_pcg_t%accumulate` (110 lines) follow `reconstructor%insert_plane_oversamp` (144
lines) step for step: one OpenMP region per batch, the strided-row `!$omp do
schedule(static,1)` over `h` with `LATENT_SAFE_STRIDE`, the Nyquist-disk bounds, the
`kbinterpol` stencil as `wx*wy*wz`, the symmetry loop, the Friedel half-plane rule, the CTF²
density deposit. What differs is inherent to the coupled model (`ncomp` volumes and
`ncomp(ncomp+1)/2` pair densities per sample). No algorithmic optimisation of the insertion
is planned.

**E-step loop.** `fit_estep_pass` reads planes in batches of `MAXIMGBATCHSZ`, runs
`!$omp parallel do schedule(dynamic)` over the particles of the batch with thread-private
work arrays, then the batched M-step insertion: the pattern of `simple_matcher_3Drec`.
Nothing further to batch at the loop level.

**E-step cost, structurally.** Per fit and iteration: the bank, `ndir = clamp(nptcls/40,
1000, 4000)` directions × `ncomp+1` volumes on the rings `rhyb+1..band`, plus the ring
Grams. Per particle: (a) the exact low-`k` part over the `npos` lattice positions of the
disk `h²+k² ≤ rhyb(rhyb+1)` (about half the band disk's area at the 0.72 split), one
stencil per position shared by the `ncomp+1` volumes; (b) the ring part, one fused KB
gather of the data plane, one `sgemv` (`nsamp2 × (ncomp+1)`), two `dgemv` over the `nk`
ring Grams; then an `ncomp × ncomp` posterior solve. No allocation inside the particle
loop. What only measurement settles (phase 0, section 4.5): the (a)/(b) split; the bank
build against the particle loop at small particle counts; read/prep against the cache.

## 2. Contract

Preserved for the user:

1. `simple_exec prg=flex_pca` with its current parameters and the shared-memory,
   distributed-master and worker execution modes; `inpl_refine=no` (default) is the
   current behaviour.
2. The delivered products: embedding cache, eigenvolumes, latent readouts and figures,
   state weights (`flex_weights_state_NNN.bin`), state maps (`vol_flex`), the
   `ptcl3D/state` labels and the out-segment registrations; with `inpl_refine=yes` also
   the refined `ptcl3D` poses.
3. The scientific method: the projected-model PPCA/MCFA probe fit with the hybrid E-step
   of section 1.1, the coupled gridding or PCG M-step, the paired mod-4 halves with the
   cross-fit FSC, the latent deconvolution, state placement, two-gate merge and per-state
   reconstruction.
4. `rec_backend` and `rec_states_backend` keep their separate meanings; `maxits_pcg=0`
   keeps its gridding/identity meaning.

Not preserved:

1. Byte identity with `master` outputs. Results must agree with the simulation truth and
   with the pre-refactor observation within the tolerances and margins of section 4. The
   PCG migration (phase 5) is expected to move numbers at the solver-tolerance level.
2. Artifact formats: part files, caches, the weights store and the new pose parts may
   change layout and version freely; no reader accepts an older layout (none does today).
   A run restarts on the new binary; resume reads only caches the same binary wrote.
3. The `SIMPLE_COV_*` environment keys and every behaviour reachable only through one.
4. The unreachable pure-Cartesian former and the pinned split/oversampling overrides.

Maintainer rules apply throughout (`CLAUDE.md`, `AGENTS.md`): edits by read-modify-write,
no compilation by the agent, the maintainer builds and runs the gate after each phase,
commits are the maintainer's.

## 3. Decisions the owner must make before the phases that need them

| # | Decision | Needed by | Proposed |
|---|---|---|---|
| D1 | CTest process budget: `SIMPLE_CTEST_BUDGET` is 32 and "only goes down"; the workflow gate `flex_pca_blobs` is one more highlevel process. | phase 0 | Raise to 33 with this note as the recorded reason, or retire one existing highlevel entry; the owner chooses. |
| D2 | `polarft_calc`: a thread-safe per-reference memoize (`memoize_ref(iref, ithr)`) so the exhaustive stage can use the FFT path on a per-thread reference. | phase 7 | Not required: phase 7 uses the direct per-rotation evaluation (section 6.3, stage 2). Add it later only if phase 7's timing says the direct path dominates. |
| D3 | Previous-iteration latent for the `inpl_refine` prediction: persisted artifact, two-pass current iteration, or mean-only. | phase 6 | Two-pass (section 6.1). |
| D4 | The target dependency DAG and consolidation map of sections 5.1–5.3. | phase 3 | As drawn; file count 45 → 33. |

## 4. The numerical test net

Three tiers, following `doc/policies/test_environment_policy.md`. The net is built on
`master` first (phase 0) and runs after every phase; a phase lands only when every tier is
green.

### 4.1 Terminology (policy section 2.1)

Three different things, kept apart in every test and in the ledger:

- **Expected value from an independent oracle**: simulation truth, a closed form, a brute
  force written in the test, a numpy emulation whose constants are pinned with a comment
  naming the script. Its **tolerance is derived from the arithmetic** (single precision,
  KB interpolation order, quadrature, a statistical standard error).
- **Observation on the pre-refactor code**: a number the phase 0 run produced. It is not an
  expected value. It goes in the ledger as the baseline, and a test may assert a
  **regression margin** against it only where no independent oracle exists (end-to-end
  metrics), with the margin and its reason stated next to the number.
- **Floor**: the threshold a test asserts. It is one of the two above, labelled as such.

### 4.2 Fast tier: `unit_heterogeneity` (per-estimator, independent oracles)

Added to `simple_flex_pca_tester` and `simple_flex_pcg_tester`, each sub-suite well under a
second on one thread:

- Posterior solve: for a random 4×4 SPD `G`, `Gamma` and contrast `a`, mean and covariance
  against the closed form `(a²G + Gamma⁻¹)⁻¹ a b`; 1e-10 relative in double. Negative
  control: `a=0` returns the prior.
- The E-step former (hybrid): for one plane and a basis of two synthetic volumes at a pose
  placed exactly on a bank direction (snap zero), the former's `G`, `b`, `c`, `e_mm`
  against a brute-force central-section sum over the whole band written in the test
  (direct KB gather from the unexpanded volumes, direct half-plane inner product). The
  low-`k` contribution alone equals the oracle restricted to `h²+k² ≤ rhyb(rhyb+1)` at
  round-off. Linearity oracle: a noise-free particle built as `μ + Σ c_q u_q` must give
  `b = G c` and `myv = e_mm + cᵀc₀` to the interpolation tolerance. Tolerances from the KB
  order and the ring quadrature, derived in the test. (The earlier "90° rotation gives an
  orthogonal `b`" control is withdrawn: not generally true.)
- Coupled M-step on a toy that is a sub-second oracle: two 8³ blob volumes (512 unknowns
  per half), 40 projections at fixed-seed poses with known `E[zz']`; gridding and PCG
  M-steps recover the basis to a subspace principal angle under a declared tolerance, and
  the PCG solve agrees with a dense solve of the same normal equations (numpy emulation,
  constants pinned). The 16³ version moves to the library tier.
- PCG operator: symmetry `<x,Hy> = <Hx,y>` and `<x,Hx> ≥ 0` on random vectors (the masked
  unregularised operator is positive semidefinite, not definite); positive definiteness
  asserted only with the ridge applied. Kept as the one white-box self-test.
- Cross-fit FSC: two half bases related by a signed permutation (what the matcher supports,
  `crossfsc.f90:36`) give the series of the identity; a random basis gives the noise floor.
- Probe-part codec: write/reduce/read round trip of a random `flex_probe_part` is exact.
- Plane cache (with phase 4): changed selection → rebuild; changed project → rebuild;
  worker adoption of a cache whose selection does not cover the partition → refused; two
  sessions in one process in either order see a pristine store.
- PCG solver engine (with phase 5): on a dense SPD operator against a direct solve; on a
  diagonal operator convergence in one step; stop reasons and outcome fields pinned.
- `inpl_refine` kernels (with phase 7): bank-frame ↔ particle-frame pose conversion
  round-trips; the search-bank prediction `U_0 + Σ ẑ_q U_q` equals the polar projection of
  the combined volume to KB tolerance; a particle simulated from a known prediction with a
  known in-plane rotation and shift increment is recovered by the three stages to a
  declared precision; an unrelated particle is left at its seed by the acceptance guards.
- Kept: embedding cache I/O, auto settings, population floor, kernel bandwidth, state
  weights, deconvolution.

### 4.3 Library tier: `lib_heterogeneity` (in-process end-to-end on a phantom)

A new tester next to the code it tests, `simple_flex_pca_application_tester.f90`, one
sub-suite `flex PCA two-state phantom`, hermetic, deterministic, minutes, registered in
`suites_lib_heterogeneity` and in the area's `suite=` help:

- Phantom: an asymmetric sum of Gaussian blobs (`add_gaussian_blob`), box 64, conformers A
  and B differing by one blob moved 8 pixels; 2,000 particles at fixed-seed uniform
  orientations, CTF and noise at a declared SNR, labels alternating by row. Truth: the maps,
  the labels, A−B, the poses.
- Runs `flex_pca_application%run` on an in-memory project for the four backend arms.
- Asserts against **simulation truth** with arithmetic-derived tolerances where the metric
  has them: label assignment accuracy after permutation; leading eigenvolume correlation
  with A−B (sign-free); each state map docked to its truth map (`compare_to_truth`); weight
  sum rules; a one-state phantom gives no bimodality and the merge folds to one state.
  Where the metric is an end-to-end number without a closed form, the floor is the phase 0
  observation minus a stated regression margin, labelled as such.
- Asserts the invariants a refactor breaks silently: the two halves agree with the single
  fit's subspace; the in-process `run_worker` arm over two partitions reproduces the
  shared-memory arm. **This arm tests partition reduction only**; it is not distributed
  coverage (no worker processes, serialisation, project reload or part transfer), which is
  the workflow gate's job.
- `inpl_refine` arm (phase 7): poses jittered by known angles (σ 3°) and shifts (σ 2 px);
  `inpl_refine=yes` recovers the pose RMS error to a declared precision and does not lower
  the truth metrics beyond a stated margin; `inpl_refine=no` on the jittered project is the
  negative control; `inpl_refine=no` on the clean project is A/B identical to the
  pre-feature code at 1e-5.
- Records the E-step split (section 4.5) for the ledger.

### 4.4 Workflow gate: `flex_pca_blobs` (highlevel, nightly)

A commander in `simple_commanders_test_highlevel_*` (policy section 4.4), with its case in
the router, its program in the test UI, and the CMake entry with its `COMMAND` arguments:
`simple_add_test(flex_pca_blobs highlevel 3600 8 RUN_SERIAL COMMAND $<TARGET_FILE:simple_test_exec> test=flex_pca_blobs nthr=8 ...)`,
subject to D1. It simulates the two-state set with `simulate_particles` from the two
phantom maps, imports both into one project with the truth labels, runs
`simple_exec prg=flex_pca` **distributed (`nparts=2`, real worker processes)** and
shared-memory, with `inpl_refine=no` and (from phase 7) `=yes` on a jittered copy, writes
`metrics.tsv` through `simple_test_gate`, and fails through the gate. This is where
distributed behaviour is covered: part files, pose parts, project reload between rounds,
master merge, shared-vs-distributed agreement.

### 4.5 Instrumentation (phase 0, before any number is compared)

One definition per bucket, applied to both the per-pass line and the per-thread sums:
`read`, `prep`, `bank`, `exact_lowk`, `ring`, `solve`, `insert`, and with phase 7
`search`, each reported as **thread-seconds** (sum over threads of the time inside the
bucket) and the pass as **wall seconds**. The ledger carries both. No optimisation
decision is taken from today's mixed timers.

## 5. Target architecture

### 5.1 The dependency DAG

The current `use` graph among the scoped modules has these edges that any merge must
respect (extracted from the sources at the baseline; arrows point at the imported module):

`application → {basis, delivery_3d, embed, embedding_io, fit_driver, merge, mstep, plane_cache, planes, project_gateway, rec3D, records, rounds, run_types, stages, state_service, targets, weights, probe_fit}`;
`fit_driver → {basis, embed, pairmerge, records, rounds, run_types, stages, util, probe_fit}`;
`probe_fit → {fit_types, records, rounds, run_types}` and its submodules → `{basis, crossfsc, plane_cache, planes, polar, posterior, mstep, pcg, latent_ops, artifacts, stages, util}`;
`fit_types → {crossfsc, mstep, polar, records, run_types}`; `posterior → fit_types`;
`basis → {artifacts, fit_types, mstep, pcg, plane_cache, planes, records, rounds, run_types, util, latent_ops}`;
`embed → {artifacts, basis, deconv, fit_types, plane_cache, planes, posterior, records, rounds, run_types, stages, util, latent_ops}`;
`pairmerge → {basis, crossfsc, mstep, pcg, records, util, probe_fit, latent_ops}`;
`mstep → {pcg, run_types, latent_ops}`; `pcg → latent_ops`; `polar → latent_ops`; `crossfsc → latent_ops`;
`planes → {plane_cache, latent_ops}`; `latent_ops → plane_cache`;
`rec3D → {rounds, run_types, stages, state_delivery, state_parts, states_backend, states_gridding, states_pcg}`;
`state_delivery → {project_gateway, states_backend, util}`; `states_* → {rounds, run_types, state_parts, states_backend, pcg | latent_ops}`;
`weights → {gmm, rec3D, records, rounds, run_types, targets, util}`; `gmm → {targets, util}`;
`merge → {pcg, rec3D, records, run_types}`; `state_service → {deconv, embedding_io, records, weights}`;
`project_gateway → {rounds, run_types, weights_state, weights_file}`; `delivery_3d → {plot, records, umap}`;
`run_types → records`; `rounds → stages`; `strategy → {application, artifacts, rounds, stages}`.

Target layering (each module imports only lower layers; contract modules stay small and
low, as the review asked):

| Layer | Modules (target names) | Role |
|---|---|---|
| L0 contracts | `simple_defs_flex` (in `src/defs`), `flex_pca_session` (today's records + the session of run_types), `flex_pca_rounds` (rounds + stages + artifacts: the distribution contract, the stage ids, part naming) | value records, identifiers, names |
| L1 leaves | `flex_plane_store` (cache + resident store, data only), `flex_latent_ops`, `flex_coupled_operator`, `flex_polar_bank`, `flex_pca_crossfsc`, `flex_pca_deconv`, `flex_pca_latent_clustering` (gmm + targets), `flex_pca_util`, `flex_pca_figures` (plot + umap), `flex_pca_io` (project_gateway + embedding_io + state_parts + weights_state: codecs and project reads/writes) | numerics and I/O with no upward imports |
| L2 fit | `flex_pca_fit_types` (payload), `flex_pca_mstep`, `flex_pca_posterior`, `flex_probe_fit` (+ submodules `_estep`, `_update` with the crossfsc submodule folded in, `_engine`), `flex_pca_basis`, `flex_pca_embed`, `flex_pca_pairmerge`, `flex_pca_fit_driver`, `flex_pose_refiner` (new, phase 7) | the probe fit |
| L3 states | `flex_pca_states` (states_backend + gridding + pcg + state_delivery + rec3D), `flex_pca_weights` (domain only), `flex_pca_merge` | state inference and reconstruction |
| L4 services | `flex_pca_delivery` (delivery_3d + state_service, which receives `cv_select_bandwidths`' orchestration) | application services |
| L5 | `flex_pca_application` | the application |
| L6 | strategy, commander | unchanged |

Ownership moves that make the layering true (phase 3, before any merge):

- `cv_select_bandwidths` splits: the bandwidth scoring (numerical) stays in `weights`; the
  loop that sets `params%outvol`, runs `reconstruct_flex_weighted_states`, reads and
  deletes the maps moves to the state service (L4), which already owns state placement.
  `weights → rec3D` disappears.
- `planes_batch_load` (prep orchestration calling `latent_ops`) moves out of the store into
  `latent_ops`, next to `prep_imgs4projected_model`; the store (L1) holds planes and the
  cache and imports nothing from flex. `latent_ops ↔ planes` cycle disappears.
- `fit_types` keeps only the payload and fit state records; its imports of `crossfsc`,
  `polar`, `mstep` become components typed through those L1/L2 modules in dependency
  order (crossfsc and polar are L1, mstep is L2 below fit_types in the file order).
- `state_delivery → project_gateway` is fine once gateway is in L1 `flex_pca_io`.
- `run_types` loses its settings half (phase 2); the session moves down to L0.

### 5.2 Contracts (phase 4)

- **Plane cache**: header = version, provenance (project file path digest and its
  modification stamp), geometry (box, box_crop, smpd), selection (count, hash, highest
  row written). Master validates every field; a worker adopts only when the header
  matches and its partition's rows are all below the highest row written. `kill` resets
  every component including the probe flags and the file name; a second session in the
  same process starts pristine. Tests in 4.2.
- **Session**: `flex_run_session` owns every product of the run and the four
  command-line facts; it is allocated by the application and released in its `kill`.
- **Payload**: `flex_probe_part`, one codec, version 1, borrowed from and restored to the
  fit without copies; the pose part of section 6.4 is a second payload with its own codec
  in the same module.
- **Artifacts**: every part file name is a function in `flex_pca_rounds` of (kind, fit,
  part, iteration); every part is published atomically (`.tmp` then rename); the master's
  completeness check lists what it expects before it reads.
- **Lifecycle**: every stateful type has `new`/`kill`, large components allocatable,
  no module-level variables in the scoped tree.

### 5.3 Consolidation map (phase 8, last)

Cycle-free by construction (every merge is within one layer and the merged module imports
only lower layers), cohesion-driven (one responsibility per file), 45 → 33 scoped files:

| Target | Today's | ~lines |
|---|---|---:|
| `simple_flex_pca_session` | records, run_types (session only) | 300 |
| `simple_flex_pca_rounds` | rounds, stages, artifacts | 240 |
| `simple_flex_plane_store` | planes (store), plane_cache | 400 |
| `simple_flex_latent_ops` | reconstructor_latent_ops + `planes_batch_load` | 1,050 |
| `simple_flex_coupled_operator` | pcg after the solver extraction | 1,900 |
| `simple_flex_polar_bank` | polar (as a type) | 580 |
| `simple_flex_pca_crossfsc` | crossfsc (as a type) | 450 |
| `simple_flex_pca_deconv` | deconv | 660 |
| `simple_flex_pca_latent_clustering` | gmm, targets | 1,340 |
| `simple_flex_pca_util` | util | 310 |
| `simple_flex_pca_figures` | plot, umap | 500 |
| `simple_flex_pca_io` | project_gateway, embedding_io, state_parts, weights_state | 1,400 |
| `simple_flex_pca_fit_types` | fit_types | 440 |
| `simple_flex_pca_mstep` | mstep | 340 |
| `simple_flex_pca_posterior` | posterior | 500 |
| `simple_flex_probe_fit` + `_estep`, `_update`, `_engine` | probe_fit, estep, update + crossfsc submodule, engine | 300 / 830 / 840 / 960 |
| `simple_flex_pca_basis` | basis | 1,210 |
| `simple_flex_pca_embed` | embed | 980 |
| `simple_flex_pca_pairmerge` | pairmerge | 910 |
| `simple_flex_pca_fit_driver` | fit_driver | 360 |
| `simple_flex_pose_refiner` | new (phase 7) | ~600 |
| `simple_flex_pca_states` | states_backend, states_gridding, states_pcg, state_delivery, rec3D | 1,180 |
| `simple_flex_pca_weights` | weights (domain only) | 350 |
| `simple_flex_pca_merge` | merge | 640 |
| `simple_flex_pca_delivery` | delivery_3d, state_service + the cv loop | 900 |
| `simple_flex_pca_application` | application | 540 |
| testers | `flex_pca_tester`, `flex_pcg_tester`, `flex_pca_application_tester` (new) | 350+ |

Not merged, on purpose: `embed` and `pairmerge` (different stages), `crossfsc` into `basis`
(`fit_types` imports crossfsc, so the merge would create `fit_types → basis → fit_types`),
`weights_state` into `weights` (store codec versus domain). The import cleanup (41 umbrella
imports → `only:` lists of what each module uses) is part of phase 8.

### 5.4 The PCG vector contract (phase 5)

`reconstructor_pcg` works on rank-3 volumes (`x(box,box,box)`), the coupled operator on
rank-4 component arrays; deferred bindings cannot vary dummy rank. The shared engine
therefore works on **rank-1 contiguous vectors**:

```fortran
type, abstract :: pcg_operator
contains
    procedure(op_size_i),    deferred :: size      ! n
    procedure(op_apply_i),   deferred :: apply     ! y(1:n) = H x(1:n)
    procedure(op_precond_i), deferred :: precond   ! z(1:n) = M^-1 r(1:n)
    procedure(op_dot_i),     deferred :: dot       ! <a,b> in double, the operator's own weighting
end type
```

with `real, contiguous, intent(in) :: x(:)` dummies. The engine (`simple_pcg_solver`:
residual replacement every `RESID_REPLACE`, the `xtol`/`rtol` stops, cold restart, the
outcome record) owns `r`, `p`, `z`, `Hp` as rank-1 of length `n` and never copies the
client's vector: the client passes a rank-1 pointer remap of its own storage
(`x_flat(1:n) => x`, legal for a contiguous target), and inside its `apply` the operator
remaps the rank-1 dummies back to its natural shape by the same pointer remap (zero copy).
`reconstructor_pcg` migrates first, A/B-checked by the `pcg_recon` and `lib_reconstruction`
gates; the coupled operator second. The engine's own tests (4.2) land before either
client moves.

## 6. `inpl_refine`: optional in-plane pose refinement in the E-step

Request (Hans): per particle and EM iteration, refine the in-plane rotation and shifts
before the latent solve, as the probabilistic polar refinement does: (1) shifts only
(shift-first), (2) exhaustive rotation at the found shift, (3) a joint `(sx, sy, rot_frac)`
polish (`inpl_cont`). `inpl_refine=(yes|no)`, default `no` = the current path. Decided: every
EM iteration; refined poses written back to the run's project copy; fastest-running
implementation. The E-step itself is not changed by the feature.

### 6.1 Where the latent and contrast come from (D3): two passes in the current iteration

`model%z` is zeroed at the top of every iteration (`engine:744`), distributed workers are
fresh processes with freshly constructed fits (`fit_driver`, one relaunch per iteration),
and no row-complete latent or contrast artifact exists. The plan therefore does not use
"the previous iteration's latent". Per particle, within the batch:

1. **Pass 1**: the E-step former at the stored pose gives `ẑ` and the contrast `a` under
   the current basis (the same computation the default path runs; its outputs are kept in
   the thread's scratch, not committed).
2. **Search**: the prediction `a·(U_0 + Σ ẑ_q U_q)` at the particle's bank direction, the
   three stages of 6.3, the accepted or rejected pose.
3. **Pass 2**: the E-step former at the refined pose; its outputs are what the batch rows,
   the posterior statistics and the M-step insertion use.

Cost: one extra former evaluation per particle per iteration (a gather and three BLAS
calls; section 1.2), no artifact, no cross-iteration state, and shared/distributed parity by
construction because everything a particle needs is in its batch. Mean-only (`ẑ=0`, skip
pass 1) is the same code behind a constant and is the fallback of 6.5.

### 6.2 The search bank

Flex's bank samples ring `r` at `nint(π r)` angles; `polarft_calc` samples every ring at
`pftsz`. Rotation by index shift exists only on the pftc grid, so a search-only bank is
built once per fit and iteration on the pftc grid: `ncomp+1` volumes × `ndir_search`
directions through the projector's `fproject_polar` (refine3D's reference preparation),
stored as a plain flex array `(pftsz, nk, 0:ncomp, ndir_search, eo)`. Memory at
`box_crop=64` (`pftsz≈98`, `nk≈30`, single-precision complex, even/odd): 0.5 GB at
`ndir_search=1000` for `ncomp=10`; the cap is a constant in `simple_defs_flex` (default
1000; the snap is the same approximation the E-step makes). The prediction is one `cgemv`
into a per-thread pftc reference slot.

### 6.3 The search context: exactly which pftc APIs, in which frames

- **Ownership.** `pftc_shsrch_grad` keeps a `class(builder), pointer` and evaluates
  `b_ptr%pftc`. Flex already passes `build` everywhere and `builder` carries a
  `polarft_calc` (`simple_builder.f90:28`). `flex_pose_refiner` therefore constructs
  `build%pftc` for the run (`nrefs = nthr` per-thread reference slots, `MAXIMGBATCHSZ`
  particle slots, `kfromto` = the fit band) and the optimizers take `build`; no interface
  refactor of the optimizer. With two fits in one process the slots are per thread, so one
  instance serves both.
- **Particles per batch**, the 3D matcher's contract (`simple_matcher_ptcl_batch`,
  `build_batch_particles3D`): the prepared image carries the **stored shift already
  applied** (`norm_noise_fft_clip_shift`), so every pftc solve is for a **shift increment
  in the cropped-pixel frame**, seeded at zero, never at the stored shift; `set_eo`;
  `create_polar_absctfmats(build%spproj, 'ptcl3D')` for the batch (the Euclidean objective
  needs the polar CTF matrices); sigma2 from flex's per-particle noise model through
  `assign_sigma2_noise`; `memoize_ptcls` once per batch, outside the OpenMP region.
- **Thread scratch.** The refiner owns, per thread: the pftc reference slot, a
  `pftc_shsrch_grad` built with `new_fixed` (stage 1), one built with `new_joint` (stage
  3), and the prediction buffer. Evaluation inside the particle loop touches only the
  thread's slot and the particle's slot, as the matcher does.
- **Stage 1, shifts only**: `new_fixed` at the stored discrete angle (the alternating
  constructor keeps updating the angle whenever `opt_angle` is set, so it is not used
  here); `minimize` from the zero increment within `±trs`.
- **Stage 2, exhaustive rotation at the found increment**: the FFT path
  (`gen_objfun_vals`) needs reference memos, and `memoize_refs` refuses to run inside an
  OpenMP region; a per-thread reference set inside the particle loop cannot be memoized
  as pftc stands. Stage 2 therefore evaluates the objective **directly per rotation**
  (`gen_euclid_for_rot_8`/`gen_corr_for_rot_8` at each of the `nrots` rotations, the
  explicit-rotation path that `calc_frc` uses): `nrots × pftsz × nk` complex
  multiply-adds per particle, about 0.6 M at `box_crop=64`, well under the former's own
  cost. D2 reopens the FFT path only if measurement says so.
- **Stage 3, joint polish**: `new_joint` with the angular window `±2` cells around the
  stage 2 winner and the shift box `±trs`, `minimize_joint`, both acceptance guards as in
  refine3D (material improvement over the seed; a bound-pinned solution is non-improving).
- **Frames and signs.** Rotations found by pftc are relative to the bank direction's
  frame; the particle's `e3` is recovered by the inverse of `polar_relative_inplane`.
  Shift increments are in cropped pixels in the particle's frame (not the reference
  frame; refine3D rotates stored shifts into the native frame before its joint solve,
  `strategy3D_srch:367`); they are converted to native pixels by the crop factor and added
  to the stored shift, with refine3D's sign convention, when the pose is written. The
  prepared plane used by pass 2 and by the insertion receives the increment as a Fourier
  phase ramp.
- **Acceptance.** An accepted solve commits `(e3, x, y)` to `orientations(i)`; a
  non-improving or invalid solve leaves the stored pose; the batch log reports the refined
  fraction and the mean `|Δe3|`, `|Δshift|` from discrete cells (statistics parity).

### 6.4 Distribution: the pose part

Workers publish per part. One pose part per (fit, part, iteration):
`flex_pca_poses_fit<f>_part<p>_iter<it>.bin`, atomic publication, header (magic, version,
fit, part, iteration, row count), rows `(project row, e3, x, y, accepted flag, score)` in
the worker's partition order. The master, before merging: expects exactly the parts of its
partition plan, refuses a missing part, refuses duplicate rows across parts and rows
outside the plan, then writes the accepted poses into its project copy's `ptcl3D` and
rewrites the project before launching the next round (the same write it already does for
`ptcl3D/state` after the state stage; phase 7(e) verifies it is in place for the probe
rounds). Shared-memory runs commit directly. Parity between the two is a workflow-gate
metric (4.4).

### 6.5 Fallbacks and the statistical note

Mean-only prediction (`ẑ=0`) if the search bank's memory is prohibitive at large
`box_crop`. The pose is refined against the plug-in prediction and the posterior
recomputed at the refined pose (coordinate ascent on the joint likelihood over pose and
z), as refine3D's `inpl_cont` does for (pose | reference); not the marginal likelihood. An
early lock to the wrong conformer is limited by the acceptance guards and by pass 1's
latent being the current iteration's (not a stale one); the jittered-phantom arm's "truth
metrics not lowered" assertion is the guard in the test net.

**Parameters.** `inpl_refine` (`yes|no`, default `no`) through the normal lifecycle
(`simple_parameters.f90`, `_parse.f90`, the heterogeneity UI, the commander); `trs`,
`maxits_sh` reused; `objfun` fixed to euclid inside; validation refuses `inpl_refine=yes`
without the sigma2 input.

## 7. Optimisation policy

No phase changes loop order, precision or arithmetic for speed on its own initiative.
After phase 0's instrumented measurement (4.5) on the phantom and, in a worktree, at a
realistic box and particle count, a conditional phase (5b) may address **one** bucket
that dominates, with a fast-tier test pinning the new kernel to the old at round-off and
the thread-seconds before and after in the ledger; a step that does not pay is reverted.
Otherwise 5b is skipped and the ledger says so.

## 8. Non-goals

- No change to the estimators, defaults or scientific behaviour of the default path.
- No change to the E-step: not its hybrid split at `0.72·band`, not its exact low-`k`
  part, not its polar grid or noise identity.
- No new abstractions beyond the solver interface, the bank type and the refiner type.
- No reader for old artifacts; no migration tool.
- No marginal-likelihood pose objective; no projection-direction search.
- No move of `cls_expansion`. No GPU or device path.

## 9. Phases (the review's order)

Each phase is one or two commits by the maintainer, lands only with every tier green, and
is recorded in the ledger.

| # | Phase | Scope | Exit |
|---|---|---|---|
| 0 | Test net on `master` | Instrumentation (4.5); 4.2 oracles that exist today, 4.3 tester, 4.4 gate and its registrations (D1); observations and oracle-derived tolerances recorded separately | all tiers green on `master`; ledger row 0 |
| 1 | Proven-dead code only | The section 1.1 dead list (pure-Cartesian former, pinned overrides, ECM, dead procedures), the second codec, stage ids 1–4, duplicate constants. Kept untouched: the hybrid E-step and all its fields | tiers green; A/B 1e-5 |
| 2 | Configuration into `parameters` | Environment reader and all non-default branches gone; `flex_run_settings` removed; leaves stop taking `parameters`/`builder` | tiers green; A/B 1e-5; `rg SIMPLE_COV_ src` empty |
| 3 | Dependency DAG and ownership moves (D4) | The 5.1 moves; `simple_defs_flex`; module-level `use` lists checked against the DAG by a script that fails on an upward edge | tiers green; A/B 1e-5; the script passes |
| 4 | Contracts and objects | 5.2: plane cache (with its tests), session, payload, artifacts, lifecycle; `flex_plane_store`, `flex_polar_bank`, `flex_pca_crossfsc`, the weights store as types; no module-level variables | tiers green; A/B 1e-5; cache tests pass |
| 5 | PCG vector contract | 5.4: engine + its tests, then `reconstructor_pcg`, then the coupled operator | tiers and reconstruction gates green at their floors; before/after metrics recorded |
| 5b | Conditional optimisation | Section 7 | thread-seconds before/after; or "skipped" |
| 6 | `inpl_refine` design closure | D2/D3 confirmed; the pftc prerequisites checked in a tester (`new_fixed` holds the angle; direct per-rotation evaluation matches `gen_objfun_vals` on a memoized reference outside OpenMP; `create_polar_absctfmats` contract); the frame conversions pinned by tests before any flex code | tests pass; section 6 updated to what was verified |
| 7 | `inpl_refine` implementation | (a) parameter plumbing, default no, A/B identical; (b) `flex_pose_refiner`: search bank, pftc instance, thread scratch; (c) two-pass E-step, three stages, guards, frames; (d) shared-memory write-back and the log line; (e) pose parts and the master merge (6.4); (f) the 4.3/4.4 arms | all tiers green incl. new arms; pose RMS, refined fraction, bank GB, search thread-seconds in the ledger |
| 8 | Physical consolidation and import cleanup | 5.3 map; umbrella imports → `only:`; README table; CMake re-configure note | tiers green; A/B 1e-5; `ls` matches the table; DAG script passes |
| 9 | Style | Indentation scars, headers, comments, line length, `code_base_map.md` | tiers green; `git diff -w --ignore-blank-lines` empty for code |
| 10 | Smaller reuse | Cholesky family into `simple_linalg`; cross-fit FSC on `image%fsc` where shells coincide; k-means on `utils/clustering` if the inputs match | tiers green; A/B 1e-5 |

## 10. Risks

- The E-step former's oracle (phase 0) is the only independent check of the E-step; it
  lands before phase 1 touches anything near it.
- Phase 2: a branch the default path reaches under some input (`merge` with
  `preimage_auto=yes`) is default behaviour and stays; phase 2 lists, per key, what is
  kept.
- Phase 3's DAG script is what keeps phase 8 honest; without it a merge can reintroduce an
  upward edge that compiles only by file order.
- Phase 5 moves numbers; floors, not identity, decide; a floor is never loosened without
  a written reason.
- Phase 7: the search bank's memory (6.2, capped); the pftc behaviours relied on are
  checked by phase 6's tester before flex code depends on them; the two-pass E-step
  doubles one bucket, measured in the ledger.
- Tolerance creep: a floor widened to pass is a finding (policy 2.1).

## 11. Ledger

| Phase | Commit | Date | Tiers | A/B latents | Buckets (thread-s / wall-s per 4.5) | Notes |
|---:|---|---|---|---|---|---|
| 0 | | | | n/a | master: __ ; worktree realistic: __ | observations: assignment acc __, eigvol corr __, state FSC0.143 __, shmem/distr latent corr __; oracle tolerances: former __, M-step toy __ |
| 1 | | | | | | lines removed __ |
| 2 | | | | | | keys removed 20; branches kept: __ |
| 3 | | | | | | DAG script: __ |
| 4 | | | | | | cache tests __ |
| 5 | | | | n/a | | before/after metrics __ |
| 5b | | | | n/a | | run or skipped |
| 6 | | | | n/a | | pftc prerequisites verified: __ |
| 7a–7f | | | | default path | search __ | pose RMS __ (jitter __), refined fraction __, bank GB __ |
| 8 | | | | | | files 45 → __ |
| 9 | | | | | | |
| 10 | | | | | | |

## 12. What the review of 2026-10-06 corrected

1. Baseline and scope: 45 files / 21,422 lines scoped; the class-expansion file excluded
   from every count; the consolidation target restated as 45 → 33 cohesive, cycle-free
   modules instead of 46 → 23.
2. The hybrid E-step was misclassified as unreachable; it is the production E-step and is
   kept untouched. Only proven-dead code is removed.
3. The first consolidation map created `application ↔ probe_fit`, `states ↔ delivery`,
   `basis ↔ fit_types` and `io ↔ weights` cycles; the plan now draws the DAG, moves
   ownership first (5.1) and merges last (5.3).
4. The PCG abstraction had no viable vector contract; 5.4 specifies the rank-1 contiguous
   contract with pointer remaps and tests the engine before migrating either client.
5. `inpl_refine` assumed a previous-iteration latent that does not exist (reset per
   iteration, fresh worker fits, no artifact); 6.1 chooses the two-pass design. The pftc
   assumptions were wrong (stored shift already applied in the batch prep; `sh_rot` is not
   shift-only; the optimizer evaluates `builder%pftc`; memoization is forbidden inside
   OpenMP; CTF matrices are part of the batch contract); 6.3 specifies the search context
   against the actual APIs. 6.4 specifies the per-part pose artifact.
6. Plane-cache correctness (selection-blind adoption, stale lifecycle) is now a contract
   with tests, not an ownership move.
7. Test design: observation, oracle and regression margin separated (4.1); the invalid
   90° control withdrawn; the M-step toy shrunk to a sub-second oracle; PSD instead of PD;
   signed permutations for the cross-fit test; the in-process worker arm labelled as
   partition-reduction coverage only.
8. Registration: the `COMMAND` argument, the CTest budget decision (D1), the next-to-code
   application tester and the synchronised suite/UI registration.
9. `cv_select_bandwidths` orchestration leaves the weights domain module; the umbrella
   import cleanup is a phase item.
10. Instrumentation: thread-seconds and wall seconds defined per bucket before any timing
    drives a decision.
