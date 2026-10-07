# FLEX refactoring plan: reduction, reuse, and a numerical test net

Date: 2026-10-06. Revision 39. Status: the maintainer reports the phase-0 test net,
shared/distributed workflow gate, phase 1, phase 2 and phase 3a green. The log values for the same-frame
map/eigenvolume, workflow-parity and timing observations have not been supplied, so their ledger
cells remain open rather than being reconstructed. Phase 3b is source-complete: bandwidth-CV
reconstruction/I/O orchestration belongs to the state service while the weights module retains the
grid, scoring and adoption numerics; shared cross-layer constants live in `simple_defs_flex`; and
the fast gate runs a FLEX layer/cycle checker. Phase 4a is source-complete: the plane cache now
validates project provenance, geometry, exact master selection and worker coverage; its lifecycle
resets all session state; and the fast FLEX suite contains the required contract tests. Phase 4b
is source-complete: one application-owned `flex_plane_store` now owns the
resident-plane data and its allocatable disk-cache object; the former module-level resident
and cache state is gone, and the store is passed explicitly through the mean, fit, E-step and
embedding paths. Phase 4c is source-complete: the application owns an allocatable session;
fits own allocatable polar banks and PCG environments; cross-FSC records, files and contexts
have explicit lifecycles; rounds own the artifact catalog passed to every codec; and the
weights-store type owns validated deliveries. The scoped production tree retains no mutable
run state at module scope. Static checks are green; the maintainer compile and tier rerun are
pending.

Phase 5 is source-complete: `simple_pcg_solver` owns the rank-1 Krylov recurrence,
outcome and explicit policy options; reconstruction and the coupled FLEX operator expose
zero-copy rank remaps through small adapters while retaining their own operators,
preconditioners, windowing and retry policy. The optimiser unit suite now contains the
dense-SPD, exact-diagonal and cold-restart checks. Phase 5b is skipped because no retained
before/after thread-second measurements exist; the optimisation policy forbids speculative
changes without that evidence. A maintainer compile corrected the reconstruction adapter's
explicit-shape target declaration: it is inherently contiguous and cannot carry the Fortran
`CONTIGUOUS` attribute. Phases 6--8 are source-complete. Every module and submodule remains in
its own source file; every core umbrella import in production FLEX code is now an explicit
`only:` list, and the DAG checker enforces both the file boundary and the dependency layers.
A maintainer compile exposed seven former transitive OpenMP
runtime imports; their owning module units now import the exact `omp_lib` procedures directly.
The next compile exposed duplicate use association in the four probe-fit submodules; names supplied
by the parent are now taken only by host association, while submodule-local dependencies retain
narrow imports. The pair-merge deflation helper's fit dummy is `intent(inout)`, matching its calls
through the fit-owned envelope. The FLEX commander's narrowed API import explicitly includes the
`simple_abspath` dependency used while resolving a project-owned consensus volume.
The style pass updated physical-file descriptions, this README
and the generated source maps without reformatting procedure bodies. The exact duplicated
real-SPD/complex-RHS Cholesky solve now lives in `simple_linalg`; cross-FSC already uses
`image%fsc`, and the weighted deterministic latent k-means has no contract-compatible shared
replacement. Static checks are green; compilation and runtime tiers remain pending under the
repository policy.

This plan is deliberately limited to refactoring and its numerical safety net. New scientific
features belong in separate design work and are not part of any phase, contract, test arm, module
target, or artifact change below.

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

Validation level: static source inspection (`rg`/`awk`/python over the tree and the extracted
module `use` graph), plus the maintainer's targeted compile and test runs. Verified measurements
and outstanding observations live only in the ledger. Timing statements remain structural until
phase 0 measures them with the instrumentation of section 4.5.

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
  `flex_weights_commit`, `_merge_local_ranges`, `_prepare_update`, `_range_path`; the unfused
  `polar_sample_particle`, `polar_sample_at_pose`,
  `polar_dir_neighbours`, `polar_apply_shift` and their unused noise-ring/half-label grid storage;
  `invert_lower`, `punit`;
  `flex_pcg_t%get_cdim`, `%get_pijk`, `%list_signature`. Removing them changes no number.
- **Configuration outside `parameters`**: `flex_run_settings` reads 20 `SIMPLE_COV_*`
  environment keys. Every behaviour reachable only through a non-default key value is
  unsupported after this round.
- **The baseline had two live codecs for one payload**: separate dense single-fit and packed
  paired formats both serialised `flex_probe_part`; stage ids 4–7 kept 1–3 unused. Phase 1b
  replaces them with one version-13 packed codec for one- or two-fit containers and stage ids
  1–4; old part files are deliberately rejected.
- **The baseline duplicated constants**: `MERGED_PC_FBODY`, `MERGED_META`,
  `MERGED_EIG_FNAME`, `PAIRED_MANIFEST`; the FLEX relative PCG lambda; the 0.143 FSC
  signal threshold; and three `reconstructor_pcg` stopping/preconditioner constants.
  Phase 1b gives each set one owning module without changing its value.
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
   distributed-master and worker execution modes.
2. The delivered products: embedding cache, eigenvolumes, latent readouts and figures,
   state weights (`flex_weights_state_NNN.bin`), state maps (`vol_flex`), the
   `ptcl3D/state` labels and the out-segment registrations.
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
2. Artifact formats: part files, caches, and the weights store may
   change layout and version freely; no reader accepts an older layout (none does today).
   A run restarts on the new binary; resume reads only caches the same binary wrote.
3. The `SIMPLE_COV_*` environment keys and every behaviour reachable only through one.
4. The unreachable pure-Cartesian former and the pinned split/oversampling overrides.

Maintainer rules apply throughout (`CLAUDE.md`, `AGENTS.md`): edits by read-modify-write,
no compilation by the agent, the maintainer builds and runs the gate after each phase,
commits are the maintainer's.

## 3. Owner rulings

| # | Decision | Needed by | Ruling |
|---|---|---|---|
| D1 | CTest process budget for the additional `flex_pca_blobs` workflow gate. | phase 0 | Approved: add the gate and raise `SIMPLE_CTEST_BUDGET` from 32 to 33 in the same registration change. |
| D2 | The target dependency DAG and file-boundary rule of sections 5.1–5.3. | phase 3 | Dependency map accepted. Every module and submodule stays in its own file; physical consolidation is not permitted. With the application tester, the scope contains 46 source files. |

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
  against a brute-force sum over the identical polar rows written in the test, with a
  separately labelled polar-quadrature versus Cartesian-lattice observation. The deterministic
  fixture is a sum of off-centre narrow blobs with variances 0.8--1.1, and it asserts that the
  mean and both basis diagonals retain resolved power in the polar-ring region (a direct DFT of
  the central sections puts 16--32% of their lattice power between the hybrid gate and `BAND`).
  A fixture-design floor requires every tested diagonal to retain at least 10% of its whole-band
  power on those rings; this protects the energetic fixture rather than merely excluding round-off.
  The ring part is
  isolated as the former result minus the Cartesian low-`k` oracle; `G`, `b`, `c`, `e_mm`, and
  `myv` are each normalized by the matching ring-only oracle rather than by a low-`k` aggregate.
  The phase-0 fixture gathers from `cmat_exp` with an independently written normalized KB loop;
  it guards the ring algebra, indexing, half-plane and Friedel paths, but deliberately does
  not claim an independent check of expansion or the interpolation model. The
  low-`k` contribution alone equals the oracle restricted to `h²+k² ≤ rhyb(rhyb+1)` at
  round-off. Linearity oracle: the Cartesian sufficient statistics of a noise-free particle built
  as `μ + Σ c_q u_q` must give `b = G c` and `myv = e_mm + cᵀc₀` at the accumulation bound.
  The production hybrid uses a Cartesian plane resampled to polar points for the particle but
  direct 3-D volume samples for the bank, so its corresponding identity defect is an interpolation-
  path observation rather than an exact algebraic assertion; its measured values live in the
  ledger. (The earlier "90° rotation gives an orthogonal `b`" control is withdrawn: not generally
  true.)
- Coupled M-step on a compact toy: two supported 8³ blob volumes, 40 projections at fixed-seed
  poses and deterministic latents with known `E[zz']=zz'`. Gridding and PCG must retain both
  truth directions (minimum principal-angle cosine 0.90, a rank-collapse guard on this coarse
  three-tap/40-view fixture). The PCG hard support leaves 64 coupled real unknowns; the test
  assembles the finalized operator by unit-vector probes and solves that matrix independently
  by double-precision Cholesky. A before/after operator probe asserts that the zero-vector
  preconditioner preparation installs a positive diagonal Tikhonov lift. The relative-residual
  bound is twice the single-precision accumulation estimate, leaving order/compiler margin while
  remaining a strict backward-error check. The independently solved dense system remains a reported
  forward-error observation: its condition-amplified bound is too loose to be a meaningful
  assertion on this fixture. The separate PCG
  operator fixture remains the independent check of the probed operator.
  The 16³ recovery version moves to the library tier.
- PCG operator: symmetry `<x,Hy> = <Hx,y>` and `<x,Hx> ≥ 0` on two deterministic formula
  vectors (the masked unregularised operator is positive semidefinite, not definite); positive
  definiteness asserted only after `prep_floor` derives the operator's Tikhonov term from the
  matching gridding density. Kept as the one white-box self-test.
- Cross-fit FSC: two half bases related by a signed permutation (what the matcher supports,
  `crossfsc.f90:36`) give the series of the identity; a random basis gives the noise floor.
- Probe-part codec: write/reduce/read round trip of a random `flex_probe_part` is exact.
- Plane cache (with phase 4): changed selection → rebuild; changed project → rebuild;
  worker adoption of a cache whose selection does not cover the partition → refused; two
  sessions in one process in either order see a pristine store.
- PCG solver engine (with phase 5): on a dense SPD operator against a direct solve; on a
  diagonal operator convergence in one step; stop reasons and outcome fields pinned.
- Kept: embedding cache I/O, auto settings, population floor, kernel bandwidth, state
  weights, deconvolution.

### 4.3 Library tier: `lib_heterogeneity` (in-process end-to-end on a phantom)

A new tester next to the application it tests, `run/simple_flex_pca_application_tester.f90`, one
sub-suite `flex PCA two-state phantom`, hermetic, deterministic, minutes, registered in
`suites_lib_heterogeneity` and in the area's `suite=` help:

- Phantom: an asymmetric sum of Gaussian blobs (`add_gaussian_blob`), box 64, conformers A
  and B differing by one blob moved 8 pixels; 2,000 particles at fixed-seed uniform
  orientations, CTF and noise at a declared SNR, labels alternating by row. Truth: the maps,
  the labels, A−B, the poses.
- Builds one project-backed simulation and runs `commander_flex_pca%execute` in shared memory for
  the four existing combinations of `rec_backend={gridding,pcg}` and
  `rec_states_backend={gridding,pcg}`. The project is part of the fixture because publication of
  the weight store, maps and hard labels is part of the application contract. Calling the public
  commander covers its defaults, canonical-sigma preflight, selection setup and strategy lifecycle
  rather than duplicating those steps in the fixture. The native and covariance boxes are both 64,
  so the known plane-cache lifecycle fault cannot couple successive arms before phase 4 fixes it.
  The fixture explicitly sets `min_neff=100`: the production default is 2,000, which would prune
  every split of a 2,000-particle phantom rather than test recovery.
- Asserts against **simulation truth** with arithmetic-derived tolerances where the metric
  has them: label assignment accuracy after optimal state-to-truth recoding; leading eigenvolume
  correlation with A−B (sign-free); two distinct delivered maps optimally matched to the two
  truth maps in their shared simulation frame (`compare_to_truth` directly, without docking),
  with each matched map required to prefer its own conformer; weight store invariants,
  maximum-weight hard labels and project-label agreement.
  Where the metric is an end-to-end number without a closed form, the floor is the phase 0
  observation minus a stated regression margin, labelled as such.
- The initial run asserts the independent checks available before calibration: publication retains
  at least two states and no more than the explicit provision ceiling of three; optimal binary
  recoding of the delivered state ids clears the 0.60 chance-separation floor; and the published
  stores agree. With at most three ids there are only `2^3` recodings, so Hoeffding's union bound
  under random truth labels gives `P(best accuracy >= 0.60) < 8 exp(-2·2000·0.10²) < 4e-17`.
  The per-state truth counts and full state-by-truth correlation matrix are logged for each arm so
  every ledger observation can be regenerated from test output. Eigenvolume and same-frame map
  correlations are phase-0 observations; regression floors are added only after the observations
  and margins are recorded.
- No one-state application arm is added: the current publication contract rejects fewer than two
  states, so requiring a collapse to one would first require a production feature change. No
  synthetic in-process `run_worker` arm is added either: `run_worker` is a stage/artifact endpoint,
  not an in-memory partition interface. Real partition reduction, serialization, project reload
  and shared-versus-distributed agreement remain the workflow gate's responsibility.
- Records the E-step split (section 4.5) for the ledger.

### 4.4 Workflow gate: `flex_pca_blobs` (highlevel, nightly)

A commander in `simple_commanders_test_highlevel_*` (policy section 4.4), with its case in
the router, its program in the test UI, and this CMake entry:
`simple_add_test(flex_pca_blobs highlevel 3600 8 RUN_SERIAL COMMAND $<TARGET_FILE:simple_test_exec> test=flex_pca_blobs nthr=4)`.
The entry reserves eight processors because its distributed arm starts two four-thread workers.
Under ruling D1 the same registration raises `SIMPLE_CTEST_BUDGET` to 33. The gate reuses the
library fixture's exact truth, particle stack, project recipe and application command line, then
runs the public `flex_pca` commander in shared memory and **distributed (`nparts=2`, real local
worker processes)** PCG/PCG modes. Each arm independently asserts the truth-label and same-frame
map-ordering oracles and publication invariants. It writes `metrics.tsv` through
`simple_test_gate` and fails through the gate. The first run reports, without premature floors,
the sign/permutation-invariant latent correlation, truth-recoded label agreement, delivered state-
count difference and corresponding truth-matched map correlations. This is where distributed
behaviour is covered: part files, project reload between rounds, master merge, and
shared-versus-distributed agreement.

### 4.5 Instrumentation (phase 0, before any number is compared)

One definition per bucket, applied to both the per-pass line and the per-thread sums:
`read`, `prep`, `bank`, `exact_lowk`, `ring`, `solve`, and `insert`, each reported as
**thread-seconds** (sum over threads of the time inside the
bucket) and the pass as **wall seconds**. The ledger carries both. No optimisation
decision is taken from today's mixed timers.

## 5. Target architecture

### 5.1 The dependency DAG

The current `use` graph among the scoped modules has these edges that any ownership move must
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

Target layering (each module imports only lower or same-layer modules, with the
same-layer graph acyclic; these are logical groupings, and every module or
submodule remains in its own source file):

| Layer | Logical grouping | Role |
|---|---|---|
| L0 contracts | `simple_defs_flex` (in `src/defs`), `simple_flex_pca_records`, `simple_flex_pca_run_types`, `simple_flex_pca_stages`, `simple_flex_pca_artifacts`, `simple_flex_pca_rounds` | value records, identifiers, names |
| L1 leaves | `simple_flex_pca_plane_cache`, `simple_flex_pca_planes`, `simple_flex_reconstructor_latent_ops`, `simple_flex_pca_pcg`, `simple_flex_pca_polar`, `simple_flex_pca_crossfsc`, `simple_flex_pca_deconv`, `simple_flex_pca_gmm`, `simple_flex_pca_targets`, `simple_flex_pca_util`, `simple_flex_pca_plot`, `simple_umap`, `simple_flex_pca_project_gateway`, `simple_flex_pca_embedding_io`, `simple_flex_pca_state_parts`, `simple_flex_weights_state`, `simple_flex_weights_file` | numerics and I/O with no upward imports |
| L2 fit | `simple_flex_pca_fit_types`, `simple_flex_pca_mstep`, `simple_flex_pca_posterior`, `simple_flex_probe_fit` and its four submodules, `simple_flex_pca_basis`, `simple_flex_pca_embed`, `simple_flex_pca_pairmerge`, `simple_flex_pca_fit_driver` | the probe fit |
| L3 states | `simple_flex_pca_states_backend`, `simple_flex_pca_states_gridding`, `simple_flex_pca_states_pcg`, `simple_flex_pca_state_delivery`, `simple_flex_pca_rec3D`, `simple_flex_pca_weights`, `simple_flex_pca_merge` | state inference and reconstruction |
| L4 services | `simple_flex_pca_delivery_3d`, `simple_flex_pca_state_service` (which receives `cv_select_bandwidths`' orchestration) | application services |
| L5 | `flex_pca_application` | the application |
| L6 | strategy, commander | unchanged |

Ownership moves that make the layering true (phase 3, before import cleanup):

- `cv_select_bandwidths` splits: the bandwidth scoring (numerical) stays in `weights`; the
  loop that sets `params%outvol`, runs `reconstruct_flex_weighted_states`, reads and
  deletes the maps moves to the state service (L4), which already owns state placement.
  `weights → rec3D` disappears.
- `planes_batch_load` (prep orchestration calling `latent_ops`) moves out of the resident
  store into `latent_ops`, next to `prep_imgs4projected_model`; in phase 3 the resident
  store and disk cache are separate L1 leaves. In phase 4 the application-owned
  `flex_plane_store` composes the cache through the acyclic same-layer `planes →
  plane_cache` edge, and `latent_ops` receives that one object explicitly.
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
- **Plane store**: one allocatable `flex_plane_store` component belongs to each application
  process. It owns the resident row store and one allocatable `flex_pca_plane_cache`, is
  passed explicitly to every particle-read path, and releases both in application `kill`.
  Neither source module retains run state at module scope.
- **Session**: `flex_run_session` owns every product of the run and the four
  command-line facts; it is allocated by the application and released in its `kill`.
- **Payload**: `flex_probe_part`, one codec, version 13, borrowed from and restored to the
  fit without copies.
- **Artifacts**: every part file name is a function in `flex_pca_rounds` of (kind, fit,
  part, iteration); every part is published atomically (`.tmp` then rename); the master's
  completeness check lists what it expects before it reads.
- **Lifecycle**: every stateful type has `new`/`kill`, large components allocatable,
  no module-level variables in the scoped tree.

### 5.3 File boundary and import cleanup (phase 6, last)

Every module and submodule has its own `.f90` source. The 45-file baseline gains only the planned
application tester, leaving 46 scoped source files. `scripts/check_flex_dag.py` rejects a source
that declares more than one module or submodule, in addition to rejecting upward edges, cycles,
unmapped modules and broad core imports. The repository-wide `scripts/check_descr.py` enforces
the same one-unit boundary for all SIMPLE Fortran sources, and the fast gate runs both checks.

Phase 6 keeps the established file boundaries for records/run session, stages/artifacts/rounds,
cache/plane store, update/cross-FSC submodules, targets/GMM, plotting/UMAP, persistence modules,
state reconstruction roles and delivery/state service. Their responsibilities remain separate;
the dependency cleanup is expressed through ownership and explicit interfaces, not physical
co-location. Production FLEX imports of `simple_core_module_api` all use explicit `only:` lists;
the strategy boundary is subject to the same gate.

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

## 6. Optimisation policy

No phase changes loop order, precision or arithmetic for speed on its own initiative.
After phase 0's instrumented measurement (4.5) on the phantom and, in a worktree, at a
realistic box and particle count, a conditional phase (5b) may address **one** bucket
that dominates, with a fast-tier test pinning the new kernel to the old at round-off and
the thread-seconds before and after in the ledger; a step that does not pay is reverted.
Otherwise 5b is skipped and the ledger says so.

## 7. Non-goals

- No change to the estimators, defaults or scientific behaviour of the default path.
- No change to the active E-step arithmetic: not its hybrid split at `0.72·band`, its exact
  low-`k` part, or the polar grid's geometry, ordering and weights.
- No new abstractions beyond the solver interface and the bank type.
- No reader for old artifacts; no migration tool.
- No new scientific features.
- No move of `cls_expansion`. No GPU or device path.

## 8. Phases (the review's order)

Each phase is one or two commits by the maintainer, lands only with every tier green, and
is recorded in the ledger.

| # | Phase | Scope | Exit |
|---|---|---|---|
| 0 | Test net on `master` | Instrumentation (4.5); 4.2 oracles that exist today, 4.3 tester, 4.4 gate and its registrations (D1); observations and oracle-derived tolerances recorded separately | all tiers green on `master`; ledger row 0 |
| 1 | Proven-dead code only | The section 1.1 dead list (pure-Cartesian former, pinned overrides, ECM, dead procedures), the second codec, stage ids 1–4, duplicate constants. Kept untouched: the active hybrid E-step arithmetic and data | tiers green; A/B 1e-5 |
| 2 | Configuration into `parameters` | Environment reader and all non-default branches gone; `flex_run_settings` removed; leaves stop taking `parameters`/`builder` | tiers green; A/B 1e-5; `rg SIMPLE_COV_ src` empty |
| 3 | Dependency DAG and ownership moves (D2) | The 5.1 moves; `simple_defs_flex`; module-level `use` lists checked against the DAG by a script that fails on an upward edge | tiers green; A/B 1e-5; the script passes |
| 4 | Contracts and objects | 5.2: plane cache (with its tests), session, payload, artifacts, lifecycle; `flex_plane_store`, `flex_polar_bank`, `flex_pca_crossfsc`, the weights store as types; no module-level variables | tiers green; A/B 1e-5; cache tests pass |
| 5 | PCG vector contract | 5.4: engine + its tests, then `reconstructor_pcg`, then the coupled operator | tiers and reconstruction gates green at their floors; before/after metrics recorded |
| 5b | Conditional optimisation | Section 6 | thread-seconds before/after; or "skipped" |
| 6 | File-boundary enforcement and import cleanup | 5.3 rule; umbrella imports → `only:`; README layout; CMake re-configure note | 46 scoped files; one module/submodule per file; no broad core import; DAG script passes; runtime tiers pending |
| 7 | Style | Indentation scars, headers, comments, line length, `code_base_map.md` | tiers green; `git diff -w --ignore-blank-lines` empty for code |
| 8 | Smaller reuse | Exact duplicated Cholesky solve into `simple_linalg`; use existing `image%fsc`; reuse clustering only if its weighted deterministic contract matches | shared solve + direct unit oracle; mismatched candidates documented; runtime tiers pending |

## 9. Risks

- The E-step former's oracle (phase 0) is the only independent check of the E-step; it
  lands before phase 1 touches anything near it.
- Phase 2: a branch the default path reaches under some input (`merge` with
  `preimage_auto=yes`) is default behaviour and stays; phase 2 lists, per key, what is
  kept.
- Phase 4 closes the former `flex_env_init` runtime precondition: the application constructs one
  canonical `flex_pcg_environment` per process, resident fits own exact copies of its already-
  resampled support, and the state merge borrows the application object. A future caller must
  obtain an initialized environment object explicitly.
- Phase 3's DAG script keeps phase 6 honest: it rejects upward edges and any attempt to put
  multiple module units in one source file.
- Phase 5 moves numbers; floors, not identity, decide; a floor is never loosened without
  a written reason.
- The energetic phase-0 fixture exposes a production interpolation-path mismatch: the particle
  reaches polar points by Cartesian-plane resampling while the bank samples the 3-D volumes
  directly. The resulting linearity defect can bias latent estimates even for a noise-free
  particle. This refactoring does not change the E-step; follow-up scientific work must measure
  the band-edge bias on real data and decide whether both paths should use one interpolation.
- The two-gate merge retained an intermediate representative on a discrete two-mode phantom, but
  it held only 5 of 2,000 hard labels: a nearly empty extra state, not a substantial split of one
  conformation. Label separation remains strong, but a downstream multi-state 3-D refinement could
  still start from a redundant near-empty class. This refactoring does not change merge policy;
  follow-up scientific work must test whether that state persists on realistic data and whether
  the merge needs a stronger distinctness criterion.
- Tolerance creep: a floor widened to pass is a finding (policy 2.1).

## 10. Ledger

| Phase | Commit | Date | Tiers | A/B latents | Buckets (thread-s / wall-s per 4.5) | Notes |
|---:|---|---|---|---|---|---|
| 0 | working tree | 2026-10-06 | maintainer reports all current tiers and the shared/distributed workflow gate green; application phantom completed all four arms in 85.14 s and retained three provisioned states per arm | green; numeric parity observation not supplied | master: __ ; worktree realistic: __ | hybrid same-sample 3.140e-7, ring fraction 1.887e-1 (fixture floor 0.10), ring error 1.341e-6, low-k 3.975e-7, polar linearity 1.911e-2 (~10% of ring contribution), polar/lattice 3.572e-2 (~19%); M-step ridge 2.247e-3/2.176e-3, subspace cosines 0.9346/0.9600, dense error observation 6.638e-3, backward residual 1.904e-4, condition 692.5; flex PCG operator 7/7; application state-to-truth accuracy 0.9995 in all arms, hard-label counts `A=(1,5,994)`, `B=(1000,0,0)` over states 1:3; remaining: transcribe same-frame map/eigenvolume, workflow-parity and timing observations from a retained verbose log |
| 1 | working tree | 2026-10-07 | maintainer reports all tests green | green | unchanged by contract | Phase 1a removed the proven-dead list. Phase 1b replaces the single/paired probe-part codecs with one version-13 packed container codec, compacts worker stage ids to 1--4, centralizes the duplicated constants, removes the retired weight-range codec, and corrects the workflow gate's `get_vol` output arguments. Old part files intentionally fail the version check. |
| 2a | working tree | 2026-10-07 | maintainer reports all tests green | green | unchanged by contract | Removed all 20 environment keys and `flex_run_settings`; removed environment-only basis composition, optional merge policies, state-filter overrides, probe/calibration caps, deflation opt-outs and PCG overrides. Kept the reachable defaults: pairing 1, all-particle probe, 20k calibration, background+dilation deflation, default PCG settings, filtered state delivery, and complete-linkage two-gate merge under `preimage_auto=yes`. `rg SIMPLE_COV_ src` is empty. |
| 2b | working tree | 2026-10-07 | maintainer reports 14/14 tests green | green | unchanged by contract | Removed all six whole-`parameters` signatures and the `simple_parameters` dependency from the coupled-PCG leaf. Its support/window/operator helpers now receive explicit file, geometry, mask and band inputs. The shared utility leaf now receives `outvol` and the project orientation field directly, eliminating both `parameters` and `builder`; the projected-model band helper takes only `box_crop`, `smpd_crop` and `lp`; the latent-target module's stale `parameters` import is gone. |
| 2c | working tree | 2026-10-07 | maintainer reports all tests green | green | unchanged by contract | Closed the phase-2 review residue: state delivery has only the fixed per-state 0.143 low-pass path; merge uses one fixed 0.98 gate and one complete-linkage statistic; `pair_map_ratio` has no dead band-cap argument; the mod-4 split has no alternate pairing argument while persisted provenance remains 1; the dead `signal_subspace` helper is gone and `noise_and_projection` is private. Corrected the merge-loop indentation. |
| 3a | working tree | 2026-10-07 | maintainer reports all tests green | green | unchanged by contract | Moved `planes_batch_load` into latent ops beside particle preparation. The resident-plane module now exposes only enable/fetch/store/kill and imports no other FLEX module, `parameters`, `builder` or matcher. The copy/read/prep timing boundaries and resident-store counters are preserved. Regenerated graph is acyclic and has only `latent_ops -> {planes, plane_cache}`. |
| 3b | working tree | 2026-10-07 | static checks green; maintainer compile/tier rerun pending | pending | unchanged by contract | Moved the bandwidth-CV reconstruction, half-selection, map I/O and cleanup loop from `weights` to `state_service`; `weights` retains the bin-grid, opposite-half scoring and minimum-error adoption numerics and no longer imports `builder`, `parameters`, images, rounds or reconstruction. Added L0 `simple_defs_flex` for the three shared cross-layer constants. `scripts/check_flex_dag.py` covers 42 production modules and 137 edges, rejects multi-unit sources, upward imports and cycles, ignores tester modules by suffix, and runs before the fast gate. The review closure documents the lower-or-acyclic-same-layer rule and removes the redundant half-map filter copy/branch. |
| 4a | working tree | 2026-10-07 | static checks green; maintainer compile/tier rerun pending | pending | unchanged by contract | Plane-cache version 2 records the absolute project-path digest and modification stamp, project row count, geometry, exact master selection hash/count and highest completed row. Masters require every field; workers require matching provenance/geometry and a partition covered by the completed cache. `flex_pca_plane_cache%kill` resets the open unit, buffers, geometry, availability/probe flags and filename. The fast FLEX suite checks changed selection, changed path, in-place project replacement, worker coverage and both session orders. |
| 4b | working tree | 2026-10-07 | static checks green; maintainer compile/tier rerun pending | pending | unchanged by contract | Replaced the resident-plane and disk-cache singletons with one application-owned `flex_plane_store`. Its resident rows and allocatable `flex_pca_plane_cache` are private object state, passed explicitly through mean scaling, paired/single fit, E-step and embedding paths, and released together by application `kill`. Only probe, polish and embed workers adopt it; state workers avoid the unused 130--320 MB batch buffer. No particle-preparation arithmetic or cache format changed. |
| 4c | working tree | 2026-10-07 | static checks green; maintainer compile/tier rerun pending | pending | unchanged by contract | Completed object ownership: the application allocates its session and one canonical support environment; each fit allocates its polar bank and owns an exact copy of the already-resampled support, while the state merge borrows the application environment. This reduces one paired master run from repeated mask reads/resamples to one without changing support values. Cross-FSC file/record/context state has type-bound lifecycle; rounds own the artifact catalog threaded through every part codec; the validated weights set is a `flex_weights_store`. The pair-merge deflation helper accepts its fit as `intent(inout)`, as required by the fit-owned envelope lifecycle. Focused singleton/free-API scans are empty and the scoped production tree has no mutable module state. |
| 4 | working tree | 2026-10-07 | source-complete; static gates green; maintainer compile/tier rerun pending | pending | unchanged by contract | Phase 4 contracts are implemented. Cache contract tests are present; runtime result pending. |
| 5 | working tree | 2026-10-07 | source-complete; static gates green; maintainer rebuild/reconstruction/FLEX tier rerun pending | n/a | before/after runtime metrics unavailable | Added the rank-1 `simple_pcg_solver` engine and shared outcome. Reconstruction and coupled-FLEX adapters remap their native rank-3/rank-4 storage without copying the client solution; their operators, preconditioners, window conversions and client-specific policies remain local. The reconstruction adapter's explicit-shape target relies on its inherent contiguity rather than the illegal `CONTIGUOUS` attribute; the coupled adapter's assumed-shape declarations retain that attribute. A requested cold restart is always logged, preserving the former FLEX diagnostic at default iteration verbosity. Optimiser tests cover a dense SPD oracle, an exact one-step diagonal solve, pinned outcome fields and a deterministic cold retry. |
| 5b | working tree | 2026-10-07 | skipped per section 6 | n/a | no retained phase-5 before/after timing evidence | No speculative optimisation was attempted. |
| 6 | working tree | 2026-10-07 | source-complete; static gates green; maintainer rebuild/tier rerun pending | pending | unchanged by contract | Retained one source file per module or submodule: 45 baseline files plus the application tester gives 46 scoped files. Every production FLEX/strategy core-API import has an explicit `only:` list. Maintainer compiles exposed seven OpenMP runtime dependencies formerly inherited through that umbrella, duplicate use association in the four probe-fit submodules, and one omitted `simple_abspath` commander dependency; each owner now imports its exact local dependencies, and submodules rely on host association for names supplied by their parent. The DAG checker rejects multi-unit sources as well as broad core imports, upward edges and cycles. `src/main` uses a recursive glob without `CONFIGURE_DEPENDS`, so the maintainer must reconfigure once after the restored file set before building. |
| 7 | working tree | 2026-10-07 | static style/documentation checks green; maintainer compile/tier rerun pending | n/a | n/a | Updated per-file descriptors, README layout, generated indexes and code map. Fixed the index generator so `module procedure/subroutine/function` declarations cannot create false module entries. Removed the three labelled `CONTINUE` no-op targets from the refactored scope without changing their branch destinations or reindenting the routines. A maintainer compile exposed and removed a stray period after the `MEAN_SCALE_FNAME` identifier. The embed worker now materializes the polymorphic `rounds%part_fname` result before passing it onward, matching the safe pattern at every other call site. Procedure bodies and maintainer-fixed indentation were otherwise not reformatted. `git diff --check` and descriptor checks pass. |
| 8 | working tree | 2026-10-07 | source-complete; static gates green; maintainer compile/tier rerun pending | pending | unchanged by contract | Replaced the two byte-equivalent real-SPD/complex-RHS Cholesky solvers with `simple_linalg::solve_real_spd_complex` and added a direct SPD/non-SPD unit oracle. Cross-fit shell calculations already call `image%fsc`. The latent k-means remains local because its reliability-weighted metric, deterministic farthest-point seeding and empty-cluster recovery have no shared contract-compatible implementation. |

## 11. What the review of 2026-10-06 corrected

1. Baseline and scope: 45 files / 21,422 lines scoped; the class-expansion file excluded
   from every count. The application tester raises the scope to 46 files. Each module and
   submodule remains in its own source file; the 42 production modules remain cycle-free.
2. The hybrid E-step was misclassified as unreachable; it is the production E-step and is
   kept untouched. Only proven-dead code is removed.
3. The first proposed consolidation map created `application ↔ probe_fit`, `states ↔ delivery`,
   `basis ↔ fit_types` and `io ↔ weights` cycles. Physical consolidation was subsequently
   rejected; the plan now draws the DAG and fixes ownership while retaining every file boundary.
4. The PCG abstraction had no viable vector contract; 5.4 specifies the rank-1 contiguous
   contract with pointer remaps and tests the engine before migrating either client.
5. Plane-cache correctness (selection-blind adoption, stale lifecycle) is now a contract
   with tests, not an ownership move.
6. Test design: observation, oracle and regression margin separated (4.1); the invalid
   90° control withdrawn; the M-step toy shrunk to a sub-second oracle; PSD instead of PD;
   signed permutations for the cross-fit test. The later application-contract review removed
   the proposed synthetic in-process worker arm because no such API exists; real partition
   reduction remains in the workflow gate.
7. Registration: the `COMMAND` argument, the CTest budget decision (D1), the next-to-code
   application tester and the synchronised suite/UI registration.
8. `cv_select_bandwidths` orchestration leaves the weights domain module; the umbrella
   import cleanup is a phase item.
9. Instrumentation: thread-seconds and wall seconds defined per bucket before any timing
    drives a decision.
10. The smooth Gaussian former fixture carried negligible power on the polar rings, so its apparent
    accuracy only measured low-`k` dominance. Phase 0 now uses narrow blobs, normalizes each
    ring statistic by its own ring oracle, guards the retained ring-power fraction, and treats
    polar-versus-Cartesian quadrature as an observation. The energetic rerun exposed the distinct
    Cartesian-plane-to-polar and direct-volume-to-polar interpolation paths, so the exact
    synthesis identity is asserted on the Cartesian oracle while the production polar identity
    defect is recorded as a latent-bias risk requiring separate scientific investigation. The
    M-step toy now probes its Tikhonov lift explicitly and guards backward rather than
    condition-amplified forward error.
11. The proposed application fixture assumed unsupported behavior. A small phantom must lower
    `min_neff` explicitly; a one-state collapse cannot be published by the current project contract;
    and `run_worker` is not an in-memory partition API. The library test now stays on the supported
    four shared-memory backend combinations, while the real two-part comparison remains in the
    workflow gate. This avoids adding production features merely to support the refactoring test.
12. The first complete application run showed that a two-mode phantom need not collapse a provision
    ceiling of three to exactly two delivered maps: production retained the two endpoints and an
    intermediate representative. The invalid exact-count and two-id permutation checks are replaced
    by a bounded delivered-state count, an optimal state-to-truth recoding with a union-bound chance
    floor, and same-frame ordering checks on two distinct truth-matched maps. The tester now enters
    through the public commander, and logs the full truth-count and correlation matrices.
13. The phase-1b review removed the final unused weight-range codec and stale imports/comments,
    routed both worker fit counts through the same open/write/close path, and corrected the workflow
    gate to receive `get_vol` sampling and box outputs in variables of the right meaning and type.
    The phantom's third state is recorded as nearly empty (5 of 2,000 hard labels), not as a split
    conformation.
14. Phase 2 separates supported typed configuration from command-line provenance. Values that
    affect production numerics come from `parameters` or the owning module's documented default;
    the application alone asks whether `infile`, `vol1`, `pindfile`, or `npreimages` was explicitly
    supplied. The former environment-only basis-composition and policy variants are deleted rather
    than retained as dormant default-only branches.
15. The phase-2 closure review verified all 20 former environment defaults against the baseline and
    found residual branches that could no longer be reached: alternate state filters, a second merge
    statistic and adaptive-threshold scaffolding, pairing 3, and the old composition subspace helper.
    They are removed; the fixed pairing identifier remains in artifacts as provenance. The review
    also made the envelope initialization precondition explicit in the risks above.
16. SIMPLE's modern Fortran style does not use `CONTINUE` statements. The newly added state-target
    branch target and the two retained embedding targets now label their following executable
    statements directly; no control flow or surrounding indentation changed.
17. SIMPLE keeps one module or submodule per source file. The phase-6 co-location was reversed,
    all 15 source paths were restored with their current implementations, and the FLEX DAG gate
    now rejects any multi-unit source.
18. The final ownership review removed repeated support-mask reads/resampling by giving the
    application one canonical environment and cloning its exact image into resident fits; state
    workers no longer adopt an unused plane cache; the last direct polymorphic filename-result
    call is materialized; and FLEX cold-restart messages remain visible at default verbosity.
19. The one-module-or-submodule-per-file rule is intentionally repository-wide for first-party
    Fortran, not only a FLEX constraint; `check_descr.py` remains in the fast gate. The FLEX DAG
    parser now scans imports in contained procedures and recognizes intrinsic/non-intrinsic,
    continuation and rename forms, so alternate syntax cannot evade layer or broad-core checks.
