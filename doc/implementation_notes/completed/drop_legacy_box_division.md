# Dropping the legacy box division from reconstruction output

**Status:** in progress on branch `drop_legacy_box_division` (started
2026-08-22). Steps 1 (instrumentation, §5.1) and 2 (dual-backend test, §5.2)
are implemented. **The premise of §2 is refuted by the step-1 baseline and by
the projector code (§0 below); step 3 is revised. Decision 2026-08-22: do
the deapodization fix (§5.3) first; retiring the ÷box/×box pair together
(§5.3a) was then done on 2026-08-22 as well; both verified (§5.3).**

## 0. Revision (2026-08-22): the division is one half of a matched pair

Two findings after steps 1 and 2:

1. **Real-data gridding baseline** (streptavidin abinitio3D, `rec_backend=
   gridding`, 86 iterations in `stage1_rec_backend_gridding/6_abinitio3D`):
   ref/ptcl amplitude ratio 0.24–0.55 in the lowest band, falling with
   resolution to 0.006–0.3 at the matching limit; euclid `v` q05/q50/q95 ≈
   0.84/0.89/0.95 in stage 1 and 0.70/0.77/0.89 at the end of stage 2; no
   invalid particles. This is the healthy regime of §3, not the
   "refs ≪ particles by ~box" regime §2 predicts for the ÷box convention.
2. **The projector multiplies by the original box.** `simple_projector::
   expand_cmat(orig_box)` sets `factor = real(orig_box)` and stores
   `cmat_exp = factor × cmat`; every reference path goes through it
   (`read_mask_filter_refvols` → `expand_cmat(params%box)` →
   `vol_pad2ref_pfts_opt`, and the batch projectors; `simulate_particles`
   and `reproject` likewise). Under SIMPLE's `fft = (1/N)Σ` convention the
   2D transform of a projection equals `box ×` the central slice of the 3D
   transform (one extra summed dimension), so a volume whose Fourier
   coefficients sit at the particle-coefficient scale (the plain data
   quotient) must be divided by `box` to be stored as a volume, and its slice
   multiplied by `box` to become a reprojection. `reconstructor_eo%
   mag_correction = box` and `expand_cmat`'s `factor = box` are that pair.
   They cancel on the reference path, so the map convention has NO effect on
   the euclid objective; the §1 statement that "nothing else in the chain
   introduces a box-sized factor" overlooked the projection side.

Consequences:

- The PCG crash mechanism stands, but its reading changes: PCG's undivided
  solution was projected with the ×box factor still applied, so references
  came out box× too large. `pcg_mag_correction` is therefore the correct
  convention for a volume that will be projected by `expand_cmat`, not a
  kludge, and the gridding↔PCG handoff is scale-continuous with it (step-2
  test: in-band ratio 1.09, FSC 0.95–0.99).
- **Step 3 as originally written (drop the division only) is wrong**: it
  would put every gridding reference box× too large, i.e. reproduce the
  crash in the production path. If the division is to be retired, the
  projector's multiplication must be retired with it; the change is then
  net-neutral for refinement (same references, same sigma2, same `v`) and
  only changes map values on disk by ×box. That is a cosmetic/robustness
  refactor (it removes the trap a new backend can fall into), to be decided
  on its own merits. The acceptance criterion is equality, not improvement:
  EUCLID DIAG and `test_rec3D_backends` output identical before and after.
- The gridding deapodization mismatch of §2.1 is unaffected by this
  revision and remains a real, separate item (§5.3 of the revised plan).
- §2's explanation of the 2D class averages being "already at the
  data-quotient scale" is consistent: there is no slice theorem in 2D, so
  `cavger` carries neither factor; the 2D control's ratio ≈ 1.1 at low
  resolution vs 3D's ≈ 0.5 reflects that and the 3D filters/masks, not a
  convention error.
- Open: why the synthetic 3D case (§5.1, SNR 0.1 simulated particles) shows
  ratio 0.08 where real data shows 0.5. Both use the same projector, so it is
  not the convention; likely `simulate_particles`' SNR definition and/or the
  static `lp=10`. A 2D control on the real particles would settle what the
  "expected signal" ratio is for that dataset. Until step 3 lands, the PCG backend mirrors the convention
(`pcg_mag_correction` in `simple_rec3D_pcg_strategy.f90`) so the two backends
are interchangeable.

## 1. What the convention is

`reconstructor_eo` divides every real-space map it produces by the ORIGINAL
box size before writing it:

```
simple_reconstructor_eo.f90:28   real :: mag_correction=1.  !< scaling factor to correct for slice insertion, cropping & padding
simple_reconstructor_eo.f90:119  self%mag_correction = real(self%p_ptr%box) ! consistent with the current scheme
simple_reconstructor_eo.f90:522/532/571/582/601/602   call even/odd%div(self%mag_correction)   (sampl_dens_correct_eos)
simple_reconstructor_eo.f90:770  call self%eosum%div(self%mag_correction)                        (sampl_dens_correct_sum)
```

Note that it is `params%box`, not `box_crop`: since Fourier cropping is
amplitude-neutral in SIMPLE's transform convention (`fft` = (1/N)Σ, `ifft` =
Σ), dividing by the original box keeps the output scale independent of the
crop factor. Nothing else in the reconstruction chain introduces a box-sized
factor: the padded 2D plane is rescaled by `pf²` back to the native
coefficient convention in `insert_plane_oversamp`, the KB weights are
normalized per axis and enter numerator and `rho` alike, `sampl_dens_correct`
is a pointwise quotient, and the inverse transform is on the native lattice.
So the map gridding writes is

    map_gridding = deapod( ifft( Σ w·ctf·y/σ² / Σ w·ctf²/σ² ) ) / box

i.e. the plain data quotient, divided by 128 for a 128-pixel box. The
division is a retained historical normalization ("consistent with the current
scheme"), not a property of the reconstruction.

The same division is applied by hand in:

```
src/main/flex/simple_flex_pca_columns.f90:1243   call mean_rec%div(real(params%box))
src/main/flex/simple_flex_pca_columns.f90:2054   call work%div(real(params%box))
src/main/flex/simple_flex_pca_columns.f90:3688   call Yeven(q)%div(real(params%box))
src/main/flex/simple_flex_pca_columns.f90:3694   call Yodd(q)%div(real(params%box))
src/main/flex/simple_flex_pca_rec3D.f90:265      call cur(state)%div(real(params%box))
src/main/commanders/simple/simple_commanders_volops.f90:605   call noisevol%div(real(params%box))
src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90   pcg_mag_correction / apply_output_convention
```

The 2D class-average restore (`simple_classaverager_restore.f90`) does NOT
carry an analogous division (it divides by the class population only), so
class averages are already at the data-quotient scale.

## 2. How it was found, and what it costs

The PCG backend solves the normal equations in the solver's native
convention; all its internal factors (`padsc` in `b`, `padsc²` in `Khat`, the
padded-lattice interpolation, the normalized KB weights, the deapodization
bracket on both `H` and `b`) cancel, and its solution is the plain data
quotient — without the `/box`. Measured on the streptavidin run
(`pcg_euclid_crash_investigation.md` (consolidated into `pcg_priors.md`) §15): gridding stage-1 lp map L2 0.99 vs
PCG iteration-41 lp map L2 194; per-shell ratios 160–260 over the shared band,
i.e. 128 × ~1.5, the remainder being the different alignment state of the two
maps (the ratio rises with resolution) and second-order KB-deconvolution
differences. This ×128 step is what overflowed the euclid objective at the
gridding→PCG handoff.

The cost of the convention is in the euclid objective. The reference
reprojection is, physically, the expected signal in a noise-normalized
particle, and at the data-quotient scale that is what it is: in the
restart run driven from a PCG-scale volume, polar reference amplitudes were
0.0055 against particle amplitudes 0.012 — a reference at the signal fraction
of the particle, as a least-squares objective assumes. At the gridding scale
the references come out ~128× below that, so every euclid refinement runs in
a refs ≪ particles regime: `v = 1 + crvec/norm` sits just above 1, the
probability contrasts are compressed by ~128² in `crvec`, and the group
sigma2 estimates reduce to particle power because the reference barely
enters the residual. It works — SIMPLE has refined in this regime for years —
but it is not the regime the objective was written for, and it is the reason
a scale-consistent backend could not be dropped in.

### 2.1 The gridding deapodization is a padded-era formula applied to a native-lattice reconstruction

The gridding deapodization lives in `volassemble`: `restore_state_from_parts`
multiplies the merged volume (`restore_merged_volume`) and the
nonuniform-filter source halves by `gridcorr_img = prep3D_inv_instrfun4mul(
ldim, ldim_pd, smpd_crop)` (`simple_gridding.f90`), the inverse of the
continuous KB instrument function `kbinterpol%instr` evaluated at
`pad_sc = 1/ldim_croppd` — i.e. for a window living on the 2x PADDED lattice
("mirror original intent"). The reconstructor, however, now inserts on the
NATIVE lattice (`reconstructor_eo%new`: `ldim = [box,box,box]`; the
`clip_inplace([box,box,box])` after `ifft` in `sampl_dens_correct_eos` is a
no-op left over from the padded era) with a KB window of half-width 1.5
NATIVE voxels and per-axis-normalized 3-point weights. Its true real-space
envelope is therefore the transform of that normalized stencil with period
`box`, and the applied correction is the one for period `2*box`:

| native offset r from centre (box 88) | true envelope E_grid (period 88) | applied correction envelope (instr, period 176) | residual attenuation R_grid |
|---|---|---|---|
| 10 | 0.924 | 0.982 | 0.94 |
| 20 | 0.735 | 0.928 | 0.79 |
| 30 | 0.523 | 0.843 | 0.62 |
| 40 | 0.394 | 0.736 | 0.53 |

(PCG's own envelope, built by `build_kb_envelope_1d` on its padded lattice,
is 0.98/0.92/0.84/0.73 at the same offsets and is removed exactly by the
solver's bracket; numerically it coincides with the continuous `instr` at
`1/boxpd`, which is why the padded-era formula was right in its day.)

So gridding maps are under-deapodized: density fades toward the box edge by
up to ~2x (R_grid ≈ 0.53 on axis at r = 40, ≈ 0.55 at the corner of the
masked sphere). For a small particle like streptavidin (radius ≲ 20 px) the
effect inside the molecule is ~1.1–1.2; for a box-filling particle it is a
real radial amplitude error in every gridding map and every reprojection.
The even/odd half-map files written by `sampl_dens_correct_eos` receive no
correction at all; the FSC (halves vs halves) is unaffected, but the
non-`lpset` reference path (`read_mask_filter_refvols` reads the halves when
`lp` is not set) matches against apodized references.

Measured on same-alignment maps (gridding stage-2 snapshot, iteration 40, vs
PCG iteration 41; `envelopes.py`/`compare2.py` in the investigation
scratchpad, recorded in `pcg_euclid_crash_investigation.md` (consolidated into `pcg_priors.md`) §15.4): L2 ratio
inside the molecule (r < 24) 175–177; the box convention alone predicts 128;
box convention × gridding's envelope deficit predicts 140–156; the remaining
×1.2–1.26 rises with spatial frequency (shell ratio 136 at shell 2 → 184 at
shell 9 after the /128) and with proximity to the core, consistent with the
shift search switching on at stage 3 (`trs = 0` in stages 1–2) sharpening the
iteration-41 map — a property of the maps, not of the conventions.

The correct gridding deapodization is the discrete normalized-stencil
transform with period `box` — exactly what `build_kb_envelope_1d` computes
for the PCG operator — not the continuous `instr` at any scale (the
continuous function disagrees with the discrete stencil by ~2x at the box
edge: 0.21 vs 0.38 at r = 43 for period 88; at the padded scale the two
coincide, 0.701 vs 0.702 at the native box edge, which is why the 2D path
below is fine). Fixing it belongs with the convention change (§5.3): it
changes every gridding map's radial profile.

Same formula, same fix, in the flex reconstruction paths:
`simple_flex_pca_columns.f90:1244` and `:2000` call
`prep3D_inv_instrfun4mul([box_crop]*3, OSMPL_PAD_FAC*[box_crop]*3, ...)` on
reconstructor-object output (native lattice).

The 2D class-average restore is CONSISTENT as it is: `cavger` accumulates on
the padded lattice (`simple_classaverager_restore.f90:~945`,
`loc = matmul([hp,kp],mat)` with `hp = h*OSMPL_PAD_FAC`), so
`prep2D_inv_instrfun4mul(ldim_crop, ldim_croppd, ...)` at `1/ldim_croppd` is
the right envelope for its window. No change needed there.

Known approximation shared by PCG (today) and gridding (after the fix): the
envelope is the transform of the stencil at fractional offset 0, whereas the
per-sample normalization (`apod_mat_3d_fast`, `kb_apod_vecs_3d_fast`) makes
the effective stencil offset-dependent (at offset 0.5 it is [0, 0.5, 0.5]).
The true envelope is the offset average; the difference is second order and
identical for both backends, so it cannot create a backend discrepancy. It
is noted here so nobody later mistakes it for one.

### 2.2 Adjacent, not required for backend equivalence: the projection side

Reference reprojection (`read_mask_filter_reproject_refvols` →
`pad_fft` → `expand_cmat` → `vol_pad2ref_pfts_opt`, and `simple_projector`)
KB-gathers from the padded volume's FFT without pre-dividing the real-space
volume by the gather envelope (no use of `kbinterpol%instr` or an inverse
instrument function exists outside `simple_gridding`, which only the
reconstruction side calls). The interpolated slices are therefore those of
the volume multiplied by the padded-lattice envelope (0.98/0.92/0.84/0.73 at
r = 10/20/30/40 for box 88): references of a box-filling particle are
attenuated toward the periphery by up to ~25%. This is identical for both
backends and for the class averages, so it does not enter the backend
comparison; it does enter the absolute ref/particle ratio of §3 for large
particles. The fix, if wanted, is to multiply the prepared reference volume
by the inverse padded-lattice envelope (the same `build_kb_envelope_1d`
construction with period `box_croppd`) before `pad_fft`.


## 3. Target convention

Every reconstruction backend writes maps at the **data-quotient scale**: the
Fourier coefficients of the map are the CTF-corrected, sampling-density-
normalized (and deapodized) average of the particle Fourier coefficients in
SIMPLE's (1/N)Σ convention. Operationally: `mag_correction = 1`, and the
PCG output path applies no factor.

The verifiable consequence, and the acceptance test: for a converged
refinement, the polar reference reprojections must have amplitudes of the
order of the signal fraction of the noise-normalized particle polar
transforms (ref/ptcl ~ 0.1–0.5 per shell inside the matching band, falling
with resolution), and the euclid `v` values must sit clearly below 1. Since
`v = Σ_k (k/σ²_k) Σ_p |ptcl − CTF·ref|² / Σ_k (k/σ²_k) Σ_p |ptcl|²` (§5.1), a
reference that explains particle variance gives `v < 1` (`v = 0` for a perfect
noise-free match, `v = 1` for a zero reference, `v = 2` for an uncorrelated
reference of equal power); `v ≈ 1.000 across the table` is the refs ≪ particles
regime, `v > 1` means the reference adds more power than it explains, and the
23.03 invalidation threshold is >20 particle powers of residual. (An earlier
version of this section said "typically 1–3"; that was written before the
identity above was established and is wrong.)

## 4. Consumers of absolute map amplitude — audit

| consumer | scale behaviour | action |
|---|---|---|
| euclid objective (`gen_prob_likelihood_euclid_val`, `crvec`/`norm`) | scale-sensitive by design | the point of the change; validate per §3 |
| group sigma2 estimation (`calc_group_sigmas`) | residual now contains the reference | values drop toward true noise power; the bootstrap (`calc_pspec`) is particle-based and unaffected |
| prob-table invalidation (`euclid_dist_from_crvec`, `v > 23.03`) | absolute threshold on `v` | safe at the correct scale; the threshold only bites at >100× mis-scaling |
| cc objective | invariant | none |
| ML regularization (`add_invtausq2rho`, PCG `ml_prior`) | built from FSC and rho only | none |
| reference noise regularization (`regularize_ref_with_noise` → `add_gauran(snr)`) | noise sdev = sqrt(var/snr), relative | none (verified) |
| automasking / NU evidence envelope | thresholds on map statistics — verify whether any absolute level is used | audit `simple_image_msk%automask3D`, `nu_evidence_envelope` |
| nonuniform filter (`simple_nu_filter`) | filter-bank selection on correlation-type criteria | verify; expected scale-free |
| postprocess (B-factor, low-pass, masking) | scale-free | none |
| `volops` `noisevol%div(real(params%box))` | simulation noise calibrated to the map convention | change together with the maps |
| flex PCA reconstruction paths (`simple_flex_pca_columns/rec3D`) | consume maps in the map convention | change together |
| `abinitio3D_cavgs` (volume vs class averages) | class averages are at the data-quotient scale already and correctly deapodized (§2.1) | check that the cavgs↔reprojection comparison does not compensate the old factor anywhere |
| trailing reconstruction accumulators / PCG raw (B,D) files | pre-division sufficient statistics | unaffected |
| volumes on disk from earlier runs used as starting references | at the old scale (÷box) | first iteration runs refs ≪ particles, recovers after one sigma2 update; document, optionally warn when the ref/ptcl ratio is ≪ 1 |
| external viewers (isosurface thresholds) | map values ×box | cosmetic; release note |

## 5. Plan

1. **Instrument first (no behaviour change).** Log, once per iteration, the
   polar reference/particle amplitude ratio (median over particles, per shell
   band) and the quantiles of euclid `v` — the same quantities the temporary
   `EUCLID DIAG` reported during the investigation, made permanent and
   compact. Record the baseline on gridding as it is today.

   *Implemented (2026-08-22).* `gen_sigma_contrib` (`simple_polarft_corr.f90`)
   takes optional outputs `ref_pow(k)` (mean |CTF·ref|² per polar component),
   `ptcl_pow(k)` (mean |ptcl|²) and `v`, all at the assigned orientation it
   already builds for the sigma2 residual. `v` is computed directly as
   `Σ_k (k/σ²_k) Σ_p |ptcl − CTF·ref|² / Σ_k (k/σ²_k) Σ_p |ptcl|²`, which is
   algebraically what `euclid_dist_from_crvec` returns as `1 + crvec/norm`
   (the `2·nrots` in `norm` cancels the unnormalized FFT round trip and the
   conjugate-symmetric half), so it is valid in prob-align mode too, where the
   memoized FFT tables are not populated. `euclid_sigma2%calc_sigma2` stores,
   per particle, the amplitude ratio `sqrt(Σ_band ref_pow / Σ_band ptcl_pow)`
   in `NDIAG_BANDS = 4` equal bands of the search range `params%kfromto`, and
   `v`; `write_sigma2` (end of every 2D and 3D iteration) calls
   `report_euclid_diag`, which prints two lines on part 1 only:

   ```
   >>> EUCLID DIAG ITER 12 NPTCLS 4321 KFROMTO 1-29 REF/PTCL AMP (q50) k[1-7]: 3.1E-01 k[8-14]: ... k[22-29]: ...
   >>> EUCLID DIAG ITER 12 V q05: 1.0210 q50: 1.1503 q95: 1.6400 max: 2.31 THRES: 23.03 NINVALID: 0
   ```

   `NINVALID` counts particles whose `v` exceeds the `-log(TINY)` threshold at
   which `euclid_dist_from_crvec` returns `huge`. Only particles visited by
   `calc_sigma2` in the iteration (i.e. updated, state ≠ 0) enter the
   quantiles. The 2D report is the control: class averages are already at the
   data-quotient scale (§1), so `cluster2D` shows the regime the 3D path
   should reach after step 3. Cost: one extra per-shell sum per particle per
   iteration, negligible.

   First reading (synthetic check, 2026-08-22; 400 particles simulated from
   the streptavidin `rec_final_state01.mrc` at SNR 0.1 with CTF, `refine3D`
   gridding backend, `lp=10`, 3 iterations): ref/ptcl amplitude ratio
   ≈ 0.08 and flat over shells 2–14; `v` q05/q50/q95 ≈ 0.85/0.87/0.89, max
   0.92, no invalid particles. This is a self-consistent synthetic case (the
   particles were simulated from the reference), not the gridding baseline
   on real data; that baseline is still to be recorded. The 2D control on
   the same particles (`abinitio2D`, `ncls=10`, 24 iterations) shows the
   regime the plan expects of correctly scaled references: ratio ≈ 1.1 in
   the lowest band falling to ≈ 0.05–0.1 at the matching limit, `v`
   q05/q50/q95 ≈ 0.4–0.8/0.55–0.98/0.7–1.2 with a real spread across
   particles (iteration 1, random classes: ratio 0.008, `v` 0.986). The 3D
   gridding case's flat ratio 0.08 with `v` compressed to 0.84–0.89 is the
   refs ≪ particles regime of §2: `v ≈ 1 − 2·0.08·corr + 0.08²`.

   Testing on this branch: the distributed workers are launched from
   `$SIMPLE_PATH/bin/simple_private_exec`, and `build/bin` is only refreshed
   by `make install` (install prefix = build dir). For this branch's build to
   be the one that runs, do `cd build && make -j && make install` and point
   `SIMPLE_PATH` at `pcg_integration_bug/SIMPLE/build` in the shell that
   launches the run; the shell-memory master relays the worker's report into
   its own log, so the lines appear in the main output with `nparts=1`.
2. **Same-inputs dual-backend test.** A test commander that reconstructs one
   fixed set of particles/orientations with both backends and compares shell
   profiles (`shellspec`-style). With `mag_correction` removed and the
   gridding deapodization corrected the ratio must be ≈1 across the band and
   flat in radius; with only `mag_correction` removed the radial profile of
   the ratio must follow 1/R_grid of §2.1.

   *Implemented (2026-08-22).* `simple_test_exec test=rec3D_backends
   projfile=<run>.simple pgrp=.. mskdiam=.. nthr=..` (`exec_test_rec3D_backends`
   in `simple_commanders_test_highlevel.f90`). Run it inside a `refine3D` run
   directory: it executes the production `reconstruct3D` commander twice in
   the cwd (gridding, then pcg; same project, orientations and
   `sigma2_it_N.star`), keeps the merged maps as
   `recvol_stateXX_gridding.mrc` / `recvol_stateXX_pcg.mrc`, and prints (a) a
   shell table — mean Fourier amplitude per shell for both, ratio pcg/gridding,
   FSC between the two maps; (b) a radial table — mean |ρ| in 16 radial bins,
   ratio, ratio normalised to the centre bin; (c) a summary — the agreement
   band (contiguous shells from k=2 with FSC(gridding,pcg) > 0.5), the median
   shell ratio inside it, and the min/max of the normalised radial ratio over
   bins fully inside 0.85×mask radius. No thresholds are enforced; it is the
   measurement for step 3.

   Caveat for the radial test: the PCG solver solves on a spherical support
   (`pcgop%set_mask(params%msk_crop)`), so its map is zero outside the mask
   radius and soft at the edge, while gridding's is not; only bins well inside
   the mask are comparable, hence the 0.85 factor. Bin 1 (r < box/32) holds
   few voxels and is noisy.

   First reading (same synthetic project as §5.1, box 128, `mskdiam=60`,
   PCG 2 iterations): agreement band k = 2–38 (3.6 Å), median in-band
   amplitude ratio pcg/gridding **1.09** (1.01–1.12 for k = 3–30), FSC
   between backends 0.95–0.99 in band — the convention mirror works. Beyond
   k ≈ 31 gridding amplitudes fall much faster than PCG's and the ratio
   explodes past k ≈ 40 with FSC ≈ 0: the PCG beyond-band excess of §6, now
   visible in a same-inputs comparison. Radial: normalised ratio 0.78–0.81
   flat over r = 4–24 px (centre bin 1.0, 18 % above), i.e. no rise with
   radius on a particle of radius ≲ 26 px in a 128 box — the §2.1 prediction
   (R_grid ≈ 0.8 at r/box ≈ 0.19) is not resolved by this case; a
   box-filling particle is needed to test it.
3. **[REVISED, see §0] Retire the ÷box / ×box pair together, or leave both.**
   If retired: `mag_correction` → 1 (six `div` sites), `pcg_mag_correction`/
   `apply_output_convention` and the two warm-start multiplications deleted,
   the flex and `volops` divisions of §1 removed, AND `expand_cmat`'s
   `factor` → 1 (all callers pass `box` or `ldim(1)`; audit `simulate_
   particles`, `reproject`, `symanalyzer`, `volpft_corrcalc`, `volinterp` for
   any consumer that relied on the product). Acceptance: EUCLID DIAG and
   `test_rec3D_backends` tables identical before/after (values on disk ×box).
   The deapodization fix below is independent of this choice.

   **Deapodization fix — implemented and verified against ground truth
   (2026-08-22).** `simple_gridding` now owns `kb_stencil_envelope_1d`
   (exact 1-D transform of the normalized origin stencil, period n),
   `kb_stencil_inv_envelope_1d`, `deapodize3D_inplace` and
   `prep3D_inv_kbenvelope4mul`; `prep3D_inv_instrfun4mul` is gone. The PCG
   reconstructor's `build_kb_envelope_1d` delegates to the shared routine
   (period `boxpd`; `test=pcg_recon` stages 1–9 still pass). For gridding the
   correction lives in `reconstructor_eo` (`deapodize`, period `box_crop`,
   built in `new`) and is applied to the merged map in
   `sampl_dens_correct_sum` and to BOTH half-maps (`_unfil` and filtered) in
   `sampl_dens_correct_eos`, so halves and merged share one convention and
   the non-`lpset` reference path reads deapodized halves. `volassemble`
   (`restore_state_from_parts`) no longer multiplies the files it reads;
   the two flex sites and `simple_flex_pca_rec3D` use
   `prep3D_inv_kbenvelope4mul`. The 2D `prep2D_inv_instrfun4mul` is unchanged.

   Verification (`test=rec3D_backends ... vol1=<truth> lp=10 hp=50`, ground
   truth = the volume the particles were simulated from; per-radial-shell
   least-squares scale recon/truth after background removal, normalised to
   r = 4–8 px; the test also re-deapodizes the gridding map analytically
   with the legacy padded-period instrument function for a same-run
   before/after):

   | r (px) | 0–4 | 8–12 | 12–16 | 16–20 | 20–24 |
   |---|---|---|---|---|---|
   | gridding, new envelope | 0.97 | 1.03 | 0.99 | 1.01 | 0.95 |
   | gridding, legacy envelope | 0.96 | 0.94 | 0.86 | 0.89 | 0.81 |
   | PCG | 0.93 | 1.01 | 0.92 | 0.96 | 0.86 |

   The legacy maps fade 14–19 % by r = 12–24 px (box 128); the corrected
   maps are flat within ±5 %. The offset-averaged effective kernel was also
   computed and differs from the origin stencil by < 0.5 % (0.965 vs 0.963
   at r = 10), confirming the "second order" note above.

   Two further findings from the same test: (i) without `hp` the gridding
   map appears to fade even after the fix (0.91/0.88/0.79/0.71) — that is
   NOT an envelope but a pre-existing low-frequency artifact: gridding's
   amplitude in Fourier shells k = 1–2 (137 / 69 Å) is ~25 % above the
   truth-relative level of every other shell (ratio 15.8 / 16.4 vs ≈ 12.8),
   PCG's is uniform (12.1 / 12.9 vs ≈ 12.5); with FSC to truth ≥ 0.99 in
   those shells it is an amplitude excess, not a phase error. It shows up as
   a broad positive central blob with a negative ring at r = 32–44 px and
   mimics a radial fade in any |ρ|-type comparison. Open item, §6.
   (ii) `EUCLID DIAG` is unchanged by the fix within noise (synthetic:
   ratio 0.083 → 0.081, `v` 0.86 → 0.86) — expected, since the envelope
   enters reprojections only through the map's radial profile.

   **§5.3a — the pair retired (implemented 2026-08-22).** `expand_cmat` no
   longer takes a box argument and no longer scales (`cmat_exp` holds the
   volume's own coefficients); `reconstructor_eo%mag_correction` and its
   seven `div` sites are gone; the PCG strategy's `pcg_mag_correction`,
   `apply_output_convention` and the two warm-start multiplications are gone;
   the five flex divisions and `volops`' `noisevol` division are gone; all
   `expand_cmat(box)` callers in `src` and `production/tests` updated.
   `report_euclid_diag` now warns when the low-band reference/particle ratio
   is < 0.02 with `v` ≈ 1 (an old-convention starting volume). Acceptance
   (equality) on the synthetic project: pre-change build on an old-convention
   map vs post-change build on the same map × 128, 2 euclid iterations —
   EUCLID DIAG 1.190/1.070/0.970/0.903 `v` 0.2088 vs 1.190/1.070/0.972/0.905
   `v` 0.2088 (iteration 2 likewise); `test=rec3D_backends` on the neutral
   phantom: every ratio and FSC identical, amplitudes × 128; `pcg_recon`,
   unit tests and the refine3D recovery suite pass.

   Release notes for the convention change:
   - Reconstructed maps (gridding and PCG, merged and halves, flex) are now
     `box` × larger in absolute value than before; isosurface thresholds
     scale accordingly.
   - A starting volume from an older run reprojects `box` × too weak; the
     first sigma2 update absorbs it, but the first iteration aligns against
     near-zero references and, if the resulting FSC collapses, the
     FSC-filtered references stay near zero thereafter (observed on the
     synthetic project: overlap 0.1, cFAR 3e-4, refs 6e-4 for all
     iterations). Scale old volumes by the box size before use; the EUCLID
     DIAG warning flags the condition.
   - `rotvol`/`symmetrize_map`/`reproject`/`simulate_particles` outputs were
     `box` × the input volume's scale (nothing compensated; harmless in
     `abinitio3D` because the symmetrized map is only used to map
     orientations before re-reconstruction); they are now scale-preserving.
     Verified: `reproject` gain 22.7 (old build, = 128 × 0.177) → 0.177.

   *Original text, kept for the record:* **Remove the division and fix the gridding deapodization.**
   `mag_correction` → 1 (or delete the member and its six `div` sites); delete
   `pcg_mag_correction`/`apply_output_convention` and the two warm-start
   multiplications in the PCG strategy; the flex and `volops` sites in §1.
   Replace `prep3D_inv_instrfun4mul`'s continuous `instr` at `1/ldim_croppd`
   by the discrete normalized-stencil envelope with period `box` (hoist
   `build_kb_envelope_1d` from `simple_reconstructor_pcg` into
   `simple_gridding` and use it for both backends), apply it to the
   half-map files as well so halves and merged share one convention, and
   use it at the two flex sites (§2.1). The 2D `prep2D_inv_instrfun4mul` is
   already consistent with its padded-lattice accumulation and stays.
   Optionally, in the same pass, add the projection-side pre-compensation of
   §2.2.
4. **Audit the scale-sensitive consumers** in §4 (automask, NU filter,
   `abinitio3D_cavgs`), fixing any absolute level tuned under the old scale.
5. **Validate.** Streptavidin (10335) plus at least one other dataset:
   abinitio3D, refine3D (gridding and PCG), refine3D_auto, abinitio3D_cavgs.
   Acceptance: §3 ratios and `v` ranges; FSC/resolution and orientation
   overlap not worse than baseline; sigma2 magnitudes consistent with
   particle noise power; gridding↔PCG switch mid-workflow with no scale
   discontinuity (the `PCG SCALE CONTINUITY`-type guard must not be needed
   and is not reinstated).
   *Stage-2 results (2026-08-22, streptavidin, commit `285f55d2` = step 2
   tool but BEFORE the deapodization fix; in `stage2_tests/`):*

   - `abinitio3D rec_backend=pcg` (PCG from stage 3, iteration 41 on) vs the
     stage-1 gridding run: low-band ref/ptcl ratio and `v` q50 at iterations
     40 / 41 / 60 / 72 / 85 / 100 — gridding 0.458/0.545/0.527/0.523/0.51/—
     and 0.83/0.79/0.77/0.77/0.85/—; PCG 0.462/0.546/0.623/0.638/0.60/0.52
     and 0.83/0.80/0.77/0.76/0.84/0.94. No discontinuity at the backend
     handoff (iteration 41 identical), no invalid particles in 220 iterations,
     `v` identical to within 0.01; PCG references come out ~20 % stronger in
     the lowest band during stages 3–4 (see the low-shell excess below).
   - `abinitio2D` control on the same particles: low-band ratio 0.44–0.57
     from iteration 4 on, `v` q50 0.78–0.94 — the same level as 3D
     (0.46–0.66), which settles §6's second question: the polar
     normalizations put the expected signal at ≈ 0.5 of the particle
     amplitude at low resolution for this dataset, for cavgs and
     reprojections alike.
   - `test=rec3D_backends` inside the gridding run dir (pre-fix, and before
     the comparison soft-masked both maps): FSC(gridding,pcg) 0.96–0.99 for
     k = 4–40, in-band amplitude ratio pcg/gridding 0.80–0.87 (partly the
     unmasked gridding map's outside-mask power), **pcg/gridding 1.7–2.0 at
     k = 1–3** (137–46 Å) — the low-shell discrepancy between backends is
     real on real data too, with the opposite sign to the synthetic case
     relative to truth. Radial pcg/gridding 1.03 → 0.94 over r = 4–24 px,
     centre bin 0.50 (gridding's central blob). To be rerun with the fixed
     build.

   *Stage-3 result (2026-08-22, same run dir, commit `11838829` = with the
   deapodization fix; `stage3_gridding_vs_pcg.txt`):* PCG map identical to
   stage 2 (s3/s2 = 1.00 ± 0.01 at every shell, as it must be). Gridding
   amplitudes vs stage 2: ×1.66 / 1.25 / 0.97 at k = 1–3, ×1.09–1.19 at
   k = 4–15, ×1.00 at k = 20–30 — the envelope's gain inside the 37 px mask,
   as predicted. Consequently pcg/gridding in band moved from 0.83–0.87 to
   0.75–0.86 (median 0.835 → 0.835 over the agreement band, which now
   extends to k = 48); FSC(gridding,pcg) 0.97–0.996 for k ≥ 5; shells 1–3
   still pcg/gridding 1.10 / 1.60 / 1.76. On the neutral synthetic truth the
   two backends agree within ±5 % in band, so the real-data 20–25 % in-band
   deficit of PCG relative to gridding is data/settings-dependent, not a
   convention or envelope effect. Prime suspect: convergence — the PCG
   solves report `ITS=2 RESID=4.6E-01` (46 % residual after the 2-iteration
   budget on N = 2905 particles), whereas the synthetic solves converge;
   second suspect: the two ML regularizations (gridding's FSC-tau in rho vs
   PCG's `KIND=ml` prior). Next real-data runs: `maxits_pcg=20 rtol=1e-3`
   and `ml_reg=no`, one at a time.

   *Stage-4 results (2026-08-22, real data, same run dir, deapodization fix
   in, convention pair still present; `stage4_test1.txt`, `stage4_test2.txt`):*
   `maxits_pcg=8 rtol=1e-3` → PCG `RESID` 0.46 → 0.18, in-band median
   pcg/gridding **0.835 → 1.11**: the in-band deficit was PCG's 2-iteration
   budget. `ml_reg=no` (2 iterations) → pcg/gridding 0.55 in band, FSC
   between backends negative at k ≤ 3 and < 0.85 below k = 9: without its ML
   prior the 2-iteration PCG solve is badly under-converged; the prior was
   compensating. The PCG-side action (§6) is therefore convergence control
   in refinement (iterations/tolerance), not the amplitude convention.

6. **Compatibility.** Release note on the ×box change of map values; a
   warning when a starting reference reprojects ≪ the particle scale.

## 6. Open items it does not settle

- ~~Gridding low-shell amplitude excess~~ — RESOLVED 2026-08-22, not a
  gridding defect. The "+25 % at k = 1–2" of §5.3 (ii) was an artifact of
  the truth volume's provenance: `rec_final_state01.mrc` of the first
  successful run is a PCG product (stages ≥ 3 and the final reconstruction
  ran `rec_backend=pcg`), so its lowest shells carry PCG's own low-k
  treatment and PCG "reproducing" it exactly proved nothing. With a neutral
  truth (asymmetric atom-cluster phantom via `pdb2mrc`, 400 particles
  simulated from it at SNR 0.1 with CTF, reconstructed at the exact
  simulated orientations with `objfun=cc`, `mskdiam=60`), both backends
  have FSC ≈ 1.000 to the truth and an amplitude ratio to the truth that is
  uniform in k (gridding 0.0159–0.0176, PCG 0.0156–0.0182 over k = 1–16);
  the residual backend difference at k ≤ 3 is ≈ −8 % (gridding) / +5 %
  (PCG) relative to their mid-band levels. The same run gives the cleanest
  confirmation of the deapodization fix: radial LS scale recon/truth with
  no high-pass, r = 0–24 px: gridding new envelope 1.01/1.00/0.99/0.99/
  0.98/0.97, legacy 1.01/1.00/0.97/0.94/0.92/0.90, PCG 0.99–1.02. Along
  the way two measurement pitfalls were found and fixed in the tool: the
  Fourier-shell comparison now removes each map's outside-mask background
  before masking (a constant solvent level masked by a sphere lands in
  shells 1–3), and `ml_reg=no` was shown NOT to change the low shells.
  Still open on real data: the pre-fix `rec3D_backends` in
  `stage2_tests` showed pcg/gridding 1.7–2.0 at k = 1–3 with the gridding
  map unmasked; to be re-measured with the current tool (masked, background
  removed) and the fixed build before drawing conclusions.
- The PCG solver's beyond-band behaviour under near-flat bootstrap sigma2
  (`pcg_euclid_crash_investigation.md` (consolidated into `pcg_priors.md`) §13, `PCG BEYOND-BAND EXCESS`
  diagnostic) is a separate solver-level question.
- Whether the 2D polar pipeline's own normalizations (`pftc` particle and
  reference preparation) place "expected signal" exactly at the data-quotient
  scale, or at a fixed multiple of it, should be confirmed with the §5.1
  instrumentation before reading the §3 ratios as absolute.

## 7. For later consideration: PCG convergence control in refinement

Decisions (2026-08-23): no back-compatibility for old (÷box) starting
volumes is required — the EUCLID DIAG warning is the only support.
Eight PCG iterations per refinement iteration (the stage-4 setting that
moved the in-band pcg/gridding ratio from 0.835 to 1.11) is not
affordable; `maxits_pcg=4` is being timed.

### 7.1 Warm-starting the solve from the previous iteration's half maps

The cheaper lever is to start each solve from the previous iteration's
map instead of `x = 0` (`simple_rec3D_pcg_strategy.f90`, the
`allocate(x, source=0.0)` branch of `reduce_solve_state_half`). The
within-iteration warm start (unregularized → ML solve) already exists;
the cross-iteration one does not.

**Box changes.** SIMPLE changes the degree of down-sampling only at
constant field of view (`box*smpd == box_crop*smpd_crop`, enforced at
the PCG entry). The previous map and the new one therefore sample the
same continuous density on different lattices, and the transform between
them is a Fourier-space zero-pad (box grows) or clip (box shrinks) with
**no amplitude factor**: under `fft = (1/N)Σ` the Fourier coefficients
are per-voxel means, invariant to resampling at constant FOV, and the
particle data are at the same scale for the same reason. This holds only
now that the ÷box convention is gone — under the old convention a box
change would additionally have required ×(box_new/box_old). The
operation exists as `image%read_and_crop` (`fft → pad_inplace/clip_inplace
→ ifft`, no scaling), already used by the PCG trailing bootstrap on the
previous `_even`/`_odd` half maps.

Grown shells start at zero, which is the cold start they would have had
anyway (the previous iteration carried no information beyond its Nyquist
and the FSC/ML prior is zero there). CG then spends its iterations where
the residual is — the new band and whatever moved — so the warm start is
never worse than cold and is best exactly where stage 4 found the
deficit (in-band amplitude).

Implementation rules:

1. Per half from its own half: even from previous `_even`, odd from
   `_odd`, never from the merged map (FSC independence).
2. Re-apply the support mask after padding: the solver solves inside
   `msk_crop`; Fourier padding rings slightly outside it, so zero the
   outside again (`msk_crop` scales with the box at constant FOV, same
   physical sphere).
3. The warm start must be at the data-quotient scale; a user-supplied
   first-iteration volume from an old build is caught by the EUCLID DIAG
   warning, no extra machinery.

Caveat: CG restarts its Krylov space at every solve, so a warm start buys
a smaller initial error, not continued convergence. The gain is largest
late in refinement (small orientation changes) and smallest at
abinitio3D stage handoffs. Expectation: with warm start, the residual
reached at `maxits_pcg=4` cold should be reachable at 2–3 iterations.

Acceptance: `RESID` at the final iteration and the in-band pcg/gridding
ratio from `test=rec3D_backends` at equal `maxits_pcg`, cold vs warm;
EUCLID DIAG and final resolution unchanged or better; FSC between
halves not inflated (compare half-map FSC cold vs warm at the same
iteration).

## 8. Post-release regression: stage-1 starting volumes were calibrated for the retired factor

Found 2026-08-23 from the first 10-run benchmark of the committed
convention change (gridding 8/10 correct maps, PCG(maxits_pcg=4) 9/10 —
both below the expected 10/10 on streptavidin). The EUCLID DIAG warning
fired on every stage-1 iteration of every run: low-band ref/ptcl
1.4-1.8E-3 (healthy: 0.25-0.55), v = 0.999, for ALL of stage 1
(iterations 1-20), snapping to healthy at the first stage-2 iteration.
Direct low-k amplitude measurement of the run volumes showed
`startvol_state01.mrc` at the SAME absolute scale in the old and new
builds (~5E-5 at k=1-6), while all reconstructed stage maps carry the
expected ×box in the new build.

Cause: three code sites fabricate or normalize starting volumes to an
ABSOLUTE scale that was calibrated for the era when `expand_cmat`
multiplied references by the original box:

1. `generate_random_volumes` (`simple_abinitio_utils.f90`) — noise
   startvol N(0, 5/box), with the comment "scaled to suit the
   euclid/sigma2 alignment scheme". Old effective reference scale:
   (5/box_crop)*box_orig. New code projected it as-is → references
   ~box_orig too small → stage-1 alignment effectively reference-free
   (v=0.999) → reduced abinitio3D reliability. (Stage 1 holds its
   references at the startvol scale throughout — pre-existing behavior,
   identical in the control run, where the DIAG ratio is flat 0.25 for
   iterations 1-20.)
2. `normalize_input_volumes` (`simple_abinitio_utils.f90`) — user-input
   volumes normalized to foreground stdev 1/box.
3. `prepare_external_init_vol` (`simple_commanders_refine3D.f90`) —
   e/o input volumes, same 1/box normalization (3 sites).

Fix (same day): all three now target the data-quotient reference scale —
the legacy value times the retired projection factor `params%box`:
noise b = 5*box_orig/box_crop; input normalization foreground stdev =
box_orig/ldim (=1 for native-box inputs).

Note the acceptance tests of §5.3a could not catch this: they exercised
reconstruction-produced volumes (which carry the new scale) and
externally ×box-scaled inputs, not the fabricated/normalized startvol
paths. The EUCLID DIAG warning caught it in the first real benchmark —
exactly the failure mode it was written for, though its text
("the first sigma2 update recovers") understated the impact: sigma2
adapts, but the reference amplitudes stay wrong for the whole lp-set
stage, degrading stage-1 alignment.

Validation required: rerun the 10x abinitio3D benchmark (both backends);
acceptance = no EUCLID DIAG warning, stage-1 low-band ref/ptcl ~0.25
(control level), 10/10 correct maps at control-level execution time.

## 9. Review response (2026-08-23)

Responses to the static branch review (all six findings addressed in code):

- **4.1 (P1, old maps collapse refinement)** — resolved by AUTOMATIC
  RESCALING at reference preparation (user's choice): a data-quotient
  reference has foreground sigma of order 1, an old-convention map sits
  ~box lower; `read_mask_filter_refvols` now measures the foreground
  sigma of every reference volume it reads and multiplies by the box
  size when sigma < 10/box (an order of magnitude of margin each way;
  idempotent by construction), logging loudly. External inputs were
  already normalized to the data-quotient scale by S8's fix. The EUCLID
  DIAG warning remains as a should-never-fire backstop and no longer
  claims that sigma2 recovery makes the condition benign.
- **4.2 (P2, Flex deapodization)** — `realize_hermitian_volume` now
  consumes its correction image (it was a dummy), and the coupled
  M-step deapodizes both half-volumes on the native lattice before FSC,
  Wiener merge, filtering and masking. Also fixed a padded-pointer
  whole-array energy sum (rmat sliced to ldim).
- **4.3 (P2, conventional-gridding helper)** — the helper now finalizes
  like production: no /BOX, native-lattice KB deapodization via
  `kb_stencil_inv_envelope_1d` + `deapodize3D_inplace`.
- **4.4 (P2, ungated test)** — `test=rec3D_backends` is now gated:
  agreement-band width >= 10 shells, >= 5 valid shells, median in-band
  amplitude ratio in [0.67, 1.5], median in-band FSC >= 0.9, >= 3 radial
  bins with normalised ratio range within [0.5, 2.0]; ground-truth mode
  additionally gates the gridding LS profile flatness in [0.92, 1.08]
  (legacy deapodization fades to ~0.90 and must fail) and median
  FSC(truth, gridding) >= 0.8. Violations print FAIL lines and end in
  THROW_HARD. Thresholds derive from the neutral-phantom fixture and the
  streptavidin reference runs recorded above; the loose in-band ratio
  bound intentionally passes the known PCG convergence dependence
  (0.835 at 2 iterations, 1.16 at 4) while failing any box-factor
  (x88) or deapodization mutation.
- **4.5 (P2, hot-path cost)** — new parameter `euclid_diag` (yes|no,
  default no; registered for refine3D, reconstruct3D, abinitio3D). With
  it off, `calc_sigma2` allocates only the sigma contribution and calls
  `gen_sigma_contrib` without the optional power/objective arguments, so
  the particle path executes the pre-branch calculation shape. Diag
  arrays are not allocated. Enable with `euclid_diag=yes` for transition
  validation runs.
- **4.6 (P3, part-1-only report)** — the report now appends
  "[PART 1/N ONLY]" when nparts > 1; aggregation across parts is
  deliberately not implemented.

Outstanding from the review's validation sequence: run the corrected
halfset/matrix production tests, the gated rec3D_backends fixture with
its two intentional-failure mutations, an old-map restart control
against a box-scaled control, the Flex fixtures, and the refine3D
recovery suite. Benchmarks should add `euclid_diag=yes`.

### 9.1 Gate verification (2026-08-23, local, neutral-phantom fixture)

- Base run: PASS. Gated band k=2-14 (agreement band capped at lp=10),
  median in-band pcg/gridding 0.935, radial and truth-flatness gates met.
- One gate-design correction was needed on the first run: the agreement
  band (backend-vs-backend FSC > 0.5) extends to Nyquist because the two
  backends correlate on shared noise far beyond the data band, where
  PCG's beyond-band excess (S6, tracked separately) dominates (median
  ratio 1.70 over k=2-64). The gates therefore act on the agreement band
  CAPPED AT THE DATA BAND (lp when given).
- Mutation 1 (restore /box on reconstructor output): FAIL as required —
  median gated-band ratio 119.7.
- Mutation 2 (omit gridding deapodization): FAIL as required — the
  gridding/truth LS profile hits 0.9 at radial bin 5, outside
  [0.92, 1.08], while the shell ratio stays 1.01 (radial gate carries
  the detection, as designed).
- Clean rebuild afterwards: PASS again; the mutations were reverted and
  the working tree verified clean of them.

## 10. Second benchmark round and the FSC/deapodization order fix

Benchmark round 2 (2026-08-23, commit 03844c75, euclid_diag=yes):

- Stage-1 startvol fix confirmed on real data: stage-1 ref/ptcl 0.25
  (gridding) / 0.21 (pcg) from iteration 1, matching the pre-refactor
  control; no DIAG warnings; no auto-rescales triggered.
- PCG (maxits_pcg=4): 5/5 correct maps, best 3.19 A, +23% wall time
  (1465 vs 1192 s); maps judged better than gridding by the user.
- Gridding: 7/8 correct, best 3.43 A — but the pre-refactor control
  reached 3.19 A: a gridding-only resolution regression.
- The gated rec3D_backends PASSes on real data (gated band k=2-45,
  in-band ratio 1.15).

Cause of the gridding regression: `sampl_dens_correct_eos` applied the
new deapodization to the half maps BEFORE computing the half-map FSC.
Legacy computed the FSC on uncorrected halves (the old grid correction
lived in the commander, after assembly). The inverse envelope's gain
rises toward the volume edge, so rim noise inside the FSC mask is
up-weighted (~1.4x at the mask radius for box 128), which pushes the
FSC crossing down ~3 shells (3.19 -> 3.43) and — because this FSC also
drives ml_reg and the reference filter — degrades the refinement itself.
PCG was unaffected (its envelope lives inside the solver).

Fix: the halves written to disk (final and _unfil) are deapodized
copies, while the in-memory halves used for the FSC/cFAR stay
undeapodized — restoring the legacy FSC estimate exactly. Verified: the
gated neutral fixture still PASSes. Validation: rerun the gridding
benchmark; acceptance is 3.19 A best resolution at the control failure
rate.

Aside from round 2: the PCG ml-replay RESID printed 0.72 without rtol
vs 0.31 with rtol=1e-3 at the same 4 iterations (STOP=fixed_iterations
vs maxits) — worth a look when tuning the §6 convergence-control item.

S10 validation (2026-08-23, user): with the FSC-order fix, gridding is
back at 3.191 A (RESOLUTION @ FSC=0.143 AVG/MIN/MAX 3.191) — the
pre-refactor control value. Reliability over a full batch to be
confirmed with the ongoing runs.

## 11. Opt-in PCG Laplacian (biharmonic) prior

Implemented 2026-08-23 as the first quadratic prior from the
"free, no algorithm change" family. R(x) = lambda ||nabla^2 x||^2 is a
Fourier multiplier ~ |k|^4, i.e. a diagonal on the same lattice as
Khat/ml_prior, so it enters BOTH the kernel operator and the
preconditioner (same padsc^2 conventions as the ML prior; disabled
during kernel-scale calibration like the other regularizers; kernel op
mode required, matrix-free throws).

Parameterization: `pcg_lambda_lap` (default 0 = off; UI: reconstruct3D,
refine3D, test=rec3D_backends). The added diagonal is
lambda_lap * data_scale * (k/k_nyq_native)^4 — lambda_lap is the
prior-to-data ratio at the native Nyquist shell; at nu = 0.5 the prior
is 6% of nominal, so the data band is essentially untouched by design.

First measurement (neutral fixture, objfun=cc, DEFAULT maxits_pcg=2,
lambda_lap = 1.0 vs off):

- Beyond-band excess ELIMINATED: pcg amplitude at k=44..64 falls from a
  rising noise floor (1.6E-3..2.5E-3, tracking 1.7-2.7x gridding) to a
  smooth roll-off (2.0E-4..3.1E-5, 0.23..0.03x gridding) — 10-80x
  suppression, while shells k<=16 change < 2%.
- UNEXPECTED WIN: the in-band deficit at 2 iterations largely
  disappears — gated-band median pcg/gridding 0.840 (off) -> 0.9986
  (on). The prior in the preconditioner improves conditioning enough
  that the 2-iteration solve is nearly converged in band; this is the
  same deficit that previously needed maxits_pcg=4-8.
- All rec3D_backends gates PASS with the prior on and off.

Suggested real-data protocol (reconstruct3D on streptavidin, opt-in):
sweep pcg_lambda_lap = 0.3 / 1.0 / 3.0 at maxits_pcg = 2 AND 4; read
the shell table beyond-band rows, the gated-band ratio, RESID, and
FSC=0.143. If in-band convergence at 2 iterations holds on real data,
the prior may remove the need for maxits_pcg=4 entirely (its +23%
wall-time cost).

Note for the test tool: pcg-only keys (pcg_lambda_lap, pcg_lambda_rel)
are deleted from the gridding leg's command line — they hard-error under
rec_backend/=pcg by parameter validation.

### 11.1 First real-data Laplacian reading (bgal, box 256) and gate refinements

bgal (d2, mskdiam 180, maxits_pcg=2, pcg_lambda_lap=0.3, no lp):

- The prior behaves as the Wiener term it is: attenuation rho/(rho+
  lambda*nu^4*S) — ratio ~1.0 to k~35, crossing the falling shell SNR at
  k~70-90 (-25% at 5.4 A, -50% at 4.2 A), beyond-band tail crushed
  (0.03-0.10 at k>=105). On the fixture the crossing sat beyond the
  band edge; on real data 0.3 puts it mid-band. Sweep DOWN:
  0.01/0.03/0.1; target = crossing at/past the FSC band edge. Pass lp=
  so the gate acts in-band.
- bgal exposes the PCG centre-bin deficit dramatically (0.31 vs
  streptavidin's 0.83 at 2 iterations; k=1-2 excess 2.0/1.6) — the S6
  low-shell item scales with box. Gate changes: the radial-flatness
  ratio is now normalised to the IN-MASK MEDIAN (normalising to the
  deficient centre bin inflated all bins ~3x -> false FAIL), the centre
  bin is excluded from the range and printed as its own diagnostic
  line, and the range gate applies to bins 2+.
- Post-mortem of a self-inflicted test bug: while making these edits
  the mskrad computation was moved after the ground-truth block that
  uses it, so truth comparisons and backend spectra ran with an
  uninitialised mask radius (truth FSC appeared to collapse to
  0.87-0.97; reconstruction was fine). Fixed by restoring the early
  computation; the fixture PASSes with the same 0.84 in-band ratio as
  the pre-edit baseline. Lesson recorded: any reordering inside the
  test must respect mskrad/nrb_msk being consumed by the truth block.

### 11.2 bgal sweep verdict: the Laplacian is a convergence aid, not a spectral tool

Sweep (box 256, maxits_pcg=2, lp=3.6, lambda_lap = 0 / 0.01 / 0.03 / 0.1,
plus 0.1 at maxits 4). Key numbers:

- ML-replay RESID at 2 its: 6.53 (lambda=0!) -> 1.88 (0.01) -> 1.13
  (0.03) -> 0.86 (0.1); at 4 its with 0.1: 0.25. Base solves: 0.20 at
  2 its, 0.048 at 4. The production output IS the ML replay, so at
  lambda=0 and 2 iterations bgal ships a barely-updated warm start.
- The centre-bin deficit tracks ML convergence, not lambda: 0.73
  (lambda=0, map ~ base warm start), 0.34 (lambda>0, transitional),
  0.50 (4 its). The S6 "centre/low-shell" item IS under-convergence of
  the ML replay, amplified by box size (streptavidin box 128: 0.83).
- k=1-2 excess likewise: worsens at 2 its with lambda (3.0/2.1),
  recovers at 4 its (1.6/1.6).
- In-band spectral cost: even lambda=0.01 attenuates k>=80 (0.58 vs
  0.77 at k=80 for lambda=0), and at 4 its/0.1 the converged
  over-regularized solution emerges (0.30 at k=85). The nu^4 prior
  anchored to the low-band data scale cannot separate "in band" from
  "beyond band": rho falls ~3 orders across the band while nu^4 spans
  ~1.5 from half-band to edge — the usable lambda window is empty on
  real data. Band-selective suppression is the ML/FSC prior's job.

Verdict: keep pcg_lambda_lap opt-in as a solver/conditioning tool; do
NOT adopt for production spectral shaping. The production PCG path
forward is ML-replay convergence: cross-iteration warm start from the
PREVIOUS ML map (S7) so the solve starts near the regularized solution
instead of the unregularized base (the RESID 6.5 gap), plus a
box-scaled iteration budget. Discriminating experiment: maxits_pcg=4
and 8 at lambda=0 on bgal — if RESID and the centre bin recover with
iterations alone, convergence control closes both S6 items without any
new prior.

### 11.3 Both experimental lambda paths removed (2026-08-24)

Decision (user): too many paths — `pcg_lambda_lap` and the user-facing
`pcg_lambda_rel` are REMOVED. What they established is recorded here so
the experiments need not be repeated:

- `pcg_lambda_lap` (biharmonic smoothness, S11-11.2): a nu^4 prior
  anchored to the data scale cannot separate in-band from beyond-band
  on real data (rho falls ~3 orders across the band, nu^4 spans ~1.5);
  its real benefit was CONDITIONING the ML replay (RESID 6.5 -> 0.86 at
  2 its on bgal). That benefit should be captured instead by the S7
  warm start from the previous ML map and a box-scaled iteration
  budget, not by a spectral prior.
- `pcg_lambda_rel` (relative Tikhonov CLI): never used in production.
  The INTERNAL mechanism `reconstructor_pcg%set_lambda_relative`
  remains (used by the pcg mass/reduction test commanders with
  MASS_LAMBDA_REL); only the CLI parameter, its parsing/validation, the
  strategy and abinitio3D forwarding, and the UI registrations were
  removed. The plain absolute lambda (PCG_LAMBDA) is unchanged.

Removal verified: full sweep clean (only the internal setter and the
'pcg_lambda_effective' diagnostic label remain), build clean, gated
neutral fixture PASS at the baseline in-band ratio (0.843).

## 12. Failed-run forensics (streptavidin, ref_taper batch)

Two failed gridding runs (RESTART5/6) analyzed against healthy references:

- Both ran with ref_taper=yes (stage-1 ref/ptcl 0.127/0.167 = 0.25 x the
  volume-averaged KB envelope — the taper's measurement footprint; DIAG
  levels are not comparable to untapered runs).
- The taper reached ALL alignment paths (matcher AND the prob_tab/
  prob_align reprojection model, which is materialized through
  read_mask_filter_reproject_refvols) — and did not prevent these
  failures.
- Failure mode: a self-consistent WRONG C1 fold locked in before/around
  symmetrization; the D2 axis itself is essentially correct (C1-vs-
  symmetrized correlation 0.85/0.86 vs 0.88-0.93 in healthy runs); runs
  then refine to plausible statistics (3.4-3.6 A, v 0.963) on the wrong
  structure. No warnings, no scale anomalies, no late-stage collapse.
- This is the classic stochastic ab initio failure (bad early
  trajectory), not a detectable defect of the new convention.

MISSING CONTROL: the true master baseline failure rate on the same data
was never measured in this investigation (the "should never fail" prior
is experiential). Current evidence: 2/10 + 1/9 + 1/10 across three
datasets on the branch. Decisive experiment: 10x abinitio3D on MASTER
(same data, box, seed policy) — if master also drops ~1/10, the branch
is exonerated on reliability and the remaining delta is the cosmetic
periphery brightness of correctly-deapodized maps; if master is 0/10,
the residual suspect list is the early-stage trajectory sensitivity
(stage 1-3) and needs a successful-vs-failed trajectory comparison
within one batch.

### 12.1 Reliability question CLOSED (2026-08-25)

The missing control was run: MASTER gives 1/10 failures on the same
streptavidin data. The branch rates (2/10, 1/9, 1/10) are statistically
indistinguishable from the intrinsic stochastic ab initio failure rate
— the branch has NO reliability regression. Consequences:

- `ref_taper` is REMOVED (hypothesis dead: the legacy KB reference
  taper is not what kept master reliable; master fails at the same rate
  without any of the branch changes). Findings preserved in S12.
- The "overfitted look" of exports stands as what S10-12 measured it to
  be: correctly-deapodized flat density showing rim/periphery content
  the legacy KB fade used to suppress at any given display threshold —
  cosmetic, with FSC, DIAG, and failure rate all at baseline. If the
  display preference matters, it belongs in postprocessing policy, not
  in the density convention.
- Remaining open items are PCG-only: the S7 ML-replay warm start
  (RESID 6.5 at 2 its on box 256 is the real S6 root cause) and a
  box-scaled iteration budget.

## 13. Cross-iteration ML warm start (branch pcg_ml_warm_start)

Implemented 2026-08-25 per S7/S11.2: both ML-replay solve sites in
simple_rec3D_pcg_strategy.f90 (shared and distributed) now warm-start
from the PREVIOUS refinement iteration's ML half map when one exists on
disk (params%vols(state) + '_even'/'_odd'), via read_and_crop
(factor-free constant-FOV pad/clip under the data-quotient convention)
with the soft support mask re-applied after resampling. Each half uses
strictly its own previous half (FSC independence). The first-iteration
noise starting volume is excluded by its 'startvol' workflow-contract
name; with no usable previous half the solve keeps the base-solution
warm start (previous behavior, e.g. reconstruct3D standalone — so
test=rec3D_backends is unaffected by design and cannot exercise this).

Expected effect: closes the warm-start gap measured on bgal (ML RESID
6.5 at 2 iterations from the base solution) because the previous ML map
already carries the regularized beyond-band roll-off; the centre-bin
and k=1-2 anomalies (S6) should shrink with it.

Acceptance: a refine3D continuation with rec_backend=pcg maxits_pcg=2 on
bgal — '>>> PCG ML WARM START' lines present from iteration 1, ML RESID
per iteration well below the cold 6.5 and decreasing, centre-bin
diagnostic (post-run rec3D_backends) toward 1, resolution not degraded
vs maxits_pcg=4 cold. Compare against master (no warm start) at
identical settings.
