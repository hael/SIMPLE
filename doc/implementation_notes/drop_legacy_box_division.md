# Dropping the legacy box division from reconstruction output

**Status:** in progress on branch `drop_legacy_box_division` (started
2026-08-22). Step 1 (instrumentation, §5.1) is implemented; steps 2–6 are
not started. Until step 3 lands, the PCG backend mirrors the convention
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
(`pcg_euclid_crash_investigation.md` §15): gridding stage-1 lp map L2 0.99 vs
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
scratchpad, recorded in `pcg_euclid_crash_investigation.md` §15.4): L2 ratio
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
3. **Remove the division and fix the gridding deapodization.**
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
6. **Compatibility.** Release note on the ×box change of map values; a
   warning when a starting reference reprojects ≪ the particle scale.

## 6. Open items it does not settle

- The PCG solver's beyond-band behaviour under near-flat bootstrap sigma2
  (`pcg_euclid_crash_investigation.md` §13, `PCG BEYOND-BAND EXCESS`
  diagnostic) is a separate solver-level question.
- Whether the 2D polar pipeline's own normalizations (`pftc` particle and
  reference preparation) place "expected signal" exactly at the data-quotient
  scale, or at a fixed multiple of it, should be confirmed with the §5.1
  instrumentation before reading the §3 ratios as absolute.
