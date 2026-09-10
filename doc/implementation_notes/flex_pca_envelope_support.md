# flex_pca: sampling-distance control, reconstruction backends, and the envelope mask as solve support

## Status

**Planning only (2026-09-10, revised after review). No code is proposed for immediate implementation.**

Scope, as decided:

- `flex_pca` contains no masking logic. The mask is an input: a real-space
  `[0,1]` volume with a smooth falloff, supplied through `pcg_mskfile`. The
  program resamples and validates it, nothing more.
- Before anything else, the fixed `box_crop=64` working box is replaced by a
  sampling-distance parameter that controls down-sampling, in the way
  `refine3D` and `abinitio3D` do it.
- Then `rec_backend=gridding` (the current implementation) and
  `rec_backend=pcg` (new) are supported in the same pattern as those
  workflows, reusing `reconstructor_pcg` wherever it applies, first with the
  spherical support that is used today.
- Only then does `pcg_mskfile` give the envelope-constrained analysis, and it
  is available on the PCG backend only, as it is everywhere else
  (`parameters_phases.f90:897`).

Code references are to `0eee68726`.

## 1. Order of work

```text
step 0   smpd_target-driven working box (replaces box_crop=64)
step 1   rec_backend=gridding|pcg in flex_pca, spherical support, reconstructor_pcg reused
step 2   pcg_mskfile -> soft envelope as the solve support (PCG only)
```

Each step is gated on its own and none needs the next; step 0 is a
prerequisite for both others because the box decides the memory and the
resolution of everything below it.

## 2. Step 0: a sampling distance instead of a fixed box

### 2.1 Today

`exec_flex_pca` sets `box_crop=64` unless given
(`simple_commanders_flex_pca.f90:43`), then `derive_flex_pca_band` derives
`lp = 2.5 * smpd_crop` (`COV_LP_OVER_NYQUIST`) and `box_rec = box`
(lines 133-170); `parameters` derives
`smpd_crop = box/box_crop * smpd` (`parameters_phases.f90:326`). The
particle cache (`ptcl_cache_ensure`), the projected-model image prep, the
basis reconstructors, the `vol1` mean and the probe-basis volumes
(`read_and_crop` at `box_crop`/`smpd_crop`) all follow `box_crop`. The
covariance is therefore fitted at whatever sampling `box/64` happens to be:
~3.5 A on a 200 A particle in a 256 box at 0.9 A, ~6 A in a 400 box -- the
resolution of the analysis is an accident of the box.

### 2.2 Proposed

Follow `refine3D` (`simple_commanders_refine3D.f90:187-215`): a target
sampling distance, `autoscale(box, smpd, smpd_target, box_crop, smpd_crop,
scale, minbox)` from `simple_magic_boxes.f90:104` (magic box sizes, never
finer than the data, floor at 64), and `box_crop`/`smpd_crop` set on the
command line for the workers. Points specific to `flex_pca`:

- **Default.** Helix-level detail needs the working Nyquist at or below
  ~4.5 A (alpha-helical pitch 5.4 A, strand separation 4.8 A), i.e.
  `smpd_target <= 2.25 A`. A default of **2.2 A** is proposed; the exact
  value is an open question (§6), and `refine3D`'s 1.3 A is too fine for
  the covariance stage at present accumulator sizes.
- **`lp` and `smpd_target` are one decision.** Today `lp` is derived from
  the box; with a sampling distance both directions are possible and only
  one should be primary. Proposed: `smpd_target` is primary and
  `lp = COV_LP_OVER_NYQUIST * smpd_crop` stays derived; an explicit `lp`
  instead sets `smpd_target = lp / COV_LP_OVER_NYQUIST` before
  `autoscale`, so a user asking for a 6 A band gets a 2.4 A lattice rather
  than a 6 A band on a 3.5 A lattice. An explicit `box_crop` remains an
  override for tests.
- **`box_rec` stays the native box.** The delivered state maps are not
  capped at the covariance Nyquist; unchanged.
- **Memory scales as `box_crop^3`.** The coupled M-step accumulators
  (`rho_cross_exp`, `npairs x` expanded grid) and, on the PCG backend, the
  pair kernels (§3.3) grow with the cube of the box. Going from 64 to 96
  or 128 is x3.4 or x8. The existing `COV_ATHR_BUDGET` /
  `cov_dim_budget` machinery budgets the reduced solve's dimension; the
  same budgeting has to cover the accumulators, and the log should print
  the footprint at start-up.
- **Provenance.** The embedding cache (`infile` resume) is tied to the
  particle selection; it must also be tied to `box_crop`/`smpd_crop`, or a
  cache built at one sampling will be resumed at another. Same for the
  probe-basis files and `flex_pca_probe.txt`.

## 3. Step 1: `rec_backend` in flex_pca

### 3.1 Which reconstructions

`flex_pca` reconstructs in two places (the mean is an input, `vol1`;
`COV_MEAN_FROM_DATA=.false.`):

1. **State maps** (`simple_flex_pca_rec3D.f90`): kernel-weighted
   reconstruction per state and half at `box_rec`, currently gridding with
   a shell-relative density floor (`floor_rho`), `zero_background` and
   `mask3D_soft` at the broadest sphere in the box (lines 495-502). This is
   the direct analogue of `reconstruct3D`.
2. **The M-step basis solve** (`simple_flex_pca_em_iter.f90:2394-2432`):
   `solve_coupled_basis_exp` (`simple_flex_projected_latent_model.f90:446`)
   divides the `k` right-hand sides per Fourier voxel by the `k x k`
   gridding density `rho_cross_exp`, then inverse-KB envelope, half-set
   FSC Wiener `2F/(1+F)`, band-limit, `mask3D_soft(msk_crop)`.

Both get the backend switch; the gridding paths stay as they are.

### 3.2 Pattern

`rec_backend` is already a validated global parameter
(`parameters_phases.f90:779-782`), `pcg_mskfile` already requires it to be
`pcg` (line 897), and the PCG strategy's preconditions apply unchanged
(`pcgop=kernel`, `mskdiam` set, `projrec=no`;
`simple_rec3D_pcg_strategy.f90:613-632`). `flex_pca` has its own
master/worker stage protocol (`simple_flex_pca_distr.f90`,
`PCA_STAGE_STATES`), so the selector lives in the flex layer rather than in
`create_rec3D_strategy`: one switch on `params%rec_backend` at the two
sites above, and the `>>> ... REC3D EXECUTION (gridding|kernel PCG)` style
log line so a run says which backend it used.

### 3.3 States on `reconstructor_pcg`

Reused wholesale; the kernel weights enter through the noise model, not
through a new weight array.

**Mechanism.** In `prepare_fused_planes`
(`simple_reconstructor_pcg.f90:2158-2222`) each particle contributes the
right-hand-side plane `weighted(h,k) = conj(CTF.shift) y(h,k) / sig2(shell,i)`
and the density plane `absT2(h,k) = CTF^2 / sig2(shell,i)`; `absT2_plane`
(2087-2145) builds the same density for the kernel path. Everything else in
the solve -- the preconditioner `1/(rho+floor)`, the Gram kernel `Khat`, the
data scale behind a relative `lambda` and the `P_tau` scale, and the raw
accumulator artifacts the workers write -- derives from those two
accumulators (`finalize_density_accum`, `finalize_khat`,
`update_lambda_from_density`). The weighted least-squares objective
`sum_i w_i ||Sigma_i^{-1/2}(A_i x - y_i)||^2` is therefore obtained exactly
by passing `prep_particles` the array
`sig2(0:R,i) = sigma2_noise(:,pinds(i)) / w_i(s)`: the weight multiplies the
RHS term and `|T_i|^2` identically, as the normal equations require, and
propagates to the preconditioner, the kernel, the ridge and the raw
artifacts without touching any of them. The rec3D PCG strategy already
builds that array per particle from the sigma store
(`simple_rec3D_pcg_strategy.f90:712-718`; under `cc` it is a unit array and
the division applies just the same), so the change is confined to the flex
caller: one division per particle and shell before `prep_particles`, per
state and per half, the even/odd split coming from the disjoint halves of
the weight table as `mask_state_weights_by_half` does today.

**Details.**

- Particles with `w_i(s) = 0` are removed from `selection`/`pinds` for
  that state rather than given an infinite sigma; they contribute nothing
  and dropping them shrinks the particle pass, which matters for a state
  holding a small fraction of the data. A floor on `w` (of order 1e-3)
  keeps `sig2` finite in single precision; particles below it are dropped
  too.
- The ridge is set in relative mode (`set_lambda_relative`), not the
  absolute `PCG_LAMBDA` the refinement strategy passes to `new`: with
  weights in `[0,1]` the effective particle count of a sparse state can be
  a small fraction of `N`, and an absolute 1e-3 would be a materially
  stronger prior on low-occupancy states than on populated ones, whereas
  the relative form scales with the weighted `D` and treats every state
  alike.
- The shell-relative preconditioner floor `RHO_FLOOR_FRAC` scales with the
  weights automatically; it takes over the role of flex's `floor_rho`
  without biasing the estimate.
- Rejected alternative: an explicit per-particle weight array inside the
  reconstructor, multiplied in at the two plane routines. Same numbers,
  more code, a second scaling path beside sigma; `add_raw_accum_weighted`
  is no help since it weights whole artifacts (trailing-reconstruction
  continuation), not particles. It would only be worth it if the solver
  sidecars had to report unweighted sigma statistics, and the data scale
  of the weighted problem is the one that describes the solve.
- The sigma route works for the states because the weight is a scalar per
  particle. It does not carry over to the M-step, where the per-particle
  weight is the `k x k` matrix `E[z_i z_i']` and the pair-weighted kernels
  need the explicit generalization of `accumulate_absT2` (§3.4).

**Flow.** Per state and half: accumulate at `box_rec`, workers
`write_raw_accum` (keyed by state and half), master `add_raw_accum`, base
solve cold with at least `FINAL_PCG_MAXITS_FLOOR` iterations, spherical
support via `set_mask` at the `msk` radius, shipped map `window * u`,
provenance sidecar. This replaces `floor_rho`, `zero_background` and the
broadest-sphere mask. The merge gate
(`simple_flex_pca_merge.f90:392-417`) then compares support-constrained
states; its 0.98 similarity threshold was calibrated on sphere-masked
gridding maps and is re-calibrated on the synthetic fixture.

### 3.4 The M-step on the PCG operator

The coupled normal system has `(q,r)` blocks
`H_qr = sum_i M_i(q,r) a_i^2 P_i^H C_i^2 P_i / sigma^2`, with
`M_i = E[z_i z_i']`: the single-volume operator with a scalar weight per
particle and pair. The present per-voxel `k x k` solve on the gridding
density is the block-Jacobi preconditioner of that system, so the PCG
version keeps `solve_coupled_basis_exp` in exactly that role and adds:

- pair-weighted Gram kernels `Khat_qr`, `k(k+1)/2` of them, from the same
  particle pass that accumulates the right-hand sides (the weighted
  generalization of `accumulate_absT2`); real and symmetric, packed
  half-grid, ~4 MB each at box 64, x3.4 at 96, x8 at 128 -- budgeted with
  the pair count (§2.2);
- the support `P (x) I` on the hard domain `window > 0`, deapodization
  inside the operator (replacing the post-solve `mstep_gridcorr`
  multiply), shipped basis `window * U`;
- CG to a relative-residual tolerance, not a fixed count: at these box
  sizes an iteration is `2 x npairs` small FFTs and costs nothing next to
  the accumulation, and the EM fixed point should be the least-squares
  fixed point rather than an iteration-count artefact.

Even and odd bases are solved separately on the same support, as now, and
the per-component half-set FSC Wiener filter and band-limit stay
post-solve in the first version, followed by a re-window (the filter leaks
a little across the support edge). Moving the per-component FSC precision
inside the operator as `P_tau` is a possible second version, not part of
step 1. The marginal-likelihood diagnostic `-2logL/N` logged per probe
iteration (`em_iter.f90:2367-2393`) is the arbiter between the backends,
since a half-set FSC under a common window is optimistic on both.

### 3.5 Unchanged

E-step, embedding, kernel weights, targets and the merge logic. The
initial basis (`simple_flex_pca_em_basis.f90:88`) and the probe-basis load
(`simple_flex_pca_em_fit.f90:505`) multiply by the same window the solve
uses -- with the spherical support that is the sphere they use today, so
step 1 changes nothing there.

## 4. Step 2: `pcg_mskfile` as the support

### 4.1 Input contract

A real-space `[0,1]` volume with a smooth falloff at the project's native
box and sampling (`automask3D_stateNN.mrc` from a `refine3D_auto`
`automsk=yes` run, or any user mask). `flex_pca` does not threshold,
dilate, soften or otherwise edit it. It validates: cubic, box equal to the
project box, values in `[0,1]` after resampling, non-empty. Whether the
envelope is generous enough to contain the moving parts is the user's
responsibility when preparing the mask; the dilation and edge controls
live in `automask3D` and the refinement policies, not here.

### 4.2 Resampling without Fourier-crop artefacts

`read_and_crop` (`simple_image_geom.f90:673`) is an `fft`,
`clip_inplace`, `ifft`: a brick-wall low-pass. On a mask that is
harmless when the falloff is wide compared with the target voxel, and
ringing when it is not: overshoot above 1, undershoot below 0, and --
the case that matters for a hard support -- positive ripples *outside* the
mask that would extend the domain `window > 0` into solvent. The mask is
needed on two lattices, `box_crop`/`smpd_crop` for the basis and
`box_rec`/`smpd_rec` for the states, and the same routine serves both:

1. Fourier crop to the target box (`read_and_crop` as is).
2. Record `min` and `max` before clamping and the fraction of voxels
   outside `[0,1]`; these are the ringing diagnostic and go to the log.
3. Clamp to `[0,1]`.
4. Floor: values below a small threshold are set to exactly 0 so the
   hard domain is the intended envelope and not its ripple; the threshold
   is a constant of the mask contract (the PCG solver already treats
   `window < PCG_SUPPORT_DIV_MIN = 0.1` as outside for its warm-start
   division, `simple_reconstructor_pcg.f90:597-606`; the floor here should
   be no larger than that and is a decision to make on the fixture).
5. Report occupancy (fraction of the box, and of the `msk` sphere) on the
   native and the target lattice; a large change is the sign that the
   falloff was too narrow for the target sampling.

The falloff width that keeps the Fourier route clean is a few target
voxels; at 2.2 A that is ~7-10 A, which the 7.5 A `ENVMSKWIDTH_A_MIN`
dilation with its cosine skirt already satisfies for masks that come from
`automask3D`. If the diagnostics of step 2 show ringing on user masks in
practice, the fallback is real-space (trilinear) down-sampling of the mask,
which cannot ring and only blurs the edge slightly; it should not be the
default because the Fourier route is what every other volume in the run
goes through.

### 4.3 Where the support enters

Exactly where the sphere enters in step 1, with `set_mask_volume` in place
of `set_mask`: the state solves at `box_rec`, the coupled M-step at
`box_crop`, the initial-basis and probe-basis windows, and the merge gate
(which then compares windowed states and needs no mask of its own). The
`vol1` mean from a PCG `automsk=yes` refinement is already windowed by the
density envelope; passing that same envelope file as `pcg_mskfile` makes
mean and basis share one support, which is the recommended use. On the
gridding backend `pcg_mskfile` is rejected by the existing validation, so
the envelope-constrained analysis is a PCG-only mode, as in `refine3D`.

## 5. Validation

- **Step 0:** the same run at `smpd_target` 3.5 / 2.6 / 2.2 A on the
  synthetic fixture (`production/tests/simple_test_flex_pca.f90`) and on
  one real dataset; ground-truth basis capture, `-2logL/N`, wall time and
  peak RSS per setting. The default is the coarsest setting on the capture
  plateau.
- **Step 1:** gridding vs PCG at the same sampling and spherical support:
  per-state unmasked FSC, basis capture, `-2logL/N`, reproducible
  dimension count and the principal-angle history (the 0.97 / 0.02
  stopping constants were tuned on gridding bases and may move), the
  merge-gate threshold, and a `refine3D_states` run seeded from each set of
  states.
- **Step 2:** the resampling diagnostics of §4.2 on `automask3D` masks and
  on deliberately sharp masks; sphere vs envelope support on the fixture
  and on a membrane protein (PfCRT), watching the leading-axis background
  ratio that the `DEFLATE_BG` comment records at 27x on real data and the
  micelle's share of the leading components.

## 6. Open questions

- Default `smpd_target`: 2.2 A is proposed on the helix argument; the
  fixture sweep may argue for a coarser default with the fine value as an
  option.
- Is `smpd_target` primary and `lp` derived (proposed), or the reverse?
- The floor below which a resampled mask value counts as outside the
  hard support (§4.2 item 4).
- Hard support (`window > 0`) at 2.2 A sampling, or the soft `P = window`
  alternative that `install_support` rejects for refinement? The
  refinement argument (a solver-state-dependent mixture of `P*u` and `u`)
  applies here too, so hard is the default proposal.
- Whether the per-component Wiener step ever moves inside the operator
  as `P_tau` (a second version of §3.4, not part of the plan).
