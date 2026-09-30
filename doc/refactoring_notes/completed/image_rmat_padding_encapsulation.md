# The padded real-space array of `image`: zero padding, and a box-shaped pointer

Date: 2026-09-25

Status: steps 0, 1 and 2 implemented and tested on Linux, 2026-09-30, on master `f9bbf72ce`
(uncommitted; section 8 is the record). Open: the fast gate on macOS Debug with bounds checking
(section 7). Step 3, removing the pointer, is out of scope (Hans, 2026-09-25). Reviewed 2026-09-30
against the code before the work: the survey holds (158 calls, 29 files, 63 statements); the review
added the routines that enter real space without a transform, the whole-array writes through names
other than `self`, the reductions inside the image class, the padding query for the tester and the
bounds-checking note of step 2. The implementation recounted the survey and corrected sections 1,
3, 4 and 5 where a count or a classification was wrong; section 8.7 lists the corrections.

Validation level: sections 1 to 7 were written from static inspection of the image class and of
every caller of `get_rmat_ptr`. Section 8 is measured: Release and Debug (bounds checking) builds
with GCC 15.2 and FFTW 3.3.5 on Oracle Linux 8.10, the fast gate, the library tier, four highlevel
gates and a real-data `flex_pca` run.

## 1. The problem

An image is one FFTW buffer seen two ways (`simple_image_core.f90`, `new`): `cmat(fdim(n1), n2, n3)`
in Fourier space and `rmat(2*fdim(n1), n2, n3)` in real space, where `fdim(n1) = n1/2 + 1`. The real
array therefore has padding rows `n1+1 .. 2*fdim(n1)` in its first dimension, two for an even box and
one for an odd box, which the in-place transforms need and which are not part of the image.

`get_rmat_ptr` (flagged `! VIOLATES ENCAPSULATION` in `simple_image.f90`) hands out the whole buffer.
It is called 158 times from 29 files outside the image submodules; nothing outside uses a `cmat`
pointer. Of the pointer uses, about 260 index or slice it within `ldim`, which is safe. About 120 use it
as a whole array or pass it on, and those fail in two ways:

- **Non-conforming expressions** against a box-sized array. Invalid Fortran; bounds checking stops on
  it, and without bounds checking it can pass by luck. The `cif2mrc` tester did this and passed on
  macOS, whose Debug build had no bounds checking (restored 2026-09-25); policy section 4.6 records the
  trap.
- **Reductions and updates over the padded array itself.** These conform, so no checker sees them,
  and they are right only if the padding holds zeros. Nothing guarantees that:
  - FFTW leaves the padding of an in-place complex-to-real output undefined, and `bwd_ft`
    (`simple_image_fft.f90`) does not clear it; nor does `ifft_mask_pad_fft`, which also leaves the
    image in real space.
  - Routines that enter real space by filling the box leave in the padding whatever the buffer held.
    `read`, `read_raw_mrc` and `read_single_mrc_image` (`simple_image_io.f90`) fill only the box, so
    in a read, `fft`, read loop the padding holds Fourier coefficients of the previous image. `ran`,
    `gauran` and the `gauimg` family write the box and set `ft = .false.` without clearing the rest
    (`corners`, `square` and `window` zero the whole buffer first). The image submodules have 33
    `ft = .false.` statements in all.
  - The image class writes into the padding in real space. Of the 63 whole-array `self%rmat = ...`
    statements in the image submodules, the scalar broadcasts put values there: assigning a scalar
    (`img = 1.`, `rmat = realin`), adding or subtracting a constant, `norm` (subtract the mean, divide
    by the standard deviation), and dividing one image by another (`self1%rmat/self2%rmat`: 0/0, a
    NaN, in zero padding). So do maps and masked assignments that send a zero to a nonzero value:
    `bin_inv` (0 becomes 1), `binarize_1` with a threshold at or below zero, `ring` with an inner
    radius at or below zero.
  - The count of 63 is of statements on `self`. Another 28 whole-array writes go through other names
    (`self_out` 14, `self_in` 2, `img_pad` 2, `thread_safe_tmp_imgs` 2, and one each of `img_out`,
    `tmp`, `pc`, `response_freq`, `img_filt`, `self_prev`, `even_prev`, `odd_prev`; recounted
    2026-09-30, the review had 26). Of those, `pad` (`self_out%rmat = backgr`) and `collage`
    (`img_out%rmat = background`) end in real space with a value in the padding. The scalar fills of
    the `*_pad_fft` family are followed by a forward transform, which ignores the padding.
  - Neither count includes a section that is whole in the first dimension only. `cendist` adds the
    squared distances with `self%rmat(:,i,1) = self%rmat(:,i,1) + ...`, so the padding gets a partial
    distance (`disc` already bounds its mask for that reason; `ring` did not, and its `npix` counted
    padding pixels: 359 for 317 at a radius of 10 in a box of 64).
  - Five reductions inside the image class read the padded array: the foreground count of
    `binarize_1` (`count(self_in%rmat > thres)`, wrong for a negative threshold even when the padding
    is zero), `npix` of `ring` (wrong when the padding was set to 1, above), the range checks of
    `remove_edge` and `one_at_edge`, `contains_nans` and `checkimg4nans` (`check4nans3D(self%rmat)`;
    all but the first are right only when the padding is zero).
  - The outside reductions that matter are in flex PCA: the Gram matrix of `em_basis`
    (`sum(real(rmat_i,dp)*real(rmat_j,dp))` over volumes straight out of an `ifft`, a gridding
    correction and a window), the norms, projections and deflations of `em_mstep` and
    `em_pairmerge`, the Gram and Gram-Schmidt of `em_solve`, the cosines of `em_crossfsc`.
    `realize_hermitian_volume` slices its energy to `ldim`; the others do not. Whether their results
    are affected depends on what the padding holds, which step 0 measures.

## 2. The contract

- In real space, the padding of an image holds zeros, and the image class keeps it so (step 1).
- Outside the image class, a real-space pointer is the box: `get_rmat_ptr` returns
  `self%rmat(:ldim(1),:ldim(2),:ldim(3))` (step 2). The whole buffer is available only through an
  explicitly named accessor, for the few callers that need the FFT layout.
- Inside the image class, a real-space operation that would write a nonzero value into the padding
  acts on the box. The test is the value at zero: an elemental map `f` with `f(0) /= 0`, a scalar
  assignment, or a masked assignment whose mask is true at zero. A reduction (`count`, `any`,
  `minval`, `maxval`) acts on the box as well, since zero is not neutral for it.
- After step 2 the first rule protects only code inside the image class and the caller of the padded
  accessor; outside, the box pointer makes every reduction a box reduction whatever the padding holds.

## 3. Step 0: measure first

Before any change, so that what steps 1 and 2 change is known:

- A type-bound query on `image`, for example `max_abs_padding()`, returning the largest magnitude in
  the padding (NaN if any value there is NaN). The tester uses it, so it needs no pointer to the
  padded buffer, before or after step 2.
- Image tester, as assertions that the padding is zero: after `fft` then `ifft`; for `self` after
  `ifft_mask_pad_fft`; after read, `fft`, read; after `ran` and `gauimg` on an image that was in
  Fourier space; after a scalar assignment, adding a constant, `norm` and an image division; after
  `bin_inv`, `binarize` with a threshold of zero and `ring` with an inner radius of zero; after `pad`
  with a background and `collage`. Even and odd boxes, 2D and 3D, where the operation takes them:
  `fft` and `ifft` stop on an odd box (`shift_phorig` assumes even dimensions), so the round trip is
  tested on even boxes, and `ifft_mask_pad_fft`, which is 2D and shifts the origin itself, on an
  even and an odd box. Expected to fail now; step 1 makes them pass and they stay as its pin. This
  list is the pin; it is not extended to every operation of the class.
- The padding values after `ifft`, recorded here for both platforms: FFTW does not specify them, so
  they can differ with the platform, the plan and the box parity.
- Flex: for the representative volumes of one flex case, the `em_basis` Gram over the padded arrays
  and over the box. Record the relative difference here. Zero means the flex numbers of that case on
  that platform do not move in step 1, and no more than that, for the reason above; anything else is
  a defect that step 1 fixes, and the before/after flex numbers are recorded with it. (No test suite
  runs the EM engine: `unit_heterogeneity` and `lib_heterogeneity` test the latent-model utilities and
  the PCG operator. The case is a real-data `flex_pca` run; section 8.)

## 4. Step 1: the image class keeps the padding zero (about a day)

- **Back to real space.** A private `zero_padding` helper (2 x n2 x n3 stores at most), called after
  the complex-to-real transform wherever the image stays in real space: `bwd_ft`, and
  `ifft_mask_pad_fft` for `self`. Its `self_out` ends in Fourier space (the forward transform at the
  end of the routine) and needs nothing, as does `ifft_mask_fft`, since the forward transform ignores
  the padding.
- **Into real space without a transform.** `zero_padding` after the box is filled in `read`,
  `read_raw_mrc` and `read_single_mrc_image`, and in `ran`, `gauran` and the `gauimg` family. The
  three readers clear the padding only when the image is flagged real: `read` takes the flag from
  the file mode, and the two raw readers keep the caller's flag and also read Fourier data (the
  gridding accumulators of the reconstructor, the Fourier stacks of `simple_discrete_stack_io`),
  where the padding rows are the last Fourier column. Each of the 33 `ft = .false.` statements of
  the image submodules is checked against this list; the check added the real-space outputs of
  `clip` and `rtsq`, which fill the box of an output that may have held Fourier data. `new`,
  `set_rmat`, `zero_and_unflag_ft`, `unserialize`, `corners`, `square`, `window`, `window_slim`,
  `pad_mirr` and `ft2img` already zero the whole buffer.
- **`set_ft` stays a flag setter.** It does not clear the padding: a caller that flags a buffer of
  Fourier data as real would lose the last Fourier column. Its three `.false.` callers are covered
  by what follows them: `window` zeroes the whole output (`simple_particle_extractor.f90`, twice), and
  the micrograph of `simple_reextract_strategy.f90` is read again, which clears the padding after
  this step.
- **Real-space writes.** The 63 whole-array `self%rmat` statements and the 26 through other names:
  the ones that can put a nonzero value into the padding in real space act on the box slice. By the
  test of section 2 these are the scalar assignments, adding or subtracting a constant, the mean and
  standard-deviation normalisations in `simple_image_norm.f90` and `simple_image_calc.f90`, the image
  division, `bin_inv`, the masked assignments of `binarize_1` and `ring`, the background fill of
  `pad` and the background fills of `collage`, and the first-dimension sections of `cendist`.
  Operations that keep zero padding zero (sums,
  differences and products of images, copies, multiplication or division by a constant) stay as they
  are. The Fourier-space branches of the same routines and the scalar fills of the `*_pad_fft` family
  are untouched.
- **Reductions.** The foreground count of `binarize_1`, `npix` of `ring`, the range checks of
  `remove_edge` and `one_at_edge`, `contains_nans` and `checkimg4nans` act on the box.
- **Tests.** The step-0 image-tester assertions pass. The fast gate on Linux and macOS Debug with bounds
  checking, and the image tester once in a Debug build with `--traps` (the image division no longer
  divides zero by zero in the padding). `lib_heterogeneity` nightly and one flex run, compared with
  step 0 and recorded here.

## 5. Step 2: the pointer is the box (about a day)

- `get_rmat_ptr` returns the box view, and the `! VIOLATES ENCAPSULATION` flag goes. A new, explicitly
  named accessor for the whole buffer (for example `get_rmat_ptr_padded`) serves the callers that need
  the FFT layout. From the survey that is one caller: `simple_denoise_project_strategy.f90`, which hands
  the buffer and the Fourier-space flag to `imgfile%wmrcSlices`. `simple_stack_io.f90` and
  `simple_flex_gpu.f90` slice to `ldim` and work with the box view. The image tester reads the
  padding through the query of step 0 and is not a caller of the padded accessor.
- **The callers, file by file.** Explicit indexing within `ldim` needs no change; whole-array
  reductions and updates become box reductions by construction.

  | Callers | Calls | What they do with the pointer |
  |---|---:|---|
  | `simple_nu_filter_apply`, `_state`, `_envmask`, `_stats`; `simple_pcg_solvent_sidecar`; `simple_flex_pca_em_compose`; `simple_gridding`; `simple_motion_align_nano`, `_hybrid`; `single_tseries_tracker`; `simple_ctf_estimate_cost` | 43 | index within `ldim`: no change |
  | `simple_flex_pca_em_basis`, `_em_mstep`, `_em_pairmerge`, `_em_solve`, `_em_crossfsc`, `simple_flex_pca_util` | 43 | whole-array Gram, norm, projection, deflation and zeroing: become box operations; results as after step 1. Two of the 43 only cleared the padding after an `ifft` (`rm_dfl(params%box_crop+1:,:,:) = 0.` in `_em_mstep` and `_em_pairmerge`): empty sections on a box pointer, removed with their calls |
  | `simple_flex_pca_pcg` | 21 | mostly indexed; passes the pointer to `occupancy` and `window_product` |
  | `simple_calpha_finder` | 12 | indexed; passes it to `suppress_neighborhood`; `zero_score_border` writes -1 through sections that are whole in the first dimension (`values(:,1:border,:)`), into the padding before this step, box by construction after |
  | `simple_nanoparticle` | 2 | passes it to `calc_isotropic_disp` and `calc_anisotropic_disp` |
  | `simple_segmentation`, `simple_ctf_estimate_fit` | 8 | indexed; whole-array `where` clipping (box by construction) |
  | `simple_stack_io`, `simple_flex_gpu` | 4 | slice to `ldim`: no change |
  | `simple_denoise_project_strategy` | 1 | the padded accessor (Fourier data to `wmrcSlices`, which itself writes the first `ldim(1)` rows only) |
  | testers and `simple_commanders_test_highlevel` | 24 | whole-array fills, `background_mean`, `ls_scale_profile`: box by construction |

  Not in the 158: `simple_image_bin.f90`, a separate module in the image directory, calls
  `get_rmat_ptr` five times and indexes within its own `bldim`; no change.

- **Contiguity.** The box view is not contiguous (the second-dimension stride is `2*fdim(n1)`). An
  assumed-shape dummy takes it as it is; an explicit-shape or assumed-size dummy gets a copy in and out,
  which is correct but costs time in a hot loop. Checked 2026-09-30: the seven procedures of the table
  that receive the pointer (`occupancy`, `window_product`, `suppress_neighborhood`,
  `calc_isotropic_disp`, `calc_anisotropic_disp`, `background_mean`, `ls_scale_profile`) have
  assumed-shape or pointer dummies, no pointer that takes it is declared `contiguous`, and none is
  passed to `c_loc`, so no dummy needs to change; repeat the check when the step is done. OpenMP loops
  that index the pointer element by element are unaffected: the inner loop runs over the first,
  contiguous dimension.
- **Bounds checking.** An index of `ldim(1)+1` or `ldim(1)+2` in the first dimension is inside the
  padded pointer today and outside the box pointer. A Debug build then stops on it; a Release build
  reads the same memory as before. A new bounds failure after this step is an existing read of the
  padding (a neighbour read at the box edge, for example) and is fixed in the caller, not by going
  back to the padded accessor.
- **Tests.** The fast gate on both platforms; `lib_heterogeneity`, `lib_single` and
  `lib_reconstruction` nightly; the times of `unit_heterogeneity` and `lib_heterogeneity` before and
  after, to catch a copy-in regression. The policy's section 4.6 trap entry is then history (the
  pointer is the box) and is rewritten to say so.

## 6. Not in scope

- Removing the pointer (step 3): weeks of work, moving about forty operations, flex-specific linear
  algebra among them, into the image class, with copies in the performance-critical loops of
  `nu_filter`, `calpha_finder` and flex. Steps 1 and 2 give the safety at a fraction of the cost.
- Padded FFTW buffers owned by other classes (`simple_polarft_corr`, `simple_classaverager_core`):
  they are not exposed through `image`; review them separately if needed.
- `get_rmat` and `get_rmat_sub` already return copies of the box.

## 7. Done when

- Outside the `simple_image` module and its submodules, no code holds a pointer to the padded buffer
  except through the named accessor (a grep for `get_rmat_ptr_padded` lists the one caller of
  section 5).
- The image tester asserts zero padding, through the padding query, after the operations listed in
  section 3.
- The fast gate passes on Linux and on macOS Debug with bounds checking; the library suites pass
  nightly; the flex results before and after, and the step-0 measurement, are recorded below.

## 8. Record

Steps 0, 1 and 2 were done on 2026-09-30 on master `f9bbf72ce`, on Linux (Oracle Linux 8.10, Intel
Arrow Lake, GCC 15.2, FFTW 3.3.5), in Release (`-O3`) and Debug (`-fcheck=all`) builds. The macOS
Debug gate of section 7 is open. Nothing is committed; the record below is what was measured.

### 8.1 Decisions

- **The padding query is `max_abs_padding`** (`simple_image_checks.f90`, a pure type-bound function):
  the name says what it returns, the largest magnitude in the padding rows, NaN when one of them is
  NaN. It is meaningful in real space; in Fourier space those rows are the last Fourier column.
- **The padded accessor is `get_rmat_ptr_padded`**: the name of section 7's done criterion, and it
  reads as `get_rmat_ptr` with the difference stated.
- **`zero_padding`** is a private, pure type-bound procedure in the zeroing section of
  `simple_image_ops.f90` (`self%rmat(self%ldim(1)+1:,:,:) = 0.`).
- **The step-0 flex measurement used throw-away copies of the tree**, never the work tree: the
  step-0 sources with an optional tag on `get_rmat_ptr` that counts, per call site of
  `src/main/flex`, the real-space images whose padding is not zero, plus the box Gram next to the
  padded Gram in `orthonormalize_representatives`. The same instrumentation on the step-1 sources,
  with a backtrace wherever a real-space image with nonzero padding is transformed, copied, written,
  pointed at or killed, checked step 1. None of it is in the diff.

### 8.2 Step 0: the padding before any change

Image tester, 220 assertions in the `image` sub-suite, 78 of them failing before step 1, all padding
assertions, the same 78 in Release and Debug. Largest magnitude in the padding, for boxes 64x64,
63x63, 24x24x24 and 23x23x23 (N(0,1) noise where a content is needed):

| operation | padding found (64 / 63 / 24 / 23) | assertions failing |
|---|---|---:|
| `fft`, `ifft` (even boxes) | 0.326 / - / 0.657 / - | 2 |
| `ifft_mask_pad_fft`, `self` (2D) | 0.326 / 0.187 | 2 |
| `read`, `fft_noshift`, `read`, MRC | 0.037 / 0.022 / 0.018 / 0.020, the Fourier data of the image before | 4 |
| the same through SPIDER | 0: that read path zeroes the whole array first | 0 |
| `ran`, `gauran`, `gauimg` on an image in Fourier space | 0.006 to 0.037 | 12 |
| scalar assignment, `add`, `subtr`, `+` with a constant | the constant: 1, 3, 2, 5 | 16 |
| image division | 2.5, the quotient of the two paddings (NaN when both are zero) | 4 |
| `norm`, `norm` to a mean and deviation | 1.5, 3.5 | 8 |
| `bin_inv`, `binarize` at zero (in place, into an output), `ring` | 1 | 16 |
| `npix` of `ring`, radius 10 | 359 for 317 / 336 for 316 / 4803 for 4169 / 4540 for 4224 | 4 |
| `pad` with a background, `collage` | the background: 5, 128 | 6 |
| `contains_nans` with a NaN planted in the padding | true | 4 |

The values after `ifft` are FFTW's and belong to this platform, plan and box; the Release and Debug
builds agree to six digits. `fft` and `ifft` stop on an odd box (`shift_phorig`), so an odd box has
padding after a transform only through `ifft_mask_pad_fft`.

Flex, a real-data `flex_pca` run (bgal, 5513 particles, box 256 cropped to 64, 10 components,
`nthr=24`, 57 s; no suite runs the EM engine):

- The `em_basis` Gram over the padded arrays equals the Gram over the box exactly: relative
  Frobenius difference 0 for all 12 Grams of the run (14 and 10 representatives). The
  representatives have zero padding: the deapodisation multiplies it by the zero padding of the
  correction image.
- Of the 47 `src/main/flex` call sites the run reaches, 11 see an image with nonzero padding. Nine
  are harmless: the PCG scratch images (`simple_flex_pca_pcg`, up to 9e7 in the padding) are sliced
  to the box, `set_rmat` takes the box, and the templates built on a copy of the consensus are
  zeroed before they are filled. Two are not: the consensus norm of the mean-shaped deflation in `em_mstep` and
  `em_pairmerge`, `mnorm_dfl = sqrt(sum(real(rv_dfl,dp)**2))`, runs over the padding that
  `read_and_crop` left (0.154 at most): 265.19509 for 265.19184 over the box, 1.2e-5 too large. It
  scales the threshold below which a deflation template is dropped (`1e-3 * mnorm_dfl`) and a log
  line. The outputs of the run do not change when the norm takes the box value (8.5).

### 8.3 Step 1: what changed in the image class

- `zero_padding` after the complex-to-real transform in `bwd_ft` and, for `self`, in
  `ifft_mask_pad_fft`; after the box is filled in `read`, `read_raw_mrc`, `read_single_mrc_image`
  (real images only, section 4), `ran`, `gauran`, `gauimg_1`, `gauimg_2`, `gauimg2D`, `gauimg3D`,
  the scalar assignment, and the real-space outputs of `clip` and `rtsq`.
- On the box instead of the whole buffer: the scalar assignment; `+` with a constant, `add` and
  `subtr` of a constant (their Fourier branches unchanged); the image division; the mean
  subtractions of `norm`, `norm_ext`, `norm_noise`, `norm_within`, `cure` and `prenorm4real_corr_3`
  and the offset of `norm`; the constants of `remove_neg`, `zero_background`, `zero_edgeavg` and
  `zero_env_background`; `bin_inv`; the masked assignments of `binarize_1` and `ring`; the
  background fills of `pad` and `collage`; the first-dimension sections of `cendist`.
- Reductions on the box: the foreground count of `binarize_1`, `npix` of `ring`, the range checks of
  `remove_edge` and `one_at_edge`, `contains_nans`, `checkimg4nans`.
- `reshape2cube` is deleted with its binding and interface: it has no caller in `src` or
  `production`, and its whole-array assignment between buffers of different shape does not conform.
- `set_ft` is unchanged.

All 220 padding and image assertions pass after step 1, in Release, Debug and a Debug build with
`SIMPLE_DEBUG_TRAPS` (the image division no longer divides in the padding). The step-1 sources with
the backtrace instrumentation of 8.1 went through the fast gate, the library tier, `pcg_recon`,
`single_atoms_stats`, `single_workflow`, `simulated_workflow_1jxy` and the flex run: no real-space
image with nonzero padding anywhere (the same instrumentation on the step-0 sources reports them
from the first `ifft` of the image tester on). In the flex run all 47 call sites see zero padding
and `mnorm_dfl` is the box value.

### 8.4 Step 2: what changed at the callers

- `get_rmat_ptr` points at `self%rmat(:ldim(1),:ldim(2),:ldim(3))`; the `! VIOLATES ENCAPSULATION`
  flag is gone. `get_rmat_ptr_padded` has one caller, `write_diffmap_stack_image` of
  `simple_denoise_project_strategy.f90`.
- Two caller edits besides that one: the explicit padding clears of `em_mstep` and `em_pairmerge`
  (section 5) are removed with their `get_rmat_ptr` calls, so the 158 calls are now 156 of
  `get_rmat_ptr` and one of `get_rmat_ptr_padded`; and the comment of the `cif2mrc` tester that
  called the pointer padded is removed. No other caller needed a change.
- **No bounds failure was exposed.** Debug with bounds checking: the fast gate, the library tier,
  `pcg_recon`, `single_atoms_stats`, `single_workflow` (35 min) and a real-data `flex_pca` run
  (208 s) pass. `simulated_workflow_1jxy` stops in Debug at `simple_abinitio_utils.f90:509`
  (`res(:nbest)` with `nbest = 3` and one 2D class), in a file this work does not touch and
  identically on a Debug build of the base tree; with that one section bounded in a scratch copy,
  the gate passes in Debug on the step-2 sources. It is reported, not fixed here.
- **Contiguity**, checked again: every procedure that receives the pointer has an assumed-shape or
  pointer dummy (`occupancy`, `window_product`, `suppress_neighborhood`, `calc_isotropic_disp`,
  `calc_anisotropic_disp`, `background_mean`, `ls_scale_profile`, `pick`, `set_rmat`,
  `wmrcSlices`), no `contiguous` pointer takes it and none goes to `c_loc`. The Debug runs print
  gfortran's "array temporary was created" warning at the same four places of the fast gate before
  and after the step. One of them involves the pointer (`simple_segmentation.f90:671`, a section
  of it passed to the explicit-shape dummy of `neigh_8`); it was a section of the padded array,
  and a copy, before the step as well.
- `doc/policies/test_environment_policy.md`, section 4.6: the trap entry now says the pointer is the
  box.

### 8.5 Flex and the other numbers, before and after

The real-data `flex_pca` run (91 output files: 66 maps, 20 text tables, 5 binary files):

| comparison | identical files | largest relative difference, text tables | maps |
|---|---:|---:|---|
| base, run 1 against run 2 | 85 | 2.6e-15 | identical |
| base against step 1 | 85 | 2.3e-15 | identical |
| base against step 2 | 85 | 1.7e-15 | identical |
| step 1 against step 2 | 85 | 2.7e-15 | identical |

The six files that differ are the four probe tables, at double-precision round-off, and two binary
files (the cross-FSC series and the embedding cache); the eigenvalue tables, the coordinates of
all 5513 particles and every map are byte-identical. The 72 metric lines of `unit_heterogeneity` and the 55 of `lib_heterogeneity`
(noise calibration, deconvolution likelihoods and errors, the operator, adjoint and solve-sweep
figures of the PCG test) are identical in the base, step-1 and step-2 builds.

The three highlevel gates that run threaded pipelines are not reproducible run to run on either
build (their first threaded stage already differs between two runs of the base build), so a number
of the work build is compared with the spread of the base build:

| gate (Release) | base, two runs | step 2 |
|---|---|---|
| `single_atoms_stats`: recall, precision, RMS position error, density correlation | 1.000, 1.000, 0.015 A, 1.000 (both) | the same (two runs; step 1: the same) |
| `single_atoms_stats`: relative RMS difference of the simulated-density map between two runs | 3.8e-4 | 6.6e-4 against base, 7.6e-4 between two work runs |
| `simulated_workflow_1jxy`: whole-volume correlation, FSC = 0.143 | 0.850 at 5.9 A; 0.801 at 7.6 A | 0.847 at 6.4 A |
| `single_workflow`: whole-volume correlation, FSC = 0.143 | 0.415 at 4.8 A; 0.460 at 3.6 A | 0.463 at 3.6 A |
| `single_workflow`: correlation of the final maps of two runs | 0.775 | 0.767 and 0.869 against the two base runs |

Every gate passes its floors in every run; the work build is inside the spread of the base build.

### 8.6 Tests

| tier | base | step 0 | step 1 | step 2 |
|---|---|---|---|---|
| fast gate, Release | - | 12 of 13 (`unit_image`: the 78) | 13 of 13 | 13 of 13 |
| fast gate, Debug | - | 12 of 13 (the same 78) | 13 of 13 | 13 of 13 |
| `unit_image`, Debug with traps | - | - | pass | pass |
| library, Release | 5 of 5 | - | 5 of 5 | 5 of 5 |
| library, Debug | - | - | - | 5 of 5 |
| `pcg_recon`, `single_atoms_stats`, `single_workflow`, Release | - | - | - | pass |
| the same three, Debug | - | - | - | pass |
| `simulated_workflow_1jxy`, Release | - | - | - | pass |
| `simulated_workflow_1jxy`, Debug | fails (8.4) | - | - | fails as base; passes with the guard |

`unit_image` has 783 assertions after step 2 (755 after step 1): the `image` sub-suite went from 76
to 248. `lib_reconstruction` passes on the base and on the work build (43 assertions).

Caller coverage under bounds checking, from a Debug build of the step-2 sources with a tag on each of
the 161 call sites (156, and 5 in `simple_image_bin.f90`); the number of call lines each run executes:

| caller file | call lines | fast gate | library | pcg_recon | single_atoms_stats | single_workflow | simulated_workflow_1jxy | flex_pca run | executed by any |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `simple_stack_io.f90` | 3 | 3 | 2 | - | - | 2 | 2 | - | 3 |
| `simple_stack_io_tester.f90` | 13 | 13 | - | - | - | - | - | - | 13 |
| `simple_commanders_test_highlevel.f90` | 9 | - | - | - | - | - | - | - | 0 |
| `simple_ctf_estimate_cost.f90` | 2 | - | - | - | - | - | 2 | - | 2 |
| `simple_ctf_estimate_fit.f90` | 3 | - | - | - | - | - | 3 | - | 3 |
| `simple_flex_gpu.f90` | 1 | - | - | - | - | - | - | - | 0 |
| `simple_flex_pca_em_basis.f90` | 10 | - | - | - | - | - | - | 7 | 7 |
| `simple_flex_pca_em_compose.f90` | 4 | - | - | - | - | - | - | - | 0 |
| `simple_flex_pca_em_crossfsc.f90` | 2 | - | - | - | - | - | - | 2 | 2 |
| `simple_flex_pca_em_mstep.f90` | 9 | - | - | - | - | - | - | 9 | 9 |
| `simple_flex_pca_em_pairmerge.f90` | 12 | - | - | - | - | - | - | 12 | 12 |
| `simple_flex_pca_em_solve.f90` | 7 | - | - | - | - | - | - | 7 | 7 |
| `simple_flex_pca_pcg.f90` | 21 | 7 | 7 | - | - | - | - | 7 | 9 |
| `simple_flex_pca_util.f90` | 1 | - | - | - | - | - | - | 1 | 1 |
| `simple_image_bin.f90` | 5 | 4 | 3 | - | 3 | 3 | 2 | - | 4 |
| `simple_segmentation.f90` | 5 | 5 | 1 | - | 1 | 1 | 1 | - | 5 |
| `simple_gridding.f90` | 1 | 1 | 1 | - | - | 1 | 1 | - | 1 |
| `simple_motion_align_hybrid.f90` | 1 | - | - | - | - | - | 1 | - | 1 |
| `simple_motion_align_nano.f90` | 1 | - | - | - | - | - | - | - | 0 |
| `simple_calpha_finder.f90` | 12 | 11 | - | - | - | - | - | - | 11 |
| `simple_calpha_finder_tester.f90` | 1 | 1 | - | - | - | - | - | - | 1 |
| `simple_cif2mrc_tester.f90` | 1 | 1 | - | - | - | - | - | - | 1 |
| `simple_nanoparticle.f90` | 2 | - | - | - | 2 | 1 | - | - | 2 |
| `single_tseries_tracker.f90` | 1 | - | - | - | - | - | - | - | 0 |
| `simple_nu_filter_apply.f90` | 16 | - | - | - | - | - | 6 | - | 6 |
| `simple_nu_filter_envmask.f90` | 2 | - | - | - | - | - | - | - | 0 |
| `simple_nu_filter_state.f90` | 6 | - | - | - | - | - | 2 | - | 2 |
| `simple_nu_filter_stats.f90` | 1 | - | - | - | - | - | 1 | - | 1 |
| `simple_denoise_project_strategy.f90` | 1 | - | - | - | - | - | - | - | 0 |
| `simple_pcg_solvent_sidecar.f90` | 8 | - | - | - | - | - | - | - | 0 |

The `flex_pca run` column is the real-data run, not a suite: no suite reaches the EM callers. The
two long gates were traced on the scratch copy that bounds the section of 8.4.

102 of the 161 call lines are executed. The other 59 were read line by line: each indexes or
slices within `ldim`, or is a whole-array expression that is a box expression now. Eight files have
no executed call: `simple_commanders_test_highlevel.f90` (all nine calls are in `rec3D_backends`,
which needs a user's project and is not a CTest entry), `simple_flex_gpu.f90` (GPU only; slices to
the box), `simple_flex_pca_em_compose.f90`, `simple_motion_align_nano.f90`,
`single_tseries_tracker.f90`, `simple_nu_filter_envmask.f90`, `simple_pcg_solvent_sidecar.f90` and
`simple_denoise_project_strategy.f90` (the padded accessor).

### 8.7 Timings and corrections

Wall time of the CTest entry, Release, median of three runs with nothing else running:

| entry | base | step 2 | ratio |
|---|---:|---:|---:|
| `unit_image` | 0.316 s | 0.329 s | 1.04 |
| `unit_heterogeneity` | 1.293 s | 1.295 s | 1.00 |
| `lib_heterogeneity` | 2.948 s | 2.996 s | 1.02 |

`unit_image` runs 172 more assertions than on the base build (the `image` sub-suite takes 0.06 s for
0.04 s), which is the difference: nine more runs give 0.309 s and 0.328 s, and the five heavier
sub-suites whose tests did not change (`masks`, `nano_mask`, `binary_image`, `segmentation`,
`volume_shape`), timed alone over nine runs, sum to 0.309 s on base and 0.307 s on the work build.
No copy-in regression.

Corrections made to sections 1 to 5 during the work:

- Section 1: 28 whole-array writes through other names, not 26; five reductions over the padded
  array, not four (`checkimg4nans`); `cendist` and its effect on `ring`.
- Section 3: `fft` and `ifft` take even boxes only; no suite runs the flex EM engine, so the flex
  case is a real-data run.
- Section 4: the readers clear the padding of real images only; `clip` and `rtsq` added to the
  routines that enter real space; `cendist` added to the real-space writes; `checkimg4nans` added to
  the reductions.
- Section 5: `zero_score_border` of `simple_calpha_finder` wrote into the padding; two of the flex
  calls only cleared the padding and are gone; `wmrcSlices` writes the first `ldim(1)` rows also of
  Fourier data, so the padded accessor of the denoise caller is the layout it was given, not a
  need of the file format.
