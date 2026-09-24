# The padded real-space array of `image`: zero padding, and a box-shaped pointer

Date: 2026-09-25

Status: planned, not started (Hans, 2026-09-25: steps 1 and 2 below; step 3, removing the pointer,
is out of scope).

Validation level: static inspection of the image class and of every caller of `get_rmat_ptr`. Nothing
was compiled or run; the effect on flex PCA results is unknown until step 0 measures it.

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
  - The image class writes into the padding in real space. Of the 63 whole-array `self%rmat = ...`
    statements in the image submodules, the scalar broadcasts put values there: assigning a scalar
    (`img = 1.`, `rmat = realin`), adding or subtracting a constant, `norm` (subtract the mean, divide
    by the standard deviation), and dividing one image by another (`self1%rmat/self2%rmat`: 0/0, a
    NaN, in zero padding).
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
- Inside the image class, a real-space operation that writes a value into the padding acts on the box.

## 3. Step 0: measure first

Before any change, so that what steps 1 and 2 change is known:

- Image tester: the padding after `fft` then `ifft`, after `norm`, after adding a constant and after an
  image division, as assertions that it is zero. Expected to fail now; step 1 makes them pass and they
  stay as its pin.
- Flex: for the representative volumes of one `unit_heterogeneity` flex case, the `em_basis` Gram over
  the padded arrays and over the box. Record the relative difference here. Zero means the flex numbers
  do not move in step 1; anything else is a defect that step 1 fixes, and the before/after flex numbers
  are recorded with it.

## 4. Step 1: the image class keeps the padding zero (about half a day)

- **Back to real space.** A private `zero_padding` helper (2 x n2 x n3 stores at most), called after
  the complex-to-real transform wherever the image stays in real space: `bwd_ft` and
  `ifft_mask_pad_fft` (and its `self_out`, which it fills in real space). `ifft_mask_fft` returns to
  Fourier space and needs nothing, since the forward transform ignores the padding.
- **Real-space writes.** The 63 whole-array `self%rmat` statements: the ones that can put a value into
  the padding in real space (scalar assignment, adding or subtracting a constant, the mean and
  standard-deviation normalisations in `simple_image_norm.f90` and `simple_image_calc.f90`, the image
  division) act on the box slice. Operations that keep zero padding zero (sums, differences and
  products of images, copies) stay as they are. The Fourier-space branches of the same routines are
  untouched.
- **Other ways into real space.** Check `read` into an image that held Fourier data and every
  `set_ft(.false.)` after a write into `cmat`; `new`, `set_rmat` and `zero_and_unflag_ft` already
  zero the whole buffer.
- **Tests.** The step-0 image-tester assertions pass. The fast gate on Linux and macOS Debug with bounds
  checking. `lib_heterogeneity` nightly and one flex run, compared with step 0 and recorded here.

## 5. Step 2: the pointer is the box (about a day)

- `get_rmat_ptr` returns the box view, and the `! VIOLATES ENCAPSULATION` flag goes. A new, explicitly
  named accessor for the whole buffer (for example `get_rmat_ptr_padded`) serves the callers that need
  the FFT layout. From the survey that is one caller: `simple_denoise_project_strategy.f90`, which hands
  the buffer and the Fourier-space flag to `imgfile%wmrcSlices`. `simple_stack_io.f90` and
  `simple_flex_gpu.f90` slice to `ldim` and work with the box view.
- **The callers, file by file.** Explicit indexing within `ldim` needs no change; whole-array
  reductions and updates become box reductions by construction.

  | Callers | Calls | What they do with the pointer |
  |---|---:|---|
  | `simple_nu_filter_apply`, `_state`, `_envmask`, `_stats`; `simple_pcg_solvent_sidecar`; `simple_flex_pca_em_compose`; `simple_gridding`; `simple_motion_align_nano`, `_hybrid`; `single_tseries_tracker`; `simple_ctf_estimate_cost` | 43 | index within `ldim`: no change |
  | `simple_flex_pca_em_basis`, `_em_mstep`, `_em_pairmerge`, `_em_solve`, `_em_crossfsc`, `simple_flex_pca_util` | 43 | whole-array Gram, norm, projection, deflation and zeroing: become box operations; results as after step 1 |
  | `simple_flex_pca_pcg` | 21 | mostly indexed; passes the pointer to `occupancy` and `window_product` |
  | `simple_calpha_finder` | 12 | indexed; passes it to `suppress_neighborhood` |
  | `simple_nanoparticle` | 2 | passes it to `calc_isotropic_disp` and `calc_anisotropic_disp` |
  | `simple_segmentation`, `simple_ctf_estimate_fit` | 8 | indexed; whole-array `where` clipping (box by construction) |
  | `simple_stack_io`, `simple_flex_gpu` | 4 | slice to `ldim`: no change |
  | `simple_denoise_project_strategy` | 1 | the padded accessor (Fourier data to `wmrcSlices`) |
  | testers and `simple_commanders_test_highlevel` | 24 | whole-array fills, `background_mean`, `ls_scale_profile`: box by construction |

- **Contiguity.** The box view is not contiguous (the second-dimension stride is `2*fdim(n1)`). An
  assumed-shape dummy takes it as it is; an explicit-shape or assumed-size dummy gets a copy in and out,
  which is correct but costs time in a hot loop. Check the dummies of the procedures in the table that
  receive the pointer, and make the hot ones assumed-shape. OpenMP loops that index the pointer
  element by element are unaffected: the inner loop runs over the first, contiguous dimension.
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

- Outside the image submodules, no code holds a pointer to the padded buffer except through the named
  accessor (a grep for `get_rmat_ptr_padded` lists the callers of section 5).
- The image tester asserts zero padding after every operation that ends in real space.
- The fast gate passes on Linux and on macOS Debug with bounds checking; the library suites pass
  nightly; the flex results before and after, and the step-0 measurement, are recorded below.

## 8. Record

(Update as the steps land: the step-0 measurements, the flex numbers before and after, the dummies
changed to assumed-shape, and the timings.)
