# image: Fourier layout, addressing, memos and threading in detail

Reference for `simple-main-image`; read when writing a loop over Fourier components, a new memo, a
parallel caller or a new section cutter.

## Buffer and layout (`simple_image_core.f90`, `new`)

- One FFTW allocation: `cmat(fdim(n1), n2, n3)` and `rmat(2*fdim(n1), n2, n3)`, `fdim(n) = n/2 + 1`.
  The real view has 2 padding rows (even n1) or 1 (odd n1) in its first dimension.
- `new` on an existing image of the same ldim keeps the buffer, remakes the plans, zeroes it, sets
  `ft = .false.` and takes `wthreads` from the argument (default on). A different ldim kills first.
- `get_rmat` and `get_cmat` copy; `get_rmat_ptr` is the box (lower bounds 1, not contiguous: pass it to
  assumed-shape or pointer dummies); `get_rmat_ptr_padded` is the whole buffer for FFT-layout callers;
  `get_cmat_ptr` contract in `simple_image_access.f90` (in place only, `set_ft` yourself, drop before
  `kill`/`new`).
- `ifft` and `ifft_mask_pad_fft` clear the real padding; routines that fill only the box (`read`, `ran`,
  `gauran`, ...) were made to clear it too (`image_rmat_padding_encapsulation.md`). A new routine that
  leaves Fourier space without a transform must do the same.

## Logical and physical addresses (`simple_ftiter.f90`)

- Stored half: h = 0..n1/2 at physical h+1. k (and l in 3D) from -n/2 to n/2-1 (even n) or
  -(n-1)/2..(n-1)/2 (odd n); negative k at physical k+1+n. 2D images have l = 0 only.
- `loop_lims(2)`: the stored half. `loop_lims(3)`: every logical address, the redundant h < 0 half
  included (h from -n1/2). `loop_lims(1, lp)`: band-limited for correlations, h from the high-pass
  index, k and l symmetric, at least 3.
- `comp_addr_phys(h,k,l)` with h < 0 returns the address of (-h,-k,-l); only `get_fcomp`, `set_fcomp`,
  `get_fcomp2D` and the gathers conjugate. `comp_addr_logi` inverts it on the stored half.
- Resolution: `get_lp(1, ind)` = dsteps/ind with dsteps = n1*smpd (even n1); `get_find` inverts it;
  `get_nyq` = `get_filtsz` = n1/2.
- Shift constant: `shconst(i)` = 2 pi/n for every n (odd included; the same bits as before for even n);
  only `fft_noshift` paths reach an odd box (`fft` refuses odd boxes). The polar calculator's argument
  transfer (`polarft_core`) and the class-sum shift (`cavg_sums%shift_class`) use the same 2 pi/n.

## Phase origin and shifts

- `fft`: `shift_phorig` (quadrant swap; THROW_HARD for odd dims), r2c, divide by product(ldim).
  `ifft`: c2r, clear the padding, swap back. `fft_noshift`: no swap, any box, no inverse of its own.
- `shift`/`oshift`: multiply by exp(+i sum(h_i s_i shconst_i)), so out(x) = in(x + s)
  (`simple_image_tester.f90`, `test_shifts`). `shift` transforms on demand; `shift2Dserial` and
  `shift2D_demoivre` assume Fourier space.

## pad and clip (`simple_image_geom.f90`)

- Fourier state: `pad` copies the stored half into the larger grid, weighted by `hannw` unless
  `antialiasing=.false.`, and sets the output's smpd to smpd*n_in/n_out; `clip` crops it, smpd coarser.
- Real space: `pad` centres the box at (n_out-n_in)/2+1 (optional `backgr` fill); `clip` cuts the
  centre. Equal ldim: `pad` calls `copy` (reallocates, takes the source's `wthreads`).
- The output must already exist with the target ldim; the opposite size relation is a silent no-op.
  `pad_inplace`, `clip_inplace` and `pad_mirr(ldim)` reallocate self (keeping its `wthreads`).
- The Fourier crop round trip is exact in mathematics, not in single precision.
- `pad_fft` pads in real space and transforms (the projector's input); the `*_clip_shift*` and
  `*_mask_*pad_fft` routines are fused particle-preparation kernels used by `strategies/search`.

## Process-wide memos

| memo | made by | read by | guard |
| --- | --- | --- | --- |
| mask coordinates | `memoize_mask_coords` | `mask2D_*`, `mask3D_*`; the fused `mask_fft`, `mask_pad_fft`, `ifft_mask_fft`, `ifft_mask_pad_fft` | THROW_HARD inside a region; `mask2D_*`/`mask3D_*` re-memoize on a box change outside a region and stop inside; the fused kernels do not check |
| power-spectrum geometry | `memoize_powspec_coords(mskrad)` | `norm_noise_mask_fft_powspec` | the same, keyed on box and mskrad |
| CTF Fourier maps | `memoize_ft_maps(ldim, smpd)` (`simple_memoize_ft_maps.f90`), `forget_ft_maps` | `ctf2img`, `apply_ctf`, `apply_ctf_wpad` (in real space it pads to 2*ldim and needs the maps of that box), the `*_ctf_flip` kernels, `ft_map_ctf_kernel` | THROW_HARD inside a region; `ctf2img` and `apply_ctf` stop on another box (`ft_maps_match`; `prepimg4align_cart` checks by `ft_map_get_farray_shape`), the `*_ctf_flip` kernels do not check |
| polarizer geometry | `memoize4polarize(pdim)` (KB alpha 1, unpadded image; kept for corrmat, `strategy2D_utils`, graphene subtraction), `memoize4polarize_oversamp(pdim)` (KBALPHA, image padded by OSMPL_PAD_FAC, result times pf^2) | `polarize`, `polarize_oversamp` (given the caller's `pdim`) | THROW_HARD inside a region; the readers stop unless box, `pdim` and mode match (`polar_memo_valid`); one geometry per process |
| per-thread work images | `construct_thread_safe_tmp_imgs(nthr)` | `calc_shiftcen_serial` (index `omp_get_thread_num()+1`) | THROW_HARD inside a region; a box mismatch stops, an index beyond nthr does not |
| ft_expanded transfer matrix | `ftexp_transfmat_init(img, lp)` | `ft_expanded%shift`, the shift-search costs | none |

`unmemoize_mask_coords` and `unmemoize_powspec_coords` free the first two (testers leave no state).

## Polar grid

pdim = [pftsz, kfrom, kto]. Rotation i = 1..pftsz at angle (i-1)*pi/pftsz; sample (h,k) =
(k*sin, -k*cos). `polarft_calc%new` builds the same grid over 2*pftsz rotations
(`pftc/simple_polarft_core.f90`); the projector's polar batch takes its coordinates from
`pftc%get_polar_coords`. `pftc%get_pdim_srch` ends at kfromto(2), `get_pdim_interp` at interpklim.

## Projector

- Order: `pad_fft` to box_croppd, `expand_cmat`, gather (read-only, shared by every thread), then `kill`
  and `kill_expanded` (`kill` keeps `cmat_exp`). `cmat_exp` is a copy: re-expand after a change.
- `expand_cmat`: FT 3D volume only (THROW_HARD); copies the stored half into a cyclic, Friedel-complete
  halo array `cmat_exp` with KB(KBWINSZ, KBALPHA).
- `interp_fcomp(loc)`: KB gather at loc in the volume's own units (`fproject`, `fproject_serial`:
  Cartesian planes). `interp_fcomp_oversamp(loc)`: loc in native units, multiplied by OSMPL_PAD_FAC,
  result times pf^3 (`fproject_polar_batch*`: the polar reference bank). All gathers are pure.
- `fproject_polar_batch_opt` evaluates weights inline and holds at most 3 taps (THROW_HARD);
  it compiles for GPU offload only with `USE_OPENMP_OFFLOAD`.

## Threads

- Plans: `FFTW_ESTIMATE` by default, inside an unnamed `omp critical` in `new`, `set_wthreads` and
  `kill`. Threaded plans only with `wthreads`, more than one plan thread (`nthr_glob` or
  `fft_nthreads`) and a dimension of at least 200 (or ldim(3) of at least 100).
- `fft_measure=.true.` plans by measurement within 10 s and keeps the wisdom for the process; the PCG
  operators and the wave projector use it from box 200 and never under a fixed seed (measured plans
  defeat exact replay).
- A `wthreads` image opens parallel regions in `fft`, `copy_fast`, the quadrant swap, `bp` and some
  filters; per-thread images are made with `wthreads=.false.` (`alloc_imgarr` defaults to it). Some
  methods open a region whatever the flag (`shift`, `fproject`, `add_workshare`).
- Variants for use inside a parallel loop: `shift2Dserial`, `rtsq_serial`, `fproject_serial`,
  `apply_filter_serial`, `insert_lowres_serial`, `calc_shiftcen_serial`, `corr_unnorm_serial`
  (ft_expanded). The other methods open their own parallel region.
- Copy policy (`test_copy_and_inplace_threading`): `copy` and `=` clone the source's `wthreads`;
  `copy_fast`, `pad_mirr`, `real_space_filter`, `upsample_square_background` keep the receiver's.

## Smaller contracts

- `image_bin`: `bimat` (integer) beside rmat; `set_imat` syncs bimat from rmat (or sets it and updates
  rmat), `update_img_rmat` the reverse; methods call `set_imat` only while `bimat_is_set` is false.
  `kill` is overridden to free bimat.
- `fplane_grid`: shell table of the padded plane box plus 1/sigma2 and 1/sqrt(sigma2) resampled to its
  shells; `euclid_sigma2%prep_fplane_grid` keeps it while `matches(box_croppd, kfrom)`, and
  `read_global` kills it when the spectra change.
- `wave_field`: its own buffer and plans, FFT order (no centre swap, no `ft` flag); one field per thread,
  never copied (intrinsic assignment aliases the buffer). Scales: `wave_model_policy.md` section 3.1.
- `ft_expanded%shift(shvec, self_out)` writes the shifted copy into `self_out`; there is no in-place form.
