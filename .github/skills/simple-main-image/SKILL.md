---
name: simple-main-image
description: Use when working in SIMPLE's src/main/image or src/main/image_processing subsystems, including the image type (one FFTW buffer seen as rmat/cmat, the Fourier-state flag, FFT normalisation and phase origin, Fourier addressing, pad/clip, copies, FFTW plans and threads, process-wide memos), its submodules (FFT, filters, masks, CTF application, geometry, I/O), image_bin and image_msk, the polarizer and projector, ft_expanded for motion correction, the reconstruction plane grid, the wave field and imaging, and the image-array, stack-file, micrograph and segmentation helpers.
---

# `src/main/image`: the image type and its Fourier conventions

Owns `image`: one FFTW buffer seen as `rmat` (real, with padding rows) or `cmat` (the stored Hermitian
half), the flag `ft`, the transforms and the operations on one image (filters, masks, CTF application,
geometry, I/O), plus the polarizer (`simple_image_polar`) and projector (`simple_projector*`) that cut
sections for `pftc`, `ft_expanded` (motion), `fplane_grid` (planes of `gen_fplane4rec`) and the wave
field and imaging (`doc/policies/3D/wave_model_policy.md`). It does not prepare particles
(`strategies/search`, `prepimg4align*`), reconstruct (`volume/`), model the CTF (`ctf/`) or score
(`pftc/`, `cftc/`). `src/main/image_processing` holds procedures on images (arrays, stack-file loops,
micrograph and segmentation helpers): no type, reads a project, never writes one; these rules apply.
Details and the memo table: `references/fourier-and-memos.md`.

## Rules

Silent ones first: nothing reports these.

- Know the state; `set_ft` only flips the flag. Methods transform on demand and restore (`shift`,
  `corr`, `bp`), stop on the wrong state (`window`, `pad_mirr`, `expand_cmat`), or trust it unchecked
  (`shift2Dserial`, `apply_ctf`, `polarize*`, `mul_cmat`; `add` reads the operand in self's state).
- `fft` divides by N = product(ldim), `ifft` does not; `fft` swaps quadrants first so the real-space
  origin is the centre pixel ldim/2+1, even boxes only (THROW_HARD). `fft_noshift` takes the origin at
  (1,1,1), any box; `ifft` always swaps back. `comp_addr_phys` maps h < 0 to the Friedel mate without
  conjugating: use `get_fcomp`/`set_fcomp`, which conjugate, or conjugate yourself.
- `shift(s)`, `shift2Dserial(s)`: out(x) = in(x + s), the content moves by -s (phase exp(+2 pi i h.s/N)).
  A particle is centred on its stored shift with `-shift`; cftc uses the same sign.
- In real space the padding rows hold zeros (the class's reductions rely on it): outside the class use
  `get_rmat_ptr` (the box, not contiguous); `get_cmat_ptr` is in place only, then `set_ft` by hand.
- `pad`/`clip` act in the current state: in Fourier space they resample (the output's smpd changes;
  `pad` applies a Hann anti-aliasing weight by default), in real space they embed or cut the box. The
  output must exist with the target ldim; the wrong direction does nothing.
- Copies: scalar `=` and `copy` reallocate via `new` and take the source's `wthreads`; `copy_fast` needs
  equal ldim (unchecked) and keeps the receiver's. Array assignment and `allocate(source=)` copy the
  pointers: two owners of one buffer. No finaliser or elemental assignment, and none is to be added:
  `kill` every image (`dealloc_imgarr` for arrays); operator results and the `image()` constructor
  leak a buffer and two plans each.
- Process-wide memos (mask coordinates, `memoize_ft_maps` for every CTF application, the polarizer's
  geometry, `construct_thread_safe_tmp_imgs`, `ftexp_transfmat_init`): make them serially, before the
  parallel loop, for the box in use (`apply_ctf_wpad` on a real-space image works at twice the box).
  `apply_ctf` and `ctf2img` stop on maps of another box or sampling (`ft_maps_match`); some readers still do not
  check (the `*_ctf_flip` kernels).
- Polarizer: one geometry per process, the last `memoize4polarize`/`_oversamp` wins (both forms stay);
  `polarize*` take the caller's `pdim` and stop unless box, `pdim` and mode match the memo, and the
  memo stops inside an OpenMP region: memoize with the pdim given to `polarize_*_pft` and pass the
  matching `oversamp`. Its angles duplicate `polarft_calc%polar` (`pftc/simple_polarft_core.f90`): change both.
- Projector order and `cmat_exp` (a copy: re-expand after a change): the reference.
- `image_bin` syncs `bimat` from rmat only once: after changing rmat call `set_imat`.
- Threads: `new`, `copy`, `=`, `set_wthreads`, `kill` plan or destroy inside an unnamed `omp critical`.
  In a parallel loop use per-thread images made before it (`wthreads=.false.`, sized by the team),
  `copy_fast` and the `_serial` variants. Distinct images transform concurrently; share one image only
  if every thread just reads it in its current state (on-demand transforms write). Reproducibility
  limits (Fourier crop, `fft_measure`): the reference.

## Where to start

- `simple_image.f90`: the type, every binding (grouped under the name of its submodule) and interface.
- `simple_ftiter.f90` (addressing, loop limits); `simple_image_tester.f90` (the conventions as assertions).
- Callers that show the order: `strategies/search/simple_matcher_refvol_utils.f90` (projector),
  `simple_matcher_ptcl_batch.f90` (polarizer and CTF maps around the particle loop).

## Tests

- `build/bin/simple_test_exec test=unit_image` (about 1 s; `suite=` image, masks, binary_image,
  fourier_iterator, ctf, segmentation, ...). Elsewhere: `unit_numerics suite=shift_search` (ft_expanded),
  `unit_reconstruction suite=wave_field`, `wave_imaging`, `observation_noise`.
- Polarizer or projector: `unit_pftc_align2D3D suite=polar_correlation`, `unit_cart_align3D
  suite=cart_calculator`, `test=simulate_particles` (11 s; images against their reprojections).
- Alignment or maps move: `test=simulated_workflow`, a robustness profile (one run shows breakage, not
  improvement). New OpenMP code: a forced team against the serial result.

## See also

`doc/policies/KB_Interpolation_Policy.md`, `doc/refactoring_notes/completed/image_rmat_padding_encapsulation.md`;
skills `simple-main-interp`, `simple-main-pftc`, `simple-main-cftc`, `simple-main-sigma2`, `simple-modern-fortran`.
