---
name: simple-main-interp
description: Use when working in SIMPLE's src/main/interp subsystem, including the Kaiser-Bessel kernel (KBWINSZ, KBALPHA, OSMPL_PAD_FAC, beta, the fast polynomial and device forms, normalized 2D/3D stencils and their gradients, the instrument function), gridding deapodization (the stencil envelope at the insertion period), window functions and square windows, and the neutral Cartesian Fourier lattice layer (centred embed and crop, native plane extraction, packed KB and B-spline gathers) used by reconstruction, PCG, class averaging, polarization and Cartesian alignment.
---

# `src/main/interp`: interpolation kernels and the Cartesian Fourier lattice

Owns the Kaiser-Bessel kernel (`kbinterpol`), the gridding corrections that undo it
(`simple_gridding`), window functions (`winfuns`, `edges_sqwins`) and the neutral Cartesian lattice
layer (`simple_cartesian_fourier`) shared by the PCG code (`volume/`), the record buffer
(`strategies/search`) and the Cartesian calculator (`cftc/`). It evaluates kernels and single gathers;
loops over a plane or bank live with their owners (polarizer and projector in `image/`, reconstructors
in `volume/`, class sums in `class/`). Rule set: `doc/policies/KB_Interpolation_Policy.md` (padded
units, pf^2/pf^3 scales, gather/splat duality); read it before changing a kernel, padding or scale.

## Rules

Silent ones first: nothing reports these.

- One kernel: `kbinterpol(KBWINSZ, KBALPHA)`, KBALPHA = OSMPL_PAD_FAC (`src/defs/simple_defs.f90`), in the
  projector, the oversampled polarizer, both reconstructors, the class averager, cftc's taper and
  `ref_weight` and every deapodization envelope. A producer and its consumer must agree on kernel and
  padding, so change the constants, never one construction site. Exceptions: `memoize4polarize` (alpha
  1, unpadded image) and the wave projector's lattice (half-width `winsz`, its envelope from the same
  kernel object).
- Stencil weights are normalized per axis to sum one (`apod_mat_*`, the polarizer, the batch projector):
  a gather is a weighted mean. The amplitude of a padded transform is restored by pf^2 or pf^3 at the
  call site, not in the kernel.
- The fast forms are the default kernel's I0 series (coefficients for KBWINSZ=1.5, KBALPHA=2): with other
  parameters `apod_fast`/`apod_fast_device` fall back to the Bessel form, `apod_fast_value_deriv` stops
  (`error stop`), `apod_kb15_a2` (the offload reconstructor) keeps the default kernel without a word, and
  `fproject_polar_batch_opt` stops beyond 3 taps.
- Deapodize at the period of the lattice the splat used: the gridding reconstructor at the native box
  (`kb_stencil_inv_envelope_1d`), PCG and the wave projector at the padded box cropped to the native
  one (`kb_stencil_centered_crop_inv_envelope_1d`), on the real-space volume (`deapodize3D_inplace`
  stops in Fourier space). A wrong period gives a radial amplitude error and no message. Class
  averages use the continuous instrument function instead (`prep2D_inv_instrfun4mul`): ruled to move
  to the stencil envelope (check a solve2D run against the current build first); copy it nowhere.
- Packed gathers read the image `cmat` layout (stored half h >= 0, k and l wrapped) and conjugate taps
  with h < 0 themselves; `loc` is in padded logical units. The wrap table (cyclic index over the padded
  limits plus the halo, `cyci_1d` on `loop_lims(3)`) is passed with its lower bound `wrap_lower`.
  `gather_packed_window` checks nothing (a tap off the table reads out of bounds); the B-spline gather
  returns `inside=.false.` and zeros, and its caller decides (cftc stops).
- The three-tap KB stencil jumps where `nint(loc)` changes (`test_stencil_switch`); a gather whose
  derivative must be continuous (cftc's LM) uses the quadratic B-spline, whose on-grid response
  (1/8, 3/4, 1/8: `bspline3_stencil_envelope_1d`) cftc divides out through `ref_weight`. Change the two
  together.
- `center_embed_real3d`, `center_crop_real3d` and the centred-crop envelope `error stop` (not
  THROW_HARD) on a non-cubic array or an odd size difference.
- `kbinterpol` is a value type with pure methods: safe in any thread. `prep2D_inv_instrfun4mul` and
  `prep3D_inv_kbenvelope4mul` return an `image` by value (a buffer and plans per call): call them once,
  outside parallel regions.

## Couplings

- KB constants: every user above, `KB_BETA_KB15_A2`, the polynomial coefficients, the kbinterpol
  tester's closed-form values and `doc/policies/KB_Interpolation_Policy.md`.
- Gather signatures: `volume/simple_pcg_lattice.f90`, `volume/simple_reconstructor_pcg.f90`,
  `cftc/simple_cartft_calc.f90`, `strategies/search/simple_matcher_ptcl_io.f90`
  (`extract_native_fourier_plane`).

## Tests

- `build/bin/simple_test_exec test=unit_numerics suite=kaiser_bessel_kernel` (the kernel against
  closed-form values; the area about 2 s) and `test=unit_cart_align3D suite=cartesian_fourier`
  (embed/crop, envelopes, packed gathers, the B-spline oracle, the stencil switch; the area about 5 s).
- Consumers: `unit_reconstruction` (`suite=rec3d_backend`, `weighted_insertion`, `pcg_lattice`),
  `unit_pftc_align2D3D suite=polar_correlation`, `unit_cart_align3D suite=cart_calculator`; library
  `lib_reconstruction` and `pcg_recon`; `test=simulate_particles` (11 s, reprojection);
  `test=simulated_workflow` for maps (a robustness profile: one run shows breakage, not improvement).

## See also

`doc/implementation_notes/completed/drop_legacy_box_division.md` (the deapodization envelope);
skills `simple-main-image`, `simple-main-volume`, `simple-main-cftc`.
