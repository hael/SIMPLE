---
name: simple-main-cftc
description: Use when working in SIMPLE's src/main/cftc subsystem, the Cartesian Fourier calculator of continuous pose refinement (refine=cont, the pose_cont polish), including cartft_calc (padded reference transforms per state and half, ref_weight, particle slots of a batch, sigma2 whitening, the B-spline gather, predict, objective, gradient and Gauss-Newton normal terms, shell residuals, the reference-volume files) and cartft_pose_opt (the bounded Levenberg-Marquardt transaction of one particle's five-parameter pose, CARTFT_* statuses, the work record).
---

# `src/main/cftc`: the Cartesian Fourier calculator

Owns `cartft_calc` (`build%cftc`), the Cartesian twin of `pftc`: the padded reference transforms of
every state and half, the prepared particle slots of one batch, and a pose's objective with its
gradient and Gauss-Newton terms; and `cartft_pose_opt`, the bounded LM transaction of one particle's
pose (three right-increment rotation components, two shifts). It scores and solves one particle; it
does not choose particles or states, prepare references or particles, write scores, the project or
sigma2. Those are `strategies/search`: references in `simple_matcher_refvol_utils.f90`, slots in
`prep_cart_batch` (`simple_matcher_ptcl_batch.f90`) with `prepimg4align_cart`
(`simple_matcher_2Dprep.f90`), the transaction and its statuses in `simple_strategy3D_cont.f90`; the
shell residual becomes sigma2 in `sigma2/simple_euclid_sigma2.f90`. Contract: refine3D_policy.md,
"Continuous Cartesian pose policy", and KB policy section 6a.

## Rules

Silent ones first: nothing reports these.

- Order in a pass: `read` (`new`, then `set_ref` per state and half from `cart_refvols_{even,odd}.bin`;
  returns the header's kfromto, which the pass adopts) -> `new_ptcls(batchsz_max)` -> per batch:
  `set_sigma2` (euclid) -> `set_ptcl` per slot -> `set_ptcl_inds(pinds)` -> per particle
  `get_ptcl_slot`, `cartft_pose_opt%refine` -> `kill`. Materializer: `new(nstates, box_crop, 1)` ->
  `set_refvol` per state and half -> `write` -> `kill`.
- A Euclidean slot is whitened by the spectra it was prepared with: `set_sigma2` marks every Euclidean
  slot stale, and any use of a stale slot stops (`ptcl_is_valid`, so also the pose transaction and the
  sigma contribution), so prepare slots after setting the spectra; cc slots stay. A degenerate particle
  (no power) is still an invalid slot without a stop. Euclid without any spectra stops (THROW_HARD).
- `iptcl` of `set_ptcl`, the evaluators and `refine` is a slot (1..nptcls), not a particle index; map it
  with `set_ptcl_inds`/`get_ptcl_slot`. `new_ptcls` drops the map, and a slot of 0 or an invalid slot
  makes `refine` return `CARTFT_INVALID_PREPARATION`, not an error.
- The observation is the polar particle up to the soft mask (noise-normalized, Fourier-clipped, shifted
  by minus the stored shift, phase-flipped for CTFFLAG_YES), times `get_ptcl_taper` in real space
  before its transform (`prepimg4align_cart`). Without the taper `ref_weight` no longer matches and the
  objective is biased; a shift added in `set_ptcl` is applied twice. A change to `prepimg4align` is
  mirrored there, or the two representations stop sharing one sigma2.
- Units and signs: shifts in cropped pixels, the model shift an increment on the stored one (`refine`
  converts it by box/box_crop); phase exp(+2 pi i (h s1 + k s2)/box), image's sign; rotation R exp([w]x),
  updated on the right; sample at padf*([h,k,0] R), as `fproject`; a pixel is in shell nint(|(h,k)|)
  inside |(h,k)| <= box/2, the polar ring membership.
- Halves: `refs(1,:)` even, `refs(2,:)` odd, picked by `iseven` (eo 0); the sigma2 column is eo+1
  (`sigma2_half`). Coefficients carry padf^3 and `ref_weight` (KB stencil over B-spline response at the
  padded period); KB constants change it and the taper together (`simple-main-interp`).
- Reference file: the prepared real-space volumes after a header [format_version, box_crop, nstates,
  kfrom, kto] and smpd_crop; writer, reader and `cart_refvols_header_compatible` change together, with
  CART_REFVOLS_FORMAT_VERSION bumped.
- Threads: the evaluators take the calculator `intent(in)` and run concurrently on different slots;
  `set_ptcl` writes only its slot (the batch fills slots in a parallel loop); `new`, `set_ref`,
  `set_sigma2`, `new_ptcls`, `set_ptcl_inds`, `read` and `write` are serial. One `cartft_pose_opt` per
  particle search object; `cartft_pose_work` is summed per thread.
- Failures: a slot without power or with an unusable sigma2 shell is invalid (`ptcl_is_valid`), a
  per-particle outcome; evaluating it, an absent reference or a gather off the wrap table stops. A solve
  that does not lower the objective returns the seed bit for bit with a typed status.

## Where to start

- `simple_cartft_calc.f90`: types, lifecycle, `fill_ptcl`, `pose_normal_terms` (the half-plane production
  kernel); the `*_full_disk` kernels are the tests' reference.
- `simple_cartft_pose_opt.f90`: statuses, `refine` -> `refine_pose` -> `refine_joint`; `refine_shift`
  runs only on the `cont_route=shift_then_joint` test route (an `internal:` key, ruled to be retired).

## Tests

- `build/bin/simple_test_exec test=unit_cart_align3D` (about 5 s; `suite=` cartesian_fourier,
  cart_calculator, pose_optimizer, pose_strategy, pose_statistics and the local-start suites).
- `test=lib_cart_align3D` (about 30 s): 1 000 simulated 1JYX particles from perturbed poses, cc and
  euclid at 8 and 4 A, gated on the truth.
- Behaviour change: `test=cont_refine3D_1jxy` (highlevel, a median 390 s at 8 threads; median rotation
  error at most 1.5 degrees, shift 0.2 px, 95 % in the basin; its profile passes every run, so one
  failure is a regression). Speed: `test=cont_refine3D_bench` (timings, work counts, no thresholds).

## See also

`doc/refactoring_notes/completed/pose_cont_refactoring.md` (the C*/O* contract items cited in comments),
`doc/implementation_notes/completed/cont_lm_acceleration.md`; skills `simple-main-pftc`,
`simple-main-interp`, `simple-refine3d`.
