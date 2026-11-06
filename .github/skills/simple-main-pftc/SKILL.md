---
name: simple-main-pftc
description: Use when working in SIMPLE's src/main/pftc subsystem, including the polarft_calc type (polar Fourier references and particles, memoised FFT terms, cc and Euclidean objectives and their shift gradients, the rotation samplers and keyed posterior draws of the probabilistic searches, per-thread work buffers, reference I/O), the polar shift maps, and the gradient shift search pftc_shsrch_grad used by the 2D and 3D alignment workflows.
---

# `src/main/pftc`: the polar Fourier calculator

Owns `polarft_calc` (`build%pftc`): polar references and particles, their memoised transforms, the
cc and Euclidean objectives with their shift gradients, the in-plane draws, the polar shift maps and
the shift search of one pair (`pftc_shsrch_grad`). It scores pairs; it does not decide which
references a particle sees, how particles are sampled or what goes into a probability table: that
is `strategies/search`. Sigma2 is estimated in `sigma2/`; polar sections are cut by `image/`
(`simple_image_polar`, `simple_projector*`).

## Rules

Silent ones first: nothing reports these.

- Memoise before scoring, outside any parallel region. Order: `new` -> references
  (`set_ref_pft`/`polarize_ref_pft` per half, or even then `cp_even2odd_ref`) -> `memoize_refs` ->
  particles (`reallocate_ptcls`, `set_ptcl_pft`/`polarize_ptcl_pft`) -> `memoize_ptcls` -> score ->
  `kill`. Every particle write bumps a generation that `memoize_ptcls` records; the scorers that read
  the memo stop on a stale one (`check_ptcl_memo`; the `_for_rot_8` kernels read the particles, not the
  memo, and do not check). `memoize_*` inside an OpenMP region stops (`THROW_HARD`).
- When the sigma2 spectra change, call `memoize_sigma2_weights` before scoring (`prep_sigmas_objfun`
  does after `read_global`): the particle self-term `wsqsums_ptcls` is refreshed only there. The
  Euclidean scorers stop without weights or on spectra newer than them (`check_euclid_weights`: the
  owner's `spectra_gen`, given to `assign_sigma2_noise` by `esig%new`). A
  normalisation that restores the residual by the "+1" identity must use `wsqsum_ptcl_now`, or the
  loss is unbounded below and exp(-loss) scores exceed 1.
- Mutable work arrays are per OpenMP thread, indexed by `omp_get_thread_num() + 1` and sized by the team
  in use at `new` (`omp_get_max_threads()`): `heap_vars(ithr)`, the FFTW buffers `cmat2_many`,
  `cmat_joint_many`, `crmat_many`, `crmat1_many`, `crvec1`, `drvec`, and the shift-map work `smwork`
  (sized like `crvec1`; a team that grew after `new` stops). Keep the thread index local to the call. New scratch goes into these buffers,
  never into a module variable, a per-call allocation in the score loop or an internal procedure
  (heap trampoline).
- Load before evaluate, on one thread. `load_ref_aln_8(iref, iptcl)` copies the reference into the
  thread's buffer; `gen_corr_for_rot_8_aln` and `gen_corr_grad_for_rot_8_aln` evaluate against that
  copy and stop on an index or half mismatch; `pftc_shsrch_grad%minimize` does the load. Shift maps
  likewise: `gen_shiftmap`, then `shiftmap_cross`, `shiftmap_objfun_grid` or `shiftmap_draw` on the
  same thread. The loaded state spans a whole minimisation: no nested parallel region or untied task
  inside it.
- References are stored per half (`pfts_refs(:,:,iref,ieo)`) and picked by the particle's `iseven`;
  sigma2 is two spectra, one per half.
- Shift matrices span only the searched shells `kfromto(1):kfromto(2)` (`gen_shmat4aln_8`; the
  single-precision `gen_shmat4aln` is converted from it on those shells only).
- The probability tables score a reference with the posterior in-plane draw (`gen_posterior_val`,
  `gen_posterior_val_bank`: the Euclidean row over every rotation, `draw_posterior_row` at
  `ptcl_posterior_temp`, objfun=euclid only). Its variate is keyed on the calculator's scope, the
  particle and the reference (`RND_FAM_INPL`), so a pair draws the same rotation on any thread, batch or
  part. The table sets the scope in each fill (`set_rnd_scope(prob_rnd_scope(...))`), `kill` drops it,
  and a draw without one stops. `shiftmap_draw`'s disk draw is keyed the same way (`RND_FAM_DISK`);
  production calls only `shiftmap_mode`, which draws nothing. `gen_euclid_dists(_bank)` return the row
  without the draw (bench, testers). The SHC samplers (`sample_bounded_dists`, `sample_power_dists` in
  `strategies/search/probabilistic/simple_eul_prob_tab_utils.f90`, used by the classic `*_smpl` and
  `neigh` routes) take distance arrays filled in the thread's buffers and draw from the thread's `ran3`
  stream.
- `gen_many_euclids_gpu` and `gen_many2many_euclids_cufft` compile only with `USE_OPENMP_OFFLOAD`
  (off by default) and run only on a GPU build: say so when a change touches them.

## Where to start

- `simple_polarft_calc.f90`: the type, `heap_vars` and every interface. Bodies are in submodules:
  `_core` (new/kill, setters), `_memo` (memoised transforms, sigma2 weights, shifted-particle bank),
  `_corr` (objectives, gradients, samplers, `_aln` evaluators), `_geom` (shift matrices, rotation,
  mirroring), `_shiftmap` (a pair's shift maps over every rotation), `_access`, `_ctf`, `_ops_io`.
- `simple_polar_shiftmap.f90`: one pair's shift-dependent cross term over a disk of shifts
  (Fourier-Bessel tables, driven by `_shiftmap`).
- `simple_pftc_shsrch_grad.f90`: the L-BFGS-B shift (or joint shift and in-plane) search of one pair;
  its `kill` deallocates the optimiser.
- Callers that show the order: `strategies/search/simple_matcher_pftc_prep.f90` (2D references),
  `simple_matcher_ptcl_batch.f90` (particle batches, sigma2 weights), `simple_shiftmap_srch.f90`.

## Tests

- `build/bin/simple_test_exec test=unit_pftc_align2D3D` (about 12 s; `suite=` runs one of: polar
  correlation, polar shift map, continuous in-plane, refine3D in-plane state, the 2D and 3D
  probability tables, sigma2 state and bootstrap, ...). The testers are here and in `strategies/search`.
- New OpenMP code: repeated runs at many threads, not only the one-thread identity; for the tables'
  draws `test=polar_matching_bench suite=fills` (threads, batches and parts give identical entries).
- Kernel speed or layout: `test=polar_matching_bench` (timings and identity checks, a dump of the
  production kernels' distance vectors for comparing builds).
- A change in scores moves alignment: run `test=simulated_workflow`. Its map correlation is a
  robustness profile (1-2 misses in 20 runs is normal), so one run shows breakage, not improvement.

## See also

`doc/policies/3D/refine3D_policy.md` (shift search of the tables, continuous in-plane),
`doc/policies/2D/solve2D_policy.md`, `doc/refactoring_notes/completed/polar_shift_search_performance_refactoring.md`;
skills `simple-refine3d`, `simple-solve2d`.
