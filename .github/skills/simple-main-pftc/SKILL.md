---
name: simple-main-pftc
description: Use when working in SIMPLE's src/main/pftc subsystem, including the polarft_calc type (polar Fourier references and particles, memoised FFT terms, cc and Euclidean objectives and their shift gradients, the rotation samplers of the probabilistic searches, per-thread work buffers, reference I/O) and the gradient shift search pftc_shsrch_grad used by the 2D and 3D alignment workflows.
---

# SIMPLE `src/main/pftc`

This folder is central to polar-Fourier search workflows.

## Read First

- `simple_polarft_calc.f90`: the `polarft_calc` type, the per-thread `heap_vars` type and every interface
- `simple_polarft_core.f90`: `new`/`kill`, setters, the per-thread allocations
- `simple_polarft_corr.f90`: objectives, gradients, the sampler front ends, the `_aln` evaluators
- `simple_polarft_geom.f90`: shift matrices (`gen_shmat*`), rotation, mirroring
- `simple_polarft_memo.f90`: memoised particle/reference FFTs, square sums, sigma2 weights, the shifted-particle bank
- `simple_polarft_access.f90`, `simple_polarft_ctf.f90`, `simple_polarft_ops_io.f90`
- `simple_pftc_shsrch_grad.f90`: the L-BFGS-B shift (or joint shift and in-plane) search of one reference/particle pair

## Structure

- The type and its interfaces are declared in `simple_polarft_calc.f90`; the bodies live in its submodules
- Mutable work arrays are per OpenMP thread, indexed by `omp_get_thread_num() + 1`: `heap_vars(ithr)` (shift matrices, double-precision reference copies, rotation scores, distances and cdf; `params%nthr` of them, from `new`) and the FFTW buffers (`cmat2_many`, `cmat_joint_many`, `crmat_many`, `crmat1_many`, `crvec1`, `drvec`; `nthr_glob` of them, from the memo workspace). Keep the thread index local to the call; never share these buffers across threads
- References are stored per half (`pfts_refs(:,:,iref,ieo)`), picked by the particle's `iseven`; sigma2 is two global spectra, one per half (`assign_sigma2_noise`, `memoize_sigma2_weights`)
- The rotation samplers (`sample_bounded_dists`, `sample_likelihood_dists`, `sample_power_dists` in `strategies/search/probabilistic/simple_eul_prob_tab_utils.f90`) take arrays of distances filled in the thread's buffers; do not reintroduce internal procedures (heap trampolines) or per-call allocations in the score loop
- The alignment shift matrix (`gen_shmat4aln_8`) spans only the searched shells `kfromto(1):kfromto(2)`; the single-precision one (`gen_shmat4aln`) is converted from it on those shells only

## Load-Before-Evaluate Contract

- A shift search copies its reference once per minimisation into the calling thread's buffer: `load_ref_aln_8(iref, iptcl)` records the reference index and half. `gen_corr_for_rot_8_aln` and `gen_corr_grad_for_rot_8_aln` evaluate against that copy and stop (`THROW_HARD`) on an index or half mismatch
- `pftc_shsrch_grad%minimize` does the load; a new caller of the `_aln` evaluators must load first, on the same thread
- The loaded state lives for a whole minimisation, so a shift search runs on one thread: no nested parallel region or untied task inside it

## Connections

- Polar sections come from `image/` (`simple_image_polar`, `simple_projector`, `simple_projector_pft_batch`); CTF from `ctf/`
- Driven by `strategies/search` (2D/3D matchers, strategy objects, probability tables) through the builder (`build%pftc`)
- Tests: `simple_test_exec test=unit_pftc_align2D3D`. New OpenMP code here also needs repeated runs at many threads, not only one-thread identity
- Policy: `doc/policies/3D/refine3D_policy.md` (shift search of the tables, continuous in-plane), `doc/policies/2D/solve2D_policy.md`

## Working Rule

Changes here often affect alignment performance and correctness simultaneously, so treat caching, per-thread buffers and geometry logic as part of the design, not incidental detail.
