---
name: simple-main-opt
description: Use when working in SIMPLE's src/main/opt subsystem, including the optimizer abstraction, the factory and opt_spec (cost and gradient callbacks, limits, tolerances, the fun_self context), L-BFGS-B, Nelder-Mead simplex, differential evolution, the rank-1 preconditioned conjugate-gradient engine of the PCG reconstructions, and the nonlinear conjugate-gradient engine of reconstruct3D_wave.
---

# `src/main/opt`: minimisers and the conjugate-gradient engines

Owns the minimisers behind `opt_spec` and `opt_factory` (`lbfgsb`, `simplex`, `de`) and two engines on
a client's operator: `pcg_solve` (linear PCG) and `nlcg_solve` (nonlinear CG). It owns no objective:
the cost, its gradient and the meaning and scale of `x` belong to the client (shift searches in `pftc/`
and `image/simple_ftexp_shsrch`, CTF fitting, `volume/simple_volpft_symsrch`, `quality/`; PCG:
`volume/simple_reconstructor_pcg`, `simple_pcg_atomic_prior`, `flex/fit/simple_flex_pca_pcg`; NLCG:
`strategies/parallelization/simple_rec3D_wave_strategy`). A wrong result is usually the client's cost
or gradient: check those first.

## Rules

Silent ones first: nothing reports these.

- Order: `ospec%specify(str_opt, ndim, ...)` -> cost pointers (`set_costfun`, `set_gcostfun`,
  `set_fdfcostfun` or the `_8` set) -> `opt_factory%new(ospec, ptr)` -> per solve `ospec%x` (start)
  and limits -> `ptr%minimize(ospec, fun_self, lowest_cost)`; the result is `ospec%x`. `specify` starts
  with `kill`: it nulls the cost pointers and the callback and zeroes `x`.
- Every optimiser minimises: return a cost (a negated correlation), flip `lowest_cost` back.
- `ftol`/`gtol` differ by optimiser. Simplex and DE stop when the fractional cost range is below
  `ftol`. In `lbfgsb` they are the line search's sufficient-decrease and curvature constants; it stops on
  `factr` (relative decrease in machine epsilons) and `pgtol` (projected gradient), and `maxits` counts
  cost-and-gradient evaluations, not iterations. A client states its convergence through `factr`/`pgtol`
  and leaves the line search at its standard constants (`ftol=1e-3`, `gtol=0.9`, as the CTF fit does),
  unless a measurement says otherwise: the polar shift search keeps a strict line search (`ftol=gtol=1e-5`,
  stated at its `specify`), since the standard one fails there (more evaluations, line-search failures,
  larger shift errors; plan agent_friendly, phase 4). `opt_spec`'s defaults are `ftol=gtol=1e-5`.
- `x` is the start, the origin included; `ospec%nostart = .true.` (reset by `specify`) says there is
  none: simplex then draws a random start within the limits, DE keeps its random population. `lbfgsb` starts at
  `ospec%x` (copied into `x_8`; setting only `x_8` has no effect).
- Bounds: `lbfgsb` needs `limits` (`THROW_HARD`) and keeps both; DE clamps trials to them; simplex uses
  `limits` (or `limits_init`) only to place its vertices, the amoeba step is unbounded and the result
  can leave them. `specify('simplex')` stops without limits; DE reads them unchecked.
- Precision: `lbfgsb` works in double and calls `fdfcostfun_8`, else `costfun_8` + `gcostfun_8`, else
  the single pair on `x` rounded to real. Give cost and gradient in one precision: a mixed pair passes
  the `lbfgsb` check, then every evaluation returns f = 0, g = 0 with a printed error only. Simplex and
  DE take only the single-precision `costfun`.
- DE takes F and CR from `str_mode` (`specify('de')` sets `multimodal`) and the population size from
  `specify(npop=)`, else from the mode; the start, then the members given by `set_inipop`, replace its
  first random members (all members are drawn first, so the random stream does not depend on them).
- `fun_self` is the client's object (`class(*)`), resolved in the cost by `select type` with a
  `class default` `THROW_HARD`. `opt_callback` runs at each accepted `lbfgsb` iterate, `x` updated, and
  may change the objective (the discrete in-plane angle of `pftc_shsrch_grad`).
- Threads: a spec and an optimiser hold all their state (L-BFGS-B's saved locals live in the instance):
  one pair per thread, never shared (`minimize` writes `x`, counters and workspace). Simplex and DE draw
  random numbers: results depend on the seed and the thread.
- `opt_factory%new` allocates a new optimiser per call; the caller owns the pointer: `kill`,
  `deallocate` (as `image/simple_ftexp_shsrch`, `pftc_shsrch_grad`, `ctf_estimate_cost4Dcont` and
  `volpft_symsrch` do in their `kill`); a `kill` that only nullifies leaks the optimiser, and nothing reports it.
- `pcg_solve` needs `apply` and `precond` symmetric positive definite. Non-positive curvature returns
  `stop_reason='indefinite'`, `converged=.false.`, no error (one cold retry if `cold_restart` and the
  start was warm): the caller checks. `rtol <= 0` runs `maxits`. Contract:
  `doc/policies/3D/reconstruct3D_pcg_policy.md`.
- `nlcg_solve` (contract: `doc/policies/3D/wave_model_policy.md` 6.1): `indefinite` and `no_descent`
  return without error; `line_losses` compare only with each other within one iteration, not with
  `eval`; `project` must satisfy P^2 = P. Both engines stop (`THROW_HARD`) on a size mismatch and
  when the first dot(r,z) <= 0 (preconditioner not positive definite).

## Where to start

- `simple_opt_spec.f90` (spec, callback interfaces, `eval_*`), `simple_opt_factory.f90`.
  `simple_opt_lbfgsb.f90` is a legacy port (`goto`, state in `isave`/`dsave`): don't copy its style.
- Callers that show the order: `image/simple_ftexp_shsrch.f90` (`ftexp_shsrch_new`),
  `pftc/simple_pftc_shsrch_grad.f90` (pointers set after `specify`), `volume/simple_volpft_symsrch.f90`
  (one spec and optimiser per thread).

## Tests

- `build/bin/simple_test_exec test=unit_numerics suite=optimisers` (spec bookkeeping, L-BFGS-B, DE,
  simplex, PCG with cold restart) and `suite=nonlinear_conjugate_gradient`; all of `unit_numerics`
  takes about 2 s. Clients: `suite=shift_search` (`ftexp_shsrch`), `test=unit_pftc_align2D3D`.
- PCG engine: `test=pcg_recon`, `test=lib_reconstruction` (its wave reconstruction takes about five
  minutes), `unit_reconstruction`, `unit_heterogeneity`.
- A change in `lbfgsb` or the spec moves every shift search: run `test=simulated_workflow`, a robustness
  profile (1-2 misses in 20 runs is normal), so one run shows breakage, not improvement.

## See also

Skills `simple-main-pftc`, `simple-main-volume`, `simple-modern-fortran`.
