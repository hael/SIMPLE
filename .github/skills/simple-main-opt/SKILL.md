---
name: simple-main-opt
description: Use when working in SIMPLE's src/main/opt subsystem, including the optimizer abstraction, the factory and spec, L-BFGS-B, Nelder-Mead simplex, differential evolution, and the rank-1 preconditioned conjugate-gradient engine used by the PCG reconstructions.
---

# SIMPLE `src/main/opt`

This folder owns general optimization infrastructure.

## Read First

- `simple_optimizer.f90`: the abstract `optimizer`
- `simple_opt_factory.f90`: `simplex`, `de` or `lbfgsb` from `opt_spec%str_opt`
- `simple_opt_spec.f90`: bounds, starting point, cost and gradient callbacks
- `simple_opt_lbfgsb.f90`
- `simple_opt_simplex.f90`, `simple_opt_subs.f90` (the amoeba step)
- `simple_opt_de.f90`
- `simple_pcg_solver.f90`: rank-1 PCG engine with client-supplied operator, preconditioner and inner product (`volume/simple_reconstructor_pcg`, `simple_pcg_atomic_prior`, `flex/fit/simple_flex_pca_pcg`)

## Connections

- Used from the shift searches (`pftc/simple_pftc_shsrch_grad`, `image/simple_ftexp_shsrch`), patched motion correction, CTF fitting, symmetry search (`volume/simple_volpft_symsrch`), and class-average quality learning
- Tests: `unit_numerics` ("optimisers")

## Working Rule

When debugging optimizer behavior, distinguish optimizer mechanics from the objective/gradient supplier living elsewhere.
