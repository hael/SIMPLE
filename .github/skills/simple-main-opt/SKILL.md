---
name: simple-main-opt
description: Use when working in SIMPLE's src/main/opt subsystem, including optimizer abstractions, BFGS, L-BFGS-B, simplex, particle swarm, differential evolution, conjugate-gradient variants, and optimization factories/specs used throughout the platform.
---

# SIMPLE `src/main/opt`

This folder owns general optimization infrastructure.

## Read First

- `simple_optimizer.f90`
- `simple_opt_factory.f90`
- `simple_opt_spec.f90`
- `simple_opt_lbfgsb.f90`
- `simple_opt_bfgs.f90`
- `simple_opt_simplex.f90`

## Connections

- Used from motion correction, image alignment, search strategies, and numerical utilities

## Working Rule

When debugging optimizer behavior, distinguish optimizer mechanics from the objective/gradient supplier living elsewhere.
