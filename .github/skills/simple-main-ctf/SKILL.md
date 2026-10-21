---
name: simple-main-ctf
description: Use when working in SIMPLE's src/main/ctf subsystem, including the microscope contrast transfer function model, anisotropic/patched CTF estimation, fitting cost functions, and iterative estimation workflows.
---

# SIMPLE `src/main/ctf`

This folder owns CTF modeling and estimation.

## Read First

- `simple_ctf.f90`: the `ctf` type (the transfer function model)
- `simple_ctf_estimate_cost.f90`
- `simple_ctf_estimate_fit.f90`: patched (anisotropic) estimation
- `simple_ctf_estimate_iter.f90`: one estimation iteration per micrograph
- `doc/algorithms/ctf_estimation.md`, `doc/policies/phase_shift_ctf_policy.md`, `doc/policies/microscope_parameters_policy.md`

## Connections

- Estimation is driven by `strategies/parallelization` (`simple_ctf_estimate_strategy.f90`, `simple_preprocess_strategy.f90`, also run by the stream's preprocessing jobs)
- Applying a CTF to images lives in `image/` (`simple_image_ctf.f90`); the polar CTF matrices in `pftc/`; `project/` holds the parameters
- Tests: `unit_image` ("CTF")

## Working Rule

Separate the microscope model from estimation-loop policy when making changes; many bugs come from blending those concerns.
