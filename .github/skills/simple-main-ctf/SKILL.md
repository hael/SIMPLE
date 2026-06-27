---
name: simple-main-ctf
description: Use when working in SIMPLE's src/main/ctf subsystem, including microscope contrast transfer function models, anisotropic/patched CTF estimation, fitting cost functions, and iterative estimation workflows.
---

# SIMPLE `src/main/ctf`

This folder owns CTF modeling and estimation.

## Read First

- `simple_ctf.f90`
- `simple_ctf_estimate_cost.f90`
- `simple_ctf_estimate_fit.f90`
- `simple_ctf_estimate_iter.f90`

## Connections

- Preprocessing and motion-correction workflows often feed or consume CTF results
- `image/`, `motion/`, `project/`, and `exec/commanders` are common neighboring folders

## Working Rule

Separate the microscope model from estimation-loop policy when making changes; many bugs come from blending those concerns.
