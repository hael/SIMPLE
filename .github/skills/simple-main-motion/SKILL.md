---
name: simple-main-motion
description: Use when working in SIMPLE's src/main/motion subsystem, including hybrid global/local motion correction, patched alignment, iterative correction, nano-specific motion paths, and helper utilities for movie alignment workflows.
---

# SIMPLE `src/main/motion`

This folder owns movie motion correction.

## Read First

- `doc/algorithms/motion_correction.md`
- `simple_motion_correct.f90`
- `simple_motion_correct_iter.f90`
- `simple_motion_align_hybrid.f90`
- `simple_motion_patched.f90`
- `simple_motion_correct_utils.f90`

## Connections

- Commonly interacts with `image/`, `ctf/`, `project/`, `stream/`, and preprocessing commanders/strategies

## Working Rule

Keep the global-alignment stage, patch-alignment stage, and deformation-model policy conceptually separate.
