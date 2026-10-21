---
name: simple-main-motion
description: Use when working in SIMPLE's src/main/motion subsystem, including hybrid global/local motion correction, patched alignment, the polynomial motion model, iterative correction, gain-reference analysis (gain flip), movie denoising, nano-specific motion paths, and helper utilities for movie alignment workflows.
---

# SIMPLE `src/main/motion`

This folder owns movie motion correction.

## Read First

- `doc/algorithms/motion_correction.md`, `doc/policies/motion_correction_policy.md`, `doc/policies/motion_gain_analysis_policy.md`
- `simple_motion_correct.f90`
- `simple_motion_correct_iter.f90`
- `simple_motion_align_hybrid.f90` (global registration on `image/simple_ftexp_shsrch`)
- `simple_motion_patched.f90`
- `simple_motion_model.f90`: the polynomial motion model (`refine_motion_model`, particle extraction)
- `simple_motion_gain_analysis.f90`, `simple_motion_gain_helpers.f90`: gain-flip analysis
- `simple_motion_correct_utils.f90`

## Connections

- Driven by `strategies/parallelization` (`simple_motion_correct_strategy.f90`, `simple_preprocess_strategy.f90`, `simple_refine_motion_model_strategy.f90`) and the stream's preprocessing stage (gain flip)
- `simple_motion_align_nano.f90` serves the SINGLE time-series path (`nano/single_tseries_tracker.f90`)
- Commonly interacts with `image/`, `ctf/`, `project/`, `preprocess/`, and `stream/`
- Tests: `simple_test_exec test=unit_motion`

## Working Rule

Keep the global-alignment stage, patch-alignment stage, and deformation-model policy conceptually separate. Alignment references come from the movie itself, never from a volume or particle reference.
