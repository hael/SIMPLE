---
name: simple-main-pick
description: Use when working in SIMPLE's src/main/pick subsystem, including segmentation-based (diameter) picking, batch picking with one diameter-bin and box decision, reference-based picking and its correlation kernel, the picker iterator, and picker utilities used by particle-extraction and stream workflows.
---

# SIMPLE `src/main/pick`

This folder owns picking logic.

## Read First

- `doc/algorithms/particle_picking.md`; the stream's picker is described in `doc/policies/stream/reference_generation_policy.md`, and `doc/policies/segdiampick_mics_policy.md` covers `segdiampick_mics`, the `mini_stream` validation program's estimate
- `simple_picksegdiam.f90`: segmentation-based picking with diameter estimation (the stream's version)
- `simple_segdiam_bin_picker.f90`: picks micrograph batches with one diameter-bin selection and box size, decided by the first batch
- `simple_pickref.f90`, `simple_pickref_corr_batch.f90`: reference-based picking and its BLAS-backed correlation
- `simple_picker_iter.f90`
- `simple_picker_utils.f90`: reference, segmentation and Gaussian picking entry points

## Connections

- Driven by `strategies/parallelization/simple_pick_strategy.f90` and `commanders/simple/simple_commanders_pick.f90`; in the stream, p03 picks in process (`simple_stream_stage_initial_analysis`, `segdiam_bin_picker`), p04 through `pick_extract` jobs; the stream's diameter bins and box come from `segdiam_bin_picker`, and `simple_mini_stream_utils` (`segdiampick_mics`) serves only the `mini_stream` validation program
- Works with `project/`, `image/` (`image_bin`), and preprocessing

## Working Rule

Separate reference generation, segmentation policy, and extraction/update plumbing when tracing behavior.
