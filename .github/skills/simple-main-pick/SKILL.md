---
name: simple-main-pick
description: Use when working in SIMPLE's src/main/pick subsystem, including segmentation-based picking, picking references, iterative picking helpers, and picker utilities used by particle-extraction workflows.
---

# SIMPLE `src/main/pick`

This folder owns picking logic.

## Read First

- `simple_pickseg.f90`
- `simple_picksegdiam.f90`
- `simple_pickref.f90`
- `simple_picker_iter.f90`
- `simple_picker_utils.f90`

## Connections

- Works with `project/`, `image/`, preprocessing, and stream workflows

## Working Rule

Separate reference generation, segmentation policy, and extraction/update plumbing when tracing behavior.
