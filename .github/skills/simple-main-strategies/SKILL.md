---
name: simple-main-strategies
description: Use when working in SIMPLE's src/main/strategies subsystem, including parallelization strategies, 2D and 3D search strategies, matchers, batching, reconstruction helpers, and the execution-policy layer between commanders and low-level scientific objects.
---

# SIMPLE `src/main/strategies`

This folder is the execution-policy and search-engine layer.

## Read First

- `parallelization/simple_refine3D_strategy.f90`
- `parallelization/simple_cluster2D_strategy.f90`
- `parallelization/simple_preprocess_strategy.f90`
- `search/simple_strategy2D_matcher.f90`
- `search/simple_strategy3D_matcher.f90`
- `search/simple_matcher_3Drec.f90`

## Subareas

- `parallelization/` handles shared-memory vs distributed policy and workflow iteration
- `search/` holds 2D/3D alignment strategy families, matchers, batching, and search helpers

## Working Rule

This folder is usually the right home for “how a workflow runs” but not “what the raw data structure means.”
