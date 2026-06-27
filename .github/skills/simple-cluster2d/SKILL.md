---
name: simple-cluster2d
description: Use when working on SIMPLE's cluster2D subsystem, including CLI/UI entrypoints, 2D search strategies, class averaging, class statistics, project updates, rejection/streaming helpers, and the interaction between commanders, strategies, class objects, and stream workflows.
---

# SIMPLE Cluster2D

`cluster2D` also crosses several folders. Treat it as a workflow, not a single file.

## Start Here

1. `src/main/ui/simple/simple_ui_cluster2D.f90`
2. `src/main/exec/simple_exec_cluster2D.f90`
3. `src/main/commanders/simple/simple_commanders_cluster2D.f90`
4. `src/main/strategies/parallelization/simple_cluster2D_strategy.f90`
5. `src/main/strategies/search/simple_strategy2D_matcher.f90`
6. `src/main/class/simple_classaverager.f90`
7. `src/main/class/simple_class_frcs.f90`

## Ownership Map

- CLI/UI surface: `ui`, `exec`
- Command-layer workflow: `simple_commanders_cluster2D.f90`
- Parallel/distributed policy: `strategies/parallelization/simple_cluster2D_strategy.f90`
- Core 2D matching/search behavior: `strategies/search/simple_strategy2D*.f90`
- Class products and restore paths: `class/`
- Project/orientation persistence: `project/`, `ori/`
- Streaming variants and rejection logic: `stream/`

## Nearby Folders That Matter

- `image/` for image transforms and FFT-backed operations
- `pftc/` for polar Fourier calculations used by search
- `project/` for particle/class metadata persistence
- `stream/` for online/mini-stream cluster2D paths

## Focused Companion Skill

For probabilistic 2D sampling, online class-average restoration, distributed
partial class-average sums, `sampled`/`updatecnt`, or even/odd consumption from
prior 2D workflows, read
`.github/skills/simple-2d-classification-restoration/SKILL.md`.

## Common Traps

- `cluster2D` logic is spread across many `simple_strategy2D_*` modules.
- Streaming `cluster2D` code is related but not identical to batch execution.
- Class averaging, selection, rejection, and project updates often live in different files than the search loop itself.
