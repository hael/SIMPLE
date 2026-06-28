---
name: simple-abinitio2d
description: Use when working on SIMPLE's abinitio2D workflow and its cluster2D stages, including 2D search strategies, probabilistic 2D sampling, class-average restoration, distributed partial class-average assembly, online matcher particle I/O, project updates, and streaming/rejection helpers.
---

# SIMPLE Abinitio2D

Use this skill for the orchestrating `abinitio2D` workflow and the `cluster2D`
stages it drives. Treat `cluster2D` as the core stage implementation, not as a
separate top-level mental model.

## Read First

1. `doc/policies/abinitio2D_policy.md`
2. `src/main/ui/simple/simple_ui_cluster2D.f90`
3. `src/main/exec/simple_exec_cluster2D.f90`
4. `src/main/commanders/simple/simple_commanders_cluster2D.f90`
5. `src/main/strategies/parallelization/simple_cluster2D_strategy.f90`
6. `src/main/strategies/search/simple_strategy2D_matcher.f90`
7. `src/main/class/simple_classaverager.f90`
8. `src/main/class/simple_class_frcs.f90`

If fractional update or probabilistic sampling is involved, also read
`doc/policies/importance_sampling_fractional_update_policy.md`.

## Ownership Map

- `abinitio2D`: orchestrating workflow policy and staged `cluster2D` calls.
- `cluster2D`: 2D classification stage implementation.
- `simple_strategy2D_matcher`: online particle alignment, orientation/class updates, and worker partial sums.
- `simple_matcher_ptcl_batch`: batch image loading and preprocessing.
- `class/` and `simple_commanders_mkcavgs.f90`: class-average products and explicit assembly stages.
- `project/` and `ori/`: particle/class metadata persistence.
- `stream/`: online and mini-stream 2D variants.

## Working Rules

- Keep command/controller code orchestration-focused and matcher code focused
  on particle-domain alignment and online partial-sum production.
- Preserve the probabilistic sample-once-and-reuse contract:
  `prob_align2D` chooses the sampled subset, and `prob_tab2D` plus the matcher
  reproduce that subset rather than resampling independently.
- Preserve `sampled` as the current-round marker and `updatecnt` as cumulative
  update history when 2D workflows share the fractional-update machinery.
- Preserve the online single-read particle I/O contract. The 2D matcher should
  read each particle batch once, keep the already-read raw images for
  class-average restoration, and restore/update class-average sums from those
  images after assignment in the same batch.
- Treat explicit offline or terminal class-average assembly commands as
  separate workflow stages. They may perform their own reads, but that does not
  relax the matcher worker's online single-read contract.
- For workflows that take a project file, set `mkdir=yes` before
  `params%new(cline)` unless the command is already a child command. Child
  workflow calls from inside the execution directory should use `mkdir=no` and
  the copied local project file.
- Treat previous 2D `eo` values as input project state. Do not regenerate them
  with `partition_eo` in a workflow whose purpose is to consume prior 2D class
  assignments; validate `eo` with `get_eo()` equal to `0` or `1`.

## Nearby Skills

- For stream microchunk rejection, read `.github/skills/simple-microchunk-rejection/SKILL.md`.
- For shared fractional-update contracts, read `.github/skills/simple-frac-update-trailing/SKILL.md`.

## Common Traps

- `cluster2D` logic is spread across many `simple_strategy2D_*` modules.
- Streaming `cluster2D` code is related but not identical to batch execution.
- Class averaging, selection, rejection, and project updates often live in
  different files than the search loop itself.
