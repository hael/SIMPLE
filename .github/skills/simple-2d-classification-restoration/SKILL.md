---
name: simple-2d-classification-restoration
description: Use when working on SIMPLE's abinitio2D/cluster2D classification, probabilistic 2D sampling, class-average restoration, distributed partial class-average assembly, and online matcher particle I/O.
---

# SIMPLE 2D Classification and Restoration

Use this skill for SIMPLE's 2D particle-classification and class-average
restoration workflow, especially when changes touch `abinitio2D`, `cluster2D`,
probabilistic 2D alignment, distributed partial class-average sums, or online
matcher image I/O.

## Quick Start

1. Read `doc/policies/abinitio2D_policy.md` for the current public workflow policy.
2. If fractional update or probabilistic sampling is involved, also read
   `doc/policies/importance_sampling_fractional_update_policy.md`.
3. Inspect the owner layer before editing:
   - `simple_commanders_cluster2D.f90` and ab initio 2D command/controller code
     for orchestration and stage policy
   - `simple_strategy2D_matcher.f90` for online particle alignment,
     orientation/class updates, and worker partial sums
   - `simple_matcher_ptcl_batch.f90` for batch image loading and preprocessing
   - classaverager modules and `simple_commanders_mkcavgs.f90` for explicit
     class-average assembly

## Working Rules

- Keep command/controller code orchestration-focused and matcher code focused
  on particle-domain alignment and online partial-sum production.
- Preserve the probabilistic sample-once-and-reuse contract:
  `prob_align2D` chooses the sampled subset, and `prob_tab2D` plus the matcher
  reproduce that subset rather than resampling independently.
- Preserve `sampled` as the current-round marker and `updatecnt` as cumulative
  update history when 2D workflows share the fractional-update machinery.
- Preserve the online single-read particle I/O contract.
  The 2D matcher should read each particle batch once, keep the already-read raw
  images for class-average restoration, and restore/update class-average sums
  from those images after assignment in the same batch.
- Do not split online class-average restoration into a second full particle pass
  that re-reads image stacks to lower peak memory unless the repository policy
  explicitly changes.
- Treat explicit offline or terminal class-average assembly commands as
  separate workflow stages. They may perform their own reads, but that does not
  relax the matcher worker's online single-read contract.

## SIMPLE Project And Even/Odd Conventions

- For user-facing SIMPLE workflows that take a project file, follow the standard
  execution-directory path: set `mkdir=yes` before `params%new(cline)` unless
  the command is already a child command, let `params%new` create and enter the
  execution directory, and modify the copied local `params%projfile`.
- Keep the copied project filename stable. The execution directory identifies
  the producing workflow; do not invent workflow-specific project basenames
  unless an existing parent/temporary-project pattern explicitly requires one.
- For child workflow calls from inside that execution directory, set `mkdir=no`
  and pass the copied local `params%projfile` to avoid nested directories and
  stale parent-project paths.
- Treat `eo` from previous 2D workflows as input project state. Do not
  regenerate it with `partition_eo` in a workflow whose purpose is to consume
  prior 2D class assignments; fail clearly if the copied project lacks valid
  `eo`.
- Validate particle `eo` with `get_eo()` equal to `0` or `1`. Do not use
  `isthere('eo')` as the validity test, because `eo=0` is a valid even
  assignment but can be reported as absent for particle orientation parameters.

## Language To Prefer

- "online class-average restoration reuses matcher batch images"
- "probabilistic 2D pre-alignment followed by hard assignment"
- "distributed partial class-average sums"
- "sample once, then reproduce the same subset"
- "offline assembly is separate from online matcher restoration"
