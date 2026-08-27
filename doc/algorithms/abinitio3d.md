# Staged Ab Initio 3D Reconstruction

## Problem

Infer initial particle orientations, translations, state labels, and 3D maps
without a trusted high-resolution reference. The workflow must delay fine
frequencies, preserve physical scale across crops, and separate stage policy
from the base refinement estimator.

## Controller and inherited estimator

`abinitio3D` prepares project state and emits one concrete `refine3D` command
per stage. It does not implement another matcher or reconstructor. Every stage
inherits reference preparation, probabilistic pre-alignment, hard assignment,
partial reconstruction, FSC, and volume assembly from base `refine3D`.

Fresh particle starts clear prior 3D alignment while preserving transferred 2D
shifts. Selected particles become active; orientations are randomized unless a
supported input-volume or class-average initialization supplies them.

## Frequency and crop schedule

The full particle workflow has up to eight stages. Each stage independently
receives a matching low-pass limit, crop, sampling distance, translation limit,
search mode, sampling target, and iteration cap. Stage 1 is never finer than
20 A. Later values march toward `lpstop` while satisfying

```text
box * smpd = box_crop * smpd_crop.
```

When one downscaled particle cache is requested, all eligible stages use the
final planned crop so the cache key is stable; the low-pass limit still changes
per stage.

## Search schedule

Single-state and docked runs start with stochastic probabilistic-neighborhood
search (`prob_neigh_mode=shc`), enter full probabilistic search in the middle,
and use local probabilistic neighborhoods late. Independent multi-state starts
use direct stochastic hill climbing for the first two stages, then `prob`, and
only enter neighborhood stages when the user extends the default schedule.

Early stages use fewer projection directions (`nspace=500`, then 1000) and
coarser bands. Optional symmetry-axis search raises or lowers from `pgrp_start`
to the target `pgrp`; in multi-state runs each state finds and applies its own
axis transform.

## Fractional update

The controller converts a fixed `nsample` target to an outer update fraction,
capped by the workflow maximum. If `nsample/N_active > 0.9`, the global
full-sampling switch removes fractional controls and disables trailing.
Otherwise each stage records the realized subset and volume assembly trails
state/half accumulators where policy enables it.

## Multi-volume modes

- `single` evolves one state.
- `independent` initializes and refines all states independently. Its default
  five-stage, 6 A endpoint is inspection-first: it stops before late
  neighborhood, NU, trailing, and automask stages but still writes final maps.
- `docked` builds one shared model through stage 5, then creates a balanced
  persistent cohort, randomizes it into requested states, reconstructs split
  halfmaps, and refines those states with shared geometric neighborhoods.

The docked split starts a new `sampled/updatecnt` epoch. Ordinary post-split
fractional sampling remains inside the explicit cohort; a separate terminal
coverage path handles particles outside it before final reconstruction.

## Filtering and final reconstruction

Late stages may enable static nonuniform filtering and, when requested,
lag-by-one NU-evidence reference masking. `abinitio3D` does not enable NU
high-resolution shell extension. Its FSC is diagnostic rather than a claim of
gold-standard independence because current ab initio halves are not fully
independent.

At the end of a full or default independent schedule, a fresh all-particle
reconstruction runs at original sampling. It does not inherit staged search or
trailing controls. Docked mode first proves every active particle received a
post-split update.

## Implementation

- Orchestration: `src/main/commanders/simple/simple_commanders_abinitio.f90`.
- Stage policy: `src/main/simple_abinitio_controller.f90`.
- Base estimator: `src/main/commanders/simple/simple_commanders_refine3D.f90`.
- Policy: `doc/policies/abinitio3D_policy.md`.

