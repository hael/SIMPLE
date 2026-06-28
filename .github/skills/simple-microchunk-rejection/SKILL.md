---
name: simple-microchunk-rejection
description: Use when working on SIMPLE's streaming microchunk class-average rejection, simple_microchunked2D lifecycle sentinels, pass-1/pass-2/refchunk/match generation, skip_pass_2 behavior, simple_cluster2D_rejector criteria, or boundaries with model_cavgs_rejection.
---

# SIMPLE Microchunk Rejection

Use this skill for `src/main/stream/simple_microchunked2D.f90`, `src/main/stream/simple_cluster2D_rejector.f90`, and stream integration questions involving `model_cavgs_rejection`.

## Read First

- `doc/microchunk_and_rejection/microchunk_rejection_policy.md`
- `doc/microchunk_and_rejection/microchunk_rejection_model_integration.md`
- `doc/microchunk_and_rejection/model_cavgs_rejection.md` when comparing the rule engine with the learned model backend.

## Live Ownership

`simple_microchunked2D` owns stream lifecycle, sentinel reconstruction, tier generation, match finalization, and combined outputs. `simple_cluster2D_rejector` owns the live scalar rejection rules. `src/main/cavg_quality` and `model_cavgs_rejection` are separate unless a task explicitly integrates them.

## Sentinel Contract

Restart state is reconstructed from sentinel files:

- `ABINITIO2D_FINISHED`: chunk job finished.
- `REJECTION_FINISHED`: rejection and particle cleanup finished.
- `REJECTION_FAILED`: rejection could not be performed safely.
- `COMPLETE`: chunk was consumed or finalized.

On import, failed chunks count as rejection-complete and complete. Do not consume a source chunk into the next tier until its project file for the next tier has been written successfully.

## Rejection Rules

The live rejector applies cumulative class-average rejection in this order:

1. population, with tier-specific fractions;
2. resolution, rejecting `res > 40.0`;
3. foreground mask geometry after low-pass/Otsu/connected-component cleanup;
4. local variance, including unconditional rejection of both-zero scores.

After rejection, rejected classes map to particle `state=0`, and rejected particles must have `class=0` and `class_match=0`.

## Tier Flow

The normal ladder is pass-1 -> pass-2 -> refchunk -> match. Pass-1 chunks become complete only after pass-2 projects are written. Pass-2 chunks feed the refchunk but are marked complete only after match chunks are generated. Match finalization updates accepted/rejected counters, writes `COMPLETE`, and copies project files to the completion directory.

When `skip_pass_2` is active, pass-2 import and generation are skipped. Rejection-complete pass-1 chunks feed the refchunk directly, and remaining pass-1 chunks become match chunks directly after the reference stack and box are available.

## Model Backend Boundary

`model_cavgs_rejection` is a learned feature-vector backend with hard validity gates, normalized features, model weights, clustering, and threshold controls. It currently exposes `chunk_default_v2`, `chunk_lp4`, `pool_default_v2`, `microchunk_p1`, and `microchunk_p2` presets. Do not replace stream rejection decisions with model output unless the task also preserves stream sentinels, selected/rejected stacks, project state propagation, and match/finalization semantics.
