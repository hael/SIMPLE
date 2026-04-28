# Abinitio2D Policy

## 1. Purpose and Scope

This document defines the current architectural policy for `abinitio2D` and the `cluster2D` workflow it drives.

It mirrors the recent `refine3D` cleanup where the same design pressure exists:

- keep stage policy separate from execution mechanics
- keep particle/class assignment work separate from class-average assembly
- preserve shared-memory and distributed parity
- make sampled-update and probabilistic handoffs explicit
- treat class-average files, assignment files, FRCs, and partial sums as workflow contracts rather than incidental scratch files

The 2D workflow is intentionally Cartesian. `polar=yes` is not supported for `abinitio2D` or `cluster2D`.

## 2. Architectural Policy

`abinitio2D` is a staged 2D classification workflow:

1. set run defaults and read project state
2. determine stage geometry, low-pass limits, sampling policy, and cluster2D command lines
3. initialize references when needed
4. run staged `cluster2D` iterations
5. optionally run probabilistic pre-alignment in later stages
6. update particle class, in-plane, shift, sampled, and update-count state
7. restore class averages through shared-memory or distributed class-average pathways
8. run a final fill-in assignment pass for active particles that were never updated
9. generate final class averages, FRC metadata, and ranked outputs

The main policy boundary is:

- particle-domain work owns particle sampling, probabilistic assignment tables, search, class assignment, shift/in-plane updates, and partition-local outputs
- class-average assembly/restoration owns class-average sums, even/odd outputs, merged class averages, FRC/class documents, and project output metadata

## 3. Ownership Policy

`simple_commanders_abinitio2D.f90` owns:

- the `abinitio2D` entry point
- top-level defaults
- run orchestration across stages
- initial reference handling
- final fill-in dispatch
- final class-average generation/ranking

This layer should stay thin enough that stage rules are readable elsewhere.

`simple_abinitio2D_controller.f90` owns:

- stage counts and stage constants
- low-pass limit helpers
- stage-local `cluster2D` command construction
- search-mode policy by stage
- sampled-update policy, including `NPTCLS2SAMPLE_2D` and `nsample` override handling
- the rule that stage 1 may sample particles but does not fractionally restore previous class averages

`simple_cluster2D_strategy.f90` owns:

- shared-memory versus distributed execution selection
- iteration control inside one `cluster2D` invocation
- scheduler interaction
- probabilistic pre-alignment dispatch
- distributed worker scheduling
- distributed class-average assembly dispatch
- convergence and run-finalization bookkeeping

`simple_strategy2D_matcher.f90` owns:

- particle-domain alignment/search
- reproduction of the probabilistic sampled subset when `prob_align2D` is active
- strategy-object selection
- sigma updates during Euclidean search
- writing orientation updates
- writing distributed partial class-average sums when running as a worker

`simple_commanders_mkcavgs.f90` and the classaverager modules own:

- explicit class-average assembly from partial sums
- merged/even/odd class-average output
- class-document generation
- class FRC output and project output metadata

## 4. Sampling and Fractional-Update Policy

`abinitio2D` uses a fixed run-local target sample size:

- default: `NPTCLS2SAMPLE_2D = 200000`
- override: `nsample=<integer>`

The effective update fraction is:

```text
update_frac_2D = min(1.0, real(min(nptcls_eff, nsample_target_2D)) / real(nptcls_eff))
```

where `nptcls_eff` is the number of active particles with `state > 0`.

Stage policy:

- stage 1 uses a random sampled subset but disables fractional carry-over of previous class-average sums
- stages 2 and later use sampled update with fractional class-average restoration when the sample is smaller than the active set
- probabilistic stages preserve sample-once-and-reuse: `prob_align2D` chooses the subset, and `prob_tab2D`/`cluster2D_exec` reproduce that subset rather than resampling
- final fill-in is an assignment-only pass for active particles that still have `updatecnt == 0`

The desired restoration model is class-local: each class average should carry forward previous sums according to the realized sampled fraction for that class. The current implementation has moved toward this policy; changes in this area should preserve class-local semantics where available and avoid reintroducing a single ambiguous global owner for sampled-update state.

## 5. Iteration Semantics

For `cluster2D`:

- `startit` is the stage/invocation start
- `which_iter` is the current iteration
- `extr_iter` tracks the 2D extrapolation/search schedule
- `endit` is written after an invocation finishes and is consumed by the next stage setup

Do not collapse these counters into one another. Child command lines, including probabilistic pre-alignment and fill-in, must preserve the distinction between stage start and current iteration.

## 6. Artifact and Handoff Policy

Stable 2D workflow artifacts include:

- `assignment_part*.dat` and `assignment.dat`
- `dist_part*.dat` and `dist.dat`
- `cavgs_even_part*.mrc`, `cavgs_odd_part*.mrc`
- `ctfsqsums_even_part*.mrc`, `ctfsqsums_odd_part*.mrc`
- `cavgs_iterNNN.mrc`, `cavgs_iterNNN_even.mrc`, `cavgs_iterNNN_odd.mrc`
- `FRCS_FILE`
- `sigma2` iteration files
- `ptcl2D`, `cls2D`, `cls3D`, and `out` project segments

Partition-local probabilistic assignment/dist files are per-iteration artifacts and should be removed before the next distributed iteration writes new ones. Class-average partial sums are different when fractional restoration is active: they are the carry-over input for the next iteration and must be preserved until the worker has read and updated them.

## 7. Review Checklist

For any `abinitio2D` or `cluster2D` change, check:

- Does the command layer remain mostly orchestration?
- Is stage policy in the controller rather than scattered through matcher or strategy code?
- Does probabilistic 2D sample once and then reproduce the same subset?
- Do shared-memory and distributed paths use the same scientific workflow?
- Are class-average assembly/restoration responsibilities explicit?
- Are stale distributed handoffs removed without deleting fractional class-average carry-over inputs?
- Are `startit`, `which_iter`, `extr_iter`, and `endit` semantics preserved?
- Does fill-in remain assignment-only unless the policy is explicitly changed?
- Does the change preserve Cartesian-only `abinitio2D`?

## 8. Rules to Preserve During Refactors

- Do not reintroduce `polar=yes` into `abinitio2D`.
- Do not bury stage-policy tables in the matcher.
- Do not let probabilistic pre-alignment and matcher update sample different particle subsets.
- Do not make distributed-only class-average assembly semantics diverge from shared-memory scientific behavior.
- Do not treat final fill-in as a normal class-average restoration stage.
- Do not reuse stale assignment files as valid current-iteration inputs.
- Do not delete class-average partial sums at the start of a fractional-update iteration; workers need them as previous-sum carry-over.
