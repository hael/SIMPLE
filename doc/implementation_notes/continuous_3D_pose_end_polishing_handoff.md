# Continuous 3-D PCG pose-polishing scientific handoff

**Status:** truth-controlled scientific gate failed; research decision required
**Implementation default:** `pcg_pose_polish=no`
**SPEC:** [continuous_3D_pose_end_polishing_spec.md](continuous_3D_pose_end_polishing_spec.md)
**PLAN:** [continuous_3D_pose_end_polishing_plan.md](continuous_3D_pose_end_polishing_plan.md)

## What exists

The experimental component is owned by a top-level `reconstruct3D` command.
With `rec_backend=pcg`, `combine_eo=no`, and `pcg_pose_polish=yes`, it performs:

1. a fixed two-iteration even/odd PCG reconstruction;
2. one joint rotation-and-shift polish against fixed same-half maps;
3. persistence of accepted complete poses; and
4. one final fixed two-iteration PCG reconstruction.

The numerical derivative, bounded LM, half ownership, symmetry/gauge,
persistence, operator-contract, and regression tests have passed their recorded
Oracle gates. The corrected production objective applies the inverse KB
envelope exactly once and uses the PCG workspace shell range.

Major source anchors:

- top-level base/polish/final sequence:
  `src/main/commanders/simple/simple_commanders_rec.f90:36` and `:69`;
- activation and invalid-route contract:
  `src/main/params/simple_parameters_phases.f90:767`;
- distributed master ownership:
  `src/main/strategies/parallelization/simple_rec3D_strategy.f90:306`;
- fixed-map preparation and inverse envelope:
  `src/main/strategies/search/simple_pcg_pose_polisher.f90:343` and `:347`;
- particle-level OpenMP LM loop and deterministic summary reduction:
  `src/main/strategies/search/simple_pcg_pose_polisher.f90:158` and `:164`;
- truth runner and predeclared analyzer thresholds:
  `production/tests/continuous_3D_pcg_pose_validation/run_truth_diagnostic.sh:1`
  and `analyze_truth_matrix.py:26`.

The first 64-thread truth run exposed a serial pose-polishing loop after the
parallel PCG reconstruction. The corrected implementation uses the active
OpenMP thread count, parallelizes independent particle LM systems, and reports
batch progress. Oracle compilation, one-thread/four-thread equivalence, the
ten-case mother suite, a complete 64-worker truth matrix, and the focused
500-particle-batch runtime gate have now completed. Large-N scaling remains
unmeasured.

## Important workflow limitation

This is not integrated into the normal `abinitio3D` or `refine3D_auto`
workflow. The standard scripts under `SIMPLE_data_testing` call those workflow
programs but do not call `prg=reconstruct3D` directly. The public parameter is
accepted only by `reconstruct3D`, and distributed workers remove it so the
master owns one polish.

Consequently, a direct reconstruction A/B tests the isolated component. It
does not show that a normal SIMPLE refinement run uses polishing. Frozen
beta-gal and HPC direct-command A/B tests are deferred until the intended
workflow boundary is decided.

## Final Oracle experiment

Purpose: test whether the corrected production operator keeps exact truth poses
stationary and recovers controlled perturbations.

Oracle directory `continuous_3D_pose_truth_diagnostic_20260818_114404`
completed all 16 arms. Every enabled arm reported `REQUESTED 64 ACTIVE 64
BATCH 500` and completed in `68.5`--`70.2` seconds. The harness and 64-thread
execution passed, but the scientific analyzer failed. Clean exact rotation
drifted to `0.0038695` rad;
clean perturbed rotation improved by only about 17 percent; noisy exact FSC
declined; and noisy perturbed truth-average FSC declined by `0.00630986`.
Nominal full/FSC05/FSC0143 noisy arms were identical, so `lp` did not create
distinct pose-objective shell policies. Review the retained
`analysis/truth_matrix.md` and PLAN Phases 10F--10G before proposing another
numerical change. Its analysis JSON is byte-identical to the earlier
64-particle-batch run, so batching improved runtime without changing the
scientific outcome.

## Manual code-review conclusion

The core production path was reviewed manually after the batch-500 gate. The
review checked the command sequence, PCG preprocessing, fixed half-map and
inverse-envelope preparation, CTF/noise transfer, right-increment rotation
derivative, positive shift phase, five-parameter Gauss--Newton system, LM gain
ratio, rollback, persistence, and final reconstruction. No basic sign,
handedness, weighting, or pose-persistence defect was found.

This makes the failure less likely to be a misunderstanding about the coded
derivative. It does not prove the proposed method is scientifically wrong. The
tested same-half reconstructed reference is biased toward the particles it
polishes, and the independent simulator includes a finite-box operation that
the declared PCG operator intentionally excludes. The optimizer can therefore
lower its fixed-map objective while truth pose or FSC worsens.

Three narrower implementation gaps remain part of the handoff:

- the `lp` arguments did not set the pose workspace shell range, so the three
  nominal frequency arms are not three frequency experiments;
- stencil-switch margins are measured but not used as a step safeguard; and
- the production truth matrix covers `ml_reg=no`, `pgrp=c1`, and no held-out
  acceptance signal. The common-FSC regularized reference and non-`c1` routes
  require separate evidence if a revised SPEC selects them.

These gaps do not justify an agent-selected patch. They depend on the reference,
frequency, regularization, symmetry, and acceptance policies requested below.

## `reconstruct3D` caller map

`commander_rec3D` is called by the public `reconstruct3D` dispatcher, by
`refine3D_auto`, `refine3D_multi`, the base refinement strategy,
`abinitio3D`, `abinitio3D_cavgs`, bootstrap reconstruction, random
reconstruction, and nano-trajectory reconstruction. The source anchors are:

- public dispatch:
  `production/simple_exec.f90:53`,
  `src/main/exec/simple_exec_refine3D.f90:39`, and
  `src/main/commanders/simple/simple_commanders_rec.f90:36`;
- distributed worker dispatch:
  `production/simple_private_exec_driver.f90:322`;
- `refine3D_auto` and `refine3D_multi`:
  `src/main/commanders/simple/simple_commanders_refine3D.f90:200`, `:216`,
  `:670`, and `:990`;
- missing starting-volume reconstruction:
  `src/main/strategies/parallelization/simple_refine3D_strategy.f90:1033`;
- ab initio reconstruction, symmetry, split, and final helpers:
  `src/main/simple_abinitio_utils.f90:507`, `:584`, and `:744`, called from
  `src/main/commanders/simple/simple_commanders_abinitio.f90:244`, `:283`,
  `:821`, `:893`, `:896`, `:912`, and `:924`;
- bootstrap and random reconstruction:
  `src/main/commanders/simple/simple_commanders_rec.f90:116`, `:122`, and
  `:431`; and
- nano-trajectory reconstruction:
  `src/main/commanders/single/single_commanders_nano3D.f90:332` and `:374`.

No higher-level caller currently propagates `pcg_pose_polish=yes`; all retain
the default `no`. Broad propagation would affect startup, split, symmetry,
bootstrap, random, trajectory, and final reconstructions. It would not mean
one well-defined end-polishing pass. Hans must choose the owning caller and
stopping protocol first.

Run from the Oracle validation checkout:

```bash
bash production/tests/continuous_3D_pcg_pose_validation/run_truth_diagnostic.sh \
  --truth-volume "$HOME/Projects/frozen_noisy_betagal_20260813_153030/frozen_fixture/truth_1JYX_box144.mrc" \
  --truth-oris "$HOME/Projects/frozen_noisy_betagal_20260813_153030/frozen_fixture/betagal_truth_oris.txt" \
  --noisy-stack "$HOME/Projects/frozen_noisy_betagal_20260813_153030/frozen_fixture/betagal_noisy.mrcs" \
  --box 144 \
  --smpd 1.3 \
  --mskdiam 120 \
  --pgrp c1 \
  --kv 300 \
  --cs 2.7 \
  --fraca 0.1 \
  --lp-fsc05 4.255 \
  --lp-fsc0143 3.671 \
  --nthr 64 \
  --output-parent "$HOME/Projects"
```

The runner continues through all 16 arms and creates one timestamped directory.
Download the complete directory. Read `analysis/truth_matrix.md` first, then
check `STATUS.txt`, `failures.txt`, `analysis/truth_matrix.json`,
`run_config.env`, the source and fixture manifests, and the individual arm
logs. Do not accept or reject the method from process exit alone.

## Predeclared decisions

- Clean exact poses must remain within `0.001745` rad rotation RMS and `0.05`
  pixel shift RMS of the disabled control.
- Clean perturbed rotation and shift RMS errors must each decrease by at least
  50 percent.
- At least one noisy frequency band must keep exact poses stable and reduce
  both perturbed-pose RMS errors by at least 10 percent.
- No evaluated half-map or truth-map FSC area can decline by more than `0.005`.
- Invalid numerical terminal outcomes must be zero.

Interpret the result as follows:

- A harness failure is not scientific evidence. Correct only the harness and
  repeat the same frozen matrix.
- A clean-exact failure means the production fixed-map objective is not
  stationary at the known pose. Do not integrate the feature into refinement.
- Clean stability with failed clean recovery points to local pose numerics,
  scaling, or solver policy.
- Clean cases passing while every noisy band fails points to overfitting or a
  missing frequency/regularization policy.
- Passing the complete gate qualifies the isolated component for workflow
  design work; it does not by itself authorize production refinement use.

## Decisions requested from Hans and the refinement researchers

1. Is the intended product a single terminal polish or a repeated polish after
   selected reconstructions inside `abinitio3D` or `refine3D_auto`?
2. If repeated, which reconstruction owns it: every iteration, only combined
   iterations, only the final iteration, or an explicit new stage?
3. Which half-map and frequency policy must remain fixed during each particle
   update, and when may a reconstructed map be refreshed?
4. What convergence, rollback, and restart contract applies when accepted poses
   feed the next refinement iteration?
5. What dataset-level acceptance metric is required in addition to the
   truth-controlled pose and FSC gates?

No workflow integration should start until these answers are written into a
revised FINAL SPEC. Keep the feature disabled by default.
