# Hard-stream-seeded L-BFGS-B experiment for `abinitio2D`

## Status

Closed on 2026-08-11. The experiment was implemented through workflow testing,
reviewed, rejected, and removed before any implementation commit or push.

The cleanup was performed on `master` at `c8045dd41`, with `origin/master` at
`2f151c1ad`. All production, parameter, UI, numerical, and plan-specific test
changes were restored exactly to their `HEAD` versions. The only tracked changes
retained from the work are this completed record and two independent SGD fixture
argument fixes described under **Cleanup implementation**.

This document is a postmortem and prevention reference, not an authorization to
reintroduce `hard_stream`.

## Final decision

Do not ship or continue the `hard_stream=yes` feature.

The proposed optimization was intended to eliminate a wrapper-level all-angle
rescan before joint L-BFGS-B. Repository review after Phase 4 established that
the current production `inpl_cont=yes` path had already eliminated that rescan:

```text
production conventional path
    strategy2D_srch%refine_selected_continuously
    -> minimize_joint(..., irot_in=self%best_rot)
    -> local seeded L-BFGS-B

experimental hard_stream path
    strategy2D_srch%refine_selected_continuously
    -> minimize_joint_seeded_lbfgsb(self%best_rot, ...)
    -> minimize_joint(..., irot_in=irot_seed)
    -> the same local seeded L-BFGS-B
```

The experimental wrapper changed the result representation at its boundary, but
did not remove work from the production route. The Phase 2 evaluation-count
reduction compared the wrapper against a test-only all-angle invocation, not
against the actual conventional production call.

The HPC workflows therefore supply non-regression evidence only. They cannot
support a runtime or resolution-improvement claim.

## Scope of the concluded experiment

The implementation was deliberately narrowed to ordinary non-SGD
`abinitio2D`:

- conventional arm: `abinitio2D objfun=euclid inpl_cont=yes`;
- opt-in arm: `abinitio2D objfun=euclid inpl_cont=yes hard_stream=yes`;
- `hard_stream=no` or `inpl_cont=no` was intended to preserve the original path;
- `abinitio2D_sgd` and its `inpl_cont=yes` rejection were out of scope;
- no SoftMax, top-K, probability-table, stage-policy, terminal-pass, or SGD
  behavior change was authorized;
- no L-BFGS-B internals were to be changed.

The intended numerical transaction was:

```text
selected hard class, integer rotation, and native shift
    -> local L-BFGS-B over (sx, sy, rho)
    -> raw loss and validity result
    -> coupled shift, nearest inpl, fractional e3, and exp(-L) score
```

That transaction remains useful as a description of the existing local joint
solver. It did not justify a second public route.

## Original hypothesis

The plan began with this performance hypothesis:

1. the conventional continuous wrapper rescans all `nrots` angles at the
   selected shift;
2. the preceding discrete search already owns the correct integer winner;
3. supplying that winner directly would remove `nrots` objective evaluations;
4. retaining the same L-BFGS-B solver would preserve quality while reducing
   runtime.

Items 2 and 4 were reasonable. Item 1 was false for the production route under
development. The plan noticed that `minimize_joint` supported an authoritative
`irot_in` mode, but failed to complete the final caller trace and confirm that
`strategy2D_srch` already used that mode. That missed call-site fact invalidated
the performance premise before implementation started.

## Work performed

### Phase 0: freeze current behavior

The existing route and supporting tests were recorded before production
integration.

Observed Oracle/Linux evidence:

- UI visibility passed 151 assertions;
- `abinitio2D_sgd` dispatch passed all 15 assertions;
- the complete six-case SGD base suite passed;
- the corrected-gradient test passed with maximum error `5.29e-5` against a
  `1e-2` tolerance;
- integer-angle parity passed with worst error `1.38e-7` against a `1e-2`
  tolerance;
- the route fixture was corrected so its hybrid-objective case used the valid
  default-off continuous state.

Two SGD fixture helpers were also corrected so a mother-suite `case=` selector
is not parsed as an old-style SIMPLE program argument. This fix is independent
of `hard_stream` and is retained after cleanup.

### Phase 1: isolated seeded API

An experimental `minimize_joint_seeded_lbfgsb` entry point was added to
`simple_pftc_shsrch_grad`.

Its intended contract was:

- immutable integer seed plus native-frame shift;
- the existing joint raw-Euclidean callback and L-BFGS-B object;
- no all-angle evaluator inside the seeded entry point;
- raw merit `-L` at the new boundary;
- separate output rotation and fractional rotation;
- explicit valid, improved, and no-improvement outcomes;
- invalid evaluation returned no committable pose;
- valid no-improvement retained the supplied seed exactly.

Observed Oracle/Linux evidence:

- invalid integer seed behavior passed;
- valid no-improvement retained the seed;
- improved cases returned comparable final raw losses;
- periodic windows crossing the first and last rotation cells passed;
- seed 288 could return unwrapped fractional index 289 and canonical bin 1.

The implementation delegated to the established local
`minimize_joint(..., irot_in=irot_seed)` path, which preserved the mature solver
but also meant it was not a new production algorithm.

### Phase 2: focused A/B mechanics

The route fixture was expanded to three references, three native shifts per
reference, and periodic boundaries at both grid ends.

Observed Oracle/Linux results:

- 9 resolved A/B cases passed;
- 6 periodic-edge cases passed;
- resolved hard and rescan winners selected equivalent raw losses;
- final raw losses, rotations, outcomes, and gradient counts agreed;
- the artificial Case A/B objective counts differed by exactly 288;
- objective median/tail were `304/321` for Case A and `16/33` for Case B;
- gradient median/tail were `16/33` for both;
- CPU median/tail were approximately `0.011107/0.029149` seconds for Case A
  and `0.010463/0.028500` seconds for Case B.

This proved that omitting an all-angle scan saves its objective evaluations when
such a scan exists. It did not prove a production saving, because the production
conventional call already supplied `irot_in` and skipped the scan.

The test accidentally made the hypothetical baseline executable by changing
`minimize_joint` so `irot_in` became optional and then calling it without the
seed. That was a test-created baseline and an inappropriate production-API
change. Both were removed during cleanup.

### Phase 3: production opt-in

The experimental production integration added:

- a binary `hard_stream` parameter with default `no`;
- parser registration and binary validation;
- developer UI metadata;
- propagation to ordinary `abinitio2D` child stages;
- a separate `hard_stream_active` strategy flag;
- the two-flag production branch;
- raw-loss to conventional-score conversion;
- UI and production-route diagnostics.

Observed Oracle/Linux evidence from `phase3_20260810_143632`:

- UI visibility passed 158 assertions;
- SGD dispatch remained 15/15;
- Stage 1 validation ended normally;
- route identity ended normally;
- the explicit route diagnostic reported 33 objective and 33 gradient
  evaluations and a positive conventional score.

Later review found three integration defects:

1. The seeded branch converted `objective_final` with
   `exp(-objective_final)` instead of the established
   `exp(-max(0, objective_final))`. Because the joint evaluator permits small
   negative series undershoots, the experimental branch could persist a score
   greater than one.
2. `hard_stream` was added to a descriptor shared by `abinitio2D` and
   `abinitio2D_sgd`, so the supposedly non-SGD parameter leaked onto the SGD UI
   surface even though the SGD numerical route remained disabled.
3. The new `hard_stream_active` lifecycle state was not reset in `kill`.

These defects made the implementation unacceptable independently of the failed
performance premise. They were removed by restoring the affected files rather
than extending an abandoned feature.

### Phase 4: workflow evidence

The established apoferritin and beta-galactosidase HPC preparation scripts were
adapted to run the conventional and explicit opt-in commands from isolated
copies of one extracted project.

Observed HPC evidence from run stamp `20260810_150234`:

| Dataset | Conventional | Hard-seeded | Best resolution | Result |
| --- | ---: | ---: | ---: | --- |
| Apoferritin | 128.1 s | 129.9 s | 4.03 A / 4.03 A | normal stop in both arms |
| Beta-galactosidase | 133.9 s | 132.8 s | 7.23 A / 7.23 A | normal stop in both arms |

Apoferritin changed by about `+1.4%`; beta-galactosidase changed by about
`-0.8%`. Those differences are consistent with noise and show no useful runtime
effect. Best resolution was unchanged.

The two arms copied the same extracted input, but separate SIMPLE processes
seeded their stochastic trajectories independently. The runs were therefore
input-matched, not fully random-trajectory-matched. This is adequate to record
that both modes completed with comparable output, but not to attribute small
quality or timing differences to the option.

The HPC scripts were also improved so future runs use one
`Exp-<timestamp>` directory containing preparation outputs, both comparison
projects, run logs, summary CSV, and SLURM `.out`/`.err` files. That operational
improvement is independent of the rejected feature and is retained.

## What the final review established

### 1. The production premise was invalid

The first mandatory review question should have been:

```text
What exact production call is replaced, and what exact work does it perform?
```

The answer would have shown that the conventional production caller already
passed its authoritative integer seed. The experiment should have stopped
before Phase 1.

### 2. An API capability is not evidence of production use

The presence of an all-angle initializer or selector in a numerical owner does
not prove that the active workflow invokes it. Performance plans must trace the
complete active call path from command to numerical call, including actual
optional arguments and constructor mode.

### 3. Tests must not manufacture the claimed baseline

The Phase 2 comparison modified a production API to make `irot_in` optional,
then omitted it to create Case A. Exact counter differences against that case
were real but irrelevant to current production. A test-only route cannot be the
baseline for a production speed claim.

### 4. Boundary representations require independent review

Even a wrapper around the same solver can introduce scientific defects. Raw
loss `L`, raw merit `-L`, and stored score `exp(-L)` are different contracts.
The established clamp protecting the `[0,1]` score range was lost in the
experimental strategy conversion.

### 5. A shared UI descriptor can violate scope without changing numerics

Adding a parameter to shared metadata exposed it to a command that was
explicitly out of scope. UI, parser, controller, strategy, and numerical
ownership must all be traced separately.

### 6. A matched checkpoint is not a matched stochastic trajectory

Copies of one input project control data and schedule, but independent random
initialization still prevents causal interpretation of small A/B differences.
Future stochastic comparisons must preserve random state or use enough paired
replicates to estimate run-to-run variance.

## Cleanup implementation

The rejected experiment was removed by comparing every changed tracked file to
`HEAD` and restoring the plan-scoped files exactly. No replacement algorithm or
additional production code was written.

Restored exactly to `HEAD`:

- `production/tests/simple_test_continuous_inplane_rotation2D_route_identity.f90`;
- `production/tests/simple_test_ui_visibility.f90`;
- `src/main/params/simple_parameters.f90`;
- `src/main/params/simple_parameters_parse.f90`;
- `src/main/params/simple_parameters_phases.f90`;
- `src/main/pftc/simple_pftc_shsrch_grad.f90`;
- `src/main/simple_abinitio2D_controller.f90`;
- `src/main/strategies/search/simple_strategy2D_srch.f90`;
- `src/main/ui/simple/simple_ui_cluster2D.f90`.

This restoration removes:

- the `hard_stream` parameter, parser registration, validation, UI metadata,
  controller propagation, strategy state, and production branch;
- `minimize_joint_seeded_lbfgsb`;
- the experimental optional-`irot_in` all-angle route;
- the unused `minimize_joint_rounded` addition from this worktree;
- Phase 1/2 seeded API and artificial-rescan tests;
- Phase 3 route and UI tests.

Retained unchanged during cleanup:

- `production/tests/simple_sgd_base_suite_baseline_test.f90`;
- `production/tests/simple_sgd_base_suite_v2_test.f90`;
- `.codex/hpc_apof_prepare.sh`;
- `.codex/hpc_betagal_prepare.sh`;
- `.codex/hpc_summarize_sgd_runs.py`;
- downloaded Oracle/HPC logs and workflow projects outside the repository.

The two tracked fixture files retain only the independent `case=` mother-suite
argument handling. The three `.codex` files retain the consolidated experiment
directory and Python-compatibility improvements. None depends on `hard_stream`
for ordinary repository operation; their historical command labels may be
repurposed in a future, separately reviewed A/B experiment.

No implementation commit was created and nothing was pushed. Any later commit
containing the two fixture fixes or this postmortem must be explicitly requested
and must not be described as a seeded-L-BFGS-B performance implementation.

## Verification after cleanup

The cleanup acceptance checks are intentionally lightweight because the
production and focused-test sources are restored byte-for-byte to `HEAD`:

1. `git diff HEAD` shows no change in any restored file.
2. No tracked production or test source contains `hard_stream` or
   `minimize_joint_seeded_lbfgsb`.
3. No unmerged paths or conflict markers remain.
4. `git diff --check` passes.
5. The remaining tracked diff contains only this postmortem and the two retained
   SGD fixture fixes.

Compilation and runtime reruns are unnecessary for the exact `HEAD`
restorations. The retained fixture fixes already passed in the observed full
Oracle/Linux SGD base suite. No new executable behavior was added by cleanup.

## Mandatory prevention rules for future performance plans

Before implementing a performance feature:

1. **Name the exact production caller and callee.** Record the full call with
   actual keyword arguments, constructor mode, stage, and command route.
2. **Measure the production baseline before creating an alternate API.** A
   counter or profile must prove the suspected work occurs in production.
3. **Compare against the production call, not a public API possibility.** Do
   not omit an argument or enable an unused mode merely to create Case A.
4. **Make the predicted delta explicit.** State which function calls and how
   many objective, gradient, I/O, allocation, or synchronization operations
   will disappear.
5. **Stop after Phase 0 if the predicted work is absent.** Do not continue
   because a proposed API is technically implementable.
6. **Keep scientific representations explicit.** Document lower/higher-is-
   better direction and every conversion among raw loss, merit, correlation,
   probability weight, and persisted score.
7. **Trace UI ownership independently.** Shared descriptors and child command
   propagation must not expose an experiment to excluded commands.
8. **Use production-context tests before workflow timing.** The test must use
   the same constructor and arguments as production and prove the changed call
   count directly.
9. **Use matched stochastic evidence.** Preserve checkpoint, particles,
   references, schedule, and random state, or run enough paired replicates to
   characterize variance.
10. **Define the stop decision before implementation.** No measured production
    saving means removal, not reinterpretation as a feature.

## Conditions for reopening related work

Do not reopen `hard_stream` under the same hypothesis.

A future optimizer or direct-solver plan may proceed only after a fresh
production profile identifies a real cost and the plan names a distinct call
sequence that removes or reduces it. The new baseline must be the current local
`minimize_joint(..., irot_in=authoritative_seed)` route. Replacing L-BFGS-B is
justified only if a candidate solver achieves comparable raw losses and improves
measured production runtime, evaluation count, robustness, or maintainability.

The existing continuous in-plane refinement and its authoritative-seed behavior
remain the baseline; this concluded experiment does not invalidate them.
