# Specification: table-free class stability for `abinitio2D_sgd`

## Status

IN REVIEW.

This specification is the complete start-to-end contract for table-free class
stability in the opt-in `abinitio2D_sgd` workflow. It governs
[`abinitio2d_sgd_class_stability_plan.md`](abinitio2d_sgd_class_stability_plan.md).
After requester approval, its status becomes `FINAL` and its requirements are
frozen. Implementation details, phase checkpoints, and test evidence belong in
the PLAN. A material change to this contract requires a major requirement
review before implementation continues.

## Request and context

The opt-in `abinitio2D_sgd` stream assigns one hard class/rotation winner per
sampled particle and then refines that winner's shift. Exhaustive competition
can create a winner-takes-all feedback loop: a temporarily weak class receives
fewer particles, its average becomes noisier, and it loses still more support.

SIMPLE defines a 2D class with fewer than two particles as nonviable for the
next search. The streamed class-average safeguard preserves only a class with
zero assigned particles. The search and restoration viability boundaries
therefore do not match.

## Problem statement

The current hard exhaustive stream can destroy a previously meaningful class
during one stochastic update. Every sampled particle compares against every
currently viable class, and a class that falls to population zero or one is
not searchable on the next iteration. Zero-only image preservation does not
protect a singleton, and preservation alone does not retain the minimum
particle support required for future competition.

The class-stability mechanism must reduce this collapse without restoring
probability tables, assigning unrelated particles to weak classes, or changing
standard SIMPLE behavior.

## Initial state / entry point

Before this class-stability work:

1. `simple_exec prg=abinitio2D` is the conventional production workflow and
   contains no SGD-specific class-stability behavior.
2. `simple_exec prg=abinitio2D_sgd` is an opt-in development workflow.
3. Every active streamed iteration draws a fresh random particle subset using
   the existing SGD fraction, 60 percent by default, capped by SIMPLE's
   existing `NSAMPLE_DEFAULT_2D=200000` limit.
4. Each sampled particle exhaustively evaluates every viable class and the
   complete discrete rotation grid with the raw Euclidean objective.
5. One hard class/rotation winner is retained and the existing direct-gradient
   shift minimizer refines only that winner.
6. Classes with population below two are nonviable for the next search.
7. Streamed restoration preserves a previous class average only at population
   zero, not population one.
8. No bounded class frontier or two-member support floor influences production
   assignments.

## Desired outcome

Provide an explicit, table-free class-stability policy for `abinitio2D_sgd`
that limits destructive global class competition while retaining one hard
winner per sampled particle. Preserve previously viable class models and at
most the minimum two incumbent particles needed to keep those classes viable.

The mechanism must be deterministic for fixed input state and seed, bounded in
memory, safe under OpenMP scheduling, and independently selectable from the
existing exhaustive stream.

## Final endpoint / required outcomes

The completed work has all of the following properties:

1. `sgd_class_policy=exhaustive|frontier` is the only public class-stability
   selector.
2. `exhaustive` is the default and preserves the pre-change streamed behavior.
3. `frontier` is accepted only by `abinitio2D_sgd` and affects only iterations
   where the SGD stream is already active.
4. The frontier is deterministic, bounded, duplicate-free, retains a valid
   previously viable incumbent, draws alternatives only from viable classes,
   and periodically performs restart-stable full exploration.
5. Candidate restriction changes only class competition. The raw Euclidean
   class/rotation evaluator, hard-winner semantics, and direct shift minimizer
   remain the same in exhaustive and frontier modes.
6. After provisional winners are known, every previously viable class with a
   valid previous model retains up to the two lowest-penalty particles from
   its own incumbents when needed to reach population two.
7. The support floor never recruits a particle from an unrelated incumbent
   class and never fabricates support when fewer than two eligible incumbents
   remain.
8. A previous valid class model is preserved when final support is zero or
   one, while the stored population remains truthful.
9. Non-sampled particles are not silently reassigned, and every sampled
   particle still has exactly one final hard assignment.
10. Detailed class-stability reporting remains behind the existing
    default-off `sgd_diagnostic` switch.
11. Focused, production-context, and matched scientific tests demonstrate
    better class survival without unacceptable class-quality or runtime
    regression.

## Why this matters

Preserving meaningful class diversity is more important than improving only
the single best reported class. A fast hard-assignment stream is not useful if
its stochastic updates eliminate plausible class averages and make the final
class distribution scientifically misleading.

The requested policy aims to retain streaming memory and runtime advantages
while making class loss conservative, observable, and testable.

## Scope

- Freeze and measure the existing exhaustive stream.
- Reuse the existing 60-percent stochastic sample and 200,000-particle cap.
- Centralize and test the existing two-particle viability boundary.
- Add default-off sampling, support, migration, margin, and timing diagnostics.
- Provide a pure deterministic class-frontier constructor.
- Add the development-only `sgd_class_policy` selector with default
  `exhaustive`.
- Integrate the frontier only into active `abinitio2D_sgd` streamed searches.
- Add a two-member floor using only a class's own eligible incumbents.
- Preserve a previous valid model at final populations zero and one.
- Add focused policy tests and a deterministic multi-iteration production-
  context fixture.
- Run matched exhaustive/frontier comparisons on beta-galactosidase and
  apoferritin from identical preprocessing checkpoints and fixed seeds.

## Out of scope

- Any change to standard `abinitio2D` behavior or UI.
- Enabling the frontier by default.
- SoftMax, top-K tables, likelihood tables, `eulprob_tab2D`, or probabilistic
  global assignment.
- Recruiting unrelated particles to rescue a weak class.
- Forced equal class populations or unconditional empty-class bonuses.
- A support-aware loss penalty in this version.
- Class-balanced input sampling presented as the collapse solution.
- Changes to the raw Euclidean objective, discrete rotation search, direct
  shift optimizer, stage schedule, or terminal conventional pass.
- Changes to continuous in-plane refinement, L-BFGS-B, or a replacement pose
  optimizer.
- Making candidate-count or exploration-cadence controls public in the first
  version.

## Opt-in interface contract

The interface is:

```text
sgd_class_policy=exhaustive|frontier
```

- Default: `exhaustive`.
- Owner: the development command `abinitio2D_sgd`.
- `abinitio2D` must not display, accept, set, or forward the key.
- `frontier` is inactive wherever the existing SGD stream is inactive.
- Unsupported values fail typed command validation with a meaningful error.
- Candidate-count and full-exploration cadence remain owned implementation
  constants or constructor inputs for the first version. They are frozen in
  the PLAN before matched scientific runs and are not user-facing controls.

## Acceptance criteria

### Compatibility and route identity

1. Omitting `sgd_class_policy` is equivalent to
   `sgd_class_policy=exhaustive`.
2. Exhaustive mode reproduces the existing viable-class candidate order,
   hard winners, stored score semantics, and direct-shift calls for identical
   input state.
3. Standard `abinitio2D`, SGD-inactive iterations, and the terminal
   conventional pass do not construct or call the frontier helper.
4. No candidate table or particle-by-class allocation is introduced.

### Frontier behavior

5. A frontier is bounded and duplicate-free.
6. A valid, previously viable incumbent is included; an invalid or never-
   viable incumbent is not specially retained.
7. All non-incumbent candidates are currently viable classes.
8. Fixed project state, seed, iteration, and particle index produce identical
   candidates independent of thread scheduling and call order.
9. Periodic full exploration is deterministic and restart-stable.
10. Both policies use the same raw Euclidean evaluator and direct shift
    minimizer and retain one final hard winner per sampled particle.

### Two-member support and model preservation

11. Populations zero and one are under-supported; population two is viable.
12. Only a previously viable class with a valid previous model is eligible for
    protected anchors.
13. Only that class's own incumbent particles are eligible, and the smallest
    finite incumbent-loss penalties are retained.
14. No more than two anchors per class are needed, final assignments remain
    unique, and final populations are recomputed exactly.
15. An impossible floor reports class death and does not recruit an unrelated
    particle.
16. A valid previous class image remains finite and unchanged when final
    population is zero or one; the recorded population is not falsified.

### Regression and scientific outcome

17. Focused sampling, frontier, support-floor, restoration, and route-identity
    tests pass, followed by the complete SGD regression suite.
18. A deterministic multi-iteration fixture proves that non-sampled particles
    are unchanged, sampled particles have one final winner, protected classes
    receive exactly their eligible lowest-penalty anchors, and accepted shifts
    remain finite, improving, and bounded.
19. Matched experiments use beta-galactosidase and apoferritin, at least three
    predeclared fixed seeds, identical preprocessing checkpoints, and identical
    stage schedules for exhaustive and frontier modes.
20. For every matched arm, frontier mode must not finish with fewer viable
    classes than exhaustive mode. The median final viable-class count across
    seeds must be strictly higher for both datasets.
21. For each dataset, frontier mode's across-seed median ranked-class
    resolution and median resolution among classes below 9 Angstrom must each
    be no more than 0.5 Angstrom worse than exhaustive mode.
22. For each dataset, frontier mode's median runtime must not exceed exhaustive
    mode by more than 10 percent.
23. Every matched arm completes normally and produces a consistent readable
    project. One improved best class alone does not satisfy this contract.

## Constraints and invariants

- Existing unrelated worktree changes must remain untouched.
- `abinitio2D_sgd` remains development-labelled and opt-in.
- The existing SGD stage policy owns when streaming is active.
- The existing fresh mini-batch uses the current SGD fraction and
  `NSAMPLE_DEFAULT_2D`; do not use or modify legacy `update_frac` for this work.
- No mutable random-number generator is shared across OpenMP threads.
- Frontier state is linear in sampled particles plus classes, not particles
  times classes.
- Detailed diagnostics are default-off and must not determine numerical
  behavior.
- Linux/BOX or HPC results are claimed only from observed logs.

## Assumptions

- The production matcher can supply the incumbent class, whether that class
  was previously viable, and whether its previous model is valid without
  inspecting image data inside the frontier helper.
- Raw incumbent and provisional-winner losses can be retained for sampled
  particles without constructing a class table.
- Three fixed seeds per dataset provide the minimum matched evidence for this
  development-stage decision; broader scientific validation may still be
  required before any future default change.

## Open questions

No implementation-blocking question remains in this proposed contract.

Requester approval is required for the opt-in interface and the numerical
matched-experiment thresholds before the status can become `FINAL` and Phase 2
production integration can begin.
