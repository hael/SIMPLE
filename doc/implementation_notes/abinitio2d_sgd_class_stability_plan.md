# Table-free class stability for `abinitio2D_sgd`

## Status

Active implementation plan. This document supersedes the former direct-joint
optimizer plan as the next development priority for `abinitio2D_sgd`.

Repository baseline reviewed on 2026-08-12: local `master` at `dc71602e8`.
The branch is one commit ahead of `origin/master` because of unrelated
continuous-3D work. The worktree also contains unrelated continuous-3D files
and one continuous-in-plane route-identity test edit; none is part of this
plan.

The earlier hard-stream-seeded L-BFGS-B experiment is closed. Production
already supplied the selected rotation seed to the established joint polisher,
so the proposed extra seeding layer did not remove work from the real caller.
That experiment is not a prerequisite for this plan.

The completed development-commander plan is archived under
`doc/implementation_notes/completed/`. The standard `abinitio2D` command must
remain free of SGD-specific UI, checkpoint, scheduling, and search behavior.

## Decision summary

The next SGD problem is class stability, not replacement of the mature inner
pose optimizer.

1. Keep the current table-free hard stream and its validated direct shift
   refinement as the starting point.
2. Keep the existing fresh random particle mini-batch: 60 percent by default,
   capped by SIMPLE's existing `NSAMPLE_DEFAULT_2D=200000` limit. Do not add a
   second sample cap and do not reuse or modify legacy `update_frac`.
3. Replace the zero-only safeguard with SIMPLE's actual minimum viable class
   support: two particles. `prep_strategy2D_glob` converts populations below
   two to zero before search, so a one-particle class is already dead to the
   next matcher iteration.
4. Protect a previously viable class with its two best incumbent particles.
   This is a minimum support floor, not an empty-class rescue: it applies only
   to a class that previously had at least two members and a valid model.
5. Add a table-free class candidate frontier so one noisy iteration cannot
   make every particle abandon a previously meaningful class merely because
   all classes compete globally on every update.
6. Keep one final hard winner per sampled particle. Do not restore SoftMax, top-K
   likelihood tables, `eulprob_tab2D`, or the obsolete table SGD path.
7. Keep continuous in-plane L-BFGS-B as the reference local polisher. The
   development command currently rejects `inpl_cont=yes`; integrating that
   independent feature is not part of the class-stability experiment.
8. Do not introduce Adam, AdamW, or another stateful optimizer for the
   three-variable per-particle pose solve. That small bounded problem has exact
   gradients, and persistent moments or weight decay would add state and bias
   without addressing class collapse.
9. Retain the proposed direct joint projected solver only as deferred research.
   It should be reconsidered only if a focused comparison shows it matches
   L-BFGS-B final loss and pose quality with meaningfully fewer evaluations or
   lower workflow runtime.

## Current verified behavior

### Activation and schedule

`abinitio2D_sgd` is an opt-in development commander. Its existing
`sgd_stage4_mode` policy decides which iterations use the stream. This plan
does not change stage boundaries, warmup iterations, restart behavior, or the
terminal conventional pass. A class-stability policy is evaluated only when
`params%l_sgd_streaming_active` is already true.

Stages 1 and 2 remain conventional. Existing stage-policy tests remain the
authority for the exact `on` and `alternate` schedules.

### Particle mini-batch

On every active SGD iteration,
`sample_ptcls4update2D` counts active particles and computes

\[
f_{\mathrm{batch}} = \min\left(
    \mathtt{sgd\_update\_frac},
    \frac{200000}{N_{\mathrm{active}}}
\right).
\]

It then calls `sample4update_rnd`, so the default batch is

\[
N_{\mathrm{batch}} = \min(0.6N_{\mathrm{active}}, 200000).
\]

The sampled set is freshly drawn on each active iteration. This is the
stochastic outer loop. The SGD path must continue to use this SIMPLE-owned
sampling API rather than add another sampler or cap.

### Assignment and pose refinement

For every sampled particle, the current stream:

1. enters `strategy2D_greedy`;
2. reads class populations through `prep_strategy2D_glob`, which rewrites
   populations below two to zero, and skips those classes;
3. evaluates all rotation bins of every remaining class with
   `gen_raw_euclid_vals`;
4. selects the single minimum raw Euclidean class/rotation state; and
5. refines shift for that fixed winner with the bounded analytical-gradient
   direct minimizer.

In mathematical form,

\[
(c_i,r_i)=\arg\min_{c\in C_{\mathrm{viable}},r}L_i(c,r,s_i),
\]

followed by bounded refinement of `s_i` for the selected `(c_i,r_i)`.

There is no candidate table, SoftMax, top-K posterior, or global
class-frontier assignment in this path.

### Existing class-preservation gap

Before streamed restoration clears the accumulators, the class averager backs
up the previous class images. If a class receives exactly zero particles in
the current update and a previous model exists, `cavger_restore_cavgs` copies
the previous image instead of writing a zero image.

That boundary does not match search eligibility. `prep_strategy2D_glob`
contains the explicit rule

```fortran
where( s2D%cls_pops < 2 )
    ! ignoring classes with one particle
    s2D%cls_pops = 0
end where
```

Therefore a one-particle class is reconstructed as a singleton but is treated
as dead by the next search. The protection threshold must be `pop < 2`, not
`pop == 0`.

This safeguard is necessary but incomplete:

- it preserves the model only at population zero;
- it does not preserve a population-one model at the actual death boundary;
- it does not preserve support;
- current greedy search skips both zero- and one-particle classes after the
  population normalization; and
- a preserved class can therefore remain unavailable to particles after its
  support falls below two.

The existing focused test proves only the old `pop == 0` predicate and the
60-percent sample count. It does not cover population one or prove that an
end-to-end class retains the two members required by the production search.

## Problem statement

The hard exhaustive stream has a winner-takes-all feedback loop:

```text
slightly lower support
        -> noisier current average
        -> higher raw losses for more particles
        -> still lower support
        -> population below two
        -> normalized to dead for search
```

Preserving only an empty average misses the actual death boundary and does not
break the feedback loop. Independent and matched preprocessing experiments
showed that the stream can improve best-class resolution and runtime while
retaining fewer useful/ranked classes than the conventional path. A better
best class is not sufficient if biologically meaningful classes disappear.

The objective of this plan is therefore:

> Improve survival and recoverability of previously viable classes in the
> stochastic hard-assignment stream without restoring the global probability
> table or forcing particles into scientifically implausible classes.

## Terminology

- **Minimum viable support**: two particles, matching the existing
  `s2D%cls_pops < 2` search rule.
- **Viable class**: current population is at least two.
- **Under-supported class**: current population is zero or one.
- **Previously viable class**: the preceding assignment had at least two
  members and a valid previous class model exists.
- **Incumbent class**: the particle's class before the current SGD update.
- **Candidate frontier**: the bounded set of classes evaluated for one
  particle on a non-full-scan iteration.
- **Sticky candidate**: the incumbent class is guaranteed to be in that
  particle's frontier. It remains a candidate; it is not automatically the
  winner.
- **Sticky support anchors**: the two best incumbent particles retained for a
  previously viable class when its provisional new support would fall below
  two. Anchors are part of v1.
- **Full exploration**: an iteration that evaluates all currently viable
  classes. It does not offer every unrelated particle every under-supported
  class.

## Non-negotiable compatibility contract

1. Standard `simple_exec prg=abinitio2D` is byte-for-behavior unchanged.
2. Default `abinitio2D_sgd` behavior remains the current exhaustive hard
   stream until a new class-stability policy is explicitly selected.
3. The new policy is development-only, explicit opt-in, and default disabled.
4. Stages and iterations where streaming SGD is inactive remain unchanged.
5. The terminal conventional pass remains unchanged.
6. The particle mini-batch continues to use `sgd_update_frac` and the existing
   200,000 SIMPLE cap.
7. The raw Euclidean objective, its sign, the class/rotation winner rule, and
   the direct shift optimizer are not redefined.
8. Each sampled particle still produces exactly one hard class/rotation/shift
   assignment.
9. No SoftMax, top-K table, likelihood table, global candidate table, or
   second baseline run is added to production.
10. A never-viable or invalid class is not filled merely to satisfy a
    target class count.
11. The two-member floor may retain only particles whose incumbent class was
    that same previously viable class. It never imports an unrelated particle
    solely to revive a class.
12. All detailed reporting remains behind `sgd_diagnostic`, which is default
    off.

## Proposed v1: incumbent-aware frontier with a two-member support floor

### Public policy

Add one development-only policy selector after the isolated policy has passed
focused tests:

```text
sgd_class_policy=exhaustive|frontier
```

The default is `exhaustive`, which is the current behavior. `frontier` is
accepted only by `abinitio2D_sgd`; it is not shown by standard `abinitio2D`.

Do not add candidate-count or cadence UI parameters in the first code slice.
Keep those values in the owned policy object/test fixture until measurements
show that users need to tune them. This avoids creating permanent unused SGD
arguments while the algorithm is still experimental.

### Candidate construction

For sampled particle `i` on an ordinary frontier iteration, build

\[
C_i = \{c_i^{\mathrm{old}}\} \cup R_i,
\]

where:

- `c_i_old` is included when it was previously viable and its previous model
  is valid, even when its provisional current population is below two;
- `R_i` is a bounded, duplicate-free sample of currently viable classes;
- the incumbent is not counted twice when it is also viable; and
- if no incumbent model is valid, the frontier consists only of viable
  alternatives.

Candidate generation must be deterministic for a fixed project state,
iteration, user seed, and particle index, independent of OpenMP scheduling.
Do not share a mutable random-number generator across matcher threads. A
stateless per-particle permutation/hash or an equivalent thread-safe existing
SIMPLE API is required.

For each candidate class, v1 continues to evaluate the complete current
rotation grid with `gen_raw_euclid_vals`. This deliberately changes only the
class competition. A local rotation neighborhood may be studied later, after
class behavior is understood.

The winner remains

\[
(c_i,r_i)=\arg\min_{c\in C_i,r}L_i(c,r,s_i).
\]

The existing direct shift refinement then operates on that one winner.

### Two-member support floor

The frontier produces provisional hard winners. Before assignments are
committed, calculate the resulting class populations. For every previously
viable class whose provisional population is below two:

1. consider only particles whose incumbent class was that class;
2. rank the incumbents that proposed leaving by the finite raw-loss penalty

   \[
   \Delta_i = L_i(c_i^{\mathrm{old}})-L_i(c_i^{\mathrm{winner}});
   \]

3. retain the smallest-penalty incumbents needed to reach population two; and
4. commit one final hard assignment per particle.

This is a constrained hard assignment, not a probability or rescue bonus. It
never recruits a particle from an unrelated incumbent class. A class that was
never viable receives no anchors. If fewer than two active incumbent particles
remain available, the floor is impossible; the class is reported as dead and
no unrelated particle is fabricated to satisfy it.

The implementation must not build a particle-by-class table. It needs only
the provisional winner, incumbent loss, winner loss, and a bounded two-entry
selection per class. The resulting memory is linear in sampled particles plus
classes.

Class-average restoration must use the same viability rule. If an SGD class
has population below two and a previous viable model exists, retain the
previous model rather than constructing a singleton or zero image. The stored
population remains truthful; model preservation and support accounting are
separate facts.

### Exploration safety

Restricting every iteration forever can freeze a bad assignment. The policy
therefore needs periodic full exploration over all viable classes. The
cadence must:

- use global/restart-stable iteration state;
- remain independent of OpenMP thread count;
- be visible in diagnostics; and
- be frozen before matched scientific runs.

Under-supported classes are not globally advertised during full exploration.
They remain available only to particles for which they are a previously
viable incumbent. This is a conservative protection rule: the class can retain
its own best incumbents without receiving unrelated particles or a global
rescue bonus.

The first implementation must expose the cadence as an explicit constructor
argument to the policy helper and tests. Whether it becomes a user parameter
is a later decision based on sensitivity measurements.

### Why this is analogous to `prob_snhc`

The conventional probabilistic neighborhood machinery reduces candidate
competition through a sparse candidate set and later performs global
probabilistic assignment. The proposed frontier reuses only the useful
locality idea:

- evaluate a bounded changing subset rather than every class;
- always retain a relevant incumbent candidate; and
- periodically widen exploration.

It intentionally does not reproduce the probability table, SoftMax weights,
top-K storage, or global assignment pass. It is therefore a table-free
analogue, not a reimplementation of `prob_snhc`.

## Phase plan

### Phase 0: freeze and measure the current exhaustive stream

Before production behavior changes:

1. Record the exact source baseline and stage schedule.
2. Verify the 60-percent/200,000 sampling contract in the production caller.
3. Record and test the existing `< 2` search-viability rule, including the
   mismatch with the current zero-only restoration predicate.
4. Add diagnostic summaries for:
   - sampled particle count and fraction;
   - active, viable, under-supported, zero-population, and previously viable
     classes;
   - per-class population before and after the iteration;
   - particles that stayed in or left their incumbent class;
   - class reactivation count;
   - raw winner and runner-up loss margin without storing a top-K table; and
   - time spent in candidate evaluation and shift refinement.
5. Run the existing focused suite and one current exhaustive-stream workflow.

No new class policy is enabled in this phase.

### Phase 1: isolated candidate-policy helper

Create one owner-aligned helper in the 2D search strategy subsystem. It owns
candidate indices only; it does not evaluate images, losses, gradients, or
class averages.

Pure/focused tests must prove:

- incumbent inclusion;
- under-supported incumbent inclusion only when it was previously viable and
  its previous model is valid;
- bounded size and no duplicates;
- alternatives are viable classes;
- deterministic results for fixed seed/iteration/particle;
- different iterations change the alternatives;
- thread-order independence;
- full-exploration cadence and restart stability; and
- exact exhaustive behavior when the policy is disabled.

### Phase 2: provisional matcher integration behind explicit opt-in

Integrate the helper only into the active SGD branch of
`strategy2D_greedy`/its owning search state.

Required behavior:

1. `sgd_class_policy=exhaustive` follows the existing loop unchanged.
2. `sgd_class_policy=frontier` evaluates the incumbent-aware frontier.
3. Both modes call the same `gen_raw_euclid_vals` and the same direct shift
   minimizer.
4. All accepted orientations and stored score semantics remain unchanged.
5. Standard and SGD-inactive routes do not construct or call the helper.
6. Diagnostics identify frontier size, incumbent winner/stay/leave outcome,
   and full-scan iterations.

Do not commit assignments yet through a test-only shortcut. The next phase
owns the global two-member support constraint before final orientation writes.

### Phase 3: two-member floor and low-support restoration

Add the bounded per-class anchor selection after provisional hard assignment
and before final orientation/class-average updates.

Focused tests must prove:

- population zero and population one are both under-supported;
- population two is viable;
- only a previously viable class can receive anchors;
- only that class's own incumbents are eligible;
- the smallest finite raw-loss penalties are retained;
- at most two anchors are needed per class;
- final assignments remain unique and populations are recomputed exactly;
- an impossible floor reports class death rather than recruiting an unrelated
  particle; and
- class-average restoration retains the previous model for both population
  zero and population one while storing the truthful population.

### Phase 4: production-context stability fixture

Extend the SGD regression suite with a deterministic multi-class fixture that
runs multiple stochastic updates. It must include:

- at least one strong class;
- at least one weaker but valid class;
- a class that falls below two members under exhaustive competition;
- a preserved prior average for that class; and
- particles whose incumbent loss remains scientifically competitive.

Compare exhaustive and frontier policies from the same initial particle
states, references, noise, and random seed. Report every iteration, then fail
only after the complete comparison so all diagnostics are retained.

The fixture must prove:

- each sampled particle has one hard winner;
- non-sampled particles are not silently reassigned;
- preserved class images are finite and unchanged at populations zero and
  one;
- an under-supported previously viable incumbent remains a legitimate
  candidate;
- the final protected class has exactly the two best incumbent anchors when
  provisional support is below two;
- class support never increases through an unrelated rescue assignment;
- accepted shift steps remain finite, improving, and bounded; and
- default exhaustive results remain unchanged.

### Phase 5: matched workflow experiment

Fork every arm from the same preprocessing/project checkpoint and preserve the
same stage schedule. Use at least beta-galactosidase and apoferritin and more
than one fixed random seed.

Compare:

1. conventional `abinitio2D` reference;
2. current `abinitio2D_sgd` exhaustive hard stream; and
3. `abinitio2D_sgd` incumbent-aware frontier.

The conventional arm is scientific reference evidence, not a runtime
dependency of production SGD.

Record:

- normal completion and runtime;
- sampled particles per iteration;
- active, ranked, viable, under-supported, and zero-population classes per
  iteration;
- number of anchored classes and retained anchor particles;
- class-survival curve and reactivations;
- minimum, median, maximum, and entropy/effective class population;
- best resolution;
- median resolution of ranked classes;
- classes below 9 Angstrom and their median resolution;
- raw objective and winner-margin distributions;
- incumbent stay/leave rate;
- direct-shift attempts, accepted steps, final loss, and bound compliance; and
- output project consistency.

Freeze acceptance thresholds before reading final experiment results. The
frontier is acceptable only if it improves class survival consistently without
a material loss in class quality or an unacceptable runtime increase. A
single improved best class is not sufficient.

### Phase 6: decide whether stronger support control is necessary

The two-member floor is already part of v1 because it matches SIMPLE's actual
class-viability rule. If the frontier plus this floor is sufficient, stop. Do
not add more machinery.

If previously meaningful classes still collapse or remain permanently pinned
at two without improving, evaluate one stronger mechanism at a time.

#### Option: capped support-aware score adjustment

Only if anchors are insufficient, test a bounded class-support term such as

\[
L'_i(c)=L_i(c)+\lambda\,\mathrm{clip}\left(
\log\frac{n_c+\epsilon}{\bar n+\epsilon},-b,b
\right).
\]

For an under-supported class the logarithm is negative and lowers the adjusted
assignment loss. The unadjusted raw loss must still be retained for scientific
reporting and shift optimization.

This option changes the assignment objective and therefore requires stronger
justification, scale calibration, and matched validation. Its first experiment
must compare against the fixed two-member floor without silently changing the
anchor rule.

#### Rejected first-line mechanisms

- An unconditional bonus for every empty or singleton class.
- Forced equal class populations.
- A separate conventional baseline run during production.
- Global top-K/SoftMax tables.
- Class-balanced particle sampling presented as a collapse solution: it can
  balance the input batch by old labels, but particles can still all migrate
  after scoring.

## Deferred optimizer track

The former direct-joint optimizer plan is not the current work item.

Current evidence supports this boundary:

- the established joint L-BFGS-B route already consumes a discrete seed and
  performs bounded local continuous polishing;
- the attempted extra seeded wrapper was redundant in the production caller;
- `inpl_cont=yes` focused tests established finite analytical gradients,
  integer-route parity, accurate synthetic pose recovery, and operational
  polishing behavior; and
- a separate direct projected solver has not demonstrated a workflow benefit.

Future optimizer research may benchmark a bounded block-scaled projected-
gradient method with Armijo backtracking against L-BFGS-B on the same particle,
class, rotation seed, shift seed, bounds, and objective. It should advance only
if it achieves comparable final loss and pose error with fewer evaluations or
lower runtime across realistic data. Coordinate equality is not required in a
shallow local objective; scientific pose quality and final loss are.

This future track must remain orthogonal to class-frontier development. Do not
change the inner optimizer and the outer class policy in the same experiment.

## File ownership for the planned implementation

Expected owners, subject to a fresh call-trace before each phase:

- `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90`:
  preserve the existing SGD mini-batch contract; no new cap.
- `src/main/strategies/search/simple_strategy2D_matcher.f90`:
  active-stream orchestration, provisional/final assignment boundary, anchor
  application, and aggregate diagnostics.
- `src/main/strategies/search/simple_strategy2D_greedy.f90`:
  consume the class frontier while retaining the existing raw-loss and shift
  calls.
- a small owner-aligned 2D search policy module or existing allocation/state
  module: deterministic candidate generation and bounded two-anchor selection.
- `src/main/strategies/search/simple_strategy2D_alloc.f90`:
  authority for the existing population-below-two viability rule; replace the
  literal threshold with one shared 2D viability constant only if all current
  callers and tests are audited together.
- `src/main/class/simple_classaverager_restore.f90`:
  align restoration with minimum viable support by preserving the previous
  valid model at populations zero and one.
- `src/main/params` and `src/main/ui/simple/simple_ui_cluster2D.f90`:
  add only the final explicit development policy after the internal helper is
  validated.
- `production/tests/simple_sgd_base_suite_sampling_restore_test.f90` and the
  split SGD suite modules: pure and production-context regressions.

Do not modify standard `simple_commanders_abinitio2D.f90` to implement this
policy. The development wrapper and existing internal delegation remain the
activation boundary.

## Validation and workflow rules

- Follow `.codex/windows-msys-install-workflow.md` and
  `.codex/compile_windows.sh` for the local Windows/MSYS build, with at most 12
  build jobs.
- Use targeted source checks before asking for a long build.
- The user will run `simple_test_sgd_base_suite` and scientific workflows on
  the Oracle Linux box.
- Do not claim Linux/BOX validation without observed output.
- A shared extracted project is necessary but not sufficient for a
  deterministic A/B test. Record and control the assignment/random seed in
  each arm.
- Checker scripts must prove attempted work, finite losses, class-support
  behavior, and output consistency, not merely find an activation marker.
- Diagnostics must not change default numerical behavior and must remain off
  in production unless explicitly requested.

## Commit boundaries

Keep implementation commits independently reviewable:

1. baseline diagnostics, shared minimum-support contract, and
   production-context low-support restoration test;
2. isolated deterministic class-frontier helper and pure tests;
3. opt-in provisional matcher/greedy integration and route-identity tests;
4. bounded two-member anchor selection and final-assignment tests;
5. deterministic multi-iteration class-stability fixture;
6. workflow scripts/reporting, if changed;
7. optional support-aware scoring only after a separate approved decision.

Do not mix direct-joint optimizer work, standard `abinitio2D` changes,
continuous-3D work, or unrelated test development into these commits.

## Completion criteria

This plan is complete only when:

1. default standard and development paths retain their established behavior;
2. the explicit frontier mode is deterministic, table-free, and produces one
   hard winner per sampled particle;
3. fresh mini-batch sampling and the existing 200,000 cap remain unchanged;
4. the population-below-two death rule is represented by one tested policy
   contract;
5. previously viable classes retain their previous model at populations zero
   and one;
6. the two-member floor retains only the best same-class incumbents and never
   fabricates an unrelated rescue assignment;
7. matched multi-seed workflows show improved class survival on at least two
   datasets without material class-quality or runtime regression;
8. diagnostics are default-off and all focused tests pass; and
9. the result is reviewed before any default is changed.

Until those gates pass, `sgd_class_policy=frontier` remains an explicitly
development-only experiment and `exhaustive` remains the default stream
policy.
