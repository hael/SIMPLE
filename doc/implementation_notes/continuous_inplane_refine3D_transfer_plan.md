# Continuous in-plane rotation transfer plan for `refine3D`

> **Current policy (2026-08-07).** This file is a historical implementation
> record. The earlier `refine=shc`-only gate, rejection of probabilistic modes,
> and invalid-result callback fallback described in later phase notes are
> superseded. The authoritative policy is
> `doc/policies/refine3D_policy.md`: `inpl_cont=no` uses the legacy callback;
> eligible `inpl_cont=yes` replaces every callback use with joint
> `(sx,sy,rotind_frac)` optimization, including probabilistic candidate
> profiling, and never falls back to the callback. Each joint solve discards
> incoming in-plane state, selects the best discrete cell once at the supplied
> native shift, and then optimizes within plus or minus two cells. The continuous
> route does not use the legacy 5-by-5 shift/all-angle coarse initializer or a
> persisted fractional restart seed.

## Project summary and quick start

This project adds an opt-in continuous in-plane refinement step to the
standalone `refine3D` workflow. The existing search still selects the discrete
state and projection direction. For the selected winner, the new route jointly
refines the two image-plane shifts and in-plane angle,
$(s_x,s_y,\theta)$, with the raw Euclidean objective. The historical
grid-angle route remains the default.

Enable the feature by supplying `inpl_cont=yes` together with
`objfun=euclid`:

```bash
simple_exec \
  prg=refine3D \
  projfile=/absolute/path/input.simple \
  mkdir=yes \
  objfun=euclid \
  inpl_cont=yes \
  maxits=1 \
  mskdiam=180 \
  nparts=10 \
  nthr=8
```

Omit `inpl_cont`, or set `inpl_cont=no`, to use the unchanged legacy path.
The opt-in route requires raw Euclidean matching: `objfun=euclid`,
`objfun_den=no`, `ptcl_src=raw`, and `projrec=no`. There is no `refine=shc`
restriction. Deterministic, neighborhood, evaluation, and probabilistic
matcher routes follow the same callback-replacement policy wherever they
perform pose search; modes with no pose search invoke neither optimizer.

Probability tables keep their established artifact format. They carry the
rounded in-plane index, score, and rounded-frame shift while the fractional
coordinate remains transient. After the final hard assignment, refine3D
uses the stored rounded index only to recover the native shift, discards all
prior in-plane state, selects the best discrete cell once at that shift, and
reruns the joint optimizer. It persists fractional `e3`, integer `inpl`, shift,
and score from one pose. Valid non-improving work commits the newly selected
grid seed; invalid work retains the incoming candidate. Neither invokes the
callback.

The tests are organized under `production/tests`. Build and invoke the mother
test from an already configured build tree:

```bash
repo=/absolute/path/hael_SIMPLE
build="$repo/build"

cmake --build "$build" \
  --target simple_test_continuous_inplane_refine3D \
  --parallel 12

cd "$build"
simple_test_continuous_inplane_refine3D
```

The parameter-free command runs the self-contained policy and search-state
groups and reports the project-backed baseline and synthetic recovery groups as
skipped. Supply the project, expected baseline, and volume inputs to run all
seven groups:

```bash
simple_test_continuous_inplane_refine3D \
  projfile=/absolute/path/refine3D_output.simple \
  snapshot="$PWD/refine3D_current_baseline.tsv" \
  expected_snapshot=/absolute/path/refine3D_legacy_baseline.tsv \
  vol1=/absolute/path/reference_volume.mrc
```

Create the expected baseline once from the accepted legacy project before
running the comparison:

```bash
simple_test_continuous_inplane_refine3D \
  case=baseline \
  projfile=/absolute/path/legacy_refine3D_output.simple \
  snapshot="$PWD/refine3D_legacy_baseline.tsv"
```

An individual group can be rerun with `case=<name>`, for example
`case=policy_suite`, `case=direct_route`, `case=synthetic_recovery`, or
`case=metadata_project`. The files named
`simple_continuous_inplane_refine3D_*.f90` are helper modules compiled into the
mother executable; they are not separate test executables.

---

> **Historical record boundary.** The remaining sections preserve the staged
> transfer history and validation observations. Any statements below that
> limit activation to `refine=shc`, reject `refine=prob*`, or prescribe a
> legacy callback fallback are descriptions of superseded intermediate policy.

## Plan Scope

This plan transfers the validated 2D continuous in-plane refinement design from
`abinitio2D` to the standalone `refine3D` workflow.

The scope is deliberately limited to:

```text
simple_exec -> prg=refine3D
```

This phase does not modify `abinitio3D`, `abinitio3D_cavgs`, SGD, or the CC
objective. The default behavior must remain unchanged.

## Final completion audit

Core feature status on 2026-08-07: **complete and verified on Oracle Linux**.
The subsequent P2/P3 hardening corrections are also complete and verified by
the Oracle acceptance run recorded in the peer-review follow-up section.

Phases 3D-0 through 3D-5 are implemented. The policy, search-state, joint-state,
direct-route, metadata, restart, default-off, opt-in, and multi-iteration
real-data checks described below have passed on Oracle Linux. A diagnostic
one-iteration run also accounted for all 3,081 active particles: 2,839 joint
refinements improved the objective, 242 were finite but non-improving, and none
entered the legacy fallback.

The later initialization correction supersedes the restart-seed and
selected-angle retention behavior measured in that audit. For the common
single-state search, `xy_first` remains the native shift seed, but the joint
route now discards the incoming integer/fractional angle and selects the best
discrete cell once at `xy_first`. A finite non-improving solve stores that new
grid cell, score, and rotated shift with false continuous validity. Persisted
fractional metadata is not reused as a restart seed, and invalid numerical work
retains the incoming candidate without requesting the legacy callback. The
focused `direct_route` and route-identity contracts cover the new behavior.

The planned truth-controlled 3D recovery fixture is implemented as the named
`synthetic_recovery` case. It fixes the 3D projection direction,
constructs a particle with a known off-grid in-plane angle and known image-plane
shift, and checks finite evaluation, non-worsening raw Euclidean loss, improved
angular recovery, and bounded shift recovery through the production PFTC joint
optimizer. Its first Oracle run exposed a fixture mismatch: starting the joint
optimizer at zero shift reduced the objective but converged to an
angle/translation-coupled local minimum. Production `refine3D` does not start
there. It first runs the legacy alternating angle/shift preselection search and
uses the resulting native-frame `xy_first` to score the discrete candidates.
The fixture now reproduces that production seed with an explicit five-pixel
translation range before invoking the joint optimizer. The corrected fixture,
`direct_route`, and all seven mother-suite groups subsequently passed on Oracle
Linux.

The final policy audit found that `inpl_cont=yes` still accepted the non-search
`refine=eval` and `refine=sigma` modes. `eval` could reach the continuous pose
update even though it should only evaluate a pose; `sigma` accepted the option
without performing refinement. Both modes are now rejected in the commander
command-shape gate and the typed parameter-consistency gate. Focused
`policy_eval` and `policy_sigma` rejection cases passed on Oracle Linux together
with the eligible default/off/on contract and every earlier rejection case.
The first post-audit smoke attempt exposed one additional no-improvement-route
defect. `minimize_joint` preserves its legacy API by returning `irot=0` when a
valid evaluation finds no strict improvement. The new seed-retention branch
then used that output value in `get_rot`, producing an out-of-range failure on
multiple distributed partitions. `refine_selected_continuously` now preserves
the selected grid index before calling the optimizer and uses that saved index
only to rotate and retain the non-improving seed. The optimizer API and its
legacy callers remain unchanged.

The corrected Oracle smoke completed and accounted for all 3,081 active
particles: 2,830 improved, 251 produced a finite no-improvement result, and
zero entered fallback. The totals balanced exactly, proving that the corrected
no-improvement branch executed in production. Output
`27_refine3D/betagal.simple` then passed the metadata contract with 3,081 active
particles, 1,930 continuous angles, and zero invalid indices, non-finite values,
or integer/angle mismatches on the 440-angle grid.

Those outputs completed the original feature acceptance. The peer-review
hardening added afterward must pass its targeted Oracle checks before the
current worktree is committed. Another long multi-iteration run is not
required unless the short smoke checks expose a regression.

## Peer-review follow-up: P2 and P3 findings

The 2026-08-07 peer review of commit `6e8b0d521` found no P0 or P1
correctness blocker. It identified the following hardening and maintainability
work. The detailed review is recorded in the Obsidian report
`302 Executive peer review and code walkthrough - continuous in-plane refine3D.md`.

### P2 findings to resolve first

1. **Default-off continuous-state overhead.** The non-probabilistic 3D search
   currently allocates, resets, and writes fractional in-plane coordinates and
   validity flags even when `inpl_cont=no`. The numerical legacy result is
   unchanged, but the default route pays continuous-only memory and hot-path
   write costs. Allocate and update these arrays only for the opt-in route.
2. **Incomplete joint-result validity check.** `minimize_joint` verifies that
   its initial and final costs are finite, but it does not verify that the
   returned optimizer coordinate vector is finite and inside its active
   bounds before using `nint` and returning metadata. Invalid coordinates must
   set `evaluation_valid=.false.` and enter the existing legacy fallback.
3. **Enabled deterministic-mode surface exceeds production validation.** The
   gate rejects known unsupported modes but otherwise admits deterministic
   modes beyond the distributed `refine=shc` route validated with biological
   data. Until those modes and shared-memory execution are tested, restrict the
   public opt-in contract to `refine=shc` and reject other explicit modes with
   a clear message.
4. **Baseline helper is not a comparator.** The baseline test validates a
   project and replaces a TSV snapshot, but it does not compare the snapshot
   with an expected result. Add an explicit expected-snapshot comparison mode
   with exact integer checks and documented floating-point tolerances.

### P3 findings and remediation

1. The public commander and typed parameter layer duplicated the compatibility
   policy, while the focused negative tests directly exercised only the
   commander gate. **Resolved locally:** both layers now call the same pure
   `refine3D_inpl_cont_policy_error` decision function owned by
   `simple_parameters`; the commander only gathers effective command values,
   and typed worker validation supplies resolved parameter values.
2. The reviewed commander refinement-mode block called `value%kill` on fatal
   branches and then contained a later test of the same value. The P2
   refinement-mode allowlist removed this pattern incidentally; it remains
   recorded here as a resolved peer-review observation.
3. Rejection tests used fixed log filenames in the build directory. Concurrent
   mother-suite runs could overwrite one another. **Resolved locally:** every
   rejection log now includes the mother-suite process ID and case label.
4. Historical pending statements remained in this plan after final Oracle
   acceptance, and the definition of done inherited an `abinitio2D` coarse
   stage split that standalone `refine3D` does not have. **Resolved locally:**
   completion language now records both accepted core behavior and completed
   peer-review validation, and the final criteria describe per-iteration
   `refine3D` activation rather than early and late stages.

P2 remediation is tracked as part of the feature before the implementation is
described as hardened beyond the validated distributed `refine=shc` route.

### P2 remediation status

Implementation status on 2026-08-07: all four P2 corrections are implemented
and verified by Oracle Linux compilation, focused contracts, the complete
mother suite, and matched production smoke settings.

1. **Default-off continuous-state overhead: implemented.**
   `simple_strategy3D_alloc.f90` now allocates the fractional-coordinate and
   validity arrays only for `inpl_cont=yes`. `store_solution` preserves the
   historical score, shift, and integer-angle writes when those optional
   arrays are absent. The search-state contract now exercises this default-off
   storage route and verifies that it does not create continuous-only state.
2. **Joint-result coordinate validity: implemented.**
   `simple_pftc_shsrch_grad.f90` now requires every returned optimizer
   coordinate to be finite and inside the active L-BFGS-B limits, allowing a
   small floating-point tolerance. A violation sets
   `evaluation_valid=.false.` and therefore uses the established legacy
   fallback in `simple_strategy3D_srch.f90`.
3. **Deterministic-mode surface: implemented.** Both the public commander
   gate and typed parameter validation now permit the opt-in route only with
   base `refine=shc`. The strategy repeats this condition as defense in depth.
   `policy_neigh` was added alongside the existing `eval`, `sigma`, and
   probabilistic rejection cases.
4. **Baseline comparison: implemented.** The baseline helper retains an
   explicit generate-only mode and adds `expected_snapshot`, `snapshot_atol`,
   and `snapshot_rtol`. Comparison mode checks the header, row count, integer
   fields, finite floating values, and floating values within the documented
   absolute-plus-relative tolerance. The mother suite runs this group only
   when both the project and expected snapshot are supplied.

Observed acceptance:

```text
policy suite: all eligible and nine rejection cases passed
search state: default-off allocation and candidate storage passed
baseline: generate-only and exact comparison passed for 3,081 active particles
mother suite: 7/7 groups passed, zero skipped, zero failed
default-off smoke: 3,081 active, 0 continuous, 0 invalid/non-finite/mismatch
opt-in smoke: 3,081 attempted, 2,848 improved, 233 no-improvement, 0 fallback
opt-in metadata: 1,937 continuous, 0 invalid/non-finite/mismatch
```

The second P3 finding was removed incidentally by the P2 mode allowlist: the
commander now reads each command value before `value%kill` and never inspects a
killed `string`. All four P3 corrections passed the same Oracle acceptance
boundary as the P2 work.

## Existing 2D contract to preserve

The public option is:

```text
inpl_cont=no|yes
```

`no` is the default. It preserves the historical grid-angle path.

`yes` is an explicit opt-in. It applies only to an eligible Euclidean search.
The search first selects a discrete candidate. It then jointly refines

$$
\mathbf{x} = (s_x, s_y, \theta).
$$

Here $s_x$ and $s_y$ are the two-dimensional particle-image shifts used by
the polar Fourier matching objective. They are not the three components of a
free-space 3D translation.

The integer rotation index remains available for legacy consumers. The
continuous angle is retained separately in orientation metadata.

The 3D transfer must preserve this same sequence:

```text
discrete candidate selection
    -> optional continuous joint refinement
    -> legacy integer metadata plus continuous angle metadata
```

## Current `refine3D` path

```text
simple_exec
  -> simple_exec_refine3D.f90
  -> simple_commanders_refine3D.f90
  -> simple_refine3D_strategy.f90
  -> simple_strategy3D_matcher.f90
  -> simple_strategy3D_* search strategy
  -> simple_strategy3D_srch.f90
  -> simple_pftc_shsrch_grad.f90
  -> simple_polarft_corr.f90
```

The numerical integration point is `simple_strategy3D_srch.f90`.
The commander and parallel strategy should carry policy and lifecycle state;
they should not implement the angular mathematics.

The current 3D search already performs continuous optimization of the
two-dimensional image shift through `pftc_shsrch_grad`, but it stores only an
integer in-plane rotation index. The transfer therefore extends the existing
3D search result rather than replacing the 3D search architecture.

The orientation object also exposes `get_3Dshift()`, and it can contain a
`z` field. That does not mean that `refine3D` currently optimizes $s_z$ in the
particle-matching loop. The current search calls `get_2Dshift()` and passes a
two-component shift to the PFTC objective. A viewing-direction translation is
not an independently identifiable image shift in this projection-matching
step. It must not be added to the continuous optimizer unless a separate 3D
objective and a clear metadata/update contract are designed.

## Exact transfer boundary

The transfer follows the same candidate-selection boundary as `abinitio2D`.

In 2D, class identity and the initial angle index are selected discretely. The
selected class is kept fixed, while the angle and two image shifts are jointly
refined as

$$
\mathbf{x}_{2D}=(s_x,s_y,\theta).
$$

In `refine3D`, state identity, projection direction, and the initial in-plane
angle index are selected discretely. The selected state and projection
direction remain fixed, while the in-plane angle and two image shifts are
jointly refined as

$$
\mathbf{x}_{3D}=(s_x,s_y,\gamma).
$$

The existing PFTC reference is already a two-dimensional reprojection for the
selected 3D projection direction. Therefore this three-variable polish is the
direct, testable transfer of the validated 2D algorithm.

Full continuous viewing-direction or five-parameter pose refinement is outside
this plan. It would require a new arbitrary-orientation central-slice evaluator
and derivatives with respect to the 3D rotation. Those are not part of the
current PFTC reference-table contract and are not needed for continuous
in-plane rotation in `refine3D`.

## Phase 3D-0: establish the legacy baseline

Implementation status: complete and verified on Oracle Linux. No `refine3D`
production file is changed by this phase.

The verified project was the interrupted real-data beta-galactosidase
`refine3D_auto` checkpoint after two fully assembled iterations and the third
distributed matcher pass. The baseline reader reported:

```text
total particles:       4231
active particles:      3081
invalid indices:          0
non-finite records:       0
correlation sum:       1.1696115581989288E+03
squared-shift sum:     6.8470346561560500E+05
```

The run ended with `SIMPLE_TEST_CONT_INPL_REFINE3D_BASELINE NORMAL STOP`.

Before changing code, run `refine3D` with the current checkout and record:

- normal completion;
- output orientation metadata;
- integer in-plane indices;
- shifts and objective values;
- refinement-mode and iteration behavior;
- shared-memory and distributed execution behavior, where available.

The baseline must use the default with no continuous option supplied.

Acceptance condition: the current `refine3D` result is reproducible and the
existing grid-only path is understood before the opt-in path is introduced.

The registered mother test is:

```text
simple_test_continuous_inplane_refine3D
```

A parameter-free mother-suite run skips the project-dependent baseline with a
clear notice and continues with every self-contained test group. Supplying both
`projfile` and `expected_snapshot` enables baseline comparison and includes it
in the pass/fail summary. A focused `case=baseline` invocation without
`expected_snapshot` generates and validates a candidate snapshot but clearly
reports `MODE: GENERATE_ONLY` rather than claiming a comparison.

Phase-specific test implementations use the non-registered dependency naming
pattern `simple_continuous_inplane_refine3D_*.f90`. The Phase 3D-0 dependency
is `simple_continuous_inplane_refine3D_baseline.f90`. New phase tests will be
added to the mother test as named child cases.

Generate the expected snapshot from an accepted legacy `refine3D` invocation:

```text
simple_test_continuous_inplane_refine3D \
  "case=baseline" \
  "projfile=/absolute/path/to/refine3D_output.simple" \
  "snapshot=refine3D_legacy_baseline.tsv"
```

Compare a current default-off project against that expected snapshot with:

```text
simple_test_continuous_inplane_refine3D \
  "case=baseline" \
  "projfile=/absolute/path/to/current_refine3D_output.simple" \
  "snapshot=refine3D_current_baseline.tsv" \
  "expected_snapshot=refine3D_legacy_baseline.tsv"
```

The snapshot records one deterministic row per active `ptcl3D` orientation:

```text
particle, state, projection index, in-plane index,
Euler angles, two-dimensional shift, stored objective
```

The checker rejects missing or non-positive legacy indices and non-finite pose
or objective values. It also prints aggregate correlation and squared-shift
sums for a quick run-to-run comparison. In comparison mode it requires exact
particle, state, projection, and integer in-plane indices. Euler angles,
shifts, and scores use configurable absolute and relative tolerances, both
defaulting to $10^{-6}$.

The normal workflow log remains the authority for iteration count, refinement
mode, normal completion, and shared-memory or distributed execution behavior.

## Phase 3D-1: add the policy gate

Implementation status: complete and verified on Oracle Linux. No continuous
`refine3D` numerical path is enabled by this phase.

The original focused policy suite passed the eligible default/off/on contract
and eight rejection cases. The P2 hardening update adds `policy_neigh`, making
the current rejection matrix nine cases:

```text
policy
policy_bad_value
policy_cc
policy_hybrid
policy_denoised
policy_projrec
policy_eval
policy_sigma
policy_probabilistic
policy_neigh
```

The run ended with `Continuous in-plane refine3D policy suite: PASS`.

Add `inpl_cont=no|yes` to the `refine3D` interface and parameter propagation.
The default must be `no`.

The gate must reject unsupported combinations before search setup. Initially:

```text
inpl_cont=yes and objfun=cc      -> reject
inpl_cont=yes and objfun=euclid -> eligible for further checks
```

The gate must also preserve the current `refine3D` contracts for denoised or
hybrid objectives, probabilistic assignment modes, and `projrec=yes`. Each
route must be declared supported or excluded according to whether it provides
the same raw Euclidean objective and metadata semantics as the joint polish.

The initial hardened gate deliberately supports only raw-Euclidean
`refine=shc` matching with `projrec=no`. It rejects `objfun=cc`,
`objfun_den=yes`, `ptcl_src=den`, `projrec=yes`, and every explicit refinement
mode other than `shc` when `inpl_cont=yes`. Omitting `inpl_cont`, or setting it
to `no`, does not apply these restrictions.

The focused Phase 3D-1 test does not require a project:

```text
simple_test_continuous_inplane_refine3D case=policy_suite
```

Acceptance conditions:

- omitted `inpl_cont` behaves exactly like `inpl_cont=no`;
- explicit `inpl_cont=no` behaves exactly like the omitted option;
- unsupported objective combinations fail clearly;
- no `inpl_cont` value changes 2D or SGD behavior.

## Phase 3D-2: extend the 3D search state

Implementation status: complete and verified on Oracle Linux. The new state is
not consumed by orientation metadata or any numerical refinement path in this
phase.

The focused search-state test verified grid seeding, real `store_solution`
integration, rejection preservation, improving-candidate replacement, and
continuous-validity reset. It reported:

```text
REFINE3D_INPL_CONT_STATE GRID_SEED/STORE/REPLACE: PASS
```

Extend the per-particle and per-reference search state used by
`simple_strategy3D_srch` so that a selected candidate can carry:

- the legacy integer in-plane index;
- the continuous rotation coordinate;
- the refined shifts;
- the reference/projection identity associated with the candidate;
- a validity flag for the continuous result.

The continuous angle must remain attached to the selected 3D reference and
particle. It must not be stored as a global or shared value.

The integer index remains the fallback whenever refinement is rejected,
non-finite, non-improving, or outside the bounded search interval.

The implementation adds a fractional grid coordinate and validity flag beside
the existing per-reference integer in-plane index. Every legacy
`store_solution` call seeds the fractional coordinate from the integer index
and clears validity. Allocation, per-thread reset, and teardown include both
new arrays.

The focused Phase 3D-2 test does not require a project:

```text
simple_test_continuous_inplane_refine3D case=search_state
```

## Phase 3D-3: add joint refinement after candidate selection

Implementation status: complete and verified on Oracle Linux with the
focused joint-state test and a one-iteration distributed production-context
opt-in `refine3D` run.

The implementation constructs the existing raw-Euclidean joint PFTC optimizer
only when `inpl_cont=yes` passes the Phase 3D-1 policy. Immediately before the
final orientation is assigned, it fixes the selected state and projection
reference and refines only $(s_x,s_y,\theta)$. The accepted result updates the
selected reference's shifts, nearest legacy integer index, continuous grid
coordinate, validity flag, and score. The existing `minimize_joint` acceptance
contract rejects non-finite or non-improving continuous updates. The
post-audit correction additionally ensures that a non-improving grid candidate
retains the same selected-angle shift used to compute its discrete score.

At this phase the accepted continuous coordinate remains in the 3D search
state. Final project metadata still uses the legacy integer angle until Phase
3D-5 adds the explicit continuous-angle assignment and restart contract.

The focused state-commit case is:

```text
simple_test_continuous_inplane_refine3D case=joint_state
```

The observed focused result was:

```text
REFINE3D_INPL_CONT_JOINT SELECTED_REFERENCE/STORE: PASS
```

The parameter-free mother suite also passed all three self-contained groups.
It skipped the baseline with a clear missing-`projfile` notice and reported
three groups run, one skipped, three passed, and zero failures. The full suite
with the beta-galactosidase project passed all four groups.

The first distributed production smoke run exposed a command-propagation
defect before matching began. The parent `inpl_cont=yes` key was copied into an
internal `reconstruct3D objfun=cc` startup command, where the global objective
gate rejected it. The refine3D strategy now removes matcher-only arguments
from reconstruction, sigma, power-spectrum, probabilistic, assembly, and
postprocessing child commands while retaining `inpl_cont=yes` in the parent
and scheduled matcher workers.

After that correction, the policy test passed and the beta-galactosidase smoke
run completed one distributed `refine3D` iteration with
`objfun=euclid inpl_cont=yes`, `nparts=10`, and 3,081 active particles. The run
completed sigma estimation, matching, volume assembly, and postprocessing. It
reported 100% of active particles sampled and updated, then ended with:

```text
**** SIMPLE_REFINE3D NORMAL STOP ****
```

This establishes production execution and child-workflow isolation for Phase
3D-3. Continuous-angle project metadata and restart persistence remain Phase
3D-5 work.

The 3D search should follow the 2D control flow:

```text
1. Search the 3D projection/reference candidates on the existing grid.
2. Select the best discrete candidate.
3. Search the two image-plane shifts using the existing 3D objective path.
4. If inpl_cont=yes and the activation policy permits it, seed theta from the
   selected integer in-plane index.
5. Evaluate the joint objective and gradient for (sx, sy, theta).
6. Apply bounded L-BFGS-B refinement to $(s_x,s_y,\theta)$.
7. Accept only finite, improving results.
8. Retain the integer index and store the continuous angle separately.
9. Fall back completely to the legacy candidate on rejection.
```

The existing `minimize_joint` interface should be reused only if its objective
and derivative semantics are valid for the 3D reference projection. If the
3D objective requires a different evaluator, add that evaluator in the PFTC
owner module rather than duplicating numerical logic in the matcher.

## Phase 3D-4: activation policy

Implementation status: complete and accepted. The focused policy suite and
both opt-in and default-off one-iteration production workflows passed on
Oracle Linux.

The implementation does not introduce an `abinitio2D`-style stage number.
Omitted `inpl_cont` and explicit `inpl_cont=no` leave
`strategy3D_srch%continuous_active` false. Explicit `inpl_cont=yes` constructs
and invokes the joint optimizer only for raw-Euclidean `refine=shc` matching.
State selection, projection-direction selection, symmetry handling, and the
initial integer in-plane candidate remain discrete. Unsupported CC, hybrid,
denoised, projection-reconstruction, and probabilistic combinations are
rejected before strategy setup. Matcher-only arguments are removed from all
non-matcher child workflows.

Base `prg=refine3D` owns an iteration loop, not the staged coarse-to-fine
schedule used by `abinitio2D`. This phase must not invent a stage number or
reuse the `abinitio2D` stage-3 rule.

The intended policy is:

```text
inpl_cont omitted or no -> historical refine3D path
inpl_cont=yes and raw-Euclidean refine=shc search
                        -> discrete search followed by joint polish
unsupported objective or search mode -> reject or retain legacy path by an
                                        explicit documented rule
```

Probabilistic pre-alignment, hard assignment, symmetry selection, state
selection, and global projection-direction exploration remain discrete. A
probabilistic candidate may receive deterministic polishing only after its
assignment contract and stored objective semantics are verified separately.

The default-off run omitted `inpl_cont`, processed all 3,081 active particles,
completed reconstruction, matching, assembly, and postprocessing, and ended
with `SIMPLE_REFINE3D NORMAL STOP`. Its total runtime was 414.8 seconds,
compared with 420.1 seconds for the opt-in run. The observed 5.3-second or
approximately 1.3% difference does not establish a performance regression and
is small enough to include normal run-to-run variation. A later scientific
off/on comparison should still use matched copies of the same checkpoint and
schedule.

## Phase 3D-4.5: direct joint refinement with legacy fallback

Implementation status: complete and verified on Oracle Linux. Focused routing,
the parameter-free mother suite, and matched-source opt-in/default-off
production workflows passed. No commit has been created.

The opt-in route must replace, rather than follow, the legacy post-selection
alternating optimizer. The pre-selection shift seed remains part of the
existing discrete candidate search. After that search, the control flow is:

```text
global discrete state/projection/angle selection
    -> select the final discrete winner
    -> directly refine (sx,sy,theta) for that fixed state/projection
    -> accept a finite improving joint result
    -> retain the discrete winner when a valid joint run finds no improvement
    -> run the legacy post-selection alternating optimizer only when the
       joint evaluation is numerically invalid
```

The joint route is considered numerically invalid when it cannot produce
finite initial and final objective values or when its required seed is invalid.
A finite joint run that finds no strict improvement is complete and successful:
it retains the discrete score, integer angle, and grid-validity state, commits
the shift that was used to score that winner, and does not invoke the legacy
optimizer. This no-improvement policy may be revisited after advisor review.

The default-off route must continue to call the original post-selection
alternating optimizer without any new branch or changed numerical input. The
opt-in fallback must call that same legacy implementation with the untouched
discrete winner. It must not duplicate the legacy algorithm.

Acceptance conditions:

- `inpl_cont=no` follows the historical post-selection path;
- `inpl_cont=yes` does not run the legacy post-selection optimizer before the
  joint attempt;
- a finite improving joint result is committed;
- a finite non-improving joint result retains the discrete seed without
  invoking legacy refinement;
- an invalid joint evaluation invokes the legacy optimizer once for the fixed
  final reference;
- state and projection identity remain fixed throughout joint refinement and
  fallback.

The implementation keeps the existing pre-selection shift search because its
result is used to score the discrete reference/angle candidates. For
`inpl_cont=yes`, the ordinary `inpl_srch` and `inpl_srch_peaks` calls return
before invoking the legacy post-selection optimizer. The final winner then
enters `refine_selected_continuously` directly. Its joint seed is the native
pre-selection shift when discrete scoring used that shift; otherwise it is the
selected candidate's stored reference-frame shift converted back to the native
particle frame.

`minimize_joint` now reports numerical validity separately from strict
improvement. A valid non-improving result preserves the discrete candidate and
couples it to the shift used by its score. An invalid result calls the unchanged
legacy `inpl_srch` routine with an explicit fallback override for the fixed
final reference. Default-off calls continue through the ordinary legacy route.

The focused routing contract is registered through the mother suite as:

```text
simple_test_continuous_inplane_refine3D case=direct_route
```

It checks default-off legacy retention, opt-in legacy bypass, valid
no-improvement completion without fallback, score/shift consistency, and
invalid-result fallback selection. After the post-audit update, the expected
focused result is:

```text
REFINE3D_INPL_CONT_DIRECT BYPASS/NO_IMPROVEMENT/SCORE_SHIFT/FALLBACK: PASS
```

The parameter-free mother suite also passed its four self-contained groups:
policy, search state, joint state, and direct routing. It skipped only the
project-dependent baseline and reported five groups scheduled, four run, one
skipped, four passed, and zero failures.

Production-context validation used two one-iteration beta-galactosidase
`refine3D` executions from the same source project and with identical workflow
arguments except for the opt-in key. The `inpl_cont=yes` run completed normally
in 405.2 seconds. The default-off run, with `inpl_cont` omitted, completed
normally in 411.1 seconds. Both ended with:

```text
**** SIMPLE_REFINE3D NORMAL STOP ****
```

The opt-in observation was 5.9 seconds, or approximately 1.4%, faster. One
pair of runs is insufficient to establish a performance improvement because
search ordering and system load can vary. It does show that replacing the
legacy post-selection optimizer did not introduce an obvious runtime penalty
in this workflow. Repeated profiling may be performed later without blocking
Phase 3D-5.

## Phase 3D-5: metadata and restart behavior

Implementation status: complete and verified on Oracle Linux. Focused state,
opt-in metadata, default-off metadata, and opt-in restart metadata checks
passed. No commit has been created.

After each accepted update, verify that the 3D orientation metadata contains:

- the correct projection direction;
- the correct state assignment;
- the legacy integer in-plane index;
- the continuous angle in the project convention;
- the refined shifts.

Restarted `refine3D` runs must preserve the continuous angle when the opt-in
mode is enabled. Default-off runs must remain grid-aligned and restart exactly
as before.

The durable representation keeps two coupled values:

- `inpl` remains the nearest integer rotation-grid index for legacy consumers
  within the producing iteration;
- Euler `e3` stores the accepted continuous in-plane angle when
  `inpl_cont=yes` is active.

The default-off assignment is unchanged:

```text
inpl = selected integer grid index
e3   = 360 - grid_rotation(inpl)
```

The opt-in assignment converts the accepted fractional rotation coordinate to
the project Euler convention while retaining its nearest integer index. On
restart, `refine3D` may rebuild the PFTC at a different crop or resolution, so
the previous integer `inpl` is not authoritative on the new grid. Both the
current integer seed and fractional coordinate are reconstructed from `e3`
using the current PFTC spacing. The coordinate is unwrapped to the periodic
copy nearest that current integer seed, so a coordinate near the end of one
period can safely enter local bounds around index `1` through its equivalent
coordinate near `1`.

If the restarted particle selects the same state, projection, and integer
in-plane index, the recovered continuous coordinate becomes the joint seed. A
valid joint evaluation that finds no improvement retains this continuous
restart pose. A changed winner starts from its newly selected grid seed. An
invalid joint evaluation still enters the Phase 3D-4.5 legacy fallback.

The self-contained metadata and restart contract is registered as:

```text
simple_test_continuous_inplane_refine3D case=metadata_state
```

The project-backed contract is available as `case=metadata_project`. It checks
active particles for finite pose values, valid integer indices, nearest-index
agreement between `inpl` and `e3`, grid-only metadata in default-off output,
and continuous metadata in opt-in output. The producing grid must be supplied
explicitly as `dang=...` or `nrots=...`; it must not be inferred from another
project because separate refine3D runs may use different PFTC grids.

The first project-backed opt-in check incorrectly inferred `dang=1.25` from a
separate default-off project. Against the opt-in project this reported 1,113
invalid indices and 1,954 nearest-index mismatches among 3,081 active
particles. A second diagnostic run with an assumed 576-angle grid showed that
both projects actually contained indices from 1 through 440. The producing
grid was therefore confirmed as `nrots=440`, with `dang=0.818182` degrees.
These failures exposed validator configuration errors; they were not accepted
as production metadata failures. The validator now requires the producing
run's explicit `dang` or `nrots`, and restart no longer trusts an integer index
produced on a potentially different PFTC grid.

Project-backed metadata validation then passed on Oracle Linux. The opt-in
project contained 3,081 active particles, 1,224 continuous angles, and zero
invalid indices, non-finite values, or integer/angle mismatches. The matched
default-off project contained 3,081 active particles, zero continuous angles,
and zero invalid indices, non-finite values, or mismatches. Both tests ended
with `SIMPLE_TEST_CONT_INPL_REFINE3D_METADATA NORMAL STOP`. The remaining
Phase 3D-5 acceptance step was an opt-in restart followed by the same metadata
validation on its output.

The restart workflow completed normally and produced
`21_refine3D/betagal.simple`. Its metadata validation found 3,081 active
particles, 1,652 continuous angles, and zero invalid indices, non-finite
values, or integer/angle mismatches on the 440-angle grid. It ended with
`SIMPLE_TEST_CONT_INPL_REFINE3D_METADATA NORMAL STOP`. This verifies that a
project containing continuous angles can restart through the opt-in route and
write a consistent continuous representation again.

The opt-in restart workflow completed in 472.3 seconds. A matched default-off
control started from the same `18_refine3D/betagal.simple` input and completed
normally in 460.3 seconds. The 12.0-second difference is approximately 2.6%.
Together with the earlier matched observation in which opt-in was about 1.4%
faster, one pair in each direction does not establish either a speedup or a
regression. It does rule out the initially suspected approximately 60-second
feature penalty. Longer-run timing should be interpreted from per-iteration
measurements and repeated matched inputs.

A longer opt-in restart was requested with `maxits=6` from
`18_refine3D/betagal.simple`. It converged after iteration 2 and therefore ran
iteration 3 as the final combined iteration rather than consuming all six
allowed iterations. The workflow ended with `SIMPLE_REFINE3D NORMAL STOP` in
1,255.3 seconds. Orientation overlap increased from 0.987 to 0.997 and then
1.000. Mean in-plane update distance decreased from 0.085 to 0.020 and then
0.006 degrees, while mean shift increment decreased from 0.062 to 0.018 and
then 0.006. The FSC-area score increased from the starting 0.6330 to 0.6468,
then 0.6521, and remained stable at 0.6520 in the final iteration. Reported
resolution stabilized at 3.739 A.

The resulting `23_refine3D/betagal.simple` metadata also passed: 3,081 active
particles, 2,070 continuous angles, and zero invalid indices, non-finite
values, or integer/angle mismatches on the 440-angle grid. The only numerical
runtime warning was the previously observed array temporary in
`fit_straight_line`; rotational-storage particle I/O notices also appeared and
may affect wall time but are not continuous-refinement failures.

The independent fresh-source run requested with `maxits=6` also completed.
It converged after iteration 3 rather than consuming all six allowed
iterations and ended with `SIMPLE_REFINE3D NORMAL STOP` in 1,211.6 seconds.
Output `24_refine3D/betagal.simple` passed the metadata contract with 3,081
active particles, 2,399 continuous angles, and zero invalid indices,
non-finite values, or integer/angle mismatches on the 440-angle grid.

The workflow now has opt-in aggregate route diagnostics for production smoke
tests. Each distributed matcher partition reports the number of selected
particles that attempted continuous refinement and the mutually exclusive
terminal outcomes: improved, finite no-improvement, or legacy fallback. The
matcher verifies that the three outcome counts sum to the attempted count.
These counters are collected per OpenMP thread and reduced after matching, so
they do not add per-particle output or synchronization. Collection and output
are disabled by default and are enabled only with
`SIMPLE_INPL_DIAGNOSTIC=yes`; `inpl_cont=yes` must also be active. Ordinary
runs do not allocate or update the diagnostic counters.

## Validation gates

### Gate A: numerical unit tests

Status: implemented through the shared continuous-angle numerical tests and
the refine3D policy, state, joint-store, direct-route, metadata-state, recovery,
and score/shift consistency contracts. All seven mother-suite groups passed on
Oracle Linux before the final `eval`/`sigma` policy additions. The updated
policy suite and current mother executable require one final Oracle run.

Add focused tests for:

- route identity and objective consistency;
- finite-difference gradient agreement;
- periodic angle behavior near $0^\circ$ and $360^\circ$;
- bounds and fallback behavior;
- no-op behavior when `inpl_cont=no`.

### Gate B: 3D synthetic recovery

Status: complete and verified on Oracle Linux as
`simple_continuous_inplane_refine3D_recovery.f90`. The corrected focused case
and the seven-group mother suite passed. The mother test exposes it as
`case=synthetic_recovery` and requires `vol1`.

Use a controlled 3D fixture with known projection direction, in-plane angle,
and shifts. Compare:

```text
grid-only refine3D
continuous refine3D
```

The fixture first reproduces production's legacy preselection search to obtain
the native-frame `xy_first` shift seed. It then scores the grid with that seed
and starts joint refinement from the same pose. The continuous route must not
worsen the accepted objective and must reduce the angular error for a truth
between grid points. The recovered image-plane shift RMS error must remain below
0.25 pixels.

Run the focused fixture on Oracle Linux with:

```text
simple_test_continuous_inplane_refine3D \
  case=synthetic_recovery \
  vol1=/absolute/path/reference_volume.mrc
```

### Gate C: real-data refine3D default-off regression

Status: complete. Matched-source default-off workflows completed normally and
their project metadata remained grid-aligned with zero invalid, non-finite, or
integer/angle mismatch counts.

Run the existing biological `refine3D` workflow without `inpl_cont`. Confirm:

- normal completion;
- unchanged metadata contract;
- no continuous angle values are written;
- no CC or SGD path is entered.

### Gate D: real-data opt-in regression

Status: complete. The no-improvement score/shift inconsistency described in the
final audit is corrected and covered by both the focused direct-route contract
and a production diagnostic smoke. The final run accounted for all 3,081 active
particles with 2,830 improvements, 251 finite no-improvement outcomes, zero
fallbacks, and balanced totals. Its output metadata had 1,930 continuous angles
and zero invalid, non-finite, or mismatch counts.

Run the same `refine3D` workflow with `objfun=euclid inpl_cont=yes`.
Compare resolution, runtime, objective values, orientation stability, and
restart behavior with the default-off run.

### Gate E: rejection tests

Status: complete. Unsupported value, CC, hybrid, denoised, projection-
reconstruction, evaluation-only, sigma-only, and probabilistic combinations
were rejected by the focused policy suite without changing omitted/default-off
or eligible opt-in behavior.

Verify that unsupported combinations fail clearly, especially:

```text
objfun=cc inpl_cont=yes
refine=eval inpl_cont=yes
refine=sigma inpl_cont=yes
```

This gate must not alter the normal CC path when `inpl_cont` is omitted or set
to `no`.

## Implemented file ownership

The implementation follows these ownership boundaries:

- `simple_ui_refine3D` for the opt-in user key;
- `simple_parameters` and `simple_parameters_phases` for defaults and gates;
- `simple_commanders_refine3D` for policy propagation;
- `simple_refine3D_strategy` for lifecycle propagation and stripping
  matcher-only arguments from child workflows;
- `simple_strategy3D_srch` for candidate state and joint refinement;
- `simple_pftc_shsrch_grad` for the bounded joint API and result validation;
- the existing `simple_polarft_calc` and `simple_polarft_corr` Euclidean
  objective and derivative implementation, without a duplicate 3D evaluator;
- focused `production/tests` helper modules and one registered mother test for
  the 3D acceptance gates.

## Definition of done

The original eight feature-completion items are implemented and verified. The
seven-group mother suite, expanded policy suite, synthetic recovery fixture,
default-off and opt-in workflows, restart and multi-iteration checks, final
production smoke, and metadata validation passed on Oracle Linux. The P2/P3
hardening changes documented above also passed their targeted Oracle rerun.

The transfer is complete only when:

1. `inpl_cont=no` is the default and preserves the historical `refine3D`
   path.
2. `inpl_cont=yes` is explicitly accepted only for supported Euclidean paths.
3. Standalone `refine3D` does not invent an `abinitio2D`-style stage split;
   every default-off matcher iteration remains on the historical grid route.
4. Every eligible opt-in `refine3D` matcher iteration performs bounded joint
   $(s_x,s_y,\theta)$ refinement after discrete winner selection. This phase
   does not add $s_z$.
5. Integer and continuous orientation metadata remain consistent.
6. Restarted runs preserve the selected representation.
7. CC behavior is unchanged and unsupported opt-in combinations are rejected.
8. Focused unit, synthetic, default-off, opt-in, and real-data tests pass.
