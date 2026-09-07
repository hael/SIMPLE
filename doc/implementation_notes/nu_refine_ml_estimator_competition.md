# ML-estimator competition in adaptive nonuniform refinement

## Status

Design proposal, 2026-09-06. No implementation is approved by this note.

This is the living design record for using the ML-regularized half-map pair in
`nu_refine=yes` volume-domain nonuniform filtering. It contains the proposed
contract, implementation shape, review findings, experiment order, and
validation criteria. If the design is implemented, update this note rather
than creating a parallel specification or plan.

The current production policy remains authoritative until an experiment is
accepted: the unregularized half pair builds and extends the adaptive NU bank,
and the ML-regularized auxiliary pair is not supplied when `nu_refine=yes`.
Static `nu_refine=no` behavior is outside this proposal and remains unchanged.

## 1. Context

Both reconstruction backends now use the same gridding-style, post-hoc NU
competition. The unregularized even/odd pair supplies the discrete low-pass
bank and the derived `_nu_filt` references. With `nu_refine=no`, an
ML-regularized pair can replace the finest static candidate when its effective
resolution is finer than that candidate. With `nu_refine=yes`, the adaptive
shell walk currently owns resolution extension and excludes that auxiliary
replacement.

The adaptive route therefore does not use a potentially useful second
estimator. Simply adding the complete ML pair as one terminal resolution label
would be ambiguous: ML regularization has a smooth spectral transfer rather
than a hard cutoff, and the resulting label would mix two separate questions:

1. what local resolution is supported by the data;
2. which estimator is preferable at that supported resolution.

The proposed design separates those decisions.

## 2. Decision

Represent the adaptive NU state as two fields:

```text
z(v) = (k(v), s(v))

k(v) = selected local cutoff / resolution label
s(v) = raw or ML estimator at that same cutoff
```

The unregularized half pair exclusively determines `k(v)`. After that field is
fixed, the unregularized and ML-regularized pairs compete at matched cutoffs to
determine `s(v)`.

This factorization is the recommended first experiment and prospective
production design. Allowing ML to advance `k(v)` is a separate, later
experiment with stronger gates; it must not be folded into the first
implementation.

## 3. Required scientific contract

### 3.1 Resolution ownership

- Build the static candidate bank from the current unregularized even/odd pair.
- Run the existing `nu_refine=yes` shell walk from that pair only.
- Preserve the current frontier definition, unary-win fraction, minimum
  absolute support, bank thinning, memory cap, ordered-label cleanup, and
  persisted high-resolution depth.
- Keep the unregularized FSC/cFAR as the resolution authority.
- Derive the NU local-resolution map and matching low-pass handoff only from
  `k(v)`. The estimator source must not change their meaning.

### 3.2 Matched-cutoff estimator competition

For every retained cutoff `k` selected somewhere in `k(v)`, construct two
candidate pairs:

```text
raw(k) = low-pass_k(unregularized even/odd)
ML(k)  = low-pass_k(ML-regularized even/odd)
```

Evaluate both with the same radially whitened, symmetric cross-half NU
objective used by the ordinary candidate bank. Compare the estimators only on
voxels whose fixed resolution field selected `k`, with the usual
candidate-scale normalized objective smoothing.

The ML estimator may be selected only when it improves on the raw estimator by
a declared nonzero margin and has sufficient spatial support. Raw wins ties.
The first experiment must not rely on a zero-margin winner-takes-all rule,
because adding a second candidate creates a multiple-comparison advantage and
can produce unstable source speckle without demonstrating a material gain.

The source field may receive a lightweight binary spatial regularizer or
hysteresis after the resolution field has been fixed. It must never feed back
into the ordered resolution-label optimization in this first experiment.

### 3.3 Half independence and synthesis

- Use an ML-regularized even/odd pair, not a merged ML map.
- A selected raw voxel in the even output comes from the raw even candidate;
  a selected ML voxel comes from the ML even candidate. Apply the analogous
  rule independently to the odd output.
- Never copy phase-bearing information from one half into the other.
- Continue to form the merged `_nu_filt` map only as the average of the final
  independently synthesized even/odd outputs.
- Keep the base and ML reconstruction products unchanged. `_nu_filt` remains
  a derived matching reference, not a replacement reconstruction artifact.

### 3.4 Masks and support

- Preserve spherical `mskdiam` support as the NU objective domain.
- With `automsk=yes`, preserve the current preliminary static-bank evidence
  envelope as the coarsest-label background constraint on `k(v)`.
- Do not use the density envelope, the NU-evidence envelope, or the complement
  of either as the ML source-selection domain.
- Outside the active NU support, preserve existing output behavior and record
  no ML source assignment.

## 4. Proposed processing sequence

For each state and iteration:

1. Obtain the unregularized and ML-regularized even/odd pairs after the same
   trailing blend and geometry preparation required by the shared NU path.
2. Build the static NU bank from the unregularized pair without an auxiliary
   replacement.
3. Generate and arm the NU-evidence background constraint when requested.
4. Optimize the static cutoff field.
5. Run the current raw-only high-resolution shell extension and ordered-label
   cleanup, producing the fixed `k(v)` field.
6. Construct matched-cutoff raw and ML candidates only for retained cutoffs
   that are populated in `k(v)`.
7. Evaluate their cross-half objectives and determine the guarded source field
   `s(v)` without changing `k(v)`.
8. Synthesize independent even/odd NU outputs from `(k(v), s(v))`.
9. Write the merged `_nu_filt` and the existing `_nu_locres` diagnostic.
10. Record the matching low-pass from `k(v)` and report source-selection
    diagnostics separately.

The common execution site remains `simple_nu_state_filter`; backend-specific
commanders and PCG operators must not acquire independent versions of this
policy.

## 5. Candidate representation and memory

The estimator source is not a second ordered resolution label. In particular,
`raw(k)` and `ML(k)` share the same cutoff coordinate. This avoids presenting
the ordered-label Potts model with two spectrally different candidates at an
invented distance from one another.

Do not double the persistent full candidate bank by retaining an ML sibling for
every label. The initial implementation should process populated cutoffs one at
a time, retain only the information required for final source selection, and
release temporary full-volume buffers promptly. The implementation plan must
account for:

- mask-packed raw-versus-ML score differences or source labels;
- temporary matched-cutoff filtered images;
- reuse or regeneration of existing raw cache products;
- cleanup on normal and error-prone paths;
- the existing NU candidate-bank memory cap.

If recomputing ML candidates is too expensive, a measured cache design may be
proposed later. Disk or memory caching is an implementation choice, not part of
the scientific contract.

## 6. Diagnostics

The first implementation should report enough information to distinguish a
useful estimator choice from candidate-count bias:

- raw and ML selected voxel counts and percentages by cutoff;
- score-improvement distribution for ML-selected voxels;
- number of assignments rejected by the margin or support gate;
- spatial fragmentation or boundary measure of the source field;
- iteration-to-iteration ML occupancy by state;
- selected-cutoff histogram and high-resolution depth, unchanged in meaning;
- matching low-pass trajectory;
- half-map FSC and beyond-band spectral excess.

An optional development-only source map may encode raw versus ML selection.
It must not be called a local-resolution map or become a matcher input.

## 7. Experiment order and gates

### Arm A: raw-only control

The current `nu_refine=yes` path. This establishes the resolution trajectory,
matching bandwidth, map quality, runtime, and memory baseline.

### Arm B: fixed-resolution estimator competition

Implement Sections 2--6. This is the recommended experiment.

Arm B passes only if:

- synthetic truth metrics improve or remain neutral overall and do not regress
  materially in any established difficult region;
- real-data half-map FSC is not spuriously inflated;
- matching stability and resolution progression are equal or better than Arm A;
- ML selection is spatially coherent and associated with a measurable
  cross-half objective improvement rather than near-tie churn;
- beyond-band excess and abrupt reference-amplitude changes do not worsen;
- the source choice remains stable enough across adjacent iterations to avoid
  a new reference fixed-point oscillation;
- peak memory and runtime remain acceptable for both reconstruction backends.

Failure should first be adjudicated by increasing the selection margin or
support/hysteresis guard. Do not allow the ML candidate to alter the resolution
field as a workaround for a failing fixed-resolution competition.

### Arm C: optional ML-assisted resolution extension

Attempt only after Arm B passes and only if the results suggest that useful ML
detail is systematically truncated by the raw-owned endpoint.

After the raw shell walk stops, evaluate raw and ML challengers at the same
single next shell. An ML-assisted shell may be accepted only if:

- the ML challenger clears the existing frontier percentage and absolute seed
  requirements;
- it beats the raw challenger by a calibrated margin;
- the raw challenger is not strongly rejected at those voxels;
- the first experiment is capped at one shell beyond the raw endpoint;
- support persists on the following iteration.

Mark such a shell and its matching-band contribution explicitly as
ML-assisted. Do not silently merge it with raw-supported shell history.
Candidate-count calibration is mandatory because the minimum of raw and ML
challenger costs has a built-in selection advantage over a single incumbent.

Arm C must beat both Arms A and B on synthetic truth and truth-free real-data
stability before it can be considered for production. Otherwise retain Arm B
and keep raw data as the sole resolution-extension authority.

## 8. Validation matrix

Use at least:

- the existing synthetic truth-judged reconstruction fixture;
- PfCRT, because matching-band restrictions have previously caused a
  high-resolution refinement regression there;
- bgal or the current representative PCG real-data case;
- both gridding and PCG reconstruction backends through the shared NU path;
- `automsk=no` and `automsk=yes`;
- restart/continuation with persisted NU high-resolution depth;
- plain `nonuniform` and, where applicable, `nonuniform_lpset` reference
  topology.

Compare Arms A and B at identical inputs, iteration count, matching settings,
and reconstruction convergence. If Arm C is attempted, include both preceding
arms in the same adjudication; do not compare it only with raw-only control.

## 9. Ownership map for a future implementation

- `simple_nu_state_filter`: shared orchestration, input lifecycle, experiment
  sequencing, output writing, and matching-LP handoff.
- `simple_nu_filter`: estimator-competition state and public lifecycle.
- NU filter submodules: matched-cutoff candidate generation, objective
  comparison, guarded source selection, synthesis, diagnostics, and cleanup.
- reconstruction strategies/commanders: provide the unregularized and
  ML-regularized half pairs to the shared path; no scientific selection logic.
- matcher reference utilities: unchanged consumption of `_nu_filt` products.
- policy documents: update only after an experiment is accepted.

Shared-memory and distributed execution may differ in transport and launch,
but must use the same scientific competition and artifact semantics.

## 10. Non-goals

- Reintroducing the removed in-solve NU prior or any NU precision in PCG.
- Replacing the global ML `P_tau` regularized reconstruction.
- Letting the ML pair redefine FSC/cFAR or unregularized resolution authority.
- Feeding a merged map into either half-map candidate path.
- Multiplying matching references by density or NU-evidence envelopes.
- Changing static `nu_refine=no` auxiliary replacement behavior.
- Making the source map a new public tuning surface before the experiment has
  established a robust automatic rule.

## 11. Open decisions before implementation

1. Define and calibrate the raw-versus-ML score margin.
2. Choose the minimum spatial support and whether a binary Potts step,
   connected-support rule, or temporal hysteresis is the simplest stable
   source regularizer.
3. Decide whether matched-cutoff ML candidates are generated only for populated
   labels or also for a narrow neighborhood around each selected label.
4. Define the development-only source-map encoding and restart diagnostics.
5. Establish runtime and peak-memory budgets for Arm B before selecting a cache
   strategy.
6. Define quantitative pass/fail tolerances for the validation metrics before
   examining Arm B results.

