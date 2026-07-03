# Refactoring Plan: Pairwise Logistic Class-Average Rejection Model

Date: 2026-07-03

## Summary

The current `model_cavgs_rejection` framework uses hard validity gates followed
by a learned monotone linear score over normalized class-average features. This
works well for broad quality signals, but it cannot naturally express
conditional evidence such as:

```text
low corr_frc_proxy is usually bad
unless support, center/edge signal, and population look like a valid top view
```

A pairwise logistic model would keep the existing feature extraction, hard
gates, training tables, and apply/analyze/evaluate command structure, but replace
the scalar weighted sum with a probability:

```text
P(accept | z_features, selected pairwise feature interactions)
```

This is a smaller and safer step than a full Bayesian network. It can express
important feature-feature dependencies while remaining fast, inspectable, and
trainable from the current analysis-table data.

## Current Baseline

Relevant files:

- `src/main/cavg_quality/simple_cavg_quality_types.f90`
- `src/main/cavg_quality/simple_cavg_quality_feats.f90`
- `src/main/cavg_quality/simple_cavg_quality_model.f90`
- `src/main/cavg_quality/simple_cavg_quality_learn.f90`
- `src/main/cavg_quality/simple_cavg_quality_analysis.f90`
- `doc/microchunk_and_rejection/model_cavgs_rejection.md`

The current learned artifact contains:

- a feature-family policy;
- non-negative feature weights;
- clustering/Otsu/min-acceptance threshold controls.

The model score is effectively:

```text
score = sum_i weight_i * z_feature_i
```

The classifier then chooses a dataset-relative boundary using score clustering,
optional Otsu thresholding, optional cluster rescue, and optional minimum
accepted fraction enforcement.

## Proposed Model

Add a new model family, tentatively named `pairwise_logistic`, with linear and
selected pairwise terms:

```text
logit(P(accept)) =
    beta0
  + sum_i beta_i * z_i
  + sum_(i,j) beta_ij * z_i * z_j
```

The predicted probability is:

```text
P(accept) = 1 / (1 + exp(-logit))
```

The most important interactions to consider first are not all possible pairs.
They should be a curated, small set motivated by observed failure modes:

- `corr_frc_proxy * log_center_edge_snr`
- `corr_frc_proxy * cc_area_frac`
- `corr_frc_proxy * log_pop`
- `corr_frc_proxy * neg_log_locvar_fg`
- `cc_area_frac * log_center_edge_snr`
- `log_locvar_fg * log_locvar_bg`
- `neg_log_locvar_fg * cc_area_frac`
- `presence * cc_area_frac`

This keeps the model compact and limits overfitting. With 12 base features, a
full pairwise expansion would add 66 interaction terms; that is probably too
wide for the current number of independent datasets.

## Model Artifact Changes

Extend `cavg_quality_model_spec` to represent two model families:

```text
model_family = linear_score | pairwise_logistic
```

For `pairwise_logistic`, the model file should include:

- `intercept`
- `feature_policy`
- `linear_coefficients`
- `interaction_terms`
- `interaction_coefficients`
- `prob_threshold`
- `regularization_lambda`
- optional `calibration_temperature`

Example model-file shape:

```text
model_version=9
name=chunk100mics_logistic
model_family=pairwise_logistic
feature_policy=microchunk_plus_score_signal
intercept=...
linear_coefficients=...
interaction_terms=corr_frc_proxy:log_center_edge_snr,corr_frc_proxy:cc_area_frac
interaction_coefficients=...
prob_threshold=...
regularization_lambda=...
```

The existing `model_version=8` linear model reader should remain supported for
backward compatibility until the old presets are intentionally retired.

## Training Changes

Training should move from AUC-derived non-negative weights to regularized
logistic regression over manually labelled trainable rows.

Recommended objective:

```text
minimize weighted_binary_cross_entropy + lambda * L2(coefficients)
```

Use dataset-aware weighting so large datasets do not dominate:

- give each dataset equal total weight;
- within each dataset, balance manual good and manual bad rows when both exist;
- keep the existing special handling for good-only or bad-only datasets in
  evaluation diagnostics, but avoid letting one-class datasets define
  interaction coefficients by themselves.

Regularization is not optional. The model should search a small grid of
`lambda` values and choose by leave-one-dataset-out or grouped cross-validation.

Suggested first grid:

```text
lambda = 0.001, 0.003, 0.01, 0.03, 0.10, 0.30, 1.0
prob_threshold = 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65
interaction_policy = none, score_conditionals, compact_conditionals, all_curated
```

## Existing Optimization Infrastructure

SIMPLE already has generic scalar-objective optimization infrastructure under
`src/main/opt`. The pairwise logistic implementation should try to reuse this
before introducing a new optimizer.

Relevant infrastructure:

- `simple_opt_spec.f90`: callback-based optimizer specification, including
  single- and double-precision cost, gradient, and combined cost/gradient
  interfaces.
- `simple_opt_factory.f90`: factory for constructing optimizer objects from
  `opt_spec`.
- `simple_opt_lbfgsb.f90`: bounded L-BFGS-B implementation.
- `simple_opt_bfgs2.f90`: unconstrained BFGS-style implementation.

Existing production-style examples include:

- `simple_pftc_shsrch_grad.f90`, which uses `lbfgsb` with double-precision
  cost, gradient, and combined `fdf` callbacks;
- `simple_opt_image_weights.f90`, which wraps an optimization state object and
  passes it through `class(*)` callback dispatch.

Recommended first optimizer path:

1. Build a compact logistic-training state object holding the design matrix,
   labels, dataset weights, regularization strength, and interaction metadata.
2. Implement a double-precision `fdf` callback that returns weighted binary
   cross entropy plus L2 penalty and its analytic gradient.
3. Use `opt_spec%specify('lbfgsb', ...)` with broad coefficient bounds, for
   example `[-20, 20]`, to prevent runaway coefficients on separable data.
4. Register `costfun_8`, `gcostfun_8`, and `fdfcostfun_8`, following the
   existing production callback pattern.
5. Use `bfgs2` only as a secondary option if bounded coefficients prove
   unnecessary and the unconstrained path is numerically stable.

The bounded `lbfgsb` route is a better first fit than a custom optimizer because
logistic training is smooth, differentiable, and low-dimensional. The required
new work is mostly the objective wrapper, dataset weighting, regularization, and
diagnostics, not the numerical minimizer itself.

The logistic `fdf` callback should:

- consume the existing standardized normalized features;
- compute logits and probabilities over all trainable rows;
- accumulate gradient for intercept, linear coefficients, and interaction
  coefficients;
- add L2 regularization to all coefficients except the intercept;
- clamp logits or use a stable sigmoid/cross-entropy formulation to avoid
  overflow;
- return both objective and gradient in one pass over the training rows.

No stochastic machinery is required for the first version. A small custom
batch-gradient optimizer should be kept as a fallback only if the existing
optimizer callback plumbing becomes a larger maintenance burden than expected.

## Apply-Time Changes

Feature extraction and hard gates remain unchanged.

For non-hard-rejected rows:

1. build the active base feature vector;
2. build the active pairwise terms;
3. compute logistic probability;
4. accept when `P(accept) >= prob_threshold`;
5. write probability as the quality score.

The classifier can still report:

- `quality_score` as probability;
- `quality_cluster` as a simple probability bin or `1/2` accept/reject label;
- `soft_decision='pairwise_logistic'`;
- `soft_reason='probability_threshold'`.

The existing k-medoids/Otsu/min-acceptance machinery should initially remain
available only for `linear_score` models. Mixing adaptive score clustering with
probability calibration would make the first version harder to interpret.

## Diagnostics

The learn and evaluate reports should add:

- selected model family;
- coefficient table sorted by absolute coefficient;
- interaction-term table with feature names and coefficients;
- per-dataset log loss and balanced accuracy;
- probability-threshold scan;
- calibration bins, for example predicted probability deciles versus empirical
  accept fraction;
- leave-one-dataset-out summary;
- interaction ablation summary.

Important diagnostic questions:

- Does the model improve MotAB without increasing false positives on PepT2,
  SLC, RYPER, or other difficult datasets?
- Are the largest coefficients stable under leave-one-dataset-out?
- Does the model rely on a single dataset-specific interaction?
- Does the probability threshold sit on a search-grid edge?

## Implementation Phases

### Phase 1: Model Representation

- Done: added `model_family` to `cavg_quality_model_spec` and
  `cavg_quality_model`.
- Done: added intercept, linear coefficients, interaction term indices,
  coefficients, probability threshold, and regularization metadata.
- Done: updated model read/write, analysis comments, and fixed-model summaries
  for the new family.
- Done: kept `linear_score` as the default behavior for existing models.

### Phase 2: Logistic Scoring

- Done: added logistic probability evaluation to
  `simple_cavg_quality_model.f90`.
- Done: routed `classify` by `model_family`.
- Done: hard-rejected rows remain final and get probability zero in the
  pairwise-logistic path.

### Phase 3: Training

- Done: added an explicit logistic candidate path in
  `simple_cavg_quality_learn.f90`, selected with
  `quality_mode=learn quality_model=logistic`.
- Done: added a logistic-training state object with design matrix, labels,
  weights, regularization strength, and interaction metadata.
- Done: implemented a double-precision `fdf` callback for weighted binary cross entropy
  plus L2 penalty.
- Done: uses existing `lbfgsb` infrastructure with broad coefficient bounds as the
  first optimizer path.
- Done: searches over feature policy, regularization, and probability threshold.
- Reuse existing analysis-table reading and dataset-role diagnostics.
- Pending: add leave-one-dataset-out diagnostics and interaction-policy
  ablations once enough training data accumulates.

### Phase 4: Reporting and Validation

- Add coefficient, interaction, threshold, calibration, and leave-one-dataset-out
  diagnostics.
- Compare against the current `chunk100mics` model and the `chunk_hard_reject`
  approximation on the same training file tables.
- Add a held-out or leave-one-dataset-out report section focused on known
  failure modes such as MotAB top views and PepT2 false positives.

### Phase 5: Promotion Policy

- Promote only if the pairwise model improves dataset-level generalization, not
  just aggregate class-level metrics.
- Keep the existing linear model selectable until enough production runs confirm
  the probabilistic model is stable.
- Do not replace the streaming `cluster2D_rejector` path unless the stream
  lifecycle, sentinels, stack writing, and particle-state propagation are
  explicitly integrated.

## Expected Advantages

- Captures conditional evidence without hardcoded rescue rules.
- Gives a direct probability-like score for reports and downstream thresholds.
- Can model cases where a feature is useful only in the presence of another
  feature.
- Should handle MotAB-style top views more naturally than a monotone linear
  score.
- Keeps apply-time cost negligible compared with image feature extraction.

## Risks

- The number of independent datasets is small relative to the number of possible
  interactions.
- A full pairwise expansion can overfit specimen-specific quirks.
- Probability calibration may look good in aggregate while failing on one
  specimen family.
- Interactions involving stored score may learn data-collection or refinement
  artifacts rather than intrinsic class-average quality.
- Reporting becomes more complex because the model is less visually obvious
  than a weighted feature sum.

## Computational Cost

Apply-time cost is small:

- base features already exist;
- curated pairwise terms are simple multiplications;
- one logistic function per trainable class is negligible.

Training cost is moderate:

- for hundreds to a few thousand class averages and a curated interaction set,
  `lbfgsb` training should be seconds to minutes;
- leave-one-dataset-out and grid search multiply that cost but remain offline;
- feature extraction remains the expensive part when working from projects, but
  training from analysis tables avoids re-extracting images.

The main cost is development and validation, not runtime.

## Development Effort

Estimated effort for a production-quality first version:

- model representation and I/O: 1-2 days;
- logistic scoring path: 1 day;
- logistic training wrapper using `lbfgsb`: 2-4 days;
- regularization, grouped weighting, and grid search: 2-4 days;
- diagnostics and reports: 2-4 days;
- validation on current and held-out data: 2-5 days.

Total: roughly 1.5-3 focused weeks. Existing optimizer infrastructure should
reduce numerical-optimization work, but report quality and dataset-level
validation remain the dominant effort.

## Recommended First Cut

Start with a constrained `pairwise_logistic` implementation:

- use all current normalized base features selected by feature policy;
- allow only curated interaction policies, not arbitrary all-pairs by default;
- train with L2 regularization and dataset-balanced weights using `lbfgsb`;
- use a fixed probability threshold selected by grouped validation;
- preserve the existing linear model path and hard-gate options.

The first success criterion should not be class-level aggregate accuracy. It
should be improved dataset-level behavior: rescue MotAB top views without
materially increasing false positives on PepT2, SLC, RYPER, or similar datasets.
