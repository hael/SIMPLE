# model_cavgs_rejection

## Scope

`model_cavgs_rejection` is the SIMPLE class-average rejection-model program. It evaluates `cls2D` class averages, applies hard validity rejects, extracts a fixed scalar feature bank, normalizes those features within the dataset, and applies a named model to partition class averages into accepted and rejected sets. The model can also abstain from additional soft rejection and report that the final decision was hard preselection only.

This is a learned feature-vector model. It should be read alongside the streaming microchunk rejector, which is a cumulative rule engine. The two systems share several image-processing primitives, but they use them differently: microchunking applies fixed sequential criteria and rejects a class as soon as any criterion fires; `model_cavgs_rejection` uses hard rejects only for validity failures, then learns a weighted scalar quality score and dataset-specific thresholding behavior for the remaining class averages.

The implementation lives in `src/main/cavg_quality`. The command entry point is `exec_model_cavgs_rejection` in `src/main/commanders/simple/simple_commanders_cavgs.f90`.

## Implementation Files

- `simple_cavg_quality_types.f90`: shared constants and derived types.
- `simple_cavg_quality_feats.f90`: feature definitions, feature extraction, hard rejects, and robust normalization.
- `simple_cavg_quality_stats.f90`: binary metrics, AUC, distance-matrix normalization, and statistical helpers.
- `simple_cavg_quality_model.f90`: built-in presets, model file I/O, scoring, clustering, thresholding, and promotion snippets.
- `simple_cavg_quality_analysis.f90`: evaluation, analysis reports, and feature tables.
- `simple_cavg_quality_learn.f90`: analysis-table reader, feature-policy search, model search, and learn reports.

## Microchunk Comparison

The stream path in `simple_microchunked2D` currently calls `cluster2D_rejector`. That engine is deliberately deterministic and rule based. `model_cavgs_rejection` is a separate backend that uses a learned feature vector over class averages.

| Aspect | Microchunk rejection | `model_cavgs_rejection` |
| --- | --- | --- |
| Decision type | Cumulative rule engine. | Learned feature-vector model with hard validity gates. |
| Primary owner | `simple_cluster2D_rejector` called from `simple_microchunked2D`. | `src/main/cavg_quality` called by `model_cavgs_rejection`. |
| Unit of evidence | Individual scalar rules applied one after another. | Normalized feature vector plus model weights, clustering, and score thresholding. |
| Rejection semantics | A class rejected by any rule remains rejected. | Hard rejects are final; remaining classes are scored and partitioned by the model. |
| Threshold source | Fixed constants, with tier-specific population and local-variance overrides. | Built-in or learned model specification; learn mode chooses feature policy, weights, and threshold controls from analysis datasets. |
| Adaptation to dataset | Robust local-variance z-scores are dataset-relative, but thresholds are fixed. | Robust feature normalization, k-medoids distances, Otsu thresholds, and selected model controls are dataset-relative. |
| Training data | None. | Manual selections from `quality_mode=analyze` files. |
| Output role | Online stream cleanup and particle deselection. | Batch/chunk/pool class-average selection, analysis, training, and model promotion. |

## Evidence Comparison

The model feature bank keeps the microchunk-style image-processing evidence, optional stored-score and signal evidence, and the local-variance/support evidence needed by optional overfit hard rejection. The overfit local-variance features are part of the learned model feature vector and are included in every learn-mode feature policy. Texture descriptors were removed after they failed to justify the extra feature and extraction complexity. The current `chunk100mics` preset uses the `microchunk_plus_signal` policy: microchunk plus overfit-family features, center/edge signal, presence, and fuzzy-ball signal evidence.

| Evidence | Microchunk rule engine | `model_cavgs_rejection` feature-vector model | Current chunk role |
| --- | --- | --- | --- |
| Class population | Hard rule: reject below a tier-specific fraction of total population. | `log_pop` feature plus hard reject below `0.0035` of total population. | Active feature and hard gate. |
| Resolution | Hard rule: reject when `res > 40.0`. | `neg_log_res` feature plus hard reject when `res > 40.0`. | Active feature and hard gate. |
| Foreground centering | Otsu foreground connected components; reject if any centroid lies outside the mask radius. | `centered` feature is negative normalized centroid displacement; centroid outside the mask is a hard reject. | Active feature and hard gate. |
| Foreground outside mask | Reject when the largest valid foreground component has more than 10 pixels outside the mask disc. | `cc_area_frac` measures largest component area relative to the mask disc; outside-mask excess is a hard reject. | Active feature and hard gate. |
| Full-image component cleanup | Connected components spanning the full image are removed before mask tests. | Same concept is applied during feature extraction before geometry features and hard mask rejection. | Hard-gate preprocessing. |
| Foreground local variance | Robust-z rule on foreground/background local variance after 10 A low-pass and Otsu masking. | `log_locvar_fg` feature; only both variances non-positive is a hard reject. | Active feature. |
| Background local variance | Robust-z rule paired with foreground local variance. | `log_locvar_bg` feature; only both variances non-positive is a hard reject. | Active feature. |
| Stored score evidence | Not used by the rule engine. | `corr_frc_proxy` reads stored class correlation or FRC-like score when present. | Inactive in the current chunk preset. |
| Center/edge signal | Not used by the rule engine. | `log_center_edge_snr` feature in the `signal` family. | Active feature. |
| Presence | Not used by the rule engine. | `presence` feature in the `signal` family. | Active feature. |
| Overfit local variance | Not used by the rule engine. | `neg_log_locvar_fg`, `neg_log_locvar_bg`, and `log_locvar_fg_bg_ratio` provide low/localized variance evidence. | Active features in learned policies. |
| Overfit band-pass localization | Not used by the rule engine. | `log_bp40_100_center_edge_var` measures center/edge variance after 100 to 40 A band-pass; raw values below `log(2)` match Joe's fuzzy-ball heuristic. | Active feature in learned policies. |
| Fuzzy-ball signal | Not used by the rule engine. | `fuzzy_ball_signal` combines foreground texture, central presence, and 100 to 40 A center/edge localization. Low values mark the fuzzy-ball pattern as continuous learned evidence, not as an apply-time hard gate. | Active feature in learned policies. |
| Invalid pixels | Not a named microchunk rule. | Image values containing invalid pixels are hard rejected. | Hard gate. |

## Hard-Reject Comparison

Hard rejects in `model_cavgs_rejection` are intended to mirror the non-negotiable validity parts of microchunking while leaving ordinary quality variation to the learned score. The main difference is local variance: microchunking uses local variance as a direct rejection rule, while `model_cavgs_rejection` keeps local variance as learned scalar evidence except for degenerate zero-variance cases.

| Criterion | Microchunk rule engine | `model_cavgs_rejection` |
| --- | --- | --- |
| Population | Rejects `pop < ceiling(sum(pop) * fraction)`. Fractions are tier-specific. | Rejects `pop <= 0` and `pop < ceiling(sum(pop) * 0.0035)`. |
| Resolution | Rejects `res > 40.0`. | Rejects `res > 40.0`. |
| Empty or mismatched class-average stack | Stream chunk is marked `REJECTION_FAILED` and `COMPLETE`. | Apply/analyze require a valid project, `cls2D` entries, matching class-average stack, and `mskdiam`; invalid inputs stop the command. |
| Invalid pixels | No separate named rule. | Hard rejected. |
| No valid foreground component | Rejected after full-image connected components are pruned. | Hard rejected after the same style of foreground-component pruning. |
| Foreground centroid outside mask | Rejected. | Hard rejected. |
| Largest foreground component outside mask | Rejected when outside pixels exceed 10. | Hard rejected when outside pixels exceed 10. |
| Local variance exactly degenerate | Rejected when both inside/outside scores are near zero. | Hard rejected when both foreground/background local variances are non-positive. |
| Local variance low but not degenerate | Rejected by fixed robust-z thresholds. | Not a hard reject; encoded as `log_locvar_fg` and `log_locvar_bg` features. |

Microchunk tier thresholds:

| Tier | Population fraction | Local-variance strong threshold | Local-variance weak threshold |
| --- | ---: | ---: | ---: |
| Pass 1 | `0.0050` | `-0.5` | `-0.1` |
| Pass 2 | `0.0035` | `-1.0` | `-1.0` |
| Reference chunk | `0.0025` | `-2.0` | `-2.0` |
| Match chunk | `0.0025` | `-2.0` | `-2.0` |

`model_cavgs_rejection` uses a fixed population hard-reject fraction of `0.0035` for the current feature extractor, independent of the selected model preset.

## Learning Model

The rejection model is a low-dimensional, monotone feature-vector classifier with explicit hard constraints. In learning-theory terms, the hard rejects are prior validity constraints: they remove classes the model is not asked to rescue or fit. The remaining rows form a supervised training set from manual selections.

The standard linear learned artifact has three parts:

- a feature policy, which selects a cumulative family set;
- non-negative feature weights, which define a linear scalar quality score;
- thresholding controls, which govern how the per-dataset score boundary is chosen after k-medoids clustering and optional Otsu thresholding.

`model_family=logistic` trains a pairwise logistic artifact instead. Its fitted parameters are an intercept, linear feature coefficients, pairwise feature-interaction coefficients, a probability threshold, and a regularization strength. The coefficients are fit from the training analysis files only. On two-class trainable datasets, the logistic loss gives manually selected classes moderately higher total weight than manually deselected classes so the probability surface learns rejection evidence without becoming an overly aggressive rejector. Manually deselected examples with the overfit/support-poor or band-pass-localization signature get an additional training-time loss multiplier, so fuzzy-ball examples influence the fit without becoming an apply-time hard reject. Learn reports also expose an internal overfit-focus loss scale for experiments; it is currently neutral.

The apply-time classifier is partly dataset-adaptive. Features are robustly normalized inside the current dataset. K-medoids and Otsu operate on the current score/feature distribution. The learned parameters provide the inductive bias: which features are active, how the score is weighted, and how aggressively the cluster-derived boundary is shifted or replaced.

The learning objective is empirical risk minimization over manually annotated `quality_mode=analyze` runs. Feature-weight candidates are learned only from datasets where both manual states remain present after hard rejects, because those are the datasets with a learnable soft boundary. Candidate models are scored over all scoreable datasets using only non-hard-rejected rows. Two-class trainable datasets contribute specificity with a small false-negative-rate tolerance over the trainable manually good rows. Datasets where at least 20% of the trainable class averages are manually bad rows matching the fuzzy-ball signature get an additional hinge/quadratic objective penalty when the accepted fuzzy-ball-signature rate is above the tolerated rate. The signature is poor support plus either low local-variance evidence or low 100 to 40 A band-pass center/edge localization. This keeps selected-class protection global while making the objective care specifically about fuzzy-ball leakage in small-specimen chunks. The final macro score is a fixed robust score: half the mean dataset score and half the mean score over the lower tail of datasets. During pairwise-logistic fitting, manually bad rows with the same fuzzy-ball signature get extra loss weight. Datasets where only manually good rows remain after hard rejects contribute recall, so they teach the model not to reject good classes that passed the hard gates. Datasets where only manually bad rows remain contribute specificity only when no manually good classes were removed by hard rejects. If manually good classes were entirely removed by hard rejects, the dataset is reported as hard-gate blocked and skipped for soft-model scoring. Hard-rejected rows remain visible in diagnostics, including counts of manually good classes lost to hard gates, but they do not participate in fitting the learned boundary.

Feature-weight search is ab initio for the supplied training data. The learner first constructs an AUC-derived seed:

```text
weight(feature) = max(0, pooled_auc(feature, manual_state) - 0.5)
```

The weights are normalized after the active feature policy is applied. Each feature policy contributes one deterministic AUC-derived weight vector to the full classifier search. There is no base-weight blending in the current learner. This means each training run derives the scalar model from the current analysis files rather than nudging an inherited default.

Feature-family policies act as a small structural model-selection layer. For the linear learner, the training data derive the AUC feature weights and select the feature policy and thresholding controls. For the logistic learner, the training data fit the coefficients and select the feature policy, regularization lambda, and probability threshold.

Use `filetab=` to define the target workflow scenario being trained. For example, early streaming chunk models should be fit from chunk analyses, while late pooled-refinement models should be fit from pool analyses. To inspect cross-scenario behavior without influencing the model, train the target model and then run `quality_mode=evaluate filetab=... infile=...` on the other scenario.

| Policy | Active families | Purpose |
| --- | --- | --- |
| `microchunk` | `microchunk`, `overfit` | Tests whether scalar analogues of the stream rejector plus overfit local-variance evidence are sufficient. |
| `microchunk_plus_score` | `microchunk`, `overfit`, `score` | Adds stored class correlation/FRC-like evidence. |
| `microchunk_plus_signal` | `microchunk`, `overfit`, `signal` | Adds center/edge and presence evidence without stored score evidence. |
| `microchunk_plus_score_signal` | `microchunk`, `overfit`, `score`, `signal` | Uses the full compact quality feature set. |

Learn mode reports feature signal, feature-drop diagnostics, and leave-one-dataset-out feature-policy diagnostics so that the selected model can be interpreted as a compact scientific statement about which evidence family and weight structure generalized on the supplied training set.

## Command Modes

`quality_mode=apply` computes quality features, applies the selected model, writes `cavgs_quality_features.txt`, writes `quality_selected_cavgs.mrc` and `quality_rejected_cavgs.mrc`, maps the selected class averages into particle states, annotates `cls2D` with `quality`, `quality_cluster`, and `accept`, annotates `cls3D` when its class count matches `cls2D`, optionally prunes particles with `prune=yes`, and writes the project.

`quality_mode=analyze` computes the same model output but treats the existing `cls2D` state as the manual reference. It writes `cavgs_quality_analysis.txt` and the selected/rejected stacks. The project selection is left unchanged.

`quality_mode=learn` reads a training file table of `cavgs_quality_analysis.txt` files from `filetab=` and searches for a model specification from a neutral `abinitio_learn_base` foundation. `model_family=linear|logistic` selects the model family to train; if omitted, learn mode uses `logistic`. It does not accept `quality_model` or `infile` as a seed. It writes a learned model file controlled by `fname=` and writes `cavgs_quality_learn_report.txt`.

`quality_mode=evaluate` applies the selected fixed model without refitting. With `filetab=`, it evaluates one or more saved `cavgs_quality_analysis.txt` files. Without `filetab=`, it evaluates a single project directly using the existing `cls2D` state as the manual reference, like analyze mode. It writes `cavgs_quality_evaluate_report.txt`, or the report path controlled by `fname=`.

`quality_mode=promote` reads a model file from `infile=` and writes a Fortran promotion snippet controlled by `fname=`.

`apply` and `analyze` require `projfile` and `mskdiam`. `learn` requires `filetab`. `evaluate` requires either `filetab` or `projfile` plus `mskdiam`. `promote` requires `infile`. The commander sets `oritype=cls2D`, defaults `mkdir=yes`, and defaults `prune=no`.

For `apply`, `analyze`, and `evaluate`, the command uses `chunk100mics` unless `quality_model` or `infile` is supplied. Use `overfit_hard_reject=yes` when a fixed support/local-variance hard gate should be layered on top of the standard hard gates. Use `chunk_hard_reject=yes` when the fixed legacy chunk-quality hard gate should be layered on top of the standard hard gates without loading or applying a model.

## Model Selection

`quality_model` selects a built-in preset outside learn mode. The promoted built-ins are:

- `chunk100mics`: default chunk/stream-style pairwise logistic model trained from `/Users/elmlundho/cavgs_quality/chunk100mic_training_data_v3`.
- `chunk100mics_linear`: interpretable linear chunk/stream-style model trained from `/Users/elmlundho/cavgs_quality/chunk100mic_training_data`.
- `pool`: late pooled-refinement pairwise logistic model trained from `/Users/elmlundho/cavgs_quality/pool_training`.

When `infile` is supplied, the model file is treated as a complete model and wins over the built-in preset.

Linear model files use `model_version=8`, pairwise logistic model files use `model_version=9`, and both use explicit key-value fields:

- `name`
- `feature_policy`
- `feature_weights`
- `boundary_margin`
- `min_score_separation`
- `otsu_min_offset`
- `otsu_max_offset`
- `cluster_rescue_margin`
- `min_accept_frac`
- `use_lowsep_otsu`
- `use_otsu_window`
- `use_cluster_rescue`
- `enforce_min_accept_frac`

Unknown model-file keys are rejected.

## Feature Bank

`CAVG_QUALITY_NFEATS` is 14. All feature directions are `higher_is_better`. The model consumes normalized `z_*` values; raw values are written for diagnostics.

| Index | Feature | Family | Source | Meaning |
| --- | --- | --- | --- | --- |
| 1 | `log_pop` | `microchunk` | `cls2D` `pop` | Log class population. |
| 2 | `neg_log_res` | `microchunk` | `cls2D` `res` | Negative log class resolution estimate. |
| 3 | `centered` | `microchunk` | Otsu foreground connected components | Negative normalized centroid displacement of segmented signal. |
| 4 | `log_locvar_fg` | `microchunk` | Foreground-mask local variance | Log local variance in the foreground Otsu mask. |
| 5 | `log_locvar_bg` | `microchunk` | Background-mask local variance | Log local variance outside the foreground Otsu mask. |
| 6 | `corr_frc_proxy` | `score` | `cls2D` `corr` when present | Stored class correlation or FRC-like quality proxy. |
| 7 | `log_center_edge_snr` | `signal` | Center/edge variance statistics | Log central variance relative to edge variance. |
| 8 | `cc_area_frac` | `microchunk` | Foreground connected-component area | Largest foreground component area divided by expected mask area. |
| 9 | `presence` | `signal` | Center/border statistics | Central signal excess relative to border background noise. |
| 10 | `neg_log_locvar_fg` | `overfit` | Foreground-mask local variance | Negative log local variance in the foreground Otsu mask; higher values indicate lower foreground variation. |
| 11 | `neg_log_locvar_bg` | `overfit` | Background-mask local variance | Negative log local variance outside the foreground Otsu mask; higher values indicate quieter background. |
| 12 | `log_locvar_fg_bg_ratio` | `overfit` | Foreground/background local variance | Log foreground/background local variance ratio; higher values indicate support-localized detail. |
| 13 | `log_bp40_100_center_edge_var` | `overfit` | Band-pass center/edge variance | Log center/edge variance ratio after 100 to 40 A band-pass; low raw values flag overfit fuzzy balls. |
| 14 | `fuzzy_ball_signal` | `overfit` | Combined overfit evidence | Sum of `log_locvar_fg`, `presence`, and `log_bp40_100_center_edge_var`; low values indicate fuzzy-ball overfit evidence. |

Feature families are used by learn mode. Every learned policy includes `microchunk` and `overfit`; the policy variants append optional evidence families:

- `microchunk`
- `microchunk_plus_score`
- `microchunk_plus_signal`
- `microchunk_plus_score_signal`

Features outside the selected policy are encoded by zero weights.

## Extraction Constants

- `FOREGROUND_SEG_LP = 30.0`
- `SIGNAL_METRIC_LP = 10.0`
- `OVERFIT_SIGNAL_BP_HP = 100.0`
- `OVERFIT_SIGNAL_BP_LP = 40.0`
- `CAVG_RES_HARD_REJECT_A = 40.0`
- `POP_FRACTION_HARD_REJECT = 0.0035`
- `LOCVAR_WINDOW = 10`
- `MASK_HARD_OUTSIDE_PIXELS = 10`
- `LOG_EPS = 1.0e-12`
- `CLIP_Z = 4.0`

## Hard Rejects

Hard rejects are validity gates applied before model fitting, clustering, and training-score calculation.

A class average is hard rejected when any of these conditions hold:

- `pop <= 0`;
- `pop < ceiling(sum(pop) * 0.0035)`;
- `res > 40.0`;
- the image contains invalid pixels;
- foreground segmentation has no valid component after full-image connected components are pruned;
- a segmented foreground component centroid lies outside the mask radius;
- the largest foreground component has more than 10 pixels outside the mask disc;
- foreground and background local variance are both non-positive.

Hard-rejected classes receive rejected state directly. Their normalized features are set to `-CLIP_Z`, their model scores are set to `-CLIP_Z`, and they remain visible in analysis and learning reports.

### Optional Overfit Hard Gate

`overfit_hard_reject=yes` is a no-model application path for project-backed `apply`, `analyze`, and `evaluate` runs. It applies the standard hard gates first, computes normalized features over the remaining rows, and then applies a fixed overfit hard rule. It does not read `infile`, does not call `model%classify`, and does not use k-medoids, Otsu thresholds, rescue, or minimum-acceptance model logic.

The default hard reject is:

```text
(z_neg_log_locvar_fg > 0.0 AND z_cc_area_frac < 0.5)
```

On the FlipQR overfit-training table, the support/local-variance rule rejects 18/26 trainable overfit-labeled classes while retaining 31/33 trainable good classes. Rows rejected by the overfit hard rule are marked as hard rejects in the final decision, but their normalized features remain inspectable in the output tables because those values explain the hard-gate decision.

### Optional Chunk Hard Gate

`chunk_hard_reject=yes` is a no-model application path for project-backed `apply`, `analyze`, and `evaluate` runs. It applies the standard hard gates first, computes normalized features over the remaining rows, and then applies a fixed legacy chunk-quality hard rule. It does not read `infile`, does not call `model%classify`, and does not use k-medoids, Otsu thresholds, model cluster-rescue, or minimum-acceptance model logic. It is mutually exclusive with `overfit_hard_reject=yes`.

The fixed chunk hard rule rejects a trainable class when any of these dataset-normalized gates fire:

```text
z_corr_frc_proxy      < -0.224682
z_neg_log_locvar_fg   < -2.073230
z_cc_area_frac        < -2.398173
z_log_center_edge_snr < -3.274551
z_log_pop             < -3.513258
```

The stored-score gate has one guarded exception for MotAB-style top views: if `z_corr_frc_proxy` is the only failing gate, the class is rescued when independent image-derived evidence is strong enough. The rescue conditions are:

```text
z_corr_frc_proxy > -1.0
AND z_log_center_edge_snr > 0.75
AND z_cc_area_frac > -2.4
AND z_log_pop > 0.45
```

or:

```text
z_corr_frc_proxy > -1.5
AND z_neg_log_locvar_fg > 3.0
AND z_cc_area_frac > 0.5
AND z_log_pop > -1.0
AND z_log_center_edge_snr > -1.5
```

These thresholds were derived from `/Users/elmlundho/cavgs_quality/chunk100mic_training_data` as a compact hard-gate approximation to the earlier learned chunk operating point. They have not been re-fit to the current `chunk100mics` v3 logistic preset. Rows rejected by the chunk hard rule are marked as hard rejects in the final decision, while their normalized feature values remain inspectable in the output tables.

## Normalization

`normalize_cavg_quality_features` computes per-feature robust statistics over non-hard-rejected rows only. Each feature is centered by the median and scaled by the Gaussian MAD. Values are clipped to `[-CLIP_Z, CLIP_Z]`. If a feature has no robust spread, its normalized value is set to zero for trainable rows.

## Built-In Presets

`chunk100mics` uses feature policy `microchunk_plus_signal` and the pairwise logistic family, with low-separation Otsu, Otsu windowing, cluster rescue, and minimum accepted fraction enforcement disabled. It was trained from `/Users/elmlundho/cavgs_quality/chunk100mic_training_data_v3` with the chunk learn objective in effect at the time of promotion.

```text
model_family        pairwise_logistic
prob_threshold      4.500000E-01
regularization      1.000000E-03
feature_weights     uniform over the active microchunk_plus_signal vector with zero weight on stored score evidence
```

```text
boundary_margin         0.00
min_score_separation    0.05
otsu_min_offset         0.00
otsu_max_offset         0.00
cluster_rescue_margin   0.20
min_accept_frac         0.00
use_lowsep_otsu         false
use_otsu_window         false
use_cluster_rescue      false
enforce_min_accept_frac false
```

`pool` uses feature policy `microchunk_plus_score_signal` and the pairwise logistic family, but is tuned for late pooled-refinement data from `/Users/elmlundho/cavgs_quality/pool_training`.

```text
model_family        pairwise_logistic
prob_threshold      3.000000E-01
regularization      1.000000E-04
feature_weights     uniform over the legacy score/signal vector with zero weight on newly added features
```

`chunk100mics_linear` preserves the previous linear score-and-threshold model as an interpretability tool. Its feature weights can be read directly as non-negative contributions to the normalized scalar quality score.

```text
log_pop             9.978756E-02
neg_log_res         1.167914E-01
centered            3.642511E-02
log_locvar_fg       1.329548E-01
log_locvar_bg       1.402481E-01
corr_frc_proxy      1.610645E-01
log_center_edge_snr 6.630784E-02
cc_area_frac        6.981037E-02
presence            1.294257E-01
neg_log_locvar_fg   0.000000E+00
neg_log_locvar_bg   0.000000E+00
log_locvar_fg_bg_ratio 4.718454E-02
log_bp40_100_center_edge_var 0.000000E+00
fuzzy_ball_signal   0.000000E+00
```

```text
boundary_margin         0.30
min_score_separation    0.05
otsu_min_offset         0.05
otsu_max_offset         0.40
cluster_rescue_margin   0.20
min_accept_frac         0.60
use_lowsep_otsu         false
use_otsu_window         true
use_cluster_rescue      false
enforce_min_accept_frac true
```

Linear model weights are clipped to non-negative values and normalized to sum to one when a linear model is initialized or read.

## Classification

For non-hard-rejected class averages, `model_family=linear_score` uses the historical score-and-cluster path:

1. Compute a linear quality score with `matmul(normalized_features, weights)`.
2. Build a pairwise Euclidean distance matrix in normalized feature space using nonzero-weight features.
3. Robustly normalize the distance matrix.
4. Run two-cluster k-medoids clustering.
5. Choose the good cluster as the cluster with the higher mean model score. If the score means tie, choose the larger cluster; if cluster sizes also tie, choose the cluster with the higher medoid score.
6. Place the raw threshold halfway between the good-cluster and bad-cluster mean scores.
7. If cluster separation is below `min_score_separation`, use the Otsu score threshold only when `use_lowsep_otsu` is enabled and Otsu separation is sufficient. Otherwise, accept all trainable class averages as a single cluster.
8. If cluster separation is sufficient, start from `raw_threshold - boundary_margin`.
9. If `use_otsu_window` is enabled, replace the threshold with Otsu when the Otsu threshold is at least `otsu_min_offset` and at most `otsu_max_offset` above the raw threshold, and Otsu separation is sufficient.
10. If `use_cluster_rescue` is enabled, also accept good-cluster members within `cluster_rescue_margin` below the threshold.
11. If `enforce_min_accept_frac` is enabled, accept enough top-scoring trainable rows to satisfy `min_accept_frac`.

For `model_family=pairwise_logistic`, the model computes:

```text
prob_accept = sigmoid(intercept + sum(linear_coefficients * features) + sum(interaction_coefficients * feature_i * feature_j))
```

Non-hard-rejected rows are accepted when `prob_accept >= prob_threshold`. Standard hard gates still run before the logistic model.

If there are fewer than four trainable rows, the distance matrix is degenerate, clustering fails, or the low-separation branch has no acceptable Otsu threshold, all non-hard-rejected rows are accepted as a single cluster.

`threshold_offset` is reported as `raw_threshold - threshold`. A positive value means the effective threshold is lower than the cluster midpoint.

The classifier reports its soft decision explicitly:

- `soft_decision=soft_threshold`: a learned score threshold was used to reject at least one non-hard-rejected class average.
- `soft_decision=hard_only`: no additional soft rejection was made after hard preselection.

Common hard-only reasons are `too_few_trainable`, `flat_feature_distances`, `invalid_two_cluster_result`, `low_score_separation`, `soft_threshold_accepts_all`, and `no_trainable_after_hard`. This is the real-world abstention path: if the post-hard feature distribution does not contain a credible low-quality partition, the model stops after hard preselection.

## Analysis Output

`cavgs_quality_analysis.txt` starts with `# model_cavgs_rejection_analysis_version=7`. It includes:

- dataset and model identity;
- selected/rejected counts;
- `soft_decision`, `soft_reason`, and the number of trainable rows rejected by the soft model;
- confusion metrics against the manual `cls2D` state;
- score AUC;
- raw threshold, threshold offset, and effective threshold;
- best learn-score and F1 thresholds against the manual reference;
- model fields as comments;
- feature inventory;
- per-feature AUC, robust separation, current weight, and suggested weight;
- threshold scan;
- one row per class average containing state, hard-reject flag, quality cluster, score, manual state, raw features, and normalized features.

`cavgs_quality_features.txt` uses the same per-class row format without the analysis summary comments.

## Learning

`quality_mode=learn` reads analysis files generated by the current analysis table format. Required columns are:

- `dataset_id`
- `hard_reject`
- `manual_state`
- all required `z_*` feature columns for the 14-feature bank

Hard-rejected rows are kept in reports but excluded from feature-weight estimation, candidate scoring, feature-signal diagnostics, feature-drop diagnostics, and feature-policy diagnostics.

Learn mode assigns each dataset an automatic role:

- `balanced`: both manual good and manual bad rows remain after hard rejects; used for feature weights and scored by specificity with a small false-negative-rate tolerance.
- `trainable_good_only`: only manual good rows remain after hard rejects; not used for feature weights, scored by guarded recall.
- `trainable_bad_only`: only manual bad rows remain after hard rejects and no manually good rows were hard rejected; not used for feature weights, scored by specificity.
- `skip_hard_gate_blocked`: only manual bad rows remain because all manually good rows were hard rejected; reported but skipped for soft-model fitting and scoring.

Learn mode derives candidate feature weights from the training data. The starting candidate is:

```text
weight(feature) = max(0, pooled_auc(feature, manual_state) - 0.5)
```

The weights are normalized to sum to one. If all active weights are zero for a policy, uniform weights are assigned over the active policy features. The full classifier search chooses among the AUC-derived candidate vectors for the available feature policies.

Balanced datasets use specificity with a small false-negative-rate tolerance over trainable manually good rows. Overfit-focused datasets also subtract a hinge/quadratic penalty for the accepted manually bad fuzzy-ball-signature rate above the tolerated rate. The final macro score blends the mean dataset score with a fractional lower tail of the dataset scores. This makes threshold selection care about datasets where bad fuzzy-ball classes still leak through, while still penalizing candidates that reject too large a fraction of manually selected classes. Good-only datasets use a recall guard in the learn score. Raw recall is still reported in the dataset table, but the guarded learn score applies an additional shortfall penalty below the good-only recall floor. This prevents a candidate model from improving balanced datasets by rejecting a large fraction of classes from datasets that contain only trainable manual positives.

Learn mode searches:

- feature policies: `microchunk`, `microchunk_plus_score`, `microchunk_plus_signal`, `microchunk_plus_score_signal`;
- feature weights: one AUC-derived candidate for each feature policy;
- `min_score_separation`: `0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50`;
- `boundary_margin`: `-0.60, -0.50, -0.40, -0.30, -0.25, -0.15, -0.05, 0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80`,
  `0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00`;
- `use_lowsep_otsu`: `false, true`;
- `use_otsu_window`: `false, true`;
- `otsu_min_offset`: `0.05, 0.10, 0.15, 0.25, 0.35, 0.40, 0.45, 0.50, 0.60` when the Otsu window is enabled;
- `otsu_max_offset`: `0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.75, 0.90` when the Otsu window is enabled;
- `use_cluster_rescue`: `false, true`;
- `enforce_min_accept_frac`: `false, true`;
- `min_accept_frac`: `0.50, 0.60, 0.65, 0.70, 0.80, 0.85, 0.90, 0.925, 0.95, 0.975, 1.00` when `enforce_min_accept_frac` is true.

The training algorithm is unified and ab initio with respect to the rejection model. It uses a neutral foundation with zero starting weights, no rescue, no minimum-acceptance enforcement, and no Otsu rescue/window behavior. This foundation is only a tie-break anchor and a source of inactive-field defaults; it is not a chunk or pool preset. Chunk, pool, and stage-specific behavior should be represented by the learned model parameters or by promoted model definitions, not by separate training procedures.

Each candidate is evaluated by running the full classifier on every scoreable training dataset and averaging the role-specific learn score. If multiple candidates tie, the selected candidate is the one closest to the neutral foundation in the searched threshold controls and policy flags.

The learn report includes the search grid, `macro_learn_score`, suggested weights, selected model, dataset-role diagnostics, Otsu ablation diagnostics, feature-signal diagnostics, feature-drop diagnostics, leave-one-dataset-out feature-policy diagnostics, top candidates, best ties, and per-dataset confusion metrics.

`quality_mode=evaluate` uses the same trainable-row scoring semantics as learn mode, but it does not derive weights, search thresholds, or write a model file. This is intended for held-out validation. It can score saved analysis tables through `filetab`, or a single manually selected project directly through `projfile` and `mskdiam`. The evaluate report includes the fixed model settings, `macro_evaluate_score`, dataset-role diagnostics, Otsu ablation diagnostics, and per-dataset confusion metrics.

## Promotion

`quality_mode=promote` writes a Fortran snippet for adding a learned model as a built-in preset in `simple_cavg_quality_model.f90`. The snippet also lists the UI/help files that need the model name added when a new built-in is introduced.

## Example Commands

Analyze a manually selected project:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=analyze \
  projfile=my_project.simple \
  mskdiam=180 \
  mkdir=yes
```

Apply the default chunk model:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=apply \
  projfile=my_project.simple \
  mskdiam=180 \
  mkdir=yes
```

Apply only the standard hard gates plus the fixed overfit hard gate:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=apply \
  overfit_hard_reject=yes \
  projfile=my_project.simple \
  mskdiam=180 \
  mkdir=yes
```

Apply only the standard hard gates plus the fixed chunk hard gates:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=apply \
  chunk_hard_reject=yes \
  projfile=my_project.simple \
  mskdiam=180 \
  mkdir=yes
```

Train from a file table of analysis outputs:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=learn \
  model_family=logistic \
  filetab=analysis_files.txt \
  fname=cavgs_quality_model_learned.txt \
  mkdir=yes
```

Evaluate a fixed built-in model on held-out analysis outputs:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=evaluate \
  quality_model=chunk100mics \
  filetab=holdout_analysis_files.txt \
  fname=cavgs_quality_evaluate_report.txt \
  mkdir=yes
```

Evaluate a single manually selected project directly:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=evaluate \
  quality_model=chunk100mics \
  projfile=my_project.simple \
  mskdiam=180 \
  fname=cavgs_quality_evaluate_report.txt \
  mkdir=yes
```

Apply a learned model file:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=apply \
  infile=cavgs_quality_model_learned.txt \
  projfile=my_project.simple \
  mskdiam=180 \
  mkdir=yes
```

Promote a learned model:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=promote \
  infile=cavgs_quality_model_learned.txt \
  fname=cavgs_quality_model_learned_builtin_code.txt \
  mkdir=yes
```
