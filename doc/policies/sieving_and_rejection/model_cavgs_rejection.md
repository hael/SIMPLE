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

The model feature bank keeps the microchunk-style image-processing evidence, optional stored-score and signal evidence, and the local-variance/support evidence needed to detect fuzzy-ball overfitting in learned models. The overfit local-variance features are part of the learned model feature vector and are included in every learn-mode feature policy. Texture descriptors were removed after they failed to justify the extra feature and extraction complexity. The current `chunk100mics` preset uses the `microchunk_plus_score_signal` policy: microchunk plus stored-score, overfit-family, center/edge signal, presence, and fuzzy-ball signal evidence.

| Evidence | Microchunk rule engine | `model_cavgs_rejection` feature-vector model | Current chunk role |
| --- | --- | --- | --- |
| Class population | Hard rule: reject below a tier-specific fraction of total population. | `log_pop` feature plus context-specific population hard gates. | Active feature and hard gate. |
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
| Overfit band-pass localization | Not used by the rule engine. | `log_bp40_100_center_edge_var` measures center/edge variance after 100 to 40 A band-pass; low raw values can be context-specific hard gates, while ordinary low values remain learned evidence. | Active feature in learned policies plus context-specific hard floors. |
| Fuzzy-ball signal | Not used by the rule engine. | `fuzzy_ball_signal` combines foreground texture, central presence, and 100 to 40 A center/edge localization. Low values mark the fuzzy-ball pattern as continuous learned evidence, not as an apply-time hard gate. | Active feature in learned policies. |
| Invalid pixels | Not a named microchunk rule. | Image values containing invalid pixels are hard rejected. | Hard gate. |

## Hard-Reject Comparison

Hard rejects in `model_cavgs_rejection` are selected by `quality_context=chunk|pool|sieve`, with `chunk` as the default.

The contexts represent three workflow phases:

- `sieve`: very small 2D runs on small particle chunks. This phase uses conservative hard gates only and deliberately avoids learned rejection.
- `chunk`: the next 2D stage, where sieve-cleaned particles are processed in larger chunks, typically 10-30k particles. This phase uses hard gates plus the chunk logistic model.
- `pool`: the final 2D phase before 3D analysis, where highly clean particle sets from chunk modeling have been merged. This phase uses its own pool model.

Chunk and pool share non-negotiable validity gates. Chunk adds early-streaming cleanup gates for undersupported and fuzzy-ball-like class averages. Pool adds its own final pre-3D cleanup gates for low population, low band-pass localization, and poor nominal resolution. Sieve is a separate compatibility-analysis route with its own hard-gate policy and deliberately does not inherit the shared chunk/pool validity gates.

| Criterion | Microchunk rule engine | `model_cavgs_rejection` |
| --- | --- | --- |
| Population | Rejects `pop < ceiling(sum(pop) * fraction)`. Fractions are tier-specific. | Shared: reject `pop <= 0`. Chunk/sieve: reject `pop < ceiling(sum(pop) * 0.0035)`. Pool: reject `pop < ceiling(sum(pop) * 5.0e-4)`. |
| Resolution | Rejects `res > 40.0`. | Chunk/pool shared validity rejects `res > 40.0`; pool additionally rejects `res > 25.0`. |
| Empty or mismatched class-average stack | Stream chunk is marked `REJECTION_FAILED` and `COMPLETE`. | Apply/analyze require a valid project, `cls2D` entries, matching class-average stack, and `mskdiam`; invalid inputs stop the command. |
| Invalid pixels | No separate named rule. | Hard rejected. |
| No valid foreground component | Rejected after full-image connected components are pruned. | Hard rejected after the same style of foreground-component pruning. |
| Foreground centroid outside mask | Rejected. | Hard rejected. |
| Largest foreground component outside mask | Rejected when outside pixels exceed 10. | Hard rejected when outside pixels exceed 10. |
| Local variance exactly degenerate | Rejected when both inside/outside scores are near zero. | Hard rejected when both foreground/background local variances are non-positive. |
| Local variance low but not degenerate | Rejected by fixed robust-z thresholds. | Chunk only: reject an extreme absolute foreground local-variance floor; otherwise encoded as learned evidence. Pool does not apply this extra gate. |
| Band-pass center/edge variance extremely low | Not used by the rule engine. | Chunk/sieve hard reject when raw `bp_center_edge_var < 1.5`; pool hard rejects when raw `bp_center_edge_var < 10.0`; less extreme values remain learned evidence. |

Microchunk tier thresholds:

| Tier | Population fraction | Local-variance strong threshold | Local-variance weak threshold |
| --- | ---: | ---: | ---: |
| Pass 1 | `0.0050` | `-0.5` | `-0.1` |
| Pass 2 | `0.0035` | `-1.0` | `-1.0` |
| Reference chunk | `0.0025` | `-2.0` | `-2.0` |
| Match chunk | `0.0025` | `-2.0` | `-2.0` |

`model_cavgs_rejection` uses context-specific population hard-reject fractions: `0.0035` for chunk/sieve and `5.0e-4` for pool.

## Learning Model

The rejection model is a low-dimensional, monotone feature-vector classifier with explicit hard constraints. In learning-theory terms, the hard rejects are prior validity constraints: they remove classes the model is not asked to rescue or fit. The remaining rows form a supervised training set from manual selections.

The legacy linear apply-only artifact has three parts:

- a feature policy, which selects a cumulative family set;
- non-negative feature weights, which define a linear scalar quality score;
- thresholding controls, which govern how the per-dataset score boundary is chosen after k-medoids clustering and optional Otsu thresholding.

New training fits a relational pairwise-logistic artifact. Its fitted parameters are an intercept, linear feature coefficients, pairwise base-feature interaction coefficients, the required relational coefficient, a probability threshold, and a regularization strength. The coefficients are fit from the canonical training files only. On two-class trainable datasets, the logistic loss gives manually selected classes moderately higher total weight than manually deselected classes so the probability surface learns rejection evidence without becoming an overly aggressive rejector. Manually deselected examples with the overfit/support-poor or band-pass-localization signature get an additional training-time loss multiplier, so fuzzy-ball examples influence the fit without becoming an apply-time hard reject. Learn reports also expose an internal overfit-focus loss scale for experiments; it is currently neutral.

All features are robustly normalized inside the current dataset. Legacy linear
artifacts retain their k-medoids/Otsu apply path. Relational logistic artifacts
apply their learned probability surface and fixed probability threshold.

The learning objective is empirical risk minimization over manually annotated `quality_mode=analyze` runs. Feature-weight candidates are learned only from datasets where both manual states remain present after hard rejects, because those are the datasets with a learnable soft boundary. Candidate models are scored over all scoreable datasets using only non-hard-rejected rows. Two-class trainable datasets contribute specificity with a small false-negative-rate tolerance over the trainable manually good rows. Datasets where at least 20% of the trainable class averages are manually bad rows matching the fuzzy-ball signature get an additional hinge/quadratic objective penalty when the accepted fuzzy-ball-signature rate is above the tolerated rate. The signature is poor support plus either low local-variance evidence or low 100 to 40 A band-pass center/edge localization. This keeps selected-class protection global while making the objective care specifically about fuzzy-ball leakage in small-specimen chunks. The final macro score is a fixed robust score: half the mean dataset score and half the mean score over the lower tail of datasets. During pairwise-logistic fitting, manually bad rows with the same fuzzy-ball signature get extra loss weight. Datasets where only manually good rows remain after hard rejects contribute recall, so they teach the model not to reject good classes that passed the hard gates. Datasets where only manually bad rows remain contribute specificity only when no manually good classes were removed by hard rejects. If manually good classes were entirely removed by hard rejects, the dataset is reported as hard-gate blocked and skipped for soft-model scoring. Hard-rejected rows remain visible in diagnostics, including counts of manually good classes lost to hard gates, but they do not participate in fitting the learned boundary.

Feature-family policies act as a small structural model-selection layer. The
relational logistic learner fits coefficients independently for each policy
and selects the feature policy, regularization lambda, and probability
threshold from the supplied training data.

Use `filetab=` to define the target workflow scenario being trained. For example, early streaming chunk models should be fit from chunk analyses, while late pooled-refinement models should be fit from pool analyses. To inspect cross-scenario behavior without influencing the model, train the target model and then run `quality_mode=evaluate filetab=... infile=...` on the other scenario.

| Policy | Active families | Purpose |
| --- | --- | --- |
| `microchunk` | `microchunk`, `overfit` | Tests whether scalar analogues of the stream rejector plus overfit local-variance evidence are sufficient. |
| `microchunk_plus_score` | `microchunk`, `overfit`, `score` | Adds stored class correlation/FRC-like evidence. |
| `microchunk_plus_signal` | `microchunk`, `overfit`, `signal` | Adds center/edge and presence evidence without stored score evidence. |
| `microchunk_plus_score_signal` | `microchunk`, `overfit`, `score`, `signal` | Uses the full compact quality feature set. |

Learn mode reports feature signal, feature-drop diagnostics, and leave-one-dataset-out feature-policy diagnostics so that the selected model can be interpreted as a compact scientific statement about which evidence family and weight structure generalized on the supplied training set.

## Command Modes

`quality_mode=apply` computes quality features, applies the selected model, writes `cavgs_quality_features.txt`, writes `quality_selected_cavgs.mrc`, `quality_rejected_cavgs.mrc`, `hard_gate_rejections.mrc`, `quality_ranked_cavgs.mrc`, and `quality_ranked_cavgs.txt`, maps the selected class averages into particle states, annotates `cls2D` with `quality`, `quality_cluster`, and `accept`, annotates `cls3D` when its class count matches `cls2D`, optionally prunes particles with `prune=yes`, and writes the project.

`quality_mode=analyze` computes the same model output but treats the existing `cls2D` state as the manual reference. It writes the single canonical learner input `cavgs_quality_training.txt`, the selected/rejected stacks, `hard_gate_rejections.mrc`, and the score-ranked class-average stack plus rank table. The project selection is left unchanged.

`quality_mode=learn` reads a training file table of `cavgs_quality_training.txt` files from `filetab=` and fits the relational logistic model from a neutral `abinitio_learn_base` foundation. Every input must declare `relational_feature_schema=corr_knn_signal_v1` and contain the normalized CC-neighbour signal-statistics feature; missing or mixed relational schemas are rejected. `quality_context=chunk|pool` labels the learned version-10 model context. `quality_context=sieve` is intentionally rejected because sieve is a hard-gates-only screening phase. Learn mode does not accept `quality_model` or `infile` as a seed. It writes a learned model file controlled by `fname=` and writes `cavgs_quality_learn_report.txt`.

`quality_mode=evaluate` applies the selected fixed model without refitting. With `filetab=`, it evaluates one or more saved `cavgs_quality_training.txt` files. Without `filetab=`, it evaluates a single project directly using the existing `cls2D` state as the manual reference, like analyze mode. It writes `cavgs_quality_evaluate_report.txt`, or the report path controlled by `fname=`.

`quality_mode=promote` reads a model file from `infile=` and writes a Fortran promotion snippet controlled by `fname=`.

`apply` and `analyze` require `projfile` and `mskdiam`. `learn` requires `filetab`. `evaluate` requires either `filetab` or `projfile` plus `mskdiam`. `promote` requires `infile`. The commander sets `oritype=cls2D`, defaults `mkdir=yes`, and defaults `prune=no`.

For `apply`, `analyze`, and project-backed `evaluate`, the command uses `chunk100mics` unless `quality_model` or `infile` is supplied. If `quality_context` is omitted, the command uses the loaded model context and also falls back to model-name inference from `pool`, `sieve`, or `chunk` substrings. Project-backed runs with `quality_context=sieve` use the hard-gates-only evaluator and skip learned model scoring. Saved-analysis `evaluate filetab=...` runs do not write image stacks because they do not load a project or class-average stack.

## Model Selection

`quality_model` selects the built-in preset outside learn mode. The promoted built-in is:

- `chunk100mics`: default chunk/stream-style version-10 pairwise logistic model trained from `/Users/elmlundho/model_cavgs_rejection/chunk_training5`; it includes `corr_knn_signal_v1` relational evidence.

When `infile` is supplied, the model file is treated as a complete model and wins over the built-in preset.

Relational logistic model files use `model_version=10` and explicit key-value fields:

- `name`
- `context`
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

Resolution is intentionally not controlled by a separate training switch. If nominal resolution is unreliable for a particular scenario, the regularized linear or pairwise-logistic fit is expected to down-weight it or use it conditionally through interactions with other features.

Features outside the selected policy are encoded by zero weights.

## Extraction Constants

- `FOREGROUND_SEG_LP = 30.0`
- `SIGNAL_METRIC_LP = 10.0`
- `OVERFIT_SIGNAL_BP_HP = 100.0`
- `OVERFIT_SIGNAL_BP_LP = 40.0`
- `CAVG_RES_HARD_REJECT_A = 40.0`
- `POP_FRACTION_HARD_REJECT = 0.0035`
- `BP_CENTER_EDGE_VAR_HARD_REJECT_MIN = 1.5`
- `CHUNK_LOCVAR_FG_HARD_REJECT_MAX = exp(-4.5)`
- `POOL_RES_HARD_REJECT_A = 25.0`
- `POOL_POP_FRACTION_HARD_REJECT = 5.0e-4`
- `POOL_BP_CENTER_EDGE_VAR_HARD_REJECT_MIN = 10.0`
- `LOCVAR_WINDOW = 10`
- `MASK_HARD_OUTSIDE_PIXELS = 10`
- `LOG_EPS = 1.0e-12`
- `CLIP_Z = 4.0`

## Hard Rejects

Hard rejects are context-sensitive gates. For `chunk` and `pool`, they run before model fitting, clustering, and training-score calculation. For `sieve`, they are the whole rejection policy.

A class average is hard rejected in `chunk` and `pool` when any shared validity condition holds:

- `pop <= 0`;
- `res > 40.0`;
- the image contains invalid pixels;
- foreground segmentation has no valid component after full-image connected components are pruned;
- a segmented foreground component centroid lies outside the mask radius;
- the largest foreground component has more than 10 pixels outside the mask disc;
- foreground and background local variance are both non-positive.

For `quality_context=chunk`, these additional early-streaming gates also apply:

- `pop < ceiling(sum(pop) * 0.0035)`;
- raw foreground local variance is below `exp(-4.5)`;
- raw 100 to 40 A band-pass center/edge variance is below `1.5`.

For `quality_context=pool`, the chunk local-variance floor is not applied. Pool instead adds stricter final pre-3D cleanup gates learned from the current pool training set:

- `res > 25.0`;
- `pop < ceiling(sum(pop) * 5.0e-4)`;
- raw 100 to 40 A band-pass center/edge variance is below `10.0`.

For `quality_context=sieve`, the route intentionally uses a separate hard-gate policy for compatibility analysis. It does not inherit the shared chunk/pool validity gates, and it does not run a learned model. The sieve-specific gates are:

- `pop < ceiling(sum(pop) * 0.0035)`;
- raw 100 to 40 A band-pass center/edge variance is below `1.5`.

Hard-rejected classes receive rejected state directly. Their normalized features are set to `-CLIP_Z`, their model scores are set to `-CLIP_Z`, and they remain visible in analysis and learning reports.

The band-pass center/edge floor is deliberately conservative. In the v4 chunk training tables, `bp_center_edge_var < 1.5` preserved all manually selected classes while still removing additional manually rejected classes beyond population and resolution. More aggressive fuzzy-ball evidence should remain in the normalized feature vector and learned model, not in the default pre-training hard rejects.

Project-backed `apply`, `analyze`, and `evaluate` runs also write `hard_gate_rejections.mrc`. For `chunk` and `pool`, this is a bookkeeping stack containing classes removed before the model stage. For `sieve`, it is the rejection decision stack because the context is hard-gates-only.

## Normalization

`normalize_cavg_quality_features` computes per-feature robust statistics over non-hard-rejected rows only. Each feature is centered by the median and scaled by the Gaussian MAD. Values are clipped to `[-CLIP_Z, CLIP_Z]`. If a feature has no robust spread, its normalized value is set to zero for trainable rows.

## Built-In Presets

`chunk100mics` uses feature policy `microchunk_plus_score_signal` and the relational logistic model. It was trained from `/Users/elmlundho/model_cavgs_rejection/chunk_training5`.

```text
prob_threshold      3.500000E-01
regularization      1.000000E-03
feature_weights     uniform over all 14 microchunk_plus_score_signal features
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

On the refreshed v4 chunk-training table, this promoted preset scored `macro_evaluate_score=0.50847`, improving over the previous built-in chunk preset (`-0.80859`). The main gain was selected-class protection: soft-classification totals moved from `tp=297, fp=95, tn=151, fn=33` to `tp=324, fp=83, tn=163, fn=6`.

## Classification

For each non-hard-rejected class average, the relational logistic model computes:

```text
prob_accept = sigmoid(intercept + sum(linear_coefficients * features)
                    + sum(interaction_coefficients * feature_i * feature_j)
                    + relational_coefficient * relational_feature)
```

Every supported artifact declares the relational schema. Non-hard-rejected rows are accepted when `prob_accept >= prob_threshold`.
Standard hard gates still run before the logistic model.

If there are fewer than four trainable rows, the distance matrix is degenerate, clustering fails, or the low-separation branch has no acceptable Otsu threshold, all non-hard-rejected rows are accepted as a single cluster.

`threshold_offset` is reported as `raw_threshold - threshold`. A positive value means the effective threshold is lower than the cluster midpoint.

The classifier reports its soft decision explicitly:

- `soft_decision=soft_threshold`: a learned score threshold was used to reject at least one non-hard-rejected class average.
- `soft_decision=hard_only`: no additional soft rejection was made after hard preselection.

Common hard-only reasons are `too_few_trainable`, `flat_feature_distances`, `invalid_two_cluster_result`, `low_score_separation`, `soft_threshold_accepts_all`, and `no_trainable_after_hard`. This is the real-world abstention path: if the post-hard feature distribution does not contain a credible low-quality partition, the model stops after hard preselection.

## Analysis Output

`cavgs_quality_training.txt` starts with
`# cavg_quality_training_version=1`. It includes:

- dataset and model identity;
- selected/rejected counts;
- `soft_decision`, `soft_reason`, and the number of trainable rows rejected by the soft model;
- confusion metrics against the manual `cls2D` state;
- relational schema and provider policy when supported by the model;
- model fields as comments;
- feature inventory;
- per-feature AUC, robust separation, current weight, and suggested weight;
- threshold scan;
- one row per class average containing state, hard-reject flag, quality score,
  manual state, raw/normalized base features, and raw/normalized promoted
  relational feature.

`cavgs_quality_features.txt` uses the same per-class row format without the analysis summary comments.

## Learning

`quality_mode=learn` reads canonical files generated by the current analyze format. Required columns are:

- `dataset_id`
- `hard_reject`
- `manual_state`
- all required `z_*` feature columns for the 14-feature bank

Every file must declare `relational_feature_schema=corr_knn_signal_v1`, contain
`z_signal_stats_anchor_topk_mean`, and agree on `k`, high-pass, low-pass, and
shift-range policy. Analyze guarantees this even when its scoring model is an
older base-only artifact.

Hard-rejected rows are kept in reports but excluded from feature-weight estimation, candidate scoring, feature-signal diagnostics, feature-drop diagnostics, and feature-policy diagnostics.

Learn mode assigns each dataset an automatic role:

- `balanced`: both manual good and manual bad rows remain after hard rejects; used for feature weights and scored by specificity with a small false-negative-rate tolerance.
- `trainable_good_only`: only manual good rows remain after hard rejects; not used for feature weights, scored by guarded recall.
- `trainable_bad_only`: only manual bad rows remain after hard rejects and no manually good rows were hard rejected; not used for feature weights, scored by specificity.
- `skip_hard_gate_blocked`: only manual bad rows remain because all manually good rows were hard rejected; reported but skipped for soft-model fitting and scoring.

Balanced datasets use specificity with a small false-negative-rate tolerance over trainable manually good rows. Overfit-focused datasets also subtract a hinge/quadratic penalty for the accepted manually bad fuzzy-ball-signature rate above the tolerated rate. The final macro score blends the mean dataset score with a fractional lower tail of the dataset scores. This makes threshold selection care about datasets where bad fuzzy-ball classes still leak through, while still penalizing candidates that reject too large a fraction of manually selected classes. Good-only datasets use a recall guard in the learn score. Raw recall is still reported in the dataset table, but the guarded learn score applies an additional shortfall penalty below the good-only recall floor. This prevents a candidate model from improving balanced datasets by rejecting a large fraction of classes from datasets that contain only trainable manual positives.

The learner fits an intercept, selected base-feature coefficients, all pairwise
interactions among those selected base features, and one relational
coefficient. The relational feature remains a standalone term; it is not
expanded into another interaction bank. Learn mode searches:

- feature policies: `microchunk`, `microchunk_plus_score`, `microchunk_plus_signal`, `microchunk_plus_score_signal`;
- L2 regularization: `0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1`;
- probability thresholds from `0.02` through `0.90` on the fixed grid reported
  by the learner.

The fit is ab initio with respect to the rejection model. Chunk, pool, and
stage-specific behaviour is represented by the learned artifact context and
coefficients, not by separate fitting procedures. A relational pool-context
fit therefore produces a version-10 artifact even though the existing built-in
pool preset remains a version-9-compatible base-only model.

Each candidate is evaluated on every scoreable training dataset using the
role-specific learn score. Objective value breaks score ties.

The learn report includes the search grid, `macro_learn_score`, fitted base and
relational coefficients, dataset-role diagnostics, and per-dataset confusion
metrics.

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

Train from a file table of canonical outputs:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=learn \
  filetab=cavgs_quality_training_filetab.txt \
  fname=cavgs_quality_model_learned.txt \
  mkdir=yes
```

Evaluate a fixed built-in model on independent canonical outputs:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=evaluate \
  quality_model=chunk100mics \
  filetab=holdout_training_files.txt \
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
