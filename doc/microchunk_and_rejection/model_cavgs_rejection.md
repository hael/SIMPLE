# model_cavgs_rejection

## Scope

`model_cavgs_rejection` is the SIMPLE class-average quality selection program. It evaluates `cls2D` class averages, applies hard validity rejects, extracts a fixed scalar feature bank, normalizes those features within the dataset, and applies a named quality model to partition class averages into accepted and rejected sets. The model can also abstain from additional soft rejection and report that the final decision was hard preselection only.

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

The model feature bank keeps the microchunk-style image-processing evidence as a feature family and appends optional evidence families. The current `chunk_default_v2` preset uses the `microchunk_plus_score` policy: the microchunk-family features plus `corr_frc_proxy`.

| Evidence | Microchunk rule engine | `model_cavgs_rejection` feature-vector model | Current chunk role |
| --- | --- | --- | --- |
| Class population | Hard rule: reject below a tier-specific fraction of total population. | `log_pop` feature plus hard reject below `0.0035` of total population. | Active feature and hard gate. |
| Resolution | Hard rule: reject when `res > 40.0`. | `neg_log_res` feature plus hard reject when `res > 40.0`. | Active feature and hard gate. |
| Foreground centering | Otsu foreground connected components; reject if any centroid lies outside the mask radius. | `centered` feature is negative normalized centroid displacement; centroid outside the mask is a hard reject. | Active feature and hard gate. |
| Foreground outside mask | Reject when the largest valid foreground component has more than 10 pixels outside the mask disc. | `cc_area_frac` measures largest component area relative to the mask disc; outside-mask excess is a hard reject. | Active feature and hard gate. |
| Full-image component cleanup | Connected components spanning the full image are removed before mask tests. | Same concept is applied during feature extraction before geometry features and hard mask rejection. | Hard-gate preprocessing. |
| Foreground local variance | Robust-z rule on foreground/background local variance after 10 A low-pass and Otsu masking. | `log_locvar_fg` feature; only both variances non-positive is a hard reject. | Active feature. |
| Background local variance | Robust-z rule paired with foreground local variance. | `log_locvar_bg` feature; only both variances non-positive is a hard reject. | Active feature. |
| Stored score evidence | Not used by the rule engine. | `corr_frc_proxy` reads stored class correlation or FRC-like score when present. | Active feature. |
| Center/edge signal | Not used by the rule engine. | `log_center_edge_snr` feature in the `signal` family. | Available, zero weight in current chunk default. |
| Presence | Not used by the rule engine. | `presence` feature in the `signal` family. | Available, zero weight in current chunk default. |
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

`model_cavgs_rejection` uses a fixed population hard-reject fraction of `0.0035` for the current feature extractor, independent of model context.

## Learning Model

The rejection model is a low-dimensional, monotone feature-vector classifier with explicit hard constraints. In learning-theory terms, the hard rejects are prior validity constraints: they remove classes the model is not asked to rescue or fit. The remaining rows form a supervised training set from manual selections.

The learned artifact has three parts:

- a feature policy, which selects a cumulative family set;
- non-negative feature weights, which define a linear scalar quality score;
- thresholding controls, which govern how the per-dataset score boundary is chosen after k-medoids clustering and optional Otsu thresholding.

The apply-time classifier is partly dataset-adaptive. Features are robustly normalized inside the current dataset. K-medoids and Otsu operate on the current score/feature distribution. The learned parameters provide the inductive bias: which features are active, how the score is weighted, and how aggressively the cluster-derived boundary is shifted or replaced.

The learning objective is empirical risk minimization over manually annotated `quality_mode=analyze` runs. Feature weights are learned only from datasets where both manual states remain present after hard rejects, because those are the datasets with a learnable soft boundary. Candidate models are scored over all scoreable datasets using only non-hard-rejected rows. Two-class trainable datasets contribute balanced accuracy. Datasets where only manually good rows remain after hard rejects contribute recall, so they teach the model not to reject good classes that passed the hard gates. Datasets where only manually bad rows remain contribute specificity only when no manually good classes were removed by hard rejects. If manually good classes were entirely removed by hard rejects, the dataset is reported as hard-gate blocked and skipped for soft-model scoring. Hard-rejected rows remain visible in diagnostics, including counts of manually good classes lost to hard gates, but they do not participate in fitting the learned boundary.

Feature weights are learned ab initio from the supplied training data:

```text
weight(feature) = max(0, pooled_auc(feature, manual_state) - 0.5)
```

The weights are normalized after the active feature policy is applied. There is no base-weight blending in the current learner. This means each training run derives the scalar model from the current analysis files rather than nudging an inherited default.

Feature-family policies act as a small structural model-selection layer:

| Policy | Active families | Purpose |
| --- | --- | --- |
| `microchunk` | `microchunk` | Tests whether scalar analogues of the stream rejector evidence are sufficient. |
| `microchunk_plus_score` | `microchunk`, `score` | Adds stored class correlation/FRC-like evidence. This is the current chunk default. |
| `microchunk_plus_signal` | `microchunk`, `signal` | Adds center/edge and presence evidence without stored score evidence. |
| `microchunk_plus_score_signal` | `microchunk`, `score`, `signal` | Uses the full current feature bank. |

Learn mode reports feature signal, feature-drop diagnostics, and leave-one-dataset-out feature-policy diagnostics so that the selected model can be interpreted as a compact scientific statement about which evidence family generalized on the supplied training set.

## Command Modes

`quality_mode=apply` computes quality features, applies the selected model, writes `cavgs_quality_features.txt`, writes `quality_selected_cavgs.mrc` and `quality_rejected_cavgs.mrc`, maps the selected class averages into particle states, annotates `cls2D` with `quality`, `quality_cluster`, and `accept`, annotates `cls3D` when its class count matches `cls2D`, optionally prunes particles with `prune=yes`, and writes the project.

`quality_mode=analyze` computes the same model output but treats the existing `cls2D` state as the manual reference. It writes `cavgs_quality_analysis.txt` and the selected/rejected stacks. The project selection is left unchanged.

`quality_mode=learn` reads a file table of `cavgs_quality_analysis.txt` files and searches for a model specification. It writes a learned model file controlled by `fname=` and writes `cavgs_quality_learn_report.txt`.

`quality_mode=promote` reads a model file from `infile=` and writes a Fortran promotion snippet controlled by `fname=`.

`apply` and `analyze` require `projfile` and `mskdiam`. `learn` requires `filetab`. `promote` requires `infile`. The commander sets `oritype=cls2D`, defaults `mkdir=yes`, and defaults `prune=no`.

## Model Selection

`quality_model` selects a built-in preset. The default is `chunk_default_v2`. The available built-ins are:

- `chunk_default_v2`: default chunk/stream-style model.
- `pool_default_v1`: pool/batch-style model with minimum accepted fraction enforcement.
- `pool_exp`: experimental pool/batch-style model learned from the widened pool search grid.

When `infile` is supplied, the model file is treated as a complete model and wins over the built-in preset.

Model files use `model_version=5` and explicit key-value fields:

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

`CAVG_QUALITY_NFEATS` is 9. All feature directions are `higher_is_better`. The model consumes normalized `z_*` values; raw values are written for diagnostics.

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

Feature families are used by learn mode. Every policy starts with the `microchunk` family and then appends optional evidence families:

- `microchunk`
- `microchunk_plus_score`
- `microchunk_plus_signal`
- `microchunk_plus_score_signal`

Features outside the selected policy are encoded by zero weights.

## Extraction Constants

- `FOREGROUND_SEG_LP = 30.0`
- `SIGNAL_METRIC_LP = 10.0`
- `CAVG_RES_HARD_REJECT_A = 40.0`
- `POP_FRACTION_HARD_REJECT = 0.0035`
- `LOCVAR_WINDOW = 10`
- `MASK_HARD_OUTSIDE_PIXELS = 10`
- `LOG_EPS = 1.0e-12`
- `CLIP_Z = 4.0`

## Hard Rejects

Hard rejects are validity gates applied before model fitting, clustering, and training-score calculation. A class average is hard rejected when any of these conditions hold:

- `pop <= 0`;
- `pop < ceiling(sum(pop) * 0.0035)`;
- `res > 40.0`;
- the image contains invalid pixels;
- foreground segmentation has no valid component after full-image connected components are pruned;
- a segmented foreground component centroid lies outside the mask radius;
- the largest foreground component has more than 10 pixels outside the mask disc;
- foreground and background local variance are both non-positive.

Hard-rejected classes receive rejected state directly. Their normalized features are set to `-CLIP_Z`, their model scores are set to `-CLIP_Z`, and they remain visible in analysis and learning reports.

## Normalization

`normalize_cavg_quality_features` computes per-feature robust statistics over non-hard-rejected rows only. Each feature is centered by the median and scaled by the Gaussian MAD. Values are clipped to `[-CLIP_Z, CLIP_Z]`. If a feature has no robust spread, its normalized value is set to zero for trainable rows.

## Built-In Presets

`chunk_default_v2` has context `chunk`, feature policy `microchunk_plus_score`, cluster rescue disabled, and minimum accepted fraction disabled.

```text
log_pop             1.355581E-01
neg_log_res         1.808395E-01
centered            4.942806E-02
log_locvar_fg       1.635655E-01
log_locvar_bg       1.726983E-01
corr_frc_proxy      2.108072E-01
log_center_edge_snr 0.000000E+00
cc_area_frac        8.710331E-02
presence            0.000000E+00
```

```text
boundary_margin         0.15
min_score_separation    0.15
otsu_min_offset         0.35
otsu_max_offset         0.50
cluster_rescue_margin   0.20
min_accept_frac         0.00
use_lowsep_otsu         true
use_otsu_window         true
use_cluster_rescue      false
enforce_min_accept_frac false
```

`pool_default_v1` has context `pool`, feature policy `microchunk_plus_score_signal`, cluster rescue enabled, and minimum accepted fraction enforcement enabled.

```text
log_pop             3.953488E-01
neg_log_res         2.093023E-01
centered            0.000000E+00
log_locvar_fg       1.860465E-01
log_locvar_bg       2.093023E-01
corr_frc_proxy      0.000000E+00
log_center_edge_snr 0.000000E+00
cc_area_frac        0.000000E+00
presence            0.000000E+00
```

```text
boundary_margin         0.05
min_score_separation    0.15
otsu_min_offset         0.25
otsu_max_offset         0.50
cluster_rescue_margin   0.20
min_accept_frac         0.65
use_lowsep_otsu         false
use_otsu_window         false
use_cluster_rescue      true
enforce_min_accept_frac true
```

`pool_exp` has context `pool`, feature policy `microchunk_plus_signal`, cluster rescue enabled, low-separation Otsu enabled, and minimum accepted fraction enforcement enabled. It is an experimental validation preset, not the pool default.

```text
log_pop             2.029138E-01
neg_log_res         2.802879E-01
centered            4.939443E-02
log_locvar_fg       8.316658E-02
log_locvar_bg       9.741970E-02
corr_frc_proxy      0.000000E+00
log_center_edge_snr 1.498684E-01
cc_area_frac        4.935930E-02
presence            8.758996E-02
```

```text
boundary_margin         0.80
min_score_separation    0.20
otsu_min_offset         0.25
otsu_max_offset         0.50
cluster_rescue_margin   0.20
min_accept_frac         0.90
use_lowsep_otsu         true
use_otsu_window         false
use_cluster_rescue      true
enforce_min_accept_frac true
```

All model weights are clipped to non-negative values and normalized to sum to one when a model is initialized or read.

## Classification

For non-hard-rejected class averages:

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

If there are fewer than four trainable rows, the distance matrix is degenerate, clustering fails, or the low-separation branch has no acceptable Otsu threshold, all non-hard-rejected rows are accepted as a single cluster.

`threshold_offset` is reported as `raw_threshold - threshold`. A positive value means the effective threshold is lower than the cluster midpoint.

The classifier reports its soft decision explicitly:

- `soft_decision=soft_threshold`: a learned score threshold was used to reject at least one non-hard-rejected class average.
- `soft_decision=hard_only`: no additional soft rejection was made after hard preselection.

Common hard-only reasons are `too_few_trainable`, `flat_feature_distances`, `invalid_two_cluster_result`, `low_score_separation`, `soft_threshold_accepts_all`, and `no_trainable_after_hard`. This is the real-world abstention path: if the post-hard feature distribution does not contain a credible low-quality partition, the model stops after hard preselection.

## Analysis Output

`cavgs_quality_analysis.txt` starts with `# model_cavgs_rejection_analysis_version=4`. It includes:

- dataset and model identity;
- selected/rejected counts;
- `soft_decision`, `soft_reason`, and the number of trainable rows rejected by the soft model;
- confusion metrics against the manual `cls2D` state;
- score AUC;
- raw threshold, threshold offset, and effective threshold;
- best balanced-accuracy and F1 thresholds against the manual reference;
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
- all `z_*` feature columns for the 9-feature bank

Hard-rejected rows are kept in reports but excluded from feature-weight estimation, candidate scoring, feature-signal diagnostics, feature-drop diagnostics, and feature-policy diagnostics.

Learn mode assigns each dataset an automatic role:

- `balanced`: both manual good and manual bad rows remain after hard rejects; used for feature weights and scored by balanced accuracy.
- `trainable_good_only`: only manual good rows remain after hard rejects; not used for feature weights, scored by recall.
- `trainable_bad_only`: only manual bad rows remain after hard rejects and no manually good rows were hard rejected; not used for feature weights, scored by specificity.
- `skip_hard_gate_blocked`: only manual bad rows remain because all manually good rows were hard rejected; reported but skipped for soft-model fitting and scoring.

Learn mode derives feature weights from the training data:

```text
weight(feature) = max(0, pooled_auc(feature, manual_state) - 0.5)
```

The weights are normalized to sum to one. If all active weights are zero for a policy, uniform weights are assigned over the active policy features.

Learn mode searches:

- feature policies: `microchunk`, `microchunk_plus_score`, `microchunk_plus_signal`, `microchunk_plus_score_signal`;
- `min_score_separation`: `0.05, 0.10, 0.15, 0.20, 0.30`;
- chunk `boundary_margin`: `-0.60, -0.50, -0.40, -0.30, -0.25, -0.15, -0.05, 0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50`;
- pool `boundary_margin`: `-0.60, -0.50, -0.40, -0.30, -0.25, -0.15, -0.05, 0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80`;
- `use_lowsep_otsu`: `false, true`;
- `use_otsu_window`: `false, true`;
- `otsu_min_offset`: `0.05, 0.10, 0.15, 0.25, 0.35` when the Otsu window is enabled;
- `otsu_max_offset`: `0.40, 0.50, 0.65` when the Otsu window is enabled;
- `min_accept_frac`: `0.50, 0.60, 0.65, 0.70, 0.80, 0.85, 0.90` for pool-context training.

Each candidate is evaluated by running the full classifier on every scoreable training dataset and averaging the role-specific learn score. If multiple candidates tie, the selected candidate is the one closest to the starting model in the searched threshold controls.

The learn report includes the search grid, `macro_learn_score`, suggested weights, selected model, dataset-role diagnostics, Otsu ablation diagnostics, feature-signal diagnostics, feature-drop diagnostics, leave-one-dataset-out feature-policy diagnostics, top candidates, best ties, and per-dataset confusion metrics.

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

Train from a file table of analysis outputs:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=learn \
  filetab=analysis_files.txt \
  fname=cavgs_quality_model_chunk_learned.txt \
  mkdir=yes
```

Apply a learned model file:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=apply \
  infile=cavgs_quality_model_chunk_learned.txt \
  projfile=my_project.simple \
  mskdiam=180 \
  mkdir=yes
```

Promote a learned model:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=promote \
  infile=cavgs_quality_model_chunk_learned.txt \
  fname=cavgs_quality_model_chunk_builtin_code.txt \
  mkdir=yes
```
