# model_cavgs_rejection

## Scope

`model_cavgs_rejection` is the SIMPLE class-average quality selection program. It evaluates `cls2D` class averages, applies hard validity rejects, extracts a fixed scalar feature bank, normalizes those features within the dataset, and applies a named quality model to partition class averages into accepted and rejected sets.

The implementation lives in `src/main/cavg_quality`. The command entry point is `exec_model_cavgs_rejection` in `src/main/commanders/simple/simple_commanders_cavgs.f90`.

## Implementation Files

- `simple_cavg_quality_types.f90`: shared constants and derived types.
- `simple_cavg_quality_feats.f90`: feature definitions, feature extraction, hard rejects, and robust normalization.
- `simple_cavg_quality_stats.f90`: binary metrics, AUC, distance-matrix normalization, and statistical helpers.
- `simple_cavg_quality_model.f90`: built-in presets, model file I/O, scoring, clustering, thresholding, and promotion snippets.
- `simple_cavg_quality_analysis.f90`: evaluation, analysis reports, and feature tables.
- `simple_cavg_quality_learn.f90`: analysis-table reader, feature-policy search, model search, and learn reports.

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

## Analysis Output

`cavgs_quality_analysis.txt` starts with `# model_cavgs_rejection_analysis_version=4`. It includes:

- dataset and model identity;
- selected/rejected counts;
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

Learn mode derives feature weights from the training data:

```text
weight(feature) = max(0, pooled_auc(feature, manual_state) - 0.5)
```

The weights are normalized to sum to one. If all active weights are zero for a policy, uniform weights are assigned over the active policy features.

Learn mode searches:

- feature policies: `microchunk`, `microchunk_plus_score`, `microchunk_plus_signal`, `microchunk_plus_score_signal`;
- `min_score_separation`: `0.05, 0.10, 0.15, 0.20, 0.30`;
- `boundary_margin`: `-0.60, -0.50, -0.40, -0.30, -0.25, -0.15, -0.05, 0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50`;
- `use_lowsep_otsu`: `false, true`;
- `use_otsu_window`: `false, true`;
- `otsu_min_offset`: `0.05, 0.10, 0.15, 0.25, 0.35` when the Otsu window is enabled;
- `otsu_max_offset`: `0.40, 0.50, 0.65` when the Otsu window is enabled;
- `min_accept_frac`: `0.50, 0.60, 0.65, 0.70, 0.80` for pool-context training.

Each candidate is evaluated by running the full classifier on every training dataset and averaging balanced accuracy over datasets with at least one non-hard-rejected row. If multiple candidates tie, the selected candidate is the one closest to the starting model in the searched threshold controls.

The learn report includes the search grid, suggested weights, selected model, search diagnostics, Otsu ablation diagnostics, feature-signal diagnostics, feature-drop diagnostics, leave-one-dataset-out feature-policy diagnostics, top candidates, best ties, and per-dataset confusion metrics.

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
