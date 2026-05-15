# Class-Average Quality Selection in SIMPLE

## Overview

`model_cavgs_rejection` is the SIMPLE class-average quality selection program. It evaluates the `cls2D` class averages in a project, applies conservative hard validity rejects, extracts a fixed scalar feature bank, normalizes those features within the dataset, and applies a named quality model to partition class averages into accepted and rejected sets.

The class-average quality backend lives under `src/main/cavg_quality` and is also usable as a library component from other SIMPLE workflows.

## Main Files

- `simple_cavg_quality_types.f90`: shared constants and derived types.
- `simple_cavg_quality_feats.f90`: feature definitions, feature extraction, hard-reject checks, and robust normalization.
- `simple_cavg_quality_stats.f90`: binary metrics and statistical helpers.
- `simple_cavg_quality_model.f90`: built-in model presets, model file I/O, scoring, clustering, thresholding, and promotion-code generation.
- `simple_cavg_quality_analysis.f90`: `apply` and `analyze` evaluation/reporting.
- `simple_cavg_quality_learn.f90`: analysis-table reader and model search.

The command entry point is `exec_model_cavgs_rejection` in `src/main/commanders/simple/simple_commanders_cavgs.f90`.

## Quality Modes

### `quality_mode=apply`

`apply` computes quality features, applies the selected model, writes selected/rejected class-average stacks, annotates the project, maps class selection into particles, optionally prunes, and writes the updated project.

Main outputs:

- `cavgs_quality_features.txt`
- `quality_selected_cavgs.mrc`
- `quality_rejected_cavgs.mrc`
- updated SIMPLE project

### `quality_mode=analyze`

`analyze` runs the selected model and treats the existing `cls2D` state as the manual reference. It writes a self-contained analysis file and diagnostic stacks without changing the project selection.

Main outputs:

- `cavgs_quality_analysis.txt`
- `quality_selected_cavgs.mrc`
- `quality_rejected_cavgs.mrc`

### `quality_mode=learn`

`learn` reads a table of `cavgs_quality_analysis.txt` files and searches for a model specification that reproduces the manual selections in those files. Hard-rejected rows remain in the report but are excluded from feature-weight estimation and model scoring.

Main outputs:

- learned model file, controlled by `fname=`
- `cavgs_quality_learn_report.txt`

### `quality_mode=promote`

`promote` reads a learned model file and writes a Fortran snippet for adding that model as a built-in preset.

Main output:

- promotion-code text file, controlled by `fname=`

Use a `.txt` suffix for the promotion snippet.

## Model Definition

A quality model contains:

- `name`
- `context`
- `weights(CAVG_QUALITY_NFEATS)`
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

`quality_model` selects a built-in preset. When no `quality_model` is provided, SIMPLE uses `chunk_default_v2`.

Built-in presets:

- `chunk_default_v2`: default chunk/stream model.
- `pool_default_v1`: recall-preserving model for larger pooled/batch datasets.

The active chunk model uses `boundary_margin=0.15`, `min_score_separation=0.15`, low-separation Otsu enabled, Otsu-window threshold replacement enabled with offsets `0.35` to `0.50`, cluster rescue disabled, and minimum accepted fraction disabled.

Model files are parsed as explicit key-value specifications. Unknown keys are rejected.

## Feature Bank

The fixed feature bank has 9 scalar features:

| Index | Feature | Direction | Meaning |
| --- | --- | --- | --- |
| 1 | `log_pop` | higher is better | Log class population. |
| 2 | `neg_log_res` | higher is better | Negative log class resolution estimate. |
| 3 | `centered` | higher is better | Negative normalized centroid displacement of segmented signal. |
| 4 | `log_locvar_fg` | higher is better | Log local variance in the foreground Otsu mask. |
| 5 | `log_locvar_bg` | higher is better | Log local variance in the complementary background mask. |
| 6 | `corr_frc_proxy` | higher is better | Stored class correlation or FRC-like score when present. |
| 7 | `log_center_edge_snr` | higher is better | Log central variance relative to edge variance. |
| 8 | `cc_area_frac` | higher is better | Main connected-component area divided by expected circular mask area. |
| 9 | `presence` | higher is better | Central signal excess relative to border background noise. |

The model uses normalized `z_*` feature values. Raw values are written for diagnostics.

## Hard Rejects

Hard rejects are validity gates applied before model fitting and clustering. A class average is hard rejected when any of these conditions hold:

- class population is non-positive;
- stored resolution is worse than `CAVG_RES_HARD_REJECT_A`;
- image values contain invalid pixels;
- foreground segmentation produces no valid component;
- the segmented foreground centroid lies outside the mask radius;
- the largest foreground component has more than `MASK_HARD_OUTSIDE_PIXELS` pixels outside the mask disc;
- foreground and background local variance are both non-positive.

Hard-rejected classes receive rejected state directly. The learned model operates only on non-hard-rejected class averages.

## Normalization

For each feature, non-hard-rejected rows are normalized with robust within-dataset statistics. Values are clipped by `CLIP_Z`. Hard-rejected rows are retained in the output tables but are excluded from robust-statistic estimation and training metrics.

## Classification

For non-hard-rejected class averages:

1. The model computes a linear score from normalized features and model weights.
2. Pairwise distances are computed in normalized feature space using the nonzero-weight features.
3. The distance matrix is robustly normalized.
4. Two-cluster k-medoids clustering is applied.
5. The cluster with higher mean model score is labeled as the good cluster.
6. A raw threshold is placed halfway between the good-cluster and bad-cluster mean scores.
7. Threshold controls determine the final score threshold:
   - `boundary_margin` shifts the cluster-derived threshold;
   - `use_lowsep_otsu` enables Otsu thresholding when cluster separation is low;
   - `use_otsu_window` enables Otsu replacement inside the configured offset window;
   - `use_cluster_rescue` accepts good-cluster members within `cluster_rescue_margin`;
   - `enforce_min_accept_frac` accepts the top-scoring classes required by `min_accept_frac`.

If there are fewer than four trainable class averages, the distance matrix is degenerate, clustering fails, or cluster separation is too low without an acceptable low-separation Otsu threshold, the model accepts all non-hard-rejected class averages as a single cluster.

## Learning

`quality_mode=learn` reads analysis files generated by the current feature bank. Each row supplies normalized features, the hard-reject flag, and the manual class state.
The analysis table must include the current header row; learn mode maps columns by name.

The learner derives feature weights from the training data:

```text
weight(feature) = max(0, pooled_auc(feature, manual_state) - 0.5)
```

Weights are normalized to sum to one. When all feature weights are zero, uniform weights are used.

The learner searches:

- `min_score_separation`
- `boundary_margin`
- `use_lowsep_otsu`
- `use_otsu_window`
- `otsu_min_offset`
- `otsu_max_offset`
- `min_accept_frac` for pool-context models

Each candidate is evaluated with the full classifier on every training dataset. Scoring uses macro balanced accuracy over non-hard-rejected rows.

The learn report records:

- search grid values;
- suggested weights;
- search diagnostics;
- Otsu ablation diagnostics;
- per-feature signal diagnostics;
- per-feature drop diagnostics;
- top candidate rows;
- best-tie rows;
- per-dataset confusion metrics.

## Example Commands

Analyze a manually selected project:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=analyze \
  projfile=my_project.simple \
  mkdir=yes
```

Apply the default chunk model:

```bash
simple_exec prg=model_cavgs_rejection \
  quality_mode=apply \
  projfile=my_project.simple \
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
