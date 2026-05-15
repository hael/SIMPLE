# Class-Average Quality Selection in SIMPLE

## Overview

`cluster_cavgs_quality` is the SIMPLE class-average quality selection program. It replaces a field of scalar rejection thresholds with an explicit, interpretable multivariate model. The program evaluates every 2D class average, applies hard validity rejects for non-negotiable failures, measures a fixed inventory of quality features, normalizes those features within the dataset, and applies an instantiable quality model that partitions the remaining class averages into accepted and rejected sets.

The current implementation is intentionally not a black-box classifier. A model is a named `linear_boundary` specification with feature weights, threshold controls, and policy flags for chunk-like or pool-like behavior. Learned models can be used directly from a model file, and a separate promotion mode can emit the Fortran code needed to add a validated learned model as a built-in library option.

## Main Files

The current implementation lives under `src/main/cavg_quality`:

- `simple_cavg_quality_types.f90`: shared derived types and constants.
- `simple_cavg_quality_feats.f90`: feature inventory and feature extraction.
- `simple_cavg_quality_stats.f90`: binary metrics and statistical helpers.
- `simple_cavg_quality_model.f90`: instantiable decision model, model I/O, model promotion-code generation.
- `simple_cavg_quality_analysis.f90`: apply/analyze evaluation and reporting.
- `simple_cavg_quality_learn.f90`: training-table reader and model search.

The command entry point is `exec_cluster_cavgs_quality` in `src/main/commanders/simple/simple_commanders_cavgs.f90`.

## Quality Modes

### `quality_mode=apply`

`apply` is the production selection mode. It reads a SIMPLE project, computes quality features for `cls2D`, applies the selected model, writes accepted/rejected class-average stacks, annotates the project, maps the class selection into particles, optionally prunes, and writes the updated project.

Important outputs:

- `cavgs_quality_features.txt`
- `quality_selected_cavgs.mrc`
- `quality_rejected_cavgs.mrc`
- updated SIMPLE project

### `quality_mode=analyze`

`analyze` runs the same model but treats the current `cls2D` state as the manual reference. It writes a self-contained analysis table and diagnostics, but leaves the project selection unchanged.

Important outputs:

- `cavgs_quality_analysis.txt`
- `quality_selected_cavgs.mrc`
- `quality_rejected_cavgs.mrc`

The analysis file contains per-class model output, raw features, normalized `z_*` features, manual states, model parameters, feature summaries, and threshold scans.

### `quality_mode=learn`

`learn` reads a file table of `cavgs_quality_analysis.txt` files and searches for a new model specification that best reproduces the manual selections in those files.

Important outputs:

- learned model file, controlled by `fname=`
- `cavgs_quality_learn_report.txt`

The learned model file can be passed directly to `apply` or `analyze` through `infile=`.

### `quality_mode=promote`

`promote` reads a learned model file and emits a text snippet with the Fortran code needed to make that model a built-in preset in `simple_cavg_quality_model.f90`.

Important output:

- promotion-code text file, controlled by `fname=`

Use a `.txt` suffix for the promotion snippet. `fname=...f90` is treated as a file-format-bearing SIMPLE argument and is rejected by the generic file-format checker.

## Model Definition

The current model family is `linear_boundary`. A model contains:

- `name`
- `family`
- `context`
- `feature_policy`
- `rejection_type`
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

`quality_model` is the primary selector for built-in presets. `rejection_type` is retained as a legacy hint for choosing a default model when no explicit model is requested.

The current built-in presets are:

- `chunk_default_v2`
- `chunk_default_v1`
- `pool_default_v1`

`chunk_default_v2` is the default chunk/stream operating point promoted from the representative batch7chunk learning cycle. Its current feature bank excludes the former scalar histogram-neighborhood feature, pairwise histogram/spectrum matrices, and unstable spectrum/ice/foreground-background ratio diagnostics. The promoted feature policy is `base12_pruned_plus_histvar`, which leaves `mask_inside`, `single_component`, `cc_area_frac`, `presence`, and `log_contrast` at zero weight while keeping masked histogram variance available to the scalar score and clustering distance. `chunk_default_v1` is retained as the legacy chunk preset. `pool` is the more recall-preserving operating point for larger pooled or batch sets.

## Feature Space

The feature inventory is maintained in `simple_cavg_quality_feats.f90`. Feature extraction produces raw per-class values and a normalized feature matrix. The model uses the normalized `z_*` features, not the raw feature values.

The current feature set includes population/resolution support, foreground geometry, local signal variance, class statistics, center-edge signal, bounded histogram entropy, two connected-component shape diagnostics, central presence, full-image contrast, and masked histogram variance. A former scalar histogram-neighborhood feature, pairwise histogram/spectrum distance matrices, `spectrum_dynrange`, `neg_ice_score`, and `log_fg_bg_locvar_ratio` were removed because representative chunk training showed they were weak, unstable, or redundant for the default scalar model.

Foreground geometry also contributes hard validity rejects before the learned model is fit. After Otsu segmentation and connected-component cleanup, a class is hard rejected if no valid foreground component remains, if any foreground component centroid lies outside the mask radius, or if the largest foreground component has more than 10 pixels outside the mask disc. The scalar geometry features remain in the analysis table as diagnostics, but catastrophic outside-mask density is not left for the model to rescue.

The feature policy is encoded by zeroing excluded weights. Nonzero-weight features contribute to both the linear score and the k-medoids feature-space distance; hard rejections remain outside the learned policy and are applied before clustering.

Feature definitions should remain centralized in `simple_cavg_quality_feats.f90`, so future feature additions are visible in one inventory rather than being scattered through model code.

## Classification Method

For one dataset, the model classifies as follows:

1. Compute the linear quality score for every class:

   ```text
   score = sum(z_feature_i * weight_i)
   ```

2. Assign hard-rejected classes a low sentinel score and remove them from the fitted partition. Hard rejects include empty/dead classes and foreground geometry failures outside the expected mask.

3. Build a Euclidean distance matrix between non-hard-rejected normalized feature vectors, using only features whose model weight is nonzero.

4. Normalize that feature distance matrix.

5. Run two-cluster k-medoids on the distance matrix.

6. Identify the good cluster as the cluster with the higher mean linear score.

7. Set the raw score threshold halfway between the good-cluster and bad-cluster mean scores.

8. Apply model policy:

    - If cluster score separation is too low, use low-separation Otsu only when enabled and sufficiently separated.
    - Otherwise apply `boundary_margin`.
    - Optionally replace the threshold with Otsu when `use_otsu_window` permits it.
    - Optionally rescue good-cluster classes close to the threshold.
    - Optionally enforce a minimum accepted fraction for pool-like models.

The final state is `state > 0` for accepted and `state <= 0` for rejected.

## Accept-All Fallback

Some datasets should not be split. The current mechanism is an implicit single-cluster fallback. The model accepts all non-hard-rejected classes when the data do not support a reliable partition.

This happens when:

- fewer than four non-hard-rejected classes are available;
- the feature distance matrix is flat after normalization;
- k-medoids fails to produce two valid nonempty clusters;
- the score separation is below `min_score_separation` and no acceptable low-separation Otsu threshold is available.

Hard rejects are still final. A class with `hard_reject=.true.` remains rejected even when the fitted model falls back to accepting the non-hard-rejected classes.

## Training Method

Training is an explicit grid search over the existing interpretable model. It is not a black-box optimizer.

Input is a file table of `cavgs_quality_analysis.txt` files. For each analysis file, the learner reads:

- the normalized `z_*` feature columns by name;
- `manual_state`, where `manual_state > 0` is the reference accepted class;
- `hard_reject`.

The learner first derives suggested feature weights. For each feature, all training rows across all datasets are pooled, feature AUC is computed against the manual labels, and the suggested feature weight is:

```text
max(0, AUC - 0.5)
```

The suggested weights are normalized to sum to one. Features with AUC below chance are not inverted; they receive zero weight. If all suggested weights are zero, the learner falls back to the default weights.

The learner then searches a fixed grid. Candidate weights are interpolated between the base model weights and the suggested weights:

```text
candidate_weights = (1 - alpha) * base_weights + alpha * suggested_weights
```

The current grid searches:

- feature-policy candidate;
- feature-weight interpolation alpha;
- minimum score separation;
- boundary margin;
- pool minimum accepted fraction, only for pool models.

Each candidate is evaluated by running the full model classifier on every training dataset. The score is macro balanced accuracy: balanced accuracy is computed per dataset and then averaged across datasets, so each dataset contributes equally.

The best candidate becomes the learned model and is written to the requested model file. The learn report also records the full search grid, the selected feature policy, the top candidates, all best-score ties, suggested weights, and per-dataset confusion metrics.

The learn report includes a `search_diagnostic` section. These records flag two classes of follow-up:

- searched parameters whose winning value is on the edge of the current grid, such as `boundary_margin`, `min_score_separation`, or pool `min_accept_frac`;
- model policy controls that are currently inherited from the base preset rather than searched, such as Otsu-window settings, low-separation Otsu, cluster rescue, and minimum-accept enforcement.

Warnings in this section do not mean the learned model is invalid. They mean the training set is touching one of the current search boundaries or inherited policy assumptions, so the next validation round should decide whether to broaden the grid or promote that policy control into the searched model space.

The report also includes feature-screen diagnostics for deciding what to try next:

- `feature_signal` gives pooled AUC, mean/min/max per-dataset AUC, inversion count, base weight, suggested weight, and learned weight for every scalar feature.
- `feature_drop` reruns the learned model with one scalar feature weight set to zero; negative `delta_vs_learned` means the feature helped the learned model on the training set.
- `feature_group_lodo` performs leave-one-dataset-out scalar scoring for fixed feature groups. These rows are a separability screen, not a promoted model: weights are learned from all other datasets, then the held-out dataset is scored by AUC and by an optimistic best-threshold balanced accuracy.

## Basic Instruction Manual

### 1. Generate Analysis Files

Run `analyze` on each dataset with a trusted manual selection already encoded in `cls2D` state.

```bash
simple_exec prg=cluster_cavgs_quality quality_mode=analyze \
  projfile=my_project.simple \
  mskdiam=180 \
  quality_model=chunk_default_v2
```

For pool-like datasets, use:

```bash
simple_exec prg=cluster_cavgs_quality quality_mode=analyze \
  projfile=my_project.simple \
  mskdiam=180 \
  quality_model=pool_default_v1
```

Keep the resulting `cavgs_quality_analysis.txt` files from each dataset. Rename or store them in separate run directories so they are not overwritten.

### 2. Build a Training File Table

Create a text file with one analysis file path per line:

```text
/path/to/dataset_001/cavgs_quality_analysis.txt
/path/to/dataset_002/cavgs_quality_analysis.txt
/path/to/dataset_003/cavgs_quality_analysis.txt
```

Use separate file tables for chunk and pool training unless deliberately testing cross-context behavior.

### 3. Learn a New Model

For chunk:

```bash
simple_exec prg=cluster_cavgs_quality quality_mode=learn \
  quality_model=chunk_default_v2 \
  filetab=cavgs_quality_chunk_analyses.txt \
  fname=cavgs_quality_model_chunk_learned.txt
```

For pool:

```bash
simple_exec prg=cluster_cavgs_quality quality_mode=learn \
  quality_model=pool_default_v1 \
  filetab=cavgs_quality_pool_analyses.txt \
  fname=cavgs_quality_model_pool_learned.txt
```

Inspect `cavgs_quality_learn_report.txt`. Pay attention to:

- `macro_balanced_accuracy`;
- `best_tie_count`;
- `top_candidate` rows;
- per-dataset false positives and false negatives;
- `hard_rejected_manual_good`;
- `feature_signal` rows with low minimum AUC or inverted datasets;
- `feature_drop` rows with large negative or positive `delta_vs_learned`;
- `feature_group_lodo` rows, especially mean/min AUC and the weakest held-out dataset.

Many tied models mean the current training set does not constrain the searched parameters strongly enough. That is not automatically bad, but it should be understood before promoting a model.

### 4. Test the Learned Model

Use the learned model file directly through `infile=`:

```bash
simple_exec prg=cluster_cavgs_quality quality_mode=analyze \
  projfile=my_project.simple \
  mskdiam=180 \
  infile=cavgs_quality_model_chunk_learned.txt
```

Compare the new analysis file and accepted/rejected stacks against the original manual selection and against the previous default model.

### 5. Apply the Learned Model

After validation, apply the learned model:

```bash
simple_exec prg=cluster_cavgs_quality quality_mode=apply \
  projfile=my_project.simple \
  mskdiam=180 \
  infile=cavgs_quality_model_chunk_learned.txt
```

This updates the project selection. Use `analyze` first when validating a model.

### 6. Generate Built-In Promotion Code

When a learned model is good enough to become a named library option, run:

```bash
simple_exec prg=cluster_cavgs_quality quality_mode=promote \
  infile=cavgs_quality_model_chunk_learned.txt \
  fname=cavgs_quality_model_chunk_builtin_code.txt
```

The generated text file contains:

- a `CAVG_QUALITY_MODEL_*` constant;
- the required `BUILTIN_MODEL_NAMES` addition;
- the `builtin_spec` case;
- the full `*_model_spec()` function;
- notes about updating the UI/help option lists.

Review the code, paste it into `simple_cavg_quality_model.f90`, update the user-visible `quality_model` option lists, rebuild, and then test the promoted model by name:

```bash
simple_exec prg=cluster_cavgs_quality quality_mode=analyze \
  projfile=my_project.simple \
  mskdiam=180 \
  quality_model=chunk_default_v2
```

## Recommended Validation Practice

Start with a small number of high-confidence datasets. Inspect false positives and false negatives visually using the selected/rejected MRC stacks. Then add harder datasets, especially cases where the correct behavior is to reject almost nothing.

For every candidate learned model, check:

- whether hard rejects remove any manual-good classes;
- whether chunk and pool data need different parameters;
- whether any newly added scalar feature is stable across training sets;
- whether the best model is unique or one of many ties;
- whether balanced accuracy hides unacceptable false negatives;
- whether the model behaves sensibly on a dataset not used for learning.

## Current Limitations and Future Directions

The current learner searches only a limited set of model parameters. It inherits behavior flags such as Otsu use and cluster rescue from the base model. That is deliberate for now. The learn report flags when those inherited controls are active on datasets that still have false positives or false negatives, which is the cue to consider searching those flags explicitly.

The accept-all behavior is currently implicit through the single-cluster fallback and `min_score_separation`. A future model version may make accept-all behavior an explicit learned operating mode.

Feature weights are currently suggested by independent AUC scoring and then blended with the base weights. A linear SVM or regularized logistic regression could be added later as an alternative way to propose multivariate weights, while still keeping the explicit SIMPLE decision model around the learned coefficients.
