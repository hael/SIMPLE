# Class-Average Quality Selection in SIMPLE

## Overview

`model_cavgs_rejection` is the SIMPLE class-average quality selection program and the validation harness for the class-average quality library. The program evaluates every 2D class average, applies hard validity rejects for non-negotiable failures, measures a fixed inventory of scalar quality features, normalizes those features within the dataset, and applies an instantiable quality model that partitions the remaining class averages into accepted and rejected sets.

The current implementation is intentionally not a black-box classifier. A model is a named `linear_boundary` specification with feature weights, threshold controls, policy flags, and a context label such as `chunk` or `pool`. Learned models can be used directly from a model file, and a separate promotion mode can emit the Fortran code needed to add a validated learned model as a built-in library option.

The current policy is scalar-feature vector learning. Pairwise histogram/spectrum matrices and relational clustering-style quality evidence are not part of the default model path. Clustering is used only as a thresholding helper in normalized scalar-feature space; the quality evidence itself comes from interpretable per-class scalar metrics and their learned weights.

There is one user-facing selector: `quality_model`. If no `quality_model` is provided, SIMPLE initializes `CAVG_QUALITY_MODEL_CHUNK_DEFAULT`, currently `chunk_default_v2`. There is no separate `rejection_type` parameter and no user-visible `default` model name; chunk/pool behavior is encoded by the selected model itself.

## Main Files

The current implementation lives under `src/main/cavg_quality`:

- `simple_cavg_quality_types.f90`: shared derived types and constants.
- `simple_cavg_quality_feats.f90`: feature inventory and feature extraction.
- `simple_cavg_quality_stats.f90`: binary metrics and statistical helpers.
- `simple_cavg_quality_model.f90`: instantiable decision model, model I/O, model promotion-code generation.
- `simple_cavg_quality_analysis.f90`: apply/analyze evaluation and reporting.
- `simple_cavg_quality_learn.f90`: training-table reader and model search.

The command entry point is `exec_model_cavgs_rejection` in `src/main/commanders/simple/simple_commanders_cavgs.f90`.

Other SIMPLE workflows should treat `simple_cavg_quality` as a library backend. The intended streaming integration is for `simple_microchunked2D` to use an internal logical flag, for example `USE_CAVG_QUALITY_BACKEND`, to route class-average rejection through the quality backend. When enabled, that backend should instantiate `cavg_quality_model`, initialize `CAVG_QUALITY_MODEL_CHUNK_DEFAULT`, call `evaluate_cavg_quality`, and map the resulting states and scores through the standard project fields for class-average selection.

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

`quality_model` is the selector for built-in presets. When no explicit `quality_model` is provided, SIMPLE uses `chunk_default_v2`; pool-style behavior is selected by requesting a pool model such as `pool_default_v1` or by passing a model file through `infile=`.

The current built-in presets are:

- `chunk_default_v2`
- `chunk_default_v1`
- `pool_default_v1`

`chunk_default_v2` is the default chunk/stream operating point. The current promoted feature policy is `full_geom_pruned`, which leaves `mask_inside`, `single_component`, and `cc_area_frac` at zero weight because foreground geometry is handled by the hard validity gate or is not stable enough as soft model evidence. Resolution remains active model evidence through the nonzero-weight `neg_log_res` feature.

The current chunk model uses `boundary_margin=0.05`, `min_score_separation=0.15`, low-separation Otsu enabled, Otsu-window threshold replacement disabled, cluster rescue disabled, and minimum accepted fraction disabled. `chunk_default_v1` is an alternate chunk preset. `pool_default_v1` is the more recall-preserving operating point for larger pooled or batch sets and enables minimum-accepted-fraction behavior.

## Feature Space

The feature inventory is maintained in `simple_cavg_quality_feats.f90`. Feature extraction produces raw per-class values and a normalized feature matrix. The model uses the normalized `z_*` features, not the raw feature values.

The current feature set has 16 scalar features:

- `log_pop`
- `neg_log_res`
- `mask_inside`
- `centered`
- `log_locvar_fg`
- `log_locvar_bg`
- `single_component`
- `corr_frc_proxy`
- `log_center_edge_snr`
- `cc_area_frac`
- `interior_curvature`
- `presence`
- `log_detail_bg_snr`
- `log_detail_signal_ratio`
- `detail_coverage`
- `detail_edge_density`

The current model path is scalar-only. Pairwise histogram/spectrum distance matrices, relational quality matrices, histogram-neighborhood density, and global intensity-softness descriptors are not part of the active feature bank.

The internal-detail diagnostics measure organized in-mask molecular texture. They use one broad 20-6 A band-pass per class average and summarize detail relative to outside-mask background, detail relative to total in-mask signal, in-mask detail coverage, band-passed edge density, and an `interior_curvature` score. `interior_curvature` measures 4-neighbor Laplacian RMS in a shrunken central mask after softly downweighting pixels dominated by first-derivative edge response. This is meant to penalize smooth ice blobs whose signal is mostly a sharp contour, while preserving signal from true internal extrema/detail. These are scalar diagnostics rather than pairwise distances.

Resolution has a dual role. The continuous `neg_log_res` feature is part of the learned model and contributes to both the linear score and the feature-space distance whenever its learned weight is nonzero. Only catastrophic resolution failure is outside the learned boundary: a class with stored `cls2D` resolution worse than `40 A` is hard rejected before fitting. Foreground geometry is handled more conservatively; after Otsu segmentation and connected-component cleanup, a class is hard rejected if no valid foreground component remains, if any foreground component centroid lies outside the mask radius, or if the largest foreground component has more than 10 pixels outside the mask disc. The scalar geometry features remain in the analysis table as diagnostics, but catastrophic low-resolution or outside-mask failures are not left for the model to rescue.

The feature policy is encoded by zeroing excluded weights. Nonzero-weight features contribute to both the linear score and the k-medoids feature-space distance; hard rejections remain outside the learned policy and are applied before clustering. The current default keeps the backend scalar-only: pairwise histogram and spectrum matrices are not part of the model search or classification path.

Feature definitions should remain centralized in `simple_cavg_quality_feats.f90`, so future feature additions are visible in one inventory rather than being scattered through model code.

## Classification Method

For one dataset, the model classifies as follows:

1. Compute the linear quality score for every class:

   ```text
   score = sum(z_feature_i * weight_i)
   ```

2. Assign hard-rejected classes a low sentinel score and remove them from the fitted partition. Hard rejects include empty/dead classes, `cls2D res > 40 A`, and foreground geometry failures outside the expected mask.

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

Training is an explicit grid search over the existing interpretable model. It is not a black-box optimizer. The model selected by `quality_model` or `infile=` is the base model for learning; if neither is supplied, learn mode starts from `chunk_default_v2`.

Input is a file table of `cavgs_quality_analysis.txt` files. For each analysis file, the learner reads:

- the normalized `z_*` feature columns by name;
- `manual_state`, where `manual_state > 0` is the reference accepted class;
- `hard_reject`. When raw diagnostic columns are present, learn mode can also infer the current hard resolution and geometry vetoes directly from those columns.

Hard-rejected rows are not part of model fitting. They are excluded from suggested feature-weight AUCs, candidate-model scoring, feature-drop diagnostics, and leave-one-dataset-out feature-group screens. They remain in the report as hard-rejection diagnostics, especially `hard_rejected_manual_good`, because those cases indicate that the upstream validity gate may need attention rather than that the learned model should move its boundary.

The learner first derives suggested feature weights. For each feature, all non-hard-rejected training rows across all datasets are pooled, feature AUC is computed against the manual labels, and the suggested feature weight is:

```text
max(0, AUC - 0.5)
```

The suggested weights are normalized to sum to one. Features with AUC below chance are not inverted; they receive zero weight. If all suggested weights are zero, the feature-policy candidates use uniform weights over their active features.

The learner then searches a fixed grid. Candidate weights are not blended with the base model. For each feature-policy candidate, the learner starts from the data-derived suggested weights, zeroes the excluded features, and renormalizes the active weights:

```text
candidate_weights = normalize(policy_mask * suggested_weights)
```

If every active feature has chance-level AUC and the masked weight sum is zero, the candidate receives uniform weights over the active policy features. This is a neutral fallback, not base-model weight reuse.

The current grid searches:

- feature-policy candidate;
- minimum score separation;
- boundary margin;
- low-separation Otsu thresholding on/off;
- Otsu-window thresholding on/off;
- Otsu-window minimum and maximum offsets when Otsu-window thresholding is enabled;
- pool minimum accepted fraction, only when the base model context is `pool`.

The feature-policy candidates are hand-curated masks over the scalar feature bank. They include base scalar subsets, full-feature variants, and geometry-pruned variants such as `full_geom_pruned`.

Each candidate is evaluated by running the full model classifier on every training dataset, then scoring only the non-hard-rejected rows. The score is macro balanced accuracy: balanced accuracy is computed per dataset and then averaged across datasets with trainable rows, so each informative dataset contributes equally.

The best candidate becomes the learned model and is written to the requested model file. If multiple candidates tie exactly, the learner chooses the tied threshold-policy variant closest to the base model's stable interior settings; this prevents grid-order artifacts such as always selecting the lowest `min_score_separation` when the whole searched range performs identically. The learn report also records the full search grid, the selected feature policy, the top candidates, all best-score ties, suggested weights, and per-dataset confusion metrics. The per-dataset rows include both `n_classes` and `n_trainable`; the reported confusion metrics are computed over `n_trainable`, while hard-rejected classes are summarized separately.

The learn report includes a `search_diagnostic` section. These records flag two classes of follow-up:

- searched parameters whose winning value is on the edge of the current grid, such as `boundary_margin`, `min_score_separation`, or pool `min_accept_frac`;
- model policy controls that are still inherited from the base preset rather than searched, such as cluster rescue and chunk/pool minimum-accept enforcement.

Warnings in this section do not mean the learned model is invalid. They mean the training set is touching one of the current search boundaries or inherited policy assumptions, so the next validation round should decide whether to broaden the grid or promote that policy control into the searched model space.

The report also includes feature-screen diagnostics for deciding what to try next:

- `otsu_ablation` compares the learned model against the same learned model with both Otsu threshold routes disabled, dataset by dataset.
- `feature_signal` gives pooled AUC, mean/min/max per-dataset AUC, inversion count, base weight, suggested weight, and learned weight for every scalar feature, using only non-hard-rejected rows.
- `feature_drop` reruns the learned model with one scalar feature weight set to zero; negative `delta_vs_learned` means the feature helped the learned model on the non-hard-rejected training set.
- `feature_group_lodo` performs leave-one-dataset-out scalar scoring for fixed feature groups. These rows are a separability screen, not a promoted model: weights are learned from non-hard-rejected rows in all other datasets, then the non-hard-rejected rows of the held-out dataset are scored by AUC and by an optimistic best-threshold balanced accuracy.

## Basic Instruction Manual

### 1. Generate Analysis Files

Run `analyze` on each dataset with a trusted manual selection already encoded in `cls2D` state. For chunk/stream data, `quality_model` may be omitted because `chunk_default_v2` is the default; specifying it explicitly can be useful for reproducible command logs.

```bash
simple_exec prg=model_cavgs_rejection quality_mode=analyze \
  projfile=my_project.simple \
  mskdiam=180 \
  quality_model=chunk_default_v2
```

For pool-like datasets, use:

```bash
simple_exec prg=model_cavgs_rejection quality_mode=analyze \
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

Use separate file tables for chunk and pool training unless deliberately testing cross-context behavior. The model context is carried by the model, not by a separate rejection-type parameter.

### 3. Learn a New Model

For chunk:

```bash
simple_exec prg=model_cavgs_rejection quality_mode=learn \
  quality_model=chunk_default_v2 \
  filetab=cavgs_quality_chunk_analyses.txt \
  fname=cavgs_quality_model_chunk_learned.txt
```

For pool:

```bash
simple_exec prg=model_cavgs_rejection quality_mode=learn \
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
simple_exec prg=model_cavgs_rejection quality_mode=analyze \
  projfile=my_project.simple \
  mskdiam=180 \
  infile=cavgs_quality_model_chunk_learned.txt
```

Compare the new analysis file and accepted/rejected stacks against the manual selection and the selected built-in model.

### 5. Apply the Learned Model

After validation, apply the learned model:

```bash
simple_exec prg=model_cavgs_rejection quality_mode=apply \
  projfile=my_project.simple \
  mskdiam=180 \
  infile=cavgs_quality_model_chunk_learned.txt
```

This updates the project selection. Use `analyze` first when validating a model.

### 6. Generate Built-In Promotion Code

When a learned model is good enough to become a named library option, run:

```bash
simple_exec prg=model_cavgs_rejection quality_mode=promote \
  infile=cavgs_quality_model_chunk_learned.txt \
  fname=cavgs_quality_model_chunk_builtin_code.txt
```

The generated text file contains:

- a `CAVG_QUALITY_MODEL_*` constant;
- the required `BUILTIN_MODEL_NAMES` addition;
- the `builtin_spec` case;
- the full `*_model_spec()` function;
- notes about updating the UI/help option lists.

Review the code, paste it into `simple_cavg_quality_model.f90`, update `BUILTIN_MODEL_NAMES` and the user-visible `quality_model` option lists, rebuild, and then test the promoted model by name:

```bash
simple_exec prg=model_cavgs_rejection quality_mode=analyze \
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

For streaming integration, validate the promoted chunk model in two passes:

- first through `model_cavgs_rejection quality_mode=analyze`, where trusted manual selections can be compared directly;
- then through the microchunk backend guarded by the internal `USE_CAVG_QUALITY_BACKEND` flag, where the important checks are unchanged sentinel behavior, identical project-state propagation, and acceptable selected/rejected class-average stacks.

## Current Limitations and Future Directions

The current learner searches only a limited set of model parameters. It learns feature weights from the training set ab initio, but still inherits behavior flags such as Otsu use and cluster rescue from the base model. That is deliberate for now. The learn report flags when those inherited controls are active on datasets that still have false positives or false negatives, which is the cue to consider searching those flags explicitly.

The accept-all behavior is currently implicit through the single-cluster fallback and `min_score_separation`. A future model version may make accept-all behavior an explicit learned operating mode.

Feature weights are currently suggested by independent AUC scoring on non-hard-rejected rows. A linear SVM or regularized logistic regression could be added later as an alternative way to propose multivariate weights, while still keeping the explicit SIMPLE decision model around the learned coefficients.
