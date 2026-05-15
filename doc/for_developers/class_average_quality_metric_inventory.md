# Class-Average Quality Metric Inventory

**Date:** May 14, 2026
**Status:** Developer inventory / future feature planning
**Scope:** Existing SIMPLE routines that can support scalar features for `cluster_cavgs_quality`.

## Purpose

This note inventories image-processing, class-average, and clustering routines already present in SIMPLE that can be reused when extending the `cluster_cavgs_quality` feature bank. The goal is to keep quality selection interpretable and reproducible: every feature should have a clear source, a documented direction, and a diagnostic raw value in the analysis table.

The current implementation treats class-average quality as a feature-vector and model-selection problem:

1. Extract raw scalar descriptors for each 2D class average.
2. Apply only conservative hard rejections for pathological cases.
3. Normalize non-hard-rejected classes within the dataset by robust statistics.
4. Use an explicit quality model to score and cluster the class averages.
5. Preserve the raw features, normalized features, model parameters, and manual reference labels in `cavgs_quality_analysis.txt`.
6. Use `quality_mode=learn` to tune interpretable model parameters against trusted manual selections.

## Current Implementation

The refactored code lives under `src/main/cavg_quality`:

- `simple_cavg_quality_types.f90`: shared constants and derived types.
- `simple_cavg_quality_feats.f90`: feature definitions, feature extraction, robust normalization, and histogram distance extraction.
- `simple_cavg_quality_stats.f90`: binary metrics, AUC, robust summaries, and distance-matrix normalization helpers.
- `simple_cavg_quality_model.f90`: built-in model presets, model file I/O, linear scoring, mixed-distance clustering, threshold policy, and promotion-code generation.
- `simple_cavg_quality_analysis.f90`: `apply`/`analyze` evaluation and reporting.
- `simple_cavg_quality_learn.f90`: analysis-table reader and grid search over interpretable model parameters.

The command entry point is `exec_cluster_cavgs_quality` in `src/main/commanders/simple/simple_commanders_cavgs.f90`.

The current application modes are:

- `quality_mode=apply`: production selection, writes selected/rejected stacks and updates the project.
- `quality_mode=analyze`: runs the model against an existing manual selection and writes a self-contained analysis file.
- `quality_mode=learn`: reads a file table of analysis files and learns a model file.
- `quality_mode=promote`: converts a validated learned model file into a Fortran built-in preset snippet.

## Current Feature Bank

The active inventory is centralized in `simple_cavg_quality_feats.f90` through `FEATURE_DEFS` and `CAVG_QUALITY_NFEATS`. Every feature is emitted as both a raw value and a normalized `z_*` value in analysis output.

| Index | Feature | Source | Direction | Notes |
| --- | --- | --- | --- | --- |
| 1 | `log_pop` | `cls2D` `pop` | higher is better | Log population support. |
| 2 | `neg_log_res` | `cls2D` `res` | higher is better | Lower Angstrom resolution maps to a larger score; `res > 40 A` is hard rejected before modeling. |
| 3 | `mask_inside` | Otsu foreground CC vs. mask disc | higher is better | Negative outside-mask fraction of the main segmented component. Extreme failures are hard rejected before modeling. |
| 4 | `centered` | Otsu CC mass centers | higher is better | Negative normalized centroid displacement. Component centroids outside the mask are hard rejected before modeling. |
| 5 | `log_locvar_fg` | `image%loc_var_masked` | higher is better | Foreground local variance after low-pass filtering and Otsu masking. |
| 6 | `log_locvar_bg` | `image%loc_var_masked` | higher is better | Complementary background local variance. |
| 7 | `single_component` | `image_bin%find_ccs` | higher is better | Negative distance from exactly one valid component. |
| 8 | `corr_frc_proxy` | optional `cls2D` `corr` | higher is better | Uses stored correlation/FRC-like metadata when present. |
| 9 | `log_center_edge_snr` | `image%center_edge_snr` | higher is better | Central variance relative to edge variance. |
| 10 | `hist_entropy` | `histogram%entropy` | diagnostic | Bounded entropy from masked intensity histograms. |
| 11 | `cc_area_frac` | `image_bin%size_ccs` plus Otsu CC mask | diagnostic | Main connected-component area divided by the expected circular mask area. |
| 12 | `cc_diameter_norm` | `image_bin%diameter_cc` | diagnostic | Main connected-component diameter divided by the expected mask diameter. |
| 13 | `presence` | `image%presence` | higher is better | Central signal excess relative to border background noise. |
| 14 | `log_contrast` | `image%contrast` | diagnostic | Log full-image intensity standard deviation after edge background subtraction. |
| 15 | `log_hist_variance` | `histogram%variance` | diagnostic | Log masked intensity-histogram variance. |

Pairwise histogram and rotational-spectrum matrices are no longer part of the default model infrastructure. The active histogram features are scalar entropy and variance only. The earlier `hist_knn` scalar was removed because nearest-neighbor histogram density did not have a stable quality direction across datasets.

The representative chunk runs also removed `spectrum_dynrange`, `neg_ice_score`, and `log_fg_bg_locvar_ratio` from the feature bank. Those signals were either weak, redundant with retained scalar features, or inverted on some datasets.

Resolution and mask-geometry validity are deliberately outside the learned model. The feature extractor hard rejects class averages when `cls2D res > 40 A`, when the Otsu/connected-component pass finds no valid foreground component, when any foreground-component centroid lies outside the mask radius, or when more than 10 pixels of the largest component lie outside the mask disc. This mirrors the proven microchunk rejection concept without depending on the stream/microchunk modules.

## Current Model Controls

The built-in models are complete `linear_boundary` presets:

- `chunk_default_v2`: current default chunk/stream behavior, promoted from the representative batch8chunk learning cycle with the `all_features_no_mask_single` feature policy.
- `chunk_default_v1`: legacy chunk/stream behavior retained for comparison.
- `pool_default_v1`: recall-preserving behavior for larger pooled/batch datasets.

The important model controls are:

- `feature_policy`: learn-mode candidate name describing which scalar features are allowed into the score and clustering distance; the model file encodes this by zeroing excluded weights. The current chunk default zeros `mask_inside` and `single_component` because those soft geometry flags duplicate the hard connected-component rejection.
- `feature_weights`: nonnegative linear score coefficients, normalized to sum to one.
- `boundary_margin`: shifts the decision threshold; more negative rejects more aggressively, more positive preserves more borderline classes.
- `min_score_separation`: below this separation the model treats the partition as unreliable unless low-separation Otsu is enabled and acceptable.
- `otsu_min_offset` and `otsu_max_offset`: constrain when an Otsu score threshold may replace the cluster-derived threshold.
- `cluster_rescue_margin`: lets pool-like models rescue good-cluster members close to threshold.
- `min_accept_frac`: enforces a minimum accepted fraction for pool-like models.

The current learner searches feature-policy candidates, feature-weight interpolation, minimum score separation, boundary margin, and pool minimum accepted fraction. The candidate policies include an `all_features_no_mask_single` option that keeps the newer scalar diagnostics while excluding the two soft geometry flags most duplicated by hard geometry rejection. It deliberately inherits the higher-level policy flags from the base model.

## Reusable Existing Routines

### Microchunk Rejector

Source:

- `src/main/stream/simple_cluster2D_rejector.f90`

Relevant routines:

- `reject_pop`
- `reject_res`
- `reject_mask`
- `reject_local_variance`

These routines remain useful as historical, operationally tested rejection criteria. They are now better viewed as feature sources or sanity checks than as the final quality selector. The refactored quality app already mirrors the strongest pieces: population, resolution, mask geometry, connected components, and local variance.

### Older Non-Junk Gate

Source:

- `src/main/strategies/search/simple_strategy2D_utils.f90`

Relevant routine:

- `flag_non_junk_cavgs`

This path combines spectral dynamic range, density inside/outside a mask, center offset, and minimum population. It is useful as a source of candidate features and for understanding prior behavior, but its hard gate should not be copied wholesale into the new model.

### Cluster-Average Similarity Matrices

Source:

- `src/main/strategies/search/simple_strategy2D_utils.f90`

Relevant routines and fields:

- `calc_sigstats_dmats`
- `calc_cc_and_res_dmats`
- `calc_cluster_cavgs_dmat`
- `align_and_score_cavg_clusters`
- `clust_info%homogeneity`
- `clust_info%resscore`
- `clust_info%clustscore`
- `clust_info%jointscore`
- `clust_info%icescore`

These routines are most valuable when quality is better judged by consistency with other class averages than by single-image statistics. They are also more expensive because some paths perform class matching or cluster alignment. They are deferred experiments, not part of the scalar default model.

### FRC and Resolution Infrastructure

Sources:

- `src/main/class/simple_class_frcs.f90`
- `src/main/class/simple_classaverager_restore.f90`
- `src/main/class/simple_classaverager_core.f90`
- `src/utils/filter/simple_estimate_ssnr.f90`
- `src/utils/filter/simple_opt_filter.f90`

Relevant routines:

- `class_frcs%frc_getter`
- `class_frcs%avg_frc_getter`
- `class_frcs%estimate_res`
- `class_frcs%estimate_find_for_eoavg`
- `class_frcs%estimate_lp_for_align`
- `fsc2ssnr`
- `fsc2optlp`
- `estimate_lplim`
- `estimate_lplims2D`

Stored `cls2D` metadata should be preferred when available. Recomputing even/odd FRC-like summaries during quality selection would be more invasive and should be justified by validation, but FRC-band summaries are scientifically plausible future features.

### Shape Ranking

Source:

- `src/main/commanders/simple/simple_commanders_cavgs.f90`

Relevant routine:

- `exec_shape_rank_cavgs`

This path uses `automask2D`, integrated masked intensity, object diameter, and center shifts to rank class averages. Several of those outputs are attractive quality features if they can be computed without making quality selection too slow or too dependent on automask parameters.

### Image-Level Primitives

Sources:

- `src/main/image/simple_image.f90`
- `src/main/image/simple_image_calc.f90`
- `src/main/image/simple_image_msk.f90`
- `src/main/image/simple_image_bin.f90`
- `src/main/simple_pspecs.f90`
- `src/utils/math/simple_histogram.f90`
- `src/utils/math/simple_stat.f90`

Useful primitives include:

- `image%loc_var_masked`, `image%variance`, `image%skew`, `image%kurt`, `image%GLCM`
- `image%spectrum`, `image%power_spectrum`, `image%guinier`, `image%guinier_bfac`, `image%frc_pspec`
- `image%fcomps_below_noise_power_stats`
- `image%center_edge_snr`, `image%presence`, `image%contrast`, `image%contains_nans`
- `image%corr`, `image%corr_shifted`, `image%real_corr`, `image%phase_corr`, `image%radial_cc`, `image%sqeuclid`
- `image%calc_ice_score`
- `image_bin%find_ccs`, `image_bin%size_ccs`, `image_bin%diameter_cc`, `image_bin%masscen_cc`
- `automask2D`, `density_inoutside_mask`, `image%density_inoutside`
- `histogram%variance`, `histogram%skew`, `histogram%entropy`, `histogram%TVD`, `histogram%JSD`, `histogram%HD`

## Feasible Future Feature Extensions

The best future additions are cheap, interpretable, and likely to generalize across datasets. They should be added one or two at a time, then validated through `quality_mode=analyze` and `quality_mode=learn` before changing built-in weights.

### Low-Risk Scalar Features

These are feasible first additions because they reuse existing per-image routines and do not require pairwise alignment.

| Candidate | Existing source | Why it may help | Caution |
| --- | --- | --- | --- |
| `contrast` variants | `image%contrast` | Detects extremely weak or washed-out averages. | `log_contrast` is now active as a diagnostic; consider variants only if validation shows a more stable transform is needed. |
| `presence` variants | `image%presence` | Direct particle-presence proxy. | `presence` is now active; inspect scale and sign across stream and pool datasets before assigning much model weight. |
| `hist_variance` variants | `histogram%variance` | Complements local variance with global intensity spread. | `log_hist_variance` is now active and correlated with contrast and local variance. |
| `fg_bg_locvar_ratio` variants | `image%loc_var_masked` | Tests whether foreground texture exceeds background texture. | Removed from the default bank after inverted chunk behavior; only revisit with new evidence. |
| `nan_or_bad_pixel_flag` | `image%contains_nans`, min/max checks | Hard diagnostic for corrupt inputs. | Should probably remain a hard reject or warning, not a learned soft feature. |

### Medium-Risk Spectral Features

These are plausible and cheap enough, but our experience with `spectrum_dynrange` says they should start as diagnostics.

| Candidate | Existing source | Why it may help | Caution |
| --- | --- | --- | --- |
| `bandpower_mid_low_ratio` | `image%spectrum` or `image%power_spectrum` | Separates real mid-frequency structure from low-frequency blobs. | Ring artifacts can mimic high quality. |
| `bandpower_high_mid_ratio` | `image%spectrum` | Catches over-sharpened/overfitted averages. | Needs ice/ring-aware calibration. |
| `guinier_slope` or `bfac_proxy` | `image%guinier`, `image%guinier_bfac` | Captures spectral falloff shape. | Sensitive to preprocessing and box/mask choices. |
| `below_noise_fraction` | `image%fcomps_below_noise_power_stats` | Identifies weak averages with little usable Fourier signal. | Needs a stable noise definition for class averages. |
| `frc_pspec_summary` | `image%frc_pspec` | May approximate consistency without full even/odd class FRC infrastructure. | Must verify the input assumptions for class-average images. |

Do not reintroduce spectrum-derived defaults casually. Representative chunk runs did not justify `spectrum_dynrange`, `neg_ice_score`, or rotational-spectrum pairwise distances.

### Deferred High-Cost Relational Features

These are attractive for difficult cases, especially when junk is statistically similar to good classes, but they add quadratic cost and sometimes alignment cost. They should not be reintroduced as pairwise matrices in the default model.

| Candidate | Existing source | Why it may help | Caution |
| --- | --- | --- | --- |
| `sig_knn` | `calc_sigstats_dmats` | Could summarize signal-statistics neighborhoods as scalar support. | Scalar neighborhood density can invert across datasets. |
| `cc_knn` | `calc_cc_and_res_dmats` | Good classes should correlate with nearby good classes after in-plane search. | Expensive for stream chunks; consider optional analyze/learn experiments. |
| `res_knn` | `calc_cc_and_res_dmats` | Pairwise FRC/resolution consistency may catch overfit junk. | Depends on robust pairwise matching. |
| `mixed_support_knn` | `calc_cluster_cavgs_dmat` | Reuses Joe's signal/correlation/resolution blend as scalar support. | Do not wire the matrix itself into the default model. |
| `cluster_jointscore` | `align_and_score_cavg_clusters` | Converts cluster homogeneity/resolution into interpretable scores. | More invasive and too heavy for a default streaming path until timed. |

### FRC-Derived Features

These are scientifically appealing if the relevant metadata is already available from class-average restoration.

| Candidate | Existing source | Why it may help | Caution |
| --- | --- | --- | --- |
| `frc_band_mean` | `class_frcs%frc_getter` | Uses more of the FRC curve than one resolution cutoff. | Avoid recomputing FRCs in quality selection if they are not already persisted. |
| `frc_drop_index` | `class_frcs%frc_getter` | Captures where class consistency collapses. | Resolution estimates already encode part of this. |
| `ssnr_band_mean` | `fsc2ssnr` | Converts FRC/FSC into a signal-to-noise style summary. | Needs robust handling near invalid correlation values. |
| `optlp_margin` | `estimate_lplims2D`, `fsc2optlp` | Measures whether usable signal extends beyond the alignment low-pass. | More useful in batch/pool contexts than tiny stream chunks. |

### Shape and Automask Features

These may help when class averages contain strong non-particle objects.

| Candidate | Existing source | Why it may help | Caution |
| --- | --- | --- | --- |
| `automask_area_frac` | `automask2D` | Particle-like classes should occupy a plausible mask area. | Automask adds parameters and runtime. |
| `automask_diameter_norm` | `automask2D` | Flags tiny specks and over-large blobs. | May reject unusual but real projections. |
| `automask_shift_norm` | `automask2D` | Independent centering signal. | Correlated with current `centered` feature. |
| `masked_integrated_intensity` | `exec_shape_rank_cavgs` pattern | Catches strong coherent particle density. | Contrast/sign conventions must be consistent. |

## Suggested Extension Strategy

Do not add a large feature batch in one pass. The current learner can evaluate all features, but feature interpretation becomes harder if many correlated features arrive together.

A practical sequence is:

1. Re-run `quality_mode=analyze` with the active cheap diagnostics: `hist_entropy`, `cc_area_frac`, `cc_diameter_norm`, `presence`, `log_contrast`, and `log_hist_variance`.
2. Inspect the learn report `feature_signal`, `feature_drop`, and `feature_group_lodo` rows. A feature should have stable per-dataset AUC, improve the one-feature-drop test, or improve leave-one-dataset-out group behavior before it is trusted as default-model signal.
3. Promote only features that improve holdout behavior, especially on difficult chunk datasets, without damaging easy/manual-trusted cases.
4. Revisit spectral additions only if the scalar feature bank cannot solve a representative failure mode; keep them diagnostic until they prove they do not simply rediscover ice/ring artifacts.
5. Consider automask or shape-ranking features only if the cheap scalar diagnostics are still insufficient.
6. Treat alignment-derived relational features as optional scalar experiments only; do not reintroduce pairwise matrices into the default model.

For any new feature:

- Add its index and definition in `simple_cavg_quality_feats.f90`.
- Keep the sign convention `higher_is_better` unless the feature is explicitly diagnostic.
- Write both raw and normalized values through the existing analysis path.
- Add it to model files through the existing `feature_weights` machinery.
- Re-run learning from analysis files; do not hand-tune weights first.
- Validate separately for `chunk` and `pool` contexts.

## Features to Avoid for Now

Some possible descriptors are available but should not be first-line additions:

- Raw unmasked image mean or total intensity: too dependent on normalization and background convention.
- Many texture features from `GLCM`: potentially useful, but easy to overfit and harder to interpret.
- Absolute thresholds on spectrum or ice scores: dataset and sampling dependent; prefer normalized features and learned weights.
- Full pairwise image alignment in the default stream path: potentially informative but likely too expensive unless timing proves otherwise.
- A black-box classifier that bypasses the current model file and analysis/reporting flow: this would make failures harder to debug.

## Implementation Notes

- Prefer stored `cls2D` metadata (`pop`, `res`, `corr`) before recomputing expensive quantities.
- Keep feature extraction in `simple_cavg_quality_feats.f90`; commanders should orchestrate, not own quality logic.
- Keep raw values in analysis files even if a feature has zero model weight.
- Use per-dataset robust normalization over non-hard-rejected classes.
- Treat hard rejection as exceptional and conservative. Most discrimination should happen through the model.
- Treat feature-order changes as breaking model-file changes: model files contain one weight per `CAVG_QUALITY_NFEATS` entry.
- Prefer appending candidate diagnostics over inserting them into the established feature order, and treat any feature-order change as a deliberate model-version break.
- When changing `CAVG_QUALITY_NFEATS`, update promotion snippets, docs, and any tests or analysis parsers that assume the feature count.
- Benchmark any relational or alignment-derived scalar feature on stream chunks before considering it for default use.
