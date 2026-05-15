# Class-Average Quality Metric Inventory

**Date:** May 14, 2026
**Status:** Developer policy / future feature planning
**Scope:** Existing SIMPLE routines that can support scalar features for `model_cavgs_rejection`.

## Purpose

This note inventories image-processing and class-average routines already present in SIMPLE that can support the `model_cavgs_rejection` scalar feature bank. The goal is to keep quality selection interpretable and reproducible: every feature should have a clear source, a documented direction, a feature family, and a diagnostic raw value in the analysis table.

The current implementation treats class-average quality as a feature-vector and model-selection problem:

1. Extract raw scalar descriptors for each 2D class average.
2. Apply only conservative hard rejections for pathological cases.
3. Normalize non-hard-rejected classes within the dataset by robust statistics.
4. Use an explicit quality model to score and cluster the class averages.
5. Preserve the raw features, normalized features, model parameters, and manual reference labels in `cavgs_quality_analysis.txt`.
6. Use `quality_mode=learn` to tune interpretable model parameters against trusted manual selections.

The current feature policy is deliberately scalar. Pairwise histogram/spectrum matrices and relational clustering-style quality metrics are not part of the default path. The supported direction is vectorized learning over normalized scalar quality metrics, with clustering used only to help choose an automatic decision threshold inside that scalar feature space.

## Current Implementation

The code lives under `src/main/cavg_quality`:

- `simple_cavg_quality_types.f90`: shared constants and derived types.
- `simple_cavg_quality_feats.f90`: feature definitions, feature-family metadata, feature extraction, and robust normalization.
- `simple_cavg_quality_stats.f90`: binary metrics, AUC, robust summaries, and distance-matrix normalization helpers.
- `simple_cavg_quality_model.f90`: built-in model presets, model file I/O, linear scoring, scalar feature-space clustering, threshold policy, and promotion-code generation.
- `simple_cavg_quality_analysis.f90`: `apply`/`analyze` evaluation and reporting.
- `simple_cavg_quality_learn.f90`: analysis-table reader and grid search over interpretable model parameters.

The command entry point is `exec_model_cavgs_rejection` in `src/main/commanders/simple/simple_commanders_cavgs.f90`.

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
| 10 | `cc_area_frac` | `image_bin%size_ccs` plus Otsu CC mask | diagnostic | Main connected-component area divided by the expected circular mask area. |
| 11 | `detail_texture` | 20-6 A band-pass texture | higher is better | Central-mask local residual-energy heterogeneity plus radial-residual detail after subtracting the radial mean. |
| 12 | `presence` | `image%presence` | higher is better | Central signal excess relative to border background noise. |
| 13 | `log_detail_bg_snr` | 20-6 A band-pass | diagnostic | Log in-mask detail RMS relative to outside-mask detail RMS. |
| 14 | `log_detail_signal_ratio` | 20-6 A band-pass | diagnostic | Log in-mask detail RMS relative to in-mask total signal RMS. |
| 15 | `detail_coverage` | 20-6 A band-pass | diagnostic | Fraction of in-mask pixels with detail above background and image-scale thresholds. |
| 16 | `detail_edge_density` | 20-6 A band-pass | diagnostic | Fraction of in-mask pixels with band-passed gradient above background and image-scale thresholds. |

Pairwise histogram and rotational-spectrum matrices are outside the model infrastructure. Global intensity-softness descriptors are not part of the feature bank.

The model policy favors scalar features with stable quality direction under robust per-dataset normalization.

The internal-detail features are zero-weighted in the current built-in defaults. They are meant to answer a specific validation question: whether smooth blob-like class averages that pass the scalar geometry checks are separable by direct measurements of organized in-mask molecular texture.

Resolution is both a hard validity check and an active model feature. The feature extractor hard rejects class averages when `cls2D res > 40 A`, but non-catastrophic resolution variation remains in the learned model through `neg_log_res`. Mask-geometry validity is handled as a hard gate when the Otsu/connected-component pass finds no valid foreground component, when any foreground-component centroid lies outside the mask radius, or when more than 10 pixels of the largest component lie outside the mask disc. These hard rejections are implemented in the quality library and do not depend on stream/microchunk modules.

## Current Model Controls

The built-in models are complete `linear_boundary` presets:

- `chunk_default_v2`: current default chunk/stream behavior, guarded by the `full_geom_pruned` feature policy.
- `chunk_default_v1`: alternate chunk/stream preset.
- `pool_default_v1`: recall-preserving behavior for larger pooled/batch datasets.

The important model controls are:

- `feature_policy`: learn-mode candidate name describing which scalar features are allowed into the score and clustering distance; the model file encodes this by zeroing excluded weights. The current chunk default zeros `mask_inside` and `single_component` because those soft geometry flags duplicate the hard connected-component rejection, and zeros `cc_area_frac` because its direction is not stable enough for the chunk default.
- `feature_weights`: nonnegative linear score coefficients, normalized to sum to one.
- `boundary_margin`: shifts the decision threshold; more negative rejects more aggressively, more positive preserves more borderline classes.
- `min_score_separation`: below this separation the model treats the partition as unreliable unless low-separation Otsu is enabled and acceptable.
- `otsu_min_offset` and `otsu_max_offset`: constrain when an Otsu score threshold may replace the cluster-derived threshold.
- `cluster_rescue_margin`: lets pool-like models rescue good-cluster members close to threshold.
- `min_accept_frac`: enforces a minimum accepted fraction for pool-like models.

The current learner searches feature-policy candidates, feature weights derived from the training data, minimum score separation, boundary margin, Otsu threshold controls, and pool minimum accepted fraction. The policy grid includes base, full, and geometry-pruned variants, including `full_geom_pruned`.

## Reusable Existing Routines

### Microchunk Rejector

Source:

- `src/main/stream/simple_cluster2D_rejector.f90`

Relevant routines:

- `reject_pop`
- `reject_res`
- `reject_mask`
- `reject_local_variance`

These routines define operationally useful rejection criteria and candidate feature sources. The quality app uses the strongest scalar pieces directly: population, resolution, mask geometry, connected components, and local variance.

### Non-Junk Gate Utilities

Source:

- `src/main/strategies/search/simple_strategy2D_utils.f90`

Relevant routine:

- `flag_non_junk_cavgs`

This path combines density inside/outside a mask, center offset, and minimum population. Treat it as a source of possible scalar diagnostics, not as a quality-selection policy.

### Relational Similarity Routines

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

These routines produce relational or pairwise class-average evidence. They are not part of the current quality policy. Any future use should reduce their output to cheap scalar diagnostics and must outperform the current scalar feature-vector learner on held-out datasets before becoming model-eligible.

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

This path uses `automask2D`, integrated masked intensity, and center shifts to rank class averages. Several of those outputs are attractive quality features if they can be computed without making quality selection too slow or too dependent on automask parameters.

### Image-Level Primitives

Sources:

- `src/main/image/simple_image.f90`
- `src/main/image/simple_image_calc.f90`
- `src/main/image/simple_image_msk.f90`
- `src/main/image/simple_image_bin.f90`
- `src/main/simple_pspecs.f90`
- `src/utils/math/simple_stat.f90`

Useful primitives include:

- `image%loc_var_masked`, `image%variance`, `image%skew`, `image%kurt`, `image%GLCM`
- `image%spectrum`, `image%power_spectrum`, `image%guinier`, `image%guinier_bfac`, `image%frc_pspec`
- `image%fcomps_below_noise_power_stats`
- `image%center_edge_snr`, `image%presence`, `image%contains_nans`
- `image%corr`, `image%corr_shifted`, `image%real_corr`, `image%phase_corr`, `image%radial_cc`, `image%sqeuclid`
- `image_bin%find_ccs`, `image_bin%size_ccs`, `image_bin%masscen_cc`
- `automask2D`, `density_inoutside_mask`, `image%density_inoutside`

## Feasible Future Feature Extensions

The best future additions are cheap, interpretable, and likely to generalize across datasets. They should be added one or two at a time, then validated through `quality_mode=analyze` and `quality_mode=learn` before changing built-in weights.

### Low-Risk Scalar Features

These are feasible first additions because they reuse existing per-image routines and do not require pairwise alignment.

| Candidate | Existing source | Why it may help | Caution |
| --- | --- | --- | --- |
| `presence` variants | `image%presence` | Direct particle-presence proxy. | Inspect scale and sign across stream and pool datasets before assigning much model weight. |
| `nan_or_bad_pixel_flag` | `image%contains_nans`, min/max checks | Hard diagnostic for corrupt inputs. | Should probably remain a hard reject or warning, not a learned soft feature. |

### Medium-Risk Spectral Features

These are plausible and cheap enough, but they should start as diagnostics.

| Candidate | Existing source | Why it may help | Caution |
| --- | --- | --- | --- |
| `bandpower_mid_low_ratio` | `image%spectrum` or `image%power_spectrum` | Separates real mid-frequency structure from low-frequency blobs. | Ring artifacts can mimic high quality. |
| `bandpower_high_mid_ratio` | `image%spectrum` | Catches over-sharpened/overfitted averages. | Needs artifact-aware calibration. |
| `guinier_slope` or `bfac_proxy` | `image%guinier`, `image%guinier_bfac` | Captures spectral falloff shape. | Sensitive to preprocessing and box/mask choices. |
| `below_noise_fraction` | `image%fcomps_below_noise_power_stats` | Identifies weak averages with little usable Fourier signal. | Needs a stable noise definition for class averages. |
| `frc_pspec_summary` | `image%frc_pspec` | May approximate consistency without full even/odd class FRC infrastructure. | Must verify the input assumptions for class-average images. |

Do not add spectrum-derived defaults casually. The current policy does not include rotational-spectrum pairwise distances.

### Relational Features Outside Current Policy

These routines are available in the code base, but pairwise or alignment-derived quality evidence is outside the current model policy. The default path should remain scalar, dataset-normalized, and learnable. If any relational idea is reconsidered, it should be converted to one cheap scalar diagnostic and validated as a normal feature-family candidate.

| Candidate | Existing source | Why it may help | Caution |
| --- | --- | --- | --- |
| `sig_knn` | `calc_sigstats_dmats` | Could summarize signal-statistics neighborhoods as scalar support. | Scalar neighborhood density can invert across datasets. |
| `cc_knn` | `calc_cc_and_res_dmats` | Good classes should correlate with nearby good classes after in-plane search. | Expensive for stream chunks; consider optional analyze/learn diagnostics only. |
| `res_knn` | `calc_cc_and_res_dmats` | Pairwise FRC/resolution consistency may catch overfit junk. | Depends on robust pairwise matching. |
| `mixed_support_knn` | `calc_cluster_cavgs_dmat` | Summarizes signal/correlation/resolution neighborhood support. | Do not wire the matrix itself into the default model. |
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
| `automask_shift_norm` | `automask2D` | Independent centering signal. | Correlated with current `centered` feature. |
| `masked_integrated_intensity` | `exec_shape_rank_cavgs` pattern | Catches strong coherent particle density. | Contrast/sign conventions must be consistent. |

## Suggested Extension Strategy

Do not add a large feature batch in one pass. The current learner can evaluate all features, but feature interpretation becomes harder if many correlated features arrive together.

A practical sequence is:

1. Re-run `quality_mode=analyze` with the active cheap diagnostics: `cc_area_frac`, `detail_texture`, `presence`, and the internal-detail features.
2. Inspect the learn report `feature_signal`, `feature_drop`, and `feature_group_lodo` rows. A feature should have stable per-dataset AUC, improve the one-feature-drop test, or improve leave-one-dataset-out group behavior before it is trusted as default-model signal.
3. Promote only features that improve holdout behavior, especially on difficult chunk datasets, without damaging easy/manual-trusted cases.
4. Revisit spectral additions only if the scalar feature bank cannot solve a well-defined failure mode; keep them diagnostic until they prove they do not simply rediscover periodic artifacts.
5. Consider automask or shape-ranking features only if the cheap scalar diagnostics are still insufficient.
6. Treat alignment-derived relational features as optional scalar diagnostics only; do not introduce pairwise matrices into the default model.

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
- Absolute thresholds on spectral scores: dataset and sampling dependent; prefer normalized features and learned weights.
- Full pairwise image alignment in the default stream path: outside the current scalar feature-vector policy and too expensive for routine stream use.
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
