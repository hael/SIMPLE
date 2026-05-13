# Class-Average Quality Metric Inventory

**Date:** May 12, 2026  
**Status:** Developer inventory / metric design input  
**Scope:** Existing SIMPLE routines that can be reused to derive scalar features for separating good and bad cryo-EM 2D class averages.

## Purpose

This note inventories image-processing and class-average routines already present in SIMPLE that can be reused to build scalar quality metrics for 2D class averages. The goal is to avoid inventing a new image-quality stack when the library already has useful signal, mask, resolution, connected-component, spectral, and clustering primitives.

The strongest path is to treat class-average quality as a small feature-vector problem:

1. Extract robust scalar descriptors from each class average.
2. Keep raw values and normalized values.
3. Apply hard rejection for pathological cases.
4. Combine the remaining descriptors with robust normalization, simple weights, or a small clustering step.
5. Persist both the final quality score and the raw feature table so that bad decisions can be diagnosed.

## Current Worktree Capability

The current worktree contains a directly relevant work-in-progress module:

- `src/main/strategies/search/simple_cavg_quality.f90`
  - `evaluate_cavg_quality`
  - `extract_cavg_quality_features`
  - `normalize_cavg_quality_features`
  - `cluster_cavg_quality`
  - `cavg_quality_result`

This module already implements an eight-feature class-average quality vector:

| Feature | Raw source | Quality direction |
| --- | --- | --- |
| `log_pop` | `cls2D` population | larger is better |
| `neg_log_res` | `cls2D` resolution | larger is better, because lower resolution in Angstrom is better |
| `mask_inside` | largest Otsu foreground component outside mask disc | larger is better after sign flip |
| `centered` | foreground component centroid offset normalized by mask radius | larger is better after sign flip |
| `log_locvar_fg` | masked foreground local variance | larger is usually better |
| `log_locvar_bg` | masked background local variance | larger is usually better |
| `spectrum_dynrange` | square-root spectrum at high-pass and low-pass indices | larger is usually better |
| `single_component` | connected-component count penalty | closer to one component is better |

The current worktree also has a commander hook in `src/main/commanders/simple/simple_commanders_cavgs.f90`:

- `exec_cluster_cavgs_quality`
  - reads class averages from a project
  - calls `evaluate_cavg_quality`
  - writes `cavgs_quality_features.txt`
  - writes selected and rejected stacks
  - annotates `cls2D` fields `quality`, `quality_cluster`, and `accept`
  - uses an internal signed boundary margin; positive values retain more
    borderline classes, while negative values raise the boundary and reject more
    junk
  - keeps `spectrum_dynrange` in the feature table for diagnostics, but does
    not currently give it positive score weight because high values can reflect
    ring or ice artifacts rather than good class-average signal

This is already very close to a reusable scalar metric pipeline. The remaining inventory below lists additional existing routines that can either support this module or expand it.

## Ready-Made Class-Average Signals

### Cluster2D Rejection

Source:

- `src/main/stream/simple_cluster2D_rejector.f90`

Reusable routines:

- `reject_pop`
- `reject_res`
- `reject_mask`
- `reject_local_variance`

These routines already encode practical good/bad criteria:

- population rejection using a small fraction of the total population
- resolution rejection using `RES_THRESHOLD = 40.0`
- mask rejection using Otsu foreground, connected components, and outside-mask pixels
- flat-class rejection using foreground/background local variance and robust z-scores

Metric candidates:

- population fraction
- resolution pass/fail or continuous resolution score
- outside-mask pixel count or outside-mask fraction
- largest-component centroid offset
- foreground and background local variance
- robust local-variance z-score
- hard-reject flag

### Older Non-Junk Gate

Source:

- `src/main/strategies/search/simple_strategy2D_utils.f90`

Reusable routine:

- `flag_non_junk_cavgs`

This routine combines several lightweight quality tests:

- spectral dynamic range
- density inside/outside the mask
- center offset
- minimum class population

Metric candidates:

- spectrum dynamic range between chosen spatial-frequency bands
- density ratio inside the mask
- center-offset penalty
- population threshold flag
- combined non-junk flag

### Cluster-Based Class-Average Scores

Source:

- `src/main/strategies/search/simple_strategy2D_utils.f90`
- `src/defs/simple_type_defs.f90`

Reusable routines and fields:

- `calc_sigstats_dmats`
- `calc_cc_and_res_dmats`
- `calc_cluster_cavgs_dmat`
- `align_and_score_cavg_clusters`
- `clust_info%euclid`
- `clust_info%homogeneity`
- `clust_info%resscore`
- `clust_info%clustscore`
- `clust_info%jointscore`
- `clust_info%good_bad`
- `clust_info%icescore`

These routines already build distance matrices and convert cluster properties into percent-style scores. They are useful when quality is better judged by consistency with other class averages rather than by single-image statistics alone.

Metric candidates:

- distance to medoid
- mean distance to assigned cluster
- cluster homogeneity percentile
- cluster resolution percentile
- cluster score percentile
- joint score
- good/bad cluster label

### FRC and Resolution Infrastructure

Sources:

- `src/main/class/simple_class_frcs.f90`
- `src/main/class/simple_classaverager_restore.f90`
- `src/main/class/simple_classaverager_core.f90`
- `src/utils/filter/simple_estimate_ssnr.f90`
- `src/utils/filter/simple_opt_filter.f90`

Reusable routines:

- `class_frcs%frc_getter`
- `class_frcs%avg_frc_getter`
- `class_frcs%estimate_res`
- `class_frcs%estimate_find_for_eoavg`
- `class_frcs%estimate_lp_for_align`
- `fsc2ssnr`
- `fsc2optlp`
- `estimate_lplim`
- `estimate_lplims2D`

The restoration path already computes even/odd FRCs and persists per-class `pop`, `res`, `corr`, and `w` fields. These are high-value quality signals and should be reused before recomputing anything expensive.

Metric candidates:

- FRC 0.5 resolution
- FRC 0.143 resolution
- mean FRC over a target band
- first frequency index below threshold
- estimated alignment low-pass limit
- even/odd class-average correlation
- SSNR-derived quality summary

### Shape Ranking

Source:

- `src/main/commanders/simple/simple_commanders_cavgs.f90`

Reusable routine:

- `exec_shape_rank_cavgs`

This path uses `automask2D`, integrated masked intensity, object diameter, and center shifts to rank non-junk class averages.

Metric candidates:

- automask diameter
- automask center shift
- integrated masked intensity
- shape rank
- masked area or occupancy fraction

## Reusable Image-Level Primitives

### Masking and Connected Components

Sources:

- `src/main/image/simple_image_msk.f90`
- `src/main/image/simple_image_bin.f90`
- `src/main/image/simple_image.f90`

Reusable routines:

- `automask2D`
- `density_inoutside_mask`
- `image%density_inoutside`
- `image%nforeground`
- `image%nbackground`
- `image_bin%find_ccs`
- `image_bin%size_ccs`
- `image_bin%diameter_cc`
- `image_bin%masscen_cc`
- `image_bin%order_ccs`
- `image_bin%cc2bin`

Metric candidates:

- foreground area
- background area
- foreground/background density ratio
- largest connected-component area
- connected-component count
- largest-component diameter
- largest-component centroid offset
- outside-mask fraction
- occupancy fraction within expected particle mask

### Local Signal and Texture

Source:

- `src/main/image/simple_image.f90`

Reusable routines:

- `image%loc_sdev`
- `image%avg_loc_sdev`
- `image%loc_var`
- `image%loc_var_masked`
- `image%GLCM`
- `image%variance`
- `image%skew`
- `image%kurt`
- `image%noisesdev`

Metric candidates:

- average local standard deviation
- average local variance
- foreground local variance
- background local variance
- foreground/background local-variance ratio
- image variance
- skewness
- kurtosis
- gray-level co-occurrence texture summaries

### Spectral and Frequency-Domain Signal

Sources:

- `src/main/image/simple_image.f90`
- `src/main/simple_pspecs.f90`

Reusable routines:

- `image%spectrum`
- `image%power_spectrum`
- `image%guinier_bfac`
- `image%guinier`
- `image%fsc`
- `image%frc_pspec`
- `image%fsvar`
- `image%get_res`
- `image%fcomps_below_noise_power_stats`
- `pspecs%new`

Metric candidates:

- spectrum dynamic range
- high-frequency to low-frequency power ratio
- band-limited power ratios
- Guinier slope or B-factor proxy
- estimated resolution from FRC/FSC curve
- fraction of Fourier components below noise power
- power-spectrum distance to good-class medoid

### Correlation, Similarity, and Alignment Consistency

Sources:

- `src/main/image/simple_image.f90`
- `src/main/strategies/search/simple_strategy2D_utils.f90`

Reusable routines:

- `image%corr`
- `image%corr_shifted`
- `image%real_corr`
- `image%phase_corr`
- `image%fcorr_shift`
- `image%radial_cc`
- `image%sqeuclid`
- `image%weighted_subtr_corr`
- `image%masked_subtr_corr`
- `calc_cc_and_res_dmats`
- `calc_cluster_cavgs_dmat`

Metric candidates:

- nearest-neighbor correlation
- shifted correlation to a medoid
- phase-correlation peak strength
- radial correlation summary
- masked Euclidean distance to medoid
- cluster-average similarity
- consistency of class average with nearby class averages

### Centering, Presence, Contrast, and Sanity Checks

Source:

- `src/main/image/simple_image.f90`

Reusable routines:

- `image%center_edge_snr`
- `image%presence`
- `image%contrast`
- `image%contains_nans`
- `image%rmsd`
- `image%snr`

Metric candidates:

- center-edge SNR
- particle-presence score
- contrast
- NaN flag
- RMSD to reference, medoid, or local average
- simple SNR estimate

### Ice and Artifact Scores

Source:

- `src/main/image/simple_image.f90`

Reusable routine:

- `image%calc_ice_score`

Metric candidates:

- ice score
- ice-score percentile within a run
- hard artifact flag for extreme ice-score outliers

### Histogram-Based Statistics

Source:

- `src/utils/math/simple_histogram.f90`

Reusable routines:

- `histogram%variance`
- `histogram%entropy`
- histogram moments and distance measures

Metric candidates:

- bounded entropy of intensity histogram
- histogram variance
- histogram skewness
- histogram distance to good-class distribution
- Jensen-Shannon or Hellinger distance to a reference distribution

## Recommended First Metric Vector

The current `simple_cavg_quality` feature vector is a good first implementation because it is scalar, interpretable, robustly normalized, and tied to existing class-average rejection logic. For a first production-quality pass, the following vector is a practical target:

| Feature | Existing source | Notes |
| --- | --- | --- |
| `log_pop` | `cls2D%pop` | Avoids over-penalizing large-population classes by using log scale. |
| `neg_log_res` | `cls2D%res`, `class_frcs%estimate_res` | Reuse stored resolution when available. |
| `corr` or FRC-band mean | `cls2D%corr`, `class_frcs%frc_getter` | Adds more information than a single resolution cutoff. |
| `outside_mask_frac` | `density_inoutside_mask`, connected components | Penalizes off-mask junk. |
| `centroid_offset` | `image_bin%masscen_cc` | Penalizes off-center particles. |
| `ncc_penalty` | `image_bin%find_ccs` | Penalizes fragmented foreground. |
| `locvar_fg` | `image%loc_var_masked` | Detects flat or featureless foreground. |
| `locvar_ratio` | `image%loc_var_masked` | Compares foreground and background texture. |
| `spectrum_dynrange` | `image%spectrum` | Detects weak or overfiltered averages. |
| `center_edge_snr` | `image%center_edge_snr` | Captures particle-like center signal. |
| `contrast` | `image%contrast` | Lightweight single-image contrast signal. |
| `ice_score` | `image%calc_ice_score` | Useful artifact indicator when calibrated per data set. |

The scalar score should keep each feature direction explicit. A simple robust approach is:

1. Hard-reject NaNs, empty populations, missing foreground components, and near-zero local variance.
2. Normalize each feature by median and MAD over non-rejected classes.
3. Clip z-scores to a small range such as `[-4, 4]`.
4. Combine with transparent weights.
5. Optionally split scores into good and bad by k-medoids or by an adaptive threshold.

## Suggested Composite Outputs

Useful outputs for developers and users:

- `quality`: scalar normalized score, larger means better.
- `accept`: hard selection state, `1` for accepted and `0` for rejected.
- `quality_cluster`: cluster label when clustering was used.
- `hard_reject`: diagnostic hard-rejection flag.
- `quality_reason`: optional future field for the dominant rejection reason.
- `cavgs_quality_features.txt`: raw and normalized feature table.

The current commander already writes the feature table and selected/rejected stacks. That is the right diagnostic shape.

## Implementation Notes

- Prefer reusing stored `cls2D` metadata (`pop`, `res`, `corr`) before recomputing FRCs.
- Keep feature extraction in a standalone module, not buried in a commander, so stream and project workflows can share it.
- Preserve both raw and normalized values. Raw values are essential for debugging thresholds.
- Normalize per run with robust statistics. Absolute thresholds are useful for hard rejection, but dataset-level contrast and ice behavior vary.
- Document the sign of every scalar feature. Every combined feature should have the same direction before scoring.
- Treat resolution and population as correlated features. A quality score should not become only a population score.
- Avoid training an opaque classifier until the raw feature table has been validated across several datasets.
- Keep mask radius and sampling distance explicit in every metric that depends on physical scale.

## Highest-Value Additions Beyond The Current WIP Module

The current worktree feature vector can be strengthened with these existing routines:

1. Add `image%center_edge_snr` as a direct particle-presence metric.
2. Add `image%contrast` to detect extremely weak averages.
3. Add `image%calc_ice_score`, but report it separately until calibrated.
4. Add an FRC-band summary from `class_frcs%frc_getter`, not only the final resolution estimate.
5. Add medoid-distance or nearest-neighbor consistency from `calc_cluster_cavgs_dmat` when enough classes are present.

These additions would make the scalar quality score sensitive to alignment consistency, spectral support, particle localization, artifact content, and class population without leaving the existing SIMPLE image-processing library.
