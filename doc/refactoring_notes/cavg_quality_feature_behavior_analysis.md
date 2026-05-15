# Class-Average Quality Feature Behavior Analysis

Date: 2026-05-15

## Scope

This note summarizes how the old and current class-average rejection features behave on the
available chunk-style annotated runs under `/Users/elmlundho/cavgs_quality`.

Primary data sources:

- `prev/batch1chunk`: early scalar feature bank with `hist_knn`, `spectrum_dynrange`, and `neg_ice_score`.
- `prev/batch6chunk`: expanded bank including histogram/intensity softness and `log_fg_bg_locvar_ratio`.
- `prev/batch10chunk`: 15-feature bank before the global intensity-softness features were removed.
- `batch_train1`: current 16-feature bank with internal-detail features and no global intensity-softness trio.
- `prev/test_feat`: ApoF and FlipQR feature probe used to inspect the new internal-detail features.

Metrics below use trainable rows only, with hard-rejected rows excluded. Per-feature behavior is reported
as mean per-dataset AUC and median z-score separation between manual-good and manual-bad classes.
An AUC below 0.5 means the feature is inverted for that dataset.

## Feature-Bank Evolution

The feature bank evolved in three practical phases:

| Phase | Representative batch | Added or tested behavior | Outcome |
| --- | --- | --- | --- |
| Early scalar plus pairwise-like features | `prev/batch1chunk` | `hist_knn`, `spectrum_dynrange`, `neg_ice_score` | `hist_knn` and `spectrum_dynrange` were not general. `neg_ice_score` was flat in these rows. |
| Expanded scalar and intensity softness | `prev/batch6chunk`, `prev/batch10chunk` | `hist_entropy`, `log_contrast`, `log_hist_variance`, `log_fg_bg_locvar_ratio`, `presence`, connected-component shape | Presence was useful. The global histogram/intensity-softness trio looked strong on average but was easy to fool and inconsistent. |
| Current scalar/detail bank | `batch_train1` | Removed global softness and older spectral/ice features; added `log_detail_bg_snr`, `log_detail_signal_ratio`, `detail_coverage`, `detail_edge_density` | Core scalar features remain reliable. Internal detail is dataset-specific and should be treated cautiously. |

## Current Feature Behavior

Current 16-feature bank, trainable rows in `batch_train1`:

| Feature | Mean AUC | Min AUC | Inverted datasets | Learned weight | Drop delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| `corr_frc_proxy` | 0.867 | 0.699 | 0 | 0.1385 | -0.0429 |
| `neg_log_res` | 0.860 | 0.631 | 0 | 0.1299 | -0.0117 |
| `presence` | 0.814 | 0.574 | 0 | 0.1177 | -0.0050 |
| `log_locvar_bg` | 0.797 | 0.586 | 0 | 0.1077 | -0.0213 |
| `log_locvar_fg` | 0.781 | 0.500 | 0 | 0.1029 | -0.0273 |
| `log_pop` | 0.762 | 0.292 | 1 | 0.1025 | +0.0047 |
| `log_center_edge_snr` | 0.713 | 0.375 | 2 | 0.0736 | -0.0434 |
| `cc_diameter_norm` | 0.676 | 0.500 | 0 | 0.0721 | -0.0155 |
| `log_detail_bg_snr` | 0.654 | 0.245 | 4 | 0.0582 | -0.0158 |
| `centered` | 0.649 | 0.511 | 0 | 0.0414 | -0.0126 |
| `detail_edge_density` | 0.567 | 0.167 | 5 | 0.0307 | -0.0198 |
| `detail_coverage` | 0.564 | 0.250 | 5 | 0.0247 | -0.0032 |
| `cc_area_frac` | 0.635 | 0.307 | 4 | 0.0000 | 0.0000 |
| `log_detail_signal_ratio` | 0.504 | 0.101 | 6 | 0.0000 | 0.0000 |
| `mask_inside` | 0.500 | 0.500 | 0 | 0.0000 | 0.0000 |
| `single_component` | 0.500 | 0.500 | 0 | 0.0000 | 0.0000 |

Interpretation:

- `corr_frc_proxy` is the most important scalar evidence in the current run. It has the best drop
  diagnostic and is never inverted.
- `neg_log_res` remains one of the safest features. It is not always sufficient on hard chunks, but it
  is never directionally wrong in the current representative set.
- `presence` is a very useful general image-processing feature. It is weaker on CLC7 and TMEM16A but
  does not invert.
- `log_locvar_fg` and `log_locvar_bg` remain core features. They are strong on PepT2 and many easy cases,
  but their dataset minimums show they should not dominate alone.
- `log_pop` is useful overall but can invert. It is a support prior, not a direct quality measurement.
- `log_center_edge_snr` is powerful when correct, but it inverts on Not and TRPM4. It should remain
  a secondary feature.
- `cc_diameter_norm` is a helpful weak geometry scalar. It has no inversions in the current batch.
- `mask_inside` and `single_component` are flat in trainable rows because the catastrophic geometry
  cases are now handled as hard rejects. They should stay out of the learned model unless their
  continuous definitions become informative again.

## Internal-Detail Features

The internal-detail features were motivated by visually obvious cases like featureless ice blobs versus
secondary-structure-rich averages. They do capture real information, but not in a uniformly stable way.

Strong positive examples:

- ApoF: `log_detail_bg_snr` AUC 1.000, `log_detail_signal_ratio` AUC 1.000.
- TMEM16A: all detail features are strong, with `log_detail_bg_snr` AUC 0.965 and
  `detail_edge_density` AUC 0.936.
- CLC7: `log_detail_bg_snr` is strong, AUC 0.915.

Inversions:

- PepT2: `log_detail_bg_snr` AUC 0.245, `log_detail_signal_ratio` AUC 0.101.
- RNApol: all four detail features invert.
- TRPM4: all four detail features are weak or inverted.
- FlipQR: all four detail features invert in the small two-dataset probe.

Conclusion: internal detail is real but currently not a universally safe scalar family. It can help TMEM16A
and ApoF-like failures, but it explains the PepT2 penalty in the promoted `batch_train1` model. If we want
these features in the production model, the learner needs either stronger evidence that they generalize or
a policy that can reject the whole internal-detail family when it improves only a subset of datasets.

## Old Feature Behavior

### `hist_knn`

`hist_knn` was intended to preserve the useful histogram-distance intuition from the original
`cluster_cavgs` approach. On `prev/batch1chunk` it had:

- Mean AUC: 0.532
- Inverted datasets: 5 of 10
- Strong positive datasets: ApoF, Not, partly TMEM16A
- Clear inversions: PepT2, RNApol, TRPM4, FlipQR, MotAB

This is the strongest evidence against the old pairwise histogram-distance style as a general learned
feature. It can cluster some easy cases, but its direction is not stable enough to serve as scalar model
evidence.

### `spectrum_dynrange`

On `prev/batch1chunk`, `spectrum_dynrange` had mean AUC 0.589 and inverted on RNApol, TMEM16A,
and TRPM4. In `prev/batch6chunk`, the same pattern persisted. It can be very strong on ApoF and Not,
but it fails on several stream-like chunks. Removing it from the feature bank was justified.

### `neg_ice_score`

`neg_ice_score` was flat in the tested rows, with AUC 0.500 across the representative batches. That
does not mean ice-ring logic is useless as a hard-rejection concept, but the tested scalar did not behave
as model evidence.

### Global Histogram And Intensity Softness

In `prev/batch10chunk`, before removal:

| Feature | Mean AUC | Min AUC | Inverted datasets | Comment |
| --- | ---: | ---: | ---: | --- |
| `log_hist_variance` | 0.805 | 0.505 | 0 | Strong average, but mostly a contrast/softness proxy. |
| `log_contrast` | 0.798 | 0.455 | 1 | Strong on many datasets, inverted on ApoF. |
| `hist_entropy` | 0.714 | 0.152 | 1 | Highly dataset-specific, strongly inverted on ApoF. |
| `log_fg_bg_locvar_ratio` | 0.591 in `batch6chunk` | 0.218 | 2 | Weak and unstable. |

The intensity-softness trio was tempting because it scored well on many datasets, including PepT2.
The problem is that it measures global image contrast/statistical softness, not reliable molecular
quality. It can be fooled by ice, masks, staining-like density, and processing differences. Removing it
from the model was a conservative choice even though it cost some apparent in-sample signal.

## Feature-Policy Evidence

The `batch_train1` learner selected `full_geom_pruned`, but the holdout-style group screen is more
cautious:

| Policy/group | Mean AUC | Min AUC dataset | Mean oracle balanced accuracy | Total FP | Total FN |
| --- | ---: | --- | ---: | ---: | ---: |
| `base_scalar` | 0.883 | RNApol | 0.879 | 73 | 13 |
| `base_no_soft_geom` | 0.883 | RNApol | 0.879 | 73 | 13 |
| `base_geom_pruned` | 0.875 | RNApol | 0.876 | 78 | 10 |
| `full_geom_pruned` | 0.874 | RNApol | 0.867 | 70 | 23 |
| `internal_detail4` | 0.597 | FlipQR | 0.701 | 126 | 45 |

This is important: direct training selected a policy that includes weak internal-detail weights, but the
group screen says the old 12-feature scalar family generalizes better. The promoted model is defensible
because PepT2 manual selection is less reliable and the macro score improved, but the feature evidence
does not yet prove that the detail features should be permanent production evidence.

## Dataset-Specific Notes

PepT2:

- Core scalar features are strong: `presence` 0.893, `log_locvar_bg` 0.886,
  `corr_frc_proxy` 0.863, `log_locvar_fg` 0.860.
- Internal-detail features invert badly. This is the main reason the promoted model loses PepT2 recall.
- The older global softness features also looked strong on PepT2, but they were removed because they
  are easy to fool elsewhere.

RNApol and TRPM4:

- Resolution and correlation are weaker than on easier datasets.
- Presence and connected-component shape are useful.
- Internal-detail features invert, so they are not reliable rescue evidence here.

TMEM16A:

- Detail features are unusually strong.
- The model still admits many false positives because a lot of junk is statistically close to good class
  averages. TMEM16A is the clearest case where we need either better molecular-detail features or a
  more explicit hard-rejection concept.

ApoF:

- New detail features behave exactly as hoped.
- Histogram entropy is strongly inverted, supporting the decision to avoid global histogram softness.

## Practical Conclusions

1. The stable model backbone is scalar, not pairwise: `corr_frc_proxy`, `neg_log_res`, `presence`,
   local variance, and selected weak geometry.
2. Pairwise/clustering-style histogram evidence did not generalize. The old `hist_knn` behavior is too
   often inverted.
3. Global intensity softness was useful in-sample but not trustworthy enough for production model evidence.
4. Internal-detail features capture something real, but they are not ready to be treated as generally safe.
5. Hard rejects and learned model features should stay conceptually separate. Geometry and catastrophic
   resolution failures belong in hard rejection; the model should learn among remaining trainable rows.

## Suggested Next Experiments

1. Run one training comparison where the learner can choose between `base_scalar`, `base_geom_pruned`,
   and `full_geom_pruned`, but score or tie-break using the feature-group holdout diagnostics rather than
   direct in-sample macro balanced accuracy only.
2. Add a model policy that excludes the entire `internal_detail` family, then compare it against the
   promoted model on the same fresh validation batch. This is not hand-zeroing individual features; it is
   testing whether a whole experimental family generalizes.
3. Keep the internal-detail features in the feature bank for diagnostics. They are informative on ApoF and
   TMEM16A, and they may guide a better future feature that is less sensitive to dataset-specific texture.
4. If the goal is to catch featureless ice blobs specifically, design a targeted feature that measures
   organized internal structure relative to expected molecular scale, not just high-frequency detail power.

