# Class-Average Quality Feature Behavior

## Scope

This note describes the active feature bank used by `model_cavgs_rejection`.

## Active Features

| Feature | Model role |
| --- | --- |
| `log_pop` | Measures class population support. |
| `neg_log_res` | Measures nominal class-average resolution; catastrophic resolution failure is also a hard reject. |
| `centered` | Measures foreground centering after segmentation. |
| `log_locvar_fg` | Measures foreground local variance. |
| `log_locvar_bg` | Measures background local variance. |
| `corr_frc_proxy` | Uses stored class correlation/FRC-like metadata when present. |
| `log_center_edge_snr` | Measures central signal relative to edge variance. |
| `cc_area_frac` | Measures foreground connected-component area relative to the mask. |
| `presence` | Measures central signal excess relative to border background noise. |

## Model Use

The model consumes robustly normalized `z_*` versions of these features. Learned weights define the scalar score and the feature dimensions used for k-medoids distance calculations.

Hard rejects are applied before model fitting. The model boundary is trained and evaluated only on non-hard-rejected class averages.

## Reports

`quality_mode=analyze` reports raw features, normalized features, model scores, states, feature summaries, and threshold scans.

`quality_mode=learn` reports learned weights, candidate scores, Otsu ablation diagnostics, feature-signal diagnostics, feature-drop diagnostics, and per-dataset confusion metrics.
