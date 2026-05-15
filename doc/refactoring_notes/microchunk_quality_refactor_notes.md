# Microchunk Class-Average Quality Backend

## Scope

This note documents the current integration contract between `cluster2D_microchunked`/`simple_microchunked2D` orchestration and the shared class-average quality backend in `src/main/cavg_quality`.

## Backend Selection

The microchunk layer controls backend selection with an internal logical flag:

```fortran
logical, parameter :: USE_CAVG_QUALITY_BACKEND = .true.
```

This flag is not a command-line parameter. When enabled, microchunk class-average rejection uses the shared class-average quality backend with `CAVG_QUALITY_MODEL_CHUNK_DEFAULT`.

## Backend Call

The quality backend path:

1. Instantiates `type(cavg_quality_model)`.
2. Initializes `CAVG_QUALITY_MODEL_CHUNK_DEFAULT`.
3. Calls `evaluate_cavg_quality` on the chunk class-average stack, class-average orientations, and mask diameter.
4. Reads `quality%states`, `quality%scores`, `quality%labels`, and `quality%hard_reject`.
5. Writes the accepted/rejected class states into the same project fields consumed by the existing microchunk pipeline.

The microchunk scheduler, sentinel files, project partitioning, and chunk-combination logic remain outside the quality backend.

## Model

The enabled backend uses `chunk_default_v2` unless another model is explicitly initialized by the caller.

`chunk_default_v2` uses the fixed 9-feature scalar bank:

- `log_pop`
- `neg_log_res`
- `centered`
- `log_locvar_fg`
- `log_locvar_bg`
- `corr_frc_proxy`
- `log_center_edge_snr`
- `cc_area_frac`
- `presence`

Hard validity rejects are applied before model fitting. Non-hard-rejected class averages are scored, clustered in normalized feature space, and thresholded by the model controls.

## Data Flow

The backend returns:

- `quality%states`: accepted/rejected class-average state.
- `quality%scores`: scalar model score.
- `quality%labels`: feature-space cluster label.
- `quality%medoids`: class indices of cluster medoids.
- `quality%hard_reject`: hard validity reject mask.
- `quality%threshold`: final score threshold.
- `quality%raw_threshold`: cluster-derived score threshold before model controls.
- `quality%threshold_margin`: difference between raw and final threshold.
- `quality%separation`: score separation between the two feature-space clusters.
- `quality%nclust`: number of clusters used by the classifier.

Microchunk code maps these outputs into the project fields already used for class-average selection and particle-state propagation.

## Diagnostics

`model_cavgs_rejection` remains the diagnostic and training harness for the same backend. The microchunk path uses the backend result directly and does not need to write full analysis tables during routine streaming runs.
