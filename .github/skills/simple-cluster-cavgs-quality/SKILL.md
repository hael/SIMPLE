---
name: simple-cluster-cavgs-quality
description: Use when working on SIMPLE's model_cavgs_rejection application, simple_cavg_quality feature-space selector, class-average quality selection, quality_mode=analyze outputs, cavgs_quality feature/reference tables, chunk vs pool rejection tuning, or microchunk quality-vector integration.
---

# SIMPLE `model_cavgs_rejection`

Use this skill for the class-average quality application and for analysis of `quality_mode=analyze` test outputs.

## Mental Model

`model_cavgs_rejection` is a `simple_exec` class-average processing program. It is not routed through `simple_exec_cluster2D`.

Flow:

`simple_ui_cavgproc -> simple_exec_cavgproc -> commander_model_cavgs_rejection -> simple_cavg_quality_*`

The application reads class averages from an input SIMPLE project, evaluates feature-space quality, writes selected/rejected class-average stacks, and either updates the project (`quality_mode=apply`) or compares against the current manual `cls2D` state without modifying the project (`quality_mode=analyze`).

## Key Files

- `src/main/ui/simple/simple_ui_cavgproc.f90`: UI registration for `model_cavgs_rejection`.
- `src/main/ui/simple_ui_params_common.f90`: `quality_mode` and `rejection_type` parameter definitions.
- `src/main/exec/simple_exec_cavgproc.f90`: routes `model_cavgs_rejection` to the commander.
- `src/main/commanders/simple/simple_commanders_cavgs.f90`: `exec_model_cavgs_rejection`, output writing, and project annotation.
- `src/main/cavg_quality/simple_cavg_quality_feats.f90`: feature extraction and normalization.
- `src/main/cavg_quality/simple_cavg_quality_model.f90`: built-in decision models, classification, cached decisions, and model promotion helpers.
- `src/main/cavg_quality/simple_cavg_quality_analysis.f90`: analyze-mode reporting, feature tables, threshold scans, and reference summaries.
- `src/main/cavg_quality/simple_cavg_quality_learn.f90`: learn/evaluate mode model fitting and reporting.
- `src/main/cavg_quality/simple_cavg_quality_types.f90`: shared result/model data structures.
- `doc/microchunk_and_rejection/model_cavgs_rejection.md`: current model behavior, feature bank, built-ins, reports, learning, and promotion.
- `doc/microchunk_and_rejection/microchunk_rejection_model_integration.md`: current stream/model boundary.

## Parameters

- `quality_mode=apply`: writes `cavgs_quality_features.txt`, selected/rejected stacks, annotates `cls2D` and `cls3D` with `quality`, `quality_cluster`, and `accept`, maps selection into particles, optionally prunes, and writes the project.
- `quality_mode=analyze`: treats the existing `cls2D` state as the manual reference, writes diagnostic outputs, and leaves the project selection unchanged.
- `quality_mode=learn`: learns a model from analysis files listed in `filetab`.
- `quality_mode=evaluate`: evaluates a selected model against analysis files.
- `quality_mode=promote`: emits Fortran snippets for promoting a model into built-in presets.
- `rejection_type=chunk`: stricter stream-partition operating point for high-junk chunks.
- `rejection_type=pool`: more recall-preserving operating point for larger pooled or batch sets.
- `quality_model`: built-ins currently include `chunk_default_v2`, `chunk_lp4`, `pool_default_v2`, `microchunk_p1`, and `microchunk_p2`.

## Analyze Outputs

Expected files from `quality_mode=analyze`:

- `cavgs_quality_features.txt`: per-class auto state, hard reject flag, quality cluster, score, manual state, match flag, raw features, and z-normalized features.
- `cavgs_quality_reference_summary.txt`: confusion metrics, AUC, current and best thresholds, and per-feature AUC/separation/weights.
- `cavgs_quality_reference_threshold_scan.txt`: threshold sweep with selected counts and binary metrics.
- `quality_selected_cavgs.mrc` and `quality_rejected_cavgs.mrc`: stacks for visual inspection.

When analyzing a run directory, first run:

```bash
python3 .github/skills/simple-cluster-cavgs-quality/scripts/analyze_quality_run.py /path/to/run
```

Use `--recursive` on a parent directory to summarize many run folders:

```bash
python3 .github/skills/simple-cluster-cavgs-quality/scripts/analyze_quality_run.py --recursive /path/to/test-data
```

Then inspect the reported false positives, false negatives, best thresholds, weak/inverted features, and hard-rejected manual-good classes. For chunk-branch tuning, pay special attention to whether `chunk` mode loses manual-good classes through hard rejection or through the stricter effective threshold.

## Interpretation Checklist

1. Confirm `quality_mode=analyze` was used. In analyze mode, project selection should be unchanged.
2. Compare current threshold metrics with `best_balacc_threshold` and `best_f1_threshold`.
3. Separate false negatives caused by `hard_reject=T` from false negatives caused by score thresholding.
4. Rank features by AUC and robust separation. A high current weight with poor or inverted AUC is a tuning warning.
5. Check whether false positives are accepted by population/resolution alone or by local signal and histogram support.
6. For chunk validation, compare `rejection_type=chunk` against `pool` only when the same project/manual states were used.
7. Use the MRC stacks when numeric diagnostics disagree with visual quality.

## Chunk-Branch Guidance

For `cluster2D_microchunked`, preserve the chunk lifecycle and sentinel behavior described in `.github/skills/simple-microchunk-rejection/SKILL.md`. The live stream path still uses `simple_cluster2D_rejector`; `model_cavgs_rejection` is a separate backend unless an integration task explicitly changes that boundary.

Prefer small local changes:

- share or extract project annotation/state-mapping code before changing rejection policy;
- keep `simple_cluster2D_rejector` temporarily as a comparison backend;
- avoid exposing calibration margins as casual CLI knobs;
- introduce stage-aware quality profiles only if validation data clearly requires pass/ref/match-specific behavior.
