# EMBB abinitio3D NU Filtering Summary

Parsed `600` NU-filtered iterations from `10` repeated abinitio3D runs.

## Final-Iteration Aggregate

- Mean base-bank local resolution across final maps: `15.61 A`
- Mean of final base-bank local-resolution medians: `15.06 A`
- Mean auxiliary assignment before/after Potts smoothing: `12.96%` -> `11.79%` of mask voxels
- Mean Potts smoothing changed voxels in final iterations: `37809`
- Mean discontinuous neighbor-pair rate: `1.16%`
- Mean voxels with any discontinuous neighbor: `10.35%`

Note: SIMPLE reports the NU local-resolution summary over base-bank voxels. Auxiliary ML-regularized assignments are reported separately as a percentage of the support mask.

## Final-Iteration Runs

| Run | Stage | Iter | Mean A | Median A | Finest selected base A | Dominant base LP | Aux % | Disc pair % | Disc voxel % | FSC143 A | Matching LP A |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| R01 | 8 | 129 | 15.69 | 15.06 | 4.0 | 19.7 A (42.6%) | 11.57 | 0.99 | 9.22 | 4.063 | 4.491 |
| R02 | 8 | 151 | 15.64 | 15.06 | 4.0 | 19.7 A (41.7%) | 12.09 | 1.19 | 10.65 | 4.000 | 4.491 |
| R03 | 8 | 151 | 15.48 | 15.06 | 4.0 | 19.7 A (40.0%) | 11.78 | 1.11 | 10.35 | 4.000 | 4.491 |
| R04 | 8 | 134 | 15.65 | 15.06 | 4.0 | 19.7 A (42.2%) | 11.48 | 1.30 | 11.55 | 4.000 | 4.491 |
| R05 | 8 | 143 | 15.45 | 15.06 | 4.0 | 19.7 A (39.2%) | 11.73 | 1.20 | 10.53 | 4.000 | 4.491 |
| R06 | 8 | 136 | 15.55 | 15.06 | 4.0 | 19.7 A (40.7%) | 11.63 | 1.25 | 10.72 | 4.000 | 4.491 |
| R07 | 8 | 151 | 15.52 | 15.06 | 4.0 | 19.7 A (40.1%) | 11.95 | 1.18 | 10.44 | 4.000 | 4.491 |
| R08 | 8 | 151 | 15.62 | 15.06 | 4.0 | 19.7 A (42.1%) | 11.97 | 1.28 | 11.04 | 4.000 | 4.491 |
| R09 | 8 | 151 | 15.61 | 15.06 | 4.0 | 19.7 A (40.4%) | 12.11 | 1.09 | 9.85 | 4.000 | 4.491 |
| R10 | 8 | 129 | 15.92 | 15.06 | 4.0 | 19.7 A (46.1%) | 11.60 | 1.05 | 9.19 | 4.063 | 4.491 |

## Generated Figures

- `nu_validation_dashboard.svg`: final assignment heatmap plus final continuity/local-resolution summary.
- `nu_trajectories.svg`: per-run NU-stage trajectories for median resolution, auxiliary assignment, discontinuity, and Potts changed voxels.

## Parsed Data

- `nu_iteration_metrics.csv`: one row per NU-filtered iteration.
- `nu_bank_assignments.csv`: one row per low-pass bank member per NU-filtered iteration.
- `nu_step_distribution.csv`: one row per LP-step-difference class per NU-filtered iteration.
