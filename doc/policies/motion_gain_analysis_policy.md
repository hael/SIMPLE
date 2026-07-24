# Motion Gain Analysis Policy

## 1. Purpose and Scope

This document records the current policy for SIMPLE's motion gain-orientation
analysis workflow used by test-time gain-flip evaluation.

Scope:

- `simple_motion_gain_analysis.f90`
- `simple_motion_gain_analysis_helpers.f90`
- test wiring through `simple_motion_gain_tester.f90`,
  `simple_test_search_gain_flips.f90`, and unit-test entrypoints

This is a workflow contract document, not a line-by-line implementation map.

## 2. Architectural Policy

The workflow has two distinct responsibilities that must stay separate:

1. movie ingestion and frame summation
2. gain-orientation analysis and convergence tracking

Current ownership:

- `simple_motion_gain_analysis_helpers.f90` owns reading movie stacks and
  returning a summed image plus count metadata.
- `simple_motion_gain_analysis.f90` owns gain-reference variants,
  correlation scoring, convergence metrics, and best-gain output artifacts.

Do not collapse these responsibilities back into one monolithic routine.

## 3. Data Contract

### 3.1 Helper Contract

`read_movies_and_sum_frames(movie_fnames, smpd, sum_img, n_movies, total_frames)`
returns:

- `sum_img`: sum of all frames across all input movies
- `n_movies`: number of movies consumed
- `total_frames`: total number of frames consumed

Contract rules:

- input list must be non-empty
- every movie file must exist
- EER (`fname2format == 'K'`) is rejected in this path
- all movie dimensions must match the first movie dimensions
- `sum_img` is initialized by the helper based on the first movie

### 3.2 Analyzer Contract

`analyze_if_due(self, sum_part_img, part_frames, part_movies, ran_analysis)`
consumes batch-local partials and updates analyzer-owned cumulative state.

- `sum_part_img`: summed image for the current batch/part
- `part_frames`: number of frames represented by `sum_part_img`
- `part_movies`: number of movies represented by `sum_part_img`

Analyzer-owned cumulative state:

- `self%sum_img` accumulates partial sums
- `self%total_frames` accumulates `part_frames`
- `self%total_movies` accumulates `part_movies`

Callers must pass batch-local counts, not global cumulative counts.

## 4. Analysis Scheduling Policy

Analysis eligibility is controlled by cumulative movie count:

- first analysis at `GAIN_FIRST_ANALYSIS_AT` movies (default `10`)
- subsequent analyses every `GAIN_ANALYSIS_STEP` movies (default `5`)

A call that is not due must return with `ran_analysis = .false.`.
A due call with valid `part_frames > 0` must perform analysis and set
`ran_analysis = .true.`.

## 5. Scoring and Convergence Policy

### 5.1 Scoring

- Build average image: `self%sum_img / self%total_frames`
- Apply low-pass filter at `GAIN_LP_RES_A` (current default `200 A`)
- Correlate the averaged movie image with four gain candidates:
  - unchanged
  - flip X
  - flip Y
  - flip XY
- Select best orientation by **most negative** correlation (minimum value)

### 5.2 Convergence

Convergence uses both orientation stability and smoothed metric stability:

- best-index stability (or near-tie by `GAIN_CONV_EPS_DELTA`)
- smoothed best-correlation drift threshold
- smoothed delta-r drift threshold
- streak threshold `GAIN_CONV_STREAK_REQ`

`get_converged()` exposes convergence state to callers.

## 6. Artifact Policy

Analysis writes these workflow artifacts:

- `average_all_movies.mrc`
- `best_gainref_by_corr.mrc`

These artifacts are intentional outputs for inspection and downstream usage in
this workflow. Test code may create and clean them up explicitly.

## 7. Batch Execution Policy

The current batch runner (`simple_test_search_gain_flips`) processes movies in
fixed-size batches (default `10` movies per batch), and for each batch:

1. sum frames for that batch only
2. call analyzer with batch-local `sum_part_img`, `part_frames`, `part_movies`
3. allow analyzer to decide whether analysis is due

This preserves the analyzer contract that cumulative ownership stays internal.

## 8. Testing Policy

`simple_motion_gain_tester.f90` is the canonical unit-test module for this
workflow and must cover at least:

- helper counts and summation correctness
- analyzer due/not-due scheduling behavior
- analyzer cumulative frame/movie bookkeeping
- expected analysis output artifact creation

Unit-test runners should invoke `run_all_motion_gain_tests` alongside existing
`run_all_*_tests` suites.

## 9. Refactor Guardrails

Do not:

- pass cumulative counts into `part_frames`/`part_movies`
- move cumulative ownership out of the analyzer
- merge helper and analyzer responsibilities into one routine
- change best-orientation selection away from minimum correlation without an
  explicit policy update
- silently change schedule constants/meaning without updating this document

When behavior changes are intentional, update this policy in the same change.
