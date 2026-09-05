# Automatic Gain-Reference Flip Detection

## Problem

A detector gain reference must be applied in the same orientation as the
movie frames it corrects, but acquisition software and file formats disagree
about axis conventions, so a supplied reference may be flipped in `x`, `y`,
or both relative to the data. Applying a wrongly oriented gain reference
imprints the detector's fixed-pattern response onto every micrograph. The
task is to decide among the four orientations from the first movies of a
session, once, before any of them is processed.

## Model

The raw (uncorrected) frame average of many movies is dominated at low
spatial frequency by the detector's multiplicative response pattern, which is
the reciprocal of the gain reference. At 200 A and coarser, specimen content
averages out across movies while the fixed pattern accumulates. The correctly
oriented gain reference is therefore the candidate whose low-pass version is
most *anti*-correlated with the low-pass raw average:

```text
orientation = argmin_o  corr( LP_200A(raw average),  LP_200A(flip_o(gain)) ),
o in { none, x, y, xy }.
```

Using the minimum rather than the maximum correlation is the essential
detail: the raw frames are proportional to the specimen intensity times the
detector response, and the gain reference is the inverse of that response.

## Algorithm

1. Accumulate the pixelwise sum of all frames of the first 10 usable movies
   (EER movies are skipped because they cannot be summed without
   fractionation) and divide by the frame count.
2. Low-pass the average and the four candidate orientations to 200 A and
   compute the four Pearson correlations. Record the winner and its margin
   `delta_r = second best - best`.
3. Add 10 more movies and repeat. Smooth the best correlation and the margin
   with an exponential moving average of weight 0.4 on the new value.
4. Declare convergence after two consecutive analyses in which the winning
   orientation is unchanged (or the margin is within 0.01, a near-tie) and
   both smoothed statistics have drifted by at most 0.003.
5. Stop at convergence or after 80 movies, using the most recent winner in
   the latter case. The earliest possible convergence is at 30 movies.

If the gain reference is absent, unreadable, or no complete non-EER batch
arrives within the wait limit, no flip is applied.

The selected orientation is then applied to the original, full-resolution
gain reference, which is written as a new file and used for all subsequent
motion correction. The low-passed candidates exist only for the decision.

## Rationale

- Ten movies per analysis and a 200 A cutoff make the statistic insensitive
  to any single specimen; convergence on a smoothed margin, rather than on a
  single winner, protects against a lucky early tie.
- The decision is made once because the orientation is a property of the
  acquisition setup, not of individual movies, and re-deciding per movie
  would only add noise.

## Implementation

- Scoring and convergence: `src/main/motion/simple_motion_gain_analysis.f90`.
- Frame summation: `src/main/motion/simple_motion_gain_helpers.f90`.
- Gain materialization and application: `src/main/motion/simple_motion_correct_utils.f90`.
- Stream integration: `src/main/stream/simple_stream_p01_preprocess_new.f90`.
- Policy: `doc/policies/motion_gain_analysis_policy.md`.
