# Automatic Gain-Reference Flip Detection

Automatic gain flip is an opt-in startup calibration in the streaming
preprocessing workflow. It determines whether the supplied gain reference
should be used unchanged or flipped along the X axis, Y axis, or both axes.
It runs once before normal movie preprocessing; it is not repeated for every
movie.

## Activation

The stream UI exposes the following gain-processing modes:

- `none` (default)
- `flip_auto`
- `flip_x`
- `flip_y`
- `flip_xy`
- `generate`

Selecting `flip_auto` and supplying a gain reference enables automatic
orientation detection. The default remains `none`, so existing streams do not
run this calibration unless it is explicitly requested.

## Movie collection

The preprocessing stage creates a dedicated movie watcher for gain analysis
and waits for batches of 10 usable movies. The normal preprocessing watcher is
a separate object, so movies sampled during gain analysis remain available for
subsequent stream processing.

For every batch, the gain-analysis helper reads each movie frame and adds it to
a batch sum. The analyzer owns the cumulative sum and cumulative movie and
frame counts across batches.

The streaming caller supplies batches of 10 movies and stops after convergence
or after eight batches, so at most 80 movies are analyzed.

## Orientation scoring

The analyzer constructs four versions of the supplied gain reference:

1. unchanged
2. flipped along X
3. flipped along Y
4. flipped along both X and Y

The gain candidates are low-pass filtered to 200 A. At each analysis point,
the cumulative raw-frame sum is divided by the cumulative frame count and
low-pass filtered to the same resolution. The averaged movie image is then
correlated with each candidate.

The candidate with the **most negative correlation**, meaning the minimum
correlation value, is selected. This reflects the inverse relationship between
the raw detector-response pattern and the multiplicative gain-correction
reference.

The separation between the best and second-best candidates is recorded as

```text
delta_r = second_best_correlation - best_correlation
```

The current average and selected candidate are written as
`average_all_movies.mrc` and `best_gainref_by_corr.mrc` for inspection.

## Convergence

The reusable analyzer is eligible to run first at 10 cumulative movies and
then every 5 movies. Because the streaming caller provides 10-movie batches,
the stream evaluates at 10, 20, 30, and subsequent multiples of 10 movies.

Best correlation and `delta_r` are exponentially smoothed with an update
weight of 0.4. A stable update requires:

- the same best orientation as the previous analysis, or a near-tie with
  `delta_r <= 0.01`;
- smoothed best-correlation drift no greater than `0.003`; and
- smoothed `delta_r` drift no greater than `0.003`.

Two consecutive stable updates are required for convergence. Consequently,
the earliest streaming convergence is at 30 movies. If convergence is not
reached by 80 movies, the most recently selected best orientation is still
used.

## Applying the result

The analyzer returns an orientation index, which the stream maps to `no`, `x`,
`y`, or `xy`. For a selected flip, `flip_gain` reads the original, unfiltered
gain reference, flips it, and writes a new MRC file such as
`gainref_flipX.mrc`. It then replaces the `gainref` command-line value with the
absolute path to that file.

The low-pass-filtered `best_gainref_by_corr.mrc` analysis artifact is not used
for correction. Motion correction uses the newly flipped, full-resolution gain
reference and multiplies each movie frame by it.

## EER handling

EER movies are not supported by the current automatic orientation-analysis
path. The streaming collection loop explicitly skips them:

```fortran
if( fname2format(gain_movies(imov)) == 'K' ) cycle
```

This check is in
`src/main/stream/simple_stream_p01_preprocess_new.f90` in
`auto_detect_gain_flip`. Format code `K` identifies EER.

The frame-summing helper also contains a defensive rejection in case an EER
file reaches it:

```fortran
select case(fname2format(movie_fnames(imov)))
case('K')
    THROW_HARD('EER movies are not supported by this test: '//movie_fnames(imov)%to_char())
case DEFAULT
    call find_ldim_nptcls(movie_fnames(imov), ldim_cur, nframes)
end select
```

An EER-only stream therefore cannot provide a valid 10-movie analysis batch.
After the wait limit is reached, no orientation has been selected and the
stream falls back to no flip.

## Other fallback and wait behavior

Automatic detection defaults to no flip when the gain-reference path is empty,
the gain-reference file does not exist, or the expected movie folders cannot
be detected. While collecting movies, the stage allows 180 wait cycles of two
seconds for a complete non-EER batch, approximately six minutes. Input files
must also be untouched for the watcher's 60-second reporting interval before
they are considered ready.

If at least one analysis has run but formal convergence is not reached before
the batch or wait limit, the last best-scoring orientation is used.

## Source locations

- Stream option definition:
  `src/main/ui/simple/simple_ui_stream.f90`
- Stream collection, selection, EER skip, and fallback:
  `src/main/stream/simple_stream_p01_preprocess_new.f90`
- Frame summation and defensive EER rejection:
  `src/main/motion/simple_motion_gain_helpers.f90`
- Candidate scoring and convergence:
  `src/main/motion/simple_motion_gain_analysis.f90`
- Gain materialization and command-line update:
  `src/main/motion/simple_motion_correct_utils.f90`
- Motion-correction application:
  `src/main/motion/simple_motion_correct_utils.f90`, in `correct_gain`
