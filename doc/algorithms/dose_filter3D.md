# Dose-Series Accumulation

## Problem

Every frame of a movie adds dose, and dose fades the signal: radiation damage first, and in practice also
residual motion that the motion model does not remove. Grant & Grigorieff measured how fast on one specimen
and summarised it as a critical exposure per spatial frequency; every dose weighting in common use applies
that one curve to every dataset. `dose_filter3D` measures the fading of the data at hand instead, from
half-map pairs reconstructed from consecutive windows of fractions, and accumulates the windows under it.
`reextract` can then weight the movie fractions of every particle by the same measured curve.

## Data

`fractionate_model` reconstructs state 1 from blocks of `stepf` fractions with the refined poses held
fixed and writes `dose_series.txt`: per window its fractions, its physical dose bounds and its unfiltered,
unmasked even and odd half maps ([fractionate_model](fractionate_model.md), dose-series export). The
manifest also carries the box, sampling, voltage and dose per fraction, so the program needs only the
manifest, `smpd` and `mskdiam`:

```text
simple_exec prg=dose_filter3D filetab=dose_series.txt smpd=<A> mskdiam=<A> \
    [hp=20] [lp=<A>] [bfac_model=line|free] [zerodose=no|yes] [outvol=dose_series.mrc] nthr=16
```

An overlapping last window replaces its predecessor; at least two disjoint windows are required. Doses are
in the scaled coordinate of `calc_dose_weights`: physical dose divided by 0.8 at 200 kV and by 0.64 at
100 kV.

## Model

Window `w` sees the structure attenuated, relative to window 1, by

```text
a_w(k) = c_w exp(-B_w k^2 / 4),      c_1 = 1, B_1 = 0,
```

a frequency-flat contrast `c_w` and a B-factor `B_w`. Across windows both follow straight lines in dose,

```text
B(D) = B0 + beta D,      ln c(D) = ln c0 - gamma D,
```

fitted over windows 2 to W, with an excess scatter `tau` about each. The weights use each window's own
values, so a fading that is not linear in dose is followed rather than forced onto the lines.

## Algorithm

1. **Statistics.** Both halves of every window are masked by the soft sphere of `mskdiam`, Fourier
   transformed and rescaled to unit mean power. Per Fourier shell, one pass forms the `2W x 2W` matrix of
   products of all halves, the noise variance of each half from the half difference, and the variance of
   the cross-half products.

2. **Band.** From `hp` (20 Å by default) to `lp`, or when `lp` is not given to the last shell where any
   window's own even/odd FSC exceeds 0.1. At least six shells.

3. **Attenuation, model-free.** Per shell, a rank-one weighted least-squares fit of the symmetrised
   cross-half products gives `a_w(s)` relative to window 1, with its full covariance between windows (a
   sandwich estimate: each window's noise projected on the structure is shared by every pair holding it).
   The independent count of a shell under the mask is the plane formula
   `f = 3 (sum M^2)^2 / (N^2 sum_axes sum_plane (sum_perp M^2)^2)`.

4. **One B-factor per window.** A Guinier fit of `c_w exp(-B_w k^2/4)` to `a_w(s)` over the band,
   started by a log-domain regression and refined in the linear domain by Levenberg-Marquardt with step
   acceptance. Errors come from the curvature, inflated by the correlation between shells under the mask:

   ```text
   F = (1/3) sum_axes  n_a sum M^4 / sum_perp (sum_a M^2)^2,
   ```

   1 without a mask. A window whose fit fails or whose contrast does not clear five of its errors is
   rejected: it leaves the lines and the model's weights.

5. **The lines.** Weighted lines of `B_w` and `ln c_w` in dose, in two passes: weights `1/se^2`, then
   `1/(se^2 + tau^2)` with `tau^2` the DerSimonian-Laird estimate under the trace of a weighted line. With
   `bfac_model=line` each `B_w` is shrunk towards the line,
   `B*_w = (B_w/se^2 + B(D_w)/tau^2) / (1/se^2 + 1/tau^2)`; with `free`, `B*_w = B_w`.

6. **Weights and accumulation.** Per shell, the matched filter
   `omega_w ∝ c_w exp(-B*_w k^2/4) / sigma_w^2`, normalised to sum one, applied to the halves as given
   (unmasked). The even accumulation uses only even halves and the odd only odd, so they remain a
   gold-standard pair. Past the band the edge shells' weights continue. The same is done under Grant &
   Grigorieff's attenuation (the control) and with equal weights (the flat mean).

7. **Readouts.**
   - FSC 0.5 and 0.143 of the model, the control, the flat mean and window 1 alone, unmasked and under the
     soft sphere of `mskdiam` with zero background, as reconstruct3D masks before its FSC.
   - The cross-half prediction error of each scheme: each window's odd half predicted from the even
     halves of all windows under the scheme's attenuation and vice versa, whitened by the window's noise
     (1 is a noise-free prediction, 2 the window alone). The model-minus-control difference is descriptive:
     the attenuations are fitted on both halves.
   - The critical exposure per shell with no model, `ln a_w(s)` regressed through the origin on
     `D_w - D_1`, `N_e = -1/(2 slope)`; beside it the lines' `1/(2 gamma + beta k^2/2)`, the per-window
     model through the same estimator, and Grant & Grigorieff's.
   - The critical exposure between neighbouring windows at 20, 10, 8, 6 and 4 Å; flat when the fading is
     exponential in dose.
   - The FSC of every pair of windows, the even half of one against the odd half of the other.
   - With `zerodose=yes`, the accumulated map divided per shell by its residual attenuation and continued
     along the B-factor line to zero dose.

Model values are reported against the sharpest fitted window (least `B*`): early beam-induced motion can
leave window 1 blurrier than the next, whose `B` would otherwise read negative. The shift changes no ratio
between windows, no weight and no map.

## Re-extraction under the measured curve

```text
simple_exec prg=reextract projfile=<project> extractfrommov=yes dose_weighting=measured \
    dose_curve=<outvol>_curve.txt nparts=<n> nthr=<n>
```

Fraction `f` of a movie with physical dose `d` per fraction sits at the scaled centre dose
`x_f = (f - 1/2) d / s_kV`. Between fitted windows `ln c` and `B*` are interpolated linearly in dose; before
the first window the first segment continues; past the last the lines' slopes continue, clamped at zero so
an unmeasured fraction never gains on the last measured one. Each frame's Fourier transform is weighted by

```text
omega_f(k) = c_f exp(-B_f k^2/4) / sqrt( sum_g c_g^2 exp(-B_g k^2/2) ),
```

unit squared sum at every frequency, so the summed noise spectrum stays flat, as in the Grant & Grigorieff
mode. The curve was measured on plain sums of fractions, so it also replaces the motion-correction frame
weights. The run stops if a micrograph's voltage or dose per fraction differs from the curve's, or if it has
no motion model (the integrated micrograph carries Grant & Grigorieff's weights).

## Outputs

`<outvol>` is `outvol` without its extension.

| File | Contents |
| --- | --- |
| `<outvol>_even.mrc`, `_odd.mrc`, `.mrc` | The model's accumulated halves and their mean, the delivered map: unfiltered, unsharpened, unmasked |
| `_fsc.bin`, `_fsc_masked.bin` | The FSC of the model's half pair, unmasked and masked |
| `_control_even.mrc`, `_control_odd.mrc`, `_control.mrc`, `_control_fsc.bin`, `_control_fsc_masked.bin` | The same under Grant & Grigorieff's weights |
| `_fsc.txt` | Per shell: the model's and the control's FSC, unmasked and masked |
| `_window_fsc.txt` | Per shell: the FSC of every pair of windows; `fsc_w_w` is window `w`'s own |
| `_zerodose.mrc` | With `zerodose=yes`: the zero-dose extrapolation |
| `_model.txt` | Per window: dose, `c_w`, `B_w` with its error, the line's `B(D_w)`, `B*_w`, residual `z_w`, Guinier rms, `used`; the B-factor and contrast lines |
| `_curve.txt` | The dose curve `reextract` reads |
| `_attenuation.txt` | Per shell: every window's model-free attenuation and its error |
| `_weights.txt`, `_weights_res.txt` | The model's weights per shell, and at 20, 10, 8, 6 and 4 Å beside Grant & Grigorieff's |
| `_cv.txt` | Per shell: the three cross-half prediction errors and the model-minus-control difference |
| `_decay.txt`, `_decay_steps.txt` | The critical exposures, per shell and between neighbouring windows |
| `_decay.png`, `_weights_res.png`, `_plots.gp`, `_plots_gnuplot.log` | Fading per window, critical exposure per shell and between windows, and the weights, each against Grant & Grigorieff; the gnuplot script and its log |

## Rationale

The windows share particles and poses, so they differ only in dose: whatever separates them is the fading,
measured on the specimen and the processing at hand. Ratios of cross-half products cancel the structure's
power and keep the noise of each half independent of the other's. A B-factor per window is the fewest
parameters that follow a fading of any shape in dose; the lines summarise it, extrapolate it past the last
window and, under `bfac_model=line`, steady windows measured with large errors.

Measured on the data, the fading need not follow Grant & Grigorieff. On β-galactosidase and apoferritin it
accelerates with dose, mostly in `B`, as an uncorrected displacement growing linearly with dose would; the
per-window weights follow it, so a curve measured before per-particle polishing does not hold after it.
With equal-width windows, averaging over a window's fractions cancels against window 1 when the fading is
exponential in dose, so window centres stand for fractions; under curvature the bias grows with the square
of the window width.

## Implementation

- CLI: `src/main/ui/simple/simple_ui_filter.f90` (`dose_filter3D`), `src/main/ui/simple/simple_ui_preproc.f90`
  (`reextract`: `extractfrommov`, `dose_weighting`, `dose_curve`).
- Entry point: `src/main/commanders/simple/simple_commanders_resolest.f90`, `exec_dose_filter3D`; plotting
  constants in `src/defs/simple_defs_dose.f90`.
- Manifest: `src/main/volume/simple_dose_series_manifest.f90`.
- Statistics, fits, weights, readouts and the curve export: `src/utils/filter/simple_dose_series.f90`.
- The curve file and its evaluation per fraction: `src/utils/filter/simple_dose_curve.f90`.
- Fraction weights: `image%apply_frame_bfactor_weights` (`src/main/image/simple_image_ops.f90`), applied in
  `ptcl_extractor%init_mov` (`src/main/preprocess/simple_particle_extractor.f90`) from
  `src/main/strategies/parallelization/simple_reextract_strategy.f90`.
- Tests: `unit_numerics suite=dose_series` and `suite=dose_curve`, `unit_image suite=image`,
  `unit_motion suite=particle_extractor`, `lib_reconstruction suite=dose_filter3d`; the condition matrix
  `scripts/dose_filter3D_matrix.sh` with `scripts/make_dose_series.py`.
- Design, equations and decisions: `doc/implementation_notes/planned/dose_filter3D.md`.
