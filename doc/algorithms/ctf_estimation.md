# CTF Modeling and Estimation

## Problem

Estimate, for each micrograph, the defocus pair `(dfx, dfy)`, the
astigmatism angle, and optionally the additional phase shift of a phase
plate, from the oscillation pattern (Thon rings) in the micrograph's power
spectrum. Every downstream comparison of images with references depends on
these parameters through the contrast transfer function.

## Model

The astigmatic defocus in direction `theta` is

```text
d(theta) = 1/2 [ dfx + dfy + (dfx - dfy) cos 2(theta - angast) ],
```

and with wavelength `lambda`, spherical aberration `Cs`, phase shift `phi`,
and amplitude-contrast phase `alpha`, the phase aberration and CTF at
spatial frequency `s` are

```text
chi(s, theta) = pi lambda s^2 [ d(theta) - 1/2 lambda^2 s^2 Cs ] + phi,
CTF(s, theta) = sin( chi(s, theta) + alpha ).
```

Underfocus is positive. The power spectrum of the micrograph is proportional
to `CTF^2` times the specimen and envelope power, so the fitting target is
the modulation `|CTF|`. Because `|sin|` is invariant to a sign flip of its
argument, phase windows spanning about `pi` admit sign-opposite equivalent
solutions; a warning is issued when a fitted phase lies at such a window
edge.

## Spectrum construction

The micrograph is tiled with `512`-pixel tiles at half-tile stride (50
percent overlap). Each tile is normalized, edge-averaged to zero, Fourier
transformed, and its power spectrum accumulated. The mean tile spectrum
then has its central cross damped and a slowly varying background removed:
the background at each pixel is the mean over a square window of side
`~1.4 x` the high-pass radius, computed with a summed-area table. What
remains is the oscillatory part. The spectrum is normalized to zero mean and
unit variance inside the annulus `[hp, lp]` so that the fitting objective is
a correlation rather than a power comparison. A rotational average gives a
cheaper 1D spectrum for initialization.

## Algorithm

1. **Coarse 1D search.** Evaluate 200 mean-defocus values on `[dfmin,
   dfmax]` (default 0.2 to 5 um) against the rotational average, and if phase
   fitting is enabled, a grid of phase values at each. The score is the
   centered cosine correlation between `|CTF|` and the normalized spectrum.
2. **Multi-start selection.** Divide the defocus range into six bins and
   promote the best sample from each bin. A single best start would often
   sit in a narrow basin created by ring aliasing; six starts from distinct
   defocus regions make the global optimum reachable.
3. **Global 2D search.** From each start, differential evolution over
   `(dfx, dfy, angast[, phi])` with population 136, up to 400 generations,
   within `+/- 2 max(astigtol, df_step)` of the start in defocus. The
   objective is

   ```text
   cost = - corr(P, |CTF|)  +  ((dfx - dfy) / astigtol)^2 / (2 N),
   ```

   with `N` the number of Fourier pixels in the annulus and `astigtol`
   defaulting to 0.05 um; solutions with `|dfx - dfy| > 2 astigtol` are
   rejected outright. The penalty regularizes implausible astigmatism without
   forcing equality.
4. **Local refinement.** L-BFGS-B on the same parameters with analytic
   gradients of the normalized correlation (quotient rule), angle within
   `+/- 30` degrees and phase within `+/- pi/6` of the start.
5. **Selection.** The best final correlation among the six refined solutions
   is the estimate. Axes are ordered and the angle canonicalized before
   publication.

## Diagnostics

**CTF resolution.** The spectrum is resampled to at most 1.4 A per pixel and
normalized between the third and fourth zeros. Along the mid-astigmatism
direction, two 1D profiles are formed, the observed spectrum and `|CTF|`,
both rank-normalized within each half-period so that amplitude decay does not
dominate. A sliding-window Pearson correlation between them is computed with a
window that widens with the density of extrema. Starting at 10 percent of
Nyquist, the first shell at which the correlation drops below 0.1 after three
shells above it, below 0.5 after three shells above it, or below 0.5 in more
than three of the last five shells, defines the resolution to which the
fitted model explains the observed rings.

**Ice fraction.** When the sampling admits 3.7 A, the ratio of the spectrum
amplitude at the crystalline-ice peak (3.7 A, searched within 10 shells) to
the amplitude of the first CTF maximum, each averaged over 3 shells, measures
ice contamination. If the CTF resolution is better than 3.8 A, half the CTF
peak is subtracted from the ice peak to correct for Thon-ring overlap.

**Astigmatism** is reported as `|dfx - dfy| / mean(dfx, dfy)`.

## Patch CTF

With `ctfpatch=yes`, local spectra are formed at half the tile-grid density
as Gaussian-weighted mixtures of nearby tiles (`w = exp(-(r/box)^2 / 2)`), and
each is fitted for `(dfx, dfy)` alone within 1 um of the global solution,
holding angle and phase. The patch values are then regressed onto a 10-term
polynomial in normalized position,

```text
{ 1, x, x^2, x^3, y, y^2, y^3, xy, xy^2, x^2y },
```

giving a smooth defocus surface that can be evaluated at any particle
coordinate. Tilted or bent specimens produce a defocus gradient of a few
hundred nanometers across a micrograph; the polynomial captures it without
propagating the noise of independent patch fits.

## Implementation

- Model: `src/main/ctf/simple_ctf.f90`.
- Objectives and gradients: `src/main/ctf/simple_ctf_estimate_cost.f90`.
- Spectrum, search, patch fit, diagnostics: `src/main/ctf/simple_ctf_estimate_fit.f90`.
- Per-micrograph driver: `src/main/ctf/simple_ctf_estimate_iter.f90`.
- Differential evolution: `src/main/opt/simple_opt_de.f90`.
- Policy: `doc/policies/phase_shift_ctf_policy.md`.
