# CTF Modeling and Estimation

## Problem

Estimate the microscope contrast transfer function (CTF) of each corrected
micrograph so downstream picking, class averaging, matching, and reconstruction
can interpret phase reversals and frequency-dependent attenuation.

## Microscope model

For squared spatial frequency `s^2` and azimuth `theta`, SIMPLE uses the
astigmatic defocus

```text
d(theta) = 1/2 [dfx + dfy + cos(2(theta-angast))(dfx-dfy)].
```

With electron wavelength `lambda`, spherical aberration `Cs`, additive phase
shift `phi`, and amplitude-contrast phase `alpha`, the phase aberration is

```text
chi(s,theta) = pi*lambda*s^2 [d(theta) - 1/2 lambda^2*s^2*Cs] + phi
CTF(s,theta) = sin(chi(s,theta) + alpha).
```

Underfocus is positive. `phi` is always supplied explicitly and canonicalized
to `[0,2*pi)`. The 1D initializer and final continuous power-spectrum model use
`|CTF|`, so phase windows spanning approximately `pi` admit sign-opposite
equivalent solutions; SIMPLE warns when a fitted phase lies at such a window
edge. The intermediate global 2D search uses the signed positive-part model
described below only to find continuous starts.

## Spectrum construction

The micrograph is tiled. Each tile is transformed to a power spectrum, the
central Fourier cross is damped, and a slowly varying local background is
subtracted. Averaging the prepared tiles gives the global 2D spectrum. A
rotational average supplies a cheaper 1D spectrum for initialization.

Only the requested high/low-resolution annulus participates. The prepared data
are normalized in that annulus so the fitting objective is a correlation, not
an arbitrary power-scale comparison.

## Coarse search

The initializer evaluates 200 mean-defocus values over `[dfmin,dfmax]` on the
rotational average. If phase fitting is enabled, it also evaluates the
requested phase grid. For model samples `t_j=|CTF_j|` and normalized spectrum
samples `p_j`, the 1D score is the centered cosine correlation; the minimized
cost is its negative.

The defocus range is divided into six bins. The best initializer from each bin
is promoted to a full 2D fit, preventing one narrow coarse basin from being the
only continuous start.

## Astigmatic refinement

Each promoted start first undergoes bounded differential-evolution search over
`(dfx,dfy,angast)` and, when enabled, `phi`. This coarse 2D route clamps the
signed CTF model below at a small positive floor. Its objective is

```text
cost = -corr(P, max(CTF,epsilon))
       + ((dfx-dfy)/astigtol)^2 / (2N),
```

where `N` is the number of Fourier pixels in the fitting annulus. The penalty
regularizes implausibly large astigmatism without imposing equality.

A bounded analytic-gradient refinement then optimizes the same parameters
locally against the absolute CTF model `|CTF|`, including the quotient-rule
derivative of its normalized correlation. The best final correlation among the
six promoted solutions becomes the global estimate. Axis ordering and angle
are canonicalized before publication.

## Patch CTF

With `ctfpatch=yes`, overlapping local spectra are formed as Gaussian-weighted
mixtures of nearby tiles. Each patch fits local defocus near the global
solution while holding the global microscope constants and phase policy. A
10-term spatial polynomial is fit to the patch estimates, yielding a smooth
defocus surface rather than independent noisy patch values.

## Outputs and diagnostics

The micrograph record receives `dfx`, `dfy`, mean defocus, astigmatism angle,
phase shift, fit correlation, CTF resolution, ice fraction, and astigmatism
magnitude. Diagnostic images show the observed and fitted spectra; patch mode
also writes the polynomial document used to evaluate local CTF values.

## Implementation

- Model: `src/main/ctf/simple_ctf.f90`.
- Objectives: `src/main/ctf/simple_ctf_estimate_cost.f90`.
- Spectrum, search, patch fit, and diagnostics:
  `src/main/ctf/simple_ctf_estimate_fit.f90`.
- Per-micrograph lifecycle: `src/main/ctf/simple_ctf_estimate_iter.f90`.
