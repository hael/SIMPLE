# FSC Area Score Policy

## Scope

This document defines how SIMPLE calculates the `fsc_area_score` diagnostic.
The implementation is a native half-map conical FSC area ratio, inspired by the
public cFAR description used in orientation-diagnostics tools. It is not a
claim of bitwise compatibility with any external package.

The diagnostic is implemented in Fortran. A Python plotting helper can render
the native CSV outputs in a CryoSPARC-like summary style, but SIMPLE does not
have a CryoSPARC `csfsc.csv` reader or CSV validation mode.

## Public interface

The program is exposed as:

```bash
simple_exec prg=fsc_area_score vol1=odd.mrc vol2=even.mrc smpd=<A/px>
```

The relevant controls are:

- `vol1`: first half map, conventionally the odd map
- `vol2`: second half map, conventionally the even map
- `smpd`: sampling distance in Angstrom per pixel
- `nspace`: number of cone axes, default `256`
- `athres`: cone half-angle in degrees, default `20`
- `lplim_crit`: FSC crossing threshold, default `0.143`
- `automsk`: envelope-mask control; `yes` and `tight` use automasking
- `mskdiam`: spherical-mask diameter when `automsk=no`
- `fbody`: output file body, default `fsc_area_score`
- `nthr`: thread count through the ordinary SIMPLE parameter path

## Inputs and masking

The input maps must be same-size cubic half maps. The commander reads `vol1`
and `vol2` into SIMPLE image objects using the provided `smpd`.

Masking follows the same broad policy as ordinary half-map FSC estimation:

1. If `automsk != no`, SIMPLE builds an automask from the half maps. The
   `tight` value requests the tight automask path. Both half maps have
   envelope background zeroed and are multiplied by the automask.
2. If `automsk=no`, both half maps receive the ordinary soft spherical mask
   from `mskdiam`.

The maps are Fourier transformed after masking. No explicit FFT normalization
is applied because the scale cancels in the FSC ratio.

## Direction sampling

SIMPLE samples cone axes deterministically with a Fibonacci-sphere sequence.
The number of axes is `nspace`. Cone axes are Fourier-space axes; they are not
particle viewing directions.

Each axis is treated as an unoriented line through the Fourier origin. A
Fourier sample with unit vector `k_hat` belongs to a cone around axis `d` when:

```text
abs(dot(k_hat, d)) >= cos(athres)
```

This is an antipodal double-cone selection, which is appropriate for real
half maps and Hermitian Fourier data.

## Fourier shell convention

The conical FSC loop uses SIMPLE's standard nonredundant Fourier half-space,
`fit%loop_lims(2)`, and resolves stored Fourier coefficients through
`fit%comp_addr_phys`.

The shell index is:

```text
shell = nint(sqrt(h*h + k*k + l*l))
```

The DC shell and shells outside the image filter size are skipped. This matches
the shell-binning convention used by the ordinary SIMPLE FSC calculation.

## Conical FSC calculation

For each cone axis and shell, SIMPLE accumulates over all included Fourier
samples:

```text
numerator = sum(F1 * conjg(F2))
power1    = sum(abs(F1)**2)
power2    = sum(abs(F2)**2)
count     = number of included Fourier samples
```

The conical FSC value is:

```text
cFSC(shell, axis) = real(numerator) / sqrt(power1 * power2)
```

The value is reported only when `count >= min_count` and both power sums are
positive. The current commander uses `min_count=1`. Otherwise the internal FSC
value is zero and the output curve table marks that shell/axis entry as `nan`.

## Area scoring

Each conical FSC curve is reduced to a weighted area under the curve, `wAuC`.
SIMPLE walks shells in increasing Fourier index and stops at the first valid
shell where:

```text
cFSC(shell, axis) <= lplim_crit
```

The crossing shell itself is not included in the area. If a curve never crosses
the threshold, all valid shells are included. Invalid sparse shells are skipped.

The current area weight is proportional to Fourier shell surface area:

```text
weight(shell) = shell**2
wAuC(axis)   = sum(cFSC(shell, axis) * shell**2)
```

Constant factors such as `4*pi` are not included because they cancel in the
final ratio.

## Reported score

The final FSC area score is the conical FSC area ratio:

```text
cFAR = min_axis(wAuC(axis)) / max_axis(wAuC(axis))
```

If the maximum directional area is not positive, SIMPLE reports `cFAR=0`.

The log reports the number of axes, cone half-angle, FSC threshold, minimum and
maximum weighted areas, and the final `cFAR`.

## Outputs

For an output file body `fbody`, SIMPLE writes:

- `fbody_summary.txt`: settings, `cFAR`, min/max `wAuC`, and min/max directions
- `fbody_curves.csv`: wave number, resolution, and one cFSC curve per axis
- `fbody_directions.csv`: axis coordinates, `wAuC`, threshold crossing, and
  number of included shells per axis
- `fbody_plot.png`: backend CPlot2D rendering of the directional cFSC summary

These CSV files are diagnostic outputs from the native calculation. SIMPLE does
not use them as inputs and does not validate external CSV files.

The backend plot is generated unconditionally from the in-memory
`fsc_area_score_result`. It includes the directional mean curve, +/- one
standard-deviation envelope, min/max envelope, FSC threshold, and relative
occurrence of directional threshold crossings.

The optional Python helper can be used to reproduce or customize the PNG from
the CSV diagnostics:

```bash
python3 scripts/plot_fsc_area_score.py fsc_area_score
```

It reads `fsc_area_score_summary.txt`, `fsc_area_score_curves.csv`, and
`fsc_area_score_directions.csv`, then writes `fsc_area_score_plot.png`.

## Interpretation

The score is intended to summarize directional anisotropy in half-map signal.
Lower `cFAR` means the weakest cone has much less weighted FSC area than the
strongest cone. The direction associated with a cone is a Fourier-space signal
axis, not a direct particle-orientation assignment.

Because this is a finite-shell, finite-cone diagnostic, the value depends on
masking, cone angle, directional sampling density, box size, and the threshold.
Comparisons are most meaningful when those controls are held fixed.
