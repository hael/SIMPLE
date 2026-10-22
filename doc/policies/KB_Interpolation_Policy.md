# Kaiser-Bessel Interpolation Policy

This policy states how SIMPLE interpolates on the oversampled (padded) Fourier
lattice: image polarization, volume reprojection, 3D reconstruction (plane
insertion), class averaging and the continuous Cartesian pose refiner. Read it
before changing a kernel, the padding factor, a gather or a splat, or an
amplitude scale factor.

## 1. Convolution Interpolation Framework

All interpolation operations are performed on the oversampled (padded) Fourier
lattice. Logical Fourier coordinates are converted from native units to padded
units using the oversampling factor pf.

```text
x_pd = pf · x_orig
```

Interpolation weights are evaluated in padded logical units and sampling is
performed using unit increments on the padded grid. The result of interpolation
is written into the native (unpadded) target lattice.

General interpolation formula in padded logical coordinates:

```text
F(x_pd) = Σ_{k_pd ∈ W} w(x_pd − k_pd) · F_pd(k_pd)
```

In 3D, the Kaiser-Bessel kernel remains separable:

```text
w(x,y,z) = w(x) · w(y) · w(z)
```

### Key Parameters

- Oversampling factor pf = 2.
- Kernel width W (wdim): determines support size.
- Spectral parameter β: controls concentration and stability.
- Amplitude scaling: pf² (2D) or pf³ (3D).

### Constants

The padding factor is the compile-time constant `OSMPL_PAD_FAC = 2`, and the
kernel's oversampling parameter α is `KBALPHA = real(OSMPL_PAD_FAC)`
(`src/defs/simple_defs.f90`). The kernel half width is `KBWINSZ = 1.5`, so
W = 2·`KBWINSZ` = 3 padded units and the stencil has three taps per axis
(`kb_windim`). The spectral parameter is

```text
β = π · sqrt( (W² / α²) · (α − 1/2)² − 0.8 )
```

evaluated as the constant `KB_BETA_KB15_A2` for the default W and α
(`simple_kbinterpol.f90`).

The padded box is `boxpd = 2·round2even(α·box/2)`. Class sums at a cropped box
interpolate on `box_croppd = 2·round2even(α·box_crop/2)`, which spans the same
physical extent as `boxpd`, so a Fourier index h denotes the same spatial
frequency on both grids. A spectrum on the plane grid is read at the padded
shell pf·n for native shell n.

## 2. Gather vs. Splat Duality under Oversampled Policy

Projection (gather):

```text
f = Σ w_pd(x_pd − k_pd) · F_pd(k_pd)
```

Reconstruction (splat):

```text
F_native += w_pd(x_pd − k_pd) · f_pd
```

Both operations evaluate weights in padded logical units and operate on padded
data samples. The target grid (polar or Cartesian) remains native.

## 3. Polarization of 2D Particle Images

Particle images are noise-normalized and Fourier transformed after real-space
padding by pf=2. Polar sampling coordinates are converted into padded logical
units prior to interpolation.

```text
x_pd = pf · (r cosθ)
y_pd = pf · (r sinθ)
```

Interpolation is performed directly on the padded FFT lattice.

Amplitude correction:

```text
Scale factor = pf²
```

## 4. Polar Central Section Extraction from 3D Volumes

Volumes are padded in real space, Fourier transformed, and stored in expanded
logical form. For each projection orientation R:

```text
loc_pd = pf · (R · [h,k,0])
```

Interpolation weights are evaluated in padded units at loc_pd. Samples are
gathered from the padded expanded Fourier volume.

Amplitude correction for padded 3D FFT:

```text
Scale factor = pf³
```

## 5. 3D Reconstruction (Plane Insertion)

Fourier planes are padded before insertion. Only the padded samples at
(pf·h, pf·k) are used, so the plane is stored on the native lattice:
F_plane(h,k) holds the padded sample (pf·h, pf·k). Only k ≤ 0 is stored;
Friedel symmetry gives k > 0. The plane's frlims and nyq keep their padded
meaning.

```text
comp = pwght · pf² · F_plane(h,k)
```

Non-uniform sampling location is computed in native units and converted to
padded units:

```text
loc_pd = pf · (R · [h,k,0])
```

KB weights are evaluated in padded logical units while splat updates are
applied to native expanded volume indices.

```text
Splat update: F_native(window) += comp · w_pd
```

## 6. Class Average Restoration

Class averaging follows the same oversampled-grid policy. Padded FFT samples are
interpolated using padded-grid kernel weights and accumulated onto native
class-average grids.

### 6a. Continuous Cartesian Pose Refinement (refine=cont)

The Cartesian pose refiner (`cartft_calc`) is the one gather that does not
use the KB kernel. It optimizes a continuous pose with analytic derivatives,
and the three-tap KB stencil jumps where a sample's nearest grid point changes
(the window is not zero at its edge, `test_stencil_switch`). Its gather is the
separable quadratic B-spline on the same three taps
(`gather_packed_bspline3_window_grad`): weights and slopes vanish at the taps
that enter or leave, so the value and its derivatives are continuous, and the
weights sum to one at every offset.

So that the model is the polar branch's on grid, the padded reference is
weighted per axis by the KB stencil's real-space response over the B-spline
stencil's (`ref_weight`, applied in `set_ref`). On grid, the B-spline gather of
the weighted transform then equals the KB gather of the unweighted one, and the
particle keeps the KB taper, so the sigma2 of the two representations agree.
Off grid the two kernels differ by their interpolation errors
(`test_matched_projector_boundary`, `simple_cartft_calc_tester.f90`; the gather
itself is guarded by `test_fused_bspline3_gather_oracle` and
`test_stencil_switch`, `simple_cartesian_fourier_tester.f90`).

## 7. Scaling Summary

- Padded 2D FFT → multiply by pf²
- Padded 3D FFT → multiply by pf³
- Kernel geometry defined in padded logical units
- Native target grids remain unpadded
