# Internal Technical Memorandum

Subject: Kaiser-Bessel Convolution Interpolation on Oversampled Fourier Grids

Purpose: Consolidated technical description of the oversampled-grid convolution policy used in image polarization, volume reprojection, 3D reconstruction, and class averaging.

## 1\. Convolution Interpolation Framework

All interpolation operations are performed on the oversampled (padded) Fourier lattice. Logical Fourier coordinates are converted from native units to padded units using the oversampling factor pf.

x_pd = pf · x_orig

Interpolation weights are evaluated in padded logical units and sampling is performed using unit increments on the padded grid. The result of interpolation is written into the native (unpadded) target lattice.

General interpolation formula in padded logical coordinates:

F(x*pd) = Σ*{k_pd ∈ W} w(x_pd − k_pd) · F_pd(k_pd)

In 3D, the Kaiser-Bessel kernel remains separable:

w(x,y,z) = w(x) · w(y) · w(z)

### Key Parameters

- Oversampling factor pf = 2 (current configuration).
- Kernel width W (wdim): determines support size.
- Spectral parameter β: controls concentration and stability.
- Amplitude scaling: pf² (2D) or pf³ (3D).

## 2\. Gather vs. Splat Duality under Oversampled Policy

Projection (gather):

f = Σ w_pd(x_pd − k_pd) · F_pd(k_pd)

Reconstruction (splat):

F_native += w_pd(x_pd − k_pd) · f_pd

Both operations evaluate weights in padded logical units and operate on padded data samples. The target grid (polar or Cartesian) remains native.

## 3\. Polarization of 2D Particle Images

Particle images are noise-normalized and Fourier transformed after real-space padding by pf=2. Polar sampling coordinates are converted into padded logical units prior to interpolation.

x_pd = pf · (r cosθ)

y_pd = pf · (r sinθ)

Interpolation is performed directly on the padded FFT lattice.

Amplitude correction:

Scale factor = pf²

## 4\. Polar Central Section Extraction from 3D Volumes

Volumes are padded in real space, Fourier transformed, and stored in expanded logical form. For each projection orientation R:

loc_pd = pf · (R · \[h,k,0\])

Interpolation weights are evaluated in padded units at loc_pd. Samples are gathered from the padded expanded Fourier volume.

Amplitude correction for padded 3D FFT:

Scale factor = pf³

## 5\. 3D Reconstruction (Plane Insertion)

Fourier planes are padded before insertion. Sampling from the padded plane uses indices (pf·h, pf·k).

comp = pwght · pf² · F_plane(pf·h, pf·k)

Non-uniform sampling location is computed in native units and converted to padded units:

loc_pd = pf · (R · \[h,k,0\])

KB weights are evaluated in padded logical units while splat updates are applied to native expanded volume indices.

Splat update: F_native(window) += comp · w_pd

## 6\. Class Average Restoration

Class averaging follows the same oversampled-grid policy. Padded FFT samples are interpolated using padded-grid kernel weights and accumulated onto native class-average grids.

## 7\. Scaling Summary

- Padded 2D FFT → multiply by pf²
- Padded 3D FFT → multiply by pf³
- Kernel geometry defined in padded logical units
- Native target grids remain unpadded

## 8\. Summary of Policy Change

The previous strided policy evaluated kernels in native logical units and sampled padded data using stride-based parity selection. The updated policy evaluates kernels directly in padded logical units and performs interpolation using unit increments. This provides a conceptually cleaner, symmetric gather/splat formulation and aligns with standard oversampled-grid interpolation practice.