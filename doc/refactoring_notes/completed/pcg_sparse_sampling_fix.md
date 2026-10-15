# PCG on sparse or uneven view sampling: the preconditioner floods unmeasured Fourier voxels

Status 2026-10-15: the fix removes the defect on the FePt tilt series. After two
iterations the PCG map holds 0.003% of its power beyond the Nyquist sphere (13%
before, gridding 0.006%), has no centre spike (centre maximum 311 against a 99.9th
percentile of 319), and its Nyquist shell carries 0.65 times the power of the shells
inside it with no even/odd correlation (13 times and 0.92 before). The half-map FSC
reaches 0.5 at 3.41 Å after two iterations and 3.31 Å after five (9.2 Å before;
gridding 3.4 Å). Against an independent filtered back-projection (band 8-1.5 Å) PCG
now correlates better than gridding: 0.911 and 0.919 against 0.896. The
sparse-sampling test has not been run yet.

## What goes wrong

The FePt tilt series (34 views per half, tilts from -65 to +64 degrees) reconstructs
well with gridding (FSC 0.5 at 3.4 Å). PCG gives 6.6 Å after one iteration and 9.2 Å
after two. The PCG map has a checkerboard spike at the centre, about -3100 against
atom peaks of about 90, and 10-13% of its power lies beyond the Nyquist sphere
(gridding: 0.01%). That excess sits in a thin layer at the Nyquist shell, near the
Nyquist faces in the plane of the tilt. It is the same in both half maps (correlation
0.9997 on the Nyquist plane), so the solver makes it from the low-resolution content
both halves share, not from noise. Low-resolution voxels between the measured planes
and in the missing wedge also hold 3-7 times the gridding power.

## Cause

The first PCG step is `z0 = M^-1 b`, where `M = rho + floor` and the floor is 1% of
the shell-mean density of the measured voxels. Before the Fourier divide, `b` is cut
by the hard solve support (`window > 0`) and by the crop to the native box. The cut
spreads `b` into Fourier voxels that no view measured. There `M` is the floor alone,
so the spread is amplified about 100 times. With few views most voxels of the outer
shells are unmeasured, so the gain is largest near Nyquist. The support edge is a
sphere, and a sphere filtered by a gain that depends on the radius alone focuses on
the centre, which gives the centre spike. The same gain floods the gaps and the
missing wedge at low resolution. Iterations repair this only slowly: the residual
rises at iteration 2 on FePt and in the emulation.

The emulation follows the production chain: a 3-point KB stencil on the 2x padded
lattice with periodic wrap, the full data disk including the Nyquist bins,
deapodization, the hard support with a soft output window, and the matrix-free
operator with full CG. It reproduces the FePt symptoms (spike, Nyquist layer,
residual 0.55 then 0.57). Correlation with the truth after the production budget of
two iterations, at box 64:

| View set | Gridding | PCG now | PCG with the fix |
| --- | --- | --- | --- |
| Single-axis tilt, 34 views | 0.816 | 0.477 | 0.812 |
| 60 random views | 0.893 | 0.700 | 0.890 |
| Preferred orientation, 300 views within 20 degrees | 0.587 | 0.125 | 0.573 |
| 400 random views (dense control) | 0.937 | 0.929 | 0.938 |

The power beyond Nyquist falls from 17% to 0.35% on the tilt series (gridding 0.3%).
Box 32 shows the same picture. So the defect is not specific to tomography. It also
affects cryo-EM whenever the outer shells are sparsely sampled: few views, such as
class-average starts, or preferred orientation. Dense cryo-EM data are not affected.

## Fix

A voxel whose density is below 10% of its shell's mean measured density
(`PCG_RHO_SAMPLED_FRAC`) is treated as unmeasured. Its data term leaves `M`, so
`M^-1` is `1/P_tau` with the ML prior and 0 without it. This is the rule that already
holds for modes beyond the Nyquist sphere. The floor stays for measured voxels.
Thresholds from 0.03 to 0.3 all work; 0.1 is in the middle. With 400 random views
only 0.3% of the voxels fall below it, so dense data are unchanged. The closed-form
regularized map then shrinks unmeasured voxels to zero, where it now shrinks them to
`floor/(floor+P_tau)`, which is close to zero anyway.

Alternatives tried in the emulation:

- A floor equal to the full shell mean also works, but it weakens the
  preconditioning of the measured voxels.
- Starting from the gridding solution works, but it breaks the policy that every
  solve starts from zero.
- The exact diagonal of `P H P` (the density blurred by the support's power
  spectrum) helps only from the second iteration on.

The FLEX PCG (`simple_flex_pca_pcg.f90`) builds its floor the same way. It is not
changed here.

## Shipped ML map: noise left at Nyquist

On dense bgal data (`rec3D_backends`, d2, about 2300 particles per half) the fix
changes nothing and the gate passes. The shipped ML map still had a layer of power
at Nyquist: shells 127 to 130 held 36 to 227 times the power of shells 116 to 125,
with zero even/odd correlation, so it is noise, not solver-made structure. The
base pair holds 2.8% of its power beyond the Nyquist sphere, noise that the hard
support leaks there and that no data term removes. The closed form scaled every
voxel beyond the data sphere by 1, so that noise survived while the shells inside
were shrunk about a hundredfold, and the native crop spread it back into the
outer shells. Beyond the data sphere the factor is now 0.

## Test

The `pcg_sparse_sampling` sub-suite of `lib_reconstruction`
(`simple_pcg_sparse_tester.f90`) uses an atomic phantom at box 32 whose atoms have
strong Nyquist content. It reconstructs three view sets: a single-axis tilt, a
preferred orientation, and a dense control. Each goes through the production path
(raw accumulation, kernel operator, sphere support, deapodization, two fixed
iterations) and through production gridding. On the tilt and preferred-orientation
sets, PCG must correlate with the truth within 0.05 of gridding and keep the power
beyond Nyquist under 3%. On the dense control it must stay within 0.02 of gridding.
In the emulation, the current code fails the first two sets and the fix passes all
three.

## Next

1. Build and run `simple_test_exec test=lib_reconstruction suite=pcg_sparse_sampling`,
   then `test=pcg_recon` and `test=rec3D_backends` for regressions.
2. Run one cryo-EM refine3D regression. Re-check the PfCRT gridding-versus-PCG case,
   and any `PCG BEYOND-BAND EXCESS` warnings seen before, against this cause.
