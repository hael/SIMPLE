# refine_tomo: tilt-series refinement for atomic electron tomography in SINGLE

Status: written 2026-10-15, not compiled; the ground-truth test is still to do (see "Status" at the end). Plan agreed with Hans on 2026-10-15. It turns the hand-run FePt kit
(`for_simple/run_fept_refine_cont.sh`) into a SINGLE workflow. On the FePt tilt
series that kit raised level-1 atom recall from 0.35 to 0.39 (their GENFIRE
map: 0.41) by refining the projection poses with `refine3D refine=cont`.

## What the workflow does

`single_exec prg=refine_tomo projfile=<project> mskdiam=<A>` on a project made by
`import_tilt_series`:

1. **Gridding** reconstruction from the imported poses, half maps by index
   parity (the import sets the half-set labels).
2. **Density mask**: `automask3D` on the gridding map, written as a mask
   file. A tilt series measures far less than a cryo-EM data set, so the
   solve support matters more; a density envelope is tighter than the
   `mskdiam` sphere.
3. **First reference**: masked PCG, `reconstruct3D rec_backend=pcg
   pcg_mskfile=<mask>`. The PCG strategy installs `pcg_mskfile` as the solve
   support in `reconstruct3D` and `refine3D` alike.
4. **Refinement**: `refine3D refine=cont`, by default one stage of two
   iterations at 2 A (`nstages=1`, `lpstart=lpstop=2`, `maxits=2`) under the
   standard continuous-pose convergence rule, at the native sampling:
   high-resolution tomography is not downscaled. Longer refinement drifts (see
   below). The low-pass never goes beyond the half-map FSC 0.143 resolution of
   the map it starts from. `nstages > 1` marches the low-pass geometrically from
   `lpstart` to `lpstop`, rebuilding the mask after each stage, with the
   iteration count running on across the stages as in the staged `solve3D`.
5. **Final map and poses**: a masked PCG reconstruction from the final poses,
   the poses written in SIMPLE's convention and in the GENFIRE/RESIRE
   convention. Atom detection stays a separate `detect_atoms` run.

Defaults follow the hand-run kit: objective `cc` (no sigma2), `ml_reg=no`,
`center=no` (the frame of the input poses is kept), `trs=3` pixels,
`athres_cont=5` degrees, five PCG iterations.

## Tunable convergence of refine=cont

The thresholds of the continuous-pose rule are constants in `simple_defs_conv`
(0.5 degrees, 1 A, 90% of the particles). They become `refine3D` parameters
with the constants as defaults, so no existing run changes:

- `cont_conv_rot`: the largest rotation of a stable particle (degrees);
- `cont_conv_shift`: the largest shift of a stable particle (A);
- `cont_conv_frac`: the stable fraction per state and the coverage of the
  stable streak.

`refine_tomo` passes them on and keeps the defaults: a tight 0.05 degrees and
0.1 A was never met on FePt (the poses keep moving about 0.1 degrees per
iteration), so every stage ran to `maxits` and drifted.

## Why the refinement is short

On FePt (34 projections per half, 2026-10-15) every iteration against a map
built from the current poses moved the geometry along two modes the data barely
constrain: a uniform stretch of the tilt range and a rotation of the whole
frame. FSC does not see it (both half maps stretch alike); the atom positions
against their final model do:

| Refinement | Tilt-range stretch | Mean frame tilt | Level-1 recall (precision) |
| --- | --- | --- | --- |
| 2 iterations at 2 A | 0.2% | 0.06 degrees | 0.42 (0.98) |
| 9 iterations, 2.2 to 1.6 A | 1.7% | 0.32 degrees | 0.38 (0.98) |
| 20 iterations, 4 to 1.66 A | 6.9% | 1.41 degrees | 0.22 (0.63) |

Their GENFIRE map gives 0.41. Two iterations at 2 A with the density mask is
the best map so far, so it is the default. Removing the two modes after each
stage (the goniometer fixes them, not the data) might allow longer refinement;
it is not implemented.

## import_tilt_series

`refine=cont` continues from the poses of a discrete workflow and checks each
particle's state, half-set and projection direction. `import_tilt_series`
makes such a project from an MRC stack and an angle file:

- `angconv=simple`: SIMPLE Euler angles, `e1= e2= e3=` and optional `x= y=`
  per line (an oritab);
- `angconv=tilt`: one tilt angle per line, single-axis tilt about image x;
- `angconv=genfire`: three angles per line, GENFIRE/RESIRE convention
  (`phi theta psi`, degrees).

It sets state 1, the half-set by index parity, and the projection direction
to the nearest node of a spiral grid, then adds the stack with `ctf=no`.

## The GENFIRE/RESIRE angle convention

Their rotation is `rotMAT = Rz(psi) Ry(theta) Rz(phi)` (standard active
rotations, acting on column vectors). Their central-plane point `(a,0,c)` sits
at volume frequency `rotMAT (a,0,c)`, with `c` along the first array dimension
of a projection and `a` along the second. In SIMPLE the plane point `(h,k)`
sits at `[h,k,0] R`, with `h` along image x. With image x the first array
dimension, `R(1,:) = rotMAT(:,3)` and `R(2,:) = rotMAT(:,1)`, and the volume
frames coincide, so their atom models drop into our maps after an origin and
pixel conversion. The single-axis case (`theta = psi = 0`) gives SIMPLE
`(90+phi, 90, 180)` and was verified on the FePt data (hand and frame check,
2026-10-14). The general case follows from their code and is checked by unit
tests (round trip and the single-axis anchor).

## Converting their MATLAB files

Their data come as `.mat` files (projections, volumes, angles, atom
positions), MATLAB v5 (usually compressed) or v7.3 (HDF5).
`scripts/aet_mat2simple.py` converts them without SIMPLE:

- a 3D array to an MRC stack or volume, MRC x along the first MATLAB
  dimension, optionally zero-padded to a larger box;
- an angle array to a text file in their convention (`angconv=genfire` or
  `tilt`);
- atom positions (pixels, 1-based) and types to a PDB in SIMPLE's map
  convention (A = 0-based voxel times the pixel size), with an optional origin
  shift for a padded box;
- a listing of the variables in a file.

It needs numpy only for v5 files and h5py for v7.3. On the FePt files it
reproduces the hand-made kit byte for byte: the padded projection stack, the
GENFIRE volume and their final model.

## Tests

- `unit_ori`, `orientation` sub-suite: the angle conversion (single-axis
  anchor, round trip of a grid of triplets through the rotation matrix).
- `unit_cart_align3D`, `pose statistics` sub-suite: the three parameters are
  honoured, and their defaults are the old constants.
- A ground-truth test (later): an atomic phantom projected as a tilt series
  with known angle errors; `refine_tomo` must reduce the orientation error
  against the truth. Scored on orientations, not FSC.

## Status

Written 2026-10-15, not compiled: the convergence parameters, the angle
conversion and its test, `import_tilt_series`, `refine_tomo` and the converter.
The ground-truth test is still to do.
