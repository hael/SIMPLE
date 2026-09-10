# C-alpha finder

## Contract

The first model-building slice adds `simple_exec prg=detect_calpha`. It locates
oriented C-alpha seed candidates in a cryo-EM map using a built-in expected
backbone density target.

Inputs:

- `vol1`: work map to search.
- `smpd`: sampling distance in Angstroms per voxel.
- `angstep`: coarse SO(3) angular sampling step.
- `npeaks`: maximum number of non-overlapping candidates.
- `thres`: minimum weighted normalized-correlation score.

Outputs:

- `pdbout`: candidate C-alpha coordinates in PDB format, with the score in
  the B-factor field.
- A same-stem CSV file containing coordinates, scores, and orientation
  matrices. Each matrix maps the canonical target frame into the work
  map frame.
- `outvol`: maximum score over orientations at every work-map voxel.

The initial coordinate contract follows the existing SIMPLE atom/map
convention: voxel `(1,1,1)` corresponds to `(0,0,0)` Angstrom. MRC start and
origin metadata are not interpreted by this slice. The FFT-searched work map
must have even dimensions.

## Numerical method

1. Construct an ideal local `N-CA-C` geometry using 1.458 Angstrom `N-CA` and
   1.525 Angstrom `CA-C` bonds and a 111.2 degree `N-CA-C` angle.
2. Represent the three backbone atoms with Gaussian density. The C-alpha lobe
   is slightly stronger to anchor the reported coordinate at the target center.
   Gaussian sigma is `max(0.85 Angstrom, 0.75 * smpd)` to remain sampled on
   coarser voxel grids.
3. Sample this expected density inside a 4 Angstrom sphere and use a cosine
   taper as the target reliability weight.
4. Generate a coarse SO(3) grid from Fibonacci-sphere axes and uniformly
   sampled roll angles.
5. For every orientation, rotate the analytic target and evaluate weighted
   normalized cross-correlation at every translation. Three 3-D correlations
   provide the weighted data sum, squared-data sum, and target-data sum.
6. Retain the best orientation per voxel and select score maxima with a
   2 Angstrom exclusion radius.

The normalized score fits a local scale and offset implicitly and is the first
SIMPLE analogue of Buccaneer's cryo-EM correlation mode. Chain growth,
fragment joining, sequence assignment, side-chain classification, and model
refinement are deliberately outside this slice.

## Review findings

- `image%ccf` supplies the required FFT correlation primitive and preserves
  SIMPLE's phase-origin convention. A new `image%ccf_into` variant writes into
  an already allocated image so the search reuses FFT plans and buffers rather
  than rebuilding them for every orientation.
- An isolated spherical C-alpha Gaussian would make orientation search
  redundant and provide poor discrimination. Including the expected adjacent
  N and C density preserves a useful oriented Buccaneer-style target without
  requiring `vol2` or a reference model.
- The existing atom-centered validation code assumes the same zero-origin
  coordinate convention adopted here.
- The algorithm belongs in a new model-building domain. The commander owns
  I/O and lifecycle; the domain object owns analytic target construction and
  numerical search.

## Validation criteria

- Static Fortran checks report no overlong lines or invalid multiline
  diagnostic macros.
- The UI registry contains `detect_calpha` and the execution API routes it to
  one commander.
- `simple_test_exec test=detect_calpha` builds repeated analytic backbone
  density in memory at an orientation in the search grid, searches it, and
  verifies that a reported candidate is within 1.5 Angstrom of a known
  C-alpha position.
- User validation remains: compile the affected targets, run the synthetic
  test, and evaluate precision/recall on representative SIMPLE maps.

## Built-in molecule benchmark

`simple_test_exec test=detect_calpha_molecules` is an opt-in accuracy benchmark
using SIMPLE's hard-coded 6VXX spike-protein and 1JYX beta-galactosidase
coordinates. The fixtures contain 2,916 and 4,044 protein C-alpha atoms,
respectively. For each structure the benchmark:

1. creates a tightly padded, even-dimension density map at the requested
   sampling distance;
2. runs the same analytic-target search used by `simple_exec`;
3. performs greedy one-to-one matching of score-ordered predictions to the
   ground truth within 2 Angstrom; and
4. reports truth, prediction, unique-match, missed, false-positive, recall,
   and precision counts while retaining the input map, centered model,
   candidate PDB/CSV, and score volume.

The default benchmark settings are `smpd=1.3`, `angstep=45`, and `thres=0.25`.
They can be overridden on the test command line. No accuracy threshold is a
test failure yet; the first observed run establishes the baseline. Candidate
selection sorts threshold-eligible scores once before non-maximum suppression,
avoiding the previous full-volume scan for every requested peak.

### 2026-09-10 baseline

The Release `simple_test_exec` benchmark was run with eight threads and the
default settings. Matches are greedy, one-to-one, and within 2 Angstrom.

| Structure | Truth | Top-N matches | Top-N recall/precision | Top-2N matches | Top-2N recall | Top-2N precision |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 6VXX | 2,916 | 1,753 | 0.601 | 2,755 | 0.945 | 0.472 |
| 1JYX | 4,044 | 2,582 | 0.638 | 3,482 | 0.861 | 0.431 |

Both structures reached the `2 * ntruth` output cap. At `thres=0.25`, 103,078
6VXX voxels and 126,413 1JYX voxels were score-eligible, so this threshold does
not determine the reported operating points. Follow-up calibration should use
the score-ranked outputs to choose a precision/recall tradeoff before setting a
production threshold.
