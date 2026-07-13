# Distributed `flex_eigenvol` Refactoring Plan

Date: 2026-07-10

## Executive Summary

`flex_eigenvol` is currently a shared-memory projected-PPCA workflow. The
outer EM loop lives in `src/main/strategies/parallelization/simple_flex_eigenvol_strategy.f90`,
and the particle-dependent expectation/maximization kernels live in
`src/main/pca/simple_projected_latent_model.f90`.

The distributed implementation keeps particle Fourier planes bounded to the
worker's current batch. For the M-step, each worker immediately reduces those
planes into fixed-size sufficient statistics: one complex RHS grid per latent
component and one real cross-density grid per upper-triangular component pair.
The worker writes a binary M-step statistics stream and releases its particle buffers. The
master streams and sums those grids directly into its global accumulators,
then performs the existing coupled solve and basis finalization.

This avoids the former whole-part plane bundle, whose size grew with the number
of particles. Worker M-step memory is now independent of partition size, while
the shared-memory implementation remains the numerical reference path.

## Current Shape

The present `run_flex_eigenvol_linear` flow is:

1. Validate inputs, select active `ptcl3D` particles, and initialize global
   latent coordinates.
2. Initialize the supplied mean map and output basis reconstructors.
3. For each PPCA iteration:
   - `update_basis_from_latents` loops over particles, prepares the particle
     Fourier plane, subtracts the projected mean, inserts into global basis RHS
     and `rho_cross_exp`, solves the coupled per-voxel system, and finalizes
     basis volumes.
   - `infer_latents_from_basis` loops over particles, projects the basis, solves
     each particle posterior, and updates global latent tables.
   - The strategy logs residuals and orthonormalizes `z`.
4. Run a final refit/inference pass and write basis, coordinates, and trajectory
   volumes.

The refactor splits the particle-local 2D preparation/inference work from the
master-owned 3D stitching work.

## Target Architecture

### Strategy Layer

Introduce a strategy interface in `simple_flex_eigenvol_strategy.f90`, mirroring
the `cluster2D_strategy` orchestration pattern while keeping the public command
stable:

- `flex_eigenvol_strategy`, abstract base.
- `flex_eigenvol_shmem_strategy`, the existing single-process behavior.
- `flex_eigenvol_distr_master_strategy`, which schedules/streams partition
  work and performs shared-memory stitching.
- `flex_eigenvol_worker_strategy`, which prepares only 2D/per-particle products
  for one particle range.

Strategy selection follows the existing SIMPLE convention:

```text
nparts defined and part not defined -> distributed master
part defined                         -> distributed worker
otherwise                            -> shared-memory
```

The commander should call a single strategy factory or `run_flex_eigenvol`
entrypoint rather than knowing about execution mode details.

### Numerical Layer

Refactor `simple_projected_latent_model` around worker-local sufficient
statistics and master reduction routines:

- `write_mstep_stats_part_file`: stream bounded particle batches, build
  CTF/shift-corrected residual planes, and accumulate raw sufficient statistics
  into a fixed-size `.bin` file.
- `update_basis_from_mstep_stats_part_files`: validate and stream each part's
  RHS/cross-density grids into the master accumulators, deleting each file once
  consumed.
- `projected_latent_mstep_solve`: master-only routine that runs the existing
  coupled solve and basis finalization.
- `projected_latent_estep_part`: particle-local E-step routine that projects the
  current mean/basis into 2D planes and returns posterior rows, covariance rows,
  residual energies, and local mode-second sums.
- `projected_latent_estep_reduce`: master-only routine that updates global `z`,
  `z_postcov`, residual arrays, and `mode_vars` from per-part E-step products.

The existing `update_basis_from_latents` and `infer_latents_from_basis` can
remain as shared-memory wrappers over these lower-level routines. That keeps the
single-process algorithm as the reference implementation.

## Partition Products

The M-step plane-bundle design below is retained as historical context only.
The implemented M-step interchange is the fixed-size binary sufficient-
statistics stream described above. E-step posterior products continue to use
their own compact part files.

### Historical Per-Part Plane Bundle (Superseded)

Partition output should be represented as one typed bundle per part and PPCA
iteration, not as a family of loose files. The bundle is the interchange object
between partition workers and the master assembly phase.

Proposed filename:

```text
flex_eigenvol_partNNN_iterIII.bin
```

The object should be written and read through one bundle API; the `.bin`
extension reflects that it is a self-describing binary bundle rather than a
human-readable table.

The bundle should be a self-describing binary stream with:

- A fixed header containing magic, version, endianness marker, `part`, `nparts`,
  `iter`, `ncomp`, `fromp`, `top`, `box_crop`, `smpd_crop`, `kfromto`, record
  count, and content flags.
- A section table containing offsets and sizes for each logical section.
- A record index containing particle index, selected-row index, pose/orientation
  fields needed for insertion, and per-record plane metadata.
- An M-step section containing prepared residual 2D Fourier-plane data and
  latent first/second moments.
- An E-step section containing compact posterior/residual products after the
  E-step has run.

The Fortran-side owner should be a small I/O type, separate from the numerical
solver state, for example:

```fortran
type :: flexvol_part
    integer :: version, iter, part, nparts, ncomp
    integer :: fromp, top, nrecords, box_crop, kfromto(2)
    real    :: smpd_crop
    type(bundle_section_table) :: sections
contains
    procedure :: create_mstep
    procedure :: append_estep
    procedure :: read_header
    procedure :: stream_mstep_records
    procedure :: read_estep_products
    procedure :: validate_against_params
    procedure :: close
end type flexvol_part
```

Internally, the M-step section should be stored as structure-of-arrays blocks,
not as many tiny serialized derived types. That keeps I/O contiguous and lets
the master stream a bounded block of records into `projected_latent_mstep_insert_2d`.
The record index gives offsets into packed residual-plane storage.

This keeps the filesystem contract simple:

```text
one part, one iteration, one bundle
```

The M-step worker creates the bundle with the header, index, and M-step section.
After the master solves/finalizes the basis, the E-step worker for the same part
appends or rewrites the same bundle to add the E-step section. Writes should use
a temporary filename followed by an atomic rename so the master never reads a
partial bundle.

The linear solve should not read 2D records directly. The master assembly phase
streams each bundle's M-step section, reconstructs bounded `fplane_type`/pose
records, inserts them into the master's global 3D normal equations, and then
calls the existing coupled solver. In other words:

```text
part bundle -> master 3D assembly -> global normal equations -> linear solve
```

### Historical M-step 2D Records (Superseded)

Each M-step partition provides bounded batches of 2D records inside its part
bundle. For local workstation execution, the same data structure can also be
used in memory before it is written, but the durable external-worker interface
is the single per-part bundle.

Each M-step 2D record contains:

- Particle index and orientation parameters needed for 3D insertion.
- The prepared residual Fourier plane after mean subtraction.
- Transfer/CTF weighting needed by the insertion kernel.
- `z(row,:)`.
- `z_postcov(row,:,:) + z(row,:)*z(row,:)^T`.
- Frequency limits and plane metadata such as `nyq`, `frlims`, and `kfromto`.

The master consumes records immediately:

```text
for each partition:
    open part bundle
    for each bounded 2D-record block in the M-step section:
        master inserts records into global 3D accumulators
        discard block
master solves/finalizes basis once
```

No worker allocates `basis_recs(:)` for M-step output. No worker allocates or
writes `rho_cross_exp`. The only global 3D M-step state is the master's state.

### E-step Per-Particle Products

Each E-step partition returns compact per-particle products:

```text
particle index
z(row,:)
z_postcov(row,:,:)
resid_mean_energy(row)
resid_energy(row)
sum_q(z_q**2 + postcov_qq) contribution
```

These products live in E-step-specific per-part files. The master
reduces them by particle index or selected-particle row map into the global
arrays. Workers do not write project segments or mutate global orientation
state.

## Iteration Flow

Distributed master:

1. Parse params, build the project/general toolbox, select active particles, and
   initialize global latent arrays exactly once.
2. Initialize the supplied mean map, global basis reconstructors, and global
   M-step accumulators.
3. Build a `qsys_env` only when external worker scheduling is requested. For
   single-workstation memory reduction, prefer an in-process partition loop with
   bounded 2D batches.
4. For each PPCA iteration:
   - Reset global M-step accumulators.
   - For each partition, prepare bounded M-step 2D record batches and insert
     them immediately into the master's global 3D accumulators.
   - Run the coupled M-step solve and basis finalization on the master.
   - For each partition, run the E-step and collect compact per-particle
     posterior/residual products.
   - Reduce E-step products into global `z`, `z_postcov`, residual arrays, and
     `mode_vars`.
   - Log residual and latent statistics, then orthonormalize global `z`.
5. Run final refit/inference using the same 2D-record M-step and compact E-step
   products.
6. Write final basis volumes, coordinates table, trajectory volumes, and
   project/out metadata.
7. Clean up any temporary per-part bundles.

Distributed worker:

1. Parse the normal command line with `part`, `fromp`, `top`, and `nparts`.
2. Build only the toolbox state needed for its assigned particle range.
3. For M-step work, prepare residual 2D records in bounded batches and write
   them into the part bundle. Local partition providers can also pass the same
   block structure in memory before durable bundle I/O is needed.
4. For E-step work, project the current mean/basis into particle planes and
   write compact posterior/residual products.
5. Release particle batches, padded FFT heaps, and 2D record buffers promptly.

## Memory Policy

The purpose of partitioning on a single workstation is to limit live
particle-local 2D state. Each M-step worker now owns one fixed-size set of raw
3D sufficient statistics rather than retaining a particle-sized plane bundle.

Required memory rules:

- `nparts` is a particle partition count, not permission to keep all partitions
  resident at once.
- `ncunits` or the local scheduler must cap concurrent workers. For large boxes
  or many eigenvolumes, the expected workstation setting is often `ncunits=1`.
- The master owns the solved basis reconstructors, the reduced global
  `rho_cross_exp`, and basis-finalization scratch.
- Each M-step worker owns one raw RHS grid per component, one raw cross-density
  grid per component pair, its mean projection context, and only the particle
  buffers/latent rows for the current bounded batch.
- E-step partition code may read/project the current solved basis but must not
  allocate M-step accumulator tensors.
- Master-side reduction streams one M-step statistics file in bounded z-slabs and
  deletes it immediately after successful consumption.

The default workstation recommendation is:

```text
nparts = enough to bound particle-local work
ncunits = chosen by memory budget, commonly 1 for large boxes/neigs
```

## Command-Line and UI Changes

Add `nparts` to the `flex_eigenvol` UI under compute, alongside `nthr`. Reuse
the existing `nparts`, `part`, `fromp`, `top`, `numlen`, `ncunits`, and `qsys`
machinery.

Do not add ad hoc command-line keys for worker phase selection. Use separate
private worker programs/commanders for M-step bundle generation and E-step
bundle completion.

## File and Module Changes

Expected code touch points:

- `src/main/commanders/simple/simple_commanders_flex.f90`
  - Switch from direct linear runner to strategy factory.
  - Preserve current defaults and normal stop behavior.
- `src/main/strategies/parallelization/simple_flex_eigenvol_strategy.f90`
  - Add strategy types and distributed orchestration.
  - Keep output/log helpers here.
  - Add cleanup of temporary per-part bundle artifacts.
- `src/main/pca/simple_projected_latent_model.f90`
  - Accumulate bounded M-step particle batches into fixed-size worker-local
    sufficient statistics.
  - Write, validate, and slab-reduce binary statistics files before the existing
    master-only 3D solve/finalize.
  - Add compact E-step product write/read/reduce helpers.
  - Keep existing shared-memory wrappers behavior-compatible.
- `src/main/volume/simple_reconstructor_latent_ops.f90`
  - Add a raw-statistics insertion kernel numerically equivalent to direct
    reconstructor insertion.
- `production/tests/simple_test_projected_latent_mstep_stats.f90`
  - Compare direct and sufficient-statistics insertion and round-trip the
    stats stream.
- `src/main/ui/simple/simple_ui_denoise.f90`
  - Add `nparts` to `flex_eigenvol`.
- `production/simple_private_exec_driver.f90` and/or exec API wiring
  - Expose separate private M-step and E-step worker commands if qsys jobs run
    through `simple_private_exec`.

## Numerical Equivalence Policy

The shared-memory path remains the oracle. The partitioned path should match it
within floating-point reduction-order tolerance.

Important invariants:

- The selected active particle set must be identical between modes.
- Worker `fromp/top` ranges are only scheduling partitions; global output rows
  must be keyed by particle index or selected-particle row map.
- The M-step solve must see every prepared 2D record inserted exactly once into
  the master's global RHS accumulators and global `rho_cross_exp`.
- Ridge regularization is applied only during the master coupled solve, after
  all 2D records have been inserted.
- `mode_vars` are computed from global mode-second sums after all E-step
  products are reduced.
- Latent orthonormalization stays master-only so every partition begins the next
  iteration from the same global `z`.

## Implementation Phases

### Phase 1: Extract 2D M-step Record APIs

Refactor `simple_projected_latent_model` so the current shared-memory M-step
uses a bounded 2D record preparation routine followed by a master insertion
routine. Keep execution single-process.

Validation:

- Existing `flex_eigenvol` command produces the same files and comparable logs.
- The new 2D record path is exercised even in shared-memory mode.

### Phase 2: Add Local Partition Loop

Split the selected particle set into partition ranges inside the master process.
Process each range through bounded M-step 2D batches, insert immediately into
global accumulators, and discard the buffers before moving to the next range.

Validation:

- `nparts=1` matches direct shared-memory behavior.
- `nparts>1` matches within expected reduction-order tolerance.
- Peak memory is measured and lower than a naive multi-process setup.

### Phase 3: Split and Reduce E-step Products

Make the E-step return compact per-particle posterior/residual products by
partition and reduce them into global arrays on the master.

Validation:

- Global `z`, `z_postcov`, residual statistics, and `mode_vars` match the
  current shared-memory E-step within tolerance.
- Particle-index metadata catches missing, duplicated, or misordered rows.

### Phase 4: Add External Worker Scheduling

Wire qsys/private-exec scheduling through separate private M-step and E-step
worker commands. External workers still provide only 2D/per-particle data; the
master still owns all 3D stitching and solving.

Validation:

- qsys scripts are generated with `part`, `fromp`, `top`, and `nparts` for the
  appropriate private worker command.
- Master refuses to continue if any expected bundle is missing, incomplete, or
  fails metadata validation.
- Workstation runs can cap concurrency through `ncunits`.

### Phase 5: Complete Full EM Loop and Outputs

Wire all PPCA iterations, final refit/inference, final outputs, cleanup, and UI
documentation.

Validation:

- `flex_eigenvol nparts=1` through the partitioned path matches shared-memory.
- `flex_eigenvol nparts=2/4` completes and produces stable residual trends.
- Logs clearly separate 2D preparation, master insertion, master solve, and
  E-step reduction time.

## Risks and Mitigations

- 2D bundle I/O can be large for external workers. Mitigate first with the local
  in-process partition loop; for qsys, store records in bounded internal blocks
  and stream them into the master instead of staging all records.
- The master still owns `rho_cross_exp`, which scales with
  `ncomp*(ncomp+1)/2` times the expanded Fourier grid. Mitigate with the current
  `neigs <= 16` cap, explicit logging of accumulator sizes, and early allocation
  failure messages.
- E-step workers need access to solved basis projections. Mitigate by capping
  worker concurrency and, if needed, projecting/loading basis components in
  small groups without changing the M-step ownership rule.
- Row mapping bugs can silently corrupt `z`. Store particle indices in all
  E-step products and validate missing/duplicate rows on reduce.
- Worker/master parameter drift can corrupt insertion. Include metadata in the
  bundle header and fail on mismatched `ncomp`, box, sampling, `kfromto`, or
  iteration.
- Floating-point differences are expected because insertion/reduction order
  changes. Require tolerances, not bitwise identity, except for metadata and
  shape checks.
- Private exec wiring may not currently expose `flex_eigenvol` to scheduled
  workers. Confirm this early in Phase 4.

## Settled Design Decisions

- Use separate private worker commands/programs for M-step and E-step work,
  rather than a typed `flex_phase` parameter.
- Use normal image-batch-style bounded streaming for 2D record blocks. Large
  boxes are not a special file-format concern; SIMPLE already writes large-box
  particle records routinely. The block boundary is only a memory-streaming
  convenience for bundle assembly and master insertion.
- Keep the current `neigs <= 16` cap for the first partitioned version.

## Recommended First Patch Set

Start with Phase 1 only. The review target should be:

- A 2D M-step record type in `simple_projected_latent_model`.
- A particle-local routine that prepares bounded residual Fourier-plane records.
- A master-only routine that inserts those records into global basis RHS and
  global `rho_cross_exp`.
- Existing `update_basis_from_latents` rewritten as a wrapper over the new
  prepare/insert/solve sequence.
- No qsys scheduling yet.
- A before/after run on the same tiny dataset with residual and latent-stat log
  comparison.
