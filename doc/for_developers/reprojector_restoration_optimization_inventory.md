# Reprojector and Restoration Optimization Inventory

**Date:** May 15, 2026
**Status:** Developer inventory / delegation note
**Audience:** Cyril and SIMPLE developers working on CPU/GPU performance kernels

## Purpose

This note inventories places where SIMPLE can reuse the optimization pattern introduced in the new reference reprojector:

- `src/main/pftc/simple_polarft_core.f90`
  - `vol_pad2ref_pfts_opt`
- `src/main/image/simple_projector_pft_batch.f90`
  - `fproject_polar_batch_opt`

The goal is to identify follow-up tasks that preserve current scientific behavior while reducing polymorphic calls, function calls, temporary allocation, and memory traffic inside the deepest projection and restoration loops.

## Current Optimized Pattern

The old reference-section path used this shape:

```fortran
do iproj = 1, nspace
    e_rotmat = eulspace%get_mat(iproj)
    do k = kfrom, kto
        do irot = 1, pftsz
            loc = rotated_polar_coordinate(...)
            refs_state(irot,kloc,iproj) = vol_pad%interp_fcomp_oversamp(loc)
        enddo
    enddo
enddo
```

That means every polar sample enters a type-bound interpolation call.

The optimized path instead:

1. Hoists polar coordinates into a packed `hkcoords` array.
2. Hoists all rotation matrices into `rotmats`.
3. Extracts `vol_pad%cmat_exp` and fixed interpolation metadata once.
4. Runs a collapsed OpenMP or OpenMP-offload kernel over `(iproj, kk, irot)`.
5. Computes KB weights with fixed small arrays inside the kernel.
6. Reads the expanded Fourier volume directly and writes `refs_state` directly.

The current implementation is in:

- `src/main/image/simple_projector_pft_batch.f90`
  - `fproject_polar_batch_opt`
  - inner kernel `fproject_polar_batch_opt_kernel`
- `src/main/pftc/simple_polarft_core.f90`
  - `vol_pad2ref_pfts_opt`
- `src/main/strategies/search/simple_matcher_refvol_utils.f90`
  - `read_mask_filter_reproject_refvols`

Important current limitation:

- `vol_pad2ref_pfts_opt` hard-stops for `mirr_proj=yes`.
- The kernel currently assumes `wdim <= 3` through fixed local KB weight storage.

## General Rules For These Refactors

Keep changes behavior-preserving first. Do not change search bands, weighting, CTF handling, or restoration policy while optimizing loops.

Prefer a two-step workflow for each task:

1. Add a new optimized routine alongside the old routine.
2. Validate numerical parity, then switch the caller or keep both behind a clear branch until enough benchmark coverage exists.

For CPU optimization:

- Hoist type-bound calls out of the innermost loop.
- Replace `matmul` on tiny vectors with scalar expressions.
- Replace per-sample array temporaries with fixed local scalars or small static arrays.
- Keep Fortran memory layout in mind: the first array index is contiguous.
- Avoid whole-array expressions in very hot loops if they create temporaries.

For GPU/offload readiness:

- Do not call polymorphic or type-bound procedures inside target kernels.
- Do not call normal type-bound KB methods inside target kernels. Use `apod_fast_device` or a scalar, non-polymorphic device helper.
- Prefer direct array arguments and explicit metadata over passing large derived-type objects into kernels.
- Treat scatter kernels separately from gather kernels; reconstruction gridding has race constraints that pure reprojection does not.

## Task 1: Finish The Optimized Mirrored Reference Projection

**Priority:** High

**Why:** This is the closest sibling of the new optimized reprojector and should be the cleanest win. It also removes the current unsupported branch for `mirr_proj=yes`.

**Current code:**

- `src/main/pftc/simple_polarft_core.f90`
  - `vol_pad2ref_pfts`
  - `vol_pad2ref_pfts_opt`
- `src/main/image/simple_projector_pft_batch.f90`
  - `fproject_polar_batch_mirr`
  - `fproject_polar_batch_opt`

**Current issue:**

`fproject_polar_batch_mirr` still extracts each unique section with:

```fortran
refs_state(irot,kloc,iproj) = vol_pad%interp_fcomp_oversamp(loc)
```

`vol_pad2ref_pfts_opt` currently throws for `mirr_proj=yes`.

**Suggested implementation:**

- Add `fproject_polar_batch_mirr_opt`.
- Reuse the optimized gather kernel for `iproj = 1:nspace/2`.
- Keep the existing mirror-fill logic:
  - get `m = eulspace%get_int(iproj, 'mirr')`
  - reverse/conjugate the in-plane samples as in `fproject_polar_batch_mirr`
  - respect the existing `psi` conditional.
- Route `vol_pad2ref_pfts_opt` to the mirrored optimized routine instead of throwing.

**Validation:**

- Compare old `vol_pad2ref_pfts` and new `vol_pad2ref_pfts_opt` for `mirr_proj=yes`.
- Validate even and odd reference banks.
- Check both CPU and `USE_OPENMP_OFFLOAD` builds if available.

**Acceptance criteria:**

- `mirr_proj=yes` no longer falls back or throws in the optimized materialization path.
- Max absolute and relative differences versus the old mirrored path are within expected floating-point tolerance.

## Task 2: Add A Batch Cartesian Reprojector For `reproject`

**Priority:** High to medium

**Why:** The public `reproject` path still uses one orientation at a time and calls the old Cartesian projector. It has the same gather nature as the optimized polar reprojector, so it is a good candidate.

**Current code:**

- `src/main/volume/simple_volinterp.f90`
  - `reproject`
- `src/main/image/simple_projector.f90`
  - `fproject`
  - `fproject_serial`
  - `interp_fcomp`

**Current issue:**

`reproject` calls:

```fortran
call vol_pad%fproject_serial(o2, imgs_pad(ithr))
```

`fproject_serial` calls `self%interp_fcomp(loc)` for each Cartesian Fourier-plane sample.

**Suggested implementation:**

- Add a batch Cartesian projection helper, for example near `simple_projector_pft_batch.f90` or a new `simple_projector_batch.f90`.
- Inputs should be the expanded volume matrix, orientation matrices, Fourier plane dimensions/limits, and output image/plane storage.
- Use the same scalarized interpolation style as `fproject_polar_batch_opt`.
- Consider producing Fourier planes first, then keeping the existing per-image `ifft` and `clip`.

**Validation:**

- Compare `simple_reproject` output stacks against the current path for a representative volume/orientation table.
- Include both `oritab` and `nspace` generated-orientation modes.
- Include `state` filtering and `neg=yes` at the command level, though those should be unaffected.

**Acceptance criteria:**

- Reprojection stack images match the old path within tolerance.
- Benchmark shows reduced time in projection before FFT/clip costs dominate.

## Task 3: Optimize Volume-PFT Extraction For Symmetry And Docking

**Priority:** Medium

**Why:** `volpft_corrcalc` repeatedly interpolates from reference and target volumes in tight loops. This is structurally close to the optimized reprojector and has no scatter-update races.

**Current code:**

- `src/main/volume/simple_volpft_corrcalc.f90`
  - `extract_ref`
  - `extract_target_1`
  - `extract_target_2`
- Callers include:
  - `src/main/volume/simple_volpft_symsrch.f90`
  - `src/main/volume/simple_dock_vols.f90`

**Current issue:**

The extraction routines repeatedly call:

```fortran
self%vol_ref%interp_fcomp(loc)
self%vol_target%interp_fcomp(loc)
```

inside `(ispace, k)` loops.

**Suggested implementation:**

- Add direct-array extraction kernels for reference and target V-PFT lines.
- Hoist:
  - line coordinates
  - rotation matrices
  - expanded matrix bounds
  - KB metadata
- For `extract_target_2`, handle the shift multiplier explicitly after interpolation.
- Keep this local to volume PFT code unless a reusable non-oversampled interpolation kernel naturally emerges.

**Validation:**

- Compare symmetry-axis search results before and after.
- Compare volume docking scores for fixed rotations and shifts.
- Check that `corr_mag`, `corr_1`, and shifted correlation paths remain numerically stable.

**Acceptance criteria:**

- Same best symmetry axis/docking result on reference test cases.
- Lower time in the grid-search scoring loop.

## Task 4: Split And Optimize 3D Reconstruction Gridding

**Priority:** Medium, larger design task

**Why:** Online 3D restoration spends meaningful time preparing Fourier planes and gridding them into even/odd reconstructors. The current gridding routine is already scalarized, but it still has type-bound KB calls and scatter updates.

**Current code:**

- `src/main/strategies/search/simple_matcher_3Drec.f90`
  - `prep_imgs4rec`
  - `update_rec`
  - `grid_ptcl`
- `src/main/image/simple_image_ctf.f90`
  - `gen_fplane4rec`
- `src/main/volume/simple_reconstructor.f90`
  - `insert_plane_oversamp`
- `src/main/volume/simple_reconstructor_eo.f90`
  - `grid_plane`

**Current issue:**

`insert_plane_oversamp` already follows much of the optimized style, but:

- KB weights still go through `self%kbwin%apod_fast`.
- The routine scatters into `self%cmat_exp` and `self%rho_exp`.
- The current race-avoidance relies on a stride over the `h` dimension.
- GPU offload would need a different update strategy or atomics.

**Suggested implementation:**

Treat this as two subtasks:

1. CPU cleanup:
   - replace the inner KB type-bound call with a non-polymorphic helper
   - hoist all scalar metadata out of loops
   - inspect whether the current `stride = self%wdim` race-avoidance leaves vectorization opportunities on the table
2. GPU/design experiment:
   - prototype per-thread or per-tile accumulation
   - compare atomics versus thread-private tiles plus reduction
   - keep the old CPU path until the scatter strategy is validated

**Validation:**

- Compare partial reconstruction files before and after for a fixed particle subset.
- Compare downstream `volassemble` halfmaps/FSC.
- Use existing match3D benchmark output:
  - `match3D partial reconstruction`
  - `match3D particle preparation`
- Use `volassemble` benchmark output for downstream effects.

**Acceptance criteria:**

- No reconstruction drift beyond expected floating-point differences.
- Lower `match3D partial reconstruction` time for online reconstruction.

## Task 5: Reduce Memory Traffic In 2D Class-Average Restoration

**Priority:** Medium

**Why:** This is the main analogous path for 2D restoration. It already has a fused interpolation/splat loop, but it materializes per-particle intermediate arrays and then performs a second class accumulation pass.

**Current code:**

- `src/main/strategies/search/simple_strategy2D_matcher.f90`
  - `restore_class_averages_for_batch`
- `src/main/class/simple_classaverager_restore.f90`
  - `cavger_update_sums`
  - `cavger_restore_cavgs`

**Current issue:**

`cavger_update_sums` fills:

```fortran
interp_cmats(:,:,i)
interp_rhos(:,:,i)
```

then later loops over classes and particles to add those intermediates into even/odd class sums.

**Suggested implementation:**

- Explore per-thread, per-class partial sums to avoid full per-particle `interp_cmats/interp_rhos` storage.
- Keep even/odd separation explicit.
- Preserve the single-read matcher contract: class-average restoration must reuse already-read batch images.
- Keep `cavger_restore_cavgs` unchanged initially; focus on `cavger_update_sums`.

**Validation:**

- Compare class-average partial sums and restored class averages before/after.
- Use existing `CLUSTER2D_BENCH_ITER*.txt` metrics:
  - `match2D class averaging`
  - `match2D cavg FFT/CTF/interpolation`
- Test both normal and fractional-update restoration paths.

**Acceptance criteria:**

- Same class assignment/restoration outputs within tolerance.
- Lower memory traffic and/or lower `match2D cavg FFT/CTF/interpolation` time.

## Task 6: Optimize `transform_ptcls`

**Priority:** Low to medium

**Why:** This helper has a classic nested Fourier rotation/interpolation loop with type-bound image accessors. It is not as central as online restoration, but it is easy to isolate.

**Current code:**

- `src/main/class/simple_classaverager_restore.f90`
  - `transform_ptcls`

**Current issue:**

The loop around Fourier rotation calls:

```fortran
call kbwin%apod_mat_2d_fast(...)
img(ithr)%get_cmat_at(...)
call timg(ithr)%set_cmat_at(...)
```

inside the deepest loops.

**Suggested implementation:**

- Replace accessor calls with direct `cmat` access where safe.
- Hoist FT address maps and KB metadata.
- Scalarize the rotation calculation.
- Consider a small direct 2D interpolation kernel shared with 2D restoration if it remains clean.

**Validation:**

- Compare transformed particle stacks and optional class average output.
- Exercise callers in:
  - `simple_commanders_cluster2D.f90`
  - `simple_cls_split_strategy.f90`

**Acceptance criteria:**

- Same transformed images within tolerance.
- No change in optional phase-flip and original-image output behavior.

## Task 7: Batch Polarization Of Particle And Reference PFTs

**Priority:** Medium, exploratory

**Why:** `polarize_oversamp` is already vector-friendly, but it is called image by image for particle and reference preparation. A batch version could improve cache behavior and make a future device path cleaner.

**Current code:**

- `src/main/image/simple_image_polar.f90`
  - `memoize4polarize_oversamp`
  - `polarize_oversamp`
- `src/main/pftc/simple_polarft_core.f90`
  - `polarize_ref_pft`
  - `polarize_ptcl_pft`
- `src/main/strategies/search/simple_matcher_ptcl_batch.f90`
  - `build_batch_particles3D`
  - `build_batch_particles2D`
- `src/main/strategies/search/simple_matcher_pftc_prep.f90`
  - `prep_pftc4align2D`

**Current issue:**

Particle preparation loops over images and calls `polarize_ptcl_pft` one particle at a time. Reference preparation does the same for class averages.

**Suggested implementation:**

- Add a batch polarizer that takes a stack or array of image `cmat` data and writes a contiguous PFT block.
- Reuse the existing memoized polar address/weight arrays.
- Avoid introducing hidden global state; keep memo ownership compatible with current image/PFTC setup.
- Preserve the range distinction:
  - `get_pdim_interp()` for restoration/averaging particle data
  - `get_pdim_srch()` for matching-only references

**Validation:**

- Compare `pfts_ptcls`, `pfts_refs_even`, and `pfts_refs_odd` arrays before/after.
- Check 2D and 3D matcher preparation.
- Confirm `memoize_ptcls` and CTF-matrix generation still consume the same data.

**Acceptance criteria:**

- Same PFT arrays within tolerance.
- Lower particle/reference preparation time in match2D/match3D benchmarks.

## Task 8: Remove Callback Overhead In Probabilistic Objective Sampling

**Priority:** Low for restoration, medium for probabilistic search

**Why:** This is not a restoration kernel, but it follows the same philosophy. The probabilistic sampling routines evaluate distances through a procedure callback, which is clean but expensive and unsuitable for GPU kernels.

**Current code:**

- `src/main/pftc/simple_polarft_corr.f90`
  - `gen_prob_euclid_val`
  - `gen_prob_power_euclid_val`
  - `euclid_dist_from_crvec`
- `src/main/simple_eul_prob_tab_utils.f90`
  - `sample_bounded_dist`
  - `sample_power_dist`

**Current issue:**

`sample_bounded_dist` and `sample_power_dist` call a procedure argument for each candidate rotation:

```fortran
pvec_sorted(i) = dist_fun(i)
```

**Suggested implementation:**

- Add array-based sampling variants that consume a precomputed distance vector.
- In `gen_prob_*`, compute the rotation-distance vector once from `crvec1`.
- Pass the vector to the array-based sampler.
- Keep the old callback sampler if it is still useful elsewhere.

**Validation:**

- Compare chosen distributions statistically over repeated seeded runs if deterministic identity is not expected.
- Compare deterministic pieces: sorted candidate distances, bounds, and top-ranked candidate.
- Use existing probability-table benchmark output for `prob_align` and `prob_tab`.

**Acceptance criteria:**

- Same sampled behavior under fixed RNG stream, or documented acceptable statistical equivalence if ordering changes.
- Lower probability-table fill time.

## Benchmark And Validation Checklist

For each task, collect before/after timing with the same:

- dataset
- command line
- thread count
- compiler flags
- CPU/GPU environment
- `kfromto`
- `nspace`
- `nstates`
- particle subset

Useful existing timing outputs:

- `CLUSTER2D_BENCH_ITER*.txt`
  - `match2D reference preparation`
  - `match2D particle preparation`
  - `match2D alignment search`
  - `match2D class averaging`
  - `match2D cavg FFT/CTF/interpolation`
- match3D benchmark file
  - `match3D particle preparation`
  - `match3D reference preparation`
  - `match3D orientation search`
  - `match3D partial reconstruction`
- `volassemble` benchmark file
  - `volassemble restore_eos_and_write_fsc`
  - `volassemble restore_merged_volume`
  - `volassemble trail_restored_halves`
  - `volassemble total time`
- probability-table benchmark files from `prob_align` / `prob_tab`
  - reference preparation
  - particle polar prep
  - table fill
  - assignment

Numerical checks should include:

- max absolute difference
- max relative difference where denominators are not tiny
- aggregate norms over full PFT/volume/class-average arrays
- final workflow-level outputs where possible:
  - class averages
  - partial reconstructions
  - assembled halfmaps
  - FSC/resolution summaries
  - assignment maps for probabilistic paths

## Suggested Delegation Order

1. Task 1: optimized mirrored reference projection.
2. Task 2: batch Cartesian `reproject`.
3. Task 3: optimized volume-PFT extraction for symmetry/docking.
4. Task 4: CPU cleanup of 3D reconstruction gridding, then GPU scatter experiment.
5. Task 5: 2D class-average restoration accumulation redesign.
6. Task 7: batch polarizer if preparation benchmarks remain significant.
7. Task 6 and Task 8 as opportunistic cleanup or GPU-enablement work.

This order starts with pure gather kernels, where Cyril's current optimized reprojector applies almost directly, before moving into scatter-heavy reconstruction/restoration kernels where correctness and race control dominate the design.
