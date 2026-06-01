# 3D NU filter soft label synthesis

Date: 2026-06-01

## Summary

The 2D nonuniform filter now keeps the Potts field discrete for assignment and
statistics, but synthesizes the output from a tent-smoothed continuous label
coordinate field. That is the right separation of concerns:

- the discrete `filtmap` decides which low-pass member won the robust objective
  and remains the object analyzed by Potts/continuity statistics
- the synthesis field decides how sharply the final image switches between
  adjacent filtered products

The same idea should be moved into 3D. The current 3D implementation already has
the ingredients: an ordered `filtmap`, `candidate_coords`, a mask-packed voxel
list, cached filtered halfmaps, optional auxiliary replacement volumes, and
`tent_smooth_3d`.

## Current 3D Behavior

The hard-switch synthesis lives in:

- `src/main/nu_filt/simple_nu_filter_apply.f90`

`nu_filter_vols` currently:

1. Seeds the even output with the coarsest cached even volume.
2. For each finer label, reads that cached even volume and overwrites voxels
   whose `filtmap(i,j,k) == icut`.
3. Repeats the same procedure for the odd output.
4. If `nu_aux_replacement_label > 0`, overwrites voxels assigned to that label
   with the stashed auxiliary even/odd pair.

`nu_filter_vol` has the same hard-switch structure for a single input volume,
except it generates each filtered volume on demand and currently rejects
selected auxiliary replacement labels.

This is scientifically conservative, but it makes every label boundary a
pixel-sharp boundary in the final map. Even when the Potts statistics report
very few strong jumps, a visible discontinuity can remain because the final
volume jumps directly between two different low-pass products.

## Proposed Design

Keep the existing discrete field. Add a synthesis-only continuous coordinate
field:

```fortran
real, allocatable :: synth_coord(:,:,:)
real, allocatable :: synth_tmp(:,:,:)
```

or allocate equivalent local scratch in `simple_nu_filter_apply.f90`. It does
not need to survive cleanup and should not replace `filtmap`.

The synthesis coordinate is built from `candidate_coords`, not raw integer label
indices:

```fortran
synth_coord(i,j,k) = nu_candidate_coord_for_label(int(filtmap(i,j,k)))
```

Then smooth it with a small 3D tent kernel and clamp it to the retained
coordinate range. The user-facing/default parameter should be in Angstroms, not
pixels, so different samplings get the same physical transition width. A good
first default is:

```fortran
real, parameter :: NU_SYNTH_LABEL_SMOOTH_RADIUS_A = 10.0
```

Convert it to pixels only at the smoothing call:

```fortran
radius_px = max(1, nint(NU_SYNTH_LABEL_SMOOTH_RADIUS_A / smpd))
```

Ten Angstrom is the current 2D synthesis setting. If faint seams remain in 3D,
the next thing to try is increasing the physical radius rather than hard-coding
a larger pixel count.

## Masked Normalized Smoothing

Unlike the current 2D path, the 3D path has a real NU support mask. The synthesis
coordinate should therefore be smoothed with mask normalization so outside-mask
values do not dilute the particle support.

Implementation shape:

1. Fill `synth_coord` only inside `nu_lmask`.
2. Fill outside `nu_lmask` with zero for the smoothing pass.
3. Smooth `synth_coord` with `tent_smooth_3d`.
4. Smooth a mask-weight volume with the same pixel radius, where mask voxels are
   1 and outside voxels are 0.
5. For mask voxels, divide smoothed coordinate by smoothed mask weight.
6. Outside the mask, preserve existing semantics: use the coarsest filtered
   volume.

The existing `smooth_nu_objective` path already does this pattern through
`prepare_nu_smooth_norm`. The synthesis path can either reuse that normalization
machinery with the chosen synthesis radius, or factor it into a small generic
helper so objective smoothing and synthesis smoothing are both explicit.

## Soft Synthesis

After `synth_coord` is built, synthesize each output as a weighted blend of
neighboring bank members. For the current compressed coordinate system,
triangular weights are enough:

```fortran
weight = max(0., 1. - abs(synth_coord(i,j,k) - candidate_coords(icut)))
```

The weights should be applied only inside `nu_lmask`. Outside the mask, keep the
coarsest filtered volume exactly as today.

For `nu_filter_vols`:

1. Initialize even/odd outputs from the coarsest even/odd cached volumes.
2. Set masked voxels to zero before accumulating the soft blend.
3. Loop over `icut = 1,size(cutoff_finds)`.
4. For ordinary labels, read the cached even/odd filtered volumes.
5. For `nu_label_is_aux_replacement(icut)`, use `aux_even_bank(1)` and
   `aux_odd_bank(1)` as the synthesis source.
6. Add `weight * source` to each masked output voxel.

This keeps the aux behavior compatible with the current replacement semantics:
aux still backs the finest replacement label, but boundaries into that label are
softened instead of overwritten.

For `nu_filter_vol`:

- Keep the current rejection of selected aux replacement labels unless a
  single-map auxiliary source is explicitly added later.
- Otherwise use the same synthesis coordinate and triangular weights, generating
  each filtered candidate from `vol_in_ft` as today.

## Local Resolution Map

The existing `write_nu_local_resolution_map` reports the discrete `filtmap`
frequency. That remains useful for debugging the assignment field. After soft
synthesis, it may be worth adding a synthesis local-resolution map too, either
by replacing the current output or by adding a second optional writer.

The synthesis map should use the same interpolation weights as the output:

```fortran
freq_eff = sum_i weight_i / nu_label_lowpass_limit(i)
```

This handles auxiliary replacement naturally because `nu_label_lowpass_limit`
already returns `nu_aux_replacement_resolution` for the replacement label.

## Statistics

Do not reinterpret the existing Potts/continuity statistics as soft-field
statistics. They should continue to describe the discrete `filtmap`.

Add only a small synthesis note to the logs at first:

```text
>>> NU filter synthesis: blended bank members over a 10 A masked tent-smoothed label field
```

If more diagnostics are needed later, useful synthesis-field numbers would be:

- smoothing radius in Angstroms and voxels
- min/max/mean of `synth_coord` inside `nu_lmask`
- fraction of mask voxels with a fractional coordinate, i.e. transition voxels
- optional effective-frequency histogram from the soft local-resolution map

## Refactoring Steps

1. Add a physical synthesis radius parameter in `simple_nu_filter.f90`, near the
   Potts constants:

   ```fortran
   real, parameter :: NU_SYNTH_LABEL_SMOOTH_RADIUS_A = 10.0
   ```

2. Add module interfaces for helper routines implemented in
   `simple_nu_filter_apply.f90` or a new small submodule:

   ```fortran
   module subroutine build_nu_synthesis_coord(synth_coord, work)
   module real function nu_synthesis_weight(coord, ilabel)
   module subroutine synthesize_nu_filter_vols_soft(vol_even, vol_odd)
   ```

3. In `build_nu_synthesis_coord`, use `nu_mask_vox` and
   `nu_candidate_coord_for_label`. Use masked normalized smoothing with
   `tent_smooth_3d`.

4. Replace the hard overwrite loops in `nu_filter_vols` with soft accumulation.
   Preserve current outside-mask coarsest initialization.

5. Update `nu_filter_vol` similarly, but keep the current aux restriction.

6. Update `write_nu_local_resolution_map` only if we want it to reflect the
   synthesis field. If both maps are useful, add a second writer rather than
   changing the meaning of the existing discrete map silently.

7. Add one log line in `commander_refine3D` or the NU apply path describing the
   synthesis behavior.

## Risks And Checks

- Memory: two or three full-size real scratch arrays are needed if the
  normalized smoothing is local to the apply path. This should be acceptable for
  most refine3D runs, but the implementation should release the scratch before
  returning.
- Cache I/O: soft synthesis still reads each cached filtered volume once. It may
  be possible to combine even/odd accumulation in one loop over labels to reduce
  repeated traversal, but correctness comes first.
- Aux replacement: do not blend into aux unless `nu_label_is_aux_replacement`
  is true and the aux bank is allocated.
- High-resolution extension: the smoothed coordinate must use
  `candidate_coords` after extension thinning/compaction, not stale original
  labels.
- Matching low-pass policy: keep using the discrete finest selected label for
  policy decisions. Soft synthesis should not promote the alignment LP by itself.

## Validation Plan

1. Build the project and run existing metadata/lp-stage smoke tests.
2. Run a small `nu_filt3D`/`refine3D` case with the same inputs before and after.
3. Compare the discrete label histogram and continuity report; they should be
   unchanged except for any log additions.
4. Inspect filtered maps for boundary seams near high-label islands and shell
   extension frontiers.
5. Check FSC and visual map quality against the old hard-switch output. The goal
   is not to sharpen aggressively, but to remove synthesis discontinuities while
   preserving the NU assignment policy.
