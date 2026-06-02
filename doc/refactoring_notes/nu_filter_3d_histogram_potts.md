# 3D NU filter histogram-constrained Potts regularization

Date: 2026-06-01

## Summary

The current 3D nonuniform filter uses a discrete ordered-label Potts cleanup
after the unary objective assigns `filtmap`. That cleanup reduces noisy label
fragmentation, but it can also change the label histogram freely. In practice
that means the Potts prior can either erase small high-resolution regions or
grow labels beyond the amount of support the robust unary objective actually
found.

The 2D NU filter handles this by adding a histogram term during ICM. For 3D,
the useful part is the histogram constraint itself: keep the raw unary label
fractions as a target, then let the Potts field smooth spatially while paying a
cost for deviating from that target histogram.

The goal is not to force every voxel to keep its original label. The goal is to
make the Potts prior redistribute labels spatially without changing the global
amount of each retained bank member too much.

## Current 3D Behavior

The 3D ordered-label cleanup lives in:

- `src/main/nu_filt/simple_nu_filter_potts.f90`

`refine_nu_candidate_map_ordered_labels` currently minimizes, per voxel:

```text
unary_objective(candidate) + beta * neighborhood_cost(candidate)
```

where:

- `dmats_mask(imask,icand)` is the unary objective
- `candidate_coords` defines ordered label distances, including after high-res
  shell extension thinning/compaction
- `nu_label_smooth_neighborhood_cost` penalizes large jumps to neighboring
  labels inside `nu_lmask`
- `estimate_nu_label_smooth_beta` scales the Potts weight from local unary
  gaps

The current color loop is OpenMP-parallel. That is fine for a local Potts prior,
but a histogram term has global state, so the implementation needs to be
careful.

## Proposed Energy

Add a global histogram penalty:

```text
E_hist = hist_scale * sum_l (count_l - target_l)^2
```

At each proposed voxel move from `cur` to `cand`, only two bins change, so the
incremental cost is cheap:

```fortran
cur_before  = real(counts(cur)  - target_counts(cur))
cand_before = real(counts(cand) - target_counts(cand))

delta_hist = hist_scale * ((cur_before  - 1.)**2 - cur_before**2 + &
                           (cand_before + 1.)**2 - cand_before**2)
```

The ICM candidate energy becomes:

```text
unary + beta * neighborhood + delta_hist
```

For 3D, the ordered-label Potts should continue to tolerate one-label neighbor
jumps at zero cost. Benchmarking showed that charging one-label jumps underfits
the data and slows resolution progression.

A good initial scale is the 2D policy:

```fortran
real, parameter :: NU_LABEL_HIST_BETA_FRAC = 8.0
hist_scale = NU_LABEL_HIST_BETA_FRAC * beta / real(max(1, n_nu_mask))
```

If 3D proves too resistant to spatial cleanup, reduce this to `4.0`; if labels
still collapse away from their unary-supported fractions, increase it.

## Target Histogram

The target histogram should be captured at the start of each ordered-label
cleanup call, immediately after the unconstrained unary competition between
filter-bank members:

```fortran
call count_nu_candidate_labels(candmap, n_candidates, target_counts)
counts = target_counts
```

This is important because there are two meaningful 3D cleanup sites:

1. After the initial base-bank unary assignment in `optimize_nu_cutoff_finds`.
   The target is the raw base-bank unary histogram.
2. After high-resolution shell extension in
   `refine_nu_extension_filtmap_ordered_labels`.
   The target is the post-extension, post-compaction label histogram just before
   final ordered-label cleanup.

Do not reuse stale target counts across extension compaction. Labels can be
thinned, remapped, or dropped, so the target must match the current
`candidate_coords`, `cutoff_finds`, and `dmats_mask` layout.

Auxiliary replacement labels do not need special histogram handling. In 3D, aux
backs the finest replacement label rather than appending a sidecar candidate, so
it can be counted as that label.

## Deterministic ICM Update

The implementation keeps the ordered-label ICM loop parallel by freezing global
histogram counts within each 8-color pass:

1. Compute `target_counts` from `candmap` at routine entry.
2. Set `counts = target_counts`.
3. During each color pass, evaluate candidate energies against read-only
   `counts`.
4. Accept same-color voxel changes in parallel.
5. Recount `counts` after each color pass before moving to the next color.
6. Include `nu_label_histogram_delta(...)` in each candidate energy.

The 8-color parity coloring separates voxels in the local neighborhood, so the
neighborhood term remains safe to update in parallel. Freezing the histogram
counts within a color pass is a softer approximation than live serial count
updates, but it avoids making the experimental gate much slower than the
working Potts implementation.

## Refactoring Steps

1. Add constants in `simple_nu_filter.f90` near the Potts constants:

   ```fortran
   real, parameter :: NU_LABEL_HIST_BETA_FRAC = 8.0
   ```

2. Add helpers in `simple_nu_filter_potts.f90`:

   ```fortran
   subroutine count_nu_candidate_labels(candmap, n_candidates, counts)
   real function nu_label_histogram_delta(icand, cur_icand, counts, target_counts, hist_scale)
   real function calc_nu_label_histogram_energy(counts, target_counts, hist_scale)
   ```

3. In `refine_nu_candidate_map_ordered_labels`, allocate:

   ```fortran
   integer, allocatable :: counts(:), target_counts(:)
   ```

   Then count labels at routine entry and compute `hist_scale`.

4. Add `nu_label_histogram_delta` to the candidate energy:

   ```fortran
   e = dmats_mask(imask,icand) + beta * neighborhood_cost + &
       &nu_label_histogram_delta(icand, cur_icand, counts, target_counts, hist_scale)
   ```

5. After each color pass, recount `counts` from the updated label map.

6. Update logging to report:

   - `NU_LABEL_HIST_BETA_FRAC`
   - target counts before smoothing
   - final counts after smoothing
   - histogram energy per voxel before/after
   - total label-count drift, e.g. `sum(abs(counts - target_counts))`

7. Leave `calc_nu_label_smooth_site_energy` either unchanged and add a separate
   histogram-energy line, or explicitly rename it so logs do not imply that the
   reported site energy includes the global histogram term.

## Pseudocode

```fortran
call count_nu_candidate_labels(candmap, n_candidates, target_counts)
counts = target_counts
hist_scale = NU_LABEL_HIST_BETA_FRAC * beta / real(max(1, n_nu_mask))

do iter = 1, NU_LABEL_SMOOTH_MAXITS
    nchanged = 0
    do color = 0, NU_LABEL_SMOOTH_NCOLORS - 1
        !$omp parallel do schedule(static) default(shared) &
        !$omp private(i,j,k,imask,icand,cur_icand,best_icand,n_full,nsz,e,best_e) &
        !$omp reduction(+:nchanged) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            if( nu_label_smooth_color(i,j,k) /= color ) cycle

            cur_icand = int(candmap(i,j,k))
            best_icand = cur_icand
            best_e = local_energy(cur_icand)

            do icand = 1, n_candidates
                if( icand == cur_icand ) cycle
                e = local_energy(icand) + &
                    &nu_label_histogram_delta(icand, cur_icand, counts, target_counts, hist_scale)
                if( nu_label_smooth_is_better(e, best_e) )then
                    best_e = e
                    best_icand = icand
                endif
            end do

            if( best_icand /= cur_icand )then
                candmap(i,j,k) = int(best_icand, kind=NU_LABEL_KIND)
                nchanged = nchanged + 1
            endif
        end do
        !$omp end parallel do
        call count_nu_candidate_labels(candmap, n_candidates, counts)
    end do
    if( nchanged == 0 ) exit
end do
```


## Risks

- Too strong a histogram term can prevent useful cleanup and preserve noisy
  unary speckles.
- Too weak a histogram term behaves like the current 3D Potts prior.
- Parallel live-count updates would be nondeterministic and should be avoided
  in the first implementation.
- High-resolution extension labels can be sparse; logging target/final counts
  is essential to catch accidental collapse or overgrowth.
- 2D-style support gating or weak one-label jump penalties are too conservative
  for 3D and can underfit by slowing resolution progression.

## Validation Plan

1. Run a small `nu_filt3D` or `refine3D` case with the current Potts and with
   histogram-constrained Potts.
2. Confirm raw target counts and final counts are close, but not necessarily
   identical.
3. Confirm neighbor continuity improves or stays comparable.
4. Inspect high-resolution extension labels to verify they do not collapse away
   unless the unary support was negligible.
