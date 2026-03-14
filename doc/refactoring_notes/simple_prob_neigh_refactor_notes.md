# Neighborhood probabilistic search refactor notes

These revised sources implement the structural fixes discussed in the review:

1. **Separate representative projections from neighborhood labels**
   - `subspace_inds` is treated as the representative full-space projection list.
   - `subspace_full2sub_map` is treated as the full-space `iproj -> isub` label map.
   - Neighborhood pooling uses `subspace_full2sub_map`; coarse scoring uses `subspace_inds`.

2. **Restore global assignment semantics**
   - Part jobs write sparse score tables.
   - The driver creates a metadata-only global object (`new_global`), reads all sparse part files (`read_tabs_to_glob`), then performs one global sparse `ref_assign`.
   - This matches the dense algorithm’s "score everywhere locally, assign once globally" contract.

3. **Allow state moves within a projection-direction neighborhood**
   - `build_sparse_from_mask` expands each allowed projection direction into candidates across **all active states**.
   - `proj_active_state_count(iproj)` is used to size the CSR rows correctly.

4. **Add existence-aware fallbacks**
   - Coarse subspace scoring skips inactive references.
   - If the pooled neighborhood becomes empty after existence filtering, the code falls back to:
     1. previous-projection neighborhood if valid,
     2. otherwise all projections that are active in at least one state.

5. **Additional fixes included in the proposed sources**
   - Bound-check `iproj` before indexing `ref_index_map` in the shift-first branch of `fill_tab_sparse`.
   - Kill `o_prev` inside the particle loop instead of once after the OpenMP region.
   - Make neighborhood score files use a distinct `_neigh_` suffix to avoid collisions with dense `prob_tab` runs.
   - Make `read_assignment` linear-time via a temporary `pind -> local slot` lookup table.
   - Make `build_ref_adjacency` safe to call after a prior adjacency has already been allocated.

## Integration note
The refactor assumes the builder provides:
- `subspace_inds(nspace_sub)` as the representative projection list.
- `subspace_full2sub_map(nspace)` as the full-space label map.

That builder-side field is already referenced by the revised neighborhood code.
