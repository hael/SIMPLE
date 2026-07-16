# Polar PFT Range Policy Refactor (2026-03-19)

## Purpose
This note captures the refactor decisions and code changes made to separate:
- averaging/restoration paths, which must use the interpolation limit (`interpklim`, typically Nyquist), and
- matching/search paths, which should remain bounded by the search band upper limit (`kfromto(2)`).

## Final Policy
1. Generation and storage for averaging/restoration are full interpolation range:
`k = kfromto(1):interpklim`.
2. Matching objective functions and alignment scoring remain search range:
`k = kfromto(1):kfromto(2)`.
3. On-disk PFT/CTF2 files use header dims `[pftsz, kfromto(1), interpklim, ncls]`.
4. Read paths are strict: headers must match expected interpolation-aware dimensions.

## Why This Split Exists
- Averaging and Wiener/FRC restoration require shells beyond the alignment search band.
- Matching kernels are intentionally constrained for speed and statistical stability.
- Mixing these concerns caused range mismatches and fragile I/O behavior.

## API and Naming Cleanup
- Replaced ambiguous getters with explicit intent:
- `get_pdim_srch()` for matching/search range.
- `get_pdim_interp()` for interpolation/averaging range.
- `polarft_calc%new(...)` no longer takes optional `iklim`; `interpklim` is now derived from Nyquist (`fdim(box_crop)-1`).

## Core Invariants Introduced
- Particle/reference storage used for averaging must be allocated with `...:interpklim`.
- Any routine writing class-average partial sums must serialize to `interpklim`.
- Any routine computing correlation/objective/FRC for alignment continues to use search-band interfaces.

## Affected Paths

### Averaging/Restoration Paths (Interpolation Range)
- 3D particle build path memoizes oversampled polarization with interpolation dims.
- `src/main/strategies/search/simple_matcher_2Dprep.f90`
- `build_batch_particles`: `memoize4polarize_oversamp(build%pftc%get_pdim_interp())`

- 2D update path that feeds `polar_cavger_update_sums` now memoizes with interpolation dims.
- `src/main/strategies/search/simple_strategy2D_matcher.f90`
- `build_batch_particles2D`: `memoize4polarize_oversamp(b_ptr%pftc%get_pdim_interp())`

- Averaging test paths switched to interpolation dims.
- `production/tests/simple_test_polarops.f90`
- `src/main/commanders/test/simple_commanders_test_fft.f90` (polar-ops routine)

- Restore utilities and common-lines averaging interfaces use interpolation bounds.
- `src/main/pftc/simple_polarft_calc.f90`
- `src/main/pftc/simple_polarft_ops_restore.f90`

### Matching/Search Paths (Search Range)
- In-plane matching and scoring utilities remain search-bound.
- `src/main/strategies/search/simple_strategy2D_utils.f90`
- `src/utils/simple_corrmat.f90`
- Matching-only setup in `src/main/strategies/search/simple_strategy2D_matcher.f90` (`prep_pftc4align2D`) remains `get_pdim_srch()` by design.

## Partial Sum I/O Confirmation
`polar_cavger_readwrite_partial_sums('write')` writes arrays that are interpolation-range class sums (`pfts_even/odd`, `ctf2_even/odd`) and not raw matching buffers.

- Writer entry point:
`src/main/pftc/simple_polarft_ops_io.f90` (`polar_cavger_readwrite_partial_sums`)
- Payload writers:
`write_pft_array_local`, `write_ctf2_array_local`
- Header contract:
`[self%pftsz, self%kfromto(1), self%interpklim, self%ncls]`

Upstream, these class sums are accumulated from polarized particle data via `polar_cavger_update_sums` over `kfromto(1):interpklim`.

## Strict I/O Behavior
- Compatibility/padding behavior for short legacy files was removed.
- Read helpers now enforce exact header match and fail hard on mismatches.
- Transfer helpers perform direct full-range assignment for the stored interpolation extent.

## Practical Rule for Future Changes
When touching any code that calls `memoize4polarize` or `memoize4polarize_oversamp`:
1. If outputs feed `polar_cavger_update_sums` or restoration/class-sum state, use `get_pdim_interp()`.
2. If outputs are used only for matching/objective/FRC scoring, use `get_pdim_srch()`.
3. Do not widen matching kernels unless explicitly changing search policy.

## Files Most Relevant to This Refactor
- `src/main/pftc/simple_polarft_core.f90`
- `src/main/pftc/simple_polarft_access.f90`
- `src/main/pftc/simple_polarft_calc.f90`
- `src/main/pftc/simple_polarft_ops_state.f90`
- `src/main/pftc/simple_polarft_ops_restore.f90`
- `src/main/pftc/simple_polarft_ops_io.f90`
- `src/main/strategies/search/simple_matcher_2Dprep.f90`
- `src/main/strategies/search/simple_strategy2D_matcher.f90`
- `src/main/strategies/search/simple_strategy2D_utils.f90`

## Open Follow-Up
- Run full build and regression tests to validate all strategy/test combinations after the dimension-policy split.
