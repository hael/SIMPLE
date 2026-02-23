# Euclid Sigma Refactoring: Builder-Owned Policy

## Summary
This refactor consolidates Euclidean sigma ownership into the `builder` object and removes the previous implicit/global access pattern.

The key architectural decision is:

- `euclid_sigma2` must be instantiated and accessed through `builder%esig`.

## Why this was needed
- The previous model mixed ownership styles (local objects, global-like access, and indirect state), which made ML regularization paths fragile.
- Reconstruction and polar update flows consumed sigma with assumptions about allocation timing and bounds.
- Runtime failures (e.g. `sig2arr` shape/range mismatches) were easier to trigger when sigma initialization and sigma consumption were not tied to the same owner.

## What changed (from current diff)
### 1) Builder now owns sigma explicitly
- `builder` now includes:
  - `type(euclid_sigma2) :: esig`
- Added in:
  - `src/main/simple_builder.f90`

### 2) Global sigma exposure removed
- Removed `eucl_sigma2_glob` export and pointer lifecycle.
- `simple_image_ctf` / polar code paths now receive sigma arrays explicitly as arguments.
- Updated in:
  - `src/main/simple_euclid_sigma2.f90`
  - `src/main/image/simple_image.f90`
  - `src/main/image/simple_image_ctf.f90`
  - `src/main/pftc/simple_polarft_ops_state.f90`

### 3) Call graph switched to builder sigma
- Flows now read/write sigma through `build%esig` instead of local `eucl_sigma` in matcher/reconstruction paths.
- Representative updates:
  - `src/main/strategies/simple_strategy2D_matcher.f90`
  - `src/main/strategies/simple_strategy3D_matcher.f90`
  - `src/main/strategies/simple_strategy2D3D_common.f90`
  - `src/main/commanders/simple_commanders_rec.f90`
  - `src/main/commanders/simple_commanders_refine3D.f90`
  - `src/main/commanders/simple_commanders_refine3D_inmem.f90`

## New policy (mandatory)
### Ownership rule
- `euclid_sigma2` is **builder-owned state**.
- The canonical instance is `build%esig` (or `b_ptr%esig` in strategy modules).

### Instantiation rule
- Instantiate sigma only via builder-owned object:
  - `call build%esig%new(...)`
- Read/update sigma only via builder-owned object:
  - `call build%esig%read_groups(...)`
  - `call build%esig%calc_sigma2(...)`
  - `call build%esig%write_sigma2(...)`
  - `call build%esig%kill`

### API rule
- Do not reintroduce global sigma accessors/pointers.
- Pass sigma arrays explicitly into consumers that need raw arrays.

## Migration guidance for future edits
- If a function needs sigma values, prefer either:
  - builder access (`build%esig`) if the function already has builder context, or
  - explicit sigma array argument if a lower-level helper should remain data-oriented.
- Do not create local per-routine `type(euclid_sigma2)` objects in command/strategy/reconstruction paths.
- Keep allocation/read/kill lifecycle near execution boundaries (commander or strategy orchestration), not deep in low-level numerical kernels.

## Expected benefits
- Single source of truth for sigma state.
- Fewer hidden dependencies and less order-of-calls fragility.
- Better runtime diagnostics and easier debugging when sigma is missing/misaligned.
- Clearer ownership model for ML regularization and reconstruction pipelines.

## Regression checklist
- `ml_reg=yes` path allocates `build%esig` before any sigma-weighted reconstruction call.
- No new references to removed global sigma symbols.
- No new local `euclid_sigma2` lifecycles in orchestration paths.
- All sigma-consuming helper signatures remain explicit about input arrays where applicable.
