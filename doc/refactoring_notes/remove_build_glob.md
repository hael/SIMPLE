# Remove build_glob

## Summary
This refactor removes the global builder pointer (`build_glob`) and replaces it with explicit `builder` passing and per-instance builder pointers where needed.

## Motivation
- Global state caused hidden coupling between commands and strategies.
- Shared-memory paths were especially fragile when local builders overwrote the global pointer.
- Explicit ownership makes call paths safer and easier to reason about.

## Key Changes
- Deleted `build_glob` from `simple_builder` and removed all assignments to it.
- Strategy3D now stores a `b_ptr` in `strategy3D_srch`, set by constructors.
- Commanders use local builder instances and pass them explicitly to helpers.
- `commander_base::execute_safe` no longer nullifies or restores `build_glob`.

## Behavior Notes
- Sigma generation and reconstruction now rely on the caller-owned builder instance.
- Any remaining implicit global access will fail at compile time, making regressions obvious.

## Files Touched (Representative)
- src/main/simple_builder.f90
- src/main/strategies/simple_strategy3D_srch.f90
- src/main/strategies/simple_strategy3D_utils.f90
- src/main/strategies/simple_strategy3D_matcher.f90
- src/main/commanders/simple_commander_base.f90
- src/main/commanders/simple_commanders_euclid.f90
- src/main/commanders/simple_commanders_refine3D.f90
