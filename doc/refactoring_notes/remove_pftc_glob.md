# Refactoring Document: Global `pftc_glob` Removal

## Overview
This refactoring eliminates the global `pftc_glob` variable from the SIMPLE codebase and replaces it with instance-based access through the `builder` object. This change improves object lifecycle management, particularly for distributed execution where `pftc` may be reinitialized during 3D refinement.

## Motivation
- **Object Lifecycle Management**: `pftc` (polarft_calc) is reinitialized during 3D refinement workflows when switching from initial setup to actual refinement phases
- **Dependency Injection**: Centralizes pftc instance management through the builder pattern

## Affected Modules

### Core Framework Changes

#### 1. Strategy Interface Updates
**Files Modified:**
- `simple_strategy2D.f90`
- `simple_strategy2D_srch.f90`
- `simple_strategy3D_srch.f90`

**Changes:**
- Added `import :: builder` to abstract interfaces
- Updated abstract `generic_new` interface to include `build` parameter:
  ```fortran
  subroutine generic_new( self, params, spec, build )
      import :: builder
      class(builder), intent(in) :: build
  end subroutine generic_new
  ```
- Added `b_ptr` (builder pointer) to both `strategy2D_srch` and `strategy3D_srch` types
- Constructor now stores builder reference for accessing `b_ptr%pftc`

**Impact:**
- All concrete strategy implementations must pass `build` to parent constructor
- Enables access to builder-managed pftc instance

#### 2. Strategy2D Implementation Updates
**Files Modified:**
- `simple_strategy2D_eval.f90`
- `simple_strategy2D_greedy.f90`
- `simple_strategy2D_greedy_smpl.f90`
- `simple_strategy2D_inpl.f90`
- `simple_strategy2D_inpl_smpl.f90`
- `simple_strategy2D_prob.f90`
- `simple_strategy2D_snhc.f90`
- `simple_strategy2D_snhc_smpl.f90`
- `simple_strategy2D_tseries.f90`

**Changes:**
- Updated all `new_*` constructors to accept `build` parameter
- Updated abstract interface calls from `self%s%new(params, spec)` to `self%s%new(params, spec, build)`
- Replaced all `pftc_glob%gen_objfun_vals(...)` calls with `self%s%b_ptr%pftc%gen_objfun_vals(...)`
- Replaced all `pftc_glob%gen_corr_for_rot_8(...)` calls with `self%s%b_ptr%pftc%gen_corr_for_rot_8(...)`
- Replaced all `pftc_glob%get_rot(...)` calls with `self%s%b_ptr%pftc%get_rot(...)`
- Replaced all `pftc_glob%get_roind(...)` calls with `self%s%b_ptr%pftc%get_roind(...)`

**Example Migration:**
```fortran
! Before
call pftc_glob%gen_objfun_vals(iref, self%s%iptcl, [0.,0.], corrs)

! After
call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, [0.,0.], corrs)
```

#### 3. Strategy3D Implementation Updates
**Files Modified:**
- `simple_strategy3D_greedy.f90`
- `simple_strategy3D_greedy_smpl.f90`
- `simple_strategy3D_greedy_sub.f90`
- `simple_strategy3D_shc.f90`
- `simple_strategy3D_shc_smpl.f90`
- `simple_strategy3D_snhc_smpl.f90`
- `simple_strategy3D_utils.f90`

**Changes:**
- Removed `use simple_polarft_calc, only: pftc_glob` imports
- Updated constructors to accept and use `build` parameter
- Replaced all `pftc_glob` references with `self%s%b_ptr%pftc` or `s%b_ptr%pftc`
- Added builder pointer to `strategy3D_srch` type for accessing pftc

**Key Updates in `simple_strategy3D_srch.f90`:**
```fortran
! New member variable
class(builder), pointer :: b_ptr => null()

! Constructor now captures builder reference
subroutine new( self, params, spec, build )
    self%b_ptr => build
    self%nrots = self%b_ptr%pftc%get_nrots()
end subroutine new
```

### High-Level Strategy Matchers

#### 4. `simple_strategy2D_matcher.f90`
**Changes:**
- Removed module-level `type(polarft_calc) :: pftc` declaration
- Updated public subroutine signatures:
  - `prep_strategy2D_glob(p_ptr, spproj, nrots, neigh_frac)` - added `nrots` parameter
  - `prep_strategy2D_batch(p_ptr, spproj, which_iter, ...)` - removed `pftc` parameter
  - `build_batch_particles2D(nptcls_here, pinds_here, ...)` - removed `pftc` parameter
  - `preppftc4align2D(batchsz_max, which_iter, l_stream)` - removed `pftc` parameter
  - `prep_polar_pftc4align2D(batchsz_max, which_iter, l_stream)` - removed `pftc` parameter

**Method Implementations:**
- All `pftc%method_name()` calls replaced with `b_ptr%pftc%method_name()`
- Strategy constructor calls updated: `call strategy2Dsrch(iptcl_batch)%ptr%new(p_ptr, strategy2Dspec, b_ptr)`
- Sigma2 calculation now passes pftc: `call b_ptr%esig%calc_sigma2(b_ptr%pftc, iptcl, orientation, 'class')`

#### 5. `simple_strategy3D_matcher.f90`
**Changes:**
- Removed module-level `type(polarft_calc) :: pftc` declaration
- Updated `prepare_refs_sigmas_ptcls` call to remove `pftc` parameter
- Updated `build_batch_particles` call to remove `pftc` parameter
- All `pftc%polar_cavger_*()` calls replaced with `b_ptr%pftc%polar_cavger_*()`
- Sigma2 calculation now passes pftc: `call b_ptr%esig%calc_sigma2(b_ptr%pftc, iptcl, orientation, 'proj')`

#### 6. `simple_strategy2D3D_common.f90`
**Changes:**
- `prepare_refs_sigmas_ptcls` signature: removed `pftc` parameter, now uses `build%pftc` internally
- `prepare_polar_references` signature: removed `pftc` parameter, now uses `build%pftc` internally
- `build_batch_particles` signature: removed `pftc` parameter, now uses `build%pftc` internally
- All function bodies updated to use `build%pftc` instead of passed `pftc`:
  ```fortran
  ! Before
  subroutine prepare_refs_sigmas_ptcls( params, build, pftc, cline, ... )
      call pftc%new(params, nrefs, [1,batchsz], params%kfromto)
  end subroutine
  
  ! After
  subroutine prepare_refs_sigmas_ptcls( params, build, cline, ... )
      call build%pftc%new(params, nrefs, [1,batchsz], params%kfromto)
  end subroutine
  ```
- Updated `read_groups` call: `call build%esig%read_groups(build%pftc, build%spproj_field)`
- Updated `new` call: `call build%esig%new(params, build%pftc, fname, params%box)`

#### 7. `simple_strategy2D_alloc.f90`
**Changes:**
- `prep_strategy2D_glob` now takes `nrots` parameter instead of accessing `pftc_glob%get_nrots()`
- `prep_strategy2D_batch` no longer takes `pftc` parameter
- Replaced `pftc_glob%get_nrots()` with parameter `nrots`

### Utility Modules

#### 8. `simple_strategy2D_utils.f90`
**Changes:**
- Added `use simple_builder, only: builder` import
- Functions `match_imgs2ref` and `match_imgs` now use local `type(builder) :: build` instead of `type(polarft_calc) :: pftc`
- All method calls on pftc now use `build%pftc%`
- Constructor calls updated: `call grad_shsrch_obj(ithr)%new(build, lims, ...)`

#### 9. `simple_corrmat.f90`
**Changes:**
- Updated `calc_inpl_invariant_cc_nomirr` to use `type(builder) :: build` instead of `type(polarft_calc) :: pftc`
- Constructor calls updated: `call grad_shsrch_obj(ithr)%new(build, lims, ...)`
- All method references use `build%pftc%`

### Shift Search Objects

#### 10. `simple_strategy2D_srch.f90` (pftc_shsrch_grad interface changes)
**Changes:**
- Gradient search object constructors now require `build` as first parameter:
  ```fortran
  ! Before
  call self%grad_shsrch_obj%new(lims, lims_init=lims_init, ...)
  
  ! After  
  call self%grad_shsrch_obj%new(self%b_ptr, lims, lims_init=lims_init, ...)
  ```
- All `pftc_glob` references replaced with `self%b_ptr%pftc`
- Constructor stores builder reference in strategy search object

#### 11. `simple_strategy3D_srch.f90` (pftc_shsrch_grad interface changes)
**Changes:**
- Same pattern as strategy2D_srch.f90
- All shift search constructor calls updated to pass `self%b_ptr` as first parameter
- All pftc method calls use `self%b_ptr%pftc%`

### One-Time Initialization

#### 12. `simple_strategy2D_matcher.f90` - `prep_strategy2D_glob` signature
**Before:**
```fortran
subroutine prep_strategy2D_glob( params, spproj, neigh_frac )
    ...
    s2D%snhc_smpl_ninpl = neighfrac2nsmpl(neigh_frac, pftc_glob%get_nrots())
```

**After:**
```fortran
subroutine prep_strategy2D_glob( params, spproj, nrots, neigh_frac )
    ...
    s2D%snhc_smpl_ninpl = neighfrac2nsmpl(neigh_frac, nrots)
```

**Call Sites Updated:**
- `simple_strategy2D_matcher.f90` line ~196: `call prep_strategy2D_glob(p_ptr, b_ptr%spproj, b_ptr%pftc%get_nrots(), neigh_frac)`

#### 13. `simple_strategy3D_alloc.f90`
**Changes:**
- `prep_strategy3D` allocation replaced `pftc_glob%get_nrots()` with `build%pftc%get_nrots()`
- Added `use simple_builder, only: builder` import

## Key Architectural Patterns

### 1. Builder-Based Injection Pattern
All code now follows this pattern:
```fortran
! Strategy objects store builder pointer
type strategy2D_srch
    class(builder), pointer :: b_ptr => null()
end type

! Constructor captures builder
subroutine new(self, params, spec, build)
    self%b_ptr => build
end subroutine

! Methods access pftc through builder
call self%b_ptr%pftc%gen_objfun_vals(...)
```

### 2. Gradient Search Objects
Shift search objects now require builder for initialization:
```fortran
subroutine new(self, build, lims, lims_init, ...)
    class(builder), intent(in) :: build
    ! Store builder reference internally or use it for initialization
end subroutine
```

### 3. Removal of Module-Level pftc
Pattern eliminated:
```fortran
! REMOVED - module-level declarations
type(polarft_calc) :: pftc  ! NO LONGER USED
class(builder), pointer :: b_ptr => null()  ! NOW USED
```

## Dependency Flow

```
builder (contains pftc)
    ↓
strategy_srch (stores b_ptr)
    ↓
concrete_strategy (uses self%s%b_ptr)
    ↓
pftc methods (accessed via self%s%b_ptr%pftc)
```

## Function Signature Changes Summary

| Old Signature | New Signature | Reason |
|---|---|---|
| `new(self, params, spec)` | `new(self, params, spec, build)` | Pass builder reference for pftc access |
| `prep_strategy2D_glob(params, spproj, neigh_frac)` | `prep_strategy2D_glob(params, spproj, nrots, neigh_frac)` | Eliminate pftc_glob dependency |
| `prep_strategy2D_batch(params, spproj, pftc, ...)` | `prep_strategy2D_batch(params, spproj, ...)` | Remove pftc parameter |
| `prepare_refs_sigmas_ptcls(params, build, pftc, ...)` | `prepare_refs_sigmas_ptcls(params, build, ...)` | Use build%pftc instead |
| `build_batch_particles(params, build, pftc, ...)` | `build_batch_particles(params, build, ...)` | Use build%pftc instead |

## Migration Path for New Code

When adding new strategy implementations:

1. **Accept builder in constructor:**
   ```fortran
   subroutine new_mystrat( self, params, spec, build )
       class(builder), intent(in) :: build
       call self%s%new(params, spec, build)
   end subroutine
   ```

2. **Access pftc through builder pointer:**
   ```fortran
   call self%s%b_ptr%pftc%gen_objfun_vals(...)
   call self%s%b_ptr%pftc%get_nrots()
   ```

3. **Initialize shift search objects with builder:**
   ```fortran
   call grad_shsrch_obj%new(self%b_ptr, lims, lims_init, ...)
   ```

## Files Modified: Complete List

**Strategy Core:**
- simple_strategy2D.f90
- simple_strategy2D_srch.f90
- simple_strategy3D_srch.f90

**Strategy2D Implementations (15 files):**
- simple_strategy2D_eval.f90
- simple_strategy2D_greedy.f90
- simple_strategy2D_greedy_smpl.f90
- simple_strategy2D_inpl.f90
- simple_strategy2D_inpl_smpl.f90
- simple_strategy2D_prob.f90
- simple_strategy2D_snhc.f90
- simple_strategy2D_snhc_smpl.f90
- simple_strategy2D_tseries.f90
- simple_strategy2D_matcher.f90
- simple_strategy2D_alloc.f90
- simple_strategy2D_utils.f90
- simple_eul_prob_tab2D.f90

**Strategy3D Implementations (9 files):**
- simple_strategy3D_greedy.f90
- simple_strategy3D_greedy_smpl.f90
- simple_strategy3D_greedy_sub.f90
- simple_strategy3D_shc.f90
- simple_strategy3D_shc_smpl.f90
- simple_strategy3D_snhc_smpl.f90
- simple_strategy3D_utils.f90
- simple_strategy3D_matcher.f90
- simple_strategy3D_alloc.f90

**Common Utilities:**
- simple_strategy2D3D_common.f90
- simple_corrmat.f90

**Total files modified: 36**

## Summary

This refactoring systematically eliminates the global `pftc_glob` variable by:

1. **Injecting builder references** into all strategy objects
2. **Updating constructors** to accept and propagate builder references
3. **Modifying function signatures** to eliminate pftc parameter passing
4. **Centralizing pftc access** through `build%pftc` pattern
5. **Enabling safe reinitialization** of pftc during distributed 3D refinement

The changes maintain backward compatibility for external interfaces while internally restructuring object lifetime management for improved robustness.
