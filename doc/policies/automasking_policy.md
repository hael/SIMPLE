# Automasking Policy (Refactored Multi-State Architecture)

## Scope
This document covers automasking policy in 3D workflows (refine3D, abinitio3D, rec3D), centered on:
- `automsk` (string policy input: `yes`, `tight`, `no`) 
- `l_automsk` (runtime automasking enabled flag)
- `l_envfsc` (runtime envelope-for-FSC flag)

**Key change**: Masks are now **internal and state-specific**; `mskfile` is NO LONGER a CLI parameter.

## Definitions
- **Envelope mask**: 3D binary or soft-edged mask volume used for envelope masking/background suppression.
- **State-specific mask file**: Each reconstruction state gets its own mask file: `automask3D_state{N:02d}.mrc`.
- **Automasking (`automsk`)**: Policy to generate an envelope mask from even/odd volumes (`yes`, `tight`, `no`).
- **FSC envelope masking (`envfsc`)**: Policy controlling whether envelope masking is applied before FSC calculation.
- **Circular fallback mask**: Soft spherical mask based on `msk`/`mskdiam` when no automasking policy is active.

## Current Flag Semantics (Refactored)

### `automsk`
- Declared as `character(len=5)` with values `yes|tight|no`, default `no`.
- Activates mask generation in volassemble for all multi-state modes.
- `tight` selects tighter Otsu behavior in the mask algorithm.
- NOT passed as CLI parameter to refine3D; only used in volassemble internally.

### `l_automsk`
- Set to `.true.` when `automsk != 'no'`.
- Runtime flag in volassemble deciding whether to generate state-specific masks.
- Each state gets its own mask file if enabled.

### `l_envfsc`
- Derived from `envfsc != 'no'`.
- Governs whether envelope masking is applied before FSC calculation in reconstructor.

### `l_filemsk` (DEPRECATED for CLI use)
- No longer accepted as CLI input.
- Internally set by volassemble when state-specific mask file exists/is generated.
- Used in `simple_parameters.f90` to trigger hard error if user attempts to pass `mskfile` CLI.

## Policy Precedence (New)
For FSC masking in reconstruction:
1. If `automsk != 'no'` and state-specific mask file exists or can be generated:
   - Load or generate `automask3D_state{N:02d}.mrc`
   - Apply envelope masking if `l_envfsc=true`
2. Else:
   - Use soft circular mask (fallback)

This is implemented in `reconstructor_eo::load_state_mask_or_fallback()`.

## Lifecycle / State Transition Model (Refactored)
1. **Initialization**: User sets `automsk` policy (yes/tight/no); masks are internal, not on CLI.
2. **volassemble generation event**: For each state `s`:
   - If `automsk != 'no'`, generate mask from even/odd volumes
   - Write as `automask3D_state{N:02d}.mrc` (e.g., `automask3D_state01.mrc`)
   - Set internal flag `params%l_filemsk = true` only within volassemble
3. **reconstructor discovery event**: During FSC calculation:
   - Reconstructor_eo calls `load_state_mask_or_fallback(state)`
   - Looks for state-specific file by convention (no CLI coupling)
   - Loads from disk if exists, generates on-the-fly if `automsk != 'no'` and missing
   - Falls back to circular mask if `automsk == 'no'`
4. **Operational state**: Each state has independent mask, persisted on disk for reproducibility.
5. **Downstream reuse**: Reconstructor discovers masks autonomously; no CLI parameter required.

## Multi-State Support
- **All multi-state modes now support automasking**: `single`, `independent`, `docked`, `input_oris_start`, `input_oris_fixed`
- **Per-state mask persistence**: Each state generates/loads its own mask file
- **Deterministic discovery**: Reconstructor finds masks by naming convention `automask3D_state{N:02d}.mrc`
- **No CLI coupling**: Masks are never passed on command line

## Workflow-Specific Behavior

### volassemble (commander_volassemble)
- **Owns mask generation** for all Cartesian 3D workflows
- For each state `s` where `automsk != 'no'`:
  - Generates mask from even/odd volumes
  - Writes to `automask3D_state{N:02d}.mrc` (state-specific)
  - Sets internal `params%l_filemsk = true` (not exposed to CLI)
- Uses `AUTOMASK_FBODY` parameter from `simple_defs_fname.f90`

### Reconstruct3D / FSC core
- Reconstructor_eo **discovers masks by convention** via `load_state_mask_or_fallback(state)`:
  - Checks for state-specific mask file
  - If exists: loads from disk
  - If missing but `automsk != 'no'`: generates on-the-fly and saves
  - If `automsk == 'no'`: uses circular fallback
- FSC calculation applies envelope masking only if `l_envfsc=true`
- No CLI dependency; autonomous discovery

### Refine3D
- User specifies `automsk` policy (yes/tiny/no) on CLI
- Does NOT accept `mskfile` on CLI (rejected with error)
- Internally, volassemble generates state-specific masks
- Reconstructor discovers masks by convention
- Smooth automation: single policy flag controls all mask operations

### Ab Initio 3D
- Global `l_automsk` flag set from `automsk` parameter
- Supports automasking for **all multi-state modes** (was: single-state only)
- No artificial restriction; volassemble + reconstructor handle multi-state gracefully
- Mask generation deferred to volassemble, discovery deferred to reconstructor

### Strategy / Filtering / Resolution utilities
- `l_filemsk` internal flag no longer used in CLI-facing code paths
- Reference preparation always uses circular masking (volassemble's responsibility)
- FSC paths use autonomous mask discovery in reconstructor

## CLI Interface Changes

### Before (Deprecated)
```bash
simple refine3D prg=refine3D vol1=vol.mrc mskfile=automask3D.mrc automsk=no ...
```
- `mskfile` was explicit CLI parameter
- `automsk` and `mskfile` competed for control
- Single mask file `automask3D.mrc` overwritten by multi-state workflows

### After (Current)
```bash
simple refine3D prg=refine3D vol1=vol.mrc automsk=yes ...
```
- `automsk` is the only mask policy control
- State-specific masks generated automatically by volassemble
- Masks discovered automatically by reconstructor
- `mskfile` parameter now rejected with error

## Parameter Definitions & Enforcement

- **`AUTOMASK_FBODY`**: Defined in `simple_defs_fname.f90` as `'automask3D_state'`
  - Used to construct filenames: `automask3D_state01.mrc`, `automask3D_state02.mrc`, etc.
- **`MSKVOL_FILE`**: Legacy single-file mask (retained for standalone `mask` command only)
  - No longer used in refine3D/reconstruct/abinitio workflows

**Enforcement** (`simple_parameters.f90`): If `mskfile` detected on CLI, throws hard error with message:
```
"mskfile is no longer supported on command line; masks are internal and per-state"
```

**UI cleanup**: `mskfile` input removed from UI definitions

**Strategy cleanup**: No `mskfile` in output metadata

**Reference prep**: Simplified to circular mask only

## Verification Checklist (Refactored)
- [x] `automsk=no` -> circular fallback in FSC (always)
- [x] `automsk=yes`, `envfsc=no` -> automask may be generated, not applied to FSC
- [x] `automsk=yes`, `envfsc=yes` -> automask generated and applied to FSC
- [x] Multi-state workflows: each state gets independent mask file
- [x] `mskfile` on CLI -> rejected with error
- [x] `automsk=tight` in ab initio works for all multivol modes (was: single-state only)
- [x] Mask files persist on disk for reproducibility
- [x] Reconstructor discovers masks by convention (no CLI coupling)

## Migration Guide (From Old Policy)
1. **Remove `mskfile` from your command lines** – will now error
2. **Keep `automsk` parameter** – works as before but now in volassemble+reconstructor
3. **Multi-state workflows now support `automsk`** – previously restricted to single-state
4. **Mask files are state-specific** – check for `automask3D_state01.mrc`, `automask3D_state02.mrc`, etc.
5. **No output metadata changes needed** – masks no longer appear in project output

## Source Pointers (Implementation Anchors)
**New/Refactored:**
- `src/defs/simple_defs_fname.f90` – `AUTOMASK_FBODY` parameter definition
- `src/main/volume/simple_reconstructor_eo.f90` – `load_state_mask_or_fallback()` subroutine
- `src/main/commanders/simple/simple_commanders_rec_distr.f90` – volassemble state-specific mask generation
- `src/main/commanders/simple/simple_commanders_refine3D.f90` – removed MSKVOL_FILE assignments
- `src/main/commanders/simple/simple_commanders_abinitio.f90` – multi-state automasking enabled
- `src/main/simple_parameters.f90` – mskfile CLI rejection

**Related (unchanged logic):**
- `src/main/image/simple_image_msk.f90` – mask generation algorithms
- `src/main/strategies/simple_refine3D_strategy.f90` – removed mskfile from output
- `src/main/strategies/search/simple_matcher_refvol_utils.f90` – circular-mask-only reference prep
