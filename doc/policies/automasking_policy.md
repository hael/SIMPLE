# Automasking Policy

## Scope
This document covers automasking policy in 3D workflows, centered on:
- `automsk` (string policy input)
- `l_filemsk` (runtime mask-file active flag)
- `l_envfsc` (runtime envelope-for-FSC flag)

## Definitions
- **Envelope mask file (`mskfile`)**: Explicit or generated 3D mask volume used for envelope masking.
- **Automasking (`automsk`)**: Policy to generate an envelope mask from even/odd volumes (`yes`, `tight`, `no`).
- **FSC envelope masking (`envfsc`)**: Policy controlling whether envelope masking is applied before FSC calculation.
- **Circular fallback mask**: Soft spherical mask based on `msk`/`mskdiam` when no envelope policy is active.

## Current Flag Semantics (As Implemented)

### `automsk`
- Declared as `character(len=5)` with values `yes|tight|no`, default `no`.
- Used in reconstructor FSC paths and in mask generation commands.
- `tight` selects a tighter Otsu behavior in the mask algorithm.

### `l_filemsk`
- Set to `.true.` when `mskfile` is defined and exists.
- Governs whether pre-existing/generated mask file is used in many downstream steps.
- Has effects beyond FSC (e.g., center-shift logic and reference filtering paths).

### `l_envfsc`
- Derived from `envfsc != 'no'`.
- Primarily gates whether envelope masking is applied before FSC.

## Policy Precedence (Observed)
For FSC masking in reconstruction:
1. If `l_filemsk && l_envfsc`: use envelope mask from file.
2. Else if `automsk != 'no'`: generate automask (write `automask3D.mrc`), then apply envelope masking only if `l_envfsc`.
3. Else: use soft circular mask.

This is the effective precedence in `reconstructor_eo` FSC generation.

## Lifecycle / State Transition Model
1. **Initial state**: no file mask (`l_filemsk=false`), optional automask policy set by CLI/UI.
2. **Generation event**: automask can be generated from even/odd maps.
3. **Persistence event**: generated mask is written as `automask3D.mrc` (`MSKVOL_FILE`).
4. **Operational state**: workflow promotes to file-mask mode (`mskfile` set, `l_filemsk=true`).
5. **Downstream reuse**: subsequent iterations/stages consume the file mask as primary envelope source.

## Workflow-Specific Behavior

### Reconstruct / FSC core
- Reconstructor reads `mskfile` into `envmask` when `l_filemsk` is true.
- FSC paths apply precedence listed above.
- Secondary FSC helper (`calc_fsc4sampl_dens_correct`) currently does **not** include an automask branch, only file-mask+`l_envfsc` else spherical fallback.

### Refine3D
- `l_automsk = (automsk != 'no')` for runtime decisions.
- On selected iterations, automask is generated from even/odd volumes.
- Generated mask is written and promoted to active `mskfile` (`l_filemsk=true`), then reused.

### Ab initio 3D
- Uses global `l_automsk` and stage gate (`AUTOMSK_STAGE`) to activate automasking in later stages.
- Reconstruction command path removes `automsk`/`mskfile`; then explicit automask execution may run and inject `mskfile` back into refine command line.
- Activation logic currently enables global automasking only when `automsk == 'yes'` (not `tight`).

### Strategy / Filtering / Resolution utilities
- `l_filemsk` affects centering eligibility, low-pass estimation masking, ICM mask constraints, and filter preparation.
- FSC utility commanders and 3D uniform filter also branch on `l_filemsk`.

## Known Deviations / Ambiguities
1. **`envfsc` default mismatch in comments**
   - In-code field value defaults to `no`, while nearby comments/documentation text suggests `{yes}` in at least one location.
2. **`tight` activation asymmetry**
   - `automsk!='no'` enables automask behavior.
3. **FSC helper inconsistency**
   - `calc_fsc4sampl_dens_correct` has no `automsk` branch unlike main FSC generation path.

## Proposed Normative Policy (Recommended)
1. **Single precedence contract** (all FSC-related paths)
   - `mskfile` (if active and `envfsc=yes`) > `automsk` generation > circular fallback.
2. **Activation contract for `automsk`**
   - Treat both `yes` and `tight` as “enabled”; `tight` is a mode variant, not a separate enable flag.
3. **Scope contract for `envfsc`**
   - Define `envfsc` as FSC-only gating policy; document that other filtering/regularization may still use `mskfile`.
4. **Persistence contract**
   - Generated automask should promote to `mskfile` for deterministic downstream behavior unless explicitly disabled.

## Verification Checklist (for future policy hardening)
- [ ] `automsk=no`, no `mskfile` -> circular fallback in FSC.
- [ ] `automsk=yes`, `envfsc=no` -> automask may be generated, not applied to FSC.
- [ ] `automsk=yes`, `envfsc=yes` -> automask generated and applied to FSC.
- [ ] `mskfile` + `envfsc=yes` -> file mask applied to FSC.
- [ ] `mskfile` + `envfsc=no` -> FSC not envelope-masked (fallback behavior explicit).
- [ ] `automsk=tight` in ab initio behaves consistently with non-ab initio workflows.

## Source Pointers (Implementation Anchors)
- `src/main/simple_parameters.f90`
- `src/main/volume/simple_reconstructor_eo.f90`
- `src/main/image/simple_image_msk.f90`
- `src/main/commanders/simple_commanders_refine3D.f90`
- `src/main/simple_abinitio_config.f90`
- `src/main/simple_abinitio_utils.f90`
- `src/main/commanders/simple_commanders_abinitio.f90`
- `src/main/strategies/simple_matcher_2Dprep.f90`
- `src/main/commanders/simple_commanders_resolest.f90`
- `src/defs/simple_defs_fname.f90`
