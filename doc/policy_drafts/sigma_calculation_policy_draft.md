# Sigma Calculation & ML Regularization Policy (Euclidean Refinement)

**Document Status**: Draft  
**Last Updated**: 2026-03-05  
**Scope**: 2D clustering / 3D refinement with `objfun='euclid'` and optionally with `ml_reg='yes'`  
**Audience**: Developers, refinement strategy architects, QA/test authors

---

## 1. Executive Summary

Sigma calculations in SIMPLE embody a **noise-weighted reconstruction philosophy** where particle/group-specific noise estimates guide Fourier-space weighting during 2D/3D image restoration. This policy defines:

- **When** sigmas are calculated (always with `objfun='euclid'`)
- **How** they accumulate across iterations (per-group averaging)
- **Where** they are applied (ML regularization, gated by `ml_reg='yes'`)
- **Who** owns them (builder object, single source of truth)

**Key Rule**: Sigmas are **calculated unconditionally** when Euclidean distance is the objective function, but **used for restoration only when both `ml_reg='yes'` AND `objfun='euclid'`**.

---

## 2. Definitions

### 2.1 Core Flags & Parameters

#### `objfun` (string: 'cc' | 'euclid' )
- **Purpose**: Selects the alignment/reconstruction metric
- **Enum conversion** (`simple_parameters.f90:268`):
  - `'cc'` → `OBJFUN_CC = 0`
  - `'euclid'` → `OBJFUN_EUCLID = 1`
- **Impact on sigma**: Only `objfun='euclid'` triggers sigma calculation

#### `ml_reg` (string: 'yes' | 'no')
- **Purpose**: Enable/disable ML-regularization (Wiener-filter-style weighting)
- **Default in all offline 2D clustering** (`simple_default_clines.f90:39`): **'yes'**
- **Default in refine3D_auto** (`simple_commanders_refine3D.f90:113`): **'yes'**
- **Derived boolean** (`simple_parameters.f90:1752–1756`):
  ```fortran
  self%l_ml_reg = trim(self%ml_reg).eq.'yes'
  if( self%l_ml_reg )then
      self%l_ml_reg = self%cc_objfun==OBJFUN_EUCLID  ! ← CRITICAL GATE
  endif
  ```
  **Result**: `l_ml_reg=.true.` **only if BOTH `ml_reg='yes'` AND `objfun='euclid'`**

#### `sigma_est` (string: 'group' | 'global')
- **'group'** (default): Group-based averaging (particles binned by stack/class ID, per even/odd)
- **'global'**: Single global average across all selected particles (single 'group' per even/odd)
- **Field**: `l_sigma_glob = (sigma_est == 'global')`

#### `l_filemsk`, `l_envfsc` (boolean)
- Not directly sigmas, but interact: automasking and FSC calculation use separate policies
- **Note**: Documented in sibling automasking policy

### 2.2 Sigma Data Structures

Type `euclid_sigma2` (`simple_euclid_sigma2.f90:18–44`):

| Member | Shape | Purpose |
|--------|-------|---------|
| `sigma2_noise(:,:)` | `[kfromto(1):kfromto(2), fromp:top]` | **Active** noise spectrum used in alignment & reconstruction (from group averages) |
| `sigma2_part(:,:)` | `[kfromto(1):kfromto(2), fromp:top]` | **Temporary** per-particle sigmas during iteration; aggregated into groups |
| `sigma2_groups(:,:,:)` | `[eo(2), ngroups, k]` | Group-averaged sigmas (read from `.star` files) |
| `pinds(:)`, `micinds(:)` | Integer arrays | Particle/micrograph index mappings |
| `fromp`, `top` | Scalars | Partition bounds `[fromp:top]` |
| `kfromto(2)` | Integer pair | Fourier shell range `[k_min, k_max]` |
| `binfname` | String | Binary filename for I/O |
| `exists` | Logical | Allocation flag |

**File Storage**:
- **Binary**: `SIGMA2_FBODY//int2str_pad(part,numlen)//'.dat'` = `'sigma2_part_X.dat'` (per partition)
  - Content: `real(kfromto(1):kfromto(2), fromp:top)` per-particle & group-aggregated sigmas
- **Starfile**: `SIGMA2_GROUP_FBODY//int2str(iter)//STAR_EXT` = `'sigma2_group_iter_N.star'`
  - Content: Group-averaged sigmas per even/odd, with EMDL metadata

---

## 3. Sigma Calculation Lifecycle

### 3.1 Phase 1: Initial Moise Power Estimation (Pre-Refinement)

**When triggered** (`simple_commanders_refine3D.f90:783`):
```fortran
if( .not.file_exists(sigma2_star_from_iter(params%startit)) )then
    call xcalc_pspec_distr%execute(cline_calc_pspec_distr)
endif
```
When no sigma `.star` file exists for starting iteration → distributed `calc_pspec` is invoked.

**What it computes** (`simple_commanders_euclid.f90:exec_calc_pspec`, line 103):
- Per partition, per particle: **noise power spectrum** from noise-masked FFT
- Formula (`simple_commanders_euclid.f90:153–157`):
  ```fortran
  call build%imgbatch(imatch)%norm_noise_mask_fft_powspec(build%lmsk, params%msk, pspec)
  If l_scale_update_frac:
      sigma2(:,iptcl) = pspec / (2.0 * params%update_frac)
  Else:
      sigma2(:,iptcl) = pspec / 2.0
  ```
  - Calculated for all frequencies of the input 2D particle (params%box x params%box) 
  - Particles are uniformly selected (no refinement information yet)
  - If sampling is active (`l_update_frac`), scale by reciprocal of sampling fraction (bias correction)

**Output**: Per-partition binary files `init_pspec_part_k.dat` (k = partition index)

**Group aggregation** (`simple_commanders_euclid.f90:exec_calc_pspec_assemble`, line 195):
- Read all partition files
- Group particles: by `stkind` if `sigma_est='group'`, else single global group
- Per group/eo: `group_pspecs(eo, igroup, k) = mean(pspecs(:, iptcls_in_group))`
- Subtract global power spectrum: `group_pspecs -= pspec_ave`
- Remove unphysical negative sigmas via neighbor propagation
- **Outputs**:
  - `SIGMA2_GROUP_FBODY_1.star` (iteration 1, group averages)
  - Updated `sigma2_part_k.dat` files (per partition, group-filled values)

### 3.2 Phase 2: Per-Particle Sigma During Refinement Iteration

**When triggered** (`simple_strategy3D_matcher.f90:250`):
```fortran
if ( p_ptr%cc_objfun == OBJFUN_EUCLID ) then
    call b_ptr%spproj_field%get_ori(iptcl, orientation)
    call orientation%set_shift(incr_shifts(:,iptcl_batch))
    call b_ptr%esig%calc_sigma2(b_ptr%pftc, iptcl, orientation, 'proj')
endif
```
After each particle search (alignment), sigma recalculated for that particle.

**Calculation** (`simple_euclid_sigma2.f90:calc_sigma2`, line 212):
```
For each particle after alignment:
  1. Get best reference (even or odd, based on particle's eo flag)
  2. Apply particle's shift & rotation
  3. Apply CTF correction
  4. Subtract particle FFT
  5. Compute per-shell residual: sigma2(k) = sum_over_FFT_box(|residual|^2) / (2 * pftsz),
     shell k in the [kfromto(1);kfromto(2)] range.
```

**Core formula** (`simple_polarft_corr.f90:gen_sigma_contrib`, line 547):
```fortran
pft_ref_tmp_8 = pft_ref_tmp_8 - pfts_ptcl(:,:,i)  ! residual
sigma_contrib = real(sum(|pft_ref_tmp_8|^2, dim=1) / (2.0 * pftsz))
```

**Semantics**: Sigma **measures mismatch** between refined reference and particle image
- **Lower sigma** = better correspondence = higher SNR
- **Higher sigma** = poor fit or noisy particle = lower SNR

### 3.3 Phase 3: Group Aggregation (End of Iteration)

**When triggered** (`simple_commanders_refine3D.f90:794`):
```fortran
if( trim(params%objfun).eq.'euclid' )then
    call cline_calc_group_sigmas%set('which_iter', iter)
    call xcalc_group_sigmas%execute(cline_calc_group_sigmas)
endif
```
After all particles have been searched in current iteration.

**Aggregation** (`simple_commanders_euclid.f90:exec_calc_group_sigmas`, line 357):
- Read all `sigma2_part_k.dat` files from current iteration
- Group particles: by `stkind` if `group`, else single global group
- Per group/eo: `group_pspecs(eo, igroup, k) = mean(sigma2(:, iptcls_in_group))`
- Remove negatives (same logic as Phase 1)
- **Outputs**:
  - `SIGMA2_GROUP_FBODY_iter_X.star` (group averages for iteration X)
  - Updated `sigma2_part_k.dat` files (group-filled for next iteration's usage)

**Rationale**: Group averaging **smooths noisy per-particle estimates** and **prepares σ for next iteration's reconstruction**.

### 3.4 Timeline Example

```
Iteration 0 (before refinement)
├── calc_pspec_distr (if no prior sigma.star)
│   ├── Per-partition: compute initial noise power spectra from noise-masked FFT
│   └── Assemble: group average → sigma2_group_iter_1.star
│
Iteration 1 (first refinement iteration)
├── Read sigma2_group_iter_1.star → b_ptr%esig%sigma2_noise (for alignment/reconstruction)
├── Per-particle loop:
│   └── calc_sigma2 on each particle after search → sigma2_part_k.dat
├── [End of iteration]
├── calc_group_sigmas
│   └── Average sigma2_part_k.dat → sigma2_group_iter_1.star (updated)
│
Iteration 2
├── Read sigma2_group_iter_1.star → b_ptr%esig%sigma2_noise
├── Per-particle loop:
│   └── calc_sigma2 on each particle after search → sigma2_part_k.dat
├── [End of iteration]
├── calc_group_sigmas
│   └── Average sigma2_part_k.dat → sigma2_group_iter_2.star
│
...continuation...
```

---

## 4. ML Regularization: Sigma Application in Reconstruction

### 4.1 Activation Gate

**Both conditions required**:
1. `ml_reg == 'yes'` (user opt-in)
2. `objfun == 'euclid'` (Euclidean distance metric)

**Result**: Boolean `l_ml_reg` in parameters (`simple_parameters.f90:1752–1756`)

**If either condition fails**:
- Sigmas **still calculated** (for diagnostics, logging, potential future use)
- But **not applied** to reconstruction (uniform Fourier weights used instead)

### 4.2 Weighting Formula (Wiener-Filter-Like)

Applied in Fourier-domain reconstruction (`simple_image_ctf.f90:gen_fplane4rec`, line 265–267):
```fortran
if (l_ml_reg) then
    c      = c      / sigma2_noise(shell)      ! ← Fourier component weighted
    tvalsq = tvalsq / sigma2_noise(shell)      ! ← CTF magnitude weighted
end if
```

**Interpretation**:
- **High sigma** (poor fit, large residual) → small weight → Fourier component contributes less
- **Low sigma** (good fit, small residual) → large weight → Fourier component contributes more
- **Effect**: Down-weights particles with high noise; up-weights particles with strong signal

### 4.3 Sigma Upsampling for Padded Reconstruction

Since reconstruction uses **padded images** but sigmas computed at **original size**, interpolation is required.

**Implementation** (`simple_math_ft.f90:upsample_sigma2`, line 271):
```fortran
subroutine upsample_sigma2( kstart, nyq, sigma2, nyq_out, sigma2_out)
    ! Linear interpolation between Fourier shells
    ! Boundary handling:
    !   k < kstart → sigma2(kstart)   [lower frequency extrapolation]
    !   k > nyq    → sigma2(nyq)      [higher frequency extrapolation]
```

**Called from** (`simple_image_ctf.f90:204`):
```fortran
call upsample_sigma2(kfromto(1), sigma_nyq, sig2arr, fplane%nyq, sigma2_noise)
```
Converts cropped-resolution sigmas to padded-resolution sigmas before Fourier weighting.

### 4.4 Class Averaging Integration (2D Clustering)

When `ml_reg='yes'` in 2D class averaging (`simple_classaverager_restore.f90:93–97`):
```fortran
if( p_ptr%l_ml_reg )then
    fname = SIGMA2_FBODY//int2str_pad(p_ptr%part,p_ptr%numlen)//'.dat'
    call b_ptr%esig%new(p_ptr, b_ptr%pftc, fname, p_ptr%box)
    call b_ptr%esig%read_part(  b_ptr%spproj_field)      ! per-particle sigmas
    call b_ptr%esig%read_groups(b_ptr%pftc, b_ptr%spproj_field)  ! group sigmas for alignment
endif
```

**Allocation** (`simple_classaverager_restore.f90:257–265`):
- Allocates thread-local `sigma2(1:nyq_crop, nthr_glob)` and upsampled variants
- Bounds: `sigma2_kfromto(1:2)` extracted from `b_ptr%esig%sigma2_noise`

**Usage during class averaging**:
- Padded image weighing use the same interpolation and formula as 3D (`simple_classaverager_restore.f90:324-338`)

---

## 5. Workflow-Specific Behaviors & Constraints

### 5.1 Sigma Estimation Methods

**Default: `sigma_est='group'`**
- Groups determined by **stack index** (`stkind` in particle orientation database)
- One sigma spectrum **per group per even/odd**
- Accounts for **group-specific quality variation** (e.g., different data sources, collection times)
- Recommended for: mixed datasets, multiple acquisition sessions

**Alternative: `sigma_est='global'`**
- Single group for all particles (both even & odd)
- One global sigma spectrum per even/odd
- May underestimate noise if data quality varies significantly by group
- **Used in**:
  - Nano3D multi-body pipeline (`single_commanders_nano3D.f90:260`)
  - Restricted data scenarios (single micrograph set, uniform quality)
  - Is the default strategy for all 2D clustering/3D abinitio

### 5.2 Update Fraction Scaling

When `l_update_frac=true` (partial particle sampling in early iterations):
```fortran
sigma2(:,iptcl) = pspec / (2.0 * params%update_frac)  ! ← scaled by reciprocal
```

**Reason**: Corrects bias from incomplete particle sampling
- Ensures group averages remain unbiased when only subset of particles included
- Only applied in initial Phase 1 (`calc_pspec`), not in per-iteration calculations

### 5.3 Special Cases: ML-reg Disabled or Modified

#### Case 1: ICM Filtering Stage (Ab Initio Pipeline)
- **Location** (`simple_abinitio_config.f90:28`):
  ```fortran
  integer, parameter :: ICM_STAGE = PROBREFINE_STAGE
  ! Comment: "we switch FROM ML regularization when refine=prob is switched on"
  ```
- **Reason**: ICM (Iterated Conditional Modes) regularization incompatible with ML-style Wiener weighting
- **Action**: At `ICM_STAGE`, transition to ICM; ML-reg disabled

#### Case 2: Nano3D Multi-Body Pipeline
- **Location** (`single_commanders_nano3D.f90:259`):
  ```fortran
  if( .not. cline%defined('ml_reg') ) call cline%set('ml_reg', 'no')
  ```
- **Reason**: Too few atoms for statistically valid noise estimates
- **Action**: Explicitly disables ML-reg; uses uniform weighting instead

#### Case 3: Objfun Mismatch
- If `objfun='cc'`: Sigmas not calculated; `l_ml_reg` forced false regardless of `ml_reg` setting

### 5.4 Continuation & Carryover

When `params%continue='yes'` (resuming from previous 3D refinement round):
```fortran
if( params%cc_objfun==OBJFUN_EUCLID )then
    call simple_list_files(prev_refine_path%to_char()//SIGMA2_FBODY//'*', list)
    do i=1,nfiles
        target_name = string(PATH_HERE)//basename(list(i))
        call simple_copy_file(list(i), target_name)
    end do
endif
```
All `sigma2_part_k.dat` files from previous run carried over to maintain continuity.

---

## 6. Builder Ownership Model

Per refactoring policy (`doc/refactoring_notes/euclid_build_new_policy.md`):

### 6.1 Mandatory Rules

1. **`type(euclid_sigma2)` instantiation**:
   - **Must** be instantiated via builder: `call build%esig%new(...)`
   - **Never** create local `euclid_sigma2` objects in strategy or command paths

2. **Access**:
   - All sigma access via `.`b_ptr%esig` in strategies (pointer to builder)
   - All sigma access via `build%esig` in commanders (owned object)
   - No global sigma variables or pointers

3. **Lifecycle boundaries**:
   - **Allocation**: Orchestration layer (commanders, top-level strategies)
   - **Update**: Per-iteration and per-particle during refinement
   - **Deallocation**: Orchestration layer post-iteration or on exit

4. **Explicit passing**:
   - Low-level helpers (e.g., `gen_fplane4rec`, `gen_sigma_contrib`) receive **sigma arrays as arguments**, not builder reference
   - Keeps numerical kernels data-oriented; decouples from builder dependency

### 6.2 Invariants to Maintain

- Before any sigma-weighted reconstruction call, `build%esig` must be **allocated and populated**
- `sigma2_noise` array bounds **must match** `[kfromto(1):kfromto(2), fromp:top]` for active partition
- After sigma write, partition-local `sigma2_part_k.dat` file **must be consistent** with in-memory `sigma2_noise`
- Shape mismatches (e.g., from cropped/padded dimension inconsistencies) **must be caught early** via `binfile%write` validation

---

## 7. Known Deviations & Future Work

### 7.1 FSC Helper (`calc_fsc4sampl_dens_correct`)

**Issue**: FSC calculation helper has **no automasking branch** while main FSC path uses automasking precedence rules. Sigma-based FSC weighting not yet unified.

**Status**: Documented but not yet resolved; listed for future policy harmonization.

**Recommended action**: Align FSC helper to use same masking/sigma precedence as main reconstruction.

---

## 8. Verification Checklist

- [ ] Sigmas calculated **only** when `objfun='euclid'` (verified in Phase 1 trigger and per-particle loop gate)
- [ ] Sigmas applied **only** when `l_ml_reg = (ml_reg=='yes' AND objfun=='euclid')` (verified in `gen_fplane4rec`)
- [ ] Initial power spectra generated via `calc_pspec` if no prior iteration's `.star` file exists
- [ ] Per-particle sigma updated during refinement after each particle search (per-iteration loop integration)
- [ ] Group averaging **triggered at end of every iteration** via `calc_group_sigmas` call
- [ ] Sigma upsampling applied before Fourier-domain weighting (checked in `gen_fplane4rec` logic)
- [ ] Update-fraction scaling applied only in Phase 1 initial pspec calc (verified in `exec_calc_pspec`)
- [ ] `build%esig` is builder-owned; no local `euclid_sigma2` instantiations in strategies
- [ ] Sigma filenames use constants: `SIGMA2_FBODY`, `SIGMA2_GROUP_FBODY` (audit codebase for hardcoded strings)
- [ ] Partition carryover on `continue='yes'` includes all `sigma2_part_k.dat` files (checked in continuation branch)
- [ ] ICM stage properly disables ML-reg transition (verified in ab initio config)
- [ ] Nano3D explicitly sets `ml_reg='no'` (verified in nano3D defaults)

---

## 9. Source Code Anchor Points

| Concept | File | Lines | Key Reference |
|---------|------|-------|----------------|
| Parameter gate derivation | `simple_parameters.f90` | 1752–1756 | `l_ml_reg` initialization |
| Sigma type definition | `simple_euclid_sigma2.f90` | 18–44 | Data structures, methods |
| Constructor | `simple_euclid_sigma2.f90` | 54–83 | `new()` subroutine |
| Initial pspec calc | `simple_commanders_euclid.f90` | 103–175 | `exec_calc_pspec` |
| Pspec assembly & grouping | `simple_commanders_euclid.f90` | 195–340 | `exec_calc_pspec_assemble` |
| Per-particle calc | `simple_euclid_sigma2.f90` | 212–232 | `calc_sigma2()` |
| Gen sigma contrib formula | `simple_polarft_corr.f90` | 547–583 | `gen_sigma_contrib()` |
| ML weighting applied | `simple_image_ctf.f90` | 265–267 | Wiener filter formula |
| Sigma upsampling | `simple_math_ft.f90` | 271–300 | `upsample_sigma2()` |
| Per-iteration integration | `simple_strategy3D_matcher.f90` | 250–255 | Refinement loop call |
| Iteration-end grouping | `simple_commanders_refine3D.f90` | 790–800 | Main loop trigger |
| Class averaging integration | `simple_classaverager_restore.f90` | 93–97 | 2D read phase |
| Class averaging ML weighing | `simple_classaverager_restore.f90` | 324–338 | Wiener filter |
| ICM transition | `simple_abinitio_config.f90` | 28 | Disable marker |
| Nano3D defaults | `single_commanders_nano3D.f90` | 259–260 | ML-reg explicit disable |

---

## 10. Policy Examples

### Example 1: Typical Refine3D with ML-reg

```
User command:
  simple_refine3D input.simple ...objfun=euclid ml_reg=yes...

Sequence:
  1. params%new() → l_ml_reg = true (both conditions met)
  2. No prior sigma.star → calc_pspec_distr invoked
     - Computes initial power spectra from noise-masked FFT
     - Groups by stkind (default)
     - Outputs: sigma2_group_iter_1.star
  3. Iteration 1:
     - Load sigma2_group_iter_1.star → b_ptr%esig%sigma2_noise
     - Per particle: calc_sigma2() → sigma2_part_k.dat updated
     - End iteration: calc_group_sigmas() → sigma2_group_iter_1.star re-written
     - Reconstruction: gen_fplane4rec() divides Fourier components by sigma2_noise (weighting applied)
  4. Iteration 2, 3, ... : repeat
```

### Example 2: Refine3D without ML-reg

```
User command:
  simple_refine3D input.simple ...objfun=euclid ml_reg=no...

Sequence:
  1. params%new() → l_ml_reg = false (ml_reg='no' fails gate)
  2. calc_pspec & calc_group_sigmas still run (per objfun='euclid')
  3. Per-particle calc_sigma2 still runs (per objfun='euclid')
  4. Reconstruction: gen_fplane4rec() skips sigma weighting → uniform Fourier weights
  Result: Sigmas calculated but not used; reconstruction identical to uniform weighting
```

### Example 3: CC Objective (No Sigma)

```
User command:
  simple_refine3D input.simple ...objfun=cc ml_reg=yes...

Sequence:
  1. params%new() → objfun='cc' → l_ml_reg = false (gate fails)
  2. calc_pspec_distr NOT called (no sigma calculation)
  3. Per-particle: cc_objfun != OBJFUN_EUCLID → calc_sigma2() skipped
  4. Reconstruction: uniform Fourier weights (l_ml_reg=false)
  Result: No sigma files generated; ml_reg setting ignored
```

---

## 11. Summary Table

| Aspect | Detail |
|--------|--------|
| **When calculated** | Always when `objfun='euclid'` |
| **When applied** | Only when `ml_reg='yes'` AND `objfun='euclid'` |
| **Primary input** | Per-particle Euclidean residuals (reference - particle)² |
| **Aggregation** | Group-based (by stack ID or global) per even/odd |
| **Weighting effect** | Wiener-filter: high-sigma (noisy) particles down-weighted; low-sigma (clean) particles up-weighted |
| **Ownership** | Builder object (`b_ptr%esig` in strategies, `build%esig` in commanders) |
| **Lifecycle** | Phase 1 (initial pspec) → Phase 2 (per-iteration) → Phase 3 (group average) → repeat |
| **Files** | Binary: `sigma2_part_k.dat` per partition; Starfile: `sigma2_group_iter_n.star` per iteration |
| **Special cases** | Disabled at ICM stage; nano3D uses global only; continuation carries over all partitions |
| **Known issue** | FSC helper missing automasking branch (low priority) |

---

## 12. References & Related Policies

- **Automasking Policy**: [automasking_policy_draft.md](automasking_policy_draft.md) — complementary mask-based refinement control
- **Builder Refactoring**: [euclid_build_new_policy.md](euclid_build_new_policy.md) — sigma ownership model
- **Objective Function Semantics**: `simple_type_defs.f90` (ENUM_OBJFUN), `simple_eul_prob_tab.f90` (distance conversions)
- **ML Regularization Variants**: `simple_convergence.f90` (logging), `simple_abinitio_utils.f90` (ab initio gating)

---

**Document History**:
- 2026-03-05: Initial draft (comprehensive analysis from codebase archaeology)
- 2026-03-06: Document review and minor additions
- (Future: Field test results, user feedback, deviations discovered)
