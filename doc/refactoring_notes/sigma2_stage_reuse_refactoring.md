# Sigma2 reuse across refine3D stages

## Status and decision

The expensive image-power-spectrum bootstrap should run once, before the first
refine3D stage that needs sigma2. Later stages should reuse the residual sigma2
model produced by the preceding matcher.

The central design decision is:

> Preserve sigma2 on its existing native particle-image grid. Do not regrid the
> persistent artifacts when `box_crop` or `smpd_crop` changes.

This follows the current data model and an existing, tested consumer path:

- `calc_pspec` calculates and writes spectra with the original `params%box` and
  `params%smpd`.
- `euclid_sigma2%new(..., params%box)` allocates the same persistent shell range,
  `[1, fdim(params%box)-1]`, for matching and reconstruction.
- stage-specific `box_crop` and `smpd_crop` select how much of that model is
  consumed; they do not redefine the stored model.
- `image%gen_fplane4rec` already calls
  `simple_math_ft::upsample_sigma2` to interpolate the selected sigma2 range onto
  the padded Fourier-plane grid used by 3D reconstruction.

Consequently, a crop change is neither a compatibility failure nor a reason to
invent a second interpolation algorithm at stage handoff. The reuse work is
artifact validation, particle/group rematerialization when needed, and lifecycle
control so that `calc_pspec` is skipped.

## 1. Why the current stage bootstrap is wrong

`calc_pspec` rereads particle stacks and estimates sigma2 from masked image power
spectra. This is useful as an initial bootstrap, but it is both expensive and a
weaker estimate than the orientation-conditioned residual spectra already
produced by the matcher.

At the end of an established refine3D stage, the useful state is:

- `sigma2_noise_partN.dat`: the newest per-particle residual spectra, written by
  the matcher;
- `sigma2_it_<iter>.star`: the grouped model written by
  `calc_group_sigmas` at the start of an iteration.

The partition files are preferred because the final group STAR can be one
matcher update behind them. Re-estimating from image power spectra at every stage
throws away that information and creates the long `DISTR_CALC_PSPEC` barrier.

## 2. The native sigma2 grid

### 2.1 Persistent representation

The persistent sigma2 representation is tied to the original particle images:

```text
native shell k  <->  params%box, params%smpd
storage range   =    1 : fdim(params%box)-1
```

This is visible in the current producers and consumers:

- `simple_calc_pspec_strategy.f90` builds full-size images, obtains `nyq` from
  `build%img`, and writes `[1, nyq]` records;
- `simple_euclid_sigma2.f90` constructs the sigma object with `params%box`;
- refine3D matching derives its active Fourier limits in the same original-image
  shell convention.

Changing `box_crop` and `smpd_crop` between stages therefore does not change the
meaning or dimensions of a stored sigma2 record. The staged crop normally
preserves the physical field of view, while exposing a different usable
frequency range.

### 2.2 Existing reconstruction-plane interpolation

The 3D reconstruction path already owns conversion from the selected native
sigma range to the padded reconstruction plane:

```fortran
call upsample_sigma2(kfromto(1), sigma_nyq, sig2arr, &
    &fplane%nyq, sigma2_noise)
```

This call is in `image%gen_fplane4rec` in
`src/main/image/simple_image_ctf.f90`. The implementation is
`simple_math_ft::upsample_sigma2` in `src/utils/math/simple_math_ft.f90`.
`prep_imgs4rec` invokes `gen_fplane4rec` from its existing OpenMP particle loop,
so this conversion is already parallelized at the correct granularity and must
not be moved into a serial stage-initialization pass.

Its established contract is:

- the input upper index is the unpadded current-crop Nyquist shell;
- the output upper index is the padded Fourier-plane Nyquist shell;
- output shells are mapped to the input by the Nyquist ratio;
- interpolation is linear in sigma2 itself;
- frequencies below `kstart` use `sigma2(kstart)`;
- the upper endpoint uses the last input value;
- `nyq_out >= nyq` is a precondition.

The stage-reuse implementation must not duplicate, replace, or subtly alter this
behavior.

### 2.3 Newly active shells already exist

When a later stage increases the active frequency range, those shells are not
missing from the persistent model. The partition records span the native
full-size range. Shells outside the earlier matching band retain their existing
values; once activated, the normal matcher residual update replaces them.

No stage-transition tail extrapolation is required.

### 2.4 Scope boundary

A change to the original particle-image grid (`params%box` or `params%smpd`) is
different from a stage crop change. It implies re-extraction or a different
native data representation and is outside this refactoring.

`upsample_sigma2` must not be misused as a general physical-frequency mapper: it
does not accept sampling distances and supports Nyquist-ratio upsampling, not an
arbitrary native-grid conversion. If native-grid reuse is needed later, it
requires a separate, explicitly tested design.

## 3. Required lifecycle

Both shared-memory and distributed refine3D currently run
`calc_group_sigmas` before probabilistic alignment and matcher preparation.
`calc_group_sigmas` opens the current `sigma2_noise_partN.dat` files.

Therefore, stage reuse must complete before the first iteration and before the
first `calc_group_sigmas` call:

```text
resolve predecessor sigma artifacts
validate native grid and semantic provenance
materialize current partition files
skip calc_pspec
enter normal iteration loop
    calc_group_sigmas
    probabilistic alignment
    matcher and residual sigma update
```

A group STAR alone cannot be left for the matcher to expand later: the earlier
`calc_group_sigmas` call would fail first. Shared-memory and distributed refine3D
must call one common initialization helper so their behavior cannot diverge.

## 4. Initialization classes

### `EXACT_PARTICLE_REUSE`

Use when every source partition file has:

- the expected native shell bounds;
- exactly the target partition's particle range;
- compatible particle identity/order and source semantics.

Copy or retain the files transactionally. Do not transform their spectra.

Validation is per numbered filename, not only over the union of ranges. Reject
swapped partitions, overlaps, gaps, malformed headers, and ordering mismatches.

### `REPARTITIONED_PARTICLE_REUSE`

Use when valid per-particle spectra exist but `nparts` or partition boundaries
changed.

For each target partition:

1. resolve each target particle identity to its source record;
2. copy the complete native sigma2 record unchanged;
3. write a target `sigma2_noise_partN.dat` with current particle bounds.

This is bounded stream-binary I/O. It does not read particle stacks or perform
numerical interpolation.

### `GROUP_EXPANSION_INITIALIZATION`

Use only when compatible grouped spectra exist but the newer per-particle files
are unavailable or cannot be associated safely with current particle identities.

Expand the source even/odd group spectra into every current partition according
to current `eo` and `stkind`, and write all target partition files before
iteration entry. Preserve the source shell values unchanged.

This is a lower-quality fallback than particle reuse because the final group
model may lag the final matcher update, but it remains preferable to a new image
power-spectrum bootstrap.

### `FRESH_CALC`

Run `calc_pspec` only when:

- no predecessor exists, normally the first stage;
- no predecessor artifact can be interpreted safely; or
- the user explicitly requests `sigma2_init=refresh`.

For staged abinitio after the first sigma-producing stage, use strict reuse. An
unexpected incompatibility should fail with a precise message instead of
silently launching a multi-minute `DISTR_CALC_PSPEC` calculation.

## 5. Compatibility policy

### Required

- Source and target native `params%box` and `params%smpd` match.
- Binary shell bounds match the expected full native sigma range.
- The source refers to the same particle set, ordering, and image channel (`raw`
  versus `den`).
- `oritype` and active-particle semantics are compatible.
- Per-particle reuse has an unambiguous identity mapping.
- Group expansion has valid current `stkind` and `eo` values and an unambiguous
  source-group mapping.
- Values are finite and positive after the existing sanitization policy.
- Even/odd convention remains index 1 = even and index 2 = odd.

### Never reuse-invalidating

- `box_crop` or `smpd_crop` changes;
- active `kfromto` or low-pass changes within the stored native range;
- `nspace`, `nspace_sub`, or angular search changes;
- `prob_neigh_mode` changes;
- `update_frac`, `nsample`, or `nparts` changes.

`nparts` can require rematerialization, but it does not invalidate the sigma
model.

## 6. Provenance

The current binary header contains particle and shell bounds, but not all
semantic provenance. The staged abinitio controller should pass trusted source
information to the next refine3D stage:

- predecessor path;
- predecessor final iteration;
- native `box` and `smpd`;
- source `ptcl_src`, `sigma_est`, and `oritype`;
- particle-project identity/order confirmation.

For restart or reuse outside one controller invocation, persist the necessary
metadata in the group STAR `general` block. Before adding fields, make the fast
`read_sigma2_groups` reader key-driven rather than dependent on fixed line
positions.

Crop metadata may be reported for diagnostics, but it is not part of native-grid
compatibility.

## 7. Controller and strategy changes

### Controller

After each completed stage, retain the predecessor path, returned final
iteration, and sigma provenance before constructing the next refine3D command.

Use:

- first stage: `sigma2_init=auto`;
- later stage with predecessor: `sigma2_init=reuse`.

Modes:

- `auto`: reuse a supplied predecessor; bootstrap only when no predecessor
  exists;
- `reuse`: require valid reuse and fail clearly otherwise;
- `refresh`: force `calc_pspec` for controlled recovery or benchmarking.

### Common refine3D initializer

Implement one helper used by shared-memory and distributed strategies:

```text
initialize_sigma2_for_stage
    resolve predecessor
    validate native grid and semantics
    inspect binary/group headers
    classify exact / repartition / group expansion / fresh
    materialize target files transactionally
    return class and diagnostics
```

Artifact inspection and rematerialization belong close to
`simple_euclid_sigma2` and `simple_sigma2_binfile`. Stage orchestration belongs
to refine3D strategy initialization. Numerical interpolation remains owned by
`simple_math_ft::upsample_sigma2` and its existing reconstruction consumer.

Write temporary target files and promote them only after every required file has
been produced successfully. A failed attempt must not leave a partial set that
can poison a retry.

## 8. I/O and performance requirements

- Retain the existing `sigma2_noise_partN.dat` stream-binary format.
- Retain `sigma2_it_<iter>.star` for grouped spectra.
- Prefer the newest partition files over the older group STAR.
- Process one target partition at a time when repartitioning.
- Do not read particle stacks during reuse.
- Do not perform a serial per-particle interpolation pass.
- Preserve the existing OpenMP behavior of downstream readers and matchers.
- Use temporary files plus atomic rename for promotion.

`CALCPSPEC_FINISHED` remains owned by the actual distributed `calc_pspec`
workflow. Reuse should not synthesize that marker; worker dispatch starts only
after master initialization has completed.

## 9. Diagnostics

Print one concise record per stage, for example:

```text
SIGMA2 INIT: EXACT_PARTICLE_REUSE source_iter=73 target_iter=74 \
  native_grid=128x1.30 crop=88x1.89->112x1.49 partitions=8 calc_pspec=no
```

```text
SIGMA2 INIT: REPARTITIONED_PARTICLE_REUSE source_iter=73 target_iter=74 \
  native_grid=128x1.30 records=4600 partitions=8->12 calc_pspec=no
```

```text
SIGMA2 INIT: GROUP_EXPANSION_INITIALIZATION source_iter=73 target_iter=74 \
  native_grid=128x1.30 groups=4 records=4600 partitions=8 calc_pspec=no
```

```text
SIGMA2 INIT: FRESH_CALC reason=no predecessor sigma artifacts
```

Report the initialization class, source iteration, native grid, source and target
crop for context, record count, partition counts, and whether `calc_pspec` ran.

## 10. Failure policy

Reject reuse for:

- missing or malformed files;
- native `box` or `smpd` mismatch;
- unexpected shell bounds;
- particle identity/order or image-channel mismatch;
- incompatible `oritype`, even/odd convention, or group mapping;
- non-finite or non-positive spectra that cannot be sanitized;
- ambiguous predecessor stage.

Do not reject reuse for crop, active-resolution, angular-search, update-fraction,
or partition-layout changes.

Under `sigma2_init=reuse`, rejection is a hard error with the exact reason.
Under `auto`, absence of a predecessor may bootstrap, but a malformed explicitly
provided predecessor should fail rather than silently launch `calc_pspec`.

## 11. Implementation phases

### Phase 1 — stage reuse

- Add controller handoff metadata and `sigma2_init=auto|reuse|refresh`.
- Add the common initializer to shared-memory and distributed refine3D.
- Implement strict per-file validation.
- Implement exact reuse, repartitioned particle reuse, and group expansion.
- Materialize all current partition files before `calc_group_sigmas`.
- Add transactional promotion and concise diagnostics.
- Skip `calc_pspec` when reuse succeeds.

### Phase 2 — restart robustness

- Persist semantic provenance for reuse outside one controller invocation.
- Make the fast group STAR reader key-driven before extending its general block.
- Optimize source-record lookup if repartitioning benchmarks show it matters.

No new sigma interpolation algorithm is part of either phase.

## 12. Validation plan

### Unit tests

- `upsample_sigma2` remains the numerical oracle for reconstruction-plane sigma
  interpolation:
  - equal Nyquist sizes preserve input values;
  - padded output uses the established linear mapping;
  - values below `kstart` and at the upper endpoint follow the current boundary
    rules.
- A changed `box_crop`/`smpd_crop` with unchanged native `box`/`smpd` is accepted
  without modifying stored spectra.
- Exact partition validation rejects swapped files, overlaps, gaps, malformed
  headers, and particle-order mismatches.
- Repartitioning preserves every source spectrum bit-for-bit.
- Group expansion respects current `eo`, `stkind`, and global/group policy.
- `sigma2_init=refresh` forces bootstrap.
- `sigma2_init=reuse` rejects an unusable source without launching calc_pspec.

### Integration tests

1. Stage A produces residual sigma artifacts.
2. Stage B changes `box_crop`, `smpd_crop`, `nspace`, and `nspace_sub` while
   retaining native `box` and `smpd`.
3. Assert that `DISTR_CALC_PSPEC` is skipped and stored sigma values are
   unchanged before the first target iteration.
4. Verify reconstruction-plane sigma arrays against the existing
   `gen_fplane4rec`/`upsample_sigma2` behavior at both crop sizes.
5. Repeat with different `nparts`; assert repartitioned reuse and exact particle
   identity/value preservation.
6. Remove source particle files; assert group expansion creates every target
   partition before `calc_group_sigmas`.
7. Run the same cases through shared-memory and distributed refine3D.
8. Verify strict failure for a true native-grid mismatch.

Do not compare a reused residual model numerically with a fresh `calc_pspec`
model: they are scientifically different estimates. Validate preservation and
consumer behavior instead.

### Performance acceptance criteria

- Later stages do not run `DISTR_CALC_PSPEC`.
- Exact reuse is file-validation/copy cost only.
- Repartitioning is bounded stream I/O with no image reads and no numerical
  interpolation pass.
- Initialization does not introduce a long serial master phase.
- Existing `calc_group_sigmas`, probabilistic alignment, matcher, and
  reconstruction timings remain unchanged after initialization.

## 13. Review decisions

1. **Require equal stage crop size or sampling?** No. `box_crop` and `smpd_crop`
   are consumer settings, not persistent sigma-grid compatibility gates.
2. **Where is interpolation done?** In the existing
   `image%gen_fplane4rec -> simple_math_ft::upsample_sigma2` path when padded 3D
   reconstruction planes are prepared.
3. **Add stage-handoff interpolation?** No. Stored spectra remain on the native
   `params%box`/`params%smpd` grid.
4. **What about newly active frequencies?** Their values already exist in the
   full native records and are replaced by normal matcher residual updates once
   active.
5. **Can group-only reuse leave partition files absent?** No. Current iteration
   ordering requires them before `calc_group_sigmas`.
6. **Automatic fallback after stage 1?** No expensive silent fallback. Staged
   abinitio requests strict reuse and reports incompatibility.
7. **Preserve the binary format?** Yes.
