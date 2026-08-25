# Review of `drop_legacy_box_division`

## Status

**Static branch review — changes requested before merge.**

- Review date: 2026-08-23
- Branch: `drop_legacy_box_division`
- Reviewed range: `fdf2dee3c36e33c6ef604ff617f016ab20bf8655..9b82e8760`
- Branch remote state at review: local HEAD matched
  `origin/drop_legacy_box_division`.

This document reviews all work added to `drop_legacy_box_division` after its
merge base with `master`. It is independent of the implementation record in
[`drop_legacy_box_division.md`](drop_legacy_box_division.md). The implementation
record supplies useful experimental history, but its reported builds and test
runs were not repeated as part of this review.

## 1. Recommendation

Do not merge the branch in its reviewed state.

The central new convention is internally consistent for newly reconstructed
maps in the principal gridding and PCG paths. However, the branch can silently
damage refinements that start from maps written under the old convention. Two
Flex PCA paths also omit the new deapodization step. Test coverage does not
currently enforce the behavior that the branch intends to establish.

Required disposition before merge:

1. Prevent an old-convention starting map from entering the first alignment
   iteration without conversion or an explicit compatibility decision.
2. Apply the new gridding deapodization in both affected Flex PCA paths.
3. Update the conventional-gridding test helper to the new output convention.
4. Convert `test=rec3D_backends` from a measurement tool into a gated test.
5. Remove, gate, sample, or fuse the unconditional Euclidean diagnostics in the
   particle hot path.
6. Aggregate diagnostic statistics across distributed parts before reporting
   them as iteration statistics.

## 2. Review scope

The branch contains five commits and changes 35 files, with 1,103 insertions
and 185 deletions.

| Commit | Purpose reviewed |
|---|---|
| `1adaedad5` | Add per-iteration Euclidean reference/particle scale and objective diagnostics. |
| `285f55d2f` | Add the same-input gridding-versus-PCG reconstruction comparison. |
| `118388291` | Add native-period discrete Kaiser-Bessel deapodization and share the envelope implementation with PCG. |
| `6894d7687` | Record deapodization measurements and refine the backend comparison. |
| `9b82e8760` | Remove the matched reconstruction `/box` and projection `*box` factors and add an old-map warning. |

The review traced the affected behavior through:

- projector expansion and reprojection;
- gridding half-map and merged-map finalization;
- PCG solve output and warm starts;
- distributed volume assembly and nonuniform-filter inputs;
- automasking and FSC consumers;
- Flex PCA reconstruction and basis generation;
- Euclidean sigma estimation and diagnostics;
- production test helpers and the new high-level backend comparison.

## 3. Severity convention

- **P1 — blocking:** can invalidate a production workflow or scientific result
  without a reliable safeguard.
- **P2 — important:** causes incorrect output, invalid validation evidence, or a
  material production regression in a bounded path.
- **P3 — follow-up:** behavior is misleading or incomplete but does not by
  itself corrupt the primary reconstruction result.

## 4. Findings

### 4.1 P1 — old starting maps can collapse refinement before detection

#### Evidence

The projector now expands the volume's Fourier coefficients without the former
box-size multiplier:

- [`src/main/image/simple_projector.f90:41`](../../src/main/image/simple_projector.f90#L41)
- [`src/main/image/simple_projector.f90:90`](../../src/main/image/simple_projector.f90#L90)

Maps produced before this change contain the matching historical `1/box`
reconstruction scale. Such a map therefore reprojects approximately `box`
times below the particle scale in the new code.

The new warning is evaluated only when Euclidean sigma output is written:

- [`src/main/simple_euclid_sigma2.f90:287`](../../src/main/simple_euclid_sigma2.f90#L287)
- [`src/main/simple_euclid_sigma2.f90:342`](../../src/main/simple_euclid_sigma2.f90#L342)

By then, the current iteration has already aligned particles against the weak
reference. The branch's own implementation record documents a case where the
first iteration collapsed FSC and the references remained near zero in later
iterations:

- [`doc/implementation_notes/drop_legacy_box_division.md:450`](drop_legacy_box_division.md#L450)

The warning text says that the first sigma update recovers. The recorded result
shows that this is not a safe general claim.

#### Impact

Restarting or continuing an existing project can silently produce failed
orientations, collapsed FSC, and persistently weak references. A warning after
the damaging alignment is not a compatibility mechanism.

The unconditional convention switch also conflicts with the repository rule
that changed behavior must be explicit and default-off:

- [`.codex/AGENTS.md:3`](../../.codex/AGENTS.md#L3)

#### Required correction

Detect the convention before reference preparation and alignment. Then do one
of the following:

- convert an identified old map to the new scale;
- reject it with an actionable error before alignment; or
- place the new convention behind an explicit compatibility key until project
  or map metadata can identify the convention reliably.

A statistical warning after matching may remain as a diagnostic, but it must
not be the primary safeguard.

#### Acceptance condition

An old-convention restart must either stop before the first matcher invocation
or produce first-iteration references and alignment diagnostics equivalent to
an explicitly `box`-scaled control. The behavior must also be verified for a
new-convention map so that it is not scaled twice.

### 4.2 P2 — two Flex PCA paths do not apply deapodization

#### Evidence: covariance-column representatives

`columns_to_real_representatives` creates the new inverse-envelope image and
passes it to `realize_hermitian_volume`:

- [`src/main/flex/simple_flex_pca_columns.f90:1998`](../../src/main/flex/simple_flex_pca_columns.f90#L1998)
- [`src/main/flex/simple_flex_pca_columns.f90:2019`](../../src/main/flex/simple_flex_pca_columns.f90#L2019)

The callee says that it deapodizes the volume, but `gridcorr_img` is an unused
dummy argument. The routine inverse-transforms and masks without multiplying by
the correction:

- [`src/main/flex/simple_flex_pca_columns.f90:2034`](../../src/main/flex/simple_flex_pca_columns.f90#L2034)
- [`src/main/flex/simple_flex_pca_columns.f90:2050`](../../src/main/flex/simple_flex_pca_columns.f90#L2050)

#### Evidence: coupled M-step basis volumes

The coupled M-step compresses and inverse-transforms the even and odd
reconstructors, then proceeds directly to FSC, merge, filtering, and masking:

- [`src/main/flex/simple_flex_pca_columns.f90:3681`](../../src/main/flex/simple_flex_pca_columns.f90#L3681)
- [`src/main/flex/simple_flex_pca_columns.f90:3691`](../../src/main/flex/simple_flex_pca_columns.f90#L3691)
- [`src/main/flex/simple_flex_pca_columns.f90:3698`](../../src/main/flex/simple_flex_pca_columns.f90#L3698)

No inverse-envelope correction is applied to either half-map or the merged
basis volume.

#### Impact

The affected Flex basis volumes retain the Kaiser-Bessel real-space envelope.
This introduces a systematic radial taper before FSC and orthonormalization and
makes these paths inconsistent with production gridding and the corrected
weighted-state Flex reconstruction.

#### Required correction

Apply the same native-lattice inverse envelope after inverse FFT and before FSC,
masking, energy calculation, or basis orthonormalization. Avoid leaving a dummy
correction argument that the implementation does not consume.

#### Acceptance condition

A centered synthetic volume passed through each Flex reconstruction/finalizing
path must recover the same radial profile, up to a fitted global gain, as the
corrected production gridding path. The test must fail when the correction
multiplication is deliberately removed.

### 4.3 P2 — the conventional-gridding test helper uses the retired scale

#### Evidence

The production test helper still finalizes conventional gridding with `/BOX`
and does not apply the new inverse envelope:

- [`production/tests/simple_continuous_3D_pcg_refinement_halfset_gridding.f90:68`](../../production/tests/simple_continuous_3D_pcg_refinement_halfset_gridding.f90#L68)

The extended half-set matrix uses this result as its conventional gridding
baseline and requires the best PCG raw L2 result to beat it:

- [`production/tests/simple_continuous_3D_pcg_refinement_halfset_matrix.f90:274`](../../production/tests/simple_continuous_3D_pcg_refinement_halfset_matrix.f90#L274)
- [`production/tests/simple_continuous_3D_pcg_refinement_halfset_matrix.f90:297`](../../production/tests/simple_continuous_3D_pcg_refinement_halfset_matrix.f90#L297)

#### Impact

The baseline no longer represents production gridding. Its raw amplitude is
approximately `1/BOX` of the new convention and it retains the old envelope.
The scale-sensitive PCG-versus-gridding assertion is therefore invalid and can
become artificially easy to satisfy.

#### Required correction

Finalize the helper through the same convention and deapodization contract as
production `reconstructor_eo`, preferably by reusing the production finalizer
or a shared numerical helper rather than duplicating the sequence.

#### Acceptance condition

For identical accumulators, the helper and production gridding path must agree
in shell amplitude and radial profile to a fixed numerical tolerance. The PCG
raw-error comparison must use that corrected baseline.

### 4.4 P2 — `test=rec3D_backends` measures but does not test

#### Evidence

The new command explicitly states that it enforces no thresholds:

- [`src/main/commanders/test/simple_commanders_test_highlevel.f90:2072`](../../src/main/commanders/test/simple_commanders_test_highlevel.f90#L2072)
- [`src/main/commanders/test/simple_commanders_test_highlevel.f90:2082`](../../src/main/commanders/test/simple_commanders_test_highlevel.f90#L2082)

Its final section logs the measured shell and radial ratios and the expected
behavior, but it does not issue an error when the measurements violate those
expectations:

- [`src/main/commanders/test/simple_commanders_test_highlevel.f90:2340`](../../src/main/commanders/test/simple_commanders_test_highlevel.f90#L2340)

#### Impact

A restored box factor, a missing deapodization multiplication, or a large
backend amplitude disagreement can still end in a normal test stop. The main
behavior introduced by the branch has no enforceable end-to-end regression
gate.

#### Required correction

Turn the measured properties into explicit assertions with PASS/FAIL output and
a nonzero failure path. At minimum, gate:

- the median shell amplitude ratio in the valid agreement band;
- gridding-versus-PCG FSC over the declared comparison band;
- the normalized radial-ratio range inside the supported mask; and
- the number of valid shells and radial bins used by the comparison.

Thresholds must be based on a checked-in deterministic fixture or a documented
validated reference run. A test must not pass because no bins qualified.

#### Acceptance condition

The test must pass for the validated implementation and fail under each of two
intentional mutations: restore the legacy box factor and omit gridding
deapodization.

### 4.5 P2 — unconditional diagnostics add work to the Euclidean hot path

#### Evidence

The diagnostic arrays are allocated whenever `euclid_sigma2` is constructed:

- [`src/main/simple_euclid_sigma2.f90:92`](../../src/main/simple_euclid_sigma2.f90#L92)

For every particle, `calc_sigma2` allocates three temporary shell arrays and
always requests reference power, particle power, and the normalized objective:

- [`src/main/simple_euclid_sigma2.f90:263`](../../src/main/simple_euclid_sigma2.f90#L263)
- [`src/main/simple_euclid_sigma2.f90:268`](../../src/main/simple_euclid_sigma2.f90#L268)

`gen_sigma_contrib` then performs separate reference-power and particle-power
reductions plus a weighted shell loop that recomputes residual and particle
power before the original sigma contribution reduction:

- [`src/main/pftc/simple_polarft_corr.f90:1991`](../../src/main/pftc/simple_polarft_corr.f90#L1991)
- [`src/main/pftc/simple_polarft_corr.f90:2001`](../../src/main/pftc/simple_polarft_corr.f90#L2001)
- [`src/main/pftc/simple_polarft_corr.f90:2018`](../../src/main/pftc/simple_polarft_corr.f90#L2018)

#### Impact

Temporary allocation and several additional full-band reductions now occur in
a per-particle, per-iteration production path. The diagnostic is useful during
the convention transition, but its cost is paid even when no report or
compatibility investigation was requested.

#### Required correction

Use one or more of these approaches:

- protect the diagnostics with an explicit default-off key;
- sample a bounded subset of particles;
- retain running aggregates instead of per-particle arrays; and
- fuse power and residual accumulation with the existing sigma calculation.

#### Acceptance condition

With diagnostics disabled, the particle path must execute the pre-branch
calculation shape without diagnostic allocation or extra power reductions.
With diagnostics enabled, reported values must remain equivalent to the
current definitions.

### 4.6 P3 — distributed diagnostics report only part 1

#### Evidence

The diagnostic storage is local to `fromp:top`, and reporting returns
immediately unless the current process is part 1:

- [`src/main/simple_euclid_sigma2.f90:86`](../../src/main/simple_euclid_sigma2.f90#L86)
- [`src/main/simple_euclid_sigma2.f90:317`](../../src/main/simple_euclid_sigma2.f90#L317)

There is no reduction or consolidation of diagnostic values from the other
parts before the quantiles and particle count are printed.

#### Impact

The output is presented as a once-per-iteration diagnostic, but its `NPTCLS`,
quantiles, and warning decision describe only part 1's assigned subset. This can
misrepresent a heterogeneous distributed particle population.

#### Required correction

Aggregate sufficient statistics or diagnostic samples across all parts before
part 1 reports. If aggregation is deliberately not implemented, label the
output explicitly as part-1-only and do not use it as a global compatibility
decision.

#### Acceptance condition

A partitioned deterministic fixture must report the same diagnostic result as
the equivalent single-process particle set, within the selected quantile or
sampling tolerance.

## 5. Paths that were consistent in static review

The following areas did not expose an additional blocking issue:

- The principal `reconstructor_eo` half-map and merged-map paths apply the new
  deapodization before their downstream FSC, automask, and nonuniform-filter
  consumers.
- The PCG output and warm-start paths consistently remove the retired map-scale
  conversion.
- The weighted-state Flex reconstruction applies the new inverse envelope
  before filtering and masking.
- Otsu-based automasking is relative to the image distribution and therefore
  does not introduce a new direct dependence on the removed global box factor.
- All reviewed `expand_cmat` call sites use the new argument-free API.
- Distributed nonuniform-filter source halves are read after
  `reconstructor_eo` has already applied deapodization; the removed second
  correction in that path is consistent with the new ownership boundary.

These observations do not replace runtime verification.

## 6. Validation performed

The review used read-only source inspection and lightweight static checks:

- established the branch merge base and complete commit range;
- inspected every changed file and every branch commit;
- traced changed APIs and their call sites;
- searched for remaining legacy box-scale operations and correction helpers;
- checked the primary reconstruction, projection, Flex, nonuniform, and test
  consumers; and
- ran `git diff --check` over the reviewed range; it reported no whitespace
  errors.

No compilation, linking, CMake build, or test binary execution was performed.
Repository instructions reserve those steps for the user unless explicitly
requested. Consequently, statements in the implementation record that tests
or recovery suites pass are historical evidence, not independently observed
results of this review.

## 7. Suggested validation sequence after correction

1. Run editor or language-server diagnostics on every modified Fortran file.
2. Compile the focused reconstruction, Flex, and high-level test targets.
3. Run the corrected conventional-gridding and PCG half-set tests.
4. Run the gated `rec3D_backends` fixture and verify its intentional-failure
   mutations.
5. Compare new-map and old-map restart controls through the first refinement
   iteration.
6. Run the affected Flex covariance-column and coupled-M-step fixtures and
   compare their radial profiles with production gridding.
7. Run a distributed Euclidean fixture and compare its report with the
   equivalent single-process result.
8. Measure Euclidean iteration time with diagnostics disabled and enabled.
9. Run the established refine3D recovery suite, including nonuniform filtering
   and automasking configurations.

Record the exact commands, commit, platform, result, and observed numerical
gates. Do not report Linux, BOX, or other platform results unless their output
was observed directly.
