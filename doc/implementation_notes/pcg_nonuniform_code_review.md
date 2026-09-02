# Code review: PCG reconstruction with nonuniform regularization

**Review status:** Corrections implemented in the working tree; latest runtime fix pending rebuild/validation
**Reviewed revision:** `5e9933ffe6dee71e3e464589c74fa991980799ce`
(`more automsk CSI`, 2026-09-02)  
**Scope:** `rec_backend=pcg` with
`filt_mode=nonuniform|nonuniform_lpset`, including `nu_refine`, `automsk`,
`envfsc`, reconstruction support, matching-reference preparation, and final
postprocessing.

This is a static source review. No compilation, linking, or runtime tests were
performed, in accordance with repository policy.

The findings below preserve the reviewed-revision snapshot. The corrective
implementation removes post-hoc PCG NU filtering and masking, makes the
conservative solve support independent of `automsk`, enforces the requested
`envfsc` solve matrix once a reconstruction-derived mask exists, bootstraps a
missing-reference base solve on the sphere before deriving replay support,
makes the envelope background a fixed coarsest-label Potts boundary condition
without changing the unconstrained null-readiness statistic, and carries the
evidence-derived LP handoff through both PCG execution routes.

## Required contract

SIMPLE has two envelope concepts with different ownership:

1. The **conservative density envelope** may be used as the PCG solve support.
2. The more aggressive **NU-evidence envelope** may shape the local-resolution
   field used to construct the real-space PCG precision prior, but it must
   never be installed as the PCG solve support.

Neither envelope may be multiplied into an iterative reference volume. When NU
refinement is active, the appropriate envelope changes the filter field by
forcing its background to the coarsest, lowest-resolution bank member:

- `automsk=yes`: use the NU-evidence envelope;
- `automsk=no`: use the conservative density envelope in its place.

For the PCG backend, that filter-field evidence must enter reconstruction only
as the Potts-regularized, spatially varying precision prior `Q_NU`. No
post-hoc NU filtering or post-hoc masking may alter the primary PCG maps or the
references used by the next refinement iteration.

The expected solve-support matrix is:

| `envfsc` | Base/unregularized PCG solve | Regularized PCG replay |
| --- | --- | --- |
| `yes` | conservative density envelope | conservative density envelope |
| `no` | spherical support | conservative density envelope |

In this note, `envfsc` is assumed to be the intended spelling of `envfs` in the
review request.

## Findings

### [P1] Plain `nonuniform` applies a second, post-hoc NU filter after `Q_NU`

The `Q_NU` branch in
[`filter_pcg_nonuniform_maps`](../../src/main/commanders/simple/simple_commanders_rec_distr.f90#L808)
returns early only for `nonuniform_lpset`. For plain `nonuniform`, lines
880-886 deliberately continue into `nu_filter_vols`; lines 1020-1032 then
write `_nu_filt` half maps, an averaged `_nu_filt` map, and `_nu_locres`.

This gives plain `nonuniform` two NU mechanisms:

1. an in-solve `Q_NU` precision; and
2. a post-hoc bank filter used to create matching references.

That violates the mode contract and makes the two NU topologies scientifically
different for reasons unrelated to their intended matching topology. The PCG
path should skip `nu_filter_vols` for both modes. `nonuniform_lpset` may retain
only its evidence-derived scalar LP handoff; plain `nonuniform` should use the
primary per-half `Q_NU` maps directly.

### [P1] Envelope background is encoded as “no support in any band,” not the coarsest bank member

The ordinary NU filter implements the requested behavior correctly:
[`optimize_nu_cutoff_finds`](../../src/main/nu_filt/simple_nu_filter_bank.f90#L186)
sets envelope-background voxels to candidate 1 before and after ordered-label
smoothing (lines 218-235).

The compact evidence state used by PCG does something different. In
[`build_nu_evidence_state`](../../src/main/nu_filt/simple_nu_filter_evidence.f90#L197),
the envelope-background override occurs after the Potts pass and sets
`selected_label=0`, `selected_cutoff=0`, and every `band_support` value to zero
(lines 207-218). During expansion,
[`expand_nu_evidence_band_weights`](../../src/main/nu_filt/simple_nu_filter_evidence.f90#L540)
computes `band_w = 1 - band_support` (line 562), producing maximum precision in
every band.

This is not equivalent to selecting the coarsest member. The first PCG band
covers the coarse detail at or below the first resolution boundary, as
documented by
[`ensure_nu_band_index`](../../src/main/volume/simple_reconstructor_pcg.f90#L696).
A coarsest-bank background should preserve that coarse band and penalize only
unsupported finer bands. Conceptually, its support should be
`[1, 0, 0, ...]`, yielding lack-of-evidence weights `[0, 1, 1, ...]`; the
current null encoding yields `[1, 1, 1, ...]` and suppresses all non-DC
background detail.

The background constraint should therefore be represented as the coarsest
candidate in the same Potts-regularized evidence field, not replaced with a
post-Potts null label. A focused test should assert the expanded per-band
weights for a clamped background voxel.

### [P1] The density solve-support policy is route-dependent and incorrectly gated by `automsk`

[`build_pcg_state_support`](../../src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90#L551)
returns immediately when `automsk=no` (line 568), even though `automsk` chooses
the filter-field envelope and should not disable the conservative density
support required by PCG. It also falls back to a sphere for a missing reference
or a `startvol` reference (lines 570-573).

The distributed solve partly implements the intended `envfsc` matrix:
[`reduce_solve_state_half`](../../src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90#L2335)
uses the density support when
`l_state_support .and. (l_ml_solve .or. params%l_envfsc)` (lines 2352-2357).
This produces the requested behavior only when `automsk/=no`, a compatible
lagged reference exists, and no explicit mask overrides it.

The shared route never supplies a density support at all. Its base solve calls
`set_pcg_solve_support(pcgop, params)` at line 1080 and its replay does the same
at line 1195, so both use the spherical fallback unless `pcg_mskfile` is set.

The effective matrix is therefore:

| Route/configuration | `envfsc=yes` | `envfsc=no` |
| --- | --- | --- |
| distributed, `automsk=yes`, valid prior reference | density for both solves | sphere for base; density for replay |
| distributed, `automsk=no` | sphere for both solves | sphere for both solves |
| shared | sphere for both solves | sphere for both solves |
| missing or `startvol` prior reference | sphere for both solves | sphere for both solves |

Density-support construction must be independent of `automsk` and available to
both shared and distributed owners. The base/replay choice should be made only
from `envfsc` and solve phase. Bootstrap behavior needs an explicit policy that
does not silently violate `envfsc=yes`.

### [P1] PCG skips envelope-corrected FSC based on backend name, even when the solve was spherical

[`calculate_halfmap_diagnostics`](../../src/main/volume/simple_halfmap_diagnostics.f90#L45)
sets `l_envfsc_preproc` false for every PCG run (line 68), on the assumption
that the density envelope already constrained the solve. The support defects
above make that assumption false for the shared path, `automsk=no`, bootstrap,
and missing-reference fallbacks.

Consequently, `envfsc=yes` can silently produce a spherical-support FSC while
the configuration claims density-envelope semantics. Once support routing is
fixed, avoiding a second envelope multiplication for PCG is correct. Until
then, this decision cannot be inferred from `rec_backend` alone; it must follow
the support actually installed for the base pair or fail when the promised
support is unavailable.

### [P2] PCG+NU configurations can bypass `Q_NU`

The dynamic default installs `pcg_nu_lambda_rel=0.1` only when the Euclidean ML
replay is active. Otherwise
[`simple_parameters_phases`](../../src/main/params/simple_parameters_phases.f90#L902)
warns that `Q_NU` is disabled (lines 928-930). An explicit
`pcg_nu_lambda_rel=0` is also documented as selecting the ordinary `P_tau`
control (lines 906-907), and
[`validate_nu_replay_request`](../../src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90#L583)
accepts zero.

Those controls were useful for development, but they violate the production
contract that `Q_NU` is the sole NU regularization mechanism on PCG. Production
`rec_backend=pcg` plus either NU filter mode should require a compatible
Euclidean regularized replay and active `Q_NU`, or reject the configuration.
If strength-zero and `P_tau` controls are retained, they should be isolated to
an explicit test/development route and must not trigger post-hoc NU filtering.

### [P2] Matcher input selection still treats `_nu_filt` files as the NU contract

[`prepare_refvol`](../../src/main/strategies/search/simple_matcher_refvol_utils.f90#L230)
first looks for `_nu_filt` files in every nonuniform mode (lines 242-267). If
they are absent, it treats the primary half maps as raw fallback inputs and
applies an FSC-optimal or global low-pass filter (lines 332-359).

After removing post-hoc PCG NU filtering, primary PCG half maps are not an
error fallback: they are the `Q_NU`-regularized matching references. The
matcher needs an explicit PCG/QNU path that selects those maps as authoritative
and does not apply another NU, FSC-optimal, or envelope filter. This must be
consistent in shared and distributed execution.

The nested `mask_matching_reference` helper is correct: lines 365-378 apply
only the spherical soft support and explicitly prohibit density- or
NU-envelope multiplication.

### [P2] `pcg_mskfile` can replace the conservative density support with an arbitrary mask

[`set_pcg_solve_support`](../../src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90#L522)
gives `pcg_mskfile` precedence over the state density envelope (lines 534-540).
If the rule that only the conservative density envelope may constrain a PCG
solve is absolute, this public override violates it.

Either remove the override from production NU routes, validate that it is a
conservative density mask with the required provenance, or clearly isolate it
as a development-only escape hatch. Under no circumstances should an
NU-evidence envelope be supplied through this parameter.

### [P2] Final `_pproc` output is still multiplied by the density envelope

[`postprocess_volume_from_files`](../../src/main/commanders/simple/simple_commanders_volops.f90#L391)
zeros the background and multiplies `vol_bfac` by
`automask3D_stateNN.mrc` when `envfsc` postprocessing is enabled (lines
391-410). This affects a derived `_pproc` artifact, not the primary PCG half
maps or iterative matching references.

If “no post-hoc masking in the PCG path” includes delivered derived maps, this
is a direct violation and should be disabled for PCG. If display/deposition
products are intentionally exempt, that exception should be documented and
the unmasked primary PCG map must remain the authoritative iterative and
scientific output.

### [P2] `refine3D_auto` does not actually make `envfsc` overridable

[`commander_refine3D_auto`](../../src/main/commanders/simple/simple_commanders_refine3D.f90#L93)
labels `envfsc` as an “overridable default” but unconditionally calls
`cline%set('envfsc','yes')` at line 105. Unlike the adjacent defaults, it does
not first test `cline%defined('envfsc')`.

If `cline%set` overwrites a supplied value, `envfsc=no` cannot exercise the
required base-sphere/replay-density path through `refine3D_auto`. This should
use the same guarded-default pattern as `filt_mode`, `nu_refine`, and
`automsk`.

### [P3] Existing policy and implementation notes contradict current source

The historical
[`nu_evidence_envelope_masking.md`](nu_evidence_envelope_masking.md) still says
matching references are multiplied by the NU-evidence envelope and describes
`envfsc=yes` as an unfinished hard-error path. Current source explicitly
prohibits reference-envelope multiplication and has a live density-mask
`envfsc` path. Other automasking/nonuniform policy sections contain both old
and new statements.

These documents should be reconciled after the implementation contract is
settled. Stale policy is especially risky here because several source comments
cite `pcg_priors.md` as normative behavior.

## Confirmed behavior

The following parts of the implementation agree with the intended design:

- `build_nu_replay_evidence` selects the NU-evidence envelope when
  `automsk=yes` and computes a conservative density envelope when
  `automsk=no`. It does not pass the NU-evidence envelope to
  `set_pcg_solve_support`.
- The evidence labels are produced by a Potts-regularized spatial model in
  [`build_nu_evidence_state`](../../src/main/nu_filt/simple_nu_filter_evidence.f90#L143).
- The replay attaches `Q_NU` through `set_nu_prior` in both shared and
  distributed regularized solves. `set_nu_prior` also rejects simultaneous
  `P_tau` and `Q_NU` attachment.
- [`apply_nu_precision`](../../src/main/volume/simple_reconstructor_pcg.f90#L797)
  implements a real-space, spatially varying precision: it isolates Fourier
  bands, transforms each band to real space, applies the per-voxel weight, and
  transforms the result back. The NU evidence therefore enters the PCG normal
  operator as a prior rather than as a phase-bearing target or post-hoc map
  multiplication.
- The PCG support mask is applied through the projected solve. The resulting
  support application is part of the estimator, not merely a cosmetic
  multiplication of an already-finished primary map.
- Iterative matching references currently receive only the broad spherical
  soft mask in `mask_matching_reference`; direct density- and NU-envelope
  multiplication has been removed there.

## Required validation after correction

At minimum, the corrected implementation should cover these cases with
focused route-level tests or deterministic diagnostics:

1. shared and distributed execution;
2. `nonuniform` and `nonuniform_lpset`;
3. `automsk=yes` and `automsk=no`;
4. `envfsc=yes` and `envfsc=no`;
5. bootstrap/start-volume and later-iteration references;
6. base and regularized solves, asserting the exact installed support source;
7. background evidence expansion, asserting coarsest-band preservation and
   finer-band suppression;
8. active `Q_NU`, with no post-hoc `nu_filter_vols`, no `_nu_filt` dependency,
   and no envelope multiplication of primary or matching-reference maps;
9. rejection or explicit development isolation of `pcg_nu_lambda_rel=0`,
   missing ML replay, and `pcg_mskfile` overrides.
10. a density-supported half pair that is exactly zero outside the conservative
    solve envelope but whose NU evidence is evaluated on the broader spherical
    support; its radial whitening profile and explicit null unary must remain
    finite.

## Runtime follow-up: hard-supported pair null-unary overflow

A user-side run reached the intended PCG topology: spherical bootstrap when no
reference existed, in-solve `Q_NU`, no post-hoc NU filtering, then a
lag-one-density-supported base pair on iteration 1. Evidence construction for
that second pair failed at the explicit zero-predictor unary with
`IEEE_INVALID_FLAG`/`IEEE_DIVIDE_BY_ZERO`.

The radial noise estimator had treated the exact zero/zero voxels created by
the narrower PCG solve support as zero-noise observations, even though the NU
objective is evaluated on the broader spherical domain. This can collapse a
shell scale. The null predictor is uniquely exposed because it measures the
full map amplitude; ordinary filtered candidates compare two similar maps and
can remain finite.

The correction keeps the evidence domain broad but omits exact zero/zero
boundary samples from the whitening fit, nearest-fills unsupported radial
shells as before, validates the resulting profile, and evaluates standardized
Huber residuals in double precision with a numerically safe finite saturation
well above meaningful unary differences. The focused NU-envelope test now
includes a hard-supported-pair regression over a broader NU sphere.

The latest correction has only static validation in this workspace. Rebuild
and focused/runtime reruns remain outstanding for the user.

## Follow-up: default-on `nu_refine` hardening

`refine3D_auto` intentionally retains `nu_refine=yes`; staged abinitio3D
continues to pin `nu_refine=no`. The latter is heavily used, so its PCG
evidence state is now an explicit compatibility boundary: the eight static
signal candidates keep their integer Potts coordinates and unit softmax
measure, the four-band ladder is unchanged, and the raw-finest matching
handoff is unchanged.

Only an evidence bank containing accepted extension candidates takes the new
spacing-aware representation. The static Potts coordinates remain exactly
1 through 8, while extension coordinates continue in Fourier-shell distance
normalized by the finest static interval. Signal-candidate posterior terms are
integrated with Voronoi cell widths normalized to the static bank's total
signal mass. Candidate thinning or inserting densely spaced probes therefore
cannot create support solely through label multiplicity.

PCG adaptive shell discovery is additionally capped at two Fourier shells
beyond the current evidence-pair FSC=0.143 crossing. The adaptive matching
handoff uses the established 5% assigned-support statistic and explicitly
retains the same two-shell FSC headroom, rather than allowing one extreme voxel
to set the global matching bandwidth. The cache identity includes the
`nu_refine` mode so a mode transition cannot ride incompatible evidence.

With `automsk=yes`, the NU envelope is deliberately the preliminary
static-bank boundary. It is fixed before shell challenges so adaptive evidence
cannot select its own spatial support. Accepted candidates refine `Q_NU`
inside that boundary. Documentation that previously claimed the envelope
always described the final accepted bank has been corrected.

The focused NU-envelope test now asserts the static candidate count, band
count, and absence of adaptive geometry, in addition to its existing evidence
and hard-support checks. Compilation and runtime validation remain with the
user.

## Follow-up: native-grid Q_NU calibration

The first concurrent abinitio3D run completed normally and verified the
compute-only even/odd parallel solve boundary. Its final reconstruction also
showed a 48--50% suppression jump when the box changed from the staged crop to
native sampling. That change is explained by the different downscaling/grid:
the Q_NU plant gain changes, so the reduced-grid `lambda_rel` is not a valid
native-grid operating point.

The final bootstrap now treats its first automatic-Q_NU solve as a current-grid
calibration rather than shipping it. It preserves the dataset-level suppression
target, clears the obsolete reduced-grid response during the unregularized
sigma pass, measures suppression at the normal cold-start strength on the
current grid, and performs a second regularized solve after the existing
auto-lambda controller adapts from that measurement. Postprocessing and output
registration occur only for the calibrated solve. User-pinned lambda values
remain single-pass and unchanged.
