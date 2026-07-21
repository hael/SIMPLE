# Phase-shift CTF Policy

Status: implemented; code-level validation complete; scientific data
validation pending.

Last reviewed: 2026-07-20

## 1. Purpose and scope

This document defines SIMPLE's policy for fitting, storing, propagating, and
using an additional CTF phase shift. It applies equally to conventional,
Volta-phase-plate, laser-phase-plate, and externally imported cryo-EM data.

This is the supported domain model in the current codebase. It is not an
implementation proposal and does not define a separate phase-plate
reconstruction mode.

## 2. Core domain policy

SIMPLE represents phase-plate effects as an ordinary numerical term in the CTF
model:

```text
H(k) = sin(chi_defocus,Cs(k) + phi_amplitude_contrast + phshift)
```

The following rules are normative:

1. `phshift` is the additional CTF phase and is always part of the numerical
   CTF state.
2. `fit_phshift=yes|no` controls only whether CTF estimation optimizes
   `phshift`.
3. All later CTF-dependent routines consume the stored numerical `phshift`
   regardless of how the data were collected or how that value was obtained.
4. Zero is the identity/default phase. It is a real stored value, not an
   indication that a field or acquisition mode is absent.
5. Acquisition provenance is not a CTF execution switch. SIMPLE has no
   supported `phaseplate`, `l_phaseplate`, or device-kind processing flag.
6. No compatibility path for the removed legacy phase-plate flag is required.

Imported nonzero phase shifts are therefore honored without
`fit_phshift=yes`. Enabling `fit_phshift` affects CTF fitting and persists its
result for all downstream consumers.

## 3. Units, range, and defaults

Internal and project units are radians:

| Boundary | Unit and range |
|---|---|
| `ctfparams` | canonical radians in `[0, pi)` |
| `ctfvars` | canonical radians in `[0, pi)` |
| `ori`/`oris` and SIMPLE project fields | canonical radians in `[0, pi)` |
| optimized and scalar CTF routines | canonical radians in `[0, pi)` |
| RELION `rlnPhaseShift` | degrees on input and output |
| CTF-fit UI controls | degrees |
| plaintext particle import | controlled by `phshiftunit` |

The common orientation setter canonicalizes numerical phase modulo `pi`.
Values outside the canonical interval are normalized at the domain boundary,
not repeatedly inside optimized pixel loops.

The supported fit controls and defaults are:

```text
fit_phshift=no
phshift_min=0 degrees
phshift_max=180 degrees
phshift_step=10 degrees
```

The 0--180-degree interval reflects the phase ambiguity of power-spectrum CTF
fitting and agrees with the established RELION/cisTEM approach. The grid is an
initial search, not an assertion that a laser or Volta phase must be exactly
90 degrees. The user-facing 180-degree maximum denotes the ambiguity boundary;
the estimator excludes the exact endpoint from its grid and maps the continuous
optimizer's upper limit to the largest representable value below `pi`.

## 4. In-memory CTF contract

`ctfparams` is the complete caller-facing CTF parameter bundle. It contains
microscope parameters, defocus/astigmatism, `phshift`, the CTF mode, and the
phase-fitting policy used by the estimator.

`ctfvars` is the flattened coefficient/value bundle used by optimized CTF
routines. It also contains `phshift`. Constructing `ctfvars` from a `ctf`
object requires the numerical phase explicitly; optimized routines must not
reconstruct a phase-free bundle implicitly.

Performance-sensitive CTF routines must take a required phase argument or a
required phase-bearing value object. Optional phase arguments and per-pixel
feature branches are not permitted in these paths. Conventional callers pass
zero.

## 5. Project persistence and propagation

Every row in the following project segments carries `phshift`:

- `os_mic`
- `os_stk`
- `os_ptcl2D`
- `os_ptcl3D`

Project serialization materializes a missing phase as zero. Project validation
requires the field. A stored zero survives project copying, merging,
subsetting, and STAR conversion like any other CTF value.

The numerical phase follows the same mapping boundaries as defocus and
astigmatism:

| Mapping boundary | Required behavior |
|---|---|
| estimator result -> `os_mic` | Store the fitted canonical phase unconditionally, including a fitted zero |
| `os_mic` -> extracted stack and particles | Copy phase into `os_stk`, `os_ptcl2D`, and `os_ptcl3D` |
| imported `os_stk` -> particles | Copy stack phase to each generated particle unless an authoritative per-particle value is present |
| `os_ptcl2D` <-> `os_ptcl3D` | Preserve phase in full-record and manual mappings |
| subset, merge, chunk, and distributed assembly | Preserve the source row's phase while remapping project indices |
| project field -> `ctfparams` | Return canonical phase for mic, stack, particle-2D, and particle-3D queries |
| SIMPLE -> RELION STAR | Emit `rlnPhaseShift` in degrees on every relevant row |
| RELION STAR -> SIMPLE | Convert degrees to radians; materialize zero when the column is absent |

A consolidated stack may contain particles from source stacks with different
phases. In that case the aggregate `os_stk` row carries explicit zero and the
particle rows remain authoritative. Aggregate stack metadata must never
overwrite heterogeneous particle phases.

## 6. CTF application and restoration

All scalar, diagnostic, alignment, 2D, and 3D CTF evaluation uses the same
signed phase-aware transfer function. Specialized forms derive from that
single value:

| Input/operation | Numerator or data factor | Denominator/density factor |
|---|---:|---:|
| raw particle restoration | `H * Y` | `H^2` |
| already phase-flipped particle restoration | `abs(H) * Y_flipped` | `H^2` |
| CTF disabled | `Y` | `1` |
| phase-flip operation | `sign(H)` | not applicable |

This contract applies to class-average restoration and 3D Fourier-plane
generation. Noise weighting and ML regularization retain their existing
ownership outside this identity.

There is no phase-plate-specific Wiener filter, reconstructor mode,
postprocessing branch, or FSC correction. `volassemble`, FSC calculation,
automasking, and nonuniform filtering remain device-agnostic.

## 7. Fitting and UI policy

Programs that perform CTF fitting must expose `fit_phshift` and the phase search
range in their UI/CLI contract. The shared CTF-estimation iterator owns policy
resolution so standalone, batch, streaming, and validation workflows cannot
silently diverge.

The active fitting surfaces are:

- standalone `ctf_estimate`
- batch `preprocess`
- streaming `preproc`
- `mini_stream`
- `check_refpick`

Programs that only consume stored CTF parameters do not need a fitting flag.
They use the stored phase automatically.

### 7.1 Possible future evolution: shared-phase refinement after patch fitting

The current patch-based estimator first fits one micrograph-level CTF and then
holds its phase shift and astigmatism angle fixed while fitting local `dfx` and
`dfy`. The resulting patch defocus values define the spatial defocus
polynomials. This remains the supported behavior.

A possible future evolution is a final refinement of the single shared
exposure phase after the spatial defocus surface has been determined. This may
be useful for strongly tilted or curved specimens, where averaging spectra
with different local defoci can smear CTF rings and bias the phase obtained by
the initial micrograph-level fit.

The candidate staged algorithm is:

1. fit the micrograph-level defocus, astigmatism, and optional phase as now
2. fit local patch defoci while holding the shared phase fixed
3. fit the spatial `dfx` and `dfy` polynomials
4. hold the defocus surface fixed and optimize exactly one phase against an
   aggregate objective over the valid patch spectra, evaluating each patch at
   its local defocus

One carefully bounded alternation between the defocus-surface and shared-phase
steps may be considered if validation demonstrates a material benefit. An
unbounded joint alternation is not part of this evolution because phase and
defocus are correlated CTF parameters.

Any implementation of this evolution must satisfy all of the following:

- the optimized phase remains one exposure-level value; independent per-patch
  phases and spatial phase polynomials are prohibited
- no patch-phase field is added to the project schema; the final shared phase
  is stored and propagated through the existing `phshift` contract
- phase units, canonicalization, search bounds, and downstream CTF use remain
  unchanged
- patch selection, weighting, and outlier rejection are explicit parts of the
  aggregate objective
- the defocus surface is held fixed during shared-phase optimization unless a
  separately reviewed, identifiable joint model is introduced
- any new fitting control is a typed parameter exposed consistently on all CTF
  fitting surfaces; the current staged behavior remains the default until the
  new method is scientifically validated
- optimized CTF routines continue to receive a required numerical phase and
  do not gain optional phase arguments or per-pixel policy branches

Before promotion to supported policy, validation must demonstrate improved or
equivalent phase recovery on synthetic tilted/curved micrographs without
degrading local defocus recovery, stable behavior on paired laser-on and
laser-off data, no boundary pile-up or conventional-data regression, and
equivalence to the current estimator when the defocus surface is constant.
Results must also be reproducible across thread counts and in-memory versus
distributed execution where the fitting workflow supports both.

## 8. Architectural ownership

| Concern | Owner |
|---|---|
| phase representation and canonicalization | `simple_type_defs`, `ori`, and `oris` |
| CTF equation and fitting | `src/main/ctf` and the low-level CTF math kernel |
| required phase-bearing image operations | `src/main/image` and `src/main/pftc` |
| project schema, mapping, and serialization | `src/main/project` |
| STAR unit conversion and interoperability | `src/main/star` |
| fitting controls and typed parameters | `src/main/ui` and `src/main/params` |
| distributed propagation | the corresponding execution strategy |
| 2D and 3D restoration orchestration | existing matcher, class-averager, reconstructor, and `volassemble` owners |

No module-global phase-plate state may be introduced. Workflow policy belongs
in typed parameters and existing orchestration objects; numerical phase travels
through explicit arguments and domain records.

## 9. Completed code-level validation

The implementation was built and the following focused checks passed on
2026-07-20:

- full CMake build
- `simple_test_ctf`: 12/12 scalar/optimized, canonicalization, and restoration checks
- `simple_test_sp_project`: phase presence and all four project-field mappings
- `simple_test_phshift_star`: RELION degree/radian contract
- `simple_test_phshift_policy`: all five fitting-program UI contracts
- `simple_test_project_merge`: project assembly and validation behavior
- `simple_test_projdir_accumulator`: unchanged reconstruction accumulation contract

The broad pre-existing `simple_test_starfile` executable is not currently a
clean acceptance gate because of unrelated legacy failures and nondeterministic
lifetime errors. Phase-specific STAR checks must remain independently runnable
until that harness is repaired.

## 10. Scientific validation procedure

Scientific validation is required before declaring phase-shift processing
production-validated. Perform the work on a fixed SIMPLE commit and record it
as a reproducible validation report.

### 10.1 Record the validation environment

Before processing data, record:

- SIMPLE commit, build type, compiler, compiler version, FFT library, and host
- exact command lines or exported workflow configurations
- random seeds where exposed
- input accession, file list, and checksums
- RELION and cisTEM revisions used for comparisons
- whether each particle stack is raw or already phase-flipped

Build SIMPLE and run the focused tests listed in section 9. Save their complete
logs with the validation report. A code-level test failure stops the campaign.
Run file-producing tests from separate temporary directories. The expected
executables are:

```text
build/production/simple_test_ctf
build/production/simple_test_sp_project
build/production/simple_test_phshift_star
build/production/simple_test_phshift_policy
build/production/simple_test_project_merge
build/production/simple_test_projdir_accumulator
```

Before inspecting experimental results, record the quantitative tolerances to
be used for synthetic phase, defocus, and astigmatism recovery and the criteria
for a material 2D/3D regression. Use established SIMPLE tolerances where they
exist. Report results without inventing a favorable threshold after the fact.

### 10.2 Synthetic CTF-fitting validation

Generate or select synthetic spectra with known microscope parameters,
defocus, astigmatism angle, and phase. Include at least:

- phases of 0, 45, 90, 135, and a value near 180 degrees
- isotropic and non-axis-aligned astigmatic cases
- several defocus values spanning the intended experimental range
- noiseless, moderate-noise, and difficult low-SNR cases

Fit every case with `fit_phshift=yes`. Fit the phase-zero controls once more
with `fit_phshift=no`.

Report phase error, defocus error, astigmatism-angle error, fit score, and
whether any parameter reaches a search bound. Inspect the full four-parameter
solution: an apparently correct phase is not sufficient if defocus or angle is
biased.

Acceptance requires unbiased recovery in identifiable synthetic cases, stable
behavior as noise increases, and preservation of the conventional zero-phase
control within existing numerical tolerances. Persistent bound clipping,
phase/defocus swaps, or failure of the non-axis-aligned case is blocking.

### 10.3 Initial paired-data pilot

Use one small matched pair before processing the full study. The recommended
pilot is:

| Specimen | Laser phase plate on | Laser phase plate off |
|---|---:|---:|
| Aldolase A1 | EMPIAR-13528 | EMPIAR-13527 |

Use all movies if practical or choose a documented subset with the same subset
rule in both accessions. Keep motion correction, dose weighting, CTF search
limits other than phase, picking, extraction, classification, and refinement
settings identical between paired comparisons.

For each accession, run two CTF fits:

1. baseline: `fit_phshift=no`
2. phase-aware: `fit_phshift=yes phshift_min=0 phshift_max=180 phshift_step=10`

Save, per micrograph:

- fitted phase in degrees
- defocus X/Y and astigmatism angle
- CTF fit score and estimated resolution
- whether any fitted value is at or near a bound
- a representative power-spectrum/CTF overlay

Plot phase histograms and phase versus acquisition order, defocus, fit score,
and estimated CTF resolution. These plots should reveal drift,
discontinuities, and phase/defocus coupling that a single summary statistic
would hide.

The laser-on data should show a stable, interpretable nonzero phase population
and credible overlays. The laser-off control is not required to fit exactly
zero, but phase fitting must not create strong boundary pile-up, unstable
exposure-to-exposure jumps, or implausible improvements obtained by damaging
defocus/astigmatism estimates.

### 10.4 Metadata and interchange audit on real data

Choose at least five micrographs spanning the fitted phase distribution. For
each, trace one or more particles through:

```text
os_mic -> os_stk -> os_ptcl2D -> os_ptcl3D -> get_ctfparams
```

Confirm that:

- every relevant row contains `phshift`
- values agree modulo `pi` at every mapping boundary
- a true fitted zero remains present as zero
- project write/read preserves the value
- subset, merge, and distributed project assembly preserve the value
- RELION STAR export converts radians to degrees
- STAR re-import recovers the original radians within serialization tolerance
- a STAR file without `rlnPhaseShift` imports explicit zero

Include one deliberately mixed project containing conventional and nonzero
phase particles. Confirm that particle values remain authoritative after any
consolidated-stack operation.

### 10.5 2D validation

Using the pilot particles, run otherwise identical 2D processing with the
stored fitted phases and with a validation copy in which phases are forced to
zero. Validate both raw-CTF and already-phase-flipped input contracts when both
are part of supported workflows.

Save representative picked-particle views, matched class-average galleries,
class populations, rejection statistics, available class-resolution/FRC
diagnostics, command lines, seeds, and project files.

Also compare in-memory and distributed class-average restoration on the same
small particle set. Their restored averages and CTF-squared accumulators should
agree within the normal floating-point tolerance for those execution modes.

Acceptance requires internally consistent class restoration, no contrast-sign
inversion, no phase-dependent crash or NaN, and no loss of phase metadata. For
laser-on data, forcibly zeroing a credible fitted phase should not produce a
systematically better CTF-consistent result. The laser-off control must not
show a material regression when phase fitting is enabled.

### 10.6 3D validation

Reconstruct and refine the same pilot subset with otherwise identical settings:

1. stored fitted phases
2. phases forced to zero

Exercise the applicable raw and already-phase-flipped reconstruction modes. At
least one small comparison must also cover in-memory versus distributed
assembly, even/odd half maps, and any routinely used trailing or ML-regularized
reconstruction path.

Save half maps, merged maps, masks, FSC curves, reported resolutions, map
handedness/contrast sign, and reconstruction logs. Inspect representative
density features rather than relying only on one FSC crossing.

Acceptance requires phase-aware numerator and density behavior consistent with
`H*Y/H^2` or `abs(H)*Y_flipped/H^2`, agreement between execution modes within
their established tolerance, no unexplained map inversion, and unchanged
device-agnostic postprocessing. The conventional zero-phase control must remain
consistent with the established SIMPLE result.

### 10.7 External implementation cross-check

For a small common subset, compare SIMPLE against the pinned local RELION and
cisTEM revisions. Use identical physical CTF parameters and explicit unit
conversion rather than treating serialized fields as unquestioned golden
values.

Compare calculated one-dimensional and two-dimensional CTF curves, fitted
phase/defocus/astigmatism, representative overlays, and one restored 2D average
or small 3D reconstruction where practical. Document numerical sign
conventions. Resolve discrepancies at the transfer-function level before
comparing final maps.

### 10.8 Full paired-data validation

After the A1 pilot passes, process the remaining paired data with a frozen
workflow:

| Pair | Laser phase plate on | Laser phase plate off |
|---|---:|---:|
| Aldolase A1 | EMPIAR-13528 | EMPIAR-13527 |
| Aldolase A2 | EMPIAR-13529 | EMPIAR-13526 |
| Aldolase A3 | EMPIAR-13530 | EMPIAR-13525 |
| Hemoglobin H1 | EMPIAR-13535 | EMPIAR-13533 |
| Hemoglobin H2 | EMPIAR-13534 | EMPIAR-13532 |
| Hemoglobin H3 | EMPIAR-13537 | EMPIAR-13531 |

For every accession, retain the CTF distributions, representative overlays,
particle/class diagnostics, maps, and FSC curves. The paired study is a
scientific validation signal, not a requirement that every small-sample metric
for an on dataset exceed its off partner. Investigate consistency, failure
modes, and whether the software uses the fitted phase correctly.

### 10.9 Performance check

Benchmark fused alignment preparation, class-average restoration, and 3D
Fourier-plane generation/reconstruction on the same hardware and inputs. Use
repeated runs and report medians. Compare a zero-phase baseline with the
phase-aware implementation and compare zero with nonzero phase in the same
implementation.

The optimized kernel should add only one scalar phase and no per-pixel policy
branch. A repeatable slowdown greater than approximately 5% should be
investigated before production sign-off.

### 10.10 Validation deliverables and sign-off

The validation package should contain:

- a manifest with commits, environment, data checksums, and command lines
- focused test and build logs
- synthetic truth and recovery tables
- per-micrograph CTF tables and phase-distribution plots
- selected spectrum/CTF overlays
- project-field and STAR round-trip audit tables
- matched 2D class galleries and diagnostics
- 3D half maps, merged maps, FSC curves, and processing logs
- performance measurements
- a deviations/failures log, including investigated negative results

Production scientific sign-off requires all of the following:

1. focused code tests pass on the validation build
2. synthetic fitting recovers identifiable phases without systematic coupling
3. real laser-on fits are stable and interpretable without widespread bound clipping
4. conventional controls show no unexplained regression
5. project and STAR round-trips preserve numerical phase and units
6. 2D and 3D results demonstrate that stored phase is consumed
7. distributed and in-memory execution agree within established tolerances
8. no unexplained contrast inversion, NaN, crash, or device-specific postprocessing appears
9. performance is acceptable or any regression is understood and approved

Record any failed item with the smallest reproducible project and exact command
line before modifying fit limits or other scientific settings. Do not tune the
acceptance test until a software defect has been excluded.

## 11. External basis

The policy is consistent with the mature phase-shift treatment audited in:

- RELION `f6d3ad20cde0acc32469c85bf0bcafd2be84c943`
  (`5.0.1-11-gf6d3ad20`)
- cisTEM `455de9a1e84bf96fc572ca94d9657aed25cd7a5d`
  (`1.0.0-979-g455de9a1`)

Both always include a stored numerical phase in CTF evaluation and use a
separate policy to enable fitting. RELION STAR stores phase in degrees; cisTEM
and SIMPLE use radians internally.

The scientific validation data are described by Petrov et al.,
["Laser phase plate improves structure determination of small proteins by
cryo-EM"](https://doi.org/10.1126/science.aeh0665).

## 12. Non-goals and prohibited regressions

Outside this policy unless separately reviewed:

- restoring `phaseplate_correct_fsc` or flattening low-frequency FSC values
- restoring deleted partial-Wiener or first-zero experimental modes
- adding Volta-versus-laser reconstruction branches
- inferring acquisition hardware from a numerical phase
- estimating independent per-patch phases or fitting a spatial phase surface
- estimating independent per-particle phases without supporting evidence
- moving restoration or postprocessing ownership into commanders

Future changes must not:

- reintroduce `phaseplate` or `l_phaseplate` as processing switches
- make phase arguments optional in optimized CTF routines
- omit `phshift` from `ctfparams`, `ctfvars`, or relevant project records
- condition downstream phase use on `fit_phshift`
- perform STAR degree/radian conversion inside numerical kernels
- replace heterogeneous particle phases with an aggregate stack default

## 13. Change-control checklist

Any change touching phase fitting, CTF evaluation, project mappings, STAR
translation, class restoration, or reconstruction must verify:

- zero-phase scalar and optimized results remain within current tolerance
- nonzero phase reaches every affected numerical consumer
- `ctfparams` and `ctfvars` remain consistent
- all four project fields retain phase through write/read and mapping
- RELION STAR conversion remains degrees <-> radians
- all five fitting UI surfaces remain synchronized
- raw and phase-flipped restoration retain their numerator/density contracts
- no phase-plate-specific postprocessing branch has been introduced

Update this policy and the nearest focused regression test whenever any of
these contracts changes.
