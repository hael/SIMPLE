# Microscope-Dependent Parameter Policy

This document records the current SIMPLE policy for microscope-dependent and
optics/CTF-model parameters, especially acceleration voltage (`kv`), spherical
aberration constant (`cs`), and amplitude contrast fraction (`fraca`). It
describes the behavior in the current codebase, not the desired future state.

## Scope

The microscope-parameter model used for CTF evaluation is:

- `kv`: acceleration voltage, in kV
- `cs`: spherical aberration constant, in mm
- `fraca`: amplitude contrast fraction
- `smpd`: sampling distance, in Angstroms per pixel

`fraca` is not technically a microscope hardware parameter, but it is part of
the same model identity for SIMPLE because it changes the CTF that is fitted and
applied. It must therefore travel with `kv`, `cs`, and `smpd` when datasets are
combined.

Defocus parameters (`dfx`, `dfy`, `angast`, `phshift`) are image- or
particle-dependent CTF parameters. They are closely related, but they are not
part of the microscope/CTF-model identity considered here. Numerical phase is
governed by the [phase-shift CTF policy](phase_shift_ctf_policy.md); acquisition
hardware is not a SIMPLE CTF execution flag.

## Current Data Model

The canonical in-memory CTF bundle is `ctfparams` in
`src/defs/simple_type_defs.f90`. It carries `smpd`, `kv`, `cs`, `fraca`,
defocus, phase shift, CTF flag, and phase-fitting policy.

Project metadata stores these values in several places:

- `os_mic`: movie/micrograph rows. CTF fitting reads CTF-model parameters
  from these rows.
- `os_stk`: particle-stack rows. 2D and 3D analysis read CTF-model
  parameters from these rows through `sp_project%get_ctfparams`.
- `os_ptcl2D` / `os_ptcl3D`: particle rows. These normally carry particle
  CTF values such as defocus, astigmatism angle, phase shift, and `stkind`.
- `os_optics`: optics-group rows. These can carry `ogid`, `ogname`, `smpd`,
  `kv`, `cs`, `fraca`, population, centroid, and box information.

The important current rule is that `os_optics` is not the analysis-time source
of truth. It is used for STAR import/export, stream optics assignment, plotting,
and group propagation, but CTF evaluation in 2D/3D does not resolve the
CTF-model parameters by `ogid`.

## Import Policy

### Movie and micrograph import

`import_movies` accepts scalar command-line values for `smpd`, `kv`, `cs`,
`fraca`, and `ctf`.

When importing into a project that already contains movies or micrographs,
`exec_import_movies` and `sp_project%add_movies` compare the new values against
the first existing micrograph and abort if `smpd`, `kv`, `cs`, `fraca`, CTF
status differ.

When importing integrated micrographs with `deftab`, per-micrograph defocus can
come from the table, but `smpd`, `kv`, `cs`, and `fraca` still come from the
single command-line CTF bundle and are written onto every micrograph row.

Therefore direct mixed-microscope movie or micrograph import into one project is
not currently supported by `import_movies`.

### Particle import

`import_particles` has several paths with different behavior:

- STAR import reads the STAR optics table, imports optics rows, and appends
  optics values to imported micrograph, stack, or particle rows. Command-line
  `smpd`, `kv`, `cs`, and `fraca` are intentionally ignored in this path.
- Single-stack import and per-particle `stktab` import use one scalar CTF
  bundle for all stack rows.
- Per-stack `stktab` import can preserve per-stack CTF-model parameters when
  the input metadata contains one row per stack. If `kv`, `cs`, or `fraca` are
  also provided on the command line, those command-line values are applied to
  all rows.

`import_cavgs` only imports class averages and sampling distance. It does not
model the CTF-model parameter set.

## CTF Fitting

CTF fitting is micrograph-driven. The CTF estimation strategy reads each active
micrograph row and calls `o%get_ctfvars()` on that row before fitting. The
fitting code then constructs the CTF model from that row's `smpd`, `kv`, `cs`,
and `fraca`.

This means CTF fitting can use per-micrograph CTF-model parameters if they are
already present on `os_mic`. The primary SIMPLE movie/micrograph importer,
however, prevents creating such mixed rows in a single project.

## Extraction

Particle extraction reads CTF values from each micrograph row and passes them to
`sp_project%add_stk`, which writes `smpd`, `kv`, `cs`, `fraca`, CTF status, and
the numerical phase onto the generated stack and particle rows. It also
propagates `ogid` from micrographs to particles when present.

The current extractor resets `ctfparms%smpd` to the run-level `params%smpd`
before adding the stack. That is consistent with the current same-sampling
project model and with workflows where projects have already been scaled to a
common sampling distance.

## 2D and 3D Analysis

2D and 3D CTF-aware analysis call `sp_project%get_ctfparams`. That function maps
particles to `stkind`, reads `smpd`, `kv`, `cs`, `fraca`, and CTF status from
`os_stk`, and reads defocus-related values, including the authoritative phase,
from the particle row. A stack-level query reads the stack phase instead. It
does not consult `os_optics`.

Consequences:

- 2D/3D analysis can handle different `kv` and `cs` values across stacks if the
  correct values are present on each `os_stk` row.
- 2D/3D analysis will not pick up a corrected value that exists only in
  `os_optics`.
- `ogid` is metadata for these paths unless duplicated stack-row values are also
  correct.

Project-level sampling helpers such as `sp_project%get_smpd` and parameter
derivation read the first stack or micrograph. The current workflow therefore
assumes one effective analysis sampling distance, even if `kv`, `cs`, or
`fraca` vary by stack.

## Optics Groups

SIMPLE has an `os_optics` segment and STAR-compatible optics-table support.
`assign_optics_groups` creates groups from filename-derived beam-tilt
information, optional XML beam shifts, `maxpop`, and `optics_offset`. It then
propagates `ogid` to micrographs, stacks, and particles and writes optics rows.

Current optics-group assignment does not split groups by `kv`, `cs`, `fraca`,
or `smpd`. The optics row inherits CTF-model values from a
representative image in the group. This is adequate for same-microscope
beam-tilt grouping, but it is not a safe policy for mixed microscopes or mixed
CTF-model settings unless the input grouping already prevents those settings
from being mixed.

## Project Merging

The refactored `merge_projects` path accepts a `projtab` file listing N source
projects and writes an explicit `projfile_merged` output. It does not use the
legacy in-place `projfile` plus `projfile_target` append interface.

The merger is project-field based: every input must have the same populated
data segments, but the row counts may differ. This permits stack-only,
movie/micrograph-only, particle, class/output-bearing, or optics-bearing
projects to merge when their populated fields match.

For already scaled particle projects with the same particle sampling distance
and correct per-stack CTF-model values, merging preserves the information that
2D/3D analysis currently uses. The merger does not reject projects merely
because `kv`, `cs`, or `fraca` differ across row-level CTF-model groups.

Row-level `ogid` assignments are remapped during merge even when no `os_optics`
segment is present. `os_optics` remains optional metadata and is not used to
backfill missing authoritative CTF-model values.

## Practical Current Guidance

The supported direct `import_movies` path is homogeneous in `smpd`, `kv`, `cs`,
`fraca`, and CTF flag.

For multiple microscopes that have already been scaled to the same sampling
distance, the safest current workflow is:

1. Import and process each microscope or CTF-model parameter group in its own
   project with the correct scalar `kv`, `cs`, `fraca`, and `smpd`.
2. Run CTF fitting and extraction in those projects so the generated stack rows
   carry the correct CTF-model values.
3. Merge the projects with `merge_projects projtab=<table> projfile_merged=<out>`
   after the particle stacks are at a common sampling distance and compatible
   box policy.
4. Before 2D or 3D analysis, verify that each relevant `os_stk` row contains the
   intended `smpd`, `kv`, `cs`, `fraca`, and `ctf` values, and verify that the
   particle rows carry the intended numerical `phshift`.

Do not rely on `os_optics` alone for CTF-aware 2D/3D analysis in the current
codebase.

## Current Gaps

- There is no single resolver/validator for the row-level CTF-model parameter
  set used by CTF fitting and 2D/3D analysis.
- `import_movies` has no per-file CTF-model parameter table.
- `assign_optics_groups` can group images with different CTF-model settings.
- Several consumers duplicate CTF-model parameters onto `mic` and `stk` rows.
- Analysis-time CTF lookup ignores `ogid` and `os_optics`.
- Sampling-distance defaults are project-global in practice because they are
  derived from the first stack or micrograph.
