# Project micrograph subset selection

## Contract

`simple_exec prg=select_micrographs projfile=input.simple nmics=100`
explicitly opts into a new, independent project operation. Other workflows
never call it. It writes a separate `_mic_subset.simple` sibling by default,
keeps original row/image identifiers, poses and half assignments, and never
modifies the input project or image stacks. Existing rejected rows stay rejected.
`outfile=` may name a new project elsewhere; image paths are made absolute.
The command uses `mkdir=no` to avoid execution-directory project copying.
Outputs also include per-micrograph and per-direction CSV reports.

`subset_mode=quality_coverage` uses active, updated particles in `state=1`
(or the requested state), with finite poses and positive finite `subset_score`
(`corr`, default; `corr_cart`; `none` for unweighted support). These scores
are ranking proxies, not calibrated signal fractions or frame-FRC measurements.
Symmetry-aware spiral directions (`nspace=64`, even; `pgrp=c1` by default) are
tracked separately for existing particle halves. Complete project half
provenance is required; selection never repartitions particles.

The objective is the sum of `log(1 + support/scale)` across bins, with scale
equal to the expected selected total support divided by all bins. Unlike a
hard cap, additional good particles keep contributing. Empty directions stay
in the report. This is an angular-support proxy, not a proven volume optimum.
The scale is computed once from the complete candidate pool and requested
cardinality, never fitted to a reconstruction outcome.

`subset_mode=particles` exactly selects the micrographs with the largest
currently active particle counts in the requested state, without requiring
orientations, half provenance or scores. Ties use original micrograph order.

## Implementation and review

- A thin project commander owns parameter construction, file checks, reports
  and output project identity. Numerical selection and project masks belong
  in `simple_mic_subset` next to `sp_project`.
- Quality/coverage uses deterministic greedy initialization followed by global
  branch-and-bound. `subset_max_nodes=100000` bounds search work. A complete
  search proves the objective optimum; a limited search reports an incumbent
  and a conservative upper bound, explicitly marked unproven.
- Explicit `prg=select_micrographs` is the default-off feature gate: existing
  project commands and refinement paths are unchanged.
- Validate the project-owned stack/micrograph mapping; reject absent or
  inconsistent provenance instead of inventing micrograph identities.
- Output deselects noneligible particles and unchosen mic/stack rows without
  pruning. Clear inherited output products/state weights and sigma state,
  because they describe the parent selection.
- Prior synthetic experiments showed quality/coverage can trade overall
  correlation against directional error. The command reports its proxy and
  never labels the selected subset a validated better reconstruction.

## Validation

Register a project-area unit sub-suite. Cover a closed-form greedy trap,
comparison with an independent exhaustive oracle, search-limit reporting,
count-mode optimum, invalid budgets, masks preserving inactive rows, poses,
structural labels, physical image indices and halves, and invalid mappings.
Run static registry/DAG/header/whitespace checks and language-server syntax
diagnostics. Compilation and runtime tests remain for the user for this
implementation task; earlier synthetic-test authorization is not a build request.

## Current status

Implemented (commit `86bf97213`): the `select_micrographs` command, its UI
entry and router entries, `simple_mic_subset`, and the tester
`simple_mic_subset_tester`. An early build aborted in `bind_input` while
registering the command, because the shared `pgrp` descriptor is required and
conflicts with Advanced visibility. The command now declares its own optional
symmetry input with the default `c1`. Rebuilding and runtime validation of the
latest tree remain with the maintainer.

The user's first build exposed `trim(prg)` on a derived `string`. Both the
parameter-validation guard and commander guard now use the existing overloaded
string comparisons directly. Earlier language-server checks did not detect
this type error; compiler validation remains pending.

The selector was restored from the user-work stash after the warning-only
test-budget fix. The stash is retained; profiling output and the older
stashed budget-default edit were not reapplied. The user confirmed that the
current interface should take an exact input number of micrographs (`nmics`);
automatic minimum-count estimation and lower-bound options are deferred.

Implemented commander, domain functions, UI/router/typed parameters, reports,
user guide and the `unit_project suite=micrograph_subset` sub-suite.
Registry, FLEX DAG, source headers and whitespace checks pass. New Fortran
files have clean fortls diagnostics; the parameters phase file has three
pre-existing local-shadow warnings away from this change. An independent
Python calculation confirms the greedy-trap fixture's closed form.
No compilation or Fortran test execution has been performed. Required next
validation is the user's build and `simple_test_exec test=unit_project
suite=micrograph_subset`, followed by representative real-project selections
and reconstructions. Navigation artifacts were regenerated.
