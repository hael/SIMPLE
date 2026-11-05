# Developer documentation

Developer lifecycle records are classified by directory:

- [Completed implementation notes](../implementation_notes/completed/) describe
  implemented features or finished investigations and reviews.
- [Planned implementation notes](../implementation_notes/planned/) describe
  features, fixes, and experiments that are not complete.
- [Rejected implementation notes](../implementation_notes/rejected/) preserve
  retired, superseded, or rejected designs when their rationale is still useful.
- [Completed refactoring notes](../refactoring_notes/completed/) describe landed
  structural changes or completed architecture reviews.
- [Planned refactoring notes](../refactoring_notes/planned/) describe open
  refactors, validation work, and test handovers.
- [Rejected refactoring notes](../refactoring_notes/rejected/) preserve refactors
  that should not be resumed without a new decision.

`completed` is the repository's existing name for the implemented/finished
classification. A completed design record may still mention optional research,
future optimization, or user-owned runtime validation; those do not make its
core deliverable unimplemented. Conversely, a completed review can identify
work that is tracked separately under `planned`.

Current guides and policies live outside these lifecycle folders because they
describe the software as it is, rather than the history of one change:

- [Policy index](../policies/README.md): every policy, with what to read before changing what
- [Layering and project ownership policy](../policies/layering_and_project_ownership_policy.md):
  read before adding or changing a commander, a strategy routine, a worker step, or any code that
  reads or writes a project
- [Fortran module index](../code_overview/fortran-indexes/module_index.md)
- [Architecture design guide](guides/architecture_design_dev_guide.docx)
- [Test environment policy](../policies/test_environment_policy.md)
- [Memory estimator and benchmark guide](../how2s/memory_estimator.md)
- [Developer onboarding](onboarding/)

## Work that still needs attention

Navigation reconciled 2026-11-04. The topic indexes list the active records and
their next actions; the linked notes retain the dated evidence and acceptance
criteria. This navigation update supplies no new runtime validation.

- [Implementation topics and completed investigations](../implementation_notes/README.md)
- [Refactoring topics and validation closeout](../refactoring_notes/README.md)
- [Current atomicity-prior algorithm](../algorithms/atomicity_prior.md), including
  the fitted-model mode absent from the original methods study

The delivered polar shift-search performance refactor, `refine3D_local`
workflow/ownership refactor, stream refactor, 3D snapshots, Release 4 cleanup,
unreferenced-code review, historical pose/alignment-table closeouts and dated
policy-audit reports are indexed under
[completed refactors and reviews](../refactoring_notes/README.md#completed-refactors-and-reviews).
Their recorded evidence limits are retained. The policy audit's live gap
register remains in `planned/`.

Completed implementation records and reports, including trajectory denoising,
posterior search, atomic density/B factors, resource planning, truth-suite
tooling, BLAS configuration and the nano-context investigation, are indexed
under [completed implementation evidence](../implementation_notes/README.md#completed-investigations-and-reference-material).
Their validation and research limits remain explicit; later proposals and the
testing/benchmarking backlog retain their open scope.

The earlier source-backed audit identified the following open areas. Missing
split records have been reconstructed from the surviving completed records;
their absence was not evidence that the outstanding work passed:

1. [Batch scheduler completion detection](../implementation_notes/planned/qsys_scheduler_completion_bug_report_unfinished.md):
   the batch exit-code and lost-part detection landed (`2307ea96f`); the periodic
   report of outstanding parts, the wall-clock timeout and the unused
   `qsys_watcher_diag` remain. Background:
   [the completed report](../implementation_notes/completed/qsys_scheduler_completion_bug_report_completed.md).
2. [Stream test handover](../refactoring_notes/planned/stream_area_tests_handover.md):
   remaining workflow-level coverage is explicitly unfinished. The delivered
   [stage refactor](../refactoring_notes/completed/stream_refactor.md) is closed,
   retaining its follow-ups and validation limits.
3. [Workflow truth gates and nightly runner](../refactoring_notes/planned/phase5_workflow_gates_and_nightly_runner_handover.md)
   and [standing nightly triage](../refactoring_notes/planned/nightly_suite_triage.md):
   reconcile the remaining truth checks and unattended evidence with the latest
   [triage report](../refactoring_notes/planned/nightly_suite_triage_report.md).
4. [Multi-microscope parameter follow-up](../refactoring_notes/planned/need2basis/multi_microscope_parameters.md):
   the first phase landed; the CTF-construction audit, most of the merge tests,
   and the open policy questions remain.
5. [FLEX PCA architecture validation](../refactoring_notes/planned/flex_pca_architecture_audit_and_refactoring_plan_2026_09_18_unfinished.md):
   phases 0b--11 landed
   ([record](../refactoring_notes/completed/flex_pca_architecture_audit_and_refactoring_plan_2026_09_18_completed.md));
   the phase-0a cluster baseline/replay is still an outstanding validation activity.

Splits: where a note mixed finished and open work it was split on 2026-10-20 into
`<name>_completed.md` (under `completed/`) and `<name>_unfinished.md` (under `planned/`).

Additional feature proposals, active experiments, test handovers and remaining
validation are indexed by topic in the two READMEs. Review the owning record
before starting work so an existing design is extended rather than duplicated.
