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
- [Setting up a machine for development with Claude Code](../how2s/how2set_up_claude_code_dev.md)
- [Writing a plan for an agent](../how2s/how2write_agent_plans.md): what the
  record of the runs so far says a plan needs
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

Current navigation (2026-11-17) separates active plans from retained historical
limits. Removed follow-up notes are not recreated, and their absence does not
establish a passing result:

1. [Batch scheduler completion record](../implementation_notes/completed/qsys_scheduler_completion_bug_report_completed.md):
   the batch exit-code and lost-part detection landed (`2307ea96f`). The completed
   record retains the periodic-report/wall-clock-backstop limitation; its former
   planned companion is no longer present. Recheck current code before reopening it.
2. [Testing and benchmarking](../implementation_notes/planned/testing_and_benchmarking_plan.md#6-consolidated-test-and-ci-work):
   the single live list of remaining SINGLE, stream, workflow-truth and CI/evidence
   work, consolidated on 2026-11-17. The three handovers are archived as superseded;
   their incomplete validation is not a passing result. C-alpha belongs to Ruben.
3. [Standing nightly triage](../refactoring_notes/planned/nightly_suite_triage.md):
   a plan that can be re-executed. Each run regenerates its report and updates the
   consolidated testing plan; the old report is no longer in the tree.
4. [Random seed parameter](../refactoring_notes/planned/rnd_seed_parameter.md):
   approved on 2026-10-18, with the design refreshed for 64-bit keyed RNG and
   restart metadata. Implementation remains open.
5. [Policy gap register](../refactoring_notes/planned/policy_set_audit_2026-10-21.md):
   unresolved decisions and candidates, with current checks separated from October
   claims that need revalidation and from closed history.

[Multi-microscope metadata](../refactoring_notes/planned/need2basis/multi_microscope_parameters.md)
and MRC suffix naming remain deferred in `need2basis/`; the two architecture
proposals remain in `next_release/`. [FLEX tree work](../refactoring_notes/planned/flex_tree_strategy.md)
belongs to Afan. The former FLEX follow-up notes have been removed; the completed
records retain their historical validation limits without claiming those gates passed.

Splits: where a note mixed finished and open work it was split on 2026-10-20 into
`<name>_completed.md` (under `completed/`) and `<name>_unfinished.md` (under `planned/`).

Additional feature proposals, active experiments, reusable plans and remaining
validation are indexed by topic in the two READMEs. Review the owning record
before starting work so an existing design is extended rather than duplicated.
