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

- [Fortran module index](../code_overview/fortran-indexes/module_index.md)
- [Architecture design guide](guides/architecture_design_dev_guide.docx)
- [Test environment policy](../policies/test_environment_policy.md)
- [Memory estimator and benchmark guide](../how2s/memory_estimator.md)
- [Developer onboarding](onboarding/)

## Work that still needs attention

The source-backed audit on 2026-10-05 found these concrete open areas:

1. [Batch scheduler completion detection](../implementation_notes/planned/qsys_scheduler_completion_bug_report.md):
   streaming jobs have exit-code files, but batch `generate_script_1` jobs can
   still wait forever when a worker exits before writing the success sentinel.
2. [Stream refactor follow-up](../refactoring_notes/planned/stream_refactor.md)
   and [stream test handover](../refactoring_notes/planned/stream_area_tests_handover.md):
   remaining cleanup and workflow-level coverage are explicitly unfinished.
3. [Workflow truth gates and nightly runner](../refactoring_notes/planned/phase5_workflow_gates_and_nightly_runner_handover.md):
   the simulation-truth gates and unattended nightly evidence are still open.
4. [Probabilistic-alignment peak-memory refactor](../refactoring_notes/planned/probabilistic_alignment_table_peak_memory_refactoring.md):
   implementation and equivalence/memory validation are incomplete.
5. [Multi-microscope parameter follow-up](../refactoring_notes/planned/multi_microscope_parameters.md):
   the first phase landed, but resolver work and open policy questions remain.
6. [Release 4 cleanup inventory](../refactoring_notes/planned/release4_legacy_cleanup_inventory.md)
   and [continuous-pose aftermath](../refactoring_notes/planned/pose_cont_refactor_aftermath_issues.md):
   these are the current cleanup ledgers.
7. [FLEX PCA architecture record](../refactoring_notes/completed/flex_pca_architecture_audit_and_refactoring_plan_2026_09_18.md):
   phases 0b--11 landed, but the phase-0a cluster baseline/replay is still an
   outstanding validation activity.

The remaining files under the two `planned` directories are feature proposals
or lower-priority optimizations. They should be reviewed there before starting
new work so an existing design is extended rather than duplicated.
