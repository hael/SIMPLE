# Refactoring-note status

Every refactoring note belongs to one lifecycle directory:

- [completed](completed/) — landed refactors and completed reviews
- [planned](planned/) — open refactors, validation, and test handovers
- [rejected](rejected/) — retired or rejected architectural directions

The directory is authoritative. A completed review may identify follow-up work;
that work belongs in a separate planned record rather than changing the review's
status.

The topic index below identifies the live records. Cross-cutting feature work is
indexed in [implementation notes](../implementation_notes/README.md).

Record conventions:

- Existing `<name>_completed.md` / `<name>_unfinished.md` splits preserve the
  delivered record and one live list of open work. New development keeps its
  contract, plan, review findings and validation criteria in one living note.
- [Nightly triage](planned/nightly_suite_triage.md) deliberately keeps its
  standing plan and replaceable run report together; retain that arrangement.
- `planned/next_release/` holds proposals deferred past the current release, and
  `planned/need2basis/` holds proposals that need a decision or basis before work can start;
  the `rejected/` folder is top level only (no nested `planned/rejected/`).

For the cross-cutting audit summary and the most important open work, see the
[developer documentation index](../for_developers/index.md).

## Open work by topic

Navigation reviewed 2026-11-02. Each note owns its detailed acceptance criteria;
this index does not establish a new passing test result.

| Topic | Record | Status and next action |
| --- | --- | --- |
| Nightly suite | [Standing triage](planned/nightly_suite_triage.md), [latest report](planned/nightly_suite_triage_report.md) | Reconcile the latest run's findings and maintainer rulings; preserve the recurring workflow. |
| Workflow truth/CI | [Phase 5 closeout](planned/phase5_workflow_gates_and_nightly_runner_handover.md) | Implementation delivered; remaining truth checks and complete unattended evidence need closeout. |
| SINGLE coverage | [Test handover](planned/single_area_tests_handover.md) | Most checks delivered; finish the remaining assertion gaps. |
| Stream coverage | [Test handover](planned/stream_area_tests_handover.md) | Reconcile remaining workflow coverage and observed runs. |
| Stream architecture | [Follow-up](planned/stream_refactor.md) | Stage switch-over complete; remaining ownership/helper cleanup and workflow evidence remain. |
| Stream snapshots | [3D snapshots](planned/stream_3D_snapshots_plan_2026-10-06.md) | Reconcile implementation and end-to-end/UI acceptance against the plan. |
| FLEX reuse | [SIMPLE machinery closeout](planned/flex_on_simple_machinery_unfinished.md) | Phases implemented; prescribed equivalence gates and follow-ups remain. |
| FLEX architecture | [PCA architecture closeout](planned/flex_pca_architecture_audit_and_refactoring_plan_2026_09_18_unfinished.md) | Refactor implemented; baseline capture/replay and optional cleanup remain. |
| Alignment-table memory | [Validation closeout](planned/probabilistic_alignment_table_peak_memory_refactoring_unfinished.md) | Reconcile scientific equivalence and multi-state memory measurements with later assignment-policy changes. |
| Continuous pose | [Aftermath validation V1](planned/pose_cont_refactor_aftermath_issues_unfinished.md) | Reconcile requested gates against their revisions before closing the implemented changes. |
| Local refinement | [Workflow record](planned/refine3d_local_workflow.md) | Preserve the recorded integration and equivalence limitations; the checkout-local evidence artifact is unavailable here. |
| Release cleanup | [Outstanding inventory](planned/release4_legacy_cleanup_inventory.md) | Recheck the dated open/confirm/deferred items before acting. |
| Policy debt | [Gap register](planned/policy_set_audit_2026-10-21.md), [dated supporting audits](planned/policy_set_audit_2026-10-21/) | The register owns dispositions; raw cluster reports preserve the pre-correction evidence. |
| Random seed | [Parameter proposal](planned/rnd_seed_parameter.md) | Awaiting review; production still uses the environment seed. |
| Tomography | [refine_tomo](planned/refine_tomo.md) | Reconcile compilation and ground-truth evidence with subsequent workflow work. |
| Microscope metadata | [Multi-microscope parameters](planned/need2basis/multi_microscope_parameters.md) | First phase delivered; audit remaining CTF construction/merge tests and policy decisions. |
| MRC naming | [Stack/volume suffixes](planned/need2basis/mrc_stack_volume_output_suffix_refactoring.md) | Decision/basis needed before migration. |
| Deferred architecture | [Domain-driven design](planned/next_release/staged_domain_driven_design_refactor.md), [execution-layer adoption](planned/next_release/x_execution_layer_adoption.md) | Next-release proposals. |
