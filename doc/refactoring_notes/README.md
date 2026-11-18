# Refactoring-note status

Every refactoring note belongs to one lifecycle directory:

- [completed](completed/) — landed refactors and completed reviews
- [planned](planned/) — open refactors, validation, and reusable plans
- [rejected](rejected/) — retired or rejected architectural directions

The directory is authoritative. A completed review may identify follow-up work;
active work belongs in a separate planned record rather than changing the review's
status. When the maintainer closes the delivered scope, deferred improvements,
unverified checks and release reminders can remain in its completed record as
explicit limits, without a new implementation task or passing-test claim.

The topic index below identifies the live records. Cross-cutting feature work is
indexed in [implementation notes](../implementation_notes/README.md).

Record conventions:

- Existing `<name>_completed.md` / `<name>_unfinished.md` splits preserve the
  delivered record and one live list of open work. New development keeps its
  contract, plan, review findings and validation criteria in one living note.
- [Nightly triage](planned/nightly_suite_triage.md) is a standing plan; each run
  writes a fresh report beside it.
- `planned/next_release/` holds proposals deferred past the current release, and
  `planned/need2basis/` holds proposals that need a decision or basis before work can start;
  the `rejected/` folder is top level only (no nested `planned/rejected/`).

For the cross-cutting audit summary and the most important open work, see the
[developer documentation index](../for_developers/index.md).

## Open work by topic

Navigation reviewed 2026-11-17. Eight Markdown plans remain in `planned/`
(including its four deferred proposals). Each note owns its detailed acceptance criteria;
this index does not establish a new passing test result.

| Topic | Record | Status and next action |
| --- | --- | --- |
| Nightly suite | [Standing triage](planned/nightly_suite_triage.md) | Re-executable plan; each run records fresh findings and feeds the consolidated testing plan. |
| Tests and CI | [Testing and benchmarking plan](../implementation_notes/planned/testing_and_benchmarking_plan.md#6-consolidated-test-and-ci-work) | One live list for remaining SINGLE, stream, workflow truth and CI/evidence work. C-alpha belongs to Ruben; runner closeout to Ruben, stream work to Joe/Ruben. |
| FLEX | [Tree strategy](planned/flex_tree_strategy.md) | Afan owns this work; leave its implementation and closeout to him. |
| Policy debt | [Gap register](planned/policy_set_audit_2026-10-21.md) | Unresolved items first, closed items as dated history. Limited November checks are distinguished from October claims awaiting recheck. |
| Random seed | [Approved parameter plan](planned/rnd_seed_parameter.md) | Approved 2026-10-18; updated for 64-bit keyed RNG, seed offsets and restart metadata. Implementation remains open. |
| Microscope metadata | [Multi-microscope parameters](planned/need2basis/multi_microscope_parameters.md) | Deferred in `need2basis/`; retain the delivered first phase and remaining CTF/merge/policy questions until work is reopened. |
| MRC naming | [Stack/volume suffixes](planned/need2basis/mrc_stack_volume_output_suffix_refactoring.md) | Decision/basis needed before migration. |
| Deferred architecture | [Domain-driven design](planned/next_release/staged_domain_driven_design_refactor.md), [execution-layer adoption](planned/next_release/x_execution_layer_adoption.md) | Next-release proposals. |

## Completed refactors and reviews

These records preserve dated delivery, review and supersession decisions.
Their validation limits and follow-up measurements remain explicit; archiving
supplies no new passing result.

| Topic | Record | Delivered scope and retained limits |
| --- | --- | --- |
| Test handovers | [SINGLE](completed/single_area_tests_handover.md), [stream](completed/stream_area_tests_handover.md), [Phase 5](completed/phase5_workflow_gates_and_nightly_runner_handover.md) | Superseded 2026-11-17 by the consolidated testing plan. Original evidence and unchecked validation remain in these historical records. |
| Tomography | [refine_tomo](completed/refine_tomo.md) | Closed by maintainer decision on 2026-11-17. Implementation compiled; the dedicated phantom orientation-recovery test remains a validation limit. Historical FePt evidence is retained. |
| Bootstrap class-average applications | [Retired contract](completed/class_average_bootstrap_retired.md) | `bootstrap_cavgs` and `unbootstrap_cavgs` retired on 2026-11-17; use `cls_expansion`. Their exclusive parameters and the old project-overwrite path are removed. |
| Agent-friendly code | [Plan](completed/agent_friendly_code.md), [run report](completed/agent_friendly_code_report.md) | Folder guides and defect fixes 2026-11-04; the run `agent_friendly` (2026-11-05, uncommitted for review): the command-line gate check, the small fixes and the checks of 3.1, the medium items, renames, glossary, measurements of I1, O1, O4, S4, and (phase 6, ruled after the report) O1 and O4 adopted, X1/G3 finished, the particle-row slot stop and the sieve fix. The items that change `ms_posterior2`'s files (R4/X4, C3, C4, R6, 3.7 and the phase-6 remainder) go to the follow-up run `agent_friendly2`. |
| Polar shift-search performance | [Refactoring record](completed/polar_shift_search_performance_refactoring.md) | Integrated in `69a623cb8`, with recorded Debug/Release and simulated-workflow validation. Experimental-data timing/memory and Linux/BOX scaling remain follow-up measurements. |
| Local refinement | [Workflow and architecture record](completed/refine3d_local_workflow.md) | Workflow and startup/ownership cleanup delivered in `0c4801f61` and `9c9ec1105`. Preserve the dated integration, replay-equivalence and broader scientific-validation limits; the checkout-local evidence artifact remains unavailable here. |
| Stream architecture | [Stage refactor](completed/stream_refactor.md) | Stage switch-over complete; class-average selection tester registered. Archived by maintainer decision with helper-placement/optics follow-ups and workflow/GUI validation limits retained. |
| Stream snapshots | [3D snapshots](completed/stream_3D_snapshots_plan_2026-10-06.md) | SIMPLE and NICE implementation delivered. Archived by maintainer decision with test-run evidence and the two final-project/selection checks explicitly unverified. |
| Release cleanup | [Closeout inventory](completed/release4_legacy_cleanup_inventory.md) | Implementation work closed by maintainer decision; dated confirmations, deferred Joe/GPU items and the v4.0.0 release-tag documentation reminder remain in the record. |
| Unreferenced-code review | [Closed inventory](completed/dead_code_inventory_2026-11-02.md) | Approved XML replacement and external PDB data delivered in `ef5fe84d5`; other deletion proposals withdrawn. The owner ruling remains explicit and authorizes no further removals. |
| Continuous-pose aftermath | [Implementation record](completed/pose_cont_refactor_aftermath_issues_completed.md), [validation closeout](completed/pose_cont_refactor_aftermath_validation_closeout.md) | Fixes/tests delivered; historical V1 validation closed by maintainer decision with unrecorded checks retained as limits, without a passing-run claim. |
| Alignment-table memory | [Implementation record](completed/probabilistic_alignment_table_peak_memory_refactoring_completed.md), [validation closeout](completed/probabilistic_alignment_table_peak_memory_validation_closeout.md) | Refactor delivered; historical memory/equivalence closeout archived by maintainer decision. Missing measurements and changed-policy comparison limits remain explicit. |
| Policy-set audit evidence | [Seven dated audit reports](completed/policy_set_audit_2026-10-21/) | Finished source reviews and the applied instruction-edit proposal, preserving pre-correction findings. The [live gap register](planned/policy_set_audit_2026-10-21.md) continues to own unresolved work. |
