# Implementation-note status

Every implementation note belongs to one lifecycle directory:

- [completed](completed/) — implemented features and finished investigations
- [planned](planned/) — proposed or partially implemented work
- [rejected](rejected/) — retired, superseded, or rejected designs

The directory is authoritative. Historical sections inside a note can retain
their original proposal language when that chronology is useful.

The [topic index below](#open-work-by-topic) lists the active records and where to
resume. The [developer documentation index](../for_developers/index.md) connects
these records to refactoring work, current algorithms and policies.

Record conventions:

- Existing `<name>_completed.md` / `<name>_unfinished.md` splits preserve the
  delivered record and one live list of open work. New development keeps its
  contract, plan, review findings and validation criteria in one living note.
- Dated plans and requested run reports retain their measurements and caveats.
  A completed investigation can leave research questions or maintainer-owned
  map assessment open; its classification describes the delivered scope.
- `planned/next_release/` holds proposals deferred past the current release;
  the `rejected/` folder is top level only (no nested `planned/rejected/`).

## Open work by topic

Navigation reviewed 2026-11-02. Status describes the records, with no new runtime
validation implied. Follow each note's current banner and dated evidence before
resuming an old command or experiment.

| Topic | Record | Status and next action |
| --- | --- | --- |
| Polar shift search | [Living plan](planned/polar_shift_search.md), [review](planned/polar_shift_search_review.md), [dynamic centring](planned/polar_shift_search_dynamic_centring.md) | Active implementation; reconcile the latest changes and review findings against the current source. |
| Posterior search | [Pose marginalization](planned/pose_marginalization.md) | Implemented policy with later experiments and rulings; reconcile real-data evidence and superseded controls. |
| Stage scheduling | [Data-driven progression](planned/data_driven_stage_progression.md) | Proposal; decide how stage progression should use measured information. |
| Multi-state regression | [Investigation](planned/multistate_solve3D_regression.md), [fix check](planned/multistate_solve3D_regression_report.md), [stage-5 check](planned/multistate_stage5_prob_check_report.md) | Measurements show partial recovery; diagnosis and scientific closeout remain. |
| Density model | [Living record](planned/atomic_density_model.md), [run report](planned/atomic_density_model_report.md) | Implementation delivered in two runs; unresolved model interpretation and validation remain. |
| Species discovery | [Single open-work list](planned/species_discovery_unfinished.md) | Finish the applicable sensitivity/half-map checks and promotion decisions; historical numbers are in the linked evidence record. |
| Atomic B factors | [B-factor model](planned/atom_bfactor_model.md) | Follow the map-only measurement ruling and assess the density-model interpretation. |
| Atomicity prior | [Integration closeout](planned/single_pcg_integration_unfinished.md) | Prior implemented; reconcile remaining workflow integration and validation. |
| Nanoparticle rejection | [Earlier context run](planned/nano_cavgs_rejection.md), [report](planned/nano_cavgs_rejection_report.md), [separate trajectory workflow](planned/trajectory_cavgs_rejection.md) | Earlier context delivered; the separate workflow has its own proposal and acceptance criteria. |
| Trajectory denoising | [Flex kernel](planned/trajectory_denoise_flex.md) | Implemented; reconcile the remaining validation with the retirement of the alternative modes. |
| Wave-model integration | [Wave summary](planned/wave_model_summary.md) | Operators and simulator delivered; refinement coupling remains to be planned. |
| NU auxiliary candidates | [Blush competition](planned/blush_auxiliary_competition.md) | Proposed experiment; precomputed-map controls precede native inference. |
| NU masking | [Cross-validated threshold](planned/nu_evidence_mask_cross_validated_threshold.md) | Proposal; start with the diagnostic table. |
| Distributed execution | [Coarray SPMD](planned/coarray_spmd_execution.md) | Proposal; complete measurements and lifecycle gates before implementation. |
| Threading | [Safety and budget](planned/thread_budget_and_safety.md), [BLAS configuration](planned/blas_openmp_policy.md) | Policy decisions and real-toolchain validation remain. |
| Resource planning | [Planner closeout](planned/hardware_resource_planning_unfinished.md) | External planner delivered; reconcile remaining calibration and scientific gates. |
| Performance follow-up | [Speed-week closeout](planned/single_state_speed_week_unfinished.md), [multi-state closeout](planned/multistate_speed_run_unfinished.md), [architecture proposals](planned/speed_architecture_proposals.md) | Keep applicable measurements open; dispose explicitly of acceptance criteria for retired search policies. |
| Testing and benchmarks | [Order of work](planned/testing_and_benchmarking_plan.md), [truth-suite record](planned/test_suite_truth_and_robustness.md), [report](planned/test_suite_truth_and_robustness_report.md) | Core truth tooling delivered; dataset work and remaining validation stay open. Current rules live in the test policy. |
| Batch completion | [Diagnostics and timeout](planned/qsys_scheduler_completion_bug_report_unfinished.md) | Exit-code/lost-job handling delivered; periodic reporting and hung-job backstop remain. |
| Chunked solve3D | [Emulation closeout](planned/emulate_solve3D_stream_unfinished.md) | Reconcile supported chunk/publication checks; the removed emulation command's rollback check is retired. |
| Deferred NU atom support | [C-alpha mask](planned/next_release/calpha_nu_evidence_mask.md) | Next-release proposal. |
| Deferred optics | [CTF, aberrations and Ewald](planned/next_release/ctf_refinement_aberrations_ewald_plan.md) | Next-release survey and agreed order of work. |

## Completed investigations and reference material

| Topic | Entry point | Retained follow-up |
| --- | --- | --- |
| Wave optics | [Operators and report](completed/ghostbuster_operators_report.md), [simulator](completed/ghostbuster_scattering_simulator.md) | Operator limitations and future integration remain explicit in the report and wave summary. |
| Polar performance | [Layout benchmark](completed/polar_matching_layout_bench_report.md), [Mac kernels](completed/polar_kernel_codegen_report.md), [x86 kernels](completed/polar_kernel_codegen_x86_report.md) | Adoption decisions and real-run measurements retain their original scope. Candidate sources are linked from the x86 report. |
| Platinum pipeline | [Hand-run pipeline](completed/nano_pipeline_pt_report.md), [one-program workflow](completed/trajectory_solve3D_report.md) | Maintainer map assessment and further rejection research remain recorded. |
| Atomic reconstruction | [Current prior algorithm](../algorithms/atomicity_prior.md), [dated methods](completed/atomicity_prior_methods.md), [AET comparison](completed/aet_codes_and_simple.md) | Current behavior is separated from the original experiments and outlook. |
| Species evidence | [Dated evidence and handoff](completed/species_discovery_status.md) | Next steps are maintained only in the planned species record. |

The lifecycle directories above contain the remaining historical records. Their
length alone is not a reason to discard measurements or decision rationale.
