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

Navigation reviewed 2026-11-04. Status describes the records, with no new runtime
validation implied. Follow each note's current banner and dated evidence before
resuming an old command or experiment.

| Topic | Record | Status and next action |
| --- | --- | --- |
| Polar shift search | [Laplace width](planned/polar_shift_laplace_width.md) | Proposal split from the completed shift search (2026-11-03); Hans's call on conditional or marginal width. |
| Stage scheduling | [Data-driven progression](planned/data_driven_stage_progression.md) | Proposal; decide how stage progression should use measured information. |
| Multi-state regression | [Investigation](planned/multistate_solve3D_regression.md), [fix check](planned/multistate_solve3D_regression_report.md), [stage-5 check](planned/multistate_stage5_prob_check_report.md) | Measurements show partial recovery; diagnosis and scientific closeout remain. |
| Learned pose assignment | [Pose student](planned/learned_pose_state_assignment.md) | Proposal of 2026-11-05, revised 2026-11-06; nothing built. Start with the deterministic teacher and the feature module on the Ribosembly subset (its section 11, phase 0); Hans rules the open decisions of section 12. |
| Species discovery | [Single open-work list](planned/species_discovery_unfinished.md) | Finish the applicable sensitivity/half-map checks and promotion decisions; historical numbers are in the linked evidence record. |
| Atomicity prior | [Integration closeout](planned/single_pcg_integration_unfinished.md) | Prior implemented; reconcile remaining workflow integration and validation. |
| Nanoparticle rejection | [Separate trajectory workflow](planned/trajectory_cavgs_rejection.md) | The earlier nano-context run and report are completed; the separate workflow retains its own proposal and acceptance criteria. |
| Wave-model integration | [Wave summary](planned/wave_model_summary.md) | Operators and simulator delivered; refinement coupling remains to be planned. |
| NU auxiliary candidates | [Blush competition](planned/blush_auxiliary_competition.md) | Proposed experiment; precomputed-map controls precede native inference. |
| NU masking | [Cross-validated threshold](planned/nu_evidence_mask_cross_validated_threshold.md) | Proposal; start with the diagnostic table. |
| Distributed execution | [Coarray SPMD](planned/coarray_spmd_execution.md) | Proposal; complete measurements and lifecycle gates before implementation. |
| Probabilistic alignment pipeline | [Batch pipeline](completed/prob_align_batch_pipeline.md) | Done by the run `prob_pipeline` on the Dell (2026-11-04/05): keyed draws, one batch pipeline for shared memory and the parts, the table programs retired; the Progress section holds the evidence and the measurements. |
| Threading | [Safety and budget](planned/thread_budget_and_safety.md) | Policy decisions and safety work remain planned; the delivered BLAS configuration is archived below with its qualification limits. |
| Performance follow-up | [Multi-state closeout](planned/multistate_speed_run_unfinished.md), [architecture proposals](planned/speed_architecture_proposals.md) | Keep applicable measurements open; dispose explicitly of acceptance criteria for retired search policies. The original speed-week checklist is archived below. |
| Testing and benchmarks | [Order of work](planned/testing_and_benchmarking_plan.md) | Dataset work, benchmarks and remaining extensions stay planned; the delivered truth-suite implementation and report are archived below. Current rules live in the test policy. |
| Batch completion | [Diagnostics and timeout](planned/qsys_scheduler_completion_bug_report_unfinished.md) | Exit-code/lost-job handling delivered; periodic reporting and hung-job backstop remain. |
| Chunked solve3D | [solve3D_chunks open validation](completed/solve3D_chunks_completed.md) | Three validation runs remain (its section 9): the whole-stack ribosome rerun, the boundary cases, the published project. |
| Deferred NU atom support | [C-alpha mask](planned/next_release/calpha_nu_evidence_mask.md) | Next-release proposal. |
| Deferred optics | [CTF, aberrations and Ewald](planned/next_release/ctf_refinement_aberrations_ewald_plan.md) | Next-release survey and agreed order of work. |

## Completed investigations and reference material

| Topic | Entry point | Retained follow-up |
| --- | --- | --- |
| Wave optics | [Operators and report](completed/ghostbuster_operators_report.md), [simulator](completed/ghostbuster_scattering_simulator.md) | Operator limitations and future integration remain explicit in the report and wave summary. |
| Polar performance | [Layout benchmark](completed/polar_matching_layout_bench_report.md), [Mac kernels](completed/polar_kernel_codegen_report.md), [x86 kernels](completed/polar_kernel_codegen_x86_report.md) | Adoption decisions and real-run measurements retain their original scope. Candidate sources are linked from the x86 report. |
| Speed week | [Implementation record](completed/single_state_speed_week_completed.md), [closeout](completed/single_state_speed_week_closeout.md), [run report](completed/single_state_speed_week_report.md) | Delivered work closed by maintainer decision. Remaining high-resolution profiling, optional shared-memory reuse and representative benchmarks are explicit limits; later pipeline and multi-state work keep their own plans. |
| Platinum pipeline | [Hand-run pipeline](completed/nano_pipeline_pt_report.md), [one-program workflow](completed/trajectory_solve3D_report.md) | Maintainer map assessment and further rejection research remain recorded. |
| Atomic reconstruction | [Current prior algorithm](../algorithms/atomicity_prior.md), [dated methods](completed/atomicity_prior_methods.md), [AET comparison](completed/aet_codes_and_simple.md) | Current behavior is separated from the original experiments and outlook. |
| Species evidence | [Dated evidence and handoff](completed/species_discovery_status.md) | Next steps are maintained only in the planned species record. |
| Nano rejection context | [Implementation and investigation](completed/nano_cavgs_rejection.md), [run report](completed/nano_cavgs_rejection_report.md) | Run closed; phase 4 cancelled by maintainer ruling. Findings and possible improvements are retained without reopening the delivered scope. |
| Atomic density model | [Implementation and run record](completed/atomic_density_model.md), [report](completed/atomic_density_model_report.md) | Both runs delivered; unresolved model interpretation, known failure cases and validation remain explicit limits of the completed implementation. |
| Atomic B factors | [Rendering and fitted-model record](completed/atom_bfactor_model.md) | Products and reconstruction feedback delivered through the density-model work. Early controls are historical; physical interpretation and targeted validation remain open limits. |
| Truth-suite implementation | [Completed phases 0-9](completed/test_suite_truth_and_robustness.md), [run report](completed/test_suite_truth_and_robustness_report.md) | Delivered tooling and run evidence; dataset work and remaining extensions stay in the [testing and benchmarking plan](planned/testing_and_benchmarking_plan.md). No new passing-run claim. |
| Trajectory denoising | [FLEX kernel](completed/trajectory_denoise_flex.md) | Implemented and now the only trajectory-denoising path. Targeted validation and possible improvements are retained; comparisons with retired alternatives are historical. |
| Posterior search | [Pose marginalization](completed/pose_marginalization.md) | Posterior draws, annealing, diagnostics and cone-restrained weighted insertion delivered. Superseded controls and later proposals remain dated history; calibration and real-data limits are explicit. |
| Resource planning | [Implementation record](completed/hardware_resource_planning_completed.md), [validation closeout](completed/hardware_resource_planning_validation_closeout.md) | External planner delivered. Full-workload/four-node/production qualification remains unmeasured; native JSON consumption and optional extensions are not implemented. |
| BLAS configuration | [OpenMP BLAS policy](completed/blas_openmp_policy.md) | Configure/build logic implemented. Recorded probes use stand-in libraries; real-toolchain and workload qualification remain unrecorded. |
| Polar shift search | [Living plan, closed](completed/polar_shift_search.md), [review](completed/polar_shift_search_review.md), [dynamic centring](completed/polar_shift_search_dynamic_centring.md) | Closed 2026-11-03; the items left open are listed at the top of the plan, the Laplace width is planned. |

The lifecycle directories above contain the remaining historical records. Their
length alone is not a reason to discard measurements or decision rationale.
