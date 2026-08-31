# Reference-Guided 3D Classification Policy

This document defines `classify3D_refs`, the canonical workflow for
competitively classifying particles against references with independent or
untrusted alignment provenance. No historical command alias is registered or
routed.

Related policies:

- [refine3D_policy.md](refine3D_policy.md)
- [refine3D_states_policy.md](refine3D_states_policy.md)
- [sigma_calculation_policy.md](sigma_calculation_policy.md)
- [importance_sampling_fractional_update_policy.md](importance_sampling_fractional_update_policy.md)

Primary implementation:

- `src/main/ui/simple/simple_ui_refine3D.f90`
- `src/main/exec/simple_exec_refine3D.f90`
- `src/main/commanders/simple/simple_commanders_refine3D.f90`
- `src/main/simple_external_reference_pose_initialization.f90`
- `src/main/simple_refine3D_stage_plan.f90`

## 1. Scientific Scope

`classify3D_refs` answers whether and how particles match a supplied set of
structural hypotheses. It does not assume that current particle poses or state
labels were established against those references. The workflow therefore
initializes pose/state correspondence with cross-correlation before it uses a
noise-normalized Euclidean objective.

The primary output is a hard state and pose assignment for every active
particle. Native-sampling state reconstructions and their normal artifacts are
also produced so the classification can be inspected scientifically.

## 2. Input Contract

- `nstates >= 2` is required.
- A complete external `vol1..volN` set is required.
- Partial reference sets are rejected.
- Every reference must exist and match the project particle box and sampling.
- The active particle population must be nonempty.
- Existing state labels may seed the run but are not treated as evidence that
  reference-conditioned initialization has occurred.

Same-lineage conformational refinement belongs to `refine3D_states`, not this
workflow.

## 3. Mandatory CC Pose Initialization

The workflow always calls the shared external-reference pose-initialization
service before Euclidean refinement:

1. prepare every state reference at a common 15 Å low-pass;
2. select one cohort capped at 100,000 active particles;
3. run exactly one broad `objfun=cc` matching pass;
4. keep all supplied references fixed with `volrec=no` and `trail_rec=no`;
5. commit one hard pose and state assignment per selected particle;
6. consolidate Euclidean residual sigma contributions produced at those
   assignments;
7. reconstruct data-derived checkpoint maps from the assigned cohort;
8. use those checkpoint maps, rather than the external inputs, as the
   authoritative references for Euclidean refinement.

The 100,000-particle cap is a one-pass initialization cohort, not the first
fraction of a coverage schedule. There is no repeated CC pass and no
missing-particle sweep during pose initialization.

This transition uses global grouped residual statistics. That keeps the capped
cohort sufficient even when it does not contain every acquisition group; a
group-specific bootstrap would require coverage that this one-pass contract
deliberately does not promise.

The shared implementation is
`initialize_poses_against_external_references` in
`simple_external_reference_pose_initialization`. `abinitio3D` input-volume
starts and `refine3D_auto ref_pose_init=cc` use the same service.

## 4. Euclidean Classification

After the checkpoint, the wrapper switches to `objfun=euclid` and uses full
state-pooled probabilistic matching. The effective internal policy is
`refine=prob_neigh`, `prob_neigh_mode=state`: every state is evaluated over the
same pooled candidate geometry and the committed update remains a hard state,
pose, and shift assignment.

The wrapper uses stochastic fractional updates when the active count exceeds
the normal sample target and full updates otherwise. The automatic target is
10,000 particles per state, capped at 100,000 total particles per iteration.
Balanced class sampling is used only when appropriate metadata exists.

## 5. Frequency and Reconstruction

`lpstart` and `lpstop` define one common frequency schedule for all states.
The workflow consumes `simple_refine3D_stage_plan` in short blocks, applying
the returned low-pass, crop, translation limit, and iteration range to each
base-`refine3D` call.

After the frequency march, a missing-update pass ensures every active particle
has a committed reference-conditioned assignment. Final maps are reconstructed
and postprocessed from all active particles at native sampling. External maps
are never blended into the authoritative final reconstruction.

## 6. Ownership

The wrapper owns provenance policy, the CC-to-Euclidean transition, sampling,
frequency planning, coverage, and final command construction. Search modules
own candidate scoring and pose/state updates. Reconstruction and volume-domain
modules own partial sums, assembly, half maps, FSCs, automasking, filtering,
and final state products.

## 7. Failure Policy

The command fails before Euclidean refinement when:

- the reference set is incomplete or incompatible;
- the CC pose-initialization cohort is empty;
- a checkpoint state is empty or invalid;
- sigma consolidation or checkpoint reconstruction is missing;
- the requested frequency interval is invalid;
- active-particle assignment coverage cannot be completed.

## 8. Validation

User-side validation must verify:

- one and only one fixed-reference CC pass;
- a cohort of at most 100,000 particles and a common 15 Å bandwidth;
- unchanged prepared external-reference hashes throughout that pass;
- residual sigma availability before the first Euclidean block;
- data-derived checkpoint state populations and maps;
- global probabilistic matching and monotonic frequency marching;
- final all-particle update coverage;
- native-sampling maps, half maps, FSCs, and diagnostic products;
- shared-memory and distributed parity.

Compilation and runtime tests are performed by the user. No Linux or BOX
result is recorded as passing without observed output.
