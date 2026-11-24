# Deterministic multi-state assignment gate

## Contract

Replace the automated `multistate_workflows` pipeline with one focused
`refine3D` correctness gate. The gate uses the existing three-state
compositional beta-galactosidase fixture (50/30/20 percent) at its planted
poses and labels, a residual noise estimate, and one fixed reference per state
reconstructed from that planted assignment. It runs one production
`refine=prob_neigh`, `objfun=euclid` assignment block with a full update.

The shared-memory result must:

- preserve every planted state with the existing strict assignment,
  per-state recall, and population bounds;
- leave every estimated state populated and produce one map for each;
- record the expected three-state residual-`sigma` then `prob_neigh` path;
- preserve the planted poses and shifts within the established
  `refine3D_states` ceilings; and
- discriminate the recovered maps while rejecting copies of one map and
  copies of their mean as distinct states.

The same seeded input is then run through two local workers. Its path and maps
must be present and its particle state assignments must agree exactly with the
shared-memory run.

## Implementation

1. `suite=assignment` of the retained `multistate_workflows` commander owns
   the fixture, path checker, truth metrics, and map controls.
2. Only that suite is registered as the `multistate_assignment` high-level
   CTest and nightly entry. The five long workflow routes remain available by
   hand.
3. The shared run uses four threads; the distributed comparison uses two
   workers of two threads and reserves five CPUs at its busiest point.
4. The fixed CTest budget is 46 entries, including 23 high-level entries.

## Review findings

- The retired assignment test used `objfun=cc`, searched poses from a
  noise-model-free start, and no longer represented production. The new gate
  uses the current Euclidean residual sigma2 path and starts at planted poses.
- Atomic truth maps did not provide the numerical scale/support of a normal
  refinement handoff. References reconstructed once from the planted labels
  produced the intended deterministic boundary: the gate tests refinement,
  not reference preparation or de novo discovery.
- The former `multistate_workflows` CTest coupled five long, partly de novo
  routes. Their convergence profiles remain useful evidence but are not a
  stable single-run correctness boundary.
- `simple_commanders_test_highlevel_msworkflows.f90::evaluate` already owns
  the required label, population, map, pose, shuffled-label, copied-map, and
  mean-map checks; the replacement reuses that contract.

## Validation

- `ctest --test-dir build -R '^multistate_assignment$' --output-on-failure`:
  passed in 107.65 s from the registered working directory.
- Scientific metrics: ARI, mapped accuracy, and minimum recall 1.0; exact
  populations 750/450/300; discrimination margin 0.456; median rotation error
  1.61 degrees; median shift error 0.050 pixels; all 1 500 distributed labels
  identical to shared memory; both copied-map controls passed.
- Incremental Release build: passed, with 46 configured CTest entries matching
  the budget.
- Fast gate: 15/15 passed.
- Registry (23 entries), descriptor, FLEX DAG, and `git diff --check`: passed.
- Incremental install: passed. The pre-existing optional `fortls` warning
  (`pip._internal` unavailable) did not affect the build or install.
- Linux validation passed on lb (registered CTest, 136.83 s) and on x3 under
  Slurm: one node on `cn128` (job `55694754`, 139.02 s) and two independent
  replicas on `cn117` and `cn120` (job `55694755`, 123.9/124.1 s). The x3
  clean build, fast gate, and install also passed in job `55694753` on `cn134`.
- After upstream geometric-neighbourhood changes through `c1f255d20`, lb and
  x3 each passed ten runs with `seed_offset=0..9`. lb ran serially in
  122.9--138.2 s; x3 Slurm array `55695214` ran across five nodes in
  122.1--164.2 s. All 20 runs reported accuracy 1.0, exact shared/two-worker
  state agreement, the final planted-state PASS, and normal stop. The refreshed
  x3 build (`55695094`), registered one-node run (`55695095`), and two-node run
  (`55695096`) also passed.
