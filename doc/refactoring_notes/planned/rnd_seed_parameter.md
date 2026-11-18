# Random seed through parameters and job descriptions

**Approved, not implemented.** The maintainer approved the seed key on
2026-10-18 ([testing plan, ruling 6](../../implementation_notes/planned/testing_and_benchmarking_plan.md#2-rulings-of-2026-10-18)).
Updated 2026-11-17 against master `68638be28` plus the working tree. This replaces
the October implementation sketch, which predated the current 64-bit keyed RNG
and restart contract. No RNG code or numerical behaviour changes in this update.

## Current behaviour and purpose

`simple_rnd` reads `SIMPLE_SEED` as an `integer(int64)` root. Legacy `seed_rnd`
advances the intrinsic generator's seed by `7919 * call_count`; the fixed-seed
API resets that count. `parameters%new` currently seeds before parsing inputs.
Unseeded legacy calls obtain operating-system randomness. The keyed generator
uses a cached run seed and logical coordinates (program, pass, iteration,
particle, draw family, reference and ordinal), not thread or partition identity.
Its cache must be initialized serially before parallel reads.

Seeds are already recorded in sampling/restart metadata. `ensure_run_seed`
creates and exports a root when needed; `set_run_seed` adopts a recorded root
when none is set and rejects a conflicting existing root. Consequently the
problem is an implicit environment configuration/propagation path, not an absence
of all provenance. Replace that path with the ordinary command/job interface
while retaining these RNG and restart contracts.

## Contract

1. Add the user key `rnd_seed` with **lossless signed 64-bit integer semantics**.
   Keep presence separate from value: zero is a seed, not an unset sentinel.
   Specify and test the accepted range, including how existing signed test seed
   offsets are handled. Reject malformed or out-of-range input explicitly.
2. Preserve the exact decimal value through parsing, parameter binding, command
   copies, job descriptions, child commands and persisted metadata. The current
   `cmdline` parser has default-integer/`real(dp)` numeric storage and a character
   fallback, and the parameter registry's `add_int` binds default integers.
   Neither is a complete 64-bit contract. Add a narrow exact-decimal path and
   typed 64-bit binding in those owners; never round-trip a seed through default
   integer or floating point. Audit all setters/getters and serialization paths,
   not only the initial parse. Generated argument/UI metadata must agree.
3. Keep the run-root state in `simple_rnd`, reusing its existing lifecycle APIs.
   Parsing must establish the explicitly supplied root before the corresponding
   `seed_rnd` call, without adding or losing calls. Re-entering a commander with
   the same root, or an omitted key in the same process, must not reset the
   intrinsic call counter or generate a new keyed root. Conflicting roots on a
   continuation must still fail; an intentional fresh-run/reset operation must
   be explicit. Audit whether parsing itself draws before moving initialization.
4. Preserve `seed_rnd_fixed`'s local-versus-propagating behaviour and the test
   dispatcher's `seed_offset` arithmetic. The offset is applied once to each
   fixed test base, before propagation. Keep keyed generator version, scope and
   draw coordinates unchanged; do not introduce thread/partition offsets.
5. Every fresh worker receives the established run root explicitly. Extend
   `cmdline%gen_job_descr` and audit launchers that bypass it, including the test
   dispatcher, queue backends, coarray images and persistent workers. Persistent
   processes need the existing run-boundary/reset semantics, not an unconditional
   reset on every sub-command. A missing key in a newly launched process must not
   silently lose an inherited root.
6. On restart, load the recorded root before any keyed work. Adopt it when the
   caller has supplied none; reject a different explicit root. Do not generate a
   replacement before reading restart metadata. Retain existing seed/generator
   metadata and sample-plan replay validation; add new formats only if necessary
   and define compatibility explicitly.
7. Remove the production `SIMPLE_SEED` read/write path only after every consumer
   and launcher has explicit propagation. Unset fresh production runs retain
   stochastic behaviour; record/propagate the root when a keyed run establishes
   one. Do not silently change the legacy unseeded call lifecycle into a fixed
   sequence. Update the hardware planner's remaining uses without reviving its
   retired testing role.

Sequence equivalence is an acceptance criterion to demonstrate, not a conclusion
from moving a setter. In particular, inspect overflow in the legacy call-count
advance and signed seed-offset arithmetic at the supported limits. Any necessary
change to those semantics needs an explicit compatibility decision and tests;
it must not be hidden inside the parameter migration.

## Implementation order and owners

| Area | Work |
| --- | --- |
| `src/main/simple_cmdline.f90` and its tester | Exact 64-bit parse/set/get/copy/serialization and worker job propagation; retain existing behaviour for other keys |
| `src/main/params/simple_parameters*.f90` | Seed value and presence, registry binding, initialization order and validation; regenerate argument metadata through the existing build |
| `src/utils/math/simple_rnd.f90` | Module-owned root lifecycle in place of environment transport, retaining intrinsic and keyed APIs and restart mismatch checks |
| `src/utils/simple_test_utils.f90`, `production/simple_test_exec.f90`, test-driver children | Fixed seeds, offset arithmetic, local/propagating reseeds and fresh-child initialization |
| Queue/coarray/persistent-worker launchers and job serializers | Inventory and prove explicit transport across every process/image boundary that previously inherited the environment |
| `src/main/strategies/search/simple_matcher_smpl_and_lplims.f90` and sample-plan metadata owners | Recorded-root adoption, mismatch rejection and restart replay before draws |
| Hardware planner and runner scripts, `production/CMakeLists.txt` | Remove remaining live environment configuration; keep resource prediction and test roles separate |
| Test policy, UI/help and generated indexes | Document the final public contract and regenerate navigation after implementation |

Inventory all callers of `seed_rnd`, `seed_rnd_fixed`, `simple_seed_base`,
`put_simple_seed`, `ensure_run_seed`, `set_run_seed` and all `SIMPLE_SEED` uses
before editing. Extend the existing owners and suites; do not add parallel RNG
helpers or a second plan. Keep historical run commands and seed settings in
completed notes unchanged: they describe the mechanism used for those runs.

## Validation and completion criteria

- Exact round trips for zero, ordinary seeds, a value above 2^53 (for example
  9007199254740993), supported signed/64-bit boundaries, copies, serialization,
  workers and restart metadata. Invalid syntax and overflow must fail.
- Paired before/after intrinsic sequences for repeated `parameters%new`, local
  fixed reseeding, propagating reseeding, omitted/same-root sub-commands and
  fresh children, including nonzero seed offsets applied exactly once.
- Keyed-draw identity for the same root/logical coordinates across thread and
  partition layouts; unchanged generator version and restart sample plans.
- Restart with an omitted seed adopts the recorded value; the same explicit seed
  succeeds; a different one fails before work. Cover fresh runs in reused workers.
- Check registrations and source headers, build incrementally, run the fast gate
  and affected command/parameter/RNG/project tests, then library/high-level cases
  that exercise child propagation and replay. Compare numerical metrics with the
  pre-change evidence without changing floors to hide a seed regression. Record
  platform and compiler limits; no Linux or full-nightly claim without its output.
- Search `src`, `production`, `scripts` and `.github` for live `SIMPLE_SEED` reads,
  writes and launcher dependencies. Resolve each before declaring the migration
  complete; dated historical evidence may still mention it.

Changing generators, draw identities, numerical floors or per-partition seed
policy is outside this migration. The implementation and validation ledger will
be added to this same note when the work is executed.
