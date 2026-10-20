# FLEX PCA architecture refactor: open work

Status (2026-10-20): the refactor itself is implemented. What is open is the proof that it changed
no result, plus one optional cleanup. Background, measurements, the target design and the full
phase ledger are in
[flex_pca_architecture_audit_and_refactoring_plan_2026_09_18_completed.md](../completed/flex_pca_architecture_audit_and_refactoring_plan_2026_09_18_completed.md).

## What we learned

Phases 0b to 11 landed (settings factory, application layer, composed fit state, posterior service,
one fit engine, M-step backend, state reconstruction split, argument contracts, folder layout under
`src/main/flex/{run,fit,states}`). Moves were checked procedure by procedure, and a seeded toy run
was byte-identical to `master` for the shared-memory, PCG M-step, PCG states and gridding M-step
arms. The distributed arm agrees to the last ulp. None of this used the planned fixtures, so the
"identical" gate in every phase is still unverified at the level the plan asked for.

## Open item 1: the parity harness (phase 0a)

Not done. FLEX (the `flex_pca` heterogeneity estimator, 3-D) needs cluster time to capture three of
the eight fixtures.

A fixture is a captured reference output for a fixed input, stored under the `heterogeneity` CTest
area, which later code is compared against (bit-for-bit for integer and byte payloads, to a stated
tolerance for floating-point statistics).

| Fixture | Configuration | Captured reference |
|---|---|---|
| F1 | shared-memory, Cartesian, gridding M-step, single fit, 2 iterations | E-step sufficient statistics per batch (`G`, `b`, `c`, `e_mm`, `myv`), posterior moments, contrast, likelihood |
| F2 | as F1 with polar E-step | same as F1 |
| F3 | distributed master plus 2 workers, F1 settings | probe-part bytes from each worker, master fold result, embedding sufficient-statistics parts |
| F4 | paired fit, gridding and PCG M-step (two captures) | per-iteration basis, ridge, convergence trace, paired-merge snapshot |
| F5 | resume from embedding cache | embedding cache bytes, deconvolution block, state labels |
| F6 | state reconstruction, gridding backend, shared | combined/even/odd maps, FSC, filtered maps, weight tables, manifest |
| F7 | state reconstruction, PCG backend, shared | as F6 plus per-state PCG table |
| F8 | state reconstruction, distributed | as F6 plus per-part state weights and gridding parts |

F3, F4 and F8 need cluster runs. Budget them before starting.

Binding rules:
- Capture on `master` at commit `6642aaefa` (the audit baseline, before the refactor), then replay on
  the refactored code. A capture taken after the refactor proves nothing.
- Register the harness in `unit_heterogeneity` (fast subset) and `lib_heterogeneity` (full).
- A phase passes only on evidence that ownership moved and the applicable gates hold, never because
  files became shorter.

Acceptance checks (the "identical" gates, run per area):
- E-step: Cartesian and polar paths give the same sufficient statistics, posterior moments, contrast
  and likelihood (F1, F2).
- Probe parts and distribution: same part bytes, same master fold, master-only publication (F3).
- Fit engine and M-step: same rounds, convergence trace, basis, ridge, part payloads and stopping
  per backend (F4).
- Resume: same embedding cache bytes, deconvolution block, state labels (F5).
- State reconstruction: same maps, FSC and filter products, filenames and project metadata per
  backend (F6 to F8).

Note: the testing work since September (the truth suite) may cover part of this. Check what it
already pins before capturing anything new.

## Open item 2: optional late physical splits (section 8.1 of the completed note)

Not done, and optional. Two large numerical leaf files could be split mechanically, with no change
in public behaviour, once the contracts are stable (they are):
- `simple_flex_reconstructor_latent_ops` (about 1,000 lines): projection, backprojection, coupled
  grid solve, observation preparation.
- `simple_flex_pca_pcg` (about 2,600 lines): geometry, accumulation, operator, solve, support, test
  support submodules.

Do this only after open item 1 gives a trustworthy comparison, since a split is exactly the kind of
move the fixtures are meant to guard.
