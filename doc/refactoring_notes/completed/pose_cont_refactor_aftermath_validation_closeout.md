# Continuous-pose refactoring aftermath: validation closeout

Status: closed and archived on 2026-11-04 by the maintainer's decision to retain
the missing historical validation as explicit limits rather than require
historical reruns for archival. The checks below are not claimed as passed;
no new build or runtime evidence is supplied by this closeout.

Historical record, 2026-11-02: reconstructed from V1 of the [completed aftermath record](pose_cont_refactor_aftermath_issues_completed.md).
Changes C1-C4, tests T1-T2 and documentation D1-D2 were implemented. The missing split record
must not be read as evidence that their requested validation passed.

## Unrecorded historical checks retained as limits

- Debug compilation and the fast gate, including `unit_cart_align3D`.
- The Cartesian library gate (`lib_cart_align3D`, historical E25).
- Continuous refinement (`cont_refine3D_1jxy`, historical N20), after T1 and T2.
- The 1JXY and 6VXX polar workflow guards and `solve3D_addon`; the original record explicitly
  kept the add-on guard open until a current pass was recorded.

T1 covers publication of the current particle field before file-backed fractional assembly;
T2 covers the in-plane continuous phase and the polar route before polishing. Reconcile later
run evidence against those behaviors and its exact source revision, rather than treating an
unrelated later pass as closeout. Use current registry selectors for the historical guard names.

These were the original validation requests. Any future comparison must retain
its commands, exact revision, platform, assertions and logs and account for
superseded workflow behavior. The maintainer closed this historical validation
task with those evidence gaps retained; this cleanup ran no numerical tests.
