# Continuous-pose refactoring aftermath: validation V1

Status, 2026-11-02: reconstructed from V1 of the [completed aftermath record](../completed/pose_cont_refactor_aftermath_issues_completed.md).
Changes C1-C4, tests T1-T2 and documentation D1-D2 were implemented. The missing split record
must not be read as evidence that their requested validation passed.

## Runs to reconcile and record

- Debug compilation and the fast gate, including `unit_cart_align3D`.
- The Cartesian library gate (`lib_cart_align3D`, historical E25).
- Continuous refinement (`cont_refine3D_1jxy`, historical N20), after T1 and T2.
- The 1JXY and 6VXX polar workflow guards and `solve3D_addon`; the original record explicitly
  kept the add-on guard open until a current pass was recorded.

T1 covers publication of the current particle field before file-backed fractional assembly;
T2 covers the in-plane continuous phase and the polar route before polishing. Reconcile later
run evidence against those behaviors and its exact source revision, rather than treating an
unrelated later pass as closeout. Use current registry selectors for the historical guard names.

Record commands, revision, platform, assertions and retained log locations here. Any requirement
superseded by later workflow changes needs an explicit disposition. This cleanup ran no numerical tests.
