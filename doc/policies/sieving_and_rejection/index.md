# Sieving and Rejection Policies

This folder collects policy and design documentation for class-average rejection,
size-compatibility filtering, and staged particle sieving workflows.

## Documents

- [ptcl_sieve_policy.md](ptcl_sieve_policy.md): staged coarse/fine particle-sieve policy, chunk lifecycle, submission/collection order, completion semantics, counters, and change guardrails.
- [class_compatibility_policy.md](class_compatibility_policy.md): support-model fit/infer policy for class-size compatibility (`a/b/c` axes), convergence semantics, metrics contract, and test requirements.
- [model_cavgs_rejection.md](model_cavgs_rejection.md): shared class-average quality backend, command modes, feature bank, hard rejects, built-in models, analysis output, and promotion rules.
- [distance_transform_shape_rejection_plan.md](distance_transform_shape_rejection_plan.md): rotationally invariant distance-transform shape evidence plan for model-backed rejection.

## Implementation Pointers

- `src/main/sieve/simple_ptcl_sieve.f90`: staged coarse/fine chunk orchestration and rejection flow.
- `src/main/class/simple_class_compatibility.f90`: class-size compatibility model training/inference.
- `src/main/cavg_quality`: shared class-average quality backend and model-based rejection implementation.
- `src/main/sieve/simple_ptcl_sieve_tester.f90`: policy-level tests for staged sieve lifecycle/query behavior.
- `src/main/class/simple_class_compatibility_tester.f90`: policy-level tests for support-model validity/convergence/inference behavior.
