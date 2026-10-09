# Continuous KB Gather Performance

Contract status: RETROSPECTIVE RECORD
Design/PLAN status: COMPLETED
Implementation status: COMPLETE; human-reviewed, transferred, and revalidated on local master

This note is the single portable development record for the optimization. The
numerical contract, scope, oracle, and performance gates were written and
reviewed before implementation. The original note was a concise PLAN, not a
formally approved full SPEC. This completed record does not claim retroactive
SPEC approval.

## Contract

### Problem and outcome

Continuous pose refinement spends most trial-normal-term time in rotation and
Kaiser-Bessel gathering. Reduce that cost without changing coordinates,
objectives, derivatives, optimizer policy, or accepted endpoints.

### Scope

- Specialize the fixed three-tap Cartesian Fourier gather used by the Cartesian
  calculator.
- Add an independent numerical-oracle test beside the Fourier implementation.
- Validate numerical parity, downstream calculator behavior, the complete fast
  area, and matched production performance.

### Numerical requirements

- Keep `matmul([h,k,0],rotmat)` coordinate generation unchanged.
- Keep `KBWINSZ=1.5`, `KBALPHA=2`, and three taps per axis.
- Keep stencil selection, wrap support, packed/Friedel mapping, accumulation
  order, native Fourier scaling, and fixed-cell derivatives unchanged.
- Keep the generic stencil and packed-gather routines as the independent oracle
  and for other callers.
- Reject the optimization if oracle values or derivatives exceed 32
  single-precision epsilons, downstream tests fail, matched endpoints differ,
  or the declared performance and memory gates fail.

### Ownership and integration

- `simple_cartesian_fourier` owns Cartesian Fourier stencil construction and
  packed KB gathering.
- `simple_cartft_calc` consumes that gather in continuous-pose objective and
  derivative evaluation.
- `simple_cartesian_fourier_tester` owns the independent oracle regression.
- No UI, parameter, workflow, persistence, reconstruction, or optimizer-policy
  owner changes.

### Non-goals

- Do not replace the coordinate `matmul` operation.
- Do not change the KB kernel, stencil, arithmetic ordering, or derivatives.
- Do not remove or modify generic interpolation APIs.
- Do not introduce objective-only trials, boundary-stop policy, new profiling,
  or a new command-line option.

### Human review and authority

The user directed the bounded performance work and reviewed the diff in the
isolated performance worktree before transfer. This note does not authorize a
commit or push. Normal repository commit and push authority still applies.

## Design and PLAN

### Existing pattern

The retained generic path generates a fixed 3x3x3 stencil, materializes full
weight and derivative arrays, and then gathers packed Fourier samples. It is the
numerical oracle.

### Implementation method

`gather_packed_kb3_window_grad` forms three normalized one-dimensional KB weight
and derivative vectors. It accumulates the same 27 packed Fourier samples
directly. It does not materialize the former `w(3,3,3)` and
`dw(3,3,3,3)` arrays. The Cartesian calculator calls the specialized routine.
Other callers continue to use the generic APIs.

### Intentional duplication and future cleanup

`gather_packed_window_grad` and `gather_packed_kb3_window_grad` contain the same
packed/Friedel indexing convention by design. The generic routine applies
caller-supplied weights and derivatives and serves as the independent traversal
and accumulation oracle. The specialized routine generates and consumes the
fixed three-tap KB stencil without materializing 3-D arrays. Separate
packed-gather and analytic-field tests protect the shared packed/Friedel
convention, so the fused-oracle comparison is not the only check on that logic.

A future cleanup can remove the production `gather_packed_window_grad` API if a
test-local brute-force oracle first replaces it with equal or stronger coverage.
That oracle must independently construct the stencil, map positive and negative
Fourier locations, accumulate the value and three derivatives, and cover wrap
and Friedel cases. Do not delete the oracle test with the generic routine. A
shared hot-loop lookup helper is a separate performance change and requires the
same numerical parity and matched benchmark gates before adoption.

### Test acceptance contract

| Level | Owner | Acceptance requirement |
| --- | --- | --- |
| Numerical oracle | `simple_cartesian_fourier_tester` | Fused value and three derivatives agree with the retained generic path within 32 single-precision epsilons. |
| Production consumer | `cart calculator` sub-suite | Existing projector, objective, gradient, sigma, and preparation contracts pass. |
| Area regression | `unit_cart_align3D` | Every registered fast sub-suite passes. |
| Matched performance | Frozen beta-gal comparison | Every pair reduces median worker orientation-search time; the gain exceeds baseline variation; strategy time and memory stay within declared limits; work and endpoints match. |

The oracle fixture covers positive and negative locations, packed/Friedel
access, periodic wrap edges, integer coordinates, both sides of positive and
negative half-grid stencil switches, and an out-of-support location.

Use the canonical repository selectors:

```text
simple_test_exec test=unit_cart_align3D suite=cartesian_fourier
simple_test_exec test=unit_cart_align3D suite=cart_calculator
simple_test_exec test=unit_cart_align3D
```

## Implementation record

### Source changes

- Added the fused fixed 3x3x3 KB gather to `simple_cartesian_fourier`.
- Routed the Cartesian calculator through the specialized gather.
- Added the independent fused-versus-generic oracle test.
- Added no UI, parameter, persistent file, workflow, or build target.

### Validation evidence

The post-rebase unit runs completed successfully on 9 October 2026:

- `cartesian_fourier`: 390 of 390 assertions passed; the fused-gather oracle ran.
- `cart_calculator`: 168 of 168 assertions passed.
- `unit_cart_align3D`: five of five sub-suites and 742 of 742 assertions passed.

The local-master reports identify base commit `29808bf3`, and the executed
`test_fused_kb3_gather_oracle` proves that the uncommitted optimized source was
present in the tested executable. Executable and source provenance for the
matched performance run is retained in the checkout-local diagnostic record.

Three order-balanced beta-gal pairs used fresh copies of one frozen checkpoint.
Every pair reduced median worker orientation-search time by 36.84-37.50%.
End-to-end wall time improved by 20.45-21.86%. Memory remained within the 5%
limit. Every arm completed normally with the same work count, identical endpoint
files, and the same reported cFAR and FSC thresholds.

Decision: retain the fused gather. Detailed run identities, hashes, failed
preflights, absolute paths, and measurements remain checkout-local diagnostic
evidence rather than portable repository documentation.

### Deviations and remaining review

- The development began from a concise PLAN rather than a formally approved
  full SPEC. This record makes that process limitation explicit.
- No required numerical, performance, build, test, or human-review gate remains open.
- The human-approved commit message and commit/push decision remain open.

### Future performance work

The following ideas are separate changes. Each requires its own numerical
contract, independent regression coverage, and matched performance evidence:

- Evaluate only the objective for a trial first. Calculate its gradient and
  Hessian only when later optimizer work requires them.
- Improve SIMD vectorization of the fixed 27-sample gather without changing
  its accumulation semantics.
- Reduce repeated packed-index and wrap calculations while preserving the
  packed/Friedel convention.
- Evaluate several Fourier locations as a batch to increase instruction-level
  and memory-level parallelism.

Do not combine these experiments with the accepted fused-gather change. Lower
precision, including `bfloat16`, is not proposed: the current kernel already
uses single precision, and reduced mantissa precision threatens normalized
weights, derivatives, and LM acceptance decisions without a clear hardware
benefit for this irregular gather.
