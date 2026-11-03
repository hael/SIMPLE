# Atomicity prior in PCG reconstruction

Source checked 2026-11-02. This reference describes the implemented prior modes; the
[dated methods and FePt study](../implementation_notes/completed/atomicity_prior_methods.md)
preserves the original derivation and measurements. It is not evidence of a new validation run.

## Problem

Given fixed projection poses and an atomic model, regularize the reconstructed density toward
the model while retaining a linear PCG solve on a fixed real-space support.

## Model

Let `G` contain Gaussian atom kernels on the solve support and let `P` project density into
the allowed model subspace. The prior adds `lambda_A * (I - P)` to the normal operator,
penalizing the density outside that subspace. The source supports two choices:

- **Free atom amplitudes:** `P x = G (G^T G)^-1 G^T x`. Each atom may have its own fitted B
  factor; otherwise a common kernel width is estimated from the prior-free map or supplied.
  The sparse Gram solve uses a small diagonal ridge and double-precision inner CG.
- **Fitted model:** per-atom widths, occupancy and posterior intensity define a model `m`.
  Only its overall scale is free: `P x = m * (m^T x) / (m^T m)`. This is a one-dimensional
  subspace, so the independent-amplitude interpretation of the original methods note does
  not describe this mode.

## Algorithm

The strategy first solves the prior-free even/odd pair. It reads the atom positions from
`pdbfile`, converts them to the cropped solve grid and removes duplicates or atoms outside
the support. A companion CSV with `BFAC` supplies per-atom widths; `POSTERIOR_INTENSITY`
selects the fitted-model mode, with `OCCUPANCY` used when present.

The strategy estimates the relative prior strength from the base pair's noise, signal and
explained signal fraction. A positive `lambda_atm` overrides that estimate; zero retains the
estimate. `atm_prior=no` selects the prior-free behavior. Both half operators receive the
same model subspace, and each half is solved again from zero with the atomic-prior budget.
The fitted mode uses its supplied model coefficients for the active projection. Both modes
also fit free coefficients through each half's observation operator for diagnostics; the
residual map uses the active projection in fitted mode and the half-specific models otherwise.

## Rationale and evidence

The model can stabilize poorly measured density, but both regularized halves share it.
Their agreement is consequently not independent resolution evidence. Keep the prior-free
`_base` pair and its FSC, and assess gains against ground truth, held-out images or an
independent model. The strategy writes residual/model diagnostics and `_atmprior.txt`;
the historical FePt numbers retain the scope of the original experiment.

Open workflow integration and validation are tracked in
[the integration backlog](../implementation_notes/planned/single_pcg_integration_unfinished.md).

## Implementation

- [Prior domain object](../../src/main/volume/simple_pcg_atomic_prior.f90): `new`, `build`,
  `project`, `apply`, `solve_gram` and `estimate_lambda`.
- [PCG strategy](../../src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90):
  `prepare_atomic_prior_on_pair`, `read_species_prior_fields` and `write_atomic_prior_products`.
- [PCG reconstruction policy](../policies/3D/reconstruct3D_pcg_policy.md) owns the solver contract.
