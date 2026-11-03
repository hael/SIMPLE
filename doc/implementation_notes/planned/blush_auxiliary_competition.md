# Blush half maps as auxiliary candidates in NU competition

Proposed experiment, 2026-11-02.

Include an independently Blush-denoised even/odd pair in SIMPLE's existing
nonuniform filtering competition, alongside the static low-pass bank and the
currently eligible ML-regularized pair. Let local cross-half prediction error
decide where Blush contributes. The hypothesis is that the learned estimator
can improve weak or anisotropically resolved regions while the conventional
candidates remain available wherever they predict the data better.

Each auxiliary pair carries its own effective resolution in Angstrom through
`aux_resolutions`. Use that value to handle the Blush candidate in the same
way as the ML-regularized candidate. The first experiment extends the number
of auxiliary competitors while reusing the existing resolution conventions.

Start with precomputed Blush pairs and gold-standard `refine3D`. Establish
scientific value before integrating native inference. Keep this as the single
living record of the contract, implementation decisions and results.

## 1 Candidate scoring contract

Let `E` and `O` be the unregularized base halves used by NU. Apply the same fixed
pretrained denoiser separately to them, producing `B_E = D(E)` and `B_O = D(O)`.
Blush preprocessing and spectral controls are part of `D`. Score this pair
with the existing `image%nu_objective`:

```text
C_Blush(v) = ((E(v) - B_O(v))^2 + (O(v) - B_E(v))^2) / (2 sigma_0^2)
```

Use the same global `sigma_0`, spherical support from `mskdiam`, and
mask-normalized objective smoothing as the other candidates. Every comparison
must smooth both raw cost fields at the same physical radius. Candidate
generation must preserve common dimensions, sampling, origin and density scale.
Do not fit a candidate's amplitude to the opposite half before scoring it.

The argument does not require a linear estimator. For a fixed estimator,
independent zero-mean target noise makes the expected unnormalized prediction
loss equal to the estimator's signal error plus a target-noise constant.
This motivates scoring a nonlinear denoiser against the opposite raw half.
Scoring `B_E - B_O` alone would instead reward agreement between denoised maps.

The selected label is shared by both output halves: a Blush-labelled voxel
takes its even value from `B_E` and its odd value from `B_O`. Keep the exact
candidate pair that was scored available for composition. The merged reference
is formed after composing the halves.

The statistical argument applies to a fixed candidate under the independence
assumptions. Selection reuses the noisy targets, and refinement, shared FSC
controls and masking introduce further dependence. A lower selected cost is
therefore a diagnostic, not an independent validation of the final map.

## 2 Findings from the current implementation

The source baseline for these findings is `bfef77378` plus the existing
working-tree changes on 2026-11-02. The governing production contract remains
the [NU filtering policy](../../policies/NU/nonuniform_filtering_policy.md).

| Owner | Current behavior and consequence |
| --- | --- |
| [simple_nu_state_filter](../../../src/main/volume/simple_nu_state_filter.f90) | Assembly supplies the unregularized base pair and, under `ml_reg=yes`, the ML pair. This is the shared integration point for gridding and PCG. |
| [simple_nu_filter_bank](../../../src/main/nu_filt/simple_nu_filter_bank.f90) | The array interface retains only the finest supplied auxiliary pair. It admits that pair only at or beyond the finest retained static rung. Passing two pairs does not create two competitors. |
| [simple_image_calc](../../../src/main/image/simple_image_calc.f90) | `nu_objective` already scores arbitrary candidate pairs against the raw halves with the required squared prediction error. |
| [simple_nu_filter_potts](../../../src/main/nu_filt/simple_nu_filter_potts.f90) | Pairwise penalties use candidate coordinates, but the no-finer cap uses numeric label indices. Equal-bandwidth estimators need explicit treatment. The estimated smoothing strength also depends on candidate cost margins. |
| [simple_nu_filter_apply](../../../src/main/nu_filt/simple_nu_filter_apply.f90) | Composition overlays one stored auxiliary pair at one label. Extend this mapping to all admitted auxiliaries. |
| [simple_nu_filter_state](../../../src/main/nu_filt/simple_nu_filter_state.f90) | Evidence accumulation takes the minimum over the competing costs, so adding Blush would also change the evidence envelope unless its contribution is explicitly excluded. |

Other assumptions to audit are the 24-column candidate cap, allocation and
cleanup, ML-specific reporting, `_nu_locres`, and the matching low-pass handoff.
The handoff currently uses the last label as the finest bank member,
independently of which candidates won voxels. Auxiliary storage and objective
columns add memory; measure the full setup, optimization and composition cost.

## 3 Initial experiment boundaries

Generate one Blush pair per state from that state's own base halves. Retain the
published spectral-trailing procedure initially: determine the cutoff from
pre-denoising FSC and limit the denoised contribution about two Fourier shells
below its 0.143 crossing. That cutoff supplies the effective resolution for a
strictly trailed Blush pair. Where the published procedure mixes conventional
regularized reconstruction above the network's range, record both the crossover
and the effective resolution of the resulting hybrid pair. Do not estimate
achieved resolution solely from FSC between denoised halves: shared learned
features can inflate their agreement.

The standard published network uses 1.5 A sampling and 64-cubed patches;
deployment must read the selected checkpoint's metadata rather than assume
these values for every model.

Preserve the current static bank and spherical NU support. Initially apply
the same resolution-based auxiliary admission rule to ML and Blush: the pair
must be at least as fine as the finest retained static rung. Log each pair's
resolution and admission separately. A trailed Blush pair can fail this gate;
choose at least one qualifying snapshot for the initial positive test. Broader
admission inside the ladder is a subsequent experimental variation, rather
than a reason to assign the Blush pair an artificially finer resolution.

Keep envelope generation and application based on the baseline bank without
Blush for this first comparison, including any NU evidence diagnostic. This
isolates candidate selection from a simultaneous change in masking. The
separate compact evidence/replay API currently rejects auxiliary pairs; do not
relax that contract as part of this experiment.

Keep FSC and resolution diagnostics on the particle-derived halves before
Blush application. Blush influences alignment references; final primary maps
are reconstructed conventionally from the resulting particle assignments.
Fractional-update trailing continues to carry raw accumulators. Blush spectral
trailing is a separate frequency control and never replaces that mechanism.

Begin with gridding, independent half-set references, one state and fixed
precomputed inputs. Extend to iterative refinement and then verify PCG and
distributed parity. Defer merged-reference `nonuniform_lpset`, staged `solve3D`,
and PCG solvent-prior apply-pair behavior until their contracts are addressed.
`solve3D` has no gold-standard stage, so its FSC cannot automatically receive
the same interpretation as the initial refinement experiment.

## 4 Implementation sequence

### 4.1 Generalize the auxiliary competition

Retain each admitted pair with an explicit label-to-storage mapping. Separate
the stored pair's identity from its resolution value, using the existing
`aux_resolutions` representation. Keep every eligible auxiliary, rather than
choosing only the finest one. Order the candidates by their supplied resolutions
so the existing coarse-to-fine comparisons remain meaningful.

Use the same resolution-to-smoothing and auxiliary Potts-coordinate rules as
for the ML member. Resolve ties deterministically in favor of an existing
conventional candidate. For multiple auxiliaries at the same effective
resolution, the Potts cap must permit eligible siblings regardless of storage
order. Preserve conventional-only behavior when Blush is absent, and measure
how the expanded bank changes the estimated spatial-prior strength.

Use the finest admitted resolution for the matching handoff, maintaining the
existing last-label convention through correct ordering. When Blush wins a
voxel, use its supplied resolution in the selected-bandwidth diagnostic just
as for ML. Record estimator identity separately; that scalar is the candidate's
effective resolution, not a new voxelwise resolution measurement.

### 4.2 Test precomputed pairs

Add a controlled path for supplying registered Blush half-map pairs to the
existing assembly/NU workflow. Use SIMPLE's normal typed parameter and command
metadata paths if exposed through the CLI. Validate both halves, state,
dimensions, sampling and association with the base reconstruction; reject
stale or mismatched inputs. This stage can use the original Python inference
outside SIMPLE to produce the comparison data.

Keep the composite outputs under the existing NU reference contract. Record
the candidate sources, model identity, preprocessing, spectral controls and
bandwidths in provenance. Reference publication and subsequent loading must
recognize the complete product set and avoid accidental extra filtering.

### 4.3 Add native inference if the experiment succeeds

The preferred route is `Fortran -> ISO_C_BINDING -> small C wrapper -> ONNX
Runtime`. Export the actual checkpoint once using developer-side Python;
the checkpoint contains its own model definition, weights and metadata.
Installed SIMPLE would consume a versioned model artifact through a native
inference library. ONNX export compatibility remains an engineering gate.

Preserve the full inference procedure: Fourier resampling and amplitude
scaling, local-standard-deviation input, volume and patch normalization,
padding, overlapping patch weights, spherical taper, and spectral crossover.
The two network inputs are density and local standard deviation; they are
not the two experimental half-maps. Reuse model sessions and buffers across
halves and states, control CPU threading, and assign GPU work at assembly.
Compare CPU and CUDA execution against the original inference before choosing
precision or performance shortcuts.

## 5 Validation and decision criteria

Use fixed fixtures and seeds, preserve the starting assignments and masks,
and separate candidate-generation effects from subsequent pose changes.

| Comparison | Purpose |
| --- | --- |
| Current NU bank with its eligible ML member | Scientific and numerical baseline. |
| The same bank plus one Blush pair | Primary test of local competition. |
| The same bank plus a conventional control at the Blush comparison scale | Check whether an apparent gain comes from expanding the candidate bank. |
| Blush reference without local competition, with the same spectral controls | Optional attribution of the selector's contribution. |

First compare the methods on identical saved halves. Record selection before
and after spatial smoothing, Blush occupancy in signal and solvent, paired
cost margins against the best conventional candidate at a common radius,
bandwidth maps, boundary artifacts, wall time and peak memory. Candidate
occupancy alone is not a success criterion.

Use synthetic common signal with independent noise for truth-based assessment,
including a phantom whose regions favor different estimators. Include pure
noise, an unhelpful or biased auxiliary, and duplicated-candidate controls to
expose selection bias and unintended changes to spatial regularization. Evaluate
the selected outputs against known truth and, where useful, fresh independent
noise realizations rather than only the targets used for selection.

Then run matched gold-standard refinements with conventional NU and the hybrid
bank. Assess the final conventional reconstructions using truth where
available, independent structural evidence on real data, pre-denoising FSC,
directional quality, weak-domain recovery, solvent artifacts and convergence.
Repeat representative starts. If precomputed maps are used through a live
refinement, regenerate them from each iteration's own base halves.

Implementation checks belong in the existing `unit_reconstruction` area under
the [test policy](../../policies/test_environment_policy.md): closed-form
cross-half costs; simultaneous scoring and composition of two auxiliaries;
candidate-order invariance away from declared ties; same-scale sibling
eligibility under the cap; metadata validation; and identity with Blush
disabled. Preserve the baseline evidence envelope for fixed inputs. Add
meaningful serial/thread comparisons where parallel code changes, then run
the affected high-level workflows and the normal build/fast gate.

Proceed to native integration if the hybrid improves final conventional maps
or refinement robustness over NU across representative weak-signal cases,
without systematic loss of valid density or unsupported resolution gains,
and with acceptable measured runtime and memory. Preselect datasets,
evaluation measures and acceptable tradeoffs before tuning the competition.
Lower training-pair cost or cleaner denoised maps alone do not meet this gate.

## 6 Open decisions

1. Choose the first saved half-map pairs and matched refinement datasets,
   including at least one weak-signal case and a well-resolved control.
2. Pin the effective-resolution convention if a hybrid Blush/ML pair is used.
   A strictly trailed Blush pair initially uses its known cutoff. Start with
   one Blush candidate; test additional bandwidths or admission inside the
   static ladder only if the first comparison warrants it.
3. Set the equal-resolution tie rule and spatial-prior calibration so estimator
   ordering cannot create an advantage. Resolve this before interpreting
   selection occupancy.
4. Decide how a future PCG solvent-prior apply pair should interact with a
   stored Blush candidate. The current separate scoring/apply-pair contract
   does not define that extension.

## 7 References

- Kimanius et al., [Data-driven regularization lowers the size barrier of
  cryo-EM structure determination](https://www.nature.com/articles/s41592-024-02304-8),
  Nature Methods, 2024. Pretrained iterative denoising, spectral trailing and
  final conventional reconstruction.
- RELION checkout reviewed at `afdd8167`: `src/backprojector.cpp`,
  `BackProjector::externalReconstruct`, for the Fourier-data handoff.
- relion-blush checkout reviewed at `ab48e63`: `relion_blush/util.py`,
  `get_reconstructions`, `apply_model` and `install_and_load_model`;
  `relion_blush/command_line.py`, `refine3d`, for the deployed procedure.
- [ONNX Runtime native C API](https://onnxruntime.ai/docs/get-started/with-c.html)
  and [CUDA provider](https://onnxruntime.ai/docs/execution-providers/CUDA-ExecutionProvider.html).
