# Continuous 3-D pose refinement: consolidated scientific review

**Status:** COMPLETED REVIEW WITH SUPERSEDED RECOMMENDATIONS IDENTIFIED
**Consolidated:** 2026-08-25
**History and handoff:** [continuous_3D_pose_end_polishing_history_and_handoff.md](continuous_3D_pose_end_polishing_history_and_handoff.md)
**Validation evidence:** [continuous_3D_pose_end_polishing_validation_evidence.md](continuous_3D_pose_end_polishing_validation_evidence.md)

## Purpose

This note preserves the durable scientific conclusions of the former literature review while separating them from recommendations that were later superseded by operator diagnostics and Hans's clarification.

## Durable conclusions

### Pose refinement is a local, non-convex inference problem

Continuous angular refinement is scientifically plausible and supported in the literature, but local minima and outliers are expected. A lower image residual does not establish a correct pose. The observed 43.1-degree wrong endpoint after a 99.7-percent objective reduction is consistent with this general non-convexity rather than evidence of floating-point overflow.

### The statistical objective matters

The Cartesian numerical fixture must use the intended particle-specific CTF and shell-variance weighting. “Euclidean” does not mean unweighted. Whitening by $1/\sqrt{\sigma_s^2}$ changes both residual and Jacobian contributions by frequency shell.

The existing PFTC matcher objective cannot be inserted unchanged into the Cartesian LM. Polar and Cartesian Fourier sums use different integration measures. A future PFTC derivative would need a separate derivation with the executed radial/Jacobian weighting, normalization, CTF, sigma, and interpolation order.

### Forward and adjoint operators must match the chosen model

For any iterative reconstruction model $A$, its backprojection must implement the corresponding adjoint $A^H$ closely enough for the stated normal equations and adjoint tests. This general inverse-problem principle remains valid.

It does not follow that every difference between a simulator and a production operator is a production defect. First decide which discrete forward model is authoritative, then validate its forward/adjoint pair and use a separate generator to test generalization.

### Finite boxes and real-space windows change Fourier data

Cropping or masking in real space is convolution in Fourier space and can create frequency correlations. The `REALSPACE_CLIP_NATIVE_FFT` diagnostic correctly localized a difference between the independent simulator path and the declared PCG Cartesian model.

The later conclusion is narrower than the former review: the difference proves that the two generators are not the same discrete operator. It does not by itself require adding the simulator's finite-box transform to production reconstruction or pose evaluation.

### Reference-derived pose bias is scientifically plausible

Pose estimation against a reconstructed reference can follow reference error or fitted noise, especially at weak high frequencies. Gold-standard half maps reduce correlated overfitting but do not make every local pose estimate unbiased. Cross-fitted references, independent particles, frequency/SNR weighting, and pose priors are useful diagnostic or future design tools if drift remains after the numerical operator is validated.

### Priors and frequency control can improve pose accuracy

Local polishing starts from an informative pose. A Gaussian or trust-region prior can encode that information, and frequency/SNR weighting can prevent weak high-frequency terms from dominating. These are later estimator-policy choices, not fixes to apply before the present LM implementation passes its independent mathematical oracles.

### Joint pose and volume optimization is deferred

Alternating or joint reconstruction and pose refinement is supported by published work, but it creates a feedback loop in which pose and map errors can reinforce each other. It should not be used to hide an unresolved fixed-reference forward, derivative, or statistical mismatch.

## Superseded recommendation

The former review recommended changing production PCG to the simulator-matched operator

$$
A=F_N C F_P^{-1}G_P
$$

and implementing its adjoint before interpreting reconstructed-reference drift.

That was a reasonable hypothesis after the finite-box diagnostic, but it is not the current production recommendation. Later source and operator-contract review established that the simulator's crop/native-FFT step is outside the declared PCG model. Adding it only because it makes one simulator agree would change the estimator rather than repair a confirmed defect.

The retained lesson is to use the finite-box path as one explicit independent-oracle variant, label its model boundary, and avoid claiming equivalence to the PCG operator.

## Recommended scientific sequence now

1. Independently validate the current shell-weighted Cartesian objective.
2. Independently validate all five Jacobian columns, $J^H J$, $J^H r$, damping, scaling, and one LM step.
3. Validate CTF and sigma edge cases and exact rollback.
4. Compare the optimized evaluator with a structurally separate forward oracle away from stencil boundaries.
5. Repeat over multiple boxes, shell ranges, and asymmetric volumes.
6. Only after those gates pass, write a separate contract for a standalone continuous `refine3D` validation mode.
7. Use that standalone route to decide whether frequency weighting, pose priors, cross-fitting, or another acceptance signal is needed.

## Literature retained for future design

1. Sorzano, C. O. S., et al. (2022). “On bias, variance, overfitting, gold standard and consensus in single-particle analysis by cryo-electron microscopy.” *Acta Crystallographica D* 78, 410--423. https://doi.org/10.1107/S2059798322001978
2. Scheres, S. H. W. (2012). “RELION: Implementation of a Bayesian approach to cryo-EM structure determination.” *Journal of Structural Biology* 180, 519--530. https://doi.org/10.1016/j.jsb.2012.09.006
3. Wang, L., Shkolnisky, Y., and Singer, A. (2013). “A Fourier-based approach for iterative 3D reconstruction from cryo-EM images.” arXiv:1307.5824. https://arxiv.org/abs/1307.5824
4. Zeng, G. L., and Gullberg, G. T. (2000). “Unmatched projector/backprojector pairs in an iterative reconstruction algorithm.” *IEEE Transactions on Medical Imaging* 19, 548--555. https://doi.org/10.1109/42.870265
5. Stewart, A., and Grigorieff, N. (2004). “Noise bias in the refinement of structures derived from single particles.” *Ultramicroscopy* 102, 67--84. https://doi.org/10.1016/j.ultramic.2004.08.008
6. Shaikh, T. R., Hegerl, R., and Frank, J. (2003). “An approach to examining model dependence in EM reconstructions using cross-validation.” *Journal of Structural Biology* 142, 301--310. https://doi.org/10.1016/S1047-8477(03)00029-7
7. Scheres, S. H. W., and Chen, S. (2012). “Prevention of overfitting in cryo-EM structure determination.” *Nature Methods* 9, 853--854. https://doi.org/10.1038/nmeth.2115
8. Scheres, S. H. W. (2012). “A Bayesian view on cryo-EM structure determination.” *Journal of Molecular Biology* 415, 406--418. https://doi.org/10.1016/j.jmb.2011.11.010
9. Diepeveen, W., Lellmann, J., Öktem, O., and Schönlieb, C.-B. (2023). “Regularizing orientation estimation in cryogenic electron microscopy three-dimensional map refinement through measure-based lifting over Riemannian manifolds.” *SIAM Journal on Imaging Sciences* 16, 1440--1490. https://doi.org/10.1137/22M1520773
10. Ramlaul, K., Palmer, C. M., and Nakane, T. (2020). “Mitigating local over-fitting during single particle reconstruction with SIDESPLITTER.” *Journal of Structural Biology* 211, 107545. https://doi.org/10.1016/j.jsb.2020.107545
11. Ortiz, S., et al. (2020). “Validation tests for cryo-EM maps using an independent particle set.” *Journal of Structural Biology: X* 4, 100032. https://doi.org/10.1016/j.yjsbx.2020.100032
12. Zehni, M., et al. (2020). “Joint angular refinement and reconstruction for single-particle cryo-EM.” *IEEE Transactions on Image Processing* 29, 6151--6163. https://doi.org/10.1109/TIP.2020.2984313

The references support methodology and risk analysis. They do not prove that the current SIMPLE implementation realizes the same objectives or numerical operators.
