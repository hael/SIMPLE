# Motion Correction in SIMPLE

## Overview

SIMPLE performs movie motion correction in two stages: a global translational alignment stage followed by a local deformation-modeling stage. The algorithm first reads the movie frames from the detector stack and preprocesses them by optional EER fractionation, gain correction, defective-pixel suppression, optional dose-based truncation of the frame series, and optional Fourier downscaling before alignment begins. By default, the algorithm uses the central frame as the alignment reference.

## Global Alignment

In the first stage, SIMPLE estimates a single global translational trajectory for the full field of view. It does so with a hybrid iterative registration scheme that combines a coarse semi-discrete correlation search with a continuous Fourier-space refinement.

During the discrete phase, the algorithm shifts frames according to the current translation estimates, forms a weighted reference average, and re-registers each frame to update the translations. As the iterations proceed, it tightens the alignment resolution from an initial low-pass limit to a final low-pass limit and updates frame weights from the per-frame correlations.

The algorithm then refines the global translations by continuous per-frame optimization in Fourier space against leave-one-out references. After this stage, SIMPLE applies the estimated global shifts to all frames before proceeding to local motion estimation.

## Local Patch Alignment

In the second stage, SIMPLE estimates residual nonuniform motion from a regular grid of image patches. It chooses the patch grid automatically from the physical micrograph dimensions using a target patch size of about 200 A while enforcing a minimum patch size in pixels after any movie scaling.

SIMPLE then aligns each patch stack independently with the same hybrid translational alignment engine used in the global stage, thereby obtaining a local shift trajectory for each patch across frames.

## Deformation Model

SIMPLE regularizes the patch trajectories by fitting a smooth deformation model separately for the x- and y-displacement components. Let \(x\) and \(y\) denote normalized patch-center coordinates and let \(t\) denote the frame index relative to the fixed reference frame. The algorithm expands each displacement component in the 18-term basis

\[
\{t,\ t^2,\ t^3,\ xt,\ xt^2,\ xt^3,\ x^2t,\ x^2t^2,\ x^2t^3,\ yt,\ yt^2,\ yt^3,\ y^2t,\ y^2t^2,\ y^2t^3,\ xyt,\ xyt^2,\ xyt^3\}.
\]

Accordingly, for each component \(u_d(x,y,t)\), with \(d \in \{x,y\}\),

\[
u_d(x,y,t) = \sum_{k=1}^{18} c_{d,k}\,\phi_k(x,y,t),
\]

where the coefficients \(c_{d,k}\) come from a least-squares fit to the measured patch trajectories.

## Patch-Refine Variant

In the optional `patch_refine` mode, SIMPLE inserts an additional model-guided refinement step. It first performs a robust trimmed fit to the initially estimated patch trajectories, then evaluates the fitted deformation field at each patch center to generate predicted local trajectories, and finally uses those predicted trajectories to initialize a further patch-level refinement before fitting the final deformation model.

## Model Acceptance and Fallback

SIMPLE evaluates the adequacy of the fitted deformation model by computing the root-mean-square deviation between fitted and measured patch shifts separately along x and y. If the fit exceeds the acceptance threshold, the algorithm reduces the number of patches and retries the fit. If the fit remains unsatisfactory, it discards the local model and falls back to the global translational correction.

## Final Frame Integration

When the local model passes the acceptance test, SIMPLE uses the fitted deformation field to warp the full movie frames and generate both an unweighted integrated image for CTF estimation and a frame-weighted corrected micrograph. If dose weighting is enabled, it applies dose weighting before computing the final corrected sum.

## Compact Summary

SIMPLE first estimates whole-frame translational motion by iterative frame-to-average registration using a hybrid discrete/continuous correlation optimizer. It then partitions the image into patches, estimates local translational trajectories for each patch, and fits those trajectories with a smooth polynomial deformation field defined over image coordinates and frame time. The algorithm accepts the local solution only if the fitted model reproduces the measured patch trajectories within a prescribed tolerance. Otherwise, it reduces model complexity and ultimately falls back to the global solution. When accepted, the fitted deformation field warps the original movie frames and produces both the corrected micrograph and the companion sum for CTF estimation.
