# Gabor Features for Class-Average Rejection

Date: 2026-07-14  
Status: proposal

## Opportunity

`model_cavgs_rejection` uses hard validity gates followed by a compact feature
bank and pairwise logistic regression. Cyril's new `image%gabor_filter2D`
routine may provide additional evidence about structured molecular detail.

The aim is not to add a large filter bank to the classifier. It is to derive a
few scalar features describing whether Gabor response is spatially associated
with the particle foreground.

Class averages have no 3D pose assignments. `Filter direction` refers only to
a direction in the image plane, and invariance refers to in-plane rotation of
the class average.

## Current Gabor Output

The routine combines even and odd filters into a phase-insensitive response,
keeps the maximum across sampled filter directions, normalizes each frequency
by its per-image median, and pools the frequencies into one response image.

This image shows where directional structure is strong relative to the rest of
the class average. It is not an absolute measure of texture strength and does
not retain separate filter-direction or frequency responses.

## Proposed First Features

Use the existing cleaned foreground mask to derive three scalars:

1. `log_gabor_fg_bg`: robust foreground response divided by robust background
   response. This measures whether structure is localized to the particle
   rather than diffuse noise or ice.

2. `gabor_support_enrichment`: enrichment of high-response pixels inside the
   foreground, corrected for foreground area. This tests agreement with the
   particle mask rather than a carbon edge or box-spanning contamination.

3. `gabor_fg_coverage`: fraction of foreground pixels above a robust
   background-derived threshold. This distinguishes distributed internal
   structure from a smooth fuzzy ball or a single isolated line.

These features must be interpreted jointly. Strong Gabor response alone is
not necessarily good; crystalline ice and carbon edges can also respond
strongly.

## Initial Policy

- Start with a conservative fixed bank, tentatively `40, 30, 20 A` for chunk
  data. Evaluate pool data separately before using finer scales.
- Remove the edge average, preserve physical sampling, and filter before
  applying any hard particle mask. Use the mask only to summarize the response.
- Add a separate `gabor` feature family, with zero coefficients in all existing
  built-in models.
- Do not use Gabor evidence as a hard reject.
- Do not expose development settings as user-facing parameters.

The bank should avoid frequencies that many classes cannot resolve. Since each
frequency is median-normalized before pooling, an unresolved band can otherwise
contribute normalized noise with the same weight as a useful band.

## Logistic Model

Adding three features expands the possible pairwise terms from 91 to 136. The
first fitted model should therefore use linear Gabor terms or a compact set of
motivated interactions, not every possible pair.

Useful interactions to test include Gabor evidence with nominal resolution,
foreground area, class population, presence, and the existing 100-40 A
center/edge feature. All terms should remain regularized and visible in the
model report.

## Evaluation

First, add the three features to analysis output without changing
classification. Examine per-dataset AUC, correlation with existing features,
current false positives and false negatives, dependence on population and
nominal resolution, ranked class-average stacks, and extraction time.

Then compare:

1. the existing logistic baseline;
2. baseline plus individual Gabor features;
3. baseline plus all three linear Gabor features;
4. baseline plus a small set of curated interactions.

Use identical hard-gate survivors in every comparison. Hold out complete
projects or specimens rather than random class-average rows, and evaluate chunk
and pool models separately.

Promote the feature family only if it improves held-out projects, protects
manual-good recall, improves a recognizable failure mode such as fuzzy balls or
ice, remains useful across specimen types, adds information beyond existing
features, and has acceptable extraction cost.

## Possible Later Extension

If the pooled response is useful but does not separate ice or carbon edges, the
Gabor API could expose per-filter-direction energies for in-plane-rotation-
invariant directional entropy or dominant-direction measurements. Per-frequency
outputs could support coarse-to-fine response ratios.

These richer descriptors should be considered only after the compact feature
family demonstrates held-out value.

## Recommended First Step

Add the three pooled-response features as diagnostics with zero model
coefficients and run them across the existing manually selected chunk datasets.
No production behavior should change until those results are reviewed.
