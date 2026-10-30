# Half-map local enhancement

## Contract

`simple_exec prg=enhance` consumes two independent, unprocessed 3D
half-maps on the same grid and writes four interpretation-only maps:

- the enhanced merged map requested by `outvol`;
- `<outvol>_confidence`, a [0,1] local half-map support weight;
- `<outvol>_locres`, the local FSC=0.143 resolution in Angstroms; and
- `<outvol>_envelope`, a conservative soft molecular envelope.

The input half-maps remain the sole resolution authority. None of the outputs
may be reused for FSC or resolution estimation.

## Implementation

The commander owns command-line validation, image I/O and product naming. The
`halfmap_enhancer` volume-domain object owns the numerical operation:

1. estimate the merged-map background and a conservative noise floor from the
   outer shell using median/MAD statistics, taking the larger robust scale of
   the merged background and half-map difference so correlated solvent is not
   classified as molecular support;
2. visit overlapping cubic windows, tapering each with a separable sine window;
3. regularize the local FSC, convert it to the merged-map Wiener weight
   `2 FSC / (1 + FSC)`, and stop support at the first FSC=0.143 crossing;
4. fit a bounded local Guinier B-factor only within supported shells and fold it
   into a single radial transfer function whose total sharpening gain is capped;
5. overlap-add the filtered windows with the squared taper as normalization;
6. derive a soft envelope from local signal power in excess of the independent
   half-map noise and suppress unsupported background.

The default physical window is approximately 48 Angstroms, clamped to an even
16--64 pixel cube and a 50% stride. The domain constructor accepts explicit
window and stride values for deterministic small tests; they are not user
controls in the first version.

## Review findings

This is complementary to `postprocess_nu`, not a replacement. It is a
standalone local spectral estimator for arbitrary raw half-maps, uses continuous
local FSC/Wiener weights and per-window B factors, and composes tapered windows.
`postprocess_nu` remains the refinement-aware route: its discrete filter bank,
cross-half evidence and Potts-regularized voxel labels use SIMPLE project state
and a global B factor.

The confidence product is intentionally conservative. On the fixed-seed noisy
phantom its molecular-core mean was 0.0050 and its solvent mean was 0.0001;
the enhanced-map correlation with the truth increased from 0.3897 to 0.8062.

## Validation

- `python3 scripts/check_descr.py .`: passed.
- `python3 scripts/check_test_registry.py . --verbose`: passed.
- `python3 scripts/check_flex_dag.py . --verbose`: passed.
- `git diff --check`: passed.
- Incremental Release build: passed.
- `simple_test_exec test=unit_reconstruction suite=half_map_enhancement`:
  17/17 assertions passed, including correlated-background regression.
- `simple_test_exec test=unit_ui suite=ui_visibility`: 357/357 assertions
  passed.
- Complete `unit_reconstruction`: 18 suites and 504/504 assertions passed.
- Complete fast gate: 15/15 entries passed in 40.49 seconds.
- `make -C build -j8 install`: passed, including complete UI JSON validation.

Linux/BOX validation was not run on this macOS host.

## EMDB smoke validation

The command was run directly on the archived half-map pairs for EMD-4997
(232 cubed at 1.075 A/pixel) and EMD-7099 (384 cubed at 0.817 A/pixel).  An
independent SIMPLE automasked half-map FSC gave FSC=0.143 resolutions of 4.02 A
and 3.87 A, respectively.  The deposited primary maps were used only as
external comparison maps, never as inputs to the enhancer.

EMD-4997 was the clean positive case.  In voxels with confidence at least 0.2,
correlation to the deposited map increased from 0.701 for the raw half-map
average to 0.802 after enhancement.  Outside the conservative envelope, the
standard deviation fell to 0.0022 of the raw value.  The high-confidence local
resolution median was 4.83 A.

EMD-7099 exposed two robustness issues during this validation.  Its half-maps
have strongly correlated outer-box background (outer-shell correlation about
0.78), which initially filled the inferred envelope.  The noise model now
estimates excess correlated background from the robust merged-map and
half-difference scales.  Separately applied B-factor and cap transforms also
produced unstable local amplification; these are now one radial transfer
function capped at an eightfold gain.  With both corrections, high-confidence
correlation to the deposited map was retained (0.894 raw, 0.891 enhanced),
while outside-envelope standard deviation fell to 0.032 of the raw value.  The
high-confidence local-resolution median was 4.43 A.  The broader 58.7% envelope
and correlated input background make this a qualified rather than clean
positive case.

The first file-level run also caught and fixed sibling-output naming: the
extension returned by `fname2ext` does not include a leading period.
