# Even working grid policy

Status: active. Decided 2026-11-10 after the dimension-parity review of the same date.

## Rule

Every grid SIMPLE transforms is even on each transformed axis: particle boxes (`box`), working
boxes (`box_crop`), padded boxes, class averages, volumes, micrographs and movie frames. The
real-space origin of such a grid is voxel `ldim/2+1`, the pixel `fft` moves to the Fourier phase
origin; masks, Gaussians, shape generators, PDB placement (`center_pdbcoord`), rotation (`rtsq`)
and the corner-to-centre conversion of box files (`corner + box/2`) all use that voxel.

Odd dimensions exist only in real space and only outside the workflows: an image utility may read
an odd or non-square image (`stk=` without `box=`) and binarize, mask or measure it. The first
transform stops it (`shift_phorig`).

## Where the rule stops a run

- `params%new`: an odd `box` or `box_crop` from the command line, a project or a volume header
  stops with the key named; a box read from a plain `stk=` header is left to the transform.
- Project stack adders (`add_stk`, `add_single_stk`, `add_stktab`): an odd box is refused next to the
  square check, so an imported stack or class-average set never enters a project odd.
- Extraction: an odd box from a box file stops with the hint to give an even `box=`.
- `shift_phorig` (every `fft`), `rtsq`, `memoize_mask_coords`, `cartft_calc`, `wave_field`.

Pickers and every internal rescale produce even sizes (`round2even`, `find_larger_magic_box`); the
magic-box table is all even.

## What is not supported, by decision

- Odd working grids: no `fftshift`/`ifftshift` pair, no odd Hermitian storage, no odd branches in
  rotation, resize or Fourier expansion. The odd branches of `rtsq` and `expand_ft` were deleted
  rather than fixed. Do not add odd-parity formulas; add a stop where an odd size can enter.
- Silent resizing of an odd input: padding changes the field of view, cropping discards data,
  Fourier resizing changes the sampling and pose scale. The user chooses the size
  (`prg=print_magic_boxes` lists FFT-efficient ones).

## Consequences for code

- A centre is `ldim/2+1` (1-based) or `real(ldim/2)*smpd` (physical). On an even grid
  `real(ldim)/2.` and `nint(real(ldim)/2.)` agree with it; do not migrate working expressions
  mechanically, and do not introduce `(ldim-1)/2` or `(ldim+1)/2` midpoints for a grid origin.
- Physical extents (`real(box)/2.` as a radius, `2*pi/N` periods, taper widths) are lengths, not
  origins; keep them as they are.
- Same-parity resize (`pad`, `clip`, `cartesian_fourier` embed/crop) is implied by even plus even.
