# Record conventions in detail

Reference for `simple-main-ori`.

## Shifts

`x, y` are pixels of the particle's stack (`os_stk` `box`, `smpd`) in the particle frame; an image
is centred by shifting it by minus `[x, y]`. Cropped work multiplies them by `box_crop/box` and
never writes the scaled value back; resampling a stack multiplies them (`mul_shifts`, as
`scale_project` does). STAR holds `x*smpd` in Angstrom (skill `simple-main-star`).

## Calls that reallocate (row indices go stale)

`new` (kills first), `reallocate` (grows only), `append`, `delete`, `compress`, `merge` (also kills
its 2nd argument), `symmetrize` (repeats rows, no operator).

## The south pole under c1

`sym%rot_to_asym` keeps a record no operator brings into the asymmetric unit: under c1 with e2 exactly
180 degrees (which `m2euler` returns within 0.026 degrees of the south pole) the record stays whole
(fixed 2026-11-05; it was overwritten by an uninitialised local, `symmetry` suite).
