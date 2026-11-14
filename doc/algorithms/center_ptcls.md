# Centring Particles at Scale

## Problem

A picked particle is off centre by the picking error: on the real sets measured here, 8 to 10 Å at the
median and 16 to 22 Å at the 90th percentile. Every alignment searches that error again until the
particle is centred. The polar shift map searches a disk of three prior widths for a particle not yet
centred, and a window about its stored shift once it is centred. The pose students tolerate about 5 Å:
the 3D student's cone holds its answer for 96% of particles at 6.4 Å of shift, and for 28% at 12.8 Å.

Today only a full solve2D centres particles, at about 30 iterations of reading every particle.
`center_ptcls` does it in one pass over the data, plus a small reference step.

## Algorithm

1. **The working grid.** The references' band `lp` (20 Å) sets the box by the library's rule for a stage at
   that band (`lpstages_setlims`): `smpd_target = max(2 Å, lp/3)`, a magic box of at least 88 pixels, never
   above the native box. On bgal and trpv1 it is the 88-pixel floor.
2. **The references.** Draw `nsample` (10,000) active particles at random, all of them when there are fewer.
   On that sample, refine2D makes `ncls` (10) class averages from random references. It uses the
   probabilistic dense table under the polar shift map, at `lp`, centring the class averages, for at most
   15 iterations. A rotational average is not used: it smears elongated views and centred 19% of trpv1's
   picks within 3 Å, against the class averages' 45%.
3. **The pass.** Every active particle is searched once against the references:
   - at the finer band `lpstop` (12 Å), every in-plane rotation, the dense table under the polar shift map;
   - as the first iteration of its own run, so every particle starts not yet centred and searches the disk
     of three `shift_sigma` (10 Å, so 30 Å) coarse-to-fine;
   - the posterior's shift width is measured as it goes.
4. **The record.** Only the shift, its width (`shsig`) and the move (`shincr`) go into `ptcl2D` and `ptcl3D`.
   The class and in-plane angle of the pass mean nothing downstream, and every other parameter is left as it
   was. A particle whose drawn class carried no measured width (an unrefined candidate) is not moved and
   keeps the sentinel, so a later run searches its disk.

A later refine2D or refine3D run then treats the particles as centred. Its dynamic centring searches the
window about each stored shift and takes the disk only for the particles left not centred. `refine3D_local`
keeps that state from the workflow before.

## What it achieves

The residuals below are against each set's converged full solve2D, after their common offset (a set of class
averages defines the centre up to a translation). Two converged solve2D runs agree within 3 Å for 69%
(trpv1) and 76% (bgal) of the same particles.

| set | method | median | within 3 Å | within 5 Å | within 6.4 Å | time |
| --- | --- | --- | --- | --- | --- | --- |
| trpv1, 35,645 picks | raw | 8.0 Å | 0.09 | 0.24 | 0.36 | |
| | classical centring (rotational average) | 5.7 Å | 0.18 | 0.42 | 0.57 | |
| | `center_ptcls` (prototype) | 3.3 Å | 0.45 | 0.67 | 0.76 | 38 s |
| | `center_ptcls` (the program, another sample) | 3.5 Å | 0.43 | | | 40 s |
| | full solve2D | | | | | 468 s |
| bgal, 5,513 picks | raw | 10.1 Å | 0.07 | 0.16 | 0.25 | |
| | `center_ptcls` | 5.6 to 7.7 Å | 0.22 to 0.25 | 0.45 | 0.55 | 11 s |
| | one pass against the full run's 50 classes (ceiling) | 2.1 Å | 0.57 | 0.66 | 0.70 | |

On trpv1, two thirds of the picks end within 5 Å in a twelfth of solve2D's time. bgal is a small set. Its
low-resolution references, not the pass, limit it, and its worse-matching half carries a tail (90th
percentile about 20 Å) that no method moves.

## Rationale

- **Ten references, not five:** 45% within 3 Å on trpv1, against 33%.
- **The finer pass band:** the pass at 12 Å against references at 20 Å centred more particles than either
  at a single band. References made at 15 Å were noisier and did worse.
- **Things that changed nothing:** a wider prior (13 Å), and forcing 15 reference iterations instead of
  refine2D's own stop.
- **The pass's iteration number:** a run marks particles not yet centred only at its first iteration. A pass
  numbered after the references' iterations took the stored width, zero for raw picks, and searched only
  the window: 38% against 45% with the same references.

## Implementation

The program is `center_ptcls`, in `simple_exec`. It is a controller: `commander_center_ptcls` and
`exec_center_ptcls` (`src/main/commanders/simple/simple_commanders_center_ptcls.f90`) run refine2D twice on
work copies of the project. The pass runs in `center_pass/`. The shifts are copied back by the owner in one
write of `ptcl2D` and `ptcl3D`.

The keys are existing ones with this program's defaults: `nsample`, `ncls`, `lp`, `lpstop`, `shift_sigma`.
`nparts` runs both refine2D runs distributed. The workflow test is `center_ptcls`
(`doc/policies/test_environment_policy.md`); the design and the measurements are in
`doc/implementation_notes/completed/center_ptcls.md`.
