# The 2D class student against reference matching

Status, 2026-11-14: measurement record; the 2D student stays open (Hans, 2026-11-15, section 6). Base: `06d3a4de9` (in-plane posterior weighting and sigma2
annealing in solve2D). Runs, scripts and tools: `/data2/agent_runs/refmatch_2D/` (`tools/`, `runs/`). One seed per
set: a lead, not a finding.

## 1. The question

The class student (`doc/algorithms/class_student.md`) is meant to make extending a 2D solution cheaper: it gives
each new particle a shortlist of likely classes, and alignment then tests only the shortlist. The alternative
needs no network: match every new particle against the teacher's class averages in one `refine2D` pass
(reference matching). Which route should extend a solve2D on a large data set?

Rule frozen before the runs: reference matching suffices if it costs no more per particle than the student's
prediction alone, and its class averages are no worse than the student's.

## 2. What was run

On bgal (5,513 particles, 50 classes, 2,756 training) and trpv1 (35,645 particles, 100 classes, 10,000 training),
seed 1, 24 threads, the Dell otherwise idle:

- `validate_class_student` on the current code: the solve2D teacher, the student (width 16), `center_ptcls`, the
  student's prediction, and the full search (one `refine2D refine=prob` iteration at `prob_beta` 100 over all
  classes against the teacher's averages, at the full box, on the centred particles);
- reference matching of the same unseen particles: one `refine2D refine=prob` pass at `prob_beta` 100, the 30 Å
  shift disk (`shift_sigma` 10) and `lp` 15, as the full search, in four arms: from the raw picks or from
  `center_ptcls`'s shifts; at the full box or at the teacher's 120 px working box; at the working box against the
  teacher's iteration-29 averages (its last at 120 px) or its final averages scaled to 120 px;
- class averages of each matching arm's assignment by `make_cavgs` at the full box, as the validation makes the
  student's;
- every route's class averages scored against the teacher's on one common class set (the classes holding at least
  ten particles in every route's assignment): the FRC 0.5 resolution (`tools/frc.py`, which reproduces the
  validation's own values exactly), the crossing interpolated between shells, and the mean FRC out to 15 Å.

## 3. Results

**The validation on the current code:**

| | bgal | trpv1 |
| --- | --- | --- |
| particles given a shortlist | 80.9% | 65.8% |
| shortlists holding the full search's class (goal 95%) | 72.0% | 90.9% |
| same, perfectly centred | 87.7% | 96.0% |
| shortlist length (median) | 4 of 50 | 7 of 100 |
| student's route cost (centring + prediction) | 3.24 ms | 1.65 ms |
| full search cost (full box, centred) | 1.90 ms | 2.00 ms |

**Cost per unseen particle and agreement with the full search's class:**

| route | bgal | trpv1 | same class as the full search |
| --- | --- | --- | --- |
| student, prediction alone | 0.56 ms | 0.41 ms | |
| student, centring + prediction | 3.24 ms | 1.65 ms | |
| matching, raw picks, full box | 1.88 ms | 1.68 ms | 100% / 99.8% |
| matching, raw picks, 120 px, final averages | 0.68 ms | 0.52 ms | 88.2% / 90.9% |
| matching, raw picks, 120 px, iteration 29 | 0.74 ms | 0.58 ms | 85.9% / 87.3% |
| matching, centred, 120 px, iteration 29 | 0.68 ms | 0.60 ms | 85.8% / 87.3% |

**Class averages against the teacher's** (37 common classes on bgal, 53 on trpv1; FRC 0.5 crossing interpolated /
mean FRC to 15 Å):

| route | bgal | trpv1 |
| --- | --- | --- |
| student | 14.41 Å / 0.803 | 18.16 Å / 0.740 |
| full search (centred, full box) | 14.36 Å / 0.802 | 20.09 Å / 0.735 |
| matching, raw picks, full box | 14.24 Å / 0.805 | 21.03 Å / 0.776 |
| matching, raw picks, 120 px, final averages | 14.27 Å / 0.799 | 21.67 Å / 0.773 |

## 4. What it shows

1. **Reference matching needs no centring.** From the raw picks it gives the same classes as from `center_ptcls`'s
   shifts (100% and 99.8% at the full box) and the same class averages: its 30 Å shift disk centres the particle.
2. **The working box costs little.** At 120 px matching is three times cheaper than at the full box; 88 to 91% of
   the particles keep the full box's class, and the class averages stay within this measure's noise.
3. **No route has consistently better class averages.** On bgal all are equal. On trpv1 the two measures
   disagree: the student is better at the FRC 0.5 crossing (18.2 Å against 21.7 Å), matching on the mean FRC to
   15 Å (0.77 against 0.74).
4. **Against the frozen rule:** quality met; cost not met, narrowly. The student's prediction alone is 20 to 25%
   cheaper than matching at the working box. The rule was lenient to the student: its prediction does not finish
   the job (shortlisted particles still need a search over the shortlist, the others a full search), and it
   needs centred particles. On raw picks its route costs 2.5 to 5 times more than matching.
5. **The full search got costlier with the in-plane posterior weighting:** 1.9 to 2.0 ms per particle at the full
   box, against 1.1 ms before it.

## 5. Corrections on the way

- The FRC 0.5 resolution is quantized to whole Fourier shells (256 × 1.275 / 22 = 14.84 Å on bgal), so a median over
  classes hides differences under a shell. The interpolated crossing and the band's mean FRC were added.
- First reads of 17.0 Å (bgal) and 25.9 Å (trpv1) for matching at the working box compared different class sets
  (the population filter depends on the two assignments it is given). On one common set they disappear.
- `make_cavgs` from the command line needs `nparts`; the validation runs it in process, which skips the check.
- `scale` writes into its own execution folder unless given `mkdir=no`.

## 6. Reading

On these two small sets, reference matching at the working box does the whole job (class, in-plane rotation and shift
for every particle) for 0.5 to 0.7 ms per particle, with no training and no centring, at no measurable cost in
quality, and the student's prediction is about 0.1 ms per particle cheaper on centred particles.

**What this does not decide** (Hans, 2026-11-15): costs on sets of 5,500 and 35,600 particles say little about the
data sets of millions the student is meant for, where the fixed costs (the teacher, the training, the centring's
references) are spread over far more particles and the per-particle terms are what remains. Nothing here closes
the 2D student. Its open questions are the shortlists' safety under centring error and what a better teacher and
readout give; the in-plane posterior is planned for both students in
`planned/student_inplane_posterior.md`. Reference matching stays the comparison to report beside the student.

## 7. Relation to NUViT (Gu, Heeter and Zhong, arXiv 2610.32041, October 2026)

NUViT is a vision transformer for Fourier-space data. It resamples the transform with a non-uniform FFT onto a polar
grid over half a turn, as the polar search and our HHCNN already do, and argues the same reason: rings align with
frequency shells and the low frequencies get their share. Where the two differ:

| | NUViT | HHCNN (the students' network) |
| --- | --- | --- |
| size and training | about 25M parameters, 12 transformer layers, hours on an H100, across data sets | three small convolution blocks and a template readout, about a minute on CPU, per data set |
| in-plane rotation | learned from data (absolute positional encoding) | exact on the rotation grid; the readout searches every rotation, which becomes an output |
| mirror | not handled | invariant by construction |
| representation | amplitude and cos/sin of the phase | complex values, real kernels and a modulus gate (commutes with conjugation) |
| CTF and noise | raw coefficients | the CTF cross term, shells whitened by sigma2, the first block normalised by the CTF² envelope |
| shifts | not handled (the authors note encoders are not robust to them) | phase ramps; robustness from training augmentation |
| supervision | hard labels; their baselines overfit (99% train, 1 to 50% test) | soft posterior rows of the probability table, a calibrated held-out threshold |

One idea touches a measured weakness: NUViT's amplitude channel is exactly shift-invariant. The prior evidence is
against amplitudes alone: a pose student on power spectra in place of the complex values reached 0.006 within 6° on
T1000 (0.80 with phases) and 0.05 on bgal (0.50) (`completed/class_student_2D.md`, the probes after Phase 0). An
amplitude channel beside the complex input is a different question, and the pose student is where it could pay, its
full search being far larger, and it is planned there first, in `planned/pose_student_amplitude_channel.md`.
