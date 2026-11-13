# Polar students: directions beyond the supervised assignment

Status, 2026-11-12: ideas, not planned. Base: `f2e907c24`. The supervised 2D student is
`planned/class_student_2D.md`; the 3D student is `doc/algorithms/pose_student.md`. Nothing here is to be
built before a plan of its own. Each idea lists what it would give, what it needs, the pitfalls, and the
first measurement that would decide whether to pursue it.

## 1. Rotation-invariant embeddings from the HHCNN

**The observation.** The network is rotation-equivariant along the angle. Its convolutions are circular on
the half ring with the conjugate wrap, and its gate and ring pooling act pointwise. The readout takes the
angular Fourier coefficients F_ch,band(m) of its map (orders m = 0 to 16, eight channels by eight ring
bands). An in-plane rotation by α multiplies F(m) by exp(-i m α). Invariant summaries follow:
- **power invariants** |F(m)|^2: about 1,100 numbers, cheap, but blind to the relative phase between
  orders;
- **bispectrum-type invariants** F(m1) F(m2) F*(m1 + m2) over chosen triples: complete up to rotation, and
  without bias from Gaussian noise in expectation. The cryo-EM literature uses bispectrum invariants with
  steerable PCA for class averaging (ASPIRE); here they would be taken on learned feature maps;
- **the class posterior** of a trained student: invariant by construction (the readout maximises over
  rotations), but limited to the teacher's classes.

**Three ways to obtain one.**
1. *From a trained student* (comes almost free with the 2D student). Uses: novelty detection in the
   stream (a confident-nowhere particle far from every template is a view the pool lacks, and a group of
   them is the signal to add classes); splitting a class whose members' embeddings are bimodal; merging
   classes whose member clouds overlap; a per-particle fit score beside `cavg_quality`.
2. *Fixed invariants without training.* Bispectrum invariants of the polar Fourier rings, compressed by
   PCA. Use: class proposals for millions of particles by mini-batch k-means or a nearest-neighbour graph,
   to seed solve2D instead of random references. The cost is linear: at the 3D student's 670 particles
   per second, a million particles in about 25 minutes on 24 threads.
3. *Self-supervised training of the HHCNN* (contrastive). Two augmented copies of one particle (rotation,
   shift, noise, another defocus's CTF) map to the same invariant vector, and different particles apart.
   Rotation invariance comes from the architecture; shift and CTF invariance are learned. The data sets of
   millions of particles supply the volume this needs.

**Pitfalls.**
- **The CTF.** Untrained invariants of CTF-multiplied data cluster by defocus. Option 2 needs phase
  flipping or Wiener filtering first; options 1 and 3 learn around it.
- **Shift.** The readout invariants are not shift-invariant, since the network sees complex Fourier values
  whose phases move with a shift. A modulus front end is shift-invariant but loses phase information.
- **Mirror.** 2D searches no mirrors, so invariants must keep the sign of m (power invariants are fine;
  bispectrum triples must not be symmetrised).

**First measurement.** On the fixed 2D teachers of the supervised plan: do clusters in embedding space
reproduce the teacher's classes (adjusted Rand index), and does the distance to the assigned template
predict which held-out assignments are wrong?

## 2. A shift-robust front end

The supervised student registers shifts after prediction and trains with a measured shift augmentation.
If the 2D student cannot learn under the picking error of raw particles, there are three alternatives:
- a centring step before the student (a cheap translational alignment to a rotationally averaged
  reference);
- a front end of Fourier moduli, shift-invariant but phase-blind;
- a learned equivariance to shifts, folded into option 3 above.

**First measurement.** The 3D student's shift check (`SHIFT_CHECK_MAGS`) at the picking-error percentiles
Phase 0 of the supervised plan records.

## 3. Generative models of the polar transform

A network that predicts the clean polar Fourier transform of a particle at its read-out pose could serve
two purposes:
- **denoising for alignment:** match against a denoised particle and keep the noise model for the
  likelihood;
- **shift-invariant registration by design:** align raw particles to their own generated clean version.

It needs:
- a decoder, which the HHCNN has not, and a loss that respects the Euclidean objective's noise model;
- a guard against the model's own reference bias, for example training on one half and applying it to the
  other.

**First measurement.** Whether polar alignment of denoised particles against the same references beats raw
alignment, on truth-pose simulated data at the SNR of the fixed sets.

## 4. The angle cap in the polar machinery

The student caps its network's half ring at twice the band's angular sampling limit (2 π k_max). bgal
trained twice as fast with its reads inside the seed spread, and the student lost accuracy at 1 π k_max.
The polar calculator still sizes its half ring from the mask (`magic_pftsz`), which oversamples a
low-pass band 2 to 3 times at the native box. On cropped solve3D stages it lands near 1 π k_max, where the
student failed. This was raised with Cyril on 2026-11-11.

**First measurement.** A refine3D stage at a fixed crop with the calculator's half ring at 1, 1.5 and
2 π k_max: the search's accuracy against truth poses and its time.

## 5. The continuous machinery after prediction

After the 3D student, the fastest route to a refined map is the continuous pose refinement
(`refine3D_local refine=cont`) from the student's poses and registered shifts, with no polar search. It
does not yet match the polar machinery's quality; the maintainer is working on it. The same holds in 2D
for continuous in-plane refinement after the 2D student's assignment.

**First measurement.** The student's poses and shifts into `refine3D_local refine=cont`, against the same
poses into `refine3D_local refine=prob_neigh`, at equal iterations on bgal.
