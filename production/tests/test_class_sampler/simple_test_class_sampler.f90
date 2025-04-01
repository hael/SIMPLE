program simple_test_class_sampler
include 'simple_lib.f08'
implicit none
integer, parameter :: NSCORES = 200, NSAMPLES = 75, NTIMES = 1000
real    :: scores(NSCORES), score, probs(NSCORES), counts(NSCORES)
integer :: i, inds(NSCORES), sample(NSAMPLES), smpl_inds(NTIMES)
logical :: l_sampled(NSCORES)
score = 0.9
do i = 1, NSCORES
    scores(i) = score - ran3() * 0.5
    inds(i)   = i
end do
probs     = scores / sum(scores)
smpl_inds = nmultinomal(probs, NTIMES)
counts    = 0.
do i = 1, NTIMES
    counts(smpl_inds(i)) = counts(smpl_inds(i)) + 1.
enddo
call hpsort(counts, inds)
call reverse(inds)
l_sampled = .false.
do i = 1, NSAMPLES
    l_sampled(inds(i)) = .true.
end do
print *, '#       SAMPLES:      ', count(l_sampled)
print *, 'AVG     SAMPLE SCORE: ', sum(scores, mask=l_sampled) / count(l_sampled)
print *, 'AVG NON-SAMPLE SCORE: ', sum(scores, mask=.not. l_sampled) / count(.not. l_sampled)
end program simple_test_class_sampler
