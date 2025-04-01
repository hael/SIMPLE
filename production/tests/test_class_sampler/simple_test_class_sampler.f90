program simple_test_class_sampler
include 'simple_lib.f08'
implicit none
integer, parameter :: NSCORES = 200, NSAMPLES = 75
real    :: scores(NSCORES), score, probs(NSCORES)
integer :: i, inds(NSCORES), sample(NSAMPLES)
logical :: l_sampled(NSCORES)
score = 0.9
do i = 1, NSCORES
    scores(i) = score - ran3() * 0.5
    inds(i)   = i
end do
probs     = scores / sum(scores)
sample    = nmultinomal_sampling(probs, NSAMPLES)
l_sampled = .false.
do i = 1, NSAMPLES
    l_sampled(sample(i)) = .true.
end do
print *, '#       SAMPLES:      ', count(l_sampled)
print *, 'AVG     SAMPLE SCORE: ', sum(scores, mask=l_sampled) / count(l_sampled)
print *, 'AVG NON-SAMPLE SCORE: ', sum(scores, mask=.not. l_sampled) / count(.not. l_sampled)
end program simple_test_class_sampler