program simple_test_multinomal
use simple_rnd
implicit none
real    :: pvec(7), cnts(7)
integer :: which, i
pvec(7) = 10.
pvec(6) = 20.
pvec(5) = 30.
pvec(4) = 100.
pvec(3) = 80.
pvec(2) = 20.
pvec(1) = 10.

pvec = pvec / sum(pvec)

print *, pvec

! sample the distribution and calculate frequencies
call seed_rnd
cnts = 0.
do i=1,10000
    which = multinomal( pvec )
    cnts(which) = cnts(which) + 1.0
end do
pvec = cnts / sum(cnts)

print *, pvec

end program simple_test_multinomal
