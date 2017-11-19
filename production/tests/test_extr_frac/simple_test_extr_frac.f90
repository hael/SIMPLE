program simple_test_extr_frac
use simple_defs
implicit none
integer, parameter :: NPTCLS=50000, MAXITS=40
integer :: iextr_lim, i
real    :: extr_frac
iextr_lim = ceiling(2.*log(real(NPTCLS)))
do i=1,MAXITS
    extr_frac = EXTRINITHRESH * cos(PI/2. * real(i-1)/real(iextr_lim)) ! cosine decay
    extr_frac = min(EXTRINITHRESH, max(0.0, extr_frac))                ! fraction of particles to randomize
    print *, i, extr_frac
end do
end program simple_test_extr_frac