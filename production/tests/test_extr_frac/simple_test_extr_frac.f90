program simple_test_extr_frac
use simple_defs
implicit none

integer, parameter :: NPTCLS=50000, MAXITS=40
real,    parameter :: UPDATE_FRAC = 0.1, INITHRESH = 0.8
integer :: iextr_lim, i
real    :: anneal_ratio, extr_thresh

iextr_lim = ceiling(2. * log(real(NPTCLS)) * (2. - UPDATE_FRAC))
do i=1,MAXITS
    anneal_ratio      = max(0., cos(PI / 2. * real(i - 1) / real(iextr_lim))) ! cosine decay
    extr_thresh       = INITHRESH * anneal_ratio                              ! fraction of particles
    print *, i, extr_thresh
end do
end program simple_test_extr_frac
