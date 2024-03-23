program simple_test_extr_frac
include 'simple_lib.f08'
use simple_dyn_ufrac
implicit none

integer, parameter :: NPTCLS = 2000000 ! two million
integer, parameter :: MAXITS = 40      ! upper iteration bound
real    :: update_frac
integer :: nsampl, i, nsampl_fromto(2)

do i = 1,MAXITS
    nsampl      = inv_nsampl_decay( i, MAXITS, NPTCLS)
    update_frac = real(nsampl) / real(NPTCLS)
    print *, i, nsampl, update_frac
end do

end program simple_test_extr_frac
