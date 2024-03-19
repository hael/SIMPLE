program simple_test_extr_frac
include 'simple_lib.f08'
implicit none

integer, parameter :: NPTCLS    = 2000000 ! two million
integer, parameter :: NSAMPL_LB = 50000   ! minimum 50k,  default: min(NSAMPL_LB, nint(0.5 * real(NPTCLS)))
integer, parameter :: NSAMPL_UB = 200000  ! maximum 200k, default: min(NSAMPL_UB, NPTCLS) 
integer, parameter :: MAXITS    = 40      ! upper iteration bound
real    :: update_frac
integer :: nsampl, i

do i = 1,MAXITS
    update_frac = inv_frac_decay(i, MAXITS, [NSAMPL_LB,NSAMPL_UB], NPTCLS)
    nsampl      = nint(update_frac * real(NPTCLS))
    print *, inv_nsampl_decay(i, MAXITS, [NSAMPL_LB,NSAMPL_UB]), nsampl
end do

contains

    function inv_frac_decay( i, maxits, nsampl_fromto, nptcls ) result( frac )
        integer, intent(in) :: i, maxits, nsampl_fromto(2), nptcls
        real :: frac
        frac = cos_decay(maxits - i, maxits, real(nsampl_fromto)) / real(nptcls)
    end function inv_frac_decay

    function inv_nsampl_decay( i, maxits, nsampl_fromto ) result( nsampl )
        integer, intent(in) :: i, maxits, nsampl_fromto(2)
        integer :: nsampl
        nsampl = nint(cos_decay(maxits - i, maxits, real(nsampl_fromto)))
    end function inv_nsampl_decay

    function cos_decay( i, maxits, eps_fromto ) result( eps )
        integer, intent(in) :: i, maxits
        real,    intent(in) :: eps_fromto(2)
        real :: delta, eps
        delta = 0.5 * (eps_fromto(2) - eps_fromto(1))
        eps   = eps_fromto(1) + delta + delta * cos((real(i) * PI) / real(MAXITS))
    end function cos_decay

end program simple_test_extr_frac
