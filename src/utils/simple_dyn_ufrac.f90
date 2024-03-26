module simple_dyn_ufrac
include 'simple_lib.f08'
implicit none

public :: inv_nsampl_decay
private

contains

    function inv_nsampl_decay( it, maxits, nptcls ) result( nsampl )
        integer, intent(in) :: it, maxits, nptcls
        integer :: nsampl, nsampl_fromto(2)
        nsampl_fromto = calc_nsampl_fromto(nptcls)
        nsampl = nint(inv_cos_decay(min(it,maxits), maxits, real(nsampl_fromto)))
    end function inv_nsampl_decay

    function cos_decay( i, maxits, eps_fromto ) result( eps )
        integer, intent(in) :: i, maxits
        real,    intent(in) :: eps_fromto(2)
        real :: delta, eps
        delta = 0.5 * (eps_fromto(2) - eps_fromto(1))
        eps   = eps_fromto(1) + delta + delta * cos((real(i) * PI) / real(MAXITS))
    end function cos_decay

    function inv_cos_decay( i, maxits, eps_fromto ) result( eps )
        integer, intent(in) :: i, maxits
        real,    intent(in) :: eps_fromto(2)
        real :: delta, eps
        delta = 0.5 * (eps_fromto(2) - eps_fromto(1))
        eps   = eps_fromto(1) + delta + delta * cos((real(maxits - i) * PI) / real(MAXITS))
    end function inv_cos_decay

    function calc_nsampl_fromto( nptcls ) result( nsampl_fromto )
        integer, intent(in)  :: nptcls
        integer :: nsampl_fromto(2)
        if( nptcls < NSAMPL_LB )then
            nsampl_fromto(1) = min(NSAMPL_LB, nint(UPDATE_FRAC_LB_SMALL * real(nptcls)))
            nsampl_fromto(2) = min(NSAMPL_UB, nint(UPDATE_FRAC_UB_SMALL * real(nptcls))) 
        else
            nsampl_fromto(1) = min(NSAMPL_LB, nint(UPDATE_FRAC_LB * real(nptcls)))
            nsampl_fromto(2) = min(NSAMPL_UB, nint(UPDATE_FRAC_UB * real(nptcls))) 
        endif
    end function calc_nsampl_fromto

end module simple_dyn_ufrac







