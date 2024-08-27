module simple_decay_funs
include 'simple_lib.f08'
implicit none

public :: nsampl_decay, inv_nsampl_decay, calc_nsampl_fromto, inv_cos_decay
private

real,    parameter :: UPDATE_FRAC_MIN = 0.1
real,    parameter :: UPDATE_FRAC_MAX = 0.5
integer, parameter :: NSAMPL_MIN      = 10000
integer, parameter :: NSAMPL_MAX      = 200000

contains

    function nsampl_decay( it, maxits, nptcls ) result( nsampl )
        integer, intent(in) :: it, maxits, nptcls
        integer :: nsampl, nsampl_fromto(2)
        nsampl_fromto = calc_nsampl_fromto(nptcls)
        nsampl = nint(cos_decay(min(it,maxits), maxits, real(nsampl_fromto)))
    end function nsampl_decay

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
        if( nint(UPDATE_FRAC_MAX * real(nptcls)) < NSAMPL_MIN )then
            nsampl_fromto(1) = nptcls/4
            nsampl_fromto(2) = nptcls
        else
            nsampl_fromto(1) = nint(UPDATE_FRAC_MIN * real(nptcls))
            nsampl_fromto(1) = max(nsampl_fromto(1), NSAMPL_MIN)
            nsampl_fromto(2) = nint(UPDATE_FRAC_MAX * real(nptcls))
            nsampl_fromto(2) = min(nsampl_fromto(2), NSAMPL_MAX)
        endif
    end function calc_nsampl_fromto

end module simple_decay_funs







