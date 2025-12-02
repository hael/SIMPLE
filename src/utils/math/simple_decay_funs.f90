module simple_decay_funs
include 'simple_lib.f08'
implicit none

public :: calc_update_frac, calc_update_frac_dyn, nsampl_decay, inv_nsampl_decay, calc_nsampl_fromto
public :: cos_decay, inv_cos_decay, extremal_decay2D, extremal_decay
private
#include "simple_local_flags.inc"

contains

    function calc_update_frac( nptcls, nstates, nsample_minmax ) result( update_frac )
        integer, intent(in) :: nptcls, nstates, nsample_minmax(2)
        real    :: update_frac
        integer :: nsampl, nsample_minmax_here(2)
        nsample_minmax_here    = nsample_minmax * nstates
        nsample_minmax_here(1) = min(nptcls,nsample_minmax_here(1))
        nsample_minmax_here(2) = min(nptcls,nsample_minmax_here(2))
        nsampl       = min(nsample_minmax_here(2), nint(0.5 * real(nptcls)))
        nsampl       = max(nsampl, nsample_minmax_here(1))
        nsampl       = min(nptcls, max(nsampl,nsample_minmax_here(2)))
        update_frac  = real(nsampl) / real(nptcls)
        update_frac  = min(1.0, update_frac)
    end function calc_update_frac

    function calc_update_frac_dyn( nptcls, nstates, nsample_minmax, it, maxits ) result( update_frac )
        integer, intent(in) :: nptcls, nstates, nsample_minmax(2), it, maxits
        real    :: update_frac
        integer :: nsampl, nsample_minmax_here(2)
        nsample_minmax_here    = nsample_minmax * nstates
        nsample_minmax_here(1) = min(nptcls,nsample_minmax_here(1))
        nsample_minmax_here(2) = min(nptcls,nsample_minmax_here(2))
        nsampl      = inv_nsampl_decay(it, maxits, nptcls, nsample_minmax_here)
        update_frac = real(nsampl) / real(nptcls)
        update_frac = min(1.0, update_frac)
    end function calc_update_frac_dyn

    function nsampl_decay( it, maxits, nptcls, nsample_minmax ) result( nsampl )
        integer, intent(in) :: it, maxits, nptcls, nsample_minmax(2)
        integer :: nsampl, nsampl_fromto(2)
        nsampl_fromto = calc_nsampl_fromto(nptcls, nsample_minmax)
        nsampl = nint(cos_decay(min(it,maxits), maxits, real(nsampl_fromto)))
    end function nsampl_decay

    function inv_nsampl_decay( it, maxits, nptcls, nsample_minmax ) result( nsampl )
        integer, intent(in) :: it, maxits, nptcls, nsample_minmax(2)
        integer :: nsampl, nsampl_fromto(2)
        nsampl_fromto = calc_nsampl_fromto(nptcls, nsample_minmax)
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

    function calc_nsampl_fromto( nptcls, nsample_minmax ) result( nsampl_fromto )
        integer, intent(in)  :: nptcls, nsample_minmax(2)
        integer :: nsampl_fromto(2)
        if( nint(0.5 * real(nptcls)) < nsample_minmax(1) )then
            nsampl_fromto(1) = nptcls/4
            nsampl_fromto(2) = nptcls
        else
            ! upper limit has half of the particles or NSAMPLE_MAX if overshoot
            nsampl_fromto(2) = nint(0.5 * real(nptcls))
            nsampl_fromto(2) = min(nsampl_fromto(2), nsample_minmax(2))
            nsampl_fromto(1) = nint(real(nsampl_fromto(2)) / 20.)
            nsampl_fromto(1) = max(nsampl_fromto(1), nsample_minmax(1))
        endif
    end function calc_nsampl_fromto

    pure real function extremal_decay2D( extr_iter, extr_lim )
        integer, intent(in) :: extr_iter, extr_lim
        real :: power
        ! factorial decay, -2 because first step is always greedy
        power = real(extr_iter)*MAX_EXTRLIM2D/real(extr_lim) - 2.0
        extremal_decay2D = SNHC2D_INITFRAC * (1.-SNHC2D_DECAY)**power
        extremal_decay2D = min(SNHC2D_INITFRAC, max(0.,extremal_decay2D))
    end function extremal_decay2D

    ! is cosine decay with updated bounds
    pure real function extremal_decay( it, maxits )
        integer, intent(in) :: it, maxits
        if( it <=  2 )then
            ! because first iteration is greedy
            extremal_decay = 0.5
        else if( it <= maxits )then
            ! cosine
            extremal_decay = 0.5**2 * (1.0 + cos(real(it)*PI / real(maxits)))
        else
            ! off
            extremal_decay = 0.0
        endif
    end function extremal_decay

end module simple_decay_funs







