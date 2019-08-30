module simple_nano_tseries_subs
use simple_image
use simple_ctf
include 'simple_lib.f08'
implicit none

public :: prep_nano_tseries_subs, subtr_graphene_backgr
private

integer, parameter :: NNN=8
integer, parameter :: NITERS=3
type(image)        :: subtracted(NNN)
type(ctf)          :: tfun
type(ctfparams)    :: ctfvars
logical            :: do_flip = .false.
integer(kind=kind(ENUM_WCRIT)), parameter :: WCRITS(2) = [CORRW_CRIT,RANK_INV_CRIT]

contains

    subroutine prep_nano_tseries_subs( ldim, ctfvars_in )
        integer,         intent(in) :: ldim(3)
        type(ctfparams), intent(in) :: ctfvars_in
        integer :: i
        ! create images for subtraction
        do i=1,NNN
            call subtracted(i)%new(ldim, ctfvars%smpd, wthreads=.false.)
        end do
        ! create CTF object
        ctfvars = ctfvars_in
        tfun    = ctf(ctfvars%smpd, ctfvars%kv, ctfvars%cs, ctfvars%fraca)
        ! set flip flag
        do_flip = ctfvars%ctfflag.eq.CTFFLAG_FLIP
    end subroutine prep_nano_tseries_subs

    subroutine subtr_graphene_backgr( particle, neighs, mask, corrected )
        class(image), intent(inout) :: particle, neighs(NNN), corrected
        logical,      intent(in)    :: mask(NNN)
        integer :: i, j
        real    :: ws(NNN)
        real    :: corrs(NNN)
        ! Fourier transformation
        call particle%fft
        do i=1,NNN
            if( mask(i) ) call neighs(i)%fft
        end do
        ! CTF correction
        if( do_flip )then
            call tfun%apply_serial(particle, 'flip', ctfvars)
            do i=1,NNN
                if( mask(i) ) call tfun%apply_serial(neighs(i), 'flip', ctfvars)
            end do
        endif
        ! create subtracted images and corrected average
        call corrected%zero
        ws = 1. / real(count(mask))
        do i=1,NNN
            if( mask(i) )then
                call subtracted(i)%set(particle)
                call subtracted(i)%subtr(neighs(i))
                call corrected%add(subtracted(i), ws(i))
            endif
        end do
        ! optimise weights
        do j=1,NITERS
            ! (1) calc corrs
            corrs = 0.
            do i=1,NNN
                if( mask(i) ) corrs(i) = corrected%corr_serial(subtracted(i))
            end do
            ! (2) derive weights from corrs
            ws = corrs2weights(corrs, WCRITS(1), WCRITS(2))
            ! (3) update weighted average
            call corrected%zero
            do i=1,NNN
                if( mask(i) ) call corrected%add(subtracted(i), ws(i))
            end do
        end do
        ! return in real-space & normalize
        call corrected%ifft
        call corrected%norm
    end subroutine subtr_graphene_backgr

end module simple_nano_tseries_subs
