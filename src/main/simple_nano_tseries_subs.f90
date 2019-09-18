module simple_nano_tseries_subs
use simple_image
use simple_ctf
include 'simple_lib.f08'
implicit none

public :: prep_nano_tseries_subs, calc_neigh_weights, subtr_graphene_backgr
private

integer, parameter :: NNN = 8
integer, parameter :: NITERS = 3

type(ctf)            :: tfun
type(ctfparams)      :: ctfvars
logical, allocatable :: corr_mask(:,:,:), corr_mask_inv(:,:,:)
real                 :: neigh_weights(NNN)
logical              :: do_flip = .false.
integer(kind=kind(ENUM_WCRIT)), parameter :: WCRIT = RANK_INV_CRIT

contains

    subroutine prep_nano_tseries_subs( ldim, mskrad, ctfvars_in )
        integer,         intent(in) :: ldim(3)
        real,            intent(in) :: mskrad
        type(ctfparams), intent(in) :: ctfvars_in
        type(image) :: img_msk
        integer :: i
        ! create real-space correlation masks
        call img_msk%new(ldim, ctfvars%smpd)
        img_msk = 1.0
        call img_msk%mask(mskrad, 'hard')
        corr_mask = img_msk%bin2logical()
        allocate(corr_mask_inv(ldim(1),ldim(2),ldim(3)))
        where( corr_mask )
            corr_mask_inv = .false.
        elsewhere
            corr_mask_inv = .true.
        endwhere
        call img_msk%kill
        ! create CTF object
        ctfvars = ctfvars_in
        tfun    = ctf(ctfvars%smpd, ctfvars%kv, ctfvars%cs, ctfvars%fraca)
        ! set flip flag
        do_flip = ctfvars%ctfflag.eq.CTFFLAG_FLIP
    end subroutine prep_nano_tseries_subs

    subroutine calc_neigh_weights( particle, neighs, neigh_mask )
        class(image), intent(inout) :: particle, neighs(NNN)
        logical,      intent(in)    :: neigh_mask(NNN)
        integer :: i, j
        real    :: neigh_corrs(NNN), sxx
        ! CTF correction by phase-flipping
        if( do_flip )then
            call particle%fft
            call tfun%apply_serial(particle, 'flip', ctfvars)
            call particle%ifft
            do i=1,NNN
                if( neigh_mask(i) )then
                    call neighs(i)%fft
                    call tfun%apply_serial(neighs(i), 'flip', ctfvars)
                    call neighs(i)%ifft
                endif
            end do
        endif
        ! calculate background correlations
        call particle%prenorm4real_corr(sxx, corr_mask_inv)
        neigh_corrs = -1.
        do i=1,NNN
            if( neigh_mask(i) ) neigh_corrs(i) = particle%real_corr_prenorm(neighs(i), sxx, corr_mask_inv)
        end do
        ! calculate background weights
        neigh_weights = corrs2weights(neigh_corrs, WCRIT)
    end subroutine calc_neigh_weights

    subroutine subtr_graphene_backgr( particle, neighs, corrected )
        class(image), intent(inout) :: particle, neighs(NNN), corrected





    end subroutine subtr_graphene_backgr

end module simple_nano_tseries_subs
