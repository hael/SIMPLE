module simple_symanalyzer
include 'simple_lib.f08'
use simple_volpft_symsrch
use simple_image,          only: image
use simple_projector,      only: projector
use simple_projector_hlev, only: rotvol_slim
use simple_sym,            only: sym
use simple_ori,            only: ori, m2euler
implicit none

public :: symmetrize_map
private

type sym_stats
    character(len=:), allocatable :: str
    real,             allocatable :: fsc(:)
    real :: cc
end type sym_stats

contains

    subroutine symmetrize_map( vol_in, pgrp, hp, lp, vol_out )
        class(projector), intent(inout) :: vol_in
        character(len=*), intent(in)    :: pgrp
        real,             intent(in)    :: hp, lp
        class(image),     intent(inout) :: vol_out
        type(ori)         :: symaxis, o
        type(sym)         :: symobj
        type(image)       :: rovol_pad, rovol
        type(projector)   :: vol_pad
        real, allocatable :: sym_rmats(:,:,:)
        integer           :: isym, nsym, ldim(3), boxpd, ldim_pd(3)
        real              :: rmat_symaxis(3,3), rmat(3,3), smpd
        ! make point-group object
        call symobj%new(pgrp)
        nsym = symobj%get_nsym()
        ! extract the rotation matrices for the symops
        allocate(sym_rmats(nsym,3,3))
        do isym=1,nsym
            o = symobj%get_symori(isym)
            sym_rmats(isym,:,:) = o%get_mat()
        end do
        ! init search object
        call volpft_symsrch_init(vol_in, pgrp, hp, lp)
        ! identify the symmetry axis
        call volpft_srch4symaxis(symaxis)
        ! get the rotation matrix for the symaxis
        rmat_symaxis = symaxis%get_mat()
        ! prepare for volume rotation
        ldim    = vol_in%get_ldim()
        smpd    = vol_in%get_smpd()
        boxpd   = 2 * round2even(KBALPHA * real(ldim(1) / 2))
        ldim_pd = [boxpd,boxpd,boxpd]
        call vol_in%ifft
        call vol_out%new(ldim, smpd)
        call rovol%new(ldim, smpd)
        call rovol_pad%new(ldim_pd, smpd)
        call vol_pad%new(ldim_pd, smpd)
        call vol_in%pad(vol_pad)
        call vol_pad%fft
        call vol_pad%expand_cmat(KBALPHA)
        ! rotate over symmetry related rotations and update average
        do isym=1,nsym
            rmat = matmul(sym_rmats(isym,:,:), rmat_symaxis)
            call o%set_euler(m2euler(rmat))
            call rotvol_slim(vol_pad, rovol_pad, rovol, o)
            call vol_out%add_workshare(rovol)
        end do
        call vol_out%div(real(nsym))
        ! destruct
        call vol_pad%kill_expanded
        call vol_pad%kill
        call rovol%kill
        call rovol_pad%kill
        call o%kill
        call symaxis%kill
        call symobj%kill
    end subroutine symmetrize_map

    subroutine eval_platonic_point_groups( vol_in, hp, lp )
        class(projector), intent(inout) :: vol_in
        real,             intent(in)    :: hp, lp
        type(sym_stats), allocatable :: pgrps(:)
        integer, parameter :: NGRPS = 11
        integer           :: sz
        real              :: cc_t, cc_t_anti, cc_o, cc_o_anti, cc_i, cc_i_anti
        real, allocatable :: fsc_t(:), fsc_t_anti(:), fsc_o(:), fsc_o_anti(:), fsc_i(:), fsc_i_anti(:)
        ! prepare point-group stats object
        allocate(pgrps(NGRPS))
        pgrps(1)%str  = 'c2'
        pgrps(2)%str  = 'c3'
        pgrps(3)%str  = 'c4'
        pgrps(4)%str  = 'c5'
        pgrps(5)%str  = 'd2'
        pgrps(6)%str  = 'd3'
        pgrps(7)%str  = 'd4'
        pgrps(8)%str  = 'd5'
        pgrps(9)%str  = 't'
        pgrps(10)%str = 'o'
        pgrps(11)%str = 'i'
        ! gather stats
        call eval_point_groups(vol_in, hp, lp, pgrps)
        ! calculate cc-based scores and anti-scores
        cc_t      = pgrps(1)%cc + pgrps(2)%cc  + pgrps(5)%cc  + pgrps(9)%cc / 4.
        cc_t_anti = pgrps(3)%cc + pgrps(4)%cc  + pgrps(6)%cc  + pgrps(7)%cc +&
                  & pgrps(8)%cc + pgrps(10)%cc + pgrps(11)%cc / 7.
        cc_o      = pgrps(1)%cc + pgrps(2)%cc  + pgrps(3)%cc  + pgrps(5)%cc +&
                  & pgrps(6)%cc + pgrps(7)%cc  + pgrps(9)%cc  + pgrps(10)%cc / 8.
        cc_o_anti = pgrps(4)%cc + pgrps(8)%cc  + pgrps(11)%cc / 3.
        cc_i      = pgrps(1)%cc + pgrps(2)%cc  + pgrps(4)%cc  + pgrps(5)%cc +&
                  & pgrps(6)%cc + pgrps(8)%cc  + pgrps(9)%cc  + pgrps(11)%cc / 8.
        cc_i_anti = pgrps(3)%cc + pgrps(7)%cc  + pgrps(10)%cc / 3.
        ! calculate fsc-based scores and anti-scores
        sz = size(pgrps(1)%fsc)
        allocate(fsc_t(sz), fsc_t_anti(sz), fsc_o(sz), fsc_o_anti(sz), fsc_i(sz), fsc_i_anti(sz), source=0.)
        fsc_t      = pgrps(1)%fsc + pgrps(2)%fsc  + pgrps(5)%fsc  + pgrps(9)%fsc / 4.
        fsc_t_anti = pgrps(3)%fsc + pgrps(4)%fsc  + pgrps(6)%fsc  + pgrps(7)%fsc +&
                   & pgrps(8)%fsc + pgrps(10)%fsc + pgrps(11)%fsc / 7.
        fsc_o      = pgrps(1)%fsc + pgrps(2)%fsc  + pgrps(3)%fsc  + pgrps(5)%fsc +&
                   & pgrps(6)%fsc + pgrps(7)%fsc  + pgrps(9)%fsc  + pgrps(10)%fsc / 8.
        fsc_o_anti = pgrps(4)%fsc + pgrps(8)%fsc  + pgrps(11)%fsc / 3.
        fsc_i      = pgrps(1)%fsc + pgrps(2)%fsc  + pgrps(4)%fsc  + pgrps(5)%fsc +&
                   & pgrps(6)%fsc + pgrps(8)%fsc  + pgrps(9)%fsc  + pgrps(11)%fsc / 8.
        fsc_i_anti = pgrps(3)%fsc + pgrps(7)%fsc  + pgrps(10)%fsc / 3.
    end subroutine eval_platonic_point_groups

    subroutine eval_c_and_d_point_groups( vol_in, hp, lp, cn_start, cn_stop, dihedral )
        class(projector), intent(inout) :: vol_in
        real,             intent(in)    :: hp, lp
        integer,          intent(in)    :: cn_start, cn_stop
        logical,          intent(in)    :: dihedral
        type(sym_stats), allocatable :: pgrps(:)
        integer :: ncsyms, nsyms, icsym, cnt, idsym
        ! count # symmetries
        ncsyms = cn_stop - cn_start + 1
        if( dihedral )then
            nsyms = ncsyms * 2
        else
            nsyms = ncsyms
        endif
        ! prepare point-group stats object
        allocate(pgrps(nsyms))
        cnt = 0
        do icsym=cn_start,cn_stop
            cnt  = cnt + 1
            pgrps(cnt)%str = 'c'//int2str(icsym)
        end do
        if( dihedral )then
            do idsym=cn_start,cn_stop
                cnt  = cnt + 1
                pgrps(cnt)%str = 'd'//int2str(idsym)
            end do
        endif
        ! gather stats
        call eval_point_groups(vol_in, hp, lp, pgrps)
    end subroutine eval_c_and_d_point_groups

    subroutine eval_point_groups( vol_in, hp, lp, pgrps )
        class(projector), intent(inout) :: vol_in
        real,             intent(in)    :: hp, lp
        type(sym_stats),  intent(inout) :: pgrps(:)
        type(projector) :: vol_pad
        type(image)     :: rovol_pad, rovol, vol_asym_aligned2axis, vol_sym
        type(ori)       :: symaxis
        type(sym)       :: symobj
        real            :: rmat_symaxis(3,3), smpd
        integer         :: filtsz, ldim(3), boxpd, igrp, ldim_pd(3), ngrps
        ! prepare for volume rotations
        ldim    = vol_in%get_ldim()
        smpd    = vol_in%get_smpd()
        boxpd   = 2 * round2even(KBALPHA * real(ldim(1) / 2))
        ldim_pd = [boxpd,boxpd,boxpd]
        ! make padded volume for interpolation
        call vol_pad%new(ldim_pd, smpd)
        call vol_in%pad(vol_pad)
        call vol_pad%fft
        call vol_pad%expand_cmat(KBALPHA)
        ! make outputs
        call vol_sym%new(ldim, smpd)
        call vol_asym_aligned2axis%new(ldim, smpd)
        filtsz = vol_sym%get_filtsz()
        ! intermediate vols
        call rovol%new(ldim, smpd)
        call rovol_pad%new(ldim_pd, smpd)
        ! loop over point-groups
        ngrps = size(pgrps)
        do igrp=1,ngrps
            ! make point-group object
            call symobj%new(pgrps(igrp)%str)
            ! locate the symmetry axis
            call find_symaxis(pgrps(igrp)%str)
            ! rotate input (non-symmetrized) volume to symmetry axis
            call rotvol_slim(vol_pad, rovol_pad, vol_asym_aligned2axis, symaxis)
            call vol_asym_aligned2axis%write('vol_c1_aligned_'//trim(pgrps(igrp)%str)//'.mrc')
            ! generate symmetrized volume
            call symaverage
            call vol_sym%write('vol_sym_'//trim(pgrps(igrp)%str)//'.mrc')
            ! calculate a correlation coefficient
            pgrps(igrp)%cc = vol_sym%corr(vol_asym_aligned2axis, lp_dyn=lp, hp_dyn=hp)
            ! calculate FSC
            if( allocated(pgrps(igrp)%fsc) ) deallocate(pgrps(igrp)%fsc)
            allocate(pgrps(igrp)%fsc(filtsz))
            call vol_sym%fsc(vol_asym_aligned2axis, pgrps(igrp)%fsc)
        end do
        ! destruct
        call vol_pad%kill
        call rovol_pad%kill
        call rovol%kill
        call vol_asym_aligned2axis%kill
        call vol_sym%kill
        call symaxis%kill
        call symobj%kill

        contains

            subroutine find_symaxis( pgrp )
                character(len=*), intent(in) :: pgrp
                call volpft_symsrch_init(vol_in, pgrp, hp, lp)
                call volpft_srch4symaxis(symaxis)
                call vol_in%ifft ! return in real-space
                ! get the rotation matrix for the symaxis
                rmat_symaxis = symaxis%get_mat()
            end subroutine find_symaxis

            subroutine symaverage
                real, allocatable :: sym_rmats(:,:,:)
                integer           :: isym, nsym
                type(ori)         :: o
                real              :: rmat(3,3)
                ! extract the rotation matrices for the symops
                nsym = symobj%get_nsym()
                allocate(sym_rmats(nsym,3,3))
                do isym=1,nsym
                    o = symobj%get_symori(isym)
                    sym_rmats(isym,:,:) = o%get_mat()
                end do
                ! rotate over symmetry related rotations and update vol_sym
                vol_sym = 0.
                do isym=1,nsym
                    rmat = matmul(sym_rmats(isym,:,:), rmat_symaxis)
                    call o%set_euler(m2euler(rmat))
                    call rotvol_slim(vol_pad, rovol_pad, rovol, o)
                    call vol_sym%add_workshare(rovol)
                end do
                call vol_sym%div(real(nsym))
            end subroutine symaverage

    end subroutine eval_point_groups

end module simple_symanalyzer
