module simple_symanalyzer
include 'simple_lib.f08'
use simple_volpft_symsrch
use simple_image,          only: image
use simple_projector,      only: projector
use simple_projector_hlev, only: rotvol_slim
use simple_sym,            only: sym
use simple_ori,            only: ori, m2euler
implicit none

public :: symmetrize_map, symmetry_tester
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG_HERE = .false.

type sym_stats
    character(len=:), allocatable :: str
    real,             allocatable :: fsc(:)
    real :: cc, cc_avg, score
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

    subroutine symmetry_tester( vol_in, msk, hp, lp, cn_start, cn_stop, dihedral, platonic )
        class(projector), intent(inout) :: vol_in
        real,             intent(in)    :: msk, hp, lp
        integer,          intent(in)    :: cn_start, cn_stop
        logical,          intent(in)    :: dihedral, platonic
        type(sym_stats), allocatable    :: pgrps(:)
        logical, allocatable  :: scoring_groups(:)
        integer, allocatable  :: inds(:)
        real,    allocatable  :: scores(:), res(:), zscores(:)
        character(len=3)      :: subgrp
        character(len=STDLEN) :: errmsg
        type(sym) :: symobj
        integer   :: ncsyms, nsyms, icsym, cnt, idsym, nscoring, j, ldim(3), ccn_start
        integer   :: isym, jsym, isub, nsubs, filtsz, iisym, fnr, kfromto(2)
        real      :: smpd
        ! to ensure correct input
        ccn_start  = max(2,cn_start)
        ! get info from vol_in
        res        = vol_in%get_res()
        smpd       = vol_in%get_smpd()
        filtsz     = vol_in%get_filtsz()
        ldim       = vol_in%get_ldim()
        ! set Fourier index range
        kfromto(1) = calc_fourier_index(hp, ldim(1), smpd)
        kfromto(2) = calc_fourier_index(lp, ldim(1), smpd)
        ! count # symmetries
        ncsyms = cn_stop - ccn_start + 1
        if( dihedral .or. platonic )then
            nsyms = ncsyms * 2
        else
            nsyms = ncsyms
        endif
        if( platonic )then
            if( ccn_start > 2 .or. cn_stop < 5 )then
                errmsg = 'cn range must include rotational symmetries from orders 2-5 when searching for Platonic groups. '//&
                    &'Set cn_start = 2 and cn_stop > 5 on command line; eval_point_groups'
                THROW_HARD(trim(errmsg))
            endif
            nsyms = nsyms + 3
        endif
        write(*,'(a,1x,i2,1x,a,1x,i2)') '>>> TESTING C-SYMMETRIES FROM', ccn_start, 'TO', cn_stop
        ! prepare point-group stats object
        allocate(pgrps(nsyms))
        cnt = 0
        do icsym=ccn_start,cn_stop
            cnt            = cnt + 1
            pgrps(cnt)%str = 'c'//int2str(icsym)
        end do
        if( dihedral .or. platonic )then
            do idsym=ccn_start,cn_stop
                cnt            = cnt + 1
                pgrps(cnt)%str = 'd'//int2str(idsym)
            end do
            write(*,'(a)') '>>> TESTING DIHEDRAL SYMMETRIES'
        endif
        if( platonic )then
            pgrps(cnt + 1)%str = 't'
            pgrps(cnt + 2)%str = 'o'
            pgrps(cnt + 3)%str = 'i'
            write(*,'(a)') '>>> TESTING PLATONIC SYMMETRIES'
        endif
        ! gather stats
        call eval_point_groups(vol_in, msk, hp, lp, pgrps)
        allocate(scoring_groups(nsyms), inds(nsyms), scores(nsyms))
        scoring_groups = .false.
        do isym=1,nsyms
            ! make point-group object
            if( DEBUG_HERE ) print *, 'point group considered: ', trim(pgrps(isym)%str)
            call symobj%new(pgrps(isym)%str)
            nsubs = symobj%get_nsubgrp()
            ! identify scoring groups
            ! all subgroups of the group under consideration are part of the scoring group
            scoring_groups = .false.
            do isub=1,nsubs
                subgrp = symobj%get_subgrp_descr(isub)
                do jsym=1,nsyms
                    select case(trim(subgrp))
                        case('c1','C1')
                            if( DEBUG_HERE ) print *, 'c1 case: ', trim(subgrp)
                            cycle
                        case DEFAULT
                            if( trim(subgrp) .eq. trim(pgrps(jsym)%str) )then
                                scoring_groups(jsym) = .true.
                                if( DEBUG_HERE )then
                                    if( scoring_groups(jsym) ) print *, 'scoring group: ', trim(subgrp)
                                endif
                            endif
                    end select
                end do
            end do
            pgrps(isym)%cc_avg = 0.
            do jsym=1,nsyms
                if( scoring_groups(jsym) )then
                    pgrps(isym)%cc_avg  = pgrps(isym)%cc_avg  + pgrps(jsym)%cc
                endif
            end do
            nscoring           = count(scoring_groups)
            pgrps(isym)%cc_avg = pgrps(isym)%cc_avg / real(nscoring)
            pgrps(isym)%score  = max(0.,median(pgrps(isym)%fsc(kfromto(1):kfromto(2))))
            scores(isym)       = pgrps(isym)%score
        end do
        ! calculate Z-scores
        zscores = z_scores(scores)
        ! produce ranked output
        inds = (/(isym,isym=1,nsyms)/)
        call hpsort(scores, inds)
        call reverse(inds)
        do isym=1,nsyms
            iisym = inds(isym)
            write(*,'(a,1x,i2,1x,a,1x,a,1x,a,f5.2,1x,a,1x,f5.2,1x,a,1x,f5.2)') 'RANK', isym, 'POINT-GROUP:',&
                &pgrps(iisym)%str, 'SCORE:', pgrps(iisym)%score, 'CORRELATION:', pgrps(iisym)%cc_avg, 'Z-SCORE:', zscores(iisym)
        end do
        call fopen(fnr, status='replace', file='symmetry_test_fscs.txt', action='write')
        do isym=1,nsyms
            iisym = inds(isym)
            write(fnr,'(a,1x,i2,1x,a,1x,a,1x,a,f5.2,1x,a,1x,f5.2,1x,a,1x,f5.2)') 'RANK', isym, 'POINT-GROUP:',&
                &pgrps(iisym)%str, 'SCORE:', pgrps(iisym)%score, 'CORRELATION:', pgrps(iisym)%cc_avg, 'Z-SCORE:', zscores(iisym)
            do j=1,size(res)
                write(fnr,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', pgrps(iisym)%fsc(j)
            end do
        end do
        call fclose(fnr)
    end subroutine symmetry_tester

    subroutine eval_point_groups( vol_in, msk, hp, lp, pgrps )
        use simple_vol_srch
        class(projector), intent(inout) :: vol_in
        real,             intent(in)    :: msk, hp, lp
        type(sym_stats),  intent(inout) :: pgrps(:)
        type(projector) :: vol_pad
        type(image)     :: rovol_pad, rovol, vol_asym_aligned2axis, vol_sym
        type(ori)       :: symaxis
        type(sym)       :: symobj
        real            :: rmat_symaxis(3,3), smpd, cxyz(4)
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
        filtsz = vol_in%get_filtsz()
        ! intermediate vols
        call rovol%new(ldim, smpd)
        call rovol_pad%new(ldim_pd, smpd)
        ! loop over point-groups
        ngrps = size(pgrps)
        do igrp=1,ngrps
            if( trim(pgrps(igrp)%str) .eq. 'c1' .or. trim(pgrps(igrp)%str) .eq. 'C1' )then
                THROW_HARD('cannot evaluate pgrp=c1, nonsensical; eval_point_groups')
            endif
            call progress(igrp, ngrps)
            ! make point-group object
            call symobj%new(pgrps(igrp)%str)
            ! locate the symmetry axis
            call find_symaxis(pgrps(igrp)%str)
            ! rotate input (non-symmetrized) volume to symmetry axis
            call rotvol_slim(vol_pad, rovol_pad, vol_asym_aligned2axis, symaxis)
            call vol_asym_aligned2axis%write('vol_c1_aligned2_'//trim(pgrps(igrp)%str)//'axis.mrc')
            call vol_asym_aligned2axis%mask(msk, 'soft')
            ! generate symmetrized volume
            call symaverage
            call vol_sym%write('vol_sym_'//trim(pgrps(igrp)%str)//'.mrc')
            call vol_sym%mask(msk, 'soft')
            ! correct for any small discrepancy in shift between the volumes
            call vol_asym_aligned2axis%fft
            call vol_sym%fft
            call vol_srch_init(vol_asym_aligned2axis, vol_sym, hp, lp, 2.0)
            cxyz = vol_shsrch_minimize()
            ! read back in unmasked volume and shift it before re-applying the mask
            call vol_sym%read('vol_sym_'//trim(pgrps(igrp)%str)//'.mrc')
            call vol_sym%shift(cxyz(2:4))
            call vol_sym%write('vol_sym_'//trim(pgrps(igrp)%str)//'.mrc')
            call vol_sym%mask(msk, 'soft')
            call vol_sym%fft
            ! calculate a correlation coefficient
            pgrps(igrp)%cc = vol_sym%corr(vol_asym_aligned2axis, lp_dyn=lp, hp_dyn=hp)
            ! calculate FSC
            if( allocated(pgrps(igrp)%fsc) ) deallocate(pgrps(igrp)%fsc)
            allocate(pgrps(igrp)%fsc(filtsz), source=0.)
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
