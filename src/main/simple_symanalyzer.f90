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

logical, parameter :: DEBUG_HERE        = .false.
real,    parameter :: SHSRCH_HWDTH      = 4.0
real,    parameter :: SCORE_PEAK_BOUND  = 0.9
real,    parameter :: ZSCORE_PEAK_BOUND = 1.0
integer            :: nsym       ! total number of symmetries considered
integer            :: kfromto(2) ! Fourier index range

type sym_stats
    character(len=:), allocatable :: str
    real,             allocatable :: fsc(:)
    real :: cc, score
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
        integer           :: isym, ldim(3), boxpd, ldim_pd(3)
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

    subroutine symmetry_tester( vol_in, msk, hp, lp, cn_stop, platonic, pgrp_best )
        class(projector),              intent(inout) :: vol_in
        real,                          intent(in)    :: msk, hp, lp
        integer,                       intent(in)    :: cn_stop
        logical,                       intent(in)    :: platonic
        character(len=:), allocatable, intent(out)   :: pgrp_best
        type(sym_stats), allocatable    :: pgrps(:)
        real,    allocatable  :: scores(:), res(:), zscores(:)
        real,    allocatable  :: scores_peak(:), scores_backgr(:)
        logical, allocatable  :: peak_msk(:)
        type(sym)             :: symobj
        character(len=STDLEN) :: errmsg
        integer :: ncsym, icsym, cnt, idsym, j, ldim(3), nsub, nsub_max, peak_flag
        integer :: isym, jsym, filtsz, iisym, fnr, npeaks, isym_most_likely
        real    :: smpd, kstwo_stat, prob_null
        ! get info from vol_in
        res        = vol_in%get_res()
        smpd       = vol_in%get_smpd()
        filtsz     = vol_in%get_filtsz()
        ldim       = vol_in%get_ldim()
        ! set Fourier index range
        kfromto(1) = calc_fourier_index(hp, ldim(1), smpd)
        kfromto(2) = calc_fourier_index(lp, ldim(1), smpd)
        ! count # symmetries
        ncsym      = cn_stop
        nsym       = ncsym * 2          ! because we always search for dihedral symmetries
        nsym       = nsym  - 1          ! because d1 is omitted
        if( platonic ) nsym = nsym  + 3 ! because of the Platonics (t,o,i)
        if( cn_stop < 10 )then
            errmsg = 'cn range must include rotational symmetries from orders up to 10. '//&
                &'Set cn_stop >= 10 on command line; eval_point_groups'
            THROW_HARD(trim(errmsg))
        endif
        if( platonic )then
            write(logfhandle,'(a,1x,i2,1x,a,1x,i2,1x,a)') '>>> TESTING ROTATIONAL SYMMETRIES FROM', 1, 'TO', cn_stop, ' & DIHEDRAL & PLATONIC GROUPS'
        else
            write(logfhandle,'(a,1x,i2,1x,a,1x,i2,1x,a)') '>>> TESTING ROTATIONAL SYMMETRIES FROM', 1, 'TO', cn_stop, ' & DIHEDRAL GROUPS'
        endif
        ! prepare point-group stats object
        allocate(pgrps(nsym))
        ! C-symmetries
        cnt = 0
        do icsym=1,cn_stop
            cnt            = cnt + 1
            pgrps(cnt)%str = 'c'//int2str(icsym)
        end do
        ! D-symmetries
        do idsym=2,cn_stop
            cnt            = cnt + 1
            pgrps(cnt)%str = 'd'//int2str(idsym)
        end do
        if( platonic )then
            pgrps(cnt + 1)%str = 't'
            pgrps(cnt + 2)%str = 'o'
            pgrps(cnt + 3)%str = 'i'
        endif
        ! gather stats
        allocate(scores(nsym))
        ! by definition for c1
        pgrps(1)%cc    = 1.0
        pgrps(1)%score = 1.0
        scores(1)      = 1.0
        call eval_point_groups(vol_in, msk, hp, lp, pgrps)
        ! fetch data
        scores(:)     = pgrps(:)%score
        ! calculate Z-scores
        zscores       = z_scores(scores)
        ! extract peak and background scores
        peak_msk      = zscores >= ZSCORE_PEAK_BOUND .and. scores >= SCORE_PEAK_BOUND
        npeaks        = count(peak_msk)
        if( npeaks == 0 )then
            THROW_WARN('no symmetry could be identified, npeaks == 0')
            pgrp_best = 'c1'
            return
        endif
        scores_peak   = pack(scores, mask =       peak_msk)
        scores_backgr = pack(scores, mask = .not. peak_msk)
        ! calculate Kolmogorov-Smirnov stats
        call kstwo(scores_peak, npeaks, scores_backgr, nsym - npeaks, kstwo_stat, prob_null)
        ! prob_null represents the significance level for the null hypothesis that the two data sets are drawn from the same distribution, i.e. small prob_null values show that the cumulative distribution functions of the two data sets differ significantly
        ! identify most likely point-group symmetry as highest order one among the peaks
        nsub_max = 0
        do isym=1,nsym
            if( peak_msk(isym) )then
                call symobj%new(pgrps(isym)%str)
                nsub = symobj%get_nsubgrp()
                if( nsub > nsub_max )then
                    isym_most_likely = isym
                    nsub_max = nsub
                endif
            endif
        end do
        ! output
        call fopen(fnr, status='replace', file='symmetry_test_output.txt', action='write')
        write(fnr,'(a)') '>>> RESULTS RANKED ACCORDING TO DEGREE OF SYMMETRY'
        do isym=1,nsym
            if( peak_msk(isym) )then
                peak_flag = 1
            else
                peak_flag = 0
            endif
            write(fnr,'(a,1x,a5,1x,a,f5.2,1x,a,1x,f5.2,1x,a,1x,f5.2,1x,a,1x,i1)') 'POINT-GROUP:',&
                &pgrps(isym)%str, 'SCORE:', pgrps(isym)%score, 'CORRELATION:', pgrps(isym)%cc,&
                'Z-SCORE:', zscores(isym), 'PEAK:', peak_flag
        end do
        write(fnr,'(a)') ''
        write(fnr,'(a)') '>>> MOST LIKELY POINT-GROUP DEFINED AS HIGHEST GROUP AMONG PEAKS'
        pgrp_best = pgrps(isym_most_likely)%str
        write(fnr,'(a,1x,a5,1x,a,f5.2,1x,a,1x,f5.2,1x,a,1x,f5.2)') 'POINT-GROUP:',&
            &pgrps(isym_most_likely)%str, 'SCORE:', pgrps(isym_most_likely)%score, 'CORRELATION:',&
            pgrps(isym_most_likely)%cc, 'Z-SCORE:', zscores(isym_most_likely)
        write(fnr,'(a)') ''
        write(fnr,'(a)') 'KOLMOGOROV-SMIRNOV TEST OF PEAK VS. NON-PEAK DISTRIBUTION'
        write(fnr,'(a)') 'P represents the significance level for the null hypothesis that the two sets are drawn from the same distribution'
        write(fnr,'(a)') 'A small P shows that the cumulative distribution functions of the peak vs. non-peak sets differ significantly'
        write(fnr,'(a)') 'A high K-S value, where K-S .in. [0,1] indicates the same'
        write(fnr,'(a,f5.2)') 'P   = ', prob_null
        write(fnr,'(a,f5.2)') 'K-S = ', kstwo_stat
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
        integer         :: filtsz, ldim(3), boxpd, igrp, ldim_pd(3)
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
        do igrp=2,nsym
            if( DEBUG_HERE )then
                write(logfhandle,*) 'gathering info for point-group: ', pgrps(igrp)%str
            else
                call progress(igrp, nsym)
            endif
            ! make point-group object
            call symobj%new(pgrps(igrp)%str)
            ! locate the symmetry axis
            if( DEBUG_HERE ) write(logfhandle,*) 'searching for the symmetry axis'
            call find_symaxis(pgrps(igrp)%str)
            ! rotate input (non-symmetrized) volume to symmetry axis
            if( DEBUG_HERE ) write(logfhandle,*) 'rotating input volume to symmetry axis'
            call rotvol_slim(vol_pad, rovol_pad, vol_asym_aligned2axis, symaxis)
            call vol_asym_aligned2axis%write('vol_c1_aligned2_'//trim(pgrps(igrp)%str)//'axis.mrc')
            call vol_asym_aligned2axis%mask(msk, 'soft')
            ! generate symmetrized volume
            if( DEBUG_HERE ) write(logfhandle,*) 'generating symmetrized volume'
            call symaverage
            call vol_sym%write('vol_sym_'//trim(pgrps(igrp)%str)//'.mrc')
            call vol_sym%mask(msk, 'soft')
            ! correct for any small discrepancy in shift between the volumes
            if( DEBUG_HERE ) write(logfhandle,*) 'correcting for any small discrepancy in shift between the volumes'
            call vol_asym_aligned2axis%fft
            call vol_sym%fft
            call vol_srch_init(vol_asym_aligned2axis, vol_sym, hp, lp, SHSRCH_HWDTH)
            cxyz = vol_shsrch_minimize()
            ! read back in unmasked volume and shift it before re-applying the mask
            if( DEBUG_HERE ) write(logfhandle,*) 'shifting volume'
            call vol_sym%read('vol_sym_'//trim(pgrps(igrp)%str)//'.mrc')
            call vol_sym%shift(cxyz(2:4))
            call vol_sym%write('vol_sym_'//trim(pgrps(igrp)%str)//'.mrc')
            call vol_sym%mask(msk, 'soft')
            call vol_sym%fft
            ! calculate a correlation coefficient
            if( DEBUG_HERE ) write(logfhandle,*) 'calculating correlation'
            pgrps(igrp)%cc = vol_sym%corr(vol_asym_aligned2axis, lp_dyn=lp, hp_dyn=hp)
            ! calculate FSC
            if( DEBUG_HERE ) write(logfhandle,*) 'calculating FSC'
            if( allocated(pgrps(igrp)%fsc) ) deallocate(pgrps(igrp)%fsc)
            allocate(pgrps(igrp)%fsc(filtsz), source=0.)
            call vol_sym%fsc(vol_asym_aligned2axis, pgrps(igrp)%fsc)
            ! set score (median of FSC in resolution interval)
            pgrps(igrp)%score = max(0.,median(pgrps(igrp)%fsc(kfromto(1):kfromto(2))))
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
            integer           :: isym, nsym_local
            type(ori)         :: o
            real              :: rmat(3,3)
            ! extract the rotation matrices for the symops
            nsym_local = symobj%get_nsym()
            allocate(sym_rmats(nsym_local,3,3))
            do isym=1,nsym_local
                o = symobj%get_symori(isym)
                sym_rmats(isym,:,:) = o%get_mat()
            end do
            ! rotate over symmetry related rotations and update vol_sym
            vol_sym = 0.
            do isym=1,nsym_local
                rmat = matmul(sym_rmats(isym,:,:), rmat_symaxis)
                call o%set_euler(m2euler(rmat))
                call rotvol_slim(vol_pad, rovol_pad, rovol, o)
                call vol_sym%add_workshare(rovol)
            end do
            call vol_sym%div(real(nsym_local))
        end subroutine symaverage

    end subroutine eval_point_groups

end module simple_symanalyzer
