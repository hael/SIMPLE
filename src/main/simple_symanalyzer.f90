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
integer            :: nsym ! tot # symmetries considered

type sym_stats
    character(len=:), allocatable :: str
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

    subroutine symmetry_tester( vol_in, msk, hp, lp, cn_stop, platonic )
        class(projector), intent(inout) :: vol_in
        real,             intent(in)    :: msk, hp, lp
        integer,          intent(in)    :: cn_stop
        logical,          intent(in)    :: platonic
        type(sym_stats), allocatable    :: pgrps(:)
        real,            allocatable    :: zscores(:), ccs(:), diffs_peak(:), diffs_backgr(:)
        integer,         allocatable    :: peaks(:)
        type(sym)             :: symobj
        character(len=STDLEN) :: errmsg
        integer :: ncsym, icsym, cnt, idsym, ldim(3), isym, fnr, nbackgr
        real    :: smpd, peakval, backgrval, x
        ! get info from vol_in
        smpd  = vol_in%get_smpd()
        ldim  = vol_in%get_ldim()
        ! count # symmetries
        ncsym = cn_stop
        nsym  = ncsym * 2               ! because we always search for dihedral symmetries
        nsym  = nsym  - 1               ! because d1 is omitted
        if( platonic ) nsym = nsym  + 3 ! because of the Platonics (t,o,i)
        if( platonic )then
            if( cn_stop < 5 )then
                errmsg = 'cn range must include rotational symmetries from orders up to 5 when platonic=yes. '//&
                    &'Set cn_stop >= 5 on command line; eval_point_groups'
                THROW_HARD(trim(errmsg))
            endif
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
        allocate(ccs(nsym), peaks(nsym), diffs_peak(nsym), diffs_backgr(nsym))
        ! by definition for c1
        pgrps(1)%cc = 1.0
        ccs(1)      = 1.0
        call eval_point_groups(vol_in, msk, hp, lp, pgrps)
        ! fetch data
        ccs(:)  = pgrps(:)%cc
        ! calculate Z-scores
        zscores = z_scores(ccs)
        ! identify peaks
        peakval   = maxval(zscores(2:))
        nbackgr   = count(zscores(2:) > TINY) - 1 ! -1 because peak position is omitted
        backgrval = 0.
        if( nbackgr > 0 ) backgrval = (sum(zscores(2:), mask=zscores(2:) > TINY) - peakval) / real(nbackgr)
        where( zscores > TINY )
            diffs_peak = sqrt((zscores - peakval)**2.0)
        elsewhere
            diffs_peak = huge(x)
        endwhere
        where( zscores > TINY ) diffs_backgr = sqrt((zscores - backgrval)**2.0)
        where(diffs_peak < diffs_backgr)
            peaks = 1
        elsewhere
            peaks = 0
        endwhere
        ! output
        call fopen(fnr, status='replace', file='symmetry_test_output.txt', action='write')
        write(fnr,'(a)') '>>> RESULTS RANKED ACCORDING TO DEGREE OF SYMMETRY'
        do isym=1,nsym
            write(fnr,'(a,f5.2,a,f5.2,a,i1)') 'POINT-GROUP: '//trim(pgrps(isym)%str)//' CORRELATION: ',&
            &pgrps(isym)%cc, ' Z-SCORE: ', zscores(isym), ' PEAK: ', peaks(isym)
        end do
        write(fnr,'(a)') ''
        call fclose(fnr)
    end subroutine symmetry_tester

    subroutine eval_point_groups( vol_in, msk, hp, lp, pgrps )
        use simple_vol_srch
        class(projector), intent(inout) :: vol_in
        real,             intent(in)    :: msk, hp, lp
        type(sym_stats),  intent(inout) :: pgrps(:)
        type(projector) :: vol_pad
        type(ori)       :: symaxis
        type(sym)       :: symobj
        real            :: smpd
        integer         :: ldim(3), boxpd, igrp, ldim_pd(3)
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
        ! loop over point-groups
        do igrp=2,nsym
            if( DEBUG_HERE )then
                write(logfhandle,*) 'gathering info for point-group: ', pgrps(igrp)%str
            else
                call progress(igrp, nsym)
            endif
            ! make point-group object
            call symobj%new(pgrps(igrp)%str)
            ! locate the symmetry axis and retrieve corr
            if( DEBUG_HERE ) write(logfhandle,*) 'searching for the symmetry axis'
            call volpft_symsrch_init(vol_in, pgrps(igrp)%str, hp, lp)
            call volpft_srch4symaxis(symaxis)
            pgrps(igrp)%cc    = symaxis%get('corr')
        end do
        ! destruct
        call vol_pad%kill
        call vol_in%ifft ! return in real-space
        call symaxis%kill
        call symobj%kill
    end subroutine eval_point_groups

end module simple_symanalyzer
