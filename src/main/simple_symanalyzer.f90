module simple_symanalyzer
include 'simple_lib.f08'
use simple_volpft_symsrch
use simple_image,          only: image
use simple_projector,      only: projector
use simple_projector_hlev, only: rotvol_slim, rotvol
use simple_sym,            only: sym
use simple_ori,            only: ori, m2euler
implicit none

public :: symmetrize_map, symmetry_tester, print_subgroups
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG_HERE         = .false.
logical, parameter :: WRITE_VOLUMES      = .false.
real,    parameter :: SHSRCH_HWDTH       = 4.0
real,    parameter :: ZSCORE_LBOUND_PEAK = 1.0
real,    parameter :: ZSCORE_LBOUND_C1   = 2.9
integer            :: nsym       ! total number of symmetries considered
integer            :: kfromto(2) ! Fourier index range

type sym_stats
    character(len=:), allocatable :: str
    real :: cc
end type sym_stats

contains

    subroutine symmetrize_map( vol_in, params, vol_out )
        use simple_parameters, only: parameters
        class(projector),  intent(inout) :: vol_in
        class(parameters), intent(inout) :: params
        class(image),      intent(inout) :: vol_out
        type(ori)         :: symaxis, o
        type(sym)         :: symobj
        type(image)       :: rovol_pad, rovol, vol_asym_aligned2axis
        type(projector)   :: vol_pad
        real, allocatable :: sym_rmats(:,:,:)
        integer           :: isym, ldim(3), boxpd, ldim_pd(3)
        real              :: rmat_symaxis(3,3), rmat(3,3), smpd
        ! make point-group object
        call symobj%new(params%pgrp)
        nsym = symobj%get_nsym()
        ! extract the rotation matrices for the symops
        allocate(sym_rmats(nsym,3,3))
        do isym=1,nsym
            call symobj%get_symori(isym, o)
            sym_rmats(isym,:,:) = o%get_mat()
        end do
        ! init search object
        call volpft_symsrch_init(vol_in, params%pgrp, params%hp, params%lp)
        ! identify the symmetry axis
        call volpft_srch4symaxis(symaxis)
        ! get the rotation matrix for the symaxis
        rmat_symaxis = symaxis%get_mat()
        ! prepare for volume rotation
        ldim = [params%box,params%box,params%box]
        call vol_in%new(ldim, params%smpd)
        call vol_in%read(params%vols(1))
        call vol_in%shift([params%xsh,params%ysh,params%zsh])
        ldim    = vol_in%get_ldim()
        smpd    = vol_in%get_smpd()
        boxpd   = 2 * round2even(KBALPHA * real(ldim(1) / 2))
        ldim_pd = [boxpd,boxpd,boxpd]
        call vol_in%ifft
        ! rotate asymmetric volume
        call o%set_euler(m2euler(rmat_symaxis))
        ! rotate over symmetry related rotations and update average
        call vol_out%new(ldim, smpd)
        call vol_asym_aligned2axis%new(ldim, smpd)
        call rovol%new(ldim, smpd)
        call rovol_pad%new(ldim_pd, smpd)
        call vol_pad%new(ldim_pd, smpd)
        call vol_in%pad(vol_pad)
        call vol_pad%fft
        call vol_pad%expand_cmat(KBALPHA)
        ! generate volume aligned to symaxis
        call rotvol_slim(vol_pad, rovol_pad, vol_asym_aligned2axis, symaxis)
        call vol_asym_aligned2axis%write('vol_c1_aligned2_'//trim(params%pgrp)//'axis.mrc')
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
        call vol_asym_aligned2axis%kill
        call o%kill
        call symaxis%kill
        call symobj%kill
    end subroutine symmetrize_map

    subroutine print_subgroups
        type(sym_stats), allocatable :: pgrps(:)
        integer, parameter :: cn_stop = 10
        integer   :: ncsym, nsym, cnt, icsym, idsym, isym
        type(sym) :: symobj
        ! count # symmetries
        ncsym = cn_stop
        nsym  = ncsym * 2 ! because we always search for dihedral symmetries
        nsym  = nsym  - 1 ! because d1 is omitted
        nsym  = nsym  + 3 ! because of the Platonics (t,o,i)
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
        ! platonic
        pgrps(cnt + 1)%str = 't'
        pgrps(cnt + 2)%str = 'o'
        pgrps(cnt + 3)%str = 'i'
        do isym=1,nsym
            ! make point-group object
            call symobj%new(pgrps(isym)%str)
            write(logfhandle,'(a)') 'Subroups of point-group '//trim(pgrps(isym)%str)//':'
            ! print subgroups
            write(logfhandle,'(a)') symobj%get_all_subgrps_descr()
            write(logfhandle,'(a)') ''
        end do
        call symobj%kill
    end subroutine print_subgroups

    subroutine symmetry_tester( vol_in, msk, hp, lp, cn_stop, platonic, pgrp_out )
        class(projector), intent(inout) :: vol_in
        real,             intent(in)    :: msk, hp, lp
        integer,          intent(in)    :: cn_stop
        logical,          intent(in)    :: platonic
        character(len=3), intent(out)   :: pgrp_out
        type(sym_stats), allocatable    :: pgrps(:)
        real,    allocatable  :: ccs(:), zscores(:)
        integer, allocatable  :: peaks(:)
        type(sym)             :: symobj
        character(len=STDLEN) :: errmsg
        character(len=3)      :: pgrp_sub
        integer :: ncsym, icsym, cnt, idsym, j, ldim(3), nsub, isub
        integer :: isym, jsym, filtsz, iisym, fnr, highest_pgrp_detected
        real    :: smpd
        ! get info from vol_in
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
        if( platonic )then
            if( cn_stop < 5 )then
                errmsg = 'cn range must include rotational symmetries from orders up to 5. '//&
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
        allocate(ccs(nsym), peaks(nsym))
        ! by definition for c1
        pgrps(1)%cc = 1.0
        ccs(1)      = 1.0
        call eval_point_groups(vol_in, msk, hp, lp, pgrps)
        ! fetch data
        ccs(:) = pgrps(:)%cc
        ! calculate Z-scores
        zscores = z_scores(ccs)
        ! set peaks
        ! test for assymetry
        if( zscores(1) >= ZSCORE_LBOUND_C1 )then
            peaks    = 0
            peaks(1) = 1
        else ! identify symmetry peaks
            where( zscores >= ZSCORE_LBOUND_PEAK )
                peaks = 1
            elsewhere
                peaks = 0
            endwhere
        endif
        ! output
        call fopen(fnr, status='replace', file='symmetry_test_output.txt', action='write')
        write(fnr,'(a)') '>>> RESULTS RANKED ACCORDING TO DEGREE OF SYMMETRY'
        do isym=1,nsym
            write(fnr,'(a,f5.2,a,f5.2,a,i1)') 'POINT-GROUP: '//trim(pgrps(isym)%str)//' CORRELATION: ',&
            &ccs(isym), ' Z-SCORE: ', zscores(isym), ' PEAK: ', peaks(isym)
            if( peaks(isym) == 1 ) highest_pgrp_detected = isym
        end do
        pgrp_out = trim(pgrps(highest_pgrp_detected)%str)
        if( highest_pgrp_detected == 1 )then
            write(fnr,'(a)') '>>> NO SYMMETRY DETECTED, SUGGESTED POINT-GROUP: c1'
        else
            write(fnr,'(a)') '>>> SYMMETRY DETECTED, SUGGESTED POINT-GROUP: '//trim(pgrp_out)
            call symobj%new(pgrp_out)
            nsub = symobj%get_nsubgrp()
            cnt  = 0
            do isub=1,nsub
                pgrp_sub = symobj%get_subgrp_descr(isub)
                do isym=1,nsym
                    if( trim(pgrps(isym)%str) .eq. trim(pgrp_sub) ) cnt = cnt + peaks(isym)
                enddo
            end do
            write(fnr,'(a)') int2str(cnt)//'/'//int2str(nsub)//' SUBGROUPS OF THE SUGGESTED POINT-GROUP ALSO DETECTED'
        endif
        call symobj%kill
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
            if( WRITE_VOLUMES ) call vol_asym_aligned2axis%write('vol_c1_aligned2_'//trim(pgrps(igrp)%str)//'axis.mrc')
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
            if( WRITE_VOLUMES )then
                call vol_sym%write('vol_sym_'//trim(pgrps(igrp)%str)//'.mrc')
            else
                call del_file('vol_sym_'//trim(pgrps(igrp)%str)//'.mrc')
            endif
            call vol_sym%mask(msk, 'soft')
            call vol_sym%fft
            ! calculate a correlation coefficient
            if( DEBUG_HERE ) write(logfhandle,*) 'calculating correlation'
            pgrps(igrp)%cc = vol_sym%corr(vol_asym_aligned2axis, lp_dyn=lp, hp_dyn=hp)
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
                call symobj%get_symori(isym, o)
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
            call o%kill
        end subroutine symaverage

    end subroutine eval_point_groups

end module simple_symanalyzer
