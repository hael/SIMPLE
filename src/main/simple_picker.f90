! particle picker
module simple_picker
!$ use omp_lib
!$ use omp_lib_kinds
use simple_defs
use simple_math,   only: sortmeans, round2even
use simple_image,  only: image
use simple_math,   only: euclid, hpsort   
use simple_syslib, only: alloc_errchk
implicit none

public :: init_picker, exec_picker, kill_picker
private

! PEAK STATS INDICES
integer,          parameter   :: CC2REF   = 1
integer,          parameter   :: SDEV     = 2
integer,          parameter   :: DYNRANGE = 3
integer,          parameter   :: SSCORE   = 4
! OTHER PARAMS
integer,          parameter   :: NSTAT   = 4
integer,          parameter   :: MAXKMIT = 20
real,             parameter   :: BOXFRAC = 0.5
logical,          parameter   :: WRITESHRUNKEN = .false., DOPRINT = .true.
! VARS
type(image)                   :: micrograph, mic_shrunken, mic_shrunken_refine, ptcl_target
type(image),      allocatable :: refs(:), refs_refine(:)
logical,          allocatable :: selected_peak_positions(:)
real,             allocatable :: sxx(:), sxx_refine(:), corrmat(:,:)
integer,          allocatable :: peak_positions(:,:), peak_positions_refined(:,:), refmat(:,:)
character(len=:), allocatable :: micname, refsname
real, allocatable             :: peak_stats(:,:)
character(len=STDLEN)         :: boxname
integer                       :: ldim(3), ldim_refs(3), ldim_refs_refine(3), ldim_shrink(3)
integer                       :: ldim_shrink_refine(3), ntargets, nx, ny, nx_refine, ny_refine
integer                       :: nrefs, npeaks, npeaks_sel, orig_box, lfny, nbackgr
real                          :: smpd, smpd_shrunken, smpd_shrunken_refine, corrmax, corrmin
real                          :: msk, msk_refine, lp, distthr, ndev
logical                       :: rm_outliers = .true.

contains

    subroutine init_picker( micfname, refsfname, smpd_in, lp_in, distthr_in, ndev_in )
        use simple_fileio,   only:  remove_abspath,fname_new_ext
        use simple_imgfile,  only: find_ldim_nptcls
        character(len=*),           intent(in) :: micfname, refsfname
        real,                       intent(in) :: smpd_in
        real,             optional, intent(in) :: lp_in, distthr_in, ndev_in
        type(image) :: refimg
        integer     :: alloc_stat, ifoo, iref
        allocate(micname,  source=trim(micfname), stat=alloc_stat)
        call alloc_errchk('picker;init, 1', alloc_stat)
        allocate(refsname, source=trim(refsfname), stat=alloc_stat)
        call alloc_errchk('picker;init, 2', alloc_stat)
        boxname = remove_abspath( fname_new_ext(micname,'box') )   
        smpd    = smpd_in
        lp      = 20.0
        if( present(lp_in) ) lp = lp_in
        ndev    = 2.0
        if( present(ndev_in)) ndev = ndev_in
        ! read micrograph
        call find_ldim_nptcls(micname, ldim, ifoo)
        call micrograph%new(ldim, smpd)
        call micrograph%read(micname)
        ! find reference dimensions
        call find_ldim_nptcls(refsname, ldim_refs, nrefs)
        orig_box              = ldim_refs(1)
        ! modify according to PICKER_SHRINK & PICKER_SHRINK_REFINE
        ldim_refs_refine      = ldim_refs
        ldim_refs(1)          = round2even(real(ldim_refs(1))/PICKER_SHRINK)
        ldim_refs(2)          = round2even(real(ldim_refs(2))/PICKER_SHRINK)
        ldim_refs(3)          = 1
        ldim_refs_refine(1)   = round2even(real(ldim_refs_refine(1))/PICKER_SHRINK_REFINE)
        ldim_refs_refine(2)   = round2even(real(ldim_refs_refine(2))/PICKER_SHRINK_REFINE)
        ldim_refs_refine(3)   = 1
        ldim_shrink(1)        = round2even(real(ldim(1))/PICKER_SHRINK)
        ldim_shrink(2)        = round2even(real(ldim(2))/PICKER_SHRINK)
        ldim_shrink(3)        = 1
        ldim_shrink_refine(1) = round2even(real(ldim(1))/PICKER_SHRINK_REFINE)
        ldim_shrink_refine(2) = round2even(real(ldim(2))/PICKER_SHRINK_REFINE)
        ldim_shrink_refine(3) = 1
        nx                    = ldim_shrink(1)-ldim_refs(1)
        ny                    = ldim_shrink(2)-ldim_refs(2)
        nx_refine             = ldim_shrink_refine(1)-ldim_refs_refine(1)
        ny_refine             = ldim_shrink_refine(2)-ldim_refs_refine(2)
        smpd_shrunken         = PICKER_SHRINK*smpd
        smpd_shrunken_refine  = PICKER_SHRINK_REFINE*smpd
        msk                   = real(ldim_refs(1)/2-5)
        msk_refine            = real(ldim_refs_refine(1)/2-5)
        distthr               = BOXFRAC*real(ldim_refs(1))
        if( present(distthr_in) ) distthr = distthr_in/PICKER_SHRINK
        ! read and shrink references
        allocate( refs(nrefs), refs_refine(nrefs), sxx(nrefs), sxx_refine(nrefs), stat=alloc_stat )
        call alloc_errchk( "In: simple_picker :: init_picker, 1", alloc_stat)
        do iref=1,nrefs
            call refs(iref)%new(ldim_refs, smpd_shrunken)
            call refs_refine(iref)%new(ldim_refs_refine, smpd_shrunken_refine)
            call refimg%new([orig_box,orig_box,1], smpd)
            call refimg%read(refsname, iref)
            call refimg%fwd_ft
            call refimg%clip(refs(iref))
            call refimg%clip(refs_refine(iref))
            call refs(iref)%bwd_ft
            call refs_refine(iref)%bwd_ft
            call refs(iref)%mask(msk, 'hard')
            call refs_refine(iref)%mask(msk_refine, 'hard')
            call refs(iref)%prenorm4real_corr(sxx(iref))
            call refs_refine(iref)%prenorm4real_corr(sxx_refine(iref))
        end do
        ! pre-process micrograph
        call micrograph%fwd_ft
        call micrograph%bp(0., lp)
        call mic_shrunken%new(ldim_shrink, smpd_shrunken)
        call mic_shrunken_refine%new(ldim_shrink_refine, smpd_shrunken_refine)
        call micrograph%clip(mic_shrunken)
        call micrograph%clip(mic_shrunken_refine)
        call mic_shrunken%bwd_ft
        call mic_shrunken_refine%bwd_ft
        if( WRITESHRUNKEN ) call mic_shrunken%write('shrunken.mrc')
    end subroutine init_picker

    subroutine exec_picker( boxname_out, nptcls_out )
        character(len=STDLEN), intent(out) :: boxname_out
        integer,               intent(out) :: nptcls_out
        call extract_peaks
        call distance_filter
        call refine_positions
        call gather_stats
        call one_cluster_clustering
        ! if( rm_outliers ) call remove_outliers
        nptcls_out = count(selected_peak_positions)
        ! bring back coordinates to original sampling
        peak_positions_refined = nint(PICKER_SHRINK_REFINE)*peak_positions_refined
        call write_boxfile
        boxname_out = boxname
    end subroutine exec_picker

    subroutine extract_peaks
        real    :: means(2), corrs(nrefs)
        integer :: xind, yind, alloc_stat, funit, iref, i, loc(1), ind
        integer, allocatable :: labels(:), target_positions(:,:)
        real,    allocatable :: target_corrs(:), spec(:)
        logical :: outside
        write(*,'(a)') '>>> EXTRACTING PEAKS'
        ntargets = 0
        do xind=0,nx,PICKER_OFFSET
            do yind=0,ny,PICKER_OFFSET
                ntargets = ntargets + 1
            end do
        end do
        allocate( target_corrs(ntargets), target_positions(ntargets,2),&
                  corrmat(0:nx,0:ny), refmat(0:nx,0:ny), stat=alloc_stat )
        call alloc_errchk( 'In: simple_picker :: gen_corr_peaks, 1', alloc_stat )
        target_corrs     = 0.
        target_positions = 0
        corrmat          = -1.
        refmat           = 0
        ntargets         = 0
        corrmax          = -1.
        corrmin          = 1.
        call ptcl_target%new(ldim_refs, smpd_shrunken)
        do xind=0,nx,PICKER_OFFSET
            do yind=0,ny,PICKER_OFFSET
                ntargets = ntargets + 1
                target_positions(ntargets,:) = [xind,yind]
                call mic_shrunken%window_slim([xind,yind], ldim_refs(1), ptcl_target, outside)
                !$omp parallel do schedule(static) default(shared) private(iref) proc_bind(close)
                do iref=1,nrefs
                    corrs(iref) = refs(iref)%real_corr_prenorm(ptcl_target, sxx(iref))
                end do
                !$omp end parallel do
                loc = maxloc(corrs)
                target_corrs(ntargets) = corrs(loc(1))
                corrmat(xind,yind)     = target_corrs(ntargets)
                refmat(xind,yind)      = loc(1)
                if( target_corrs(ntargets) > corrmax ) corrmax = target_corrs(ntargets)
                if( target_corrs(ntargets) < corrmin ) corrmin = target_corrs(ntargets)
            end do
        end do
        call sortmeans(target_corrs, MAXKMIT, means, labels)
        npeaks = count(labels == 2)
        allocate( peak_positions(npeaks,2), stat=alloc_stat)
        call alloc_errchk( 'In: simple_picker :: gen_corr_peaks, 2', alloc_stat )
        peak_positions = 0
        npeaks = 0
        do i=1,ntargets
            if( labels(i) == 2 )then
                npeaks = npeaks + 1
                peak_positions(npeaks,:) = target_positions(i,:)
            endif
        end do
    end subroutine extract_peaks

    subroutine distance_filter
        integer :: ipeak, jpeak, ipos(2), jpos(2), alloc_stat, loc(1)
        real    :: dist
        logical, allocatable :: mask(:)
        real,    allocatable :: corrs(:)
        write(*,'(a)') '>>> DISTANCE FILTERING'
        allocate( mask(npeaks), corrs(npeaks), selected_peak_positions(npeaks), stat=alloc_stat)
        call alloc_errchk( 'In: simple_picker :: distance_filter', alloc_stat )
        selected_peak_positions = .true.
        do ipeak=1,npeaks
            ipos = peak_positions(ipeak,:)
            mask = .false.
            !$omp parallel do schedule(static) default(shared) private(jpeak,jpos,dist) proc_bind(close)
            do jpeak=1,npeaks
                jpos = peak_positions(jpeak,:)
                dist = euclid(real(ipos),real(jpos))
                if( dist < distthr ) mask(jpeak) = .true.
                corrs(jpeak) = corrmat(jpos(1),jpos(2))
            end do
            !$omp end parallel do
            ! find best match in the neigh
            loc = maxloc(corrs, mask=mask)
            ! eliminate all but the best
            mask(loc(1)) = .false.
            where( mask )
                selected_peak_positions = .false.
            end where
        end do
        npeaks_sel = count(selected_peak_positions)
        write(*,'(a,1x,I5)') 'peak positions left after distance filtering: ', npeaks_sel
    end subroutine distance_filter

    subroutine refine_positions
        integer :: ipeak, xrange(2), yrange(2), xind, yind, ref, cnt
        real    :: corr, target_corr
        logical :: outside
        write(*,'(a)') '>>> REFINING POSITIONS & GATHERING FIRST STATS'
        allocate( peak_stats(npeaks_sel,NSTAT) )
        ! bring back coordinates to refinement sampling
        allocate( peak_positions_refined(npeaks,2), source=nint(PICKER_SHRINK/PICKER_SHRINK_REFINE)*peak_positions)
        call ptcl_target%new(ldim_refs_refine, smpd_shrunken_refine)
        cnt = 0
        do ipeak=1,npeaks
            if( selected_peak_positions(ipeak) )then
                cnt = cnt + 1
                ! best match in crude first scan
                ref = refmat(peak_positions(ipeak,1),peak_positions(ipeak,2))
                ! refinement range
                call srch_range(peak_positions_refined(ipeak,:))
                ! extract image, correlate, find peak
                corr = -1
                do xind=xrange(1),xrange(2)
                    do yind=yrange(1),yrange(2)
                        call mic_shrunken_refine%window_slim([xind,yind], ldim_refs_refine(1), ptcl_target, outside)
                        target_corr = refs_refine(ref)%real_corr_prenorm(ptcl_target, sxx_refine(ref))
                        if( target_corr > corr )then
                            peak_positions_refined(ipeak,:) = [xind,yind]
                            corr = target_corr
                        endif
                    end do
                end do
                peak_stats(cnt,CC2REF) = corr
            endif
        end do

        contains

            subroutine srch_range( pos )
                integer, intent(in) :: pos(2)
                xrange(1) = max(0,         pos(1) - PICKER_OFFSET*nint(PICKER_SHRINK/PICKER_SHRINK_REFINE))
                xrange(2) = min(nx_refine, pos(1) + PICKER_OFFSET*nint(PICKER_SHRINK/PICKER_SHRINK_REFINE))
                yrange(1) = max(0,         pos(2) - PICKER_OFFSET*nint(PICKER_SHRINK/PICKER_SHRINK_REFINE))
                yrange(2) = min(ny_refine, pos(2) + PICKER_OFFSET*nint(PICKER_SHRINK/PICKER_SHRINK_REFINE))
            end subroutine srch_range

    end subroutine refine_positions

    subroutine gather_stats
        use simple_stat, only: normalize_minmax
        integer           :: ipeak, cnt, istat
        logical           :: outside
        real              :: ave, maxv, minv
        real, allocatable :: spec(:)
        write(*,'(a)') '>>> GATHERING REMAINING STATS'
        call ptcl_target%new(ldim_refs_refine, smpd_shrunken_refine)
        cnt = 0
        do ipeak=1,npeaks
            if( selected_peak_positions(ipeak) )then
                cnt = cnt + 1
                call mic_shrunken_refine%window_slim(peak_positions_refined(ipeak,:),&
                    &ldim_refs_refine(1), ptcl_target, outside)
                spec = ptcl_target%spectrum('power')
                peak_stats(cnt,SSCORE) = sum(spec)/real(size(spec))
                call ptcl_target%stats('background', ave, peak_stats(cnt,SDEV), maxv, minv)
                peak_stats(cnt,DYNRANGE) = maxv - minv
            endif
        end do
        ! min/max normalise to get all vars on equal footing
        do istat=1,NSTAT
            call normalize_minmax(peak_stats(:,istat))
        end do
    end subroutine gather_stats

    subroutine one_cluster_clustering
        use simple_stat, only: dev_from_dmat
        real, allocatable :: dmat(:,:)
        integer           :: i_median, i, j, nnincl, cnt, ipeak
        real              :: dmed, ddev
        allocate(dmat(npeaks_sel,npeaks_sel), source=0.)
        do i=1,npeaks_sel - 1
            do j=i + 1,npeaks_sel
                dmat(i,j) = euclid(peak_stats(i,:), peak_stats(j,:))
                dmat(j,i) = dmat(i,j)
            end do
        end do
        call dev_from_dmat( dmat, i_median, ddev )
        cnt = 0
        do ipeak=1,npeaks
            if( selected_peak_positions(ipeak) )then
                cnt = cnt + 1
                if( dmat(i_median,cnt) <=  ndev * ddev )then
                    ! we are keeping this one
                else
                    ! we are removing this one
                    selected_peak_positions(ipeak) = .false.
                endif
            endif
        end do
        npeaks_sel = count(selected_peak_positions)
        write(*,'(a,1x,I5)') 'peak positions left after one cluster clustering: ', npeaks_sel
    end subroutine one_cluster_clustering

    subroutine write_boxfile
        use simple_fileio, only: fopen, fclose, fileio_errmsg
        integer :: funit, ipeak,iostat
        if(.not.fopen(funit, status='REPLACE', action='WRITE', file=boxname,iostat=iostat))&
             call fileio_errmsg('picker; write_boxfile ', iostat)
        do ipeak=1,npeaks
            if( selected_peak_positions(ipeak) )then             
                write(funit,'(I7,I7,I7,I7,I7)') peak_positions_refined(ipeak,1),&
                peak_positions_refined(ipeak,2), orig_box, orig_box, -3
            endif
        end do
        if(.not.fclose(funit,iostat=iostat))&
             call fileio_errmsg('picker; write_boxfile end', iostat)
    end subroutine write_boxfile

    subroutine kill_picker
        integer :: iref, alloc_stat
        if( allocated(micname) )then
            call micrograph%kill
            call mic_shrunken%kill
            call mic_shrunken_refine%kill
            call ptcl_target%kill
            deallocate(selected_peak_positions,sxx,sxx_refine,corrmat,peak_positions,&
                &peak_positions_refined,refmat,micname,refsname,peak_stats)
            do iref=1,nrefs
                call refs(iref)%kill
                call refs_refine(iref)%kill
            end do
            deallocate(refs, refs_refine)
        endif
    end subroutine kill_picker

end module simple_picker
