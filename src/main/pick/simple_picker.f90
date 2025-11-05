! particle picker
module simple_picker
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: init_picker, exec_picker, kill_picker
private

! PEAK STATS INDICES
integer,          parameter   :: CC2REF   = 1
integer,          parameter   :: SDEV     = 2
integer,          parameter   :: DYNRANGE = 3
integer,          parameter   :: SSCORE   = 4
! OTHER PARAMS
integer,          parameter   :: NSTAT    = 4
integer,          parameter   :: MAXKMIT  = 20
real,             parameter   :: BOXFRAC  = 0.5
! VARS
type(image)                   :: micrograph, mic_shrunken, mic_shrunken_copy, mic_shrunken_refine
type(image),      allocatable :: refs(:), refs_refine(:), ptcls(:)
logical,          allocatable :: selected_peak_positions(:), l_refs_err(:), l_refs_refine_err(:)
real,             allocatable :: corrmat(:,:), peak_stats(:,:)
integer,          allocatable :: peak_positions(:,:), peak_positions_refined(:,:), refmat(:,:)
character(len=:), allocatable :: micname
character(len=LONGSTRLEN)     :: boxname
integer                       :: ldim(3), ldim_refs(3), ldim_refs_refine(3), ldim_shrink(3)
integer                       :: ldim_shrink_refine(3), ntargets, nx, ny, nx_refine, ny_refine
integer                       :: nrefs, nmax, nmax_sel, orig_box, cnt_glob=0
real                          :: smpd, smpd_shrunken, smpd_shrunken_refine
real                          :: msk, msk_refine, lp, distthr, ndev

contains

    subroutine init_picker( micfname, pickrefs, smpd_in, lp_in, distthr_in, ndev_in, dir_out )
        character(len=*),           intent(in) :: micfname
        type(image),   allocatable, intent(in) :: pickrefs(:)
        real,                       intent(in) :: smpd_in
        real,             optional, intent(in) :: lp_in, distthr_in, ndev_in
        character(len=*), optional, intent(in) :: dir_out
        type(image)       :: refimg
        integer           :: orig_ldim_refs(3), ifoo, iref
        real              :: hp
        allocate(micname,  source=trim(micfname))
        boxname = basename( fname_new_ext(micname,'box') )
        if( present(dir_out) )boxname = trim(dir_out)//trim(boxname)
        smpd = smpd_in
        lp   = 20.0
        if( present(lp_in) ) lp = lp_in
        ndev = 2.0
        if( present(ndev_in)) ndev = ndev_in
        ! read micrograph
        call find_ldim_nptcls(micname, ldim, ifoo)
        call micrograph%new(ldim, smpd)
        call micrograph%read(micname)
        ! find reference dimensions
        orig_ldim_refs = pickrefs(1)%get_ldim()
        nrefs          = size(pickrefs)
        orig_box       = orig_ldim_refs(1)
        ! modify according to PICKER_SHRINK & PICKER_SHRINK_REFINE
        ldim_refs(1)          = round2even(real(orig_ldim_refs(1))/PICKER_SHRINK)
        ldim_refs(2)          = round2even(real(orig_ldim_refs(2))/PICKER_SHRINK)
        ldim_refs(3)          = 1
        ldim_refs_refine(1)   = round2even(real(orig_ldim_refs(1))/PICKER_SHRINK_REFINE)
        ldim_refs_refine(2)   = round2even(real(orig_ldim_refs(2))/PICKER_SHRINK_REFINE)
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
        msk                   = max(msk, real(ldim_refs(1)/2-2)) ! for tiny particles
        msk_refine            = real(ldim_refs_refine(1)/2-5)
        msk_refine            = max(PICKER_SHRINK/PICKER_SHRINK_REFINE*real(ldim_refs(1)/2-2), msk_refine)
        distthr               = BOXFRAC*real(ldim_refs(1))
        if( present(distthr_in) ) distthr = distthr_in/smpd_shrunken
        ! read and shrink references
        allocate(refs(nrefs),refs_refine(nrefs),l_refs_err(nrefs),l_refs_refine_err(nrefs))
        call refimg%new([orig_box,orig_box,1], smpd)
        do iref=1,nrefs
            call refs(iref)%new(ldim_refs, smpd_shrunken)
            call refs_refine(iref)%new(ldim_refs_refine, smpd_shrunken_refine)
            call refimg%copy_fast(pickrefs(iref))
            call refimg%fft()
            call refimg%clip(refs(iref))
            call refimg%clip(refs_refine(iref))
        end do
        call refimg%kill
        ! could be multi-threaded?
        do iref=1,nrefs
            call refs(iref)%ifft()
            call refs_refine(iref)%ifft()
            call refs(iref)%mask(msk, 'hard')
            call refs_refine(iref)%mask(msk_refine, 'hard')
            call refs(iref)%prenorm4real_corr(l_refs_err(iref))
            call refs_refine(iref)%prenorm4real_corr(l_refs_refine_err(iref))
        end do
        allocate(ptcls(nthr_glob))
        ! pre-process micrograph
        call micrograph%fft()
        call mic_shrunken%new(ldim_shrink, smpd_shrunken)
        call mic_shrunken_refine%new(ldim_shrink_refine, smpd_shrunken_refine)
        call micrograph%clip(mic_shrunken)
        call micrograph%clip(mic_shrunken_refine)
        hp = real(ldim_shrink(1) / 2) * smpd_shrunken
        call mic_shrunken%bp(hp, lp)
        hp = real(ldim_shrink_refine(1) / 2) * smpd_shrunken_refine
        call mic_shrunken_refine%bp(hp, lp)
        call mic_shrunken%ifft()
        call mic_shrunken_refine%ifft()
    end subroutine init_picker

    subroutine exec_picker( boxname_out, nptcls_out )
        character(len=LONGSTRLEN), intent(out) :: boxname_out
        integer,                   intent(out) :: nptcls_out
        call extract_peaks
        call distance_filter
        call refine_positions
        call gather_stats
        call one_cluster_clustering
        nptcls_out = count(selected_peak_positions)
        ! bring back coordinates to original sampling
        peak_positions_refined = nint(PICKER_SHRINK_REFINE)*peak_positions_refined
        call write_boxfile
        ! returns absolute path
        boxname_out = simple_abspath(boxname)
    end subroutine exec_picker

    subroutine extract_peaks
        integer, allocatable :: labels(:), target_positions(:,:)
        real,    allocatable :: target_corrs(:)
        real    :: means(2), corrmax, corr
        integer :: xind, yind, iref, i, loc, ithr, ntx, nty
        logical :: outside, l_err
        write(logfhandle,'(a)') '>>> EXTRACTING PEAKS'
        ntx = floor(real(nx+PICKER_OFFSET)/real(PICKER_OFFSET))
        nty = floor(real(ny+PICKER_OFFSET)/real(PICKER_OFFSET))
        ntargets = ntx*nty
        allocate( target_corrs(ntargets), target_positions(ntargets,2),&
                    corrmat(0:nx,0:ny), refmat(0:nx,0:ny))
        target_corrs     = 0.
        target_positions = 0
        corrmat          = -1.
        refmat           = 0
        !$omp parallel do schedule(static) default(shared) private(iref,ithr,ntargets,xind,yind,corr,corrmax,loc,outside,l_err)&
        !$omp proc_bind(close) collapse(2)
        do xind = 0,nx,PICKER_OFFSET
            do yind = 0,ny,PICKER_OFFSET
                ithr     = omp_get_thread_num() + 1
                ntargets = xind/PICKER_OFFSET * nty + (yind+PICKER_OFFSET)/PICKER_OFFSET
                target_positions(ntargets,:) = [xind,yind]
                if( .not. ptcls(ithr)%exists() ) call ptcls(ithr)%new(ldim_refs, smpd_shrunken, wthreads=.false.)
                call mic_shrunken%window_slim([xind,yind], ldim_refs(1), ptcls(ithr), outside)
                call ptcls(ithr)%prenorm4real_corr(l_err)
                corrmax = -1.0
                loc     = 0
                if( .not.l_err )then
                    do iref=1,nrefs
                        if( l_refs_err(iref) ) cycle
                        corr = refs(iref)%real_corr_prenorm(ptcls(ithr))
                        if( corr > corrmax )then
                            corrmax = corr
                            loc     = iref
                        endif
                    end do
                endif
                target_corrs(ntargets) = corrmax
                corrmat(xind,yind)     = corrmax
                refmat(xind,yind)      = loc
            enddo
        enddo
        !$omp end parallel do
        call sortmeans(target_corrs, MAXKMIT, means, labels)
        nmax = count(labels == 2)
        allocate( peak_positions(nmax,2) )
        peak_positions = 0
        nmax = 0
        do i=1,ntargets
            if( labels(i) == 2 )then
                nmax = nmax + 1
                peak_positions(nmax,:) = target_positions(i,:)
            endif
        end do
    end subroutine extract_peaks

    subroutine distance_filter
        logical, allocatable :: mask(:)
        real,    allocatable :: corrs(:)
        integer :: ipos(2), jpos(2), loc(1), ipeak, jpeak
        real    :: distsq, distthrsq
        write(logfhandle,'(a)') '>>> DISTANCE FILTERING'
        allocate( mask(nmax), corrs(nmax), selected_peak_positions(nmax) )
        selected_peak_positions = .true.
        distthrsq = distthr**2
        do ipeak=1,nmax
            ipos = peak_positions(ipeak,:)
            !$omp parallel do schedule(static) default(shared) private(jpeak,jpos,distsq) proc_bind(close)
            do jpeak=1,nmax
                jpos   = peak_positions(jpeak,:)
                distsq = real(sum((jpos-ipos)**2))
                if( distsq < distthrsq )then
                    mask(jpeak)  = .true.
                    corrs(jpeak) = corrmat(jpos(1),jpos(2))
                else
                    mask(jpeak) = .false.
                endif
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
        nmax_sel = count(selected_peak_positions)
        write(logfhandle,'(a,1x,I5)') '>>> PEAK POSITIONS AFTER DISTANCE FILTERING: ', nmax_sel
    end subroutine distance_filter

    subroutine refine_positions
        integer :: xrange(2), yrange(2), xind, yind, ref, cnt, ithr, ipeak
        real    :: corr, target_corr
        logical :: outside, l_err
        write(logfhandle,'(a)') '>>> REFINING POSITIONS & GATHERING STATS'
        allocate( peak_stats(nmax_sel,NSTAT) )
        ! bring back coordinates to refinement sampling
        allocate( peak_positions_refined(2,nmax), source=nint(PICKER_SHRINK/PICKER_SHRINK_REFINE)*transpose(peak_positions))
        do ithr = 1,nthr_glob
            call ptcls(ithr)%new(ldim_refs_refine, smpd_shrunken_refine, wthreads=.false.)
        enddo
        !$omp parallel do private(ipeak,cnt,ref,xind,yind,ithr,target_corr,corr,outside,l_err,xrange,yrange)&
        !$omp schedule(static) default(shared) proc_bind(close)
        do ipeak=1,nmax
            if( selected_peak_positions(ipeak) )then
                ithr = omp_get_thread_num() + 1
                cnt  = merge(1, count(selected_peak_positions(:ipeak-1))+1, ipeak==1)
                ! best match in crude first scan
                ref = refmat(peak_positions(ipeak,1),peak_positions(ipeak,2))
                corr = -1.
                if( .not.l_refs_refine_err(ref) )then
                    ! refinement range
                    call srch_range(peak_positions_refined(:,ipeak), xrange, yrange)
                    ! extract image, correlate, find peak
                    do xind=xrange(1),xrange(2)
                        do yind=yrange(1),yrange(2)
                            call mic_shrunken_refine%window_slim([xind,yind], ldim_refs_refine(1), ptcls(ithr), outside)
                            call ptcls(ithr)%prenorm4real_corr(l_err)
                            if( l_err ) cycle
                            target_corr = refs_refine(ref)%real_corr_prenorm(ptcls(ithr))
                            if( target_corr > corr )then
                                peak_positions_refined(:,ipeak) = [xind,yind]
                                corr = target_corr
                            endif
                        end do
                    end do
                endif
                peak_stats(cnt,CC2REF) = corr
            endif
        end do
        !$omp end parallel do
        contains

            subroutine srch_range( pos, xr, yr )
                integer, intent(in)  :: pos(2)
                integer, intent(out) :: xr(2), yr(2)
                xr(1) = max(0,         pos(1) - PICKER_OFFSET*nint(PICKER_SHRINK/PICKER_SHRINK_REFINE))
                xr(2) = min(nx_refine, pos(1) + PICKER_OFFSET*nint(PICKER_SHRINK/PICKER_SHRINK_REFINE))
                yr(1) = max(0,         pos(2) - PICKER_OFFSET*nint(PICKER_SHRINK/PICKER_SHRINK_REFINE))
                yr(2) = min(ny_refine, pos(2) + PICKER_OFFSET*nint(PICKER_SHRINK/PICKER_SHRINK_REFINE))
            end subroutine srch_range

    end subroutine refine_positions

    subroutine gather_stats
        real, allocatable :: spec(:)
        real              :: ave, maxv, minv
        integer           :: ipeak, cnt, istat, ithr
        logical           :: outside
        allocate(spec(fdim(ldim_refs_refine(1)) - 1))
        !$omp parallel do private(ithr,ipeak,cnt,spec,outside,ave,minv,maxv)&
        !$omp schedule(static) default(shared) proc_bind(close)
        do ipeak=1,nmax
            if( selected_peak_positions(ipeak) )then
                ithr = omp_get_thread_num() + 1
                cnt  = merge(1, count(selected_peak_positions(:ipeak-1))+1, ipeak==1)
                call ptcls(ithr)%set_ft(.false.)
                call mic_shrunken_refine%window_slim(peak_positions_refined(:,ipeak), ldim_refs_refine(1), ptcls(ithr), outside)
                call ptcls(ithr)%stats('background', ave, peak_stats(cnt,SDEV), maxv, minv)
                peak_stats(cnt,DYNRANGE) = maxv - minv
                call ptcls(ithr)%fft
                call ptcls(ithr)%power_spectrum(spec)
                peak_stats(cnt,SSCORE) = sum(spec)/real(size(spec))
            endif
        end do
        !$omp end parallel do
        ! min/max normalise to get all vars on equal footing
        do istat=1,NSTAT
            call normalize_minmax(peak_stats(:,istat))
        end do
    end subroutine gather_stats

    subroutine one_cluster_clustering
        real, allocatable :: dmat(:,:)
        integer           :: i_median, i, j, cnt, ipeak
        real              :: ddev
        allocate(dmat(nmax_sel,nmax_sel), source=0.)
        do i=1,nmax_sel - 1
            do j=i + 1,nmax_sel
                dmat(i,j) = euclid(peak_stats(i,:), peak_stats(j,:))
                dmat(j,i) = dmat(i,j)
            end do
        end do
        call medoid_from_dmat(dmat, i_median)
        ddev = median(dmat(i_median,:))
        cnt  = 0
        do ipeak=1,nmax
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
        nmax_sel = count(selected_peak_positions)
        write(logfhandle,'(a,1x,I5)') '>>> PEAK POSITIONS AFTER ONE-CLUSTER CLUSTERING: ', nmax_sel
    end subroutine one_cluster_clustering

    subroutine write_boxfile
        integer :: funit, ipeak,iostat
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(boxname)),iostat=iostat)
        call fileiochk('picker; write_boxfile ', iostat)
        do ipeak=1,nmax
            if( selected_peak_positions(ipeak) )then
                write(funit,'(I7,I7,I7,I7,I7)') peak_positions_refined(1,ipeak),&
                peak_positions_refined(2,ipeak), orig_box, orig_box, -3
            endif
        end do
        call fclose(funit)
    end subroutine write_boxfile

    subroutine kill_picker
        integer :: iref, ithr
        if( allocated(micname) )then
            call micrograph%kill
            call mic_shrunken%kill
            call mic_shrunken_refine%kill
            deallocate(selected_peak_positions,corrmat,l_refs_err,l_refs_refine_err)
            deallocate(peak_positions,peak_positions_refined,refmat,micname,peak_stats)
            do iref=1,nrefs
                call refs(iref)%kill
                call refs_refine(iref)%kill
            end do
            deallocate(refs, refs_refine)
            do ithr=1,nthr_glob
                call ptcls(ithr)%kill
            end do
            deallocate(ptcls)
        endif
    end subroutine kill_picker

end module simple_picker
