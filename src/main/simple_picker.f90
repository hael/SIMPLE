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
type(image)                   :: micrograph, mic_shrunken, mic_shrunken_copy, mic_shrunken_refine, ptcl_target, mic_saved
type(image),      allocatable :: refs(:), refs_refine(:)
logical,          allocatable :: selected_peak_positions(:)
real,             allocatable :: sxx(:), sxx_refine(:), corrmat(:,:), peak_stats(:,:)
integer,          allocatable :: peak_positions(:,:), peak_positions_refined(:,:), refmat(:,:)
character(len=:), allocatable :: micname, refsname
character(len=LONGSTRLEN)     :: boxname
integer                       :: ldim(3), ldim_refs(3), ldim_refs_refine(3), ldim_shrink(3)
integer                       :: ldim_shrink_refine(3), ntargets, nx, ny, nx_refine, ny_refine
integer                       :: nrefs, nmax, nmax_sel, orig_box, cnt_glob=0
real                          :: smpd, smpd_shrunken, smpd_shrunken_refine, corrmax, corrmin
real                          :: msk, msk_refine, lp, distthr, ndev

contains

    subroutine init_picker( micfname, refsfname, smpd_in, lp_in, distthr_in, ndev_in, dir_out )
        character(len=*),           intent(in) :: micfname, refsfname
        real,                       intent(in) :: smpd_in
        real,             optional, intent(in) :: lp_in, distthr_in, ndev_in
        character(len=*), optional, intent(in) :: dir_out
        type(image)       :: refimg
        integer           :: ifoo, iref
        real              :: sigma, hp
        allocate(micname,  source=trim(micfname))
        allocate(refsname, source=trim(refsfname))
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
        msk                   = max(msk, real(ldim_refs(1)/2-2)) ! for tiny particles
        msk_refine            = real(ldim_refs_refine(1)/2-5)
        msk_refine            = max(PICKER_SHRINK/PICKER_SHRINK_REFINE*real(ldim_refs(1)/2-2), msk_refine)
        distthr               = BOXFRAC*real(ldim_refs(1))
        if( present(distthr_in) ) distthr = distthr_in/smpd_shrunken
        ! read and shrink references
        allocate( refs(nrefs), refs_refine(nrefs), sxx(nrefs), sxx_refine(nrefs) )
        do iref=1,nrefs
            call refs(iref)%new(ldim_refs, smpd_shrunken)
            call refs_refine(iref)%new(ldim_refs_refine, smpd_shrunken_refine)
            call refimg%new([orig_box,orig_box,1], smpd)
            call refimg%read(refsname, iref)
            call refimg%fft()
            call refimg%clip(refs(iref))
            call refimg%clip(refs_refine(iref))
            call refs(iref)%ifft()
            call refs_refine(iref)%ifft()
            call refs(iref)%mask(msk, 'hard')
            call refs_refine(iref)%mask(msk_refine, 'hard')
            call refs(iref)%prenorm4real_corr(sxx(iref))
            call refs_refine(iref)%prenorm4real_corr(sxx_refine(iref))
        end do
        call refimg%kill
        ! pre-process micrograph
        call micrograph%fft()
        call mic_shrunken%new(ldim_shrink, smpd_shrunken)
        call mic_shrunken_refine%new(ldim_shrink_refine, smpd_shrunken_refine)
        call mic_shrunken%set_ft(.true.)
        call mic_shrunken_refine%set_ft(.true.)
        call micrograph%clip(mic_shrunken)
        call micrograph%clip(mic_shrunken_refine)
        hp = real(ldim_shrink(1) / 2) * smpd_shrunken
        call mic_shrunken%bp(hp, lp)
        hp = real(ldim_shrink_refine(1) / 2) * smpd_shrunken_refine
        mic_saved = mic_shrunken
        call mic_shrunken_refine%bp(hp, lp)
        call mic_shrunken%ifft()
        call mic_saved%ifft()
        call mic_shrunken_refine%ifft()
    end subroutine init_picker

    subroutine exec_picker( boxname_out, nptcls_out )
        character(len=LONGSTRLEN), intent(out) :: boxname_out
        integer,                   intent(out) :: nptcls_out
        real, pointer :: prmat(:,:,:)
        integer       :: i
        real          :: maxv
        call extract_peaks
        call distance_filter
        call refine_positions
        call gather_stats
        call one_cluster_clustering
        nptcls_out = count(selected_peak_positions)
        call mic_saved%kill
        ! bring back coordinates to original sampling
        peak_positions_refined = nint(PICKER_SHRINK_REFINE)*peak_positions_refined
        call write_boxfile
        ! returns absolute path
        call make_relativepath(CWD_GLOB, boxname, boxname_out)
    end subroutine exec_picker

    subroutine extract_peaks
        real    :: means(2), corrs(nrefs)
        integer :: xind, yind, iref, i, loc(1)
        integer, allocatable :: labels(:), target_positions(:,:)
        real,    allocatable :: target_corrs(:)
        logical :: outside
        write(logfhandle,'(a)') '>>> EXTRACTING PEAKS'
        ntargets = 0
        do xind=0,nx,PICKER_OFFSET
            do yind=0,ny,PICKER_OFFSET
                ntargets = ntargets + 1
            end do
        end do
        allocate( target_corrs(ntargets), target_positions(ntargets,2),&
                  corrmat(0:nx,0:ny), refmat(0:nx,0:ny))
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
        integer :: ipeak, jpeak, ipos(2), jpos(2), loc(1)
        real    :: dist
        logical, allocatable :: mask(:)
        real,    allocatable :: corrs(:)
        write(logfhandle,'(a)') '>>> DISTANCE FILTERING'
        allocate( mask(nmax), corrs(nmax), selected_peak_positions(nmax) )
        selected_peak_positions = .true.
        do ipeak=1,nmax
            ipos = peak_positions(ipeak,:)
            mask = .false.
            !$omp parallel do schedule(static) default(shared) private(jpeak,jpos,dist) proc_bind(close)
            do jpeak=1,nmax
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
        nmax_sel = count(selected_peak_positions)
        write(logfhandle,'(a,1x,I5)') 'peak positions left after distance filtering: ', nmax_sel
    end subroutine distance_filter

    subroutine refine_positions
        integer :: ipeak, xrange(2), yrange(2), xind, yind, ref, cnt
        real    :: corr, target_corr
        logical :: outside
        write(logfhandle,'(a)') '>>> REFINING POSITIONS & GATHERING FIRST STATS'
        allocate( peak_stats(nmax_sel,NSTAT) )
        ! bring back coordinates to refinement sampling
        allocate( peak_positions_refined(2,nmax), source=nint(PICKER_SHRINK/PICKER_SHRINK_REFINE)*transpose(peak_positions))
        call ptcl_target%new(ldim_refs_refine, smpd_shrunken_refine)
        cnt = 0
        do ipeak=1,nmax
            if( selected_peak_positions(ipeak) )then
                cnt = cnt + 1
                ! best match in crude first scan
                ref = refmat(peak_positions(ipeak,1),peak_positions(ipeak,2))
                ! refinement range
                call srch_range(peak_positions_refined(:,ipeak))
                ! extract image, correlate, find peak
                corr = -1
                do xind=xrange(1),xrange(2)
                    do yind=yrange(1),yrange(2)
                        call mic_shrunken_refine%window_slim([xind,yind], ldim_refs_refine(1), ptcl_target, outside)
                        target_corr = refs_refine(ref)%real_corr_prenorm(ptcl_target, sxx_refine(ref))
                        if( target_corr > corr )then
                            peak_positions_refined(:,ipeak) = [xind,yind]
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
        integer           :: ipeak, cnt, istat
        logical           :: outside
        real              :: ave, maxv, minv
        real, allocatable :: spec(:)
        write(logfhandle,'(a)') '>>> GATHERING REMAINING STATS'
        call ptcl_target%new(ldim_refs_refine, smpd_shrunken_refine)
        cnt = 0
        do ipeak=1,nmax
            if( selected_peak_positions(ipeak) )then
                cnt = cnt + 1
                call mic_shrunken_refine%window_slim(peak_positions_refined(:,ipeak),&
                    &ldim_refs_refine(1), ptcl_target, outside)
                call ptcl_target%spectrum('power', spec)
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
        call dev_from_dmat( dmat, i_median, ddev )
        cnt = 0
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
        write(logfhandle,'(a,1x,I5)') 'peak positions left after one cluster clustering: ', nmax_sel
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
        integer :: iref
        if( allocated(micname) )then
            call micrograph%kill
            call mic_shrunken%kill
            call mic_shrunken_refine%kill
            call ptcl_target%kill
            deallocate(selected_peak_positions,sxx,sxx_refine,corrmat)
            deallocate(peak_positions,peak_positions_refined,refmat,micname,refsname,peak_stats)
            do iref=1,nrefs
                call refs(iref)%kill
                call refs_refine(iref)%kill
            end do
            deallocate(refs, refs_refine)
        endif
    end subroutine kill_picker

end module simple_picker
