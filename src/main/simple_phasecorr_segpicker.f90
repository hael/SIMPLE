! particle picker
module simple_phasecorr_segpicker
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,        only: image
implicit none

public :: init_phasecorr_segpicker, exec_phasecorr_segpicker, kill_phasecorr_segpicker
private

! PEAK STATS INDICES
integer,          parameter   :: SDEV     = 2
integer,          parameter   :: DYNRANGE = 3
integer,          parameter   :: SSCORE   = 4
real,             parameter   :: PICKER_SHRINK_HERE   = 4.
! OTHER PARAMS
integer,          parameter   :: NSTAT   = 4
integer,          parameter   :: MAXKMIT = 20
real,             parameter   :: BOXFRAC = 0.5
logical,          parameter   :: DOWRITEIMGS = .true.

! VARS
type(image)                   :: micrograph, mic_shrunken
type(image),      allocatable :: refs(:)
logical,          allocatable :: selected_peak_positions(:)
real,             allocatable :: corrmat(:,:), peak_stats(:,:)
integer,          allocatable :: peak_positions(:,:)
character(len=:), allocatable :: micname
character(len=LONGSTRLEN)     :: boxname
integer                       :: ldim(3), ldim_shrink(3)
integer                       :: ntargets
integer                       :: nrefs, npeaks, npeaks_sel
integer                       :: orig_box, box_shrunken
real                          :: min_rad, max_rad
real                          :: smpd, smpd_shrunken
real                          :: lp, distthr, ndev

contains

    subroutine init_phasecorr_segpicker( micfname, minrad, maxrad, smpd_in, lp_in, distthr_in, ndev_in, dir_out )
        use simple_procimgfile, only :  clip_imgfile
        character(len=*),           intent(in) :: micfname
        real,                       intent(in) :: minrad, maxrad ! in pixels
        real,                       intent(in) :: smpd_in
        real,             optional, intent(in) :: lp_in, distthr_in, ndev_in
        character(len=*), optional, intent(in) :: dir_out
        integer           :: ifoo
        real              :: hp
        allocate(micname,  source=trim(micfname), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('picker;init, 1',alloc_stat)
        boxname = basename( fname_new_ext(micname,'box') )
        if( present(dir_out) )boxname = trim(dir_out)//trim(boxname)
        smpd = smpd_in
        orig_box = round2even(4.*maxrad)
        min_rad = minrad/(PICKER_SHRINK_HERE) !in pixels, shrunken dimensions
        max_rad = maxrad/(PICKER_SHRINK_HERE) !in pixels, shrunken dimensions
        lp   = 20.0
        if( present(lp_in) ) lp = lp_in
        ndev = 2.0
        if( present(ndev_in)) ndev = ndev_in
        ! read micrograph
        call find_ldim_nptcls(micname, ldim, ifoo)
        call micrograph%new(ldim, smpd)
        call micrograph%read(micname)
        box_shrunken = round2even(orig_box/PICKER_SHRINK_HERE)
        ! modify according to PICKER_SHRINK_HERE
        ldim_shrink(1)        = round2even(real(ldim(1))/PICKER_SHRINK_HERE)
        ldim_shrink(2)        = round2even(real(ldim(2))/PICKER_SHRINK_HERE)
        ldim_shrink(3)        = 1
        smpd_shrunken         = PICKER_SHRINK_HERE*smpd
        distthr               = BOXFRAC*real(box_shrunken) !In shrunken dimensions
        if( present(distthr_in) ) distthr = distthr_in/PICKER_SHRINK_HERE
        ! generate references
        call generate_gaussian_refs(refs,min_rad,max_rad)
        ! pre-process micrograph
        call micrograph%fft()
        call mic_shrunken%new(ldim_shrink, smpd_shrunken)
        call mic_shrunken%set_ft(.true.)
        call micrograph%clip(mic_shrunken)
        hp = real(ldim_shrink(1) / 2) * smpd_shrunken
        call mic_shrunken%bp(hp, lp)
        call mic_shrunken%ifft()
    contains
        subroutine generate_gaussian_refs(refs,min_rad,max_rad)
            type(image), allocatable, intent(inout) :: refs(:)
            real,   intent(in) :: min_rad, max_rad
            real,  allocatable :: rmat_out(:,:,:) !to use rtsq_serial
            integer, parameter :: STEP = 2
            real    :: sigma_x, sigma_y
            real    :: rad_x, rad_y
            integer :: i, j
            logical :: rotate_ref
            nrefs = 0
            do i = 0, STEP-1
                rad_x   = (min_rad + real(i*STEP))
                sigma_x = rad_x/(sqrt(2.*log(2.)))           !FWHM = 2*min_rad (see Wiki formula)
                if(rad_x < max_rad+STEP) then
                    do j = 0, STEP-1
                        rad_y   = (min_rad + real(j*STEP))
                        sigma_y = rad_y/(sqrt(2.*log(2.)))   !FWHM = 2*min_rad (see Wiki formula)
                        if(rad_y < max_rad+STEP) then
                            nrefs = nrefs + 1
                            if(rad_y/rad_x > 1.5) rotate_ref = .true.
                            if(rotate_ref) nrefs = nrefs + 2 !total 3 references
                        endif
                        rotate_ref = .false.                 !restore
                    enddo
                endif
            enddo
            if(allocated(refs)) deallocate(refs)
            allocate( refs(nrefs) )
            allocate(rmat_out(ldim_shrink(1),ldim_shrink(2),1), source = 0.)
            nrefs = 0
            do i = 0, STEP-1
                rad_x = (min_rad + real(i*STEP))
                sigma_x = rad_x/(sqrt(2.*log(2.))) ! FWHM = 2*min_rad (see Wiki formula)
                if(rad_x < max_rad+STEP) then      !+STEP tollerance
                    do j = 0, STEP-1
                        rad_y = (min_rad + real(j*STEP))
                        sigma_y = rad_y/(sqrt(2.*log(2.))) !FWHM = 2*min_rad (see Wiki formula)
                        if(rad_y < max_rad+STEP) then
                             if(rad_y/rad_x > 1.5) rotate_ref = .true.
                             nrefs = nrefs + 1
                             call refs(nrefs)%new(ldim_shrink,smpd_shrunken)
                             call refs(nrefs)%gauimg2D(sigma_x,sigma_y,(rad_x+rad_y)/2.+5.) !five pxls cutoff
                             call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                             if(rotate_ref) then
                                 call refs(nrefs)%rtsq_serial( 45., 0., 0., rmat_out )
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_shrink,smpd_shrunken)
                                 call refs(nrefs)%set_rmat(rmat_out)
                                 call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                                 call refs(nrefs)%rtsq_serial( 90., 0., 0., rmat_out )
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_shrink,smpd_shrunken)
                                 call refs(nrefs)%set_rmat(rmat_out)
                                 call refs(nrefs)%write(PATH_HERE//'_GaussianReference.mrc',nrefs)
                             endif
                        endif
                        rotate_ref = .false. !restore
                    enddo
                endif
            enddo
            deallocate(rmat_out)
        end subroutine generate_gaussian_refs
    end subroutine init_phasecorr_segpicker

    subroutine exec_phasecorr_segpicker( boxname_out, nptcls_out )
        character(len=LONGSTRLEN), intent(out) :: boxname_out
        integer,                   intent(out) :: nptcls_out
        call extract_peaks
        call distance_filter
        call gather_stats
        call one_cluster_clustering
        nptcls_out = count(selected_peak_positions)
        ! bring back coordinates to original sampling
        peak_positions = nint(PICKER_SHRINK_HERE*(real(peak_positions)))-orig_box/2
        call write_boxfile
        ! returns absolute path
        call make_relativepath(CWD_GLOB, boxname, boxname_out)
    end subroutine exec_phasecorr_segpicker

    subroutine extract_peaks
        type(image) :: mask_img
        integer, allocatable :: labels(:), target_positions(:,:)
        logical, allocatable :: mask(:,:)
        real,    allocatable :: target_corrs(:)
        real,    pointer     :: rmat_phasecorr(:,:,:)
        real    :: ave, sdev, maxv, minv
        real    :: means(2)
        integer :: xind, yind, alloc_stat, i
        write(logfhandle,'(a)') '>>> EXTRACTING PEAKS'
        call gen_phase_correlation(mic_shrunken,mask_img)
        call mic_shrunken%stats( ave=ave, sdev=sdev, maxv=maxv, minv=minv,mskimg=mask_img)
        call mask_img%kill
        call mic_shrunken%get_rmat_ptr(rmat_phasecorr)
        allocate(corrmat(1:ldim_shrink(1),1:ldim_shrink(2)))
        corrmat(1:ldim_shrink(1),1:ldim_shrink(2)) = rmat_phasecorr(1:ldim_shrink(1),1:ldim_shrink(2),1)
        call mic_shrunken%bin(ave+.8*sdev)
        if(DOWRITEIMGS) call mic_shrunken%write(PATH_HERE//basename(trim(micname))//'_shrunken_bin.mrc')
        allocate(mask(1:ldim_shrink(1), 1:ldim_shrink(2)), source = .false.)
        ntargets = 0
        !$omp parallel do collapse(2) default(shared) private(xind,yind) proc_bind(close) schedule(static) reduction(+:ntargets)
        do xind=1,ldim_shrink(1),PICKER_OFFSET
            do yind=1,ldim_shrink(2),PICKER_OFFSET
                if(rmat_phasecorr(xind,yind,1) > 0.5) then
                    ntargets = ntargets + 1
                    mask(xind,yind) = .true.
                endif
            enddo
        enddo
        !$omp end parallel do
        allocate( target_corrs(ntargets),target_positions(ntargets,2))
        ntargets = 0
        do xind=1,ldim_shrink(1),PICKER_OFFSET
            do yind=1,ldim_shrink(2),PICKER_OFFSET
                if(mask(xind,yind)) then
                    ntargets = ntargets + 1
                    target_positions(ntargets,:) = [xind,yind]
                endif
            enddo
        enddo
        target_corrs = pack(corrmat,mask)
        call sortmeans(target_corrs, MAXKMIT, means, labels)  !TO IMPROVE
        npeaks = count(labels == 2)
        allocate( peak_positions(npeaks,2),  stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk( 'In: simple_picker :: gen_corr_peaks, 2',alloc_stat)
        peak_positions = 0
        npeaks = 0
        do i=1,ntargets
            if( labels(i) == 2 )then
                npeaks = npeaks + 1
                peak_positions(npeaks,:) = target_positions(i,:)
            endif
        end do
    contains
        ! Reference generation and Phase Correlation calculation
        ! FORMULA: phasecorr = ifft2(fft2(field).*conj(fft2(reference)));
        subroutine gen_phase_correlation(field,mask)
            type(image),         intent(inout) :: field
            type(image),optional,intent(inout) :: mask
            type(image)   :: phasecorr, aux
            real, pointer :: mask_rmat(:,:,:)
            integer :: iref
            call phasecorr%new(ldim_shrink, smpd_shrunken)
            call aux%new(ldim_shrink, smpd_shrunken)
            do iref = 1, nrefs
                aux = field%phase_corr(refs(iref),border=box_shrunken/2,lp=lp) !phase correlation
                if(DOWRITEIMGS)  call aux%write(PATH_HERE//'PhaseCorrs.mrc',iref)
                if(iref > 1) then
                    call max_image(phasecorr,phasecorr,aux) !save in phasecorr the maximum value between previous phasecorr and new phasecorr
                else
                    phasecorr   = aux
                endif
            enddo
            call field%copy(phasecorr)
            call field%neg() !The correlations are inverted because the references are white particles on black backgound
            if(DOWRITEIMGS) call field%write(PATH_HERE//basename(trim(micname))//'MaxValPhaseCorr.mrc')
            if(present(mask)) then
                call mask%new(ldim_shrink, smpd_shrunken)
                call mask%get_rmat_ptr(mask_rmat)
                mask_rmat(box_shrunken/2+1:ldim_shrink(1)-box_shrunken/2,box_shrunken/2+1:ldim_shrink(2)-box_shrunken/2,1)=1.
            endif
        end subroutine gen_phase_correlation

        ! This subroutine creates an image in which each pixel value
        ! is the max value between the gray level in img1 and img2
        subroutine max_image(img,img1,img2)
            type(image), intent(inout) :: img !output image
            type(image), intent(in)    :: img1, img2
            real, pointer :: rmat(:,:,:)
            real, pointer :: rmat1(:,:,:), rmat2(:,:,:) !matrices for img1 and img2
            integer :: i, j
            call img%get_rmat_ptr(rmat)
            call img1%get_rmat_ptr(rmat1)
            call img2%get_rmat_ptr(rmat2)
            !$omp parallel do collapse(2) schedule(static) default(shared) private(i,j) proc_bind(close)
            do i = 1, ldim_shrink(1)
                do j = 1, ldim_shrink(2)
                    rmat(i,j,1) = max(rmat1(i,j,1), rmat2(i,j,1))
                 enddo
            enddo
            !$omp end parallel do
        end subroutine max_image
    end subroutine extract_peaks

    subroutine distance_filter
        integer :: ipeak, jpeak, ipos(2), jpos(2), loc(1)
        real    :: dist
        logical, allocatable :: mask(:)
        real,    allocatable :: corrs(:)
        write(logfhandle,'(a)') '>>> DISTANCE FILTERING'
        allocate( mask(npeaks), corrs(npeaks), selected_peak_positions(npeaks), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk( 'In: simple_picker :: distance_filter',alloc_stat)
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
        write(logfhandle,'(a,1x,I5)') 'peak positions left after distance filtering: ', npeaks_sel
        deallocate(mask,corrs)
    end subroutine distance_filter

    subroutine gather_stats
        integer           :: ipeak, cnt, istat
        logical           :: outside
        real              :: ave, maxv, minv
        real, allocatable :: spec(:)
        type(image)       :: ptcl_target
        write(logfhandle,'(a)') '>>> GATHERING REMAINING STATS'
        call ptcl_target%new([box_shrunken,box_shrunken,1], smpd_shrunken)
        cnt = 0
        allocate( peak_stats(npeaks_sel,NSTAT) )
        do ipeak=1,npeaks
            if( selected_peak_positions(ipeak) )then
                cnt = cnt + 1
                call mic_shrunken%window_slim(peak_positions(ipeak,:)-box_shrunken,&
                    &box_shrunken, ptcl_target, outside)
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
        call ptcl_target%kill()
        deallocate(spec)
    end subroutine gather_stats

    subroutine one_cluster_clustering
        real, allocatable :: dmat(:,:)
        integer           :: i_median, i, j, cnt, ipeak
        real              :: ddev
        allocate(dmat(npeaks_sel,npeaks_sel), source=0., stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('picker::one_cluster_clustering dmat ',alloc_stat)
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
        write(logfhandle,'(a,1x,I5)') 'peak positions left after one cluster clustering: ', npeaks_sel
        deallocate(dmat)
    end subroutine one_cluster_clustering

    subroutine write_boxfile
        integer :: funit, ipeak,iostat
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(boxname)),iostat=iostat)
        call fileiochk('phasecorr_picker; write_boxfile ', iostat)
        do ipeak=1,npeaks
            if( selected_peak_positions(ipeak) )then
                write(funit,'(I7,I7,I7,I7,I7)') peak_positions(ipeak,1)-1,&
                peak_positions(ipeak,2)-1, orig_box, orig_box, -3
            endif
        end do
        call fclose(funit,errmsg='picker; write_boxfile end')
    end subroutine write_boxfile

    subroutine kill_phasecorr_segpicker
        integer :: iref
        if( allocated(micname) )then
            call micrograph%kill
            call mic_shrunken%kill
            deallocate(selected_peak_positions,corrmat, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('phasecorr_picker kill, 1',alloc_stat)
            deallocate(peak_positions,micname,peak_stats, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('phasecorr_picker kill, 2',alloc_stat)
            do iref=1,nrefs
                call refs(iref)%kill
            end do
            deallocate(refs, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('phasecorr_picker; kill 3',alloc_stat)
        endif
    end subroutine kill_phasecorr_segpicker
end module simple_phasecorr_segpicker
