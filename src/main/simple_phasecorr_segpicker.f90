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
integer,          parameter   :: NSTAT    = 4
integer,          parameter   :: MAXKMIT  = 20
integer,          parameter   :: MIN_NCCS = 5    ! minimum number of connected components to identify after size-filtering
real,             parameter   :: BOXFRAC  = 0.5
logical,          parameter   :: DOWRITEIMGS = .false.
logical,          parameter   :: DOPRINT     = .false.

! VARS
type(image)                   :: micrograph, mic_shrunken
type(image),      allocatable :: refs(:)
logical,          allocatable :: selected_peak_positions(:)
real,             allocatable :: corrmat(:,:), peak_stats(:,:)
integer,          allocatable :: peak_positions(:,:)
character(len=:), allocatable :: micname
character(len=LONGSTRLEN)     :: boxname
character(len=3)              :: elongated  !additional info inputted by the user
integer                       :: ldim(3), ldim_shrink(3)
integer                       :: ntargets
integer                       :: nrefs, npeaks, npeaks_sel
integer                       :: orig_box, box_shrunken
real                          :: min_rad, max_rad, step_sz
real                          :: smpd, smpd_shrunken
real                          :: lp, distthr, ndev

contains

    subroutine init_phasecorr_segpicker( micfname, minrad, maxrad, stepsz, eelongated, smpd_in, lp_in, distthr_in, ndev_in, dir_out )
        use simple_procimgfile, only :  clip_imgfile
        character(len=*),           intent(in) :: micfname
        real,                       intent(in) :: minrad, maxrad, stepsz ! in A
        character(len=3),           intent(in) :: eelongated
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
        min_rad = minrad/(PICKER_SHRINK_HERE*smpd)    !in pixel, shrunken dimensions
        max_rad = maxrad/(PICKER_SHRINK_HERE*smpd)    !in pixel, shrunken dimensions
        step_sz = stepsz/(PICKER_SHRINK_HERE*smpd)    !in pixel, shrunken dimensions
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
        elongated = eelongated
        ! generate references
        if(elongated .eq. 'yes') then
            write(logfhandle, *) 'Elongated references generation'
            call generate_gaussian_refs_elongated(refs,min_rad,max_rad,step_sz)
        else
            call generate_gaussian_refs(refs,min_rad,max_rad,step_sz)
        endif
        ! pre-process micrograph
        call micrograph%fft()
        call mic_shrunken%new(ldim_shrink, smpd_shrunken)
        call mic_shrunken%set_ft(.true.)
        call micrograph%clip(mic_shrunken)
        hp = real(ldim_shrink(1) / 2) * smpd_shrunken
        call mic_shrunken%bp(hp, lp)
        call mic_shrunken%ifft()
        if(DOWRITEIMGS) call mic_shrunken%write('mic_shrunken.mrc')
    contains

        subroutine generate_gaussian_refs_elongated(refs,min_rad,max_rad,step_sz)
            type(image), allocatable, intent(inout) :: refs(:)
            real, intent(in)  :: min_rad, max_rad, step_sz
            real, allocatable :: rmat_out(:,:,:) !to use rtsq_serial
            real    :: sigma_x, sigma_y
            real    :: rad_x, rad_y
            integer :: i, j, max_nrefs
            nrefs = 0
            max_nrefs = ceiling((max_rad-min_rad)/step_sz)*4 ! possibility to rotate 3 times
            do i = 0, max_nrefs-1
                rad_x   = (min_rad + real(i)*step_sz)
                if(rad_x <= max_rad) then
                    do j = i+1, max_nrefs
                        rad_y   = (min_rad + real(j)*step_sz)
                        if(rad_y <= max_rad) then
                            if(rad_y/rad_x >= 2.) then
                              nrefs = nrefs + 4 !total 4 references (original and 3 rotations)
                            endif
                        endif
                    enddo
                endif
            enddo
            if(allocated(refs)) deallocate(refs)
            allocate( refs(nrefs) )
            allocate(rmat_out(ldim_shrink(1),ldim_shrink(2),1), source = 0.)
            nrefs = 0
            do i = 0, max_nrefs-1
                rad_x   = (min_rad + real(i)*step_sz)
                if(rad_x <= max_rad) then
                    do j = i+1, max_nrefs
                        rad_y   = (min_rad + real(j)*step_sz)
                        if(rad_y <= max_rad) then
                             if(rad_y/rad_x >= 2.) then
                                 sigma_x = rad_x/(2.*sqrt(2.*log(2.))) ! FWHM = rad (see Wiki formula)
                                 sigma_y = rad_y/(2.*sqrt(2.*log(2.))) ! FWHM = rad (see Wiki formula)
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_shrink,smpd_shrunken)
                                 call refs(nrefs)%gauimg2D(sigma_x,sigma_y,(rad_x+rad_y)/2.+5.) !five pxls cutoff
                                 call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                                 if(DOPRINT) write(logfhandle,*) 'generating ref with rads: ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                                 call refs(nrefs)%rtsq_serial( 45., 0., 0., rmat_out )
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_shrink,smpd_shrunken)
                                 call refs(nrefs)%set_rmat(rmat_out)
                                 call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                                 if(DOPRINT) write(logfhandle,*) 'rotating 45 degrees ref:  ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                                 call refs(nrefs)%rtsq_serial( 90., 0., 0., rmat_out )
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_shrink,smpd_shrunken)
                                 call refs(nrefs)%set_rmat(rmat_out)
                                 call refs(nrefs)%write(PATH_HERE//'_GaussianReference.mrc',nrefs)
                                 if(DOPRINT) write(logfhandle,*) 'rotating 90 degrees ref:  ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                                 call refs(nrefs)%rtsq_serial( 135., 0., 0., rmat_out )
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_shrink,smpd_shrunken)
                                 call refs(nrefs)%set_rmat(rmat_out)
                                 call refs(nrefs)%write(PATH_HERE//'_GaussianReference.mrc',nrefs)
                                 if(DOPRINT) write(logfhandle,*) 'rotating 135 degrees ref:  ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                              endif
                        endif
                    enddo
                endif
            enddo
            deallocate(rmat_out)
            if(nrefs == 0) then
                write(logfhandle, *) 'Min_rad and max_rad have to have ratio > 2 to generate elongated refs; generate_gaussian_refs_elongated'
                stop
            endif
            if(DOPRINT) write(logfhandle,*) 'Refs generation completed'
        end subroutine generate_gaussian_refs_elongated

        subroutine generate_gaussian_refs(refs,min_rad,max_rad,step_sz)
            type(image), allocatable, intent(inout) :: refs(:)
            real, intent(in)  :: min_rad, max_rad, step_sz
            real, allocatable :: rmat_out(:,:,:) !to use rtsq_serial
            real    :: sigma_x, sigma_y
            real    :: rad_x, rad_y
            integer :: i, j, max_nrefs
            logical :: rotate_ref
            nrefs = 0
            max_nrefs = ceiling((max_rad-min_rad)/step_sz)
            do i = 0, max_nrefs
                rad_x   = (min_rad + real(i)*step_sz)
                if(rad_x <= max_rad) then
                    do j = 0, max_nrefs
                        rad_y   = (min_rad + real(j)*step_sz)
                        if(rad_y <= max_rad) then
                            nrefs = nrefs + 1
                            if(rad_y/rad_x >= 1.5) then
                              rotate_ref = .true.
                              nrefs = nrefs + 2 !total 3 references
                            endif
                        endif
                        rotate_ref = .false.                 !restore
                    enddo
                endif
            enddo
            if(allocated(refs)) deallocate(refs)
            allocate( refs(nrefs) )
            allocate(rmat_out(ldim_shrink(1),ldim_shrink(2),1), source = 0.)
            nrefs = 0
            do i = 0, max_nrefs
                rad_x = (min_rad + real(i)*step_sz)
                if(rad_x <= max_rad) then
                    do j = 0, max_nrefs
                        rad_y = (min_rad + real(j)*step_sz)
                        if(rad_y <= max_rad) then
                             if(rad_y/rad_x >= 1.5) rotate_ref = .true.
                             sigma_x = 0.5*rad_x/(sqrt(log(2.))) !FWHM = half rad (see Wiki formula)
                             sigma_y = 0.5*rad_y/(sqrt(log(2.))) !FWHM = half rad (see Wiki formula)
                             nrefs = nrefs + 1
                             call refs(nrefs)%new(ldim_shrink,smpd_shrunken)
                             call refs(nrefs)%gauimg2D(sigma_x,sigma_y,(rad_x+rad_y)/2.+5.) !five pxls cutoff
                             call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                             if(DOPRINT) write(logfhandle,*)  'generating ref with rads: ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                             if(rotate_ref) then
                                 call refs(nrefs)%rtsq_serial( 45., 0., 0., rmat_out )
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_shrink,smpd_shrunken)
                                 call refs(nrefs)%set_rmat(rmat_out)
                                 call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                                 if(DOPRINT) write(logfhandle,*)  'rotating 45 degrees ref:  ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                                 call refs(nrefs)%rtsq_serial( 90., 0., 0., rmat_out )
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_shrink,smpd_shrunken)
                                 call refs(nrefs)%set_rmat(rmat_out)
                                 call refs(nrefs)%write(PATH_HERE//'_GaussianReference.mrc',nrefs)
                                 if(DOPRINT) write(logfhandle,*) 'rotating 90 degrees ref:  ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                             endif
                        endif
                        rotate_ref = .false. !restore
                    enddo
                endif
            enddo
            deallocate(rmat_out)
            if(DOPRINT) write(logfhandle,*)  'Refs generation completed'
        end subroutine generate_gaussian_refs
    end subroutine init_phasecorr_segpicker

    subroutine exec_phasecorr_segpicker( boxname_out, nptcls_out, center )
        use simple_timer
        character(len=LONGSTRLEN), intent(out) :: boxname_out
        integer,                   intent(out) :: nptcls_out
        character(len=3),          intent(in)  :: center
        integer(timer_int_kind) :: time
        time = tic()
        call extract_peaks
        print *, 'time of extract peaks: ', toc(time)
        if(center .eq. 'yes') then
            call distance_filter_and_center
        else
            call distance_filter
        endif
        print *, 'time of distance filtering: ', toc(time)
        call gather_stats
        print *, 'time of gather stats: ', toc(time)
        call one_cluster_clustering
        print *, 'time of one cluster clustering: ', toc(time)
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
        call mic_shrunken%stats(ave=ave, sdev=sdev, maxv=maxv, minv=minv,mskimg=mask_img)
        call mask_img%kill
        call mic_shrunken%get_rmat_ptr(rmat_phasecorr)
        allocate(corrmat(1:ldim_shrink(1),1:ldim_shrink(2)))
        corrmat(1:ldim_shrink(1),1:ldim_shrink(2)) = rmat_phasecorr(1:ldim_shrink(1),1:ldim_shrink(2),1)
        call mic_shrunken%binarize(ave+.8*sdev)
        rmat_phasecorr(1:box_shrunken/2,:,1) = 0. !set to zero the borders
        rmat_phasecorr(ldim_shrink(1)-box_shrunken/2:ldim_shrink(1),:,1) = 0. !set to zero the borders
        rmat_phasecorr(:,1:box_shrunken/2,1) = 0. !set to zero the borders
        rmat_phasecorr(:,ldim_shrink(2)-box_shrunken/2:ldim_shrink(2),1) = 0. !set to zero the borders
        if(DOWRITEIMGS) call mic_shrunken%write(PATH_HERE//basename(trim(micname))//'_shrunken_bin.mrc')
        allocate(mask(1:ldim_shrink(1),1:ldim_shrink(2)), source = .false.)
        ntargets = 0
        !$omp parallel do collapse(2) default(shared) private(xind,yind) proc_bind(close) schedule(static) reduction(+:ntargets)
        do xind=1,ldim_shrink(1)
            do yind=1,ldim_shrink(2)
                if(rmat_phasecorr(xind,yind,1) > 0.5) then
                    ntargets = ntargets + 1
                    mask(xind,yind) = .true.
                endif
            enddo
        enddo
        !$omp end parallel do
        allocate( target_corrs(ntargets),target_positions(ntargets,2))
        ntargets = 0
        do xind=1,ldim_shrink(1)
            do yind=1,ldim_shrink(2)
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
            call aux%zero_and_flag_ft()
            call phasecorr%zero_and_flag_ft()
            call field%fft
            do iref = 1, nrefs
                call refs(iref)%fft
                call aux%zero_and_flag_ft
                call field%phase_corr(refs(iref),aux,lp,border=box_shrunken/2) !phase correlation
                if(DOWRITEIMGS) call aux%write(PATH_HERE//'PhaseCorrs.mrc',iref)
                if(iref > 1) then
                    call max_image(phasecorr,phasecorr,aux) !save in phasecorr the maximum value between previous phasecorr and new phasecorr
                    call phasecorr%write('Maxphasecorr.mrc', iref)
                else
                    phasecorr   = aux
                    call phasecorr%write('Maxphasecorr.mrc', iref)
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
            !omp parallel do collapse(2) schedule(static) default(shared) private(i,j) proc_bind(close)
            do i = 1, ldim_shrink(1)
                do j = 1, ldim_shrink(2)
                    rmat(i,j,1) = max(rmat1(i,j,1), rmat2(i,j,1))
                 enddo
            enddo
            !omp end parallel do
        end subroutine max_image
    end subroutine extract_peaks

    subroutine distance_filter_and_center
        integer :: ipeak, jpeak, ipos(2), jpos(2), loc(1), cnt_peaks
        real    :: dist
        logical, allocatable :: mask(:)
        real,    allocatable :: corrs(:)
        write(logfhandle,'(a)') '>>> DISTANCE FILTERING'
        allocate( mask(npeaks), corrs(npeaks), selected_peak_positions(npeaks), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk( 'In: simple_picker :: distance_filter',alloc_stat)
        selected_peak_positions = .true.
        call center_particles(cnt_peaks)
        print *, 'after centering and cc elim we have ', cnt_peaks, 'peaks'
        print *, 'should be the same as ', count(selected_peak_positions)
        ! HEREEEE TO PUT BACK OPENMP
        !omp do collapse(2) schedule(static) default(shared) private(ipeak,jpeak) proc_bind(close)
        do ipeak = 1, cnt_peaks-1            !fix one coord
            do jpeak = ipeak+1, cnt_peaks       !fix another coord to compare
                if(selected_peak_positions(ipeak) .and. selected_peak_positions(jpeak)) then !not compare twice ,and if the particles haven t been deleted yet
                    if( euclid(real(peak_positions(ipeak,:)),real(peak_positions(jpeak,:))) <= distthr) then
                        selected_peak_positions(ipeak) = .false.
                        selected_peak_positions(jpeak) = .false.
                    endif
                endif
            enddo
        enddo
        !omp end do
        print *, 'after distance threshold we have', count(selected_peak_positions), 'peaks'
        npeaks_sel = count(selected_peak_positions)
        write(logfhandle,'(a,1x,I5)') 'peak positions left after distance filtering: ', npeaks_sel
        deallocate(mask,corrs)
    end subroutine distance_filter_and_center

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

    ! This subroutine centers the picked particles. The center
    ! is defined as the center of mass of the
    ! connected component to which the particle belongs.
    subroutine center_particles(cnt_peaks)
      use simple_binimage, only : binimage
      integer, intent(out) :: cnt_peaks
      type(binimage)       :: micrograph_bin,micrograph_cc
      logical, allocatable :: center(:) ! not to center more than once the same cc
      integer, allocatable :: imat_cc(:,:,:), imat_aux(:,:,:)
      integer :: ipeak, cc, cnt, i, j
      print *, 'Initialization centering'
      call micrograph_bin%new_bimg(ldim_shrink, smpd_shrunken)
      ! I need a copy to have a binimage instead of an image
      call micrograph_bin%copy(mic_shrunken) ! mic_shrunken is already binary
      call micrograph_bin%find_ccs(micrograph_cc)
      ! 6) cc filtering, based on the shape of the particle
      call micrograph_cc%polish_ccs([min_rad,max_rad], circular=' no',elongated=elongated, min_nccs = MIN_NCCS)
      call micrograph_cc%write_bimg('CcsElimin.mrc')
      call micrograph_cc%get_imat(imat_cc)
      allocate(imat_aux(ldim_shrink(1), ldim_shrink(2),1), source = 0)
      allocate(center(maxval(imat_cc)), source = .true.)
      cnt_peaks = 0 ! count the number of peaks after centering and cc filtering
      !$omp parallel do default(shared) private(ipeak,cc,cnt,i,j) proc_bind(close) schedule(static)
      do ipeak = 1,npeaks
          cc = imat_cc(peak_positions(ipeak,1), peak_positions(ipeak,2),1)
          if(cc < 0.5) then
             ! remove peaks that correspond to filtered ccs.
              selected_peak_positions(ipeak) = .false.
           elseif(center(cc)) then
              cnt_peaks = cnt_peaks + 1
              center(cc) = .false.
              ! center peaks, center of mass of the cc
              peak_positions(ipeak,:) = 0 ! reset
              cnt = 0
              do i = 1, ldim_shrink(1)
                    do j = 1, ldim_shrink(2)
                        if(imat_cc(i,j,1) == cc) then
                            cnt = cnt + 1
                            peak_positions(ipeak,1:2) = peak_positions(ipeak,1:2) + [i,j]
                        endif
                    enddo
                enddo
                peak_positions(ipeak,1:2) = peak_positions(ipeak,1:2)/cnt
              else
                selected_peak_positions(ipeak) = .false.
           endif
      enddo
      !$omp end parallel do
      call micrograph_bin%kill_bimg
      call micrograph_cc%kill_bimg
    end subroutine center_particles

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
                write(funit,'(I7,I7,I7,I7,I7)') peak_positions(ipeak,1),&
                peak_positions(ipeak,2), orig_box, orig_box, -3
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
