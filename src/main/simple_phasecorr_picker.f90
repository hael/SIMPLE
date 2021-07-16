! particle picker
module simple_phasecorr_picker
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,        only: image
implicit none

public :: init_phasecorr_picker_refs, init_phasecorr_picker_gauss, exec_phasecorr_picker, kill_phasecorr_picker
private

! PEAK STATS INDICES
integer,          parameter   :: SDEV     = 1
integer,          parameter   :: DYNRANGE = 2
integer,          parameter   :: SSCORE   = 3
! OTHER PARAMS
integer,          parameter   :: NSTAT   = 3
integer,          parameter   :: MAXKMIT = 20
real,             parameter   :: BOXFRAC = 0.5
logical,          parameter   :: DOWRITEIMGS = .false., DEBUG_HERE = .false.
integer,          parameter   :: PICKER_OFFSET_HERE = 3, OFFSET_HWIN = 1

! VARS
type(image)                   :: micrograph, mic_shrunken, mic_saved
type(image),      allocatable :: refs(:)
logical,          allocatable :: selected_peak_positions(:)
real,             allocatable :: corrmat(:,:), peak_stats(:,:)
integer,          allocatable :: peak_positions(:,:)
character(len=:), allocatable :: micname, refsname
character(len=LONGSTRLEN)     :: boxname
integer                       :: ldim(3), ldim_refs(3), ldim_shrink(3)
integer                       :: ntargets
integer                       :: nrefs, nmax, nmax_sel, orig_box
real                          :: smpd, smpd_shrunken
real                          :: msk, hp,lp, distthr, ndev
real                          :: min_rad, max_rad, step_sz
character(len=3)              :: elongated  !additional info inputted by the user
logical                       :: template_based = .false.


contains

    subroutine init_phasecorr_picker_refs( micfname, refsfname, smpd_in, lp_in, distthr_in, ndev_in, dir_out )
        use simple_procimgfile, only :  clip_imgfile
        character(len=*),           intent(in) :: micfname, refsfname
        real,                       intent(in) :: smpd_in
        real,             optional, intent(in) :: lp_in, distthr_in, ndev_in
        character(len=*), optional, intent(in) :: dir_out
        logical, allocatable :: lmsk(:,:,:)
        type(image)          :: refimg, mskimg
        integer              :: ifoo, iref
        real                 :: sdev_noise
        template_based = .true.
        allocate(micname,  source=trim(micfname), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('picker;init, 1',alloc_stat)
        allocate(refsname, source=trim(refsfname), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('picker;init, 2')
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
        ! modify according to PICKER_SHRINK
        ldim_refs(1)          = round2even(real(ldim_refs(1))/PICKER_SHRINK)
        ldim_refs(2)          = round2even(real(ldim_refs(2))/PICKER_SHRINK)
        ldim_refs(3)          = 1
        ldim_shrink(1)        = round2even(real(ldim(1))/PICKER_SHRINK)
        ldim_shrink(2)        = round2even(real(ldim(2))/PICKER_SHRINK)
        ldim_shrink(3)        = 1
        smpd_shrunken         = PICKER_SHRINK*smpd
        msk                   = real(ldim_refs(1)/2-5)
        msk                   = max(msk, real(ldim_refs(1)/2-2)) ! for tiny particles
        distthr               = BOXFRAC*real(ldim_refs(1))
        if( present(distthr_in) ) distthr = distthr_in / smpd_shrunken ! pixels
        ! read and shrink references
        allocate( refs(nrefs), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk( "In: simple_picker :: init_picker, refs etc. ",alloc_stat)
        call mskimg%disc([orig_box,orig_box,1], smpd, msk*real(orig_box)/real(ldim_refs(1)), lmsk)
        call mskimg%kill
        do iref=1,nrefs
            call refs(iref)%new(ldim_refs, smpd_shrunken)
            call refimg%new([orig_box,orig_box,1], smpd)
            call refimg%read(refsname, iref)
            call refimg%noise_norm(lmsk, sdev_noise)
            call refimg%mask(msk*real(orig_box)/real(ldim_refs(1)), 'soft')
            call refimg%fft()
            call refimg%clip(refs(iref))
            call refs(iref)%ifft()
            call refs(iref)%mask(msk, 'soft')
        end do
        call refimg%kill
        ! pre-process micrograph
        call micrograph%fft()
        call mic_shrunken%new(ldim_shrink, smpd_shrunken)
        call mic_shrunken%set_ft(.true.)
        call micrograph%clip(mic_shrunken)
        hp = real(ldim_shrink(1) / 2) * smpd_shrunken
        call mic_shrunken%fft
        mic_saved = mic_shrunken
        call mic_shrunken%bp(hp, lp)
        call mic_saved%ifft
      end subroutine init_phasecorr_picker_refs

      subroutine init_phasecorr_picker_gauss( micfname, minrad, maxrad, stepsz, eelongated, smpd_in, lp_in, distthr_in, ndev_in, dir_out )
          use simple_procimgfile, only :  clip_imgfile
          character(len=*),           intent(in) :: micfname
          real,                       intent(in) :: minrad, maxrad, stepsz ! in A
          character(len=3),           intent(in) :: eelongated
          real,                       intent(in) :: smpd_in
          real,             optional, intent(in) :: lp_in, distthr_in, ndev_in
          character(len=*), optional, intent(in) :: dir_out
          integer           :: ifoo
          real              :: hp
          template_based = .false.  ! reset to be sure
          allocate(micname,  source=trim(micfname), stat=alloc_stat)
          if(alloc_stat.ne.0)call allocchk('picker;init, 1',alloc_stat)
          boxname = basename( fname_new_ext(micname,'box') )
          if( present(dir_out) )boxname = trim(dir_out)//trim(boxname)
          smpd = smpd_in
          orig_box = round2even(2.7*maxrad/smpd)    !in pixel, original dimensions
          min_rad  = minrad/(PICKER_SHRINK*smpd)    !in pixel, shrunken dimensions
          max_rad  = maxrad/(PICKER_SHRINK*smpd)    !in pixel, shrunken dimensions
          step_sz  = stepsz/(PICKER_SHRINK*smpd)    !in pixel, shrunken dimensions
          lp   = 20.0
          if( present(lp_in) ) lp = lp_in
          ndev = 2.0
          if( present(ndev_in)) ndev = ndev_in
          ! read micrograph
          call find_ldim_nptcls(micname, ldim, ifoo)
          call micrograph%new(ldim, smpd)
          call micrograph%read(micname)
          ldim_refs(1:2)   = round2even(orig_box/PICKER_SHRINK)
          ldim_refs(3)     = 1
          msk              = real(ldim_refs(1)/2-5)
          msk              = max(msk, real(ldim_refs(1)/2-2)) ! for tiny particles
          ! modify according to PICKER_SHRINK_HERE
          ldim_shrink(1)   = round2even(real(ldim(1))/PICKER_SHRINK)
          ldim_shrink(2)   = round2even(real(ldim(2))/PICKER_SHRINK)
          ldim_shrink(3)   = 1
          smpd_shrunken    = PICKER_SHRINK*smpd
          distthr          = BOXFRAC*real(ldim_refs(1)) !In shrunken dimensions
          if( present(distthr_in) ) distthr = distthr_in/PICKER_SHRINK/smpd_shrunken ! pixels
          elongated = eelongated
          ! generate references
          if(elongated .eq. 'yes') then
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
          call mic_shrunken%fft
          mic_saved = mic_shrunken
          call mic_shrunken%bp(hp, lp)
          call mic_saved%ifft

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
            allocate(rmat_out(ldim_refs(1),ldim_refs(2),1), source = 0.)
            nrefs = 0
            do i = 0, max_nrefs-1
                rad_x   = (min_rad + real(i)*step_sz)
                if(rad_x <= max_rad) then
                    do j = i+1, max_nrefs
                        rad_y   = (min_rad + real(j)*step_sz)
                        if(rad_y <= max_rad) then
                             if(rad_y/rad_x >= 2.) then
                                 sigma_x = rad_x/2.!(2.*sqrt(2.*log(2.))) ! FWHM = rad (see Wiki formula)
                                 sigma_y = rad_y/2.!(2.*sqrt(2.*log(2.))) ! FWHM = rad (see Wiki formula)
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_refs,smpd_shrunken)
                                 call refs(nrefs)%gauimg2D(sigma_x,sigma_y)! cutoff is hard mask, so we do soft mask afterwards
                                 call refs(nrefs)%mask(max_rad,'soft',backgr=0.)
                                 if(DOWRITEIMGS) call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                                 if(DEBUG_HERE)  write(logfhandle,*) 'generating ref with rads: ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                                 call refs(nrefs)%rtsq_serial( 45., 0., 0., rmat_out )
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_refs,smpd_shrunken)
                                 call refs(nrefs)%set_rmat(rmat_out)
                                 if(DOWRITEIMGS) call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                                 if(DEBUG_HERE) write(logfhandle,*) 'rotating 45 degrees ref:  ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                                 call refs(nrefs)%rtsq_serial( 90., 0., 0., rmat_out )
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_refs,smpd_shrunken)
                                 call refs(nrefs)%set_rmat(rmat_out)
                                 if(DOWRITEIMGS) call refs(nrefs)%write(PATH_HERE//'_GaussianReference.mrc',nrefs)
                                 if(DEBUG_HERE) write(logfhandle,*) 'rotating 90 degrees ref:  ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                                 call refs(nrefs)%rtsq_serial( 135., 0., 0., rmat_out )
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_refs,smpd_shrunken)
                                 call refs(nrefs)%set_rmat(rmat_out)
                                 if(DOWRITEIMGS) call refs(nrefs)%write(PATH_HERE//'_GaussianReference.mrc',nrefs)
                                 if(DEBUG_HERE) write(logfhandle,*) 'rotating 135 degrees ref:  ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
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
            if(DEBUG_HERE) write(logfhandle,*) 'Elongated refs generation completed'
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
            allocate(rmat_out(ldim_refs(1),ldim_refs(2),1), source = 0.)
            nrefs = 0
            do i = 0, max_nrefs
                rad_x = (min_rad + real(i)*step_sz)
                if(rad_x <= max_rad) then
                    do j = 0, max_nrefs
                        rad_y = (min_rad + real(j)*step_sz)
                        if(rad_y <= max_rad) then
                             if(rad_y/rad_x >= 1.5) rotate_ref = .true.
                             sigma_x = rad_x/2.!(2.*sqrt(2.*log(2.))) ! FWHM = rad (see Wiki formula)
                             sigma_y = rad_y/2.!(2.*sqrt(2.*log(2.))) ! FWHM = rad (see Wiki formula)
                             nrefs = nrefs + 1
                             call refs(nrefs)%new(ldim_refs,smpd_shrunken)
                             call refs(nrefs)%gauimg2D(sigma_x,sigma_y)!,(rad_x+rad_y)/2.+5.) !five pxls cutoff
                             call refs(nrefs)%mask(max_rad,'soft',backgr=0.)
                             if(DOWRITEIMGS) call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                             if(DEBUG_HERE)  write(logfhandle,*)  'generating ref with rads: ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                             if(rotate_ref) then
                                 call refs(nrefs)%rtsq_serial( 45., 0., 0., rmat_out )
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_refs,smpd_shrunken)
                                 call refs(nrefs)%set_rmat(rmat_out)
                                 if(DOWRITEIMGS) call refs(nrefs)%write('_GaussianReference.mrc',nrefs)
                                 if(DEBUG_HERE) write(logfhandle,*)  'rotating 45 degrees ref:  ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                                 call refs(nrefs)%rtsq_serial( 90., 0., 0., rmat_out )
                                 nrefs = nrefs + 1
                                 call refs(nrefs)%new(ldim_refs,smpd_shrunken)
                                 call refs(nrefs)%set_rmat(rmat_out)
                                 if(DOWRITEIMGS) call refs(nrefs)%write(PATH_HERE//'_GaussianReference.mrc',nrefs)
                                 if(DEBUG_HERE)  write(logfhandle,*) 'rotating 90 degrees ref:  ', rad_x*smpd_shrunken, rad_y*smpd_shrunken, ' A'
                             endif
                        endif
                        rotate_ref = .false. !restore
                    enddo
                endif
            enddo
            deallocate(rmat_out)
            if(DEBUG_HERE) write(logfhandle,*)  'Refs generation completed'
        end subroutine generate_gaussian_refs
    end subroutine init_phasecorr_picker_gauss

    subroutine exec_phasecorr_picker( boxname_out, nptcls_out )
        character(len=LONGSTRLEN), intent(out) :: boxname_out
        integer,                   intent(out) :: nptcls_out
        real, pointer :: prmat(:,:,:)
        integer       :: i
        real          :: maxv
        call extract_peaks
        call distance_filter
        call gather_stats
        call one_cluster_clustering
        nptcls_out = count(selected_peak_positions)
        if(DOWRITEIMGS) then
            call mic_saved%get_rmat_ptr(prmat)
            maxv = 5.*maxval(prmat)
            do i=1,size(selected_peak_positions)
                if(.not.selected_peak_positions(i))cycle
                call mic_saved%set([peak_positions(i,1)-1,peak_positions(i,2)-1,1],maxv)
                call mic_saved%set([peak_positions(i,1)-1,peak_positions(i,2),1],maxv)
                call mic_saved%set([peak_positions(i,1)-1,peak_positions(i,2)+1,1],maxv)
                call mic_saved%set([peak_positions(i,1),  peak_positions(i,2)-1,1],maxv)
                call mic_saved%set([peak_positions(i,1),  peak_positions(i,2),1],maxv)
                call mic_saved%set([peak_positions(i,1),  peak_positions(i,2)+1,1],maxv)
                call mic_saved%set([peak_positions(i,1)+1,peak_positions(i,2)-1,1],maxv)
                call mic_saved%set([peak_positions(i,1)+1,peak_positions(i,2),1],maxv)
                call mic_saved%set([peak_positions(i,1)+1,peak_positions(i,2)+1,1],maxv)
            enddo
            call mic_saved%write(PATH_HERE//basename(trim(micname))//'_picked.mrc')
        endif
        ! bring back coordinates to original sampling
        peak_positions = nint(PICKER_SHRINK*(real(peak_positions)))-orig_box/2
        call write_boxfile
        ! returns absolute path
        call make_relativepath(CWD_GLOB, boxname, boxname_out)
    end subroutine exec_phasecorr_picker

    subroutine extract_peaks
        type(image)              :: mask_img, sdevimg, circ_mask, tmp_mic
        type(image), allocatable :: ptcls(:)
        real,    pointer         :: rmat_phasecorr(:,:,:)
        integer, allocatable     :: labels(:), target_positions(:,:),ref_inds(:,:)
        real,    allocatable     :: target_corrs(:)
        logical, allocatable     :: mask(:,:), l_mask(:,:,:)
        real                     :: means(2), ave, sdev, maxv, minv, msk_shrunken
        integer                  :: xind, yind, alloc_stat, i, border, j, l, r, u, d, iref, ithr
        logical :: outside
        write(logfhandle,'(a)') '>>> EXTRACTING PEAKS'
        write(logfhandle,'(a)') '>>> FOURIER CORRELATIONS'
        allocate(ref_inds(1:ldim_shrink(1), 1:ldim_shrink(2)), source=0)
        call gen_phase_correlation(mic_shrunken,mask_img)
        call mic_shrunken%stats( ave=ave, sdev=sdev, maxv=maxv, minv=minv,mskimg=mask_img)
        call mic_shrunken%get_rmat_ptr(rmat_phasecorr)
        allocate(corrmat(1:ldim_shrink(1),1:ldim_shrink(2)))
        corrmat(1:ldim_shrink(1),1:ldim_shrink(2)) = rmat_phasecorr(1:ldim_shrink(1),1:ldim_shrink(2),1)
        !$omp parallel do default(shared) private(xind,yind,l,r,u,d) proc_bind(close) schedule(static)
        do xind=1,ldim_shrink(1),PICKER_OFFSET_HERE
            l = max(1,xind-OFFSET_HWIN)
            r = min(xind+ OFFSET_HWIN,ldim_shrink(1))
            do yind=1,ldim_shrink(2),PICKER_OFFSET_HERE
                u = max(1,yind- OFFSET_HWIN)
                d = min(yind+OFFSET_HWIN,ldim_shrink(2))
                corrmat(xind,yind) = maxval(rmat_phasecorr(l:r,u:d,1))
            enddo
        enddo
        !$omp end parallel do
        write(logfhandle,'(a)') '>>> BINARIZATION'
        call mic_shrunken%binarize(ave+.8*sdev)
        border = max(ldim_refs(1)/2,ldim_refs(2)/2)
        rmat_phasecorr(1:border,:,1) = 0. !set to zero the borders
        rmat_phasecorr(ldim_shrink(1)-border:ldim_shrink(1),:,1) = 0. !set to zero the borders
        rmat_phasecorr(:,1:border,1) = 0. !set to zero the borders
        rmat_phasecorr(:,ldim_shrink(2)-border:ldim_shrink(2),1) = 0. !set to zero the borders
        if(DOWRITEIMGS) call mic_shrunken%write(PATH_HERE//basename(trim(micname))//'_shrunken_bin.mrc')
        allocate(mask(1:ldim_shrink(1), 1:ldim_shrink(2)), source = .false.)
        ! Select initial peaks
        ntargets = 0
        !$omp parallel do collapse(2) default(shared) private(xind,yind) proc_bind(close) schedule(static) reduction(+:ntargets)
        do xind=1,ldim_shrink(1),PICKER_OFFSET_HERE
            do yind=1,ldim_shrink(2),PICKER_OFFSET_HERE
                if(rmat_phasecorr(xind,yind,1) > 0.5) then
                    ntargets = ntargets + 1
                    mask(xind,yind) = .true.
                endif
            enddo
        enddo
        !$omp end parallel do
        allocate( target_corrs(ntargets),target_positions(ntargets,2))
        ntargets = 0
        do xind=1,ldim_shrink(1),PICKER_OFFSET_HERE
            do yind=1,ldim_shrink(2),PICKER_OFFSET_HERE
                if(mask(xind,yind)) then
                    ntargets = ntargets + 1
                    target_positions(ntargets,:) = [xind,yind]
                endif
            enddo
        enddo
        write(logfhandle,'(a)') '>>> REAL SPACE CORRELATIONS'
        allocate(ptcls(max(1,nthr_glob)))
        do i = 1, size(ptcls)
            call ptcls(i)%new(ldim_refs,smpd_shrunken)
        enddo
        do iref = 1,NREFS
            call refs(iref)%bp(hp,lp)
        enddo
        msk_shrunken = msk*real(orig_box)/real(ldim_refs(1))
        call circ_mask%disc(ldim_refs, smpd_shrunken, msk_shrunken, l_mask)
        call circ_mask%kill
        if(DOWRITEIMGS) then
            sdevimg = mic_shrunken
            call sdevimg%zero_and_unflag_ft
        endif
        tmp_mic = mic_saved
        call tmp_mic%bp(hp,lp)
        !$omp parallel do default(shared) private(i,outside,xind,yind,ithr) proc_bind(close) schedule(static)
        do i = 1, ntargets
            ithr = omp_get_thread_num()+1
            xind = target_positions(i,1)
            yind = target_positions(i,2)
            call tmp_mic%window_slim([xind,yind]-ldim_refs(1)/2-1, ldim_refs(1), ptcls(ithr), outside)
            target_corrs(i) = ptcls(ithr)%real_corr(refs(ref_inds(xind,yind)),l_mask)
            if(DOWRITEIMGS) call sdevimg%set([xind,yind,1],target_corrs(i))
        enddo
        !$omp end parallel do
        if(DOWRITEIMGS) then
            call sdevimg%write(PATH_HERE//basename(trim(micname))//'_corrs.mrc')
            call sdevimg%kill
        endif
        write(logfhandle,'(a)') '>>> PEAKS SELECTION'
        call sortmeans(target_corrs, MAXKMIT, means, labels)
        nmax = count(labels == 2)
        if(DEBUG_HERE) then
            write(logfhandle,'(a,1x,I7)') 'peak positions initially identified:         ', ntargets
            write(logfhandle,'(a,1x,I7)') 'peak positions after sortmeans:              ',nmax
         endif
        ! get peak positions
        allocate( peak_positions(nmax,2),  stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk( 'In: simple_picker :: gen_corr_peaks, 2',alloc_stat)
        peak_positions = 0
        nmax = 0
        do i=1,ntargets
            if( labels(i) == 2 )then
                nmax = nmax + 1
                peak_positions(nmax,:) = target_positions(i,:)
            endif
        end do
        ! cleanup
        call mask_img%kill
        call tmp_mic%kill
        do i = 1, size(ptcls)
            call ptcls(i)%kill
        enddo
        deallocate(ptcls)
    contains
        ! Reference generation and Phase Correlation calculation
        ! FORMULA: phasecorr = ifft2(fft2(field).*conj(fft2(reference)));
        subroutine gen_phase_correlation(field,mask)
            type(image),         intent(inout) :: field
            type(image),optional,intent(inout) :: mask
            type(image)   :: phasecorr, aux, ref_ext
            real, pointer :: mask_rmat(:,:,:)
            real    :: v
            integer :: iref, border
            border =  max(ldim_refs(1),ldim_refs(2))
            call phasecorr%new(ldim_shrink, smpd_shrunken)
            call aux%new(ldim_shrink, smpd_shrunken)
            call aux%set_ft(.true.)
            call ref_ext%new(ldim_shrink, smpd_shrunken)
            call field%fft
            do iref = 1, nrefs
                call refs(iref)%pad(ref_ext, 0.) ! zero padding
                call ref_ext%fft
                call field%phase_corr(ref_ext,aux,lp,border=max(ldim_refs(1)/2,ldim_refs(2)/2)) !phase correlation
                if(iref > 1) then
                    !$omp parallel do collapse(2) schedule(static) default(shared) private(i,j,v) proc_bind(close)
                    do i = 1, ldim_shrink(1)
                        do j = 1, ldim_shrink(2)
                            v = aux%get([i,j,1])
                            if( v > phasecorr%get([i,j,1]))then
                                call phasecorr%set([i,j,1], v)
                                ref_inds(i,j) = iref
                            endif
                         enddo
                    enddo
                    !$omp end parallel do
                else
                    phasecorr = aux
                    ref_inds  = 1
                endif
                call aux%fft
            enddo
            if(.not. template_based) call phasecorr%neg()
            call field%copy(phasecorr)
            if(DOWRITEIMGS) call field%write(PATH_HERE//basename(trim(micname))//'MaxValPhaseCorr.mrc')
            if(present(mask)) then
                call mask%new(ldim_shrink, smpd_shrunken)
                call mask%get_rmat_ptr(mask_rmat)
                mask_rmat(border+1:ldim_shrink(1)-border,border+1:ldim_shrink(2)-border,1)=1.
            endif
            call phasecorr%kill
            call aux%kill
            call ref_ext%kill
        end subroutine gen_phase_correlation
    end subroutine extract_peaks

    subroutine distance_filter
        integer :: ipeak, jpeak, ipos(2), jpos(2), loc(1)
        real    :: dist, corr, dist_sq, distthr_sq
        logical, allocatable :: mask(:)
        real,    allocatable :: corrs(:)
        write(logfhandle,'(a)') '>>> DISTANCE FILTERING'
        allocate( mask(nmax), corrs(nmax), selected_peak_positions(nmax), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk( 'In: simple_picker :: distance_filter',alloc_stat)
        selected_peak_positions = .true.
        distthr_sq = nint(distthr*distthr)
        !$omp parallel do schedule(static) default(shared) private(ipeak,jpeak,ipos,jpos,dist,mask,loc,corr) proc_bind(close)
        do ipeak=1,nmax
            ipos = peak_positions(ipeak,:)
            corr = -HUGE(corr)
            mask = .false.
            do jpeak=1,nmax
                jpos = peak_positions(jpeak,:)
                ! dist = euclid(real(ipos),real(jpos))
                dist_sq = sum((ipos-jpos)**2)
                if( dist_sq < distthr_sq ) then
                  mask(jpeak) = .true.
                  ! find best match in the neigh
                  if(corrmat(jpos(1),jpos(2)) > corr) then
                    corr   = corrmat(jpos(1),jpos(2))
                    loc(1) = jpeak
                  endif
                endif
              end do
            ! eliminate all but the best
            mask(loc(1)) = .false.
            where( mask )
                selected_peak_positions = .false.
            end where
        end do
        !$omp end parallel do
        nmax_sel = count(selected_peak_positions)
        if(DEBUG_HERE) write(logfhandle,'(a,1x,I7)') 'peak positions left after distance filtering: ', nmax_sel
    end subroutine distance_filter

    subroutine gather_stats
        integer           :: ipeak, cnt, istat
        logical           :: outside
        real              :: ave, maxv, minv
        real, allocatable :: spec(:)
        type(image)       :: ptcl_target
        write(logfhandle,'(a)') '>>> GATHERING REMAINING STATS'
        call ptcl_target%new(ldim_refs, smpd_shrunken)
        cnt = 0
        allocate( peak_stats(nmax_sel,NSTAT) )
        do ipeak=1,nmax
            if( selected_peak_positions(ipeak) )then
                cnt = cnt + 1
                call mic_saved%window_slim(peak_positions(ipeak,:)-ldim_refs(1)/2,&
                    &ldim_refs(1), ptcl_target, outside)
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
    end subroutine gather_stats


    ! In the old school picker there are 4 stats. In picker there are 3.
    ! The disregarded stat is the correlation (CC2REF). The effect it has, if put
    ! back, is that the picking becomes more permissive.
    subroutine one_cluster_clustering
        real, allocatable :: dmat(:,:)
        integer           :: i_median, i, j, cnt, ipeak
        real              :: ddev
        allocate(dmat(nmax_sel,nmax_sel), source=0., stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk('picker::one_cluster_clustering dmat ',alloc_stat)
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
        write(logfhandle,'(a,1x,I5)') 'peak positions after one cluster clustering: ', nmax_sel
    end subroutine one_cluster_clustering

    subroutine write_boxfile
        integer :: funit, ipeak,iostat
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(boxname)),iostat=iostat)
        call fileiochk('phasecorr_picker; write_boxfile ', iostat)
        do ipeak=1,nmax
            if( selected_peak_positions(ipeak) )then
                write(funit,'(I7,I7,I7,I7,I7)') peak_positions(ipeak,1),&
                peak_positions(ipeak,2), orig_box, orig_box, -3
            endif
        end do
        call fclose(funit,errmsg='picker; write_boxfile end')
    end subroutine write_boxfile

    subroutine kill_phasecorr_picker
        integer :: iref
        if( allocated(micname) )then
            call micrograph%kill
            call mic_shrunken%kill
            call mic_saved%kill
            deallocate(selected_peak_positions,corrmat, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('phasecorr_picker kill, 1',alloc_stat)
            deallocate(peak_positions,micname,peak_stats, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('phasecorr_picker kill, 2',alloc_stat)
            if(allocated(refsname)) deallocate(refsname)
            do iref=1,nrefs
                call refs(iref)%kill
            end do
            deallocate(refs, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('phasecorr_picker; kill 3',alloc_stat)
        endif
    end subroutine kill_phasecorr_picker
end module simple_phasecorr_picker

! IN distance filtering previous implementation of the loop
! do ipeak=1,nmax
!     ipos = peak_positions(ipeak,:)
!     mask = .false.
!     !$omp parallel do schedule(static) default(shared) private(jpeak,jpos,dist) proc_bind(close)
!     do jpeak=1,nmax
!         jpos = peak_positions(jpeak,:)
!         dist_sq = sum((ipos-jpos)**2)
!         if( dist_sq < distthr_sq ) mask(jpeak) = .true.
!         corrs(jpeak) = corrmat(jpos(1),jpos(2))
!     end do
!     !$omp end parallel do
!     ! find best match in the neigh
!     loc = maxloc(corrs, mask=mask)
!     ! eliminate all but the best
!     mask(loc(1)) = .false.
!     where( mask )
!         selected_peak_positions = .false.
!     end where
! end do
