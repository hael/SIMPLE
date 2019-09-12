module simple_segpicker
include 'simple_lib.f08'
use simple_image, only : image
use simple_parameters, only: params_glob
implicit none

 public :: segpicker

#include "simple_local_flags.inc"

! module global constants
real,    parameter :: SHRINK  = 4.  !HEREEEE
logical, parameter :: DEBUG_HERE = .false.
integer, parameter :: N_ROT = 18

type :: segpicker
    private
    type(image) :: img
    type(image) :: img_cc
    type(image) :: reference
    type(image) :: phasecorr
    real, allocatable :: particles_coord(:,:)
    real, allocatable :: stdev_gray_level(:)
    real    :: min_rad               = 0.
    real    :: max_rad               = 0.
    real    :: smpd                  = 0.
    real    :: smpd_shrunken         = 0.
    real    :: hp_box                = 0.
    real    :: lambda                = 0.  ! for tv denoising
    integer :: winsz                 = 0   ! for median filtering
    integer :: ldim(3)               = 0
    integer :: ldim_shrunken(3)      = 0
    integer :: n_particles           = 0
    integer :: orig_box              = 0
    character(len=STDLEN)         :: color      = ''   !color in which to draw on the mic to identify picked particles
    character(len=STDLEN)         :: pickername = ''   !fname
    character(len=STDLEN)         :: fbody      = ''   !fbody
    character(len=LONGSTRLEN)     :: boxname    = ''   !name of the .box file
    character(len=STDLEN)         :: detector   = ''

contains
    ! constructor
    procedure          :: new
    ! picking functions
    procedure          :: identify_particle_positions
    procedure          :: write_boxfile
    procedure          :: elimin_aggregation
    procedure, private :: relative_intensity_filtering
    procedure, private :: center_cc
    procedure, private :: center_mass_cc
    ! preprocess mic prior picking
    procedure          :: preprocess_mic
    ! setters/getters
    procedure :: get_n_ccs !TO REMOVEEE
    ! output
    procedure          :: print_info
    procedure          :: output_identified_particle_positions
    ! kill
    procedure          :: kill
end type segpicker

contains

    subroutine new(self, fname, min_rad, max_rad, smpd, color)
        class(segpicker),       intent(inout) :: self
        character(len=*),           intent(in)    :: fname
        real,                       intent(in)    :: min_rad
        real,                       intent(in)    :: max_rad
        real,                       intent(in)    :: smpd
        character(len=*), optional, intent(in)    :: color
        integer :: nptcls
        self%pickername = fname
        call find_ldim_nptcls(self%pickername, self%ldim, nptcls, self%smpd)
        self%smpd = smpd
        self%ldim_shrunken(1) = round2even(real(self%ldim(1))/SHRINK)
        self%ldim_shrunken(2) = round2even(real(self%ldim(2))/SHRINK)
        self%ldim_shrunken(3) = 1
        self%smpd_shrunken = self%smpd*SHRINK
        self%min_rad = min_rad
        self%max_rad = max_rad
        self%fbody = get_fbody(trim(fname), trim(fname2ext(fname)))
        ! self%boxname = basename( fname_new_ext(self%fbody,'box') )
        self%boxname =  PATH_HERE//basename(trim(self%fbody))//'.box'
        call self%img%new      (self%ldim_shrunken, self%smpd_shrunken)
        call self%img_cc%new   (self%ldim_shrunken, self%smpd_shrunken)
        call self%reference%new(self%ldim_shrunken, self%smpd_shrunken)
        call self%phasecorr%new(self%ldim_shrunken, self%smpd_shrunken)
        self%n_particles = 0
        self%lambda      = 3. !HEREEE
        self%detector    = 'bin'
        self%color       = 'white' !default
        if(present(color)) self%color = color
    end subroutine new

    subroutine preprocess_mic(self, detector)
        use simple_tvfilter, only : tvfilter
        use simple_micops
        use simple_segmentation, only: sobel, automatic_thresh_sobel
        class(segpicker), intent(inout) :: self
        character(len= *),    intent(in)    :: detector
        real, allocatable :: rmat(:,:,:)
        real, allocatable :: x(:), x_out(:)
        type(tvfilter)    :: tvf
        type(image) :: micrograph_shrunken
        integer     :: box_shrunken
        real        :: ave, sdev, maxv, minv
        real        :: thresh(1)
        real        :: sigma_x, sigma_y
        type(image) :: aux !HEREE can I get rid of it?
        type(image) :: mask_img
        ! 0) Reading and saving original micrograph
        call read_micrograph(self%pickername, smpd=self%smpd)
        ! 1) Shrink and high pass filtering
        call shrink_micrograph(SHRINK, self%ldim_shrunken, self%smpd_shrunken)
        self%hp_box =  4.*self%max_rad+2.*self%max_rad
        call set_box(int(SHRINK*(self%hp_box)), box_shrunken, micrograph_shrunken)
        call self%img%copy(micrograph_shrunken)
        ! To take care of shrinking
        self%min_rad = self%min_rad/SHRINK ! I am thingking I shouldn't multiply by the smpd cuz I am working in pxls
        self%max_rad = self%max_rad/SHRINK
        ! 2) Low pass filtering
        ! call micrograph_shrunken%bp(0.,40.) !HEREEE
        call micrograph_shrunken%bp(0., params_glob%lp)
        call micrograph_shrunken%ifft()
        ! 2.1) TV denoising !HEREEEE
        ! call tvf%new()
        ! call tvf%apply_filter(micrograph_shrunken, self%lambda)
        ! call tvf%kill
        ! if(DEBUG_HERE) call micrograph_shrunken%write_jpg(PATH_HERE//basename(trim(self%fbody))//'_tvfiltered.jpg')
        ! call micrograph_shrunken%write(PATH_HERE//basename(trim(self%fbody))//'_tvfiltered.mrc')
        ! 2.2) Negative image, to obtain a binarization with white particles
        call micrograph_shrunken%neg() !TO REMOVE IN CASE OF NEGATIVE STAINING
        ! I PADDED HEREEE
        ! 2.3) New approach with phase correlation calculation
        call gen_phase_correlation(micrograph_shrunken,maxv,minv,mask_img)
        call micrograph_shrunken%stats( ave=ave, sdev=sdev, maxv=maxv, minv=minv,mskimg=mask_img) !HEREE, it was before phase correlation calculation
        ! 3) Binarization
        if(detector .eq. 'sobel') then
            self%detector = 'sobel'
            thresh(1) = ave+.5*sdev !sobel needs lower thresh not to pick just edges
            call sobel(micrograph_shrunken,thresh)
        else if (detector .eq. 'bin') then
            ! default self%detector is bin
            call micrograph_shrunken%bin(ave+.8*sdev)
        else if (detector .eq. 'otsu') then
            self%detector = 'otsu'
            rmat = micrograph_shrunken%get_rmat()
            x = pack(rmat, .true.)
            call otsu(x,x_out)
            rmat = reshape(x_out, [self%ldim_shrunken(1),self%ldim_shrunken(2),1])
            call micrograph_shrunken%set_rmat(rmat)
            call micrograph_shrunken%erosion() !morphological erosion
            deallocate(x,x_out,rmat)
        else
            THROW_HARD('Invalid detector; preprocess_mic')
        endif
        if( DEBUG_HERE ) call micrograph_shrunken%write_jpg(PATH_HERE//basename(trim(self%fbody))//'_Bin.jpg')
        call micrograph_shrunken%write(PATH_HERE//basename(trim(self%fbody))//'_Bin.mrc')
        self%winsz = int(self%min_rad+self%max_rad)/4 !half of the avg between the dimensions of the particles
        call micrograph_shrunken%real_space_filter(self%winsz,'median') !median filtering allows easy calculation of cc
        if( DEBUG_HERE ) call micrograph_shrunken%write_jpg(PATH_HERE//basename(trim(self%fbody))//'_BinMedian.jpg')
        call micrograph_shrunken%write(PATH_HERE//basename(trim(self%fbody))//'_BinMedian.mrc')
        ! 5) Connected components (cc) identification
        call micrograph_shrunken%find_connected_comps(self%img_cc)
        ! 6) cc filtering !HEREEE
        print *, 'before polishing the ccs: ', self%get_n_ccs()
        call self%img_cc%polish_cc(self%min_rad,self%max_rad)
        print *, 'after polishing the ccs: ', self%get_n_ccs()
        if( DEBUG_HERE ) call self%img_cc%write_jpg(PATH_HERE//basename(trim(self%fbody))//'_ConnectedComponentsElimin.jpg')
        call self%img_cc%write(PATH_HERE//basename(trim(self%fbody))//'_ConnectedComponentsElimin.mrc')
        call micrograph_shrunken%kill

    contains

        ! Reference generation and Phase Correlation calculation
        ! FORMULA: phasecorr = ifft2(fft2(field).*conj(fft2(reference)));
        subroutine gen_phase_correlation(field,maxv, minv,mask_img)
            use simple_procimgfile, only :  clip_imgfile
            type(image), intent(inout) :: field
            real,        intent(in)    :: maxv, minv
            type(image), intent(out)   :: mask_img
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            real, allocatable          :: rmat_out(:,:,:) !to use rtsq_serial
            real    :: sigma_x, sigma_y
            real    :: rot_step ! reference angular rotation step
            integer :: n_ang
            logical :: rotate_ref
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            type(image) :: pickref
            type(image) :: pickref_ext
            real    :: ref_smpd
            real    :: maskrad
            real    :: shrink_factor
            integer :: NREFS, n_ref
            integer :: ref_dim(3), ldim_clip(3)
            integer :: border
            integer :: i, j, nptcls
            border = 3*nint(self%max_rad)
            maskrad= 60.
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! sigma_x = (self%min_rad)/2. ! half of the radius
            ! sigma_y = (self%max_rad)/2.
            ! rotate_ref = .false.
            ! if(sigma_y/sigma_x > 1.5) rotate_ref = .true.
            ! print *, 'sigma_x = ', sigma_x, 'sigma_y = ', sigma_y
            ! rot_step = 180./real(N_ROT) !it's gaussian, so rotate just 180 degrees
            ! call self%reference%gauimg2D(sigma_x,sigma_y)
            ! call self%reference%scale_pixels([minv,maxv]) !should I??
            ! call self%reference%write(PATH_HERE//basename(trim(self%fbody))//'_GaussianReference.mrc')
            call field%fft()
            ! allocate(rmat_out(self%ldim_shrunken(1),self%ldim_shrunken(2),1), source = 0.)
            ! call self%phasecorr%zero_and_unflag_ft()
            ! print *, 'rotate_ref = ', rotate_ref
            !!!!!!!!!!!!!!TENMPLATEEEEE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !call aux2%new(self%ldim_shrunken, self%smpd_shrunken)
            call pickref_ext%new(self%ldim_shrunken, self%smpd_shrunken)
            call aux%new        (self%ldim_shrunken, self%smpd_shrunken)
            call aux%set_ft(.true.)
            call find_ldim_nptcls('../pickrefs.mrc', ref_dim, nptcls, ref_smpd)
            NREFS = nptcls
            shrink_factor = self%smpd_shrunken/ref_smpd
            ldim_clip(1) = round2even(real(ref_dim(1))/shrink_factor)
            ldim_clip(2) = round2even(real(ref_dim(2))/shrink_factor)
            ldim_clip(3) = 1
            call clip_imgfile('../pickrefs.mrc', 'pickrefs_rightsmpd.mrc', ldim_clip, ref_smpd )
            call pickref%new(ldim_clip,self%smpd_shrunken)
            do n_ref = 1, NREFS
                call pickref%read('pickrefs_rightsmpd.mrc', n_ref)
                call pickref%norm() !normalise, DOESN'T HAVE ANY EFFECT, CAN BE REMOVED
                call pickref%pad(pickref_ext, 0.) ! zero padding
                call pickref_ext%mask(mskrad=maskrad, which='soft', backgr=0.)
                call pickref_ext%write('pickref_extended_masked.mrc', n_ref)
                call pickref_ext%fft
                call field%phase_corr(pickref_ext,aux,params_glob%lp,border=border) !correlation
                if(n_ref > 1) then
                    call max_image(self%phasecorr,self%phasecorr,aux) !save in phasecorr the maximum value between previous phasecorr and new phasecorr
                else
                    self%phasecorr = aux
                endif
                call self%phasecorr%write('MaxValPhaseCorr.mrc',n_ref)
                call aux%fft()
            enddo
            call self%phasecorr%write(PATH_HERE//basename(trim(self%fbody))//'MaxValPhaseCorr.mrc')
            call field%copy(self%phasecorr)
            call field%neg() ! to fix wrt the previous negative
            call mask_img%new(self%ldim_shrunken, self%smpd_shrunken)
            do i = border+1, self%ldim_shrunken(1) - border
                do j = border+1, self%ldim_shrunken(2) - border
                    call mask_img%set([i,j,1], 1.)
                enddo
            enddo
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! if(rotate_ref) then
            !     do n_ang = 1, N_ROT
            !         if(n_ang > 1) then !do not rotate the first time
            !             call self%reference%ifft() ! in order to rotate it has to be real (don't need to check if it's FT, it' done internally)
            !             call self%reference%rtsq_serial( rot_step, 0., 0., rmat_out )
            !             call self%reference%set_rmat(rmat_out)
            !         endif
            !         ! I think I can improve this bit, write just one command
            !         call self%reference%fft()
            !         call self%phasecorr%fft()    ! before sum they need to be in the same state
            !         aux = self%reference%conjg() ! necessary
            !         call max_image(self%phasecorr,self%phasecorr,field*aux) !save in phasecorr the maximum value between phasecorr and field*aux
            !     enddo
            ! else
            !     call self%reference%fft()
            !     call self%phasecorr%fft()    ! before sum they need to be in the same state
            !     aux = self%reference%conjg() ! necessary
            !     call self%phasecorr%copy(field*aux)
            ! endif
            ! call self%phasecorr%ifft()
            ! call self%phasecorr%write(PATH_HERE//basename(trim(self%fbody))//'MaxValPhaseCorr.mrc')
            ! call field%copy(self%phasecorr)
        end subroutine gen_phase_correlation

        ! This subroutine creates an image in which each pixel value
        ! is the max value between the gray level in img1 and img2
        subroutine max_image(img,img1,img2)
            type(image), intent(inout) :: img !output image
            type(image), intent(in)    :: img1, img2
            ! real, pointer     :: rmat(:,:,:)
            real, allocatable :: rmat(:,:,:)
            real, allocatable :: rmat1(:,:,:), rmat2(:,:,:) !matrices for img1 and img2
            integer :: i, j
            ! call img%get_rmat_ptr(rmat)
            rmat  = img%get_rmat()
            rmat1 = img1%get_rmat()
            rmat2 = img2%get_rmat()
            do i = 1, self%ldim_shrunken(1)
                do j = 1, self%ldim_shrunken(2)
                    if(rmat1(i,j,1) > 0. .and. rmat2(i,j,1) > 0.) then
                        rmat(i,j,1) = max(rmat1(i,j,1),rmat2(i,j,1))
                    elseif(rmat1(i,j,1) < 0. .and. rmat2(i,j,1) < 0.) then
                        rmat(i,j,1) = min(rmat1(i,j,1),rmat2(i,j,1))
                    else
                        rmat(i,j,1) = (rmat1(i,j,1)+rmat2(i,j,1))/2.
                    endif
                 enddo
            enddo
            call img%set_rmat(rmat)
        end subroutine max_image

        !>  \brief clip_imgfile is for clipping
        !! \param fname2clip output filename
        !! \param fname input filename
        !! \param ldim_clip clipped logical dimension
        !! \param smpd sampling distance TO FIX COMMENTSSS
        ! I AM NOT ACTUALLY USING THIS SUBROUTINE BUT I USE IT
        ! TO GENERATE THE REFERENCES WITH THE CORRECT SMPD
        ! SO I SAVE IT HERE SO I KNOW HOW IT WORKS AND IT
        ! DOESN'T GET LOST.
        subroutine clip_vol( vol_in, vol_out, ldim_clip, smpd )
            type(image), intent(inout) :: vol_in
            type(image), intent(out)   :: vol_out
            integer,     intent(in)    :: ldim_clip(3)
            real,        intent(in)    :: smpd
            type(image) :: img, img_clip
            integer     :: n, i, ldim(3)
            ldim = vol_in%get_ldim()
            ! do the work
            if( ldim_clip(1) <= ldim(1) .and. ldim_clip(2) <= ldim(2)&
                 .and. ldim_clip(3) <= ldim(3) )then
                call img_clip%new(ldim_clip,smpd)
                write(logfhandle,'(a)') '>>> CLIPPING VOL IN FOURIER SPACE'
                ! soft mask the vol. You don't know where it comes from.
                ! Need to get rid of border effects
                call img%mask(144.,'soft') ! in general, use the maximum radius + something --> 1.1*max_rad
                call img%fft() !clip in Fourier Space to change SMPD
                call img%clip(img_clip) ! FT state preserved
                call img_clip%ifft() ! Come back in real space
            else
                ! need to downscale, how??
            end if
        end subroutine clip_vol
    end subroutine preprocess_mic

    subroutine print_info(self)
        class(segpicker), intent(inout) :: self
        open(unit = 17, file = "PickerInfo.txt")
        write(unit = 17, fmt = '(a)') '>>>>>>>>>>>>>>>>>>>>PARTICLE PICKING>>>>>>>>>>>>>>>>>>'
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a,f0.0)")             'Mic Shrunken, fact ', SHRINK
        write(unit = 17, fmt = "(a,i4,tr1,i4,tr1,i4)") 'Dim before  shrink ', self%ldim
        write(unit = 17, fmt = "(a,i4,tr1,i4,tr1,i4)") 'Dim after   shrink ', self%ldim_shrunken
        write(unit = 17, fmt = "(a,f4.2)")             'Smpd before shrink ', self%smpd
        write(unit = 17, fmt = "(a,f4.2)")             'Smpd after  shrink ', self%smpd_shrunken
        write(unit = 17, fmt = "(a,a)")                'Hp box              ', trim(int2str(int(self%hp_box)))
        write(unit = 17, fmt = "(a,i4,tr1,i4)")        'Ccs size filtering ', int(self%min_rad*self%max_rad), int(2*3.14*(self%max_rad)**2)
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a)")  'SELECTED PARAMETERS '
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a,f4.2)")  'Lp filter parameter     ', params_glob%lp
        write(unit = 17, fmt = "(a,i4)")    'Winsz for median filter ', self%winsz
        write(unit = 17, fmt = "(a,f4.2)")  'TV filter parameter     ', self%lambda
        write(unit = 17, fmt = "(a,a,a,a)") 'particle dimensions     ', trim(int2str(int(self%min_rad))),' ', trim(int2str(int(self%max_rad)))
        write(unit = 17, fmt = "(a,a)")     'detector                ', self%detector
        close(17, status = "keep")
    end subroutine print_info

    ! This routine is aimed to eliminate aggregations of particles.
    ! It takes in input the list of the coordinates of the identified
    ! particles and the approximate radius of the particle.
    ! The output is a new list of coordinates where aggregate picked particles
    ! have been deleted
    ! If the dist beetween the 2 particles is < 1.5*r -> delete them. (not too stringent)
    subroutine elimin_aggregation( self, saved_coord, refined_coords )
        class(segpicker),     intent(inout)  :: self
        real,                 intent(in)     :: saved_coord(:,:)     !Coordinates of picked particles
        real, allocatable,    intent(out)    :: refined_coords(:,:)  !New coordinates of not aggregated particles
        logical, allocatable :: msk(:)
        integer              :: i, j, cnt
        !!!!!!!!!!!!!!!!
        logical :: outside
        type(image) :: imgwin_particle
        call imgwin_particle%new([self%orig_box,self%orig_box,1],self%smpd)
        !!!!!!!!!!!!!!!
        allocate(msk(self%n_particles), source = .true.) !initialise to 'keep all'
        do i = 1, self%n_particles             !fix one coord
            do j = i+1, self%n_particles       !fix another coord to compare
                if(msk(i) .and. msk(j)) then !not compare twice ,and if the particles haven t been deleted yet
                    if(  sqrt(real((saved_coord(i,1)-saved_coord(j,1))**2 &
                        &        + (saved_coord(i,2)-saved_coord(j,2))**2)) <= 1.5*self%min_rad) then!&
                        !& .and. self%min_rad < &
                        !& real((saved_coord(i,1)-saved_coord(j,1))**2 &        ! TO CHECK WHAT TO DO
                        !&    + (saved_coord(i,2)-saved_coord(j,2))**2)) then
                        msk(i) = .false.
                        msk(j) = .false.
                    endif
                endif
            enddo
        enddo
        !print *, 'before relative intensiry filtering nb of particles is ', count(msk)
        call self%relative_intensity_filtering(msk)
        ! print *, 'after  relative intensiry filtering nb of particles is ', count(msk)
        allocate(refined_coords(count(msk),2), source = 0.)
        cnt = 0
        do i = 1, self%n_particles
            if(msk(i)) then
                cnt = cnt + 1
                refined_coords(cnt, :) = saved_coord(i, :)
            endif
        enddo
        deallocate(msk)
    end subroutine elimin_aggregation


    subroutine relative_intensity_filtering(self, selected_particle_positions)
        class(segpicker),     intent(inout) :: self
        logical,              intent(inout) :: selected_particle_positions(:)
        integer, allocatable :: imat(:,:,:)
        logical, allocatable :: mask(:,:,:)
        real,    allocatable :: avg_gray_level(:), stdev_gray_level(:)
        integer :: n_cc
        real    :: avg_level
        integer :: m(1), mm(1)
        imat = nint(self%img_cc%get_rmat())
        allocate(mask(self%ldim_shrunken(1),self%ldim_shrunken(2),1), source = .false.)
        allocate(avg_gray_level(maxval(imat)),self%stdev_gray_level(maxval(imat)),source = 0.) ! SUBSTITUTE MAXVAL(IMAT) WITH THE NUMBER OF PARTICLES
        do n_cc = 1, maxval(imat) ! for each connected component
            where(imat == n_cc) mask = .true.
            call calc_avgst_intensity_particle(mask,avg_gray_level(n_cc),self%stdev_gray_level(n_cc))
            mask = .false. !restore
        enddo
        ! avg_level   = 0. !initialise
        ! stdev_level = 0.
        ! avg_level   = sum(avg_gray_level)/real(size(avg_gray_level))
        ! do n_cc = 1, size(avg_gray_level)
        !     stdev_level = stdev_level + (avg_gray_level(n_cc)-avg_level)**2
        ! enddo
        ! stdev_level = sqrt(stdev_level/real(size(avg_gray_level)-1))
        ! ! FIlter assuming gaussian distribution in the gray levels of the particles
        ! do n_cc = 1, size(selected_particle_positions)
        !     if(avg_gray_level(n_cc) > avg_level + 2.*stdev_level .or. avg_gray_level(n_cc) < avg_level - 2.*stdev_level) selected_particle_positions(n_cc) = .false.
        ! enddo
        ! print *, 'thesh up = ', avg_level + 2.*stdev_level , 'thresh down ', avg_level - 2.*stdev_level
    contains

        subroutine calc_avgst_intensity_particle(mask,avg_gray_level,stdev_gray_level)
            logical, intent(in)  :: mask(:,:,:)
            real,    intent(out) :: avg_gray_level
            real,    intent(out) :: stdev_gray_level
            real,    allocatable :: rmat(:,:,:)
            integer :: i,j
            rmat = self%img%get_rmat()
            avg_gray_level   = 0.
            stdev_gray_level = 0.
            do i = 1, self%ldim_shrunken(1)
                do j = 1, self%ldim_shrunken(2)
                    if(mask(i,j,1)) avg_gray_level = avg_gray_level + rmat(i,j,1)
                enddo
            enddo
            avg_gray_level = avg_gray_level/count(mask)
            do i = 1, self%ldim_shrunken(1)
                do j = 1, self%ldim_shrunken(2)
                    if(mask(i,j,1)) stdev_gray_level = stdev_gray_level + (avg_gray_level-rmat(i,j,1))**2
                enddo
            enddo
            stdev_gray_level = sqrt(stdev_gray_level/(real(count(mask)-1)))
        end subroutine calc_avgst_intensity_particle
    end subroutine relative_intensity_filtering


    ! ! Adapted from simple_picker.f90
    ! subroutine distance_filter(self)
    !     class(segpicker),  intent(inout) :: self
    !     integer :: i, j, ipos(2), jpos(2), loc(1)
    !     real    :: dist
    !     logical, allocatable :: mask(:)
    !     real,    allocatable :: corrs(:)
    !     write(logfhandle,'(a)') '>>> DISTANCE FILTERING'
    !     allocate( mask(npeaks), corrs(npeaks), selected_peak_positions(npeaks), stat=alloc_stat)
    !     if(alloc_stat.ne.0)call allocchk( 'In: simple_picker :: distance_filter',alloc_stat)
    !     selected_peak_positions = .true.
    !     do i=1, size(saved_coord, dim = 1)
    !         ipos = saved_coord(i,:)
    !         mask = .false.
    !         !$omp parallel do schedule(static) default(shared) private(j,jpos,dist) proc_bind(close)
    !         do j=1,npeaks
    !             jpos = saved_coord(j,:)
    !             dist = euclid(real(ipos),real(jpos))
    !             if( dist < 2.*self%min_rad ) mask(j) = .true.
    !             corrs(j) = corrmat(jpos(1),jpos(2))
    !         end do
    !         !$omp end parallel do
    !         ! find best match in the neigh
    !         loc = maxloc(corrs, mask=mask)
    !         ! eliminate all but the best
    !         mask(loc(1)) = .false.
    !         where( mask )
    !             selected_peak_positions = .false.
    !         end where
    !     end do
    !     npeaks_sel = count(selected_peak_positions)
    !     write(logfhandle,'(a,1x,I5)') 'peak positions left after distance filtering: ', npeaks_sel
    ! end subroutine distance_filter

    ! This subroutine takes in input an image, its connected components image,
    ! and extract particles. It doesn't use mass_center.
    ! notation:: cc = connected component.
    subroutine identify_particle_positions(self)
      class(segpicker), intent(inout) :: self
      integer              :: n_cc
      integer              :: cnt
      real                 :: pos(3)            !center of each cc
      real, allocatable    :: xyz_saved(:,:), xyz_norep_noagg(:,:)
      integer, allocatable :: imat(:,:,:), imat_cc(:,:,:)
      ! Initialisations
      self%orig_box = int(4.*(self%max_rad)+2.*self%max_rad) !needs to be bigger than the particle
      imat_cc = int(self%img_cc%get_rmat())
      allocate(xyz_saved(maxval(imat_cc),2), source = 0.) ! size of the # of cc (likely particles)
      allocate(imat(1:self%ldim_shrunken(1),1:self%ldim_shrunken(2),1:self%ldim_shrunken(3)), source = 0)
      ! Particle identification, extraction and centering
      where(imat_cc > 0.5) imat = 1
      do n_cc = 1, maxval(imat_cc)
          ! pos(:) = self%center_cc(n_cc)
          pos(:) = self%center_mass_cc(n_cc)
          xyz_saved(n_cc,:) = pos(:2)
      enddo
      deallocate(imat_cc)
      self%n_particles = size(xyz_saved, dim = 1) ! first estimation
      print *, 'before elimin aggregations: ', size(xyz_saved, dim=1)
      call self%elimin_aggregation(xyz_saved, xyz_norep_noagg)      ! x2 not to pick too close particles
       allocate(self%particles_coord(count(xyz_norep_noagg(:,1)> TINY),2), source = 0.) !TO IMPROVE
       ! allocate(self%particles_coord(size(xyz_saved,dim=1),2), source = xyz_saved) !TO IMPROVE
      cnt = 0 !initialise
      do n_cc = 1,  size(xyz_norep_noagg, dim=1)
          if(abs(xyz_norep_noagg(n_cc,1)) > TINY .and. abs(xyz_norep_noagg(n_cc,2))> TINY) then
              cnt = cnt + 1
              self%particles_coord(cnt,:) = xyz_norep_noagg(n_cc,:)
          endif
      enddo
      self%n_particles = size(self%particles_coord, dim = 1) !update after elim aggregations
      print *, 'after elimin aggregations: ', self%n_particles
      if(allocated(xyz_saved))       deallocate(xyz_saved)
      if(allocated(xyz_norep_noagg)) deallocate(xyz_norep_noagg)
    end subroutine identify_particle_positions

      subroutine output_identified_particle_positions(self)
          class(segpicker), intent(inout) :: self
          type(image) :: imgwin_particle
          integer :: n_cc
          integer :: cnt
          logical :: outside(self%n_particles)
          real    :: dev
          outside(:) = .false.
          cnt = 0
          call imgwin_particle%new([self%orig_box,self%orig_box,1],self%smpd)
          open(unit = 23, file = basename(trim(self%fbody))//"Stdev.txt")
          write(unit = 23, fmt = "(a,f4.2)") basename(trim(self%fbody))
          call self%img%rmsd(dev)
          write(unit = 23, fmt = "(a,f4.2)") 'Stdev ', dev
          do n_cc = 1, self%n_particles
              call self%img%window_slim(nint(self%particles_coord(n_cc,:)-self%orig_box/2), self%orig_box, imgwin_particle, outside(n_cc))
              if( .not. outside(n_cc)) then
                  cnt = cnt + 1
                  if( DEBUG_HERE )call imgwin_particle%write_jpg(PATH_HERE//basename(trim(self%fbody))//'_centered_particles.jpg', cnt)
                  call imgwin_particle%write(PATH_HERE//basename(trim(self%fbody))//'_centered_particles.mrc', cnt)
                  write(unit = 23, fmt = "(a,i4,a,f4.2)") 'Particle ', cnt, ' Stdev ',  self%stdev_gray_level(n_cc)
              endif
          enddo
          close(23)
          ! Cannot put in the same loop as before otherwise when I extract particles I can see the drawings
          do n_cc = 1, self%n_particles
              if(.not. outside(n_cc)) then
                  call self%img%draw_picked(nint(self%particles_coord(n_cc,:)),nint((self%min_rad+self%max_rad)/2.),2, self%color)
              endif
          end do
          call imgwin_particle%kill
          if( DEBUG_HERE ) call self%img%write_jpg(PATH_HERE//basename(trim(self%fbody))//'_picked_particles.jpg')
           call self%img%write(PATH_HERE//basename(trim(self%fbody))//'_picked_particles.mrc')
      end subroutine output_identified_particle_positions

    ! This function returns the index of a pixel (assuming to have a 2D)
    ! image in a connected component. The pixel identified is the one
    ! that minimizes the distance between itself and all the other pixels of
    ! the connected component. It corresponds to the geometric median.
    function center_cc(self,n_cc) result (px)
        class(segpicker),  intent(inout) :: self
        integer,               intent(in)    :: n_cc
        real        :: px(3)               !index of the central px of the cc
        integer     :: i, j, k
        integer     :: n_px                !counter
        integer     :: idx(2)
        real,    allocatable :: dist(:,:)        !to extract the window according to the px which minimize the dist
        logical, allocatable :: mask(:,:)        !to calc the min of an array along a specific dim
        logical, allocatable :: mask_dist(:)     !for sum dist calculation
        integer, allocatable :: pos(:,:)         !position of the pixels of a fixed cc
        integer, allocatable :: imat_cc(:,:,:)
        imat_cc = int(self%img_cc%get_rmat())
        where(imat_cc .ne. n_cc) imat_cc = 0
        call get_pixel_pos(imat_cc,pos)
        allocate(dist(4,   size(pos, dim = 2)), source = 0.)
        allocate(mask(4,   size(pos, dim = 2)), source = .false.)
        allocate(mask_dist(size(pos, dim = 2)), source = .true.)
        mask(4,:) = .true. !to calc the min wrt the dist
        n_px = 0
        do i = 1, self%ldim_shrunken(1)
            do j = 1, self%ldim_shrunken(2)
                do k = 1, self%ldim_shrunken(3)
                    if(imat_cc(i,j,k) > 0.5) then
                        n_px = n_px + 1
                        dist( 4, n_px) = pixels_dist([i,j,k], pos, 'sum', mask_dist)
                        dist(:3, n_px) = [real(i),real(j),real(k)]
                    endif
                enddo
            enddo
        enddo
        idx   = minloc(dist, mask)
        px(:) = dist(1:3, idx(2))
        deallocate(pos, mask, dist, mask_dist)
        if(allocated(imat_cc)) deallocate(imat_cc)
    end function center_cc

    ! This function returns the index of a pixel (assuming to have a 2D)
    ! image in a connected component. The pixel identified is the one
    ! center of mass of the cc.
    function center_mass_cc(self,n_cc) result (px)
        class(segpicker),  intent(inout) :: self
        integer,               intent(in)    :: n_cc
        real        :: px(3)               !index of the central px of the cc
        integer, allocatable :: pos(:,:)         !position of the pixels of a fixed cc
        integer, allocatable :: imat_cc(:,:,:)
        imat_cc = int(self%img_cc%get_rmat())
        where(imat_cc .ne. n_cc) imat_cc = 0
        call get_pixel_pos(imat_cc,pos)
        px(1) = sum(pos(1,:))/real(size(pos,dim = 2))
        px(2) = sum(pos(2,:))/real(size(pos,dim = 2))
        px(3) = 1.
        if(allocated(imat_cc)) deallocate(imat_cc)
    end function center_mass_cc

    function get_n_ccs(self) result(ncc)
        class(segpicker),  intent(inout) :: self
        integer :: ncc
        integer, allocatable :: imat(:,:,:)
        imat = nint(self%img_cc%get_rmat())
        ncc  = maxval(imat)
    end function get_n_ccs


    subroutine write_boxfile(self)
        class(segpicker),  intent(inout) :: self
        integer :: funit, n_cc,iostat
        call fopen(funit, status='REPLACE', action='WRITE', file=self%boxname,iostat=iostat)
        call fileiochk('picker; write_boxfile ', iostat)
        do n_cc=1,self%n_particles
            write(funit,'(I7,I7,I7,I7,I7)') int(self%particles_coord(n_cc,1)),&
            int(self%particles_coord(n_cc,2)), self%orig_box, self%orig_box, -3
        end do
        call fclose(funit,errmsg='picker; write_boxfile end')
    end subroutine write_boxfile


    subroutine kill(self)
        class(segpicker),  intent(inout) :: self
        !kill images
        call self%img%kill
        call self%img_cc%kill
        call self%reference%kill
        call self%phasecorr%kill
        if(allocated(self%particles_coord)) deallocate(self%particles_coord)
        if(allocated(self%stdev_gray_level)) deallocate(self%stdev_gray_level)
        self%min_rad          = 0.
        self%max_rad          = 0.
        self%smpd             = 0.
        self%smpd_shrunken    = 0.
        self%hp_box           = 0.
        self%winsz            = 0
        self%ldim(:)          = 0
        self%ldim_shrunken(:) = 0
        self%n_particles      = 0
        self%orig_box         = 0
        self%color            = ''
        self%pickername       = ''   !fname
        self%fbody            = ''   !fbody
        self%boxname          = ''
        self%detector         = ''
        self%lambda           = 0.
    end subroutine kill
end module simple_segpicker
