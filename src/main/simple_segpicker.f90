module simple_segpicker
include 'simple_lib.f08'
use simple_image,      only : image
use simple_binimage,   only: binimage
use simple_parameters, only: params_glob
implicit none

 public :: segpicker

#include "simple_local_flags.inc"

! module global constants
real,    parameter :: SHRINK      = 4.
logical, parameter :: DEBUG_HERE  = .false.
logical, parameter :: DOWRITEIMGS = .false.
integer, parameter :: N_ROT       = 18
integer, parameter :: MIN_NCCS    = 5    ! minimum number of connected components to identify after size-filtering

type :: segpicker
    private
    type(image)       :: mic
    type(binimage)    :: img   ! mic shrunken
    type(binimage)    :: img_cc
    real, allocatable :: particles_coord(:,:)
    real, allocatable :: avg_gray_level(:),stdev_gray_level(:)
    real    :: min_rad               = 0.
    real    :: max_rad               = 0.
    real    :: smpd                  = 0.
    real    :: smpd_shrunken         = 0.
    real    :: hp_box                = 0.
    real    :: lambda                = 0.  ! for tv denoising
    real    :: distthr               = 0.
    integer :: winsz                 = 0   ! for median filtering
    integer :: ldim(3)               = 0
    integer :: ldim_shrunken(3)      = 0
    integer :: n_particles           = 0
    integer :: orig_box              = 0
    character(len=3)              :: elongated  = ''   !user inputted particle shape
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
    ! procedure, private :: center_cc    ! geometric center
    procedure, private :: center_mass_cc ! center of mass
    ! preprocess mic prior picking
    procedure          :: preprocess_mic
    ! debugging output
    procedure          :: print_info
    procedure          :: output_identified_particle_positions
    ! kill
    procedure          :: kill
end type segpicker

contains

    subroutine new(self, fname, min_rad, max_rad, eelongated, smpd, color, distthr_in)
        class(segpicker),       intent(inout) :: self
        character(len=*),           intent(in)    :: fname
        real, optional,             intent(in)    :: distthr_in
        real,                       intent(in)    :: min_rad, max_rad
        character(len=3),           intent(in)    :: eelongated
        real,                       intent(in)    :: smpd
        character(len=*),           intent(in)    :: color
        integer :: nptcls
        self%pickername = fname
        call find_ldim_nptcls(self%pickername, self%ldim, nptcls, self%smpd)
        self%smpd = smpd
        self%ldim_shrunken(1) = round2even(real(self%ldim(1))/SHRINK)
        self%ldim_shrunken(2) = round2even(real(self%ldim(2))/SHRINK)
        self%ldim_shrunken(3) = 1
        self%smpd_shrunken = self%smpd*SHRINK
        self%min_rad = min_rad/(SHRINK*self%smpd)    !in pixel, shrunken dimensions
        self%max_rad = max_rad/(SHRINK*self%smpd)    !in pixel, shrunken dimensions
        self%distthr = 2.*self%max_rad/self%smpd     !in pixel, shrunken dimensions
        if( present(distthr_in) ) self%distthr = distthr_in/SHRINK
        self%fbody = get_fbody(trim(fname), trim(fname2ext(fname)))
        self%boxname =  basename(self%fbody)//'.box'
        call self%mic%new        (self%ldim, self%smpd)
        call self%img%new        (self%ldim_shrunken, self%smpd_shrunken)
        call self%img_cc%new_bimg(self%ldim_shrunken, self%smpd_shrunken)
        self%n_particles = 0
        self%lambda      = 3.
        self%color       = color
        ! set shape
        self%elongated = eelongated
    end subroutine new

    subroutine preprocess_mic(self, detector)
        use simple_segmentation, only : otsu_robust_fast
        use simple_tvfilter,     only : tvfilter
        use simple_micops
        class(segpicker),  intent(inout) :: self
        character(len= *), intent(in)    :: detector
        type(tvfilter) :: tvf
        type(binimage) :: micrograph_shrunken
        real           :: ave, sdev, maxv, minv
        real           :: ttthresh(3)
        real           :: sigma_x, sigma_y
        integer        :: box_shrunken
        type(image)    :: mask_img
        ! 0) Reading and saving original micrograph
        call read_micrograph(self%pickername, smpd=self%smpd, mic_out=self%mic)
        ! 1) Shrink and high pass filtering
        call shrink_micrograph(SHRINK, self%ldim_shrunken, self%smpd_shrunken)
        self%hp_box =  4.*self%max_rad+2.*self%max_rad
        call set_box(int(SHRINK*(self%hp_box)),box_shrunken,micrograph_shrunken)
        call self%img%copy(micrograph_shrunken)
        ! 2) Low pass filtering
        call micrograph_shrunken%bp(0., params_glob%lp)
        call micrograph_shrunken%ifft()
        ! 2.1) TV denoising !HEREEEE
        call tvf%new()
        call tvf%apply_filter(micrograph_shrunken, self%lambda)
        call tvf%kill
        if(DOWRITEIMGS) call micrograph_shrunken%write(PATH_HERE//basename(trim(self%fbody))//'_tvfiltered.mrc')
        ! 2.2) Negative image, to obtain a binarization with white particles
        call micrograph_shrunken%neg() !TO REMOVE IN CASE OF NEGATIVE STAINING
        call micrograph_shrunken%stats( ave=ave, sdev=sdev, maxv=maxv, minv=minv)
        ! 3) Binarization
        if (detector .eq. 'bin') then
            self%detector = 'bin'
            call micrograph_shrunken%binarize(ave+.8*sdev)
        else if (detector .eq. 'otsu') then
            self%detector = 'otsu'
            call otsu_robust_fast(micrograph_shrunken, is2D=.true., noneg=.false., thresh=ttthresh)
            call micrograph_shrunken%erode() !morphological erosion
        else
          ! default self%detector is bin
           self%detector = 'bin'
           call micrograph_shrunken%binarize(ave+.8*sdev)
        endif
        if( DOWRITEIMGS ) call micrograph_shrunken%write(PATH_HERE//basename(trim(self%fbody))//'_bin.mrc')
        self%winsz = int(self%min_rad+self%max_rad)/4 !half of the avg between the dimensions of the particles
        call micrograph_shrunken%real_space_filter(self%winsz,'median') !median filtering allows easy calculation of cc
        if(DOWRITEIMGS) call micrograph_shrunken%write(PATH_HERE//basename(trim(self%fbody))//'_bin_median.mrc')
        ! 5) Connected components (cc) identification
        call micrograph_shrunken%find_ccs(self%img_cc,update_imat=.true.)
        ! 6) cc filtering
        call self%img_cc%polish_ccs([self%min_rad,self%max_rad],circular=' no',elongated=self%elongated, min_nccs = MIN_NCCS)
        if(DOWRITEIMGS) call self%img_cc%write_bimg(PATH_HERE//basename(trim(self%fbody))//'_CcsElimin.mrc')
        call micrograph_shrunken%kill_bimg
    contains

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
        write(unit = 17, fmt = "(a,i4)")    'Winsz for median filter ', self%winsz
        write(unit = 17, fmt = "(a,f4.2)")  'TV filter parameter     ', self%lambda
        write(unit = 17, fmt = "(a,a,a,a)") 'particle dimensions     ', trim(int2str(int(self%min_rad))),' ', trim(int2str(int(self%max_rad)))
        write(unit = 17, fmt = "(a,a)")     'detector                ', self%detector
        write(unit = 17, fmt = "(a,f4.2,f4.2)")  'min/max rad in A        ', self%min_rad*(SHRINK*self%smpd), self%max_rad*(SHRINK*self%smpd)
        close(17, status = "keep")
    end subroutine print_info

    ! This routine is aimed to eliminate aggregations of particles.
    ! It takes in input the list of the coordinates of the identified
    ! particles and the approximate radius of the particle.
    ! The output is a new list of coordinates where aggregate picked particles
    ! have been deleted
    ! If the dist beetween the 2 particles is < 2.*r -> delete them. (not too stringent)
    subroutine elimin_aggregation( self, saved_coord, msk )
        class(segpicker),     intent(inout)  :: self
        real,                 intent(in)     :: saved_coord(:,:)     !Coordinates of picked particles
        logical,              intent(inout)  :: msk(:)
        integer              :: i, j, cnt
        do i = 1, self%n_particles-1           !fix one coord
            do j = i+1, self%n_particles       !fix another coord to compare
                if(msk(i) .and. msk(j)) then !not compare twice ,and if the particles haven t been deleted yet
                    if( euclid(saved_coord(i,:), saved_coord(j,:)) <= self%distthr) then
                        msk(i) = .false.
                        msk(j) = .false.
                    endif
                endif
            enddo
        enddo
        ! TO CHECK IF NECESSARY
        call self%relative_intensity_filtering(msk)
        if(DEBUG_HERE) write(logfhandle,*) 'after  relative intensity filtering nb of particles is ', count(msk)
    end subroutine elimin_aggregation

    subroutine relative_intensity_filtering(self, selected_particle_positions)
        class(segpicker),     intent(inout) :: self
        logical,              intent(inout) :: selected_particle_positions(:)
        integer, allocatable :: imat(:,:,:)
        integer :: n_cc
        logical :: mask(self%ldim_shrunken(1),self%ldim_shrunken(2),1)
        real    :: avg_stdev,stdev_level
        imat = nint(self%img_cc%get_rmat())
        allocate(self%avg_gray_level(self%n_particles),self%stdev_gray_level(self%n_particles),source = 0.)
        do n_cc = 1, self%n_particles ! for each connected component
            where(imat == n_cc)
              mask = .true.
            elsewhere
              mask = .false.
            endwhere
            call calc_avgst_intensity_particle(mask,self%avg_gray_level(n_cc),self%stdev_gray_level(n_cc))
        enddo
        stdev_level = 0.
        avg_stdev   = sum(self%stdev_gray_level)/real(size(self%stdev_gray_level))
        do n_cc = 1, self%n_particles
            stdev_level = stdev_level + (avg_stdev- self%stdev_gray_level(n_cc))**2
        enddo
        stdev_level = sqrt(stdev_level/real(self%n_particles-1))
        ! Filter assuming gaussian distribution
        do n_cc = 1, self%n_particles
            if( self%stdev_gray_level(n_cc) < avg_stdev - 2.*stdev_level) then
              selected_particle_positions(n_cc) =  .false.
            else
              selected_particle_positions(n_cc) =  .true.
            endif
        enddo
    contains
        ! PUT IT TOGETHER WITH THE PREVIOUS FUNCTION
        subroutine calc_avgst_intensity_particle(mask,avg_gray_level,stdev_gray_level)
            logical, intent(in)  :: mask(:,:,:)
            real,    intent(out) :: avg_gray_level
            real,    intent(out) :: stdev_gray_level
            real,    allocatable :: rmat(:,:,:)
            integer :: i,j
            rmat = self%img%get_rmat()
            if(count(mask) == 1) then
                avg_gray_level = sum(rmat, mask) ! value of rmat in the mask
                stdev_level    = 0.
                return
            endif
            avg_gray_level   = 0.
            stdev_gray_level = 0.
            avg_gray_level   = sum(rmat, mask)/real(count(mask))
            do i = 1, self%ldim_shrunken(1)
                do j = 1, self%ldim_shrunken(2)
                    if(mask(i,j,1)) stdev_gray_level = stdev_gray_level + (avg_gray_level-rmat(i,j,1))**2
                enddo
            enddo
            stdev_gray_level = sqrt(stdev_gray_level/(real(count(mask)-1)))
        end subroutine calc_avgst_intensity_particle
    end subroutine relative_intensity_filtering


    ! This subroutine takes in input an image, its connected components image,
    ! and extract particles. It doesn't use mass_center.
    ! notation:: cc = connected component.
    subroutine identify_particle_positions(self)
      class(segpicker), intent(inout) :: self
      integer              :: n_cc
      integer              :: cnt
      real                 :: pos(3)            !center of each cc
      real,    allocatable :: xyz_saved(:,:)
      logical, allocatable :: mask(:)
      integer, allocatable :: imat(:,:,:), imat_cc(:,:,:)
      ! Initialisations
      self%orig_box = int(4.*(self%max_rad)+2.*self%max_rad) !needs to be bigger than the particle
      call self%img_cc%get_imat(imat_cc)
      allocate(xyz_saved(maxval(imat_cc),2), source = 0.) ! size of the # of cc (likely particles)
      allocate(imat(1:self%ldim_shrunken(1),1:self%ldim_shrunken(2),1:self%ldim_shrunken(3)), source = 0)
      ! Particle identification, extraction and centering
      where(imat_cc > 0.5) imat = 1 ! binary version
      do n_cc = 1, maxval(imat_cc)
          pos(:) = self%center_mass_cc(n_cc)
          xyz_saved(n_cc,:) = pos(:2)
      enddo
      deallocate(imat_cc)
      self%n_particles = size(xyz_saved, dim = 1) ! first estimation
      if(DEBUG_HERE) write(logfhandle,*) 'before elimin aggregations: ', size(xyz_saved, dim=1)
      allocate(mask(size(xyz_saved, dim=1)), source = .true. )
      call self%elimin_aggregation(xyz_saved, mask)
      allocate(self%particles_coord(count(mask),2), source = 0.)
      cnt = 0 !initialise
      do n_cc = 1,  size(xyz_saved, dim=1)
          if(mask(n_cc)) then
              cnt = cnt + 1
              self%particles_coord(cnt,:) = xyz_saved(n_cc,:)
          endif
      enddo
      self%n_particles = size(self%particles_coord, dim = 1) !update after elim aggregations
      if(DEBUG_HERE) write(logfhandle,*) 'after elimin aggregations: ', self%n_particles
      if(allocated(xyz_saved))       deallocate(xyz_saved)
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
        write(logfhandle,*) 'n_particles inside output particle position: ', self%n_particles
        open(unit = 23, file = basename(trim(self%fbody))//"Stdev.txt")
        write(unit = 23, fmt = "(a,f4.2)") basename(trim(self%fbody))
        call imgwin_particle%new([self%orig_box,self%orig_box,1],self%smpd)
        call self%img%rmsd(dev)
        write(unit = 23, fmt = "(a,f4.2)") 'Stdev ', dev
        do n_cc = 1, self%n_particles
            call self%img%window_slim(nint(self%particles_coord(n_cc,:)-self%orig_box/2), self%orig_box, imgwin_particle, outside(n_cc))
            if( .not. outside(n_cc)) then
                cnt = cnt + 1
                if(DEBUG_HERE) call imgwin_particle%write('centered_particles.mrc', cnt)
                write(unit = 23, fmt = "(a,i4,a,f4.2)") 'Particle ', cnt, 'Avg,  Stdev ', self%avg_gray_level(n_cc), self%stdev_gray_level(n_cc)
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
        if(DOWRITEIMGS) call self%img%write(PATH_HERE//basename(trim(self%fbody))//'_picked_particles.mrc')
    end subroutine output_identified_particle_positions

    ! This function returns the index of a pixel (assuming to have a 2D)
    ! image in a connected component. The pixel identified is the one
    ! center of mass of the cc.
    function center_mass_cc(self,n_cc) result (px)
        class(segpicker),  intent(inout) :: self
        integer,           intent(in)    :: n_cc
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

    subroutine write_boxfile(self, boxfile)
        class(segpicker),          intent(inout) :: self
        character(len=LONGSTRLEN), intent(out)   :: boxfile
        integer :: funit, n_cc,iostat
        self%particles_coord = (real(SHRINK)*self%particles_coord)-self%orig_box/2.
        call fopen(funit, status='REPLACE', action='WRITE', file=self%boxname, iostat=iostat)
        call fileiochk('picker; write_boxfile ', iostat)
        do n_cc=1,self%n_particles
            write(funit,'(I7,I7,I7,I7,I7)') int(self%particles_coord(n_cc,1)),&
            int(self%particles_coord(n_cc,2)), self%orig_box, self%orig_box, -3
        end do
        call fclose(funit,errmsg='picker; write_boxfile end')
        ! returns absolute path
        call make_relativepath(CWD_GLOB, self%boxname, boxfile)
    end subroutine write_boxfile

    subroutine kill(self)
        class(segpicker),  intent(inout) :: self
        !kill images
        call self%img%kill
        call self%img_cc%kill_bimg
        if(allocated(self%particles_coord)) deallocate(self%particles_coord)
        if(allocated(self%stdev_gray_level)) deallocate(self%stdev_gray_level)
        self%min_rad          = 0.
        self%max_rad          = 0.
        self%smpd             = 0.
        self%smpd_shrunken    = 0.
        self%hp_box           = 0.
        self%distthr          = 0.
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

! This function returns the index of a pixel (assuming to have a 2D)
! image in a connected component. The pixel identified is the one
! that minimizes the distance between itself and all the other pixels of
! the connected component. It corresponds to the geometric median.
! function center_cc(self,n_cc) result (px)
!     class(segpicker),  intent(inout) :: self
!     integer,               intent(in)    :: n_cc
!     real        :: px(3)               !index of the central px of the cc
!     integer     :: i, j, k
!     integer     :: n_px                !counter
!     integer     :: idx(2)
!     real,    allocatable :: dist(:,:)        !to extract the window according to the px which minimize the dist
!     logical, allocatable :: mask(:,:)        !to calc the min of an array along a specific dim
!     logical, allocatable :: mask_dist(:)     !for sum dist calculation
!     integer, allocatable :: pos(:,:)         !position of the pixels of a fixed cc
!     integer, allocatable :: imat_cc(:,:,:)
!     imat_cc = int(self%img_cc%get_rmat())
!     where(imat_cc .ne. n_cc) imat_cc = 0
!     call get_pixel_pos(imat_cc,pos)
!     allocate(dist(4,   size(pos, dim = 2)), source = 0.)
!     allocate(mask(4,   size(pos, dim = 2)), source = .false.)
!     allocate(mask_dist(size(pos, dim = 2)), source = .true.)
!     mask(4,:) = .true. !to calc the min wrt the dist
!     n_px = 0
!     do i = 1, self%ldim_shrunken(1)
!         do j = 1, self%ldim_shrunken(2)
!             do k = 1, self%ldim_shrunken(3)
!                 if(imat_cc(i,j,k) > 0.5) then
!                     n_px = n_px + 1
!                     dist( 4, n_px) = pixels_dist([i,j,k], pos, 'sum', mask_dist)
!                     dist(:3, n_px) = [real(i),real(j),real(k)]
!                 endif
!             enddo
!         enddo
!     enddo
!     idx   = minloc(dist, mask)
!     px(:) = dist(1:3, idx(2))
!     deallocate(pos, mask, dist, mask_dist)
!     if(allocated(imat_cc)) deallocate(imat_cc)
! end function center_cc
