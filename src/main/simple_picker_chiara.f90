! USAGE: simple_private_exec prg=pick_chiara detector=bin smpd=1. part_radius=15 fname='/home/chiara/Desktop/Chiara/ParticlePICKING/PickingResults/SomeExamples/NegativeSTORIGINAL.mrc'
module simple_picker_chiara
include 'simple_lib.f08'
use simple_image, only : image
implicit none

 public :: extract_particles, preprocess_mic, picker_chiara, print_info

#include "simple_local_flags.inc"

! module global constants
real, parameter :: SHRINK  = 4.

type :: picker_chiara
    private
    type(image) :: img
    type(image) :: img_cc
    real, allocatable :: particles_coord(:,:)
    real    :: part_radius           = 0.
    real    :: smpd                  = 0.
    real    :: smpd_shrunken         = 0.
    real    :: hp_box                = 0.
    integer :: ldim(3)               = 0
    integer :: ldim_shrunken(3)      = 0
    integer :: n_particles           = 0
    character(len=STDLEN) :: pickername = ''   !fname
    character(len=STDLEN) :: fbody      = ''   !fbody
    character(len=STDLEN) :: detector   = ''
    real     :: lambda ! for tv denoising
    real     :: lp     ! low pass filtering
contains
    ! constructor
    procedure          :: new => new_picker
    ! picking functions
    procedure          :: extract_particles
    procedure          :: elimin_aggregation
    procedure          :: center_cc
    ! preprocess mic prior picking
    procedure          :: preprocess_mic
    ! setters/getters

    ! output
    procedure          :: print_info
    ! kill
    procedure          :: kill => kill_picker
end type picker_chiara

private

contains

    subroutine new_picker(self, fname, radius, smpd)
        class(picker_chiara),  intent(inout) :: self
        character(len=*),      intent(in)    :: fname
        real,                  intent(in)    :: radius
        real,                  intent(in)    :: smpd
        integer :: nptcls
        self%pickername = fname
        call find_ldim_nptcls(self%pickername, self%ldim, nptcls, self%smpd)
        self%smpd = smpd
        self%ldim_shrunken(1) = round2even(real(self%ldim(1))/SHRINK)
        self%ldim_shrunken(2) = round2even(real(self%ldim(2))/SHRINK)
        self%ldim_shrunken(3) = 1
        self%smpd_shrunken = self%smpd*SHRINK
        self%part_radius = radius
        self%fbody = get_fbody(trim(fname), trim(fname2ext(fname)))
        call self%img%new   (self%ldim_shrunken, self%smpd_shrunken)
        call self%img_cc%new(self%ldim_shrunken, self%smpd_shrunken)
        self%n_particles = 0
        self%lambda      = 5.
        self%lp          = 20.
        self%detector    = 'bin'
    end subroutine new_picker

    subroutine preprocess_mic(self, detector,lp)
        use simple_tvfilter, only : tvfilter
        use simple_micops
        use simple_segmentation, only: sobel, automatic_thresh_sobel
        class(picker_chiara), intent(inout) :: self
        character(len= *),    intent(in)    :: detector
        real, optional,       intent(in)    :: lp
        real, allocatable :: rmat(:,:,:)
        real, allocatable :: x(:), x_out(:)
        type(tvfilter)    :: tvf
        type(image) :: mic_copy
        integer     :: box_shrunken, winsz
        real        :: ave, sdev, maxv, minv
        real        :: thresh(1)
        ! 0) Reading and saving original micrograph
        call read_micrograph(self%pickername, smpd = self%smpd)
        ! 1) Shrink and high pass filtering
        call shrink_micrograph(SHRINK, self%ldim_shrunken, self%smpd_shrunken)
        self%hp_box =  4.*self%part_radius+2.*self%part_radius
        call set_box(int(SHRINK*(self%hp_box)), box_shrunken)
        ! To take care of shrinking
        self%part_radius = self%part_radius/SHRINK ! I am thingking I shouldn't multiply by the smpd cuz I am working in pxls
        call mic_copy%new(self%ldim_shrunken, self%smpd_shrunken)
        call self%img%read('shrunken_hpassfiltered.mrc')
        call mic_copy%read('shrunken_hpassfiltered.mrc')
        ! 2) Low pass filtering
        if(present(lp)) self%lp = lp
        call mic_copy%bp(0.,self%lp)
        call mic_copy%ifft()
        ! 2.1) TV denoising
        call tvf%new()
        call tvf%apply_filter(mic_copy, self%lambda)
        call tvf%kill
        call mic_copy%write(trim(self%fbody)//'_tvfiltered.mrc')
        ! 2.3) negative image, to obtain a binarization with white particles
        call mic_copy%neg() !TO REMOVE IN CASE OF NEGATIVE STAINING
        ! 3) Binarization
        call mic_copy%stats( ave, sdev, maxv, minv )
        if(detector .eq. 'sobel') then
            self%detector = 'sobel'
            thresh(1) = ave+.5*sdev !sobel needs lower thresh not to pick just edges
            call sobel(mic_copy,thresh)
        else if (detector .eq. 'bin') then
            ! default self%detector is bin
            call mic_copy%bin(ave+.8*sdev)
        else if (detector .eq. 'otsu') then
            self%detector = 'otsu'
            rmat = mic_copy%get_rmat()
            x = pack(rmat, .true.)
            call otsu(x,x_out)
            rmat = reshape(x_out, [self%ldim_shrunken(1),self%ldim_shrunken(2),1])
            call mic_copy%set_rmat(rmat)
            deallocate(x,x_out,rmat)
        else
            THROW_HARD('Invalid detector; preprocess_mic')
        endif
        call mic_copy%write(trim(self%fbody)//'_Bin.mrc')
        winsz = int(self%part_radius)/2
        call mic_copy%real_space_filter(winsz,'median') !median filtering allows easy calculation of cc
        call mic_copy%write(trim(self%fbody)//'_BinMedian.mrc')
        ! 5) Connected components (cc) identification
        call mic_copy%find_connected_comps(self%img_cc)
        ! 6) cc filtering
        call self%img_cc%polish_cc(self%part_radius)
        call self%img_cc%write(trim(self%fbody)//'_ConnectedComponentsElimin.mrc')
        call mic_copy%kill
    end subroutine preprocess_mic

    subroutine print_info(self)
        class(picker_chiara), intent(inout) :: self
        open(unit = 17, file = "PickerInfo.txt")
        write(unit = 17, fmt = '(a)') '>>>>>>>>>>>>>>>>>>>>PARTICLE PICKING>>>>>>>>>>>>>>>>>>'
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a,f0.0)")             'Mic Shrunken, fact ', SHRINK
        write(unit = 17, fmt = "(a,i4,tr1,i4,tr1,i4)") 'Dim before  shrink ', self%ldim
        write(unit = 17, fmt = "(a,i4,tr1,i4,tr1,i4)") 'Dim after   shrink ', self%ldim_shrunken
        write(unit = 17, fmt = "(a,f4.2)")             'Smpd before shrink ', self%smpd
        write(unit = 17, fmt = "(a,f4.2)")             'Smpd after  shrink ', self%smpd_shrunken
        write(unit = 17, fmt = "(a,a)")                'Hp box              ', trim(int2str(int(self%hp_box)))
        write(unit = 17, fmt = "(a,i4,tr1,i4)")        'Ccs size filtering ', int(5*self%part_radius), int(2*3.14*(self%part_radius)**2)
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a)")  'SELECTED PARAMETERS '
        write(unit = 17, fmt = '(a)') ''
        write(unit = 17, fmt = "(a,f0.0)")  'Lp filter parameter ', self%lp
        write(unit = 17, fmt = "(a,f0.0)")  'TV filter parameter ', self%lambda
        write(unit = 17, fmt = "(a,a)")     'part_radius         ', trim(int2str(int(self%part_radius)))
        write(unit = 17, fmt = "(a,a)")     'detector            ', self%detector
        close(17, status = "keep")
    end subroutine print_info

    ! This routine is aimed to eliminate aggregations of particles.
    ! It takes in input the list of the coordinates of the identified
    ! particles and the approximate radius of the particle.
    ! The output is a new list of coordinates where aggregate picked particles
    ! have been deleted
    ! If the dist beetween the 2 particles is in [r,2r] -> delete them.
    subroutine elimin_aggregation( self, saved_coord, refined_coords )
        class(picker_chiara),  intent(inout) :: self
        real,                 intent(in)     :: saved_coord(:,:)     !Coordinates of picked particles
        real, allocatable,    intent(out)    :: refined_coords(:,:)  !New coordinates of not aggregated particles
        logical, allocatable :: msk(:)
        integer              :: i, j, cnt
        allocate(msk(size(saved_coord,dim=1)), source = .true.) !initialise to 'keep all'
        cnt = 0
        do i = 1, size(saved_coord, dim = 1)             !fix one coord
            do j = i+1, size(saved_coord, dim = 1)       !fix another coord to compare
                if(msk(i) .and. msk(j)) then !not compare twice ,and if the particles haven t been deleted yet
                    if(  sqrt(real((saved_coord(i,1)-saved_coord(j,1))**2 &
                        &        + (saved_coord(i,2)-saved_coord(j,2))**2)) <= 2.*self%part_radius) then !&
                        !& .and. self%part_radius < &
                        !& real((saved_coord(i,1)-saved_coord(j,1))**2 &        ! TO CHECK WHAT TO DO
                        !&    + (saved_coord(i,2)-saved_coord(j,2))**2)) then
                        msk(i) = .false.
                        msk(j) = .false.
                        cnt = cnt + 1                    !number of deleted couples
                    endif
                endif
            enddo
        enddo
        allocate(refined_coords(size(saved_coord, dim = 1)-cnt,2), source = 0.)
        cnt = 0
        do i = 1, size(saved_coord, dim = 1)
            if(msk(i)) then
                cnt = cnt + 1
                refined_coords(cnt, :) = saved_coord(i, :)
            endif
        enddo
        deallocate(msk)
    end subroutine elimin_aggregation


  ! This subroutine takes in input an image, its connected components image,
  ! and extract particles. It doesn't use mass_center.
  ! notation:: cc = connected component.
  subroutine extract_particles(self)
      class(picker_chiara), intent(inout) :: self
      type(image)          :: imgwin_particle
      type(image)          :: img_back
      integer              :: box
      integer              :: n_cc
      integer              :: cnt
      real                 :: pos(3)            !center of each cc
      real, allocatable    :: xyz_saved(:,:), xyz_norep_noagg(:,:)
      integer, allocatable :: imat_cc(:,:,:)
      integer, allocatable :: imat(:,:,:)
      logical              :: outside
      ! Initialisations
      box     = int(4.*(self%part_radius)+2.*self%part_radius) !needs to be bigger than the particle
      imat_cc = int(self%img_cc%get_rmat())
      call imgwin_particle%new([box,box,1],self%smpd)
      allocate(xyz_saved(maxval(imat_cc),2), source = 0.) ! size of the # of cc (likely particles)
      allocate(imat(1:self%ldim_shrunken(1),1:self%ldim_shrunken(2),1:self%ldim_shrunken(3)), source = 0)
      ! Copy of the micrograph, to highlight on it the picked particles
      call img_back%copy(self%img)
      ! Particle identification, extraction and centering
      where(imat_cc > 0.5) imat = 1
      do n_cc = 1, maxval(imat_cc)
          pos(:) = self%center_cc(n_cc)
          xyz_saved(n_cc,:) = pos(:2)
      enddo
      deallocate(imat_cc)
      call self%elimin_aggregation(xyz_saved, xyz_norep_noagg)      ! x2 not to pick too close particles
      allocate(self%particles_coord(count(xyz_norep_noagg(:,1)> TINY),2), source = 0.) !TO IMPROVE
      cnt = 0 !initialise
      do n_cc = 1,  size(xyz_norep_noagg, dim=1)
          if(abs(xyz_norep_noagg(n_cc,1)) > TINY .and. abs(xyz_norep_noagg(n_cc,2))> TINY) then
              cnt = cnt + 1
              self%particles_coord(cnt,:) = xyz_norep_noagg(n_cc,:)
          endif
      enddo
      self%n_particles = size(self%particles_coord, dim = 1)
      outside = .false.
      cnt = 0
      do n_cc = 1, self%n_particles
          !if( abs(self%particles_coord(n_cc,1)) > TINY) then !useless?
              call img_back%draw_picked(nint(self%particles_coord(n_cc,:)),nint(self%part_radius),3, 'white')
              call self%img%window_slim(nint(self%particles_coord(n_cc,:)-box/2), box, imgwin_particle, outside)
              if( .not. outside) then
                  cnt = cnt + 1
                  call imgwin_particle%write(trim(self%fbody)//'_centered_particles.mrc', cnt)

              endif
          !endif
      end do
      call img_back%write(trim(self%fbody)//'_picked_particles.mrc')
      deallocate(xyz_saved,xyz_norep_noagg)
  end subroutine extract_particles

    ! This function returns the index of a pixel (assuming to have a 2D)
    ! image in a connected component. The pixel identified is the one
    ! that minimizes the distance between itself and all the other pixels of
    ! the connected component. It corresponds to the geometric median.
    function center_cc(self,n_cc) result (px)
        class(picker_chiara),  intent(inout) :: self
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

    subroutine kill_picker(self)
        class(picker_chiara),  intent(inout) :: self
        !kill images
        call self%img%kill
        call self%img_cc%kill
        if(allocated(self%particles_coord)) deallocate(self%particles_coord)
        self%part_radius      = 0.
        self%smpd             = 0.
        self%smpd_shrunken    = 0.
        self%hp_box           = 0.
        self%ldim(:)          = 0
        self%ldim_shrunken(:) = 0
        self%n_particles      = 0
        self%pickername       = ''   !fname
        self%fbody            = ''   !fbody
        self%detector         = ''
        self%lp               = 0.
        self%lambda           = 0.
    end subroutine kill_picker
end module simple_picker_chiara
