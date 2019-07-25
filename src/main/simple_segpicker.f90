module simple_segpicker
include 'simple_lib.f08'
use simple_image, only : image
use simple_parameters, only: params_glob
implicit none

 public :: segpicker

#include "simple_local_flags.inc"

! module global constants
real,    parameter :: SHRINK  = 4.
logical, parameter :: DEBUG_HERE = .true.

type :: segpicker
    private
    type(image) :: img
    type(image) :: img_cc
    real, allocatable :: particles_coord(:,:)
    real    :: min_rad               = 0.
    real    :: max_rad               = 0.
    real    :: smpd                  = 0.
    real    :: smpd_shrunken         = 0.
    real    :: hp_box                = 0.
    real    :: lambda                = 10. ! for tv denoising
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
    procedure, private :: center_cc
    procedure, private :: center_mass_cc
    ! preprocess mic prior picking
    procedure          :: preprocess_mic
    ! setters/getters

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
        self%boxname =  trim(self%fbody)//'.box'
        call self%img%new   (self%ldim_shrunken, self%smpd_shrunken)
        call self%img_cc%new(self%ldim_shrunken, self%smpd_shrunken)
        self%n_particles = 0
        self%lambda      = 5.
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
        integer     :: box_shrunken, winsz
        real        :: ave, sdev, maxv, minv
        real        :: thresh(1)
        ! 0) Reading and saving original micrograph
        call read_micrograph(self%pickername, smpd=self%smpd)
        ! 1) Shrink and high pass filtering
        call shrink_micrograph(SHRINK, self%ldim_shrunken, self%smpd_shrunken)
        self%hp_box =  4.*self%max_rad+2.*self%max_rad
        call set_box(int(SHRINK*(self%hp_box)), box_shrunken, micrograph_shrunken)
        call self%img%copy(micrograph_shrunken)
        ! To take care of shrinking
        print *, 'self%min_rad = ', self%min_rad, 'self%max_rad = ', self%max_rad
        self%min_rad = self%min_rad/SHRINK ! I am thingking I shouldn't multiply by the smpd cuz I am working in pxls
        self%max_rad = self%max_rad/SHRINK
        ! 2) Low pass filtering
        call micrograph_shrunken%bp(0., params_glob%lp)
        call micrograph_shrunken%ifft()
        ! 2.1) TV denoising
        call tvf%new()
        print *, 'tv filtering, with lambda ' , self%lambda
        call tvf%apply_filter(micrograph_shrunken, self%lambda)
        call tvf%kill
        call micrograph_shrunken%write(trim(self%fbody)//'_tvfiltered.mrc')
        ! 2.2) negative image, to obtain a binarization with white particles
        call micrograph_shrunken%neg() !TO REMOVE IN CASE OF NEGATIVE STAINING
        ! 3) Binarization
        call micrograph_shrunken%stats( ave, sdev, maxv, minv )
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
            deallocate(x,x_out,rmat)
        else
            THROW_HARD('Invalid detector; preprocess_mic')
        endif
        if( DEBUG_HERE ) call micrograph_shrunken%write(trim(self%fbody)//'_Bin.mrc')
        winsz = int(self%min_rad+self%max_rad)/4 !half of the avg between the dimensions of the particles
        call micrograph_shrunken%real_space_filter(winsz,'median') !median filtering allows easy calculation of cc
        if( DEBUG_HERE ) call micrograph_shrunken%write(trim(self%fbody)//'_BinMedian.mrc')
        ! 5) Connected components (cc) identification
        call micrograph_shrunken%find_connected_comps(self%img_cc)
        ! 6) cc filtering
        call self%img_cc%polish_cc(self%min_rad,self%max_rad)
        if( DEBUG_HERE ) call self%img_cc%write(trim(self%fbody)//'_ConnectedComponentsElimin.mrc')
        call micrograph_shrunken%kill
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
        write(unit = 17, fmt = "(a,f0.0)")  'Lp filter parameter ', params_glob%lp
        write(unit = 17, fmt = "(a,f0.0)")  'TV filter parameter ', self%lambda
        write(unit = 17, fmt = "(a,a,a,a)") 'particle dimensions ', trim(int2str(int(self%min_rad))),' ', trim(int2str(int(self%max_rad)))
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
        class(segpicker),  intent(inout) :: self
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
                        &        + (saved_coord(i,2)-saved_coord(j,2))**2)) <= 2.*self%min_rad) then!&
                        !& .and. self%min_rad < &
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
      deallocate(xyz_saved,xyz_norep_noagg)
    end subroutine identify_particle_positions

      subroutine output_identified_particle_positions(self)
          class(segpicker), intent(inout) :: self
          type(image) :: imgwin_particle
          integer :: n_cc
          integer :: cnt
          logical :: outside
          outside = .false.
          cnt = 0
          call imgwin_particle%new([self%orig_box,self%orig_box,1],self%smpd)
          do n_cc = 1, self%n_particles
              call self%img%window_slim(nint(self%particles_coord(n_cc,:)-self%orig_box/2), self%orig_box, imgwin_particle, outside)
              if( .not. outside) then
                  cnt = cnt + 1
                  if( DEBUG_HERE )call imgwin_particle%write(trim(self%fbody)//'_centered_particles.mrc', cnt)
              endif
          enddo
          do n_cc = 1, self%n_particles
              call self%img%draw_picked(nint(self%particles_coord(n_cc,:)),nint((self%min_rad+self%max_rad)/2.),2, self%color)
          end do
          call imgwin_particle%kill
          if( DEBUG_HERE ) call self%img%write(trim(self%fbody)//'_picked_particles.mrc')
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

    ! subroutine write_boxfile_seg(self)
    !     class(segpicker),  intent(inout) :: self
    !     integer :: n_cc, iostat
    !     open(121, file=self%boxname, iostat=iostat)
    !     call fileiochk('segpicker; write_boxfile_seg ', iostat)
    !     print *, 'Before loop'
    !     print *, 'self%n_particles', self%n_particles
    !     do n_cc=1,self%n_particles
    !         print *, 'n_cc = ', n_cc
    !         print *, self%particles_coord(n_cc,1),&
    !         & self%particles_coord(n_cc,2), self%orig_box, self%orig_box, -3
    !         write(unit = 121, fmt = '(I7,I7,I7,I7,I7)') self%particles_coord(n_cc,1),&
    !                                  & self%particles_coord(n_cc,2), self%orig_box, self%orig_box, -3
    !     end do
    !     print *, 'After loop'
    !     close(121)
    ! end subroutine write_boxfile_seg


    subroutine write_boxfile(self)
        class(segpicker),  intent(inout) :: self
        integer :: funit, n_cc,iostat
        print *, 'funit = ', funit
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
        if(allocated(self%particles_coord)) deallocate(self%particles_coord)
        self%min_rad          = 0.
        self%max_rad          = 0.
        self%smpd             = 0.
        self%smpd_shrunken    = 0.
        self%hp_box           = 0.
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
