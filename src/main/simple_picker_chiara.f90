module simple_picker_chiara
include 'simple_lib.f08'
use simple_image, only : image
implicit none

! public :: extract_particles, center_cc, discard_borders !maybe to remove center_cc from public
!
! private
#include "simple_local_flags.inc"
contains

    ! This function takes in input an image and gives in output another image
    ! with smaller size in which are discarded about 4% of the borders. It
    ! is meantt o deal with border effects.
    subroutine discard_borders(img_in, img_out, reduce_ldim)
        class(image),      intent(in)  :: img_in
        class(image),      intent(out) :: img_out
        integer, optional, intent(out) :: reduce_ldim(3)
        logical :: outside
        integer :: ldim(3), rreduce_ldim(3), box
        real    :: smpd
        outside = .false. !initialization
        ldim = img_in%get_ldim()
        smpd = img_in%get_smpd()
        rreduce_ldim = nint(4.*minval(ldim(:2))/100.)
        rreduce_ldim(3) = 1
        if(present(reduce_ldim)) reduce_ldim = rreduce_ldim
        box = minval(ldim(:2))- 2*rreduce_ldim(1)
        call img_out%new([box,box,1], smpd)
        call img_in%window_slim(rreduce_ldim(:2),box, img_out, outside)
        if(outside) THROW_HARD('It is outside! discard_borders')
    end subroutine discard_borders

    ! This function tells whether the new_coord of a likely particle
    ! have already been identified.
    function is_picked( new_coord, saved_coord, part_radius, saved ) result(yes_no)
        integer, intent(in) :: new_coord(2)       !Coordinates of a new particle to pick
        integer, intent(in) :: saved_coord(:,:)   !Coordinates of picked particles
        integer, intent(in) :: part_radius        !Approximate radius of the particle
        integer, intent(in), optional :: saved    !How many particles have already been saved
        logical :: yes_no
        integer :: iwind, ssaved, s(2)
        s = shape(saved_coord)
        if(s(2) /= 2) THROW_HARD('dimension error; is_picked')
        yes_no = .false.
        ssaved = s(1)
        if(present(saved)) ssaved = saved
        do iwind = 1, ssaved
            if(  sqrt(real((new_coord(1)-saved_coord(iwind,1))**2 + (new_coord(2)-saved_coord(iwind,2)) **2)) <= real(part_radius)) then
                yes_no = .true.
                return
            endif
        enddo
    end function is_picked

    ! This routine is aimed to eliminate aggregations of particles.
    ! It takes in input the list of the coordinates of the identified
    ! particles and the approximate radius of the particle.
    ! The output is a new list of coordinates where aggregate picked particles
    ! have been deleted
    !if the dist beetween the 2 particles is in [r,2r], delete them
    subroutine elimin_aggregation( saved_coord, part_radius, refined_coords )
        integer,              intent(in)  :: saved_coord(:,:)     !Coordinates of picked particles
        integer,              intent(in)  :: part_radius          !Approximate radius of the particle
        integer, allocatable, intent(out) :: refined_coords(:,:)  !New coordinates of not aggregated particles
        logical, allocatable :: msk(:)
        integer              :: i, j, cnt
        allocate(msk(size(saved_coord, dim = 1)), source = .true.) !initialise to 'keep all'
        cnt = 0
        do i = 1, size(saved_coord, dim = 1)             !fix one coord
            do j = 1, size(saved_coord, dim = 1)         !fix another coord to compare
                if(j > i .and. msk(i) .and. msk(j)) then !not compare twice ,and if the particles haven t been deleted yet
                    if(  sqrt(real((saved_coord(i,1)-saved_coord(j,1))**2 &
                        & + (saved_coord(i,2)-saved_coord(j,2))**2)) <= real(2*part_radius) &
                        & .and. real(part_radius) < &
                        & real((saved_coord(i,1)-saved_coord(j,1))**2 &
                        &    + (saved_coord(i,2)-saved_coord(j,2))**2)) then
                        msk(i) = .false.
                        msk(j) = .false.
                        cnt = cnt + 1                    !number of deleted couples
                    endif
                endif
            enddo
        enddo
        allocate(refined_coords(size(saved_coord, dim = 1)-cnt,2), source = 0)
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
  subroutine extract_particles(self, img_cc, part_radius)
      type(image),  intent(inout) :: self        !original image
      type(image),  intent(inout) :: img_cc      !connected components image
      integer,      intent(in)    :: part_radius !extimated particle radius, to discard aggregations
      type(image)          :: imgwin_particle!, imgwin_bin
      type(image)          :: self_back
      integer              :: box
      integer              :: n_cc                            !n_cc  = # cc in the input image
      integer              :: cnt_particle, cnt_likely        !counters
      integer              :: xyz(3)                          !mass center coordinates
      integer :: ldim(3)
      real    :: pos(3)                                       !central position of each cc
      real    :: ave, sdev, maxv, minv, med                   !stats
      real    :: aveB, sdevB, maxvB, minvB, medB      !stats
      ! character (len = *), parameter :: fmt_1 = "(a)" !I/O
      integer, allocatable :: xyz_saved(:,:), xyz_norep_noagg(:,:)
      real,    allocatable :: rmat(:,:,:)                     !matrix corresponding to the cc image
      integer, allocatable :: rmat_masked(:,:,:)              !to identify each cc
      logical              :: outside, discard, errout, erroutB
      ! Initialisations
      ldim = self%get_ldim()
      box  = (part_radius)*4 + 10
      rmat = img_cc%get_rmat()
      call imgwin_particle%new([box,box,1],1.)
      !call imgwin_bin%new([box,box,1],1.)
      outside = .false.
      allocate(xyz_saved(int(maxval(rmat)),2), source = 0) ! int(maxval(rmat)) is the # of cc (likely particles)
      allocate(rmat_masked(ldim(1),ldim(2),1))
      ! Copy of the micrograph, to highlight on it the picked particles
      call self_back%copy(self)
      ! Particle identification, extraction and centering
      !open(unit = 17, file = "Statistics.txt")
      cnt_likely = 0
      do n_cc = 1, int(maxval(rmat))   !fix cc
          rmat_masked = 0              !reset
          if(any(abs(int(rmat)-n_cc)< TINY)) then
              where((abs(int(rmat)-n_cc)) < TINY) rmat_masked = 1
              cnt_likely = cnt_likely + 1
              pos = center_cc(rmat_masked)
              xyz_saved(cnt_likely,:) = int(pos(:))
          endif
      enddo
      call elimin_aggregation(xyz_saved, part_radius  + 0, xyz_norep_noagg) !+ 7 not to pick too close particles
      cnt_particle = 0
      do n_cc = 1, cnt_likely
          if( abs(xyz_norep_noagg(n_cc,1)) >  0) then !useless?
              if( abs(xyz_norep_noagg(n_cc,1)) >  0) call self_back%draw_picked( xyz_norep_noagg(n_cc,:),part_radius,3)
              call   self%window_slim(xyz_norep_noagg(n_cc,:)-box/2, box, imgwin_particle, outside)
              !call img_cc%window_slim(xyz_norep_noagg(n_cc,:)-box/2, box, imgwin_bin, outside)
              if( .not. outside) then
                  cnt_particle = cnt_particle + 1
                  call imgwin_particle%write('centered_particles.mrc', cnt_particle)
                  ! call imgwin_bin%write('centered_particles_BIN.mrc',  cnt_particle)
                  ! BINARY IMAGE STATS
                  !   TO START AGAIN FROM HEREEEEE
                  ! write(unit = 17, fmt = "(a,i0,2(a,f0.0))") &
                  ! &'image=', cnt_particle,' diameter=', part_longest_dim(imgwin_bin), ' nforeground=', imgwin_bin%nforeground()
              endif
          endif
      end do
      call self_back%write('picked_particles.mrc')
      deallocate(xyz_saved, xyz_norep_noagg, rmat)
    !  close(17, status = "keep")
  end subroutine extract_particles

    ! This function returns the index of a pixel (assuming to have a 2D)
    ! image in a connected component. The pixel identified is the one
    ! that minimizes the distance between itself and all the other pixels of
    ! the connected component.
    function center_cc(masked_mat) result (px)
        integer, intent(in) :: masked_mat(:,:,:)
        integer     :: px(3)               !index of the central px of the cc
        integer     :: i, j, k
        integer     :: s(3)                !shape of the input matrix
        integer     :: n_px                !counter
        integer     :: idx(2)
        real,    allocatable :: dist(:,:)        !to extract the window according to the px which minimize the dist
        logical, allocatable :: mask(:,:)        !to calc the min of an array along a specific dim
        integer, allocatable :: pos(:,:)         !position of the pixels of a fixed cc
        s = shape(masked_mat)
        call get_pixel_pos(masked_mat, pos)
        if(size(pos, dim = 2) < 1) print *, 'SHIT'
        allocate(dist(4,size(pos, dim = 2)), source = 0.)
        allocate(mask(4,size(pos, dim = 2)), source = .false.)
        mask(4,:) = .true.
        n_px = 0
        do i = 1, s(1)
            do j = 1, s(2)
                do k = 1, s(3)
                    if(masked_mat(i,j,k) > 0.5) then
                        n_px = n_px + 1
                        dist( 4, n_px) = pixels_dist( [i,j,k], pos, 'sum')
                        dist(:3, n_px) = [real(i),real(j),real(k)]
                    endif
                enddo
            enddo
        enddo
        idx   = minloc(dist, mask)
        px(:) = int(dist(1:3, idx(2)))
        deallocate(pos, mask, dist)
    end function center_cc

  ! This subroutine stores in pos the indeces corresponding to
  ! the pixels with value > 0 in the binary matrix rmat_masked.
  subroutine get_pixel_pos(rmat_masked, pos)
      integer,              intent(in)  :: rmat_masked(:,:,:)
      integer, allocatable, intent(out) :: pos(:,:)
      integer :: s(3), i, j, k, cnt
      ! if( any(rmat_masked > 1.0001) .or. any(rmat_masked < 0. ))&
      ! &THROW_HARD('Input not binary; get_pixel_pos')
      s = shape(rmat_masked)
      allocate(pos(3, count(rmat_masked > 0.5)), source = 0)
      cnt = 0
      do i = 1, s(1)
            do j = 1, s(2)
                do k = 1, s(3)
                    if(rmat_masked(i,j,k) > 0.5) then !rmat_masked is binary
                        cnt = cnt + 1
                        pos(:,cnt) = [i,j,k]
                    endif
                enddo
            enddo
        enddo
    end subroutine get_pixel_pos

    !>   calculates the euclidean distance between one pixel and a list of other pixels.
    ! if which == 'max' then distance is the maximum value of the distance between
    !              the selected pixel and all the others
    ! if which == 'min' then distance is the minimum value of the distance between
    !              the selected pixel and all the others
    ! if which == 'sum' then distance is the sum of the distances between the
    !              selected pixel and all the others.
    function pixels_dist( px, vec, which) result( dist )
        integer, intent(in)           :: px(3)
        integer, intent(in)           :: vec(:,:)
        character(len=*),  intent(in) :: which
        logical :: mask(size(vec, dim = 2))
        real    :: dist
        integer :: i
        mask(:) = .true. !to calculation of the 'min' excluding the pixel itself, otherwise it d always be 0
        do i = 1, size(vec, dim = 2)
            if( px(1)==vec(1,i) .and. px(2)==vec(2,i) .and. px(3)==vec(3,i) ) mask(i) = .false.
        enddo
        select case(which)
        case('max')
            dist =  maxval(sqrt((px(1)-vec(1,:))**2.+(px(2)-vec(2,:))**2.+(px(3)-vec(3,:))**2.))
        case('min')
            dist =  minval(sqrt((px(1)-vec(1,:))**2.+(px(2)-vec(2,:))**2.+(px(3)-vec(3,:))**2.), mask)
        case('sum')
            dist =  maxval(sqrt((px(1)-vec(1,:))**2.+(px(2)-vec(2,:))**2.+(px(3)-vec(3,:))**2.))
        case DEFAULT
            write(logfhandle,*) 'Pixels_dist kind: ', trim(which)
            THROW_HARD('Unsupported pixels_dist kind; pixels_dist')
        end select
    end function pixels_dist

    function points_dist( px, vec, which) result( dist )
        real,    intent(in)           :: px(3)
        real, intent(in)           :: vec(:,:)
        character(len=*),  intent(in) :: which
        logical :: mask(size(vec, dim = 2))
        real    :: dist
        integer :: i
        mask(:) = .true. !to calculation of the 'min' excluding the pixel itself, otherwise it d always be 0
        do i = 1, size(vec, dim = 2)
            if(      abs(px(1)-vec(1,i)) < TINY .and. abs(px(2)-vec(2,i)) < TINY  &
            &  .and. abs(px(3)-vec(3,i)) < TINY ) mask(i) = .false.
        enddo
        select case(which)
        case('max')
            dist =  maxval(sqrt((px(1)-vec(1,:))**2.+(px(2)-vec(2,:))**2.+(px(3)-vec(3,:))**2.))
        case('min')
            dist =  minval(sqrt((px(1)-vec(1,:))**2.+(px(2)-vec(2,:))**2.+(px(3)-vec(3,:))**2.), mask)
        case('sum')
            dist =  sum(sqrt((px(1)-vec(1,:))**2.+(px(2)-vec(2,:))**2.+(px(3)-vec(3,:))**2.))
        case DEFAULT
            write(logfhandle,*) 'Points_dist kind: ', trim(which)
            THROW_HARD('Unsupported pixels_dist kind; pixels_dist')
        end select
    end function points_dist

    ! This function calculates the longest dimension of the biggest connected
    ! component in the image img. Img is supposed to be either a binary
    ! image or a connected component image.
    ! The connected component whose longest dimension is calculated is
    ! supposed to be circular.
    function part_longest_dim(img) result(diameter)
        type(image), intent(inout) :: img
        real, allocatable    :: rmat(:,:,:), dist(:)
        integer, allocatable :: rmat_masked(:,:,:), pos(:,:)
        real :: diameter
        integer :: ldim(3), i, j,k, cnt
        ldim = img%get_ldim()
        rmat = img%get_rmat()
        call img%prepare_connected_comps() !consider just one particle, not multiple
        allocate(rmat_masked(ldim(1),ldim(2),ldim(3)), source = 0)
        rmat_masked = int(img%get_rmat())
        call get_pixel_pos(rmat_masked,pos)
        cnt = 0
        allocate(dist(count(rmat_masked > 0.5))) !one for each white px
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                do k = 1, ldim(3)
                    if(rmat_masked(i,j,k) > 0.5) then  !just ones
                        cnt = cnt + 1
                        dist(cnt) = pixels_dist([i,j,k],pos, 'max')
                    endif
                enddo
            enddo
        enddo
        diameter = maxval(dist)
    end function part_longest_dim
end module simple_picker_chiara
