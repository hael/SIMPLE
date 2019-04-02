! USAGE: simple_private_exec prg=pick_chiara detector=bin smpd=1. part_radius=15 fname='/home/chiara/Desktop/Chiara/ParticlePICKING/PickingResults/SomeExamples/NegativeSTORIGINAL.mrc'
module simple_picker_chiara
include 'simple_lib.f08'
use simple_image, only : image
implicit none

 public :: extract_particles, center_cc, discard_borders, pixels_dist, get_pixel_pos, polish_cc !maybe to remove center_cc from public

 private
interface pixels_dist
    module procedure pixels_dist_1
    module procedure pixels_dist_2
end interface

#include "simple_local_flags.inc"
contains


    ! This function takes in input an image and gives in output another image
    ! with smaller size in which are discarded about 4% of the borders. It
    ! is meant to deal with border effects.
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
    ! If the dist beetween the 2 particles is in [r,2r] -> delete them.
    subroutine elimin_aggregation( saved_coord, part_radius, refined_coords )
        integer,              intent(in)  :: saved_coord(:,:)     !Coordinates of picked particles
        integer,              intent(in)  :: part_radius          !Approximate radius of the particle
        integer, allocatable, intent(out) :: refined_coords(:,:)  !New coordinates of not aggregated particles
        logical, allocatable :: msk(:)
        integer              :: i, j, cnt
        real :: ppart_radius
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ppart_radius = part_radius*4. !because of shrinking
        ppart_radius = part_radius
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(msk(size(saved_coord,dim=1)), source = .true.) !initialise to 'keep all'
        cnt = 0
        do i = 1, size(saved_coord, dim = 1)             !fix one coord
            do j = 1, size(saved_coord, dim = 1)         !fix another coord to compare
                if(j > i .and. msk(i) .and. msk(j)) then !not compare twice ,and if the particles haven t been deleted yet
                    if(  sqrt(real((saved_coord(i,1)-saved_coord(j,1))**2 &
                        & + (saved_coord(i,2)-saved_coord(j,2))**2)) <= real(2*ppart_radius) &
                        & .and. real(ppart_radius) < &
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
      integer, allocatable :: xyz_saved(:,:), xyz_norep_noagg1(:,:),xyz_norep_noagg(:,:)
      integer, allocatable :: imat(:,:,:)
      real,    allocatable :: rmat(:,:,:)                     !matrix corresponding to the cc image
      integer, allocatable :: imat_masked(:,:,:)              !to identify each cc
      logical              :: outside, discard, errout, erroutB
      real :: vec_tiny(2)
      ! Initialisations
      ldim = self%get_ldim()
      box  = (part_radius)*4 + 2*part_radius !needs to be slightly bigger than the particle
      rmat = img_cc%get_rmat()
      call imgwin_particle%new([box,box,1],1.)
      !call imgwin_bin%new([box,box,1],1.)
      outside = .false.
      allocate(xyz_saved(int(maxval(rmat)),2), source = 0) ! int(maxval(rmat)) is the # of cc (likely particles)
      allocate(imat_masked(1:ldim(1),1:ldim(2),ldim(3)))
      ! Copy of the micrograph, to highlight on it the picked particles
      call self_back%copy(self)
      ! Particle identification, extraction and centering
      !open(unit = 17, file = "Statistics.txt")
      allocate(imat(ldim(1),ldim(2),ldim(3)), source = int(rmat)) !rmat is the connected component matrix
      cnt_likely = 0
      do n_cc = 1, int(maxval(rmat))   !fix cc
          imat_masked = 0              !reset
          if(any(imat == n_cc)) then
              where(imat == n_cc) imat_masked = 1
              cnt_likely = cnt_likely + 1
              pos = center_cc(imat_masked)
              xyz_saved(cnt_likely,:) = int(pos(:2))
          endif
      enddo
      deallocate(imat)
      !call elimin_aggregation(xyz_saved, 2*part_radius, xyz_norep_noagg)
      call elimin_aggregation(xyz_saved, part_radius, xyz_norep_noagg1)      ! x2 not to pick too close particles
      vec_tiny = [TINY, TINY]
      cnt_likely = 0.
      do n_cc = 1, size(xyz_norep_noagg1, dim=1)
          if(xyz_norep_noagg1(n_cc,1) > TINY) cnt_likely = cnt_likely+1
      enddo
      allocate(xyz_norep_noagg(cnt_likely,2), source = xyz_norep_noagg1(1:cnt_likely,:))
      cnt_particle = 0
      print *, 'xyz_saved = '
      call vis_mat(xyz_saved)
      print *, 'xyz_norep_noagg = '
      call vis_mat(xyz_norep_noagg)
      do n_cc = 1, size(xyz_norep_noagg, dim = 1)
          if( abs(xyz_norep_noagg(n_cc,1)) >  0) then !useless?
              if( abs(xyz_norep_noagg(n_cc,1)) >  0) call self_back%draw_picked( xyz_norep_noagg(n_cc,:),part_radius,3, 'white')
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
    ! the connected component. It corresponds to the geometric median.
    function center_cc(masked_mat) result (px)
        integer, intent(in) :: masked_mat(:,:,:)
        integer     :: px(3)               !index of the central px of the cc
        integer     :: i, j, k
        integer     :: s(3)                !shape of the input matrix
        integer     :: n_px                !counter
        integer     :: idx(2)
        real,    allocatable :: dist(:,:)        !to extract the window according to the px which minimize the dist
        logical, allocatable :: mask(:,:)        !to calc the min of an array along a specific dim
        logical, allocatable :: mask_dist(:)     !for sum dist calculation
        integer, allocatable :: pos(:,:)         !position of the pixels of a fixed cc
        s = shape(masked_mat)
        call get_pixel_pos(masked_mat, pos)
        allocate(dist(4,   size(pos, dim = 2)), source = 0.)
        allocate(mask(4,   size(pos, dim = 2)), source = .false.)
        allocate(mask_dist(size(pos, dim = 2)), source = .true.)
        mask(4,:) = .true. !to calc the min wrt the dist
        n_px = 0
        do i = 1, s(1)
            do j = 1, s(2)
                do k = 1, s(3)
                    if(masked_mat(i,j,k) > 0.5) then
                        n_px = n_px + 1
                        dist( 4, n_px) = pixels_dist([i,j,k], pos, 'sum', mask_dist)
                        dist(:3, n_px) = [real(i),real(j),real(k)]
                    endif
                enddo
            enddo
        enddo
        idx   = minloc(dist, mask)
        px(:) = int(dist(1:3, idx(2)))
        deallocate(pos, mask, dist, mask_dist)
    end function center_cc

  ! This subroutine stores in pos the indeces corresponding to
  ! the pixels with value > 0 in the binary matrix imat_masked.
  subroutine get_pixel_pos(imat_masked, pos)
      integer,              intent(in)  :: imat_masked(:,:,:)
      integer, allocatable, intent(out) :: pos(:,:)
      integer :: s(3), i, j, k, cnt
      s = shape(imat_masked)
      if(allocated(pos)) deallocate(pos)
      allocate(pos(3,count(imat_masked(:s(1),:s(2),:s(3)) > 0.5)), source = 0)
      cnt = 0
      do i = 1, s(1)
            do j = 1, s(2)
                do k = 1, s(3)
                    if(imat_masked(i,j,k) > 0.5) then
                        cnt = cnt + 1
                        pos(:3,cnt) = [i,j,k]
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
    function pixels_dist_1( px, vec, which, mask, location) result( dist )
        integer,           intent(in)     :: px(3)
        integer,           intent(in)     :: vec(:,:)
        character(len=*),  intent(in)     :: which
        logical,           intent(inout)  :: mask(:)
        integer, optional, intent(out)    :: location(1)
        real    :: dist
        integer :: i
        if(size(mask,1) .ne. size(vec, dim = 2)) THROW_HARD('Incompatible sizes mask and input vector; pixels_dist_1')
        if((any(mask .eqv. .false.)) .and. which .eq. 'sum') THROW_WARN('Not considering mask for sum; pixels_dist_1')
        !to calculation of the 'min' excluding the pixel itself, otherwise it d always be 0
        do i = 1, size(vec, dim = 2)
            if( px(1)==vec(1,i) .and. px(2)==vec(2,i) .and. px(3)==vec(3,i) )then
                if(which .ne. 'sum') mask(i) = .false.
            endif
        enddo
        select case(which)
        case('max')
            dist =  maxval(sqrt((real(px(1)-vec(1,:)))**2+(real(px(2)-vec(2,:)))**2+(real(px(3)-vec(3,:)))**2), mask)
            if(present(location)) location =  maxloc(sqrt((real(px(1)-vec(1,:)))**2+(real(px(2)-vec(2,:)))**2+(real(px(3)-vec(3,:)))**2), mask)
        case('min')
            dist =  minval(sqrt((real(px(1)-vec(1,:)))**2+(real(px(2)-vec(2,:)))**2+(real(px(3)-vec(3,:)))**2), mask)
            if(present(location)) location =  minloc(sqrt((real(px(1)-vec(1,:)))**2+(real(px(2)-vec(2,:)))**2+(real(px(3)-vec(3,:)))**2), mask)
        case('sum')
            if(present(location))  THROW_HARD('Unsupported location parameter with sum mode; pixels_dist_1')
            dist =  sum   (sqrt((real(px(1)-vec(1,:)))**2+(real(px(2)-vec(2,:)))**2+(real(px(3)-vec(3,:)))**2))
        case DEFAULT
            write(logfhandle,*) 'Pixels_dist kind: ', trim(which)
            THROW_HARD('Unsupported pixels_dist kind; pixels_dist_1')
        end select
    end function pixels_dist_1

    function pixels_dist_2( px, vec, which, mask, location) result( dist )
        real,              intent(in)     :: px(3)
        real,              intent(in)     :: vec(:,:)
        character(len=*),  intent(in)     :: which
        logical,           intent(inout)  :: mask(:)
        integer, optional, intent(out)    :: location(1)
        real    :: dist
        integer :: i
        if(size(mask,1) .ne. size(vec, dim = 2)) THROW_HARD('Incompatible sizes mask and input vector; pixels_dist_2')
        if(any(mask .eqv. .false.) .and. which .eq. 'sum') THROW_WARN('Not considering mask for sum; pixels_dist_2')
        !to calculation of the 'min' excluding the pixel itself, otherwise it d always be 0
        do i = 1, size(vec, dim = 2)
            if(      abs(px(1)-vec(1,i)) < TINY .and. abs(px(2)-vec(2,i)) < TINY  &
            &  .and. abs(px(3)-vec(3,i)) < TINY ) mask(i) = .false.
        enddo
        select case(which)
        case('max')
            dist =  maxval(sqrt((px(1)-vec(1,:))**2+(px(2)-vec(2,:))**2+(px(3)-vec(3,:))**2), mask)
            if(present(location)) location = maxloc(sqrt((px(1)-vec(1,:))**2+(px(2)-vec(2,:))**2+(px(3)-vec(3,:))**2), mask)
        case('min')
            dist =  minval(sqrt((px(1)-vec(1,:))**2+(px(2)-vec(2,:))**2+(px(3)-vec(3,:))**2), mask)
            if(present(location)) location = minloc(sqrt((px(1)-vec(1,:))**2+(px(2)-vec(2,:))**2+(px(3)-vec(3,:))**2), mask)
        case('sum')
            if(present(location))  THROW_HARD('Unsupported location parameter with sum mode; pixels_dist_1')
            dist =  sum(sqrt((px(1)-vec(1,:))**2+(px(2)-vec(2,:))**2+(px(3)-vec(3,:))**2))
        case DEFAULT
            write(logfhandle,*) 'Pixels_dist kind: ', trim(which)
            THROW_HARD('Unsupported pixels_dist kind; pixels_dist_2')
        end select
    end function pixels_dist_2

    ! This subroutine takes in input a connected components (cc)
    ! image and eliminates some of the ccs according to thei size.
    ! The decision method consists in calculate the avg size of the ccs
    ! and their standar deviation.
    ! Elimin ccs which have size: > ave + 0.8*stdev
    !                             < ave - 0.8*stdev
    subroutine polish_cc(img_cc)
        type(image), intent(inout) :: img_cc
        integer, allocatable :: sz(:)
        real :: lt, ht !low and high thresh for ccs polising
        real :: ave, stdev ! avg and stdev of the size od the ccs
        integer :: n_cc
        sz = img_cc%size_connected_comps()
        ave = sum(sz)/size(sz)
        stdev = 0.
        do n_cc = 1, size(sz)
            stdev = stdev + (sz(n_cc)-ave)**2
         enddo
        stdev = sqrt(stdev/(size(sz)-1))
        call img_cc%elim_cc([ floor(ave-0.8*stdev) , ceiling(ave+0.8*stdev) ])
        call img_cc%order_cc()
    end subroutine polish_cc
end module simple_picker_chiara
