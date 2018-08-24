module simple_picker_chiara
include 'simple_lib.f08'
use simple_image, only : image
implicit none

public :: extract_particles

private
#include "simple_local_flags.inc"
contains

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
      if(s(2) /= 2) THROW_HARD('dimension error')
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
              if(j > i .and. msk(i) .and. msk(j)) then !not to compare twice the same couple and if the selected
                                                       !particles haven t been deleted yet
                  !if the dist beetween the 2 particles is <= 2r, delete them
                  if(  sqrt(real((saved_coord(i,1)-saved_coord(j,1))**2 &
                             & + (saved_coord(i,2)-saved_coord(j,2))**2)) <= real(2*part_radius)) then
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
  ! and extract particles using a mass centering process.
  !notation:: cc = connected component.
  subroutine extract_particles(self, img_cc, part_radius)
      type(image),  intent(inout) :: self        !original image
      type(image),  intent(inout) :: img_cc      !connected components image
      integer,      intent(in)    :: part_radius !extimated particle radius, to discard aggregations
      type(image)          :: imgwin_particle, imgwin_centered, self_back_no_aggreg
      integer              :: box
      integer              :: n_cc               ! n_cc  = # cc in the input image
      integer              :: cnt_particle, cnt_centered  !counters
      real                 :: xyz_real(3)
      integer              :: xyz(3)             ! mass center coordinates
      integer, allocatable :: pos(:)             ! position of each cc
      integer, allocatable :: xyz_saved(:,:), xyz_no_rep(:,:), xyz_norep_noagg(:,:)
      real,    allocatable :: rmat(:,:,:)
      logical              :: outside_particle, outside_centered, discard, picked
      ! Initialisations
      box  = (part_radius)*4
      rmat = img_cc%get_rmat()
      cnt_particle = 0
      cnt_centered = 0
      call imgwin_particle%new([box,box,1],1.)
      call imgwin_centered%new([box,box,1],1.)
      outside_particle  = .false.
      outside_centered  = .false.
      picked = .false.
      allocate(xyz_saved(int(maxval(rmat)),2),xyz_no_rep(int(maxval(rmat)),2), source = 0) ! int(maxval(rmat)) is the # of cc (likely particles)
      ! Some copies of the lp-micrograph, to highlight on them the picked particles
      call self_back_no_aggreg%copy(self)
      ! Particle identification, extraction and centering
      do n_cc = 1, int(maxval(rmat))   !automatically do not consider background
          if(minval(abs(rmat-n_cc)) < TINY) then
              pos = minloc(abs(rmat-n_cc)) !where is the selected cc, particle identification
              call img_cc%window_slim(pos-box/2, box, imgwin_particle, outside_particle) !particle extraction
              if(.not. outside_particle) then
                  cnt_particle = cnt_particle + 1
                  call  imgwin_particle%prepare_connected_comps(discard, 1)
                  call  imgwin_particle%morpho_closing()
                  call imgwin_particle%masscen(xyz_real) !particle centering
                  xyz = nint(xyz_real)
                  xyz_saved(cnt_particle,:) =      int(xyz(:2)) + pos(:2)
                  if(n_cc > 1) picked = is_picked( int(xyz(:2)) + pos(:2), xyz_saved, part_radius, cnt_particle-1 )
                  call self%window_slim(int(xyz(:2))+pos(:2)-box/2, box, imgwin_centered, outside_centered)
                  if( .not. picked .and. .not. discard .and. .not. outside_centered) then
                      cnt_centered = cnt_centered + 1
                      xyz_no_rep(cnt_centered,:) =int(xyz(:2))+pos(:2)
                  endif
              endif
          endif
      enddo
      call elimin_aggregation( xyz_no_rep, part_radius+5, xyz_norep_noagg ) !+5 not to pick too close particles
      cnt_particle = 0
      do n_cc = 1, cnt_centered
          if( abs(xyz_norep_noagg(n_cc,1)) >  0) then
              if( abs(xyz_norep_noagg(n_cc,1)) >  0) call self_back_no_aggreg%draw_picked( xyz_norep_noagg(n_cc,:),part_radius,3)
              call self%window_slim( xyz_norep_noagg(n_cc,:)-box/2, box, imgwin_centered, outside_centered )
              if( .not. outside_centered) then
                  cnt_particle = cnt_particle + 1
                  call imgwin_centered%write('centered_particles.mrc', cnt_particle)
              endif
          endif
      end do
      call self_back_no_aggreg%write('picked_particles.mrc')
      deallocate(xyz_saved, xyz_no_rep, xyz_norep_noagg, rmat)
  end subroutine extract_particles
end module simple_picker_chiara
