module connected_components
include 'simple_lib.f08'
use simple_image, only : image
implicit none
#include "simple_local_flags.inc"

contains

    subroutine print_mat(matrix)
        integer, intent(in) :: matrix(:,:,:)
        integer             :: j, s(3)
        s = shape(matrix)
        do j = 1, s(1)
            print *, matrix(j,:,1)
        enddo
    end subroutine print_mat

    subroutine build_ellipse(img_in,center, axis, rot)
        type(image), intent(inout) :: img_in
        real,        intent(in)    :: center(2), axis(2), rot
        real, allocatable          :: theta(:)
        integer                    :: ldim(3), i, j, k
        real :: rot_r
        if(rot < 0. .or. rot > 360. ) THROW_HARD("please insert an angle in the range [0,360]")
        rot_r = deg2rad(rot)
        ldim = img_in%get_ldim()
        if(ldim(3) /= 1) THROW_HARD("the image has to be 2D!")
        theta = (/ (deg2rad(real(i)),i=1,360) /)
        do k = 1,size(theta)
            do i = 1,ldim(1)
                do j = 1,ldim(2)
                    if(abs(real(i)-center(1)-axis(1)*cos(theta(k))*cos(rot_r)+axis(2)*sin(theta(k))*sin(rot_r))<1 .and. &
                    &  abs(real(j)-center(1)-axis(1)*cos(theta(k))*sin(rot_r)-axis(2)*sin(theta(k))*cos(rot_r))<1) then
                        call img_in%set([i,j,1], 1.)
                    end if
                enddo
            enddo
        enddo
        deallocate(theta)
    end subroutine build_ellipse

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
        do i = 1, size(saved_coord, dim = 1)            !fix one coord
            do j = 1, size(saved_coord, dim = 1)        !fix another coord to compare
                if(j > i .and. msk(i) .and. msk(j)) then !not to compare twice the same couple and if the selected particles
                                                        !haven t been deleted yet
                    if(  sqrt(real((saved_coord(i,1)-saved_coord(j,1))**2 + (saved_coord(i,2)-saved_coord(j,2))**2)) <= real(2*part_radius)) then   !if the distance beetween the 2 selected particles is <= 2r, delete them
                        msk(i) = .false.
                        msk(j) = .false.
                        cnt = cnt + 1           !number of deleted couples
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

    subroutine extract_particles(self, img_cc, part_radius)
        type(image),  intent(inout) :: self        !original image
        type(image),  intent(inout) :: img_cc      !connected components image
        integer,      intent(in)    :: part_radius !extimated particle radius, to discard aggregations
        type(image)          :: imgwin_particle, imgwin_centered, self_back_no_aggreg
        integer              :: box
        integer              :: n_cc               ! n_cc  = # cc in the input image
        integer              :: cnt_particle, cnt_centered  !counters
        integer              :: xyz(3)            !mass center coordinates
        integer, allocatable :: pos(:)            !position of each cc
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
                    xyz = imgwin_particle%masscen() !particle centering
                    xyz_saved(cnt_particle,:) =      int(xyz(:2)) + pos(:2)
                    if(n_cc > 1) picked = is_picked( int(xyz(:2)) + pos(:2), xyz_saved, part_radius, cnt_particle-1 )
                    call self%window_slim(int(xyz(:2))+pos(:2)-box/2, box, imgwin_centered, outside_centered)
                    if( .not. picked .and. .not. discard .and. .not. outside_centered) then
                        cnt_centered = cnt_centered + 1
                        call imgwin_centered%write('CenteredParticleEdge.mrc', cnt_centered)
                        xyz_no_rep(cnt_centered,:) =int(xyz(:2))+pos(:2)
                    endif
                endif
            endif
        enddo
        call elimin_aggregation( xyz_no_rep, part_radius+5, xyz_norep_noagg ) !+5 not to pick too close particles
        cnt_particle = 0
        do n_cc = 1, cnt_centered
            if( abs(xyz_norep_noagg(n_cc,1)) >  0) then
                if( abs(xyz_norep_noagg(n_cc,1)) >  0) call self_back_no_aggreg%draw_picked( xyz_norep_noagg(n_cc,:), part_radius,2)
                call self%window_slim( xyz_norep_noagg(n_cc,:)-box/2, box, imgwin_centered, outside_centered )
                if( .not. outside_centered) then
                    cnt_particle = cnt_particle + 1
                    call imgwin_centered%write('centered_particles_NOrepNOagg.mrc', cnt_particle)
                endif
            endif
        end do
        call self_back_no_aggreg%write('picked_particles_no_aggreg.mrc')
        deallocate(xyz_saved, xyz_no_rep, xyz_norep_noagg, rmat)
    end subroutine extract_particles

  ! !!!!!!!!!!!!!ADDED BY CHIARA!!!!!!!!!!!!!!
  ! ! This functions draws a white square in self centered in square_center.
  ! ! wide is the length of the side of the square.
  ! ! This function is meant to mark picked particles (visual purpose).
  ! subroutine draw_square(self, square_center, value, wide)
  !   class(image),      intent(inout) :: self
  !   integer,           intent(in)    :: square_center(2)  !coordinates of the center of the square
  !   real,              intent(in)    :: value             !intensity value of the square
  !   integer, optional, intent(in)    :: wide              !side width of the square
  !   integer :: wwide
  !   if( .not. self%is_2d() ) THROW_HARD('only for 2D images; draw_square')
  !   wwide = 7
  !   if(present(wide)) wwide = wide
  !   if(square_center(1)-wwide < 1 .or. square_center(1)+wwide > self%ldim(1) .or. &
  !   &  square_center(2)-wwide < 1 .or. square_center(2)+wwide > self%ldim(2) ) then
  !     print *, 'The square is put of the border of the image!'
  !     return  !Do not throw error, just do not draw
  !   endif
  !   self%rmat(square_center(1)-wwide:square_center(1)+wwide, square_center(2)-wwide:square_center(2)+wwide, 1) = value
  ! end subroutine draw_square


  ! subroutine draw_picked(self, part_coords, part_radius, border)
  !   type(image),intent(inout)     :: self
  !   integer,    intent(in)        :: part_coords(2)  !coordinates of the picked particle
  !   integer,    intent(in)        :: part_radius
  !   integer, optional, intent(in) :: border          !width of the border of the drawn line
  !   integer :: bborder
  !   real    :: value             !intensity value of the window
  !   integer :: wide              !side width of the window
  !   integer :: length            !length of the drawn side
  !
  !   real, allocatable :: rmat(:,:,:)
  !   integer :: ldim(3)
  !   ldim = self%get_ldim()
  !   rmat = self%get_rmat()
  !   bborder = 1
  !   if(present(border)) bborder = border
  !   value = maxval(rmat(:,:,:))+100
  !   wide = 4*part_radius
  !   length = int(part_radius/2)
  !   if(part_coords(1)-wide/2-int((bborder-1)/2) < 1 .or. part_coords(1)+wide/2+int((bborder-1)/2) > ldim(1) .or. &
  !   &  part_coords(2)-wide/2-int((bborder-1)/2) < 1 .or. part_coords(2)+wide/2+int((bborder-1)/2) > ldim(2) ) then
  !     return  !Do not throw error, just do not draw
  !   endif
  !   ! Edges of the window
  !   rmat(part_coords(1)-wide/2:part_coords(1)-wide/2+length, &
  !      & part_coords(2)-wide/2-int((bborder-1)/2):part_coords(2)-wide/2+int((bborder-1)/2), 1) = value
  !   rmat(part_coords(1)-wide/2-int((bborder-1)/2):part_coords(1)-wide/2+int((bborder-1)/2),&
  !      & part_coords(2)-wide/2:part_coords(2)-wide/2+length, 1) = value
  !   rmat(part_coords(1)+wide/2-length:part_coords(1)+wide/2,&
  !      & part_coords(2)-wide/2-int((bborder-1)/2):part_coords(2)-wide/2+int((bborder-1)/2), 1) = value
  !   rmat(part_coords(1)+wide/2-int((bborder-1)/2):part_coords(1)+wide/2+int((bborder-1)/2),&
  !      & part_coords(2)-wide/2:part_coords(2)-wide/2+length, 1) = value
  !   rmat(part_coords(1)-wide/2:part_coords(1)-wide/2+length,&
  !      & part_coords(2)+wide/2-int((bborder-1)/2):part_coords(2)+wide/2+int((bborder-1)/2), 1) = value
  !   rmat(part_coords(1)-wide/2-int((bborder-1)/2):part_coords(1)-wide/2+int((bborder-1)/2),&
  !      & part_coords(2)+wide/2-length:part_coords(2)+wide/2, 1) = value
  !   rmat(part_coords(1)+wide/2-length:part_coords(1)+wide/2,&
  !      & part_coords(2)+wide/2-int((bborder-1)/2):part_coords(2)+wide/2+int((bborder-1)/2), 1) = value
  !   rmat(part_coords(1)+wide/2-int((bborder-1)/2):part_coords(1)+wide/2+int((bborder-1)/2),&
  !   & part_coords(2)+wide/2-length:part_coords(2)+wide/2, 1) = value
  !
  !   ! Central cross
  !   rmat(part_coords(1)-length:part_coords(1)+length,  part_coords(2), 1) = value
  !   rmat(part_coords(1),  part_coords(2)-length:part_coords(2)+length, 1) = value
  !   call self%set_rmat(rmat)
  ! end subroutine draw_picked

end module connected_components

program simple_test_chiara_try
include 'simple_lib.f08'
use simple_micops
use simple_image
use simple_stackops
use simple_math
use connected_components
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
type(cmdline)        :: cline
type(parameters)     :: params
type(image)          :: mic_shrunken,  mic_bin, mic_back_no_aggreg, mic_copy
type(image)          :: imgcc, imgwin
integer              :: ldim_shrunken(3), box_shrunken,  xind, yind, min_sz, max_sz
real                 :: part_radius
real,    parameter   :: SHRINK = 4.
integer, parameter   :: BOX = 230, OFFSET = BOX/SHRINK-20, BOFFSET = 1
real                 :: smpd_shrunken, lp, sobel_thresh(1)
logical              :: outside, discard, picked
call mic_shrunken%new([1024,1024,1],1.)
call mic_shrunken%read('/home/lenovoc30/Desktop/MassCenter/try/shrunken_hpassfiltered.mrc')
part_radius = 15.
call mic_shrunken%draw_picked([100,100], int(part_radius), 2)
call mic_shrunken%write('DrownMic.mrc')
end program simple_test_chiara_try
! matrix = reshape(real([ 1,1,1,0,0,6,5, &
!                  & 1,1,0,0,6,6,6, &
!                  & 1,0,0,2,0,6,0, &
!                  & 0,0,2,2,0,0,4, &
!                  & 0,5,0,0,0,4,4, &
!                  & 0,5,5,5,0,0,0, &
!                  & 0,5,5,0,0,3,3]),[7,7,1])
! ldim = [box,box,1]
! call img_in%new(ldim,1.)
! call img_in%disc( ldim, 1.0 , 10.)
! call build_ellipse(img_in,[10.,10.], [4.,4.], 0.)
! call build_ellipse(img_in,[20.,20.], [4.,4.], 0.)
! call build_ellipse(img_in,[90.,10.], [4.,4.], 0.)
! call build_ellipse(img_in,[30.,30.], [4.,4.], 0.)
! call build_ellipse(img_in,[40.,40.], [4.,4.], 0.)
! call build_ellipse(img_in,[50.,50.], [4.,4.], 0.)
! call build_ellipse(img_in,[60.,60.], [4.,4.], 0.)
! call build_ellipse(img_in,[70.,70.], [4.,4.], 0.)
! call build_ellipse(img_in,[80.,80.], [4.,4.], 0.)
! call build_ellipse(img_in,[90.,90.], [4.,4.], 0.)
! call build_ellipse(img_in,[100.,100.], [4.,4.], 0.)
! call build_ellipse(img_in,[110.,110.], [4.,4.], 0.)
! call build_ellipse(img_in,[120.,120.], [4.,4.], 0.)
