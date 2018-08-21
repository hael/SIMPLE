module test_chiara_pick_edge
include 'simple_lib.f08'
use simple_image
implicit none
#include "simple_local_flags.inc"

contains

    ! This function tells whether the new_coord of a likely particle have already been identified.
    function is_picked( new_coord, saved_coord, part_radius, saved ) result(yes_no)
        integer, intent(in) :: new_coord(2)       !Coordinates of a new particle to pick
        integer, intent(in) :: saved_coord(:,:)   !Coordinates of picked particles
        integer, intent(in) :: part_radius        !Approximate radius of the particle
        integer, intent(in), optional :: saved    !How many particles have already been saved
        logical :: yes_no
        integer :: iwind, ssaved, s(2)
        s = shape(saved_coord)
        if(s(2) /= 2) then
            THROW_HARD('dimension error; is_picked')
        endif
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
    ! It takes in input the list of  ! where( rmat < 0. )
    !    rmat = 0.
    ! end where the coordinates of the identified
    ! particles and the approximate radius of the particle.
    ! The output is a new list of coordinates where aggregate picked_particles
    ! have been deleted
    subroutine elimin_aggregation( saved_coord, part_radius, refined_coords )
        integer,              intent(in)  :: saved_coord(:,:)     !Coordinates of picked particles
        integer,              intent(in)  :: part_radius          !Approximate radius of the particle
        integer, allocatable, intent(out) :: refined_coords(:,:)  !New coordinates of not aggregated particles
        logical, allocatable :: msk(:)
        integer              :: i, j, cnt
        allocate(msk(size(saved_coord, dim = 1)), source = .true.)
        cnt = 0
        do i = 1, size(saved_coord, dim = 1)        !fix one coord
            do j = 1, size(saved_coord, dim = 1)    !fix another coord to compare
               if(j > i .and. msk(i) .and. msk(j)) then !not to compare twice the same couple and if it hasn t already been deleted
                    if(  sqrt(real((saved_coord(i,1)-saved_coord(j,1))**2 + (saved_coord(i,2)-saved_coord(j,2))**2)) <= real(2*part_radius)) then
                        msk(i) = .false.
                        msk(j) = .false.
                        cnt = cnt + 1
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
                if( abs(xyz_norep_noagg(n_cc,1)) >  0) call self_back_no_aggreg%draw_picked( xyz_norep_noagg(n_cc,:),part_radius,3)
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

end module test_chiara_pick_edge

program simple_test_chiara_pick_neg_staining
include 'simple_lib.f08'
use simple_micops
use simple_image
use simple_stackops
use simple_math
use test_chiara_pick_edge
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
type(cmdline)      :: cline
type(parameters)   :: params
type(image)        :: mic_shrunken,  mic_bin, mic_back_no_aggreg, mic_copy
type(image)        :: imgcc, imgwin
integer            :: ldim_shrunken(3), box_shrunken,  xind, yind, min_sz, max_sz
real               :: part_radius
real,    parameter :: SHRINK = 4.
integer, parameter :: BOX = 230, OFFSET = BOX/SHRINK-20, BOFFSET = 1
real               :: smpd_shrunken, lp, sobel_thresh(1)
logical            :: outside, discard, picked
! In this new approach I don't extract any window, I just work directly on the entire
! micrograph using edge detection.

! Image processing steps: 1) Shrink and high pass filtering
!                         2) Low pass filtering
!                         3) Edge Detection
!                         4) Median Filtering
!                         5) Connected components (cc) identification
!                         6) cc filtering

! As a test use
!      fname = '/home/lenovoc30/Desktop/MassCenter/NegativeStaining/16.06.34 CCD Acquire_0000_1.mrc'
!      smpd  = 1.
!      part_radius = 15.
if( command_argument_count() < 3 )then
    write(*,'(a)',advance='no') 'simple_test_chiara_try smpd=<sampling distance(in A)> [fname = file name] [part_radius = <radius of the particle (# pixels)]'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('fname', 1)
call cline%checkvar('smpd', 2)
call cline%checkvar('part_radius', 3)
call params%new(cline)              !<read cline parameters
call read_micrograph(micfname = params%fname, smpd = params%smpd)
! 1) Shrink and high pass filtering
call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
call set_box(BOX, box_shrunken)
call mic_shrunken%new(ldim_shrunken, smpd_shrunken)
! 1b) Negative image, visual purpose to make it look like if it was cryo em
call mic_shrunken%read('shrunken_hpassfiltered.mrc')
call mic_shrunken%negative_image()
call mic_shrunken%write('shrunken_hpassfiltered.mrc')
! 2) Low pass filtering
lp = 30.
call mic_copy%copy(mic_shrunken) !work on a copy not to modify the original mic
call mic_copy%bp(0.,lp)
call mic_copy%ifft()
call mic_copy%write('LowPassFiltered.mrc')
! 3) Edge Detection
call mic_copy%automatic_thresh_sobel( mic_bin, 0.3, sobel_thresh(1))  !Automatic threshold selection
call mic_bin%write('AutomaticThresh.mrc')
! 4) Median Filtering
call mic_bin%real_space_filter(6, 'median') !median filtering allows me to calculate cc in an easy way
call mic_bin%write('MedianFiltered.mrc')
! 5) Connected components (cc) identification
call imgcc%new(ldim_shrunken, smpd_shrunken)
call mic_bin%find_connected_comps(imgcc)
! 6) cc filtering
part_radius = params%part_radius
min_sz =  7*int(part_radius)
max_sz = 50*int(part_radius)
call imgcc%elim_cc([min_sz,max_sz])
call imgcc%write('ElimSmallBigCC.mrc')
! Particle extraction
call extract_particles(mic_shrunken, imgcc, int(part_radius))
end program simple_test_chiara_pick_neg_staining
