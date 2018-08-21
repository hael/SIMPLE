module test_chiara_pick_neg_staining
include 'simple_lib.f08'
use simple_image
implicit none
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

end module test_chiara_pick_neg_staining

program simple_test_chiara_pick_neg_staining
include 'simple_lib.f08'
use simple_micops
use simple_image
use simple_stackops
use simple_math
use test_chiara_pick_neg_staining
type(image)          :: mic_shrunken, mic_bin, mic_back, mic_back_no_aggreg, mic_copy
type(image)          :: img_extrac, imgwin
integer              :: ldim_shrunken(3), box_shrunken
integer              :: n_images, ifeat, cnt, cnt1, xind, yind
integer              :: part_radius, median_filt_winsz
integer, allocatable :: coord(:,:)
integer, allocatable :: xyz_saved(:,:), xyz_no_rep(:,:), xyz_norep_noagg(:,:)
real,    parameter   :: SHRINK = 4.
integer, parameter   :: BOX = 230, OFFSET = BOX/SHRINK-20, BOFFSET = 20
real                 :: smpd_shrunken, lp, xyz(3), sobel_thresh(1)
logical              :: outside, discard, picked
! Parameters selection
lp = 30.                 !low pass filtering threshold
part_radius = 15         !extimate particle radius
median_filt_winsz = 6    !size of the window for median filter
call read_micrograph( micfname ='/home/lenovoc30/Desktop/MassCenter/NegativeStaining/16.06.34 CCD Acquire_0000_1.mrc'&
&, smpd = 1.)
! Image processing on the micrograph
! 1) shrinking and hp filtering
call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
call set_box(BOX, box_shrunken)
call mic_shrunken%new(ldim_shrunken, smpd_shrunken)   !Shrunken-high pass filtered version
call mic_shrunken%read('shrunken_hpassfiltered.mrc')
call mic_shrunken%negative_image()                    !neg staining --> similar to cryo
call mic_shrunken%write('shrunken_hpassfiltered.mrc')
! Copying the shrunken hp filtered mic, to work on it without modifying the original one
! 2) low pass filtering
call mic_copy%copy(mic_shrunken)
call mic_copy%bp(0.,lp)
call mic_copy%ifft()
! 3) Edge detection
call mic_copy%automatic_thresh_sobel( mic_bin, 0.3, sobel_thresh(1)) !Automatic threshold selection for the entire micrograph
! 4) Median filtering
call mic_bin%real_space_filter(median_filt_winsz, 'median')          !Median filtering allows easly calculation of cc
call mic_bin%write('MedianFiltered.mrc')
! Windows extraction
call extract_windows2file(mic_bin, OFFSET, 'extracted_windows.mrc', n_images, BOFFSET, coord)
call img_extrac%new ([box_shrunken,box_shrunken,1],1.)   !it will be the extracted window
call imgwin%new     ([box_shrunken,box_shrunken,1],1.)   !it will be the centered  window
picked = .false. !logical, not to pick the same particle more than once
cnt  = 0         !counter
allocate(xyz_saved(n_images,2), xyz_no_rep(n_images,2), source = 0)  !vectors of coordinates
! Copy of the lp-micrograph, to highlight on it the picked particles
call mic_back_no_aggreg%copy(mic_shrunken)
! Process of idenfifying particles and centering them
do ifeat = 1, n_images
call img_extrac%read('extracted_windows.mrc', ifeat)
xyz = img_extrac%center_edge( discard )  !center coordinates (found by mean of mass center)
xyz_saved(ifeat,:) = int(xyz(:2))+coord(ifeat,:)+box_shrunken/2 !saving coordinates, not to pick more than once
call mic_shrunken%window_slim   ( int(xyz(:2)) + coord(ifeat,:),  box_shrunken, imgwin, outside )
if(ifeat > 1) picked = is_picked( int(xyz(:2)) + coord(ifeat,:) + box_shrunken/2, xyz_saved, part_radius, ifeat-1 )
if( .not. picked .and. .not. discard .and. .not. outside) then
    cnt = cnt + 1
    call imgwin%write('centered_particles.mrc', cnt)
    xyz_no_rep(cnt,:) = int(xyz) + coord(ifeat,:) + box_shrunken/2 !coordinates of the particles with no repetition
endif
end do
call elimin_aggregation( xyz_no_rep, part_radius, xyz_norep_noagg ) !xyz_norep_noagg = coordinates, no aggregations
! Extract the particles with no aggregations from the mic
cnt1 = 0
do ifeat = 1, cnt
if( abs(xyz_norep_noagg(ifeat,1)) > 0) then
    if( abs(xyz_norep_noagg(ifeat,1)) >  0) call mic_back_no_aggreg%draw_picked( xyz_norep_noagg(ifeat,:), part_radius,2)
    call mic_shrunken%window_slim( xyz_norep_noagg(ifeat,:) - box_shrunken/2, box_shrunken, imgwin, outside )
    if( .not. outside) then
        cnt1 = cnt1 + 1
        call imgwin%write('centered_particles_NOagg.mrc', cnt1)
    endif
endif
end do
call mic_back_no_aggreg%write('picked_particles_no_aggreg.mrc')
!deallocate(coord, xyz_saved, xyz_no_rep, xyz_norep_noagg)
end program simple_test_chiara_pick_neg_staining

! To center again the picked_particles
! cnt = 0
! do ifeat = 1, cnt1
!     call img_extrac%read('centered_particles_NOrepNOagg.mrc', ifeat)
!     call img_extrac%mask( real(part_radius), 'soft') !This time I mask, that's why I improve
!     xyz = center_edge(img_extrac, discard )
!     !PRINT *, "ifeat = ", ifeat, 'xyz = ', ceiling(xyz(:2))
!     call mic_shrunken%window_slim( ceiling(xyz(:2))+xyz_norep_noagg(ifeat,:) - box_shrunken/2, box_shrunken, imgwin, outside )
!     if(.not. outside) then
!         cnt = cnt + 1
!         call imgwin%write('centered_particles_RECENTERED.mrc', cnt)
!         call mic_adjusted%draw_picked(int(xyz(:2))+xyz_norep_noagg(ifeat,:),3000., 4)
!     endif
! enddo
! call mic_adjusted%write('picked_particles_REcentered.mrc.mrc')
