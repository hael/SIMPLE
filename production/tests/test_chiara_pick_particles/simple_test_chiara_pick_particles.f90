module chiara_pick_particles_mod
include 'simple_lib.f08'
use simple_image
implicit none
#include "simple_local_flags.inc"

contains

    function center_edge( img_in, lp, thresh, bin_img, discard ) result( xyz )
        class(image),   intent(inout)      :: img_in
        real,           intent(in)         :: lp, thresh(1)
        logical, intent(out)               :: discard
        type(image)       :: tmp, bin_img, imgcc
        real              :: xyz(3)
        integer           :: ldim(3), min_sz
        real, allocatable :: rmat(:,:,:)
        min_sz = 30*30/2   !(part_radius**2)/2
        ldim = img_in%get_ldim()
        call tmp%copy(img_in)
        call tmp%bp(0., lp)
        call tmp%ifft()
        rmat = tmp%get_rmat()
        where( rmat < 0. ) rmat = 0.
        call tmp%sobel(bin_img, thresh(1))  !I use the threshold automatically selected for the whole micrograph
        !call tmp%automatic_thresh_sobel(bin_img)
        call bin_img%real_space_filter(3, 'median') !median filtering allows me to calculate cc in an easy way
        call bin_img%find_connected_comps(imgcc)
        call imgcc%prepare_connected_comps(discard, min_sz)
        xyz = imgcc%masscen()
        call tmp%kill
    end function center_edge

    ! This function tells whether the new_coord of a likely particle have already been identified.
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

end module chiara_pick_particles_mod

program simple_test_chiara_pick_particles
include 'simple_lib.f08'
use simple_micops
use simple_image
use simple_stackops
use simple_math
use chiara_pick_particles_mod
integer              :: ldim_shrunken(3), n_images, ifeat, box_shrunken, ldim(3), cnt, xind, yind, part_radius, cnt1
integer, allocatable :: coord_build(:,:), coord(:,:), xyz_saved(:,:), xyz_no_rep(:,:), xyz_norep_noagg(:,:)
real,    parameter   :: SHRINK = 4.
integer, parameter   :: BOX = 370, OFFSET = BOX/SHRINK-20, BOFFSET = 1
type(image)          :: img, img_extrac, imgwin, mic_original, mic_shrunken, bin_img, imgcc, &
& mic_bin, mic_back_no_aggreg
real                 :: smpd_shrunken, lp, xyz(3), sobel_thresh(1)
logical              :: outside, discard, picked
call read_micrograph( micfname = &
!& '/home/lenovoc30/Desktop/NLmeanDenoising/Denoised.mrc', smpd = 1.0)
& '/home/lenovoc30/Desktop/MassCenter/MicrographsMArion/micrographs/FoilHole_1574549_Data_1584496_1584497_20180703_1915-40767_intg.mrc'&
&, smpd = 1.)
call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
call set_box(BOX, box_shrunken)  !Here I can decide to add noise
call extract_boxes2file(OFFSET, 'extracted_windows.mrc', n_images, BOFFSET, coord)
call mic_shrunken%new(ldim_shrunken, smpd_shrunken)   !Shrunken-high pass filtered version
call mic_shrunken%read('shrunken_hpassfiltered.mrc')
call mic_shrunken%automatic_thresh_sobel( mic_bin, 2., sobel_thresh(1))  !Automatic selection of the threshold
call mic_bin%write('AutomaticThresh.mrc')
call img_extrac%new ([box_shrunken,box_shrunken,1],1.)   !it will be the extracted window
call imgwin%new ([box_shrunken,box_shrunken,1],1.)       !it will be the centered window
lp = 4.
part_radius = 30
picked = .false.
cnt  = 0
allocate(xyz_saved(n_images,2), xyz_no_rep(n_images,2), source = 0)
call mic_back_no_aggreg%copy(mic_shrunken)
do ifeat = 1, n_images
    call img_extrac%read('extracted_windows.mrc', ifeat)
    xyz = center_edge(img_extrac, lp, sobel_thresh(1), bin_img, discard )
    xyz_saved(ifeat,:) = int(xyz(:2))+coord(ifeat,:)+box_shrunken/2
    call bin_img%write('binary_images_sobel.mrc', ifeat)
    call mic_shrunken%window_slim( int(xyz(:2))+coord(ifeat,:), box_shrunken, imgwin, outside )
    if(ifeat > 1) picked = is_picked( int(xyz) + coord(ifeat,:) + box_shrunken/2, xyz_saved, part_radius, ifeat-1 )
    if( .not. picked .and. .not. discard .and. .not. outside) then
        cnt = cnt + 1
        xyz_no_rep(cnt,:) = int(xyz) + coord(ifeat,:) + box_shrunken/2
    endif
end do
call elimin_aggregation( xyz_no_rep, part_radius, xyz_norep_noagg )
do ifeat = 1, cnt
    if( abs(xyz_norep_noagg(ifeat,1)) > 0) call mic_back_no_aggreg%draw_picked( xyz_norep_noagg(ifeat,:),part_radius,2)
enddo
call mic_back_no_aggreg%write('picked_particles_no_aggreg.mrc')
cnt1 = 0
do ifeat = 1, cnt
    if( xyz_norep_noagg(ifeat,1) /=  0) then
        call mic_shrunken%window_slim( xyz_norep_noagg(ifeat,:) - box_shrunken/2, box_shrunken, imgwin, outside )
        if( .not. outside) then
            cnt1 = cnt1 + 1
            call imgwin%write('centered_particles_NOrepNOagg.mrc', cnt1)
        endif
    endif
end do
!To center again the picked_particles
! cnt = 0
! do ifeat = 1, cnt1
!     call img_extrac%read('centered_particles_NOrepNOagg.mrc', ifeat)
!     call img_extrac%mask( real(part_radius), 'soft') !This time I mask, that's why I improve
!     xyz = center_edge(img_extrac, lp, sobel_thresh(1), bin_img, discard )
!     call mic_shrunken%window_slim( int(xyz(:2))+xyz_norep_noagg(ifeat,:) - box_shrunken/2, box_shrunken, imgwin, outside )
!     if(.not. outside) then
!         cnt = cnt + 1
!         call imgwin%write('centered_particles_RECENTERED.mrc', cnt)
!     endif
! enddo
end program simple_test_chiara_pick_particles
