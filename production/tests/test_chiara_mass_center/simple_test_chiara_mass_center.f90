module chiara_mass_center_mod
include 'simple_lib.f08'
use simple_image, only: image
implicit none
#include "simple_local_flags.inc"

contains

    ! This function is to enumerate the labels of a connected components image.
    ! This routine assumes that the input pixel values are 0 or cnt (see enumerate_white_pixels, above)
    ! subroutine enumerate_connected_comps( img )
    !     class(image), intent(inout) :: img
    !     integer           :: i, j
    !     real, allocatable :: tmp_mat(:,:,:) !unnecessary if in simple_image module
    !     real :: tmp
    !     allocate(tmp_mat(img%ldim(1), img%ldim(2), 1), source = 0.)
    !     do i = 1, img%ldim(1)
    !         do j = 1, img%ldim(2)
    !             if( img%rmat(i,j,1) > 0.5 ) then  !rmat == 0  --> background
    !                 tmp = img%rmat(i,j,1) ! tmp is now [1,cnt]
    !                 where(img%rmat == tmp)
    !                     img%rmat = 0.     !Not to consider this cc again
    !                     tmp_mat = tmp
    !                 endwhere
    !             endif
    !         enddo
    !     enddo
    ! end subroutine enumerate_connected_comps

    !The output of this function is the crossing number of a pixel as defined in
    !'Binary digital image processing' (book).
    !The crossing number indicates the number of 4-connected components in the 8-neighbourhoods.
    !It takes in input the 8-neighbourhoods of a px in a BINARY image.
    !It doesn't check if the image is binarised, but it should be for definition.
    function crossing_number(neigh_8) result(cn)
        real, intent(in) :: neigh_8(:)
        real    :: cn         !crossing number of the pixel px
        real    :: prod, sum  !just for comfort
        integer :: i          !loop
        prod = 1.
        sum  = 0.
        do i = 1, size(neigh_8)-1  !-1 because I don't consider the pixel itself (stored in the last entrance)
            prod  = prod*neigh_8(i)
            if(i<size(neigh_8)-1) then
                sum = sum + abs(neigh_8(i+1)-neigh_8(i))
            endif
        enddo
        sum = sum + abs(neigh_8(1) - neigh_8(size(neigh_8)-1)) !first one with the last one
        sum = 0.5*sum
        cn = prod + sum
    end function crossing_number

    !This subroutine eliminates isolated white points from a binary image
    subroutine elimin_isolated_points(bin_img)
        type(image), intent(inout) :: bin_img
        real, allocatable :: neigh_8(:), rmat(:,:,:)
        integer           :: ldim(3), i, j
        real              :: cn   !crossing number
        ldim = bin_img%get_ldim()
        rmat = bin_img%get_rmat()
        do i = 1, ldim(1)
            do j = 1, ldim(2)
                if(rmat(i,j,1) /= 0.) then         !Consider just white pixels
                    neigh_8 = bin_img%calc_neigh_8([i,j,1])
                    cn = crossing_number(neigh_8)
                    if(cn == 0.) rmat(i,j,1) = 0.  !Eliminate isolated points
                    deallocate(neigh_8)
                endif
            enddo
        enddo
        call bin_img%set_rmat(rmat)
        deallocate(rmat)
    end subroutine elimin_isolated_points

    function center_edge( img_in, lp, thresh, bin_img, discard ) result( xyz )
        class(image),   intent(inout)      :: img_in
        real,           intent(in)         :: lp, thresh(1)
        logical, intent(out)               :: discard
        type(image)       :: tmp, bin_img, imgcc
        real              :: xyz(3)
        integer           :: ldim(3), min_sz
        real, allocatable :: rmat(:,:,:)
        min_sz = 250
        ldim = img_in%get_ldim()
        call tmp%copy(img_in)
        call tmp%bp(0., lp)
        call tmp%ifft()
        rmat = tmp%get_rmat()
        where( rmat < 0. ) rmat = 0.
        call tmp%sobel(bin_img, thresh(1))
        call bin_img%real_space_filter(3, 'median') !median filtering allows me to calculate cc in an easy way
        call bin_img%find_connected_comps(imgcc)
        call imgcc%prepare_connected_comps(discard, min_sz)
        xyz = imgcc%masscen()
        call tmp%kill
    end function center_edge

    subroutine print_mat(matrix)
        integer , intent(in) :: matrix(:,:)
        integer :: j, s(2)
        s = shape(matrix)
        do j = 1, s(1)
            print *, matrix(j,:)
        enddo
    end subroutine print_mat

    function is_picked( new_coord,saved_coord, part_radius, saved ) result(yes_no)
        integer, intent(in) :: new_coord(2)       !Coordinates of a new window to extract
        integer, intent(in) :: saved_coord(:,:)   !Coordinates of extracted windows
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

end module chiara_mass_center_mod

program simple_test_chiara_mass_center
include 'simple_lib.f08'
use simple_micops
use simple_image
use simple_stackops
use simple_math
use chiara_mass_center_mod
integer              :: ldim_shrunken(3), n_images, ifeat, box_shrunken, ldim(3), cnt, xind, yind, part_radius
integer, allocatable :: coord_build(:,:), coord(:,:), coord_center_edge(:,:)
real,    parameter   :: SHRINK = 2.
integer, parameter   :: BOX = 256, OFFSET = BOX/SHRINK-20, BOFFSET = 1, N_PROJ = 64
type(image)          :: img, img_extrac, imgwin, mic_original, mic_shrunken, bin_img
real                 :: smpd_shrunken, lp, xyz(3), thresh(1)
logical              :: outside, discard, picked
logical, allocatable :: mask(:)
cnt = 0
do xind=1,BOX*8,BOX
    do yind=1,BOX*8,BOX
        cnt = cnt + 1
    end do
end do
allocate(coord_build(cnt,2))
cnt = 0
do xind=1,BOX*8,BOX
    do yind=1,BOX*8,BOX
        cnt = cnt + 1
        coord_build(cnt,:) = [xind,yind]
    end do
end do
call mic_original%new([BOX*8, BOX*8,1], 1.) !building my micrograph
call img%new([BOX,BOX,1],1.)
do ifeat=1,N_PROJ
    call img%read('outstk.mrc', ifeat)
    call mic_original%add_window(img, coord_build(ifeat,:))
end do
deallocate(coord_build)
call mic_original%write('original_micrograph.mrc')
call read_micrograph( micfname = 'original_micrograph.mrc', smpd = 1.0)
call shrink_micrograph(SHRINK, ldim_shrunken, smpd_shrunken)
call set_box(BOX, box_shrunken,0.2)  !Here I can decide to add noise
part_radius = 20
call extract_boxes2file(OFFSET, 'extracted_windows.mrc', n_images, BOFFSET, coord)
call mic_shrunken%new(ldim_shrunken, smpd_shrunken)   !Shrunken-high pass filtered version
call mic_shrunken%read('shrunken_hpassfiltered.mrc')
call img_extrac%new ([box_shrunken,box_shrunken,1],1.)   !it will be the extracted window
call imgwin%new ([box_shrunken,box_shrunken,1],1.)       !it will be the centered window
lp = 4.
cnt  = 0
allocate(coord_center_edge(n_images,2), source = 0)
picked = .false. !initialization
thresh(1) = 0.99
do ifeat = 1, n_images
    call img_extrac%read('extracted_windows.mrc', ifeat)
    xyz = center_edge(img_extrac, lp, thresh(1), bin_img, discard )
    call bin_img%write('binary_images_sobel.mrc', ifeat)
    if(.not. discard) then    !it means sobel has found enough points
        call mic_shrunken%window_slim( [int(xyz(1))+coord(ifeat,1), int(xyz(2))+coord(ifeat,2)], &
                                                          &  box_shrunken, imgwin, outside )
        if(ifeat > 1) picked = is_picked([int(xyz(1))+coord(ifeat,1), int(xyz(2))+coord(ifeat,2)],coord_center_edge, part_radius, ifeat-1 )
        coord_center_edge(ifeat,:) =     [int(xyz(1))+coord(ifeat,1), int(xyz(2))+coord(ifeat,2)]
        if(.not. picked) then  !window not already selected
            if( .not. outside  )then  !not outside of the borders of the image
               cnt = cnt + 1
               call imgwin%write('centered_particles_edgeNOREP.mrc', cnt)
            endif
        endif
    endif
end do
end program simple_test_chiara_mass_center
