module simple_nano_picker_utils
include 'simple_lib.f08'
use simple_image
use simple_atoms
use simple_parameters
use simple_srch_sort_loc, only : hpsort
implicit none
#include "simple_local_flags.inc"
    
    contains 
    
    subroutine find_closest( coords_1, coords_2, length_1, length_2, distances, filename )
        integer,                    intent(in)  :: length_1, length_2
        real,                       intent(in)  :: coords_1(3,length_1), coords_2(3,length_2)
        real,                       intent(out) :: distances(length_1)
        character(len=*), optional, intent(in)  :: filename
        integer :: i, j, min_loc(1)
        real :: closest_coord(3)
        real, allocatable :: dists(:)
        real :: min_dist
    
        if( present(filename) )then
            open(unit=22, file=filename)
        else
            open(unit=22, file='combined_coords.csv')
        endif
    
        do i = 1, length_1
            allocate(dists(length_2))
            do j = 1, length_2
                dists(j) = euclid(real(coords_1(:,i)),real(coords_2(:,j)))
            enddo
            min_loc       = minloc(dists)
            min_dist      = minval(dists)
            distances(i)  = min_dist
            closest_coord = coords_2(:,min_loc(1))
            write(22,'(7(f8.3,a))') coords_1(1,i), ',', coords_1(2,i), ',', coords_1(3,i), ',', closest_coord(1), ',', closest_coord(2), ',', closest_coord(3), ',', min_dist
            deallocate(dists)
        enddo
    
        close(22)
    end subroutine find_closest
    
    subroutine write_centers( fname, coords, smpd )
        character(len=*),           intent(in)    :: fname
        real,                       intent(in)    :: coords(:,:)
        real,                       intent(in)    :: smpd
        type(atoms) :: centers_pdb
        integer     :: cc
        call centers_pdb%new(size(coords, dim = 1), dummy=.true.)
        do cc = 1, size(coords, dim = 1)
            call centers_pdb%set_name(cc,params_glob%element)
            call centers_pdb%set_element(cc,params_glob%element)
            call centers_pdb%set_coord(cc,coords(cc,:)*smpd)
            !call centers_pdb%set_beta(cc,nano%atominfo(cc)%valid_corr) ! use per atom valid corr
            call centers_pdb%set_resnum(cc,cc)
        enddo
        call centers_pdb%writepdb(fname)
    end subroutine write_centers

    subroutine threshold_img(img_in, img_out, level)
        type(image), intent(in)  :: img_in
        type(image), intent(out) :: img_out
        integer,     intent(in)  :: level
        real              :: thres
        real, allocatable :: intensities(:)
        intensities = pack(img_in%get_rmat(), mask=.true.)
        call detect_peak_thres(size(intensities), level, intensities, thres)
        call img_out%copy(img_in)
        call img_out%zero_below(thres)
        call img_out%write('post_thresholding_map.mrc')
        deallocate(intensities)
    end subroutine threshold_img

    subroutine make_intensity_mask(img_in, mask_out, level)
        type(image),          intent(in)  :: img_in
        logical, allocatable, intent(out) :: mask_out(:,:,:)
        integer,              intent(in)  :: level
        integer           :: ldim(3)
        real, allocatable :: rmat(:,:,:)
        type(image)       :: thres_img
        ldim = img_in%get_ldim()
        allocate(mask_out(ldim(1), ldim(2), ldim(3)), source=.false.)
        call threshold_img(img_in, thres_img, level)
        rmat = thres_img%get_rmat()
        where (rmat > 0)
            mask_out = .true.
        end where
    end subroutine make_intensity_mask

    function logical_box(array, pos, boxsize) result(truth_value)
        logical, intent(in)  :: array(:,:,:)
        integer, intent(in)  :: pos(3)
        integer, intent(in)  :: boxsize
        logical              :: truth_value
        logical, allocatable :: small_box(:,:,:)
        integer              :: to(3), from(3)
        allocate(small_box(boxsize,boxsize,boxsize))
        from      = pos
        to        = pos + boxsize
        small_box = array(from(1):to(1),from(2):to(2),from(3):to(3))
        if(count(small_box) > 0) then 
            truth_value = .true.
        else
            truth_value = .false.
        end if
        deallocate(small_box)
    end function logical_box
    
end module simple_nano_picker_utils