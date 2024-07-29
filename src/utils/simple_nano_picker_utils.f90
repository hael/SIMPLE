module simple_nano_picker_utils
include 'simple_lib.f08'
use simple_image
use simple_atoms
use simple_parameters
use simple_defs_atoms
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
            call centers_pdb%set_coord(cc,(coords(cc,:))*smpd)
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

    ! intensity thresholding based on Otsu's method
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

    ! intensity thresholding based on nanoparticle properties
    subroutine make_intensity_mask_2(img_in, mask_out, element, NP_diam)
        type(image),          intent(in)  :: img_in
        logical, allocatable, intent(out) :: mask_out(:,:,:)
        character(len=*),     intent(in)  :: element
        real,                 intent(in)  :: NP_diam ! approx. diameter of nanoparticle
        character(len=2)  :: el_ucase
        character(len=8)  :: crystal_system
        real              :: msksq, a, ha, x, y, z, center(3), smpd, radius, radius_vx, sphere_vol, total_vol, thres
        real              :: x1, x2, x3, y1, y2, y3, z1, z2, z3
        real, allocatable :: intensities_flat(:), rmat(:,:,:)
        integer           :: i, j, k, n, ncubes, box, ldim(3), ZZ, nvx, last_index
        print *, 'NP_diam = ', NP_diam
        msksq    = (NP_diam / 2.)**2.
        print *, 'msksq = ', msksq
        el_ucase = uppercase(trim(adjustl(element)))
        call get_lattice_params(el_ucase, crystal_system, a)
        print *, 'crystal_system = ', crystal_system
        print *, 'el_ucase = ', el_ucase
        print *, 'a = ', a
        ha       = a / 2.
        ldim     = img_in%get_ldim()
        box      = ldim(1)
        print *, 'box = ', box
        smpd     = img_in%get_smpd()
        ncubes   = floor(real(box) * smpd / a)
        ! atoms at edges
        n = 0
        center = (real([box,box,box]/2 + 1) - 1.) * smpd
        ! find number of atoms in startvol
        select case( trim(crystal_system) )
            case('fcc')
                ! count
                do i = 1, ncubes
                    x  = real(i - 1) * a
                    x1 = x + ha; x2 = x; x3 = x + ha
                    do j = 1, ncubes
                        y  = real(j - 1) * a
                        y1 = y + ha; y2 = y + ha; y3 = y
                        do k = 1, ncubes
                            z  = real(k - 1) * a
                            z1 = z; z2 = z + ha; z3 = z + ha
                            ! edge
                            if( sum(([x ,y ,z ] - center)**2.) <= msksq ) n = n+1
                            ! faces
                            if( sum(([x1,y1,z1] - center)**2.) <= msksq ) n = n+1
                            if( sum(([x2,y2,z2] - center)**2.) <= msksq ) n = n+1
                            if( sum(([x3,y3,z3] - center)**2.) <= msksq ) n = n+1
                        enddo
                    enddo
                enddo
            case('bcc')
                ! count
                do i = 1, ncubes
                    x  = real(i - 1) * a; x1 = x + ha
                    do j = 1, ncubes
                        y  = real(j - 1) * a; y1 = y + ha
                        do k= 1, ncubes
                            z  = real(k - 1) * a; z1 = z + ha
                            ! edge
                            if( sum(([x ,y ,z ] - center)**2.) <= msksq ) n = n + 1
                            ! body-centered
                            if( sum(([x1,y1,z1] - center)**2.) <= msksq ) n = n + 1
                        enddo
                    enddo
                enddo
            case('rocksalt')
                ! count
                do i = 1, ncubes
                    x  = real(i - 1) * a; x1 = x + ha
                    do j = 1, ncubes
                        y  = real(j - 1) * a; y1 = y + ha
                        do k=1,ncubes
                            z  = real(k - 1) * a; z1 = z + ha
                            if( sum(([x ,y ,z ] - center)**2.) <= msksq ) n = n + 1
                            if( sum(([x1,y1,z ] - center)**2.) <= msksq ) n = n + 1
                            if( sum(([x ,y1,z1] - center)**2.) <= msksq ) n = n + 1
                            if( sum(([x1,y ,z1] - center)**2.) <= msksq ) n = n + 1
                            if( sum(([x1,y1,z1] - center)**2.) <= msksq ) n = n + 1
                            if( sum(([x1,y ,z ] - center)**2.) <= msksq ) n = n + 1
                            if( sum(([x ,y1,z ] - center)**2.) <= msksq ) n = n + 1
                            if( sum(([x ,y ,z1] - center)**2.) <= msksq ) n = n + 1
                        enddo
                    enddo
                enddo
        end select
        print *, 'n = ', n
        if (n .eq. 0) print *, 'Error! n should be greater than 0: make_intensity_mask_2.'
        call get_element_Z_and_radius(trim(adjustl(element)), ZZ, radius)
        radius_vx  = radius / smpd
        sphere_vol = PI * (4./3.) * radius_vx**3
        total_vol  = sphere_vol * n
        print *, 'total_vol = ', total_vol
        nvx        = anint(total_vol)
        ! find intensity threshold value based on highest 'nvx' intensity values
        rmat       = img_in%get_rmat()
        intensities_flat = pack(rmat, mask=.true.)
        call hpsort(intensities_flat)                   ! sort array (low to high)
        last_index = size(intensities_flat)             ! last_index is position of largest intensity value in sorted array
        thres      = intensities_flat(last_index - nvx) ! last_index - nvx position of sorted array is threshold
        ! make intensity logical mask, with positions with intensities greater than the threshold being set to true
        allocate(mask_out(ldim(1), ldim(2), ldim(3)), source=.false.)
        where (rmat > thres)
            mask_out = .true.
        end where
        print *, 'count(mask_out) = ', count(mask_out)
        deallocate(intensities_flat)
    end subroutine make_intensity_mask_2

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

    function one_atom_valid_corr(coord, maxrad, raw_img, simatms) result(valid_corr)
        real,        intent(in) :: coord(3)
        real,        intent(in) :: maxrad
        type(image), intent(in) :: raw_img, simatms
        real              :: valid_corr ! ouput
        real, allocatable :: pixels1(:), pixels2(:)
        integer           :: winsz, npix_in, npix_out1, npix_out2, ijk(3)
        winsz     = ceiling(maxrad)
        npix_in   = (2 * winsz + 1)**3
        allocate(pixels1(npix_in), pixels2(npix_in), source=0.)
        ijk = nint(coord)
        call raw_img%win2arr_rad(ijk(1), ijk(2), ijk(3), winsz, npix_in, maxrad, npix_out1, pixels1)
        call simatms%win2arr_rad(ijk(1), ijk(2), ijk(3), winsz, npix_in, maxrad, npix_out2, pixels2)
        valid_corr = pearsn_serial(pixels1(:npix_out1),pixels2(:npix_out2))
        deallocate(pixels1,pixels2)
    end function one_atom_valid_corr
    
end module simple_nano_picker_utils