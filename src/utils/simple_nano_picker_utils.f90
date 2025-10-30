module simple_nano_picker_utils
include 'simple_lib.f08'
use simple_srch_sort_loc, only : hpsort
use simple_atoms
use simple_defs_atoms
use simple_image
use simple_nanoparticle_utils
use simple_parameters
use simple_segmentation
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

       ! input both pdbfile_* with .pdb extension
    subroutine compare_pick( pdbfile_ref, pdbfile_exp )
        character(len=*),           intent(in)    :: pdbfile_ref
        character(len=*), optional, intent(in)    :: pdbfile_exp
        real,    allocatable :: pdb_ref_coords(:,:), pdb_exp_coords(:,:), distances(:)
        integer              :: iostat, i
        call read_pdb2matrix(trim(pdbfile_ref), pdb_ref_coords)
        if( present(pdbfile_exp) )then 
            call read_pdb2matrix(trim(pdbfile_exp),pdb_exp_coords)
        else
            open(unit = 40, file='test_atomic_centers.pdb', iostat=iostat)
            if( iostat /= 0 )then
                print *, 'compare_pick: test_atomic_centers.pdb does not exist, please enter valid filename for pdbfile_exp'
                close(40)
                return
            endif
            call read_pdb2matrix('test_atomic_centers.pdb',pdb_exp_coords)
            close(40)
        endif
        allocate(distances(max(size(pdb_ref_coords,dim=2),size(pdb_exp_coords,dim=2))))
        call find_closest(pdb_ref_coords,pdb_exp_coords,size(pdb_ref_coords,dim=2),size(pdb_exp_coords,dim=2),distances)
        print *, 'AVG DISTANCE = ', sum(distances)/size(distances)
    end subroutine compare_pick
    
    subroutine write_centers( fname, coords, smpd )
        character(len=*),           intent(in)    :: fname
        real,                       intent(in)    :: coords(:,:)
        real,                       intent(in)    :: smpd
        type(atoms) :: centers_pdb
        integer     :: cc
        call centers_pdb%new(size(coords, dim = 1), dummy=.true.)
        open(unit=45,file=trim(fname)//'.csv')
        do cc = 1, size(coords, dim = 1)
            call centers_pdb%set_name(cc,' '//params_glob%element)
            call centers_pdb%set_num(cc,cc)
            call centers_pdb%set_element(cc,params_glob%element)
            call centers_pdb%set_coord(cc,(coords(cc,:)-1)*smpd)
            !call centers_pdb%set_beta(cc,nano%atominfo(cc)%valid_corr) ! use per atom valid corr
            call centers_pdb%set_resnum(cc,1)
            write(45,'(f8.4,a,f8.4,a,f8.4)') coords(cc,1), ',', coords(cc,2), ',', coords(cc,3)-1
        enddo
        close(45)
        call centers_pdb%writepdb(trim(fname)//'.pdb')
    end subroutine write_centers

    subroutine threshold_img(img_in, img_out, level, threshold_out)
        type(image),    intent(in)  :: img_in
        type(image),    intent(out) :: img_out
        integer,        intent(in)  :: level
        real, optional, intent(out) :: threshold_out
        real              :: thres
        real, allocatable :: intensities(:)
        intensities = pack(img_in%get_rmat(), mask=img_in%get_rmat() > 0.)
        call detect_peak_thres(size(intensities), level, intensities, thres)
        call img_out%copy(img_in)
        call img_out%zero_below(thres)
        call img_out%write('post_thresholding_map.mrc')
        deallocate(intensities)
        if (present(threshold_out)) threshold_out = thres
    end subroutine threshold_img

    ! intensity thresholding based on Otsu's method
    subroutine make_intensity_mask(img_in, mask_out, level, intensity_thres)
        type(image),          intent(in)  :: img_in
        logical, allocatable, intent(out) :: mask_out(:,:,:)
        integer,              intent(in)  :: level
        real,    optional,    intent(out) :: intensity_thres
        integer           :: ldim(3)
        real, allocatable :: rmat(:,:,:)
        type(image)       :: thres_img
        ldim = img_in%get_ldim()
        allocate(mask_out(ldim(1), ldim(2), ldim(3)), source=.false.)
        call threshold_img(img_in, thres_img, level, intensity_thres)
        rmat = thres_img%get_rmat()
        where (rmat > 0)
            mask_out = .true.
        end where
    end subroutine make_intensity_mask

    ! intensity thresholding based on nanoparticle properties
    subroutine make_intensity_mask_2(img_in, mask_out, element, NP_diam, level, intensity_thres)
        type(image),          intent(in)  :: img_in
        logical, allocatable, intent(out) :: mask_out(:,:,:)
        character(len=*),     intent(in)  :: element
        real,                 intent(in)  :: NP_diam ! approx. diameter of nanoparticle
        integer, optional,    intent(in)  :: level
        real,    optional,    intent(out) :: intensity_thres
        type(image)       :: img_copy
        character(len=2)  :: el_ucase
        character(len=8)  :: crystal_system
        real              :: msksq, a, ha, x, y, z, center(3), smpd, radius, radius_vx, sphere_vol, total_vol, thres
        real              :: x1, x2, x3, y1, y2, y3, z1, z2, z3
        real, allocatable :: intensities_flat(:), rmat(:,:,:)
        integer           :: i, j, k, n, ncubes, box, ldim(3), ZZ, nvx, last_index
        msksq    = (NP_diam / 2.)**2.
        el_ucase = uppercase(trim(adjustl(element)))
        call get_lattice_params(el_ucase, crystal_system, a)
        ha       = a / 2.
        ldim     = img_in%get_ldim()
        box      = ldim(1)
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
        if (n .eq. 0) THROW_HARD( 'Error! n should be greater than 0: make_intensity_mask_2.' )
        call get_element_Z_and_radius(element, ZZ, radius)
        radius_vx  = radius / smpd
        sphere_vol = PI * (4./3.) * radius_vx**3
        total_vol  = sphere_vol * n
        nvx        = anint(total_vol)
        ! find intensity threshold value based on highest 'nvx' intensity values
        rmat       = img_in%get_rmat()
        intensities_flat = pack(rmat, mask=.true.)
        call hpsort(intensities_flat)                   ! sort array (low to high)
        last_index = size(intensities_flat)             ! last_index is position of largest intensity value in sorted array
        thres      = intensities_flat(last_index - nvx) ! last_index - nvx position of sorted array is threshold
        print *, 'initial thres = ', thres
        call img_copy%copy(img_in)
        call img_copy%zero_below(thres)
        call img_copy%write('post_thresholding_map2.mrc')
        if (present(level)) then
            call make_intensity_mask(img_copy, mask_out, level, intensity_thres)
        else
            ! make intensity logical mask, with positions with intensities greater than the threshold being set to true
            allocate(mask_out(ldim(1), ldim(2), ldim(3)), source=.false.)
            where (rmat > thres)
                mask_out = .true.
            end where
            intensity_thres = thres
        end if
        deallocate(intensities_flat)
    end subroutine make_intensity_mask_2

    function logical_box(array, pos, boxsize, num_px) result(truth_value)
        logical, intent(in)           :: array(:,:,:)
        integer, intent(in)           :: pos(3)
        integer, intent(in)           :: boxsize
        integer, optional, intent(in) :: num_px
        logical              :: truth_value
        logical, allocatable :: small_box(:,:,:)
        integer              :: to(3), from(3), num_px_here
        allocate(small_box(boxsize,boxsize,boxsize))
        if (present(num_px)) then
            num_px_here = num_px
        else
            num_px_here = 0
        end if
        from      = pos
        to        = pos + boxsize
        small_box = array(from(1):to(1),from(2):to(2),from(3):to(3))
        if(count(small_box) > num_px_here) then 
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

    ! finds per-atom valid correlation for a subset of coordinates, given a raw and simulated image
    ! input coords in PIXELS (for now, can change this later)
    subroutine find_corr_by_coords(coords, maxrad, raw_img, simatms, smpd, corrs, filename)
        real,                       intent(in)  :: coords(:,:) !should be dimensions 3 x n, where n is number of positions
        real,                       intent(in)  :: maxrad
        type(image),                intent(in)  :: raw_img, simatms
        real,                       intent(in)  :: smpd
        character(len=*), optional, intent(in)  :: filename
        real,          allocatable, intent(out) :: corrs(:)
        integer           :: npos, ipos
        npos = size(coords, dim=2)
        allocate(corrs(npos))
        do ipos = 1, npos
            corrs(ipos) = one_atom_valid_corr(coords(:,ipos), maxrad, raw_img, simatms)
        end do
        if(present(filename)) then
            open(unit=80,file=filename)
            do ipos = 1, npos
                write(80,'(f8.4,a,f8.4,a,f8.4,a,f6.5)') (coords(1,ipos) - 1) * smpd, ',', (coords(2,ipos) - 1) * smpd, ',', (coords(3,ipos) - 1) * smpd, ',', corrs(ipos)
            end do
            close(80)
        end if
    end subroutine find_corr_by_coords
    
end module simple_nano_picker_utils