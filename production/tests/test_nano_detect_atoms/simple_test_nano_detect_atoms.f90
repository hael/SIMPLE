module nano_picker_utils
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

    function avg_loc_sdev_3D( img, winsz ) result( asdev )
        type(image), intent(in) :: img
        integer,     intent(in) :: winsz
        integer           :: img_ldim(3), i, j, k, ir(2), jr(2), kr(2), npix, isz, jsz, ksz
        real              :: avg, asdev
        real, allocatable :: rmat(:,:,:), sdevs(:,:,:)
        img_ldim = img%get_ldim()
        allocate( rmat(img_ldim(1),img_ldim(2),img_ldim(3)))
        allocate(sdevs(img_ldim(1),img_ldim(2),img_ldim(3)))
        rmat = img%get_rmat()
        do i = 1, img_ldim(1)
            ir(1) = max(1,           i-winsz)
            ir(2) = min(img_ldim(1), i+winsz)
            isz   = ir(2) - ir(1) + 1
            do j = 1, img_ldim(2)
                jr(1) = max(1,           j-winsz)
                jr(2) = min(img_ldim(2), j+winsz)
                jsz   = jr(2) - jr(1) + 1
                do k = 1, img_ldim(3)
                    kr(1)        = max(1,           k-winsz)
                    kr(2)        = min(img_ldim(3), k+winsz)
                    ksz          = kr(2) - kr(1) + 1
                    npix         = isz * jsz * ksz
                    avg          = sum(rmat(ir(1):ir(2),jr(1):jr(2),kr(1):kr(2))) / real(npix)
                    sdevs(i,j,k) = sqrt(sum((rmat(ir(1):ir(2),jr(1):jr(2),kr(1):kr(2))-avg)**2.0) / real(npix - 1))
                enddo
            enddo
        enddo
        asdev = sum(sdevs) / real(img_ldim(1) * img_ldim(2) * img_ldim(3))
        deallocate(rmat,sdevs)
    end function avg_loc_sdev_3D

end module nano_picker_utils

module nano_detect_atoms
    !$ use omp_lib
    !$ use omp_lib_kinds
    include 'simple_lib.f08'
    use simple_image
    use simple_binimage
    use simple_atoms
    use simple_nanoparticle
    use simple_parameters
    use nano_picker_utils
    use simple_math
    use simple_linalg
    use simple_nanoparticle_utils
    use simple_defs_atoms
    use simple_stat
    use simple_aff_prop

    implicit none
    
    type :: nano_picker
        private
        character(len=2)         :: element
        character(len=100)       :: raw_filename, pdb_filename
        integer                  :: boxsize, ldim(3), nxyz_offset(3), offset, peak_thres_level
        integer, allocatable     :: inds_offset(:,:,:), positions(:,:)
        type(image)              :: simulated_atom, nano_img, sim_img
        type(image), allocatable :: convolved_atoms(:)
        real                     :: smpd, thres, radius
        real, allocatable        :: box_scores(:,:,:), loc_sdevs(:,:,:), center_positions(:,:)
        logical                  :: is_denoised

    contains
        procedure :: new
        procedure :: simulate_atom
        procedure :: setup_iterators
        procedure :: match_boxes
        procedure :: identify_threshold
        procedure :: distance_filter
        procedure :: remove_outliers
        procedure :: find_centers
        procedure :: calc_atom_stats
        procedure :: write_pdb
        procedure :: write_boximgs
        procedure :: write_positions
        procedure :: write_NP_image
        procedure :: write_corr_dist
        procedure :: compare_pick
        procedure :: refine_threshold
        procedure :: kill

    end type nano_picker

    contains

    subroutine new( self, smpd, element, filename, peak_thres_level, denoise )
        class(nano_picker), intent(inout) :: self
        real,               intent(in)    :: smpd
        character(len=*),   intent(in)    :: element
        character(len=100), intent(in)    :: filename
        integer,            intent(in)    :: peak_thres_level
        logical, optional,  intent(in)    :: denoise
        character(len=2)         :: el_ucase
        type(nanoparticle)       :: nano
        type(parameters), target :: params
        integer                  :: Z, nptcls
        call self%kill
        self%smpd         = smpd
        self%element      = element
        self%raw_filename = filename
        self%is_denoised  = .false.
        ! retrieve nano_img from filename and find ldim
        params_glob => params
        params_glob%element = self%element
        params_glob%smpd    = self%smpd
        el_ucase            = upperCase(trim(adjustl(params_glob%element)))
        call get_element_Z_and_radius(el_ucase, Z, self%radius)
        if( Z == 0 ) THROW_HARD('Unknown element: '//el_ucase)
        self%boxsize = round2even(self%radius / self%smpd) * 2
        call find_ldim_nptcls(self%raw_filename, self%ldim, nptcls, self%smpd)
        call self%nano_img%new(self%ldim, self%smpd)
        call self%nano_img%read(filename)
        call self%simulated_atom%new([self%boxsize,self%boxsize,self%boxsize], self%smpd)
        call self%sim_img%new(self%ldim, self%smpd)

        !call self%convolved_atoms(:)%new([self%boxsize,self%boxsize,self%boxsize],self%smpd)
        !call nano%new(filename)
        !call nano%get_img(self%nano_img)
        !self%ldim             = self%nano_img%get_ldim()

        self%peak_thres_level = peak_thres_level
        ! denoise nano_img if requested
        if( present(denoise) )then
            if( denoise )then
                call phasecorr_one_atom(self%nano_img, self%nano_img, self%element)
                self%is_denoised = .true.
            endif
        endif
    end subroutine new

    subroutine simulate_atom( self )
        class(nano_picker), intent(inout) :: self
        type(atoms) :: atom
        integer     :: Z, ldim_box(3)
        logical     :: l_err_atom
        call get_element_Z_and_radius(self%element, Z, self%radius)

        call atom%new(1)
        call atom%set_element(1,trim(self%element))
        ldim_box     = [self%boxsize,self%boxsize,self%boxsize]
        call atom%set_coord(1,(self%smpd*real(ldim_box)/2.)) ! make sure atom is in center of box
        !call self%simulated_atom%new(ldim_box,self%smpd)
        call atom%convolve(self%simulated_atom, cutoff=8*self%smpd)
        call self%simulated_atom%write('simulated_atom.mrc')
        !call self%simulated_atom%prenorm4real_corr(l_err_atom)
        call atom%kill
    end subroutine simulate_atom

    subroutine setup_iterators( self, offset )
        class(nano_picker), intent(inout) :: self
        integer,            intent(in)    :: offset
        integer :: nxyz(3), xind, yind, zind, nboxes
        ! set up picking infrastructure
        self%offset      = offset
        nxyz             = self%ldim - self%boxsize
        self%nxyz_offset = 0
        nboxes           = 0
        ! find number of boxes
        do xind = 0, nxyz(1), self%offset
            self%nxyz_offset(1) = self%nxyz_offset(1) + 1
            self%nxyz_offset(2) = 0
            do yind = 0, nxyz(2), self%offset
                self%nxyz_offset(2) = self%nxyz_offset(2) + 1
                self%nxyz_offset(3) = 0
                do zind = 0, nxyz(3), self%offset
                    self%nxyz_offset(3) = self%nxyz_offset(3) + 1
                    nboxes              = nboxes + 1
                enddo
            enddo
        enddo
        ! set up positions and inds_offset
        allocate(self%positions(nboxes,3),        source = 0)
        allocate(self%center_positions(nboxes,3), source = 0.)
        allocate(self%inds_offset(self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)), source=0)
        self%nxyz_offset = 0
        nboxes           = 0
        do xind = 0, nxyz(1), self%offset
            self%nxyz_offset(1) = self%nxyz_offset(1) + 1
            self%nxyz_offset(2) = 0
            do yind = 0, nxyz(2), self%offset
                self%nxyz_offset(2) = self%nxyz_offset(2) + 1
                self%nxyz_offset(3) = 0
                do zind = 0, nxyz(3), self%offset
                    self%nxyz_offset(3)      = self%nxyz_offset(3) + 1
                    nboxes                   = nboxes + 1
                    self%positions(nboxes,:) = [xind,yind,zind]
                    self%inds_offset(self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)) = nboxes
                enddo
            enddo
        enddo
        allocate(self%box_scores(self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)), source = -1.)
        allocate(self%loc_sdevs( self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)), source = -1.)
    end subroutine setup_iterators

    subroutine match_boxes( self, circle )
        class(nano_picker), intent(inout) :: self
        logical, optional,  intent(in)    :: circle
        type(image), allocatable :: boximgs(:)
        integer                  :: xoff, yoff, zoff, pos(3), pos_center(3), ithr, nthr, winsz, npix_in, npix_out1, npix_out2
        logical                  :: l_err_box, circle_here
        real                     :: maxrad, xyz(3)
        real,        allocatable :: pixels1(:), pixels2(:)
        logical                  :: outside
        if( present(circle) )then
            circle_here = circle
        else
            circle_here = .false.
        endif
        ! construct array of boximgs
        !$ nthr = omp_get_max_threads()
        allocate(boximgs(nthr))
        do ithr = 1,nthr
            call boximgs(ithr)%new([self%boxsize,self%boxsize,self%boxsize],self%smpd)
        enddo
        if( .not. circle_here )then
            ! use entire boxes for correlation scores
            ! iterate through positions in nanoparticle image, compare to simulated atom 
            !$omp parallel do schedule(static) collapse(3) default(shared) private(xoff,yoff,zoff,ithr,pos,l_err_box) proc_bind(close)
            do xoff = 1, self%nxyz_offset(1)
                do yoff = 1, self%nxyz_offset(2)
                    do zoff = 1, self%nxyz_offset(3)
                        ithr = omp_get_thread_num() + 1
                        pos  = self%positions(self%inds_offset(xoff,yoff,zoff),:)
                        call self%nano_img%window_slim( pos, self%boxsize, boximgs(ithr), outside)
                        !call boximgs(ithr)%prenorm4real_corr(l_err_box)
                        self%box_scores(xoff,yoff,zoff) = self%simulated_atom%real_corr(boximgs(ithr))
                        self%loc_sdevs( xoff,yoff,zoff) = avg_loc_sdev_3D(boximgs(ithr),self%offset)
                    enddo 
                enddo 
            enddo
            !$omp end parallel do
        else
        ! circular correlation
            maxrad    = (self%radius * 1.5) / self%smpd ! in pixels
            winsz     = ceiling(maxrad)
            npix_in   = (2 * winsz + 1)**3
            allocate(pixels1(npix_in), pixels2(npix_in), source=0.)
            ! !$omp parallel do schedule(static) collapse(3) default(shared) private(xoff,yoff,zoff,ithr,pos,pos_center,xyz,l_err_box) proc_bind(close)
            do xoff = 1, self%nxyz_offset(1)
                do yoff = 1, self%nxyz_offset(2)
                    do zoff = 1, self%nxyz_offset(3)
                        pos  = self%positions(self%inds_offset(xoff,yoff,zoff),:)
                        !pos_center = pos + [self%boxsize/2,self%boxsize/2,self%boxsize/2]
                        ithr = omp_get_thread_num() + 1
                        call self%nano_img%window_slim( pos, self%boxsize, boximgs(ithr), outside)
                        call boximgs(ithr)%norm_minmax
                        call boximgs(ithr)%masscen(xyz)
                        pos_center = pos + anint(xyz) + [self%boxsize/2,self%boxsize/2,self%boxsize/2]
                        do 
                            if( pos_center(1)-winsz < 1 .or. pos_center(2)-winsz < 1 .or. pos_center(3)-winsz < 1 )then
                                pos_center = pos_center + [1,1,1]
                            else
                                exit
                            endif
                        enddo
                        do 
                            if( pos_center(1)+winsz > self%ldim(1) .or. pos_center(2)+winsz > self%ldim(2) .or. pos_center(3)+winsz > self%ldim(3) )then
                                pos_center = pos_center - [1,1,1]
                            else
                                exit
                            endif
                        enddo
                        call self%nano_img%win2arr_rad(      pos_center(1),  pos_center(2),  pos_center(3),  winsz, npix_in, maxrad, npix_out1, pixels1)
                        call self%simulated_atom%win2arr_rad(self%boxsize/2, self%boxsize/2, self%boxsize/2, winsz, npix_in, maxrad, npix_out2, pixels2)
                        self%box_scores(xoff,yoff,zoff) = pearsn_serial(pixels1(:npix_out1),pixels2(:npix_out2))
                        self%loc_sdevs( xoff,yoff,zoff) = avg_loc_sdev_3D(boximgs(ithr),self%offset)
                    enddo 
                enddo 
            enddo
        endif
        ! !$omp end parallel do
        ! kill boximgs
        do ithr = 1, nthr
            call boximgs(ithr)%kill
        enddo
        ! deallocate(boximgs)
    end subroutine match_boxes

    subroutine identify_threshold( self, min_thres )
        class(nano_picker), intent(inout) :: self
        real, optional,     intent(in)    :: min_thres
        real, allocatable :: tmp(:)
        ! find peak thresholding value
        if( present(min_thres) )then
            tmp = pack(self%box_scores, mask=(self%box_scores > min_thres))
        else ! idk if this is something that makes sense to do..
            tmp = pack(self%box_scores, mask=(self%box_scores > -1 + 1e-10))
        endif
        call detect_peak_thres(size(tmp), size(tmp), self%peak_thres_level, tmp, self%thres)
        print *, 'Peak threshold is ', self%thres
        deallocate(tmp)
    end subroutine identify_threshold

    subroutine distance_filter( self, dist_thres )
        class(nano_picker), intent(inout) :: self
        real,   optional,   intent(in)    :: dist_thres
        real                 :: dist_thres_here, dist
        integer              :: nbox, ibox, jbox, xoff, yoff, zoff, npeaks, ipeak, loc
        integer, allocatable :: pos_inds(:)
        real,    allocatable :: pos_scores(:)
        logical, allocatable :: mask(:), selected_pos(:)
        logical              :: is_peak
        character(len=8)     :: crystal_system
        ! distance threshold
        if( present(dist_thres) )then
            dist_thres_here = dist_thres
        else
            dist_thres_here = self%offset
        endif
        pos_inds   = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres)
        pos_scores = pack(self%box_scores(:,:,:),   mask=self%box_scores(:,:,:) >= self%thres)
        nbox       = size(pos_inds)
        allocate(mask(nbox),         source=.false.)
        allocate(selected_pos(nbox), source=.true. )
        do ibox = 1, nbox
            mask = .false.
            ! identify boxes in neighborhood
            !$omp parallel do schedule(static) default(shared) private(jbox, dist) proc_bind(close)
            do jbox = 1, nbox
                dist = euclid(real(self%positions(pos_inds(ibox),:)),real(self%positions(pos_inds(jbox),:)))
                if( dist <= dist_thres_here ) mask(jbox) = .true.
            enddo
            !$omp end parallel do
            ! find highest correlation score in neighborhood
            loc = maxloc(pos_scores, mask=mask, dim=1)
            ! eliminate all but the best
            mask(loc) = .false.
            where( mask ) selected_pos = .false.
        enddo
        npeaks = count(selected_pos)
        print *, 'NPEAKS BEFORE DISTANCE FILTER = ', nbox
        print *, 'NPEAKS AFTER DISTANCE FILTER = ', npeaks
        ! update packed arrays
        pos_inds   = pack(pos_inds,   mask=selected_pos)
        pos_scores = pack(pos_scores, mask=selected_pos)
        ! update box scores
        !$omp parallel do schedule(static) collapse(3) default(shared) private(xoff,yoff,zoff,is_peak,ipeak)
        do xoff = 1, self%nxyz_offset(1)
            do yoff = 1, self%nxyz_offset(2)
                do zoff = 1, self%nxyz_offset(3)
                    is_peak = .false.
                    do ipeak = 1,npeaks
                        if( pos_inds(ipeak) == self%inds_offset(xoff,yoff,zoff) )then
                            is_peak = .true.
                            exit
                        endif
                    enddo
                    if( .not. is_peak ) self%box_scores(xoff,yoff,zoff) = -1.
                enddo
            enddo
        enddo
        !$omp end parallel do
        deallocate(pos_inds, pos_scores, mask, selected_pos)
    end subroutine distance_filter

    subroutine remove_outliers( self, ndev )
        class(nano_picker), intent(inout) :: self
        real,               intent(in)    :: ndev
        real, allocatable :: tmp(:)
        real              :: avg, sdev, t
        integer           :: npeaks, xoff, yoff, zoff
        tmp    = pack(self%loc_sdevs, mask = self%box_scores(:,:,:) >= self%thres .and. self%loc_sdevs(:,:,:) > 0.)
        call avg_sdev(tmp, avg, sdev)
        t      = avg + ndev * sdev
        npeaks = count(tmp < t)
        print *, 'NPEAKS AFTER REMOVE OUTLIERS = ', npeaks
        ! update box scores
        !$omp parallel do schedule(static) collapse(3) default(shared) private(xoff,yoff,zoff) proc_bind(close)
        do xoff = 1, self%nxyz_offset(1)
            do yoff = 1, self%nxyz_offset(2)
                do zoff = 1, self%nxyz_offset(3)
                    if( self%loc_sdevs(xoff,yoff,zoff) < t )then
                        ! it is a peak
                    else
                        self%box_scores(xoff,yoff,zoff) = -1
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine remove_outliers

    ! no longer writes the file
    subroutine find_centers( self )
        class(nano_picker), intent(inout) :: self
        type(image), allocatable :: atms_array(:)
        integer,     allocatable :: pos_inds(:)
        real,        allocatable :: coords(:,:)
        integer                  :: nbox, iimg, pos(3)
        logical                  :: outside
        ! make array of images containing the images of identified atoms and extract coordinates of peaks
        pos_inds = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres)
        nbox     = size(pos_inds, dim=1)
        allocate(coords(nbox,3))
        allocate(atms_array(nbox))
        if( .not. self%is_denoised )then
            if (allocated(self%convolved_atoms)) deallocate(self%convolved_atoms)
            allocate(self%convolved_atoms(nbox))
            call self%simulated_atom%fft()
            do iimg = 1, nbox
                pos = self%positions(pos_inds(iimg),:)
                call atms_array(iimg)%new([self%boxsize,self%boxsize,self%boxsize],self%smpd)
                call self%convolved_atoms(iimg)%new([self%boxsize,self%boxsize,self%boxsize],self%smpd)
                call self%nano_img%window_slim( pos, self%boxsize, atms_array(iimg), outside)
                !call atms_array(iimg)%write('boximgs/boximg_'//trim(int2str(iimg))//'.mrc')
                call atms_array(iimg)%fft()
                self%convolved_atoms(iimg) = atms_array(iimg)%conjg() * self%simulated_atom
                !call self%convolved_atoms(iimg)%write('boximgs_ft/boximg_ft_'//trim(int2str(iimg))//'.mrc')
                call self%convolved_atoms(iimg)%ifft()
                !call self%convolved_atoms(iimg)%write('boximgs_convolved/boximg_conv_'//trim(int2str(iimg))//'.mrc')
                call atms_array(iimg)%ifft()
                ! want coordinates of atoms to be at the center of the images
                call self%convolved_atoms(iimg)%norm_minmax
                call self%convolved_atoms(iimg)%masscen(coords(iimg,:)) 
                coords(iimg,:) = coords(iimg,:) + real(self%convolved_atoms(iimg)%get_ldim())/2. + pos !adjust center by size and position of box
                ! update center positions for chosen boxes
                self%center_positions(pos_inds(iimg),:) = coords(iimg,:)
            enddo
            call self%simulated_atom%ifft()
        else 
            do iimg = 1, nbox
                pos = self%positions(pos_inds(iimg),:)
                call atms_array(iimg)%new([self%boxsize,self%boxsize,self%boxsize],self%smpd)
                call self%nano_img%window_slim( pos, self%boxsize, atms_array(iimg), outside)
                ! want coordinates of atoms to be at the center of the images
                call atms_array(iimg)%norm_minmax
                call atms_array(iimg)%masscen(coords(iimg,:)) 
                coords(iimg,:) = coords(iimg,:) + real([self%boxsize,self%boxsize,self%boxsize])/2. + pos !adjust center by size and position of box
                ! update center positions for chosen boxes
                self%center_positions(pos_inds(iimg),:) = coords(iimg,:)
            enddo
            if( allocated(self%convolved_atoms) ) deallocate(self%convolved_atoms)
            allocate(self%convolved_atoms(nbox))
            self%convolved_atoms = atms_array
        endif
        deallocate(atms_array)
        deallocate(coords)
        deallocate(pos_inds)
    end subroutine find_centers

    subroutine calc_atom_stats(self)
        class(nano_picker), intent(inout) :: self
        type(nanoparticle)       :: nano
        type(image), allocatable :: boximgs(:)
        type(binimage)           :: bimg
        type(parameters), target :: params
        integer, allocatable     :: pos_inds(:), imat(:,:,:), imat_adj(:,:,:)
        integer                  :: nbox, ibox, int_pos(3), edge_pos(3), x, y, z, rad, xbox, ybox, zbox, nthr, ithr
        real, allocatable        :: rmat(:,:,:), rmat_flat(:)
        real                     :: local_thres
        logical                  :: outside
        params_glob => params
        params_glob%element = self%element
        params_glob%smpd    = self%smpd
        pos_inds            = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres)
        nbox                = size(pos_inds, dim=1)
        allocate(imat(    self%ldim(1),self%ldim(2),self%ldim(3)),source=0)
        allocate(imat_adj(self%ldim(1),self%ldim(2),self%ldim(3)),source=0)
        !$ nthr = omp_get_max_threads()
        allocate(boximgs(nthr))
        do ithr = 1,nthr
            call boximgs(ithr)%new([self%boxsize,self%boxsize,self%boxsize],self%smpd)
        enddo
        rad = anint((self%radius * 1.) / self%smpd) ! need to convert to pixels, give wiggle room
        !$omp parallel do schedule(static) default(shared) private(ibox,x,y,z,xbox,ybox,zbox,rmat,rmat_flat) proc_bind(close)
        do ibox = 1, nbox
            int_pos   = anint(self%center_positions(pos_inds(ibox),:))
            edge_pos  = int_pos - [rad,rad,rad]
            ithr      = omp_get_thread_num() + 1
            call self%nano_img%window_slim(edge_pos, 2*rad, boximgs(ithr), outside)
            rmat      = boximgs(ithr)%get_rmat()
            rmat_flat = pack(rmat, mask=.true.)
            call detect_peak_thres(size(rmat_flat), 2, rmat_flat, local_thres)
            xbox      = 0
            ybox      = 0
            zbox      = 0
            do x = int_pos(1) - rad, int_pos(1) + rad
                xbox = xbox + 1 ! keeps count of position in box, while iterator keeps track of position in larger image
                do y = int_pos(2) - rad, int_pos(2) + rad
                    ybox = ybox + 1
                    do z = int_pos(3) - rad, int_pos(3) + rad
                        zbox = zbox + 1
                        if( euclid(real(int_pos),real([x,y,z])) <= rad )then
                            if( rmat(xbox,ybox,zbox) >= local_thres )then
                                imat(x,y,z) = ibox
                            endif
                        endif
                    enddo
                    zbox = 0
                enddo
                ybox = 0
            enddo
            xbox = 0
            deallocate(rmat,rmat_flat)
        enddo
        !$omp end parallel do
        call bimg%new_bimg(self%ldim,self%smpd)
        where( imat >= 1 )
            imat_adj = 1
        elsewhere
            imat_adj = 0
        end where
        call bimg%set_imat(imat_adj)
        call bimg%write_bimg('CC.mrc')
        ! create new nanoparticle object
        call nano%new(trim(self%raw_filename))
        call nano%set_atomic_coords(trim(self%pdb_filename))
        call nano%set_img(trim(self%raw_filename),'img_raw')
        call nano%fillin_atominfo(imat=imat) ! should imat be in pixels or Angstroms? i think pixels
        call nano%write_csv_files
        call nano%kill
        deallocate(imat,boximgs)
    end subroutine calc_atom_stats

    ! input filename with no extension
    subroutine write_pdb( self, filename )
        class(nano_picker),         intent(inout) :: self
        character(len=*), optional, intent(in)    :: filename
        integer, allocatable :: pos_inds(:)
        real,    allocatable :: coords(:,:)
        integer :: nbox, iimg
        real :: pos(3)
        ! make array of images containing the images of identified atoms and extract coordinates of peaks
        pos_inds = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres)
        nbox     = size(pos_inds, dim=1)
        allocate(coords(nbox,3))
        do iimg = 1, nbox
            pos = self%center_positions(pos_inds(iimg),:)
            coords(iimg,:) = pos
        enddo
        if( present(filename) )then
            call write_centers(filename,coords,self%smpd)
            self%pdb_filename = trim(filename)//'.pdb'
        else
            call write_centers('test_atomic_centers',coords,self%smpd)
            self%pdb_filename = 'test_atomic_centers.pdb'
        endif
        deallocate(coords)
        deallocate(pos_inds)
    end subroutine write_pdb

    subroutine write_boximgs( self, foldername )
        class(nano_picker),          intent(inout) :: self
        character(len=*), optional,  intent(in)    :: foldername
        integer              :: iimg, nbox
        integer, allocatable :: pos_inds(:)
        pos_inds   = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres)
        nbox       = size(pos_inds, dim=1)
        do iimg = 1, nbox
            if( present(foldername) )then 
                call self%convolved_atoms(iimg)%write(trim(adjustl(foldername))//'/boximg_'//trim(int2str(iimg))//'.mrc')
            else
                call self%convolved_atoms(iimg)%write('boximgs/boximg_'//trim(int2str(iimg))//'.mrc')
            endif
        enddo
        deallocate(pos_inds)
    end subroutine write_boximgs

    subroutine write_positions( self, filename )
        class(nano_picker), intent(inout) :: self
        character(len=*),   intent(in)    :: filename
        integer,     allocatable :: pos_inds(:)
        real,        allocatable :: coords(:,:)
        integer                  :: nbox, ipos, i, j
        pos_inds = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres)
        nbox     = size(pos_inds)
        allocate(coords(nbox,3))
        do ipos = 1, nbox
            coords(ipos,:) = self%center_positions(pos_inds(ipos),:) * self%smpd
        enddo
        open(unit=25, file=filename, status='replace', action='write')
        do i = 1, nbox
            write(25, '(I9,a)', advance='no') pos_inds(i), ','
            do j = 1, 3
                if( j /= 1 ) write(25, '(A)', advance='no') ', '
                write(25, '(F10.3)', advance='no') coords(i, j)
            enddo
            write(25, *)
        enddo
        close(25)
    end subroutine write_positions

    subroutine write_NP_image( self, ref_img_name, sim_img_name )
        class(nano_picker), intent(inout) :: self
        character(len=*),   intent(in)    :: ref_img_name, sim_img_name
        type(nanoparticle)       :: nano
        type(image)              :: sim_img
        type(parameters), target :: params
        params_glob => params
        params_glob%element = self%element
        params_glob%smpd    = self%smpd
        call self%write_pdb('simulate_NP')
        call nano%new(trim(ref_img_name))
        call nano%set_atomic_coords('simulate_NP.pdb')
        call nano%simulate_atoms(simatms=sim_img)
        call sim_img%write(sim_img_name)
        self%sim_img = sim_img
        call nano%kill
    end subroutine write_NP_image

    subroutine write_corr_dist( self, csv_name )
        class(nano_picker), intent(inout) :: self
        character(len=*),   intent(in)    :: csv_name
        real,    allocatable :: pos_scores(:), lower_half_scores(:), upper_half_scores(:)
        integer, allocatable :: pos_inds(:)
        real                 :: Q1, mid, Q3, IQR, mean
        integer              :: i
        pos_inds          = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres)
        pos_scores        = pack(self%box_scores(:,:,:),   mask=self%box_scores(:,:,:) >= self%thres)
        mid               = median(pos_scores)
        lower_half_scores = pack(pos_scores(:), pos_scores(:) < mid)
        upper_half_scores = pack(pos_scores(:), pos_scores(:) > mid)
        Q1                = median(lower_half_scores)
        Q3                = median(upper_half_scores)
        IQR               = Q3 - Q1
        mean = sum(pos_scores) / size(pos_scores)
        print *, 'SUMMARY STATISTICS OF ATOMIC CORRELATION SCORES'
        print *, 'Q1 = ', Q1
        print *, 'MEDIAN = ', mid
        print *, 'Q3 = ', Q3
        print *, 'IQR = ', IQR
        print *, 'MEAN = ', mean
        open(unit=99,file=trim(csv_name))
        do i = 1, size(pos_scores)
            write(99,'(1x,f4.3)') pos_scores(i)
        enddo
        close(99)
    end subroutine write_corr_dist

    ! input both pdbfile_* with .pdb extension
    subroutine compare_pick( self, pdbfile_ref, pdbfile_exp )
        class(nano_picker),         intent(inout) :: self
        character(len=*),           intent(in)    :: pdbfile_ref
        character(len=*), optional, intent(in)    :: pdbfile_exp
        real,    allocatable :: pdb_ref_coords(:,:), pdb_exp_coords(:,:), distances(:)
        integer, allocatable :: pos_inds(:)
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
    
    subroutine refine_threshold( self, num_thres, ref_pdb_name, ref_img_name, max_thres )
        class(nano_picker), intent(inout) :: self
        integer,            intent(in)    :: num_thres
        character(len=*),   intent(in)    :: ref_pdb_name, ref_img_name
        real,    optional,  intent(in)    :: max_thres
        type(nanoparticle)       :: nano_ref, nano_exp
        type(image)              :: ref_NP,   exp_NP
        type(parameters), target :: params
        real                     :: thresholds(num_thres), thres_corrs(num_thres), max_thres_here, step
        integer                  :: i, optimal_index(1), num_pos
        integer, allocatable     :: pos_inds(:)
        print *, 'REFINE_THRESHOLD'
        if( present(max_thres) )then
            max_thres_here = max_thres
        else
            max_thres_here = 0.5
        endif
        step = (max_thres_here - self%thres) / (num_thres - 1)
        ! set up array of potential thresholds
        do i = 1, num_thres
            thresholds(i) = self%thres + step * (i-1)
        enddo
        ! simulate nanoparticle with initial pdb file (for reference / comparison)
        ! params_glob has to be set because of the way simple_nanoparticle is set up
        params_glob => params
        params_glob%element = self%element
        params_glob%smpd    = self%smpd
        call nano_ref%new(trim(ref_img_name))
        call nano_ref%set_atomic_coords(trim(ref_pdb_name))
        call nano_ref%simulate_atoms(simatms=ref_NP)
        ! iterate through following steps:
        ! 1. remove boxes with correlations below each threshold
        ! 2. call find centers and write_pdb
        ! 3. use the resulting pdb file to simulate nanoparticle
        ! 4. calculate correlation between this simulated nanoparticle and original? simulated nanoparticle
        ! 5. save correlations in array, at end will find maximum and return corresponding threshold
        do i = 1, num_thres
            self%thres     = thresholds(i) ! need to set self%thres because it is called in multiple subroutines
            call self%find_centers
            call self%write_pdb('sim_centers')
            call nano_exp%new(trim(ref_img_name))
            call nano_exp%set_atomic_coords('sim_centers.pdb')
            call nano_exp%simulate_atoms(simatms=exp_NP)
            thres_corrs(i) = ref_NP%real_corr(exp_NP)
            !thres_corrs(i) = self%nano_img%real_corr(exp_NP) ! see if correlating with raw img helps
            call nano_exp%kill
        enddo
        optimal_index = maxloc(thres_corrs)
        self%thres    = thresholds(optimal_index(1))
        call self%find_centers ! call again to set positions to the optimal
        pos_inds      = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres)
        num_pos       = size(pos_inds)
        print *, 'OPTIMAL THRESHOLD = ', self%thres
        print *, 'OPTIMAL CORRELATION = ', thres_corrs(optimal_index(1))
        print *, 'NUMBER POSITIONS = ', num_pos
        call nano_ref%kill
    end subroutine refine_threshold

    subroutine kill( self )
        class(nano_picker), intent(inout) :: self
        integer :: i
        call self%nano_img%kill
        call self%sim_img%kill
        call self%simulated_atom%kill
        if(allocated(self%convolved_atoms))then
            do i = 1, size(self%convolved_atoms)
                call self%convolved_atoms(i)%kill
            enddo
            deallocate(self%convolved_atoms)
        endif
        if(allocated(self%positions)       ) deallocate(self%positions)
        if(allocated(self%center_positions)) deallocate(self%center_positions)
        if(allocated(self%inds_offset)     ) deallocate(self%inds_offset)
        if(allocated(self%box_scores)      ) deallocate(self%box_scores)
    end subroutine kill

end module nano_detect_atoms
    
program simple_test_nano_detect_atoms
include 'simple_lib.f08'
use nano_detect_atoms
use nano_picker_utils
use simple_nanoparticle
use simple_nanoparticle_utils
use simple_image
use simple_parameters
use simple_strings, only: int2str
implicit none
#include "simple_local_flags.inc"

type(nano_picker)  :: test_sim
type(nano_picker)  :: test_exp4
type(nanoparticle) :: nano
type(image)        :: simulated_NP
real               :: smpd, dist_thres
character(len=2)   :: element
character(len=100) :: filename_exp, filename_sim, pdbfile_ref
character(STDLEN)  :: timestr
integer            :: offset, peak_thres_level, startTime, stopTime, subStart, subStop, ldim(3)
type(parameters), target :: params

! keeping track of how long program takes
startTime        = real(time())

! Inputs
filename_exp     = 'rec_merged.mrc'
filename_sim     = 'simulated_NP.mrc'
!pdbfile_ref     = 'ATMS.pdb'
pdbfile_ref      = 'reference.pdb'
element          = 'PT'
smpd             = 0.358
offset           = 2
peak_thres_level = 2
dist_thres       = 3.

! simulate nanoparticle
! params_glob has to be set because of the way simple_nanoparticle is set up
params_glob => params
params_glob%element = element
params_glob%smpd    = smpd

call nano%new(trim(filename_exp))
call nano%set_atomic_coords(trim(pdbfile_ref))
call nano%simulate_atoms(simatms=simulated_NP)
call nano%get_ldim(ldim) 
call simulated_NP%new(ldim, smpd)
call simulated_NP%write(trim(filename_sim))

subStart = real(time())

call test_exp4%new(smpd, element, filename_exp, peak_thres_level, denoise=.false.)
call test_exp4%simulate_atom()
call test_exp4%setup_iterators(offset)
call test_exp4%match_boxes(circle=.true.)
call test_exp4%identify_threshold()
call test_exp4%find_centers()
call test_exp4%distance_filter(dist_thres)
call test_exp4%refine_threshold(100,pdbfile_ref,filename_sim,max_thres=0.75)
call test_exp4%remove_outliers(3.)
call test_exp4%write_boximgs(foldername='boximgs')
! OUTPUT FILES
call test_exp4%write_pdb('experimental_centers')
call test_exp4%compare_pick('experimental_centers.pdb',trim(pdbfile_ref))
call test_exp4%write_NP_image(filename_sim,'result.mrc')
call test_exp4%calc_atom_stats
call test_exp4%kill
subStop = real(time())
print *, 'TEST 1 RUNTIME = ', (subStop - subStart), ' s'

print *, ' '
stopTime = time()
print *, 'TOTAL RUNTIME = ', (stopTime - startTime), ' s'

end program simple_test_nano_detect_atoms
