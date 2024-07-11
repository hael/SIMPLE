module simple_nano_detect_atoms
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image
use simple_binimage
use simple_atoms
use simple_nanoparticle
use simple_parameters
use simple_nano_picker_utils
use simple_math
use simple_linalg
use simple_nanoparticle_utils
use simple_defs_atoms
use simple_stat
#include "simple_local_flags.inc"

    implicit none
    
    type :: nano_picker
        private
        character(len=2)         :: element
        character(len=100)       :: raw_filename, pdb_filename
        integer                  :: boxsize, ldim(3), nxyz_offset(3), offset, peak_thres_level
        integer, allocatable     :: inds_offset(:,:,:), positions(:,:)
        type(image)              :: simulated_atom, nano_img
        type(image), allocatable :: convolved_atoms(:)
        real                     :: smpd, thres, radius, euclid_thres, mask_radius
        real, allocatable        :: box_scores(:,:,:), loc_sdevs(:,:,:), avg_int(:,:,:), euclid_dists(:,:,:), center_positions(:,:)
        logical                  :: is_denoised, use_euclids, has_mask, wrote_pdb
        logical, allocatable     :: msk(:,:,:)

    contains
        procedure :: new
        procedure :: simulate_atom
        procedure :: setup_iterators
        procedure :: match_boxes
        procedure :: one_box_corr
        procedure :: identify_threshold
        procedure :: identify_high_scores
        procedure :: apply_threshold
        procedure :: distance_filter
        procedure :: remove_outliers_sdev
        procedure :: remove_outliers_position
        procedure :: find_centers
        procedure :: calc_atom_stats
        procedure :: calc_per_atom_corr
        procedure :: write_pdb
        procedure :: write_boximgs
        procedure :: write_positions_and_scores
        procedure :: write_NP_image
        procedure :: write_dist
        procedure :: extract_img_from_pos
        procedure :: compare_pick
        procedure :: refine_threshold
        procedure :: kill

    end type nano_picker

    contains

    subroutine new( self, smpd, element, raw_filename, peak_thres_level, offset, denoise, use_euclids, mskdiam)
        class(nano_picker), intent(inout) :: self
        real,               intent(in)    :: smpd
        character(len=*),   intent(in)    :: element
        character(len=100), intent(in)    :: raw_filename
        integer,            intent(in)    :: peak_thres_level, offset
        logical, optional,  intent(in)    :: denoise
        logical, optional,  intent(in)    :: use_euclids
        real,    optional,  intent(in)    :: mskdiam ! Angstroms
        character(len=2)         :: el_ucase
        integer                  :: Z, nptcls
        logical                  :: outside
        type(parameters), target :: params
        call self%kill
        self%smpd           = smpd
        self%element        = element
        self%raw_filename   = raw_filename
        self%offset         = offset
        self%is_denoised    = .false.
        self%use_euclids    = .false.
        self%wrote_pdb      = .false.
        self%has_mask       = .false.
        params_glob         => params
        params_glob%element = self%element
        params_glob%smpd    = self%smpd
        ! retrieve nano_img from filename and find ldim
        el_ucase            = upperCase(trim(adjustl(params_glob%element)))
        call get_element_Z_and_radius(el_ucase, Z, self%radius)
        if( Z == 0 ) THROW_HARD('Unknown element: '//el_ucase)
        ! hardcode
        self%boxsize = round2even(self%radius / self%smpd) * 2
        call find_ldim_nptcls(self%raw_filename, self%ldim, nptcls, self%smpd)
        call self%nano_img%new(self%ldim, self%smpd)
        call self%nano_img%read(self%raw_filename)
        call self%simulated_atom%new([self%boxsize,self%boxsize,self%boxsize], self%smpd)
        self%peak_thres_level = peak_thres_level
        self%thres = -0.999 ! initial value, will be updated later
        ! denoise nano_img if requested
        if( present(denoise) )then
            if( denoise )then
                call phasecorr_one_atom(self%nano_img, self%element)
                self%is_denoised = .true.
                call self%nano_img%write('denoised.mrc')
            endif
        endif
        if( present(use_euclids) ) then
            self%use_euclids = use_euclids
        end if
        if( present(mskdiam) ) then
            self%has_mask = .true.
            self%mask_radius   = (mskdiam / self%smpd) / 2
        end if
    end subroutine new

    subroutine simulate_atom( self )
        class(nano_picker), intent(inout) :: self
        type(atoms)       :: atom
        integer           :: Z, ldim_box(3)
        logical           :: l_err_atom
        real, allocatable :: rmat(:,:,:), rmat_simatm(:,:,:)
        real              :: rmat_min, rmat_max, delta
        call get_element_Z_and_radius(self%element, Z, self%radius)
        if (Z == 0) THROW_HARD('Unknown element : '//self%element)
        call atom%new(1)
        call atom%set_element(1,trim(self%element))
        ldim_box  = [self%boxsize,self%boxsize,self%boxsize]
        call atom%set_coord(1,(self%smpd*real(ldim_box)/2.)) ! make sure atom is in center of box
        call atom%convolve(self%simulated_atom, cutoff=8*self%smpd)
        ! normalize simulated atom
        rmat         = self%nano_img%get_rmat()
        rmat_min     = minval(rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        rmat_max     = maxval(rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        delta        = rmat_max - rmat_min
        call self%simulated_atom%norm_minmax
        rmat_simatm  = self%simulated_atom%get_rmat()
        rmat_simatm  = rmat_simatm * delta
        rmat_simatm  = rmat_simatm + rmat_min
        call self%simulated_atom%set_rmat(rmat_simatm, ft=.false.)
        ! output 
        call self%simulated_atom%write('simulated_atom.mrc')
        call atom%kill
        deallocate(rmat,rmat_simatm)
    end subroutine simulate_atom

    subroutine setup_iterators( self )
        class(nano_picker), intent(inout) :: self
        integer :: nxyz(3), xind, yind, zind, nboxes
        ! set up picking infrastructure
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
        allocate(self%msk(self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)), source=.true.)
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
                    if(self%has_mask) then
                        if(euclid(real(self%ldim / 2), real([xind,yind,zind])) > self%mask_radius) then
                            self%msk(self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3))=.false.
                        end if
                    end if
                enddo
            enddo
        enddo
        allocate(self%box_scores(  self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)), source = -1.)
        allocate(self%loc_sdevs(   self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)), source = -1.)
        allocate(self%avg_int(     self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)), source = -1.)
        allocate(self%euclid_dists(self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)), source = -1.)
    end subroutine setup_iterators

    subroutine match_boxes( self, circle )
        class(nano_picker), intent(inout) :: self
        logical, optional,  intent(in)    :: circle
        type(image), allocatable :: boximgs(:)
        type(image)              :: boximg, boximg_minmax
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
        if( .not. circle_here )then
            ! use entire boxes for correlation scores
            ! iterate through positions in nanoparticle image, compare to simulated atom 
            ! construct array of boximgs
            !$ nthr = omp_get_max_threads()
            allocate(boximgs(nthr))
            do ithr = 1,nthr
                call boximgs(ithr)%new([self%boxsize,self%boxsize,self%boxsize],self%smpd)
            enddo
            !$omp parallel do schedule(static) collapse(3) default(shared) private(xoff,yoff,zoff,ithr,pos,l_err_box) proc_bind(close)
            do xoff = 1, self%nxyz_offset(1)
                do yoff = 1, self%nxyz_offset(2)
                    do zoff = 1, self%nxyz_offset(3)
                        ithr = omp_get_thread_num() + 1
                        pos  = self%positions(self%inds_offset(xoff,yoff,zoff),:)
                        call self%nano_img%window_slim( pos, self%boxsize, boximgs(ithr), outside)
                        self%box_scores(  xoff,yoff,zoff)   = self%simulated_atom%real_corr(boximgs(ithr))            ! revert to
                        self%euclid_dists(xoff,yoff,zoff)   = self%simulated_atom%euclid_dist_two_imgs(boximgs(ithr)) ! new method
                        self%loc_sdevs(   xoff,yoff,zoff)   = boximgs(ithr)%avg_loc_sdev(self%offset)
                        self%avg_int(     xoff,yoff,zoff)   = boximgs(ithr)%get_avg_int()
                    enddo 
                enddo 
            enddo
            !$omp end parallel do
            ! kill boximgs
            do ithr = 1, nthr
                call boximgs(ithr)%kill
            enddo
            deallocate(boximgs)
        else
        ! circular correlation
            maxrad    = (self%radius * 1.5) / self%smpd ! in pixels
            winsz     = ceiling(maxrad)
            npix_in   = (2 * winsz + 1)**3
            ! !$omp parallel do schedule(static) collapse(3) default(shared) private(xoff,yoff,zoff,ithr,pos,pos_center,xyz,l_err_box) proc_bind(close)
            do xoff = 1, self%nxyz_offset(1)
                do yoff = 1, self%nxyz_offset(2)
                    do zoff = 1, self%nxyz_offset(3)
                        call boximg%new([self%boxsize,self%boxsize,self%boxsize],self%smpd)
                        if (allocated(pixels1)) deallocate(pixels1)
                        if (allocated(pixels2)) deallocate(pixels2)
                        allocate(pixels1(npix_in), pixels2(npix_in), source=0.)
                        pos  = self%positions(self%inds_offset(xoff,yoff,zoff),:)
                        call self%nano_img%window_slim( pos, self%boxsize, boximg, outside)
                        call boximg_minmax%copy(boximg)
                        call boximg_minmax%norm_minmax
                        call boximg_minmax%masscen(xyz)
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
                        self%box_scores(  xoff,yoff,zoff) = pearsn_serial(pixels1(:npix_out1),pixels2(:npix_out2))      
                        self%euclid_dists(xoff,yoff,zoff) = same_energy_euclid(pixels1(:npix_out1),pixels2(:npix_out2)) 
                        self%loc_sdevs(   xoff,yoff,zoff) = boximg%avg_loc_sdev(self%offset)
                        self%avg_int(     xoff,yoff,zoff) = boximg%get_avg_int()
                        call boximg%kill
                    enddo 
                enddo 
            enddo
            ! !$omp end parallel do
            deallocate(pixels1,pixels2)
        endif
    end subroutine match_boxes

    subroutine one_box_corr(self, pos, circle, corr)
        class(nano_picker), intent(inout) :: self
        integer,            intent(in)    :: pos(3)
        logical,            intent(in)    :: circle
        real,               intent(out)   :: corr
        type(image)       :: boximg
        real              :: maxrad, xyz(3)
        real, allocatable :: pixels1(:), pixels2(:)
        integer           :: pos_center(3), winsz, npix_in, npix_out1, npix_out2
        logical           :: outside
        print *, 'inside one_box_corr'
        print *, 'pos = ', pos
        call boximg%new([self%boxsize,self%boxsize,self%boxsize],self%smpd)
        call self%nano_img%window_slim( pos, self%boxsize, boximg, outside)
        call boximg%norm_minmax
        call boximg%masscen(xyz)
        print *, 'xyz = ', xyz
        pos_center = pos + anint(xyz) + [self%boxsize/2,self%boxsize/2,self%boxsize/2]
        print *, 'pos_center = ', pos_center
        if (.not. circle) then
            corr      = self%simulated_atom%real_corr(boximg)
        else
            maxrad    = (self%radius * 1.5) / self%smpd ! in pixels
            winsz     = ceiling(maxrad)
            npix_in   = (2 * winsz + 1)**3
            allocate(pixels1(npix_in), pixels2(npix_in), source=0.)
            call self%nano_img%win2arr_rad(      pos_center(1),  pos_center(2),  pos_center(3),  winsz, npix_in, maxrad, npix_out1, pixels1)
            call self%simulated_atom%win2arr_rad(self%boxsize/2, self%boxsize/2, self%boxsize/2, winsz, npix_in, maxrad, npix_out2, pixels2)
            !corr = pearsn_serial(pixels1(:npix_out1),pixels2(:npix_out2)) ! revert to
            corr = same_energy_euclid(pixels1(:npix_out1),pixels2(:npix_out2)) ! new method
            print *, 'corr = ', corr
        end if
        call boximg%kill
    end subroutine one_box_corr

    subroutine identify_threshold( self, min_thres )
        class(nano_picker), intent(inout) :: self
        real, optional,     intent(in)    :: min_thres
        real, allocatable :: tmp(:)
        ! find peak thresholding value
        if( present(min_thres) )then
            tmp = pack(self%box_scores(:,:,:), mask=(self%box_scores(:,:,:) > min_thres)  .and. self%msk(:,:,:))
        else ! idk if this is something that makes sense to do..
            tmp = pack(self%box_scores(:,:,:), mask=(self%box_scores(:,:,:) > -1 + 1e-10) .and. self%msk(:,:,:))
        endif
        call detect_peak_thres(size(tmp), size(tmp), self%peak_thres_level, tmp, self%thres)
        print *, 'Peak threshold is ', self%thres
        deallocate(tmp)
    end subroutine identify_threshold

    subroutine identify_high_scores(self, use_zscores)
        class(nano_picker), intent(inout) :: self
        logical, optional,  intent(in)    :: use_zscores
        real                 :: Q1, mid, Q3, IQR, outlier_cutoff, temp_thres
        real, allocatable    :: pos_scores(:), lower_half_scores(:), upper_half_scores(:)
        logical              :: use_zscores_here
        temp_thres = self%thres
        if (present(use_zscores)) then
            use_zscores_here = use_zscores
        else
            use_zscores_here = .false.
        end if
        if (.not. use_zscores_here) then
            ! identify lower threshold for outlier correlation scores
            pos_scores        = pack(self%box_scores(:,:,:),   mask=self%box_scores(:,:,:) >= temp_thres .and. self%msk(:,:,:))
            mid               = median(pos_scores)
            lower_half_scores = pack(pos_scores(:), pos_scores(:) < mid)
            upper_half_scores = pack(pos_scores(:), pos_scores(:) > mid)
            Q1                = median(lower_half_scores)
            Q3                = median(upper_half_scores)
            IQR               = Q3 - Q1
            outlier_cutoff    = Q3 + 1.25*(IQR)
            self%thres        = outlier_cutoff
            ! identify euclid distance threshold (this is an upper rather than lower threshold)
            deallocate(pos_scores, lower_half_scores, upper_half_scores)
            pos_scores        = pack(self%euclid_dists(:,:,:), mask=self%box_scores(:,:,:) >= temp_thres .and. self%msk(:,:,:))
            mid               = median(pos_scores)
            lower_half_scores = pack(pos_scores(:), pos_scores(:) < mid)
            upper_half_scores = pack(pos_scores(:), pos_scores(:) > mid)
            Q1                = median(lower_half_scores)
            Q3                = median(upper_half_scores)
            IQR               = Q3 - Q1
            outlier_cutoff    = Q1 - 0.5*(IQR)
            self%euclid_thres = outlier_cutoff
            deallocate(lower_half_scores, upper_half_scores, pos_scores)
        else
            ! identify lower threshold for outlier correlation scores
            pos_scores        = pack(self%box_scores(:,:,:),   mask=self%box_scores(:,:,:) >= temp_thres .and. self%msk(:,:,:))
            mid               = median(pos_scores)
            outlier_cutoff    = 3.5 * mad_gau(pos_scores,mid) + mid
            self%thres        = outlier_cutoff
            deallocate(pos_scores)
            ! identify euclid distance threshold (this is an upper rather than lower threshold)
            pos_scores        = pack(self%euclid_dists(:,:,:), mask=self%box_scores(:,:,:) >= temp_thres .and. self%msk(:,:,:))
            mid               = median(pos_scores)
            outlier_cutoff    = mid - 3.5 * mad_gau(pos_scores,mid)
            self%euclid_thres = outlier_cutoff
            deallocate(pos_scores)
        end if
        if (self%use_euclids) self%thres = -0.999
    end subroutine identify_high_scores

    subroutine apply_threshold(self)
        class(nano_picker), intent(inout) :: self
        integer, allocatable :: pos_inds(:)
        integer              :: xoff, yoff, zoff, ipeak, npeaks
        logical              :: is_peak
        if (.not. self%use_euclids) then
            print *, 'Correlation threshold = ', self%thres
            pos_inds = pack(self%inds_offset(:,:,:), mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
            npeaks   = size(pos_inds)
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
        else
            print *, 'Euclidean distance threshold = ', self%euclid_thres
            pos_inds = pack(self%inds_offset(:,:,:), mask=self%euclid_dists(:,:,:) <= self%euclid_thres .and. self%msk(:,:,:))
            npeaks   = size(pos_inds)
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
        end if
        deallocate(pos_inds)
    end subroutine apply_threshold

    subroutine distance_filter( self, dist_thres)
        class(nano_picker), intent(inout) :: self
        real,    optional,  intent(in)    :: dist_thres
        real                 :: dist_thres_here, dist
        integer              :: nbox, ibox, jbox, xoff, yoff, zoff, npeaks, ipeak, loc
        integer, allocatable :: pos_inds(:)
        real,    allocatable :: pos_scores(:)
        logical, allocatable :: mask(:), selected_pos(:)
        logical              :: is_peak
        character(len=8)     :: crystal_system
        ! distance threshold
        if( present(dist_thres) )then
            dist_thres_here  = dist_thres
        else
            dist_thres_here  = self%offset
        endif
        if (.not. self%use_euclids) then
            pos_inds   = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
            pos_scores = pack(self%box_scores(:,:,:),   mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
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
        else
            pos_inds   = pack(self%inds_offset(:,:,:),  mask=self%euclid_dists(:,:,:) <= self%euclid_thres .and. self%msk(:,:,:))
            pos_scores = pack(self%euclid_dists(:,:,:), mask=self%euclid_dists(:,:,:) <= self%euclid_thres .and. self%msk(:,:,:))
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
                loc = minloc(pos_scores, mask=mask, dim=1)
                ! eliminate all but the best
                mask(loc) = .false.
                where( mask ) selected_pos = .false.
            enddo
        end if
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

    subroutine remove_outliers_sdev( self, ndev )
        class(nano_picker), intent(inout) :: self
        real,               intent(in)    :: ndev
        real, allocatable :: tmp(:)
        real              :: avg, sdev, t
        integer           :: npeaks, xoff, yoff, zoff
        tmp    = pack(self%loc_sdevs, mask = self%box_scores(:,:,:) >= self%thres .and. self%loc_sdevs(:,:,:) > 0. .and. self%msk(:,:,:))
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
    end subroutine remove_outliers_sdev

    ! finds distance of each atom to closest neighbor and removes atoms too far away
    subroutine remove_outliers_position(self, dist_cutoff, filename)
        class(nano_picker),         intent(inout) :: self
        real,                       intent(in)    :: dist_cutoff
        character(len=*), optional, intent(in)    :: filename
        integer, allocatable :: pos_inds(:), peak_inds(:)
        integer              :: ibox, jbox, nbox, xoff, yoff, zoff, ipeak, npeaks, index
        real,    allocatable :: dists(:)
        real                 :: pos(3), pos2(3), dist, min_dist
        logical, allocatable :: mask(:)
        logical              :: is_peak
        if (.not. allocated(self%center_positions)) THROW_HARD('self%center_positions must be allocated to call remove_outliers_positions')
        pos_inds   = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
        nbox       = size(pos_inds, dim=1)
        allocate(dists(nbox))
        allocate( mask(nbox), source = .true.)
        do ibox = 1, nbox
            pos = self%center_positions(pos_inds(ibox),:)
            min_dist = 1000
            do jbox = 1, nbox
                if (ibox .eq. jbox) cycle
                pos2 = self%center_positions(pos_inds(jbox),:)
                dist = euclid(pos,pos2)
                if (dist < min_dist) min_dist = dist
            end do
            dists(ibox) = min_dist
            if (min_dist > dist_cutoff) mask(ibox) = .false.
        end do
        ! update box scores
        peak_inds = pack(pos_inds, mask=mask)
        npeaks    = size(peak_inds)
        print *, 'NPEAKS AFTER REMOVE OUTLIERS (BY POSITION) = ', npeaks
        !$omp parallel do schedule(static) collapse(3) default(shared) private(xoff,yoff,zoff,is_peak,ipeak)
        do xoff = 1, self%nxyz_offset(1)
            do yoff = 1, self%nxyz_offset(2)
                do zoff = 1, self%nxyz_offset(3)
                    is_peak = .false.
                    do ipeak = 1,npeaks
                        if( peak_inds(ipeak) == self%inds_offset(xoff,yoff,zoff) )then
                            is_peak = .true.
                            exit
                        endif
                    enddo
                    if( .not. is_peak ) self%box_scores(xoff,yoff,zoff) = -1.
                enddo
            enddo
        enddo
        !$omp end parallel do
        if (present(filename)) then
            open(unit = 99, file=trim(filename))
            write(99, '(9(a))') 'index', ',', 'x', ',', 'y', ',', 'z', ',', 'distance'
            do ibox = 1, nbox
                index = pos_inds(ibox)
                pos   = self%center_positions(pos_inds(ibox),:)
                write(99,'(I8,4(a,f8.4))') index, ',', pos(1), ',', pos(2), ',', pos(3), ',', dists(ibox)
            end do
            close(99)
        end if
        deallocate(dists, mask, pos_inds, peak_inds)
    end subroutine remove_outliers_position

    subroutine find_centers( self )
        class(nano_picker), intent(inout) :: self
        type(image), allocatable :: atms_array(:)
        integer,     allocatable :: pos_inds(:)
        real,        allocatable :: coords(:,:)
        integer                  :: nbox, iimg, pos(3)
        logical                  :: outside
        pos_inds = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
        nbox     = size(pos_inds, dim=1)
        allocate(coords(nbox,3))
        allocate(atms_array(nbox))
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
            call self%convolved_atoms(iimg)%ifft()
            call atms_array(iimg)%ifft()
            ! want coordinates of atoms to be at the center of the images
            call self%convolved_atoms(iimg)%norm_minmax
            call self%convolved_atoms(iimg)%masscen(coords(iimg,:)) 
            coords(iimg,:) = coords(iimg,:) + real(self%convolved_atoms(iimg)%get_ldim())/2. + pos !adjust center by size and position of box
            ! update center positions for chosen boxes
            self%center_positions(pos_inds(iimg),:) = coords(iimg,:)
        enddo
        call self%simulated_atom%ifft()
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
        pos_inds = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
        nbox     = size(pos_inds, dim=1)
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
            do x = max(int_pos(1) - rad, 1), int_pos(1) + rad
                xbox = xbox + 1 ! keeps count of position in box, while iterator keeps track of position in larger image
                do y = max(int_pos(2) - rad, 1), int_pos(2) + rad
                    ybox = ybox + 1
                    do z = max(int_pos(3) - rad,1), int_pos(3) + rad
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
        print *, 'FINISHED CALC_ATOM_STATS'
    end subroutine calc_atom_stats

    subroutine calc_per_atom_corr(self)
        class(nano_picker), intent(inout) :: self
        type(nanoparticle)       :: nano
        type(image)              :: simatms
        type(parameters), target :: params
        params_glob         => params
        params_glob%element = self%element
        params_glob%smpd    = self%smpd
        call nano%new(self%raw_filename)
        if (.not. self%wrote_pdb) call self%write_pdb()
        call nano%set_atomic_coords(trim(self%pdb_filename))
        call nano%simulate_atoms(simatms=simatms)
        call nano%validate_atoms(simatms=simatms)
        call nano%kill
    end subroutine

    ! input filename with no extension
    subroutine write_pdb( self, filename )
        class(nano_picker),         intent(inout) :: self
        character(len=*), optional, intent(in)    :: filename
        integer, allocatable :: pos_inds(:)
        real,    allocatable :: coords(:,:)
        integer              :: nbox, iimg
        real                 :: pos(3)
        pos_inds = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
        nbox     = size(pos_inds, dim=1)
        allocate(coords(nbox,3))
        do iimg = 1, nbox
            pos = self%positions(pos_inds(iimg),:)
            coords(iimg,:) = self%center_positions(pos_inds(iimg),:)
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
        self%wrote_pdb = .true.
    end subroutine write_pdb

    subroutine write_boximgs( self, foldername )
        class(nano_picker),          intent(inout) :: self
        character(len=*), optional,  intent(in)    :: foldername
        integer              :: iimg, nbox
        integer, allocatable :: pos_inds(:)
        pos_inds = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
        nbox     = size(pos_inds, dim=1)
        do iimg  = 1, nbox
            if( present(foldername) )then 
                call self%convolved_atoms(iimg)%write(trim(adjustl(foldername))//'/boximg_'//trim(int2str(iimg))//'.mrc')
            else
                call self%convolved_atoms(iimg)%write('boximgs/boximg_'//trim(int2str(iimg))//'.mrc')
            endif
        enddo
        deallocate(pos_inds)
    end subroutine write_boximgs

    subroutine write_positions_and_scores( self, filename, type )
        class(nano_picker), intent(inout) :: self
        character(len=*),   intent(in)    :: filename
        character(len=*),   intent(in)    :: type
        integer,     allocatable :: pos_inds(:)
        real,        allocatable :: coords(:,:), scores(:)
        integer                  :: nbox, ipos, i, j
        pos_inds = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
        nbox     = size(pos_inds)
        allocate(coords(nbox,3))
        select case(trim(type))
            case('centers')
                scores   = pack(self%box_scores(:,:,:),   mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
                do ipos  = 1, nbox
                    coords(ipos,:) = self%center_positions(pos_inds(ipos),:) * self%smpd ! using center positions and Angstroms (same as pdb)
                enddo
            case('pixels')
                scores   = pack(self%box_scores(:,:,:),   mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
                do ipos  = 1, nbox
                    coords(ipos,:) = self%positions(pos_inds(ipos),:) ! keep in pixels
                enddo
            case('intensities')
                scores   = pack(self%avg_int(:,:,:),      mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
                do ipos  = 1, nbox
                    coords(ipos,:) = self%center_positions(pos_inds(ipos),:) * self%smpd ! using center positions and Angstroms (same as pdb)
                enddo
            case('euclid')
                scores   = pack(self%euclid_dists(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
                do ipos  = 1, nbox
                    coords(ipos,:) = self%center_positions(pos_inds(ipos),:) * self%smpd ! using center positions and Angstroms (same as pdb)
                enddo
            case DEFAULT
                THROW_HARD('write_positions_and_scores: type='//trim(type)//' is unsupported')
        end select
        open(unit=25, file=filename, status='replace', action='write')
        do i = 1, nbox
            write(25, '(I9,a)', advance='no') pos_inds(i), ','
            do j = 1, 3
                if( j /= 1 ) write(25, '(A)', advance='no') ', '
                write(25, '(F10.3)', advance='no') coords(i, j)
            enddo
            write(25, '(A)', advance='no') ', '
            write(25, '(F12.5)', advance='no') scores(i)
            write(25, *)
        enddo
        close(25)
    end subroutine write_positions_and_scores

    subroutine write_NP_image( self, sim_img_name )
        class(nano_picker), intent(inout) :: self
        character(len=*),   intent(in)    :: sim_img_name
        type(nanoparticle)       :: nano
        type(image)              :: sim_img
        type(parameters), target :: params
        params_glob => params
        params_glob%element = self%element
        params_glob%smpd    = self%smpd
        call self%write_pdb('simulate_NP')
        call nano%new(trim(self%raw_filename))
        call nano%set_atomic_coords('simulate_NP.pdb')
        call nano%simulate_atoms(simatms=sim_img)
        call sim_img%write(sim_img_name)
        call nano%kill
    end subroutine write_NP_image

    subroutine write_dist( self, csv_name, type, to_write )
        class(nano_picker), intent(inout) :: self
        character(len=*),   intent(in)    :: csv_name, type
        logical, optional,  intent(in)    :: to_write
        real,    allocatable :: pos_scores(:), lower_half_scores(:), upper_half_scores(:)
        real                 :: Q1, mid, Q3, IQR, mean
        integer              :: i
        logical              :: write_here
        if (present(to_write)) then
            write_here = to_write
        else
            write_here = .false.
        end if
        select case(trim(type))
            case('corr')
                pos_scores = pack(self%box_scores(:,:,:),   mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
                if(write_here) print *, 'SUMMARY STATISTICS OF ATOMIC CORRELATION SCORES'
            case('avg_int')
                pos_scores = pack(self%avg_int(:,:,:),      mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
                if(write_here) print *, 'SUMMARY STATISTICS OF PER-ATOM AVERAGE INTENSITY'
            case('euclid')
                pos_scores = pack(self%euclid_dists(:,:,:), mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
                if(write_here) print *, 'SUMMARY STATISTICS OF SAME ENERGY EUCLIDEAN DISTANCE'
            case DEFAULT
                THROW_HARD('write_dist: type='//trim(type)//' is unsupported')
        end select
        if (write_here) then
            mid               = median(pos_scores)
            lower_half_scores = pack(pos_scores(:), pos_scores(:) < mid)
            upper_half_scores = pack(pos_scores(:), pos_scores(:) > mid)
            Q1                = median(lower_half_scores)
            Q3                = median(upper_half_scores)
            IQR               = Q3 - Q1
            mean              = sum(pos_scores) / size(pos_scores)
            print *, 'Q1 = ', Q1
            print *, 'MEDIAN = ', mid
            print *, 'Q3 = ', Q3
            print *, 'IQR = ', IQR
            print *, 'MEAN = ', mean
            deallocate(lower_half_scores, upper_half_scores)
        end if
        open(unit=99,file=trim(csv_name))
        do i = 1, size(pos_scores)
            write(99,'(1x,f11.6)') pos_scores(i)
        enddo
        close(99)
        deallocate(pos_scores)
    end subroutine write_dist

    ! input pos in pixels, not Angstroms
    subroutine extract_img_from_pos(self, pos, img)
        class(nano_picker), intent(inout) :: self
        integer,            intent(in)    :: pos(3)
        type(image),        intent(out)   :: img
        logical :: outside
        call img%new([self%boxsize,self%boxsize,self%boxsize],self%smpd)
        call self%nano_img%window_slim( pos, self%boxsize, img, outside)
    end subroutine extract_img_from_pos

    ! input both pdbfile_* with .pdb extension
    subroutine compare_pick( self, pdbfile_ref, pdbfile_exp )
        class(nano_picker),         intent(inout) :: self
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
    
    subroutine refine_threshold( self, num_thres, ref_pdb_name, max_thres )
        class(nano_picker), intent(inout) :: self
        integer,            intent(in)    :: num_thres
        character(len=*),   intent(in)    :: ref_pdb_name
        real,    optional,  intent(in)    :: max_thres
        type(nanoparticle)       :: nano_ref, nano_exp
        type(image)              :: ref_NP,   exp_NP
        type(parameters), target :: params
        real                     :: thresholds(num_thres), thres_corrs(num_thres), max_thres_here, step, avg_corr
        integer                  :: i, optimal_index(1), num_pos
        integer, allocatable     :: pos_inds(:)
        print *, 'REFINE_THRESHOLD'
        if( present(max_thres) )then
            max_thres_here = max_thres
        else
            max_thres_here = 0.75
        end if
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
        ! call nano_ref%new(trim(self%raw_filename))
        ! call nano_ref%set_atomic_coords(trim(ref_pdb_name))
        ! call nano_ref%simulate_atoms(simatms=ref_NP)
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
            call nano_exp%new(self%raw_filename)
            call nano_exp%set_atomic_coords('sim_centers.pdb')
            call nano_exp%simulate_atoms(simatms=exp_NP)
            !thres_corrs(i) = ref_NP%real_corr(exp_NP)
            thres_corrs(i) = self%nano_img%real_corr(exp_NP) ! see if correlating with raw img helps
            !call nano_exp%validate_atoms(exp_NP, avg_corr)
            !thres_corrs(i) = avg_corr
            call nano_exp%kill
        enddo
        optimal_index = maxloc(thres_corrs)
        self%thres    = thresholds(optimal_index(1))
        call self%find_centers ! call again to set positions to the optimal
        pos_inds      = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
        num_pos       = size(pos_inds)
        print *, 'OPTIMAL THRESHOLD = ', self%thres
        print *, 'OPTIMAL CORRELATION = ', thres_corrs(optimal_index(1))
        print *, 'NUMBER POSITIONS = ', num_pos
        call nano_ref%kill
    end subroutine refine_threshold

    subroutine kill( self )
        class(nano_picker), intent(inout) :: self
        integer :: i_atom
        call self%nano_img%kill
        call self%simulated_atom%kill
        if(allocated(self%convolved_atoms))then
            do i_atom = 1, size(self%convolved_atoms)
                call self%convolved_atoms(i_atom)%kill
            enddo
            deallocate(self%convolved_atoms)
        endif
        if(allocated(self%positions)       ) deallocate(self%positions)
        if(allocated(self%center_positions)) deallocate(self%center_positions)
        if(allocated(self%inds_offset)     ) deallocate(self%inds_offset)
        if(allocated(self%box_scores)      ) deallocate(self%box_scores)
        if(allocated(self%avg_int)         ) deallocate(self%avg_int)
        if(allocated(self%euclid_dists)    ) deallocate(self%euclid_dists)
    end subroutine kill

end module simple_nano_detect_atoms