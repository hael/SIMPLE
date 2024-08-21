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
        real                     :: smpd, thres, radius, mask_radius, dist_thres, intensity_thres
        real, allocatable        :: box_scores(:,:,:), loc_sdevs(:,:,:), avg_int(:,:,:), center_positions(:,:)
        logical                  :: has_mask, wrote_pdb, circle
        logical, allocatable     :: msk(:,:,:), thres_msk(:,:,:)

    contains
        procedure :: new                        ! create new nano_picker object, initialized object variables
        procedure :: exec_nano_picker           ! execute overall picking workflow
        procedure :: simulate_atom              ! simulate atom based on element identity
        procedure :: setup_iterators            ! setup picking infrastructure (indexing, initialize arrays, etc)
        procedure :: match_boxes                ! calculate correlation between image at each position and simulated atom
        procedure :: identify_threshold         ! identify threshold for correlation scores
        procedure :: identify_high_scores       ! idenitfy outliers for correlation scores (currently retired method, kept in case useful later)
        procedure :: apply_threshold            ! remove positions with a correlation below the threshold
        procedure :: distance_filter            ! remove all positions in local neighborhoods except position with highest correlation to simualted atom
        procedure :: remove_outliers_sdev       
        procedure :: remove_outliers_position   ! remove all positions that are not within a certain distance of their nearest neighbor
        procedure :: find_centers               ! find center of mass of selected positions, save as center_positions
        procedure :: calc_atom_stats            
        procedure :: calc_per_atom_corr         ! find per-atom valid correlation by comparing selected positions with experimental map
        procedure :: refine_threshold           ! refine correlation threshold by iterating over various thresholds and optimizing correlation with exp. map
        procedure :: refine_positions           ! refine selected positions to those with highest per-atom valid correlation
        procedure :: discard_atoms              ! discard atoms with low per-atom valid correlation and / or low contact score
        procedure :: whole_map_correlation      ! calculate correlation between simulated map based on selected positions and experimental map
        ! utils
        procedure :: one_box_corr  
        procedure :: write_pdb
        procedure :: write_boximgs
        procedure :: write_positions_and_scores
        procedure :: write_NP_image
        procedure :: write_dist
        procedure :: extract_img_from_pos
        procedure :: extract_img_from_center_pos
        procedure :: kill

    end type nano_picker

    contains

    subroutine new( self, smpd, element, raw_filename, peak_thres_level, offset, dist_thres, mskdiam, intensity_level, circle)
        class(nano_picker), intent(inout) :: self
        real,               intent(in)    :: smpd
        character(len=*),   intent(in)    :: element
        character(len=100), intent(in)    :: raw_filename
        integer,            intent(in)    :: peak_thres_level, offset
        real,    optional,  intent(in)    :: dist_thres
        real,    optional,  intent(in)    :: mskdiam ! Angstroms
        integer, optional,  intent(in)    :: intensity_level
        logical, optional,  intent(in)    :: circle
        character(len=2)         :: el_ucase
        integer                  :: Z, nptcls
        logical                  :: outside
        type(image)              :: img_copy
        type(parameters), target :: params
        call self%kill
        self%smpd            = smpd
        self%element         = element
        self%raw_filename    = raw_filename
        self%offset          = offset
        self%wrote_pdb       = .false.
        self%has_mask        = .false.
        self%circle          = .false.
        self%intensity_thres = -1
        params_glob          => params
        params_glob%element  = self%element
        params_glob%smpd     = self%smpd
        ! retrieve nano_img from filename and find ldim
        el_ucase            = upperCase(params_glob%element)
        call get_element_Z_and_radius(el_ucase, Z, self%radius)
        if( Z == 0 ) THROW_HARD('Unknown element: '//el_ucase)
        ! hardcode
        self%boxsize    = round2even(self%radius / self%smpd) * 2
        self%dist_thres = self%radius / self%smpd
        call find_ldim_nptcls(self%raw_filename, self%ldim, nptcls, self%smpd)
        call self%nano_img%new(self%ldim, self%smpd)
        call self%nano_img%read(self%raw_filename)
        call self%simulated_atom%new([self%boxsize,self%boxsize,self%boxsize], self%smpd)
        self%peak_thres_level = peak_thres_level
        self%thres = -0.999 ! initial value, will be updated later
        ! create thresholding mask
        if (present(intensity_level)) then
            if (    intensity_level .eq. 1) then
                call make_intensity_mask(      self%nano_img, self%thres_msk, level=2, intensity_thres=self%intensity_thres)
            else if(intensity_level .eq. 2) then
                if (present(mskdiam)) then
                    call make_intensity_mask_2(self%nano_img, self%thres_msk, self%element, mskdiam*0.75, level=1, intensity_thres=self%intensity_thres)
                else
                    THROW_HARD('ERROR: MSKDIAM must be present to use intensity level 2')
                end if
            else
                allocate(self%thres_msk(self%ldim(1),self%ldim(2),self%ldim(3)), source=.true.)
            end if 
        else
            allocate(self%thres_msk(self%ldim(1),self%ldim(2),self%ldim(3)), source=.true.)
        end if
        print *, 'INTENSITY THRES = ', self%intensity_thres
        if( present(mskdiam) ) then
            self%has_mask = .true.
            self%mask_radius   = (mskdiam / self%smpd) / 2
            call self%nano_img%mask(self%mask_radius, 'soft')
            call self%nano_img%write('masked_img.mrc')
        end if
        if( present(circle) ) self%circle = circle
        if( present(dist_thres) ) self%dist_thres = dist_thres
    end subroutine new

    ! executes overall picking workflow
    ! some things are hard-coded for now. Can be changed to be more / less customizable
    subroutine exec_nano_picker(self, corr_thres, cs_thres)
        class(nano_picker), intent(inout) :: self
        real,  optional   , intent(in)    :: corr_thres, cs_thres
        call self%simulate_atom
        call self%setup_iterators
        call self%match_boxes
        call self%identify_threshold()
        call self%apply_threshold
        call self%distance_filter
        call self%write_pdb('before_refinement')
        print *, 'before refinement'
        call self%calc_per_atom_corr
        call self%write_NP_image('simimg3.mrc')
        call self%refine_threshold(20, max_thres=0.7)
        call self%write_pdb('after_refinement')
        print *, 'after refinement'
        call self%calc_per_atom_corr
        if (.not. present(corr_thres) .and. .not. present(cs_thres))  then
            call self%discard_atoms(use_valid_corr=.true., use_cs_thres=.true.)
        else if (present(corr_thres) .and. .not. present(cs_thres))   then
            call self%discard_atoms(use_valid_corr=.true., use_cs_thres=.true., corr_thres=corr_thres)
        else if (present(cs_thres)   .and. .not. present(corr_thres)) then
            call self%discard_atoms(use_valid_corr=.true., use_cs_thres=.true., cs_thres=cs_thres)
        else 
            call self%discard_atoms(use_valid_corr=.true., use_cs_thres=.true., corr_thres=corr_thres, cs_thres=cs_thres)
        end if
        call self%write_pdb('after_discard_atoms')
        print *, 'after discard atoms'
        call self%calc_per_atom_corr
    end subroutine exec_nano_picker

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
        call atom%set_element(1,self%element)
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
        do xind = 1, nxyz(1), self%offset
            self%nxyz_offset(1) = self%nxyz_offset(1) + 1
            self%nxyz_offset(2) = 0
            do yind = 1, nxyz(2), self%offset
                self%nxyz_offset(2) = self%nxyz_offset(2) + 1
                self%nxyz_offset(3) = 0
                do zind = 1, nxyz(3), self%offset
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
        do xind = 1, nxyz(1), self%offset
            self%nxyz_offset(1) = self%nxyz_offset(1) + 1
            self%nxyz_offset(2) = 0
            do yind = 1, nxyz(2), self%offset
                self%nxyz_offset(2) = self%nxyz_offset(2) + 1
                self%nxyz_offset(3) = 0
                do zind = 1, nxyz(3), self%offset
                    self%nxyz_offset(3)      = self%nxyz_offset(3) + 1
                    nboxes                   = nboxes + 1
                    self%positions(nboxes,:) = [xind,yind,zind]
                    self%inds_offset(self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)) = nboxes
                    if(self%has_mask) then
                        if(euclid(real(self%ldim / 2), real([xind,yind,zind])) > self%mask_radius) then
                            self%msk(self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3))=.false.
                        end if
                    end if
                    if(self%msk(self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3))) then
                        self%msk(self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)) = logical_box(self%thres_msk, pos=self%positions(nboxes,:), boxsize=self%boxsize, num_px=ceiling(self%radius / self%smpd))
                    end if
                enddo
            enddo
        enddo
        allocate(self%box_scores(  self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)), source = -1.)
        allocate(self%loc_sdevs(   self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)), source = -1.)
        allocate(self%avg_int(     self%nxyz_offset(1),self%nxyz_offset(2),self%nxyz_offset(3)), source = -1.)
    end subroutine setup_iterators

    subroutine match_boxes( self )
        class(nano_picker), intent(inout) :: self
        type(image), allocatable :: boximgs(:)
        type(image)              :: boximg, boximg_minmax
        integer                  :: xoff, yoff, zoff, pos(3), pos_center(3), ithr, nthr, winsz, npix_in, npix_out1, npix_out2
        logical                  :: l_err_box, circle_here, outside
        real                     :: maxrad, xyz(3)
        real,        allocatable :: pixels1(:), pixels2(:)
        if( .not. self%circle )then
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

    subroutine one_box_corr(self, pos, circle, corr, center)
        class(nano_picker), intent(inout) :: self
        integer,            intent(in)    :: pos(3)
        logical,            intent(in)    :: circle
        real,               intent(out)   :: corr
        real, optional,     intent(out)   :: center(3)
        type(image)          :: boximg
        real                 :: maxrad, xyz(3)
        real, allocatable    :: pixels1(:), pixels2(:), boximg_rmat(:,:,:)
        integer              :: pos_center(3), winsz, npix_in, npix_out1, npix_out2
        logical              :: outside
        logical, allocatable :: boximg_mask(:,:,:)
        call boximg%new([self%boxsize,self%boxsize,self%boxsize],self%smpd)
        call self%nano_img%window_slim( pos, self%boxsize, boximg, outside)
        allocate(boximg_rmat(self%boxsize,self%boxsize,self%boxsize), source=0.)
        allocate(boximg_mask(self%boxsize,self%boxsize,self%boxsize), source=.true.)
        boximg_rmat = boximg%get_rmat()
        where (boximg_rmat < self%intensity_thres)
            boximg_mask = .false.
        end where
        call boximg%masscen_adjusted(xyz, mask_in=boximg_mask)
        pos_center = pos + anint(xyz)
        if (present(center)) then
            center = pos + xyz
        end if
        if (.not. circle) then
            corr      = self%simulated_atom%real_corr(boximg)
        else
            maxrad    = (self%radius * 1.5) / self%smpd ! in pixels
            winsz     = ceiling(maxrad)
            npix_in   = (2 * winsz + 1)**3
            allocate(pixels1(npix_in), pixels2(npix_in), source=0.)
            call self%nano_img%win2arr_rad(      pos_center(1),  pos_center(2),  pos_center(3),  winsz, npix_in, maxrad, npix_out1, pixels1)
            call self%simulated_atom%win2arr_rad(self%boxsize/2, self%boxsize/2, self%boxsize/2, winsz, npix_in, maxrad, npix_out2, pixels2)
            corr = pearsn_serial(pixels1(:npix_out1),pixels2(:npix_out2))
        end if
        call boximg%kill
        deallocate(boximg_rmat, boximg_mask)
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
            deallocate(pos_scores)
        else
            ! identify lower threshold for outlier correlation scores
            pos_scores        = pack(self%box_scores(:,:,:),   mask=self%box_scores(:,:,:) >= temp_thres .and. self%msk(:,:,:))
            mid               = median(pos_scores)
            outlier_cutoff    = 3.5 * mad_gau(pos_scores,mid) + mid
            self%thres        = outlier_cutoff
            deallocate(pos_scores)
        end if
    end subroutine identify_high_scores

    subroutine apply_threshold(self)
        class(nano_picker), intent(inout) :: self
        integer, allocatable :: pos_inds(:)
        integer              :: xoff, yoff, zoff, ipeak, npeaks
        logical              :: is_peak
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
        deallocate(pos_inds)
    end subroutine apply_threshold

    subroutine distance_filter( self )
        class(nano_picker), intent(inout) :: self
        real                 :: dist
        integer              :: nbox, ibox, jbox, xoff, yoff, zoff, npeaks, ipeak, loc
        integer, allocatable :: pos_inds(:)
        real,    allocatable :: pos_scores(:)
        logical, allocatable :: mask(:), selected_pos(:)
        logical              :: is_peak
        character(len=8)     :: crystal_system
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
                if( dist <= self%dist_thres ) mask(jbox) = .true.
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
        call self%find_centers
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
        real,        allocatable :: coords(:,:), boximg_rmat(:,:,:)
        integer                  :: nbox, iimg, pos(3)
        logical                  :: outside
        logical,     allocatable :: boximg_mask(:,:,:)
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
            self%convolved_atoms(iimg) = atms_array(iimg)
            ! want coordinates of atoms to be at the center of the images
            allocate(boximg_rmat(self%boxsize,self%boxsize,self%boxsize), source=0.)
            allocate(boximg_mask(self%boxsize,self%boxsize,self%boxsize), source=.true.)
            boximg_rmat = self%convolved_atoms(iimg)%get_rmat()
            where (boximg_rmat < self%intensity_thres)
                boximg_mask = .false.
            end where
            call self%convolved_atoms(iimg)%masscen_adjusted(coords(iimg,:),boximg_mask) 
            deallocate(boximg_rmat, boximg_mask)
            coords(iimg,:) = coords(iimg,:) + pos
            ! update center positions for chosen boxes
            self%center_positions(pos_inds(iimg),:) = coords(iimg,:)
        enddo
        call self%simulated_atom%ifft()
        self%wrote_pdb = .false. ! when centers are updated, we will want to update pdb file
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
        if (self%has_mask) then
            call nano%new(self%raw_filename, msk=self%mask_radius)
        else
            call nano%new(self%raw_filename)
        end if
        if (.not. self%wrote_pdb) call self%write_pdb()
        call nano%set_atomic_coords(trim(self%pdb_filename))
        call nano%simulate_atoms(simatms=simatms)
        call nano%validate_atoms(simatms=simatms)
        !call nano%write_centers('coordinates_in_nanoparticle_angstroms',which='valid_corr')
        call nano%kill
    end subroutine calc_per_atom_corr

    subroutine refine_threshold( self, num_thres, max_thres )
        class(nano_picker), intent(inout) :: self
        integer,            intent(in)    :: num_thres
        real,    optional,  intent(in)    :: max_thres
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
        ! iterate through thresholds and calculate correlation to experimental map
        print *, 'ITERATIONS BEGIN'
        do i = 1, num_thres
            self%thres     = thresholds(i) ! need to set self%thres because it is called in multiple subroutines
            print *, 'i = ', i, 'self%thres = ', self%thres
            call self%find_centers
            thres_corrs(i) = self%whole_map_correlation()
            print *, 'CORR = ', thres_corrs(i)
        enddo
        print *, 'ITERATIONS END'
        optimal_index = maxloc(thres_corrs)
        self%thres    = thresholds(optimal_index(1))
        call self%find_centers ! call again to set positions to the optimal
        pos_inds      = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
        num_pos       = size(pos_inds)
        print *, 'OPTIMAL THRESHOLD = ', self%thres
        print *, 'OPTIMAL CORRELATION = ', thres_corrs(optimal_index(1))
        print *, 'NUMBER POSITIONS = ', num_pos
    end subroutine refine_threshold

    subroutine refine_positions(self)
        class(nano_picker), intent(inout) :: self
        integer, allocatable :: pos_inds(:)
        integer              :: npos, ipos, pos(3), x, y, z, xoff, yoff, zoff, optimal_pos(3)
        real                 :: corr, optimal_corr, xyz(3), center(3), optimal_center(3)
        pos_inds   = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
        npos       = size(pos_inds)
        do ipos = 1, npos
            pos            = self%positions(pos_inds(ipos),:)
            optimal_pos    = pos
            call self%one_box_corr(pos, circle=.true., corr=optimal_corr)
            optimal_center = self%center_positions(pos_inds(ipos),:)
            do x = pos(1) - self%offset, pos(1) + self%offset
                do y = pos(2) - self%offset, pos(2) + self%offset
                    do z = pos(3) - self%offset, pos(3) + self%offset
                        call self%one_box_corr([x,y,z], circle=.true., corr=corr, center=center)      
                        if (corr > optimal_corr) then
                            optimal_corr   = corr
                            optimal_pos    = [x,y,z]
                            optimal_center = center
                        end if
                    end do
                end do
            end do
            !$omp parallel do schedule(static) collapse(3) default(shared) private(xoff,yoff,zoff,ipos)
            do xoff = 1, self%nxyz_offset(1)
                do yoff = 1, self%nxyz_offset(2)
                    do zoff = 1, self%nxyz_offset(3)
                        if(pos_inds(ipos) .eq. self%inds_offset(xoff,yoff,zoff)) then
                            self%center_positions(pos_inds(ipos),:) = optimal_center
                            self%positions(pos_inds(ipos),:) = optimal_pos
                        end if
                    enddo
                enddo
            enddo
            !$omp end parallel do
        end do
        deallocate(pos_inds)
        call self%find_centers
    end subroutine refine_positions

    subroutine discard_atoms(self, use_cs_thres, use_valid_corr, corr_thres, cs_thres)
        class(nano_picker), intent(inout) :: self
        logical, optional,  intent(in)    :: use_cs_thres, use_valid_corr
        real,    optional,  intent(in)    :: corr_thres, cs_thres
        logical                  :: cs_here, corr_here ! true by default, otherwise take input vals
        logical                  :: is_peak
        logical, allocatable     :: selected_pos(:)
        real                     :: corr_thres_here, corr_thres_sigma, cs_thres_here, pos(3), maxrad
        real,    allocatable     :: corrs(:), positions(:,:), coords(:,:)
        integer                  :: ibox, nbox, xoff, yoff, zoff, ipeak, npeaks
        integer, allocatable     :: pos_inds(:), contact_scores(:)
        type(nanoparticle)       :: nano 
        type(image)              :: simatms
        type(parameters), target :: params
        type(stats_struct)       :: cscore_stats
        params_glob         => params
        params_glob%element = self%element
        params_glob%smpd    = self%smpd
        if (present(use_cs_thres)) then
            cs_here = use_cs_thres
        else
            cs_here = .true.
        end if
        if (present(use_valid_corr)) then
            corr_here = use_valid_corr
        else
            corr_here = .true.
        end if
        ! filter out lowly correlated atoms
        if (corr_here) then
            call nano%new(trim(self%raw_filename))
            if (.not. self%wrote_pdb) call self%write_pdb()
            call nano%set_atomic_coords(trim(self%pdb_filename))
            call nano%simulate_atoms(simatms=simatms)
            call nano%kill
            pos_inds = pack(self%inds_offset(:,:,:), mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
            nbox     = size(pos_inds)
            allocate(selected_pos(nbox), source=.true.)
            allocate(corrs(nbox))
            maxrad  = (self%radius * 1.5) / self%smpd ! in pixels
            open(unit=44,file='discard_coords.csv')
            do ibox = 1, nbox
                pos         = self%center_positions(pos_inds(ibox),:)
                corrs(ibox) = one_atom_valid_corr(pos, maxrad, self%nano_img, simatms)
                write(44,'(1x,4(f8.3,a))') pos(1), ',', pos(2), ',', pos(3), ',', corrs(ibox)
            end do
            close(44)
            if (present(corr_thres)) then
                corr_thres_here = corr_thres
            else
                corr_thres_sigma = -2.0 ! took from nanoparticle class paramter
                corr_thres_here = min(max(robust_sigma_thres(corrs(:), corr_thres_sigma), 0.3), 0.7)
            end if
            print *, 'Valid_corr threshold: ', corr_thres_here
            do ibox = 1, nbox
                if (corrs(ibox) < corr_thres_here) selected_pos(ibox) = .false.
            end do
            pos_inds = pack(pos_inds, mask=selected_pos)
            npeaks   = size(pos_inds)
            print *, 'Number positions removed due to low valid-corr: ', nbox - npeaks
            print *, 'NPEAKS AFTER LOW VALID CORR REMOVED = ', npeaks
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
            deallocate(pos_inds, selected_pos)
            call nano%kill
            self%wrote_pdb = .false.
            deallocate(corrs)
        end if
        ! remove atoms with low contact scores (cs)
        if(cs_here) then
            pos_inds = pack(self%inds_offset(:,:,:), mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
            nbox     = size(pos_inds)
            allocate(positions(3,nbox))
            allocate(contact_scores(nbox))
            do ibox = 1, nbox
                positions(:,ibox) = (self%center_positions(pos_inds(ibox),:) - 1.) * self%smpd
            end do
            call calc_contact_scores(self%element, positions, contact_scores)
            call calc_stats(real(contact_scores), cscore_stats)
            if( present(cs_thres) )then
                cs_thres_here = cs_thres
            else
                cs_thres_here = min(5,max(3,nint(cscore_stats%avg - cscore_stats%sdev)))
                print *, 'Contact score threshold: ', cs_thres_here
            endif
            allocate(selected_pos(nbox), source=.true.)
            do ibox = 1, nbox
                if (contact_scores(ibox) < cs_thres_here) then
                    selected_pos(ibox) = .false.
                end if
            end do
            pos_inds = pack(pos_inds, mask=selected_pos)
            npeaks   = size(pos_inds)
            print *, 'Number positions removed due to low contact score: ', nbox - npeaks
            print *, 'NPEAKS AFTER LOW CONTACT SCORES REMOVED = ', npeaks
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
            deallocate(pos_inds, selected_pos)
        end if
    end subroutine discard_atoms

    function whole_map_correlation(self, compare_to_pdb) result(corr)
        class(nano_picker),         intent(inout) :: self
        character(len=*), optional, intent(in)    :: compare_to_pdb
        real                     :: corr
        type(nanoparticle)       :: nano, nano_ref
        type(image)              :: simatms, simatms_ref
        type(parameters), target :: params
        if (.not. self%wrote_pdb) call self%write_pdb()
        params_glob => params
        params_glob%element = self%element
        params_glob%smpd    = self%smpd
        call nano%new(self%raw_filename)
        call nano%set_atomic_coords(trim(self%pdb_filename))
        call nano%simulate_atoms(simatms)
        if (present(compare_to_pdb)) then
            call nano_ref%new(self%raw_filename)
            call nano_ref%set_atomic_coords(trim(compare_to_pdb))
            call nano_ref%simulate_atoms(simatms_ref)
            corr = simatms_ref%real_corr(simatms)
            call nano_ref%kill
            call simatms_ref%kill
        else
            corr = self%nano_img%real_corr(simatms)
        end if
        call nano%kill
        call simatms%kill
    end function whole_map_correlation

    ! input filename with no extension
    subroutine write_pdb( self, filename )
        class(nano_picker),         intent(inout) :: self
        character(len=*), optional, intent(in)    :: filename
        integer, allocatable     :: pos_inds(:)
        real,    allocatable     :: coords(:,:), pos_scores(:)
        integer                  :: nbox, iimg
        real                     :: pos(3)
        type(parameters), target :: params
        params_glob         => params
        params_glob%element = self%element
        params_glob%smpd    = self%smpd
        pos_inds   = pack(self%inds_offset(:,:,:),  mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
        pos_scores = pack(self%box_scores(:,:,:),   mask=self%box_scores(:,:,:) >= self%thres .and. self%msk(:,:,:))
        nbox       = size(pos_inds, dim=1)
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
                    coords(ipos,:) = ( self%center_positions(pos_inds(ipos),:) - 1 ) * self%smpd ! using center positions and Angstroms (same as pdb)
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
        if (.not. self%wrote_pdb) call self%write_pdb('simulate_NP')
        call nano%new(trim(self%raw_filename))
        call nano%set_atomic_coords(trim(self%pdb_filename))
        call nano%simulate_atoms(simatms=sim_img)
        call sim_img%write(sim_img_name)
        call nano%kill
        call sim_img%kill
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

    ! input center_pos in pixels, not Angstroms
    subroutine extract_img_from_center_pos(self, center_pos, img)
        class(nano_picker), intent(inout) :: self
        real,               intent(in)    :: center_pos(3)
        type(image),        intent(out)   :: img
        real    :: maxrad
        integer :: winsz, edge_pos(3)
        logical :: outside
        maxrad    = (self%radius * 1.5) / self%smpd ! in pixels
        winsz     = ceiling(maxrad)
        call img%new([winsz,winsz,winsz],self%smpd)
        edge_pos  = anint(center_pos - [winsz/2., winsz/2., winsz/2.])
        call self%nano_img%window_slim( edge_pos, winsz, img, outside)
        call img%mask(maxrad, 'hard')
    end subroutine extract_img_from_center_pos

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
        if(allocated(self%thres_msk)       ) deallocate(self%thres_msk)
    end subroutine kill

end module simple_nano_detect_atoms