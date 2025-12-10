module simple_pickref
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
use simple_segmentation
implicit none

public :: read_mic_raw_pickref, pickref, multiref_merge
private
#include "simple_local_flags.inc"

! class constants
real,    parameter :: NDEV_DEFAULT = 2.5
real,    parameter :: MSKDIAM2LP = 0.15, LP_LB = 30., LP_UB = 15.
integer, parameter :: OFFSET_DEFAULT = 3
logical, parameter :: L_WRITE  = .false.
logical, parameter :: L_DEBUG  = .false.

! class variables
integer     :: ldim_raw(3)
real        :: smpd_raw
type(image) :: mic_raw

! instance
type pickref
    private
    ! control parameters
    real                     :: smpd_shrink = 0., maxdiam = 0., sxx = 0., t = 0., ndev = 0.
    real                     :: dist_thres  = 0.
    integer                  :: ldim(3), ldim_box(3), ldim_raw_box(3), nboxes = 0, nboxes_ub = 0, nx = 0, ny = 0
    integer                  :: nx_offset = 0, ny_offset = 0, npeaks = 0, nrefs = 0, offset = 0, offset_ub = 0
    integer                  :: peak_thres_level = 2, nboxes_max = 0
    ! peak stats
    real                     :: smd_corr = 0., ksstat_corr = 0., prob_corr = 0.
    real                     :: a_corr_peak = 0., s_corr_peak = 0., a_corr_nonpeak = 0., s_corr_nonpeak = 0.
    ! images
    type(image)              :: mic_shrink, mic_roi, mic_copy
    type(image), allocatable :: boxrefs(:)
    ! masks
    logical,     allocatable :: l_mic_mask(:,:), l_err_refs(:)
    ! indices
    integer,     allocatable :: positions(:,:), inds_offset(:,:)
    ! scores & local standard deviations
    real,        allocatable :: box_scores(:,:), box_scores_mem(:,:), loc_sdevs(:,:), loc_sdevs_mem(:,:)
    ! flags
    logical                  :: l_black_ptcls_input = .true.
    logical                  :: l_roi               = .false.
    logical                  :: exists              = .false.
contains
    procedure          :: new
    procedure          :: refpick
    procedure, private :: setup_iterators
    procedure, private :: match_boximgs
    procedure, private :: detect_peaks
    procedure, private :: distance_filter
    procedure, private :: remove_outliers
    procedure, private :: peak_vs_nonpeak_stats
    procedure          :: get_positions
    procedure          :: get_loc_sdevs
    procedure          :: get_scores
    procedure          :: get_nboxes
    procedure          :: get_box
    procedure          :: get_maxdiam
    procedure, private :: write_boxfile
    procedure          :: report_thumb_den
    procedure          :: report_boxfile
    procedure, private :: refine_upscaled
    procedure          :: kill
end type

contains

    subroutine read_mic_raw_pickref( micname, smpd, subtr_backgr )
        class(string),     intent(in) :: micname !< micrograph file name
        real,              intent(in) :: smpd    !< sampling distance in A
        logical, optional, intent(in) :: subtr_backgr
        integer :: nframes
        ! set micrograph info
        call find_ldim_nptcls(micname, ldim_raw, nframes)
        if( ldim_raw(3) /= 1 .or. nframes /= 1 ) THROW_HARD('Only for 2D images')
        smpd_raw = smpd
        ! read micrograph
        call mic_raw%new(ldim_raw, smpd_raw)
        call mic_raw%read(micname)
        if( present(subtr_backgr) )then
            if( subtr_backgr ) call mic_raw%subtract_background(HP_BACKGR_SUBTR)
        endif
    end subroutine read_mic_raw_pickref

    subroutine refpick( self, self_refine )
        use simple_micproc, only: flag_amorphous_carbon, flag_ice
        class(pickref),           intent(inout) :: self
        class(pickref), optional, intent(inout) :: self_refine
        integer, allocatable :: pos(:,:)
        if( self%l_roi )then
            call flag_ice(mic_raw, self%l_mic_mask)
            call flag_amorphous_carbon(self%mic_roi, self%l_mic_mask)
            if( count(self%l_mic_mask) < nint(0.01*product(self%ldim)) )then
                self%npeaks = 0
                if( present(self_refine) ) self_refine%npeaks = 0
                return
            endif
        endif
        ! establish iterators
        call self%setup_iterators
        ! matching
        call self%match_boximgs
        call self%detect_peaks
        if( self%npeaks == 0 )then
            if( present(self_refine) ) self_refine%npeaks = 0
            return
        endif
        ! disabled for now. This needs to be done in a more global fashion, not per-micrograph
        ! call self%remove_outliers
        call self%peak_vs_nonpeak_stats
        if( present(self_refine) )then
            call self%get_positions(pos, self_refine%smpd_shrink)
            call self_refine%setup_iterators
            call self_refine%refine_upscaled(pos, self%smpd_shrink, self%offset)
            call self_refine%distance_filter
            deallocate(pos)
        endif
    end subroutine refpick

    subroutine new( self, pcontrast, smpd_shrink, imgs, offset, ndev, roi, nboxes_max )
        class(pickref),    intent(inout) :: self
        character(len=*),  intent(in)    :: pcontrast
        real,              intent(in)    :: smpd_shrink
        class(image),      intent(inout) :: imgs(:)
        integer, optional, intent(in)    :: offset
        real,    optional, intent(in)    :: ndev    !< # std devs for outlier detection
        logical, optional, intent(in)    :: roi
        integer, optional, intent(in)    :: nboxes_max
        real    :: scale, lp
        integer :: iref, box
        if( self%exists ) call self%kill
        self%l_roi = .false.
        if( present(roi) ) self%l_roi = roi
        self%ldim_raw_box = imgs(1)%get_ldim()
        if( self%ldim_raw_box(3) /= 1 ) THROW_HARD('box references must be 2D')
        box = self%ldim_raw_box(1)
        ! pixel size after downscaling
        self%smpd_shrink = max(smpd_raw,smpd_shrink)
        ! set shrunken logical dimensions
        scale             = smpd_raw / self%smpd_shrink
        self%ldim(1)      = round2even(real(ldim_raw(1)) * scale)
        self%ldim(2)      = round2even(real(ldim_raw(2)) * scale)
        self%ldim(3)      = 1  
        self%maxdiam      = smpd_raw * real(box)
        self%ldim_box(1)  = round2even(real(box) * scale)
        self%ldim_box(2)  = self%ldim_box(1)
        self%ldim_box(3)  = 1
        self%nx           = self%ldim(1) - box
        self%ny           = self%ldim(2) - box
        self%dist_thres   = (self%maxdiam/3.)/self%smpd_shrink
        ! take care of matching references
        self%nrefs = size(imgs)
        ! set # pixels in x/y for both box sizes
        self%nx = self%ldim(1) - self%ldim_box(1)
        self%ny = self%ldim(2) - self%ldim_box(2)
        allocate(self%l_err_refs(self%nrefs), self%boxrefs(self%nrefs))
        !$omp parallel do schedule(static) default(shared) private(iref) proc_bind(close)
        do iref = 1,self%nrefs
            call imgs(iref)%fft
            call self%boxrefs(iref)%new(self%ldim_box, self%smpd_shrink, wthreads=.false.)
            call self%boxrefs(iref)%set_ft(.true.)
            call imgs(iref)%clip(self%boxrefs(iref))
            call imgs(iref)%ifft
            call self%boxrefs(iref)%ifft
            call self%boxrefs(iref)%prenorm4real_corr(self%l_err_refs(iref))
        end do
        !$omp end parallel do
        ! set peak threshold level
        select case(trim(params_glob%particle_density))
            case('low')
                self%peak_thres_level = 1
            case('optimal')
                self%peak_thres_level = 2
            case('high')
                self%peak_thres_level = 3
            case DEFAULT
                THROW_HARD('Unsupported particle density level')
        end select
        ! set ndev
        self%ndev = NDEV_DEFAULT
        if( present(ndev) ) self%ndev = ndev
        ! set offset
        self%offset = OFFSET_DEFAULT
        if( present(offset) ) self%offset = offset
        self%offset_ub = self%offset * 2
        ! shrink micrograph
        call self%mic_shrink%new(self%ldim, self%smpd_shrink)
        call self%mic_shrink%set_ft(.true.)
        call mic_raw%mul(real(product(ldim_raw))) ! to prevent numerical underflow when performing FFT
        call mic_raw%fft
        call mic_raw%clip(self%mic_shrink)
        self%l_black_ptcls_input = .false.
        select case(trim(pcontrast))
            case('black')
                ! flip contrast (assuming black particle contrast on input)
                call self%mic_shrink%mul(-1.)
                self%l_black_ptcls_input = .true.
            case('white')
                ! nothing to do
            case DEFAULT
                THROW_HARD('uknown pcontrast parameter, use (black|white)')
        end select
        ! back to real-space
        call mic_raw%ifft
        call mic_raw%div(real(product(ldim_raw)))
        call self%mic_shrink%ifft
        if( self%l_roi )then
            ! backup of shrunken micrograph
            call self%mic_roi%copy(self%mic_shrink)
        endif
        ! low-pass filter mic_shrink
        lp = min(max(LP_LB,MSKDIAM2LP * self%maxdiam),LP_UB)
        call self%mic_copy%copy(self%mic_shrink)
        call self%mic_shrink%bp(0.,lp)
        if( allocated(self%l_mic_mask) ) deallocate(self%l_mic_mask)
        allocate(self%l_mic_mask(self%ldim(1),self%ldim(2)), source=.true.)
        if( present(nboxes_max) ) self%nboxes_max = nboxes_max
        ! flag existence
        self%exists = .true.
    end subroutine new

    subroutine setup_iterators( self )
        class(pickref), intent(inout) :: self
        integer :: xind, yind
        self%nboxes    = 0
        self%nx_offset = 0
        do xind = 0,self%nx,self%offset
            self%nx_offset = self%nx_offset + 1
            self%ny_offset = 0
            do yind = 0,self%ny,self%offset
                self%nboxes    = self%nboxes    + 1
                self%ny_offset = self%ny_offset + 1
                if( self%l_mic_mask(xind+1,yind+1) ) self%nboxes = self%nboxes + 1
            end do
        end do
        ! count # boxes, upper bound
        self%nboxes_ub = 0
        do xind = 0,self%nx,self%offset_ub
            do yind = 0,self%ny,self%offset_ub
                if( self%l_mic_mask(xind+1,yind+1) ) self%nboxes_ub = self%nboxes_ub + 1
            end do
        end do
        ! allocate and set positions and inds_offset
        if( allocated(self%positions)   ) deallocate(self%positions)
        if( allocated(self%inds_offset) ) deallocate(self%inds_offset)
        allocate(self%positions(self%nboxes,2), self%inds_offset(self%nx_offset,self%ny_offset), source=0)
        ! calculate total # boxes
        self%nboxes    = 0
        self%nx_offset = 0
        do xind = 0,self%nx,self%offset
            self%nx_offset = self%nx_offset + 1
            self%ny_offset = 0
            do yind = 0,self%ny,self%offset
                self%ny_offset = self%ny_offset + 1
                if( self%l_mic_mask(xind+1,yind+1) )then
                    self%nboxes = self%nboxes + 1
                    self%positions(self%nboxes,:) = [xind,yind]
                    self%inds_offset(self%nx_offset,self%ny_offset) = self%nboxes
                endif
            end do
        end do
        ! allocate box_scores
        if( allocated(self%box_scores) ) deallocate(self%box_scores)
        allocate(self%box_scores(self%nx_offset,self%ny_offset), source = -1.)
        ! allocate loc_sdevs
        if( allocated(self%loc_sdevs) ) deallocate(self%loc_sdevs)
        allocate(self%loc_sdevs(self%nx_offset,self%ny_offset),  source = -1.)
    end subroutine setup_iterators

    subroutine match_boximgs( self )
        class(pickref), intent(inout) :: self
        logical     :: outside, l_err
        integer     :: pos(2), ioff, joff, ithr, iref
        real        :: scores(self%nrefs)
        type(image) :: boximgs_heap(nthr_glob)
        ! construct heap
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box, self%smpd_shrink)
        end do
        !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr,outside,iref,scores,pos,l_err) proc_bind(close)
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                if( self%inds_offset(ioff,joff) == 0 )then
                    self%box_scores(ioff,joff) = -1.
                    self%loc_sdevs(ioff,joff)  = -1.
                else
                    pos = self%positions(self%inds_offset(ioff,joff),:)
                    if( self%l_mic_mask(pos(1)+1,pos(2)+1) )then
                        ithr = omp_get_thread_num() + 1
                        call self%mic_shrink%window_slim(pos, self%ldim_box(1), boximgs_heap(ithr), outside)
                        self%loc_sdevs(ioff,joff) = boximgs_heap(ithr)%avg_loc_sdev(self%offset)
                        call boximgs_heap(ithr)%prenorm4real_corr(l_err)
                        if( l_err )then
                            self%box_scores(ioff,joff) = -1.
                        else
                            do iref = 1,self%nrefs
                                if( self%l_err_refs(iref) )then
                                    scores(iref) = -1.
                                else
                                    scores(iref) = self%boxrefs(iref)%real_corr_prenorm(boximgs_heap(ithr))
                                endif
                            end do
                            self%box_scores(ioff,joff) = maxval(scores)
                        endif
                    else
                        self%box_scores(ioff,joff) = -1.
                        self%loc_sdevs(ioff,joff)  = -1.
                    endif
                endif
            end do
        end do
        !$omp end parallel do
        ! save a box_scores copy in memory
        if( allocated(self%box_scores_mem) ) deallocate(self%box_scores_mem)
        allocate(self%box_scores_mem(self%nx_offset,self%ny_offset), source=self%box_scores)
        ! save a loc_sdevs copy in memory
        if( allocated(self%loc_sdevs_mem) ) deallocate(self%loc_sdevs_mem)
        allocate(self%loc_sdevs_mem(self%nx_offset,self%ny_offset), source=self%loc_sdevs)
        ! destruct heap
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
    end subroutine match_boximgs

    subroutine detect_peaks( self )
        class(pickref), intent(inout) :: self
        real, allocatable :: tmp(:)
        integer :: n
        n = count(self%box_scores > 0.)
        if( n == 0 )then
            self%npeaks = 0
            return
        endif
        ! select peak targets > 0
        tmp = pack(self%box_scores, mask=(self%box_scores > 0.))
        ! upper bound thresholding
        call detect_peak_thres(n, self%nboxes_ub, 0, tmp, self%t)
        ! modify box_scores accordingly
        where( self%box_scores >= self%t )
            ! there's a peak
        elsewhere
            self%box_scores = -1.
        end where
        ! apply distance filter
        call self%distance_filter
        ! peak detection
        n = count(self%box_scores > 0.)
        if( n == 0 )then
            self%npeaks = 0
            return
        endif
        deallocate(tmp)
        tmp = pack(self%box_scores, mask=(self%box_scores > 0.))
        ! Otsu thres (level=1)
        call detect_peak_thres(n, 1, tmp, self%t)
        call refine_peak_thres_sortmeans(n, self%peak_thres_level, tmp, self%t)
        ! call detect_peak_thres_sortmeans(n, self%peak_thres_level, tmp, self%t)
        where( self%box_scores >= self%t )
            ! there's a peak
        elsewhere
            self%box_scores = -1.
        end where
        self%npeaks = count(self%box_scores >= self%t)
        write(logfhandle,'(a,1x,f5.2)') 'peak threshold identified:             ', self%t
        write(logfhandle,'(a,1x,I5)'  ) '# peaks detected:                      ', self%npeaks
        deallocate(tmp)
    end subroutine detect_peaks

    subroutine distance_filter( self )
        class(pickref), intent(inout) :: self
        integer, allocatable :: pos_inds(:), order(:)
        real,    allocatable :: pos_scores(:), tmp(:)
        logical, allocatable :: mask(:), selected_pos(:)
        integer :: nbox, npeaks, ibox, jbox, loc, ioff, joff, ithr, ipeak, i
        real    :: dist
        logical :: is_corr_peak
        pos_inds   = pack(self%inds_offset(:,:), mask=self%box_scores(:,:) >= self%t)
        pos_scores = pack(self%box_scores(:,:),  mask=self%box_scores(:,:) >= self%t)
        nbox       = size(pos_inds)
        allocate(mask(nbox),         source=.false.)
        allocate(selected_pos(nbox), source=.true. )
        ! This is a greedy algorithm for distance thresholding. Hence, the order in which the peaks are considered
        ! matters. Create an order from highest to lowest peak and traverse them in that order. In this way, the highest
        ! peaks are selected and their neighbors are stripped off first
        order = (/(i,i=1,nbox)/)
        tmp   = pack(self%box_scores, mask=self%box_scores(:,:) >= self%t)
        call hpsort(tmp, order)
        deallocate(tmp)
        do i = nbox, 1, -1 ! because highest peaks are last
            ibox = order(i)
            if( selected_pos(ibox) )then
                mask = .false.
                !$omp parallel do schedule(static) default(shared) private(jbox,dist) proc_bind(close)
                do jbox = 1,nbox
                    dist = euclid(real(self%positions(pos_inds(ibox),:)),real(self%positions(pos_inds(jbox),:)))
                    if( dist <= self%dist_thres ) mask(jbox) = .true.
                end do
                !$omp end parallel do
                ! find best match in the neigh
                loc = maxloc(pos_scores, mask=mask, dim=1)
                ! eliminate all but the best
                mask(loc) = .false.
                where( mask ) selected_pos = .false.
            endif
        end do
        npeaks = count(selected_pos)
        write(logfhandle,'(a,1x,I5)') '# positions before distance filtering: ', nbox
        write(logfhandle,'(a,1x,I5)') '# positions after  distance filtering: ', npeaks
        ! update packed arrays
        pos_inds   = pack(pos_inds,   mask=selected_pos)
        pos_scores = pack(pos_scores, mask=selected_pos)
        ! update box_scores
        !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr,is_corr_peak,ipeak) proc_bind(close)
        do ioff = 1,self%nx_offset
            do joff = 1,self%ny_offset
                ithr = omp_get_thread_num() + 1
                is_corr_peak = .false.
                do ipeak = 1,npeaks
                    if( pos_inds(ipeak) == self%inds_offset(ioff,joff) )then
                        is_corr_peak = .true.
                        exit
                    endif
                end do
                if( .not. is_corr_peak ) self%box_scores(ioff,joff) = -1.
            end do
        end do
        !$omp end parallel do
        self%npeaks = count(self%box_scores >= self%t)
        if( L_DEBUG ) write(logfhandle,'(a,1x,I5)') '# positions after updating box_scores: ', self%npeaks
    end subroutine distance_filter

    subroutine remove_outliers( self )
        class(pickref), intent(inout) :: self
        real, allocatable :: tmp(:)
        integer     :: npeaks
        real        :: avg, sdev, t
        tmp = pack(self%loc_sdevs, mask=self%box_scores >= self%t .and. self%loc_sdevs > 0.)
        call avg_sdev(tmp, avg, sdev)
        t = avg + self%ndev * sdev
        npeaks = count(tmp < t)
        write(logfhandle,'(a,1x,I5)') '# positions after  outlier    removal: ', npeaks
        ! update box_scores
        !$omp parallel workshare proc_bind(close)
        where( self%loc_sdevs >= t )
            ! not a peak, update box_scores
            self%box_scores = -1.
        end where
        !$omp end parallel workshare
        if( self%nboxes_max > 0 .and. count(self%box_scores >= self%t) > self%nboxes_max )then
            tmp = pack(self%box_scores, mask=(self%box_scores > -1. + TINY))
            call detect_peak_thres_for_npeaks(size(tmp), self%nboxes_max, tmp, self%t)
            deallocate(tmp)
        endif
        self%npeaks = count(self%box_scores >= self%t)
        if( L_DEBUG ) write(logfhandle,'(a,1x,I5)') '# positions after updating box_scores: ', self%npeaks
    end subroutine remove_outliers

    subroutine peak_vs_nonpeak_stats( self )
        class(pickref), intent(inout) :: self
        real,    allocatable :: corrs_peak(:), corrs_nonpeak(:)
        integer, allocatable :: pos(:,:)
        real    :: pixrad_shrink, rpos(2)
        integer :: xrange(2), yrange(2), ibox, xind, yind, nx_offset, ny_offset, mask_count
        logical :: mask_backgr(0:self%nx,0:self%ny), mask_backgr_offset(self%nx_offset,self%ny_offset)
        ! prepare background mask
        call self%get_positions(pos)
        mask_backgr   = .true. 
        pixrad_shrink = (self%maxdiam / 2.) / self%smpd_shrink
        do ibox = 1,size(pos,dim=1)
            rpos      = real(pos(ibox,:))
            xrange(1) = max(0,       nint(rpos(1) - pixrad_shrink))
            xrange(2) = min(self%nx, nint(rpos(1) + pixrad_shrink))
            yrange(1) = max(0,       nint(rpos(2) - pixrad_shrink))
            yrange(2) = min(self%ny, nint(rpos(2) + pixrad_shrink))
            do xind = xrange(1),xrange(2)
                do yind = yrange(1),yrange(2)
                    mask_backgr(xind,yind) = .false.
                end do
            end do
        end do
        ! translate background mask to offset coordinates
        mask_backgr_offset = .false.
        mask_count = 0
        nx_offset = 0
        do xind = 0,self%nx,self%offset
            nx_offset = nx_offset + 1
            ny_offset = 0
            do yind = 0,self%ny,self%offset
                ny_offset = ny_offset + 1
                if( mask_backgr(xind,yind) .and. self%box_scores_mem(nx_offset,ny_offset) > -1. + TINY )then
                    mask_backgr_offset(nx_offset,ny_offset) = .true.
                    mask_count = mask_count + 1 ! not currently detecting any background coordinates
                endif
            end do
        end do
        ! extract scores
        corrs_peak        = pack(self%box_scores_mem, mask=self%box_scores >= self%t)
        corrs_nonpeak     = pack(self%box_scores_mem, mask=mask_backgr_offset)
        ! calc stats
        call avg_sdev(corrs_peak,        self%a_corr_peak,        self%s_corr_peak)
        call avg_sdev(corrs_nonpeak,     self%a_corr_nonpeak,     self%s_corr_nonpeak)
        self%smd_corr = std_mean_diff(self%a_corr_peak, self%a_corr_nonpeak, self%s_corr_peak, self%s_corr_nonpeak)
        call kstwo(corrs_peak, size(corrs_peak), corrs_nonpeak, size(corrs_nonpeak), self%ksstat_corr, self%prob_corr)
        write(logfhandle,'(a,1x,f4.2)') 'A_PEAK        = ', self%a_corr_peak
        write(logfhandle,'(a,1x,f4.2)') 'A_NONPEAK     = ', self%a_corr_nonpeak
        write(logfhandle,'(a,1x,f4.2)') 'S_PEAK        = ', self%s_corr_peak
        write(logfhandle,'(a,1x,f4.2)') 'S_NONPEAK     = ', self%s_corr_nonpeak
        write(logfhandle,'(a,1x,f4.2)') 'SMD           = ', self%smd_corr
        write(logfhandle,'(a,1x,f4.2)') 'K-S statistic = ', self%ksstat_corr
        write(logfhandle,'(a,1x,f4.2)') 'P             = ', self%prob_corr
        if( self%smd_corr < 0.2 .and. self%prob_corr > 0.5 ) write(logfhandle,'(a)') 'peak and non-peak distributions of CORRS are similar'
    end subroutine peak_vs_nonpeak_stats

    subroutine get_positions( self, pos, smpd_new )
        class(pickref),       intent(in)    :: self
        integer, allocatable, intent(inout) :: pos(:,:)
        real, optional,       intent(in)    :: smpd_new
        integer, allocatable :: pos_inds(:)
        real    :: scale
        integer :: nbox, ibox
        pos_inds = pack(self%inds_offset(:,:), mask=self%box_scores(:,:) >= self%t)
        nbox     = size(pos_inds)
        if( allocated(pos) ) deallocate(pos)
        allocate( pos(nbox,2),  source=0 )
        if( present(smpd_new) )then
            scale = self%smpd_shrink / smpd_new
        else
            scale = 1.
        endif
        do ibox = 1,nbox
            pos(ibox,:) = nint(scale * real(self%positions(pos_inds(ibox),:)))
        end do
    end subroutine get_positions

    subroutine get_loc_sdevs( self, loc_sdevs )
        class(pickref),       intent(in)    :: self
        real,    allocatable, intent(inout) :: loc_sdevs(:)
        if( allocated(loc_sdevs) ) deallocate(loc_sdevs)
        loc_sdevs = pack(self%loc_sdevs(:,:), mask=self%box_scores(:,:) >= self%t)
    end subroutine get_loc_sdevs

    subroutine get_scores( self, scores )
        class(pickref),       intent(in)    :: self
        real,    allocatable, intent(inout) :: scores(:)
        if( allocated(scores) ) deallocate(scores)
        scores = pack(self%box_scores(:,:), mask=self%box_scores(:,:) >= self%t)
    end subroutine get_scores

    integer function get_box( self )
        class(pickref), intent(in) :: self
        get_box = self%ldim_raw_box(1)
    end function get_box

    real function get_maxdiam( self )
        class(pickref), intent(in) :: self
        get_maxdiam = self%maxdiam
    end function get_maxdiam

    pure function get_nboxes( self ) result( nboxes )
        class(pickref), intent(in) :: self
        integer :: nboxes
        nboxes = self%npeaks
    end function get_nboxes

    ! for writing boxes with picking box size
    subroutine write_boxfile( self, fname )
        class(pickref), intent(in) :: self
        class(string),  intent(in) :: fname
        integer, allocatable :: pos(:,:)
        integer :: funit, ibox, iostat
        call self%get_positions(pos)
        call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=iostat)
        call fileiochk('simple_pickref; write_boxfile ', iostat)
        do ibox = 1,size(pos,dim=1)
            write(funit,'(I7,I7,I7,I7,I7)') pos(ibox,1), pos(ibox,2), self%ldim_box(1), self%ldim_box(1), -3
        end do
        call fclose(funit)
    end subroutine write_boxfile

    subroutine report_thumb_den( self, fname_thumb_den )
        use simple_micproc,   only: tv_filter_biomol
        use simple_gui_utils, only: mic2thumb
        class(pickref), intent(inout) :: self
        class(string),  intent(in)    :: fname_thumb_den 
        call tv_filter_biomol(self%mic_copy)
        call mic2thumb(self%mic_copy, fname_thumb_den, l_neg=self%l_black_ptcls_input) ! particles black
    end subroutine report_thumb_den
    
    ! for writing boxes with arbitrary box size
    subroutine report_boxfile( self, box, smpd, fname, nptcls )
        class(pickref), intent(in)  :: self
        integer,        intent(in)  :: box
        real,           intent(in)  :: smpd
        class(string),  intent(in)  :: fname
        integer,        intent(out) :: nptcls
        logical,   parameter :: l_write_outside = .false.
        integer, allocatable :: pos(:,:), updated_pos(:,:)
        logical, allocatable :: mask(:)
        integer :: i,j, ninside, funit, iostat
        nptcls  = self%npeaks
        if( nptcls == 0 ) return
        call self%get_positions(pos, smpd_new=smpd)
        if( box == self%ldim_raw_box(1) )then
            ! nothing to adjust
            updated_pos = pos
        else
            ! transform
            updated_pos = nint(real(pos) - real(box-self%ldim_raw_box(1))/2.)
            ! to remove positions outside micrograph
            if( l_write_outside )then
                ! flag positions
                allocate(mask(nptcls),source=.true.)
                do i = 1,nptcls
                    if( any(updated_pos(i,:) < 0) ) mask(i) = .false.
                    if( updated_pos(i,1)+self%ldim_raw_box(1)-1 >= ldim_raw(1)) mask(i) = .false.
                    if( updated_pos(i,2)+self%ldim_raw_box(2)-1 >= ldim_raw(2)) mask(i) = .false.
                enddo
                ninside = count(mask)
                if( ninside == 0 )then
                    ! all outside
                    nptcls = 0
                else if( ninside == nptcls )then
                    ! all inside
                else
                    ! some outside
                    deallocate(pos)
                    allocate(pos(ninside,2),source=0)
                    j = 0
                    do i = 1,nptcls
                        if( mask(i) ) then
                            j = j+1
                            pos(j,:) = updated_pos(i,:)
                        endif
                    enddo
                    updated_pos = pos
                    nptcls      = ninside
                endif
            endif
        endif
        ! write
        if( nptcls == 0 )then
            return
        else
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=iostat)
            call fileiochk('simple_pickref; report_boxfile ', iostat)
            do i = 1,nptcls
                write(funit,'(I7,I7,I7,I7,I7)')updated_pos(i,1),updated_pos(i,2),box,box,-3
            end do
            call fclose(funit)
        endif
    end subroutine report_boxfile

    subroutine refine_upscaled( self, pos, smpd_old, offset_old )
        class(pickref), intent(inout) :: self
        integer,        intent(in)    :: pos(:,:)
        real,           intent(in)    :: smpd_old
        integer,        intent(in)    :: offset_old
        integer, allocatable :: pos_refined(:,:)
        real,    allocatable :: scores_refined(:)
        type(image) :: boximgs_heap(nthr_glob)
        integer     :: nbox, ibox, jbox, ithr, xrange(2), yrange(2), xind, yind, iref, loc(1)
        real        :: factor, rpos(2), box_score, box_score_trial, scores(self%nrefs), dists(self%nboxes)
        logical     :: outside, l_err
        if( self%offset /= 1 ) THROW_HARD('Pixel offset in pickref instance subjected to refinement must be 1')
        nbox = size(pos, dim=1)
        allocate(pos_refined(nbox,2),  source= 0 )
        allocate(scores_refined(nbox), source=-1.)
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box, self%smpd_shrink)
        end do
        factor = real(offset_old) * (smpd_old / self%smpd_shrink)
        !$omp parallel do schedule(static) default(shared) proc_bind(close)&
        !$omp private(ibox,rpos,xrange,yrange,box_score,xind,yind,ithr,outside,iref,scores,box_score_trial,l_err)
        do ibox = 1,nbox
            rpos      = real(pos(ibox,:))
            xrange(1) = max(0,       nint(rpos(1) - factor))
            xrange(2) = min(self%nx, nint(rpos(1) + factor))
            yrange(1) = max(0,       nint(rpos(2) - factor))
            yrange(2) = min(self%ny, nint(rpos(2) + factor))
            box_score = -1.
            do xind = xrange(1),xrange(2)
                do yind = yrange(1),yrange(2)
                    ithr = omp_get_thread_num() + 1
                    call self%mic_shrink%window_slim([xind,yind], self%ldim_box(1), boximgs_heap(ithr), outside)
                    call boximgs_heap(ithr)%prenorm4real_corr(l_err)
                    if( l_err )then
                        box_score_trial = -1.
                    else
                        do iref = 1,self%nrefs
                            if( self%l_err_refs(iref) )then
                                scores(iref) = -1.
                            else
                                scores(iref) = self%boxrefs(iref)%real_corr_prenorm(boximgs_heap(ithr))
                            endif
                        end do
                        box_score_trial = maxval(scores)
                    endif
                    if( box_score_trial > box_score )then
                        pos_refined(ibox,:) = [xind,yind]
                        box_score = box_score_trial
                    endif
                end do
            end do
            scores_refined(ibox) = box_score
        end do
        !$omp end parallel do
        ! translate refined positions into instance
        self%t = minval(scores_refined)
        self%box_scores(:,:) = -1.
        do ibox = 1,nbox
            !$omp parallel do schedule(static) default(shared) private(jbox) proc_bind(close)
            do jbox = 1,self%nboxes
                dists(jbox) = euclid(real(pos_refined(ibox,:)),real(self%positions(jbox,:)))
            end do
            !$omp end parallel do
            loc = minloc(dists)
            where(self%inds_offset(:,:) == loc(1)) self%box_scores(:,:) = scores_refined(ibox)
        end do
        self%npeaks = count(self%box_scores >= self%t)
        write(logfhandle,'(a,1x,I5)') '# positions after refining upscaled:   ', self%npeaks
        ! destruct
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
    end subroutine refine_upscaled

    subroutine kill( self )
        class(pickref), intent(inout) :: self
        integer :: iimg
        if( self%exists )then
            call self%mic_shrink%kill
            call self%mic_roi%kill
            call self%mic_copy%kill
            if( allocated(self%l_mic_mask)     ) deallocate(self%l_mic_mask)
            if( allocated(self%l_err_refs)     ) deallocate(self%l_err_refs)
            if( allocated(self%positions)      ) deallocate(self%positions)
            if( allocated(self%inds_offset)    ) deallocate(self%inds_offset)
            if( allocated(self%box_scores)     ) deallocate(self%box_scores)
            if( allocated(self%box_scores_mem) ) deallocate(self%box_scores_mem)
            if( allocated(self%loc_sdevs)      ) deallocate(self%loc_sdevs)
            if( allocated(self%loc_sdevs_mem)  ) deallocate(self%loc_sdevs_mem)
            if( allocated(self%boxrefs) )then
                do iimg = 1,size(self%boxrefs)
                    call self%boxrefs(iimg)%kill
                end do
                deallocate(self%boxrefs)
            endif
            self%l_roi  = .false.
            self%exists = .false.
        endif
    end subroutine kill

    ! PUBLIC UTILITY


    ! Competitive reference picking
    ! the first picker contains the positions
    subroutine multiref_merge( np, pickers, pickind )
        integer,        intent(in)    :: np
        class(pickref), intent(inout) :: pickers(np)
        integer,        intent(inout) :: pickind
        integer              :: nx, ny, ioff, joff, ipick, c, cmax
        real,    allocatable :: scores(:,:)
        integer, allocatable :: map(:,:)
        nx = pickers(1)%nx_offset
        ny = pickers(1)%ny_offset
        allocate(map(nx,ny),    source=0)
        allocate(scores(nx,ny), source=-1.0)
        do ioff = 1,nx
            do joff = 1,ny
                do ipick = 1,np
                    if( pickers(ipick)%box_scores(ioff,joff) > -1. + TINY )then
                        if( pickers(ipick)%box_scores(ioff,joff) > scores(ioff,joff) )then
                            scores(ioff,joff) = pickers(ipick)%box_scores(ioff,joff)
                            map(ioff,joff)    = ipick
                        endif
                    endif
                end do
            end do
        end do
        ! output is first in the vector
        pickers(1)%box_scores(:,:) = scores
        pickers(1)%t = minval(scores, mask=scores >= 0.)
        ! apply distance filter to merged
        call pickers(1)%distance_filter
        ! vote for reference
        where( pickers(1)%box_scores(:,:) < pickers(1)%t ) map = 0
        pickind = 0
        cmax    = -1
        do ipick = 1,np
            c = count(map==ipick)
            if( c > cmax )then
                cmax    = c
                pickind = ipick
            endif
        enddo
        ! cleanup
        deallocate(scores,map)
    end subroutine multiref_merge

end module simple_pickref
