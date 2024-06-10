module simple_pickgau
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
implicit none

public :: read_mic_raw, pickgau, gaupick_multi
private
#include "simple_local_flags.inc"

! class constants
real,    parameter :: GAUSIG = 2.5, BOX_EXP_FAC = 0.111, NDEV_DEFAULT = 2.5
real,    parameter :: MSKDIAM2LP = 0.15, lP_LB = 30., LP_UB = 15.
integer, parameter :: OFFSET_DEFAULT = 3
logical, parameter :: L_WRITE  = .false.
logical, parameter :: L_DEBUG  = .false.

! class variables
integer                       :: ldim_raw(3)
real                          :: smpd_raw
type(image)                   :: mic_raw
character(len=:), allocatable :: fbody

! instance
type pickgau
    private
    ! control parameters
    real                     :: smpd_shrink = 0., moldiam = 0., maxdiam = 0., sig = 0., sxx = 0., t = 0., ndev = 0.
    real                     :: dist_thres  = 0., refine_dist_thres = 0.
    integer                  :: ldim(3), ldim_box(3), ldim_raw_box(3), nboxes = 0, nboxes_ub = 0, nx = 0, ny = 0
    integer                  :: nx_offset = 0, ny_offset = 0, npeaks = 0, nrefs = 0, offset = 0, offset_ub = 0
    integer                  :: peak_thres_level
    ! peak stats
    real                     :: smd_corr = 0., ksstat_corr = 0., prob_corr = 0.
    real                     :: a_corr_peak = 0., s_corr_peak = 0., a_corr_nonpeak = 0., s_corr_nonpeak = 0.
    ! images
    type(image)              :: mic_shrink, gauref, mic_roi
    type(image), allocatable :: boxrefs(:)
    ! masks
    logical,     allocatable :: l_mic_mask(:,:), l_err_refs(:)
    ! indices
    integer,     allocatable :: positions(:,:), inds_offset(:,:)
    ! scores & local standard deviations
    real,        allocatable :: box_scores(:,:), box_scores_mem(:,:), loc_sdevs(:,:), loc_sdevs_mem(:,:)
    ! flags
    logical                  :: l_roi   = .false.
    logical                  :: refpick = .false.
    logical                  :: exists  = .false.
contains
    procedure          :: gaupick
    procedure          :: new_gaupicker
    procedure          :: new_gaupicker_multi
    procedure          :: new_refpicker
    procedure, private :: new
    procedure, private :: set_refs
    procedure, private :: set_gaurefs
    procedure, private :: setup_iterators
    procedure, private :: flag_ice
    procedure, private :: flag_amorphous_carbon
    procedure, private :: gaumatch_boximgs
    procedure, private :: refmatch_boximgs
    procedure, private :: detect_peaks
    procedure, private :: distance_filter
    procedure, private :: remove_outliers
    procedure, private :: peak_vs_nonpeak_stats
    procedure          :: get_positions
    procedure          :: get_loc_sdevs
    procedure          :: get_scores
    procedure          :: get_nboxes
    procedure          :: get_box
    procedure          :: get_moldiam
    procedure, private :: write_boxfile
    procedure          :: report_boxfile
    procedure, private :: refine_upscaled
    procedure          :: kill
end type

contains

    subroutine gaupick_multi( pcontrast, smpd_shrink, moldiams, boxfile_multi, offset, ndev, moldiam_opt )
        character(len=*),  intent(in)  :: pcontrast
        real,              intent(in)  :: smpd_shrink, moldiams(:)
        character(len=*),  intent(in)  :: boxfile_multi
        integer, optional, intent(in)  :: offset
        real,    optional, intent(in)  :: ndev
        real,    optional, intent(out) :: moldiam_opt
        type(pickgau), allocatable :: pickers(:)
        integer,       allocatable :: picker_map(:,:)
        real,          allocatable :: smds_corr(:)
        type(pickgau) :: picker_merged
        integer :: npickers, ipick, ioff, joff, loc(1)
        real    :: moldiam_max, dist_threshold
        logical :: l_merge
        npickers    = size(moldiams)
        moldiam_max = maxval(moldiams)
        allocate(pickers(npickers), smds_corr(npickers))
        ! multi-gaussian pick
        do ipick = 1,npickers
            call pickers(ipick)%new_gaupicker(pcontrast, smpd_shrink, moldiams(ipick), moldiams(ipick), offset, ndev)
            ! call pickers(ipick)%new_gaupicker(pcontrast, smpd_shrink, moldiams(ipick), moldiam_max, offset, ndev)
            dist_threshold = (pickers(ipick)%maxdiam/3.) / smpd_shrink
            call pickers(ipick)%gaupick(dist_thres=dist_threshold)
            smds_corr(ipick) = pickers(ipick)%smd_corr
            call pickers(ipick)%mic_shrink%kill
        end do
        call report_multipick( npickers, pickers, boxfile_multi)
        ! return molecular diameter based on maximizing standard mean distance between peak- and
        ! non-peak distributions of corr
        loc = maxloc(smds_corr)
        moldiam_opt = moldiams(loc(1))
        write(logfhandle,'(a,1x,f7.2)') 'optimal molecular diameter identified:', moldiam_opt
        ! Is the following necessary since nothing is done with output?
        l_merge = .false.
        if( l_merge )then
            ! merged pick
            call picker_merged%new_gaupicker(pcontrast, smpd_shrink, moldiam_max, moldiam_max, offset, ndev)
            call picker_merged%setup_iterators
            allocate(picker_map(picker_merged%nx_offset,picker_merged%ny_offset), source=0)
            picker_merged%box_scores = -1.
            do ioff = 1,picker_merged%nx_offset
                do joff = 1,picker_merged%ny_offset
                    do ipick = 1,npickers
                        if( pickers(ipick)%box_scores(ioff,joff) > -1. + TINY )then
                            if( pickers(ipick)%box_scores(ioff,joff) > picker_merged%box_scores(ioff,joff) )then
                                picker_merged%box_scores(ioff,joff) = pickers(ipick)%box_scores(ioff,joff)
                                picker_map(ioff,joff) = ipick
                            endif
                        endif
                    end do
                end do
            end do
            picker_merged%t = minval(picker_merged%box_scores, mask=picker_merged%box_scores >= 0.)
            ! apply distance filter to merged
            call picker_merged%distance_filter
            ! update picker_map
            do ioff = 1,picker_merged%nx_offset
                do joff = 1,picker_merged%ny_offset
                    if( picker_merged%box_scores(ioff,joff) < picker_merged%t )then
                        picker_map(ioff,joff) = 0
                    endif
                end do
            end do
            ! cleanup
            call picker_merged%kill
        endif
        ! cleanup
        do ipick = 1,npickers
            call pickers(ipick)%kill
        end do
        deallocate(pickers)
    end subroutine gaupick_multi

    subroutine gaupick( self, self_refine, dist_thres )
        class(pickgau),           intent(inout) :: self
        class(pickgau), optional, intent(inout) :: self_refine
        real,           optional, intent(in)    :: dist_thres ! distance threshold that applies to self, not self_refine
        integer, allocatable :: pos(:,:)
        if( self%l_roi )then
            call self%flag_ice
            call self%flag_amorphous_carbon
            if( count(self%l_mic_mask) < nint(0.01*product(self%ldim)) )then
                self%npeaks = 0
                if( present(self_refine) ) self_refine%npeaks = 0
                return
            endif
        endif
        ! establish iterators
        call self%setup_iterators
        ! matching
        if( self%refpick )then
            call self%refmatch_boximgs
        else
            call self%gaumatch_boximgs
        endif
        call self%detect_peaks
        if( self%npeaks == 0 )then
            if( present(self_refine) ) self_refine%npeaks = 0
            return
        endif
        ! if( L_WRITE ) call self%write_boxfile('pickgau_after_detect_peaks.box')
        call self%distance_filter(dist_thres=dist_thres)
        ! if( L_WRITE ) call self%write_boxfile('pickgau_after_distance_filter.box')
        call self%remove_outliers
        ! if( L_WRITE ) call self%write_boxfile('pickgau_after_remove_outliers.box')
        call self%peak_vs_nonpeak_stats
        if( present(self_refine) )then
            call self%get_positions(pos, self_refine%smpd_shrink)
            call self_refine%setup_iterators
            call self_refine%refine_upscaled(pos, self%smpd_shrink, self%offset)
            call self_refine%distance_filter(dist_thres=self_refine%refine_dist_thres)
            ! if( L_WRITE ) call self_refine%write_boxfile('pickgau_after_refine_upscaled.box')
            deallocate(pos)
        endif
    end subroutine gaupick

    subroutine read_mic_raw( micname, smpd, subtr_backgr )
        character(len=*), intent(in) :: micname !< micrograph file name
        real,             intent(in) :: smpd    !< sampling distance in A
        logical, optional, intent(in) :: subtr_backgr
        character(len=:), allocatable :: ext
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
        ! set fbody
        ext   = fname2ext(trim(micname))
        fbody = trim(get_fbody(basename(trim(micname)), ext))
    end subroutine read_mic_raw

    subroutine new_gaupicker( self, pcontrast, smpd_shrink, moldiam, moldiam_max, offset, ndev, roi )
        class(pickgau),    intent(inout) :: self
        character(len=*),  intent(in)    :: pcontrast
        real,              intent(in)    :: smpd_shrink, moldiam
        real,    optional, intent(in)    :: moldiam_max
        integer, optional, intent(in)    :: offset
        real,    optional, intent(in)    :: ndev    !< # std devs for outlier detection
        logical, optional, intent(in)    :: roi
        call self%new( pcontrast, smpd_shrink, moldiam, moldiam_max, offset=offset, ndev=ndev, roi=roi )
    end subroutine new_gaupicker

    subroutine new_gaupicker_multi( self, pcontrast, smpd_shrink, moldiams, offset, ndev, roi)
        class(pickgau),    intent(inout) :: self
        character(len=*),  intent(in)    :: pcontrast
        real,              intent(in)    :: smpd_shrink
        real,              intent(in)    :: moldiams(:)
        integer, optional, intent(in)    :: offset
        real,    optional, intent(in)    :: ndev    !< # std devs for outlier detection
        logical, optional, intent(in)    :: roi
        real    :: moldiam_max, maxdiam
        integer :: box
        moldiam_max = maxval(moldiams)
        maxdiam     = moldiam_max + moldiam_max * BOX_EXP_FAC
        box         = round2even(maxdiam / smpd_raw)
        call self%new(pcontrast, smpd_shrink, moldiam_max, moldiam_max, box=box, offset=offset, ndev=ndev, roi=roi)
        call self%set_gaurefs(moldiams)
    end subroutine new_gaupicker_multi

    subroutine new_refpicker( self, pcontrast, smpd_shrink, imgs, offset, ndev, roi )
        class(pickgau),    intent(inout) :: self
        character(len=*),  intent(in)    :: pcontrast
        real,              intent(in)    :: smpd_shrink
        class(image),      intent(inout) :: imgs(:)
        integer, optional, intent(in)    :: offset
        real,    optional, intent(in)    :: ndev    !< # std devs for outlier detection
        logical, optional, intent(in)    :: roi
        integer :: ldim(3)
        ldim    = imgs(1)%get_ldim()
        params_glob%moldiam = (real(ldim(1)) - 2. * COSMSKHALFWIDTH ) * imgs(1)%get_smpd() 
        call self%new( pcontrast, smpd_shrink, params_glob%moldiam, params_glob%moldiam, box=ldim(1), offset=offset, ndev=ndev, roi=roi  )
        call self%set_refs( imgs )
    end subroutine new_refpicker

    subroutine new( self, pcontrast, smpd_shrink, moldiam, moldiam_max, box, offset, ndev, roi )
        class(pickgau),    intent(inout) :: self
        character(len=*),  intent(in)    :: pcontrast
        real,              intent(in)    :: smpd_shrink, moldiam, moldiam_max
        integer, optional, intent(in)    :: box     !< reference-based picking if present
        integer, optional, intent(in)    :: offset
        real,    optional, intent(in)    :: ndev    !< # std devs for outlier detection
        logical, optional, intent(in)    :: roi
        character(len=:), allocatable   :: numstr
        integer     :: box_max
        real        :: scale, maxdiam_max, lp
        if( self%exists ) call self%kill
        self%smpd_shrink = smpd_shrink
        if( present(roi) )then
            self%l_roi = roi
        else
            self%l_roi = .false.
        endif
        ! set peak threshold level
        if( trim(params_glob%crowded) .eq. 'yes' )then
            self%peak_thres_level = 1
        else
            self%peak_thres_level = 2
        endif
        ! set shrunken logical dimensions
        scale            = smpd_raw / self%smpd_shrink
        self%ldim(1)     = round2even(real(ldim_raw(1)) * scale)
        self%ldim(2)     = round2even(real(ldim_raw(2)) * scale)
        self%ldim(3)     = 1
        ! set logical dimensions of boxes
        self%moldiam = moldiam
        if( present(box) )then ! reference-based picking
            self%maxdiam      = smpd_raw * real(box)
            self%ldim_raw_box = [box,box,1]
            self%ldim_box(1)  = round2even(real(box) * scale)
            self%ldim_box(2)  = self%ldim_box(1)
            self%ldim_box(3)  = 1
            maxdiam_max       = self%maxdiam     ! only one box size considered
            box_max           = self%ldim_box(1) ! only one box size considered
        else
            self%maxdiam         = self%moldiam + self%moldiam * BOX_EXP_FAC
            self%ldim_raw_box(1) = round2even(self%maxdiam / smpd_raw)
            self%ldim_raw_box(2) = self%ldim_raw_box(1)
            self%ldim_raw_box(3) = 1
            self%ldim_box(1)     = round2even(self%maxdiam / self%smpd_shrink)
            self%ldim_box(2)     = self%ldim_box(1)
            self%ldim_box(3)     = 1
            maxdiam_max          = moldiam_max + moldiam_max * BOX_EXP_FAC    ! iteration bounds based on maximum box size
            box_max              = round2even(maxdiam_max / self%smpd_shrink) ! iteration bounds based on maximum box size
        endif
        self%nx = self%ldim(1) - box_max
        self%ny = self%ldim(2) - box_max
        ! set Gaussian
        self%sig = ((self%maxdiam / 2.) / self%smpd_shrink) / GAUSIG
        call self%gauref%new(self%ldim_box, self%smpd_shrink)
        call self%gauref%gauimg2D(self%sig, self%sig)
        call self%gauref%prenorm4real_corr(self%sxx)
        ! set ndev
        self%ndev = NDEV_DEFAULT
        if( present(ndev) ) self%ndev = ndev
        ! set offset
        self%offset = OFFSET_DEFAULT
        if( present(offset) ) self%offset = offset
        self%offset_ub = self%offset * 2
        ! set distance threshold
        if( self%offset == 1 )then
            self%dist_thres = 4.
        else
            self%dist_thres = self%offset
        endif
        self%refine_dist_thres = (self%maxdiam/3.)/self%smpd_shrink
        ! shrink micrograph
        call self%mic_shrink%new(self%ldim, self%smpd_shrink)
        call self%mic_shrink%set_ft(.true.)
        call mic_raw%mul(real(product(ldim_raw))) ! to prevent numerical underflow when performing FFT
        call mic_raw%fft
        call mic_raw%clip(self%mic_shrink)
        select case(trim(pcontrast))
            case('black')
                ! flip contrast (assuming black particle contrast on input)
                call self%mic_shrink%mul(-1.)
            case('white')
                ! nothing to do
            case DEFAULT
                THROW_HARD('uknown pcontrast parameter, use (black|white)')
        end select
        ! back to real-space
        call mic_raw%ifft
        call mic_raw%div(real(product(ldim_raw)))
        call self%mic_shrink%ifft
        if( L_WRITE )then
            if( present(box) )then ! reference-based picking
                call self%mic_shrink%write('mic_shrink.mrc')
                call self%gauref%write('gauref.mrc')
            else
                numstr = int2str(nint(self%moldiam))
                call self%mic_shrink%write('mic_shrink_moldiam'//numstr//'.mrc')
                call self%gauref%write('gauref_moldiam'//numstr//'.mrc')
            endif
        endif
        if( self%l_roi )then
            ! backup of shrunk micrograph
            call self%mic_roi%copy(self%mic_shrink)
        endif
        ! low-pass filter mic_shrink
        lp = min(max(LP_LB,MSKDIAM2LP * self%maxdiam),LP_UB)
        call self%mic_shrink%bp(0.,lp)
        if( L_WRITE )then
            if( present(box) )then ! reference-based picking
                call self%mic_shrink%write('mic_shrink_lp.mrc')
            else
                call self%mic_shrink%write('mic_shrink_lp_moldiam'//numstr//'.mrc')
            endif
        endif
        if( allocated(self%l_mic_mask) ) deallocate(self%l_mic_mask)
        allocate(self%l_mic_mask(self%ldim(1),self%ldim(2)), source=.true.)
        ! flag existence
        self%exists = .true.
    end subroutine new

    subroutine set_refs( self, imgs )
        class(pickgau), intent(inout) :: self
        class(image),   intent(inout) :: imgs(:)
        integer     :: ldim(3), iref
        ldim       = imgs(1)%get_ldim()
        if( ldim(3) /= 1 ) THROW_HARD('box references must be 2D')
        self%nrefs = size(imgs)
        ! set # pixels in x/y for both box sizes
        self%nx    = self%ldim(1) - self%ldim_box(1)
        self%ny    = self%ldim(2) - self%ldim_box(2)
        if( allocated(self%boxrefs) )then
            do iref = 1,size(self%boxrefs)
                call self%boxrefs(iref)%kill
            end do
            deallocate(self%boxrefs)
        endif
        if( allocated(self%l_err_refs) ) deallocate(self%l_err_refs)
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
        if( L_WRITE )then
            do iref = 1,self%nrefs
                call self%boxrefs(iref)%write('boxrefs_shrink.mrc', iref)
            enddo
        endif
        self%refpick = .true.
    end subroutine set_refs

    subroutine set_gaurefs(self, moldiams)
        class(pickgau), intent(inout) :: self
        real,           intent(in)    :: moldiams(:)
        integer                       :: nmoldiams, idiam, iimg
        real                          :: maxdiam, sig
        character(len=:), allocatable :: numstr
        nmoldiams = size(moldiams)
        self%nrefs = nmoldiams ! one ref per moldiam
        if( allocated(self%boxrefs) )then
            do iimg = 1,size(self%boxrefs)
                call self%boxrefs(iimg)%kill
            end do
            deallocate(self%boxrefs)
        endif
        allocate(self%boxrefs(self%nrefs))
        if( allocated(self%l_err_refs) ) deallocate(self%l_err_refs)
        allocate(self%l_err_refs(self%nrefs))
        do idiam=1,nmoldiams
            ! define constants used in subsequent calculations (per molecular diameter)
            maxdiam = moldiams(idiam) + moldiams(idiam) * BOX_EXP_FAC
            sig = maxdiam / self%smpd_shrink / GAUSIG
            ! create new gaussian reference image for each molecular diameter 
            call self%boxrefs(idiam)%new(self%ldim_box, self%smpd_shrink)
            call self%boxrefs(idiam)%gauimg2D(sig,sig)
            call self%boxrefs(idiam)%prenorm4real_corr(self%l_err_refs(idiam))
            if( L_WRITE )then
                numstr = int2str(nint(moldiams(idiam)))
                call self%boxrefs(idiam)%write('gauref_moldiam'//numstr//'.mrc')
            endif
        end do
        self%refpick = .true.
    end subroutine set_gaurefs

    subroutine setup_iterators( self )
        class(pickgau), intent(inout) :: self
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

    subroutine flag_ice( self )
        class(pickgau), intent(inout) :: self
        real,    parameter   :: THRESHOLD = 5.
        real,    parameter   :: SMPD_ICE  = ICE_BAND1/2. - 0.15
        integer, parameter   :: BOX = 128
        type(image)          :: boximgs_heap(nthr_glob), img
        real,    allocatable :: scores(:,:)
        integer, allocatable :: counts(:,:)
        real        :: score, scale, radius
        integer     :: ldim(3),i,j,ithr,ii,jj
        logical     :: outside
        if( smpd_raw > SMPD_ICE ) return
        scale   = self%smpd_shrink / SMPD_ICE
        ldim(1) = round2even(real(self%ldim(1)) * scale)
        ldim(2) = round2even(real(self%ldim(2)) * scale)
        ldim(3) = 1
        radius  = real(BOX)/2.- COSMSKHALFWIDTH
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new([BOX,BOX,1], SMPD_ICE, wthreads=.false.)
        end do
        call img%new(ldim,SMPD_ICE)
        call img%set_ft(.true.)
        call mic_raw%fft
        call mic_raw%clip(img)
        call mic_raw%ifft
        call img%ifft
        allocate(scores(ldim(1),ldim(2)),source=0.)
        allocate(counts(ldim(1),ldim(2)),source=0)
        !$omp parallel do collapse(2) schedule(static) default(shared) private(i,j,ithr,outside,score) proc_bind(close)
        do i = 1,ldim(1)-BOX+1,BOX/2
            do j = 1,ldim(2)-BOX+1,BOX/2
                ithr = omp_get_thread_num() + 1
                call img%window_slim([i-1, j-1], BOX, boximgs_heap(ithr), outside)
                call boximgs_heap(ithr)%mask(radius,'soft')
                call boximgs_heap(ithr)%fft
                call boximgs_heap(ithr)%calc_ice_score(score)
                scores(i:i+BOX-1,j:j+BOX-1) = scores(i:i+BOX-1,j:j+BOX-1) + score
                counts(i:i+BOX-1,j:j+BOX-1) = counts(i:i+BOX-1,j:j+BOX-1) + 1
            end do
        end do
        !$omp end parallel do
        i = ldim(1)-BOX+1
        !$omp parallel do schedule(static) default(shared) private(j,ithr,outside,score) proc_bind(close)
        do j = 1,ldim(2)-BOX+1,BOX/2
            ithr = omp_get_thread_num() + 1
            call img%window_slim([i-1, j-1], BOX, boximgs_heap(ithr), outside)
            call boximgs_heap(ithr)%mask(radius,'soft')
            call boximgs_heap(ithr)%fft
            call boximgs_heap(ithr)%calc_ice_score(score)
            scores(i:i+BOX-1,j:j+BOX-1) = scores(i:i+BOX-1,j:j+BOX-1) + score
            counts(i:i+BOX-1,j:j+BOX-1) = counts(i:i+BOX-1,j:j+BOX-1) + 1
        enddo
        !$omp end parallel do
        j = ldim(2)-BOX+1
        !$omp parallel do schedule(static) default(shared) private(i,ithr,outside,score) proc_bind(close)
        do i = 1,ldim(1)-BOX+1,BOX/2
            ithr = omp_get_thread_num() + 1
            call img%window_slim([i-1, j-1], BOX, boximgs_heap(ithr), outside)
            call boximgs_heap(ithr)%mask(radius,'soft')
            call boximgs_heap(ithr)%fft
            call boximgs_heap(ithr)%calc_ice_score(score)
            scores(i:i+BOX-1,j:j+BOX-1) = scores(i:i+BOX-1,j:j+BOX-1) + score
            counts(i:i+BOX-1,j:j+BOX-1) = counts(i:i+BOX-1,j:j+BOX-1) + 1
        enddo
        !$omp end parallel do
        where( counts > 0 ) scores = scores / real(counts)
        scale = real(ldim(1)) / real(self%ldim(1))
        do i = 1,self%ldim(1)
            ii = min(ldim(1), max(1, nint(scale*real(i))))
            do j = 1,self%ldim(2)
                jj = min(ldim(2), max(1, nint(scale*real(j))))
                self%l_mic_mask(i,j) = self%l_mic_mask(i,j) .and. (scores(ii,jj) < THRESHOLD)
            enddo
        enddo
        if( .not.all(self%l_mic_mask(:,:)) )then
            ! do i = 1,ldim(1)
            !     do j = 1,ldim(2)
            !         call img%set([i,j,1], scores(i,j))
            !     end do
            ! end do
            ! call img%write(trim(self%fbody)//'_ice.mrc')
            print *,' % ice water: ',&
                &100.-100.*count(self%l_mic_mask(:,:))/real(product(self%ldim))
        endif
        call img%kill
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%kill
        end do
    end subroutine flag_ice

    subroutine flag_amorphous_carbon( self )
        use simple_histogram, only: histogram
        class(pickgau), intent(inout) :: self
        integer, parameter :: K             = 5
        integer, parameter :: BOX           = 32 ! multiple of 4
        integer, parameter :: NBINS         = 64
        real,    parameter :: TVD_THRESHOLD = 0.2
        real,    parameter :: MIN_TVD_DIFF  = 0.05
        real,    parameter :: MIC_LP        = 15.
        type(histogram),  allocatable :: hists(:,:)
        logical,          allocatable :: final_mask(:,:)
        character(len=:), allocatable :: string
        type(ran_tabu)  :: rt
        type(histogram) :: khists(K)
        type(image)     :: mic, patches(nthr_glob)
        integer :: dims(3), pad, i
        logical :: found, empty
        found  = .false.
        dims   = self%ldim
        string = trim(fbody)
        allocate(final_mask(dims(1),dims(2)),source=.true.)
        call mic%copy(self%mic_roi)
        call clustering_rejection( found )
        if( found )then
            self%l_mic_mask = self%l_mic_mask .and. final_mask
            ! accounting for center of picking boxes
            pad = nint(real(self%ldim_box(1))/2.)
            final_mask = self%l_mic_mask
            self%l_mic_mask = .false.
            self%l_mic_mask(1:dims(1)-pad,1:dims(2)-pad) = final_mask(pad+1:,pad+1:)
        endif
        print *,string//' % amorphous carbon: ',100.-100.*count(final_mask(:,:))/real(product(dims))
        call mic%kill
        call rt%kill
        do i = 1,nthr_glob
            call patches(i)%kill
        enddo
        call khists(:)%kill
        call hists(:,:)%kill
        deallocate(hists)

        contains

            subroutine clustering_rejection( found )
                logical, intent(out) :: found
                integer, allocatable :: kmeans_labels(:), tmp_labels(:), inds(:)
                type(image)          :: tmpimg, tmpimg2
                type(histogram)      :: hist1, hist2
                real(dp) :: tmp
                real     :: bin_vars(K), tmpvec(K), tvd_inter(K-1)
                real     :: kmeans_score, mean, tvd, min_tvd
                logical  :: mask(dims(1),dims(2)), bin_mask(K), outside, cluster_mask(K-1)
                integer  :: b,twob,bon2,bon4,i,j,m,n,ilb,jlb,iub,jub,nx,ny,nxy
                integer  :: bin,repeat,ii,jj,carbon,ithr
                found = .false.
                ! gradients magnitude image
                call mic%fft
                call mic%bp(0.,MIC_LP)
                call mic%ifft
                call mic%norm
                call mic%gradients_magnitude(tmpimg)
                b     = BOX
                twob  = 2*b
                bon2  = b/2
                bon4  = b/4
                nx    = ceiling(real(dims(1))/real(b))
                ny    = ceiling(real(dims(2))/real(b))
                nxy   = nx*ny
                ! histograms
                allocate(kmeans_labels(nxy), tmp_labels(nxy),hists(nx,ny))
                !$omp parallel do collapse(2) proc_bind(close) private(i,j,ithr,ii,ilb,jlb,iub,jub)&
                !$omp default(shared)
                do i = 1,nx
                    do j =1,ny
                        ithr = omp_get_thread_num() + 1
                        ii  = (i-1)*nx+j
                        ilb = max(1, (i-1)*b-bon2+1)
                        iub = min(ilb+twob-1, dims(1))
                        ilb = iub-twob+1
                        jlb = max(1, (j-1)*b-bon2+1)
                        jub = min(jlb+twob-1, dims(2))
                        jlb = jub-twob+1
                        if(.not.patches(ithr)%exists() ) call patches(ithr)%new([twob,twob,1], self%smpd_shrink, wthreads=.false.)
                        call tmpimg%window_slim([ilb-1,jlb-1],twob, patches(ithr), outside)
                        call hists(i,j)%new(patches(ithr),  NBINS, minmax=[0.,4.])
                    enddo
                enddo
                !$omp end parallel do
                ! Cluster patches into K clusters
                kmeans_score = huge(0.)
                do i = 1,K
                    call khists(i)%new(hists(1,1))
                enddo
                call seed_rnd
                do repeat = 1,300
                    rt = ran_tabu(nxy)
                    call cluster(nxy, nx, ny, K, tmp, tmp_labels)
                    if( real(tmp) < kmeans_score )then
                        kmeans_score  = real(tmp)
                        kmeans_labels = tmp_labels
                    endif
                enddo
                do bin = 1,K
                    bin_mask(bin) = count(kmeans_labels==bin) > 0
                enddo
                if( count(bin_mask) == 1 ) return ! clustering fail, do nothing
                ! histograms for all clusters
                do i = 1,K
                    call khists(i)%zero
                enddo
                bin = 0
                do i = 1,nx
                    do j =1,ny
                        bin = bin+1
                        call khists(kmeans_labels(bin))%add(hists(i,j))
                    enddo
                enddo
                ! deviation from mode
                bin_vars = -1.
                do bin = 1,K
                    if( bin_mask(bin) )then
                        mean          = khists(bin)%hmode()
                        bin_vars(bin) = khists(bin)%variance(mean=mean)
                    endif
                enddo
                ! sorting by cluster variance
                inds   = (/(i,i=1,K)/)
                tmpvec = bin_vars
                call hpsort(tmpvec, inds)
                ! debug
                if( L_DEBUG )then
                    do i = 1,K
                        call tmpimg2%copy(tmpimg)
                        j = 0
                        do ii = 1,nx
                            ilb = max(1, (ii-1)*b+1)
                            iub = min(ilb+b-1, dims(1))
                            ilb = iub-b+1
                            do jj =1,ny
                                j = j + 1
                                if( kmeans_labels(j) == inds(i) ) cycle
                                jlb = max(1, (jj-1)*b+1)
                                jub = min(jlb+b-1, dims(2))
                                jlb = jub-b+1
                                do m = ilb,iub
                                    do n = jlb,jub
                                        call tmpimg2%set([m,n,1],0.)
                                    enddo
                                enddo
                            enddo
                        enddo
                        call tmpimg2%write(trim(string)//'_clusters.mrc',i)
                    enddo
                    call mic%write(trim(string)//'_clusters.mrc',K+1)
                    call tmpimg2%kill
                endif
                ! Segmentation by agglomerative inter-cluster distance
                tvd_inter = -.1
                do i = 1,K-1
                    empty = .true.
                    call hist1%new(khists(1))
                    do j = 1,i
                        jj = inds(j)
                        if( bin_mask(jj) )then
                            call hist1%add(khists(jj))
                            empty = .false.
                        endif
                    enddo
                    if( empty ) cycle
                    empty = .true.
                    call hist2%new(khists(1))
                    do j = i+1,K
                        jj = inds(j)
                        if( bin_mask(jj) )then
                            call hist2%add(khists(jj))
                            empty= .false.
                        endif
                    enddo
                    if( empty ) cycle
                    tvd_inter(i) = hist1%tvd(hist2)
                    if( L_DEBUG )then
                        write(*,'(I6,5F9.4)') i, hist1%variance(), hist2%variance(), hist1%tvd(hist2),&
                            &khists(inds(i))%variance(), khists(inds(i))%mean()
                    endif
                enddo
                cluster_mask = tvd_inter > 0.
                carbon       = maxloc(tvd_inter,dim=1)
                tvd          = tvd_inter(carbon)
                min_tvd      = minval(tvd_inter,mask=cluster_mask)
                bin          = findloc(cluster_mask(1:K-1),.true.,dim=1) ! first non-empty cluster
                ! Decision
                found = .false.
                if( bin == 0 )then
                    tvd = 0.
                else
                    if( carbon == bin ) tvd = 0.
                    if( abs(tvd-min_tvd) < MIN_TVD_DIFF ) tvd = 0.
                endif
                found = tvd > TVD_THRESHOLD
                ! Dejection
                if( found )then
                    mask = final_mask
                    do bin = carbon+1,K
                        n = inds(bin)
                        if( .not.bin_mask(n) ) cycle
                        m     = 0
                        do i = 1,nx
                            ii = min((i-1)*b+1, dims(1)-b+1)
                            do j =1,ny
                                jj = min((j-1)*b+1, dims(2)-b+1)
                                m = m+1
                                if( kmeans_labels(m) == n ) mask(ii:ii+b-1,jj:jj+b-1) = .false.
                            enddo
                        enddo
                    enddo
                    ! erosion with half window size
                    final_mask = mask
                    do i = 1,2*nx
                        ilb = max(1, (i-1)*bon2-bon4+1)
                        iub = min(ilb+b-1, dims(1))
                        ilb = iub - b + 1 + bon4
                        iub = iub - bon4
                        do j = 1,2*ny
                            jlb = max(1, (j-1)*bon2-bon4+1)
                            jub = min(jlb+b-1, dims(2))
                            jlb = jub - b + 1 + bon4
                            jub = jub - bon4
                            n = count(mask(ilb-bon4:iub+bon4,jlb-bon4:jub+bon4))
                            if( n == 0   )cycle
                            if( n == b*b )cycle
                            final_mask(ilb:iub,jlb:jub) = .true.
                        enddo
                    enddo
                    ! dilation with half window size x 2
                    mask = final_mask
                    do i = 1,2*nx
                        ilb = max(1, (i-1)*bon2-bon4+1)
                        iub = min(ilb+b-1, dims(1))
                        ilb = iub - b + 1 + bon4
                        iub = iub - bon4
                        do j = 1,2*ny
                            jlb = max(1, (j-1)*bon2-bon4+1)
                            jub = min(jlb+b-1, dims(2))
                            jlb = jub - b + 1 + bon4
                            jub = jub - bon4
                            if( all(mask(ilb:iub,jlb:jub)) )then
                                n = count(.not.mask(ilb-bon4:iub+bon4,jlb-bon4:jub+bon4))
                                if( n > 0 )then
                                    final_mask(ilb:iub,jlb:jub) = .false.
                                endif
                            endif
                        enddo
                    enddo
                    mask = final_mask
                    do i = 1,2*nx
                        ilb = max(1, (i-1)*bon2-bon4+1)
                        iub = min(ilb+b-1, dims(1))
                        ilb = iub - b + 1 + bon4
                        iub = iub - bon4
                        do j = 1,2*ny
                            jlb = max(1, (j-1)*bon2-bon4+1)
                            jub = min(jlb+b-1, dims(2))
                            jlb = jub - b + 1 + bon4
                            jub = jub - bon4
                            if( all(mask(ilb:iub,jlb:jub)) )then
                                n = count(.not.mask(ilb-bon4:iub+bon4,jlb-bon4:jub+bon4))
                                if( n > 0 )then
                                    final_mask(ilb:iub,jlb:jub) = .false.
                                endif
                            endif
                        enddo
                    enddo
                endif
                ! cleanup
                call hist1%kill
                call hist2%kill
            end subroutine clustering_rejection

            subroutine cluster( n, nx, ny, K, overall_score, labels )
                integer,         intent(in)  :: n, nx, ny, K
                real(dp),        intent(out) :: overall_score
                integer,         intent(out) :: labels(n)
                integer, parameter :: nits = 100
                integer  :: pops(K), new_labels(n)
                real(dp) :: scores(K), kscores(K), prev_score, score
                integer  :: i,j,c,cl,it,nk
                ! initial labels
                c = 0
                do i = 1,n
                    c = c+1
                    if( c > K ) c = 1
                    labels(i) = c
                enddo
                call rt%shuffle(labels)
                ! main loop
                prev_score = huge(0.0)
                do it = 1,nits
                    scores = 0.d0
                    !$omp parallel private(i,j,c,cl,kscores) default(shared) proc_bind(close)
                    !$omp do
                    ! centers
                    do cl = 1,K
                        call khists(cl)%zero
                        pops(cl) = count(labels==cl)
                        if( pops(cl) == 0 ) cycle
                        do i = 1,nx
                            do j = 1,ny
                                c = (i-1)*ny+j
                                if( labels(c) /= cl )cycle
                                call khists(cl)%add(hists(i,j))
                            enddo
                        enddo
                    enddo
                    !$omp end do
                    !$omp do collapse(2) reduction(+:scores)
                    do i = 1,nx
                        do j = 1,ny
                            c = (i-1)*ny+j
                            do cl = 1,K
                                if( pops(cl) > 0 )then
                                    kscores(cl) = real(khists(cl)%tvd(hists(i,j)),dp)
                                else
                                    kscores(cl) = huge(0.)
                                endif
                            enddo
                            ! current parttioning
                            cl = labels(c)
                            scores(cl) = scores(cl) + kscores(cl)
                            ! new labels
                            new_labels(c) = minloc(kscores,dim=1)
                        enddo
                    enddo
                    !$omp end do
                    !$omp end parallel
                    nk     = count(pops>0)
                    score  = sum(scores/real(pops,dp),mask=pops>0) / real(nk,dp)
                    ! convergence
                    if( score < prev_score )then
                        if( abs(score - prev_score) < 1.d-6 ) exit
                    endif
                    prev_score = score
                    labels     = new_labels
                enddo
                overall_score = score
            end subroutine cluster

    end subroutine flag_amorphous_carbon

    subroutine gaumatch_boximgs( self )
        class(pickgau), intent(inout) :: self
        logical     :: outside
        integer     :: pos(2), ioff, joff, ithr
        type(image) :: boximgs_heap(nthr_glob)
        ! construct heap
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box, self%smpd_shrink)
        end do
        ! calculate box_scores
        !$omp parallel do schedule(static) collapse(2) default(shared) private(ioff,joff,ithr,outside,pos) proc_bind(close)
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
                        self%loc_sdevs(ioff,joff)  = boximgs_heap(ithr)%avg_loc_sdev(self%offset)
                        self%box_scores(ioff,joff) = self%gauref%real_corr_prenorm(boximgs_heap(ithr), self%sxx)
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
    end subroutine gaumatch_boximgs

    subroutine refmatch_boximgs( self )
        class(pickgau), intent(inout) :: self
        logical     :: outside, l_err
        integer     :: pos(2), ioff, joff, ithr, iref
        real        :: scores(self%nrefs)
        type(image) :: boximgs_heap(nthr_glob)
        ! construct heap
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box, self%smpd_shrink)
        end do
        ! calculate box_scores
        if( self%refpick )then
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
        else
            THROW_HARD('Instance not setup for reference-based picking')
        endif
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
    end subroutine refmatch_boximgs

    subroutine detect_peaks( self )
        class(pickgau), intent(inout) :: self
        real, allocatable :: tmp(:)
        integer :: n
        tmp = pack(self%box_scores, mask=(self%box_scores > -1. + TINY))
        n   = size(tmp)
        if( n == 0 )then
            self%npeaks = 0
            return
        endif
        call detect_peak_thres(n, self%nboxes_ub, self%peak_thres_level, tmp, self%t)
        deallocate(tmp)
        self%t = max(0.,self%t)
        where( self%box_scores >= self%t )
            ! there's a peak
        elsewhere
            self%box_scores = -1.
        end where
        self%npeaks = count(self%box_scores >= self%t)
        write(logfhandle,'(a,1x,f5.2)') 'peak threshold identified:             ', self%t
        write(logfhandle,'(a,1x,I5)'  ) '# peaks detected:                      ', self%npeaks
    end subroutine detect_peaks

    subroutine distance_filter( self, dist_thres )
        class(pickgau), intent(inout) :: self
        real, optional, intent(in)    :: dist_thres
        integer, allocatable :: pos_inds(:)
        real,    allocatable :: pos_scores(:)
        logical, allocatable :: mask(:), selected_pos(:)
        integer :: nbox, npeaks, ibox, jbox, loc, ioff, joff, ithr, ipeak
        real    :: dist, dthres
        logical :: is_corr_peak
        if( present(dist_thres) )then
            dthres = dist_thres
        else
            dthres = self%dist_thres
        endif
        pos_inds   = pack(self%inds_offset(:,:), mask=self%box_scores(:,:) >= self%t)
        pos_scores = pack(self%box_scores(:,:),  mask=self%box_scores(:,:) >= self%t)
        nbox       = size(pos_inds)
        allocate(mask(nbox),         source=.false.)
        allocate(selected_pos(nbox), source=.true. )
        do ibox = 1,nbox
            mask = .false.
            !$omp parallel do schedule(static) default(shared) private(jbox,dist) proc_bind(close)
            do jbox = 1,nbox
                dist = euclid(real(self%positions(pos_inds(ibox),:)),real(self%positions(pos_inds(jbox),:)))
                if( dist <= dthres ) mask(jbox) = .true.
            end do
            !$omp end parallel do
            ! find best match in the neigh
            loc = maxloc(pos_scores, mask=mask, dim=1)
            ! eliminate all but the best
            mask(loc) = .false.
            where( mask ) selected_pos = .false.
        end do
        npeaks = count(selected_pos)
        if( L_DEBUG ) write(logfhandle,'(a,1x,I5)') '# positions before distance filtering: ', nbox
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
        class(pickgau), intent(inout) :: self
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
        self%npeaks = count(self%box_scores >= self%t)
        if( L_DEBUG ) write(logfhandle,'(a,1x,I5)') '# positions after updating box_scores: ', self%npeaks
    end subroutine remove_outliers

    subroutine peak_vs_nonpeak_stats( self )
        class(pickgau), intent(inout) :: self
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
        class(pickgau),       intent(in)    :: self
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
        class(pickgau),       intent(in)    :: self
        real,    allocatable, intent(inout) :: loc_sdevs(:)
        if( allocated(loc_sdevs) ) deallocate(loc_sdevs)
        loc_sdevs = pack(self%loc_sdevs(:,:), mask=self%box_scores(:,:) >= self%t)
    end subroutine get_loc_sdevs

    subroutine get_scores( self, scores )
        class(pickgau),       intent(in)    :: self
        real,    allocatable, intent(inout) :: scores(:)
        if( allocated(scores) ) deallocate(scores)
        scores = pack(self%box_scores(:,:), mask=self%box_scores(:,:) >= self%t)
    end subroutine get_scores

    integer function get_box( self )
        class(pickgau), intent(in) :: self
        get_box = self%ldim_raw_box(1)
    end function get_box

    real function get_moldiam( self )
        class(pickgau), intent(in) :: self
        get_moldiam = self%moldiam
    end function get_moldiam

    pure function get_nboxes( self ) result( nboxes )
        class(pickgau), intent(in) :: self
        integer :: nboxes
        nboxes = self%npeaks
    end function get_nboxes

    ! for writing boxes with picking box size
    subroutine write_boxfile( self, fname )
        class(pickgau),   intent(in) :: self
        character(len=*), intent(in) :: fname
        integer, allocatable :: pos(:,:)
        integer :: funit, ibox, iostat
        call self%get_positions(pos)
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(fname)), iostat=iostat)
        call fileiochk('simple_pickgau; write_boxfile ', iostat)
        do ibox = 1,size(pos,dim=1)
            write(funit,'(I7,I7,I7,I7,I7)') pos(ibox,1), pos(ibox,2), self%ldim_box(1), self%ldim_box(1), -3
        end do
        call fclose(funit)
    end subroutine write_boxfile

    ! for writing boxes with arbitrary box size
    subroutine report_boxfile( self, box, smpd, fname, nptcls )
        class(pickgau),            intent(in)  :: self
        integer,                   intent(in)  :: box
        real,                      intent(in)  :: smpd
        character(len=LONGSTRLEN), intent(in)  :: fname
        integer,                   intent(out) :: nptcls
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
            call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(fname)), iostat=iostat)
            call fileiochk('simple_pickgau; report_boxfile ', iostat)
            do i = 1,nptcls
                write(funit,'(I7,I7,I7,I7,I7)')updated_pos(i,1),updated_pos(i,2),box,box,-3
            end do
            call fclose(funit)
        endif
    end subroutine report_boxfile

    subroutine refine_upscaled( self, pos, smpd_old, offset_old )
        class(pickgau), intent(inout) :: self
        integer,        intent(in)    :: pos(:,:)
        real,           intent(in)    :: smpd_old
        integer,        intent(in)    :: offset_old
        integer, allocatable :: pos_refined(:,:)
        real,    allocatable :: scores_refined(:)
        type(image) :: boximgs_heap(nthr_glob)
        integer     :: nbox, ibox, jbox, ithr, xrange(2), yrange(2), xind, yind, iref, loc(1)
        real        :: factor, rpos(2), box_score, box_score_trial, scores(self%nrefs), dists(self%nboxes)
        logical     :: outside, l_err
        if( self%offset /= 1 ) THROW_HARD('Pixel offset in pickgau instance subjected to refinement must be 1')
        nbox = size(pos, dim=1)
        allocate(pos_refined(nbox,2),  source= 0 )
        allocate(scores_refined(nbox), source=-1.)
        do ithr = 1,nthr_glob
            call boximgs_heap(ithr)%new(self%ldim_box, self%smpd_shrink)
        end do
        factor = real(offset_old) * (smpd_old / self%smpd_shrink)
        if( self%refpick )then
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
        else
            !$omp parallel do schedule(static) default(shared) proc_bind(close)&
            !$omp private(ibox,rpos,xrange,yrange,box_score,xind,yind,ithr,outside,box_score_trial)
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
                        box_score_trial = self%gauref%real_corr_prenorm(boximgs_heap(ithr), self%sxx)
                        if( box_score_trial > box_score )then
                            pos_refined(ibox,:) = [xind,yind]
                            box_score = box_score_trial
                        endif
                    end do
                end do
                scores_refined(ibox) = box_score
            end do
            !$omp end parallel do
        endif
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
        class(pickgau), intent(inout) :: self
        integer :: iimg
        if( self%exists )then
            call self%mic_shrink%kill
            call self%gauref%kill
            call self%mic_roi%kill
            if( allocated(self%l_mic_mask)     ) deallocate(self%l_mic_mask)
            if( allocated(self%l_err_refs)     ) deallocate(self%l_err_refs)
            if( allocated(self%positions)      ) deallocate(self%positions)
            if( allocated(self%inds_offset)    ) deallocate(self%inds_offset)
            if( allocated(self%box_scores)     ) deallocate(self%box_scores)
            if( allocated(self%box_scores_mem) ) deallocate(self%box_scores_mem)
            if( allocated(self%loc_sdevs)     ) deallocate(self%loc_sdevs)
            if( allocated(self%loc_sdevs_mem) ) deallocate(self%loc_sdevs_mem)
            if( allocated(self%boxrefs) )then
                do iimg = 1,size(self%boxrefs)
                    call self%boxrefs(iimg)%kill
                end do
                deallocate(self%boxrefs)
            endif
            self%refpick = .false.
            self%exists  = .false.
        endif
    end subroutine kill

    !
    ! Private module utilities
    !

    subroutine report_multipick( npick, picker, fname )
        integer,                   intent(in)  :: npick
        class(pickgau),            intent(in)  :: picker(npick)
        character(len=*), intent(in)  :: fname
        integer, allocatable :: pos(:,:)
        real,    allocatable :: scores(:), loc_sdevs(:)
        real    :: moldiam
        integer :: i,j, funit, iostat, npos, box, box_raw
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(fname)), iostat=iostat)
        call fileiochk('simple_pickgau; report_multipick ', iostat)
        do i = 1,npick
            npos = picker(i)%get_nboxes()
            if( npos == 0 ) cycle
            call picker(i)%get_positions(pos, smpd_new=smpd_raw)
            call picker(i)%get_scores(scores)
            call picker(i)%get_loc_sdevs(loc_sdevs)
            box_raw = picker(i)%get_box()
            moldiam = picker(i)%get_moldiam()
            ! to account tight box
            box = find_larger_magic_box( nint(params_glob%box_exp_factor*real(box_raw)) )
            if( box /= box_raw ) pos = nint(real(pos) - real(box-box_raw)/2.)
            ! write
            do j = 1,npos
                write(funit,'(I7,I7,I7,F8.1,F8.3,F12.3)') pos(j,1),pos(j,2),box,moldiam,scores(j),loc_sdevs(j)
            enddo
        enddo
        call fclose(funit)
    end subroutine report_multipick

end module simple_pickgau
