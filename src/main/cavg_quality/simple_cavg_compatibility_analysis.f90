!@descr: class-average compatibility analysis, support-model training, and rejection
module simple_cavg_compatibility_analysis
! Module overview:
! - Preprocess class averages and derive geometric/statistical descriptors.
! - Reject class averages with support and cluster-level criteria.
! - Write selected/rejected stacks and per-class rejection reasons.
use unix,                    only: c_float
use simple_defs,             only: logfhandle
use simple_error,            only: simple_exception
use simple_image,            only: image
use simple_string,           only: string
use simple_fileio,           only: swap_suffix
use simple_cmdline,          only: cmdline
use simple_gui_utils,        only: mrc2jpeg_tiled
use simple_image_bin,        only: image_bin
use simple_defs_fname,       only: METADATA_EXT, MRC_EXT, JPG_EXT
use simple_sp_project,       only: sp_project
use simple_imgarr_utils,     only: read_cavgs_into_imgarr, read_stk_into_imgarr, write_imgarr, dealloc_imgarr
use simple_segmentation,     only: otsu_img
use simple_srch_sort_loc,    only: hpsort
use simple_commanders_cavgs, only: commander_cluster_cavgs
use simple_cavg_quality_types, only: cavg_quality_result
use simple_cavg_quality_analysis, only: evaluate_cavg_quality_hard_reject

implicit none

public :: cavg_compatibility_analysis

private
#include "simple_local_flags.inc"

! Constants for cavg compatibility analysis
integer, parameter :: ANALYSIS_BOXSIZE             = 128    ! Working box size for resized class averages (pixels)
integer, parameter :: ANALYSIS_MORPH_SIZE          = 5      ! Number of dilate/erode passes for morphological closing
logical, parameter :: ANALYSIS_AUTOTUNE_SIZE_MODEL = .true. ! Enable grid-search autotuning of size-support parameters

! Size-subset support model constants.
real,    parameter :: RESCUE_EDGE_FRAC       = 0.00                          ! Extra margin for boundary rescue in size compatibility
integer, parameter :: NRELAX                 = 5                             ! Number of relaxation candidates in support-model grid
integer, parameter :: NQ                     = 3                             ! Number of lower/upper quantile candidates in grid
real,    parameter :: SUPPORT_EDGE_SOFT_FRAC = 0.0                           ! Soft edge slack for c/a bounds during support check
real,    parameter :: RELAX_GRID(NRELAX)     = [0.03, 0.05, 0.07, 0.1, 0.15] ! Relative interval expansion candidates
real,    parameter :: QLOW_GRID(NQ)          = [0.05, 0.1, 0.15]             ! Lower-quantile candidates for c-axis estimate
real,    parameter :: QHIGH_GRID(NQ)         = [0.85, 0.9, 0.95]             ! Upper-quantile candidates for a-axis estimate

! Cluster overfitting rejection constants.
real,    parameter :: LOWVAR_CLUSTER_REJECT_FRAC          = 0.60 ! Reject cluster when low-variance fraction exceeds this value
integer, parameter :: MIN_CLUSTER_REJECTION_CAVGS         = 10   ! Minimum selected class averages required to run cluster overfitting rejection.

! Rejection reason codes for image_set
integer, parameter :: REJECT_REASON_NONE                  = 0    ! Rejection code for non-rejected images.
integer, parameter :: REJECT_REASON_ZERO_VARIANCE         = 1    ! Rejection code for zero image variance.
integer, parameter :: REJECT_REASON_MASK_OUTSIDE_SUPPORT  = 2    ! Rejection code for mask outside circular support.
integer, parameter :: REJECT_REASON_SIZE_INCOMPATIBLE     = 3    ! Rejection code for size-incompatible subset.
integer, parameter :: REJECT_REASON_SUSPECTED_OVERFITTING = 4    ! Rejection code for cluster-level overfitting suspicion.
integer, parameter :: REJECT_REASON_QUALITY_HARD_REJECT   = 5    ! Rejection code for quality hard-reject gate.

type :: support_model_state
    logical              :: do_autotune = ANALYSIS_AUTOTUNE_SIZE_MODEL
    logical              :: valid       = .false.
    real                 :: axis_c      = 0.0
    real                 :: axis_b      = 0.0
    real                 :: axis_a      = 0.0
    real                 :: used_relax  = 0.07
    real                 :: used_qlo    = 0.08
    real                 :: used_qhi    = 0.93
end type support_model_state

type :: support_training_set_state
    real,    allocatable :: min_dim(:)
    real,    allocatable :: max_dim(:)
    integer, allocatable :: batch_id(:)
    integer              :: nbatches = 0
end type support_training_set_state

type :: image_set
    type(image)          :: img                          ! Processed class-average image used for analysis and outputs.
    type(image_bin)      :: mask                         ! Binary foreground mask derived from Otsu + morphology + LCC.
    integer              :: rejection_reason   = REJECT_REASON_NONE ! Integer rejection code used in reports/logging.
    real                 :: variance           = 0.0     ! Global variance of the processed class-average image.
    real                 :: local_var_in_mask  = 0.0     ! Local variance statistic measured inside the mask.
    real                 :: local_var_out_mask = 0.0     ! Local variance statistic measured outside the mask.
    real                 :: feret_max          = 0.0     ! Maximum Feret diameter of the binary mask.
    real                 :: feret_min          = 0.0     ! Minimum Feret diameter of the binary mask.
    logical              :: is_rejected        = .false. ! True when this image_set is rejected by any criterion.
end type image_set

type :: cavg_compatibility_analysis
    private
    type(image_set),     allocatable :: imagesets(:)         ! Per-class analysis state for all input class averages.
    type(support_model_state)        :: support_model        ! Support model used by infer() through apply_support_model().
    type(support_training_set_state) :: support_training_set ! Training pool used by generate_support_model().
    type(sp_project)                 :: spproj               ! Project handle used to read/write class-average metadata.
    type(string)                     :: input_stkname        ! Input class-average stack filename.
    integer                          :: nimagesets = 0       ! Number of class averages stored in imagesets.
    
contains
    procedure, public :: new
    procedure, public :: kill
    procedure, public :: kill_support_training_set
    procedure, public :: kill_support_model
    procedure, public :: infer
    procedure, public :: get_rejection_states
    generic,   public :: train => train_1, train_2
    procedure, private :: kill_imagesets
    procedure, private :: train_1
    procedure, private :: train_2
    procedure, private :: cache_support_training_set
    procedure, private :: generate_support_model
    procedure, private :: print_rejection_reasons
    procedure, private :: preprocess
    procedure, private :: write_selected_rejected_stacks
    procedure, private :: run_cluster_overfitting_rejection
end type cavg_compatibility_analysis

contains

    ! Initialize analysis object state.
    subroutine new( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        call self%kill()
    end subroutine new

    ! Release owned resources and reset object state.
    subroutine kill( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        call self%input_stkname%kill()
        call self%spproj%kill()
        call self%kill_imagesets()
        call self%kill_support_training_set()
        call self%kill_support_model()
    end subroutine kill

    ! Release processed image/mask buffers and reset image counter.
    subroutine kill_imagesets( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        integer :: i
        if( allocated(self%imagesets) )then
            do i = 1, size(self%imagesets)
                call self%imagesets(i)%img%kill()
                call self%imagesets(i)%mask%kill_bimg()
            end do
            deallocate(self%imagesets)
        end if
        self%nimagesets = 0
    end subroutine kill_imagesets

    ! Release module-level training arrays for support-model fitting.
    subroutine kill_support_training_set( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        if( allocated(self%support_training_set%min_dim)  ) deallocate(self%support_training_set%min_dim)
        if( allocated(self%support_training_set%max_dim)  ) deallocate(self%support_training_set%max_dim)
        if( allocated(self%support_training_set%batch_id) ) deallocate(self%support_training_set%batch_id)
        self%support_training_set%nbatches = 0
    end subroutine kill_support_training_set

    ! Reset module-level support model state.
    subroutine kill_support_model( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        self%support_model%axis_c     = 0.0
        self%support_model%axis_b     = 0.0
        self%support_model%axis_a     = 0.0
        self%support_model%used_relax = 0.07
        self%support_model%used_qlo   = 0.08
        self%support_model%used_qhi   = 0.93
        self%support_model%valid      = .false.
    end subroutine kill_support_model

    ! Append current non-rejected min/max Feret dimensions to the training pool.
    subroutine cache_support_training_set( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        integer :: ii, jj, nkeep, nold, ntot, ibatch
        real,    allocatable :: min_new(:), max_new(:), min_all(:), max_all(:)
        integer, allocatable :: bid_new(:), bid_all(:)
        
        if( self%nimagesets < 1 ) return

        ! Keep min/max arrays consistent if previous state was partially allocated.
        if( allocated(self%support_training_set%min_dim) .neqv. allocated(self%support_training_set%max_dim) .or. &
            allocated(self%support_training_set%min_dim) .neqv. allocated(self%support_training_set%batch_id) )then
            call self%kill_support_training_set()
        end if

        nkeep = 0
        do ii = 1, self%nimagesets
            if( .not. self%imagesets(ii)%is_rejected ) nkeep = nkeep + 1
        end do
        if( nkeep < 1 )then
            write(logfhandle,'(A)') 'cavg_compat: train_cache skipped (accepted=0)'
            return
        end if

        ibatch = self%support_training_set%nbatches + 1
        allocate(min_new(nkeep), max_new(nkeep), bid_new(nkeep))
        jj = 0
        do ii = 1, self%nimagesets
            if( self%imagesets(ii)%is_rejected ) cycle
            jj = jj + 1
            min_new(jj) = self%imagesets(ii)%feret_min
            max_new(jj) = self%imagesets(ii)%feret_max
            bid_new(jj) = ibatch
        end do

        if( .not. allocated(self%support_training_set%min_dim) )then
            allocate(self%support_training_set%min_dim(nkeep), self%support_training_set%max_dim(nkeep), self%support_training_set%batch_id(nkeep))
            self%support_training_set%min_dim = min_new
            self%support_training_set%max_dim = max_new
            self%support_training_set%batch_id = bid_new
            ntot = nkeep
        else
            nold = size(self%support_training_set%min_dim)
            ntot = nold + nkeep
            allocate(min_all(ntot), max_all(ntot), bid_all(ntot))
            min_all(1:nold) = self%support_training_set%min_dim
            max_all(1:nold) = self%support_training_set%max_dim
            bid_all(1:nold) = self%support_training_set%batch_id
            min_all(nold+1:ntot) = min_new
            max_all(nold+1:ntot) = max_new
            bid_all(nold+1:ntot) = bid_new
            call move_alloc(min_all, self%support_training_set%min_dim)
            call move_alloc(max_all, self%support_training_set%max_dim)
            call move_alloc(bid_all, self%support_training_set%batch_id)
        end if
        self%support_training_set%nbatches = ibatch

        deallocate(min_new, max_new, bid_new)

        write(logfhandle,'(A,I0,A,I0,A,I0)') 'cavg_compat: train_cache appended=', nkeep, ' run_total=', self%nimagesets, ' pool_total=', ntot
    end subroutine cache_support_training_set

    ! Run compatibility inference and write outputs.
    subroutine infer( self, spproj )
        class(cavg_compatibility_analysis), intent(inout) :: self
        type(sp_project),                      intent(in) :: spproj
        type(image),                          allocatable :: imgs(:), img_out(:), mask_out(:)
        type(cavg_quality_result)                         :: quality
        type(string)                                      :: out_imgs, out_masks
        integer                                           :: iimg, ncls, i, nrej, ncomp
        real                                              :: smpd

        if( .not. self%support_model%valid )then
            THROW_HARD('cavg_compatibility_analysis: infer requires successful train (support model is invalid)')
        end if

        call self%kill_imagesets()
        call self%spproj%kill()
        self%spproj = spproj
        call self%spproj%get_cavgs_stk(self%input_stkname, ncls, smpd)
        imgs            = read_cavgs_into_imgarr(self%spproj)
        self%nimagesets = size(imgs)
        if( self%nimagesets < 1 ) THROW_HARD('cavg_compatibility_analysis: no images found in input stack')
        allocate(self%imagesets(self%nimagesets), img_out(self%nimagesets), mask_out(self%nimagesets))
        call evaluate_cavg_quality_hard_reject(imgs, self%spproj%os_cls2D, 5000.0, quality)
        do iimg = 1, self%nimagesets
            if( quality%states(iimg) == 0 ) then
                self%imagesets(iimg)%is_rejected      = .true.
                self%imagesets(iimg)%rejection_reason = REJECT_REASON_QUALITY_HARD_REJECT
            end if
            call imgs(iimg)%set_smpd(smpd)
            call self%preprocess(iimg, imgs(iimg))
            call img_out(iimg)%copy(self%imagesets(iimg)%img)
            if( self%imagesets(iimg)%is_rejected )then
                call mask_out(iimg)%copy(self%imagesets(iimg)%img)
                call mask_out(iimg)%zero()
            else
                call mask_out(iimg)%copy(self%imagesets(iimg)%mask)
            end if
        end do

        out_imgs  = self%input_stkname//'.cavg_compatibility_imgs.mrcs'
        out_masks = self%input_stkname//'.cavg_compatibility_masks.mrcs'
        call write_imgarr(img_out,   out_imgs)
        call write_imgarr(mask_out, out_masks)

        call dealloc_imgarr(img_out)
        call dealloc_imgarr(mask_out)
        call dealloc_imgarr(imgs)

        write(logfhandle,'(A,I0)') 'cavg_compat: analyze n_images=', self%nimagesets

        call self%run_cluster_overfitting_rejection(5000.0)
        nrej = 0
        do i = 1, self%nimagesets
            if( self%imagesets(i)%is_rejected ) nrej = nrej + 1
        end do
        write(logfhandle,'(A,I0)') 'cavg_compat: reject_cluster n_rejected=', nrej
        write(logfhandle,'(A,ES14.6,A,ES14.6,A,ES14.6)') 'cavg_compat: support_model axes c/b/a=', &
            self%support_model%axis_c, '/', self%support_model%axis_b, '/', self%support_model%axis_a
        write(logfhandle,'(A,L1,A,ES14.6,A,ES14.6,A,ES14.6)') 'cavg_compat: support_model autotune=', self%support_model%do_autotune, &
            ' relax=', self%support_model%used_relax, ' qlo=', self%support_model%used_qlo, ' qhi=', self%support_model%used_qhi

        call apply_support_model(self, ncomp)
        write(logfhandle,'(A,I0)') 'cavg_compat: support_apply n_compatible=', ncomp
        nrej = 0
        do i = 1, self%nimagesets
            if( self%imagesets(i)%is_rejected ) nrej = nrej + 1
        end do
        write(logfhandle,'(A,I0)') 'cavg_compat: support_apply n_rejected=', nrej

        call self%write_selected_rejected_stacks()

        call self%print_rejection_reasons()
        write(logfhandle,'(A)') 'cavg_compat: wrote rejection reasons file=rejection_reasons.txt'
    end subroutine infer

    ! Train support model from class averages in a project.
    subroutine train_1( self, spproj )
        class(cavg_compatibility_analysis), intent(inout) :: self
        type(sp_project),                      intent(in) :: spproj
        type(image),                          allocatable :: imgs(:)
        type(cavg_quality_result)                         :: quality
        integer                                           :: iimg, ncls, i, nrej
        real                                              :: smpd
        call self%kill_imagesets()
        call self%spproj%kill()
        self%spproj = spproj
        call self%spproj%get_cavgs_stk(self%input_stkname, ncls, smpd)
        imgs            = read_cavgs_into_imgarr(self%spproj)
        self%nimagesets = size(imgs)
        if( self%nimagesets < 1 ) THROW_HARD('cavg_compatibility_analysis: no images found in input stack')
        allocate(self%imagesets(self%nimagesets))
        call evaluate_cavg_quality_hard_reject(imgs, self%spproj%os_cls2D, 5000.0, quality)
        do iimg = 1, self%nimagesets
            if( quality%states(iimg) == 0 ) then
                self%imagesets(iimg)%is_rejected      = .true.
                self%imagesets(iimg)%rejection_reason = REJECT_REASON_QUALITY_HARD_REJECT
            end if
            call imgs(iimg)%set_smpd(smpd)
            call self%preprocess(iimg, imgs(iimg))
        end do
        call dealloc_imgarr(imgs)

        write(logfhandle,'(A,I0)') 'cavg_compat: analyze n_images=', self%nimagesets

        call self%run_cluster_overfitting_rejection(5000.0)
        nrej = 0
        do i = 1, self%nimagesets
            if( self%imagesets(i)%is_rejected ) nrej = nrej + 1
        end do
        write(logfhandle,'(A,I0)') 'cavg_compat: reject_cluster n_rejected=', nrej
        call self%cache_support_training_set()
        call self%generate_support_model()
    end subroutine train_1

    ! Train support model from class averages in a stack.
    subroutine train_2( self, refs )
        class(cavg_compatibility_analysis), intent(inout) :: self
        type(string),                          intent(in) :: refs
        type(image),                          allocatable :: imgs(:)
        integer                                           :: iimg
        call self%kill_imagesets()
        call self%spproj%kill()
        imgs            = read_stk_into_imgarr(refs)
        self%nimagesets = size(imgs)
        if( self%nimagesets < 1 ) THROW_HARD('cavg_compatibility_analysis: no images found in input stack')
        allocate(self%imagesets(self%nimagesets))
        do iimg = 1, self%nimagesets
            call self%preprocess(iimg, imgs(iimg))
        end do
        call dealloc_imgarr(imgs)

        write(logfhandle,'(A,I0)') 'cavg_compat: analyze n_images=', self%nimagesets

        call self%cache_support_training_set()
        call self%generate_support_model()
    end subroutine train_2

    ! Export selection states as integers (1=accepted, 0=rejected).
    subroutine get_rejection_states( self, states )
        class(cavg_compatibility_analysis), intent(in)  :: self
        integer, allocatable,               intent(out) :: states(:)
        integer :: i

        if( allocated(states) ) deallocate(states)

        if( self%nimagesets < 1 )then
            allocate(states(0))
            return
        end if

        allocate(states(self%nimagesets), source=1)
        do i = 1, self%nimagesets
            if( self%imagesets(i)%is_rejected ) states(i) = 0
        end do
    end subroutine get_rejection_states

    ! Write selected and rejected class-average stacks and JPEG summaries.
    subroutine write_selected_rejected_stacks( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        type(image), allocatable :: selected_imgs(:), rejected_imgs(:)
        type(string)             :: selected_out, rejected_out
        type(string)             :: selected_jpg, rejected_jpg
        integer                  :: i, nsel, nrej, isel, irej

        nsel = 0
        nrej = 0
        do i = 1, self%nimagesets
            if( .not. self%imagesets(i)%is_rejected ) nsel = nsel + 1
            if( self%imagesets(i)%is_rejected ) nrej = nrej + 1
        end do

        selected_out =  swap_suffix(self%input_stkname, '_selected'//MRC_EXT, MRC_EXT)
        rejected_out =  swap_suffix(self%input_stkname, '_rejected'//MRC_EXT, MRC_EXT)

        if( nsel > 0 )then
            allocate(selected_imgs(nsel))
            isel = 0
            do i = 1, self%nimagesets
                if( self%imagesets(i)%is_rejected ) cycle
                isel = isel + 1
                call selected_imgs(isel)%copy(self%imagesets(i)%img)
            end do
            call write_imgarr(selected_imgs, selected_out)
            call dealloc_imgarr(selected_imgs)
        end if

        if( nrej > 0 )then
            allocate(rejected_imgs(nrej))
            irej = 0
            do i = 1, self%nimagesets
                if( .not. self%imagesets(i)%is_rejected ) cycle
                irej = irej + 1
                call rejected_imgs(irej)%copy(self%imagesets(i)%img)
            end do
            call write_imgarr(rejected_imgs, rejected_out)
            call dealloc_imgarr(rejected_imgs)
        end if

        write(logfhandle,'(A,I0,A,A)') 'cavg_compat: wrote stack kind=selected n=', nsel, ' file=', selected_out%to_char()
        write(logfhandle,'(A,I0,A,A)') 'cavg_compat: wrote stack kind=rejected n=', nrej, ' file=', rejected_out%to_char()

        selected_jpg = swap_suffix(selected_out, JPG_EXT, MRC_EXT)
        rejected_jpg = swap_suffix(rejected_out, JPG_EXT, MRC_EXT)
        call mrc2jpeg_tiled(selected_out, selected_jpg)
        call mrc2jpeg_tiled(rejected_out, rejected_jpg)
        write(logfhandle,'(A,A)') '>>> JPEG ', selected_jpg%to_char()
        write(logfhandle,'(A,A)') '>>> JPEG ', rejected_jpg%to_char()
    end subroutine write_selected_rejected_stacks


    ! Fit latent support-model axes from cached training dimensions.
    subroutine generate_support_model( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        integer :: nn, i, j, k, best_count
        integer :: ir, iq, jq, final_count, best_final_count
        real, allocatable :: cands(:), mins_sup(:), maxs_sup(:), w_sup(:), sample_weights(:)
        integer, allocatable :: batch_counts(:)
        integer :: nbatches
        real    :: b_cand, left, right
        real    :: spread, best_spread, score, best_score
        real    :: weighted_count, best_weighted_count, best_final_weighted_count
        real    :: relax_cur, qlo_cur, qhi_cur
        real    :: c_cur, b_cur, a_cur, c_soft, a_soft
        real    :: c_best, b_best, a_best
        logical :: in_b, in_c, in_a
        logical, allocatable :: tmp_support(:), support_cur(:)
        
        call self%kill_support_model()

        if( allocated(self%support_training_set%min_dim) .neqv. allocated(self%support_training_set%max_dim) )then
            call self%kill_support_training_set()
            write(logfhandle,'(A)') 'cavg_compat: support_model invalid training allocation state (min/max mismatch)'
            return
        end if

        if( allocated(self%support_training_set%min_dim) .and. .not. allocated(self%support_training_set%batch_id) )then
            nn = size(self%support_training_set%min_dim)
            allocate(self%support_training_set%batch_id(nn), source=1)
            self%support_training_set%nbatches = 1
            write(logfhandle,'(A)') 'cavg_compat: support_model reconstructed batch metadata from legacy training set'
        else if( .not. allocated(self%support_training_set%min_dim) .and. allocated(self%support_training_set%batch_id) )then
            call self%kill_support_training_set()
            write(logfhandle,'(A)') 'cavg_compat: support_model invalid training allocation state (batch without dims)'
            return
        end if

        if( .not. allocated(self%support_training_set%min_dim) )then
            write(logfhandle,'(A)') 'cavg_compat: support_model training set is not allocated'
            return
        end if

        if( size(self%support_training_set%min_dim) /= size(self%support_training_set%max_dim) .or. &
            size(self%support_training_set%min_dim) /= size(self%support_training_set%batch_id) )then
            call self%kill_support_training_set()
            write(logfhandle,'(A)') 'cavg_compat: support_model training set size mismatch'
            return
        end if

        nn = size(self%support_training_set%min_dim)

        if( nn <= 2 )then
            write(logfhandle,'(A)') 'cavg_compat: support_model skipped (training set too small)'
            return
        end if

        nbatches = self%support_training_set%nbatches
        if( nbatches < 1 ) nbatches = maxval(self%support_training_set%batch_id)
        if( nbatches < 1 )then
            write(logfhandle,'(A)') 'cavg_compat: support_model invalid batch metadata'
            return
        end if

        allocate(batch_counts(nbatches), source=0)
        allocate(sample_weights(nn),     source=0.0)
        do i = 1, nn
            if( self%support_training_set%batch_id(i) < 1 .or. self%support_training_set%batch_id(i) > nbatches )then
                deallocate(batch_counts, sample_weights)
                write(logfhandle,'(A)') 'cavg_compat: support_model invalid batch id in training set'
                return
            end if
            batch_counts(self%support_training_set%batch_id(i)) = batch_counts(self%support_training_set%batch_id(i)) + 1
        end do
        do i = 1, nn
            sample_weights(i) = 1.0 / real(batch_counts(self%support_training_set%batch_id(i)))
        end do

        allocate(tmp_support(nn), source=.false.)
        allocate(support_cur(nn), source=.false.)
        allocate(cands(2*nn))
        k = 0
        do i = 1, nn
            k = k + 1
            cands(k) = self%support_training_set%min_dim(i)
            k = k + 1
            cands(k) = self%support_training_set%max_dim(i)
        end do

        best_count = -1
        best_weighted_count = -1.0
        best_spread = huge(1.0)
        best_final_count = -1
        best_final_weighted_count = -1.0
        best_score = -huge(1.0)
        c_best = 0.0
        b_best = 0.0
        a_best = 0.0

        do ir = 1, NRELAX
            if( self%support_model%do_autotune )then
                relax_cur = RELAX_GRID(ir)
            else
                if( ir > 1 ) cycle
                relax_cur = 0.07
            end if

            do iq = 1, NQ
                if( self%support_model%do_autotune )then
                    qlo_cur = QLOW_GRID(iq)
                else
                    if( iq > 1 ) cycle
                    qlo_cur = 0.08
                end if

                do jq = 1, NQ
                    if( self%support_model%do_autotune )then
                        qhi_cur = QHIGH_GRID(jq)
                    else
                        if( jq > 1 ) cycle
                        qhi_cur = 0.93
                    end if
                    if( qhi_cur <= qlo_cur ) cycle

                    best_count = -1
                    best_spread = huge(1.0)
                    b_cur = 0.0
                    support_cur = .false.

                    do i = 1, k
                        b_cand = cands(i)
                        final_count = 0
                        weighted_count = 0.0
                        spread = 0.0
                        tmp_support = .false.

                        do j = 1, nn
                            left = self%support_training_set%min_dim(j) * (1.0 - relax_cur)
                            right = self%support_training_set%max_dim(j) * (1.0 + relax_cur)
                            if( b_cand >= left .and. b_cand <= right )then
                                final_count = final_count + 1
                                weighted_count = weighted_count + sample_weights(j)
                                tmp_support(j) = .true.
                                spread = spread + sample_weights(j) * abs(b_cand - 0.5 * (self%support_training_set%min_dim(j) + self%support_training_set%max_dim(j))) / &
                                    max(self%support_training_set%max_dim(j)-self%support_training_set%min_dim(j), 1.0e-6)
                            end if
                        end do

                        if( weighted_count > best_weighted_count .or. &
                            (abs(weighted_count - best_weighted_count) <= 1.0e-6 .and. spread < best_spread) )then
                            best_count = final_count
                            best_weighted_count = weighted_count
                            best_spread = spread
                            b_cur = b_cand
                            support_cur = tmp_support
                        end if
                    end do

                    if( best_count <= 0 ) cycle

                    allocate(mins_sup(best_count), maxs_sup(best_count), w_sup(best_count))
                    j = 0
                    do i = 1, nn
                        if( .not. support_cur(i) ) cycle
                        j = j + 1
                        mins_sup(j) = self%support_training_set%min_dim(i)
                        maxs_sup(j) = self%support_training_set%max_dim(i)
                        w_sup(j) = sample_weights(i)
                    end do

                    c_cur = percentile_weighted_real(mins_sup, w_sup, qlo_cur)
                    a_cur = percentile_weighted_real(maxs_sup, w_sup, qhi_cur)
                    if( c_cur > b_cur ) c_cur = minval(mins_sup)
                    if( a_cur < b_cur ) a_cur = maxval(maxs_sup)

                    final_count = 0
                    weighted_count = 0.0
                    spread = 0.0
                    tmp_support = .false.
                    do i = 1, nn
                        left = self%support_training_set%min_dim(i) * (1.0 - relax_cur)
                        right = self%support_training_set%max_dim(i) * (1.0 + relax_cur)
                        c_soft = c_cur * (1.0 - relax_cur - SUPPORT_EDGE_SOFT_FRAC)
                        a_soft = a_cur * (1.0 + relax_cur + SUPPORT_EDGE_SOFT_FRAC)
                        in_b = (b_cur >= left .and. b_cur <= right)
                        in_c = (self%support_training_set%min_dim(i) >= c_cur * (1.0 - relax_cur))
                        in_a = (self%support_training_set%max_dim(i) <= a_cur * (1.0 + relax_cur))
                        if( in_b .and. ( (in_c .and. in_a) .or. (self%support_training_set%min_dim(i) >= c_soft .and. in_a) .or. (in_c .and. self%support_training_set%max_dim(i) <= a_soft) ) )then
                            final_count = final_count + 1
                            weighted_count = weighted_count + sample_weights(i)
                            tmp_support(i) = .true.
                            spread = spread + sample_weights(i) * abs(b_cur - 0.5 * (self%support_training_set%min_dim(i) + self%support_training_set%max_dim(i))) / &
                                max(self%support_training_set%max_dim(i)-self%support_training_set%min_dim(i), 1.0e-6)
                        end if
                    end do

                    score = weighted_count - 0.01 * spread - 0.10 * relax_cur
                    if( score > best_score .or. (abs(score - best_score) <= 1.0e-6 .and. weighted_count > best_final_weighted_count) )then
                        best_score = score
                        best_final_count = final_count
                        best_final_weighted_count = weighted_count
                        c_best = c_cur
                        b_best = b_cur
                        a_best = a_cur
                        self%support_model%used_relax = relax_cur
                        self%support_model%used_qlo = qlo_cur
                        self%support_model%used_qhi = qhi_cur
                    end if

                    deallocate(mins_sup, maxs_sup, w_sup)
                end do
            end do
        end do

        if( best_final_count <= 0 )then
            deallocate(cands, tmp_support, support_cur, batch_counts, sample_weights)
            return
        end if

        self%support_model%axis_c = c_best
        self%support_model%axis_b = b_best
        self%support_model%axis_a = a_best
        self%support_model%valid = .true.

        deallocate(cands, tmp_support, support_cur, batch_counts, sample_weights)

    contains

        ! Return the nearest-rank percentile of a real-valued vector.
        real function percentile_real(arr, q) result(pval)
            real, intent(in) :: arr(:)
            real, intent(in) :: q
            real, allocatable :: sorted(:)
            integer :: nvals, idx

            nvals = size(arr)
            if( nvals <= 0 )then
                pval = 0.0
                return
            end if

            allocate(sorted(nvals))
            sorted = arr
            call hpsort(sorted)

            idx = 1 + int(q * real(max(nvals - 1, 0)))
            idx = max(1, min(nvals, idx))
            pval = sorted(idx)
            deallocate(sorted)
        end function percentile_real

        ! Return the weighted nearest-rank percentile of a real-valued vector.
        real function percentile_weighted_real(arr, warr, q) result(pval)
            real, intent(in) :: arr(:)
            real, intent(in) :: warr(:)
            real, intent(in) :: q
            real, allocatable :: sorted(:), wsorted(:)
            real :: wtot, wacc, qtarget, wtmp
            integer :: nvals, ia, ib, kk

            nvals = size(arr)
            if( nvals <= 0 )then
                pval = 0.0
                return
            end if
            if( size(warr) /= nvals )then
                pval = percentile_real(arr, q)
                return
            end if

            allocate(sorted(nvals), wsorted(nvals))
            sorted = arr
            wsorted = warr

            do ia = 1, nvals - 1
                do ib = ia + 1, nvals
                    if( sorted(ib) < sorted(ia) )then
                        pval = sorted(ia)
                        sorted(ia) = sorted(ib)
                        sorted(ib) = pval
                        wtmp = wsorted(ia)
                        wsorted(ia) = wsorted(ib)
                        wsorted(ib) = wtmp
                    end if
                end do
            end do

            wtot = sum(max(wsorted, 0.0))
            if( wtot <= 0.0 )then
                pval = percentile_real(arr, q)
                deallocate(sorted, wsorted)
                return
            end if

            qtarget = min(max(q, 0.0), 1.0) * wtot
            wacc = 0.0
            pval = sorted(nvals)
            do kk = 1, nvals
                wacc = wacc + max(wsorted(kk), 0.0)
                if( wacc >= qtarget )then
                    pval = sorted(kk)
                    exit
                end if
            end do

            deallocate(sorted, wsorted)
        end function percentile_weighted_real
        
    end subroutine generate_support_model

    ! Apply support-model bounds and mark size-based rejections.
    subroutine apply_support_model(self, ncomp)
        class(cavg_compatibility_analysis), intent(inout) :: self
        integer,                            intent(out)   :: ncomp
        integer :: ii
        logical :: compatible
        real    :: min_dim, max_dim, relax

        ncomp = 0
        relax = self%support_model%used_relax
        do ii = 1, self%nimagesets
            if( self%imagesets(ii)%is_rejected )then
                cycle
            end if

            min_dim = self%imagesets(ii)%feret_min
            max_dim = self%imagesets(ii)%feret_max
            compatible = self%support_model%axis_b >= min_dim * (1.0 - relax) .and. &
                self%support_model%axis_b <= max_dim * (1.0 + relax) .and. &
                min_dim >= self%support_model%axis_c * (1.0 - relax) .and. &
                max_dim <= self%support_model%axis_a * (1.0 + relax)

            if( .not. compatible )then
                if( min_dim >= self%support_model%axis_c * (1.0 - relax - RESCUE_EDGE_FRAC) .and. &
                    max_dim <= self%support_model%axis_a * (1.0 + relax + RESCUE_EDGE_FRAC) )then
                    compatible = .true.
                    write(logfhandle,'(A,I0,A)') 'cavg_compat: support_apply rescued idx=', ii, ' reason=within_c_to_a_envelope'
                end if
            end if

            if( compatible ) ncomp = ncomp + 1
            if( .not. compatible )then
                self%imagesets(ii)%is_rejected = .true.
                self%imagesets(ii)%rejection_reason = REJECT_REASON_SIZE_INCOMPATIBLE
                write(logfhandle,'(A,I0,A)') 'cavg_compat: support_apply rejected idx=', ii, ' reason=size_incompatible_subset'
            end if
        end do
    end subroutine apply_support_model

    ! Write rejection reason codes for rejected image sets.
    subroutine print_rejection_reasons( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        integer :: i, funit
        character(len=64) :: reason_txt

        open(newunit=funit, file='rejection_reasons.txt', status='replace', action='write')
        write(funit,'(A)') '# class_index reason_code reason_text'
        write(logfhandle,'(A)') 'cavg_compat: writing rejection reasons'

        do i = 1, self%nimagesets
            if( .not. self%imagesets(i)%is_rejected ) cycle
            select case(self%imagesets(i)%rejection_reason)
            case(REJECT_REASON_ZERO_VARIANCE)
                reason_txt = 'zero_variance'
            case(REJECT_REASON_MASK_OUTSIDE_SUPPORT)
                reason_txt = 'mask_outside_circular_support'
            case(REJECT_REASON_SIZE_INCOMPATIBLE)
                reason_txt = 'size_incompatible_subset'
            case(REJECT_REASON_SUSPECTED_OVERFITTING)
                reason_txt = 'suspected_overfitting'
            case(REJECT_REASON_QUALITY_HARD_REJECT)
                reason_txt = 'quality_hard_reject'
            case default
                reason_txt = 'unknown'
            end select
            write(funit,'(I0,1X,I0,1X,A)') i, self%imagesets(i)%rejection_reason, trim(reason_txt)
            write(logfhandle,'(A,I0,A,I0,A,A)') 'cavg_compat: reject idx=', i, ' code=', self%imagesets(i)%rejection_reason, ' reason=', trim(reason_txt)
        end do

        close(funit)
    end subroutine print_rejection_reasons

    ! Preprocess one class average and populate image-set descriptors.
    subroutine preprocess( self, i, img )
        class(cavg_compatibility_analysis), intent(inout) :: self
        type(image),                        intent(inout) :: img
        integer,                            intent(in)    :: i
        logical,                            allocatable   :: l_circ(:,:,:)
        real(kind=c_float),                 pointer       :: rmat_mask(:,:,:) => null()
        type(image)     :: circ_mask
        type(image_bin) :: cc_img
        integer, allocatable :: cc_sizes(:)
        integer :: largest_cc
        integer     :: ldim(3)
        integer     :: ldim_target(3)
        integer     :: imorph

        ! Preprocess the input image: zero edge average and bandpass filter.
        call img%zero_edgeavg()
        call img%bp(0., 10.)

        ! Frequency-domain resize to 128x128: FFT -> crop/pad -> inverse FFT.
        call img%fft()
        ldim        = img%get_ldim()
        ldim_target = [ANALYSIS_BOXSIZE, ANALYSIS_BOXSIZE, ldim(3)]
        if( ldim(1) > ldim_target(1) .or. ldim(2) > ldim_target(2) )then
            call img%clip_inplace(ldim_target)
        else if( ldim(1) < ldim_target(1) .or. ldim(2) < ldim_target(2) )then
            call img%pad_inplace(ldim_target)
        end if
        call img%ifft()
        call img%set_smpd(img%get_smpd() * real(ldim(1)) / real(ldim_target(1)))

        ! Keep a processed copy of the resized image for downstream analysis/output.
        self%imagesets(i)%img = img

        if( self%imagesets(i)%is_rejected ) return

        ! Reset per-image outputs so repeated calls on the same slot are safe.
        self%imagesets(i)%rejection_reason = REJECT_REASON_NONE
        self%imagesets(i)%local_var_in_mask = 0.0
        self%imagesets(i)%local_var_out_mask = 0.0

        ! Compute global variance and reject if zero.
        self%imagesets(i)%variance      = img%variance()
        self%imagesets(i)%is_rejected   = self%imagesets(i)%variance == 0.0
        if( self%imagesets(i)%is_rejected ) then
            self%imagesets(i)%rejection_reason = REJECT_REASON_ZERO_VARIANCE
            return
        end if
        ! Compute Otsu threshold and binary mask.
        call otsu_img(img)
        ! Apply 5 px morphological closing to the binary Otsu mask.
        call self%imagesets(i)%mask%transfer2bimg(img)
        do imorph = 1, ANALYSIS_MORPH_SIZE
            call self%imagesets(i)%mask%dilate()
        end do
        do imorph = 1, ANALYSIS_MORPH_SIZE
            call self%imagesets(i)%mask%erode()
        end do
        ! Keep only the largest connected foreground component in the mask.
        call self%imagesets(i)%mask%find_ccs(cc_img)
        cc_sizes = cc_img%size_ccs()
        if( size(cc_sizes) > 0 .and. maxval(cc_sizes) > 0 )then
            largest_cc = maxloc(cc_sizes, dim=1)
            call cc_img%cc2bin(largest_cc)
            call self%imagesets(i)%mask%copy_bimg(cc_img)
        end if
        if( allocated(cc_sizes) ) deallocate(cc_sizes)
        call cc_img%kill_bimg()
        ! Build a hard circular support mask with diameter equal to the box size.
        ldim = self%imagesets(i)%mask%get_ldim()
        call circ_mask%disc(ldim, self%imagesets(i)%mask%get_smpd(), 0.5 * real(min(ldim(1), ldim(2))), l_circ)
        call self%imagesets(i)%mask%get_rmat_ptr(rmat_mask)
        ! Any foreground outside support marks this class average as rejected.
        if( any(rmat_mask(1:ldim(1),1:ldim(2),1:ldim(3)) > 0.5 .and. .not. l_circ(1:ldim(1),1:ldim(2),1:ldim(3))) ) then
            self%imagesets(i)%is_rejected = .true.
            self%imagesets(i)%rejection_reason = REJECT_REASON_MASK_OUTSIDE_SUPPORT
            ! Early exit: free temporary support resources before returning.
            nullify(rmat_mask)
            if( allocated(l_circ) ) deallocate(l_circ)
            call circ_mask%kill()
            return
        end if
        if( allocated(l_circ) ) deallocate(l_circ)
        call circ_mask%kill()
        ! Compute Feret diameters of the mask and store them in the pointset.
        call self%imagesets(i)%mask%feret_minmax(self%imagesets(i)%feret_min, self%imagesets(i)%feret_max)
        ! Compute local variance inside and outside the mask.
        call self%imagesets(i)%img%loc_var_masked(rmat_mask(1:ldim(1), 1:ldim(2), 1), 10, &
            self%imagesets(i)%local_var_in_mask, self%imagesets(i)%local_var_out_mask)
        write(logfhandle,'(A,I0,A,ES14.6,A,ES14.6)') 'cavg_compat: preprocess idx=', i, ' var_in=', self%imagesets(i)%local_var_in_mask, &
            ' var_out=', self%imagesets(i)%local_var_out_mask
        nullify(rmat_mask)
        
    end subroutine preprocess

    ! Reject clusters with globally low local-variance statistics.
   subroutine run_cluster_overfitting_rejection( self, mskdiam, nclust_in )
        class(cavg_compatibility_analysis), intent(inout) :: self
        real,                                  intent(in) :: mskdiam
        integer,                     optional, intent(in) :: nclust_in
        integer,                              allocatable :: states(:), clusters(:)
        integer,                              allocatable :: per_cluster_lowvar(:), per_cluster_total(:)
        logical,                              allocatable :: reject_cluster(:)
        real,                                 allocatable :: vals_in(:), vals_out(:), absdev(:)
        type(string)                                      :: projname
        type(cmdline)                                     :: cline
        type(sp_project)                                  :: spproj
        type(commander_cluster_cavgs)                     :: commander
        integer                                           :: i, nclustered, nclusters, k, iclust
        real                                              :: thr_in, thr_out, med_in_all, med_out_all, mad_in_all, mad_out_all
        real                                              :: vin, vout, frac_lowvar
        projname = 'cluster_cavgs'//METADATA_EXT

        allocate(states(self%nimagesets), source=1)
        do i = 1, self%nimagesets
            if( self%imagesets(i)%is_rejected ) states(i) = 0
        end do
        call self%spproj%map_cavgs_selection(states)
        deallocate(states)
        call self%spproj%write(projname)
        call cline%set('projfile', projname)
        call cline%set('mskdiam', mskdiam)
        if( present(nclust_in) ) call cline%set('ncls', nclust_in)
        call commander%execute(cline)
        call spproj%read(projname)
        clusters = spproj%os_cls2D%get_all_asint('cluster')
        if( size(clusters) == self%nimagesets )then
            nclustered = count(clusters > 0)
            nclusters  = maxval(clusters)
            write(logfhandle,'(A,I0,A,I0)') 'run_exec_cluster_cavgs: clustered classes=', nclustered, ' total=', self%nimagesets

            k = 0
            do i = 1, self%nimagesets
                if( self%imagesets(i)%is_rejected ) cycle
                k = k + 1
            end do
            if( k > 0 .and. nclusters > 0 )then
                allocate(vals_in(k), vals_out(k))
                k = 0
                do i = 1, self%nimagesets
                    if( self%imagesets(i)%is_rejected ) cycle
                    k = k + 1
                    vals_in(k)  = log(max(self%imagesets(i)%local_var_in_mask,  1.0e-8))
                    vals_out(k) = log(max(self%imagesets(i)%local_var_out_mask, 1.0e-8))
                end do

                med_in_all = median_real_local(vals_in)
                med_out_all = median_real_local(vals_out)
                allocate(absdev(size(vals_in)))
                absdev = abs(vals_in - med_in_all)
                mad_in_all = median_real_local(absdev)
                deallocate(absdev)
                allocate(absdev(size(vals_out)))
                absdev = abs(vals_out - med_out_all)
                mad_out_all = median_real_local(absdev)
                deallocate(absdev)

                thr_in  = med_in_all  - 2.0 * max(mad_in_all,  1.0e-3)
                thr_out = med_out_all - 0.25 * max(mad_out_all, 1.0e-3)

                allocate(per_cluster_lowvar(nclusters), source=0)
                allocate(per_cluster_total(nclusters), source=0)
                allocate(reject_cluster(nclusters), source=.false.)
                do i = 1, self%nimagesets
                    if( clusters(i) < 1 .or. clusters(i) > nclusters ) cycle
                    per_cluster_total(clusters(i)) = per_cluster_total(clusters(i)) + 1
                    if( self%imagesets(i)%is_rejected ) cycle
                    vin = log(max(self%imagesets(i)%local_var_in_mask,   1.0e-8))
                    vout = log(max(self%imagesets(i)%local_var_out_mask, 1.0e-8))
                    write(logfhandle,'(A,I0,A,I0,A,ES14.6,A,ES14.6)') 'run_exec_cluster_cavgs: idx=', i, ' cluster=', clusters(i), &
                        ' var_in=', vin, ' var_out=', vout
                    if( vin <= thr_in .and. vout <= thr_out )then
                        per_cluster_lowvar(clusters(i)) = per_cluster_lowvar(clusters(i)) + 1
                    end if
                end do

                write(logfhandle,'(A,ES14.6,A,ES14.6)') 'run_exec_cluster_cavgs: localvar thresholds in/out=', thr_in, ' / ', thr_out
                do iclust = 1, nclusters
                    if( per_cluster_total(iclust) > 0 )then
                        frac_lowvar = real(per_cluster_lowvar(iclust)) / real(per_cluster_total(iclust))
                    else
                        frac_lowvar = 0.0
                    end if
                    if( frac_lowvar > LOWVAR_CLUSTER_REJECT_FRAC ) reject_cluster(iclust) = .true.
                    write(logfhandle,'(A,I0,A,I0,A,I0,A,F6.2,A,L1)') 'run_exec_cluster_cavgs: cluster=', iclust, &
                        ' lowvar_in_and_out=', per_cluster_lowvar(iclust), ' total=', per_cluster_total(iclust), &
                        ' lowvar_pct=', 100.0 * frac_lowvar, ' reject_cluster=', reject_cluster(iclust)
                end do

                do i = 1, self%nimagesets
                    if( clusters(i) >= 1 .and. clusters(i) <= nclusters )then
                        if( reject_cluster(clusters(i)) )then
                            self%imagesets(i)%is_rejected = .true.
                            if( self%imagesets(i)%rejection_reason == REJECT_REASON_NONE ) then
                                self%imagesets(i)%rejection_reason = REJECT_REASON_SUSPECTED_OVERFITTING
                            end if
                        end if
                    end if
                end do
                deallocate(reject_cluster, per_cluster_lowvar, per_cluster_total, vals_in, vals_out)
            end if
            
        else
            write(logfhandle,'(A,I0,A,I0)') 'run_exec_cluster_cavgs: cluster size mismatch cls2D=', size(clusters), ' imagesets=', self%nimagesets
        end if
        if( allocated(clusters) ) deallocate(clusters)
        
        call spproj%kill()
        call cline%kill()

    contains

        function median_real_local(arr) result(med_val)
            real, intent(in) :: arr(:)
            real :: med_val
            real, allocatable :: sorted(:)
            integer :: n, mid
            if( size(arr) == 0 )then
                med_val = 0.0
                return
            end if
            n = size(arr)
            mid = (n + 1) / 2
            allocate(sorted(n))
            sorted = arr
            call hpsort(sorted)
            if( mod(n,2) == 1 )then
                med_val = sorted(mid)
            else
                med_val = 0.5 * (sorted(mid) + sorted(mid+1))
            end if
            deallocate(sorted)
        end function median_real_local

    end subroutine run_cluster_overfitting_rejection

end module simple_cavg_compatibility_analysis
