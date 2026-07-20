!@descr: class-average compatibility analysis, support-model fitting, and size-based rejection
!@header: Builds a size-based support model from class averages, then applies
!         fitted Feret-axis bounds (c, b, a) to reject incompatible classes.
!         Includes convergence and metric reporting helpers for monitoring.

module simple_class_compatibility

use simple_defs,                  only: logfhandle
use simple_image,                 only: image
use simple_error,                 only: simple_exception
use simple_string,                only: string
use simple_image_bin,             only: image_bin
use simple_sp_project,            only: sp_project
use simple_segmentation,          only: otsu_img
use simple_string_utils,          only: int2str
use simple_imgarr_utils,          only: read_cavgs_into_imgarr, read_stk_into_imgarr, dealloc_imgarr
use simple_srch_sort_loc,         only: hpsort

implicit none

public :: class_compatibility, support_model_metrics, PREPROCESS_MORPH_SIZE

private
#include "simple_local_flags.inc"

! Preprocessing constants for class-average compatibility analysis.
integer, parameter :: PREPROCESS_BOXSIZE    = 128    ! Working box size for resized class averages (pixels)
integer, parameter :: PREPROCESS_MORPH_SIZE = 5      ! Number of dilate/erode passes for morphological closing

! Support-model search and convergence constants.
integer, parameter :: NRELAX                     = 5                              ! Number of relaxation candidates in support-model grid
integer, parameter :: NQ                         = 3                              ! Number of lower/upper quantile candidates in grid
real,    parameter :: DEFAULT_RELAX              = 0.07                           ! Default relative interval expansion for c/a axis estimation
real,    parameter :: DEFAULT_QLO                = 0.08                           ! Default quantile candidates for c/a axis estimation
real,    parameter :: DEFAULT_QHI                = 0.93                           ! Default quantile candidates for c/a axis estimation
real,    parameter :: CONV_ABS_EPS               = 0.5                            ! Absolute tolerance (pixels) for axis stability
real,    parameter :: CONV_REL_EPS               = 0.01                           ! Relative tolerance for axis stability
real,    parameter :: RESCUE_EDGE_FRAC           = 0.00                           ! Extra margin for boundary rescue in size compatibility
real,    parameter :: SUPPORT_EDGE_SOFT_FRAC_MAX = 0.10                           ! Initial soft edge slack (10%) for c/a bounds during support check
real,    parameter :: RELAX_GRID(NRELAX)         = [0.03, 0.05, 0.07, 0.10, 0.15] ! Relative interval expansion candidates
real,    parameter :: QLOW_GRID(NQ)              = [0.05, 0.10, 0.15]             ! Lower-quantile candidates for c-axis estimate
real,    parameter :: QHIGH_GRID(NQ)             = [0.85, 0.90, 0.95]             ! Upper-quantile candidates for a-axis estimate

type :: support_model_training_set
    ! Per-training-batch Feret bounds extracted from selected class averages.
    real,    allocatable :: min_dim(:)
    real,    allocatable :: max_dim(:)
end type support_model_training_set

type :: support_model_metrics
    ! Snapshot of current fitted axes and their latest update deltas.
    real    :: axis_c      = 0.0
    real    :: axis_b      = 0.0
    real    :: axis_a      = 0.0
    real    :: delta_c     = 0.0
    real    :: delta_b     = 0.0
    real    :: delta_a     = 0.0
    logical :: delta_valid = .false.
    logical :: valid       = .false.
    logical :: converged   = .false.
end type support_model_metrics

type :: support_model
    ! History of appended training batches; each batch is reweighted equally.
    type(support_model_training_set), allocatable :: training_set(:)
    logical                                       :: valid          = .false.
    logical                                       :: converged      = .false.
    logical                                       :: delta_valid    = .false.
    real                                          :: axis_c         = 0.0
    real                                          :: axis_b         = 0.0
    real                                          :: axis_a         = 0.0
    real                                          :: delta_c        = 0.0
    real                                          :: delta_b        = 0.0
    real                                          :: delta_a        = 0.0
    real                                          :: used_relax     = DEFAULT_RELAX
    real                                          :: used_qlo       = DEFAULT_QLO
    real                                          :: used_qhi       = DEFAULT_QHI
    real                                          :: used_soft_frac = SUPPORT_EDGE_SOFT_FRAC_MAX
end type support_model

type :: class_compatibility
    private
    type(support_model) :: support_model        ! Internal support model used by infer() via apply_support_model().

contains
    procedure, public  :: new
    procedure, public  :: kill
    procedure, public  :: kill_support_model
    procedure, public  :: train_1
    procedure, public  :: train_2
    procedure, public  :: infer
    procedure, public  :: converged
    procedure, public  :: get_support_model_metrics
    procedure, private :: apply_support_model
    procedure, private :: update_support_model
    procedure, private :: new_support_training_set
    procedure, private :: preprocess_img
    generic,   public  :: train => train_1, train_2
end type class_compatibility

contains

    ! ----------------------------------------------------------------
    ! lifecycle management
    ! ----------------------------------------------------------------

    ! Initialize analysis object state.
    subroutine new( self )
        class(class_compatibility), intent(inout) :: self
        call self%kill()
    end subroutine new

    ! Release owned resources and reset object state.
    subroutine kill( self )
        class(class_compatibility), intent(inout) :: self
        call self%kill_support_model()
    end subroutine kill

    ! Reset module-level support model state.
    subroutine kill_support_model( self )
        class(class_compatibility), intent(inout) :: self
        integer                                   :: i
        if( allocated(self%support_model%training_set) ) then
            do i = 1, size(self%support_model%training_set)
                if( allocated(self%support_model%training_set(i)%min_dim) ) deallocate(self%support_model%training_set(i)%min_dim)
                if( allocated(self%support_model%training_set(i)%max_dim) ) deallocate(self%support_model%training_set(i)%max_dim)
            end do
            deallocate(self%support_model%training_set)
        end if
        self%support_model%axis_c         = 0.0
        self%support_model%axis_b         = 0.0
        self%support_model%axis_a         = 0.0
        self%support_model%delta_c        = 0.0
        self%support_model%delta_b        = 0.0
        self%support_model%delta_a        = 0.0
        self%support_model%used_relax     = DEFAULT_RELAX
        self%support_model%used_qlo       = DEFAULT_QLO
        self%support_model%used_qhi       = DEFAULT_QHI
        self%support_model%used_soft_frac = SUPPORT_EDGE_SOFT_FRAC_MAX
        self%support_model%delta_valid    = .false.
        self%support_model%converged      = .false.
        self%support_model%valid          = .false.
    end subroutine kill_support_model

    ! ----------------------------------------------------------------
    ! getters
    ! ----------------------------------------------------------------

    ! Expose convergence status of the fitted support model.
    logical function converged( self ) result( l_converged )
        class(class_compatibility), intent(in) :: self
        l_converged = self%support_model%converged
    end function converged

    ! Return fitted axes and latest axis deltas for monitoring/logging.
    subroutine get_support_model_metrics( self, metrics )
        class(class_compatibility), intent(in)  :: self
        type(support_model_metrics), intent(out) :: metrics
        metrics%axis_c      = self%support_model%axis_c
        metrics%axis_b      = self%support_model%axis_b
        metrics%axis_a      = self%support_model%axis_a
        metrics%delta_c     = self%support_model%delta_c
        metrics%delta_b     = self%support_model%delta_b
        metrics%delta_a     = self%support_model%delta_a
        metrics%delta_valid = self%support_model%delta_valid
        metrics%valid       = self%support_model%valid
        metrics%converged   = self%support_model%converged
    end subroutine get_support_model_metrics

    ! ----------------------------------------------------------------
    ! model training
    ! ----------------------------------------------------------------

    ! Train support model from class averages in a project.
    subroutine train_1( self, spproj )
        class(class_compatibility),           intent(inout) :: self
        type(sp_project),                     intent(inout) :: spproj
        type(image),                          allocatable   :: imgs(:)
        integer,                              allocatable   :: states(:)
        imgs   = read_cavgs_into_imgarr(spproj)
        states = spproj%os_cls2D%get_all_asint('state')
        call self%new_support_training_set(imgs, states)
        call self%update_support_model()
        call dealloc_imgarr(imgs)
    end subroutine train_1

    ! Train support model from class averages in a stack.
    subroutine train_2( self, refs )
        class(class_compatibility),           intent(inout) :: self
        type(string),                         intent(in)    :: refs
        type(image),                          allocatable   :: imgs(:)
        integer,                              allocatable   :: states(:)
        imgs = read_stk_into_imgarr(refs)
        allocate(states(size(imgs)), source=1)
        call self%new_support_training_set(imgs, states)
        call self%update_support_model()
        call dealloc_imgarr(imgs)
    end subroutine train_2

    ! ----------------------------------------------------------------
    ! model inference
    ! ----------------------------------------------------------------

    ! Apply current support model to active classes and update project selection.
    subroutine infer( self, spproj )
        class(class_compatibility), intent(inout) :: self
        type(sp_project),           intent(inout) :: spproj
        type(image),                allocatable   :: imgs(:)
        integer,                    allocatable   :: states(:)
        integer                                   :: nset, iimg
        real                                      :: min_dim, max_dim
        imgs   = read_cavgs_into_imgarr(spproj)
        states = spproj%os_cls2D%get_all_asint('state')
        if( size(imgs) /= size(states) ) THROW_HARD('class_compatibility: training set size mismatch')
        do iimg = 1, size(imgs)
            if( states(iimg) <= 0 ) cycle
            call self%preprocess_img(imgs(iimg), min_dim, max_dim)
            call self%apply_support_model(states(iimg), min_dim, max_dim)
            if(states(iimg) <= 0) call spproj%os_cls2D%set(iimg, 'rejection_reason', string('class_compatibility: size_incompatible_subset'))
        end do
        call spproj%map_cavgs_selection(states)
        call dealloc_imgarr(imgs)
    end subroutine infer

    ! ----------------------------------------------------------------
    ! data preparation
    ! ----------------------------------------------------------------

    ! Build one training batch from active classes and append it to model history.
    subroutine new_support_training_set( self, imgs, states )
        class(class_compatibility), intent(inout) :: self
        type(image),                intent(in)    :: imgs(:)
        integer,                    intent(in)    :: states(:)
        type(support_model_training_set)          :: new_set
        type(support_model_training_set), allocatable :: tmp(:)
        integer                                   :: nset, iimg, idx
        if( size(imgs) /= size(states) ) THROW_HARD('class_compatibility: training set size mismatch')
        nset = count(states > 0)
        if( nset < 3 ) then
            write(logfhandle,'(A)') 'class_compatibility: support_model training set too small (n='//int2str(nset)//')'
            return
        end if
        allocate(new_set%min_dim(nset), new_set%max_dim(nset))
        idx = 0
        do iimg = 1, size(imgs)
            if( states(iimg) <= 0 ) cycle
            idx = idx + 1
            call self%preprocess_img(imgs(iimg), new_set%min_dim(idx), new_set%max_dim(idx))
        end do
        if( .not. allocated(self%support_model%training_set) ) then
            allocate(self%support_model%training_set(1))
            self%support_model%training_set(1) = new_set
        else
            nset = size(self%support_model%training_set)
            allocate(tmp(nset+1))
            tmp(1:nset) = self%support_model%training_set
            tmp(nset+1) = new_set
            call move_alloc(tmp, self%support_model%training_set)
        end if    
    end subroutine new_support_training_set

    ! Extract Feret min/max dimensions from a preprocessed class-average mask.
    subroutine preprocess_img( self, img, min_dim, max_dim )
        class(class_compatibility), intent(inout) :: self
        type(image),                intent(in)    :: img
        real,                       intent(out)   :: min_dim, max_dim
        integer,                    allocatable   :: cc_sizes(:)
        real,                       allocatable   :: weight_vals(:,:,:)
        type(image)                               :: img_loc
        type(image_bin)                           :: mask_loc, cc_img
        integer                                   :: imorph, nx, ny, i, j
        real                                      :: cx, cy, r2, r2max, radial_w
        integer                                   :: ldim(3), ldim_target(3)
        min_dim = 0.0
        max_dim = 0.0
        ! Preprocess the input image: zero edge average and bandpass filter.
        call img_loc%copy(img)
        call img_loc%zero_edgeavg()
        call img_loc%bp(0., 10.)

        ! Frequency-domain resize to PREPROCESS_BOXSIZE: FFT -> crop/pad -> inverse FFT.
        call img_loc%fft()
        ldim        = img_loc%get_ldim()
        ldim_target = [PREPROCESS_BOXSIZE, PREPROCESS_BOXSIZE, ldim(3)]
        if( ldim(1) > ldim_target(1) .or. ldim(2) > ldim_target(2) )then
            call img_loc%clip_inplace(ldim_target)
        else if( ldim(1) < ldim_target(1) .or. ldim(2) < ldim_target(2) )then
            call img_loc%pad_inplace(ldim_target)
        end if
        call img_loc%ifft()
        call img_loc%set_smpd(img_loc%get_smpd() * real(ldim(1)) / real(ldim_target(1)))

        ! Center-weight intensities before Otsu so dark/bright edge backgrounds
        ! are less likely to dominate the selected foreground component.
        weight_vals = img_loc%get_rmat()
        nx = size(weight_vals, 1)
        ny = size(weight_vals, 2)
        cx = 0.5 * real(nx + 1)
        cy = 0.5 * real(ny + 1)
        r2max = max(1.0, (0.5 * real(min(nx, ny)))**2)
        do j = 1, ny
            do i = 1, nx
            r2 = (real(i) - cx)**2 + (real(j) - cy)**2
            radial_w = 0.2 + 0.8 * max(0.0, 1.0 - min(1.0, r2 / r2max))
            weight_vals(i,j,1) = weight_vals(i,j,1) * radial_w
            end do
        end do
        call img_loc%set_rmat(weight_vals, .false.)
        if( allocated(weight_vals) ) deallocate(weight_vals)

        ! Compute Otsu threshold and binary mask.
        call otsu_img(img_loc)

        ! Apply 5 px morphological closing to the binary Otsu mask.
        call mask_loc%transfer2bimg(img_loc)
        do imorph = 1, PREPROCESS_MORPH_SIZE
            call mask_loc%dilate()
        end do
        do imorph = 1, PREPROCESS_MORPH_SIZE
            call mask_loc%erode()
        end do

        ! Keep only the largest connected foreground component in the mask.
        call mask_loc%find_ccs(cc_img)
        cc_sizes = cc_img%size_ccs()
        if( size(cc_sizes) > 0 .and. maxval(cc_sizes) > 0 )then
            call cc_img%cc2bin(maxloc(cc_sizes, dim=1))
            call mask_loc%copy_bimg(cc_img)
        end if

        ! Compute Feret diameters of the largest connected component in the mask.
        call mask_loc%feret_minmax(min_dim, max_dim)

        ! Cleanup temporary mask/CC buffers.
        if( allocated(cc_sizes) ) deallocate(cc_sizes)
        call cc_img%kill_bimg()
        call mask_loc%kill_bimg()
    end subroutine preprocess_img

    ! Apply support-model bounds and mark size-based rejections.
    ! A candidate is compatible when b falls in its relaxed [min,max] interval
    ! and [min,max] stays within relaxed global c/a bounds.
    subroutine apply_support_model( self, state, min_dim, max_dim )
        class(class_compatibility), intent(inout) :: self
        integer,                    intent(inout) :: state
        real,                       intent(in)    :: min_dim, max_dim
        logical :: compatible
        real    :: relax

        if( .not. self%support_model%valid ) return
        if( state <= 0 ) return

        relax = self%support_model%used_relax
        compatible = self%support_model%axis_b >= min_dim * (1.0 - relax) .and. &
            self%support_model%axis_b <= max_dim * (1.0 + relax) .and. &
            min_dim >= self%support_model%axis_c * (1.0 - relax) .and. &
            max_dim <= self%support_model%axis_a * (1.0 + relax)

        if( .not. compatible )then
            if( min_dim >= self%support_model%axis_c * (1.0 - relax - RESCUE_EDGE_FRAC) .and. &
                max_dim <= self%support_model%axis_a * (1.0 + relax + RESCUE_EDGE_FRAC) )then
                compatible = .true.
            end if
        end if

        if( .not. compatible )then
            state = 0
        end if
    end subroutine apply_support_model

    ! Fit latent support-model axes (c <= b <= a) from cached training batches.
    subroutine update_support_model( self )
        class(class_compatibility), intent(inout) :: self
        integer,                    allocatable   :: batch_id(:), batch_counts(:)
        real,                       allocatable   :: mins(:), maxs(:), cands(:), mins_sup(:), maxs_sup(:), w_sup(:), sample_weights(:)
        logical,                    allocatable   :: tmp_support(:), support_cur(:)
        integer :: nn, i, j, k, best_count
        integer :: ir, iq, jq, final_count, best_final_count
        integer :: nset, ntot, idx, nbatches
        logical :: in_b, in_c, in_a, had_prev_fit
        logical :: stable_a, stable_b, stable_c
        real    :: b_cand, left, right
        real    :: spread, best_spread, score, best_score
        real    :: weighted_count, best_weighted_count, best_final_weighted_count
        real    :: relax_cur, qlo_cur, qhi_cur
        real    :: c_cur, b_cur, a_cur, c_soft, a_soft
        real    :: c_best, b_best, a_best
        real    :: relax_best, qlo_best, qhi_best
        real    :: soft_frac, dc_norm, db_norm, da_norm, stability
        real    :: c_prev, b_prev, a_prev
        real    :: dc, db, da

        had_prev_fit = self%support_model%valid
        c_prev = self%support_model%axis_c
        b_prev = self%support_model%axis_b
        a_prev = self%support_model%axis_a

        ! Start permissive (5%) and decay with axis-update stability from prior fit.
        soft_frac = SUPPORT_EDGE_SOFT_FRAC_MAX
        if( self%support_model%delta_valid )then
            dc_norm = self%support_model%delta_c / max(CONV_ABS_EPS, CONV_REL_EPS * max(abs(c_prev), 1.0))
            db_norm = self%support_model%delta_b / max(CONV_ABS_EPS, CONV_REL_EPS * max(abs(b_prev), 1.0))
            da_norm = self%support_model%delta_a / max(CONV_ABS_EPS, CONV_REL_EPS * max(abs(a_prev), 1.0))
            stability = max(dc_norm, db_norm, da_norm)
            soft_frac = SUPPORT_EDGE_SOFT_FRAC_MAX * min(1.0, max(0.0, stability))
        end if
        if( self%support_model%converged ) soft_frac = 0.0

        if( .not. allocated(self%support_model%training_set) )then
            write(logfhandle,'(A)') 'cavg_compat: support_model training set is not allocated'
            return
        end if

        ! Flatten all valid batches into one working view while retaining batch ids
        ! so each batch can contribute with equal total weight.
        ntot = 0
        nbatches = 0
        do nset = 1, size(self%support_model%training_set)
            if( allocated(self%support_model%training_set(nset)%min_dim) .neqv. allocated(self%support_model%training_set(nset)%max_dim) )then
                write(logfhandle,'(A)') 'cavg_compat: support_model invalid training allocation state (min/max mismatch)'
                return
            end if
            if( .not. allocated(self%support_model%training_set(nset)%min_dim) ) cycle
            if( size(self%support_model%training_set(nset)%min_dim) /= size(self%support_model%training_set(nset)%max_dim) )then
                write(logfhandle,'(A)') 'cavg_compat: support_model training set size mismatch'
                return
            end if
            if( size(self%support_model%training_set(nset)%min_dim) == 0 ) cycle
            nbatches = nbatches + 1
            ntot = ntot + size(self%support_model%training_set(nset)%min_dim)
        end do

        if( ntot <= 2 )then
            write(logfhandle,'(A)') 'cavg_compat: support_model skipped (training set too small)'
            return
        end if

        if( nbatches < 1 )then
            write(logfhandle,'(A)') 'cavg_compat: support_model invalid batch metadata'
            return
        end if

        allocate(mins(ntot), maxs(ntot), batch_id(ntot))
        idx = 0
        nn = 0
        do nset = 1, size(self%support_model%training_set)
            if( .not. allocated(self%support_model%training_set(nset)%min_dim) ) cycle
            if( size(self%support_model%training_set(nset)%min_dim) == 0 ) cycle
            nn = nn + 1
            k = size(self%support_model%training_set(nset)%min_dim)
            mins(idx+1:idx+k) = self%support_model%training_set(nset)%min_dim
            maxs(idx+1:idx+k) = self%support_model%training_set(nset)%max_dim
            batch_id(idx+1:idx+k) = nn
            idx = idx + k
        end do
        nn = ntot

        allocate(batch_counts(nbatches), source=0)
        allocate(sample_weights(nn),     source=0.0)
        do i = 1, nn
            batch_counts(batch_id(i)) = batch_counts(batch_id(i)) + 1
        end do
        do i = 1, nn
            sample_weights(i) = 1.0 / real(batch_counts(batch_id(i)))
        end do

        allocate(tmp_support(nn), source=.false.)
        allocate(support_cur(nn), source=.false.)
        allocate(cands(2*nn))
        k = 0
        do i = 1, nn
            k = k + 1
            cands(k) = mins(i)
            k = k + 1
            cands(k) = maxs(i)
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
        relax_best = DEFAULT_RELAX
        qlo_best   = DEFAULT_QLO
        qhi_best   = DEFAULT_QHI

        ! Grid search over relaxation and quantile hyperparameters.
        do ir = 1, NRELAX
            relax_cur = RELAX_GRID(ir)

            do iq = 1, NQ
                qlo_cur = QLOW_GRID(iq)

                do jq = 1, NQ
                    qhi_cur = QHIGH_GRID(jq)
                    if( qhi_cur <= qlo_cur ) cycle

                    best_count = -1
                    best_weighted_count = -1.0
                    best_spread = huge(1.0)
                    b_cur = 0.0
                    support_cur = .false.

                    ! First stage: pick b that maximizes weighted in-interval support.
                    do i = 1, k
                        b_cand = cands(i)
                        final_count = 0
                        weighted_count = 0.0
                        spread = 0.0
                        tmp_support = .false.

                        do j = 1, nn
                            left = mins(j) * (1.0 - relax_cur)
                            right = maxs(j) * (1.0 + relax_cur)
                            if( b_cand >= left .and. b_cand <= right )then
                                final_count = final_count + 1
                                weighted_count = weighted_count + sample_weights(j)
                                tmp_support(j) = .true.
                                spread = spread + sample_weights(j) * abs(b_cand - 0.5 * (mins(j) + maxs(j))) / &
                                    max(maxs(j)-mins(j), 1.0e-6)
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
                        mins_sup(j) = mins(i)
                        maxs_sup(j) = maxs(i)
                        w_sup(j) = sample_weights(i)
                    end do

                    ! Second stage: estimate c/a from weighted support quantiles.
                    c_cur = percentile_weighted_real(mins_sup, w_sup, qlo_cur)
                    a_cur = percentile_weighted_real(maxs_sup, w_sup, qhi_cur)
                    if( c_cur > b_cur ) c_cur = minval(mins_sup)
                    if( a_cur < b_cur ) a_cur = maxval(maxs_sup)

                    final_count = 0
                    weighted_count = 0.0
                    spread = 0.0
                    tmp_support = .false.
                    do i = 1, nn
                        left = mins(i) * (1.0 - relax_cur)
                        right = maxs(i) * (1.0 + relax_cur)
                        c_soft = c_cur * (1.0 - relax_cur - soft_frac)
                        a_soft = a_cur * (1.0 + relax_cur + soft_frac)
                        in_b = (b_cur >= left .and. b_cur <= right)
                        in_c = (mins(i) >= c_cur * (1.0 - relax_cur))
                        in_a = (maxs(i) <= a_cur * (1.0 + relax_cur))
                        if( in_b .and. ( (in_c .and. in_a) .or. (mins(i) >= c_soft .and. in_a) .or. (in_c .and. maxs(i) <= a_soft) ) )then
                            final_count = final_count + 1
                            weighted_count = weighted_count + sample_weights(i)
                            tmp_support(i) = .true.
                            spread = spread + sample_weights(i) * abs(b_cur - 0.5 * (mins(i) + maxs(i))) / &
                                max(maxs(i)-mins(i), 1.0e-6)
                        end if
                    end do

                    ! Prefer high weighted support, then compactness, then less relaxation.
                    score = weighted_count - 0.01 * spread - 0.10 * relax_cur
                    if( score > best_score .or. (abs(score - best_score) <= 1.0e-6 .and. weighted_count > best_final_weighted_count) )then
                        best_score = score
                        best_final_count = final_count
                        best_final_weighted_count = weighted_count
                        c_best = c_cur
                        b_best = b_cur
                        a_best = a_cur
                        relax_best = relax_cur
                        qlo_best   = qlo_cur
                        qhi_best   = qhi_cur
                    end if

                    deallocate(mins_sup, maxs_sup, w_sup)
                end do
            end do
        end do

        if( best_final_count <= 0 )then
            deallocate(cands, tmp_support, support_cur, mins, maxs, batch_id, batch_counts, sample_weights)
            return
        end if

        self%support_model%axis_c = c_best
        self%support_model%axis_b = b_best
        self%support_model%axis_a = a_best
        self%support_model%used_relax = relax_best
        self%support_model%used_qlo   = qlo_best
        self%support_model%used_qhi   = qhi_best
        self%support_model%used_soft_frac = soft_frac
        self%support_model%valid = .true.

        ! Convergence requires a previous fit and stable c/b/a updates.
        dc = abs(self%support_model%axis_c - c_prev)
        db = abs(self%support_model%axis_b - b_prev)
        da = abs(self%support_model%axis_a - a_prev)
        self%support_model%delta_c = dc
        self%support_model%delta_b = db
        self%support_model%delta_a = da
        self%support_model%delta_valid = had_prev_fit

        stable_c = dc <= CONV_ABS_EPS .or. dc <= CONV_REL_EPS * max(abs(c_prev), abs(self%support_model%axis_c), 1.0)
        stable_b = db <= CONV_ABS_EPS .or. db <= CONV_REL_EPS * max(abs(b_prev), abs(self%support_model%axis_b), 1.0)
        stable_a = da <= CONV_ABS_EPS .or. da <= CONV_REL_EPS * max(abs(a_prev), abs(self%support_model%axis_a), 1.0)

        self%support_model%converged = had_prev_fit .and. stable_c .and. stable_b .and. stable_a

        deallocate(cands, tmp_support, support_cur, mins, maxs, batch_id, batch_counts, sample_weights)

    contains

        ! Return the nearest-rank percentile of a real-valued vector.
        real function percentile_real(arr, q) result(pval)
            real, intent(in) :: arr(:)
            real, intent(in) :: q
            real, allocatable :: sorted(:)
            integer :: nvals, idx_local

            nvals = size(arr)
            if( nvals <= 0 )then
                pval = 0.0
                return
            end if

            allocate(sorted(nvals))
            sorted = arr
            call hpsort(sorted)

            idx_local = 1 + int(q * real(max(nvals - 1, 0)))
            idx_local = max(1, min(nvals, idx_local))
            pval = sorted(idx_local)
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
        
    end subroutine update_support_model

end module simple_class_compatibility
