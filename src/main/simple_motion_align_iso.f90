! linearly align frames in a maximum-likelihood way, double (nested) iteration.
! convergence and frameweights can be determined by callback
module simple_motion_align_iso
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_error
use simple_ft_expanded,           only: ft_expanded
use simple_ftexp_shsrch,          only: ftexp_shsrch
use simple_image,                 only: image
use simple_parameters,            only: params_glob
use simple_opt_lbfgsb,            only: PRINT_NEVALS
use CPlot2D_wrapper_module
implicit none
public :: motion_align_iso, align_iso_callback, align_iso_fw_callback
private
#include "simple_local_flags.inc"

type :: motion_align_iso
    private
    logical, public                                   :: existence = .false.
    type(image), pointer                              :: frames(:)                         !< pointer to stack of frames
    integer                                           :: nframes_allocd = 0                !< number of frames allocated
    integer                                           :: nframes        = 0                !< number of frames
    integer                                           :: iter           = 0                !< iteration number
    integer                                           :: mitsref                           !< maximum iteration number
    integer                                           :: fixed_frame = 1                   !< fixed (non-shifted) frame
    real                                              :: trs                               !< size of box constraint
    real                                              :: ftol                              !< tolerance parameter for minimizer
    real                                              :: gtol                              !< tolerance parameter for minimizer
    integer                                           :: nimproved                         !< number of improved frames
    real                                              :: frac_improved                     !< fraction of improved frames
    real                                              :: smallshift                        !< range for initial random shift
    logical                                           :: rand_init_shifts                  !< randomize initial condition?
    real                                              :: corrfrac                          !< quotient new/old correlation
    real                                              :: hp = -1.                          !< high pass value
    real                                              :: lp = -1.                          !< low pass value
    real                                              :: corr                              !< correlation
    real                                              :: corr_prev                         !< previous correlation
    real                                              :: corr_saved                        !< best correlation
    real                                              :: shsrch_tol                        !< tolerance parameter for shsrch update
    integer                                           :: maxits = 0                        !< maximum number of iterations; 0: default
    logical                                           :: has_shsrch_tol = .false.          !< is shsrch_tol set
    integer                                           :: coord_x, coord_y                  !< x,y coordinates for patch-based alignment callback
    real,              allocatable                    :: shifts_toplot(:,:)                !< for plotting
    real,              allocatable                    :: opt_shifts(:,:)                   !< shifts identified
    real,              allocatable                    :: opt_shifts_saved(:,:)             !< best shifts
    real,              allocatable                    :: frameweights(:)                   !< array of frameweights
    real,              allocatable                    :: frameweights_saved(:)             !< array of frameweights
    real,              allocatable                    :: corrs(:)                          !< per-frame correlations
    type(ft_expanded), allocatable                    :: movie_frames_ftexp(:)             !< movie frames
    type(ft_expanded), allocatable                    :: movie_frames_ftexp_sh(:)          !< shifted movie frames
    type(ft_expanded), allocatable                    :: movie_sum_global_ftexp_threads(:) !< array of global movie sums for parallel refinement
    procedure(align_iso_callback),    pointer, nopass :: callback              => null()   !< callback to determine convergence
    procedure(align_iso_fw_callback), pointer, nopass :: frameweights_callback => null()   !< callback to compute frameweights
                                                                                           !< if not set, use default weighting
contains
    procedure, private                                :: motion_align_iso_new
    procedure, private                                :: motion_align_iso_align
    procedure, private                                :: motion_align_iso_kill
    procedure, private                                :: recenter_shifts                   !< put shifts rel. to fixed frame
    procedure, private                                :: calc_frameweights                 !< compute new frameweights
    procedure, private                                :: corrmat2weights
    procedure, private                                :: allocate_fields
    procedure, private                                :: deallocate_fields
    procedure, private                                :: create_ftexp_objs
    procedure                                         :: shift_wsum_and_calc_corrs         !< shift, sum and calculate new correlatons
    procedure                                         :: new                  => motion_align_iso_new
    procedure                                         :: align                => motion_align_iso_align
    procedure                                         :: kill                 => motion_align_iso_kill
    procedure                                         :: set_shsrch_tol       => motion_align_iso_set_shsrch_tol
    procedure                                         :: set_frames           => motion_align_iso_set_frames
    procedure                                         :: set_hp_lp            => motion_align_iso_set_hp_lp
    procedure                                         :: set_trs              => motion_align_iso_set_trs
    procedure                                         :: set_ftol_gtol        => motion_align_iso_set_ftol_gtol
    procedure                                         :: set_smallshift       => motion_align_iso_set_smallshift
    procedure                                         :: get_iter             => motion_align_iso_get_iter
    procedure                                         :: set_weights          => motion_align_iso_set_weights
    procedure                                         :: set_even_weights     => motion_align_iso_set_even_weights
    procedure                                         :: get_weights          => motion_align_iso_get_weights
    procedure                                         :: set_rand_init_shifts => motion_align_iso_set_rand_init_shifts
    procedure                                         :: get_corr             => motion_align_iso_get_corr
    procedure                                         :: get_corrs            => motion_align_iso_get_corrs
    procedure                                         :: get_opt_shifts       => motion_align_iso_get_opt_shifts
    procedure                                         :: get_shifts_toplot    => motion_align_iso_get_shifts_toplot
    procedure                                         :: set_mitsref          => motion_align_iso_set_mitsref
    procedure                                         :: set_fixed_frame      => motion_align_iso_set_fixed_frame
    procedure                                         :: get_frac_improved    => motion_align_iso_get_frac_improved
    procedure                                         :: get_nimproved        => motion_align_iso_get_nimproved
    procedure                                         :: get_corrfrac         => motion_align_iso_get_corrfrac
    procedure                                         :: set_maxits           => motion_align_iso_set_maxits
    procedure                                         :: set_coords           => motion_align_iso_set_coords
    procedure                                         :: get_coords           => motion_align_iso_get_coords
    procedure                                         :: set_callback         => motion_align_iso_set_callback
    procedure                                         :: set_frameweights_callback => motion_align_iso_set_frameweights_callback
end type motion_align_iso

abstract interface
    subroutine align_iso_callback(aptr, align_iso, converged)
        import motion_align_iso
        class(*),                intent(inout) :: aptr
        class(motion_align_iso), intent(inout) :: align_iso
        logical,                 intent(out)   :: converged
    end subroutine align_iso_callback

    subroutine align_iso_fw_callback(aptr, align_iso)
        import motion_align_iso
        class(*),                intent(inout) :: aptr
        class(motion_align_iso), intent(inout) :: align_iso
    end subroutine align_iso_fw_callback
end interface

real,    parameter :: SMALLSHIFT_DEFAULT = 2.
integer, parameter :: MITSREF_DEFAULT    = 30
real,    parameter :: NIMPROVED_TOL      = 1e-7

contains

    subroutine motion_align_iso_new( self )
        class(motion_align_iso), intent(inout) :: self
        call self%kill()
        self%trs                   = params_glob%trs
        self%ftol                  = params_glob%motion_correctftol
        self%gtol                  = params_glob%motion_correctgtol
        self%fixed_frame           = 1
        self%smallshift            = SMALLSHIFT_DEFAULT
        self%nframes_allocd        = 0
        self%rand_init_shifts      = .false.
        self%mitsref               = MITSREF_DEFAULT
        self%nimproved             = -1
        self%frac_improved         = -1.
        self%callback              => null()
        self%frameweights_callback => null()
        self%hp                    = -1.
        self%lp                    = -1.
        self%existence             = .true.
    end subroutine motion_align_iso_new

    subroutine motion_align_iso_align( self, callback_ptr )
        class(motion_align_iso),    intent(inout) :: self
        class(*),                   intent(inout) :: callback_ptr  !< callback pointer to be passed as first argument
        type(ftexp_shsrch), allocatable           :: ftexp_srch(:)
        integer :: i, iter, iframe, nimproved, maxits_saved
        real    :: cxy(3)
        logical :: callback_convgd
        real    :: hp_saved, lp_saved
        if ( .not. self%existence ) then
            THROW_HARD('not instantiated; simple_motion_align_iso: align')
        end if
        if ( self%nframes < 2 ) then
            THROW_HARD('nframes < 2; simple_motion_align_iso: align')
        end if
        if (( self%hp < 0. ) .or. ( self%lp < 0.)) then
            THROW_HARD('hp or lp < 0; simple_motion_align_iso: align')
        end if
        self%iter = 0
        call self%create_ftexp_objs
        call self%calc_frameweights( callback_ptr )
        if ( self%rand_init_shifts ) then
            do iframe = 1, self%nframes
                self%opt_shifts(iframe,1) = ran3() * 2. * self%smallshift - self%smallshift
                self%opt_shifts(iframe,2) = ran3() * 2. * self%smallshift - self%smallshift
            end do
        else
            self%opt_shifts = 0.
        end if
        call self%recenter_shifts(self%opt_shifts)
        ! generate movie sum for refinement
        call self%shift_wsum_and_calc_corrs
        self%corr = sum(self%corrs)/real(self%nframes)
        allocate( ftexp_srch(self%nframes), stat=alloc_stat )
        if (alloc_stat.ne.0) call allocchk('align 1; simple_motion_align_init')
        do iframe = 1, self%nframes
            call ftexp_srch(iframe)%new(self%movie_sum_global_ftexp_threads(iframe),&
                self%movie_frames_ftexp(iframe),&
                self%trs, motion_correct_ftol=self%ftol, motion_correct_gtol=self%gtol)
            if (self%has_shsrch_tol) then
                call ftexp_srch(iframe)%set_shsrch_tol(self%shsrch_tol)
            end if
            if (self%maxits > 0) ftexp_srch(iframe)%ospec%maxits = self%maxits
        end do
        self%corr_saved = -1.
        do iter=1,self%mitsref
            self%iter = iter
            nimproved = 0
            !$omp parallel do default(shared) private(iframe,cxy) proc_bind(close) reduction(+:nimproved)
            do iframe=1,self%nframes
                call self%movie_sum_global_ftexp_threads(iframe)%subtr(&
                    self%movie_frames_ftexp_sh(iframe), w=self%frameweights(iframe))
                PRINT_NEVALS = .false.
                cxy = ftexp_srch(iframe)%minimize(self%corrs(iframe), self%opt_shifts(iframe,:))
                if( cxy(1) - self%corrs(iframe) > NIMPROVED_TOL ) nimproved = nimproved + 1
                self%opt_shifts(iframe,:) = cxy(2:3)
                self%corrs(iframe)        = cxy(1)
            end do
            !$omp end parallel do
            self%nimproved = nimproved
            self%frac_improved = real(self%nimproved) / real(self%nframes) * 100.
            write(logfhandle,'(a,1x,f4.0)') 'This % of frames improved their alignment: ', self%frac_improved
            call self%recenter_shifts(self%opt_shifts)
            call self%calc_frameweights( callback_ptr )
            ! build new reference
            call self%shift_wsum_and_calc_corrs
            self%corr_prev = self%corr
            self%corr      = sum(self%corrs) / real(self%nframes)
            if( self%corr >= self%corr_saved ) then ! save the local optimum
                self%corr_saved         = self%corr
                self%opt_shifts_saved   = self%opt_shifts
                self%frameweights_saved = self%frameweights
            endif
            self%corrfrac = self%corr_prev / self%corr
            if (associated(self%callback)) then
                hp_saved = self%hp
                lp_saved = self%lp
                callback_convgd = .false.
                maxits_saved = self%maxits
                call self%callback(callback_ptr, self, callback_convgd)
                if (maxits_saved /= self%maxits) then
                    do iframe = 1, self%nframes
                        ftexp_srch(iframe)%ospec%maxits = self%maxits
                    end do
                end if
                if (callback_convgd) exit
                if ((abs(hp_saved-self%hp) > epsilon(hp_saved)) .or. &
                    (abs(lp_saved-self%lp) > epsilon(lp_saved))) then
                    ! need to re-make the ftexps
                    call self%create_ftexp_objs
                    call self%recenter_shifts(self%opt_shifts)
                    call self%calc_frameweights( callback_ptr )
                    call self%shift_wsum_and_calc_corrs
                    ! need to destroy all previous knowledge about correlations
                    self%corr       = sum(self%corrs) / real(self%nframes)
                    self%corr_prev  = self%corr
                    self%corr_saved = self%corr
                end if
            end if
        end do
        ! put the best local optimum back
        self%corr         = self%corr_saved
        self%opt_shifts   = self%opt_shifts_saved
        self%frameweights = self%frameweights_saved
        do i = 1, self%nframes
            call ftexp_srch(i)%kill
        end do
        deallocate(ftexp_srch)
    end subroutine motion_align_iso_align

    subroutine motion_align_iso_kill( self )
        class(motion_align_iso), intent(inout) :: self
        call self%deallocate_fields
        self%nframes_allocd        =    0
        self%callback              =>   null()
        self%frameweights_callback =>   null()
        self%existence             =    .false.
    end subroutine motion_align_iso_kill

    subroutine recenter_shifts( self, shifts )
        class(motion_align_iso), intent(inout) :: self
        real,     intent(inout) :: shifts(self%nframes,2)
        integer :: iframe
        real :: xsh, ysh
        xsh = -shifts(self%fixed_frame,1)
        ysh = -shifts(self%fixed_frame,2)
        do iframe=1,self%nframes
            shifts(iframe,1) = shifts(iframe,1) + xsh
            shifts(iframe,2) = shifts(iframe,2) + ysh
            if( abs(shifts(iframe,1)) < 1e-6 ) shifts(iframe,1) = 0.
            if( abs(shifts(iframe,2)) < 1e-6 ) shifts(iframe,2) = 0.
            self%shifts_toplot(iframe,1) = shifts(iframe, 1)
            self%shifts_toplot(iframe,2) = shifts(iframe, 2)
        end do
    end subroutine recenter_shifts

    subroutine calc_frameweights( self, callback_ptr )
        class(motion_align_iso), intent(inout) :: self
        class(*),                intent(inout) :: callback_ptr
        if ( associated(self%frameweights_callback) ) then
            call self%frameweights_callback(callback_ptr, self)
        else
            if( self%iter == 0 ) then
                call self%corrmat2weights ! initialisation
            else
                self%frameweights = corrs2weights(self%corrs) ! update
            endif
        end if
    end subroutine calc_frameweights

    subroutine corrmat2weights( self )
        class(motion_align_iso), intent(inout) :: self
        integer :: iframe, jframe
        real    :: corrmat(self%nframes,self%nframes)
        corrmat = 1. ! diagonal elements are 1
        self%corrs   = 0.
        !$omp parallel default(shared) private(iframe,jframe) proc_bind(close)
        !$omp do schedule(guided)
        do iframe=1,self%nframes-1
            do jframe=iframe+1,self%nframes
                corrmat(iframe,jframe) = self%movie_frames_ftexp_sh(iframe)%corr(self%movie_frames_ftexp_sh(jframe))
                corrmat(jframe,iframe) = corrmat(iframe,jframe)
            end do
        end do
        !$omp end do
        !$omp do schedule(static)
        do iframe=1,self%nframes
            do jframe=1,self%nframes
                if( jframe == iframe ) cycle
                self%corrs(iframe) = self%corrs(iframe)+corrmat(iframe,jframe)
            end do
            self%corrs(iframe) = self%corrs(iframe)/real(self%nframes-1)
        end do
        !$omp end do
        !$omp end parallel
        self%frameweights = corrs2weights(self%corrs)
    end subroutine corrmat2weights

    subroutine allocate_fields( self )
        class(motion_align_iso), intent(inout) :: self
        allocate(&
            self%movie_frames_ftexp(self%nframes),&
            self%movie_frames_ftexp_sh(self%nframes),&
            self%movie_sum_global_ftexp_threads(self%nframes),&
            self%shifts_toplot(self%nframes,2),&
            self%opt_shifts(self%nframes,2),&
            self%opt_shifts_saved(self%nframes,2),&
            self%corrs(self%nframes),&
            self%frameweights(self%nframes),&
            self%frameweights_saved(self%nframes), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('allocate_fields; simple_motion_align_iso')
        self%nframes_allocd = self%nframes
    end subroutine allocate_fields

    subroutine deallocate_fields( self )
        class(motion_align_iso), intent(inout) :: self
        integer :: i
        if (self%nframes_allocd > 0) then
            do i = 1, self%nframes_allocd
                call self%movie_frames_ftexp(i)%kill
                call self%movie_frames_ftexp_sh(i)%kill
                call self%movie_sum_global_ftexp_threads(i)%kill
            end do
            deallocate( self%movie_frames_ftexp, self%movie_frames_ftexp_sh, &
                self%movie_sum_global_ftexp_threads, self%opt_shifts, self%frameweights, self%corrs, &
                self%shifts_toplot )
        else
            if ( allocated(self%movie_frames_ftexp) .or. allocated(self%movie_frames_ftexp_sh) .or.&
                 allocated(self%movie_sum_global_ftexp_threads) .or. allocated(self%opt_shifts) .or.&
                 allocated(self%frameweights) .or. allocated(self%corrs) .or. allocated(self%shifts_toplot) ) then
                THROW_HARD('inconsistency; simple_motion_align_iso: deallocate_fields')
            end if
        end if
    end subroutine deallocate_fields

    subroutine create_ftexp_objs( self )
        class(motion_align_iso), intent(inout) :: self
        integer :: iframe
        do iframe=1,self%nframes
            call self%movie_frames_ftexp(iframe)%new(self%frames(iframe), self%hp, self%lp, .true.)
            call self%movie_frames_ftexp_sh(iframe)%new(self%frames(iframe), self%hp, self%lp, .false.)
            call self%movie_sum_global_ftexp_threads(iframe)%new(self%frames(iframe), self%hp, self%lp, .false.)
        end do
    end subroutine create_ftexp_objs

    subroutine shift_wsum_and_calc_corrs( self )
        class(motion_align_iso), intent(inout) :: self
        complex, allocatable :: cmat_sum(:,:,:), cmat(:,:,:)
        integer :: iframe, flims(3,2)
        real    :: shvec(3)
        ! allocate matrices for reduction
        flims = self%movie_sum_global_ftexp_threads(1)%get_flims()
        allocate(cmat(flims(1,1):flims(1,2),flims(2,1):flims(2,2),flims(3,1):flims(3,2)),&
            cmat_sum(flims(1,1):flims(1,2),flims(2,1):flims(2,2),flims(3,1):flims(3,2)), source=cmplx(0.,0.))
        ! FIRST LOOP TO OBTAIN WEIGHTED SUM
        !$omp parallel default(shared) private(iframe,shvec,cmat) proc_bind(close)
        !$omp do schedule(static) reduction(+:cmat_sum)
        do iframe=1,self%nframes
            shvec(1) = -self%opt_shifts(iframe,1)
            shvec(2) = -self%opt_shifts(iframe,2)
            shvec(3) = 0.0
            call self%movie_frames_ftexp(iframe)%shift(shvec, self%movie_frames_ftexp_sh(iframe))
            call self%movie_frames_ftexp_sh(iframe)%get_cmat(cmat)
            cmat_sum = cmat_sum + cmat * self%frameweights(iframe)
        end do
        !$omp end do
        ! SECOND LOOP TO UPDATE movie_sum_global_ftexp_threads AND CALCULATE CORRS
        !$omp do schedule(static)
        do iframe=1,self%nframes
            ! update array of sums (for future parallel exec)
            call self%movie_sum_global_ftexp_threads(iframe)%set_cmat(cmat_sum)
            ! subtract the movie frame being correlated to reduce bias
            call self%movie_sum_global_ftexp_threads(iframe)%subtr(self%movie_frames_ftexp_sh(iframe), &
                w=self%frameweights(iframe))
            ! calc corr
            self%corrs(iframe) = self%movie_sum_global_ftexp_threads(iframe)%&
                corr(self%movie_frames_ftexp_sh(iframe))
            ! add the subtracted movie frame back to the weighted sum
            call self%movie_sum_global_ftexp_threads(iframe)%add(self%movie_frames_ftexp_sh(iframe), &
                w=self%frameweights(iframe))
        end do
        !$omp end do
        !$omp end parallel
    end subroutine shift_wsum_and_calc_corrs

    subroutine motion_align_iso_set_frames( self, frames_ptr, nframes )
        class(motion_align_iso),          intent(inout) :: self
        type(image), allocatable, target, intent(in)    :: frames_ptr(:)
        integer,                          intent(in)    :: nframes
        logical :: do_alloc
        if (.not. self%existence) then
            THROW_HARD('not instantiated; simple_motion_align_iso: set_frames')
        end if
        self%frames  => frames_ptr
        self%nframes = nframes
        if (size(frames_ptr, 1) < nframes) then
            THROW_HARD('nframes > #frames provided; simple_motion_align_iso: set_frames')
        end if
        do_alloc = .true.
        if ( self%nframes_allocd == self%nframes) then
            do_alloc = .false.
        else
            if ( self%nframes_allocd > 0 ) then
                call self%deallocate_fields
            end if
        end if
        if ( do_alloc ) call self%allocate_fields
    end subroutine motion_align_iso_set_frames

    subroutine motion_align_iso_set_hp_lp( self, hp, lp )
        class(motion_align_iso), intent(inout) :: self
        real,                    intent(in)    :: hp, lp
        if (.not. self%existence) then
            THROW_HARD('not instantiated; simple_motion_align_iso: set_hp_lp')
        end if
        self%hp = hp
        self%lp = lp
        if (.not. self%nframes > 0) then
            THROW_HARD('nframes < 1; simple_motion_align_iso: set_hp_lp')
        end if
    end subroutine motion_align_iso_set_hp_lp

    subroutine motion_align_iso_set_trs( self, trs )
        class(motion_align_iso), intent(inout) :: self
        real,                    intent(in)    :: trs
        self%trs = trs
    end subroutine motion_align_iso_set_trs

    subroutine motion_align_iso_set_ftol_gtol( self, ftol, gtol )
        class(motion_align_iso), intent(inout) :: self
        real,                    intent(in)    :: ftol, gtol
        self%ftol = ftol
        self%gtol = gtol
    end subroutine motion_align_iso_set_ftol_gtol

    subroutine motion_align_iso_set_smallshift( self, smallshift )
        class(motion_align_iso), intent(inout) :: self
        real,                    intent(in)    :: smallshift
        self%smallshift = smallshift
    end subroutine motion_align_iso_set_smallshift

    function motion_align_iso_get_iter( self ) result( iter )
        class(motion_align_iso), intent(inout) :: self
        integer :: iter
        iter = self%iter
    end function motion_align_iso_get_iter

    subroutine motion_align_iso_set_weights( self, frameweights )
        class(motion_align_iso), intent(inout) :: self
        real, allocatable,       intent(in)    :: frameweights(:)
        if (size(frameweights) /= self%nframes) then
            THROW_HARD('inconsistency; simple_motion_align_iso: set_weights')
        end if
        self%frameweights(:) = frameweights(:)
    end subroutine motion_align_iso_set_weights

    subroutine motion_align_iso_set_even_weights( self )
        class(motion_align_iso), intent(inout) :: self
        self%frameweights = 1. / real(self%nframes)
    end subroutine motion_align_iso_set_even_weights

    subroutine motion_align_iso_get_weights( self, frameweights )
        class(motion_align_iso), intent(inout) :: self
        real, allocatable,       intent(out)   :: frameweights(:)
        allocate(frameweights(self%nframes), source=self%frameweights)
    end subroutine motion_align_iso_get_weights

    subroutine motion_align_iso_set_rand_init_shifts( self, rand_init_shifts )
        class(motion_align_iso), intent(inout) :: self
        logical,                 intent(in)    :: rand_init_shifts
        self%rand_init_shifts = rand_init_shifts
    end subroutine motion_align_iso_set_rand_init_shifts

    function motion_align_iso_get_corr( self ) result( corr )
        class(motion_align_iso), intent(inout) :: self
        real :: corr
        corr = self%corr
    end function motion_align_iso_get_corr

    subroutine motion_align_iso_get_corrs( self, corrs )
        class(motion_align_iso), intent(inout) :: self
        real, allocatable,       intent(out)   :: corrs(:)
        allocate( corrs(self%nframes), source=self%corrs )
    end subroutine motion_align_iso_get_corrs

    subroutine motion_align_iso_get_opt_shifts( self, opt_shifts )
        class(motion_align_iso), intent(inout) :: self
        real, allocatable,       intent(out)   :: opt_shifts(:,:)
        allocate( opt_shifts(self%nframes, 2), source=self%opt_shifts )
    end subroutine motion_align_iso_get_opt_shifts

    subroutine motion_align_iso_get_shifts_toplot( self, shifts_toplot )
        class(motion_align_iso), intent(inout) :: self
        real, allocatable,       intent(out)   :: shifts_toplot(:,:)
        allocate( shifts_toplot(self%nframes, 2), source=self%shifts_toplot )
    end subroutine motion_align_iso_get_shifts_toplot

    subroutine motion_align_iso_set_mitsref( self, mitsref )
        class(motion_align_iso), intent(inout) :: self
        integer, intent(in) :: mitsref
        self%mitsref = mitsref
    end subroutine motion_align_iso_set_mitsref

    subroutine motion_align_iso_set_fixed_frame( self, fixed_frame )
        class(motion_align_iso), intent(inout) :: self
        integer, intent(in) :: fixed_frame
        self%fixed_frame = fixed_frame
    end subroutine motion_align_iso_set_fixed_frame

    function motion_align_iso_get_frac_improved( self ) result( frac_improved )
        class(motion_align_iso), intent(inout) :: self
        real :: frac_improved
        frac_improved = self%frac_improved
    end function motion_align_iso_get_frac_improved

    function motion_align_iso_get_nimproved( self ) result( nimproved )
        class(motion_align_iso), intent(inout) :: self
        integer :: nimproved
        nimproved = self%nimproved
    end function motion_align_iso_get_nimproved

    function motion_align_iso_get_corrfrac( self ) result( corrfrac )
        class(motion_align_iso), intent(inout) :: self
        real :: corrfrac
        corrfrac = self%corrfrac
    end function motion_align_iso_get_corrfrac

    subroutine motion_align_iso_set_shsrch_tol( self, shsrch_tol )
        class(motion_align_iso), intent(inout) :: self
        real, intent(in) :: shsrch_tol
        self%shsrch_tol = shsrch_tol
        self%has_shsrch_tol = .true.
    end subroutine motion_align_iso_set_shsrch_tol

    subroutine motion_align_iso_set_maxits( self, maxits )
        class(motion_align_iso), intent(inout) :: self
        integer, intent(in) :: maxits
        self%maxits = maxits
    end subroutine motion_align_iso_set_maxits

    subroutine motion_align_iso_set_coords( self, x, y )
        class(motion_align_iso), intent(inout) :: self
        integer, intent(in) :: x, y
        self%coord_x = x
        self%coord_y = y
    end subroutine motion_align_iso_set_coords

    subroutine motion_align_iso_get_coords(self, x, y )
        class(motion_align_iso), intent(inout) :: self
        integer, intent(out) :: x, y
        x = self%coord_x
        y = self%coord_y
    end subroutine motion_align_iso_get_coords

    subroutine motion_align_iso_set_callback( self, callback )
        class(motion_align_iso), intent(inout) :: self
        procedure(align_iso_callback)          :: callback
        self%callback => callback
    end subroutine motion_align_iso_set_callback

    subroutine motion_align_iso_set_frameweights_callback( self, frameweights_callback )
        class(motion_align_iso), intent(inout) :: self
        procedure(align_iso_fw_callback)       :: frameweights_callback
        self%frameweights_callback => frameweights_callback
    end subroutine motion_align_iso_set_frameweights_callback

end module simple_motion_align_iso
