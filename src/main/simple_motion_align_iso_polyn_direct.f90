! align frames in a maximum-likelihood way, with direct optimization of polynomial coefficients
! convergence can be determined by callback
! no frameweights are being applied here
module simple_motion_align_iso_polyn_direct
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_error
use simple_ft_expanded,           only: ft_expanded
use simple_image,                 only: image
use simple_parameters,            only: params_glob
use simple_opt_lbfgsb,            only: PRINT_NEVALS
use simple_opt_spec,              only: opt_spec
use simple_optimizer,             only: optimizer
use CPlot2D_wrapper_module
implicit none
public :: motion_align_iso_polyn_direct, align_iso_callback_polyn, POLYDIM, POLYDIM2
private
#include "simple_local_flags.inc"

integer,  parameter :: POLYDIM            = 3
integer,  parameter :: POLYDIM2           = 2*3
real,     parameter :: ISO_POLYN_DIRECT_FTOL_DEF  = 1e-8
real,     parameter :: ISO_POLYN_DIRECT_GTOL_DEF  = 1e-8
real(dp), parameter :: ISO_POLYN_DIRECT_FACTR_DEF = 1d+3
real(dp), parameter :: ISO_POLYN_DIRECT_PGTOL_DEF = 1d-7
real,     parameter :: ISO_POLYN_DIRECT_TRS_DEF   = 40.

logical, parameter :: inc_or_not = .true. ! use increments ("velocities") instead of shifts

type :: motion_align_iso_polyn_direct
    private
    logical                                           :: existence = .false.
    type(image), pointer                              :: frames(:)                         !< pointer to stack of frames
    integer                                           :: nframes_allocd = 0                !< number of frames allocated
    integer                                           :: nframes        = 0                !< number of frames
    integer                                           :: iter           = 0                !< iteration number
    real                                              :: trs                               !< size of box constraint
    real                                              :: ftol                              !< tolerance parameter for minimizer
    real                                              :: gtol                              !< tolerance parameter for minimizer
    real(dp)                                          :: factr, pgtol                      !< tolerance parameters for LBFGSB convergence
    real                                              :: hp      = -1.                     !< high pass value
    real                                              :: lp      = -1.                     !< low pass value
    real                                              :: bfactor = -1.                     !< b-factor for alignment weights
    integer                                           :: coord_x, coord_y                  !< x,y coordinates for patch-based alignment callback
    real,              allocatable                    :: opt_shifts(:,:)                   !< shifts identified
    real(dp),          allocatable, public            :: shifts(:,:)                       !< for contribs functions !****TODO**** remove
    real,              allocatable                    :: corrs(:)                          !< per-frame correlations
    real,              allocatable                    :: frameweights(:)                   !< weights per frame
    real                                              :: corr                              !< total correlation
    logical                                           :: group_frames = .false.            !< whether to group frames
    integer                                           :: lp_updates   = 0                  !< # of resolution updates performed [0;3]
    type(ft_expanded), allocatable                    :: movie_frames_Ij(:)
    type(ft_expanded), allocatable                    :: movie_frames_Ij_saved(:)          !< movie frames; saved for group frames
    type(ft_expanded), allocatable                    :: movie_frames_Ij_sh(:)             !< shifted movie frames
    type(ft_expanded), allocatable                    :: movie_frames_dIj_sh(:,:)          !< gradients of shifted movie frames
    type(ft_expanded), allocatable                    :: movie_frames_R                    !< reference movie frame
    type(ft_expanded), allocatable                    :: movie_frames_dR                   !< gradient of reference movie frame
    procedure(align_iso_callback_polyn), pointer, nopass :: callback           => null()   !< callback to determine convergence
contains
    procedure, private                                :: allocate_fields
    procedure, private                                :: deallocate_fields
    procedure, private                                :: coeffs_to_shifts          ! for polyn_direct
    procedure, private                                :: dcoeff                    ! for polyn_direct
    procedure, private                                :: vec_to_shifts             ! for direct
    procedure, private                                :: vec_to_shift_incrs        ! for direct`
    procedure, private                                :: shifts_to_vec             ! for direct
    procedure, private                                :: calc_corrs
    procedure, private                                :: calc_weights
    procedure                                         :: create_ftexp_objs
    procedure                                         :: refine_direct        => motion_align_iso_polyn_direct_refine_direct
    procedure                                         :: new                  => motion_align_iso_polyn_direct_new
    procedure                                         :: align_polyn          => motion_align_iso_polyn_direct_align_polyn
    procedure                                         :: kill                 => motion_align_iso_polyn_direct_kill
    procedure                                         :: set_frames           => motion_align_iso_polyn_direct_set_frames
    procedure                                         :: set_hp_lp            => motion_align_iso_polyn_direct_set_hp_lp
    procedure                                         :: set_trs              => motion_align_iso_polyn_direct_set_trs
    procedure                                         :: set_ftol_gtol        => motion_align_iso_polyn_direct_set_ftol_gtol
    procedure                                         :: set_factr_pgtol      => motion_align_iso_polyn_direct_set_factr_pgtol
    procedure                                         :: reset_tols           => motion_align_iso_polyn_direct_reset_tols
    procedure                                         :: get_corr             => motion_align_iso_polyn_direct_get_corr
    procedure                                         :: get_corrs            => motion_align_iso_polyn_direct_get_corrs
    procedure                                         :: get_opt_shifts       => motion_align_iso_polyn_direct_get_opt_shifts
    procedure                                         :: get_shifts_toplot    => motion_align_iso_polyn_direct_get_shifts_toplot
    procedure                                         :: set_coords           => motion_align_iso_polyn_direct_set_coords
    procedure                                         :: get_coords           => motion_align_iso_polyn_direct_get_coords
    procedure                                         :: set_callback         => motion_align_iso_polyn_direct_set_callback
    procedure                                         :: get_weights          => motion_align_iso_polyn_direct_get_weights
    procedure                                         :: set_group_frames     => motion_align_iso_polyn_direct_set_group_frames
    procedure                                         :: set_bfactor          => motion_align_iso_polyn_direct_set_bfactor
    procedure                                         :: motion_align_iso_polyn_direct_cost
    procedure                                         :: motion_align_iso_polyn_direct_gcost
    procedure                                         :: motion_align_iso_polyn_direct_fdf
    procedure                                         :: motion_align_iso_direct_cost
    procedure                                         :: motion_align_iso_direct_gcost
    procedure                                         :: motion_align_iso_direct_fdf
    procedure                                         :: motion_align_iso_contribs_cost
    procedure                                         :: motion_align_iso_contribs_gcost
    procedure                                         :: motion_align_iso_contribs_fdf
end type motion_align_iso_polyn_direct

abstract interface
    subroutine align_iso_callback_polyn(aptr, align_iso, converged)
        import motion_align_iso_polyn_direct
        class(*),                             intent(inout) :: aptr
        class(motion_align_iso_polyn_direct), intent(inout) :: align_iso
        logical,                              intent(out)   :: converged
    end subroutine align_iso_callback_polyn
end interface

contains

    subroutine motion_align_iso_polyn_direct_new( self )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        call self%kill()
        self%trs                   = ISO_POLYN_DIRECT_TRS_DEF
        self%ftol                  = ISO_POLYN_DIRECT_FTOL_DEF
        self%gtol                  = ISO_POLYN_DIRECT_GTOL_DEF
        self%factr                 = ISO_POLYN_DIRECT_FACTR_DEF
        self%pgtol                 = ISO_POLYN_DIRECT_PGTOL_DEF
        self%nframes_allocd        = 0
        self%callback              => null()
        self%hp                    = -1.
        self%lp                    = -1.
        self%bfactor               = -1.
        self%lp_updates            = 0
        self%group_frames          = .false.
        self%existence             = .true.
    end subroutine motion_align_iso_polyn_direct_new

    subroutine motion_align_iso_polyn_direct_align_polyn( self, callback_ptr )
        use simple_opt_factory, only: opt_factory
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        class(*),                             intent(inout) :: callback_ptr  !< callback pointer to be passed as first argument
        integer  :: iter
        real     :: cxy(3)
        logical  :: callback_convgd
        real     :: opt_lims(2*(self%nframes-1),2)
        real     :: lowest_cost
        logical  :: convergd
        real(dp) :: tmp_shifts(self%nframes, 2)
        type(opt_factory)         :: ofac
        type(opt_spec)            :: ospec                     !< optimizer specification object
        class(optimizer), pointer :: nlopt                     !< pointer to nonlinear optimizer
        write (*,*) ' ^^^^^^^^^^^^^^^^^^^^^^^^ align_polyn ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
        nlopt => null()
        if ( .not. self%existence ) then
            THROW_HARD('not instantiated; simple_motion_align_iso_polyn_direct: align')
        end if
        if ( self%nframes < 2 ) then
            THROW_HARD('nframes < 2; simple_motion_align_iso_polyn_direct: align')
        end if
        if (( self%hp < 0. ) .or. ( self%lp < 0.)) then
            THROW_HARD('hp or lp < 0; simple_motion_align_iso_polyn_direct: align')
        end if
        self%iter       = 0
        self%lp_updates = 0
        call self%create_ftexp_objs
        self%opt_shifts = 0.
        opt_lims(:,1) = -self%trs
        opt_lims(:,2) =  self%trs
        call ospec%specify('lbfgsb', (self%nframes-1)*2, ftol=self%ftol, gtol=self%gtol, limits=opt_lims, factr=self%factr, pgtol=self%pgtol)
        write (*,*) 'specifying: ftol=', self%ftol, 'gtol=', self%gtol, 'factr=', self%factr, 'pgtol=', self%pgtol
        call ofac%new(ospec, nlopt)
        call ospec%set_costfun_8(motion_align_iso_direct_cost_wrapper)
        call ospec%set_gcostfun_8(motion_align_iso_direct_gcost_wrapper)
        call ospec%set_fdfcostfun_8(motion_align_iso_direct_fdf_wrapper)
        !call ospec%set_costfun_8(motion_align_iso_polyn_direct_cost_wrapper)
        !call ospec%set_gcostfun_8(motion_align_iso_polyn_direct_gcost_wrapper)
        !call ospec%set_fdfcostfun_8(motion_align_iso_polyn_direct_fdf_wrapper)
        if (alloc_stat.ne.0) call allocchk('align 1; motion_align_iso_polyn_direct_align')
        convergd = .false.
        ospec%x_8 = 0._dp
        ospec%x   = 0.
        do while (.not. convergd)
            call nlopt%minimize(ospec, self, lowest_cost)
            iter = iter + 1
            call self%vec_to_shifts(ospec%x_8, tmp_shifts)!call self%coeffs_to_shifts(ospec%x_8, tmp_shifts)
            self%opt_shifts = - real(tmp_shifts)
            self%corr       = - real(self%motion_align_iso_direct_cost(ospec%x_8) / real(self%nframes,dp))
            if (associated(self%callback)) then
                call self%callback(callback_ptr, self, convergd)
            else
                convergd = .true.
            end if
            if (.not. convergd) then
                call self%create_ftexp_objs    ! lp should have changed, otherwise it should have converged
                self%lp_updates = self%lp_updates + 1
            end if
            ospec%ftol = self%ftol
            ospec%gtol = self%gtol
            ospec%pgtol = self%pgtol
            ospec%factr = self%factr
        end do
        call self%calc_corrs
        call self%calc_weights
        call nlopt%kill
        deallocate(nlopt)
    end subroutine motion_align_iso_polyn_direct_align_polyn

    subroutine motion_align_iso_polyn_direct_kill( self )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        call self%deallocate_fields
        self%nframes_allocd        =    0
        self%callback              =>   null()
        self%existence             =    .false.
    end subroutine motion_align_iso_polyn_direct_kill

    subroutine allocate_fields( self )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        allocate(&
            self%movie_frames_Ij(self%nframes),       &
            self%movie_frames_Ij_saved(self%nframes), &
            self%movie_frames_Ij_sh(self%nframes),    &
            self%movie_frames_dIj_sh(self%nframes,2), &
            self%movie_frames_R,                      &
            self%movie_frames_dR,                     &
            self%opt_shifts(self%nframes,2),          &
            self%shifts(self%nframes, 2),             &
            self%corrs(self%nframes),                 &
            self%frameweights(self%nframes)         )
        if(alloc_stat.ne.0)call allocchk('allocate_fields; simple_motion_align_iso_polyn_direct')
        self%nframes_allocd = self%nframes
    end subroutine allocate_fields

    subroutine deallocate_fields( self )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        integer :: i, j
        if (self%nframes_allocd > 0) then
            do i = 1, self%nframes_allocd
                call self%movie_frames_Ij(i)%kill
                call self%movie_frames_Ij_saved(i)%kill
                call self%movie_frames_Ij_sh(i)%kill
                do j = 1, 2
                    call self%movie_frames_dIj_sh(i,j)%kill
                end do
            end do
            call self%movie_frames_R%kill
            call self%movie_frames_dR%kill
            deallocate( self%movie_frames_Ij, self%movie_frames_Ij_saved, self%movie_frames_Ij_sh, &
                self%movie_frames_dIj_sh, self%movie_frames_R,                                     &
                self%movie_frames_dR, self%opt_shifts, self%shifts, self%corrs, self%frameweights )
        else
            if ( allocated(self%movie_frames_Ij)     .or. allocated(self%movie_frames_Ij_sh) .or. &
                 allocated(self%movie_frames_dIj_sh) .or. allocated(self%movie_frames_R)     .or. &
                 allocated(self%movie_frames_dR)     .or. allocated(self%opt_shifts)         .or. &
                 allocated(self%shifts)              .or. allocated(self%corrs)              .or. &
                 allocated(self%movie_frames_Ij_saved)                                       .or. &
                 allocated(self%frameweights)                                                      ) then
                THROW_HARD('inconsistency; simple_motion_align_iso_polyn_direct: deallocate_fields')
            end if
        end if
    end subroutine deallocate_fields

    subroutine vec_to_shifts( self, vec, ashifts )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp),                             intent(in)    :: vec(:)
        real(dp),                             intent(out)   :: ashifts(:,:)
        integer  :: j
        ashifts(1,:) = 0.
        do j = 1, self%nframes-1
            if (inc_or_not) then
                ashifts(j+1,1) = ashifts(j,1) + vec(2*(j-1)+1)
                ashifts(j+1,2) = ashifts(j,2) + vec(2*(j-1)+2)
            else
                ashifts(j+1,1) = vec(2*(j-1)+1)
                ashifts(j+1,2) = vec(2*(j-1)+2)
            end if
        end do
    end subroutine vec_to_shifts

    subroutine vec_to_shift_incrs( self, vec, ashift_incrs )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp),                             intent(in)    :: vec(:)
        real(dp),                             intent(out)   :: ashift_incrs(:,:)
        integer  :: j
        ashift_incrs(1,:) = 0.
        do j = 1, self%nframes-1
            ashift_incrs(j+1,1) = vec(2*(j-1)+1)
            ashift_incrs(j+1,2) = vec(2*(j-1)+2)
        end do
    end subroutine vec_to_shift_incrs


    subroutine shifts_to_vec( self, ashifts, vec )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real,                                 intent(in)    :: ashifts(:,:)
        real,                                 intent(out)   :: vec(:)
        integer :: j
        ! we should assume that ashifts(1,:) = 0, but just in case, we subtract ashifts(1,:) anyway
        do j = 1, self%nframes-1
            if (inc_or_not) then
                vec(               j) = ashifts(j+1,1) - ashifts(j,1)
                vec(self%nframes-1+j) = ashifts(j+1,2) - ashifts(j,2)
            else
                vec(               j) = ashifts(j+1,1) - ashifts(1,1)
                vec(self%nframes-1+j) = ashifts(j+1,2) - ashifts(1,2)
            end if
            write (*,*) 'shifts_to_vec: accessing elements ', j, 'and', self%nframes-1+j
        end do
    end subroutine shifts_to_vec

    subroutine create_ftexp_objs( self )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        integer :: j, xy
        real    :: w, sumw
        !$omp parallel do default(shared) private(j,xy) schedule(static) proc_bind(close)
        do j=1,self%nframes
            if( self%group_frames )then
                call self%movie_frames_Ij_saved(j)%new(self%frames(j), self%hp, self%lp, .true., bfac=self%bfactor)
                call self%movie_frames_Ij      (j)%new(self%frames(j), self%hp, self%lp, .false.)
            else
                call self%movie_frames_Ij(j)%new(self%frames(j), self%hp, self%lp, .true., bfac=self%bfactor)
                call self%movie_frames_Ij(j)%normalize_mat
            end if
            call self%movie_frames_Ij_sh(j)%new(self%frames(j), self%hp, self%lp, .false.)
            do xy = 1,2
                call self%movie_frames_dIj_sh(j, xy)%new(self%frames(j), self%hp, self%lp, .false.)
            end do
        end do
        !$omp end parallel do
        call self%movie_frames_R%new(self%frames(1),  self%hp, self%lp, .false.)
        call self%movie_frames_dR%new(self%frames(1), self%hp, self%lp, .false.)
        if( self%group_frames )then
            w = 0.5 * exp(-real(self%lp_updates))
            do j=1,self%nframes
                sumw = 1.
                if( j == 1 )then
                    sumw = sumw + 1.2*w
                    call add_shifted_frame(1,2,    w/sumw)
                    call add_shifted_frame(1,3,0.2*w/sumw)
                elseif( j == self%nframes )then
                    sumw = sumw + 1.2*w
                    call add_shifted_frame(self%nframes,self%nframes-1,    w/sumw)
                    call add_shifted_frame(self%nframes,self%nframes-2,0.2*w/sumw)
                else
                    sumw = sumw + 2.*w
                    call add_shifted_frame(j,j-1,w/sumw)
                    call add_shifted_frame(j,j+1,w/sumw)
                end if
                call self%movie_frames_Ij(j)%add(self%movie_frames_Ij_saved(j), w=1./sumw)
                call self%movie_frames_Ij(j)%normalize_mat
            end do
        end if

    contains

        subroutine add_shifted_frame(ii, jj, w_here)
            integer, intent(in) :: ii,jj
            real,    intent(in) :: w_here
            real :: shvec(3)
            shvec(1:2) = self%opt_shifts(ii,:) - self%opt_shifts(jj,:)
            shvec(3)   = 0.
             call self%movie_frames_Ij_saved(jj)%shift(-shvec, self%movie_frames_Ij_sh(ii))
             ! shift the SAVED one by -shvec and store
             call self%movie_frames_Ij(ii)%add(self%movie_frames_Ij_sh(ii),w=w_here)
            ! add to the one to be update
        end subroutine add_shifted_frame

    end subroutine create_ftexp_objs

    subroutine coeffs_to_shifts( self, vec, ashifts )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp),                             intent(in)    :: vec(:)
        real(dp),                             intent(out)   :: ashifts(:,:)
        real(dp) :: acoeffs(POLYDIM,2)
        real(dp) :: tmp, t
        integer  :: j, xy
        acoeffs(:,1) = vec(        1:  POLYDIM)
        acoeffs(:,2) = vec(POLYDIM+1:2*POLYDIM)
        do j = 1, self%nframes
            t = real(j - 1, dp)
            do xy = 1, 2
                tmp =       t    * acoeffs(1,xy)
                tmp = tmp + t**2 * acoeffs(2,xy)
                tmp = tmp + t**3 * acoeffs(3,xy)
                ashifts(j,xy) = tmp
            end do
        end do
    end subroutine coeffs_to_shifts

    function dcoeff(self, j, p) result( res )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        integer,                              intent(in)    :: j, p
        real(dp) :: res
        real(dp) :: t
        t = real(j-1, dp)
        if (p == 1) then
            res = t
        else if (p == 2) then
            res = t**2
        else
            res = t**3
        end if
    end function dcoeff

    subroutine calc_corrs( self )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real    :: shvec(3)
        integer :: j
        real    :: N_R
        ! shifts frames
        do j = 1, self%nframes
            shvec(  3) = 0._dp
            shvec(1:2) = - self%opt_shifts(j, :)
            call self%movie_frames_Ij(j)%shift( shvec, self%movie_frames_Ij_sh(j) )
        end do
        ! R = sum_j I_j
        call self%movie_frames_R%zero
        do j = 1, self%nframes
            call self%movie_frames_R%add( self%movie_frames_Ij_sh(j) )
        end do
        N_R = sqrt(self%movie_frames_R%corr_unnorm(self%movie_frames_R))
        do j = 1, self%nframes
            self%corrs(j) = self%movie_frames_R%corr_unnorm(self%movie_frames_Ij_sh(j)) / N_R
        end do
    end subroutine calc_corrs

    subroutine calc_weights( self )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        if( params_glob%l_rankw )then
            self%frameweights = corrs2weights(self%corrs, params_glob%ccw_crit, params_glob%rankw_crit)
        else
            self%frameweights = corrs2weights(self%corrs, params_glob%ccw_crit)
        endif
    end subroutine calc_weights

    subroutine motion_align_iso_polyn_direct_set_frames( self, frames_ptr, nframes )
        class(motion_align_iso_polyn_direct),          intent(inout) :: self
        type(image), allocatable, target,              intent(in)    :: frames_ptr(:)
        integer,                                       intent(in)    :: nframes
        logical :: do_alloc
        if (.not. self%existence) then
            THROW_HARD('not instantiated; simple_motion_align_iso_polyn_direct: set_frames')
        end if
        self%frames  => frames_ptr
        self%nframes = nframes
        if (size(frames_ptr, 1) < nframes) then
            THROW_HARD('nframes > #frames provided; simple_motion_align_iso_polyn_direct: set_frames')
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
    end subroutine motion_align_iso_polyn_direct_set_frames

    subroutine motion_align_iso_polyn_direct_set_hp_lp( self, hp, lp )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real,                                 intent(in)    :: hp, lp
        if (.not. self%existence) then
            THROW_HARD('not instantiated; simple_motion_align_iso_polyn_direct: set_hp_lp')
        end if
        self%hp = hp
        self%lp = lp
        if (.not. self%nframes > 0) then
            THROW_HARD('nframes < 1; simple_motion_align_iso_polyn_direct: set_hp_lp')
        end if
    end subroutine motion_align_iso_polyn_direct_set_hp_lp

    subroutine motion_align_iso_polyn_direct_set_trs( self, trs )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real,                                 intent(in)    :: trs
        self%trs = trs
    end subroutine motion_align_iso_polyn_direct_set_trs

    subroutine motion_align_iso_polyn_direct_set_ftol_gtol( self, ftol, gtol )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real,                                 intent(in)    :: ftol, gtol
        self%ftol = ftol
        self%gtol = gtol
    end subroutine motion_align_iso_polyn_direct_set_ftol_gtol

    subroutine motion_align_iso_polyn_direct_set_factr_pgtol( self, factr, pgtol )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp),                             intent(in)    :: factr, pgtol
        self%factr = factr
        self%pgtol = pgtol
    end subroutine motion_align_iso_polyn_direct_set_factr_pgtol

    subroutine motion_align_iso_polyn_direct_reset_tols( self )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        self%ftol = ISO_POLYN_DIRECT_FTOL_DEF
        self%gtol = ISO_POLYN_DIRECT_GTOL_DEF
        self%factr = ISO_POLYN_DIRECT_FACTR_DEF
        self%pgtol = ISO_POLYN_DIRECT_PGTOL_DEF
    end subroutine motion_align_iso_polyn_direct_reset_tols

    function motion_align_iso_polyn_direct_get_corr( self ) result( corr )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real :: corr
        corr = self%corr
    end function motion_align_iso_polyn_direct_get_corr

    subroutine motion_align_iso_polyn_direct_get_corrs( self, corrs )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real, allocatable,                    intent(out)   :: corrs(:)
        ! have to calculate corrs here
        allocate( corrs(self%nframes), source=self%corrs )
    end subroutine motion_align_iso_polyn_direct_get_corrs

    subroutine motion_align_iso_polyn_direct_get_opt_shifts( self, opt_shifts )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real, allocatable,                    intent(out)   :: opt_shifts(:,:)
        allocate( opt_shifts(self%nframes, 2), source=self%opt_shifts )
    end subroutine motion_align_iso_polyn_direct_get_opt_shifts

    subroutine motion_align_iso_polyn_direct_get_shifts_toplot( self, shifts_toplot )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real, allocatable,                    intent(out)   :: shifts_toplot(:,:)
        allocate( shifts_toplot(self%nframes, 2), source=self%opt_shifts )
    end subroutine motion_align_iso_polyn_direct_get_shifts_toplot


    subroutine motion_align_iso_polyn_direct_set_coords( self, x, y )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        integer,                              intent(in) :: x, y
        self%coord_x = x
        self%coord_y = y
    end subroutine motion_align_iso_polyn_direct_set_coords

    subroutine motion_align_iso_polyn_direct_get_coords(self, x, y )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        integer,                              intent(out) :: x, y
        x = self%coord_x
        y = self%coord_y
    end subroutine motion_align_iso_polyn_direct_get_coords

    subroutine motion_align_iso_polyn_direct_set_callback( self, callback )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        procedure(align_iso_callback_polyn) :: callback
        self%callback => callback
    end subroutine motion_align_iso_polyn_direct_set_callback

    subroutine motion_align_iso_polyn_direct_get_weights( self, frameweights )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real, allocatable,                    intent(out)   :: frameweights(:)
        allocate(frameweights(self%nframes), source=self%frameweights)
    end subroutine motion_align_iso_polyn_direct_get_weights

    subroutine motion_align_iso_polyn_direct_set_group_frames( self, group_frames )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        logical,                              intent(in)    :: group_frames
        self%group_frames = group_frames
    end subroutine motion_align_iso_polyn_direct_set_group_frames

    subroutine motion_align_iso_polyn_direct_set_bfactor( self, bfac )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real,                                 intent(in)    :: bfac
        self%bfactor = bfac
    end subroutine motion_align_iso_polyn_direct_set_bfactor

    subroutine motion_align_iso_polyn_direct_refine_direct( self )
        use simple_opt_factory, only: opt_factory
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp)                  :: tmp_shifts(self%nframes,2)
        type(opt_factory)         :: ofac
        type(opt_spec)            :: ospec                     !< optimizer specification object
        class(optimizer), pointer :: nlopt                     !< pointer to nonlinear optimizer
        real :: opt_lims(2*(self%nframes-1),2)
        real :: lowest_cost
        write (*,*) ' ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ refine_direct ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
        opt_lims(:,1) = - self%trs
        opt_lims(:,2) =   self%trs
        write (*,*) 'self%nframes: ', self%nframes, '(self%nframes-1)*2: ', (self%nframes-1)*2
        call ospec%specify('lbfgsb', (self%nframes-1)*2, ftol=self%ftol, gtol=self%gtol, limits=opt_lims, factr=self%factr, pgtol=self%pgtol)
        call ofac%new(ospec, nlopt)
        call ospec%set_costfun_8(motion_align_iso_direct_cost_wrapper)
        call ospec%set_gcostfun_8(motion_align_iso_direct_gcost_wrapper)
        call ospec%set_fdfcostfun_8(motion_align_iso_direct_fdf_wrapper)
        call self%shifts_to_vec(self%opt_shifts, ospec%x)
        ospec%x   = - ospec%x
        ospec%x_8 = real(ospec%x, dp)
        ! minimize
        call nlopt%minimize(ospec, self, lowest_cost)
        call self%vec_to_shifts(ospec%x_8, tmp_shifts)
        self%opt_shifts = -real(tmp_shifts)
        call self%calc_corrs
        call self%calc_weights
        ! cleanup
        call nlopt%kill
        deallocate(nlopt)
    end subroutine motion_align_iso_polyn_direct_refine_direct

    function motion_align_iso_polyn_direct_cost( self, vec ) result( cost )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp),                             intent(in)    :: vec(POLYDIM2)
        real(dp) :: shifts(self%nframes,2), shvec(3)
        real(dp) :: cost
        integer :: j
        ! determine shifts
        call self%coeffs_to_shifts(vec,shifts)
        ! shifts frames
        !$omp parallel do default(shared) private(j,shvec) schedule(static) proc_bind(close)
        do j = 1, self%nframes
            shvec(1:2) = shifts(j,:)
            shvec(3) = 0.
            call self%movie_frames_Ij(j)%shift( real(shvec), self%movie_frames_Ij_sh(j) )
        end do
        !$omp end parallel do
        ! R = sum_j I_j
        call self%movie_frames_R%zero
        do j = 1, self%nframes                    !! TODO: try openmp reduction (accuracy lost significant?)
            call self%movie_frames_R%add( self%movie_frames_Ij_sh(j) )
        end do
        ! corr = sqrt( R * R )
        cost = dsqrt(self%movie_frames_R%corr_unnorm(self%movie_frames_R))
        cost = - cost
    end function motion_align_iso_polyn_direct_cost

    function motion_align_iso_polyn_direct_cost_wrapper( self, vec, D ) result( cost )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(dp),     intent(in)    :: vec(D)
        real(dp) :: cost
        select type(self)
        class is (motion_align_iso_polyn_direct)
            cost = self%motion_align_iso_polyn_direct_cost( vec )
        class DEFAULT
            THROW_HARD('unknown type; simple_motion_align_iso_polyn_direct: cost_wrapper')
        end select
    end function motion_align_iso_polyn_direct_cost_wrapper

    subroutine motion_align_iso_polyn_direct_gcost( self, vec, grad )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp),                             intent(in)    :: vec(POLYDIM2)
        real(dp),                             intent(out)   :: grad(POLYDIM2)
        real(dp) :: shifts(self%nframes,2), shvec(3)
        real(dp) :: cost, w, atmp
        integer  :: j, p, xy
        ! determine shifts
        call self%coeffs_to_shifts(vec,shifts)
        ! shifts frames
        !$omp parallel do default(shared) private(j,shvec) schedule(static) proc_bind(close)
        do j = 1, self%nframes
            shvec(3) = 0._dp
            shvec(1:2) = shifts(j,:)
            call self%movie_frames_Ij(j)%shift( real(shvec), self%movie_frames_Ij_sh(j) )
            call self%movie_frames_Ij(j)%gen_grad( shvec, self%movie_frames_dIj_sh(j,1), &
                                                          self%movie_frames_dIj_sh(j,2) )
        end do
        !$omp end parallel do
        ! R = sum_j I_j
        call self%movie_frames_R%zero
        do j = 1, self%nframes
            call self%movie_frames_R%add( self%movie_frames_Ij_sh(j) )
        end do
        ! corr = sqrt( R * R )
        cost = sqrt(self%movie_frames_R%corr_unnorm(self%movie_frames_R))
        do p = 1, POLYDIM
            do xy = 1, 2
                call self%movie_frames_dR%zero
                do j = 1, self%nframes
                    w = self%dcoeff(j,p)
                    call self%movie_frames_dR%add( self%movie_frames_dIj_sh(j, xy), real(w) )
                end do
                atmp = self%movie_frames_R%corr_unnorm(self%movie_frames_dR) / cost
                grad((xy-1)*POLYDIM+p) = atmp
            end do
        end do
        grad = - grad
    end subroutine motion_align_iso_polyn_direct_gcost

    subroutine motion_align_iso_polyn_direct_gcost_wrapper( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        grad = 0.d0
        select type(self)
        class is (motion_align_iso_polyn_direct)
            call self%motion_align_iso_polyn_direct_gcost( vec, grad )
        class DEFAULT
            THROW_HARD('unknown type; simple_motion_align_iso_polyn_direct: gcost_wrapper')
        end select
    end subroutine motion_align_iso_polyn_direct_gcost_wrapper

    subroutine motion_align_iso_polyn_direct_fdf( self, vec, f, grad )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp),                             intent(in)    :: vec(POLYDIM2)
        real(dp),                             intent(out)   :: f
        real(dp),                             intent(out)   :: grad(POLYDIM2)
        real(dp) :: shifts(self%nframes,2), shvec(3)
        real(dp) :: w, atmp
        integer  :: j, p, xy
        ! determine shifts
        write (*,*) 'fdf, polyn_direct, vec=', vec
        call self%coeffs_to_shifts(vec,shifts)
        ! shifts frames
        !$omp parallel do default(shared) private(j,shvec) schedule(static) proc_bind(close)
        do j = 1, self%nframes
            shvec(  3) = 0._dp
            shvec(1:2) = shifts(j,:)
            call self%movie_frames_Ij(j)%shift( real(shvec), self%movie_frames_Ij_sh(j) )
            call self%movie_frames_Ij(j)%gen_grad( shvec, self%movie_frames_dIj_sh(j,1), &
                                                          self%movie_frames_dIj_sh(j,2) )
        end do
        !$omp end parallel do
        ! R = sum_j I_j
        call self%movie_frames_R%zero
        do j = 1, self%nframes                         !! TODO: try openmp reduction (accuracy lost significant?)
            call self%movie_frames_R%add( self%movie_frames_Ij_sh(j) )
        end do
        ! corr = sqrt( R * R )
        f = sqrt(self%movie_frames_R%corr_unnorm(self%movie_frames_R))
        do p = 1, POLYDIM
            do xy = 1, 2
                call self%movie_frames_dR%zero
                do j = 1, self%nframes                 !! TODO: try openmp reduction (accuracy lost significant?)
                    w = self%dcoeff(j,p)
                    call self%movie_frames_dR%add( self%movie_frames_dIj_sh(j, xy), real(w) )
                end do
                atmp = self%movie_frames_R%corr_unnorm(self%movie_frames_dR) / f
                grad((xy-1)*POLYDIM+p) = atmp
            end do
        end do
        f    = - f
        grad = - grad
    end subroutine motion_align_iso_polyn_direct_fdf

    subroutine motion_align_iso_polyn_direct_fdf_wrapper( self, vec, f, grad, D )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
        f    = 0.d0
        grad = 0.d0
        select type(self)
        class is (motion_align_iso_polyn_direct)
            call self%motion_align_iso_polyn_direct_fdf( vec, f, grad )
        class DEFAULT
            THROW_HARD('unknown type; simple_motion_align_iso_polyn_direct: fdf_wrapper')
        end select
    end subroutine motion_align_iso_polyn_direct_fdf_wrapper

    function motion_align_iso_direct_cost( self, vec ) result( cost )
        class(motion_align_iso_polyn_direct), intent(inout)    :: self
        real(dp),                       intent(in) :: vec((self%nframes)*2)
        real(dp) :: shifts(self%nframes,2), shvec(3)
        real(dp) :: cost
        integer  :: j
        ! determine shifts
        call self%vec_to_shifts(vec, shifts)
        ! shifts frames
        !$omp parallel do default(shared) private(j,shvec) schedule(static) proc_bind(close)
        do j = 1, self%nframes
            shvec(1:2) = shifts(j,:)
            shvec(3) = 0.
            call self%movie_frames_Ij(j)%shift( real(shvec), self%movie_frames_Ij_sh(j) )
        end do
        !$omp end parallel do
        ! R = sum_j I_j
        call self%movie_frames_R%zero
        do j = 1, self%nframes
            call self%movie_frames_R%add( self%movie_frames_Ij_sh(j) )
        end do
        ! corr = sqrt( R * R )
        cost = sqrt(self%movie_frames_R%corr_unnorm(self%movie_frames_R))
        cost = - cost
    end function motion_align_iso_direct_cost

    function motion_align_iso_direct_cost_wrapper( self, vec, D ) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(in)    :: vec(D)
        real(dp) :: cost
        select type(self)
            class is (motion_align_iso_polyn_direct)
                cost = self%motion_align_iso_direct_cost( vec )
            class DEFAULT
                THROW_HARD('unknown type; motion_align_iso_direct: cost_wrapper')
        end select
    end function motion_align_iso_direct_cost_wrapper

    subroutine motion_align_iso_direct_gcost( self, vec, grad )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp),                       intent(in)    :: vec((self%nframes-1)*2)
        real(dp),                       intent(out)   :: grad((self%nframes-1)*2)
        real(dp) :: f
        call self%motion_align_iso_direct_fdf(vec, f, grad)
    end subroutine motion_align_iso_direct_gcost

    subroutine motion_align_iso_direct_gcost_wrapper( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        grad = 0.d0
        select type(self)
        class is (motion_align_iso_polyn_direct)
            call self%motion_align_iso_direct_gcost( vec, grad )
        class DEFAULT
            THROW_HARD('unknown type; motion_align_iso_direct: gcost_wrapper')
        end select
    end subroutine motion_align_iso_direct_gcost_wrapper

    subroutine motion_align_iso_direct_fdf( self, vec, f, grad )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp),                       intent(in)    :: vec((self%nframes-1)*2)
        real(dp),                       intent(out)   :: f
        real(dp),                       intent(out)   :: grad((self%nframes-1)*2)
        real(dp) :: grad_tmp((self%nframes-1)*2)
        real(dp) :: shifts(self%nframes,2), shvec(3)
        real(dp) :: atmp
        integer  :: j, p, xy
        write (*,*) 'fdf, direct, vec=', vec
        ! determine shifts
        call self%vec_to_shifts(vec, shifts)
        ! shifts frames
        !$omp parallel do default(shared) private(j,shvec) schedule(static) proc_bind(close)
        do j = 1, self%nframes
            shvec(3) = 0._dp
            shvec(1:2) = shifts(j,:)
            call self%movie_frames_Ij(j)%shift( real(shvec), self%movie_frames_Ij_sh(j) )
            call self%movie_frames_Ij(j)%gen_grad( shvec, self%movie_frames_dIj_sh(j,1), &
                self%movie_frames_dIj_sh(j,2) )
        end do
        !$omp end parallel do
        ! R = sum_j I_j
        call self%movie_frames_R%zero
        do j = 1, self%nframes
            call self%movie_frames_R%add( self%movie_frames_Ij_sh(j) )
        end do
        ! corr = sqrt( R * R )
        f = sqrt(self%movie_frames_R%corr_unnorm(self%movie_frames_R))
        !$omp parallel do collapse(2) default(shared) private(xy,j,atmp) schedule(static) proc_bind(close)
        do xy = 1, 2
            do j = 1, self%nframes-1 !parallelize here
                atmp = self%movie_frames_R%corr_unnorm(self%movie_frames_dIj_sh(j+1, xy)) / f
                grad_tmp(2*(j-1)+xy) = atmp
            end do
        end do
        !$omp end parallel do
        f    = - f
        grad_tmp = - grad_tmp
        if (inc_or_not) then
            grad = 0._dp
            do xy = 1, 2
                do j = self%nframes-1, 1, -1
                    if (j < self%nframes-1) then
                        grad(2*(j-1)+xy) = grad(2*(j)+xy)
                    end if
                    grad(2*(j-1)+xy) = grad(2*(j-1)+xy) + grad_tmp(2*(j-1)+xy)
                end do
            end do
        else
            grad  = grad_tmp
        end if
    end subroutine motion_align_iso_direct_fdf

    subroutine motion_align_iso_direct_fdf_wrapper( self, vec, f, grad, D )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
        f    = 0.d0
        grad = 0.d0
        select type(self)
            class is (motion_align_iso_polyn_direct)
                call self%motion_align_iso_direct_fdf( vec, f, grad )
            class DEFAULT
                THROW_HARD('unknown type; motion_align_iso_direct: fdf_wrapper')
        end select
    end subroutine motion_align_iso_direct_fdf_wrapper

    function motion_align_iso_contribs_cost( self ) result( cost )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp) :: shvec(3)
        real(dp) :: cost
        integer  :: j
        ! shifts will have been correctly determined by owning motion_patched object
        ! shifts frames
        do j = 1, self%nframes
            shvec(1:2) = self%shifts(j,:)
            shvec(  3) = 0.
            call self%movie_frames_Ij(j)%shift( real(shvec), self%movie_frames_Ij_sh(j) )
        end do
        ! R = sum_j I_j
        call self%movie_frames_R%zero
        do j = 1, self%nframes
            call self%movie_frames_R%add( self%movie_frames_Ij_sh(j) )
        end do
        ! corr = sqrt( R * R )
        cost = - real(sqrt(self%movie_frames_R%corr_unnorm(self%movie_frames_R)),dp)
    end function motion_align_iso_contribs_cost

    subroutine motion_align_iso_contribs_gcost( self, grad )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp),                             intent(out)   :: grad(POLYDIM2)
        real(dp) :: shvec(3)
        real(dp) :: cost, w, atmp, t
        integer  :: j, p, xy
        ! shifts will have been correctly determined by owning motion_patched object
        ! shifts frames
        do j = 1, self%nframes
            shvec(3) = 0._dp
            shvec(1:2) = self%shifts(j,:)
            call self%movie_frames_Ij(j)%shift( real(shvec), self%movie_frames_Ij_sh(j) )
            call self%movie_frames_Ij(j)%gen_grad( shvec, self%movie_frames_dIj_sh(j,1), &
                                                          self%movie_frames_dIj_sh(j,2) )
        end do
        ! R = sum_j I_j
        call self%movie_frames_R%zero
        do j = 1, self%nframes
            call self%movie_frames_R%add( self%movie_frames_Ij_sh(j) )
        end do
        ! corr = sqrt( R * R )
        cost = sqrt(self%movie_frames_R%corr_unnorm(self%movie_frames_R))
        do p = 1, POLYDIM
            do xy = 1, 2
                call self%movie_frames_dR%zero
                do j = 1, self%nframes
                    t = real(j - 1, dp)
                    w = t**p
                    call self%movie_frames_dR%add( self%movie_frames_dIj_sh(j, xy), real(w) )
                end do
                atmp = self%movie_frames_R%corr_unnorm(self%movie_frames_dR) / cost
                grad((xy-1)*POLYDIM+p) = atmp
            end do
        end do
        grad = - grad
    end subroutine motion_align_iso_contribs_gcost

    subroutine motion_align_iso_contribs_fdf( self, f, grad )
        class(motion_align_iso_polyn_direct), intent(inout) :: self
        real(dp),                               intent(out)   :: f
        real(dp),                               intent(out)   :: grad(POLYDIM2)
        real(dp) :: shvec(3)
        real(dp) :: w, t
        integer  :: j, p, xy
        ! shifts will have been correctly determined by owning motion_patched object
        ! shifts frames
!ss        !$omp parallel do default(shared) private(j) schedule(static) proc_bind(close)
        do j = 1, self%nframes
            shvec(3) = 0._dp
            shvec(1:2) = self%shifts(j,:)
            call self%movie_frames_Ij(j)%shift( real(shvec), self%movie_frames_Ij_sh(j) )
            call self%movie_frames_Ij(j)%gen_grad( shvec, self%movie_frames_dIj_sh(j,1), &
                                                          self%movie_frames_dIj_sh(j,2) )
        end do
!ss        !$omp end parallel do
        ! R = sum_j I_j
        call self%movie_frames_R%zero
        do j = 1, self%nframes
            call self%movie_frames_R%add( self%movie_frames_Ij_sh(j) )
        end do
        ! corr = sqrt( R * R )
        f = sqrt(self%movie_frames_R%corr_unnorm(self%movie_frames_R))
        do p = 1, POLYDIM
            do xy = 1, 2
                call self%movie_frames_dR%zero
                do j = 1, self%nframes
                    t = real(j - 1, dp)
                    w = t**p
                    call self%movie_frames_dR%add( self%movie_frames_dIj_sh(j, xy), real(w) )
                end do
                grad((xy-1)*POLYDIM+p) = self%movie_frames_R%corr_unnorm(self%movie_frames_dR) / f
            end do
        end do
        f    = - f
        grad = - grad
    end subroutine motion_align_iso_contribs_fdf

end module simple_motion_align_iso_polyn_direct
