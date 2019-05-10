! patched-based anisotropic motion correction through direct optimization of polynomial coefficients
module simple_motion_patched_direct
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,      only: params_glob
use simple_opt_lbfgsb,       only: PRINT_NEVALS
use simple_opt_factory,      only: opt_factory
use simple_opt_spec,         only: opt_spec
use simple_optimizer,        only: optimizer
use simple_image,            only: image
use simple_ft_expanded,      only: ft_expanded
use simple_ftexp_shsrch,     only: ftexp_shsrch
use CPlot2D_wrapper_module
implicit none
private
public :: motion_patched_direct
#include "simple_local_flags.inc"

! module global constants
integer, parameter :: NX_PATCHED           = 5   ! number of patches in x-direction
integer, parameter :: NY_PATCHED           = 5   !       "      "       y-direction
integer, parameter :: X_OVERLAP            = 0   ! number of overlapping pixels per patch in x-direction
integer, parameter :: Y_OVERLAP            = 0   !       "      "        "         "         y-direction
real,    parameter :: TOL                  = 1e-6 !< tolerance parameter
real,    parameter :: TRS_DEFAULT          = 7.
integer, parameter :: PATCH_PDIM           = 18  ! dimension of fitted polynomial
integer, parameter :: PATCH_PDIM_DBL       = 18*2
integer, parameter :: MITSREF_DEFAULT_DIR  = 70
integer, parameter :: MAXITS_DEFAULT_DIR   = 400
real,    parameter :: NIMPROVED_TOL        = 1e-7
real,    parameter :: COEFFS_BOUND_DEFAULT = 1.  ! TODO: make this more sensible; make this a vector and according to polynomial term order

type :: rmat_ptr_type
    real, pointer :: rmat_ptr(:,:,:)
end type rmat_ptr_type

type :: stack_type
    type(image), allocatable :: stack(:)
end type stack_type

type :: ftexp_stack_type
    type(ft_expanded), allocatable :: stack(:)
end type ftexp_stack_type

type :: shsrch_stack_type
    type(ftexp_shsrch), allocatable :: stack(:)
end type shsrch_stack_type

type :: cmat_type
    complex, allocatable :: mat(:,:,:)
end type cmat_type

type :: motion_patched_direct
    private
    logical                              :: existence
    type(stack_type),        allocatable :: frame_patches(:,:)
    real,                    allocatable :: shifts_patches(:,:,:,:)
    character(len=:),        allocatable :: shift_fname
    real,                    allocatable :: global_shifts(:,:)
    real,                    allocatable :: frameweights(:)
    logical                              :: has_global_shifts
    logical                              :: has_frameweights    = .false.
    integer                              :: nframes
    integer                              :: ldim(3)       ! size of entire frame, reference
    integer                              :: ldim_patch(3) ! size of one patch
    integer                              :: lims_patches (NX_PATCHED,NY_PATCHED,2,2) ! corners of the patches
    real                                 :: patch_centers(NX_PATCHED,NY_PATCHED,2)
    real                                 :: motion_correct_ftol
    real                                 :: motion_correct_gtol
    real(dp)                             :: poly_coeffs    (PATCH_PDIM,2)  ! coefficients of fitted polynomial
    real(dp)                             :: poly_coeffs_tmp(PATCH_PDIM,2)  ! coefficients of fitted polynomial
    real                                 :: trs
    real                                 :: hp
    real                                 :: lp
    integer                              :: updateres
    real                                 :: resstep
    type(ftexp_stack_type),  allocatable :: movie_frames_ftexp(:,:)             !< movie frames
    type(ftexp_stack_type),  allocatable :: movie_frames_ftexp_sh(:,:)          !< shifted movie frames
    type(ftexp_stack_type),  allocatable :: movie_sum_global_ftexp_threads(:,:) !< array of global movie sums for parallel refinement
    type(shsrch_stack_type), allocatable :: ftexp_shsrch_stack(:,:)
    type(opt_spec)                       :: ospec                     !< optimizer specification object
    class(optimizer),           pointer  :: nlopt     => null()       !< pointer to nonlinear optimizer
    real                                 :: coeffs_bound
    integer                              :: iter
    integer                              :: maxits
    real                                 :: corr
    real                                 :: corr_saved
    real                                 :: corr_prev
    integer                              :: mitsref
    real                                 :: frac_improved
    integer                              :: nimproved
    real,                    allocatable :: corrs(:)
    real(dp)                             :: opt_poly_saved(PATCH_PDIM,2)
    real,                    allocatable :: frameweights_saved(:)
    real                                 :: corrfrac
contains
    procedure, private                   :: allocate_fields
    procedure, private                   :: deallocate_fields
    procedure, private                   :: set_size_frames_ref
    procedure, private                   :: set_patches
    procedure, private                   :: det_shifts
    procedure, private                   :: get_local_shift
    procedure, private                   :: apply_polytransfo
    procedure, private                   :: plot_shifts
    procedure, private                   :: create_ftexp_objs
    procedure, private                   :: create_shsrch_objs
    procedure, private                   :: calc_frameweights
    procedure, private                   :: shift_wsum_and_calc_corrs
    procedure, private                   :: gen_sum_global_ftexp_threads
    procedure, private                   :: determine_nimproved
    procedure, private                   :: determine_convergence
    procedure, private                   :: shifts_patchcenters_from_poly
    procedure, private                   :: deriv_patch_poly
    procedure, private                   :: motion_patched_direct_set_frameweights
    procedure, private                   :: motion_patched_direct_set_mitsref
    procedure                            :: set_frameweights => motion_patched_direct_set_frameweights
    procedure                            :: set_mitsref      => motion_patched_direct_set_mitsref
    procedure                            :: new              => motion_patched_direct_new
    procedure                            :: correct          => motion_patched_direct_correct
    procedure                            :: kill             => motion_patched_direct_kill
    procedure                            :: motion_patched_costfun_8
    procedure                            :: motion_patched_gcostfun_8
    procedure                            :: motion_patched_fdfcostfun_8
end type motion_patched_direct


contains

    ! Polynomial for patch motion
    function patch_poly(p, n) result(res)
        real(dp), intent(in) :: p(:)
        integer,  intent(in) :: n
        real(dp) :: res(n)
        real(dp) :: x, y, t
        x = p(1)
        y = p(2)
        t = p(3)
        res(    1) = t
        res(    2) = t**2
        res(    3) = t**3
        res( 4: 6) = x * res( 1: 3)  ! x   * {t,t^2,t^3}
        res( 7: 9) = x * res( 4: 6)  ! x^2 * {t,t^2,t^3}
        res(10:12) = y * res( 1: 3)  ! y   * {t,t^2,t^3}
        res(13:15) = y * res(10:12)  ! y^2 * {t,t^2,t^3}
        res(16:18) = y * res( 4: 6)  ! x*y * {t,t^2,t^3}
    end function patch_poly

    function apply_patch_poly(c, x, y, t) result(res_sp)
        real(dp), intent(in) :: c(PATCH_PDIM), x, y, t
        real(sp) :: res_sp
        real(dp) :: res
        real(dp) :: x2, y2, xy, t2, t3
        x2 = x * x
        y2 = y * y
        xy = x * y
        t2 = t * t
        t3 = t2 * t
        res = 0._dp
        res = res + c( 1) * t      + c( 2) * t2      + c( 3) * t3
        res = res + c( 4) * t * x  + c( 5) * t2 * x  + c( 6) * t3 * x
        res = res + c( 7) * t * x2 + c( 8) * t2 * x2 + c( 9) * t3 * x2
        res = res + c(10) * t * y  + c(11) * t2 * y  + c(12) * t3 * y
        res = res + c(13) * t * y2 + c(14) * t2 * y2 + c(15) * t3 * y2
        res = res + c(16) * t * xy + c(17) * t2 * xy + c(18) * t3 * xy
        res_sp = real(res)
    end function apply_patch_poly

    subroutine deriv_patch_poly(self, x0, y0, iframe, deriv)
        ! TODO: memoize this
        class(motion_patched_direct), intent(inout) :: self
        real(dp), intent(in)  :: x0, y0
        integer,  intent(in)  :: iframe
        real(dp), intent(out) :: deriv(PATCH_PDIM)
        real(dp) :: t, x, y, x2, y2, xy, t2, t3
        t = real(iframe - 1)
        x = x0 / real(self%ldim(1)) - 0.5
        y = y0 / real(self%ldim(2)) - 0.5
        x2 = x * x
        y2 = y * y
        xy = x * y
        t2 = t * t
        t3 = t2 * t
        deriv( 1) = t
        deriv( 2) = t2
        deriv( 3) = t3
        deriv( 4) = t  * x
        deriv( 5) = t2 * x
        deriv( 6) = t3 * x
        deriv( 7) = t  * x2
        deriv( 8) = t2 * x2
        deriv( 9) = t3 * x2
        deriv(10) = t  * y
        deriv(11) = t2 * y
        deriv(12) = t3 * y
        deriv(13) = t  * y2
        deriv(14) = t2 * y2
        deriv(15) = t3 * y2
        deriv(16) = t  * xy
        deriv(17) = t2 * xy
        deriv(18) = t3 * xy
    end subroutine deriv_patch_poly

    subroutine plot_shifts(self)
        class(motion_patched_direct), intent(inout) :: self
        real, parameter       :: SCALE = 40.
        real                  :: shift_scale
        type(str4arr)         :: title
        type(CPlot2D_type)    :: plot2D
        type(CDataSet_type)   :: dataSetStart, dataSet        !!!!!! todo: we don't need this
        type(CDataSet_type)   :: fit
        type(CDataSet_type)   :: patch_start
        type(CDataPoint_type) :: point2, p_obs, p_fit, point
        real                  :: xcenter, ycenter
        integer               :: ipx, ipy, iframe, j
        real                  :: loc_shift(2)
        shift_scale = SCALE
        call CPlot2D__new(plot2D, self%shift_fname)
        call CPlot2D__SetXAxisSize(plot2D, 600._c_double)
        call CPlot2D__SetYAxisSize(plot2D, 600._c_double)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetFlipY(plot2D, C_TRUE)
        if (self%has_global_shifts) then
            call CDataSet__new(dataSet)
            call CDataSet__SetDrawMarker(dataSet, C_FALSE)
            call CDataSet__SetDatasetColor(dataSet, 0.0_c_double, 0.0_c_double, 1.0_c_double)
            xcenter = real(self%ldim(1))/2.
            ycenter = real(self%ldim(2))/2.
            do j = 1, self%nframes
                call CDataPoint__new2(&
                    real(xcenter + SCALE * self%global_shifts(j, 1), c_double), &
                    real(ycenter + SCALE * self%global_shifts(j, 2), c_double), &
                    point)
                call CDataSet__AddDataPoint(dataSet, point)
                call CDataPoint__delete(point)
            end do
            call CPlot2D__AddDataSet(plot2D, dataset)
            call CDataSet__delete(dataset)
        end if
        call CDataSet__new(dataSetStart)
        call CDataSet__SetDrawMarker(dataSetStart, C_TRUE)
        call CDataSet__SetMarkerSize(dataSetStart,5._c_double)
        call CDataSet__SetDatasetColor(dataSetStart, 1.0_c_double,0.0_c_double,0.0_c_double)
        call CDataPoint__new2(real(xcenter,c_double), real(ycenter,c_double), point2)
        call CDataSet__AddDataPoint(dataSetStart, point2)
        call CPlot2D__AddDataSet(plot2D, dataSetStart)
        do ipx = 1, NX_PATCHED
            do ipy = 1, NY_PATCHED
                call CDataSet__new(patch_start)
                call CDataSet__SetDrawMarker(patch_start,C_TRUE)
                call CDataSet__SetMarkerSize(patch_start,5.0_c_double)
                call CDataSet__SetDatasetColor(patch_start,1.0_c_double,0.0_c_double,0.0_c_double)
                call CDataSet__new(fit)
                call CDataSet__SetDrawMarker(fit, C_FALSE)
                call CDataSet__SetDatasetColor(fit, 0.0_c_double,0.0_c_double,0.0_c_double)
                do iframe = 1, self%nframes
                    call self%get_local_shift(iframe, self%patch_centers(ipx, ipy, 1), &
                        self%patch_centers(ipy, ipy, 2), loc_shift)
                    !loc_shift = -1. * loc_shift
                    call CDataPoint__new2(&
                        real(self%patch_centers(ipx, ipy, 1) + SCALE * loc_shift(1), c_double), &
                        real(self%patch_centers(ipx, ipy, 2) + SCALE * loc_shift(2) ,c_double), &
                        p_fit)
                    call CDataSet__AddDataPoint(fit, p_fit)
                    if (iframe == 1) then
                        call CDataSet__AddDataPoint(patch_start, p_fit)
                    end if
                    call CDataPoint__delete(p_fit)
                end do
                call CPlot2D__AddDataSet(plot2D, fit)
                call CPlot2D__AddDataSet(plot2D, patch_start)
                call CDataSet__delete(patch_start)
                call CDataSet__delete(fit)
            end do
        end do
        title%str = 'X (in pixels; trajectory scaled by ' // trim(real2str(SHIFT_SCALE)) // ')' // C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%str)
        title%str(1:1) = 'Y'
        call CPlot2D__SetYAxisTitle(plot2D, title%str)
        call CPlot2D__OutputPostScriptPlot(plot2D, self%shift_fname)
        call CPlot2D__delete(plot2D)
    end subroutine plot_shifts

    subroutine get_local_shift( self, iframe, x, y, shift )
        class(motion_patched_direct), intent(inout) :: self
        integer,               intent(in)  :: iframe
        real,                  intent(in)  :: x, y
        real,                  intent(out) :: shift(2)
        real :: t, xx, yy
        t  = real(iframe - 1)
        xx = x / real(self%ldim(1)) - 0.5
        yy = y / real(self%ldim(2)) - 0.5
        shift(1) = apply_patch_poly(self%poly_coeffs(:,1),real(xx,dp),real(yy,dp),real(t,dp))
        shift(2) = apply_patch_poly(self%poly_coeffs(:,2),real(xx,dp),real(yy,dp),real(t,dp))
    end subroutine get_local_shift

    subroutine apply_polytransfo( self, frames, frames_output )
        class(motion_patched_direct),    intent(inout) :: self
        type(image), allocatable, intent(inout) :: frames(:)
        type(image), allocatable, intent(inout) :: frames_output(:)
        integer :: i, j, iframe
        real    :: x, y, t
        real    :: x_trafo, y_trafo
        type(rmat_ptr_type) :: rmat_ins(self%nframes), rmat_outs(self%nframes)
        do iframe = 1, self%nframes
            call frames_output(iframe)%new(self%ldim, params_glob%smpd)
            if (frames(iframe)%is_ft()) call frames(iframe)%ifft()
            call frames(iframe)%get_rmat_ptr(rmat_ins(iframe)%rmat_ptr)
            call frames_output(iframe)%get_rmat_ptr(rmat_outs(iframe)%rmat_ptr)
        end do
        ! foo !$omp parallel do collapse(3) default(shared) private(iframe,j,i,x,y,t,x_trafo,y_trafo) proc_bind(close) schedule(static)
        do iframe = 1, self%nframes
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    t = real(iframe - 1)
                    x = real(i) / real(self%ldim(1)) - 0.5
                    y = real(j) / real(self%ldim(2)) - 0.5
                    x_trafo = real(i) - apply_patch_poly(self%poly_coeffs(:,1),real(x,dp),real(y,dp),real(t,dp))
                    y_trafo = real(j) - apply_patch_poly(self%poly_coeffs(:,2),real(x,dp),real(y,dp),real(t,dp))
                    rmat_outs(iframe)%rmat_ptr(i,j,1) = interp_bilin(x_trafo, y_trafo, iframe, rmat_ins )
                end do
            end do
        end do
        ! foo !$omp end parallel do
    contains

        pure real function interp_bilin( xval, yval, iiframe, rmat_ins2 )
            real,                intent(in)  :: xval, yval
            integer,             intent(in)  :: iiframe
            type(rmat_ptr_type), intent(in) :: rmat_ins2(self%nframes)
            integer  :: x1_h,  x2_h,  y1_h,  y2_h
            real     :: y1, y2, y3, y4, t, u
            logical :: outside
            outside = .false.
            x1_h = floor(xval)
            x2_h = x1_h + 1
            if( x1_h<1 .or. x2_h<1 )then
                x1_h    = 1
                outside = .true.
            endif
            if( x1_h>self%ldim(1) .or. x2_h>self%ldim(1) )then
                x1_h    = self%ldim(1)
                outside = .true.
            endif
            y1_h = floor(yval)
            y2_h = y1_h + 1
            if( y1_h<1 .or. y2_h<1 )then
                y1_h    = 1
                outside = .true.
            endif
            if( y1_h>self%ldim(2) .or. y2_h>self%ldim(2) )then
                y1_h    = self%ldim(2)
                outside = .true.
            endif
            if( outside )then
                interp_bilin = rmat_ins2(iiframe)%rmat_ptr(x1_h, y1_h, 1)
                return
            endif
            y1 = rmat_ins2(iiframe)%rmat_ptr(x1_h, y1_h, 1)
            y2 = rmat_ins2(iiframe)%rmat_ptr(x2_h, y1_h, 1)
            y3 = rmat_ins2(iiframe)%rmat_ptr(x2_h, y2_h, 1)
            y4 = rmat_ins2(iiframe)%rmat_ptr(x1_h, y2_h, 1)
            t   = xval - x1_h
            u   = yval - y1_h
            interp_bilin =  (1. - t) * (1. - u) * y1 + &
                        &t  * (1. - u) * y2 + &
                        &t  *       u  * y3 + &
                        &(1. - t) * u  * y4
        end function interp_bilin
    end subroutine apply_polytransfo

    subroutine allocate_fields( self )
        class(motion_patched_direct), intent(inout) :: self
        integer :: alloc_stat
        logical :: do_allocate
        integer :: i, j
        do_allocate = .true.
        if (allocated(self%shifts_patches)) then
            if (size(self%shifts_patches, dim=1) < self%nframes) then
                do_allocate = .false.
            else
                call self%deallocate_fields()
            end if
        end if
        if (do_allocate) then
            allocate(self%shifts_patches           (2, self%nframes, NX_PATCHED, NY_PATCHED),&
                self%frame_patches                 (                 NX_PATCHED, NY_PATCHED),&
                self%movie_frames_ftexp            (                 NX_PATCHED, NY_PATCHED),&
                self%movie_frames_ftexp_sh         (                 NX_PATCHED, NY_PATCHED),&
                self%movie_sum_global_ftexp_threads(                 NX_PATCHED, NY_PATCHED),&
                self%ftexp_shsrch_stack            (                 NX_PATCHED, NY_PATCHED),&
                self%corrs                         (   self%nframes                        ),&
                self%frameweights_saved            (   self%nframes                        ),&
                stat=alloc_stat )
            if (alloc_stat /= 0) call allocchk('allocate_fields 1; simple_motion_patched_direct')
            do i = 1, NX_PATCHED
                do j = 1, NY_PATCHED
                    allocate( self%frame_patches           (i, j)%stack(self%nframes),&
                        self%movie_frames_ftexp            (i, j)%stack(self%nframes),&
                        self%movie_frames_ftexp_sh         (i, j)%stack(self%nframes),&
                        self%movie_sum_global_ftexp_threads(i, j)%stack(self%nframes),&
                        self%ftexp_shsrch_stack            (i, j)%stack(self%nframes),&
                        stat=alloc_stat )
                        if (alloc_stat /= 0) call allocchk('allocate_fields 2; simple_motion_patched_direct')
                end do
            end do
        end if
    end subroutine allocate_fields

    subroutine deallocate_fields( self )
        class(motion_patched_direct), intent(inout) :: self
        integer :: iframe, i, j
        if (allocated(self%shifts_patches)    ) deallocate(self%shifts_patches    )
        if (allocated(self%shift_fname)       ) deallocate(self%shift_fname       )
        if (allocated(self%global_shifts)     ) deallocate(self%global_shifts     )
        if (allocated(self%frameweights)      ) deallocate(self%frameweights      )
        if (allocated(self%frameweights_saved)) deallocate(self%frameweights_saved)
        if (allocated(self%corrs)             ) deallocate(self%corrs             )
        if (allocated(self%frame_patches)) then
            do j = 1, NY_PATCHED
                do i = 1, NX_PATCHED
                    do iframe = 1, self%nframes
                        call self%frame_patches(i,j)%stack(iframe)%kill
                        call self%movie_frames_ftexp(i, j)%stack(iframe)%kill
                        call self%movie_frames_ftexp_sh(i, j)%stack(iframe)%kill
                        call self%movie_sum_global_ftexp_threads(i, j)%stack(iframe)%kill
                        call self%ftexp_shsrch_stack(i, j)%stack(iframe)%kill
                    end do
                    deallocate(self%frame_patches(i,j)%stack,            &
                        self%movie_frames_ftexp(i, j)%stack,             &
                        self%movie_frames_ftexp_sh(i, j)%stack,          &
                        self%movie_sum_global_ftexp_threads(i, j)%stack, &
                        self%ftexp_shsrch_stack(i, j)%stack              &
                        )
                end do
            end do
            deallocate(self%frame_patches, self%movie_frames_ftexp, self%movie_frames_ftexp_sh,&
                self%movie_sum_global_ftexp_threads, self%ftexp_shsrch_stack)
        end if
    end subroutine deallocate_fields

    subroutine set_size_frames_ref( self )
        class(motion_patched_direct), intent(inout) :: self
        integer :: ldim_nooverlap(2)
        integer :: i, j
        self%ldim_patch(1) = round2even(real(self%ldim(1)) / real(NX_PATCHED)) + 2 * X_OVERLAP
        self%ldim_patch(2) = round2even(real(self%ldim(2)) / real(NY_PATCHED)) + 2 * Y_OVERLAP
        self%ldim_patch(3) = 1
        ldim_nooverlap(1) = round2even(real(self%ldim(1)) / real(NX_PATCHED))
        ldim_nooverlap(2) = round2even(real(self%ldim(2)) / real(NY_PATCHED))
        do j = 1, NY_PATCHED
            do i = 1, NX_PATCHED
                self%lims_patches(i,j,1,1) = (i-1) * ldim_nooverlap(1) - X_OVERLAP + 1
                self%lims_patches(i,j,1,2) =  i    * ldim_nooverlap(1) + X_OVERLAP
                self%lims_patches(i,j,2,1) = (j-1) * ldim_nooverlap(2) - Y_OVERLAP + 1
                self%lims_patches(i,j,2,2) =  j    * ldim_nooverlap(2) + Y_OVERLAP
                self%patch_centers(i,j,1)  = real(self%lims_patches(i,j,1,1) + self%lims_patches(i,j,1,2)) / 2.
                self%patch_centers(i,j,2)  = real(self%lims_patches(i,j,2,1) + self%lims_patches(i,j,2,2)) / 2.
            end do
        end do
    end subroutine set_size_frames_ref

    subroutine set_patches( self, stack, patches )
        class(motion_patched_direct),          intent(inout) :: self
        type(image),       allocatable, intent(inout) :: stack(:)
        type(stack_type),  allocatable, intent(inout) :: patches(:,:)
        real, allocatable :: res(:)
        integer :: i, j, iframe, k, l, kk, ll
        integer :: ip, jp           ! ip, jp: i_patch, j_patch
        integer :: lims_patch(2,2)
        type(rmat_ptr_type) :: rmat_ptrs(self%nframes)
        real, pointer :: rmat_patch(:,:,:)
        do j = 1, NY_PATCHED
            do i = 1, NX_PATCHED
                do iframe=1,self%nframes
                    call patches(i,j)%stack(iframe)%new(self%ldim_patch, params_glob%smpd)
                end do
            end do
        end do
        do iframe=1,self%nframes
            call stack(iframe)%get_rmat_ptr(rmat_ptrs(iframe)%rmat_ptr)
        end do
        ! foo !$omp parallel do collapse(3) default(shared) private(iframe,j,i,lims_patch,rmat_patch,k,l,kk,ll,ip,jp) proc_bind(close) schedule(static)
        do iframe=1,self%nframes
            do j = 1, NY_PATCHED
                do i = 1, NX_PATCHED
                    call patches(i,j)%stack(iframe)%get_rmat_ptr(rmat_patch)
                    lims_patch(:,:) = self%lims_patches(i,j,:,:)
                    do k = lims_patch(1,1), lims_patch(1,2)
                        kk = k
                        if (kk < 1) then
                            kk = kk + self%ldim(1)
                        else if (kk > self%ldim(1)) then
                            kk = kk - self%ldim(1)
                        end if
                        ip = k - lims_patch(1,1) + 1
                        do l = lims_patch(2,1), lims_patch(2,2)
                            ll = l
                            if (ll < 1) then
                                ll = ll + self%ldim(2)
                            else if (ll > self%ldim(2)) then
                                ll = ll - self%ldim(2)
                            end if
                            jp = l - lims_patch(2,1) + 1
                            ! now copy the value
                            rmat_patch(ip,jp,1) = rmat_ptrs(iframe)%rmat_ptr(kk,ll,1)
                        end do
                    end do
                end do
            end do
        end do
        ! foo !$omp end parallel do
        ! updates high-pass according to new dimensions
        res = self%frame_patches(1,1)%stack(1)%get_res()
        self%hp = min(self%hp,res(1))
        deallocate(res)
        write(logfhandle,'(A,F6.1)')'>>> PATCH HIGH-PASS: ',self%hp
    end subroutine set_patches

    subroutine det_shifts( self )
        class(motion_patched_direct), target, intent(inout) :: self
        integer           :: iframe, ipx, ipy
        integer           :: alloc_stat
        type(opt_factory) :: ofac
        real              :: opt_lims(PATCH_PDIM_DBL,2)
        integer           :: nimproved
        real              :: corr_new
        logical           :: convgd
        real              :: hp_saved, lp_saved
        integer           :: maxits_saved
        integer           :: iter
        self%shifts_patches = 0.
        if (alloc_stat /= 0) call allocchk('det_shifts 1; simple_motion_patched_direct')
        opt_lims(:,1) = - self%coeffs_bound
        opt_lims(:,2) =   self%coeffs_bound
        call self%ospec%specify('lbfgsb', PATCH_PDIM_DBL, ftol=self%motion_correct_ftol, &
            gtol=self%motion_correct_gtol, limits=opt_lims)
        call self%ospec%set_costfun_8   ( motion_patched_costfun_wrapper    )
        call self%ospec%set_gcostfun_8  ( motion_patched_gcostfun_wrapper   )
        call self%ospec%set_fdfcostfun_8( motion_patched_fdfcostfun_wrapper )
        ! generate optimizer object with the factory
        if( associated(self%nlopt) )then
            call self%nlopt%kill
            deallocate(self%nlopt)
        end if
        call ofac%new(self%ospec, self%nlopt)
        !initialize point
        self%ospec%x   = 0.
        self%ospec%x_8 = 0._dp
        self%iter      = 0
        call self%create_ftexp_objs
        call self%create_shsrch_objs
        call self%calc_frameweights
        call self%shifts_patchcenters_from_poly
        ! generate movie sum for refinement
        call self%shift_wsum_and_calc_corrs
        self%corr = sum(self%corrs)/real(self%nframes)
        call self%gen_sum_global_ftexp_threads
        do ipx = 1, NX_PATCHED
            do ipy = 1, NY_PATCHED
                do iframe = 1, self%nframes
                    call self%ftexp_shsrch_stack(ipx, ipy)%stack(iframe)%set_dims_and_alloc()
                    call self%ftexp_shsrch_stack(ipx, ipy)%stack(iframe)%calc_tmp_cmat12()
                end do
            end do
        end do
        if (self%maxits > 0) then
                self%ospec%maxits = self%maxits
            end if
            self%corr_saved = -1.
            do iter=1,self%mitsref
                self%iter = iter
                nimproved = 0
                call self%nlopt%minimize(self%ospec, self, corr_new)
                self%poly_coeffs(:,1) = self%ospec%x(           1:PATCH_PDIM    )
                self%poly_coeffs(:,2) = self%ospec%x(PATCH_PDIM+1:PATCH_PDIM_DBL)
                call self%shifts_patchcenters_from_poly
                call self%determine_nimproved ! also: corrs
                self%frac_improved = real(self%nimproved) / real(self%nframes)  * 100.
                write(logfhandle,'(a,1x,f4.0)') 'This % of frames improved their alignment: ', self%frac_improved
                call self%calc_frameweights
                ! build new reference
                call self%shift_wsum_and_calc_corrs
                call self%gen_sum_global_ftexp_threads
                self%corr_prev = self%corr
                self%corr      = sum(self%corrs) / real(self%nframes)
                if( self%corr >= self%corr_saved ) then ! save the local optimum
                    self%corr_saved         = self%corr
                    self%opt_poly_saved     = self%poly_coeffs
                    self%frameweights_saved = self%frameweights
                endif
                self%corrfrac = self%corr_prev / self%corr
                hp_saved = self%hp
                lp_saved = self%lp
                convgd = .false.
                maxits_saved = self%maxits
                call self%determine_convergence(convgd)
                if (maxits_saved /= self%maxits) then
                    self%ospec%maxits = self%maxits
                end if
                if (convgd) exit
                if ((abs(hp_saved-self%hp) > epsilon(hp_saved)) .or. &
                    (abs(lp_saved-self%lp) > epsilon(lp_saved))) then
                    ! need to re-make the ftexps
                    call self%create_ftexp_objs
                    call self%create_shsrch_objs
                    call self%calc_frameweights
                    call self%shift_wsum_and_calc_corrs
                    call self%gen_sum_global_ftexp_threads
                    do ipx = 1, NX_PATCHED
                        do ipy = 1, NY_PATCHED
                            do iframe = 1, self%nframes
                                call self%ftexp_shsrch_stack(ipx, ipy)%stack(iframe)%set_dims_and_alloc()
                            end do
                        end do
                    end do
                    ! need to destroy all previous knowledge about correlations
                    self%corr       = sum(self%corrs) / real(self%nframes)
                    self%corr_prev  = self%corr
                    self%corr_saved = self%corr
                end if
                do ipx = 1, NX_PATCHED
                    do ipy = 1, NY_PATCHED
                        do iframe = 1, self%nframes
                            call self%ftexp_shsrch_stack(ipx, ipy)%stack(iframe)%calc_tmp_cmat12()
                        end do
                    end do
                end do
            end do
            ! put the best local optimum back
            self%corr         = self%corr_saved
            self%poly_coeffs  = self%opt_poly_saved
            self%frameweights = self%frameweights_saved
    end subroutine det_shifts

    subroutine gen_sum_global_ftexp_threads( self )
        ! calculate the references by subtracting the individual frames (to reduce self-bias)
        class(motion_patched_direct), intent(inout) :: self
        integer :: ipx, ipy, iframe
        !$omp parallel do collapse(3) default(shared) private(ipx,ipy,iframe) proc_bind(close) schedule(static)
        do ipx = 1, NX_PATCHED
            do ipy = 1, NY_PATCHED
                do iframe=1,self%nframes
                    call self%movie_sum_global_ftexp_threads(ipx, ipy)%stack(iframe)%subtr(&
                        self%movie_frames_ftexp_sh(ipx, ipy)%stack(iframe), w=self%frameweights(iframe))
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine gen_sum_global_ftexp_threads

    subroutine determine_convergence( self, converged )
        class(motion_patched_direct), intent(inout) :: self
        logical,                      intent(out)   :: converged
        logical :: didupdateres
        didupdateres  = .false.
        select case(self%updateres)
        case(0)
            call update_res( 0.96, 40., self%updateres )
        case(1)
            call update_res( 0.97, 30., self%updateres )
        case(2)
            call update_res( 0.98, 20., self%updateres )
        case DEFAULT
            ! nothing to do
        end select
        if ( self%updateres > 2 .and. .not. didupdateres ) then ! at least one iteration with new lim
            if( self%nimproved == 0 .and. self%iter     > 2      )  converged = .true.
            if( self%iter      > 10 .and. self%corrfrac > 0.9999 )  converged = .true.
        else
            converged = .false.
        end if
    contains
        subroutine update_res( thres_corrfrac, thres_frac_improved, which_update )
            real,    intent(in) :: thres_corrfrac, thres_frac_improved
            integer, intent(in) :: which_update
            if( self%corrfrac > thres_corrfrac .and. self%frac_improved <= thres_frac_improved&
                .and. self%updateres == which_update )then
                self%lp = self%lp - self%resstep
                write(logfhandle,'(a,1x,f7.4)') '>>> LOW-PASS LIMIT UPDATED TO:', self%lp
                ! need to indicate that we updated resolution limit
                self%updateres  = self%updateres + 1
                ! indicate that reslim was updated
                didupdateres = .true.
            endif
        end subroutine update_res
    end subroutine determine_convergence

    subroutine shifts_patchcenters_from_poly( self )
        ! input: self%poly_coeffs; output: self%shifts_patches
        class(motion_patched_direct), intent(inout) :: self
        integer :: iframe, ipx, ipy
        real :: loc_shift(2)
        do iframe = 1, self%nframes
            do ipx = 1, NX_PATCHED
                do ipy = 1, NY_PATCHED
                    call self%get_local_shift(iframe, self%patch_centers(ipx, ipy, 1), &
                        self%patch_centers(ipx, ipy, 2), loc_shift)
                    self%shifts_patches(:,iframe,ipx,ipy) = loc_shift(:)
                end do
            end do
        end do
    end subroutine shifts_patchcenters_from_poly

    subroutine motion_patched_direct_set_frameweights( self, frameweights )
        class(motion_patched_direct), intent(inout) :: self
        real, allocatable, intent(in) :: frameweights(:)
        integer :: nlen
        nlen = size(frameweights)
        if (allocated(self%frameweights)) deallocate(self%frameweights)
        allocate(self%frameweights(nlen), source=frameweights)
        self%has_frameweights = .true.
    end subroutine motion_patched_direct_set_frameweights

    subroutine motion_patched_direct_set_mitsref( self, mitsref )
        class(motion_patched_direct), intent(inout) :: self
        integer,                      intent(in)    :: mitsref
        self%mitsref = mitsref
    end subroutine motion_patched_direct_set_mitsref

    subroutine motion_patched_direct_new( self, motion_correct_ftol, motion_correct_gtol, trs )
        class(motion_patched_direct), intent(inout) :: self
        real, optional,        intent(in)    :: motion_correct_ftol, motion_correct_gtol
        real, optional,        intent(in)    :: trs
        call self%kill()
        if (present(motion_correct_ftol)) then
            self%motion_correct_ftol = motion_correct_ftol
        else
            self%motion_correct_ftol = TOL
        end if
        if (present(motion_correct_gtol)) then
            self%motion_correct_gtol = motion_correct_gtol
        else
            self%motion_correct_gtol = TOL
        end if
        self%motion_correct_ftol = 1e-8
        self%motion_correct_gtol = 1e-8
        if (present(trs)) then
            self%trs = trs
        else
            self%trs = TRS_DEFAULT
        end if
        self%coeffs_bound  = COEFFS_BOUND_DEFAULT
        self%nimproved     = -1
        self%frac_improved = -1.
        self%hp            = -1.
        self%lp            = params_glob%lpstart
        self%updateres     = 0
        self%mitsref       = MITSREF_DEFAULT_DIR
        self%maxits        = MAXITS_DEFAULT_DIR
        self%existence = .true.
    end subroutine motion_patched_direct_new

    subroutine motion_patched_direct_correct( self, hp, resstep, frames, frames_output, shift_fname, &
        global_shifts )
        class(motion_patched_direct),  intent(inout) :: self
        real,                          intent(in)    :: hp, resstep
        type(image),      allocatable, intent(inout) :: frames(:)
        type(image),      allocatable, intent(inout) :: frames_output(:)
        character(len=:), allocatable, intent(in)    :: shift_fname
        real, optional,   allocatable, intent(in)    :: global_shifts(:,:)
        integer :: ldim_frames(3)
        integer :: i
        self%hp = hp
        self%lp = params_glob%lpstart
        self%resstep = resstep
        self%updateres = 0
        self%shift_fname = shift_fname // C_NULL_CHAR
        if (allocated(self%global_shifts)) deallocate(self%global_shifts)
        if (present(global_shifts)) then
            allocate(self%global_shifts(size(global_shifts, 1), size(global_shifts, 2)))
            self%global_shifts = global_shifts
            self%has_global_shifts = .true.
        else
            self%has_global_shifts = .false.
        end if
        self%nframes = size(frames,dim=1)
        self%ldim   = frames(1)%get_ldim()
        do i = 1,self%nframes
            ldim_frames = frames(i)%get_ldim()
            if (any(ldim_frames(1:2) /= self%ldim(1:2))) then
                THROW_HARD('error in motion_patched_direct_correct: frame dimensions do not match reference dimension; simple_motion_patched_direct')
            end if
        end do
        call self%allocate_fields()
        call self%set_size_frames_ref()
        ! divide the reference into patches
        call self%set_patches(frames, self%frame_patches)
        ! determine shifts for patches
        call self%det_shifts()
        ! apply transformation
        call self%apply_polytransfo(frames, frames_output)
        ! report visual results
        call self%plot_shifts()
    end subroutine motion_patched_direct_correct

    subroutine motion_patched_direct_kill( self )
        class(motion_patched_direct), intent(inout) :: self
        call self%deallocate_fields()
        if (allocated(self%frameweights)      ) deallocate(self%frameweights      )
        if (allocated(self%frameweights_saved)) deallocate(self%frameweights_saved)
        self%has_frameweights = .false.
        self%existence        = .false.
    end subroutine motion_patched_direct_kill

    subroutine create_ftexp_objs( self )
        class(motion_patched_direct), target, intent(inout) :: self
        type(image), pointer :: image_ptr
        integer :: iframe, ipx, ipy
        do iframe=1,self%nframes
            do ipx = 1, NX_PATCHED
                do ipy = 1, NY_PATCHED
                    image_ptr => self%frame_patches(ipx,ipy)%stack(iframe)
                    call self%movie_frames_ftexp(ipx, ipy)%            stack(iframe)%new(image_ptr, &
                        self%hp, self%lp, .true. )
                    call self%movie_frames_ftexp_sh(ipx, ipy)%         stack(iframe)%new(image_ptr, &
                        self%hp, self%lp, .false.)
                    call self%movie_sum_global_ftexp_threads(ipx, ipy)%stack(iframe)%new(image_ptr, &
                        self%hp, self%lp, .false.)
                end do
            end do
        end do
    end subroutine create_ftexp_objs

    subroutine create_shsrch_objs( self )
        class(motion_patched_direct), target, intent(inout) :: self
        integer :: iframe, ipx, ipy
        do iframe=1,self%nframes
            do ipx = 1, NX_PATCHED
                do ipy = 1, NY_PATCHED
                    call self%ftexp_shsrch_stack(ipx, ipy)%stack(iframe)%new(&
                        self%movie_sum_global_ftexp_threads(ipx, ipy)%stack(iframe),&
                        self%movie_frames_ftexp(ipx, ipy)%            stack(iframe),&
                        params_glob%trs&
                    )
                end do
            end do
        end do
    end subroutine create_shsrch_objs

    subroutine calc_frameweights( self )
        class(motion_patched_direct), intent(inout) :: self
        ! do nothing
    end subroutine calc_frameweights

    subroutine shift_wsum_and_calc_corrs( self )
        class(motion_patched_direct), target, intent(inout) :: self
        type(cmat_type), allocatable, target :: cmat_sums(:,:), cmats(:,:)
        complex, pointer :: cmat(:,:,:), cmat_sum(:,:,:)
        integer :: iframe, flims(3,2), ipx, ipy
        real    :: shvec(3)
        ! allocate matrices for reduction
        allocate(cmat_sums(NX_PATCHED, NY_PATCHED), cmats(NX_PATCHED, NY_PATCHED))
        do ipx = 1, NX_PATCHED
            do ipy = 1, NY_PATCHED
                flims = self%movie_sum_global_ftexp_threads(ipx,ipy)%stack(1)%get_flims()
                allocate(cmats(ipx,ipy)%mat(flims(1,1):flims(1,2),flims(2,1):flims(2,2),flims(3,1):flims(3,2)),&
                    cmat_sums(ipx,ipy)%mat(flims(1,1):flims(1,2),flims(2,1):flims(2,2),flims(3,1):flims(3,2)),&
                    source=cmplx(0.,0.))
            end do
        end do
        self%corrs = 0.
        do ipx = 1, NX_PATCHED
            do ipy = 1, NY_PATCHED
                cmat     => cmats(ipx,ipy)%mat
                cmat_sum => cmat_sums(ipx,ipy)%mat
                ! FIRST LOOP TO OBTAIN WEIGHTED SUM
                ! foo !$omp parallel default(shared) private(iframe,shvec,cmat) proc_bind(close)
                ! foo !$omp do schedule(static) reduction(+:cmat_sum)
                do iframe=1,self%nframes
                    shvec(1:2) = self%shifts_patches(:, iframe, ipx, ipy)
                    shvec(3)   = 0.
                    call self%movie_frames_ftexp(ipx,ipy)%stack(iframe)%shift(-shvec, &
                        self%movie_frames_ftexp_sh(ipx,ipy)%stack(iframe))
                    call self%movie_frames_ftexp_sh(ipx,ipy)%stack(iframe)%get_cmat(cmat)
                    cmat_sum = cmat_sum + cmat * self%frameweights(iframe)
                end do
                ! foo !$omp end do

                ! SECOND LOOP TO UPDATE movie_sum_global_ftexp_threads AND CALCULATE CORRS
                ! foo !$omp do schedule(static)
                do iframe=1,self%nframes
                    ! update array of sums (for future parallel exec)
                    call self%movie_sum_global_ftexp_threads(ipx,ipy)%stack(iframe)%set_cmat(cmat_sum)
                    ! subtract the movie frame being correlated to reduce bias
                    call self%movie_sum_global_ftexp_threads(ipx,ipy)%stack(iframe)%subtr(&
                        self%movie_frames_ftexp_sh(ipx,ipy)%stack(iframe), &
                        w=self%frameweights(iframe))
                    ! calc corr
                    self%corrs(iframe) = self%corrs(iframe) + self%movie_sum_global_ftexp_threads(ipx,ipy)%&
                        stack(iframe)%corr(self%movie_frames_ftexp_sh(ipx,ipy)%stack(iframe))
                    ! add the subtracted movie frame back to the weighted sum
                    call self%movie_sum_global_ftexp_threads(ipx,ipy)%stack(iframe)%add(&
                        self%movie_frames_ftexp_sh(ipx,ipy)%stack(iframe), &
                        w=self%frameweights(iframe))
                end do
                ! foo !$omp end do
                ! foo !$omp end parallel
            end do
        end do
        do ipx = 1, NX_PATCHED
            do ipy = 1, NY_PATCHED
                deallocate(cmats(ipx,ipy)%mat, cmat_sums(ipx,ipy)%mat)
            end do
        end do
        deallocate(cmats, cmat_sums)
    end subroutine shift_wsum_and_calc_corrs

    subroutine determine_nimproved( self )
        ! NOTE: could also consider individual patched instead of frames
        class(motion_patched_direct), intent(inout) :: self
        integer :: ipx, ipy, iframe
        real    :: shvec(3)
        real    :: corr_tmp, corr_tmp2
        integer :: nimproved
        nimproved = 0
        do iframe = 1, self%nframes
            corr_tmp = 0.
            do ipx = 1, NX_PATCHED
                do ipy = 1, NY_PATCHED
                    shvec(1:2) = self%shifts_patches(:, iframe, ipx, ipy)
                    shvec(3)   = 0.
                    call self%movie_frames_ftexp(ipx,ipy)%stack(iframe)%shift(-shvec, &
                        self%movie_frames_ftexp_sh(ipx,ipy)%stack(iframe))
                    corr_tmp2 = self%movie_sum_global_ftexp_threads(ipx,ipy)%&
                        stack(iframe)%corr(self%movie_frames_ftexp_sh(ipx,ipy)%stack(iframe))
                    corr_tmp  = corr_tmp + corr_tmp2
                end do
            end do
            if (corr_tmp - self%corrs(iframe) > NIMPROVED_TOL) nimproved = nimproved + 1
        end do
        self%nimproved = nimproved
    end subroutine determine_nimproved

    function motion_patched_costfun_8( self, vec ) result( cost )
        class(motion_patched_direct), intent(inout) :: self
        real(dp),                     intent(in)    :: vec(:)
        real(dp) :: cost
        real(dp) :: shvec(2)
        real(dp) :: cost_tmp
        integer  :: iframe, ipx, ipy
        self%poly_coeffs(:, 1) = vec(           1:  PATCH_PDIM)
        self%poly_coeffs(:, 2) = vec(PATCH_PDIM+1:2*PATCH_PDIM)
        call self%shifts_patchcenters_from_poly
        cost = 0._dp
        !$omp parallel do collapse(3) default(shared) private(ipx,ipy,iframe,shvec,cost_tmp) reduction(+:cost) proc_bind(close) schedule(static)
        do ipx = 1, NX_PATCHED
            do ipy = 1, NY_PATCHED
                do iframe = 1, self%nframes
                    shvec = self%shifts_patches(:,iframe,ipx,ipy)
                    cost_tmp = self%ftexp_shsrch_stack(ipx,ipy)%stack(iframe)%corr_shifted_cost_8(-shvec)
                    cost = cost + cost_tmp
                end do
            end do
        end do
        !$omp end parallel do
    end function motion_patched_costfun_8

    function motion_patched_costfun_wrapper( self, vec, D ) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(in)    :: vec(D)
        real(dp) :: cost
        select type(self)
        class is(motion_patched_direct)
            cost = -self%motion_patched_costfun_8( vec )
        class DEFAULT
            THROW_HARD('unknown type; motion_patched_costfun_wrapper')
        end select
    end function motion_patched_costfun_wrapper

    subroutine motion_patched_gcostfun_8( self, vec, grad )
        class(motion_patched_direct), intent(inout) :: self
        real(dp),                     intent(in)    :: vec(:)
        real(dp),                     intent(out)   :: grad(:)
        real(dp) :: grad_tmp(2)
        real(dp) :: shvec(2)
        real(dp) :: pgrad_tmp(PATCH_PDIM)
        real(dp) :: patch_cntr(2)
        integer  :: iframe, ipx, ipy
        self%poly_coeffs(:, 1) = vec(           1:  PATCH_PDIM)
        self%poly_coeffs(:, 2) = vec(PATCH_PDIM+1:2*PATCH_PDIM)
        call self%shifts_patchcenters_from_poly
        grad = 0._dp
        !$omp parallel do collapse(3) default(shared) private(ipx,ipy,iframe,shvec,grad_tmp,patch_cntr,pgrad_tmp) reduction(+:grad) proc_bind(close) schedule(static)
        do ipx = 1, NX_PATCHED
            do ipy = 1, NY_PATCHED
                do iframe = 1, self%nframes
                    shvec      = self%shifts_patches(:,iframe,ipx,ipy)
                    call self%ftexp_shsrch_stack(ipx,ipy)%stack(iframe)%corr_gshifted_cost_8(-shvec, grad_tmp)
                    patch_cntr = self%patch_centers(ipx, ipy, :)
                    call self%deriv_patch_poly(patch_cntr(1), patch_cntr(2), iframe, pgrad_tmp)
                    grad(           1:  PATCH_PDIM) = grad(           1:  PATCH_PDIM) + grad_tmp(1) * pgrad_tmp
                    grad(PATCH_PDIM+1:2*PATCH_PDIM) = grad(PATCH_PDIM+1:2*PATCH_PDIM) + grad_tmp(2) * pgrad_tmp
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine motion_patched_gcostfun_8

    subroutine motion_patched_gcostfun_wrapper( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        grad = 0.d0
        select type(self)
        class is (motion_patched_direct)
            call self%motion_patched_gcostfun_8( vec, grad )
        class DEFAULT
            THROW_HARD('unknown type; motion_patched_gcostfun_wrapper')
        end select
    end subroutine motion_patched_gcostfun_wrapper

    subroutine motion_patched_fdfcostfun_8( self, vec, f, grad )
        class(motion_patched_direct), intent(inout) :: self
        real(dp),                     intent(in)    :: vec(:)
        real(dp),                     intent(out)   :: f
        real(dp),                     intent(out)   :: grad(:)
        real(dp) :: f_tmp, grad_tmp(2)
        real(dp) :: shvec(2)
        real(dp) :: pgrad_tmp(PATCH_PDIM)
        real(dp) :: patch_cntr(2)
        integer  :: iframe, ipx, ipy
        self%poly_coeffs(:, 1) = vec(           1:  PATCH_PDIM)
        self%poly_coeffs(:, 2) = vec(PATCH_PDIM+1:2*PATCH_PDIM)
        call self%shifts_patchcenters_from_poly
        f    = 0._dp
        grad = 0._dp
        !$omp parallel do collapse(3) default(shared) private(ipx,ipy,iframe,shvec,f_tmp,grad_tmp,patch_cntr,pgrad_tmp) reduction(+:f,grad) proc_bind(close) schedule(static)
        do ipx = 1, NX_PATCHED
            do ipy = 1, NY_PATCHED
                do iframe = 1, self%nframes
                    shvec = self%shifts_patches(:,iframe,ipx,ipy)
                    call self%ftexp_shsrch_stack(ipx,ipy)%stack(iframe)%corr_fdfshifted_cost_8&
                        (-shvec, f_tmp, grad_tmp)
                    f    = f + f_tmp
                    patch_cntr = self%patch_centers(ipx, ipy, :)
                    call self%deriv_patch_poly(patch_cntr(1), patch_cntr(2), iframe, pgrad_tmp)
                    grad(           1:  PATCH_PDIM) = grad(           1:  PATCH_PDIM) + grad_tmp(1) * pgrad_tmp
                    grad(PATCH_PDIM+1:2*PATCH_PDIM) = grad(PATCH_PDIM+1:2*PATCH_PDIM) + grad_tmp(2) * pgrad_tmp
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine motion_patched_fdfcostfun_8

    subroutine motion_patched_fdfcostfun_wrapper( self, vec, f, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: f, grad(D)
        f    = 0.d0
        grad = 0.d0
        select type(self)
        class is (motion_patched_direct)
            call self%motion_patched_fdfcostfun_8( vec, f, grad )
            f    = f    * (-1._dp)
        class DEFAULT
            THROW_HARD('unknown type; motion_patched_fdfcostfun_wrapper')
        end select
    end subroutine motion_patched_fdfcostfun_wrapper


end module simple_motion_patched_direct
