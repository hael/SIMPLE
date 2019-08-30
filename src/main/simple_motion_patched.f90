! patched-based anisotropic motion correction
module simple_motion_patched
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,                    only: params_glob
use simple_opt_factory,                   only: opt_factory
use simple_opt_spec,                      only: opt_spec
use simple_optimizer,                     only: optimizer
use simple_image,                         only: image, imstack_type
use simple_ft_expanded,                   only: ft_expanded, ftexp_transfmat_init, ftexp_transfmat_kill
use simple_motion_align_hybrid,           only: motion_align_hybrid
use simple_motion_align_iso_polyn_direct, only: motion_align_iso_polyn_direct, POLYDIM, POLYDIM2
use CPlot2D_wrapper_module
implicit none
private
public :: motion_patched
#include "simple_local_flags.inc"

! module global constants
real,    parameter :: TOL            = 1e-6 !< tolerance parameter
real,    parameter :: TRS_DEFAULT    = 5.
integer, parameter :: PATCH_PDIM     = 18   ! dimension of fitted polynomial
real,    parameter :: DIRECT_FTOL    = 1e-7
real,    parameter :: DIRECT_GTOL    = 1e-7

type :: rmat_ptr_type
   real, pointer :: rmat_ptr(:,:,:)
end type rmat_ptr_type

type :: motion_patched
    private
    logical                             :: existence
    type(imstack_type),                  allocatable :: frame_patches(:,:)
    type(motion_align_iso_polyn_direct), allocatable :: align_iso_polyn_direct(:,:)
    real,                   allocatable :: shifts_patches(:,:,:,:)
    real,                   allocatable :: shifts_patches_for_fit(:,:,:,:)
    real,                   allocatable :: lp(:,:)
    real,                   allocatable :: global_shifts(:,:)       ! isotropic solution
    real,                   allocatable :: frameweights(:)
    real,                   allocatable :: patch_centers(:,:,:)     ! patches centers
    integer,                allocatable :: updateres(:,:)
    integer,                allocatable :: lims_patches(:,:,:,:)    ! patches corners
    character(len=:),       allocatable :: shift_fname
    integer                             :: nframes
    integer                             :: ldim(3)                  ! size of entire frame, reference
    integer                             :: ldim_patch(3)            ! size of one patch
    integer                             :: fixed_frame        = 1   ! frame of reference for alignment
    integer                             :: interp_fixed_frame = 1   ! frame of reference for interpolation
    real                                :: smpd               = 0.
    real                                :: bfactor            = -1.
    real                                :: motion_correct_ftol
    real                                :: motion_correct_gtol
    real                                :: motion_patched_direct_ftol
    real                                :: motion_patched_direct_gtol
    real(dp)                            :: poly_coeffs(PATCH_PDIM,2)  ! coefficients of fitted polynomial
    real                                :: poly_chisq(2)              ! polynomial fit goodness of fit
    real                                :: trs
    real, public                        :: hp
    real                                :: resstep
    logical                             :: has_global_shifts
    logical                             :: has_frameweights  = .false.
    logical                             :: fitshifts         = .false.
contains
    procedure                           :: new
    procedure, private                  :: allocate_fields
    procedure, private                  :: deallocate_fields
    procedure, private                  :: set_size_frames_ref
    procedure, private                  :: gen_patch
    procedure, private                  :: set_patches
    procedure, private                  :: det_shifts
    procedure, private                  :: det_shifts_polyn
    procedure, private                  :: det_shifts_direct
    procedure, private                  :: fit_polynomial
    procedure, private                  :: get_local_shift
    procedure, private                  :: apply_polytransfo
    procedure, private                  :: plot_shifts
    procedure, private                  :: motion_patched_polyn_callback
    procedure, private                  :: pix2polycoords
    procedure, private                  :: cleanup_polyn
    procedure, private                  :: grad_contrib_to_full_grad
    procedure, private                  :: motion_patched_direct_cost
    procedure, private                  :: motion_patched_direct_gcost
    procedure, private                  :: motion_patched_direct_fdf
    procedure, private                  :: coeffs_to_shifts
    procedure                           :: set_frameweights
    procedure                           :: set_fitshifts
    procedure                           :: set_fixed_frame
    procedure                           :: set_interp_fixed_frame
    procedure                           :: set_bfactor
    procedure                           :: get_poly4star
    procedure                           :: get_polyfit_chisq
    procedure                           :: polytransfo
    procedure                           :: correct         => motion_patched_correct
    procedure                           :: correct_polyn   => motion_patched_correct_polyn
    procedure                           :: kill            => motion_patched_kill
end type motion_patched

contains

    subroutine new( self, motion_correct_ftol, motion_correct_gtol, trs )
        class(motion_patched), intent(inout) :: self
        real, optional,        intent(in)    :: motion_correct_ftol, motion_correct_gtol
        real, optional,        intent(in)    :: trs
        call self%kill()
        self%motion_correct_ftol = TOL
        if (present(motion_correct_ftol)) self%motion_correct_ftol = motion_correct_ftol
        self%motion_correct_gtol = TOL
        if (present(motion_correct_gtol)) self%motion_correct_gtol = motion_correct_gtol
        self%trs = TRS_DEFAULT
        if (present(trs)) self%trs = trs
        allocate(self%lp(params_glob%nxpatch,params_glob%nypatch),&
                &self%updateres(params_glob%nxpatch,params_glob%nypatch),&
                &self%patch_centers(params_glob%nxpatch,params_glob%nypatch,2),&
                &self%lims_patches(params_glob%nxpatch,params_glob%nypatch,2,2))
        self%updateres   = 0
        self%lp          = -1.
        self%hp          = -1.
        self%resstep     = -1.
        self%fixed_frame = 1
        self%bfactor     = -1.
        self%motion_patched_direct_ftol = DIRECT_FTOL
        self%motion_patched_direct_gtol = DIRECT_GTOL
        self%existence = .true.
    end subroutine new

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

    pure function apply_patch_poly(c, x, y, t) result(res_sp)
        real(dp), intent(in) :: c(PATCH_PDIM), x, y, t
        real(sp) :: res_sp
        real(dp) :: res
        real(dp) :: x2, y2, xy, t2, t3
        x2 = x * x
        y2 = y * y
        xy = x * y
        t2 = t * t
        t3 = t2 * t
        res =       c( 1) * t      + c( 2) * t2      + c( 3) * t3
        res = res + c( 4) * t * x  + c( 5) * t2 * x  + c( 6) * t3 * x
        res = res + c( 7) * t * x2 + c( 8) * t2 * x2 + c( 9) * t3 * x2
        res = res + c(10) * t * y  + c(11) * t2 * y  + c(12) * t3 * y
        res = res + c(13) * t * y2 + c(14) * t2 * y2 + c(15) * t3 * y2
        res = res + c(16) * t * xy + c(17) * t2 * xy + c(18) * t3 * xy
        res_sp = real(res)
    end function apply_patch_poly

    subroutine fit_polynomial( self )
        class(motion_patched), intent(inout) :: self
        real(dp) :: yx(self%nframes*params_glob%nxpatch*params_glob%nypatch)      ! along x
        real(dp) :: yy(self%nframes*params_glob%nxpatch*params_glob%nypatch)      ! along y
        real(dp) :: x(3,self%nframes*params_glob%nxpatch*params_glob%nypatch)     ! x,y,t
        real(dp) :: sig(self%nframes*params_glob%nxpatch*params_glob%nypatch)
        real(dp) :: v(PATCH_PDIM,PATCH_PDIM), w(PATCH_PDIM), chisq
        integer  :: idx, iframe, i, j
        self%poly_chisq = 0.
        sig = 1.d0
        idx = 0
        do iframe = 1, self%nframes
            do i = 1, params_glob%nxpatch
                do j = 1, params_glob%nypatch
                    idx     = idx+1
                    yx(idx) = real(self%shifts_patches_for_fit(1,iframe,i,j),dp)
                    yy(idx) = real(self%shifts_patches_for_fit(2,iframe,i,j),dp)
                    call self%pix2polycoords(real(self%patch_centers(i,j,1),dp),&
                                            &real(self%patch_centers(i,j,2),dp),x(1,idx),x(2,idx))
                    x(3,idx) = real(iframe-self%fixed_frame,dp)
                end do
            end do
            if( self%frameweights(iframe) < 1.e-6 ) sig(iframe) = 1.d6
        end do
        call svd_multifit(x,yx,sig,self%poly_coeffs(:,1),v,w,chisq,patch_poly)
        self%poly_chisq(1) = real(chisq)
        call svd_multifit(x,yy,sig,self%poly_coeffs(:,2),v,w,chisq,patch_poly)
        self%poly_chisq(2) = real(chisq)
        self%poly_chisq    = sqrt( self%poly_chisq / real(count(self%frameweights>1.e-6)) ) ! average difference in pixels
    end subroutine fit_polynomial

    subroutine plot_shifts(self)
        class(motion_patched), intent(inout) :: self
        real, parameter       :: SCALE = 40.
        type(str4arr)         :: title
        type(CPlot2D_type)    :: plot2D
        type(CDataSet_type)   :: dataSetStart, dataSet, fit, obs, patch_start
        type(CDataPoint_type) :: point2, p_obs, p_fit, point
        integer               :: ipx,ipy, iframe, j
        real                  :: shifts(self%nframes,2), loc_shift(2), ref_shift(2)
        real                  :: xcenter,ycenter, cx,cy
        call CPlot2D__new(plot2D, self%shift_fname)
        call CPlot2D__SetXAxisSize(plot2D, 600._c_double)
        call CPlot2D__SetYAxisSize(plot2D, 600._c_double)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetFlipY(plot2D, C_TRUE)
        if (self%has_global_shifts) then
            ! centering to first frame for display only
            shifts      = self%global_shifts
            shifts(:,1) = shifts(:,1) - shifts(1,1)
            shifts(:,2) = shifts(:,2) - shifts(1,2)
            ! plot
            call CDataSet__new(dataSet)
            call CDataSet__SetDrawMarker(dataSet, C_FALSE)
            call CDataSet__SetDatasetColor(dataSet, 0.0_c_double, 0.0_c_double, 1.0_c_double)
            xcenter = real(self%ldim(1))/2.
            ycenter = real(self%ldim(2))/2.
            do j = 1, self%nframes
                call CDataPoint__new2( real(xcenter + SCALE * shifts(j, 1), c_double), &
                                      &real(ycenter + SCALE * shifts(j, 2), c_double), point)
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
        do ipx = 1, params_glob%nxpatch
            do ipy = 1, params_glob%nypatch
                cx = self%patch_centers(ipx, ipy, 1)
                cy = self%patch_centers(ipx, ipy, 2)
                ! centering to first frame for display only
                shifts(:,1) = self%shifts_patches_for_fit(1,:,ipx,ipy)
                shifts(:,2) = self%shifts_patches_for_fit(2,:,ipx,ipy)
                shifts(:,1) = shifts(:,1) - shifts(1,1)
                shifts(:,2) = shifts(:,2) - shifts(1,2)
                call self%get_local_shift(self%fixed_frame,cx,cy, ref_shift)
                ref_shift = ref_shift - shifts(self%fixed_frame,:)
                ! plot
                call CDataSet__new(patch_start)
                call CDataSet__SetDrawMarker(patch_start,C_TRUE)
                call CDataSet__SetMarkerSize(patch_start,5.0_c_double)
                call CDataSet__SetDatasetColor(patch_start,1.0_c_double,0.0_c_double,0.0_c_double)
                call CDataSet__new(fit)
                call CDataSet__new(obs)
                call CDataSet__SetDrawMarker(fit, C_FALSE)
                call CDataSet__SetDatasetColor(fit, 0.0_c_double,0.0_c_double,0.0_c_double)
                call CDataSet__SetDrawMarker(obs, C_FALSE)
                call CDataSet__SetDatasetColor(obs, 0.5_c_double,0.5_c_double,0.5_c_double)
                call CDataPoint__new2(real(cx,c_double), real(cy, c_double), p_fit)
                call CDataSet__AddDataPoint(patch_start, p_fit)
                call CDataPoint__delete(p_fit)
                do iframe = 1, self%nframes
                    if( self%frameweights(iframe) < 1.e-6 ) cycle
                    call CDataPoint__new2(real(cx + SCALE*shifts(iframe,1), c_double),&
                                         &real(cy + SCALE*shifts(iframe,2), c_double), p_obs)
                    call CDataSet__AddDataPoint(obs, p_obs)
                    call self%get_local_shift(iframe,cx,cy, loc_shift)
                    loc_shift = loc_shift - ref_shift
                    call CDataPoint__new2(real(cx + SCALE*loc_shift(1), c_double),&
                                         &real(cy + SCALE*loc_shift(2), c_double), p_fit)
                    call CDataSet__AddDataPoint(fit, p_fit)
                    call CDataPoint__delete(p_fit)
                    call CDataPoint__delete(p_obs)
                end do
                call CPlot2D__AddDataSet(plot2D, obs)
                call CPlot2D__AddDataSet(plot2D, fit)
                call CPlot2D__AddDataSet(plot2D, patch_start)
                call CDataSet__delete(patch_start)
                call CDataSet__delete(fit)
                call CDataSet__delete(obs)
            end do
        end do
        title%str = 'X (in pixels; trajectory scaled by '//trim(int2str(nint(SCALE)))//')'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%str)
        title%str(1:1) = 'Y'
        call CPlot2D__SetYAxisTitle(plot2D, title%str)
        call CPlot2D__OutputPostScriptPlot(plot2D, self%shift_fname)
        call CPlot2D__delete(plot2D)
    end subroutine plot_shifts

    ! produces shifts for 'polishing' close to relion 3 convention
    subroutine get_poly4star( self, polycoeffs )
        class(motion_patched), intent(inout) :: self
        real(dp), allocatable, intent(inout) :: polycoeffs(:)
        real(dp) :: yx(self%nframes*params_glob%nxpatch*params_glob%nypatch)      ! along x
        real(dp) :: yy(self%nframes*params_glob%nxpatch*params_glob%nypatch)      ! along y
        real(dp) :: x(3,self%nframes*params_glob%nxpatch*params_glob%nypatch)     ! x,y,t
        real(dp) :: sig(self%nframes*params_glob%nxpatch*params_glob%nypatch)
        real(dp) :: v(PATCH_PDIM,PATCH_PDIM), w(PATCH_PDIM), chisq
        integer  :: idx, iframe, i, j
        if( allocated(polycoeffs) ) deallocate(polycoeffs)
        allocate(polycoeffs(2*PATCH_PDIM),source=0.d0)
        idx = 0
        do iframe = 1, self%nframes
            do i = 1, params_glob%nxpatch
                do j = 1, params_glob%nypatch
                    idx     = idx+1
                    yx(idx) = real(self%shifts_patches_for_fit(1,iframe,i,j)-self%shifts_patches_for_fit(1,1,i,j),dp)
                    yy(idx) = real(self%shifts_patches_for_fit(2,iframe,i,j)-self%shifts_patches_for_fit(2,1,i,j),dp)
                    call self%pix2polycoords(real(self%patch_centers(i,j,1),dp),&
                                            &real(self%patch_centers(i,j,2),dp),x(1,idx),x(2,idx))
                    x(3,idx) = real(iframe-1,dp)
                end do
            end do
        end do
        sig = 1.d0
        call svd_multifit(x,yx,sig,polycoeffs(           1:  PATCH_PDIM),v,w,chisq,patch_poly)
        call svd_multifit(x,yy,sig,polycoeffs(PATCH_PDIM+1:2*PATCH_PDIM),v,w,chisq,patch_poly)
    end subroutine get_poly4star

    function get_polyfit_chisq( self )result( chisq )
        class(motion_patched), intent(in) :: self
        real :: chisq(2)
        chisq = self%poly_chisq
    end function get_polyfit_chisq

    elemental subroutine pix2polycoords( self, xin, yin, x, y )
        class(motion_patched), intent(in)  :: self
        real(dp),              intent(in)  :: xin, yin
        real(dp),              intent(out) :: x, y
        x = (xin-1.d0) / real(self%ldim(1)-1,dp) - 0.5d0
        y = (yin-1.d0) / real(self%ldim(2)-1,dp) - 0.5d0
    end subroutine pix2polycoords

    pure subroutine get_local_shift( self, iframe, x, y, shift )
        class(motion_patched), intent(inout) :: self
        integer,               intent(in)  :: iframe
        real,                  intent(in)  :: x, y
        real,                  intent(out) :: shift(2)
        real(dp) :: t, xx, yy
        t  = real(iframe-self%fixed_frame, dp)
        call self%pix2polycoords(real(x,dp),real(y,dp), xx,yy)
        shift(1) = apply_patch_poly(self%poly_coeffs(:,1), xx,yy,t)
        shift(2) = apply_patch_poly(self%poly_coeffs(:,2), xx,yy,t)
    end subroutine get_local_shift

    !>  Per frame real space polynomial interpolation
    subroutine polytransfo( self, iframe, frame, frame_output )
        class(motion_patched), intent(inout) :: self
        integer,               intent(in)    :: iframe
        type(image),           intent(inout) :: frame, frame_output
        integer  :: i, j
        real     :: coords(2), sh(2), shinterp(2)
        type(rmat_ptr_type) :: rmat_in, rmat_out
        call frame_output%new(self%ldim, self%smpd, wthreads=.false.)
        call frame%ifft()
        call frame%get_rmat_ptr(rmat_in%rmat_ptr)
        call frame_output%get_rmat_ptr(rmat_out%rmat_ptr)
        do i = 1, self%ldim(1)
            do j = 1, self%ldim(2)
                call self%get_local_shift(iframe, real(i), real(j), sh)
                call self%get_local_shift(self%interp_fixed_frame, real(i), real(j), shinterp)
                coords = real([i,j]) - sh + shinterp
                rmat_out%rmat_ptr(i,j,1) = interp_bilin(coords(1), coords(2))
            end do
        end do

    contains

        pure real function interp_bilin( xval, yval )
            real, intent(in) :: xval, yval
            integer  :: x1_h,  x2_h,  y1_h,  y2_h
            real     :: y1, y2, y3, y4, t, u
            logical  :: outside
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
                interp_bilin = rmat_in%rmat_ptr(x1_h, y1_h, 1)
                return
            endif
            y1 = rmat_in%rmat_ptr(x1_h, y1_h, 1)
            y2 = rmat_in%rmat_ptr(x2_h, y1_h, 1)
            y3 = rmat_in%rmat_ptr(x2_h, y2_h, 1)
            y4 = rmat_in%rmat_ptr(x1_h, y2_h, 1)
            t   = xval - real(x1_h)
            u   = yval - real(y1_h)
            interp_bilin =  (1. - t) * (1. - u) * y1 + &
                                 &t  * (1. - u) * y2 + &
                                 &t  *       u  * y3 + &
                           &(1. - t) *       u  * y4
        end function interp_bilin

    end subroutine polytransfo

    subroutine apply_polytransfo( self, frames, frames_output )
        class(motion_patched),    intent(inout) :: self
        type(image), allocatable, intent(inout) :: frames(:)
        type(image), allocatable, intent(inout) :: frames_output(:)
        integer  :: i, j, iframe
        real     :: x, y, sh(2), shinterp(2)
        type(rmat_ptr_type) :: rmat_ins(self%nframes), rmat_outs(self%nframes)
        !$omp parallel default(shared) private(iframe,j,i,sh,shinterp,x,y) proc_bind(close)
        !$omp do schedule(static)
        do iframe = 1, self%nframes
            call frames_output(iframe)%new(self%ldim, self%smpd, wthreads=.false.)
            call frames(iframe)%ifft()
            call frames(iframe)%get_rmat_ptr(rmat_ins(iframe)%rmat_ptr)
            call frames_output(iframe)%get_rmat_ptr(rmat_outs(iframe)%rmat_ptr)
        end do
        !$omp end do
        !$omp do collapse(3) schedule(static)
        do iframe = 1, self%nframes
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    call self%get_local_shift(iframe, real(i), real(j), sh)
                    call self%get_local_shift(self%interp_fixed_frame, real(i), real(j), shinterp)
                    x = real(i) - sh(1) + shinterp(1)
                    y = real(j) - sh(2) + shinterp(2)
                    rmat_outs(iframe)%rmat_ptr(i,j,1) = interp_bilin(x, y, iframe)
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
    contains

        pure real function interp_bilin( xval, yval, iiframe )
            real,                intent(in)  :: xval, yval
            integer,             intent(in)  :: iiframe
            integer  :: x1_h,  x2_h,  y1_h,  y2_h
            real     :: y1, y2, y3, y4, t, u
            logical  :: outside
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
                interp_bilin = rmat_ins(iiframe)%rmat_ptr(x1_h, y1_h, 1)
                return
            endif
            y1 = rmat_ins(iiframe)%rmat_ptr(x1_h, y1_h, 1)
            y2 = rmat_ins(iiframe)%rmat_ptr(x2_h, y1_h, 1)
            y3 = rmat_ins(iiframe)%rmat_ptr(x2_h, y2_h, 1)
            y4 = rmat_ins(iiframe)%rmat_ptr(x1_h, y2_h, 1)
            t   = xval - real(x1_h)
            u   = yval - real(y1_h)
            interp_bilin =  (1. - t) * (1. - u) * y1 + &
                        &t  * (1. - u) * y2 + &
                        &t  *       u  * y3 + &
                        &(1. - t) * u  * y4
        end function interp_bilin

    end subroutine apply_polytransfo

    subroutine allocate_fields( self )
        class(motion_patched), intent(inout) :: self
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
            allocate(self%shifts_patches   (2, self%nframes, params_glob%nxpatch, params_glob%nypatch),&
                self%shifts_patches_for_fit(2, self%nframes, params_glob%nxpatch, params_glob%nypatch),&
                self%frame_patches(params_glob%nxpatch, params_glob%nypatch), stat=alloc_stat )
            if (alloc_stat /= 0) call allocchk('allocate_fields 1; simple_motion_patched')
            self%shifts_patches         = 0.
            self%shifts_patches_for_fit = 0.
            do i = 1, params_glob%nxpatch
                do j = 1, params_glob%nypatch
                    allocate( self%frame_patches(i, j)%stack(self%nframes), stat=alloc_stat )
                    if (alloc_stat /= 0) call allocchk('allocate_fields 2; simple_motion_patched')
                end do
            end do
        end if
    end subroutine allocate_fields

    subroutine deallocate_fields( self )
        class(motion_patched), intent(inout) :: self
        integer :: iframe, i, j
        if (allocated(self%shifts_patches)) deallocate(self%shifts_patches)
        if (allocated(self%shifts_patches_for_fit)) deallocate(self%shifts_patches_for_fit)
        if (allocated(self%frame_patches)) then
            do j = 1, params_glob%nypatch
                do i = 1, params_glob%nxpatch
                    do iframe = 1, self%nframes
                        call self%frame_patches(i,j)%stack(iframe)%kill()
                    end do
                    deallocate(self%frame_patches(i,j)%stack)
                end do
            end do
            deallocate(self%frame_patches)
        end if
    end subroutine deallocate_fields

    subroutine set_size_frames_ref( self )
        class(motion_patched), intent(inout) :: self
        integer :: i,j
        real    :: cen, dist
        self%ldim_patch(1) = round2even(real(self%ldim(1)) / real(params_glob%nxpatch))
        self%ldim_patch(2) = round2even(real(self%ldim(2)) / real(params_glob%nypatch))
        self%ldim_patch(3) = 1
        ! along X
        ! limits & center first patches
        self%lims_patches(1,:,1,1) = 1
        self%lims_patches(1,:,1,2) = self%ldim_patch(1)
        self%patch_centers(1,:,1)  = sum(self%lims_patches(1,:,1,1:2),dim=2) / 2.
        ! limits & center last patches
        self%lims_patches(params_glob%nxpatch,:,1,1) = self%ldim(1)-self%ldim_patch(1)+1
        self%lims_patches(params_glob%nxpatch,:,1,2) = self%ldim(1)
        self%patch_centers(params_glob%nxpatch,:,1)  = sum(self%lims_patches(params_glob%nxpatch,:,1,1:2),dim=2) / 2.
        ! adjust other patch centers to be evenly spread
        dist = real(self%patch_centers(params_glob%nxpatch,1,1)-self%patch_centers(1,1,1)+1) / real(params_glob%nxpatch-1)
        do i=2,params_glob%nxpatch-1
            cen = self%patch_centers(1,1,1) + real(i-1)*dist
            self%lims_patches(i,:,1,1) = ceiling(cen) - self%ldim_patch(1)/2
            self%lims_patches(i,:,1,2) = self%lims_patches(i,:,1,1) + self%ldim_patch(1) - 1
            self%patch_centers(i,:,1)  = sum(self%lims_patches(i,:,1,1:2),dim=2) / 2.
        enddo
        ! along Y
        self%lims_patches(:,1,2,1) = 1
        self%lims_patches(:,1,2,2) = self%ldim_patch(2)
        self%patch_centers(:,1,2)  = sum(self%lims_patches(:,1,2,1:2),dim=2) / 2.
        self%lims_patches(:,params_glob%nypatch,2,1) = self%ldim(2)-self%ldim_patch(2)+1
        self%lims_patches(:,params_glob%nypatch,2,2) = self%ldim(2)
        self%patch_centers(:,params_glob%nypatch,2)  = sum(self%lims_patches(:,params_glob%nypatch,2,1:2),dim=2) / 2.
        dist = real(self%patch_centers(1,params_glob%nypatch,2)-self%patch_centers(1,1,2)+1) / real(params_glob%nypatch-1)
        do j=2,params_glob%nypatch-1
            cen = self%patch_centers(1,1,2) + real(j-1)*dist
            self%lims_patches(:,j,2,1) = ceiling(cen) - self%ldim_patch(2)/2
            self%lims_patches(:,j,2,2) = self%lims_patches(:,j,2,1) + self%ldim_patch(2) - 1
            self%patch_centers(:,j,2)  = sum(self%lims_patches(:,j,2,1:2),dim=2) /2.
        enddo
    end subroutine set_size_frames_ref

    subroutine set_patches( self, stack )
        class(motion_patched),          intent(inout) :: self
        type(image),       allocatable, intent(inout) :: stack(:)
        real, allocatable   :: res(:)
        integer             :: i, j, iframe, k, l, kk, ll
        integer             :: ip, jp           ! ip, jp: i_patch, j_patch
        integer             :: lims_patch(2,2)
        type(rmat_ptr_type) :: rmat_ptrs(self%nframes)
        real, pointer       :: rmat_patch(:,:,:)
        ! init
        do j = 1, params_glob%nypatch
            do i = 1, params_glob%nxpatch
                do iframe=1,self%nframes
                    call self%frame_patches(i,j)%stack(iframe)%new(self%ldim_patch, self%smpd)
                end do
            end do
        end do
        ! initialize transfer matrix to correct dimensions
        call ftexp_transfmat_init(self%frame_patches(1,1)%stack(1))
        ! fill patches
        do iframe=1,self%nframes
            call stack(iframe)%get_rmat_ptr(rmat_ptrs(iframe)%rmat_ptr)
        end do
        !$omp parallel do collapse(3) default(shared) private(iframe,j,i,lims_patch,rmat_patch,k,l,kk,ll,ip,jp) proc_bind(close) schedule(static)
        do iframe=1,self%nframes
            do j = 1, params_glob%nypatch
                do i = 1, params_glob%nxpatch
                    call self%frame_patches(i,j)%stack(iframe)%get_rmat_ptr(rmat_patch)
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
        !$omp end parallel do
        ! updates high-pass according to new dimensions
        res = self%frame_patches(1,1)%stack(1)%get_res()
        self%hp = min(self%hp,res(1))
        deallocate(res)
        write(logfhandle,'(A,F6.1)')'>>> PATCH HIGH-PASS: ',self%hp
    end subroutine set_patches

    subroutine gen_patch( self, stack, pi, pj )
        class(motion_patched),          intent(inout) :: self
        type(image),       allocatable, intent(inout) :: stack(:)
        integer,                        intent(in)    :: pi, pj
        integer             :: iframe, k, l, kk, ll
        integer             :: ip, jp           ! ip, jp: i_patch, j_patch
        real, pointer       :: prmat_patch(:,:,:), prmat_frame(:,:,:)
        do iframe=1,self%nframes
            ! init
            call self%frame_patches(pi,pj)%stack(iframe)%new(self%ldim_patch, self%smpd)
            call self%frame_patches(pi,pj)%stack(iframe)%get_rmat_ptr(prmat_patch)
            call stack(iframe)%get_rmat_ptr(prmat_frame)
            ! copy
            do k = self%lims_patches(pi,pj,1,1), self%lims_patches(pi,pj,1,2)
                kk = k
                if (kk < 1) then
                    kk = kk + self%ldim(1)
                else if (kk > self%ldim(1)) then
                    kk = kk - self%ldim(1)
                end if
                ip = k - self%lims_patches(pi,pj,1,1) + 1
                do l = self%lims_patches(pi,pj,2,1), self%lims_patches(pi,pj,2,2)
                    ll = l
                    if (ll < 1) then
                        ll = ll + self%ldim(2)
                    else if (ll > self%ldim(2)) then
                        ll = ll - self%ldim(2)
                    end if
                    jp = l - self%lims_patches(pi,pj,2,1) + 1
                    ! now copy the value
                    prmat_patch(ip,jp,1) = prmat_frame(kk,ll,1)
                end do
            end do
        end do
    end subroutine gen_patch

    subroutine det_shifts( self, frames )
        class(motion_patched), target, intent(inout) :: self
        type(image),      allocatable, intent(inout) :: frames(:)
        type(motion_align_hybrid), allocatable :: align_hybrid(:,:)
        real, allocatable :: opt_shifts(:,:), res(:)
        real              :: corr_avg
        integer           :: iframe, i, j, alloc_stat
        self%shifts_patches = 0.
        allocate(align_hybrid(params_glob%nxpatch, params_glob%nypatch), stat=alloc_stat )
        if (alloc_stat /= 0) call allocchk('det_shifts 1; simple_motion_patched')
        corr_avg = 0.
        ! initialize transfer matrix to correct dimensions
        call self%frame_patches(1,1)%stack(1)%new(self%ldim_patch, self%smpd)
        call ftexp_transfmat_init(self%frame_patches(1,1)%stack(1))
        res = self%frame_patches(1,1)%stack(1)%get_res()
        self%hp = min(self%hp,res(1))
        write(logfhandle,'(A,F6.1)')'>>> PATCH HIGH-PASS: ',self%hp
        !$omp parallel do collapse(2) default(shared) private(i,j,iframe,opt_shifts)&
        !$omp proc_bind(close) schedule(static) reduction(+:corr_avg)
        do i = 1,params_glob%nxpatch
            do j = 1,params_glob%nypatch
                ! init
                self%lp(i,j) = (params_glob%lpstart+params_glob%lpstop)/2.
                call self%gen_patch(frames,i,j)
                call align_hybrid(i,j)%new(self%frame_patches(i,j)%stack)
                call align_hybrid(i,j)%set_group_frames(trim(params_glob%groupframes).eq.'yes')
                call align_hybrid(i,j)%set_rand_init_shifts(.true.)
                call align_hybrid(i,j)%set_reslims(self%hp, self%lp(i,j), params_glob%lpstop)
                call align_hybrid(i,j)%set_trs(params_glob%scale*params_glob%trs)
                call align_hybrid(i,j)%set_shsrch_tol(TOL)
                call align_hybrid(i,j)%set_coords(i,j)
                call align_hybrid(i,j)%set_fitshifts(self%fitshifts)
                call align_hybrid(i,j)%set_fixed_frame(self%fixed_frame)
                call align_hybrid(i,j)%set_bfactor(self%bfactor)
                call align_hybrid(i,j)%set_downscale(.not.(trim(params_glob%groupframes).eq.'yes'))
                ! align
                call align_hybrid(i,j)%align(frameweights=self%frameweights)
                ! fetch info
                corr_avg = corr_avg + align_hybrid(i,j)%get_corr()
                call align_hybrid(i,j)%get_opt_shifts(opt_shifts)
                ! making sure the shifts are in reference to fixed_frame
                do iframe = 1, self%nframes
                    self%shifts_patches(:,iframe,i,j) = opt_shifts(iframe,:) - opt_shifts(self%fixed_frame,:)
                end do
                ! cleanup
                call align_hybrid(i,j)%kill
                do iframe=1,self%nframes
                    call self%frame_patches(i,j)%stack(iframe)%kill
                end do
            end do
        end do
        !$omp end parallel do
        self%shifts_patches_for_fit = self%shifts_patches
        corr_avg = corr_avg / real(params_glob%nxpatch*params_glob%nypatch)
        write(logfhandle,'(A,F6.3)')'>>> AVERAGE PATCH & FRAMES CORRELATION: ', corr_avg
        deallocate(align_hybrid,res)
    end subroutine det_shifts

    subroutine det_shifts_polyn( self )
        class(motion_patched), target, intent(inout) :: self
        real, allocatable :: opt_shifts(:,:)
        real              :: corr_avg
        integer           :: iframe, i, j, alloc_stat
        self%shifts_patches = 0.
        allocate( self%align_iso_polyn_direct(params_glob%nxpatch, params_glob%nypatch), stat=alloc_stat )
        if (alloc_stat /= 0) call allocchk('det_shifts 1; simple_motion_patched_polyn')
        !$omp parallel do collapse(2) default(shared) private(j,i) proc_bind(close) schedule(static)
        do i = 1, params_glob%nxpatch
            do j = 1, params_glob%nypatch
                call self%align_iso_polyn_direct(i,j)%new
                call self%align_iso_polyn_direct(i,j)%set_frames(self%frame_patches(i,j)%stack, self%nframes)
                call self%align_iso_polyn_direct(i,j)%set_hp_lp(self%hp, self%lp(i,j))
                call self%align_iso_polyn_direct(i,j)%set_trs(self%trs)
                call self%align_iso_polyn_direct(i,j)%set_ftol_gtol(TOL, TOL)
                call self%align_iso_polyn_direct(i,j)%set_coords(i,j)
                call self%align_iso_polyn_direct(i,j)%set_callback(motion_patched_polyn_callback_wrapper)
                call self%align_iso_polyn_direct(i,j)%align_polyn(self)
            end do
        end do
        !$omp end parallel do
        corr_avg = 0.
        do i = 1, params_glob%nxpatch
            do j = 1, params_glob%nypatch
                corr_avg = corr_avg + self%align_iso_polyn_direct(i,j)%get_corr()
            enddo
        enddo
        corr_avg = corr_avg / real(params_glob%nxpatch*params_glob%nypatch)
        write(logfhandle,'(A,F6.3)')'>>> AVERAGE PATCH & FRAMES CORRELATION: ', corr_avg
        ! Set the first shift to 0.
        do i = 1, params_glob%nxpatch
            do j = 1, params_glob%nypatch
                call self%align_iso_polyn_direct(i,j)%get_opt_shifts(opt_shifts)
                do iframe = 1, self%nframes
                    self%shifts_patches(:,iframe,i,j) = opt_shifts(iframe,:)
                end do
            end do
        end do
        do iframe = self%nframes, 1, -1
            self%shifts_patches(:,iframe,:,:) = self%shifts_patches(:,iframe,:,:)-self%shifts_patches(:,1,:,:)
        enddo
        do iframe = 1, self%nframes
            self%shifts_patches_for_fit(1,iframe,:,:) = self%shifts_patches(1,iframe,:,:) + 0.5*self%shifts_patches(1,1,:,:)
            self%shifts_patches_for_fit(2,iframe,:,:) = self%shifts_patches(2,iframe,:,:) + 0.5*self%shifts_patches(2,1,:,:)
        enddo
        ! no cleanup yet
    end subroutine det_shifts_polyn

    subroutine det_shifts_direct( self )
        class(motion_patched), target, intent(inout) :: self
        type(opt_factory)         :: ofac
        type(opt_spec)            :: ospec
        class(optimizer), pointer :: nlopt
        real                      :: opt_lims(PATCH_PDIM*2, 2), lowest_cost
        opt_lims(:,1) = -self%trs
        opt_lims(:,2) =  self%trs
        call ospec%specify('lbfgsb', PATCH_PDIM*2, ftol=self%motion_patched_direct_ftol, &
            gtol=self%motion_patched_direct_gtol, limits=opt_lims, maxits=800)
        call ospec%set_costfun_8(patched_direct_cost_wrapper)
        call ospec%set_gcostfun_8(patched_direct_gcost_wrapper)
        call ospec%set_fdfcostfun_8(patched_direct_fdf_wrapper)
        call ofac%new(ospec, nlopt)
        ospec%x_8(           1:PATCH_PDIM  ) = self%poly_coeffs(:,1)
        ospec%x_8(PATCH_PDIM+1:PATCH_PDIM*2) = self%poly_coeffs(:,2)
        ospec%x = real(ospec%x_8)
        call nlopt%minimize(ospec, self, lowest_cost)
        self%poly_coeffs(:,1) = ospec%x_8(           1:PATCH_PDIM  )
        self%poly_coeffs(:,2) = ospec%x_8(PATCH_PDIM+1:PATCH_PDIM*2)
        nlopt => null()
    end subroutine det_shifts_direct

    subroutine cleanup_polyn( self )
        class(motion_patched), intent(inout) :: self
        integer :: i, j
        do i = 1, params_glob%nxpatch
            do j = 1, params_glob%nypatch
                call self%align_iso_polyn_direct(i,j)%kill
            end do
        end do
        deallocate(self%align_iso_polyn_direct)
    end subroutine cleanup_polyn

    subroutine grad_contrib_to_full_grad( self, grad_contrib, x, y, grad_full )
        class(motion_patched), intent(inout) :: self
        real(dp),                     intent(in)    :: grad_contrib(POLYDIM2)
        real(dp),                     intent(in)    :: x, y
        real(dp),                     intent(out)   :: grad_full(PATCH_PDIM*2)
        real(dp) :: x2, y2, xy
        ! X:    (c01+c04 x+c07 x^2+c10 y+c13 y^2+c16 xy)*t   + (c02+c05 x+c08 x^2+c11 y+c14 y^2+c17 xy)*t^2
        !     + (c03+c06 x+c09 x^2+c12 y+c15 y^2+c18 xy)*t^3
        ! Y:    (c19+c22 x+c25 x^2+c28 y+c31 y^2+c34 xy)*t   + (c20+c23 x+c26 x^2+c29 y+c32 y^2+c35 xy)*t^2
        !     + (c21+c24 x+c27 x^2+c30 y+c33 y^2+c36 xy)*t^3
        x2 = x**2
        y2 = y**2
        xy = x*y
        grad_full( 1) =      grad_contrib(1)
        grad_full( 4) = x  * grad_contrib(1)
        grad_full( 7) = x2 * grad_contrib(1)
        grad_full(10) = y  * grad_contrib(1)
        grad_full(13) = y2 * grad_contrib(1)
        grad_full(16) = xy * grad_contrib(1)
        grad_full( 2) =      grad_contrib(2)
        grad_full( 5) = x  * grad_contrib(2)
        grad_full( 8) = x2 * grad_contrib(2)
        grad_full(11) = y  * grad_contrib(2)
        grad_full(14) = y2 * grad_contrib(2)
        grad_full(17) = xy * grad_contrib(2)
        grad_full( 3) =      grad_contrib(3)
        grad_full( 6) = x  * grad_contrib(3)
        grad_full( 9) = x2 * grad_contrib(3)
        grad_full(12) = y  * grad_contrib(3)
        grad_full(15) = y2 * grad_contrib(3)
        grad_full(18) = xy * grad_contrib(3)
        grad_full(19) =      grad_contrib(4)
        grad_full(22) = x  * grad_contrib(4)
        grad_full(25) = x2 * grad_contrib(4)
        grad_full(28) = y  * grad_contrib(4)
        grad_full(31) = y2 * grad_contrib(4)
        grad_full(34) = xy * grad_contrib(4)
        grad_full(20) =      grad_contrib(5)
        grad_full(23) = x  * grad_contrib(5)
        grad_full(26) = x2 * grad_contrib(5)
        grad_full(29) = y  * grad_contrib(5)
        grad_full(32) = y2 * grad_contrib(5)
        grad_full(35) = xy * grad_contrib(5)
        grad_full(21) =      grad_contrib(6)
        grad_full(24) = x  * grad_contrib(6)
        grad_full(27) = x2 * grad_contrib(6)
        grad_full(30) = y  * grad_contrib(6)
        grad_full(33) = y2 * grad_contrib(6)
        grad_full(36) = xy * grad_contrib(6)
    end subroutine grad_contrib_to_full_grad

    subroutine set_frameweights( self, frameweights )
        class(motion_patched), intent(inout) :: self
        real, allocatable,     intent(in) :: frameweights(:)
        integer :: nlen
        nlen = size(frameweights)
        if (allocated(self%frameweights)) deallocate(self%frameweights)
        allocate(self%frameweights(nlen), source=frameweights)
        self%has_frameweights = .true.
    end subroutine set_frameweights

    subroutine set_fitshifts( self, fitshifts )
        class(motion_patched), intent(inout) :: self
        logical,               intent(in)    :: fitshifts
        self%fitshifts = fitshifts
    end subroutine set_fitshifts

    subroutine set_fixed_frame( self, fixed_frame )
        class(motion_patched), intent(inout) :: self
        integer,               intent(in)    :: fixed_frame
        self%fixed_frame = fixed_frame
    end subroutine set_fixed_frame

    subroutine set_interp_fixed_frame( self, fixed_frame )
        class(motion_patched), intent(inout) :: self
        integer,               intent(in)    :: fixed_frame
        self%interp_fixed_frame = fixed_frame
    end subroutine set_interp_fixed_frame

    subroutine set_bfactor( self, bfac )
        class(motion_patched), intent(inout) :: self
        real,                  intent(in)    :: bfac
        self%bfactor = bfac
    end subroutine set_bfactor

    ! EXECUTION ROUTINES

    subroutine motion_patched_correct( self, hp, resstep, frames, shift_fname, global_shifts )
        class(motion_patched),           intent(inout) :: self
        real,                            intent(in)    :: hp, resstep
        type(image),        allocatable, intent(inout) :: frames(:)
        character(len=:),   allocatable, intent(in)    :: shift_fname
        real,     optional, allocatable, intent(in)    :: global_shifts(:,:)
        integer :: ldim_frames(3)
        integer :: i
        ! prep
        self%hp          = hp
        self%lp          = params_glob%lpstart
        self%resstep     = resstep
        self%updateres   = 0
        self%shift_fname = shift_fname // C_NULL_CHAR
        if (allocated(self%global_shifts)) deallocate(self%global_shifts)
        self%has_global_shifts = .false.
        if (present(global_shifts)) then
            allocate(self%global_shifts(size(global_shifts, 1), size(global_shifts, 2)))
            self%global_shifts     = global_shifts
            self%has_global_shifts = .true.
        end if
        self%nframes = size(frames,dim=1)
        self%ldim    = frames(1)%get_ldim()
        self%smpd    = frames(1)%get_smpd()
        do i = 1,self%nframes
            ldim_frames = frames(i)%get_ldim()
            if (any(ldim_frames(1:2) /= self%ldim(1:2))) then
                THROW_HARD('error in motion_patched_correct: frame dimensions do not match reference dimension; simple_motion_patched')
            end if
        end do
        call self%allocate_fields()
        ! determines patch geometry
        call self%set_size_frames_ref()
        ! determine shifts for patches
        call self%det_shifts(frames)
        ! fit the polynomial model against determined shifts
        call self%fit_polynomial()
        ! report visual results
        call self%plot_shifts()
    end subroutine motion_patched_correct

    subroutine motion_patched_correct_polyn( self, hp, resstep, frames, shift_fname, refine_direct, global_shifts )
        class(motion_patched),           intent(inout) :: self
        real,                            intent(in)    :: hp, resstep
        type(image),        allocatable, intent(inout) :: frames(:)
        character(len=:),   allocatable, intent(in)    :: shift_fname
        logical,                         intent(in)    :: refine_direct
        real,     optional, allocatable, intent(in)    :: global_shifts(:,:)
        integer :: ldim_frames(3)
        integer :: i
        write (*,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ patched polyn ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
        ! prep
        self%hp          = hp
        self%lp          = params_glob%lpstart
        self%resstep     = resstep
        self%updateres   = 0
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
                THROW_HARD('error in motion_patched_correct: frame dimensions do not match reference dimension; simple_motion_patched')
            end if
        end do
        call self%allocate_fields()
        call self%set_size_frames_ref()
        ! divide the reference into patches & updates high-pass accordingly
        call self%set_patches(frames)
        ! determine shifts for patches
        call self%det_shifts_polyn()
        ! fit the polynomial model against determined shifts
        call self%fit_polynomial()
        write (*,*) '^^^^^^^^^^^^^^^^^^^^^^^ fitted polynomial; poly_coeffs(:, 1)=', self%poly_coeffs(:,1), 'poly_coeffs(:, 2)=', self%poly_coeffs(:,2)
        if ( refine_direct ) then
            write (*,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ patched direct ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ '
            call self%det_shifts_direct()
        end if
        ! clean up align_iso_polyn_direct objects
        call self%cleanup_polyn()
        ! report visual results
        call self%plot_shifts()
    end subroutine motion_patched_correct_polyn

    subroutine motion_patched_kill( self )
        class(motion_patched), intent(inout) :: self
        call self%deallocate_fields()
        if (allocated(self%frameweights)) deallocate(self%frameweights)
        if (allocated(self%updateres)) deallocate(self%updateres)
        if (allocated(self%lp)) deallocate(self%lp)
        if (allocated(self%patch_centers)) deallocate(self%patch_centers)
        if (allocated(self%lims_patches)) deallocate(self%lims_patches)
        self%has_frameweights = .false.
        self%existence        = .false.
        call ftexp_transfmat_kill
    end subroutine motion_patched_kill

    subroutine motion_patched_polyn_callback(self, align_iso_polyn_direct, converged)
        class(motion_patched),                intent(inout) :: self
        class(motion_align_iso_polyn_direct), intent(inout) :: align_iso_polyn_direct
        logical,                              intent(out)   :: converged
        integer :: i, j
        logical :: didupdateres
        didupdateres = .false.
        call align_iso_polyn_direct%get_coords(i, j)
        select case(self%updateres(i,j))
        case(0)
            call update_res( self%updateres(i,j) )
        case(1)
            call update_res( self%updateres(i,j) )
        case(2)
            call update_res( self%updateres(i,j) )
            call align_iso_polyn_direct%set_factr_pgtol(1d+6, 1d-6)
        case DEFAULT
            ! nothing to do
        end select
        if( self%updateres(i,j) > 2 .and. .not. didupdateres)then ! at least one iteration with new lim
            call align_iso_polyn_direct%reset_tols
            converged = .true.
        else
            converged = .false.
        end if

    contains

        subroutine update_res( which_update )
            integer, intent(in) :: which_update
            self%lp(i,j) = self%lp(i,j) - self%resstep
            call align_iso_polyn_direct%set_hp_lp(self%hp, self%lp(i,j))
            write(logfhandle,'(a,1x,f7.4)') '>>> LOW-PASS LIMIT UPDATED TO:', self%lp(i,j)
            ! need to indicate that we updated resolution limit
            self%updateres(i,j)  = self%updateres(i,j) + 1
            ! indicate that reslim was updated
            didupdateres = .true.
        end subroutine update_res

    end subroutine motion_patched_polyn_callback

    subroutine motion_patched_polyn_callback_wrapper(aptr, align_iso_polyn_direct, converged)
        class(*),                             intent(inout) :: aptr
        class(motion_align_iso_polyn_direct), intent(inout) :: align_iso_polyn_direct
        logical,                              intent(out)   :: converged
        select type(aptr)
        class is (motion_patched)
            call aptr%motion_patched_polyn_callback(align_iso_polyn_direct, converged)
        class default
            THROW_HARD('error in motion_patched_polyn_callback_wrapper: unknown type; simple_motion_patched_polyn')
        end select
    end subroutine motion_patched_polyn_callback_wrapper

    subroutine coeffs_to_shifts( self )
        class(motion_patched), intent(inout) :: self
        integer  :: xi, yi, j
        real(dp) :: ashift(2)
        real(dp) :: x, y, t
        do xi = 1, params_glob%nxpatch
            do yi = 1, params_glob%nypatch
                call self%pix2polycoords(real(self%patch_centers(xi,yi,1),dp), real(self%patch_centers(xi,yi,2),dp), x, y) ! TODO: memoize
                do j = 1, self%nframes
                    t = real(j - 1, dp)
                    ashift(1) = apply_patch_poly(self%poly_coeffs(:,1), x, &
                        y, t)
                    ashift(2) = apply_patch_poly(self%poly_coeffs(:,2), x, &
                        y, t)
                    self%align_iso_polyn_direct(xi,yi)%shifts(j,:) = - ashift(:)
                end do
            end do
        end do
    end subroutine coeffs_to_shifts

    function motion_patched_direct_cost( self, vec ) result( r )
        class(motion_patched), intent(inout) :: self
        real(dp),              intent(in)    :: vec(PATCH_PDIM*2)
        real(dp) :: r
        integer  :: xi, yi, j
        real(dp) :: x, y
        real(dp) :: val, t
        real(dp) :: ashift(2)
        r = 0.d0
        ! convert polynomial coefficients into shifts
        !$omp parallel do collapse(2) default(shared) private(x,y) proc_bind(close) schedule(static)
        do xi = 1, params_glob%nxpatch
            do yi = 1, params_glob%nypatch
                call self%pix2polycoords(real(self%patch_centers(xi,yi,1),dp), real(self%patch_centers(xi,yi,2),dp), x, y) ! TODO: memoize
                do j = 1, self%nframes
                    t = real(j - 1, dp)
                    ashift(1) = apply_patch_poly(vec(           1:PATCH_PDIM  ), x, &
                        y, t)
                    ashift(2) = apply_patch_poly(vec(PATCH_PDIM+1:PATCH_PDIM*2), x, &
                        y, t)
                    self%align_iso_polyn_direct(xi,yi)%shifts(j,:) = - ashift(:)
                end do
                r = r + self%align_iso_polyn_direct(xi,yi)%motion_align_iso_contribs_cost()
            end do
        end do
        !$omp end parallel do
    end function motion_patched_direct_cost

    function patched_direct_cost_wrapper( self, vec, D ) result( cost )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(kind=8), intent(in)    :: vec(D)
        real(kind=8) :: cost
        select type(self)
            class is (motion_patched)
                cost = self%motion_patched_direct_cost( vec )
            class DEFAULT
                THROW_HARD('unknown type; patched_direct_cost_wrapper')
        end select
    end function patched_direct_cost_wrapper

    subroutine motion_patched_direct_gcost( self, vec, grad )
        class(motion_patched), intent(inout) :: self
        real(dp),                     intent(in)    :: vec(PATCH_PDIM*2)
        real(dp),                     intent(out)   :: grad(PATCH_PDIM*2)
        real(dp) :: r
        integer  :: xi, yi, j
        real(dp) :: x, y
        real(dp) :: t
        real(dp) :: ashift(2)
        real(dp) :: grad_contrib(POLYDIM2), grad_full_tmp(PATCH_PDIM*2)
        grad(:) = 0.d0
        ! convert polynomial coefficients into shifts
        do xi = 1, params_glob%nxpatch
            do yi = 1, params_glob%nypatch
                do j = 1, self%nframes
                    t = real(j - 1, dp)
                    call self%pix2polycoords(real(self%patch_centers(xi,yi,1),dp), real(self%patch_centers(xi,yi,2),dp), x, y) ! TODO: memoize
                    ashift(1) = apply_patch_poly(vec(           1:PATCH_PDIM  ), x, y, t)
                    ashift(2) = apply_patch_poly(vec(PATCH_PDIM+1:PATCH_PDIM*2), x, y, t)
                    self%align_iso_polyn_direct(xi,yi)%shifts(j,:) = - ashift(:)
                end do
                ! calculate gradient contribution (partial)
                call self%align_iso_polyn_direct(xi,yi)%motion_align_iso_contribs_gcost(grad_contrib)
                call self%grad_contrib_to_full_grad( grad_contrib, x, y, grad_full_tmp )
                grad(:) = grad(:) - grad_full_tmp(:)
            end do
        end do
    end subroutine motion_patched_direct_gcost

    subroutine patched_direct_gcost_wrapper( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        grad = 0.d0
        select type(self)
            class is (motion_patched)
                call self%motion_patched_direct_gcost( vec, grad )
            class DEFAULT
                THROW_HARD('unknown type; patched_direct_gcost_wrapper')
        end select
    end subroutine patched_direct_gcost_wrapper

    !< cost function for minimizer, f and gradient
    subroutine motion_patched_direct_fdf( self, vec, f, grad )
        class(motion_patched), intent(inout) :: self
        real(dp),                     intent(in)    :: vec (PATCH_PDIM*2)
        real(dp),                     intent(out)   :: grad(PATCH_PDIM*2), f
        integer  :: xi, yi, j
        real(dp) :: x, y
        real(dp) :: ashift(2), t, ftmp, grad_contrib(POLYDIM2), grad_full_tmp(PATCH_PDIM*2)
        write (*,*) 'fdf, patched direct, vec=', vec
        f       = 0.d0
        grad(:) = 0.d0
        ! convert polynomial coefficients into shifts
        !$omp parallel do collapse(2) default(shared) private(xi,yi,j,t,x,y,ashift,ftmp,grad_contrib,grad_full_tmp) reduction(+:f,grad) proc_bind(close) schedule(static)
        do xi = 1, params_glob%nxpatch
            do yi = 1, params_glob%nypatch
                do j = 1, self%nframes
                    t = real(j - 1, dp)
                    call self%pix2polycoords(real(self%patch_centers(xi,yi,1),dp), real(self%patch_centers(xi,yi,2),dp), x, y) ! TODO: memoize
                    ashift(1) = apply_patch_poly(vec(           1:PATCH_PDIM  ), x, y, t)
                    ashift(2) = apply_patch_poly(vec(PATCH_PDIM+1:PATCH_PDIM*2), x, y, t)
                    self%align_iso_polyn_direct(xi,yi)%shifts(j,:) = - ashift(:)
                end do
                call self%align_iso_polyn_direct(xi,yi)%motion_align_iso_contribs_fdf(ftmp, grad_contrib)
                f = f + ftmp
                call self%grad_contrib_to_full_grad( grad_contrib, x, y, grad_full_tmp )
                grad(:) = grad(:) + (- grad_full_tmp(:))
            end do
        end do
        !$omp end parallel do
    end subroutine motion_patched_direct_fdf

    subroutine patched_direct_fdf_wrapper( self, vec, f, grad, D )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
        f    = 0.d0
        grad = 0.d0
        select type(self)
            class is (motion_patched)
                call self%motion_patched_direct_fdf( vec, f, grad )
            class DEFAULT
                THROW_HARD('unknown type; patched_direct_fdf_wrapper')
        end select
    end subroutine patched_direct_fdf_wrapper


end module simple_motion_patched
