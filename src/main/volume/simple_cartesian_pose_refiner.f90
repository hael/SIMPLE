!@descr: Fixed-volume Cartesian Fourier particle-pose refinement numerics
module simple_cartesian_pose_refiner
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_image, only: image
use simple_ctf, only: ctf
use simple_cartesian_fourier, only: center_embed_real3d, gather_packed_window_grad
use simple_gridding, only: kb_stencil_centered_crop_inv_envelope_1d
implicit none

public :: cartesian_pose_refiner, cartesian_pose_data
public :: POSE_DATA_VALID, POSE_DATA_INVALID_NOISE_RANGE
public :: SHIFT_LM_ACCEPTED_IMPROVEMENT, SHIFT_LM_FINITE_NO_IMPROVEMENT
public :: SHIFT_LM_NO_RELIABLE_UPDATE, SHIFT_LM_STEP_BOUND_REJECTED
public :: SHIFT_LM_INVALID_NUMERICS, SHIFT_LM_ITERATION_LIMIT
public :: POSE_LM_ACCEPTED_IMPROVEMENT, POSE_LM_FINITE_NO_IMPROVEMENT
public :: POSE_LM_NO_RELIABLE_UPDATE, POSE_LM_STEP_BOUND_REJECTED
public :: POSE_LM_INVALID_NUMERICS, POSE_LM_ITERATION_LIMIT
public :: right_increment_rotation
private
#include "simple_local_flags.inc"

integer, parameter :: SHIFT_LM_ACCEPTED_IMPROVEMENT = 1
integer, parameter :: SHIFT_LM_FINITE_NO_IMPROVEMENT = 2
integer, parameter :: SHIFT_LM_NO_RELIABLE_UPDATE = 3
integer, parameter :: SHIFT_LM_STEP_BOUND_REJECTED = 4
integer, parameter :: SHIFT_LM_INVALID_NUMERICS = 5
integer, parameter :: SHIFT_LM_ITERATION_LIMIT = 6
integer, parameter :: POSE_LM_ACCEPTED_IMPROVEMENT = SHIFT_LM_ACCEPTED_IMPROVEMENT
integer, parameter :: POSE_LM_FINITE_NO_IMPROVEMENT = SHIFT_LM_FINITE_NO_IMPROVEMENT
integer, parameter :: POSE_LM_NO_RELIABLE_UPDATE = SHIFT_LM_NO_RELIABLE_UPDATE
integer, parameter :: POSE_LM_STEP_BOUND_REJECTED = SHIFT_LM_STEP_BOUND_REJECTED
integer, parameter :: POSE_LM_INVALID_NUMERICS = SHIFT_LM_INVALID_NUMERICS
integer, parameter :: POSE_LM_ITERATION_LIMIT = SHIFT_LM_ITERATION_LIMIT
real(dp), parameter :: POSE_NUMERIC_FLOOR = epsilon(1._dp)**2
integer, parameter :: POSE_DATA_VALID = 0
integer, parameter :: POSE_DATA_INVALID_NOISE_RANGE = 1

!> One shift-free, noise-whitened particle observation for a fixed reference.
type :: cartesian_pose_data
    private
    complex, allocatable :: observed(:,:), transfer(:,:)
    integer :: shell_range(2) = 0
    logical :: valid = .false.
contains
    procedure :: is_valid => pose_data_is_valid
    procedure :: get_shell_range => pose_data_get_shell_range
end type cartesian_pose_data

type :: cartesian_pose_refiner
    private
    integer          :: box     = 0
    integer          :: boxpd   = 0
    integer          :: padf    = 1
    integer          :: iwinsz  = 0
    integer          :: wdim    = 0
    integer          :: lims2(2,2) = 0
    integer          :: sqhp    = 0
    integer          :: sqlp    = 0
    real             :: padsc   = 1.0
    type(kbinterpol) :: kbwin
    integer, allocatable :: wrap(:)
    complex, allocatable :: cmat(:,:,:)
    logical :: exists = .false.
contains
    procedure :: new => new_pose_refiner
    procedure :: new_prepared_test => new_pose_refiner_prepared_test
    procedure :: prepare_particle => prepare_pose_particle
    procedure :: prepared_objective_gradient => prepared_pose_objective_gradient
    procedure :: refine_prepared_pose_lm
    procedure :: kill => kill_fourier_workspace
    procedure :: get_lims2 => get_fourier_workspace_lims2
    procedure :: set_shell_range => set_fourier_workspace_shell_range
    procedure :: sample_with_grad => sample_fourier_with_grad
    procedure :: shift_residual
    procedure :: shift_jvp
    procedure :: shift_jhz
    procedure :: shift_objective_gradient
    procedure :: refine_shift_lm
    procedure :: rotation_jvp
    procedure :: pose_objective_gradient
    procedure :: refine_pose_lm
    procedure, private :: shift_normal_terms
    procedure, private :: pose_normal_terms
    procedure :: count_stencil_switches
end type cartesian_pose_refiner

contains

    pure logical function pose_data_is_valid(self) result(valid)
        class(cartesian_pose_data), intent(in) :: self
        valid = self%valid
    end function pose_data_is_valid

    pure function pose_data_get_shell_range(self) result(shell_range)
        class(cartesian_pose_data), intent(in) :: self
        integer :: shell_range(2)
        shell_range = self%shell_range
    end function pose_data_get_shell_range

    !> Construct one immutable Fourier reference from a physical real-space volume.
    subroutine new_pose_refiner(self, volume)
        class(cartesian_pose_refiner), intent(inout) :: self
        real, intent(in) :: volume(:,:,:)
        call load_pose_reference(self,volume,.true.)
    end subroutine new_pose_refiner

    !> Test-only constructor for diagnostics that compare prepared-volume models.
    subroutine new_pose_refiner_prepared_test(self, volume)
        class(cartesian_pose_refiner), intent(inout) :: self
        real, intent(in) :: volume(:,:,:)
        call load_pose_reference(self,volume,.false.)
    end subroutine new_pose_refiner_prepared_test

    subroutine load_pose_reference(self, volume, apply_inverse_envelope)
        class(cartesian_pose_refiner), intent(inout) :: self
        real, intent(in) :: volume(:,:,:)
        logical, intent(in) :: apply_inverse_envelope
        type(image) :: padded_image
        real, allocatable :: prepared(:,:,:), inv1d(:)
        integer :: box, lims3(3,2), wlims(2), lo, hi, i, j, k

        call self%kill
        box = size(volume,1)
        if( box < 2 .or. mod(box,2) /= 0 .or. size(volume,2) /= box .or. size(volume,3) /= box ) &
            &error stop 'cartesian pose reference requires an even cubic volume'
        self%box = box
        self%padf = OSMPL_PAD_FAC
        self%boxpd = self%padf*box
        self%padsc = real(self%padf)**3
        self%kbwin = kbinterpol(KBWINSZ,KBALPHA)
        self%iwinsz = ceiling(self%kbwin%get_winsz()-0.5)
        self%wdim = 2*self%iwinsz+1
        self%lims2(1,:) = [-box/2,box/2]
        self%lims2(2,:) = [-box/2,box/2]
        self%sqhp = 0
        self%sqlp = (box/2)**2
        allocate(prepared,source=volume)
        if( apply_inverse_envelope )then
            call kb_stencil_centered_crop_inv_envelope_1d(self%kbwin,self%boxpd,box,inv1d)
            do k = 1, box
                do j = 1, box
                    do i = 1, box
                        prepared(i,j,k) = prepared(i,j,k)*inv1d(i)*inv1d(j)*inv1d(k)
                    enddo
                enddo
            enddo
            deallocate(inv1d)
        endif
        call padded_image%new([self%boxpd,self%boxpd,self%boxpd],1.0)
        call padded_image%set_rmat(center_embed_real3d(prepared,self%boxpd),.false.)
        call padded_image%fft()
        lims3 = padded_image%loop_lims(3)
        wlims = lims3(2,:)
        lo = wlims(1)-self%iwinsz-1
        hi = wlims(2)+self%iwinsz+1
        allocate(self%wrap(lo:hi))
        do i = lo, hi
            self%wrap(i) = cyci_1d(wlims,i)
        enddo
        self%cmat = padded_image%get_cmat()
        self%exists = .true.
        call padded_image%kill
        deallocate(prepared)
    end subroutine load_pose_reference

    !> Prepare the shift-free CTF/noise transfer and whitened observation.
    subroutine prepare_pose_particle(self, raw_observed, ctfparms, sigma2, requested_range, data, status)
        class(cartesian_pose_refiner), intent(in) :: self
        complex, intent(in) :: raw_observed(self%lims2(1,1):self%lims2(1,2), &
            &self%lims2(2,1):self%lims2(2,2))
        type(ctfparams), intent(in) :: ctfparms
        real, intent(in) :: sigma2(0:)
        integer, intent(in) :: requested_range(2)
        type(cartesian_pose_data), intent(out) :: data
        integer, intent(out) :: status
        type(ctf) :: tfun
        type(ctfvars) :: ctfvals
        real :: cval, sigma, s2, cterm, df, phsh
        real :: wl, half_wl2_cs, sum_df, diff_df, angast, accc, phc
        integer :: h, k, shell, lower_shell, upper_shell
        logical :: use_ctf, phase_flip

        status = POSE_DATA_INVALID_NOISE_RANGE
        data%valid = .false.
        data%shell_range = 0
        lower_shell = max(0,requested_range(1))
        upper_shell = min(requested_range(2),ubound(sigma2,1),self%box/2)
        if( upper_shell < lower_shell ) return
        do shell = lower_shell, upper_shell
            if( .not. ieee_is_finite(sigma2(shell)) .or. sigma2(shell) <= 0.0 ) return
        enddo
        allocate(data%observed(self%lims2(1,1):self%lims2(1,2), &
            &self%lims2(2,1):self%lims2(2,2)),source=cmplx(0.,0.))
        allocate(data%transfer(self%lims2(1,1):self%lims2(1,2), &
            &self%lims2(2,1):self%lims2(2,2)),source=cmplx(0.,0.))
        use_ctf = ctfparms%ctfflag /= CTFFLAG_NO
        phase_flip = ctfparms%ctfflag == CTFFLAG_FLIP
        if( use_ctf )then
            tfun = ctf(ctfparms%smpd,ctfparms%kv,ctfparms%cs,ctfparms%fraca)
            call tfun%init(ctfparms%dfx,ctfparms%dfy,ctfparms%angast)
            ctfvals = tfun%get_ctfvars(ctfparms%phshift)
            wl = ctfvals%wl
            half_wl2_cs = 0.5*wl*wl*ctfvals%cs
            sum_df = ctfvals%dfx+ctfvals%dfy
            diff_df = ctfvals%dfx-ctfvals%dfy
            angast = ctfvals%angast
            accc = ctfvals%amp_contr_const
            phc = ctfvals%phshift
        endif
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                shell = nint(sqrt(real(h*h+k*k)))
                if( shell < lower_shell .or. shell > upper_shell ) cycle
                sigma = sigma2(shell)
                cval = 1.0
                if( use_ctf )then
                    s2 = real(h*h+k*k)/(real(self%box)*ctfparms%smpd)**2
                    cterm = cos(2.0*(atan2(real(k),real(h))-angast))
                    df = 0.5*(sum_df+cterm*diff_df)
                    phsh = PI*wl*s2*(df-half_wl2_cs*s2)
                    cval = sin(phsh+phc+accc)
                    if( phase_flip ) cval = abs(cval)
                endif
                data%transfer(h,k) = cval/sqrt(sigma)
                data%observed(h,k) = raw_observed(h,k)/sqrt(sigma)
            enddo
        enddo
        data%shell_range = [lower_shell,upper_shell]
        data%valid = .true.
        status = POSE_DATA_VALID
    end subroutine prepare_pose_particle

    subroutine prepared_pose_objective_gradient(self, rotmat, shift, data, objective, gradient)
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        type(cartesian_pose_data), intent(in) :: data
        real(dp), intent(out) :: objective, gradient(5)
        if( .not. data%valid ) error stop 'prepared pose objective requires valid particle data'
        ! Reuse the immutable Fourier lattice; only the particle shell mask varies.
        call self%pose_objective_gradient(rotmat,shift,data%observed,objective,gradient, &
            &data%transfer,data%shell_range)
    end subroutine prepared_pose_objective_gradient

    subroutine refine_prepared_pose_lm(self, rotmat, shift, data, rotation_scale, max_iterations, &
        &accepted_objectives, naccepted, status, nattempted, max_rotation_step, max_shift_step, &
        &nstencil_switches)
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(inout) :: rotmat(3,3), shift(2)
        type(cartesian_pose_data), intent(in) :: data
        real(dp), intent(in) :: rotation_scale
        integer, intent(in) :: max_iterations
        real(dp), intent(out) :: accepted_objectives(0:)
        integer, intent(out) :: naccepted, status, nattempted, nstencil_switches
        real(dp), intent(out) :: max_rotation_step, max_shift_step
        if( .not. data%valid ) error stop 'prepared pose LM requires valid particle data'
        ! Shell limits are local arguments, so no padded Fourier reference is copied.
        call self%refine_pose_lm(rotmat,data%observed,shift,rotation_scale,max_iterations, &
            &accepted_objectives,naccepted,status,nattempted,max_rotation_step,max_shift_step, &
            &nstencil_switches,transfer=data%transfer,shell_range=data%shell_range)
    end subroutine refine_prepared_pose_lm

    pure function right_increment_rotation( rotmat, omega ) result(updated_rotmat)
        real(dp), intent(in) :: rotmat(3,3), omega(3)
        real(dp) :: updated_rotmat(3,3), skew(3,3), exp_skew(3,3)
        real(dp) :: identity(3,3), theta2, theta4, sinc_theta, cosc_theta

        identity = 0._dp
        identity(1,1) = 1._dp
        identity(2,2) = 1._dp
        identity(3,3) = 1._dp
        ! [omega]x u = omega x u.
        skew = reshape([0._dp,omega(3),-omega(2), &
            &-omega(3),0._dp,omega(1),omega(2),-omega(1),0._dp],[3,3])
        theta2 = dot_product(omega,omega)
        if( theta2 < 1.e-8_dp )then
            theta4 = theta2*theta2
            sinc_theta = 1._dp-theta2/6._dp+theta4/120._dp
            cosc_theta = 0.5_dp-theta2/24._dp+theta4/720._dp
        else
            sinc_theta = sin(sqrt(theta2))/sqrt(theta2)
            cosc_theta = (1._dp-cos(sqrt(theta2)))/theta2
        endif
        exp_skew = identity+sinc_theta*skew+cosc_theta*matmul(skew,skew)
        updated_rotmat = matmul(rotmat,exp_skew)
    end function right_increment_rotation

    subroutine kill_fourier_workspace( self )
        class(cartesian_pose_refiner), intent(inout) :: self
        if( allocated(self%wrap) ) deallocate(self%wrap)
        if( allocated(self%cmat) ) deallocate(self%cmat)
        self%box    = 0
        self%boxpd  = 0
        self%padf   = 1
        self%iwinsz = 0
        self%wdim   = 0
        self%lims2  = 0
        self%sqhp   = 0
        self%sqlp   = 0
        self%padsc  = 1.0
        self%exists = .false.
    end subroutine kill_fourier_workspace

    !> Return the native packed-plane bounds used by observations and transfers.
    pure function get_fourier_workspace_lims2( self ) result(lims2)
        class(cartesian_pose_refiner), intent(in) :: self
        integer :: lims2(2,2)
        lims2 = self%lims2
    end function get_fourier_workspace_lims2

    !> Restrict the local objective to an inclusive native Fourier-shell range.
    subroutine set_fourier_workspace_shell_range( self, kfromto )
        class(cartesian_pose_refiner), intent(inout) :: self
        integer, intent(in) :: kfromto(2)
        integer :: khi, klo
        if( .not. self%exists ) error stop 'set_shell_range called on an empty Fourier workspace'
        klo = max(0,kfromto(1))
        khi = min(self%box/2,kfromto(2))
        if( khi < klo ) error stop 'set_shell_range requires an ordered nonempty range'
        self%sqhp = klo*klo
        self%sqlp = khi*khi
    end subroutine set_fourier_workspace_shell_range

    !>  \brief  Samples the packed Fourier snapshot and its three fixed-cell
    !!          spatial derivatives at one oversampled-lattice coordinate.
    pure subroutine sample_fourier_with_grad( self, loc, value, dvalue_dloc, switch_margin )
        class(cartesian_pose_refiner), intent(in)  :: self
        real(sp),                     intent(in)  :: loc(3)
        complex,                      intent(out) :: value, dvalue_dloc(3)
        real(sp),                     intent(out) :: switch_margin(3)
        real(sp) :: w(self%wdim,self%wdim,self%wdim)
        real(sp) :: dw(self%wdim,self%wdim,self%wdim,3)
        integer  :: i0(3)
        if( .not. self%exists ) error stop 'sample_with_grad called on an empty Fourier workspace'
        ! Build w and dw/dloc on the same fixed interpolation stencil.
        call self%kbwin%apod_mat_3d_fast_grad(loc, self%iwinsz, self%wdim, i0, switch_margin, w, dw)
        if( any(i0 < lbound(self%wrap,1)) .or. &
            &any(i0 + self%wdim - 1 > ubound(self%wrap,1)) )then
            error stop 'sample_with_grad location lies outside the periodic wrap table'
        endif
        call gather_packed_window_grad(self%cmat,lbound(self%wrap,1),self%wrap, &
            &i0,w,dw,value,dvalue_dloc)
        ! Apply native Fourier scaling to the value and all derivatives.
        value       = self%padsc * value
        dvalue_dloc = self%padsc * dvalue_dloc
    end subroutine sample_fourier_with_grad

    !>  \brief  Fixed-volume, CTF-free, unit-noise residual for a shifted
    !!          Fourier projection: r = S(t) G(R)V - y.
    subroutine shift_residual( self, rotmat, shift, observed, residual, objective )
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        complex, intent(in)  :: observed(self%lims2(1,1):self%lims2(1,2),&
                                          &self%lims2(2,1):self%lims2(2,2))
        complex, intent(out) :: residual(self%lims2(1,1):self%lims2(1,2),&
                                          &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(out) :: objective
        complex :: value, dvalue_dloc(3), phase
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: arg
        integer :: h, k
        if( .not. self%exists ) error stop 'shift_residual called on an empty Fourier workspace'
        residual = cmplx(0.,0.)
        objective = 0._dp
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k < self%sqhp .or. h*h + k*k > self%sqlp ) cycle
                loc = real(self%padf,sp) * real(matmul(real([h,k,0],dp),rotmat),sp)
                call self%sample_with_grad(loc, value, dvalue_dloc, switch_margin)
                arg = 2._dp * real(PI,dp) * &
                    &(real(h,dp)*shift(1) + real(k,dp)*shift(2)) / real(self%box,dp)
                phase = cmplx(cos(arg),sin(arg),kind=sp)
                residual(h,k) = phase * value - observed(h,k)
                objective = objective + 0.5_dp * real(conjg(cmplx(residual(h,k),kind=dp)) * &
                    &cmplx(residual(h,k),kind=dp),dp)
            enddo
        enddo
    end subroutine shift_residual

    !>  \brief  Directional derivative of the CTF-free, unit-noise residual
    !!          with respect to the two real image-shift parameters.
    subroutine shift_jvp( self, rotmat, shift, direction, jv )
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2), direction(2)
        complex, intent(out) :: jv(self%lims2(1,1):self%lims2(1,2),&
                                    &self%lims2(2,1):self%lims2(2,2))
        complex :: value, dvalue_dloc(3), phase
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: arg, directional_frequency
        integer :: h, k
        if( .not. self%exists ) error stop 'shift_jvp called on an empty Fourier workspace'
        jv = cmplx(0.,0.)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k < self%sqhp .or. h*h + k*k > self%sqlp ) cycle
                loc = real(self%padf,sp) * real(matmul(real([h,k,0],dp),rotmat),sp)
                call self%sample_with_grad(loc, value, dvalue_dloc, switch_margin)
                arg = 2._dp * real(PI,dp) * &
                    &(real(h,dp)*shift(1) + real(k,dp)*shift(2)) / real(self%box,dp)
                directional_frequency = 2._dp * real(PI,dp) * &
                    &(real(h,dp)*direction(1) + real(k,dp)*direction(2)) / real(self%box,dp)
                phase = cmplx(cos(arg),sin(arg),kind=sp)
                jv(h,k) = cmplx(cmplx(0._dp,directional_frequency,kind=dp) * &
                    &cmplx(phase*value,kind=dp),kind=sp)
            enddo
        enddo
    end subroutine shift_jvp

    !>  \brief  Real-parameter adjoint of the two shift-Jacobian columns.
    subroutine shift_jhz( self, rotmat, shift, z, jhz )
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        complex, intent(in) :: z(self%lims2(1,1):self%lims2(1,2),&
                                  &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(out) :: jhz(2)
        complex :: value, dvalue_dloc(3), phase
        complex(dp) :: jacobian_value
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: arg, frequency
        integer :: axis, h, k
        if( .not. self%exists ) error stop 'shift_jhz called on an empty Fourier workspace'
        jhz = 0._dp
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k < self%sqhp .or. h*h + k*k > self%sqlp ) cycle
                loc = real(self%padf,sp) * real(matmul(real([h,k,0],dp),rotmat),sp)
                call self%sample_with_grad(loc, value, dvalue_dloc, switch_margin)
                arg = 2._dp * real(PI,dp) * &
                    &(real(h,dp)*shift(1) + real(k,dp)*shift(2)) / real(self%box,dp)
                phase = cmplx(cos(arg),sin(arg),kind=sp)
                do axis = 1, 2
                    frequency = 2._dp * real(PI,dp) * real(merge(h,k,axis==1),dp) / real(self%box,dp)
                    jacobian_value = cmplx(0._dp,frequency,kind=dp) * &
                        &cmplx(phase*value,kind=dp)
                    jhz(axis) = jhz(axis) + real(conjg(jacobian_value) * &
                        &cmplx(z(h,k),kind=dp),dp)
                enddo
            enddo
        enddo
    end subroutine shift_jhz

    !>  \brief  Fused shift objective, gradient and Gauss-Newton block. This
    !!          avoids materializing derivative planes inside local refinement.
    subroutine shift_normal_terms( self, rotmat, shift, observed, objective, gradient, hessian, transfer )
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        complex, intent(in) :: observed(self%lims2(1,1):self%lims2(1,2),&
                                         &self%lims2(2,1):self%lims2(2,2))
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2),&
                                                   &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(out) :: objective, gradient(2), hessian(2,2)
        complex :: value, dvalue_dloc(3), phase
        complex(dp) :: model, residual, jacobian(2)
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: arg, frequency(2)
        integer :: axis, h, jaxis, k
        if( .not. self%exists ) error stop 'shift_normal_terms called on an empty Fourier workspace'
        objective = 0._dp
        gradient = 0._dp
        hessian = 0._dp
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k < self%sqhp .or. h*h + k*k > self%sqlp ) cycle
                loc = real(self%padf,sp) * real(matmul(real([h,k,0],dp),rotmat),sp)
                call self%sample_with_grad(loc, value, dvalue_dloc, switch_margin)
                arg = 2._dp * real(PI,dp) * &
                    &(real(h,dp)*shift(1) + real(k,dp)*shift(2)) / real(self%box,dp)
                phase = cmplx(cos(arg),sin(arg),kind=sp)
                model = cmplx(phase*value,kind=dp)
                if( present(transfer) ) model = model * cmplx(transfer(h,k),kind=dp)
                residual = model - cmplx(observed(h,k),kind=dp)
                frequency = 2._dp * real(PI,dp) * real([h,k],dp) / real(self%box,dp)
                jacobian = cmplx(0._dp,frequency,kind=dp) * model
                objective = objective + 0.5_dp*real(conjg(residual)*residual,dp)
                do axis = 1, 2
                    gradient(axis) = gradient(axis) + real(conjg(jacobian(axis))*residual,dp)
                    do jaxis = 1, 2
                        hessian(axis,jaxis) = hessian(axis,jaxis) + &
                            &real(conjg(jacobian(axis))*jacobian(jaxis),dp)
                    enddo
                enddo
            enddo
        enddo
    end subroutine shift_normal_terms

    !> Return the shift objective and gradient without exposing the computed
    !! two-by-two Gauss-Newton block.
    subroutine shift_objective_gradient( self, rotmat, shift, observed, objective, gradient, transfer )
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        complex, intent(in) :: observed(self%lims2(1,1):self%lims2(1,2), &
                                         &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(out) :: objective, gradient(2)
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2), &
                                                   &self%lims2(2,1):self%lims2(2,2))
        real(dp) :: hessian(2,2)
        if( present(transfer) )then
            call self%shift_normal_terms(rotmat,shift,observed,objective,gradient,hessian,transfer)
        else
            call self%shift_normal_terms(rotmat,shift,observed,objective,gradient,hessian)
        endif
    end subroutine shift_objective_gradient

    !>  \brief  Directional derivative of the Fourier gather for a right
    !!          tangent-space rotation. The transfer plane includes CTF and
    !!          whitening when it is supplied by the production caller.
    subroutine rotation_jvp( self, rotmat, shift, direction, jv, transfer )
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2), direction(3)
        complex, intent(out) :: jv(self%lims2(1,1):self%lims2(1,2),&
                                    &self%lims2(2,1):self%lims2(2,2))
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2),&
                                                   &self%lims2(2,1):self%lims2(2,2))
        complex :: value, dvalue_dloc(3), phase
        complex(dp) :: derivative
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: args, dloc(3)
        integer :: h, k

        if( .not. self%exists ) error stop 'rotation_jvp called on an empty Fourier workspace'
        jv = cmplx(0.,0.)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h+k*k < self%sqhp .or. h*h+k*k > self%sqlp ) cycle
                loc = real(self%padf,sp)*real(matmul(real([h,k,0],dp),rotmat),sp)
                call self%sample_with_grad(loc,value,dvalue_dloc,switch_margin)
                ! For row-vector gathers, d(loc)/d(epsilon) = loc x direction.
                dloc = [real(loc(2),dp)*direction(3)-real(loc(3),dp)*direction(2), &
                    &real(loc(3),dp)*direction(1)-real(loc(1),dp)*direction(3), &
                    &real(loc(1),dp)*direction(2)-real(loc(2),dp)*direction(1)]
                args = 2._dp*real(PI,dp)*(real(h,dp)*shift(1)+real(k,dp)*shift(2))/real(self%box,dp)
                phase = cmplx(cos(args),sin(args),kind=sp)
                derivative = cmplx(phase,kind=dp)*sum(cmplx(dvalue_dloc,kind=dp)*dloc)
                if( present(transfer) ) derivative = derivative*cmplx(transfer(h,k),kind=dp)
                jv(h,k) = cmplx(derivative,kind=sp)
            enddo
        enddo
    end subroutine rotation_jvp

    !>  \brief  Fused objective, five-vector gradient and Gauss-Newton block
    !!          for three right-rotation coordinates and two pixel shifts.
    subroutine pose_normal_terms( self, rotmat, shift, observed, objective, gradient, &
        &hessian, min_switch_margin, transfer, shell_range )
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        complex, intent(in) :: observed(self%lims2(1,1):self%lims2(1,2),&
                                         &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(out) :: objective, gradient(5), hessian(5,5), min_switch_margin
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2),&
                                                   &self%lims2(2,1):self%lims2(2,2))
        integer, optional, intent(in) :: shell_range(2)
        complex :: value, dvalue_dloc(3), phase
        complex(dp) :: weighted_phase, model, residual, jacobian(5)
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: arg, dloc(3,3), frequency(2)
        integer :: active_sqhp, active_sqlp, axis, h, jaxis, k

        if( .not. self%exists ) error stop 'pose_normal_terms called on an empty Fourier workspace'
        objective = 0._dp
        gradient = 0._dp
        hessian = 0._dp
        min_switch_margin = huge(0._dp)
        active_sqhp = self%sqhp
        active_sqlp = self%sqlp
        if( present(shell_range) )then
            if( shell_range(1) < 0 .or. shell_range(2) > self%box/2 .or. &
                &shell_range(2) < shell_range(1) ) &
                &error stop 'pose shell range lies outside the native Fourier disk'
            active_sqhp = shell_range(1)*shell_range(1)
            active_sqlp = shell_range(2)*shell_range(2)
        endif
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h+k*k < active_sqhp .or. h*h+k*k > active_sqlp ) cycle
                loc = real(self%padf,sp)*real(matmul(real([h,k,0],dp),rotmat),sp)
                call self%sample_with_grad(loc,value,dvalue_dloc,switch_margin)
                min_switch_margin = min(min_switch_margin,real(minval(switch_margin),dp))
                ! Columns are loc x e1, loc x e2 and loc x e3.
                dloc(:,1) = [0._dp,real(loc(3),dp),-real(loc(2),dp)]
                dloc(:,2) = [-real(loc(3),dp),0._dp,real(loc(1),dp)]
                dloc(:,3) = [real(loc(2),dp),-real(loc(1),dp),0._dp]
                arg = 2._dp*real(PI,dp)*(real(h,dp)*shift(1)+real(k,dp)*shift(2))/real(self%box,dp)
                phase = cmplx(cos(arg),sin(arg),kind=sp)
                weighted_phase = cmplx(phase,kind=dp)
                if( present(transfer) ) weighted_phase = weighted_phase*cmplx(transfer(h,k),kind=dp)
                model = weighted_phase*cmplx(value,kind=dp)
                residual = model-cmplx(observed(h,k),kind=dp)
                do axis = 1, 3
                    jacobian(axis) = weighted_phase*sum(cmplx(dvalue_dloc,kind=dp)*dloc(:,axis))
                enddo
                frequency = 2._dp*real(PI,dp)*real([h,k],dp)/real(self%box,dp)
                jacobian(4:5) = cmplx(0._dp,frequency,kind=dp)*model
                objective = objective+0.5_dp*real(conjg(residual)*residual,dp)
                do axis = 1, 5
                    gradient(axis) = gradient(axis)+real(conjg(jacobian(axis))*residual,dp)
                    do jaxis = 1, 5
                        hessian(axis,jaxis) = hessian(axis,jaxis)+ &
                            &real(conjg(jacobian(axis))*jacobian(jaxis),dp)
                    enddo
                enddo
            enddo
        enddo
        if( min_switch_margin == huge(0._dp) ) min_switch_margin = 0._dp
    end subroutine pose_normal_terms

    !> Return the joint pose objective and five-vector gradient without exposing
    !! the computed Gauss-Newton block or stencil margin.
    subroutine pose_objective_gradient( self, rotmat, shift, observed, objective, gradient, &
        &transfer, shell_range )
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        complex, intent(in) :: observed(self%lims2(1,1):self%lims2(1,2),&
                                         &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(out) :: objective, gradient(5)
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2),&
                                                   &self%lims2(2,1):self%lims2(2,2))
        integer, optional, intent(in) :: shell_range(2)
        real(dp) :: hessian(5,5), min_switch_margin

        call self%pose_normal_terms(rotmat,shift,observed,objective,gradient,hessian, &
            &min_switch_margin,transfer,shell_range)
    end subroutine pose_objective_gradient

    !>  \brief  Counts active Fourier samples whose nearest-grid interpolation
    !!          stencil changes between two rotation matrices.
    function count_stencil_switches( self, rotmat, trial_rotmat, shell_range ) result(nswitches)
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), trial_rotmat(3,3)
        integer, optional, intent(in) :: shell_range(2)
        integer :: nswitches
        real(dp) :: loc(3), trial_loc(3)
        integer :: active_sqhp, active_sqlp, h, k

        if( .not. self%exists ) error stop 'count_stencil_switches called on an empty Fourier workspace'
        nswitches = 0
        active_sqhp = self%sqhp
        active_sqlp = self%sqlp
        if( present(shell_range) )then
            if( shell_range(1) < 0 .or. shell_range(2) > self%box/2 .or. &
                &shell_range(2) < shell_range(1) ) &
                &error stop 'stencil-switch shell range lies outside the native Fourier disk'
            active_sqhp = shell_range(1)*shell_range(1)
            active_sqlp = shell_range(2)*shell_range(2)
        endif
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h+k*k < active_sqhp .or. h*h+k*k > active_sqlp ) cycle
                loc = real(self%padf,dp)*matmul(real([h,k,0],dp),rotmat)
                trial_loc = real(self%padf,dp)*matmul(real([h,k,0],dp),trial_rotmat)
                if( any(nint(loc) /= nint(trial_loc)) ) nswitches = nswitches+1
            enddo
        enddo
    end function count_stencil_switches

    !>  \brief  Damped two-parameter Gauss-Newton refinement. Only accepted,
    !!          fully recomputed objective values are appended to the trace.
    subroutine refine_shift_lm( self, rotmat, observed, shift, max_iterations, &
        &accepted_objectives, naccepted, status, nattempted, max_trial_step, transfer )
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3)
        complex, intent(in) :: observed(self%lims2(1,1):self%lims2(1,2),&
                                         &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(inout) :: shift(2)
        integer, intent(in) :: max_iterations
        real(dp), intent(out) :: accepted_objectives(0:)
        integer, intent(out) :: naccepted
        integer, intent(out), optional :: status, nattempted
        real(dp), intent(out), optional :: max_trial_step
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2),&
                                                   &self%lims2(2,1):self%lims2(2,2))
        real(dp) :: gradient(2), hessian(2,2), trial_gradient(2), trial_hessian(2,2)
        real(dp) :: solve_matrix(2,2), diagonal(2), direction(2), trial_shift(2)
        real(dp) :: objective, trial_objective, mu, det, predicted, actual, ratio, maxdiag
        real(dp) :: discriminant, lambda_max, lambda_min, step_norm, relative_reduction
        integer :: axis, iteration, outcome, attempted
        logical :: bounded_trial
        if( ubound(accepted_objectives,1) < max_iterations )then
            error stop 'refine_shift_lm objective trace is shorter than max_iterations+1'
        endif
        if( present(transfer) )then
            call self%shift_normal_terms(rotmat,shift,observed,objective,gradient,hessian,transfer)
        else
            call self%shift_normal_terms(rotmat,shift,observed,objective,gradient,hessian)
        endif
        accepted_objectives = huge(0._dp)
        accepted_objectives(0) = objective
        naccepted = 0
        attempted = 0
        bounded_trial = .false.
        outcome = SHIFT_LM_ITERATION_LIMIT
        if( present(max_trial_step) ) max_trial_step = 0._dp
        if( .not. ieee_is_finite(objective) .or. any(.not. ieee_is_finite(gradient)) .or. &
            &any(.not. ieee_is_finite(hessian)) )then
            outcome = SHIFT_LM_INVALID_NUMERICS
            if( present(status) ) status = outcome
            if( present(nattempted) ) nattempted = attempted
            return
        endif
        mu = 1.e-3_dp
        do iteration = 1, max_iterations
            maxdiag = max(maxval([(hessian(axis,axis),axis=1,2)]),1._dp)
            ! Eigenvalues diagnose whether both shift directions are observable.
            discriminant = sqrt(max(0._dp,(hessian(1,1)-hessian(2,2))**2 + &
                &4._dp*hessian(1,2)*hessian(2,1)))
            lambda_max = 0.5_dp*(hessian(1,1)+hessian(2,2)+discriminant)
            lambda_min = 0.5_dp*(hessian(1,1)+hessian(2,2)-discriminant)
            if( lambda_max <= sqrt(epsilon(1._dp))*maxdiag .or. &
                &lambda_min <= sqrt(epsilon(1._dp))*max(lambda_max,1._dp) )then
                outcome = SHIFT_LM_NO_RELIABLE_UPDATE
                exit
            endif
            if( sqrt(dot_product(gradient,gradient)) < 1.e-8_dp )then
                outcome = merge(SHIFT_LM_ACCEPTED_IMPROVEMENT,SHIFT_LM_FINITE_NO_IMPROVEMENT,naccepted>0)
                exit
            endif
            do axis = 1, 2
                diagonal(axis) = max(hessian(axis,axis),sqrt(epsilon(1._dp))*maxdiag,epsilon(1._dp))
            enddo
            solve_matrix = hessian
            solve_matrix(1,1) = solve_matrix(1,1) + mu*diagonal(1)
            solve_matrix(2,2) = solve_matrix(2,2) + mu*diagonal(2)
            det = solve_matrix(1,1)*solve_matrix(2,2) - solve_matrix(1,2)*solve_matrix(2,1)
            if( abs(det) <= epsilon(1._dp)*maxdiag*maxdiag )then
                outcome = SHIFT_LM_NO_RELIABLE_UPDATE
                exit
            endif
            direction(1) = (-solve_matrix(2,2)*gradient(1) + solve_matrix(1,2)*gradient(2)) / det
            direction(2) = ( solve_matrix(2,1)*gradient(1) - solve_matrix(1,1)*gradient(2)) / det
            if( any(.not. ieee_is_finite(direction)) )then
                outcome = SHIFT_LM_INVALID_NUMERICS
                exit
            endif
            ! Shift coordinates are pixels; cap every trial displacement at one pixel.
            step_norm = sqrt(dot_product(direction,direction))
            if( step_norm > 1._dp )then
                direction = direction/step_norm
                bounded_trial = .true.
            endif
            step_norm = min(step_norm,1._dp)
            if( present(max_trial_step) ) max_trial_step = max(max_trial_step,step_norm)
            predicted = -dot_product(gradient,direction) - 0.5_dp * &
                &dot_product(direction,matmul(hessian,direction))
            if( .not. ieee_is_finite(predicted) )then
                outcome = SHIFT_LM_INVALID_NUMERICS
                exit
            elseif( predicted <= 0._dp )then
                mu = 4._dp * mu
                cycle
            endif
            trial_shift = shift + direction
            if( present(transfer) )then
                call self%shift_normal_terms(rotmat,trial_shift,observed,trial_objective,&
                    &trial_gradient,trial_hessian,transfer)
            else
                call self%shift_normal_terms(rotmat,trial_shift,observed,trial_objective,&
                    &trial_gradient,trial_hessian)
            endif
            attempted = attempted + 1
            if( .not. ieee_is_finite(trial_objective) .or. any(.not. ieee_is_finite(trial_gradient)) .or. &
                &any(.not. ieee_is_finite(trial_hessian)) )then
                mu = 4._dp * mu
                outcome = SHIFT_LM_INVALID_NUMERICS
                cycle
            endif
            actual = objective - trial_objective
            ratio = actual / predicted
            if( actual > 0._dp .and. ratio >= 0.25_dp )then
                relative_reduction = actual/max(abs(objective),1._dp)
                shift = trial_shift
                objective = trial_objective
                gradient = trial_gradient
                hessian = trial_hessian
                naccepted = naccepted + 1
                accepted_objectives(naccepted) = objective
                if( ratio > 0.75_dp ) mu = max(mu/2._dp,epsilon(1._dp))
                outcome = SHIFT_LM_ACCEPTED_IMPROVEMENT
                if( step_norm < 1.e-8_dp .or. relative_reduction < 1.e-10_dp ) exit
            else
                mu = 4._dp * mu
            endif
        enddo
        if( outcome == SHIFT_LM_ITERATION_LIMIT .and. naccepted == 0 .and. bounded_trial ) &
            &outcome = SHIFT_LM_STEP_BOUND_REJECTED
        if( present(status) ) status = outcome
        if( present(nattempted) ) nattempted = attempted
    end subroutine refine_shift_lm

    !>  \brief  Scaled, bounded five-parameter LM refinement for a right
    !!          rotation increment and two image shifts.
    subroutine refine_pose_lm( self, rotmat, observed, shift, rotation_scale, max_iterations, &
        &accepted_objectives, naccepted, status, nattempted, max_rotation_step, &
        &max_shift_step, nstencil_switches, transfer, accepted_rotmats, accepted_shifts, &
        &active_parameters, anchor_rotmat, anchor_shift, max_total_rotation, max_total_shift, &
        &shell_range )
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(inout) :: rotmat(3,3), shift(2)
        complex, intent(in) :: observed(self%lims2(1,1):self%lims2(1,2),&
                                         &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(in) :: rotation_scale
        integer, intent(in) :: max_iterations
        real(dp), intent(out) :: accepted_objectives(0:)
        integer, intent(out) :: naccepted, status, nattempted
        real(dp), intent(out) :: max_rotation_step, max_shift_step
        integer, intent(out) :: nstencil_switches
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2),&
                                                   &self%lims2(2,1):self%lims2(2,2))
        real(dp), optional, intent(out) :: accepted_rotmats(:,:,0:), accepted_shifts(:,0:)
        logical, optional, intent(in) :: active_parameters(5)
        real(dp), optional, intent(in) :: anchor_rotmat(3,3), anchor_shift(2)
        real(dp), optional, intent(in) :: max_total_rotation, max_total_shift
        integer, optional, intent(in) :: shell_range(2)
        real(dp) :: gradient(5), hessian(5,5), trial_gradient(5), trial_hessian(5,5)
        real(dp) :: scaled_gradient(5), scaled_hessian(5,5), solve_matrix(5,5)
        real(dp) :: diagonal(5), coordinate_scale(5), scaled_direction(5), direction(5)
        real(dp) :: trial_rotmat(3,3), trial_shift(2), objective, trial_objective
        real(dp) :: mu, predicted, actual, ratio, rotation_norm, shift_norm, hessian_scale
        real(dp) :: relative_reduction, min_switch_margin, trial_switch_margin
        real(dp) :: cumulative_rotation, cumulative_shift, sine_half
        integer :: axis, jaxis, iteration, trial_switches
        logical :: active(5), bounded_trial, cumulative_guard, reliable

        if( max_iterations < 1 ) error stop 'refine_pose_lm requires at least one LM iteration'
        if( rotation_scale <= 0._dp .or. .not. ieee_is_finite(rotation_scale) ) &
            &error stop 'refine_pose_lm requires a positive finite rotation scale'
        if( ubound(accepted_objectives,1) < max_iterations ) &
            &error stop 'refine_pose_lm objective trace is shorter than max_iterations+1'
        if( present(accepted_rotmats) )then
            if( size(accepted_rotmats,1) /= 3 .or. size(accepted_rotmats,2) /= 3 .or. &
                &ubound(accepted_rotmats,3) < max_iterations ) &
                &error stop 'refine_pose_lm rotation trace has invalid dimensions'
        endif
        if( present(accepted_shifts) )then
            if( size(accepted_shifts,1) /= 2 .or. ubound(accepted_shifts,2) < max_iterations ) &
                &error stop 'refine_pose_lm shift trace has invalid dimensions'
        endif
        active = .true.
        if( present(active_parameters) ) active = active_parameters
        if( .not. any(active) ) error stop 'refine_pose_lm requires one active parameter'
        cumulative_guard = present(anchor_rotmat) .and. present(anchor_shift) .and. &
            &present(max_total_rotation) .and. present(max_total_shift)
        if( cumulative_guard .neqv. (present(anchor_rotmat) .or. present(anchor_shift) .or. &
            &present(max_total_rotation) .or. present(max_total_shift)) ) &
            &error stop 'refine_pose_lm cumulative guard requires all four arguments'
        if( cumulative_guard )then
            if( max_total_rotation <= 0._dp .or. max_total_shift <= 0._dp ) &
                &error stop 'refine_pose_lm cumulative bounds must be positive'
        endif
        call self%pose_normal_terms(rotmat,shift,observed,objective,gradient,hessian, &
            &min_switch_margin,transfer,shell_range)
        accepted_objectives = huge(0._dp)
        accepted_objectives(0) = objective
        if( present(accepted_rotmats) ) accepted_rotmats(:,:,0) = rotmat
        if( present(accepted_shifts) ) accepted_shifts(:,0) = shift
        naccepted = 0
        nattempted = 0
        nstencil_switches = 0
        max_rotation_step = 0._dp
        max_shift_step = 0._dp
        bounded_trial = .false.
        status = POSE_LM_ITERATION_LIMIT
        if( .not. ieee_is_finite(objective) .or. any(.not. ieee_is_finite(gradient)) .or. &
            &any(.not. ieee_is_finite(hessian)) )then
            status = POSE_LM_INVALID_NUMERICS
            return
        endif

        ! Dimensionless variables balance radians against pixel shifts.
        coordinate_scale = [rotation_scale,rotation_scale,rotation_scale,1._dp,1._dp]
        do axis = 1, 5
            scaled_gradient(axis) = coordinate_scale(axis)*gradient(axis)
            do jaxis = 1, 5
                scaled_hessian(axis,jaxis) = &
                    &coordinate_scale(axis)*hessian(axis,jaxis)*coordinate_scale(jaxis)
            enddo
        enddo
        call apply_pose_parameter_mask(scaled_gradient,scaled_hessian,active)
        mu = 1.e-3_dp
        do iteration = 1, max_iterations
            ! Damping must not hide an unidentifiable five-parameter block.
            call solve_pose_cholesky(scaled_hessian,-scaled_gradient,scaled_direction,reliable)
            if( .not. reliable )then
                status = merge(POSE_LM_ACCEPTED_IMPROVEMENT,POSE_LM_NO_RELIABLE_UPDATE,naccepted>0)
                exit
            endif
            if( sqrt(dot_product(scaled_gradient,scaled_gradient)) < 1.e-8_dp )then
                status = merge(POSE_LM_ACCEPTED_IMPROVEMENT,POSE_LM_FINITE_NO_IMPROVEMENT,naccepted>0)
                exit
            endif
            hessian_scale = max(maxval(abs(scaled_hessian)),POSE_NUMERIC_FLOOR)
            do axis = 1, 5
                diagonal(axis) = max(scaled_hessian(axis,axis), &
                    &sqrt(epsilon(1._dp))*hessian_scale,POSE_NUMERIC_FLOOR)
            enddo
            solve_matrix = scaled_hessian
            do axis = 1, 5
                solve_matrix(axis,axis) = solve_matrix(axis,axis)+mu*diagonal(axis)
            enddo
            call solve_pose_cholesky(solve_matrix,-scaled_gradient,scaled_direction,reliable)
            if( .not. reliable )then
                status = POSE_LM_NO_RELIABLE_UPDATE
                exit
            endif
            direction = coordinate_scale*scaled_direction
            if( any(.not. ieee_is_finite(direction)) )then
                status = POSE_LM_INVALID_NUMERICS
                exit
            endif
            ! Bound rotations in radians and shifts in pixels independently.
            rotation_norm = sqrt(dot_product(direction(1:3),direction(1:3)))
            if( rotation_norm > rotation_scale )then
                direction(1:3) = direction(1:3)*(rotation_scale/rotation_norm)
                bounded_trial = .true.
            endif
            shift_norm = sqrt(dot_product(direction(4:5),direction(4:5)))
            if( shift_norm > 1._dp )then
                direction(4:5) = direction(4:5)/shift_norm
                bounded_trial = .true.
            endif
            rotation_norm = sqrt(dot_product(direction(1:3),direction(1:3)))
            shift_norm = sqrt(dot_product(direction(4:5),direction(4:5)))
            max_rotation_step = max(max_rotation_step,rotation_norm)
            max_shift_step = max(max_shift_step,shift_norm)
            ! Quadratic LM model: predicted = -g^T d - 1/2 d^T H d.
            predicted = -dot_product(gradient,direction)-0.5_dp* &
                &dot_product(direction,matmul(hessian,direction))
            if( .not. ieee_is_finite(predicted) )then
                status = POSE_LM_INVALID_NUMERICS
                exit
            elseif( predicted <= 0._dp )then
                mu = 4._dp*mu
                cycle
            endif
            trial_rotmat = right_increment_rotation(rotmat,direction(1:3))
            trial_shift = shift+direction(4:5)
            nattempted = nattempted+1
            if( cumulative_guard )then
                sine_half = sqrt(sum((trial_rotmat-anchor_rotmat)**2))/(2._dp*sqrt(2._dp))
                cumulative_rotation = 2._dp*asin(max(0._dp,min(1._dp,sine_half)))
                cumulative_shift = sqrt(sum((trial_shift-anchor_shift)**2))
                if( cumulative_rotation > max_total_rotation+10._dp*epsilon(1._dp) .or. &
                    &cumulative_shift > max_total_shift+10._dp*epsilon(1._dp) )then
                    mu = 4._dp*mu
                    bounded_trial = .true.
                    cycle
                endif
            endif
            trial_switches = self%count_stencil_switches(rotmat,trial_rotmat,shell_range)
            nstencil_switches = nstencil_switches+trial_switches
            call self%pose_normal_terms(trial_rotmat,trial_shift,observed,trial_objective, &
                &trial_gradient,trial_hessian,trial_switch_margin,transfer,shell_range)
            if( .not. ieee_is_finite(trial_objective) .or. any(.not. ieee_is_finite(trial_gradient)) .or. &
                &any(.not. ieee_is_finite(trial_hessian)) )then
                mu = 4._dp*mu
                status = POSE_LM_INVALID_NUMERICS
                cycle
            endif
            ! Gain ratio compares the recomputed reduction with the local model.
            actual = objective-trial_objective
            ratio = actual/predicted
            if( actual > 0._dp .and. ratio >= 0.25_dp )then
                relative_reduction = actual/max(abs(objective),1._dp)
                rotmat = trial_rotmat
                shift = trial_shift
                objective = trial_objective
                gradient = trial_gradient
                hessian = trial_hessian
                do axis = 1, 5
                    scaled_gradient(axis) = coordinate_scale(axis)*gradient(axis)
                    do jaxis = 1, 5
                        scaled_hessian(axis,jaxis) = &
                            &coordinate_scale(axis)*hessian(axis,jaxis)*coordinate_scale(jaxis)
                    enddo
                enddo
                call apply_pose_parameter_mask(scaled_gradient,scaled_hessian,active)
                naccepted = naccepted+1
                accepted_objectives(naccepted) = objective
                if( present(accepted_rotmats) ) accepted_rotmats(:,:,naccepted) = rotmat
                if( present(accepted_shifts) ) accepted_shifts(:,naccepted) = shift
                if( ratio > 0.75_dp ) mu = max(mu/2._dp,epsilon(1._dp))
                status = POSE_LM_ACCEPTED_IMPROVEMENT
                if( max(rotation_norm,shift_norm) < 1.e-8_dp .or. relative_reduction < 1.e-10_dp ) exit
            else
                mu = 4._dp*mu
            endif
        enddo
        if( status == POSE_LM_ITERATION_LIMIT .and. naccepted == 0 .and. bounded_trial ) &
            &status = POSE_LM_STEP_BOUND_REJECTED
    end subroutine refine_pose_lm

    !> Freeze inactive pose coordinates while retaining one five-vector LM path.
    pure subroutine apply_pose_parameter_mask( gradient, hessian, active )
        real(dp), intent(inout) :: gradient(5), hessian(5,5)
        logical, intent(in) :: active(5)
        integer :: axis

        do axis = 1, 5
            if( active(axis) ) cycle
            gradient(axis) = 0._dp
            hessian(axis,:) = 0._dp
            hessian(:,axis) = 0._dp
            hessian(axis,axis) = 1._dp
        enddo
    end subroutine apply_pose_parameter_mask

    !>  \brief  Cholesky solve with a relative pivot test for a 5-by-5
    !!          symmetric positive-definite pose block.
    pure subroutine solve_pose_cholesky( matrix, rhs, solution, reliable )
        real(dp), intent(in) :: matrix(5,5), rhs(5)
        real(dp), intent(out) :: solution(5)
        logical, intent(out) :: reliable
        real(dp) :: lower(5,5), intermediate(5), pivot, pivot_floor, matrix_scale
        integer :: i, j

        solution = 0._dp
        intermediate = 0._dp
        lower = 0._dp
        reliable = .false.
        if( any(.not. ieee_is_finite(matrix)) .or. any(.not. ieee_is_finite(rhs)) ) return
        matrix_scale = maxval(abs(matrix))
        if( matrix_scale <= POSE_NUMERIC_FLOOR ) return
        pivot_floor = sqrt(epsilon(1._dp))*matrix_scale
        do i = 1, 5
            do j = 1, i-1
                lower(i,j) = (matrix(i,j)-dot_product(lower(i,1:j-1),lower(j,1:j-1)))/lower(j,j)
            enddo
            pivot = matrix(i,i)-dot_product(lower(i,1:i-1),lower(i,1:i-1))
            if( .not. ieee_is_finite(pivot) .or. pivot <= pivot_floor ) return
            lower(i,i) = sqrt(pivot)
        enddo
        do i = 1, 5
            intermediate(i) = (rhs(i)-dot_product(lower(i,1:i-1),intermediate(1:i-1)))/lower(i,i)
        enddo
        do i = 5, 1, -1
            solution(i) = (intermediate(i)-dot_product(lower(i+1:5,i),solution(i+1:5)))/lower(i,i)
        enddo
        reliable = all(ieee_is_finite(solution))
    end subroutine solve_pose_cholesky

end module simple_cartesian_pose_refiner
