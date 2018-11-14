! parametric image registration using a biquadratic polynomial with L-BFGS-B in real-space (used in motion_correct)
module simple_motion_anisocor
include 'simple_lib.f08'
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_image,       only: image
implicit none
private
public :: motion_anisocor, POLY_DIM   !, test_motion_anisocor
#include "simple_local_flags.inc"

real,    parameter                  :: TOL      = 1e-4      !< tolerance parameter
integer, parameter                  :: MAXITS   = 30        !< maximum number of iterations
integer, parameter                  :: POLY_DIM = 12        !< dimensionality of polynomial dome model

type :: motion_anisocor
    private
    type(opt_spec)                  :: ospec                !< optimizer specification object
    class(optimizer),   pointer     :: nlopt      => null() !< pointer to nonlinear optimizer
    real                            :: maxHWshift = 0.      !< maximum half-width of shift
    class(image),       pointer     :: reference  => null() !< reference image ptr
    class(image),       pointer     :: frame      => null() !< particle image ptr
    integer                         :: ldim(2)              !< dimensions of reference, particle
    real(kind=c_float), pointer     :: rmat_ref   (:,:,:)   !< reference matrix
    real(kind=c_float), pointer     :: rmat       (:,:,:)   !< particle  matrix
!    real(kind=c_float), allocatable :: rmat_grad  (:,:,:)   !< gradient of particle matrix
    real(kind=c_float), pointer     :: rmat_T     (:,:,:)   !< transformed (interpolated) particle matrix
    real(kind=c_float), allocatable :: rmat_T_grad(:,:,:)   !< gradient of transformed particle matrix
    real                            :: motion_correctftol
    real                            :: motion_correctgtol
contains
    procedure, private              :: eval_fdf
    procedure, private              :: biquad_poly_df
    procedure, private              :: biquad_poly
    procedure, private              :: calc_mat_tld
    procedure, private              :: alloc_hlpmats
    ! procedure, private              :: interp_bilin
    procedure, private              :: interp_bilin_fdf
    procedure                       :: interp_bilin2
    procedure, private              :: interp_bilin2_fdf
    procedure, private              :: calc_gradient_cross_deriv
    procedure, private              :: bcuint
    procedure, private              :: bcucof
    procedure                       :: new      => motion_anisocor_new
    procedure                       :: kill     => motion_anisocor_kill
    procedure                       :: minimize => motion_anisocor_minimize
end type motion_anisocor

contains

    subroutine alloc_hlpmats( self )
        class(motion_anisocor), intent(inout) :: self
!!$        logical                               :: do_allocate
!!$        do_allocate = .false.
!!$        if (.not.(allocated(self%rmat_T_grad))) then
!!$            do_allocate = .true.
!!$        else
!!$            if ((size(self%rmat_T, 1) .ne. self%ldim(1)) .or. &
!!$                (size(self%rmat_T, 2) .ne. self%ldim(2))) do_allocate = .true.
!!$        end if
!!$        if (do_allocate) then
!!$            if (allocated(self%rmat_T))      deallocate(self%rmat_T)
!!$            if (allocated(self%rmat_T_grad)) deallocate(self%rmat_T_grad)
!!$            allocate(self%rmat_T     (self%ldim(1),self%ldim(2)),   &
!!$                self%rmat_T_grad(self%ldim(1),self%ldim(2),2)     )
!!$        end if
    end subroutine alloc_hlpmats

    ! ! bilinear interpolation
    ! function interp_bilin( self, x ) result( val )
    !     class(motion_anisocor), intent(inout) :: self
    !     real,                   intent(in)    :: x(2)
    !     real                                  :: val
    !     integer                               :: x1,x2,y1,y2
    !     real                                  :: NP1, NP2, NP3, NP4
    !     real                                  :: PW1, PW2, PW3, PW4
    !     ! Any values out of acceptable range
    !     if ((x(1) < 1.).or.(x(1) > real(self%ldim(1))).or.(x(2) < 1.).or.(x(2) > real(self%ldim(2)))) then
    !         val = 0.
    !         return
    !     end if
    !     x1 = min(floor(x(1)), self%ldim(1) - 1)
    !     x2 = x1 + 1
    !     y1 = min(floor(x(2)), self%ldim(2) - 1)
    !     y2 = y1 + 1
    !     ! Neighboring Pixels
    !     NP1 = self%rmat(y1,x1,1)
    !     NP2 = self%rmat(y1,x2,1)
    !     NP3 = self%rmat(y2,x1,1)
    !     NP4 = self%rmat(y2,x2,1)
    !     ! Pixel Weights
    !     PW1 = (y2-x(2))*(x2-x(1))
    !     PW2 = (y2-x(2))*(x(1)-x1)
    !     PW3 = (x2-x(1))*(x(2)-y1)
    !     PW4 = (x(2)-y1)*(x(1)-x1)
    !     val = PW1 * NP1 + PW2 * NP2 + PW3 * NP3 + PW4 * NP4
    ! end function interp_bilin
    !
    ! function interp_bilin2( self, x ) result( val )
    !     class(motion_anisocor), intent(inout) :: self
    !     real,                   intent(in)    :: x(2)
    !     real                                  :: val
    !     logical                               :: x1_valid, x2_valid, y1_valid, y2_valid
    !     integer                               :: x1_h, x2_h, y1_h, y2_h
    !     real                                  :: y1, y2, y3, y4, t, u
    !
    !     x1_h = floor(x(1))
    !     x2_h = x1_h + 1
    !     y1_h = floor(x(2))
    !     y2_h = y1_h + 1
    !
    !     x1_valid = ((x1_h > 0).and.(x1_h <= size(self%rmat, 1)))
    !     x2_valid = ((x2_h > 0).and.(x2_h <= size(self%rmat, 1)))
    !     y1_valid = ((y1_h > 0).and.(y1_h <= size(self%rmat, 2)))
    !     y2_valid = ((y2_h > 0).and.(y2_h <= size(self%rmat, 2)))
    !
    !     if (x1_valid .and. y1_valid) then
    !         y1 = self%rmat(x1_h, y1_h, 1)
    !     else
    !         y1 = 0.
    !     end if
    !
    !     if (x2_valid .and. y1_valid) then
    !         y2 = self%rmat(x2_h, y1_h, 1)
    !     else
    !         y2 = 0.
    !     end if
    !
    !     if (x2_valid .and. y2_valid) then
    !         y3 = self%rmat(x2_h, y2_h, 1)
    !     else
    !         y3 = 0.
    !     end if
    !
    !     if (x1_valid .and. y2_valid) then
    !         y4 = self%rmat(x1_h, y2_h, 1)
    !     else
    !         y4 = 0.
    !     end if
    !
    !
    !     t    = x(1) - real(x1_h);
    !     u    = x(2) - real(y1_h);
    !
    !     val  =  (1. - t) * (1. - u) * y1 + &
    !                   t  * (1. - u) * y2 + &
    !                   t  *       u  * y3 + &
    !             (1. - t) *       u  * y4
    ! end function interp_bilin2

    ! bilinear interpolation with gradient gradient
    subroutine interp_bilin_fdf( self, x, f, grad )
        class(motion_anisocor), intent(inout) :: self
        real,                   intent(in)    :: x(2)
        real,                   intent(out)   :: f
        real,                   intent(out)   :: grad(2)
        integer                               :: x1,x2,y1,y2
        real                                  :: NP1, NP2, NP3, NP4
        real                                  :: PW1, PW2, PW3, PW4
        real                                  :: PW1dx, PW2dx, PW3dx, PW4dx
        real                                  :: PW1dy, PW2dy, PW3dy, PW4dy
        ! Any values out of acceptable range
        if ((x(1) < 1.).or.(x(1) > real(self%ldim(1))).or.(x(2) < 1.).or.(x(2) > real(self%ldim(2)))) then
            f       = 0.
            grad(1) = 0.
            grad(2) = 0.
            return
        end if
        x1 = min(floor(x(1)), self%ldim(1) - 1)
        x2 = x1 + 1
        y1 = min(floor(x(2)), self%ldim(2) - 1)
        y2 = y1 + 1
        ! Neighboring Pixels
        NP1 = self%rmat(y1,x1,1)
        NP2 = self%rmat(y1,x2,1)
        NP3 = self%rmat(y2,x1,1)
        NP4 = self%rmat(y2,x2,1)
        ! Pixel Weights
        PW1 = (y2-x(2))*(x2-x(1))
        PW2 = (y2-x(2))*(x(1)-x1)
        PW3 = (x2-x(1))*(x(2)-y1)
        PW4 = (x(2)-y1)*(x(1)-x1)
        ! Pixel Weights, derivatives
        PW1dx = -(y2-x(2))
        PW2dx =  (y2-x(2))
        PW3dx = -(x(2)-y1)
        PW4dx =  (x(2)-y1)
        PW1dy = -(x2-x(1))
        PW2dy = -(x(1)-x1)
        PW3dy =  (x2-x(1))
        PW4dy =  (x(1)-x1)
        f       = PW1   * NP1 + PW2   * NP2 + PW3   * NP3 + PW4   * NP4
        grad(1) = PW1dx * NP1 + PW2dx * NP2 + PW3dx * NP3 + PW4dx * NP4
        grad(2) = PW1dy * NP1 + PW2dy * NP2 + PW3dy * NP3 + PW4dy * NP4
    end subroutine interp_bilin_fdf

    subroutine interp_bilin2_fdf( self, x, val, grad )
        class(motion_anisocor), intent(inout) :: self
        real,                   intent(in)    :: x(2)
        real,                   intent(out)   :: val, grad(2)
        logical                               :: x1_valid, x2_valid, y1_valid, y2_valid
        integer                               :: x1_h, x2_h, y1_h, y2_h
        real                                  :: y1, y2, y3, y4, t, u

        x1_h = floor(x(1))
        x2_h = x1_h + 1
        y1_h = floor(x(2))
        y2_h = y1_h + 1

        x1_valid = ((x1_h > 0) .and. (x1_h <= self%ldim(1)))
        x2_valid = ((x2_h > 0) .and. (x2_h <= self%ldim(1)))
        y1_valid = ((y1_h > 0) .and. (y1_h <= self%ldim(2)))
        y2_valid = ((y2_h > 0) .and. (y2_h <= self%ldim(2)))

        if (x1_valid .and. y1_valid) then
            y1 = self%rmat(x1_h, y1_h, 1)
        else
            y1 = 0.
        end if

        if (x2_valid .and. y1_valid) then
            y2 = self%rmat(x2_h, y1_h, 1)
        else
            y2 = 0.
        end if

        if (x2_valid .and. y2_valid) then
            y3 = self%rmat(x2_h, y2_h, 1)
        else
            y3 = 0.
        end if

        if (x1_valid .and. y2_valid) then
            y4 = self%rmat(x1_h, y2_h, 1)
        else
            y4 = 0.
        end if


        t    = x(1) - x1_h
        u    = x(2) - y1_h

        val  =  (1. - t) * (1. - u) * y1 + &
                      t  * (1. - u) * y2 + &
                      t  *       u  * y3 + &
                (1. - t) *       u  * y4

        grad(1) =   - (1. - u) * y1 + &
                      (1. - u) * y2 + &
                            u  * y3 - &
                            u  * y4
        grad(2) =   - (1. - t) * y1 - &
                            t  * y2 + &
                            t  * y3 + &
                      (1. - t) * y4
    end subroutine interp_bilin2_fdf

    subroutine interp_bilin2( self, x, rmat, val )
        class(motion_anisocor), intent(inout) :: self
        real,                   intent(in)    :: x(2), rmat(:,:)
        real,                   intent(out)   :: val
        logical                               :: x1_valid, x2_valid, y1_valid, y2_valid
        integer                               :: x1_h, x2_h, y1_h, y2_h
        real                                  :: y1, y2, y3, y4, t, u

        x1_h = floor(x(1))
        x2_h = x1_h + 1
        y1_h = floor(x(2))
        y2_h = y1_h + 1

        x1_valid = ((x1_h > 0) .and. (x1_h <= self%ldim(1)))
        x2_valid = ((x2_h > 0) .and. (x2_h <= self%ldim(1)))
        y1_valid = ((y1_h > 0) .and. (y1_h <= self%ldim(2)))
        y2_valid = ((y2_h > 0) .and. (y2_h <= self%ldim(2)))

        if (x1_valid .and. y1_valid) then
            y1 = rmat(x1_h, y1_h)
        else
            y1 = 0.
        end if

        if (x2_valid .and. y1_valid) then
            y2 = rmat(x2_h, y1_h)
        else
            y2 = 0.
        end if

        if (x2_valid .and. y2_valid) then
            y3 = rmat(x2_h, y2_h)
        else
            y3 = 0.
        end if

        if (x1_valid .and. y2_valid) then
            y4 = rmat(x1_h, y2_h)
        else
            y4 = 0.
        end if


        t    = x(1) - x1_h
        u    = x(2) - y1_h

        val  =  (1. - t) * (1. - u) * y1 + &
                      t  * (1. - u) * y2 + &
                      t  *       u  * y3 + &
                (1. - t) *       u  * y4

    end subroutine interp_bilin2

    subroutine bcuint(self, y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2, ansy, &
        &             ansy1, ansy2)
        ! TODO: check for cases of idx outside of range
        ! Bicubic interpolation within a grid square.  Input quantities are y,y1,y2,y12
        ! (as described in bcucof ); x1l and x1u,  the  lower  and  upper  coordinates
        ! of  the  grid  square  in  the  1-direction; x2l and x2u likewise  for  the  2-direction;
        ! and x1,x2,  the  coordinates  of  the desired  point  for  the  interpolation.  !
        ! The  interpolated  function  value  is  returned  as ansy, and the interpolated gradient values as
        ! ansy1 and ansy2.  This  routine calls bcucof.
        class(motion_anisocor), intent(inout) :: self
        real                                  :: ansy,ansy1,ansy2,x1,x1l,x1u,x2,x2l,x2u,y(4),y1(4), &
            &                                    y12(4),y2(4)
        integer                               :: i
        real                                  :: t,u,c(4,4)
        call self%bcucof(y, y1, y2, y12, x1u - x1l, x2u - x2l, c) ! Get the c's.
        if ((x1u == x1l) .or. (x2u == x2l)) THROW_HARD('bad input in bcuint')
        t     = (x1 - x1l) / (x1u - x1l)           ! Equation (3.6.4).
        u     = (x2 - x2l) / (x2u - x2l)
        ansy  = 0.
        ansy2 = 0.
        ansy1 = 0.
        do i = 4, 1, -1
            ! Equation (3.6.6).
            ansy  = t * ansy  + ((c(i,4)*u  +   c(i,3))*u+c(i,2)) * u + c(i,1)
            ansy2 = t * ansy2 +  (3.*c(i,4)*u+2.*c(i,3)         ) * u + c(i,2)
            ansy1 = u * ansy1 +  (3.*c(4,i)*t+2.*c(3,i)         ) * t + c(2,i)
        enddo
        ansy1 = ansy1 / (x1u - x1l)
        ansy2 = ansy2 / (x2u - x2l)
    end subroutine bcuint

    subroutine bcucof(self, y, y1, y2, y12, d1, d2, c)
        ! Given  arrays y,y1,y2, and y12, each of length 4, containing the function,
        ! gradients, and cross derivative at the four grid points of a rectangular
        ! grid cell (numbered counterclockwise from  the  lower  left),  and  given
        ! d1 and d2,  the  length  of  the  grid  cell  in  the  1-  and  2-directions,
        ! this  routine  returns  the  table c(1:4,1:4) that  is  used  by  routine
        ! bcuint for bicubic  interpolation.
        class(motion_anisocor), intent(inout) :: self
        real                                  :: d1, d2, c(4,4), y(4), y1(4), y12(4), y2(4)
        integer                               :: i, j, k, l
        real                                  :: d1d2, xx, cl(16), wt(16,16), x(16)
        SAVE wt
        DATA wt/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4 &
            ,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4             &
            ,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2    &
            ,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2             &
            ,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2        &
            ,10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2             &
            ,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1        &
            ,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
        d1d2 = d1 * d2
        do i = 1, 4 ! Pack a  temporary vector x.
            x(i)    = y(i)
            x(i+4)  = y1(i)  * d1
            x(i+8)  = y2(i)  * d2
            x(i+12) = y12(i) * d1d2
        enddo
        do i = 1, 16 ! Matrix  multiply  by  the  stored table.
            xx = 0.
            do k = 1,16
                xx = xx + wt(i,k) * x(k)
            enddo
            cl(i) = xx
        enddo
        l = 0
        do i = 1, 4 ! Unpack  the  result  into  the output  table.
            do j=1,4
                l      = l + 1
                c(i,j) = cl(l)
            enddo
        enddo
    end subroutine bcucof

    subroutine calc_gradient_cross_deriv( self )
        class(motion_anisocor), intent(inout) :: self


    end subroutine calc_gradient_cross_deriv

    subroutine calc_aniso_shifted( self, a, img )
        class(motion_anisocor), intent(inout) :: self
        real,                   intent(in)    :: a(POLY_DIM)
        class(image),           intent(inout) :: img
        real, allocatable :: rmat_tmp(:,:)
        real, pointer     :: rmat_ptr(:,:,:)
        integer           :: i, j, ldim(3)
        real              :: x_tld(2), grad(2), val
        ldim = img%get_ldim()
        if( .not. all(ldim(1:2) .eq. self%ldim(1:2)) )then
            THROW_HARD('in calc_aniso_shifted: dimension do not match; simple_motion_anisocor')
        endif
        allocate(rmat_tmp(self%ldim(1),self%ldim(2)))
        call img%get_rmat_ptr(rmat_ptr)
        rmat_tmp(1:self%ldim(1),1:self%ldim(2)) = rmat_ptr(1:self%ldim(1),1:self%ldim(2),1)
        do i = 1, self%ldim(1)
            do j = 1, self%ldim(2)
                call self%biquad_poly(i, j, a, x_tld)
                call self%interp_bilin2_fdf(x_tld, val, grad)
                rmat_tmp(i,j) = val
            end do
        end do
        deallocate(rmat_tmp)
    end

    subroutine calc_mat_tld( self, a )
        class(motion_anisocor), intent(inout) :: self
        real,                   intent(in)    :: a(POLY_DIM)
        integer                               :: i, j
        real                                  :: x_tld(2)
        real                                  :: val, grad(2)
        do i = 1, self%ldim(1)
            do j = 1, self%ldim(2)
                call self%biquad_poly(i, j, a, x_tld)
                call self%interp_bilin2_fdf(x_tld, val, grad)
                self%rmat_T     (i,j,1) = val
                self%rmat_T_grad(i,j,1) = grad(1)
                self%rmat_T_grad(i,j,2) = grad(2)
            end do
        end do
    end subroutine calc_mat_tld

    subroutine biquad_poly( self, i, j, a, x_tld)
        class(motion_anisocor), intent(inout) :: self
        integer,                intent(in)    :: i, j
        real,                   intent(in)    :: a(POLY_DIM)
        real,                   intent(out)   :: x_tld(2)
        real                                  :: x, y, xm, ym, xms, yms, z1, z2
        real                                  :: Nx, Ny
        Nx = self%ldim(1)
        Ny = self%ldim(2)
        x  = real(i)
        y  = real(j)
        ! map x,y to [-1:1]
        xm = (2. * x - Nx - 1.) / (Nx - 1.)
        ym = (2. * y - Ny - 1.) / (Ny - 1.)
        xms = xm * xm
        yms = ym * ym
        ! Polynomial:
        ! Tx = a1 + a2 * x + a3 * y + a4  * x^2 + a5  * y^2 + a6  * x*y
        ! Ty = a7 + a8 * x + a9 * y + a10 * x^2 + a11 * y^2 + a12 * x*y
        z1 = real(xm + a(1) + a(2) * xm + a(3) * ym + a( 4) * xms + a( 5) * yms + a( 6) * xm*ym)
        z2 = real(ym + a(7) + a(8) * xm + a(9) * ym + a(10) * xms + a(11) * yms + a(12) * xm*ym)
        ! remap to [1:N]
        x_tld(1) = ((Nx - 1.) * z1 + Nx + 1.) / 2.
        x_tld(2) = ((Nx - 1.) * z2 + Ny + 1.) / 2.
    end subroutine biquad_poly

    subroutine biquad_poly_df( self, i, j, a, poly_grad )
        class(motion_anisocor), intent(inout) :: self
        integer,                intent(in)    :: i, j
        real,                   intent(in)    :: a(POLY_DIM)
        real                                  :: poly_grad(POLY_DIM)
        real                                  :: x, y, xm, ym, xms, yms
        real                                  :: Nx, Ny
        Nx  = self%ldim(1)
        Ny  = self%ldim(2)
        x   = real(i)
        y   = real(j)
        ! map x,y to [-1:1]
        xm  = (2. * x - Nx - 1.) / (Nx - 1.)
        ym  = (2. * y - Ny - 1.) / (Ny - 1.)
        xms = xm * xm
        yms = ym * ym
        ! Polynomial:
        ! Tx = a1 + a2 * x + a3 * y + a4  * x^2 + a5  * y^2 + a6  * x*y
        ! Ty = a7 + a8 * x + a9 * y + a10 * x^2 + a11 * y^2 + a12 * x*y
        poly_grad( 1)    = 1.
        poly_grad( 2)    = xm
        poly_grad( 3)    = ym
        poly_grad( 4)    = xms
        poly_grad( 5)    = yms
        poly_grad( 6)    = xm*ym
        poly_grad( 7)    = 1.
        poly_grad( 8)    = xm
        poly_grad( 9)    = ym
        poly_grad(10)    = xms
        poly_grad(11)    = yms
        poly_grad(12)    = xm*ym
        poly_grad( 1: 6) = poly_grad( 1: 6) * (Nx - 1.) / 2.
        poly_grad( 7:12) = poly_grad( 7:12) * (Ny - 1.) / 2.
    end subroutine biquad_poly_df

    subroutine eval_fdf( self, a, f, grad )
        class(motion_anisocor), intent(inout) :: self
        real,                   intent(in)    :: a(POLY_DIM)
        real                                  :: f, grad(POLY_DIM)
        real                                  :: N, D_I, D_R, grad_tmp1, grad_tmp2
        integer                               :: i, j, r
        real                                  :: poly_dfs(POLY_DIM)
        call self%calc_mat_tld(a)
        N   = 0.
        D_I = 0.
        D_R = 0.
        do i = 1, self%ldim(1)
            do j = 1, self%ldim(2)
                N   = N   + self%rmat_ref(i,j,1) * self%rmat_T  (i,j,1)
                D_I = D_I + self%rmat_ref(i,j,1) * self%rmat_ref(i,j,1)
                D_R = D_R + self%rmat_T  (i,j,1) * self%rmat_T  (i,j,1)
            end do
        end do
        f = N / sqrt(D_I * D_R)
        ! Computing gradient
        grad = 0.
        do r = 1,POLY_DIM
            grad_tmp1 = 0.
            grad_tmp2 = 0.
            do i = 1,self%ldim(1)
                do j = 1,self%ldim(2)
                    call self%biquad_poly_df( i, j, a, poly_dfs )
                    if ( r <= 6 ) then
                        grad_tmp1 = grad_tmp1 + self%rmat_ref(i,j,1) * poly_dfs(r) * self%rmat_T_grad(i,j,1)
                        grad_tmp2 = grad_tmp2 + self%rmat_T  (i,j,1) * poly_dfs(r) * self%rmat_T_grad(i,j,1)
                    else
                        grad_tmp1 = grad_tmp1 + self%rmat_ref(i,j,1) * poly_dfs(r) * self%rmat_T_grad(i,j,2)
                        grad_tmp2 = grad_tmp2 + self%rmat_T  (i,j,1) * poly_dfs(r) * self%rmat_T_grad(i,j,2)
                    end if
                end do
            end do
            grad(r) = 1. / sqrt(D_I*D_R) * (grad_tmp1 - f * sqrt(D_I / D_R) * grad_tmp2)
        end do
    end subroutine eval_fdf

    !> Initialise  ftexp_shsrch
    subroutine motion_anisocor_new( self, ref, frame, aniso_shifted_frame, motion_correct_ftol, motion_correct_gtol )
        class(motion_anisocor), intent(inout) :: self
        class(image), target, intent(in) :: ref, frame, aniso_shifted_frame
        real,       optional, intent(in) :: motion_correct_ftol, motion_correct_gtol
        type(opt_factory)                :: ofac                 !< optimizer factory
        integer                          :: i
        real                             :: lims(POLY_DIM,2)
        integer                          :: ldim_here1(3), ldim_here2(3)
        call self%kill()
        ldim_here1 = ref  %get_ldim()
        ldim_here2 = frame%get_ldim()
        if ((ldim_here1(1) /= ldim_here2(1)) .or. (ldim_here1(2) /= ldim_here2(2))) then
            THROW_HARD('motion_anisocor_new: dimensions of particle and reference do not match; simple_motion_anisocor')
        end if
        self%reference  => ref
        call ref%get_rmat_ptr(self%rmat_ref)
        self%frame      => frame
        call frame%get_rmat_ptr(self%rmat)
        call aniso_shifted_frame%get_rmat_ptr(self%rmat_T)
        self%ldim       =  ref%get_ldim()
        allocate( self%rmat_T_grad ( self%ldim(1), self%ldim(2), 2 ) )
        self%maxHWshift = 0.1
        if( present(motion_correct_ftol) )then
            self%motion_correctftol = motion_correct_ftol
        else
            self%motion_correctftol = TOL
        end if
        if( present(motion_correct_gtol) )then
            self%motion_correctgtol = motion_correct_gtol
        else
            self%motion_correctgtol = TOL
        end if
        do i = 1,POLY_DIM
            lims(i,1) = - self%maxHWshift
            lims(1,2) =   self%maxHWshift
        end do
        call self%ospec%specify('lbfgsb', POLY_DIM, ftol=self%motion_correctftol, gtol=self%motion_correctgtol, limits=lims)
        call self%ospec%set_costfun_8   (motion_anisocor_cost_8)
        call self%ospec%set_gcostfun_8  (motion_anisocor_gcost_8)
        call self%ospec%set_fdfcostfun_8(motion_anisocor_fdfcost_8)
        ! generate optimizer object with the factory
        if( associated(self%nlopt) )then
            call self%nlopt%kill()
            deallocate(self%nlopt)
        end if
        call ofac%new(self%ospec, self%nlopt)
    end subroutine motion_anisocor_new

    !> Cost function, double precision
    function motion_anisocor_cost_8( fun_self, vec, D ) result(cost)
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(in)    :: vec(D)
        real(kind=8)                :: cost
        real                        :: f, grad(POLY_DIM)
        select type(fun_self)
        class is (motion_anisocor)
            call fun_self%eval_fdf(real(vec), f, grad)
            cost = -real(f, kind=dp)
        class default
            THROW_HARD('error in motion_anisocor_cost_8: unknown type; simple_motion_anisocor')
        end select
    end function motion_anisocor_cost_8

    !> Gradient function, double precision
    subroutine motion_anisocor_gcost_8( fun_self, vec, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(dp),     intent(inout) :: vec( D )
        real(dp),     intent(out)   :: grad( D )
        real                        :: f, spgrad(POLY_DIM)
        select type(fun_self)
        class is (motion_anisocor)
            call fun_self%eval_fdf(real(vec), f, spgrad)
            grad = -real(spgrad, kind=dp)
        class default
            THROW_HARD('error in motion_anisocor_gcost_8: unknown type; simple_motion_anisocor')
        end select
    end subroutine motion_anisocor_gcost_8

    !> Gradient & cost function, double precision
    subroutine motion_anisocor_fdfcost_8( fun_self, vec, f, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
        real                        :: spf, spgrad(POLY_DIM)
        select type(fun_self)
        class is (motion_anisocor)
            call fun_self%eval_fdf(real(vec), spf, spgrad)
            f    = -real(f,      kind=dp)
            grad = -real(spgrad, kind=dp)
        class default
            THROW_HARD('error in motion_anisocor_fdfcost_8: unknown type; simple_motion_anisocor')
        end select
    end subroutine motion_anisocor_fdfcost_8

    !> Main search routine
    function motion_anisocor_minimize( self ) result( cxy )
        class(motion_anisocor), intent(inout) :: self
        real :: cxy(POLY_DIM+1) ! corr + model
        class(*), pointer :: fun_self => null()
        integer           :: i
        self%ospec%x = 0.
        call self%nlopt%minimize(self%ospec, self, cxy(1))
        cxy(2:) = self%ospec%x


!!$        lims(1,1) = - maxHWshift
!!$        lims(1,2) =   maxHWshift
!!$        lims(2,1) = - maxHWshift
!!$        lims(2,2) =   maxHWshift
!!$        if( present(prev_shift) )then
!!$            lims(1,:) = lims(1,:) + prev_shift(1)
!!$            lims(2,:) = lims(2,:) + prev_shift(2)
!!$        endif
!!$        call ospec%specify('lbfgsb', 2, ftol=motion_correctftol, gtol=motion_correctgtol, limits=lims) ! nrestarts=nrestarts
!!$        call ospec%set_costfun_8(ftexp_shsrch_cost_8)
!!$        call ospec%set_gcostfun_8(ftexp_shsrch_gcost_8)
!!$        call ospec%set_fdfcostfun_8(ftexp_shsrch_fdfcost_8)
!!$        ! generate optimizer object with the factory
!!$        if( associated(nlopt) )then
!!$            call nlopt%kill
!!$            deallocate(nlopt)
!!$        end if
!!$        call ofac%new(ospec, nlopt)
!!$        if( present(prev_shift) )ospec%x = prev_shift
!!$        ! set initial solution to previous shift
!!$        call nlopt%minimize(ospec, fun_self, cxy(1))
!!$        call reference%corr_normalize(particle, cxy(1))
!!$        call ft_exp_reset_tmp_pointers
!!$        cxy(1)  = -cxy(1) ! correlation
!!$        cxy(2:) = ospec%x ! shift
!!$        if( present(prev_corr) )then
!!$            if( abs(cxy(1)-prev_corr) <= TOL )then
!!$                cxy(1)  = prev_corr
!!$                if( present(prev_shift) ) cxy(2:) = prev_shift
!!$            endif
!!$        endif
    end function motion_anisocor_minimize

    subroutine motion_anisocor_kill( self )
        class(motion_anisocor), intent(inout) :: self
!        if (allocated(self%rmat_grad)) deallocate(self%rmat_grad)
        if (allocated(self%rmat_T_grad)) deallocate(self%rmat_T_grad)
        if (associated(self%nlopt)) then
            call self%nlopt%kill()
            self%nlopt => null()
        end if
        call self%ospec%kill()
        self%maxHWshift  = 0.
        self%reference   => null()
        self%frame       => null()
        self%ldim        = 0
    end subroutine motion_anisocor_kill

!!$    subroutine test_motion_anisocor
!!$        type(image)       :: img_ref, img_ptcl
!!$        type(ft_expanded) :: ftexp_ref, ftexp_ptcl
!!$        real, parameter   :: TRS=5.0, lp=6., hp=100.
!!$        real              :: cxy(3), x, y, lims(2,2)
!!$        integer           :: i
!!$        lims(:,1) = -TRS
!!$        lims(:,2) = TRS
!!$        call img_ref%new([100,100,1], 2.)
!!$        call img_ptcl%new([100,100,1], 2.)
!!$        call img_ref%square(20)
!!$        call img_ref%fft()
!!$        call ftexp_ref%new(img_ref, hp, lp)
!!$        call ftexp_ptcl%new(img_ptcl, hp, lp)
!!$        call ftexp_shsrch_init(ftexp_ref, ftexp_ptcl, trs)
!!$        do i=1,100
!!$            x = ran3()*2*TRS-TRS
!!$            y = ran3()*2*TRS-TRS
!!$            img_ptcl = img_ref
!!$            call img_ptcl%shift([x,y,0.])
!!$            call ftexp_ptcl%new(img_ptcl, hp, lp)
!!$            cxy = ftexp_shsrch_minimize()
!!$            if( cxy(1) < 0.995 )then
!!$                THROW_HARD('shift alignment failed; test_ftexp_shsrch')
!!$            endif
!!$        end do
!!$        write(*,'(a)') 'SIMPLE_ftexp_shsrch_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
!!$    end subroutine test_motion_anisocor

end module simple_motion_anisocor
