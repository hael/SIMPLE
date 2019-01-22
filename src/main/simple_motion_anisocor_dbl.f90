! parametric image registration using a biquadratic polynomial with L-BFGS-B in real-space (used in motion_correct), double precision
module simple_motion_anisocor_dbl
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_opt_factory,     only: opt_factory
use simple_opt_spec,        only: opt_spec
use simple_optimizer,       only: optimizer
use simple_image,           only: image
use simple_motion_anisocor, only: motion_anisocor, POLY_DIM, TOL
implicit none
private
public :: motion_anisocor_dbl
#include "simple_local_flags.inc"

type, extends(motion_anisocor)       :: motion_anisocor_dbl
    private
    real(kind=c_double), allocatable :: rmat_ref      (:,:)   !< reference matrix
    real(kind=c_double), allocatable :: rmat          (:,:)   !< particle  matrix
    real(kind=c_double), allocatable :: rmat_T        (:,:)   !< transformed (interpolated) particle matrix
    real(kind=c_double), allocatable :: rmat_T_grad   (:,:,:) !< gradient of transformed particle matrix
    real(kind=c_double), allocatable :: T_coords      (:,:,:) !< transformed coordinates
    real(kind=c_double), allocatable :: T_coords_xm   (:,:,:) !< transformed coordinates, for use in derivative of obj. function
    real(kind=c_double), allocatable :: T_coords_xm_sq(:,:,:) !< square of T_coords_xm, memoized for speedup
    real(kind=c_double), allocatable :: T_coords_xm_cr(:,:)   !< cross-product of square of T_coords_xm, memoized for speedup
    real(kind=c_double), allocatable :: T_coords_out  (:,:,:) !< transformed coordinates for output image
    real(dp)                         :: sc_alpha_x, sc_alpha_y!< scaling factors, used for gradients
    real(dp)                         :: sc_beta_x, sc_beta_y
    real(dp)                         :: sc_gamma_x, sc_gamma_y
    real(dp)                         :: gamma
contains
    procedure                        :: eval_fdf_foo     ! for debug purpopses
    procedure, private               :: eval_fdf
    procedure, private               :: eval_f_noregu    ! for debug purposes
    procedure, private               :: calc_mat_tld
    procedure, private               :: interp_bilin
    procedure, private               :: interp_bilin_fdf
    procedure, private               :: interp_bilin_out
    procedure, private               :: calc_T_coords
    procedure, private               :: calc_grad_coords    ! specifically for gradient of regularization term
    procedure, private               :: calc_T_coords_only
    procedure, private               :: calc_T_coords_out
    procedure                        :: calc_T_out => motion_anisocor_dbl_calc_T_out
    procedure                        :: new        => motion_anisocor_dbl_new
    procedure                        :: kill       => motion_anisocor_dbl_kill
    procedure                        :: minimize   => motion_anisocor_dbl_minimize
end type motion_anisocor_dbl

contains

    function interp_bilin( self, x ) result(val)
        class(motion_anisocor_dbl), intent(inout) :: self
        real(dp),                   intent(in)    :: x(2)
        real                                  :: val
        logical  :: x1_valid, x2_valid, y1_valid, y2_valid
        integer  :: x1_h, x2_h, y1_h, y2_h
        real(dp) :: y1, y2, y3, y4, t, u
        ! if outside of image
        if ((x(1) < 1._dp) .or. (x(1) >= real(self%ldim(1),dp)) .or. &
            (x(2) < 1._dp) .or. (x(2) >= real(self%ldim(2),dp))) then
            val  = 0._dp
            return
        end if
        x1_h = floor(x(1))
        x2_h = x1_h + 1
        y1_h = floor(x(2))
        y2_h = y1_h + 1
        y1 = self%rmat(x1_h, y1_h)
        y2 = self%rmat(x2_h, y1_h)
        y3 = self%rmat(x2_h, y2_h)
        y4 = self%rmat(x1_h, y2_h)
        t    = x(1) - x1_h
        u    = x(2) - y1_h
        val  =  (1._dp - t) * (1._dp - u) * y1 + &
                         t  * (1._dp - u) * y2 + &
                         t  *          u  * y3 + &
                (1._dp - t) *          u  * y4
    end function interp_bilin

    subroutine interp_bilin_fdf( self, x, val, grad )
        class(motion_anisocor_dbl), intent(inout) :: self
        real(dp),                   intent(in)    :: x(2)
        real(dp),                   intent(out)   :: val, grad(2)
        logical  :: x1_valid, x2_valid, y1_valid, y2_valid
        integer  :: x1_h, x2_h, y1_h, y2_h
        real(dp) :: y1, y2, y3, y4, t, u
        ! if outside of image
        if ((x(1) < 1._dp) .or. (x(1) >= real(self%ldim(1),dp)) .or. &
            (x(2) < 1._dp) .or. (x(2) >= real(self%ldim(2),dp))) then
            val  = 0._dp
            grad = 0._dp
            return
        end if
        x1_h = floor(x(1))
        x2_h = x1_h + 1
        y1_h = floor(x(2))
        y2_h = y1_h + 1
        y1 = self%rmat(x1_h, y1_h)
        y2 = self%rmat(x2_h, y1_h)
        y3 = self%rmat(x2_h, y2_h)
        y4 = self%rmat(x1_h, y2_h)
        t    = x(1) - x1_h
        u    = x(2) - y1_h
        val     =   (1._dp - t) * (1._dp - u) * y1 + &
                             t  * (1._dp - u) * y2 + &
                             t  *          u  * y3 + &
                    (1._dp - t) *          u  * y4
        grad(1) = -               (1._dp - u) * y1 + &
                                  (1._dp - u) * y2 + &
                                           u  * y3 - &
                                           u  * y4
        grad(2) = - (1._dp - t) *               y1 - &
                             t  *               y2 + &
                             t  *               y3 + &
                    (1._dp - t) *               y4
    end subroutine interp_bilin_fdf

    function interp_bilin_out( self, x, rmat_in ) result(val)
        class(motion_anisocor_dbl), intent(inout) :: self
        real(dp),                   intent(in)    :: x(2)
        real, pointer,              intent(in)    :: rmat_in(:,:,:)
        real(dp) :: val
        logical  :: x1_valid, x2_valid, y1_valid, y2_valid
        integer  :: x1_h, x2_h, y1_h, y2_h
        real(dp) :: y1, y2, y3, y4, t, u
        ! if outside of image
        if ((x(1) < 1._dp) .or. (x(1) >= self%ldim_out(1)) .or. (x(2) < 1._dp) .or. (x(2) >= self%ldim_out(2))) then
            val  = 0._dp
            return
        end if
        x1_h = floor(x(1))
        x2_h = x1_h + 1
        y1_h = floor(x(2))
        y2_h = y1_h + 1
        y1 = rmat_in(x1_h, y1_h, 1)
        y2 = rmat_in(x2_h, y1_h, 1)
        y3 = rmat_in(x2_h, y2_h, 1)
        y4 = rmat_in(x1_h, y2_h, 1)
        t    = x(1) - x1_h
        u    = x(2) - y1_h
        val  =  (1._dp - t) * (1._dp - u) * y1 + &
                         t  * (1._dp - u) * y2 + &
                         t  *          u  * y3 + &
                (1._dp - t) *          u  * y4
    end function interp_bilin_out

    subroutine calc_T_coords( self, a )
        class(motion_anisocor_dbl), intent(inout) :: self
        real(dp),                   intent(in)    :: a(POLY_DIM)
        integer  :: i, j
        real(dp) :: Nx, Ny
        real(dp) :: x, y, xm, ym, xms, yms, xym, z1, z2
        Nx = real(self%ldim(1), dp)
        Ny = real(self%ldim(2), dp)
        do j = 1, self%ldim(2)
            do i = 1, self%ldim(1)
                x   = real(i, dp)
                y   = real(j, dp)
                ! map x,y to [-1:1]
                xm  = (2._dp * x - Nx - 1._dp) / (Nx - 1._dp)
                ym  = (2._dp * y - Ny - 1._dp) / (Ny - 1._dp)
                xms = xm * xm
                yms = ym * ym
                xym = xm * ym
                self%T_coords_xm   (i,j,1) = xm
                self%T_coords_xm   (i,j,2) = ym
                self%T_coords_xm_sq(i,j,1) = xms
                self%T_coords_xm_sq(i,j,2) = yms
                self%T_coords_xm_cr(i,j)   = xym
                ! Polynomial:
                ! Tx = a1 + a2 * x + a3 * y + a4  * x^2 + a5  * y^2 + a6  * x*y
                ! Ty = a7 + a8 * x + a9 * y + a10 * x^2 + a11 * y^2 + a12 * x*y
                z1 = xm + a(1) + a(2) * xm + a(3) * ym + a( 4) * xms + a( 5) * yms + a( 6) * xym
                z2 = ym + a(7) + a(8) * xm + a(9) * ym + a(10) * xms + a(11) * yms + a(12) * xym
                ! remap to [1:N]
                self%T_coords      (1,i,j) = ((Nx - 1._dp) * z1 + Nx + 1._dp) / 2._dp
                self%T_coords      (2,i,j) = ((Ny - 1._dp) * z2 + Ny + 1._dp) / 2._dp
            end do
        end do
    end subroutine calc_T_coords

    subroutine calc_grad_coords( self, a, grad )                 ! specifically for gradient of regularization term
        class(motion_anisocor_dbl), intent(inout) :: self
        real(dp),                   intent(in)    :: a(POLY_DIM)
        real(dp),                   intent(out)   :: grad(POLY_DIM)
        integer  :: i, j
        real(dp) :: Nx, Ny
        real(dp) :: x, y, xm, ym, xms, yms, xym, z1, z2
        real(dp) :: T_coords(2)
        grad = 0.
        Nx = real(self%ldim(1), dp)
        Ny = real(self%ldim(2), dp)
        do j = 1, self%ldim(2)
            do i = 1, self%ldim(1)
                x   = real(i, dp)
                y   = real(j, dp)
                ! map x,y to [-1:1]
                xm  = (2._dp * x - Nx - 1._dp) / (Nx - 1._dp)
                ym  = (2._dp * y - Ny - 1._dp) / (Ny - 1._dp)
                xms = xm * xm
                yms = ym * ym
                xym = xm * ym
                self%T_coords_xm   (i,j,1) = xm
                self%T_coords_xm   (i,j,2) = ym
                self%T_coords_xm_sq(i,j,1) = xms
                self%T_coords_xm_sq(i,j,2) = yms
                self%T_coords_xm_cr(i,j)   = xym
                ! Polynomial:
                ! Tx = a1 + a2 * x + a3 * y + a4  * x^2 + a5  * y^2 + a6  * x*y
                ! Ty = a7 + a8 * x + a9 * y + a10 * x^2 + a11 * y^2 + a12 * x*y
                z1 = xm + a(1) + a(2) * xm + a(3) * ym + a( 4) * xms + a( 5) * yms + a( 6) * xym
                z2 = ym + a(7) + a(8) * xm + a(9) * ym + a(10) * xms + a(11) * yms + a(12) * xym
                ! remap to [1:N]
                T_coords(1) = ((Nx - 1._dp) * z1 + Nx + 1._dp) / 2._dp
                T_coords(2) = ((Ny - 1._dp) * z2 + Ny + 1._dp) / 2._dp
                grad(1) = grad(1) - 2. * (x - T_coords(1)) * (Nx - 1._dp) / 2._dp
                grad(2) = grad(2) - 2. * (x - T_coords(1)) * (Nx - 1._dp) / 2._dp * xm
                grad(3) = grad(3) - 2. * (x - T_coords(1)) * (Nx - 1._dp) / 2._dp * ym
                grad(4) = grad(4) - 2. * (x - T_coords(1)) * (Nx - 1._dp) / 2._dp * xms
                grad(5) = grad(5) - 2. * (x - T_coords(1)) * (Nx - 1._dp) / 2._dp * yms
                grad(6) = grad(6) - 2. * (x - T_coords(1)) * (Nx - 1._dp) / 2._dp * xym
                grad(7) = grad(7) - 2. * (y - T_coords(2)) * (Ny - 1._dp) / 2._dp 
                grad(8) = grad(8) - 2. * (y - T_coords(2)) * (Ny - 1._dp) / 2._dp * xm
                grad(9) = grad(9) - 2. * (y - T_coords(2)) * (Ny - 1._dp) / 2._dp * ym
                grad(10) = grad(10) - 2. * (y - T_coords(2)) * (Ny - 1._dp) / 2._dp * xms
                grad(11) = grad(11) - 2. * (y - T_coords(2)) * (Ny - 1._dp) / 2._dp * yms
                grad(12) = grad(12) - 2. * (y - T_coords(2)) * (Ny - 1._dp) / 2._dp * xym                          
            end do
        end do
    end subroutine calc_grad_coords                          ! specifically for gradient of regularization term
 
    subroutine calc_T_coords_only( self, a )
        class(motion_anisocor_dbl), intent(inout) :: self
        real(dp),                   intent(in)    :: a(POLY_DIM)
        integer  :: i, j
        real(dp) :: Nx, Ny
        real(dp) :: x, y, xm, ym, xms, yms, z1, z2
        Nx = real(self%ldim(1), dp)
        Ny = real(self%ldim(2), dp)
        do j = 1, self%ldim(2)
            do i = 1, self%ldim(1)
                x   = real(i, dp)
                y   = real(j, dp)
                ! map x,y to [-1:1]
                xm  = (2._dp * x - Nx - 1._dp) / (Nx - 1._dp)
                ym  = (2._dp * y - Ny - 1._dp) / (Ny - 1._dp)
                xms = xm * xm
                yms = ym * ym
                ! Polynomial:
                ! Tx = a1 + a2 * x + a3 * y + a4  * x^2 + a5  * y^2 + a6  * x*y
                ! Ty = a7 + a8 * x + a9 * y + a10 * x^2 + a11 * y^2 + a12 * x*y
                z1 = xm + a(1) + a(2) * xm + a(3) * ym + a( 4) * xms + a( 5) * yms + a( 6) * xm*ym
                z2 = ym + a(7) + a(8) * xm + a(9) * ym + a(10) * xms + a(11) * yms + a(12) * xm*ym
                ! remap to [1:N]
                self%T_coords(1,i,j) = ((Nx - 1._dp) * z1 + Nx + 1._dp) / 2._dp
                self%T_coords(2,i,j) = ((Ny - 1._dp) * z2 + Ny + 1._dp) / 2._dp
            end do
        end do
    end subroutine calc_T_coords_only

    subroutine calc_T_coords_out( self, a )
        class(motion_anisocor_dbl), intent(inout) :: self
        real(dp),                   intent(in)    :: a(POLY_DIM)
        integer  :: i, j
        real(dp) :: Nx, Ny
        real(dp) :: x, y, xm, ym, xms, yms, z1, z2
        Nx = real(self%ldim_out(1), dp)
        Ny = real(self%ldim_out(2), dp)
        do j = 1, self%ldim_out(2)
            do i = 1, self%ldim_out(1)
                x   = real(i, dp)
                y   = real(j, dp)
                ! map x,y to [-1:1]
                xm  = (2._dp * x - Nx - 1._dp) / (Nx - 1._dp)
                ym  = (2._dp * y - Ny - 1._dp) / (Ny - 1._dp)
                xms = xm * xm
                yms = ym * ym
                ! Polynomial:
                ! Tx = a1 + a2 * x + a3 * y + a4  * x^2 + a5  * y^2 + a6  * x*y
                ! Ty = a7 + a8 * x + a9 * y + a10 * x^2 + a11 * y^2 + a12 * x*y
                z1 = xm + a(1) + a(2) * xm + a(3) * ym + a( 4) * xms + a( 5) * yms + a( 6) * xm*ym
                z2 = ym + a(7) + a(8) * xm + a(9) * ym + a(10) * xms + a(11) * yms + a(12) * xm*ym
                ! remap to [1:N]
                self%T_coords_out(1,i,j) = ((Nx - 1._dp) * z1 + Nx + 1._dp) / 2._dp
                self%T_coords_out(2,i,j) = ((Ny - 1._dp) * z2 + Ny + 1._dp) / 2._dp
            end do
        end do
    end subroutine calc_T_coords_out

    subroutine motion_anisocor_dbl_calc_T_out( self, a, frame_in, frame_out )
        class(motion_anisocor_dbl), intent(inout) :: self
        real(dp),                   intent(in)    :: a(12)
        class(image), target,       intent(in)    :: frame_in
        class(image), target,       intent(inout) :: frame_out
        integer       :: i,j
        real(dp)      :: x_tld(2)
        integer       :: ldim_out_here1(3), ldim_out_here2(3)
        logical       :: do_alloc
        real, pointer :: rmat_in_ptr(:,:,:), rmat_out_ptr(:,:,:)
        ldim_out_here1 = frame_in%get_ldim()
        ldim_out_here2 = frame_out%get_ldim()
        if( any(ldim_out_here1(1:2) .ne. ldim_out_here2(1:2)) )then
            THROW_HARD('calc_T_out: dimensions of particle and reference do not match; simple_motion_anisocor_dbl')
        end if
        self%ldim_out(1) = ldim_out_here1(1)
        self%ldim_out(2) = ldim_out_here1(2)
        call frame_in %get_rmat_ptr( rmat_in_ptr  )
        call frame_out%get_rmat_ptr( rmat_out_ptr )
        do_alloc = .false.
        if (.not. allocated(self%T_coords_out)) then
            do_alloc = .true.
        else
            if ((size(self%T_coords_out, 2) .ne. self%ldim_out(1)).or.(size(self%T_coords_out, 3) .ne. self%ldim_out(2))) then
                do_alloc = .true.
            end if
        end if
        if (do_alloc) then
            if (allocated(self%T_coords_out)) deallocate(self%T_coords_out)
            allocate(self%T_coords_out(2, self%ldim_out(1), self%ldim_out(2)))
        end if
        call self%calc_T_coords_out(a)
        do j = 1, self%ldim_out(2)
            do i = 1, self%ldim_out(1)
                x_tld = self%T_coords_out(:,i,j)
                rmat_out_ptr(i,j,1) = self%interp_bilin_out(x_tld, rmat_in_ptr)
            end do
        end do
    end subroutine motion_anisocor_dbl_calc_T_out

    subroutine calc_mat_tld( self )
        class(motion_anisocor_dbl), intent(inout) :: self
        integer  :: i, j
        real(dp) :: x_tld(2)
        real(dp) :: val, grad(2)
        do j = 1, self%ldim(2)
            do i = 1, self%ldim(1)
                x_tld = self%T_coords(:,i,j)
                call self%interp_bilin_fdf(x_tld, val, grad)
                self%rmat_T     (i,j)   = val
                self%rmat_T_grad(i,j,1) = grad(1)
                self%rmat_T_grad(i,j,2) = grad(2)
            end do
        end do
    end subroutine calc_mat_tld

    subroutine eval_fdf_foo( self, ref, frame, a, f, grad )  ! for debug purpopse
        class(motion_anisocor_dbl), intent(inout) :: self
        class(image), target,       intent(in)    :: ref, frame
        real(dp),                   intent(in)    :: a(POLY_DIM)
        real(dp),                   intent(out)   :: f, grad(POLY_DIM)
        class(*), pointer :: fun_self => null()
        integer           :: i, j
        real, pointer     :: rmat_ref_ptr(:,:,:), rmat_ptr(:,:,:)
        integer           :: ldim_here1(3), ldim_here2(3)
        logical           :: do_alloc
        ldim_here1 = ref  %get_ldim()
        ldim_here2 = frame%get_ldim()
        if( any(ldim_here1(1:2) .ne. ldim_here2(1:2)) )then
            THROW_HARD('eval_fdf_foo: dimensions of particle and reference do not match; simple_motion_anisocor_dbl')
        end if
        self%ldim(1) = ldim_here1(1)
        self%ldim(2) = ldim_here1(2)
        self%reference => ref
        call ref%get_rmat_ptr( rmat_ref_ptr )
        do_alloc = .false.
        if (.not. allocated(self%rmat_ref)) then
            do_alloc = .true.
        else
            if ((size(self%rmat_ref, 1) .ne. self%ldim(1)).or.(size(self%rmat_ref, 2) .ne. self%ldim(2))) then
                do_alloc = .true.
            end if
        end if
        if (do_alloc) then
            if (allocated(self%rmat_ref      )) deallocate(self%rmat_ref      )
            if (allocated(self%rmat          )) deallocate(self%rmat          )
            if (allocated(self%rmat_T        )) deallocate(self%rmat_T        )
            if (allocated(self%rmat_T_grad   )) deallocate(self%rmat_T_grad   )
            if (allocated(self%T_coords      )) deallocate(self%T_coords      )
            if (allocated(self%T_coords_xm   )) deallocate(self%T_coords_xm   )
            if (allocated(self%T_coords_xm_sq)) deallocate(self%T_coords_xm_sq)
            if (allocated(self%T_coords_xm_cr)) deallocate(self%T_coords_xm_cr)
            allocate( self%rmat_ref(self%ldim(1), self%ldim(2)) )
            allocate( self%rmat(self%ldim(1), self%ldim(2)) )
            allocate( self%rmat_T         (    self%ldim(1),     self%ldim(2)       ) )
            allocate( self%rmat_T_grad    (    self%ldim(1),     self%ldim(2),    2 ) )
            allocate( self%T_coords       ( 2, self%ldim(1),     self%ldim(2)       ) )
            allocate( self%T_coords_xm    (    self%ldim(1),     self%ldim(2),    2 ) )
            allocate( self%T_coords_xm_sq (    self%ldim(1),     self%ldim(2),    2 ) )
            allocate( self%T_coords_xm_cr (    self%ldim(1),     self%ldim(2)       ) )
            self%sc_alpha_x = 2._dp / (real(self%ldim(1), dp) - 1._dp)
            self%sc_alpha_y = 2._dp / (real(self%ldim(2), dp) - 1._dp)
        end if
        call ref%get_rmat_ptr( rmat_ref_ptr )
        self%rmat_ref(1:self%ldim(1), 1:self%ldim(2)) = rmat_ref_ptr(1:self%ldim(1), 1:self%ldim(2), 1)
        self%frame  => frame
        call frame%get_rmat_ptr( rmat_ptr )
        self%rmat(1:self%ldim(1), 1:self%ldim(2))     = rmat_ptr(1:self%ldim(1), 1:self%ldim(2), 1)
        call self%eval_fdf(a, f, grad)
        
    end subroutine eval_fdf_foo    ! for debug purpopse

    
    
    subroutine eval_fdf( self, a, f, grad )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(motion_anisocor_dbl), intent(inout) :: self
        real(dp),                   intent(in)    :: a(POLY_DIM)
        real(dp),                   intent(out)   :: f, grad(POLY_DIM)
        real(dp) :: N, D_I, D_R, grad_tmp1, grad_tmp2, regu_term
        integer  :: i, j, r
        real(dp) :: poly_dfs(POLY_DIM)
        real(dp), allocatable :: atmp(:,:)
        real(dp) :: regu_grad(12)
!!$        integer :: ithr
!!$        ithr = omp_get_thread_num() + 1
!!$        if( ithr == 1 )then
!!$            write(logfhandle,*) 'eval_fdf, a=', a
!!$        end if
        call self%calc_T_coords(a)
        call self%calc_mat_tld()
        N   = sum( self%rmat_ref(:,:) * self%rmat_T  (:,:) )
        D_R = sum( self%rmat_ref(:,:) * self%rmat_ref(:,:) )
        D_I = sum( self%rmat_T  (:,:) * self%rmat_T  (:,:) )
        f   = N / sqrt(D_I * D_R)
        regu_term = 0.
        do j = 1, self%ldim(2)
            do i = 1, self%ldim(1)
                regu_term = regu_term + (real(i, dp) - self%T_coords(1,i,j))**2 + (real(j, dp) - self%T_coords(2,i,j))**2
            end do
        end do
        f = f - self%gamma * regu_term / real(self%ldim(1)*self%ldim(2),dp)
                
        ! Computing gradient
        ! Polynomial:
        ! Tx = a1 + a2 * x + a3 * y + a4  * x^2 + a5  * y^2 + a6  * x*y
        ! Ty = a7 + a8 * x + a9 * y + a10 * x^2 + a11 * y^2 + a12 * x*y
        ! poly_grad( 1)    = 1.
        ! poly_grad( 2)    = xm
        ! poly_grad( 3)    = ym
        ! poly_grad( 4)    = xms
        ! poly_grad( 5)    = yms
        ! poly_grad( 6)    = xm*ym
        ! poly_grad( 7)    = 1.
        ! poly_grad( 8)    = xm
        ! poly_grad( 9)    = ym
        ! poly_grad(10)    = xms
        ! poly_grad(11)    = yms
        ! poly_grad(12)    = xm*ym
        ! poly_grad( 1: 6) = poly_grad( 1: 6) * (Nx - 1.) / 2.
        ! poly_grad( 7:12) = poly_grad( 7:12) * (Ny - 1.) / 2.
        ! d/da1: poly_dfs=1
        grad_tmp1 = sum( self%rmat_ref(:,:) *                              self%rmat_T_grad(:,:,1) )
        grad_tmp2 = sum( self%rmat_T  (:,:) *                              self%rmat_T_grad(:,:,1) )
        grad(1)   = grad_tmp1 - N *  grad_tmp2 / D_I
        ! d/da2: poly_dfs=xm
        grad_tmp1 = sum( self%rmat_ref(:,:) * self%T_coords_xm(:,:,1)    * self%rmat_T_grad(:,:,1) )
        grad_tmp2 = sum( self%rmat_T  (:,:) * self%T_coords_xm(:,:,1)    * self%rmat_T_grad(:,:,1) )
        grad(2)   = grad_tmp1 - N * grad_tmp2 / D_I
        ! d/da3: poly_dfs=ym
        grad_tmp1 = sum( self%rmat_ref(:,:) * self%T_coords_xm(:,:,2)    * self%rmat_T_grad(:,:,1) )
        grad_tmp2 = sum( self%rmat_T  (:,:) * self%T_coords_xm(:,:,2)    * self%rmat_T_grad(:,:,1) )
        grad(3)   = grad_tmp1 - N * grad_tmp2 / D_I
        ! d/da4: poly_dfs=xm^2
        grad_tmp1 = sum( self%rmat_ref(:,:) * self%T_coords_xm_sq(:,:,1) * self%rmat_T_grad(:,:,1) )
        grad_tmp2 = sum( self%rmat_T  (:,:) * self%T_coords_xm_sq(:,:,1) * self%rmat_T_grad(:,:,1) )
        grad(4)   = grad_tmp1 - N * grad_tmp2 / D_I
        ! d/da5: poly_dfs=ym^2
        grad_tmp1 = sum( self%rmat_ref(:,:) * self%T_coords_xm_sq(:,:,2) * self%rmat_T_grad(:,:,1) )
        grad_tmp2 = sum( self%rmat_T  (:,:) * self%T_coords_xm_sq(:,:,2) * self%rmat_T_grad(:,:,1) )
        grad(5)   = grad_tmp1 - N * grad_tmp2 / D_I
        ! d/da6: poly_dfs=xm*ym
        grad_tmp1 = sum( self%rmat_ref(:,:) * self%T_coords_xm_cr(:,:)   * self%rmat_T_grad(:,:,1) )
        grad_tmp2 = sum( self%rmat_T  (:,:) * self%T_coords_xm_cr(:,:)   * self%rmat_T_grad(:,:,1) )
        grad(6)   = grad_tmp1 - N * grad_tmp2 / D_I
        ! d/da7: poly_dfs=1
        grad_tmp1 = sum( self%rmat_ref(:,:) *                              self%rmat_T_grad(:,:,2) )
        grad_tmp2 = sum( self%rmat_T  (:,:) *                              self%rmat_T_grad(:,:,2) )
        grad(7)   = grad_tmp1 - N * grad_tmp2 / D_I
        ! d/da8: poly_dfs=xm
        grad_tmp1 = sum( self%rmat_ref(:,:) * self%T_coords_xm(:,:,1)    * self%rmat_T_grad(:,:,2) )
        grad_tmp2 = sum( self%rmat_T  (:,:) * self%T_coords_xm(:,:,1)    * self%rmat_T_grad(:,:,2) )
        grad(8)   = grad_tmp1 - N * grad_tmp2 / D_I
        ! d/da9: poly_dfs=ym
        grad_tmp1 = sum( self%rmat_ref(:,:) * self%T_coords_xm(:,:,2)    * self%rmat_T_grad(:,:,2) )
        grad_tmp2 = sum( self%rmat_T  (:,:) * self%T_coords_xm(:,:,2)    * self%rmat_T_grad(:,:,2) )
        grad(9)   = grad_tmp1 - N * grad_tmp2 / D_I
        ! d/da10: poly_dfs=xm^2
        grad_tmp1 = sum( self%rmat_ref(:,:) * self%T_coords_xm_sq(:,:,1) * self%rmat_T_grad(:,:,2) )
        grad_tmp2 = sum( self%rmat_T  (:,:) * self%T_coords_xm_sq(:,:,1) * self%rmat_T_grad(:,:,2) )
        grad(10)  = grad_tmp1 - N * grad_tmp2 / D_I
        ! d/da11: poly_dfs=ym^2
        grad_tmp1 = sum( self%rmat_ref(:,:) * self%T_coords_xm_sq(:,:,2) * self%rmat_T_grad(:,:,2) )
        grad_tmp2 = sum( self%rmat_T  (:,:) * self%T_coords_xm_sq(:,:,2) * self%rmat_T_grad(:,:,2) )
        grad(11)  = grad_tmp1 - N * grad_tmp2 / D_I
        ! d/da12: poly_dfs=xm*ym
        grad_tmp1 = sum( self%rmat_ref(:,:) * self%T_coords_xm_cr(:,:)   * self%rmat_T_grad(:,:,2) )
        grad_tmp2 = sum( self%rmat_T  (:,:) * self%T_coords_xm_cr(:,:)   * self%rmat_T_grad(:,:,2) )
        grad(12)  = grad_tmp1 - N * grad_tmp2 / D_I
        grad(1: 6) = grad(1: 6) / self%sc_alpha_x
        grad(7:12) = grad(7:12) / self%sc_alpha_y
        grad = grad / sqrt(D_I*D_R)
        call self%calc_grad_coords(a, regu_grad)
        grad = grad - self%gamma * regu_grad / real(self%ldim(1)*self%ldim(2),dp)
    end subroutine eval_fdf
    
    subroutine eval_f_noregu( self, a, f, regu )      ! for debug purposes
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(motion_anisocor_dbl), intent(inout) :: self
        real(dp),                   intent(in)    :: a(POLY_DIM)
        real(dp),                   intent(out)   :: f, regu
        real(dp) :: N, D_I, D_R, grad_tmp1, grad_tmp2
        integer  :: i, j, r
        real(dp) :: poly_dfs(POLY_DIM)
        real(dp), allocatable :: atmp(:,:)
        real(dp) :: regu_term
        integer :: ithr
!!$        integer :: ithr
!!$        ithr = omp_get_thread_num() + 1
!!$        if( ithr == 1 )then
!!$            write(logfhandle,*) 'eval_fdf, a=', a
!!$        end if
        call self%calc_T_coords(a)
        call self%calc_mat_tld()
        N   = sum( self%rmat_ref(:,:) * self%rmat_T  (:,:) )
        D_R = sum( self%rmat_ref(:,:) * self%rmat_ref(:,:) )
        D_I = sum( self%rmat_T  (:,:) * self%rmat_T  (:,:) )
        f   = N / sqrt(D_I * D_R)
        regu_term = 0.
        
        ithr = omp_get_thread_num() + 1
        if ((.false.).and.(ithr == 1)) then
            open(123, file='acoords')
            write (123,*) 'a=', a
            write (123,*) 'Nx=', self%ldim(1)
            write (123,*) 'Ny=', self%ldim(2)
            write (123,*) 'T_coords_1=[...'
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    write (123,'(A)', advance='no') trim(dbl2str(self%T_coords(1,i,j)))
                    if (j < self%ldim(2)) write (123,'(A)', advance='no') ', '                    
                end do
                if (i < self%ldim(1)) write (123,'(A)') '; ...'
            end do
            write (213,*) '];'
            write (123,*) 'T_coords_2=[...'
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    write (123,'(A)', advance='no') trim(dbl2str(self%T_coords(2,i,j)))
                    if (j < self%ldim(2)) write (123,'(A)', advance='no') ', '                    
                end do
                if (i < self%ldim(1)) write (123,'(A)') '; ...'
            end do
            write (213,*) '];'
            close(123)
            stop 1
        end if
        do j = 1, self%ldim(2)
            do i = 1, self%ldim(1)
                regu_term = regu_term + (real(i, dp) - self%T_coords(1,i,j))**2 + (real(j, dp) - self%T_coords(2,i,j))**2
            end do
        end do
        write (*,*) 'regu_term_pre:', regu_term
        regu_term = regu_term / real(self%ldim(1)*self%ldim(2),dp)
        write (*,*) 'regu_term_post:', regu_term
        regu = regu_term * self%gamma
        write (*,*) 'regu:', regu
    end subroutine eval_f_noregu     ! for debug purposes


    !> Initialise  ftexp_shsrch
    subroutine motion_anisocor_dbl_new( self, motion_correct_ftol, motion_correct_gtol )
        class(motion_anisocor_dbl), intent(inout) :: self
        real,         optional,     intent(in)    :: motion_correct_ftol, motion_correct_gtol
        type(opt_factory) :: ofac !< optimizer factory
        integer           :: i
        real              :: lims(POLY_DIM,2)
        call self%kill()
        self%gamma      = 0.01
        self%maxHWshift = 0.2
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
        lims(:,1) = - self%maxHWshift
        lims(:,2) =   self%maxHWshift
        call self%ospec%specify('lbfgsb', POLY_DIM, ftol=self%motion_correctftol, gtol=self%motion_correctgtol, limits=lims)
        call self%ospec%set_costfun_8   (motion_anisocor_dbl_cost_8)
        call self%ospec%set_gcostfun_8  (motion_anisocor_dbl_gcost_8)
        call self%ospec%set_fdfcostfun_8(motion_anisocor_dbl_fdfcost_8)
        ! generate optimizer object with the factory
        if( associated(self%nlopt) )then
            call self%nlopt%kill()
            deallocate(self%nlopt)
        end if
        self%ospec%x   = 0.
        self%ospec%x_8 = 0._dp
        call ofac%new(self%ospec, self%nlopt)
    end subroutine motion_anisocor_dbl_new

    !> Cost function, double precision
    function motion_anisocor_dbl_cost_8( fun_self, vec, D ) result(cost)
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(in)    :: vec(D)
        real(kind=8)                :: cost
        real(kind=8)                :: f, grad(POLY_DIM)
        select type(fun_self)
        class is (motion_anisocor_dbl)
            call fun_self%eval_fdf(vec, f, grad)
            cost = -f
        class default
            THROW_HARD('error in motion_anisocor_dbl_cost_8: unknown type; simple_motion_anisocor_dbl')
        end select
    end function motion_anisocor_dbl_cost_8

    !> Gradient function, double precision
    subroutine motion_anisocor_dbl_gcost_8( fun_self, vec, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(dp),     intent(inout) :: vec( D )
        real(dp),     intent(out)   :: grad( D )
        real(dp)                    :: f, agrad( D )
        select type(fun_self)
        class is (motion_anisocor_dbl)
            call fun_self%eval_fdf(vec, f, agrad)
            grad = -agrad
        class default
            THROW_HARD('error in motion_anisocor_dbl_gcost_8: unknown type; simple_motion_anisocor_dbl')
        end select
    end subroutine motion_anisocor_dbl_gcost_8

    !> Gradient & cost function, double precision
    subroutine motion_anisocor_dbl_fdfcost_8( fun_self, vec, f, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
        real(kind=8)                :: agrad( D )
        select type(fun_self)
        class is (motion_anisocor_dbl)
            call fun_self%eval_fdf(vec, f, agrad)
            f    = -f
            grad = -agrad
        class default
            THROW_HARD('error in motion_anisocor_dbl_fdfcost_8: unknown type; simple_motion_anisocor_dbl')
        end select
    end subroutine motion_anisocor_dbl_fdfcost_8

    !> Main search routine
    function motion_anisocor_dbl_minimize( self, ref, frame, corr, regu ) result(cxy)
        class(motion_anisocor_dbl), intent(inout) :: self
        class(image), target,       intent(in)    :: ref, frame
        real(dp),                   intent(out)   :: corr, regu
        real(dp)          :: cxy(POLY_DIM+1) ! corr + model
        real              :: cxy1
        class(*), pointer :: fun_self => null()
        integer           :: i, j
        real, pointer     :: rmat_ref_ptr(:,:,:), rmat_ptr(:,:,:)
        integer           :: ldim_here1(3), ldim_here2(3)
        logical           :: do_alloc
        integer           :: irestart
        write (*,*) 'gamma=', self%gamma ; call flush(6)
        ldim_here1 = ref  %get_ldim()
        ldim_here2 = frame%get_ldim()
        if( any(ldim_here1(1:2) .ne. ldim_here2(1:2)) )then
            THROW_HARD('motion_anisocor_dbl_minimize: dimensions of particle and reference do not match; simple_motion_anisocor_dbl')
        end if
        self%ldim(1) = ldim_here1(1)
        self%ldim(2) = ldim_here1(2)
        self%reference => ref
        call ref%get_rmat_ptr( rmat_ref_ptr )
        do_alloc = .false.
        if (.not. allocated(self%rmat_ref)) then
            do_alloc = .true.
        else
            if ((size(self%rmat_ref, 1) .ne. self%ldim(1)).or.(size(self%rmat_ref, 2) .ne. self%ldim(2))) then
                do_alloc = .true.
            end if
        end if
        if (do_alloc) then
            if (allocated(self%rmat_ref      )) deallocate(self%rmat_ref      )
            if (allocated(self%rmat          )) deallocate(self%rmat          )
            if (allocated(self%rmat_T        )) deallocate(self%rmat_T        )
            if (allocated(self%rmat_T_grad   )) deallocate(self%rmat_T_grad   )
            if (allocated(self%T_coords      )) deallocate(self%T_coords      )
            if (allocated(self%T_coords_xm   )) deallocate(self%T_coords_xm   )
            if (allocated(self%T_coords_xm_sq)) deallocate(self%T_coords_xm_sq)
            if (allocated(self%T_coords_xm_cr)) deallocate(self%T_coords_xm_cr)
            allocate( self%rmat_ref(self%ldim(1), self%ldim(2)) )
            allocate( self%rmat(self%ldim(1), self%ldim(2)) )
            allocate( self%rmat_T         (    self%ldim(1),     self%ldim(2)       ) )
            allocate( self%rmat_T_grad    (    self%ldim(1),     self%ldim(2),    2 ) )
            allocate( self%T_coords       ( 2, self%ldim(1),     self%ldim(2)       ) )
            allocate( self%T_coords_xm    (    self%ldim(1),     self%ldim(2),    2 ) )
            allocate( self%T_coords_xm_sq (    self%ldim(1),     self%ldim(2),    2 ) )
            allocate( self%T_coords_xm_cr (    self%ldim(1),     self%ldim(2)       ) )
            self%sc_alpha_x = 2._dp / (real(self%ldim(1), dp) - 1._dp)
            self%sc_alpha_y = 2._dp / (real(self%ldim(2), dp) - 1._dp)
        end if
        call ref%get_rmat_ptr( rmat_ref_ptr )
        self%rmat_ref(1:self%ldim(1), 1:self%ldim(2)) = rmat_ref_ptr(1:self%ldim(1), 1:self%ldim(2), 1)
        self%frame  => frame
        call frame%get_rmat_ptr( rmat_ptr )
        self%rmat(1:self%ldim(1), 1:self%ldim(2))     = rmat_ptr(1:self%ldim(1), 1:self%ldim(2), 1)
        self%ospec%x   = randn(POLY_DIM) * self%maxHWshift * 0.001 ! randomized (yet small) starting point
        self%ospec%x_8 = real(self%ospec%x, dp)
        self%ospec%maxits = 400
        irestart = 0
        self%ospec%converged = .false.
        do while ((irestart < 5).and.(.not. self%ospec%converged))
            irestart = irestart + 1
            call self%nlopt%minimize(self%ospec, self, cxy1)
        end do
        if (irestart > 1) then            
            if (.not. self%ospec%converged) then
                write (*,*) 'QQQQQQQQQQQQQQQQQQQQQQQQQQ not converged!'
            else
                write (*,*) 'converged with irestart=', irestart
            end if
        end if
        if (self%ospec%converged) then
            cxy(1)  = real(cxy1)
            cxy(2:) = self%ospec%x_8
            call self%eval_f_noregu(self%ospec%x_8, corr, regu)
        else
            cxy(2:) = 0.
            call self%eval_f_noregu(cxy(2:), cxy(1), regu)
            corr = cxy(1)
        end if
    end function motion_anisocor_dbl_minimize

    subroutine motion_anisocor_dbl_kill( self )
        class(motion_anisocor_dbl), intent(inout) :: self
        if (allocated(self%rmat       )) deallocate(self%rmat       )
        if (allocated(self%rmat_ref   )) deallocate(self%rmat_ref   )
        if (allocated(self%rmat_T     )) deallocate(self%rmat_T     )
        if (allocated(self%rmat_T_grad)) deallocate(self%rmat_T_grad)
        if (associated(self%nlopt)) then
            call self%nlopt%kill()
            self%nlopt => null()
        end if
        if (allocated(self%T_coords      )) deallocate(self%T_coords      )
        if (allocated(self%T_coords_xm   )) deallocate(self%T_coords_xm   )
        if (allocated(self%T_coords_xm_sq)) deallocate(self%T_coords_xm_sq)
        if (allocated(self%T_coords_xm_cr)) deallocate(self%T_coords_xm_cr)
        if (allocated(self%T_coords_out  )) deallocate(self%T_coords_out  )
        call self%ospec%kill()
        self%maxHWshift  = 0.
        self%reference   => null()
        self%frame       => null()
        self%ldim        = 0
    end subroutine motion_anisocor_dbl_kill

end module simple_motion_anisocor_dbl
