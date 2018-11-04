! parametric image registration using a biquadratic polynomial with L-BFGS-B in real-space (used in motion_correct)
module simple_motion_anisocor
include 'simple_lib.f08'
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_image,       only: image
implicit none
private
public :: motion_anisocor_init, motion_anisocor_set_ptrs, motion_anisocor_minimize   !, test_motion_anisocor
#include "simple_local_flags.inc"

real,    parameter                  :: TOL=1e-4           !< tolerance parameter
integer, parameter                  :: MAXITS=30          !< maximum number of iterations

type :: motion_anisocor
    type(opt_factory)               :: ofac               !< optimizer factory
    type(opt_spec)                  :: ospec              !< optimizer specification object
    class(optimizer),   pointer     :: nlopt=>null()      !< pointer to nonlinear optimizer
    real                            :: maxHWshift = 0.    !< maximum half-width of shift
    integer,            allocatable :: ldim(:)   
    class(image),       pointer     :: reference=>null()  !< reference pft
    class(image),       pointer     :: particle =>null()  !< particle pft
    real(kind=c_float), pointer     :: rmat_ref   (:,:,:) !< reference matrix
    real(kind=c_float), pointer     :: rmat       (:,:,:) !< particle  matrix        
    real(kind=c_float), allocatable :: rmat_T     (:,:)   !< transformed (interpolated) particle matrix
    real(kind=c_float), allocatable :: rmat_T_grad(:,:,:) !< gradient of transformed particle matrix    
contains
    procedure                       :: motion_anisocor_init
    procedure                       :: motion_anisocor_set_ptrs
    procedure                       :: motion_anisocor_minimize
    procedure                       :: eval_fdf
    procedure                       :: biquad_poly_df
    procedure                       :: biquad_poly
    procedure                       :: calc_mat_tld
    procedure                       :: alloc_hlpmats
    procedure                       :: interp_bilin
    procedure                       :: interp_bilin_fdf
end type motion_anisocor

contains

    subroutine alloc_hlpmats( self )
        class(motion_anisocor), intent(inout) :: self
        logical                               :: do_allocate
        do_allocate = .false.
        if (.not.(allocated(self%rmat_T))) then
            do_allocate = .true.
        else
            if ((size(self%rmat_T, 1) .ne. self%ldim(1)) .or. &
                (size(self%rmat_T, 2) .ne. self%ldim(2))) do_allocate = .true.
        end if
        if (do_allocate) then
            if (allocated(self%rmat_T))      deallocate(self%rmat_T)
            if (allocated(self%rmat_T_grad)) deallocate(self%rmat_T_grad)
            allocate(self%rmat_T     (self%ldim(1),self%ldim(2)),   &
                     self%rmat_T_grad(self%ldim(1),self%ldim(2),2)     )
        end if
    end subroutine alloc_hlpmats

    ! bilinear interpolation
    function interp_bilin( self, x ) result( val )
        class(motion_anisocor), intent(inout) :: self
        real,                   intent(in)    :: x(2)
        real                                  :: val
        integer                               :: x1,x2,y1,y2
        real                                  :: NP1, NP2, NP3, NP4
        real                                  :: PW1, PW2, PW3, PW4
        ! Any values out of acceptable range
        if ((x(1) < 1.).or.(x(1) > real(self%ldim(1))).or.(x(2) < 1.).or.(x(2) > real(self%ldim(2)))) then
            val = 0.
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
        val = PW1 * NP1 + PW2 * NP2 + PW3 * NP3 + PW4 * NP4
    end function interp_bilin

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

    subroutine calc_mat_tld( self, a )
        class(motion_anisocor), intent(inout) :: self
        real(dp),     intent(in)    :: a(12)
        integer                     :: i, j
        real                        :: x_tld(2)
        real                        :: val, grad(2)
        do i = 1, self%ldim(1)
            do j = 1, self%ldim(2)
                call self%biquad_poly(i, j, a, x_tld)
                call self%interp_bilin_fdf(x_tld, val, grad)
                self%rmat_T     (i,j)   = val
                self%rmat_T_grad(i,j,1) = grad(1)
                self%rmat_T_grad(i,j,2) = grad(2)
            end do
        end do        
    end subroutine calc_mat_tld

    subroutine biquad_poly( self, i, j, a, x_tld)
        class(motion_anisocor), intent(inout) :: self
        integer,      intent(in)    :: i, j
        real(dp),     intent(in)    :: a(12)
        real,         intent(out)   :: x_tld(2)
        real                        :: x, y, xm, ym, xms, yms, z1, z2
        real                        :: Nx, Ny
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

    subroutine biquad_poly_df( self, i, j, a, grad )
        class(motion_anisocor), intent(inout) :: self
        integer,      intent(in)    :: i, j
        real(dp),     intent(in)    :: a(12)
        real                        :: grad(12)
        real                        :: x, y, xm, ym, xms, yms
        real                        :: Nx, Ny
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
        grad( 1) = 1.
        grad( 2) = xm
        grad( 3) = ym
        grad( 4) = xms
        grad( 5) = yms
        grad( 6) = xm*ym
        grad( 7) = 1.
        grad( 8) = xm
        grad( 9) = ym
        grad(10) = xms
        grad(11) = yms
        grad(12) = xm*ym
        grad( 1: 6) = grad( 1: 6) * (Nx - 1.) / 2.
        grad( 7:12) = grad( 7:12) * (Ny - 1.) / 2.
    end subroutine biquad_poly_df
    
    subroutine eval_fdf( self, a, f_dp, grad )
        class(motion_anisocor), intent(inout) :: self
        real(dp),               intent(in)    :: a(12)
        real(dp),               intent(out)   :: f_dp
        real(dp),               intent(out)   :: grad(12)
        real                                  :: f
        real                                  :: N, D_I, D_R        
        integer                               :: i, j
        integer                               :: r
        real                                  :: grad_tmp1
        real                                  :: grad_tmp2
        real                                  :: poly_dfs(12)
        call self%calc_mat_tld(a)        
        N   = 0.
        D_I = 0.
        D_R = 0.
        do i = 1, self%ldim(1)
            do j = 1, self%ldim(2)
                f   = f   + self%rmat  (i,j,1) * self%rmat_T(i,j  )
                D_I = D_I + self%rmat  (i,j,1) * self%rmat  (i,j,1)
                D_R = D_R + self%rmat_T(i,j  ) * self%rmat_T(i,j  )                
            end do
        end do
        f = N / (D_I * D_R)
        ! Computing gradient
        grad_tmp1 = 0.
        grad_tmp2 = 0.
        do r = 1, 12
            do i = 1, self%ldim(1)
                do j = 1, self%ldim(2)
                    call self%biquad_poly(i, j, a, poly_dfs)   ! TODO: take this out of the r-loop
                    if ( r <= 6) then                        
                        grad_tmp1 = grad_tmp1 + self%rmat  (i,j,1) * poly_dfs(r) * self%rmat_T_grad(i,j,1)
                        grad_tmp2 = grad_tmp2 + self%rmat_T(i,j  ) * poly_dfs(r) * self%rmat_T_grad(i,j,1)
                    else
                        grad_tmp1 = grad_tmp1 + self%rmat  (i,j,1) * poly_dfs(r) * self%rmat_T_grad(i,j,2)
                        grad_tmp2 = grad_tmp2 + self%rmat_T(i,j  ) * poly_dfs(r) * self%rmat_T_grad(i,j,2)                        
                    end if
                end do
            end do
            grad(r) = 1./(D_I*D_R)*(grad_tmp1 - f * D_I * 2. * grad_tmp2)
        end do        
    end subroutine eval_fdf

    
    !> Initialise  ftexp_shsrch
    subroutine motion_anisocor_init( self, ref, ptcl, trs, motion_correct_ftol, motion_correct_gtol )
        class(motion_anisocor), intent(inout) :: self
        class(image), target, intent(in) :: ref, ptcl
        real,                 intent(in) :: trs
        real,       optional, intent(in) :: motion_correct_ftol, motion_correct_gtol
!!$        if( present(motion_correct_ftol) )then
!!$            motion_correctftol = motion_correct_ftol
!!$        else
!!$            motion_correctftol = TOL
!!$        end if
!!$        if( present(motion_correct_gtol) )then
!!$            motion_correctgtol = motion_correct_gtol
!!$        else
!!$            motion_correctgtol = TOL 
!!$        end if
!!$        reference  => ref
!!$        particle   => ptcl
!!$        maxHWshift =  trs
    end subroutine motion_anisocor_init

    subroutine motion_anisocor_set_ptrs( self, ref, ptcl )
        class(motion_anisocor), intent(inout)      :: self
        class(image),           intent(in), target :: ref, ptcl
        self%reference  => ref
        self%particle   => ptcl
        self%ldim       =  ref%get_ldim()
    end subroutine motion_anisocor_set_ptrs

    !> Cost function, double precision
    function motion_anisocor_cost_8( fun_self, vec, D ) result(cost)
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(in)    :: vec(D)
        real(kind=8)                :: cost
!!$        cost = -reference%corr_shifted_8(particle, -vec)
    end function motion_anisocor_cost_8

    !> Gradient function, double precision
    subroutine motion_anisocor_gcost_8( fun_self, vec, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(dp), intent(inout) :: vec( D )
        real(dp), intent(out)   :: grad( D )
!!$        call reference%corr_gshifted_8(particle, -vec, grad)
    end subroutine motion_anisocor_gcost_8

    !> Gradient & cost function, double precision
    subroutine motion_anisocor_fdfcost_8( fun_self, vec, f, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
!!$        call reference%corr_fdfshifted_8(particle, -vec, f, grad)
!!$        f    = f    * (-1.0_dp)
    end subroutine motion_anisocor_fdfcost_8

    !> Main search routine
    function motion_anisocor_minimize( self, prev_corr, prev_shift ) result( cxy )
        class(motion_anisocor), intent(inout) :: self
        real, optional, intent(in) :: prev_corr, prev_shift(2)
        real :: cxy(3), lims(2,2) ! max_shift, maxlim
        class(*), pointer :: fun_self => null()
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
