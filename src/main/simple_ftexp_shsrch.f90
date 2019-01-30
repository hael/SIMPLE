! shift search with L-BFGS-B using expanded Fourier transforms (used in motion_correct)
module simple_ftexp_shsrch
include 'simple_lib.f08'
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_ft_expanded, only: ft_expanded

implicit none

public :: ftexp_shsrch, test_ftexp_shsrch
private
#include "simple_local_flags.inc"

real,    parameter :: TOL    = 1e-4 !< tolerance parameter
integer, parameter :: MAXITS = 30   !< maximum number of iterations

complex(dp), parameter   :: J = CMPLX(0.0_dp, 1.0_dp, kind=dp)
real(dp),    parameter   :: denom = 0.00075_dp ! denominator for rescaling of cost function

type :: ftexp_shsrch
    type(opt_spec)              :: ospec                     !< optimizer specification object
    class(optimizer),   pointer :: nlopt     => null()       !< pointer to nonlinear optimizer
    class(ft_expanded), pointer :: reference => null()       !< reference ft_exp
    class(ft_expanded), pointer :: particle  => null()       !< particle ft_exp
    real                        :: maxHWshift         = 0.   !< maximum half-width of shift
    real                        :: motion_correctftol = 1e-4 !< function error tolerance
    real                        :: motion_correctgtol = 1e-4 !< gradient error tolerance
    integer                     :: lims(3,2)                 !< physical limits for the Fourier transform
    integer                     :: flims(3,2)                !< shifted limits
    integer                     :: ldim(3)                   !< logical dimension
    integer                     :: flims_nyq(3,2)            !< nyquist limits
    real                        :: ptcl_smpd                 !< sampling distance
    real(dp),    allocatable    :: ftexp_tmpmat_re_2d(:,:)   !< temporary matrix for shift search
    real(dp),    allocatable    :: ftexp_tmpmat_im_2d(:,:)   !< temporary matrix for shift search
    complex(dp), allocatable    :: ftexp_tmp_cmat12(:,:)     !< temporary matrix for shift search
    logical                     :: existence
contains
    procedure          :: new            => ftexp_shsrch_new
    procedure          :: minimize       => ftexp_shsrch_minimize
    procedure          :: corr_shifted_8 => ftexp_shsrch_corr_shifted_8
    procedure          :: kill           => ftexp_shsrch_kill
    procedure, private :: set_dims_and_alloc                 !< set dimensions from images and allocate tmp matrices
    procedure, private :: corr_shifted_cost_8                !< cost function for minimizer, f only
    procedure, private :: corr_gshifted_cost_8               !< cost function for minimizer, gradient only
    procedure, private :: corr_fdfshifted_cost_8             !< cost function for minimizer, f and gradient
    procedure, private :: calc_tmpmat_re                     !< calculate tmp matrix for cost function
    procedure, private :: calc_tmpmat_im                     !< calculate tmp matrix for cost function
    procedure, private :: calc_tmpmat_re_im                  !< calculate tmp matrix for cost function
    procedure, private :: calc_tmp_cmat12                    !< calculate tmp matrix for cost function
end type ftexp_shsrch

contains

    !> Initialise  ftexp_shsrch
    subroutine ftexp_shsrch_new( self, ref, ptcl, trs, motion_correct_ftol, motion_correct_gtol )
        use simple_opt_factory, only: opt_factory
        class(ftexp_shsrch),        intent(inout) :: self
        class(ft_expanded), target, intent(in)    :: ref, ptcl
        real,                       intent(in)    :: trs
        real,             optional, intent(in)    :: motion_correct_ftol, motion_correct_gtol
        type(opt_factory) :: ofac
        real              :: opt_lims(2,2)
        call self%kill()
        self%reference  => ref
        self%particle   => ptcl
        self%maxHWshift =  trs
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
        opt_lims(1,1) = - self%maxHWshift
        opt_lims(1,2) =   self%maxHWshift
        opt_lims(2,1) = - self%maxHWshift
        opt_lims(2,2) =   self%maxHWshift
        call self%ospec%specify('lbfgsb', 2, ftol=self%motion_correctftol, gtol=self%motion_correctgtol, limits=opt_lims)
        call self%ospec%set_costfun_8(ftexp_shsrch_cost_8)
        call self%ospec%set_gcostfun_8(ftexp_shsrch_gcost_8)
        call self%ospec%set_fdfcostfun_8(ftexp_shsrch_fdfcost_8)
        ! generate optimizer object with the factory
        if( associated(self%nlopt) )then
            call self%nlopt%kill
            deallocate(self%nlopt)
        end if
        call ofac%new(self%ospec, self%nlopt)
        self%existence = .true.
    end subroutine ftexp_shsrch_new

    !> Main search routine
    function ftexp_shsrch_minimize( self, prev_corr, prev_shift ) result( cxy )
        class(ftexp_shsrch), intent(inout) :: self
        real, optional,      intent(in) :: prev_corr, prev_shift(2)
        real :: cxy(3)
        self%ospec%limits(1,1) = - self%maxHWshift
        self%ospec%limits(1,2) =   self%maxHWshift
        self%ospec%limits(2,1) = - self%maxHWshift
        self%ospec%limits(2,2) =   self%maxHWshift
        if( present(prev_shift) )then
            self%ospec%limits(1,:) = self%ospec%limits(1,:) + prev_shift(1)
            self%ospec%limits(2,:) = self%ospec%limits(2,:) + prev_shift(2)
        endif
        if( present(prev_shift) ) self%ospec%x = prev_shift
        call self%set_dims_and_alloc()
        call self%calc_tmp_cmat12()
        ! set initial solution to previous shift
        call self%nlopt%minimize(self%ospec, self, cxy(1))
        call self%reference%corr_normalize(self%particle, cxy(1))
        cxy(1)  = -cxy(1) ! correlation
        cxy(2:) = self%ospec%x ! shift
        if( present(prev_corr) )then
            if( abs(cxy(1)-prev_corr) <= TOL )then
                cxy(1)  = prev_corr
                if( present(prev_shift) ) cxy(2:) = prev_shift
            endif
        endif
    end function ftexp_shsrch_minimize

    subroutine ftexp_shsrch_kill( self )
        class(ftexp_shsrch), intent(inout) :: self
        if ( self%existence ) then
            call self%ospec%kill
            if( associated( self%nlopt ) )then
                call self%nlopt%kill
                nullify(self%nlopt)
            end if
            if ( associated( self%reference ) ) self%reference => null()
            if ( associated( self%particle )  ) self%particle  => null()
            if ( allocated( self%ftexp_tmpmat_re_2d ) ) deallocate( self%ftexp_tmpmat_re_2d )
            if ( allocated( self%ftexp_tmpmat_im_2d ) ) deallocate( self%ftexp_tmpmat_im_2d )
            if ( allocated( self%ftexp_tmp_cmat12   ) ) deallocate( self%ftexp_tmp_cmat12   )
            self%existence = .false.
        end if
    end subroutine ftexp_shsrch_kill

    !< set dimensions from images and allocate tmp matrices
    subroutine set_dims_and_alloc( self )
        class(ftexp_shsrch), intent(inout) :: self
        integer :: ref_flims (3,2), ref_ldim (3), ref_lims (3,2), ref_flims_nyq (3,2)
        integer :: ptcl_flims(3,2), ptcl_ldim(3), ptcl_lims(3,2), ptcl_flims_nyq(3,2)
        logical :: do_alloc
        ref_flims  = self%reference%get_flims()
        ptcl_flims = self%particle %get_flims()
        if (any( ref_flims /= ptcl_flims ) ) then
            THROW_HARD('set_dims_and_alloc: ptcl and ref have inconsistens flims; simple_ftexp_shsrch')
        end if
        ref_ldim  = self%reference%get_ldim()
        ptcl_ldim = self%particle %get_ldim()
        if ( any( ref_ldim /= ptcl_ldim ) ) then
            THROW_HARD('set_dims_and_alloc: ptcl and ref have inconsistens ldim; simple_ftexp_shsrch')
        end if
        ref_lims  = self%reference%get_lims()
        ptcl_lims = self%particle %get_lims()
        if ( any( ref_lims /= ptcl_lims ) ) then
            THROW_HARD('set_dims_and_alloc: ptcl and ref have inconsistens lims; simple_ftexp_shsrch')
        end if
        ref_flims_nyq  = self%reference%get_flims_nyq()
        ptcl_flims_nyq = self%particle %get_flims_nyq()
        if ( any( ref_flims_nyq /= ptcl_flims_nyq ) ) then
            THROW_HARD('set_dims_and_alloc: ptcl and ref have inconsistens flims_nyq; simple_ftexp_shsrch')
        end if
        self%flims     = ref_flims
        self%ldim      = ref_ldim
        self%lims      = ref_lims
        self%flims_nyq = ref_flims_nyq
        do_alloc = .true.
        if ( allocated( self%ftexp_tmpmat_re_2d ) ) then
            if ( ( ubound( self%ftexp_tmpmat_re_2d, 1 ) == ref_flims_nyq(1,2) ) .and. &
                 ( ubound( self%ftexp_tmpmat_re_2d, 2 ) == ref_flims_nyq(2,2) ) ) then
                do_alloc = .false.
            else
                deallocate( self%ftexp_tmpmat_re_2d )
                deallocate( self%ftexp_tmpmat_im_2d )
                deallocate( self%ftexp_tmp_cmat12 )
            end if
        end if
        if ( do_alloc ) then
            allocate( self%ftexp_tmpmat_re_2d( 1:ref_flims_nyq(1,2), 1:ref_flims_nyq(2,2) ),&
                      self%ftexp_tmpmat_im_2d( 1:ref_flims_nyq(1,2), 1:ref_flims_nyq(2,2) ),&
                      self%ftexp_tmp_cmat12  ( 1:ref_flims_nyq(1,2), 1:ref_flims_nyq(2,2) ),&
                      stat=alloc_stat ) ! note: we could allocate 1:flims(1,2), 1:flims(2,2) elements to save memory
            if (alloc_stat /= 0) call allocchk('In: set_dims_and_alloc; simple_ftexp_shsrch')
        end if
    end subroutine set_dims_and_alloc

    !> Cost function, double precision
    function ftexp_shsrch_cost_8( self, vec, D ) result( cost )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(kind=8), intent(in)    :: vec(D)
        real(kind=8) :: cost
        select type(self)
            class is (ftexp_shsrch)
                !cost = -self%reference%corr_shifted_8(self%particle, -vec)
                cost = -self%corr_shifted_cost_8( -vec )
            class DEFAULT
                THROW_HARD('unknown type; ftexp_shsrch_cost_8')
        end select
    end function ftexp_shsrch_cost_8

    !< cost function for minimizer, f only
    function corr_shifted_cost_8( self, shvec ) result( r )
        class(ftexp_shsrch), intent(inout) :: self
        real(dp),            intent(in)    :: shvec(2)
        real(dp) :: r
        call self%calc_tmpmat_re( shvec )
        r = sum(self%ftexp_tmpmat_re_2d(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2))) / denom
    end function corr_shifted_cost_8

    !< cost function for minimizer, gradient only
    subroutine corr_gshifted_cost_8( self, shvec, grad )
        class(ftexp_shsrch), intent(inout) :: self
        real(dp),            intent(in)    :: shvec(2)
        real(dp),            intent(out)   :: grad(2)
        real, pointer :: transfmat_here(:,:,:,:)
        call self%particle%get_transfmat_ptr( transfmat_here )
        call self%calc_tmpmat_im( shvec )
        grad(1) = sum(self%ftexp_tmpmat_im_2d(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2)) *&
                   &transfmat_here(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2),1,1))
        grad(2) = sum(self%ftexp_tmpmat_im_2d(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2)) *&
                   &transfmat_here(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2),1,2))
        grad = grad / denom
    end subroutine corr_gshifted_cost_8

    !< cost function for minimizer, f and gradient
    subroutine corr_fdfshifted_cost_8( self, shvec, f, grad )
        class(ftexp_shsrch), intent(inout) :: self
        real(dp),            intent(in)    :: shvec(2)
        real(dp),            intent(out)   :: grad(2), f
        real, pointer :: transfmat_here(:,:,:,:)
        call self%particle%get_transfmat_ptr( transfmat_here )
        call self%calc_tmpmat_re_im( shvec )
        f       = sum(self%ftexp_tmpmat_re_2d(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2))) / denom
        grad(1) = sum(self%ftexp_tmpmat_im_2d(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2)) *&
                   &transfmat_here(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2),1,1))
        grad(2) = sum(self%ftexp_tmpmat_im_2d(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2)) *&
                   &transfmat_here(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2),1,2))
        grad    = grad / denom
    end subroutine corr_fdfshifted_cost_8

    !< calculate tmp matrix for cost function
    subroutine calc_tmpmat_re( self, shvec )
        class(ftexp_shsrch), intent(inout) :: self
        real(dp),            intent(in)    :: shvec(2)
        real(dp) :: arg
        integer  :: hind,kind
        real, pointer :: transfmat_here(:,:,:,:)
        call self%particle%get_transfmat_ptr(transfmat_here)
        do kind=self%flims(2,1),self%flims(2,2)
            do hind=self%flims(1,1),self%flims(1,2)
                arg  = dot_product(shvec(:), transfmat_here(hind,kind,1,1:2))
                self%ftexp_tmpmat_re_2d(hind,kind) = real(self%ftexp_tmp_cmat12(hind,kind) * exp(-J * arg),kind=dp)
            end do
        end do
    end subroutine calc_tmpmat_re

    !< calculate tmp matrix for cost function
    subroutine calc_tmpmat_im( self, shvec )
        class(ftexp_shsrch), intent(inout) :: self
        real(dp),            intent(in)    :: shvec(2)
        real(dp) :: arg
        integer  :: hind,kind
        real, pointer :: transfmat_here(:,:,:,:)
        call self%particle%get_transfmat_ptr(transfmat_here)
        do kind=self%flims(2,1),self%flims(2,2)
            do hind=self%flims(1,1),self%flims(1,2)
                arg  = dot_product(shvec(:), transfmat_here(hind,kind,1,1:2))
                self%ftexp_tmpmat_im_2d(hind,kind) = aimag(self%ftexp_tmp_cmat12(hind,kind) * exp(-J * arg))
            end do
        end do
    end subroutine calc_tmpmat_im

    !< calculate tmp matrix for cost function
    subroutine calc_tmpmat_re_im( self, shvec )
        class(ftexp_shsrch), intent(inout) :: self
        real(dp),            intent(in)    :: shvec(2)
        real(dp)    :: arg
        complex(dp) :: tmp
        integer     :: hind,kind
        real, pointer :: transfmat_here(:,:,:,:)
        call self%particle%get_transfmat_ptr(transfmat_here)
        do kind=self%flims(2,1),self%flims(2,2)
            do hind=self%flims(1,1),self%flims(1,2)
                arg  = dot_product(shvec(:), transfmat_here(hind,kind,1,1:2))
                tmp  = self%ftexp_tmp_cmat12(hind,kind) * exp(-J * arg)
                self%ftexp_tmpmat_re_2d(hind,kind) = real(tmp,kind=dp)
                self%ftexp_tmpmat_im_2d(hind,kind) = aimag(tmp)
            end do
        end do
    end subroutine calc_tmpmat_re_im

    !< calculate tmp matrix for cost function
    subroutine calc_tmp_cmat12( self )
        class(ftexp_shsrch), intent(inout) :: self
        complex, pointer :: cmat1_ptr(:,:,:), cmat2_ptr(:,:,:)
        call self%reference%get_cmat_ptr(cmat1_ptr)
        call self%particle %get_cmat_ptr(cmat2_ptr)
        self%ftexp_tmp_cmat12(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2)) = &
                  cmat1_ptr(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2),1) * &
            conjg(cmat2_ptr(self%flims(1,1):self%flims(1,2),self%flims(2,1):self%flims(2,2),1))
    end subroutine calc_tmp_cmat12

    function ftexp_shsrch_corr_shifted_8( self, shvec ) result( r )
        class(ftexp_shsrch), intent(inout) :: self
        real(dp),            intent(in)    :: shvec(2)
        real(dp) :: r
        r = self%corr_shifted_cost_8( shvec )
        call self%particle%corr_normalize( self%reference, r )
    end function ftexp_shsrch_corr_shifted_8

    !> Gradient function, double precision
    subroutine ftexp_shsrch_gcost_8( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        grad = 0.d0
        select type(self)
            class is (ftexp_shsrch)
                !call self%reference%corr_gshifted_8(self%particle, -vec, grad)
                call self%corr_gshifted_cost_8( -vec, grad )
            class DEFAULT
                THROW_HARD('unknown type; ftexp_shsrch_gcost_8')
        end select
    end subroutine ftexp_shsrch_gcost_8

    !> Gradient & cost function, double precision
    subroutine ftexp_shsrch_fdfcost_8( self, vec, f, grad, D )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
        f    = 0.d0
        grad = 0.d0
        select type(self)
            class is (ftexp_shsrch)
                !call self%reference%corr_fdfshifted_8(self%particle, -vec, f, grad)
                call self%corr_fdfshifted_cost_8( -vec, f, grad )
                f = f * (-1.0_dp)
            class DEFAULT
                THROW_HARD('unknown type; ftexp_shsrch_fdfcost_8')
        end select
    end subroutine ftexp_shsrch_fdfcost_8

    subroutine test_ftexp_shsrch
        use simple_image, only: image
        type(image)       :: img_ref, img_ptcl
        type(ft_expanded) :: ftexp_ref, ftexp_ptcl
        type(ftexp_shsrch)  :: ftexp_srch
        real, parameter   :: TRS=5.0, lp=6., hp=100.
        real              :: cxy(3), x, y, lims(2,2)
        integer           :: i
        lims(:,1) = -TRS
        lims(:,2) = TRS
        call img_ref%new([100,100,1], 2.)
        call img_ptcl%new([100,100,1], 2.)
        call img_ref%square(20)
        call img_ref%fft()
        call ftexp_ref%new(img_ref, hp, lp, .true.)
        call ftexp_ptcl%new(img_ptcl, hp, lp, .true.)
        call ftexp_srch%new(ftexp_ref, ftexp_ptcl, trs)
        do i=1,100
            x = ran3()*2*TRS-TRS
            y = ran3()*2*TRS-TRS
            img_ptcl = img_ref
            call img_ptcl%shift([x,y,0.])
            call ftexp_ptcl%new(img_ptcl, hp, lp, .true.)
            cxy = ftexp_srch%minimize()
            if( cxy(1) < 0.995 )then
                THROW_HARD('shift alignment failed; test_ftexp_shsrch')
            endif
        end do
        write(logfhandle,'(a)') 'SIMPLE_ftexp_shsrch_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_ftexp_shsrch

end module simple_ftexp_shsrch
