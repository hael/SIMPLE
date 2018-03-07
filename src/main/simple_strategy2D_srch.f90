! common strategy2D methods and type specification for polymorphic strategy2D object creation are delegated to this class
module simple_strategy2D_srch
#include "simple_lib.f08"
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad ! gradient-based angle and shift search
use simple_oris,              only: oris
use simple_params,            only: params
use simple_strategy2D_alloc   ! use all in there
use simple_timer              ! use all in there
implicit none

public :: strategy2D_srch, strategy2D_spec
private

logical, parameter :: DEBUG   = .false.

type strategy2D_spec
    class(params),           pointer :: pp         => null()
    class(polarft_corrcalc), pointer :: ppftcc     => null()
    class(oris),             pointer :: pa         => null()
    integer,                 pointer :: nnmat(:,:) => null()
    integer :: iptcl=0, iptcl_map=0
    logical :: do_extr=.false.
    real    :: corr_bound=0.
end type strategy2D_spec

type strategy2D_srch
    class(polarft_corrcalc), pointer :: pftcc_ptr => null()  !< pointer to pftcc (corrcalc) object
    class(oris),             pointer :: a_ptr     => null()  !< pointer to b%a (primary particle orientation table)
    type(pftcc_shsrch_grad)          :: grad_shsrch_obj      !< origin shift search object, L-BFGS with gradient
    integer                          :: nrefs         =  0   !< number of references
    integer                          :: nrots         =  0   !< number of in-plane rotations in polar representation
    integer                          :: nrefs_eval    =  0   !< nr of references evaluated
    integer                          :: prev_class    =  0   !< previous class index
    integer                          :: best_class    =  0   !< best class index found by search
    integer                          :: prev_rot      =  0   !< previous in-plane rotation index
    integer                          :: best_rot      =  0   !< best in-plane rotation found by search
    integer                          :: nthr          =  0   !< number of threads
    integer                          :: fromp         =  1   !< from particle index
    integer                          :: top           =  1   !< to particle index
    integer                          :: nnn           =  0   !< # nearest neighbors
    integer                          :: iptcl         =  0   !< global particle index
    integer                          :: iptcl_map     =  0   !< index in pre-allocated arrays
    real                             :: trs           =  0.  !< shift range parameter [-trs,trs]
    real                             :: prev_shvec(2) =  0.  !< previous origin shift vector
    real                             :: best_shvec(2) =  0.  !< best shift vector found by search
    real                             :: prev_corr     = -1.  !< previous best correlation
    real                             :: best_corr     = -1.  !< best corr found by search
    real                             :: specscore     =  0.  !< spectral score
    logical                          :: dyncls  = .true.     !< whether to turn on dynamic class update (use of low population threshold)
    logical                          :: doshift = .true.     !< origin shift search indicator
  contains
    procedure :: new
    procedure :: prep4srch
    procedure :: inpl_srch
    procedure :: calc_corr
    procedure :: fit_bfac
    procedure :: store_solution
    procedure :: kill
end type strategy2D_srch

contains

    subroutine new( self, spec )
        use simple_params, only: params
        class(strategy2D_srch), intent(inout) :: self
        class(strategy2D_spec), intent(in)    :: spec
        integer, parameter :: MAXITS = 60
        real :: lims(2,2), lims_init(2,2)
        ! set constants
        self%pftcc_ptr  => spec%ppftcc
        self%a_ptr      => spec%pa
        self%iptcl      =  spec%iptcl
        self%iptcl_map  =  spec%iptcl_map
        self%nrefs      =  spec%pp%ncls
        self%nrots      =  round2even(twopi*real(spec%pp%ring2))
        self%nrefs_eval =  0
        self%trs        =  spec%pp%trs
        self%doshift    =  spec%pp%l_doshift
        self%nthr       =  spec%pp%nthr
        self%fromp      =  spec%pp%fromp
        self%top        =  spec%pp%top
        self%nnn        =  spec%pp%nnn
        self%dyncls     =  (spec%pp%dyncls.eq.'yes')
        ! construct composites
        lims(:,1)       = -spec%pp%trs
        lims(:,2)       =  spec%pp%trs
        lims_init(:,1)  = -SHC_INPL_TRSHWDTH
        lims_init(:,2)  =  SHC_INPL_TRSHWDTH
        call self%grad_shsrch_obj%new(self%pftcc_ptr, lims, lims_init=lims_init, maxits=MAXITS)
        if( DEBUG ) print *, '>>> strategy2D_srch::CONSTRUCTED NEW SIMPLE_strategy2D_srch OBJECT'
    end subroutine new

    subroutine prep4srch( self )
        class(strategy2D_srch), intent(inout) :: self
        real :: corrs(self%pftcc_ptr%get_nrots()), bfac
        ! find previous discrete alignment parameters
        self%prev_class = nint(self%a_ptr%get(self%iptcl,'class')) ! class index
        if( self%dyncls )then
            if( self%prev_class > 0 )then
                ! reassignement to a class with higher population
                do while( cls_pops(self%prev_class) <= MINCLSPOPLIM )
                   self%prev_class = irnd_uni(self%nrefs)
                enddo
            endif
        endif
        self%prev_rot   = self%pftcc_ptr%get_roind(360.-self%a_ptr%e3get(self%iptcl))     ! in-plane angle index
        self%prev_shvec = [self%a_ptr%get(self%iptcl,'x'),self%a_ptr%get(self%iptcl,'y')] ! shift vector
        ! set best to previous best by default
        self%best_class = self%prev_class
        self%best_rot   = self%prev_rot
        ! calculate previous best corr (treshold for better)
        if( self%prev_class > 0 )then
            call self%pftcc_ptr%gencorrs(self%prev_class, self%iptcl, corrs)
            self%prev_corr  = max(0., corrs(self%prev_rot))
            self%best_corr  = self%prev_corr
        else
            self%prev_class = irnd_uni(self%nrefs)
            self%prev_corr  = 0.
            self%best_corr  = 0.
        endif
        ! calculate spectral score
        self%specscore = self%pftcc_ptr%specscore(self%prev_class, self%iptcl, self%prev_rot)
        ! B-factor memoization
        if( self%pftcc_ptr%objfun_is_ccres() )then
            if( .not. self%a_ptr%isthere(self%iptcl, 'bfac') )then
                bfac = self%pftcc_ptr%fit_bfac(self%prev_class, self%iptcl, self%prev_rot, [0.,0.])
                call self%pftcc_ptr%memoize_bfac(self%iptcl, bfac)
            endif
        endif
        if( DEBUG ) print *, '>>> strategy2D_srch::PREPARED FOR SIMPLE_strategy2D_srch'
    end subroutine prep4srch

    subroutine inpl_srch( self )
        class(strategy2D_srch), intent(inout) :: self
        real, allocatable :: cxy(:)
        integer           :: irot
        self%best_shvec = [0.,0.]
        if( self%doshift )then
            ! BFGS
            call self%grad_shsrch_obj%set_indices(self%best_class, self%iptcl)
            cxy = self%grad_shsrch_obj%minimize(irot=irot)
            if( irot > 0 )then
                self%best_corr  = cxy(1)
                self%best_rot   = irot
                self%best_shvec = cxy(2:3)
            endif
        endif
        if( DEBUG ) write(*,'(A)') '>>> strategy2D_srch::FINISHED SHIFT SEARCH'
    end subroutine inpl_srch

    subroutine calc_corr( self )
        class(strategy2D_srch),   intent(inout) :: self
        integer :: iref, prev_roind, state
        real    :: cc
        state = self%a_ptr%get_state(self%iptcl)
        if( state > 0 )then
            self%prev_class = nint(self%a_ptr%get(self%iptcl,'class'))                   ! class index
            if( self%prev_class > 0 )then
                prev_roind = self%pftcc_ptr%get_roind(360.-self%a_ptr%e3get(self%iptcl)) ! rotation index
                cc = self%pftcc_ptr%gencorr_cc_for_rot(self%prev_class, self%iptcl, [0.,0.], prev_roind)
                call self%a_ptr%set(self%iptcl,'corr', cc)
            endif
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        if( DEBUG ) print *, '>>> strategy2D_srch::FINISHED CALC_CORR'
    end subroutine calc_corr

    subroutine fit_bfac( self )
        class(strategy2D_srch), intent(inout) :: self
        real :: bfac
        if( self%pftcc_ptr%objfun_is_ccres() )then
            bfac = self%pftcc_ptr%fit_bfac(self%best_class, self%iptcl, self%best_rot, self%best_shvec)
            call self%a_ptr%set(self%iptcl, 'bfac',  bfac )
        endif
    end subroutine fit_bfac

    subroutine store_solution( self )
        use simple_ori,  only: ori
        class(strategy2D_srch), intent(in) :: self
        real :: dist, mat(2,2), u(2), x1(2), x2(2)
        real :: e3, mi_class, mi_inpl, mi_joint
        ! get in-plane angle
        e3   = 360. - self%pftcc_ptr%get_rot(self%best_rot) ! change sgn to fit convention
        ! calculate in-plane rot dist (radians)
        u(1) = 0.
        u(2) = 1.
        mat  = rotmat2d(e3)
        x1   = matmul(u,mat)
        mat  = rotmat2d(self%a_ptr%e3get(self%iptcl))
        x2   = matmul(u,mat)
        dist = myacos(dot_product(x1,x2))
        ! calculate overlap between distributions
        mi_class = 0.
        mi_inpl  = 0.
        mi_joint = 0.
        if( self%prev_class == self%best_class )then
            mi_class = mi_class + 1.
            mi_joint = mi_joint + 1.
        endif
        if( self%prev_rot ==  self%best_rot )then
            mi_inpl  = mi_inpl  + 1.
            mi_joint = mi_joint + 1.
        endif
        mi_joint = mi_joint / 2.
        ! update parameters
        call self%a_ptr%e3set(self%iptcl,e3)
        call self%a_ptr%set_shift(self%iptcl, self%prev_shvec + self%best_shvec)
        call self%a_ptr%set(self%iptcl, 'inpl',       real(self%best_rot))
        call self%a_ptr%set(self%iptcl, 'class',      real(self%best_class))
        call self%a_ptr%set(self%iptcl, 'corr',       self%best_corr)
        call self%a_ptr%set(self%iptcl, 'specscore',  self%specscore)
        call self%a_ptr%set(self%iptcl, 'dist_inpl',  rad2deg(dist))
        call self%a_ptr%set(self%iptcl, 'mi_class',   mi_class)
        call self%a_ptr%set(self%iptcl, 'mi_inpl',    mi_inpl)
        call self%a_ptr%set(self%iptcl, 'mi_joint',   mi_joint)
        call self%a_ptr%set(self%iptcl, 'frac',       100.*(real(self%nrefs_eval)/real(self%nrefs)))
        if( DEBUG ) print *, '>>> strategy2D_srch::GOT BEST ORI'
    end subroutine store_solution

    subroutine kill( self )
        class(strategy2D_srch),  intent(inout) :: self !< instance
        call self%grad_shsrch_obj%kill
    end subroutine kill

end module simple_strategy2D_srch
