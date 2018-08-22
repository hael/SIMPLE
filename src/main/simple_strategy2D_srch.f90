! common strategy2D methods and type specification for polymorphic strategy2D object creation are delegated to this class
module simple_strategy2D_srch
include 'simple_lib.f08'
use simple_polarft_corrcalc,  only: pftcc_glob
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad ! gradient-based angle and shift search
use simple_oris,              only: oris
use simple_parameters,        only: params_glob
use simple_builder,           only: build_glob
use simple_strategy2D_alloc   ! s2D singleton
implicit none

public :: strategy2D_srch, strategy2D_spec
private

#include "simple_local_flags.inc"

type strategy2D_spec
    real    :: stoch_bound = 0.
    integer :: iptcl       = 0
    integer :: iptcl_map   = 0
end type strategy2D_spec

type strategy2D_srch
    type(pftcc_shsrch_grad) :: grad_shsrch_obj      !< origin shift search object, L-BFGS with gradient
    integer                 :: nrefs         =  0   !< number of references
    integer                 :: nrots         =  0   !< number of in-plane rotations in polar representation
    integer                 :: nrefs_eval    =  0   !< nr of references evaluated
    integer                 :: prev_class    =  0   !< previous class index
    integer                 :: best_class    =  0   !< best class index found by search
    integer                 :: prev_rot      =  0   !< previous in-plane rotation index
    integer                 :: best_rot      =  0   !< best in-plane rotation found by search
    integer                 :: nthr          =  0   !< number of threads
    integer                 :: fromp         =  1   !< from particle index
    integer                 :: top           =  1   !< to particle index
    integer                 :: nnn           =  0   !< # nearest neighbors
    integer                 :: iptcl         =  0   !< global particle index
    integer                 :: iptcl_map     =  0   !< index in pre-allocated arrays
    real                    :: trs           =  0.  !< shift range parameter [-trs,trs]
    real                    :: prev_shvec(2) =  0.  !< previous origin shift vector
    real                    :: best_shvec(2) =  0.  !< best shift vector found by search
    real                    :: prev_corr     = -1.  !< previous best correlation
    real                    :: best_corr     = -1.  !< best corr found by search
    real                    :: prev_bfac     =  0.  !< previous b-factor
    real                    :: specscore     =  0.  !< spectral score
    logical                 :: doshift    = .true.  !< origin shift search indicator
  contains
    procedure :: new
    procedure :: prep4srch
    procedure :: inpl_srch
    procedure :: store_solution
    procedure :: kill
end type strategy2D_srch

contains

    subroutine new( self, spec )
        class(strategy2D_srch), intent(inout) :: self
        class(strategy2D_spec), intent(in)    :: spec
        integer, parameter :: MAXITS = 60
        real :: lims(2,2), lims_init(2,2)
        ! set constants
        self%iptcl      =  spec%iptcl
        self%iptcl_map  =  spec%iptcl_map
        self%nrefs      =  params_glob%ncls
        self%nrots      =  round2even(twopi*real(params_glob%ring2))
        self%nrefs_eval =  0
        self%trs        =  params_glob%trs
        self%doshift    =  params_glob%l_doshift
        self%nthr       =  params_glob%nthr
        self%fromp      =  params_glob%fromp
        self%top        =  params_glob%top
        self%nnn        =  params_glob%nnn
        ! construct composites
        lims(:,1)       = -params_glob%trs
        lims(:,2)       =  params_glob%trs
        lims_init(:,1)  = -SHC_INPL_TRSHWDTH
        lims_init(:,2)  =  SHC_INPL_TRSHWDTH
        call self%grad_shsrch_obj%new(lims, lims_init=lims_init, maxits=MAXITS)
     !   DebugPrint  '>>> strategy2D_srch::CONSTRUCTED NEW SIMPLE_strategy2D_srch OBJECT'
    end subroutine new

    subroutine prep4srch( self )
        class(strategy2D_srch), intent(inout) :: self
        real :: corrs(pftcc_glob%get_nrots())
        self%nrefs_eval = 0
        ! find previous discrete alignment parameters
        self%prev_class = nint(build_glob%spproj_field%get(self%iptcl,'class'))                ! class index
        self%prev_rot   = pftcc_glob%get_roind(360.-build_glob%spproj_field%e3get(self%iptcl)) ! in-plane angle index
        self%prev_shvec = build_glob%spproj_field%get_2Dshift(self%iptcl)                      ! shift vector
        ! set best to previous best by default
        self%best_class = self%prev_class
        self%best_rot   = self%prev_rot
        ! calculate previous best corr (treshold for better) & b-factor
        if( self%prev_class > 0 )then
            self%prev_bfac = pftcc_glob%fit_bfac(self%prev_class, self%iptcl, self%prev_rot, [0.,0.])
            if(params_glob%cc_objfun == OBJFUN_RES)then
                ! prior to correlation calculation
                call pftcc_glob%memoize_bfac(self%iptcl, self%prev_bfac)
            endif
            call pftcc_glob%gencorrs(self%prev_class, self%iptcl, corrs)
            self%prev_corr  = max(0., corrs(self%prev_rot))
            self%best_corr  = self%prev_corr
        else
            self%prev_class = irnd_uni(self%nrefs)
            self%prev_corr  = 0.
            self%best_corr  = 0.
            self%prev_bfac  = 500. ! default b-factor
            if(params_glob%cc_objfun == OBJFUN_RES)then
                call pftcc_glob%memoize_bfac(self%iptcl, self%prev_bfac)
            endif
        endif
        call build_glob%spproj_field%set(self%iptcl, 'bfac', self%prev_bfac)
        ! calculate spectral score
        self%specscore = pftcc_glob%specscore(self%prev_class, self%iptcl, self%prev_rot)
        if( DEBUG ) print *, '>>> strategy2D_srch::PREPARED FOR SIMPLE_strategy2D_srch'
    end subroutine prep4srch

    subroutine inpl_srch( self )
        class(strategy2D_srch), intent(inout) :: self
        real              :: cxy(3)
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
      !  DebugPrint '>>> strategy2D_srch::FINISHED SHIFT SEARCH'
    end subroutine inpl_srch

    subroutine store_solution( self, nrefs )
        use simple_ori,  only: ori
        class(strategy2D_srch), intent(in) :: self
        integer,      optional, intent(in) :: nrefs
        real :: dist, mat(2,2), u(2), x1(2), x2(2)
        real :: e3, mi_class, mi_inpl, mi_joint, frac
        ! get in-plane angle
        e3   = 360. - pftcc_glob%get_rot(self%best_rot) ! change sgn to fit convention
        ! calculate in-plane rot dist (radians)
        u(1) = 0.
        u(2) = 1.
        mat  = rotmat2d(e3)
        x1   = matmul(u,mat)
        mat  = rotmat2d(build_glob%spproj_field%e3get(self%iptcl))
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
        ! search psace explored
        if( present(nrefs) )then
            frac = 100.*(real(self%nrefs_eval)/real(nrefs))
        else
            frac = 100.*(real(self%nrefs_eval)/real(self%nrefs))
        endif
        ! update parameters
        call build_glob%spproj_field%e3set(self%iptcl,e3)
        call build_glob%spproj_field%set_shift(self%iptcl, self%prev_shvec + self%best_shvec)
        call build_glob%spproj_field%set(self%iptcl, 'inpl',       real(self%best_rot))
        call build_glob%spproj_field%set(self%iptcl, 'class',      real(self%best_class))
        call build_glob%spproj_field%set(self%iptcl, 'corr',       self%best_corr)
        call build_glob%spproj_field%set(self%iptcl, 'specscore',  self%specscore)
        call build_glob%spproj_field%set(self%iptcl, 'dist_inpl',  rad2deg(dist))
        call build_glob%spproj_field%set(self%iptcl, 'mi_class',   mi_class)
        call build_glob%spproj_field%set(self%iptcl, 'mi_inpl',    mi_inpl)
        call build_glob%spproj_field%set(self%iptcl, 'mi_joint',   mi_joint)
        call build_glob%spproj_field%set(self%iptcl, 'frac',       frac)
      !  DebugPrint  '>>> strategy2D_srch::GOT BEST ORI'
    end subroutine store_solution

    subroutine kill( self )
        class(strategy2D_srch),  intent(inout) :: self !< instance
        call self%grad_shsrch_obj%kill
    end subroutine kill

end module simple_strategy2D_srch
