! common strategy2D methods and type specification for polymorphic strategy2D object creation are delegated to this class
module simple_strategy2D_srch
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_polarft_calc,  only: pftc_glob
use simple_pftc_shsrch_grad, only: pftc_shsrch_grad ! gradient-based angle and shift search
use simple_parameters,        only: params_glob
use simple_builder,           only: build_glob
use simple_eul_prob_tab2D,    only: eul_prob_tab2D
use simple_strategy2D_alloc   ! s2D singleton
implicit none

public :: strategy2D_srch, strategy2D_spec
private

#include "simple_local_flags.inc"

type strategy2D_spec
    real    :: stoch_bound = 0.
    integer :: iptcl       = 0  ! global particle index
    integer :: iptcl_batch = 0  ! maps to index in batch
    integer :: iptcl_map   = 0  ! index in all contiguous batches
end type strategy2D_spec

type strategy2D_srch
    type(pftc_shsrch_grad) :: grad_shsrch_obj        !< origin shift search object, L-BFGS with gradient
    type(pftc_shsrch_grad) :: grad_shsrch_obj2       !< origin shift search object, L-BFGS with gradient, no call back
    type(pftc_shsrch_grad) :: grad_shsrch_first_obj  !< origin shift search object, L-BFGS with gradient, used for initial shift search on previous ref
    integer                 :: nrefs           =  0   !< number of references
    integer                 :: nrots           =  0   !< number of in-plane rotations in polar representation
    integer                 :: nrefs_eval      =  0   !< nr of references evaluated
    integer                 :: prev_class      =  0   !< previous class index
    integer                 :: best_class      =  0   !< best class index found by search
    integer                 :: best_rot        =  0   !< best in-plane rotation found by search
    integer                 :: prev_rot        =  0   !< previous in-plane rotation found by search
    integer                 :: iptcl           =  0   !< global particle index
    integer                 :: iptcl_batch     =  0   !< index in pre-allocated batch array
    integer                 :: iptcl_map       =  0   !< index in all contiguous batches combined
    integer                 :: ithr            =  0   !< current thread
    real                    :: prev_shvec(2)   =  0.  !< previous origin shift vector
    real                    :: best_shvec(2)   =  0.  !< best shift vector found by search
    real                    :: xy_first(2)     =  0.  !< initial shifts identified by searching the previous best reference
    real                    :: xy_first_rot(2) =  0.  !< initial shifts identified by searching the previous best reference, rotated
    real                    :: prev_corr       = -1.  !< previous best correlation
    real                    :: best_corr       = -1.  !< best corr found by search
    real                    :: trs             =  0.  !< shift boundary
    logical                 :: l_sh_first      = .false. !< Whether to search the shifts on previous best reference
  contains
    procedure :: new
    procedure :: prep4srch
    procedure :: inpl_srch_first
    procedure :: inpl_srch
    procedure :: store_solution
    procedure :: kill
end type strategy2D_srch

contains

    subroutine new( self, spec )
        class(strategy2D_srch), intent(inout) :: self
        class(strategy2D_spec), intent(in)    :: spec
        real :: lims(2,2), lims_init(2,2)
        call self%kill
        ! set constants
        self%iptcl       = spec%iptcl
        self%iptcl_batch = spec%iptcl_batch
        self%iptcl_map   = spec%iptcl_map
        self%nrefs       = params_glob%ncls
        self%nrots       = pftc_glob%get_nrots()
        self%nrefs_eval  = 0
        ! construct composites
        self%trs        =  params_glob%trs
        lims(:,1)       = -params_glob%trs
        lims(:,2)       =  params_glob%trs
        lims_init(:,1)  = -SHC_INPL_TRSHWDTH
        lims_init(:,2)  =  SHC_INPL_TRSHWDTH
        if( trim(params_glob%tseries).eq.'yes' )then
            ! shift only search
            call self%grad_shsrch_obj%new(lims, lims_init=lims_init,&
            maxits=params_glob%maxits_sh, opt_angle=.false.)
            call self%grad_shsrch_first_obj%new(lims, lims_init=lims_init,&
            maxits=params_glob%maxits_sh, opt_angle=.false., coarse_init=.true.)
        else
            call self%grad_shsrch_obj%new(lims, lims_init=lims_init,&
            maxits=params_glob%maxits_sh)
            call self%grad_shsrch_first_obj%new(lims, lims_init=lims_init,&
            maxits=params_glob%maxits_sh, coarse_init=.true.)
        endif
        call self%grad_shsrch_obj2%new(lims, lims_init=lims_init, maxits=params_glob%maxits_sh, opt_angle=.false.)
    end subroutine new

    subroutine prep4srch( self )
        class(strategy2D_srch), intent(inout) :: self
        real    :: corrs(pftc_glob%get_nrots())
        self%nrefs_eval = 0
        self%ithr       = omp_get_thread_num() + 1
        ! find previous discrete alignment parameters
        self%prev_class = nint(build_glob%spproj_field%get(self%iptcl,'class'))                ! class index
        self%prev_rot   = pftc_glob%get_roind(360.-build_glob%spproj_field%e3get(self%iptcl)) ! in-plane angle index
        self%prev_shvec = build_glob%spproj_field%get_2Dshift(self%iptcl)                      ! shift vector
        self%best_shvec = 0.
        if( self%prev_class > 0 )then
            if( s2D%cls_pops(self%prev_class) > 0 )then
                ! all done
            else
                ! for limiting cases
                self%prev_class = irnd_uni(self%nrefs)
                do while( s2D%cls_pops(self%prev_class) <= 0 )
                    self%prev_class = irnd_uni(self%nrefs)
                enddo
            endif
        else
            ! initialization
            self%prev_class = irnd_uni(self%nrefs)
            do while( s2D%cls_pops(self%prev_class) <= 0 )
                self%prev_class = irnd_uni(self%nrefs)
            enddo
        endif
        ! set best to previous best by default
        self%best_class = self%prev_class
        self%best_rot   = self%prev_rot
        ! calculate previous best corr (treshold for better)
        call pftc_glob%gen_objfun_vals(self%prev_class, self%iptcl, [0.,0.], corrs)
        if( params_glob%cc_objfun == OBJFUN_CC )then
            self%prev_corr  = max(0., corrs(self%prev_rot))
        else
            self%prev_corr  = corrs(self%prev_rot)
        endif
        self%best_corr = self%prev_corr
        ! whether to search shifts first
        self%l_sh_first   = s2D%do_inplsrch(self%iptcl_batch) .and. params_glob%l_sh_first
        self%xy_first     =  0.
        self%xy_first_rot =  0.
    end subroutine prep4srch

    subroutine inpl_srch_first( self )
        class(strategy2D_srch), intent(inout) :: self
        real    :: cxy(3), rotmat(2,2)
        integer :: irot
        self%best_shvec = [0.,0.]
        if( .not. self%l_sh_first ) return
        ! BFGS
        irot = 0
        call self%grad_shsrch_first_obj%set_indices(self%prev_class, self%iptcl)
        if( .not.self%grad_shsrch_first_obj%does_opt_angle() )then
            ! shift-only optimization
            irot = self%prev_rot
        endif
        cxy = self%grad_shsrch_first_obj%minimize(irot=irot, sh_rot=.false.)
        if( irot == 0 ) cxy(2:3) = 0.
        self%xy_first = cxy(2:3)
        self%xy_first_rot = 0.
        if( irot > 0 )then
            ! rotate the shift vector to the frame of reference
            call rotmat2d(pftc_glob%get_rot(irot), rotmat)
            self%xy_first_rot = matmul(cxy(2:3), rotmat)
            ! update best
            self%best_corr  = cxy(1)
            self%best_rot   = irot
            self%best_shvec = self%xy_first_rot
        endif
    end subroutine inpl_srch_first

    subroutine inpl_srch( self )
        class(strategy2D_srch), intent(inout) :: self
        real    :: cxy(3)
        integer :: irot
        irot = 0
        self%best_shvec = [0.,0.]
        if( s2D%do_inplsrch(self%iptcl_batch) )then
            ! BFGS
            call self%grad_shsrch_obj%set_indices(self%best_class, self%iptcl)
            if( .not.self%grad_shsrch_obj%does_opt_angle() )then
                ! shift-only optimization
                irot = self%best_rot
            endif
            if( self%l_sh_first )then
                cxy = self%grad_shsrch_obj%minimize(irot=irot, xy_in=self%xy_first)
            else
                cxy = self%grad_shsrch_obj%minimize(irot=irot)
            endif
            if( irot > 0 )then
                self%best_corr  = cxy(1)
                self%best_rot   = irot
                self%best_shvec = cxy(2:3)
            endif
        endif
    end subroutine inpl_srch

    subroutine store_solution( self, nrefs, w_in )
        class(strategy2D_srch), intent(in) :: self
        integer,      optional, intent(in) :: nrefs
        real,         optional, intent(in) :: w_in
        real :: dist, mat(2,2), u(2), x1(2), x2(2)
        real :: e3, mi_class, frac, w
        ! get in-plane angle
        e3   = 360. - pftc_glob%get_rot(self%best_rot) ! change sgn to fit convention
        ! calculate in-plane rot dist (radians)
        u(1) = 0.
        u(2) = 1.
        call rotmat2d(e3, mat)
        x1   = matmul(u,mat)
        call rotmat2d(build_glob%spproj_field%e3get(self%iptcl), mat)
        x2   = matmul(u,mat)
        dist = myacos(dot_product(x1,x2))
        ! calculate overlap between distributions
        mi_class = 0.
        if( self%prev_class == self%best_class ) mi_class = 1.
        ! search space explored
        if( present(nrefs) )then
            frac = 100.*(real(self%nrefs_eval)/real(nrefs))
        else
            frac = 100.*(real(self%nrefs_eval)/real(self%nrefs))
        endif
        ! weight
        w = 1.0
        if( present(w_in) ) w = w_in
        ! update parameters
        call build_glob%spproj_field%e3set(self%iptcl,e3)
        call build_glob%spproj_field%set_shift(self%iptcl, self%prev_shvec + self%best_shvec)
        call build_glob%spproj_field%set(self%iptcl, 'shincarg',   arg(self%best_shvec))
        call build_glob%spproj_field%set(self%iptcl, 'inpl',       real(self%best_rot))
        call build_glob%spproj_field%set(self%iptcl, 'class',      real(self%best_class))
        call build_glob%spproj_field%set(self%iptcl, 'corr',       self%best_corr)
        call build_glob%spproj_field%set(self%iptcl, 'dist_inpl',  rad2deg(dist))
        call build_glob%spproj_field%set(self%iptcl, 'mi_class',   mi_class)
        call build_glob%spproj_field%set(self%iptcl, 'frac',       frac)
        call build_glob%spproj_field%set(self%iptcl, 'w',          w)
    end subroutine store_solution

    subroutine kill( self )
        class(strategy2D_srch), intent(inout) :: self
        call self%grad_shsrch_obj%kill
        call self%grad_shsrch_obj2%kill
        call self%grad_shsrch_first_obj%kill
    end subroutine kill

end module simple_strategy2D_srch
