! common strategy3D methods and type specification for polymorphic strategy3D object creation are delegated to this class
module simple_strategy3D_srch
include 'simple_lib.f08'
use simple_pftcc_shsrch_grad,  only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
use simple_polarft_corrcalc,   only: pftcc_glob, polarft_corrcalc
use simple_cartft_corrcalc,    only: cftcc_glob
use simple_cftcc_shsrch_grad,  only: cftcc_shsrch_grad
use simple_parameters,         only: params_glob
use simple_builder,            only: build_glob
use simple_strategy3D_alloc    ! singleton s3D
implicit none

public :: strategy3D_srch, strategy3D_spec, set_ptcl_stats, eval_ptcl
private
#include "simple_local_flags.inc"

type strategy3D_spec
    integer, pointer :: symmat(:,:) => null()
    integer :: iptcl=0, szsn=0
    logical :: do_extr=.false.
    real    :: extr_score_thresh=0.
end type strategy3D_spec

type strategy3D_srch
    type(pftcc_shsrch_grad) :: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
    type(cftcc_shsrch_grad) :: cart_shsrch_obj           !< origin shift search object in cartesian, L-BFGS with gradient
    type(ori)               :: o_prev                    !< previous orientation, used in continuous search
    type(oris)              :: eulspace                  !< copy of build_glob%eulspace, needed for peak detection
    integer                 :: iptcl         = 0         !< global particle index
    integer                 :: ithr          = 0         !< thread index
    integer                 :: nrefs         = 0         !< total # references (nstates*nprojs)
    integer                 :: nsample       = 0         !< # of continuous 3D rotational orientations to sample uniformly
    integer                 :: nsample_neigh = 0         !< # of continuous 3D rotational orientations to sample Gaussian
    integer                 :: nsample_trs   = 0         !< # of continuous origin shifts (2D) to sample Gaussian
    integer                 :: nstates       = 0         !< # states
    integer                 :: nprojs        = 0         !< # projections
    integer                 :: nrots         = 0         !< # in-plane rotations
    integer                 :: nsym          = 0         !< symmetry order
    integer                 :: nnn           = 0         !< # nearest neighbors
    integer                 :: nbetter       = 0         !< # better orientations identified
    integer                 :: nrefs_eval    = 0         !< # references evaluated
    integer                 :: prev_roind    = 0         !< previous in-plane rotation index
    integer                 :: prev_state    = 0         !< previous state index
    integer                 :: class         = 0         !< 2D class index
    integer                 :: prev_ref      = 0         !< previous reference index
    integer                 :: prev_proj     = 0         !< previous projection direction index
    real                    :: athres        = 15.       !< angular treshold (refine=neighc) for neighborhood continuous Cartesian search
    real                    :: prev_corr     = 1.        !< previous best correlation
    real                    :: prev_shvec(2) = 0.        !< previous origin shift vector
    logical                 :: l_neigh       = .false.   !< neighbourhood refinement flag
    logical                 :: l_greedy      = .false.   !< greedy        refinement flag
    logical                 :: l_ptclw       = .false.   !< whether to calculate particle weight
    logical                 :: doshift       = .true.    !< 2 indicate whether 2 serch shifts
    logical                 :: exists        = .false.   !< 2 indicate existence
    logical                 :: ref_only      = .false.   !< in regularization: updating the reference oris (stored in builder%reg_oris), not the particle oris
  contains
    procedure          :: new
    procedure          :: prep4srch
    procedure          :: prep4_cont_srch
    procedure, private :: inpl_srch_1
    procedure, private :: inpl_srch_2
    generic            :: inpl_srch => inpl_srch_1, inpl_srch_2
    procedure          :: store_solution
    procedure          :: shift_srch_cart
    procedure          :: kill
end type strategy3D_srch

contains

    ! class method: set_ptcl_stats for filling in stats
    ! of particles not part of update fraction
    subroutine set_ptcl_stats( pftcc, iptcl )
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        integer   :: prev_state, prev_roind, prev_proj, prev_ref
        type(ori) :: o_prev
        real      :: specscore
        prev_state = build_glob%spproj_field%get_state(iptcl)      ! state index
        if( prev_state == 0 )return
        ! previous parameters
        call build_glob%spproj_field%get_ori(iptcl, o_prev)
        prev_roind = pftcc%get_roind(360.-o_prev%e3get())          ! in-plane angle index
        prev_proj  = build_glob%eulspace%find_closest_proj(o_prev) ! previous projection direction
        prev_ref   = (prev_state-1)*params_glob%nspace + prev_proj ! previous reference
        ! calc specscore
        specscore = pftcc%specscore(prev_ref, iptcl, prev_roind)
        ! update spproj_field
        call build_glob%spproj_field%set(iptcl, 'specscore',  specscore)
        call o_prev%kill
    end subroutine set_ptcl_stats

    ! class method: evaluation of stats and objective function
    subroutine eval_ptcl( pftcc, iptcl )
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        integer   :: prev_state, prev_roind, prev_proj, prev_ref
        type(ori) :: o_prev
        real      :: specscore, corr, corrs(pftcc_glob%get_nrots())
        prev_state = build_glob%spproj_field%get_state(iptcl)      ! state index
        if( prev_state == 0 )return
        ! previous parameters
        call build_glob%spproj_field%get_ori(iptcl, o_prev)
        prev_roind = pftcc%get_roind(360.-o_prev%e3get())          ! in-plane angle index
        prev_proj  = build_glob%eulspace%find_closest_proj(o_prev) ! previous projection direction
        prev_ref   = (prev_state-1)*params_glob%nspace + prev_proj ! previous reference
        ! calc specscore
        specscore = pftcc%specscore(prev_ref, iptcl, prev_roind)
        ! prep corr
        call pftcc_glob%gencorrs(prev_ref, iptcl, corrs)
        if( params_glob%cc_objfun == OBJFUN_EUCLID )then
            corr = maxval(corrs)
        else
            corr = max(0.,maxval(corrs))
        endif
        ! update spproj_field
        call build_glob%spproj_field%set(iptcl, 'proj',      real(prev_proj))
        call build_glob%spproj_field%set(iptcl, 'corr',      corr)
        call build_glob%spproj_field%set(iptcl, 'specscore', specscore)
        call build_glob%spproj_field%set(iptcl, 'w',         1.)
        call o_prev%kill
    end subroutine eval_ptcl

    subroutine new( self, spec )
        class(strategy3D_srch), intent(inout) :: self
        class(strategy3D_spec), intent(in)    :: spec
        integer, parameter :: MAXITS = 60
        real    :: lims(2,2), lims_init(2,2)
        ! set constants
        self%iptcl         = spec%iptcl
        self%nstates       = params_glob%nstates
        self%nprojs        = params_glob%nspace
        self%nrefs         = self%nprojs*self%nstates
        self%nsample       = params_glob%nsample
        self%nsample_neigh = params_glob%nsample_neigh
        self%nsample_trs   = params_glob%nsample_trs
        self%athres        = params_glob%athres
        self%nbetter       = 0
        self%nrefs_eval    = 0
        self%nsym          = build_glob%pgrpsyms%get_nsym()
        self%doshift       = params_glob%l_doshift
        self%l_neigh       = str_has_substr(params_glob%refine, 'neigh')
        self%l_greedy      = str_has_substr(params_glob%refine, 'greedy')
        self%l_ptclw       = trim(params_glob%ptclw).eq.'yes'
        lims(:,1)          = -params_glob%trs
        lims(:,2)          =  params_glob%trs
        lims_init(:,1)     = -SHC_INPL_TRSHWDTH
        lims_init(:,2)     =  SHC_INPL_TRSHWDTH
        call self%o_prev%new(.true.)
        call self%eulspace%copy(build_glob%eulspace, is_ptcl=.false.)
        ! create in-plane search objects
        if( params_glob%l_cartesian )then
            self%nrots = params_glob%nrots
            call self%cart_shsrch_obj%new(lims, lims_init=lims_init, shbarrier=params_glob%shbarrier, maxits=MAXITS)
        else
            self%nrots = pftcc_glob%get_nrots()
            call self%grad_shsrch_obj%new(lims, lims_init=lims_init,&
                &shbarrier=params_glob%shbarrier, maxits=MAXITS, opt_angle=.true.)
        endif
        self%exists = .true.
    end subroutine new

    subroutine prep4srch( self )
        class(strategy3D_srch), intent(inout) :: self
        type(ori) :: o_prev
        real      :: corrs(self%nrots), corr
        ! previous parameters
        call build_glob%spproj_field%get_ori(self%iptcl, o_prev)            ! previous ori
        self%prev_state = o_prev%get_state()                                ! state index
        self%class      = o_prev%get_class()                                ! 2D class index
        self%prev_roind = pftcc_glob%get_roind(360.-o_prev%e3get())         ! in-plane angle index
        self%prev_shvec = o_prev%get_2Dshift()                              ! shift vector
        self%prev_proj  = build_glob%eulspace%find_closest_proj(o_prev)     ! previous projection direction
        call build_glob%spproj_field%set(self%iptcl, 'proj', real(self%prev_proj))
        self%prev_ref = (self%prev_state-1)*self%nprojs + self%prev_proj    ! previous reference
        ! init threaded search arrays
        call prep_strategy3D_thread(self%ithr)
        ! search order
        call s3D%rts(self%ithr)%ne_ran_iarr(s3D%srch_order(self%ithr,:))
        call put_last(self%prev_ref, s3D%srch_order(self%ithr,:))
        ! sanity check
        if( self%prev_state > 0 )then
            if( self%prev_state > self%nstates )          THROW_HARD('previous best state outside boundary; prep4srch')
            if( .not. s3D%state_exists(self%prev_state) ) THROW_HARD('empty previous state; prep4srch')
        endif
        ! prep corr
        call pftcc_glob%gencorrs(self%prev_ref, self%iptcl, corrs)
        if( params_glob%cc_objfun == OBJFUN_EUCLID )then
            corr = maxval(corrs)
        else
            corr = max(0.,maxval(corrs))
        endif
        self%prev_corr = corr
        call o_prev%kill
    end subroutine prep4srch

    subroutine prep4_cont_srch( self )
        class(strategy3D_srch), intent(inout) :: self
        ! previous parameters
        call build_glob%spproj_field%get_ori(self%iptcl, self%o_prev) ! previous ori
        self%prev_state = self%o_prev%get_state()                     ! state index
        self%prev_shvec = self%o_prev%get_2Dshift()                   ! shift vector
        ! sanity check
        if( self%prev_state > 0 )then
            if( self%prev_state > self%nstates )          THROW_HARD('previous best state outside boundary; prep4_cont_srch')
            if( .not. s3D%state_exists(self%prev_state) ) THROW_HARD('empty previous state; prep4_cont_srch')
        endif
        ! prep corr
        self%prev_corr = cftcc_glob%project_and_correlate(self%iptcl, self%o_prev)
        call self%o_prev%set('corr', self%prev_corr)
    end subroutine prep4_cont_srch

    function shift_srch_cart( self, nevals, shvec ) result( cxy )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(inout) :: nevals(2)
        real, optional,         intent(in)    :: shvec(2)
        real :: cxy(3)
        call self%cart_shsrch_obj%set_pind( self%iptcl )
        cxy = self%cart_shsrch_obj%minimize(nevals, shvec)
    end function shift_srch_cart

    subroutine inpl_srch_1( self )
        class(strategy3D_srch), intent(inout) :: self
        real      :: cxy(3)
        integer   :: ref, irot, loc(1)
        if( self%doshift )then
            ! BFGS over shifts with in-plane rot exhaustive callback
            loc = maxloc(s3D%proj_space_corrs(self%ithr,:))
            ref = loc(1)
            call self%grad_shsrch_obj%set_indices(ref, self%iptcl)
            cxy = self%grad_shsrch_obj%minimize(irot=irot)
            if( irot > 0 )then
                ! irot > 0 guarantees improvement found, update solution
                s3D%proj_space_euls(3,ref,self%ithr)  = 360. - pftcc_glob%get_rot(irot)
                s3D%proj_space_corrs(self%ithr,ref)   = cxy(1)
                s3D%proj_space_shift(:,ref,self%ithr) = cxy(2:3)
                s3D%proj_space_mask(ref,self%ithr)    = .true.
            endif
        endif
    end subroutine inpl_srch_1

    subroutine inpl_srch_2( self, iptcl )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(inout) :: iptcl
        type(ori) :: o
        real      :: cxy(3)
        integer   :: ref, irot
        ! class constrained refinement
        call build_glob%spproj_field%get_ori(iptcl, o)
        ref  = nint(build_glob%spproj_field%get(iptcl, 'class'))
        ! BFGS over shifts with in-plane rot exhaustive callback
        call self%grad_shsrch_obj%set_indices(ref, self%iptcl)
        cxy = self%grad_shsrch_obj%minimize(irot=irot)
        if( irot > 0 )then
            ! irot > 0 guarantees improvement found, update solution
            s3D%proj_space_euls(3,ref,self%ithr)  = 360. - pftcc_glob%get_rot(irot)
            s3D%proj_space_corrs(self%ithr,ref)   = cxy(1)
            s3D%proj_space_shift(:,ref,self%ithr) = cxy(2:3)
            s3D%proj_space_mask(ref,self%ithr)    = .true.
        else
            iptcl = 0
        endif
    end subroutine inpl_srch_2

    subroutine store_solution( self, ref, inpl_ind, corr )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: ref, inpl_ind
        real,                   intent(in)    :: corr
        s3D%proj_space_inplinds(self%ithr,ref) = inpl_ind
        s3D%proj_space_euls(3,ref,self%ithr)   = 360. - pftcc_glob%get_rot(inpl_ind)
        s3D%proj_space_corrs(self%ithr,ref)    = corr
        s3D%proj_space_mask(ref,self%ithr)     = .true.
    end subroutine store_solution

    subroutine kill( self )
        class(strategy3D_srch), intent(inout) :: self
        call self%grad_shsrch_obj%kill
        call self%o_prev%kill
        call self%eulspace%kill
    end subroutine kill

end module simple_strategy3D_srch
