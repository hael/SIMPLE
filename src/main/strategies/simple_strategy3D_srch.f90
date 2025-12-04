! common strategy3D methods and type specification for polymorphic strategy3D object creation are delegated to this class
module simple_strategy3D_srch
include 'simple_lib.f08'
use simple_pftc_shsrch_grad,  only: pftc_shsrch_grad  ! gradient-based in-plane angle and shift search
use simple_polarft_calc,   only: pftc_glob, polarft_calc
use simple_parameters,         only: params_glob
use simple_builder,            only: build_glob
use simple_eul_prob_tab,       only: eul_prob_tab
use simple_strategy3D_alloc    ! singleton s3D
implicit none

public :: strategy3D_srch, strategy3D_spec
private
#include "simple_local_flags.inc"

type strategy3D_spec
    type(eul_prob_tab), pointer :: eulprob_obj_part
    integer :: iptcl=0, iptcl_map=0
end type strategy3D_spec

type strategy3D_srch
    character(len=:), allocatable :: refine                !< 3D refinement flag
    type(pftc_shsrch_grad) :: grad_shsrch_obj             !< origin shift search object, L-BFGS with gradient
    type(pftc_shsrch_grad) :: grad_shsrch_obj2            !<
    type(pftc_shsrch_grad) :: grad_shsrch_first_obj       !< origin shift search object, L-BFGS with gradient, used for initial shift search on previous ref
    type(ori)               :: o_prev                      !< previous orientation
    type(oris)              :: opeaks                      !< peak orientations to consider for refinement
    integer                 :: iptcl           = 0         !< global particle index
    integer                 :: iptcl_map       = 0         !< map particle index
    integer                 :: ithr            = 0         !< thread index
    integer                 :: nrefs           = 0         !< total # references (nstates*nprojs)
    integer                 :: nrefs_sub       = 0         !< total # references (nstates*nprojs), subspace
    integer                 :: npeaks          = 0         !< # peak subspace orientations to consider
    integer                 :: npeaks_inpl     = 0         !< # # multi-neighborhood peaks to refine with L-BFGS
    integer                 :: nstates         = 0         !< # states
    integer                 :: nprojs          = 0         !< # projections
    integer                 :: nprojs_sub      = 0         !< # projections, subspace
    integer                 :: nrots           = 0         !< # in-plane rotations
    integer                 :: nsym            = 0         !< symmetry order
    integer                 :: nnn             = 0         !< # nearest neighbors
    integer                 :: nbetter         = 0         !< # better orientations identified
    integer                 :: nrefs_eval      = 0         !< # references evaluated
    integer                 :: ntrs_eval       = 0         !< # shifts evaluated
    integer                 :: prev_roind      = 0         !< previous in-plane rotation index
    integer                 :: prev_state      = 0         !< previous state index
    integer                 :: prev_ref        = 0         !< previous reference index
    integer                 :: prev_proj       = 0         !< previous projection direction index
    real                    :: athres          = 10.       !< angular treshold (refine=neighc) for neighborhood continuous Cartesian search
    real                    :: prev_corr       = 1.        !< previous best correlation
    real                    :: prev_shvec(2)   = 0.        !< previous origin shift vector
    real                    :: xy_first(2)     = 0.        !< initial shifts identified by searching the previous best reference
    real                    :: xy_first_rot(2) = 0.        !< initial shifts identified by searching the previous best reference, rotated
    logical                 :: l_neigh         = .false.   !< neighbourhood refinement flag
    logical                 :: l_greedy        = .false.   !< greedy        refinement flag
    logical                 :: doshift         = .true.    !< 2 indicate whether 2 serch shifts
    logical                 :: exists          = .false.   !< 2 indicate existence
  contains
    procedure :: new
    procedure :: prep4srch
    procedure :: inpl_srch_first
    procedure :: inpl_srch
    procedure :: inpl_srch_peaks
    procedure :: store_solution
    procedure :: kill
end type strategy3D_srch

contains

    subroutine new( self, spec )
        class(strategy3D_srch), intent(inout) :: self
        class(strategy3D_spec), intent(in)    :: spec
        real :: lims(2,2), lims_init(2,2)
        ! set constants
        self%iptcl         = spec%iptcl
        self%iptcl_map     = spec%iptcl_map
        self%nstates       = params_glob%nstates
        self%nprojs        = params_glob%nspace
        self%nprojs_sub    = params_glob%nspace_sub
        self%nrefs         = self%nprojs     * self%nstates
        self%nrefs_sub     = self%nprojs_sub * self%nstates
        self%npeaks        = params_glob%npeaks
        self%npeaks_inpl   = params_glob%npeaks_inpl
        self%athres        = params_glob%athres
        self%nbetter       = 0
        self%nrefs_eval    = 0
        self%ntrs_eval     = 0
        self%nsym          = build_glob%pgrpsyms%get_nsym()
        self%doshift       = params_glob%l_doshift
        self%l_neigh       = params_glob%l_neigh
        self%refine        = trim(params_glob%refine)
        self%l_greedy      = str_has_substr(params_glob%refine, 'greedy')
        lims(:,1)          = -params_glob%trs
        lims(:,2)          =  params_glob%trs
        lims_init(:,1)     = -SHC_INPL_TRSHWDTH
        lims_init(:,2)     =  SHC_INPL_TRSHWDTH
        call self%o_prev%new(is_ptcl=.true.)
        call self%opeaks%new(self%npeaks, is_ptcl=.true.)
        ! create in-plane search objects
        self%nrots = pftc_glob%get_nrots()
        call self%grad_shsrch_obj%new(lims, lims_init=lims_init, shbarrier=params_glob%shbarrier,&
        &maxits=params_glob%maxits_sh, opt_angle=.true.)
        call self%grad_shsrch_obj2%new(lims, lims_init=lims_init, maxits=params_glob%maxits_sh,&
        &shbarrier=params_glob%shbarrier, opt_angle=.false.)
        call self%grad_shsrch_first_obj%new(lims, lims_init=lims_init, shbarrier=params_glob%shbarrier,&
        &maxits=params_glob%maxits_sh, opt_angle=.true., coarse_init=.true.)
        self%exists = .true.
    end subroutine new

    subroutine prep4srch( self )
        class(strategy3D_srch), intent(inout) :: self
        real      :: corrs(self%nrots), corr
        integer   :: ipeak, tmp_inds(self%nrefs_sub), iref_sub, prev_proj_sub
        ! previous parameters
        call build_glob%spproj_field%get_ori(self%iptcl, self%o_prev)        ! previous ori
        self%prev_state = self%o_prev%get_state()                            ! state index
        self%prev_roind = pftc_glob%get_roind(360.-self%o_prev%e3get())     ! in-plane angle index
        self%prev_shvec = self%o_prev%get_2Dshift()                          ! shift vector
        self%prev_proj  = build_glob%eulspace%find_closest_proj(self%o_prev) ! previous projection direction
        self%prev_ref   = (self%prev_state-1)*self%nprojs + self%prev_proj   ! previous reference
        call build_glob%spproj_field%set(self%iptcl, 'proj', real(self%prev_proj))
        ! copy self%o_prev to opeaks to transfer paticle-dependent parameters
        do ipeak = 1, self%npeaks
            call self%opeaks%set_ori(ipeak, self%o_prev)
        end do
        ! init threaded search arrays
        call prep_strategy3D_thread(self%ithr)
        ! search order
        ! -- > full space
        call s3D%rts(     self%ithr)%ne_ran_iarr(s3D%srch_order(:,self%ithr))
        call s3D%rts_inpl(self%ithr)%ne_ran_iarr(s3D%inpl_order(:,self%ithr))
        call put_last(self%prev_ref, s3D%srch_order(:,self%ithr))
        ! --> subspace
        if( self%l_neigh )then
            call s3D%rts_sub(self%ithr)%ne_ran_iarr(tmp_inds)
            ! --> do the mapping
            do iref_sub = 1, self%nrefs_sub
                ! index for build_glob%subspace_inds needs to be a projection direction
                prev_proj_sub = build_glob%subspace_inds(&
                &tmp_inds(iref_sub) - (self%prev_state - 1) * params_glob%nspace_sub)
                ! but then we need to turn it back into a reference index in the full search space
                s3D%srch_order_sub(iref_sub,self%ithr) = (self%prev_state - 1) * self%nprojs + prev_proj_sub
            enddo
            ! --> check if we have prev_ref in subspace, put last
            if( any(s3D%srch_order_sub(:,self%ithr) == self%prev_ref ) )then
                call put_last(self%prev_ref, s3D%srch_order_sub(:,self%ithr))
            endif
        endif
        ! sanity check
        if( self%prev_state > 0 )then
            if( self%prev_state > self%nstates )          THROW_HARD('previous best state outside boundary; prep4srch')
            if( .not. s3D%state_exists(self%prev_state) ) THROW_HARD('empty previous state; prep4srch')
        endif
        ! prep corr
        call pftc_glob%gen_objfun_vals(self%prev_ref, self%iptcl, [0.,0.], corrs)
        corr = max(0.,maxval(corrs))
        self%prev_corr = corr
    end subroutine prep4srch

    subroutine inpl_srch_first( self )
        class(strategy3D_srch), intent(inout) :: self
        real    :: cxy(3), rotmat(2,2)
        integer :: irot
        if( .not. self%doshift           ) return
        if( .not. params_glob%l_sh_first ) return
        ! BFGS over shifts with in-plane rot exhaustive callback
        irot = self%prev_roind
        call self%grad_shsrch_first_obj%set_indices(self%prev_ref, self%iptcl)
        cxy = self%grad_shsrch_first_obj%minimize(irot=irot, sh_rot=.false.)
        if( irot == 0 ) cxy(2:3) = 0.
        self%xy_first = cxy(2:3)
        self%xy_first_rot = 0.
        if( irot > 0 )then
            ! rotate the shift vector to the frame of reference
            call rotmat2d(pftc_glob%get_rot(irot), rotmat)
            self%xy_first_rot = matmul(cxy(2:3), rotmat)
        endif
    end subroutine inpl_srch_first

    subroutine inpl_srch( self, ref, irot_in )
        class(strategy3D_srch), intent(inout) :: self
        integer, optional,      intent(in)    :: ref
        integer, optional,      intent(inout) :: irot_in
        real    :: cxy(3)
        integer :: iref, irot, loc(1)
        if( .not. self%doshift ) return
        if( present(ref) )then
            iref = ref
        else
            loc  = maxloc(s3D%proj_space_corrs(:,self%ithr))
            iref = loc(1)
        endif
        irot = s3D%proj_space_inplinds(iref,self%ithr)
        if( present(irot_in) ) irot = irot_in
        ! BFGS over shifts with in-plane rot exhaustive callback
        call self%grad_shsrch_obj%set_indices(iref, self%iptcl)
        if( params_glob%l_sh_first )then
            cxy = self%grad_shsrch_obj%minimize(irot=irot, xy_in=self%xy_first)
        else
            cxy = self%grad_shsrch_obj%minimize(irot=irot)
        endif
        if( irot > 0 ) call self%store_solution(iref, irot, cxy(1), sh=cxy(2:3))
        if( present(irot_in) ) irot_in = irot
    end subroutine inpl_srch

    subroutine inpl_srch_peaks( self )
        class(strategy3D_srch), intent(inout) :: self
        real    :: cxy(3)
        integer :: refs(self%npeaks_inpl), irot, ipeak
        if( .not. self%doshift ) return
        ! BFGS over shifts with in-plane rot exhaustive callback
        refs = maxnloc(s3D%proj_space_corrs(:,self%ithr), self%npeaks_inpl)
        do ipeak = 1, self%npeaks_inpl
            call self%grad_shsrch_obj%set_indices(refs(ipeak), self%iptcl)
            if( params_glob%l_sh_first )then
                cxy = self%grad_shsrch_obj%minimize(irot=irot, xy_in=self%xy_first)
            else
                cxy = self%grad_shsrch_obj%minimize(irot=irot)
            endif
            if( irot > 0 )then
                ! irot > 0 guarantees improvement found, update solution
                call self%store_solution(refs(ipeak), irot, cxy(1), sh=cxy(2:3))
            endif
        enddo
    end subroutine inpl_srch_peaks

    subroutine store_solution( self, ref, inpl_ind, corr, sh, w )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: ref, inpl_ind
        real,                   intent(in)    :: corr
        real,       optional,   intent(in)    :: sh(2)
        real,       optional,   intent(in)    :: w
        if( present(sh) ) s3D%proj_space_shift(:,ref,self%ithr) = sh
        if( present(w)  ) s3D%proj_space_w(      ref,self%ithr) = w
        s3D%proj_space_inplinds(ref,self%ithr) = inpl_ind
        s3D%proj_space_euls(  3,ref,self%ithr) = 360. - pftc_glob%get_rot(inpl_ind)
        s3D%proj_space_corrs(   ref,self%ithr) = corr
    end subroutine store_solution

    subroutine kill( self )
        class(strategy3D_srch), intent(inout) :: self
        call self%grad_shsrch_obj%kill
        call self%grad_shsrch_obj2%kill
        call self%grad_shsrch_first_obj%kill
        call self%opeaks%kill
        call self%o_prev%kill
        if( allocated(self%refine) ) deallocate(self%refine)
    end subroutine kill

end module simple_strategy3D_srch
