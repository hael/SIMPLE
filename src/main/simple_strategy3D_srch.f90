! common strategy3D methods and type specification for polymorphic strategy3D object creation are delegated to this class
module simple_strategy3D_srch
#include "simple_lib.f08"
use simple_params,            only: params
use simple_oris,              only: oris
use simple_ori,               only: ori
use simple_sym,               only: sym
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_pftcc_shsrch,      only: pftcc_shsrch       ! simplex-based angle and shift search
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based angle and shift search
use simple_strategy3D_alloc   ! use all in there
implicit none

public :: strategy3D_srch, strategy3D_spec
private

logical, parameter :: DEBUG = .false.

type strategy3D_spec
    class(params),           pointer :: pp            => null()
    class(polarft_corrcalc), pointer :: ppftcc        => null()
    class(oris),             pointer :: pa            => null()
    class(sym),              pointer :: pse           => null()
    integer,                 pointer :: nnmat(:,:)    => null()
    integer,                 pointer :: grid_projs(:) => null()
    integer,                 pointer :: symmat(:,:)   => null()
    integer :: iptcl=0, iptcl_map=0, szsn=0
    logical :: do_extr=.false.
    real    :: corr_thresh=0.
end type strategy3D_spec

type strategy3D_srch
    class(polarft_corrcalc), pointer :: pftcc_ptr => null()       !< corrcalc object
    class(oris),             pointer :: a_ptr     => null()       !< b%a (primary particle orientation table)
    class(sym),              pointer :: se_ptr    => null()       !< b%se (symmetry elements)
    type(pftcc_shsrch)               :: shsrch_obj                !< origin shift search object
    type(pftcc_shsrch_grad)          :: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
    integer, allocatable             :: nnvec(:)                  !< nearest neighbours indices
    integer                          :: iptcl         = 0         !< global particle index
    integer                          :: iptcl_map     = 0         !< index in pre-allocated 2D arrays
    integer                          :: nrefs         = 0         !< total # references (nstates*nprojs)
    integer                          :: nnnrefs       = 0         !< total # neighboring references (nstates*nnn)
    integer                          :: nstates       = 0         !< # states
    integer                          :: nprojs        = 0         !< # projections
    integer                          :: nrots         = 0         !< # in-plane rotations in polar representation
    integer                          :: npeaks        = 0         !< # peaks (nonzero orientation weights)
    integer                          :: npeaks_eff    = 0         !< effective # peaks
    integer                          :: npeaks_grid   = 0         !< # peaks after coarse search
    integer                          :: nbetter       = 0         !< # better orientations identified
    integer                          :: nrefs_eval    = 0         !< # references evaluated
    integer                          :: nnn_static    = 0         !< # nearest neighbors (static)
    integer                          :: nnn           = 0         !< # nearest neighbors (dynamic)
    integer                          :: prev_roind    = 0         !< previous in-plane rotation index
    integer                          :: prev_state    = 0         !< previous state index
    integer                          :: prev_ref      = 0         !< previous reference index
    integer                          :: kstop_grid    = 0         !< Frequency limit of first coarse grid search
    integer                          :: nsym          = 0         !< symmetry order
    real                             :: prev_corr     = 1.        !< previous best correlation
    real                             :: specscore     = 0.        !< spectral score
    real                             :: prev_shvec(2) = 0.        !< previous origin shift vector
    character(len=STDLEN)            :: opt           = ''        !< optimizer flag
    logical                          :: neigh         = .false.   !< nearest neighbour refinement flag
    logical                          :: doshift       = .true.    !< 2 indicate whether 2 serch shifts
    logical                          :: exists        = .false.   !< 2 indicate existence
  contains
    procedure :: new
    procedure :: prep4srch
    procedure :: greedy_subspace_srch
    procedure :: inpl_srch
    procedure :: calc_corr
    procedure :: store_solution
    procedure :: kill
end type strategy3D_srch

contains

    subroutine new( self, spec )
        class(strategy3D_srch), intent(inout) :: self
        class(strategy3D_spec), intent(in)    :: spec
        integer :: nstates_eff
        real    :: lims(2,2), lims_init(2,2)
        ! set constants
        self%pftcc_ptr  => spec%ppftcc
        self%a_ptr      => spec%pa
        self%se_ptr     => spec%pse
        self%iptcl      =  spec%iptcl
        self%iptcl_map  =  spec%iptcl_map
        self%nstates    =  spec%pp%nstates
        self%nprojs     =  spec%pp%nspace
        self%nrefs      =  self%nprojs*self%nstates
        self%nrots      =  round2even(twopi*real(spec%pp%ring2))
        self%npeaks     =  spec%pp%npeaks
        self%nbetter    =  0
        self%nrefs_eval =  0
        self%nsym       =  self%se_ptr%get_nsym()
        self%doshift    =  spec%pp%l_doshift
        self%neigh      =  spec%pp%neigh == 'yes'
        self%nnn_static =  spec%pp%nnn
        self%nnn        =  spec%pp%nnn
        self%nnnrefs    =  self%nnn*self%nstates
        self%kstop_grid =  spec%pp%kstop_grid
        self%opt        =  spec%pp%opt
        ! multiple states
        if( self%nstates == 1 )then
            self%npeaks_grid = GRIDNPEAKS
        else
            ! number of populated states
            nstates_eff = count(state_exists)
            select case(trim(spec%pp%refine))
                case('cluster')
                    self%npeaks_grid = 1
                case DEFAULT
                    ! "-(nstates_eff-1)" because all states share the same previous orientation
                    self%npeaks_grid = GRIDNPEAKS * nstates_eff - (nstates_eff - 1)
            end select
        endif
        if( spec%pp%neigh.eq.'yes' )then
            self%npeaks_grid = min(self%npeaks_grid,self%nnnrefs)
        else
            self%npeaks_grid = min(self%npeaks_grid,self%nrefs)
        endif
        ! create in-plane search objects
        lims(:,1)      = -spec%pp%trs
        lims(:,2)      =  spec%pp%trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        call self%shsrch_obj%new(self%pftcc_ptr, lims, lims_init=lims_init,&
            &shbarrier=spec%pp%shbarrier, nrestarts=3, maxits=60)
        call self%grad_shsrch_obj%new(self%pftcc_ptr, lims, lims_init=lims_init,&
            &shbarrier=spec%pp%shbarrier, maxits=60)
        self%exists = .true.
        if( DEBUG ) print *, '>>> STRATEGY3D_SRCH :: CONSTRUCTED NEW STRATEGY3D_SRCH OBJECT'
    end subroutine new

    subroutine prep4srch( self, nnmat, target_projs )
        use simple_combinatorics, only: merge_into_disjoint_set
        class(strategy3D_srch), intent(inout) :: self
        integer, optional,   intent(in)    :: nnmat(self%nprojs,self%nnn_static), target_projs(self%npeaks_grid)
        real      :: corrs(self%nrots)
        type(ori) :: o_prev
        real      :: corr, bfac
        if( self%neigh )then
            if( .not. present(nnmat) )&
            &stop 'need optional nnmat to be present for refine=neigh modes :: prep4srch (strategy3D_srch)'
            if( .not. present(target_projs) )&
            &stop 'need optional target_projs to be present for refine=neigh modes :: prep4srch (strategy3D_srch)'
        endif
        o_prev          = self%a_ptr%get_ori(self%iptcl)
        self%prev_state = o_prev%get_state()                                          ! state index
        self%prev_roind = self%pftcc_ptr%get_roind(360.-o_prev%e3get())               ! in-plane angle index
        self%prev_shvec = o_prev%get_2Dshift()                                        ! shift vector
        self%prev_ref   = (self%prev_state-1)*self%nprojs + prev_proj(self%iptcl_map) ! previous reference
        if( self%prev_state > 0 )then
            if( self%prev_state > self%nstates ) stop 'previous best state outside boundary; prep4srch; simple_strategy3D_srch'
            if( .not. state_exists(self%prev_state) ) stop 'empty previous state; prep4srch; simple_strategy3D_srch'
        endif
        if( self%neigh )then
            ! disjoint nearest neighbour set
            self%nnvec = merge_into_disjoint_set(self%nprojs, self%nnn_static, nnmat, target_projs)
        endif
        ! calc specscore
        self%specscore = self%pftcc_ptr%specscore(self%prev_ref, self%iptcl, self%prev_roind)
        ! B-factor memoization
        if( self%pftcc_ptr%objfun_is_ccres() )then
            if( .not. self%a_ptr%isthere(self%iptcl, 'bfac') )then
                bfac = self%pftcc_ptr%fit_bfac(self%prev_ref, self%iptcl, self%prev_roind, [0.,0.])
                call self%pftcc_ptr%memoize_bfac(self%iptcl, bfac)
            endif
        endif
        ! prep corr
        call self%pftcc_ptr%gencorrs(self%prev_ref, self%iptcl, corrs)
        corr = max(0.,maxval(corrs))
        if( corr - 1.0 > 1.0e-5 .or. .not. is_a_number(corr) )then
            print *, 'FLOATING POINT EXCEPTION ALARM; simple_strategy3D_srch :: prep4srch'
            print *, 'corr > 1. or isNaN'
            print *, 'corr = ', corr
            if( corr > 1. )               corr = 1.
            if( .not. is_a_number(corr) ) corr = 0.
            call o_prev%print_ori()
        endif
        self%prev_corr = corr
        if( DEBUG ) print *, '>>> STRATEGY3D_SRCH :: PREPARED FOR SIMPLE_STRATEGY3D_SRCH'
    end subroutine prep4srch

    subroutine greedy_subspace_srch( self, grid_projs, target_projs )
        class(strategy3D_srch), intent(inout) :: self
        integer,             intent(in)    :: grid_projs(:)
        integer,             intent(inout) :: target_projs(:)
        real      :: inpl_corrs(self%nrots), corrs(self%nrefs), inpl_corr
        integer   :: iref, isample, nrefs, ntargets, cnt, istate
        integer   :: state_cnt(self%nstates), iref_state, inpl_ind(1)
        if( self%a_ptr%get_state(self%iptcl) > 0 )then
            ! initialize
            target_projs = 0
            proj_space_corrs(self%iptcl_map,:) = -1.
            nrefs = size(grid_projs)
            ! search
            do isample = 1, nrefs
                do istate = 1, self%nstates
                    iref = grid_projs(isample)           ! set the projdir reference index
                    iref = (istate-1)*self%nprojs + iref ! set the state reference index
                    call per_ref_srch                    ! actual search
                end do
            end do
            ! sort in correlation projection direction space
            corrs = proj_space_corrs(self%iptcl_map,:)
            call hpsort(corrs, proj_space_inds(self%iptcl_map,:))
            ! return target points
            ntargets = size(target_projs)
            cnt      = 1
            target_projs( cnt ) = prev_proj(self%iptcl_map) ! previous always part of the targets
            if( self%nstates == 1 )then
                ! Single state
                do isample=self%nrefs,self%nrefs - ntargets + 1,-1
                    if( target_projs(1) == proj_space_inds(self%iptcl_map,isample) )then
                        ! direction is already in target set
                    else
                        cnt = cnt + 1
                        target_projs(cnt) = proj_space_inds(self%iptcl_map,isample)
                        if( cnt == ntargets ) exit
                    endif
                end do
            else
                ! Multiples states
                state_cnt = 1                                            ! previous always part of the targets
                do isample = self%nrefs, 1, -1
                    if( cnt >= self%npeaks_grid )exit                    ! all that we need
                    iref_state = proj_space_inds(self%iptcl_map,isample) ! reference index to multi-state space
                    istate     = ceiling(real(iref_state)/real(self%nprojs))
                    iref       = iref_state - (istate-1)*self%nprojs     ! reference index to single state space
                    if( .not.state_exists(istate) )cycle
                    if( any(target_projs == iref) )cycle                 ! direction is already set
                    if( state_cnt(istate) >= GRIDNPEAKS )cycle           ! state is already filled
                    cnt = cnt + 1
                    target_projs(cnt) = iref
                    state_cnt(istate) = state_cnt(istate) + 1
                end do
            endif
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        if( DEBUG ) print *, '>>> STRATEGY3D_SRCH :: FINISHED GREEDY SUBSPACE SEARCH'

        contains

            subroutine per_ref_srch
                if( state_exists(istate) )then
                    call self%pftcc_ptr%gencorrs(iref, self%iptcl, self%kstop_grid, inpl_corrs) ! In-plane correlations
                    inpl_ind   = maxloc(inpl_corrs)                   ! greedy in-plane
                    inpl_corr  = inpl_corrs(inpl_ind(1))              ! max in plane correlation
                    proj_space_corrs(self%iptcl_map,iref) = inpl_corr ! stash in-plane correlation for sorting
                    proj_space_inds(self%iptcl_map,iref)  = iref      ! stash the index for sorting
                endif
            end subroutine per_ref_srch

    end subroutine greedy_subspace_srch

    subroutine inpl_srch( self )
        class(strategy3D_srch), intent(inout) :: self
        real, allocatable :: cxy(:)
        integer :: i, ref, irot
        if( self%doshift )then
            select case(trim(self%opt))
                case('bfgs')
                    do i=self%nrefs,self%nrefs-self%npeaks+1,-1
                        ref = proj_space_inds(self%iptcl_map, i)
                        call self%grad_shsrch_obj%set_indices(ref, self%iptcl)
                        cxy = self%grad_shsrch_obj%minimize(irot=irot)
                        if( irot > 0 )then
                            ! irot > 0 guarantees improvement found
                            proj_space_euls(self%iptcl_map, ref,3) = 360. - self%pftcc_ptr%get_rot(irot)
                            proj_space_corrs(self%iptcl_map,ref)   = cxy(1)
                            proj_space_shift(self%iptcl_map,ref,:) = cxy(2:3)
                        endif
                    end do
                case DEFAULT
                    do i=self%nrefs,self%nrefs-self%npeaks+1,-1
                        ref = proj_space_inds(self%iptcl_map, i)
                        call self%shsrch_obj%set_indices(ref, self%iptcl)
                        cxy = self%shsrch_obj%minimize(irot=irot)
                        if( irot > 0 )then
                            ! irot > 0 guarantees improvement found
                            proj_space_euls(self%iptcl_map, ref,3) = 360. - self%pftcc_ptr%get_rot(irot)
                            proj_space_corrs(self%iptcl_map,ref)   = cxy(1)
                            proj_space_shift(self%iptcl_map,ref,:) = cxy(2:3)
                        endif
                    end do
            end select
        endif
        if( DEBUG ) print *, '>>> STRATEGY3D_SRCH :: FINISHED INPL SEARCH'
    end subroutine inpl_srch

    subroutine calc_corr( self )
        class(strategy3D_srch), intent(inout) :: self
        integer :: iref
        real    :: cc
        if( prev_states(self%iptcl_map) > 0 )then
            self%prev_roind = self%pftcc_ptr%get_roind(360.-self%a_ptr%e3get(self%iptcl))
            iref = (prev_states(self%iptcl_map)-1) * self%nprojs + prev_proj(self%iptcl_map)
            cc = self%pftcc_ptr%gencorr_cc_for_rot(iref, self%iptcl, [0.,0.], self%prev_roind)
            call self%a_ptr%set(self%iptcl,'corr', cc)
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        if( DEBUG ) print *, '>>> STRATEGY3D_SRCH :: FINISHED CALC_CORR'
    end subroutine calc_corr

    subroutine store_solution( self, ind, ref, inpl_ind, corr )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: ind, ref, inpl_ind
        real,                   intent(in)    :: corr
        proj_space_inds(self%iptcl_map,ind)   = ref
        proj_space_euls(self%iptcl_map,ref,3) = 360. - self%pftcc_ptr%get_rot(inpl_ind)
        proj_space_corrs(self%iptcl_map,ref)  = corr
    end subroutine store_solution

    subroutine kill( self )
        class(strategy3D_srch), intent(inout) :: self
        if(allocated(self%nnvec))deallocate(self%nnvec)
        call self%grad_shsrch_obj%kill
    end subroutine kill

end module simple_strategy3D_srch
