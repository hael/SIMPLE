! common strategy3D methods and type specification for polymorphic strategy3D object creation are delegated to this class
module simple_strategy3D_srch
include 'simple_lib.f08'
use simple_oris,              only: oris
use simple_ori,               only: ori
use simple_sym,               only: sym
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
use simple_strategy3D_alloc   ! singleton s3D
use simple_polarft_corrcalc,  only: pftcc_glob
use simple_parameters,        only: params_glob
use simple_builder,           only: build_glob
implicit none

public :: strategy3D_srch, strategy3D_spec
private

#include "simple_local_flags.inc"

type strategy3D_spec
    integer, pointer :: grid_projs(:) => null()
    integer, pointer :: symmat(:,:)   => null()
    integer :: iptcl=0, iptcl_map=0, szsn=0
    logical :: do_extr=.false.
    real    :: corr_thresh=0.
end type strategy3D_spec

type strategy3D_srch
    type(pftcc_shsrch_grad) :: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
    integer, allocatable    :: nnvec(:)                  !< nearest neighbours indices
    integer                 :: iptcl         = 0         !< global particle index
    integer                 :: iptcl_map     = 0         !< index in pre-allocated 2D arrays
    integer                 :: kstop_grid    = 0         !< Frequency limit of first coarse grid search
    integer                 :: nrefs         = 0         !< total # references (nstates*nprojs)
    integer                 :: nnnrefs       = 0         !< total # neighboring references (nstates*nnn)
    integer                 :: nstates       = 0         !< # states
    integer                 :: nprojs        = 0         !< # projections
    integer                 :: nrots         = 0         !< # in-plane rotations in polar representation
    integer                 :: npeaks        = 0         !< # peaks (nonzero orientation weights)
    integer                 :: npeaks_eff    = 0         !< effective # peaks
    integer                 :: npeaks_grid   = 0         !< # peaks after coarse search
    integer                 :: nsym          = 0         !< symmetry order
    integer                 :: nbetter       = 0         !< # better orientations identified
    integer                 :: nrefs_eval    = 0         !< # references evaluated
    integer                 :: nnn_static    = 0         !< # nearest neighbors (static)
    integer                 :: nnn           = 0         !< # nearest neighbors (dynamic)
    integer                 :: prev_roind    = 0         !< previous in-plane rotation index
    integer                 :: prev_state    = 0         !< previous state index
    integer                 :: prev_ref      = 0         !< previous reference index
    real                    :: prev_corr     = 1.        !< previous best correlation
    real                    :: specscore     = 0.        !< spectral score
    real                    :: prev_shvec(2) = 0.        !< previous origin shift vector
    logical                 :: neigh         = .false.   !< nearest neighbour refinement flag
    logical                 :: doshift       = .true.    !< 2 indicate whether 2 serch shifts
    logical                 :: exists        = .false.   !< 2 indicate existence
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

    subroutine new( self, spec, npeaks )
        class(strategy3D_srch), intent(inout) :: self
        class(strategy3D_spec), intent(in)    :: spec
        integer,                intent(in)    :: npeaks
        integer, parameter :: MAXITS = 60
        integer :: nstates_eff
        real    :: lims(2,2), lims_init(2,2)
        ! set constants
        self%iptcl      =  spec%iptcl
        self%iptcl_map  =  spec%iptcl_map
        self%nstates    =  params_glob%nstates
        self%nprojs     =  params_glob%nspace
        self%nrefs      =  self%nprojs*self%nstates
        self%nrots      =  round2even(twopi*real(params_glob%ring2))
        self%npeaks     =  npeaks
        self%nbetter    =  0
        self%nrefs_eval =  0
        self%nsym       =  build_glob%pgrpsyms%get_nsym()
        self%doshift    =  params_glob%l_doshift
        self%neigh      =  params_glob%neigh == 'yes'
        self%nnn_static =  params_glob%nnn
        self%nnn        =  params_glob%nnn
        self%nnnrefs    =  self%nnn*self%nstates
        self%kstop_grid =  params_glob%kstop_grid
        ! multiple states
        if( self%nstates == 1 )then
            self%npeaks_grid = GRIDNPEAKS
        else
            ! number of populated states
            nstates_eff = count(s3D%state_exists)
            select case(trim(params_glob%refine))
            case('cluster','clustersym','clusterdev')
                    self%npeaks_grid = 1
                case DEFAULT
                    ! "-(nstates_eff-1)" because all states share the same previous orientation
                    self%npeaks_grid = GRIDNPEAKS * nstates_eff - (nstates_eff - 1)
            end select
        endif
        if( self%neigh )then
            self%npeaks_grid = min(self%npeaks_grid,self%nnnrefs)
        else
            self%npeaks_grid = min(self%npeaks_grid,self%nrefs)
        endif
        ! create in-plane search objects
        lims(:,1)      = -params_glob%trs
        lims(:,2)      =  params_glob%trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        call self%grad_shsrch_obj%new(lims, lims_init=lims_init,&
            &shbarrier=params_glob%shbarrier, maxits=MAXITS)
        self%exists = .true.
        DebugPrint  '>>> STRATEGY3D_SRCH :: CONSTRUCTED NEW STRATEGY3D_SRCH OBJECT'
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
        o_prev          = build_glob%spproj_field%get_ori(self%iptcl)
        self%prev_state = o_prev%get_state()                                          ! state index
        self%prev_roind = pftcc_glob%get_roind(360.-o_prev%e3get())               ! in-plane angle index
        self%prev_shvec = o_prev%get_2Dshift()                                        ! shift vector
        self%prev_ref   = (self%prev_state-1)*self%nprojs + s3D%prev_proj(self%iptcl_map) ! previous reference
        if( self%prev_state > 0 )then
            if( self%prev_state > self%nstates ) stop 'previous best state outside boundary; prep4srch; simple_strategy3D_srch'
            if( .not. s3D%state_exists(self%prev_state) ) stop 'empty previous state; prep4srch; simple_strategy3D_srch'
        endif
        if( self%neigh )then
            ! disjoint nearest neighbour set
            self%nnvec = merge_into_disjoint_set(self%nprojs, self%nnn_static, nnmat, target_projs)
        endif
        ! B-factor memoization
        if( pftcc_glob%objfun_is_ccres() )then
            bfac = pftcc_glob%fit_bfac(self%prev_ref, self%iptcl, self%prev_roind, [0.,0.])
            call pftcc_glob%memoize_bfac(self%iptcl, bfac)
            call build_glob%spproj_field%set(self%iptcl, 'bfac', bfac)
        endif
        ! calc specscore
        self%specscore = pftcc_glob%specscore(self%prev_ref, self%iptcl, self%prev_roind)
        ! prep corr
        call pftcc_glob%gencorrs(self%prev_ref, self%iptcl, corrs)
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
        DebugPrint  '>>> STRATEGY3D_SRCH :: PREPARED FOR SIMPLE_STRATEGY3D_SRCH'
    end subroutine prep4srch

    subroutine greedy_subspace_srch( self, grid_projs, target_projs )
        class(strategy3D_srch), intent(inout) :: self
        integer,             intent(in)    :: grid_projs(:)
        integer,             intent(inout) :: target_projs(:)
        real      :: inpl_corrs(self%nrots), corrs(self%nrefs), inpl_corr
        integer   :: iref, isample, nrefs, ntargets, cnt, istate
        integer   :: state_cnt(self%nstates), iref_state, inpl_ind(1)
        if( build_glob%spproj_field%get_state(self%iptcl) > 0 )then
            ! initialize
            target_projs = 0
            s3D%proj_space_corrs(self%iptcl_map,:) = -1.
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
            corrs = s3D%proj_space_corrs(self%iptcl_map,:)
            call hpsort(corrs, s3D%proj_space_inds(self%iptcl_map,:))
            ! return target points
            ntargets = size(target_projs)
            cnt      = 1
            target_projs( cnt ) = s3D%prev_proj(self%iptcl_map) ! previous always part of the targets
            if( self%nstates == 1 )then
                ! Single state
                do isample=self%nrefs,self%nrefs - ntargets + 1,-1
                    if( target_projs(1) == s3D%proj_space_inds(self%iptcl_map,isample) )then
                        ! direction is already in target set
                    else
                        cnt = cnt + 1
                        target_projs(cnt) = s3D%proj_space_inds(self%iptcl_map,isample)
                        if( cnt == ntargets ) exit
                    endif
                end do
            else
                ! Multiples states
                state_cnt = 1                                            ! previous always part of the targets
                do isample = self%nrefs, 1, -1
                    if( cnt >= self%npeaks_grid )exit                    ! all that we need
                    iref_state = s3D%proj_space_inds(self%iptcl_map,isample) ! reference index to multi-state space
                    istate     = ceiling(real(iref_state)/real(self%nprojs))
                    iref       = iref_state - (istate-1)*self%nprojs     ! reference index to single state space
                    if( .not.s3D%state_exists(istate) )cycle
                    if( any(target_projs == iref) )cycle                 ! direction is already set
                    if( state_cnt(istate) >= GRIDNPEAKS )cycle           ! state is already filled
                    cnt = cnt + 1
                    target_projs(cnt) = iref
                    state_cnt(istate) = state_cnt(istate) + 1
                end do
            endif
        else
            call build_glob%spproj_field%reject(self%iptcl)
        endif
        DebugPrint  '>>> STRATEGY3D_SRCH :: FINISHED GREEDY SUBSPACE SEARCH'

        contains

            subroutine per_ref_srch
                if( s3D%state_exists(istate) )then
                    call pftcc_glob%gencorrs(iref, self%iptcl, self%kstop_grid, inpl_corrs) ! In-plane correlations
                    inpl_ind   = maxloc(inpl_corrs)                   ! greedy in-plane
                    inpl_corr  = inpl_corrs(inpl_ind(1))              ! max in plane correlation
                    s3D%proj_space_corrs(self%iptcl_map,iref) = inpl_corr ! stash in-plane correlation for sorting
                    s3D%proj_space_inds(self%iptcl_map,iref)  = iref      ! stash the index for sorting
                endif
            end subroutine per_ref_srch

    end subroutine greedy_subspace_srch

    subroutine inpl_srch( self )
        class(strategy3D_srch), intent(inout) :: self
        real    :: cxy(3)
        integer :: i, ref, irot
        if( self%doshift )then
            ! BFGS
            do i=self%nrefs,self%nrefs-self%npeaks+1,-1
                ref = s3D%proj_space_inds(self%iptcl_map, i)
                call self%grad_shsrch_obj%set_indices(ref, self%iptcl)
                cxy = self%grad_shsrch_obj%minimize(irot=irot)
                if( irot > 0 )then
                    ! irot > 0 guarantees improvement found
                    s3D%proj_space_euls(self%iptcl_map, ref,3) = 360. - pftcc_glob%get_rot(irot)
                    s3D%proj_space_corrs(self%iptcl_map,ref)   = cxy(1)
                    s3D%proj_space_shift(self%iptcl_map,ref,:) = cxy(2:3)
                endif
            end do
        endif
        DebugPrint  '>>> STRATEGY3D_SRCH :: FINISHED INPL SEARCH'
    end subroutine inpl_srch

    subroutine calc_corr( self )
        class(strategy3D_srch), intent(inout) :: self
        integer :: iref
        real    :: cc
        self%prev_state = build_glob%spproj_field%get_state(self%iptcl)
        if( self%prev_state > 0 )then
            self%prev_roind = pftcc_glob%get_roind(360.-build_glob%spproj_field%e3get(self%iptcl))
            iref = (self%prev_state-1) * self%nprojs + s3D%prev_proj(self%iptcl_map)
            cc = pftcc_glob%gencorr_cc_for_rot(iref, self%iptcl, [0.,0.], self%prev_roind)
            call build_glob%spproj_field%set(self%iptcl,'corr', cc)
        else
            call build_glob%spproj_field%reject(self%iptcl)
        endif
        DebugPrint  '>>> STRATEGY3D_SRCH :: FINISHED CALC_CORR'
    end subroutine calc_corr

    subroutine store_solution( self, ind, ref, inpl_ind, corr )
        class(strategy3D_srch), intent(inout) :: self
        integer,                intent(in)    :: ind, ref, inpl_ind
        real,                   intent(in)    :: corr
        s3D%proj_space_inds(self%iptcl_map,ind)   = ref
        s3D%proj_space_euls(self%iptcl_map,ref,3) = 360. - pftcc_glob%get_rot(inpl_ind)
        s3D%proj_space_corrs(self%iptcl_map,ref)  = corr
    end subroutine store_solution

    subroutine kill( self )
        class(strategy3D_srch), intent(inout) :: self
        if(allocated(self%nnvec))deallocate(self%nnvec)
        call self%grad_shsrch_obj%kill
    end subroutine kill

end module simple_strategy3D_srch
