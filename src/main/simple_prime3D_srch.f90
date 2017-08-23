! PRIME3D stochastic search routines
module simple_prime3D_srch
use simple_defs
use simple_syslib
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_strings,          only: str_has_substr, int2str_pad
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_pftcc_shsrch,     only: pftcc_shsrch
use simple_pftcc_inplsrch,   only: pftcc_inplsrch
use simple_math              ! use all in there
implicit none

public :: prime3D_srch
private
#include "simple_local_flags.inc"

real,    parameter :: FACTWEIGHTS_THRESH = 0.001        !< threshold for factorial weights
integer, parameter :: INPL_EXPAND_FAC    = 3            !< expansion factor for in-plane weighting

!> struct for prime3d params
type prime3D_srch
    private
    type(oris)              :: o_refs                   !< projection directions search space
    type(oris)              :: o_peaks                  !< orientations of best npeaks oris
    type(pftcc_shsrch)      :: shsrch_obj               !< origin shift search object
    type(pftcc_inplsrch)    :: inplsrch_obj             !< in-plane search object
    integer                 :: nrefs          = 0       !< total # references (nstates*nprojs)
    integer                 :: nnnrefs        = 0       !< total # neighboring references (nstates*nnn)
    integer                 :: nstates        = 0       !< # states
    integer                 :: nprojs         = 0       !< # projections (same as # oris in o_refs)
    integer                 :: nrots          = 0       !< # in-plane rotations in polar representation
    integer                 :: npeaks         = 0       !< # peaks (nonzero orientation weights)
    integer                 :: npeaks_grid    = 0       !< # peaks after coarse search
    integer                 :: nbetter        = 0       !< # better orientations identified
    integer                 :: nrefs_eval     = 0       !< # references evaluated
    integer                 :: nnn_static     = 0       !< # nearest neighbors (static)
    integer                 :: nnn            = 0       !< # nearest neighbors (dynamic)
    integer                 :: nsym           = 0       !< symmetry order
    integer                 :: prev_roind     = 0       !< previous in-plane rotation index
    integer                 :: prev_state     = 0       !< previous state index
    integer                 :: prev_ref       = 0       !< previous reference index
    integer                 :: prev_proj      = 0       !< previous projection index
    integer                 :: kstop_grid     = 0       !< Frequency limit of first coarse grid search
    real                    :: prev_corr      = 1.      !< previous best correlation
    real                    :: specscore      = 0.      !< spectral score
    real                    :: prev_shvec(2)  = 0.      !< previous origin shift vector
    real                    :: lims(2,2)      = 0.      !< shift search range limit
    real                    :: athres         = 0.      !< angular threshold
    real                    :: dfx            = 0.      !< ctf x-defocus
    real                    :: dfy            = 0.      !< ctf y-defocus
    real                    :: angast         = 0.      !< ctf astigmatism
    integer, allocatable    :: proj_space_inds(:)       !< projection space index array
    integer, allocatable    :: srch_order(:)            !< stochastic search order
    logical, allocatable    :: state_exists(:)          !< indicates whether each state is populated
    character(len=STDLEN)   :: refine         = ''      !< refinement flag
    character(len=STDLEN)   :: ctf            = ''      !< ctf flag
    character(len=STDLEN)   :: shbarr         = ''      !< shift barrier flag
    character(len=STDLEN)   :: pgrp           = 'c1'    !< point-group symmetry
    character(len=STDLEN)   :: eo             = 'no'    !< wheteher to weigh distances
    logical                 :: doshift        = .true.  !< 2 indicate whether 2 serch shifts
    logical                 :: greedy_inpl    = .true.  !< 2 indicate whether in-plane search is greedy or not
    logical                 :: exists         = .false. !< 2 indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SEARCH ROUTINES
    procedure          :: exec_prime3D_srch
    procedure          :: exec_prime3D_srch_het
    procedure, private :: dev_srch
    procedure, private :: greedy_srch
    procedure, private :: greedy_subspace_srch
    procedure, private :: stochastic_srch
    procedure, private :: stochastic_srch_shc
    procedure, private :: stochastic_srch_snhc
    procedure, private :: stochastic_srch_het
    procedure          :: inpl_srch
    ! PREPARATION ROUTINES
    procedure          :: prep4srch
    procedure          :: prep_reforis
    ! CALCULATORS
    procedure          :: prep_npeaks_oris
    procedure          :: stochastic_weights
    procedure, private :: calc_specscore

    ! GETTERS & SETTERS
    procedure          :: update_best
    procedure          :: get_ori
    procedure          :: get_oris
    procedure          :: store_solution
    procedure          :: set_o_peaks
    procedure          :: get_o_peaks
    procedure          :: get_o_refs
    procedure          :: get_srch_order
    procedure          :: get_ntotrefs
    procedure          :: get_nrefs
    procedure          :: get_nrots
    procedure          :: get_prevstate
    procedure          :: get_prevcorr
    procedure          :: get_prevroind
    procedure          :: get_prevshvec
    procedure          :: get_prevref
    procedure          :: get_angast
    procedure          :: get_dfx
    procedure          :: get_dfy
    procedure          :: set_ctf
    procedure          :: set_greedy_inpl
    ! MEMORY MANAGEMENT
    procedure          :: online_allocate
    procedure          :: online_destruct
    ! DESTRUCTOR
    procedure          :: kill
end type prime3D_srch

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine new( self, a, p, pftcc )
        use simple_params, only: params
        class(prime3D_srch),     intent(inout) :: self     !< instance
        class(oris),             intent(inout) :: a        !< oris
        class(params),           intent(in)    :: p        !< parameters
        class(polarft_corrcalc), intent(inout) :: pftcc    !< correlator
        integer  :: alloc_stat, nstates_eff
        ! destroy possibly pre-existing instance
        call self%kill
        ! set constants
        self%nstates     =  p%nstates
        self%nprojs      =  p%nspace
        self%nrefs       =  self%nprojs*self%nstates
        self%nrots       =  round2even(twopi*real(p%ring2))
        self%npeaks      =  p%npeaks
        self%nbetter     =  0
        self%nrefs_eval  =  0
        self%athres      =  p%athres
        self%doshift     =  p%doshift
        self%refine      =  p%refine
        self%ctf         =  p%ctf
        self%nnn_static  =  p%nnn
        self%nnn         =  p%nnn
        self%nnnrefs     =  self%nnn*self%nstates
        self%shbarr      =  p%shbarrier
        self%pgrp        =  p%pgrp
        self%kstop_grid  =  p%kstop_grid
        self%greedy_inpl = .true.
        self%eo          =  p%eo
        if( str_has_substr(self%refine,'shc') )then
            if( self%npeaks > 1 ) stop 'npeaks must be equal to 1 with refine=shc|shcneigh'
            self%greedy_inpl = .false.
        endif
        ! construct composites
        if( p%oritab.ne.'' )then
            self%state_exists = a%states_exist(self%nstates)
        else
            allocate(self%state_exists(self%nstates), stat=alloc_stat)   
            call alloc_errchk('In: new; simple_prime3D_srch, 1', alloc_stat)
            self%state_exists = .true.
        endif
        ! multiple states
        if( self%nstates == 1 )then
            self%npeaks_grid = GRIDNPEAKS
        else
            ! number of populated states
            nstates_eff = count(self%state_exists)
            select case(trim(p%refine))
                case('het')
                    self%npeaks_grid = 1
                case('exp')
                    self%npeaks_grid = GRIDNPEAKS
                case DEFAULT
                    ! "-(nstates_eff-1)" because all states share the same previous orientation
                    self%npeaks_grid = GRIDNPEAKS*nstates_eff - (nstates_eff-1)
            end select
        endif
        if(str_has_substr(self%refine,'neigh'))then
            self%npeaks_grid = min(self%npeaks_grid,self%nnnrefs)
        else
            self%npeaks_grid = min(self%npeaks_grid,self%nrefs)
        endif
        ! generate oris oject in which the best npeaks refs will be stored
        call self%o_peaks%new(self%npeaks)
        ! updates option to search shift
        self%doshift = p%doshift
        if( self%doshift )then
            self%lims(:,1) = -p%trs
            self%lims(:,2) =  p%trs
        endif
        self%exists = .true.
        ! create in-plane search objects
        call self%shsrch_obj%new(pftcc, self%lims, shbarrier=self%shbarr)
        call self%inplsrch_obj%new(pftcc, self%lims, shbarrier=self%shbarr)
        DebugPrint '>>> PRIME3D_SRCH::CONSTRUCTED NEW SIMPLE_PRIME3D_SRCH OBJECT'
    end subroutine new

    ! SEARCH ROUTINES
    
    !>  \brief  exec_prime3D_srch is a master prime search routine
    !! \param pftcc polarft corrcalc search storage
    !! \param iptcl particle index
    !! \param a,e search orientation
    !! \param lp low-pass cutoff freq
    !! \param greedy are we processing using greedy
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs grid projections
    !! \param szsn optional size for snhc refinement
    !!
    subroutine exec_prime3D_srch( self, pftcc, iptcl, a, e, lp, greedy, nnmat, grid_projs, szsn )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(oris),             intent(inout) :: a, e
        real,                    intent(in)    :: lp
        logical, optional,       intent(in)    :: greedy
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:), szsn
        logical :: ggreedy
        ggreedy = .false.
        if( present(greedy) ) ggreedy = greedy
        call self%online_allocate
        if( trim(self%refine).eq.'exp' )then
            ! experimental refinement mode
            call self%dev_srch(pftcc,  iptcl, a, e, lp, nnmat)
            !call self%o_peaks%write('ptcl_'//int2str_pad(iptcl,3)//'.txt')
        else
            ! classic refinement modes
            if( ggreedy )then
                call self%greedy_srch(pftcc,  iptcl, a, e, lp, nnmat, grid_projs)
            else if( self%refine.eq.'snhc' )then
                if( .not. present(szsn) )then
                    stop 'refine=snhc mode needs optional input szsn; simple_prime3D_srch :: exec_prime3D_srch'
                endif
                call self%stochastic_srch_snhc( pftcc, iptcl, a, e, lp, szsn )
            else if( str_has_substr(self%refine,'shc') )then
                call self%stochastic_srch_shc( pftcc, iptcl, a, e, lp, nnmat, grid_projs )
            else
                call self%stochastic_srch( pftcc, iptcl, a, e, lp, nnmat, grid_projs )
            endif
        endif
        call self%online_destruct
        DebugPrint '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SRCH'
    end subroutine exec_prime3D_srch

    !>  \brief state labeler
    subroutine exec_prime3D_srch_het( self, pftcc, iptcl, a, e, extr_bound, statecnt )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(oris),             intent(inout) :: a, e
        real,                    intent(in)    :: extr_bound
        integer,                 intent(inout) :: statecnt(self%nstates)
        call self%stochastic_srch_het(pftcc, iptcl, a, e, extr_bound, statecnt)
        call self%online_destruct
        DebugPrint '>>> PRIME3D_SRCH::EXECUTED PRIME3D_HET_SRCH'
    end subroutine exec_prime3D_srch_het

    !>  \brief  Dev search
    !! \param pftcc polarft corrcalc search storage
    !! \param iptcl particle index
    !! \param a,e search orientation
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs grid projections
    subroutine dev_srch( self, pftcc, iptcl, a, e, lp, nnmat )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(oris),             intent(inout) :: a, e
        real,                    intent(in)    :: lp
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn_static)
        real      :: projspace_corrs(self%nrefs), wcorr
        integer   :: iref, isample
        ! execute search
        if( nint(a%get(iptcl,'state')) > 0 )then
            ! initialize
            call self%prep4srch(pftcc, iptcl, a, e, lp, nnmat=nnmat)
            self%nbetter         = 0
            self%nrefs_eval      = 0
            self%proj_space_inds = 0
            projspace_corrs      = -1.
            ! search
            do isample = 1, self%nnn
                iref = self%srch_order(isample)        ! set the stochastic reference index
                if( iref == self%prev_ref ) cycle      ! previous best considered last
                call per_ref_srch(iref)                ! actual search
                if( self%nbetter >= self%npeaks ) exit ! exit condition
            end do
            if( self%nbetter < self%npeaks )then
                call per_ref_srch(self%prev_ref )      ! evaluate previous best ref last
            endif
            ! sort in correlation projection direction space
            call hpsort(self%nrefs, projspace_corrs, self%proj_space_inds) 
            call self%inpl_srch(pftcc, iptcl) ! search shifts
            ! prepare weights and orientations
            call self%prep_npeaks_oris
            call self%stochastic_weights(wcorr)
            call self%update_best(pftcc, iptcl, a)
            call a%set(iptcl, 'corr', wcorr)
        else
            call a%reject(iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED DEV SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind
                corrs     = pftcc%gencorrs(iref, iptcl) ! In-plane correlations
                loc       = maxloc(corrs)               ! greedy in-plane
                inpl_ind  = loc(1)                      ! in-plane angle index
                inpl_corr = corrs(inpl_ind)             ! max in plane correlation
                call self%store_solution(pftcc, iref, iref, inpl_ind, inpl_corr)
                projspace_corrs( iref ) = inpl_corr ! stash in-plane correlation for sorting                   
                ! update nbetter to keep track of how many improving solutions we have identified
                if( self%npeaks == 1 )then
                    if( inpl_corr > self%prev_corr ) self%nbetter = self%nbetter + 1
                else
                    if( inpl_corr >= self%prev_corr ) self%nbetter = self%nbetter + 1
                endif
                ! keep track of how many references we are evaluating
                self%nrefs_eval = self%nrefs_eval + 1
            end subroutine per_ref_srch

    end subroutine dev_srch

    !>  \brief  greedy hill-climbing
    !! \param pftcc polarft corrcalc search storage
    !! \param iptcl particle index
    !! \param a,e search orientation
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs grid projections
    subroutine greedy_srch( self, pftcc, iptcl, a, e, lp, nnmat, grid_projs )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(oris),             intent(inout) :: a, e
        real,                    intent(in)    :: lp
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:)
        real    :: projspace_corrs(self%nrefs),wcorr
        integer :: iref,isample,nrefs,target_projs(self%npeaks_grid)
        if( nint(a%get(iptcl,'state')) > 0 )then
            if( str_has_substr(self%refine, 'neigh') )then
                ! for neighbour modes we do a coarse grid search first
                if( .not. present(grid_projs) ) stop 'need optional grid_projs 4 subspace srch; prime3D_srch :: greedy_srch'
                call self%greedy_subspace_srch(pftcc, iptcl, a, e, grid_projs, target_projs)
                ! initialize
                call self%prep4srch(pftcc, iptcl, a, e, lp, nnmat, target_projs)
                nrefs = self%nnnrefs
            else
                ! initialize
                call self%prep4srch(pftcc, iptcl, a, e, lp)
                nrefs = self%nrefs
            endif
            self%nbetter         = 0
            self%nrefs_eval      = 0
            self%proj_space_inds = 0
            projspace_corrs      = -1.
            ! search
            do isample=1,nrefs
                iref = self%srch_order(isample) ! set the reference index
                call per_ref_srch(iref)         ! actual search
            end do
            ! in greedy mode, we evaluate all refs
            self%nrefs_eval = nrefs
            ! sort in correlation projection direction space
            call hpsort(self%nrefs, projspace_corrs, self%proj_space_inds)
            ! take care of the in-planes
            call self%inpl_srch(pftcc, iptcl)
            ! prepare weights & orientation
            call self%prep_npeaks_oris
            call self%stochastic_weights(wcorr)
            call self%update_best(pftcc, iptcl, a)
            call a%set(iptcl, 'corr', wcorr)
        else
            call a%reject(iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED GREEDY SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind, state
                state = 1
                if( self%nstates > 1 ) state = nint( self%o_refs%get(iref, 'state') )
                if( self%state_exists(state) )then
                    corrs     = pftcc%gencorrs(iref, iptcl)        ! In-plane correlations
                    loc       = maxloc(corrs)                      ! greedy in-plane
                    inpl_ind  = loc(1)                             ! in-plane angle index
                    inpl_corr = corrs(inpl_ind)                    ! max in plane correlation
                    call self%store_solution(pftcc, iref, iref, inpl_ind, inpl_corr)
                    projspace_corrs( iref ) = inpl_corr            ! stash in-plane correlation for sorting
                endif
            end subroutine per_ref_srch

    end subroutine greedy_srch

    !>  \brief  greedy hill-climbing (4 initialisation)
    !! \param pftcc polarft corrcalc search storage
    !! \param iptcl particle index
    !! \param a,e search orientation
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs grid projections
    !! \param target_projs target projections
    subroutine greedy_subspace_srch( self, pftcc, iptcl, a, e, grid_projs, target_projs )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(oris),             intent(inout) :: a, e
        integer,                 intent(in)    :: grid_projs(:)
        integer,                 intent(inout) :: target_projs(:)
        real      :: projspace_corrs(self%nrefs)
        integer   :: iref, isample, nrefs, ntargets, cnt, prev_proj, istate
        integer   :: state_cnt(self%nstates), iref_state
        type(ori) :: o_prev
        if( nint(a%get(iptcl,'state')) > 0 )then
            o_prev    = a%get_ori(iptcl)
            prev_proj = e%find_closest_proj(o_prev,1)
            ! initialize
            self%proj_space_inds = 0
            target_projs         = 0
            projspace_corrs      = -1.
            nrefs = size(grid_projs)
            ! search
            do isample = 1, nrefs
                do istate = 1, self%nstates
                    iref = grid_projs(isample)           ! set the projdir reference index
                    iref = (istate-1)*self%nprojs + iref ! set the state reference index
                    call per_ref_srch(iref)              ! actual search
                end do
            end do
            ! sort in correlation projection direction space
            call hpsort(self%nrefs, projspace_corrs, self%proj_space_inds)
            ! return target points
            ntargets = size(target_projs)
            cnt = 1
            target_projs( cnt ) = prev_proj ! previous always part of the targets
            if( self%nstates == 1 )then
                ! Single state
                do isample=self%nrefs,self%nrefs - ntargets + 1,-1
                    if( target_projs(1) == self%proj_space_inds(isample) )then
                        ! direction is already in target set
                    else
                        cnt = cnt + 1
                        target_projs(cnt) = self%proj_space_inds(isample)
                        if( cnt == ntargets ) exit
                    endif
                end do
            else
                ! Multiples states
                state_cnt = 1                                           ! previous always part of the targets
                do isample = self%nrefs, 1, -1
                    if( cnt >= self%npeaks_grid )exit                   ! all that we need
                    iref_state = self%proj_space_inds(isample)          ! reference index to multi-state space
                    istate     = ceiling(real(iref_state)/real(self%nprojs))
                    iref       = iref_state - (istate-1)*self%nprojs    ! reference index to single state space
                    if( .not.self%state_exists(istate) )cycle
                    if( any(target_projs == iref) )cycle                ! direction is already set
                    if( state_cnt(istate) >= GRIDNPEAKS )cycle          ! state is already filled
                    cnt = cnt + 1
                    target_projs(cnt) = iref
                    state_cnt(istate) = state_cnt(istate) + 1
                end do
            endif
        else
            call a%reject(iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED GREEDY SUBSPACE SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind
                if( self%state_exists(istate) )then
                    corrs     = pftcc%gencorrs(iref, iptcl, self%kstop_grid) ! In-plane correlations
                    loc       = maxloc(corrs)                                ! greedy in-plane
                    inpl_ind  = loc(1)                                       ! in-plane angle index
                    inpl_corr = corrs(inpl_ind)                              ! max in plane correlation
                    projspace_corrs( iref )      = inpl_corr                 ! stash in-plane correlation for sorting
                    self%proj_space_inds( iref ) = iref                      ! stash the index for sorting
                endif
            end subroutine per_ref_srch

    end subroutine greedy_subspace_srch

    !>  \brief  executes the stochastic soft orientation search
    !! \param pftcc polarft corrcalc search storage
    !! \param iptcl particle index
    !! \param a,e search orientation
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs grid projections
    subroutine stochastic_srch( self, pftcc, iptcl, a, e, lp, nnmat, grid_projs )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(oris),             intent(inout) :: a, e
        real,                    intent(in)    :: lp
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:)
        real      :: projspace_corrs(self%nrefs),wcorr
        integer   :: iref,isample,nrefs,target_projs(self%npeaks_grid)
        ! execute search
        if( nint(a%get(iptcl,'state')) > 0 )then
            if( str_has_substr(self%refine, 'neigh') )then
                ! for neighbour modes we do a coarse grid search first
                if( .not. present(grid_projs) ) stop 'need optional grid_projs 4 subspace srch; prime3D_srch :: stochastic_srch'
                call self%greedy_subspace_srch(pftcc, iptcl, a, e, grid_projs, target_projs)
                ! initialize
                call self%prep4srch(pftcc, iptcl, a, e, lp, nnmat, target_projs)
                nrefs = self%nnnrefs
            else
                ! initialize
                call self%prep4srch(pftcc, iptcl, a, e, lp)
                nrefs = self%nrefs
            endif
            ! initialize, ctd
            self%nbetter         = 0
            self%nrefs_eval      = 0
            self%proj_space_inds = 0
            projspace_corrs      = -1. 
            do isample=1,nrefs
                iref = self%srch_order(isample)        ! set the stochastic reference index
                if( iref == self%prev_ref ) cycle      ! previous best considered last
                call per_ref_srch(iref)                ! actual search
                if( self%nbetter >= self%npeaks ) exit ! exit condition
            end do
            if( self%nbetter < self%npeaks )then
                call per_ref_srch(self%prev_ref )      ! evaluate previous best ref last
            endif
            ! sort in correlation projection direction space
            call hpsort(self%nrefs, projspace_corrs, self%proj_space_inds) 
            call self%inpl_srch(pftcc, iptcl) ! search shifts
            ! prepare weights and orientations
            call self%prep_npeaks_oris
            call self%stochastic_weights(wcorr)
            call self%update_best(pftcc, iptcl, a)
            call a%set(iptcl, 'corr', wcorr)
        else
            call a%reject(iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind, state
                state = 1
                if( self%nstates > 1 ) state = nint( self%o_refs%get(iref, 'state') )
                if( self%state_exists(state) )then
                    corrs     = pftcc%gencorrs(iref, iptcl) ! In-plane correlations
                    loc       = maxloc(corrs)   ! greedy in-plane
                    inpl_ind  = loc(1)          ! in-plane angle index
                    inpl_corr = corrs(inpl_ind) ! max in plane correlation
                    call self%store_solution(pftcc, iref, iref, inpl_ind, inpl_corr)
                    projspace_corrs( iref ) = inpl_corr ! stash in-plane correlation for sorting                   
                    ! update nbetter to keep track of how many improving solutions we have identified
                    if( self%npeaks == 1 )then
                        if( inpl_corr > self%prev_corr ) self%nbetter = self%nbetter + 1
                    else
                        if( inpl_corr >= self%prev_corr ) self%nbetter = self%nbetter + 1
                    endif
                    ! keep track of how many references we are evaluating
                    self%nrefs_eval = self%nrefs_eval + 1
                endif
            end subroutine per_ref_srch

    end subroutine stochastic_srch

    !>  \brief  executes the stochastic hard orientation search using pure stochastic hill climbing
    !!          (no probabilistic weighting + stochastic search of in-plane angles)
    !! \param pftcc polarft corrcalc search storage
    !! \param iptcl particle index
    !! \param a,e search orientation
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs grid projections
    subroutine stochastic_srch_shc( self, pftcc, iptcl, a, e, lp, nnmat, grid_projs )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(oris),             intent(inout) :: a, e
        real,                    intent(in)    :: lp
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:)
        real    :: inpl_corr,corrs(self%nrots)
        integer :: iref,isample,nrefs,inpl_ind,loc(1),target_projs(self%npeaks_grid)
        logical :: found_better
        if( nint(a%get(iptcl,'state')) > 0 )then
            if( str_has_substr(self%refine, 'neigh') )then
                ! for neighbour modes we do a coarse grid search first
                if( .not. present(grid_projs) ) stop 'need optional grid_projs 4 subspace srch; prime3D_srch :: stochastic_srch_shc'
                call self%greedy_subspace_srch(pftcc, iptcl, a, e, grid_projs, target_projs)
                ! initialize
                call self%prep4srch(pftcc, iptcl, a, e, lp, nnmat, target_projs)
                nrefs = self%nnnrefs
            else
                ! initialize
                call self%prep4srch(pftcc, iptcl, a, e, lp)
                nrefs = self%nrefs
            endif
            ! initialize
            found_better         = .false.
            self%nrefs_eval      = 0
            self%proj_space_inds = 0
            ! search
            do isample=1,nrefs
                iref = self%srch_order( isample ) ! stochastic reference index
                if( iref == self%prev_ref ) cycle ! previous best considered last
                call per_ref_srch( iref )         ! actual search
                if( inpl_ind > 0 )then
                    ! correlation-improving reference found, store it
                    call self%store_solution(pftcc, self%nrefs, iref, inpl_ind, inpl_corr)
                    found_better = .true.         ! flag for found solution
                    exit
                endif
            end do
            if( .not. found_better )then
                iref = self%prev_ref                          
                call per_ref_srch( iref )
                if( inpl_ind == 0 )then
                    ! be greedy
                    loc       = maxloc(corrs)
                    inpl_ind  = loc(1)
                    inpl_corr = corrs(inpl_ind)
                endif
                call self%store_solution(pftcc, self%nrefs, iref, inpl_ind, inpl_corr)
            endif
            ! search shifts
            call self%inpl_srch(pftcc, iptcl)
            ! output
            call self%prep_npeaks_oris
            call self%update_best(pftcc, iptcl, a)
        else
            call a%reject(iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH (REFINE=SHC|SHCNEIGH)'

        contains
        
            subroutine per_ref_srch( iref )
                use simple_rnd, only: shcloc
                integer, intent(in) :: iref
                integer :: state
                inpl_ind  = 0                                            ! init in-plane index 
                inpl_corr = 0.                                           ! init correlation
                state     = 1
                if( self%nstates > 1 ) state = nint( self%o_refs%get(iref, 'state') )
                if( self%state_exists( state ) )then
                    corrs    = pftcc%gencorrs( iref, iptcl )             ! in-plane correlations 
                    inpl_ind = shcloc(self%nrots, corrs, self%prev_corr) ! first improving in-plane index
                    if( inpl_ind > 0 ) inpl_corr = corrs( inpl_ind )     ! improving correlation found
                    self%nrefs_eval = self%nrefs_eval + 1                ! updates fractional search space
                endif
            end subroutine per_ref_srch

    end subroutine stochastic_srch_shc

    !>  \brief  stochastic neighborhood hill-climbing
    !! \param pftcc polarft corrcalc search storage
    !! \param iptcl particle index
    !! \param a,e search orientation
    !! \param lp low-pass cutoff freq
    !! \param szsn size ref evals
    subroutine stochastic_srch_snhc( self, pftcc, iptcl, a, e, lp, szsn )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(oris),             intent(inout) :: a, e
        real,                    intent(in)    :: lp
        integer,                 intent(in)    :: szsn
        real    :: projspace_corrs(self%nrefs), wcorr
        integer :: iref, isample
        if( nint(a%get(iptcl,'state')) > 0 )then
            ! initialize
            call self%prep4srch(pftcc, iptcl, a, e, lp)
            self%nbetter         = 0
            self%nrefs_eval      = 0
            self%proj_space_inds = 0
            projspace_corrs      = -1.
            ! search
            do isample=1,szsn
                iref = self%srch_order(isample) ! set the stochastic reference index
                call per_ref_srch(iref)         ! actual search
            end do
            self%nrefs_eval = szsn
            ! sort in correlation projection direction space
            call hpsort(self%nrefs, projspace_corrs, self%proj_space_inds) 
            ! output
            call self%prep_npeaks_oris
            call self%stochastic_weights(wcorr)
            call self%update_best(pftcc, iptcl, a)
        else
            call a%reject(iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED EXTREMAL SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind, state
                state = 1
                if( self%nstates > 1 ) state = nint( self%o_refs%get(iref, 'state') )
                if( self%state_exists(state) )then       
                    corrs     = pftcc%gencorrs(iref, iptcl) ! In-plane correlations
                    loc       = maxloc(corrs)               ! greedy in-plane
                    inpl_ind  = loc(1)                      ! in-plane angle index
                    inpl_corr = corrs(inpl_ind)             ! max in plane correlation
                    call self%store_solution(pftcc, iref, iref, inpl_ind, inpl_corr)
                    projspace_corrs( iref ) = inpl_corr     ! stash in-plane correlation for sorting
                endif
            end subroutine per_ref_srch

    end subroutine stochastic_srch_snhc

    !> stochastic search het
    !! \param pftcc polarft corrcalc search storage
    !! \param iptcl particle index
    !! \param a,e search orientation
    !! \param lp low-pass cutoff freq
    !! \param statecnt state-variable matrix
    !! \param extr_bound corr threshold
    subroutine stochastic_srch_het( self, pftcc, iptcl, a, e, extr_bound, statecnt)
        use simple_rnd, only: shcloc, irnd_uni
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(oris),             intent(inout) :: a, e
        real,                    intent(in)    :: extr_bound
        integer,                 intent(inout) :: statecnt(self%nstates)
        type(ori) :: o
        integer   :: iref, state
        real      :: corr, mi_state, frac, corrs(self%nstates)
        self%prev_state = nint(a%get(iptcl, 'state'))
        if( self%prev_state > self%nstates ) stop 'previous best state outside boundary; stochastic_srch_het; simple_prime3D_srch'
        if( .not. self%state_exists(self%prev_state) ) stop 'empty previous state; stochastic_srch_het; simple_prime3D_srch'
        if( self%prev_state > 0 )then
            ! initialize
            o = a%get_ori(iptcl)
            self%prev_roind = pftcc%get_roind(360.-o%e3get())
            self%prev_proj  = e%find_closest_proj(o,1)
            if(a%get(iptcl,'corr') < extr_bound)then
                ! state randomization
                statecnt(self%prev_state) = statecnt(self%prev_state) + 1
                self%nrefs_eval = 1
                state = irnd_uni(self%nstates)
                do while(state == self%prev_state .or. .not.self%state_exists(state))
                    state = irnd_uni(self%nstates)
                enddo
                iref = (state - 1) * self%nprojs + self%prev_proj
                corr = pftcc%corr(iref, iptcl, self%prev_roind)
            else
                ! SHC
                corrs = -1.
                do state = 1, self%nstates
                    if( .not.self%state_exists(state) )cycle
                    iref = (state-1) * self%nprojs + self%prev_proj 
                    corrs(state) = pftcc%corr(iref, iptcl, self%prev_roind)
                enddo
                self%prev_corr = corrs(self%prev_state)
                state          = shcloc(self%nstates, corrs, self%prev_corr)
                if(state == 0)state = self%prev_state ! numerical stability
                corr            = corrs(state)
                self%nrefs_eval = count(corrs <= self%prev_corr)
            endif
            ! updates orientation
            frac = 100.*real(self%nrefs_eval) / real(self%nstates)
            call o%set('frac', frac)
            call o%set('state', real(state))
            call o%set('corr', corr)
            call o%set('mi_proj', 1.)
            call o%set('mi_inpl', 1.)
            if( self%prev_state .ne. state )then
                mi_state = 0.
            else
                mi_state = 1.
            endif
            call o%set('mi_state', mi_state)
            call o%set('mi_joint', mi_state)
            call o%set('w', 1.)
            ! specscore
            iref = (state - 1) * self%nprojs + self%prev_proj
            call self%calc_specscore(pftcc, iref, iptcl, self%prev_roind)
            call o%set('specscore', self%specscore)
            ! updates orientations objects
            call self%o_peaks%set_ori(1, o)
            call a%set_ori(iptcl, o)
        else
            call a%reject(iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED HET SEARCH'
    end subroutine stochastic_srch_het

    !>  \brief  executes the in-plane search for discrete mode
    !! \param pftcc polarft corrcalc search storage
    !! \param iptcl particle index
    subroutine inpl_srch( self, pftcc, iptcl )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        real, allocatable :: cxy(:), crxy(:)
        type(ori) :: o
        real      :: cc
        integer   :: i, ref, inpl_ind
        if( self%doshift )then
            do i=self%nrefs,self%nrefs-self%npeaks+1,-1
                ref      = self%proj_space_inds( i )
                o        = self%o_refs%get_ori( ref )
                cc       = o%get('corr')
                inpl_ind = pftcc%get_roind( 360.-o%e3get() )
                if( self%greedy_inpl )then
                    call self%inplsrch_obj%set_indices(ref, iptcl)
                    crxy = self%inplsrch_obj%minimize(irot=inpl_ind)
                    if( crxy(1) >= cc )then
                        call o%set( 'corr', crxy(1) )
                        call o%e3set( 360.-crxy(2) )
                        call o%set_shift( crxy(3:4) )
                        call self%o_refs%set_ori( ref, o )
                    endif
                else
                    call self%shsrch_obj%set_indices(ref, iptcl, inpl_ind)
                    cxy = self%shsrch_obj%minimize()
                    if( cxy(1) >= cc )then
                        call o%set('corr', cxy(1))
                        call o%set_shift( cxy(2:3) )
                        call self%o_refs%set_ori( ref, o )
                    endif
                endif
            end do
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED INPL SEARCH'
    end subroutine inpl_srch

    !>  \brief  prepares reference indices for the search & fetches ctf
    !! \param pftcc polarft corrcalc search storage
    !! \param iptcl particle index
    !! \param a,e search orientation
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param target_projs grid projections
    subroutine prep4srch( self, pftcc, iptcl, a, e, lp, nnmat, target_projs )
        use simple_combinatorics, only: merge_into_disjoint_set
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(oris),             intent(inout) :: a, e
        real,                    intent(in)    :: lp
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn_static), target_projs(self%npeaks_grid)
        integer, allocatable :: nnvec(:)
        type(ori) :: o_prev
        real      :: cc_t_min_1, corr
        if( str_has_substr(self%refine,'neigh') )then
            if( .not. present(nnmat) )&
            &stop 'need optional nnmat to be present for refine=neigh modes :: prep4srch (prime3D_srch)'
            if( .not. present(target_projs) )&
            &stop 'need optional target_projs to be present for refine=neigh modes :: prep4srch (prime3D_srch)'
        endif
        o_prev          = a%get_ori(iptcl)
        self%prev_state = nint(o_prev%get('state'))                                    ! state index            
        self%prev_roind = pftcc%get_roind(360.-o_prev%e3get())                         ! in-plane angle index
        self%prev_shvec = o_prev%get_shift()                                           ! shift vector
        self%prev_proj  = e%find_closest_proj(o_prev,1)                                ! projection direction
        if( self%prev_state > self%nstates ) stop 'previous best state outside boundary; prep4srch; simple_prime3D_srch'
        if( .not. self%state_exists(self%prev_state) ) stop 'empty previous state; prep4srch; simple_prime3D_srch'
        select case( self%refine )
            case( 'no','shc','snhc','greedy' )                                         ! DISCRETE CASE
                call self%prep_reforis(e)                                              ! search space & order prep
                self%prev_ref = self%o_refs%find_closest_proj(o_prev, self%prev_state) ! find closest ori with same state
            case( 'neigh','shcneigh', 'greedyneigh' )                                  ! DISCRETE CASE WITH NEIGHBOURHOOD
                nnvec = merge_into_disjoint_set(self%nprojs, self%nnn_static, nnmat, target_projs) ! disjoint nearest neighbour set
                call self%prep_reforis(e, nnvec=nnvec)                                 ! search space & order prep
                self%prev_ref = self%o_refs%find_closest_proj(o_prev, self%prev_state) ! find closest ori with same state
            case( 'het' )
                self%prev_ref = (self%prev_state-1)*self%nprojs+self%prev_proj
            case( 'exp' )
                call self%prep_reforis(e, nnvec=nnmat(self%prev_proj,:))
                self%prev_ref = self%prev_proj ! because oris only cover one state
            case DEFAULT
                stop 'Unknown refinement mode; simple_prime3D_srch; prep4srch'
        end select
        ! prep corr
        if( self%refine .ne. 'het' )then
            corr = max( 0., pftcc%corr(self%prev_ref, iptcl, self%prev_roind) )
            if( corr - 1.0 > 1.0e-5 .or. .not. is_a_number(corr) )then
                print *, 'FLOATING POINT EXCEPTION ALARM; simple_prime3D_srch :: prep4srch'
                print *, 'corr > 1. or isNaN'
                print *, 'corr = ', corr
                if( corr > 1. )               corr = 1.
                if( .not. is_a_number(corr) ) corr = 0.
                call o_prev%print_ori()
            endif
            if( self%refine.eq.'no' .and. self%nstates==1 )then
                ! moving average for single state only
                cc_t_min_1 = -1.
                if( o_prev%isthere('corr') ) cc_t_min_1 = o_prev%get('corr')
                if( o_prev%isthere('lp') )then
                    if( abs(o_prev%get('lp') - lp) < 0.5 )then   ! previous and present correlations comparable
                        if( cc_t_min_1 > 0. )then
                            corr = 0.5 * corr + 0.5 * cc_t_min_1 ! diversifying limit
                        endif
                    endif
                endif
            endif
            self%prev_corr = corr
        endif
        DebugPrint '>>> PRIME3D_SRCH::PREPARED FOR SIMPLE_PRIME3D_SRCH'
    end subroutine prep4srch

    !>  \brief  prepares the search space (ref oris) & search order per particle
    !! \param e search orientation
    !! \param nnvec nearest neighbour matrix
    subroutine prep_reforis( self, e, nnvec )
        use simple_ran_tabu, only: ran_tabu
        class(prime3D_srch),  intent(inout) :: self
        class(oris),          intent(inout) :: e
        integer,    optional, intent(in)    :: nnvec(:)
        type(ran_tabu) :: rt
        integer        :: i, cnt, istate, iproj
        type(ori)      :: o
        if( trim(self%refine).eq.'exp' )then
            ! dev mode: search prev_state only with neighbourhood
            if(.not.present(nnvec))stop 'need neighbours for refine=exp; simple_prime3d_srch%prep_reforis'
            self%nnn        = size(nnvec)
            self%nnnrefs    = self%nnn
            self%srch_order = nnvec + (self%prev_state-1)*self%nprojs
            rt = ran_tabu(self%nnn)
            call rt%shuffle(self%srch_order)
        else
            ! dynamic update of number of nearest neighbours
            if( present(nnvec) )then
                self%nnn     = size(nnvec)
                self%nnnrefs =  self%nnn*self%nstates
            endif
            ! on exit all the oris are clean and only the out-of-planes, 
            ! state & proj fields are present
            if( str_has_substr(self%refine, 'neigh') )then ! local refinement
                allocate(self%srch_order(self%nnnrefs), source=0)
                rt = ran_tabu(self%nnnrefs)
            else
                allocate(self%srch_order(self%nrefs), source=0)
                rt = ran_tabu(self%nrefs)
            endif
            if( present(nnvec) )then
                do istate=0,self%nstates-1 ! concatenate nearest neighbor per state...
                    i = istate*self%nnn+1
                    self%srch_order(i:i+self%nnn-1) = nnvec + istate*self%nprojs
                enddo
                call rt%shuffle( self%srch_order ) ! ...& wizz it up
            else
                ! refine=no|shc
                call rt%ne_ran_iarr( self%srch_order )
            endif
        endif
        if( any(self%srch_order == 0) ) stop 'Invalid index in srch_order; simple_prime3d_srch::prep_ref_oris'
        ! prepare discrete reforis
        ! The underlying principle here is that o_refs is congruent with pftcc
        call self%o_refs%new( self%nrefs )          ! init references object
        cnt = 0
        do istate=1,self%nstates
            do iproj=1,self%nprojs
                cnt = cnt + 1
                o = e%get_ori( iproj )
                call o%set( 'state', real(istate) ) ! Updates state
                call o%set( 'proj', real(iproj) )   ! Updates proj
                call self%o_refs%set_ori( cnt,o )
            enddo
        enddo
        call rt%kill
    end subroutine prep_reforis

    ! CALCULATORS

    !>  \brief retrieves and preps npeaks orientations for reconstruction
    subroutine prep_npeaks_oris( self )
        use simple_sym, only: sym
        class(prime3D_srch),   intent(inout) :: self
        type(ori)  :: o, o_best, o_new
        type(oris) :: sym_os, o_peaks
        type(sym)  :: se
        real, allocatable :: corrs(:)
        real       :: euls(3), shvec(2)
        real       :: corr, frac, ang_sdev, dist
        integer    :: ipeak, cnt, ref, state, proj, loc(1)
        integer    :: neff_states ! number of effective (populated) states
        ! empty states
        neff_states = 1
        if(self%nstates > 1) neff_states = count(self%state_exists)
        ! init npeaks
        call o_peaks%new(self%npeaks)
        do ipeak = 1, self%npeaks
            ! get ipeak-th ori
            cnt = self%nrefs - self%npeaks + ipeak
            ref = self%proj_space_inds( cnt )
            if( ref < 1 .or. ref > self%nrefs )then
                print *, 'ref: ', ref
                stop 'ref index out of bound; simple_prime3D_srch::prep_npeaks_oris'
            endif
            o = self%o_refs%get_ori( ref )
            ! grab info
            state = nint( o%get('state') )
            if( .not. self%state_exists(state) )then
                print *,'empty state:',state,' ; simple_prime3D_srch::prep_npeaks_oris'
                stop
            endif
            proj  = nint( o%get('proj') )
            corr  = o%get('corr')
            if( .not. is_a_number(corr) ) stop 'correlation is NaN in simple_prime3D_srch::prep_npeaks_oris'
            euls  = o%get_euler()
            ! add shift
            shvec = self%prev_shvec
            if( self%doshift )shvec = shvec + o%get_shift()
            where( abs(shvec) < 1e-6 ) shvec = 0.
            ! copy info to new ori
            call o_new%new
            call o_new%set_euler( euls )  
            call o_new%set_shift( shvec )
            call o_new%set('state', real(state))
            call o_new%set('proj',  real(proj) )
            call o_new%set('corr',  corr       )
            ! stashes in self
            call o_peaks%set_ori( ipeak, o_new )
        enddo
        ! other variables
        if( str_has_substr(self%refine, 'neigh') )then
            frac = 100.*real(self%nrefs_eval) / real(self%nnn * neff_states)
        else if( trim(self%refine).eq.'exp' )then
            frac = 100.*real(self%nrefs_eval) / real(self%nnn) ! 1 state searched          
        else
            frac = 100.*real(self%nrefs_eval) / real(self%nprojs * neff_states)
        endif
        call o_peaks%set_all2single('frac',    frac )
        call o_peaks%set_all2single('mi_hard', 0.   )
        call o_peaks%set_all2single('dist',    0.   )
        ! angular standard deviation
        ang_sdev = 0.
        if( trim(self%pgrp).eq.'c1' )then
            ang_sdev = o_peaks%ang_sdev(self%refine, self%nstates, self%npeaks)
        else
            if( self%npeaks > 2 )then
                corrs  = o_peaks%get_all('corr')
                loc    = maxloc(corrs)
                o_best = o_peaks%get_ori(loc(1))
                call se%new(trim(self%pgrp))
                sym_os = o_peaks
                do ipeak = 1, self%npeaks
                    if(ipeak == loc(1))cycle
                    o = o_peaks%get_ori(ipeak)
                    call se%sym_euldist(o_best, o, dist)
                    call sym_os%set_ori(ipeak, o)
                enddo
                ang_sdev = sym_os%ang_sdev(self%refine, self%nstates, self%npeaks)
            endif         
        endif
        call o_peaks%set_all2single('sdev', ang_sdev)
        ! ctf parameters
        if( self%ctf.ne.'no' )then
            call o_peaks%set_all2single('dfx',   self%dfx   )
            call o_peaks%set_all2single('dfy',   self%dfy   )
            call o_peaks%set_all2single('angast',self%angast)
        endif
        ! the end
        self%o_peaks = o_peaks
        DebugPrint '>>> PRIME3D_SRCH::EXECUTED PREP_NPEAKS_ORIS'
    end subroutine prep_npeaks_oris

    !>  \brief  determines and updates stochastic weights
    !! \param pftcc polarft corrcalc search storage
    !! \param iptcl particle index
    !! \param e search orientation
    !! \param wcorr weights
    subroutine stochastic_weights( self, wcorr )
        class(prime3D_srch),     intent(inout) :: self
        real,                    intent(out)   :: wcorr
        real,    allocatable :: corrs(:), ws(:), logws(:)
        integer, allocatable :: order(:) 
        logical, allocatable :: included(:)
        integer              :: ipeak, alloc_stat
        if( self%npeaks == 1 )then
            call self%o_peaks%set(1,'ow',1.0)
            wcorr = self%o_peaks%get(1,'corr')
            return
        endif
        ! calculate the exponential of the negative distances
        ! so that when diff==0 the weights are maximum and when
        ! diff==corrmax the weights are minimum
        corrs = self%o_peaks%get_all('corr')
        ws    = exp(-(1.-corrs))
        logws = log(ws)
        order = (/(ipeak,ipeak=1,self%npeaks)/)
        call hpsort(self%npeaks, logws, order)
        call reverse(order)
        call reverse(logws)
        forall(ipeak=1:self%npeaks) ws(order(ipeak)) = exp(sum(logws(:ipeak))) 
        ! thresholding of the weights
        allocate(included(self%npeaks), source=.true.)
        included = (ws>=FACTWEIGHTS_THRESH)
        where( .not.included ) ws = 0.
        ! weighted corr
        wcorr = sum(ws*corrs,mask=included) / sum(ws,mask=included)
        ! update npeaks individual weights
        call self%o_peaks%set_all('ow', ws)
    end subroutine stochastic_weights

    !>  \brief  calculates the sprectral score
    subroutine calc_specscore( self, pftcc, iref, iptcl, roind )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iref, iptcl, roind
        real, allocatable :: frc(:)
        frc = pftcc%genfrc(iref, iptcl, roind)
        self%specscore = max(0., median_nocopy(frc))
    end subroutine calc_specscore

    ! GETTERS & SETTERS
    
    !>  \brief  to get the best orientation
    !! \param pftcc polarft corrcalc search storage
    !! \param iptcl particle index
    !! \param a search orientation
    subroutine update_best( self, pftcc, iptcl, a )
        use simple_sym,  only: sym
        use simple_math, only: myacos, rotmat2d, rad2deg
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(oris),             intent(inout) :: a
        type(sym)         :: se
        type(ori)         :: o_new, o_old, o_new_copy
        real, allocatable :: corrs(:)
        real              :: euldist, mi_joint, mi_proj, mi_inpl, mi_state, dist_inpl
        integer           :: roind, state, best_loc(1), iref
        o_old    = a%get_ori(iptcl)
        corrs    = self%o_peaks%get_all('corr')
        best_loc = maxloc(corrs) 
        o_new    = self%o_peaks%get_ori(best_loc(1))
        ! angular distances
        o_new_copy = o_new
        call se%new(trim(self%pgrp))
        call se%sym_euldist( o_old, o_new_copy, euldist )
        dist_inpl = rad2deg( o_new_copy.inplrotdist.o_old )
        call se%kill
        ! specscore & new parameters
        state = nint( o_new%get('state') )
        if( .not. self%state_exists(state) )then
            print *,'Empty state in simple_prime3d_srch; update_best'
            stop
        endif
        roind = pftcc%get_roind( 360.-o_new%e3get() )
        iref  = self%o_refs%find_closest_proj(o_new, state)
        call self%calc_specscore(pftcc, iref, iptcl, roind)
        ! calculate overlap between distributions
        mi_proj  = 0.
        mi_inpl  = 0.
        mi_state = 0.
        mi_joint = 0.
        if( euldist < 0.5 )then
            mi_proj = mi_proj + 1.
            mi_joint = mi_joint + 1.
        endif
        if( self%prev_roind == roind )then
            mi_inpl  = mi_inpl  + 1.
            mi_joint = mi_joint + 1.
        endif
        if( self%nstates > 1 )then
            if( self%prev_state == state )then
                mi_state = mi_state + 1.
                mi_joint = mi_joint + 1.
            endif
            mi_joint = mi_joint/3.
        else
            mi_joint = mi_joint/2.
        endif
        ! set the overlaps
        call a%set(iptcl, 'mi_proj',  mi_proj )
        call a%set(iptcl, 'mi_inpl',  mi_inpl )
        call a%set(iptcl, 'mi_state', mi_state)
        call a%set(iptcl, 'mi_joint', mi_joint)
        ! set the distances before we update the orientation
        call a%set(iptcl, 'dist', 0.5*euldist + 0.5*o_old%get('dist'))
        call a%set(iptcl, 'dist_inpl', dist_inpl)
        ! all the other stuff
        call a%set_euler(iptcl, o_new%get_euler()    )
        call a%set_shift(iptcl, o_new%get_shift()    )
        call a%set(iptcl, 'state', real(state)       )
        call a%set(iptcl, 'frac',  o_new%get('frac') )
        call a%set(iptcl, 'corr',  o_new%get('corr') )
        call a%set(iptcl, 'specscore', self%specscore)
        call a%set(iptcl, 'ow',    o_new%get('ow')   )
        call a%set(iptcl, 'mirr',  0.                )
        call a%set(iptcl, 'proj',  o_new%get('proj') )
        call a%set(iptcl, 'sdev',  o_new%get('sdev') )
        ! stash and return
        o_new = a%get_ori(iptcl)
        call self%o_peaks%set_ori(best_loc(1), o_new)
        DebugPrint '>>> PRIME3D_SRCH::GOT BEST ORI'
    end subroutine update_best

    !>  \brief  to get one orientation
    subroutine get_ori( self, ipeak, o2update )
        class(prime3D_srch), intent(inout) :: self
        integer,             intent(in)    :: ipeak    !< which peak
        class(ori),          intent(inout) :: o2update !< search orientation
        real      :: euls(3), x, y, rstate, ow
        integer   :: ind, npeaks
        if( str_has_substr(self%refine,'shc') .and. ipeak /= 1 ) stop 'get_ori not for shc-modes; simple_prime3D_srch'
        npeaks = self%o_peaks%get_noris()
        if( ipeak < 1 .or. ipeak > npeaks ) stop 'Invalid index in simple_prime3D_srch::get_ori'
        ind    = npeaks - ipeak + 1
        euls   = self%o_peaks%get_euler( ind )
        x      = self%o_peaks%get( ind, 'x'    )
        y      = self%o_peaks%get( ind, 'y'    )
        rstate = self%o_peaks%get( ind, 'state')
        if( .not. self%state_exists(nint(rstate)) )then
            print *,'Empty state in simple_prime3d_srch; get_ori'
            stop
        endif
        ow = self%o_peaks%get( ind, 'ow' )
        call o2update%set_euler( euls )
        call o2update%set( 'x',     x      )
        call o2update%set( 'y',     y      )
        call o2update%set( 'state', rstate )
        call o2update%set( 'ow',    ow     )
    end subroutine get_ori

    !>  \brief  to produce oris object for probabilistic reconstruction
    subroutine get_oris( self, os, o )
        class(prime3D_srch), intent(inout) :: self
        class(oris),         intent(inout) :: os !< search orientation list
        class(ori),          intent(inout) :: o  !< search orientation
        type(ori) :: o_new
        integer   :: ipeak, npeaks
        npeaks = self%o_peaks%get_noris()
        call os%new( npeaks )
        do ipeak=1,npeaks
            o_new = o
            call self%get_ori(ipeak, o_new)
            call os%set_ori(ipeak, o_new)
        enddo
    end subroutine get_oris

    !>  \brief store_solution stash solution in o_refs
    !! \param pftcc  polarft corrcalc search storage
    !! \param ind proj index
    !! \param ref ref index
    !! \param inpl_ind in-plane index
    !! \param corr solution
    subroutine store_solution( self, pftcc, ind, ref, inpl_ind, corr )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: ind, ref, inpl_ind
        real,                    intent(in)    :: corr
        real :: e3
        self%proj_space_inds( ind ) = ref          ! stash reference similarly to other search modes
        e3 = 360. - pftcc%get_rot( inpl_ind )      ! psi
        call self%o_refs%set( ref, 'e3',   e3   )  ! stash psi
        call self%o_refs%set( ref, 'corr', corr )  ! stash correlation
        if( self%npeaks==1 ) call self%o_refs%set( ref, 'ow', 1. ) ! set reconstruction weight
    end subroutine store_solution

    ! GETTERS FOR TESTING

    !>  \brief  is for creating the search order
    function get_srch_order( self )result( inds )
        class(prime3D_srch), intent(inout) :: self
        integer,allocatable :: inds(:)
        integer             :: alloc_stat
        allocate( inds(size(self%srch_order)), stat=alloc_stat )
        call alloc_errchk( 'simple_prime3D_srch::get_srch_order', alloc_stat)
        inds(:) = self%srch_order(:)
    end function get_srch_order

    !> \brief  for setting o_peaks
    subroutine set_o_peaks( self, oris_in )
        class(prime3D_srch), intent(inout) :: self
        class(oris),         intent(in)    :: oris_in
        self%o_peaks = oris_in
        DebugPrint '>>> PRIME3D_SRCH::EXECUTED SET_O_PEAKS'
    end subroutine set_o_peaks

    !>  \brief  is for getting o_peaks
    function get_o_peaks( self )result( out_os )
        class(prime3D_srch), intent(inout) :: self
        type(oris) :: out_os
        out_os = self%o_peaks
    end function get_o_peaks

    !>  \brief  is for getting the n-last reference orientations
    function get_o_refs( self, n )result( out_os )
        class(prime3D_srch), intent(inout) :: self
        integer, optional,   intent(in)    :: n !<  nth-last reference
        type(oris) :: out_os
        integer    :: i, cnt, in
        if( present(n) )then
            if( n < 1 )then
                stop 'invalid index in prime3D_srch::get_o_refs, 1'
            elseif( n>self%nrefs)then
                stop 'invalid index in prime3D_srch::get_o_refs, 3'
            endif
            in = n
        else
            in = self%nrefs
        endif
        out_os = oris(n)
        cnt = 0
        do i=self%nrefs,self%nrefs-n+1,-1
            cnt = cnt+1
            call out_os%set_ori( i, self%o_refs%get_ori(i) )
        enddo
    end function get_o_refs

    function get_ntotrefs( self )result( nrefs )
        class(prime3D_srch), intent(inout) :: self
        integer :: nrefs
        nrefs = self%nrefs
    end function get_ntotrefs

    function get_nrefs( self )result( nrefs )
        class(prime3D_srch), intent(inout) :: self
        integer :: nrefs
        if( str_has_substr(self%refine,'neigh') )then
            nrefs = self%nnn_static * self%nstates
        else
            nrefs = self%nrefs
        endif
    end function get_nrefs

    function get_nrots( self )result( nrots )
        class(prime3D_srch), intent(inout) :: self
        integer :: nrots
        nrots = self%nrots
    end function get_nrots

    function get_prevshvec( self )result( vec )
        class(prime3D_srch), intent(inout) :: self
        real :: vec(2)
        vec = self%prev_shvec
    end function get_prevshvec

    function get_prevcorr( self )result( corr )
        class(prime3D_srch), intent(inout) :: self
        real :: corr
        corr = self%prev_corr
    end function get_prevcorr

    function get_prevroind( self )result( roind )
        class(prime3D_srch), intent(inout) :: self
        integer :: roind
        roind = self%prev_roind
    end function get_prevroind

    function get_prevstate( self )result( state )
        class(prime3D_srch), intent(inout) :: self
        integer :: state
        state = self%prev_state
    end function get_prevstate

    function get_prevref( self )result( ref )
        class(prime3D_srch), intent(inout) :: self
        integer :: ref
        ref = self%prev_ref
    end function get_prevref

    function get_angast( self )result( ang )
        class(prime3D_srch), intent(inout) :: self
        real :: ang
        ang = self%angast
    end function get_angast

    function get_dfx( self )result( dfx )
        class(prime3D_srch), intent(inout) :: self
        real :: dfx
        dfx = self%dfx
    end function get_dfx

    function get_dfy( self )result( dfy )
        class(prime3D_srch), intent(inout) :: self
        real :: dfy
        dfy = self%dfy
    end function get_dfy

    subroutine set_ctf( self, ctf )
        class(prime3D_srch), intent(inout) :: self
        character(len=*), intent(in) :: ctf
        self%ctf = ctf
    end subroutine set_ctf

    subroutine set_greedy_inpl( self, log )
        class(prime3D_srch), intent(inout) :: self
        logical, intent(in) :: log
        self%greedy_inpl = log
    end subroutine set_greedy_inpl

    ! MEMORY MANAGEMENT

    subroutine online_allocate( self )
        class(prime3D_srch), intent(inout) :: self
        integer :: alloc_stat
        allocate(self%proj_space_inds(self%nrefs), stat=alloc_stat)  
        call alloc_errchk('In: prime3D_srch_allocate; simple_prime3D_srch', alloc_stat)
        self%proj_space_inds = 0
    end subroutine online_allocate

    subroutine online_destruct( self )
        class(prime3D_srch), intent(inout) :: self
        if( allocated(self%proj_space_inds)) deallocate(self%proj_space_inds)
        if( allocated(self%srch_order )    ) deallocate(self%srch_order     )
        call self%o_refs%kill
    end subroutine online_destruct

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(prime3D_srch), intent(inout) :: self !< instance
        if( self%exists )then
            call self%o_peaks%kill
            call self%online_destruct
            if( allocated(self%state_exists) ) deallocate(self%state_exists)
            self%exists = .false.
        endif        
    end subroutine kill

end module simple_prime3D_srch
