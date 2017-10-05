! PRIME3D stochastic search routines
module simple_prime3D_srch
#include "simple_lib.f08"
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_shc_inplane,      only: shc_inplane
use simple_pftcc_shsrch,     only: pftcc_shsrch
use simple_pftcc_inplsrch,   only: pftcc_inplsrch
implicit none

public :: prime3D_srch
private
#include "simple_local_flags.inc"

real,    parameter :: FACTWEIGHTS_THRESH = 0.001 !< threshold for factorial weights
real,    parameter :: E3HALFWINSZ        = 60.   !< in-plane angle half window size

!> struct for prime3d params
type prime3D_srch
    private
    class(polarft_corrcalc), pointer :: pftcc_ptr => null()      !< pointer to pftcc (corrcalc) object
    class(oris),             pointer :: a_ptr     => null()      !< pointer to b%a (primary particle orientation table)
    class(oris),             pointer :: e_ptr     => null()      !< pointer to b%e (reference orientations)
    type(oris)                       :: o_refs                   !< projection directions search space
    type(oris)                       :: o_peaks                  !< orientations of best npeaks oris
    type(shc_inplane)                :: shcgrid                  !< in-plane grid search object
    type(pftcc_shsrch)               :: shsrch_obj               !< origin shift search object
    type(pftcc_inplsrch)             :: inplsrch_obj             !< in-plane search object
    integer                          :: iptcl          = 0       !< global particle index
    integer                          :: nrefs          = 0       !< total # references (nstates*nprojs)
    integer                          :: nnnrefs        = 0       !< total # neighboring references (nstates*nnn)
    integer                          :: nstates        = 0       !< # states
    integer                          :: nprojs         = 0       !< # projections (same as # oris in o_refs)
    integer                          :: nrots          = 0       !< # in-plane rotations in polar representation
    integer                          :: npeaks         = 0       !< # peaks (nonzero orientation weights)
    integer                          :: npeaks_grid    = 0       !< # peaks after coarse search
    integer                          :: nbetter        = 0       !< # better orientations identified
    integer                          :: nrefs_eval     = 0       !< # references evaluated
    integer                          :: nnn_static     = 0       !< # nearest neighbors (static)
    integer                          :: nnn            = 0       !< # nearest neighbors (dynamic)
    integer                          :: ntn            = 1       !< # of time neighbors
    integer                          :: nsym           = 0       !< symmetry order
    integer                          :: prev_roind     = 0       !< previous in-plane rotation index
    integer                          :: prev_state     = 0       !< previous state index
    integer                          :: prev_ref       = 0       !< previous reference index
    integer                          :: prev_proj      = 0       !< previous projection index
    integer                          :: kstop_grid     = 0       !< Frequency limit of first coarse grid search
    real                             :: prev_corr      = 1.      !< previous best correlation
    real                             :: specscore      = 0.      !< spectral score
    real                             :: prev_shvec(2)  = 0.      !< previous origin shift vector
    real                             :: athres         = 0.      !< angular threshold
    real                             :: dfx            = 0.      !< ctf x-defocus
    real                             :: dfy            = 0.      !< ctf y-defocus
    real                             :: angast         = 0.      !< ctf astigmatism
    integer, allocatable             :: proj_space_inds(:)       !< projection space index array
    integer, allocatable             :: srch_order(:)            !< stochastic search order
    logical, allocatable             :: state_exists(:)          !< indicates whether each state is populated
    character(len=STDLEN)            :: refine         = ''      !< refinement flag
    character(len=STDLEN)            :: ctf            = ''      !< ctf flag
    character(len=STDLEN)            :: shbarr         = ''      !< shift barrier flag
    character(len=STDLEN)            :: pgrp           = 'c1'    !< point-group symmetry
    logical                          :: dev            = .false. !< development flag
    logical                          :: doshift        = .true.  !< 2 indicate whether 2 serch shifts
    logical                          :: greedy_inpl    = .true.  !< 2 indicate whether in-plane search is greedy or not
    logical                          :: inpl_neigh     = .false. !< 2 indicate whether to use implane neighbourhood search
    logical                          :: exists         = .false. !< 2 indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SEARCH ROUTINES
    procedure          :: exec_prime3D_srch
    procedure          :: exec_prime3D_srch_het
    procedure, private :: greedy_srch
    procedure, private :: greedy_subspace_srch
    procedure, private :: stochastic_srch
    procedure, private :: stochastic_srch_shc
    procedure, private :: stochastic_srch_snhc
    procedure, private :: stochastic_srch_het
    procedure, private :: stochastic_srch_bystate
    procedure          :: inpl_srch
    procedure, private :: inpl_grid_srch
    procedure, private :: inpl_grid_srch_exhaustive
    ! PREPARATION ROUTINES
    procedure          :: prep4srch
    procedure          :: prep_reforis
    ! CALCULATORS
    procedure          :: prep_npeaks_oris
    procedure          :: stochastic_weights
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
    procedure          :: get_npeaks
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
    subroutine new( self, iptcl, pftcc, a, e, p )
        use simple_params, only: params
        class(prime3D_srch),             intent(inout) :: self   !< instance
        integer,                         intent(in)    :: iptcl  !< global particle index
        class(polarft_corrcalc), target, intent(inout) :: pftcc  !< correlator
        class(oris),             target, intent(in)    :: a      !< primary particle orientation table
        class(oris),             target, intent(in)    :: e      !< reference orientations
        class(params),                   intent(in)    :: p      !< parameters
        integer :: alloc_stat, nstates_eff
        real    :: lims(2,2), lims_init(2,2)
        ! destroy possibly pre-existing instance
        call self%kill
        ! set constants
        self%pftcc_ptr   => pftcc
        self%a_ptr       => a
        self%e_ptr       => e
        self%iptcl       =  iptcl
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
        if( str_has_substr(self%refine,'shc') )then
            if( self%npeaks > 1 ) stop 'npeaks must be equal to 1 with refine=shc|shcneigh'
            self%greedy_inpl = .false.
        endif
        self%dev = .false.
        if( p%dev .eq. 'yes' ) self%dev = .true.
        ! construct composites
        if( p%oritab.ne.'' )then
            self%state_exists = self%a_ptr%states_exist(self%nstates)
        else
            allocate(self%state_exists(self%nstates), stat=alloc_stat)
            allocchk('In: new; simple_prime3D_srch, 1')
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
                case('states')
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
        ! in-plane neighbourhood
        select case(trim(self%refine))
        case('tseries')
            self%inpl_neigh = .true.
        case DEFAULT
            self%inpl_neigh = .false.
        end select
        ! generate oris oject in which the best npeaks refs will be stored
        call self%o_peaks%new(self%npeaks)
        ! updates option to search shift
        self%doshift = p%doshift
        ! create in-plane search objects
        lims(:,1)      = -p%trs
        lims(:,2)      =  p%trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        call self%shcgrid%new
        call self%shsrch_obj%new(  pftcc, lims, lims_init=lims_init, shbarrier=self%shbarr)
        call self%inplsrch_obj%new(pftcc, lims, lims_init=lims_init, shbarrier=self%shbarr)
        self%exists = .true.
        DebugPrint '>>> PRIME3D_SRCH::CONSTRUCTED NEW SIMPLE_PRIME3D_SRCH OBJECT'
    end subroutine new

    ! SEARCH ROUTINES

    !>  \brief  exec_prime3D_srch is a master prime search routine
    !! \param lp low-pass cutoff freq
    !! \param greedy greedy search flag
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs projection indices for grid search
    !! \param szsn optional size for snhc refinement
    subroutine exec_prime3D_srch( self, lp, greedy, nnmat, grid_projs, szsn )
        class(prime3D_srch), intent(inout) :: self
        real,                intent(in)    :: lp
        logical, optional,   intent(in)    :: greedy
        integer, optional,   intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:), szsn
        logical :: ggreedy
        ggreedy = .false.
        if( present(greedy) ) ggreedy = greedy
        call self%online_allocate
        if( ggreedy )then
            if( trim(self%refine).eq.'states' )then
                call self%stochastic_srch_bystate(lp, nnmat)
            else
                call self%greedy_srch(lp, nnmat, grid_projs)
            endif
        else if( self%refine.eq.'snhc' )then
            if( .not. present(szsn) )then
                stop 'refine=snhc mode needs optional input szsn; simple_prime3D_srch :: exec_prime3D_srch'
            endif
            call self%stochastic_srch_snhc(lp, szsn)
        else if( str_has_substr(self%refine,'shc') )then
            call self%stochastic_srch_shc(lp, nnmat, grid_projs)
        else
            call self%stochastic_srch(lp, nnmat, grid_projs)
        endif
        call self%online_destruct
        DebugPrint '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SRCH'
    end subroutine exec_prime3D_srch

    !>  \brief state labeler
    subroutine exec_prime3D_srch_het( self, extr_bound, statecnt )
        class(prime3D_srch), intent(inout) :: self
        real,                intent(in)    :: extr_bound
        integer,             intent(inout) :: statecnt(self%nstates)
        call self%stochastic_srch_het(extr_bound, statecnt)
        call self%online_destruct
        DebugPrint '>>> PRIME3D_SRCH::EXECUTED PRIME3D_HET_SRCH'
    end subroutine exec_prime3D_srch_het

    !>  \brief  Individual stochastic search by state
    !! \param pftcc polarft corrcalc search storage
    !! \param a,e search orientation
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs grid projections
    subroutine stochastic_srch_bystate( self, lp, nnmat )
        class(prime3D_srch),     intent(inout) :: self
        real,                    intent(in)    :: lp
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn_static)
        real      :: projspace_corrs(self%nrefs), wcorr
        integer   :: iref, isample
        ! execute search
        if( nint(self%a_ptr%get(self%iptcl,'state')) > 0 )then
            ! initialize
            call self%prep4srch(lp, nnmat=nnmat)
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
            call self%inpl_srch ! search shifts
            ! prepare weights and orientations
            call self%prep_npeaks_oris
            call self%stochastic_weights(wcorr)
            call self%update_best
            call self%a_ptr%set(self%iptcl, 'corr', wcorr)
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED STOCHASTIC BYSTATE SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind, state
                state = nint(self%o_refs%get(iref, 'state'))
                if(state .ne. self%prev_state )then
                    ! checker that will need to go
                    print *,self%iptcl, self%prev_ref, self%prev_proj, self%prev_state, iref, state
                    stop 'srch order error; simple_prime3d_srch; stochastic_srch_bystate'
                endif
                corrs     = self%pftcc_ptr%gencorrs(iref, self%iptcl) ! In-plane correlations
                loc       = maxloc(corrs)               ! greedy in-plane
                inpl_ind  = loc(1)                      ! in-plane angle index
                inpl_corr = corrs(inpl_ind)             ! max in plane correlation
                call self%store_solution(iref, iref, inpl_ind, inpl_corr)
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

    end subroutine stochastic_srch_bystate

    !>  \brief  greedy hill-climbing
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs projection indices for grid search
    subroutine greedy_srch( self, lp, nnmat, grid_projs )
        class(prime3D_srch), intent(inout) :: self
        real,                intent(in)    :: lp
        integer, optional,   intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:)
        real    :: projspace_corrs(self%nrefs),wcorr
        integer :: iref,isample,nrefs,target_projs(self%npeaks_grid)
        if( nint(self%a_ptr%get(self%iptcl,'state')) > 0 )then
            if( str_has_substr(self%refine, 'neigh') )then
                ! for neighbour modes we do a coarse grid search first
                if( .not. present(grid_projs) ) stop 'need optional grid_projs 4 subspace srch; prime3D_srch :: greedy_srch'
                call self%greedy_subspace_srch(grid_projs, target_projs)
                ! initialize
                call self%prep4srch(lp, nnmat, target_projs)
                nrefs = self%nnnrefs
            else
                ! initialize
                call self%prep4srch(lp)
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
            call self%inpl_srch
            ! prepare weights & orientation
            call self%prep_npeaks_oris
            call self%stochastic_weights(wcorr)
            call self%update_best
            call self%a_ptr%set(self%iptcl, 'corr', wcorr)
        else
            call self%a_ptr%reject(self%iptcl)
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
                    corrs     = self%pftcc_ptr%gencorrs(iref, self%iptcl) ! In-plane correlations
                    loc       = maxloc(corrs)                        ! greedy in-plane
                    inpl_ind  = loc(1)                               ! in-plane angle index
                    inpl_corr = corrs(inpl_ind)                      ! max in plane correlation
                    call self%store_solution(iref, iref, inpl_ind, inpl_corr)
                    projspace_corrs( iref ) = inpl_corr              ! stash in-plane correlation for sorting
                endif
            end subroutine per_ref_srch

    end subroutine greedy_srch

    !>  \brief  greedy hill-climbing (4 initialisation)
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs projection indices for grid search
    !! \param target_projs target projections (output)
    subroutine greedy_subspace_srch( self, grid_projs, target_projs )
        class(prime3D_srch), intent(inout) :: self
        integer,             intent(in)    :: grid_projs(:)
        integer,             intent(inout) :: target_projs(:)
        real      :: projspace_corrs(self%nrefs)
        integer   :: iref, isample, nrefs, ntargets, cnt, prev_proj, istate
        integer   :: state_cnt(self%nstates), iref_state
        type(ori) :: o_prev
        if( nint(self%a_ptr%get(self%iptcl,'state')) > 0 )then
            o_prev    = self%a_ptr%get_ori(self%iptcl)
            prev_proj = self%e_ptr%find_closest_proj(o_prev,1)
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
            call self%a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED GREEDY SUBSPACE SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind
                if( self%state_exists(istate) )then
                    corrs     = self%pftcc_ptr%gencorrs(iref, self%iptcl, self%kstop_grid) ! In-plane correlations
                    loc       = maxloc(corrs)                                         ! greedy in-plane
                    inpl_ind  = loc(1)                                                ! in-plane angle index
                    inpl_corr = corrs(inpl_ind)                                       ! max in plane correlation
                    projspace_corrs( iref )      = inpl_corr                          ! stash in-plane correlation for sorting
                    self%proj_space_inds( iref ) = iref                               ! stash the index for sorting
                endif
            end subroutine per_ref_srch

    end subroutine greedy_subspace_srch

    !>  \brief  executes the stochastic soft orientation search
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs projection indices for grid search
    subroutine stochastic_srch( self, lp, nnmat, grid_projs )
        class(prime3D_srch), intent(inout) :: self
        real,                intent(in)    :: lp
        integer, optional,   intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:)
        integer, allocatable :: roind_vec(:)    ! in-plane neighbourhood
        real                 :: projspace_corrs(self%nrefs),wcorr,e3_prev
        integer              :: iref,isample,nrefs,target_projs(self%npeaks_grid)
        ! execute search
        if( nint(self%a_ptr%get(self%iptcl,'state')) > 0 )then
            if( str_has_substr(self%refine, 'neigh') )then
                ! for neighbour modes we do a coarse grid search first
                if( .not. present(grid_projs) ) stop 'need optional grid_projs 4 subspace srch; prime3D_srch :: stochastic_srch'
                call self%greedy_subspace_srch(grid_projs, target_projs)
                ! initialize
                call self%prep4srch(lp, nnmat, target_projs)
                nrefs = self%nnnrefs
            else
                ! initialize
                call self%prep4srch(lp)
                nrefs = self%nrefs
            endif
            if( self%inpl_neigh )then
                ! in-plane neighbourhood
                e3_prev = 360. - self%pftcc_ptr%get_rot(self%prev_roind)
                roind_vec = self%pftcc_ptr%get_win_roind(e3_prev, E3HALFWINSZ)
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
            call self%inpl_srch ! search shifts
            ! prepare weights and orientations
            call self%prep_npeaks_oris
            call self%stochastic_weights(wcorr)
            call self%update_best
            call self%a_ptr%set(self%iptcl, 'corr', wcorr)
        else
            call self%a_ptr%reject(self%iptcl)
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
                    ! In-plane correlations
                    if( self%inpl_neigh )then
                        corrs = self%pftcc_ptr%gencorrs(iref, self%iptcl, roind_vec=roind_vec)
                    else
                        corrs = self%pftcc_ptr%gencorrs(iref, self%iptcl)
                    endif
                    loc       = maxloc(corrs)   ! greedy in-plane
                    inpl_ind  = loc(1)          ! in-plane angle index
                    inpl_corr = corrs(inpl_ind) ! max in plane correlation
                    call self%store_solution(iref, iref, inpl_ind, inpl_corr)
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
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs projection indices for grid search
    subroutine stochastic_srch_shc( self, lp, nnmat, grid_projs )
        class(prime3D_srch), intent(inout) :: self
        real,                intent(in)    :: lp
        integer, optional,   intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:)
        real    :: inpl_corr,corrs(self%nrots)
        integer :: iref,isample,nrefs,inpl_ind,loc(1),target_projs(self%npeaks_grid)
        logical :: found_better
        if( nint(self%a_ptr%get(self%iptcl,'state')) > 0 )then
            if( str_has_substr(self%refine, 'neigh') )then
                ! for neighbour modes we do a coarse grid search first
                if( .not. present(grid_projs) ) stop 'need optional grid_projs 4 subspace srch; prime3D_srch :: stochastic_srch_shc'
                call self%greedy_subspace_srch(grid_projs, target_projs)
                ! initialize
                call self%prep4srch(lp, nnmat, target_projs)
                nrefs = self%nnnrefs
            else
                ! initialize
                call self%prep4srch(lp)
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
                    call self%store_solution(self%nrefs, iref, inpl_ind, inpl_corr)
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
                call self%store_solution(self%nrefs, iref, inpl_ind, inpl_corr)
            endif
            ! search shifts
            call self%inpl_srch
            ! output
            call self%prep_npeaks_oris
            call self%update_best
        else
            call self%a_ptr%reject(self%iptcl)
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
                    corrs    = self%pftcc_ptr%gencorrs( iref, self%iptcl ) ! in-plane correlations
                    inpl_ind = shcloc(self%nrots, corrs, self%prev_corr) ! first improving in-plane index
                    if( inpl_ind > 0 ) inpl_corr = corrs( inpl_ind )     ! improving correlation found
                    self%nrefs_eval = self%nrefs_eval + 1                ! updates fractional search space
                endif
            end subroutine per_ref_srch

    end subroutine stochastic_srch_shc

    !>  \brief  stochastic neighborhood hill-climbing
    !! \param lp low-pass cutoff freq
    !! \param szsn size ref evals
    subroutine stochastic_srch_snhc( self, lp, szsn )
        class(prime3D_srch), intent(inout) :: self
        real,                intent(in)    :: lp
        integer,             intent(in)    :: szsn
        real    :: projspace_corrs(self%nrefs), wcorr
        integer :: iref, isample
        if( nint(self%a_ptr%get(self%iptcl,'state')) > 0 )then
            ! initialize
            call self%prep4srch(lp)
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
            call self%update_best
        else
            call self%a_ptr%reject(self%iptcl)
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
                    corrs     = self%pftcc_ptr%gencorrs(iref, self%iptcl) ! In-plane correlations
                    loc       = maxloc(corrs)               ! greedy in-plane
                    inpl_ind  = loc(1)                      ! in-plane angle index
                    inpl_corr = corrs(inpl_ind)             ! max in plane correlation
                    call self%store_solution(iref, iref, inpl_ind, inpl_corr)
                    projspace_corrs( iref ) = inpl_corr     ! stash in-plane correlation for sorting
                endif
            end subroutine per_ref_srch

    end subroutine stochastic_srch_snhc

    !> stochastic search het
    !! \param lp low-pass cutoff freq
    !! \param statecnt state counter array
    !! \param extr_bound corr threshold
    subroutine stochastic_srch_het( self, extr_bound, statecnt)
        use simple_rnd, only: shcloc, irnd_uni
        class(prime3D_srch), intent(inout) :: self
        real,                intent(in)    :: extr_bound
        integer,             intent(inout) :: statecnt(self%nstates)
        type(ori) :: o
        integer   :: iref, state
        real      :: corr, mi_state, frac, corrs(self%nstates)
        self%prev_state = nint(self%a_ptr%get(self%iptcl, 'state'))
        if( self%prev_state > self%nstates ) stop 'previous best state outside boundary; stochastic_srch_het; simple_prime3D_srch'
        if( self%prev_state > 0 )then
            if( .not. self%state_exists(self%prev_state) ) stop 'empty previous state; stochastic_srch_het; simple_prime3D_srch'
            ! initialize
            o = self%a_ptr%get_ori(self%iptcl)
            self%prev_roind = self%pftcc_ptr%get_roind(360.-o%e3get())
            self%prev_proj  = self%e_ptr%find_closest_proj(o,1)
            if( self%a_ptr%get(self%iptcl,'corr') < extr_bound)then
                ! state randomization
                statecnt(self%prev_state) = statecnt(self%prev_state) + 1
                self%nrefs_eval = 1
                state = irnd_uni(self%nstates)
                do while(state == self%prev_state .or. .not.self%state_exists(state))
                    state = irnd_uni(self%nstates)
                enddo
                iref = (state - 1) * self%nprojs + self%prev_proj
                corr = self%pftcc_ptr%corr(iref, self%iptcl, self%prev_roind)
            else
                ! SHC
                corrs = -1.
                do state = 1, self%nstates
                    if( .not.self%state_exists(state) )cycle
                    iref = (state-1) * self%nprojs + self%prev_proj
                    corrs(state) = self%pftcc_ptr%corr(iref, self%iptcl, self%prev_roind)
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
            ! updates orientations objects
            call self%o_peaks%set_ori(1, o)
            call self%a_ptr%set_ori(self%iptcl, o)
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED HET SEARCH'
    end subroutine stochastic_srch_het

    !>  \brief  executes the in-plane search. This improved routine took HCN
    !!          from stall @ 5.4 to close to convergence @ 4.5 after 10 iters
    !!          refine=no with 1000 referecnes
    subroutine inpl_srch( self )
        class(prime3D_srch), intent(inout) :: self
        integer, allocatable :: inpl_inds(:)
        real,    allocatable :: cxy(:), crxy(:), shvecs(:,:)
        type(ori) :: o
        real      :: cc, e3
        integer   :: i, ref
        if( self%doshift )then
            call self%inpl_grid_srch(inpl_inds, shvecs)
            do i=self%nrefs,self%nrefs-self%npeaks+1,-1
                ref = self%proj_space_inds( i )
                o   = self%o_refs%get_ori( ref )
                cc  = self%o_refs%get( ref, 'corr' )
                if( self%greedy_inpl )then
                    call self%inplsrch_obj%set_indices(ref, self%iptcl)
                    crxy = self%inplsrch_obj%minimize(irot=inpl_inds(i), shvec=shvecs(i,:))
                    if( crxy(1) >= cc )then
                        call o%set( 'corr', crxy(1) )
                        call o%e3set( 360.-crxy(2) )
                        call o%set_shift( crxy(3:4) )
                        call self%o_refs%set_ori( ref, o )
                    endif
                else
                    call self%shsrch_obj%set_indices(ref, self%iptcl, inpl_inds(i))
                    cxy = self%shsrch_obj%minimize(shvec=shvecs(i,:))
                    if( cxy(1) >= cc )then
                        e3 = 360. - self%pftcc_ptr%get_rot(inpl_inds(i)) ! psi
                        call o%e3set(e3)                                 ! stash psi
                        call o%set('corr', cxy(1))
                        call o%set_shift( cxy(2:3) )
                        call self%o_refs%set_ori( ref, o )
                    endif
                endif
            end do
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED INPL SEARCH'
    end subroutine inpl_srch

    !>  \brief  discrete stochastic hill-climbing for the in-plane search
    subroutine inpl_grid_srch( self, inpl_inds, shvecs )
        class(prime3D_srch),    intent(inout) :: self
        integer, allocatable,   intent(out)   :: inpl_inds(:)
        real,    allocatable,   intent(out)   :: shvecs(:,:)
        integer           :: i, ref, istop, prev_rot
        if( allocated(inpl_inds) ) deallocate(inpl_inds)
        if( allocated(shvecs)    ) deallocate(shvecs)
        istop = self%nrefs - self%npeaks + 1
        allocate( shvecs(istop:self%nrefs,2), inpl_inds(istop:self%nrefs) )
        do i=self%nrefs,istop,-1
            ref      = self%proj_space_inds(i)
            prev_rot = self%pftcc_ptr%get_roind( 360.-self%o_refs%e3get(ref) )
            call self%shcgrid%srch( self%pftcc_ptr, ref, self%iptcl, self%nrots, prev_rot, inpl_inds(i), shvecs(i,:))
        end do
        DebugPrint '>>> PRIME3D_SRCH::FINISHED INPL GRID SEARCH'
    end subroutine inpl_grid_srch

    !>  \brief  implemented to test wheter we had an in-plane search deficiency
    !!          test on HCN took it from 5.4 A to 4.8 A with Nspace=1000 and refine=no
    !!          but it is too slow to be feasible. Implemented the above routine for
    !!          SHC-based grid search in the plane followed by more focused
    !!          continuous simplex refinement
    subroutine inpl_grid_srch_exhaustive( self, inpl_inds, shvecs )
        class(prime3D_srch),    intent(inout) :: self
        integer, allocatable,   intent(out)   :: inpl_inds(:)
        real,    allocatable,   intent(out)   :: shvecs(:,:)
        real    :: cc, cc_best, xsh, ysh
        integer :: i, j, jrot, ref, inpl_ind, istop, n_inpl_changes, n_trs_changes
        if( allocated(inpl_inds) ) deallocate(inpl_inds)
        if( allocated(shvecs)    ) deallocate(shvecs)
        istop = self%nrefs - self%npeaks + 1
        allocate( shvecs(istop:self%nrefs,2), inpl_inds(istop:self%nrefs) )
        n_inpl_changes = 0
        do i=self%nrefs,istop,-1
            ref      = self%proj_space_inds( i )
            inpl_ind = self%pftcc_ptr%get_roind( 360.-self%o_refs%e3get(ref) )
            cc_best  = -1.0
            do j=inpl_ind-SHC_INPL_INPLHWDTH,inpl_ind+SHC_INPL_INPLHWDTH
                jrot = cyci_1d([1,self%nrots], j)
                xsh  = -SHC_INPL_TRSHWDTH
                do while( xsh <= SHC_INPL_TRSHWDTH )
                    ysh = -SHC_INPL_TRSHWDTH
                    do while( ysh <= SHC_INPL_TRSHWDTH )
                        cc = self%pftcc_ptr%corr(ref, self%iptcl, jrot, [xsh,ysh])
                        if( cc > cc_best )then
                            cc_best      = cc
                            inpl_inds(i) = jrot
                            shvecs(i,:)  = [xsh,ysh]
                        endif
                        ysh = ysh + SHC_INPL_TRSSTEPSZ
                    end do
                    xsh = xsh + SHC_INPL_TRSSTEPSZ
                end do
            end do
            if( inpl_inds(i) /= inpl_ind )    n_inpl_changes = n_inpl_changes + 1
            if( sum(abs(shvecs(i,:))) > 0.1 ) n_trs_changes  = n_trs_changes  + 1
        end do
        call self%a_ptr%set(self%iptcl, 'inpl_changes', real(n_inpl_changes)/real(self%nrefs - istop + 1))
        call self%a_ptr%set(self%iptcl, 'trs_changes',  real(n_inpl_changes)/real(self%nrefs - istop + 1))
        DebugPrint '>>> PRIME3D_SRCH::FINISHED INPL GRID SEARCH'
    end subroutine inpl_grid_srch_exhaustive

    !>  \brief  prepares reference indices for the search & fetches ctf
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param target_projs projection indices for grid search
    subroutine prep4srch( self, lp, nnmat, target_projs )
        use simple_combinatorics, only: merge_into_disjoint_set
        class(prime3D_srch), intent(inout) :: self
        real,                intent(in)    :: lp
        integer, optional,   intent(in)    :: nnmat(self%nprojs,self%nnn_static), target_projs(self%npeaks_grid)
        integer, allocatable :: nnvec(:)
        real,    allocatable :: frc(:)
        type(ori) :: o_prev
        real      :: cc_t_min_1, corr
        if( str_has_substr(self%refine,'neigh') )then
            if( .not. present(nnmat) )&
            &stop 'need optional nnmat to be present for refine=neigh modes :: prep4srch (prime3D_srch)'
            if( .not. present(target_projs) )&
            &stop 'need optional target_projs to be present for refine=neigh modes :: prep4srch (prime3D_srch)'
        endif
        o_prev          = self%a_ptr%get_ori(self%iptcl)
        self%prev_state = nint(o_prev%get('state'))                                    ! state index
        self%prev_roind = self%pftcc_ptr%get_roind(360.-o_prev%e3get())                ! in-plane angle index
        self%prev_shvec = o_prev%get_2Dshift()                                         ! shift vector
        self%prev_proj  = self%e_ptr%find_closest_proj(o_prev,1)                       ! projection direction
        if( self%prev_state > self%nstates ) stop 'previous best state outside boundary; prep4srch; simple_prime3D_srch'
        if( self%prev_state > 0 )then
            if( .not. self%state_exists(self%prev_state) ) stop 'empty previous state; prep4srch; simple_prime3D_srch'
        endif
        select case( self%refine )
            case( 'no','shc','snhc','greedy','tseries' )                               ! DISCRETE CASE
                call self%prep_reforis                                                 ! search space & order prep
                self%prev_ref = self%o_refs%find_closest_proj(o_prev, self%prev_state) ! find closest ori with same state
            case( 'neigh','shcneigh', 'greedyneigh' )                                  ! DISCRETE CASE WITH NEIGHBOURHOOD
                nnvec = merge_into_disjoint_set(self%nprojs, self%nnn_static, nnmat, target_projs) ! disjoint nearest neighbour set
                call self%prep_reforis(nnvec=nnvec)                                    ! search space & order prep
                self%prev_ref = self%o_refs%find_closest_proj(o_prev, self%prev_state) ! find closest ori with same state
            case( 'het' )
                self%prev_ref = (self%prev_state-1)*self%nprojs + self%prev_proj
            case( 'states' )
                call self%prep_reforis(nnvec=nnmat(self%prev_proj,:))
                self%prev_ref = self%o_refs%find_closest_proj(o_prev, self%prev_state)
            case DEFAULT
                stop 'Unknown refinement mode; simple_prime3D_srch; prep4srch'
        end select
        ! calc specscore
        frc = self%pftcc_ptr%genfrc(self%prev_ref, self%iptcl, self%prev_roind)
        self%specscore = max(0., median_nocopy(frc))
        ! prep corr
        if( self%refine .ne. 'het' )then
            corr = max( 0., self%pftcc_ptr%corr(self%prev_ref, self%iptcl, self%prev_roind) )
            if( corr - 1.0 > 1.0e-5 .or. .not. is_a_number(corr) )then
                print *, 'FLOATING POINT EXCEPTION ALARM; simple_prime3D_srch :: prep4srch'
                print *, 'corr > 1. or isNaN'
                print *, 'corr = ', corr
                if( corr > 1. )               corr = 1.
                if( .not. is_a_number(corr) ) corr = 0.
                call o_prev%print_ori()
            endif
            self%prev_corr = corr
        endif
        DebugPrint '>>> PRIME3D_SRCH::PREPARED FOR SIMPLE_PRIME3D_SRCH'
    end subroutine prep4srch

    !>  \brief  prepares the search space (ref oris) & search order per particle
    !! \param nnvec nearest neighbour array
    subroutine prep_reforis( self, nnvec )
        use simple_ran_tabu, only: ran_tabu
        class(prime3D_srch),  intent(inout) :: self
        integer,    optional, intent(in)    :: nnvec(:)
        type(ran_tabu) :: rt
        integer        :: i, cnt, istate, iproj
        type(ori)      :: o
        ! prepare discrete reforis
        ! The underlying principle here is that o_refs is congruent with pftcc
        call self%o_refs%new( self%nrefs )          ! init references object
        cnt = 0
        do istate=1,self%nstates
            do iproj=1,self%nprojs
                cnt = cnt + 1
                o = self%e_ptr%get_ori( iproj )
                call o%set( 'state', real(istate) ) ! Updates state
                call o%set( 'proj', real(iproj) )   ! Updates proj
                call self%o_refs%set_ori( cnt,o )
            enddo
        enddo
        ! dynamic update of number of nearest neighbours
        if( present(nnvec) )then
            self%nnn = size(nnvec)
            if( trim(self%refine) .eq. 'states')then
                self%nnnrefs = self%nnn
            else
                self%nnnrefs =  self%nnn*self%nstates
            endif
        endif
        if( str_has_substr(self%refine, 'neigh') )then ! local refinement
            allocate(self%srch_order(self%nnnrefs), source=0)
            rt = ran_tabu(self%nnnrefs)
        else if( trim(self%refine).eq.'states' )then
            allocate(self%srch_order(self%nnn), source=0)
            rt = ran_tabu(self%nnn)
        else
            allocate(self%srch_order(self%nrefs), source=0)
            if( trim(self%refine).ne.'tseries' )rt = ran_tabu(self%nrefs)
        endif
        if( present(nnvec) )then
            if( trim(self%refine).eq.'states' )then
                self%srch_order = nnvec + (self%prev_state-1)*self%nprojs
            else
                do istate=0,self%nstates-1 ! concatenate nearest neighbor per state...
                    i = istate*self%nnn+1
                    self%srch_order(i:i+self%nnn-1) = nnvec + istate*self%nprojs
                enddo
            endif
            call rt%shuffle( self%srch_order ) ! ...& wizz it up
        else
            select case( trim(self%refine) )
            case('no','shc','snhc','greedy')
                call rt%ne_ran_iarr( self%srch_order )
            case('tseries')
                call prep_tseries_srchorder
            case DEFAULT
                stop 'Unknown refinement mode; simple_prime3d_srch%prep_reforis'
            end select
        endif
        if( any(self%srch_order == 0) ) stop 'Invalid index in srch_order; simple_prime3d_srch::prep_ref_oris'
        ! cleanup
        call rt%kill

        contains

            subroutine prep_tseries_srchorder
                type(oris)        :: trefs
                type(ori)         :: oref, o, o_mirr, optcl
                real, allocatable :: tdists(:)
                real       :: somereal, dists(self%nrefs)
                integer    :: n_trefs, itref, iref, istart, iend
                ! does not make sense to have multiple states support
                optcl = self%a_ptr%get_ori(self%iptcl)
                ! particle range
                istart  = max(1, self%iptcl-self%ntn)
                iend    = min(self%a_ptr%get_noris(), self%iptcl+self%ntn)
                n_trefs = iend-istart+1
                call trefs%new( n_trefs )
                ! particles range object
                cnt = 0
                do itref = istart, iend
                    cnt    = cnt + 1
                    if( itref.eq.self%iptcl )then
                        o = optcl
                    else
                        o      = self%a_ptr%get_ori(itref)
                        o_mirr = o
                        call o_mirr%mirror2d
                        if((optcl.euldist.o) > (optcl.euldist.o_mirr))then
                            o = o_mirr
                        else
                            ! all good
                        endif
                    endif
                    call trefs%set_ori(cnt, o)
                enddo
                ! particles-to-references distances
                do iref = 1, self%nrefs
                    self%srch_order(iref) = iref
                    if(iref .eq. self%prev_ref)then
                        dists(iref) = huge(somereal) ! guarantees previous best last
                    else
                        oref = self%o_refs%get_ori(iref)
                        call trefs%calc_euldists(oref, tdists)
                        dists(iref) = sum(tdists)
                    endif
                enddo
                ! sorting
                call hpsort(self%nrefs, dists, self%srch_order)
            end subroutine prep_tseries_srchorder

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
                print *, 'empty state: ', state
                stop 'simple_prime3D_srch::prep_npeaks_oris'
            endif
            proj  = nint( o%get('proj') )
            corr  = o%get('corr')
            if( .not. is_a_number(corr) ) stop 'correlation is NaN in simple_prime3D_srch::prep_npeaks_oris'
            euls  = o%get_euler()
            ! add shift
            shvec = self%prev_shvec
            if( self%doshift )shvec = shvec + o%get_2Dshift()
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
        else if( trim(self%refine).eq.'states' )then
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
    !! \param wcorr weights
    subroutine stochastic_weights( self, wcorr )
        class(prime3D_srch),     intent(inout) :: self
        real,                    intent(out)   :: wcorr
        real, allocatable :: corrs(:)
        real              :: ws(self%npeaks), logws(self%npeaks)
        integer           :: order(self%npeaks), ipeak
        logical           :: included(self%npeaks)
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
        included = (ws >= FACTWEIGHTS_THRESH)
        where( .not.included ) ws = 0.
        ! weighted corr
        wcorr = sum(ws*corrs,mask=included) / sum(ws,mask=included)
        ! update npeaks individual weights
        call self%o_peaks%set_all('ow', ws)
    end subroutine stochastic_weights

    ! GETTERS & SETTERS

    !>  \brief  to get the best orientation
    subroutine update_best( self )
        use simple_sym,  only: sym
        class(prime3D_srch), intent(inout) :: self
        type(sym)         :: se
        type(ori)         :: o_new, o_old, o_new_copy
        real, allocatable :: corrs(:)
        real              :: euldist, mi_joint, mi_proj, mi_inpl, mi_state, dist_inpl
        integer           :: roind, state, best_loc(1)
        o_old    = self%a_ptr%get_ori(self%iptcl)
        corrs    = self%o_peaks%get_all('corr')
        best_loc = maxloc(corrs)
        o_new    = self%o_peaks%get_ori(best_loc(1))
        ! angular distances
        o_new_copy = o_new
        call se%new(trim(self%pgrp))
        call se%sym_euldist( o_old, o_new_copy, euldist )
        dist_inpl = rad2deg( o_new_copy.inplrotdist.o_old )
        call se%kill
        state = nint( o_new%get('state') )
        if( .not. self%state_exists(state) )then
            print *, 'empty state: ', state
            stop 'simple_prime3d_srch; update_best'
        endif
        roind = self%pftcc_ptr%get_roind( 360.-o_new%e3get() )
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
        call self%a_ptr%set(self%iptcl, 'mi_proj',  mi_proj )
        call self%a_ptr%set(self%iptcl, 'mi_inpl',  mi_inpl )
        call self%a_ptr%set(self%iptcl, 'mi_state', mi_state)
        call self%a_ptr%set(self%iptcl, 'mi_joint', mi_joint)
        ! set the distances before we update the orientation
        call self%a_ptr%set(self%iptcl, 'dist', 0.5*euldist + 0.5*o_old%get('dist'))
        call self%a_ptr%set(self%iptcl, 'dist_inpl', dist_inpl)
        ! all the other stuff
        call self%a_ptr%set_euler(self%iptcl, o_new%get_euler()    )
        call self%a_ptr%set_shift(self%iptcl, o_new%get_2Dshift()    )
        call self%a_ptr%set(self%iptcl, 'state', real(state)       )
        call self%a_ptr%set(self%iptcl, 'frac',  o_new%get('frac') )
        call self%a_ptr%set(self%iptcl, 'corr',  o_new%get('corr') )
        call self%a_ptr%set(self%iptcl, 'specscore', self%specscore)
        call self%a_ptr%set(self%iptcl, 'ow',    o_new%get('ow')   )
        call self%a_ptr%set(self%iptcl, 'mirr',  0.                )
        call self%a_ptr%set(self%iptcl, 'proj',  o_new%get('proj') )
        call self%a_ptr%set(self%iptcl, 'sdev',  o_new%get('sdev') )
        ! stash and return
        o_new = self%a_ptr%get_ori(self%iptcl)
        call self%o_peaks%set_ori(best_loc(1), o_new)
        DebugPrint '>>> PRIME3D_SRCH::GOT BEST ORI'
    end subroutine update_best

    !>  \brief  to get one orientation
    subroutine get_ori( self, ipeak, o2update )
        class(prime3D_srch), intent(inout) :: self
        integer,             intent(in)    :: ipeak    !< which peak
        class(ori),          intent(inout) :: o2update !< search orientation
        real      :: euls(3), x, y, rstate, rproj, corr, ow
        integer   :: npeaks
        if( str_has_substr(self%refine,'shc') .and. ipeak /= 1 ) stop 'get_ori not for shc-modes; simple_prime3D_srch'
        npeaks = self%o_peaks%get_noris()
        if( ipeak < 1 .or. ipeak > npeaks ) stop 'Invalid index in simple_prime3D_srch::get_ori'
        euls   = self%o_peaks%get_euler( ipeak )
        x      = self%o_peaks%get( ipeak, 'x'    )
        y      = self%o_peaks%get( ipeak, 'y'    )
        rstate = self%o_peaks%get( ipeak, 'state')
        if( allocated(self%state_exists) )then
            if( .not. self%state_exists(nint(rstate)) )then
                print *, 'empty state: ', nint(rstate)
                stop 'simple_prime3d_srch; get_ori'
            endif
        endif
        rproj = self%o_peaks%get( ipeak, 'proj' )
        corr  = self%o_peaks%get( ipeak, 'corr' )
        ow    = self%o_peaks%get( ipeak, 'ow'   )
        call o2update%set_euler( euls )
        call o2update%set_shift([x, y])
        call o2update%set( 'state', rstate )
        call o2update%set( 'proj',  rproj  )
        call o2update%set( 'corr',  corr   )
        call o2update%set( 'ow',    ow     )
    end subroutine get_ori

    !>  \brief  to produce oris object for probabilistic reconstruction
    subroutine get_oris( self, os, o_in )
        class(prime3D_srch), intent(inout) :: self
        class(oris),         intent(out)   :: os    !< search orientation list
        class(ori),          intent(in)    :: o_in  !< search orientation
        type(ori) :: o
        integer   :: ipeak!, npeaks
        !npeaks = self%o_peaks%get_noris()
        call os%new( self%npeaks )
        do ipeak=1,self%npeaks
            o = o_in
            call self%get_ori(ipeak, o)
            call os%set_ori(ipeak, o)
        enddo
    end subroutine get_oris

    !>  \brief store_solution stash solution in o_refs
    !! \param ind proj index
    !! \param ref ref index
    !! \param inpl_ind in-plane index
    !! \param corr solution
    subroutine store_solution( self, ind, ref, inpl_ind, corr )
        class(prime3D_srch),     intent(inout) :: self
        integer,                 intent(in)    :: ind, ref, inpl_ind
        real,                    intent(in)    :: corr
        real :: e3
        self%proj_space_inds( ind ) = ref              ! stash reference similarly to other search modes
        e3 = 360. - self%pftcc_ptr%get_rot( inpl_ind ) ! psi
        call self%o_refs%set( ref, 'e3',   e3   )      ! stash psi
        call self%o_refs%set( ref, 'corr', corr )      ! stash correlation
        if( self%npeaks==1 ) call self%o_refs%set( ref, 'ow', 1. ) ! set reconstruction weight
    end subroutine store_solution

    ! GETTERS FOR TESTING

    !>  \brief  is for getting the search order
    function get_srch_order( self )result( inds )
        class(prime3D_srch), intent(inout) :: self
        integer,allocatable :: inds(:)
        allocate( inds(size(self%srch_order)), stat=alloc_stat )
        allocchk( 'simple_prime3D_srch::get_srch_order')
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

    integer function get_ntotrefs( self )
        class(prime3D_srch), intent(inout) :: self
        get_ntotrefs = self%nrefs
    end function get_ntotrefs

    integer function get_nrefs( self )
        class(prime3D_srch), intent(inout) :: self
        if( str_has_substr(self%refine,'neigh') )then
            get_nrefs = self%nnn_static * self%nstates
        else
            get_nrefs = self%nrefs
        endif
    end function get_nrefs

    integer function get_nrots( self )
        class(prime3D_srch), intent(inout) :: self
        get_nrots = self%nrots
    end function get_nrots

    integer function get_npeaks( self )
        class(prime3D_srch), intent(inout) :: self
        get_npeaks = self%npeaks
    end function get_npeaks

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
        allocate(self%proj_space_inds(self%nrefs), stat=alloc_stat)
        allocchk('In: prime3D_srch_allocate; simple_prime3D_srch')
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
            self%a_ptr => null()
            self%e_ptr => null()
            call self%o_peaks%kill
            call self%online_destruct
            call self%shcgrid%kill
            if( allocated(self%state_exists) ) deallocate(self%state_exists)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_prime3D_srch
