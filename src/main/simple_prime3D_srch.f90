! PRIME3D stochastic search routines
module simple_prime3D_srch
#include "simple_lib.f08"
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_sym,              only: sym
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_pftcc_shsrch,     only: pftcc_shsrch
implicit none

public :: prime3D_srch, associate_prime3D_srch_ptrs, nullify_prime3D_srch_ptrs
private
#include "simple_local_flags.inc"

real, parameter :: FACTWEIGHTS_THRESH = 0.001 !< threshold for factorial weights

! pointers to objects shared across threads
class(polarft_corrcalc), pointer :: pftcc_ptr                => null()  !< corrcalc) object
class(oris),             pointer :: a_ptr                    => null()  !< b%a (primary particle orientation table)
class(oris),             pointer :: e_ptr                    => null()  !< b%e (reference orientations)
class(sym),              pointer :: se_ptr                   => null()  !< b%se (symmetry elements)
type(oris),              pointer :: o_refs_ptr(:)            => null()  !< for storing search space and solutiuons
type(oris),              pointer :: o_peaks_ptr(:)           => null()  !< corrr peak oris for soft assignment
integer,                 pointer :: srch_order_ptr(:,:)      => null()  !< reference search order
integer,                 pointer :: proj_space_inds_ptr(:,:) => null()  !< projection indices
logical,                 pointer :: state_exists_ptr(:)      => null()  !< state existence flags

type prime3D_srch
    private
    type(pftcc_shsrch)               :: shsrch_obj               !< origin shift search object
    integer                          :: iptcl          = 0       !< global particle index
    integer                          :: nrefs          = 0       !< total # references (nstates*nprojs)
    integer                          :: nnnrefs        = 0       !< total # neighboring references (nstates*nnn)
    integer                          :: nstates        = 0       !< # states
    integer                          :: nprojs         = 0       !< # projections
    integer                          :: nrots          = 0       !< # in-plane rotations in polar representation
    integer                          :: npeaks         = 0       !< # peaks (nonzero orientation weights)
    integer                          :: npeaks_grid    = 0       !< # peaks after coarse search
    integer                          :: nbetter        = 0       !< # better orientations identified
    integer                          :: nrefs_eval     = 0       !< # references evaluated
    integer                          :: nnn_static     = 0       !< # nearest neighbors (static)
    integer                          :: nnn            = 0       !< # nearest neighbors (dynamic)
    integer                          :: prev_roind     = 0       !< previous in-plane rotation index
    integer                          :: prev_state     = 0       !< previous state index
    integer                          :: prev_ref       = 0       !< previous reference index
    integer                          :: prev_proj      = 0       !< previous projection index
    integer                          :: kstop_grid     = 0       !< Frequency limit of first coarse grid search
    real                             :: prev_corr      = 1.      !< previous best correlation
    real                             :: specscore      = 0.      !< spectral score
    real                             :: prev_shvec(2)  = 0.      !< previous origin shift vector
    character(len=STDLEN)            :: refine         = ''      !< refinement flag
    character(len=STDLEN)            :: shbarr         = ''      !< shift barrier flag
    logical                          :: doshift        = .true.  !< 2 indicate whether 2 serch shifts
    logical                          :: greedy_inpl    = .true.  !< 2 indicate whether in-plane search is greedy or not
    logical                          :: exists         = .false. !< 2 indicate existence
  contains
    procedure          :: new
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
    procedure          :: prep4srch
    procedure          :: prep_npeaks_oris_and_weights
    procedure          :: store_solution
end type prime3D_srch

contains

    subroutine associate_prime3D_srch_ptrs(pftcc, a, e, se, o_refs, o_peaks, srch_order, proj_space_inds, state_exists) 
        class(polarft_corrcalc), target, intent(in) :: pftcc                !< corrcalc) object
        class(oris),             target, intent(in) :: a                    !< b%a (primary particle orientation table)
        class(oris),             target, intent(in) :: e                    !< b%e (reference orientations)
        class(sym),              target, intent(in) :: se                   !< b%se (symmetry elements)
        type(oris),              target, intent(in) :: o_refs(:)            !< for storing search space and solutiuons
        type(oris),              target, intent(in) :: o_peaks(:)           !< corrr peak oris for soft assignment
        integer,                 target, intent(in) :: srch_order(:,:)      !< reference search order
        integer,                 target, intent(in) :: proj_space_inds(:,:) !< projection indices
        logical,                 target, intent(in) :: state_exists(:)      !< state existence flags
        pftcc_ptr           => pftcc
        a_ptr               => a
        e_ptr               => e
        se_ptr              => se
        o_refs_ptr          => o_refs
        o_peaks_ptr         => o_peaks
        srch_order_ptr      => srch_order
        proj_space_inds_ptr => proj_space_inds
        state_exists_ptr    => state_exists
    end subroutine associate_prime3D_srch_ptrs

    subroutine new( self, iptcl, p )
        use simple_params, only: params
        class(prime3D_srch), intent(inout) :: self   !< instance
        integer,             intent(in)    :: iptcl  !< global particle index
        class(params),       intent(in)    :: p      !< parameters
        integer :: alloc_stat, nstates_eff
        real    :: lims(2,2), lims_init(2,2)
        ! set constants
        self%iptcl       =  iptcl
        self%nstates     =  p%nstates
        self%nprojs      =  p%nspace
        self%nrefs       =  self%nprojs*self%nstates
        self%nrots       =  round2even(twopi*real(p%ring2))
        self%npeaks      =  p%npeaks
        self%nbetter     =  0
        self%nrefs_eval  =  0
        self%doshift     =  p%doshift
        self%refine      =  p%refine
        self%nnn_static  =  p%nnn
        self%nnn         =  p%nnn
        self%nnnrefs     =  self%nnn*self%nstates
        self%shbarr      =  p%shbarrier
        self%kstop_grid  =  p%kstop_grid
        self%greedy_inpl = .true.
        if( str_has_substr(self%refine,'shc') )then
            if( self%npeaks > 1 ) stop 'npeaks must be equal to 1 with refine=shc|shcneigh'
            self%greedy_inpl = .false.
        endif
        ! multiple states
        if( self%nstates == 1 )then
            self%npeaks_grid = GRIDNPEAKS
        else
            ! number of populated states
            nstates_eff = count(state_exists_ptr)
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
        ! updates option to search shift
        self%doshift = p%doshift
        ! create in-plane search objects
        lims(:,1)      = -p%trs
        lims(:,2)      =  p%trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        call self%shsrch_obj%new(pftcc_ptr, lims, lims_init=lims_init, shbarrier=self%shbarr, nrestarts=3, maxits=60)
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
        DebugPrint '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SRCH'
    end subroutine exec_prime3D_srch

    !>  \brief state labeler
    subroutine exec_prime3D_srch_het( self, extr_bound, statecnt )
        class(prime3D_srch), intent(inout) :: self
        real,                intent(in)    :: extr_bound
        integer,             intent(inout) :: statecnt(self%nstates)
        call self%stochastic_srch_het(extr_bound, statecnt)
        DebugPrint '>>> PRIME3D_SRCH::EXECUTED PRIME3D_HET_SRCH'
    end subroutine exec_prime3D_srch_het

    !>  \brief  Individual stochastic search by state
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    subroutine stochastic_srch_bystate( self, lp, nnmat )
        class(prime3D_srch),     intent(inout) :: self
        real,                    intent(in)    :: lp
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn_static)
        real      :: projspace_corrs(self%nrefs), wcorr
        integer   :: iref, isample
        ! execute search
        if( nint(a_ptr%get(self%iptcl,'state')) > 0 )then
            ! initialize
            call self%prep4srch(lp, nnmat=nnmat)
            self%nbetter         = 0
            self%nrefs_eval      = 0
            projspace_corrs      = -1.
            ! search
            do isample = 1, self%nnn
                iref = srch_order_ptr(self%iptcl, isample)        ! set the stochastic reference index
                call per_ref_srch(iref)                ! actual search
                if( self%nbetter >= self%npeaks ) exit ! exit condition
            end do
            ! sort in correlation projection direction space
            call hpsort(self%nrefs, projspace_corrs, proj_space_inds_ptr(self%iptcl,:))
            call self%inpl_srch ! search shifts
            ! prepare weights and orientations
            call self%prep_npeaks_oris_and_weights
        else
            call a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED STOCHASTIC BYSTATE SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind, state
                state = nint(o_refs_ptr(self%iptcl)%get(iref, 'state'))
                if(state .ne. self%prev_state )then
                    ! checker that will need to go
                    print *,self%iptcl, self%prev_ref, self%prev_proj, self%prev_state, iref, state
                    stop 'srch order error; simple_prime3d_srch; stochastic_srch_bystate'
                endif
                corrs     = pftcc_ptr%gencorrs(iref, self%iptcl) ! In-plane correlations
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
        if( nint(a_ptr%get(self%iptcl,'state')) > 0 )then
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
            projspace_corrs      = -1.
            ! search
            do isample=1,nrefs
                iref = srch_order_ptr(self%iptcl,isample) ! set the reference index
                call per_ref_srch(iref)         ! actual search
            end do
            ! in greedy mode, we evaluate all refs
            self%nrefs_eval = nrefs
            ! sort in correlation projection direction space
            call hpsort(self%nrefs, projspace_corrs, proj_space_inds_ptr(self%iptcl,:))
            ! take care of the in-planes
            call self%inpl_srch
            ! prepare weights & orientation
            call self%prep_npeaks_oris_and_weights
        else
            call a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED GREEDY SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind, state
                state = 1
                if( self%nstates > 1 ) state = nint( o_refs_ptr(self%iptcl)%get(iref, 'state') )
                if( state_exists_ptr(state) )then
                    corrs     = pftcc_ptr%gencorrs(iref, self%iptcl) ! In-plane correlations
                    loc       = maxloc(corrs)                        ! greedy in-plane
                    inpl_ind  = loc(1)                               ! in-plane angle index
                    inpl_corr = corrs(inpl_ind)                      ! max in plane correlation
                    call self%store_solution(iref, iref, inpl_ind, inpl_corr)
                    projspace_corrs( iref ) = inpl_corr              ! stash in-plane correlation for sorting
                endif
            end subroutine per_ref_srch

    end subroutine greedy_srch

    !>  \brief  greedy hill-climbing (4 initialisation)
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
        if( nint(a_ptr%get(self%iptcl,'state')) > 0 )then
            prev_proj = e_ptr%find_closest_proj(a_ptr%get_ori(self%iptcl),1)
            ! initialize
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
            call hpsort(self%nrefs, projspace_corrs, proj_space_inds_ptr(self%iptcl,:))
            ! return target points
            ntargets = size(target_projs)
            cnt = 1
            target_projs( cnt ) = prev_proj ! previous always part of the targets
            if( self%nstates == 1 )then
                ! Single state
                do isample=self%nrefs,self%nrefs - ntargets + 1,-1
                    if( target_projs(1) == proj_space_inds_ptr(self%iptcl,isample) )then
                        ! direction is already in target set
                    else
                        cnt = cnt + 1
                        target_projs(cnt) = proj_space_inds_ptr(self%iptcl,isample)
                        if( cnt == ntargets ) exit
                    endif
                end do
            else
                ! Multiples states
                state_cnt = 1                                           ! previous always part of the targets
                do isample = self%nrefs, 1, -1
                    if( cnt >= self%npeaks_grid )exit                   ! all that we need
                    iref_state = proj_space_inds_ptr(self%iptcl,isample)    ! reference index to multi-state space
                    istate     = ceiling(real(iref_state)/real(self%nprojs))
                    iref       = iref_state - (istate-1)*self%nprojs    ! reference index to single state space
                    if( .not.state_exists_ptr(istate) )cycle
                    if( any(target_projs == iref) )cycle                ! direction is already set
                    if( state_cnt(istate) >= GRIDNPEAKS )cycle          ! state is already filled
                    cnt = cnt + 1
                    target_projs(cnt) = iref
                    state_cnt(istate) = state_cnt(istate) + 1
                end do
            endif
        else
            call a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED GREEDY SUBSPACE SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind
                if( state_exists_ptr(istate) )then
                    corrs     = pftcc_ptr%gencorrs(iref, self%iptcl, self%kstop_grid) ! In-plane correlations
                    loc       = maxloc(corrs)                                         ! greedy in-plane
                    inpl_ind  = loc(1)                                                ! in-plane angle index
                    inpl_corr = corrs(inpl_ind)                                       ! max in plane correlation
                    projspace_corrs( iref )      = inpl_corr                          ! stash in-plane correlation for sorting
                    proj_space_inds_ptr(self%iptcl,iref)  = iref                               ! stash the index for sorting
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
        real                 :: projspace_corrs(self%nrefs),wcorr,e3_prev
        integer              :: iref,isample,nrefs,target_projs(self%npeaks_grid)
        ! execute search
        if( nint(a_ptr%get(self%iptcl,'state')) > 0 )then
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
            ! initialize, ctd
            self%nbetter         = 0
            self%nrefs_eval      = 0
            projspace_corrs      = -1.
            do isample=1,nrefs
                iref = srch_order_ptr(self%iptcl,isample)        ! set the stochastic reference index
                call per_ref_srch(iref)                ! actual search
                if( self%nbetter >= self%npeaks ) exit ! exit condition
            end do
            ! sort in correlation projection direction space
            call hpsort(self%nrefs, projspace_corrs, proj_space_inds_ptr(self%iptcl,:))
            call self%inpl_srch ! search shifts
            ! prepare weights and orientations
            call self%prep_npeaks_oris_and_weights
        else
            call a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind, state
                state = 1
                if( self%nstates > 1 ) state = nint( o_refs_ptr(self%iptcl)%get(iref, 'state') )
                if( state_exists_ptr(state) )then
                    ! In-plane correlations
                    corrs     = pftcc_ptr%gencorrs(iref, self%iptcl)
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
        if( nint(a_ptr%get(self%iptcl,'state')) > 0 )then
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
            ! search
            do isample=1,nrefs
                iref = srch_order_ptr(self%iptcl, isample ) ! stochastic reference index
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
            call self%prep_npeaks_oris_and_weights
        else
            call a_ptr%reject(self%iptcl)
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
                if( self%nstates > 1 ) state = nint( o_refs_ptr(self%iptcl)%get(iref, 'state') )
                if( state_exists_ptr( state ) )then
                    corrs    = pftcc_ptr%gencorrs( iref, self%iptcl ) ! in-plane correlations
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
        if( nint(a_ptr%get(self%iptcl,'state')) > 0 )then
            ! initialize
            call self%prep4srch(lp)
            self%nbetter         = 0
            self%nrefs_eval      = 0
            projspace_corrs      = -1.
            ! search
            do isample=1,szsn
                iref = srch_order_ptr(self%iptcl,isample) ! set the stochastic reference index
                call per_ref_srch(iref)         ! actual search
            end do
            self%nrefs_eval = szsn
            ! sort in correlation projection direction space
            call hpsort(self%nrefs, projspace_corrs, proj_space_inds_ptr(self%iptcl,:))
            ! output
            call self%prep_npeaks_oris_and_weights
        else
            call a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED EXTREMAL SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind, state
                state = 1
                if( self%nstates > 1 ) state = nint( o_refs_ptr(self%iptcl)%get(iref, 'state') )
                if( state_exists_ptr(state) )then
                    corrs     = pftcc_ptr%gencorrs(iref, self%iptcl) ! In-plane correlations
                    loc       = maxloc(corrs)               ! greedy in-plane
                    inpl_ind  = loc(1)                      ! in-plane angle index
                    inpl_corr = corrs(inpl_ind)             ! max in plane correlation
                    call self%store_solution(iref, iref, inpl_ind, inpl_corr)
                    projspace_corrs( iref ) = inpl_corr     ! stash in-plane correlation for sorting
                endif
            end subroutine per_ref_srch

    end subroutine stochastic_srch_snhc

    !> stochastic search het
    !! \param extr_bound corr threshold
    !! \param statecnt state counter array
    subroutine stochastic_srch_het( self, extr_bound, statecnt)
        use simple_rnd, only: shcloc, irnd_uni
        class(prime3D_srch), intent(inout) :: self
        real,                intent(in)    :: extr_bound
        integer,             intent(inout) :: statecnt(self%nstates)
        type(ori) :: o
        integer   :: iref, state
        real      :: corrs_inpl(pftcc_ptr%get_nrots())
        real      :: corr, mi_state, frac, corrs_state(self%nstates)
        self%prev_state = nint(a_ptr%get(self%iptcl, 'state'))
        if( self%prev_state > self%nstates ) stop 'previous best state outside boundary; stochastic_srch_het; simple_prime3D_srch'
        if( self%prev_state > 0 )then
            if( .not. state_exists_ptr(self%prev_state) ) stop 'empty previous state; stochastic_srch_het; simple_prime3D_srch'
            ! initialize
            o = a_ptr%get_ori(self%iptcl)
            self%prev_roind = pftcc_ptr%get_roind(360.-o%e3get())
            self%prev_proj  = e_ptr%find_closest_proj(o,1)
            if( a_ptr%get(self%iptcl,'corr') < extr_bound)then
                ! state randomization
                statecnt(self%prev_state) = statecnt(self%prev_state) + 1
                self%nrefs_eval = 1
                state = irnd_uni(self%nstates)
                do while(state == self%prev_state .or. .not.state_exists_ptr(state))
                    state = irnd_uni(self%nstates)
                enddo
                iref  = (state - 1) * self%nprojs + self%prev_proj
                corrs_inpl = pftcc_ptr%gencorrs(iref, self%iptcl)
                corr       = corrs_inpl(self%prev_roind)
            else
                ! SHC
                corrs_state = -1.
                do state = 1, self%nstates
                    if( .not.state_exists_ptr(state) )cycle
                    iref = (state-1) * self%nprojs + self%prev_proj
                    corrs_inpl = pftcc_ptr%gencorrs(iref, self%iptcl)
                    corrs_state(state) = corrs_inpl(self%prev_roind)
                enddo
                self%prev_corr = corrs_state(self%prev_state)
                state          = shcloc(self%nstates, corrs_state, self%prev_corr)
                if(state == 0)state = self%prev_state ! numerical stability
                corr            = corrs_state(state)
                self%nrefs_eval = count(corrs_state <= self%prev_corr)
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
            call o_peaks_ptr(self%iptcl)%set_ori(1, o)
            call a_ptr%set_ori(self%iptcl, o)
        else
            call a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED HET SEARCH'
    end subroutine stochastic_srch_het

    !>  \brief  executes the in-plane search. This improved routine took HCN
    !!          from stall @ 5.4 to close to convergence @ 4.5 after 10 iters
    !!          refine=no with 1000 referecnes
    subroutine inpl_srch( self )
        class(prime3D_srch), intent(inout) :: self
        real,    allocatable :: cxy(:), crxy(:)
        type(ori) :: o
        real      :: cc, e3
        integer   :: i, ref, irot
        if( self%doshift )then
            do i=self%nrefs,self%nrefs-self%npeaks+1,-1
                ref = proj_space_inds_ptr(self%iptcl, i )
                o   = o_refs_ptr(self%iptcl)%get_ori( ref )
                cc  = o_refs_ptr(self%iptcl)%get( ref, 'corr' )
                call self%shsrch_obj%set_indices(ref, self%iptcl)
                cxy = self%shsrch_obj%minimize(irot=irot)
                if( irot > 0 )then
                    e3 = 360. - pftcc_ptr%get_rot(irot) ! psi
                    call o%e3set(e3)                         ! stash psi
                    call o%set('corr', cxy(1))
                    call o%set_shift( cxy(2:3) )
                    call o_refs_ptr(self%iptcl)%set_ori( ref, o )
                endif
            end do
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED INPL SEARCH'
    end subroutine inpl_srch

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
        real      :: corrs(pftcc_ptr%get_nrots())
        type(ori) :: o_prev
        real      :: cc_t_min_1, corr
        if( str_has_substr(self%refine,'neigh') )then
            if( .not. present(nnmat) )&
            &stop 'need optional nnmat to be present for refine=neigh modes :: prep4srch (prime3D_srch)'
            if( .not. present(target_projs) )&
            &stop 'need optional target_projs to be present for refine=neigh modes :: prep4srch (prime3D_srch)'
        endif
        o_prev          = a_ptr%get_ori(self%iptcl)
        self%prev_state = nint(o_prev%get('state'))                                    ! state index
        self%prev_roind = pftcc_ptr%get_roind(360.-o_prev%e3get())                ! in-plane angle index
        self%prev_shvec = o_prev%get_2Dshift()                                         ! shift vector
        self%prev_proj  = e_ptr%find_closest_proj(o_prev,1)                       ! projection direction
        if( self%prev_state > self%nstates ) stop 'previous best state outside boundary; prep4srch; simple_prime3D_srch'
        if( self%prev_state > 0 )then
            if( .not. state_exists_ptr(self%prev_state) ) stop 'empty previous state; prep4srch; simple_prime3D_srch'
        endif
        select case( self%refine )
            case( 'no','shc','snhc','greedy')                               ! DISCRETE CASE
                self%prev_ref = o_refs_ptr(self%iptcl)%find_closest_proj(o_prev, self%prev_state) ! find closest ori with same state
            case( 'neigh','shcneigh', 'greedyneigh' )                                  ! DISCRETE CASE WITH NEIGHBOURHOOD
                nnvec = merge_into_disjoint_set(self%nprojs, self%nnn_static, nnmat, target_projs) ! disjoint nearest neighbour set
                self%prev_ref = o_refs_ptr(self%iptcl)%find_closest_proj(o_prev, self%prev_state) ! find closest ori with same state
            case( 'het' )
                self%prev_ref = (self%prev_state-1)*self%nprojs + self%prev_proj
            case( 'states' )
                self%prev_ref = o_refs_ptr(self%iptcl)%find_closest_proj(o_prev, self%prev_state)
            case DEFAULT
                stop 'Unknown refinement mode; simple_prime3D_srch; prep4srch'
        end select
        ! calc specscore
        self%specscore = pftcc_ptr%specscore(self%prev_ref, self%iptcl, self%prev_roind)
        ! prep corr
        if( self%refine .ne. 'het' )then
            corrs = pftcc_ptr%gencorrs(self%prev_ref, self%iptcl)
            corr  = max(0.,maxval(corrs))
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

    !>  \brief retrieves and preps npeaks orientations for reconstruction
    subroutine prep_npeaks_oris_and_weights( self )
        class(prime3D_srch),   intent(inout) :: self
        type(ori)  :: o
        type(oris) :: sym_os
        real       :: shvec(2), corrs(self%npeaks), ws(self%npeaks), logws(self%npeaks)
        real       :: frac, ang_sdev, dist, inpl_dist, euldist, mi_joint
        real       :: mi_proj, mi_inpl, mi_state, dist_inpl, wcorr
        integer    :: ipeak, cnt, ref, state, best_loc(1), order(self%npeaks), loc(1)
        integer    :: roind, neff_states
        logical    :: included(self%npeaks)
        ! empty states
        neff_states = 1
        if(self%nstates > 1)neff_states = count(state_exists_ptr)
        ! init npeaks
        do ipeak = 1, self%npeaks
            cnt = self%nrefs - self%npeaks + ipeak
            ref = proj_space_inds_ptr(self%iptcl, cnt )
            if( ref < 1 .or. ref > self%nrefs )then
                print *, 'ref: ', ref
                stop 'ref index out of bound; simple_prime3D_srch::prep_npeaks_oris'
            endif
            state = nint( o_refs_ptr(self%iptcl)%get(ref, 'state') )
            if( .not. state_exists_ptr(state) )then
                print *, 'empty state: ', state
                stop 'simple_prime3D_srch::prep_npeaks_oris'
            endif
            ! add shift
            shvec = self%prev_shvec
            if( self%doshift )shvec = shvec + o_refs_ptr(self%iptcl)%get_2Dshift(ipeak)
            where( abs(shvec) < 1e-6 ) shvec = 0.
            ! transfer to solution set
            corrs(ipeak) = o_refs_ptr(self%iptcl)%get(ref, 'corr')
            call o_peaks_ptr(self%iptcl)%set(ipeak, 'state', real(state))
            call o_peaks_ptr(self%iptcl)%set(ipeak, 'proj', o_refs_ptr(self%iptcl)%get(ref, 'proj'))
            call o_peaks_ptr(self%iptcl)%set(ipeak, 'corr', corrs(ipeak))
            call o_peaks_ptr(self%iptcl)%set_euler(ipeak,   o_refs_ptr(self%iptcl)%get_euler(ref))
            call o_peaks_ptr(self%iptcl)%set_shift(ipeak, shvec)
        enddo
        best_loc = maxloc(corrs)
        ! stochastic weights
        if( self%npeaks == 1 )then
            call o_peaks_ptr(self%iptcl)%set(1,'ow',1.0)
            wcorr = o_peaks_ptr(self%iptcl)%get(1,'corr')
        else
            ! calculate the exponential of the negative distances
            ! so that when diff==0 the weights are maximum and when
            ! diff==corrmax the weights are minimum
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
            call o_peaks_ptr(self%iptcl)%set_all('ow', ws)
        endif
        ! angular standard deviation
        ang_sdev = 0.
        if( trim(se_ptr%get_pgrp()).eq.'c1' )then
            ang_sdev = o_peaks_ptr(self%iptcl)%ang_sdev(self%refine, self%nstates, self%npeaks)
        else
            if( self%npeaks > 2 )then
                loc    = maxloc(corrs)
                sym_os = o_peaks_ptr(self%iptcl)
                do ipeak = 1, self%npeaks
                    if(ipeak == loc(1))cycle
                    o = o_peaks_ptr(self%iptcl)%get_ori(ipeak)
                    call se_ptr%sym_dists( o_peaks_ptr(self%iptcl)%get_ori(best_loc(1)), o, dist, inpl_dist)
                    call sym_os%set_ori(ipeak, o)
                enddo
                ang_sdev = sym_os%ang_sdev(self%refine, self%nstates, self%npeaks)
            endif
        endif
        ! Update the best orientation
        ! angular distances
        call se_ptr%sym_dists( a_ptr%get_ori(self%iptcl),&
            &o_peaks_ptr(self%iptcl)%get_ori(best_loc(1)), euldist, dist_inpl )
        ! convergence parameters
        state = nint( o_peaks_ptr(self%iptcl)%get(best_loc(1), 'state') )
        if( .not. state_exists_ptr(state) )then
            print *, 'empty state: ', state
            stop 'simple_prime3d_srch; update_best'
        endif
        roind = pftcc_ptr%get_roind( 360.-o_peaks_ptr(self%iptcl)%e3get(best_loc(1)) )
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
        ! fraction search space
        if( str_has_substr(self%refine, 'neigh') )then
            frac = 100.*real(self%nrefs_eval) / real(self%nnn * neff_states)
        else if( trim(self%refine).eq.'states' )then
            frac = 100.*real(self%nrefs_eval) / real(self%nnn) ! 1 state searched
        else
            frac = 100.*real(self%nrefs_eval) / real(self%nprojs * neff_states)
        endif
        ! set the overlaps
        call a_ptr%set(self%iptcl, 'mi_proj',  mi_proj )
        call a_ptr%set(self%iptcl, 'mi_inpl',  mi_inpl )
        call a_ptr%set(self%iptcl, 'mi_state', mi_state)
        call a_ptr%set(self%iptcl, 'mi_joint', mi_joint)
        ! set the distances before we update the orientation
        call a_ptr%set(self%iptcl, 'dist', 0.5*euldist + 0.5*a_ptr%get(self%iptcl,'dist'))
        call a_ptr%set(self%iptcl, 'dist_inpl', dist_inpl)
        ! all the other stuff
        call a_ptr%set_euler(self%iptcl, o_peaks_ptr(self%iptcl)%get_euler(best_loc(1)))
        call a_ptr%set_shift(self%iptcl, o_peaks_ptr(self%iptcl)%get_2Dshift(best_loc(1)))
        call a_ptr%set(self%iptcl, 'state', real(state))
        call a_ptr%set(self%iptcl, 'frac', frac )
        call a_ptr%set(self%iptcl, 'corr', wcorr )
        call a_ptr%set(self%iptcl, 'specscore', self%specscore)
        call a_ptr%set(self%iptcl, 'ow',    o_peaks_ptr(self%iptcl)%get(best_loc(1),'ow')   )
        call a_ptr%set(self%iptcl, 'proj',  o_peaks_ptr(self%iptcl)%get(best_loc(1),'proj') )
        call a_ptr%set(self%iptcl, 'sdev',  ang_sdev )
        DebugPrint '>>> PRIME3D_SRCH::EXECUTED PREP_NPEAKS_ORIS'
    end subroutine prep_npeaks_oris_and_weights

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
        proj_space_inds_ptr(self%iptcl,ind) = ref          ! stash reference similarly to other search modes
        e3 = 360. - pftcc_ptr%get_rot( inpl_ind ) ! psi
        call o_refs_ptr(self%iptcl)%set( ref, 'e3',   e3   )      ! stash psi
        call o_refs_ptr(self%iptcl)%set( ref, 'corr', corr )      ! stash correlation
        if( self%npeaks==1 ) call o_refs_ptr(self%iptcl)%set( ref, 'ow', 1. ) ! set reconstruction weight
    end subroutine store_solution

    subroutine nullify_prime3D_srch_ptrs
        pftcc_ptr           => null()
        a_ptr               => null()
        e_ptr               => null()
        se_ptr              => null()
        o_refs_ptr          => null()
        o_peaks_ptr         => null()
        srch_order_ptr      => null()
        proj_space_inds_ptr => null()
        state_exists_ptr    => null()
    end subroutine nullify_prime3D_srch_ptrs

end module simple_prime3D_srch
