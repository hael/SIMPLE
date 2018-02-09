! PRIME3D stochastic search routines
module simple_prime3D_srch
#include "simple_lib.f08"
use simple_oris,              only: oris
use simple_ori,               only: ori
use simple_sym,               only: sym
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_pftcc_shsrch,      only: pftcc_shsrch       ! simplex-based angle and shift search
use simple_pftcc_grad_shsrch, only: pftcc_grad_shsrch  ! gradient-based angle and shift search
implicit none

public :: cleanprime3D_srch, prep4prime3D_srch
public :: prime3D_srch, o_peaks
private

real,    parameter :: SOFTMAXW_THRESH = 0.01 !< threshold for softmax weights
logical, parameter :: DEBUG = .false.

! allocatables for prime3D_srch are class variables to improve caching and reduce alloc overheads
type(oris), allocatable :: o_peaks(:)                             !< solution objects
real,       allocatable :: proj_space_euls(:,:,:)                 !< euler angles
real,       allocatable :: proj_space_shift(:,:,:)                !< shift vector
real,       allocatable :: proj_space_corrs(:,:)                  !< correlations vs. reference orientations
integer,    allocatable :: proj_space_inds(:,:)                   !< stochastic index of reference orientations
integer,    allocatable :: proj_space_state(:,:)                  !< reference orientations state
integer,    allocatable :: proj_space_proj(:,:)                   !< reference orientations projection direction (1 state assumed)
integer,    allocatable :: prev_proj(:)                           !< particle previous reference projection direction
integer,    allocatable :: prev_states(:)                         !< particle previous state
integer,    allocatable :: srch_order(:,:)                        !< stochastic search index
logical,    allocatable :: state_exists(:)                        !< indicates state existence
logical,    pointer     :: ptcl_mask_ptr(:) => null()             !< pointer to particle mask

type prime3D_srch
    private
    class(polarft_corrcalc), pointer :: pftcc_ptr => null()       !< corrcalc object
    class(oris),             pointer :: a_ptr     => null()       !< b%a (primary particle orientation table)
    class(sym),              pointer :: se_ptr    => null()       !< b%se (symmetry elements)
    type(pftcc_shsrch)               :: shsrch_obj                !< origin shift search object
    type(pftcc_grad_shsrch)          :: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
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
    character(len=STDLEN)            :: refine        = ''        !< refinement flag
    character(len=STDLEN)            :: opt           = ''        !< optimizer flag
    logical                          :: doshift       = .true.    !< 2 indicate whether 2 serch shifts
    logical                          :: exists        = .false.   !< 2 indicate existence
  contains
    procedure          :: new
    procedure          :: exec_prime3D_srch
    procedure          :: exec_prime3D_srch_het
    procedure          :: calc_corr
    procedure, private :: greedy_srch
    procedure, private :: greedy_subspace_srch
    procedure, private :: stochastic_srch
    procedure, private :: stochastic_srch_shc
    procedure, private :: stochastic_srch_snhc
    procedure, private :: stochastic_srch_het
    procedure, private :: inpl_srch
    procedure, private :: prep4srch
    procedure, private :: prep_npeaks_oris_and_weights
    procedure, private :: store_solution
    procedure          :: kill
end type prime3D_srch

contains

    subroutine cleanprime3D_srch
        if( allocated(proj_space_euls)  ) deallocate(proj_space_euls)
        if( allocated(proj_space_shift) ) deallocate(proj_space_shift)
        if( allocated(proj_space_corrs) ) deallocate(proj_space_corrs)
        if( allocated(proj_space_inds)  ) deallocate(proj_space_inds)
        if( allocated(proj_space_state) ) deallocate(proj_space_state)
        if( allocated(proj_space_proj)  ) deallocate(proj_space_proj)
        if( allocated(prev_proj)        ) deallocate(prev_proj)
        if( allocated(srch_order)       ) deallocate(srch_order)
        if( allocated(state_exists)     ) deallocate(state_exists)
        if( allocated(prev_states)      ) deallocate(prev_states)
    end subroutine cleanprime3D_srch

    subroutine prep4prime3D_srch( b, p, ptcl_mask )
        use simple_ran_tabu, only: ran_tabu
        use simple_params,   only: params
        use simple_build,    only: build
        class(build),    intent(inout) :: b
        class(params),   intent(inout) :: p
        logical, target, intent(in)    :: ptcl_mask(p%fromp:p%top)
        integer,        allocatable :: pinds(:)
        type(ran_tabu), allocatable :: rts(:)
        type(ran_tabu) :: rt
        integer        :: i, istate, iproj, iptcl, prev_state
        integer        :: nnnrefs, cnt, prev_ref, nrefs, nptcls
        ! clean all class arrays & types
        call cleanprime3D_srch()
        if( allocated(o_peaks) )then
            do i=p%fromp,p%top
                call o_peaks(i)%kill
            enddo
            deallocate(o_peaks)
        endif
        ! parameters
        nrefs  = p%nspace * p%nstates
        nptcls = count(ptcl_mask)
        ! particle index mapping
        allocate(pinds(p%fromp:p%top))
        pinds = 0
        cnt   = 0
        do i=p%fromp,p%top
            if( ptcl_mask(i) )then
                cnt = cnt + 1
                pinds(i) = cnt
            endif
        end do
        ! set patcl_mask pointer
        ptcl_mask_ptr => ptcl_mask
        ! shared-memory arrays
        allocate(proj_space_euls(nptcls,nrefs,3), proj_space_shift(nptcls,nrefs,2),&
            &proj_space_corrs(nptcls,nrefs), proj_space_inds(nptcls,nrefs),&
            &proj_space_state(nptcls,nrefs), proj_space_proj(nptcls,nrefs),&
            &prev_proj(nptcls) )
        ! states existence
        if( p%oritab.ne.'' )then
            state_exists = b%a%states_exist(p%nstates)
        else
            allocate(state_exists(p%nstates), source=.true.)
        endif
        ! The shared memory used in a parallel section should be initialised
        ! with a (redundant) parallel section, because of how pages are organised.
        ! Memory otherwise becomes associated with the single thread used for
        ! allocation, causing load imbalance. This will reduce cache misses.
        !$omp parallel default(shared) private(i,iptcl,cnt,istate,iproj) proc_bind(close)
        !$omp do schedule(static)
        do i=1,nptcls
            proj_space_shift(i,:,:) = 0.
            proj_space_corrs(i,:)   = -1.
            proj_space_inds( i,:)   = 0
            proj_space_state(i,:)   = 0
            proj_space_proj( i,:)   = 0
            ! reference projection directions
            cnt = 0
            do istate=1,p%nstates
                do iproj=1,p%nspace
                    cnt = cnt + 1
                    proj_space_state(i,cnt)   = istate
                    proj_space_proj(i, cnt)   = iproj
                    proj_space_euls(i, cnt,:) = b%e%get_euler(iproj)
                enddo
            enddo
        end do
        !$omp end do nowait
        !$omp do schedule(static)
        do iptcl = p%fromp, p%top
            if( pinds(iptcl) > 0 ) prev_proj(pinds(iptcl)) = b%e%find_closest_proj(b%a%get_ori(iptcl))
        enddo
        !$omp end do nowait
        !$omp end parallel
        ! projection direction peaks, eo & CTF transfer
        select case(trim(p%refine))
            case('het','hetsym')
               ! nothing to do
            case DEFAULT
                allocate(o_peaks(p%fromp:p%top))
                do iptcl = p%fromp, p%top
                    if( ptcl_mask(iptcl) )then
                        ! ORIENTATION PEAKS
                        call o_peaks(iptcl)%new_clean(p%npeaks)
                        ! transfer CTF params
                        if( p%ctf.ne.'no' )then
                            call o_peaks(iptcl)%set_all2single('kv',    b%a%get(iptcl,'kv')   )
                            call o_peaks(iptcl)%set_all2single('cs',    b%a%get(iptcl,'cs')   )
                            call o_peaks(iptcl)%set_all2single('fraca', b%a%get(iptcl,'fraca'))
                            call o_peaks(iptcl)%set_all2single('dfx',   b%a%get(iptcl,'dfx')  )
                            if( p%tfplan%mode .eq. 'astig' )then
                                call o_peaks(iptcl)%set_all2single('dfy',    b%a%get(iptcl,'dfy')   )
                                call o_peaks(iptcl)%set_all2single('angast', b%a%get(iptcl,'angast'))
                            else
                                call o_peaks(iptcl)%set_all2single('dfy',    b%a%get(iptcl,'dfx'))
                                call o_peaks(iptcl)%set_all2single('angast', 0.)
                            endif
                            if( p%tfplan%l_phaseplate )then
                                call o_peaks(iptcl)%set_all2single('phshift', b%a%get(iptcl,'phshift'))
                            endif
                        endif
                        ! transfer eo flag
                        call o_peaks(iptcl)%set_all2single('eo', b%a%get(iptcl,'eo'))
                    endif
                enddo
        end select
        ! refine mode specific allocations and initialisations
        select case( trim(p%refine) )
            case( 'het', 'hetsym' )
                allocate(prev_states(nptcls), source=0)
                !$omp parallel do default(shared) private(i,iptcl) schedule(static) proc_bind(close)
                do iptcl = p%fromp, p%top
                    i = pinds(iptcl)
                    if( i > 0 )then
                        prev_states(i) = b%a%get_state(iptcl)
                        if( prev_states(i) > 0 )then
                             if(.not. state_exists(prev_states(i)) )&
                            &stop 'empty previous state; prep4prime3D_srch; simple_prime3D_srch'
                            if( prev_states(i) > p%nstates )&
                            &stop 'previous best state outside boundary; prep4prime3D_srch; simple_prime3D_srch'
                        endif
                    endif
                end do
                !$omp end parallel do
            case( 'neigh','shcneigh', 'greedyneigh' )
                nnnrefs =  p%nnn * p%nstates
                allocate(srch_order(nptcls, nnnrefs), rts(nptcls))
                do i=1,nptcls
                    rts(i) = ran_tabu(nnnrefs)
                end do
                !$omp parallel do default(shared) private(iptcl,prev_state,prev_ref,istate,i) schedule(static) proc_bind(close)
                do iptcl = p%fromp, p%top
                    if( pinds(iptcl) > 0 )then
                        prev_state = b%a%get_state(iptcl)
                        prev_ref   = (prev_state - 1) * p%nspace + prev_proj(pinds(iptcl))
                        do istate = 0, p%nstates - 1
                            i = istate * p%nnn + 1
                            srch_order(pinds(iptcl),i:i + p%nnn - 1) = b%nnmat(prev_proj(pinds(iptcl)),:) + istate * p%nspace
                        enddo
                        call rts(pinds(iptcl))%shuffle(srch_order(pinds(iptcl),:))
                        call put_last(prev_ref, srch_order(pinds(iptcl),:))
                    endif
                enddo
                !$omp end parallel do
                do i=1,nptcls
                    call rts(i)%kill
                end do
            case('no','shc','snhc','greedy')
                allocate(srch_order(nptcls,nrefs), rts(nptcls))
                do i=1,nptcls
                    rts(i) = ran_tabu(nrefs)
                end do
                !$omp parallel do default(shared) private(iptcl,prev_state,prev_ref,istate,i) schedule(static) proc_bind(close)
                do iptcl = p%fromp, p%top
                    if( pinds(iptcl) > 0 )then
                        call rts(pinds(iptcl))%ne_ran_iarr(srch_order(pinds(iptcl),:))
                        prev_state = b%a%get_state(iptcl)
                        prev_ref   = (prev_state - 1) * p%nspace + prev_proj(pinds(iptcl))
                        call put_last(prev_ref, srch_order(pinds(iptcl),:))
                    endif
                enddo
                !$omp end parallel do
                do i=1,nptcls
                    call rts(i)%kill
                end do
            case DEFAULT
                stop 'Unknown refinement mode; simple_hadamard3D_matcher; prep4primesrch3D'
        end select
        call rt%kill
        if( allocated(srch_order) )then
            if( any(srch_order == 0) ) stop 'Invalid index in srch_order; simple_hadamard3D_matcher :: prep4primesrch3D'
        endif
    end subroutine prep4prime3D_srch

    subroutine new( self, iptcl, iptcl_map, p, pftcc, a, se )
        use simple_params, only: params
        class(prime3D_srch),             intent(inout) :: self      !< instance
        integer,                         intent(in)    :: iptcl     !< global particle index
        integer,                         intent(in)    :: iptcl_map !< index in pre-allocated 2D arrays
        class(params),                   intent(in)    :: p         !< parameters
        class(polarft_corrcalc), target, intent(in)    :: pftcc     !< corrcalc object
        class(oris),             target, intent(in)    :: a         !< b%a (primary particle orientation table)
        class(sym),              target, intent(in)    :: se        !< b%se (symmetry elements)
        integer :: nstates_eff
        real    :: lims(2,2), lims_init(2,2)
        ! set constants
        self%pftcc_ptr  => pftcc
        self%a_ptr      => a
        self%se_ptr     => se
        self%iptcl      =  iptcl
        self%iptcl_map  =  iptcl_map
        self%nstates    =  p%nstates
        self%nprojs     =  p%nspace
        self%nrefs      =  self%nprojs*self%nstates
        self%nrots      =  round2even(twopi*real(p%ring2))
        self%npeaks     =  p%npeaks
        self%nbetter    =  0
        self%nrefs_eval =  0
        self%nsym       =  se%get_nsym()
        self%doshift    =  p%l_doshift
        self%refine     =  p%refine
        self%nnn_static =  p%nnn
        self%nnn        =  p%nnn
        self%nnnrefs    =  self%nnn*self%nstates
        self%kstop_grid =  p%kstop_grid
        self%opt        =  p%opt
        if( str_has_substr(self%refine,'shc') )then
            if( self%npeaks > 1 ) stop 'npeaks must be equal to 1 with refine=shc|shcneigh'
        endif
        ! multiple states
        if( self%nstates == 1 )then
            self%npeaks_grid = GRIDNPEAKS
        else
            ! number of populated states
            nstates_eff = count(state_exists)
            select case(trim(p%refine))
            case('het','hetsym')
                    self%npeaks_grid = 1
                case DEFAULT
                    ! "-(nstates_eff-1)" because all states share the same previous orientation
                    self%npeaks_grid = GRIDNPEAKS * nstates_eff - (nstates_eff - 1)
            end select
        endif
        if(str_has_substr(self%refine,'neigh'))then
            self%npeaks_grid = min(self%npeaks_grid,self%nnnrefs)
        else
            self%npeaks_grid = min(self%npeaks_grid,self%nrefs)
        endif
        ! updates option to search shift
        self%doshift   = p%l_doshift
        ! create in-plane search objects
        lims(:,1)      = -p%trs
        lims(:,2)      =  p%trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        call self%shsrch_obj%new(self%pftcc_ptr, lims, lims_init=lims_init,&
            &shbarrier=p%shbarrier, nrestarts=3, maxits=60)
        call self%grad_shsrch_obj%new(self%pftcc_ptr, lims, lims_init=lims_init,&
            &shbarrier=p%shbarrier, maxits=60)
        self%exists = .true.
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::CONSTRUCTED NEW SIMPLE_PRIME3D_SRCH OBJECT'
    end subroutine new

    ! SEARCH ROUTINES

    !>  \brief  exec_prime3D_srch is a master prime search routine
    !! \param lp low-pass cutoff freq
    !! \param greedy greedy search flag
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs projection indices for grid search
    !! \param szsn optional size for snhc refinement
    subroutine exec_prime3D_srch( self, greedy, nnmat, grid_projs, szsn )
        class(prime3D_srch), intent(inout) :: self
        logical, optional,   intent(in)    :: greedy
        integer, optional,   intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:), szsn
        logical :: ggreedy
        ggreedy = .false.
        if( present(greedy) ) ggreedy = greedy
        if( ggreedy )then
            call self%greedy_srch(nnmat, grid_projs)
        else if( self%refine.eq.'snhc' )then
            if( .not. present(szsn) )then
                stop 'refine=snhc mode needs optional input szsn; simple_prime3D_srch :: exec_prime3D_srch'
            endif
            call self%stochastic_srch_snhc(szsn)
        else if( str_has_substr(self%refine,'shc') )then
            call self%stochastic_srch_shc(nnmat, grid_projs)
        else
            call self%stochastic_srch(nnmat, grid_projs)
        endif
        if(allocated(self%nnvec))deallocate(self%nnvec)
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SRCH'
    end subroutine exec_prime3D_srch

    !>  \brief state labeler
    subroutine exec_prime3D_srch_het( self, corr_thresh, do_extr_opt, counts, symmat )
        class(prime3D_srch),        intent(inout) :: self
        real,                       intent(in)    :: corr_thresh
        logical,                    intent(in)    :: do_extr_opt
        integer,                    intent(inout) :: counts(2)
        integer,          optional, intent(in)    :: symmat(self%nprojs, self%nsym)
        call self%stochastic_srch_het( corr_thresh, do_extr_opt, counts, symmat )
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SRCH_HET'
    end subroutine exec_prime3D_srch_het

    !>  \brief  greedy hill-climbing
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs projection indices for grid search
    subroutine greedy_srch( self, nnmat, grid_projs )
        class(prime3D_srch), intent(inout) :: self
        integer, optional,   intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:)
        integer :: iref,isample,nrefs,target_projs(self%npeaks_grid), inpl_ind(1), state
        real    :: corrs(self%nrefs), inpl_corrs(self%nrots), inpl_corr
        if( self%a_ptr%get_state(self%iptcl) > 0 )then
            if( str_has_substr(self%refine, 'neigh') )then
                ! for neighbour modes we do a coarse grid search first
                if( .not. present(grid_projs) ) stop 'need optional grid_projs 4 subspace srch; prime3D_srch :: greedy_srch'
                call self%greedy_subspace_srch(grid_projs, target_projs)
                ! initialize
                call self%prep4srch(nnmat, target_projs)
                nrefs = self%nnnrefs
            else
                ! initialize
                call self%prep4srch()
                nrefs = self%nrefs
            endif
            self%nbetter    = 0
            self%nrefs_eval = 0
            proj_space_corrs(self%iptcl_map,:) = -1.
            ! search
            do isample=1,nrefs
                iref = srch_order(self%iptcl_map,isample) ! set the reference index
                call per_ref_srch                         ! actual search
            end do
            ! in greedy mode, we evaluate all refs
            self%nrefs_eval = nrefs
            ! sort in correlation projection direction space
            corrs = proj_space_corrs(self%iptcl_map,:)
            call hpsort(corrs, proj_space_inds(self%iptcl_map,:))
            ! take care of the in-planes
            call self%inpl_srch
            ! prepare weights & orientation
            call self%prep_npeaks_oris_and_weights
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::FINISHED GREEDY SEARCH'

        contains

            subroutine per_ref_srch
                state = proj_space_state(self%iptcl_map,iref)
                if( state_exists(state) )then
                    call self%pftcc_ptr%gencorrs(iref, self%iptcl, inpl_corrs) ! In-plane correlations
                    inpl_ind   = maxloc(inpl_corrs)                            ! greedy in-plane index
                    inpl_corr  = inpl_corrs(inpl_ind(1))                       ! max in plane correlation
                    call self%store_solution(iref, iref, inpl_ind(1), inpl_corr)
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
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::FINISHED GREEDY SUBSPACE SEARCH'

        contains

            subroutine per_ref_srch
                if( state_exists(istate) )then
                    call self%pftcc_ptr%gencorrs(iref, self%iptcl, self%kstop_grid, inpl_corrs) ! In-plane correlations
                    inpl_ind   = maxloc(inpl_corrs)               ! greedy in-plane
                    inpl_corr  = inpl_corrs(inpl_ind(1))          ! max in plane correlation
                    proj_space_corrs(self%iptcl_map,iref) = inpl_corr ! stash in-plane correlation for sorting
                    proj_space_inds(self%iptcl_map,iref)  = iref      ! stash the index for sorting
                endif
            end subroutine per_ref_srch

    end subroutine greedy_subspace_srch

    !>  \brief  executes the stochastic soft orientation search
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs projection indices for grid search
    subroutine stochastic_srch( self, nnmat, grid_projs )
        class(prime3D_srch), intent(inout) :: self
        integer, optional,   intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:)
        integer :: iref,isample,nrefs,target_projs(self%npeaks_grid), loc(1), inpl_ind
        real    :: corrs(self%nrefs), inpl_corrs(self%nrots), inpl_corr
        ! execute search
        if( self%a_ptr%get_state(self%iptcl) > 0 )then
            if( str_has_substr(self%refine, 'neigh') )then
                ! for neighbour modes we do a coarse grid search first
                if( .not. present(grid_projs) ) stop 'need optional grid_projs 4 subspace srch; prime3D_srch :: stochastic_srch'
                call self%greedy_subspace_srch(grid_projs, target_projs)
                ! initialize
                call self%prep4srch(nnmat, target_projs)
                nrefs = self%nnnrefs
            else
                ! initialize
                call self%prep4srch()
                nrefs = self%nrefs
            endif
            ! initialize, ctd
            self%nbetter    = 0
            self%nrefs_eval = 0
            proj_space_corrs(self%iptcl_map,:) = -1.
            do isample=1,nrefs
                iref = srch_order(self%iptcl_map,isample)  ! set the stochastic reference index
                call per_ref_srch                          ! actual search
                if( self%nbetter >= self%npeaks ) exit     ! exit condition
            end do
            ! sort in correlation projection direction space
            corrs = proj_space_corrs(self%iptcl_map,:)
            call hpsort(corrs, proj_space_inds(self%iptcl_map,:))
            call self%inpl_srch ! search shifts
            ! prepare weights and orientations
            call self%prep_npeaks_oris_and_weights
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH'

        contains

            subroutine per_ref_srch
                if( state_exists( proj_space_state(self%iptcl_map,iref) ) )then
                    ! In-plane correlations
                    call self%pftcc_ptr%gencorrs(iref, self%iptcl, inpl_corrs)
                    loc       = maxloc(inpl_corrs)   ! greedy in-plane
                    inpl_ind  = loc(1)               ! in-plane angle index
                    inpl_corr = inpl_corrs(inpl_ind) ! max in plane correlation
                    call self%store_solution(iref, iref, inpl_ind, inpl_corr)
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
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs projection indices for grid search
    subroutine stochastic_srch_shc( self, nnmat, grid_projs )
        class(prime3D_srch), intent(inout) :: self
        integer, optional,   intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:)
        real    :: inpl_corr, inpl_corrs(self%nrots)
        integer :: iref,isample,nrefs,inpl_ind,loc(1),target_projs(self%npeaks_grid)
        logical :: found_better
        if( self%a_ptr%get_state(self%iptcl) > 0 )then
            if( str_has_substr(self%refine, 'neigh') )then
                ! for neighbour modes we do a coarse grid search first
                if( .not. present(grid_projs) ) stop 'need optional grid_projs 4 subspace srch; prime3D_srch :: stochastic_srch_shc'
                call self%greedy_subspace_srch(grid_projs, target_projs)
                ! initialize
                call self%prep4srch(nnmat, target_projs)
                nrefs = self%nnnrefs
            else
                ! initialize
                call self%prep4srch()
                nrefs = self%nrefs
            endif
            ! initialize
            found_better    = .false.
            self%nrefs_eval = 0
            ! search
            do isample=1,nrefs
                iref = srch_order(self%iptcl_map, isample ) ! stochastic reference index
                call per_ref_srch                           ! actual search
                if( inpl_ind > 0 )then
                    ! correlation-improving reference found, store it
                    call self%store_solution(self%nrefs, iref, inpl_ind, inpl_corr)
                    found_better = .true. ! flag for found solution
                    exit
                endif
            end do
            if( .not. found_better )then
                iref = self%prev_ref
                call per_ref_srch
                if( inpl_ind == 0 )then
                    ! be greedy
                    loc       = maxloc(inpl_corrs)
                    inpl_ind  = loc(1)
                    inpl_corr = inpl_corrs(inpl_ind)
                endif
                call self%store_solution(self%nrefs, iref, inpl_ind, inpl_corr)
            endif
            ! search shifts
            call self%inpl_srch
            ! output
            call self%prep_npeaks_oris_and_weights
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH (REFINE=SHC|SHCNEIGH)'

        contains

            subroutine per_ref_srch
                use simple_rnd, only: shcloc
                inpl_ind  = 0                                                    ! init in-plane index
                inpl_corr = 0.                                                   ! init correlation
                if( state_exists( proj_space_state(self%iptcl_map,iref) ) )then
                    call self%pftcc_ptr%gencorrs( iref, self%iptcl, inpl_corrs ) ! in-plane correlations
                    inpl_ind = shcloc(self%nrots, inpl_corrs, self%prev_corr)    ! first improving in-plane index
                    if( inpl_ind > 0 ) inpl_corr = inpl_corrs( inpl_ind )        ! improving correlation found
                    self%nrefs_eval = self%nrefs_eval + 1                        ! updates fractional search space
                endif
            end subroutine per_ref_srch

    end subroutine stochastic_srch_shc

    !>  \brief  stochastic neighborhood hill-climbing
    !! \param szsn size ref evals
    subroutine stochastic_srch_snhc( self, szsn )
        class(prime3D_srch), intent(inout) :: self
        integer,             intent(in)    :: szsn
        integer :: iref, isample, loc(1), inpl_ind
        real    :: corrs(self%nrefs), inpl_corrs(self%nrots), inpl_corr
        if( self%a_ptr%get_state(self%iptcl) > 0 )then
            ! initialize
            call self%prep4srch()
            self%nbetter    = 0
            self%nrefs_eval = 0
            proj_space_corrs(self%iptcl_map,:) = -1.
            ! search
            do isample=1,szsn
                iref = srch_order(self%iptcl_map,isample) ! set the stochastic reference index
                call per_ref_srch                         ! actual search
            end do
            self%nrefs_eval = szsn
            ! sort in correlation projection direction space
            corrs = proj_space_corrs(self%iptcl_map,:)
            call hpsort(corrs, proj_space_inds(self%iptcl_map,:))
            ! output
            call self%prep_npeaks_oris_and_weights
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::FINISHED EXTREMAL SEARCH'

        contains

            subroutine per_ref_srch
                if( state_exists( proj_space_state(self%iptcl_map,iref) ) )then
                    call self%pftcc_ptr%gencorrs(iref, self%iptcl, inpl_corrs) ! In-plane correlations
                    loc        = maxloc(inpl_corrs)               ! greedy in-plane
                    inpl_ind   = loc(1)                           ! in-plane angle index
                    inpl_corr  = inpl_corrs(inpl_ind)             ! max in plane correlation
                    call self%store_solution(iref, iref, inpl_ind, inpl_corr)
                endif
            end subroutine per_ref_srch

    end subroutine stochastic_srch_snhc

    !> stochastic search heterogeneity
    subroutine stochastic_srch_het( self, corr_thresh, do_extr_opt, counts, symmat )
        use simple_rnd,      only: shcloc, irnd_uni
        class(prime3D_srch),   intent(inout) :: self
        real,                  intent(in)    :: corr_thresh
        logical,               intent(in)    :: do_extr_opt
        integer,               intent(inout) :: counts(2)
        integer,     optional, intent(in)    :: symmat(self%nprojs, self%nsym)
        integer :: sym_projs(self%nstates), loc(1), iproj, iref, isym, state
        real    :: corrs(self%nstates), corrs_sym(self%nsym), corrs_inpl(self%nrots)
        real    :: shvec(2), corr, mi_state, frac, mi_inpl, mi_proj
        logical :: hetsym
        if( prev_states(self%iptcl_map) > 0 )then
            hetsym = present(symmat)
            self%prev_roind = self%pftcc_ptr%get_roind(360.-self%a_ptr%e3get(self%iptcl))
            self%prev_corr  = self%a_ptr%get(self%iptcl, 'corr')
            ! evaluate all correlations
            corrs = -1.
            do state = 1, self%nstates
                if( .not.state_exists(state) )cycle
                if( hetsym )then
                    ! greedy in symmetric unit
                    do isym = 1, self%nsym
                        iproj = symmat(prev_proj(self%iptcl_map), isym)
                        iref  = (state-1) * self%nprojs + iproj
                        call self%pftcc_ptr%gencorrs(iref, self%iptcl, corrs_inpl)
                        corrs_sym(isym) = corrs_inpl(self%prev_roind)
                    enddo
                    loc              = maxloc(corrs_sym)
                    isym             = loc(1)
                    corrs(state)     = corrs_sym(isym)
                    sym_projs(state) = symmat(prev_proj(self%iptcl_map), isym)
                else
                    iref = (state-1) * self%nprojs + prev_proj(self%iptcl_map)
                    call self%pftcc_ptr%gencorrs(iref, self%iptcl, corrs_inpl)
                    ! corrs(state) = corrs_inpl(self%prev_roind)
                    ! this is to somewhat take into account alignment errors
                    corrs(state) = maxval(corrs_inpl)
                endif
            enddo
            ! make moves
            mi_state = 0.
            mi_inpl  = 1.
            mi_proj  = 1.
            if( self%prev_corr < corr_thresh )then
                ! state randomization
                counts(1) = counts(1) + 1
                state = irnd_uni(self%nstates)
                do while(state == prev_states(self%iptcl_map) .or. .not.state_exists(state))
                    state = irnd_uni(self%nstates)
                enddo
                corr            = corrs(state)
                self%nrefs_eval = 1
            else
                ! SHC state optimization
                counts(2) = counts(2) + 1
                self%prev_corr  = corrs(prev_states(self%iptcl_map))
                state           = shcloc(self%nstates, corrs, self%prev_corr)
                corr            = corrs(state)
                self%nrefs_eval = count(corrs <= self%prev_corr)
                if( prev_states(self%iptcl_map) .eq. state ) mi_state = 1.
                if( do_extr_opt )then
                    ! extremal optimization
                    if( hetsym )then
                        ! greedy in symmetric unit
                        if( sym_projs(state) .ne. prev_proj(self%iptcl_map) )then
                            mi_proj = 0.
                            iref    = (state-1)*self%nprojs + sym_projs(state)
                            call self%a_ptr%e1set(self%iptcl, proj_space_euls(self%iptcl_map, iref, 1))
                            call self%a_ptr%e2set(self%iptcl, proj_space_euls(self%iptcl_map, iref, 2))
                        endif
                    endif
                else
                    ! greedy moves after extremal optimization complete
                    if( mi_state > 0.5 )then
                        if( hetsym )then
                            ! greedy in symmetric units
                            iproj   = sym_projs(state)
                            if( iproj .ne. prev_proj(self%iptcl_map) )then
                                mi_proj = 0.
                                iref    = (state-1) * self%nprojs + iproj
                                call self%a_ptr%e1set(self%iptcl, proj_space_euls(self%iptcl_map, iref, 1))
                                call self%a_ptr%e2set(self%iptcl, proj_space_euls(self%iptcl_map, iref, 2))
                            endif
                        else
                            ! greedy in plane
                            iref = (state-1)*self%nprojs + prev_proj(self%iptcl_map)
                            proj_space_corrs(self%iptcl_map,iref)      = corr
                            proj_space_inds(self%iptcl_map,self%nrefs) = iref
                            call self%inpl_srch
                            if( proj_space_corrs(self%iptcl_map,iref) > corr )then
                                corr  = proj_space_corrs(self%iptcl_map,iref)
                                shvec = self%a_ptr%get_2Dshift(self%iptcl) + proj_space_shift(self%iptcl_map,iref,:)
                                call self%a_ptr%set_shift(self%iptcl, shvec)
                                call self%a_ptr%e3set(self%iptcl, proj_space_euls(self%iptcl_map, iref, 3))
                            endif
                            if( self%prev_roind .ne. self%pftcc_ptr%get_roind(360.-proj_space_euls(self%iptcl_map, iref, 3)) )then
                                mi_inpl = 0.
                            endif
                        endif
                    endif
                endif
            endif
            ! updates peaks and orientation orientation
            frac = 100.*real(self%nrefs_eval) / real(self%nstates)
            call self%a_ptr%set(self%iptcl,'frac',     frac)
            call self%a_ptr%set(self%iptcl,'state',    real(state))
            call self%a_ptr%set(self%iptcl,'corr',     corr)
            call self%a_ptr%set(self%iptcl,'mi_proj',  mi_proj)
            call self%a_ptr%set(self%iptcl,'mi_inpl',  mi_inpl)
            call self%a_ptr%set(self%iptcl,'mi_state', mi_state)
            call self%a_ptr%set(self%iptcl,'mi_joint', (mi_state+mi_inpl)/2.)
            call self%a_ptr%set(self%iptcl,'w',        1.)
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::FINISHED HET SEARCH'
    end subroutine stochastic_srch_het

    !> evaluates correlation
    subroutine calc_corr( self )
        class(prime3D_srch),   intent(inout) :: self
        integer :: iref
        real    :: corrs_inpl(self%nrots)
        if( prev_states(self%iptcl_map) > 0 )then
            self%prev_roind = self%pftcc_ptr%get_roind(360.-self%a_ptr%e3get(self%iptcl))
            iref = (prev_states(self%iptcl_map)-1) * self%nprojs + prev_proj(self%iptcl_map)
            call self%pftcc_ptr%gencorrs(iref, self%iptcl, corrs_inpl)
            call self%a_ptr%set(self%iptcl,'corr', corrs_inpl(self%prev_roind))
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::FINISHED calc_corr SEARCH'
    end subroutine calc_corr

    !>  \brief  executes the in-plane search. This improved routine took HCN
    !!          from stall @ 5.4 to close to convergence @ 4.5 after 10 iters
    !!          refine=no with 1000 references
    subroutine inpl_srch( self )
        class(prime3D_srch), intent(inout) :: self
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
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::FINISHED INPL SEARCH'
    end subroutine inpl_srch

    !>  \brief  prepares reference indices for the search & fetches ctf
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param target_projs projection indices for grid search
    subroutine prep4srch( self, nnmat, target_projs )
        use simple_combinatorics, only: merge_into_disjoint_set
        class(prime3D_srch), intent(inout) :: self
        integer, optional,   intent(in)    :: nnmat(self%nprojs,self%nnn_static), target_projs(self%npeaks_grid)
        real      :: corrs(self%nrots)
        type(ori) :: o_prev
        real      :: corr
        if( str_has_substr(self%refine,'neigh') )then
            if( .not. present(nnmat) )&
            &stop 'need optional nnmat to be present for refine=neigh modes :: prep4srch (prime3D_srch)'
            if( .not. present(target_projs) )&
            &stop 'need optional target_projs to be present for refine=neigh modes :: prep4srch (prime3D_srch)'
        endif
        o_prev          = self%a_ptr%get_ori(self%iptcl)
        self%prev_state = o_prev%get_state()                                          ! state index
        self%prev_roind = self%pftcc_ptr%get_roind(360.-o_prev%e3get())               ! in-plane angle index
        self%prev_shvec = o_prev%get_2Dshift()                                        ! shift vector
        self%prev_ref   = (self%prev_state-1)*self%nprojs + prev_proj(self%iptcl_map) ! previous reference
        if( self%prev_state > 0 )then
            if( self%prev_state > self%nstates ) stop 'previous best state outside boundary; prep4srch; simple_prime3D_srch'
            if( .not. state_exists(self%prev_state) ) stop 'empty previous state; prep4srch; simple_prime3D_srch'
        endif
        select case( self%refine )
            case( 'neigh','shcneigh', 'greedyneigh' )
                ! disjoint nearest neighbour set
                self%nnvec = merge_into_disjoint_set(self%nprojs, self%nnn_static, nnmat, target_projs)
            case( 'no','shc','snhc','greedy','yes','het' )
                ! all good
            case DEFAULT
                stop 'Unknown refinement mode; simple_prime3D_srch; prep4srch'
        end select
        ! calc specscore
        self%specscore = self%pftcc_ptr%specscore(self%prev_ref, self%iptcl, self%prev_roind)
        ! prep corr
        if( self%refine .ne. 'het' )then
            call self%pftcc_ptr%gencorrs(self%prev_ref, self%iptcl, corrs)
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
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::PREPARED FOR SIMPLE_PRIME3D_SRCH'
    end subroutine prep4srch

    !>  \brief retrieves and preps npeaks orientations for reconstruction
    subroutine prep_npeaks_oris_and_weights( self )
        class(prime3D_srch),   intent(inout) :: self
        type(ori)  :: osym
        type(oris) :: sym_os
        real       :: shvec(2), corrs(self%npeaks), ws(self%npeaks), dists(self%npeaks)
        real       :: arg4softmax(self%npeaks), state_ws(self%nstates)
        real       :: mi_proj, mi_inpl, mi_state, dist_inpl, wcorr, realfrac, frac
        real       :: ang_sdev, dist, inpl_dist, euldist, mi_joint, bfac
        integer    :: best_loc(1), loc(1), kfromto(2), states(self%npeaks)
        integer    :: s, ipeak, cnt, ref, state, roind, neff_states
        logical    :: included(self%npeaks)
        ! empty states
        neff_states = count(state_exists)
        ! init npeaks
        do ipeak = 1, self%npeaks
            cnt = self%nrefs - self%npeaks + ipeak
            ref = proj_space_inds(self%iptcl_map, cnt)
            if( ref < 1 .or. ref > self%nrefs )then
                print *, 'ref: ', ref
                stop 'ref index out of bound; simple_prime3D_srch::prep_npeaks_oris'
            endif
            state = proj_space_state(self%iptcl_map,ref)
            if( .not. state_exists(state) )then
                print *, 'empty state: ', state
                stop 'simple_prime3D_srch::prep_npeaks_oris'
            endif
            ! add shift
            shvec = self%prev_shvec
            if( self%doshift )shvec = shvec + proj_space_shift(self%iptcl_map,ref,1:2)
            where( abs(shvec) < 1e-6 ) shvec = 0.
            ! transfer to solution set
            corrs(ipeak) = proj_space_corrs(self%iptcl_map,ref)
            if( corrs(ipeak) < 0. ) corrs(ipeak) = 0.
            call o_peaks(self%iptcl)%set(ipeak, 'state', real(state))
            call o_peaks(self%iptcl)%set(ipeak, 'proj',  real(proj_space_proj(self%iptcl_map,ref)))
            call o_peaks(self%iptcl)%set(ipeak, 'corr',  corrs(ipeak))
            call o_peaks(self%iptcl)%set_euler(ipeak, proj_space_euls(self%iptcl_map,ref,1:3))
            call o_peaks(self%iptcl)%set_shift(ipeak, shvec)
        enddo
        best_loc = maxloc(corrs)
        ! stochastic weights
        if( self%npeaks == 1 )then
            call o_peaks(self%iptcl)%set(1,'ow',1.0)
            wcorr = o_peaks(self%iptcl)%get(1,'corr')
        else
            ! convert correlations to distances
            dists = 1.0 - corrs
            ! scale distances with TAU
            dists = dists / TAU
            ! argument for softmax function is negative distances
            arg4softmax = -dists
            ! subtract maxval of negative distances for numerical stability
            arg4softmax = arg4softmax - maxval(arg4softmax)
            ! calculate softmax weights
            ws = exp(arg4softmax)
            ws = ws / sum(ws)
            ! threshold weights
            included = (ws >= SOFTMAXW_THRESH)
            self%npeaks_eff = count(included)
            where( .not. included ) ws = 0.
            ! weighted corr
            wcorr = sum(ws*corrs,mask=included)
            ! update npeaks individual weights
            call o_peaks(self%iptcl)%set_all('ow', ws)
        endif
        if( self%npeaks > 1 .and. self%nstates > 1 )then
            ! states weights
            do ipeak = 1, self%npeaks
                corrs(ipeak)  = o_peaks(self%iptcl)%get(ipeak,'corr')
                states(ipeak) = nint(o_peaks(self%iptcl)%get(ipeak,'state'))
            enddo
            ! greedy state assignment
            do s = 1, self%nstates
                state_ws(s) = sum(ws,mask=(states==s))
            enddo
            loc      = maxloc(state_ws)
            state    = loc(1)
            ! in-state re-weighing
            included = included .and. (states==state)
            where( .not.included )ws = 0.
            ws       = ws / sum(ws, mask=included)
            best_loc = maxloc(ws)
            ! weighted corr
            wcorr    = sum(ws*corrs, mask=included)
            ! update npeaks individual weights
            call o_peaks(self%iptcl)%set_all('ow', ws)
        endif
        ! B factors
        do ipeak = 1, self%npeaks
            if( ws(ipeak) > TINY .or. self%npeaks==1 )then
                cnt   = self%nrefs - self%npeaks + ipeak
                ref   = proj_space_inds(self%iptcl_map, cnt)
                shvec = 0.
                if( self%doshift )shvec = proj_space_shift(self%iptcl_map, ref, 1:2)
                roind = self%pftcc_ptr%get_roind(360. - proj_space_euls(self%iptcl_map, ref, 3))
                bfac  = self%pftcc_ptr%calc_bfac(ref, self%iptcl, roind, shvec)
                call o_peaks(self%iptcl)%set(ipeak, 'bfac', bfac)
            else
                call o_peaks(self%iptcl)%set(ipeak, 'bfac', 0.)
            endif
        enddo
        ! angular standard deviation
        ang_sdev = 0.
        if( trim(self%se_ptr%get_pgrp()).eq.'c1' )then
            ang_sdev = o_peaks(self%iptcl)%ang_sdev(self%refine, self%nstates, self%npeaks)
        else
            if( self%npeaks > 2 )then
                sym_os = o_peaks(self%iptcl)
                do ipeak = 1, self%npeaks
                    if( ipeak == best_loc(1) )cycle
                    call self%se_ptr%sym_dists( o_peaks(self%iptcl)%get_ori(best_loc(1)),&
                        &o_peaks(self%iptcl)%get_ori(ipeak), osym, dist, inpl_dist)
                    call sym_os%set_ori(ipeak, osym)
                enddo
                ang_sdev = sym_os%ang_sdev(self%refine, self%nstates, self%npeaks)
            endif
        endif
        ! angular distances
        call self%se_ptr%sym_dists( self%a_ptr%get_ori(self%iptcl),&
            &o_peaks(self%iptcl)%get_ori(best_loc(1)), osym, euldist, dist_inpl)
        ! convergence parameters
        roind    = self%pftcc_ptr%get_roind(360. - o_peaks(self%iptcl)%e3get(best_loc(1)))
        mi_proj  = 0.
        mi_inpl  = 0.
        mi_joint = 0.
        if( euldist < 0.5 )then
            mi_proj = 1.
            mi_joint = mi_joint + 1.
        endif
        if( self%prev_roind == roind )then
            mi_inpl  = 1.
            mi_joint = mi_joint + 1.
        endif
        ! states convergence
        mi_state = 0.
        state = nint( o_peaks(self%iptcl)%get(best_loc(1), 'state') )
        if( .not. state_exists(state) )then
            print *, 'empty state: ', state
            stop 'simple_prime3d_srch; update_best'
        endif
        if(self%nstates > 1)then
            if( self%prev_state == state )then
                mi_state = 1.
                mi_joint = mi_joint + mi_state
            endif
            mi_joint = mi_joint/3.
        else
            mi_joint = mi_joint/2.
        endif
        ! fraction search space
        if( str_has_substr(self%refine, 'neigh') )then
            frac = 100.*real(self%nrefs_eval) / real(self%nnn * neff_states)
        else
            frac = 100.*real(self%nrefs_eval) / real(self%nprojs * neff_states)
        endif
        ! set the overlaps
        call self%a_ptr%set(self%iptcl, 'mi_proj',  mi_proj )
        call self%a_ptr%set(self%iptcl, 'mi_inpl',  mi_inpl )
        call self%a_ptr%set(self%iptcl, 'mi_state', mi_state)
        call self%a_ptr%set(self%iptcl, 'mi_joint', mi_joint)
        ! set the distances before we update the orientation
        call self%a_ptr%set(self%iptcl, 'dist', 0.5*euldist + 0.5*self%a_ptr%get(self%iptcl,'dist'))
        call self%a_ptr%set(self%iptcl, 'dist_inpl', dist_inpl)
        ! all the other stuff
        call self%a_ptr%set_euler(self%iptcl, o_peaks(self%iptcl)%get_euler(best_loc(1)))
        call self%a_ptr%set_shift(self%iptcl, o_peaks(self%iptcl)%get_2Dshift(best_loc(1)))
        call self%a_ptr%set(self%iptcl, 'state', real(state))
        call self%a_ptr%set(self%iptcl, 'frac', frac )
        call self%a_ptr%set(self%iptcl, 'corr', wcorr )
        call self%a_ptr%set(self%iptcl, 'specscore', self%specscore)
        call self%a_ptr%set(self%iptcl, 'ow',    o_peaks(self%iptcl)%get(best_loc(1),'ow')   )
        call self%a_ptr%set(self%iptcl, 'proj',  o_peaks(self%iptcl)%get(best_loc(1),'proj') )
        call self%a_ptr%set(self%iptcl, 'sdev',  ang_sdev )
        call self%a_ptr%set(self%iptcl, 'npeaks', real(self%npeaks_eff) )
        if( DEBUG ) print *,  '>>> PRIME3D_SRCH::EXECUTED PREP_NPEAKS_ORIS'
    end subroutine prep_npeaks_oris_and_weights

    subroutine store_solution( self, ind, ref, inpl_ind, corr )
        class(prime3D_srch),     intent(inout) :: self
        integer,                 intent(in)    :: ind, ref, inpl_ind
        real,                    intent(in)    :: corr
        proj_space_inds(self%iptcl_map,ind)   = ref
        proj_space_euls(self%iptcl_map,ref,3) = 360. - self%pftcc_ptr%get_rot(inpl_ind)
        proj_space_corrs(self%iptcl_map,ref)  = corr
    end subroutine store_solution

    subroutine kill( self )
        class(prime3D_srch),  intent(inout) :: self !< instance
        if(allocated(self%nnvec))deallocate(self%nnvec)
        call self%grad_shsrch_obj%kill
    end subroutine kill

end module simple_prime3D_srch
