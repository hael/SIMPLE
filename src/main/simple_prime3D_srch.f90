! PRIME3D stochastic search routines
module simple_prime3D_srch
#include "simple_lib.f08"
use simple_oris,                  only: oris
use simple_ori,                   only: ori
use simple_sym,                   only: sym
use simple_polarft_corrcalc,      only: polarft_corrcalc
use simple_pftcc_shsrch,          only: pftcc_shsrch       ! simplex-based angle and shift search
use simple_pftcc_grad_shsrch,     only: pftcc_grad_shsrch  ! gradient-based angle and shift search
implicit none

public :: cleanprime3D_srch, prep4prime3D_srch, prime3D_srch_extr
public :: prime3D_srch, o_peaks
private
#include "simple_local_flags.inc"

real, parameter :: FACTWEIGHTS_THRESH = 0.001 !< threshold for factorial weights

! allocatables for prime3D_srch are class variables to improve caching and reduce alloc overheads 
type(oris), allocatable :: o_peaks(:)                            !< solution objects
real,       allocatable :: proj_space_euls(:,:,:)                !< euler angles
real,       allocatable :: proj_space_shift(:,:,:)               !< shift vector
real,       allocatable :: proj_space_corrs(:,:)                 !< correlations vs. reference orientations
integer,    allocatable :: proj_space_inds(:,:)                  !< stochastic index of reference orientations
integer,    allocatable :: proj_space_state(:,:)                 !< reference orientations state
integer,    allocatable :: proj_space_proj(:,:)                  !< reference orientations projection direction (1 state assumed)
integer,    allocatable :: prev_proj(:)                          !< particle previous reference projection direction
integer,    allocatable :: prev_states(:)                        !< particle previous state
integer,    allocatable :: srch_order(:,:)                       !< stochastic search index
logical,    allocatable :: state_exists(:)                       !< indicates state existence
real,       allocatable :: het_corrs(:,:)                        !< per particle state correlations

type prime3D_srch
    private
    class(polarft_corrcalc), pointer :: pftcc_ptr      => null() !< corrcalc object
    class(oris),             pointer :: a_ptr          => null() !< b%a (primary particle orientation table)
    class(sym),              pointer :: se_ptr         => null() !< b%se (symmetry elements)
    type(pftcc_shsrch)               :: shsrch_obj               !< origin shift search object
!    type(pftcc_grad_shsrch)          :: shsrch_obj               !< activate for gradient-based search object
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
    procedure          :: kill
    procedure          :: exec_prime3D_srch
    procedure          :: exec_prime3D_srch_het
    procedure, private :: greedy_srch
    procedure, private :: greedy_subspace_srch
    procedure, private :: stochastic_srch
    procedure, private :: stochastic_srch_shc
    procedure, private :: stochastic_srch_snhc
    procedure, private :: stochastic_srch_het
    procedure, private :: stochastic_srch_bystate
    procedure, private :: inpl_srch
    procedure, private :: prep4srch
    procedure, private :: prep_npeaks_oris_and_weights
    procedure, private :: store_solution
end type prime3D_srch

contains

    subroutine cleanprime3D_srch( p )
        use simple_params, only: params
        class(params),  intent(in) :: p
        if( allocated(proj_space_euls)  ) deallocate(proj_space_euls)
        if( allocated(proj_space_shift) ) deallocate(proj_space_shift)
        if( allocated(proj_space_corrs) ) deallocate(proj_space_corrs)
        if( allocated(proj_space_inds)  ) deallocate(proj_space_inds)
        if( allocated(proj_space_state) ) deallocate(proj_space_state)
        if( allocated(proj_space_proj)  ) deallocate(proj_space_proj)
        if( allocated(prev_proj)        ) deallocate(prev_proj)
        if( allocated(srch_order)       ) deallocate(srch_order)
        if( allocated(state_exists)     ) deallocate(state_exists)
        if( allocated(het_corrs)        ) deallocate(het_corrs)
        if( allocated(prev_states)      ) deallocate(prev_states)
    end subroutine cleanprime3D_srch

    subroutine prep4prime3D_srch( b, p )
        use simple_ran_tabu, only: ran_tabu
        use simple_params,   only: params
        use simple_build,    only: build
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        type(ran_tabu) :: rt
        real           :: euldists(p%nspace)
        integer        :: i, istate, iproj, iptcl, prev_state, nnnrefs, cnt, prev_ref, nrefs
        ! clean all class arrays
        call cleanprime3D_srch(p)
        if( allocated(o_peaks) )then
            do i=p%fromp,p%top
                call o_peaks(i)%kill
            enddo
            deallocate(o_peaks)
        endif
        ! reference projection directions
        nrefs = p%nstates * p%nspace
        allocate(proj_space_euls( p%fromp:p%top, nrefs,3), source=0.)
        allocate(proj_space_shift(p%fromp:p%top, nrefs,2), source=0.)
        allocate(proj_space_corrs(p%fromp:p%top, nrefs), source=-1.)
        allocate(proj_space_inds( p%fromp:p%top, nrefs), source=0)
        allocate(proj_space_state(p%fromp:p%top, nrefs), source=0)
        allocate(proj_space_proj( p%fromp:p%top, nrefs), source=0)
        do iptcl = p%fromp, p%top
            cnt = 0
            do istate=1,p%nstates
                do iproj=1,p%nspace
                    cnt = cnt + 1
                    proj_space_state(iptcl,cnt)   = istate
                    proj_space_proj(iptcl, cnt)   = iproj
                    proj_space_euls(iptcl, cnt,:) = b%e%get_euler(iproj)
                enddo
            enddo
        enddo
        ! projection direction peaks, eo & CTF transfer
        if( p%refine.ne.'het')then
            allocate(o_peaks(p%fromp:p%top))
            do iptcl = p%fromp, p%top
                call o_peaks(iptcl)%new_clean(p%npeaks)
                ! transfer CTF params
                if( p%ctf.ne.'no' )then
                    call o_peaks(iptcl)%set_all2single('kv', b%a%get(iptcl,'kv'))
                    call o_peaks(iptcl)%set_all2single('cs', b%a%get(iptcl,'cs'))
                    call o_peaks(iptcl)%set_all2single('fraca', b%a%get(iptcl,'fraca'))
                    call o_peaks(iptcl)%set_all2single('dfx', b%a%get(iptcl,'dfx'))
                    if( p%tfplan%mode .eq. 'astig' )then
                        call o_peaks(iptcl)%set_all2single('dfy', b%a%get(iptcl,'dfy'))
                        call o_peaks(iptcl)%set_all2single('angast', b%a%get(iptcl,'angast'))
                    else
                        call o_peaks(iptcl)%set_all2single('dfy', b%a%get(iptcl,'dfx'))
                        call o_peaks(iptcl)%set_all2single('angast', 0.)
                    endif
                    if( p%tfplan%l_phaseplate )then
                        call o_peaks(iptcl)%set_all2single('phshift', b%a%get(iptcl,'phshift'))
                    endif
                endif
                ! transfer eo flag
                call o_peaks(iptcl)%set_all2single('eo', b%a%get(iptcl,'eo'))
            enddo
        endif
        ! search order & previous projection direction
        allocate(prev_proj(p%fromp:p%top), source=0)
            !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
            do iptcl = p%fromp, p%top
                prev_proj(iptcl) = b%e%find_closest_proj(b%a%get_ori(iptcl))
            enddo
            !$omp end parallel do
        select case( trim(p%refine) )
            case( 'het' )
                allocate(het_corrs(p%fromp:p%top,p%nstates),source=-1.)
                allocate(prev_states(p%fromp:p%top),source=0)
            case( 'states' )
                allocate(srch_order(p%fromp:p%top, p%nnn), source=0)
                rt = ran_tabu(p%nnn)
                do iptcl = p%fromp, p%top
                    prev_state          = nint(b%a%get(iptcl, 'state'))
                    srch_order(iptcl,:) = b%nnmat(iptcl,:) + (prev_state - 1) * p%nspace
                    call rt%reset
                    prev_ref = (prev_state - 1) * p%nspace + prev_proj(iptcl)
                    call put_last(prev_ref, srch_order(iptcl,:))
                enddo
            case( 'neigh','shcneigh', 'greedyneigh' )
                nnnrefs =  p%nnn * p%nstates
                allocate(srch_order(p%fromp:p%top, nnnrefs), source=0)
                rt = ran_tabu(nnnrefs)
                do iptcl = p%fromp, p%top   
                    prev_state = nint(b%a%get(iptcl, 'state'))
                    prev_ref   = (prev_state - 1)*p%nspace + prev_proj(iptcl)
                    do istate = 0, p%nstates - 1
                        i = istate * p%nnn + 1
                        srch_order(iptcl,i:i + p%nnn - 1) = b%nnmat(prev_proj(iptcl),:) + istate * p%nspace
                    enddo
                    call rt%reset
                    call rt%shuffle(srch_order(iptcl,:))
                    call put_last(prev_ref, srch_order(iptcl,:))
                enddo
            case( 'yes' )
                allocate(srch_order(p%fromp:p%top,p%nspace*p%nstates), source=0)
                do iptcl = p%fromp, p%top
                    prev_state = nint(b%a%get(iptcl, 'state'))
                    call b%e%calc_euldists(b%a%get_ori(iptcl), euldists)
                    srch_order(iptcl,:p%nspace) = (/(i,i=1,p%nspace)/)
                    call hpsort( p%nspace, euldists, srch_order(iptcl,:p%nspace) )
                    prev_ref = (prev_state - 1) * p%nspace + srch_order(iptcl,1)
                    do istate = 0, p%nstates - 1
                        i = istate * p%nspace + 1
                        srch_order(iptcl,i:i + p%nspace - 1) = srch_order(iptcl,:p%nspace) + istate * p%nspace
                    enddo
                    call put_last(prev_ref, srch_order(iptcl,:))
                enddo
            case('no','shc','snhc','greedy')
                allocate(srch_order(p%fromp:p%top,nrefs), source=0)
                rt = ran_tabu(nrefs)
                do iptcl = p%fromp, p%top
                    call rt%reset
                    call rt%ne_ran_iarr( srch_order(iptcl,:) )
                    prev_state = nint(b%a%get(iptcl, 'state'))
                    prev_ref   = (prev_state - 1) * p%nspace + prev_proj(iptcl)
                    call put_last(prev_ref, srch_order(iptcl,:))
                enddo
            case DEFAULT
                stop 'Unknown refinement mode; simple_hadamard3D_matcher; prep4primesrch3D'
        end select
        call rt%kill
        if( allocated(srch_order) )then
            if( any(srch_order == 0) ) stop 'Invalid index in srch_order; simple_hadamard3D_matcher :: prep4primesrch3D'
        endif
        ! states existence
        if( p%oritab.ne.'' )then
            state_exists = b%a%states_exist(p%nstates)
        else
            allocate(state_exists(p%nstates), source=.true.)
        endif
    end subroutine prep4prime3D_srch

    subroutine prime3D_srch_extr( a, p, do_extr )
        use simple_params, only: params
        class(oris),   intent(inout) :: a
        class(params),  intent(in)   :: p
        logical,        intent(in)   :: do_extr
        real,    allocatable :: curr_corrs(:)
        integer, allocatable :: curr_inds(:)
        real    :: extr_frac, prev_corr, frac, mi_state, corr_thresh, avg_rndcorr, avg_shccorr, avg_currcorr
        real    :: avg_prevrndcorr, avg_prevshccorr
        integer :: iextr_lim, iptcl, n_incl, n_rnd, cnt, ind, istate, nrefs_eval, prev_state, zero_pop
        if( .not.do_extr )return
        zero_pop  = a%get_pop( 0, 'state', consider_w=.false.)
        iextr_lim = ceiling(2.*log( real(p%nptcls-zero_pop) ))
        n_incl    = count(prev_states(p%fromp:p%top)>0)      ! # non-zero state particles in part
        extr_frac = EXTRINITHRESH * cos(PI/2. * real(p%extr_iter-1)/real(iextr_lim)) ! cosine decay
        extr_frac = min(EXTRINITHRESH, max(0.0, extr_frac)) ! fraction of particles to randomize
        n_rnd     = nint(real(n_incl)*extr_frac)            ! # of particles to randomize
        ! collate and sort correlations
        allocate(curr_corrs(n_incl), curr_inds(n_incl))
        cnt = 0
        do iptcl = p%fromp, p%top
            if(prev_states(iptcl) > 0)then
                cnt = cnt + 1 
                curr_corrs(cnt) = het_corrs(iptcl,prev_states(iptcl))
                curr_inds(cnt)  = iptcl
            endif
        enddo
        avg_currcorr = sum(curr_corrs) / real(n_incl)
        if( n_rnd == 0 )then
            corr_thresh = 0.
        else
            call hpsort(n_incl, curr_corrs, curr_inds)
        endif
        ! assign moves
        avg_rndcorr     = 0.
        avg_shccorr     = 0.
        avg_prevrndcorr = 0.
        avg_prevshccorr = 0.
        do ind = 1, n_incl
            iptcl      = curr_inds(ind)
            prev_corr  = het_corrs(iptcl,prev_states(iptcl))
            prev_state = prev_states(iptcl)
            if( ind <= n_rnd )then
                ! stochastic diversifying move
                istate = irnd_uni(p%nstates)
                do while(istate == prev_state .or. .not.state_exists(istate))
                    istate = irnd_uni(p%nstates)
                enddo
                nrefs_eval      = 1
                avg_prevrndcorr = avg_prevrndcorr + prev_corr
                avg_rndcorr     = avg_rndcorr + het_corrs(iptcl, istate)
                if( ind == n_rnd )corr_thresh = het_corrs(iptcl, prev_state)
            else
                ! SHC move
                istate          = shcloc(p%nstates, het_corrs(iptcl,:), prev_corr)
                nrefs_eval      = count(het_corrs(iptcl,:) <= prev_corr)
                avg_prevshccorr = avg_prevshccorr + prev_corr
                avg_shccorr     = avg_shccorr + het_corrs(iptcl, istate)
            endif
            ! orientation update
            frac = 100.*real(nrefs_eval) / real(p%nstates)
            if( prev_state .ne. istate )then
                mi_state = 0.
            else
                mi_state = 1.
            endif      
            call a%set(iptcl, 'state',    real(istate))
            call a%set(iptcl, 'corr',     het_corrs(iptcl, istate))
            call a%set(iptcl, 'frac',     frac)
            call a%set(iptcl, 'mi_state', mi_state)
            call a%set(iptcl, 'mi_joint', mi_state)
            call a%set(iptcl, 'mi_proj',  1.)
            call a%set(iptcl, 'mi_inpl',  1.)
            call a%set(iptcl, 'w',        1.)
            call a%set(iptcl, 'ow',       1.)
        enddo
        write(*,'(A,F8.3)') '>>> PARTICLE RANDOMIZATION(%): ', 100.*extr_frac
        write(*,'(A,F8.3)') '>>> CORRELATION THRESHOLD    : ', corr_thresh
        write(*,'(A,I8)')   '>>> # RANDOMIZED PARTICLES   : ', n_rnd
        write(*,'(A,I8)')   '>>> # SHC        PARTICLES   : ', n_incl-n_rnd
        write(*,'(A,F8.3)') '>>> AVG PREV     CORRELATION : ', avg_currcorr
        if(n_rnd > 0)then
            avg_rndcorr     = avg_rndcorr / real(n_rnd)
            avg_prevrndcorr = avg_prevrndcorr / real(n_rnd)
            write(*,'(A,F8.3)') '>>> AVG PREV RND CORRELATION : ', avg_prevrndcorr
            write(*,'(A,F8.3)') '>>> AVG      RND CORRELATION : ', avg_rndcorr
        endif
        avg_shccorr     = avg_shccorr / real(n_incl-n_rnd)
        avg_prevshccorr = avg_prevshccorr / real(n_incl-n_rnd)
        write(*,'(A,F8.3)') '>>> AVG PREV SHC CORRELATION : ', avg_prevshccorr
        write(*,'(A,F8.3)') '>>> AVG      SHC CORRELATION : ', avg_shccorr
        deallocate(curr_corrs, curr_inds)
    end subroutine prime3D_srch_extr

    subroutine new( self, iptcl, p, pftcc, a, se )
        use simple_params, only: params
        class(prime3D_srch),             intent(inout) :: self   !< instance
        integer,                         intent(in)    :: iptcl  !< global particle index
        class(params),                   intent(in)    :: p      !< parameters
        class(polarft_corrcalc), target, intent(in)    :: pftcc  !< corrcalc object
        class(oris),             target, intent(in)    :: a      !< b%a (primary particle orientation table)
        class(sym),              target, intent(in)    :: se     !< b%se (symmetry elements)
        integer :: alloc_stat, nstates_eff
        real    :: lims(2,2), lims_init(2,2)
        ! set constants
        self%pftcc_ptr   => pftcc
        self%a_ptr       => a
        self%se_ptr      => se
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
            nstates_eff = count(state_exists)
            select case(trim(p%refine))
                case('het')
                    self%npeaks_grid = 1
                case('states')
                    self%npeaks_grid = GRIDNPEAKS
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
        self%doshift   = p%doshift
        ! create in-plane search objects
        lims(:,1)      = -p%trs
        lims(:,2)      =  p%trs
        lims_init(:,1) = -SHC_INPL_TRSHWDTH
        lims_init(:,2) =  SHC_INPL_TRSHWDTH
        call self%shsrch_obj%new(self%pftcc_ptr, lims, lims_init=lims_init,&
            &shbarrier=self%shbarr, nrestarts=3, maxits=60)
!!$  activate below for gradient-based shift search
!!$        call self%shsrch_obj%new(self%pftcc_ptr, lims, lims_init=lims_init,&
!!$            &shbarrier=self%shbarr, nrestarts=1, maxits=60)
        self%exists = .true.
        DebugPrint '>>> PRIME3D_SRCH::CONSTRUCTED NEW SIMPLE_PRIME3D_SRCH OBJECT'
    end subroutine new

    subroutine kill( self )
        class(prime3D_srch),  intent(inout) :: self   !< instance
        call self%shsrch_obj%kill
    end subroutine kill

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
    subroutine exec_prime3D_srch_het( self, do_extr )
        class(prime3D_srch), intent(inout) :: self
        logical,             intent(in)    :: do_extr
        call self%stochastic_srch_het( do_extr )
        DebugPrint '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SRCH_HET'
    end subroutine exec_prime3D_srch_het

    !>  \brief  Individual stochastic search by state
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    subroutine stochastic_srch_bystate( self, lp, nnmat )
        class(prime3D_srch),     intent(inout) :: self
        real,                    intent(in)    :: lp
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn_static)
        integer :: iref, isample, inpl_ind(1)
        real    :: inpl_corrs(self%nrots), corrs(self%nrefs), inpl_corr
        ! execute search
        if( nint(self%a_ptr%get(self%iptcl,'state')) > 0 )then
            ! initialize
            call self%prep4srch(lp, nnmat=nnmat)
            self%nbetter    = 0
            self%nrefs_eval = 0
            proj_space_corrs(self%iptcl,:) = -1.
            ! search
            do isample = 1, self%nnn
                iref = srch_order(self%iptcl, isample) ! set the stochastic reference index
                call per_ref_srch                      ! actual search
                if( self%nbetter >= self%npeaks ) exit ! exit condition
            end do
            ! sort in correlation projection direction space
            corrs = proj_space_corrs(self%iptcl,:)
            call hpsort(self%nrefs, corrs, proj_space_inds(self%iptcl,:))
            call self%inpl_srch ! search shifts
            ! prepare weights and orientations
            call self%prep_npeaks_oris_and_weights
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED STOCHASTIC BYSTATE SEARCH'

        contains

            subroutine per_ref_srch
                if( proj_space_state(self%iptcl,iref) .ne. self%prev_state )then
                    ! checker that will need to go
                    print *,self%iptcl, self%prev_ref, prev_proj(self%iptcl), self%prev_state, iref, proj_space_state(self%iptcl,iref)
                    stop 'srch order error; simple_prime3d_srch; stochastic_srch_bystate'
                endif
                call self%pftcc_ptr%gencorrs(iref, self%iptcl, inpl_corrs) ! In-plane correlations
                inpl_ind   = maxloc(inpl_corrs)                            ! greedy in-plane index
                inpl_corr  = inpl_corrs(inpl_ind(1))                       ! max in plane correlation
                call self%store_solution(iref, iref, inpl_ind(1), inpl_corr)
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
        integer :: iref,isample,nrefs,target_projs(self%npeaks_grid), inpl_ind(1), state
        real    :: corrs(self%nrefs), inpl_corrs(self%nrots), inpl_corr
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
            self%nbetter    = 0
            self%nrefs_eval = 0
            proj_space_corrs(self%iptcl,:) = -1.
            ! search
            do isample=1,nrefs
                iref = srch_order(self%iptcl,isample) ! set the reference index
                call per_ref_srch                     ! actual search
            end do
            ! in greedy mode, we evaluate all refs
            self%nrefs_eval = nrefs
            ! sort in correlation projection direction space
            corrs = proj_space_corrs(self%iptcl,:)
            call hpsort(self%nrefs, corrs, proj_space_inds(self%iptcl,:))
            ! take care of the in-planes
            call self%inpl_srch
            ! prepare weights & orientation
            call self%prep_npeaks_oris_and_weights
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED GREEDY SEARCH'

        contains

            subroutine per_ref_srch
                state = proj_space_state(self%iptcl,iref)
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
        integer   :: state_cnt(self%nstates), iref_state, loc(1), inpl_ind
        if( nint(self%a_ptr%get(self%iptcl,'state')) > 0 )then
            ! initialize
            target_projs = 0
            proj_space_corrs(self%iptcl,:) = -1.
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
            corrs = proj_space_corrs(self%iptcl,:)
            call hpsort(self%nrefs, corrs, proj_space_inds(self%iptcl,:))
            ! return target points
            ntargets = size(target_projs)
            cnt = 1
            target_projs( cnt ) = prev_proj(self%iptcl) ! previous always part of the targets
            if( self%nstates == 1 )then
                ! Single state
                do isample=self%nrefs,self%nrefs - ntargets + 1,-1
                    if( target_projs(1) == proj_space_inds(self%iptcl,isample) )then
                        ! direction is already in target set
                    else
                        cnt = cnt + 1
                        target_projs(cnt) = proj_space_inds(self%iptcl,isample)
                        if( cnt == ntargets ) exit
                    endif
                end do
            else
                ! Multiples states
                state_cnt = 1                                        ! previous always part of the targets
                do isample = self%nrefs, 1, -1
                    if( cnt >= self%npeaks_grid )exit                ! all that we need
                    iref_state = proj_space_inds(self%iptcl,isample) ! reference index to multi-state space
                    istate     = ceiling(real(iref_state)/real(self%nprojs))
                    iref       = iref_state - (istate-1)*self%nprojs ! reference index to single state space
                    if( .not.state_exists(istate) )cycle
                    if( any(target_projs == iref) )cycle             ! direction is already set
                    if( state_cnt(istate) >= GRIDNPEAKS )cycle       ! state is already filled
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

            subroutine per_ref_srch
                if( state_exists(istate) )then
                    call self%pftcc_ptr%gencorrs(iref, self%iptcl, self%kstop_grid, inpl_corrs) ! In-plane correlations
                    loc        = maxloc(inpl_corrs)               ! greedy in-plane
                    inpl_ind   = loc(1)                           ! in-plane angle index
                    inpl_corr  = inpl_corrs(inpl_ind)             ! max in plane correlation
                    proj_space_corrs(self%iptcl,iref) = inpl_corr ! stash in-plane correlation for sorting
                    proj_space_inds(self%iptcl,iref)  = iref      ! stash the index for sorting
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
        integer :: iref,isample,nrefs,target_projs(self%npeaks_grid), loc(1), inpl_ind
        real    :: corrs(self%nrefs), inpl_corrs(self%nrots), inpl_corr
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
            ! initialize, ctd
            self%nbetter    = 0
            self%nrefs_eval = 0
            proj_space_corrs(self%iptcl,:) = -1.
            do isample=1,nrefs
                iref = srch_order(self%iptcl,isample)  ! set the stochastic reference index
                call per_ref_srch                      ! actual search
                if( self%nbetter >= self%npeaks ) exit ! exit condition
            end do
            ! sort in correlation projection direction space
            corrs = proj_space_corrs(self%iptcl,:)
            call hpsort(self%nrefs, corrs, proj_space_inds(self%iptcl,:))
            call self%inpl_srch ! search shifts
            ! prepare weights and orientations
            call self%prep_npeaks_oris_and_weights
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH'

        contains

            subroutine per_ref_srch
                if( state_exists( proj_space_state(self%iptcl,iref) ) )then
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
    !! \param lp low-pass cutoff freq
    !! \param nnmat nearest neighbour matrix
    !! \param grid_projs projection indices for grid search
    subroutine stochastic_srch_shc( self, lp, nnmat, grid_projs )
        class(prime3D_srch), intent(inout) :: self
        real,                intent(in)    :: lp
        integer, optional,   intent(in)    :: nnmat(self%nprojs,self%nnn_static), grid_projs(:)
        real    :: inpl_corr, inpl_corrs(self%nrots)
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
            found_better    = .false.
            self%nrefs_eval = 0
            ! search
            do isample=1,nrefs
                iref = srch_order(self%iptcl, isample ) ! stochastic reference index
                call per_ref_srch                       ! actual search
                if( inpl_ind > 0 )then
                    ! correlation-improving reference found, store it
                    call self%store_solution(self%nrefs, iref, inpl_ind, inpl_corr)
                    found_better = .true.               ! flag for found solution
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
        DebugPrint '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH (REFINE=SHC|SHCNEIGH)'

        contains

            subroutine per_ref_srch
                use simple_rnd, only: shcloc
                inpl_ind  = 0                                                    ! init in-plane index
                inpl_corr = 0.                                                   ! init correlation
                if( state_exists( proj_space_state(self%iptcl,iref) ) )then
                    call self%pftcc_ptr%gencorrs( iref, self%iptcl, inpl_corrs ) ! in-plane correlations
                    inpl_ind = shcloc(self%nrots, inpl_corrs, self%prev_corr)    ! first improving in-plane index
                    if( inpl_ind > 0 ) inpl_corr = inpl_corrs( inpl_ind )        ! improving correlation found
                    self%nrefs_eval = self%nrefs_eval + 1                        ! updates fractional search space
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
        integer :: iref, isample, loc(1), inpl_ind
        real    :: corrs(self%nrefs), inpl_corrs(self%nrots), inpl_corr
        if( nint(self%a_ptr%get(self%iptcl,'state')) > 0 )then
            ! initialize
            call self%prep4srch(lp)
            self%nbetter    = 0
            self%nrefs_eval = 0
            proj_space_corrs(self%iptcl,:) = -1.
            ! search
            do isample=1,szsn
                iref = srch_order(self%iptcl,isample) ! set the stochastic reference index
                call per_ref_srch                     ! actual search
            end do
            self%nrefs_eval = szsn
            ! sort in correlation projection direction space
            corrs = proj_space_corrs(self%iptcl,:)
            call hpsort(self%nrefs, corrs, proj_space_inds(self%iptcl,:))
            ! output
            call self%prep_npeaks_oris_and_weights
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED EXTREMAL SEARCH'

        contains

            subroutine per_ref_srch
                if( state_exists( proj_space_state(self%iptcl,iref) ) )then
                    call self%pftcc_ptr%gencorrs(iref, self%iptcl, inpl_corrs) ! In-plane correlations
                    loc        = maxloc(inpl_corrs)               ! greedy in-plane
                    inpl_ind   = loc(1)                           ! in-plane angle index
                    inpl_corr  = inpl_corrs(inpl_ind)             ! max in plane correlation
                    call self%store_solution(iref, iref, inpl_ind, inpl_corr)
                endif
            end subroutine per_ref_srch

    end subroutine stochastic_srch_snhc

    !> stochastic search het
    !! \param extr_bound corr threshold
    !! \param statecnt state counter array
    subroutine stochastic_srch_het( self, do_extr )
        use simple_rnd, only: shcloc, irnd_uni
        class(prime3D_srch), intent(inout) :: self
        logical,             intent(in)    :: do_extr
        integer :: iref, state
        real    :: corr, mi_state, mi_inpl, frac
        real    :: corrs_state(self%nstates), corrs_inpl(self%nrots), shvec(2)
        prev_states(self%iptcl) = nint(self%a_ptr%get(self%iptcl, 'state'))
        if( prev_states(self%iptcl) > 0 )then
            if( prev_states(self%iptcl) > self%nstates )&
                &stop 'previous best state outside boundary; stochastic_srch_het; simple_prime3D_srch'
            if( .not. state_exists(prev_states(self%iptcl)) )&
                &stop 'empty previous state; stochastic_srch_het; simple_prime3D_srch'
            ! all states correlations
            self%prev_roind = self%pftcc_ptr%get_roind(360.-self%a_ptr%e3get(self%iptcl))
            corrs_state = -1.
            do state = 1, self%nstates
                if( .not.state_exists(state) )cycle
                iref = (state-1) * self%nprojs + prev_proj(self%iptcl)
                call self%pftcc_ptr%gencorrs(iref, self%iptcl, corrs_inpl)
                corrs_state(state) = corrs_inpl(self%prev_roind)
            enddo
            if( do_extr )then
                ! populates state correlations matrix
                het_corrs(self%iptcl,:) = corrs_state(:)
            else
                ! SHC out of plane
                self%prev_corr  = corrs_state(prev_states(self%iptcl))
                state           = shcloc(self%nstates, corrs_state, self%prev_corr)
                corr            = corrs_state(state)
                self%nrefs_eval = count(corrs_state <= self%prev_corr)
                ! state           = maxloc(corrs_state) ! greedy in state
                ! self%nrefs_eval = self%nstates
                ! greedy inpl move
                iref = (state-1)*self%nprojs + prev_proj(self%iptcl)
                proj_space_corrs(self%iptcl,iref)      = corr
                proj_space_inds(self%iptcl,self%nrefs) = iref
                call self%inpl_srch
                if( proj_space_corrs(self%iptcl,iref) > corr )then
                    corr  = proj_space_corrs(self%iptcl,iref)
                    shvec = self%a_ptr%get_2Dshift(self%iptcl) + proj_space_shift(self%iptcl,iref,:)
                    call self%a_ptr%e3set(self%iptcl, proj_space_euls(self%iptcl, iref, 3))
                    call self%a_ptr%set_shift(self%iptcl, shvec)
                endif
                ! updates peaks and orientation orientation
                frac    = 100.*real(self%nrefs_eval) / real(self%nstates)
                mi_inpl = 1.
                if( self%prev_roind .ne. self%pftcc_ptr%get_roind(360.-proj_space_euls(self%iptcl, iref, 3)) )then
                    mi_inpl = 0.
                endif
                call self%a_ptr%set(self%iptcl,'frac', frac)
                call self%a_ptr%set(self%iptcl,'state', real(state))
                call self%a_ptr%set(self%iptcl,'corr', corr)
                call self%a_ptr%set(self%iptcl,'mi_proj', 1.)
                call self%a_ptr%set(self%iptcl,'mi_inpl', mi_inpl)
                if( prev_states(self%iptcl) .ne. state )then
                    mi_state = 0.
                else
                    mi_state = 1.
                endif
                call self%a_ptr%set(self%iptcl,'mi_state', mi_state)
                call self%a_ptr%set(self%iptcl,'mi_joint', (mi_state+mi_inpl)/2.)
                call self%a_ptr%set(self%iptcl,'w', 1.)
            endif
        else
            call self%a_ptr%reject(self%iptcl)
        endif
        DebugPrint '>>> PRIME3D_SRCH::FINISHED HET SEARCH'
    end subroutine stochastic_srch_het

    !>  \brief  executes the in-plane search. This improved routine took HCN
    !!          from stall @ 5.4 to close to convergence @ 4.5 after 10 iters
    !!          refine=no with 1000 references
    subroutine inpl_srch( self )
        class(prime3D_srch), intent(inout) :: self
        real, allocatable :: cxy(:), crxy(:)
        integer :: i, ref, irot
        if( self%doshift )then
            do i=self%nrefs,self%nrefs-self%npeaks+1,-1
                ref = proj_space_inds(self%iptcl, i)
                call self%shsrch_obj%set_indices(ref, self%iptcl)
                cxy = self%shsrch_obj%minimize(irot=irot)
                if( irot > 0 )then
                    ! irot > 0 guarantees improvement found
                    proj_space_euls(self%iptcl, ref,3) = 360. - self%pftcc_ptr%get_rot(irot)
                    proj_space_corrs(self%iptcl,ref)   = cxy(1)
                    proj_space_shift(self%iptcl,ref,:) = cxy(2:3)
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
        real      :: corrs(self%nrots)
        type(ori) :: o_prev
        real      :: cc_t_min_1, corr
        if( str_has_substr(self%refine,'neigh') )then
            if( .not. present(nnmat) )&
            &stop 'need optional nnmat to be present for refine=neigh modes :: prep4srch (prime3D_srch)'
            if( .not. present(target_projs) )&
            &stop 'need optional target_projs to be present for refine=neigh modes :: prep4srch (prime3D_srch)'
        endif
        o_prev          = self%a_ptr%get_ori(self%iptcl)
        self%prev_state = nint(o_prev%get('state'))                        ! state index
        self%prev_roind = self%pftcc_ptr%get_roind(360.-o_prev%e3get())    ! in-plane angle index
        self%prev_shvec = o_prev%get_2Dshift()                             ! shift vector
        self%prev_ref   = (self%prev_state-1)*self%nprojs + prev_proj(self%iptcl) ! previous reference
        if( self%prev_state > 0 )then
            if( self%prev_state > self%nstates ) stop 'previous best state outside boundary; prep4srch; simple_prime3D_srch'
            if( .not. state_exists(self%prev_state) ) stop 'empty previous state; prep4srch; simple_prime3D_srch'
        endif
        select case( self%refine )
            case( 'no','shc','snhc','greedy','yes')
                ! 
            case( 'neigh','shcneigh', 'greedyneigh' )
                ! disjoint nearest neighbour set
                nnvec = merge_into_disjoint_set(self%nprojs, self%nnn_static, nnmat, target_projs)
            case( 'het' )
                !
            case( 'states' )
                !
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
        neff_states = count(state_exists)
        ! init npeaks
        do ipeak = 1, self%npeaks
            cnt = self%nrefs - self%npeaks + ipeak
            ref = proj_space_inds(self%iptcl, cnt )
            if( ref < 1 .or. ref > self%nrefs )then
                print *, 'ref: ', ref
                stop 'ref index out of bound; simple_prime3D_srch::prep_npeaks_oris'
            endif
            state = proj_space_state(self%iptcl,ref)
            if( .not. state_exists(state) )then
                print *, 'empty state: ', state
                stop 'simple_prime3D_srch::prep_npeaks_oris'
            endif
            ! add shift
            shvec = self%prev_shvec
            if( self%doshift )shvec = shvec + proj_space_shift(self%iptcl,ref,1:2)
            where( abs(shvec) < 1e-6 ) shvec = 0.
            ! transfer to solution set
            corrs(ipeak) = proj_space_corrs(self%iptcl,ref)
            call o_peaks(self%iptcl)%set(ipeak, 'state', real(state))
            call o_peaks(self%iptcl)%set(ipeak, 'proj',  real(proj_space_proj(self%iptcl,ref)))
            call o_peaks(self%iptcl)%set(ipeak, 'corr',  corrs(ipeak))
            call o_peaks(self%iptcl)%set_euler(ipeak, proj_space_euls(self%iptcl,ref,1:3))
            call o_peaks(self%iptcl)%set_shift(ipeak, shvec)
        enddo
        best_loc = maxloc(corrs)
        ! stochastic weights
        if( self%npeaks == 1 )then
            call o_peaks(self%iptcl)%set(1,'ow',1.0)
            wcorr = o_peaks(self%iptcl)%get(1,'corr')
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
            ws = ws / sum(ws,mask=included)
            ! weighted corr
            wcorr = sum(ws*corrs,mask=included)
            ! update npeaks individual weights
            call o_peaks(self%iptcl)%set_all('ow', ws)
        endif
        ! angular standard deviation
        ang_sdev = 0.
        if( trim(self%se_ptr%get_pgrp()).eq.'c1' )then
            ang_sdev = o_peaks(self%iptcl)%ang_sdev(self%refine, self%nstates, self%npeaks)
        else
            if( self%npeaks > 2 )then
                loc    = maxloc(corrs)
                sym_os = o_peaks(self%iptcl)
                do ipeak = 1, self%npeaks
                    if(ipeak == loc(1))cycle
                    o = o_peaks(self%iptcl)%get_ori(ipeak)
                    call self%se_ptr%sym_dists( o_peaks(self%iptcl)%get_ori(best_loc(1)), o, dist, inpl_dist)
                    call sym_os%set_ori(ipeak, o)
                enddo
                ang_sdev = sym_os%ang_sdev(self%refine, self%nstates, self%npeaks)
            endif
        endif
        ! Update the best orientation
        ! angular distances
        call self%se_ptr%sym_dists( self%a_ptr%get_ori(self%iptcl),&
            &o_peaks(self%iptcl)%get_ori(best_loc(1)), euldist, dist_inpl )
        ! convergence parameters
        state = nint( o_peaks(self%iptcl)%get(best_loc(1), 'state') )
        if( .not. state_exists(state) )then
            print *, 'empty state: ', state
            stop 'simple_prime3d_srch; update_best'
        endif
        roind = self%pftcc_ptr%get_roind( 360.-o_peaks(self%iptcl)%e3get(best_loc(1)) )
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
        DebugPrint '>>> PRIME3D_SRCH::EXECUTED PREP_NPEAKS_ORIS'
    end subroutine prep_npeaks_oris_and_weights

    subroutine store_solution( self, ind, ref, inpl_ind, corr )
        class(prime3D_srch),     intent(inout) :: self
        integer,                 intent(in)    :: ind, ref, inpl_ind
        real,                    intent(in)    :: corr
        proj_space_inds(self%iptcl,ind)   = ref
        proj_space_euls(self%iptcl,ref,3) = 360. - self%pftcc_ptr%get_rot(inpl_ind)
        proj_space_corrs(self%iptcl,ref)  = corr
    end subroutine store_solution

end module simple_prime3D_srch
