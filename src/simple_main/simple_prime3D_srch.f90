module simple_prime3D_srch
use simple_defs
use simple_prime_srch,       only: prime_srch
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_strings,          only: str_has_substr
use simple_sym,              only: sym
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_pftcc_shsrch,     only: pftcc_shsrch
use simple_pftcc_inplsrch,   only: pftcc_inplsrch
use simple_math             ! use all in there
implicit none

public :: prime3D_srch
private

logical, parameter :: DEBUG = .false.

type prime3D_srch
    private
    type(oris)              :: e                       !< copy oforientations of references in proj_space
    type(prime_srch)        :: srch_common             !< functionalities common to primesrch2D/3D
    type(oris)              :: o_refs                  !< projection directions search space
    type(oris)              :: o_peaks                 !< orientations of best npeaks oris
    type(sym)               :: se                      !< symmetry object
    type(pftcc_shsrch)      :: shsrch_obj              !< origin shift search object
    type(pftcc_inplsrch)    :: inplsrch_obj            !< in-plane search object
    integer                 :: nrefs          = 0      !< total number of references (nstates*nprojs)
    integer                 :: nnnrefs        = 0      !< total number of neighboring references (nstates*nnn)
    integer                 :: nstates        = 0      !< number of states
    integer                 :: nprojs         = 0      !< number of projections (same as number of oris in o_refs)
    integer                 :: nrots          = 0      !< number of in-plane rotations in polar representation
    integer                 :: npeaks         = 0      !< number of peaks (nonzero orientation weights)
    integer                 :: nbetter        = 0      !< nr of better orientations identified
    integer                 :: nrefs_eval     = 0      !< nr of references evaluated
    integer                 :: nnn            = 0      !< number of nearest neighbors
    integer                 :: nsym           = 0      !< symmetry order
    integer                 :: prev_roind     = 0      !< previous in-plane rotation index
    integer                 :: prev_state     = 0      !< previous state index
    integer                 :: prev_ref       = 0      !< previous reference index
    integer                 :: prev_proj      = 0      !< previous projection index
    real                    :: prev_corr      = 1.     !< previous best correlation
    real                    :: prev_shvec(2)  = 0.     !< previous origin shift vector
    real                    :: lims(2,2)      = 0.     !< shift search range limit
    real                    :: eullims(3,2)   = 0.     !< asymetric unit bounds
    real                    :: athres         = 0.     !< angular threshold
    real                    :: dfx            = 0.     !< ctf x-defocus
    real                    :: dfy            = 0.     !< ctf y-defocus
    real                    :: angast         = 0.     !< ctf astigmatism
    integer, allocatable    :: proj_space_inds(:)      !< projection space index array
    integer, allocatable    :: srch_order(:)           !< stochastic search order
    integer, allocatable    :: nnmat_sym(:,:)          !< nearest neighbor matrix for adaptive symmetry refinement
    logical, allocatable    :: state_exists(:)         !< indicates whether each state is populated
    character(len=STDLEN)   :: refine        = ''      !< refinement flag
    character(len=STDLEN)   :: ctf           = ''      !< ctf flag
    character(len=STDLEN)   :: shbarr        = ''      !< shift barrier flag
    logical                 :: doshift       = .true.  !< origin shift search indicator
    logical                 :: greedy_inpl   = .true.  !< indicator for whether in-plane search is greedy or not
    logical                 :: exists        = .false. !< 2 indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! CALCULATORS
    procedure          :: prep_npeaks_oris
    procedure          :: sort_shifted_peaks
    ! GETTERS
    procedure          :: get_ori_best
    procedure          :: get_ori
    procedure          :: get_oris
    procedure          :: ang_sdev
    ! PREPARATION ROUTINES
    procedure          :: prep4srch
    procedure          :: prep_corr4srch
    procedure          :: prep_reforis
    ! SEARCH ROUTINES
    procedure          :: exec_prime3D_srch
    procedure          :: exec_prime3D_srch_shc
    procedure          :: exec_prime3D_srch_het
    procedure, private :: stochastic_srch
    procedure, private :: stochastic_srch_shc
    procedure, private :: stochastic_srch_het
    procedure, private :: gen_symnnmat
    procedure          :: inpl_srch
    procedure, private :: greedy_srch
    procedure          :: stochastic_weights
    ! GETTERS & SETTERS
    procedure          :: store_solution
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
    ! DESTRUCTOR
    procedure          :: kill
end type prime3D_srch

contains

    !>  \brief  is a constructor
    subroutine new( self, a, e, p, pftcc )
        use simple_params, only: params
        class(prime3D_srch),             intent(inout) :: self  !< instance
        class(oris),                     intent(inout) :: a     !< ptcls oris
        class(oris),                     intent(in)    :: e     !< reference oris
        class(params),                   intent(in)    :: p     !< parameters
        class(polarft_corrcalc), target, intent(inout) :: pftcc !< correlator
        real,    allocatable :: corrs(:)
        integer, allocatable :: inds(:)
        integer  :: alloc_stat, s, ind, thresh_ind
        ! destroy possibly pre-existing instance
        call self%kill
        ! set constants
        self%e           =  e
        self%nstates     =  p%nstates
        self%nprojs      =  p%nspace
        self%nrefs       =  self%nprojs*self%nstates
        self%nrots       =  round2even(twopi*real(p%ring2))
        self%srch_common =  prime_srch(p)
        self%npeaks      =  p%npeaks
        self%nbetter     =  0
        self%nrefs_eval  =  0
        self%athres      =  p%athres
        self%doshift     =  p%doshift
        self%refine      =  p%refine
        self%ctf         =  p%ctf
        self%nnn         =  p%nnn
        self%nnnrefs     =  self%nnn*self%nstates
        self%shbarr      =  p%shbarrier
        self%greedy_inpl = .true.
        if( str_has_substr(self%refine,'shc') )then
            if( self%npeaks > 1 ) stop 'npeaks must be equal to 1 with refine=shc|shcneigh'
            self%greedy_inpl = .false.
        endif
        ! construct composites
        allocate(self%state_exists(self%nstates), self%proj_space_inds(self%nrefs), stat=alloc_stat)   
        call alloc_err('In: new; simple_prime3D_srch, 1', alloc_stat)
        self%proj_space_inds = 0
        if( p%oritab.ne.'' )then
            self%state_exists = .false.
            do s=1,self%nstates
                if( a%get_statepop(s) > 0 ) self%state_exists(s) = .true.
            enddo
        else
            self%state_exists = .true.
        endif
        ! generate oris oject in which the best npeaks refs will be stored
        call self%o_peaks%new(self%npeaks)
        ! symmetry-related variables
        call self%se%new( p%pgrp )
        self%nsym    = self%se%get_nsym()
        self%eullims = self%se%srchrange()
        if( str_has_substr(self%refine,'adasym') ) call self%gen_symnnmat
        ! updates option to search shift
        self%doshift = p%doshift
        if( self%doshift )then
            self%lims(:,1) = -p%trs
            self%lims(:,2) =  p%trs
        endif
        self%exists = .true.
        ! create in-plane search objects
        call self%shsrch_obj%new(pftcc, self%lims, self%shbarr)
        call self%inplsrch_obj%new(pftcc, self%lims, self%shbarr)
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::CONSTRUCTED NEW SIMPLE_PRIME3D_SRCH OBJECT'
    end subroutine new
    
    ! CALCULATORS

    !> \brief  for sorting already shifted npeaks oris
    subroutine sort_shifted_peaks( self )
        class(prime3D_srch), intent(inout) :: self
        type(oris)        :: o_sorted
        type(ori)         :: o
        real, allocatable :: corrs(:)
        integer           :: i, inds(self%npeaks)
        if( self%o_peaks%get_noris() /= self%npeaks ) stop 'invalid number of oris in simple_prime3D_srch :: sort_shifted_peaks'
        call o_sorted%new( self%npeaks )
        corrs = self%o_peaks%get_all('corr')
        inds  = (/(i,i=1,self%npeaks)/)
        call hpsort(self%npeaks, corrs, inds)
        do i=1,self%npeaks
            o = self%o_peaks%get_ori( inds(i) )
            call o_sorted%set_ori( i,o )
        enddo
        self%o_peaks = o_sorted
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED SORT_SHIFTED_PEAKS'
    end subroutine sort_shifted_peaks

    ! GETTERS & SETTERS
    
    !>  \brief  to get the best orientation
    subroutine get_ori_best( self, o2update )
        use simple_math, only: myacos, rotmat2d, rad2deg
        class(prime3D_srch), intent(inout) :: self
        class(ori),          intent(inout) :: o2update
        type(ori) :: o
        real      :: euldist, mi_joint, mi_proj, mi_inpl, mi_state
        real      :: mat(2,2), u(2), x1(2), x2(2)
        integer   :: roind,state
        ! make unit vector
        u(1) = 0.
        u(2) = 1.
        ! calculate previous vec
        mat  = rotmat2d(o2update%e3get())
        x1   = matmul(u,mat)
        ! get new orientation
        o = self%o_peaks%get_ori( self%npeaks )
        ! calculate new vec
        mat     = rotmat2d(o%e3get())
        x2      = matmul(u,mat)
        euldist = rad2deg( o2update.euldist.o )
        roind   = self%srch_common%roind( 360.-o%e3get() )
        state   = nint( o%get('state') )
        if( .not.self%state_exists(state) )then
            print *,'Empty state in simple_prime3d_srch; get_ori_best'
            stop
        endif
        ! calculate overlap between distributions
        mi_proj  = 0.
        mi_inpl  = 0.
        mi_state = 0.
        mi_joint = 0.
        if( euldist < 0.1 )then
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
        call o2update%set('mi_proj',  mi_proj)
        call o2update%set('mi_inpl',  mi_inpl)
        call o2update%set('mi_state', mi_state)
        call o2update%set('mi_joint', mi_joint)
        ! set the distances before we update the orientation
        call o2update%set('dist', 0.5*euldist + 0.5*o2update%get('dist'))
        call o2update%set('dist_inpl', rad2deg(myacos(dot_product(x1,x2))))
        ! all the other stuff
        call o2update%set_euler( o%get_euler()    )
        call o2update%set_shift( o%get_shift()    )
        call o2update%set( 'state', real(state)   )
        call o2update%set( 'frac',  o%get('frac') )
        call o2update%set( 'corr',  o%get('corr') )
        call o2update%set( 'ow',    o%get('ow')   )
        call o2update%set( 'mirr',  0.            )
        call o2update%set( 'proj',  o%get('proj') )
        call o2update%set( 'sdev',  o%get('sdev') )
        ! stash and return
        call self%o_peaks%set_ori(self%npeaks, o2update)
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::GOT BEST ORI'
    end subroutine get_ori_best

    !>  \brief  to produce oris object for probabilistic reconstruction
    subroutine get_oris( self, os, o )
        class(prime3D_srch), intent(inout) :: self
        class(oris),         intent(inout) :: os
        class(ori),          intent(inout) :: o
        type(ori) :: o_new
        integer   :: ipeak
        call os%new( self%npeaks )
        do ipeak=1,self%npeaks
            o_new = o
            call self%get_ori(ipeak, o_new)
            call os%set_ori(ipeak, o_new)
        enddo
    end subroutine get_oris

    !>  \brief  to get one orientation
    subroutine get_ori( self, ipeak, o2update )
        class(prime3D_srch), intent(inout) :: self
        integer,             intent(in)    :: ipeak
        class(ori),          intent(inout) :: o2update
        type(ori) :: o
        real      :: euls(3), x, y, rstate, ow
        integer   :: ind
        if( str_has_substr(self%refine,'shc') .and. ipeak /= 1 ) stop 'get_ori not for shc-modes; simple_prime3D_srch'
        if( ipeak < 1 .or. ipeak > self%npeaks ) stop 'Invalid index in simple_prime3D_srch::get_ori'
        ind    = self%npeaks - ipeak + 1
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

    !>  \brief  standard deviation
    function ang_sdev( self )result( sdev )
        class(prime3D_srch), intent(inout) :: self
        real    :: sdev
        integer :: nstates, state, pop
        sdev = 0.
        if( self%npeaks < 3 .or. self%refine.eq.'shc' .or. self%refine.eq.'shcneigh' ) return
        nstates = 0
        do state=1,self%nstates
            pop = self%o_peaks%get_statepop( state )
            if( pop > 0 )then
                sdev = sdev + ang_sdev_state( state )
                nstates = nstates + 1
            endif
        enddo
        sdev = sdev / real( nstates )
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::CALCULATED ANG_SDEV'
    
        contains
            
            function ang_sdev_state( istate )result( isdev )
                use simple_stat, only: moment
                integer, intent(in)  :: istate
                type(ori)            :: o_best, o
                type(oris)           :: os
                real, allocatable    :: dists(:), ws(:)
                integer, allocatable :: inds(:)
                real                 :: ave, isdev, var
                integer              :: loc(1), alloc_stat, i, ind, n, cnt
                logical              :: err
                isdev = 0.
                inds = self%o_peaks%get_state( istate )
                n = size(inds)
                if( n < 3 )return ! because one is excluded in the next step & moment needs at least 2 objs
                call os%new( n )
                allocate( ws(n), dists(n-1), stat=alloc_stat )
                call alloc_err( 'ang_sdev_state; simple_prime3D_srch', alloc_stat )
                ws    = 0.
                dists = 0.
                ! get best ori
                do i=1,n
                    ind = inds(i)
                    ws(i) = self%o_peaks%get(ind,'ow')
                    call os%set_ori( i, self%o_peaks%get_ori(ind) )
                enddo
                loc    = maxloc( ws )
                o_best = os%get_ori( loc(1) )
                ! build distance vector
                cnt = 0
                do i=1,n
                    if( i==loc(1) )cycle
                    cnt = cnt+1
                    o = os%get_ori( i )
                    dists( cnt ) = rad2deg( o.euldist.o_best )
                enddo
                call moment(dists, ave, isdev, var, err)
                deallocate( ws, dists, inds )
                call os%kill
            end function ang_sdev_state

    end function ang_sdev
    
    ! PREPARATION ROUTINES

    !>  \brief  prepares reference indices for the search & fetches ctf
    subroutine prep4srch( self, o_prev, nnmat )
        class(prime3D_srch),        intent(inout) :: self
        class(ori),       optional, intent(inout) :: o_prev
        integer,          optional, intent(in)    :: nnmat(self%nprojs,self%nnn)
        if(str_has_substr(self%refine,'neigh') .and. .not. present(nnmat) )&
        &stop 'need optional nnmat to be present for refine=neigh modes :: prep4srch (prime3D_srch)'
        if( (str_has_substr(self%refine,'neigh') .or. str_has_substr(self%refine,'shc')).and. .not. present(o_prev) )&
        &stop 'need optional o_prev to be present refine=shc modes :: prep4srch (prime3D_srch)'
        ! Default values
        self%prev_roind  = 1
        self%prev_state  = 1
        self%prev_ref    = 1
        self%prev_proj   = 1
        self%prev_shvec  = 0.
        ! sets refence alignemnt parameters
        if( present(o_prev) )then
            ! PREVIOUS ALIGNMENT PARAMETERS
            self%prev_state = nint(o_prev%get('state'))                                      ! state index            
            self%prev_roind = self%srch_common%roind(360.-o_prev%e3get())                    ! in-plane angle index
            self%prev_shvec = o_prev%get_shift()                                             ! shift vector
            self%prev_proj  = self%e%find_closest_proj(o_prev,1)                             ! projection direction
            if( self%prev_state>self%nstates ) stop 'previous best state outside boundary; prep4srch; simple_prime3D_srch'
            if( .not.self%state_exists(self%prev_state) )stop 'empty previous state; prep4srch; simple_prime3D_srch'
            select case( self%refine )
                case( 'no', 'shc','adasym')                                                  ! DISCRETE CASE
                    call self%prep_reforis                                                   ! search space & order prep
                    self%prev_ref = self%o_refs%find_closest_proj(o_prev, self%prev_state)   ! find closest ori with same state
                case( 'neigh', 'shcneigh' )                                                  ! DISCRETE CASE WITH NEIGHBOURHOOD
                    call self%prep_reforis( nnvec=nnmat(self%prev_proj,:) )                  ! search space & order prep
                    self%prev_ref = self%o_refs%find_closest_proj(o_prev, self%prev_state)   ! find closest ori with same state
                case( 'het' )
                    self%prev_ref = (self%prev_state-1)*self%nprojs+self%prev_proj
                case DEFAULT
                    stop 'Unknown refinement mode; simple_prime3D_srch; prep4srch'
            end select
        else
            call self%prep_reforis
        endif
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::PREPARED FOR SIMPLE_PRIME3D_SRCH'
    end subroutine prep4srch

    !>  \brief  prepares the search space (ref oris) & search order per particle
    subroutine prep_reforis( self, o_prev, nnvec )
        use simple_strings,  only: int2str_pad
        use simple_ran_tabu, only: ran_tabu
        class(prime3D_srch), intent(inout) :: self
        type(ori), optional, intent(inout) :: o_prev
        integer,   optional, intent(in)    :: nnvec(self%nnn)
        type(ran_tabu)       :: irt, rt
        integer, allocatable :: isrch_order(:)
        integer              :: i, istate, n, istart, iend, start, end, nprojs
        ! on exit all the oris are clean and only the out-of-planes, 
        ! state & proj fields are present
        select case( self%refine )
            case( 'no', 'shc', 'neigh', 'shcneigh', 'adasym' )
                if( allocated(self%srch_order) ) deallocate(self%srch_order)
                if( str_has_substr(self%refine, 'neigh') )then ! local refinement
                    allocate(self%srch_order(self%nnnrefs))
                    rt = ran_tabu(self%nnnrefs)
                else
                    allocate(self%srch_order(self%nrefs))
                    rt = ran_tabu(self%nrefs)
                endif
                self%srch_order = 0
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
                if( any(self%srch_order == 0) ) stop 'Invalid index in srch_order; simple_prime3d_srch::prep_ref_oris'
                call prep_discrete_reforis
                call rt%kill
            case DEFAULT
                stop 'Unknown refinement method; prime3D_srch; prep_reforis'                
        end select

        contains

            !>  Discrete search space: multi-state spiral taken from the builder's
            subroutine prep_discrete_reforis
                type(ori)  :: o
                integer    :: cnt, istate, iproj
                call self%o_refs%new( self%nrefs )          ! init references object
                cnt = 0
                do istate=1,self%nstates
                    do iproj=1,self%nprojs
                        cnt = cnt + 1
                        o = self%e%get_ori( iproj )
                        call o%set( 'state', real(istate) ) ! Updates state
                        call o%set( 'proj', real(iproj) )   ! Updates proj
                        call self%o_refs%set_ori( cnt,o )
                    enddo
                enddo
            end subroutine prep_discrete_reforis

    end subroutine prep_reforis

    !>  \brief  prepares correlation target (previous best) for the search
    subroutine prep_corr4srch( self, pftcc, iptcl, lp, o_prev )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        real,                    intent(in)    :: lp
        class(ori), optional,    intent(inout) :: o_prev
        real :: cc_t_min_1, corr
        corr = 1.
        if( present(o_prev) )then
            corr = max( 0., pftcc%corr(self%prev_ref, iptcl, self%prev_roind) )
            if( corr > 1. .or. .not. is_a_number(corr) )then
                stop 'Invalid correlation value in simple_prime3d_srch::prep_corr4srch'
            endif
            if( (self%refine.eq.'no' .or. self%refine.eq.'adasym') .and. self%nstates==1 )then
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
        endif
        self%prev_corr = corr
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::PREPARED CORRELATION FOR SIMPLE_PRIME3D_SRCH'
    end subroutine prep_corr4srch

    !>  \brief retrieves and preps npeaks orientations for reconstruction
    subroutine prep_npeaks_oris( self )
        class(prime3D_srch), intent(inout) :: self
        type(ori) :: o, o_new
        real      :: euls(3), shvec(2)
        real      :: corr, ow, frac
        integer   :: ipeak, cnt, ref, state, proj
        integer   :: neff_states ! number of effective (non-empty) states
        ! empty states
        neff_states = 1
        if( self%nstates > 1 ) neff_states = count( self%state_exists )
        ! init
        call self%o_peaks%new( self%npeaks )
        do ipeak=1,self%npeaks
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
            ow    = o%get('ow')
            euls  = o%get_euler()
            ! add shift
            shvec = self%prev_shvec
            if( self%doshift )shvec = shvec + o%get_shift()
            if( abs(shvec(1)) < 1e-6 ) shvec(1) = 0.
            if( abs(shvec(2)) < 1e-6 ) shvec(2) = 0.
            ! copy info to new ori
            call o_new%new
            call o_new%set_euler( euls )  
            call o_new%set_shift( shvec )
            call o_new%set('state', real(state))
            call o_new%set('proj',  real(proj) )
            call o_new%set('corr',  corr       )
            call o_new%set('ow',    ow         )
            ! stashes in self
            call self%o_peaks%set_ori( ipeak, o_new )
        enddo
        ! other variables
        if( str_has_substr(self%refine, 'neigh') )then
            frac = 100.*real(self%nrefs_eval) / real(self%nnn * neff_states)
        else
            frac = 100.*real(self%nrefs_eval) / real(self%nprojs * neff_states)
        endif
        call self%o_peaks%set_all2single('frac',   frac            )
        call self%o_peaks%set_all2single('mi_hard',0.              )
        call self%o_peaks%set_all2single('dist',   0.              )
        call self%o_peaks%set_all2single('sdev',   self%ang_sdev() )
        ! ctf parameters
        if( self%ctf.ne.'no' )then
            call self%o_peaks%set_all2single('dfx',   self%dfx   )
            call self%o_peaks%set_all2single('dfy',   self%dfy   )
            call self%o_peaks%set_all2single('angast',self%angast)
        endif
        ! DEBUG
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED PREP_NPEAKS_ORIS'
    end subroutine prep_npeaks_oris

    ! SEARCH ROUTINES
    
    !>  \brief a master prime search routine
    subroutine exec_prime3D_srch( self, pftcc, iptcl, orientation, lp, greedy, nnmat )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(ori),              intent(inout) :: orientation
        real,                    intent(in)    :: lp
        logical, optional,       intent(in)    :: greedy
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn)
        logical :: ggreedy
        ggreedy = .false.
        if( present(greedy) ) ggreedy = greedy
        if( ggreedy )then
            call self%greedy_srch( pftcc,  iptcl, orientation, lp, nnmat )
        else
            call self%stochastic_srch( pftcc, iptcl, orientation, lp, nnmat )
        endif
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SRCH'
    end subroutine exec_prime3D_srch

    !>  \brief a master prime search routine
    subroutine exec_prime3D_srch_shc( self, pftcc, iptcl, orientation, lp, greedy, nnmat )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(ori),              intent(inout) :: orientation
        real,                    intent(in)    :: lp
        logical, optional,       intent(in)    :: greedy
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn)
        logical :: ggreedy
        ggreedy = .false.
        if( present(greedy) ) ggreedy = greedy
        if( ggreedy )then
            call self%greedy_srch( pftcc, iptcl, orientation, lp, nnmat )
        else
            call self%stochastic_srch_shc( pftcc, iptcl, orientation, lp, nnmat )
        endif
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SHC_SRCH'
    end subroutine exec_prime3D_srch_shc

    !>  \brief state labeler
    subroutine exec_prime3D_srch_het( self, pftcc, iptcl, orientation, extr_bound, statecnt )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(ori),              intent(inout) :: orientation
        real,                    intent(in)    :: extr_bound
        integer,                 intent(inout) :: statecnt(self%nstates)
        call self%stochastic_srch_het( pftcc, iptcl, orientation, extr_bound, statecnt )
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED PRIME3D_HET_SRCH'
    end subroutine exec_prime3D_srch_het

    !>  \brief  greedy hill-climbing (4 initialisation)
    subroutine greedy_srch( self, pftcc, iptcl, orientation, lp, nnmat )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(ori),              intent(inout) :: orientation
        real,                    intent(in)    :: lp
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn)
        real    :: projspace_corrs(self%nrefs),wcorr
        integer :: iref,isample,nrefs
        if( nint(orientation%get('state')) > 0 )then
            ! initialize
            call self%prep4srch(orientation, nnmat)
            call self%prep_corr4srch(pftcc, iptcl, lp, orientation)
            self%nbetter         = 0
            self%nrefs_eval      = 0
            self%proj_space_inds = 0
            projspace_corrs      = -1.
            nrefs = self%nrefs
            if( self%refine.eq.'neigh' ) nrefs = self%nnnrefs 
            ! search
            do isample=1,nrefs
                iref = self%srch_order(isample) ! set the stochastic reference index
                call per_ref_srch(iref)         ! actual search
            end do
            ! in greedy mode, we evaluate all refs
            self%nrefs_eval = nrefs
            ! sort in correlation projection direction space
            call hpsort(self%nrefs, projspace_corrs, self%proj_space_inds) 
            call self%inpl_srch(iptcl) ! search shifts
            ! prepare weights and orientations
            call self%prep_npeaks_oris
            call self%stochastic_weights(wcorr)
            if( self%doshift ) call self%sort_shifted_peaks
            call self%get_ori_best(orientation)
            call orientation%set('corr', wcorr)
        else
            call orientation%reject
        endif
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED GREEDY SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind, state
                state = 1
                if( self%nstates > 1 ) state = nint( self%o_refs%get(iref, 'state') )
                if( self%state_exists(state) )then
                    corrs     = pftcc%gencorrs_serial(iref, iptcl) ! In-plane correlations
                    loc       = maxloc(corrs)                      ! greedy in-plane
                    inpl_ind  = loc(1)                             ! in-plane angle index
                    inpl_corr = corrs(inpl_ind)                    ! max in plane correlation
                    call self%store_solution( iref, iref, inpl_ind, inpl_corr )
                    projspace_corrs( iref ) = inpl_corr            ! stash in-plane correlation for sorting
                endif
            end subroutine per_ref_srch

    end subroutine greedy_srch

    !>  \brief  executes the stochastic soft orientation search
    subroutine stochastic_srch( self, pftcc, iptcl, orientation, lp, nnmat )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(ori),              intent(inout) :: orientation
        real,                    intent(in)    :: lp
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn)
        real      :: projspace_corrs(self%nrefs),wcorr
        integer   :: iref,isample,nrefs
        logical   :: adasym_refine
        adasym_refine = .false.
        if( str_has_substr(self%refine,'adasym') ) adasym_refine = .true.
        ! execute search
        if( nint(orientation%get('state')) > 0 )then
            ! initialize
            call self%prep4srch(orientation, nnmat)
            call self%prep_corr4srch(pftcc, iptcl, lp, orientation)
            self%nbetter         = 0
            self%nrefs_eval      = 0
            self%proj_space_inds = 0
            projspace_corrs      = -1.
            nrefs = self%nrefs
            if( self%refine.eq.'neigh' ) nrefs = self%nnnrefs
            if( adasym_refine )then
                do isample=1,nrefs
                    iref = self%srch_order(isample)         ! set the stochastic reference index
                    if( iref == self%prev_ref ) cycle       ! previous best considered last
                    call per_ref_srch_adasym(iref)          ! actual search
                    if( self%nbetter >= self%npeaks ) exit  ! exit condition
                end do
                if( self%nbetter < self%npeaks )then
                    call per_ref_srch_adasym(self%prev_ref) ! evaluate previous best ref last
                endif
            else
                do isample=1,nrefs
                    iref = self%srch_order(isample)         ! set the stochastic reference index
                    if( iref == self%prev_ref ) cycle       ! previous best considered last
                    call per_ref_srch(iref)                 ! actual search
                    if( self%nbetter >= self%npeaks ) exit  ! exit condition
                end do
                if( self%nbetter < self%npeaks )then
                    call per_ref_srch(self%prev_ref )       ! evaluate previous best ref last
                endif
            endif
            ! sort in correlation projection direction space
            call hpsort(self%nrefs, projspace_corrs, self%proj_space_inds) 
            call self%inpl_srch(iptcl) ! search shifts
            ! prepare weights and orientations
            call self%prep_npeaks_oris
            call self%stochastic_weights(wcorr)
            if( self%doshift ) call self%sort_shifted_peaks
            call self%get_ori_best(orientation)
            call orientation%set('corr', wcorr)
        else
            call orientation%reject
        endif
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH'

        contains

            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real    :: corrs(self%nrots), inpl_corr
                integer :: loc(1), inpl_ind, state
                state = 1
                if( self%nstates > 1 ) state = nint( self%o_refs%get(iref, 'state') )
                if( self%state_exists(state) )then
                    corrs     = pftcc%gencorrs_serial(iref, iptcl) ! In-plane correlations
                    loc       = maxloc(corrs)                      ! greedy in-plane
                    inpl_ind  = loc(1)                             ! in-plane angle index
                    inpl_corr = corrs(inpl_ind)                    ! max in plane correlation
                    call self%store_solution( iref, iref, inpl_ind, inpl_corr )
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

            subroutine per_ref_srch_adasym( iref )
                integer, intent(in) :: iref
                real                :: corrs(self%nrots), inpl_corr
                integer             :: loc(1), inpl_ind, state, iref_sym, isym
                state = 1
                if( self%nstates>1 ) state = nint( self%o_refs%get(iref, 'state') )
                if( self%state_exists( state ) )then
                    corrs     = pftcc%gencorrs_serial(iref, iptcl) ! In-plane correlations
                    loc       = maxloc(corrs)                      ! greedy in-plane
                    inpl_ind  = loc(1)                             ! in-plane angle index
                    inpl_corr = corrs(loc(1))                      ! max in plane correlation
                    call self%store_solution( iref, iref, inpl_ind, inpl_corr )
                    projspace_corrs( iref ) = inpl_corr            ! stash in-plane correlation for sorting
                    ! update nbetter to keep track of how many improving solutions we have identified
                    if( self%npeaks == 1 )then
                        if( inpl_corr > self%prev_corr ) self%nbetter = self%nbetter + 1
                    else
                        if( inpl_corr >= self%prev_corr ) self%nbetter = self%nbetter + 1
                    endif
                    ! keep track of how many references we are evaluating
                    self%nrefs_eval = self%nrefs_eval + 1
                    if( inpl_corr >= self%prev_corr )then
                        ! evaluate symmetry-related orientations
                        do isym=1,self%nsym - 1
                            iref_sym = self%nnmat_sym(iref,isym)
                            if( self%proj_space_inds(iref_sym) /= 0 ) cycle
                            corrs     = pftcc%gencorrs_serial(iref_sym, iptcl) ! In-plane correlations
                            loc       = maxloc(corrs)                          ! greedy in-plane
                            inpl_ind  = loc(1)                                 ! in-plane angle index
                            inpl_corr = corrs(loc(1))                          ! max in plane correlation
                            call self%store_solution( iref_sym, iref_sym, inpl_ind, inpl_corr )
                            projspace_corrs( iref_sym ) = inpl_corr            ! stash in-plane correlation for sorting
                        end do
                    endif
                endif
            end subroutine per_ref_srch_adasym

    end subroutine stochastic_srch

    !>  \brief  executes the stochastic hard orientation search using pure stochastic hill climbing
    !!          (no probabilistic weighting + stochastic search of in-plane angles)
    subroutine stochastic_srch_shc( self, pftcc, iptcl, orientation, lp, nnmat )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(ori),              intent(inout) :: orientation
        real,                    intent(in)    :: lp
        integer, optional,       intent(in)    :: nnmat(self%nprojs,self%nnn)
        real      :: inpl_corr,corrs(self%nrots)
        integer   :: iref,isample,nrefs,inpl_ind,loc(1)
        logical   :: found_better
        if( nint(orientation%get('state')) > 0 )then
            ! initialize
            call self%prep4srch( orientation, nnmat )
            call self%prep_corr4srch( pftcc, iptcl, lp, orientation )
            found_better         = .false.
            self%nrefs_eval      = 0
            self%proj_space_inds = 0
            nrefs                = self%nrefs
            if( self%refine.eq.'shcneigh' ) nrefs = self%nnnrefs               
            ! search
            do isample=1,nrefs
                iref = self%srch_order( isample ) ! stochastic reference index
                if( iref == self%prev_ref ) cycle ! previous best considered last
                call per_ref_srch( iref )         ! actual search
                if( inpl_ind > 0 )then
                    ! correlation-improving reference found, store it
                    call self%store_solution( self%nrefs, iref, inpl_ind, inpl_corr )
                    found_better = .true.        ! flag for found solution
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
            call self%inpl_srch(iptcl)
            ! output
            call self%prep_npeaks_oris
            call self%get_ori_best(orientation)
        else
            call orientation%reject
        endif
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH (REFINE=SHC|SHCNEIGH)'

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
                    corrs    = pftcc%gencorrs_serial( iref, iptcl )      ! in-plane correlations 
                    inpl_ind = shcloc(self%nrots, corrs, self%prev_corr) ! first improving in-plane index
                    if( inpl_ind > 0 ) inpl_corr = corrs( inpl_ind )     ! improving correlation found
                    self%nrefs_eval = self%nrefs_eval + 1                ! updates fractional search space
                endif
            end subroutine per_ref_srch

    end subroutine stochastic_srch_shc
    
    subroutine stochastic_srch_het( self, pftcc, iptcl, orientation, extr_bound, statecnt )
        use simple_rnd, only: shcloc, irnd_uni
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        class(ori),              intent(inout) :: orientation
        real,                    intent(in)    :: extr_bound
        integer,                 intent(inout) :: statecnt(self%nstates)
        integer :: istate,iref,state,refs(self%nstates),loc(1)
        real    :: corr,mi_state,frac,corrs(self%nstates)
        if( nint(orientation%get('state')) > 0 )then
            ! initialize
            call self%prep4srch(orientation)
            if( orientation%get('corr') < extr_bound )then
                ! state randomization
                statecnt(self%prev_state) = statecnt(self%prev_state) + 1
                self%nrefs_eval = 1
                state = irnd_uni(self%nstates)
                do while( state == self%prev_state )
                    state = irnd_uni(self%nstates)
                enddo
                iref = (state - 1) * self%nprojs + self%prev_proj
                corr = pftcc%corr(iref, iptcl, self%prev_roind)
            else
                ! SHC
                refs = 0
                do istate=1,self%nstates
                    refs(istate) =  (istate-1)*self%nprojs + self%prev_proj 
                    if( .not. self%state_exists(istate) ) refs(istate) = 0
                enddo
                corrs           = pftcc%corrs(refs, self%nstates, iptcl, self%prev_roind)
                self%prev_corr  = corrs(self%prev_state)
                loc             = shcloc(self%nstates, corrs, self%prev_corr)
                state           = loc(1)
                corr            = corrs(state)
                self%nrefs_eval = count(corrs <= self%prev_corr)
            endif
            frac = 100.*real(self%nrefs_eval)/real(self%nstates)
            call orientation%set('frac', frac)
            call orientation%set('state', real(state))
            call orientation%set('corr', corr)
            call orientation%set('mi_proj', 1.)
            call orientation%set('mi_inpl',  1.)
            if( self%prev_state .ne. state )then
                mi_state = 0.
            else
                mi_state = 1.
            endif
            call orientation%set('mi_state', mi_state)
            call orientation%set('mi_joint', mi_state)
            call orientation%set('w', 1.)
            call self%o_peaks%set_ori(1,orientation)
            call self%prep_npeaks_oris
            call self%get_ori_best(orientation)
        else
            call orientation%reject
        endif
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED HET SEARCH'
    end subroutine stochastic_srch_het
    
    !>  \brief  executes the in-plane search for discrete mode
    subroutine inpl_srch( self, iptcl )
        class(prime3D_srch), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        real, allocatable :: cxy(:), crxy(:)
        type(ori) :: o
        real      :: cc
        integer   :: i, ref, inpl_ind
        if( self%doshift )then
            do i=self%nrefs,self%nrefs-self%npeaks+1,-1
                ref      = self%proj_space_inds( i )
                o        = self%o_refs%get_ori( ref )
                cc       = o%get('corr')
                inpl_ind = self%srch_common%roind( 360.-o%e3get() )
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
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED INPL SEARCH'
    end subroutine inpl_srch

    subroutine gen_symnnmat( self )
        class(prime3D_srch), intent(inout) :: self
        type(ori) :: oref,osym
        integer :: i, isym,iref
        if( allocated(self%nnmat_sym) )deallocate(self%nnmat_sym)
        allocate( self%nnmat_sym(self%nrefs,self%nsym-1) )
        do iref = 1,self%nrefs
            oref = self%e%get_ori( iref )
            do isym=2,self%nsym
                osym = self%se%apply(oref, isym)
                self%nnmat_sym(iref,isym-1) = self%e%find_closest_proj( osym )
            enddo
        enddo
    end subroutine gen_symnnmat

    subroutine stochastic_weights( self, wcorr )
        class(prime3D_srch), intent(inout) :: self
        real,                intent(out)   :: wcorr
        real             :: ws(self%npeaks)
        real,allocatable :: corrs(:)
        integer          :: i
        if( self%npeaks == 1 )then
            call self%o_peaks%set(1,'ow',1.0)
            return
        endif
        if( self%o_peaks%get_noris() /= self%npeaks ) stop 'invalid number of references in simple_prime3D_srch::stochastic_weights'
        ws    = 0.
        wcorr = 0.
        ! get unnormalised correlations
        corrs = self%o_peaks%get_all('corr')
        ! calculate normalised weights and weighted corr
        where( corrs > TINY ) ws = exp(corrs) ! ignore invalid corrs
        ws    = ws/sum(ws)
        wcorr = sum(ws*corrs) 
        ! update npeaks individual weights
        do i=1,self%npeaks
            call self%o_peaks%set(i,'ow',ws(i))
        enddo
        deallocate(corrs)
        if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH :: CALCULATED STOCHASTIC WEIGHTS PRIME3D_SRCH'        
    end subroutine stochastic_weights

    subroutine store_solution( self, ind, ref, inpl_ind, corr )
        class(prime3D_srch), intent(inout) :: self
        integer,             intent(in)    :: ind, ref, inpl_ind
        real,                intent(in)    :: corr
        real :: e3
        self%proj_space_inds( ind ) = ref          ! stash reference similarly to other search modes
        e3 = 360.-self%srch_common%rot( inpl_ind ) ! psi
        call self%o_refs%set( ref, 'e3',   e3   )  ! stash psi
        call self%o_refs%set( ref, 'corr', corr )  ! stash correlation
        if( self%npeaks==1 ) call self%o_refs%set( ref, 'ow', 1. ) ! set reconstruction weight
    end subroutine store_solution

    ! GETTERS FOR TESTING

    !>  \brief  is for getting the search order
    function get_srch_order( self )result( inds )
        class(prime3D_srch), intent(inout) :: self
        integer,allocatable :: inds(:)
        integer             :: alloc_stat
        allocate( inds(size(self%srch_order)), stat=alloc_stat )
        call alloc_err( 'simple_prime3D_srch::get_srch_order', alloc_stat)
        inds(:) = self%srch_order(:)
    end function get_srch_order

    !>  \brief  is for getting o_peaks
    function get_o_peaks( self, iptcl )result( out_os )
        class(prime3D_srch), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        type(oris) :: out_os
        out_os = self%o_peaks
    end function get_o_peaks

    !>  \brief  is for getting the n-last reference orientations
    function get_o_refs( self, n )result( out_os )
        class(prime3D_srch), intent(inout) :: self
        integer, optional,   intent(in)    :: n
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
            nrefs = self%nnnrefs
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

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(prime3D_srch), intent(inout) :: self !< instance
        integer :: i, istart, istop
        if( self%exists )then
            call self%e%kill
            call self%srch_common%kill
            call self%o_refs%kill
            call self%o_peaks%kill 
            call self%se%kill
            if( allocated( self%proj_space_inds) ) deallocate( self%proj_space_inds )
            if( allocated( self%srch_order )     ) deallocate( self%srch_order)
            if( allocated( self%nnmat_sym)       ) deallocate(self%nnmat_sym)
            if( allocated( self%state_exists )   ) deallocate(self%state_exists)
            self%exists = .false.
        endif        
    end subroutine kill

    !>  \brief  executes the stochastic soft orientation search
    ! subroutine stochastic_srch( self, pftcc, iptcl )
    !     class(prime3D_srch),     intent(inout) :: self
    !     class(polarft_corrcalc), intent(inout) :: pftcc
    !     integer,                 intent(in)    :: iptcl
    !     real    :: projspace_corrs( self%nrefs )
    !     integer :: ref, isample, nrefs
    !     ! initialize
    !     self%nbetter         = 0
    !     self%nrefs_eval      = 0
    !     self%proj_space_inds = 0
    !     projspace_corrs      = -1.
    !     ! search
    !     nrefs = self%nrefs
    !     if( self%refine.eq.'neigh' )nrefs=self%nnnrefs
    !     do isample=1,nrefs
    !         ref = self%srch_order( isample )                             ! set the stochastic reference index
    !         if( ref==self%prev_ref )cycle                                ! previous best considered last
    !         if( ref < 1 .or. ref > self%nrefs ) stop 'ref index out of bound; simple_prime3D_srch::stochastic_srch'
    !         call per_ref_srch( ref )                                     ! actual search
    !         if( self%nbetter >= self%npeaks ) exit                       ! exit condition
    !     end do
    !     if( self%nbetter<self%npeaks )call per_ref_srch( self%prev_ref ) ! evaluate previous best ref last
    !     call hpsort(self%nrefs, projspace_corrs, self%proj_space_inds)   ! sort in correlation projection direction space
    !     call self%inpl_srch( iptcl )                                     ! search shifts
    !     if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH'

    !     contains

    !         subroutine per_ref_srch( iref )
    !             integer, intent(in) :: iref
    !             real                :: corrs(self%nrots), inpl_corr
    !             integer             :: loc(1), inpl_ind, state
    !             state = 1
    !             if( self%nstates>1 )state = nint( self%o_refs%get(iref, 'state') )
    !             if( self%state_exists( state ) )then
    !                 corrs     = pftcc%gencorrs(iref, iptcl) ! In-plane correlations
    !                 loc       = maxloc(corrs)               ! greedy in-plane
    !                 inpl_ind  = loc(1)                      ! in-plane angle index
    !                 inpl_corr = corrs(inpl_ind)             ! max in plane correlation
    !                 call self%store_solution( iref, iref, inpl_ind, inpl_corr )
    !                 projspace_corrs( iref ) = inpl_corr     ! stash in-plane correlation for sorting
    !                 ! update nbetter to keep track of how many improving solutions we have identified
    !                 if( self%npeaks == 1 )then
    !                     if( inpl_corr > self%prev_corr ) self%nbetter = self%nbetter+1
    !                 else
    !                     if( inpl_corr >= self%prev_corr ) self%nbetter = self%nbetter+1
    !                 endif
    !                 ! keep track of how many references we are evaluating
    !                 self%nrefs_eval = self%nrefs_eval+1
    !             endif
    !         end subroutine per_ref_srch

    ! end subroutine stochastic_srch

    !>  \brief a master prime search routine
    ! subroutine exec_prime3D_srch( self, pftcc, iptcl, lp, o, nnmat )
    !     class(prime3D_srch),     intent(inout) :: self
    !     class(polarft_corrcalc), intent(inout) :: pftcc
    !     integer,                 intent(in)    :: iptcl
    !     real,                    intent(in)    :: lp
    !     class(ori), optional,    intent(inout) :: o
    !     integer,    optional,    intent(in)    :: nnmat(self%nprojs,self%nnn)
    !     real :: wcorr
    !     call self%prep4srch( o, nnmat )
    !     call self%prep_corr4srch( pftcc, iptcl, lp, o )
    !     call self%prep_inpl_srch( pftcc )
    !     if( str_has_substr(self%refine,'adasym') )then
    !         call self%stochastic_srch_adasym( pftcc, iptcl )
    !     else
    !         call self%stochastic_srch( pftcc, iptcl )
    !     endif
    !     call self%prep_npeaks_oris
    !     call self%stochastic_weights( wcorr, self%o_peaks )
    !     if( present(o)   ) call o%set('corr', wcorr)
    !     if( self%doshift ) call self%sort_shifted_npeaks( self%o_peaks )
    !     if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SRCH'
    ! end subroutine exec_prime3D_srch

    !>  \brief a master shift/in-plane prime search routine
    ! subroutine exec_prime3D_inpl_srch( self, pftcc, iptcl, lp, o, greedy )
    !     class(prime3D_srch),     intent(inout) :: self
    !     class(polarft_corrcalc), intent(inout) :: pftcc
    !     integer,                 intent(in)    :: iptcl
    !     real,                    intent(in)    :: lp
    !     class(ori),              intent(inout) :: o
    !     logical,    optional,    intent(in)    :: greedy
    !     logical :: doshift, dogreedy, greedy_inpl, found_better
    !     if( .not.self%doshift )stop 'In-plane search not activated yet; simple_prime3D_srch%exec_prime3D_inpl_srch'
    !     if( self%npeaks>1 )stop 'In-plane only permitted with npeaks=1; simple_prime3D_srch%exec_prime3D_inpl_srch'
    !     ! Whether to perform shift or inpl search
    !     greedy_inpl = self%greedy_inpl
    !     dogreedy = .false.
    !     if( present(greedy) )dogreedy = greedy
    !     self%greedy_inpl = dogreedy
    !     ! Search
    !     call self%prep4srch( o )
    !     call self%prep_corr4srch( pftcc, iptcl, lp, o )
    !     call self%prep_inpl_srch( pftcc )
    !     call self%stochastic_srch_inpl( pftcc, iptcl, o )
    !     call self%prep_npeaks_oris(iptcl)
    !     call o%set('corr', self%os_npeaks(iptcl)%get( self%npeaks,'corr'))
    !     self%greedy_inpl = greedy_inpl ! restores setting
    !     if( DEBUG ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED PRIME3D_INPL_SRCH'
    ! end subroutine exec_prime3D_inpl_srch

    ! subroutine stochastic_srch_inpl( self, pftcc, iptcl, o )
    !     use simple_rnd,      only: ran3
    !     use simple_ran_tabu, only: ran_tabu
    !     class(prime3D_srch),     intent(inout) :: self
    !     class(polarft_corrcalc), intent(inout) :: pftcc
    !     class(ori),              intent(inout) :: o
    !     integer,                 intent(in)    :: iptcl
    !     type(ori)      :: o_tmp
    !     type(ran_tabu) :: rt
    !     real           :: target_corr, corr, e3, delta, p, rnd
    !     integer        :: i, istate, ref, nrefs, isrch_order(self%nstates), loc(1)
    !     logical        :: found_better
    !     ! init
    !     nrefs = self%nrefs
    !     if( str_has_substr(self%refine,'neigh') )nrefs = self%nnn
    !     found_better         = .false.
    !     self%nrefs_eval      = 0
    !     self%proj_space_inds = 0
    !     !e3 = 360.-self%srch_common%rot( self%prev_roind )       ! psi
    !     if( self%nstates==1 )then
    !         ! single state
    !         call self%store_solution( self%nrefs, self%prev_ref, self%prev_roind, self%prev_corr ) 
    !         self%nrefs_eval = nrefs                              ! updates frac
    !     else
    !         ! multi-state, first improving
    !         rt = ran_tabu( self%nstates )
    !         call rt%ne_ran_iarr( isrch_order )
    !         call rt%kill
    !         do i=1,self%nstates
    !             istate = isrch_order(i)                          ! randomized state search order
    !             if( istate==self%prev_state )cycle               ! previous state last
    !             if( .not.self%state_exists(istate) )cycle        ! empty state
    !             self%nrefs_eval = self%nrefs_eval+1              ! updates frac
    !             ref  = (istate-1)*self%nprojs + self%prev_proj   ! reference index to state projection direction
    !             corr = pftcc%corr( ref, iptcl, self%prev_roind ) ! state correlation
    !             if( corr >= self%prev_corr )then
    !                 found_better = .true.                        ! first improvement exit condition
    !                 exit
    !             endif
    !         enddo
    !         if( .not. found_better )then
    !             ref  = self%prev_ref                             ! defaults to previous state
    !             corr = self%prev_corr
    !         endif
    !         call self%store_solution( self%nrefs, ref, self%prev_roind, corr )       
    !         self%nrefs_eval = nint( real(nrefs*self%nrefs_eval)/real(self%nstates) ) ! updates frac
    !     endif
    !     call self%inpl_srch( iptcl )                             ! performs in-plane search
    ! end subroutine stochastic_srch_inpl

end module simple_prime3D_srch
