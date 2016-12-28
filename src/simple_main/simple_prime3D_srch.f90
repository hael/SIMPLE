module simple_prime3D_srch
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_prime_srch,       only: prime_srch
use simple_ran_tabu,         only: ran_tabu
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_params,           only: params
use simple_pftcc_shsrch      ! singleton
use simple_pftcc_inplsrch    ! singleton
use simple_defs              ! singleton
use simple_jiffys            ! singleton
use simple_math              ! singleton
implicit none

public :: prime3D_srch
private

logical, parameter :: debug = .false.

type prime3D_srch
    private
    class(params), pointer :: pp => null()            !< pointer to parameters singleton
    class(oris),   pointer :: pe => null()            !< pointer to orientations of references in proj_space
    type(prime_srch)       :: srch_common             !< functionalities common to primesrch2D/3D
    type(oris)             :: o_refs                  !< projection directions search space
    type(oris)             :: o_npeaks                !< orientations of best npeaks oris
    type(ran_tabu)         :: rt                      !< object for random number generation
    integer                :: nrefs          = 0      !< total number of references (nstates*nprojs)
    integer                :: nnnrefs        = 0      !< total number of neighboring references (nstates*nnn)
    integer                :: nstates        = 0      !< number of states
    integer                :: nprojs         = 0      !< number of projections (same as number of oris in o_refs)
    integer                :: nrots          = 0      !< number of in-plane rotations in polar representation
    integer                :: npeaks         = 0      !< number of peaks (nonzero orientation weights)
    integer                :: nbetter        = 0      !< nr of better orientations identified
    integer                :: nrefs_eval     = 0      !< nr of references evaluated
    integer                :: nnn            = 0      !< number of nearest neighbors
    integer                :: nsym           = 0      !< symmetry order
    integer                :: prev_roind     = 0      !< previous in-plane rotation index
    integer                :: prev_state     = 0      !< previous state index
    integer                :: prev_ref       = 0      !< previous reference index
    integer                :: prev_proj      = 0      !< previous projection index
    real                   :: prev_corr      = 1.     !< previous best correlation
    real                   :: prev_shvec(2)  = 0.     !< previous origin shift vector
    real                   :: lims(2,2)      = 0.     !< shift search range limit
    real                   :: eullims(3,2)   = 0.     !< asymetric unit bounds
    real                   :: athres         = 0.     !< angular threshold
    real                   :: dfx            = 0.     !< ctf x-defocus
    real                   :: dfy            = 0.     !< ctf y-defocus
    real                   :: angast         = 0.     !< ctf astigmatism
    integer, allocatable   :: proj_space_inds(:)      !< projection space index array
    integer, allocatable   :: srch_order(:)           !< stochastic search order
    real,    allocatable   :: prev_corrs(:)           !< previous particle-reference correlations for GPU-based search
    character(len=STDLEN)  :: refine        = ''      !< refinement flag
    character(len=STDLEN)  :: opt           = ''      !< optimiser flag
    character(len=STDLEN)  :: ctf           = ''      !< ctf flag
    character(len=STDLEN)  :: shbarr        = 'yes'   !< shift barrier flag
    logical                :: use_cpu       = .true.  !< indicates if CPU logics is being used
    logical                :: bench_gpu     = .false. !< indicates if GPU logics should be tested
    logical                :: use_gpu       = .false. !< indicates if GPU logics should be used
    logical                :: doshift       = .true.  !< origin shift search indicator
    logical                :: greedy_inpl   = .true.  !< indicator for whether in-plane search is greedy or not
    logical                :: exists        = .false. !< 2 indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! CALCULATORS
    procedure, private :: calc_qcont_corr
    procedure          :: prep_npeaks_oris
    procedure          :: sort_shifted_npeaks
    ! GETTERS
    procedure          :: get_ori_best
    procedure          :: get_ori
    procedure          :: ang_sdev
    ! PREPARATION ROUTINES
    procedure          :: prep4srch
    procedure          :: prep_corr4srch
    procedure          :: prep_reforis
    procedure, private :: prep_inpl_srch
    procedure          :: prep_ctfparms
    procedure          :: calc_corrs
    ! SEARCH ROUTINES
    procedure          :: exec_prime3D_srch
    procedure          :: exec_prime3D_qcont_srch
    procedure          :: exec_prime3D_shc_srch
    procedure          :: exec_prime3D_inpl_srch
    procedure, private :: stochastic_srch
    procedure, private :: stochastic_srch_qcont
    procedure, private :: stochastic_srch_shc
    procedure, private :: stochastic_srch_inpl
    procedure          :: inpl_srch
    procedure          :: inpl_srch_qcont
    procedure, private :: greedy_srch
    procedure          :: stochastic_weights
    ! GETTERS & SETTERS
    procedure          :: store_solution
    procedure          :: get_o_npeaks
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
    procedure :: kill
end type prime3D_srch

contains

    !>  \brief  is a constructor
    subroutine new( self, a, e, p )
        use simple_sym,   only: sym
        class(prime3D_srch),   intent(inout) :: self  !< instance
        class(params), target, intent(in)    :: p     !< parameters
        class(oris),   target, intent(in)    :: e     !< reference oris
        class(oris),           intent(inout) :: a     !< ptcls oris
        type(sym) :: symobj
        integer   :: alloc_stat
        ! destroy possibly pre-existing instance
        call self%kill
        ! set constants
        self%pp            => p 
        self%pe            => e
        self%nstates       =  p%nstates
        self%nprojs        =  p%nspace
        self%nrefs         =  self%nprojs*self%nstates
        self%nrots         =  round2even(twopi*real(p%ring2))
        self%srch_common   =  prime_srch(p, self%nrefs, self%nrots)
        self%npeaks        =  p%npeaks
        self%nbetter       =  0
        self%nrefs_eval    =  0
        self%athres        =  p%athres
        self%doshift       =  p%doshift
        self%refine        =  p%refine
        self%opt           =  p%opt
        self%ctf           =  p%ctf
        self%nnn           =  p%nnn
        self%nnnrefs       =  self%nnn*self%nstates
        self%shbarr        =  p%shbarrier
        self%bench_gpu     =  .false.
        self%use_gpu       =  .false.
        self%use_cpu       =  .true.
        if( p%bench_gpu .eq. 'yes' .or. p%use_gpu .eq. 'yes' )then
            if( p%bench_gpu .eq. 'yes' ) self%bench_gpu = .true.
            if( p%use_gpu   .eq. 'yes' ) self%use_gpu   = .true.
            if( p%top-p%fromp+1 /= self%nrefs )&
            &stop 'the particle chunk is not correctly balanced for GPU execution!'
            self%use_cpu = .false.
        endif
        if( str_has_substr(self%refine,'neigh') )then
            if( self%bench_gpu .or. self%use_gpu ) stop 'Neigh modes not implemented for GPU; simple_prime3D :: stochastic_srch'
        endif
        self%greedy_inpl = .true.
        if( str_has_substr(self%refine,'shc') )then
            if( self%npeaks > 1 ) stop 'npeaks must be equal to 1 with refine=shc|shcneigh'
            self%greedy_inpl = .false.
        endif
        ! construct composites
        allocate(self%proj_space_inds(self%nrefs), stat=alloc_stat)
        call alloc_err('In: new; simple_prime3D_srch, 1', alloc_stat)
        self%proj_space_inds = 0
        if( str_has_substr(self%refine, 'neigh') )then ! local refinement
            allocate(self%srch_order(self%nnnrefs), stat=alloc_stat)
            call alloc_err('In: new; simple_prime3D_srch, 2', alloc_stat)
            self%srch_order = 0
            self%rt = ran_tabu( self%nnnrefs )
        else
            allocate(self%srch_order(self%nrefs), stat=alloc_stat)
            call alloc_err('In: new; simple_prime3D_srch, 3', alloc_stat)
            self%srch_order = 0
            self%rt = ran_tabu(self%nrefs)
        endif
        ! generate oris oject in which the best npeaks refs will be stored
        call self%o_npeaks%new( self%npeaks )
        ! symmetry-related variables
        call symobj%new( self%pp%pgrp )
        self%nsym    = symobj%get_nsym()
        self%eullims = symobj%srchrange()
        call symobj%kill
        ! Updates option to search shift
        self%doshift = self%pp%doshift
        if( self%doshift )then
            self%lims(:,1) = -self%pp%trs
            self%lims(:,2) =  self%pp%trs
        endif
        ! the end
        self%exists = .true.
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::CONSTRUCTED NEW SIMPLE_PRIME3D_SRCH OBJECT'
    end subroutine new
    
    ! CALCULATORS

    !>  \brief  calculates correlation to previous volume in previous orientation
    function calc_qcont_corr( self, iptcl, o, pftcc, proj, tfun, refvols )result( corr )
        use simple_image,     only: image
        use simple_projector, only: projector
        use simple_ctf,       only: ctf
        class(prime3D_srch),     intent(inout) :: self
        class(image),            intent(inout) :: refvols( self%nstates )
        class(projector),        intent(inout) :: proj
        class(polarft_corrcalc), intent(inout) :: pftcc
        class(ori),              intent(inout) :: o
        class(ctf),              intent(inout) :: tfun
        integer,                 intent(in)    :: iptcl
        real    :: corr
        integer :: state, inpl_ind
        state    = nint( o%get('state') )
        inpl_ind = self%srch_common%roind( 360. )
        call proj%fproject_polar( 1, refvols(state), o, self%pp, pftcc ) ! on-the-fly polar projection
        if( self%ctf .ne. 'no' )then
            call pftcc%apply_ctf(self%pp%smpd, tfun, self%dfx, self%dfy, self%angast)
        endif
        corr = pftcc%corr( 1, iptcl, inpl_ind )
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED CALC_QCONT_CORR'
    end function calc_qcont_corr

    !> \brief  for sorting already shifted npeaks oris
    subroutine sort_shifted_npeaks( self, os )
        class(prime3D_srch), intent(inout) :: self
        type(oris),          intent(inout) :: os
        type(oris):: o_sorted
        type(ori) :: o
        real, allocatable :: corrs(:)
        integer :: i, inds(self%npeaks)
        if( os%get_noris() /= self%npeaks ) stop 'invalid number of oris in simple_prime3D_srch::sort_shifted_npeaks'
        call o_sorted%new( self%npeaks )
        corrs = os%get_all('corr')
        inds  = (/(i,i=1,self%npeaks)/)
        call hpsort(self%npeaks, corrs, inds)
        do i=1,self%npeaks
            o = os%get_ori( inds(i) )
            call o_sorted%set_ori( i,o )
        enddo
        os = o_sorted
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED SORT_SHIFTED_NPEAKS'
    end subroutine sort_shifted_npeaks

    ! GETTERS & SETTERS
    
    !>  \brief  to get the best orientation from the discrete space
    subroutine get_ori_best( self, o_prev )
        use simple_math, only: myacos, rotmat2d, rad2deg
        class(prime3D_srch), intent(inout) :: self
        class(ori),          intent(inout) :: o_prev
        type(ori) :: o
        real      :: euldist, mi_joint, mi_class, mi_inpl, mi_state
        real      :: mat(2,2), u(2), x1(2), x2(2)
        integer   :: roind,state
        ! make unit vector
        u(1) = 0.
        u(2) = 1.
        ! calculate previous vec
        mat  = rotmat2d(o_prev%e3get())
        x1   = matmul(u,mat)
        ! get new orientation
        o = self%o_npeaks%get_ori( self%npeaks )
        ! calculate new vec
        mat     = rotmat2d(o%e3get())
        x2      = matmul(u,mat)
        euldist = rad2deg( o_prev.euldist.o )
        roind   = self%srch_common%roind( 360.-o%e3get() )
        state   = nint( o%get('state') )
        ! calculate overlap between distributions
        mi_class = 0.
        mi_inpl  = 0.
        mi_state = 0.
        mi_joint = 0.
        if( euldist < 0.1 )then
            mi_class = mi_class + 1.
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
        call o_prev%set('mi_class', mi_class)
        call o_prev%set('mi_inpl',  mi_inpl)
        call o_prev%set('mi_state', mi_state)
        call o_prev%set('mi_joint', mi_joint)
        ! set the distances before we update the orientation
        call o_prev%set('dist', 0.5*euldist + 0.5*o_prev%get('dist'))
        call o_prev%set('dist_inpl', rad2deg(myacos(dot_product(x1,x2))))
        ! all the other stuff
        call o_prev%set_euler( o%get_euler() )
        call o_prev%set_shift( o%get_shift() )
        call o_prev%set( 'state',   real(state) )
        call o_prev%set( 'frac',    o%get('frac') )
        call o_prev%set( 'corr',    o%get('corr') )
        call o_prev%set( 'ow',      o%get('ow') )
        call o_prev%set( 'mirr',    0. )
        call o_prev%set( 'class',   o%get('class') )
        call o_prev%set( 'sdev',    o%get('sdev') )
        ! stash and return
        call self%o_npeaks%set_ori( self%npeaks,o_prev )
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::GOT BEST ORI'
    end subroutine get_ori_best

    !>  \brief  to get one orientation from the discrete space
    subroutine get_ori( self, i, o )
        class(prime3D_srch), intent(inout) :: self
        type(ori),           intent(inout) :: o
        integer,             intent(in)    :: i
        if( str_has_substr(self%refine,'shc') .and. i /= 1 )stop 'get_ori not for shc-modes; simple_prime3D_srch'
        if( i < 1 .or. i > self%npeaks )stop 'Invalid index in simple_prime3D_srch::get_ori'
        call o%copy( self%o_npeaks%get_ori( self%npeaks-i+1 ) )
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
            pop = self%o_npeaks%get_statepop( state )
            if( pop > 0 )then
                sdev = sdev + ang_sdev_state( state )
                nstates = nstates + 1
            endif
        enddo
        sdev = sdev / real( nstates )
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::CALCULATED ANG_SDEV'
    
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
                inds = self%o_npeaks%get_state( istate )
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
                    ws(i) = self%o_npeaks%get( ind,'ow')
                    call os%set_ori( i, self%o_npeaks%get_ori(ind) )
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

    !>  \brief  prepares reference indices for the search on CPU & fetches ctf
    subroutine prep4srch( self, o_prev, nnmat )
        class(prime3D_srch),        intent(inout) :: self
        class(ori),       optional, intent(inout) :: o_prev
        integer,          optional, intent(in)    :: nnmat(self%nprojs,self%nnn)
        if( (self%refine.eq.'neigh' .or. self%refine.eq.'shcneigh') .and. .not.present(nnmat) )&
        & stop 'o_prev and nnmat must be provided with refine=shc|shcneigh'
        ! Default values
        self%prev_roind  = 1
        self%prev_state  = 1
        self%prev_ref    = 1
        self%prev_proj   = 1
        self%prev_shvec  = 0.
        ! sets refence alignemnt parameters
        if( present(o_prev) )then
            ! PREVIOUS ALIGNMENT PARAMETERS
            self%prev_state = nint(o_prev%get('state'))                                 ! state index
            self%prev_roind = self%srch_common%roind( 360.-o_prev%e3get() )             ! in-plane angle index
            self%prev_shvec = o_prev%get_shift()                                        ! shift vector
            self%prev_proj  = self%pe%find_closest_proj( o_prev, 1 )                    ! projection direction
            if( self%prev_state>self%nstates ) stop 'previous best state outside boundary; prep4srch; simple_prime3D_srch'
            select case( self%refine )
                case( 'no', 'shc', 'shift')                                             ! DISCRETE CASE
                    call self%prep_reforis                                              ! search space & order prep
                case( 'neigh', 'shcneigh' )                                             ! DISCRETE CASE WITH NEIGHBOURHOOD
                    call self%prep_reforis( nnvec=nnmat(self%prev_proj,:) )             ! search space & order prep
                case( 'qcont', 'qcontneigh' )                                           ! QUASI-CONTINUOUS CASE
                    call self%prep_reforis( o_prev=o_prev )                             ! search space & order prep
                    !call self%o_refs%set_ori( self%prev_ref, o_prev )  ! necessary???  ! replaces closest ref with previous best
                case DEFAULT
                    stop 'Unknown refinement mode; simple_prime3D_srch; prep4srch'
            end select
            self%prev_ref = self%o_refs%find_closest_proj( o_prev, self%prev_state )    ! find closest ori with same state
         else
            call self%prep_reforis
         endif
         if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::PREPARED FOR SIMPLE_PRIME3D_SRCH'
    end subroutine prep4srch

    !>  \brief  prepares the search space (ref oris) & search order per particle
    subroutine prep_reforis( self, o_prev, nnvec )
        use simple_jiffys, only: int2str_pad
        class(prime3D_srch), intent(inout) :: self
        type(ori), optional, intent(inout) :: o_prev
        integer,   optional, intent(in)    :: nnvec(self%nnn)
        type(ran_tabu)       :: irt
        integer, allocatable :: isrch_order(:)
        integer              :: i, istate, n, istart, iend, start, end, nprojs
        ! on exit all the oris are clean and only the out-of-planes, 
        ! state & class fields are present
        select case( self%refine )
            case( 'no', 'shc', 'neigh', 'shcneigh' )
                ! DISCRETE
                ! Search order
                call self%rt%reset
                self%srch_order = 0
                if( present(nnvec) )then
                    ! refine=neigh|shcneigh
                    do istate=0,self%nstates-1 ! concatenate nearest neighbor per state...
                        i = istate*self%nnn+1
                        self%srch_order(i:i+self%nnn-1) = nnvec + istate*self%nprojs
                    enddo
                    call self%rt%shuffle( self%srch_order ) ! ...& wizz it up
                else
                    ! refine=no|shc
                    call self%rt%ne_ran_iarr( self%srch_order )
                endif
                ! PUT PREVIOUS STATE LAST - UNTESTED
                ! if( self%nstates==1 )then
                !     ! Single state
                !     if( present(nnvec) )then
                !         ! refine=neigh|shcneigh
                !         do istate=0,self%nstates-1                  ! concatenate nearest neighbor per state...
                !             i = istate*self%nnn+1
                !             self%srch_order(i:i+self%nnn-1) = nnvec + istate*self%nprojs
                !         enddo
                !         call self%rt%shuffle( self%srch_order )     ! ...& wizz it up                        
                !     else
                !         ! refine=no|shc
                !         call self%rt%ne_ran_iarr( self%srch_order ) ! no neihgborhood
                !     endif
                ! else
                !     ! Multi-state
                !     ! other states
                !     nprojs = self%nprojs
                !     if( present(nnvec) )nprojs = self%nnn
                !     n =  (self%nstates-1)*nprojs
                !     allocate( isrch_order(n) )
                !     irt = ran_tabu( n )
                !     start = 1 - nprojs
                !     end   = 0
                !     do istate=1,self%nstates
                !         if( istate==self%prev_state )cycle
                !         start  = start+nprojs                 ! cumulative index to isrch_order
                !         end    = start+nprojs-1
                !         if( present(nnvec) )then
                !             isrch_order( start:end ) = nnvec + (istate-1)*self%nprojs
                !         else
                !             istart = (istate-1)*self%nprojs+1       ! absolute index to references
                !             iend   = istart+self%nprojs-1
                !             isrch_order( start:end ) = (/ (i,i=istart,iend) /)
                !         endif
                !     enddo
                !     call irt%shuffle( isrch_order )
                !     self%srch_order(1:n) = isrch_order(:)       ! updates global srch_order for other states
                !     deallocate( isrch_order )
                !     call irt%kill
                !     ! previous state
                !     allocate( isrch_order( nprojs ) )
                !     irt = ran_tabu( nprojs )
                !     if( present(nnvec) )then
                !         isrch_order = nnvec + (self%prev_state-1)*self%nprojs
                !     else
                !         istart = (self%prev_state-1)*self%nprojs + 1
                !         iend   = self%prev_state * self%nprojs
                !         isrch_order = (/ (i,i=istart,iend) /)                        
                !     endif
                !     call irt%shuffle( isrch_order )
                !     self%srch_order(n+1:) = isrch_order         ! updates global srch_order for previous state
                !     deallocate( isrch_order )
                !     call irt%kill
                ! endif
                if( any(self%srch_order == 0))stop 'Invalid index in srch_order; simple_prime3d_srch::prep_ref_oris'
                ! Search space
                call prep_discrete_reforis
            case( 'shift' )
                call prep_discrete_reforis      
            case( 'qcont' )
                ! QUASI-CONTINUOUS
                ! TODO: PREVIOUS STATE LAST
                call self%o_refs%new( self%nrefs )
                call self%o_refs%rnd_proj_space( self%nrefs, eullims=self%eullims ) ! stochastic search space
                if( self%nstates > 1 )call self%o_refs%rnd_states( self%nstates )   ! randomise states
            case( 'qcontneigh' )
                ! QUASI-CONTINUOUS & NEIHGBORHOOD
                ! TODO: PREVIOUS STATE LAST
                if( .not.present(o_prev) )stop 'Previous orientation must be provided; simple_prime3D_srch;prep_ref_oris'
                call self%o_refs%new( self%nnnrefs )
                call self%o_refs%rnd_proj_space( self%nnnrefs, o_prev, self%athres, self%eullims ) ! neighborhood stochastic search space
                if( self%nstates > 1 )call self%o_refs%rnd_states( self%nstates )                  ! randomise states
            case DEFAULT
                stop 'Unknown refinement method; prime3D_srch; prep_reforis'                
        end select

        contains

            !>  Discrete search space: multi-state spiral similar to the builder's
            subroutine prep_discrete_reforis
                type(ori)  :: o
                integer    :: cnt, istate, iproj
                call self%o_refs%new( self%nrefs )              ! init references object
                cnt = 0
                do istate=1,self%nstates
                    do iproj=1,self%nprojs
                        cnt = cnt + 1
                        o = self%pe%get_ori( iproj )
                        call o%set( 'state', real(istate) )     ! Updates state
                        call o%set( 'class', real(iproj) )      ! Updates class (proj dir)
                        call self%o_refs%set_ori( cnt,o )
                    enddo
                enddo
            end subroutine prep_discrete_reforis

    end subroutine prep_reforis

    !>  \brief  prepares correlation target (previous best) for the search
    subroutine prep_corr4srch( self, pftcc, iptcl, lp, o_prev, corr_t )
        class(prime3D_srch),        intent(inout) :: self
        class(polarft_corrcalc),    intent(inout) :: pftcc
        integer,                    intent(in)    :: iptcl
        real,                       intent(in)    :: lp
        class(ori),       optional, intent(inout) :: o_prev
        real,             optional, intent(in)    :: corr_t
        real      :: cc_t, cc_t_min_1
        logical   :: calc_corr, did_set_prev_corr
        !
        integer :: nnn, inpl_ind
        real    :: wprev_corr, sumw, angthresh, e3
        real    :: corrs(self%nrots), ws(self%nrots), inpl_dist
        !
        calc_corr = .true.
        if( present(corr_t) )calc_corr = .false.
        if( (self%refine.eq.'qcont'.or.self%refine.eq.'qcontneigh' ).and. &
        & .not.present(corr_t) )stop 'corr to be provided with qcont-/neigh; prep_corr4srch; simple_prime3D_srch'
        ! DEFAULT VALUE
        self%prev_corr = 1.
        if( present(o_prev) )then
            ! CALCULATE PREVIOUS BEST CORRELATION (treshold for better)
            did_set_prev_corr = .false.
            if( calc_corr )then
                cc_t = max( 0., pftcc%corr(self%prev_ref, iptcl, self%prev_roind) ) ! correlation with prev best ori
            else
                cc_t = corr_t
            endif
            if( cc_t > 1. .or. cc_t < -1. .or. isnan(cc_t) )then
                stop 'Invalid correlation value in simple_prime3d_srch::prep_corr4srch'
            endif
            self%prev_corr = cc_t                                                        ! default
            if( (self%refine.eq.'no') .and. self%nstates==1 .and. self%ctf.eq.'no' )then ! moving average for single state only
                cc_t_min_1 = -1.
                if( o_prev%isthere('corr') ) cc_t_min_1 = o_prev%get('corr')
                if( o_prev%isthere('lp') )then
                    if( abs(o_prev%get('lp')-lp) < 0.5 )then                             ! previous and present correlations comparable
                        if( cc_t_min_1 > 0. )then
                            self%prev_corr    = 0.5*cc_t+0.5*cc_t_min_1                  ! diversifying limit
                            did_set_prev_corr = .true.
                        endif
                    endif
                endif
                if( .not. did_set_prev_corr ) self%prev_corr = cc_t                      ! greedy limit
            endif
        endif
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::PREPARED CORRELATION FOR SIMPLE_PRIME3D_SRCH'
    end subroutine prep_corr4srch

    !>  \brief  initialises pftcc in-plane search
    subroutine prep_inpl_srch( self, pftcc)
        class(prime3D_srch),        intent(inout) :: self
        class(polarft_corrcalc),    intent(inout) :: pftcc
        ! Init shift search
        if( self%doshift )then
            ! we are using simplex with barrier constraint
            if( self%greedy_inpl )then
                call pftcc_inplsrch_init( pftcc, self%lims, self%shbarr )
            else
                call pftcc_shsrch_init(   pftcc, self%lims, self%shbarr )
            endif
        endif
    end subroutine prep_inpl_srch

    !>  \brief  fetches ctf params for online application (quasi-continuous only)
    subroutine prep_ctfparms( self, o )
        class(prime3D_srch), intent(inout) :: self
        type(ori), optional, intent(inout) :: o
        self%dfx    = 0.
        self%dfy    = 0.
        self%angast = 0.
        if( present(o) )then
            if( self%ctf.ne.'no' )then
                self%dfx = o%get('dfx')
                if( self%pp%tfplan%mode.eq.'astig' )then ! astigmatic CTF
                    self%dfy    = o%get('dfy')
                    self%angast = o%get('angast')
                else if( self%pp%tfplan%mode.eq.'noastig' )then
                    self%dfy    = self%dfx
                    self%angast = 0.
                else
                    stop 'Unsupported ctf mode; prep_ctfparms; simple_prime3D_srch'
                endif
            endif
        endif
    end subroutine prep_ctfparms

    !>  \brief retrieves and preps npeaks orientations for reconstruction
    subroutine prep_npeaks_oris( self )
        class(prime3D_srch), intent(inout) :: self
        type(ori) :: o, o_new
        real      :: euls(3), shvec(2)
        real      :: corr, ow, frac
        integer   :: ipeak, cnt, ref, state, class
        ! Init
        call self%o_npeaks%new( self%npeaks )
        do ipeak=1,self%npeaks
            ! get ipeak-th ori
            cnt = self%nrefs-self%npeaks+ipeak
            ref = self%proj_space_inds( cnt )
            if( ref < 1 .or. ref > self%nrefs ) stop 'ref index out of bound; simple_prime3D_srch::prep_npeaks_oris'
            o   = self%o_refs%get_ori( ref )
            ! grab info
            state = nint( o%get('state') )
            class = nint( o%get('class') )
            corr  = o%get('corr')
            if( isnan(corr) ) stop 'correlation is NaN in simple_prime3D_srch::prep_npeaks_oris'
            ow    = o%get('ow')
            euls  = o%get_euler()
            ! Add shift
            shvec = self%prev_shvec
            if( self%doshift )shvec = shvec + o%get_shift()
            if( abs(shvec(1)) < 1e-6 ) shvec(1) = 0.
            if( abs(shvec(2)) < 1e-6 ) shvec(2) = 0.
            ! copy info to new ori
            call o_new%new
            call o_new%set_euler( euls )  
            call o_new%set_shift( shvec )
            call o_new%set('state', real(state))
            call o_new%set('class', real(class))
            call o_new%set('corr',  corr)
            call o_new%set('ow',    ow)
            ! stashes in self
            call self%o_npeaks%set_ori( ipeak, o_new )
        enddo
        ! other variables
        if( str_has_substr(self%refine, 'neigh') )then
            frac = 100.*real(self%nrefs_eval) / real(self%nnnrefs)
        else
            frac = 100.*real(self%nrefs_eval) / real(self%nrefs)
        endif            
        call self%o_npeaks%set_all2single('frac',   frac)
        call self%o_npeaks%set_all2single('mi_hard',0.)
        call self%o_npeaks%set_all2single('dist',   0.)
        call self%o_npeaks%set_all2single('sdev', self%ang_sdev() )
        ! ctf parameters
        if( self%ctf.ne.'no' )then
            call self%o_npeaks%set_all2single('dfx',   self%dfx)
            call self%o_npeaks%set_all2single('dfy',   self%dfy)
            call self%o_npeaks%set_all2single('angast',self%angast)
        endif
        ! debug
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED PREP_NPEAKS_ORIS'
    end subroutine prep_npeaks_oris

    !>  \brief  prepares the matrices for PRIME3D search
    subroutine calc_corrs( self, pftcc, mode )
        class(prime3D_srch),        intent(inout) :: self
        class(polarft_corrcalc),    intent(inout) :: pftcc
        character(len=*), optional, intent(in)    :: mode
        if( self%refine .eq. 'shc' )then
            call self%srch_common%calc_corrs(pftcc, self%refine, mode, self%prev_corrs)
        else
           call self%srch_common%calc_corrs(pftcc, self%refine, mode)
        endif        
    end subroutine calc_corrs

    ! SEARCH ROUTINES
    
    !>  \brief a master prime search routine 4 CPU
    subroutine exec_prime3D_srch( self, pftcc, iptcl, lp, o, nnmat, cnt_glob )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        real,                    intent(in)    :: lp
        class(ori), optional,    intent(inout) :: o
        integer,    optional,    intent(in)    :: nnmat(:,:)
        integer,    optional,    intent(in)    :: cnt_glob
        real :: wcorr
        call self%prep4srch( o, nnmat )
        call self%prep_corr4srch( pftcc, iptcl, lp, o )
        call self%prep_inpl_srch( pftcc )
        call self%stochastic_srch( pftcc, iptcl, cnt_glob )
        call self%prep_npeaks_oris
        call self%stochastic_weights( wcorr, self%o_npeaks )
        if( present(o)   ) call o%set('corr', wcorr)
        if( self%doshift ) call self%sort_shifted_npeaks( self%o_npeaks )
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SRCH'
    end subroutine exec_prime3D_srch
    
    !>  \brief a master prime search routine 4 CPU
    subroutine exec_prime3D_qcont_srch( self, refvols, proj, pftcc, iptcl, lp, o, athres )
        use simple_image,     only: image
        use simple_projector, only: projector
        use simple_ctf,       only: ctf
        class(prime3D_srch),     intent(inout) :: self
        class(image),            intent(inout) :: refvols(self%nstates)
        class(projector),        intent(inout) :: proj
        class(ori),              intent(inout) :: o
        class(polarft_corrcalc), intent(inout) :: pftcc
        real,                    intent(in)    :: lp
        integer,                 intent(in)    :: iptcl
        real, optional,          intent(in)    :: athres
        type(ctf) :: tfun
        real      :: cc, wcorr
        if( present(athres) )then
            if( athres <= 0.)then
                write(*,*)'Invalid athres value: ',athres,' in simple_prime3D_srch; exec_prime3D_qcont_srch'
                stop
            endif
            self%athres = athres
        endif
        ! we here need to re-create the CTF object as kV/cs/fraca are now per-particle params
        ! that these parameters are part of the doc is checked in the params class
        tfun = ctf(self%pp%smpd, o%get('kv'), o%get('cs'), o%get('fraca'))
        call self%prep4srch( o )
        call self%prep_ctfparms( o )
        cc = self%calc_qcont_corr( iptcl, o, pftcc, proj, tfun, refvols)
        call self%prep_corr4srch( pftcc, iptcl, lp, o, cc )
        call self%prep_inpl_srch( pftcc )
        call self%stochastic_srch_qcont(refvols, proj, tfun, pftcc, iptcl)
        call self%prep_npeaks_oris
        call self%stochastic_weights(wcorr, self%o_npeaks)
        call o%set('corr', wcorr)
        if( self%doshift )call self%sort_shifted_npeaks( self%o_npeaks )
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::EXECUTED PRIME3D_SRCH'
    end subroutine exec_prime3D_qcont_srch

    !>  \brief a master prime search routine 4 CPU
    subroutine exec_prime3D_shc_srch( self, pftcc, iptcl, lp, o, nnmat, cnt_glob )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        real,                    intent(in)    :: lp
        integer,                 intent(in)    :: iptcl
        class(ori), optional,    intent(inout) :: o
        integer,    optional,    intent(in)    :: nnmat(self%nprojs,self%nnn)
        integer,    optional,    intent(in)    :: cnt_glob
        call self%prep4srch( o, nnmat )
        call self%prep_corr4srch( pftcc, iptcl, lp, o )
        call self%prep_inpl_srch( pftcc )
        if( present(o) )then
            call self%stochastic_srch_shc( pftcc, iptcl, cnt_glob )
        else
            call self%greedy_srch( pftcc, iptcl, cnt_glob )
        endif
        call self%prep_npeaks_oris
        if( present(o) ) call o%set('corr', self%o_npeaks%get( self%npeaks,'corr'))
        if( debug ) write(*,'(A)') '>>> PRIME3D_SHC_SRCH::EXECUTED PRIME3D_SRCH'
    end subroutine exec_prime3D_shc_srch

    !>  \brief a master shift/in-plane prime search routine 4 CPU
    subroutine exec_prime3D_inpl_srch( self, pftcc, iptcl, lp, o, greedy )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        class(ori),              intent(inout) :: o
        real,                    intent(in)    :: lp
        integer,                 intent(in)    :: iptcl
        logical,    optional,    intent(in)    :: greedy
        logical :: doshift, dogreedy, greedy_inpl, found_better
        if( .not.self%doshift )stop 'In-plane search not activated yet; simple_prime3D_srch%exec_prime3D_inpl_srch'
        if( self%npeaks>1 )stop 'In-plane only permitted with npeaks=1; simple_prime3D_srch%exec_prime3D_inpl_srch'
        ! Whether to perform shift or inpl search
        greedy_inpl = self%greedy_inpl
        dogreedy = .false.
        if( present(greedy) )dogreedy = greedy
        self%greedy_inpl = dogreedy
        ! Search
        call self%prep4srch( o )
        call self%prep_corr4srch( pftcc, iptcl, lp, o )
        call self%prep_inpl_srch( pftcc )
        call self%stochastic_srch_inpl( pftcc, iptcl, o )
        call self%prep_npeaks_oris
        call o%set('corr', self%o_npeaks%get( self%npeaks,'corr'))
        self%greedy_inpl = greedy_inpl ! restores setting
        if( debug ) write(*,'(A)') '>>> PRIME3D_SHC_SRCH::EXECUTED PRIME3D_INPL_SRCH'
    end subroutine exec_prime3D_inpl_srch

    !>  \brief  executes the stochastic soft orientation search on CPU
    subroutine stochastic_srch( self, pftcc, iptcl, cnt_glob )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        integer, optional,       intent(in)    :: cnt_glob
        real    :: projspace_corrs( self%nrefs )
        integer :: ref, isample, nrefs
        if( present(cnt_glob) )then
            if( self%use_cpu ) stop 'cnt_glob dummy not required for CPU execution; simple_prime3D :: stochastic_srch'
        endif
        ! initialize
        self%nbetter         = 0
        self%nrefs_eval      = 0
        self%proj_space_inds = 0
        projspace_corrs      = -1.
        ! search
        nrefs = self%nrefs
        if( self%refine.eq.'neigh' )nrefs=self%nnnrefs
        do isample=1,nrefs
            ref = self%srch_order( isample )                             ! set the stochastic reference index
            if( ref==self%prev_ref )cycle                                ! previous best considered last
            if( ref < 1 .or. ref > self%nrefs ) stop 'ref index out of bound; simple_prime3D_srch::stochastic_srch'
            call per_ref_srch( ref )                                     ! actual search
            if( self%nbetter >= self%npeaks ) exit                       ! exit condition
        end do
        if( self%nbetter<self%npeaks )call per_ref_srch( self%prev_ref ) ! evaluate previous best ref last
        call hpsort(self%nrefs, projspace_corrs, self%proj_space_inds)   ! sort in correlation projection direction space
        call self%inpl_srch( iptcl )                                     ! search shifts
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH'

        contains

            ! align using best improving
            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                real                :: corrs(self%nrots), inpl_corr
                integer             :: loc(1), inpl_ind
                if( self%use_cpu )then
                    corrs     = pftcc%gencorrs(iref, iptcl)             ! In-plane correlations
                    loc       = maxloc(corrs)                           ! greedy in-plane
                    inpl_ind  = loc(1)                                  ! in-plane angle index
                    inpl_corr = corrs(loc(1))                           ! max in plane correlation
                else
                    inpl_ind  = self%srch_common%inpl(cnt_glob,iref)
                    inpl_corr = self%srch_common%corr(cnt_glob,iref)
                endif
                call self%store_solution( iref, iref, inpl_ind, inpl_corr )
                projspace_corrs( iref ) = inpl_corr                     ! stash in-plane correlation for sorting
                ! update nbetter to keep track of how many improving solutions we have identified
                if( self%npeaks == 1 )then
                    if( inpl_corr > self%prev_corr ) self%nbetter = self%nbetter+1
                else
                    if( inpl_corr >= self%prev_corr ) self%nbetter = self%nbetter+1
                endif
                ! keep track of how many references we are evaluating
                self%nrefs_eval = self%nrefs_eval+1
            end subroutine per_ref_srch

    end subroutine stochastic_srch

    !>  \brief  executes the stochastic quasi-continuous soft orientation search on CPU
    subroutine stochastic_srch_qcont( self, refvols, proj, tfun, pftcc, iptcl )
        use simple_image,     only: image
        use simple_projector, only: projector
        use simple_ctf,       only: ctf
        class(prime3D_srch),     intent(inout) :: self
        class(image),            intent(inout) :: refvols(self%nstates)
        class(projector),        intent(inout) :: proj
        class(ctf),              intent(inout) :: tfun
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        integer :: iref, nrefs
        real    :: projspace_corrs(self%nrefs)
        ! INIT
        self%nbetter         = 0
        self%nrefs_eval      = 0
        self%proj_space_inds = 0
        projspace_corrs = -1.
        ! ANGULAR SEARCH
        nrefs = self%nrefs
        if( self%refine.eq.'qcontneigh' )nrefs=self%nnnrefs
        do iref=1,nrefs
            if( iref==self%prev_ref )cycle                                ! skips previous reference
            call per_ref_srch( iref )                                     ! do thy thing
            if( self%nbetter>=self%npeaks ) exit   ! exit condition
        end do
        if( self%nbetter<self%npeaks )call per_ref_srch( self%prev_ref )  ! evaluate previous reference last
        call hpsort(self%nrefs, projspace_corrs, self%proj_space_inds)    ! sort in projection direction space
        ! SHIFT SEARCH
        call self%inpl_srch_qcont(refvols, proj, tfun, pftcc, iptcl)
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED QUASI-CONTINUOUS STOCHASTIC SEARCH'

        contains
            
            ! per random reference calculation
            subroutine per_ref_srch( iref )
                integer, intent(in) :: iref
                type(ori)           :: o_ref
                real                :: corrs(self%nrots), inpl_corr
                integer             :: state, loc(1), inpl_ind
                ! Info is stashed on orientation object except for index and correlation
                ! Out-of-plane
                o_ref = self%o_refs%get_ori( iref )                                       ! get reference ori
                state = nint( o_ref%get('state') )
                call proj%fproject_polar( 1, refvols( state ), o_ref, self%pp, pftcc )    ! online generation of polar central section
                if( self%ctf .ne. 'no' )then
                    call pftcc%apply_ctf(self%pp%smpd, tfun, self%dfx, self%dfy, self%angast)
                endif
                ! In-plane
                corrs     = pftcc%gencorrs(1, iptcl)                                      ! in-plane correlations
                loc       = maxloc(corrs)                                                 ! greedy in-plane
                inpl_ind  = loc(1)                                                        ! in-plane index
                inpl_corr = corrs( inpl_ind )
                call self%store_solution( iref, iref, inpl_ind, inpl_corr )
                projspace_corrs( iref ) = inpl_corr                                       ! stash correlation for sorting
                ! keeps track of how many improving solutions identified
                if( self%npeaks == 1 )then
                    if( inpl_corr > self%prev_corr )self%nbetter = self%nbetter+1
                else
                    if( inpl_corr >= self%prev_corr )self%nbetter = self%nbetter+1
                endif
                ! keeps track of how many references have been evaluated
                self%nrefs_eval = self%nrefs_eval+1
            end subroutine per_ref_srch

    end subroutine stochastic_srch_qcont
    
    !>  \brief  executes the stochastic hard orientation search using pure stochastic hill climbing
    !!          (no probabilistic weighting + stochastic search of in-plane angles) on CPU
    subroutine stochastic_srch_shc( self, pftcc, iptcl, cnt_glob )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        integer, optional,       intent(in)    :: cnt_glob
        real      :: inpl_corr
        real      :: state_inpl_corrs(self%nstates), state_probs(self%nstates), corrs(self%nrots)
        integer   :: state_counts(self%nstates),  state_refs(self%nstates), state_inpl_inds(self%nstates)
        integer   :: isample, ref, inpl_ind, nrefs, s, loc(1)
        logical   :: found_better
        if( self%refine.eq.'shcneigh' .or. self%nstates>1 )then
            if( self%bench_gpu .or. self%use_gpu )&
            & stop 'Neigh and multistates modes not implemented for GPU; simple_prime3D :: stochastic_srch'
        endif
        if( present(cnt_glob) )then
            if( self%use_cpu ) stop 'cnt_glob dummy not required for CPU execution; simple_prime3D :: stochastic_srch_shc'
        endif
        ! initialize
        self%proj_space_inds = 0
        self%nrefs_eval = 0
        nrefs = self%nrefs
        if( self%refine.eq.'shcneigh' )nrefs=self%nnnrefs
        ! search
        found_better = .false.                                  ! whether an improving reference has been found
        do isample=1,nrefs
            ref = self%srch_order( isample )                    ! stochastic reference index
            if( ref==self%prev_ref )cycle                       ! previous best considered last
            if( ref < 1 .or. ref > self%nrefs )stop 'ref index out of bound; simple_prime3D_srch::stochastic_srch_shc'
            call per_ref_srch( ref )                            ! actual search
            if( inpl_ind>0 )then
                ! correlation-improving reference found
                call self%store_solution( self%nrefs, ref, inpl_ind, inpl_corr) ! store solution
                found_better = .true.                           ! flag for found solution
                exit                                            ! exit criterion
            endif
        enddo
        if( .not.found_better )then
            ref = self%prev_ref                                 ! if none found defaults to previous best...
            call per_ref_srch( ref )                            ! ...and do calculations
            if( inpl_ind == 0 )then
                ! be greedy
                loc       = maxloc(corrs)
                inpl_ind  = loc(1)
                inpl_corr = corrs(inpl_ind)
            endif
            call self%store_solution(self%nrefs, ref, inpl_ind, inpl_corr) ! store solution
        endif
        ! search shifts
        call self%inpl_srch( iptcl )
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED STOCHASTIC SEARCH (REFINE=SHC|SHCNEIGH)'
        
        contains
        
            ! align using SHC
            subroutine per_ref_srch( iref )
                use simple_rnd, only: shcloc
                integer, intent(in) :: iref
                inpl_ind  = 0                                             ! init in-plane index 
                inpl_corr = 0.                                            ! init correlation
                if( self%use_cpu )then
                    corrs     = pftcc%gencorrs( iref, iptcl )             ! in-plane correlations 
                    inpl_ind  = shcloc(self%nrots, corrs, self%prev_corr) ! first improving in-plane index
                    if( inpl_ind > 0 )inpl_corr = corrs( inpl_ind )       ! improving correlation found
                else
                    inpl_ind  = self%srch_common%inpl( cnt_glob, iref )
                    inpl_corr = self%srch_common%corr( cnt_glob, iref )
                    ! Multi-states TODO
                endif
                self%nrefs_eval = self%nrefs_eval+1                       ! updates fractional search space
            end subroutine per_ref_srch

    end subroutine stochastic_srch_shc

    subroutine stochastic_srch_inpl( self, pftcc, iptcl, o )
        use simple_rnd, only: ran3
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        class(ori),              intent(inout) :: o
        integer,                 intent(in)    :: iptcl
        type(ori)      :: o_tmp
        type(ran_tabu) :: rt
        real           :: target_corr, corr, e3, delta, p, rnd
        integer        :: i, istate, ref, nrefs, isrch_order(self%nstates), loc(1)
        logical        :: found_better
        ! init
        nrefs = self%nrefs
        if( str_has_substr(self%refine,'neigh') )nrefs = self%nnn
        found_better         = .false.
        self%nrefs_eval      = 0
        self%proj_space_inds = 0
        !e3 = 360.-self%srch_common%rot( self%prev_roind )               ! psi
        if( self%nstates==1 )then
            ! single state
            call self%store_solution( self%nrefs, self%prev_ref, self%prev_roind, self%prev_corr ) 
            self%nrefs_eval = nrefs                                     ! updates frac
        else
            ! multi-state, first improving
            rt = ran_tabu( self%nstates )
            call rt%ne_ran_iarr( isrch_order )
            call rt%kill
            do i=1,self%nstates
                istate = isrch_order(i)                                 ! randomized state search order
                if( istate==self%prev_state )cycle                      ! previous state last
                self%nrefs_eval = self%nrefs_eval+1                     ! updates frac
                ref  = (istate-1)*self%nprojs + self%prev_proj          ! reference index to state projection direction
                corr = pftcc%corr( ref, iptcl, self%prev_roind )        ! state correlation
                if( corr >= self%prev_corr )then
                    found_better = .true.                               ! first improvement exit condition
                    exit
                endif
            enddo
            if( .not. found_better )then
                ref  = self%prev_ref                                    ! defaults to previous state
                corr = self%prev_corr
                self%nrefs_eval = self%nstates
            endif
            call self%store_solution( self%nrefs, ref, self%prev_roind, corr )       
            self%nrefs_eval = nint( real(nrefs*self%nrefs_eval)/real(self%nstates) ) ! updates frac
        endif
        call self%inpl_srch( iptcl )                                    ! performs in-plane search
    end subroutine stochastic_srch_inpl


    !>  \brief  greedy hill climbing (4 initialisation)
    subroutine greedy_srch( self, pftcc, iptcl, cnt_glob )
        class(prime3D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        integer, optional,       intent(in)    :: cnt_glob
        type(ori) :: o_best
        integer   :: i, ref, inpl_ind, loc(1)
        real      :: corrs(self%nrots), inpl_corr, e3, corr_best
        if( present(cnt_glob) )then
            if( self%use_cpu ) stop 'cnt_glob dummy not required for CPU execution; simple_prime3D :: greedy_srch'
        endif
        ! initialize
        self%nrefs_eval  = 0
        inpl_ind         = 0
        inpl_corr        = 0
        corr_best        = -1.
        ! systematic search
        do i=1,self%nrefs
            ref = i                                               ! set reference index
            ! align
            if( self%use_cpu )then
                corrs     = pftcc%gencorrs(ref, iptcl)            ! in-plane correlations
                loc       = maxloc(corrs)                         ! best improving correlation index
                inpl_ind  = loc(1)
                inpl_corr = corrs(inpl_ind)                       ! best improving correlation
            else
                inpl_ind  = self%srch_common%inpl(cnt_glob,ref)
                inpl_corr = self%srch_common%corr(cnt_glob,ref)
            endif
            self%nrefs_eval = i                                   ! fractional search space
            if( inpl_corr > corr_best )then
                self%proj_space_inds( self%nrefs ) = ref
                corr_best = inpl_corr
                o_best    = self%o_refs%get_ori( ref )            ! get reference
                e3        = 360.-self%srch_common%rot( inpl_ind ) ! psi
                call o_best%e3set( e3 )                           ! stash psi
                call o_best%set( 'corr', corr_best)               ! stash best in-plane correlation
                call o_best%set( 'ow', 1.)                        ! sets recontruction weight
                call self%o_refs%set_ori( ref, o_best )           ! stash new best reference
            endif
        end do
        call self%inpl_srch(iptcl)                               ! search shifts
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED GREEDY SEARCH 4 INITI (REFINE=SHC)'
    end subroutine greedy_srch
    
    !>  \brief  executes the in-plane search for discrete mode
    subroutine inpl_srch( self, iptcl )
        class(prime3D_srch), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        type(ori) :: o
        real      :: cc, cxy(3), crxy(4)
        integer   :: i, ref, inpl_ind
        if( self%doshift )then
            do i=self%nrefs,self%nrefs-self%npeaks+1,-1
                ref      = self%proj_space_inds( i )
                o        = self%o_refs%get_ori( ref )
                cc       = o%get('corr')
                inpl_ind = self%srch_common%roind( 360.-o%e3get() )
                if( self%greedy_inpl )then
                    call pftcc_inplsrch_set_indices(ref, iptcl)
                    crxy = pftcc_inplsrch_minimize(inpl_ind)
                    if( crxy(1) >= cc )then
                        call o%set( 'corr', crxy(1) )
                        call o%e3set( 360.-crxy(2) )
                        call o%set_shift( crxy(3:4) )
                        call self%o_refs%set_ori( ref, o )
                    endif
                else
                    call pftcc_shsrch_set_indices(ref, iptcl, inpl_ind )
                    cxy = pftcc_shsrch_minimize()
                    if( cxy(1) >= cc )then
                        call o%set('corr', cxy(1))
                        call o%set_shift( cxy(2:3) )
                        call self%o_refs%set_ori( ref, o )
                    endif
                endif
            end do
        endif
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH::FINISHED INPL SEARCH'
    end subroutine inpl_srch

    !>  \brief  executes the in-plane search for quasi-continuous mode
    subroutine inpl_srch_qcont( self, refvols, proj, tfun, pftcc, iptcl )
        use simple_image,     only: image
        use simple_projector, only: projector
        use simple_ctf,       only: ctf
        class(prime3D_srch),     intent(inout) :: self
        class(image),            intent(inout) :: refvols(self%nstates)
        class(projector),        intent(inout) :: proj
        class(polarft_corrcalc), intent(inout) :: pftcc
        class(ctf),              intent(inout) :: tfun
        integer,                 intent(in)    :: iptcl
        type(ori) :: o, o_zero_e3
        integer   :: i, ref, state, inpl_ind
        real      :: cxy(3), crxy(4)
        if( self%doshift )then
            ! inpl_ind = self%srch_common%roind( 360. ) ! the in-plane angle is taken care of at the projection level
            do i=self%nrefs,self%nrefs-self%npeaks+1,-1
                ref      = self%proj_space_inds(i)
                o        = self%o_refs%get_ori(ref)
                inpl_ind = self%srch_common%roind( o%e3get() )
                state    = nint( o%get('state') )
                ! online generation of polar central section
                o_zero_e3 = o
                call o_zero_e3%e3set(0.)
                call proj%fproject_polar(1, refvols(state), o_zero_e3, self%pp, pftcc )
                if( self%ctf .ne. 'no' )then
                    call pftcc%apply_ctf(self%pp%smpd, tfun, self%dfx, self%dfy, self%angast)
                endif
                ! in-plane search
                if( self%greedy_inpl )then
                    call pftcc_inplsrch_set_indices(1, iptcl)
                    crxy = pftcc_inplsrch_minimize(inpl_ind)
                    call o%set( 'corr', crxy(1) )
                    call o%e3set( 360.-crxy(2) )
                    call o%set_shift( crxy(3:4) )
                    call self%o_refs%set_ori( ref, o )
                else
                    call pftcc_shsrch_set_indices(1, iptcl, inpl_ind )
                    cxy = pftcc_shsrch_minimize()
                    call o%set( 'corr', cxy(1) )
                    call o%set_shift( cxy(2:3) )
                    call self%o_refs%set_ori( ref, o )
                endif
            end do
        endif
    end subroutine inpl_srch_qcont

    subroutine stochastic_weights( self, wcorr, os )
        class(prime3D_srch), intent(inout) :: self
        type(oris),          intent(inout) :: os
        real,                intent(out)   :: wcorr
        real             :: ws(self%npeaks)
        real,allocatable :: corrs(:)
        integer          :: i
        if( os%get_noris() /= self%npeaks )stop 'invalid number of references in simple_prime3D_srch::stochastic_weights'
        ws    = 0.
        wcorr = 0.
        ! get unnormalised correlations
        corrs = os%get_all('corr')
        ! calculate normalised weights and weighted corr
        where( corrs > TINY ) ws = exp(corrs) ! ignore invalid corrs
        ws    = ws/sum(ws)
        wcorr = sum(ws*corrs) 
        ! update npeaks individual weights
        do i=1,self%npeaks
            call os%set(i,'ow',ws(i))
        enddo
        deallocate(corrs)
        if( debug ) write(*,'(A)') '>>> PRIME3D_SRCH :: CALCULATED STOCHASTIC WEIGHTS PRIME3D_SRCH'        
    end subroutine stochastic_weights

    ! stores a solution
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

    !>  \brief  is for getting o_npeaks
    function get_o_npeaks( self )result( out_os )
        class(prime3D_srch), intent(inout) :: self
        type(oris) :: out_os
        out_os = self%o_npeaks
    end function get_o_npeaks

    !>  \brief  is for getting the n-last reference orientations
    function get_o_refs( self, n )result( out_os )
        class(prime3D_srch), intent(inout) :: self
        integer, optional,   intent(in)    :: n
        type(oris) :: out_os
        integer    :: i, cnt, in
        if( present(n) )then
            if( n < 1 )then
                stop 'invalid index in prime3D_srch::get_o_refs, 1'
            elseif( self%refine.eq.'qcontneigh' .and. n>self%nnnrefs )then
                stop 'invalid index in prime3D_srch::get_o_refs, 2'
            elseif( n>self%nrefs)then
                stop 'invalid index in prime3D_srch::get_o_refs, 3'
            endif
            in = n
        else
            if( self%refine.eq.'qcontneigh' )then
                in = self%nnnrefs
            else
                in = self%nrefs
            endif
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
        if( self%exists )then
            self%pe => null()
            self%pp => null()
            if( allocated( self%proj_space_inds) ) deallocate(self%proj_space_inds)
            if( allocated( self%srch_order )     ) deallocate(self%srch_order)
            if( allocated( self%prev_corrs )     ) deallocate(self%prev_corrs)
            call self%srch_common%kill
            call self%rt%kill
            call self%o_refs%kill
            call self%o_npeaks%kill
            self%exists = .false.
        endif        
    end subroutine kill

end module simple_prime3D_srch
