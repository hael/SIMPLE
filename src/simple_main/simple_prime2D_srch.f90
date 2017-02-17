module simple_prime2D_srch
use simple_defs
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_prime_srch,       only: prime_srch
use simple_ran_tabu,         only: ran_tabu
use simple_strings,          only: str_has_substr
use simple_math              ! use all in there
use simple_pftcc_shsrch      ! use all in there
use simple_pftcc_inplsrch    ! use all in there
implicit none

public :: prime2D_srch
private

logical, parameter :: debug = .false.

type reference
    integer :: class = 1     !< class index
    real    :: inpl_corr     !< best in-plane correlation
    integer :: inpl_ind      !< best in-plane index
    real    :: shvec(2) = 0. !< shift vector
end type reference

type prime2D_srch
    private
    type(prime_srch)      :: srch_common            !< functionalities common to primesrch2D/3D
    integer               :: nrefs      = 0         !< number of references
    integer               :: nrots      = 0         !< number of in-plane rotations in polar representation
    integer               :: nnn        = 0         !< number of nearest neighbours
    integer               :: nrefs_eval = 0         !< nr of references evaluated
    integer               :: class      = 0         !< previous class index
    integer               :: ref        = 0         !< previous reference index (same as class)
    integer               :: rot        = 0         !< previous in-plane rotation index
    real                  :: trs        = 0.        !< shift range parameter [-trs,trs]
    real                  :: shvec(2)   = 0.        !< previous origin shift vector
    real                  :: corr_prev  = 1.        !< previous best correlation
    type(reference)       :: best_ref               !< best reference found
    integer, allocatable  :: srch_order(:)          !< stochastic search order
    real, allocatable     :: prev_corrs(:)          !< array of previous correlations required for GPU exec
    type(ran_tabu)        :: rt                     !< object for random number generation
    character(len=STDLEN) :: refine                 !< refinement flag
    logical               :: use_cpu     = .true.   !< indicates if CPU logics is being used
    logical               :: bench_gpu   = .false.  !< indicates if GPU logics should be tested
    logical               :: use_gpu     = .false.  !< indicates if GPU logics should be used
    logical               :: doinpl      = .true.   !< inplane params search indicator
    logical               :: doshift     = .true.   !< origin shift search indicator
    logical               :: exists      = .false.  !< 2 indicate existence
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! GETTERS
    procedure :: get_nrots
    procedure :: get_corr
    procedure :: get_inpl
    procedure :: get_cls
    procedure :: get_prev_corrs
    procedure :: get_roind
    ! SETTERS
    procedure :: set_use_cpu
    ! PREPARATION ROUTINES
    procedure :: prep4srch
    procedure :: prepcorrs4gpusrch
    procedure :: calc_corrs
    ! SEARCH ROUTINES
    procedure :: exec_prime2D_srch
    procedure :: stochastic_srch_shc
    procedure :: greedy_srch
    procedure :: shift_srch
    ! DESTRUCTOR
    procedure :: kill
end type prime2D_srch

contains

    ! CONSTRUCTOR
    
    !>  \brief  is a constructor
    subroutine new( self, p, ncls )
        use simple_params, only: params
        class(prime2D_srch), intent(inout) :: self !< instance
        class(params),       intent(in)    :: p    !< parameters
        integer, optional,   intent(in)    :: ncls !< overriding dummy (ncls)
        integer :: alloc_stat, i
        real    :: dang
        ! destroy possibly pre-existing instance
        call self%kill
        ! set constants
        if( present(ncls) )then
            self%nrefs = ncls
        else
            self%nrefs = p%ncls
        endif
        self%nrots      = round2even(twopi*real(p%ring2))
        self%refine     = p%refine
        self%nnn        = p%nnn
        self%nrefs_eval = 0
        self%trs        = p%trs
        self%doshift    = p%doshift
        self%doinpl     = .true.
        if( p%srch_inpl .eq. 'no') self%doinpl = .false.
        if( p%bench_gpu .eq. 'yes' .or. p%use_gpu .eq. 'yes' )then
            if( p%bench_gpu .eq. 'yes' ) self%bench_gpu = .true.
            if( p%use_gpu   .eq. 'yes' ) self%use_gpu = .true.
            if( p%top-p%fromp+1 /= self%nrefs )&
            &stop 'the particle chunk is not correctly balanced for GPU execution!'
            self%use_cpu = .false.
        endif
        ! construct composites
        self%srch_common = prime_srch(p, self%nrefs, self%nrots)
        if( str_has_substr(self%refine,'neigh') )then
            if( self%bench_gpu .or. self%use_gpu )&
            stop 'Neigh modes not implemented for GPU; simple_prime2D_srch :: new'
            allocate(self%srch_order(self%nnn), stat=alloc_stat)
            call alloc_err('In: new; simple_prime2d_srch, 3', alloc_stat)
            self%srch_order = 0
            self%rt = ran_tabu(self%nnn)
        else
            allocate(self%srch_order(self%nrefs), stat=alloc_stat)
            call alloc_err('In: new; simple_prime2d_srch, 4', alloc_stat)
            self%srch_order = 0
            self%rt = ran_tabu(self%nrefs)
        endif
        ! initialize best
        self%best_ref%class     = 0
        self%best_ref%inpl_corr = -1.
        self%best_ref%inpl_ind  = 0
        ! the instance now exists
        self%exists = .true.
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::CONSTRUCTED NEW SIMPLE_prime2D_srch OBJECT'
    end subroutine new

    ! GETTERS

    !>  \brief nrots getter
    pure integer function get_nrots( self )
        class(prime2D_srch), intent(in) :: self
        get_nrots = self%nrots
    end function get_nrots

    !>  \brief correlation getter
    real function get_corr( self, cnt_glob, iref )
        class(prime2D_srch), intent(in) :: self
        integer,             intent(in) :: cnt_glob, iref
        get_corr = self%srch_common%corr(cnt_glob, iref)
    end function get_corr
    
    !>  \brief in-plane index getter
    integer function get_inpl( self, cnt_glob, iref )
        class(prime2D_srch), intent(in) :: self
        integer,             intent(in) :: cnt_glob, iref
        get_inpl = self%srch_common%inpl(cnt_glob, iref)
    end function get_inpl
    
    !>  \brief  to get the class
    subroutine get_cls( self, o )
        use simple_math, only: myacos, rotmat2d, rad2deg
        use simple_ori,  only: ori
        class(prime2D_srch), intent(in)    :: self
        class(ori),          intent(inout) :: o
        real    :: euls(3), mi_class, mi_inpl, mi_joint
        integer :: class, rot
        real    :: x, y, mat(2,2), u(2), x1(2), x2(2)
        ! make unit vector
        u(1)     = 0.
        u(2)     = 1.
        ! calculate previous vec
        mat      = rotmat2d(o%e3get())
        x1       = matmul(u,mat)
        ! get new indices
        class    = self%best_ref%class
        rot      = self%best_ref%inpl_ind
        ! get in-plane angle
        euls     = 0.
        euls(3)  = 360.-self%srch_common%rot(rot) ! change sgn to fit convention
        if( euls(3) == 360. ) euls(3) = 0.
        call o%set_euler(euls)
        ! calculate new vec & distance (in degrees)
        mat      = rotmat2d(o%e3get())
        x2       = matmul(u,mat)
        ! calculate overlap between distributions
        mi_class = 0.
        mi_inpl  = 0.
        mi_joint = 0.
        if( self%class == class )then
            mi_class = mi_class + 1.
            mi_joint = mi_joint + 1.
        endif
        if( self%rot == rot )then
            mi_inpl  = mi_inpl  + 1.
            mi_joint = mi_joint + 1.
        endif 
        mi_joint = mi_joint/2.
        ! set parameters
        x = self%shvec(1)
        y = self%shvec(2)
        if( self%doshift )then
            ! shifts must be obtained by vector addition
            x = x+self%best_ref%shvec(1)
            y = y+self%best_ref%shvec(2)
        endif
        call o%set('x',         x)
        call o%set('y',         y)
        call o%set('class',     real(class))
        call o%set('corr',      self%best_ref%inpl_corr)
        call o%set('dist_inpl', rad2deg(myacos(dot_product(x1,x2))))
        call o%set('mi_class',  mi_class)
        call o%set('mi_inpl',   mi_inpl)
        call o%set('mi_joint',  mi_joint)
        if( str_has_substr(self%refine,'neigh') )then
            call o%set('frac', 100.*(real(self%nrefs_eval)/real(self%nnn)))
        else
            call o%set('frac', 100.*(real(self%nrefs_eval)/real(self%nrefs)))
        endif        
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::GOT BEST ORI'
    end subroutine get_cls

    !>  \brief  is for getting the previous correlations
    function get_prev_corrs( self ) result( prev_corrs )
        class(prime2D_srch), intent(inout) :: self
        real, allocatable :: prev_corrs(:)
        if( allocated(self%prev_corrs) )then
            allocate(prev_corrs(size(self%prev_corrs)), source=self%prev_corrs)
        else
            stop 'self%prev_corrs not allocated; simple_prime2D_srch :: get_prev_corrs'
        endif
    end function get_prev_corrs

    !>  \brief returns the in-plane rotational index for the rot in-plane angle
    integer function get_roind( self, rot )
        class(prime2D_srch), intent(in) :: self
        real,                intent(in) :: rot
        get_roind = self%srch_common%roind(rot)
    end function get_roind

    !>  \brief setter for the use_cpu flag (used by the tester routine)
    subroutine set_use_cpu( self, flag )
        class(prime2D_srch), intent(inout) :: self
        logical,             intent(in)    :: flag
        self%use_cpu = flag
    end subroutine set_use_cpu

    ! PREPARATION ROUTINES

    !>  \brief  prepares for the search
    subroutine prep4srch( self, pftcc, iptcl, lp, o_prev, nnmat )
        use simple_ori, only: ori
        class(prime2D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        real,                    intent(in)    :: lp
        class(ori), optional,    intent(inout) :: o_prev
        integer,    optional,    intent(in)    :: nnmat(self%nrefs,self%nnn)
        real :: lims(2,2)
        if( str_has_substr(self%refine,'neigh') .and. .not.present(nnmat) )&
        & stop 'nnmat must be provided with refine=neigh modes'
        if( self%refine .eq. 'no' )then
            ! make random reference direction order here
            call self%rt%reset
            call self%rt%ne_ran_iarr(self%srch_order)
        endif
        ! initialize in-plane search classes
        lims(1,1) = -self%trs
        lims(1,2) =  self%trs
        lims(2,1) = -self%trs
        lims(2,2) =  self%trs
        ! ininialise 
        call pftcc_shsrch_init(   pftcc, lims )
        call pftcc_inplsrch_init( pftcc, lims )
        if( present(o_prev) )then
            ! find previous discrete alignment parameters
            self%class = nint(o_prev%get('class'))                    ! class index
            self%rot   = self%srch_common%roind(360.-o_prev%e3get())  ! in-plane angle index
            self%ref   = self%class                                   ! reference index
            self%shvec = [o_prev%get('x'),o_prev%get('y')]            ! shift vector
            ! set best to previous best by default
            self%best_ref%class    = self%class         
            self%best_ref%inpl_ind = self%rot 
            ! calculate previous best corr (treshold for better)
            self%corr_prev = pftcc%corr(self%ref, iptcl, self%rot)
            self%best_ref%inpl_corr = self%corr_prev
            if( str_has_substr(self%refine,'neigh') )then
                ! the srch_order depends on the previous class index
                ! so the order needs to be established here
                self%srch_order = nnmat(self%class,:)
                if( any(self%srch_order == 0) ) stop 'Invalid index in srch_order; simple_prime2D_srch :: prep4srch'
                call self%rt%shuffle( self%srch_order )               ! wizz it up
            endif
            ! put prev_best last to avoid cycling
            call put_last(self%class, self%srch_order)
        else
            self%class     = 1
            self%rot       = 1
            self%ref       = 1
            self%shvec     = 0.
            self%corr_prev = 1.
        endif
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::PREPARED FOR SIMPLE_prime2D_srch'
    end subroutine prep4srch

    !>  \brief  in preparation for the GPU search we need to calculate all previous corrs
    subroutine prepcorrs4gpusrch( self, pftcc, a, pfromto )
        use simple_oris, only: oris
        class(prime2D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        class(oris),             intent(inout) :: a
        integer,                 intent(in)    :: pfromto(2)
        real    :: corrs_local(pfromto(1):pfromto(2))
        integer :: iptcl, ncorrs, alloc_stat, rot, ref
        !$omp parallel do schedule(auto) default(shared) private(iptcl,rot,ref)
        do iptcl=pfromto(1),pfromto(2)
            ! find previous discrete alignment parameters
            rot = self%srch_common%roind(360.-a%e3get(iptcl)) ! in-plane angle index
            ref = nint(a%get(iptcl, 'class'))                 ! reference index
            ! calculate previous best corr (treshold for better)
            corrs_local(iptcl) = pftcc%corr(ref, iptcl, rot)
        end do 
        !$omp end parallel do
        if( allocated(self%prev_corrs) ) deallocate(self%prev_corrs)
        ncorrs = pfromto(2)-pfromto(1)+1
        allocate(self%prev_corrs(ncorrs), source=corrs_local(pfromto(1):pfromto(2)), stat=alloc_stat)
        call alloc_err("In: simple_prime2D_srch :: prep4gpusrch", alloc_stat)
    end subroutine prepcorrs4gpusrch
    
    !>  \brief  prepares the matrices for PRIME2D search
    subroutine calc_corrs( self, pftcc, a, pfromto, mode )
        use simple_oris, only: oris
        class(prime2D_srch),        intent(inout) :: self
        class(polarft_corrcalc),    intent(inout) :: pftcc
        class(oris),      optional, intent(inout) :: a
        integer,          optional, intent(in)    :: pfromto(2)
        character(len=*), optional, intent(in)    :: mode
        if( present(a) )then
            if( present(pfromto) )then
                call self%prepcorrs4gpusrch(pftcc, a, pfromto)
                call self%srch_common%calc_corrs(pftcc, 'shc', mode, self%prev_corrs)
            else
                stop 'need optional pfromto (particle index range) present; simple_prime2D_srch :: calc_corrs'
            endif
        else
            call self%srch_common%calc_corrs(pftcc, 'no', mode)
        endif    
    end subroutine calc_corrs

    ! SEARCH ROUTINES
    
    !>  \brief a master prime search routine
    subroutine exec_prime2D_srch( self, pftcc, iptcl, lp, o, nnmat, cnt_glob )
        use simple_ori, only: ori
        class(prime2D_srch),      intent(inout) :: self
         class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                  intent(in)    :: iptcl
        real,                     intent(in)    :: lp
        class(ori), optional,     intent(inout) :: o
        integer,    optional,     intent(in)    :: nnmat(:,:)
        integer,    optional,     intent(in)    :: cnt_glob
        call self%prep4srch(pftcc, iptcl, lp, o, nnmat)
        if( .not.present(o) .or. self%refine.eq.'greedy' )then
            call self%greedy_srch(pftcc, iptcl, cnt_glob)
        else
            call self%stochastic_srch_shc(pftcc, iptcl, cnt_glob)
        endif
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::EXECUTED PRIME2D_SRCH'
    end subroutine exec_prime2D_srch
    
    !>  \brief  executes the greedy rotational search
    subroutine greedy_srch( self, pftcc, iptcl, cnt_glob )
        class(prime2D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        integer, optional,       intent(in)    :: cnt_glob
        integer :: ref, loc(1), inpl_ind, endit, i
        real    :: cxy(3), corrs(self%nrots), inpl_corr
        if( present(cnt_glob) )then
            if( self%use_cpu ) stop 'cnt_glob dummy not required for CPU execution; simple_prime2D :: greedy_srch'
        endif
        self%corr_prev = -1.
        endit = self%nrefs
        if( str_has_substr(self%refine,'neigh') ) endit = self%nnn
        do i=1,endit
            ref = self%srch_order(i)
            if( self%use_cpu )then
                corrs     = pftcc%gencorrs(ref, iptcl)
                loc       = maxloc(corrs)
                inpl_ind  = loc(1)
                inpl_corr = corrs(inpl_ind)
            else
                inpl_ind  = self%srch_common%inpl(cnt_glob,ref)
                inpl_corr = self%srch_common%corr(cnt_glob,ref)
            endif
            if( inpl_corr >= self%corr_prev )then
                ! update the class
                self%best_ref%class     = ref
                ! update the correlations
                self%best_ref%inpl_corr = inpl_corr
                self%corr_prev          = self%best_ref%inpl_corr
                ! update the in-plane angle
                self%best_ref%inpl_ind  = inpl_ind
            endif
        end do
        ! we always evaluate all references using the greedy approach
        self%nrefs_eval = self%nrefs
        if( str_has_substr(self%refine,'neigh') ) self%nrefs_eval = self%nnn
        ! search in-plane
        call self%shift_srch(iptcl)
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::FINISHED GREEDY SEARCH'
    end subroutine greedy_srch
    
    !>  \brief  executes the stochastic rotational search
    subroutine stochastic_srch_shc( self, pftcc, iptcl, cnt_glob )
        use simple_rnd, only: shcloc
        class(prime2D_srch),     intent(inout) :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer,                 intent(in)    :: iptcl
        integer, optional,       intent(in)    :: cnt_glob
        integer :: i, ref, inpl_ind, loc(1), endit
        real    :: corr_new, cxy(3), corrs(self%nrots), inpl_corr
        logical :: found_better
        if( present(cnt_glob) )then
            if( self%use_cpu ) stop 'cnt_glob dummy not required for CPU execution; simple_prime2D :: stochastic_srch_shc'
        endif
        ! initialize
        self%nrefs_eval = 0
        found_better    = .false.
        endit           = self%nrefs
        if( str_has_substr(self%refine,'neigh') ) endit = self%nnn
        ! search
        do i=1,endit
            ref = self%srch_order(i)
            if( self%doinpl )then
                ! keep track of how many references we are evaluating
                self%nrefs_eval = self%nrefs_eval+1
                if( self%use_cpu )then
                    corrs     = pftcc%gencorrs(ref, iptcl)
                    inpl_ind  = shcloc(self%nrots, corrs, self%corr_prev)
                    inpl_corr = 0.
                    if( inpl_ind > 0 ) inpl_corr = corrs(inpl_ind)
                else
                    inpl_ind  = self%srch_common%inpl( cnt_glob, ref )
                    inpl_corr = self%srch_common%corr( cnt_glob, ref )
                endif
                if( inpl_ind > 0 )then
                    ! update the class
                    self%best_ref%class = ref
                    ! update the correlation
                    self%best_ref%inpl_corr = inpl_corr
                    ! update the in-plane angle
                    self%best_ref%inpl_ind  = inpl_ind
                    ! indicate that we found a better solution
                    found_better = .true.
                    exit ! first-improvement heuristic
                endif
            else
                if( .not. self%use_cpu )then
                    stop 'this mode of search (srch_inpl=no) not implemented 4 GPU; simple_prime2D :: stochastic_srch_shc'
                endif
                corr_new = pftcc%corr(ref, iptcl, self%rot)
                ! keep track of how many references we are evaluating
                self%nrefs_eval = self%nrefs_eval+1
                if( corr_new > self%corr_prev )then
                    ! update the class
                    self%best_ref%class     = ref
                    ! update the correlation
                    self%best_ref%inpl_corr = corr_new
                    ! update of the in-plane angle not needed but put in for consistency
                    self%best_ref%inpl_ind  = self%rot
                    ! indicate that we found a better solution
                    found_better = .true.
                    exit ! first-improvement heuristic
                endif
            endif
        end do
        if( found_better )then
            ! best ref has already been updated
        else
            ! keep the old parameters
            self%best_ref%class     = self%class 
            self%best_ref%inpl_corr = self%corr_prev
            self%best_ref%inpl_ind  = self%rot
        endif
        if( self%doinpl )then
            call self%shift_srch(iptcl)
        else
            self%best_ref%shvec = [0.,0.]
        endif
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::FINISHED STOCHASTIC SEARCH'
    end subroutine stochastic_srch_shc

    !>  \brief  executes the in-plane search over one reference
    subroutine shift_srch( self, iptcl )
        class(prime2D_srch), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        real :: cxy(3)
        if( self%doshift )then
            call pftcc_shsrch_set_indices(self%best_ref%class, iptcl, self%best_ref%inpl_ind)
            cxy = pftcc_shsrch_minimize()
            self%best_ref%inpl_corr = cxy(1)
            self%best_ref%shvec     = cxy(2:3)
        else
            self%best_ref%shvec = [0.,0.]
        endif
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::FINISHED SHIFT SEARCH'
    end subroutine shift_srch

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(prime2D_srch), intent(inout) :: self !< instance
        if( self%exists )then
            if( allocated(self%prev_corrs) ) deallocate(self%prev_corrs)
            deallocate(self%srch_order)
            call self%srch_common%kill
            call self%rt%kill
            self%exists = .false.
        endif
    end subroutine kill

end module simple_prime2D_srch
