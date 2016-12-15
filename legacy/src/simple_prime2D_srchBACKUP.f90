module simple_prime2D_srch
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_oris,             only: oris
use simple_ran_tabu,         only: ran_tabu
use simple_jiffys,           only: alloc_err
use simple_math,             only: hpsort, rad2Deg, round2even, put_last
use simple_pftcc_shsrch      ! singleton
use simple_defs              ! singleton
implicit none

public :: prime2D_srch
private

logical, parameter :: debug=.false.

type reference
    integer :: class = 1     !< class index
    real    :: inpl_corr     !< best in-plane correlation
    integer :: inpl_ind      !< best in-plane index
    real    :: shvec(2) = 0. !< shift vector
end type

type prime2D_srch
    private
    integer               :: nrefs      = 0     !< number of references
    integer               :: nrots      = 0     !< number of in-plane rotations in polar representation
    integer               :: nrefs_eval = 0     !< nr of references evaluated
    integer               :: class      = 0     !< previous class index
    integer               :: ref        = 0     !< previous reference index (same as class)
    integer               :: rot        = 0     !< previous in-plane rotation index
    real                  :: trs        = 0.    !< shift range parameter [-trs,trs]
    real                  :: shvec(2)   = 0.    !< previous origin shift vector
    real                  :: corr_prev  = 1.    !< previous best correlation
    type(reference)       :: best_ref           !< best reference found
    integer, allocatable  :: srch_order(:)      !< stochastoc search order
    real,    allocatable  :: angtab(:)          !< list of in-plane angles
    type(ran_tabu)        :: rt                 !< object for random number generation
    character(len=STDLEN) :: refine             !< refinement flag
    logical               :: doinpl  = .true.   !< inplane params search indicator
    logical               :: doshift = .true.   !< origin shift search indicator
    logical               :: exists  = .false.  !< 2 indicate existence
  contains
    ! CONSTRUCTOR & INITIALIZER
    procedure :: new
    ! CALCULATORS
    procedure, private :: calc_roind
    procedure, private :: calc_rot
    ! GETTERS & SETTERS
    procedure :: get_cls
    ! SEARCH ROUTINES & WEIGHT CALCULATION
    procedure :: exec_prime2D_srch
    procedure, private :: prep4srch
    procedure, private :: stochastic_srch
    procedure, private :: greedy_srch
    ! DESTRUCTOR
    procedure :: kill
end type

contains
    
    !>  \brief  is a constructor
    subroutine new( self, p )
        use simple_params, only: params
        class(prime2D_srch), intent(inout) :: self !< instance
        class(params), intent(in)         :: p    !< parameters
        integer :: alloc_stat, cnt, i
        real    :: dang
        ! destroy possibly pre-existing instance
        call self%kill
        ! set constants
        self%nrefs      = p%ncls
        self%nrots      = round2even(twopi*real(p%ring2))
        self%refine     = p%refine
        self%nrefs_eval = 0
        self%trs        = p%trs
        self%doshift    = p%doshift
        self%doinpl     = .true.
        if( p%srch_inpl .eq. 'no') self%doinpl = .false.        
        ! construct composites
        allocate( self%angtab(self%nrots), stat=alloc_stat )
        call alloc_err('In: new; simple_prime2D_srch, 1', alloc_stat)        
        allocate(self%srch_order(self%nrefs), stat=alloc_stat)
        call alloc_err('In: new; simple_prime2d_srch, 3', alloc_stat)
        self%srch_order = 0
        self%rt = ran_tabu(self%nrefs)
        ! generate the list of in-plane angles and indices
        dang = twopi/real(self%nrots)
        do i=1,self%nrots
            self%angtab(i) = rad2Deg((i-1)*dang)
        end do
        ! initialize best
        self%best_ref%class     = 0
        self%best_ref%inpl_corr = -1.
        self%best_ref%inpl_ind  = 0
        ! the instance now exists
        self%exists = .true.
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::CONSTRUCTED NEW SIMPLE_prime2D_srch OBJECT'
    end subroutine
    
    ! CALCULATORS
    
    !>  \brief calculate the in-plane rotational index for the rot in-plane angle
    function calc_roind( self, rot ) result( ind )
        class(prime2D_srch), intent(in) :: self
        real, intent(in)               :: rot
        integer :: ind, i, loc(1)
        real    :: dists(self%nrots)
        do i=1,self%nrots
            dists(i) = sqrt((self%angtab(i)-rot)**2.)
        end do
        loc = minloc(dists)
        ind = loc(1)
    end function
    
    !>  \brief calculate the in-plane angle of in-plane rotational index rot
    function calc_rot( self, ind ) result( rot )
        class(prime2D_srch), intent(in) :: self
        integer, intent(in) :: ind
        real :: rot
        if( ind < 1 .or. ind > self%nrots )then
            write(*,*) 'rotational index is: ', ind, ' which is out of range; calc_rot; simple_prime2D_srch'
            stop
        endif
        rot = self%angtab(ind)
    end function
    
    ! GETTERS & SETTERS
    
    !>  \brief  to get the class
    subroutine get_cls( self, o )
        use simple_ori,  only: ori
        class(prime2D_srch), intent(in) :: self
        class(ori), intent(inout)       :: o
        type(ori)                       :: o_copy
        real                            :: euls(3), mi, frac
        integer                         :: class, rot
        real                            :: x, y, corr
        ! copy the input orientation
        o_copy = o
        ! get new indices
        class = self%best_ref%class
        rot   = self%best_ref%inpl_ind
        ! get in-plane angle
        euls = 0.
        euls(3) = 360.-self%calc_rot(rot) ! change sgn to fit convention
        if( euls(3) == 360. ) euls(3) = 0.
        call o%set_euler(euls)
        ! calculate overlap between distributions
        mi = 0.
        if( self%class == class ) mi = mi+1
        if( self%rot   == rot   ) mi = mi+1
        mi = mi/2.
        ! set parameters
        if( self%doshift )then
            ! shifts must be obtained by vector addition
            x = self%shvec(1)+self%best_ref%shvec(1)
            y = self%shvec(2)+self%best_ref%shvec(2)
        else
            x = 0.
            y = 0.
        endif
        call o%set('x',x)
        call o%set('y',y)
        corr = self%best_ref%inpl_corr
        frac = 100.*(real(self%nrefs_eval)/real(self%nrefs))
        call o%set_list(['class  ','corr   ','mi_hard','frac   '], [real(class),corr,mi,frac])
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::GOT BEST ORI'
    end subroutine
    
    ! SEARCH ROUTINES
    
    !>  \brief a master prime search routine
    subroutine exec_prime2D_srch( self, pftcc, iptcl, lp, o )
        use simple_ori, only: ori
        class(prime2D_srch), intent(inout)     :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        real, intent(in)                       :: lp
        class(ori), intent(inout), optional    :: o
        call self%prep4srch(pftcc, iptcl, lp, o)
        if( .not. present(o) )then
            call self%greedy_srch(pftcc, iptcl)
        else if( self%refine .eq. 'greedy' )then
            call self%greedy_srch(pftcc, iptcl)
        else
            call self%stochastic_srch(pftcc, iptcl)
        endif
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::EXECUTED PRIME2D_SRCH'
    end subroutine
    
    !>  \brief  prepares for the search
    subroutine prep4srch( self, pftcc, iptcl, lp, o_prev )
        use simple_ori, only: ori
        class(prime2D_srch), intent(inout)     :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        real, intent(in)                       :: lp
        class(ori), intent(inout), optional    :: o_prev
        real    :: lims(2,2)
        ! make random reference direction order
        call self%rt%reset
        call self%rt%ne_ran_iarr(self%srch_order)
        ! initialize origin shift search object
        lims(1,1) = -self%trs
        lims(1,2) = self%trs
        lims(2,1) = -self%trs
        lims(2,2) = self%trs
        ! we are using simplex with barrier constraint by default
        call pftcc_shsrch_init(pftcc, 'simplex', lims, 3)
        if( present(o_prev) )then
            ! find previous discrete alignment parameters
            self%class = nint(o_prev%get('class'))             ! class index
            self%rot   = self%calc_roind(360.-o_prev%e3get())  ! in-plane angle index
            self%ref   = self%class                            ! reference index
            if( self%doshift )then
                self%shvec = [o_prev%get('x'),o_prev%get('y')] ! shift vector
            else
                self%shvec = 0.
            endif
            ! calculate previous best corr (treshold for better)
            self%corr_prev = pftcc%corr(self%ref, iptcl, self%rot)
            ! put prev_best last to avoid cycling
            call put_last(self%class , self%srch_order)
        else
            self%class     = 1
            self%rot       = 1
            self%ref       = 1
            self%shvec     = 0.
            self%corr_prev = 1.
        endif
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::PREPARED FOR SIMPLE_prime2D_srch'
    end subroutine
    
    !>  \brief  executes the greedy rotational search
    subroutine greedy_srch( self, pftcc, iptcl )
        class(prime2D_srch), intent(inout)     :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        integer :: ref, rot_new, loc(1)
        real    :: corr_new, cxy(3), corrs(self%nrots)
        self%corr_prev = -1.
        do ref=1,self%nrefs
            corrs = pftcc%gencorrs(ref, iptcl)
            loc = maxloc(corrs)
            if( corrs(loc(1)) >= self%corr_prev )then
                ! update the class
                self%best_ref%class     = ref
                ! update the correlations
                self%best_ref%inpl_corr = corrs(loc(1))
                self%corr_prev = self%best_ref%inpl_corr
                ! update the in-plane angle
                self%best_ref%inpl_ind  = loc(1)
            endif
        end do
        ! we always evaluate all references using the greedy approach
        self%nrefs_eval = self%nrefs
        ! search shifts
        if( self%doshift )then
            call pftcc_shsrch_set_indices(self%best_ref%class, iptcl, self%best_ref%inpl_ind)
            cxy = pftcc_shsrch_minimize()
            self%best_ref%shvec = cxy(2:3)
            self%best_ref%inpl_corr = cxy(1)
        else
            self%best_ref%shvec = [0.,0.]
        endif
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::FINISHED GREEDY SEARCH'
    end subroutine
    
    !>  \brief  executes the stochastic rotational search
    subroutine stochastic_srch( self, pftcc, iptcl )
        use simple_rnd, only: shcloc
        class(prime2D_srch), intent(inout)     :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        integer :: i, ref, rot_new, this ! loc(1)
        real    :: corr_new, cxy(3), corrs(self%nrots), diffs(self%nrots)
        logical :: found_better, available(self%nrots)
        ! initialize
        self%nrefs_eval = 0
        ! search
        found_better = .false.
        do i=1,self%nrefs
            ref = self%srch_order(i)
            if( self%doinpl )then
                corrs = pftcc%gencorrs(ref, iptcl)
                ! keep track of how many references we are evaluating
                self%nrefs_eval = self%nrefs_eval+1
                this = shcloc(corrs, self%corr_prev)
                if( this > 0 )then
                    ! update the class
                    self%best_ref%class = ref
                    ! update the correlation
                    self%best_ref%inpl_corr = corrs(this)
                    ! update the in-plane angle
                    self%best_ref%inpl_ind  = this
                    ! indicate that we found a better solution
                    found_better = .true.
                    exit ! first-improvement heuristic
                endif
            else
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
            ! search shifts
            if( self%doshift )then
                call pftcc_shsrch_set_indices(self%best_ref%class, iptcl, self%best_ref%inpl_ind)
                cxy = pftcc_shsrch_minimize()
                self%best_ref%shvec = cxy(2:3)
                self%best_ref%inpl_corr = cxy(1)
            else
                self%best_ref%shvec = [0.,0.]
            endif
        else
            self%best_ref%shvec = [0.,0.]
        endif
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::FINISHED STOCHASTIC SEARCH'
    end subroutine
    
    !>  \brief  executes shift-only search
    subroutine shift_srch( self, pftcc, iptcl )
        class(prime2D_srch), intent(inout)     :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        real    :: cxy(3)
        call pftcc_shsrch_set_indices(self%class, iptcl, self%rot)
        cxy = pftcc_shsrch_minimize()
        self%best_ref%class     = self%class
        self%best_ref%inpl_corr = cxy(1)
        self%best_ref%inpl_ind  = self%rot
        self%best_ref%shvec     = cxy(2:3)
        if( debug ) write(*,'(A)') '>>> PRIME2D_SRCH::FINISHED SHIFT SEARCH'
    end subroutine

    !>  \brief  is a destructor
    subroutine kill( self )
        class(prime2D_srch), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%angtab, self%srch_order)
            call self%rt%kill
            self%exists = .false.
        endif
    end subroutine

end module simple_prime2D_srch