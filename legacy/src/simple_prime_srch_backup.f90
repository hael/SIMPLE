module simple_prime_srch
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_oris,             only: oris
use simple_ran_tabu,         only: ran_tabu
use simple_jiffys,           only: alloc_err
use simple_math,             only: hpsort, rad2deg, round2even, put_last
use simple_pftcc_shsrch      ! singleton
use simple_defs              ! singleton
implicit none

public :: prime_srch
private

logical, parameter :: debug=.false.

type reference
    integer :: proj  = 1     !< projection direction index
    integer :: state = 1     !< state of reference (conformational or compositional)
    real    :: inpl_corr     !< best in-plane correlation
    integer :: inpl_ind      !< best in-plane index
    real    :: shvec(2) = 0. !< shift vector
    real    :: w = 0.        !< orientation weight
end type

type prime_srch
    private
    integer                      :: nrefs      = 0        !< total number of references (nstates*nprojs)
    integer                      :: nstates    = 0        !< number of states
    integer                      :: nprojs     = 0        !< number of projections (same as number of oris in o_refs)
    integer                      :: nrots      = 0        !< number of in-plane rotations in polar representation
    integer                      :: npeaks     = 0        !< number of peaks (nonzero orientation weights)
    integer                      :: nbetter    = 0        !< nr of better orientations identified
    integer                      :: nrefs_eval = 0        !< nr of references evaluated
    integer                      :: proj       = 0        !< previous projection index
    integer                      :: rot        = 0        !< previous in-plane rotation index
    integer                      :: state      = 0        !< previous state index
    integer                      :: ref        = 0        !< previous reference index
    integer                      :: nnn        = 0        !< number of nearest neighbors
    integer                      :: fromp      = 0        !< from particle index
    real                         :: atres      = 0.       !< angular treshold
    real                         :: trs        = 0.       !< shift range parameter [-trs,trs]
    real                         :: shvec(2)   = 0.       !< previous origin shift vector
    real                         :: corr_prev  = 1.       !< previous best correlation
    real                         :: eps        = 0.5      !< learning rate for the state probabilities
    class(oris), pointer         :: o_refs => null()      !< pointer to orientations of references in proj_space
    type(reference), allocatable :: proj_space(:)         !< array of references (the space of projection directions)
    integer, allocatable         :: refinds(:,:)          !< table of reference indices
    integer, allocatable         :: srch_order(:)         !< stochastoc search order
    integer, allocatable         :: proj_space_inds(:)    !< projection space index array
    integer, allocatable         :: inplmat(:,:)          !< in-plane rotation index matrix
    real,    allocatable         :: proj_space_corrs(:)   !< projection space correlations (in-plane maxima)
    real,    allocatable         :: angtab(:)             !< list of in-plane angles
    real,    allocatable         :: state_probs(:)        !< state probabilities
    real,    allocatable         :: state_probs_prev(:)   !< previous state probabilities
    real,    allocatable         :: corrmat(:,:)          !< correlation matrix
    type(ran_tabu)               :: rt                    !< object for random number generation
    character(len=STDLEN)        :: refine                !< refinement flag
    character(len=STDLEN)        :: diversify             !< diversification flag
    character(len=STDLEN)        :: test                  !< test flag
    logical                      :: doshift = .true.      !< origin shift search indicator
    logical                      :: exists  = .false.     !< 2 indicate existence
  contains
    ! CONSTRUCTOR & INITIALIZER
    procedure :: new
    procedure :: init
    ! CALCULATORS
    procedure, private :: calc_roind
    procedure, private :: calc_rot
    ! GETTERS & SETTERS
    procedure :: get_ori_best
    procedure :: get_state
    procedure :: get_proj
    procedure :: get_class
    procedure :: get_ori
    procedure :: get_euler
    procedure :: ang_sdev
    ! SEARCH ROUTINES & WEIGHT CALCULATION
    procedure :: exec_prime_srch
    procedure :: exec_shift_srch
    procedure, private :: prep4srch
    procedure :: gencorrs_all
    procedure, private :: stochastic_srch
    procedure, private :: shift_srch
    procedure, private :: stochastic_weights
    ! DESTRUCTOR
    procedure :: kill
end type

contains

    !>  \brief  is a constructor
    subroutine new( self, e, p )
        use simple_params, only: params
        class(prime_srch), intent(inout) :: self !< instance
        class(params), intent(in)        :: p    !< parameters
        class(oris), intent(in), target  :: e    !< Euler angles
        integer :: alloc_stat, cnt, i, s
        real    :: dang, ang, x
        ! destroy possibly pre-existing instance
        call self%kill
        ! set constants
        self%nrefs      = p%nspace*p%nstates
        self%nstates    = p%nstates
        self%nprojs     = p%nspace
        self%nrots      = round2even(twopi*real(p%ring2))
        self%npeaks     = p%npeaks
        self%atres      = p%tres
        self%nbetter    = 0
        self%nrefs_eval = 0
        self%trs        = p%trs
        self%o_refs     => e
        self%doshift    = p%doshift
        self%refine     = p%refine
        self%diversify  = p%diversify
        self%nnn        = p%nnn
        self%test       = p%test
        self%fromp      = p%fromp
        ! construct composites
        allocate( self%proj_space(self%nrefs), self%refinds(self%nstates,self%nprojs),&
                  self%proj_space_inds(self%nrefs), self%proj_space_corrs(self%nrefs),&
                  self%state_probs(self%nstates), self%state_probs_prev(self%nstates),&
                  self%angtab(self%nrots), stat=alloc_stat )
        call alloc_err('In: new; simple_prime_srch, 1', alloc_stat)
        self%proj_space_inds  = 0
        self%proj_space_corrs = -1.
        self%state_probs      = 1./real(self%nstates)
        self%state_probs_prev = 1./real(self%nstates)
        if( self%refine .eq. 'neigh' )then ! local refinement
            allocate(self%srch_order(self%nnn), stat=alloc_stat)
            call alloc_err('In: new; simple_prime_srch, 2', alloc_stat)
            self%srch_order = 0
            self%rt = ran_tabu(self%nnn)
        else
            allocate(self%srch_order(self%nrefs), stat=alloc_stat)
            call alloc_err('In: new; simple_prime_srch, 3', alloc_stat)
            self%srch_order = 0
            self%rt = ran_tabu(self%nrefs)
        endif
        if( self%test .eq. 'yes' )then
            if( p%top-p%fromp+1 == self%nrefs )then
                allocate( self%corrmat(self%nrefs,self%nrefs), self%inplmat(self%nrefs,self%nrefs), stat=alloc_stat )
                call alloc_err('In: new; simple_prime_srch, 4', alloc_stat)
                self%corrmat = 0.
                self%inplmat = 0
            else
                stop 'nonconforming number of particles for this mode of execution; new; simple_prime_srch'
            endif
        endif
        ! generate the list of in-plane angles and indices
        dang = twopi/real(self%nrots)
        do i=1,self%nrots
            self%angtab(i) = rad2deg((i-1)*dang)
        end do
        ! initialize structures
        call self%init
        self%exists = .true.
        if( debug ) write(*,'(A)') '>>> PRIME_SRCH::CONSTRUCTED NEW SIMPLE_PRIME_SRCH OBJECT'
    end subroutine
    
    !>  \brief  is an initializer
    subroutine init( self )
        class(prime_srch), intent(inout) :: self !< instance
        integer :: cnt, i, s
        ! create index transformation & output data (correlation) storage structure
        cnt = 0
        do s=1,self%nstates
            do i=1,self%nprojs
                cnt = cnt+1
                self%proj_space(cnt)%state = s
                self%proj_space(cnt)%proj  = i
                self%refinds(s,i)          = cnt
                self%proj_space(cnt)%inpl_corr = -1.
                self%proj_space(cnt)%inpl_ind  = 0
            end do
        end do
        if( debug ) write(*,'(A)') '>>> PRIME_SRCH::INITIALISED SIMPLE_PRIME_SRCH OBJECT'
    end subroutine
    
    ! CALCULATORS
    
    !>  \brief calculate the in-plane rotational index for the rot in-plane angle
    function calc_roind( self, rot ) result( ind )
        class(prime_srch), intent(in) :: self
        real, intent(in)              :: rot
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
        class(prime_srch), intent(in) :: self
        integer, intent(in) :: ind
        real :: rot
        if( ind < 1 .or. ind > self%nrots )then
            write(*,*) 'rotational index is: ', ind, ' which is out of range; calc_rot; simple_prime_srch'
            stop
        endif
        rot = self%angtab(ind)
    end function
    
    ! GETTERS & SETTERS
    
    !>  \brief  to get the best orientation from the discrete space
    subroutine get_ori_best( self, o )
        use simple_ori,  only: ori
        class(prime_srch), intent(in) :: self
        class(ori), intent(inout)     :: o
        type(ori)                     :: o_copy
        real                          :: euls(3), mi, frac, dist
        character(len=STDLEN)         :: dig, key
        integer                       :: proj_new, rot_new, state_new, class, s, ref
        real                          :: x, y, corr, w
        ! copy the input orientation
        o_copy = o
        ! set the reference index
        ref = self%proj_space_inds(self%nrefs)
        if( self%refine .eq. 'shift' )then
            ! shifts must be obtained by vector addition
            x = self%shvec(1)+self%proj_space(ref)%shvec(1)
            y = self%shvec(2)+self%proj_space(ref)%shvec(2)
            call o%set('x',x)
            call o%set('y',y)
            corr = self%proj_space(ref)%inpl_corr
            call o%set('corr', corr)
        else
            ! get new indices
            proj_new  = self%proj_space(self%proj_space_inds(self%nrefs))%proj
            rot_new   = self%proj_space(self%proj_space_inds(self%nrefs))%inpl_ind
            state_new = self%proj_space(self%proj_space_inds(self%nrefs))%state
            euls(1)   = self%o_refs%e1get(proj_new)
            euls(2)   = self%o_refs%e2get(proj_new)
            if( rot_new < 1 .or. rot_new > self%nrots )then
                write(*,*) 'rot_new  in simple_prime_srch::get_ori is: ', rot_new, ' which is out of range'
                stop
            endif
            euls(3) = 360.-self%calc_rot(rot_new) ! change sgn to fit convention
            if( euls(3) == 360. ) euls(3) = 0.
            class = self%proj_space_inds(self%nrefs)
            ! calculate overlap between distributions
            mi = 0.
            if( self%proj  == proj_new  ) mi = mi+1
            if( self%rot   == rot_new   ) mi = mi+1
            if( self%state == state_new ) mi = mi+1
            mi = mi/3.
            ! set parameters
            call o%set_euler(euls)
            if( self%doshift )then
                ! shifts must be obtained by vector addition
                x = self%shvec(1)+self%proj_space(self%proj_space_inds(self%nrefs))%shvec(1)
                y = self%shvec(2)+self%proj_space(self%proj_space_inds(self%nrefs))%shvec(2)
            else
                x = 0.
                y = 0.
            endif
            call o%set('x',x)
            call o%set('y',y)
            corr = self%proj_space(self%proj_space_inds(self%nrefs))%inpl_corr
            if( self%refine .eq. 'neigh' )then
                frac = 100.*(real(self%nrefs_eval)/real(self%nnn))
            else
                frac = 100.*(real(self%nrefs_eval)/real(self%nrefs))
            endif
            dist = 0.5*rad2deg(o_copy.euldist.o)+0.5*o_copy%get('dist')
            w    = self%proj_space(self%proj_space_inds(self%nrefs))%w
            call o%set_list(['state  ','class  ','corr   ','dist   ','mi_hard','frac   ', 'ow     '],&
            [real(state_new),real(class),corr,dist,mi,frac,w])
            ! output multinomal state distribution
            if( self%nstates > 1 )then
                do s=1,self%nstates
                    write(dig,*) s
                    key = 'sprob'//trim(adjustl(dig))
                    call o%set(trim(adjustl(key)), self%state_probs(s))
                end do
            endif
        endif
        if( debug ) write(*,'(A)') '>>> PRIME_SRCH::GOT BEST ORI'
    end subroutine
    
    !>  \brief  get the state corresponding to class
    function get_state( self, class ) result( s )
        class(prime_srch), intent(inout) :: self
        integer, intent(in)              :: class
        integer                          :: s
        s = self%proj_space(class)%state
    end function
    
    !>  \brief  get the projection direction corresponding to class
    function get_proj( self, class ) result( p )
        class(prime_srch), intent(inout) :: self
        integer, intent(in)              :: class
        integer                          :: p
        p = self%proj_space(class)%proj
    end function
    
    !>  \brief  get the projection direction corresponding to class
    function get_class( self, state, proj ) result( class )
        class(prime_srch), intent(inout) :: self
        integer, intent(in)              :: state, proj
        integer                          :: class
        class = self%refinds(state,proj)
    end function
    
    !>  \brief  to get one orientation from the discrete space
    subroutine get_ori( self, i, o )
        use simple_ori, only: ori
        class(prime_srch), intent(inout) :: self
        integer, intent(in)              :: i
        type(ori), intent(inout)         :: o
        real                             :: euls(3), x, y
        integer                          :: ref
        ref = self%proj_space_inds(self%nrefs-i+1)
        if( self%doshift )then
            ! shifts must be obtained by vector addition
            x = self%shvec(1)+self%proj_space(ref)%shvec(1)
            y = self%shvec(2)+self%proj_space(ref)%shvec(2)
            if( abs(x) < 1e-6 ) x = 0.
            if( abs(y) < 1e-6 ) y = 0.
            call o%set('x',x)
            call o%set('y',y)
        else
            call o%set('x',0.)
            call o%set('y',0.)
        endif
        euls(1) = self%o_refs%e1get(self%proj_space(ref)%proj)
        euls(2) = self%o_refs%e2get(self%proj_space(ref)%proj)
        euls(3) = 360.-self%calc_rot(self%proj_space(ref)%inpl_ind) ! change sgn to fit convention
        call o%set_euler(euls)
        call o%set('state',real(self%proj_space(ref)%state))
        call o%set('corr',self%proj_space(ref)%inpl_corr)
        call o%set('ow',self%proj_space(ref)%w)
    end subroutine
    
    !>  \brief  to get one orientation from the discrete space
    subroutine get_euler( self, i, o )
        use simple_ori, only: ori
        class(prime_srch), intent(inout) :: self
        integer, intent(in)              :: i
        type(ori), intent(inout)         :: o
        real                             :: euls(3)
        integer                          :: ref
        ref = self%proj_space_inds(self%nrefs-i+1)
        euls(1) = self%o_refs%e1get(self%proj_space(ref)%proj)
        euls(2) = self%o_refs%e2get(self%proj_space(ref)%proj)
        euls(3) = 360.-self%calc_rot(self%proj_space(ref)%inpl_ind) ! change sgn to fit convention
        call o%set_euler(euls)
        call o%set('state',real(self%proj_space(ref)%state))
        call o%set('corr',self%proj_space(ref)%inpl_corr)
    end subroutine
    
    !>  \brief  standard deviation
    function ang_sdev( self ) result( sdev )
        use simple_ori, only: ori
        class(prime_srch), intent(inout) :: self
        real                             :: sdev
        integer                          :: nstates, s, i, cnt, alloc_stat, ref
        type(reference)                  :: best_ones(self%npeaks)
        cnt = 0
        do i=self%nrefs,self%nrefs-self%npeaks+1,-1
            cnt = cnt+1
            ref = self%proj_space_inds(i)
            best_ones(cnt) = self%proj_space(ref)   
        end do
        nstates = maxval(best_ones%state)
        sdev = 0.
        do s=1,nstates
            sdev = sdev+ang_sdev_state()
        end do
        sdev = sdev/real(nstates)
        cnt = 0
        if( debug ) write(*,'(A)') '>>> PRIME_SRCH::CALCULATED ANG_SDEV'
    
        contains

            function ang_sdev_state( ) result( sdev )
                use simple_stat, only: moment
                real                        :: ave, sdev, var
                type(ori)                   :: o_best, o
                integer                     :: loc(1), alloc_stat, i, cnt
                real, allocatable           :: dists(:)
                logical                     :: err
                real                        :: w
                call o_best%new
                call o%new
                ! count the number of occurances of state
                cnt = 0
                do i=1,size(best_ones)
                    if( best_ones(i)%state == s ) cnt = cnt+1
                end do
                if( cnt <= 3 ) then ! because one is excluded in the next step & moment needs at least 2 objs
                    sdev = 0.
                    return
                endif
                allocate( dists(cnt-1), stat=alloc_stat )
                call alloc_err( 'ang_sdev_state; simple_prime_srch', alloc_stat )
                ! copy data
                cnt = 0
                w = 0.
                loc(1) = 0
                do i=1,size(best_ones)
                    if( best_ones(i)%state == s )then
                        cnt = cnt+1
                        if( best_ones(i)%w > w )then
                            loc(1) = i
                            w = best_ones(i)%w
                        endif
                    endif
                end do
                if( loc(1) == 0 )then
                    sdev = 0.
                    return
                else
                    call self%get_euler(loc(1), o_best)
                    cnt = 0
                    do i=1,size(best_ones)
                        if( i == loc(1) ) cycle
                        if( best_ones(i)%state == s )then
                            cnt = cnt+1
                            call self%get_euler(i, o)
                            dists(cnt) = rad2deg(o_best.euldist.o) 
                        endif
                    end do
                    call moment(dists(:cnt), ave, sdev, var, err)
                endif
                deallocate(dists)
                call o_best%kill
                call o%kill
            end function
            
    end function

    ! SEARCH ROUTINES & WEIGHT CALCULATION
    
    !>  \brief a master prime search routine
    subroutine exec_prime_srch( self, pftcc, iptcl, lp, wcorr, o, nnmat )
        use simple_ori, only: ori
        class(prime_srch), intent(inout)       :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        real, intent(in)                       :: lp
        real, intent(out)                      :: wcorr
        class(ori), intent(inout), optional    :: o
        integer, intent(in), optional          :: nnmat(self%nrefs,self%nnn)
        call self%prep4srch(pftcc, iptcl, lp, o)
        call self%stochastic_srch(pftcc, iptcl, nnmat)
        call self%stochastic_weights(wcorr)
        if( debug ) write(*,'(A)') '>>> PRIME_SRCH::EXECUTED PRIME_SRCH'
    end subroutine
    
    !>  \brief a master prime search routine
    subroutine exec_shift_srch( self, pftcc, iptcl, lp, o )
        use simple_ori, only: ori
        class(prime_srch), intent(inout)       :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        real, intent(in)                       :: lp
        class(ori), intent(inout)              :: o
        call self%prep4srch(pftcc, iptcl, lp, o)
        call self%shift_srch(pftcc, iptcl)
        if( debug ) write(*,'(A)') '>>> PRIME_SRCH::EXECUTED SHIFT_SRCH'
    end subroutine

    !>  \brief  prepares for the search
    subroutine prep4srch( self, pftcc, iptcl, lp, o_prev )
        use simple_ori, only: ori
        class(prime_srch), intent(inout)       :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        real, intent(in)                       :: lp
        class(ori), intent(inout), optional    :: o_prev
        character(len=STDLEN) :: dig, key  
        integer :: s
        logical :: did_set_stateprobs, did_set_corr_prev
        real    :: lims(2,2), cc_t, cc_t_min_1
        ! make 'random' projection direction order
        call self%rt%reset
        call self%rt%ne_ran_iarr(self%srch_order)
        ! initialize origin shift search object
        lims(1,1) = -self%trs
        lims(1,2) = self%trs
        lims(2,1) = -self%trs
        lims(2,2) = self%trs
        ! we are using simplex with barrier constraint by default
        ! need to look into the out of bounds problems with BFGS before making it the standard
        call pftcc_shsrch_init(pftcc, 'simplex', lims, 1)
        if( present(o_prev) )then
            ! find previous discrete alignment parameters
            self%proj  = self%o_refs%find_closest_proj(o_prev) ! projection index
            self%rot   = self%calc_roind(360.-o_prev%e3get())  ! in-plane angle index
            self%state = nint(o_prev%get('state'))             ! state index
            self%ref   = self%refinds(self%state,self%proj)    ! reference index
            if( self%doshift )then
                self%shvec = [o_prev%get('x'),o_prev%get('y')] ! shift vector
            else
                self%shvec = 0.
            endif
            if( self%state > self%nstates )then
                stop 'previous best state outside boundary; prep4srch; simple_prime_srch'
            endif
            ! put prev_best last to avoid cycling
            call put_last(self%ref, self%srch_order)
            ! calculate previous best corr (treshold for better)
            did_set_corr_prev = .false.
            cc_t = pftcc%corr(self%ref, iptcl, self%rot)
            if( self%diversify .eq. 'yes' )then        
                cc_t_min_1 = -1.
                if( o_prev%isthere('corr') ) cc_t_min_1 = o_prev%get('corr')
                if( o_prev%isthere('lp') )then
                    if( abs(o_prev%get('lp')-lp) < 0.5 )then ! previous and present correlation values are comparable
                        if( cc_t_min_1 > 0. )then
                            self%corr_prev = 0.5*cc_t+0.5*cc_t_min_1 ! moving average correlation (diversifying limit)
                            did_set_corr_prev = .true.
                        endif
                    endif
                endif
                if( .not. did_set_corr_prev ) self%corr_prev = cc_t ! greedy limit
            else
                self%corr_prev = cc_t
            endif
        else
            self%proj      = 1
            self%rot       = 1
            self%state     = 1
            self%ref       = 1
            self%shvec     = 0.
            self%corr_prev = 1.
        endif
        ! check for presence of state probabilities and 
        ! extract the multimodal distribution
        did_set_stateprobs = .false.
        if( present(o_prev) )then
            if( o_prev%isthere('sprob1') )then
                do s=1,self%nstates
                    write(dig,*) s
                    key = 'sprob'//trim(adjustl(dig))
                    self%state_probs_prev(s) = o_prev%get(key)
                end do
                did_set_stateprobs = .true.
            endif
        endif
        if( .not. did_set_stateprobs ) self%state_probs_prev = 1./real(self%nstates)
        if( debug ) write(*,'(A)') '>>> PRIME_SRCH::PREPARED FOR SIMPLE_PRIME_SRCH'
    end subroutine
    
    !>  \brief  generates all correlations
    subroutine gencorrs_all( self, pftcc )
        class(prime_srch), intent(inout)       :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        call pftcc%expand_dim
        call pftcc%gencorrs_all(self%corrmat, self%inplmat)
    end subroutine
    
    !>  \brief  executes the stochastic rotational search
    subroutine stochastic_srch( self, pftcc, iptcl, nnmat )
        class(prime_srch), intent(inout)       :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        integer, intent(in), optional          :: nnmat(self%nrefs,self%nnn)
        integer :: i, j, loc(1), ref, endit, ptcl_ind
        real    :: corr, cxy(3), corrs(self%nrots)
        ! initialize
        self%nbetter          = 0
        self%proj_space_corrs = -1.
        self%nrefs_eval       = 0
        self%proj_space_inds  = 0
        ! search
        endit = self%nrefs
        if( present(nnmat) ) endit = self%nnn
        do i=1,endit
            ! set the stochastic reference index
            if( present(nnmat) )then
                ref = nnmat(self%proj,self%srch_order(i))
                if( allocated(self%corrmat) )then
                    stop 'cannot use neigh for this mode of execution; stochastic_srch; simple_prime_srch'
                endif
            else  
                ref = self%srch_order(i)
            endif
            ! stash the stochastic reference index
            self%proj_space_inds(ref) = ref
            ! calculate all in-plane correlations
            if( allocated(self%corrmat) )then
                ptcl_ind = iptcl-self%fromp+1
                self%proj_space(ref)%inpl_ind  = self%inplmat(ptcl_ind,ref)
                self%proj_space(ref)%inpl_corr = self%corrmat(ptcl_ind,ref)
            else
                ! we don't use shifts here because the particle images are pre-shifted
                ! when the pftcc object is filled-up
                corrs = pftcc%gencorrs(ref, iptcl)
                loc   = maxloc(corrs)
                self%proj_space(ref)%inpl_ind  = loc(1)
                self%proj_space(ref)%inpl_corr = corrs(loc(1))
            endif
            ! update the correlation
            self%proj_space_corrs(ref) = self%proj_space(ref)%inpl_corr
            ! update nbetter to keep track of how many improving solutions we have identified
            if( self%npeaks == 1 )then
                if( self%proj_space_corrs(ref) > self%corr_prev ) self%nbetter = self%nbetter+1
            else
                if( self%proj_space_corrs(ref) >= self%corr_prev ) self%nbetter = self%nbetter+1
            endif
            ! keep track of how many references we are evaluating
            self%nrefs_eval = self%nrefs_eval+1
            ! exit condition
            if( self%nbetter == self%npeaks ) exit
        end do
        ! sort in projection direction space (the max inpl corrs)
        call hpsort(self%nrefs, self%proj_space_corrs, self%proj_space_inds)
        ! search shifts
        do i=self%nrefs,self%nrefs-self%npeaks+1,-1
            ref = self%proj_space_inds(i)
            if( self%doshift )then
                call pftcc_shsrch_set_indices(ref, iptcl, self%proj_space(ref)%inpl_ind)
                cxy = pftcc_shsrch_minimize()
                self%proj_space(ref)%shvec = cxy(2:3)
                self%proj_space(ref)%inpl_corr = cxy(1)
            else
                self%proj_space(ref)%shvec = [0.,0.]
            endif
        end do
        if( debug ) write(*,'(A)') '>>> PRIME_SRCH::FINISHED STOCHASTIC SEARCH'
    end subroutine
    
    !>  \brief  executes shift-only search
    subroutine shift_srch( self, pftcc, iptcl )
        class(prime_srch), intent(inout)       :: self
        class(polarft_corrcalc), intent(inout) :: pftcc
        integer, intent(in)                    :: iptcl
        integer :: i, j, loc(1), ref
        real    :: cxy(3)
        do i=1,self%nrefs
            self%proj_space_inds(i) = i
        end do
        call pftcc_shsrch_set_indices(self%ref, iptcl, self%rot)
        cxy = pftcc_shsrch_minimize()
        self%proj_space(self%proj_space_inds(self%nrefs))%shvec = cxy(2:3)
        self%proj_space(self%proj_space_inds(self%nrefs))%inpl_corr = cxy(1)
        if( debug ) write(*,'(A)') '>>> PRIME_SRCH::FINISHED SHIFT SEARCH'
    end subroutine
    
    !>  \brief  calculates stochastic weights
    subroutine stochastic_weights( self, wcorr )
        class(prime_srch), intent(inout) :: self
        real, intent(out)                :: wcorr
        real                             :: wsum, spsum
        integer                          :: i, j, nbetter_inpl, alloc_stat, ref
        ! update parameters
        self%state_probs = 0.
        do i=self%nrefs,self%nrefs-self%npeaks+1,-1
            ref = self%proj_space_inds(i)
            self%proj_space(ref)%w = self%proj_space(ref)%inpl_corr
            self%state_probs(self%proj_space(ref)%state)= &
            self%state_probs(self%proj_space(ref)%state)+1.
        end do
        ! normalize the state probs
        if( self%nstates == 1 )then
            self%state_probs(1) = 1.
        else
            spsum = sum(self%state_probs)
            self%state_probs = self%state_probs/spsum
            self%state_probs = (1.-self%eps)*self%state_probs_prev+self%eps*(self%state_probs)
        endif
        ! update weights
        wsum = 0.
        do i=self%nrefs,self%nrefs-self%npeaks+1,-1
            ref = self%proj_space_inds(i)
            call weight_update
        end do
        ! finalize weights & calculate weighted correlation
        wcorr = 0.     
        do i=self%nrefs,self%nrefs-self%npeaks+1,-1
            ref = self%proj_space_inds(i)
            self%proj_space(ref)%w = self%proj_space(ref)%w/wsum
            wcorr = wcorr+self%proj_space(ref)%inpl_corr*self%proj_space(ref)%w
        end do
        if( debug ) write(*,'(A)') '>>> PRIME_SRCH :: CALCULATED STOCHASTIC WEIGHTS PRIME_SRCH'
        
        contains
            
            subroutine weight_update
                if( self%proj_space(ref)%inpl_corr > 0. )then
                    self%proj_space(ref)%w =&
                    exp(self%proj_space(ref)%inpl_corr*&
                    self%state_probs(self%proj_space(ref)%state))
                else
                    self%proj_space(ref)%w = 0.
                endif
                wsum = wsum+self%proj_space(ref)%w
            end subroutine
        
    end subroutine
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(prime_srch), intent(inout) :: self !< instance
        if( self%exists )then
            self%o_refs => null()
            deallocate(self%proj_space, self%refinds, self%proj_space_inds,&
            self%proj_space_corrs, self%state_probs, self%state_probs_prev,&
            self%angtab, self%srch_order)
            call self%rt%kill
            if( allocated(self%corrmat) ) deallocate(self%corrmat)
            if( allocated(self%inplmat) ) deallocate(self%inplmat)
            self%exists = .false.
        endif        
    end subroutine

end module simple_prime_srch
