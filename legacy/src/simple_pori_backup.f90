!==Class simple_pori
!
! simple_pori handles orientational information in the probabilistic setting (OASIS). pori extends ori to include
! functionality for online mean and variance estimation as well as importance sampling via a Monte Carlo procedure.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution or 
! modification is regulated by the GNU General Public License. *Author:* Hans Elmlund, 2009-05-26
!
!==Changes are documented below
!* new class, HE 2014-10-05
!
module simple_pori
use simple_ori,        only: ori
use simple_online_var, only: online_var
use simple_sll,        only: sll
use simple_opt_spec,   only: opt_spec
use simple_defs        ! singleton
implicit none

public :: pori, test_pori
private

logical :: report=.false.

type :: pori
    private
    type(ori), pointer            :: optr=>null()     !< pointer to ori obj
    integer                       :: ngaup=0          !< number of Gaussian parameters (5-7) depending on defocus refinement or not
    type(online_var), allocatable :: vars(:)          !< objects for online mean and variance estimation from all local minima
                                                      !! encountered throughout the history of the search
    type(sll)                     :: locopts          !< list of local optimas identified during this round of search
    type(opt_spec)                :: spec_linmin      !< specification of line minimizer
    type(ori), pointer            :: ori_glob=>null() !< pointer to global ori in onflyft_srch, needed to update state                                 
    real, allocatable             :: sdevs(:)         !< the stdevs of the prior Gaussian distributions
    real, allocatable             :: probs(:)         !< the multinomal distribution of probabilities for 1ofK state encoding
    integer(long), allocatable    :: statecnts(:)     !< state counts throughout the history of the search (used to calc state probs)
    real                          :: yb=1.            !< best cost
  contains
    ! CONSTRUCTORS
    procedure          :: new                  
    procedure, private :: serialize                    
    procedure, private :: unserialize                  
    ! GETTERS/SETTERS/ASSIGNERS
    procedure          :: exists     
    procedure          :: get_nbetter
    procedure          :: get_sdevs
    procedure, private :: extract_prev
    procedure, private :: add_minimum                  
    procedure, private :: get_minimum
    procedure          :: minimum2ori                 
    procedure, private :: print_minima                 
    ! I/O
    procedure          :: get_recsz                    
    procedure          :: write                   
    procedure          :: read                  
    ! OASIS IMPORTANCE SAMPLING ROUTINES
    procedure          :: oasis_init                           
    procedure          :: oasis_mc                     
    procedure, private :: update_gaup                          
    procedure, private :: update_probs          
    procedure, private :: update_ori_part           
    ! NLP REFINEMENT
    procedure          :: refine_minimum             
    ! ORISPREAD CALCULATOR
    procedure          :: calc_orispread
    ! WEIGHT CALCULATOR
    procedure          :: calc_weights                            
    ! TOURNAMENT SELECTION 
    procedure          :: binary_tournament
    ! DESTRUCTOR
    procedure, private :: test_alloc_stat              
    procedure          :: kill                  
end type

real, allocatable :: globarr(:)
type(ori)         :: o_old

! SELF%LOCOPTS CONTAINS:
! iarr(1)=state
! rarr(1)=e1
! rarr(2)=e2
! rarr(3)=e3
! rarr(4)=x
! rarr(5)=y
! rarr(self%ngaup-1)=dfx
! rarr(self%ngaup)=dfy
! rarr(self%ngaup+1)=cost
! rarr(self%ngaup+2)=weight

contains

    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    subroutine new( self, ospec, o, ori_glob )
        use simple_jiffys,   only: alloc_err
        use simple_opt_spec, only: opt_spec 
        class(pori), intent(inout)    :: self
        type(ori), intent(in), target :: o, ori_glob
        class(opt_spec), intent(in)   :: ospec
        logical :: cyclic(10)
        integer :: alloc_stat
        call self%kill
        ! set pointer to ori
        self%optr => o
        ! set the number of Gaussian parameters
        self%ngaup = ospec%ndim
        if( self%ngaup < 5 )then
            write(*,*) 'nonconforming ngaup parameter:', self%ngaup, ' in: new; simple_pori'
            stop
        endif
        ! allocate arrays
        if( allocated(globarr) ) deallocate(globarr)
        allocate( globarr(self%ngaup+2), self%vars(self%ngaup), self%sdevs(self%ngaup), stat=alloc_stat)
        call alloc_err('In: new_pori; simple_pori, 1', alloc_stat)
        self%sdevs = 0.
        ! make local minima container
        call self%locopts%new
        if( ospec%nstates > 1 )then ! allocate stateprobs and counters
            allocate(self%probs(ospec%nstates), self%statecnts(ospec%nstates), stat=alloc_stat)
            call alloc_err('In: new_pori; simple_pori, 2', alloc_stat)
            self%probs     = 1./real(ospec%nstates)
            self%statecnts = 0
        endif
        cyclic = .false.
        cyclic(1:3) = .true.
        call self%spec_linmin%specify('linmin', self%ngaup, maxits=ospec%nsample,&
        limits=ospec%limits, cyclic=cyclic(:self%ngaup), ftol=FTOL, nstates=ospec%nstates)
        ! set pointer to global ori
        self%ori_glob => ori_glob 
        if( report ) write(*,*) 'new_pori created'
    end subroutine
    
    !>  \brief  for serialization of a pori object
    function serialize( self ) result( arr )
        use simple_jiffys, only: alloc_err
        ! ngaup             => 1                    => 1
        ! nstates           => 1                    => 2
        ! nsample           => 1                    => 3
        ! limits            => 2*self%ngaup         => 3+1->3+2*ngaup
        ! vars              => 4*self%ngaup         => 3+2*ngaup+1->3+6*ngaup
        ! sdevs             => self%ngaup           => 3+6*ngaup+1->3+7*ngaup 
        ! if( nstates > 1 ) => size(self%probs)     => 3+7*ngaup+1->3+7*ngaup+nstates
        ! if( nstates > 1 ) => size(self%statecnts) => 3+7*ngaup+nstates+1->3+7*ngaup+2*nstates
        class(pori), intent(in) :: self
        real(dp), allocatable   :: arr(:)
        integer :: n, nstates, alloc_stat, cnt, i
        ! nstates:
        if( allocated(self%probs) )then
            nstates = size(self%probs)
        else
            nstates = 0
        endif
        ! ntot:
        n = 3+7*self%ngaup+2*nstates
        if( n < 38 )then
            write(*,*) 'nonconforming ngaup param:', self%ngaup
            write(*,*) 'object subjected to serialization was probably never created; serialize; simple_pori'
            stop
        endif
        allocate( arr(n), stat=alloc_stat )
        call alloc_err('In: serialize; simple_pori', alloc_stat)
        ! ngaup:
        arr(1) = dble(self%ngaup)
        ! nstates:
        arr(2) = dble(nstates)
        ! nsample:
        arr(3) = dble(self%spec_linmin%maxits)
        ! limits:
        cnt = 0
        do i=3+1,3+2*self%ngaup-1,2
            cnt = cnt+1
            arr(i:i+1) = self%spec_linmin%limits(cnt,:)
        end do        
        ! vars:
        cnt = 0
        do i=3+2*self%ngaup+1,3+6*self%ngaup,4
            cnt = cnt+1
            arr(i:i+3) = self%vars(cnt)%serialize()
        end do
        ! sdevs:
        arr(3+6*self%ngaup+1:3+7*self%ngaup) = self%sdevs
        if( nstates > 1 )then
            ! probs:
            arr(3+7*self%ngaup+1:3+7*self%ngaup+nstates)           = dble(self%probs)
            ! statecnts:
            arr(3+7*self%ngaup+nstates+1:3+7*self%ngaup+2*nstates) = dble(self%statecnts)
        endif
         if( report ) write(*,*) 'pori serialized'
    end function
    
    !>  \brief  for re-creation of a pori object from its serialization
    subroutine unserialize( self, arr, o, ori_glob )
        use simple_jiffys,   only: alloc_err
        use simple_opt_spec, only: opt_spec
        class(pori), intent(inout)    :: self
        real(dp), intent(inout)       :: arr(:)
        type(ori), intent(in), target :: o, ori_glob
        real, allocatable             :: limits(:,:)
        type(opt_spec)                :: ospec
        integer :: n, ngaup, i, cnt, nstates, nsample, alloc_stat
        call self%kill
        ngaup   = nint(arr(1))
        nstates = nint(arr(2))
        nsample = nint(arr(3))
        allocate( limits(ngaup,2), stat=alloc_stat )
        call alloc_err('In: unserialize; simple_pori', alloc_stat)
        ! limits:
        cnt = 0
        do i=3+1,3+2*ngaup,2
            cnt = cnt+1
            limits(cnt,:)= real(arr(i:i+1))
        end do
        call ospec%specify('linmin', ngaup, nstates=nstates, nsample=nsample, limits=limits)
        call self%new(ospec, o, ori_glob)
        deallocate(limits)
        ! ntot:
        n = 3+7*self%ngaup+2*nstates
        if( size(arr) /= n )then
            write(*,*) 'nonconforming arr size, expected: ', n, 'got: ', size(arr)
            stop
        endif
        ! vars:
        cnt = 0
        do i=3+2*self%ngaup+1,3+6*self%ngaup,4
            cnt = cnt+1
            call self%vars(cnt)%unserialize(arr(i:i+3))
        end do
        ! sdevs:
        self%sdevs = real(arr(3+6*self%ngaup+1:3+7*self%ngaup))
        if( nstates > 1 )then
            ! probs:
            self%probs     = real(arr(3+7*self%ngaup+1:3+7*self%ngaup+nstates))
            ! statecnts:
            self%statecnts = nint(arr(3+7*self%ngaup+nstates+1:3+7*self%ngaup+2*nstates))
        endif
         if( report ) write(*,*) 'pori unserialized'
    end subroutine
    
    ! GETTERS/SETTERS/ASSIGNERS
    
    !>  \brief  2 check existence
    function exists( self ) result( t )
        class(pori), intent(in) :: self
        logical :: t
        t = allocated(self%vars)
    end function
    
    !>  \brief  returns the number of local optimas identified
    !!          if only one (the prevous one) has been identified
    !!          the assignment should be hard (after further powell
    !!          refinement)
    function get_nbetter( self ) result( nbetter )
        class(pori), intent(in) :: self
        integer :: nbetter
        nbetter = self%locopts%size()
    end function
    
    !>  \brief  for getting the sdevs
    function get_sdevs( self ) result( sdevs )
        class(pori), intent(in) :: self
        real :: sdevs(self%ngaup)
        sdevs = self%sdevs
    end function
    
    !>  \brief  for extracting an array representation of the previous solution
    function extract_prev( self ) result( x )
        use simple_jiffys, only: alloc_err
        class(pori), intent(inout) :: self
        real, allocatable          :: x(:)
        integer                    :: i, alloc_stat, s
        allocate( x(self%ngaup), stat=alloc_stat )
        call alloc_err("In: extract_prev; simple_pori", alloc_stat)
        ! store the previous solution
        o_old = self%optr
        ! extract the previous solution
        do i=1,self%ngaup
            if( i == 1 )then
                x(i) = self%optr%e1get()
            else if( i == 2 )then
                x(i) = self%optr%e2get()
            else if( i == 3 )then
                x(i) = self%optr%e3get()
            else if( i == 4 )then
                if( self%ngaup > 5 )then
                    x(i) = -self%optr%get('x')
                else
                    x(i) = 0.
                endif 
            else if( i == 5 )then
                if( self%ngaup > 5 )then
                    x(i) = -self%optr%get('y')
                else
                    x(i) = 0.
                endif 
            else if( i == 6 )then
                x(i) = self%optr%get('dfx')
            else if( i == 7 )then
                x(i) = self%optr%get('dfy')
            endif
        end do
        call self%optr%set('nbetter', 0.) ! init
        if( report ) write(*,*) 'previous ori params extracted'
    end function
    
    !>  \brief  for adding a new local minimum to the list
    subroutine add_minimum( self, x, yb, w, state )
        class(pori), intent(inout)    :: self
        real, intent(in)              :: x(self%ngaup), yb
        real, intent(in), optional    :: w
        integer, intent(in), optional :: state
        integer :: i, sz
        if( report )then
            sz = size(globarr)
            if( sz < 7 )then
                write(*,*) 'globarr not allocated correctly!, size: ', sz
                write(*,*) 'ngaup: ', self%ngaup
                stop
            endif
        endif
        do i=1,self%ngaup
            globarr(i) = x(i)
        end do
        globarr(self%ngaup+1) = yb
        if( present(w) )then
            globarr(self%ngaup+2) = w
        else
            globarr(self%ngaup+2) = 0. 
        endif
        if( present(state) )then
            call self%locopts%add(rarr=globarr, iarr=[state])
        else
            call self%locopts%add(rarr=globarr)
        endif
        if( report ) write(*,*) 'local minimum added'
    end subroutine
    
    !>  \brief  for printing the local minima
    subroutine print_minima( self )
        class(pori), intent(in) :: self
        call self%locopts%print()
    end subroutine
    
    !>  \brief  for getting a local minimum from the list
    subroutine get_minimum( self, x, yb, w, state, which )
        class(pori), intent(in)        :: self
        real, intent(out)              :: x(self%ngaup), yb
        real, intent(out), optional    :: w
        integer, intent(out), optional :: state
        integer, intent(in), optional  :: which
        integer, allocatable           :: iarr(:)
        integer :: pos, sz
        if( present(which) )then
            if( which >= 1 .and. which <= self%locopts%size() )then
                pos = which
            else
                write(*,*) 'which=', which
                write(*,*) 'locopts%size=',  self%locopts%size()
                stop 'which out of range; get_minimum; simple_pori'
            endif
        else
            pos = self%locopts%size()
        endif
        if( present(state) )then
            call self%locopts%get(pos, rarr=globarr, iarr=iarr)
            state = iarr(1)
        else
            call self%locopts%get(pos, rarr=globarr)
        endif
        if( report )then
            sz = size(globarr)
            if( sz < self%ngaup+2 )then
                write(*,*) 'nonconforming array size of sll extracted globarr:', sz
                stop
            endif 
        endif
        x  = globarr(1:self%ngaup)
        yb = globarr(self%ngaup+1)
        if( present(w) ) w = globarr(self%ngaup+2)
        if( report ) write(*,*) 'local minimum gotten'
    end subroutine
    
    !>  \brief  for getting a local minimum from the list
    subroutine minimum2ori( self, which, o )
        class(pori), intent(in)   :: self
        integer, intent(in)       :: which
        class(ori), intent(inout) :: o
        real    :: x(self%ngaup), yb, w
        integer :: state
        state = 1
        if( allocated(self%probs) )then
            call self%get_minimum(x, yb, w, state, which )
        else
            call self%get_minimum(x, yb, w=w, which=which )
        endif
        if( .not. o%exists() ) call o%new ! create ori 4 output 
        call o%set('corr', -yb)  ! update the correlation
        call o%set('ow', w)      ! update the orientation weight
        call o%set_euler(x(1:3)) ! update the Euler angle
        ! shifts must be obtained by vector addition
        call o%set('x',o_old%get('x')-x(4)) ! update x-shift, revshsgn to fit shift convention
        call o%set('y',o_old%get('y')-x(5)) ! update y-shift, revshsgn to fit shift convention
        call o%set('state',real(state))     ! update state
    end subroutine
    
    ! I/O
    
    !>  \brief  for writing a pori object to file
    subroutine write( self, fname, i )
        use simple_jiffys, only: get_fileunit, fopen_err
        class(pori), intent(in)      :: self
        character(len=*), intent(in) :: fname
        integer, intent(in)          :: i
        real(dp), allocatable        :: arr(:)
        integer                      :: recsz, fnr, file_stat, elemsz
        logical                      :: here
        character(len=3)             :: ostat
        real(dp)                     :: tmp
        arr = self%serialize()
        inquire(iolength=recsz) arr
        inquire(iolength=elemsz) tmp
        fnr = get_fileunit()
        inquire(file=fname, exist=here)
        ostat = 'new'
        if( here ) ostat = 'old'
        open(unit=fnr, status=ostat, action='write', file=fname,&
        access='direct', form='unformatted', recl=recsz, iostat=file_stat)
        call fopen_err( 'In: write; simple_pori', file_stat )
        write(fnr,rec=i) arr
        close(fnr)
    end subroutine

    !>  \brief  for re-creating a pori object from file 
    subroutine read( self, fname, i, o, ori_glob )
        use simple_jiffys, only: get_fileunit, fopen_err, calc_pori_recsz, alloc_err
        class(pori), intent(inout)    :: self
        character(len=*), intent(in)  :: fname
        integer, intent(in)           :: i
        type(ori), intent(in), target :: o, ori_glob
        real(dp), allocatable         :: arr(:)
        integer                       :: recsz, fnr, file_stat, elemsz, alloc_stat
        logical                       :: here
        real(dp)                      :: tmp
        if( allocated(self%probs) )then
            recsz = calc_pori_recsz(self%ngaup,size(self%probs))
        else
            recsz = calc_pori_recsz(self%ngaup,0)
        endif
        recsz = self%get_recsz()
        inquire(iolength=elemsz) tmp
        allocate( arr(recsz/elemsz), stat=alloc_stat )
        call alloc_err('In: read; simple_pori', alloc_stat)
        fnr = get_fileunit()
        inquire(file=fname, exist=here)
        if( here )then
            open(unit=fnr, status='old', action='read', file=fname,&
            access='direct', form='unformatted', recl=recsz, iostat=file_stat)
            call fopen_err( 'In: write; simple_pori', file_stat )
            read(fnr,rec=i) arr
            close(fnr)
            call self%unserialize(arr, o, ori_glob)
        else
            write(*,*) 'file: ', fname, 'does not exist! read; simple_pori'
            stop 
        endif
    end subroutine
    
    !>  \brief  for inquiring the record size of a pori object
    function get_recsz( self ) result( recsz )
         class(pori), intent(inout) :: self
         real(dp), allocatable      :: arr(:)
         integer                    :: recsz
         arr = serialize(self)
         inquire(iolength=recsz) arr
         deallocate(arr)
    end function

    ! OASIS IMPORTANCE SAMPLING ROUTINES
    
    !>  \brief  initialization step for the OASIS minimization
    subroutine oasis_init( self, ospec, sdevs, probs )
        use simple_rnd, only: ran3
        class(pori), intent(inout)  :: self
        class(opt_spec), intent(in) :: ospec
        real, intent(in), optional  :: sdevs(:), probs(:)
        real, allocatable :: x(:)
        logical :: arezero(self%ngaup)
        integer :: i
        ! extract the previous solution
        x = self%extract_prev()
        ! test if best point is set
        arezero = .false.
        do i=1,self%ngaup
            if( x(i) == 0. ) arezero(i) = .true.
        end do
        ! generate initial vector
        if( all(arezero) )then
            do i=1,self%ngaup
                ! initialize each variable by randomized bounds
                x(i) = ospec%limits(i,1)+ran3()*(ospec%limits(i,2)-ospec%limits(i,1))
            end do
            ! update the parameters in ori object
            call self%update_ori_part(x)
        endif
        if( sum(self%sdevs) > 0. )then
            ! sdevs are already read
        else            
            ! initialize by setting all the sdevs to the half the interval
            do i=1,self%ngaup
                self%sdevs(i) = (ospec%limits(i,2)-ospec%limits(i,1))/2.
                if( i > 3 )then
                    if( self%sdevs(i) < 1. ) self%sdevs(i) = 1. 
                endif
            end do
        endif
        if( present(sdevs) )then
            do i=1,size(sdevs)
                self%sdevs(i) = sdevs(i)
            end do
        endif
        if( allocated(self%probs) )then
            if( present(probs) )then
                if( size(self%probs) == size(probs) )then
                    self%probs = probs
                else
                    stop 'size mismatch between probs arrays; oasis_init; simple_pori'
                endif
            else
                if( sum(self%probs) > 0. )then
                    ! probs are already read
                else
                    self%probs = 1./real(size(self%probs))
                endif
            endif
        endif
        if( report ) write(*,*) 'oasis initialized'
    end subroutine
    
    !>  \brief  importance sampling step 4 the OASIS minimization
    subroutine oasis_mc( self, ospec, refine )
        use simple_rnd,      only: gasdev, multinomal
        use simple_opt_subs, only: linmin, check_and_correct_vec
        class(pori), intent(inout)   :: self
        class(opt_spec), intent(in)  :: ospec
        character(len=*), intent(in) :: refine
        real, allocatable :: x(:)
        real    :: y, y2, x_old(self%ngaup)
        logical :: corrected, found_better
        integer :: i, j, state, state_best
        ! extract the previous solution
        x = self%extract_prev()
        ! initialize best-so-far cost
        self%yb = ospec%costfun(x, self%ngaup)
        if( report ) write(*,*) 'cost of previous solution:', self%yb
        ! add the 'old' solution to the list of local minima
        call self%locopts%new
        if( allocated(self%probs) )then
            call self%add_minimum(x, self%yb, state=nint(self%optr%get('state')))
        else  
            call self%add_minimum(x, self%yb)
        endif
        ! associate costfun in line minimizer spec
        call self%spec_linmin%set_costfun(ospec%costfun)
        ! initialize found-better-indicator
        found_better =.false.
        do i=1,self%spec_linmin%maxits                                ! [1,nsample]
            x_old = x                                                 ! store old parameter vector
            do j=1,self%ngaup                                         ! change the components using Gaussian random sampling
                x(j) = gasdev(x_old(j), self%sdevs(j))                ! no limits, since limits give infinite loop at times
            end do
            state = 1
            if( allocated(self%probs) )then                           ! need to sample multinomal state distribution
                state = multinomal(self%probs)                        ! this is the new state
            endif
            call self%ori_glob%set('state', real(state))              ! update ori_glob           
            y = ospec%costfun(x, self%ngaup)                          ! score the new solution vector         
            if( y <= self%yb )then                                    ! update the model if a better solution is found
                found_better = .true.                                 ! indicate better solution found
                if( report ) write(*,*) 'mc found better:', y 
                self%yb = y                                           ! update the best cost
                if( allocated(self%probs) )then        
                    self%statecnts(state) = self%statecnts(state)+1   ! update the state counters
                    state_best = state                                ! set best state 
                endif
                if( refine .ne. 'het' .and. refine .ne. 'stoch' )then
                    self%spec_linmin%x  = x                           ! set solution in linmin spec
                    self%spec_linmin%xi = x-x_old                     ! estimate the search direction
                    call linmin(self%spec_linmin,y2)                  ! line minimize
                    if( y2 < self%yb )then                            ! better solution found
                        self%yb = y2                                  ! update the best cost
                        x = self%spec_linmin%x                        ! update solution
                        if( report )then
                            write(*,*) 'linmin found better:', y2
                        endif
                    endif
                endif
                call self%update_gaup(x)                              ! update Gaussian part of the model
                call self%update_probs                                ! update the multinomal part of the model
                ! correct if cyclic
                call check_and_correct_vec(self%spec_linmin, x, corrected)
                if( corrected )then ! shift the Euler means if x was corrected
                    do j=1,3
                        call self%vars(j)%reset_mean(x(j))
                    end do
                endif
                ! add local minimum to the list of local minima
                if( allocated(self%probs) )then
                    call self%add_minimum(x, self%yb, state=state)
                else
                    call self%add_minimum(x, self%yb)
                endif
            else
                x = x_old ! put back the old point
            endif
        end do
        ! update the parameters in ori object
        if( found_better )then
            if( allocated(self%probs) )then
                call self%update_ori_part(x, state_best)
            else
                call self%update_ori_part(x)
            endif
        endif
    end subroutine
    
    !>  \brief  is 4 updating the Gaussian part of the model
    subroutine update_gaup( self, x )
        class(pori), intent(inout) :: self
        real, intent(in)           :: x(self%ngaup)
        integer :: j
        real    :: var, sdev
        do j=1,self%ngaup                ! [1,ngaup]
            call self%vars(j)%add(x(j))  ! add local minimum to the Gaussian part of the probabilistic model
            var = self%vars(j)%get_var() ! extract variance
            if( var /= 0. )then
                sdev = sqrt(var) ! update the standard deviations
                if( j < 4 ) then
                    self%sdevs(j) = max(MINEULSDEV,sdev) ! Euler angle standard deviations
                else
                    self%sdevs(j) = max(MINTRSSDEV,sdev) ! shift standard deviations
                endif
            endif
        end do
        if( report ) write(*,*) 'updated gaup'
    end subroutine
    
    !>  \brief  is 4 updating the state probs given statecnts
    !!          A learning is required for the initial stages
    !!          because there will be large fluctuations in 
    !!          self%probs   
    subroutine update_probs( self )
        class(pori), intent(inout) :: self
        integer(long)              :: sumcounts
        integer                    :: nstates, s
        real, parameter            :: eps=0.5
        if( allocated(self%probs) )then
            nstates = size(self%probs)
            sumcounts = sum(self%statecnts)
            do s=1,nstates
                self%probs(s) = (1.-eps)*self%probs(s)+eps*real(dble(self%statecnts(s))/dble(sumcounts))
            end do
        endif
        if( report ) write(*,*) 'updated probs'
    end subroutine
    
    !>  \brief  is 4 updating the ori part of the object
    subroutine update_ori_part( self, x, state )
        use simple_math, only: rad2deg
        class(pori), intent(inout)    :: self
        real, intent(in)              :: x(self%ngaup)
        integer, intent(in), optional :: state
        real    :: dist
        integer :: state_prev
        call self%optr%set('corr', -self%yb) ! update the best cost
        call self%optr%set_euler(x(1:3))
        ! shifts must be obtained by vector addition
        call self%optr%set('x',o_old%get('x')-x(4)) ! update x-shift, revshsgn to fit shift convention
        call self%optr%set('y',o_old%get('y')-x(5)) ! update y-shift, revshsgn to fit shift convention
        if( present(state) )then
            call self%optr%set('state',real(state))
            state_prev = nint(o_old%get('state'))
            if( state_prev == state )then
                call self%optr%set('mi_hard', 1.)
            else
                call self%optr%set('mi_hard', 0.)
            endif
        endif
        if( self%ngaup > 5 )then
            call self%optr%set('dfx',x(6))
            if( self%ngaup > 6 ) call self%optr%set('dfy',x(7))
        endif
        call self%optr%set('nbetter', 1.) ! 2 indicate that we found better best (MUST BE ZEROED AT INIT)
        dist = 0.5*rad2deg(o_old.euldist.self%optr)+0.5*o_old%get('dist')
        call self%optr%set('dist',dist)
        if( report ) write(*,*) 'updated ori_part'
    end subroutine
    
    ! NLP-REFINEMENT PROCEDURE
    
    !>  \brief  for refining a local minimum further using NLP
    subroutine refine_minimum( self, opt, ospec, which )
        use simple_optimizer, only: optimizer
        use simple_ori,       only: ori
        use simple_opt_subs,  only: check_and_correct_vec
        class(pori), intent(inout)      :: self
        class(optimizer), intent(inout) :: opt 
        class(opt_spec), intent(inout)  :: ospec
        integer, intent(in), optional   :: which
        integer :: pos, j, s
        real    :: yb, y
        logical :: corrected
        if( present(which) )then
            if( which >= 1 .and. which <= self%locopts%size() )then
                pos = which
            else
                write(*,*) 'which=', which
                write(*,*) 'locopts%size=',  self%locopts%size()
                stop 'which out of range; refine_minimum; simple_pori'
            endif
        else
            pos = self%locopts%size()
        endif
        if( allocated(self%probs) )then
            call self%get_minimum(ospec%x, yb, state=s, which=pos)
            call self%ori_glob%set('state', real(s))
        else
            call self%get_minimum(ospec%x, yb, which=pos)
        endif
        ! minimize
        call opt%minimize(ospec, y)
        if( y < yb )then ! update the local minimum
            if( report ) write(*,*) 'NLP found better:', y
            ! first update the Gaussian model
            call self%update_gaup(ospec%x)
            ! then, correct if cyclic
            call check_and_correct_vec(ospec, ospec%x, corrected)
            ! replace the local minimum
            call self%locopts%set(pos, rarr=[ospec%x,yb,0.])
            if( y < self%yb )then ! if better than previous best
                if( report ) write(*,*) 'NLP found better best:', y
                ! update best cost
                self%yb = y
                if( corrected )then ! shift the Euler means if x was corrected
                    do j=1,3
                        call self%vars(j)%reset_mean(ospec%x(j))
                    end do
                endif
                call self%update_ori_part(ospec%x)
            endif
        endif
        if( report ) write(*,*) 'refined minimum'
    end subroutine
    
    ! ORIENTATION SPREAD CALCULATOR
    
    !>  \brief  produce a measure of the orientation spread
    !!          should be calculated before binary tournament
    subroutine calc_orispread( self )
        use simple_math, only: rad2deg
        class(pori), intent(inout) :: self
        type(ori) :: o1, o2
        real      :: adist
        integer   :: nbetter, i
        adist = 0.
        nbetter = self%get_nbetter()
        if( nbetter > 1 )then
            call self%minimum2ori(nbetter, o1)
            do i=nbetter-1,1,-1
                call self%minimum2ori(i, o2)
                adist = adist+rad2deg(o1.euldist.o2)
            end do
            adist = adist/real(nbetter-1)
        endif
        call self%optr%set('sdev', adist)
        if( report ) write(*,*) 'did calculate orispread'
    end subroutine
    
    ! WEIGHT CALCULATOR
    
    !>  \brief  4 calculating orientation weights
    !!          this should be done before tournament decimation of the population
    subroutine calc_weights( self )
        use simple_stat, only: normalize_net
        class(pori), intent(inout) :: self
        integer                    :: loptsz
        real, allocatable          :: weights(:)
        integer, allocatable       :: state(:)
        real                       :: wsum
        integer                    :: i
        loptsz = self%locopts%size()
        if( loptsz == 1 )then
            call set_weight(1, 1.)
        else
            call extract_corrs ! initial weights are corr vals
            call normalize_net(weights, 'sigm') ! [0,1] normalization
            ! calculate orientation weights that are  weighted
            ! by the state probabilities (by multiplying them)
            wsum = 0.
            do i=1,loptsz
                if( allocated(self%probs) )then
                    call self%locopts%get(i, iarr=state)
                    weights(i) = exp(weights(i)*self%probs(state(1)))
                else
                    weights(i) = exp(weights(i))
                endif
                wsum = wsum+weights(i)
            end do
            ! normalize the weights and update the list structure
            do i=1,loptsz
                call set_weight(i, weights(i)/wsum)
            end do
            deallocate(weights)
        endif
        if( report ) write(*,*) 'weights calculated'
        
        contains
        
            !>  \brief  for setting the orientation weight in the list structure
            subroutine set_weight( i, w )
                integer, intent(in) :: i
                real, intent(in)    :: w
                call self%locopts%get(i, rarr=globarr)
                if( size(globarr) <  self%ngaup+2)then
                    stop 'nonconformming globarr size; set_weight; calc_weights; simple_pori'
                endif 
                globarr(self%ngaup+2) = w                
                call self%locopts%set(i, rarr=globarr)
            end subroutine
            
            !>  \brief  for extracting the correlations of the local minima
            subroutine extract_corrs
                use simple_jiffys, only: alloc_err
                integer :: alloc_stat, i
                allocate( weights(loptsz), stat=alloc_stat )
                call alloc_err('In: extract_corrs; calc_weights; simple_pori', alloc_stat)
                do i=1,loptsz
                    call self%locopts%get(i, rarr=globarr)
                     weights(i) = -globarr(self%ngaup+1)
                end do
            end subroutine
            
    end subroutine
    
    ! TOURNAMENT SELECTION OF LOCAL MINIMA
    
    !>  \brief  decimate the population of local optima by random binary 
    !!          tournament (always include the best though)
    subroutine binary_tournament( self )
        use simple_rnd, only: irnd_uni_pair
        class(pori), intent(inout) :: self
        type(sll)                  :: winners
        integer                    :: rndpair(2), loptsz
        integer, allocatable       :: iarr(:)
        call winners%new ! list of winners
        loptsz = self%locopts%size() 
        if( loptsz <= IMPORTANT )then
            ! include all with reversed order
            do
                loptsz = self%locopts%size()
                if( loptsz > 0 )then
                    call move2winners(loptsz)
                else
                    exit
                endif
            end do
        else
            ! always include the best
            call move2winners(loptsz)
            ! now, random binary tournament
            do
                loptsz = self%locopts%size()
                if( loptsz > 2 )then
                    rndpair = irnd_uni_pair(loptsz)
                    call move2winners(find_winner(rndpair))
                else if( loptsz == 2 )then
                    rndpair = [1,2]
                    call move2winners(find_winner(rndpair))
                else
                    call move2winners(1)
                endif
                if( winners%size() == IMPORTANT ) exit
            end do
        endif
        ! make the locopts list a replica of winners
        self%locopts = winners
        if( report ) write(*,*) 'binary tournament finished'
        
        contains
            
            !>  \brief  for moving node i of locopts 2 winners
            subroutine move2winners( i )
                integer, intent(in) :: i
                if( allocated(self%probs) )then ! include state (iarr)
                    call self%locopts%get(i, iarr=iarr, rarr=globarr)
                    call winners%add(iarr=iarr, rarr=globarr)
                else
                    call self%locopts%get(i, rarr=globarr)
                    call winners%add(rarr=globarr)
                endif
                call self%locopts%del(i)
            end subroutine
            
            !>  \brief  for finding the winner in a pair
            function find_winner( pair ) result( winner )
                integer, intent(in) :: pair(2)
                integer :: winner
                real    :: cost1, cost2
                call self%locopts%get(pair(1), rarr=globarr)
                cost1 = globarr(self%ngaup+1)
                call self%locopts%get(pair(2), rarr=globarr)
                cost2 = globarr(self%ngaup+1) 
                if( cost1 <= cost2 )then
                    winner = pair(1)
                else
                    winner = pair(2)
                endif
            end function

    end subroutine

    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine test_alloc_stat( self )
        class(pori), intent(inout) :: self
        if( allocated(globarr) )then
            print *, 'globarr   allocated'
        else
            print *, 'globarr   deallocated'
        endif
        if( allocated(self%vars) )then
            print *, 'vars      allocated'
        else
            print *, 'vars      deallocated'
        endif
        if( allocated(self%sdevs) )then
            print *, 'sdevs     allocated'
        else
            print *, 'sdevs     deallocated'
        endif 
        if( allocated(self%probs) )then
            print *, 'probs     allocated'
        else
            print *, 'probs     deallocated'
        endif
        if( allocated(self%statecnts) )then
            print *, 'statecnts allocated'
        else
            print *, 'statecnts deallocated'
        endif
        if( associated(self%ori_glob) )then
            print *, 'ori_glob  associated'
        else
            print *, 'ori_glob  unassociated'
        endif
    end subroutine
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(pori), intent(inout) :: self
        self%optr => null()
        self%ngaup = 0
        if( allocated(globarr) )        deallocate(globarr)
        if( allocated(self%vars) )      deallocate(self%vars)
        call self%locopts%kill
        call self%spec_linmin%kill
        self%ori_glob => null()
        if( allocated(self%sdevs) )     deallocate(self%sdevs)
        if( allocated(self%probs) )     deallocate(self%probs)
        if( allocated(self%statecnts) ) deallocate(self%statecnts)
    end subroutine
    
    ! UNIT TEST
    
    !>  \brief  is the unit test for the pori class
    subroutine test_pori
    ! this test only tests the Euler part of pori,
    ! the rest is tested in the oris class
        use simple_stat,     only: pearsn
        use simple_math,     only: rad2deg, euclid
        use simple_opt_spec, only: opt_spec
        use simple_jiffys,   only: calc_pori_recsz
        type(ori)             :: ori_glob
        type(pori)            :: e1, e2
        type(opt_spec)        :: ospec
        real(dp), allocatable :: arr(:)
        logical               :: passed
        real                  :: limits(5,2), x(5), yb
        integer               :: recsz, s
        character(len=STDLEN) :: fname='poris.bin'
        limits        = 0.
        limits(1:3,2) = 360.
        limits(2,2)   = 180.
        limits(4:5,1) = -2.
        limits(4:5,2) = 2.
        call ospec%specify('linmin', 5, nstates=2, limits=limits)        
        call e1%new(ospec, ori_glob, ori_glob) 
        call e2%new(ospec, ori_glob, ori_glob)      
        e1%sdevs     = [1.,2.,3.,4.,5.]
        e1%statecnts = [10,20]
        e2%sdevs     = [1.,2.,3.,4.,5.]
        e2%statecnts = [10,20]
        if( .not. e2%exists()) stop 'e2 does not exist!'
        if( .not. e1%exists()) stop 'e1 does not exist!'
        passed = .false.
        if( euclid(e1%sdevs,e2%sdevs)                     < 1e-9 ) passed = .true.
        if( euclid(real(e1%statecnts),real(e2%statecnts)) < 1e-9 ) passed = .true.
        if( .not. passed ) stop 'pori assigner corrupt!'
        arr = e1%serialize()
        inquire(iolength=recsz) arr
        if( e1%get_recsz() /= calc_pori_recsz(5, 2) ) stop 'theoretical/calculated recsz mismatch'
        call e2%unserialize(arr, ori_glob, ori_glob)
        if( euclid(e1%sdevs,e2%sdevs)                     < 1e-9 ) passed = .true.
        if( euclid(real(e1%statecnts),real(e2%statecnts)) < 1e-9 ) passed = .true.
        if( .not. passed ) stop 'pori serialization corrupt!'
        write(*,'(a)') '**info(simple_pori_unit_test: testing local minima container'
        call e1%add_minimum([1.,0.,0.,0.,0.], -0.1, state=1)
        call e1%add_minimum([0.,1.,0.,0.,0.], -0.2, state=2)
        call e1%add_minimum([0.,0.,1.,0.,0.], -0.3, state=1)
        call e1%add_minimum([0.,0.,0.,1.,0.], -0.4, state=2)
        call e1%add_minimum([0.,0.,0.,0.,1.], -0.5, state=1)
        call e1%add_minimum([0.,0.,0.,1.,1.], -0.6, state=2)
        call e1%add_minimum([0.,0.,1.,1.,1.], -0.7, state=1)
        call e1%add_minimum([0.,1.,1.,1.,1.], -0.8, state=2)
        call e1%add_minimum([1.,1.,1.,1.,1.], -0.9, state=1)
        call e1%add_minimum([0.,0.,0.,0.,0.], -1.0, state=2)
        call e1%print_minima
        print *, '***********'
        call e1%get_minimum(x, yb, which=1, state=s)
        print *, x, yb, s
        call e1%get_minimum(x, yb, which=2, state=s)
        print *, x, yb, s
        call e1%get_minimum(x, yb, which=3, state=s)
        print *, x, yb, s
        call e1%get_minimum(x, yb, which=4, state=s)
        print *, x, yb, s
        call e1%get_minimum(x, yb, which=5, state=s)
        print *, x, yb, s
        passed = .false.
        if( e1%get_nbetter() == 10 ) passed = .true.
        if( .not. passed ) stop 'local minima container corrupt!'
        call e1%calc_weights
        call e1%print_minima
        print *, '***********BINARY TORNAMENT FOLLOWS'
        call e1%binary_tournament
        call e1%print_minima
        write(*,'(a)') '**info(simple_pori_unit_test: testing I/O'
        call e1%write(fname,1)
        call e2%read(fname, 1, ori_glob, ori_glob)        
        if( euclid(e1%sdevs,e2%sdevs)                     < 1e-9 ) passed = .true.
        if( euclid(real(e1%statecnts),real(e2%statecnts)) < 1e-9 ) passed = .true.
        if( .not. passed ) stop 'I/O corrupt!'
        call e1%test_alloc_stat
        call e1%kill
        call e1%test_alloc_stat
        call e2%kill
        if( e2%exists()) stop 'e2 does exist (despite killed)!'
        if( e1%exists()) stop 'e1 does exist (despite killed)!'
        write(*,'(a)') 'SIMPLE_PORI_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine
    
end module simple_pori
