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
use simple_defs        ! singleton
use simple_ori,        only: ori
use simple_online_var, only: online_var
use simple_sll,        only: sll
! use simple_opt_spec,   only: opt_spec
use simple_math
implicit none

public :: pori, test_pori
private

logical :: debug=.false.

type :: pori
    private
    type(ori), pointer            :: optr=>null()     !< pointer to ori obj
    integer                       :: ngaup=0          !< number of Gaussian parameters (5-7) depending on defocus refinement or not
    type(online_var), allocatable :: vars(:)          !< objects for online mean and variance estimation from all local minima
                                                      !! encountered throughout the history of the search
    type(sll)                     :: locopts          !< list of local optimas identified during this round of search                         
    real, allocatable             :: sdevs(:)         !< the stdevs of the prior Gaussian distributions
                                                      !! the means of the priors are not needed, they are in vars
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
    procedure          :: optr2arr
    procedure          :: add_minimum                  
    procedure          :: get_minimum
    procedure          :: minimum2ori                 
    procedure          :: print_minima                 
    ! I/O
    procedure          :: get_recsz                    
    procedure          :: write                   
    procedure          :: read                  
    ! MODEL UPDATE               
    procedure          :: update_gaup                          
    procedure          :: update_probs          
    ! ORISPREAD CALCULATOR
    procedure          :: calc_orispread
    ! WEIGHT CALCULATOR
    procedure          :: calc_weights                            
    ! TOURNAMENT SELECTION 
    procedure          :: binary_tournament
    ! DESTRUCTOR
    procedure          :: kill                  
end type

real, allocatable :: globarr(:)

! SELF%LOCOPTS CONTAINS:
! iarr(1)=state
! rarr(1)=e1
! rarr(2)=e2
! rarr(3)=e3
! rarr(4)=x
! rarr(5)=y
! rarr(self%ngaup+1)=cost
! rarr(self%ngaup+2)=weight

contains

    ! CONSTRUCTORS
    
    !>  \brief  is a constructor
    subroutine new( self, o, ndim, nstates )
        class(pori), intent(inout)    :: self
        type(ori), intent(in), target :: o
        integer, intent(in)           :: ndim, nstates
        integer :: alloc_stat
        call self%kill
        ! set pointer to ori
        self%optr => o
        ! set the number of Gaussian parameters
        self%ngaup = ndim
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
        if( nstates > 1 )then ! allocate stateprobs and counters
            allocate(self%probs(nstates), self%statecnts(nstates), stat=alloc_stat)
            call alloc_err('In: new_pori; simple_pori, 2', alloc_stat)
            self%probs     = 1./real(nstates)
            self%statecnts = 0
        endif
        if( debug ) write(*,*) 'new_pori created'
    end subroutine
    
    !>  \brief  for serialization of a pori object
    function serialize( self ) result( arr )
        ! ngaup             => 1                    => 1
        ! nstates           => 1                    => 2
        ! vars              => 4*self%ngaup         => 3->2+4*self%ngaup
        ! sdevs             => self%ngaup           => 2+4*self%ngaup+1->2+5*self%ngaup
        ! if( nstates > 1 ) => size(self%probs)     => 2+5*self%ngaup+1->2+5*self%ngaup+nstates
        ! if( nstates > 1 ) => size(self%statecnts) => 2+5*self%ngaup+nstates+1->2+5*self%ngaup+2*nstates
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
        n = 2+5*self%ngaup+2*nstates
        if( n < 27 )then
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
        ! vars:
        cnt = 0
        if( debug ) write(*,*) 'simple_pori::serialize, online_vars range: ', 3, 2+4*self%ngaup
        do i=3,2+4*self%ngaup,4
            cnt = cnt+1
            arr(i:i+3) = self%vars(cnt)%serialize()
        end do
        ! sdevs:
        if( debug ) write(*,*) 'simple_pori::serialize, sdevs range: ', 2+4*self%ngaup+1, 2+5*self%ngaup
        arr(2+4*self%ngaup+1:2+5*self%ngaup) = self%sdevs
        if( nstates > 1 )then
            ! probs:
            if( debug ) write(*,*) 'simple_pori::serialize, probs range: ', 2+5*self%ngaup+1, 2+5*self%ngaup+nstates
            arr(2+5*self%ngaup+1:2+5*self%ngaup+nstates) = dble(self%probs)
            ! statecnts:
            if( debug ) write(*,*) 'simple_pori::serialize, statecnts range: ', 2+5*self%ngaup+nstates+1, 2+5*self%ngaup+2*nstates
            arr(2+5*self%ngaup+nstates+1:2+5*self%ngaup+2*nstates) = dble(self%statecnts)
        endif
        if( debug ) write(*,*) 'pori serialized'
    end function
    
    !>  \brief  for re-creation of a pori object from its serialization
    subroutine unserialize( self, arr, o )
        use simple_opt_spec, only: opt_spec
        ! ngaup             => 1                    => 1
        ! nstates           => 1                    => 2
        ! vars              => 4*self%ngaup         => 3->2+4*self%ngaup
        ! sdevs             => self%ngaup           => 2+4*self%ngaup+1->2+5*self%ngaup
        ! if( nstates > 1 ) => size(self%probs)     => 2+5*self%ngaup+1->2+5*self%ngaup+nstates
        ! if( nstates > 1 ) => size(self%statecnts) => 2+5*self%ngaup+nstates+1->2+5*self%ngaup+2*nstates
        class(pori), intent(inout)    :: self
        real(dp), intent(inout)       :: arr(:)
        type(ori), intent(in), target :: o
        integer :: n, ngaup, i, cnt, nstates, alloc_stat
        ! ntot:
        n = 2+5*self%ngaup+2*size(self%probs)
        call self%kill
        ngaup   = nint(arr(1))
        nstates = nint(arr(2))
        if( debug ) write(*,*) 'simple_pori::unserialize, ntot: ', n
        if( size(arr) /= n )then
            write(*,*) 'nonconforming arr size simple_pori::unserialize, expected: ', n, 'got: ', size(arr)
            stop
        endif
        ! create object:
        call self%new(o, ngaup, nstates)
        ! online_vars:
        cnt = 0
        if( debug ) write(*,*) 'simple_pori::unserialize, online_vars range: ', 3, 2+4*self%ngaup
        do i=3,2+4*self%ngaup,4
            cnt = cnt+1
            call self%vars(cnt)%unserialize(arr(i:i+3))
        end do
        ! sdevs:
        if( debug ) write(*,*) 'simple_pori::unserialize, sdevs range: ', 2+4*self%ngaup+1, 2+5*self%ngaup
        self%sdevs = real(arr(2+4*self%ngaup+1:2+5*self%ngaup))
        if( nstates > 1 )then
            ! probs:
            if( debug ) write(*,*) 'simple_pori::unserialize, probs range: ', 2+5*self%ngaup+1, 2+5*self%ngaup+nstates
            self%probs = real(arr(2+5*self%ngaup+1:2+5*self%ngaup+nstates))
            ! statecnts:
            if( debug ) write(*,*) 'simple_pori::unserialize, statecnts range: ', 2+5*self%ngaup+nstates+1, 2+5*self%ngaup+2*nstates
            self%statecnts = nint(arr(2+5*self%ngaup+nstates+1:2+5*self%ngaup+2*nstates))
        endif
         if( debug ) write(*,*) 'pori unserialized'
    end subroutine
    
    ! GETTERS/SETTERS/ASSIGNERS
    
    !>  \brief  2 check existence
    function exists( self ) result( t )
        class(pori), intent(in) :: self
        logical :: t
        t = allocated(self%vars)
    end function
    
    !>  \brief  returns the number of local optimas identified
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
    function optr2arr( self ) result( x )
        class(pori), intent(inout) :: self
        real, allocatable          :: x(:)
        integer                    :: i, alloc_stat
        allocate( x(self%ngaup), stat=alloc_stat )
        call alloc_err("In: optr2arr; simple_pori", alloc_stat)
        ! extract the previous solution
        do i=1,self%ngaup
            if( i == 1 )then
                x(i) = self%optr%e1get()
            else if( i == 2 )then
                x(i) = self%optr%e2get()
            else if( i == 3 )then
                x(i) = self%optr%e3get()
            else if( i == 4 )then
                x(i) = -self%optr%get('x') 
            else if( i == 5 )then
                x(i) = -self%optr%get('y')
            endif
        end do
        if( debug ) write(*,*) 'previous ori params extracted'
    end function
    
    !>  \brief  for adding a new local minimum to the list
    subroutine add_minimum( self, x, yb, w, state )
        class(pori), intent(inout)    :: self
        real, intent(in)              :: x(self%ngaup), yb
        real, intent(in), optional    :: w
        integer, intent(in), optional :: state
        integer :: i, sz
        if( debug )then
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
        if( debug ) write(*,*) 'local minimum added'
    end subroutine
    
    !>  \brief  for printing the local minima
    subroutine print_minima( self )
        class(pori), intent(in) :: self
        call self%locopts%print()
    end subroutine
    
    !>  \brief  for getting a local minimum from the list
    subroutine get_minimum( self, which, x, yb, w, state )
        class(pori), intent(in)        :: self
        integer, intent(in)            :: which
        real, intent(out)              :: x(self%ngaup), yb
        real, intent(out), optional    :: w
        integer, intent(out), optional :: state
        integer, allocatable           :: iarr(:)
        integer :: pos, sz
        if( which >= 1 .and. which <= self%locopts%size() )then
            pos = which
        else
            write(*,*) 'which=', which
            write(*,*) 'locopts%size=',  self%locopts%size()
            stop 'which out of range; get_minimum; simple_pori'
        endif
        if( present(state) )then
            call self%locopts%get(pos, rarr=globarr, iarr=iarr)
            state = iarr(1)
        else
            call self%locopts%get(pos, rarr=globarr)
        endif
        if( debug )then
            sz = size(globarr)
            if( sz < self%ngaup+2 )then
                write(*,*) 'nonconforming array size of sll extracted globarr:', sz
                stop
            endif 
        endif
        x  = globarr(1:self%ngaup)
        yb = globarr(self%ngaup+1)
        if( present(w) ) w = globarr(self%ngaup+2)
        if( debug ) write(*,*) 'local minimum gotten'
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
            call self%get_minimum(which, x, yb, w, state)
        else
            call self%get_minimum(which, x, yb, w)
        endif
        if( .not. o%exists() ) call o%new ! create ori 4 output 
        call o%set('corr', -yb)  ! update the correlation
        call o%set('ow', w)      ! update the orientation weight
        call o%set_euler(x(1:3)) ! update the Euler angle
        ! shifts must be obtained by vector addition
        call o%set('x',x(4))            ! update x-shift
        call o%set('y',x(5))            ! update y-shift
        call o%set('state',real(state)) ! update state
    end subroutine
    
    ! I/O
    
    !>  \brief  for writing a pori object to file
    subroutine write( self, fname, i )
        use simple_filehandling, only: fopen,fclose, err_msg
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
        inquire(file=fname, exist=here)
        ostat = 'new'
        if( here ) ostat = 'old'
        if(.not.fopen(fnr, status=ostat, action='write', file=fname,&
             access='direct', form='unformatted', recl=recsz, iostat=file_stat))&
             call err_msg( 'In: write; simple_pori', file_stat )
        write(fnr,rec=i) arr
        if(.not.fclose(fnr,file_stat))&
             call err_msg( 'In: write; simple_pori', file_stat )
    end subroutine

    !>  \brief  for re-creating a pori object from file 
    subroutine read( self, fname, i, o )
        use simple_jiffys, only: get_fileunit, fopen_err, alloc_err
        class(pori), intent(inout)    :: self
        character(len=*), intent(in)  :: fname
        integer, intent(in)           :: i
        type(ori), intent(in), target :: o
        real(dp), allocatable         :: arr(:)
        integer                       :: recsz, fnr, file_stat, elemsz, alloc_stat
        logical                       :: here
        real(dp)                      :: tmp
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
            call self%unserialize(arr, o)
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
        if( debug ) write(*,*) 'updated gaup'
    end subroutine
    
    !>  \brief  is 4 updating the state probs given statecnts
    !!          A learning rate is required for the initial stages
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
        if( debug ) write(*,*) 'updated probs'
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
        if( debug ) write(*,*) 'did calculate orispread'
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
        if( debug ) write(*,*) 'weights calculated'
        
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
                use simple_fileio, only: alloc_err
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
        if( debug ) write(*,*) 'binary tournament finished'
        
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
    subroutine kill( self )
        class(pori), intent(inout) :: self
        self%optr => null()
        self%ngaup = 0
        if( allocated(globarr) )   deallocate(globarr)
        if( allocated(self%vars) ) deallocate(self%vars)
        call self%locopts%kill
        if( allocated(self%sdevs) )     deallocate(self%sdevs)
        if( allocated(self%probs) )     deallocate(self%probs)
        if( allocated(self%statecnts) ) deallocate(self%statecnts)
    end subroutine
    
    ! UNIT TEST
    
    !>  \brief  is the unit test for the pori class
    subroutine test_pori
        use simple_stat,     only: pearsn
        use simple_math,     only: rad2deg, euclid
        use simple_jiffys,   only: calc_pori_recsz
        type(ori)             :: ori_glob
        type(pori)            :: e1, e2
        real(dp), allocatable :: arr(:)
        logical               :: passed
        real                  :: x(5), yb
        integer               :: recsz, s
        character(len=STDLEN) :: fname='poris.bin'
        call e1%new(ori_glob, 5, 2) 
        call e2%new(ori_glob, 5, 2)      
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
        call e2%unserialize(arr, ori_glob)
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
        call e1%get_minimum(1, x, yb, state=s)
        print *, x, yb, s
        call e1%get_minimum(2, x, yb, state=s)
        print *, x, yb, s
        call e1%get_minimum(3, x, yb, state=s)
        print *, x, yb, s
        call e1%get_minimum(4, x, yb, state=s)
        print *, x, yb, s
        call e1%get_minimum(5, x, yb, state=s)
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
        call e2%read(fname, 1, ori_glob)        
        if( euclid(e1%sdevs,e2%sdevs)                     < 1e-9 ) passed = .true.
        if( euclid(real(e1%statecnts),real(e2%statecnts)) < 1e-9 ) passed = .true.
        if( .not. passed ) stop 'I/O corrupt!'
        call e1%kill
        call e2%kill
        if( e2%exists()) stop 'e2 does exist (despite killed)!'
        if( e1%exists()) stop 'e1 does exist (despite killed)!'
        write(*,'(a)') 'SIMPLE_PORI_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine
    
end module simple_pori
