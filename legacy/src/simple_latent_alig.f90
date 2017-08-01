!>  \brief  SIMPLE class that defines a continuous latent variable model over refinement (alignment) parameters
module simple_latent_alig
use simple_defs ! singleton
implicit none

public :: latent_alig
private

type :: latent_alig
    private
    integer              :: nsample        !< number of samples per iteration during search 
    integer              :: ngaup          !< the number of Gaussian parameters
    real, allocatable    :: gauparams(:,:) !< the parameter vectors
    real, allocatable    :: means(:)       !< the means of the prior Gaussian distributions
    real, allocatable    :: best_sol(:)    !< the best solution found so far
    real, allocatable    :: sdevs(:)       !< the stdevs of the prior Gaussian distributions
    real, allocatable    :: probs(:)       !< the multinomal distribution of probabilities for 1ofK state encoding
    real, allocatable    :: probs_new(:)   !< for updating probs
    real, allocatable    :: weights(:)     !< the orientation weights
    real, allocatable    :: scores(:)      !< the scores of the individual parameter vec:s
    real, allocatable    :: limits(:,:)    !< variable limits
    integer, allocatable :: states(:)      !< the discrete state variables
    integer, allocatable :: order(:)       !< the order of solution according to score
    logical, allocatable :: cyclic(:)      !< to indicate which variables are cyclic
    real                 :: glob_wscore=0. !< weighted global score
    real                 :: best_score=-1. !< best score
    real                 :: prev_score=-1. !< previous score
    real                 :: eps=0.5        !< learning rate for the state probabilities
    integer              :: nbetter=0      !< the number of better solutions found
    integer              :: nstates=1      !< the number of states
    integer              :: important=6    !< number of solutions considered important
    logical              :: exists=.false. !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! SETTERS
    procedure          :: zero_all
    procedure          :: zero_nbetter
    procedure          :: zero_probs
    procedure, private :: set_1
    procedure, private :: set_2
    procedure, private :: set_3
    generic            :: set => set_1, set_2, set_3
    procedure          :: update_means
    procedure          :: replace
    ! GETTERS
    procedure          :: get_nbetter
    procedure, private :: get_solution_1
    procedure, private :: get_solution_2
    generic            :: get_solution => get_solution_1, get_solution_2
    procedure          :: get_stateprobs
    ! MODEL UPDATE & SAMPLING
    procedure          :: sort
    procedure          :: update
    procedure, private :: sample_1
    procedure, private :: sample_2
    generic            :: sample => sample_1, sample_2
    procedure, private :: sample_uni_1
    procedure, private :: sample_uni_2
    generic            :: sample_uni => sample_uni_1, sample_uni_2
    ! DESTRUCTOR
    procedure          :: kill
end type

contains

    ! CONSTRUCTOR
    
    !>  \brief  is a constructor
    subroutine new( self, p )
        use simple_params, only: params
        use simple_jiffys, only: alloc_err
        class(latent_alig), intent(inout) :: self
        class(params), intent(inout)      :: p
        integer                           :: alloc_stat
        call self%kill
        self%ngaup     = p%ndim
        self%nsample   = p%nsample
        self%nstates   = p%nstates
        self%important = IMPORTANT
        if( p%refine .eq. 'no' ) self%important = p%npeaks
        allocate( self%gauparams(self%nsample,self%ngaup), self%means(self%ngaup),&
        self%best_sol(self%ngaup), self%sdevs(self%ngaup), self%probs(self%nstates),&
        self%probs_new(self%nstates), self%scores(self%nsample), self%states(self%nsample),&
        self%weights(self%nsample), self%order(self%nsample), self%cyclic(self%ngaup),&
        self%limits(self%ngaup,2), stat=alloc_stat )
        call alloc_err("In: new_1; simple_latent_alig", alloc_stat)
        self%gauparams = 0.
        self%means     = 0.
        self%best_sol  = 0.
        self%sdevs     = 0.
        self%probs     = 1./real(self%nstates)
        self%probs_new = self%probs
        self%scores    = 0.
        self%states    = 1
        self%weights   = 0.
        self%order     = 0
        self%cyclic    = p%cyclic
        self%limits    = p%optlims(:p%ndim,:)
        self%exists    = .true.
    end subroutine
    
    ! SETTERS
    
    !>  \brief  zeroes all counters
    subroutine zero_all( self )
        class(latent_alig), intent(inout) :: self
        integer :: i
        call self%zero_nbetter
        call self%zero_probs
        self%order(:self%nsample)=(/(i,i=1,self%nsample)/)
        self%scores = -1.
    end subroutine
    
    !>  \brief  zeroes the nbetter counter variable
    subroutine zero_nbetter( self )
        class(latent_alig), intent(inout) :: self
        self%nbetter = 0
        self%best_score = -1.
    end subroutine
    
    !>  \brief  zeroes the state probabilities
    subroutine zero_probs( self )
        class(latent_alig), intent(inout) :: self
        self%probs_new = 0.
    end subroutine
    
    !>  \brief  is a setter
    subroutine set_1( self, means, sdevs, probs )
        class(latent_alig), intent(inout) :: self
        real, intent(in)                  :: means(self%ngaup), sdevs(self%ngaup)
        real, intent(in), optional        :: probs(self%nstates)
        self%means = means
        self%sdevs = sdevs
        if( present(probs) )self%probs = probs
    end subroutine

    !>  \brief  is a setter
    subroutine set_2( self, gauparams, score, state, full )
        class(latent_alig), intent(inout) :: self
        real, intent(in)                  :: gauparams(self%ngaup), score
        integer, intent(in), optional     :: state
        logical, intent(out), optional    :: full
        logical                           :: update_probs_new
        if( present(full) ) full =.false.
        self%nbetter = self%nbetter+1
        if( self%nbetter > self%nsample )then
            self%nbetter = self%nsample
            if( present(full) )then
                full =.true.
            else
                write(*,'(a)') 'WARNING, latent_alig object is full!'
            endif
            return
        endif
        ! set orientation params and score
        self%gauparams(self%nbetter,:) = gauparams
        self%scores(self%nbetter) = score
        ! only update probs_new if the score is improving
        update_probs_new = .false.
        if( self%best_score < 0. )then
            self%prev_score = score
            update_probs_new = .true.
        else
            if( score >= self%prev_score ) update_probs_new = .true.
        endif
        ! update the states and probs_new arrays
        if( present(state) )then
            self%states(self%nbetter) = state
            if( update_probs_new ) self%probs_new(state) = self%probs_new(state)+1.
        else
            self%states(self%nbetter) = 1
        endif
        ! update the best solution
        if( score > self%best_score )then
            self%best_score = score
            self%best_sol = gauparams
        endif
    end subroutine
    
    !>  \brief  is a setter
    subroutine set_3( self, o, full )
        use simple_ori, only: ori
        class(latent_alig), intent(inout) :: self
        class(ori), intent(inout)         :: o
        logical, intent(out), optional    :: full
        real                              :: gauparams(self%ngaup), corr
        integer                           :: state
        gauparams(1:3) = o%get_euler()
        gauparams(4)   = o%get('x')
        gauparams(5)   = o%get('y')
        if( self%ngaup > 5 )then
            gauparams(6) = o%get('dfx')
            if( self%ngaup == 7 ) gauparams(7) = o%get('dfy')
        endif
        state = 1
        if( self%nstates > 1 ) state = nint(o%get('state'))
        corr = o%get('corr')
        call self%set_2(gauparams, corr, state, full)
    end subroutine
    
    !>  \brief  is 4 updating the means of the latent variable model
    subroutine update_means( self, o )
        use simple_ori, only: ori
        class(latent_alig), intent(inout) :: self
        class(ori), intent(inout)         :: o
        real :: euls(3), x, y
        euls = o%get_euler()
        x    = o%get('x')
        y    = o%get('y')
        self%means(1:3) = euls
        self%means(4)   = x
        self%means(5)   = y
        if( self%ngaup > 5 )then
            self%means(6) = o%get('dfx')
            if( self%ngaup == 7 ) self%means(7) = o%get('dfy')
        endif
    end subroutine
    
    !>  \brief  is a setter
    subroutine replace( self, ind, o )
        use simple_ori, only: ori
        class(latent_alig), intent(inout) :: self
        integer, intent(in)               :: ind
        class(ori), intent(inout)         :: o
        real                              :: gauparams(self%ngaup), corr
        logical                           :: full
        if( ind > self%nbetter )then
            call self%set(o, full)
            if( full ) stop 'Weird! Object full in replace; simple_latent_alig'
        endif
        gauparams(1:3) = o%get_euler()
        gauparams(4)   = o%get('x')
        gauparams(5)   = o%get('y')
        if( self%ngaup > 5 )then
            gauparams(6) = o%get('dfx')
            if( self%ngaup == 7 ) gauparams(7) = o%get('dfy')
        endif
        self%gauparams(ind,:) = gauparams
        self%scores(ind)      = o%get('corr')
        if( self%scores(ind) > self%best_score )then
            self%best_score   = self%scores(ind)
            self%best_sol     = gauparams
        endif
    end subroutine
    
    ! GETTERS
    
    !>  \brief  is a getter 
    function get_nbetter( self ) result( nbetter )
        class(latent_alig), intent(in) :: self
        integer :: nbetter
        nbetter = min(self%important,self%nbetter)
    end function

    !>  \brief  is a getter 
    subroutine get_solution_1( self, i, gauparams, state, score, w )
        class(latent_alig), intent(in) :: self
        integer, intent(in)  :: i
        real, intent(out)    :: gauparams(self%ngaup)
        integer, intent(out) :: state
        real, intent(out)    :: score, w
        integer              :: j
        if( self%order(i) == 0 )then
            write(*,'(a)') 'solutions in latent alignment model has not been ordered'
            write(*,'(a)') 'or the index is out of range:', i, 'nbetter:', self%nbetter
            stop
        endif
        gauparams = self%gauparams(self%order(i),:)
        state = 1
        if( self%nstates > 1 ) state = self%states(self%order(i)) 
        score = self%scores(self%order(i))
        w = self%weights(self%order(i))
        do j=1,self%ngaup
            if( isnan(gauparams(j)) ) write(*,*) 'parameter:', j, 'is NaN; get_solution; latent_alig'
        end do
        if( isnan(score) ) write(*,*) 'score is NaN; get_solution; latent_alig'
        if( isnan(w) )     write(*,*) 'w is NaN; get_solution; latent_alig'
    end subroutine
    
    !>  \brief  is a getter 
    subroutine get_solution_2( self, i, o, score, w )
        use simple_ori, only: ori
        class(latent_alig), intent(in) :: self
        integer, intent(in)            :: i
        class(ori), intent(inout)      :: o
        real, intent(out)              :: score, w
        real                           :: gauparams(self%ngaup)
        integer                        :: state
        call self%get_solution_1(i, gauparams, state, score, w)
        call o%set_euler(gauparams(1:3))
        call o%set('x', gauparams(4))
        call o%set('y', gauparams(5))
        call o%set('state', real(state))
        if( self%ngaup > 5 )then
            call o%set('dfx', gauparams(6))
            if( self%ngaup == 7 ) call o%set('dfy', gauparams(7))
        endif
    end subroutine
    
    !>  \brief  is for getting the multimodal state distribution
    subroutine get_stateprobs( self, o )
        use simple_ori, only: ori
        class(latent_alig), intent(in) :: self
        class(ori), intent(inout)      :: o
        character(len=STDLEN)          :: dig, key
        integer                        :: s, loc(1)
        if( self%nstates > 1 )then
            do s=1,self%nstates
                write(dig,*) s
                key = 'sprob'//trim(adjustl(dig))
                call o%set(trim(adjustl(key)), self%probs(s))
            end do
            loc = maxloc(self%probs)
            call o%set('state', real(loc(1)))
        endif
    end subroutine
        
    ! MODEL UPDATE & SAMPLING
    
    !>  \brief  order the solutions
    subroutine sort( self )
        use simple_math, only: hpsort, reverse
        class(latent_alig), intent(inout) :: self
        real :: scores_copy(self%nbetter)
        integer :: i
        scores_copy = self%scores(:self%nbetter)
        self%order(:self%nbetter)=(/(i,i=1,self%nbetter)/)
        call hpsort(self%nbetter,scores_copy,self%order(:self%nbetter))
        call reverse(self%order(:self%nbetter)) ! best first
    end subroutine
    
    !>  \brief  updates the statistical model
    subroutine update( self )
        use simple_stat, only: normalize_net
        class(latent_alig), intent(inout) :: self
        real    :: wsum, scores_important(self%nbetter)
        integer :: i, en
        if( self%nbetter > 1 )then        
            ! update the model
            self%means = self%best_sol
            ! calculate the state probabilities
            if( self%nstates > 1 )then
                wsum = sum(self%probs_new)
                self%probs_new = self%probs_new/wsum
                self%probs = (1.-self%eps)*self%probs+self%eps*self%probs_new
            else
                self%probs = 1.
            endif
            ! order the solutions
            call self%sort
            ! fish out the important ones & normalize
            en = min(self%important,self%nbetter)
            do i=1,en
                scores_important(i) = self%scores(self%order(i))
            end do
            call normalize_net(scores_important(:en), 'sigm')
            ! calculate the orientation weights, weighted by the state probabilities (by multiplying them)
            wsum = 0.
            self%weights = 0.
            do i=1,en
                self%weights(self%order(i)) = exp(scores_important(i)*self%probs(self%states(self%order(i))))
                wsum = wsum+self%weights(self%order(i))
            end do
            ! normalize the weights and calculate the weighted global score
            self%glob_wscore = 0.
            do i=1,en
                self%weights(self%order(i)) = self%weights(self%order(i))/wsum
                self%glob_wscore = self%glob_wscore+scores_important(i)*self%weights(self%order(i))
            end do 
        else if( self%nbetter == 1 )then
            self%weights(1) = 1.
            self%glob_wscore = self%scores(1)
            self%order(1) = 1
        endif
    end subroutine
    
    !>  \brief  samples the statistical model
    subroutine sample_1( self, spec, gauparams, state )
        use simple_opt_subs, only: check_and_correct_vec
        use simple_opt_spec, only: opt_spec
        use simple_rnd,      only: gasdev, multinomal
        class(latent_alig), intent(inout) :: self
        class(opt_spec), intent(in)       :: spec
        real, intent(out)                 :: gauparams(self%ngaup)
        integer, intent(out), optional    :: state
        integer                           :: i
        do i=1,self%ngaup
            gauparams(i) = gasdev(self%means(i), self%sdevs(i))
        end do
        call check_and_correct_vec(spec, gauparams)
        if( present(state) )then
            state = 1
            if( self%nstates > 1 ) state = multinomal(self%probs)
        endif
    end subroutine
    
    !>  \brief  samples the statistical model
    subroutine sample_2( self, spec, o )
        use simple_opt_spec, only: opt_spec
        use simple_ori,      only: ori
        class(latent_alig), intent(inout) :: self
        class(opt_spec), intent(in)       :: spec
        class(ori), intent(inout)         :: o
        integer                           :: i, state
        real                              :: gauparams(self%ngaup)
        call self%sample_1(spec,gauparams,state)
        call o%set_euler(gauparams(1:3))
        call o%set('x',gauparams(4))
        call o%set('y',gauparams(5))
        call o%set('state',real(state))
        if( self%ngaup > 5 )then
            call o%set('dfx', gauparams(6))
            if( self%ngaup == 7 ) call o%set('dfy', gauparams(7))
        endif
    end subroutine
    
    !>  \brief  samples the statistical model
    subroutine sample_uni_1( self, gauparams, state )
        use simple_rnd,      only: ran3, gasdev, multinomal
        class(latent_alig), intent(inout) :: self
        real, intent(out)                 :: gauparams(self%ngaup)
        integer, intent(out), optional    :: state
        integer                           :: i
        gauparams(1) = 360.*ran3()
        gauparams(2) = 180.*ran3()
        gauparams(3) = 360.*ran3()
        if( self%ngaup > 3 )then
            do i=4,self%ngaup
                gauparams(i) = gasdev(self%means(i), self%sdevs(i))
            end do
        endif
        if( present(state) )then
            state = 1
            if( self%nstates > 1 ) state = multinomal(self%probs)
        endif
    end subroutine
    
    !>  \brief  samples the statistical model
    subroutine sample_uni_2( self, o )
        use simple_ori, only: ori
        class(latent_alig), intent(inout) :: self
        class(ori), intent(inout)         :: o
        integer                           :: i, state
        real                              :: gauparams(self%ngaup)
        call self%sample_uni_1(gauparams,state)
        call o%set_euler(gauparams(1:3))
        call o%set('x',gauparams(4))
        call o%set('y',gauparams(5))
        call o%set('state',real(state))
        if( self%ngaup > 5 )then
            call o%set('dfx', gauparams(6))
            if( self%ngaup == 7 ) call o%set('dfy', gauparams(7))
        endif
    end subroutine
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(latent_alig), intent(inout) :: self
        if( self%exists )then
            deallocate( self%gauparams, self%means, self%best_sol,&
            self%sdevs, self%probs, self%probs_new, self%scores,&
            self%states, self%weights, self%order, self%cyclic,&
            self%limits )
            self%exists = .false. 
        endif
    end subroutine
    
end module simple_latent_alig
