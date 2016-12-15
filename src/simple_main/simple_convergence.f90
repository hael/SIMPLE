module simple_convergence
use simple_oris,    only: oris
use simple_params,  only: params
use simple_cmdline, only: cmdline
implicit none

! origin shift related parameters
real, parameter :: MSK_FRAC     = 0.06
real, parameter :: MINSHIFT     = 2.0
real, parameter :: MAXSHIFT     = 6.0
real, parameter :: FRAC_SH_LIM  = 90.0
! convergence related parameters
real, parameter :: MI_CLASS_LIM_3D = 0.80
real, parameter :: MI_CLASS_LIM_2D = 0.98
real, parameter :: MI_STATE_LIM    = 0.98
real, parameter :: MI_INPL_LIM     = 0.98
real, parameter :: FRAC_LIM        = 98.0

type convergence
    private
    class(oris),    pointer :: bap    => null() !< pointer to alignment oris object (a) part of build (b)
    class(params),  pointer :: pp     => null() !< pointer parameters singleton object
    class(cmdline), pointer :: pcline => null() !< pointer to command line object
    real :: corr      = 0.                      !< average correlation
    real :: dist      = 0.                      !< average angular distance
    real :: dist_inpl = 0.                      !< average in-plane angular distance
    real :: frac      = 0.                      !< average fraction of search space scanned
    real :: mi_joint  = 0.                      !< joint parameter distribution overlap
    real :: mi_class  = 0.                      !< class parameter distribution overlap
    real :: mi_inpl   = 0.                      !< in-plane parameter distribution overlap 
    real :: mi_state  = 0.                      !< state parameter distribution overlap
    real :: sdev      = 0.                      !< angular standard deviation of model
  contains
    procedure :: check_conv2D
    procedure :: check_conv3D
    procedure :: get
end type convergence

interface convergence
    module procedure constructor
end interface convergence

contains

    !>  \brief  is a constructor
    function constructor( ba, p, cline ) result( self )
        class(oris),    target, intent(in) :: ba    !< alignment oris object (a) part of build (b)
        class(params),  target, intent(in) :: p     !< parameters singleton object
        class(cmdline), target, intent(in) :: cline !< parameters singleton object
        type(convergence) :: self
        self%bap    => ba
        self%pp     => p
        self%pcline => cline 
    end function constructor
    
    !>  \brief  checks convergence for the 2D case
    function check_conv2D( self, ncls ) result( converged )
        class(convergence), intent(inout) :: self
        integer, optional,  intent(in)    :: ncls
        logical :: converged
        integer :: nncls
        if( present(ncls) )then
            nncls = ncls
        else
            nncls = self%pp%ncls
        endif
        self%corr      = self%bap%get_avg('corr')
        self%dist_inpl = self%bap%get_avg('dist_inpl')
        self%frac      = self%bap%get_avg('frac')
        self%mi_joint  = self%bap%get_avg('mi_joint')
        self%mi_class  = self%bap%get_avg('mi_class')
        self%mi_inpl   = self%bap%get_avg('mi_inpl')
        write(*,'(A,1X,F7.4)') '>>> JOINT    DISTRIBUTION OVERLAP:     ', self%mi_joint
        write(*,'(A,1X,F7.4)') '>>> CLASS    DISTRIBUTION OVERLAP:     ', self%mi_class
        write(*,'(A,1X,F7.4)') '>>> IN-PLANE DISTRIBUTION OVERLAP:     ', self%mi_inpl
        write(*,'(A,1X,F7.1)') '>>> AVERAGE IN-PLANE ANGULAR DISTANCE: ', self%dist_inpl
        write(*,'(A,1X,F7.1)') '>>> PERCENTAGE OF SEARCH SPACE SCANNED:', self%frac
        write(*,'(A,1X,F7.4)') '>>> CORRELATION:                       ', self%corr
        ! dynamic shift search range update
        if( self%frac >= FRAC_SH_LIM )then
            if( .not. self%pcline%defined('trs') .or. self%pp%trs <  MINSHIFT )then
                ! determine shift bounds
                self%pp%trs = MSK_FRAC*real(self%pp%msk)
                self%pp%trs = max(MINSHIFT,self%pp%trs)
                self%pp%trs = min(MAXSHIFT,self%pp%trs)
                ! set shift search flag
                self%pp%doshift = .true.
            endif
        endif
        ! determine convergence
        if( nncls > 1 )then        
            if( self%mi_class > MI_CLASS_LIM_2D .and. self%frac > FRAC_LIM )then
                write(*,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(*,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
        else
            if( self%mi_inpl > MI_INPL_LIM .or. self%dist_inpl < 0.5 )then
                write(*,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(*,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
        endif
    end function check_conv2D
    
    !>  \brief  checks convergence for the 3D case
    function check_conv3D( self, update_res ) result( converged )
        use simple_math, only: rad2deg
        class(convergence), intent(inout) :: self
        logical, optional,  intent(inout) :: update_res
        real, allocatable :: state_mi_joint(:), statepops(:)
        real              :: min_state_mi_joint
        logical           :: converged
        integer           :: iptcl, istate
        self%corr      = self%bap%get_avg('corr')
        self%dist      = self%bap%get_avg('dist')
        self%dist_inpl = self%bap%get_avg('dist_inpl')
        self%frac      = self%bap%get_avg('frac')
        self%mi_joint  = self%bap%get_avg('mi_joint')
        self%mi_class  = self%bap%get_avg('mi_class')
        self%mi_inpl   = self%bap%get_avg('mi_inpl')
        self%mi_state  = self%bap%get_avg('mi_state')
        self%sdev      = self%bap%get_avg('sdev')
        if( self%pp%athres==0. )then
            ! required for distributed mode
            self%pp%athres = rad2deg(atan(max(self%pp%fny,self%pp%lp)/(self%pp%moldiam/2.)))
        endif
        write(*,'(A,1X,F7.1)') '>>> ANGLE OF FEASIBLE REGION:          ', self%pp%athres
        write(*,'(A,1X,F7.4)') '>>> JOINT    DISTRIBUTION OVERLAP:     ', self%mi_joint
        write(*,'(A,1X,F7.4)') '>>> CLASS    DISTRIBUTION OVERLAP:     ', self%mi_class
        write(*,'(A,1X,F7.4)') '>>> IN-PLANE DISTRIBUTION OVERLAP:     ', self%mi_inpl
        if( self%pp%nstates > 1 )&
        write(*,'(A,1X,F7.4)') '>>> STATE DISTRIBUTION OVERLAP:        ', self%mi_state
        write(*,'(A,1X,F7.1)') '>>> AVERAGE ANGULAR DISTANCE BTW ORIS: ', self%dist
        write(*,'(A,1X,F7.1)') '>>> AVERAGE IN-PLANE ANGULAR DISTANCE: ', self%dist_inpl
        write(*,'(A,1X,F7.1)') '>>> PERCENTAGE OF SEARCH SPACE SCANNED:', self%frac
        write(*,'(A,1X,F7.4)') '>>> CORRELATION:                       ', self%corr
        write(*,'(A,1X,F7.2)') '>>> ANGULAR SDEV OF MODEL:             ', self%sdev
        ! automatic resolution stepping
        if( present(update_res) )then
            if( update_res )then
                ! the previous round updated the resolution limit, so
                ! don't update this round
                update_res = .false.
            else
                update_res = .false.
                if(       self%pp%dynlp .eq. 'yes'  .and. &
                    .not. self%pcline%defined('lp') .and. &
                          self%dist <= self%pp%athres/2. )then
                    update_res = .true.          
                endif
            endif
            if( update_res )then
                write(*,'(A)') '>>> UPDATE LOW-PASS LIMIT: .YES.'
            else
                write(*,'(A)') '>>> UPDATE LOW-PASS LIMIT: .NO.'
            endif
        else
            if( self%pcline%defined('find') .and. self%dist <= self%pp%athres/2. )then
                write(*,'(A)') '>>> UPDATE LOW-PASS LIMIT: .YES.'
            else
                write(*,'(A)') '>>> UPDATE LOW-PASS LIMIT: .NO.'
            endif
        endif
        ! dynamic shift search range update
        if( self%frac >= FRAC_SH_LIM )then
            if( .not. self%pcline%defined('trs') .or. self%pp%trs <  MINSHIFT )then
                ! determine shift bounds
                self%pp%trs = MSK_FRAC*real(self%pp%msk)
                self%pp%trs = max(MINSHIFT,self%pp%trs)
                self%pp%trs = min(MAXSHIFT,self%pp%trs)
                ! set shift search flag
                self%pp%doshift = .true.
            endif
        endif
        ! determine convergence
        if( self%pp%nstates == 1 )then
            if( self%dist < self%pp%athres/5. .and.&
                self%frac > FRAC_LIM          .and.&
                self%mi_class > MI_CLASS_LIM_3D )then
                write(*,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(*,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
        else
            ! provides convergence stats for multiple states
            ! by calculating mi_joint for individual states
            allocate( state_mi_joint(self%pp%nstates), statepops(self%pp%nstates) )
            state_mi_joint = 0.
            statepops      = 0.
            do iptcl=1,self%bap%get_noris()
                istate                 = nint(self%bap%get(iptcl,'state'))
                if( istate==0 )cycle
                ! it doesn't make sense to include the state overlap here
                ! as the overall state overlap is already provided above
                state_mi_joint(istate) = state_mi_joint(istate) + self%bap%get(iptcl,'mi_class')
                state_mi_joint(istate) = state_mi_joint(istate) + self%bap%get(iptcl,'mi_inpl')
                ! 2.0 because we include two mi-values
                statepops(istate)      = statepops(istate) + 2.0
            end do
            ! normalise the overlap
            state_mi_joint     = state_mi_joint/statepops
            ! bring back the correct statepops
            statepops          = statepops/2.0
            ! the minumum overlap is in charge of convergence
            min_state_mi_joint = minval(state_mi_joint)
            ! print the overlaps and pops for the different states
            do istate=1,self%pp%nstates
                write(*,'(A,1X,I3,1X,A,1X,F7.4,1X,A,1X,I5)') '>>> STATE', istate,&
                'DISTRIBUTION OVERLAP:', state_mi_joint(istate), 'POPULATION:', nint(statepops(istate))
            end do
            if( min_state_mi_joint > MI_STATE_LIM      .and.&
                self%mi_state      > MI_STATE_LIM      .and.&
                self%dist          < self%pp%athres/5. .and.&
                self%frac          > FRAC_LIM                )then
                write(*,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(*,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
            deallocate( state_mi_joint, statepops )
        endif
    end function check_conv3D

    !>  \brief  is a getter
    real function get( self, which )
        class(convergence), intent(in) :: self
        character(len=*),   intent(in) :: which
        select case(which)
            case('corr')
                get = self%corr
            case('dist')
                get = self%dist
            case('dist_inpl')
                get = self%dist_inpl
            case('frac')
                get = self%frac
            case('mi_joint')
                get = self%mi_joint
            case('mi_class')
                get = self%mi_class
            case('mi_inpl')
                get = self%mi_inpl
            case('mi_state')
                get = self%mi_state
            case('sdev')
                get = self%sdev
            case DEFAULT
        end select
    end function get

end module simple_convergence
