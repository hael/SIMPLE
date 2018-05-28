! for checking convergence
module simple_convergence
include 'simple_lib.f08'
use simple_oris,       only: oris
use simple_parameters, only: params_glob
use simple_builder,    only: build_glob
use simple_cmdline,    only: cmdline
implicit none

public :: convergence
private

type convergence
    private
    real :: corr      = 0. !< average correlation
    real :: dist      = 0. !< average angular distance
    real :: dist_inpl = 0. !< average in-plane angular distance
    real :: npeaks    = 0. !< average # peaks
    real :: frac      = 0. !< average fraction of search space scanned
    real :: mi_joint  = 0. !< joint parameter distribution overlap
    real :: mi_class  = 0. !< class parameter distribution overlap
    real :: mi_proj   = 0. !< projection parameter distribution overlap
    real :: mi_inpl   = 0. !< in-plane parameter distribution overlap
    real :: mi_state  = 0. !< state parameter distribution overlap
    real :: sdev      = 0. !< angular standard deviation of model
    real :: bfac      = 0. !< average per-particle B-factor
  contains
    procedure :: check_conv2D
    procedure :: check_conv3D
    procedure :: check_conv_cluster
    procedure :: get
end type convergence

contains

    function check_conv2D( self, cline, ncls ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        integer,            intent(in)    :: ncls
        real,    allocatable :: updatecnts(:)
        logical, allocatable :: mask(:)
        real    :: avg_updatecnt
        logical :: converged, update_frac
        update_frac   = .false.
        avg_updatecnt = 0.
        if( build_glob%spproj_field%isthere('updatecnt') )then
            updatecnts    = build_glob%spproj_field%get_all('updatecnt')
            avg_updatecnt = sum(updatecnts) / real(size(updatecnts))
            update_frac   = count(updatecnts > 0.5) > 0 ! for the case of greedy step with frac_update on
        endif
        if( update_frac )then
            ! fractional particle update
            allocate(mask(size(updatecnts)), source=updatecnts > 0.5)
            self%corr      = build_glob%spproj_field%get_avg('corr',      mask=mask)
            self%dist_inpl = build_glob%spproj_field%get_avg('dist_inpl', mask=mask)
            self%frac      = build_glob%spproj_field%get_avg('frac',      mask=mask)
            self%bfac      = build_glob%spproj_field%get_avg('bfac',      mask=mask)
            self%mi_joint  = build_glob%spproj_field%get_avg('mi_joint',  mask=mask)
            self%mi_inpl   = build_glob%spproj_field%get_avg('mi_inpl',   mask=mask)
            self%mi_class  = build_glob%spproj_field%get_avg('mi_class',  mask=mask)
        else
            self%corr      = build_glob%spproj_field%get_avg('corr')
            self%dist_inpl = build_glob%spproj_field%get_avg('dist_inpl')
            self%frac      = build_glob%spproj_field%get_avg('frac')
            self%bfac      = build_glob%spproj_field%get_avg('bfac')
            self%mi_joint  = build_glob%spproj_field%get_avg('mi_joint')
            self%mi_class  = build_glob%spproj_field%get_avg('mi_class')
            self%mi_inpl   = build_glob%spproj_field%get_avg('mi_inpl')
        endif
        write(*,'(A,1X,F7.4)') '>>> JOINT    DISTRIBUTION OVERLAP:     ', self%mi_joint
        write(*,'(A,1X,F7.4)') '>>> CLASS    DISTRIBUTION OVERLAP:     ', self%mi_class
        write(*,'(A,1X,F7.4)') '>>> IN-PLANE DISTRIBUTION OVERLAP:     ', self%mi_inpl
        write(*,'(A,1X,F7.1)') '>>> AVERAGE # PARTICLE UPDATES:        ', avg_updatecnt
        write(*,'(A,1X,F7.1)') '>>> AVERAGE IN-PLANE ANGULAR DISTANCE: ', self%dist_inpl
        write(*,'(A,1X,F7.1)') '>>> AVERAGE PER-PARTICLE B-FACTOR:     ', self%bfac
        write(*,'(A,1X,F7.1)') '>>> PERCENTAGE OF SEARCH SPACE SCANNED:', self%frac
        write(*,'(A,1X,F7.4)') '>>> CORRELATION:                       ', self%corr
        ! dynamic shift search range update
        if( self%frac >= FRAC_SH_LIM )then
            if( .not. cline%defined('trs') .or. params_glob%trs <  MINSHIFT )then
                ! determine shift bounds
                params_glob%trs = MSK_FRAC*real(params_glob%msk)
                params_glob%trs = max(MINSHIFT,params_glob%trs)
                params_glob%trs = min(MAXSHIFT,params_glob%trs)
                ! set shift search flag
                params_glob%l_doshift = .true.
            endif
        endif
        ! determine convergence
        if( ncls > 1 )then
            converged = .false.
            if( params_glob%l_frac_update )then
                if( self%mi_joint > MI_CLASS_LIM_2D_FRAC .and. self%frac > FRAC_LIM_FRAC )converged = .true.
            else
                if( self%mi_class > MI_CLASS_LIM_2D .and. self%frac > FRAC_LIM )converged = .true.
            endif
            if( converged )then
                write(*,'(A)') '>>> CONVERGED: .YES.'
            else
                write(*,'(A)') '>>> CONVERGED: .NO.'
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

    function check_conv3D( self, cline ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        real,    allocatable :: state_mi_joint(:), statepops(:), updatecnts(:)
        logical, allocatable :: mask(:)
        real    :: min_state_mi_joint, avg_updatecnt
        logical :: converged, update_frac
        integer :: iptcl, istate
        update_frac   = .false.
        avg_updatecnt = 0.
        if( build_glob%spproj_field%isthere('updatecnt') )then
            updatecnts    = build_glob%spproj_field%get_all('updatecnt')
            avg_updatecnt = sum(updatecnts) / real(size(updatecnts))
            update_frac   = count(updatecnts > 0.5) > 0
        endif
        if( update_frac )then
            ! fractional particle update
            allocate(mask(size(updatecnts)), source=updatecnts > 0.5)
            self%corr      = build_glob%spproj_field%get_avg('corr',      mask=mask)
            self%dist      = build_glob%spproj_field%get_avg('dist',      mask=mask)
            self%dist_inpl = build_glob%spproj_field%get_avg('dist_inpl', mask=mask)
            self%npeaks    = build_glob%spproj_field%get_avg('npeaks',    mask=mask)
            self%frac      = build_glob%spproj_field%get_avg('frac',      mask=mask)
            self%mi_joint  = build_glob%spproj_field%get_avg('mi_joint',  mask=mask)
            self%mi_proj   = build_glob%spproj_field%get_avg('mi_proj',   mask=mask)
            self%mi_inpl   = build_glob%spproj_field%get_avg('mi_inpl',   mask=mask)
            self%mi_state  = build_glob%spproj_field%get_avg('mi_state',  mask=mask)
            self%sdev      = build_glob%spproj_field%get_avg('sdev',      mask=mask)
            self%bfac      = build_glob%spproj_field%get_avg('bfac',      mask=mask)
        else
            self%corr      = build_glob%spproj_field%get_avg('corr')
            self%dist      = build_glob%spproj_field%get_avg('dist')
            self%dist_inpl = build_glob%spproj_field%get_avg('dist_inpl')
            self%npeaks    = build_glob%spproj_field%get_avg('npeaks')
            self%frac      = build_glob%spproj_field%get_avg('frac')
            self%mi_joint  = build_glob%spproj_field%get_avg('mi_joint')
            self%mi_proj   = build_glob%spproj_field%get_avg('mi_proj')
            self%mi_inpl   = build_glob%spproj_field%get_avg('mi_inpl')
            self%mi_state  = build_glob%spproj_field%get_avg('mi_state')
            self%sdev      = build_glob%spproj_field%get_avg('sdev')
            self%bfac      = build_glob%spproj_field%get_avg('bfac')
        endif
        write(*,'(A,1X,F7.4)') '>>> JOINT    DISTRIBUTION OVERLAP:     ', self%mi_joint
        write(*,'(A,1X,F7.4)') '>>> PROJ     DISTRIBUTION OVERLAP:     ', self%mi_proj
        write(*,'(A,1X,F7.4)') '>>> IN-PLANE DISTRIBUTION OVERLAP:     ', self%mi_inpl
        write(*,'(A,1X,F7.1)') '>>> AVERAGE # PARTICLE UPDATES:        ', avg_updatecnt
        if( params_glob%nstates > 1 )&
        &write(*,'(A,1X,F7.4)') '>>> STATE DISTRIBUTION OVERLAP:        ', self%mi_state
        write(*,'(A,1X,F7.1)') '>>> AVERAGE ANGULAR DISTANCE BTW ORIS: ', self%dist
        write(*,'(A,1X,F7.1)') '>>> AVERAGE IN-PLANE ANGULAR DISTANCE: ', self%dist_inpl
        write(*,'(A,1X,F7.1)') '>>> AVERAGE # PEAKS:                   ', self%npeaks
        write(*,'(A,1X,F7.1)') '>>> AVERAGE PER-PARTICLE B-FACTOR:     ', self%bfac
        write(*,'(A,1X,F7.1)') '>>> PERCENTAGE OF SEARCH SPACE SCANNED:', self%frac
        write(*,'(A,1X,F7.4)') '>>> CORRELATION:                       ', self%corr
        write(*,'(A,1X,F7.2)') '>>> ANGULAR SDEV OF MODEL:             ', self%sdev
        ! dynamic shift search range update
        if( self%frac >= FRAC_SH_LIM )then
            if( .not. cline%defined('trs') .or. &
                & params_glob%trs <  MINSHIFT )then
                ! determine shift bounds
                params_glob%trs = MSK_FRAC*real(params_glob%msk)
                params_glob%trs = max(MINSHIFT,params_glob%trs)
                params_glob%trs = min(MAXSHIFT,params_glob%trs)
                ! set shift search flag
                params_glob%l_doshift = .true.
            endif
        endif
        ! determine convergence
        if( params_glob%nstates == 1 )then
            if( self%frac    > FRAC_LIM       .and.&
                self%mi_proj > MI_CLASS_LIM_3D )then
                write(*,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(*,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
        else
            ! provides convergence stats for multiple states
            ! by calculating mi_joint for individual states
            allocate( state_mi_joint(params_glob%nstates), statepops(params_glob%nstates) )
            state_mi_joint = 0.
            statepops      = 0.
            do iptcl=1,build_glob%spproj_field%get_noris()
                istate = nint(build_glob%spproj_field%get(iptcl,'state'))
                if( istate==0 )cycle
                ! it doesn't make sense to include the state overlap here
                ! as the overall state overlap is already provided above
                state_mi_joint(istate) = state_mi_joint(istate) + build_glob%spproj_field%get(iptcl,'mi_proj')
                state_mi_joint(istate) = state_mi_joint(istate) + build_glob%spproj_field%get(iptcl,'mi_inpl')
                ! 2.0 because we include two mi-values
                statepops(istate)      = statepops(istate) + 2.0
            end do
            ! normalise the overlap
            forall( istate=1:params_glob%nstates, statepops(istate)>0. )&
               &state_mi_joint(istate) = state_mi_joint(istate)/statepops(istate)
            ! bring back the correct statepops
            statepops = statepops/2.0
            ! the minumum overlap is in charge of convergence
            min_state_mi_joint = minval(state_mi_joint, MASK=statepops>0.)
            ! print the overlaps and pops for the different states
            do istate=1,params_glob%nstates
                write(*,'(A,1X,I3,1X,A,1X,F7.4,1X,A,1X,I8)') '>>> STATE', istate,&
                'DISTRIBUTION OVERLAP:', state_mi_joint(istate), 'POPULATION:', nint(statepops(istate))
            end do
            if( min_state_mi_joint > MI_STATE_LIM .and.&
                self%mi_state      > MI_STATE_LIM .and.&
                self%frac          > FRAC_LIM      )then
                write(*,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(*,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
            deallocate( state_mi_joint, statepops )
        endif
    end function check_conv3D

    function check_conv_cluster( self, cline ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        real, allocatable :: statepops(:)
        logical           :: converged
        integer           :: iptcl, istate
        self%frac      = build_glob%spproj_field%get_avg('frac')
        self%mi_state  = build_glob%spproj_field%get_avg('mi_state')
        write(*,'(A,1X,F7.4)') '>>> STATE DISTRIBUTION OVERLAP:        ', self%mi_state
        write(*,'(A,1X,F7.1)') '>>> PERCENTAGE OF SEARCH SPACE SCANNED:', self%frac
        ! provides convergence stats for multiple states
        ! by calculating mi_joint for individual states
        allocate( statepops(params_glob%nstates) )
        statepops = 0.
        do iptcl=1,build_glob%spproj_field%get_noris()
            istate = nint(build_glob%spproj_field%get(iptcl,'state'))
            if( istate==0 )cycle
            statepops(istate) = statepops(istate) + 1.0
        end do
        if( build_glob%spproj_field%isthere('bfac') )then
            self%bfac = build_glob%spproj_field%get_avg('bfac')
            write(*,'(A,1X,F6.1)') '>>> AVERAGE PER-PARTICLE B-FACTOR: ', self%bfac
        endif
        self%corr = build_glob%spproj_field%get_avg('corr')
        write(*,'(A,1X,F7.4)') '>>> CORRELATION                  :', self%corr
        ! print the overlaps and pops for the different states
        do istate=1,params_glob%nstates
            write(*,'(A,I2,1X,A,1X,I8)') '>>> STATE ',istate,'POPULATION:', nint(statepops(istate))
        end do
        if( self%mi_state > HET_MI_STATE_LIM .and.&
            self%frac     > HET_FRAC_LIM     )then
            write(*,'(A)') '>>> CONVERGED: .YES.'
            converged = .true.
        else
            write(*,'(A)') '>>> CONVERGED: .NO.'
            converged = .false.
        endif
        deallocate( statepops )
    end function check_conv_cluster

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
            case('mi_proj')
                get = self%mi_proj
            case('mi_inpl')
                get = self%mi_inpl
            case('mi_state')
                get = self%mi_state
            case('sdev')
                get = self%sdev
            case('bfac')
                get = self%bfac
            case DEFAULT
        end select
    end function get

end module simple_convergence
