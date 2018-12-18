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
    real :: mi_class  = 0. !< class parameter distribution overlap
    real :: mi_proj   = 0. !< projection parameter distribution overlap
    real :: mi_state  = 0. !< state parameter distribution overlap
    real :: spread    = 0. !< angular spread
    real :: shwmean   = 0. !< shift increment, weighted mean
    real :: shwstdev  = 0. !< shift increment, weighted std deviation
    real :: bfac      = 0. !< average per-particle B-factor (search)
    real :: bfac_rec  = 0. !< average per-particle B-factor (rec)
  contains
    procedure :: check_conv2D
    procedure :: check_conv3D
    procedure :: check_conv_cluster
    procedure :: get
end type convergence

contains

    function check_conv2D( self, cline, ncls, msk ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        integer,            intent(in)    :: ncls
        real,               intent(in)    :: msk
        real,    allocatable :: tmp_arr(:)
        logical, allocatable :: mask(:)
        real    :: avg_updatecnt
        logical :: converged
        ! generate mask
        allocate(mask(build_glob%spproj_field%get_noris()))
        tmp_arr = build_glob%spproj_field%get_all('updatecnt')
        mask    = tmp_arr > 0.5
        tmp_arr = build_glob%spproj_field%get_all('state')
        mask    = mask .and. tmp_arr > 0.5
        deallocate(tmp_arr)
        ! stats
        avg_updatecnt  = build_glob%spproj_field%get_avg('updatecnt', mask=mask)
        self%corr      = build_glob%spproj_field%get_avg('corr',      mask=mask)
        self%dist_inpl = build_glob%spproj_field%get_avg('dist_inpl', mask=mask)
        self%frac      = build_glob%spproj_field%get_avg('frac',      mask=mask)
        self%bfac      = build_glob%spproj_field%get_avg('bfac',      mask=mask)
        self%mi_class  = build_glob%spproj_field%get_avg('mi_class',  mask=mask)
        write(logfhandle,'(A,1X,F7.4)') '>>> CLASS OVERLAP:            ', self%mi_class
        write(logfhandle,'(A,1X,F7.1)') '>>> AVG # PARTICLE UPDATES:   ', avg_updatecnt
        write(logfhandle,'(A,1X,F7.1)') '>>> AVG IN-PLANE DIST (DEG):  ', self%dist_inpl
        write(logfhandle,'(A,1X,F7.1)') '>>> AVG PER-PARTICLE B-FACTOR:', self%bfac
        write(logfhandle,'(A,1X,F7.1)') '>>> % SEARCH SPACE SCANNED:   ', self%frac
        write(logfhandle,'(A,1X,F7.4)') '>>> CORRELATION:              ', self%corr
        ! dynamic shift search range update
        if( self%frac >= FRAC_SH_LIM )then
            if( .not. cline%defined('trs') .or. params_glob%trs <  MINSHIFT )then
                ! determine shift bounds
                params_glob%trs = MSK_FRAC*msk
                params_glob%trs = max(MINSHIFT,params_glob%trs)
                params_glob%trs = min(MAXSHIFT,params_glob%trs)
                ! set shift search flag
                params_glob%l_doshift = .true.
            endif
        endif
        ! determine convergence
        if( ncls > 1 )then
            converged = .false.
            if( (params_glob%l_frac_update) .or. (params_glob%stream.eq.'yes') )then
                if( self%mi_class > MI_CLASS_LIM_2D_FRAC .and. self%frac > FRAC_LIM_FRAC )converged = .true.
            else
                if( self%mi_class > MI_CLASS_LIM_2D .and. self%frac > FRAC_LIM )converged = .true.
            endif
            if( converged )then
                write(logfhandle,'(A)') '>>> CONVERGED: .YES.'
            else
                write(logfhandle,'(A)') '>>> CONVERGED: .NO.'
            endif
        else
            if( self%dist_inpl < 0.5 )then
                write(logfhandle,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(logfhandle,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
        endif
    end function check_conv2D

    function check_conv3D( self, cline, msk ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        real,               intent(in)    :: msk
        real,    allocatable :: state_mi_joint(:), statepops(:), updatecnts(:)
        logical, allocatable :: mask(:)
        real    :: min_state_mi_joint, avg_updatecnt
        logical :: converged, bfac_rec_there
        integer :: iptcl, istate
        updatecnts = build_glob%spproj_field%get_all('updatecnt')
        avg_updatecnt = sum(updatecnts) / size(updatecnts)
        allocate(mask(size(updatecnts)), source=updatecnts > 0.5)
        self%corr      = build_glob%spproj_field%get_avg('corr',      mask=mask)
        self%dist      = build_glob%spproj_field%get_avg('dist',      mask=mask)
        self%dist_inpl = build_glob%spproj_field%get_avg('dist_inpl', mask=mask)
        self%npeaks    = build_glob%spproj_field%get_avg('npeaks',    mask=mask)
        self%frac      = build_glob%spproj_field%get_avg('frac',      mask=mask)
        self%mi_proj   = build_glob%spproj_field%get_avg('mi_proj',   mask=mask)
        self%mi_state  = build_glob%spproj_field%get_avg('mi_state',  mask=mask)
        self%spread    = build_glob%spproj_field%get_avg('spread',    mask=mask)
        self%shwmean   = build_glob%spproj_field%get_avg('shwmean',   mask=mask)
        self%shwstdev  = build_glob%spproj_field%get_avg('shwstdev',  mask=mask)
        self%bfac      = build_glob%spproj_field%get_avg('bfac')     ! always updated for all ptcls with states > 0
        self%bfac_rec  = 0.
        bfac_rec_there = build_glob%spproj_field%isthere('bfac_rec')
        if( bfac_rec_there ) self%bfac_rec = build_glob%spproj_field%get_avg('bfac_rec')
        write(logfhandle,'(A,1X,F7.4)') '>>> ORIENTATION OVERLAP:               ', self%mi_proj
        write(logfhandle,'(A,1X,F7.1)') '>>> AVG # PARTICLE UPDATES:            ', avg_updatecnt
        if( params_glob%nstates > 1 )then
        write(logfhandle,'(A,1X,F7.4)') '>>> STATE OVERLAP:                     ', self%mi_state
        endif
        write(logfhandle,'(A,1X,F7.1)') '>>> AVG DIST BTW BEST ORIS (DEG):      ', self%dist
        write(logfhandle,'(A,1X,F7.1)') '>>> AVG IN-PLANE DIST      (DEG):      ', self%dist_inpl
        write(logfhandle,'(A,1X,F7.1)') '>>> AVG # PEAKS:                       ', self%npeaks
        write(logfhandle,'(A,1X,F7.1)') '>>> AVG PER-PARTICLE B-FACTOR (SEARCH):', self%bfac
        if( bfac_rec_there )then
        write(logfhandle,'(A,1X,F7.1)') '>>> AVG PER-PARTICLE B-FACTOR (REC):   ', self%bfac_rec
        endif
        write(logfhandle,'(A,1X,F7.1)') '>>> % SEARCH SPACE SCANNED:            ', self%frac
        write(logfhandle,'(A,1X,F7.4)') '>>> CORRELATION:                       ', self%corr
        write(logfhandle,'(A,1X,F7.2)') '>>> ANGULAR SPREAD (DEG):              ', self%spread
        write(logfhandle,'(A,1X,F7.2)') '>>> AVG WEIGHTED SHIFT INCREMENT:      ', self%shwmean
        write(logfhandle,'(A,1X,F7.2)') '>>> AVG WEIGHTED SHIFT INCR STDEV:     ', self%shwstdev
        ! dynamic shift search range update
        if( self%frac >= FRAC_SH_LIM )then
            if( .not. cline%defined('trs') .or. &
                & params_glob%trs <  MINSHIFT )then
                ! determine shift bounds
                params_glob%trs = MSK_FRAC*msk
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
                write(logfhandle,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(logfhandle,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
        else
            ! provides convergence stats for multiple states
            ! by calculating mi_joint for individual states
            allocate( state_mi_joint(params_glob%nstates), statepops(params_glob%nstates) )
            state_mi_joint = 0.
            statepops      = 0.
            do iptcl=1,build_glob%spproj_field%get_noris()
                istate = build_glob%spproj_field%get_state(iptcl)
                if( istate==0 )cycle
                state_mi_joint(istate) = state_mi_joint(istate) + build_glob%spproj_field%get(iptcl,'mi_proj')
                statepops(istate)      = statepops(istate) + 1.
            end do
            ! normalise the overlap
            forall( istate=1:params_glob%nstates, statepops(istate)>0. )&
               &state_mi_joint(istate) = state_mi_joint(istate)/statepops(istate)
            ! the minumum overlap is in charge of convergence
            min_state_mi_joint = minval(state_mi_joint, mask=statepops>0.)
            ! print the overlaps and pops for the different states
            do istate=1,params_glob%nstates
                write(logfhandle,'(A,1X,I3,1X,A,1X,F7.4,1X,A,1X,I8)') '>>> STATE', istate,&
                'JOINT DISTRIBUTION OVERLAP:', state_mi_joint(istate), 'POPULATION:', nint(statepops(istate))
            end do
            if( min_state_mi_joint > MI_STATE_JOINT_LIM .and.&
                self%mi_state      > MI_STATE_LIM .and.&
                self%frac          > FRAC_LIM      )then
                write(logfhandle,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(logfhandle,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
            deallocate( state_mi_joint, statepops )
        endif
        ! destruct
        deallocate(mask, updatecnts)
    end function check_conv3D

    function check_conv_cluster( self, cline ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        real, allocatable :: statepops(:)
        logical           :: converged
        integer           :: iptcl, istate
        self%frac      = build_glob%spproj_field%get_avg('frac')
        self%mi_state  = build_glob%spproj_field%get_avg('mi_state')
        write(logfhandle,'(A,1X,F7.4)') '>>> STATE OVERLAP:                ', self%mi_state
        write(logfhandle,'(A,1X,F7.1)') '>>> % SEARCH SPACE SCANNED:       ', self%frac
        ! provides convergence stats for multiple states
        ! by calculating mi_joint for individual states
        allocate( statepops(params_glob%nstates) )
        statepops = 0.
        do iptcl=1,build_glob%spproj_field%get_noris()
            istate = build_glob%spproj_field%get_state(iptcl)
            if( istate==0 )cycle
            statepops(istate) = statepops(istate) + 1.0
        end do
        if( build_glob%spproj_field%isthere('bfac') )then
            self%bfac = build_glob%spproj_field%get_avg('bfac')
            write(logfhandle,'(A,1X,F6.1)') '>>> AVERAGE PER-PARTICLE B-FACTOR:', self%bfac
        endif
        self%corr = build_glob%spproj_field%get_avg('corr')
        write(logfhandle,'(A,1X,F7.4)') '>>> CORRELATION                  :', self%corr
        ! print the overlaps and pops for the different states
        do istate=1,params_glob%nstates
            write(logfhandle,'(A,I2,1X,A,1X,I8)') '>>> STATE ',istate,'POPULATION:', nint(statepops(istate))
        end do
        if( self%mi_state > HET_MI_STATE_LIM .and.&
            self%frac     > HET_FRAC_LIM     )then
            write(logfhandle,'(A)') '>>> CONVERGED: .YES.'
            converged = .true.
        else
            write(logfhandle,'(A)') '>>> CONVERGED: .NO.'
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
            case('mi_class')
                get = self%mi_class
            case('mi_proj')
                get = self%mi_proj
            case('mi_state')
                get = self%mi_state
            case('spread')
                get = self%spread
            case('bfac')
                get = self%bfac
            case('bfac_rec')
                get = self%bfac_rec
            case DEFAULT
                get = 0.
        end select
    end function get

end module simple_convergence
