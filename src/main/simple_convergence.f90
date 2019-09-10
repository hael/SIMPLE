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
    type(stats_struct) :: corr      !< correlation stats
    type(stats_struct) :: dist      !< angular distance stats
    type(stats_struct) :: dist_inpl !< in-plane angular distance stats
    type(stats_struct) :: npeaks    !< # peaks stats
    type(stats_struct) :: frac_srch !< fraction of search space scanned stats
    type(stats_struct) :: specscore !< spectral score stats
    type(stats_struct) :: spread    !< angular spread stats
    type(stats_struct) :: shwmean   !< shift increment, weighted mean stats
    type(stats_struct) :: shwstdev  !< shift increment, weighted std deviation stats
    type(stats_struct) :: pw        !< particle weight stats
    type(stats_struct) :: ow        !< max orientation weight stats
    type(oris)         :: ostats    !< centralize stats for writing
    real :: mi_class = 0.           !< class parameter distribution overlap
    real :: mi_proj  = 0.           !< projection parameter distribution overlap
    real :: mi_state = 0.           !< state parameter distribution overlap
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
        real,    allocatable :: updatecnts(:)
        logical, allocatable :: mask(:)
        real    :: avg_updatecnt
        logical :: converged
        601 format(A,1X,F8.3)
        604 format(A,1X,F8.3,1X,F8.3,1X,F8.3,1X,F8.3)
        updatecnts    = build_glob%spproj_field%get_all('updatecnt')
        avg_updatecnt = sum(updatecnts) / size(updatecnts)
        allocate(mask(size(updatecnts)), source=updatecnts > 0.5)
        call build_glob%spproj_field%stats('corr',      self%corr,      mask=mask)
        call build_glob%spproj_field%stats('specscore', self%specscore, mask=mask)
        call build_glob%spproj_field%stats('dist_inpl', self%dist_inpl, mask=mask)
        call build_glob%spproj_field%stats('frac',      self%frac_srch, mask=mask)
        call build_glob%spproj_field%stats('w',         self%pw,        mask=mask)
        self%mi_class  = build_glob%spproj_field%get_avg('mi_class',  mask=mask)
        write(logfhandle,601) '>>> CLASS OVERLAP:                          ', self%mi_class
        write(logfhandle,601) '>>> # PARTICLE UPDATES     AVG:             ', avg_updatecnt
        write(logfhandle,604) '>>> IN-PLANE DIST (DEG)    AVG/SDEV/MIN/MAX:',&
        &self%dist_inpl%avg, self%dist_inpl%sdev, self%dist_inpl%minv, self%dist_inpl%maxv
        write(logfhandle,604) '>>> PARTICLE WEIGHT        AVG/SDEV/MIN/MAX:',&
        &self%pw%avg, self%pw%sdev, self%pw%minv, self%pw%maxv
        write(logfhandle,604) '>>> % SEARCH SPACE SCANNED AVG/SDEV/MIN/MAX:',&
        &self%frac_srch%avg, self%frac_srch%sdev, self%frac_srch%minv, self%frac_srch%maxv
        write(logfhandle,604) '>>> CORRELATION            AVG/SDEV/MIN/MAX:',&
        &self%corr%avg, self%corr%sdev, self%corr%minv, self%corr%maxv
        write(logfhandle,604) '>>> SPECSCORE              AVG/SDEV/MIN/MAX:',&
        &self%specscore%avg, self%specscore%sdev, self%specscore%minv, self%specscore%maxv
        ! dynamic shift search range update
        if( self%frac_srch%avg >= FRAC_SH_LIM )then
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
                if( self%mi_class > MI_CLASS_LIM_2D_FRAC .and. self%frac_srch%avg > FRAC_LIM_FRAC )converged = .true.
            else
                if( self%mi_class > MI_CLASS_LIM_2D .and. self%frac_srch%avg > FRAC_LIM )converged = .true.
            endif
            if( converged )then
                write(logfhandle,'(A)') '>>> CONVERGED: .YES.'
            else
                write(logfhandle,'(A)') '>>> CONVERGED: .NO.'
            endif
        else
            if( self%dist_inpl%avg < 0.5 )then
                write(logfhandle,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(logfhandle,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
        endif
        ! stats
        call self%ostats%new(1)
        call self%ostats%set(1,'CLASS_OVERLAP',self%mi_class)
        call self%ostats%set(1,'PARTICLE_UPDATES',avg_updatecnt)
        call self%ostats%set(1,'IN-PLANE_DIST',self%dist_inpl%avg)
        call self%ostats%set(1,'PARTICLE_WEIGHT',self%pw%avg)
        call self%ostats%set(1,'SEARCH_SPACE_SCANNED',self%frac_srch%avg)
        call self%ostats%set(1,'CORRELATION',self%corr%avg)
        call self%ostats%set(1,'SPECSCORE',self%specscore%avg)
        call self%ostats%write(STATS_FILE)
        call self%ostats%kill
    end function check_conv2D

    function check_conv3D( self, cline, msk ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        real,               intent(in)    :: msk
        real,    allocatable :: state_mi_joint(:), statepops(:), updatecnts(:)
        logical, allocatable :: mask(:)
        real    :: min_state_mi_joint, avg_updatecnt
        logical :: converged
        integer :: iptcl, istate
        601 format(A,1X,F8.3)
        604 format(A,1X,F8.3,1X,F8.3,1X,F8.3,1X,F8.3)
        updatecnts = build_glob%spproj_field%get_all('updatecnt')
        avg_updatecnt = sum(updatecnts) / size(updatecnts)
        allocate(mask(size(updatecnts)), source=updatecnts > 0.5)
        call build_glob%spproj_field%stats('corr',      self%corr,      mask=mask)
        call build_glob%spproj_field%stats('specscore', self%specscore, mask=mask)
        call build_glob%spproj_field%stats('dist',      self%dist,      mask=mask)
        call build_glob%spproj_field%stats('dist_inpl', self%dist_inpl, mask=mask)
        call build_glob%spproj_field%stats('npeaks',    self%npeaks,    mask=mask)
        call build_glob%spproj_field%stats('frac',      self%frac_srch, mask=mask)
        call build_glob%spproj_field%stats('w',         self%pw,        mask=mask)
        call build_glob%spproj_field%stats('ow',        self%ow,        mask=mask)
        call build_glob%spproj_field%stats('spread',    self%spread,    mask=mask)
        call build_glob%spproj_field%stats('shwmean',   self%shwmean,   mask=mask)
        call build_glob%spproj_field%stats('shwstdev',  self%shwstdev,  mask=mask)
        self%mi_proj   = build_glob%spproj_field%get_avg('mi_proj',   mask=mask)
        self%mi_state  = build_glob%spproj_field%get_avg('mi_state',  mask=mask)
        write(logfhandle,601) '>>> ORIENTATION OVERLAP:                      ', self%mi_proj
        if( params_glob%nstates > 1 )then
        write(logfhandle,601) '>>> STATE OVERLAP:                            ', self%mi_state
        endif
        write(logfhandle,601) '>>> # PARTICLE UPDATES       AVG:             ', avg_updatecnt
        write(logfhandle,604) '>>> DIST BTW BEST ORIS (DEG) AVG/SDEV/MIN/MAX:',&
        &self%dist%avg, self%dist%sdev, self%dist%minv, self%dist%maxv
        write(logfhandle,604) '>>> IN-PLANE DIST      (DEG) AVG/SDEV/MIN/MAX:',&
        &self%dist_inpl%avg, self%dist_inpl%sdev, self%dist_inpl%minv, self%dist_inpl%maxv
        write(logfhandle,604) '>>> ANGULAR SPREAD     (DEG) AVG/SDEV/MIN/MAX:',&
        &self%spread%avg, self%spread%sdev, self%spread%minv, self%spread%maxv
        write(logfhandle,604) '>>> # PEAKS                  AVG/SDEV/MIN/MAX:',&
        &self%npeaks%avg, self%npeaks%sdev, self%npeaks%minv, self%npeaks%maxv
        write(logfhandle,604) '>>> PARTICLE WEIGHT          AVG/SDEV/MIN/MAX:',&
        &self%pw%avg, self%pw%sdev, self%pw%minv, self%pw%maxv
        write(logfhandle,604) '>>> ORIENTATION WEIGHT MAX   AVG/SDEV/MIN/MAX:',&
        &self%ow%avg, self%ow%sdev, self%ow%minv, self%ow%maxv
        write(logfhandle,604) '>>> % SEARCH SPACE SCANNED   AVG/SDEV/MIN/MAX:',&
        &self%frac_srch%avg, self%frac_srch%sdev, self%frac_srch%minv, self%frac_srch%maxv
        write(logfhandle,604) '>>> CORRELATION              AVG/SDEV/MIN/MAX:',&
        &self%corr%avg, self%corr%sdev, self%corr%minv, self%corr%maxv
        write(logfhandle,604) '>>> SPECSCORE                AVG/SDEV/MIN/MAX:',&
        &self%specscore%avg, self%specscore%sdev, self%specscore%minv, self%specscore%maxv
        write(logfhandle,604) '>>> AVG  SHIFT INCR PEAKS    AVG/SDEV/MIN/MAX:',&
        &self%shwmean%avg, self%shwmean%sdev, self%shwmean%minv, self%shwmean%maxv
        write(logfhandle,604) '>>> SDEV SHIFT INCR PEAKS    AVG/SDEV/MIN/MAX:',&
        &self%shwstdev%avg, self%shwstdev%sdev, self%shwstdev%minv, self%shwstdev%maxv
        ! dynamic shift search range update
        if( self%frac_srch%avg >= FRAC_SH_LIM )then
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
            if( self%frac_srch%avg > FRAC_LIM .and.&
                self%mi_proj  > MI_CLASS_LIM_3D )then
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
                self%frac_srch%avg      > FRAC_LIM      )then
                write(logfhandle,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(logfhandle,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
            deallocate( state_mi_joint, statepops )
        endif
        ! stats
        call self%ostats%new(1)
        call self%ostats%set(1,'ORIENTATION_OVERLAP',self%mi_proj)
        if( params_glob%nstates > 1 ) call self%ostats%set(1,'STATE_OVERLAP', self%mi_state)
        call self%ostats%set(1,'PARTICLE_UPDATES',avg_updatecnt)
        call self%ostats%set(1,'DIST_BTW_BEST_ORIS',self%dist%avg)
        call self%ostats%set(1,'IN-PLANE_DIST',self%dist_inpl%avg)
        call self%ostats%set(1,'ANGULAR_SPREAD',self%spread%avg)
        call self%ostats%set(1,'#_PEAKS',self%npeaks%avg)
        call self%ostats%set(1,'PARTICLE_WEIGHT',self%pw%avg)
        call self%ostats%set(1,'ORIENTATION_WEIGHT_MAX',self%ow%avg)
        call self%ostats%set(1,'SEARCH_SPACE_SCANNED',self%frac_srch%avg)
        call self%ostats%set(1,'CORRELATION',self%corr%avg)
        call self%ostats%set(1,'SPECSCORE',self%specscore%avg)
        call self%ostats%set(1,'AVG_SHIFT_INCR_PEAKS',self%shwmean%avg)
        call self%ostats%set(1,'SDEV_SHIFT_INCR_PEAKS',self%shwstdev%avg)
        call self%ostats%write(STATS_FILE)
        ! destruct
        deallocate(mask, updatecnts)
        call self%ostats%kill
    end function check_conv3D

    function check_conv_cluster( self, cline ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        integer, allocatable :: statepops(:)
        logical :: converged
        integer :: istate
        601 format(A,1X,F8.3)
        604 format(A,1X,F8.3,1X,F8.3,1X,F8.3,1X,F8.3)
        call build_glob%spproj_field%stats('frac', self%frac_srch)
        self%mi_state  = build_glob%spproj_field%get_avg('mi_state')
        write(logfhandle,601) '>>> STATE OVERLAP:                          ', self%mi_state
        write(logfhandle,604) '>>> % SEARCH SPACE SCANNED AVG/SDEV/MIN/MAX:',&
        &self%frac_srch%avg, self%frac_srch%sdev, self%frac_srch%minv, self%frac_srch%maxv
        ! provides convergence stats for multiple states
        ! by calculating mi_joint for individual states
        call build_glob%spproj_field%get_pops(statepops,'state')
        call build_glob%spproj_field%stats('corr', self%corr)
        write(logfhandle,604) '>>> CORRELATION            AVG/SDEV/MIN/MAX:',&
        &self%corr%avg, self%corr%sdev, self%corr%minv, self%corr%maxv
        ! print the overlaps and pops for the different states
        do istate=1,params_glob%nstates
            write(logfhandle,'(A,I2,1X,A,1X,I8)') '>>> STATE ',istate,'POPULATION:', statepops(istate)
        end do
        if( self%mi_state > HET_MI_STATE_LIM .and.&
            self%frac_srch%avg > HET_FRAC_LIM     )then
            write(logfhandle,'(A)') '>>> CONVERGED: .YES.'
            converged = .true.
        else
            write(logfhandle,'(A)') '>>> CONVERGED: .NO.'
            converged = .false.
        endif
        ! stats
        call self%ostats%new(1)
        call self%ostats%set(1,'STATE_OVERLAP',       self%mi_state)
        call self%ostats%set(1,'SEARCH_SPACE_SCANNED',self%frac_srch%avg)
        call self%ostats%set(1,'CORRELATION',         self%corr%avg)
        do istate=1,params_glob%nstates
            call self%ostats%set(1,'STATE_POPULATION_'//int2str(istate), real(statepops(istate)))
        enddo
        call self%ostats%write(STATS_FILE)
        ! cleanup
        call self%ostats%kill
        deallocate( statepops )
    end function check_conv_cluster

    real function get( self, which )
        class(convergence), intent(in) :: self
        character(len=*),   intent(in) :: which
        get = 0.
        select case(which)
            case('corr')
                get = self%corr%avg
            case('dist')
                get = self%dist%avg
            case('dist_inpl')
                get = self%dist_inpl%avg
            case('frac_srch')
                get = self%frac_srch%avg
            case('mi_class')
                get = self%mi_class
            case('mi_proj')
                get = self%mi_proj
            case('mi_state')
                get = self%mi_state
            case('spread')
                get = self%spread%avg
        end select
    end function get

end module simple_convergence
