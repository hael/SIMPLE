! for checking convergence
module simple_convergence
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_builder,    only: build_glob
use simple_cmdline,    only: cmdline
use simple_progress
implicit none

public :: convergence
private

type convergence
    private
    type(stats_struct) :: corr       !< correlation stats
    type(stats_struct) :: dist       !< angular distance stats
    type(stats_struct) :: dist_inpl  !< in-plane angular distance stats
    type(stats_struct) :: dist_peaks !< average angular distance between peaks
    type(stats_struct) :: frac_srch  !< fraction of search space scanned stats
    type(stats_struct) :: frac_sh    !< fraction of search space scanned stats, shifts
    type(stats_struct) :: shincarg   !< shift increment
    type(stats_struct) :: pw         !< particle weight stats
    type(stats_struct) :: nevals     !< # cost function evaluations
    type(stats_struct) :: ngevals    !< # gradient evaluations
    type(stats_struct) :: better     !< improvement statistics
    type(stats_struct) :: better_l   !< improvement statistics, LBFGS-B
    type(stats_struct) :: npeaks     !< peak statistics
    type(stats_struct) :: cc_peak    !< cc peak statistics
    type(stats_struct) :: cc_nonpeak !< cc non-peak statistics
    type(oris)         :: ostats     !< centralize stats for writing
    real :: mi_class = 0.            !< class parameter distribution overlap
    real :: mi_proj  = 0.            !< projection parameter distribution overlap
    real :: mi_state = 0.            !< state parameter distribution overlap
    real :: progress = 0.            !< progress estimation
  contains
    procedure :: check_conv2D
    procedure :: check_conv3D
    procedure :: check_conv3Dc
    procedure :: check_conv_cluster
    procedure :: get
end type convergence

contains

    function check_conv2D( self, cline, os, ncls, msk ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(oris),        intent(inout) :: os
        integer,            intent(in)    :: ncls
        real,               intent(in)    :: msk
        real,    allocatable :: updatecnts(:), states(:), corrs(:), pws(:)
        logical, allocatable :: mask(:)
        real    :: avg_updatecnt, overlap_lim, fracsrch_lim, corr_t, percen_nonzero_pw
        logical :: converged, chk4conv
        601 format(A,1X,F8.3)
        604 format(A,1X,F8.3,1X,F8.3,1X,F8.3,1X,F8.3)
        states            = build_glob%spproj_field%get_all('state')
        corrs             = build_glob%spproj_field%get_all('corr')
        updatecnts        = build_glob%spproj_field%get_all('updatecnt')
        avg_updatecnt     = sum(updatecnts) / real(count(states > 0.5))
        allocate(mask(size(updatecnts)), source=updatecnts > 0.5 .and. states > 0.5)
        pws               = build_glob%spproj_field%get_all('w')
        percen_nonzero_pw = (real(count(mask .and. (pws > TINY))) / real(count(mask))) * 100.
        call os%stats('corr',      self%corr,      mask=mask)
        call os%stats('dist_inpl', self%dist_inpl, mask=mask)
        call os%stats('frac',      self%frac_srch, mask=mask)
        call os%stats('shincarg',  self%shincarg,  mask=mask)
        call os%stats('w',         self%pw,        mask=mask)
        self%mi_class = os%get_avg('mi_class',    mask=mask)
        corr_t        = self%corr%avg - 2. * self%corr%sdev
        ! overlaps and particle updates
        write(logfhandle,601) '>>> CLASS OVERLAP:                          ', self%mi_class
        write(logfhandle,601) '>>> # PARTICLE UPDATES     AVG:             ', avg_updatecnt
        ! dists and % search space
        write(logfhandle,604) '>>> IN-PLANE DIST    (DEG) AVG/SDEV/MIN/MAX:', self%dist_inpl%avg, self%dist_inpl%sdev, self%dist_inpl%minv, self%dist_inpl%maxv
        write(logfhandle,604) '>>> SHIFT INCR ARG         AVG/SDEV/MIN/MAX:', self%shincarg%avg, self%shincarg%sdev, self%shincarg%minv, self%shincarg%maxv
        write(logfhandle,604) '>>> % SEARCH SPACE SCANNED AVG/SDEV/MIN/MAX:', self%frac_srch%avg, self%frac_srch%sdev, self%frac_srch%minv, self%frac_srch%maxv
        ! correlation & particle weights
        write(logfhandle,604) '>>> CORRELATION            AVG/SDEV/MIN/MAX:', self%corr%avg, self%corr%sdev, self%corr%minv, self%corr%maxv
        write(logfhandle,601) '>>> % PARTICLES   CC > CC_AVG - 2 * CC_SDEV:', 100. * real(count(corrs > corr_t .and. mask)) / real(count(mask))
        write(logfhandle,601) '>>> CORRELATION                   THRESHOLD:', corr_t
        write(logfhandle,604) '>>> PARTICLE WEIGHT        AVG/SDEV/MIN/MAX:', self%pw%avg, self%pw%sdev, self%pw%minv, self%pw%maxv
        write(logfhandle,601) '>>> % PARTICLES WITH NONZERO WEIGHT         ', percen_nonzero_pw
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
        converged = .false.
        chk4conv  = .true.
        if( cline%defined('converge') )then
            if( cline%get_carg('converge') .eq. 'no' )then
                ! never converge
                chk4conv = .false.
            else
                ! to indicate that we need to check for convergence
                chk4conv = .true.
            endif
        endif
        if( chk4conv )then
            ! determine convergence
            if( ncls > 1 )then
                converged = .false.
                ! set limits for convergence
                if( (params_glob%l_frac_update) .or. (params_glob%stream.eq.'yes') )then
                    overlap_lim  = OVERLAP_2D_FRAC
                    fracsrch_lim = FRACSRCHSPACE_FRAC
                else if( trim(params_glob%tseries) .eq. 'yes' )then
                    overlap_lim  = OVERLAP_2D_NANO
                else
                    overlap_lim  = OVERLAP_2D
                    fracsrch_lim = FRACSRCHSPACE_2D
                endif
                ! override if present on command line
                if( cline%defined('overlap')  ) overlap_lim  = cline%get_rarg('overlap')
                if( cline%defined('fracsrch') ) fracsrch_lim = cline%get_rarg('fracsrch')
                ! test for convergence
                if( (params_glob%l_frac_update) .or. (params_glob%stream.eq.'yes') )then
                    converged = ( self%mi_class > overlap_lim .and. self%frac_srch%avg > fracsrch_lim )
                    self%progress = progress_estimate_2D(real(params_glob%which_iter), self%mi_class, overlap_lim, self%frac_srch%avg, fracsrch_lim, 0.0, 0.0)
                else if( trim(params_glob%tseries) .eq. 'yes' )then
                    converged = self%mi_class > overlap_lim
                    self%progress = progress_estimate_2D(real(params_glob%which_iter), self%mi_class, overlap_lim, 0.0, 0.0, 0.0, 0.0)
                else
                    converged = ( self%mi_class > overlap_lim .and. self%frac_srch%avg > fracsrch_lim )
                    self%progress = progress_estimate_2D(real(params_glob%which_iter), self%mi_class, overlap_lim, self%frac_srch%avg, fracsrch_lim, 0.0, 0.0)
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
        endif
        ! stats
        call self%ostats%new(1, is_ptcl=.false.)
        call self%ostats%set(1,'ITERATION',real(params_glob%which_iter))
        call self%ostats%set(1,'CLASS_OVERLAP',self%mi_class)
        call self%ostats%set(1,'PARTICLE_UPDATES',avg_updatecnt)
        call self%ostats%set(1,'IN-PLANE_DIST',self%dist_inpl%avg)
        call self%ostats%set(1,'SEARCH_SPACE_SCANNED',self%frac_srch%avg)
        call self%ostats%set(1,'CORRELATION',self%corr%avg)
        call self%ostats%write(STATS_FILE)
        ! destruct
        deallocate(mask, updatecnts, states, corrs, pws)
        call self%ostats%kill
    end function check_conv2D

    function check_conv3D( self, cline, msk ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        real,               intent(in)    :: msk
        real,    allocatable :: state_mi_joint(:), statepops(:), updatecnts(:), pws(:), states(:), corrs(:)
        logical, allocatable :: mask(:)
        real    :: min_state_mi_joint, avg_updatecnt, percen_nonzero_pw, overlap_lim, fracsrch_lim, corr_t
        logical :: converged
        integer :: iptcl, istate
        601 format(A,1X,F12.3)
        604 format(A,1X,F12.3,1X,F12.3,1X,F12.3,1X,F12.3)
        states        = build_glob%spproj_field%get_all('state')
        corrs         = build_glob%spproj_field%get_all('corr')

        updatecnts    = build_glob%spproj_field%get_all('updatecnt')
        avg_updatecnt = sum(updatecnts) / real(count(states > 0.5))
        allocate(mask(size(updatecnts)), source=updatecnts > 0.5 .and. states > 0.5)
        pws = build_glob%spproj_field%get_all('w')
        percen_nonzero_pw = (real(count(mask .and. (pws > TINY))) / real(count(mask))) * 100.
        call build_glob%spproj_field%stats('corr',       self%corr,       mask=mask)
        call build_glob%spproj_field%stats('npeaks',     self%npeaks,     mask=mask)
        if( self%npeaks%avg > 1e-6 )then
        call build_glob%spproj_field%stats('cc_peak',    self%cc_peak,    mask=mask, nozero=.true.)
        call build_glob%spproj_field%stats('cc_nonpeak', self%cc_nonpeak, mask=mask, nozero=.true.)
        call build_glob%spproj_field%stats('dist_peaks', self%dist_peaks, mask=mask, nozero=.true.)
        endif
        call build_glob%spproj_field%stats('dist',       self%dist,       mask=mask)
        call build_glob%spproj_field%stats('dist_inpl',  self%dist_inpl,  mask=mask)
        call build_glob%spproj_field%stats('frac',       self%frac_srch,  mask=mask)
        call build_glob%spproj_field%stats('w',          self%pw,         mask=mask)
        call build_glob%spproj_field%stats('shincarg',   self%shincarg,   mask=mask)
        self%mi_proj   = build_glob%spproj_field%get_avg('mi_proj',   mask=mask)
        self%mi_state  = build_glob%spproj_field%get_avg('mi_state',  mask=mask)
        corr_t         = self%corr%avg - 2. * self%corr%sdev
        ! overlaps and particle updates
        write(logfhandle,601) '>>> ORIENTATION OVERLAP:                      ', self%mi_proj
        if( params_glob%nstates > 1 )then
        write(logfhandle,601) '>>> STATE OVERLAP:                            ', self%mi_state
        endif
        write(logfhandle,601) '>>> # PARTICLE UPDATES       AVG:             ', avg_updatecnt
        ! dists and % search space
        write(logfhandle,604) '>>> DIST BTW BEST ORIS (DEG) AVG/SDEV/MIN/MAX:', self%dist%avg, self%dist%sdev, self%dist%minv, self%dist%maxv
        write(logfhandle,604) '>>> IN-PLANE DIST      (DEG) AVG/SDEV/MIN/MAX:', self%dist_inpl%avg, self%dist_inpl%sdev, self%dist_inpl%minv, self%dist_inpl%maxv
        if( self%npeaks%avg > 1e-6 )then
        write(logfhandle,604) '>>> # PROJECTION PEAKS       AVG/SDEV/MIN/MAX:', self%npeaks%avg, self%npeaks%sdev, self%npeaks%minv, self%npeaks%maxv
        write(logfhandle,604) '>>> PEAK DIST          (DEG) AVG/SDEV/MIN/MAX:', self%dist_peaks%avg, self%dist_peaks%sdev, self%dist_peaks%minv, self%dist_peaks%maxv
        endif
        write(logfhandle,604) '>>> SHIFT INCR ARG           AVG/SDEV/MIN/MAX:', self%shincarg%avg, self%shincarg%sdev, self%shincarg%minv, self%shincarg%maxv
        write(logfhandle,604) '>>> % SEARCH SPACE SCANNED   AVG/SDEV/MIN/MAX:', self%frac_srch%avg, self%frac_srch%sdev, self%frac_srch%minv, self%frac_srch%maxv
        ! correlation & particle weights
        write(logfhandle,604) '>>> CORRELATION              AVG/SDEV/MIN/MAX:', self%corr%avg, self%corr%sdev, self%corr%minv, self%corr%maxv
        if( self%npeaks%avg > 1e-6 )then
        write(logfhandle,604) '>>> CORRELATION, PEAK        AVG/SDEV/MIN/MAX:', self%cc_peak%avg, self%cc_peak%sdev, self%cc_peak%minv, self%cc_peak%maxv
        write(logfhandle,604) '>>> CORRELATION, NONPEAK     AVG/SDEV/MIN/MAX:', self%cc_nonpeak%avg, self%cc_nonpeak%sdev, self%cc_nonpeak%minv, self%cc_nonpeak%maxv
        endif
        write(logfhandle,601) '>>> % PARTICLES     CC > CC_AVG - 2 * CC_SDEV:', 100. * real(count(corrs > corr_t .and. mask)) / real(count(mask))
        write(logfhandle,601) '>>> CORRELATION                     THRESHOLD:', corr_t
        write(logfhandle,604) '>>> PARTICLE WEIGHT          AVG/SDEV/MIN/MAX:', self%pw%avg, self%pw%sdev, self%pw%minv, self%pw%maxv
        write(logfhandle,601) '>>> % PARTICLES WITH NONZERO WEIGHT           ', percen_nonzero_pw
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
        ! set limits for convergence
        overlap_lim  = OVERLAP_3D
        fracsrch_lim = FRACSRCHSPACE_3D
        ! override if present on command line
        if( cline%defined('overlap')  ) overlap_lim  = cline%get_rarg('overlap')
        if( cline%defined('fracsrch') ) fracsrch_lim = cline%get_rarg('fracsrch')
        ! determine convergence
        if( params_glob%nstates == 1 )then
            if( self%frac_srch%avg > fracsrch_lim .and. self%mi_proj  > overlap_lim )then
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
            if( min_state_mi_joint > OVERLAP_STATE_JOINT .and.&
                self%mi_state      > OVERLAP_STATE       .and.&
                self%frac_srch%avg > fracsrch_lim        )then
                write(logfhandle,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(logfhandle,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
            deallocate( state_mi_joint, statepops )
        endif
        ! stats
        call self%ostats%new(1, is_ptcl=.false.)
        call self%ostats%set(1,'ORIENTATION_OVERLAP',self%mi_proj)
        if( params_glob%nstates > 1 ) call self%ostats%set(1,'STATE_OVERLAP', self%mi_state)
        call self%ostats%set(1,'PARTICLE_UPDATES',avg_updatecnt)
        call self%ostats%set(1,'DIST_BTW_BEST_ORIS',self%dist%avg)
        call self%ostats%set(1,'IN-PLANE_DIST',self%dist_inpl%avg)
        call self%ostats%set(1,'PARTICLE_WEIGHT',self%pw%avg)
        call self%ostats%set(1,'SEARCH_SPACE_SCANNED',self%frac_srch%avg)
        call self%ostats%set(1,'CORRELATION',self%corr%avg)
        call self%ostats%set(1,'SHIFT_INCR_ARG',self%shincarg%avg)
        call self%ostats%write(STATS_FILE)
        ! destruct
        deallocate(mask, updatecnts, pws, states, corrs)
        call self%ostats%kill
    end function check_conv3D

    function check_conv3Dc( self, cline, msk ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        real,               intent(in)    :: msk
        real,    allocatable :: updatecnts(:), pws(:), states(:), corrs(:)
        logical, allocatable :: mask(:)
        real    :: avg_updatecnt, percen_nonzero_pw, overlap_lim, fracsrch_lim, corr_t
        logical :: converged
        integer :: iptcl
        601 format(A,1X,F12.3)
        604 format(A,1X,F12.3,1X,F12.3,1X,F12.3,1X,F12.3)
        states        = build_glob%spproj_field%get_all('state')
        corrs         = build_glob%spproj_field%get_all('corr')
        updatecnts    = build_glob%spproj_field%get_all('updatecnt')
        avg_updatecnt = sum(updatecnts) / real(count(states > 0.5))
        allocate(mask(size(updatecnts)), source=updatecnts > 0.5 .and. states > 0.5)
        pws = build_glob%spproj_field%get_all('w')
        percen_nonzero_pw = (real(count(mask .and. (pws > TINY))) / real(count(mask))) * 100.
        call build_glob%spproj_field%stats('corr',       self%corr,       mask=mask)
        call build_glob%spproj_field%stats('npeaks',     self%npeaks,     mask=mask)
        if( self%npeaks%avg > 1e-6 )then
        call build_glob%spproj_field%stats('cc_peak',    self%cc_peak,    mask=mask, nozero=.true.)
        call build_glob%spproj_field%stats('cc_nonpeak', self%cc_nonpeak, mask=mask, nozero=.true.)
        call build_glob%spproj_field%stats('dist_peaks', self%dist_peaks, mask=mask, nozero=.true.)
        endif
        call build_glob%spproj_field%stats('dist',       self%dist,       mask=mask)
        call build_glob%spproj_field%stats('dist_inpl',  self%dist_inpl,  mask=mask)
        call build_glob%spproj_field%stats('frac',       self%frac_srch,  mask=mask)
        call build_glob%spproj_field%stats('frac_sh',    self%frac_sh,    mask=mask)
        call build_glob%spproj_field%stats('w',          self%pw,         mask=mask)
        call build_glob%spproj_field%stats('shincarg',   self%shincarg,   mask=mask)
        call build_glob%spproj_field%stats('nevals',     self%nevals,     mask=mask, nozero=.true.)
        call build_glob%spproj_field%stats('ngevals',    self%ngevals,    mask=mask, nozero=.true.)
        call build_glob%spproj_field%stats('better',     self%better,     mask=mask)
        call build_glob%spproj_field%stats('better_l',   self%better_l,   mask=mask)
        corr_t = self%corr%avg - 2. * self%corr%sdev
        ! particle updates
        write(logfhandle,601) '>>> # PARTICLE UPDATES       AVG:             ', avg_updatecnt
        ! dists and % search space
        write(logfhandle,604) '>>> DIST BTW BEST ORIS (DEG) AVG/SDEV/MIN/MAX:', self%dist%avg, self%dist%sdev, self%dist%minv, self%dist%maxv
        write(logfhandle,604) '>>> IN-PLANE DIST      (DEG) AVG/SDEV/MIN/MAX:', self%dist_inpl%avg, self%dist_inpl%sdev, self%dist_inpl%minv, self%dist_inpl%maxv
        if( self%npeaks%avg > 1e-6 )then
        write(logfhandle,604) '>>> # PROJECTION PEAKS       AVG/SDEV/MIN/MAX:', self%npeaks%avg, self%npeaks%sdev, self%npeaks%minv, self%npeaks%maxv
        write(logfhandle,604) '>>> PEAK DIST          (DEG) AVG/SDEV/MIN/MAX:', self%dist_peaks%avg, self%dist_peaks%sdev, self%dist_peaks%minv, self%dist_peaks%maxv
        endif
        write(logfhandle,604) '>>> SHIFT INCR ARG           AVG/SDEV/MIN/MAX:', self%shincarg%avg, self%shincarg%sdev, self%shincarg%minv, self%shincarg%maxv
        write(logfhandle,604) '>>> % EULER SPACE SCANNED    AVG/SDEV/MIN/MAX:', self%frac_srch%avg, self%frac_srch%sdev, self%frac_srch%minv, self%frac_srch%maxv
        write(logfhandle,604) '>>> % SHIFT SPACE SCANNED    AVG/SDEV/MIN/MAX:', self%frac_sh%avg, self%frac_sh%sdev, self%frac_sh%minv, self%frac_sh%maxv
        write(logfhandle,604) '>>> % IMPROVED SOLUTIONS     AVG/SDEV/MIN/MAX:', 100.*self%better%avg, 100.*self%better%sdev, 100.*self%better%minv, 100.*self%better%maxv
        write(logfhandle,604) '>>> % IMPROVED LBFGS-B       AVG/SDEV/MIN/MAX:', 100.*self%better_l%avg, 100.*self%better_l%sdev, 100.*self%better_l%minv, 100.*self%better_l%maxv
        ! correlation & particle weights
        write(logfhandle,604) '>>> CORRELATION              AVG/SDEV/MIN/MAX:', self%corr%avg, self%corr%sdev, self%corr%minv, self%corr%maxv
        if( self%npeaks%avg > 1e-6 )then
        write(logfhandle,604) '>>> CORRELATION, PEAK        AVG/SDEV/MIN/MAX:', self%cc_peak%avg, self%cc_peak%sdev, self%cc_peak%minv, self%cc_peak%maxv
        write(logfhandle,604) '>>> CORRELATION, NONPEAK     AVG/SDEV/MIN/MAX:', self%cc_nonpeak%avg, self%cc_nonpeak%sdev, self%cc_nonpeak%minv, self%cc_nonpeak%maxv
        endif
        write(logfhandle,601) '>>> % PARTICLES     CC > CC_AVG - 2 * CC_SDEV:', 100. * real(count(corrs > corr_t .and. mask)) / real(count(mask))
        write(logfhandle,601) '>>> CORRELATION                     THRESHOLD:', corr_t
        write(logfhandle,604) '>>> PARTICLE WEIGHT          AVG/SDEV/MIN/MAX:', self%pw%avg, self%pw%sdev, self%pw%minv, self%pw%maxv
        write(logfhandle,601) '>>> % PARTICLES WITH NONZERO WEIGHT           ', percen_nonzero_pw
        ! cost and gradient evals
        write(logfhandle,604) '>>> # COST FUN EVALS         AVG/SDEV/MIN/MAX:', self%nevals%avg, self%nevals%sdev, self%nevals%minv, self%nevals%maxv
        write(logfhandle,604) '>>> # GRADIENT EVALS         AVG/SDEV/MIN/MAX:', self%ngevals%avg, self%ngevals%sdev, self%ngevals%minv, self%ngevals%maxv

        ! determine convergence
        converged = .false.
        if( self%better%avg * 100. < 5. )then ! less than 5 % improved solutions
            converged = .true.
        endif
        ! stats
        call self%ostats%new(1, is_ptcl=.false.)
        call self%ostats%set(1,'PARTICLE_UPDATES',avg_updatecnt)
        call self%ostats%set(1,'DIST_BTW_BEST_ORIS',self%dist%avg)
        call self%ostats%set(1,'IN-PLANE_DIST',self%dist_inpl%avg)
        call self%ostats%set(1,'PARTICLE_WEIGHT',self%pw%avg)
        call self%ostats%set(1,'CORRELATION',self%corr%avg)
        call self%ostats%set(1,'SHIFT_INCR_ARG',self%shincarg%avg)
        call self%ostats%write(STATS_FILE)
        ! destruct
        deallocate(mask, updatecnts, pws, corrs)
        call self%ostats%kill
    end function check_conv3Dc

    function check_conv_cluster( self, cline ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        integer, allocatable :: statepops(:)
        logical :: converged
        integer :: istate

        converged = .false.

        ! 601 format(A,1X,F8.3)
        ! 604 format(A,1X,F8.3,1X,F8.3,1X,F8.3,1X,F8.3)
        ! call build_glob%spproj_field%stats('frac', self%frac_srch)
        ! self%mi_state  = build_glob%spproj_field%get_avg('mi_state')
        ! write(logfhandle,601) '>>> STATE OVERLAP:                          ', self%mi_state
        ! write(logfhandle,604) '>>> % SEARCH SPACE SCANNED AVG/SDEV/MIN/MAX:',&
        ! &self%frac_srch%avg, self%frac_srch%sdev, self%frac_srch%minv, self%frac_srch%maxv
        ! ! provides convergence stats for multiple states
        ! ! by calculating mi_joint for individual states
        ! call build_glob%spproj_field%get_pops(statepops,'state')
        ! call build_glob%spproj_field%stats('corr', self%corr)
        ! write(logfhandle,604) '>>> CORRELATION            AVG/SDEV/MIN/MAX:',&
        ! &self%corr%avg, self%corr%sdev, self%corr%minv, self%corr%maxv
        ! ! print the overlaps and pops for the different states
        ! do istate=1,params_glob%nstates
        !     write(logfhandle,'(A,I2,1X,A,1X,I8)') '>>> STATE ',istate,'POPULATION:', statepops(istate)
        ! end do
        ! if( self%mi_state > OVERLAP_STATE_HET .and.&
        !     self%frac_srch%avg > FRACSRCHSPACE_HET     )then
        !     write(logfhandle,'(A)') '>>> CONVERGED: .YES.'
        !     converged = .true.
        ! else
        !     write(logfhandle,'(A)') '>>> CONVERGED: .NO.'
        !     converged = .false.
        ! endif
        ! ! stats
        ! call self%ostats%new(1, is_ptcl=.false.)
        ! call self%ostats%set(1,'STATE_OVERLAP',       self%mi_state)
        ! call self%ostats%set(1,'SEARCH_SPACE_SCANNED',self%frac_srch%avg)
        ! call self%ostats%set(1,'CORRELATION',         self%corr%avg)
        ! do istate=1,params_glob%nstates
        !     call self%ostats%set(1,'STATE_POPULATION_'//int2str(istate), real(statepops(istate)))
        ! enddo
        ! call self%ostats%write(STATS_FILE)
        ! ! cleanup
        ! call self%ostats%kill
        ! deallocate( statepops )
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
            case('progress')
                get = self%progress
        end select
    end function get

end module simple_convergence
