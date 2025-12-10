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
    type(stats_struct) :: score      !< objective function stats
    type(stats_struct) :: dist       !< angular distance stats
    type(stats_struct) :: dist_inpl  !< in-plane angular distance stats
    type(stats_struct) :: frac_srch  !< fraction of search space scanned stats
    type(stats_struct) :: shincarg   !< shift increment
    type(stats_struct) :: pw         !< particle weights
    type(stats_struct) :: lp         !< low-pass limit
    type(stats_struct) :: lp_est     !< low-pass limit, estimated
    type(stats_struct) :: res        !< resolution @ FSC=0.143
    integer :: iteration   = 0       !< current interation
    real    :: mi_class    = 0.      !< class parameter distribution overlap
    real    :: mi_proj     = 0.      !< projection parameter distribution overlap
    real    :: mi_state    = 0.      !< state parameter distribution overlap
    real    :: frac_greedy = 0.      !< fraction of greedy searches
    real    :: progress    = 0.      !< progress estimation
  contains
    procedure :: read
    procedure :: check_conv2D
    procedure :: check_conv3D
    procedure :: append_stats
    procedure :: plot_projdirs
    procedure :: get
end type convergence

contains

    subroutine read( self, l_err )
        class(convergence), intent(inout) :: self
        logical,            intent(out)   :: l_err
        type(oris) :: ostats
        l_err = .false.
        call ostats%new(1, is_ptcl=.false.)
        if( file_exists(STATS_FILE) )then
            call ostats%read(string(STATS_FILE))
            self%iteration     = nint(ostats%get(1,'ITERATION'))
            self%score%avg     = ostats%get(1,'SCORE')
            self%frac_srch%avg = ostats%get(1,'SEARCH_SPACE_SCANNED')
            self%mi_proj       = ostats%get(1,'ORIENTATION_OVERLAP')
            self%mi_state      = nint(ostats%get(1,'STATE_OVERLAP'))
            self%dist%avg      = ostats%get(1,'DIST_BTW_BEST_ORIS')
            self%dist_inpl%avg = ostats%get(1,'IN-PLANE_DIST')
            self%shincarg%avg  = ostats%get(1,'SHIFT_INCR_ARG')
        else
            l_err = .true.
        endif
        call ostats%kill
    end subroutine read

    function check_conv2D( self, cline, os, ncls, msk ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(oris),        intent(inout) :: os
        integer,            intent(in)    :: ncls
        real,               intent(in)    :: msk
        type(oris)           :: ostats
        real,    allocatable :: updatecnts(:), states(:), scores(:), sampled(:)
        logical, allocatable :: mask(:)
        integer :: n, nptcls
        real    :: overlap_lim, fracsrch_lim
        real    :: percen_sampled, percen_updated, percen_avg, sampled_lb
        logical :: converged, chk4conv
        601 format(A,1X,F12.3)
        604 format(A,1X,F12.3,1X,F12.3,1X,F12.3,1X,F12.3)
        states         = os%get_all('state')
        scores         = os%get_all('corr')
        updatecnts     = os%get_all('updatecnt')
        sampled        = os%get_all('sampled')
        n              = size(states)
        nptcls         = count(states > 0.5)
        sampled_lb     = maxval(sampled) - 0.5
        percen_sampled = (real(count(sampled    > sampled_lb .and. states > 0.5)) / real(nptcls)) * 100.
        percen_updated = (real(count(updatecnts > 0.5        .and. states > 0.5)) / real(nptcls)) * 100.
        percen_avg     = percen_sampled
        if( params_glob%l_update_frac )then
            allocate(mask(n), source=sampled    > sampled_lb .and. states > 0.5)
        else
            allocate(mask(n), source=updatecnts > 0.5 .and. states > 0.5)
        endif
        call os%stats('corr',      self%score,     mask=mask)
        call os%stats('dist_inpl', self%dist_inpl, mask=mask)
        call os%stats('frac',      self%frac_srch, mask=mask)
        call os%stats('shincarg',  self%shincarg,  mask=mask)
        call os%stats('lp',        self%lp,        mask=mask)
        call os%stats('lp_est',    self%lp_est,    mask=mask)
        call os%stats('res',       self%res,       mask=mask)
        self%mi_class = os%get_avg('mi_class',     mask=mask)
        ! overlaps and particle updates
        write(logfhandle,601) '>>> CLASS OVERLAP:                            ', self%mi_class
        write(logfhandle,601) '>>> % PARTICLES SAMPLED THIS ITERATION        ', percen_sampled
        write(logfhandle,601) '>>> % PARTICLES UPDATED SO FAR                ', percen_updated
        write(logfhandle,601) '>>> % PARTICLES USED FOR AVERAGING            ', percen_avg
        ! dists and % search space
        write(logfhandle,604) '>>> IN-PLANE DIST    (DEG)   AVG/SDEV/MIN/MAX:', self%dist_inpl%avg, self%dist_inpl%sdev, self%dist_inpl%minv, self%dist_inpl%maxv
        write(logfhandle,604) '>>> SHIFT INCR ARG           AVG/SDEV/MIN/MAX:', self%shincarg%avg,  self%shincarg%sdev,  self%shincarg%minv,  self%shincarg%maxv
        write(logfhandle,604) '>>> % SEARCH SPACE SCANNED   AVG/SDEV/MIN/MAX:', self%frac_srch%avg, self%frac_srch%sdev, self%frac_srch%minv, self%frac_srch%maxv
        write(logfhandle,604) '>>> MATCHING  LOW-PASS LIMIT AVG/SDEV/MIN/MAX:', self%lp%avg,        self%lp%sdev,        self%lp%minv,        self%lp%maxv
        write(logfhandle,604) '>>> ESTIMATED LOW-PASS LIMIT AVG/SDEV/MIN/MAX:', self%lp_est%avg,    self%lp_est%sdev,    self%lp_est%minv,    self%lp_est%maxv
        write(logfhandle,604) '>>> RESOLUTION @ FSC=0.143   AVG/SDEV/MIN/MAX:', self%res%avg,       self%res%sdev,       self%res%minv,       self%res%maxv
        ! score
        write(logfhandle,604) '>>> SCORE [0,1]              AVG/SDEV/MIN/MAX:', self%score%avg, self%score%sdev, self%score%minv, self%score%maxv
        ! dynamic shift search range update
        if( self%frac_srch%avg >= FRAC_SH_LIM )then
            if( cline%defined('trs') .or. params_glob%trs >  MINSHIFT)then
                ! no change: bound was explicitely defined or is large enough
            else
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
                if( (params_glob%l_update_frac) .or. (params_glob%stream.eq.'yes') )then
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
                if( trim(params_glob%refine).eq.'inpl' )then
                    converged = self%dist_inpl%avg < 0.5
                else if( (params_glob%l_update_frac) .or. (params_glob%stream.eq.'yes') )then
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
        call ostats%new(1, is_ptcl=.false.)
        call ostats%set(1,'ITERATION',                params_glob%which_iter)
        call ostats%set(1,'CLASS_OVERLAP',            self%mi_class)
        call ostats%set(1,'PERCEN_PARTICLES_SAMPLED', percen_sampled)
        call ostats%set(1,'PERCEN_PARTICLES_UPDATED', percen_updated)
        call ostats%set(1,'PERCEN_PARTICLES_AVERAGED',percen_avg)
        call ostats%set(1,'IN-PLANE_DIST',            self%dist_inpl%avg)
        call ostats%set(1,'SEARCH_SPACE_SCANNED',     self%frac_srch%avg)
        call ostats%set(1,'SCORE',                    self%score%avg)
        call ostats%write(string(STATS_FILE))
        ! destruct
        deallocate(mask, updatecnts, states, scores, sampled)
        call ostats%kill
    end function check_conv2D

    function check_conv3D( self, cline, msk ) result( converged )
        class(convergence), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        real,               intent(in)    :: msk
        type(oris)           :: ostats
        real,    allocatable :: state_mi_joint(:), statepops(:), updatecnts(:), states(:), scores(:), sampled(:)
        logical, allocatable :: mask(:)
        real    :: min_state_mi_joint, overlap_lim, fracsrch_lim, trail_rec_ufrac
        real    :: percen_sampled, percen_updated, percen_avg, sampled_lb
        logical :: converged
        integer :: iptcl, istate, n, nptcls, ucnt
        601 format(A,1X,F12.3)
        604 format(A,1X,F12.3,1X,F12.3,1X,F12.3,1X,F12.3)
        607 format(A,1X,F4.2)
        609 format(A)
        states         = build_glob%spproj_field%get_all('state')
        scores         = build_glob%spproj_field%get_all('corr')
        updatecnts     = build_glob%spproj_field%get_all('updatecnt')
        sampled        = build_glob%spproj_field%get_all('sampled')
        n              = size(states)
        nptcls         = count(states > 0.5)
        sampled_lb     = maxval(sampled) - 0.5
        percen_sampled = (real(count(sampled    > sampled_lb .and. states > 0.5)) / real(nptcls)) * 100.
        percen_updated = (real(count(updatecnts > 0.5        .and. states > 0.5)) / real(nptcls)) * 100.
        percen_avg     = percen_sampled
        if( params_glob%l_update_frac )then
            allocate(mask(n), source=sampled    > sampled_lb .and. states > 0.5)
        else
            allocate(mask(n), source=updatecnts > 0.5 .and. states > 0.5)
        endif
        call build_glob%spproj_field%stats('corr',       self%score,      mask=mask)
        call build_glob%spproj_field%stats('dist',       self%dist,       mask=mask)
        call build_glob%spproj_field%stats('dist_inpl',  self%dist_inpl,  mask=mask)
        call build_glob%spproj_field%stats('frac',       self%frac_srch,  mask=mask)
        call build_glob%spproj_field%stats('shincarg',   self%shincarg,   mask=mask)
        call build_glob%spproj_field%stats('w',          self%pw,         mask=mask)
        call build_glob%spproj_field%stats('lp',         self%lp,         mask=mask)
        call build_glob%spproj_field%stats('lp_est',     self%lp_est,     mask=mask)
        call build_glob%spproj_field%stats('res',        self%res,        mask=mask)
        self%mi_proj     = build_glob%spproj_field%get_avg('mi_proj',     mask=mask)
        self%mi_state    = build_glob%spproj_field%get_avg('mi_state',    mask=mask)
        self%frac_greedy = build_glob%spproj_field%get_avg('frac_greedy', mask=mask)
        ! overlaps and particle updates
        write(logfhandle,601) '>>> ORIENTATION OVERLAP:                      ', self%mi_proj
        if( params_glob%nstates > 1 )then
        write(logfhandle,601) '>>> STATE OVERLAP:                            ', self%mi_state
        endif
        write(logfhandle,601) '>>> % PARTICLES SAMPLED THIS ITERATION        ', percen_sampled
        write(logfhandle,601) '>>> % PARTICLES UPDATED SO FAR                ', percen_updated
        write(logfhandle,601) '>>> % PARTICLES USED FOR UPDATING THE MODEL   ', percen_avg
        write(logfhandle,601) '>>> % GREEDY SEARCHES                         ', self%frac_greedy * 100.
        ! dists and % search space
        write(logfhandle,604) '>>> DIST BTW BEST ORIS (DEG) AVG/SDEV/MIN/MAX:', self%dist%avg,      self%dist%sdev,      self%dist%minv,      self%dist%maxv
        write(logfhandle,604) '>>> IN-PLANE DIST      (DEG) AVG/SDEV/MIN/MAX:', self%dist_inpl%avg, self%dist_inpl%sdev, self%dist_inpl%minv, self%dist_inpl%maxv
        write(logfhandle,604) '>>> SHIFT INCR ARG           AVG/SDEV/MIN/MAX:', self%shincarg%avg,  self%shincarg%sdev,  self%shincarg%minv,  self%shincarg%maxv
        write(logfhandle,604) '>>> % SEARCH SPACE SCANNED   AVG/SDEV/MIN/MAX:', self%frac_srch%avg, self%frac_srch%sdev, self%frac_srch%minv, self%frac_srch%maxv
        write(logfhandle,604) '>>> PARTICLE WEIGHTS         AVG/SDEV/MIN/MAX:', self%pw%avg,        self%pw%sdev,        self%pw%minv,        self%pw%maxv
        write(logfhandle,604) '>>> MATCHING  LOW-PASS LIMIT AVG/SDEV/MIN/MAX:', self%lp%avg,        self%lp%sdev,        self%lp%minv,        self%lp%maxv
        write(logfhandle,604) '>>> ESTIMATED LOW-PASS LIMIT AVG/SDEV/MIN/MAX:', self%lp_est%avg,    self%lp_est%sdev,    self%lp_est%minv,    self%lp_est%maxv
        write(logfhandle,604) '>>> RESOLUTION @ FSC=0.143   AVG/SDEV/MIN/MAX:', self%res%avg,       self%res%sdev,       self%res%minv,       self%res%maxv
        ! score
        write(logfhandle,604) '>>> SCORE [0,1]              AVG/SDEV/MIN/MAX:', self%score%avg, self%score%sdev, self%score%minv, self%score%maxv
        write(logfhandle,609) '>>> REFINEMENT MODE IS '//trim(params_glob%refine)
        if( trim(params_glob%gauref).eq.'yes' )then
        write(logfhandle,607) '>>> GAU REGULARIZATION IS ON'
        else
        write(logfhandle,609) '>>> GAU REGULARIZATION IS OFF'
        endif
        if( params_glob%l_ml_reg )then
        write(logfhandle,607) '>>> ML  REGULARIZATION IS ON, TAU:    ', params_glob%tau
        else
        write(logfhandle,609) '>>> ML  REGULARIZATION IS OFF'
        endif
        if( params_glob%l_icm )then
        write(logfhandle,607) '>>> ICM REGULARIZATION IS ON, LAMBDA: ', params_glob%lambda
        else
        write(logfhandle,609) '>>> ICM REGULARIZATION IS OFF'
        endif
        if( params_glob%l_update_frac )then
        if( cline%defined('ufrac_trec') )then
        write(logfhandle,607) '>>> TRAILING REC UPDATE FRACTION:     ', params_glob%ufrac_trec
        else
        trail_rec_ufrac = real(count(mask)) / real(count(updatecnts > 0.5 .and. states > 0.5))
        write(logfhandle,607) '>>> TRAILING REC UPDATE FRACTION:     ', trail_rec_ufrac
        endif
        if( params_glob%l_fillin .and. mod(params_glob%which_iter,5) == 0 )then
        write(logfhandle,609) '>>> FILLIN PARTICLE SAMPLING WAS ON'
        else
        write(logfhandle,609) '>>> FILLIN PARTICLE SAMPLING WAS OFF'
        endif
        endif
        ! dynamic shift search range update
        if( self%frac_srch%avg >= FRAC_SH_LIM )then
            if( .not. cline%defined('trs') .or. &
                & params_glob%trs <  MINSHIFT )then
                if( cline%defined('trs') .and. params_glob%trs<0.01 )then
                    ! bound was defined, no update
                else
                    ! determine shift bounds
                    params_glob%trs = MSK_FRAC*msk
                    params_glob%trs = max(MINSHIFT,params_glob%trs)
                    params_glob%trs = min(MAXSHIFT,params_glob%trs)
                    ! set shift search flag
                    params_glob%l_doshift = .true.
                endif
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
                ucnt   = build_glob%spproj_field%get_updatecnt(iptcl)
                istate = build_glob%spproj_field%get_state(iptcl)
                if( istate==0 .or. ucnt==0 ) cycle
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
            if( min_state_mi_joint > OVERLAP_STATE_JOINT .and. self%frac_srch%avg > fracsrch_lim )then
                write(logfhandle,'(A)') '>>> CONVERGED: .YES.'
                converged = .true.
            else
                write(logfhandle,'(A)') '>>> CONVERGED: .NO.'
                converged = .false.
            endif
            deallocate( state_mi_joint, statepops )
        endif
        ! stats
        call ostats%new(1, is_ptcl=.false.)
        call ostats%set(1,'ITERATION',                   params_glob%which_iter)
        call ostats%set(1,'ORIENTATION_OVERLAP',         self%mi_proj)
        if( params_glob%nstates > 1 ) call ostats%set(1,'STATE_OVERLAP', self%mi_state)
        call ostats%set(1,'PERCEN_PARTICLES_SAMPLED',    percen_sampled)
        call ostats%set(1,'PERCEN_PARTICLES_UPDATED',    percen_updated)
        call ostats%set(1,'PERCEN_PARTICLES_AVERAGED',   percen_avg)
        call ostats%set(1,'PERCEN_GREEDY_SEARCHES',      self%frac_greedy * 100.)
        call ostats%set(1,'DIST_BTW_BEST_ORIS',          self%dist%avg)
        call ostats%set(1,'IN-PLANE_DIST',               self%dist_inpl%avg)
        call ostats%set(1,'SHIFT_INCR_ARG',              self%shincarg%avg)
        call ostats%set(1,'PERCEN_SEARCH_SPACE_SCANNED', self%frac_srch%avg)
        call ostats%set(1,'LP_MATCHING',                 self%lp%avg)
        call ostats%set(1,'LP_ESTIMATED',                self%lp_est%avg)
        call ostats%set(1,'RESOLUTION',                  self%res%avg)
        call ostats%set(1,'SCORE',                       self%score%avg)
        if( params_glob%l_ml_reg )then
            call ostats%set(1,'ML_REGULARIZATION',                        1.0)
            call ostats%set(1,'ML_REGULARIZATION_TAU',        params_glob%tau)
        else
            call ostats%set(1,'ML_REGULARIZATION',                        0.0)
        endif
        if( params_glob%l_icm )then
            call ostats%set(1,'ICM_REGULARIZATION',                       1.0)
            call ostats%set(1,'ICM_REGULARIZATION_LAMBDA', params_glob%lambda)
        endif
        if( params_glob%l_update_frac )then
            call ostats%set(1,'TRAIL_REC_UPDATE_FRACTION',    trail_rec_ufrac)
        endif
        
        call ostats%write(string(STATS_FILE))
        call self%append_stats(ostats)
        call self%plot_projdirs(mask)
        ! destruct
        if( allocated(state_mi_joint) ) deallocate(state_mi_joint)
        if( allocated(statepops)      ) deallocate(statepops)
        deallocate(mask, updatecnts, states, scores, sampled)
        call ostats%kill
    end function check_conv3D

    subroutine append_stats( self, ostats )
        use CPlot2D_wrapper_module, only: plot2D
        class(convergence), intent(in) :: self
        class(oris),        intent(in) :: ostats
        type(oris)        :: os_prev, os
        type(string)      :: fname
        real, allocatable :: iter(:), inpl_dist(:), proj_dist(:)
        real, allocatable :: score(:), proj_overlap(:)
        integer    :: i,nl
        if( trim(params_glob%iterstats).ne.'yes' ) return
        fname = ITERSTATS_FILE
        if( file_exists(fname) )then
            nl = nlines(fname )
            call os_prev%new(nl,is_ptcl=.false.)
            call os_prev%read(fname)
            call os%new(nl+1,is_ptcl=.false.)
            do i =1,nl
                call os%transfer_ori(i,os_prev,i)
            enddo
            nl = nl+1
            call os_prev%kill
        else
            call os%new(1,is_ptcl=.false.)
            nl = 1
        endif
        call os%transfer_ori(nl,ostats,1)
        call os%write(fname)
        if( nl > 1 )then
            iter         = os%get_all('ITERATION')
            inpl_dist    = os%get_all('IN-PLANE_DIST')
            score        = os%get_all('SCORE')
            call plot2D(nl,iter,inpl_dist,   'iter_inpl_dist',   line=.true.,xtitle='Iterations',ytitle='Average in-plane distance (degs)')
            call plot2D(nl,iter,score,       'iter_score',       line=.true.,xtitle='Iterations',ytitle='Average score')
            if( os%isthere('ORIENTATION_OVERLAP') )then
                proj_overlap = os%get_all('ORIENTATION_OVERLAP')
                proj_dist    = os%get_all('DIST_BTW_BEST_ORIS')
                call plot2D(nl,iter,proj_dist,   'iter_proj_dist',   line=.true.,xtitle='Iterations',ytitle='Average orientation distance (degs)')
                call plot2D(nl,iter,proj_overlap,'iter_proj_overlap',line=.true.,xtitle='Iterations',ytitle='Average orientation overlap')
            endif
        endif
        call os%kill
        call fname%kill
    end subroutine append_stats

    subroutine plot_projdirs( self, ptcl_mask )
        use CPlot2D_wrapper_module
        class(convergence),   intent(in) :: self
        logical, allocatable, intent(in) :: ptcl_mask(:)
        type(string)                  :: title
        type(CPlot2D_type)            :: figure
        type(CDataSet_type)           :: center, axis
        type(CDataPoint_type)         :: p
        type(oris)                    :: os
        type(string)                  :: fname_eps, fname_pdf, ps2pdf_cmd
        real,             allocatable :: phi(:), psi(:), logpops(:)
        integer,          allocatable :: pops(:), projs(:), inds(:)
        real(dp) :: color, x,y, sz
        integer  :: iptcl, nptcls, nprojs, proj, l, iostat, ind
        if( trim(params_glob%iterstats).ne.'yes' ) return
        nptcls = size(ptcl_mask)
        projs  = nint(build_glob%spproj_field%get_all('proj'))
        nprojs = max(params_glob%nspace,maxval(projs,mask=ptcl_mask))
        ! gather populations & euler angles
        allocate(phi(nprojs),psi(nprojs),pops(nprojs), logpops(nprojs))
        pops = 0
        phi  = -2.
        psi  = -2.
        do iptcl = 1,nptcls
            if( ptcl_mask(iptcl) )then
                proj = projs(iptcl)
                pops(proj) = pops(proj) + 1
                if( phi(proj) < -1. )then
                    phi(proj) = build_glob%spproj_field%e1get(iptcl)
                    psi(proj) = build_glob%spproj_field%e2get(iptcl)
                endif
            endif
        enddo
        where( pops == 0 )
            phi = 0.
            psi = 0.
        end where
        ! raw populations output
        call os%new(nprojs,is_ptcl=.false.)
        do proj = 1,nprojs
            call os%set(proj,'pop',pops(proj))
            if( pops(proj) == 0 ) cycle
            call os%set(proj,'e1',phi(proj))
            call os%set(proj,'e2',psi(proj))
        enddo
        call os%write(string('projdir_pops.txt'))
        call os%kill
        ! sorting
        where( pops==0 )
            logpops = 0
        else where
            logpops = log10(real(pops)+0.1)
        end where
        logpops = logpops / maxval(logpops)
        inds = (/(proj,proj=1,nprojs)/)
        call hpsort(logpops,inds)
        ! Plot
        call CPlot2D__new(figure, 'Projection Directions Distribution'//C_NULL_CHAR)
        call CPlot2D__SetDrawXAxisGridLines(figure, C_FALSE)
        call CPlot2D__SetDrawYAxisGridLines(figure, C_FALSE)
        call CPlot2D__SetXAxisSize(figure, 400._c_double)
        call CPlot2D__SetYAxisSize(figure, 400._c_double)
        title = 'Phi'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(figure, title%to_char())
        title = 'Psi'//C_NULL_CHAR
        call CPlot2D__SetYAxisTitle(figure, title%to_char())
        call CPlot2D__SetDrawLegend(figure, C_FALSE)
        ! axes
        call CDataSet__new(axis)
        call CDataSet__SetDrawMarker(axis, C_FALSE)
        call CDataSet__SetDrawLine(axis, C_TRUE)
        call CDataSet__SetDatasetColor(axis, 0.d0,0.d0,0.d0)
        call CDataSet_addpoint(axis,-180., 0.)
        call CDataSet_addpoint(axis, 180., 0.)
        call CPlot2D__AddDataSet(figure, axis)
        call CDataSet__delete(axis)
        call CDataSet__new(axis)
        call CDataSet__SetDrawMarker(axis, C_FALSE)
        call CDataSet__SetDrawLine(axis, C_TRUE)
        call CDataSet__SetDatasetColor(axis, 0.d0,0.d0,0.d0)
        call CDataSet_addpoint(axis, 0., 0.)
        call CDataSet_addpoint(axis, 0., 180.)
        call CPlot2D__AddDataSet(figure, axis)
        call CDataSet__delete(axis)
        ! orientations
        do ind = 1,nprojs
            proj = inds(ind)
            if( pops(proj) == 0 ) cycle
            sz    = 8.d0 * real(logpops(ind),dp)
            color = 1.d0 - real(logpops(ind),dp)
            x     = real(phi(proj)-180.,dp)
            y     = real(psi(proj),dp)
            call CDataSet__new(center)
            call CDataSet__SetDrawMarker(center, C_TRUE)
            call CDataSet__SetMarkerSize(center, sz)
            call CDataSet__SetDatasetColor(center,color,color,1.d0)
            call CDataPoint__new2(x, y, p)
            call CDataSet__AddDataPoint(center, p)
            call CDataPoint__delete(p)
            call CPlot2D__AddDataSet(figure, center)
            call CDataSet__delete(center)
        enddo
        ! write
        fname_eps = 'iter_projdir_'//int2str_pad(params_glob%which_iter,3)//'.eps'//C_NULL_CHAR
        call CPlot2D__OutputPostScriptPlot(figure, fname_eps%to_char())
        call CPlot2D__delete(figure)
        l = fname_eps%strlen_trim()
        fname_eps = fname_eps%to_char([1,l-1]) ! removing trailing C NULL character
        fname_pdf = 'iter_projdir_'//int2str_pad(params_glob%which_iter,3)//'.pdf'
        ps2pdf_cmd = 'gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=512 -dDEVICEHEIGHTPOINTS=512 -sOutputFile='&
            //fname_pdf%to_char()//' '//fname_eps%to_char()
        call exec_cmdline(ps2pdf_cmd, suppress_errors=.true., exitstat=iostat)
        if( iostat == 0 ) call del_file(fname_eps)
    end subroutine plot_projdirs

    real function get( self, which )
        class(convergence), intent(in) :: self
        character(len=*),   intent(in) :: which
        get = 0.
        select case(which)
            case('iter')
                get = real(self%iteration)
            case('corr')
                get = self%score%avg
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
