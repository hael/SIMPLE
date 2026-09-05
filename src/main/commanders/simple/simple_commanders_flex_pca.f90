!@descr: projection-aware covariance heterogeneity commander
module simple_commanders_flex_pca
use simple_commanders_api
implicit none

#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_flex_pca
  contains
    procedure :: execute => exec_flex_pca
end type commander_flex_pca

contains

    subroutine exec_flex_pca( self, cline )
        use simple_flex_pca_model, only: run_flex_pca
        use simple_flex_pca_distr, only: flex_pca_distr_init, flex_pca_distr_kill
        use simple_ptcl_cache,     only: ptcl_cache_ensure
        class(commander_flex_pca), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        ! a worker runs in the master's directory and must not descend into one of its own.
        ! NOT merge('no ','yes',..): merge pads the shorter branch, yielding 'no ' with a trailing blank.
        if( .not.cline%defined('mkdir') )then
            if( cline%defined('part') )then
                call cline%set('mkdir','no')
            else
                call cline%set('mkdir','yes')
            endif
        endif
        if( .not.cline%defined('oritype') )     call cline%set('oritype','ptcl3D')
        if( .not.cline%defined('nstates') )     call cline%set('nstates',1)
        ! npreimages is a CEILING, not a target: the two-gate merge collapses indistinct states below
        ! it. Under preimage_auto the key is deliberately LEFT UNDEFINED when the user did not pin
        ! one, because run_flex_pca reads cline%defined('npreimages') to decide whether it may raise
        ! the ceiling to the auto value -- injecting a default here would make that test always true.
        if( .not.cline%defined('npreimages') .and. .not.flex_pca_auto_states(cline) ) &
            &call cline%set('npreimages',16)
        if( .not.cline%defined('neigs') )       call cline%set('neigs',16)
        ! an UPPER BOUND, not an iteration count: the probe stops itself on COV_PROBE_CONV
        if( .not.cline%defined('n_probe_iters') ) call cline%set('n_probe_iters',5)
        if( .not.cline%defined('box_crop') )    call cline%set('box_crop',64)
        call derive_flex_pca_band(cline)
        if( .not.cline%defined('ptcl_src') )    call cline%set('ptcl_src','raw')
        if( .not.cline%defined('objfun') )      call cline%set('objfun','euclid')
        if( .not.cline%defined('outvol') )      call cline%set('outvol','flex_pca_state_001.mrc')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        call ensure_canonical_sigma_state(params, build, cline)
        ! cache=yes: build the downscaled particle cache ONCE on the master. flex re-reads the same
        ! particles on every EM iteration and relaunches its workers each round, so an in-process
        ! cache cannot help; the on-disk one is adopted by each worker (they run in the master's
        ! directory, so the cache path resolves identically).
        if( .not. cline%defined('part') ) call ptcl_cache_ensure(params, build, cline)
        call flex_pca_distr_init(params,cline)
        call run_flex_pca(params,build,cline)
        call flex_pca_distr_kill(params)
        call build%kill_general_tbox
        if( cline%defined('part') )then
            call simple_end('**** SIMPLE_FLEX_PCA WORKER NORMAL STOP ****',print_simple=.false.)
        else
            call simple_end('**** SIMPLE_FLEX_PCA NORMAL STOP ****')
        endif
    end subroutine exec_flex_pca

    subroutine ensure_canonical_sigma_state( params, build, cline )
        use, intrinsic :: iso_fortran_env, only: int64
        use simple_commanders_euclid, only: commander_calc_pspec
        use simple_sigma2_state, only: sigma2_state_project_layout_digest, sigma2_state_validate_identity
        use simple_sigma2_state_file, only: sigma2_state_validate_file, SIGMA2_GROUP_GLOBAL, &
            &SIGMA2_GROUP_STACK, SIGMA2_STATE_COMMITTED
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(commander_calc_pspec) :: xcalc_pspec
        type(cmdline) :: cline_pspec
        type(string) :: state_path
        integer(int64) :: layout_digest
        integer :: iptcl, ngroups, status
        logical :: found, rebuild
        character(len=STDLEN) :: message
        if( .not. params%l_sigma_canonical ) return
        if( params%cc_objfun /= OBJFUN_EUCLID ) return
        rebuild = .true.
        call build%spproj%get_sigma2_state_path(state_path, found)
        if( found )then
            call sigma2_state_validate_file(state_path%to_char(), status, message, deep=.true.)
            if( status == 0 )then
                layout_digest = sigma2_state_project_layout_digest(build%spproj, build%spproj_field)
                if( params%l_sigma_glob )then
                    call sigma2_state_validate_identity(state_path%to_char(), params%box, params%smpd, 1, &
                        &fdim(params%box)-1, params%nptcls, layout_digest, status, message, &
                        &expected_state=SIGMA2_STATE_COMMITTED, expected_grouping=SIGMA2_GROUP_GLOBAL, &
                        &expected_ngroups=1)
                else
                    ngroups = 0
                    do iptcl = 1, params%nptcls
                        if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                        ngroups = max(ngroups, build%spproj_field%get_int(iptcl, 'stkind'))
                    enddo
                    call sigma2_state_validate_identity(state_path%to_char(), params%box, params%smpd, 1, &
                        &fdim(params%box)-1, params%nptcls, layout_digest, status, message, &
                        &expected_state=SIGMA2_STATE_COMMITTED, expected_grouping=SIGMA2_GROUP_STACK, &
                        &expected_ngroups=ngroups)
                endif
                rebuild = status /= 0
            endif
        endif
        if( rebuild )then
            write(logfhandle,'(A)') '>>> FLEX_PCA SIGMA: initializing missing or stale canonical state from particle power'
            cline_pspec = cline
            call cline_pspec%set('prg', 'calc_pspec')
            call cline_pspec%set('mkdir', 'no')
            call cline_pspec%delete('part')
            call xcalc_pspec%execute(cline_pspec)
            call build%spproj%read_segment('projinfo', params%projfile)
            call cline_pspec%kill
        endif
        call state_path%kill
    end subroutine ensure_canonical_sigma_state

    !> Whether this run lets the data set the state count. Read straight off the cmdline because it
    !! is needed before params%new.
    logical function flex_pca_auto_states( cline ) result( auto )
        class(cmdline), intent(inout) :: cline
        type(string) :: val
        auto = .false.
        if( .not. cline%defined('preimage_auto') ) return
        val  = cline%get_carg('preimage_auto')
        auto = trim(val%to_char()) == 'yes'
        call val%kill
    end function flex_pca_auto_states

    ! Resolve lp and box_rec from project geometry; neither overrides an explicit command-line value.
    ! Called on the master before params%new so the workers inherit the resolved numbers.
    subroutine derive_flex_pca_band( cline )
        class(cmdline), intent(inout) :: cline
        type(sp_project) :: spproj
        type(string)     :: projfile
        real             :: smpd, smpd_crop, lp_here
        integer          :: box, box_crop
        ! 2*smpd_crop is the working Nyquist; the 1.25 safety factor keeps the band clear of it
        real, parameter  :: COV_LP_OVER_NYQUIST = 2.5
        if( cline%defined('lp') .and. cline%defined('box_rec') ) return
        if( .not. cline%defined('projfile') ) return
        projfile = cline%get_carg('projfile')
        if( .not. file_exists(projfile) ) return
        call spproj%read_segment('stk', projfile)
        box  = spproj%get_box()
        smpd = spproj%get_smpd()
        call spproj%kill
        if( box < 1 .or. smpd <= 0. ) return
        box_crop = box
        if( cline%defined('box_crop') ) box_crop = cline%get_iarg('box_crop')
        if( box_crop < 1 ) box_crop = box
        if( .not. cline%defined('lp') )then
            smpd_crop = real(box) / real(box_crop) * smpd
            lp_here   = COV_LP_OVER_NYQUIST * smpd_crop
            call cline%set('lp', lp_here)
            write(logfhandle,'(A,F7.3,A,F6.3,A,I0,A)') '>>> FLEX_PCA derived lp=', lp_here, &
                &' A from smpd_crop=', smpd_crop, ' A (box ', box, ' -> box_crop)'
        endif
        ! the covariance is fitted at box_crop, but the state maps need not inherit that limit
        if( .not. cline%defined('box_rec') )then
            call cline%set('box_rec', box)
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA derived box_rec=', box, &
                &' (native box; state maps not capped at the covariance Nyquist)'
        endif
        call projfile%kill
    end subroutine derive_flex_pca_band

end module simple_commanders_flex_pca
