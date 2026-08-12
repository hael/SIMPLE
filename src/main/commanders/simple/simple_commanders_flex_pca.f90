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
        if( .not.cline%defined('npreimages') )  call cline%set('npreimages',7)
        if( .not.cline%defined('neigs') )       call cline%set('neigs',16)
        if( .not.cline%defined('ncols') )       call cline%set('ncols',128)
        ! an UPPER BOUND, not an iteration count: the probe stops itself on COV_PROBE_CONV
        if( .not.cline%defined('n_probe_iters') ) call cline%set('n_probe_iters',5)
        if( .not.cline%defined('box_crop') )    call cline%set('box_crop',64)
        call derive_flex_pca_band(cline)
        if( .not.cline%defined('ptcl_src') )    call cline%set('ptcl_src','raw')
        if( .not.cline%defined('objfun') )      call cline%set('objfun','euclid')
        if( .not.cline%defined('outvol') )      call cline%set('outvol','flex_pca_state_001.mrc')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
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
