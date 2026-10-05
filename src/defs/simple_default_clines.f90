!@descr: setting command line defaults
module simple_default_clines
use simple_defs
use simple_cmdline,       only: cmdline
use simple_estimate_ssnr, only: mskdiam2lplimits
implicit none

! the automask2D defaults, also for callers that set up their masking without a command line
integer, parameter :: AUTOMASK2D_NGROW  = 3
real,    parameter :: AUTOMASK2D_WINSZ  = 5.
real,    parameter :: AUTOMASK2D_AMSKLP = 20.
integer, parameter :: AUTOMASK2D_EDGE   = 6

contains

    subroutine set_automask2D_defaults( cline )
        class(cmdline), intent(inout) :: cline
        if( .not. cline%defined('ngrow')  ) call cline%set('ngrow',  AUTOMASK2D_NGROW)
        if( .not. cline%defined('winsz')  ) call cline%set('winsz',  AUTOMASK2D_WINSZ)
        if( .not. cline%defined('amsklp') ) call cline%set('amsklp', AUTOMASK2D_AMSKLP)
        if( .not. cline%defined('edge')   ) call cline%set('edge',   AUTOMASK2D_EDGE)
    end subroutine set_automask2D_defaults

    subroutine set_refine2D_defaults( cline )
        class(cmdline), intent(inout) :: cline
        real :: mskdiam, lpstart, lpstop, lpcen
        mskdiam = cline%get_rarg('mskdiam')
        call mskdiam2lplimits(cline%get_rarg('mskdiam'), lpstart, lpstop, lpcen)
        if( .not. cline%defined('mkdir')        ) call cline%set('mkdir',        'yes')
        if( .not. cline%defined('oritype')      ) call cline%set('oritype',   'ptcl2D')
        if( .not. cline%defined('lpstart')      ) call cline%set('lpstart',    lpstart)
        if( .not. cline%defined('lpstop')       ) call cline%set('lpstop',      lpstop)
        if( .not. cline%defined('cenlp')        ) call cline%set('cenlp',        lpcen)
        if( .not. cline%defined('maxits')       ) call cline%set('maxits',          30)
        if( .not. cline%defined('autoscale')    ) call cline%set('autoscale',    'yes')
        if( .not. cline%defined('cls_init')     ) call cline%set('cls_init',    'ptcl')
        if( .not. cline%defined('center_type')  ) call cline%set('center_type', 'mass')
        if( .not. cline%defined('refine')       ) call cline%set('refine',   'snhc_smpl')
        if( .not. cline%defined('extr_lim')     ) call cline%set('extr_lim', MAX_EXTRLIM2D)
        if( .not. cline%defined('restore_cavgs')) call cline%set('restore_cavgs','yes')
        ! 2D objective function section
        if( .not. cline%defined('objfun')       ) call cline%set('objfun',    'euclid')
        if( .not. cline%defined('ml_reg')       ) call cline%set('ml_reg',       'yes')
        call set_automask2D_defaults( cline )
    end subroutine set_refine2D_defaults

end module simple_default_clines
