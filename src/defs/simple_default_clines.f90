module simple_default_clines
use simple_cmdline,       only: cmdline
use simple_estimate_ssnr, only: mskdiam2lplimits
implicit none

contains

    subroutine set_automask2D_defaults( cline )
        class(cmdline), intent(inout) :: cline
        if( .not. cline%defined('ngrow')  ) call cline%set('ngrow',   3.)
        if( .not. cline%defined('winsz')  ) call cline%set('winsz',   5.)
        if( .not. cline%defined('amsklp') ) call cline%set('amsklp', 20.)
        if( .not. cline%defined('edge')   ) call cline%set('edge',    6.)
    end subroutine set_automask2D_defaults

    subroutine set_cluster2D_defaults( cline )
        class(cmdline), intent(inout) :: cline
        real :: mskdiam, lpstart, lpstop, lpcen
        mskdiam = cline%get_rarg('mskdiam')
        call mskdiam2lplimits(cline%get_rarg('mskdiam'), lpstart, lpstop, lpcen)
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('oritype')   ) call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('lpstart')   ) call cline%set('lpstart',  lpstart)
        if( .not. cline%defined('lpstop')    ) call cline%set('lpstop',    lpstop)
        if( .not. cline%defined('cenlp')     ) call cline%set('cenlp',      lpcen)
        if( .not. cline%defined('maxits')    ) call cline%set('maxits',      30. )
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale',  'yes')
        if( .not.cline%defined('refine')     ) call cline%set('refine',    'snhc')
        if( .not. cline%defined('wiener')    ) call cline%set('wiener',    'full')
        if( .not. cline%defined('ptclw')     ) call cline%set('ptclw',       'no')
        ! 2D objective function section
        if( .not. cline%defined('objfun')    ) call cline%set('objfun',  'euclid')
        if( .not. cline%defined('ml_reg')    ) call cline%set('ml_reg',      'no')
        call set_automask2D_defaults( cline )
    end subroutine set_cluster2D_defaults

end module simple_default_clines
