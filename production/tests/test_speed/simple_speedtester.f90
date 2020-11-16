module simple_speedtester
use simple_image,              only: image
use simple_commander_refine3D, only: refine3D_commander_distr
use simple_commander_volops,   only: reproject_commander
use simple_cmdline,            only: cmdline
implicit none

integer, parameter            :: BOX=200, SQRAD=40, NSPACE=50
real,    parameter            :: SMPD=1.1, MSK=50., TRS=5.0
type(image)                   :: cube, square
type(reproject_commander)     :: xproject
type(refine3D_commander_distr) :: xprime3D
type(cmdline)                 :: cline_project, cline_refine3D

contains

    subroutine exec_speedtest( cline, be_verbose )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
        call setup_testenv( cline )
    end subroutine exec_speedtest

    subroutine setup_testenv( cline )
        class(cmdline), intent(inout) :: cline
        real    :: nthr
        integer :: ldim(3)
        ldim = [BOX,BOX,BOX]
        call cube%new(ldim, smpd)
        call cube%square(SQRAD)
        call cube%write('speedtester_cube.mrc')
        nthr = cline%get_rarg('nthr')
        call cline_project%set('vol1',   'speedtester_cube.mrc')
        call cline_project%set('outstk', 'speedtester_cubes.mrc')
        call cline_project%set('nspace',  real(NSPACE))
        call cline_project%set('nthr',    nthr)
        call xproject%execute(cline_project)
        call cline_refine3D%set('prg',    'prime3D')
        call cline_refine3D%set('stk',    'speedtester_cubes.mrc')
        call cline_refine3D%set('smpd',    real(SMPD))
        call cline_refine3D%set('msk',     real(MSK))
        call cline_refine3D%set('ctf',    'no')
        call cline_refine3D%set('pgrp',   'c1')
        call cline_refine3D%set('nparts',  1.0)
        call cline_refine3D%set('vol1',   'speedtester_cube.mrc')
        call cline_refine3D%set('trs',     TRS)
        call cline_refine3D%set('lp',      6.0)
        call cline_refine3D%set('nthr',    nthr)
        call cline_refine3D%set('eo',      'no')
        call cline_refine3D%set('maxits',  1.0)
        call xprime3D%execute(cline_refine3D)
    end subroutine setup_testenv

end module simple_speedtester
