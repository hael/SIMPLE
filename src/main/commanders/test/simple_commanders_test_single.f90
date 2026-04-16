!@descr: tests for single 
module simple_commanders_test_single
use simple_commanders_api
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_atoms_stats
  contains
    procedure :: execute      => exec_test_atoms_stats
end type commander_test_atoms_stats

type, extends(commander_base) :: commander_test_detect_atoms
  contains
    procedure :: execute      => exec_test_detect_atoms
end type commander_test_detect_atoms

type, extends(commander_base) :: commander_test_simulate_nanoparticle
  contains
    procedure :: execute      => exec_test_simulate_nanoparticle
end type commander_test_simulate_nanoparticle

type, extends(commander_base) :: commander_test_single_workflow
  contains
    procedure :: execute      => exec_test_single_workflow
end type commander_test_single_workflow

integer, parameter :: BOX          = 160
integer, parameter :: MOLDIAM      = 20
integer, parameter :: NTHR         = 24

contains

subroutine exec_test_atoms_stats( self, cline )
    use simple_commanders_atoms, only: commander_detect_atoms
    use simple_commanders_sim,   only: commander_simulate_nanoparticle
    use simple_commanders_atoms, only: commander_atoms_stats
    class(commander_test_atoms_stats), intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    type(cmdline)                         :: cline_sim, cline_detat, cline_atstats
    type(parameters)                      :: params
    type(commander_simulate_nanoparticle) :: xsim_nptcl
    type(commander_detect_atoms)          :: xdetat
    type(commander_atoms_stats)           :: xatstats
    write(logfhandle,'(a)') '>>> TEST_ATOMS_STATS:'
    call params%new(cline)
    call cline_sim%set('prg',      'simulate_nanoparticle')
    call cline_sim%set('box',                          BOX)
    call cline_sim%set('smpd',                 params%smpd)
    call cline_sim%set('moldiam',                  MOLDIAM)
    call cline_sim%set('element',           params%element)
    call cline_sim%set('nthr',                        NTHR)
    call xsim_nptcl%execute(cline_sim)
    call cline_detat%set('prg',             'detect_atoms')
    call cline_detat%set('vol1',              'outvol.mrc')
    call cline_detat%set('smpd',               params%smpd)
    call cline_detat%set('element',         params%element)
    call cline_detat%set('nthr',                      NTHR)
    call xdetat%execute(cline_detat)
    call cline_atstats%set('prg',             'atoms_stats')
    call cline_atstats%set('vol1',             'outvol.mrc')
    call cline_atstats%set('vol2',          'outvol_CC.mrc')
    call cline_atstats%set('pdbfile',     'outvol_ATMS.pdb')
    call cline_atstats%set('smpd',              params%smpd)
    call cline_atstats%set('element',        params%element)
    call cline_atstats%set('nthr',                     NTHR)
    call xatstats%execute(cline_atstats)
    call simple_end('**** SIMPLE_TEST_ATOMS_STATS NORMAL STOP ****')
end subroutine exec_test_atoms_stats

subroutine exec_test_detect_atoms( self, cline )
    use simple_commanders_atoms, only: commander_detect_atoms
    use simple_commanders_sim,   only: commander_simulate_nanoparticle
    class(commander_test_detect_atoms), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    type(cmdline)                         :: cline_sim, cline_detat
    type(parameters)                      :: params
    type(commander_simulate_nanoparticle) :: xsim_nptcl
    type(commander_detect_atoms)          :: xdetat
    write(logfhandle,'(a)') '>>> TEST_DETECT_ATOMS:'
    call params%new(cline_sim)
    call cline_sim%set('prg',      'simulate_nanoparticle')
    call cline_sim%set('box',                          BOX)
    call cline_sim%set('smpd',                 params%smpd)
    call cline_sim%set('moldiam',                  MOLDIAM)
    call cline_sim%set('element',           params%element)
    call cline_sim%set('nthr',                        NTHR)
    call xsim_nptcl%execute(cline_sim)
    call cline_detat%set('prg',             'detect_atoms')
    call cline_detat%set('vol1',              'outvol.mrc')
    call cline_detat%set('smpd',               params%smpd)
    call cline_detat%set('element',         params%element)
    call cline_detat%set('nthr',                      NTHR)
    call xdetat%execute(cline_detat)
    call simple_end('**** SIMPLE_TEST_DETECT_ATOMS NORMAL STOP ****')
end subroutine exec_test_detect_atoms

subroutine exec_test_simulate_nanoparticle( self, cline )
    use simple_commanders_sim, only: commander_simulate_nanoparticle
    class(commander_test_simulate_nanoparticle), intent(inout) :: self
    class(cmdline),                              intent(inout) :: cline
    type(cmdline)                         :: cline_sim
    type(parameters)                      :: params
    type(commander_simulate_nanoparticle) :: xsim_nptcl
    write(logfhandle,'(a)') '>>> TEST_SIMULATE_NANOPARTICLE:'
    call params%new(cline)
    call cline_sim%set('prg',      'simulate_nanoparticle')
    call cline_sim%set('box',                          BOX)
    call cline_sim%set('smpd',                 params%smpd)
    call cline_sim%set('moldiam',                  MOLDIAM)
    call cline_sim%set('element',           params%element)
    call cline_sim%set('nthr',                        NTHR)
    call xsim_nptcl%execute(cline_sim)
    call simple_end('**** SIMPLE_TEST_SIMULATE_NANOPARTICLE NORMAL STOP ****')
end subroutine exec_test_simulate_nanoparticle

subroutine exec_test_single_workflow( self, cline )
    use single_commanders_nano2D,       only: commander_analysis2D_nano
    !use simple_commanders_sim,          only: commander_simulate_nanoparticle, commander_simulate_particles
    use simple_commanders_sim,          only: commander_simulate_nanoparticle
    use simple_commanders_reproject,    only: commander_reproject
    use simple_commanders_project_ptcl, only: commander_import_particles
    use simple_commanders_project_core, only: commander_new_project
    use single_commanders_nano3D,       only: commander_autorefine3D_nano
    class(commander_test_single_workflow), intent(inout) :: self
    class(cmdline),                        intent(inout) :: cline
    ! single_exec prg=tseries_import filetab=filetab.txt cs=0 fraca=0.3 kv=300 smpd=0.358
    ! nohup single_exec prg=tseries_motion_correct nparts=6 nthr=12 > MCORR &
    ! nohup single_exec prg=tseries_track_particles boxfile=positions_all.box fbody=tracked_ptcl nthr=32 > TRACK_MCORR &
    ! single_exec prg=import_particles cs=-0.01 fraca=0.4 kv=300 smpd=0.358 stk=../subtr_den_cong.mrc ctf=no
    ! single_exec prg=analysis2D_nano element=Pt nthr=16
    ! single_exec prg=map_cavgs_selection stk2=2_analysis2D_nano/selected.spi prune=yes
    ! nohup single_exec prg=autorefine3D_nano vol1=startvol.mrc element=Pt smpd=0.358 pgrp=c1 lp=1.5 mskdiam=50 nthr=16 projfile=3_selection/reproc.simple > AUTOREFINE &
    type(cmdline)                         :: cline_sim, cline_simptcls, cline_nproj,cline_imptcls, cline_an2Dnano, cline_aref3Dnano
    type(parameters)                      :: params
    type(commander_simulate_nanoparticle) :: xsim_nptcl
    type(commander_reproject)             :: xsim_ptcls
    type(commander_import_particles)      :: ximptcls
    type(commander_analysis2D_nano)       :: xan2Dnano
    type(commander_new_project)           :: xnproj
    type(commander_autorefine3D_nano)     :: xaref3Dnano
    type(string)                          :: projname, projfile
    integer, parameter                    :: NPTCLS_SIM = 10000
    write(logfhandle,'(a)') '>>> TEST_SINGLE_WORKFLOW:'
    projname='test_single_workflow'
    call params%new(cline)
    call cline_sim%set('prg',               'simulate_nanoparticle')
    call cline_sim%set('box',                                   BOX)
    call cline_sim%set('smpd',                          params%smpd)
    call cline_sim%set('moldiam',                           MOLDIAM)
    call cline_sim%set('element',                    params%element)
    call cline_sim%set('nthr',                                 NTHR)
    call xsim_nptcl%execute(cline_sim)
    call cline_simptcls%set('prg',                      'reproject')
    call cline_simptcls%set('pgrp',                            'c1')
    call cline_simptcls%set('vol1',                    'outvol.mrc')
    call cline_simptcls%set('smpd',                     params%smpd)
    call cline_simptcls%set('nspace',                    NPTCLS_SIM)
    call cline_simptcls%set('mskdiam',                           20)
    call cline_simptcls%set('outstk',     'simulated_particles.mrc')
    call cline_simptcls%set('nthr',                            NTHR)
    ! call cline_simptcls%set('prg',             'simulate_particles')
    ! call cline_simptcls%set('vol1',                    'outvol.mrc')
    ! call cline_simptcls%set('smpd',                     params%smpd)
    ! call cline_simptcls%set('cs',                             -0.01)
    ! call cline_simptcls%set('kv',                               300)
    ! call cline_simptcls%set('fraca',                            0.4)
    ! call cline_simptcls%set('nthr',                              24)
    ! call cline_simptcls%set('nptcls',                    NPTCLS_SIM)
    ! call cline_simptcls%set('pgrp',                            'c1')
    ! call cline_simptcls%set('snr',                               1.)
    ! call cline_simptcls%set('ctf',                             'no')
    ! call cline_simptcls%set('dferr',                            0.0)
    ! call cline_simptcls%set('sherr',                            0.0)
    ! call cline_simptcls%set('even',                            'on')
    ! call cline_simptcls%set('defocus',                          0.0)
    call xsim_ptcls%execute(cline_simptcls)
    projfile = projname%to_char()//'.simple'
    call cline_nproj%set('prg',                       'new_project')
    call cline_nproj%set('projname',             projname%to_char())
    call xnproj%execute(cline_nproj)
    call cline_imptcls%set('prg',                'import_particles')
    call cline_imptcls%set('projfile',           projfile%to_char())
    call cline_imptcls%set('stk',      '../simulated_particles.mrc')
    call cline_imptcls%set('smpd',                      params%smpd)
    call cline_imptcls%set('cs  ',                            -0.01)
    call cline_imptcls%set('fraca',                             0.4)
    call cline_imptcls%set('kv',                                300)
    call cline_imptcls%set('ctf',                              'no')
    call ximptcls%execute(cline_imptcls)
    call cline_an2Dnano%set('prg',                'analysis2D_nano')
    call cline_an2Dnano%set('projfile',          projfile%to_char())
    call cline_an2Dnano%set('element',               params%element)
    call cline_an2Dnano%set('nthr',                            NTHR)
    call xan2Dnano%execute(cline_an2Dnano)
    call cline_aref3Dnano%set('prg',            'autorefine3D_nano')
    call cline_aref3Dnano%set('projfile',        projfile%to_char())
    call cline_aref3Dnano%set('vol1',               '../outvol.mrc')
    call cline_aref3Dnano%set('smpd',                   params%smpd)
    call cline_aref3Dnano%set('element',             params%element)
    call cline_aref3Dnano%set('nthr',                          NTHR)
    call cline_aref3Dnano%set('pgrp',                          'c1')
    call cline_aref3Dnano%set('lp',                             1.0)  
    call cline_aref3Dnano%set('mskdiam',                         25)
    call cline_aref3Dnano%set('mskdiam',                         50)
    call xaref3Dnano%execute(cline_aref3Dnano)
    call simple_end('**** SIMPLE_TEST_SINGLE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_single_workflow

end module simple_commanders_test_single
