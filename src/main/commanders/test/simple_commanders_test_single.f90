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

type, extends(commander_base) :: commander_test_simulate_nanoprticle
  contains
    procedure :: execute      => exec_test_simulate_nanoprticle
end type commander_test_simulate_nanoprticle

type, extends(commander_base) :: commander_test_single_workflow
  contains
    procedure :: execute      => exec_test_single_workflow
end type commander_test_single_workflow

contains

subroutine exec_test_atoms_stats( self, cline )
    class(commander_test_atoms_stats), intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_ATOMS_STATS NORMAL STOP ****')
end subroutine exec_test_atoms_stats

subroutine exec_test_detect_atoms( self, cline )
    class(commander_test_detect_atoms), intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_DETECT_ATOMS NORMAL STOP ****')
end subroutine exec_test_detect_atoms

subroutine exec_test_simulate_nanoprticle( self, cline )
    class(commander_test_simulate_nanoprticle), intent(inout) :: self
    class(cmdline),                            intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_SIMULATE_NANOPARTICLE NORMAL STOP ****')
end subroutine exec_test_simulate_nanoprticle

subroutine exec_test_single_workflow( self, cline )
    class(commander_test_single_workflow), intent(inout) :: self
    class(cmdline),                        intent(inout) :: cline
    ! single_exec prg=tseries_import filetab=filetab.txt cs=0 fraca=0.3 kv=300 smpd=0.358
    ! nohup single_exec prg=tseries_motion_correct nparts=6 nthr=12 > MCORR &
    ! nohup single_exec prg=tseries_track_particles boxfile=positions_all.box fbody=tracked_ptcl nthr=32 > TRACK_MCORR &
    ! single_exec prg=import_particles cs=-0.01 fraca=0.4 kv=300 smpd=0.358 stk=../subtr_den_cong.mrc ctf=no
    ! single_exec prg=analysis2D_nano element=Pt nthr=16
    ! single_exec prg=map_cavgs_selection stk2=2_analysis2D_nano/selected.spi prune=yes
    ! nohup single_exec prg=autorefine3D_nano vol1=startvol.mrc element=Pt smpd=0.358 pgrp=c1 lp=1.5 mskdiam=50 nthr=16 projfile=3_selection/reproc.simple > AUTOREFINE &
    call simple_end('**** SIMPLE_TEST_SINGLE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_single_workflow

end module simple_commanders_test_single
