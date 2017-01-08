!==Class simple_commander_hlev_wflows
!
! This class contains commanders responsible for execution of high-level workflows in SIMPLE. This class provides 
! the glue between the reciver (main reciever is simple_distr_exec) and the abstract action, which is simply execute 
! (defined by the base class: simple_commander_base).
!
! The code is hlevibuted with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Hans Elmlund 2017
!
module simple_commander_hlev_wflows
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_params,         only: params
use simple_commander_distr_wflows
implicit none

public :: ini3D_from_cavgs_commander
private

type, extends(commander_base) :: ini3D_from_cavgs_commander
  contains
    procedure :: execute      => exec_ini3D_from_cavgs
end type ini3D_from_cavgs_commander

contains

    ! GENERATE INITIAL 3D MODEL FROM CLASS AVERAGES
    ! (1) prime3D with pgrp=c1 and 

    subroutine exec_ini3D_from_cavgs( self, cline )
        use simple_commander_comlin, only: symsrch_commander
        use simple_commander_volops, only: projvol_commander
        ! constants
        logical, parameter            :: DEBUG=.false.
        real,    parameter            :: LPLIMS(2) = [20.,10.]
        integer, parameter            :: MAXITS1=50, MAXITS2=100
        ! distributed commanders
        type(prime3D_distr_commander) :: xprime3D_distr
        ! shared-mem commanders
        type(symsrch_commander)       :: xsymsrch
        type(projvol_commander)       :: xprojvol
        ! command lines
        type(cmdline)                 :: cline_prime3D_init
        type(cmdline)                 :: cline_prime3D_refine1
        type(cmdline)                 :: cline_prime3D_refine2
        type(cmdline)                 :: cline_symsrch
        type(cmdline)                 :: cline_projvol
        ! other variables
        type(params)                  :: p_master

        




    end subroutine exec_ini3D_from_cavgs

end module simple_commander_hlev_wflows
