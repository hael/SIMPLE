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
use simple_defs
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_commander_base, only: commander_base
use simple_commander_distr_wflows ! use all in there
use simple_filehandling           ! use all in there
use simple_jiffys                 ! use all in there
implicit none

public :: ini3D_from_cavgs_commander
private

type, extends(commander_base) :: ini3D_from_cavgs_commander
  contains
    procedure :: execute      => exec_ini3D_from_cavgs
end type ini3D_from_cavgs_commander

contains

    ! GENERATE INITIAL 3D MODEL FROM CLASS AVERAGES

    subroutine exec_ini3D_from_cavgs( self, cline )
        use simple_commander_comlin, only: symsrch_commander
        use simple_commander_volops, only: projvol_commander
        use simple_commander_rec,    only: recvol_commander
        use simple_strings,          only: int2str_pad, str2int
        class(ini3D_from_cavgs_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! constants
        logical,           parameter  :: DEBUG=.false.
        real,              parameter  :: LPLIMS(2) = [20.,10.], CENLP=50.
        integer,           parameter  :: MAXITS_INIT=50, MAXITS_REFINE=100
        integer,           parameter  :: STATE=1, NPROJS_SYMSRCH=50
        character(len=32), parameter  :: ITERFBODY = 'prime3Ddoc_'
        character(len=32), parameter  :: VOLFBODY  = 'recvol_state'
        ! distributed commanders
        type(prime3D_distr_commander) :: xprime3D_distr
        ! shared-mem commanders
        type(symsrch_commander)       :: xsymsrch
        type(recvol_commander)        :: xrecvol
        type(projvol_commander)       :: xprojvol
        ! command lines
        type(cmdline)                 :: cline_prime3D_init
        type(cmdline)                 :: cline_prime3D_refine1
        type(cmdline)                 :: cline_prime3D_refine2
        type(cmdline)                 :: cline_symsrch
        type(cmdline)                 :: cline_recvol
        type(cmdline)                 :: cline_projvol
        ! other variables
        type(params)                  :: p_master
        real                          :: iter
        character(len=2)              :: str_state
        character(len=STDLEN)         :: oritab, vol_iter
        logical                       :: srch4symaxis
        integer                       :: io_stat, pgrp_nr

        ! set cline defaults
        call cline%set('eo', 'no')
        ! make master parameters
        p_master = params(cline, checkdistr=.false.)
        ! set global state string
        str_state = int2str_pad(STATE,2)
        ! decide wether to search for the symmetry axis or put the point-group in from the start
        srch4symaxis = .false.
        if( p_master%pgrp .ne. 'c1' )then
            if(  p_master%pgrp(1:1).eq.'c'  .or. p_master%pgrp(1:1).eq.'C'&
            .or. p_master%pgrp(1:2).eq.'d2' .or. p_master%pgrp(1:2).eq.'D2' )then
                srch4symaxis = .true.
            endif
        endif

        ! prepare command lines from prototype master
        cline_prime3D_init    = cline
        cline_prime3D_refine1 = cline
        cline_prime3D_refine2 = cline
        cline_symsrch         = cline
        cline_recvol          = cline
        cline_projvol         = cline

        ! initialise command line parameters
        ! (1) PRIME3D_INIT
        call cline_prime3D_init%set('prg', 'prime3D')
        call cline_prime3D_init%set('ctf', 'no')
        call cline_prime3D_init%set('maxits', real(MAXITS_INIT))
        call cline_prime3D_init%set('dynlp', 'yes') ! better be explicit about the dynlp
        ! (2) PRIME3D REFINE STEP 1
        call cline_prime3D_refine1%set('prg', 'prime3D')
        call cline_prime3D_refine1%set('ctf', 'no')
        call cline_prime3D_refine1%set('maxits', real(MAXITS_REFINE))
        call cline_prime3D_refine1%set('dynlp', 'no') ! better be explicit about the dynlp
        call cline_prime3D_refine1%set('lp', LPLIMS(1))
        call cline_prime3D_refine1%set('refine', 'shc')
        ! (3) SYMMETRY AXIS SEARCH
        if( srch4symaxis )then
            ! need to replace original point-group flag with c1
            call cline_prime3D_init%set('pgrp', 'c1') 
            call cline_prime3D_refine1%set('pgrp', 'c1')
            ! symsrch
            call cline_symsrch%set('prg', 'symsrch')
            call cline_symsrch%delete('stk') ! volumetric symsrch
            call cline_symsrch%set('nptcls', real(NPROJS_SYMSRCH))
            call cline_symsrch%set('nspace', real(NPROJS_SYMSRCH))
            call cline_symsrch%set('cenlp', CENLP)
            call cline_symsrch%set('outfile', 'symdoc.txt')
            call cline_symsrch%set('lp', LPLIMS(1))
            if( cline%defined('nthr_master') )then
                call cline_symsrch%set('nthr', real(p_master%nthr_master))
            endif
            ! (3.5) RECONSTRUCT SYMMETRISED VOLUME
            call cline_recvol%set('prg', 'recvol')
            call cline_recvol%set('trs', 5.) ! to assure that shifts are being used
            call cline_recvol%set('ctf', 'no')
            call cline_recvol%set('oritab', 'symdoc.txt')
            if( cline%defined('nthr_master') )then
                call cline_recvol%set('nthr', real(p_master%nthr_master))
            endif
            ! 2nd refinement step now uses the symmetrised vol and doc
            call cline_prime3D_refine2%set('oritab', 'symdoc.txt')
            call cline_prime3D_refine2%set('vol1', 'rec_sym'//p_master%ext)
        endif
        ! (4) PRIME3D REFINE STEP 2
        call cline_prime3D_refine2%set('prg', 'prime3D')
        call cline_prime3D_refine2%set('ctf', 'no')
        call cline_prime3D_refine2%set('maxits', real(MAXITS_REFINE))
        call cline_prime3D_refine2%set('dynlp', 'no') ! better be explicit about the dynlp
        call cline_prime3D_refine2%set('lp', LPLIMS(2))
        call cline_prime3D_refine2%set('refine', 'shc')
        ! (5) RE-PROJECT VOLUME
        call cline_projvol%set('prg', 'projvol')
        call cline_projvol%set('outstk', 'reprojs'//p_master%ext)

        ! execute commanders
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> INITIAL 3D MODEL GENERATION WITH PRIME3D'
        write(*,'(A)') '>>>'
        call xprime3D_distr%execute(cline_prime3D_init)
        iter = cline_prime3D_init%get_rarg('endit')
        call set_iter_dependencies
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> FIRST PRIME3D REFINEMENT STEP, LP=20, REFINE=SHC'
        write(*,'(A)') '>>>'
        call cline_prime3D_refine1%set('startit', iter + 1.0)
        call cline_prime3D_refine1%set('oritab', trim(oritab))
        call cline_prime3D_refine1%set('vol1', trim(vol_iter))
        call xprime3D_distr%execute(cline_prime3D_refine1)
        iter = cline_prime3D_refine1%get_rarg('endit')
        call set_iter_dependencies
        if( srch4symaxis )then
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> SYMMETRY AXIS SEARCH'
            write(*,'(A)') '>>>'
            call cline_symsrch%set('oritab', trim(oritab))
            call cline_symsrch%set('vol1', trim(vol_iter))
            call xsymsrch%execute(cline_symsrch)
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> RECONSTRUCTION OF SYMMETRISED VOLUME'
            write(*,'(A)') '>>>'
            call xrecvol%execute(cline_recvol)
            call rename(trim(volfbody)//trim(str_state)//p_master%ext, 'rec_sym'//p_master%ext)  
        else
            ! 2nd refinement step needs to use iter dependent vol/oritab
            call cline_prime3D_refine2%set('oritab', trim(oritab))
            call cline_prime3D_refine2%set('vol1', trim(vol_iter))
        endif
        call cline_prime3D_refine2%set('startit', iter + 1.0)
        call xprime3D_distr%execute(cline_prime3D_refine2)
        iter = cline_prime3D_refine2%get_rarg('endit')
        call set_iter_dependencies
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> RE-PROJECTION OF THE FINAL VOLUME'
        write(*,'(A)') '>>>'
        call cline_projvol%set('vol1', trim(vol_iter))
        call cline_projvol%set('oritab', trim(oritab))
        call xprojvol%execute(cline_projvol)

        ! end gracefully
        call simple_end('**** SIMPLE_INI3D_FROM_CAVGS NORMAL STOP ****')

        contains

            subroutine set_iter_dependencies
                character(len=3) :: str_iter
                str_iter = int2str_pad(nint(iter),3)
                oritab   = trim(ITERFBODY)//trim(str_iter)//'.txt'
                vol_iter = trim(VOLFBODY)//trim(str_state)//'_iter'//trim(str_iter)//p_master%ext                
            end subroutine set_iter_dependencies

    end subroutine exec_ini3D_from_cavgs

end module simple_commander_hlev_wflows
