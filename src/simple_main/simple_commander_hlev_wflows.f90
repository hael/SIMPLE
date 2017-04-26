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
use simple_qsys_env,       only: qsys_env
use simple_commander_distr_wflows ! use all in there
use simple_filehandling           ! use all in there
use simple_jiffys                 ! use all in there
implicit none

public :: ini3D_from_cavgs_commander
public :: het_ensemble_commander
private

type, extends(commander_base) :: ini3D_from_cavgs_commander
  contains
    procedure :: execute      => exec_ini3D_from_cavgs
end type ini3D_from_cavgs_commander
type, extends(commander_base) :: het_ensemble_commander
  contains
    procedure :: execute      => exec_het_ensemble
end type het_ensemble_commander

contains

    ! GENERATE INITIAL 3D MODEL FROM CLASS AVERAGES

    subroutine exec_ini3D_from_cavgs( self, cline )
        use simple_commander_volops,  only: projvol_commander
        use simple_commander_rec,     only: recvol_commander
        use simple_strings,           only: int2str_pad, str2int
        use simple_scaler,            only: scaler
        use simple_oris,              only: oris
        class(ini3D_from_cavgs_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! constants
        logical,               parameter :: DEBUG=.false.
        real,                  parameter :: LPLIMS(2)=[20.,10.] ! default low-pass limits
        real,                  parameter :: CENLP=30.           ! consistency with prime3D
        real,                  parameter :: SMPD_TARGET=3.3     ! 4 auto-scaling
        integer,               parameter :: MAXITS_INIT=30, MAXITS_REFINE=80
        integer,               parameter :: STATE=1, NPROJS_SYMSRCH=100
        character(len=32),     parameter :: ITERFBODY     = 'prime3Ddoc_'
        character(len=32),     parameter :: VOLFBODY      = 'recvol_state'
        character(len=STDLEN), parameter :: STKSCALEDBODY = 'stk_sc_ini3D_from_cavgs'
        ! distributed commanders
        type(prime3D_distr_commander) :: xprime3D_distr
        type(symsrch_distr_commander) :: xsymsrch_distr
        ! shared-mem commanders
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
        type(scaler)                  :: scobj
        type(qsys_env)                :: qenv
        type(params)                  :: p_master
        type(oris)                    :: os
        real                          :: iter, scale, smpd_sc, msk_sc, native_msk, native_smpd
        character(len=2)              :: str_state
        character(len=STDLEN)         :: oritab, vol_iter
        logical                       :: srch4symaxis
        integer                       :: box_sc
        ! set cline defaults
        call cline%set('eo', 'no')
        ! make master parameters
        p_master = params(cline, checkdistr=.false.)
        ! setup the environment for distributed execution
        call qenv%new(p_master)
        ! set global state string
        str_state = int2str_pad(STATE,2)
        ! delete possibly pre-existing stack_parts
        call del_files('stack_part', p_master%nparts, ext=p_master%ext)
        ! decide wether to search for the symmetry axis or put the point-group in from the start
        srch4symaxis = .false.
        if( p_master%pgrp .ne. 'c1' )then
            if(  p_master%pgrp(1:1).eq.'c'  .or. p_master%pgrp(1:1).eq.'C'&
            .or. p_master%pgrp(1:2).eq.'d2' .or. p_master%pgrp(1:2).eq.'D2' )then
                srch4symaxis = .true.
            endif
        endif
        ! auto-scaling prep
        call scobj%init(p_master, cline, SMPD_TARGET, STKSCALEDBODY)
        ! prepare command lines from prototype master
        cline_prime3D_init    = cline
        cline_prime3D_refine1 = cline
        cline_prime3D_refine2 = cline
        cline_symsrch         = cline
        cline_recvol          = cline
        cline_projvol         = cline
        ! initialise command line parameters
        ! (2) PRIME3D_INIT
        call cline_prime3D_init%set('prg',    'prime3D')
        call cline_prime3D_init%set('ctf',    'no')
        call cline_prime3D_init%set('maxits', real(MAXITS_INIT))
        call cline_prime3D_init%set('dynlp',  'yes') ! better be explicit about the dynlp
        ! (3) PRIME3D REFINE STEP 1
        call cline_prime3D_refine1%set('prg',    'prime3D')
        call cline_prime3D_refine1%set('ctf',    'no')
        call cline_prime3D_refine1%set('maxits', real(MAXITS_REFINE))
        call cline_prime3D_refine1%set('dynlp',  'no') ! better be explicit about the dynlp
        call cline_prime3D_refine1%set('refine', 'shc')
        ! (4) SYMMETRY AXIS SEARCH
        if( srch4symaxis )then
            call scobj%update_smpd_msk(cline_symsrch, 'scaled')
            call scobj%update_stk_smpd_msk(cline_recvol, 'scaled')
            ! need to replace original point-group flag with c1
            call cline_prime3D_init%set('pgrp', 'c1') 
            call cline_prime3D_refine1%set('pgrp', 'c1')
            ! symsrch
            call cline_symsrch%set('prg', 'symsrch')
            call cline_symsrch%delete('stk')  ! volumetric symsrch
            call cline_symsrch%set('nptcls',  real(NPROJS_SYMSRCH))
            call cline_symsrch%set('nspace',  real(NPROJS_SYMSRCH))
            call cline_symsrch%set('cenlp',   CENLP)
            call cline_symsrch%set('outfile', 'symdoc.txt')
            ! (4.5) RECONSTRUCT SYMMETRISED VOLUME
            call cline_recvol%set('prg', 'recvol')
            call cline_recvol%set('trs',  5.) ! to assure that shifts are being used
            call cline_recvol%set('ctf',  'no')
            call cline_recvol%set('oritab', 'symdoc.txt')
            ! 2nd refinement step now uses the symmetrised vol and doc
            call cline_prime3D_refine2%set('oritab', 'symdoc.txt')
            call cline_prime3D_refine2%set('vol1', 'rec_sym'//p_master%ext)
        endif
        ! (5) PRIME3D REFINE STEP 2
        call cline_prime3D_refine2%set('prg', 'prime3D')
        call cline_prime3D_refine2%set('ctf', 'no')
        call cline_prime3D_refine2%set('maxits', real(MAXITS_REFINE))
        call cline_prime3D_refine2%set('dynlp', 'no') ! better be explicit about the dynlp
        call cline_prime3D_refine2%set('lp', LPLIMS(2))
        call cline_prime3D_refine2%set('refine', 'shc')
        ! (6) RE-PROJECT VOLUME
        call cline_projvol%set('prg', 'projvol')
        call cline_projvol%set('outstk', 'reprojs'//p_master%ext)
        call cline_projvol%delete('stk')
        call scobj%update_smpd_msk(cline_projvol, 'native')
        ! scale class averages
        call scobj%scale_exec()
        ! execute commanders
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> INITIAL 3D MODEL GENERATION WITH PRIME3D'
        write(*,'(A)') '>>>'
        call xprime3D_distr%execute(cline_prime3D_init)
        iter = cline_prime3D_init%get_rarg('endit')
        call set_iter_dependencies
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> FIRST PRIME3D REFINEMENT STEP, REFINE=SHC'
        write(*,'(A)') '>>>'
        call cline_prime3D_refine1%set('startit', iter + 1.0)
        call cline_prime3D_refine1%set('oritab', trim(oritab))
        call cline_prime3D_refine1%set('vol1', trim(vol_iter))
        call update_lp(cline_prime3D_refine1, 1)
        call xprime3D_distr%execute(cline_prime3D_refine1)
        iter = cline_prime3D_refine1%get_rarg('endit')
        call set_iter_dependencies
        if( srch4symaxis )then
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> SYMMETRY AXIS SEARCH'
            write(*,'(A)') '>>>'
            call cline_symsrch%set('oritab', trim(oritab))
            call cline_symsrch%set('vol1', trim(vol_iter))
            call update_lp(cline_symsrch, 1)
            call xsymsrch_distr%execute(cline_symsrch)
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> 3D RECONSTRUCTION OF SYMMETRISED VOLUME'
            write(*,'(A)') '>>>'
            call xrecvol%execute(cline_recvol)
            call rename(trim(volfbody)//trim(str_state)//p_master%ext, 'rec_sym'//p_master%ext)
        else
            ! 2nd refinement step needs to use iter dependent vol/oritab
            call cline_prime3D_refine2%set('oritab', trim(oritab))
            call cline_prime3D_refine2%set('vol1', trim(vol_iter))
        endif
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> SECOND PRIME3D REFINEMENT STEP, REFINE=SHC'
        write(*,'(A)') '>>>'
        call cline_prime3D_refine2%set('startit', iter + 1.0)
        call update_lp(cline_prime3D_refine2, 2)
        call xprime3D_distr%execute(cline_prime3D_refine2)
        iter = cline_prime3D_refine2%get_rarg('endit')
        call set_iter_dependencies
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> 3D RECONSTRUCTION AT NATIVE SAMPLING'
        write(*,'(A)') '>>>'
        ! modulate shifts
        call os%new(p_master%nptcls)
        call os%read(oritab)
        call os%mul_shifts(1./scobj%get_scaled_var('scale'))
        call os%write(oritab)
        ! prepare recvol command line
        call scobj%update_stk_smpd_msk(cline_recvol, 'native')
        call cline_recvol%set('oritab', trim(oritab))
        ! re-reconstruct volume
        call xrecvol%execute(cline_recvol)
        call rename(trim(volfbody)//trim(str_state)//p_master%ext, 'rec_final'//p_master%ext)
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> RE-PROJECTION OF THE FINAL VOLUME'
        write(*,'(A)') '>>>'
        call cline_projvol%set('vol1', 'rec_final'//p_master%ext)
        call cline_projvol%set('oritab', trim(oritab))
        call xprojvol%execute(cline_projvol)
        ! end gracefully
        call del_file(trim(STKSCALEDBODY)//p_master%ext)
        call simple_end('**** SIMPLE_INI3D_FROM_CAVGS NORMAL STOP ****')

        contains

            subroutine set_iter_dependencies
                character(len=3) :: str_iter
                str_iter = int2str_pad(nint(iter),3)
                oritab   = trim(ITERFBODY)//trim(str_iter)//'.txt'
                vol_iter = trim(VOLFBODY)//trim(str_state)//'_iter'//trim(str_iter)//p_master%ext                
            end subroutine set_iter_dependencies

            subroutine update_lp( cl, refine_step )
                class(cmdline), intent(inout) :: cl
                integer,        intent(in)    :: refine_step
                real :: lpstop, lplim_new
                if( cline%defined('lp') )then
                    call cl%set('lp', p_master%lp)
                else
                    if( cline_prime3D_init%defined('lpstop') )then
                        lpstop = cline_prime3D_init%get_rarg('lpstop')
                        call cl%set('lp', min(LPLIMS(refine_step),lpstop))
                    else
                        call cl%set('lp', LPLIMS(refine_step))
                    endif
                endif
                lplim_new = cl%get_rarg('lp')
                write(*,'(A,1X,F6.2)') '>>> UPDATED LOW-PASS LIMIT TO:', lplim_new
                write(*,'(A)') '>>>'
            end subroutine update_lp

    end subroutine exec_ini3D_from_cavgs

    ! ENSEMBLE HETEROGEINITY ANALYSIS

    subroutine exec_het_ensemble( self, cline )
        use simple_commander_rec,     only: recvol_commander
        use simple_strings,           only: int2str_pad
        use simple_oris,              only: oris
        use simple_combinatorics,     only: diverse_labeling, shc_aggregation
        use simple_filehandling,      only: file_exists, del_file
        class(het_ensemble_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        ! constants
        logical,               parameter :: DEBUG=.false.
        integer,               parameter :: MAXITS_INIT=50, NREPEATS=5
        character(len=32),     parameter :: HETFBODY    = 'hetrep_'
        character(len=32),     parameter :: REPEATFBODY = 'hetdoc_'
        character(len=32),     parameter :: VOLFBODY    = 'recvol_state'
        ! distributed commanders
        type(prime3D_distr_commander) :: xprime3D_distr
        type(recvol_distr_commander)  :: xrecvol_distr        
        ! shared-mem commanders
        !
        ! command lines
        type(cmdline)                 :: cline_prime3D_master
        type(cmdline)                 :: cline_prime3D
        type(cmdline)                 :: cline_recvol_distr
        ! other variables
        type(params)                  :: p_master
        type(oris)                    :: os, rep_os
        integer, allocatable          :: labels(:,:), labels_incl(:,:), consensus(:)
        logical, allocatable          :: included(:)
        character(len=STDLEN)         :: oritab, vol1, vol2, fname
        integer                       :: irepeat, state, iter, n_incl, it
        ! set cline defaults
        call cline%set('eo', 'no')
        if(nint(cline%get_rarg('nstates')) <= 1)stop 'Non-sensical NSTATES argument for heterogeinity analysis!'

        ! make master parameters
        p_master = params(cline, checkdistr=.false.)

        ! delete possibly pre-existing stack & pft parts
        call del_files('stack_part', p_master%nparts, ext=p_master%ext)
        call del_files('ppfts_memoized_part', p_master%nparts, ext='.bin')

        ! prepare command lines from prototype master
        cline_recvol_distr = cline
        call cline_recvol_distr%set('prg', 'recvol')
        call cline_recvol_distr%set('eo', 'yes')
        call cline_recvol_distr%delete('lp')
        cline_prime3D_master = cline
        call cline_prime3D_master%set('prg', 'prime3D')
        call cline_prime3D_master%set('startit', 1.)
        call cline_prime3D_master%set('maxits', real(MAXITS_INIT))
        call cline_prime3D_master%set('refine', 'het')
        call cline_prime3D_master%set('dynlp', 'no')
        call cline_prime3D_master%set('lp', p_master%lp) 

        ! GENERATE DIVERSE INITIAL LABELS
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> GENERATING DIVERSE LABELING'
        write(*,'(A)') '>>>'
        call os%new(p_master%nptcls)
        call os%read(p_master%oritab)
        labels   = diverse_labeling( p_master%nptcls, p_master%nstates, NREPEATS)
        included = os%included()
        n_incl   = count(included)
        do irepeat=1,NREPEATS
            where( .not.included )labels(irepeat,:) = 0
        enddo

        ! GENERATE CANDIDATE SOLUTIONS
        do irepeat = 1,NREPEATS
            write(*,'(A)')    '>>>'
            write(*,'(A,I3)') '>>> PRIME3D REPEAT ', irepeat
            write(*,'(A)')    '>>>'
            ! GENERATE ORIENTATIONS
            rep_os = os
            call rep_os%set_all('state', real(labels(irepeat,:)))
            oritab = trim(REPEATFBODY)//'init_rep'//int2str_pad(irepeat,2)//'.txt'
            call rep_os%write(trim(oritab))
            ! RUN PRIME3D
            cline_prime3D = cline_prime3D_master
            call cline_prime3D%set('oritab',oritab)
            call xprime3D_distr%execute(cline_prime3D)
            ! HARVEST OUTCOME
            iter   = nint(cline_prime3D%get_rarg('endit'))
            oritab = 'prime3Ddoc_'//int2str_pad(iter,3)//'.txt'
            call rep_os%read(trim(oritab))
            ! updates labels
            labels(irepeat,:) = nint(rep_os%get_all('state'))
            ! stash candidate solution
            fname = trim(REPEATFBODY)//'rep'//int2str_pad(irepeat,2)//'.txt'
            call rename(trim(oritab), trim(fname))
            ! STASH FINAL VOLUMES
            call stash_volumes
            ! CLEANUP
            call prime3d_cleanup
        enddo

        ! GENERATE CONSENSUS DOCUMENT
        oritab = trim(REPEATFBODY)//'consensus.txt'
        write(*,'(A)')   '>>>'
        write(*,'(A,A)') '>>> GENERATING ENSEMBLE SOLUTION: ', trim(oritab)
        write(*,'(A)')   '>>>'
        allocate(labels_incl(NREPEATS,n_incl), consensus(n_incl))
        do irepeat=1,NREPEATS
            labels_incl(irepeat,:) = pack(labels(irepeat,:), mask=included)
        enddo
        labels(1,:) = 0
        call shc_aggregation(NREPEATS, n_incl, labels_incl, consensus)
        call os%set_all('state', real(unpack(consensus, included, labels(1,:))) )
        call os%write(trim(oritab))
        call os%kill

        ! FINAL RECONSTRUCTION
        call cline_recvol_distr%set('oritab', trim(oritab))
        call xrecvol_distr%execute(cline_recvol_distr)

        ! end gracefully
        call simple_end('**** SIMPLE_HET_ENSEMBLE NORMAL STOP ****')        
        contains

            subroutine prime3d_cleanup
                ! delete starting volumes
                do state=1,p_master%nstates
                    vol1 = 'startvol_state'//int2str_pad(state,2)//p_master%ext
                    if(file_exists(vol1))call del_file(vol1)
                enddo
                ! delete iterations volumes & documents
                do it=1,iter-1
                    do state=1,p_master%nstates
                        vol1 = trim(VOLFBODY)//int2str_pad(state,2)//'_iter'//int2str_pad(it,3)//p_master%ext
                        if(file_exists(vol1))call del_file(vol1)
                        vol1 = trim(VOLFBODY)//int2str_pad(state,2)//'_iter'//int2str_pad(it,3)//'pproc'//p_master%ext
                        if(file_exists(vol1))call del_file(vol1)
                    enddo
                    oritab = 'prime3Ddoc_'//int2str_pad(it,3)//'.txt'
                    if(file_exists(oritab))call del_file(oritab)
                enddo
                ! delete restart documents
                do it=1,iter
                    fname = 'prime3D_restart_iter'//int2str_pad(it,3)//'.txt'
                    if(file_exists(fname))call del_file(fname)
                enddo
            end subroutine prime3d_cleanup

            subroutine stash_volumes
                ! renames final volumes of each repeat for all states
                do state=1,p_master%nstates
                    vol1 = trim(VOLFBODY)//int2str_pad(state,2)//'_iter'//int2str_pad(iter,3)//p_master%ext
                    vol2 = trim(HETFBODY)//int2str_pad(irepeat,2)//'_'//trim(VOLFBODY)//int2str_pad(state,2)//p_master%ext
                    call rename(trim(vol1), trim(vol2))
                    vol1 = trim(VOLFBODY)//int2str_pad(state,2)//'_iter'//int2str_pad(iter,3)//'pproc'//p_master%ext
                    vol2 = trim(HETFBODY)//int2str_pad(irepeat,2)//'_'//trim(VOLFBODY)//int2str_pad(state,2)//'pproc'//p_master%ext
                    call rename(trim(vol1), trim(vol2))
                enddo
            end subroutine stash_volumes

    end subroutine exec_het_ensemble

end module simple_commander_hlev_wflows
