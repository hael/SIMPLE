! concrete commander: high-level workflows
module simple_commander_hlev_wflows
#include "simple_lib.f08"
use simple_defs_fname
use simple_cmdline,               only: cmdline
use simple_params,                only: params
use simple_commander_base,        only: commander_base
use simple_qsys_env,              only: qsys_env
use simple_oris,                  only: oris
use simple_scaler,                only: scaler
use simple_strings,               only: int2str_pad, str2int
use simple_binoris_io,            only: binread_oritab, binwrite_oritab
use simple_commander_distr_wflows ! use all in there
use simple_commander_distr        ! use all in there
implicit none

public :: prime2D_autoscale_commander
public :: ini3D_from_cavgs_commander
public :: auto_refine3D_commander
public :: het_commander
public :: het_refine_commander
private

type, extends(commander_base) :: prime2D_autoscale_commander
  contains
    procedure :: execute      => exec_prime2D_autoscale
end type prime2D_autoscale_commander
type, extends(commander_base) :: ini3D_from_cavgs_commander
  contains
    procedure :: execute      => exec_ini3D_from_cavgs
end type ini3D_from_cavgs_commander
type, extends(commander_base) :: auto_refine3D_commander
  contains
    procedure :: execute      => exec_auto_refine3D
end type auto_refine3D_commander
type, extends(commander_base) :: het_commander
  contains
    procedure :: execute      => exec_het
end type het_commander
type, extends(commander_base) :: het_refine_commander
  contains
    procedure :: execute      => exec_het_refine
end type het_refine_commander

contains

    !> for distributed PRIME2D with two-stage autoscaling
    subroutine exec_prime2D_autoscale( self, cline )
        use simple_commander_prime2D, only: rank_cavgs_commander
        class(prime2D_autoscale_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        ! constants
        integer,           parameter :: MAXITS_STAGE1   = 10
        ! commanders
        type(split_commander)           :: xsplit
        type(makecavgs_distr_commander) :: xmakecavgs
        type(prime2D_distr_commander)   :: xprime2D_distr
        type(rank_cavgs_commander)      :: xrank_cavgs
        ! command lines
        type(cmdline) :: cline_prime2D_stage1
        type(cmdline) :: cline_prime2D_stage2
        type(cmdline) :: cline_makecavgs
        type(cmdline) :: cline_rank_cavgs
        ! other variables
        character(len=STDLEN) :: scaled_stktab
        character(len=STDLEN) :: finaldoc, finalcavgs, finalcavgs_ranked
        type(oris)            :: os
        type(scaler)          :: scobj
        type(params)          :: p_master
        real                  :: scale_stage1, scale_stage2
        integer               :: nparts, last_iter_stage1, last_iter_stage2
        p_master = params(cline, del_scaled=.true.)
        nparts   = p_master%nparts
        if( .not. cline%defined('stktab') )then
            ! split stack
            call xsplit%execute(cline)
        endif
        if( p_master%l_autoscale )then
            ! auto-scaling prep (cline is modified by scobj%init)
            call scobj%init(p_master, cline, p_master%box, p_master%smpd_targets2D(1))
            scale_stage1 = scobj%get_scaled_var('scale')
            ! scale images in parallel
            call scobj%scale_distr_exec
            if( cline%defined('stktab') )then
                ! updates command lines
                scaled_stktab = add2fbody(p_master%stktab, METADATEXT, SCALE_SUFFIX)
                ! update stktab
                call p_master%stkhandle%add_scale_tag
                call p_master%stkhandle%write_stktab(trim(scaled_stktab))
            endif
            ! execute stage 1
            cline_prime2D_stage1 = cline
            if( cline%defined('stktab') )then
                call cline_prime2D_stage1%set('stktab', trim(scaled_stktab))
            endif
            call cline_prime2D_stage1%delete('automsk') ! delete possible automsk flag from stage 1
            call cline_prime2D_stage1%set('maxits', real(MAXITS_STAGE1))
            call xprime2D_distr%execute(cline_prime2D_stage1)
            last_iter_stage1 = nint(cline_prime2D_stage1%get_rarg('endit'))
            finaldoc         = trim(PRIME2D_ITER_FBODY)//int2str_pad(last_iter_stage1,3)//METADATEXT
            finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage1,3)//p_master%ext
            ! prepare stage 2 input -- re-scale
            call scobj%uninit(cline) ! puts back the old command line
            call scobj%init(p_master, cline, p_master%box, p_master%smpd_targets2D(2))
            scale_stage2 = scobj%get_scaled_var('scale')
            call scobj%scale_distr_exec
            ! prepare stage 2 input -- shift modulation
            call os%new(p_master%nptcls)
            call binread_oritab(finaldoc, os, [1,p_master%nptcls])
            call os%mul_shifts(scale_stage2/scale_stage1)
            call binwrite_oritab(finaldoc, os, [1,p_master%nptcls])
            ! prepare stage 2 input -- command line
            cline_prime2D_stage2 = cline
            ! if automsk .eq. yes, we need to replace it with cavg
            if( p_master%automsk .eq. 'yes' )then
                call cline_prime2D_stage2%set('automsk', 'cavg')
            endif
            if( cline%defined('stktab') )then
                call cline_prime2D_stage2%set('stktab', trim(scaled_stktab))
            endif
            call cline_prime2D_stage2%delete('deftab')
            call cline_prime2D_stage2%set('oritab',  trim(finaldoc))
            call cline_prime2D_stage2%set('startit', real(MAXITS_STAGE1 + 1))
            call xprime2D_distr%execute(cline_prime2D_stage2)
            last_iter_stage2 = nint(cline_prime2D_stage2%get_rarg('endit'))
            finaldoc         = trim(PRIME2D_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//METADATEXT
            finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//p_master%ext
            if( cline%defined('stktab') )then
                ! delete downscaled stack parts & stktab (we are done with them)
                call del_file(trim(scaled_stktab))
                call p_master%stkhandle%del_stktab_files
                ! put back original stktab
                call p_master%stkhandle%del_scale_tag
            else
                ! delete downscaled stack parts (we are done with them)
                call del_files(trim(STKPARTFBODY), p_master%nparts, ext=p_master%ext, suffix='_sc')
            endif
            ! re-generate class averages at original sampling
            call scobj%uninit(cline) ! puts back the old command line
            call binread_oritab(finaldoc, os, [1,p_master%nptcls])
            call os%mul_shifts(1./scale_stage2)
            call binwrite_oritab(finaldoc, os, [1,p_master%nptcls])
            cline_makecavgs = cline
            call cline_makecavgs%delete('autoscale')
            call cline_makecavgs%delete('balance')
            call cline_makecavgs%delete('chunksz')
            if( p_master%l_chunk_distr )then
                call cline_makecavgs%delete('ncls')
            endif
            call cline_makecavgs%set('prg',    'makecavgs')
            call cline_makecavgs%set('oritab', trim(finaldoc))
            call cline_makecavgs%set('nparts', real(nparts))
            call cline_makecavgs%set('refs',   trim(finalcavgs))
            call xmakecavgs%execute(cline_makecavgs)
        else ! no auto-scaling
            call xprime2D_distr%execute(cline)
            last_iter_stage2 = nint(cline%get_rarg('endit'))
            finaldoc         = trim(PRIME2D_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//METADATEXT
            finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//p_master%ext
        endif
        ! ranking
        finalcavgs_ranked = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//'_ranked'//p_master%ext
        call cline_rank_cavgs%set('oritab',   trim(finaldoc))
        call cline_rank_cavgs%set('stk',      trim(finalcavgs))
        call cline_rank_cavgs%set('classdoc', 'classdoc_'//int2str_pad(last_iter_stage2,3)//'.txt')
        call cline_rank_cavgs%set('outstk',   trim(finalcavgs_ranked))
        call xrank_cavgs%execute( cline_rank_cavgs )
        ! cleanup
        call del_file('prime2D_startdoc'//METADATEXT)
        call del_file('start2Drefs'//p_master%ext)
        ! end gracefully
        call simple_end('**** SIMPLE_PRIME2D NORMAL STOP ****')
    end subroutine exec_prime2D_autoscale

    !> for generation of an initial 3d model from class averages
    subroutine exec_auto_refine3D( self, cline )
        use simple_defs_conv
        !use simple_commander_rec, only: recvol_commander
        class(auto_refine3D_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        ! constants
        integer,               parameter :: MAXITS_REFINE1=20, MAXITS_REFINE2=20
        integer,               parameter :: STATE=1
        character(len=32),     parameter :: ITERFBODY     = 'prime3Ddoc_'
        ! distributed commanders
        type(prime3D_distr_commander) :: xprime3D_distr
        ! shared-mem commanders
        !type(recvol_commander)        :: xrecvol
        ! command lines
        type(cmdline)         :: cline_prime3D_1, cline_prime3D_2, cline_prime3D_3
        type(cmdline)         :: cline_recvol
        ! other variables
        type(params)          :: p_master
        real                  :: iter, trs_lim, greedytrs_lim
        character(len=2)      :: str_state
        character(len=3)      :: str_iter
        character(len=STDLEN) :: vol_iter, oritab, vol_pproc
        ! set cline defaults
        call cline%set('eo', 'yes')
        call cline%set('dynlp', 'no') ! better be explicit about the dynlp
        ! make master parameters
        p_master = params(cline)
        ! delete possibly pre-existing stack_parts
        call del_files(STKPARTFBODY, p_master%nparts, ext=p_master%ext)
        ! set global state string
        str_state = int2str_pad(STATE,2)
        ! prepare command lines from prototype master
        cline_prime3D_1 = cline
        cline_prime3D_2 = cline
        cline_prime3D_3 = cline
        ! determines shift limits
        if( cline%defined('trs') )then
            ! all good
        else
            trs_lim = MSK_FRAC*real(p_master%msk)
            trs_lim = max(MINSHIFT, trs_lim)
            trs_lim = min(MAXSHIFT, trs_lim)
        endif
        greedytrs_lim = min(MAXSHIFT, 2.*trs_lim)
        ! REFINEMENT FROM MODEL - STEP 1
        call cline_prime3D_1%set('prg', 'prime3D')
        call cline_prime3D_1%set('refine', 'greedy')
        call cline_prime3D_1%set('trs', greedytrs_lim)
        ! REFINEMENT FROM MODEL - STEP 2
        call cline_prime3D_2%set('prg', 'prime3D')
        call cline_prime3D_2%set('refine', 'no')
        call cline_prime3D_2%set('trs', trs_lim)
        ! REFINEMENT FROM MODEL - STEP 3
        call cline_prime3D_3%set('prg', 'prime3D')
        call cline_prime3D_3%set('refine', 'greedyneigh')
        call cline_prime3D_3%set('trs', trs_lim)
        ! execute commanders
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> PRIME3D REFINEMENT STEP 1'
        write(*,'(A)') '>>>'
        call cline_prime3D_1%set('startit', 1.)
        call cline_prime3D_1%set('maxits', 1.)
        call xprime3D_distr%execute(cline_prime3D_1)
        iter = cline_prime3D_1%get_rarg('endit')
        call set_iter_dependencies
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> PRIME3D REFINEMENT STEP 2'
        write(*,'(A)') '>>>'
        iter = iter + 1
        call cline_prime3D_2%set('startit', real(iter))
        call cline_prime3D_2%set('maxits', real(MAXITS_REFINE1)+iter-1)
        call cline_prime3D_2%set('oritab', oritab)
        call cline_prime3D_2%set('vol1', vol_iter)
        call xprime3D_distr%execute(cline_prime3D_2)
        iter = cline_prime3D_2%get_rarg('endit')
        call set_iter_dependencies
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> PRIME3D REFINEMENT STEP 3'
        write(*,'(A)') '>>>'
        iter = iter + 1
        call cline_prime3D_2%set('maxits', real(MAXITS_REFINE2)+iter-1)
        call cline_prime3D_3%set('startit', real(iter))
        call cline_prime3D_3%set('oritab', oritab)
        call cline_prime3D_3%set('vol1', vol_iter)
        call xprime3D_distr%execute(cline_prime3D_3)
        iter = cline_prime3D_3%get_rarg('endit')
        call set_iter_dependencies
        vol_pproc = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//'_pproc'//p_master%ext
        call simple_rename(trim(vol_iter), 'rec_final'//p_master%ext)
        call simple_rename(trim(vol_pproc), 'rec_final_pproc'//p_master%ext)
        call simple_rename(trim(oritab), 'prime3Ddoc_final'//METADATEXT)
        ! delete stack parts (we are done with them)
        call del_files(trim(STKPARTFBODY), p_master%nparts, ext=p_master%ext)
        ! end gracefully
        call simple_end('**** SIMPLE_AUTO_REFINE3D NORMAL STOP ****')

        contains

            subroutine set_iter_dependencies
                str_iter = int2str_pad(nint(iter),3)
                oritab   = trim(ITERFBODY)//trim(str_iter)//METADATEXT
                vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//p_master%ext
            end subroutine set_iter_dependencies

    end subroutine exec_auto_refine3D

    !> for generation of an initial 3d model from class averages
    subroutine exec_ini3D_from_cavgs( self, cline )
        use simple_commander_volops, only: projvol_commander
        use simple_commander_rec,    only: recvol_commander
        class(ini3D_from_cavgs_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! constants
        real,                  parameter :: CENLP=30.           !< consistency with prime3D
        integer,               parameter :: MAXITS_SNHC=30, MAXITS_INIT=15, MAXITS_REFINE=40
        integer,               parameter :: STATE=1, NPROJS_SYMSRCH=50, NPEAKS_REFINE=6
        integer,               parameter :: NSPACE_SNHC = 1000, NSPACE_DEFAULT= 2500
        character(len=32),     parameter :: ITERFBODY     = 'prime3Ddoc_'
        character(len=STDLEN), parameter :: STKSCALEDBODY = 'stk_sc_ini3D_from_cavgs'
        ! distributed commanders
        type(prime3D_distr_commander) :: xprime3D_distr
        type(symsrch_distr_commander) :: xsymsrch_distr
        ! shared-mem commanders
        type(recvol_commander)        :: xrecvol
        type(projvol_commander)       :: xprojvol
        ! command lines
        type(cmdline)         :: cline_prime3D_snhc
        type(cmdline)         :: cline_prime3D_init
        type(cmdline)         :: cline_prime3D_refine
        type(cmdline)         :: cline_symsrch
        type(cmdline)         :: cline_recvol
        type(cmdline)         :: cline_projvol
        ! other variables
        type(scaler)          :: scobj
        type(params)          :: p_master
        type(oris)            :: os
        real                  :: iter, smpd_target, lplims(2)
        character(len=2)      :: str_state
        character(len=STDLEN) :: vol_iter, oritab
        logical               :: srch4symaxis, doautoscale
        ! set cline defaults
        call cline%set('eo', 'no')
        ! auto-scaling prep
        doautoscale = (cline%get_carg('autoscale').eq.'yes')
        ! now, remove autoscale flag from command line, since no scaled partial stacks 
        ! will be produced (this program used shared-mem paralllelisation of scale)
        call cline%delete('autoscale')
        ! delete possibly pre-existing stack_parts
        call del_files(STKPARTFBODY, p_master%nparts, ext=p_master%ext)
        call del_files(STKPARTFBODY, p_master%nparts, ext=p_master%ext, suffix='_sc')
        ! make master parameters
        p_master = params(cline)
        ! set global state string
        str_state = int2str_pad(STATE,2)
        ! decide wether to search for the symmetry axis or put the point-group in from the start
        ! if the point-group is considered known, it is put in from the start
        srch4symaxis = .false.
        if( p_master%pgrp_known .eq. 'no' )then
            if( p_master%pgrp .ne. 'c1' )then
                if(  p_master%pgrp(1:1).eq.'c'  .or. p_master%pgrp(1:1).eq.'C'&
                .or. p_master%pgrp(1:2).eq.'d2' .or. p_master%pgrp(1:2).eq.'D2' )then
                    srch4symaxis = .true.
                endif
            endif
        endif
        ! set lplims
        ! default
        lplims(1) = 20.
        lplims(2) = 10.
        ! passed
        if( cline%defined('lpstart') ) lplims(1) = p_master%lpstart
        if( cline%defined('lpstop')  ) lplims(2) = p_master%lpstop
        ! init scaler
        smpd_target = p_master%smpd
        if( doautoscale )then
            smpd_target = lplims(2)*LP2SMPDFAC
            call scobj%init(p_master, cline, p_master%box, smpd_target, STKSCALEDBODY)
        endif
        ! prepare command lines from prototype master
        cline_prime3D_snhc   = cline
        cline_prime3D_init   = cline
        cline_prime3D_refine = cline
        cline_symsrch        = cline
        cline_recvol         = cline
        cline_projvol        = cline
        ! recvol & projvol are not distributed executions, so remove the nparts flag
        call cline_recvol%delete('nparts')
        call cline_projvol%delete('nparts')
        ! initialise command line parameters
        ! (1) INITIALIZATION BY STOCHASTIC NEIGHBORHOOD HILL-CLIMBING
        call cline_prime3D_snhc%delete('update_frac') ! no fractional update in first phase
        call cline_prime3D_snhc%set('prg',    'prime3D')
        call cline_prime3D_snhc%set('ctf',    'no')
        call cline_prime3D_snhc%set('maxits', real(MAXITS_SNHC))
        call cline_prime3D_snhc%set('refine', 'snhc')
        call cline_prime3D_snhc%set('dynlp',  'no') ! better be explicit about the dynlp
        call cline_prime3D_snhc%set('lp',     lplims(1))
        call cline_prime3D_snhc%set('nspace', real(NSPACE_SNHC))
        ! (2) PRIME3D_INIT
        call cline_prime3D_init%set('prg',    'prime3D')
        call cline_prime3D_init%set('ctf',    'no')
        call cline_prime3D_init%set('maxits', real(MAXITS_INIT))
        call cline_prime3D_init%set('vol1',   trim(SNHCVOL)//trim(str_state)//p_master%ext)
        call cline_prime3D_init%set('oritab', SNHCDOC)
        call cline_prime3D_init%set('dynlp',  'no') ! better be explicit about the dynlp
        call cline_prime3D_init%set('lp',     lplims(1))
        if( .not. cline_prime3D_init%defined('nspace') )then
            call cline_prime3D_init%set('nspace', real(NSPACE_DEFAULT))
        endif
        ! (3) SYMMETRY AXIS SEARCH
        if( srch4symaxis )then
            ! need to replace original point-group flag with c1
            call cline_prime3D_snhc%set('pgrp', 'c1')
            call cline_prime3D_init%set('pgrp', 'c1')
            ! symsrch
            call cline_symsrch%set('prg', 'symsrch')
            call cline_symsrch%delete('stk')  ! volumetric symsrch
            call cline_symsrch%set('nptcls',  real(NPROJS_SYMSRCH))
            call cline_symsrch%set('nspace',  real(NPROJS_SYMSRCH))
            call cline_symsrch%set('cenlp',   CENLP)
            call cline_symsrch%set('outfile', 'symdoc'//METADATEXT)
            call cline_symsrch%set('lp',      lplims(2))
            ! (4.5) RECONSTRUCT SYMMETRISED VOLUME
            call cline_recvol%set('prg', 'recvol')
            call cline_recvol%set('trs',  5.) ! to assure that shifts are being used
            call cline_recvol%set('ctf',  'no')
            call cline_recvol%set('oritab', 'symdoc'//METADATEXT)
            ! refinement step now uses the symmetrised vol and doc
            call cline_prime3D_refine%set('oritab', 'symdoc'//METADATEXT)
            call cline_prime3D_refine%set('vol1', 'rec_sym'//p_master%ext)
        endif
        ! (4) PRIME3D REFINE STEP
        call cline_prime3D_refine%set('prg', 'prime3D')
        call cline_prime3D_refine%set('ctf', 'no')
        call cline_prime3D_refine%set('maxits', real(MAXITS_REFINE))
        call cline_prime3D_refine%set('refine', 'no')
        call cline_prime3D_refine%set('npeaks', real(NPEAKS_REFINE))
        call cline_prime3D_refine%set('dynlp', 'no') ! better be explicit about the dynlp
        call cline_prime3D_refine%set('lp', lplims(2))
        if( .not. cline_prime3D_refine%defined('nspace') )then
            call cline_prime3D_refine%set('nspace', real(NSPACE_DEFAULT))
        endif
        ! (5) RE-PROJECT VOLUME
        call cline_projvol%set('prg', 'projvol')
        call cline_projvol%set('outstk', 'reprojs'//p_master%ext)
        call cline_projvol%delete('stk')
        if( doautoscale )then
            call scobj%update_smpd_msk(cline_projvol, 'original')
            ! scale class averages
            call scobj%scale_exec
        endif
        ! execute commanders
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> INITIALIZATION WITH STOCHASTIC NEIGHBORHOOD HILL-CLIMBING'
        write(*,'(A)') '>>>'
        call xprime3D_distr%execute(cline_prime3D_snhc)
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> INITIAL 3D MODEL GENERATION WITH PRIME3D'
        write(*,'(A)') '>>>'
        call xprime3D_distr%execute(cline_prime3D_init)
        iter = cline_prime3D_init%get_rarg('endit')
        call set_iter_dependencies
        if( srch4symaxis )then
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> SYMMETRY AXIS SEARCH'
            write(*,'(A)') '>>>'
            call cline_symsrch%set('oritab', trim(oritab))
            call cline_symsrch%set('vol1', trim(vol_iter))
            call xsymsrch_distr%execute(cline_symsrch)
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> 3D RECONSTRUCTION OF SYMMETRISED VOLUME'
            write(*,'(A)') '>>>'
            call xrecvol%execute(cline_recvol)
            call simple_rename(trim(VOL_FBODY)//trim(str_state)//p_master%ext, 'rec_sym'//p_master%ext)
        else
            ! refinement step needs to use iter dependent vol/oritab
            call cline_prime3D_refine%set('oritab', trim(oritab))
            call cline_prime3D_refine%set('vol1', trim(vol_iter))
        endif
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> PRIME3D REFINEMENT STEP'
        write(*,'(A)') '>>>'
        call cline_prime3D_refine%set('startit', iter + 1.0)
        call xprime3D_distr%execute(cline_prime3D_refine)
        iter = cline_prime3D_refine%get_rarg('endit')
        call set_iter_dependencies
        ! delete stack parts (we are done with them)
        call del_files(trim(STKPARTFBODY), p_master%nparts, ext=p_master%ext)
        if( doautoscale )then
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> 3D RECONSTRUCTION AT ORIGINAL SAMPLING'
            write(*,'(A)') '>>>'
            ! modulate shifts
            call os%new(p_master%nptcls)
            call binread_oritab(oritab, os, [1,p_master%nptcls])
            call os%mul_shifts(1./scobj%get_scaled_var('scale'))
            call binwrite_oritab(oritab, os, [1,p_master%nptcls])
            ! prepare recvol command line
            call scobj%update_stk_smpd_msk(cline_recvol, 'original')
            call cline_recvol%set('oritab', trim(oritab))
            ! re-reconstruct volume
            call xrecvol%execute(cline_recvol)
            call simple_rename(trim(VOL_FBODY)//trim(str_state)//p_master%ext, 'rec_final'//p_master%ext)
        else
            call simple_rename(trim(vol_iter), 'rec_final'//p_master%ext)
        endif
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
                oritab   = trim(ITERFBODY)//trim(str_iter)//METADATEXT
                vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//p_master%ext
            end subroutine set_iter_dependencies

    end subroutine exec_ini3D_from_cavgs

    !> for heterogeinity analysis
    subroutine exec_het( self, cline )
        use simple_defs_conv
        use simple_commander_rec,    only: recvol_commander
        use simple_commander_volops, only: postproc_vol_commander
        use simple_sym,              only: sym
        class(het_commander), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        ! constants
        integer,            parameter :: MAXITS   = 50
        ! distributed commanders
        type(prime3D_distr_commander) :: xprime3D_distr
        type(recvol_distr_commander)  :: xrecvol_distr
        ! shared-mem commanders
        type(postproc_vol_commander)  :: xpostproc_vol
        ! command lines
        type(cmdline)                 :: cline_prime3D
        type(cmdline)                 :: cline_recvol_distr
        type(cmdline)                 :: cline_postproc_vol
        ! other variables
        type(sym)                     :: symop
        integer,     allocatable      :: labels(:)
        logical,     allocatable      :: included(:)
        type(params)                  :: p_master
        type(oris)                    :: os, os_states
        character(len=STDLEN)         :: oritab, vol, str_state
        real                          :: trs
        integer                       :: state, iter, n_incl
        call seed_rnd
        ! sanity check
        if(nint(cline%get_rarg('nstates')) <= 1)&
            &stop 'Non-sensical NSTATES argument for heterogeneity analysis!'
        ! make master parameters
        p_master = params(cline)
        if( p_master%eo .eq. 'no' .and. .not. cline%defined('lp') )&
            &stop 'need lp input when eo .eq. no; het'
        ! prepare command lines from prototype
        call cline%delete('refine')
        cline_prime3D      = cline
        cline_postproc_vol = cline ! eo always eq yes
        cline_recvol_distr = cline ! eo always eq yes
        call cline_prime3D%set('prg', 'prime3D')
        call cline_prime3D%set('maxits', real(MAXITS))
        if( trim(p_master%refine) .eq. 'sym' )then
            call cline_prime3D%set('refine', 'hetsym')
        else
            call cline_prime3D%set('refine', 'het')
        endif
        call cline_prime3D%set('dynlp', 'no')
        call cline_prime3D%set('pproc', 'no')
        call cline_prime3D%delete('oritab2')
        call cline_recvol_distr%set('prg', 'recvol')
        call cline_postproc_vol%set('prg', 'postproc_vol')
        call cline_recvol_distr%delete('lp')
        call cline_postproc_vol%delete('lp')
        call cline_recvol_distr%set('eo','yes')
        call cline_postproc_vol%set('eo','yes')
        ! works out shift lmits for in-plane search
        if( cline%defined('trs') )then
            ! all good
        else
            trs = MSK_FRAC*real(p_master%msk)
            trs = min(MAXSHIFT, max(MINSHIFT, trs))
            call cline_prime3D%set('trs',trs)
        endif
        ! generate diverse initial labels
        oritab = 'hetdoc_init'//METADATEXT
        call os%new(p_master%nptcls)
        call binread_oritab(p_master%oritab, os, [1,p_master%nptcls])
        if( p_master%eo.eq.'no' )then
            ! updates e/o flags
            call os%set_all2single('eo', -1.)
        else
            if( os%get_nevenodd() == 0 ) call os%partition_eo
        endif
        if( trim(p_master%refine) .eq. 'sym' )then
            ! randomize projection directions with respect to symmetry
            symop = sym(p_master%pgrp)
            call symop%symrandomize(os)
            call symop%kill
        endif
        if( cline%defined('oritab2') )then
            ! this is  to force initialisation (4 testing)
            call os_states%new(p_master%nptcls)
            call binread_oritab(p_master%oritab2, os_states, [1,p_master%nptcls])
            call os%set_all('state', os_states%get_all('state'))
            call os_states%kill
        else if( .not. cline%defined('startit') )then
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> GENERATING DIVERSE LABELING'
            call diverse_labeling(os, p_master%nstates, labels, corr_ranked=.true.)
            call os%set_all('state', real(labels))
        else
            ! starting from a previous solution
            labels = os%get_all('state')
        endif
        ! to accomodate state=0s in oritab input
        included = os%included()
        n_incl   = count(included)
        where( .not. included ) labels = 0
        call os%set_all('state', real(labels))
        call binwrite_oritab(trim(oritab), os, [1,p_master%nptcls])
        deallocate(labels, included)
        call cline_prime3D%set('oritab', trim(oritab))
        ! run prime3d
        write(*,'(A)')    '>>>'
        write(*,'(A,I3)') '>>> PRIME3D'
        write(*,'(A)')    '>>>'
        call xprime3D_distr%execute(cline_prime3D)
        ! final distributed reconstruction to obtain resolution estimate when eo .eq. 'no'
        if( p_master%eo .eq. 'no' )then
            iter   = nint(cline_prime3D%get_rarg('endit'))
            oritab = PRIME3D_ITER_FBODY//int2str_pad(iter,3)//METADATEXT
            call cline_recvol_distr%set('oritab', trim(oritab))
            call xrecvol_distr%execute(cline_recvol_distr)
            do state = 1, p_master%nstates
                str_state  = int2str_pad(state, 2)
                call cline_postproc_vol%set('vol1', VOL_FBODY//trim(str_state)//p_master%ext)
                call cline_postproc_vol%set('fsc', FSC_FBODY//trim(str_state)//BIN_EXT)
                call cline_postproc_vol%set('vol_filt', ANISOLP_FBODY//trim(str_state)//p_master%ext)
                call xpostproc_vol%execute(cline_postproc_vol)
            enddo
        endif
        ! kill oris
        call os%kill
        ! end gracefully
        call simple_end('**** SIMPLE_HET NORMAL STOP ****')

        contains

            subroutine diverse_labeling( os, nlabels, config_diverse, corr_ranked )
                use simple_oris, only: oris
                type(oris),           intent(inout) :: os
                integer,              intent(in)    :: nlabels
                integer, allocatable, intent(out)   :: config_diverse(:)
                logical, optional,    intent(in)    :: corr_ranked
                integer, allocatable :: tmp(:), states(:), order(:)
                real,    allocatable :: corrs(:)
                type(ran_tabu)       :: rt
                integer              :: iptcl, nptcls,nonzero_nptcls, alloc_stat, cnt, s, ind
                logical              :: corr_ranked_here
                corr_ranked_here = .false.
                if(present(corr_ranked)) corr_ranked_here = corr_ranked
                if(.not.os%isthere('corr')) corr_ranked_here = .false.
                nptcls = os%get_noris()
                states = nint(os%get_all('state'))
                nonzero_nptcls = count(states>0)
                if( .not.corr_ranked_here )then
                    allocate(config_diverse(nptcls), tmp(nonzero_nptcls), stat=alloc_stat )
                    if(alloc_stat /= 0) call alloc_errchk('In: commander_hlev_wflows::diverse_labeling ', alloc_stat)
                    rt = ran_tabu(nonzero_nptcls)
                    call rt%balanced(nlabels, tmp)
                    cnt = 0
                    do iptcl=1,nptcls
                        if(states(iptcl)>0 .and. states(iptcl)<=nlabels)then
                            cnt = cnt + 1
                            config_diverse(iptcl) = tmp(cnt)
                        else
                            config_diverse(iptcl) = states(iptcl)
                        endif
                    enddo
                    deallocate(tmp,states)
                else
                    allocate(config_diverse(nptcls), order(nptcls), tmp(nlabels), source=0, stat=alloc_stat )
                    corrs = os%get_all('corr')
                    tmp   = (/(s,s=1,nlabels)/)
                    where( states<=0 ) corrs = -1.
                    order = (/(iptcl,iptcl=1,nptcls)/)
                    call hpsort( corrs, order )
                    call reverse(order)
                    call reverse(corrs)
                    rt = ran_tabu(nlabels)
                    do iptcl = 1, nonzero_nptcls+nlabels, nlabels
                        call rt%reset
                        call rt%shuffle(tmp)
                        do s = 1, nlabels
                            ind = iptcl + s - 1
                            if(ind > nptcls)exit
                            config_diverse(order(ind)) = tmp(s)
                        enddo
                    enddo
                    where( states<=0 ) config_diverse = 0
                    deallocate(states,corrs,order,tmp)
                endif
                call rt%kill
            end subroutine diverse_labeling

    end subroutine exec_het

    !> multi-particle refinement after het
    subroutine exec_het_refine( self, cline )
        use simple_defs_conv
        use simple_commander_rec,    only: recvol_commander
        use simple_commander_volops, only: postproc_vol_commander
        class(het_refine_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        ! constants
        character(len=20),  parameter :: INIT_FBODY  = 'hetinit_refine_state'
        character(len=19),  parameter :: FINAL_FBODY = 'hetdoc_refine_state'
        character(len=17),  parameter :: FINAL_DOC   = 'hetdoc_refine'//METADATEXT
        ! distributed commanders
        type(prime3D_distr_commander) :: xprime3D_distr
        type(recvol_distr_commander)  :: xrecvol_distr
        ! shared-mem commanders
        type(postproc_vol_commander)  :: xpostproc_vol
        ! command lines
        type(cmdline)                 :: cline_prime3D_master
        type(cmdline)                 :: cline_prime3D
        type(cmdline)                 :: cline_recvol_distr
        type(cmdline)                 :: cline_postproc_vol
        ! other variables
        integer,               allocatable :: state_pops(:), states(:)
        character(len=STDLEN), allocatable :: init_docs(:), final_docs(:)
        logical,               allocatable :: l_hasmskvols(:), l_hasvols(:)
        type(params)          :: p_master
        type(oris)            :: os_master, os_state
        character(len=STDLEN) :: oritab, vol, fname
        character(len=9)      :: dir
        character(len=2)      :: str_state
        integer               :: state, iptcl, iter
        logical               :: l_singlestate, error

        ! make master parameters
        p_master      = params(cline)
        l_singlestate = cline%defined('state')
        error         = .false.

        ! sanity checks
        if( p_master%eo .eq. 'no' .and. .not. cline%defined('lp') )&
            &stop 'need lp input when eo .eq. no; het_refine'
        if(.not.file_exists(p_master%oritab))then
            print *,'Document ',trim(p_master%oritab),' does not exist!'
            stop
        endif

        ! general prep
        p_master%nptcls = nlines(p_master%oritab)
        call os_master%new(p_master%nptcls)
        call binread_oritab(p_master%oritab, os_master, [1,p_master%nptcls])
        p_master%nstates = os_master%get_n('state')
        if(p_master%nstates < 2 .and. .not.l_singlestate)then
            print *, 'Non-sensical number of states argument for heterogeneity refinemnt: ',p_master%nstates
            stop
        endif
        allocate(l_hasmskvols(p_master%nstates), l_hasvols(p_master%nstates), source=.false.)
        allocate(init_docs(p_master%nstates), final_docs(p_master%nstates))
        call os_master%get_pops(state_pops, 'state', consider_w=.false.)
        do state = 1, p_master%nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. state.ne.p_master%state )cycle
            str_state = int2str_pad(state,2)
            dir = 'state_'//str_state//'/'
            call mkdir(dir)
            ! preps individual documents
            os_state = os_master
            states   = nint(os_state%get_all('state'))
            where( states .ne. state )
                states = 0
            else where
                states = 1
            end where
            call os_state%set_all('state', real(states))
            deallocate(states)
            init_docs(state) = INIT_FBODY//str_state//METADATEXT
            call binwrite_oritab(init_docs(state), os_state, [1,p_master%nptcls])
            final_docs(state) = FINAL_FBODY//str_state//METADATEXT            
            ! check & move volumes
            l_hasvols(state) = trim(p_master%vols(state)) .ne. ''
            if( l_hasvols(state) )then
                if( p_master%eo .ne. 'no' )then
                    ! fsc
                    fname = FSC_FBODY//str_state//BIN_EXT
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                        error = .true.
                    else
                        call rename(fname, dir//trim(fname))
                    endif
                    ! FRC
                    fname = FRCS_FBODY//str_state//BIN_EXT
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                    else
                        call rename(fname, dir//trim(fname))
                    endif
                    ! aniso
                    fname = ANISOLP_FBODY//str_state//p_master%ext
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                    else
                        call rename(fname, dir//trim(fname))
                    endif
                    ! e/o
                    fname = add2fbody(trim(p_master%vols(state)), p_master%ext, '_even')
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                    else
                        call rename(fname, dir//trim(fname))
                    endif
                    fname = add2fbody(trim(p_master%vols(state)), p_master%ext, '_odd')
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                    else
                        call rename(fname, dir//trim(fname))
                    endif
                endif
                ! volume
                if( .not.file_exists(p_master%vols(state)) )then
                    print *, 'File missing: ', p_master%vols(state)
                    error = .true.
                else
                    fname = trim(p_master%vols(state))
                    p_master%vols(state) = dir//trim(fname) ! name change
                    call rename(fname, p_master%vols(state))
                endif
            endif
            ! mask volume
            l_hasmskvols(state) = trim(p_master%mskvols(state)) .ne. ''
            if( l_hasmskvols(state) )then
                if( .not.file_exists(p_master%mskvols(state)) )then
                    print *, 'File missing: ', trim(fname)
                    error = .true.
                else
                    fname = trim(p_master%mskvols(state))
                    p_master%mskvols(state) = dir//trim(fname)  ! name change
                    call exec_cmdline('cp '//trim(fname)//' '//trim(p_master%mskvols(state)))
                endif
            endif
        enddo
        if( error ) stop
        if( cline%defined('vollist') .and. count(l_hasvols).eq.0 )stop 'Missing volume(s)'
        if( cline%defined('msklist') .and. count(l_hasmskvols).eq.0 )stop 'Missing mask volume(s)'

        ! preps command-lines
        call cline%delete('vollist')
        call cline%delete('nstates')
        call cline%delete('msklist')
        call cline%delete('mskfile')
        do state = 1, p_master%nstates
            call cline%delete('vol'//int2str_pad(state,1))
        enddo
        cline_prime3D_master = cline
        cline_recvol_distr   = cline
        cline_postproc_vol   = cline
        call cline_prime3D_master%set('prg', 'prime3D')
        call cline_prime3D_master%set('dynlp', 'no')
        call cline_recvol_distr%set('prg', 'recvol')
        call cline_recvol_distr%set('oritab', trim(FINAL_DOC))
        call cline_recvol_distr%set('nstates', trim(int2str(p_master%nstates)))
        call cline_recvol_distr%set('eo', 'yes')
        call cline_postproc_vol%delete('lp')

        ! Main loop
        do state = 1, p_master%nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. state.ne.p_master%state )cycle
            str_state = int2str_pad(state,2)
            dir       = 'state_'//str_state//'/'
            write(*,'(A)')   '>>>'
            write(*,'(A,I2)')'>>> REFINING STATE: ', state
            write(*,'(A)')   '>>>'
            ! prime3D prep
            cline_prime3D = cline_prime3D_master
            call cline_prime3D%set('oritab', trim(init_docs(state)))
            if( l_hasvols(state) )then
                call cline_prime3D%set('vol1', trim(p_master%vols(state)))
                if( p_master%eo.ne.'no' )then
                    fname = dir//FSC_FBODY//str_state//BIN_EXT
                    call exec_cmdline('cp '//trim(fname)//' fsc_state01.bin')
                    fname = dir//FRCS_FBODY//str_state//BIN_EXT
                    if(file_exists(fname))call exec_cmdline('cp '//trim(fname)//' frcs_state01.bin')
                    fname = dir//ANISOLP_FBODY//str_state//p_master%ext               
                    if(file_exists(fname))call exec_cmdline('cp '//trim(fname)//' aniso_optlp_state01.mrc')
                endif
            endif
            if( l_hasmskvols(state) )call cline_prime3D%set('mskfile', trim(p_master%mskvols(state)))
            ! run prime3D
            call xprime3D_distr%execute(cline_prime3D)
            ! harvest outcome
            iter   = nint(cline_prime3D%get_rarg('endit'))
            oritab = PRIME3D_ITER_FBODY//int2str_pad(iter,3)//METADATEXT
            ! stash
            call os_state%new(p_master%nptcls)
            call binread_oritab(oritab, os_state, [1,p_master%nptcls])
            do iptcl = 1,p_master%nptcls
                if( nint(os_master%get(iptcl, 'state')) .ne. 1 )cycle
                call os_master%set_ori(iptcl, os_state%get_ori(iptcl))
                call os_master%set(iptcl, 'state', real(state))
            enddo
            call os_state%kill  
            call prime3d_cleanup
        enddo

        ! final document
        call binwrite_oritab(FINAL_DOC, os_master, [1,p_master%nptcls])

        ! final reconstruction
        if( p_master%eo .eq.'no' )then
            call xrecvol_distr%execute(cline_recvol_distr)
            do state = 1, p_master%nstates
                if( state_pops(state) == 0 )cycle
                if( l_singlestate .and. state.ne.p_master%state )cycle
                str_state = int2str_pad(state, 2)
                if( l_hasmskvols(state) )call cline_postproc_vol%set('mskfile', trim(p_master%mskvols(state)))
                vol = 'recvol_state'//trim(str_state)//p_master%ext
                call cline_postproc_vol%set('vol1', trim(vol))
                call cline_postproc_vol%set('fsc', FSC_FBODY//trim(str_state)//BIN_EXT)
                call cline_postproc_vol%set('vol_filt', ANISOLP_FBODY//trim(str_state)//p_master%ext)
                call xpostproc_vol%execute(cline_postproc_vol)
            enddo
        endif

        ! end gracefully
        call simple_end('**** SIMPLE_HET_REFINE NORMAL STOP ****')
        contains

            ! stash docs, volumes , etc.
            subroutine prime3d_cleanup
                character(len=STDLEN) :: src, dist
                character(len=2), parameter :: one = '01'
                character(len=3) :: str_iter
                integer          :: it
                dir = 'state_'//str_state//'/'
                call mkdir(dir)
                do it = 1, iter
                    str_iter = int2str_pad(it,3)
                    ! volumes
                    src  = VOL_FBODY//one//'_iter'//str_iter//p_master%ext
                    dist = dir//VOL_FBODY//one//'_iter'//str_iter//p_master%ext
                    call rename(src, dist)
                    src  = VOL_FBODY//one//'_iter'//str_iter//'_pproc'//p_master%ext
                    dist = dir//VOL_FBODY//one//'_iter'//str_iter//'_pproc'//p_master%ext
                    call rename(src, dist)
                    ! e/o
                    if( p_master%eo.ne.'no')then
                        src  = VOL_FBODY//one//'_iter'//str_iter//'_even'//p_master%ext
                        dist = dir//VOL_FBODY//one//'_iter'//str_iter//'_even'//p_master%ext
                        call rename(src, dist)
                        src  = VOL_FBODY//one//'_iter'//str_iter//'_odd'//p_master%ext
                        dist = dir//VOL_FBODY//one//'_iter'//str_iter//'_odd'//p_master%ext
                        call rename(src, dist)
                        src = 'RESOLUTION_STATE'//str_state//'_ITER'//str_iter
                        call rename(src, dir//src)
                    endif
                    ! orientation document
                    src = PRIME3D_ITER_FBODY//str_iter//METADATEXT
                    if( file_exists(src) )call rename(src, dir//src)
                enddo
                ! resolution measures
                if( p_master%eo.ne.'no')then
                    src  = FSC_FBODY//one//BIN_EXT
                    if( file_exists(src) )call rename(src, dir//src)
                    src  = FRCS_ITER_FBODY//one//BIN_EXT
                    if( file_exists(src) )call rename(src, dir//src)
                    src  = ANISOLP_FBODY//one//p_master%ext
                    if( file_exists(src) )call rename(src, dir//src)
                endif
            end subroutine prime3d_cleanup

    end subroutine exec_het_refine

end module simple_commander_hlev_wflows
