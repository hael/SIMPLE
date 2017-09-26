! concrete commander: high-level workflows
module simple_commander_hlev_wflows
use simple_cmdline,               only: cmdline
use simple_params,                only: params
use simple_commander_base,        only: commander_base
use simple_qsys_env,              only: qsys_env
use simple_commander_distr_wflows ! use all in there
use simple_fileio                 ! use all in there
use simple_commander_distr        ! use all in there
use simple_jiffys                 ! use all in there
use simple_binoris_io             ! use all in there
use simple_defs                   ! use all in there
implicit none

public :: prime2D_autoscale_commander
public :: ini3D_from_cavgs_commander
public :: het_ensemble_commander
public :: cga_hres_sel_commander
private

type, extends(commander_base) :: prime2D_autoscale_commander
  contains
    procedure :: execute      => exec_prime2D_autoscale
end type prime2D_autoscale_commander
type, extends(commander_base) :: ini3D_from_cavgs_commander
  contains
    procedure :: execute      => exec_ini3D_from_cavgs
end type ini3D_from_cavgs_commander
type, extends(commander_base) :: het_ensemble_commander
  contains
    procedure :: execute      => exec_het_ensemble
end type het_ensemble_commander
type, extends(commander_base) :: cga_hres_sel_commander
  contains
    procedure :: execute      => exec_cga_hres_sel
end type cga_hres_sel_commander

contains

    !> for distributed PRIME2D with two-stage autoscaling
    subroutine exec_prime2D_autoscale( self, cline )
        use simple_scaler,            only: scaler
        use simple_oris,              only: oris
        use simple_commander_prime2D, only: rank_cavgs_commander
        class(prime2D_autoscale_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        ! constants
        integer,           parameter :: MAXITS_STAGE1   = 10
        character(len=32), parameter :: CAVGS_ITERFBODY = 'cavgs_iter'
        character(len=32), parameter :: FINALDOC        = 'prime2Ddoc_final'//METADATEXT
        ! commanders
        type(split_commander)                        :: xsplit
        type(makecavgs_distr_commander)              :: xmakecavgs
        type(prime2D_distr_commander),       target  :: xprime2D_distr
        type(prime2D_chunk_distr_commander), target  :: xprime2D_chunk_distr
        class(commander_base),               pointer :: xprime2D => null()
        type(rank_cavgs_commander)                   :: xrank_cavgs
        ! command lines
        type(cmdline) :: cline_prime2D_stage1
        type(cmdline) :: cline_prime2D_stage2
        type(cmdline) :: cline_makecavgs
        type(cmdline) :: cline_rank_cavgs
        ! other variables
        type(oris)            :: os
        type(scaler)          :: scobj
        type(params)          :: p_master
        character(len=STDLEN) :: refs
        real                  :: scale_stage1, scale_stage2
        integer               :: nparts
        ! make master parameters
        p_master = params(cline, checkdistr=.false.)
        ! set pointer to the right commander
        if( cline%defined('chunksz') )then
            nparts = nint(real(p_master%nptcls)/real(p_master%chunksz))
            xprime2D => xprime2D_chunk_distr
            call cline%set('nparts', real(nparts))
        else
            nparts = p_master%nparts
            xprime2D => xprime2D_distr
        endif
        ! split stack
        call xsplit%execute(cline)
        if( p_master%l_autoscale )then
            ! auto-scaling prep (cline is modified by scobj%init)
            call scobj%init(p_master, cline, p_master%smpd_targets2D(1))
            scale_stage1 = scobj%get_scaled_var('scale')
            ! scale images in parallel
            call scobj%scale_distr_exec
            ! execute stage 1
            cline_prime2D_stage1 = cline
            call cline_prime2D_stage1%delete('automsk') ! deletes possible automsk flag from stage 1
            call cline_prime2D_stage1%set('maxits', real(MAXITS_STAGE1))
            call xprime2D%execute(cline_prime2D_stage1)
            ! prepare stage 2 input -- re-scale
            call scobj%uninit(cline) ! puts back the old command line
            call scobj%init(p_master, cline, p_master%smpd_targets2D(2))
            scale_stage2 = scobj%get_scaled_var('scale')
            call scobj%scale_distr_exec
            ! prepare stage 2 input -- shift modulation
            call os%new(p_master%nptcls)
            call binread_oritab(FINALDOC, os, [1,p_master%nptcls])
            call os%mul_shifts(scale_stage2/scale_stage1)
            call binwrite_oritab(FINALDOC, os, [1,p_master%nptcls])
            ! prepare stage 2 input -- command line
            cline_prime2D_stage2 = cline
            ! if automsk .eq. yes, we need to replace it with cavg
            if( p_master%automsk .eq. 'yes' )then
                call cline_prime2D_stage2%set('automsk', 'cavg')
            endif
            call cline_prime2D_stage2%delete('deftab')
            call cline_prime2D_stage2%set('oritab',  trim(FINALDOC))
            call cline_prime2D_stage2%set('startit', real(MAXITS_STAGE1 + 1))
            call xprime2D%execute(cline_prime2D_stage2)
            ! delete downscaled stack parts (we are done with them)
            call del_files(trim(STKPARTFBODY_SC), p_master%nparts, ext=p_master%ext)
            ! re-generate class averages at original sampling
            call scobj%uninit(cline) ! puts back the old command line
            call binread_oritab(FINALDOC, os, [1,p_master%nptcls])
            call os%mul_shifts(1./scale_stage2)
            call binwrite_oritab(FINALDOC, os, [1,p_master%nptcls])
            cline_makecavgs = cline
            call cline_makecavgs%delete('autoscale')
            call cline_makecavgs%delete('balance')
            call cline_makecavgs%delete('chunksz')
            if( p_master%l_chunk_distr )then
                call cline_makecavgs%delete('ncls')
            endif
            call cline_makecavgs%set('prg',     'makecavgs')
            call cline_makecavgs%set('oritab',  trim(FINALDOC))
            call cline_makecavgs%set('nparts',  real(nparts))
            call cline_makecavgs%set('refs',    'cavgs_final'//p_master%ext)
            call xmakecavgs%execute(cline_makecavgs)
        else
            call xprime2D%execute(cline)
        endif
        ! ranking
        call cline_rank_cavgs%set('oritab',   trim(FINALDOC))
        call cline_rank_cavgs%set('stk',      'cavgs_final'//p_master%ext)
        call cline_rank_cavgs%set('classdoc', 'classdoc.txt')
        call cline_rank_cavgs%set('outstk',   'cavgs_final_ranked'//p_master%ext)
        call xrank_cavgs%execute( cline_rank_cavgs )
        ! cleanup
        call del_file('prime2D_startdoc'//METADATEXT)
        call del_file('start2Drefs'//p_master%ext)
        call del_files(STKPARTFBODY_SC, p_master%nparts, ext=p_master%ext)
        ! end gracefully
        call simple_end('**** SIMPLE_PRIME2D NORMAL STOP ****')
    end subroutine exec_prime2D_autoscale

    !> for generation of an initial 3d model from class averages
    subroutine exec_ini3D_from_cavgs( self, cline )
        use simple_commander_volops, only: projvol_commander
        use simple_commander_rec,    only: recvol_commander
        use simple_strings,          only: int2str_pad, str2int
        use simple_scaler,           only: scaler
        use simple_oris,             only: oris
        class(ini3D_from_cavgs_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! constants
        real,                  parameter :: LPLIMS(2)=[20.,10.] !< default low-pass limits
        real,                  parameter :: CENLP=30.           !< consistency with prime3D
        integer,               parameter :: MAXITS_SNHC=30, MAXITS_INIT=15, MAXITS_REFINE=40
        integer,               parameter :: STATE=1, NPROJS_SYMSRCH=50, NPEAKS_REFINE=6
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
        real                  :: iter, lpstop, smpd_target, smpd
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
        call del_files(STKPARTFBODY_SC, p_master%nparts, ext=p_master%ext)
        ! make master parameters
        p_master = params(cline, checkdistr=.false.)
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
        smpd_target = p_master%smpd
        if( doautoscale )then
            if( cline%defined('lp') )then
                smpd_target = p_master%lp*LP2SMPDFAC
            else if( cline%defined('lpstop') )then
                smpd_target = min(LPLIMS(2),p_master%lpstop)*LP2SMPDFAC
            else
                smpd_target = LPLIMS(2)*LP2SMPDFAC
            endif
            call scobj%init(p_master, cline, smpd_target, STKSCALEDBODY)
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
        call cline_prime3D_snhc%set('prg',    'prime3D')
        call cline_prime3D_snhc%set('ctf',    'no')
        call cline_prime3D_snhc%set('maxits', real(MAXITS_SNHC))
        call cline_prime3D_snhc%set('dynlp',  'yes') ! better be explicit about the dynlp
        call cline_prime3D_snhc%set('refine', 'snhc')
        ! (2) PRIME3D_INIT
        call cline_prime3D_init%set('prg',    'prime3D')
        call cline_prime3D_init%set('ctf',    'no')
        call cline_prime3D_init%set('maxits', real(MAXITS_INIT))
        call cline_prime3D_init%set('dynlp',  'no') ! better be explicit about the dynlp
        call cline_prime3D_init%set('vol1',   trim(SNHCVOL)//trim(str_state)//p_master%ext)
        call cline_prime3D_init%set('oritab', SNHCDOC)
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
        call cline_prime3D_refine%set('dynlp', 'no') ! better be explicit about the dynlp
        call cline_prime3D_refine%set('lp', LPLIMS(2))
        call cline_prime3D_refine%set('refine', 'no')
        call cline_prime3D_refine%set('npeaks', real(NPEAKS_REFINE))
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
        call update_lp(cline_prime3D_init, 1)
        call xprime3D_distr%execute(cline_prime3D_init)
        iter = cline_prime3D_init%get_rarg('endit')
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
            ! refinement step needs to use iter dependent vol/oritab
            call cline_prime3D_refine%set('oritab', trim(oritab))
            call cline_prime3D_refine%set('vol1', trim(vol_iter))
        endif
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> PRIME3D REFINEMENT STEP'
        write(*,'(A)') '>>>'
        call cline_prime3D_refine%set('startit', iter + 1.0)
        call update_lp(cline_prime3D_refine, 2)
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
            call rename(trim(volfbody)//trim(str_state)//p_master%ext, 'rec_final'//p_master%ext)
        else
            call rename(trim(vol_iter), 'rec_final'//p_master%ext)
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

    !> for ensemble heterogeinity analysis
    subroutine exec_het_ensemble( self, cline )
        use simple_commander_rec,    only: recvol_commander
        use simple_commander_volops, only: postproc_vol_commander
        use simple_strings,          only: int2str_pad
        use simple_oris,             only: oris
        use simple_combinatorics,    only: diverse_labeling, shc_aggregation
        use simple_fileio,           only: file_exists, del_file
        class(het_ensemble_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        ! constants
        integer,            parameter :: MAXITS_INIT = 50
        character(len=32),  parameter :: HETFBODY    = 'hetrep_'
        character(len=32),  parameter :: REPEATFBODY = 'hetdoc_'
        character(len=32),  parameter :: VOLFBODY    = 'recvol_state'
        ! distributed commanders
        type(prime3D_distr_commander) :: xprime3D_distr
        type(recvol_distr_commander)  :: xrecvol_distr
        ! shared-mem commanders
        type(postproc_vol_commander)  :: xpostproc_vol
        ! command lines
        type(cmdline)                 :: cline_prime3D_master
        type(cmdline)                 :: cline_prime3D
        type(cmdline)                 :: cline_recvol_distr
        type(cmdline)                 :: cline_postproc_vol, cline_postproc_repvol
        ! other variables
        real,                  allocatable :: rep_corrs(:)
        integer,               allocatable :: labels(:,:), labels_incl(:,:), consensus(:)
        character(len=STDLEN), allocatable :: init_docs(:), final_docs(:)
        logical,               allocatable :: included(:)
        type(params)          :: p_master
        type(oris)            :: os
        character(len=STDLEN) :: oritab, vol, str_state
        integer               :: irepeat, state, iter, n_incl, it, alloc_stat
        integer               :: best_loc(1)

        ! sanity check
        if(nint(cline%get_rarg('nstates')) <= 1)&
            &stop 'Non-sensical NSTATES argument for heterogeneity analysis!'

        ! make master parameters
        p_master = params(cline, checkdistr=.false.)

        if( p_master%eo .eq. 'no' .and. .not. cline%defined('lp') )&
            &stop 'need lp input when eo .ne. no; het_ensemble'

        allocate(init_docs(p_master%nrepeats), final_docs(p_master%nrepeats), rep_corrs(p_master%nrepeats), stat=alloc_stat)
        do irepeat=1,p_master%nrepeats
            oritab              = trim(REPEATFBODY)//'init_rep'//int2str_pad(irepeat,2)//METADATEXT
            init_docs(irepeat)  = trim(oritab)
            oritab              = trim(REPEATFBODY)//'rep'//int2str_pad(irepeat,2)//METADATEXT
            final_docs(irepeat) = trim(oritab)
            call del_file(init_docs(irepeat))
            call del_file(final_docs(irepeat))
        enddo

        ! delete possibly pre-existing scaled stack parts & pft parts
        call del_files(STKPARTFBODY_SC, p_master%nparts, ext=p_master%ext)
        call del_files('ppfts_memoized_part', p_master%nparts, ext='.bin')

        ! prepare command lines from prototype
        cline_recvol_distr = cline
        call cline_recvol_distr%set('prg', 'recvol')
        if( p_master%eo .eq. 'no' ) call cline_recvol_distr%set('eo', 'yes')
        call cline_recvol_distr%delete('lp')
        cline_postproc_vol = cline
        call cline_postproc_vol%set('prg', 'postproc_vol')
        if( p_master%eo .eq. 'no' )  call cline_postproc_vol%set('eo', 'yes')
        call cline_postproc_vol%delete('lp')
        cline_postproc_repvol = cline
        call cline_postproc_repvol%set('prg', 'postproc_vol')
        cline_prime3D_master = cline
        call cline_prime3D_master%set('prg', 'prime3D')
        call cline_prime3D_master%set('startit', 1.)
        call cline_prime3D_master%set('maxits', real(MAXITS_INIT))
        call cline_prime3D_master%set('refine', 'het')
        call cline_prime3D_master%set('dynlp', 'no')
        call cline_prime3D_master%set('pproc', 'no')

        ! generate diverse initial labels
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> GENERATING DIVERSE LABELING'
        call os%new(p_master%nptcls)
        call binread_oritab(p_master%oritab, os, [1,p_master%nptcls])
        labels   = diverse_labeling(p_master%nptcls, p_master%nstates, p_master%nrepeats)
        included = os%included()
        n_incl   = count(included)
        do irepeat=1,p_master%nrepeats
            where( .not.included ) labels(irepeat,:) = 0
            call os%set_all('state', real(labels(irepeat,:)))
            call binwrite_oritab(trim(init_docs(irepeat)), os, [1,p_master%nptcls])
        enddo

        ! generate candidate solutions
        do irepeat = 1,p_master%nrepeats
            write(*,'(A)')    '>>>'
            write(*,'(A,I3)') '>>> PRIME3D REPEAT ', irepeat
            write(*,'(A)')    '>>>'
            ! run prime3d
            cline_prime3D = cline_prime3D_master
            call cline_prime3D%set('oritab', init_docs(irepeat))
            call xprime3D_distr%execute(cline_prime3D)
            ! harvest outcome
            iter   = nint(cline_prime3D%get_rarg('endit'))
            oritab = 'prime3Ddoc_'//int2str_pad(iter,3)//METADATEXT
            call rename(trim(oritab), trim(final_docs(irepeat)))
            call binread_oritab(trim(final_docs(irepeat)), os, [1,p_master%nptcls])
            ! updates labels & correlations
            labels(irepeat,:)  = nint(os%get_all('state'))
            rep_corrs(irepeat) = sum(os%get_all('corr'), mask=included) / real(n_incl)
            ! process final volumes
            call pproc_volumes
            ! cleanup
            call prime3d_cleanup
        enddo
        best_loc = maxloc(rep_corrs)

        ! generate consensus document
        oritab = trim(REPEATFBODY)//'consensus'//METADATEXT
        call del_file(oritab)
        write(*,'(A)')   '>>>'
        write(*,'(A,A)') '>>> GENERATING ENSEMBLE SOLUTION: ', trim(oritab)
        write(*,'(A)')   '>>>'
        if( p_master%nrepeats == 1 )then
            ! done
        else 
            ! aggregate solutions
            call binread_oritab(p_master%oritab, os, [1,p_master%nptcls])
            allocate(labels_incl(p_master%nrepeats,n_incl), consensus(n_incl), stat=alloc_stat)
            do irepeat=1,p_master%nrepeats
                labels_incl(irepeat,:) = pack(labels(irepeat,:), mask=included)
            enddo
            labels(1,:) = 0
            call shc_aggregation(p_master%nrepeats, n_incl, labels_incl, consensus, best_loc(1))
            call os%set_all('state', real(unpack(consensus, included, labels(1,:))) )
            deallocate(labels_incl, consensus, stat=alloc_stat)
        endif
        ! output
        call binwrite_oritab(trim(oritab), os, [1,p_master%nptcls])
        ! cleanup
        call os%kill
        deallocate(init_docs, final_docs, labels, included, rep_corrs)

        ! final distributed reconstruction
        call cline_recvol_distr%set('oritab', trim(oritab))
        call xrecvol_distr%execute(cline_recvol_distr)
        ! post-process
        do state = 1, p_master%nstates
            str_state = int2str_pad(state, 2)
            vol = 'recvol_state'//trim(str_state)//p_master%ext
            call cline_postproc_vol%set('vol1', trim(vol))
            call update_pproc_cline(cline_postproc_vol, str_state)
            call xpostproc_vol%execute(cline_postproc_vol)
        enddo

        ! end gracefully
        call simple_end('**** SIMPLE_HET_ENSEMBLE NORMAL STOP ****')
        contains

            subroutine prime3d_cleanup
                character(len=STDLEN) :: fname, vol
                ! delete starting volumes
                do state = 1, p_master%nstates
                    vol = 'startvol_state'//int2str_pad(state,2)//p_master%ext
                    if(file_exists(vol))call del_file(vol)
                enddo
                ! delete iterations volumes & documents
                do it = 1, iter-1
                    do state = 1, p_master%nstates
                        vol = trim(VOLFBODY)//int2str_pad(state,2)//'_iter'//int2str_pad(it,3)//p_master%ext
                        if(file_exists(vol))call del_file(vol)
                    enddo
                    oritab = 'prime3Ddoc_'//int2str_pad(it,3)//METADATEXT
                    if(file_exists(oritab))call del_file(oritab)
                enddo
            end subroutine prime3d_cleanup

            ! renames & post-process final volumes of a repeat for all states
            subroutine pproc_volumes
                character(len=STDLEN) :: srcvol, destvol
                do state=1,p_master%nstates
                    str_state = int2str_pad(state, 2)
                    srcvol  = trim(VOLFBODY)//int2str_pad(state,2)//'_iter'//int2str_pad(iter,3)//p_master%ext
                    destvol = trim(HETFBODY)//int2str_pad(irepeat,2)//'_'//trim(VOLFBODY)//int2str_pad(state,2)//p_master%ext
                    if(file_exists(destvol))call del_file(destvol)
                    call rename(trim(srcvol), trim(destvol))
                    call cline_postproc_repvol%set('vol1', trim(destvol))
                    call update_pproc_cline(cline_postproc_repvol, str_state)
                    call xpostproc_vol%execute(cline_postproc_repvol)
                enddo
            end subroutine pproc_volumes

            subroutine update_pproc_cline( cl, str_state )
                class(cmdline),   intent(inout) :: cl
                character(len=2), intent(in)    :: str_state
                character(len=STDLEN) :: fsc_file, optlp_file
                fsc_file   = 'fsc_state'//trim(str_state)//'.bin'
                optlp_file = 'aniso_optlp_state'//trim(str_state)//p_master%ext
                if( file_exists(fsc_file) .and. p_master%eo .eq. 'aniso' )then
                    if( .not. file_exists(optlp_file) )then
                        write(*,*) 'eo=aniso but file: ', trim(optlp_file)
                        stop 'is not in cwd as required; commander_hlev_wflows :: exec_het_ensemble'
                    endif
                    call cl%delete('lp')
                    call cl%set('fsc', trim(fsc_file))
                    call cl%set('vol_filt', trim(optlp_file))
                else if( file_exists(fsc_file) .and. p_master%eo .eq. 'yes' )then
                    call cl%delete('lp')
                    call cl%set('fsc', trim(fsc_file))
                else
                    call cl%delete('fsc')
                    call cl%set('lp', p_master%lp)
                endif
            end subroutine update_pproc_cline

    end subroutine exec_het_ensemble

    subroutine exec_cga_hres_sel( self, cline )
        use simple_commander_rec, only: recvol_commander
        use simple_rnd,           only: branarr
        use simple_oris,          only: oris
        use simple_math,          only: fsc1_ge_fsc2
        class(cga_hres_sel_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        ! constants
        integer,           parameter :: MAXITS        = 50
        character(len=32), parameter :: VOLFBODY      = 'recvol_state'
        character(len=32), parameter :: ORIS_MODIFIED = 'oris_from_cga_hres_sel'//METADATEXT
        character(len=32), parameter :: FSC_FILE      = 'fsc_state01.bin'
        ! distributed commanders
        type(recvol_distr_commander) :: xrecvol_distr
        ! command lines
        type(cmdline)                :: cline_recvol_distr
        type(cmdline)                :: cline_postproc_vol
        ! other variables
        integer                      :: iter, iptcl
        real                         :: eps, frac_srch_space, param_overlap, rstate
        real,            allocatable :: probs(:), probs_conv(:), fsc1(:), fsc2(:)
        integer, target, allocatable :: S1(:), S2(:)
        integer,         allocatable :: Sbest(:), Sprev_best(:)
        integer, pointer             :: Swinner(:) => null(), Sloser(:) => null()
        type(params)                 :: p_master
        type(oris)                   :: os, os_mod
        ! default for now
        call cline%set('eo', 'yes')
        ! default always
        if( .not. cline%defined('maxits') ) call cline%set('maxits', real(MAXITS))
        if( .not. cline%defined('oritab') ) stop 'need oritab input; commander_hlev :: exec_cga_hres_sel'
        ! make master parameters
        p_master = params(cline, checkdistr=.false.)
        ! set learning rate
        if( .not. cline%defined('eps') ) p_master%eps = 5.0 / real(2.0 * p_master%maxits)
        ! prep recvol cline
        cline_recvol_distr = cline
        call cline_recvol_distr%set('prg',    'recvol')
        call cline_recvol_distr%set('oritab', ORIS_MODIFIED)
        ! prepare oris
        call os%new(p_master%nptcls)
        call binread_oritab(p_master%oritab, os, [1,p_master%nptcls])
        os_mod = os
        ! allocate and initialise arrays
        allocate(probs(p_master%nptcls), probs_conv(p_master%nptcls), source=0.5)
        allocate(Sbest(p_master%nptcls), Sprev_best(p_master%nptcls), source=1)
        do iter=1,p_master%maxits
            write(*,'(A)')   '>>>'
            write(*,'(A,I6)')'>>> ITERATION ', iter
            write(*,'(A)')   '>>>' 
            ! make two Bernoulli samples from probs
            S1   = branarr(probs)
            S2   = branarr(probs)
            ! distributed generation of FSC functions
            fsc1 = gen_fsc(S1, iter)
            fsc2 = gen_fsc(S2, iter)
            ! stash previous best
            Sprev_best = Sbest
            ! compare FSCs and set winner/loser pointers
            if( fsc1_ge_fsc2(fsc1,fsc2) )then
                Sbest   =  S1
                Swinner => S1
                Sloser  => S2
            else
                Sbest   =  S2
                Swinner => S2
                Sloser  => S1
            endif
            ! update probabilistic model
            where( Swinner /= Sloser .and. Swinner == 1 )
                probs = probs + p_master%eps
            elsewhere( Swinner /= Sloser .and. Swinner == 0 )
                probs = probs - p_master%eps
            end where
            ! check convergence
            where( probs >= 0.5 )
                probs_conv = probs
            else where
                probs_conv = 1.0 - probs
            end where
            frac_srch_space = (sum(probs_conv)/real(p_master%nptcls)) * 100.0
            param_overlap   = real(count(Sbest == Sprev_best)) / real(p_master%nptcls)
            write(*,'(A,1X,F7.4)') '>>> PARAMETER DISTRIBUTION OVERLAP:    ', param_overlap
            write(*,'(A,1X,F7.1)') '>>> PERCENTAGE OF SEARCH SPACE SCANNED:', frac_srch_space
        end do
        ! generate final rec + FSC
        fsc1 = gen_fsc(Sbest)
        ! report the final solution as a state 1/0 labeling
        do iptcl=1,p_master%nptcls
            rstate = os%get(iptcl, 'state')
            if( rstate > 0.5 .and. Sbest(iptcl) == 0 )then
                call os%set(iptcl, 'state', 0.)
            endif
        end do
        call binwrite_oritab('cga_hres_seldoc_final'//METADATEXT, os, [1,p_master%nptcls])
        ! end gracefully
        call simple_end('**** SIMPLE_CGA_HRES_SEL NORMAL STOP ****')

        contains

            function gen_fsc( S, iter ) result( fsc )
                integer,           intent(in) :: S(p_master%nptcls)
                integer, optional, intent(in) :: iter
                real, allocatable   :: fsc(:)
                call del_file(FSC_FILE)
                ! set reconstruction weights according to Bernoulli sample
                ! in this way any previous state 1/0 labeling will still be considered
                call os_mod%set_all('w', real(S))
                call binwrite_oritab(ORIS_MODIFIED, os_mod, [1,p_master%nptcls])
                ! reconstruct
                call xrecvol_distr%execute(cline_recvol_distr)
                ! get FSC values
                if( file_exists(FSC_FILE) )then
                    fsc = file2rarr(FSC_FILE)
                else
                    write(*,*) 'ERROR, file: ', trim(FSC_FILE)
                    stop 'does not exist in cwd; commander_hlev_wflows :: exec_cga_hres_sel'
                endif
                if( present(iter) )then
                    call del_file('RESOLUTION_STATE01')
                endif
            end function gen_fsc


    end subroutine exec_cga_hres_sel

end module simple_commander_hlev_wflows
