! concrete commander: high-level workflows
module simple_commander_hlev_wflows
include 'simple_lib.f08'
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_sp_project,     only: sp_project
use simple_parameters,     only: parameters
implicit none

public :: cluster2D_autoscale_commander
public :: initial_3Dmodel_commander
public :: cluster3D_commander
public :: cluster3D_refine_commander
private

type, extends(commander_base) :: cluster2D_autoscale_commander
  contains
    procedure :: execute      => exec_cluster2D_autoscale
end type cluster2D_autoscale_commander
type, extends(commander_base) :: initial_3Dmodel_commander
  contains
    procedure :: execute      => exec_initial_3Dmodel
end type initial_3Dmodel_commander
type, extends(commander_base) :: cluster3D_commander
  contains
    procedure :: execute      => exec_cluster3D
end type cluster3D_commander
type, extends(commander_base) :: cluster3D_refine_commander
  contains
    procedure :: execute      => exec_cluster3D_refine
end type cluster3D_refine_commander

contains

    !> for distributed CLUSTER2D with two-stage autoscaling
    subroutine exec_cluster2D_autoscale( self, cline )
        use simple_commander_distr_wflows, only: make_cavgs_distr_commander,cluster2D_distr_commander,scale_project_distr_commander
        use simple_commander_cluster2D,    only: rank_cavgs_commander
        use simple_commander_imgproc,      only: scale_commander
        class(cluster2D_autoscale_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        ! constants
        integer, parameter :: MAXITS_STAGE1      = 10
        integer, parameter :: MAXITS_STAGE1_EXTR = 15
        ! commanders
        type(make_cavgs_distr_commander)    :: xmake_cavgs
        type(cluster2D_distr_commander)     :: xcluster2D_distr
        type(rank_cavgs_commander)          :: xrank_cavgs
        type(scale_commander)               :: xscale
        type(scale_project_distr_commander) :: xscale_distr
        ! command lines
        type(cmdline) :: cline_cluster2D_stage1
        type(cmdline) :: cline_cluster2D_stage2
        type(cmdline) :: cline_scalerefs, cline_scale1, cline_scale2
        type(cmdline) :: cline_make_cavgs
        type(cmdline) :: cline_rank_cavgs
        ! other variables
        type(parameters)              :: params
        type(sp_project)              :: spproj, spproj_sc
        character(len=:), allocatable :: projfile_sc, stk
        character(len=LONGSTRLEN)     :: finalcavgs, finalcavgs_ranked, refs_sc
        real     :: scale_stage1, scale_stage2
        integer  :: istk, nparts, last_iter_stage1, last_iter_stage2
        logical  :: scaling
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        nparts = params%nparts
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! read project file
        call spproj%read(params%projfile)
        ! sanity checks
        if( spproj%get_nptcls() == 0 )then
            write(*,*)'No particles found in project file: ', trim(params%projfile), 'simple_commander_hlev_wflows::exec_cluster2D_autoscale'
            stop 'No particles found in project file: simple_commander_hlev_wflows::exec_cluster2D_autoscale'
        endif
        ! delete any previous solution
        if( .not. spproj%is_virgin_field(params%oritype) )then
            ! removes previous cluster2D solution (states are preserved)
            call spproj%os_ptcl2D%delete_2Dclustering
            call spproj%write_segment_inside(params%oritype)
        endif
        ! refinement flag
        if(.not.cline%defined('refine')) call cline%set('refine','snhc')
        ! splitting
        call spproj%split_stk(params%nparts, (params%mkdir.eq.'yes'), dir=PATH_PARENT)
        ! general options planning
        if( params%l_autoscale )then
            ! this workflow executes two stages of CLUSTER2D
            ! Stage 1: high down-scaling for fast execution, hybrid extremal/SHC optimisation for
            !          improved population distribution of clusters, no incremental learning,
            !          objective function is standard cross-correlation (cc)
            cline_cluster2D_stage1 = cline
            if( cline%get_carg('objfun').eq.'euclid' )then
                call cline_cluster2D_stage1%set('objfun', 'euclid')
            else
                call cline_cluster2D_stage1%set('objfun', 'cc')
            endif
            if( params%l_frac_update )then
                call cline_cluster2D_stage1%delete('update_frac') ! no incremental learning in stage 1
                call cline_cluster2D_stage1%set('maxits', real(MAXITS_STAGE1_EXTR))
            else
                call cline_cluster2D_stage1%set('maxits', real(MAXITS_STAGE1))
            endif
            ! Scaling
            call spproj%scale_projfile(params%smpd_targets2D(1), projfile_sc,&
                &cline_cluster2D_stage1, cline_scale1, dir=trim(STKPARTSDIR))
            call spproj%kill
            scale_stage1 = cline_scale1%get_rarg('scale')
            scaling      = trim(projfile_sc) /= trim(params%projfile)
            if( scaling )then
                call simple_mkdir(trim(STKPARTSDIR),errmsg="commander_hlev_wflows :: exec_cluster2D_autoscale;  ")
                call xscale_distr%execute( cline_scale1 )
                ! scale references
                if( cline%defined('refs') )then
                    call cline_scalerefs%set('stk', trim(params%refs))
                    refs_sc = 'refs'//trim(SCALE_SUFFIX)//params%ext
                    call cline_scalerefs%set('outstk', trim(refs_sc))
                    call cline_scalerefs%set('smpd', params%smpd)
                    call cline_scalerefs%set('newbox', cline_scale1%get_rarg('newbox'))
                    call xscale%execute(cline_scalerefs)
                    call cline_cluster2D_stage1%set('refs',trim(refs_sc))
                endif
            endif
            ! execution
            call cline_cluster2D_stage1%set('projfile', trim(projfile_sc))
            call xcluster2D_distr%execute(cline_cluster2D_stage1)
            last_iter_stage1 = nint(cline_cluster2D_stage1%get_rarg('endit'))
            ! update original project
            if( scaling )then
                call spproj_sc%read_segment( 'ptcl2D', projfile_sc )
                call spproj_sc%os_ptcl2D%mul_shifts( 1./scale_stage1 )
                call spproj%read( params%projfile )
                spproj%os_ptcl2D = spproj_sc%os_ptcl2D
                call spproj%write_segment_inside('ptcl2D')
                call spproj%kill()
                ! clean stacks
                call spproj_sc%read_segment( 'stk', projfile_sc )
                do istk=1,spproj_sc%os_stk%get_noris()
                    call spproj_sc%os_stk%getter(istk, 'stk', stk)
                    call del_file(trim(stk))
                enddo
                call spproj_sc%kill()
                call del_file(trim(projfile_sc))
            endif
            deallocate(projfile_sc)
            ! Stage 2: refinement stage, less down-scaling, no extremal updates, incremental
            !          learning for acceleration, objective function is resolution weighted
            !          cross-correlation with automtic fitting of B-factors
            cline_cluster2D_stage2 = cline
            if( cline%get_carg('objfun').eq.'euclid' )then
                call cline_cluster2D_stage2%set('objfun', 'euclid')
            else
                call cline_cluster2D_stage2%set('objfun', 'ccres')
            endif
            call cline_cluster2D_stage2%delete('refs')
            call cline_cluster2D_stage2%set('startit', real(last_iter_stage1 + 1))
            if( cline%defined('update_frac') )then
                call cline_cluster2D_stage2%set('update_frac', params%update_frac)
            endif
            ! Scaling
            call spproj%read( params%projfile )
            call spproj%scale_projfile( params%smpd_targets2D(2), projfile_sc,&
                &cline_cluster2D_stage2, cline_scale2, dir=trim(STKPARTSDIR))
            call spproj%kill
            scale_stage2 = cline_scale2%get_rarg('scale')
            scaling      = trim(projfile_sc) /= params%projfile
            if( scaling ) call xscale_distr%execute( cline_scale2 )
            ! execution
            call cline_cluster2D_stage2%set('projfile', trim(projfile_sc))
            call xcluster2D_distr%execute(cline_cluster2D_stage2)
            last_iter_stage2 = nint(cline_cluster2D_stage2%get_rarg('endit'))
            finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//params%ext
            ! Updates project and references
            if( scaling )then
                ! shift modulation
                call spproj_sc%read_segment( 'ptcl2D', projfile_sc )
                call spproj_sc%os_ptcl2D%mul_shifts( 1./scale_stage2 )
                call spproj%read( params%projfile )
                spproj%os_ptcl2D = spproj_sc%os_ptcl2D
                call spproj%write_segment_inside('ptcl2D')
                call spproj%kill()
                ! clean stacks
                call spproj_sc%read_segment( 'stk', projfile_sc )
                do istk=1,spproj_sc%os_stk%get_noris()
                    call spproj_sc%os_stk%getter(istk, 'stk', stk)
                    call del_file(trim(stk))
                enddo
                call spproj_sc%kill()
                call del_file(trim(projfile_sc))
                ! original scale references
                cline_make_cavgs = cline ! ncls is transferred here
                call cline_make_cavgs%delete('autoscale')
                call cline_make_cavgs%delete('balance')
                call cline_make_cavgs%set('prg',      'make_cavgs')
                call cline_make_cavgs%set('projfile', params%projfile)
                call cline_make_cavgs%set('nparts',   real(nparts))
                call cline_make_cavgs%set('refs',     trim(finalcavgs))
                call xmake_cavgs%execute(cline_make_cavgs)
            endif
        else
            ! no auto-scaling
            call xcluster2D_distr%execute(cline)
            last_iter_stage2 = nint(cline%get_rarg('endit'))
            finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//params%ext
        endif
        ! adding cavgs & FRCs to project & perform match filtering
        call spproj%read( params%projfile )
        call spproj%add_frcs2os_out( trim(FRCS_FILE), 'frc2D')
        call spproj%add_cavgs2os_out( trim(finalcavgs), spproj%get_smpd())
        call spproj%write_segment_inside('out')
        call spproj%kill()
        ! ranking
        finalcavgs_ranked = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//'_ranked'//params%ext
        call cline_rank_cavgs%set('projfile', trim(params%projfile))
        call cline_rank_cavgs%set('stk',      trim(finalcavgs))
        call cline_rank_cavgs%set('outstk',   trim(finalcavgs_ranked))
        call xrank_cavgs%execute( cline_rank_cavgs )
        ! cleanup
        call del_file('start2Drefs'//params%ext)
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
    end subroutine exec_cluster2D_autoscale

    !> for generation of an initial 3d model from class averages
    subroutine exec_initial_3Dmodel( self, cline )
        use simple_commander_distr_wflows, only: refine3D_distr_commander, scale_project_distr_commander
        use simple_oris,                   only: oris
        use simple_commander_volops,       only: reproject_commander, symaxis_search_commander
        use simple_parameters,             only: params_glob
        use simple_qsys_env,               only: qsys_env
        class(initial_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! constants
        real,    parameter :: CENLP=30. !< consistency with refine3D
        integer, parameter :: STATE=1, MAXITS_SNHC=30, MAXITS_INIT=15, MAXITS_REFINE=40
        integer, parameter :: NSPACE_SNHC=1000, NSPACE_INIT=1000, NSPACE_REFINE= 2500
        character(len=STDLEN), parameter :: ORIG_WORK_PROJFILE = 'initial_3Dmodel_tmpproj.simple'
        ! distributed commanders
        type(refine3D_distr_commander)      :: xrefine3D_distr
        type(scale_project_distr_commander) :: xscale_distr
        ! shared-mem commanders
        type(symaxis_search_commander) :: xsymsrch
        type(reproject_commander)      :: xreproject
        ! command lines
        type(cmdline) :: cline_refine3D_snhc
        type(cmdline) :: cline_refine3D_init
        type(cmdline) :: cline_refine3D_refine
        type(cmdline) :: cline_symsrch
        type(cmdline) :: cline_reconstruct3D
        type(cmdline) :: cline_reproject
        type(cmdline) :: cline_scale
        ! other variables
        character(len=:), allocatable :: stk, orig_stk, frcs_fname
        character(len=:), allocatable :: WORK_PROJFILE
        real,             allocatable :: res(:)
        integer,          allocatable :: states(:)
        type(qsys_env)        :: qenv
        type(parameters)      :: params
        type(ctfparams)       :: ctfvars ! ctf=no by default
        type(sp_project)      :: spproj, work_proj1, work_proj2
        type(oris)            :: os
        character(len=2)      :: str_state
        character(len=STDLEN) :: vol_iter
        real                  :: iter, smpd_target, lplims(2), msk, scale_factor, orig_msk
        integer               :: icls, ncavgs, orig_box, box, istk, status
        logical               :: srch4symaxis, do_autoscale, do_eo
        ! hard set oritype
        call cline%set('oritype', 'out') ! because cavgs are part of out segment
        ! auto-scaling prep
        do_autoscale = (cline%get_carg('autoscale').eq.'yes')
        ! now, remove autoscale flag from command line, since no scaled partial stacks
        ! will be produced (this program used shared-mem paralllelisation of scale)
        call cline%delete('autoscale')
        ! make master parameters
        do_eo = trim(cline%get_carg('eo')) .eq. 'yes'
        call cline%set('eo','no')
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! from now on we are in the ptcl3D segment, final report is in the cls3D segment
        call cline%set('oritype', 'ptcl3D')
        ! set global state string
        str_state = int2str_pad(STATE,2)
        ! decide wether to search for the symmetry axis
        srch4symaxis = params%pgrp .ne. 'c1'
        ! read project & update sampling distance
        call spproj%read(params%projfile)
        ! retrieve cavgs stack & FRCS info
        call spproj%get_cavgs_stk(stk, ncavgs, ctfvars%smpd)
        if( .not.spproj%os_cls2D%isthere('state') )then
            ! start from import
            allocate(states(ncavgs), source=1)
        else
            ! start from previous 2D
            states = nint(spproj%os_cls2D%get_all('state'))
        endif
        if( count(states==0) .eq. ncavgs )then
            write(*,*) 'No class averages detected in project file: ',trim(params%projfile), &
                '; simple_commander_hlev_wflows::initial_3Dmodel'
            stop 'No class averages detected in project file ; simple_commander_hlev_wflows::initial_3Dmodel'
        endif
        ! sanity check
        if( do_eo )then
            call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
            if( .not.file_exists(frcs_fname) )then
                write(*,*)'The project file does not contain the required information for e/o alignment, use eo=no instead'
                stop 'The project file does not contain the required information for e/o alignment'
            endif
        endif
        ! init
        params%smpd = ctfvars%smpd
        orig_stk    = stk
        ! set lplims
        lplims(1) = 20.
        lplims(2) = 10.
        if( cline%defined('lpstop') ) lplims(2) = params%lpstop
        if( do_eo )then
            res       = spproj%os_cls2D%get_all('res')
            lplims(1) = max(median_nocopy(res), lplims(2))
            deallocate(res)
        else
            if( cline%defined('lpstart') ) lplims(1) = params%lpstart
        endif
        ! prepare a temporary project file for the class average processing
        allocate(WORK_PROJFILE, source=trim(ORIG_WORK_PROJFILE))
        call del_file(WORK_PROJFILE)
        work_proj1%projinfo = spproj%projinfo
        work_proj1%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris()  > 0 ) work_proj1%jobproc = spproj%jobproc
        call work_proj1%add_stk(trim(stk), ctfvars)
        call work_proj1%os_ptcl3D%set_all('state', real(states)) ! takes care of states
        ! name change
        call work_proj1%projinfo%delete_entry('projname')
        call work_proj1%projinfo%delete_entry('projfile')
        call cline%set('projfile', trim(WORK_PROJFILE))
        call cline%set('projname', trim(get_fbody(trim(WORK_PROJFILE),trim('simple'))))
        call work_proj1%update_projinfo(cline)
        ! split
        if(params%nparts == 1 )then
            call work_proj1%write()
        else
            call work_proj1%split_stk(params%nparts, .false.)
        endif
        ! down-scale
        orig_box     = work_proj1%get_box()
        orig_msk     = params%msk
        smpd_target  = max(params%smpd, lplims(2)*LP2SMPDFAC)
        do_autoscale = do_autoscale .and. smpd_target > work_proj1%get_smpd()
        if( do_autoscale )then
            deallocate(WORK_PROJFILE)
            call simple_mkdir(STKPARTSDIR,errmsg="commander_hlev_wflows :: exec_initial_3Dmodel;  ")
            call work_proj1%scale_projfile(smpd_target, WORK_PROJFILE, cline, cline_scale, dir=trim(STKPARTSDIR))
            scale_factor = cline_scale%get_rarg('scale')
            box          = nint(cline_scale%get_rarg('newbox'))
            msk          = cline%get_rarg('msk')
            call xscale_distr%execute( cline_scale )
        else
            box         = orig_box
            msk         = orig_msk
        endif
        ! prepare command lines from prototype
        ! projects names are subject to change depending on scaling and are updated individually
        call cline%delete('projname')
        call cline%delete('projfile')
        cline_reconstruct3D   = cline
        cline_refine3D_refine = cline
        cline_reproject       = cline
        cline_refine3D_snhc   = cline
        cline_refine3D_init   = cline
        cline_symsrch         = cline
        ! In shnc & stage 1 the objective function is always standard cross-correlation,
        ! in stage 2 it follows optional user input and defaults to ccres
        call cline_refine3D_snhc%set('objfun', 'cc')
        call cline_refine3D_init%set('objfun', 'cc')
        if( .not.cline%defined('objfun') )call cline_refine3D_refine%set('objfun', 'ccres')
        ! reconstruct3D & project are not distributed executions, so remove the nparts flag
        call cline_reconstruct3D%delete('nparts')
        call cline_reproject%delete('nparts')
        ! initialise command line parameters
        ! (1) INITIALIZATION BY STOCHASTIC NEIGHBORHOOD HILL-CLIMBING
        call cline_refine3D_snhc%set('projfile', trim(WORK_PROJFILE))
        call cline_refine3D_snhc%set('msk',      msk)
        call cline_refine3D_snhc%set('box',      real(box))
        call cline_refine3D_snhc%set('prg',    'refine3D')
        call cline_refine3D_snhc%set('refine',  'snhc')
        call cline_refine3D_snhc%set('lp',      lplims(1))
        call cline_refine3D_snhc%set('nspace',  real(NSPACE_SNHC))
        call cline_refine3D_snhc%set('maxits',  real(MAXITS_SNHC))
        call cline_refine3D_snhc%delete('update_frac') ! no fractional update in first phase
        ! (2) REFINE3D_INIT
        call cline_refine3D_init%set('prg',      'refine3D')
        call cline_refine3D_init%set('projfile', trim(WORK_PROJFILE))
        call cline_refine3D_init%set('msk',      msk)
        call cline_refine3D_init%set('box',      real(box))
        call cline_refine3D_init%set('maxits',   real(MAXITS_INIT))
        call cline_refine3D_init%set('vol1',     trim(SNHCVOL)//trim(str_state)//params%ext)
        call cline_refine3D_init%set('lp',       lplims(1))
        if( .not. cline_refine3D_init%defined('nspace') )then
            call cline_refine3D_init%set('nspace', real(NSPACE_INIT))
        endif
        ! (3) SYMMETRY AXIS SEARCH
        if( srch4symaxis )then
            ! need to replace original point-group flag with c1
            call cline_refine3D_snhc%set('pgrp', 'c1')
            call cline_refine3D_init%set('pgrp', 'c1')
            ! symsrch
            call qenv%new(1, exec_bin='simple_exec')
            call cline_symsrch%set('prg',     'symaxis_search') ! needed for cluster exec
            call cline_symsrch%set('msk',      msk)
            call cline_symsrch%set('smpd',     work_proj1%get_smpd())
            call cline_symsrch%set('projfile', trim(WORK_PROJFILE))
            call cline_symsrch%set('cenlp',    CENLP)
            call cline_symsrch%set('hp',       params%hp)
            call cline_symsrch%set('lp',       lplims(1))
            call cline_symsrch%set('oritype',  'ptcl3D')
        endif
        ! (4) REFINE3D REFINE STEP
        call cline_refine3D_refine%set('prg',      'refine3D')
        call cline_refine3D_refine%set('projfile', trim(ORIG_WORK_PROJFILE))
        call cline_refine3D_refine%set('msk',      orig_msk)
        call cline_refine3D_refine%set('box',      real(orig_box))
        call cline_refine3D_refine%set('maxits',   real(MAXITS_REFINE))
        call cline_refine3D_refine%set('refine',   'single')
        if( do_eo )then
            call cline_refine3D_refine%delete('lp')
            call cline_refine3D_refine%set('eo',         'yes')
            call cline_refine3D_refine%set('lplim_crit', 0.5)
            call cline_refine3D_refine%set('lpstop',     lplims(2))
            call cline_refine3D_refine%set('clsfrcs',   'yes')
        else
            call cline_refine3D_refine%set('lp', lplims(2))
        endif
        if( .not. cline_refine3D_refine%defined('nspace') )then
            call cline_refine3D_refine%set('nspace', real(NSPACE_REFINE))
        endif
        ! (5) RE-PROJECT VOLUME
        call cline_reproject%set('prg',    'reproject')
        call cline_reproject%set('outstk', 'reprojs'//params%ext)
        call cline_reproject%set('smpd',    params%smpd)
        call cline_reproject%set('msk',     orig_msk)
        call cline_reproject%set('box',     real(orig_box))
        ! execute commanders
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> INITIALIZATION WITH STOCHASTIC NEIGHBORHOOD HILL-CLIMBING'
        write(*,'(A,F6.1,A)') '>>> LOW-PASS LIMIT FOR ALIGNMENT: ', lplims(1),' ANGSTROMS'
        write(*,'(A)') '>>>'
        call xrefine3D_distr%execute(cline_refine3D_snhc)
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> INITIAL 3D MODEL GENERATION WITH REFINE3D'
        write(*,'(A)') '>>>'
        call xrefine3D_distr%execute(cline_refine3D_init)
        iter = cline_refine3D_init%get_rarg('endit')
        call set_iter_dependencies
        if( srch4symaxis )then
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> SYMMETRY AXIS SEARCH'
            write(*,'(A)') '>>>'
            call cline_symsrch%set('vol1', trim(vol_iter))
            if( qenv%get_qsys() .eq. 'local' )then
                call xsymsrch%execute(cline_symsrch)
            else
                call qenv%exec_simple_prg_in_queue(cline_symsrch, 'SYMAXIS_SEARCH', 'SYMAXIS_SEARCH_FINISHED')
            endif
            call del_file('SYMAXIS_SEARCH_FINISHED')
        endif
        ! prep refinement stage
        call work_proj1%read_segment('ptcl3D', trim(WORK_PROJFILE))
        os = work_proj1%os_ptcl3D
        ! modulate shifts
        if( do_autoscale )call os%mul_shifts( 1./scale_factor )
        ! clean stacks & project file
        call work_proj1%read_segment('stk', trim(WORK_PROJFILE))
        do istk=1,work_proj1%os_stk%get_noris()
            call work_proj1%os_stk%getter(istk, 'stk', stk)
            call del_file(trim(stk))
        enddo
        call work_proj1%kill()
        call del_file(trim(WORK_PROJFILE))
        deallocate(WORK_PROJFILE)
        ! re-create project
        call del_file(ORIG_WORK_PROJFILE)
        ctfvars%smpd        = params%smpd
        work_proj2%projinfo = spproj%projinfo
        work_proj2%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris()  > 0 ) work_proj2%jobproc = spproj%jobproc
        if( do_eo )then
            call prep_eo_stks_refine
            params_glob%eo     = 'yes'
            params_glob%nptcls = work_proj2%get_nptcls()
        else
            call work_proj2%add_stk(trim(orig_stk), ctfvars)
            work_proj2%os_ptcl3D = os
            call work_proj2%os_ptcl3D%set_all('state', real(states))
        endif
        call os%kill
        ! naming
        call work_proj2%projinfo%delete_entry('projname')
        call work_proj2%projinfo%delete_entry('projfile')
        call cline%set('projfile', trim(ORIG_WORK_PROJFILE))
        call cline%set('projname', trim(get_fbody(trim(ORIG_WORK_PROJFILE),trim('simple'))))
        call work_proj2%update_projinfo(cline)
        ! split
        if( do_eo )then
            call work_proj2%write()
        else
            if(params%nparts == 1)then
                call work_proj2%write()
            else
                call work_proj2%split_stk(params%nparts, .false.)
            endif
        endif
        call work_proj2%kill
        ! refinement stage
        if( do_autoscale )then
            ! refine3d_init will take care of reconstruction at original scale
        else
            if( do_eo )then
                ! refine3d_init will take care of reconstruction
            else
                call cline_refine3D_refine%set('vol1', trim(vol_iter))
            endif
        endif
        write(*,'(A)') '>>>'
        if( do_autoscale )then
            write(*,'(A)') '>>> PROBABILISTIC REFINEMENT AT ORIGINAL SAMPLING'
        else
            write(*,'(A)') '>>> PROBABILISTIC REFINEMENT'
        endif
        write(*,'(A)') '>>>'
        call cline_refine3D_refine%set('startit', iter + 1.)
        call xrefine3D_distr%execute(cline_refine3D_refine)
        iter = cline_refine3D_refine%get_rarg('endit')
        call set_iter_dependencies
        ! rename final volumes
        status = simple_rename(vol_iter, 'rec_final'//params%ext)
        status = simple_rename(add2fbody(vol_iter,params%ext,PPROC_SUFFIX), 'rec_final_pproc'//params%ext)
        ! update the original project cls3D segment with orientations
        call work_proj2%read_segment('ptcl3D', trim(ORIG_WORK_PROJFILE))
        call work_proj2%os_ptcl3D%delete_entry('stkind')
        call work_proj2%os_ptcl3D%delete_entry('eo')
        params_glob%nptcls = ncavgs
        if( do_eo )then
            ! setting the even cavgs oris here, there must be a better solution
            ! to check map2ptcls
            call spproj%os_cls3D%new(ncavgs)
            do icls=1,ncavgs
                call spproj%os_cls3D%set_ori(icls, work_proj2%os_ptcl3D%get_ori(icls))
            enddo
            call conv_eo(work_proj2%os_ptcl3D)
        else
            spproj%os_cls3D = work_proj2%os_ptcl3D
        endif
        call work_proj2%kill
        ! map the orientation parameters obtained for the clusters back to the particles
        call spproj%map2ptcls
        ! add rec_final to os_out
        call spproj%add_vol2os_out('rec_final'//params%ext, params%smpd, 1, 'vol_cavg')
        ! write results (this needs to be a full write as multiple segments are updated)
        call spproj%write()
        ! reprojections
        call spproj%os_cls3D%write('final_oris.txt')
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> RE-PROJECTION OF THE FINAL VOLUME'
        write(*,'(A)') '>>>'
        call cline_reproject%set('vol1',   'rec_final_pproc'//params%ext)
        call cline_reproject%set('oritab', 'final_oris.txt')
        call xreproject%execute(cline_reproject)
        ! end gracefully
        call spproj%kill
        call del_file(ORIG_WORK_PROJFILE)
        call simple_end('**** SIMPLE_INITIAL_3DMODEL NORMAL STOP ****')

        contains

            subroutine set_iter_dependencies
                character(len=3) :: str_iter
                str_iter = int2str_pad(nint(iter),3)
                vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//params%ext
            end subroutine set_iter_dependencies

            subroutine prep_eo_stks_refine
                use simple_ori, only: ori
                type(ori)                     :: o, o_even, o_odd
                character(len=:), allocatable :: eostk, ext
                integer :: even_ind, odd_ind, state, icls
                call os%delete_entry('lp')
                call cline_refine3D_refine%set('frcs',trim(frcs_fname))
                ! add stks
                ext   = '.'//fname2ext( stk )
                eostk = add2fbody(trim(orig_stk), trim(ext), '_even')
                call work_proj2%add_stk(eostk, ctfvars)
                eostk = add2fbody(trim(orig_stk), trim(ext), '_odd')
                call work_proj2%add_stk(eostk, ctfvars)
                ! update orientations parameters
                do icls=1,ncavgs
                    even_ind = icls
                    odd_ind  = ncavgs+icls
                    o        = os%get_ori(icls)
                    state    = os%get_state(icls)
                    call o%set('class', real(icls)) ! for mapping frcs in 3D
                    call o%set('state', real(state))
                    ! even
                    o_even = o
                    call o_even%set('eo', 0.)
                    call o_even%set('stkind', work_proj2%os_ptcl3D%get(even_ind,'stkind'))
                    call work_proj2%os_ptcl3D%set_ori(even_ind, o_even)
                    ! odd
                    o_odd = o
                    call o_odd%set('eo', 1.)
                    call o_odd%set('stkind', work_proj2%os_ptcl3D%get(odd_ind,'stkind'))
                    call work_proj2%os_ptcl3D%set_ori(odd_ind, o_odd)
                enddo
                ! cleanup
                deallocate(eostk, ext)
                call o%kill
                call o_even%kill
                call o_odd%kill
            end subroutine prep_eo_stks_refine

            subroutine conv_eo( os )
                use simple_sym, only: sym
                use simple_ori, only: ori
                class(oris), intent(inout) :: os
                type(sym) :: se
                type(ori) :: o_odd, o_even
                real      :: avg_euldist, euldist
                integer   :: icls, ncls
                call se%new(params%pgrp)
                avg_euldist = 0.
                ncls = 0
                do icls=1,os%get_noris()/2
                    o_even = os%get_ori(icls)
                    if( o_even%get_state() == 0 )cycle
                    ncls    = ncls + 1
                    o_odd   = os%get_ori(ncavgs+icls)
                    euldist = rad2deg(o_odd.euldist.o_even)
                    if( se%get_nsym() > 1 )then
                        call o_odd%mirror2d
                        call se%rot_to_asym(o_odd)
                        euldist = min(rad2deg(o_odd.euldist.o_even), euldist)
                    endif
                    avg_euldist = avg_euldist + euldist
                enddo
                avg_euldist = avg_euldist/real(ncls)
                write(*,'(A)')'>>>'
                write(*,'(A,F6.1)')'>>> EVEN/ODD AVERAGE ANGULAR DISTANCE: ', avg_euldist
            end subroutine conv_eo

    end subroutine exec_initial_3Dmodel

    !> for heterogeinity analysis
    subroutine exec_cluster3D( self, cline )
        use simple_oris,             only: oris
        use simple_sym,              only: sym
        use simple_commander_volops, only: postprocess_commander
        use simple_cluster_seed,     only: gen_labelling
        use simple_commander_distr_wflows, only: refine3D_distr_commander, reconstruct3D_distr_commander
        class(cluster3D_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        ! constants
        integer,          parameter :: MAXITS1 = 50
        integer,          parameter :: MAXITS2 = 40
        character(len=*), parameter :: one     = '01'
        ! distributed commanders
        type(refine3D_distr_commander)      :: xrefine3D_distr
        type(reconstruct3D_distr_commander) :: xreconstruct3D_distr
        ! command lines
        type(cmdline) :: cline_refine3D1, cline_refine3D2
        type(cmdline) :: cline_reconstruct3D_distr, cline_reconstruct3D_mixed_distr
        ! other variables
        type(parameters)      :: params
        type(sym)             :: symop
        type(sp_project)      :: spproj
        type(oris)            :: os
        integer,  allocatable :: labels(:)
        real                  :: trs
        integer               :: iter, startit, rename_stat
        logical               :: write_proj
        ! sanity check
        if(nint(cline%get_rarg('nstates')) <= 1)stop 'Non-sensical NSTATES argument for heterogeneity analysis!'
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! make master parameters
        call params%new(cline)
        if( params%eo .eq. 'no' .and. .not. cline%defined('lp') )stop 'need lp input when eo .eq. no; cluster3D'
        ! set mkdir to no
        call cline%set('mkdir', 'no')
        ! prepare command lines from prototype
        call cline%delete('refine')
        if(.not.cline%defined('lplim_crit'))call cline%set('lplim_crit', 0.5)
        cline_refine3D1                 = cline ! first stage, extremal optimization
        cline_refine3D2                 = cline ! second stage, stochastic refinement
        cline_reconstruct3D_distr       = cline ! eo always eq yes, for resolution only
        cline_reconstruct3D_mixed_distr = cline ! eo always eq yes, for resolution only
        ! first stage
        call cline_refine3D1%set('prg', 'refine3D')
        call cline_refine3D1%set('maxits', real(MAXITS1))
        call cline_refine3D1%delete('neigh')
        select case(trim(params%refine))
            case('sym')
                call cline_refine3D1%set('refine', 'clustersym')
            case('cluster_snhc')
                call cline_refine3D1%set('refine', 'cluster_snhc')
                call cline_refine3D1%set('neigh', 'yes')
                call cline_refine3D1%set('nnn',   10.)
            case DEFAULT
                call cline_refine3D1%set('refine', trim(params%refine))
        end select
        call cline_refine3D1%delete('update_frac')  ! no update frac for extremal optimization
        ! second stage
        call cline_refine3D2%set('prg',    'refine3D')
        call cline_refine3D2%set('refine', 'multi')
        if(.not.cline%defined('neigh'))then
            call cline_refine3D2%set('neigh',  'yes')
            call cline_refine3D2%set('nnn',    real(max(30,nint(0.05*real(params%nspace)))))
        endif
        if( .not.cline%defined('update_frac') )call cline_refine3D2%set('update_frac', 0.5)
        ! reconstructions
        call cline_reconstruct3D_distr%set('prg', 'reconstruct3D')
        call cline_reconstruct3D_distr%delete('lp')
        call cline_reconstruct3D_distr%set('eo','yes')
        call cline_reconstruct3D_mixed_distr%set('prg', 'reconstruct3D')
        call cline_reconstruct3D_mixed_distr%delete('lp')
        call cline_reconstruct3D_mixed_distr%set('nstates', 1.)
        if( trim(params%refine) .eq. 'sym' ) call cline_reconstruct3D_mixed_distr%set('pgrp', 'c1')
        if( cline%defined('trs') )then
            ! all good
        else
            ! works out shift limits for in-plane search
            trs = MSK_FRAC*real(params%msk)
            trs = min(MAXSHIFT, max(MINSHIFT, trs))
            call cline_refine3D1%set('trs',trs)
            call cline_refine3D2%set('trs',trs)
        endif

        ! PREP
        call spproj%read(params%projfile )
        ! splitting
        call spproj%split_stk(params%nparts, (params%mkdir.eq.'yes'), dir=PATH_PARENT)
        ! wipe previous states
        os     = spproj%os_ptcl3D
        labels = nint(os%get_all('states'))
        if( any(labels > 1) )then
            where(labels > 0) labels = 1
            call os%set_all('state', real(labels))
        endif
        deallocate(labels)
        ! e/o partition
        if( params%eo.eq.'no' )then
            call os%set_all2single('eo', -1.)
        else
            if( os%get_nevenodd() == 0 ) call os%partition_eo
        endif
        if( trim(params%refine) .eq. 'sym' )then
            ! randomize projection directions with respect to symmetry
            symop = sym(params%pgrp)
            call symop%symrandomize(os)
            call symop%kill
        endif

        ! MIXED MODEL RECONSTRUCTION
        ! retrieve mixed model Fourier components, normalization matrix, FSC & anisotropic filter
        if( params%eo .ne. 'no' )then
            spproj%os_ptcl3D = os
            call spproj%write
            call xreconstruct3D_distr%execute( cline_reconstruct3D_mixed_distr )
            rename_stat = simple_rename(trim(VOL_FBODY)//one//params%ext, trim(CLUSTER3D_VOL)//params%ext)
            rename_stat = simple_rename(trim(VOL_FBODY)//one//'_even'//params%ext, trim(CLUSTER3D_VOL)//'_even'//params%ext)
            rename_stat = simple_rename(trim(VOL_FBODY)//one//'_odd'//params%ext,  trim(CLUSTER3D_VOL)//'_odd'//params%ext)
            rename_stat = simple_rename(trim(FSC_FBODY)//one//BIN_EXT, trim(CLUSTER3D_FSC))
            rename_stat = simple_rename(trim(FRCS_FBODY)//one//BIN_EXT, trim(CLUSTER3D_FRCS))
            rename_stat = simple_rename(trim(ANISOLP_FBODY)//one//params%ext, trim(CLUSTER3D_ANISOLP)//params%ext)
        endif

        ! randomize state labels
        write(*,'(A)') '>>>'
        call gen_labelling(os, params%nstates, 'squared_uniform')
        call os%write('cluster3D_init.txt') ! analysis purpose only
        ! writes for refine3D
        spproj%os_ptcl3D = os
        call spproj%write
        call spproj%kill
        call os%kill

        ! STAGE1: extremal optimization, frozen orientation parameters
        write(*,'(A)')    '>>>'
        write(*,'(A,I3)') '>>> 3D CLUSTERING - STAGE 1'
        write(*,'(A)')    '>>>'
        call xrefine3D_distr%execute(cline_refine3D1)
        iter = nint(cline_refine3D1%get_rarg('endit'))
        ! for analysis purpose only
        call spproj%read_segment(params%oritype, params%projfile)
        call spproj%os_ptcl3D%write('cluster3D_stage1.txt')
        call spproj%kill

        ! STAGE2: soft multi-states refinement
        startit = iter + 1
        call cline_refine3D2%set('startit', real(startit))
        call cline_refine3D2%set('maxits',  real(startit + MAXITS2))
        write(*,'(A)')    '>>>'
        write(*,'(A,I3)') '>>> 3D CLUSTERING - STAGE 2'
        write(*,'(A)')    '>>>'
        call xrefine3D_distr%execute(cline_refine3D2)

        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER3D NORMAL STOP ****')
    end subroutine exec_cluster3D

    !> multi-particle refinement after cluster3D
    subroutine exec_cluster3D_refine( self, cline )
        use simple_oris,                   only: oris
        use simple_parameters,             only: params_glob
        use simple_commander_distr_wflows, only: refine3D_distr_commander
        class(cluster3D_refine_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! constants
        integer,                 parameter :: MAXITS = 40
        ! distributed commanders
        type(refine3D_distr_commander)     :: xrefine3D_distr
        ! command lines
        type(cmdline),         allocatable :: cline_refine3D(:)
        ! other variables
        integer,               allocatable :: state_pops(:), states(:), master_states(:)
        character(len=STDLEN), allocatable :: dirs(:), projfiles(:)
        character(len=:),      allocatable :: projname
        type(parameters)         :: params
        type(sp_project)         :: spproj, spproj_master
        integer                  :: state, iptcl, cnt, nstates, single_state
        logical                  :: l_singlestate, error
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call params%new(cline)
        l_singlestate = cline%defined('state')
        if( l_singlestate )then
            single_state = nint(cline%get_rarg('state'))
        else
            single_state = 0
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')

        ! sanity checks
        if( params%eo .eq. 'no' .and. .not. cline%defined('lp') )&
            &stop 'need lp input when eo .eq. no; cluster3D_refine'
        if( .not.cline%defined('maxits') )call cline%set('maxits',real(MAXITS))

        ! prep individual project files
        call spproj_master%read(params%projfile)
        master_states  = nint(spproj_master%os_ptcl3D%get_all('state'))
        nstates        = maxval(master_states)
        params%nstates = nstates
        if( params%nstates==1 )then
            write(*,*) 'Non-sensical number of states for heterogeneity refinement: ',params%nstates
            stop 'Non-sensical number of states'
        endif
        call spproj_master%os_ptcl3D%get_pops(state_pops, 'state', consider_w=.false.)
        if( state_pops(params%state) == 0 )then
            write(*,*) 'Empty state to refine: ', params%state
            stop 'Empty state to refine'
        endif
        allocate(projfiles(params%nstates), dirs(params%nstates), cline_refine3D(params%nstates))
        ! states are lost from the project after this loop and stored in master_states
        cnt = 0
        do state = 1, params%nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. single_state.ne.state )cycle
            cnt = cnt + 1
            ! name & directory
            projname         = 'state_'//trim(int2str_pad(state,2))
            projfiles(state) = trim(projname)//trim(METADATA_EXT)
            dirs(state)      = trim(int2str(state))//'_refine3D'
            ! command line
            cline_refine3D(state) = cline
            call cline_refine3D(state)%set('prg',     'refine3D')
            call cline_refine3D(state)%set('projname',trim(projname))
            call cline_refine3D(state)%set('projfile',trim(projfiles(state)))
            call cline_refine3D(state)%set('mkdir',   'yes')
            call cline_refine3D(state)%set('refine',  'single')
            call cline_refine3D(state)%delete('state')
            call cline_refine3D(state)%delete('nstates')
            ! states
            states = master_states
            where(states /= state) states = 0
            where(states /= 0)     states = 1
            call spproj_master%os_ptcl3D%set_all('state', real(states))
            ! write
            call spproj_master%update_projinfo(cline_refine3D(state))
            call spproj_master%write(projfiles(state))
            deallocate(states)
        enddo
        call spproj_master%update_projinfo(cline) ! restores name

        ! Execute individual refine3D jobs
        do state = 1, nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. state.ne.single_state )cycle
            write(*,'(A)')   '>>>'
            write(*,'(A,I2,A,A)')'>>> REFINING STATE: ', state, ' IN DIRECTORY: ', trim(dirs(state))
            write(*,'(A)')   '>>>'
            params_glob%projname = 'state_'//trim(int2str_pad(state,2))
            params_glob%projfile = projfiles(state)
            params_glob%nstates = 1
            params_glob%state   = 1
            call xrefine3D_distr%execute(cline_refine3D(state))
            call simple_chdir(PATH_PARENT,errmsg="commander_hlev_wflows :: exec_cluster3D_refine;")
            ! renames volumes and updates in os_out
            call stash_state(state)
        enddo
        ! restores original values
        params_glob%projname = params%projname
        params_glob%projfile = params%projfile
        params_glob%nstates  = nstates

        ! consolidates new orientations parameters
        do state = 1, nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. state.ne.single_state )cycle
            ! transfer orientations
            call spproj%read_segment(params%oritype, filepath(dirs(state),projfiles(state)))
            do iptcl=1,params%nptcls
                if( master_states(iptcl)==state )then
                    call spproj_master%os_ptcl3D%set_ori(iptcl, spproj%os_ptcl3D%get_ori(iptcl))
                endif
            enddo
            call spproj%kill
        enddo
        ! reset states
        call spproj_master%os_ptcl3D%set_all('state', real(master_states)) ! restores original states
        call spproj_master%write

        ! cleanup
        call spproj_master%kill
        deallocate(master_states, dirs, projfiles)
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER3D_REFINE NORMAL STOP ****')

        contains

            ! stash docs, volumes , etc.
            subroutine stash_state(s)
                integer, intent(in)   :: s
                character(len=2), parameter :: one = '01'
                character(len=STDLEN) :: src, dest
                character(len=2)      :: str_state
                character(len=8)      :: str_iter
                integer               :: it, final_it, stat
                final_it  = nint(cline_refine3D(s)%get_rarg('endit'))
                str_state = int2str_pad(s,2)
                do it = 1,final_it
                    str_iter = '_iter'//int2str_pad(it,3)
                    ! volume
                    src  = filepath( dirs(s),trim(VOL_FBODY)//one//str_iter//params%ext)
                    dest = trim(VOL_FBODY)//str_state//str_iter//params%ext
                    stat = simple_rename(src, dest)
                    ! post_processed volume
                    src  = filepath( dirs(s), trim(VOL_FBODY)//one//str_iter//trim(PPROC_SUFFIX)//params%ext)
                    dest = trim(VOL_FBODY)//str_state//str_iter//trim(PPROC_SUFFIX)//params%ext
                    stat = simple_rename(src, dest)
                    ! e/o
                    if( params%eo.ne.'no')then
                        ! e/o
                        src  = filepath( dirs(s), trim(VOL_FBODY)//one//str_iter//'_even'//params%ext)
                        dest = trim(VOL_FBODY)//str_state//str_iter//'_even'//params%ext
                        stat = simple_rename(src, dest)
                        src  = filepath( dirs(s), trim(VOL_FBODY)//one//str_iter//'_odd'//params%ext)
                        dest = trim(VOL_FBODY)//str_state//str_iter//'_odd'//params%ext
                        stat = simple_rename(src, dest)
                        ! FSC
                        str_iter = '_ITER'//int2str_pad(it,3)
                        src  = filepath( dirs(s), 'RESOLUTION_STATE'//one//str_iter)
                        dest = 'RESOLUTION_STATE'//str_state//str_iter
                        stat = simple_rename(src, dest)
                    endif
                enddo
                ! updates os_out
                str_iter = '_iter'//int2str_pad(final_it,3)
                src = trim(VOL_FBODY)//str_state//str_iter//params%ext
                call spproj_master%add_vol2os_out(trim(src), params%smpd, s, 'vol')
                if( params%eo.ne.'no')then
                    src  = filepath( dirs(s), trim(FSC_FBODY)//one//BIN_EXT)
                    dest = trim(FSC_FBODY)//str_state//BIN_EXT
                    stat = simple_rename(src, dest)
                    call spproj_master%add_fsc2os_out(trim(dest), s, params%box)
                    src  = filepath( dirs(s), trim(ANISOLP_FBODY)//one//params%ext)
                    dest = trim(ANISOLP_FBODY)//str_state//params%ext
                    stat = simple_rename(src, dest)
                    call  spproj_master%add_vol2os_out(trim(dest), params%smpd, s, 'vol_filt')
                endif
            end subroutine stash_state

    end subroutine exec_cluster3D_refine

end module simple_commander_hlev_wflows
