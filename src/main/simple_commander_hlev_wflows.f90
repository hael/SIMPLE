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
        ! delete any previous solution
        if( .not. spproj%is_virgin_field(params%oritype) )then
            ! removes previous cluster2D solution (states are preserved)
            call spproj%os_ptcl2D%delete_2Dclustering
            call spproj%write_segment_inside(params%oritype)
        endif
        ! stack splitting
        call spproj%split_stk(params%nparts, dir='..')
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
                call simple_mkdir(STKPARTSDIR)
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
            if( scaling )call xscale_distr%execute( cline_scale2 )
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
        use simple_commander_volops,       only: reproject_commander, symsrch_commander
        use simple_commander_rec,          only: reconstruct3D_commander
        use simple_parameters,             only: params_glob
        use simple_qsys_env,               only: qsys_env
        class(initial_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! constants
        real,    parameter :: CENLP=30. !< consistency with refine3D
        integer, parameter :: MAXITS_SNHC=30, MAXITS_INIT=15, MAXITS_REFINE=40
        integer, parameter :: STATE=1, MAXITS_SNHC_RESTART=3
        integer, parameter :: NSPACE_SNHC=1000, NSPACE_INIT=1000, NSPACE_REFINE= 2500
        integer, parameter :: NRESTARTS=5
        integer, parameter :: NPEAKS_INIT=6, NPEAKS_REFINE=1
        character(len=STDLEN), parameter :: ORIG_WORK_PROJFILE = 'initial_3Dmodel_tmpproj.simple'
        ! distributed commanders
        type(refine3D_distr_commander)      :: xrefine3D_distr
        type(scale_project_distr_commander) :: xscale_distr
        ! shared-mem commanders
        type(symsrch_commander)       :: xsymsrch
        type(reconstruct3D_commander) :: xreconstruct3D
        type(reproject_commander)     :: xreproject
        ! command lines
        type(cmdline) :: cline_refine3D_snhc_restart
        type(cmdline) :: cline_refine3D_snhc
        type(cmdline) :: cline_refine3D_init
        type(cmdline) :: cline_refine3D_refine
        type(cmdline) :: cline_symsrch
        type(cmdline) :: cline_reconstruct3D
        type(cmdline) :: cline_reproject
        type(cmdline) :: cline_scale
        ! other variables
        character(len=:), allocatable :: stk, orig_stk
        character(len=:), allocatable :: WORK_PROJFILE
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
        ! set lplims
        lplims(1) = 20.
        lplims(2) = 10.
        if( cline%defined('lpstart') ) lplims(1) = params%lpstart
        if( cline%defined('lpstop')  ) lplims(2) = params%lpstop
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
        params%smpd = ctfvars%smpd
        orig_stk    = stk
        if( do_eo )call prep_eo_stks_init
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
            call work_proj1%split_stk(params%nparts)
        endif
        ! down-scale
        orig_box     = work_proj1%get_box()
        orig_msk     = params%msk
        smpd_target  = max(params%smpd, lplims(2)*LP2SMPDFAC)
        do_autoscale = do_autoscale .and. smpd_target > work_proj1%get_smpd()
        if( do_autoscale )then
            deallocate(WORK_PROJFILE)
            call simple_mkdir(STKPARTSDIR)
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
        cline_reconstruct3D         = cline
        cline_refine3D_refine       = cline
        cline_reproject             = cline
        cline_refine3D_snhc_restart = cline
        cline_refine3D_init         = cline
        cline_symsrch               = cline
        ! In shnc & stage 1 the objective function is always standard cross-correlation,
        ! in stage 2 it follows optional user input and defaults to ccres
        call cline_refine3D_snhc_restart%set('objfun', 'cc')
        call cline_refine3D_init%set('objfun', 'cc')
        if( .not.cline%defined('objfun') )call cline_refine3D_refine%set('objfun', 'ccres')
        ! reconstruct3D & project are not distributed executions, so remove the nparts flag
        call cline_reconstruct3D%delete('nparts')
        call cline_reproject%delete('nparts')
        ! initialise command line parameters
        ! (1) INITIALIZATION BY STOCHASTIC NEIGHBORHOOD HILL-CLIMBING
        call cline_refine3D_snhc_restart%set('projfile', trim(WORK_PROJFILE))
        call cline_refine3D_snhc_restart%set('msk',      msk)
        call cline_refine3D_snhc_restart%set('box',      real(box))
        call cline_refine3D_snhc_restart%delete('update_frac') ! no fractional update in first phase
        call cline_refine3D_snhc_restart%set('prg',    'refine3D')
        call cline_refine3D_snhc_restart%set('maxits',  real(MAXITS_SNHC_RESTART))
        call cline_refine3D_snhc_restart%set('refine',  'snhc')
        call cline_refine3D_snhc_restart%set('lp',      lplims(1))
        call cline_refine3D_snhc_restart%set('nspace',  real(NSPACE_SNHC))
        cline_refine3D_snhc = cline_refine3D_snhc_restart
        call cline_refine3D_snhc%set('maxits', real(MAXITS_SNHC))
        ! (2) REFINE3D_INIT
        call cline_refine3D_init%set('prg',      'refine3D')
        call cline_refine3D_init%set('projfile', trim(WORK_PROJFILE))
        call cline_refine3D_init%set('msk',      msk)
        call cline_refine3D_init%set('box',      real(box))
        call cline_refine3D_init%set('maxits',   real(MAXITS_INIT))
        call cline_refine3D_init%set('vol1',     trim(SNHCVOL)//trim(str_state)//params%ext)
        call cline_refine3D_init%set('lp',       lplims(1))
        call cline_refine3D_init%set('npeaks',   real(NPEAKS_INIT))
        if( .not. cline_refine3D_init%defined('nspace') )then
            call cline_refine3D_init%set('nspace', real(NSPACE_INIT))
        endif
        ! (3) SYMMETRY AXIS SEARCH
        if( srch4symaxis )then
            ! need to replace original point-group flag with c1
            call cline_refine3D_snhc%set('pgrp', 'c1')
            call cline_refine3D_init%set('pgrp', 'c1')
            ! symsrch
            call qenv%new(exec_bin='simple_exec')
            call cline_symsrch%set('msk', msk)
            call cline_symsrch%set('smpd', work_proj1%get_smpd())
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
        call cline_refine3D_refine%set('npeaks',   real(NPEAKS_REFINE))
        if( do_eo )then
            call cline_refine3D_refine%delete('lp')
            call cline_refine3D_refine%set('eo', 'yes')
            call cline_refine3D_refine%set('lplim_crit', 0.5)
            call cline_refine3D_refine%set('lpstop',lplims(2))
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
                call qenv%exec_simple_prg_in_queue(cline_symsrch, 'SYMSRCH', 'SYMSRCH_FINISHED')
            endif
            call del_file('SYMSRCH_FINISHED')
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
        endif
        call work_proj2%os_ptcl3D%set_all('state', real(states)) ! takes care of states
        ! naming
        call work_proj2%projinfo%delete_entry('projname')
        call work_proj2%projinfo%delete_entry('projfile')
        call cline%set('projfile', trim(ORIG_WORK_PROJFILE))
        call cline%set('projname', trim(get_fbody(trim(ORIG_WORK_PROJFILE),trim('simple'))))
        call work_proj2%update_projinfo(cline)
        ! split
        if( do_eo )then
            if(params%nparts <= 2)then
                call work_proj2%write()
            else
                call work_proj2%split_stk(params%nparts)
            endif
        else
            if(params%nparts == 1)then
                call work_proj2%write()
            else
                call work_proj2%split_stk(params%nparts)
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
        call cline_reproject%set('vol1',     'rec_final'//params%ext)
        call cline_reproject%set('oritab',   'final_oris.txt')
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

            subroutine prep_eo_stks_init
                use simple_commander_imgproc, only: filter_commander
                type(filter_commander) :: xfilter
                type(cmdline)          :: cline_filter
                character(len=:), allocatable :: frcs_fname, ext, stk_filt
                ext = '.'//fname2ext( stk )
                call spproj%get_frcs(frcs_fname, 'frc2D')
                if( trim(frcs_fname).eq.'' )then
                    write(*,*)'The project file does not contain the required information for e/o alignment'
                    stop 'The project file does not contain the required information for e/o alignment'
                endif
                stk_filt = add2fbody(trim(stk), trim(ext), '_filt')
                call cline_filter%set('prg',  'filter')
                call cline_filter%set('nthr', real(params%nthr))
                call cline_filter%set('smpd', params%smpd)
                call cline_filter%set('frcs', trim(frcs_fname))
                call cline_filter%set('stk',   trim(stk))
                call cline_filter%set('outstk',trim(stk_filt))
                call xfilter%execute(cline_filter)
                deallocate(stk)
                stk = trim(stk_filt)
            end subroutine prep_eo_stks_init

            subroutine prep_eo_stks_refine
                character(len=:), allocatable :: eostk, ext
                ext = '.'//fname2ext( stk )
                call os%delete_entry('lp')
                ! even
                call os%set_all2single('eo', 0.)
                call os%set_all2single('stkind', 1.)
                eostk = add2fbody(trim(orig_stk), trim(ext), '_even')
                call work_proj2%add_stk(eostk, ctfvars)
                do icls=1,ncavgs
                    call work_proj2%os_ptcl3D%set_ori(icls, os%get_ori(icls))
                enddo
                deallocate(eostk)
                ! odd
                call os%set_all2single('eo', 1.)
                call os%set_all2single('stkind', 2.)
                eostk = add2fbody(trim(orig_stk), trim(ext), '_odd')
                call work_proj2%add_stk(eostk, ctfvars)
                do icls=1,ncavgs
                    call work_proj2%os_ptcl3D%set_ori(ncavgs+icls, os%get_ori(icls))
                enddo
            end subroutine prep_eo_stks_refine

    end subroutine exec_initial_3Dmodel

    !> for heterogeinity analysis
    subroutine exec_cluster3D( self, cline )
        use simple_oris,             only: oris
        use simple_sym,              only: sym
        use simple_commander_volops, only: postprocess_commander
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
        ! shared-mem commanders
        type(postprocess_commander)         :: xpostprocess
        ! command lines
        type(cmdline) :: cline_refine3D1, cline_refine3D2
        type(cmdline) :: cline_reconstruct3D_distr, cline_reconstruct3D_mixed_distr, cline_postprocess
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
        cline_refine3D1                 = cline ! first stage, extremal optimization
        cline_refine3D2                 = cline ! second stage, stochastic refinement
        cline_postprocess               = cline ! eo always eq yes, for resolution only
        cline_reconstruct3D_distr       = cline ! eo always eq yes, for resolution only
        cline_reconstruct3D_mixed_distr = cline ! eo always eq yes, for resolution only
        ! first stage
        call cline_refine3D1%set('prg', 'refine3D')
        call cline_refine3D1%set('maxits', real(MAXITS1))
        select case(trim(params%refine))
            case('sym')
                call cline_refine3D1%set('refine', 'clustersym')
            case DEFAULT
                if( .not.cline%defined('refine') )then
                    call cline_refine3D1%set('refine', 'cluster')
                else
                    call cline_refine3D1%set('refine', trim(params%refine))
                endif
        end select
        call cline_refine3D1%delete('neigh')
        call cline_refine3D1%delete('update_frac')  ! no update frac for extremal optimization
        ! second stage
        call cline_refine3D2%set('prg', 'refine3D')
        call cline_refine3D2%set('refine', 'multi')
        if( .not.cline%defined('update_frac') )call cline_refine3D2%set('update_frac', 0.5)
        ! reconstructions
        call cline_reconstruct3D_distr%set('prg', 'reconstruct3D')
        call cline_reconstruct3D_distr%delete('lp')
        call cline_reconstruct3D_distr%set('eo','yes')
        call cline_reconstruct3D_mixed_distr%set('prg', 'reconstruct3D')
        call cline_reconstruct3D_mixed_distr%delete('lp')
        call cline_reconstruct3D_mixed_distr%set('nstates', 1.)
        if( trim(params%refine) .eq. 'sym' ) call cline_reconstruct3D_mixed_distr%set('pgrp', 'c1')
        call cline_postprocess%set('prg', 'postprocess')
        call cline_postprocess%delete('lp')
        call cline_postprocess%set('eo','yes')
        call cline_postprocess%set('mkdir','no')
        if( cline%defined('trs') )then
            ! all good
        else
            ! works out shift lmits for in-plane search
            trs = MSK_FRAC*real(params%msk)
            trs = min(MAXSHIFT, max(MINSHIFT, trs))
            call cline_refine3D1%set('trs',trs)
            call cline_refine3D2%set('trs',trs)
        endif

        ! MIXED MODEL RECONSTRUCTION
        ! retrieve mixed model Fourier components, normalization matrix, FSC & anisotropic filter
        if( params%eo .ne. 'no' )then
            if( trim(params%refine) .eq. 'sym' ) call cline_reconstruct3D_mixed_distr%set('pgrp', 'c1')
            call xreconstruct3D_distr%execute( cline_reconstruct3D_mixed_distr )
            rename_stat = simple_rename(trim(VOL_FBODY)//one//params%ext, trim(CLUSTER3D_VOL)//params%ext)
            rename_stat = simple_rename(trim(VOL_FBODY)//one//'_even'//params%ext, trim(CLUSTER3D_VOL)//'_even'//params%ext)
            rename_stat = simple_rename(trim(VOL_FBODY)//one//'_odd'//params%ext,  trim(CLUSTER3D_VOL)//'_odd'//params%ext)
            rename_stat = simple_rename(trim(FSC_FBODY)//one//BIN_EXT, trim(CLUSTER3D_FSC))
            rename_stat = simple_rename(trim(FRCS_FBODY)//one//BIN_EXT, trim(CLUSTER3D_FRCS))
            rename_stat = simple_rename(trim(ANISOLP_FBODY)//one//params%ext, trim(CLUSTER3D_ANISOLP)//params%ext)
        endif

        ! PREP
        call spproj%read(params%projfile )
        os = spproj%os_ptcl3D
        ! wipe previous states
        labels = os%get_all('states')
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
        write(*,'(A)') '>>> GENERATING DIVERSE LABELING'
        call diverse_labeling(os, params%nstates, labels, corr_ranked=.true.)
        call os%set_all('state', real(labels))
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
        iter   = nint(cline_refine3D2%get_rarg('endit'))
        ! stage 2 reconstruction to obtain resolution estimate when eo .eq. 'no'
        ! if( params%eo .eq. 'no' )then
        !     call xreconstruct3D_distr%execute(cline_reconstruct3D_distr)
        !     do state = 1, params%nstates
        !         str_state  = int2str_pad(state, 2)
        !         call cline_postprocess%set('state', real(state))
        !         call cline_postprocess%set('fsc', trim(FSC_FBODY)//trim(str_state)//BIN_EXT)
        !         call cline_postprocess%set('vol_filt', trim(ANISOLP_FBODY)//trim(str_state)//params%ext)
        !         call xpostprocess%execute(cline_postprocess)
        !     enddo
        ! endif

        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER3D NORMAL STOP ****')

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
                    if(alloc_stat /= 0) call allocchk('In: commander_hlev_wflows::diverse_labeling ', alloc_stat)
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
                    allocate(config_diverse(nptcls), source=0, stat=alloc_stat )
                    allocate(order(nptcls),          source=0, stat=alloc_stat )
                    allocate(tmp(nlabels),           source=0, stat=alloc_stat )
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

    end subroutine exec_cluster3D

    !> multi-particle refinement after cluster3D
    subroutine exec_cluster3D_refine( self, cline )
        use simple_binoris_io,       only: binread_nlines, binread_oritab, binwrite_oritab
        use simple_oris,             only: oris
        use simple_commander_volops, only: postprocess_commander
        use simple_commander_distr_wflows, only: refine3D_distr_commander, reconstruct3D_distr_commander
        class(cluster3D_refine_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! constants
        integer,            parameter :: MAXITS = 40
        character(len=:), allocatable :: INIT_FBODY
        character(len=:), allocatable :: FINAL_FBODY
        character(len=:), allocatable :: FINAL_DOC
        ! distributed commanders
        type(refine3D_distr_commander)      :: xrefine3D_distr
        type(reconstruct3D_distr_commander) :: xreconstruct3D_distr
        ! shared-mem commanders
        type(postprocess_commander)   :: xpostprocess
        ! command lines
        type(cmdline)                 :: cline_refine3D_master
        type(cmdline)                 :: cline_refine3D
        type(cmdline)                 :: cline_reconstruct3D_distr
        type(cmdline)                 :: cline_postprocess
        ! other variables
        integer,               allocatable :: state_pops(:), states(:)
        character(len=STDLEN), allocatable :: init_docs(:), final_docs(:)
        logical,               allocatable :: l_hasmskvols(:), l_hasvols(:)
        type(parameters)         :: params
        type(sp_project)         :: spproj_master
        type(sp_project), target :: spproj_state
        class(oris), pointer     :: os_master => null(), os_state => null()
        character(len=STDLEN)    :: oritab, vol, fname
        character(len=9)         :: dir
        character(len=2)         :: str_state
        integer                  :: state, iptcl, iter
        logical                  :: l_singlestate, error
        integer                  :: rename_stat
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call params%new(cline)
        l_singlestate = cline%defined('state')
        error         = .false.
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')

        ! filenaming strings allocation
        allocate(INIT_FBODY, source='cluster3Dinit_refine_state')
        allocate(FINAL_FBODY, source='cluster3Ddoc_refine_state')
        allocate(FINAL_DOC, source='cluster3Ddoc_refine'//trim(METADATA_EXT))

        ! sanity checks
        if( params%eo .eq. 'no' .and. .not. cline%defined('lp') )&
            &stop 'need lp input when eo .eq. no; cluster3D_refine'
        if(.not.file_exists(params%oritab))then
            print *,'Document ',trim(params%oritab),' does not exist!'
            stop
        endif

        ! general prep
        params%nptcls = binread_nlines( params%oritab)
        call spproj_master%new_seg_with_ptr(params%nptcls, params%oritype, os_master)
        call binread_oritab(params%oritab, spproj_master, os_master, [1,params%nptcls])
        params%nstates = os_master%get_n('state')
        if(params%nstates < 2 .and. .not.l_singlestate)then
            print *, 'Non-sensical number of states argument for heterogeneity refinement: ',params%nstates
            stop
        endif
        allocate(l_hasmskvols(params%nstates), source=.false.)
        allocate(l_hasvols(params%nstates),    source=.false.)
        allocate(init_docs(params%nstates), final_docs(params%nstates))
        call os_master%get_pops(state_pops, 'state', consider_w=.false.)
        do state = 1, params%nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. state.ne.params%state )cycle
            str_state = int2str_pad(state,2)
            dir = 'state_'//str_state//'/'
            call simple_mkdir(dir)
            ! preps individual documents
            spproj_state%os_ptcl3D = spproj_master%os_ptcl3D
            os_state => spproj_state%os_ptcl3D
            states   = nint(os_state%get_all('state'))
            where( states .ne. state )
                states = 0
            else where
                states = 1
            end where
            call os_state%set_all('state', real(states))
            deallocate(states)
            init_docs(state) = trim(INIT_FBODY)//str_state//trim(METADATA_EXT)
            call binwrite_oritab(init_docs(state), spproj_state, os_state, [1,params%nptcls])
            final_docs(state) = trim(FINAL_FBODY)//str_state//trim(METADATA_EXT)
            ! check & move volumes
            l_hasvols(state) = trim(params%vols(state)) .ne. ''
            if( l_hasvols(state) )then
                if( params%eo .ne. 'no' )then
                    ! fsc
                    fname = trim(FSC_FBODY)//str_state//BIN_EXT
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                        error = .true.
                    else
                        rename_stat = simple_rename(fname, dir//trim(fname))
                    endif
                    ! FRC
                    fname = trim(FRCS_FBODY)//str_state//BIN_EXT
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                    else
                        rename_stat = simple_rename(fname, dir//trim(fname))
                    endif
                    ! aniso
                    fname = trim(ANISOLP_FBODY)//str_state//params%ext
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                    else
                        rename_stat = simple_rename(fname, dir//trim(fname))
                    endif
                    ! e/o
                    fname = add2fbody(trim(params%vols(state)), params%ext, '_even')
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                    else
                        rename_stat = simple_rename(fname, dir//trim(fname))
                    endif
                    fname = add2fbody(trim(params%vols(state)), params%ext, '_odd')
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                    else
                        rename_stat = simple_rename(fname, dir//trim(fname))
                    endif
                endif
                ! volume
                if( .not.file_exists(params%vols(state)) )then
                    print *, 'File missing: ', params%vols(state)
                    error = .true.
                else
                    fname = trim(params%vols(state))
                    params%vols(state) = dir//trim(fname) ! name change
                    rename_stat = simple_rename(fname, params%vols(state))
                endif
            endif
            ! mask volume
            l_hasmskvols(state) = trim(params%mskvols(state)) .ne. ''
            if( l_hasmskvols(state) )then
                if( .not.file_exists(params%mskvols(state)) )then
                    print *, 'File missing: ', trim(fname)
                    error = .true.
                else
                    fname = trim(params%mskvols(state))
                    params%mskvols(state) = dir//trim(fname)  ! name change
                    call simple_copy_file(trim(fname), trim(params%mskvols(state)))
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
        do state = 1, params%nstates
            call cline%delete('vol'//int2str_pad(state,1))
        enddo
        cline_refine3D_master     = cline
        cline_reconstruct3D_distr = cline
        cline_postprocess         = cline
        call cline_refine3D_master%set('prg', 'refine3D')
        if( .not.cline%defined('maxits') ) call cline_refine3D_master%set('maxits', real(MAXITS))
        call cline_reconstruct3D_distr%set('prg', 'reconstruct3D')
        call cline_reconstruct3D_distr%set('oritab', trim(FINAL_DOC))
        call cline_reconstruct3D_distr%set('nstates', trim(int2str(params%nstates)))
        call cline_reconstruct3D_distr%set('eo', 'yes')
        call cline_postprocess%set('eo', 'yes')
        call cline_postprocess%delete('lp')

        ! Main loop
        do state = 1, params%nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. state.ne.params%state )cycle
            str_state = int2str_pad(state,2)
            dir       = 'state_'//str_state//'/'
            write(*,'(A)')   '>>>'
            write(*,'(A,I2)')'>>> REFINING STATE: ', state
            write(*,'(A)')   '>>>'
            ! prep
            cline_refine3D = cline_refine3D_master
            call cline_refine3D%set('oritab', trim(init_docs(state)))
            if( l_hasvols(state) )then
                call cline_refine3D%set('vol1', trim(params%vols(state)))
                if( params%eo.ne.'no' )then
                    fname = dir//trim(FSC_FBODY)//str_state//BIN_EXT
                    call simple_copy_file(trim(fname),'fsc_state01.bin')
                    fname = dir//trim(FRCS_FBODY)//str_state//BIN_EXT
                    if(file_exists(fname))call simple_copy_file(trim(fname),'frcs_state01.bin')
                    fname = dir//trim(ANISOLP_FBODY)//str_state//params%ext
                    if(file_exists(fname))call simple_copy_file(trim(fname),'aniso_optlp_state01.mrc')
                endif
            endif
            if( l_hasmskvols(state) )call cline_refine3D%set('mskfile', trim(params%mskvols(state)))
            ! run
            call xrefine3D_distr%execute(cline_refine3D)
            ! harvest outcome
            iter   = nint(cline_refine3D%get_rarg('endit'))
            oritab = trim(REFINE3D_ITER_FBODY)//int2str_pad(iter,3)//trim(METADATA_EXT)
            ! stash
            call binread_oritab(oritab, spproj_state, os_state, [1,params%nptcls])
            do iptcl = 1,params%nptcls
                if( nint(os_master%get(iptcl, 'state')) .ne. 1 )cycle
                call os_master%set_ori(iptcl, os_state%get_ori(iptcl))
                call os_master%set(iptcl, 'state', real(state))
            enddo
            call os_state%kill
            call refine3D_cleanup
        enddo

        ! final document
        call binwrite_oritab(FINAL_DOC, spproj_master, os_master, [1,params%nptcls])

        ! final reconstruction
        if( params%eo .eq.'no' )then
            call xreconstruct3D_distr%execute(cline_reconstruct3D_distr)
            do state = 1, params%nstates
                if( state_pops(state) == 0 )cycle
                if( l_singlestate .and. state.ne.params%state )cycle
                str_state = int2str_pad(state, 2)
                if( l_hasmskvols(state) )call cline_postprocess%set('mskfile', trim(params%mskvols(state)))
                vol = 'recvol_state'//trim(str_state)//params%ext
                call cline_postprocess%set('state', real(state))
                call cline_postprocess%set('fsc', trim(FSC_FBODY)//trim(str_state)//BIN_EXT)
                call cline_postprocess%set('vol_filt', trim(ANISOLP_FBODY)//trim(str_state)//params%ext)
                call xpostprocess%execute(cline_postprocess)
            enddo
        endif

        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER3D_REFINE NORMAL STOP ****')
        contains

            ! stash docs, volumes , etc.
            subroutine refine3D_cleanup
                character(len=STDLEN) :: src, dist
                character(len=*), parameter :: one = '01'
                character(len=3) :: str_iter
                integer          :: it, rename_stat
                dir = 'state_'//str_state//'/'
                call simple_mkdir(dir)
                do it = 1, iter
                    str_iter = int2str_pad(it,3)
                    ! volumes
                    src  = trim(VOL_FBODY)//one//'_iter'//str_iter//params%ext
                    dist = dir//trim(VOL_FBODY)//one//'_iter'//str_iter//params%ext
                    rename_stat = simple_rename(src, dist)
                    src  = trim(VOL_FBODY)//one//'_iter'//str_iter//trim(PPROC_SUFFIX)//params%ext
                    dist = dir//trim(VOL_FBODY)//one//'_iter'//str_iter//trim(PPROC_SUFFIX)//params%ext
                    rename_stat = simple_rename(src, dist)
                    ! e/o
                    if( params%eo.ne.'no')then
                        src  = trim(VOL_FBODY)//one//'_iter'//str_iter//'_even'//params%ext
                        dist = dir//trim(VOL_FBODY)//one//'_iter'//str_iter//'_even'//params%ext
                        rename_stat = simple_rename(src, dist)
                        src  = trim(VOL_FBODY)//one//'_iter'//str_iter//'_odd'//params%ext
                        dist = dir//trim(VOL_FBODY)//one//'_iter'//str_iter//'_odd'//params%ext
                        rename_stat = simple_rename(src, dist)
                        src = 'RESOLUTION_STATE'//one//'_ITER'//str_iter
                        rename_stat = simple_rename(src, dir//src)
                    endif
                    ! orientation document
                    src = trim(REFINE3D_ITER_FBODY)//str_iter//trim(METADATA_EXT)
                    if( file_exists(src) ) rename_stat = simple_rename(src, dir//src)
                enddo
                ! resolution measures
                if( params%eo.ne.'no')then
                    src  = trim(FSC_FBODY)//one//BIN_EXT
                    if( file_exists(src) ) rename_stat = simple_rename(src, dir//src)
                    src  = trim(FRCS_FBODY)//one//BIN_EXT
                    if( file_exists(src) ) rename_stat = simple_rename(src, dir//src)
                    src  = trim(ANISOLP_FBODY)//one//params%ext
                    if( file_exists(src) ) rename_stat = simple_rename(src, dir//src)
                endif
            end subroutine refine3D_cleanup

    end subroutine exec_cluster3D_refine

end module simple_commander_hlev_wflows
