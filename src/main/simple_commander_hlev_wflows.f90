! concrete commander: high-level workflows
module simple_commander_hlev_wflows
include 'simple_lib.f08'
use simple_binoris_io
use simple_cmdline,               only: cmdline
use simple_params,                only: params
use simple_commander_base,        only: commander_base
use simple_qsys_env,              only: qsys_env
use simple_oris,                  only: oris
use simple_sp_project,            only: sp_project
use simple_commander_distr_wflows ! use all in there
use simple_commander_distr        ! use all in there
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
        use simple_commander_cluster2D, only: rank_cavgs_commander
        use simple_commander_imgproc,   only: scale_commander
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
        type(sp_project)      :: spproj, spproj_sc
        type(params)          :: p_master
        character(len=:), allocatable :: projfile_sc, stk
        character(len=STDLEN) :: finalcavgs, finalcavgs_ranked, refs_sc
        real                  :: scale_stage1, scale_stage2
        integer               :: istk, nparts, last_iter_stage1, last_iter_stage2
        logical               :: scaling
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        ! parameters
        p_master = params(cline, del_scaled=.true.)
        ! delete mkdir flag from cline if present (to avoid nested directory structure)
        call cline%delete('mkdir')
        nparts   = p_master%nparts
        if( p_master%l_autoscale )then
            call cline%delete('objfun') ! stage dependent objective function
            ! SPLITTING
            call spproj%read(p_master%projfile )
            call spproj%split_stk(p_master%nparts)
            ! this workflow executes two stages of CLUSTER2D
            ! Stage 1: high down-scaling for fast execution, hybrid extremal/SHC optimisation for
            !          improved population distribution of clusters, no incremental learning,
            !          objective function is standard cross-correlation (cc)
            cline_cluster2D_stage1 = cline
            call cline_cluster2D_stage1%set('objfun', 'cc')
            call cline_cluster2D_stage1%delete('automsk')
            if( p_master%l_frac_update )then
                call cline_cluster2D_stage1%delete('update_frac') ! no incremental learning in stage 1
                call cline_cluster2D_stage1%set('maxits', real(MAXITS_STAGE1_EXTR))
            else
                call cline_cluster2D_stage1%set('maxits', real(MAXITS_STAGE1))
            endif
            ! Scaling
            call spproj%scale_projfile(p_master%smpd_targets2D(1), projfile_sc,&
                &cline_cluster2D_stage1, cline_scale1)
            call spproj%kill
            scale_stage1 = cline_scale1%get_rarg('scale')
            scaling      = trim(projfile_sc) /= trim(p_master%projfile)
            if( scaling )then
                call xscale_distr%execute( cline_scale1 )
                ! scale references
                if( cline%defined('refs') )then
                    call cline_scalerefs%set('stk', cline%get_carg('refs'))
                    refs_sc = trim(cline_scalerefs%get_carg('refs')) //trim(SCALE_SUFFIX)//p_master%ext
                    call cline_scalerefs%set('outstk', trim(refs_sc))
                    call cline_scalerefs%set('smpd', cline%get_rarg('smpd'))
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
                call spproj%read( p_master%projfile )
                spproj%os_ptcl2D = spproj_sc%os_ptcl2D
                call spproj%write()
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
            call cline_cluster2D_stage2%set('objfun', 'ccres')
            if( p_master%automsk .eq. 'yes' )call cline_cluster2D_stage2%set('automsk', 'cavg')
            call cline_cluster2D_stage2%set('startit', real(last_iter_stage1 + 1))
            if( cline%defined('update_frac') )then
                call cline_cluster2D_stage2%set('update_frac', p_master%update_frac)
            endif
            ! Scaling
            call spproj%read( p_master%projfile )
            call spproj%scale_projfile( p_master%smpd_targets2D(2), projfile_sc,&
                &cline_cluster2D_stage2, cline_scale2)
            call spproj%kill
            scale_stage2 = cline_scale2%get_rarg('scale')
            scaling      = trim(projfile_sc) /= p_master%projfile
            if( scaling ) call xscale_distr%execute( cline_scale2 )
            ! execution
            call cline_cluster2D_stage2%set('projfile', trim(projfile_sc))
            call xcluster2D_distr%execute(cline_cluster2D_stage2)
            last_iter_stage2 = nint(cline_cluster2D_stage2%get_rarg('endit'))
            finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//p_master%ext
            ! Updates project and references
            if( scaling )then
                ! shift modulation
                call spproj_sc%read_segment( 'ptcl2D', projfile_sc )
                call spproj_sc%os_ptcl2D%mul_shifts( 1./scale_stage2 )
                call spproj%read( p_master%projfile )
                spproj%os_ptcl2D = spproj_sc%os_ptcl2D
                call spproj%write()
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
                cline_make_cavgs = cline
                call cline_make_cavgs%delete('autoscale')
                call cline_make_cavgs%delete('balance')
                call cline_make_cavgs%set('prg',      'make_cavgs')
                call cline_make_cavgs%set('projfile', p_master%projfile)
                call cline_make_cavgs%set('nparts',   real(nparts))
                call cline_make_cavgs%set('refs',     trim(finalcavgs))
                call xmake_cavgs%execute(cline_make_cavgs)
            endif
        else
            ! no auto-scaling
            call xcluster2D_distr%execute(cline)
            last_iter_stage2 = nint(cline%get_rarg('endit'))
            finalcavgs       = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//p_master%ext
        endif
        ! adding cavgs to project
        call spproj%read( p_master%projfile )
        call spproj%add_cavgs2os_out( trim(finalcavgs), spproj%get_smpd())
        call spproj%write()
        call spproj%kill()
        ! ranking
        finalcavgs_ranked = trim(CAVGS_ITER_FBODY)//int2str_pad(last_iter_stage2,3)//'_ranked'//p_master%ext
        call cline_rank_cavgs%set('projfile', trim(p_master%projfile))
        call cline_rank_cavgs%set('stk',      trim(finalcavgs))
        call cline_rank_cavgs%set('outstk',   trim(finalcavgs_ranked))
        call cline_rank_cavgs%printline
        call xrank_cavgs%execute( cline_rank_cavgs )
        ! cleanup
        call del_file('start2Drefs'//p_master%ext)
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
    end subroutine exec_cluster2D_autoscale

    !> for generation of an initial 3d model from class averages
    subroutine exec_initial_3Dmodel( self, cline )
        use simple_commander_volops,  only: project_commander
        use simple_commander_rec,     only: reconstruct3D_commander
        class(initial_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! constants
        real,             parameter :: CENLP=30.           !< consistency with prime3D
        integer,          parameter :: MAXITS_SNHC=30, MAXITS_INIT=15, MAXITS_REFINE=40
        integer,          parameter :: STATE=1, NPROJS_SYMSRCH=50
        integer,          parameter :: NSPACE_SNHC = 1000, NSPACE_DEFAULT= 2500
        !character(len=*), parameter :: STKSCALEDBODY = 'stk_sc_initial_3Dmodel'
        ! distributed commanders
        type(prime3D_distr_commander)       :: xprime3D_distr
        type(symsrch_distr_commander)       :: xsymsrch_distr
        type(scale_project_distr_commander) :: xscale_distr
        ! shared-mem commanders
        type(reconstruct3D_commander) :: xreconstruct3D
        type(project_commander)       :: xproject
        ! command lines
        type(cmdline)         :: cline_refine3D_snhc
        type(cmdline)         :: cline_refine3D_init
        type(cmdline)         :: cline_refine3D_refine
        type(cmdline)         :: cline_symsrch
        type(cmdline)         :: cline_reconstruct3D
        type(cmdline)         :: cline_project
        type(cmdline)         :: cline_scale
        ! other variables
        type(ctfparams)       :: ctfvars ! ctf=no by default
        type(params)          :: p_master
        type(sp_project)      :: spproj, work_proj
        type(oris)            :: os
        real                  :: iter, smpd_target, lplims(2), orig_msk, msk, scale_factor
        integer               :: ncavgs, orig_box, box, istk
        character(len=2)      :: str_state
        character(len=:), allocatable :: projfile, stk, imgkind, WORK_PROJFILE
        character(len=STDLEN) :: vol_iter
        integer               :: status
        logical               :: srch4symaxis, doautoscale
        ! set cline defaults
        call cline%set('eo', 'no')
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! auto-scaling prep
        doautoscale = (cline%get_carg('autoscale').eq.'yes')
        ! now, remove autoscale flag from command line, since no scaled partial stacks
        ! will be produced (this program used shared-mem paralllelisation of scale)
        call cline%delete('autoscale')
        ! make master parameters
        call cline%set('oritype', 'out')
        p_master = params(cline)
        allocate(WORK_PROJFILE, source='initial_3Dmodel_tmproj.simple')
        call del_file(WORK_PROJFILE)
        call cline%set('oritype', 'ptcl3D')
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
        ! First wee need a dummy project for refine3D
        ! fetch cavgs
        call spproj%read(p_master%projfile)
        ncavgs       = nint(spproj%os_out%get(1,'nptcls'))
        ctfvars%smpd = p_master%smpd
        if( spproj%os_out%isthere('imgkind') )then
            call spproj%os_out%getter(1, 'imgkind', imgkind)
            if( trim(imgkind).ne.'cavg' )then
                stop 'IMGKIND should be CAVG; simple_commander_hlev_wflows :: initial_3Dmodel'
            endif
        else
            stop 'IMGKIND should be informed; simple_commander_hlev_wflows :: initial_3Dmodel'
        endif
        if( spproj%os_out%isthere('stk') )then
            call spproj%os_out%getter(1, 'stk', stk)
        else
            stop 'No class-average stack found; simple_commander_hlev_wflows :: initial_3Dmodel'
        endif
        ! fetch general project info
        work_proj%projinfo = spproj%projinfo
        work_proj%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris()>0 ) work_proj%jobproc = spproj%jobproc
        call os%new_clean( ncavgs )
        call work_proj%add_single_stk(trim(stk), ctfvars, os)
        call spproj%kill
        call os%kill
        ! name change
        call work_proj%projinfo%delete_entry('projname')
        call work_proj%projinfo%delete_entry('projfile')
        call cline%set('projfile', trim(WORK_PROJFILE))
        call cline%set('projname', trim(get_fbody(trim(WORK_PROJFILE),trim('simple'))))
        call work_proj%update_projinfo(cline)
        ! split
        call work_proj%split_stk(p_master%nparts)
        ! Scaling
        orig_box    = work_proj%get_box()
        orig_msk    = p_master%msk
        smpd_target = max(p_master%smpd, lplims(2)*LP2SMPDFAC)
        doautoscale = doautoscale .and. smpd_target > work_proj%get_smpd()
        projfile    = WORK_PROJFILE
        if( doautoscale )then
            deallocate(WORK_PROJFILE)
            call work_proj%scale_projfile(smpd_target, WORK_PROJFILE, cline, cline_scale)
            call cline_scale%set('nparts',1.) ! shared memory
            scale_factor = cline_scale%get_rarg('scale')
            box          = nint(cline_scale%get_rarg('newbox'))
            msk          = cline%get_rarg('msk')
            call xscale_distr%execute( cline_scale )
        else
            box         = orig_box
            msk         = orig_msk
            smpd_target = p_master%smpd
            projfile    = p_master%projfile
        endif
        ! prepare command lines from prototype for original scaling
        cline_reconstruct3D   = cline
        cline_project         = cline
        ! updates scaling parameters
        call cline%set('projfile', trim(WORK_PROJFILE))
        call cline%set('msk',      msk)
        call cline%set('box',      real(box))
        ! prepare command lines from prototype
        cline_refine3D_snhc   = cline
        cline_refine3D_init   = cline
        cline_refine3D_refine = cline
        cline_symsrch         = cline
        ! reconstruct3D & project are not distributed executions, so remove the nparts flag
        call cline_reconstruct3D%delete('nparts')
        call cline_project%delete('nparts')
        ! projects names are subject to change depending on scaling
        call cline_reconstruct3D%delete('projname')
        call cline_reconstruct3D%delete('projfile')
        call cline_project%delete('projname')
        call cline_project%delete('projfile')
        ! initialise command line parameters
        ! (1) INITIALIZATION BY STOCHASTIC NEIGHBORHOOD HILL-CLIMBING
        call cline_refine3D_snhc%delete('update_frac') ! no fractional update in first phase
        call cline_refine3D_snhc%set('prg',    'refine3D')
        call cline_refine3D_snhc%set('maxits', real(MAXITS_SNHC))
        call cline_refine3D_snhc%set('refine', 'snhc')
        call cline_refine3D_snhc%set('dynlp',  'no') ! better be explicit about the dynlp
        call cline_refine3D_snhc%set('lp',     lplims(1))
        call cline_refine3D_snhc%set('nspace', real(NSPACE_SNHC))
        call cline_refine3D_snhc%set('objfun', 'cc')
        ! (2) refine3D_init
        call cline_refine3D_init%set('prg',    'refine3D')
        call cline_refine3D_init%set('maxits', real(MAXITS_INIT))
        call cline_refine3D_init%set('vol1',   trim(SNHCVOL)//trim(str_state)//p_master%ext)
        call cline_refine3D_init%set('dynlp',  'no') ! better be explicit about the dynlp
        call cline_refine3D_init%set('lp',     lplims(1))
        if( .not. cline_refine3D_init%defined('nspace') )then
            call cline_refine3D_init%set('nspace', real(NSPACE_DEFAULT))
        endif
        if( .not. cline_refine3D_init%defined('objfun') )then
            call cline_refine3D_init%set('objfun', 'ccres')
        endif
        ! (3) SYMMETRY AXIS SEARCH
        if( srch4symaxis )then
            ! need to replace original point-group flag with c1
            call cline_refine3D_snhc%set('pgrp', 'c1')
            call cline_refine3D_init%set('pgrp', 'c1')
            ! symsrch
            call cline_symsrch%set('prg', 'symsrch')
            call cline_symsrch%delete('stk')  ! volumetric symsrch
            call cline_symsrch%set('nptcls',  real(NPROJS_SYMSRCH))
            call cline_symsrch%set('nspace',  real(NPROJS_SYMSRCH))
            call cline_symsrch%set('cenlp',   CENLP)
            call cline_symsrch%set('outfile', 'symdoc'//trim(METADATA_EXT))
            call cline_symsrch%set('lp',      lplims(2))
            ! (4.5) RECONSTRUCT SYMMETRISED VOLUME
            call cline_reconstruct3D%set('prg',      'reconstruct3D')
            call cline_reconstruct3D%set('projfile', trim(WORK_PROJFILE))
            call cline_reconstruct3D%set('trs',      5.) ! to assure that shifts are being used
            call cline_reconstruct3D%set('oritab',   'symdoc'//trim(METADATA_EXT))
            ! refinement step now uses the symmetrised vol and doc
            ! call cline_refine3D_refine%set('oritab', 'symdoc'//trim(METADATA_EXT)) !! TO UPDATE
            call cline_refine3D_refine%set('vol1',   'rec_sym'//p_master%ext)
        endif
        ! (4) PRIME3D REFINE STEP
        call cline_refine3D_refine%set('prg', 'refine3D')
        call cline_refine3D_refine%set('maxits', real(MAXITS_REFINE))
        call cline_refine3D_refine%set('refine', 'single')
        call cline_refine3D_refine%set('dynlp', 'no') ! better be explicit about the dynlp
        call cline_refine3D_refine%set('lp', lplims(2))
        if( .not. cline_refine3D_refine%defined('nspace') )then
            call cline_refine3D_refine%set('nspace', real(NSPACE_DEFAULT))
        endif
        if( .not. cline_refine3D_refine%defined('objfun') )then
            call cline_refine3D_refine%set('objfun', 'ccres')
        endif
        ! (5) RE-PROJECT VOLUME
        call cline_project%set('prg',    'project')
        call cline_project%set('outstk', 'reprojs'//p_master%ext)
        call cline_project%set('smpd',    p_master%smpd)
        call cline_project%set('msk',     orig_msk)
        call cline_project%set('box',     real(orig_box))
        ! execute commanders
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> INITIALIZATION WITH STOCHASTIC NEIGHBORHOOD HILL-CLIMBING'
        write(*,'(A)') '>>>'
        call xprime3D_distr%execute(cline_refine3D_snhc)
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> INITIAL 3D MODEL GENERATION WITH PRIME3D'
        write(*,'(A)') '>>>'
        call xprime3D_distr%execute(cline_refine3D_init)
        iter = cline_refine3D_init%get_rarg('endit')
        call set_iter_dependencies
        if( srch4symaxis )then
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> SYMMETRY AXIS SEARCH'
            write(*,'(A)') '>>>'
            !call cline_symsrch%set('oritab', trim(oritab)) ! TO UPDATE
            call cline_symsrch%set('vol1', trim(vol_iter))
            call xsymsrch_distr%execute(cline_symsrch)
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> 3D RECONSTRUCTION OF SYMMETRISED VOLUME'
            write(*,'(A)') '>>>'
            call xreconstruct3D%execute(cline_reconstruct3D)
            status = simple_rename(trim(VOL_FBODY)//trim(str_state)//p_master%ext, 'rec_sym'//p_master%ext)
        else
            ! refinement step needs to use iter dependent vol/oritab
            call cline_refine3D_refine%set('vol1', trim(vol_iter))
        endif
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> PRIME3D REFINEMENT STEP'
        write(*,'(A)') '>>>'
        call cline_refine3D_refine%set('startit', iter + 1.0)
        call xprime3D_distr%execute(cline_refine3D_refine)
        iter = cline_refine3D_refine%get_rarg('endit')
        call set_iter_dependencies
        call work_proj%read_segment( 'ptcl3D', WORK_PROJFILE )
        os = work_proj%os_ptcl3D
        if( doautoscale )then
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> 3D RECONSTRUCTION AT ORIGINAL SAMPLING'
            write(*,'(A)') '>>>'
            ! modulate shifts
            call os%mul_shifts( 1./scale_factor )
            ! clean stacks
            call work_proj%read_segment( 'stk', WORK_PROJFILE )
            do istk=1,work_proj%os_stk%get_noris()
                call work_proj%os_stk%getter(istk, 'stk', stk)
                call del_file(trim(stk))
            enddo
            call work_proj%kill()
            call del_file(trim(WORK_PROJFILE))
            deallocate(WORK_PROJFILE)
            ! updates original scale project
            WORK_PROJFILE = trim(projfile)
            call work_proj%read(WORK_PROJFILE)
            work_proj%os_ptcl3D = os
            call work_proj%write()
            ! re-reconstruct volume
            call cline_reconstruct3D%set('projfile',WORK_PROJFILE)
            call xreconstruct3D%execute(cline_reconstruct3D)
            status = simple_rename(trim(VOL_FBODY)//trim(str_state)//p_master%ext, 'rec_final'//p_master%ext)
        else
            status = simple_rename(trim(vol_iter), 'rec_final'//p_master%ext)
        endif
        call work_proj%kill()
        ! update the original project cls3D segment with orientations
        call spproj%read(p_master%projfile)
        call os%delete_entry('stkind')
        spproj%os_cls3D = os
        call spproj%write()
        ! reprojections
        call os%write('final_oris.txt')
        write(*,'(A)') '>>>'
        write(*,'(A)') '>>> RE-PROJECTION OF THE FINAL VOLUME'
        write(*,'(A)') '>>>'
        call cline_project%set('projfile', trim(WORK_PROJFILE))
        call cline_project%set('vol1',     'rec_final'//p_master%ext)
        call cline_project%set('oritab',   'final_oris.txt')
        call xproject%execute(cline_project)
        ! end gracefully
        call del_file(WORK_PROJFILE)
        call simple_end('**** SIMPLE_INITIAL_3DMODEL NORMAL STOP ****')

        contains

            subroutine set_iter_dependencies
                character(len=3) :: str_iter
                str_iter = int2str_pad(nint(iter),3)
                !oritab   = trim(REFINE3D_ITER_FBODY)//trim(str_iter)//trim(METADATA_EXT)
                vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//p_master%ext
            end subroutine set_iter_dependencies

    end subroutine exec_initial_3Dmodel

    !> for heterogeinity analysis
    subroutine exec_cluster3D( self, cline )
        use simple_defs_conv
        use simple_commander_rec,    only: reconstruct3D_commander
        use simple_commander_volops, only: postprocess_commander
        use simple_sym,              only: sym
        class(cluster3D_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        ! constants
        integer,          parameter :: MAXITS1 = 50
        integer,          parameter :: MAXITS2 = 40
        character(len=*), parameter :: one     = '01'
        ! distributed commanders
        type(prime3D_distr_commander)        :: xprime3D_distr
        type(reconstruct3D_distr_commander)  :: xreconstruct3D_distr
        ! shared-mem commanders
        type(postprocess_commander) :: xpostprocess
        ! command lines
        type(cmdline)               :: cline_refine3D1, cline_refine3D2
        type(cmdline)               :: cline_reconstruct3D_distr, cline_reconstruct3D_mixed_distr
        type(cmdline)               :: cline_postprocess
        ! other variables
        type(sym)                     :: symop
        integer,     allocatable      :: labels(:)
        logical,     allocatable      :: included(:)
        type(params)                  :: p_master
        type(sp_project)              :: spproj, spproj_states
        class(oris), pointer          :: os => null(), os_states => null()
        character(len=STDLEN)         :: oritab, str_state
        character(len=:), allocatable :: recname, rhoname, part_str
        real                          :: trs
        integer                       :: state, iter, n_incl, startit
        integer                       :: rename_stat, ipart
        call seed_rnd
        ! sanity check
        if(nint(cline%get_rarg('nstates')) <= 1)&
            &stop 'Non-sensical NSTATES argument for heterogeneity analysis!'

        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')

        ! make master parameters
        p_master = params(cline)
        if( p_master%eo .eq. 'no' .and. .not. cline%defined('lp') )&
            &stop 'need lp input when eo .eq. no; cluster3D'

        ! prepare command lines from prototype
        call cline%delete('refine')
        cline_refine3D1                 = cline
        cline_refine3D2                 = cline
        cline_postprocess               = cline ! eo always eq yes
        cline_reconstruct3D_distr       = cline ! eo always eq yes
        cline_reconstruct3D_mixed_distr = cline ! eo always eq yes
        call cline_refine3D1%set('prg', 'refine3D')
        call cline_refine3D2%set('prg', 'refine3D')
        call cline_refine3D1%set('maxits', real(MAXITS1))
        select case(trim(p_master%refine))
            case('sym')
                call cline_refine3D1%set('refine', 'clustersym')
            case DEFAULT
                if( .not.cline%defined('refine') )then
                    call cline_refine3D1%set('refine', 'cluster')
                else
                    call cline_refine3D1%set('refine', trim(p_master%refine))
                endif
        end select
        call cline_refine3D2%set('refine', 'multi')
        call cline_refine3D1%set('dynlp', 'no')
        call cline_refine3D1%delete('oritab2')
        call cline_refine3D1%delete('update_frac')
        call cline_refine3D1%set('frcs', 'cluster3D_frcs.bin') !! mixed model FRCs
        call cline_refine3D2%set('dynlp', 'no')
        call cline_refine3D2%delete('oritab2')
        if( .not.cline%defined('update_frac') )call cline_refine3D2%set('update_frac', 0.5)
        call cline_reconstruct3D_distr%set('prg', 'reconstruct3D')
        call cline_reconstruct3D_mixed_distr%set('prg', 'reconstruct3D')
        call cline_postprocess%set('prg', 'postprocess')
        call cline_reconstruct3D_distr%delete('lp')
        call cline_reconstruct3D_mixed_distr%delete('lp')
        call cline_postprocess%delete('lp')
        call cline_reconstruct3D_distr%set('eo','yes')
        call cline_reconstruct3D_mixed_distr%set('nstates', 1.)
        call cline_postprocess%set('eo','yes')
        if( cline%defined('trs') )then
            ! all good
        else
            ! works out shift lmits for in-plane search
            trs = MSK_FRAC*real(p_master%msk)
            trs = min(MAXSHIFT, max(MINSHIFT, trs))
            call cline_refine3D1%set('trs',trs)
            call cline_refine3D2%set('trs',trs)
        endif

        ! generate diverse initial labels & orientations
        oritab = 'cluster3Ddoc_init'//trim(METADATA_EXT)
        call spproj%new_seg_with_ptr(p_master%nptcls, p_master%oritype, os)
        call binread_oritab(p_master%oritab, spproj, os, [1,p_master%nptcls])
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
            call binwrite_oritab(trim('symrnd_'//oritab), spproj, os, [1,p_master%nptcls])
        endif
        if( cline%defined('oritab2') )then
            ! this is  to force initialisation (4 testing)
            call spproj_states%new_seg_with_ptr(p_master%nptcls, p_master%oritype, os_states)
            call binread_oritab(p_master%oritab2, spproj_states, os_states, [1,p_master%nptcls])
            labels = nint(os_states%get_all('state'))
            call os%set_all('state', real(labels))
            call os_states%kill
        else if( .not. cline%defined('startit') )then
            write(*,'(A)') '>>>'
            write(*,'(A)') '>>> GENERATING DIVERSE LABELING'
            call diverse_labeling(os, p_master%nstates, labels, corr_ranked=.true.)
            call os%set_all('state', real(labels))
        else
            ! starting from a previous solution
            labels = nint(os%get_all('state'))
        endif

        ! to accomodate state=0s in oritab input
        included = os%included()
        n_incl   = count(included)
        where( .not. included ) labels = 0
        call os%set_all('state', real(labels))
        call binwrite_oritab(trim(oritab), spproj, os, [1,p_master%nptcls])
        call cline_refine3D1%set('oritab', trim(oritab))
        deallocate(labels, included)

        ! retrieve mixed model Fourier components, normalization matrix, FSC & anisotropic filter
        if( trim(p_master%refine) .eq. 'sym' )then
            call cline%set('oritab', trim('symrnd_'//oritab))
            call cline%set('pgrp', 'c1')
        endif
        call xreconstruct3D_distr%execute( cline_reconstruct3D_mixed_distr )
        rename_stat = simple_rename(trim(VOL_FBODY)//one//p_master%ext, trim(CLUSTER3D_VOL)//p_master%ext)
        if( p_master%eo .ne. 'no' )then
            rename_stat = simple_rename(trim(VOL_FBODY)//one//'_even'//p_master%ext, trim(CLUSTER3D_VOL)//'_even'//p_master%ext)
            rename_stat = simple_rename(trim(VOL_FBODY)//one//'_odd'//p_master%ext,  trim(CLUSTER3D_VOL)//'_odd'//p_master%ext)
            rename_stat = simple_rename(trim(FSC_FBODY)//one//BIN_EXT, trim(CLUSTER3D_FSC))
            rename_stat = simple_rename(trim(FRCS_FBODY)//one//BIN_EXT, trim(CLUSTER3D_FRCS))
            rename_stat = simple_rename(trim(ANISOLP_FBODY)//one//p_master%ext, trim(CLUSTER3D_ANISOLP)//p_master%ext)
        endif
        ! THis is experimental to keep FCs & RHO matrices from mixed model
        ! do ipart = 1, p_master%nparts
        !     allocate(part_str, source=int2str_pad(ipart,p_master%numlen))
        !     allocate(recname, source=trim(VOL_FBODY)//one//'_part'//part_str//p_master%ext)
        !     allocate(rhoname, source='rho_'//trim(VOL_FBODY)//one//'_part'//part_str//p_master%ext)
        !     rename_stat = rename(recname, 'cluster3D_'//trim(recname))
        !     rename_stat = rename(rhoname, 'cluster3D_'//trim(rhoname))
        !     deallocate(part_str,recname,rhoname)
        ! enddo

        ! STAGE1: frozen orientation parameters
        write(*,'(A)')    '>>>'
        write(*,'(A,I3)') '>>> PRIME3D - STAGE 1'
        write(*,'(A)')    '>>>'
        call xprime3D_distr%execute(cline_refine3D1)
        iter   = nint(cline_refine3D1%get_rarg('endit'))
        oritab = trim(REFINE3D_ITER_FBODY)//int2str_pad(iter,3)//trim(METADATA_EXT)
        call binread_oritab(trim(oritab), spproj, os, [1,p_master%nptcls])
        oritab = 'cluster3Ddoc_stage1'//trim(METADATA_EXT)
        call binwrite_oritab(trim(oritab), spproj, os, [1,p_master%nptcls])

        ! ! stage 1 reconstruction to obtain resolution estimate when eo .eq. 'no'
        ! if( p_master%eo .eq. 'no' )then
        !     call cline_reconstruct3D_distr%set('oritab', trim(oritab))
        !     call xreconstruct3D_distr%execute(cline_reconstruct3D_distr)
        !     do state = 1, p_master%nstates
        !         str_state  = int2str_pad(state, 2)
        !         call cline_postprocess%set('vol1', trim(VOL_FBODY)//trim(str_state)//p_master%ext)
        !         call cline_postprocess%set('fsc', trim(FSC_FBODY)//trim(str_state)//BIN_EXT)
        !         call cline_postprocess%set('vol_filt', trim(ANISOLP_FBODY)//trim(str_state)//p_master%ext)
        !         call xpostprocess%execute(cline_postprocess)
        !     enddo
        ! endif

        ! STAGE2: soft multi-states refinement
        startit = iter + 1
        call cline_refine3D2%set('startit', real(startit))
        call cline_refine3D2%set('maxits',  real(startit + MAXITS2))
        call cline_refine3D2%set('oritab',  trim(oritab))
        write(*,'(A)')    '>>>'
        write(*,'(A,I3)') '>>> PRIME3D - STAGE 2'
        write(*,'(A)')    '>>>'
        call xprime3D_distr%execute(cline_refine3D2)
        iter   = nint(cline_refine3D2%get_rarg('endit'))
        oritab = trim(REFINE3D_ITER_FBODY)//int2str_pad(iter,3)//trim(METADATA_EXT)
        call binread_oritab(trim(oritab), spproj, os, [1,p_master%nptcls])
        oritab = 'cluster3Ddoc_stage2'//trim(METADATA_EXT)
        call binwrite_oritab(trim(oritab), spproj, os, [1,p_master%nptcls])

        ! stage 2 reconstruction to obtain resolution estimate when eo .eq. 'no'
        if( p_master%eo .eq. 'no' )then
            call cline_reconstruct3D_distr%set('oritab', trim(oritab))
            call xreconstruct3D_distr%execute(cline_reconstruct3D_distr)
            do state = 1, p_master%nstates
                str_state  = int2str_pad(state, 2)
                call cline_postprocess%set('vol1', trim(VOL_FBODY)//trim(str_state)//p_master%ext)
                call cline_postprocess%set('fsc', trim(FSC_FBODY)//trim(str_state)//BIN_EXT)
                call cline_postprocess%set('vol_filt', trim(ANISOLP_FBODY)//trim(str_state)//p_master%ext)
                call xpostprocess%execute(cline_postprocess)
            enddo
        endif

        ! cleanup
        call os%kill

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
        use simple_defs_conv
        use simple_commander_rec,    only: reconstruct3D_commander
        use simple_commander_volops, only: postprocess_commander
        class(cluster3D_refine_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! constants
        integer,            parameter :: MAXITS = 40
        character(len=:), allocatable :: INIT_FBODY
        character(len=:), allocatable :: FINAL_FBODY
        character(len=:), allocatable :: FINAL_DOC
        ! distributed commanders
        type(prime3D_distr_commander)       :: xprime3D_distr
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
        type(params)             :: p_master
        type(sp_project)         :: spproj_master
        type(sp_project), target :: spproj_state
        class(oris), pointer     :: os_master => null(), os_state => null()
        character(len=STDLEN)    :: oritab, vol, fname
        character(len=9)         :: dir
        character(len=2)         :: str_state
        integer                  :: state, iptcl, iter
        logical                  :: l_singlestate, error
        integer                  :: rename_stat

        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')

        ! make master parameters
        p_master      = params(cline)
        l_singlestate = cline%defined('state')
        error         = .false.

        ! filenaming strings allocation
        allocate(INIT_FBODY, source='cluster3Dinit_refine_state')
        allocate(FINAL_FBODY, source='cluster3Ddoc_refine_state')
        allocate(FINAL_DOC, source='cluster3Ddoc_refine'//trim(METADATA_EXT))

        ! sanity checks
        if( p_master%eo .eq. 'no' .and. .not. cline%defined('lp') )&
            &stop 'need lp input when eo .eq. no; cluster3D_refine'
        if(.not.file_exists(p_master%oritab))then
            print *,'Document ',trim(p_master%oritab),' does not exist!'
            stop
        endif

        ! general prep
        p_master%nptcls = binread_nlines(p_master, p_master%oritab)
        call spproj_master%new_seg_with_ptr(p_master%nptcls, p_master%oritype, os_master)
        call binread_oritab(p_master%oritab, spproj_master, os_master, [1,p_master%nptcls])
        p_master%nstates = os_master%get_n('state')
        if(p_master%nstates < 2 .and. .not.l_singlestate)then
            print *, 'Non-sensical number of states argument for heterogeneity refinement: ',p_master%nstates
            stop
        endif
        allocate(l_hasmskvols(p_master%nstates), source=.false.)
        allocate(l_hasvols(p_master%nstates),    source=.false.)
        allocate(init_docs(p_master%nstates), final_docs(p_master%nstates))
        call os_master%get_pops(state_pops, 'state', consider_w=.false.)
        do state = 1, p_master%nstates
            if( state_pops(state) == 0 )cycle
            if( l_singlestate .and. state.ne.p_master%state )cycle
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
            call binwrite_oritab(init_docs(state), spproj_state, os_state, [1,p_master%nptcls])
            final_docs(state) = trim(FINAL_FBODY)//str_state//trim(METADATA_EXT)
            ! check & move volumes
            l_hasvols(state) = trim(p_master%vols(state)) .ne. ''
            if( l_hasvols(state) )then
                if( p_master%eo .ne. 'no' )then
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
                    fname = trim(ANISOLP_FBODY)//str_state//p_master%ext
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                    else
                        rename_stat = simple_rename(fname, dir//trim(fname))
                    endif
                    ! e/o
                    fname = add2fbody(trim(p_master%vols(state)), p_master%ext, '_even')
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                    else
                        rename_stat = simple_rename(fname, dir//trim(fname))
                    endif
                    fname = add2fbody(trim(p_master%vols(state)), p_master%ext, '_odd')
                    if( .not.file_exists(fname) )then
                        print *, 'File missing: ', trim(fname)
                    else
                        rename_stat = simple_rename(fname, dir//trim(fname))
                    endif
                endif
                ! volume
                if( .not.file_exists(p_master%vols(state)) )then
                    print *, 'File missing: ', p_master%vols(state)
                    error = .true.
                else
                    fname = trim(p_master%vols(state))
                    p_master%vols(state) = dir//trim(fname) ! name change
                    rename_stat = simple_rename(fname, p_master%vols(state))
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
                    call simple_copy_file(trim(fname), trim(p_master%mskvols(state)))
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
        cline_refine3D_master     = cline
        cline_reconstruct3D_distr = cline
        cline_postprocess         = cline
        call cline_refine3D_master%set('prg', 'refine3D')
        call cline_refine3D_master%set('dynlp', 'no')
        if( .not.cline%defined('maxits') ) call cline_refine3D_master%set('maxits', real(MAXITS))
        call cline_reconstruct3D_distr%set('prg', 'reconstruct3D')
        call cline_reconstruct3D_distr%set('oritab', trim(FINAL_DOC))
        call cline_reconstruct3D_distr%set('nstates', trim(int2str(p_master%nstates)))
        call cline_reconstruct3D_distr%set('eo', 'yes')
        call cline_postprocess%delete('lp')

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
            cline_refine3D = cline_refine3D_master
            call cline_refine3D%set('oritab', trim(init_docs(state)))
            if( l_hasvols(state) )then
                call cline_refine3D%set('vol1', trim(p_master%vols(state)))
                if( p_master%eo.ne.'no' )then
                    fname = dir//trim(FSC_FBODY)//str_state//BIN_EXT
                    call simple_copy_file(trim(fname),'fsc_state01.bin')
                    fname = dir//trim(FRCS_FBODY)//str_state//BIN_EXT
                    if(file_exists(fname))call simple_copy_file(trim(fname),'frcs_state01.bin')
                    fname = dir//trim(ANISOLP_FBODY)//str_state//p_master%ext
                    if(file_exists(fname))call simple_copy_file(trim(fname),'aniso_optlp_state01.mrc')
                endif
            endif
            if( l_hasmskvols(state) )call cline_refine3D%set('mskfile', trim(p_master%mskvols(state)))
            ! run prime3D
            call xprime3D_distr%execute(cline_refine3D)
            ! harvest outcome
            iter   = nint(cline_refine3D%get_rarg('endit'))
            oritab = trim(REFINE3D_ITER_FBODY)//int2str_pad(iter,3)//trim(METADATA_EXT)
            ! stash
            call binread_oritab(oritab, spproj_state, os_state, [1,p_master%nptcls])
            do iptcl = 1,p_master%nptcls
                if( nint(os_master%get(iptcl, 'state')) .ne. 1 )cycle
                call os_master%set_ori(iptcl, os_state%get_ori(iptcl))
                call os_master%set(iptcl, 'state', real(state))
            enddo
            call os_state%kill
            call prime3d_cleanup
        enddo

        ! final document
        call binwrite_oritab(FINAL_DOC, spproj_master, os_master, [1,p_master%nptcls])

        ! final reconstruction
        if( p_master%eo .eq.'no' )then
            call xreconstruct3D_distr%execute(cline_reconstruct3D_distr)
            do state = 1, p_master%nstates
                if( state_pops(state) == 0 )cycle
                if( l_singlestate .and. state.ne.p_master%state )cycle
                str_state = int2str_pad(state, 2)
                if( l_hasmskvols(state) )call cline_postprocess%set('mskfile', trim(p_master%mskvols(state)))
                vol = 'recvol_state'//trim(str_state)//p_master%ext
                call cline_postprocess%set('vol1', trim(vol))
                call cline_postprocess%set('fsc', trim(FSC_FBODY)//trim(str_state)//BIN_EXT)
                call cline_postprocess%set('vol_filt', trim(ANISOLP_FBODY)//trim(str_state)//p_master%ext)
                call xpostprocess%execute(cline_postprocess)
            enddo
        endif

        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER3D_REFINE NORMAL STOP ****')
        contains

            ! stash docs, volumes , etc.
            subroutine prime3d_cleanup
                character(len=STDLEN) :: src, dist
                character(len=*), parameter :: one = '01'
                character(len=3) :: str_iter
                integer          :: it, rename_stat
                dir = 'state_'//str_state//'/'
                call simple_mkdir(dir)
                do it = 1, iter
                    str_iter = int2str_pad(it,3)
                    ! volumes
                    src  = trim(VOL_FBODY)//one//'_iter'//str_iter//p_master%ext
                    dist = dir//trim(VOL_FBODY)//one//'_iter'//str_iter//p_master%ext
                    rename_stat = simple_rename(src, dist)
                    src  = trim(VOL_FBODY)//one//'_iter'//str_iter//'_pproc'//p_master%ext
                    dist = dir//trim(VOL_FBODY)//one//'_iter'//str_iter//'_pproc'//p_master%ext
                    rename_stat = simple_rename(src, dist)
                    ! e/o
                    if( p_master%eo.ne.'no')then
                        src  = trim(VOL_FBODY)//one//'_iter'//str_iter//'_even'//p_master%ext
                        dist = dir//trim(VOL_FBODY)//one//'_iter'//str_iter//'_even'//p_master%ext
                        rename_stat = simple_rename(src, dist)
                        src  = trim(VOL_FBODY)//one//'_iter'//str_iter//'_odd'//p_master%ext
                        dist = dir//trim(VOL_FBODY)//one//'_iter'//str_iter//'_odd'//p_master%ext
                        rename_stat = simple_rename(src, dist)
                        src = 'RESOLUTION_STATE'//one//'_ITER'//str_iter
                        rename_stat = simple_rename(src, dir//src)
                    endif
                    ! orientation document
                    src = trim(REFINE3D_ITER_FBODY)//str_iter//trim(METADATA_EXT)
                    if( file_exists(src) ) rename_stat = simple_rename(src, dir//src)
                enddo
                ! resolution measures
                if( p_master%eo.ne.'no')then
                    src  = trim(FSC_FBODY)//one//BIN_EXT
                    if( file_exists(src) ) rename_stat = simple_rename(src, dir//src)
                    src  = trim(FRCS_FBODY)//one//BIN_EXT
                    if( file_exists(src) ) rename_stat = simple_rename(src, dir//src)
                    src  = trim(ANISOLP_FBODY)//one//p_master%ext
                    if( file_exists(src) ) rename_stat = simple_rename(src, dir//src)
                endif
            end subroutine prime3d_cleanup

    end subroutine exec_cluster3D_refine

end module simple_commander_hlev_wflows
