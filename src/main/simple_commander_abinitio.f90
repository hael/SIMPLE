! concrete commander: high-level workflows
module simple_commander_abinitio
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters
use simple_sp_project,     only: sp_project
use simple_qsys_funs
implicit none

public :: initial_3Dmodel_commander_hlev
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: initial_3Dmodel_commander_hlev
  contains
    procedure :: execute      => exec_initial_3Dmodel
end type initial_3Dmodel_commander_hlev

contains

    !> for generation of an initial 3d model from class averages
    subroutine exec_initial_3Dmodel( self, cline )
        use simple_commander_rec,      only: reconstruct3D_commander_distr
        use simple_commander_project,  only: scale_project_commander_distr
        use simple_commander_refine3D, only: refine3D_commander_distr
        use simple_oris,               only: oris
        use simple_ori,                only: ori
        use simple_image,              only: image
        use simple_commander_volops,   only: reproject_commander, symaxis_search_commander, postprocess_commander
        use simple_parameters,         only: params_glob
        use simple_qsys_env,           only: qsys_env
        use simple_sym,                only: sym
        class(initial_3Dmodel_commander_hlev), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        ! constants
        real,                  parameter :: SCALEFAC2_TARGET = 0.5
        real,                  parameter :: CENLP=30. !< consistency with refine3D
        integer,               parameter :: MAXITS_SNHC=30, MAXITS_INIT=15, MAXITS_REFINE=40
        integer,               parameter :: NSPACE_SNHC=1000, NSPACE_INIT=1000, NSPACE_REFINE=2500
        character(len=STDLEN), parameter :: ORIG_WORK_PROJFILE   = 'initial_3Dmodel_tmpproj.simple'
        character(len=STDLEN), parameter :: REC_FBODY            = 'rec_final'
        character(len=STDLEN), parameter :: REC_PPROC_FBODY      = trim(REC_FBODY)//trim(PPROC_SUFFIX)
        character(len=STDLEN), parameter :: REC_PPROC_MIRR_FBODY = trim(REC_PPROC_FBODY)//trim(MIRR_SUFFIX)
        character(len=2) :: str_state
        ! distributed commanders
        type(refine3D_commander_distr)      :: xrefine3D_distr
        type(scale_project_commander_distr) :: xscale_distr
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        ! shared-mem commanders
        type(symaxis_search_commander) :: xsymsrch
        type(reproject_commander)      :: xreproject
        type(postprocess_commander)    :: xpostprocess
        ! command lines
        type(cmdline) :: cline_refine3D_snhc, cline_refine3D_init, cline_refine3D_refine
        type(cmdline) :: cline_symsrch
        type(cmdline) :: cline_reconstruct3D, cline_postprocess
        type(cmdline) :: cline_reproject
        type(cmdline) :: cline_scale1, cline_scale2
        ! other
        character(len=:), allocatable :: stk, orig_stk, frcs_fname
        character(len=:), allocatable :: WORK_PROJFILE
        real,             allocatable :: res(:), tmp_rarr(:)
        integer,          allocatable :: states(:), tmp_iarr(:)
        type(qsys_env)        :: qenv
        type(parameters)      :: params
        type(ctfparams)       :: ctfvars ! ctf=yes by default
        type(sp_project)      :: spproj, work_proj1, work_proj2
        type(oris)            :: os
        type(ori)             :: o_tmp
        type(sym)             :: se1,se2
        type(image)           :: img, vol
        character(len=STDLEN) :: vol_iter, pgrp_init, pgrp_refine
        real                  :: iter, smpd_target, lplims(2), msk, orig_msk, orig_smpd
        real                  :: scale_factor1, scale_factor2
        integer               :: icls, ncavgs, orig_box, box, istk, status, cnt
        logical               :: srch4symaxis, do_autoscale, symran_before_refine, l_lpset
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('ptclw') )     call cline%set('ptclw',      'no')
        ! hard set oritype
        call cline%set('oritype', 'out') ! because cavgs are part of out segment
        ! class averages, so no CTF!!
        ctfvars%ctfflag = CTFFLAG_NO
        ! auto-scaling prep
        do_autoscale = (cline%get_carg('autoscale').eq.'yes')
        if( cline%defined('vol1') ) do_autoscale = .false.
        ! now, remove autoscale flag from command line, since no scaled partial stacks
        ! will be produced (this program used shared-mem paralllelisation of scale)
        call cline%delete('autoscale')
        ! whether to perform perform ab-initio reconstruction with e/o class averages
        l_lpset = cline%defined('lpstart') .and. cline%defined('lpstop')
        ! make master parameters
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! from now on we are in the ptcl3D segment, final report is in the cls3D segment
        call cline%set('oritype', 'ptcl3D')
        ! state string
        str_state = int2str_pad(1,2)
        ! decide wether to search for the symmetry axis
        pgrp_init    = trim(params%pgrp_start)
        pgrp_refine  = trim(params%pgrp)
        srch4symaxis = trim(pgrp_refine) .ne. trim(pgrp_init)
        symran_before_refine = .false.
        if( pgrp_init.ne.'c1' .or. pgrp_refine.ne.'c1' )then
            se1 = sym(pgrp_init)
            se2 = sym(pgrp_refine)
            if(se1%get_nsym() > se2%get_nsym())then
                ! ensure se2 is a subgroup of se1
                if( .not. se1%has_subgrp(pgrp_refine) )&
                    &THROW_HARD('Incompatible symmetry groups; simple_commander_abinitio')
                ! set flag for symmetry randomisation before refinmement
                ! in case we are moving from a higher to lower group
                symran_before_refine = .true.
            else if(se2%get_nsym() > se1%get_nsym())then
                ! ensure se1 is a subgroup of se2
                if( .not. se2%has_subgrp(pgrp_init) )&
                    &THROW_HARD('Incompatible symmetry groups; simple_commander_abinitio')
            endif
        endif
        ! read project & update sampling distance
        call spproj%read(params%projfile)
        ! retrieve cavgs stack & FRCS info
        call spproj%get_cavgs_stk(stk, ncavgs, orig_smpd)
        ctfvars%smpd = orig_smpd
        params%smpd  = orig_smpd
        orig_stk     = stk
        if( .not.spproj%os_cls2D%isthere('state') )then
            ! start from import
            allocate(states(ncavgs), source=1)
        else
            ! start from previous 2D
            states = nint(spproj%os_cls2D%get_all('state'))
        endif
        if( count(states==0) .eq. ncavgs )then
            THROW_HARD('no class averages detected in project file: '//trim(params%projfile)//'; initial_3Dmodel')
        endif
        ! SANITY CHECKS
        ! e/o
        if( l_lpset )then
            ! no filtering
        else
            call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
            if( .not.file_exists(frcs_fname) )then
                THROW_HARD('the project file does not contain enough information for e/o alignment, use a low-pass instead: LPSTART/LPSTOP')
            endif
        endif
        ! set lplims
        lplims(1) = 20.
        lplims(2) = 8.
        if( l_lpset )then
            lplims(1) = params%lpstart
            lplims(2) = params%lpstop
        else
            if( cline%defined('lpstart') )then
                lplims(1) = params%lpstart
            else
                tmp_rarr  = spproj%os_cls2D%get_all('res')
                tmp_iarr  = nint(spproj%os_cls2D%get_all('state'))
                res       = pack(tmp_rarr, mask=(tmp_iarr>0))
                lplims(1) = max(median_nocopy(res), lplims(2))
                deallocate(res, tmp_iarr, tmp_rarr)
            endif
        endif
        ! prepare a temporary project file for the class average processing
        allocate(WORK_PROJFILE, source=trim(ORIG_WORK_PROJFILE))
        call del_file(WORK_PROJFILE)
        work_proj1%projinfo  = spproj%projinfo
        work_proj1%compenv   = spproj%compenv
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
        orig_box      = work_proj1%get_box()
        orig_msk      = params%msk
        smpd_target   = max(params%smpd, lplims(2)*LP2SMPDFAC)
        do_autoscale  = do_autoscale .and. smpd_target > work_proj1%get_smpd()
        scale_factor1 = 1.
        if( do_autoscale )then
            deallocate(WORK_PROJFILE)
            call simple_mkdir(STKPARTSDIR,errmsg="commander_hlev_wflows :: exec_initial_3Dmodel;  ")
            call work_proj1%scale_projfile(smpd_target, WORK_PROJFILE, cline, cline_scale1, dir=trim(STKPARTSDIR))
            scale_factor1 = cline_scale1%get_rarg('scale')
            box           = nint(cline_scale1%get_rarg('newbox'))
            msk           = cline%get_rarg('msk')
            call cline_scale1%delete('smpd')
            call xscale_distr%execute( cline_scale1 )
        else
            box = orig_box
            msk = orig_msk
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
        ! in stage 2 it follows optional user input and defaults to cc
        call cline_refine3D_snhc%set('objfun', 'cc')
        call cline_refine3D_init%set('objfun', 'cc')
        ! reconstruct3D & project are not distributed executions, so remove the nparts flag
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
        call cline_refine3D_snhc%set('match_filt', 'no')
        call cline_refine3D_snhc%set('ptclw',      'no')  ! no soft particle weights in first phase
        call cline_refine3D_snhc%delete('update_frac') ! no fractional update in first phase
        ! (2) REFINE3D_INIT
        call cline_refine3D_init%set('prg',      'refine3D')
        call cline_refine3D_init%set('projfile', trim(WORK_PROJFILE))
        call cline_refine3D_init%set('msk',      msk)
        call cline_refine3D_init%set('box',      real(box))
        call cline_refine3D_init%set('maxits',   real(MAXITS_INIT))
        call cline_refine3D_init%set('vol1',     trim(SNHCVOL)//trim(str_state)//params%ext)
        call cline_refine3D_init%set('lp',       lplims(1))
        call cline_refine3D_init%set('match_filt','no')
        call cline_refine3D_init%set('ptclw',     'no')  ! no soft particle weights in init phase
        if( .not. cline_refine3D_init%defined('nspace') )then
            call cline_refine3D_init%set('nspace', real(NSPACE_INIT))
        endif
        ! (3) SYMMETRY AXIS SEARCH
        if( srch4symaxis )then
            ! need to replace original point-group flag with c1/pgrp_start
            call cline_refine3D_snhc%set('pgrp', trim(pgrp_init))
            call cline_refine3D_init%set('pgrp', trim(pgrp_init))
            ! symsrch
            call qenv%new(1, exec_bin='simple_exec')
            call cline_symsrch%set('prg',     'symaxis_search') ! needed for cluster exec
            call cline_symsrch%set('pgrp',     trim(pgrp_refine))
            call cline_symsrch%set('msk',      msk)
            call cline_symsrch%set('smpd',     work_proj1%get_smpd())
            call cline_symsrch%set('projfile', trim(WORK_PROJFILE))
            if( .not. cline_symsrch%defined('cenlp') ) call cline_symsrch%set('cenlp', CENLP)
            call cline_symsrch%set('hp',       params%hp)
            call cline_symsrch%set('lp',       lplims(1))
            call cline_symsrch%set('oritype',  'ptcl3D')
        endif
        ! (4) REFINE3D REFINE STEP
        call cline_refine3D_refine%set('prg',      'refine3D')
        call cline_refine3D_refine%set('pgrp',     trim(pgrp_refine))
        call cline_refine3D_refine%set('maxits',   real(MAXITS_REFINE))
        call cline_refine3D_refine%set('refine',   'single')
        call cline_refine3D_refine%set('trs',      real(MINSHIFT)) ! activates shift search
        if( l_lpset )then
            call cline_refine3D_refine%set('lp', lplims(2))
        else
            call cline_refine3D_refine%delete('lp')
            call cline_refine3D_refine%set('lplim_crit',  0.5)
            call cline_refine3D_refine%set('lpstop',      lplims(2))
            call cline_refine3D_refine%set('clsfrcs',    'yes')
            call cline_refine3D_refine%set('match_filt', 'yes')
        endif
        if( .not. cline_refine3D_refine%defined('nspace') )then
            call cline_refine3D_refine%set('nspace', real(NSPACE_REFINE))
        endif
        ! (5) RE-CONSTRUCT & RE-PROJECT VOLUME
        call cline_reconstruct3D%set('prg',     'reconstruct3D')
        call cline_reconstruct3D%set('msk',      orig_msk)
        call cline_reconstruct3D%set('box',      real(orig_box))
        call cline_reconstruct3D%set('projfile', ORIG_WORK_PROJFILE)
        call cline_postprocess%set('prg',       'postprocess')
        call cline_postprocess%set('projfile',   ORIG_WORK_PROJFILE)
        call cline_postprocess%set('mkdir',      'no')
        if( l_lpset )then
            call cline_postprocess%set('lp', lplims(2))
        else
            call cline_postprocess%delete('lp')
        endif
        call cline_reproject%set('prg',   'reproject')
        call cline_reproject%set('pgrp',   trim(pgrp_refine))
        call cline_reproject%set('outstk','reprojs'//params%ext)
        call cline_reproject%set('smpd',   params%smpd)
        call cline_reproject%set('msk',    orig_msk)
        call cline_reproject%set('box',    real(orig_box))
        ! execute commanders
        if( .not. cline%defined('vol1') )then
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> INITIALIZATION WITH STOCHASTIC NEIGHBORHOOD HILL-CLIMBING'
            write(logfhandle,'(A,F6.1,A)') '>>> LOW-PASS LIMIT FOR ALIGNMENT: ', lplims(1),' ANGSTROMS'
            write(logfhandle,'(A)') '>>>'
            call xrefine3D_distr%execute(cline_refine3D_snhc)
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> INITIAL 3D MODEL GENERATION WITH REFINE3D'
            write(logfhandle,'(A)') '>>>'
            call xrefine3D_distr%execute(cline_refine3D_init)
            iter     = cline_refine3D_init%get_rarg('endit')
            vol_iter = trim(VOL_FBODY)//trim(str_state)//params%ext
            if( symran_before_refine )then
                call work_proj1%read_segment('ptcl3D', trim(WORK_PROJFILE))
                call se1%symrandomize(work_proj1%os_ptcl3D)
                call work_proj1%write_segment_inside('ptcl3D', trim(WORK_PROJFILE))
            endif
            if( srch4symaxis )then
                write(logfhandle,'(A)') '>>>'
                write(logfhandle,'(A)') '>>> SYMMETRY AXIS SEARCH'
                write(logfhandle,'(A)') '>>>'
                call cline_symsrch%set('vol1', trim(vol_iter))
                if( qenv%get_qsys() .eq. 'local' )then
                    call xsymsrch%execute(cline_symsrch)
                else
                    call qenv%exec_simple_prg_in_queue(cline_symsrch, 'SYMAXIS_SEARCH_FINISHED')
                endif
                call del_file('SYMAXIS_SEARCH_FINISHED')
            endif
        else
            iter = 0
            call cline_refine3D_refine%set('lp', lplims(2))
        endif
        ! prep refinement stage
        call work_proj1%read_segment('ptcl3D', trim(WORK_PROJFILE))
        os = work_proj1%os_ptcl3D
        ! modulate shifts
        if( do_autoscale )then
            call os%mul_shifts( 1./scale_factor1 )
            ! clean stacks & project file & o_peaks on disc
            call work_proj1%read_segment('stk', trim(WORK_PROJFILE))
            do istk=1,work_proj1%os_stk%get_noris()
                call work_proj1%os_stk%getter(istk, 'stk', stk)
                call del_file(trim(stk))
            enddo
        endif
        call work_proj1%kill()
        call del_file(WORK_PROJFILE)
        deallocate(WORK_PROJFILE)
        call del_files(O_PEAKS_FBODY, params_glob%nparts, ext=BIN_EXT)
        ! re-create project
        call del_file(ORIG_WORK_PROJFILE)
        work_proj2%projinfo = spproj%projinfo
        work_proj2%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris()  > 0 ) work_proj2%jobproc = spproj%jobproc
        if( l_lpset )then
            call work_proj2%add_stk(trim(orig_stk), ctfvars)
            work_proj2%os_ptcl3D = os
            call work_proj2%os_ptcl3D%set_all('state', real(states))
        else
            call prep_eo_stks_refine
            params_glob%nptcls = work_proj2%get_nptcls()
        endif
        call os%kill
        ! renaming
        allocate(WORK_PROJFILE, source=trim(ORIG_WORK_PROJFILE))
        call work_proj2%projinfo%delete_entry('projname')
        call work_proj2%projinfo%delete_entry('projfile')
        call cline%set('projfile', trim(WORK_PROJFILE))
        call cline%set('projname', trim(get_fbody(trim(WORK_PROJFILE),trim('simple'))))
        call work_proj2%update_projinfo(cline)
        call work_proj2%write
        ! split
        if( l_lpset )then
            if(params%nparts == 1)then
                ! all good
            else
                call work_proj2%split_stk(params%nparts)
            endif
        endif
        ! refinement scaling
        scale_factor2 = 1.0
        if( do_autoscale )then
            if( scale_factor1 < SCALEFAC2_TARGET )then
                smpd_target = orig_smpd / SCALEFAC2_TARGET
                call cline%set('msk',orig_msk)
                call work_proj2%scale_projfile(smpd_target, WORK_PROJFILE, cline, cline_scale2, dir=trim(STKPARTSDIR))
                scale_factor2 = cline_scale2%get_rarg('scale')
                box = nint(cline_scale2%get_rarg('newbox'))
                msk = cline%get_rarg('msk')
                call cline_scale2%delete('smpd') !!
                call xscale_distr%execute( cline_scale2 )
                call work_proj2%os_ptcl3D%mul_shifts(scale_factor2)
                call work_proj2%write
                if( .not.l_lpset ) call rescale_2Dfilter
            else
                do_autoscale = .false.
                box = orig_box
                msk = orig_msk
            endif
        endif
        call cline_refine3D_refine%set('msk', msk)
        call cline_refine3D_refine%set('box', real(box))
        call cline_refine3D_refine%set('projfile', WORK_PROJFILE)
        ! refinement stage
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> PROBABILISTIC REFINEMENT'
        write(logfhandle,'(A)') '>>>'
        call cline_refine3D_refine%set('startit', iter + 1.)
        call xrefine3D_distr%execute(cline_refine3D_refine)
        iter = cline_refine3D_refine%get_rarg('endit')
        ! updates shifts & deals with final volume
        call work_proj2%read_segment('ptcl3D', WORK_PROJFILE)
        if( do_autoscale )then
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
            write(logfhandle,'(A)') '>>>'
            ! modulates shifts
            os = work_proj2%os_ptcl3D
            call os%mul_shifts(1./scale_factor2)
            call work_proj2%kill
            call work_proj2%read_segment('ptcl3D', ORIG_WORK_PROJFILE)
            work_proj2%os_ptcl3D = os
            call work_proj2%write_segment_inside('ptcl3D', ORIG_WORK_PROJFILE)
            ! reconstruction
            call xreconstruct3D_distr%execute(cline_reconstruct3D)
            vol_iter = trim(VOL_FBODY)//trim(str_state)//params%ext
            ! because postprocess only updates project file when mkdir=yes
            call work_proj2%read_segment('out', ORIG_WORK_PROJFILE)
            call work_proj2%add_vol2os_out(vol_iter, params%smpd, 1, 'vol')
            if( .not.l_lpset )then
                call work_proj2%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), 1, orig_box)
                call work_proj2%add_vol2os_out(ANISOLP_FBODY//str_state//params%ext, orig_smpd, 1, 'vol_filt', box=orig_box)
            endif
            call work_proj2%write_segment_inside('out',ORIG_WORK_PROJFILE)
            call xpostprocess%execute(cline_postprocess)
            call os%kill
        else
            iter     = cline_refine3D_refine%get_rarg('endit')
            vol_iter = trim(VOL_FBODY)//trim(str_state)//params%ext
            call vol%new([orig_box,orig_box,orig_box],orig_smpd)
            call vol%read(vol_iter)
            call vol%mirror('x')
            call vol%write(add2fbody(vol_iter,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX)))
            call vol%kill
        endif
        status = simple_rename(vol_iter, trim(REC_FBODY)//params%ext)
        status = simple_rename(add2fbody(vol_iter,params%ext,PPROC_SUFFIX),&
            &trim(REC_PPROC_FBODY)//params%ext)
        status = simple_rename(add2fbody(vol_iter,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX)),&
            &trim(REC_PPROC_MIRR_FBODY)//params%ext)
        ! updates original cls3D segment
        call work_proj2%os_ptcl3D%delete_entry('stkind')
        call work_proj2%os_ptcl3D%delete_entry('eo')
        params_glob%nptcls = ncavgs
        if( l_lpset )then
            spproj%os_cls3D = work_proj2%os_ptcl3D
        else
            call spproj%os_cls3D%new(ncavgs)
            do icls=1,ncavgs
                call work_proj2%os_ptcl3D%get_ori(icls, o_tmp)
                call spproj%os_cls3D%set_ori(icls, o_tmp)
            enddo
            call conv_eo(work_proj2%os_ptcl3D)
        endif
        call work_proj2%kill
        ! revert splitting
        call spproj%os_cls3D%set_all2single('stkind',1.)
        ! map the orientation parameters obtained for the clusters back to the particles
        call spproj%map2ptcls
        ! add rec_final to os_out
        call spproj%add_vol2os_out(trim(REC_FBODY)//params%ext, params%smpd, 1, 'vol_cavg')
        ! write results (this needs to be a full write as multiple segments are updated)
        call spproj%write()
        ! reprojections
        call spproj%os_cls3D%write('final_oris.txt')
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RE-PROJECTION OF THE FINAL VOLUME'
        write(logfhandle,'(A)') '>>>'
        call cline_reproject%set('vol1',   trim(REC_PPROC_FBODY)//params%ext)
        call cline_reproject%set('oritab', 'final_oris.txt')
        call xreproject%execute(cline_reproject)
        ! write alternated stack
        call img%new([orig_box,orig_box,1], orig_smpd)
        cnt = -1
        do icls=1,ncavgs
            cnt = cnt + 2
            call img%read(orig_stk,icls)
            call img%norm
            call img%write('cavgs_reprojs.mrc',cnt)
            call img%read('reprojs.mrc',icls)
            call img%norm
            call img%write('cavgs_reprojs.mrc',cnt+1)
        enddo
        ! end gracefully
        call se1%kill
        call se2%kill
        call img%kill
        call spproj%kill
        call o_tmp%kill
        if( allocated(WORK_PROJFILE) ) call del_file(WORK_PROJFILE)
        call del_file(ORIG_WORK_PROJFILE)
        call simple_rmdir(STKPARTSDIR)
        call simple_end('**** SIMPLE_INITIAL_3DMODEL NORMAL STOP ****')

        contains

            subroutine prep_eo_stks_refine
                use simple_ori, only: ori
                type(ori)                     :: o, o_even, o_odd
                character(len=:), allocatable :: eostk, ext
                integer :: even_ind, odd_ind, state, icls
                call os%delete_entry('lp')
                call cline_refine3D_refine%set('frcs',frcs_fname)
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
                    call os%get_ori(icls, o)
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

            subroutine rescale_2Dfilter
                use simple_projection_frcs, only: projection_frcs
                type(projection_frcs) :: projfrcs, projfrcs_sc
                call projfrcs%read(frcs_fname)
                call projfrcs%downsample(box, projfrcs_sc)
                frcs_fname = trim(FRCS_FILE)
                call projfrcs_sc%write(frcs_fname)
                call cline_refine3D_refine%set('frcs',frcs_fname)
                call projfrcs%kill
                call projfrcs_sc%kill
            end subroutine rescale_2Dfilter

            subroutine conv_eo( os )
                use simple_ori, only: ori
                class(oris), intent(inout) :: os
                type(sym) :: se
                type(ori) :: o_odd, o_even
                real      :: avg_euldist, euldist
                integer   :: icls, ncls
                call se%new(pgrp_refine)
                avg_euldist = 0.
                ncls = 0
                do icls=1,os%get_noris()/2
                    call os%get_ori(icls, o_even)
                    if( o_even%get_state() == 0 )cycle
                    ncls    = ncls + 1
                    call os%get_ori(ncavgs+icls, o_odd)
                    euldist = rad2deg(o_odd.euldist.o_even)
                    if( se%get_nsym() > 1 )then
                        call o_odd%mirror2d
                        call se%rot_to_asym(o_odd)
                        euldist = min(rad2deg(o_odd.euldist.o_even), euldist)
                    endif
                    avg_euldist = avg_euldist + euldist
                enddo
                avg_euldist = avg_euldist/real(ncls)
                write(logfhandle,'(A)')'>>>'
                write(logfhandle,'(A,F6.1)')'>>> EVEN/ODD AVERAGE ANGULAR DISTANCE: ', avg_euldist
            end subroutine conv_eo

    end subroutine exec_initial_3Dmodel

end module simple_commander_abinitio
