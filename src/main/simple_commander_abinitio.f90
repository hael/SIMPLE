! concrete commander: high-level workflows
module simple_commander_abinitio
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters, params_glob
use simple_sp_project,         only: sp_project
use simple_stack_io,           only: stack_io
use simple_qsys_env,           only: qsys_env
use simple_commander_base,     only: commander_base
use simple_commander_volops,   only: reproject_commander, symaxis_search_commander, postprocess_commander
use simple_commander_rec,      only: reconstruct3D_commander, reconstruct3D_commander_distr
use simple_commander_refine3D, only: refine3D_commander, refine3D_commander_distr
use simple_commander_project,  only: scale_project_commander_distr
use simple_commander_imgproc,  only: scale_commander
use simple_procimgstk,         only: shift_imgfile
use simple_oris,               only: oris
use simple_ori,                only: ori
use simple_image,              only: image
use simple_sym,                only: sym
use simple_builder,            only: builder
use simple_opt_filter,         only: opt_2D_filter_sub
use simple_masker,             only: automask2D
use simple_qsys_funs
use simple_estimate_ssnr
implicit none

public :: initial_3Dmodel_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: initial_3Dmodel_commander
    contains
    procedure :: execute => exec_initial_3Dmodel
end type initial_3Dmodel_commander

contains

    !> for generation of an initial 3d model from class averages
    subroutine exec_initial_3Dmodel( self, cline )
        class(initial_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! constants
        real,                  parameter :: SCALEFAC2_TARGET = 0.5, AMSKLP3D = 15.
        real,                  parameter :: CENLP_DEFAULT = 30.
        integer,               parameter :: MAXITS_SNHC=20, MAXITS_INIT=15, MAXITS_REFINE=40, WINSZ_AUTOMSK = 5
        integer,               parameter :: NSPACE_SNHC=1000, NSPACE_INIT=1000, NSPACE_REFINE=2500
        character(len=STDLEN), parameter :: ORIG_work_projfile   = 'initial_3Dmodel_tmpproj.simple'
        character(len=STDLEN), parameter :: REC_FBODY            = 'rec_final'
        character(len=STDLEN), parameter :: REC_PPROC_FBODY      = trim(REC_FBODY)//trim(PPROC_SUFFIX)
        character(len=STDLEN), parameter :: REC_PPROC_MIRR_FBODY = trim(REC_PPROC_FBODY)//trim(MIRR_SUFFIX)
        ! distributed commanders
        type(refine3D_commander_distr)      :: xrefine3D_distr
        type(scale_project_commander_distr) :: xscale
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        ! shared-mem commanders
        type(refine3D_commander)            :: xrefine3D
        type(reconstruct3D_commander)       :: xreconstruct3D
        type(symaxis_search_commander)      :: xsymsrch
        type(reproject_commander)           :: xreproject
        type(postprocess_commander)         :: xpostprocess
        type(scale_commander)               :: xscale_msk
        ! command lines
        type(cmdline) :: cline_refine3D_snhc, cline_refine3D_init, cline_refine3D_refine
        type(cmdline) :: cline_symsrch
        type(cmdline) :: cline_reconstruct3D, cline_postprocess
        type(cmdline) :: cline_reproject
        type(cmdline) :: cline_scale1, cline_scale2, cline_scale_msk
        ! other
        character(len=:), allocatable :: stk, orig_stk, frcs_fname, shifted_stk, stk_even, stk_odd, ext
        character(len=:), allocatable :: work_projfile
        real,             allocatable :: res(:), tmp_rarr(:), diams(:)
        integer,          allocatable :: states(:), tmp_iarr(:)
        type(image),      allocatable :: cavgs_eo(:,:), masks(:)
        class(parameters), pointer    :: params_ptr => null()
        character(len=2)      :: str_state
        type(qsys_env)        :: qenv
        type(parameters)      :: params
        type(ctfparams)       :: ctfvars ! ctf=yes by default
        type(sp_project)      :: spproj, work_proj1, work_proj2
        type(oris)            :: os
        type(ori)             :: o_tmp
        type(sym)             :: se1,se2
        type(image)           :: img, vol
        type(stack_io)        :: stkio_r, stkio_r2, stkio_w
        character(len=STDLEN) :: vol_iter, pgrp_init, pgrp_refine, vol_iter_pproc, vol_iter_pproc_mirr
        real                  :: iter, smpd_target, lplims(2), orig_smpd, cenlp
        real                  :: scale_factor1, scale_factor2, lp3(3)
        integer               :: icls, ncavgs, orig_box, box, istk, cnt, ifoo, ldim(3)
        logical               :: srch4symaxis, do_autoscale, symran_before_refine, l_lpset, l_shmem
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('ptclw')     ) call cline%set('ptclw',      'no')
        if( .not. cline%defined('overlap')   ) call cline%set('overlap',     0.8)
        if( .not. cline%defined('fracsrch')  ) call cline%set('fracsrch',    0.9)
        if( .not. cline%defined('envfsc')    ) call cline%set('envfsc',     'no')
        if( .not. cline%defined('ngrow')     ) call cline%set('ngrow',        3.)
        if( .not. cline%defined('amsklp')    ) call cline%set('amsklp',      20.)
        if( .not. cline%defined('edge')      ) call cline%set('edge',         6.)
        ! call set_automask2D_defaults(cline)
        ! set shared-memory flag
        if( cline%defined('nparts') )then
            if( nint(cline%get_rarg('nparts')) == 1 )then
                l_shmem = .true.
                call cline%delete('nparts')
            else
                l_shmem = .false.
            endif
        else
            l_shmem = .true.
        endif
        ! hard set oritype
        call cline%set('oritype', 'out') ! because cavgs are part of out segment
        ! hard set bfactor
        call cline%set('bfac', 0.)       ! because initial models should not be sharpened
        ! class averages, so no CTF
        ctfvars%ctfflag = CTFFLAG_NO
        ! auto-scaling prep
        do_autoscale = (cline%get_carg('autoscale').eq.'yes')
        ! remove autoscale flag from command line, since no scaled partial stacks
        ! will be produced (this program always uses shared-mem parallelisation of scale)
        call cline%delete('autoscale')
        ! whether to perform perform ab-initio reconstruction with e/o class averages
        l_lpset = cline%defined('lpstart') .and. cline%defined('lpstop')
        ! make master parameters
        call params%new(cline)
        ! take care of automask flag
        if( cline%defined('automsk') ) call cline%delete('automsk')
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
                if( .not. se1%has_subgrp(pgrp_refine) ) THROW_HARD('Incompatible symmetry groups; simple_commander_abinitio')
                ! set flag for symmetry randomisation before refinmement
                ! in case we are moving from a higher to lower group
                symran_before_refine = .true.
            else if( se2%get_nsym() > se1%get_nsym() )then
                ! ensure se1 is a subgroup of se2
                if( .not. se2%has_subgrp(pgrp_init) ) THROW_HARD('Incompatible symmetry groups; simple_commander_abinitio')
            endif
        endif
        ! read project & update sampling distance
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        ! retrieve cavgs stack & FRCS info
        call spproj%get_cavgs_stk(stk, ncavgs, orig_smpd)
        ext = '.'//fname2ext( stk )
        ! e/o
        if( l_lpset )then
            ! no filtering
        else
            call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
            if( .not.file_exists(frcs_fname) )then
                THROW_HARD('the project file does not contain enough information for e/o alignment, use a low-pass instead: LPSTART/LPSTOP')
            endif
            ! update params
            params%frcs = trim(frcs_fname)
        endif
        if( params%l_nonuniform )then
            ! retrieve even/odd stack names
            stk_even = add2fbody(stk, ext, '_even')
            stk_odd  = add2fbody(stk, ext, '_odd')
            if( .not. file_exists(trim(stk_even)) ) THROW_HARD('Even stack '//trim(stk_even)//' does not exist!')
            if( .not. file_exists(trim(stk_odd))  ) THROW_HARD('Odd stack '//trim(stk_even)//' does not exist!')
            ! read even/odd class averages into array
            allocate( cavgs_eo(ncavgs,2), masks(ncavgs) )
            call find_ldim_nptcls(stk_even, ldim, ifoo)
            ldim(3) = 1 ! class averages
            call stkio_r%open(stk_odd,   orig_smpd, 'read', bufsz=500)
            call stkio_r2%open(stk_even, orig_smpd, 'read', bufsz=500)
            do icls = 1, ncavgs
                call cavgs_eo(icls,1)%new(ldim, orig_smpd)
                call cavgs_eo(icls,2)%new(ldim, orig_smpd)
                call stkio_r%read(icls, cavgs_eo(icls,1))
                call stkio_r2%read(icls, cavgs_eo(icls,2))
                call masks(icls)%copy(cavgs_eo(icls,1))
                call masks(icls)%add(cavgs_eo(icls,2))
                call masks(icls)%mul(0.5)
            end do
            call stkio_r%close
            call stkio_r2%close
            ! nonuniform filtering
            if( params%l_automsk )then
                call automask2D(masks, params%ngrow, WINSZ_AUTOMSK, params%edge, diams)
                call opt_2D_filter_sub(cavgs_eo(:,2), cavgs_eo(:,1), masks)
            else
                call opt_2D_filter_sub(cavgs_eo(:,2), cavgs_eo(:,1))
            endif
            ! write even filtered cavgs
            stk_even = 'cavgs_nonuniform_even.mrc'
            call stkio_w%open(stk_even, orig_smpd, 'write', box=ldim(1), is_ft=.false., bufsz=500)
            do icls = 1, ncavgs
                call stkio_w%write(icls, cavgs_eo(icls,2))
            end do
            call stkio_w%close
            ! write odd filtered cavgs
            stk_odd = 'cavgs_nonuniform_odd.mrc'
            call stkio_w%open(stk_odd, orig_smpd, 'write', box=ldim(1), is_ft=.false., bufsz=500)
            do icls = 1, ncavgs
                call stkio_w%write(icls, cavgs_eo(icls,1))
            end do
            call stkio_w%close
            ! write merged filtered cavgs
            stk = 'cavgs_nonuniform.mrc'
            call stkio_w%open(stk, orig_smpd, 'write', box=ldim(1), is_ft=.false., bufsz=500)
            do icls = 1, ncavgs
                call cavgs_eo(icls,1)%add(cavgs_eo(icls,2))
                call cavgs_eo(icls,1)%mul(0.5)
                call stkio_w%write(icls, cavgs_eo(icls,1))
            end do
            call stkio_w%close
            do icls = 1, ncavgs
                call cavgs_eo(icls,1)%kill
                call cavgs_eo(icls,2)%kill
                call masks(icls)%kill
            end do
            deallocate(cavgs_eo, masks)
        else
            stk_even = add2fbody(trim(stk), trim(ext), '_even')
            stk_odd  = add2fbody(trim(stk), trim(ext), '_odd')
        endif
        params%amsklp = AMSKLP3D
        call cline%set('amsklp', AMSKLP3D)
        ctfvars%smpd = orig_smpd
        params%smpd  = orig_smpd
        orig_stk     = stk
        shifted_stk  = basename(add2fbody(stk, ext, '_shifted'))
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
        ! set lplims
        call mskdiam2lplimits(params%mskdiam, lplims(1), lplims(2), cenlp)
        if( .not. cline%defined('cenlp') ) params_glob%cenlp = cenlp
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
                call hpsort(res)
                ! new way
                lplims(2) = max(median_nocopy(res(:3)), lplims(2)) ! low-pass limit is median of three best (as in 2D)
                deallocate(res, tmp_iarr, tmp_rarr)
            endif
        endif
        write(logfhandle,'(A,F5.1)') '>>> DID SET STARTING  LOW-PASS LIMIT (IN A) TO: ', lplims(1)
        write(logfhandle,'(A,F5.1)') '>>> DID SET HARD      LOW-PASS LIMIT (IN A) TO: ', lplims(2)
        write(logfhandle,'(A,F5.1)') '>>> DID SET CENTERING LOW-PASS LIMIT (IN A) TO: ', params_glob%cenlp
        ! prepare a temporary project file for the class average processing
        allocate(work_projfile, source=trim(ORIG_work_projfile))
        call del_file(work_projfile)
        work_proj1%projinfo  = spproj%projinfo
        work_proj1%compenv   = spproj%compenv
        if( spproj%jobproc%get_noris()  > 0 ) work_proj1%jobproc = spproj%jobproc
        call work_proj1%add_stk(trim(stk), ctfvars)
        call work_proj1%os_ptcl3D%set_all('state', real(states)) ! takes care of states
        ! name change
        call work_proj1%projinfo%delete_entry('projname')
        call work_proj1%projinfo%delete_entry('projfile')
        call cline%set('projfile', trim(work_projfile))
        call cline%set('projname', trim(get_fbody(trim(work_projfile),trim('simple'))))
        call work_proj1%update_projinfo(cline)
        call work_proj1%write()
        ! split
        if( .not. l_shmem ) call work_proj1%split_stk(params%nparts)
        ! down-scale
        orig_box      = work_proj1%get_box()
        smpd_target   = max(params%smpd, lplims(2)*LP2SMPDFAC)
        do_autoscale  = do_autoscale .and. smpd_target > work_proj1%get_smpd()
        scale_factor1 = 1.
        if( do_autoscale )then
            deallocate(work_projfile)
            call simple_mkdir(STKPARTSDIR,errmsg="commander_hlev_wflows :: exec_initial_3Dmodel;  ")
            call work_proj1%scale_projfile(smpd_target, work_projfile, cline, cline_scale1, dir=trim(STKPARTSDIR))
            scale_factor1 = cline_scale1%get_rarg('scale')
            box           = nint(cline_scale1%get_rarg('newbox'))
            call cline_scale1%delete('smpd')
            call xscale%execute( cline_scale1 )
        else
            box = orig_box
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
        ! re-project is never a distributed executions, so remove the nparts flag if there
        call cline_reproject%delete('nparts')
        ! initialise command line parameters
        ! (1) INITIALIZATION BY STOCHASTIC NEIGHBORHOOD HILL-CLIMBING
        call cline_refine3D_snhc%set('projfile',   trim(work_projfile))
        call cline_refine3D_snhc%set('box',        real(box))
        call cline_refine3D_snhc%set('prg',        'refine3D')
        call cline_refine3D_snhc%set('refine',     'snhc')
        call cline_refine3D_snhc%set('lp',         lplims(1))
        call cline_refine3D_snhc%set('nspace',     real(NSPACE_SNHC))
        call cline_refine3D_snhc%set('maxits',     real(MAXITS_SNHC))
        call cline_refine3D_snhc%set('match_filt', 'no')
        call cline_refine3D_snhc%set('ptclw',      'no')  ! no soft particle weights in first phase
        call cline_refine3D_snhc%set('silence_fsc','yes') ! no FSC plot printing in snhc phase
        call cline_refine3D_snhc%set('lp_iters',    0.)   ! low-pass limited resolution, no e/o
        call cline_refine3D_snhc%delete('frac')           ! no rejections in first phase
        ! (2) REFINE3D_INIT
        call cline_refine3D_init%set('projfile', trim(work_projfile))
        call cline_refine3D_init%set('box',      real(box))
        call cline_refine3D_init%set('prg',      'refine3D')
        call cline_refine3D_init%set('lp',       lplims(1))
        call cline_refine3D_init%set('lp_iters', 0.)   ! low-pass limited resolution, no e/o
        if( .not. cline_refine3D_init%defined('nspace') )then
            call cline_refine3D_init%set('nspace', real(NSPACE_INIT))
        endif
        call cline_refine3D_init%set('maxits',   real(MAXITS_INIT))
        call cline_refine3D_init%set('match_filt','no')
        call cline_refine3D_init%set('ptclw',     'no')   ! no soft particle weights in init phase
        call cline_refine3D_init%set('silence_fsc','yes') ! no FSC plot printing in 2nd phase
        call cline_refine3D_init%set('vol1',     trim(SNHCVOL)//trim(str_state)//ext)
        call cline_refine3D_init%delete('frac')           ! no rejections in 2nd phase
        ! (3) SYMMETRY AXIS SEARCH
        if( srch4symaxis )then
            ! need to replace original point-group flag with c1/pgrp_start
            call cline_refine3D_snhc%set('pgrp', trim(pgrp_init))
            call cline_refine3D_init%set('pgrp', trim(pgrp_init))
            ! symsrch
            call qenv%new(1, exec_bin='simple_exec')
            call cline_symsrch%set('prg',     'symaxis_search') ! needed for cluster exec
            call cline_symsrch%set('pgrp',     trim(pgrp_refine))
            call cline_symsrch%set('smpd',     work_proj1%get_smpd())
            call cline_symsrch%set('projfile', trim(work_projfile))
            if( .not. cline_symsrch%defined('cenlp') ) call cline_symsrch%set('cenlp', CENLP)
            call cline_symsrch%set('hp',       params%hp)
            call cline_symsrch%set('lp',       lplims(1))
            call cline_symsrch%set('oritype',  'ptcl3D')
        endif
        ! (4) REFINE3D REFINE STEP
        call cline_refine3D_refine%set('prg',      'refine3D')
        call cline_refine3D_refine%set('pgrp',     trim(pgrp_refine))
        call cline_refine3D_refine%set('maxits',   real(MAXITS_REFINE))
        call cline_refine3D_refine%set('refine',   'shc')
        call cline_refine3D_refine%set('trs',      real(MINSHIFT)) ! activates shift search
        if( l_lpset )then
            call cline_refine3D_refine%set('lp', lplims(2))
            call cline_refine3D_refine%set('lp_iters', 0.)   ! low-pass limited resolution, no e/o
        else
            call cline_refine3D_refine%delete('lp')
            call cline_refine3D_refine%set('lp_iters', 0.)   ! no lp, e/o only
            call cline_refine3D_refine%set('lpstop',      lplims(2))
            call cline_refine3D_refine%set('clsfrcs',    'yes')
            call cline_refine3D_refine%set('match_filt', 'yes')
        endif
        if( .not. cline_refine3D_refine%defined('nspace') )then
            call cline_refine3D_refine%set('nspace', real(NSPACE_REFINE))
        endif
        call cline_refine3D_refine%set('nonuniform', 'no') ! done in 2D
        if( params%l_automsk )then
            call cline_refine3D_refine%set('automsk', trim(params%automsk))
            call cline_refine3D_refine%set('amsklp', AMSKLP3D)
        endif
        ! (5) RE-CONSTRUCT & RE-PROJECT VOLUME
        call cline_reconstruct3D%set('prg',     'reconstruct3D')
        call cline_reconstruct3D%set('box',      real(orig_box))
        call cline_reconstruct3D%set('projfile', ORIG_work_projfile)
        call cline_postprocess%set('prg',       'postprocess')
        call cline_postprocess%set('projfile',   ORIG_work_projfile)
        call cline_postprocess%set('mkdir',      'no')
        call cline_postprocess%set('bfac',       0.)
        if( l_lpset )then
            call cline_postprocess%set('lp', lplims(2))
        else
            call cline_postprocess%delete('lp')
        endif
        if( params%l_automsk )then
            call cline_postprocess%set('automsk', trim(params%automsk))
        endif
        call cline_reproject%set('prg',     'reproject')
        call cline_reproject%set('pgrp',    trim(pgrp_refine))
        call cline_reproject%set('outstk',  'reprojs'//ext)
        call cline_reproject%set('smpd',    params%smpd)
        call cline_reproject%set('box',     real(orig_box))
        ! execute commanders
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> INITIALIZATION WITH STOCHASTIC NEIGHBORHOOD HILL-CLIMBING'
        write(logfhandle,'(A,F6.1,A)') '>>> LOW-PASS LIMIT FOR ALIGNMENT: ', lplims(1),' ANGSTROMS'
        write(logfhandle,'(A)') '>>>'
        if( l_shmem )then
            call rec(cline_refine3D_snhc, l_rnd=.true.)
            params_ptr  => params_glob
            params_glob => null()
            call xrefine3D%execute(cline_refine3D_snhc)
            params_glob => params_ptr
            params_ptr  => null()
        else
            call xrefine3D_distr%execute(cline_refine3D_snhc)
        endif
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> INITIAL 3D MODEL GENERATION WITH REFINE3D'
        write(logfhandle,'(A)') '>>>'
        if( l_shmem )then
            params_ptr  => params_glob
            params_glob => null()
            call xrefine3D%execute(cline_refine3D_init)
            params_glob => params_ptr
            params_ptr  => null()
        else
            call xrefine3D_distr%execute(cline_refine3D_init)
        endif
        iter     = cline_refine3D_init%get_rarg('endit')
        vol_iter = trim(VOL_FBODY)//trim(str_state)//ext
        if( .not. file_exists(vol_iter) ) THROW_HARD('input volume to symmetry axis search does not exist')
        if( symran_before_refine )then
            call work_proj1%read_segment('ptcl3D', trim(work_projfile))
            call se1%symrandomize(work_proj1%os_ptcl3D)
            call work_proj1%write_segment_inside('ptcl3D', trim(work_projfile))
        endif
        if( srch4symaxis )then
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> SYMMETRY AXIS SEARCH'
            write(logfhandle,'(A)') '>>>'
            call cline_symsrch%set('vol1', trim(vol_iter))
            if( l_shmem .or. qenv%get_qsys() .eq. 'local' )then
                call xsymsrch%execute(cline_symsrch)
            else
                call qenv%exec_simple_prg_in_queue(cline_symsrch, 'SYMAXIS_SEARCH_FINISHED')
            endif
            call del_file('SYMAXIS_SEARCH_FINISHED')
        endif
        ! prep refinement stage
        call work_proj1%read_segment('ptcl3D', trim(work_projfile))
        os = work_proj1%os_ptcl3D
        ! modulate shifts
        if( do_autoscale )then
            call os%mul_shifts( 1./scale_factor1 )
            ! clean stacks & project file on disc
            call work_proj1%read_segment('stk', trim(work_projfile))
            do istk=1,work_proj1%os_stk%get_noris()
                call work_proj1%os_stk%getter(istk, 'stk', stk)
                call del_file(trim(stk))
            enddo
        endif
        call work_proj1%kill()
        call del_file(work_projfile)
        deallocate(work_projfile)
        ! re-create project
        call del_file(ORIG_work_projfile)
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
        allocate(work_projfile, source=trim(ORIG_work_projfile))
        call work_proj2%projinfo%delete_entry('projname')
        call work_proj2%projinfo%delete_entry('projfile')
        call cline%set('projfile', trim(work_projfile))
        call cline%set('projname', trim(get_fbody(trim(work_projfile),trim('simple'))))
        call work_proj2%update_projinfo(cline)
        call work_proj2%write
        ! split
        if( l_lpset )then
            if( l_shmem )then
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
                call work_proj2%scale_projfile(smpd_target, work_projfile, cline, cline_scale2, dir=trim(STKPARTSDIR))
                scale_factor2 = cline_scale2%get_rarg('scale')
                box = nint(cline_scale2%get_rarg('newbox'))
                call cline_scale2%delete('smpd')
                call xscale%execute( cline_scale2 )
                call work_proj2%os_ptcl3D%mul_shifts(scale_factor2)
                call work_proj2%write
                if( .not.l_lpset ) call rescale_2Dfilter
            else
                do_autoscale = .false.
                box = orig_box
            endif
        endif
        call cline_refine3D_refine%set('box', real(box))
        call cline_refine3D_refine%set('projfile', work_projfile)
        ! refinement stage
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> REFINEMENT'
        write(logfhandle,'(A)') '>>>'
        call cline_refine3D_refine%set('startit', iter + 1.)
        if( l_shmem )then
            call rec(cline_refine3D_refine, l_rnd=.false.)
            params_ptr  => params_glob
            params_glob => null()
            call xrefine3D%execute(cline_refine3D_refine)
            params_glob => params_ptr
            params_ptr  => null()
        else
            call xrefine3D_distr%execute(cline_refine3D_refine)
        endif
        iter = cline_refine3D_refine%get_rarg('endit')
        ! updates shifts & deals with final volume
        call work_proj2%read_segment('ptcl3D', work_projfile)
        if( do_autoscale )then
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
            write(logfhandle,'(A)') '>>>'
            if( params%l_automsk )then
                ! scale the mask
                if( .not. file_exists('automask'//trim(ext)) ) THROW_HARD('file '//'automask'//trim(ext)//' does not exist')
                call cline_scale_msk%set('smpd',   smpd_target)
                call cline_scale_msk%set('vol1',   'automask'//trim(ext))
                call cline_scale_msk%set('newbox', real(orig_box))
                call cline_scale_msk%set('outvol', 'automask_scaled'//trim(ext))
                call cline_scale_msk%set('mkdir',  'no')
                call cline_scale_msk%set('nthr',   real(params%nthr))
                call xscale_msk%execute(cline_scale_msk)
                call del_file('automask'//trim(ext))
                call simple_rename('automask_scaled'//trim(ext), 'automask'//trim(ext))
                call cline_reconstruct3D%set('mskfile', 'automask'//trim(ext))
            endif
            ! modulates shifts
            os = work_proj2%os_ptcl3D
            call os%mul_shifts(1./scale_factor2)
            call work_proj2%kill
            call work_proj2%read_segment('ptcl3D', ORIG_work_projfile)
            work_proj2%os_ptcl3D = os
            call work_proj2%write_segment_inside('ptcl3D', ORIG_work_projfile)
            ! reconstruction
            if( l_shmem )then
                params_ptr  => params_glob
                params_glob => null()
                call xreconstruct3D%execute(cline_reconstruct3D)
                params_glob => params_ptr
                params_ptr  => null()
            else
                call xreconstruct3D_distr%execute(cline_reconstruct3D)
            endif
            vol_iter = trim(VOL_FBODY)//trim(str_state)//ext
            ! because postprocess only updates project file when mkdir=yes
            call work_proj2%read_segment('out', ORIG_work_projfile)
            call work_proj2%add_vol2os_out(vol_iter, params%smpd, 1, 'vol')
            if( .not.l_lpset )then
                call work_proj2%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), 1, orig_box)
            endif
            call work_proj2%write_segment_inside('out',ORIG_work_projfile)
            call xpostprocess%execute(cline_postprocess)
            call os%kill
        else
            vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(nint(iter),3)//ext
            call vol%new([orig_box,orig_box,orig_box],orig_smpd)
            call vol%read(vol_iter)
            call vol%mirror('x')
            call vol%write(add2fbody(vol_iter,ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX)))
            call vol%kill
        endif
        vol_iter_pproc      = add2fbody(vol_iter,ext,PPROC_SUFFIX)
        vol_iter_pproc_mirr = add2fbody(vol_iter,ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX))
        if( file_exists(vol_iter)            ) call simple_rename(vol_iter,            trim(REC_FBODY)//ext)
        if( file_exists(vol_iter_pproc)      ) call simple_rename(vol_iter_pproc,      trim(REC_PPROC_FBODY)//ext)
        if( file_exists(vol_iter_pproc_mirr) ) call simple_rename(vol_iter_pproc_mirr, trim(REC_PPROC_MIRR_FBODY)//ext)
        ! updates original cls3D segment
        call work_proj2%os_ptcl3D%delete_entry('stkind')
        call work_proj2%os_ptcl3D%delete_entry('eo')
        params_glob%nptcls = ncavgs
        if( l_lpset )then
            spproj%os_cls3D = work_proj2%os_ptcl3D
        else
            call spproj%os_cls3D%new(ncavgs, is_ptcl=.false.)
            do icls=1,ncavgs
                call work_proj2%os_ptcl3D%get_ori(icls, o_tmp)
                call spproj%os_cls3D%set_ori(icls, o_tmp)
            enddo
            call conv_eo(work_proj2%os_ptcl3D)
        endif
        call work_proj2%kill
        ! revert splitting
        call spproj%os_cls3D%set_all2single('stkind', 1.)
        ! map the orientation parameters obtained for the clusters back to the particles
        call spproj%map2ptcls
        ! add rec_final to os_out
        call spproj%add_vol2os_out(trim(REC_FBODY)//ext, params%smpd, 1, 'vol_cavg')
        ! reprojections
        call spproj%os_cls3D%write('final_oris.txt')
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RE-PROJECTION OF THE FINAL VOLUME'
        write(logfhandle,'(A)') '>>>'
        call cline_reproject%set('vol1',   trim(REC_PPROC_FBODY)//ext)
        call cline_reproject%set('oritab', 'final_oris.txt')
        call xreproject%execute(cline_reproject)
        ! write alternated stack
        call img%new([orig_box,orig_box,1], orig_smpd)
        call stkio_r%open(orig_stk,            params%smpd, 'read',                                 bufsz=500)
        call stkio_r2%open('reprojs.mrc',      params%smpd, 'read',                                 bufsz=500)
        call stkio_w%open('cavgs_reprojs.mrc', params%smpd, 'write', box=params%box, is_ft=.false., bufsz=500)
        cnt = -1
        do icls=1,ncavgs
            cnt = cnt + 2
            call stkio_r%read(icls, img)
            call img%norm
            call stkio_w%write(cnt, img)
            call stkio_r2%read(icls, img)
            call img%norm
            call stkio_w%write(cnt + 1, img)
        enddo
        call stkio_r%close
        call stkio_r2%close
        call stkio_w%close
        ! produce shifted stack
        call shift_imgfile(orig_stk, shifted_stk, spproj%os_cls3D, params%smpd)
        ! add shifted stack to project
        call spproj%add_cavgs2os_out(simple_abspath(shifted_stk), params%smpd, 'cavg_shifted')
        ! write results (this needs to be a full write as multiple segments are updated)
        call spproj%write()
        ! end gracefully
        call se1%kill
        call se2%kill
        call img%kill
        call spproj%kill
        call o_tmp%kill
        if( allocated(work_projfile) ) call del_file(work_projfile)
        if( allocated(res)      ) deallocate(res)
        if( allocated(tmp_rarr) ) deallocate(tmp_rarr)
        if( allocated(diams)    ) deallocate(diams)
        if( allocated(states)   ) deallocate(states)
        if( allocated(tmp_iarr) ) deallocate(tmp_iarr)
        call del_file(ORIG_work_projfile)
        call simple_rmdir(STKPARTSDIR)
        call simple_end('**** SIMPLE_INITIAL_3DMODEL NORMAL STOP ****')

        contains

            subroutine rec( cline, l_rnd )
                class(cmdline), intent(inout) :: cline
                logical,        intent(in)    :: l_rnd
                type(parameters) :: params
                type(builder)    :: build
                call build%init_params_and_build_spproj(cline, params)
                if( l_rnd )then
                    call build%spproj%os_ptcl3D%rnd_oris
                    call build%spproj_field%zero_shifts
                    call build%spproj%write_segment_inside('ptcl3D', params%projfile)
                endif
                call cline%set('mkdir', 'no') ! to avoid nested dirs
                params_ptr  => params_glob
                params_glob => null()
                call xreconstruct3D%execute(cline)
                params_glob => params_ptr
                params_ptr  => null()
                call build%spproj_field%kill
                call simple_rename('recvol_state01_even.mrc', 'startvol_even.mrc')
                call simple_rename('recvol_state01_odd.mrc',  'startvol_odd.mrc')
                call simple_rename('recvol_state01.mrc', 'startvol.mrc')
                call cline%set('vol1', 'startvol.mrc')
            end subroutine rec

            subroutine prep_eo_stks_refine
                use simple_ori, only: ori
                type(ori) :: o, o_even, o_odd
                integer   :: even_ind, odd_ind, state, icls
                call os%delete_entry('lp')
                call cline_refine3D_refine%set('frcs',frcs_fname)
                ! add stks
                call work_proj2%add_stk(stk_even, ctfvars)
                call work_proj2%add_stk(stk_odd,  ctfvars)
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
                call o%kill
                call o_even%kill
                call o_odd%kill
            end subroutine prep_eo_stks_refine

            subroutine rescale_2Dfilter
                use simple_class_frcs, only: class_frcs
                type(class_frcs) :: clsfrcs, clsfrcs_sc
                call clsfrcs%read(frcs_fname)
                call clsfrcs%downsample(box, clsfrcs_sc)
                frcs_fname = trim(FRCS_FILE)
                call clsfrcs_sc%write(frcs_fname)
                call cline_refine3D_refine%set('frcs',frcs_fname)
                call clsfrcs%kill
                call clsfrcs_sc%kill
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
