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
use simple_procimgstk,         only: shift_imgfile
use simple_image,              only: image
use simple_builder,            only: builder
use simple_commander_euclid
use simple_euclid_sigma2
use simple_qsys_funs
implicit none

public :: initial_3Dmodel_commander, abinitio_3Dmodel_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: initial_3Dmodel_commander
    contains
    procedure :: execute => exec_initial_3Dmodel
end type initial_3Dmodel_commander

type, extends(commander_base) :: abinitio_3Dmodel_commander
    contains
    procedure :: execute => exec_abinitio_3Dmodel
end type abinitio_3Dmodel_commander

character(len=STDLEN), parameter :: REC_FBODY            = 'rec_final'
character(len=STDLEN), parameter :: REC_PPROC_FBODY      = trim(REC_FBODY)//trim(PPROC_SUFFIX)
character(len=STDLEN), parameter :: REC_PPROC_MIRR_FBODY = trim(REC_PPROC_FBODY)//trim(MIRR_SUFFIX)

contains

    !> for generation of an initial 3d model from class averages
    subroutine exec_initial_3Dmodel( self, cline )
        class(initial_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! constants
        real,                  parameter :: SCALEFAC2_TARGET = 0.5
        real,                  parameter :: CENLP_DEFAULT = 30.
        integer,               parameter :: MAXITS_SNHC=25, MAXITS_INIT=15, MAXITS_REFINE=20
        integer,               parameter :: NSPACE_SNHC=1000, NSPACE_INIT=1000, NSPACE_REFINE=2500
        integer,               parameter :: MINBOX = 64
        character(len=STDLEN), parameter :: ORIG_work_projfile   = 'initial_3Dmodel_tmpproj.simple'
        ! distributed commanders
        type(refine3D_commander_distr)      :: xrefine3D_distr
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        type(calc_pspec_commander_distr)    :: xcalc_pspec_distr
        ! shared-mem commanders
        type(refine3D_commander)            :: xrefine3D
        type(reconstruct3D_commander)       :: xreconstruct3D
        type(symaxis_search_commander)      :: xsymsrch
        type(reproject_commander)           :: xreproject
        type(postprocess_commander)         :: xpostprocess
        ! command lines
        type(cmdline) :: cline_refine3D_snhc, cline_refine3D_init, cline_refine3D_refine
        type(cmdline) :: cline_symsrch
        type(cmdline) :: cline_reconstruct3D, cline_postprocess
        type(cmdline) :: cline_reproject, cline_calc_pspec
        ! other
        character(len=:), allocatable :: stk, orig_stk, frcs_fname, shifted_stk, stk_even, stk_odd, ext
        real,             allocatable :: res(:), tmp_rarr(:), diams(:)
        integer,          allocatable :: states(:), tmp_iarr(:)
        character(len=2)      :: str_state
        type(qsys_env)        :: qenv
        type(parameters)      :: params
        type(ctfparams)       :: ctfvars
        type(sp_project)      :: spproj, work_proj1, work_proj2
        type(sym)             :: se1,se2
        type(image)           :: img, vol
        type(stack_io)        :: stkio_r, stkio_r2, stkio_w
        character(len=STDLEN) :: vol_iter, pgrp_init, pgrp_refine, vol_iter_pproc, vol_iter_pproc_mirr
        character(len=STDLEN) :: sigma2_fname, sigma2_fname_sc, orig_objfun
        real                  :: smpd_target, lplims(2), cenlp
        real                  :: scale_factor1, scale_factor2
        integer               :: icls, ncavgs, cnt, iter, ipart
        logical               :: srch4symaxis, do_autoscale, symran_before_refine, l_lpset, l_shmem, l_euclid
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('amsklp')    ) call cline%set('amsklp',      15.)
        if( .not. cline%defined('envfsc')    ) call cline%set('envfsc',     'no')
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('overlap')   ) call cline%set('overlap',    0.98)
        if( .not. cline%defined('fracsrch')  ) call cline%set('fracsrch',   0.95)
        if( .not. cline%defined('objfun')    ) call cline%set('objfun',     'cc')
        call cline%set('nonuniform', 'no') ! for speed
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
        ! whether to perform perform ab-initio reconstruction with e/o class averages
        l_lpset = cline%defined('lpstart') .and. cline%defined('lpstop')
        ! make master parameters
        call params%new(cline)
        select case(params%cc_objfun)
            case(OBJFUN_EUCLID)
                l_euclid = .true.
            case DEFAULT
                l_euclid = .false.
        end select
        orig_objfun = trim(params%objfun)
        call cline%delete('autoscale')
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! from now on we are in the ptcl3D segment, final report is in the cls3D segment
        call cline%set('oritype', 'ptcl3D')
        ! state string
        str_state    = int2str_pad(1,2)
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
        call spproj%get_cavgs_stk(stk, ncavgs, params%smpd)
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
        stk_even = add2fbody(trim(stk), trim(ext), '_even')
        stk_odd  = add2fbody(trim(stk), trim(ext), '_odd')
        ctfvars%smpd = params%smpd
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
        ! prepare a temporary project file
        work_proj1%projinfo = spproj%projinfo
        work_proj1%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris()  > 0 ) work_proj1%jobproc = spproj%jobproc
        call work_proj1%add_stk(trim(stk), ctfvars)
        call work_proj1%os_ptcl3D%set_all('state', real(states)) ! takes care of states
        ! name change
        call work_proj1%projinfo%delete_entry('projname')
        call work_proj1%projinfo%delete_entry('projfile')
        call cline%set('projfile', trim(ORIG_work_projfile))
        call cline%set('projname', trim(get_fbody(trim(ORIG_work_projfile),trim('simple'))))
        call work_proj1%update_projinfo(cline)
        call work_proj1%write()
        ! split
        if( .not. l_shmem ) call work_proj1%split_stk(params%nparts)
        ! Cropped dimensions
        do_autoscale  = .false.
        scale_factor1 = 1.0
        params%box    = work_proj1%get_box()
        smpd_target   = max(params%smpd, lplims(2)*LP2SMPDFAC)
        params%smpd_crop = params%smpd
        params%box_crop  = params%box
        params%msk_crop  = params%msk
        if( params%l_autoscale )then
            call autoscale(params%box, params%smpd, smpd_target, params%box_crop, params%smpd_crop, scale_factor1, minbox=MINBOX)
            do_autoscale = params%box_crop < params%box
            if( do_autoscale )then
                params%msk_crop = params%msk * scale_factor1
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
            endif
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
        ! re-project is never a distributed execution, so remove the nparts flag if there
        call cline_reproject%delete('nparts')
        ! initialise command line parameters
        ! (1) INITIALIZATION BY STOCHASTIC NEIGHBORHOOD HILL-CLIMBING
        call cline_refine3D_snhc%set('projfile',   trim(ORIG_work_projfile))
        call cline_refine3D_snhc%set('box_crop',   real(params%box_crop))
        call cline_refine3D_snhc%set('smpd_crop',  params%smpd_crop)
        call cline_refine3D_snhc%set('prg',        'refine3D')
        call cline_refine3D_snhc%set('refine',     'snhc')
        call cline_refine3D_snhc%set('lp',         lplims(1))
        call cline_refine3D_snhc%set('nspace',     real(NSPACE_SNHC))
        call cline_refine3D_snhc%set('maxits',     real(MAXITS_SNHC))
        call cline_refine3D_snhc%set('silence_fsc','yes')              ! no FSC plot printing in snhc phase
        call cline_refine3D_snhc%set('lp_iters',    0.)                ! low-pass limited resolution, no e/o
        call cline_refine3D_snhc%delete('frac')                        ! no rejections in first phase
        ! (2) REFINE3D_INIT
        call cline_refine3D_init%set('prg',      'refine3D')
        call cline_refine3D_init%set('projfile', trim(ORIG_work_projfile))
        call cline_refine3D_init%set('box_crop', real(params%box_crop))
        call cline_refine3D_init%set('smpd_crop',params%smpd_crop)
        call cline_refine3D_init%set('refine',   'neigh')
        call cline_refine3D_init%set('lp',       lplims(1))
        call cline_refine3D_init%set('lp_iters', 0.)                   ! low-pass limited resolution, no e/o
        if( .not. cline_refine3D_init%defined('nspace') )then
            call cline_refine3D_init%set('nspace',real(NSPACE_INIT))
        endif
        call cline_refine3D_init%set('maxits',     real(MAXITS_INIT))
        call cline_refine3D_init%set('silence_fsc','yes') ! no FSC plot printing in 2nd phase
        call cline_refine3D_init%set('vol1',       trim(SNHCVOL)//trim(str_state)//ext)
        call cline_refine3D_init%delete('frac')
        ! (3) SYMMETRY AXIS SEARCH
        if( srch4symaxis )then
            ! need to replace original point-group flag with c1/pgrp_start
            call cline_refine3D_snhc%set('pgrp', trim(pgrp_init))
            call cline_refine3D_init%set('pgrp', trim(pgrp_init))
            ! symsrch
            call qenv%new(1, exec_bin='simple_exec')
            call cline_symsrch%set('prg',     'symaxis_search') ! needed for cluster exec
            call cline_symsrch%set('pgrp',     trim(pgrp_refine))
            call cline_symsrch%set('smpd',     params%smpd_crop)
            call cline_symsrch%set('box',      real(params%box_crop))
            call cline_symsrch%set('projfile', trim(ORIG_work_projfile))
            if( .not. cline_symsrch%defined('cenlp') ) call cline_symsrch%set('cenlp', CENLP)
            call cline_symsrch%set('hp',       params%hp)
            call cline_symsrch%set('lp',       lplims(1))
            call cline_symsrch%set('oritype',  'ptcl3D')
        endif
        ! (4) REFINE3D REFINE STEP
        call cline_refine3D_refine%set('prg',    'refine3D')
        call cline_refine3D_refine%set('projfile',   trim(ORIG_work_projfile))
        call cline_refine3D_refine%set('pgrp',   trim(pgrp_refine))
        call cline_refine3D_refine%set('maxits', real(MAXITS_REFINE))
        call cline_refine3D_refine%set('refine', 'neigh')
        call cline_refine3D_refine%set('nspace', real(NSPACE_REFINE))
        if( l_euclid )then
            call cline_refine3D_refine%set('objfun',     orig_objfun)
            call cline_refine3D_refine%set('sigma_est',  'global')
            call cline_reconstruct3D%set('sigma_est',    'global')
            call cline_reconstruct3D%set('needs_sigma',  'yes')
            if( .not. cline%defined('ml_reg') )then
                call cline_refine3D_refine%set('ml_reg', 'yes')
                call cline_reconstruct3D%set('ml_reg',   'yes')
            endif
        else
            call cline_refine3D_refine%set('objfun',     'cc')
            call cline_reconstruct3D%set('needs_sigma',  'no')
            call cline_reconstruct3D%set('ml_reg',       'no')
        endif
        if( l_lpset )then
            call cline_refine3D_refine%set('lp', lplims(2))
            call cline_refine3D_refine%set('lp_iters', real(MAXITS_REFINE)) ! low-pass limited resolution, no e/o
        else
            call cline_refine3D_refine%delete('lp')
            call cline_refine3D_refine%set('lp_iters',      0.)             ! no lp, e/o only
            call cline_refine3D_refine%set('lpstop',      lplims(2))
        endif
        ! (5) RE-CONSTRUCT & RE-PROJECT VOLUME
        call cline_reconstruct3D%set('prg',     'reconstruct3D')
        call cline_reconstruct3D%set('box',      real(params%box))
        call cline_reconstruct3D%set('projfile', ORIG_work_projfile)
        call cline_postprocess%set('prg',       'postprocess')
        call cline_postprocess%set('projfile',   ORIG_work_projfile)
        call cline_postprocess%set('mkdir',      'no')
        call cline_postprocess%delete('bfac') ! sharpen final map
        if( l_lpset )then
            call cline_postprocess%set('lp', lplims(2))
        else
            call cline_postprocess%delete('lp')
        endif
        call cline_reproject%set('prg',    'reproject')
        call cline_reproject%set('pgrp',   trim(pgrp_refine))
        call cline_reproject%set('outstk', 'reprojs'//ext)
        call cline_reproject%set('smpd',   params%smpd)
        call cline_reproject%set('box',    real(params%box))
        ! execute commanders
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> INITIALIZATION WITH STOCHASTIC NEIGHBORHOOD HILL-CLIMBING'
        write(logfhandle,'(A,F6.1,A)') '>>> LOW-PASS LIMIT FOR ALIGNMENT: ', lplims(1),' ANGSTROMS'
        write(logfhandle,'(A)') '>>>'
        if( l_shmem )then
            call rec(cline_refine3D_snhc, l_rnd=.true.)
            call xrefine3D%execute_shmem(cline_refine3D_snhc)
        else
            call xrefine3D_distr%execute(cline_refine3D_snhc)
        endif
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> INITIAL 3D MODEL GENERATION WITH REFINE3D'
        write(logfhandle,'(A)') '>>>'
        if( l_shmem )then
            call xrefine3D%execute_shmem(cline_refine3D_init)
        else
            call xrefine3D_distr%execute(cline_refine3D_init)
        endif
        iter     = nint(cline_refine3D_init%get_rarg('endit'))
        vol_iter = trim(VOL_FBODY)//trim(str_state)//ext
        if( .not. file_exists(vol_iter) ) THROW_HARD('input volume to symmetry axis search does not exist')
        if( symran_before_refine )then
            call work_proj1%read_segment('ptcl3D', ORIG_work_projfile)
            call se1%symrandomize(work_proj1%os_ptcl3D)
            call work_proj1%write_segment_inside('ptcl3D', ORIG_work_projfile)
        endif
        if( srch4symaxis )then
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> SYMMETRY AXIS SEARCH'
            write(logfhandle,'(A)') '>>>'
            call cline_symsrch%set('vol1', trim(vol_iter))
            if( l_shmem .or. qenv%get_qsys() .eq. 'local' )then
                call xsymsrch%execute_shmem(cline_symsrch)
            else
                call qenv%exec_simple_prg_in_queue(cline_symsrch, 'SYMAXIS_SEARCH_FINISHED')
            endif
            call del_file('SYMAXIS_SEARCH_FINISHED')
        endif
        ! Prep refinement stage
        call work_proj1%read_segment('ptcl3D', ORIG_work_projfile)
        work_proj2%projinfo = spproj%projinfo
        work_proj2%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris()  > 0 ) work_proj2%jobproc = spproj%jobproc
        if( l_lpset )then
            call work_proj2%add_stk(trim(orig_stk), ctfvars)
            work_proj2%os_ptcl3D = work_proj1%os_ptcl3D
            call work_proj2%os_ptcl3D%set_all('state', real(states))
        else
            call prep_eo_stks_refine
            params_glob%nptcls = work_proj2%get_nptcls()
        endif
        call work_proj1%kill
        call del_file(ORIG_work_projfile)
        call work_proj2%projinfo%delete_entry('projname')
        call work_proj2%projinfo%delete_entry('projfile')
        call cline%set('projfile', trim(ORIG_work_projfile))
        call cline%set('projname', trim(get_fbody(trim(ORIG_work_projfile),trim('simple'))))
        call work_proj2%update_projinfo(cline)
        call work_proj2%write(ORIG_work_projfile)
        ! split
        if( l_lpset )then
            if( l_shmem )then
                ! all good
            else
                call work_proj2%split_stk(params%nparts)
            endif
        endif
        ! refinement scaling
        do_autoscale     = .false.
        scale_factor2    = 1.0
        smpd_target      = max(params%smpd, params%smpd/SCALEFAC2_TARGET)
        params%smpd_crop = params%smpd
        params%box_crop  = params%box
        params%msk_crop  = params%msk
        if( params%l_autoscale )then
            call autoscale(params%box, params%smpd, smpd_target, params%box_crop, params%smpd_crop, scale_factor2, minbox=MINBOX)
            do_autoscale = params%box_crop < params%box
            if( do_autoscale )then
                params%msk_crop = params%msk * scale_factor2
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
            endif
        endif
        ! refinement stage
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> REFINEMENT'
        write(logfhandle,'(A)') '>>>'
        iter = iter + 1
        call cline_refine3D_refine%set('startit',  real(iter))
        call cline_refine3D_refine%set('box_crop', real(params%box_crop))
        call cline_refine3D_refine%set('smpd_crop',params%smpd_crop)
        if( l_euclid )then
            ! sigma2 at original sampling
            cline_calc_pspec = cline
            call cline_calc_pspec%set('prg',      'calc_pspec' )
            call cline_calc_pspec%set('projfile', ORIG_work_projfile)
            call cline_calc_pspec%set('sigma_est','global')
            call cline_calc_pspec%set('box',   real(params%box))
            call cline_calc_pspec%set('smpd',  params%smpd)
            if( l_shmem )then
                call cline_calc_pspec%set('which_iter', real(iter))
                call xcalc_pspec_distr%execute_shmem( cline_calc_pspec )
                call cline_refine3D_refine%set('which_iter', real(iter))
            else
                call cline_calc_pspec%set('which_iter', 0.)
                call xcalc_pspec_distr%execute( cline_calc_pspec )
                do ipart = 1,params%nparts
                    call del_file('init_pspec_part'//int2str(ipart)//'.dat')
                    call del_file('sigma2_noise_part'//int2str(ipart)//'.dat')
                enddo
            endif
        endif
        if( l_shmem )then
            call rec(cline_refine3D_refine, l_rnd=.false.)
            call xrefine3D%execute_shmem(cline_refine3D_refine)
        else
            call xrefine3D_distr%execute(cline_refine3D_refine)
        endif
        iter = nint(cline_refine3D_refine%get_rarg('endit'))
        if( l_euclid )then
            call cline_reconstruct3D%set('which_iter',real(iter))
            if( l_shmem )then
                ! done
            else
                if( do_autoscale )then
                    sigma2_fname    = trim(SIGMA2_GROUP_FBODY)//'0.star'
                    sigma2_fname_sc = trim(SIGMA2_GROUP_FBODY)//trim(int2str(iter))//'.star'
                    call fill_sigma2_before_nyq(sigma2_fname, sigma2_fname_sc)
                endif
                do ipart = 1,params%nparts
                    call del_file('init_pspec_part'//int2str(ipart)//'.dat')
                    call del_file('sigma2_noise_part'//int2str(ipart)//'.dat')
                enddo
            endif
        endif
        ! deals with final volume
        if( do_autoscale )then
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
            write(logfhandle,'(A)') '>>>'
            call cline_reconstruct3D%set('box',   real(params%box))
            call cline_reconstruct3D%set('smpd',  params%smpd)
            ! reconstruction
            call cline_reconstruct3D%set('ml_reg', 'no')
            if( l_shmem )then
                call xreconstruct3D%execute_shmem(cline_reconstruct3D)
            else
                call xreconstruct3D_distr%execute(cline_reconstruct3D)
            endif
            vol_iter = trim(VOL_FBODY)//trim(str_state)//ext
            ! because postprocess only updates project file when mkdir=yes
            call work_proj2%read_segment('out', ORIG_work_projfile)
            call work_proj2%add_vol2os_out(vol_iter, params%smpd, 1, 'vol')
            if( .not.l_lpset )then
                call work_proj2%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), 1, params%box)
            endif
            call work_proj2%write_segment_inside('out',ORIG_work_projfile)
            call xpostprocess%execute(cline_postprocess)
        else
            vol_iter = trim(VOL_FBODY)//trim(str_state)//ext
            call vol%new([params%box,params%box,params%box],params%smpd)
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
        call work_proj2%read_segment('ptcl3D', ORIG_work_projfile)
        call work_proj2%os_ptcl3D%delete_entry('stkind')
        call work_proj2%os_ptcl3D%delete_entry('eo')
        params_glob%nptcls = ncavgs
        if( l_lpset )then
            spproj%os_cls3D = work_proj2%os_ptcl3D
        else
            call spproj%os_cls3D%new(ncavgs, is_ptcl=.false.)
            do icls=1,ncavgs
                call spproj%os_cls3D%transfer_ori(icls, work_proj2%os_ptcl3D, icls)
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
        call img%new([params%box,params%box,1], params%smpd)
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
                call xreconstruct3D%execute_shmem(cline)
                call build%spproj_field%kill
                if( l_rnd )then
                    call simple_rename('snhc_recvol_state01_even.mrc', 'startvol_even.mrc')
                    call simple_rename('snhc_recvol_state01_odd.mrc',  'startvol_odd.mrc')
                    call simple_rename('snhc_recvol_state01.mrc',      'startvol.mrc')
                else
                    call simple_rename('recvol_state01_even.mrc', 'startvol_even.mrc')
                    call simple_rename('recvol_state01_odd.mrc',  'startvol_odd.mrc')
                    call simple_rename('recvol_state01.mrc',      'startvol.mrc')
                endif
                call cline%set('vol1', 'startvol.mrc')
            end subroutine rec

            subroutine prep_eo_stks_refine
                type(ori) :: o, o_even, o_odd
                integer   :: even_ind, odd_ind, state, icls
                call cline_refine3D_refine%set('frcs',frcs_fname)
                ! add stks
                call work_proj2%add_stk(stk_even, ctfvars)
                call work_proj2%add_stk(stk_odd,  ctfvars)
                ! update orientations parameters
                do icls=1,ncavgs
                    even_ind = icls
                    odd_ind  = ncavgs+icls
                    call work_proj1%os_ptcl3D%get_ori(icls, o)
                    state    = o%get_state()
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
                call work_proj2%os_ptcl3D%delete_entry('lp')
                call o%kill
                call o_even%kill
                call o_odd%kill
            end subroutine prep_eo_stks_refine

            subroutine conv_eo( os )
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

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio_3Dmodel( self, cline )
        use simple_convergence, only: convergence
        use simple_class_frcs,  only: class_frcs
        class(abinitio_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        real,    parameter :: SCALEFAC      = 0.667
        real,    parameter :: CENLP_DEFAULT = 30.
        real,    parameter :: LP_DEFAULT    = 6.
        real,    parameter :: LPSTART_DEFAULT = 30., LPSTOP_DEFAULT=LP_DEFAULT
        integer, parameter :: MINBOX  = 88
        integer, parameter :: NSTAGES = 5
        integer, parameter :: MAXITS1=100, MAXITS2=30, MAXITS_SHORT1=15, MAXITS_SHORT2=25
        integer, parameter :: MAXITS_GLOB_LP_SET = MAXITS1 + MAXITS2 * 2
        integer, parameter :: MAXITS_GLOB = MAXITS_SHORT1 * (NSTAGES - 2) + MAXITS_SHORT2 * 2 + MAXITS2
        integer, parameter :: NSPACE1=500, NSPACE2=1000, NSPACE3=2000
        integer, parameter :: SYMSEARCH_ITER = 3
        ! commanders
        type(refine3D_commander_distr)      :: xrefine3D_distr
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        type(postprocess_commander)         :: xpostprocess
        type(symaxis_search_commander)      :: xsymsrch
        type(reconstruct3D_commander)       :: xreconstruct3D
        ! command lines
        type(cmdline)                 :: cline_refine3D, cline_reconstruct3D
        type(cmdline)                 :: cline_postprocess, cline_symsrch
        ! other
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(convergence)             :: conv
        type(sym)                     :: se1, se2
        type(class_frcs)              :: clsfrcs 
        type(image)                   :: final_vol, reprojs
        logical,          allocatable :: mask(:)
        integer,          allocatable :: pinds(:)
        character(len=:), allocatable :: vol_type, str_state, vol, vol_pproc, vol_pproc_mirr, frcs_fname
        real    :: smpd_target, lp_target, scale, trslim, cenlp, symlp, dummy, auto_lpstart, auto_lpstop
        integer :: iter, it, prev_box_crop, maxits, noris, nsamples
        logical :: l_autoscale, l_lpset, l_err, l_srch4symaxis, l_symran, l_sym, l_lpstop_set
        logical :: l_lpstart_set
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',        'yes')
        if( .not. cline%defined('refine')     ) call cline%set('refine',      'prob')
        if( .not. cline%defined('autoscale')  ) call cline%set('autoscale',    'yes')
        if( .not. cline%defined('ml_reg')     ) call cline%set('ml_reg',        'no')
        if( .not. cline%defined('sigma_est')  ) call cline%set('sigma_est', 'global')
        if( .not. cline%defined('prob_sh')    ) call cline%set('prob_sh',       'no')
        if( .not. cline%defined('prob_athres')) call cline%set('prob_athres',    10.)
        if( .not. cline%defined('center')     ) call cline%set('center',        'no')
        if( .not. cline%defined('objfun')     ) call cline%set('objfun',    'euclid')
        if( .not. cline%defined('oritype')    ) call cline%set('oritype',   'ptcl3D')
        if( .not. cline%defined('pgrp')       ) call cline%set('pgrp',          'c1')
        if( .not. cline%defined('pgrp_start') ) call cline%set('pgrp_start',    'c1')
        ! resolution limit strategy
        l_lpset       = .false.
        l_lpstop_set  = cline%defined('lpstop')
        l_lpstart_set = cline%defined('lpstart')
        if( cline%defined('lp') )then
            if( l_lpstart_set .or. l_lpstop_set )then
                THROW_HARD('One of LP or LPSTART & LPSTOP must be defined!')
            endif
            l_lpset = .true.
        else
            if( .not.l_lpstart_set ) call cline%set('lpstart',LPSTART_DEFAULT)
            if( .not.l_lpstop_set  ) call cline%set('lpstop', LPSTOP_DEFAULT)
        endif
        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call cline%delete('autoscale')
        call cline%delete('lpstart')
        call cline%delete('lpstop')
        call cline%delete('lp')
        str_state = int2str_pad(1,2)
        if( l_lpset )then
            params%lpstop  = params%lp
            params%lpstart = params%lp
        endif
        ! read project
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
        if( .not. cline%defined('vol1') )then
            ! randomize projection directions
            select case(trim(params%oritype))
                case('cls3D')
                    if( spproj%os_cls3D%get_noris() < 1 )then
                        THROW_HARD('Class averages could not be found in the project')
                    endif
                    vol_type = 'vol_cavg'
                    call spproj%os_cls3D%rnd_oris
                    call spproj%os_cls3D%set_all2single('stkind', 1.)
                case('ptcl3D')
                    if( spproj%os_ptcl3D%get_noris() < 1 )then
                        THROW_HARD('Particles could not be found in the project')
                    endif
                    vol_type = 'vol'
                    call spproj%os_ptcl3D%rnd_oris
                case DEFAULT
                    THROW_HARD('Unsupported ORITYPE; exec_abinitio_3Dmodel')
            end select
            call spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        ! ! Automated resolution limits, not active yet
        ! call mskdiam2lplimits(params%mskdiam, auto_lpstart, auto_lpstop, cenlp)
        ! if( .not. cline%defined('cenlp') )then
        !     params%cenlp = cenlp
        !     call cline%set('cenlp', params%cenlp)
        ! endif
        ! if( .not.l_lpset )then
        !     if( l_lpstart_set .and. l_lpstop_set )then
        !         ! all set
        !     else
        !         if( .not.l_lpstart_set ) params%lpstart = 15.
        !         if( .not.l_lpstop_set  )then
        !             call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
        !             if( file_exists(frcs_fname) )then
        !                 call clsfrcs%read(frcs_fname)
        !                 auto_lpstop = clsfrcs%estimate_lp_for_align()
        !                 auto_lpstop = min(max(auto_lpstop,6.),8.)
        !                 call clsfrcs%kill
        !             else
        !                 auto_lpstop = 8.
        !             endif
        !             params%lpstop  = auto_lpstop
        !         endif
        !     endif    
        ! endif
        ! symlp = params%lpstart
        write(logfhandle,'(A,F5.1)') '>>> STARTING RESOLUTION LIMIT (IN A): ', params%lpstart
        write(logfhandle,'(A,F5.1)') '>>> HARD     RESOLUTION LIMIT (IN A): ', params%lpstop
        if( trim(params%center).eq.'yes' )then
            write(logfhandle,'(A,F5.1)') '>>> CENTERING  LOW-PASS LIMIT (IN A): ', params%cenlp
        endif
        ! centering & symmetry resolution limit
        call mskdiam2lplimits(params%mskdiam, symlp, dummy, cenlp)
        if( .not. cline%defined('cenlp') )then
            params%cenlp = cenlp
            call cline%set('cenlp', params%cenlp)
        endif
        ! symmetry
        if( l_lpset )then
            ! from mskdiam2lplimits lpstart above
        else
            symlp = (params%lpstart+params%lpstop)/2.
        endif
        l_srch4symaxis = trim(params%pgrp) .ne. trim(params%pgrp_start)
        l_symran       = .false.
        l_sym          = l_srch4symaxis
        if( params%pgrp_start.ne.'c1' .or. params%pgrp.ne.'c1' )then
            se1 = sym(params%pgrp_start)
            se2 = sym(params%pgrp)
            if(se1%get_nsym() > se2%get_nsym())then
                ! ensure se2 is a subgroup of se1
                if( .not. se1%has_subgrp(params%pgrp) )THROW_HARD('Incompatible symmetry groups; exec_abinitio_3Dmodel')
                ! set flag for symmetry randomisation
                ! in case we are moving from a higher to lower group
                l_symran = .true.
            else if( se2%get_nsym() > se1%get_nsym() )then
                ! ensure se1 is a subgroup of se2
                if( .not. se2%has_subgrp(params%pgrp_start) )THROW_HARD('Incompatible symmetry groups; exec_abinitio_3Dmodel')
            endif
        endif
        ! dimensions defaults
        params%box       = spproj%get_box()
        params%smpd_crop = params%smpd
        params%box_crop  = params%box
        l_autoscale      = .false.
        call spproj%kill ! not needed anymore
        ! command-lines
        cline_refine3D      = cline
        cline_reconstruct3D = cline
        cline_postprocess   = cline
        cline_symsrch       = cline
        call cline_refine3D%set('prg',                'refine3D')
        call cline_refine3D%set('projfile',      params%projfile)
        call cline_refine3D%set('center',                   'no')
        call cline_refine3D%set('pgrp',        params%pgrp_start)
        call cline_reconstruct3D%set('prg',      'reconstruct3D')
        call cline_reconstruct3D%set('box',     real(params%box))
        call cline_reconstruct3D%set('projfile', params%projfile)
        call cline_reconstruct3D%set('ml_reg',              'no')
        call cline_reconstruct3D%set('needs_sigma',         'no')
        call cline_reconstruct3D%set('objfun',              'cc')
        call cline_reconstruct3D%set('pgrp',         params%pgrp)
        call cline_postprocess%set('prg',          'postprocess')
        call cline_postprocess%set('projfile',   params%projfile)
        call cline_postprocess%set('imgkind',           vol_type)
        if( l_srch4symaxis )then
            call cline_symsrch%set('prg',     'symaxis_search') ! needed for cluster exec
            call cline_symsrch%set('pgrp',     params%pgrp)
            call cline_symsrch%set('projfile', params%projfile)
            call cline_symsrch%set('hp',       params%hp)
            call cline_symsrch%set('center',   'yes')
        endif
        ! Frequency marching
        iter = 0
        do it = 1,NSTAGES
            write(logfhandle,'(A)')'>>>'
            prev_box_crop = params%box_crop
            ! resolution limit
            if( it == 1 )then
                params%lp = params%lpstart
            else
                params%lp = max(params%lpstop, params%lpstop+(params%lp-params%lpstop)/2.)
            endif
            write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ',it,' WITH LP =',params%lp
            ! dimensions
            prev_box_crop = params%box_crop
            if( params%l_autoscale )then
                l_autoscale = .false.
                lp_target   = params%lp * SCALEFAC
                smpd_target = max(params%smpd, lp_target/2.)
                call autoscale(params%box, params%smpd, smpd_target, params%box_crop, params%smpd_crop, scale, minbox=MINBOX)
                trslim      = max(2.0, AHELIX_WIDTH / params%smpd_crop / 2.0)
                l_autoscale = params%box_crop < params%box
            else
                trslim      = max(2.0, AHELIX_WIDTH / params%smpd / 2.0)
            endif
            if( l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
            endif
            if( it > 1 )then
                call cline_refine3D%set('center', params%center)
                if( prev_box_crop == params%box_crop )then
                    call cline_refine3D%set('continue',  'yes')
                else
                    ! allows reconstruction to correct dimension
                    call cline_refine3D%delete('continue')
                    call cline_refine3D%delete('vol1')
                endif
            endif
            if( it == SYMSEARCH_ITER+1 .and. (l_srch4symaxis .or. l_symran) )then
                ! allows reconstruction to desired point-group
                call cline_refine3D%delete('continue')
                call cline_refine3D%delete('vol1')
                call cline_refine3D%set('pgrp', params%pgrp)
                l_srch4symaxis = .false.
                l_symran       = .false.
            endif
            ! stage updates
            call cline_refine3D%set('box_crop', params%box_crop)
            call cline_refine3D%set('lp',       params%lp)
            call cline_refine3D%set('startit',  iter+1)
            ! # of iterations
            maxits = MAXITS_SHORT1
            if( it >= NSTAGES-1 ) maxits = MAXITS_SHORT2
            call cline_refine3D%set('maxits',      maxits)
            call cline_refine3D%set('maxits_glob', MAXITS_GLOB)
            call cline_refine3D%set('lp_iters',    maxits)
            ! projection directions & shift
            if( it < NSTAGES )then
                call cline_refine3D%set('nspace', NSPACE1)
                call cline_refine3D%set('trs',      0.)
            else
                call cline_refine3D%set('nspace', NSPACE2)
                call cline_refine3D%set('trs',    trslim)
            end if
            call cline_refine3D%set('ml_reg',    'no')       ! since ml_reg seems to interfer with finding a good lowres shape            
            call cline_refine3D%set('it_history', min(3,it)) ! particles sampled in most recent past iterations also included in rec
            ! execution
            call exec_refine3D(iter)
            ! symmetrization
            if( it == SYMSEARCH_ITER ) call symmetrize
        enddo
        ! Final stage
        write(logfhandle,'(A)')'>>>'
        write(logfhandle,'(A,F6.2)')'>>> FINAL STAGE WITH LP=',params%lpstop
        prev_box_crop = params%box_crop
        if( params%l_autoscale )then
            l_autoscale = .false.
            lp_target   = params%lpstop * SCALEFAC
            smpd_target = max(params%smpd, lp_target/2.)
            call autoscale(params%box, params%smpd, smpd_target, params%box_crop, params%smpd_crop, scale, minbox=MINBOX)
            trslim      = max(MINSHIFT, AHELIX_WIDTH / params%smpd_crop)
            l_autoscale = params%box_crop < params%box
        else
            trslim      = max(MINSHIFT, AHELIX_WIDTH / params%smpd)
        endif
        if( l_autoscale )then
            write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
        endif
        ! Dealing with reconstruction
        if( params%l_ml_reg )then
            ! need to reconstruct
            call cline_refine3D%set('ml_reg', 'yes')
            call cline_refine3D%delete('continue')
            call cline_refine3D%delete('vol1')
        else
            call cline_refine3D%set('ml_reg', 'no')
            if( prev_box_crop == params%box_crop )then
                call cline_refine3D%set('continue',  'yes')
            else
                ! allows reconstruction to correct dimension
                call cline_refine3D%delete('continue')
                call cline_refine3D%delete('vol1')
            endif
        endif
        ! Final step
        call cline_refine3D%set('box_crop',   params%box_crop)
        call cline_refine3D%set('lp',         params%lpstop)
        call cline_refine3D%set('trs',        trslim)
        call cline_refine3D%set('startit',    iter+1)
        call cline_refine3D%set('maxits',     MAXITS2)
        call cline_refine3D%set('lp_iters',   MAXITS2)
        call cline_refine3D%set('nspace',     NSPACE3)
        call cline_refine3D%set('it_history',       3) ! particles sampled in three most recent past iterations also included in rec
        call cline_refine3D%set('mov_avg_vol',  'yes')
        ! execution
        call exec_refine3D(iter)
        ! for visualization
        call final_vol%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop)
        call final_vol%read(trim(VOL_FBODY)//trim(str_state)//trim(params%ext))
        call final_vol%generate_orthogonal_reprojs(reprojs)
        call reprojs%write_jpg('orthogonal_reprojs.jpg')
        call final_vol%kill
        call reprojs%kill
        ! Final reconstruction at original scale
        vol = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
        ! no ML-filtering
        call cline_reconstruct3D%set('ml_reg',      'no')
        call cline_reconstruct3D%set('needs_sigma', 'no')
        call cline_reconstruct3D%set('objfun',      'cc')
        ! no fractional or stochastic updates
        call cline_reconstruct3D%delete('update_frac')
        call cline_reconstruct3D%delete('stoch_update')
        call cline_reconstruct3D%delete('it_history')
        ! reconstruction
        if( l_autoscale )then
            write(logfhandle,'(A)') '>>>'
            write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
            write(logfhandle,'(A)') '>>>'
            params%box_crop = params%box
            call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
            call spproj%read_segment('out',params%projfile)
            call spproj%add_vol2os_out(vol, params%smpd, 1, vol_type)
            if( trim(params%oritype).eq.'ptcl3D' )then
                call spproj%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), 1, params%box)
            endif
            call spproj%write_segment_inside('out',params%projfile)
            ! post-processing
            call cline_postprocess%set('lp', params%lpstop)
            call xpostprocess%execute(cline_postprocess)
        else
            call final_vol%new([params%box,params%box,params%box],params%smpd)
            call final_vol%read(vol)
            call final_vol%mirror('x')
            call final_vol%write(add2fbody(vol,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX)))
            call final_vol%kill
        endif
        vol_pproc      = add2fbody(vol,params%ext,PPROC_SUFFIX)
        vol_pproc_mirr = add2fbody(vol,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX))
        if( file_exists(vol)            ) call simple_rename(vol,            trim(REC_FBODY)//trim(params%ext))
        if( file_exists(vol_pproc)      ) call simple_rename(vol_pproc,      trim(REC_PPROC_FBODY)//trim(params%ext))
        if( file_exists(vol_pproc_mirr) ) call simple_rename(vol_pproc_mirr, trim(REC_PPROC_MIRR_FBODY)//trim(params%ext))
        ! transfer cls3D parameters to particles
        if( trim(params%oritype) .eq. 'cls3D' )then
            call spproj%read_segment('cls3D', params%projfile)
            call spproj%read_segment('ptcl3D',params%projfile)
            call spproj%map2ptcls
            call spproj%write_segment_inside('ptcl3D',params%projfile)
            call spproj%os_cls3D%delete_entry('stkind')
            call spproj%write_segment_inside('cls3D',params%projfile)
        endif
        ! cleanup
        call spproj%kill
        call qsys_cleanup
        call simple_end('**** SIMPLE_ABINITIO_3DMODEL NORMAL STOP ****')
        contains

            subroutine exec_refine3D( iter )
                integer,          intent(out) :: iter
                character(len=:), allocatable :: stage
                call cline_refine3D%delete('endit')
                call xrefine3D_distr%execute_shmem(cline_refine3D)
                call conv%read(l_err)
                iter = nint(conv%get('iter'))
                call del_files(DIST_FBODY,      params_glob%nparts,ext='.dat')
                call del_files(ASSIGNMENT_FBODY,params_glob%nparts,ext='.dat')
                call del_file(trim(DIST_FBODY)      //'.dat')
                call del_file(trim(ASSIGNMENT_FBODY)//'.dat')
                if( it <= NSTAGES )then
                    stage = '_stage_'//int2str(it)
                    vol   = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
                    vol_pproc = add2fbody(vol,params%ext,PPROC_SUFFIX)
                    if( file_exists(vol)      ) call simple_copy_file(vol,       add2fbody(vol,params%ext,stage))
                    if( file_exists(vol_pproc)) call simple_copy_file(vol_pproc, add2fbody(vol_pproc,params%ext,stage))
                endif
            end subroutine exec_refine3D

            subroutine symmetrize()
                if( l_symran )then
                    call spproj%read_segment(params%oritype, params%projfile)
                    select case(trim(params%oritype))
                    case('cls3D')
                        call se1%symrandomize(spproj%os_cls3D)
                    case('ptcl3D')
                        call se1%symrandomize(spproj%os_ptcl3D)
                    end select
                    call spproj%write_segment_inside(params%oritype, params%projfile)
                endif
                if( l_srch4symaxis )then
                    write(logfhandle,'(A)') '>>>'
                    write(logfhandle,'(A)') '>>> SYMMETRY AXIS SEARCH'
                    write(logfhandle,'(A)') '>>>'
                    symlp = max(symlp, params%lp)
                    call cline_symsrch%set('lp',       symlp)
                    call cline_symsrch%set('box_crop', params%box_crop)
                    call cline_symsrch%set('vol1',     trim(VOL_FBODY)//trim(str_state)//trim(params%ext))
                    call xsymsrch%execute_shmem(cline_symsrch)
                    call del_file('SYMAXIS_SEARCH_FINISHED')
                endif
            end subroutine symmetrize

    end subroutine exec_abinitio_3Dmodel

end module simple_commander_abinitio
