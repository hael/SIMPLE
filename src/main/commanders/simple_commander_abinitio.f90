! concrete commander: high-level workflows
module simple_commander_abinitio
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters, params_glob
use simple_sp_project,         only: sp_project
use simple_stack_io,           only: stack_io
use simple_qsys_env,           only: qsys_env
use simple_commander_base,     only: commander_base
use simple_commander_volops,   only: reproject_commander, symaxis_search_commander, postprocess_commander, symmetrize_map_commander
use simple_commander_rec,      only: reconstruct3D_commander, reconstruct3D_commander_distr
use simple_commander_refine3D, only: refine3D_commander, refine3D_commander_distr
use simple_procimgstk,         only: shift_imgfile
use simple_image,              only: image
use simple_builder,            only: builder
use simple_class_frcs,         only: class_frcs
use simple_commander_euclid
use simple_euclid_sigma2
use simple_qsys_funs
implicit none

public :: initial_3Dmodel_commander, abinitio_3Dmodel_commander
public :: abinitio_3Dmodel2_commander
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

type, extends(commander_base) :: abinitio_3Dmodel2_commander
    contains
    procedure :: execute => exec_abinitio_3Dmodel2
end type abinitio_3Dmodel2_commander

character(len=STDLEN), parameter :: REC_FBODY            = 'rec_final'
character(len=STDLEN), parameter :: REC_PPROC_FBODY      = trim(REC_FBODY)//trim(PPROC_SUFFIX)
character(len=STDLEN), parameter :: REC_PPROC_MIRR_FBODY = trim(REC_PPROC_FBODY)//trim(MIRR_SUFFIX)
real,                  parameter :: LPSTART_LB=10., LPSTART_DEFAULT=20., LPSTOP_LB=6.
real,                  parameter :: CENLP_DEFAULT=30.
real,                  parameter :: LPSYMSRCH_LB=12.
logical                          :: l_srch4symaxis=.false., l_symran=.false., l_sym=.false.
type(sym)                        :: se1,se2

contains


    subroutine set_symmetry_class_vars( params )
        class(parameters), intent(in) :: params
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
    end subroutine set_symmetry_class_vars

    function calc_lplim_final_stage( spproj, nbest ) result( lplim )
        class(sp_project), intent(in) :: spproj
        integer,           intent(in) :: nbest
        real,    allocatable :: res(:), tmp_rarr(:)
        integer, allocatable :: states(:), tmp_iarr(:)
        real :: lplim
        tmp_rarr  = spproj%os_cls2D%get_all('res')
        tmp_iarr  = nint(spproj%os_cls2D%get_all('state'))
        res       = pack(tmp_rarr, mask=(tmp_iarr>0))
        call hpsort(res)
        lplim = median_nocopy(res(:nbest))
        deallocate(tmp_rarr, tmp_iarr, res)
    end function calc_lplim_final_stage

    !> for generation of an initial 3d model from class averages
    subroutine exec_initial_3Dmodel( self, cline )
        class(initial_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! constants
        integer,               parameter :: NSTAGES         = 8
        integer,               parameter :: PHASES(3)       = [2,6,8]
        integer,               parameter :: MAXITS(3)       = [20,15,10]
        integer,               parameter :: MAXITS_GLOB(3)  = [2*20,4*15,2*10]
        integer,               parameter :: NSPACE(3)       = [500,1000,2500]
        integer,               parameter :: SYMSRCH_STAGE   = 3
        character(len=STDLEN), parameter :: work_projfile   = 'initial_3Dmodel_tmpproj.simple'
        ! distributed commanders
        type(calc_pspec_commander_distr) :: xcalc_pspec_distr
        ! shared-mem commanders
        type(refine3D_commander)         :: xrefine3D
        type(reconstruct3D_commander)    :: xreconstruct3D
        type(symmetrize_map_commander)   :: xsymmap
        type(reproject_commander)        :: xreproject
        type(postprocess_commander)      :: xpostprocess
        ! command lines
        type(cmdline) :: cline_refine3D
        type(cmdline) :: cline_symmap
        type(cmdline) :: cline_reconstruct3D, cline_postprocess
        type(cmdline) :: cline_reproject, cline_calc_pspec
        ! other
        character(len=:),  allocatable :: stk, stkpath, orig_stk, frcs_fname, shifted_stk, stk_even, stk_odd, ext, vol_str, vol_name
        integer,           allocatable :: states(:)
        real,              allocatable :: frcs_avg(:)
        type(lp_crop_inf), allocatable :: lpinfo(:)
        character(len=2)      :: str_state
        type(ori)             :: o, o_even, o_odd
        type(parameters)      :: params
        type(ctfparams)       :: ctfvars
        type(sp_project)      :: spproj, work_proj
        type(image)           :: img, vol
        type(stack_io)        :: stkio_r, stkio_r2, stkio_w
        type(class_frcs)      :: clsfrcs
        character(len=STDLEN) :: vol_iter, vol_iter_pproc, vol_iter_pproc_mirr, frckind, vol_sym
        real                  :: lpsym, lpfinal
        integer               :: icls, ncavgs, cnt, iter, ipart, even_ind, odd_ind, state, filtsz, istage
        call cline%set('objfun',    'euclid') ! use noise normalized Euclidean distances from the start
        call cline%set('sigma_est', 'global') ! obviously
        call cline%set('oritype',      'out') ! because cavgs are part of out segment
        call cline%set('bfac',            0.) ! because initial models should not be sharpened
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('overlap')     ) call cline%set('overlap',     0.95) ! needed to prevent premature convergence
        if( .not. cline%defined('prob_athres') ) call cline%set('prob_athres',  90.) ! reduces # failed runs on trpv1 from 4->2/10
        if( .not. cline%defined('cenlp')       ) call cline%set('cenlp', CENLP_DEFAULT)
        if( .not. cline%defined('imgkind')     ) call cline%set('imgkind',   'cavg') ! whether to use classes generated from 2D/3D
        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir',       'no')   ! to avoid nested directory structure
        call cline%set('oritype', 'ptcl3D')   ! from now on we are in the ptcl3D segment, final report is in the cls3D segment
        ! state string
        str_state = int2str_pad(1,2)
        ! symmetry
        call set_symmetry_class_vars(params)
        ! read project
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
        ! whether to use classes generated from 2D or 3D
        select case(trim(params%imgkind))
            case('cavg')
                frckind = 'frc2D'
                states  = nint(spproj%os_cls2D%get_all('state'))
            case('cavg3D')
                frckind = 'frc3D'
                states  = nint(spproj%os_cls3D%get_all('state'))
            case DEFAULT
                THROW_HARD('Unsupported IMGKIND!')
        end select
        call cline%delete('imgkind') ! no interference down the line
        ! retrieve cavgs stack info
        call spproj%get_cavgs_stk(stk, ncavgs, params%smpd, imgkind=params%imgkind, stkpath=stkpath)
        if(.not. file_exists(stk)) stk = trim(stkpath) // '/' // trim(stk)
        if(.not. file_exists(stk)) THROW_HARD('cavgs stk does not exist; simple_commander_abinitio')
        orig_stk        = stk
        ext             = '.'//fname2ext(stk)
        stk_even        = add2fbody(trim(stk), trim(ext), '_even')
        stk_odd         = add2fbody(trim(stk), trim(ext), '_odd')
        if( .not. file_exists(stk_even) ) THROW_HARD('Even cavgs stk: '//trim(stk_even)//' does not exist!')
        if( .not. file_exists(stk_odd)  ) THROW_HARD('Odd cavgs stk: '//trim(stk_odd)//' does not exist!')
        ctfvars%ctfflag = CTFFLAG_NO
        ctfvars%smpd    = params%smpd
        shifted_stk     = basename(add2fbody(stk, ext, '_shifted'))
        if( count(states==0) .eq. ncavgs )then
            THROW_HARD('no class averages detected in project file: '//trim(params%projfile)//'; initial_3Dmodel')
        endif
        ! retrieve FRC info
        call spproj%get_frcs(frcs_fname, frckind, fail=.false.)
        if( .not.file_exists(frcs_fname) )then
            ! 08/24 This is a backwards compatibility patch to account for error in metadata
            ! on exit of streaming related to GUI directory structure (now fixed and cf above get_cavgs_stk).
            ! Will need to harmonize (move to absolute path?).
            frcs_fname = trim(stkpath)//'/'//trim(frcs_fname)
            if( .not.file_exists(frcs_fname) )then
                THROW_HARD('the project file does not contain an FRCs file, which is required')
            endif
        endif
        ! work out low-pass limits and downscaling parameters
        params%frcs = trim(frcs_fname)
        call clsfrcs%read(frcs_fname)
        filtsz = clsfrcs%get_filtsz()
        allocate(frcs_avg(filtsz), source=0.)
        call clsfrcs%avg_frc_getter(frcs_avg, states)
        allocate(lpinfo(NSTAGES))
        lpfinal = max(LPSTOP_LB,calc_lplim_final_stage(spproj,3))
        call lpstages(params%box, NSTAGES, frcs_avg, params%smpd, LPSTART_LB, LPSTART_DEFAULT, lpfinal, lpinfo, verbose=.true.)
        ! prepare a temporary project file
        work_proj%projinfo = spproj%projinfo
        work_proj%compenv  = spproj%compenv
        if( spproj%jobproc%get_noris()  > 0 ) work_proj%jobproc = spproj%jobproc
        ! name change
        call work_proj%projinfo%delete_entry('projname')
        call work_proj%projinfo%delete_entry('projfile')
        call cline%set('projfile', trim(work_projfile))
        call cline%set('projname', trim(get_fbody(trim(work_projfile),trim('simple'))))
        call work_proj%update_projinfo(cline)
        ! add stks to temporary project
        call work_proj%add_stk(stk_even, ctfvars)
        call work_proj%add_stk(stk_odd,  ctfvars)
        ! update orientations parameters
        do icls=1,ncavgs
            even_ind = icls
            odd_ind  = ncavgs + icls
            call work_proj%os_ptcl3D%get_ori(icls, o)
            call o%set('class', real(icls))
            call o%set('state', real(states(icls)))
            ! even
            o_even = o
            call o_even%set('eo', 0.)
            call o_even%set('stkind', work_proj%os_ptcl3D%get(even_ind,'stkind'))
            call work_proj%os_ptcl3D%set_ori(even_ind, o_even)
            ! odd
            o_odd = o
            call o_odd%set('eo', 1.)
            call o_odd%set('stkind', work_proj%os_ptcl3D%get(odd_ind,'stkind'))
            call work_proj%os_ptcl3D%set_ori(odd_ind, o_odd)
        enddo
        params_glob%nptcls = work_proj%get_nptcls()
        call work_proj%write()
        ! prepare command lines from prototype
        cline_reconstruct3D = cline
        cline_refine3D      = cline
        cline_reproject     = cline
        cline_symmap        = cline
        ! map symmetrization
        if( l_srch4symaxis )then
            call cline_symmap%set('prg',     'symaxis_search')
            call cline_symmap%set('pgrp',   trim(params%pgrp))
            call cline_symmap%set('projfile', params%projfile)
            if( .not. cline_symmap%defined('cenlp') ) call cline_symmap%set('cenlp', CENLP_DEFAULT)
            call cline_symmap%set('hp',             params%hp)
        endif
        ! re-reconstruct & re-project volume
        call cline_reconstruct3D%set('prg',       'reconstruct3D')
        call cline_reconstruct3D%set('box',      real(params%box))
        call cline_reconstruct3D%set('projfile',    work_projfile)
        call cline_reconstruct3D%set('needs_sigma',         'yes')
        call cline_postprocess%set('prg',           'postprocess')
        call cline_postprocess%set('projfile',      work_projfile)
        call cline_postprocess%set('mkdir',                  'no')
        call cline_postprocess%delete('bfac') ! sharpen final map
        call cline_reproject%set('prg',               'reproject')
        call cline_reproject%set('pgrp',        trim(params%pgrp))
        call cline_reproject%set('outstk',         'reprojs'//ext)
        call cline_reproject%set('smpd',              params%smpd)
        call cline_reproject%set('box',          real(params%box))
        ! Frequency marching
        do istage = 1, NSTAGES
            write(logfhandle,'(A)')'>>>'
            write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ', istage,' WITH LP =', lpinfo(istage)%lp
            call set_cline_refine3D(istage)
            if( istage == 1 ) call rndstart(cline_refine3D)
            if( lpinfo(istage)%l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',lpinfo(istage)%box_crop
            endif
            call xrefine3D%execute_shmem(cline_refine3D)
            call work_proj%read_segment('ptcl3D', work_projfile)

            ! Symmetrization
            if( istage == SYMSRCH_STAGE )then
                if( l_symran )then
                    call se1%symrandomize(spproj%os_ptcl3D)
                    call spproj%write_segment_inside('ptcl3D', work_projfile)
                endif
                if( l_srch4symaxis )then
                    vol_iter = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
                    if( .not. file_exists(vol_iter) ) THROW_HARD('input volume to map symmetrization does not exist')
                    call cline_symmap%set('vol1',               trim(vol_iter))
                    call cline_symmap%set('smpd',     lpinfo(istage)%smpd_crop)
                    call cline_symmap%set('box', real(lpinfo(istage)%box_crop))
                    vol_sym  = 'symmetrized_map'//trim(params%ext)
                    call cline_symmap%set('outvol', trim(vol_sym))
                    lpsym = max(LPSYMSRCH_LB,lpinfo(SYMSRCH_STAGE)%lp)
                    call cline_symmap%set('lp', lpsym)
                    write(logfhandle,'(A,F5.1)') '>>> DID SET MAP SYMMETRIZATION LOW-PASS LIMIT (IN A) TO: ', lpsym
                    write(logfhandle,'(A)') '>>>'
                    write(logfhandle,'(A)') '>>> MAP SYMMETRIZATION'
                    write(logfhandle,'(A)') '>>>'
                    call xsymmap%execute_shmem(cline_symmap)
                    call del_file('SYMAXIS_SEARCH_FINISHED')
                    call simple_copy_file(vol_sym, vol_iter)
                endif
            endif
        end do
        ! sigma2 at original sampling
        cline_calc_pspec = cline
        call cline_calc_pspec%set('prg',      'calc_pspec' )
        call cline_calc_pspec%set('projfile', work_projfile)
        call cline_calc_pspec%set('box',   real(params%box))
        call cline_calc_pspec%set('smpd',       params%smpd)
        call cline_calc_pspec%set('which_iter',  real(iter))
        call xcalc_pspec_distr%execute_shmem(cline_calc_pspec)
        iter = nint(cline_refine3D%get_rarg('endit'))
        call cline_reconstruct3D%set('which_iter',real(iter))
        ! deal with final volume
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
        write(logfhandle,'(A)') '>>>'
        call cline_reconstruct3D%set('box',  real(params%box))
        call cline_reconstruct3D%set('smpd', params%smpd)
        ! reconstruction
        call xreconstruct3D%execute_shmem(cline_reconstruct3D)
        vol_iter = trim(VOL_FBODY)//trim(str_state)//ext
        ! because postprocess only updates project file when mkdir=yes
        call work_proj%read_segment('out', work_projfile)
        call work_proj%add_vol2os_out(vol_iter, params%smpd, 1, 'vol')
        call work_proj%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), 1, params%box)
        call work_proj%write_segment_inside('out', work_projfile)
        call xpostprocess%execute(cline_postprocess)
        vol_iter_pproc      = add2fbody(vol_iter,ext,PPROC_SUFFIX)
        vol_iter_pproc_mirr = add2fbody(vol_iter,ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX))
        if( file_exists(vol_iter)            ) call simple_rename(vol_iter,            trim(REC_FBODY)//ext)
        if( file_exists(vol_iter_pproc)      ) call simple_rename(vol_iter_pproc,      trim(REC_PPROC_FBODY)//ext)
        if( file_exists(vol_iter_pproc_mirr) ) call simple_rename(vol_iter_pproc_mirr, trim(REC_PPROC_MIRR_FBODY)//ext)
        ! updates original cls3D segment
        call work_proj%read_segment('ptcl3D', work_projfile)
        call work_proj%os_ptcl3D%delete_entry('stkind')
        call work_proj%os_ptcl3D%delete_entry('eo')
        params_glob%nptcls = ncavgs
        call spproj%os_cls3D%new(ncavgs, is_ptcl=.false.)
        do icls=1,ncavgs
            call spproj%os_cls3D%transfer_ori(icls, work_proj%os_ptcl3D, icls)
        enddo
        call conv_eo(work_proj%os_ptcl3D)
        ! revert splitting
        call spproj%os_cls3D%set_all2single('stkind', 1.)
        ! map the orientation parameters obtained for the clusters back to the particles
        select case(trim(params%imgkind))
        case('cavg')
            call spproj%map2ptcls
        case('cavg3D')
            ! For internal testing, particle parameters are not updated
        end select
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
        call o%kill
        call o_even%kill
        call o_odd%kill
        call clsfrcs%kill
        call work_proj%kill
        call del_file(work_projfile)
        call simple_rmdir(STKPARTSDIR)
        call simple_end('**** SIMPLE_INITIAL_3DMODEL NORMAL STOP ****')

        contains

            subroutine rndstart( cline )
                class(cmdline), intent(inout) :: cline
                call work_proj%os_ptcl3D%rnd_oris
                call work_proj%os_ptcl3D%zero_shifts
                call work_proj%write_segment_inside('ptcl3D', work_projfile)
                call cline%set('mkdir', 'no') ! to avoid nested dirs
                call cline%set('objfun', 'cc')
                call cline%set('silence_fsc', 'yes')
                call xreconstruct3D%execute_shmem(cline)
                call cline%set('objfun', trim(params%objfun))
                call simple_copy_file('recvol_state01_even.mrc', 'startvol_even_unfil.mrc')
                call simple_copy_file('recvol_state01_odd.mrc',  'startvol_odd_unfil.mrc')
                call simple_rename(   'recvol_state01_even.mrc', 'startvol_even.mrc')
                call simple_rename(   'recvol_state01_odd.mrc',  'startvol_odd.mrc')
                call simple_rename(   'recvol_state01.mrc',      'startvol.mrc')
                call cline%set('vol1', 'startvol.mrc')
            end subroutine rndstart

            subroutine set_cline_refine3D( istage )
                integer, intent(in) :: istage
                character(len=:), allocatable :: silence_fsc, sh_first, prob_sh, ml_reg, refine, icm
                integer :: iphase, s
                real    :: trs, rnspace, rmaxits, rmaxits_glob, riter, snr_noise_reg
                if( istage > 1 )then ! use the previous volume rather than the random starting volume
                    do s = 1, params%nstates
                        vol_str   = 'vol'//trim(int2str(s))
                        call cline_refine3D%delete(vol_str)
                        str_state = int2str_pad(s,2)
                        vol_name  = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
                        call cline_refine3D%set(vol_str, vol_name)
                    enddo
                endif
                ! iteration number bookkeeping
                if( cline_refine3D%defined('endit') )then
                    riter = cline_refine3D%get_rarg('endit')
                else
                    riter = 0.
                endif
                riter = riter + 1.0
                ! phase logics
                if(      istage <= PHASES(1) )then
                    iphase = 1
                else if( istage <= PHASES(2) )then
                    iphase = 2
                else if( istage <= PHASES(3) )then
                    iphase = 3
                else 
                    THROW_HARD('Invalid istage index')
                endif
                ! phase control parameters
                select case(iphase)
                    case(1)
                        refine        = 'shc_smpl'
                        rnspace       = real(NSPACE(1))
                        rmaxits       = real(MAXITS(1))
                        rmaxits_glob  = real(MAXITS_GLOB(1))
                        silence_fsc   = 'yes'
                        trs           = 0.
                        snr_noise_reg = 2.0
                        sh_first      = 'no'
                        prob_sh       = 'no'
                        ml_reg        = 'no'
                        icm           = 'no'
                    case(2)
                        refine        = 'shc_smpl'
                        rnspace       = real(NSPACE(2))
                        rmaxits       = real(MAXITS(2))
                        rmaxits_glob  = real(MAXITS_GLOB(2))
                        silence_fsc   = 'yes'
                        trs           = lpinfo(istage)%trslim
                        snr_noise_reg = 4.0
                        sh_first      = 'yes'
                        prob_sh       = 'no'
                        ml_reg        = 'yes'
                        icm           = 'no'
                    case(3)
                        refine        = 'prob'
                        rnspace       = real(NSPACE(3))
                        rmaxits       = real(MAXITS(3))
                        rmaxits_glob  = real(MAXITS_GLOB(3))
                        silence_fsc   = 'no'
                        trs           = lpinfo(istage)%trslim
                        snr_noise_reg = 6.0
                        sh_first      = 'yes'
                        prob_sh       = 'yes'
                        ml_reg        = 'yes'
                        icm           = 'yes'
                end select
                ! symmetry
                if( l_srch4symaxis )then
                    if( istage <= SYMSRCH_STAGE )then
                        ! need to replace original point-group flag with c1/pgrp_start
                        call cline_refine3D%set('pgrp', trim(params%pgrp_start))
                    else
                        call cline_refine3D%set('pgrp', trim(params%pgrp))
                    endif
                endif
                ! command line update
                call cline_refine3D%set('prg',                     'refine3D')
                call cline_refine3D%set('startit',                      riter)
                call cline_refine3D%set('which_iter',                   riter)
                call cline_refine3D%set('refine',                      refine)
                call cline_refine3D%set('lp',               lpinfo(istage)%lp)
                call cline_refine3D%set('smpd_crop', lpinfo(istage)%smpd_crop)
                call cline_refine3D%set('box_crop',   lpinfo(istage)%box_crop)
                call cline_refine3D%set('nspace',                     rnspace)
                call cline_refine3D%set('maxits',                     rmaxits)
                call cline_refine3D%set('maxits_glob',           rmaxits_glob)
                call cline_refine3D%set('silence_fsc',            silence_fsc)
                call cline_refine3D%set('trs',                            trs)
                call cline_refine3D%set('snr_noise_reg',        snr_noise_reg)
                call cline_refine3D%set('sh_first',                  sh_first)
                call cline_refine3D%set('prob_sh',                    prob_sh)
                call cline_refine3D%set('ml_reg',                      ml_reg)
                call cline_refine3D%set('icm',                            icm)
            end subroutine set_cline_refine3D
    
            subroutine conv_eo( os )
                class(oris), intent(inout) :: os
                type(sym) :: se
                type(ori) :: o_odd, o_even
                real      :: avg_euldist, euldist
                integer   :: icls, ncls
                call se%new(params%pgrp)
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
        class(abinitio_3Dmodel_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        ! constants
        integer, parameter :: NSTAGES         = 8
        integer, parameter :: PHASES(3)       = [2,6,8]
        integer, parameter :: MAXITS(3)       = [20,15,10]
        integer, parameter :: MAXITS_GLOB(3)  = [2*20,4*15,2*10]
        integer, parameter :: NSPACE(3)       = [500,1000,2500]
        integer, parameter :: SYMSRCH_STAGE   = 3
        ! commanders
        type(refine3D_commander_distr)      :: xrefine3D_distr
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        type(postprocess_commander)         :: xpostprocess
        type(symmetrize_map_commander)      :: xsymmap
        ! command lines
        type(cmdline)                       :: cline_refine3D, cline_reconstruct3D
        type(cmdline)                       :: cline_postprocess, cline_symmap
        ! other
        character(len=:),  allocatable :: frcs_fname, vol_type, str_state, vol, vol_pproc, vol_pproc_mirr, frckind, stkpath, stk, imgkind
        integer,           allocatable :: states(:)
        real,              allocatable :: frcs_avg(:)
        type(lp_crop_inf), allocatable :: lpinfo(:)
        type(parameters)               :: params
        type(sp_project)               :: spproj
        type(convergence)              :: conv
        type(class_frcs)               :: clsfrcs
        type(image)                    :: final_vol, reprojs, noisevol
        character(len=LONGSTRLEN)      :: vol_str
        character(len=STDLEN)          :: vol_iter, vol_sym
        real    :: lpsym, lpfinal, smpd
        integer :: istage, s, filtsz, ncavgs, iter, state, box_crop
        logical :: l_err
        call cline%set('objfun',    'euclid') ! use noise normalized Euclidean distances from the start
        call cline%set('sigma_est', 'global') ! obviously
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',        'yes')
        if( .not. cline%defined('overlap')     ) call cline%set('overlap',       0.95) ! needed to prevent premature convergence
        if( .not. cline%defined('prob_athres') ) call cline%set('prob_athres',    10.)
        ! if( .not. cline%defined('stoch_update') ) call cline%set('stoch_update', 'yes') ! off 4 now
        call cline%set('stoch_update', 'no')
        if( .not. cline%defined('center')      ) call cline%set('center',        'no')
        if( .not. cline%defined('cenlp')       ) call cline%set('cenlp', CENLP_DEFAULT)
        if( .not. cline%defined('oritype')     ) call cline%set('oritype',   'ptcl3D')
        if( .not. cline%defined('pgrp')        ) call cline%set('pgrp',          'c1')
        if( .not. cline%defined('pgrp_start')  ) call cline%set('pgrp_start',    'c1')
        if( .not. cline%defined('ptclw')       ) call cline%set('ptclw',         'no')
        ! make master parameters
        if( cline%defined('update_frac') ) call cline%delete('stoch_update')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! state string
        str_state = int2str_pad(1,2)
        ! symmetry
        call set_symmetry_class_vars(params)
        ! read project
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
        if( .not. cline%defined('vol1') )then
            ! randomize projection directions
            select case(trim(params%oritype))
                case('ptcl3D')
                    if( spproj%os_ptcl3D%get_noris() < 1 )then
                        THROW_HARD('Particles could not be found in the project')
                    endif
                    vol_type = 'vol'
                    call spproj%os_ptcl3D%rnd_oris
                    ! call spproj%os_ptcl3D%set_all2single('w',1.) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                case DEFAULT
                    THROW_HARD('Unsupported ORITYPE; exec_abinitio_3Dmodel')
            end select
            call spproj%write_segment_inside(params%oritype, params%projfile)
            ! create noise starting volume(s)
            call noisevol%new([params%box,params%box,params%box], params%smpd)
            do s = 1, params%nstates
                call noisevol%ran()
                vol = 'startvol_state'//int2str_pad(s,2)//'.mrc'
                call cline%set('vol'//int2str(s), vol)
                params%vols(s) = vol
                call noisevol%write(vol)
                call noisevol%ran()
                vol = 'startvol_state'//int2str_pad(s,2)//'_even.mrc'
                call noisevol%write(vol)
                vol = 'startvol_state'//int2str_pad(s,2)//'_even_unfil.mrc'
                call noisevol%write(vol)
                call noisevol%ran()
                vol = 'startvol_state'//int2str_pad(s,2)//'_odd.mrc'
                call noisevol%write(vol)
                vol = 'startvol_state'//int2str_pad(s,2)//'_odd_unfil.mrc'
                call noisevol%write(vol)
            end do
            call noisevol%kill
        endif
        ! retrieve FRC info
        call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
        if( .not.file_exists(frcs_fname) )then
            ! 08/24 This is a backwards compatibility patch to account for error in metadata
            ! on exit of streaming related to GUI directory structure (now fixed and cf above get_cavgs_stk).
            ! Will need to harmonize (move to absolute path?).
            call spproj%get_cavgs_stk(stk, ncavgs, smpd, imgkind=imgkind, stkpath=stkpath)
            frcs_fname = trim(stkpath)//'/'//trim(frcs_fname)
            if( .not.file_exists(frcs_fname) )then
                THROW_HARD('the project file does not contain an FRCs file, which is required')
            endif
        endif
        ! work out low-pass limits and downscaling parameters
        params%frcs = trim(frcs_fname)
        call clsfrcs%read(frcs_fname)
        filtsz = clsfrcs%get_filtsz()
        allocate(frcs_avg(filtsz), source=0.)
        states  = nint(spproj%os_cls2D%get_all('state'))
        call clsfrcs%avg_frc_getter(frcs_avg, states)
        allocate(lpinfo(NSTAGES))
        lpfinal = max(LPSTOP_LB,calc_lplim_final_stage(spproj,3))
        call lpstages(params%box, NSTAGES, frcs_avg, params%smpd, LPSTART_LB, LPSTART_DEFAULT, lpfinal, lpinfo, verbose=.true.)
        ! dimensions defaults
        params%box          = spproj%get_box()
        params%smpd_crop    = params%smpd
        params%box_crop     = params%box
        ! prepare command lines from prototype
        cline_refine3D      = cline
        cline_reconstruct3D = cline
        cline_postprocess   = cline
        cline_symmap        = cline
        ! map symmetrization
        if( l_srch4symaxis )then
            call cline_symmap%set('prg',     'symaxis_search') ! needed for cluster exec
            call cline_symmap%set('pgrp',   trim(params%pgrp))
            call cline_symmap%set('projfile', params%projfile)
            if( .not. cline_symmap%defined('cenlp') ) call cline_symmap%set('cenlp', CENLP_DEFAULT)
            call cline_symmap%set('hp',             params%hp)
        endif
        ! re-reconstruction & re-projection
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
        ! Frequency marching
        do istage = 1, NSTAGES
            write(logfhandle,'(A)')'>>>'
            write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ', istage,' WITH LP =', lpinfo(istage)%lp
            call set_cline_refine3D(istage)
            if( lpinfo(istage)%l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',lpinfo(istage)%box_crop
            endif
            ! Execution
            call exec_refine3D(iter)
            ! Symmetrization
            if( istage == SYMSRCH_STAGE )then
                if( l_symran )then
                    call se1%symrandomize(spproj%os_ptcl3D)
                    call spproj%write_segment_inside('ptcl3D', params%projfile)
                endif
                if( l_srch4symaxis )then
                    vol_iter = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
                    if( .not. file_exists(vol_iter) ) THROW_HARD('input volume to map symmetrization does not exist')
                    call cline_symmap%set('vol1',               trim(vol_iter))
                    call cline_symmap%set('smpd',     lpinfo(istage)%smpd_crop)
                    call cline_symmap%set('box', real(lpinfo(istage)%box_crop))
                    vol_sym  = 'symmetrized_map'//trim(params%ext)
                    call cline_symmap%set('outvol', trim(vol_sym))
                    lpsym = max(LPSYMSRCH_LB,lpinfo(SYMSRCH_STAGE)%lp)
                    call cline_symmap%set('lp', lpsym)
                    write(logfhandle,'(A,F5.1)') '>>> DID SET MAP SYMMETRIZATION LOW-PASS LIMIT (IN A) TO: ', lpsym
                    write(logfhandle,'(A)') '>>>'
                    write(logfhandle,'(A)') '>>> MAP SYMMETRIZATION'
                    write(logfhandle,'(A)') '>>>'
                    call xsymmap%execute_shmem(cline_symmap)
                    call del_file('SYMAXIS_SEARCH_FINISHED')
                    call simple_copy_file(vol_sym, vol_iter)
                endif
            endif
        enddo
        ! for visualization
        do state = 1, params%nstates
            str_state = int2str_pad(state,2)
            box_crop = lpinfo(NSTAGES)%box_crop 
            call final_vol%new([box_crop,box_crop,box_crop],lpinfo(NSTAGES)%smpd_crop)
            call final_vol%read(trim(VOL_FBODY)//trim(str_state)//trim(params%ext))
            call final_vol%generate_orthogonal_reprojs(reprojs)
            call reprojs%write_jpg('orthogonal_reprojs_state'//trim(str_state)//'.jpg')
            call final_vol%kill
            call reprojs%kill
        enddo
        ! Final reconstruction at original scale
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
        write(logfhandle,'(A)') '>>>'
        ! no ML-filtering
        call cline_reconstruct3D%set('ml_reg',      'no')
        call cline_reconstruct3D%set('needs_sigma', 'no')
        call cline_reconstruct3D%set('objfun',      'cc')
        ! no fractional or stochastic updates
        call cline_reconstruct3D%delete('update_frac')
        call cline_reconstruct3D%delete('stoch_update')
        ! individual particles reconstruction
        call cline_reconstruct3D%set('projrec', 'no')
        ! reconstruction
        call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
        call spproj%read_segment('out',params%projfile)
        do state = 1, params%nstates
            str_state = int2str_pad(state,2)
            vol       = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
            call spproj%add_vol2os_out(vol, params%smpd, state, vol_type)
            if( trim(params%oritype).eq.'ptcl3D' )then
                call spproj%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), state, params%box)
            endif
        enddo
        call spproj%write_segment_inside('out',params%projfile)
        ! post-processing
        do state = 1, params%nstates
            call cline_postprocess%delete('lp') ! so as to obtain optimal filtration
            call cline_postprocess%set('state', real(state))
            call xpostprocess%execute(cline_postprocess)
        enddo
        do state = 1, params%nstates
            str_state      = int2str_pad(state,2)
            vol            = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
            vol_pproc      = add2fbody(vol,params%ext,PPROC_SUFFIX)
            vol_pproc_mirr = add2fbody(vol,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX))
            if( file_exists(vol)            ) call simple_rename(vol,            trim(REC_FBODY)           //trim(str_state)//trim(params%ext))
            if( file_exists(vol_pproc)      ) call simple_rename(vol_pproc,      trim(REC_PPROC_FBODY)     //trim(str_state)//trim(params%ext))
            if( file_exists(vol_pproc_mirr) ) call simple_rename(vol_pproc_mirr, trim(REC_PPROC_MIRR_FBODY)//trim(str_state)//trim(params%ext))
        enddo
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
                if( istage <= NSTAGES )then
                    stage = '_stage_'//int2str(istage)
                    do state = 1, params%nstates
                        str_state = int2str_pad(state,2)
                        vol       = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
                        vol_pproc = add2fbody(vol,params%ext,PPROC_SUFFIX)
                        if( file_exists(vol)      ) call simple_copy_file(vol,       add2fbody(vol,      params%ext,stage))
                        if( file_exists(vol_pproc)) call simple_copy_file(vol_pproc, add2fbody(vol_pproc,params%ext,stage))
                    enddo
                endif
            end subroutine exec_refine3D

            subroutine set_cline_refine3D( istage )
                integer, intent(in) :: istage
                character(len=:), allocatable :: silence_fsc, sh_first, prob_sh, ml_reg, refine, icm
                integer :: iphase, s
                real    :: trs, rnspace, rmaxits, rmaxits_glob, riter
                if( istage > 1 )then ! use the previous volume rather than the noise starting volume
                    do s = 1, params%nstates
                        vol_str   = 'vol'//trim(int2str(s))
                        call cline_refine3D%delete(vol_str)
                        str_state = int2str_pad(s,2)
                        vol       = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
                        call cline_refine3D%set(vol_str, vol)
                    enddo
                endif
                ! iteration number bookkeeping
                if( cline_refine3D%defined('endit') )then
                    riter = cline_refine3D%get_rarg('endit')
                else
                    riter = 0.
                endif
                riter = riter + 1.0
                ! phase logics
                if(      istage <= PHASES(1) )then
                    iphase = 1
                else if( istage <= PHASES(2) )then
                    iphase = 2
                else if( istage <= PHASES(3) )then
                    iphase = 3
                else 
                    THROW_HARD('Invalid istage index')
                endif
                ! phase control parameters
                select case(iphase)
                    case(1)
                        refine        = 'shc_smpl'
                        rnspace       = real(NSPACE(1))
                        rmaxits       = real(MAXITS(1))
                        rmaxits_glob  = real(MAXITS_GLOB(1))
                        silence_fsc   = 'yes'
                        trs           = 0.
                        sh_first      = 'no'
                        prob_sh       = 'no'
                        ml_reg        = 'no'
                        icm           = 'no'
                    case(2)
                        refine        = 'shc_smpl'
                        rnspace       = real(NSPACE(2))
                        rmaxits       = real(MAXITS(2))
                        rmaxits_glob  = real(MAXITS_GLOB(2))
                        silence_fsc   = 'yes'
                        trs           = lpinfo(istage)%trslim
                        sh_first      = 'yes'
                        prob_sh       = 'no'
                        ml_reg        = 'yes'
                        icm           = 'no'
                    case(3)
                        refine        = 'prob'
                        rnspace       = real(NSPACE(3))
                        rmaxits       = real(MAXITS(3))
                        rmaxits_glob  = real(MAXITS_GLOB(3))
                        silence_fsc   = 'no'
                        trs           = lpinfo(istage)%trslim
                        sh_first      = 'yes'
                        prob_sh       = 'yes'
                        ml_reg        = 'yes'
                        icm           = 'yes'
                end select
                ! symmetry
                if( l_srch4symaxis )then
                    if( istage <= SYMSRCH_STAGE )then
                        ! need to replace original point-group flag with c1/pgrp_start
                        call cline_refine3D%set('pgrp', trim(params%pgrp_start))
                    else
                        call cline_refine3D%set('pgrp', trim(params%pgrp))
                    endif
                endif
                ! command line update
                call cline_refine3D%set('prg',                     'refine3D')
                call cline_refine3D%set('startit',                      riter)
                call cline_refine3D%set('which_iter',                   riter)
                call cline_refine3D%set('refine',                      refine)
                call cline_refine3D%set('lp',               lpinfo(istage)%lp)
                call cline_refine3D%set('smpd_crop', lpinfo(istage)%smpd_crop)
                call cline_refine3D%set('box_crop',   lpinfo(istage)%box_crop)
                call cline_refine3D%set('nspace',                     rnspace)
                call cline_refine3D%set('maxits',                     rmaxits)
                call cline_refine3D%set('maxits_glob',           rmaxits_glob)
                call cline_refine3D%set('silence_fsc',            silence_fsc)
                call cline_refine3D%set('trs',                            trs)
                call cline_refine3D%set('sh_first',                  sh_first)
                call cline_refine3D%set('prob_sh',                    prob_sh)
                call cline_refine3D%set('ml_reg',                      ml_reg)
                call cline_refine3D%set('icm',                            icm)                
            end subroutine set_cline_refine3D

    end subroutine exec_abinitio_3Dmodel

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio_3Dmodel2( self, cline )
        use simple_convergence, only: convergence
        use simple_fsc,         only: plot_fsc
        class(abinitio_3Dmodel2_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        real,    parameter :: SCALEFAC        = 0.667
        real,    parameter :: CENLP_DEFAULT   = 30.
        real,    parameter :: LP_DEFAULT      = 6.
        real,    parameter :: LPSTART_DEFAULT = 30., LPSTOP_DEFAULT=LP_DEFAULT
        integer, parameter :: NPARTS  = 4
        integer, parameter :: MINBOX  = 64
        integer, parameter :: NSTAGES_DEFAULT = 22
        integer, parameter :: MAXITS_SHORT = 5
        integer, parameter :: NSPACE1 = 500, NSPACE2 = 1000, NSPACE3 = 2000
        integer, parameter :: SYMSEARCH_DEFAULT = 5
        integer, parameter :: MLREG_ITER        = 1
        integer, parameter :: SHIFT_STAGE_DEFAULT = NSTAGES_DEFAULT-5 ! in [1;NSTAGES+1]
        ! commanders
        type(refine3D_commander_distr)      :: xrefine3D_distr
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        type(postprocess_commander)         :: xpostprocess
        type(symaxis_search_commander)      :: xsymsrch
        type(calc_pspec_commander_distr)    :: xcalc_pspec_distr
        ! command lines
        type(cmdline)                 :: cline_refine3D, cline_reconstruct3D, cline_reconstruct3D_mlreg
        type(cmdline)                 :: cline_postprocess, cline_symsrch, cline_calc_pspec_distr
        ! other
        type(parameters)              :: params
        type(sp_project)              :: spproj, spproj_part
        type(convergence)             :: conv
        type(sym)                     :: se1, se2
        type(class_frcs)              :: clsfrcs
        type(image)                   :: vol_even, vol_odd, reprojs, tmpvol, vol
        type(qsys_env)                :: qenv
        real,             allocatable :: fsc(:), res(:)
        character(len=:), allocatable :: str_state, vol_pproc, vol_pproc_mirr
        character(len=:), allocatable :: stack_name, dir, fsc_fname
        integer,          allocatable :: states(:), tmp(:), iters(:), prev_iters(:)
        character(len=STDLEN), allocatable :: completion_fnames(:)
        character(len=LONGSTRLEN)     :: vol_str
        real    :: lps(NSTAGES_DEFAULT), smpds(NSTAGES_DEFAULT), trs(NSTAGES_DEFAULT)
        integer :: boxs(NSTAGES_DEFAULT)
        real    :: smpd_target, lp_target, scale, trslim, cenlp, symlp, dummy, msk
        integer :: it, prev_box_crop, maxits, nptcls_sel, filtsz
        integer :: nstages, symsearch_iter, istk, part, iter, nptcls_part, i,j, cnt
        logical :: l_autoscale, l_lpset, l_err, l_srch4symaxis, l_symran, l_sym, l_lpstop_set
        logical :: l_lpstart_set
        call cline%set('oritype',      'ptcl3D')
        call cline%set('ml_reg',       'yes')
        call cline%set('stoch_update', 'no')
        call cline%set('icm',          'no')
        if( .not. cline%defined('mkdir')        ) call cline%set('mkdir',        'yes')
        if( .not. cline%defined('refine')       ) call cline%set('refine',      'prob')
        if( .not. cline%defined('autoscale')    ) call cline%set('autoscale',    'yes')
        if( .not. cline%defined('sigma_est')    ) call cline%set('sigma_est', 'global')
        if( .not. cline%defined('prob_sh')      ) call cline%set('prob_sh',      'yes')
        if( .not. cline%defined('prob_athres')  ) call cline%set('prob_athres',    10.)
        if( .not. cline%defined('center')       ) call cline%set('center',        'no')
        if( .not. cline%defined('objfun')       ) call cline%set('objfun',    'euclid')
        if( .not. cline%defined('oritype')      ) call cline%set('oritype',   'ptcl3D')
        if( .not. cline%defined('pgrp')         ) call cline%set('pgrp',          'c1')
        if( .not. cline%defined('pgrp_start')   ) call cline%set('pgrp_start',    'c1')
        if( .not. cline%defined('shift_stage')  ) call cline%set('shift_stage', SHIFT_STAGE_DEFAULT)
        if( .not. cline%defined('ptclw')        ) call cline%set('ptclw',         'no')
        if( .not. cline%defined('nparts')       ) call cline%set('nparts',      NPARTS)
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
        if( cline%defined('update_frac') ) call cline%delete('stoch_update')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call cline%delete('autoscale')
        call cline%delete('lpstart')
        call cline%delete('lpstop')
        call cline%delete('lp')
        call cline%delete('shift_stage')
        allocate(completion_fnames(params%nparts),iters(params%nparts),prev_iters(params%nparts))
        call qenv%new(1)
        str_state = int2str_pad(1,2)
        ! stages specific parameters
        nstages        = NSTAGES_DEFAULT
        symsearch_iter = SYMSEARCH_DEFAULT
        if( l_lpset )then
            params%lpstop  = params%lp
            params%lpstart = params%lp
        endif
        if( params%shift_stage < 1 .or. params%shift_stage > nstages+1 )then
            params%shift_stage = min(nstages+1,max(1,params%shift_stage))
            THROW_WARN('SHIFT_STAGE out of range, defaulting to: '//int2str(params%shift_stage))
        endif
        ! read project
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
        if( .not. cline%defined('vol1') )then
            ! randomize projection directions
            if( spproj%os_ptcl3D%get_noris() < 1 )then
                THROW_HARD('Particles could not be found in the project')
            endif
            call spproj%os_ptcl3D%rnd_oris
            call spproj%os_ptcl3D%set_all2single('w',1.)
            states = nint(spproj%os_ptcl3D%get_all('state'))
            call spproj%write_segment_inside(params%oritype, params%projfile)
        endif
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
        ! command-lines
        cline_refine3D            = cline
        cline_reconstruct3D       = cline
        cline_postprocess         = cline
        cline_symsrch             = cline
        cline_reconstruct3D_mlreg = cline_reconstruct3D
        cline_calc_pspec_distr    = cline
        call cline_refine3D%set('prg',                'refine3D')
        call cline_refine3D%set('projfile',      params%projfile)
        call cline_refine3D%set('pgrp',        params%pgrp_start)
        call cline_refine3D%set('maxits',                    999)
        call cline_refine3D%delete('nparts')
        call cline_reconstruct3D%set('prg',      'reconstruct3D')
        call cline_reconstruct3D%set('box',     real(params%box))
        call cline_reconstruct3D%set('projfile', params%projfile)
        call cline_reconstruct3D%set('ml_reg',              'no')
        call cline_reconstruct3D%set('needs_sigma',         'no')
        call cline_reconstruct3D%set('objfun',              'cc')
        call cline_reconstruct3D%set('pgrp',   params%pgrp_start)
        call cline_postprocess%set('prg',          'postprocess')
        call cline_postprocess%set('projfile',   params%projfile)
        call cline_postprocess%set('imgkind',              'vol')
        if( l_srch4symaxis )then
            call cline_symsrch%set('prg',     'symaxis_search') ! needed for cluster exec
            call cline_symsrch%set('pgrp',     params%pgrp)
            call cline_symsrch%set('projfile', params%projfile)
            call cline_symsrch%set('hp',       params%hp)
            call cline_symsrch%set('center',   'yes')
        endif
        call cline_reconstruct3D_mlreg%set('prg',         'reconstruct3D')
        call cline_reconstruct3D_mlreg%set('objfun',      'euclid')
        call cline_reconstruct3D_mlreg%set('needs_sigma', 'yes')
        call cline_reconstruct3D_mlreg%set('sigma_est',   params%sigma_est)
        call cline_reconstruct3D_mlreg%set('ml_reg',      'yes')
        call cline_calc_pspec_distr%set('prg',      'calc_pspec')
        ! Frequency marching plan
        lps(1) = params%lpstart
        do it = 2,nstages-1
            lps(it) = params%lpstop + (params%lpstart-params%lpstop) * real(nstages-it) / real(nstages)
        enddo
        lps(nstages) = params%lpstop
        if( l_lpset )then
            ! from mskdiam2lplimits lpstart above
        else
            symlp = (params%lpstart+params%lpstop)/2.
        endif
        ! dimensions
        do it = 1,nstages
            lp_target   = lps(it) * SCALEFAC
            smpd_target = max(params%smpd, lp_target/2.)
            call autoscale(params%box, params%smpd, smpd_target, boxs(it), smpds(it), scale, minbox=MINBOX)
            if( it < params%shift_stage )then
                trs(it) = 0.
            else
                trs(it) = max(2.0, AHELIX_WIDTH / smpds(it) / 2.0)
            endif
        enddo
        ! random reconstruction
        params%smpd_crop = smpds(1)
        params%box_crop  = boxs(1)
        call cline_reconstruct3D%set('smpd_crop', params%smpd_crop)
        call cline_reconstruct3D%set('box_crop',  params%box_crop)
        if( params%l_ml_reg .and. MLREG_ITER==1 )then
            call xcalc_pspec_distr%execute_shmem(cline_calc_pspec_distr)
            call cline_reconstruct3D%set('which_iter', 1)
            call cline_reconstruct3D%set('ml_reg',     'yes')
            call cline_reconstruct3D%set('needs_sigma','yes')
            call cline_reconstruct3D%set('objfun',     'euclid')
        endif
        call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
        call spproj%read_segment('ptcl3D', params%projfile)
        call spproj%os_ptcl3D%set_all2single('updatecnt',0.)
        call spproj%write_segment_inside('ptcl3D', params%projfile)
        call spproj%read_segment('out', params%projfile)
        call spproj%add_vol2os_out('recvol_state01.mrc', params%smpd_crop, 1, 'vol')
        call spproj%write_segment_inside('out', params%projfile)
        ! updating stack names to absolute path
        call spproj%read_segment('stk', params%projfile)
        do istk = 1,spproj%os_stk%get_noris()
            stack_name = trim(spproj%get_stkname(istk))
            stack_name = simple_abspath(stack_name, check_exists=.false.)
            call spproj%os_stk%set(istk, 'stk', stack_name)
        enddo
        call spproj%write_segment_inside('stk', params%projfile)
        call spproj%read_segment('ptcl2D', params%projfile)
        call spproj%read_segment('ptcl3D', params%projfile)
        ! directory structure
        do part = 1,params%nparts
            dir = int2str(part)//'/'
            call simple_mkdir(dir)
        enddo
        ! Parts partitioning
        call cline_refine3D%delete('projfile')
        call cline_refine3D%set('projname', get_fbody(basename(params%projfile), 'simple'))
        call cline_refine3D%set('projfile', basename(params%projfile))
        nptcls_sel  = count(states==1)
        nptcls_part = ceiling(real(nptcls_sel)/real(params%nparts))
        tmp = states
        j   = 0
        do part = 1,params%nparts
            spproj_part%os_stk    = spproj%os_stk
            spproj_part%os_ptcl2D = spproj%os_ptcl2D
            spproj_part%os_ptcl3D = spproj%os_ptcl3D
            spproj_part%projinfo  = spproj%projinfo
            spproj_part%compenv   = spproj%compenv
            tmp = states
            if( j > 0 ) tmp(1:j) = 0
            cnt = 0
            do i = j+1,params%nptcls
                if( states(i) == 1)then
                    cnt       = cnt+1
                    states(i) = part
                    if( cnt == nptcls_part )then
                        j = i
                        exit
                    endif
                endif
            enddo
            if( part < params%nparts ) tmp(j+1:) = 0
            call spproj_part%os_ptcl2D%set_all('state', real(tmp))
            call spproj_part%os_ptcl3D%set_all('state', real(tmp))
            call spproj_part%prune_particles
            call chdir(int2str(part)//'/')
            call spproj_part%update_projinfo(cline_refine3D)
            call spproj_part%write
            call chdir('../')
            completion_fnames(part) = int2str(part)//'/'//trim(JOB_FINISHED_FBODY)
            call spproj_part%kill
        enddo
        deallocate(tmp)
        ! Stages loop
        iters(:) = 0
        do it = 1,nstages
            params%smpd_crop = smpds(it)
            params%box_crop  = boxs(it)
            params%lp        = lps(it)
            params%trs       = trs(it)
            write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ',it,' WITH LP =',params%lp
            write(logfhandle,'(A,I3)')        '>>> CROPPED IMAGE SIZE: ',params%box_crop
            call cline_refine3D%set('smpd_crop', params%smpd_crop)
            call cline_refine3D%set('box_crop',  params%box_crop)
            call cline_refine3D%set('lp',        params%lp)
            if( it == 1 )then
                call cline_refine3D%set('vol1', '../recvol_state01.mrc')
            else
                call cline_refine3D%set('vol1', '../recvol_state01_stage'//int2str_pad(it-1,2)//'.mrc')
            endif
            ! # of iterations
            call cline_refine3D%set('maxits', MAXITS_SHORT)
            ! projection directions & shift
            if( it < params%shift_stage )then
                call cline_refine3D%set('nspace', NSPACE1)
            else
                call cline_refine3D%set('nspace', NSPACE2)
            end if
            if( it >= nstages-2 ) call cline_refine3D%set('nspace', NSPACE3)
            call cline_refine3D%set('trs', params%trs)
            ! execution
            do part = 1,params%nparts
                call exec_refine3D(part)
            enddo
            ! waiting
            call qsys_watcher(completion_fnames)
            ! convergence, volume averaging, padding & cleanup
            prev_iters = iters
            call vol%new([params%box_crop, params%box_crop,params%box_crop],params%smpd_crop)
            call vol_even%new([params%box_crop, params%box_crop,params%box_crop],params%smpd_crop)
            call vol_odd%new([params%box_crop, params%box_crop,params%box_crop],params%smpd_crop)
            call tmpvol%new([params%box_crop, params%box_crop,params%box_crop],params%smpd_crop)
            do part = 1,params%nparts
                call chdir(int2str(part)//'/')
                ! convergence parameters
                call conv%read(l_err)
                iters(part) = nint(conv%get('iter'))
                write(logfhandle,'(A,I3,A,F7.3,A,F7.3)')'>>> PART ',part,'; PROJ OVERLAP: ',&
                    &conv%get('mi_proj'),'; SCORE: ',conv%get('score')
                ! volumes
                if( part == 1)then
                    call vol%read('recvol_state01_iter'//int2str_pad(iters(part),3)//'.mrc')
                    call vol_even%read('recvol_state01_iter'//int2str_pad(iters(part),3)//'_even.mrc')
                    call vol_odd%read('recvol_state01_iter'//int2str_pad(iters(part),3)//'_odd.mrc')
                else
                    call tmpvol%read('recvol_state01_iter'//int2str_pad(iters(part),3)//'.mrc')
                    call vol%add(tmpvol)
                    call tmpvol%read('recvol_state01_iter'//int2str_pad(iters(part),3)//'_even.mrc')
                    call vol_even%add(tmpvol)
                    call tmpvol%read('recvol_state01_iter'//int2str_pad(iters(part),3)//'_odd.mrc')
                    call vol_odd%add(tmpvol)
                endif
                ! cleanup
                call qsys_cleanup
                call del_files(DIST_FBODY,      1,ext='.dat')
                call del_files(ASSIGNMENT_FBODY,1,ext='.dat')
                call del_file(trim(DIST_FBODY)      //'.dat')
                call del_file(trim(ASSIGNMENT_FBODY)//'.dat')
                call del_file(JOB_FINISHED_FBODY)
                call del_file(trim(FSC_FBODY)//int2str_pad(1,2)//BIN_EXT)
                do i = prev_iters(part)+1,iters(part)-1
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//PPROC_SUFFIX//params%ext)
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//LP_SUFFIX//params%ext)
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//'_even'//params%ext)
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//'_odd'//params%ext)
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//'_even_unfil'//params%ext)
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//'_odd_unfil'//params%ext)
                    call del_file(trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(i,3)//params%ext)
                enddo
                call chdir('../')
            enddo
            ! averaging & fsc
            if( it < NSTAGES_DEFAULT )then
                ! Volume & FSC will be padded on the fly at the next refine3D run
                call vol%div(real(params%nparts))
                call vol_even%div(real(params%nparts))
                call vol_odd%div(real(params%nparts))
                call vol%write(trim(VOL_FBODY)//trim(str_state)//'_stage'//int2str_pad(it,2)//'.mrc')
                filtsz = fdim(params%box_crop) - 1
                msk    = real(params%box_crop / 2) - COSMSKHALFWIDTH - 1.
                allocate(fsc(filtsz),source=0.)
                call vol_even%mask(msk, 'soft', backgr=0.)
                call vol_odd%mask(msk, 'soft', backgr=0.)
                call vol_even%fft()
                call vol_odd%fft()
                call vol_even%fsc(vol_odd, fsc)
                fsc_fname = trim(FSC_FBODY)//int2str_pad(1,2)//BIN_EXT
                call arr2file(fsc, fsc_fname)
                call cline_refine3D%set('fsc', '../'//trim(fsc_fname))
                res = get_resarr(params%box_crop, params%smpd_crop)
                call plot_fsc(size(fsc), fsc, res, params%smpd_crop, 'fsc_stage_'//int2str_pad(it,2))
                deallocate(fsc,res)
                call vol%kill
                call vol_even%kill
                call vol_odd%kill
                call tmpvol%kill
            endif
            ! symmetrization
            if( it == SYMSEARCH_ITER-1 )then
                call consolidate_alnparms
                call cline_symsrch%set('vol1', trim(VOL_FBODY)//trim(str_state)//'_stage'//int2str_pad(it,2)//'.mrc')
                call symmetrize
                call cline_refine3D%set('pgrp', params%pgrp)
                call cline_reconstruct3D%set('pgrp', params%pgrp)
                l_srch4symaxis = .false.
                l_symran       = .false.
                ! transfer symmetrized parameters
                call spproj%read_segment('ptcl3D',params%projfile)
                do part = 1,params%nparts
                    call spproj_part%read_segment('ptcl3D', int2str(part)//'/'//basename(params%projfile))
                    j = 0
                    do i = 1,params%nptcls
                        if( states(i) /= part ) cycle
                        j = j+1
                        call spproj_part%os_ptcl3D%transfer_3Dparams(j, spproj%os_ptcl3D, i)
                    enddo
                    call spproj_part%write_segment_inside('ptcl3D',int2str(part)//'/'//basename(params%projfile))
                enddo
                call spproj_part%kill
            endif
        enddo
        ! gathering alignment parameters
        call consolidate_alnparms
        ! final reconstruction
        write(logfhandle,'(A)') '>>>'
        write(logfhandle,'(A)') '>>> RECONSTRUCTION AT ORIGINAL SAMPLING'
        write(logfhandle,'(A)') '>>>'
        ! no ML-filtering
        call cline_reconstruct3D%set('ml_reg',      'no')
        call cline_reconstruct3D%set('needs_sigma', 'no')
        call cline_reconstruct3D%set('objfun',      'cc')
        call cline_reconstruct3D%delete('smpd_crop')
        call cline_reconstruct3D%delete('box_crop')
        call xreconstruct3D_distr%execute_shmem(cline_reconstruct3D)
        vol_str = trim(VOL_FBODY)//trim(str_state)//trim(params%ext)
        call spproj%read_segment('out',params%projfile)
        call spproj%add_vol2os_out(vol_str, params%smpd, 1, 'vol')
        call spproj%add_fsc2os_out(FSC_FBODY//str_state//trim(BIN_EXT), 1, params%box)
        call spproj%write_segment_inside('out',params%projfile)
        ! post-processing
        call cline_postprocess%delete('lp')
        call cline_postprocess%set('state', 1)
        call xpostprocess%execute(cline_postprocess)
        vol_pproc      = add2fbody(vol_str,params%ext,PPROC_SUFFIX)
        vol_pproc_mirr = add2fbody(vol_str,params%ext,trim(PPROC_SUFFIX)//trim(MIRR_SUFFIX))
        call simple_rename(vol_str, trim(REC_FBODY)//trim(str_state)//trim(params%ext))
        if(file_exists(vol_pproc)     ) call simple_rename(vol_pproc,      trim(REC_PPROC_FBODY)     //trim(str_state)//trim(params%ext))
        if(file_exists(vol_pproc_mirr)) call simple_rename(vol_pproc_mirr, trim(REC_PPROC_MIRR_FBODY)//trim(str_state)//trim(params%ext))
        ! cleanup
        call spproj%kill
        call qsys_cleanup
        call simple_end('**** SIMPLE_ABINITIO_3DMODEL2 NORMAL STOP ****')
        contains

            subroutine exec_refine3D( part )
                integer,          intent(in)  :: part
                character(len=XLONGSTRLEN) :: cwd
                call cline_refine3D%set('startit', iters(part)+1)
                dir = int2str(part)//'/'
                call chdir(dir)
                call simple_getcwd(cwd)
                cwd_glob = trim(cwd)
                call qenv%new(1)
                call qenv%exec_simple_prg_in_queue_async(cline_refine3D, './refine3D', 'log_refine3D')
                call chdir('../')
                call simple_getcwd(cwd_glob)
            end subroutine exec_refine3D

            subroutine symmetrize()
                use simple_projector_hlev, only: rotvol_slim
                use simple_projector,      only: projector
                type(projector) :: vol_pad
                type(image)     :: rovol_pad, rovol
                type(ori)       :: o
                real    :: symaxis_rmat(3,3), symop_rmat(3,3), rmat(3,3)
                integer :: ldim_pd(3), boxpd,isym, nsym
                if( l_symran )then
                    call spproj%read_segment(params%oritype, params%projfile)
                    call se1%symrandomize(spproj%os_ptcl3D)
                    call spproj%write_segment_inside(params%oritype, params%projfile)
                endif
                if( l_srch4symaxis )then
                    write(logfhandle,'(A)') '>>>'
                    write(logfhandle,'(A)') '>>> SYMMETRY AXIS SEARCH'
                    write(logfhandle,'(A)') '>>>'
                    symlp = max(symlp, params%lp)
                    call cline_symsrch%set('lp',       symlp)
                    call cline_symsrch%set('box_crop', params%box_crop)
                    call xsymsrch%execute_shmem(cline_symsrch)
                    call del_file('SYMAXIS_SEARCH_FINISHED')
                    ! symmetrize volume
                    call vol%new([params%box_crop, params%box_crop,params%box_crop],params%smpd_crop)
                    call rovol%new([params%box_crop, params%box_crop,params%box_crop],params%smpd_crop)
                    boxpd   = 2 * round2even(KBALPHA * real(params%box_crop))
                    ldim_pd = [boxpd,boxpd,boxpd]
                    call rovol_pad%new(ldim_pd, params%smpd_crop)
                    call vol_pad%new(ldim_pd, params%smpd_crop)
                    call vol%read('vol_aligned2_'//trim(params%pgrp)//'axis'//params%ext)
                    call vol%pad(vol_pad)
                    call vol_pad%fft
                    call vol_pad%expand_cmat(KBALPHA)
                    nsym = se2%get_nsym()
                    do isym =2,nsym
                        call se2%get_symori(isym, o)
                        call rotvol_slim(vol_pad, rovol_pad, rovol, o)
                        call vol%add_workshare(rovol)
                    end do
                    call vol%div(real(nsym))
                    call vol%write(trim(VOL_FBODY)//trim(str_state)//'_stage'//int2str_pad(it,2)//'.mrc')
                    call o%kill
                    call rovol%kill
                    call rovol_pad%kill
                    call vol%kill
                    call vol_pad%kill
                endif
            end subroutine symmetrize

            subroutine consolidate_alnparms
                integer :: i,j,part
                do part = 1,params%nparts
                    call spproj_part%read_segment('ptcl3D', int2str(part)//'/'//basename(params%projfile))
                    j = 0
                    do i = 1,params%nptcls
                        if( states(i) /= part ) cycle
                        j = j+1
                        call spproj%os_ptcl3D%transfer_3Dparams(i, spproj_part%os_ptcl3D, j)
                    enddo
                enddo
                call spproj_part%kill
                call spproj%write_segment_inside('ptcl3D',params%projfile)
            end subroutine consolidate_alnparms

    end subroutine exec_abinitio_3Dmodel2

end module simple_commander_abinitio
