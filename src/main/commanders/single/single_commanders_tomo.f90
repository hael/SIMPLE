!@descr: SINGLE tilt-series commanders for atomic electron tomography: import with angle conventions, staged pose refinement
! import_tilt_series makes a refine=cont-ready project from an MRC stack and an angle file (SIMPLE, single-axis tilt or
! GENFIRE/RESIRE convention); refine_tomo runs gridding, a density mask, masked PCG and a short refine=cont pass at a
! fixed low-pass. Plan: doc/refactoring_notes/planned/refine_tomo.md
module single_commanders_tomo
use simple_commanders_api
implicit none

public :: commander_import_tilt_series, commander_refine_tomo
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_import_tilt_series
  contains
    procedure :: execute => exec_import_tilt_series
end type commander_import_tilt_series

type, extends(commander_base) :: commander_refine_tomo
  contains
    procedure :: execute => exec_refine_tomo
end type commander_refine_tomo

contains

    !> A tilt series as a project refine=cont continues from: every projection active, its half-set by index
    !! parity and its projection direction the nearest node of an nspace spiral; poses from angtab in angconv
    subroutine exec_import_tilt_series( self, cline )
        class(commander_import_tilt_series), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters)  :: params
        type(sp_project)  :: spproj
        type(oris)        :: os, grid
        type(ctfparams)   :: ctfvars
        real, allocatable :: angs(:,:), shifts(:,:)
        integer           :: ldim(3), nimgs, i
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('angtab') ) THROW_HARD('angtab (per-projection angles) is required; import_tilt_series')
        call cline%set('ctf', 'no')
        call params%new(cline)
        call find_ldim_nptcls(params%stk, ldim, nimgs)
        call read_tilt_angles(params%angtab, params%angconv, angs, shifts)
        if( size(angs,2) /= nimgs ) THROW_HARD('angtab and stk hold different numbers of projections; import_tilt_series')
        call os%new(nimgs, is_ptcl=.true.)
        do i = 1, nimgs
            call os%set_euler(i, angs(:,i))
            call os%set_shift(i, shifts(:,i))
            call os%set_state(i, 1)
            call os%set(i, 'eo', mod(i-1, 2))
        end do
        call grid%new(params%nspace, is_ptcl=.false.)
        call grid%spiral()
        call os%set_projs(grid)
        ctfvars%ctfflag = CTFFLAG_NO
        ctfvars%smpd    = params%smpd
        ctfvars%kv      = params%kv
        ctfvars%cs      = params%cs
        ctfvars%fraca   = params%fraca
        call spproj%read(params%projfile)
        call spproj%add_single_stk(params%stk, ctfvars, os)
        call spproj%write
        write(logfhandle,'(A,I0,A)') '>>> IMPORT_TILT_SERIES: ', nimgs, ' projections, angles in the '//trim(params%angconv)// &
            &' convention'
        call spproj%kill
        call os%kill
        call grid%kill
        call simple_end('**** IMPORT_TILT_SERIES NORMAL STOP ****')
    end subroutine exec_import_tilt_series

    !> Per-projection SIMPLE Euler angles and shifts (pixels) from a text file: simple, an oritab (e1 e2 e3, x y);
    !! tilt, one tilt angle per line (single-axis about image x) with optional x y; genfire, phi theta psi in the
    !! GENFIRE/RESIRE convention with optional x y. Comment and blank lines are skipped in the numeric formats
    subroutine read_tilt_angles( fname, angconv, angs, shifts )
        use simple_ori_utils, only: genfire2euler
        class(string),     intent(in)  :: fname
        character(len=*),  intent(in)  :: angconv
        real, allocatable, intent(out) :: angs(:,:), shifts(:,:)
        integer, parameter :: TXT_READ = 1 ! nrtxtfile's OPEN_TO_READ
        type(nrtxtfile)   :: tab
        type(oris)        :: os
        real, allocatable :: line(:)
        integer           :: i, n, nrecs
        select case(trim(angconv))
            case('simple')
                n = nlines(fname)
                call os%new(n, is_ptcl=.true.)
                call os%read(fname)
                allocate(angs(3,n), shifts(2,n))
                do i = 1, n
                    angs(:,i)   = os%get_euler(i)
                    shifts(:,i) = os%get_2Dshift(i)
                end do
                call os%kill
            case('tilt', 'genfire')
                call tab%new(fname, TXT_READ)
                n     = tab%get_ndatalines()
                nrecs = tab%get_nrecs_per_line()
                if( trim(angconv) == 'tilt' )then
                    if( nrecs /= 1 .and. nrecs /= 3 ) THROW_HARD('angconv=tilt reads tilt [x y] per line; read_tilt_angles')
                else
                    if( nrecs /= 3 .and. nrecs /= 5 ) THROW_HARD('angconv=genfire reads phi theta psi [x y] per line; read_tilt_angles')
                endif
                allocate(line(nrecs), angs(3,n), shifts(2,n), source=0.)
                do i = 1, n
                    call tab%readNextDataLine(line)
                    if( trim(angconv) == 'tilt' )then
                        angs(:,i) = genfire2euler([line(1), 0., 0.])
                        if( nrecs == 3 ) shifts(:,i) = line(2:3)
                    else
                        angs(:,i) = genfire2euler(line(1:3))
                        if( nrecs == 5 ) shifts(:,i) = line(4:5)
                    endif
                end do
                call tab%kill
            case default
                THROW_HARD('angconv must be simple, tilt or genfire; read_tilt_angles')
        end select
        if( n < 1 ) THROW_HARD('angtab holds no projections; read_tilt_angles')
    end subroutine read_tilt_angles

    !> Gridding from the imported poses, a density mask, masked PCG as the first reference, then refine=cont: by default
    !! one stage of two iterations at 2 A, since longer refinement drifts along the tilt-scale and frame-rotation modes
    !! (nstages > 1 marches the low-pass from lpstart to lpstop); the final map and the poses in both conventions
    subroutine exec_refine_tomo( self, cline )
        use simple_commanders_rec,      only: commander_rec3D
        use simple_commanders_refine3D, only: commander_refine3D
        use simple_ori_utils,           only: euler2genfire
        use simple_refine3D_fnames,     only: refine3D_fsc_fname, refine3D_state_halfvol_fname, refine3D_state_vol_fname
        class(commander_refine_tomo), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        character(len=*), parameter :: MASK_FNAME         = 'tomo_mask.mrc'
        character(len=*), parameter :: FINAL_FBODY        = 'tomo_final'
        real,             parameter :: LPSTART_DEFAULT    = 2.
        real,             parameter :: LPSTOP_DEFAULT     = 2.
        real,             parameter :: AMSKLP_DEFAULT     = 6.
        real,             parameter :: TRS_DEFAULT        = 3.
        real,             parameter :: ATHRES_DEFAULT     = 5.
        real,             parameter :: MASK_GROW_A        = 2.   ! envelope dilation and soft edge (A)
        integer,          parameter :: NSTAGES_DEFAULT    = 1
        integer,          parameter :: MAXITS_DEFAULT     = 2
        integer,          parameter :: MAXITS_PCG_DEFAULT = 5
        type(parameters)         :: params
        type(commander_rec3D)    :: xrec3D
        type(commander_refine3D) :: xrefine3D
        type(cmdline)            :: cline_ref
        type(string)             :: stage_vol
        real    :: lp_sched, lp_stage, res0143
        integer :: istage, iter
        ! the workflow's defaults; every one can be given on the command line
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',         LPSTART_DEFAULT)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',          LPSTOP_DEFAULT)
        if( .not. cline%defined('nstages')         ) call cline%set('nstages',         NSTAGES_DEFAULT)
        if( .not. cline%defined('maxits')          ) call cline%set('maxits',          MAXITS_DEFAULT)
        if( .not. cline%defined('maxits_pcg')      ) call cline%set('maxits_pcg',      MAXITS_PCG_DEFAULT)
        if( .not. cline%defined('amsklp')          ) call cline%set('amsklp',          AMSKLP_DEFAULT)
        if( .not. cline%defined('trs')             ) call cline%set('trs',             TRS_DEFAULT)
        if( .not. cline%defined('athres_cont')     ) call cline%set('athres_cont',     ATHRES_DEFAULT)
        call cline%set('oritype', 'ptcl3D')
        call cline%set('mkdir',   'yes')
        call cline%delete('nparts')   ! a shared-memory workflow
        call params%new(cline)        ! the execution directory and its copy of the project
        call cline%set('mkdir',   'no')
        if( params%lpstart < params%lpstop ) THROW_HARD('lpstart must be at or above lpstop (A); refine_tomo')
        if( params%nstages < 1 )             THROW_HARD('nstages must be at least 1; refine_tomo')
        ! a tilt series is sampled at sub-angstrom pixels: grow the envelope by a fixed distance, not by pixels
        if( .not. cline%defined('binwidth') ) params%binwidth = max(params%binwidth, nint(MASK_GROW_A / params%smpd_crop))
        if( .not. cline%defined('edge')     ) params%edge     = max(params%edge,     nint(MASK_GROW_A / params%smpd_crop))
        ! the first reference: gridding, its density envelope, masked PCG
        write(logfhandle,'(A)') '>>> REFINE_TOMO: GRIDDING FROM THE INPUT POSES'
        call reconstruct('gridding', .false.)
        call make_mask(refine3D_state_vol_fname(1))
        write(logfhandle,'(A)') '>>> REFINE_TOMO: MASKED PCG, THE FIRST REFERENCE'
        call reconstruct('pcg', .true.)
        stage_vol = 'tomo_stage00.mrc'
        call simple_copy_file(refine3D_state_vol_fname(1), stage_vol)
        ! pose refinement at the native sampling; with several stages the iteration count runs on across them
        ! (startit = the previous stage's endit + 1), as in the staged solve3D
        iter = 1
        do istage = 1, params%nstages
            lp_sched = scheduled_lp(istage)
            res0143  = map_res0143()
            lp_stage = max(lp_sched, res0143)
            write(logfhandle,'(A,I0,A,F6.2,A,F6.2,A,F6.2,A)') '>>> REFINE_TOMO: STAGE ', istage, ', LOW-PASS ', lp_stage, &
                &' A (SCHEDULED ', lp_sched, ' A, MAP FSC=0.143 AT ', res0143, ' A)'
            call cline_ref%kill
            call cline_ref%set('prg',             'refine3D')
            call cline_ref%set('projfile',        params%projfile)
            call cline_ref%set('vol1',            stage_vol)
            call cline_ref%set('refine',          'cont')
            call cline_ref%set('oritype',         'ptcl3D')
            call cline_ref%set('objfun',          'cc')
            call cline_ref%set('ml_reg',          'no')
            call cline_ref%set('rec_backend',     'pcg')
            call cline_ref%set('pcg_mskfile',     MASK_FNAME)
            call cline_ref%set('maxits_pcg',      params%maxits_pcg)
            call cline_ref%set('mskdiam',         params%mskdiam)
            call cline_ref%set('pgrp',            'c1')
            call cline_ref%set('lp',              lp_stage)
            call cline_ref%set('trs',             params%trs)
            call cline_ref%set('athres_cont',     params%athres_cont)
            call cline_ref%set('cont_conv_rot',   params%cont_conv_rot)
            call cline_ref%set('cont_conv_shift', params%cont_conv_shift)
            call cline_ref%set('cont_conv_frac',  params%cont_conv_frac)
            call cline_ref%set('maxits',          params%maxits)
            call cline_ref%set('startit',         iter)
            call cline_ref%set('which_iter',      iter)
            if( params%minits > 0 ) call cline_ref%set('minits', params%minits)
            call cline_ref%set('center',          'no')
            call cline_ref%set('nthr',            params%nthr)
            call cline_ref%set('mkdir',           'no')
            call xrefine3D%execute(cline_ref)
            if( .not. file_exists(refine3D_state_vol_fname(1)) ) THROW_HARD('refine3D left no map; refine_tomo')
            if( .not. cline_ref%defined('endit') ) THROW_HARD('refine3D reported no final iteration; refine_tomo')
            iter = cline_ref%get_iarg('endit') + 1
            stage_vol = 'tomo_stage'//int2str_pad(istage, 2)//'.mrc'
            call simple_copy_file(refine3D_state_vol_fname(1), stage_vol)
            call make_mask(stage_vol)
        end do
        ! the final map from the final poses, and the poses in both conventions
        write(logfhandle,'(A)') '>>> REFINE_TOMO: FINAL MASKED PCG MAP'
        call reconstruct('pcg', .true.)
        call simple_copy_file(refine3D_state_vol_fname(1),             string(FINAL_FBODY//'.mrc'))
        call simple_copy_file(refine3D_state_halfvol_fname(1, 'even'), string(FINAL_FBODY//'_even.mrc'))
        call simple_copy_file(refine3D_state_halfvol_fname(1, 'odd'),  string(FINAL_FBODY//'_odd.mrc'))
        call write_poses
        res0143 = map_res0143()
        write(logfhandle,'(A,F6.2,A)') '>>> REFINE_TOMO: FINAL HALF-MAP FSC=0.143 AT ', res0143, ' A'
        call cline_ref%kill
        call stage_vol%kill
        call simple_end('**** REFINE_TOMO NORMAL STOP ****')

    contains

        !> reconstruct3D in this directory from the project's current poses, with the mask as PCG support
        subroutine reconstruct( backend, l_mask )
            character(len=*), intent(in) :: backend
            logical,          intent(in) :: l_mask
            type(cmdline) :: cline_rec
            call cline_rec%set('prg',         'reconstruct3D')
            call cline_rec%set('projfile',    params%projfile)
            call cline_rec%set('rec_backend', backend)
            call cline_rec%set('mskdiam',     params%mskdiam)
            call cline_rec%set('pgrp',        'c1')
            call cline_rec%set('objfun',      'cc')
            call cline_rec%set('ml_reg',      'no')
            call cline_rec%set('maxits_pcg',  params%maxits_pcg)
            call cline_rec%set('nthr',        params%nthr)
            call cline_rec%set('mkdir',       'no')
            if( l_mask ) call cline_rec%set('pcg_mskfile', MASK_FNAME)
            call xrec3D%execute(cline_rec)
            call cline_rec%kill
        end subroutine reconstruct

        !> density envelope of a map (automask3D at amsklp, grown by binwidth, soft edge), written as the PCG support
        subroutine make_mask( vol_fname )
            class(string), intent(in) :: vol_fname
            type(image)     :: vol
            type(image_msk) :: msk
            integer         :: ldim(3), nvols
            call find_ldim_nptcls(vol_fname, ldim, nvols)
            call vol%new(ldim, params%smpd_crop)
            call vol%read(vol_fname)
            call msk%automask3D(params, vol, .false.)
            call msk%write(string(MASK_FNAME))
            call msk%kill
            call vol%kill
        end subroutine make_mask

        !> low-pass limit of a stage: geometric from lpstart to lpstop over nstages
        real function scheduled_lp( istage )
            integer, intent(in) :: istage
            if( params%nstages == 1 )then
                scheduled_lp = params%lpstop
            else
                scheduled_lp = params%lpstart * (params%lpstop / params%lpstart)**(real(istage - 1) / real(params%nstages - 1))
            endif
        end function scheduled_lp

        !> FSC 0.143 resolution of the last reconstruction in this directory (0 without one)
        real function map_res0143()
            real, allocatable :: fsc(:), res(:)
            real    :: res05
            integer :: n
            map_res0143 = 0.
            if( .not. file_exists(refine3D_fsc_fname(1)) ) return
            fsc = file2rarr(refine3D_fsc_fname(1))
            res = get_resarr(params%box_crop, params%smpd_crop)
            n   = min(size(fsc), size(res))
            if( n < 2 ) return
            call get_resolution(fsc(1:n), res(1:n), res05, map_res0143)
        end function map_res0143

        !> the project's final poses as an oritab and in the GENFIRE/RESIRE convention (re-importable, angconv=genfire)
        subroutine write_poses
            type(sp_project) :: spproj
            type(string)     :: fname
            real             :: ang(3), sh(2)
            integer          :: i, funit, io_stat
            call spproj%read_segment('ptcl3D', params%projfile)
            call spproj%os_ptcl3D%write(string(FINAL_FBODY//'_poses_simple.txt'))
            fname = FINAL_FBODY//'_poses_genfire.txt'
            call fopen(funit, fname, 'replace', 'unknown', iostat=io_stat, form='formatted')
            call fileiochk('refine_tomo fopen failed '//fname%to_char(), io_stat)
            write(funit,'(A)') '# phi theta psi (degrees, GENFIRE/RESIRE convention) x y (SIMPLE shifts, pixels), stack order'
            do i = 1, spproj%os_ptcl3D%get_noris()
                ang = euler2genfire(spproj%os_ptcl3D%get_euler(i))
                sh  = spproj%os_ptcl3D%get_2Dshift(i)
                write(funit,'(5F12.5)') ang, sh
            end do
            call fclose(funit)
            call fname%kill
            call spproj%kill
        end subroutine write_poses

    end subroutine exec_refine_tomo

end module single_commanders_tomo
