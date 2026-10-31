!@descr: the trajectory_solve3D commander: a tracked nanoparticle's graphene-subtracted frames to its time-integrated map in one program
module single_commanders_trajectory_solve3D
use simple_commanders_api
use simple_commanders_reproject, only: commander_reproject
use single_commanders_nano3D,    only: commander_refine3D_nano, commander_autorefine3D_nano
use simple_refine3D_fnames,      only: refine3D_state_vol_fbody
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_trajectory_solve3D
  contains
    procedure :: execute      => exec_trajectory_solve3D
end type commander_trajectory_solve3D

contains

    !> trajectory_solve3D: the graphene-subtracted frames of a tracked particle to its time-integrated map, one stage
    !! per numbered subdirectory of the execution directory, in the order of the maintainer's rulings 6 and 7
    !! (trajectory_solve3D.md, Progress): the raw frames imported and centred (center2D_nano, ten chunk averages), the
    !! one diameter estimate on the centred averages, the start model, refine2D_nano on the raw frames from the
    !! centering's shifts, the two coupled gates on its class averages (the lattice gate at LAT_PEAK_MIN_RAW, the model
    !! gate against the start model; a class is kept when both keep it), the prune, then only the kept frames denoised
    !! (extracted, flex, swapped into the project) and autorefine3D_nano. Nothing is rejected on evidence the denoiser
    !! has touched. Every stage reads the project file of the stage before and writes its own copy; the stages'
    !! records (wall time, frames in and out, numbers, done flag) are trajectory_solve3D.txt, which a rerun with
    !! mkdir=no in the same directory reads to skip the stages already done.
    !! Project: reads the stage before's file, writes the stage's own copy (the model gate in this routine: cls2D,
    !! cls3D, ptcl2D, ptcl3D through map_cavgs_selection, the prune and the 2D shifts mapped to ptcl3D)
    subroutine exec_trajectory_solve3D( self, cline )
        use single_commanders_trajectory, only: commander_trajectory_denoise, commander_trajectory_import, &
            &commander_trajectory_extract_substk, commander_trajectory_swap_stack
        use single_commanders_nano2D,     only: commander_refine2D_nano, commander_center2D_nano
        use simple_commanders_imgproc,    only: commander_estimate_diam
        use simple_commanders_cavgs,      only: commander_model_cavgs_rejection
        use simple_commanders_sim,        only: commander_simulate_nanoparticle
        use simple_cavg_quality_lattice,  only: K_MODEL, LAT_PEAK_MIN_RAW, model_gate_states
        use simple_refine3D_fnames,       only: refine3D_reproj_model_fname
        class(commander_trajectory_solve3D), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        integer,          parameter :: NSTAGES = 11
        integer,          parameter :: S_IMPORT = 1, S_CENTER2D = 2, S_DIAM = 3, S_STARTVOL = 4, S_REFINE2D = 5
        integer,          parameter :: S_GATE1 = 6, S_GATE2 = 7, S_EXTRACT = 8, S_DENOISE = 9, S_SWAP = 10
        integer,          parameter :: S_AUTOREFINE3D = 11
        character(len=*), parameter :: STAGE_PRGS(NSTAGES) = [character(len=26) :: 'trajectory_import', &
            &'center2D_nano', 'estimate_diam', 'simulate_nanoparticle', 'refine2D_nano', 'model_cavgs_rejection', &
            &'refine3D_nano', 'trajectory_extract_substk', 'trajectory_denoise', 'trajectory_swap_stack', &
            &'autorefine3D_nano']
        character(len=*), parameter :: RECORDS     = 'trajectory_solve3D.txt'
        character(len=*), parameter :: KEPT        = 'kept_frames'//STK_EXT
        character(len=*), parameter :: DENOISED    = 'trajectory_denoised'//STK_EXT
        character(len=*), parameter :: STARTVOL    = 'startvol.mrc'
        character(len=*), parameter :: GATE_TABLE  = 'model_gate_classes.tsv'
        character(len=*), parameter :: GATE_REJ    = 'model_gate_rejected_pairs.mrc'
        character(len=*), parameter :: GATE_KEPT   = 'model_gate_kept_pairs.mrc'
        character(len=*), parameter :: GATE_REPROJ = 'reprojs_startvol.mrc'
        integer,          parameter :: NCHUNKS     = 10      !< the centering's time-chunk averages
        integer,          parameter :: NKEPT_PAIRS = 20      !< kept classes shown beside the rejected ones
        real,             parameter :: LP_EST_DIAM = 3.      !< as analysis2D_nano
        real,             parameter :: MSK_FAC     = 1.5     !< mskdiam = MSK_FAC x the largest diameter
        type(commander_trajectory_import)     :: ximport
        type(commander_center2D_nano)         :: xcenter2D_nano
        type(commander_estimate_diam)         :: xest_diam
        type(commander_simulate_nanoparticle) :: xsim_np
        type(commander_refine2D_nano)         :: xrefine2D_nano
        type(commander_model_cavgs_rejection) :: xgate1
        type(commander_refine3D_nano)         :: xrefine3D_nano
        type(commander_reproject)             :: xreproject
        type(commander_trajectory_extract_substk) :: xextract
        type(commander_trajectory_denoise)    :: xdenoise
        type(commander_trajectory_swap_stack) :: xswap
        type(commander_autorefine3D_nano)     :: xautorefine3D_nano
        type(parameters) :: params
        type(oris)       :: recs
        type(string)     :: exec_dir, stk, projname, stage_dirs(NSTAGES), stage_proj(NSTAGES), map_fname
        integer(timer_int_kind) :: t0
        integer :: istage, ldim(3), nframes, box
        real    :: mskdiam, moldiam, lat_thres
        ! the program's defaults and checks, on the command line before the one params%new (the parameter class
        ! holds no program-specific checks)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('lp')      ) call cline%set('lp',      1.5)
        if( .not. cline%defined('element') ) call cline%set('element', 'Pt')
        call cline%set('oritype', 'ptcl2D')
        if( cline%get_rarg('smpd') <= 0. ) THROW_HARD('trajectory_solve3D requires a positive smpd')
        if( cline%get_rarg('lp')   <= 0. ) THROW_HARD('trajectory_solve3D requires a positive lp')
        if( cline%defined('mskdiam') )then
            if( cline%get_rarg('mskdiam') <= 0. ) THROW_HARD('trajectory_solve3D: mskdiam must be positive')
        endif
        if( cline%defined('moldiam') )then
            if( cline%get_rarg('moldiam') <= 0. ) THROW_HARD('trajectory_solve3D: moldiam must be positive')
        endif
        if( cline%defined('thres') )then
            if( cline%get_rarg('thres') <= 0. ) THROW_HARD('trajectory_solve3D: thres must be a positive log10 ratio')
        endif
        if( cline%defined('model_gate_nsd') )then
            if( cline%get_rarg('model_gate_nsd') <= 0. ) THROW_HARD('model_gate_nsd must be a positive number of robust spreads')
        endif
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call simple_getcwd(exec_dir)
        stk = simple_abspath(params%stk)
        if( .not. file_exists(stk) ) THROW_HARD('trajectory_solve3D: cannot find stk '//stk%to_char())
        call find_ldim_nptcls(stk, ldim, nframes)
        box       = ldim(1)
        projname  = basename(params%projfile)
        lat_thres = LAT_PEAK_MIN_RAW
        if( cline%defined('thres') ) lat_thres = params%thres
        do istage = 1, NSTAGES
            stage_dirs(istage) = filepath(exec_dir, int2str(istage)//'_'//trim(STAGE_PRGS(istage)))
            stage_proj(istage) = filepath(stage_dirs(istage), projname)
        enddo
        call read_records
        ! 1. import of the raw frames
        if( begin_stage(S_IMPORT, .true., params%projfile) )then
            block
                type(cmdline) :: cl
                call cl%set('prg',      'trajectory_import')
                call cl%set('mkdir',    'no')
                call cl%set('projfile', stage_proj(S_IMPORT))
                call cl%set('stk',      stk)
                call cl%set('smpd',     params%smpd)
                call ximport%execute(cl)
            end block
            call record_frames(S_IMPORT, nframes)
            call end_stage(S_IMPORT)
        endif
        ! 2. centering of the raw frames on ten time-chunk averages (in-plane alignment, shifts carried to ptcl3D)
        if( begin_stage(S_CENTER2D, .true., stage_proj(S_IMPORT)) )then
            block
                type(cmdline) :: cl
                call cl%set('prg',      'center2D_nano')
                call cl%set('mkdir',    'no')
                call cl%set('projfile', stage_proj(S_CENTER2D))
                call cl%set('ncls',     NCHUNKS)
                call cl%set('nthr',     params%nthr)
                call xcenter2D_nano%execute(cl)
            end block
            call record_shifts(S_CENTER2D, 'ptcl2D')
            call record_frames(S_CENTER2D, nframes)
            call end_stage(S_CENTER2D)
        endif
        ! 3. the one diameter estimate of the workflow, on the centred chunk averages: mskdiam and moldiam, unless given
        if( begin_stage(S_DIAM, .false.) )then
            block
                type(sp_project) :: spproj
                type(string)     :: cavgs
                integer          :: ncls
                real             :: smpd
                call spproj%read_segment('out', stage_proj(S_CENTER2D))
                call spproj%get_cavgs_stk(cavgs, ncls, smpd)
                call estimate_diameters(S_DIAM, cavgs)
                call spproj%kill
            end block
            call end_stage(S_DIAM)
        endif
        mskdiam = recs%get(S_DIAM, 'mskdiam')
        moldiam = recs%get(S_DIAM, 'moldiam')
        ! 4. the start model
        if( begin_stage(S_STARTVOL, .false.) )then
            call simulate_startvol(S_STARTVOL, moldiam)
            call recs%set(S_STARTVOL, 'moldiam', moldiam)
            call end_stage(S_STARTVOL)
        endif
        ! 5. refine2D_nano on the raw frames with the centred mask, from the centering's shifts (center=no, as
        ! analysis2D_nano runs it after centering, so that it does not accumulate shifts of its own)
        if( begin_stage(S_REFINE2D, .true., stage_proj(S_CENTER2D)) )then
            block
                type(cmdline) :: cl
                call cl%set('prg',      'refine2D_nano')
                call cl%set('mkdir',    'no')
                call cl%set('projfile', stage_proj(S_REFINE2D))
                call cl%set('mskdiam',  mskdiam)
                call cl%set('center',   'no')
                call cl%set('nthr',     params%nthr)
                call xrefine2D_nano%execute(cl)
            end block
            call record_2D(S_REFINE2D)
            call record_shifts(S_REFINE2D, 'ptcl2D')
            call record_frames(S_REFINE2D, nframes)
            call end_stage(S_REFINE2D)
        endif
        ! 6. the lattice gate on the raw class averages: class and particle states, no prune
        if( begin_stage(S_GATE1, .true., stage_proj(S_REFINE2D)) )then
            block
                type(cmdline) :: cl
                call cl%set('prg',           'model_cavgs_rejection')
                call cl%set('mkdir',         'no')
                call cl%set('projfile',      stage_proj(S_GATE1))
                call cl%set('quality_model', 'nano')
                call cl%set('mskdiam',       mskdiam)
                call cl%set('thres',         lat_thres)
                call cl%set('prune',         'no')
                call cl%set('nthr',          params%nthr)
                call xgate1%execute(cl)
            end block
            call record_gate1(S_GATE1)
            call end_stage(S_GATE1)
        endif
        ! 7. the model gate, coupled with the lattice gate: the prune, the original frame indices and the 2D shifts
        if( begin_stage(S_GATE2, .true., stage_proj(S_GATE1)) )then
            call model_gate
            call end_stage(S_GATE2)
        endif
        ! 8. the kept raw frames, extracted in project order
        if( begin_stage(S_EXTRACT, .false.) )then
            block
                type(cmdline) :: cl
                call cl%set('prg',      'trajectory_extract_substk')
                call cl%set('mkdir',    'no')
                call cl%set('projfile', stage_proj(S_GATE2))
                call cl%set('state',    1)
                call cl%set('outstk',   KEPT)
                call xextract%execute(cl)
            end block
            call require_file(S_EXTRACT, filepath(stage_dirs(S_EXTRACT), KEPT))
            call record_frames(S_EXTRACT, recs%get_int(S_GATE2, 'frames_out'))
            call end_stage(S_EXTRACT)
        endif
        ! 9. flex denoising of the kept frames only (its defaults unless the denoiser's keys are given)
        if( begin_stage(S_DENOISE, .false.) )then
            block
                type(cmdline) :: cl
                call cl%set('prg',   'trajectory_denoise')
                call cl%set('mkdir', 'no')
                call cl%set('stk',   filepath(stage_dirs(S_EXTRACT), KEPT))
                call cl%set('smpd',  params%smpd)
                call cl%set('nthr',  params%nthr)
                call cl%copy_arg(cline, 'pca_mode')
                call cl%copy_arg(cline, 'halfwin')
                call cl%copy_arg(cline, 'neigs')
                call cl%copy_arg(cline, 'eo_pool')
                call xdenoise%execute(cl)
            end block
            call require_file(S_DENOISE, filepath(stage_dirs(S_DENOISE), DENOISED))
            call record_frames(S_DENOISE, recs%get_int(S_GATE2, 'frames_out'))
            call end_stage(S_DENOISE)
        endif
        ! 10. the denoised kept frames put in the project's place
        if( begin_stage(S_SWAP, .true., stage_proj(S_GATE2)) )then
            block
                type(cmdline) :: cl
                call cl%set('prg',      'trajectory_swap_stack')
                call cl%set('mkdir',    'no')
                call cl%set('projfile', stage_proj(S_SWAP))
                call cl%set('stk',      filepath(stage_dirs(S_DENOISE), DENOISED))
                call xswap%execute(cl)
            end block
            call record_frames(S_SWAP, recs%get_int(S_GATE2, 'frames_out'))
            call end_stage(S_SWAP)
        endif
        ! 11. autorefine3D_nano on the denoised kept frames, without the atomicity prior, icm=no, the centred mask
        if( begin_stage(S_AUTOREFINE3D, .true., stage_proj(S_SWAP)) )then
            block
                type(cmdline) :: cl
                call cl%set('prg',      'autorefine3D_nano')
                call cl%set('mkdir',    'no')
                call cl%set('projfile', stage_proj(S_AUTOREFINE3D))
                call cl%set('vol1',     filepath(stage_dirs(S_STARTVOL), STARTVOL))
                call cl%set('smpd',     params%smpd)
                call cl%set('pgrp',     'c1')
                call cl%set('mskdiam',  mskdiam)
                call cl%set('lp',       params%lp)
                call cl%set('icm',      'no')
                call cl%set('nthr',     params%nthr)
                call cl%copy_arg(cline, 'maxits')
                call cl%copy_arg(cline, 'maxits_between')
                call cl%copy_arg(cline, 'nspace')
                call xautorefine3D_nano%execute(cl)
            end block
            map_fname = final_map()
            call recs%set(S_AUTOREFINE3D, 'map', map_fname)
            call record_frames(S_AUTOREFINE3D, recs%get_int(S_GATE2, 'frames_out'))
            call end_stage(S_AUTOREFINE3D)
        endif
        map_fname = recs%get_str(S_AUTOREFINE3D, 'map')
        call recs%set(NSTAGES + 1, 'stage', 'map')
        call recs%set(NSTAGES + 1, 'map',   map_fname)
        call recs%write(string(RECORDS))
        write(logfhandle,'(A)') '>>> TRAJECTORY_SOLVE3D MAP: '//map_fname%to_char()
        call recs%kill
        call simple_end('**** TRAJECTORY_SOLVE3D NORMAL STOP ****')

      contains

        !> the records of an earlier run in this directory (mkdir=no), or fresh ones
        subroutine read_records
            integer :: i
            call recs%new(NSTAGES + 1, is_ptcl=.false.)
            if( file_exists(RECORDS) )then
                call recs%read(string(RECORDS), [1, min(NSTAGES + 1, nlines(string(RECORDS)))])
                write(logfhandle,'(A)') '>>> TRAJECTORY_SOLVE3D: resuming from the records in '//RECORDS
            endif
            do i = 1, NSTAGES
                call recs%set(i, 'stage', trim(STAGE_PRGS(i)))
                if( .not. recs%isthere(i, 'done') ) call recs%set(i, 'done', 0)
            enddo
        end subroutine read_records

        !> .true. when the stage is to run: an earlier attempt's directory removed, the stage directory made (with a
        !! copy of from_proj when it reads a project) and entered, and every later stage marked as not done; .false.
        !! when its record says it is done
        logical function begin_stage( istg, l_proj, from_proj )
            integer,                 intent(in) :: istg
            logical,                 intent(in) :: l_proj
            class(string), optional, intent(in) :: from_proj
            integer :: j
            begin_stage = recs%get_int(istg, 'done') /= 1
            if( .not. begin_stage )then
                write(logfhandle,'(A)') '>>> TRAJECTORY_SOLVE3D STAGE '//int2str(istg)//' '//trim(STAGE_PRGS(istg))//&
                    &': done, skipped'
                return
            endif
            write(logfhandle,'(A)') '>>> TRAJECTORY_SOLVE3D STAGE '//int2str(istg)//' '//trim(STAGE_PRGS(istg))
            do j = istg + 1, NSTAGES
                call recs%set(j, 'done', 0)
            enddo
            call chdir_checked(exec_dir)
            call simple_rmdir_checked(stage_dirs(istg))
            call simple_mkdir(stage_dirs(istg))
            if( l_proj ) call simple_copy_file(from_proj, stage_proj(istg))
            call chdir_checked(stage_dirs(istg))
            t0 = tic()
        end function begin_stage

        !> the stage's wall time and done flag in the records, rewritten; back to the execution directory
        subroutine end_stage( istg )
            integer, intent(in) :: istg
            call chdir_checked(exec_dir)
            call recs%set(istg, 'wall_s', real(toc(t0)))
            call recs%set(istg, 'done',   1)
            call recs%write(string(RECORDS))
        end subroutine end_stage

        subroutine chdir_checked( dir )
            class(string), intent(in) :: dir
            integer :: status
            call simple_chdir(dir, status)
            if( status /= 0 ) THROW_HARD('trajectory_solve3D: cannot enter '//dir%to_char())
        end subroutine chdir_checked

        !> an unfinished stage's directory is removed before it runs again
        subroutine simple_rmdir_checked( dir )
            class(string), intent(in) :: dir
            if( dir_exists(dir) ) call exec_cmdline('rm -rf '//dir%to_char())
        end subroutine simple_rmdir_checked

        subroutine require_file( istg, fname )
            integer,       intent(in) :: istg
            class(string), intent(in) :: fname
            type(string) :: errmsg
            if( file_exists(fname) ) return
            errmsg = 'trajectory_solve3D stage '//trim(STAGE_PRGS(istg))//' did not produce '//fname%to_char()
            THROW_HARD(errmsg%to_char())
        end subroutine require_file

        !> estimate_diam on a stack of class averages (low-pass 3 A, as analysis2D_nano), taken at the 0.1 A it reports:
        !! mskdiam = MSK_FAC x the largest diameter, moldiam the smallest; given keys win
        subroutine estimate_diameters( istg, cavgs )
            integer,       intent(in) :: istg
            class(string), intent(in) :: cavgs
            type(cmdline) :: cl
            real :: dmin, dmax, mskdiam, moldiam
            call cl%set('prg',     'estimate_diam')
            call cl%set('mkdir',   'no')
            call cl%set('stk',     cavgs)
            call cl%set('smpd',    params%smpd)
            call cl%set('mskdiam', 0.)
            call cl%set('lp',      LP_EST_DIAM)
            call cl%set('nthr',    params%nthr)
            call xest_diam%execute(cl)
            dmin    = anint(10. * cl%get_rarg('min_diam')) / 10.
            dmax    = anint(10. * cl%get_rarg('max_diam')) / 10.
            mskdiam = anint(10. * MSK_FAC * dmax) / 10.
            moldiam = dmin
            if( cline%defined('mskdiam') ) mskdiam = params%mskdiam
            if( cline%defined('moldiam') ) moldiam = params%moldiam
            call recs%set(istg, 'diam_min', dmin)
            call recs%set(istg, 'diam_max', dmax)
            call recs%set(istg, 'mskdiam',  mskdiam)
            call recs%set(istg, 'moldiam',  moldiam)
            if( mskdiam > real(box) * params%smpd )then
                ! the parameter layer falls back on its box default and warns, in every stage that reads mskdiam
                call recs%set(istg, 'mask_exceeds_box', 1)
                write(logfhandle,'(A,F7.1,A)') '>>> TRAJECTORY_SOLVE3D: mskdiam ', mskdiam, &
                    &' A exceeds the box; the box default applies'
            else
                call recs%set(istg, 'mask_exceeds_box', 0)
            endif
            write(logfhandle,'(A,2F7.1)') '>>> TRAJECTORY_SOLVE3D MSKDIAM / MOLDIAM (A): ', mskdiam, moldiam
            call cl%kill
        end subroutine estimate_diameters

        subroutine simulate_startvol( istg, moldiam )
            integer, intent(in) :: istg
            real,    intent(in) :: moldiam
            type(cmdline) :: cl
            call cl%set('prg',     'simulate_nanoparticle')
            call cl%set('mkdir',   'no')
            call cl%set('outvol',  STARTVOL)
            call cl%set('smpd',    params%smpd)
            call cl%set('element', params%element)
            call cl%set('moldiam', moldiam)
            call cl%set('box',     box)
            call cl%set('nthr',    params%nthr)
            call xsim_np%execute(cl)
            call cl%kill
            call require_file(istg, filepath(stage_dirs(istg), STARTVOL))
        end subroutine simulate_startvol

        !> the 2D classes: their number with frames and the fraction of consecutive frames that share a class
        subroutine record_2D( istg )
            integer, intent(in) :: istg
            type(sp_project)     :: spproj
            integer, allocatable :: cls(:), pops(:)
            integer :: n
            call spproj%read_segment('ptcl2D', stage_proj(istg))
            call spproj%read_segment('cls2D',  stage_proj(istg))
            cls  = spproj%os_ptcl2D%get_all_asint('class')
            pops = spproj%os_cls2D%get_all_asint('pop')
            n    = size(cls)
            call recs%set(istg, 'classes', count(pops > 0))
            if( n > 1 ) call recs%set(istg, 'consec_same_class', real(count(cls(2:n) == cls(1:n-1))) / real(n - 1))
            call spproj%kill
        end subroutine record_2D

        !> gate 1: classes accepted and rejected, frames kept
        subroutine record_gate1( istg )
            integer, intent(in) :: istg
            type(sp_project)     :: spproj
            integer, allocatable :: cstates(:), pstates(:)
            call spproj%read_segment('cls2D',  stage_proj(istg))
            call spproj%read_segment('ptcl2D', stage_proj(istg))
            cstates = spproj%os_cls2D%get_all_asint('state')
            pstates = spproj%os_ptcl2D%get_all_asint('state')
            call recs%set(istg, 'cls_accepted', count(cstates > 0))
            call recs%set(istg, 'cls_rejected', count(cstates == 0))
            call recs%set(istg, 'thres',        lat_thres)
            call record_frames(istg, count(pstates > 0))
            call spproj%kill
        end subroutine record_gate1

        !> the frames a stage takes (the last earlier stage's output, the trajectory for the first) and gives on
        subroutine record_frames( istg, nout )
            integer, intent(in) :: istg, nout
            integer :: j, nin
            nin = nframes
            do j = istg - 1, 1, -1
                if( recs%isthere(j, 'frames_out') )then
                    nin = recs%get_int(j, 'frames_out')
                    exit
                endif
            enddo
            call recs%set(istg, 'frames_in',  nin)
            call recs%set(istg, 'frames_out', nout)
        end subroutine record_frames

        !> the frames' in-plane shifts of a field (pixels): median and largest magnitude
        subroutine record_shifts( istg, field )
            integer,          intent(in) :: istg
            character(len=*), intent(in) :: field
            type(sp_project)  :: spproj
            class(oris), pointer :: os
            real, allocatable :: x(:), y(:), r(:)
            call spproj%read_segment(field, stage_proj(istg))
            call spproj%ptr2oritype(field, os)
            x = os%get_all('x')
            y = os%get_all('y')
            r = sqrt(x**2 + y**2)
            if( size(r) > 0 )then
                call recs%set(istg, 'shift_median_px', median(r))
                call recs%set(istg, 'shift_max_px',    maxval(r))
            endif
            call spproj%kill
        end subroutine record_shifts

        !> gate 2 (this stage's directory): refine3D_nano oritype=cls3D maxits=1 against the start model, as
        !! autorefine3D_nano's own final pass, aligns gate 1's survivors (the classes gate 1 rejected carry state 0 in
        !! cls3D and are not searched); model_gate_states cuts them on their correlation; the states go to the classes
        !! and the particles (map_cavgs_selection) and the project is pruned. Products: the class table, the class
        !! average / reprojection pairs of the rejected classes and of NKEPT_PAIRS kept ones
        subroutine model_gate
            type(sp_project)     :: spproj
            type(cmdline)        :: cl
            type(string)         :: cavgs, errmsg
            real,    allocatable :: corrs(:), peaks(:)
            integer, allocatable :: states1(:), states2(:), pops(:), pstates(:)
            real    :: nsd, med, sd, cut, smpd, euls(3)
            integer :: ncls, icls, funit, io_stat, nframes_in, iptcl
            call cl%set('prg',      'refine3D_nano')
            call cl%set('mkdir',    'no')
            call cl%set('projfile', stage_proj(S_GATE2))
            call cl%set('vol1',     filepath(stage_dirs(S_STARTVOL), STARTVOL))
            call cl%set('smpd',     params%smpd)
            call cl%set('pgrp',     'c1')
            call cl%set('mskdiam',  mskdiam)
            call cl%set('oritype',  'cls3D')
            call cl%set('maxits',   1)
            call cl%set('lp',       1.0)
            call cl%set('trs',      10.0)
            call cl%set('icm',      'no')
            call cl%set('volrec',   'no')
            call cl%set('nthr',     params%nthr)
            call xrefine3D_nano%execute(cl)
            call cl%kill
            call del_file(refine3D_reproj_model_fname('even'))
            call del_file(refine3D_reproj_model_fname('odd'))
            call spproj%read(stage_proj(S_GATE2))
            call spproj%get_cavgs_stk(cavgs, ncls, smpd)
            if( spproj%os_cls3D%get_noris() /= ncls .or. spproj%os_cls2D%get_noris() /= ncls )then
                errmsg = 'trajectory_solve3D stage model gate: cls2D, cls3D and the class averages differ in number'
                THROW_HARD(errmsg%to_char())
            endif
            corrs   = spproj%os_cls3D%get_all('corr')
            states1 = spproj%os_cls3D%get_all_asint('state')
            pops    = spproj%os_cls2D%get_all_asint('pop')
            ! the raw lattice peak model_cavgs_rejection records on every scored class; 0 for the hard-rejected ones
            allocate(peaks(ncls), source=0.)
            do icls = 1, ncls
                if( spproj%os_cls2D%isthere(icls, 'lattice_peak') ) peaks(icls) = spproj%os_cls2D%get(icls, 'lattice_peak')
            enddo
            nsd     = K_MODEL
            if( cline%defined('model_gate_nsd') ) nsd = params%model_gate_nsd
            allocate(states2(ncls))
            call model_gate_states(corrs, states1, nsd, states2, med, sd, cut)
            if( count(states2 > 0) == 0 )then
                errmsg = 'trajectory_solve3D stage model gate: no class average left after the two gates'
                THROW_HARD(errmsg%to_char())
            endif
            ! the class table
            call fopen(funit, string(GATE_TABLE), 'replace', 'unknown', iostat=io_stat, form='formatted')
            call fileiochk('trajectory_solve3D fopen failed '//GATE_TABLE, io_stat)
            write(funit,'(A)') 'class'//char(9)//'population'//char(9)//'lattice_peak'//char(9)//'gate1_state'//char(9)//&
                &'corr'//char(9)//'e1'//char(9)//'e2'//char(9)//'e3'//char(9)//'x'//char(9)//'y'//char(9)//'gate2_state'
            do icls = 1, ncls
                euls = spproj%os_cls3D%get_euler(icls)
                write(funit,'(I0,A,I0,A,F8.4,A,I0,A,F8.4,5(A,F9.3),A,I0)') icls, char(9), pops(icls), char(9), &
                    &peaks(icls), char(9), states1(icls), char(9), corrs(icls), char(9), euls(1), char(9), euls(2), &
                    &char(9), euls(3), char(9), spproj%os_cls3D%get(icls, 'x'), char(9), spproj%os_cls3D%get(icls, 'y'), &
                    &char(9), states2(icls)
            enddo
            call fclose(funit)
            call write_gate_pairs(spproj, cavgs, ncls, states1, states2, corrs)
            ! the decision to the classes and the particles; every particle keeps its frame of the trajectory (frame,
            ! its stack index here), which the denoised stack swapped in later renumbers; then the prune and the
            ! refined 2D shifts of the kept frames to ptcl3D
            pstates    = spproj%os_ptcl2D%get_all_asint('state')
            nframes_in = count(pstates > 0)
            call spproj%map_cavgs_selection(states2)
            do iptcl = 1, spproj%os_ptcl2D%get_noris()
                call spproj%os_ptcl2D%set(iptcl, 'frame', spproj%os_ptcl2D%get_int(iptcl, 'indstk'))
                call spproj%os_ptcl3D%set(iptcl, 'frame', spproj%os_ptcl2D%get_int(iptcl, 'indstk'))
            enddo
            pstates = spproj%os_ptcl2D%get_all_asint('state')
            call recs%set(S_GATE2, 'nsd',             nsd)
            call recs%set(S_GATE2, 'corr_median',     med)
            call recs%set(S_GATE2, 'corr_spread',     sd)
            call recs%set(S_GATE2, 'cut',             cut)
            call recs%set(S_GATE2, 'cls_in',          count(states1 > 0))
            call recs%set(S_GATE2, 'cls_rejected',    count(states1 > 0 .and. states2 == 0))
            call recs%set(S_GATE2, 'frames_rejected', nframes_in - count(pstates > 0))
            call recs%set(S_GATE2, 'kept_frac',       real(count(pstates > 0)) / real(nframes))
            call record_frames(S_GATE2, count(pstates > 0))
            write(logfhandle,'(A,3F8.4)') '>>> MODEL GATE CORRELATION MEDIAN / SPREAD / CUT: ', med, sd, cut
            write(logfhandle,'(A,I0,A,I0)') '>>> MODEL GATE CLASSES REJECTED: ', count(states1 > 0 .and. states2 == 0), &
                &' OF ', count(states1 > 0)
            call spproj%prune_particles
            call spproj%map2Dshifts23D
            call spproj%write(stage_proj(S_GATE2))
            call spproj%kill
            call cavgs%kill
        end subroutine model_gate

        !> the class average / reprojection pairs for the maintainer's eye: every class gate 2 rejected, and
        !! NKEPT_PAIRS of the ones it kept spread over their correlation rank, lowest first
        subroutine write_gate_pairs( spproj, cavgs, ncls, states1, states2, corrs )
            class(sp_project), intent(inout) :: spproj
            class(string),     intent(in)    :: cavgs
            integer,           intent(in)    :: ncls, states1(:), states2(:)
            real,              intent(in)    :: corrs(:)
            type(cmdline)        :: cl
            type(image)          :: img
            integer, allocatable :: kept(:), order(:)
            integer :: icls, i, nkept, nshow, cnt
            call spproj%os_cls3D%write(string('cavgs_oris.txt'))
            call cl%set('prg',    'reproject')
            call cl%set('vol1',   filepath(stage_dirs(S_STARTVOL), STARTVOL))
            call cl%set('outstk', GATE_REPROJ)
            call cl%set('smpd',   params%smpd)
            call cl%set('oritab', 'cavgs_oris.txt')
            call cl%set('pgrp',   'c1')
            call cl%set('nthr',   params%nthr)
            call xreproject%execute(cl)
            call cl%kill
            call img%new([box,box,1], params%smpd)
            cnt = 0
            do icls = 1, ncls
                if( states1(icls) > 0 .and. states2(icls) == 0 ) call write_pair(img, cavgs, GATE_REJ, icls, cnt)
            enddo
            kept  = pack([(icls, icls = 1, ncls)], states2 > 0)
            nkept = size(kept)
            if( nkept > 0 )then
                order = kept(sort_order(corrs(kept)))
                cnt = 0
                nshow = min(NKEPT_PAIRS, nkept)
                do i = 1, nshow
                    call write_pair(img, cavgs, GATE_KEPT, order(1 + ((i - 1) * (nkept - 1)) / max(1, nshow - 1)), cnt)
                enddo
            endif
            call img%kill
            call del_file('cavgs_oris.txt')
        end subroutine write_gate_pairs

        !> one class average and its reprojection, each normalised, appended to fname
        subroutine write_pair( img, cavgs, fname, icls, cnt )
            class(image),     intent(inout) :: img
            class(string),    intent(in)    :: cavgs
            character(len=*), intent(in)    :: fname
            integer,          intent(in)    :: icls
            integer,          intent(inout) :: cnt
            call img%read(cavgs, icls)
            call img%norm
            cnt = cnt + 1
            call img%write(string(fname), cnt)
            call img%read(string(GATE_REPROJ), icls)
            call img%norm
            cnt = cnt + 1
            call img%write(string(fname), cnt)
        end subroutine write_pair

        !> ascending order of a real vector (indices into it)
        function sort_order( v ) result( inds )
            real, intent(in)     :: v(:)
            integer, allocatable :: inds(:)
            real, allocatable    :: vcopy(:)
            integer :: i
            vcopy = v
            inds  = [(i, i = 1, size(v))]
            call hpsort(vcopy, inds)
        end function sort_order

        !> the final map of autorefine3D_nano: final_results/<fbody>_iter<NNN>.mrc of its last round
        type(string) function final_map() result( fname )
            type(string), allocatable :: maps(:)
            type(string) :: fbody, pattern
            integer :: i
            fbody   = refine3D_state_vol_fbody(1)
            pattern = filepath(filepath(stage_dirs(S_AUTOREFINE3D), 'final_results'), fbody%to_char()//'_iter???'//MRC_EXT)
            call simple_list_files(pattern%to_char(), maps)
            if( .not. allocated(maps) ) THROW_HARD('trajectory_solve3D stage autorefine3D_nano: no final map')
            if( size(maps) == 0 ) THROW_HARD('trajectory_solve3D stage autorefine3D_nano: no final map')
            fname = maps(1)
            do i = 2, size(maps)
                if( maps(i)%to_char() > fname%to_char() ) fname = maps(i)
            enddo
        end function final_map

    end subroutine exec_trajectory_solve3D

end module single_commanders_trajectory_solve3D
