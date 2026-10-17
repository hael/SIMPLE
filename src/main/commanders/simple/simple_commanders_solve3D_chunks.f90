!@descr: solve3D_chunks commander, which runs solve3D followed by solve3D_addon over balanced particle chunks
! The runs are consecutive whole stacks with balanced sizes; the step projects, forwarded keys and step
! statistics are leaf utilities of this module. Contract: doc/policies/3D/solve3D_chunks_policy.md.
module simple_commanders_solve3D_chunks
use simple_commanders_api
use simple_commanders_solve3D,       only: commander_solve3D, commander_solve3D_addon
use simple_solve3D_manifest,         only: solve3D_manifest, solve3D_stage_record, manifest_records_input
use simple_solve3D_addon_report,     only: solve3D_addon_report, ADDON_REPORT_FNAME
use simple_solve3D_chunks_plan, only: chunks_report, chunks_step, &
    &partition_selected_balanced, chunk_nsample, prepare_chunks_source, select_chunks_step, &
    &CHUNKS_REPORT_FNAME
use simple_refine3D_fnames,          only: refine3D_fsc_fname
use simple_ui,                       only: get_prg_ptr
use simple_ui_program,               only: ui_program
implicit none

public :: commander_solve3D_chunks
private
#include "simple_local_flags.inc"

character(len=*), parameter :: SOURCE_PROJ_NAME    = 'source.simple'
character(len=*), parameter :: BASE_PROJ_NAME      = 'base.simple'

type, extends(commander_base) :: commander_solve3D_chunks
    contains
    procedure :: execute => exec_solve3D_chunks
end type commander_solve3D_chunks

contains

    subroutine exec_solve3D_chunks( self, cline )
        use simple_solve3D_utils, only: solve3D_symsrch_stage
        class(commander_solve3D_chunks), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(commander_solve3D)       :: xsolve3D
        type(commander_solve3D_addon) :: xsolve3D_addon
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(chunks_report)        :: report
        type(chunks_step)          :: step
        type(cmdline)                 :: cline_entry, cline_base, cline_add
        type(chash)                   :: job_descr
        type(ui_program), pointer     :: ui_addon => null(), ui_workflow => null()
        type(string),     allocatable :: keys(:)
        type(string)                  :: edir, source_in, source_fname, frozen, result_proj, step_proj, cmd_text, report_fname
        type(string)                  :: pgrp, pgrp_start
        integer,          allocatable :: part(:), part_sizes(:), stkind(:)
        logical,          allocatable :: selected(:)
        character(len=:), allocatable :: key
        character(len=STDLEN)         :: msg
        integer(8) :: t0, t1, rate
        integer    :: i, k, s, status, nchunks, nselected, nfrozen_sel, nstates, nsample, ninside, min_stage
        ! ---- validate the public command line before setting the internal sampling target or writing anything
        call get_prg_ptr(string('solve3D_chunks'), ui_workflow)
        if( .not. associated(ui_workflow) ) THROW_HARD('the solve3D_chunks user interface is not registered')
        if( .not. cline%defined('projfile') ) THROW_HARD('solve3D_chunks requires projfile')
        keys = cline%get_keys()
        do i = 1, size(keys)
            key = trim(keys(i)%to_char())
            if( key == 'prg' .or. key == 'projfile' .or. key == 'mkdir' ) cycle
            if( ui_workflow%accepts(key) ) cycle
            msg = key//' is not an input of solve3D_chunks'
            THROW_HARD(trim(msg))
        enddo
        if( cline%defined('mkdir') )then
            if( .not. (cline%get_carg('mkdir') == 'yes') ) THROW_HARD('solve3D_chunks needs mkdir=yes: its steps live in its execution directory')
        endif
        ! solve3D samples nsample*nstates particles per iteration; nsample keeps that total fixed for any nstates
        nstates = 1
        if( cline%defined('nstates') ) nstates = cline%get_iarg('nstates')
        if( nstates < 1 ) THROW_HARD('solve3D_chunks needs nstates >= 1')
        ! a base run that stops before the symmetry-search stage registers no manifest (one state)
        ! or cannot be continued by solve3D_addon: refuse it before the base run, not after
        if( cline%defined('nstages') )then
            pgrp       = 'c1'
            pgrp_start = 'c1'
            if( cline%defined('pgrp') )       pgrp       = cline%get_carg('pgrp')
            if( cline%defined('pgrp_start') ) pgrp_start = cline%get_carg('pgrp_start')
            min_stage = solve3D_symsrch_stage(lowercase(pgrp%to_char()) /= lowercase(pgrp_start%to_char()))
            if( cline%get_iarg('nstages') < min_stage )then
                write(msg,'(A,I0,A,I0)') 'solve3D_chunks needs nstages >= ', min_stage, '; got ', &
                    &cline%get_iarg('nstages')
                THROW_HARD(trim(msg))
            endif
        endif
        nsample = chunk_nsample(nstates)
        call cline%set('nsample', nsample)
        ! workflow defaults: cohort-only diagnostics and NU filtering with LP-set matching; gridding
        ! reconstruction always (the program offers no PCG controls); balance follows solve3D's default
        if( .not. cline%defined('addon_diag') ) call cline%set('addon_diag', 'yes')
        if( .not. cline%defined('filt_mode')  ) call cline%set('filt_mode',  'nonuniform_lpset')
        call cline%set('rec_backend', 'gridding')
        call cline%set('mkdir', 'yes')
        call get_prg_ptr(string('solve3D_addon'), ui_addon)
        if( .not. associated(ui_addon) ) THROW_HARD('the solve3D_addon user interface is not registered')
        ! keep the effective command line for the nested runs and report
        cline_entry = cline
        call cline_entry%set('mkdir', 'yes')
        call cline_entry%gen_job_descr(job_descr)
        cmd_text = job_descr%chash2str()
        call job_descr%kill
        source_in = simple_abspath(cline%get_carg('projfile'))
        ! ---- execution directory: the input project is copied into it
        call params%new(cline)
        call simple_getcwd(edir)
        if( simple_abspath(params%projfile) == source_in )then
            THROW_HARD('solve3D_chunks did not enter its own execution directory; its steps would overwrite the input project')
        endif
        ! ---- split active particles into balanced capped runs
        call spproj%read(params%projfile)
        if( spproj%os_ptcl2D%get_noris() < 1 ) THROW_HARD('the project has no particles (ptcl2D)')
        if( spproj%os_ptcl3D%get_noris() /= spproj%os_ptcl2D%get_noris() )then
            THROW_HARD('ptcl3D and ptcl2D of the project differ in size: solve3D_chunks needs one particle index space')
        endif
        allocate(selected(spproj%os_ptcl2D%get_noris()), source=.false.)
        do i = 1, size(selected)
            selected(i) = spproj%os_ptcl2D%get_state(i) > 0
        enddo
        nselected = count(selected)
        ! the runs are consecutive whole stacks in arrival order, as the stream forms them
        allocate(stkind(size(selected)), source=0)
        do i = 1, size(selected)
            stkind(i) = spproj%os_ptcl2D%get_int(i, 'stkind')
        enddo
        call partition_selected_balanced(selected, stkind, part, nchunks, ninside, status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        allocate(part_sizes(0:nchunks))
        do k = 0, nchunks
            part_sizes(k) = count(part == k)
        enddo
        write(logfhandle,'(A,I0,A,I0)') '>>> SOLVE3D_CHUNKS SELECTED/RUNS: ', nselected, '/', nchunks + 1
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> SOLVE3D_CHUNKS NSAMPLE/STATES/SAMPLED PER ITERATION: ', nsample, '/', &
            &nstates, '/', nsample * nstates
        write(logfhandle,'(A,I0)') '>>> SOLVE3D_CHUNKS RUNS: WHOLE STACKS IN ARRIVAL ORDER; CUTS INSIDE A STACK: ', ninside
        do k = 0, nchunks
            write(logfhandle,'(A,I0,A,I0)') '>>> SOLVE3D_CHUNKS RUN ', k, ' PARTICLES: ', part_sizes(k)
        enddo
        ! ---- the source copy: no 3D solution, no run registration
        call prepare_chunks_source(spproj)
        source_fname = edir%to_char()//'/'//SOURCE_PROJ_NAME
        call spproj%projinfo%set(1, 'projfile', source_fname%to_char())
        call spproj%write(source_fname)
        call spproj%kill
        call report%new_chunks(cmd_text%to_char(), source_in%to_char(), part_sizes, nstates, nsample, ninside)
        report_fname = edir%to_char()//'/'//CHUNKS_REPORT_FNAME
        ! ---- the base run
        write(logfhandle,'(A)') '>>> SOLVE3D_CHUNKS BASE RUN'
        step_proj = edir%to_char()//'/'//BASE_PROJ_NAME
        call write_solve3D_step_project(source_fname, step_proj, part, 0)
        cline_base = cline_entry
        call cline_base%set('prg', 'solve3D')
        call cline_base%set('projfile', step_proj%to_char())
        call cline_base%delete('addon_diag')
        call system_clock(t0, rate)
        call xsolve3D%execute(cline_base)
        call system_clock(t1)
        call simple_chdir(edir)
        result_proj = simple_abspath(cline_base%get_carg('projfile'))
        step = chunks_step()
        step%kind     = 'base'
        step%ipart    = 0
        step%nfrozen  = 0
        step%nadded   = part_sizes(0)
        step%seconds  = real(t1 - t0) / real(rate)
        step%projfile = result_proj%to_char()
        call collect_solve3D_step_stats(result_proj, step, .false.)
        call report%add_step(step)
        call report%write(report_fname)
        frozen      = result_proj
        nfrozen_sel = step%nadded
        call cline_base%kill
        ! ---- every add-on is exactly one planned run and is adopted unconditionally
        do k = 1, nchunks
            write(logfhandle,'(A,I0,A,I0)') '>>> SOLVE3D_CHUNKS ADDON STEP ', k, ' OF ', nchunks
            step_proj = edir%to_char()//'/addon_'//int2str_pad(k, 2)//METADATA_EXT
            call write_solve3D_step_project(source_fname, step_proj, part, k)
            call cline_add%kill
            call cline_add%set('prg',             'solve3D_addon')
            call cline_add%set('projfile',        step_proj%to_char())
            call cline_add%set('projfile_frozen', frozen%to_char())
            call cline_add%set('mkdir',           'yes')
            call forward_keys(cline_add)
            call system_clock(t0, rate)
            call xsolve3D_addon%execute(cline_add)
            call system_clock(t1)
            call simple_chdir(edir)
            step = chunks_step()
            step%kind     = 'addon'
            step%ipart    = k
            step%nfrozen  = nfrozen_sel
            step%nadded   = part_sizes(k)
            step%seconds  = real(t1 - t0) / real(rate)
            step%projfile = step_proj%to_char()
            call collect_solve3D_step_stats(step_proj, step, .true.)
            if( step%l_regressed )then
                msg = ''
                do s = 1, size(step%verdict)
                    if( trim(step%verdict(s)) == 'REGRESSED' ) msg = trim(msg)//' '//int2str(s)
                enddo
                write(logfhandle,'(A,I0,A,A,A)') '>>> SOLVE3D_CHUNKS ADDON STEP ', k, ' REGRESSED IN STATE(S):', &
                    &trim(msg), '; ADOPTED ALL THE SAME'
            endif
            frozen      = step_proj
            nfrozen_sel = nfrozen_sel + step%nadded
            call report%add_step(step)
            call report%write(report_fname)
        enddo
        ! ---- publish the final result over the execution directory's project, where simple_exec records the job
        call publish_solve3D_result(frozen, params%projfile)
        write(logfhandle,'(A,A)') '>>> SOLVE3D_CHUNKS REPORT WRITTEN: ', report_fname%to_char()
        write(logfhandle,'(A,A)') '>>> SOLVE3D_CHUNKS FINAL RESULT: ', frozen%to_char()
        write(logfhandle,'(A,A)') '>>> SOLVE3D_CHUNKS PUBLISHED THE FINAL RESULT TO ', params%projfile%to_char()
        call report%kill
        call cline_add%kill
        call cline_entry%kill
        call simple_touch(TASK_FINISHED)
        call simple_end('**** SIMPLE_SOLVE3D_CHUNKS NORMAL STOP ****')

    contains

        !> Forward accepted, non-manifest inputs; the internal nsample is replayed from the base manifest.
        subroutine forward_keys( cline_dst )
            class(cmdline), intent(inout) :: cline_dst
            type(string), allocatable :: ekeys(:)
            character(len=:), allocatable :: ekey
            integer :: j
            ekeys = cline_entry%get_keys()
            do j = 1, size(ekeys)
                ekey = trim(ekeys(j)%to_char())
                select case(ekey)
                    case('prg', 'projfile', 'projfile_frozen', 'mkdir', 'nsample')
                        cycle
                end select
                call copy_forwardable_addon_arg(cline_dst, cline_entry, ui_addon, ekey)
            enddo
        end subroutine forward_keys

    end subroutine exec_solve3D_chunks

    !> Write one cumulative step project from the prepared source.
    subroutine write_solve3D_step_project( source_fname, fname, part, kmax )
        class(string), intent(in) :: source_fname, fname
        integer,       intent(in) :: part(:), kmax
        type(sp_project) :: proj
        type(string)     :: pname
        call proj%read(source_fname)
        call select_chunks_step(proj, part, kmax)
        call proj%projinfo%set(1, 'projfile', fname%to_char())
        pname = get_fbody(basename(fname), string('simple'))
        call proj%projinfo%set(1, 'projname', pname%to_char())
        call proj%write(fname)
        call proj%kill
    end subroutine write_solve3D_step_project

    !> Publish a finished project over dst, as solve3D_addon publishes its result: written beside dst
    !! under a temporary name, then renamed over it. The run manifest is registered by absolute path,
    !! since a bare name resolves against the project file's own directory.
    subroutine publish_solve3D_result( src, dst )
        class(string), intent(in) :: src, dst
        type(sp_project)       :: proj
        type(solve3D_manifest) :: man
        type(string)           :: manifest_abs, pname, tmpname
        character(len=STDLEN)  :: emsg
        integer :: stat
        call proj%read(src)
        call man%read_registered(proj, src, stat, emsg)
        if( stat /= 0 ) THROW_HARD('the final result registers no solve3D manifest: '//trim(emsg))
        manifest_abs = man%get_fname()
        call man%register(proj, manifest_abs%to_char())
        call proj%projinfo%set(1, 'projfile', dst%to_char())
        pname = get_fbody(basename(dst), string('simple'))
        call proj%projinfo%set(1, 'projname', pname%to_char())
        tmpname = stemname(dst)//'/solve3D_chunks_publish_'//basename(dst)
        call proj%write(tmpname)
        call simple_rename(tmpname, dst, overwrite=.true.)
        call man%kill
        call proj%kill
    end subroutine publish_solve3D_result

    !> Copy one add-on input only when the add-on accepts it and the base manifest does not record it.
    subroutine copy_forwardable_addon_arg( cline_dst, cline_src, ui_addon, key )
        class(cmdline),    intent(inout) :: cline_dst
        class(cmdline),    intent(in)    :: cline_src
        class(ui_program), intent(in)    :: ui_addon
        character(len=*),  intent(in)    :: key
        if( .not. ui_addon%accepts(key) ) return
        if( manifest_records_input(key) ) return
        call cline_dst%copy_arg(cline_src, key)
    end subroutine copy_forwardable_addon_arg

    !> Read the last stage and limits from the manifest, result FSCs and, for an add-on, its comparison report.
    subroutine collect_solve3D_step_stats( projfile, stp, l_addon )
        class(string),        intent(in)    :: projfile
        type(chunks_step), intent(inout) :: stp
        logical,              intent(in)    :: l_addon
        type(sp_project)           :: proj
        type(solve3D_manifest)     :: man
        type(solve3D_stage_record) :: stage
        type(solve3D_addon_report) :: addon_report
        type(string)               :: run_dir, fsc_name, fsc_path, rep_fname
        real, allocatable :: fsc(:), res(:)
        real    :: fsc05, fsc0143
        integer :: s, ns, box_fsc, stat
        character(len=STDLEN) :: emsg
        call proj%read(projfile)
        call man%read_registered(proj, projfile, stat, emsg)
        if( stat /= 0 ) THROW_HARD('the step result registers no solve3D manifest ('//trim(emsg)//'); a base run needs nstages >= 3')
        run_dir = stemname(man%get_fname())
        ns      = man%get_nstates()
        stp%last_stage = man%get_last_stage()
        if( stp%last_stage >= 1 .and. stp%last_stage <= man%get_nstages() )then
            stage        = man%get_stage(stp%last_stage)
            stp%lp       = stage%lp_emitted
            stp%box_crop = stage%box_crop
        endif
        allocate(stp%res0143(ns), stp%res05(ns), stp%corr(ns), stp%res_cohort(ns), source=0.)
        allocate(stp%dshell(ns), source=0)
        allocate(stp%verdict(ns))
        stp%verdict = 'NOT_COMPARED'
        do s = 1, ns
            fsc_path = ''
            if( proj%isthere_in_osout('fsc', s) )then
                call proj%get_fsc(s, fsc_name, box_fsc)
                fsc_path = fsc_name
                if( .not. file_exists(fsc_path) ) fsc_path = stemname(projfile)//'/'//fsc_name%to_char()
            else
                box_fsc  = man%get_box()
                fsc_name = refine3D_fsc_fname(s)
                fsc_path = run_dir%to_char()//'/'//fsc_name%to_char()
            endif
            if( file_exists(fsc_path) )then
                fsc = file2rarr(fsc_path)
                res = get_resarr(box_fsc, proj%get_smpd())
                if( size(fsc) == size(res) )then
                    call get_resolution(fsc, res, fsc05, fsc0143)
                    stp%res0143(s) = fsc0143
                    stp%res05(s)   = fsc05
                endif
            endif
        enddo
        if( l_addon )then
            rep_fname = run_dir%to_char()//'/'//ADDON_REPORT_FNAME
            if( file_exists(rep_fname) )then
                call addon_report%read(rep_fname)
                do s = 1, min(ns, addon_report%get_nstates())
                    stp%verdict(s)    = addon_report%get_verdict(s)
                    stp%dshell(s)     = addon_report%get_dshell(s)
                    stp%corr(s)       = addon_report%get_corr(s)
                    stp%res_cohort(s) = addon_report%get_cohort_res0143(s)
                enddo
                stp%l_regressed = addon_report%any_regressed()
                call addon_report%kill
            else
                THROW_WARN('the add-on run wrote no report: '//rep_fname%to_char())
            endif
        endif
        call man%kill
        call proj%kill
    end subroutine collect_solve3D_step_stats

end module simple_commanders_solve3D_chunks
