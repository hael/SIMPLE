!@descr: commander of emulate_solve3D_stream, the offline replay of the stream's solve3D and solve3D_addon cycle on an existing project
! A base solve3D on the first nptcls_base selected particles, then one solve3D_addon per chunk of nptcls_addon
! selected particles, each executed commander-style from the emulation's execution directory. The commander
! only partitions, builds the step projects, routes the keys, reads the verdicts and writes the report; the
! solve3D logic stays in commander_solve3D and commander_solve3D_addon.
! Contract: doc/implementation_notes/planned/emulate_solve3D_stream.md.
module simple_commanders_solve3D_stream_emulation
use simple_commanders_api
use simple_commanders_solve3D,       only: commander_solve3D, commander_solve3D_addon
use simple_solve3D_manifest,         only: solve3D_manifest, solve3D_stage_record, manifest_records_input
use simple_solve3D_addon_report,     only: solve3D_addon_report, ADDON_REPORT_FNAME
use simple_solve3D_stream_emulation
use simple_refine3D_fnames,          only: refine3D_fsc_fname
use simple_ui,                       only: get_prg_ptr
use simple_ui_program,               only: ui_program
implicit none

public :: commander_emulate_solve3D_stream
private
#include "simple_local_flags.inc"

integer,          parameter :: MIN_PTCLS_PER_STATE = 5       !< the stream's minimum per state (base and cohort)
character(len=*), parameter :: SOURCE_PROJ_NAME    = 'source.simple'
character(len=*), parameter :: BASE_PROJ_NAME      = 'base.simple'

type, extends(commander_base) :: commander_emulate_solve3D_stream
    contains
    procedure :: execute => exec_emulate_solve3D_stream
end type commander_emulate_solve3D_stream

contains

    subroutine exec_emulate_solve3D_stream( self, cline )
        class(commander_emulate_solve3D_stream), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        type(commander_solve3D)       :: xsolve3D
        type(commander_solve3D_addon) :: xsolve3D_addon
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(emulation_report)        :: report
        type(emulation_step)          :: step
        type(cmdline)                 :: cline_entry, cline_base, cline_add
        type(chash)                   :: job_descr
        type(ui_program), pointer     :: ui_addon => null()
        type(string),     allocatable :: keys(:)
        type(string)                  :: edir, source_in, source_fname, frozen, result_proj, step_proj, cmd_text, report_fname
        integer,          allocatable :: part(:)
        logical,          allocatable :: selected(:)
        character(len=:), allocatable :: key
        character(len=STDLEN)         :: msg
        integer(8) :: t0, t1, rate
        integer    :: i, k, status, nchunks, nbase, naddon, nmin, nselected, nfrozen_sel, nstates
        logical    :: l_rollback
        ! ---- the command line: refusals and required inputs, before anything is written
        if( .not. cline%defined('projfile') )     THROW_HARD('emulate_solve3D_stream requires projfile')
        if( .not. cline%defined('nptcls_base') )  THROW_HARD('emulate_solve3D_stream requires nptcls_base')
        if( .not. cline%defined('nptcls_addon') ) THROW_HARD('emulate_solve3D_stream requires nptcls_addon')
        keys = cline%get_keys()
        do i = 1, size(keys)
            key = trim(keys(i)%to_char())
            select case(key)
                case('vol1', 'cavg_ini', 'cavg_ini_ext', 'state', 'projfile_frozen')
                    THROW_HARD(key//' is not supported by emulate_solve3D_stream: it is not an input of the stream''s base run')
            end select
        enddo
        if( cline%defined('mkdir') )then
            if( .not. (cline%get_carg('mkdir') == 'yes') ) THROW_HARD('emulate_solve3D_stream needs mkdir=yes: its steps live in its execution directory')
        endif
        call get_prg_ptr(string('solve3D_addon'), ui_addon)
        if( .not. associated(ui_addon) ) THROW_HARD('the solve3D_addon user interface is not registered')
        ! the command line as given, kept for the steps and the report
        cline_entry = cline
        call cline_entry%set('mkdir', 'yes')
        call cline_entry%gen_job_descr(job_descr)
        cmd_text = job_descr%chash2str()
        call job_descr%kill
        source_in = simple_abspath(cline%get_carg('projfile'))
        ! ---- execution directory: the input project is copied into it
        call params%new(cline)
        call simple_getcwd(edir)
        nbase      = params%nptcls_base
        naddon     = params%nptcls_addon
        nstates    = params%nstates
        nmin       = MIN_PTCLS_PER_STATE * nstates
        l_rollback = trim(params%rollback) == 'yes'
        ! ---- the selection and its partition
        call spproj%read(params%projfile)
        if( spproj%os_ptcl2D%get_noris() < 1 ) THROW_HARD('the project has no particles (ptcl2D)')
        if( spproj%os_ptcl3D%get_noris() /= spproj%os_ptcl2D%get_noris() )then
            THROW_HARD('ptcl3D and ptcl2D of the project differ in size: the emulation needs one particle index space')
        endif
        allocate(selected(spproj%os_ptcl2D%get_noris()), source=.false.)
        do i = 1, size(selected)
            selected(i) = spproj%os_ptcl2D%get_state(i) > 0
        enddo
        nselected = count(selected)
        call partition_selected(selected, nbase, naddon, nmin, part, nchunks, status, msg)
        if( status /= 0 ) THROW_HARD(trim(msg))
        write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> EMULATE_SOLVE3D_STREAM SELECTED/BASE/ADDON CHUNK/ADDON STEPS: ', &
            &nselected, '/', nbase, '/', naddon, '/', nchunks
        ! ---- the source copy: no 3D solution, no run registration
        call prepare_emulation_source(spproj)
        source_fname = edir%to_char()//'/'//SOURCE_PROJ_NAME
        call spproj%projinfo%set(1, 'projfile', source_fname%to_char())
        call spproj%write(source_fname)
        call spproj%kill
        call report%new(cmd_text%to_char(), source_in%to_char(), nselected, nbase, naddon, nchunks, nstates, l_rollback)
        report_fname = edir%to_char()//'/'//EMULATION_REPORT_FNAME
        ! ---- the base run
        write(logfhandle,'(A)') '>>> EMULATE_SOLVE3D_STREAM BASE RUN'
        step_proj = edir%to_char()//'/'//BASE_PROJ_NAME
        call write_step_project(step_proj, 0)
        cline_base = cline_entry
        call cline_base%set('prg', 'solve3D')
        call cline_base%set('projfile', step_proj%to_char())
        call cline_base%delete('nptcls_base')
        call cline_base%delete('nptcls_addon')
        call cline_base%delete('rollback')
        call cline_base%delete('addon_diag')
        call system_clock(t0, rate)
        call xsolve3D%execute(cline_base)
        call system_clock(t1)
        call simple_chdir(edir)
        result_proj = simple_abspath(cline_base%get_carg('projfile'))
        step = emulation_step()
        step%kind     = 'base'
        step%ipart    = 0
        step%nfrozen  = 0
        step%nadded   = count(part == 0)
        step%seconds  = real(t1 - t0) / real(rate)
        step%projfile = result_proj%to_char()
        call collect_stats(result_proj, step, .false.)
        call report%add_step(step)
        call report%write(report_fname)
        frozen      = result_proj
        nfrozen_sel = step%nadded
        call cline_base%kill
        ! ---- the add-on steps
        do k = 1, nchunks
            write(logfhandle,'(A,I0,A,I0)') '>>> EMULATE_SOLVE3D_STREAM ADDON STEP ', k, ' OF ', nchunks
            step_proj = edir%to_char()//'/addon_'//int2str_pad(k, 2)//METADATA_EXT
            call write_step_project(step_proj, k)
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
            step = emulation_step()
            step%kind     = 'addon'
            step%ipart    = k
            step%nfrozen  = nfrozen_sel
            step%nadded   = count(part >= 1 .and. part <= k) + count(part == 0) - nfrozen_sel
            step%seconds  = real(t1 - t0) / real(rate)
            step%projfile = step_proj%to_char()
            call collect_stats(step_proj, step, .true.)
            if( l_rollback .and. step%l_regressed )then
                step%l_adopted = .false.
                write(logfhandle,'(A,I0,A)') '>>> EMULATE_SOLVE3D_STREAM ADDON STEP ', k, &
                    &' REGRESSED: NOT ADOPTED, THE PREVIOUS RESULT STAYS FROZEN AND THE CHUNK JOINS THE NEXT STEP'
            else
                frozen      = step_proj
                nfrozen_sel = nfrozen_sel + step%nadded
            endif
            call report%add_step(step)
            call report%write(report_fname)
        enddo
        ! ---- cleanup
        write(logfhandle,'(A,A)') '>>> EMULATE_SOLVE3D_STREAM REPORT WRITTEN: ', report_fname%to_char()
        write(logfhandle,'(A,A)') '>>> EMULATE_SOLVE3D_STREAM FINAL RESULT: ', frozen%to_char()
        call report%kill
        call cline_add%kill
        call cline_entry%kill
        call simple_touch(TASK_FINISHED)
        call simple_end('**** SIMPLE_EMULATE_SOLVE3D_STREAM NORMAL STOP ****')

    contains

        !> The project of the step that has chunks 0..kmax: the prepared source with later chunks deselected
        subroutine write_step_project( fname, kmax )
            class(string), intent(in) :: fname
            integer,       intent(in) :: kmax
            type(sp_project) :: proj
            type(string)     :: pname
            call proj%read(source_fname)
            call select_emulation_step(proj, part, kmax)
            call proj%projinfo%set(1, 'projfile', fname%to_char())
            pname = get_fbody(basename(fname), string('simple'))
            call proj%projinfo%set(1, 'projname', pname%to_char())
            call proj%write(fname)
            call proj%kill
        end subroutine write_step_project

        !> The add-on command line takes the keys the solve3D_addon interface accepts and the base run's manifest
        !! does not record (those reach the add-on through the manifest); the emulation's own keys stay behind
        subroutine forward_keys( cline_dst )
            class(cmdline), intent(inout) :: cline_dst
            type(string), allocatable :: ekeys(:)
            character(len=:), allocatable :: ekey
            integer :: j
            ekeys = cline_entry%get_keys()
            do j = 1, size(ekeys)
                ekey = trim(ekeys(j)%to_char())
                select case(ekey)
                    case('prg', 'projfile', 'projfile_frozen', 'mkdir', 'nptcls_base', 'nptcls_addon', 'rollback')
                        cycle
                end select
                if( .not. ui_addon%accepts(ekey) ) cycle
                if( manifest_records_input(ekey) ) cycle
                call cline_dst%copy_arg(cline_entry, ekey)
            enddo
        end subroutine forward_keys

        !> The last stage and limits from the run manifest, the per-state resolutions from the registered FSCs
        !! and, for an add-on, the verdicts from the add-on report in its run directory
        subroutine collect_stats( projfile, stp, l_addon )
            class(string),        intent(in)    :: projfile
            type(emulation_step), intent(inout) :: stp
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
                        stp%verdict(s)     = addon_report%get_verdict(s)
                        stp%dshell(s)      = addon_report%get_dshell(s)
                        stp%corr(s)        = addon_report%get_corr(s)
                        stp%res_cohort(s)  = addon_report%get_cohort_res0143(s)
                    enddo
                    stp%l_regressed = addon_report%any_regressed()
                    call addon_report%kill
                else
                    THROW_WARN('the add-on run wrote no report: '//rep_fname%to_char())
                endif
            endif
            call man%kill
            call proj%kill
        end subroutine collect_stats

    end subroutine exec_emulate_solve3D_stream

end module simple_commanders_solve3D_stream_emulation
