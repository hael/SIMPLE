module simple_commander_abinitio3D_stream
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: params_glob
use simple_sp_project,         only: sp_project
use simple_qsys_env,           only: qsys_env
use simple_guistats,           only: guistats
use simple_oris,               only: oris
use simple_stream_utils
use simple_qsys_funs
implicit none

public :: request_snapshot, submit_one_process, check_processes

private
#include "simple_local_flags.inc"

type(procrecord), allocatable :: procrecords(:)         ! records of all jobs submitted
integer(kind=dp)              :: time_last_submission
integer                       :: njobs_glob    = 0      ! global number of abinitio3D runs currently running

contains

    ! request snapshot generation to 2D process
    ! nptcls: # of particles that will be used in next run of abinitio3D
    subroutine request_snapshot( nptcls )
        integer, intent(in) :: nptcls
        type(oris) :: os
        call os%new(1, is_ptcl=.false.)
        call os%set(1,'nptcls', nptcls)
        call os%write(trim(params_glob%dir_target)//'/'//trim(SNAPSHOT_REQUEST))
        call os%kill
        write(logfhandle,'(A,A)')'>>> REQUESTED ONE SNAPSHOT: ',cast_time_char(simple_gettime())
    end subroutine request_snapshot

    ! submits one abinitio3D job to the queue
    subroutine submit_one_process( cline, project_fname )
        type(cmdline),     intent(in) :: cline
        character(len=*),  intent(in) :: project_fname
        type(sp_project)              :: spproj
        type(chash)                   :: job_descr
        type(qsys_env)                :: qenv
        type(cmdline)                 :: cline_abinitio3D
        character(len=:), allocatable :: projfile, distrexec, log, direxec, id, string, path
        character(len=XLONGSTRLEN)    :: cwd, cwd_old
        if( .not.file_exists(project_fname) )then
            THROW_WARN('Could not find project: '//trim(project_fname))
            return
        endif
        ! updates command line
        ! TO UPDATE PARAMETERS FROM USERS HERE?
        cline_abinitio3D = cline
        call cline_abinitio3D%set('prg',     'abinitio3D')
        call cline_abinitio3D%set('oritype', 'ptcl3D')
        call cline_abinitio3D%set('mkdir',   'no')
        call cline_abinitio3D%set('stream',  'yes')
        call cline_abinitio3D%delete('dir_target')
        call cline_abinitio3D%delete('walltime')
        ! copy & rename 2D snapshot
        njobs_glob = njobs_glob + 1
        id         = int2str_pad(njobs_glob, params_glob%numlen)
        direxec    = 'run_'//id//'/'
        projfile   = 'abinitio_'  //id//trim(METADATA_EXT)
        distrexec  = 'execscript_'//id
        log        = 'execlog_'   //id
        call simple_mkdir(direxec)
        call append_procrecord(procrecords, direxec, projfile)
        call cline_abinitio3D%set('projfile', projfile)
        ! update project
        cwd_old = trim(cwd_glob)
        call chdir(direxec)
        call simple_getcwd(cwd)
        cwd_glob = trim(cwd)
        call spproj%read(project_fname)
        call spproj%update_projinfo(cline_abinitio3D)
        call spproj%os_ptcl3D%delete_2Dclustering
        ! copy frcs
        path   = trim(stemname(project_fname))
        string = trim(basename_safe(project_fname))
        string = trim(get_fbody(string, METADATA_EXT, separator=.false.))
        string = trim(path)//trim(string)//'_'//trim(FRCS_FILE)
        if( .not.file_exists(string) )then
            THROW_WARN('Cannot find '//string//'. Skipping')
            return
        endif
        call simple_copy_file(string, FRCS_FILE)
        call spproj%add_frcs2os_out(FRCS_FILE, 'frc2D')
        ! submit
        call spproj%write(projfile)
        call cline_abinitio3D%gen_job_descr(job_descr)
        call qenv%new(params_glob%nparts,exec_bin='simple_exec')
        call qenv%exec_simple_prg_in_queue_async(cline_abinitio3D, distrexec, log)
        call chdir(cwd_old)
        cwd_glob = trim(cwd_old)
        ! tidy
        time_last_submission = time8()
        call cline_abinitio3D%kill
        call qenv%kill
        write(logfhandle,'(A,A)')'>>> SUBMITTED ONE JOB: ',cast_time_char(simple_gettime())
    end subroutine submit_one_process

    subroutine check_processes( ncomplete )
        integer, intent(inout) :: ncomplete
        integer :: i
        ncomplete = 0
        if( .not.allocated(procrecords) )return
        do i = 1,size(procrecords)
            if( procrecords(i)%complete ) cycle
            procrecords(i)%complete = file_exists(trim(procrecords(i)%folder)//&
                &trim(ABINITIO3D_FINISHED))
            if( procrecords(i)%complete )then
                procrecords(i)%included = .false.
                ncomplete = ncomplete + 1
            endif
        enddo
    end subroutine check_processes

    ! number completed jobs
    integer function get_nrecords_completed()
        if( allocated(procrecords) )then
            get_nrecords_completed = count(procrecords%complete)
        else
            get_nrecords_completed = 0
        endif
    end function get_nrecords_completed

    ! total number of jobs
    integer function get_nrecords()
        if( allocated(procrecords) )then
            get_nrecords = size(procrecords)
        else
            get_nrecords = 0
        endif
    end function get_nrecords

    ! number of jobs currently submitted
    integer function get_nrunning()
        if( allocated(procrecords) )then
            get_nrunning = size(procrecords) - count(procrecords%complete) 
        else
            get_nrunning = 0
        endif
    end function get_nrunning

end module simple_commander_abinitio3D_stream
