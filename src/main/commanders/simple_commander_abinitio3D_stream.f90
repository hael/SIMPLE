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

public :: request_snapshot, add_projects2jobslist, submit_jobs, check_processes
public :: analysis, get_nrecs_submitted

private
#include "simple_local_flags.inc"

type(procrecord), allocatable :: procrecs(:)         ! records of all jobs
integer(kind=dp)              :: time_last_submission

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

    ! Append each project file to the processing queue, no submission performed
    subroutine add_projects2jobslist( projects )
        character(len=LONGSTRLEN), allocatable, intent(in) :: projects(:)
        type(sp_project)              :: spproj
        character(len=:), allocatable :: projname, projfile, direxec, id, string, path, frc
        character(len=LONGSTRLEN)     :: project_fname
        integer :: i, jobid
        if( .not.allocated(projects) )return
        do i = 1,size(projects)
            project_fname = trim(projects(i))
            if( .not.file_exists(project_fname) )then
                THROW_WARN('Could not find project: '//trim(project_fname))
                cycle
            endif
            ! prep record & folder
            jobid    = get_nrecs() + 1
            id       = int2str_pad(jobid, params_glob%numlen)
            direxec  = 'run_'//id//'/'              ! relative path
            projname = 'abinitio_'  //id
            projfile = projname//trim(METADATA_EXT)
            call simple_mkdir(direxec)
            direxec  = simple_abspath(direxec)//'/' ! absolute path & trailing /
            ! update project
            call spproj%read(project_fname)
            if( spproj%projinfo%get_noris() /= 1) THROW_HARD('Critical error add_projects2jobslist')
            call spproj%projinfo%set(1, 'projname', projname)
            call spproj%projinfo%set(1, 'projfile', projfile)
            call spproj%projinfo%set(1, 'cwd',      direxec)
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
            frc = direxec//trim(FRCS_FILE)
            call simple_copy_file(string, frc)
            call spproj%add_frcs2os_out(frc, 'frc2D') ! need to pad??
            ! write project
            call spproj%write(direxec//projfile)
            ! new record
            call append_procrecord(procrecs, id, direxec, projfile)
            ! tidy
            call spproj%kill
            write(logfhandle,'(A,A)')'>>> ADDED ONE JOB: ',procrecs(get_nrecs())%folder
            ! cleanup after copy?
        enddo
    end subroutine add_projects2jobslist

    ! submits all abinitio3D jobs to the queue
    subroutine submit_jobs( cline )
        type(cmdline),     intent(in) :: cline
        integer :: i, nsubmitted
        nsubmitted = get_nrecs_submitted()
        if( nsubmitted >= params_glob%maxjobs ) return ! all resources used
        if( nsubmitted == get_nrecs() )         return ! all jobs submitted
        do i = 1,get_nrecs()
            if( procrecs(i)%submitted ) cycle
            if( get_nrecs_submitted() >= params_glob%maxjobs ) exit
            call submit_one_process(cline, procrecs(i))
        enddo
    end subroutine submit_jobs

    ! submits one abinitio3D job to the queue
    subroutine submit_one_process( cline, record )
        type(cmdline),     intent(in)    :: cline
        class(procrecord), intent(inout) :: record
        type(chash)                   :: job_descr
        type(qsys_env)                :: qenv
        type(cmdline)                 :: cline_abinitio3D
        character(len=:), allocatable :: distrexec, log
        character(len=XLONGSTRLEN)    :: cwd, cwd_old
        ! updates command line
        cline_abinitio3D = cline
        call cline_abinitio3D%set('projfile', record%projfile)
        call cline_abinitio3D%set('prg',     'abinitio3D')
        call cline_abinitio3D%set('oritype', 'ptcl3D')
        call cline_abinitio3D%set('mkdir',   'no')
        call cline_abinitio3D%set('stream',  'yes')
        call cline_abinitio3D%delete('dir_target')
        call cline_abinitio3D%delete('walltime') !! TODO
        call cline_abinitio3D%gen_job_descr(job_descr)
        ! submit
        distrexec = 'execscript_'//record%id
        log       = 'execlog_'   //record%id
        cwd_old = trim(cwd_glob)
        call chdir(record%folder)
        call simple_getcwd(cwd)
        cwd_glob = trim(cwd)
        call qenv%new(params_glob%nparts, exec_bin='simple_exec')
        call qenv%exec_simple_prg_in_queue_async(cline_abinitio3D, distrexec, log)
        call chdir(cwd_old)
        cwd_glob = trim(cwd_old)
        time_last_submission = time8()
        ! update global record
        record%submitted = .true.
        ! tidy
        call cline_abinitio3D%kill
        call qenv%kill
        call job_descr%kill
        write(logfhandle,'(A,A)')'>>> SUBMITTED ONE JOB: ',cast_time_char(simple_gettime())
    end subroutine submit_one_process

    subroutine check_processes( ncomplete )
        integer, intent(inout) :: ncomplete
        integer :: i
        ncomplete = 0
        if( .not.allocated(procrecs) )return
        do i = 1,size(procrecs)
            if( procrecs(i)%complete ) cycle
            if( procrecs(i)%submitted )then
                procrecs(i)%complete = file_exists(trim(procrecs(i)%folder)//&
                    &trim(ABINITIO3D_FINISHED))
                if( procrecs(i)%complete )then
                    procrecs(i)%included = .false.
                    ncomplete = ncomplete + 1
                    write(logfhandle,'(A,A)')'>>> ONE JOB COMPLETED AT: ',trim(procrecs(i)%folder)
                endif
            endif
        enddo
    end subroutine check_processes

    subroutine analysis()
        type(procrecord), allocatable :: prevrecs(:)    ! already analyzed
        type(procrecord), allocatable :: newrecs(:)     ! to analyze
        type(procrecord), allocatable :: cohort(:)      ! all
        type(sp_project)              :: spproj
        character(len=:), allocatable :: projfile
        real    :: smpd
        integer :: i, j, n_new, n_prev, ncohort, box
        n_new = get_nrecs2postprocess()
        if( n_new == 0 ) return
        n_prev = get_nrecs_completed() - n_new
        if( n_prev > 0 ) prevrecs = pack(procrecs, mask=procrecs(:)%included)
        allocate(newrecs(n_new))
        ! Individual analysis
        j = 0
        do i = 1,size(procrecs)
            if( .not.procrecs(i)%complete ) cycle   ! not finished or submitted
            if(      procrecs(i)%included ) cycle   ! previously done
            projfile = procrecs(i)%folder//procrecs(i)%projfile
            ! parameters
            call spproj%read_segment('ptcl3D', projfile)
            ! parameters analysis here?
            ! unfiltered volume
            call spproj%read_segment('out', projfile)
            call spproj%get_vol('vol', 1, procrecs(i)%volume, smpd, box)
            call spproj%kill
            ! volume analysis here
            ! will not go through this loop again
            procrecs(i)%included = .true.
            ! book-keeping
            j = j + 1
            newrecs(j) = procrecs(i)
        enddo
        ! Cohort analysis
        cohort  = pack(procrecs, mask=procrecs(:)%included)
        ncohort = size(cohort)
        ! TODO
    end subroutine analysis

    ! number completed jobs
    integer function get_nrecs_completed()
        get_nrecs_completed = 0
        if( allocated(procrecs) ) get_nrecs_completed = count(procrecs%complete)
    end function get_nrecs_completed

    ! total number of jobs
    integer function get_nrecs()
        get_nrecs = 0
        if( allocated(procrecs) ) get_nrecs = size(procrecs)
    end function get_nrecs

    ! number of jobs currently submitted
    integer function get_nrecs_submitted()
        get_nrecs_submitted = 0
        if( allocated(procrecs) ) get_nrecs_submitted = count(procrecs%submitted)
    end function get_nrecs_submitted

    ! number of jobs to post-process/analyze
    integer function get_nrecs2postprocess()
        get_nrecs2postprocess = 0
        if( allocated(procrecs) ) get_nrecs2postprocess = count(procrecs%complete .and.(.not.procrecs%included))
    end function get_nrecs2postprocess

end module simple_commander_abinitio3D_stream
