module simple_commander_abinitio3D_stream
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: params_glob
use simple_sp_project,         only: sp_project
use simple_qsys_env,           only: qsys_env
use simple_oris,               only: oris
use simple_stream_utils
implicit none

public :: request_snapshot, add_projects2jobslist, submit_jobs, check_processes
public :: analysis, update_user_params3D, restart_cleanup, get_nrecs_submitted

private
#include "simple_local_flags.inc"

character(len=STDLEN), parameter :: USER_PARAMS = 'stream3D_user_params.txt'
character(len=STDLEN), parameter :: EXECDIR     = 'run_'

type(procrecord),    allocatable :: procrecs(:)             ! records of all jobs
integer                          :: nvols_glob = 0          ! Total number of volumes generated (runs completed)
integer(kind=dp)                 :: time_last_submission

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
        use simple_procimgstk, only: scale_imgfile
        use simple_class_frcs
        character(len=LONGSTRLEN), allocatable, intent(in) :: projects(:)
        character(len=LONGSTRLEN), allocatable :: files(:)
        character(len=:),          allocatable :: projname, projfile, direxec, id, string
        character(len=:),          allocatable :: path, frc, cavgs, tmpl, tmp1, tmp2
        type(sp_project)          :: spproj
        type(class_frcs)          :: frcs
        character(len=LONGSTRLEN) :: project_fname
        real    :: smpd, smpd_sc
        integer :: i,j,n, jobid, box
        logical :: found
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
            direxec  = trim(EXECDIR)//id//'/'       ! relative path
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
            ! fetch frcs & classes
            path = trim(stemname(project_fname))
            tmpl = trim(basename_safe(project_fname))
            tmpl = trim(get_fbody(tmpl, METADATA_EXT, separator=.false.))
            smpd = spproj%get_smpd()
            box  = spproj%get_box()
            call spproj%get_cavgs_stk(tmp1, n, smpd_sc)
            ! call find_ldim_nptcls(string, ldim_sc, n, smpd=smpd_sc)
            ! classes, pad to original dimensions as images were scaled during alignment
            found  = .false.
            string = trim(tmpl)//'_'//trim(CAVGS_ITER_FBODY)
            call simple_list_files_regexp(path, '\.mrc$|\.mrcs$', files)
            if( allocated(files) )then
                do j = 1,size(files)
                    if( str_has_substr(files(j), string) )then
                        if( str_has_substr(files(j),'_even') .or.&
                            &str_has_substr(files(j),'_odd') )cycle
                        string = trim(files(j))
                        found = .true.
                        exit
                    endif
                enddo
                deallocate(files)
            endif
            if( .not.found )then
                THROW_WARN('Cannot find classes at: '//string//'. Skipping')
                call spproj%kill
                return
            endif
            j     = index(string,CAVGS_ITER_FBODY)
            cavgs = trim(direxec)//string(j:len_trim(string))
            call scale_imgfile(string, cavgs, smpd_sc, [box,box,1], smpd)
            call spproj%add_cavgs2os_out(cavgs, smpd, 'cavg')
            tmp1 = add2fbody(string, params_glob%ext, '_even')
            tmp2 = add2fbody(cavgs,  params_glob%ext, '_even')
            call scale_imgfile(tmp1, tmp2, smpd_sc, [box,box,1], smpd)
            tmp1 = add2fbody(string, params_glob%ext, '_odd')
            tmp2 = add2fbody(cavgs,  params_glob%ext, '_odd')
            call scale_imgfile(tmp1, tmp2, smpd_sc, [box,box,1], smpd)
            ! frcs, pad to original dimensions as images were scaled during alignment
            string = trim(path)//trim(tmpl)//'_'//trim(FRCS_FILE)
            if( .not.file_exists(string) )then
                THROW_WARN('Cannot find '//string//'. Skipping')
                call spproj%kill
                return
            endif
            frc  = direxec//trim(FRCS_FILE)
            call frcs%read(string)
            call frcs%pad(smpd, box)
            call frcs%write(frc)
            call spproj%add_frcs2os_out(frc, 'frc2D')
            ! write project
            call spproj%write(direxec//projfile)
            ! new record
            call append_procrecord(procrecs, id, direxec, projfile)
            ! tidy
            call spproj%kill
            call frcs%kill
            write(logfhandle,'(A,A)')'>>> ADDED ONE JOB: ',procrecs(get_nrecs())%folder
        enddo
    end subroutine add_projects2jobslist

    ! submits all abinitio3D jobs to the queue
    subroutine submit_jobs( cline )
        type(cmdline),     intent(in) :: cline
        integer :: i, nsubmitted
        nsubmitted = get_nrecs_submitted()
        if( nsubmitted >= params_glob%maxnruns ) return ! all resources used
        if( nsubmitted == get_nrecs() )         return ! all jobs submitted
        do i = 1,get_nrecs()
            if( procrecs(i)%submitted ) cycle
            if( get_nrecs_submitted() >= params_glob%maxnruns ) exit
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
        if(cline_abinitio3D%defined('niceserver')) call cline_abinitio3D%delete('niceserver')
        if(cline_abinitio3D%defined('niceprocid')) call cline_abinitio3D%delete('niceprocid')
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

    subroutine analysis( master_spproj )
        class(sp_project), intent(inout) :: master_spproj
        type(procrecord), allocatable :: prevrecs(:)    ! already analyzed
        type(procrecord), allocatable :: newrecs(:)     ! to analyze
        type(procrecord), allocatable :: cohort(:)      ! all
        type(sp_project)              :: spproj
        character(len=:), allocatable :: projfile
        real    :: smpd
        integer :: i, j, n_new, n_prev, ncohort, box, nptcls
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
            nptcls = spproj%os_ptcl3D%get_noris(consider_state=.true.)
            call spproj%kill
            ! volume analysis here
            ! will not go through this loop again
            procrecs(i)%included = .true.
            ! book-keeping
            j = j + 1
            newrecs(j) = procrecs(i)
            ! add to global project as state=nvols_glob
            nvols_glob = nvols_glob+1
            call master_spproj%add_vol2os_out(procrecs(i)%volume, smpd, nvols_glob, 'vol', box, pop=nptcls)
        enddo
        ! Cohort analysis
        cohort  = pack(procrecs, mask=procrecs(:)%included)
        ncohort = size(cohort)
        ! TODO
        ! write global project
        call master_spproj%write_segment_inside('out',params_glob%projfile)
        ! cleanup
        call kill_procrecords(prevrecs)
        call kill_procrecords(newrecs)
        call kill_procrecords(cohort)
    end subroutine analysis

    !> Updates current parameters with user input
    subroutine update_user_params3D( updated )
        logical, intent(out)  :: updated
        type(oris)            :: os
        ! character(len=STDLEN) :: val
        integer               :: nptcls
        updated = .false.
        call os%new(1, is_ptcl=.false.)
        if( file_exists(USER_PARAMS) )then
            call os%read(USER_PARAMS)
            if( os%isthere(1,'nptcls') )then
                nptcls = os%get_int(1,'nptcls')
                if( nptcls > 0 )then
                    params_glob%nptcls = nptcls
                    write(logfhandle,'(A,I8)')'>>> # OF PTCLS / SNAPSHOT UPDATED TO:',params_glob%nptcls
                    updated = .true.
                endif
            endif
            call del_file(USER_PARAMS)
        endif
        call os%kill
    end subroutine update_user_params3D

    subroutine restart_cleanup
        character(len=LONGSTRLEN), allocatable :: files(:)
        character(len=STDLEN),     allocatable :: folders(:)
        integer :: i
        call del_file(USER_PARAMS)
        call del_file(TERM_STREAM)
        call simple_list_files_regexp('.', '\.mrc$|\.mrcs$|\.txt$|\.star$|\.eps$|\.jpeg$|\.jpg$|\.dat$|\.bin$', files)
        if( allocated(files) )then
            do i = 1,size(files)
                call del_file(files(i))
            enddo
            deallocate(files)
        endif
        folders = simple_list_dirs('.')
        if( allocated(folders) )then
            do i = 1,size(folders)
                if( str_has_substr(folders(i),trim(EXECDIR)) ) call simple_rmdir(folders(i))
            enddo
            deallocate(folders)
        endif
    end subroutine restart_cleanup

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
        if( allocated(procrecs) )then
            get_nrecs_submitted = count(procrecs%submitted) - count(procrecs%complete)
        endif
    end function get_nrecs_submitted

    ! number of jobs to post-process/analyze
    integer function get_nrecs2postprocess()
        get_nrecs2postprocess = 0
        if( allocated(procrecs) ) get_nrecs2postprocess = count(procrecs%complete .and.(.not.procrecs%included))
    end function get_nrecs2postprocess

end module simple_commander_abinitio3D_stream
