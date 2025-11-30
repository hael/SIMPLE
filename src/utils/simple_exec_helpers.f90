module simple_exec_helpers
include 'simple_lib.f08'
use simple_qsys_env,   only: qsys_env
use simple_parameters, only: parameters, params_glob
use simple_cmdline,    only: cmdline
use simple_sp_project, only: sp_project
implicit none

public :: restarted_exec, async_exec, gen_exec_cmd, script_exec, update_job_descriptions_in_project
public :: copy_project_file_to_root_dir, set_master_num_threads, set_shmem_flag
private
#include "simple_local_flags.inc"

contains

    subroutine restarted_exec( cline, prg, executable )
        class(cmdline), intent(inout) :: cline
        class(string),  intent(in)    :: prg, executable
        type(string) :: cmd
        type(chash)  :: job_descr
        integer      :: nrestarts, i, istart
        if( .not. cline%defined('nrestarts') )then
            THROW_HARD('nrestarts needs to be defined on command line for restarted_exec')
        else
            nrestarts = cline%get_iarg('nrestarts')
            call cline%delete('nrestarts')
            istart = 1
            if( cline%defined('istart') )then
                istart = cline%get_iarg('istart')
            endif
        endif
        if( .not. cline%defined('projfile') )then
            THROW_HARD('projfile needs to be defined on command line for restarted_exec')
        endif
        call cline%set('prg', prg)
        call cline%gen_job_descr(job_descr)
        do i = istart, nrestarts
            ! compose the command line
            cmd = executable//' '//job_descr%chash2str()//' > '//uppercase(prg%to_char())//'_OUTPUT_RESTART'//int2str(i)
            ! execute
            call exec_cmdline(cmd)
        end do
        call job_descr%kill
    end subroutine restarted_exec

    subroutine async_exec( cline, executable, output )
        class(cmdline), intent(inout) :: cline
        class(string),  intent(in)    :: executable, output
        type(string) :: cmd
        type(chash)  :: job_descr
        if( .not. cline%defined('projfile') )then
            THROW_HARD('projfile needs to be defined on command line')
        endif
        call cline%gen_job_descr(job_descr)
        ! compose the command line
        cmd = string('nohup ')//executable//' '//job_descr%chash2str()//' > '//output//' &'
        ! execute asynchronously
        call exec_cmdline(cmd, waitflag=.false., suppress_errors=.true.)
        call job_descr%kill
    end subroutine async_exec

    function gen_exec_cmd( cline, executable, output ) result( cmd )
        class(cmdline),          intent(inout) :: cline
        class(string),           intent(in)    :: executable
        class(string), optional, intent(in)    :: output
        type(string)  :: cmd
        type(chash)   :: job_descr
        call cline%gen_job_descr(job_descr)
        if( present(output) )then
            cmd = executable//' '//job_descr%chash2str()//' > '//output
        else
            cmd = executable//' '//job_descr%chash2str()
        endif
        call job_descr%kill
    end function gen_exec_cmd

    subroutine script_exec( cline, prg, executable )
        class(cmdline), intent(inout) :: cline
        class(string),  intent(in)    :: prg, executable 
        type(string)     :: projfile
        type(qsys_env)   :: qenv
        type(parameters) :: params
        ! generate script for queue submission?
        if( cline%defined('script') )then
            if( cline%get_carg('script').eq.'yes' )then
                if( .not. cline%defined('projfile') ) THROW_HARD('script-based execution route requires a project file')
                projfile = cline%get_carg('projfile')
                call cline%delete('script')
                call cline%set('prg', prg)
                call cline%set('mkdir', 'no')
                call params%new(cline)
                call cline%delete('mkdir')
                call cline%delete('projfile')
                call cline%set('projfile', projfile)
                call qenv%new(1, exec_bin=executable)
                if( cline%defined('tag') )then
                    call qenv%gen_script(cline, prg//'_script'//'_'//trim(params%tag), string(uppercase(prg%to_char())//'_OUTPUT'//'_'//trim(params%tag)))
                else
                    call qenv%gen_script(cline, prg//'_script', string(uppercase(prg%to_char())//'_OUTPUT'))
                endif
                call exit
            endif
        endif
    end subroutine script_exec

    subroutine update_job_descriptions_in_project( cline )
        class(cmdline), intent(in) :: cline
        type(string)     :: exec, name
        type(chash)      :: job_descr
        type(sp_project) :: spproj
        logical          :: did_update
        if( .not. associated(params_glob)         ) return
        if( .not. associated(params_glob%ptr2prg) ) return
        exec = params_glob%ptr2prg%get_executable()
        if( exec%has_substr('private') ) return
        name = params_glob%ptr2prg%get_name()
        if( name%has_substr('print')   ) return
        if( name%has_substr('info')    ) return
        if( name%has_substr('report')  ) return
        call cline%gen_job_descr(job_descr, name)
        if( file_exists(params_glob%projfile) )then
            call spproj%read_non_data_segments(params_glob%projfile)
            call spproj%append_job_descr2jobproc(params_glob%exec_dir, job_descr, did_update)
            if( did_update ) call spproj%write_non_data_segments(params_glob%projfile)
        endif
    end subroutine update_job_descriptions_in_project

    subroutine copy_project_file_to_root_dir( cline )
        class(cmdline), intent(in)    :: cline
        type(string) :: exec
        logical :: tests(3)
        if( .not. associated(params_glob)         ) return
        if( .not. associated(params_glob%ptr2prg) ) return
        exec = params_glob%ptr2prg%get_executable()
        if( exec%has_substr('private') ) return
        tests(1) = trim(params_glob%mkdir) .eq. 'yes'
        tests(2) = params_glob%ptr2prg%requires_sp_project()
        tests(3) = file_exists(params_glob%projfile)
        if( all(tests) ) call simple_copy_file(params_glob%projfile, string('../')//params_glob%projfile)
    end subroutine copy_project_file_to_root_dir

    ! deals with # multiprocessing threads of the master process in distributed execution
    subroutine set_master_num_threads( nthr, str )
        integer,       intent(inout) :: nthr
        class(string), intent(in)    :: str
        character(len=STDLEN) :: nthr_str
        integer               :: envlen, iostat,nthr_here
        call get_environment_variable('SLURM_CPUS_PER_TASK', nthr_str, envlen)
        if( envlen > 0 )then
            nthr_here = str2int(trim(nthr_str), iostat)
        else
            nthr_here = omp_get_max_threads()
            nthr_here = min(NTHR_SHMEM_MAX, nthr_here)
        end if
        call omp_set_num_threads(nthr_here)
        write(logfhandle,'(A,A,A,I6)')'>>> # SHARED-MEMORY THREADS USED BY ', str%to_char(),' MASTER PROCESS: ', nthr_here
        nthr = nthr_here
    end subroutine set_master_num_threads

    ! deals with logical flag for shared memory execution
    logical function set_shmem_flag( cline )
        class(cmdline), intent(inout) :: cline
        if( cline%defined('nparts') )then
            set_shmem_flag = cline%get_iarg('nparts') == 1
        else
            set_shmem_flag = .true.
        endif
        if( set_shmem_flag ) call cline%delete('nparts')
    end function

end module simple_exec_helpers
