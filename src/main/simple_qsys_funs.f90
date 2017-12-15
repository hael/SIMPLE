! batch-processing manager - functions
module simple_qsys_funs
#include "simple_lib.f08"
implicit none

interface qsys_watcher
    module procedure qsys_watcher_1
    module procedure qsys_watcher_2
end interface

integer, parameter :: SHORTTIME = 3

contains

    subroutine qsys_cleanup( p )
        use simple_params, only: params
        class(params), intent(in) :: p
        ! character(len=:), allocatable :: rec_base_str, rho_base_str, rec_base_str_part, rho_base_str_part
        integer, parameter :: NUMLEN_STATE = 2, NUMLEN_ITER = 3
        integer :: istate, iter
        ! single files
        call del_file('FOO')
        call del_file('fort.0')
        call del_file('ftab_from_sys_find_last_fname.txt')
        call del_file('VOLASSEMBLE')
        call del_file('CAVGASSEMBLE')
        call del_file('SYMSRCH')
        ! state numbered files
        call del_files('VOLASSEMBLE_FINISHED_STATE', p%nstates, numlen=NUMLEN_STATE)
        call del_files('simple_script_state',        p%nstates, numlen=NUMLEN_STATE)
        ! part numbered files
        call del_files('OUT',                        p%nparts)
        call del_files('algndoc_',                   p%nparts, ext=METADATEXT)
        call del_files('unidoc_',                    p%nparts, ext='.txt')
        call del_files('cavgs_even_part',            p%nparts, ext=p%ext )
        call del_files('cavgs_odd_part',             p%nparts, ext=p%ext )
        call del_files('JOB_FINISHED_',              p%nparts)
        call del_files('ctfsqsums_even_part',        p%nparts, ext=p%ext )
        call del_files('ctfsqsums_odd_part',         p%nparts, ext=p%ext )
        call del_files('distr_simple_script_',       p%nparts)
        call del_files('ctffind_output_part',        p%nparts, ext='.txt')
        call del_files('ctffind_ctrl_file_part',     p%nparts, ext='.txt')
        ! state and part numbered files
        do istate=1,MAXS
            ! allocate(rec_base_str, source='recvol_state'//int2str_pad(istate,NUMLEN_STATE))
            ! allocate(rho_base_str, source='rho_'//rec_base_str)
            ! allocate(rec_base_str_part, source=rec_base_str//'_part')
            ! allocate(rho_base_str_part, source='rho_'//rec_base_str_part)
            ! call del_files(rec_base_str_part, p%nparts, ext=p%ext)
            ! call del_files(rho_base_str_part, p%nparts, ext=p%ext)
            ! call del_files(rec_base_str_part, p%nparts, ext=p%ext, suffix='_even')
            ! call del_files(rec_base_str_part, p%nparts, ext=p%ext, suffix='_odd')
            ! call del_files(rho_base_str_part, p%nparts, ext=p%ext, suffix='_even')
            ! call del_files(rho_base_str_part, p%nparts, ext=p%ext, suffix='_odd')
            ! call del_file(rho_base_str//'_even'//p%ext)
            ! call del_file(rho_base_str//'_odd'//p%ext)
            ! deallocate(rec_base_str,rho_base_str,rec_base_str_part,rho_base_str_part)
        end do
        ! syscalls to del files
        ! call sys_del_files('VOLASSEMBLE_STATE', '')
    end subroutine qsys_cleanup

    subroutine terminate_if_prg_in_cwd( prg )
        character(len=*), intent(in) :: prg
        integer                      :: file_stat, funit, ios, file_rec, i, pos
        character(len=STDLEN)        :: pid, pid_dir, pwd, cmd
        cmd = 'pgrep '//trim(prg)//' > pid_from_fortran.txt'
        call exec_cmdline(cmd)
        ! read in the pid
        call fopen(funit, status='old', action='read', file='pid_from_fortran.txt', iostat=file_stat)
        call fileio_errmsg("qsys_funs; terminate_if_prg_in_cwd fopen failed", file_stat)
        read(funit,'(a)') pid
        call fclose(funit,errmsg="qsys_funs; terminate_if_prg_in_cwd fclose failed")
        ! get a string that contains the directory it is running in with lsof -p
        cmd = 'lsof -p '//trim(pid)//' | grep cwd > pid_dir_from_fortran.txt'
        call exec_cmdline(cmd)
        call fopen(funit, status='old', action='read', file='pid_dir_from_fortran.txt', iostat=file_stat)
        call fileio_errmsg("qsys_funs; terminate_if_prg_in_cwd fopen failed pid_dir_from_fortran.txt", file_stat)
        read(funit,'(a)', iostat=file_stat) pid_dir
        if(file_stat/=0)call fileio_errmsg("qsys_funs; terminate_if_prg_in_cwd read failed", file_stat)
        call fclose(funit,errmsg="qsys_funs; terminate_if_prg_in_cwd fclose failed")
        ! get pwd
        pwd = simple_getenv('PWD')
        ! assert
        pos = index(pid_dir, '/') ! start of directory
        if( trim(pid_dir(pos:)) .eq. trim(pwd) )then
            write(*,*) 'Another process of your program: ', trim(prg)
            write(*,*) 'is already running in the cwd: ', trim(pid_dir(pos:)), ' ABORTING...'
            stop
        endif
        ! remove temporary files
        call del_file('pid_from_fortran.txt')
        call del_file('pid_dir_from_fortran.txt')
    end subroutine terminate_if_prg_in_cwd

    subroutine parse_env_file( env )
        use simple_chash, only: chash
        type(chash), intent(inout)    :: env
        character(len=STDLEN)         :: env_file, qsys_name
        character(len=:), allocatable :: simple_path, simple_path_env
        integer :: i,funit, nl, file_stat
        env_file = './simple_distr_config.env'
        if( .not. file_exists(trim(env_file)) ) call autogen_env_file(env_file)
        call env%read( trim(env_file) )
        ! User keys
        if( .not. env%isthere('simple_path') )stop 'Path to SIMPLE directory is required in simple_distr_config.env (simple_path)'
        simple_path = simple_getenv('SIMPLE_PATH')
        if( allocated(simple_path) )then
            simple_path_env = env%get('simple_path')
            if( str_has_substr(simple_path, simple_path_env) .or. str_has_substr(simple_path_env, simple_path) )then
                ! all ok
            else
                write(*,*) 'WARNING! simple absolute paths in shell and simple_distr_config.env do not agree'
                write(*,*) 'In env:   ', trim(simple_path_env)
                write(*,*) 'In shell: ', trim(simple_path)
            endif
        endif
        ! Queue keys
        if( .not. env%isthere('qsys_name') ) &
            & stop 'The type of the queuing system is required in simple_distr_config.env: qsys_name=<local|slurm|pbs|sge>'
        qsys_name = env%get('qsys_name')
        if( qsys_name.ne.'local' .and. qsys_name.ne.'pbs' .and. qsys_name.ne.'slurm' .and. qsys_name.ne.'sge' ) &
            & stop 'Invalid qsys_name in simple_distr_config.env'
        if( qsys_name.eq.'slurm' )then
            if( .not. env%isthere('qsys_partition') )stop 'qsys_partition field is required in simple_distr_config.env'
        endif
        ! Job keys
        if( qsys_name.ne.'local' )then
            if( .not. env%isthere('job_memory_per_task') )stop 'Job memory is required in simple_distr_config.env (job_memory)'
            if( .not. env%isthere('job_name') )             call env%push('job_name','simple_job')
            if( .not. env%isthere('job_time') )             call env%push('job_time', '0-2:0:0')
            if( .not. env%isthere('job_ntasks') )           call env%push('job_ntasks', '1')
            if( .not. env%isthere('job_cpus_per_task') )    call env%push('job_cpus_per_task', '1')
            if( .not. env%isthere('job_ntasks_per_socket') )call env%push('job_ntasks_per_socket', '1')
        endif
    end subroutine parse_env_file

    subroutine autogen_env_file( env_file )
        use simple_syslib, only: simple_getenv
        character(len=*), intent(in)  :: env_file
        integer                       :: funit, file_stat
        character(len=:), allocatable :: simple_path, simple_email, simple_qsys
        simple_path = simple_getenv('SIMPLE_PATH')
        if( .not. allocated(simple_path) )&
        stop 'need SIMPLE_PATH env var to auto-generate simple_distr_config.env; simple_qsys_funs :: autogen_env_file'
        simple_qsys = simple_getenv('SIMPLE_QSYS')
        if( .not. allocated(simple_qsys) )&
        stop 'need SIMPLE_QSYS env var to auto-generate simple_distr_config.env; simple_qsys_funs :: autogen_env_file'
        simple_email = simple_getenv('SIMPLE_EMAIL')
        if( .not. allocated(simple_email) ) allocate(simple_email, source='me.myself\uni.edu')

        call fopen(funit, status='replace', action='write', file=trim(env_file), iostat=file_stat)
        call fileio_errmsg("autogen_env_file qsys_funs ", file_stat)
        write(funit,'(a)') '# CONFIGURATION FILE FOR DISTRIBUTED SIMPLE EXECUTION'
        write(funit,'(a)') ''
        write(funit,'(a)') '# ABSOLUTE PATH TO SIMPLE ROOT DIRECTORY'
        write(funit,'(a)') 'simple_path           = '//simple_path
        write(funit,'(a)') ''
        write(funit,'(a)') '# ESTIMATED TIME PER IMAGE (IN SECONDS)'
        write(funit,'(a)')  'time_per_image        = 400'
        write(funit,'(a)') ''
        write(funit,'(a)') '# USER DETAILS'
        write(funit,'(a)') 'user_account          ='
        write(funit,'(a)') 'user_email            = '//simple_email
        write(funit,'(a)') 'user_project          ='
        write(funit,'(a)') ''
        write(funit,'(a)') '# QSYS DETAILS (qsys_name=<local|slurm|pbs>)'
        write(funit,'(a)') 'qsys_name             = '//simple_qsys
        write(funit,'(a)') 'qsys_partition        ='
        write(funit,'(a)') 'qsys_qos              ='
        write(funit,'(a)') 'qsys_reservation      ='
        write(funit,'(a)') ''
        write(funit,'(a)') '# JOB DETAILS'
        write(funit,'(a)') 'job_ntasks            = 1'
        write(funit,'(a)') 'job_memory_per_task   = 32000'
        write(funit,'(a)') 'job_name              = default'
        write(funit,'(a)') 'job_ntasks_per_socket = 1'
        call fclose(funit,errmsg="autogen_env_file qsys_funs ")
    end subroutine autogen_env_file

    !>  Writes the JOB_FINISHED_* file to mark end of computing unit job
    subroutine qsys_job_finished( p, source )
        ! generation of this file marks completion of the partition
        ! this file is empty 4 now but may contain run stats etc.
        use simple_params, only: params
        class(params),    intent(in) :: p
        character(len=*), intent(in) :: source
        if( p%l_chunk_distr )then
            call simple_touch('JOB_FINISHED_'//int2str_pad(p%chunk,p%numlen), errmsg="qsys_job_finished")
        endif
        if( p%l_distr_exec )then
            call simple_touch('JOB_FINISHED_'//int2str_pad(p%part,p%numlen), errmsg="qsys_job_finished")
        endif
        if( p%stream.eq.'yes' )then
            call simple_touch('JOB_FINISHED_'//int2str_pad(p%part,p%numlen), errmsg="qsys_job_finished ")
        endif
    end subroutine qsys_job_finished

    !>  returns when the inputted file exists in cwd
    subroutine qsys_watcher_1( fname, wtime )
        character(len=*),  intent(in) :: fname
        integer, optional, intent(in) :: wtime
        integer :: wwtime
        logical :: there
        wwtime = SHORTTIME
        if( present(wtime) ) wwtime = wtime
        do
            there = file_exists(trim(fname))
            if( there ) exit
            call simple_sleep(wwtime)
        end do
    end subroutine qsys_watcher_1

    !>  returns when the inputted files exist in cwd
    subroutine qsys_watcher_2( fnames, wtime )
        character(len=STDLEN), intent(in)    :: fnames(:)
        integer, optional,     intent(in)    :: wtime
        integer, parameter   :: MAXITS=20000
        integer              :: wwtime, nfiles, ifile, nlogic, i
        logical              :: doreturn, fexists
        wwtime = SHORTTIME
        if( present(wtime) ) wwtime = wtime
        nfiles = size(fnames)
        do i=1,MAXITS
            doreturn = .true.
            do ifile=1,nfiles
                fexists = file_exists(trim(fnames(ifile)))
                if( .not. fexists ) doreturn = .false.
            end do
            if( doreturn )then
                return
            else
                call simple_sleep(wwtime)
            endif
        end do
    end subroutine qsys_watcher_2

    subroutine exec_simple_prg( exec_bin, cline )
        use simple_cmdline, only: cmdline
        use simple_chash,   only: chash
        character(len=*),  intent(in) :: exec_bin
        class(cmdline),    intent(in) :: cline
        type(chash) :: job_descr
        character(len=1024) :: exec_str
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        exec_str = trim(exec_bin)//' '//job_descr%chash2str()
        write(*,'(a)') '>>> EXECUTING COMMAND:'
        write(*,'(a)') trim(exec_str)
        call exec_cmdline(exec_str)
    end subroutine exec_simple_prg

end module simple_qsys_funs
