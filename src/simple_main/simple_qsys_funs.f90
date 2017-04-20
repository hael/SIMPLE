module simple_qsys_funs
use simple_defs     
use simple_syscalls ! use all in there
use simple_strings  ! use all in there
implicit none

interface qsys_watcher
    module procedure qsys_watcher_1
    module procedure qsys_watcher_2
end interface

integer, parameter :: SHORTTIME = 3

contains

    subroutine qsys_cleanup( p )
        use simple_params,       only: params
        use simple_filehandling, only: del_file, del_files
        class(params),     intent(in) :: p
        character(len=:), allocatable :: rec_base_str, rho_base_str
        integer, parameter :: NUMLEN_STATE = 2, NUMLEN_ITER = 3
        integer :: istate        
        ! single files
        call del_file('FOO')
        call del_file('fort.0')
        call del_file('qsys_submit_jobs')
        call del_file('simple_script_single')
        call del_file('ftab_from_sys_find_last_fname.txt')
        call del_file('VOLASSEMBLE')
        call del_file('CAVGASSEMBLE')
        call del_file('SYMSRCH')
        ! part numbered files
        call del_files('OUT',                    p%nparts)
        call del_files('algndoc_',               p%nparts, ext='.txt')
        call del_files('unidoc_',                p%nparts, ext='.txt')
        call del_files('cavgs_part',             p%nparts, ext=p%ext )
        call del_files('JOB_FINISHED_',          p%nparts)
        call del_files('ctfsqsums_part',         p%nparts, ext=p%ext )
        call del_files('distr_simple_script_',   p%nparts)
        call del_files('ctffind_output_part',    p%nparts, ext='.txt')
        call del_files('ctffind_ctrl_file_part', p%nparts, ext='.txt')
        ! state and part numbered files
        do istate=1,MAXS
            allocate(rec_base_str, source='recvol_state'//int2str_pad(istate,NUMLEN_STATE)//'_part')
            allocate(rho_base_str, source='rho_'//rec_base_str)
            call del_files(rec_base_str, p%nparts, ext=p%ext)
            call del_files(rho_base_str, p%nparts, ext=p%ext)
            call del_files(rec_base_str, p%nparts, ext=p%ext, suffix='_even')
            call del_files(rec_base_str, p%nparts, ext=p%ext, suffix='_odd')
            call del_files(rho_base_str, p%nparts, ext=p%ext, suffix='_even')
            call del_files(rho_base_str, p%nparts, ext=p%ext, suffix='_odd')
            deallocate(rec_base_str,rho_base_str)
        end do
    end subroutine qsys_cleanup

    function stack_is_split( stkext, npart ) result( is_split )
        character(len=4),  intent(in) :: stkext
        integer,           intent(in) :: npart  
        character(len=:), allocatable :: stack_part_fname
        logical, allocatable          :: stack_parts_exist(:) 
        integer :: ipart, numlen
        logical :: is_split
        allocate( stack_parts_exist(npart) )
        numlen = len(int2str(npart))
        do ipart=1,npart
            allocate(stack_part_fname, source='stack_part'//int2str_pad(ipart,numlen)//stkext)
            stack_parts_exist(ipart) = file_exists(stack_part_fname)
            deallocate(stack_part_fname)
        end do
        is_split = all(stack_parts_exist)
    end function stack_is_split

    subroutine stack_parts_of_correct_sizes( stkext, parts, box )
        character(len=4),  intent(in) :: stkext
        integer,           intent(in) :: parts(:,:)
        integer,           intent(in) :: box
        character(len=:), allocatable :: stack_part_fname
        integer :: ipart, numlen, sz_correct, sz, ldim(3), npart
        npart  = size(parts,1)
        numlen = len(int2str(npart))
        do ipart=1,npart
            sz_correct = parts(ipart,2)-parts(ipart,1)+1
            allocate(stack_part_fname, source='stack_part'//int2str_pad(ipart,numlen)//stkext)
            call find_ldim_nptcls(stack_part_fname, ldim, sz)
            if( sz /= sz_correct )then
                write(*,*) 'size of ', stack_part_fname, ' is ', sz, ' not ', sz_correct, 'as expected'
                stop 'simple_qsys_funs :: stack_parts_of_correct_sizes'
            endif
            if( ldim(1) == box .and. ldim(2) == box )then
                ! dimension ok
            else
                write(*,*) 'ldim of ', stack_part_fname, ' is ', [ldim(1),ldim(2),1], ' not ', [box,box,1], 'as expected'
                stop 'simple_qsys_funs :: stack_parts_of_correct_sizes'
            endif
            deallocate(stack_part_fname)
        end do
    end subroutine stack_parts_of_correct_sizes

    subroutine terminate_if_prg_in_cwd( prg )
        character(len=*), intent(in) :: prg
        integer                      :: file_stat, funit, ios, file_rec, i, pos
        character(len=STDLEN)        :: pid, pid_dir, pwd, cmd
        cmd = 'pgrep '//trim(prg)//' > pid_from_fortran.txt'
        call exec_cmdline(cmd)
        ! read in the pid
        funit = get_fileunit()
        open(unit=funit, status='old', action='read', file='pid_from_fortran.txt', iostat=file_stat)
        read(funit,'(a)') pid
        close(funit)
        ! get a string that contains the directory it is running in with lsof -p
        cmd = 'lsof -p '//trim(pid)//' | grep cwd > pid_dir_from_fortran.txt'
        call exec_cmdline(cmd)
        open(unit=funit, status='old', action='read', file='pid_dir_from_fortran.txt', iostat=file_stat)
        read(funit,'(a)') pid_dir
        close(funit)
        ! get pwd
        call get_environment_variable('PWD', pwd)
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
        simple_path = sys_get_env_var('SIMPLE_PATH')
        if( allocated(simple_path) )then
            simple_path_env = env%get('simple_path')
            if( str_has_substr(simple_path, simple_path_env) .or. str_has_substr(simple_path_env, simple_path) )then
                ! all ok
            else
                write(*,*) 'WARNING! simple absolute paths in shell and simple_distr_config.env do not agree'
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
        use simple_syscalls, only: sys_get_env_var
        character(len=*), intent(in)  :: env_file
        integer                       :: funit, file_stat
        character(len=:), allocatable :: simple_path, simple_email, simple_qsys
        simple_path = sys_get_env_var('SIMPLE_PATH')
        if( .not. allocated(simple_path) )&
        stop 'need SIMPLE_PATH env var to auto-generate simple_distr_config.env; simple_qsys_funs :: autogen_env_file'
        simple_qsys = sys_get_env_var('SIMPLE_QSYS')
        if( .not. allocated(simple_qsys) )&
        stop 'need SIMPLE_QSYS env var to auto-generate simple_distr_config.env; simple_qsys_funs :: autogen_env_file'
        simple_email = sys_get_env_var('SIMPLE_EMAIL')
        if( .not. allocated(simple_email) ) allocate(simple_email, source='me.myself@uni.edu')
        funit = get_fileunit()
        open(unit=funit, status='replace', action='write', file=trim(env_file), iostat=file_stat)
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
        close(funit)
    end subroutine autogen_env_file

    !>  Writes the JOB_FINISHED_* file to mark end of computing unit job 
    subroutine qsys_job_finished( p, source )
        ! generation of this file marks completion of the partition
        ! this file is empty 4 now but may contain run stats etc.
        use simple_params, only: params
        class(params),    intent(in) :: p
        character(len=*), intent(in) :: source 
        integer :: fnr, file_stat
        if( p%l_chunk_distr )then
            fnr = get_fileunit()
            open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%chunk,p%numlen),&
            &STATUS='REPLACE', action='WRITE', iostat=file_stat)
            call fopen_err( source, file_stat )
            close( unit=fnr )
        endif
        if( p%l_distr_exec )then
            fnr = get_fileunit()
            open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen),&
            &STATUS='REPLACE', action='WRITE', iostat=file_stat)
            call fopen_err( source, file_stat )
            close( unit=fnr )
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
    subroutine qsys_watcher_2( fnames, files_exist, wtime )
        character(len=STDLEN), intent(in)    :: fnames(:)
        logical,               intent(inout) :: files_exist(:)
        integer, optional,     intent(in)    :: wtime
        integer, parameter :: MAXITS=1000
        integer            :: wwtime, nfiles, ifile, nlogic, i
        logical            :: fexists, doreturn
        wwtime = SHORTTIME
        if( present(wtime) ) wwtime = wtime
        nfiles = size(fnames)
        nlogic = size(files_exist)
        if( nlogic .ne. nfiles ) stop 'nonconforming arrays in simple_qsys_funs :: qsys_watcher_2' 
        do i=1,MAXITS
            doreturn = .false.
            do ifile=1,nfiles
                fexists = file_exists(trim(fnames(ifile)))
                if( fexists .neqv. files_exist(ifile) )then
                    files_exist(ifile) = fexists
                    doreturn = .true.
                endif
            end do
            if( all(files_exist) .or. doreturn )then
                return
            else
                call sleep(wtime)
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
