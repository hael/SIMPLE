module simple_qsys_funs
use simple_jiffys   ! singleton
use simple_defs     ! singleton
use simple_syscalls ! singleton
implicit none

contains

    subroutine qsys_cleanup_all
        character(len=STDLEN) :: cleanup_exec_cmd
        integer, parameter    :: SHORTTIME = 5
        call sleep(SHORTTIME)
        cleanup_exec_cmd = 'rm -f FOO'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f OUT*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f algndoc_*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f distr_script_*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f rho*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f fort.0'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f primedoc_*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f recvol*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f JOB_FINISHED_*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f errfile.*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f outfile.*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f ctfsqsums_part*'
        call exec_cmdline(cleanup_exec_cmd)
    end subroutine qsys_cleanup_all

    subroutine qsys_cleanup_iter()
        character(len=STDLEN) :: cleanup_exec_cmd
        integer, parameter    :: SHORTTIME = 5
        call sleep(SHORTTIME)
        cleanup_exec_cmd = 'rm -f FOO'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f OUT*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f algndoc_*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f distr_*script_*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f rho*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f fort.0'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f recvol_state*_part*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f JOB_FINISHED_*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f errfile.*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f outfile.*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f ctfsqsums_part*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f ctfsqsums_part*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f cavgs_part*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f *_recvol_state*_even*'
        call exec_cmdline(cleanup_exec_cmd)
        cleanup_exec_cmd = 'rm -f *_recvol_state*_odd*'
        call exec_cmdline(cleanup_exec_cmd)
    end subroutine qsys_cleanup_iter

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

    subroutine stack_parts_of_correct_sizes( stkext, parts )
        character(len=4),  intent(in) :: stkext
        integer,           intent(in) :: parts(:,:)
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
        call del_txtfile('pid_from_fortran.txt')
        call del_txtfile('pid_dir_from_fortran.txt')
    end subroutine terminate_if_prg_in_cwd

    subroutine parse_env_file( env )
        use simple_chash, only: chash
        type(chash), intent(inout)  :: env
        character(len=STDLEN) :: env_file, qsys_name
        integer               :: i,funit, nl, file_stat
        env_file = './simple_distr_config.env'
        if( .not. file_exists(trim(env_file)) ) call autogen_env_file(env_file)
        call env%read( trim(env_file) )
        ! User keys
        if( .not. env%isthere('simple_path') )stop 'Path to SIMPLE directory is required in simple_distr_config.env (simple_path)'        
        ! Queue keys
        if( .not. env%isthere('qsys_name') ) &
            & stop 'The type of the queuing system is required in simple_distr_config.env: qsys_name=<local|slurm|pbs>'
        qsys_name = env%get('qsys_name')
        if( qsys_name.ne.'local' .and. qsys_name.ne.'pbs' .and. qsys_name.ne.'slurm' ) &
            & stop 'Invalid qsys_name in simple_distr_config.env'
        if( qsys_name.ne.'local' )then
            if( .not. env%isthere('qsys_partition') )stop 'qsys_partition field is required in simple_distr_config.env'        
        endif
        ! Job keys
        if( qsys_name.ne.'local' )then
            if( .not. env%isthere('job_memory_per_task') )stop 'Job memory is required in simple_distr_config.env (job_memory)'
            if( .not. env%isthere('job_name') )             call env%push('job_name','simple_job')
            if( .not. env%isthere('job_time') )             call env%push('job_time', '0-23:59:0') ! TODO: use time_per_image
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

    subroutine setup_qsys_env( p, qsys_fac, myqsys, parts, qscripts, myq_descr )
        use simple_params,       only: params
        use simple_qsys_factory, only: qsys_factory
        use simple_qsys_base,    only: qsys_base
        use simple_qsys_ctrl,    only: qsys_ctrl
        use simple_chash,        only: chash
        use simple_map_reduce   ! singleton
        class(params),             intent(in)  :: p          !< parameters
        class(qsys_factory),       intent(out) :: qsys_fac   !< qsystem factory instance
        class(qsys_base), pointer, intent(out) :: myqsys     !< pointer to constructed object
        integer, allocatable,      intent(out) :: parts(:,:) !< indices of partitions
        class(qsys_ctrl),          intent(out) :: qscripts   !< qsys controller
        class(chash),              intent(out) :: myq_descr  !< user-provided queue and job specifics from simple_distr_config.env
        character(len=:), allocatable :: qsnam
        character(len=STDLEN)         :: simplepath_exec
        integer, parameter            :: MAXNKEYS=30
        if( .not. allocated(parts) )then
            ! generate partitions
            select case(p%split_mode)
                case('even')
                    parts = split_nobjs_even(p%nptcls, p%nparts)
                case('chunk')
                    parts = split_nobjs_in_chunks(p%nptcls, p%chunksz)
                case DEFAULT
                    write(*,*) 'split_mode: ', trim(p%split_mode)
                    stop 'Unsupported split_mode'
            end select
        endif
        ! PREPARE QUEUE DEPENDENT VARIABLES
        ! retrieve environment variables from file
        call myq_descr%new(MAXNKEYS)
        call parse_env_file(myq_descr) ! parse .env file
        qsnam = myq_descr%get('qsys_name')
        call qsys_fac%new(qsnam, myqsys)
        ! create the user specific qsys and qsys controller (script generator)
        simplepath_exec = trim(myq_descr%get('simple_path'))//'/bin/simple_exec'
        call qscripts%new( simplepath_exec, myqsys, parts, [1,p%nparts], p%ncunits )
        call myq_descr%set('job_cpus_per_task', int2str(p%nthr))   ! overrides env file
        call myq_descr%set('job_nparts',        int2str(p%nparts)) ! overrides env file
        deallocate(qsnam)
    end subroutine setup_qsys_env

    !>  Writes the JOB_FINISHED_* file to mark end of computing unit job 
    subroutine qsys_job_finished( p, source )
        ! generation of this file marks completion of the partition
        ! this file is empty 4 now but may contain run stats etc.
        use simple_params,       only: params
        class(params),         intent(in)    :: p
        character(len=STDLEN), intent(in) :: source 
        integer :: fnr, file_stat
        if( p%l_distr_exec )then
            fnr = get_fileunit()
            open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen),&
            &STATUS='REPLACE', action='WRITE', iostat=file_stat)
            call fopen_err( source, file_stat )
            close( unit=fnr )
        endif
    end subroutine qsys_job_finished

    subroutine copy_bin_files( source, destination, nstates )
        use simple_syscalls, only: sys_cp
        character(len=*), intent(in) :: source, destination
        integer, intent(in)          :: nstates
        integer                      :: s
        character(len=STDLEN) :: from, to, file, str_state
        do s=1,nstates
            str_state = 'state'//int2str_pad(s,2)
            from = trim(source)//trim(str_state)//'.bin'
            to   = trim(destination)//'/'
            call sys_cp( trim(from), trim(to) )
        enddo
    end subroutine copy_bin_files

end module simple_qsys_funs