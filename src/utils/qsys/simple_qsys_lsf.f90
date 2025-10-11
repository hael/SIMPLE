! batch-processing manager - LSF
module simple_qsys_lsf
include 'simple_lib.f08'
use simple_qsys_base, only: qsys_base
implicit none

public :: qsys_lsf
private

integer, parameter    :: MAXENVITEMS=100
character(len=STDLEN) :: stderrout

type, extends(qsys_base) :: qsys_lsf
    private
    type(chash) :: env   !< defines the LSF environment
contains
    procedure :: new               => new_lsf_env
    procedure :: submit_cmd        => get_lsf_submit_cmd
    procedure :: write_instr       => write_lsf_header
    procedure :: write_array_instr => write_lsf_array_header
    procedure :: kill              => kill_lsf_env
end type qsys_lsf

contains

    !> constructor
    subroutine new_lsf_env( self )
        class(qsys_lsf), intent(inout) :: self
        call self%env%new(MAXENVITEMS)
        ! ### USER PARAMETERS
        call self%env%push('user_account',          '#BSUB -P')   ! project/account
        call self%env%push('user_email',            '#BSUB -u')   ! email for notifications
        ! ### QSYS PARAMETERS
        call self%env%push('qsys_partition',        '#BSUB -q')   ! queue
        call self%env%push('qsys_qos',              '#BSUB -sla') ! site-specific; SLA/QoS if enabled
        call self%env%push('qsys_reservation',      '#BSUB -U')   ! (best-effort mapping; site-specific)
        call self%env%push('qsys_submit_cmd',       'bsub')
        ! ### JOB PARAMETERS
        call self%env%push('job_name',              '#BSUB -J')
        call self%env%push('job_ntasks',            '#BSUB -n')             ! total slots
        call self%env%push('job_ntasks_per_socket', '#BSUB -R "span[ptile') ! completed in write_* (see below)
        call self%env%push('job_cpus_per_task',     '#BSUB -R "span[ptile') ! completed in write_* (see below)
        call self%env%push('job_memory_per_task',   '#BSUB -R "rusage[mem') ! completed in write_* (see below)
        call self%env%push('job_time',              '#BSUB -W')
        ! standard error & output folder
        stderrout = PATH_HERE//trim(STDERROUT_DIR)
        call simple_mkdir(trim(stderrout),errmsg="qsys_lsf::new lsf env")
    end subroutine new_lsf_env

    !> getter
    function get_lsf_submit_cmd( self ) result( cmd )
        class(qsys_lsf), intent(in) :: self
        character(len=:), allocatable :: cmd
        cmd = self%env%get('qsys_submit_cmd')
    end function get_lsf_submit_cmd

    !> write header instructions (single job)
    subroutine write_lsf_header( self, q_descr, fhandle )
        class(qsys_lsf),   intent(in) :: self
        class(chash),      intent(in) :: q_descr
        integer, optional, intent(in) :: fhandle
        character(len=:), allocatable :: key, bsub_cmd, bsub_val
        integer  :: i, which
        logical  :: write2file
        logical  :: have_ptile, have_mem
        character(len=:), allocatable :: ptile_str, mem_str
        write2file = .false.; if( present(fhandle) ) write2file = .true.
        have_ptile = .false.; have_mem = .false.
        do i=1,q_descr%size_of()
            key   = q_descr%get_key(i)
            which = self%env%lookup(key)
            if( which > 0 )then
                bsub_cmd = self%env%get(which)
                bsub_val = q_descr%get(i)
                select case (trim(key))
                    case ('job_cpus_per_task')
                        ! map to span[ptile=<val>] (approximation of CPU-per-task semantics)
                        ptile_str  = trim(bsub_cmd)//'='//trim(bsub_val)//']"'
                        have_ptile = .true.
                        if( write2file )then
                            write(fhandle,'(a)') '#BSUB -R "span[ptile='//trim(bsub_val)//']"'
                        else
                            write(logfhandle,'(a)') '#BSUB -R "span[ptile='//trim(bsub_val)//']"'
                        endif
                    case ('job_ntasks_per_socket')
                        ! no clean LSF equivalent; best-effort: constrain ptile as a hint
                        if( write2file )then
                            write(fhandle,'(a)') '#BSUB -R "span[ptile='//trim(bsub_val)//']"'
                        else
                            write(logfhandle,'(a)') '#BSUB -R "span[ptile='//trim(bsub_val)//']"'
                        endif
                    case ('job_memory_per_task')
                        ! mem in MB → rusage[mem=...]
                        mem_str  = '#BSUB -R "rusage[mem='//trim(bsub_val)//']"'
                        have_mem = .true.
                        if( write2file )then
                            write(fhandle,'(a)') trim(mem_str)
                        else
                            write(logfhandle,'(a)') trim(mem_str)
                        endif
                    case DEFAULT
                        ! standard space-separated BSUB syntax
                        if( write2file )then
                            write(fhandle,'(a)') trim(bsub_cmd)//' '//trim(bsub_val)
                        else
                            write(logfhandle,'(a)') trim(bsub_cmd)//' '//trim(bsub_val)
                        endif
                end select
                deallocate(bsub_cmd,bsub_val)
            endif
            deallocate(key)
        end do
        ! defaults: stdout/err & email-at-end (no exact "FAIL only" equivalent)
        if( write2file )then
            write(fhandle,'(a)') '#BSUB -o ' // trim(stderrout)//'outfile.%J'
            write(fhandle,'(a)') '#BSUB -e ' // trim(stderrout)//'errfile.%J'
            write(fhandle,'(a)') '#BSUB -N'
        else
            write(logfhandle,'(a)') '#BSUB -o ' // trim(stderrout)//'outfile.%J'
            write(logfhandle,'(a)') '#BSUB -e ' // trim(stderrout)//'errfile.%J'
            write(logfhandle,'(a)') '#BSUB -N'
        endif
    end subroutine write_lsf_header

    !> write array header instructions
    subroutine write_lsf_array_header( self, q_descr, parts_fromto, fhandle, nactive )
        class(qsys_lsf),   intent(in) :: self
        class(chash),      intent(in) :: q_descr
        integer,           intent(in) :: parts_fromto(2)
        integer, optional, intent(in) :: fhandle, nactive
        character(len=:), allocatable :: key, bsub_cmd, bsub_val, jobname
        integer  :: i, which
        logical  :: write2file
        write2file = .false.; if( present(fhandle) ) write2file = .true.
        jobname = 'ARRAY'
        ! first pass: write normal options & capture job name
        do i=1,q_descr%size_of()
            key   = q_descr%get_key(i)
            which = self%env%lookup(key)
            if( which > 0 )then
                bsub_cmd = self%env%get(which)
                bsub_val = q_descr%get(i)
                if( trim(key) == 'job_name' ) then
                    if( allocated(jobname) ) deallocate(jobname)
                    jobname = trim(bsub_val)
                    ! do not emit a plain -J here; array -J will be written below
                else
                    select case (trim(key))
                        case ('job_cpus_per_task','job_ntasks_per_socket')
                            if( write2file )then
                                write(fhandle,'(a)') '#BSUB -R "span[ptile='//trim(bsub_val)//']"'
                            else
                                write(logfhandle,'(a)') '#BSUB -R "span[ptile='//trim(bsub_val)//']"'
                            endif
                        case ('job_memory_per_task')
                            if( write2file )then
                                write(fhandle,'(a)') '#BSUB -R "rusage[mem='//trim(bsub_val)//']"'
                            else
                                write(logfhandle,'(a)') '#BSUB -R "rusage[mem='//trim(bsub_val)//']"'
                            endif
                        case DEFAULT
                            if( write2file )then
                                write(fhandle,'(a)') trim(bsub_cmd)//' '//trim(bsub_val)
                            else
                                write(logfhandle,'(a)') trim(bsub_cmd)//' '//trim(bsub_val)
                            endif
                    end select
                endif
                deallocate(bsub_cmd,bsub_val)
            endif
            deallocate(key)
        end do
        ! defaults & array specifics (LSF arrays are 1-based; no index shift)
        if( write2file )then
            write(fhandle,'(a)') '#BSUB -o ' // trim(stderrout)//'outfile.%J_%I'
            write(fhandle,'(a)') '#BSUB -e ' // trim(stderrout)//'errfile.%J_%I'
            write(fhandle,'(a)') '#BSUB -N'
            if( present(nactive) )then
                write(fhandle,'(a)') '#BSUB -J "'//trim(jobname)//'['// &
                int2str(parts_fromto(1))//'-'//int2str(parts_fromto(2))//']%'//int2str(nactive)//'"'
            else
                write(fhandle,'(a)') '#BSUB -J "'//trim(jobname)//'['// &
                int2str(parts_fromto(1))//'-'//int2str(parts_fromto(2))//']"'
            endif
            write(fhandle,'(a)') 'echo $LSB_JOBID    > LSF_JOB_ID'
            write(fhandle,'(a)') 'echo $LSB_JOBINDEX > LSF_JOB_INDEX'
        else
            write(logfhandle,'(a)') '#BSUB -o ' // trim(stderrout)//'outfile.%J_%I'
            write(logfhandle,'(a)') '#BSUB -e ' // trim(stderrout)//'errfile.%J_%I'
            write(logfhandle,'(a)') '#BSUB -N'
            if( present(nactive) )then
                write(logfhandle,'(a)') '#BSUB -J "'//trim(jobname)//'['// &
                int2str(parts_fromto(1))//'-'//int2str(parts_fromto(2))//']%'//int2str(nactive)//'"'
            else
                write(logfhandle,'(a)') '#BSUB -J "'//trim(jobname)//'['// &
                int2str(parts_fromto(1))//'-'//int2str(parts_fromto(2))//']"'
            endif
            write(logfhandle,'(a)') 'echo $LSB_JOBID    > LSF_JOB_ID'
            write(logfhandle,'(a)') 'echo $LSB_JOBINDEX > LSF_JOB_INDEX'
        endif
    end subroutine write_lsf_array_header

    !> destructor
    subroutine kill_lsf_env( self )
        class(qsys_lsf), intent(inout) :: self
        call self%env%kill
    end subroutine kill_lsf_env

end module simple_qsys_lsf
