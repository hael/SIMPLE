! batch-processing manager - PBS
module simple_qsys_pbs
include 'simple_lib.f08'
use simple_qsys_base, only: qsys_base
implicit none

public :: qsys_pbs
private

integer, parameter :: MAXENVITEMS=100

type, extends(qsys_base) :: qsys_pbs
    private
    type(chash) :: env !< defines the PBS environment
  contains
    procedure :: new               => new_pbs_env
    procedure :: submit_cmd        => get_pbs_submit_cmd
    procedure :: write_instr       => write_pbs_header
    procedure :: write_array_instr => write_pbs_array_header
    procedure :: kill              => kill_pbs_env
end type qsys_pbs

contains

    !> \brief  is a constructor
    subroutine new_pbs_env( self )
        class(qsys_pbs), intent(inout) :: self
        ! make the container
        call self%env%new(MAXENVITEMS)
        ! define the environment:
        ! ### USER PARAMETERS
        call self%env%push('user_account',          '#PBS -A')
        call self%env%push('user_email',            '#PBS -M')
        ! ### QSYS PARAMETERS
        call self%env%push('qsys_partition',        '#PBS -q')
        call self%env%push('qsys_qos',              '#PBS -l qos')
        !call self%env%push('qsys_reservation',      '#PBS --reservation') ??
        call self%env%push('qsys_submit_cmd',       'qsub')
        ! ### JOB PARAMETERS
        call self%env%push('job_name',              '#PBS -N')
        !call self%env%push('job_ntasks',            '')                   ! unnecessary
        !call self%env%push('job_ntasks_per_socket', '')                   ! unnecessary
        call self%env%push('job_cpus_per_task',     '#PBS -l nodes=1:ppn') ! overridden by p%nthr
        call self%env%push('job_memory_per_task',   '#PBS -l mem')
        call self%env%push('job_time',              '#PBS -l walltime')
    end subroutine new_pbs_env

    !> \brief  is a getter
    function get_pbs_submit_cmd( self ) result( cmd )
        class(qsys_pbs),   intent(in) :: self
        character(len=:), allocatable :: cmd
        cmd = self%env%get('qsys_submit_cmd')
    end function get_pbs_submit_cmd

    !> \brief  writes the header instructions
    subroutine write_pbs_header( self, q_descr, fhandle )
        class(qsys_pbs),   intent(in) :: self
        class(chash),      intent(in) :: q_descr
        integer, optional, intent(in) :: fhandle
        character(len=:), allocatable :: key, pbs_cmd, pbs_val
        real    :: rval
        integer :: i, which
        logical :: write2file
        write2file = .false.
        if( present(fhandle) ) write2file = .true.
        do i=1,q_descr%size_of()
            if(allocated(key))deallocate(key)
            key     = q_descr%get_key(i)
            which   = self%env%lookup(key)
            if(which == 0)cycle
            pbs_cmd = self%env%get(which)
            pbs_val = q_descr%get(i)
            select case(key)
                case('job_time', 'job_cpus_per_task', 'qsys_qos')
                    call write_formatted(pbs_cmd, pbs_val, '=')
                case('job_memory_per_task')
                    ! memory in kilobytes
                    rval    = str2real(pbs_val) / 1024.
                    pbs_val = trim(int2str(ceiling(rval)))//'kb'
                    call write_formatted(pbs_cmd, pbs_val, '=')
                case DEFAULT
                    call write_formatted(pbs_cmd, pbs_val)
            end select
            if(allocated(pbs_cmd))deallocate(pbs_cmd)
            if(allocated(pbs_val))deallocate(pbs_val)
        end do
        if(allocated(key))deallocate(key)
        ! write default instructions
        if( write2file )then
            write(fhandle,'(a)') '#PBS -V'  ! environment copy
            write(fhandle,'(a)') '#PBS -o outfile.%j'
            write(fhandle,'(a)') '#PBS -e errfile.%j'
            write(fhandle,'(a)') '#PBS -m a'
        else
            write(logfhandle,'(a)') '#PBS -V'  ! environment copy
            write(logfhandle,'(a)') '#PBS -o outfile.%j'
            write(logfhandle,'(a)') '#PBS -e errfile.%j'
            write(logfhandle,'(a)') '#PBS -m a'
        endif

        contains

            subroutine write_formatted(cmd, val, delimiter)
                character(len=*),           intent(in) :: cmd, val
                character(len=*), optional, intent(in) :: delimiter
                if(present(delimiter))then
                    if( write2file )then
                        write(fhandle,'(a)') trim(cmd)//trim(delimiter)//trim(val)
                    else
                        write(logfhandle,'(a)') trim(cmd)//trim(delimiter)//trim(val)
                    endif
                else
                    if( write2file )then
                        write(fhandle,'(a)') trim(cmd)//' '//trim(val)
                    else
                        write(logfhandle,'(a)') trim(cmd)//' '//trim(val)
                    endif
                endif
            end subroutine write_formatted

    end subroutine write_pbs_header

    !> \brief  writes the array header instructions
    subroutine write_pbs_array_header( self, q_descr, nparts, fhandle, nactive )
        class(qsys_pbs),   intent(in) :: self
        class(chash),      intent(in) :: q_descr
        integer,           intent(in) :: nparts
        integer, optional, intent(in) :: fhandle, nactive
    end subroutine write_pbs_array_header

    !> \brief  is a destructor
    subroutine kill_pbs_env( self )
        class(qsys_pbs), intent(inout) :: self
        call self%env%kill
    end subroutine kill_pbs_env

end module simple_qsys_pbs
