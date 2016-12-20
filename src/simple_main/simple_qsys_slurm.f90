module simple_qsys_slurm
use simple_qsys_base, only: qsys_base
use simple_chash,     only: chash
implicit none

public :: qsys_slurm
private

integer, parameter :: MAXENVITEMS=100

type, extends(qsys_base) :: qsys_slurm
    private
    type(chash) :: env !< defines the SLURM environment
  contains
    procedure :: new         => new_slurm_env
    procedure :: submit_cmd  => get_slurm_submit_cmd
    procedure :: write_instr => write_slurm_header
    procedure :: kill        => kill_slurm_env
end type qsys_slurm

contains
    
    !> \brief  is a constructor
    subroutine new_slurm_env( self )
        class(qsys_slurm), intent(inout) :: self
        ! make the container
        call self%env%new(MAXENVITEMS)
        ! define the environment:
        ! ### USER PARAMETERS
        call self%env%push('user_account',          '#SBATCH --account')
        call self%env%push('user_email',            '#SBATCH --mail-user')
        ! ### QSYS PARAMETERS
        call self%env%push('qsys_partition',        '#SBATCH --partition')
        call self%env%push('qsys_qos',              '#SBATCH --qos')
        call self%env%push('qsys_reservation',      '#SBATCH --reservation')
        call self%env%push('qsys_submit_cmd',       'sbatch')
        ! ### JOB PARAMETERS
        call self%env%push('job_name',              '#SBATCH --job-name')
        call self%env%push('job_ntasks',            '#SBATCH --ntasks')
        call self%env%push('job_ntasks_per_socket', '#SBATCH --ntasks-per-socket')
        call self%env%push('job_cpus_per_task',     '#SBATCH --cpus-per-task') ! overridden by p%nthr
        ! call self%env%push('job_gpus_per_task',     '#SBATCH --gpus-per-task') ! not a thing
        call self%env%push('job_memory_per_task',   '#SBATCH --mem')
        call self%env%push('job_time',              '#SBATCH --time')
    end subroutine new_slurm_env

    !> \brief  is a getter
    function get_slurm_submit_cmd( self ) result( cmd )
        class(qsys_slurm), intent(in) :: self
        character(len=:), allocatable :: cmd
        cmd = self%env%get('qsys_submit_cmd')
    end function get_slurm_submit_cmd

    !> \brief  writes the header instructions
    subroutine write_slurm_header( self, job_descr, fhandle )
        class(qsys_slurm), intent(in) :: self
        class(chash),      intent(in) :: job_descr
        integer, optional, intent(in) :: fhandle
        character(len=:), allocatable :: key, sbatch_cmd, sbatch_val
        integer :: i, which
        logical :: write2file
        write2file = .false.
        if( present(fhandle) ) write2file = .true.
        do i=1,job_descr%size_of_chash()
            key   = job_descr%get_key(i)
            which = self%env%lookup(key)
            if( which > 0 )then
                sbatch_cmd = self%env%get(which)
                sbatch_val = job_descr%get(i)
                if( write2file )then
                    write(fhandle,'(a)') sbatch_cmd//'='//sbatch_val
                else
                    write(*,'(a)') sbatch_cmd//'='//sbatch_val
                endif
                deallocate(sbatch_cmd,sbatch_val)
            endif
            deallocate(key)
        end do
        ! write default instructions
        if( write2file )then
            write(fhandle,'(a)') '#SBATCH --output=outfile.%j'
            write(fhandle,'(a)') '#SBATCH --error=errfile.%j'
            write(fhandle,'(a)') '#SBATCH --mail-type=FAIL'
        else
            write(*,'(a)') '#SBATCH --output=outfile.%j'
            write(*,'(a)') '#SBATCH --error=errfile.%j'
            write(*,'(a)') '#SBATCH --mail-type=FAIL'
        endif
    end subroutine write_slurm_header
    
    !> \brief  is a destructor
    subroutine kill_slurm_env( self )
        class(qsys_slurm), intent(inout) :: self
        call self%env%kill
    end subroutine kill_slurm_env

end module simple_qsys_slurm
