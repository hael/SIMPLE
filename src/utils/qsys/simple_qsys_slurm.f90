! batch-processing manager - SLURM
module simple_qsys_slurm
include 'simple_lib.f08'
use simple_qsys_base, only: qsys_base
implicit none

public :: qsys_slurm
private

integer, parameter    :: MAXENVITEMS=100
character(len=STDLEN) :: stderrout

type, extends(qsys_base) :: qsys_slurm
    private
    type(chash) :: env !< defines the SLURM environment
  contains
    procedure :: new               => new_slurm_env
    procedure :: submit_cmd        => get_slurm_submit_cmd
    procedure :: write_instr       => write_slurm_header
    procedure :: write_array_instr => write_slurm_array_header
    procedure :: kill              => kill_slurm_env
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
        call self%env%push('job_memory_per_task',   '#SBATCH --mem')
        call self%env%push('job_time',              '#SBATCH --time')
        ! standard error & output folder
        stderrout = PATH_HERE//trim(STDERROUT_DIR)
        call simple_mkdir(trim(stderrout),errmsg="qsys_slurm::new slurm env")
    end subroutine new_slurm_env

    !> \brief  is a getter
    function get_slurm_submit_cmd( self ) result( cmd )
        class(qsys_slurm), intent(in) :: self
        character(len=:), allocatable :: cmd
        cmd = self%env%get('qsys_submit_cmd')
    end function get_slurm_submit_cmd

    !> \brief  writes the header instructions
    subroutine write_slurm_header( self, q_descr, fhandle )
        class(qsys_slurm), intent(in) :: self
        class(chash),      intent(in) :: q_descr
        integer, optional, intent(in) :: fhandle
        character(len=:), allocatable :: key, sbatch_cmd, sbatch_val
        integer :: i, which
        logical :: write2file
        write2file = .false.
        if( present(fhandle) ) write2file = .true.
        do i=1,q_descr%size_of()
            key   = q_descr%get_key(i)
            which = self%env%lookup(key)
            if( which > 0 )then
                sbatch_cmd = self%env%get(which)
                sbatch_val = q_descr%get(i)
                if( write2file )then
                    write(fhandle,'(a)') sbatch_cmd//'='//sbatch_val
                else
                    write(logfhandle,'(a)') sbatch_cmd//'='//sbatch_val
                endif
                deallocate(sbatch_cmd,sbatch_val)
            endif
            deallocate(key)
        end do
        ! write default instructions
        if( write2file )then
            write(fhandle,'(a)') '#SBATCH --output='//trim(stderrout)//'outfile.%j'
            write(fhandle,'(a)') '#SBATCH --error='//trim(stderrout)//'errfile.%j'
            write(fhandle,'(a)') '#SBATCH --mail-type=FAIL'
        else
            write(logfhandle,'(a)') '#SBATCH --output='//trim(stderrout)//'outfile.%j'
            write(logfhandle,'(a)') '#SBATCH --error='//trim(stderrout)//'errfile.%j'
            write(logfhandle,'(a)') '#SBATCH --mail-type=FAIL'
        endif
    end subroutine write_slurm_header

    !> \brief  writes the array header instructions
    subroutine write_slurm_array_header( self, q_descr, parts_fromto, fhandle, nactive )
        class(qsys_slurm), intent(in) :: self
        class(chash),      intent(in) :: q_descr
        integer,           intent(in) :: parts_fromto(2)
        integer, optional, intent(in) :: fhandle, nactive
        character(len=:), allocatable :: key, sbatch_cmd, sbatch_val
        integer :: i, which
        logical :: write2file
        write2file = .false.
        if( present(fhandle) ) write2file = .true.
        do i=1,q_descr%size_of()
            key   = q_descr%get_key(i)
            which = self%env%lookup(key)
            if( which > 0 )then
                sbatch_cmd = self%env%get(which)
                sbatch_val = q_descr%get(i)
                if( write2file )then
                    write(fhandle,'(a)') sbatch_cmd//'='//sbatch_val
                else
                    write(logfhandle,'(a)') sbatch_cmd//'='//sbatch_val
                endif
                deallocate(sbatch_cmd,sbatch_val)
            endif
            deallocate(key)
        end do
        ! write default instructions
        if( write2file )then
            write(fhandle,'(a)') '#SBATCH --output='//trim(stderrout)//'outfile.%A_%a'
            write(fhandle,'(a)') '#SBATCH --error='//trim(stderrout)//'errfile.%A_%a'
            write(fhandle,'(a)') '#SBATCH --mail-type=FAIL'
            if( present(nactive) )then
                ! subtract 1 from array indexes as bash arrays are zero indexed
                write(fhandle,'(a)') '#SBATCH --array='//int2str(parts_fromto(1) - 1)//'-'//int2str(parts_fromto(2) - 1)//'%'//int2str(nactive)
            else
                ! subtract 1 from array indexes as bash arrays are zero indexed
                write(fhandle,'(a)') '#SBATCH --array='//int2str(parts_fromto(1) - 1)//'-'//int2str(parts_fromto(2) - 1)
            endif
            write(fhandle,'(a)') 'echo $SLURM_ARRAY_JOB_ID > SLURM_ARRAY_JOB_ID'
        else
            write(logfhandle,'(a)') '#SBATCH --output='//trim(stderrout)//'outfile.%A_%a'
            write(logfhandle,'(a)') '#SBATCH --error='//trim(stderrout)//'errfile.%A_%a'
            write(logfhandle,'(a)') '#SBATCH --mail-type=FAIL'
            if( present(nactive) )then
                ! subtract 1 from array indexes as bash arrays are zero indexed
                write(logfhandle,'(a)') '#SBATCH --array='//int2str(parts_fromto(1) - 1)//'-'//int2str(parts_fromto(2) - 1)//'%'//int2str(nactive)
            else
                ! subtract 1 from array indexes as bash arrays are zero indexed
                write(logfhandle,'(a)') '#SBATCH --array='//int2str(parts_fromto(1) - 1)//'-'//int2str(parts_fromto(2) - 1)
            endif
            write(logfhandle,'(a)') 'echo $SLURM_ARRAY_JOB_ID > SLURM_ARRAY_JOB_ID'
        endif
    end subroutine write_slurm_array_header

    !> \brief  is a destructor
    subroutine kill_slurm_env( self )
        class(qsys_slurm), intent(inout) :: self
        call self%env%kill
    end subroutine kill_slurm_env

end module simple_qsys_slurm
