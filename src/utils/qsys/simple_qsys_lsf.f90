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
        ! ### QSYS PARAMETERS
        call self%env%push('qsys_partition',        '#BSUB -q')   ! queue name
        call self%env%push('qsys_submit_cmd',       'bsub <')
        ! ### JOB PARAMETERS
        call self%env%push('job_name',              '#BSUB -J')
        call self%env%push('job_cpus_per_task',     '#BSUB -n')
        call self%env%push('job_ntasks_per_socket', '#BSUB -R "span[hosts=') ! completed in write_* (see below)
        call self%env%push('job_ntasks',            '#BSUB -R "span[hosts=') ! completed in write_* (see below)
        call self%env%push('job_memory_per_task',   '#BSUB -R "rusage[mem=') ! completed in write_* (see below)
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

     !> \brief  writes the header instructions
    subroutine write_lsf_header( self, q_descr, fhandle )
        class(qsys_lsf),   intent(in) :: self
        class(chash),      intent(in) :: q_descr
        integer, optional, intent(in) :: fhandle
        character(len=:), allocatable :: key, bsub_cmd, bsub_val, tmpstr
        integer :: i, which
        logical :: write2file
        write2file = .false.
        if( present(fhandle) ) write2file = .true.
        do i=1,q_descr%size_of()
            key   = q_descr%get_key(i)
            which = self%env%lookup(key)
            if( which > 0 )then
                bsub_cmd = self%env%get(which)
                bsub_val = q_descr%get(i)
                if( trim(key) == 'job_name' ) then
                    tmpstr   = '"'//trim(bsub_val)//'"'
                    bsub_val = tmpstr
                end if
                if( trim(key) == 'job_ntasks_per_socket' .or. trim(key) == 'job_memory_per_task' ) then
                    tmpstr   = trim(bsub_val)//']"'
                    bsub_val = tmpstr
                endif
                if( write2file )then
                    write(fhandle,'(a)') bsub_cmd//' '//bsub_val
                else
                    write(logfhandle,'(a)') bsub_cmd//' '//bsub_val
                endif
                deallocate(bsub_cmd,bsub_val)
            endif
            deallocate(key)
        end do
        ! write default instructions
        if( write2file )then
            write(fhandle,'(a)') '#BSUB -oo ' //trim(stderrout)//'outfile.%J'
            write(fhandle,'(a)') '#BSUB -eo ' //trim(stderrout)//'errfile.%J'
        else
            write(logfhandle,'(a)') '#BSUB -oo ' //trim(stderrout)//'outfile.%J'
            write(logfhandle,'(a)') '#BSUB -eo ' //trim(stderrout)//'errfile.%J'
        endif
    end subroutine write_lsf_header

    subroutine write_lsf_array_header( self, q_descr, parts_fromto, fhandle, nactive )
        class(qsys_lsf),   intent(in) :: self
        class(chash),      intent(in) :: q_descr
        integer,           intent(in) :: parts_fromto(2)
        integer, optional, intent(in) :: fhandle, nactive
    end subroutine write_lsf_array_header

    !> destructor
    subroutine kill_lsf_env( self )
        class(qsys_lsf), intent(inout) :: self
        call self%env%kill
    end subroutine kill_lsf_env

end module simple_qsys_lsf
