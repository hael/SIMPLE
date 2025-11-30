! batch-processing manager - SunGrid Engine
module simple_qsys_sge
use simple_qsys_base, only: qsys_base
include 'simple_lib.f08'
implicit none

public :: qsys_sge
private

integer, parameter :: MAXENVITEMS=100

type, extends(qsys_base) :: qsys_sge
    private
    type(chash) :: env !< defines the SGE environment
  contains
    procedure :: new               => new_sge_env
    procedure :: submit_cmd        => get_sge_submit_cmd
    procedure :: write_instr       => write_sge_header
    procedure :: write_array_instr => write_sge_array_header
    procedure :: kill              => kill_sge_env
end type qsys_sge

contains

    !> \brief  is a constructor
    subroutine new_sge_env( self )
        class(qsys_sge), intent(inout) :: self
        ! make the container
        call self%env%new(MAXENVITEMS)
        ! define the environment:
        ! ### USER PARAMETERS
        call self%env%push('user_account',           '#$ -A')
        call self%env%push('user_email',             '#$ -M')
        call self%env%push('user_project',           '#$ -P')
        ! ### QSYS PARAMETERS
        call self%env%push('qsys_queue',             '#$ -q')
        call self%env%push('qsys_submit_cmd',         'qsub')
        ! ### JOB PARAMETERS
        call self%env%push('job_name',               '#$ -N')
        call self%env%push('job_addon_line',         '')
    end subroutine new_sge_env

    !> \brief  is a getter
    function get_sge_submit_cmd( self ) result( cmd )
        class(qsys_sge), intent(in) :: self
        type(string) :: cmd
        cmd = self%env%get('qsys_submit_cmd')
    end function get_sge_submit_cmd

    !> \brief  writes the header instructions
    subroutine write_sge_header( self, q_descr, fhandle )
        class(qsys_sge), intent(in)   :: self
        class(chash),      intent(in) :: q_descr
        integer, optional, intent(in) :: fhandle
        type(string) :: key, qsub_cmd, qsub_val, bind2socket, addon_line
        integer :: i, which
        logical :: write2file
        write2file = .false.
        if( present(fhandle) ) write2file = .true.
        do i=1,q_descr%size_of()
            key   = q_descr%get_key(i)
            which = self%env%lookup(key%to_char())
            if( which > 0 )then
                qsub_cmd = self%env%get(which)
                qsub_val = q_descr%get(i)
                if( key%to_char().eq. 'job_addon_line' )then
                    call addon_line%kill
                    addon_line = qsub_val
                    cycle
                else
                    if( write2file )then
                        write(fhandle,'(a)') qsub_cmd%to_char()//' '//qsub_val%to_char()
                    else
                        write(logfhandle,'(a)') qsub_cmd%to_char()//' '//qsub_val%to_char()
                    endif
                endif
                if(key%to_char() .eq. 'job_cpus_per_task')then
                    call bind2socket%kill
                    bind2socket = '#$ -binding linear:'//qsub_val%to_char()
                endif
                call qsub_cmd%kill
                call qsub_val%kill
            endif
            call key%kill
        end do
        ! write default instructions
        if( write2file )then
            write(fhandle,'(a)') '#$ -V'
            write(fhandle,'(a)') '#$ -S /bin/bash'
            write(fhandle,'(a)') '#$ -cwd'
            write(fhandle,'(a)') '#$ -o outfile -j y'
            if( bind2socket%is_allocated() )then
                 write(fhandle,'(a)') bind2socket%to_char()
            endif
            if( addon_line%is_allocated() )then
                 write(fhandle, '(a)') addon_line%to_char()
            endif
        else
            write(logfhandle,'(a)') '#$ -V'
            write(logfhandle,'(a)') '#$ -S /bin/bash'
            write(logfhandle,'(a)') '#$ -cwd'
            write(logfhandle,'(a)') '#$ -o outfile -j y'
            if( bind2socket%is_allocated() )then
                 write(fhandle,'(a)') bind2socket%to_char()
            endif
            if( addon_line%is_allocated() )then
                 write(fhandle, '(a)') addon_line%to_char()
            endif
        endif
    end subroutine write_sge_header

    !> \brief  writes the array header instructions
    subroutine write_sge_array_header( self, q_descr, parts_fromto, fhandle, nactive )
        class(qsys_sge),   intent(in) :: self
        class(chash),      intent(in) :: q_descr
        integer,           intent(in) :: parts_fromto(2)
        integer, optional, intent(in) :: fhandle, nactive
    end subroutine write_sge_array_header

    !> \brief  is a destructor
    subroutine kill_sge_env( self )
        class(qsys_sge), intent(inout) :: self
        call self%env%kill
    end subroutine kill_sge_env

end module simple_qsys_sge
