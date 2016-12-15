module simple_qsys_local
use simple_qsys_base, only: qsys_base
use simple_chash,     only: chash
implicit none

public :: qsys_local
private

integer, parameter :: MAXENVITEMS=100

type, extends(qsys_base) :: qsys_local
    private
    type(chash) :: env !< defines the local environment
  contains
    procedure :: new         => new_local_env
    procedure :: submit_cmd  => get_local_submit_cmd
    procedure :: write_instr => write_local_header
    procedure :: kill        => kill_local_env
end type qsys_local

contains
    
    !> \brief  is a constructor
    subroutine new_local_env( self )
        class(qsys_local), intent(inout) :: self
        ! make the container
        call self%env%new(MAXENVITEMS)
        ! define the environment:
        call self%env%push('qsys_submit_cmd', 'nohup')
    end subroutine new_local_env
    
    !> \brief  is a getter
    function get_local_submit_cmd( self ) result( cmd )
        class(qsys_local), intent(in) :: self
        character(len=:), allocatable :: cmd
        cmd = self%env%get('qsys_submit_cmd')
    end function get_local_submit_cmd

    !> \brief  writes the header instructions
    subroutine write_local_header( self, job_descr, fhandle )
        class(qsys_local), intent(in) :: self
        class(chash),      intent(in) :: job_descr
        integer, optional, intent(in) :: fhandle
    end subroutine write_local_header
    
    !> \brief  is a destructor
    subroutine kill_local_env( self )
        class(qsys_local), intent(inout) :: self
        call self%env%kill
    end subroutine kill_local_env

end module simple_qsys_local
