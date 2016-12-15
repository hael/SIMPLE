module simple_qsys_base
implicit none

public :: qsys_base
private

type, abstract :: qsys_base
  contains
    procedure(generic_new),         deferred :: new
    procedure(generic_submit_cmd),  deferred :: submit_cmd
    procedure(generic_write_instr), deferred :: write_instr
    procedure(generic_kill),        deferred :: kill
end type qsys_base

abstract interface

    !>  \brief  constructor
    subroutine generic_new( self )
        import :: qsys_base
        class(qsys_base), intent(inout) :: self
    end subroutine generic_new
    
    !>  \brief  getter that returns the submit command of the qsys
    function generic_submit_cmd( self ) result ( submit_cmd )
        import :: qsys_base
        class(qsys_base), intent(in)  :: self
        character(len=:), allocatable :: submit_cmd
    end function generic_submit_cmd
    
    !>  \brief  writes a header instruction for the submit script
    subroutine generic_write_instr( self, job_descr, fhandle )
        use simple_chash, only: chash
        import :: qsys_base
        class(qsys_base),  intent(in) :: self
        class(chash),      intent(in) :: job_descr
        integer, optional, intent(in) :: fhandle   
    end subroutine generic_write_instr

    !>  \brief  destructor
    subroutine generic_kill( self )
        import :: qsys_base
        class(qsys_base), intent(inout) :: self
    end subroutine generic_kill
    
end interface

end module simple_qsys_base
