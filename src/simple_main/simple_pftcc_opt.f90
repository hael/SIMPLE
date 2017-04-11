module simple_pftcc_opt
implicit none

public :: pftcc_opt
private

type, abstract :: pftcc_opt
  contains
    procedure(generic_new),          deferred :: new
    procedure(generic_set_indices),  deferred :: set_indices
    procedure(generic_costfun),      deferred :: costfun
    procedure(generic_minimize),     deferred :: minimize
    procedure(generic_get_nevals),   deferred :: get_nevals
end type

abstract interface

    !>  \brief  is a constructor
    subroutine generic_new( self, pftcc, lims, shbarrier, nrestarts) 
        use simple_polarft_corrcalc, only: polarft_corrcalc
        import :: pftcc_opt
        class(pftcc_opt),                intent(inout) :: self
        class(polarft_corrcalc), target, intent(in)    :: pftcc
        real,                            intent(in)    :: lims(:,:)
        character(len=*), optional,      intent(in)    :: shbarrier
        integer,          optional,      intent(in)    :: nrestarts
    end subroutine generic_new

    !>  \brief  is a setter
    subroutine generic_set_indices( self, ref, ptcl, rot )
        import :: pftcc_opt
        class(pftcc_opt),  intent(inout) :: self
        integer,           intent(in)    :: ref, ptcl
        integer, optional, intent(in)    :: rot
    end subroutine generic_set_indices

    !>  \brief  is the cost function
    function generic_costfun( self, vec, D ) result( cost )
        import :: pftcc_opt
        class(pftcc_opt), intent(inout) :: self
        integer,          intent(in)    :: D
        real,             intent(in)    :: vec(D)
        real :: cost
    end function generic_costfun

    !> \brief  minimization of the costfunction
    function generic_minimize( self, irot, shvec, rxy ) result( crxy )
        import :: pftcc_opt
        class(pftcc_opt),  intent(inout) :: self
        integer, optional, intent(in)    :: irot
        real,    optional, intent(in)    :: shvec(:)
        real,    optional, intent(in)    :: rxy(:)
        real, allocatable  :: crxy(:)
    end function generic_minimize

    !> \brief  getter for number of function evaluations
    function generic_get_nevals( self ) result( nevals )
        import :: pftcc_opt
        class(pftcc_opt), intent(inout) :: self
        integer :: nevals
    end function generic_get_nevals

end interface

end module simple_pftcc_opt
