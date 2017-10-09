! abstract polar Fourier cross-correlation optimisation class
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
end type  pftcc_opt

abstract interface

    !>  \brief  is a constructor
    subroutine generic_new( self, pftcc, lims, lims_init, shbarrier, nrestarts, maxits) 
        use simple_projector,        only: projector
        use simple_polarft_corrcalc, only: polarft_corrcalc
        import :: pftcc_opt
        class(pftcc_opt),                   intent(inout) :: self
        class(polarft_corrcalc),    target, intent(in)    :: pftcc
        real,                               intent(in)    :: lims(:,:)
        real,             optional,         intent(in)    :: lims_init(:,:)
        character(len=*), optional,         intent(in)    :: shbarrier
        integer,          optional,         intent(in)    :: nrestarts, maxits
    end subroutine generic_new

    !>  \brief  is a virtual setter
    subroutine generic_set_indices( self, ref, ptcl )
        import :: pftcc_opt
        class(pftcc_opt),  intent(inout) :: self
        integer,           intent(in)    :: ref, ptcl
    end subroutine generic_set_indices

    !>  \brief  is the cost function
    function generic_costfun( self, vec, D ) result( cost )
        import :: pftcc_opt
        class(pftcc_opt), intent(inout) :: self
        integer,          intent(in)    :: D
        real,             intent(in)    :: vec(D)
        real :: cost
    end function generic_costfun

    !> \brief  Virtual minimization of the costfunction
    function generic_minimize( self, irot ) result( cxy )
        import :: pftcc_opt
        class(pftcc_opt), intent(inout) :: self
        integer,          intent(out)   :: irot
        real, allocatable  :: cxy(:)
    end function generic_minimize

end interface

end module simple_pftcc_opt
