! abstract polar Fourier cross-correlation optimisation class
module simple_pftcc_opt
implicit none

public :: pftcc_opt
private

type, abstract :: pftcc_opt
  contains
    procedure(generic_new),          deferred :: new
    procedure(generic_set_indices),  deferred :: set_indices
    procedure(generic_set_inipop),   deferred :: set_inipop
    procedure(generic_costfun),      deferred :: costfun
    procedure(generic_minimize),     deferred :: minimize
    procedure(generic_get_nevals),   deferred :: get_nevals
    procedure(generic_get_peaks),    deferred :: get_peaks
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
    subroutine generic_set_indices( self, ref, ptcl, rot, state )
        import :: pftcc_opt
        class(pftcc_opt),  intent(inout) :: self
        integer,           intent(in)    :: ref, ptcl
        integer, optional, intent(in)    :: rot, state
    end subroutine generic_set_indices

    !>  \brief  is a setter
    subroutine generic_set_inipop( self, inipop )
        import :: pftcc_opt
        class(pftcc_opt),  intent(inout) :: self
        real,              intent(in)    :: inipop(:,:)
    end subroutine generic_set_inipop

    !>  \brief  is the cost function
    function generic_costfun( self, vec, D ) result( cost )
        import :: pftcc_opt
        class(pftcc_opt), intent(inout) :: self
        integer,          intent(in)    :: D
        real,             intent(in)    :: vec(D)
        real :: cost
    end function generic_costfun

    !> \brief  Virtual minimization of the costfunction
    function generic_minimize( self, irot, shvec, rxy, fromto ) result( crxy )
        import :: pftcc_opt
        class(pftcc_opt),  intent(inout) :: self
        integer, optional, intent(inout) :: irot
        real,    optional, intent(in)    :: shvec(:)
        real,    optional, intent(in)    :: rxy(:)
        integer, optional, intent(in)    :: fromto(2)
        real, allocatable  :: crxy(:)
    end function generic_minimize

    !> \brief  Virtual getter for number of function evaluations
    function generic_get_nevals( self ) result( nevals )
        import :: pftcc_opt
        class(pftcc_opt), intent(inout) :: self
        integer :: nevals
    end function generic_get_nevals

    !> \brief  getter for number of function evaluations
    subroutine generic_get_peaks( self, peaks )
        import :: pftcc_opt
        class(pftcc_opt), intent(inout) :: self
        real, allocatable, intent(out)  :: peaks(:,:)
    end subroutine generic_get_peaks

end interface

end module simple_pftcc_opt
