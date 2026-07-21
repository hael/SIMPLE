!@descr: reusable result container for a fitted projected latent volume model
module simple_projected_latent_result
use simple_core_module_api
implicit none

public :: projected_latent_fit_result
private

type :: projected_latent_fit_result
    integer :: nptcls = 0
    integer :: ncomp  = 0
    integer,  allocatable :: pinds(:)
    real(dp), allocatable :: z(:,:)
    real(dp), allocatable :: resid_energy(:), resid_mean_energy(:)
    real(dp), allocatable :: eigvals(:)
contains
    procedure :: kill => kill_projected_latent_fit_result
end type projected_latent_fit_result

contains

    subroutine kill_projected_latent_fit_result( self )
        class(projected_latent_fit_result), intent(inout) :: self
        if( allocated(self%pinds) ) deallocate(self%pinds)
        if( allocated(self%z) ) deallocate(self%z)
        if( allocated(self%resid_energy) ) deallocate(self%resid_energy)
        if( allocated(self%resid_mean_energy) ) deallocate(self%resid_mean_energy)
        if( allocated(self%eigvals) ) deallocate(self%eigvals)
        self%nptcls = 0
        self%ncomp  = 0
    end subroutine kill_projected_latent_fit_result

end module simple_projected_latent_result
