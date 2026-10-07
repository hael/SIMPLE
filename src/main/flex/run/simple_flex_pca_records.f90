!@descr: flex_pca value records shared across fitting, state inference and delivery
module simple_flex_pca_records
use simple_core_module_api, only: dp, simple_exception
use simple_image,           only: image
use simple_reconstructor,   only: reconstructor
implicit none

public :: flex_selection, flex_fit_model, flex_latent, flex_state_set, flex_latent_readout
private

!> the particle selection: project rows in ascending order
type :: flex_selection
    integer, allocatable :: pinds(:)
    integer :: nptcls = 0
  contains
    procedure :: kill => selection_kill
end type flex_selection

!> the model handles: mean, basis, prior variances, rank, noise level, and the latent point estimates of the probe subsample
type :: flex_fit_model
    type(reconstructor), allocatable :: basis_recs(:)
    real(dp),            allocatable :: eigvals(:)
    real(dp) :: sig2_eff = 0.d0  !< calibrated whitened-noise level (raw)
    real(dp) :: sig2     = 0.d0  !< floored copy the solves consume
    integer :: ncomp    = 0  !< can shrink at the M-step swap; per-fit
    type(reconstructor) :: mean_rec  !< per-fit mean copy carrying the per-fit mean scale
    real(dp), allocatable :: z(:,:)  !< (npp,ncomp at stage entry); never reallocated
  contains
    procedure :: kill => flex_fit_model_kill
end type flex_fit_model

!> the latent product of one embedding pass: per-particle latents, contrast, posterior precision,
!! residual energies, the even/odd half solutions, the per-component reliability and the prior
type :: flex_latent
    real(dp), allocatable :: z(:,:)                !< (nptcls,ncomp)
    real(dp), allocatable :: contrast(:)           !< (nptcls)
    real(dp), allocatable :: precision(:,:,:)      !< (ncomp,ncomp,nptcls) posterior precision per particle
    real(dp), allocatable :: resid_energy(:), resid_mean_energy(:)
    real(dp), allocatable :: zhalf(:,:,:)          !< (nptcls,ncomp,2) even/odd Fourier-half solutions
    real(dp), allocatable :: comp_rho(:)           !< (ncomp) per-component cross-half reliability
    real(dp), allocatable :: prior_precision(:)    !< (ncomp)
  contains
    procedure :: kill => latent_kill
end type flex_latent

!> the state set: kernel weights per state, their targets, bandwidths, effective counts, hard
!! labels, and the distances/floors the bandwidth selection rebuilds the weights from
type :: flex_state_set
    integer :: nstates = 0
    real,     allocatable :: weights(:,:)          !< (nptcls,nstates)
    real,     allocatable :: half_weights(:,:)
    real,     allocatable :: targets(:,:)          !< (ncomp,nstates)
    real,     allocatable :: bandwidths(:), neff(:)
    integer,  allocatable :: labels(:)
    real(dp), allocatable :: kdist(:,:), kfloor(:)
  contains
    procedure :: kill => state_set_kill
end type flex_state_set

type :: flex_latent_readout
    real,    allocatable :: umap_xy(:,:)
    integer, allocatable :: umap_pind(:)
  contains
    procedure :: kill => readout_kill
end type flex_latent_readout

contains

    subroutine selection_kill( self )
        class(flex_selection), intent(inout) :: self
        if( allocated(self%pinds) ) deallocate(self%pinds)
        self%nptcls = 0
    end subroutine selection_kill

    subroutine flex_fit_model_kill( self )
        class(flex_fit_model), intent(inout) :: self
        integer :: q
        if( allocated(self%basis_recs) )then
            do q = 1, size(self%basis_recs)
                call self%basis_recs(q)%dealloc_rho; call self%basis_recs(q)%kill
            end do
            deallocate(self%basis_recs)
        endif
        if( allocated(self%eigvals) ) deallocate(self%eigvals)
        call self%mean_rec%dealloc_rho
        call self%mean_rec%kill
        if( allocated(self%z) ) deallocate(self%z)
        self%ncomp = 0; self%sig2_eff = 0.d0; self%sig2 = 0.d0
    end subroutine flex_fit_model_kill

    subroutine latent_kill( self )
        class(flex_latent), intent(inout) :: self
        if( allocated(self%z) )                 deallocate(self%z)
        if( allocated(self%contrast) )          deallocate(self%contrast)
        if( allocated(self%precision) )         deallocate(self%precision)
        if( allocated(self%resid_energy) )      deallocate(self%resid_energy)
        if( allocated(self%resid_mean_energy) ) deallocate(self%resid_mean_energy)
        if( allocated(self%zhalf) )             deallocate(self%zhalf)
        if( allocated(self%comp_rho) )          deallocate(self%comp_rho)
        if( allocated(self%prior_precision) )   deallocate(self%prior_precision)
    end subroutine latent_kill

    subroutine state_set_kill( self )
        class(flex_state_set), intent(inout) :: self
        if( allocated(self%weights) )      deallocate(self%weights)
        if( allocated(self%half_weights) ) deallocate(self%half_weights)
        if( allocated(self%targets) )      deallocate(self%targets)
        if( allocated(self%bandwidths) )   deallocate(self%bandwidths)
        if( allocated(self%neff) )         deallocate(self%neff)
        if( allocated(self%labels) )       deallocate(self%labels)
        if( allocated(self%kdist) )        deallocate(self%kdist)
        if( allocated(self%kfloor) )       deallocate(self%kfloor)
        self%nstates = 0
    end subroutine state_set_kill

    subroutine readout_kill( self )
        class(flex_latent_readout), intent(inout) :: self
        if( allocated(self%umap_xy) )   deallocate(self%umap_xy)
        if( allocated(self%umap_pind) ) deallocate(self%umap_pind)
    end subroutine readout_kill

end module simple_flex_pca_records
