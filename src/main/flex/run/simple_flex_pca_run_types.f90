!@descr: flex_pca run session: values that outlive one application phase
module simple_flex_pca_run_types
use simple_core_module_api,  only: dp, simple_exception
use simple_reconstructor,    only: reconstructor
use simple_flex_pca_records, only: flex_selection, flex_fit_model, flex_latent, flex_state_set, flex_latent_readout
implicit none
private
#include "simple_local_flags.inc"

public :: flex_run_session

!> Everything a run owns between its phases. Allocated by the phases as they always were; freed
!! once, in one place, by `kill` (safe on a partially built session: every free is guarded).
type :: flex_run_session
    type(flex_selection)      :: sel      !< the particle selection
    type(flex_fit_model)      :: model    !< mean, basis, prior variances, rank, noise level
    type(flex_latent)         :: latent   !< the embedding and its statistics
    type(flex_state_set)      :: states   !< the delivered state set
    type(flex_latent_readout) :: readout  !< UMAP readout coordinates for the state figure
    logical :: sigma_loaded = .false., l_resume = .false., l_paired_states = .false.
    integer, allocatable :: deconv_labels(:)
    logical :: l_deconv_applied = .false., l_deconv_adopted = .false.
    integer :: min_neff = 0, state_axis = 0, col_sep = 1, neigs_req = 0, nkern = 0
    real(dp), allocatable :: pviews(:,:)
    logical :: l_pop_floor = .false., l_merged = .false., l_state_rec = .true.
  contains
    procedure :: clamp_state_axis => session_clamp_state_axis
    procedure :: kill => session_kill
end type flex_run_session

contains

    !> Free everything the session may hold: model handles first (a resume never built them), then
    !! the embedding, the state table and the readout.
    !> the state axis can never exceed the rank or the kernel dimension; applied after every rank change
    subroutine session_clamp_state_axis( self )
        class(flex_run_session), intent(inout) :: self
        if( self%state_axis > 0 ) self%state_axis = min(self%state_axis, min(self%model%ncomp, self%nkern))
    end subroutine session_clamp_state_axis

    subroutine session_kill( self )
        class(flex_run_session), intent(inout) :: self
        call self%sel%kill
        call self%model%kill
        call self%latent%kill
        call self%states%kill
        call self%readout%kill
        if( allocated(self%deconv_labels) ) deallocate(self%deconv_labels)
        if( allocated(self%pviews) )        deallocate(self%pviews)
    end subroutine session_kill

end module simple_flex_pca_run_types
