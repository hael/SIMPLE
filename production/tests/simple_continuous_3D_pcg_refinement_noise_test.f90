module continuous_3D_pcg_refinement_noise_test
use continuous_3D_pcg_refinement_noise_gauran_test, only: run_gauran_replacement
use continuous_3D_pcg_refinement_noise_observation_test, only: run_observation_noise
use continuous_3D_pcg_refinement_noise_volume_test, only: run_added_volume_noise
implicit none
private
public :: run_volume_noise

contains

!> Run both volume-domain and observation-domain noise contracts.
subroutine run_volume_noise()
    call run_gauran_replacement()
    call run_added_volume_noise()
    call run_observation_noise()
    write(*,'(a)') 'CONTINUOUS_3D_PCG_NOISE: PASS'
end subroutine run_volume_noise

end module continuous_3D_pcg_refinement_noise_test
