!@descr: validates flex projection-aware latent-model I/O, Fourier operations, and canonicalization
program simple_test_flex_projected_latent_model
use simple_core_module_api
use simple_flex_projected_latent_model, only: test_projected_latent_mstep_stats_io, &
    &test_projected_latent_canonicalization
use simple_flex_reconstructor_latent_ops, only: test_coupled_batch_accumulation, test_cartesian_projection_contract
use simple_flex_diffmap_rec3D,            only: test_flex_local_linear_preimage
use simple_flex_diffmap_preimage,         only: test_flex_preimage_bandwidth_decoupling
implicit none

call test_projected_latent_mstep_stats_io
call test_coupled_batch_accumulation
call test_cartesian_projection_contract
call test_projected_latent_canonicalization
call test_flex_local_linear_preimage
call test_flex_preimage_bandwidth_decoupling
call simple_end('**** SIMPLE_FLEX_PROJECTED_LATENT_MODEL TEST NORMAL STOP ****')

end program simple_test_flex_projected_latent_model
