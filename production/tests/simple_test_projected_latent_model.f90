!@descr: validates projected EM-PCA sufficient-statistics I/O and final canonicalization
program simple_test_projected_latent_model
use simple_core_module_api
use simple_projected_latent_model, only: test_projected_latent_mstep_stats_io, &
    &test_projected_latent_canonicalization
use simple_reconstructor_latent_ops, only: test_coupled_batch_accumulation, test_cartesian_projection_contract
implicit none

call test_projected_latent_mstep_stats_io
call test_coupled_batch_accumulation
call test_cartesian_projection_contract
call test_projected_latent_canonicalization
call simple_end('**** SIMPLE_PROJECTED_LATENT_MODEL TEST NORMAL STOP ****')

end program simple_test_projected_latent_model
