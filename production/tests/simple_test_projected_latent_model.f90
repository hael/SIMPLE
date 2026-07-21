!@descr: validates projected EM-PCA sufficient-statistics I/O and final canonicalization
program simple_test_projected_latent_model
use simple_core_module_api
use simple_projected_latent_model, only: test_projected_latent_mstep_stats_io, &
    &test_projected_latent_canonicalization
implicit none

call test_projected_latent_mstep_stats_io
call test_projected_latent_canonicalization
call simple_end('**** SIMPLE_PROJECTED_LATENT_MODEL TEST NORMAL STOP ****')

end program simple_test_projected_latent_model
